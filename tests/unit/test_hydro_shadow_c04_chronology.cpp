// C04-02 full-physics chronology fixtures (operator level).
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// The smooth source operator is the PRODUCTION physics::BaryonicSourceStage
// over the production physics::ThermochemistryOperator (fixed hybrid network,
// composition sidecar), replayed on a PRIVATE single-region SimulationState
// exactly as SourceRuntime::ShadowSmoothReplay does. The hydro operator is the
// production HydroCoreSolver (CTU/PPM). The fixture builds the production
// chronology S(dt/2) H(dt) S(dt/2) + captured event forcing on an 8^3
// periodic patch and checks that the shadow:
//   * replays S1 bit-identically (replay identity = 0),
//   * completes the refined Strang chronology (kCoupledFullPhysicsStrang),
//   * transports composition privately (production sidecar untouched),
//   * fails closed with a SPECIFIC counted reason for every unsupported case.
// The workflow-level wiring (SourceRuntime observer, real T01 dispatch) is
// covered by tests/integration/test_hydro_shadow_c04_workflow_strang.cpp.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <optional>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/baryonic_source_stage.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"

namespace {

using namespace cosmosim;
namespace h = cosmosim::hydro;
namespace p = cosmosim::physics;

constexpr double k_gamma = 5.0 / 3.0;
constexpr std::size_t k_cells = 512U;
constexpr double k_interval_seconds = 3.15576e13;  // ~1 Myr accepted interval
constexpr double k_velocity_si_per_code = 1.0e3;   // km/s code velocity

p::ThermochemistryUnitBridge units() {
  p::ThermochemistryUnitBridge u;
  u.density_code_to_cgs = 1.6726e-24;
  u.velocity_code_to_cm_s = 1.0e5;
  u.length_code_to_cm = 3.0857e21;
  u.mass_code_to_g = 1.989e33;
  return u;
}

// Test-side twin of SourceRuntime::ShadowSmoothReplay (same procedure, same
// production operator); status_override injects the fail-closed paths.
class OperatorReplay final : public h::HydroShadowSmoothSourceReplay {
 public:
  OperatorReplay(const p::ThermochemistryOperator& chemistry, p::BaryonicSourceStageConfig config)
      : m_chemistry(chemistry), m_config(config) {}
  std::optional<h::HydroShadowReplayStatus> status_override;
  mutable std::uint64_t calls = 0U;

  [[nodiscard]] h::HydroShadowReplayStatus advanceSmooth(
      h::HydroShadowPrivateRegion& region, double f0, double f1) const override {
    ++calls;
    if (status_override.has_value()) return *status_override;
    const std::size_t n = region.cell_count;
    core::SimulationState mini;
    mini.cells.resize(n);
    mini.gas_cells.resize(n);
    mini.gas_composition = region.composition;
    const double volume = region.cell_volume_comoving;
    const double factor = h::temperatureKelvinPerInternalEnergyCode(k_gamma, k_velocity_si_per_code);
    for (std::size_t i = 0; i < n; ++i) {
      const double rho = region.mass_density[i];
      mini.gas_cells.gas_cell_id[i] = region.gas_cell_ids[i];
      mini.gas_cells.density_code[i] = rho;
      mini.cells.mass_code[i] = rho * volume;
      mini.gas_cells.velocity_x_peculiar[i] = region.momentum_x[i] / rho;
      mini.gas_cells.velocity_y_peculiar[i] = region.momentum_y[i] / rho;
      mini.gas_cells.velocity_z_peculiar[i] = region.momentum_z[i] / rho;
      mini.gas_cells.internal_energy_code[i] = region.specific_internal_energy[i];
      mini.gas_cells.metal_mass_code[i] = region.metal_density[i] * volume;
      mini.gas_cells.temperature_code[i] = factor * region.specific_internal_energy[i];
    }
    p::SourceInterval interval{.physical_seconds = (f1 - f0) * k_interval_seconds,
                               .evaluation_scale_factor = 1.0, .accepted_step_index = 1U};
    p::BaryonicSourceStage stage(m_chemistry, m_config);
    try {
      stage.prepare(mini, interval, {});
    } catch (const p::SourceIntervalIncompleteError&) {
      return h::HydroShadowReplayStatus::kIntervalIncomplete;
    }
    (void)stage.commit(mini);
    for (std::size_t i = 0; i < n; ++i) {
      const double rho = region.mass_density[i];
      region.specific_internal_energy[i] = mini.gas_cells.internal_energy_code[i];
      region.metal_density[i] = mini.gas_cells.metal_mass_code[i] / volume;
      const double m2 = region.momentum_x[i] * region.momentum_x[i] + region.momentum_y[i] * region.momentum_y[i] +
          region.momentum_z[i] * region.momentum_z[i];
      region.total_energy[i] = rho * region.specific_internal_energy[i] + 0.5 * m2 / rho;
    }
    region.composition = std::move(mini.gas_composition);
    return h::HydroShadowReplayStatus::kCompleted;
  }
  [[nodiscard]] std::string_view modelLabel() const noexcept override { return "test_operator_replay"; }

 private:
  const p::ThermochemistryOperator& m_chemistry;
  p::BaryonicSourceStageConfig m_config;
};

struct NullSink final : h::HydroFaceMassFluxSink {
  std::uint64_t records = 0U;
  void recordFaceMassTransfer(std::size_t, std::size_t, std::size_t, double) override { ++records; }
};

struct Fixture {
  p::CompositionLayout layout = p::CompositionLayout::make(
      {.network = p::CompositionNetworkKind::kHybridLight, .dust = p::CompositionDustKind::kTwoReservoirV1});
  h::HydroPatchGeometry geometry;
  core::SimulationState production;  // canonical rows (composition authority)
  std::vector<std::uint32_t> rows;
  std::vector<std::size_t> cells;
  std::vector<std::size_t> faces;
  double volume = 0.0;
};

[[nodiscard]] Fixture makeFixture() {
  Fixture f;
  f.geometry = h::makeCartesianPatchGeometry({.nx = 8U, .ny = 8U, .nz = 8U,
                                              .cell_width_x_comoving = 0.125,
                                              .cell_width_y_comoving = 0.125,
                                              .cell_width_z_comoving = 0.125});
  h::appendCartesianBoundaryGhostFaces(f.geometry, h::HydroBoundaryKind::kPeriodic);
  (void)h::completeReconstructionStencils(f.geometry);
  (void)h::completePeriodicReconstructionStencils(f.geometry);
  f.volume = f.geometry.cell_volume_comoving;
  f.production.resizeCells(k_cells);
  for (std::size_t r = 0; r < k_cells; ++r) {
    const auto ijk = f.geometry.cellIjk(r);
    const double rho = 1.0 + 0.02 * std::sin(0.785398 * static_cast<double>(ijk[0] + 2U * ijk[1]));
    f.production.gas_cells.gas_cell_id[r] = 5000U + r;
    f.production.gas_cells.density_code[r] = rho;
    f.production.cells.mass_code[r] = rho * f.volume;
    f.production.gas_cells.velocity_x_peculiar[r] = 0.3;
    f.production.gas_cells.velocity_y_peculiar[r] = -0.2;
    f.production.gas_cells.internal_energy_code[r] = 200.0 + 5.0 * static_cast<double>(ijk[2]);
    f.production.gas_cells.metal_mass_code[r] = 0.0134 * f.production.cells.mass_code[r];
    f.production.gas_cells.temperature_code[r] =
        h::temperatureKelvinPerInternalEnergyCode(k_gamma, k_velocity_si_per_code) *
        f.production.gas_cells.internal_energy_code[r];
  }
  p::enableGasComposition(f.production, f.layout, {});
  f.rows.resize(k_cells);
  std::iota(f.rows.begin(), f.rows.end(), 0U);
  f.cells.resize(k_cells);
  std::iota(f.cells.begin(), f.cells.end(), 0U);
  f.faces.resize(f.geometry.faces.size());
  std::iota(f.faces.begin(), f.faces.end(), 0U);
  return f;
}

[[nodiscard]] h::HydroShadowPrivateRegion snapshot(const core::SimulationState& state, double volume) {
  h::HydroShadowPrivateRegion region;
  region.resize(state.cells.size());
  region.cell_volume_comoving = volume;
  region.composition = state.gas_composition.emptyLike(state.cells.size(), false);
  for (std::size_t r = 0; r < state.cells.size(); ++r) {
    const double rho = state.gas_cells.density_code[r];
    region.gas_cell_ids[r] = state.gas_cells.gas_cell_id[r];
    region.mass_density[r] = rho;
    region.momentum_x[r] = rho * state.gas_cells.velocity_x_peculiar[r];
    region.momentum_y[r] = rho * state.gas_cells.velocity_y_peculiar[r];
    region.momentum_z[r] = rho * state.gas_cells.velocity_z_peculiar[r];
    region.specific_internal_energy[r] = state.gas_cells.internal_energy_code[r];
    region.metal_density[r] = rho * state.gas_cells.metal_mass_code[r] / state.cells.mass_code[r];
    const double v2 = state.gas_cells.velocity_x_peculiar[r] * state.gas_cells.velocity_x_peculiar[r] +
        state.gas_cells.velocity_y_peculiar[r] * state.gas_cells.velocity_y_peculiar[r] +
        state.gas_cells.velocity_z_peculiar[r] * state.gas_cells.velocity_z_peculiar[r];
    region.total_energy[r] = rho * (region.specific_internal_energy[r] + 0.5 * v2);
    region.composition.unpackRow(r, state.gas_composition.packRow(r), region.gas_cell_ids[r]);
    region.composition.mean_molecular_weight[r] = state.gas_composition.mean_molecular_weight[r];
  }
  return region;
}

[[nodiscard]] h::HydroConservedStateSoa soaFrom(const h::HydroShadowPrivateRegion& region,
                                                const h::HydroPatchGeometry& geometry) {
  h::HydroConservedStateSoa soa(geometry.totalCellStorageCount());
  for (std::size_t r = 0; r < region.cell_count; ++r) {
    soa.storeCell(r, h::HydroConservedState{.mass_density_comoving = region.mass_density[r],
                                            .momentum_density_x_comoving = region.momentum_x[r],
                                            .momentum_density_y_comoving = region.momentum_y[r],
                                            .momentum_density_z_comoving = region.momentum_z[r],
                                            .total_energy_density_comoving = region.total_energy[r],
                                            .metal_mass_density_comoving = region.metal_density[r]});
  }
  h::fillHydroBoundaryGhostCells(soa, geometry, k_gamma);
  return soa;
}

void regionFromSoa(const h::HydroConservedStateSoa& soa, h::HydroShadowPrivateRegion& region) {
  for (std::size_t r = 0; r < region.cell_count; ++r) {
    const auto c = soa.loadCell(r);
    region.mass_density[r] = c.mass_density_comoving;
    region.momentum_x[r] = c.momentum_density_x_comoving;
    region.momentum_y[r] = c.momentum_density_y_comoving;
    region.momentum_z[r] = c.momentum_density_z_comoving;
    region.metal_density[r] = c.metal_mass_density_comoving;
    const double m2 = c.momentum_density_x_comoving * c.momentum_density_x_comoving +
        c.momentum_density_y_comoving * c.momentum_density_y_comoving +
        c.momentum_density_z_comoving * c.momentum_density_z_comoving;
    region.specific_internal_energy[r] = (c.total_energy_density_comoving - 0.5 * m2 / c.mass_density_comoving) /
        c.mass_density_comoving;
    region.total_energy[r] = c.total_energy_density_comoving;
  }
}

struct Chronology {
  std::array<h::HydroShadowPrivateRegion, 7> points;
  h::HydroConservedStateSoa hydro_input;
  h::HydroConservedStateSoa hydro_output;
};

// Production chronology built with the SAME operators: S1 = replay over
// [0, 1/2], H = one production CTU/PPM step, S2 = replay over [1/2, 1],
// discrete events = an immutable forcing increment on two cells.
[[nodiscard]] Chronology productionChronology(const Fixture& f, const OperatorReplay& replay,
                                              const h::HydroUpdateContext& update) {
  Chronology c;
  c.points[0] = snapshot(f.production, f.volume);
  c.points[1] = c.points[0];  // no element diffusion / SGS shear in this deck
  c.points[2] = c.points[1];
  assert(replay.advanceSmooth(c.points[2], 0.0, 0.5) == h::HydroShadowReplayStatus::kCompleted);
  c.hydro_input = soaFrom(c.points[2], f.geometry);
  c.hydro_output = c.hydro_input;
  h::HydroNumericsPolicy numerics;
  numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  const h::HydroCoreSolver solver(k_gamma, numerics);
  const auto courant = h::hydroTransportCourantFactors(update.dt_code, 1.0, 0.125, 0.125, 0.125);
  const h::HydroReconstructionPolicy policy{.dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
                                            .adiabatic_index = k_gamma};
  const h::PpmReconstruction ppm(policy);
  const h::HllcRiemannSolver riemann;
  h::HydroScratchBuffers scratch;
  h::HydroSourceContext source_context{.update = update};
  solver.advancePatchActiveSetBatchedWithScratch(
      c.hydro_output, f.geometry, {.active_cells = f.cells, .active_faces = f.faces}, update, ppm, riemann, {},
      source_context, scratch, nullptr, solver.requiredWorkspaceFeatures({}));
  c.points[3] = c.points[2];
  regionFromSoa(c.hydro_output, c.points[3]);
  c.points[4] = c.points[3];
  assert(replay.advanceSmooth(c.points[4], 0.5, 1.0) == h::HydroShadowReplayStatus::kCompleted);
  c.points[5] = c.points[4];
  c.points[6] = c.points[5];
  for (const std::size_t cell : {std::size_t{7U}, std::size_t{300U}}) {
    c.points[6].specific_internal_energy[cell] *= 1.25;  // SN thermal deposit (captured)
    c.points[6].total_energy[cell] += c.points[6].mass_density[cell] * 0.25 * c.points[5].specific_internal_energy[cell];
  }
  return c;
}

[[nodiscard]] h::HydroShadowPatchCapture captureSpec(const Fixture& f, const Chronology& c,
                                                     const h::HydroUpdateContext& update,
                                                     const h::HydroSourceContext& source_context,
                                                     std::span<const std::size_t> active) {
  h::HydroShadowPatchCapture spec;
  spec.identity = h::HydroShadowCandidate{.patch_id = 1U, .level = 0,
                                          .region_key = h::hydroShadowPhysicalRegionKey(0, {0, 0, 0}, {8U, 8U, 8U})};
  spec.before = &c.hydro_input;
  spec.geometry = &f.geometry;
  spec.active = {.active_cells = active, .active_faces = f.faces};
  spec.update = update;
  spec.production_numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  spec.production_reconstruction = h::HydroReconstructionScheme::kPiecewiseParabolic;
  spec.production_policy = h::HydroReconstructionPolicy{.adiabatic_index = k_gamma};
  spec.production_riemann = h::HydroRiemannSolverKind::kHllc;
  spec.source_context = &source_context;
  spec.composition = &f.production.gas_composition;
  spec.composition_real_rows = f.rows;
  spec.temperature_k_per_internal_energy_code =
      h::temperatureKelvinPerInternalEnergyCode(k_gamma, k_velocity_si_per_code);
  return spec;
}

void recordAll(h::HydroShadowReference& shadow, std::uint64_t key, const Chronology& c,
               std::optional<std::size_t> skip = std::nullopt) {
  for (std::size_t point = 0; point < c.points.size(); ++point) {
    if (skip.has_value() && *skip == point) continue;
    shadow.recordChronologyPoint(key, static_cast<h::HydroShadowChronologyPoint>(point), c.points[point]);
  }
}

[[nodiscard]] std::vector<double> compositionLanes(const core::SimulationState& state) {
  std::vector<double> lanes;
  for (std::size_t l = 0; l < state.gas_composition.laneCount(); ++l) {
    const auto lane = state.gas_composition.lane(l);
    lanes.insert(lanes.end(), lane.begin(), lane.end());
  }
  return lanes;
}

struct Run {
  h::HydroShadowReport report;
  std::vector<h::HydroShadowRegionResult> results;
};

[[nodiscard]] Run runCase(const Fixture& f, const OperatorReplay& replay, const Chronology& c,
                          std::span<const std::size_t> active, std::optional<std::size_t> skip_point,
                          bool corrupt_ids, bool provide_replay) {
  const h::HydroUpdateContext update{.dt_code = 1.0e-4, .time_begin_code = 0.0, .time_end_code = 1.0e-4};
  NullSink sink;
  h::HydroSourceContext source_context{.update = update};
  source_context.face_mass_flux_sink = &sink;  // composition transport active in production
  h::HydroShadowReference shadow({.enabled = true, .cadence = 1U, .max_patches_per_rank = 1U,
                                  .max_storage_cells_per_rank = 1U << 16U,
                                  .max_workspace_bytes_per_rank = 1ULL << 30U});
  shadow.beginEpoch(0U);
  auto spec = captureSpec(f, c, update, source_context, active);
  shadow.setGlobalCandidates({spec.identity});
  const std::uint64_t key = h::hydroShadowEffectiveKey(spec.identity);
  std::vector<std::uint64_t> ids(k_cells);
  for (std::size_t r = 0; r < k_cells; ++r) ids[r] = 5000U + r;
  assert(shadow.openChronology(spec.identity, ids, f.production.gas_composition.laneCount()));
  auto sample = shadow.capture(spec);
  assert(sample);
  shadow.observeProduction(std::move(sample), c.hydro_output);
  shadow.finalizePendingEndpoints();
  if (corrupt_ids) {
    auto changed = c.points[0];
    changed.gas_cell_ids[3] = 999999U;  // a regrid/migration changed the region
    shadow.recordChronologyPoint(key, h::HydroShadowChronologyPoint::kPreHalfBegin, changed);
  } else {
    recordAll(shadow, key, c, skip_point);
  }
  shadow.finalizeChronologies(provide_replay ? &replay : nullptr);
  return Run{shadow.report(), std::vector<h::HydroShadowRegionResult>(shadow.regionResults().begin(),
                                                                       shadow.regionResults().end())};
}

void checkCoupledChronology() {
  const Fixture f = makeFixture();
  p::ThermochemistryOperatorConfig config;
  config.solver.backend = p::ThermochemistryBackend::kFixedHybrid;
  config.units = units();
  config.solver.allow_unqualified_rates = true;
  const p::ThermochemistryOperator chemistry(config, f.layout);
  const OperatorReplay replay(chemistry, p::BaryonicSourceStageConfig{});
  const h::HydroUpdateContext update{.dt_code = 1.0e-4, .time_begin_code = 0.0, .time_end_code = 1.0e-4};
  const Chronology c = productionChronology(f, replay, update);
  const auto production_composition = compositionLanes(f.production);
  const Run run = runCase(f, replay, c, f.cells, std::nullopt, false, true);
  assert(run.report.chronology_started == 1U);
  assert(run.report.completed_chronology_regions == 1U);
  // Same operator, same data: the private S1 reproduces production exactly.
  assert(run.report.chronology_replay_identity_linf == 0.0);
  assert(run.report.completed == 1U);
  assert(run.report.completed_composition_patches == 0U);  // no post-commit snapshot path here
  bool coupled = false;
  for (const auto& result : run.results) {
    if (result.comparison_class == h::HydroShadowComparisonClass::kCoupledFullPhysicsStrang) {
      coupled = true;
      assert(result.compared_cells == k_cells);
      assert(result.compared_composition_cells == k_cells);
      for (const auto& metric : result.components) assert(std::isfinite(metric.weighted_abs_difference));
    }
  }
  assert(coupled);
  // The private composition never aliases the production sidecar.
  assert(compositionLanes(f.production) == production_composition);
  // Four private smooth calls + one identity probe per chronology.
  assert(replay.calls >= 5U);
}

void checkFailClosedReasons() {
  const Fixture f = makeFixture();
  p::ThermochemistryOperatorConfig config;
  config.solver.backend = p::ThermochemistryBackend::kFixedHybrid;
  config.units = units();
  config.solver.allow_unqualified_rates = true;
  const p::ThermochemistryOperator chemistry(config, f.layout);
  OperatorReplay replay(chemistry, p::BaryonicSourceStageConfig{});
  const h::HydroUpdateContext update{.dt_code = 1.0e-4, .time_begin_code = 0.0, .time_end_code = 1.0e-4};
  const Chronology c = productionChronology(f, replay, update);
  // Missing chronology point (e.g. a dispatch boundary not observed).
  assert(runCase(f, replay, c, f.cells, std::size_t{5U}, false, true).report.chronology_skipped_missing_point == 1U);
  // Region identity changed between points (regrid/migration).
  assert(runCase(f, replay, c, f.cells, std::nullopt, true, true).report.chronology_skipped_region_changed == 1U);
  // No smooth-source owner/replay lent.
  assert(runCase(f, replay, c, f.cells, std::nullopt, false, false).report.chronology_skipped_no_replay_operator == 1U);
  // Private integrator would need a retry: refuse, never truncate.
  replay.status_override = h::HydroShadowReplayStatus::kIntervalIncomplete;
  assert(runCase(f, replay, c, f.cells, std::nullopt, false, true).report.chronology_skipped_replay_incomplete == 1U);
  replay.status_override = h::HydroShadowReplayStatus::kModelUnavailable;
  assert(runCase(f, replay, c, f.cells, std::nullopt, false, true).report.chronology_skipped_replay_unavailable == 1U);
  replay.status_override.reset();
  // Partial scheduler activity: production S advances only the active rows,
  // so the private replay is restricted to them (no refusal any more). The
  // production chronology of this fixture was all-active, so the comparison
  // completes with a nonzero (honestly reported) disagreement.
  std::vector<std::size_t> half(f.cells.begin(), f.cells.begin() + 256);
  const Run partial = runCase(f, replay, c, half, std::nullopt, false, true);
  assert(partial.report.chronology_skipped_partial_activity == 0U);
  assert(partial.report.completed_chronology_regions == 1U);
  assert(partial.report.completed == 1U);
  assert(partial.report.completed_partial_active_patches == 1U);
}

// A sample whose production face sink is installed but whose composition row
// map is missing is refused as external passives (never silently compared).
void checkCompositionMapRequired() {
  const Fixture f = makeFixture();
  const h::HydroUpdateContext update{.dt_code = 1.0e-4, .time_begin_code = 0.0, .time_end_code = 1.0e-4};
  NullSink sink;
  h::HydroSourceContext source_context{.update = update};
  source_context.face_mass_flux_sink = &sink;
  Chronology c;
  c.hydro_input = soaFrom(snapshot(f.production, f.volume), f.geometry);
  h::HydroShadowReference shadow({.enabled = true, .cadence = 1U});
  shadow.beginEpoch(0U);
  auto spec = captureSpec(f, c, update, source_context, f.cells);
  spec.composition_real_rows = {};
  shadow.setGlobalCandidates({spec.identity});
  assert(!shadow.capture(spec));
  assert(shadow.report().skipped_external_passives == 1U);
  assert(sink.records == 0U);
}

}  // namespace

int main() {
  checkCompositionMapRequired();
  checkCoupledChronology();
  checkFailClosedReasons();
  return 0;
}
