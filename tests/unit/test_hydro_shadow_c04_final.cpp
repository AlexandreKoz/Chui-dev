// C04 final source-closure fixtures (operator level).
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// Covers the contracts added by the final C04 source campaign:
//   * C04-05 containment: capture / observeProduction / setGlobalCandidates
//     are noexcept, so an optional diagnostic cannot abort a production sweep.
//   * C04-02/03 remote donors: an 8^3 CTU/PPM sample whose x-lower ghost
//     slots are REMOTE same-level imports (no local composition row) is HELD
//     when the caller promises the production donor-fraction table, is not
//     compared by finalizePendingEndpoints(false), and after
//     provideRemoteDonorFractions() resolves every remote slot with no
//     receiver-fraction fallback. Resolving the same donors through LOCAL rows
//     needs no hold and yields the identical private composition metric.
//   * Without the table (or with a malformed one) the counted receiver-
//     fraction fallback is used, never a silent substitution.
//   * C04-02 parity: the private composition is reconciled with the private
//     hydro mass exactly like CompositionFluxAccumulator::commit (finite,
//     round-off-sized reconciliation for a floor-free deck).
//   * Rejected attempts: held work is counted by unpublishedResultCount()
//     before beginEpoch() discards it (the runtime carries that count).
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"
#include "cosmosim/physics/baryon_composition.hpp"

namespace {

using namespace cosmosim;
namespace h = cosmosim::hydro;
namespace p = cosmosim::physics;

constexpr double k_gamma = 5.0 / 3.0;
constexpr std::size_t k_cells = 512U;
constexpr std::uint64_t k_remote_id_base = 900000U;

struct NullSink final : h::HydroFaceMassFluxSink {
  void recordFaceMassTransfer(std::size_t, std::size_t, std::size_t, double) override {}
};

struct Fixture {
  p::CompositionLayout layout = p::CompositionLayout::make(
      {.network = p::CompositionNetworkKind::kHybridLight, .dust = p::CompositionDustKind::kTwoReservoirV1});
  h::HydroPatchGeometry geometry;
  core::SimulationState production;
  std::vector<std::uint32_t> rows;
  std::vector<std::size_t> cells;
  std::vector<std::size_t> faces;
  h::HydroConservedStateSoa before;
  h::HydroConservedStateSoa after;
  std::vector<h::HydroShadowGhostProvenance> provenance;
  std::vector<std::uint64_t> ghost_ids;       // remote origin id per marked slot, else 0
  std::vector<std::uint32_t> local_ghost_rows;  // same donors through local rows
  std::vector<std::size_t> marked_slots;
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
  const double volume = f.geometry.cell_volume_comoving;
  f.production.resizeCells(k_cells);
  f.before = h::HydroConservedStateSoa(f.geometry.totalCellStorageCount());
  for (std::size_t r = 0; r < k_cells; ++r) {
    const auto ijk = f.geometry.cellIjk(r);
    const double rho = 1.0 + 0.05 * std::sin(0.785398 * static_cast<double>(ijk[0] + ijk[1]));
    const double vx = 0.3;
    const double e = 200.0 + 3.0 * static_cast<double>(ijk[2]);
    f.production.gas_cells.gas_cell_id[r] = 7000U + r;
    f.production.gas_cells.density_code[r] = rho;
    f.production.cells.mass_code[r] = rho * volume;
    f.production.gas_cells.velocity_x_peculiar[r] = vx;
    f.production.gas_cells.internal_energy_code[r] = e;
    f.production.gas_cells.metal_mass_code[r] = 0.0134 * rho * volume;
    f.before.storeCell(r, h::HydroConservedState{.mass_density_comoving = rho,
                                                 .momentum_density_x_comoving = rho * vx,
                                                 .total_energy_density_comoving = rho * (e + 0.5 * vx * vx),
                                                 .metal_mass_density_comoving = 0.0134 * rho});
  }
  p::enableGasComposition(f.production, f.layout, {});
  h::fillHydroBoundaryGhostCells(f.before, f.geometry, k_gamma);
  f.rows.resize(k_cells);
  std::iota(f.rows.begin(), f.rows.end(), 0U);
  f.cells.resize(k_cells);
  std::iota(f.cells.begin(), f.cells.end(), 0U);
  f.faces.resize(f.geometry.faces.size());
  std::iota(f.faces.begin(), f.faces.end(), 0U);
  // Production result of the same interval (one CTU/PPM step).
  f.after = f.before;
  h::HydroNumericsPolicy numerics;
  numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  const h::HydroCoreSolver solver(k_gamma, numerics);
  const h::HydroUpdateContext update{.dt_code = 1.0e-3, .time_begin_code = 0.0, .time_end_code = 1.0e-3};
  const auto courant = h::hydroTransportCourantFactors(update.dt_code, 1.0, 0.125, 0.125, 0.125);
  const h::HydroReconstructionPolicy policy{.dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
                                            .adiabatic_index = k_gamma};
  const h::PpmReconstruction ppm(policy);
  const h::HllcRiemannSolver riemann;
  h::HydroScratchBuffers scratch;
  h::HydroSourceContext source_context{.update = update};
  solver.advancePatchActiveSetBatchedWithScratch(
      f.after, f.geometry, {.active_cells = f.cells, .active_faces = f.faces}, update, ppm, riemann, {},
      source_context, scratch, nullptr, solver.requiredWorkspaceFeatures({}));
  // x-lower face slots are declared REMOTE same-level imports.
  const std::size_t slots = f.geometry.ghost_cells.size();
  f.provenance.resize(slots);
  f.ghost_ids.assign(slots, 0U);
  f.local_ghost_rows.assign(slots, core::kInvalidGasCellRow);
  for (std::size_t slot = 0; slot < slots; ++slot) {
    f.provenance[slot] = h::hydroShadowDefaultGhostProvenance(f.geometry, slot);
    const auto& ghost = f.geometry.ghost_cells[slot];
    if (ghost.axis == h::HydroFaceAxis::kX && ghost.side == h::HydroFaceSide::kLower &&
        ghost.transverse_axis == 3U && ghost.source_real_cell < k_cells) {
      f.provenance[slot] = h::HydroShadowGhostProvenance::kImportedSameLevelRemote;
      f.ghost_ids[slot] = k_remote_id_base + slot;
      f.local_ghost_rows[slot] = static_cast<std::uint32_t>(ghost.source_real_cell);
      f.marked_slots.push_back(slot);
    }
  }
  assert(!f.marked_slots.empty());
  return f;
}

// Production donor table for the marked slots: ids sorted, lane fractions of
// the periodic source cell (what the donor-fraction exchange would deliver).
void donorTable(const Fixture& f, std::vector<std::uint64_t>& ids, std::vector<double>& fractions) {
  const auto& sidecar = f.production.gas_composition;
  const auto& d = sidecar.descriptor();
  const std::size_t lanes = sidecar.laneCount();
  for (const std::size_t slot : f.marked_slots) {
    const std::size_t row = f.geometry.ghost_cells[slot].source_real_cell;
    double mass = 0.0;
    for (std::uint32_t l = 0; l < d.element_lane_count; ++l) mass += sidecar.value(d.element_lane_begin + l, row);
    for (std::uint32_t l = 0; l < d.dust_lane_count; ++l) mass += sidecar.value(d.dust_lane_begin + l, row);
    assert(mass > 0.0);
    ids.push_back(f.ghost_ids[slot]);
    for (std::size_t l = 0; l < lanes; ++l) fractions.push_back(sidecar.value(l, row) / mass);
  }
  assert(std::is_sorted(ids.begin(), ids.end()));
}

[[nodiscard]] h::HydroShadowPatchCapture captureSpec(const Fixture& f, const h::HydroSourceContext& source_context,
                                                     bool remote_ids, bool local_rows, bool promise) {
  h::HydroShadowPatchCapture spec;
  spec.identity = h::HydroShadowCandidate{.patch_id = 3U, .level = 0,
                                          .region_key = h::hydroShadowPhysicalRegionKey(0, {8, 0, 0}, {8U, 8U, 8U})};
  spec.before = &f.before;
  spec.geometry = &f.geometry;
  spec.active = {.active_cells = f.cells, .active_faces = f.faces};
  spec.update = source_context.update;
  spec.production_numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  spec.production_reconstruction = h::HydroReconstructionScheme::kPiecewiseParabolic;
  spec.production_policy = h::HydroReconstructionPolicy{.adiabatic_index = k_gamma};
  spec.production_riemann = h::HydroRiemannSolverKind::kHllc;
  spec.source_context = &source_context;
  spec.ghost_provenance = f.provenance;
  spec.composition = &f.production.gas_composition;
  spec.composition_real_rows = f.rows;
  if (local_rows) spec.composition_ghost_rows = f.local_ghost_rows;
  if (remote_ids) spec.composition_ghost_gas_cell_ids = f.ghost_ids;
  spec.remote_donor_fractions_will_be_provided = promise;
  return spec;
}

struct Outcome {
  h::HydroShadowReport report;
  std::uint64_t completed_before_table = 0U;
  std::uint64_t unpublished_before_table = 0U;
  h::HydroShadowMetric composition{};
};

enum class Table { kNone, kValid, kMalformed };

[[nodiscard]] Outcome run(const Fixture& f, bool remote_ids, bool local_rows, bool promise, Table table) {
  const h::HydroUpdateContext update{.dt_code = 1.0e-3, .time_begin_code = 0.0, .time_end_code = 1.0e-3};
  NullSink production_sink;
  h::HydroSourceContext source_context{.update = update};
  source_context.face_mass_flux_sink = &production_sink;
  h::HydroShadowReference shadow({.enabled = true, .cadence = 1U, .max_patches_per_rank = 1U,
                                  .max_storage_cells_per_rank = 1U << 16U,
                                  .max_workspace_bytes_per_rank = 1ULL << 30U});
  shadow.beginEpoch(0U);
  const auto spec = captureSpec(f, source_context, remote_ids, local_rows, promise);
  shadow.setGlobalCandidates({spec.identity});
  auto sample = shadow.capture(spec);
  assert(sample);
  shadow.observeProduction(std::move(sample), f.after);
  // The synchronized runtime finalizes endpoints with the hold active.
  shadow.finalizePendingEndpoints(false);
  Outcome outcome;
  outcome.completed_before_table = shadow.report().completed;
  outcome.unpublished_before_table = shadow.unpublishedResultCount();
  if (table != Table::kNone) {
    std::vector<std::uint64_t> ids;
    std::vector<double> fractions;
    donorTable(f, ids, fractions);
    const std::size_t lanes = f.production.gas_composition.laneCount();
    if (table == Table::kMalformed) fractions.resize(fractions.size() / 2U);
    shadow.provideRemoteDonorFractions(ids, fractions, lanes);
  }
  shadow.finalizePendingEndpoints();
  // Post-commit composition check against the step-start production rows
  // (any common snapshot: only equality between cases is asserted).
  h::HydroShadowPrivateRegion snapshot;
  snapshot.resize(k_cells);
  snapshot.cell_volume_comoving = f.geometry.cell_volume_comoving;
  snapshot.composition = f.production.gas_composition.emptyLike(k_cells, false);
  for (std::size_t r = 0; r < k_cells; ++r) {
    snapshot.gas_cell_ids[r] = f.production.gas_cells.gas_cell_id[r];
    snapshot.composition.unpackRow(r, f.production.gas_composition.packRow(r), snapshot.gas_cell_ids[r]);
  }
  shadow.compareCompositionAfterCommit(h::hydroShadowEffectiveKey(spec.identity), snapshot, update.time_end_code);
  for (const auto& result : shadow.regionResults()) {
    if (result.compared_composition_cells != 0U) outcome.composition = result.composition;
  }
  outcome.report = shadow.report();
  return outcome;
}

void checkNoexceptContainment() {
  h::HydroShadowReference shadow({.enabled = true, .cadence = 1U});
  const h::HydroShadowPatchCapture spec;
  static_assert(noexcept(shadow.capture(spec)));
  static_assert(noexcept(shadow.setGlobalCandidates({})));
  static_assert(noexcept(shadow.finalizePendingEndpoints()));
  static_assert(noexcept(shadow.provideRemoteDonorFractions({}, {}, 0U)));
  static_assert(noexcept(shadow.unpublishedResultCount()));
  // A spec without geometry is a counted refusal, never an exception.
  shadow.beginEpoch(0U);
  shadow.setGlobalCandidates({spec.identity});
  assert(!shadow.capture(spec));
  assert(shadow.report().skipped_geometry == 1U);
}

void checkRemoteDonorsHeldAndResolved() {
  const Fixture f = makeFixture();
  const Outcome resolved = run(f, true, false, true, Table::kValid);
  // Held: nothing compared before the production table arrived.
  assert(resolved.completed_before_table == 0U);
  assert(resolved.unpublished_before_table == 1U);
  assert(resolved.report.completed == 1U);
  assert(resolved.report.composition_remote_donor_slots_resolved == f.marked_slots.size());
  assert(resolved.report.composition_remote_donor_faces == 0U);
  // Production-rule reconciliation is small for a floor-free deck (the
  // initial composition mass equals the cell mass; only transport round-off
  // and receiver-fraction-free donors remain).
  assert(std::isfinite(resolved.report.composition_reconciled_mass_code));
  assert(std::abs(resolved.report.composition_reconciled_mass_code) <
         1.0e-6 * static_cast<double>(k_cells) * f.geometry.cell_volume_comoving);

  // Same donors through LOCAL rows: no hold, identical private composition.
  const Outcome local = run(f, false, true, true, Table::kNone);
  assert(local.completed_before_table == 1U);
  assert(local.report.composition_remote_donor_slots_resolved == 0U);
  assert(local.report.composition_remote_donor_faces == 0U);
  assert(local.composition.weighted_abs_difference == resolved.composition.weighted_abs_difference);
  assert(local.composition.linf_difference == resolved.composition.linf_difference);
}

void checkCountedFallbackWithoutTable() {
  const Fixture f = makeFixture();
  // No promise: compared immediately with the counted receiver fallback.
  const Outcome unpromised = run(f, true, false, false, Table::kNone);
  assert(unpromised.completed_before_table == 1U);
  assert(unpromised.report.composition_remote_donor_faces != 0U);
  assert(unpromised.report.composition_remote_donor_slots_resolved == 0U);
  // Promise broken by a malformed table: released, compared, fallback counted.
  const Outcome malformed = run(f, true, false, true, Table::kMalformed);
  assert(malformed.completed_before_table == 0U);
  assert(malformed.report.completed == 1U);
  assert(malformed.report.composition_remote_donor_faces != 0U);
  // Promise never fulfilled: the default (forced) finalize still compares.
  const Outcome unfulfilled = run(f, true, false, true, Table::kNone);
  assert(unfulfilled.report.completed == 1U);
  assert(unfulfilled.report.composition_remote_donor_faces != 0U);
}

void checkHeldWorkCountedOnRejection() {
  const Fixture f = makeFixture();
  const h::HydroUpdateContext update{.dt_code = 1.0e-3, .time_begin_code = 0.0, .time_end_code = 1.0e-3};
  NullSink production_sink;
  h::HydroSourceContext source_context{.update = update};
  source_context.face_mass_flux_sink = &production_sink;
  h::HydroShadowReference shadow({.enabled = true, .cadence = 1U, .max_patches_per_rank = 1U,
                                  .max_storage_cells_per_rank = 1U << 16U,
                                  .max_workspace_bytes_per_rank = 1ULL << 30U});
  shadow.beginEpoch(4U);
  const auto spec = captureSpec(f, source_context, true, false, true);
  shadow.setGlobalCandidates({spec.identity});
  auto sample = shadow.capture(spec);
  assert(sample);
  shadow.observeProduction(std::move(sample), f.after);
  shadow.finalizePendingEndpoints(false);
  assert(shadow.unpublishedResultCount() == 1U);
  // Rejected attempt: the retry opens the same step again.
  shadow.beginEpoch(4U);
  assert(shadow.unpublishedResultCount() == 0U);
  assert(shadow.report().completed == 0U);
  assert(shadow.regionResults().empty());
}

}  // namespace

int main() {
  checkNoexceptContainment();
  checkRemoteDonorsHeldAndResolved();
  checkCountedFallbackWithoutTable();
  checkHeldWorkCountedOnRejection();
  return 0;
}
