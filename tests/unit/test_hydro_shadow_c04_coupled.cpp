// C04 source fixtures. NOT COMPILED OR RUN in this handoff.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <memory>
#include <vector>

#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"
#include "cosmosim/physics/cooling_heating.hpp"

namespace {
namespace h = cosmosim::hydro;
namespace p = cosmosim::physics;

void checkGlobalTopKIndependentOfRankAndOrder() {
  std::vector<h::HydroShadowCandidate> physical;
  for (std::uint64_t id = 100U; id < 180U; ++id) {
    physical.push_back({id, static_cast<int>(id % 4U)});
  }
  const auto expected = h::rankHydroShadowCandidates(physical, 3U, 37U, 7987U);
  for (unsigned ranks : {1U, 2U, 3U, 7U, 16U}) {
    std::vector<h::HydroShadowCandidate> union_top;
    for (unsigned rank = 0; rank < ranks; ++rank) {
      std::vector<h::HydroShadowCandidate> local;
      for (const auto& candidate : physical) {
        // Change the ownership mapping: physical identities do not change.
        if ((candidate.patch_id * 23U + ranks) % ranks == rank) local.push_back(candidate);
      }
      std::reverse(local.begin(), local.end());
      const auto top = h::rankHydroShadowCandidates(local, 3U, 37U, 7987U);
      union_top.insert(union_top.end(), top.begin(), top.end());
    }
    const auto gathered = h::rankHydroShadowCandidates(union_top, 3U, 37U, 7987U);
    assert(gathered.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      assert(gathered[i].patch_id == expected[i].patch_id);
      assert(gathered[i].level == expected[i].level);
    }
  }
}

struct TestSource final : h::HydroSourceTerm {
  explicit TestSource(std::uint64_t* calls = nullptr) : production_calls(calls) {}
  std::uint64_t* production_calls = nullptr;
  [[nodiscard]] std::unique_ptr<h::HydroSourceTerm> cloneForShadow() const override {
    // The private trajectory has its own counter and never writes the source
    // object's production-side diagnostic/ledger.
    return std::make_unique<TestSource>();
  }
  [[nodiscard]] h::HydroConservedState sourceForCell(
      std::size_t, const h::HydroConservedState&,
      const h::HydroPrimitiveState&, const h::HydroSourceContext&) const override {
    if (production_calls != nullptr) ++*production_calls;
    h::HydroConservedState source{};
    source.total_energy_density_comoving = 1.0e-5;
    return source;
  }
};

class TestClosure final : public h::HydroThermodynamicClosure {
 public:
  [[nodiscard]] h::HydroThermodynamicClosureResult evaluate(
      std::size_t, const h::HydroConservedState&,
      const h::HydroPrimitiveState& ideal,
      double, double) const override {
    return {.pressure_comoving = ideal.pressure_comoving,
            .signal_speed_squared_code = ideal.signal_speed_squared_code,
            .target_specific_internal_energy_code = ideal.specific_internal_energy_code,
            .uses_effective_ism = true, .valid = true};
  }
};

h::HydroPatchGeometry makeClosedPatch(h::HydroBoundaryKind boundary) {
  auto geometry = h::makeCartesianPatchGeometry({
      .nx = 8U, .ny = 8U, .nz = 8U,
      .cell_width_x_comoving = 0.125,
      .cell_width_y_comoving = 0.125,
      .cell_width_z_comoving = 0.125});
  h::appendCartesianBoundaryGhostFaces(geometry, boundary);
  (void)h::completeReconstructionStencils(geometry);
  if (boundary == h::HydroBoundaryKind::kPeriodic) {
    (void)h::completePeriodicReconstructionStencils(geometry);
  }
  return geometry;
}

void checkSmallClosedPatchWithActiveSourcesAndCosmology(h::HydroBoundaryKind boundary) {
  constexpr double gamma = 5.0 / 3.0;
  auto geometry = makeClosedPatch(boundary);
  h::HydroConservedStateSoa state(geometry.totalCellStorageCount());
  state.enableAuxiliaryInternalEnergy(true);
  for (std::size_t i = 0; i < geometry.cellCount(); ++i) {
    const auto coords = geometry.cellIjk(i);
    h::HydroPrimitiveState primitive{};
    primitive.rho_comoving = 1.0 + 0.001 * static_cast<double>(coords[0] + coords[1]);
    primitive.vel_x_peculiar = 0.012;
    primitive.vel_y_peculiar = 0.008;
    primitive.vel_z_peculiar = -0.002;
    primitive.pressure_comoving = 1.0;
    primitive.metallicity_mass_fraction = 0.01;
    state.storeCell(i, h::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
    state.auxiliaryInternalEnergyDensityComoving()[i] = primitive.pressure_comoving / (gamma - 1.0);
  }
  h::fillHydroBoundaryGhostCells(state, geometry, gamma);
  h::fillHydroBoundaryGhostAuxiliaryInternalEnergy(state, geometry);
  std::vector<std::size_t> cells(geometry.cellCount());
  std::vector<std::size_t> faces(geometry.faces.size());
  std::iota(cells.begin(), cells.end(), 0U);
  std::iota(faces.begin(), faces.end(), 0U);
  const h::HydroActiveSetView active{.active_cells = cells, .active_faces = faces};
  const h::HydroUpdateContext update{
      .dt_code = 1.0e-5, .time_begin_code = 0.0, .time_end_code = 1.0e-5,
      .scale_factor_begin = 0.7, .scale_factor = 0.70001, .scale_factor_end = 0.70002,
      .hubble_rate_begin_code = 0.1, .hubble_rate_code = 0.09999,
      .hubble_rate_end_code = 0.09998, .comoving_coordinates = true};
  const auto courant = h::hydroTransportCourantFactors(update.dt_code, update.scale_factor,
      geometry.cell_width_x_comoving, geometry.cell_width_y_comoving, geometry.cell_width_z_comoving);
  const h::HydroReconstructionPolicy reconstruction_policy{
      .limiter = h::HydroSlopeLimiter::kVanLeer,
      .dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
      .rho_floor = 1.0e-9, .pressure_floor = 1.0e-9,
      .enable_muscl_hancock_predictor = true, .adiabatic_index = gamma};
  h::HydroNumericsPolicy numerics;
  numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  numerics.dual_energy.mode = h::HydroDualEnergyMode::kThermalFraction;
  numerics.dual_energy.thermal_fraction_switch = 0.003;
  h::HydroCoreSolver solver(gamma, numerics);
  h::PpmReconstruction ppm(reconstruction_policy);
  h::HllcRiemannSolver riemann;
  h::HydroScratchBuffers scratch;
  std::uint64_t physical_calls = 0;
  TestSource source(&physical_calls);
  const std::vector<const h::HydroSourceTerm*> sources{&source};
  TestClosure closure;
  h::HydroSourceContext source_context{};
  source_context.update = update;
  source_context.thermodynamic_closure = &closure;
  h::HydroShadowReference monitor({.enabled = true, .cadence = 1,
      .max_patches_per_rank = 1, .max_storage_cells_per_rank = 4096,
      .max_workspace_bytes_per_rank = 64U * 1024U * 1024U});
  monitor.beginEpoch(5);
  monitor.setGlobalCandidates({{88U, 1}});
  auto sample = monitor.capture(88U, 1, state, geometry, active, update, numerics,
      h::HydroReconstructionScheme::kPiecewiseParabolic, reconstruction_policy,
      h::HydroRiemannSolverKind::kHllc, sources, source_context);
  assert(sample);
  solver.advancePatchActiveSetBatchedWithScratch(state, geometry, active, update, ppm,
      riemann, sources, source_context, scratch, nullptr,
      solver.requiredWorkspaceFeatures({}));
  const auto physical_count = physical_calls;
  const auto produced = state;
  monitor.compare(std::move(sample), state, geometry, active, update, numerics,
                  h::HydroReconstructionScheme::kPiecewiseParabolic);
  assert(physical_calls == physical_count); // shadow did not touch production source
  assert(monitor.report().completed == 1U);
  assert(monitor.report().compared_cells == 512U);
  assert(monitor.report().compared_auxiliary_cells == 512U);
  assert(monitor.report().completed_closed_boundary_patches == 1U);
  assert(monitor.report().completed_source_active_patches == 1U);
  assert(monitor.report().completed_evolving_comoving_patches == 1U);
  assert(monitor.report().reference_substeps == 2U);
  for (std::size_t i = 0; i < state.size(); ++i) {
    assert(state.loadCell(i).total_energy_density_comoving == produced.loadCell(i).total_energy_density_comoving);
  }
}

h::HydroConservedStateSoa uniformState(const h::HydroPatchGeometry& geometry, double gamma) {
  h::HydroConservedStateSoa state(geometry.totalCellStorageCount());
  for (std::size_t row = 0; row < state.size(); ++row) {
    h::HydroPrimitiveState primitive{};
    primitive.rho_comoving = 1.0 + 1.0e-3 * static_cast<double>(row % 7U);
    primitive.vel_x_peculiar = 0.01;
    primitive.pressure_comoving = 1.0;
    primitive.metallicity_mass_fraction = 0.01;
    state.storeCell(row, h::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
  return state;
}

// C04-03 (rev 6): an imported ghost without a promised endpoint is no longer a
// blanket refusal. Its slot keeps step-start authority and every real cell in
// its halfstep dependency cone is excluded (bounded projection, counted); a
// patch whose EVERY cell lies in such cones is refused as skipped_temporal_ghost.
void checkImportedGhostWithoutEndpointIsConeExcluded() {
  constexpr double gamma = 5.0 / 3.0;
  auto geometry = makeClosedPatch(h::HydroBoundaryKind::kReflective);
  geometry.ghost_cells.front().mutation_rights = h::HydroGhostMutationRights::kReadOnlyImported;
  geometry.ghost_cells.front().boundary_kind = h::HydroBoundaryKind::kImportedMpi;
  auto before = uniformState(geometry, gamma);
  std::vector<std::size_t> cells(geometry.cellCount());
  std::vector<std::size_t> faces(geometry.faces.size());
  std::iota(cells.begin(), cells.end(), 0U);
  std::iota(faces.begin(), faces.end(), 0U);
  h::HydroNumericsPolicy numerics;
  numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  h::HydroUpdateContext update{.dt_code = 1.0e-5, .time_begin_code = 0.0, .time_end_code = 1.0e-5};
  h::HydroShadowReference monitor({.enabled = true, .cadence = 1});
  monitor.beginEpoch(0);
  monitor.setGlobalCandidates({{27U, 1}});
  const h::HydroReconstructionPolicy policy{.adiabatic_index = gamma};
  auto sample = monitor.capture(27U, 1, before, geometry,
      {.active_cells = cells, .active_faces = faces}, update,
      numerics, h::HydroReconstructionScheme::kPiecewiseParabolic, policy,
      h::HydroRiemannSolverKind::kHllc, {}, h::HydroSourceContext{});
  assert(sample);
  const auto produced = before;
  monitor.compare(std::move(sample), produced, geometry, {.active_cells = cells, .active_faces = faces},
                  update, numerics, h::HydroReconstructionScheme::kPiecewiseParabolic);
  assert(monitor.report().completed == 1U);
  assert(monitor.report().cone_excluded_cells > 0U);
  assert(monitor.report().compared_cells + monitor.report().cone_excluded_cells == 512U);
  assert(monitor.report().ghost_slots_step_start_only == 1U);
  assert(monitor.report().completed_by_class[static_cast<std::size_t>(
             h::HydroShadowComparisonClass::kBoundedProjection)] == 1U);

  // Every slot imported, no endpoint: every cell of 8^3 lies in a cone.
  auto imported = makeClosedPatch(h::HydroBoundaryKind::kReflective);
  for (auto& ghost : imported.ghost_cells) {
    ghost.mutation_rights = h::HydroGhostMutationRights::kReadOnlyImported;
    ghost.boundary_kind = h::HydroBoundaryKind::kImportedMpi;
  }
  auto state = uniformState(imported, gamma);
  h::HydroShadowReference refusing({.enabled = true, .cadence = 1});
  refusing.beginEpoch(0);
  refusing.setGlobalCandidates({{28U, 1}});
  assert(!refusing.capture(28U, 1, state, imported, {.active_cells = cells, .active_faces = faces}, update,
                           numerics, h::HydroReconstructionScheme::kPiecewiseParabolic, policy,
                           h::HydroRiemannSolverKind::kHllc, {}, h::HydroSourceContext{}));
  assert(refusing.report().skipped_temporal_ghost == 1U);
}

void checkCoolingReadOnlyProxyIsIndependent() {
  // The production cooling rate table stays read-only; the shadow must not
  // allocate a second persistent table or own the production event ledger.
  p::CoolingModelConfig config;
  p::CoolingHeatingSource source{p::CoolingRateProvider(config),
                                 p::CoolingSourceIntegrator(1.0e-8)};
  const auto cloned = source.cloneForShadow();
  assert(cloned);
}
} // namespace

int main() {
  checkGlobalTopKIndependentOfRankAndOrder();
  checkSmallClosedPatchWithActiveSourcesAndCosmology(h::HydroBoundaryKind::kPeriodic);
  checkSmallClosedPatchWithActiveSourcesAndCosmology(h::HydroBoundaryKind::kReflective);
  checkImportedGhostWithoutEndpointIsConeExcluded();
  checkCoolingReadOnlyProxyIsIndependent();
  return 0;
}
