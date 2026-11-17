// C04 source fixture — NOT RUN (user owns compilation and qualification).
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"

namespace {
namespace h = cosmosim::hydro;

void checkShadowOneStep(h::HydroIntegrationScheme production_mode) {
  constexpr double gamma = 1.4;
  h::HydroPatchGeometry geometry = h::makeCartesianPatchGeometry({
      .nx = 32, .ny = 1, .nz = 1,
      .cell_width_x_comoving = 1.0 / 32.0,
      .cell_width_y_comoving = 1.0,
      .cell_width_z_comoving = 1.0});
  h::HydroConservedStateSoa state(geometry.totalCellStorageCount());
  for (std::size_t row = 0; row < state.size(); ++row) {
    const double x = static_cast<double>(row) / 32.0;
    h::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.0 + 0.01 * std::sin(6.283185307179586 * x);
    primitive.vel_x_peculiar = 0.1;
    primitive.pressure_comoving = 1.0;
    primitive.metallicity_mass_fraction = 0.01;
    state.storeCell(row, h::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
  std::vector<std::size_t> cells(geometry.cellCount());
  std::iota(cells.begin(), cells.end(), 0U);
  std::vector<std::size_t> faces(geometry.faces.size());
  std::iota(faces.begin(), faces.end(), 0U);
  const h::HydroActiveSetView active{.active_cells = cells, .active_faces = faces};
  h::HydroUpdateContext update;
  update.dt_code = 1.0e-4;
  update.time_end_code = update.dt_code;
  h::HydroSourceContext source_context;
  source_context.update = update;
  h::HydroNumericsPolicy numerics;
  numerics.integration_scheme = production_mode;
  h::HydroShadowReference monitor({.enabled = true, .cadence = 1,
      .max_patches_per_rank = 1, .max_storage_cells_per_rank = 1024,
      .max_workspace_bytes_per_rank = 16U * 1024U * 1024U});
  monitor.beginEpoch(0);
  const auto reference_scheme = production_mode == h::HydroIntegrationScheme::kCornerTransportUpwind
          ? h::HydroReconstructionScheme::kPiecewiseParabolic
          : h::HydroReconstructionScheme::kPiecewiseLinear;
  monitor.setGlobalCandidates({{0U, 0}});
  const h::HydroReconstructionPolicy reference_input_policy{.adiabatic_index = gamma};
  auto sample = monitor.capture(0U, 0, state, geometry, active, update, numerics,
      reference_scheme, reference_input_policy, h::HydroRiemannSolverKind::kHllc, {}, h::HydroSourceContext{});
  assert(sample);
  h::HydroCoreSolver solver(gamma, numerics);
  const auto cfl = h::hydroTransportCourantFactors(update.dt_code, 1.0,
      geometry.cell_width_x_comoving, geometry.cell_width_y_comoving,
      geometry.cell_width_z_comoving);
  h::HydroReconstructionPolicy policy;
  policy.dt_over_dx_code = cfl[0];
  policy.dt_over_cell_width_code = cfl;
  policy.adiabatic_index = gamma;
  h::MusclHancockReconstruction plm(policy);
  h::PpmReconstruction ppm(policy);
  h::HllcRiemannSolver riemann;
  h::HydroScratchBuffers scratch;
  const auto batch = solver.requiredWorkspaceFeatures(h::HydroActiveBatchPolicy{});
  solver.advancePatchActiveSetBatchedWithScratch(state, geometry, active, update,
      production_mode == h::HydroIntegrationScheme::kCornerTransportUpwind
          ? static_cast<const h::HydroReconstruction&>(ppm)
          : static_cast<const h::HydroReconstruction&>(plm),
      riemann, {}, source_context, scratch, nullptr, batch);
  const auto physical_state = state;
  monitor.compare(std::move(sample), state, geometry, active, update, numerics,
      production_mode == h::HydroIntegrationScheme::kCornerTransportUpwind
          ? h::HydroReconstructionScheme::kPiecewiseParabolic
          : h::HydroReconstructionScheme::kPiecewiseLinear);
  assert(monitor.report().selected == 1U);
  assert(monitor.report().completed == 1U);
  assert(monitor.report().compared_cells > 0U);
  assert(monitor.report().reference_substeps ==
         (production_mode == h::HydroIntegrationScheme::kCornerTransportUpwind ? 2U : 1U));
  // In particular no observer invocation changes any canonical physics lane.
  for (std::size_t row = 0; row < state.size(); ++row) {
    assert(state.loadCell(row).mass_density_comoving == physical_state.loadCell(row).mass_density_comoving);
    assert(state.loadCell(row).total_energy_density_comoving ==
           physical_state.loadCell(row).total_energy_density_comoving);
  }
}

void checkSamplingBudgetsAndCadence() {
  h::HydroPatchGeometry geometry = h::makeCartesianPatchGeometry({
      .nx = 32, .ny = 1, .nz = 1, .cell_width_x_comoving = 1.0,
      .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  h::HydroConservedStateSoa state(geometry.totalCellStorageCount());
  std::vector<std::size_t> cells(geometry.cellCount());
  std::iota(cells.begin(), cells.end(), 0U);
  std::vector<std::size_t> faces(geometry.faces.size());
  std::iota(faces.begin(), faces.end(), 0U);
  const h::HydroActiveSetView active{.active_cells = cells, .active_faces = faces};
  h::HydroUpdateContext update;
  h::HydroNumericsPolicy numerics;
  h::HydroShadowReference limited({.enabled = true, .cadence = 2,
      .max_patches_per_rank = 1, .max_storage_cells_per_rank = 16,
      .max_workspace_bytes_per_rank = 1024});
  limited.beginEpoch(1);
  limited.setGlobalCandidates({{0U, 0}});
  assert(!limited.capture(0, 0, state, geometry, active, update, numerics, h::HydroReconstructionScheme::kPiecewiseLinear, {}, h::HydroRiemannSolverKind::kHllc, {}, h::HydroSourceContext{}));
  assert(limited.report().skipped_cadence == 1U);
  limited.beginEpoch(0);
  limited.setGlobalCandidates({{0U, 0}});
  assert(!limited.capture(0, 0, state, geometry, active, update, numerics, h::HydroReconstructionScheme::kPiecewiseLinear, {}, h::HydroRiemannSolverKind::kHllc, {}, h::HydroSourceContext{}));
  assert(limited.report().skipped_resource == 1U);
  // Budget refusal and zero candidates leave production input untouched.
  assert(state.size() == geometry.totalCellStorageCount());
}

void checkSkippedWithoutStateMutation() {
  h::HydroPatchGeometry geometry = h::makeCartesianPatchGeometry({
      .nx = 8, .ny = 1, .nz = 1, .cell_width_x_comoving = 1.0,
      .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  h::HydroConservedStateSoa state(geometry.totalCellStorageCount());
  std::vector<std::size_t> cells(geometry.cellCount());
  std::iota(cells.begin(), cells.end(), 0U);
  std::vector<std::size_t> faces(geometry.faces.size());
  std::iota(faces.begin(), faces.end(), 0U);
  const h::HydroActiveSetView active{.active_cells = cells, .active_faces = faces};
  h::HydroShadowReference monitor({.enabled = true, .cadence = 1});
  monitor.beginEpoch(0);
  monitor.setGlobalCandidates({{0U, 0}});
  auto sample = monitor.capture(0U, 0, state, geometry, active,
                                h::HydroUpdateContext{}, h::HydroNumericsPolicy{},
                                h::HydroReconstructionScheme::kPiecewiseLinear, {}, h::HydroRiemannSolverKind::kHllc, {}, h::HydroSourceContext{});
  assert(!sample && monitor.report().skipped_geometry == 1U);
  h::HydroShadowReference disabled({.enabled = false});
  disabled.beginEpoch(0);
  assert(!disabled.capture(0U, 0, state, geometry, active,
                           h::HydroUpdateContext{}, h::HydroNumericsPolicy{},
                           h::HydroReconstructionScheme::kPiecewiseLinear, {}, h::HydroRiemannSolverKind::kHllc, {}, h::HydroSourceContext{}));
  assert(disabled.report().eligible == 0U);
}
}  // namespace

int main() {
  checkShadowOneStep(h::HydroIntegrationScheme::kMusclHancock);
  checkShadowOneStep(h::HydroIntegrationScheme::kCornerTransportUpwind);
  checkSkippedWithoutStateMutation();
  checkSamplingBudgetsAndCadence();
  return 0;
}
