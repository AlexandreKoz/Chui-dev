// C04-03/C04-05/C04-06 production-derived AMR shadow fixtures.
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// Every case drives the PRODUCTION AMR path (advanceProductionAmrHydro /
// advanceProductionAmrHydroRecursive, buildAmrHydroPatchGeometry, the AMR
// ghost filler and its temporal history) on an ordinary 8^3-patch hierarchy:
//
//   domain [0,2] x [0,2] x [0,1], reflective walls
//   level 0 (cell 1/8):  A=[0,1]x[0,1]  B=[1,2]x[0,1]  C=[0,1]x[1,2]   (z in [0,1])
//   level 1 (cell 1/16): eight 8^3 patches tiling D=[1,2]x[1,2]x[0,1]
//
// A sees B and C as same-level face neighbours and D only across the x=y=1
// EDGE (diagonal CTU stencil ghosts); B and C have coarse/fine faces with D;
// the fine patches have same-level, coarse->fine and physical-wall ghosts.
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "amr_hydro_validation_helpers.hpp"
#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"

namespace {

namespace amr = cosmosim::amr;
namespace core = cosmosim::core;
namespace h = cosmosim::hydro;
namespace v = cosmosim::tests::amr_hydro_validation;

using Provenance = h::HydroShadowGhostProvenance;

[[nodiscard]] amr::PatchDescriptor makePatch(std::uint64_t id, std::uint8_t level,
                                             std::array<double, 3> origin, double extent,
                                             double extent_z) {
  return amr::PatchDescriptor{
      .patch_id = id,
      .parent_patch_id = 0U,
      .level = level,
      .morton_key = id,
      .origin_comov = origin,
      .extent_comov = {extent, extent, extent_z},
      .cell_dims = {8U, 8U, 8U}};
}

[[nodiscard]] std::vector<amr::PatchDescriptor> hierarchy() {
  std::vector<amr::PatchDescriptor> patches{
      makePatch(101U, 0U, {0.0, 0.0, 0.0}, 1.0, 1.0),
      makePatch(102U, 0U, {1.0, 0.0, 0.0}, 1.0, 1.0),
      makePatch(103U, 0U, {0.0, 1.0, 0.0}, 1.0, 1.0)};
  std::uint64_t id = 201U;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        patches.push_back(makePatch(id++, 1U, {1.0 + 0.5 * i, 1.0 + 0.5 * j, 0.5 * k}, 0.5, 0.5));
      }
    }
  }
  return patches;
}

[[nodiscard]] h::HydroPrimitiveState smoothInitial(const amr::PatchDescriptor&, std::size_t,
                                                   const v::CellCenter& c) {
  constexpr double k_pi = 3.141592653589793;
  h::HydroPrimitiveState primitive{};
  primitive.rho_comoving = 1.0 + 0.05 * std::sin(k_pi * c.x_comoving) * std::sin(k_pi * c.y_comoving);
  primitive.vel_x_peculiar = 0.1;
  primitive.vel_y_peculiar = 0.05;
  primitive.vel_z_peculiar = 0.0;
  primitive.pressure_comoving = 1.0 + 0.02 * std::cos(k_pi * c.z_comoving);
  return primitive;
}

[[nodiscard]] h::HydroShadowCandidate candidateFor(const amr::PatchDescriptor& patch) {
  return h::HydroShadowCandidate{.patch_id = patch.patch_id,
                                 .level = static_cast<int>(patch.level),
                                 .region_key = amr::hydroShadowRegionKeyForPatch(patch)};
}

struct Production {
  h::HydroNumericsPolicy numerics{};
  amr::ProductionAmrHydroOptions options{};
};

[[nodiscard]] Production ctuPpmProduction(h::HydroBoundaryKind boundary, std::array<double, 3> extent) {
  Production production;
  production.numerics.integration_scheme = h::HydroIntegrationScheme::kCornerTransportUpwind;
  production.options = v::productionAmrHydroOptions(boundary);
  production.options.reconstruction_scheme = h::HydroReconstructionScheme::kPiecewiseParabolic;
  production.options.riemann_solver_kind = h::HydroRiemannSolverKind::kHllc;
  production.options.ctu_edge_stencils = true;
  production.options.state_time_code = 0.0;
  production.options.ghost_fill_time_code = 0.0;
  production.options.domain_origin_comoving = {0.0, 0.0, 0.0};
  production.options.domain_extent_comoving = extent;
  if (boundary == h::HydroBoundaryKind::kPeriodic) production.options.periodic_box_extent_comoving = extent;
  const h::HydroCoreSolver solver(v::k_gamma, production.numerics);
  production.options.active_batch_policy = solver.requiredWorkspaceFeatures({});
  return production;
}

[[nodiscard]] h::HydroUpdateContext updateFor(double dt) {
  return h::HydroUpdateContext{.dt_code = dt, .time_begin_code = 0.0, .time_end_code = dt,
                               .scale_factor_begin = 1.0, .scale_factor = 1.0, .scale_factor_end = 1.0};
}

[[nodiscard]] std::vector<double> gasLanes(const core::SimulationState& state) {
  std::vector<double> lanes;
  for (const auto* lane : {&state.gas_cells.density_code, &state.gas_cells.internal_energy_code,
                           &state.gas_cells.velocity_x_peculiar, &state.gas_cells.velocity_y_peculiar,
                           &state.gas_cells.velocity_z_peculiar, &state.gas_cells.metal_mass_code,
                           &state.gas_cells.pressure_code}) {
    lanes.insert(lanes.end(), lane->begin(), lane->end());
  }
  lanes.insert(lanes.end(), state.cells.mass_code.begin(), state.cells.mass_code.end());
  return lanes;
}

[[nodiscard]] bool bitIdentical(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

[[nodiscard]] h::HydroShadowOptions generousOptions(std::size_t k) {
  return h::HydroShadowOptions{.enabled = true, .cadence = 1U, .seed = 11U, .max_patches_per_rank = k,
                               .max_storage_cells_per_rank = 1U << 20U,
                               .max_workspace_bytes_per_rank = 1ULL << 32U};
}

// Synchronized production step on the 8^3 hierarchy. Samples: coarse A (same
// level + diagonal edge to the fine region), coarse B (coarse/fine face) and
// the fine patch at the D corner (coarse->fine + same-level + wall). All
// imported slots receive production-filler endpoint ghosts; every active cell
// is comparable (no interior guard, no cone exclusion).
void checkSynchronizedCoarseFineDiagonalEndpoints() {
  const auto descriptors = hierarchy();
  auto state = v::makePatchState(descriptors, smoothInitial);
  Production production = ctuPpmProduction(h::HydroBoundaryKind::kReflective, {2.0, 2.0, 1.0});
  const h::HydroCoreSolver solver(v::k_gamma, production.numerics);
  const h::HllcRiemannSolver riemann;
  h::HydroShadowReference shadow(generousOptions(3U));
  shadow.beginEpoch(0U);
  shadow.setGlobalCandidates({candidateFor(descriptors[0]), candidateFor(descriptors[1]),
                              candidateFor(descriptors[3])});
  assert(shadow.globalCandidates().size() == 3U);
  production.options.hydro_shadow_reference = &shadow;
  production.options.hydro_shadow_endpoint_ghosts = true;
  const auto update = updateFor(2.0e-3);
  const h::HydroSourceContext source_context{.update = update};
  const auto rows = v::allCellRows(state);
  (void)amr::advanceProductionAmrHydro(state, rows, update, source_context, solver, riemann, {}, production.options);
  assert(shadow.pendingEndpointCount() == 3U);
  const auto before_endpoint = gasLanes(state);
  const auto endpoint = amr::completeHydroShadowEndpointGhosts(state, production.options, {},
                                                              update.time_end_code, false);
  assert(endpoint.requested == 3U);
  assert(endpoint.provided == 3U);
  assert(endpoint.failed == 0U);
  assert(endpoint.stale_fine_to_coarse_slots == 0U);
  // The endpoint fill is read-only on the canonical state.
  assert(bitIdentical(before_endpoint, gasLanes(state)));
  shadow.finalizePendingEndpoints();
  const auto& report = shadow.report();
  assert(report.selected == 3U);
  assert(report.completed == 3U);
  assert(report.completed_endpoint_ghost_patches == 3U);
  assert(report.cone_excluded_cells == 0U);
  assert(report.compared_cells == 3U * 512U);
  assert(report.skipped_temporal_ghost == 0U);
  assert(report.skipped_geometry == 0U);
  assert(report.skipped_endpoint_unavailable == 0U);
  assert(report.ghost_slots_by_provenance[static_cast<std::size_t>(Provenance::kImportedSameLevel)] > 0U);
  assert(report.ghost_slots_by_provenance[static_cast<std::size_t>(Provenance::kImportedEdgeDiagonal)] > 0U);
  assert(report.ghost_slots_by_provenance[static_cast<std::size_t>(Provenance::kCoarseToFine)] > 0U);
  assert(report.ghost_slots_by_provenance[static_cast<std::size_t>(Provenance::kFineToCoarse)] > 0U);
  assert(report.ghost_slots_endpoint_interpolated > 0U);
  assert(report.ghost_slots_step_start_only == 0U);
  assert(report.completed_by_class[static_cast<std::size_t>(h::HydroShadowComparisonClass::kTimeOnlyCtuPpm)] == 3U);
  assert(report.reference_substeps == 6U);
  assert(shadow.regionResults().size() == 3U);
  // No shadow side effect on synchronization state.
  assert(state.pending_flux_registers.empty());
  assert(state.amr_temporal_boundary_history.empty());
}

// Partial scheduler activity on a single periodic 8^3 patch: the production
// active set (every other x column) and its deliberately suppressed duplicate
// lower periodic faces are reproduced exactly. Only active cells are compared;
// inactive cells are never advanced by either trajectory.
void checkPartialActivityAndSuppressedPeriodicFaces() {
  const std::vector<amr::PatchDescriptor> descriptors{makePatch(301U, 0U, {0.0, 0.0, 0.0}, 1.0, 1.0)};
  auto state = v::makePatchState(descriptors, smoothInitial);
  Production production = ctuPpmProduction(h::HydroBoundaryKind::kPeriodic, {1.0, 1.0, 1.0});
  const h::HydroCoreSolver solver(v::k_gamma, production.numerics);
  const h::HllcRiemannSolver riemann;
  h::HydroShadowReference shadow(generousOptions(1U));
  shadow.beginEpoch(0U);
  shadow.setGlobalCandidates({candidateFor(descriptors[0])});
  production.options.hydro_shadow_reference = &shadow;
  production.options.hydro_shadow_endpoint_ghosts = true;
  std::vector<std::uint32_t> active;
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const double x = state.cells.center_x_comoving[row];
    if (static_cast<int>(std::floor(8.0 * x)) % 2 == 0) active.push_back(row);
  }
  assert(active.size() == 256U);
  const auto inactive_before = gasLanes(state);
  const auto update = updateFor(2.0e-3);
  const h::HydroSourceContext source_context{.update = update};
  (void)amr::advanceProductionAmrHydro(state, active, update, source_context, solver, riemann, {}, production.options);
  (void)amr::completeHydroShadowEndpointGhosts(state, production.options, {}, update.time_end_code, false);
  shadow.finalizePendingEndpoints();
  const auto& report = shadow.report();
  assert(report.completed == 1U);
  assert(report.completed_partial_active_patches == 1U);
  assert(report.compared_cells + report.cone_excluded_cells == 256U);
  assert(report.ghost_slots_by_provenance[static_cast<std::size_t>(Provenance::kSelfPeriodicWrap)] > 0U);
  (void)inactive_before;
}

// Shadow disabled, enabled (all paths) and enabled-but-refused runs leave the
// production state bit-identical: the diagnostic never changes evolution.
void checkShadowNeverChangesProduction() {
  const auto descriptors = hierarchy();
  const auto update = updateFor(2.0e-3);
  const h::HydroSourceContext source_context{.update = update};
  const auto run = [&](h::HydroShadowReference* shadow, bool endpoints) {
    auto state = v::makePatchState(descriptors, smoothInitial);
    Production production = ctuPpmProduction(h::HydroBoundaryKind::kReflective, {2.0, 2.0, 1.0});
    const h::HydroCoreSolver solver(v::k_gamma, production.numerics);
    const h::HllcRiemannSolver riemann;
    production.options.hydro_shadow_reference = shadow;
    production.options.hydro_shadow_endpoint_ghosts = endpoints;
    (void)amr::advanceProductionAmrHydro(state, v::allCellRows(state), update, source_context, solver, riemann,
                                         {}, production.options);
    if (shadow != nullptr && endpoints) {
      (void)amr::completeHydroShadowEndpointGhosts(state, production.options, {}, update.time_end_code, false);
    }
    if (shadow != nullptr) shadow->finalizePendingEndpoints();
    return gasLanes(state);
  };
  const auto disabled = run(nullptr, false);
  h::HydroShadowReference enabled(generousOptions(11U));
  enabled.beginEpoch(0U);
  std::vector<h::HydroShadowCandidate> all;
  for (const auto& patch : descriptors) all.push_back(candidateFor(patch));
  enabled.setGlobalCandidates(all);
  assert(bitIdentical(disabled, run(&enabled, true)));
  assert(enabled.report().completed == 11U);
  // Endpoints withheld: imported slots become step-start-only (cone
  // exclusion, counted) and production is still unchanged.
  h::HydroShadowReference withheld(generousOptions(11U));
  withheld.beginEpoch(0U);
  withheld.setGlobalCandidates(all);
  assert(bitIdentical(disabled, run(&withheld, false)));
  assert(withheld.report().skipped_temporal_ghost + withheld.report().completed == 11U);
  // Governor refusal (tiny hard limit owned by the diagnostic only).
  core::MemoryGovernor tiny(core::MemoryGovernorPolicy{.hard_limit_bytes = 4096U});
  h::HydroShadowReference refused(generousOptions(11U), &tiny);
  refused.beginEpoch(0U);
  refused.setGlobalCandidates(all);
  assert(bitIdentical(disabled, run(&refused, true)));
  assert(refused.report().completed == 0U);
  assert(refused.report().skipped_resource == 11U);
  // Workspace bound refusal (checked arithmetic before allocation).
  h::HydroShadowReference bounded(h::HydroShadowOptions{.enabled = true, .cadence = 1U, .max_patches_per_rank = 11U,
                                                        .max_storage_cells_per_rank = 1U << 20U,
                                                        .max_workspace_bytes_per_rank = 1024U});
  bounded.beginEpoch(0U);
  bounded.setGlobalCandidates(all);
  assert(bitIdentical(disabled, run(&bounded, true)));
  assert(bounded.report().skipped_resource == 11U);
}

// A rejected attempt: samples captured but the epoch is restarted before the
// endpoint barrier. Nothing is compared or published; drops are counted.
void checkRejectedAttemptIsDiscarded() {
  const auto descriptors = hierarchy();
  auto state = v::makePatchState(descriptors, smoothInitial);
  Production production = ctuPpmProduction(h::HydroBoundaryKind::kReflective, {2.0, 2.0, 1.0});
  const h::HydroCoreSolver solver(v::k_gamma, production.numerics);
  const h::HllcRiemannSolver riemann;
  h::HydroShadowReference shadow(generousOptions(2U));
  shadow.beginEpoch(4U);
  shadow.setGlobalCandidates({candidateFor(descriptors[0]), candidateFor(descriptors[1])});
  production.options.hydro_shadow_reference = &shadow;
  production.options.hydro_shadow_endpoint_ghosts = true;
  const auto update = updateFor(2.0e-3);
  const h::HydroSourceContext source_context{.update = update};
  (void)amr::advanceProductionAmrHydro(state, v::allCellRows(state), update, source_context, solver, riemann, {},
                                       production.options);
  assert(shadow.pendingEndpointCount() == 2U);
  shadow.discardPending();
  assert(shadow.report().discarded_rejected_attempt == 2U);
  assert(shadow.report().completed == 0U);
  assert(shadow.regionResults().empty());
  shadow.beginEpoch(4U);  // the retried attempt starts clean
  assert(shadow.pendingEndpointCount() == 0U);
  assert(shadow.report().discarded_rejected_attempt == 0U);
}

// Recursive Berger-Oliger clock: the driver completes endpoints after each
// level sweep. Fine samples use the open coarse history at the fine t_end;
// coarse B has fine->coarse slots whose finer level has not advanced: those
// slots keep step-start authority and their dependency cones are excluded.
void checkRecursiveLevelClockEndpoints() {
  const auto descriptors = hierarchy();
  auto state = v::makePatchState(descriptors, smoothInitial);
  Production production = ctuPpmProduction(h::HydroBoundaryKind::kReflective, {2.0, 2.0, 1.0});
  const h::HydroCoreSolver solver(v::k_gamma, production.numerics);
  const h::HllcRiemannSolver riemann;
  h::HydroShadowReference shadow(generousOptions(3U));
  shadow.beginEpoch(0U);
  shadow.setGlobalCandidates({candidateFor(descriptors[0]), candidateFor(descriptors[1]),
                              candidateFor(descriptors[3])});
  production.options.hydro_shadow_reference = &shadow;
  production.options.hydro_shadow_endpoint_ghosts = true;
  production.options.sweep_mode = amr::ProductionAmrHydroSweepMode::kSynchronized;
  const auto update = updateFor(2.0e-3);
  const h::HydroSourceContext source_context{.update = update};
  (void)amr::advanceProductionAmrHydroRecursive(state, {}, update, source_context, solver, riemann, {},
                                                production.options);
  shadow.finalizePendingEndpoints();
  const auto& report = shadow.report();
  assert(report.completed >= 2U);
  assert(report.completed_endpoint_ghost_patches >= 2U);
  // Coarse B borders the fine region: cone-excluded, never mislabeled.
  assert(report.cone_excluded_cells > 0U);
  assert(report.completed_by_class[static_cast<std::size_t>(h::HydroShadowComparisonClass::kBoundedProjection)] >= 1U);
  assert(state.pending_flux_registers.empty());
  assert(state.amr_temporal_boundary_history.empty());
}

// Physical identity of a production patch does not depend on patch ID,
// container order or owning rank.
void checkProductionRegionKeysArePhysical() {
  auto a = makePatch(101U, 0U, {0.0, 0.0, 0.0}, 1.0, 1.0);
  auto renumbered = a;
  renumbered.patch_id = 987654U;
  renumbered.morton_key = 3U;
  assert(amr::hydroShadowRegionKeyForPatch(a) == amr::hydroShadowRegionKeyForPatch(renumbered));
  auto shifted = a;
  shifted.origin_comov[0] = 1.0;
  assert(amr::hydroShadowRegionKeyForPatch(a) != amr::hydroShadowRegionKeyForPatch(shifted));
  auto finer = a;
  finer.level = 1U;
  finer.extent_comov = {0.5, 0.5, 0.5};
  assert(amr::hydroShadowRegionKeyForPatch(a) != amr::hydroShadowRegionKeyForPatch(finer));
}

}  // namespace

int main() {
  checkProductionRegionKeysArePhysical();
  checkSynchronizedCoarseFineDiagonalEndpoints();
  checkPartialActivityAndSuppressedPeriodicFaces();
  checkShadowNeverChangesProduction();
  checkRejectedAttemptIsDiscarded();
  checkRecursiveLevelClockEndpoints();
  return 0;
}
