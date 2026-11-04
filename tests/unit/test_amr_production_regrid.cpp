// CHUI Part 3 Checkpoint B — automatic production regrid and distributed
// composition transport.
// STATUS: written as source only. NOT RUN — USER-RESERVED.
//
// Invariants (exact up to roundoff unless stated):
//  1. Refine -> derefine of a patch carrying a composition sidecar (hybrid
//     species, two dust reservoirs, SGS energy) conserves mass, momentum,
//     total energy and EVERY composition lane total; each child inherits its
//     parent's lane mass fractions (piecewise-constant prolongation); the
//     derefined parent's lanes equal the sum of its children. Relative 1e-12.
//  2. The planner refines on the density/pressure jump indicator, derefines a
//     uniform octet, honours max_level, and blocks a refinement that would
//     break the 2:1 face balance with a (remote) coarser neighbour.
//  3. A regrid refused by memory admission leaves the hierarchy and the
//     composition sidecar bit-identical (no partial mutation).
//  4. Rank-striped id windows are disjoint and above the global maxima.
//  5. Distributed composition accumulator: a remote ghost donor is requested
//     by stable id, single-rank prepare() fails closed without its fractions,
//     inflow uses exactly the remote fractions, and a fine-side mirrored
//     transfer onto a REMOTE coarse cell becomes an outbound correction equal
//     and opposite to the local change (zero net lane change across ranks).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_production_regrid.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/composition_transport.hpp"
#include "../support/test_require.hpp"

namespace {

namespace amr = cosmosim::amr;
namespace core = cosmosim::core;
namespace hy = cosmosim::hydro;
namespace ph = cosmosim::physics;

constexpr double k_gamma = 5.0 / 3.0;

core::SimulationState makeRootState(bool jump) {
  core::SimulationState state;
  state.resizeCells(8U);
  state.resizePatches(1U);
  const amr::PatchDescriptor root{.patch_id = 1U, .parent_patch_id = 0U, .level = 0U, .morton_key = 1U,
                                  .origin_comov = {0.0, 0.0, 0.0}, .extent_comov = {1.0, 1.0, 1.0},
                                  .cell_dims = {2U, 2U, 2U}};
  amr::writePatchDescriptorToStateRow(state, 0U, root);
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] = 8U;
  state.patches.owning_rank[0] = 0;
  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t row = 0; row < 8U; ++row) {
    const std::size_t i = row % 2U;
    const std::size_t j = (row / 2U) % 2U;
    const std::size_t k = row / 4U;
    const double rho = jump && i == 1U ? 4.0 : 1.0;
    const double pressure = jump && i == 1U ? 2.0 : 1.0;
    state.cells.center_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(i);
    state.cells.center_y_comoving[row] = 0.25 + 0.5 * static_cast<double>(j);
    state.cells.center_z_comoving[row] = 0.25 + 0.5 * static_cast<double>(k);
    state.cells.patch_index[row] = 0U;
    state.cells.mass_code[row] = rho * 0.125;
    state.gas_cells.gas_cell_id[row] = 100U + row;
    state.gas_cells.density_code[row] = rho;
    state.gas_cells.pressure_code[row] = pressure;
    state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
    state.gas_cells.velocity_x_peculiar[row] = 0.1 * static_cast<double>(row);
    state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
    state.gas_cells.metal_mass_code[row] = 0.01 * (1.0 + static_cast<double>(row)) * state.cells.mass_code[row];
    records.push_back(core::GasCellIdentityRecord{.gas_cell_id = 100U + row, .parent_particle_id = std::nullopt,
                                                  .owning_patch_id = 1U, .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  ph::enableGasComposition(state, ph::CompositionLayout::make({
      .network = ph::CompositionNetworkKind::kHybridLight,
      .dust = ph::CompositionDustKind::kTwoReservoirV1,
      .sgs_energy = true}), {});
  return state;
}

amr::ProductionAmrHydroOptions regridOptions(const core::SimulationState& state, core::MemoryGovernor& governor) {
  governor.setBaselineOwnedBytes(core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state)));
  amr::ProductionAmrHydroOptions options;
  options.adiabatic_index = k_gamma;
  options.regrid_memory_governor = &governor;
  return options;
}

std::vector<double> laneTotals(const core::SimulationState& state) {
  std::vector<double> totals(state.gas_composition.laneCount(), 0.0);
  for (std::size_t lane = 0; lane < totals.size(); ++lane) {
    for (const double v : state.gas_composition.lane(lane)) totals[lane] += v;
  }
  return totals;
}

bool close(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max({1.0e-300, std::abs(a), std::abs(b)});
}

amr::ProductionAmrRegridPolicy jumpPolicy(amr::AmrRegridPassKind pass) {
  amr::ProductionAmrRegridPolicy policy;
  policy.criteria.use_mass_threshold = false;
  policy.criteria.use_jeans_condition = false;
  policy.criteria.use_particle_count = false;
  policy.criteria.use_gradient_indicator = true;
  policy.criteria.gradient_threshold = 0.5;
  policy.criteria.derefine_hysteresis = 0.5;
  policy.max_level = 3U;
  policy.pass = pass;
  return policy;
}

void testRefineDerefineConservesHydroAndComposition() {
  core::SimulationState state = makeRootState(true);
  const std::vector<double> lanes_before = laneTotals(state);
  // Give one row a distinct chemistry tag and dwell (restriction rule check).
  state.gas_composition.chemistry_model_tag[3] = 2U;
  state.gas_composition.model_dwell_years[3] = 5.0;
  core::MemoryGovernor governor;
  const auto plan = amr::planProductionAmrRegrid(state, {}, jumpPolicy(amr::AmrRegridPassKind::kRefineOnly));
  CHUI_REQUIRE(plan.operations.size() == 1U);
  CHUI_REQUIRE(plan.operations[0].kind == amr::AmrRegridOperationKind::kRefine);
  const auto window = amr::rankStripedAmrRegridIdWindow(1U, 107U, 0, 1, plan.required_patch_ids, plan.required_gas_cell_ids);
  const auto summary = amr::applyProductionAmrRegridPlan(state, plan, window, regridOptions(state, governor));
  CHUI_REQUIRE(summary.applied_refinements == 1U);
  CHUI_REQUIRE(state.patches.size() == 8U && state.cells.size() == 64U);
  CHUI_REQUIRE(std::abs(summary.mass_change_code) < 1.0e-13);
  CHUI_REQUIRE(std::abs(summary.momentum_x_change_code) < 1.0e-13);
  CHUI_REQUIRE(std::abs(summary.total_energy_change_code) < 1.0e-12);
  const std::vector<double> lanes_refined = laneTotals(state);
  for (std::size_t lane = 0; lane < lanes_before.size(); ++lane) {
    CHUI_REQUIRE(close(lanes_refined[lane], lanes_before[lane], 1.0e-12));
  }
  state.gas_composition.requireAlignedWith(state.gas_cells.gas_cell_id, "test refined alignment");

  // Uniformize so that every child allows derefinement, keeping lane totals.
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    state.gas_cells.density_code[row] = 2.0;
    state.gas_cells.pressure_code[row] = 1.0;
  }
  const auto deplan = amr::planProductionAmrRegrid(state, {}, jumpPolicy(amr::AmrRegridPassKind::kDerefineOnly));
  CHUI_REQUIRE(deplan.operations.size() == 1U);
  CHUI_REQUIRE(deplan.operations[0].kind == amr::AmrRegridOperationKind::kDerefine);
  CHUI_REQUIRE(deplan.operations[0].patch.patch_id == 1U);
  std::uint64_t max_cell = 0U;
  for (const auto id : state.gas_cells.gas_cell_id) max_cell = std::max(max_cell, id);
  const auto dewindow = amr::rankStripedAmrRegridIdWindow(9U, max_cell, 0, 1, 0U, deplan.required_gas_cell_ids);
  core::MemoryGovernor governor2;
  (void)amr::applyProductionAmrRegridPlan(state, deplan, dewindow, regridOptions(state, governor2));
  CHUI_REQUIRE(state.patches.size() == 1U && state.cells.size() == 8U);
  const std::vector<double> lanes_after = laneTotals(state);
  for (std::size_t lane = 0; lane < lanes_before.size(); ++lane) {
    CHUI_REQUIRE(close(lanes_after[lane], lanes_before[lane], 1.0e-12));
  }
  // Restriction metadata: most detailed tag survives, mu becomes unknown.
  bool saw_detailed = false;
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    saw_detailed = saw_detailed || state.gas_composition.chemistry_model_tag[row] == 2U;
    CHUI_REQUIRE(std::isnan(state.gas_composition.mean_molecular_weight[row]));
  }
  CHUI_REQUIRE(saw_detailed);
}

void testMemoryRefusalLeavesStateUntouched() {
  core::SimulationState state = makeRootState(true);
  const std::vector<double> lanes_before = laneTotals(state);
  const std::vector<std::uint64_t> ids_before(state.gas_cells.gas_cell_id.begin(), state.gas_cells.gas_cell_id.end());
  const auto plan = amr::planProductionAmrRegrid(state, {}, jumpPolicy(amr::AmrRegridPassKind::kRefineOnly));
  const std::uint64_t baseline = core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state));
  core::MemoryGovernor tight(core::MemoryGovernorPolicy{.hard_limit_bytes = baseline + 64U});
  const auto window = amr::rankStripedAmrRegridIdWindow(1U, 107U, 0, 1, plan.required_patch_ids, plan.required_gas_cell_ids);
  bool refused = false;
  try {
    (void)amr::applyProductionAmrRegridPlan(state, plan, window, regridOptions(state, tight));
  } catch (const core::MemoryAdmissionError&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
  CHUI_REQUIRE(state.cells.size() == 8U && state.patches.size() == 1U);
  CHUI_REQUIRE(std::equal(ids_before.begin(), ids_before.end(), state.gas_cells.gas_cell_id.begin()));
  const std::vector<double> lanes_after = laneTotals(state);
  for (std::size_t lane = 0; lane < lanes_before.size(); ++lane) {
    CHUI_REQUIRE(lanes_after[lane] == lanes_before[lane]);
  }
}

void testBalanceAndMaxLevel() {
  core::SimulationState state = makeRootState(true);
  // A remote coarser neighbour (level 0) next to a local level-1 patch blocks
  // refinement of the local patch (would create a 2-level jump).
  amr::writePatchDescriptorToStateRow(state, 0U, amr::PatchDescriptor{
      .patch_id = 1U, .parent_patch_id = 0U, .level = 1U, .morton_key = 1U,
      .origin_comov = {0.0, 0.0, 0.0}, .extent_comov = {1.0, 1.0, 1.0}, .cell_dims = {2U, 2U, 2U}});
  const amr::PatchDescriptor remote{.patch_id = 77U, .parent_patch_id = 0U, .level = 0U, .morton_key = 9U,
                                    .origin_comov = {1.0, 0.0, 0.0}, .extent_comov = {2.0, 2.0, 2.0},
                                    .cell_dims = {2U, 2U, 2U}};
  const auto blocked = amr::planProductionAmrRegrid(
      state, std::span<const amr::PatchDescriptor>(&remote, 1U), jumpPolicy(amr::AmrRegridPassKind::kRefineOnly));
  CHUI_REQUIRE(blocked.operations.empty() && blocked.blocked_by_balance == 1U);
  auto capped = jumpPolicy(amr::AmrRegridPassKind::kRefineOnly);
  capped.max_level = 1U;
  const auto at_max = amr::planProductionAmrRegrid(state, {}, capped);
  CHUI_REQUIRE(at_max.operations.empty() && at_max.blocked_by_max_level == 1U);
}

void testRankStripedWindowsAreDisjoint() {
  const auto w0 = amr::rankStripedAmrRegridIdWindow(100U, 5000U, 0, 4, 16U, 512U);
  const auto w3 = amr::rankStripedAmrRegridIdWindow(100U, 5000U, 3, 4, 16U, 512U);
  CHUI_REQUIRE(w0.first_patch_id == 101U && w0.first_gas_cell_id == 5001U);
  CHUI_REQUIRE(w3.first_patch_id == 101U + 3U * 16U && w3.first_gas_cell_id == 5001U + 3U * 512U);
}

// Minimal two-cell geometry with one imported ghost (slot 0) beyond cell 1.
hy::HydroPatchGeometry twoCellGeometry(hy::HydroFluxRegisterFaceRole role, std::uint64_t coarse_id) {
  hy::HydroPatchGeometry g;
  g.nx = 2U; g.ny = 1U; g.nz = 1U;
  g.cell_width_x_comoving = 0.5; g.cell_width_y_comoving = 1.0; g.cell_width_z_comoving = 1.0;
  g.cell_volume_comoving = 0.5;
  g.faces.push_back(hy::HydroFace{.owner_cell = 0, .neighbor_cell = 1, .area_comoving = 1.0});
  g.faces.push_back(hy::HydroFace{.owner_cell = 1, .neighbor_cell = 2, .ghost_cell_slot = 0, .area_comoving = 1.0});
  g.ghost_cells.push_back(hy::HydroGhostCell{.owner_real_cell = 1, .source_real_cell = 1, .ghost_cell = 2,
                                            .ghost_slot = 0, .boundary_kind = hy::HydroBoundaryKind::kImportedMpi});
  g.flux_register_faces.resize(2);
  g.flux_register_faces[1].role = role;
  g.flux_register_faces[1].coarse_gas_cell_id = coarse_id;
  return g;
}

void testDistributedCompositionAccumulator() {
  core::SimulationState state = makeRootState(false);
  const std::size_t lanes = state.gas_composition.laneCount();
  std::vector<std::uint32_t> dense{0U, 1U};
  const std::vector<std::uint64_t> remote_origin{9999U};  // not local
  // (a) inflow from a remote same-level ghost donor.
  {
    const auto g = twoCellGeometry(hy::HydroFluxRegisterFaceRole::kNone, 0U);
    ph::CompositionFluxAccumulator accumulator(state);
    accumulator.beginAdvance(g, dense, remote_origin);
    accumulator.recordFaceMassTransfer(1, 1, hy::k_invalid_cell_index, -1.0e-3);  // ghost -> cell 1
    accumulator.endAdvance();
    CHUI_REQUIRE(accumulator.requiredRemoteDonorIds() == std::vector<std::uint64_t>{9999U});
    bool refused = false;
    try { accumulator.prepare(state); } catch (const std::runtime_error&) { refused = true; }
    CHUI_REQUIRE(refused);
  }
  {
    const auto g = twoCellGeometry(hy::HydroFluxRegisterFaceRole::kNone, 0U);
    ph::CompositionFluxAccumulator accumulator(state);
    accumulator.beginAdvance(g, dense, remote_origin);
    accumulator.recordFaceMassTransfer(1, 1, hy::k_invalid_cell_index, -1.0e-3);
    accumulator.endAdvance();
    std::vector<double> fractions(lanes, 0.0);
    fractions[0] = 1.0;  // pure lane-0 donor: inflow must add only to lane 0
    const std::vector<std::uint64_t> ids{9999U};
    accumulator.setRemoteDonorFractions(ids, fractions, lanes);
    const double lane0_before = state.gas_composition.value(0U, 1U);
    const double lane1_before = state.gas_composition.value(1U, 1U);
    accumulator.prepareDeferred(state);
    CHUI_REQUIRE(accumulator.outboundRemoteCorrections().empty());
    accumulator.validatePrepared(state);
    state.cells.mass_code[1] += 1.0e-3;
    (void)accumulator.commit(state);
    // Reconciliation rescales by (m + dm)/(m + dm) == 1 here.
    CHUI_REQUIRE(close(state.gas_composition.value(0U, 1U), lane0_before + 1.0e-3, 1.0e-12));
    CHUI_REQUIRE(close(state.gas_composition.value(1U, 1U), lane1_before, 1.0e-12));
  }
  // (b) fine-side mirrored transfer onto a REMOTE coarse cell 8888.
  {
    core::SimulationState fresh = makeRootState(false);
    const auto g = twoCellGeometry(hy::HydroFluxRegisterFaceRole::kFine, 8888U);
    ph::CompositionFluxAccumulator accumulator(fresh);
    accumulator.beginAdvance(g, dense, remote_origin);
    const double dm = 2.0e-4;  // fine cell 1 -> remote coarse cell
    accumulator.recordFaceMassTransfer(1, 1, hy::k_invalid_cell_index, dm);
    accumulator.endAdvance();
    const std::vector<std::uint64_t> needed = accumulator.requiredRemoteDonorIds();
    CHUI_REQUIRE(std::find(needed.begin(), needed.end(), 8888U) != needed.end());
    accumulator.setRemoteDonorFractions({}, {}, lanes);  // donor is local (outflow)
    std::vector<double> local_before(lanes);
    for (std::size_t l = 0; l < lanes; ++l) local_before[l] = fresh.gas_composition.value(l, 1U);
    accumulator.prepareDeferred(fresh);
    const auto& outbound = accumulator.outboundRemoteCorrections();
    CHUI_REQUIRE(outbound.size() == 1U && outbound[0].gas_cell_id == 8888U);
    accumulator.validatePrepared(fresh);
    fresh.cells.mass_code[1] -= dm;
    (void)accumulator.commit(fresh);
    for (std::size_t l = 0; l < lanes; ++l) {
      const double local_change = fresh.gas_composition.value(l, 1U) - local_before[l];
      CHUI_REQUIRE(std::abs(local_change + outbound[0].lane_delta_code[l]) <=
                   1.0e-12 * std::max(1.0e-300, std::abs(outbound[0].lane_delta_code[l])) + 1.0e-18);
    }
  }
}

// R1 (A03): physics-aware indicators in the production planner.
//  * metallicity, species and vorticity interfaces are OPTIONAL requests;
//  * the cooling length is MANDATORY and blocks derefinement inside the
//    hysteresis band; a zero budget truncates and REPORTS a mandatory request;
//  * each indicator is inactive (no request) for uniform data.
amr::ProductionAmrRegridPolicy physicsPolicy(amr::AmrRegridPassKind pass) {
  amr::ProductionAmrRegridPolicy policy = jumpPolicy(pass);
  policy.criteria.use_gradient_indicator = false;
  return policy;
}

void makeUniform(core::SimulationState& state) {
  for (std::uint32_t row = 0; row < 8U; ++row) {
    state.gas_cells.metal_mass_code[row] = 0.01 * state.cells.mass_code[row];
    state.gas_cells.velocity_x_peculiar[row] = 0.0;
  }
}

void testPhysicsIndicators() {
  {  // metallicity interface (row-dependent Z in makeRootState)
    core::SimulationState state = makeRootState(false);
    auto policy = physicsPolicy(amr::AmrRegridPassKind::kRefineOnly);
    policy.physics.metallicity_jump = 0.5;
    const auto plan = amr::planProductionAmrRegrid(state, {}, policy);
    CHUI_REQUIRE(plan.operations.size() == 1U && plan.optional_refine_requests == 1U &&
                 plan.mandatory_refine_requests == 0U);
    makeUniform(state);
    CHUI_REQUIRE(amr::planProductionAmrRegrid(state, {}, policy).operations.empty());
  }
  {  // vorticity: shear v_y(x) gives curl_z = dv_y/dx
    core::SimulationState state = makeRootState(false);
    makeUniform(state);
    auto policy = physicsPolicy(amr::AmrRegridPassKind::kRefineOnly);
    policy.physics.vorticity = 0.1;
    CHUI_REQUIRE(amr::planProductionAmrRegrid(state, {}, policy).operations.empty());
    for (std::uint32_t row = 0; row < 8U; ++row) {
      state.gas_cells.velocity_y_peculiar[row] = (row % 2U == 1U) ? 0.5 : 0.0;  // |w| dx / c_s ~ 0.39
    }
    const auto plan = amr::planProductionAmrRegrid(state, {}, policy);
    CHUI_REQUIRE(plan.operations.size() == 1U && plan.optional_refine_requests == 1U);
  }
  {  // species interface on the hybrid composition sidecar
    core::SimulationState state = makeRootState(false);
    makeUniform(state);
    auto policy = physicsPolicy(amr::AmrRegridPassKind::kRefineOnly);
    // Initial species fractions differ by <~0.06 (metallicity-dependent H/He
    // split from makeRootState); 0.15 separates that from the 0.3 jump below.
    policy.physics.species_jump = 0.15;
    CHUI_REQUIRE(amr::planProductionAmrRegrid(state, {}, policy).operations.empty());
    const auto& d = state.gas_composition.descriptor();
    CHUI_REQUIRE(d.species_lane_count > 0U);
    double gas = 0.0;
    for (std::uint32_t l = 0; l < d.element_lane_count; ++l) gas += state.gas_composition.value(d.element_lane_begin + l, 0U);
    state.gas_composition.lane(d.species_lane_begin)[0] += 0.3 * gas;  // indicator only; not conservation data
    CHUI_REQUIRE(amr::planProductionAmrRegrid(state, {}, policy).optional_refine_requests == 1U);
  }
  {  // cooling length (mandatory) and its derefinement guard
    core::SimulationState state = makeRootState(false);
    makeUniform(state);
    std::vector<double> t_cool(8U, std::numeric_limits<double>::infinity());
    auto policy = physicsPolicy(amr::AmrRegridPassKind::kRefineOnly);
    policy.physics.cooling_length_cells = 4.0;
    policy.physics.cooling_time_code = t_cool;
    CHUI_REQUIRE(amr::planProductionAmrRegrid(state, {}, policy).operations.empty());
    // c_s t_cool = 0.5 dx  => resolution 4 dx / (0.5 dx) = 8 > 1: refine.
    const double cs = state.gas_cells.sound_speed_code[3];
    t_cool[3] = 0.5 * 0.5 / cs;
    policy.physics.cooling_time_code = t_cool;
    const auto plan = amr::planProductionAmrRegrid(state, {}, policy);
    CHUI_REQUIRE(plan.operations.size() == 1U && plan.mandatory_refine_requests == 1U);
    policy.max_operations = 0U;
    const auto starved = amr::planProductionAmrRegrid(state, {}, policy);
    CHUI_REQUIRE(starved.operations.empty() && starved.mandatory_truncated_by_budget == 1U);
    // Missing per-row cooling times are an error, not a silent no-op.
    policy.physics.cooling_time_code = {};
    bool refused = false;
    try { (void)amr::planProductionAmrRegrid(state, {}, policy); } catch (const std::invalid_argument&) { refused = true; }
    CHUI_REQUIRE(refused);
  }
}

}  // namespace

int main() {
  testPhysicsIndicators();
  testRefineDerefineConservesHydroAndComposition();
  testMemoryRefusalLeavesStateUntouched();
  testBalanceAndMaxLevel();
  testRankStripedWindowsAreDisjoint();
  testDistributedCompositionAccumulator();
  return 0;
}
