#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory_resource>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/governed_scratch_arena.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_gamma = 1.4;
constexpr double k_tol = 1.0e-9;

[[nodiscard]] cosmosim::amr::ProductionAmrHydroOptions regridOptions(
    const cosmosim::core::SimulationState& state,
    cosmosim::core::MemoryGovernor& governor) {
  governor.setBaselineOwnedBytes(cosmosim::core::memoryReportBaselineOwnedBytes(
      cosmosim::core::collectSimulationMemoryReport(state)));
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.regrid_memory_governor = &governor;
  return options;
}

struct Totals {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double total_energy = 0.0;
};

void setPatch(
    cosmosim::core::SimulationState& state,
    std::size_t patch_index,
    const cosmosim::amr::PatchDescriptor& descriptor,
    std::uint32_t first_cell,
    std::uint32_t cell_count) {
  cosmosim::amr::writePatchDescriptorToStateRow(state, patch_index, descriptor);
  state.patches.first_cell[patch_index] = first_cell;
  state.patches.cell_count[patch_index] = cell_count;
  state.patches.owning_rank[patch_index] = 0;
}

[[nodiscard]] cosmosim::hydro::HydroPrimitiveState primitiveForRow(
    const cosmosim::core::SimulationState& state,
    std::uint32_t row) {
  return cosmosim::hydro::HydroPrimitiveState{
      .rho_comoving = state.gas_cells.density_code[row],
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = state.gas_cells.pressure_code[row]};
}

[[nodiscard]] Totals totalState(
    const cosmosim::core::SimulationState& state,
    const std::vector<cosmosim::amr::PatchDescriptor>& descriptors) {
  Totals totals;
  for (const cosmosim::amr::PatchDescriptor& patch : descriptors) {
    const double volume = patch.extent_comov[0] * patch.extent_comov[1] * patch.extent_comov[2] /
        static_cast<double>(static_cast<std::size_t>(patch.cell_dims[0]) * patch.cell_dims[1] * patch.cell_dims[2]);
    const std::vector<std::uint32_t> rows = state.gas_cell_identity.rowsForPatch(patch.patch_id);
    for (const std::uint32_t row : rows) {
      const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
          primitiveForRow(state, row), k_gamma);
      totals.mass += conserved.mass_density_comoving * volume;
      totals.momentum_x += conserved.momentum_density_x_comoving * volume;
      totals.momentum_y += conserved.momentum_density_y_comoving * volume;
      totals.momentum_z += conserved.momentum_density_z_comoving * volume;
      totals.total_energy += conserved.total_energy_density_comoving * volume;
    }
  }
  return totals;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    double x,
    double rho,
    double velocity_x,
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.time_bin[row] = 0;
  state.cells.mass_code[row] = rho;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.velocity_x_peculiar[row] = velocity_x;
  state.gas_cells.velocity_y_peculiar[row] = 0.0;
  state.gas_cells.velocity_z_peculiar[row] = 0.0;
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  records.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = std::nullopt,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::core::SimulationState makeCoarseFineState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(4);
  state.resizePatches(2);
  setPatch(state, 0, cosmosim::amr::PatchDescriptor{
      .patch_id = 101,
      .parent_patch_id = 0,
      .level = 0,
      .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 0, 2);
  setPatch(state, 1, cosmosim::amr::PatchDescriptor{
      .patch_id = 201,
      .parent_patch_id = 101,
      .level = 1,
      .morton_key = 201,
      .origin_comov = {1.0, 0.0, 0.0},
      .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 2, 2);

  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(4);
  setCell(state, 0, 0.25, 1.0, 0.20, 1.0, 0, 101, 9001, records);
  setCell(state, 1, 0.75, 0.9, 0.35, 0.9, 0, 101, 9002, records);
  setCell(state, 2, 1.125, 0.8, -0.25, 0.8, 1, 201, 9101, records);
  setCell(state, 3, 1.375, 0.7, -0.10, 0.7, 1, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

void testPreparedGhostWorkspaceAdmission() {
  using namespace cosmosim;
  const std::uint64_t bound = amr::amrPreparedGhostWorkspaceBytes(2U, 12U);
  assert(bound == 12U * sizeof(hydro::HydroConservedState) +
      3U * sizeof(std::size_t) + 256U);
  core::MemoryGovernor tight(core::MemoryGovernorPolicy{.hard_limit_bytes = bound - 1U});
  bool rejected = false;
  try {
    core::GovernedScratchArena arena(&tight, core::MemoryClass::kPhaseResident,
                                    bound, "amr.test.ghost_snapshots");
  } catch (const core::MemoryAdmissionError&) { rejected = true; }
  assert(rejected);
  assert(tight.snapshot().committed_bytes == 0U);
  core::MemoryGovernor governor(core::MemoryGovernorPolicy{.hard_limit_bytes = bound});
  {
    core::GovernedScratchArena arena(&governor, core::MemoryClass::kPhaseResident,
                                    bound, "amr.test.ghost_snapshots");
    std::pmr::vector<std::size_t> offsets(arena.resource());
    std::pmr::vector<hydro::HydroConservedState> values(arena.resource());
    offsets.resize(3U);
    values.resize(12U);
    assert(values.capacity() * sizeof(hydro::HydroConservedState) +
           offsets.capacity() * sizeof(std::size_t) <= bound);
    assert(governor.snapshot().committed_bytes == bound);
  }
  assert(governor.snapshot().committed_bytes == 0U);
  bool overflow = false;
  try { (void)amr::amrPreparedGhostWorkspaceBytes(1U, std::numeric_limits<std::size_t>::max()); }
  catch (const std::overflow_error&) { overflow = true; }
  assert(overflow);
}

void testProductionAmrGeometryScatterAndRefluxPath() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));

  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  assert(descriptors.size() == 2U);
  const Totals before = totalState(state, descriptors);

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{
      .dt_code = 1.0e-4,
      .scale_factor = 1.0,
      .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  const std::vector<std::uint32_t> active_rows{0, 1, 2, 3};
  const cosmosim::amr::ProductionAmrHydroOptions options{
      .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
      .adiabatic_index = k_gamma,
      .density_floor = 1.0e-12,
      .pressure_floor = 1.0e-12};

  const auto diagnostics = cosmosim::amr::advanceProductionAmrHydro(
      state,
      active_rows,
      update,
      source_context,
      solver,
      riemann,
      {},
      options);

  assert(diagnostics.patch_count == 2U);
  assert(diagnostics.advanced_patch_count == 2U);
  assert(diagnostics.ghost_fill.same_level_ghosts_filled + diagnostics.ghost_fill.coarse_to_fine_ghosts_filled + diagnostics.ghost_fill.fine_to_coarse_ghosts_filled > 0U);
  assert(diagnostics.flux_register_entry_count > 0U);
  assert(diagnostics.reflux.corrected_cells > 0U);
  assert(diagnostics.reflux.corrected_mass_code >= 0.0);
  assert(diagnostics.reflux.corrected_momentum_x_code >= 0.0);
  assert(diagnostics.reflux.corrected_total_energy_code >= 0.0);
  assert(diagnostics.reflux.corrected_energy_code >= 0.0);

  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    assert(state.gas_cells.gas_cell_id[row] != 0U);
    assert(state.gas_cells.density_code[row] > 0.0);
    assert(state.gas_cells.pressure_code[row] > 0.0);
    assert(std::isfinite(state.gas_cells.velocity_x_peculiar[row]));
    assert(std::isfinite(state.gas_cells.internal_energy_code[row]));
  }

  const Totals after = totalState(state, cosmosim::amr::buildProductionAmrPatchDescriptors(state));
  assert(std::isfinite(after.mass));
  assert(std::isfinite(after.total_energy));
  assert(std::abs(after.mass - before.mass) < 1.0e-2);
}

[[nodiscard]] cosmosim::core::SimulationState makeRefineState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(8);
  state.resizePatches(1);
  setPatch(state, 0, cosmosim::amr::PatchDescriptor{
      .patch_id = 301,
      .parent_patch_id = 0,
      .level = 0,
      .morton_key = 301,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 2, 2}}, 0, 8);
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(8);
  std::uint32_t row = 0;
  for (std::uint32_t k = 0; k < 2; ++k) {
    for (std::uint32_t j = 0; j < 2; ++j) {
      for (std::uint32_t i = 0; i < 2; ++i) {
        const double rho = 1.0 + 0.05 * static_cast<double>(row);
        const double pressure = 1.0 + 0.02 * static_cast<double>(row);
        state.cells.center_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(i);
        state.cells.center_y_comoving[row] = 0.25 + 0.5 * static_cast<double>(j);
        state.cells.center_z_comoving[row] = 0.25 + 0.5 * static_cast<double>(k);
        state.cells.patch_index[row] = 0;
        state.cells.time_bin[row] = 0;
        state.cells.mass_code[row] = rho * 0.125;
        state.gas_cells.gas_cell_id[row] = 9201 + row;
        state.gas_cells.parent_particle_id[row] = 0;
        state.gas_cells.density_code[row] = rho;
        state.gas_cells.pressure_code[row] = pressure;
        state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
        state.gas_cells.velocity_x_peculiar[row] = 0.1 * static_cast<double>(i);
        state.gas_cells.velocity_y_peculiar[row] = -0.05 * static_cast<double>(j);
        state.gas_cells.velocity_z_peculiar[row] = 0.02 * static_cast<double>(k);
        state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
        state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
        records.push_back(cosmosim::core::GasCellIdentityRecord{
            .gas_cell_id = 9201 + row,
            .parent_particle_id = std::nullopt,
            .owning_patch_id = 301,
            .local_cell_row = row});
        ++row;
      }
    }
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

void assertClose(double lhs, double rhs) {
  assert(std::abs(lhs - rhs) < k_tol);
}

void rebuildIdentityFromCurrentRows(cosmosim::core::SimulationState& state) {
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const std::uint32_t patch_index = state.cells.patch_index[row];
    records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = state.gas_cells.gas_cell_id[row],
        .parent_particle_id = std::nullopt,
        .owning_patch_id = state.patches.patch_id[patch_index],
        .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
}

void swapGasRows(cosmosim::core::SimulationState& state, std::uint32_t a, std::uint32_t b) {
  using std::swap;
  swap(state.cells.center_x_comoving[a], state.cells.center_x_comoving[b]);
  swap(state.cells.center_y_comoving[a], state.cells.center_y_comoving[b]);
  swap(state.cells.center_z_comoving[a], state.cells.center_z_comoving[b]);
  swap(state.cells.mass_code[a], state.cells.mass_code[b]);
  swap(state.cells.time_bin[a], state.cells.time_bin[b]);
  swap(state.cells.patch_index[a], state.cells.patch_index[b]);
  swap(state.gas_cells.gas_cell_id[a], state.gas_cells.gas_cell_id[b]);
  swap(state.gas_cells.parent_particle_id[a], state.gas_cells.parent_particle_id[b]);
  swap(state.gas_cells.velocity_x_peculiar[a], state.gas_cells.velocity_x_peculiar[b]);
  swap(state.gas_cells.velocity_y_peculiar[a], state.gas_cells.velocity_y_peculiar[b]);
  swap(state.gas_cells.velocity_z_peculiar[a], state.gas_cells.velocity_z_peculiar[b]);
  swap(state.gas_cells.density_code[a], state.gas_cells.density_code[b]);
  swap(state.gas_cells.pressure_code[a], state.gas_cells.pressure_code[b]);
  swap(state.gas_cells.internal_energy_code[a], state.gas_cells.internal_energy_code[b]);
  swap(state.gas_cells.temperature_code[a], state.gas_cells.temperature_code[b]);
  swap(state.gas_cells.sound_speed_code[a], state.gas_cells.sound_speed_code[b]);
  rebuildIdentityFromCurrentRows(state);
}

void testPatchLocalMappingSurvivesRowReorder() {
  cosmosim::core::SimulationState state = makeRefineState();
  swapGasRows(state, 0, 7);
  const auto patch = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto geometry = cosmosim::amr::buildAmrHydroPatchGeometry(state, patch);
  assert(geometry.real_cells.front().local_cell_row == 7U);
  assert(geometry.real_cells.front().gas_cell_id == state.gas_cells.gas_cell_id[7]);
  auto conserved = cosmosim::amr::loadAmrHydroConservedState(state, geometry, k_gamma);
  auto cell0 = conserved.loadCell(0);
  cell0.mass_density_comoving *= 1.125;
  conserved.storeCell(0, cell0);
  const std::uint64_t target_gas_cell_id = geometry.real_cells.front().gas_cell_id;
  cosmosim::amr::scatterAmrHydroConservedState(state, geometry, conserved, k_gamma);
  const auto row = state.gas_cell_identity.rowForGasCellId(target_gas_cell_id);
  assert(row.has_value());
  assert(*row == 7U);
  assert(state.gas_cells.density_code[*row] > 1.0);
}

void testProductionRegridRejectsIdCollisions() {
  cosmosim::core::SimulationState state = makeRefineState();
  cosmosim::core::MemoryGovernor regrid_governor;
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  bool patch_collision = false;
  try {
    (void)cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 301, 20000, regridOptions(state, regrid_governor));
  } catch (const std::invalid_argument&) {
    patch_collision = true;
  }
  assert(patch_collision);

  bool gas_collision = false;
  try {
    (void)cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 400, 9201, regridOptions(state, regrid_governor));
  } catch (const std::invalid_argument&) {
    gas_collision = true;
  }
  assert(gas_collision);
}

void testProductionRefineUsesPatchLocalGeometryAfterRowReorder() {
  cosmosim::core::SimulationState state = makeRefineState();
  cosmosim::core::MemoryGovernor regrid_governor;
  swapGasRows(state, 0, 7);
  const double patch_cell_zero_density = state.gas_cells.density_code[7];
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 400, 20000, regridOptions(state, regrid_governor));
  assert(refine.refined_patch_count == 1U);
  assertClose(refine.conserved_mass_before, refine.conserved_mass_after);
  assertClose(refine.conserved_momentum_x_before, refine.conserved_momentum_x_after);
  assertClose(refine.conserved_momentum_y_before, refine.conserved_momentum_y_after);
  assertClose(refine.conserved_momentum_z_before, refine.conserved_momentum_z_after);
  assertClose(refine.conserved_total_energy_before, refine.conserved_total_energy_after);

  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const auto child_it = std::find_if(
      descriptors.begin(),
      descriptors.end(),
      [](const cosmosim::amr::PatchDescriptor& patch) {
        return patch.parent_patch_id == 301U && patch.origin_comov[0] == 0.0 &&
            patch.origin_comov[1] == 0.0 && patch.origin_comov[2] == 0.0;
      });
  assert(child_it != descriptors.end());
  for (const std::uint32_t row : state.gas_cell_identity.rowsForPatch(child_it->patch_id)) {
    assertClose(state.gas_cells.density_code[row], patch_cell_zero_density);
  }
}

void testProductionDerefineUsesPatchLocalGeometryAfterChildRowReorder() {
  cosmosim::core::SimulationState state = makeRefineState();
  cosmosim::core::MemoryGovernor regrid_governor;
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 400, 20000, regridOptions(state, regrid_governor));
  assert(refine.refined_patch_count == 1U);

  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    const std::uint32_t first = state.patches.first_cell[patch_index];
    const std::uint32_t count = state.patches.cell_count[patch_index];
    if (count >= 8U) {
      swapGasRows(state, first, static_cast<std::uint32_t>(first + count - 1U));
    }
  }

  const auto derefine = cosmosim::amr::derefineProductionPatchInSimulationState(state, parent, 30000, regridOptions(state, regrid_governor));
  assert(derefine.derefined_patch_count == 1U);
  assertClose(derefine.conserved_mass_before, derefine.conserved_mass_after);
  assertClose(derefine.conserved_momentum_x_before, derefine.conserved_momentum_x_after);
  assertClose(derefine.conserved_momentum_y_before, derefine.conserved_momentum_y_after);
  assertClose(derefine.conserved_momentum_z_before, derefine.conserved_momentum_z_after);
  assertClose(derefine.conserved_total_energy_before, derefine.conserved_total_energy_after);
}

void testPartialActiveSetSkipsIncompleteReflux() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  const cosmosim::amr::ProductionAmrHydroOptions options{
      .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
      .adiabatic_index = k_gamma,
      .density_floor = 1.0e-12,
      .pressure_floor = 1.0e-12};
  const std::vector<std::uint32_t> active_rows{0, 1};
  const auto diagnostics = cosmosim::amr::advanceProductionAmrHydro(
      state, active_rows, update, source_context, solver, riemann, {}, options);
  assert(diagnostics.flux_register_entry_count > 0U);
  assert(diagnostics.reflux.corrected_cells == 0U);
  assert(diagnostics.reflux.skipped_incomplete_register_count > 0U);
}

void testProductionRefineAndDerefineConserveSimulationState() {
  cosmosim::core::SimulationState state = makeRefineState();
  cosmosim::core::MemoryGovernor regrid_governor;
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 400, 20000, regridOptions(state, regrid_governor));
  assert(refine.refined_patch_count == 1U);
  assert(refine.created_gas_cell_count == 64U);
  assert(refine.retired_gas_cell_count == 8U);
  assert(state.patches.size() == 8U);
  assert(state.cells.size() == 64U);
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  assertClose(refine.conserved_mass_before, refine.conserved_mass_after);
  assertClose(refine.conserved_momentum_x_before, refine.conserved_momentum_x_after);
  assertClose(refine.conserved_momentum_y_before, refine.conserved_momentum_y_after);
  assertClose(refine.conserved_momentum_z_before, refine.conserved_momentum_z_after);
  assertClose(refine.conserved_total_energy_before, refine.conserved_total_energy_after);

  const auto derefine = cosmosim::amr::derefineProductionPatchInSimulationState(state, parent, 30000, regridOptions(state, regrid_governor));
  assert(derefine.derefined_patch_count == 1U);
  assert(derefine.created_gas_cell_count == 8U);
  assert(derefine.retired_gas_cell_count == 64U);
  assert(state.patches.size() == 1U);
  assert(state.cells.size() == 8U);
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  assertClose(derefine.conserved_mass_before, derefine.conserved_mass_after);
  assertClose(derefine.conserved_momentum_x_before, derefine.conserved_momentum_x_after);
  assertClose(derefine.conserved_momentum_y_before, derefine.conserved_momentum_y_after);
  assertClose(derefine.conserved_momentum_z_before, derefine.conserved_momentum_z_after);
  assertClose(derefine.conserved_total_energy_before, derefine.conserved_total_energy_after);
}

void testProductionRegridGovernorRejectionIsAtomic() {
  cosmosim::core::SimulationState state = makeRefineState();
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto cells_before = state.cells.mass_code;
  const auto gas_ids_before = state.gas_cells.gas_cell_id;
  const auto patch_ids_before = state.patches.patch_id;
  const auto identity_generation_before = state.gasCellIdentityGeneration();
  const std::uint64_t baseline = cosmosim::core::memoryReportBaselineOwnedBytes(
      cosmosim::core::collectSimulationMemoryReport(state));
  cosmosim::core::MemoryGovernor governor(
      cosmosim::core::MemoryGovernorPolicy{.hard_limit_bytes = baseline + 1U});
  governor.setBaselineOwnedBytes(baseline);
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.regrid_memory_governor = &governor;
  bool rejected = false;
  try {
    (void)cosmosim::amr::refineProductionPatchInSimulationState(
        state, parent, 400, 20000, options);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
  assert(state.cells.mass_code == cells_before);
  assert(state.gas_cells.gas_cell_id == gas_ids_before);
  assert(state.patches.patch_id == patch_ids_before);
  assert(state.gasCellIdentityGeneration() == identity_generation_before);
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  assert(governor.snapshot().rejection_count == 1U);

  // Repeat the same atomic-rejection contract for derefinement after a
  // successful governed refine has constructed the exact child octet.
  cosmosim::core::MemoryGovernor refine_governor;
  (void)cosmosim::amr::refineProductionPatchInSimulationState(
      state, parent, 400, 20000, regridOptions(state, refine_governor));
  const auto derefine_cells_before = state.cells.mass_code;
  const auto derefine_gas_ids_before = state.gas_cells.gas_cell_id;
  const auto derefine_patch_ids_before = state.patches.patch_id;
  const auto derefine_identity_generation_before = state.gasCellIdentityGeneration();
  const std::uint64_t derefine_baseline = cosmosim::core::memoryReportBaselineOwnedBytes(
      cosmosim::core::collectSimulationMemoryReport(state));
  cosmosim::core::MemoryGovernor derefine_governor(
      cosmosim::core::MemoryGovernorPolicy{.hard_limit_bytes = derefine_baseline + 1U});
  derefine_governor.setBaselineOwnedBytes(derefine_baseline);
  cosmosim::amr::ProductionAmrHydroOptions derefine_options;
  derefine_options.regrid_memory_governor = &derefine_governor;
  rejected = false;
  try {
    (void)cosmosim::amr::derefineProductionPatchInSimulationState(
        state, parent, 30000, derefine_options);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
  assert(state.cells.mass_code == derefine_cells_before);
  assert(state.gas_cells.gas_cell_id == derefine_gas_ids_before);
  assert(state.patches.patch_id == derefine_patch_ids_before);
  assert(state.gasCellIdentityGeneration() == derefine_identity_generation_before);
  assert(derefine_governor.snapshot().rejection_count == 1U);
}

void testProductionRegridRejectsPendingRefluxBeforeMutation() {
  cosmosim::core::SimulationState state = makeRefineState();
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto patch_ids_before = state.patches.patch_id;
  const auto gas_ids_before = state.gas_cells.gas_cell_id;
  cosmosim::core::PendingFluxRegisterRecord pending;
  pending.register_key = 77U;
  pending.coarse_patch_id = parent.patch_id;
  pending.coarse_gas_cell_id = state.gas_cells.gas_cell_id.front();
  state.pending_flux_registers.upsertByRegisterKey(pending);
  bool rejected = false;
  try {
    (void)cosmosim::amr::refineProductionPatchInSimulationState(
        state, parent, 400, 20000);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
  assert(state.pending_flux_registers.size() == 1U);
  assert(state.patches.patch_id == patch_ids_before);
  assert(state.gas_cells.gas_cell_id == gas_ids_before);
}

void testProductionRegridReservationDoesNotScaleWithUnrelatedParticles() {
  auto run = [](std::size_t particle_count) {
    cosmosim::core::SimulationState state = makeRefineState();
    state.resizeParticles(particle_count);
    for (std::size_t i = 0; i < particle_count; ++i) {
      state.particle_sidecar.particle_id[i] = 1000000U + i;
      state.particle_sidecar.species_tag[i] =
          static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    }
    state.rebuildSpeciesIndex();
    const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
    cosmosim::core::MemoryGovernor governor;
    const auto options = regridOptions(state, governor);
    (void)cosmosim::amr::refineProductionPatchInSimulationState(
        state, parent, 400, 20000, options);
    return governor.snapshot().peak_reserved_bytes;
  };
  const std::uint64_t without_particles = run(0U);
  const std::uint64_t with_particles = run(4096U);
  assert(without_particles > 0U);
  assert(with_particles == without_particles);
}

}  // namespace

// CHUI Part 1 stabilization (NOT RUN — USER-RESERVED): refinement writes the
// derived thermodynamic caches of every child through the runtime-wide
// policy: T = f e (not the legacy T = e in code units), P = (gamma-1) rho e,
// c_s = sqrt(gamma P / rho), consistent with the canonical e of the child.
void testRefinementWritesConsistentThermodynamicCaches() {
  cosmosim::core::SimulationState state = makeRefineState();
  cosmosim::core::MemoryGovernor regrid_governor;
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  cosmosim::amr::ProductionAmrHydroOptions options = regridOptions(state, regrid_governor);
  options.adiabatic_index = k_gamma;
  constexpr double k_factor = 1234.5;  // synthetic T per code specific energy
  options.thermodynamic_cache.temperature_k_per_internal_energy_code = k_factor;
  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 400, 20000, options);
  assert(refine.refined_patch_count == 1U);
  std::size_t child_rows = 0;
  for (std::size_t row = 0; row < state.gas_cells.size(); ++row) {
    if (state.gas_cells.gas_cell_id[row] < 20000U) {
      continue;  // only rows written by the refinement (child ids start at 20000)
    }
    ++child_rows;
    const double rho = state.gas_cells.density_code[row];
    const double e = state.gas_cells.internal_energy_code[row];
    assert(std::abs(state.gas_cells.temperature_code[row] - k_factor * e) <= 1.0e-12 * k_factor * e);
    assert(std::abs(state.gas_cells.pressure_code[row] - (k_gamma - 1.0) * rho * e) <= 1.0e-12 * rho * e);
    assert(std::abs(state.gas_cells.sound_speed_code[row] -
                    std::sqrt(k_gamma * (k_gamma - 1.0) * e)) <= 1.0e-12 * std::sqrt(e));
  }
  assert(child_rows > 0U);
}

int main() {
  testRefinementWritesConsistentThermodynamicCaches();
  testPreparedGhostWorkspaceAdmission();
  testProductionAmrGeometryScatterAndRefluxPath();
  testPatchLocalMappingSurvivesRowReorder();
  testProductionRegridRejectsIdCollisions();
  testProductionRefineUsesPatchLocalGeometryAfterRowReorder();
  testProductionDerefineUsesPatchLocalGeometryAfterChildRowReorder();
  testPartialActiveSetSkipsIncompleteReflux();
  testProductionRefineAndDerefineConserveSimulationState();
  testProductionRegridGovernorRejectionIsAtomic();
  testProductionRegridRejectsPendingRefluxBeforeMutation();
  testProductionRegridReservationDoesNotScaleWithUnrelatedParticles();
  return 0;
}
