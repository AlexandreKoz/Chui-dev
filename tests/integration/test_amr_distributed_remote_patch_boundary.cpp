#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_gamma = 1.4;

void setPatch(
    cosmosim::core::SimulationState& state,
    std::size_t patch_index,
    const cosmosim::amr::PatchDescriptor& descriptor,
    std::uint32_t first_cell,
    std::uint32_t cell_count,
    std::uint32_t owner_rank) {
  cosmosim::amr::writePatchDescriptorToStateRow(state, patch_index, descriptor);
  state.patches.first_cell[patch_index] = first_cell;
  state.patches.cell_count[patch_index] = cell_count;
  state.patches.owning_rank[patch_index] = owner_rank;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    double x,
    double rho,
    double vx,
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.mass_code[row] = rho;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0U;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.velocity_x_peculiar[row] = vx;
  state.gas_cells.velocity_y_peculiar[row] = 0.0;
  state.gas_cells.velocity_z_peculiar[row] = 0.0;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  records.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = std::nullopt,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::amr::PatchDescriptor coarsePatch() {
  return cosmosim::amr::PatchDescriptor{
      .patch_id = 101,
      .level = 0,
      .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}};
}

[[nodiscard]] cosmosim::amr::PatchDescriptor finePatch() {
  return cosmosim::amr::PatchDescriptor{
      .patch_id = 201,
      .parent_patch_id = 101,
      .level = 1,
      .morton_key = 201,
      .origin_comov = {1.0, 0.0, 0.0},
      .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}};
}

[[nodiscard]] cosmosim::core::SimulationState makeFineRankState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(2);
  state.resizePatches(1);
  setPatch(state, 0, finePatch(), 0, 2, 1);
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  setCell(state, 0, 1.125, 0.8, -0.25, 0.8, 0, 201, 9101, records);
  setCell(state, 1, 1.375, 0.7, -0.10, 0.7, 0, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

[[nodiscard]] cosmosim::amr::DistributedAmrRemotePatch makeRemoteCoarsePatch() {
  cosmosim::amr::DistributedAmrRemotePatch remote;
  remote.patch = coarsePatch();
  remote.owner_rank = 0;
  remote.ghost_hydro_epoch = 7;
  remote.expected_ghost_hydro_epoch = 7;
  // M2A.2 sparse remote contract: only the coarse patch face touching the
  // local fine patch exists in memory. Offset 0 has no placeholder state.
  remote.boundary_cells.push_back(cosmosim::amr::AmrHydroSparseRemoteCell{
      .patch_local_cell = 1U,
      .gas_cell_id = 9002U,
      .conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
          {.rho_comoving = 0.9, .vel_x_peculiar = 0.35, .pressure_comoving = 0.9}, k_gamma),
      .internal_energy_code = 0.9 / ((k_gamma - 1.0) * 0.9)});
  return remote;
}

void testRemoteCoarseGhostProducesRemoteRefluxRegister() {
  cosmosim::core::SimulationState state = makeFineRankState();
  const cosmosim::amr::DistributedAmrRemotePatch remote = makeRemoteCoarsePatch();
  std::vector<cosmosim::amr::FluxRegisterEntry> outbound;
  cosmosim::hydro::HydroCoreSolver solver(
      k_gamma,
      cosmosim::hydro::HydroNumericsPolicy{
          .dual_energy = cosmosim::hydro::HydroDualEnergyPolicy{
              .mode = cosmosim::hydro::HydroDualEnergyMode::kThermalFraction,
              .thermal_fraction_switch = 1.0e-3}});
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  const std::vector<std::uint32_t> active_rows{0, 1};
  const auto diagnostics = cosmosim::amr::advanceDistributedProductionAmrHydro(
      state,
      active_rows,
      update,
      source_context,
      solver,
      riemann,
      {},
      cosmosim::amr::ProductionAmrHydroOptions{
          .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
          .adiabatic_index = k_gamma,
          .density_floor = 1.0e-12,
          .pressure_floor = 1.0e-12,
          .state_time_code = 0.0,
          .ghost_fill_time_code = 0.0},
      cosmosim::amr::DistributedAmrHydroExchange{
          .local_rank = 1,
          .ghost_hydro_epoch = 7,
          .expected_ghost_hydro_epoch = 7,
          .remote_patches = std::span<const cosmosim::amr::DistributedAmrRemotePatch>(&remote, 1),
          .outbound_remote_flux_registers = &outbound});

  assert(diagnostics.advanced_patch_count == 1U);
  assert(diagnostics.ghost_fill.coarse_to_fine_ghosts_filled > 0U);
  // H02: the imported coarse owner supplied canonical specific internal
  // energy, so dual-energy ghost seeding must not fall back to E-K.
  assert(diagnostics.dual_energy_unresolved_ghosts == 0U);
  assert(!outbound.empty());
  bool saw_remote_coarse_target = false;
  for (const auto& entry : outbound) {
    if (entry.coarse_patch_id == 101 && entry.coarse_gas_cell_id == 9002 &&
        entry.fine_face_count > 0U && entry.coarse_face_count == 0U) {
      saw_remote_coarse_target = true;
    }
  }
  assert(saw_remote_coarse_target);
}

void testSparseRemoteResidencyTracksInterfaceAreaNotPatchVolume() {
  constexpr std::size_t face_cells = 64U * 64U;
  auto make_sparse_patch = [](std::uint64_t patch_id, std::uint16_t depth) {
    cosmosim::amr::DistributedAmrRemotePatch remote;
    remote.patch = cosmosim::amr::PatchDescriptor{
        .patch_id = patch_id,
        .level = 0,
        .morton_key = patch_id,
        .origin_comov = {0.0, 0.0, 0.0},
        .extent_comov = {1.0, 1.0, 1.0},
        .cell_dims = {depth, 64U, 64U}};
    remote.owner_rank = 0;
    remote.boundary_cells.reserve(face_cells);
    for (std::uint32_t k = 0; k < 64U; ++k) {
      for (std::uint32_t j = 0; j < 64U; ++j) {
        const std::uint32_t offset =
            static_cast<std::uint32_t>(depth - 1U) +
            static_cast<std::uint32_t>(depth) * (j + 64U * k);
        remote.boundary_cells.push_back(cosmosim::amr::AmrHydroSparseRemoteCell{
            .patch_local_cell = offset,
            .gas_cell_id = 100000U + static_cast<std::uint64_t>(j + 64U * k),
            .conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
                {.rho_comoving = 1.0, .vel_x_peculiar = 0.0, .pressure_comoving = 1.0},
                k_gamma)});
      }
    }
    return remote;
  };

  auto shallow = make_sparse_patch(301U, 4U);
  auto deep = make_sparse_patch(302U, 128U);
  assert(shallow.boundary_cells.size() == face_cells);
  assert(deep.boundary_cells.size() == face_cells);
  assert(shallow.boundary_cells.capacity() == deep.boundary_cells.capacity());
  const std::size_t shallow_bytes = shallow.boundary_cells.capacity() *
      sizeof(cosmosim::amr::AmrHydroSparseRemoteCell);
  const std::size_t deep_bytes = deep.boundary_cells.capacity() *
      sizeof(cosmosim::amr::AmrHydroSparseRemoteCell);
  assert(shallow_bytes == deep_bytes);
  assert(static_cast<std::uint64_t>(shallow.patch.cell_dims[0]) * 64U * 64U !=
         static_cast<std::uint64_t>(deep.patch.cell_dims[0]) * 64U * 64U);
}

void testMissingSparseRemoteCellFailsDeterministically() {
  cosmosim::core::SimulationState state = makeFineRankState();
  cosmosim::amr::DistributedAmrRemotePatch remote = makeRemoteCoarsePatch();
  remote.boundary_cells.clear();
  std::vector<cosmosim::amr::FluxRegisterEntry> outbound;
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0};
  bool threw = false;
  try {
    (void)cosmosim::amr::advanceDistributedProductionAmrHydro(
        state,
        std::vector<std::uint32_t>{0U, 1U},
        update,
        {.update = update},
        solver,
        riemann,
        {},
        cosmosim::amr::ProductionAmrHydroOptions{
            .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
            .adiabatic_index = k_gamma,
            .density_floor = 1.0e-12,
            .pressure_floor = 1.0e-12,
            .state_time_code = 0.0,
            .ghost_fill_time_code = 0.0},
        cosmosim::amr::DistributedAmrHydroExchange{
            .local_rank = 1,
            .ghost_hydro_epoch = 7,
            .expected_ghost_hydro_epoch = 7,
            .remote_patches = std::span<const cosmosim::amr::DistributedAmrRemotePatch>(&remote, 1),
            .outbound_remote_flux_registers = &outbound});
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void testSparseRemoteCapacityStableAcrossRepeatedHydroSteps() {
  cosmosim::amr::DistributedAmrRemotePatch remote = makeRemoteCoarsePatch();
  const std::size_t capacity_before = remote.boundary_cells.capacity();
  for (int repeat = 0; repeat < 3; ++repeat) {
    cosmosim::core::SimulationState state = makeFineRankState();
    std::vector<cosmosim::amr::FluxRegisterEntry> outbound;
    cosmosim::hydro::HydroCoreSolver solver(k_gamma);
    cosmosim::hydro::HlleRiemannSolver riemann;
    const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0};
    (void)cosmosim::amr::advanceDistributedProductionAmrHydro(
        state,
        std::vector<std::uint32_t>{0U, 1U},
        update,
        {.update = update},
        solver,
        riemann,
        {},
        cosmosim::amr::ProductionAmrHydroOptions{
            .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
            .adiabatic_index = k_gamma,
            .density_floor = 1.0e-12,
            .pressure_floor = 1.0e-12,
            .state_time_code = 0.0,
            .ghost_fill_time_code = 0.0},
        cosmosim::amr::DistributedAmrHydroExchange{
            .local_rank = 1,
            .ghost_hydro_epoch = 7,
            .expected_ghost_hydro_epoch = 7,
            .remote_patches = std::span<const cosmosim::amr::DistributedAmrRemotePatch>(&remote, 1),
            .outbound_remote_flux_registers = &outbound});
    assert(remote.boundary_cells.capacity() == capacity_before);
  }
}

void testDistributedSubcyclingIsGated() {
  cosmosim::core::SimulationState state = makeFineRankState();
  const cosmosim::amr::DistributedAmrRemotePatch remote = makeRemoteCoarsePatch();
  std::vector<cosmosim::amr::FluxRegisterEntry> outbound;
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0};
  bool threw = false;
  try {
    (void)cosmosim::amr::advanceDistributedProductionAmrHydro(
        state,
        {},
        update,
        {.update = update},
        solver,
        riemann,
        {},
        cosmosim::amr::ProductionAmrHydroOptions{
            .adiabatic_index = k_gamma,
            .sweep_mode = cosmosim::amr::ProductionAmrHydroSweepMode::kLocalSubcycling},
        cosmosim::amr::DistributedAmrHydroExchange{
            .local_rank = 1,
            .remote_patches = std::span<const cosmosim::amr::DistributedAmrRemotePatch>(&remote, 1),
            .outbound_remote_flux_registers = &outbound});
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  testRemoteCoarseGhostProducesRemoteRefluxRegister();
  testSparseRemoteResidencyTracksInterfaceAreaNotPatchVolume();
  testMissingSparseRemoteCellFailsDeterministically();
  testSparseRemoteCapacityStableAcrossRepeatedHydroSteps();
  testDistributedSubcyclingIsGated();
  return 0;
}
