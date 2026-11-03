// AMR reflux thermodynamic-cache and dual-energy regression tests
// (Part 1 stabilization).
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants:
//  * with the runtime policy (temperature factor f, effective-EOS closure),
//    a refluxed coarse cell leaves reflux with T = f e and with P / c_s from
//    the closure (no legacy T = e in code units, no ideal-gas P for an
//    effective-EOS cell);
//  * the legacy path (zero policy) is unchanged for focused callers;
//  * dual energy at reflux: a cold hypersonic coarse cell whose corrected
//    E - K is below eta E keeps its pre-reflux specific internal energy, the
//    thermal correction is ledgered, and mass/momentum are still corrected.
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "../support/test_require.hpp"

namespace {

namespace hy = cosmosim::hydro;
constexpr double k_gamma = 1.4;

void setPatch(cosmosim::core::SimulationState& state, std::size_t patch_index,
              const cosmosim::amr::PatchDescriptor& descriptor, std::uint32_t first_cell, std::uint32_t cell_count) {
  cosmosim::amr::writePatchDescriptorToStateRow(state, patch_index, descriptor);
  state.patches.first_cell[patch_index] = first_cell;
  state.patches.cell_count[patch_index] = cell_count;
  state.patches.owning_rank[patch_index] = 0;
}

void setCell(cosmosim::core::SimulationState& state, std::uint32_t row, double x, double rho, double vx,
             double pressure, std::uint32_t patch_index, std::uint64_t patch_id, std::uint64_t id,
             std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.mass_code[row] = rho;
  state.gas_cells.gas_cell_id[row] = id;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.velocity_x_peculiar[row] = vx;
  state.gas_cells.temperature_code[row] = 0.0;
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  records.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = id, .parent_particle_id = std::nullopt, .owning_patch_id = patch_id, .local_cell_row = row});
}

cosmosim::core::SimulationState makeState(double coarse_velocity, double coarse_pressure) {
  cosmosim::core::SimulationState state;
  state.resizeCells(4);
  state.resizePatches(2);
  setPatch(state, 0, cosmosim::amr::PatchDescriptor{.patch_id = 101, .level = 0, .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0}, .extent_comov = {1.0, 1.0, 1.0}, .cell_dims = {2, 1, 1}}, 0, 2);
  setPatch(state, 1, cosmosim::amr::PatchDescriptor{.patch_id = 201, .parent_patch_id = 101, .level = 1,
      .morton_key = 201, .origin_comov = {1.0, 0.0, 0.0}, .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 2, 2);
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  setCell(state, 0, 0.25, 1.0, 0.2, 1.0, 0, 101, 9001, records);
  setCell(state, 1, 0.75, 1.0, coarse_velocity, coarse_pressure, 0, 101, 9002, records);
  setCell(state, 2, 1.125, 0.8, -0.25, 0.8, 1, 201, 9101, records);
  setCell(state, 3, 1.375, 0.7, -0.10, 0.7, 1, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

cosmosim::amr::FluxRegisterEntry completeRegister() {
  return cosmosim::amr::FluxRegisterEntry{
      .register_key = 199, .coarse_patch_id = 101, .coarse_gas_cell_id = 9002, .coarse_cell_index = 1, .level = 0,
      .coarse_face_flux_code = cosmosim::amr::ConservedState{.mass_code = 0.10, .momentum_x_code = 0.20,
          .total_energy_code = 0.30},
      .fine_face_flux_code = cosmosim::amr::ConservedState{.mass_code = 0.40, .momentum_x_code = 0.70,
          .total_energy_code = 0.90},
      .face_area_comov = 0.5, .coarse_area_comov = 0.5, .fine_area_comov = 0.5, .dt_code = 1.0e-5,
      .coarse_face_count = 1, .fine_face_count = 1};
}

// Effective-EOS test double: fixed P = 3 rho, c^2 = 3 for every cell.
class FixedClosure final : public hy::HydroThermodynamicClosure {
 public:
  [[nodiscard]] hy::HydroThermodynamicClosureResult evaluate(
      std::size_t, const hy::HydroConservedState& conserved, const hy::HydroPrimitiveState&, double,
      double) const override {
    return {.pressure_comoving = 3.0 * conserved.mass_density_comoving, .signal_speed_squared_code = 3.0,
            .target_specific_internal_energy_code = 1.0, .uses_effective_ism = true, .valid = true};
  }
};

void testRefluxWritesPolicyConsistentCaches() {
  cosmosim::core::SimulationState state = makeState(0.35, 0.9);
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const FixedClosure closure;
  constexpr double k_factor = 50.0;
  const std::vector<cosmosim::amr::FluxRegisterEntry> entries{completeRegister()};
  const auto diagnostics = cosmosim::amr::applyFluxRegistersToSimulationState(
      state, entries, descriptors, k_gamma,
      cosmosim::amr::AmrThermodynamicCachePolicy{.temperature_k_per_internal_energy_code = k_factor,
                                                 .closure = &closure});
  CHUI_REQUIRE(diagnostics.corrected_cells == 1U);
  const double rho = state.gas_cells.density_code[1];
  const double e = state.gas_cells.internal_energy_code[1];
  CHUI_REQUIRE(std::abs(state.gas_cells.temperature_code[1] - k_factor * e) <= 1.0e-12 * k_factor * e);
  CHUI_REQUIRE(std::abs(state.gas_cells.pressure_code[1] - 3.0 * rho) <= 1.0e-12 * rho);
  CHUI_REQUIRE(std::abs(state.gas_cells.sound_speed_code[1] - std::sqrt(3.0)) <= 1.0e-12);
}

void testLegacyPathUnchanged() {
  cosmosim::core::SimulationState state = makeState(0.35, 0.9);
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const std::vector<cosmosim::amr::FluxRegisterEntry> entries{completeRegister()};
  (void)cosmosim::amr::applyFluxRegistersToSimulationState(state, entries, descriptors, k_gamma);
  CHUI_REQUIRE(state.gas_cells.temperature_code[1] == state.gas_cells.internal_energy_code[1]);
}

void testDualEnergyRefluxKeepsColdThermalState() {
  // Coarse cell: rho = 1, u = 100 (K = 5e3), P = 1e-8 (rho e = 2.5e-8).
  cosmosim::core::SimulationState with_dual = makeState(100.0, 1.0e-8);
  cosmosim::core::SimulationState without_dual = makeState(100.0, 1.0e-8);
  const double e_before = with_dual.gas_cells.internal_energy_code[1];
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(with_dual);
  const std::vector<cosmosim::amr::FluxRegisterEntry> entries{completeRegister()};
  const auto dual = cosmosim::amr::applyFluxRegistersToSimulationState(
      with_dual, entries, descriptors, k_gamma,
      cosmosim::amr::AmrThermodynamicCachePolicy{.dual_energy = true, .dual_energy_eta = 1.0e-3});
  (void)cosmosim::amr::applyFluxRegistersToSimulationState(without_dual, entries, descriptors, k_gamma);
  CHUI_REQUIRE(dual.dual_energy_reflux_cells == 1U);
  CHUI_REQUIRE(std::abs(with_dual.gas_cells.internal_energy_code[1] - e_before) <= 1.0e-12 * e_before);
  CHUI_REQUIRE(std::isfinite(dual.dual_energy_reflux_thermal_correction_code));
  // Mass and momentum corrections are identical with and without dual energy.
  CHUI_REQUIRE(with_dual.gas_cells.density_code[1] == without_dual.gas_cells.density_code[1]);
  CHUI_REQUIRE(with_dual.gas_cells.velocity_x_peculiar[1] == without_dual.gas_cells.velocity_x_peculiar[1]);
}

}  // namespace

int main() {
  testRefluxWritesPolicyConsistentCaches();
  testLegacyPathUnchanged();
  testDualEnergyRefluxKeepsColdThermalState();
  return 0;
}
