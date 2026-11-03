// Audit F4: post-mutation thermodynamic refresh before the next CFL decision.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants:
//  * heating a static ideal-gas cell by a factor 100 in specific internal
//    energy raises P by 100, c_s by 10 and T by 100 after the refresh, and the
//    local hydro CFL timestep shrinks by exactly 10 (static gas, CFL = C dx/c);
//  * the epoch contract: notifyGasThermodynamicMutation() makes the caches
//    stale, the refresh makes them fresh again;
//  * the refresh never changes authoritative density / internal energy /
//    velocity (cannot create or destroy mass or energy);
//  * an effective-EOS closure branch transition (cold gas on the EOS heated
//    off it) is honoured: the refreshed P/c_s come from the ideal branch;
//  * out-of-range rows are rejected before any write.
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "../support/test_require.hpp"

namespace {

namespace hy = cosmosim::hydro;
constexpr double k_gamma = 5.0 / 3.0;

void seedStaticGas(cosmosim::core::SimulationState& state, std::size_t count, double u) {
  state.resizeCells(count);
  for (std::size_t i = 0; i < count; ++i) {
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.gas_cell_id[i] = 1000U + i;
    state.gas_cells.density_code[i] = 1.0;
    state.gas_cells.internal_energy_code[i] = u;
    state.gas_cells.pressure_code[i] = (k_gamma - 1.0) * u;
    state.gas_cells.sound_speed_code[i] = std::sqrt(k_gamma * (k_gamma - 1.0) * u);
    state.gas_cells.temperature_code[i] = 1.0e4;
    state.gas_cells.metal_mass_code[i] = 0.01;
  }
}

double staticCflDt(double sound_speed) {
  return cosmosim::core::computeDirectionalCflTimeStep(
      cosmosim::core::DirectionalCflTimeStepInput{
          .cell_width_axis_code = {1.0, 1.0, 1.0},
          .velocity_axis_code = {0.0, 0.0, 0.0},
          .sound_speed_code = sound_speed,
          .coordinate_frame = cosmosim::core::CoordinateFrame::kPhysical,
          .scale_factor = 1.0},
      0.4);
}

void testHeatingByHundredTightensCflByTen() {
  cosmosim::core::SimulationState state;
  seedStaticGas(state, 3U, 2.0);
  const double dt_before = staticCflDt(state.gas_cells.sound_speed_code[1]);
  CHUI_REQUIRE(state.gasThermodynamicCachesFresh());

  // Feedback-like mutation of cell 1: u x 100.
  const double u_before = state.gas_cells.internal_energy_code[1];
  state.gas_cells.internal_energy_code[1] *= 100.0;
  state.notifyGasThermodynamicMutation();
  CHUI_REQUIRE(!state.gasThermodynamicCachesFresh());

  const std::vector<hy::GasThermodynamicTouch> touched{{.cell_row = 1U, .internal_energy_before_code = u_before}};
  const auto report = hy::refreshGasThermodynamics(
      state, touched, hy::GasThermodynamicRefreshPolicy{.adiabatic_index = k_gamma});
  CHUI_REQUIRE(state.gasThermodynamicCachesFresh());
  CHUI_REQUIRE(report.refreshed_cells == 1U);
  CHUI_REQUIRE(std::abs(report.max_sound_speed_increase_ratio - 10.0) < 1.0e-12);
  CHUI_REQUIRE(std::abs(state.gas_cells.pressure_code[1] - (k_gamma - 1.0) * 200.0) < 1.0e-12);
  CHUI_REQUIRE(std::abs(state.gas_cells.temperature_code[1] - 1.0e6) < 1.0e-6);
  const double dt_after = staticCflDt(state.gas_cells.sound_speed_code[1]);
  CHUI_REQUIRE(std::abs(dt_before / dt_after - 10.0) < 1.0e-12);
  // Authoritative inputs untouched; untouched neighbours unchanged.
  CHUI_REQUIRE(state.gas_cells.internal_energy_code[1] == 200.0);
  CHUI_REQUIRE(state.gas_cells.density_code[1] == 1.0);
  CHUI_REQUIRE(std::abs(state.gas_cells.temperature_code[0] - 1.0e4) < 1.0e-9);
}

// Test double: "effective EOS" with a fixed pressure while u stays below a
// threshold; above it the gas leaves the EOS and the ideal state is kept.
class ThresholdClosure final : public hy::HydroThermodynamicClosure {
 public:
  [[nodiscard]] hy::HydroThermodynamicClosureResult evaluate(
      std::size_t, const hy::HydroConservedState& conserved, const hy::HydroPrimitiveState& ideal,
      double, double) const override {
    if (ideal.specific_internal_energy_code < 10.0) {
      return hy::HydroThermodynamicClosureResult{
          .pressure_comoving = 7.0 * conserved.mass_density_comoving,
          .signal_speed_squared_code = 7.0,
          .target_specific_internal_energy_code = 10.0,
          .uses_effective_ism = true,
          .valid = true};
    }
    return hy::HydroThermodynamicClosureResult{
        .pressure_comoving = ideal.pressure_comoving,
        .signal_speed_squared_code = ideal.signal_speed_squared_code,
        .target_specific_internal_energy_code = 10.0,
        .uses_effective_ism = false,
        .valid = true};
  }
};

void testEffectiveBranchTransition() {
  cosmosim::core::SimulationState state;
  seedStaticGas(state, 2U, 1.0);
  const ThresholdClosure closure;
  const std::vector<hy::GasThermodynamicTouch> cold{{.cell_row = 0U, .internal_energy_before_code = 1.0}};
  auto report = hy::refreshGasThermodynamics(
      state, cold, hy::GasThermodynamicRefreshPolicy{.adiabatic_index = k_gamma, .closure = &closure});
  CHUI_REQUIRE(report.effective_closure_cells == 1U);
  CHUI_REQUIRE(std::abs(state.gas_cells.pressure_code[0] - 7.0) < 1.0e-12);
  CHUI_REQUIRE(std::abs(state.gas_cells.sound_speed_code[0] - std::sqrt(7.0)) < 1.0e-12);

  state.gas_cells.internal_energy_code[0] = 100.0;  // heated off the EOS
  state.notifyGasThermodynamicMutation();
  const std::vector<hy::GasThermodynamicTouch> hot{{.cell_row = 0U, .internal_energy_before_code = 1.0}};
  report = hy::refreshGasThermodynamics(
      state, hot, hy::GasThermodynamicRefreshPolicy{.adiabatic_index = k_gamma, .closure = &closure});
  CHUI_REQUIRE(report.effective_closure_cells == 0U);
  CHUI_REQUIRE(std::abs(state.gas_cells.pressure_code[0] - (k_gamma - 1.0) * 100.0) < 1.0e-12);
}

// Runtime-wide temperature convention: T = factor * u with
// factor = (gamma-1) mu m_p / k_B * v_unit^2 (mu = 0.59). For a km/s velocity
// unit and u = 1 (km/s)^2: T = (2/3) 0.59 1.6726e-27 1e6 / 1.3806e-23 K.
void testTemperatureConventionFactor() {
  const double factor = hy::temperatureKelvinPerInternalEnergyCode(k_gamma, 1.0e3);
  const double expected = (2.0 / 3.0) * 0.59 * 1.67262192369e-27 * 1.0e6 / 1.380649e-23;
  CHUI_REQUIRE(std::abs(factor - expected) <= 1.0e-12 * expected);
  cosmosim::core::SimulationState state;
  seedStaticGas(state, 1U, 3.0);
  state.notifyGasThermodynamicMutation();
  const std::vector<hy::GasThermodynamicTouch> touched{{.cell_row = 0U, .internal_energy_before_code = 3.0}};
  (void)hy::refreshGasThermodynamics(
      state, touched,
      hy::GasThermodynamicRefreshPolicy{.adiabatic_index = k_gamma, .temperature_k_per_internal_energy_code = factor});
  CHUI_REQUIRE(std::abs(state.gas_cells.temperature_code[0] - 3.0 * factor) <= 1.0e-12 * 3.0 * factor);
}

void testOutOfRangeRowRejectedBeforeWrites() {
  cosmosim::core::SimulationState state;
  seedStaticGas(state, 2U, 1.0);
  state.gas_cells.internal_energy_code[0] = 50.0;
  state.notifyGasThermodynamicMutation();
  const double pressure_before = state.gas_cells.pressure_code[0];
  const std::vector<hy::GasThermodynamicTouch> touched{
      {.cell_row = 0U, .internal_energy_before_code = 1.0},
      {.cell_row = 5U, .internal_energy_before_code = 1.0}};
  bool rejected = false;
  try {
    (void)hy::refreshGasThermodynamics(state, touched, hy::GasThermodynamicRefreshPolicy{});
  } catch (const std::out_of_range&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);
  CHUI_REQUIRE(state.gas_cells.pressure_code[0] == pressure_before);
  CHUI_REQUIRE(!state.gasThermodynamicCachesFresh());
}

}  // namespace

int main() {
  testHeatingByHundredTightensCflByTen();
  testEffectiveBranchTransition();
  testTemperatureConventionFactor();
  testOutOfRangeRowRejectedBeforeWrites();
  return 0;
}
