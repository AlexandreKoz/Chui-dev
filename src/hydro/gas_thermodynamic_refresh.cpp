#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cosmosim::hydro {
namespace {

constexpr double k_small = 1.0e-30;

void validatePolicy(const GasThermodynamicRefreshPolicy& policy) {
  if (!(policy.adiabatic_index > 1.0) || !std::isfinite(policy.adiabatic_index)) {
    throw std::invalid_argument("gas thermodynamic refresh requires adiabatic_index > 1");
  }
  if (!(policy.density_floor_code >= 0.0) || !(policy.pressure_floor_code >= 0.0)) {
    throw std::invalid_argument("gas thermodynamic refresh requires non-negative floors");
  }
  if (!(policy.scale_factor > 0.0) || !std::isfinite(policy.scale_factor)) {
    throw std::invalid_argument("gas thermodynamic refresh requires finite scale_factor > 0");
  }
}

struct RefreshedCell {
  double pressure = 0.0;
  double sound_speed = 0.0;
  bool effective = false;
  bool floored = false;
};

[[nodiscard]] RefreshedCell evaluateCell(
    const core::SimulationState& state,
    std::size_t row,
    const GasThermodynamicRefreshPolicy& policy) {
  const auto& gas = state.gas_cells;
  const double rho = std::max(gas.density_code[row], std::max(policy.density_floor_code, k_small));
  const double u = std::max(gas.internal_energy_code[row], 0.0);
  const double vx = gas.velocity_x_peculiar[row];
  const double vy = gas.velocity_y_peculiar[row];
  const double vz = gas.velocity_z_peculiar[row];
  RefreshedCell result;
  result.pressure = (policy.adiabatic_index - 1.0) * rho * u;
  double c_squared = policy.adiabatic_index * result.pressure / rho;
  if (policy.closure != nullptr) {
    HydroPrimitiveState ideal;
    ideal.rho_comoving = rho;
    ideal.vel_x_peculiar = vx;
    ideal.vel_y_peculiar = vy;
    ideal.vel_z_peculiar = vz;
    ideal.pressure_comoving = std::max(result.pressure, k_small);
    ideal.specific_internal_energy_code = std::max(u, k_small);
    ideal.signal_speed_squared_code = std::max(c_squared, k_small);
    const double cell_mass = row < state.cells.mass_code.size() ? state.cells.mass_code[row] : 0.0;
    ideal.metallicity_mass_fraction = cell_mass > 0.0
        ? std::clamp(gas.metal_mass_code[row] / cell_mass, 0.0, 1.0)
        : 0.0;
    HydroConservedState conserved;
    conserved.mass_density_comoving = rho;
    conserved.momentum_density_x_comoving = rho * vx;
    conserved.momentum_density_y_comoving = rho * vy;
    conserved.momentum_density_z_comoving = rho * vz;
    conserved.total_energy_density_comoving = rho * u + 0.5 * rho * (vx * vx + vy * vy + vz * vz);
    conserved.metal_mass_density_comoving = rho * ideal.metallicity_mass_fraction;
    const HydroThermodynamicClosureResult closure = policy.closure->evaluate(
        row, conserved, ideal, policy.scale_factor, policy.redshift);
    if (closure.valid && closure.pressure_comoving > 0.0 && closure.signal_speed_squared_code > 0.0) {
      result.pressure = closure.pressure_comoving;
      c_squared = closure.signal_speed_squared_code;
      result.effective = closure.uses_effective_ism;
    }
  }
  if (result.pressure < policy.pressure_floor_code) {
    result.pressure = policy.pressure_floor_code;
    c_squared = std::max(c_squared, policy.adiabatic_index * result.pressure / rho);
    result.floored = true;
  }
  result.sound_speed = std::sqrt(std::max(c_squared, 0.0));
  return result;
}

void accumulate(GasThermodynamicRefreshReport& report, const RefreshedCell& cell, double old_sound_speed) {
  ++report.refreshed_cells;
  if (cell.effective) {
    ++report.effective_closure_cells;
  }
  if (cell.floored) {
    ++report.pressure_floor_cells;
  }
  if (old_sound_speed > 0.0 && std::isfinite(old_sound_speed)) {
    report.max_sound_speed_increase_ratio =
        std::max(report.max_sound_speed_increase_ratio, cell.sound_speed / old_sound_speed);
  }
  report.max_sound_speed_code = std::max(report.max_sound_speed_code, cell.sound_speed);
}

}  // namespace

double temperatureKelvinPerInternalEnergyCode(
    double adiabatic_index, double velocity_si_per_code, double mean_molecular_weight) {
  constexpr double k_proton_mass_kg = 1.67262192369e-27;
  constexpr double k_boltzmann_j_per_k = 1.380649e-23;
  if (!(adiabatic_index > 1.0) || !(velocity_si_per_code > 0.0) || !(mean_molecular_weight > 0.0)) {
    throw std::invalid_argument(
        "temperature conversion requires gamma > 1, a positive velocity unit and mu > 0");
  }
  return (adiabatic_index - 1.0) * mean_molecular_weight * k_proton_mass_kg /
      k_boltzmann_j_per_k * velocity_si_per_code * velocity_si_per_code;
}

namespace {
// CHUI Part 2: composition-dependent temperature. When the gas composition
// sidecar carries a valid mean molecular weight for the row (written by the
// thermochemistry commit, in units of the atomic mass unit), T uses that mu
// instead of the Part 1 reference mu = 0.59; otherwise the reference factor
// applies (plain hydro / legacy runs are unchanged).
[[nodiscard]] double rowTemperatureFactor(const core::SimulationState& state, std::size_t row, double base) {
  const auto& composition = state.gas_composition;
  if (base > 0.0 && composition.enabled() && row < composition.rowCount()) {
    const double mu_amu = composition.mean_molecular_weight[row];
    if (std::isfinite(mu_amu) && mu_amu > 0.0) {
      constexpr double k_amu_over_proton_mass = 1.66053906660e-24 / 1.67262192369e-24;
      return base * mu_amu * k_amu_over_proton_mass / k_part1_reference_mean_molecular_weight;
    }
  }
  return base;
}
}  // namespace

GasThermodynamicRefreshReport refreshGasThermodynamics(
    core::SimulationState& state,
    std::span<const GasThermodynamicTouch> touched_cells,
    const GasThermodynamicRefreshPolicy& policy) {
  validatePolicy(policy);
  const std::size_t gas_count = state.gas_cells.size();
  for (const GasThermodynamicTouch& touch : touched_cells) {
    if (touch.cell_row >= gas_count) {
      throw std::out_of_range("gas thermodynamic refresh touch row is outside GasCellSidecar");
    }
  }
  GasThermodynamicRefreshReport report;
  auto& gas = state.gas_cells;
  for (const GasThermodynamicTouch& touch : touched_cells) {
    const std::size_t row = touch.cell_row;
    const RefreshedCell cell = evaluateCell(state, row, policy);
    accumulate(report, cell, gas.sound_speed_code[row]);
    gas.pressure_code[row] = cell.pressure;
    gas.sound_speed_code[row] = cell.sound_speed;
    const double u_new = gas.internal_energy_code[row];
    if (policy.temperature_k_per_internal_energy_code > 0.0) {
      gas.temperature_code[row] =
          rowTemperatureFactor(state, row, policy.temperature_k_per_internal_energy_code) * std::max(u_new, 0.0);
    } else if (touch.internal_energy_before_code > 0.0 && gas.temperature_code[row] > 0.0 && u_new > 0.0) {
      gas.temperature_code[row] *= u_new / touch.internal_energy_before_code;
    } else {
      ++report.temperature_unscaled_cells;
    }
  }
  state.markGasThermodynamicsRefreshed();
  return report;
}

bool refreshGasThermodynamicCell(
    core::SimulationState& state,
    std::size_t row,
    const GasThermodynamicRefreshPolicy& policy) {
  if (row >= state.gas_cells.size()) {
    throw std::out_of_range("gas thermodynamic cell refresh row is outside GasCellSidecar");
  }
  const RefreshedCell cell = evaluateCell(state, row, policy);
  auto& gas = state.gas_cells;
  gas.pressure_code[row] = cell.pressure;
  gas.sound_speed_code[row] = cell.sound_speed;
  if (policy.temperature_k_per_internal_energy_code > 0.0) {
    gas.temperature_code[row] = rowTemperatureFactor(state, row, policy.temperature_k_per_internal_energy_code) *
        std::max(gas.internal_energy_code[row], 0.0);
  }
  return cell.effective;
}

GasThermodynamicRefreshReport refreshAllGasThermodynamics(
    core::SimulationState& state,
    const GasThermodynamicRefreshPolicy& policy) {
  validatePolicy(policy);
  GasThermodynamicRefreshReport report;
  auto& gas = state.gas_cells;
  for (std::size_t row = 0; row < gas.size(); ++row) {
    const RefreshedCell cell = evaluateCell(state, row, policy);
    accumulate(report, cell, gas.sound_speed_code[row]);
    gas.pressure_code[row] = cell.pressure;
    gas.sound_speed_code[row] = cell.sound_speed;
    if (policy.temperature_k_per_internal_energy_code > 0.0) {
      gas.temperature_code[row] = rowTemperatureFactor(state, row, policy.temperature_k_per_internal_energy_code) *
          std::max(gas.internal_energy_code[row], 0.0);
    } else {
      ++report.temperature_unscaled_cells;
    }
  }
  state.markGasThermodynamicsRefreshed();
  return report;
}

}  // namespace cosmosim::hydro
