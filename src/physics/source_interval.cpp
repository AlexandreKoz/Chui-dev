#include "cosmosim/physics/source_interval.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cosmosim::physics {

bool SourceInterval::valid() const noexcept {
  return std::isfinite(physical_seconds) && physical_seconds > 0.0 &&
      std::isfinite(evaluation_scale_factor) && evaluation_scale_factor > 0.0 &&
      std::isfinite(evaluation_redshift) && evaluation_redshift >= 0.0 &&
      time_end_code >= time_begin_code;
}

std::string_view sourceAdvanceStatusLabel(SourceAdvanceStatus status) noexcept {
  switch (status) {
    case SourceAdvanceStatus::kCompleted: return "completed";
    case SourceAdvanceStatus::kRetryShorterInterval: return "retry_shorter_interval";
    case SourceAdvanceStatus::kRejectedInvalid: return "rejected_invalid";
  }
  return "unknown";
}

SourceLedger& SourceLedger::operator+=(const SourceLedger& rhs) {
  radiative_loss_erg += rhs.radiative_loss_erg;
  heating_gain_erg += rhs.heating_gain_erg;
  chemical_binding_delta_erg += rhs.chemical_binding_delta_erg;
  floor_injection_erg += rhs.floor_injection_erg;
  thermal_energy_delta_erg += rhs.thermal_energy_delta_erg;
  element_mass_residual_code = std::max(element_mass_residual_code, rhs.element_mass_residual_code);
  charge_residual = std::max(charge_residual, rhs.charge_residual);
  metal_ion_potential_delta_erg += rhs.metal_ion_potential_delta_erg;
  metal_ion_co_transfer_nuclei += rhs.metal_ion_co_transfer_nuclei;
  metal_ion_normalization_residual = std::max(metal_ion_normalization_residual, rhs.metal_ion_normalization_residual);
  metal_ion_inventory_residual = std::max(metal_ion_inventory_residual, rhs.metal_ion_inventory_residual);
  dust_to_gas_transfer_code += rhs.dust_to_gas_transfer_code;
  gas_to_dust_transfer_code += rhs.gas_to_dust_transfer_code;
  cells += rhs.cells;
  accepted_substeps += rhs.accepted_substeps;
  rejected_substeps += rhs.rejected_substeps;
  newton_iterations += rhs.newton_iterations;
  floor_activations += rhs.floor_activations;
  metal_table_clamps += rhs.metal_table_clamps;
  max_error_estimate = std::max(max_error_estimate, rhs.max_error_estimate);
  max_newton_residual = std::max(max_newton_residual, rhs.max_newton_residual);
  return *this;
}

bool sourceIntervalClosed(long double completed_seconds, long double requested_seconds) noexcept {
  if (!(requested_seconds > 0.0L) || !std::isfinite(static_cast<double>(completed_seconds))) {
    return false;
  }
  const long double tolerance =
      4.0L * static_cast<long double>(std::numeric_limits<double>::epsilon()) * requested_seconds;
  return std::fabs(completed_seconds - requested_seconds) <= tolerance;
}

}  // namespace cosmosim::physics
