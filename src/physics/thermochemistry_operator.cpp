#include "cosmosim/physics/thermochemistry_operator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

ThermochemistryOperator::ThermochemistryOperator(ThermochemistryOperatorConfig config, CompositionLayout layout)
    : m_config(config), m_layout(std::move(layout)), m_solver(config.solver) {
  const std::size_t network_species = m_solver.speciesCount();
  if (network_species > m_layout.descriptor().species_lane_count) {
    throw std::invalid_argument(
        "ThermochemistryOperator: composition layout '" + m_layout.descriptor().layout_name +
        "' has fewer species lanes than the selected network requires");
  }
  if (!(m_config.units.density_code_to_cgs > 0.0) || !(m_config.units.velocity_code_to_cm_s > 0.0) ||
      !(m_config.units.length_code_to_cm > 0.0) || !(m_config.units.mass_code_to_g > 0.0)) {
    throw std::invalid_argument("ThermochemistryOperator: unit bridge must be positive");
  }
}

ThermochemistryCellInput ThermochemistryOperator::cellInput(const core::SimulationState& state, std::uint32_t row,
                                                            const SourceInterval& interval) const {
  const auto& sidecar = state.gas_composition;
  ThermochemistryCellInput input;
  const double a = interval.evaluation_scale_factor;
  const double rho_code = state.gas_cells.density_code[row];
  const double frame = m_config.units.comoving ? 1.0 / (a * a * a) : 1.0;
  input.mass_density_cgs = rho_code * frame * m_config.units.density_code_to_cgs;
  input.specific_internal_energy_erg_g =
      state.gas_cells.internal_energy_code[row] * m_config.units.velocity_code_to_cm_s *
      m_config.units.velocity_code_to_cm_s;
  input.redshift = interval.evaluation_redshift;
  const double mass = compositionGasMass(sidecar, row);
  if (!(mass > 0.0)) {
    throw std::runtime_error("thermochemistry: composition row " + std::to_string(row) + " has no mass");
  }
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    input.gas_element_mass_fraction[e] =
        sidecar.value(m_layout.elementLane(static_cast<ChemicalElement>(e)), row) / mass;
  }
  double dust_mass = 0.0;
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    dust_mass += sidecar.value(m_layout.dustLane(d), row);
  }
  if (sidecar.descriptor().dust_lane_count > 0U) {
    input.dust_to_gas_relative_mw = (dust_mass / mass) / m_config.milky_way_dust_to_gas;
  } else if (m_config.metallicity_scaled_dust_proxy) {
    double z = 0.0;
    for (std::size_t e = 2; e < k_chemical_element_count; ++e) z += input.gas_element_mass_fraction[e];
    input.dust_to_gas_relative_mw = z / solarMetalMassFraction();
  }
  for (std::size_t s = 0; s < m_solver.speciesCount(); ++s) {
    const double species_mass = sidecar.value(m_layout.speciesLane(s), row);
    const double amu = chemicalSpeciesInfo(static_cast<ChemicalSpecies>(s)).mass_amu;
    input.species_cm3[s] = species_mass / mass * input.mass_density_cgs / (amu * k_amu_g);
  }
  return input;
}

void ThermochemistryOperator::prepare(const core::SimulationState& state, const SourceInterval& interval,
                                      std::span<const std::uint32_t> rows) {
  m_prepared = false;
  m_staged.clear();
  m_ledger = SourceLedger{};
  m_skipped_effective = 0U;
  if (!interval.valid()) {
    throw std::invalid_argument("thermochemistry: invalid source interval");
  }
  const auto& sidecar = state.gas_composition;
  if (!sidecar.enabled() || !(sidecar.descriptor() == m_layout.descriptor())) {
    throw SourceModelUnavailableError(
        "thermochemistry requires the gas composition sidecar with layout '" + m_layout.descriptor().layout_name + "'");
  }
  sidecar.requireAlignedWith(state.gas_cells.gas_cell_id, "ThermochemistryOperator::prepare");
  const std::size_t row_count = rows.empty() ? state.cells.size() : rows.size();
  m_staged.reserve(row_count);
  double min_retry = std::numeric_limits<double>::infinity();
  std::string first_failure;
  std::size_t failures = 0;
  for (std::size_t i = 0; i < row_count; ++i) {
    const std::uint32_t row = rows.empty() ? static_cast<std::uint32_t>(i) : rows[i];
    if (row >= state.cells.size()) throw std::out_of_range("thermochemistry: active row outside gas state");
    if (sidecar.ism_regime_tag[row] == k_ism_regime_effective) {
      ++m_skipped_effective;
      continue;
    }
    const ThermochemistryCellInput input = cellInput(state, row, interval);
    ThermochemistryCellOutput output;
    SourceLedger per_volume;
    m_solver.advanceCell(input, interval.physical_seconds, output, per_volume);
    if (output.status != SourceAdvanceStatus::kCompleted) {
      ++failures;
      if (first_failure.empty()) {
        first_failure = "row " + std::to_string(row) + " (gas_cell_id " +
            std::to_string(state.gas_cells.gas_cell_id[row]) + "): " + output.diagnostic;
      }
      min_retry = std::min(min_retry, output.suggested_retry_seconds);
      continue;
    }
    // Volume ledger: proper cell volume = (mass / rho_code) a^3 L^3.
    const double a = interval.evaluation_scale_factor;
    const double volume_code = state.cells.mass_code[row] / std::max(state.gas_cells.density_code[row], 1.0e-300);
    const double volume_cm3 = volume_code * (m_config.units.comoving ? a * a * a : 1.0) *
        std::pow(m_config.units.length_code_to_cm, 3.0);
    SourceLedger scaled = per_volume;
    scaled.radiative_loss_erg *= volume_cm3;
    scaled.heating_gain_erg *= volume_cm3;
    scaled.chemical_binding_delta_erg *= volume_cm3;
    scaled.floor_injection_erg *= volume_cm3;
    scaled.thermal_energy_delta_erg *= volume_cm3;
    m_ledger += scaled;

    ThermochemistryStagedRow staged;
    staged.row = row;
    staged.internal_energy_code = output.specific_internal_energy_erg_g /
        (m_config.units.velocity_code_to_cm_s * m_config.units.velocity_code_to_cm_s);
    staged.mean_molecular_weight = output.mean_molecular_weight;
    staged.temperature_k = output.temperature_k;
    const double mass = compositionGasMass(sidecar, row);
    for (std::size_t s = 0; s < m_solver.speciesCount(); ++s) {
      const double amu = chemicalSpeciesInfo(static_cast<ChemicalSpecies>(s)).mass_amu;
      staged.species_mass_code[s] = output.species_cm3[s] * amu * k_amu_g / input.mass_density_cgs * mass;
    }
    if (!std::isfinite(staged.internal_energy_code) || !(staged.internal_energy_code > 0.0)) {
      throw std::runtime_error("thermochemistry produced a non-positive internal energy");
    }
    m_staged.push_back(staged);
  }
  if (failures > 0U) {
    m_staged.clear();
    throw SourceIntervalIncompleteError(
        "thermochemistry could not integrate the accepted interval in " + std::to_string(failures) +
            " cell(s); first: " + first_failure,
        std::isfinite(min_retry) ? min_retry : 0.0);
  }
  m_prepared = true;
}

std::vector<std::uint32_t> ThermochemistryOperator::commit(core::SimulationState& state) {
  if (!m_prepared) throw std::logic_error("ThermochemistryOperator::commit without successful prepare");
  auto& sidecar = state.gas_composition;
  std::vector<std::uint32_t> touched;
  touched.reserve(m_staged.size());
  for (const ThermochemistryStagedRow& staged : m_staged) {
    state.gas_cells.internal_energy_code[staged.row] = staged.internal_energy_code;
    sidecar.mean_molecular_weight[staged.row] = staged.mean_molecular_weight;
    for (std::size_t s = 0; s < m_solver.speciesCount(); ++s) {
      sidecar.setValue(m_layout.speciesLane(s), staged.row, staged.species_mass_code[s]);
    }
    touched.push_back(staged.row);
  }
  m_staged.clear();
  m_prepared = false;
  return touched;
}

}  // namespace cosmosim::physics
