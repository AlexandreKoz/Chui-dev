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
  if (m_config.adaptive.enabled) {
    if (m_config.solver.backend != ThermochemistryBackend::kFixedDetailed) {
      throw std::invalid_argument("adaptive chemistry fidelity requires the fixed_detailed backend");
    }
    ThermochemistryConfig hybrid = m_config.solver;
    hybrid.backend = ThermochemistryBackend::kFixedHybrid;
    // Demoted rows keep their ladder distribution frozen (like the C+/CO
    // reservoir); the hybrid model owns no metal-ion state.
    hybrid.metal_ion_table = nullptr;
    hybrid.metal_ion_elements = {};
    m_hybrid_solver = std::make_unique<ThermochemistrySolver>(hybrid);
  }
  if (m_config.shadow.sample_fraction > 0.0 && m_config.shadow.max_cells > 0U &&
      m_config.solver.backend == ThermochemistryBackend::kFixedDetailed) {
    ThermochemistryConfig reference = m_config.solver;
    reference.relative_tolerance = std::max(1.0e-10, 0.1 * reference.relative_tolerance);
    reference.max_substeps = std::max<std::uint32_t>(reference.max_substeps, 4U * reference.max_substeps);
    m_reference_solver = std::make_unique<ThermochemistrySolver>(reference);
  }
  if (m_layout.ionStageCount() != m_config.solver.metal_ion_elements.stageCount()) {
    throw std::invalid_argument(
        "ThermochemistryOperator: composition layout '" + m_layout.descriptor().layout_name +
        "' carries " + std::to_string(m_layout.ionStageCount()) + " metal-ion stage lanes but the solver selects '" +
        m_config.solver.metal_ion_elements.label() + "'");
  }
  if (!(m_config.units.density_code_to_cgs > 0.0) || !(m_config.units.velocity_code_to_cm_s > 0.0) ||
      !(m_config.units.length_code_to_cm > 0.0) || !(m_config.units.mass_code_to_g > 0.0)) {
    throw std::invalid_argument("ThermochemistryOperator: unit bridge must be positive");
  }
}

ThermochemistryCellInput ThermochemistryOperator::cellInput(const core::SimulationState& state, std::uint32_t row,
                                                            const SourceInterval& interval) const {
  const std::uint8_t tag = state.gas_composition.chemistry_model_tag.size() > row
      ? state.gas_composition.chemistry_model_tag[row]
      : k_chemistry_tag_detailed;
  ThermochemistryCellInput input = cellInputForTag(state, row, interval, tag);
  if (rowCarriesMetalIons(tag)) {
    MetalIonRowState ions;
    MetalIonAdvanceReport report;
    attachMetalIonState(state, row, rowTemperatureK(state, row, input), input, ions, report);
  }
  return input;
}

bool ThermochemistryOperator::rowCarriesMetalIons(std::uint8_t chemistry_tag) const noexcept {
  return m_solver.metalIonSolver() != nullptr &&
      (m_hybrid_solver == nullptr || chemistry_tag == k_chemistry_tag_detailed);
}

double ThermochemistryOperator::rowTemperatureK(const core::SimulationState& state, std::uint32_t row,
                                                const ThermochemistryCellInput& input) const {
  double mu = state.gas_composition.mean_molecular_weight[row];
  if (!(mu > 0.0) || !std::isfinite(mu)) mu = 1.22;
  return std::max((m_config.solver.adiabatic_index - 1.0) * mu * k_amu_g * input.specific_internal_energy_erg_g /
                      k_boltzmann_erg_per_k,
                  m_config.solver.temperature_floor_k);
}

void ThermochemistryOperator::attachMetalIonState(const core::SimulationState& state, std::uint32_t row,
                                                  double temperature_k, ThermochemistryCellInput& input,
                                                  MetalIonRowState& ions, MetalIonAdvanceReport& report) const {
  const MetalIonLadderSolver* solver = m_solver.metalIonSolver();
  if (solver == nullptr) throw std::logic_error("metal-ion state requested without a metal-ion solver");
  const auto& sidecar = state.gas_composition;
  const MetalIonElementSet set = m_config.solver.metal_ion_elements;
  const double n_co = m_solver.speciesCount() > static_cast<std::size_t>(ChemicalSpecies::kCO)
      ? input.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kCO)]
      : 0.0;
  MetalIonEnvironment env;
  env.temperature_k = temperature_k;
  env.electron_cm3 = m_solver.totalElectronDensity(input, input.species_cm3);
  env.redshift = input.redshift;
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    const auto element = static_cast<MetalIonElement>(e);
    if (!set.contains(element)) continue;
    const ChemicalElement chemical = chemicalElementOf(element);
    // Authoritative gas-phase nuclei: element lane (dust excluded) minus the
    // nuclei bound in CO (C and O). The ladder lanes carry the DISTRIBUTION;
    // nuclei added/removed by dust, feedback or star formation inherit the
    // current distribution (conservative proportional mapping).
    double nuclei = input.gas_element_mass_fraction[static_cast<std::size_t>(chemical)] * input.mass_density_cgs /
        (chemicalElementInfo(chemical).atomic_mass_amu * k_amu_g);
    if (chemical == ChemicalElement::kC || chemical == ChemicalElement::kO) nuclei -= n_co;
    ions.nuclei_cm3[e] = std::max(nuclei, 0.0);
    const std::size_t offset = set.stageOffset(element);
    const std::size_t stages = metalIonStageCount(element);
    double sum = 0.0;
    for (std::size_t k = 0; k < stages; ++k) {
      const double lane_mass = std::max(sidecar.value(m_layout.ionStageLane(offset + k), row), 0.0);
      ions.fractions[offset + k] = lane_mass;
      sum += lane_mass;
    }
    const std::span<double> fractions(ions.fractions.data() + offset, stages);
    if (sum > 0.0 && std::isfinite(sum)) {
      for (double& x : fractions) x /= sum;
    } else {
      // Row without ladder history (first detailed interval, promoted
      // adaptive row, pre-ion restart): explicit equilibrium initialization.
      solver->equilibrium(element, fractions, env, report);
      ions.initialized_from_equilibrium = true;
    }
  }
  refreshMetalIonInput(input, ions);
}

void ThermochemistryOperator::refreshMetalIonInput(ThermochemistryCellInput& input,
                                                   const MetalIonRowState& ions) const {
  const MetalIonElementSet set = m_config.solver.metal_ion_elements;
  input.metal_ion_stage_cm3.fill(0.0);
  double electrons = 0.0;
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    const auto element = static_cast<MetalIonElement>(e);
    if (!set.contains(element)) continue;
    const std::size_t offset = set.stageOffset(element);
    const std::size_t stages = metalIonStageCount(element);
    for (std::size_t k = 0; k < stages; ++k) {
      input.metal_ion_stage_cm3[offset + k] = ions.nuclei_cm3[e] * ions.fractions[offset + k];
    }
    electrons += ions.nuclei_cm3[e] *
        metalIonElectronsPerNucleus(std::span<const double>(ions.fractions.data() + offset, stages));
    if (element == MetalIonElement::kC &&
        m_solver.speciesCount() > static_cast<std::size_t>(ChemicalSpecies::kCII)) {
      // The NL97 network already counts one electron per C+ particle.
      electrons -= input.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kCII)];
    }
  }
  input.metal_ion_electron_cm3 = electrons;
}

ThermochemistryCellInput ThermochemistryOperator::cellInputForTag(const core::SimulationState& state,
                                                                  std::uint32_t row, const SourceInterval& interval,
                                                                  std::uint8_t /*chemistry_tag*/) const {
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

ThermochemistryRowResult ThermochemistryOperator::advanceRow(const core::SimulationState& state,
                                                             const SourceInterval& interval, std::uint32_t row,
                                                             std::uint8_t chemistry_tag,
                                                             ChemistryShadowReport& shadow) const {
  ThermochemistryRowResult result;
  const auto& sidecar = state.gas_composition;
  ThermochemistryCellInput input = cellInputForTag(state, row, interval, chemistry_tag);
  const bool demoted = m_hybrid_solver != nullptr && chemistry_tag == k_chemistry_tag_hybrid;
  const ThermochemistrySolver& solver = demoted ? *m_hybrid_solver : m_solver;
  // C02 operator split inside the source interval (Strang in the ladder):
  // ions(dt/2) at the step-start (T, n_e) -> light network + energy over dt
  // with the ion distribution frozen (ion losses evaluated implicitly at the
  // Newton temperature) -> ions(dt/2) at the end (T, n_e).
  const bool ions = !demoted && rowCarriesMetalIons(chemistry_tag);
  MetalIonRowState ion_state;
  MetalIonAdvanceReport ion_report;
  const MetalIonLadderSolver* ion_solver = m_solver.metalIonSolver();
  const MetalIonElementSet ion_set = m_config.solver.metal_ion_elements;
  const auto advance_ions = [&](ThermochemistryCellInput& target, double temperature_k,
                                double electron_cm3) -> bool {
    MetalIonEnvironment env;
    env.temperature_k = temperature_k;
    env.electron_cm3 = electron_cm3;
    env.redshift = input.redshift;
    for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
      const auto element = static_cast<MetalIonElement>(e);
      if (!ion_set.contains(element)) continue;
      const std::span<double> fractions(ion_state.fractions.data() + ion_set.stageOffset(element),
                                        metalIonStageCount(element));
      if (ion_solver->advance(element, fractions, env, 0.5 * interval.physical_seconds, ion_report) !=
          SourceAdvanceStatus::kCompleted) {
        return false;
      }
    }
    refreshMetalIonInput(target, ion_state);
    return true;
  };
  const auto ion_failure = [&](const char* phase) {
    result.status = SourceAdvanceStatus::kRetryShorterInterval;
    result.suggested_retry_seconds = 0.5 * interval.physical_seconds;
    result.diagnostic = "row " + std::to_string(row) + " (gas_cell_id " +
        std::to_string(state.gas_cells.gas_cell_id[row]) + "): metal-ion ladder " + phase +
        " half did not converge within the substep cap";
    return result;
  };
  if (ions) {
    try {
      attachMetalIonState(state, row, rowTemperatureK(state, row, input), input, ion_state, ion_report);
      if (!advance_ions(input, rowTemperatureK(state, row, input),
                        m_solver.totalElectronDensity(input, input.species_cm3))) {
        return ion_failure("first");
      }
    } catch (const std::range_error& error) {
      result.status = SourceAdvanceStatus::kRejectedInvalid;
      result.diagnostic = "row " + std::to_string(row) + ": metal-ion dataset range: " + error.what();
      return result;
    }
  }
  ThermochemistryCellOutput output;
  SourceLedger per_volume;
  solver.advanceCell(input, interval.physical_seconds, output, per_volume);
  if (output.status != SourceAdvanceStatus::kCompleted) {
    result.status = output.status;
    result.suggested_retry_seconds = output.suggested_retry_seconds;
    result.diagnostic = "row " + std::to_string(row) + " (gas_cell_id " +
        std::to_string(state.gas_cells.gas_cell_id[row]) + "): " + output.diagnostic;
    return result;
  }
  if (ions) {
    // The light network's end species update the NL97 C+ correction before
    // the closing half; the end electron density is the coupled total. The
    // step-start `input` stays untouched for the shadow reference below.
    ThermochemistryCellInput end_input = input;
    end_input.species_cm3 = output.species_cm3;
    refreshMetalIonInput(end_input, ion_state);
    try {
      if (!advance_ions(end_input, std::max(output.temperature_k, m_config.solver.temperature_floor_k),
                        m_solver.totalElectronDensity(end_input, output.species_cm3))) {
        return ion_failure("closing");
      }
    } catch (const std::range_error& error) {
      result.status = SourceAdvanceStatus::kRejectedInvalid;
      result.diagnostic = "row " + std::to_string(row) + ": metal-ion dataset range: " + error.what();
      return result;
    }
    per_volume.accepted_substeps += ion_report.accepted_substeps;
    per_volume.rejected_substeps += ion_report.rejected_substeps;
    per_volume.max_error_estimate = std::max(per_volume.max_error_estimate, ion_report.max_error_estimate);
    per_volume.charge_residual = std::max(per_volume.charge_residual, ion_report.normalization_residual);
    per_volume.metal_table_clamps += ion_report.range_clamps;
  }
  if (m_reference_solver != nullptr && shadow.samples < m_config.shadow.max_cells) {
    // Deterministic sample keyed by stable id and accepted step.
    const std::uint64_t key = core::fnv1a64(std::to_string(state.gas_cells.gas_cell_id[row]) + ":" +
                                            std::to_string(interval.accepted_step_index));
    const double u01 = static_cast<double>(key >> 11U) * (1.0 / 9007199254740992.0);
    if (u01 < m_config.shadow.sample_fraction) {
      ThermochemistryCellOutput reference;
      SourceLedger reference_ledger;
      m_reference_solver->advanceCell(input, interval.physical_seconds, reference, reference_ledger);
      if (reference.status == SourceAdvanceStatus::kCompleted) {
        ++shadow.samples;
        NetworkState production = output.species_cm3;
        if (demoted) {
          for (std::size_t sp = solver.speciesCount(); sp < m_solver.speciesCount(); ++sp) {
            production[sp] = input.species_cm3[sp];  // frozen reservoir
          }
        }
        double species = 0.0;
        for (std::size_t sp = 0; sp < m_solver.speciesCount(); ++sp) {
          const double scale = std::max({std::abs(reference.species_cm3[sp]), std::abs(production[sp]),
                                         1.0e-10 * std::max(input.species_cm3[0], 1.0e-300)});
          species = std::max(species, std::abs(reference.species_cm3[sp] - production[sp]) / scale);
        }
        shadow.max_relative_species_difference = std::max(shadow.max_relative_species_difference, species);
        shadow.max_relative_electron_difference = std::max(
            shadow.max_relative_electron_difference,
            std::abs(reference.electron_per_hydrogen - output.electron_per_hydrogen) /
                std::max(reference.electron_per_hydrogen, 1.0e-12));
        shadow.max_relative_energy_difference = std::max(
            shadow.max_relative_energy_difference,
            std::abs(reference.specific_internal_energy_erg_g - output.specific_internal_energy_erg_g) /
                reference.specific_internal_energy_erg_g);
        const double ratio =
            static_cast<double>(reference_ledger.accepted_substeps + reference_ledger.rejected_substeps) /
            std::max<double>(1.0, static_cast<double>(per_volume.accepted_substeps + per_volume.rejected_substeps));
        shadow.mean_substep_cost_ratio += (ratio - shadow.mean_substep_cost_ratio) / static_cast<double>(shadow.samples);
      }
    }
  }
  // Volume ledger: proper cell volume = (mass / rho_code) a^3 L^3.
  const double a = interval.evaluation_scale_factor;
  const double volume_code = state.cells.mass_code[row] / std::max(state.gas_cells.density_code[row], 1.0e-300);
  const double volume_cm3 = volume_code * (m_config.units.comoving ? a * a * a : 1.0) *
      std::pow(m_config.units.length_code_to_cm, 3.0);
  result.ledger = per_volume;
  result.ledger.radiative_loss_erg *= volume_cm3;
  result.ledger.heating_gain_erg *= volume_cm3;
  result.ledger.chemical_binding_delta_erg *= volume_cm3;
  result.ledger.floor_injection_erg *= volume_cm3;
  result.ledger.thermal_energy_delta_erg *= volume_cm3;

  ThermochemistryStagedRow& staged = result.staged;
  staged.row = row;
  staged.internal_energy_code = output.specific_internal_energy_erg_g /
      (m_config.units.velocity_code_to_cm_s * m_config.units.velocity_code_to_cm_s);
  staged.mean_molecular_weight = output.mean_molecular_weight;
  staged.temperature_k = output.temperature_k;
  const double mass = compositionGasMass(sidecar, row);
  for (std::size_t s = 0; s < m_solver.speciesCount(); ++s) {
    if (demoted && s >= solver.speciesCount()) {
      staged.species_mass_code[s] = sidecar.value(m_layout.speciesLane(s), row);  // frozen reservoir
      continue;
    }
    const double amu = chemicalSpeciesInfo(static_cast<ChemicalSpecies>(s)).mass_amu;
    staged.species_mass_code[s] = output.species_cm3[s] * amu * k_amu_g / input.mass_density_cgs * mass;
  }
  if (ions) {
    // Stage lanes = gas-phase nuclei mass of each stage, in code mass, from
    // the authoritative nuclei density (element lane minus CO) at this row.
    staged.ion_stage_count = static_cast<std::uint32_t>(ion_set.stageCount());
    staged.ion_state_initialized_from_equilibrium = ion_state.initialized_from_equilibrium;
    for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
      const auto element = static_cast<MetalIonElement>(e);
      if (!ion_set.contains(element)) continue;
      const double amu = chemicalElementInfo(chemicalElementOf(element)).atomic_mass_amu;
      const std::size_t offset = ion_set.stageOffset(element);
      for (std::size_t k = 0; k < metalIonStageCount(element); ++k) {
        staged.ion_stage_mass_code[offset + k] = ion_state.nuclei_cm3[e] * ion_state.fractions[offset + k] *
            amu * k_amu_g / input.mass_density_cgs * mass;
      }
    }
  }
  if (!std::isfinite(staged.internal_energy_code) || !(staged.internal_energy_code > 0.0)) {
    result.status = SourceAdvanceStatus::kRejectedInvalid;
    result.diagnostic = "thermochemistry produced a non-positive internal energy at row " + std::to_string(row);
  }
  return result;
}

void ThermochemistryOperator::setGlobalPromotionBias(double bias) {
  if (!std::isfinite(bias) || bias < 0.0) {
    throw std::invalid_argument("chemistry promotion bias must be finite and >= 0");
  }
  m_promotion_bias = bias;
}

double nextChemistryPromotionBias(double bias, double global_max_shadow_difference,
                                  std::uint64_t global_samples) noexcept {
  if (global_samples == 0U || !std::isfinite(global_max_shadow_difference)) return bias;
  if (global_max_shadow_difference > 0.1) return std::min(bias + 0.5, 4.0);
  if (global_max_shadow_difference < 0.01) return std::max(bias - 0.25, 0.0);
  return bias;
}

double ThermochemistryOperator::coolingTimeSeconds(
    const core::SimulationState& state, std::uint32_t row,
    const SourceInterval& interval) const {
  const ThermochemistryCellInput input = cellInput(state, row, interval);
  double mu = state.gas_composition.mean_molecular_weight[row];
  if (!(mu > 0.0) || !std::isfinite(mu)) mu = 1.22;
  const double temperature_k = (m_config.solver.adiabatic_index - 1.0) * mu * k_amu_g *
      input.specific_internal_energy_erg_g / k_boltzmann_erg_per_k;
  if (!(temperature_k > 0.0) || !(input.mass_density_cgs > 0.0) ||
      !(input.specific_internal_energy_erg_g > 0.0)) {
    return std::numeric_limits<double>::infinity();
  }
  const CoolingHeatingTerms terms =
      m_solver.coolingHeating(input, input.species_cm3, temperature_k);
  const double loss =
      terms.radiative_cooling + terms.chemical_cooling - terms.heating;
  const double energy_density =
      input.mass_density_cgs * input.specific_internal_energy_erg_g;
  if (!(loss > 0.0) || !(energy_density > 0.0) || !std::isfinite(loss)) {
    return std::numeric_limits<double>::infinity();
  }
  const double cooling_time = energy_density / loss;
  return std::isfinite(cooling_time) && cooling_time > 0.0
      ? cooling_time : std::numeric_limits<double>::infinity();
}

ChemistryFidelityDecision ThermochemistryOperator::decideChemistryFidelity(
    const core::SimulationState& state, const SourceInterval& interval, std::uint32_t row) const {
  const auto& sidecar = state.gas_composition;
  ChemistryFidelityDecision decision;
  decision.tag = sidecar.chemistry_model_tag[row];
  decision.dwell_years = sidecar.model_dwell_years[row];
  const ThermochemistryCellInput input = cellInput(state, row, interval);
  double mu = sidecar.mean_molecular_weight[row];
  if (!(mu > 0.0) || !std::isfinite(mu)) mu = 1.22;
  const double t = (m_config.solver.adiabatic_index - 1.0) * mu * k_amu_g *
      input.specific_internal_energy_erg_g / k_boltzmann_erg_per_k;
  // A03/R3: every Part 2 backend publishes a net cooling time from the same
  // rate set that owns its source update. This diagnostic is independent of
  // adaptive fidelity so mandatory cooling-length AMR has a production data
  // path in both fixed-reference and adaptive modes.
  double t_cool = std::numeric_limits<double>::infinity();
  if (interval.physical_seconds > 0.0 && t > 0.0) {
    const CoolingHeatingTerms terms = m_solver.coolingHeating(input, input.species_cm3, t);
    const double loss = terms.radiative_cooling + terms.chemical_cooling - terms.heating;  // erg cm^-3 s^-1
    const double u = input.mass_density_cgs * input.specific_internal_energy_erg_g;     // erg cm^-3
    if (loss > 0.0 && u > 0.0) t_cool = u / loss;
  }
  decision.cooling_time_seconds = t_cool;
  if (!m_config.adaptive.enabled) return decision;
  ChemistryAdaptiveConfig c = m_config.adaptive;
  {
    const double s = 1.0 + m_promotion_bias;
    c.cold_temperature_k *= s;
    c.dense_n_h_cm3 /= s;
    c.molecular_fraction_threshold /= s;
    c.co_carbon_fraction_threshold /= s;
    c.cooling_time_ratio *= s;
    c.error_estimate_threshold /= s;
    c.residual_fraction_threshold /= s;
  }
  const double h = c.hysteresis_fraction;
  const double n_h = input.gas_element_mass_fraction[0] * input.mass_density_cgs /
      (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
  const double f_h2 = n_h > 0.0 ? 2.0 * input.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kH2)] / n_h : 0.0;
  const double n_c = input.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kCII)] +
      input.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kCO)];
  const double f_co = n_c > 0.0 ? input.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kCO)] / n_c : 0.0;
  const double stiff_horizon = c.cooling_time_ratio * interval.physical_seconds;
  const bool stiff = c.cooling_time_ratio > 0.0 && t_cool < stiff_horizon;
  const bool not_stiff = c.cooling_time_ratio <= 0.0 || t_cool > stiff_horizon * (1.0 + h);
  const double previous_error = sidecar.chemistry_last_error_estimate[row];
  const double residual_fraction = std::max(
      sidecar.chemistry_last_element_fraction_residual[row],
      sidecar.chemistry_last_charge_residual[row]);
  const bool error_history =
      c.error_estimate_threshold > 0.0 &&
      previous_error > c.error_estimate_threshold;
  const bool residual_history =
      c.residual_fraction_threshold > 0.0 &&
      residual_fraction > c.residual_fraction_threshold;
  const bool error_history_clear =
      (c.error_estimate_threshold <= 0.0 ||
       previous_error < c.error_estimate_threshold * (1.0 - h)) &&
      (c.residual_fraction_threshold <= 0.0 ||
       residual_fraction < c.residual_fraction_threshold * (1.0 - h));
  const std::uint8_t current = decision.tag;
  std::uint8_t desired = current;
  if (current == k_chemistry_tag_hybrid) {
    const bool cold_dense = t < c.cold_temperature_k * (1.0 - h) && n_h > c.dense_n_h_cm3 * (1.0 + h);
    if (cold_dense || stiff || error_history || residual_history ||
        f_h2 > c.molecular_fraction_threshold * (1.0 + h) ||
        f_co > c.co_carbon_fraction_threshold * (1.0 + h)) {
      desired = k_chemistry_tag_detailed;
      decision.cooling_time_trigger = stiff;
      decision.error_history_trigger = error_history || residual_history;
    }
  } else {
    const bool clear = not_stiff && error_history_clear &&
        t > c.cold_temperature_k * (1.0 + h) && n_h < c.dense_n_h_cm3 * (1.0 - h) &&
        f_h2 < c.molecular_fraction_threshold * (1.0 - h) && f_co < c.co_carbon_fraction_threshold * (1.0 - h);
    if (clear) desired = k_chemistry_tag_hybrid;
  }
  decision.dwell_years += interval.physical_seconds / k_seconds_per_year;
  if (desired != current) {
    // Promotion is never delayed (the more complete model is always safe);
    // demotion waits for the dwell time.
    if (desired == k_chemistry_tag_hybrid && decision.dwell_years < c.minimum_dwell_years) {
      decision.held_by_dwell = true;
    } else {
      decision.tag = desired;
      decision.dwell_years = 0.0;
      decision.promoted = desired == k_chemistry_tag_detailed;
      decision.demoted = desired == k_chemistry_tag_hybrid;
    }
  }
  return decision;
}

}  // namespace cosmosim::physics
