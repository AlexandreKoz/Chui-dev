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
    // Patch 53 (C02/C03): a demoted row has no ladder, so the losses of the
    // ladder elements (C, O, Si, Fe) on that row can only come from the
    // element-resolved metal table (the hybrid solver queries it without the
    // ladder exclusion). Without such a table those losses would silently
    // vanish on every demoted row (C/O keep only their cold fine-structure
    // terms; Si/Fe keep nothing). Refuse that combination instead.
    const bool element_table = m_config.solver.metal_mode == MetalCoolingMode::kElementTable ||
        m_config.solver.metal_mode == MetalCoolingMode::kElementTableAndFineStructure;
    if (!m_config.solver.metal_ion_elements.empty() && !element_table) {
      throw SourceModelUnavailableError(
          "adaptive chemistry fidelity with metal-ion ladders '" + m_config.solver.metal_ion_elements.label() +
          "' requires an element-resolved metal cooling table (metal_cooling_mode = element_table or "
          "element_table_and_fine_structure): demoted hybrid rows carry no ion ladder and would otherwise lose the "
          "cooling of these elements. Configure the table or disable chemistry_adaptive_selection");
    }
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

const ThermochemistrySolver& ThermochemistryOperator::fidelityModelForRow(const core::SimulationState& state,
                                                                         std::uint32_t row) const noexcept {
  const std::uint8_t tag = state.gas_composition.chemistry_model_tag.size() > row
      ? state.gas_composition.chemistry_model_tag[row]
      : k_chemistry_tag_detailed;
  return (m_hybrid_solver != nullptr && tag == k_chemistry_tag_hybrid) ? *m_hybrid_solver : m_solver;
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
  // P1-49-02: the single authoritative local radiation context (shielded
  // photo rates identical to the light network and the ion photoheating).
  MetalIonEnvironment env = m_solver.metalIonEnvironment(input, input.species_cm3, temperature_k);
  std::array<bool, k_metal_ion_element_count> from_equilibrium{};
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    const auto element = static_cast<MetalIonElement>(e);
    if (!set.contains(element)) continue;
    const ChemicalElement chemical = chemicalElementOf(element);
    // Authoritative gas-phase nuclei: element lane (dust excluded) minus the
    // nuclei bound in CO (C and O). The ladder lanes carry the DISTRIBUTION;
    // nuclei added/removed by dust, feedback or star formation inherit the
    // current distribution (conservative proportional mapping; the stage
    // reconciles the stored stage masses at commit, see
    // reconcileMetalIonStageLanes).
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
      from_equilibrium[e] = true;
    }
  }
  refreshMetalIonInput(input, ions);
  if (ions.initialized_from_equilibrium) {
    // The equilibrium ladder changes n_e, which feeds back into the
    // collisional/recombination balance: bounded fixed-point on the coupled
    // electron density (no state outside this row is touched).
    for (int iteration = 0; iteration < 8; ++iteration) {
      const double ne_before = env.electron_cm3;
      env = m_solver.metalIonEnvironment(input, input.species_cm3, temperature_k);
      for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
        if (!from_equilibrium[e]) continue;
        const auto element = static_cast<MetalIonElement>(e);
        solver->equilibrium(element,
                            std::span<double>(ions.fractions.data() + set.stageOffset(element),
                                              metalIonStageCount(element)),
                            env, report);
      }
      refreshMetalIonInput(input, ions);
      const double ne_after = m_solver.totalElectronDensity(input, input.species_cm3);
      if (std::abs(ne_after - ne_before) <= 1.0e-8 * std::max(ne_after, 1.0e-300)) break;
    }
  }
}

double ThermochemistryOperator::ionPotentialDensity(const MetalIonRowState& ions) const {
  const MetalIonLadderSolver* solver = m_solver.metalIonSolver();
  if (solver == nullptr) return 0.0;
  const MetalIonElementSet set = m_config.solver.metal_ion_elements;
  double energy = 0.0;
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    const auto element = static_cast<MetalIonElement>(e);
    if (!set.contains(element)) continue;
    energy += ions.nuclei_cm3[e] *
        solver->table().storedIonizationEnergyPerNucleusErg(
            element, std::span<const double>(ions.fractions.data() + set.stageOffset(element),
                                             metalIonStageCount(element)));
  }
  return energy;
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

namespace {

// P1-49-03: conservative atomic-reservoir <-> CO transfer of one element's
// ladder. The light network's CO chemistry (NL97 R20 C+ + O -> CO, R21
// CO + gamma -> C+ + O) moves nuclei between CO and a fixed atomic channel
// stage: C II (stage 1) for carbon, O I (stage 0) for oxygen. Formation
// removes nuclei from that stage first and only the remainder (when the
// channel stage is exhausted, i.e. the ladder and the NL97 product
// convention disagree) proportionally from all stages; destruction returns
// nuclei to the channel stage. Afterwards nuclei == target exactly and the
// fractions are renormalized (no creation/destruction of nuclei).
// Returns |delta| moved [cm^-3].
double transferLadderNuclei(std::span<double> fractions, double& nuclei, double target, std::size_t channel) {
  target = std::max(target, 0.0);
  const double delta = target - nuclei;
  if (delta == 0.0) return 0.0;
  std::array<double, 27> n{};
  for (std::size_t k = 0; k < fractions.size(); ++k) n[k] = std::max(nuclei * fractions[k], 0.0);
  if (delta > 0.0) {
    n[channel] += delta;
  } else {
    double remove = -delta;
    const double take = std::min(remove, n[channel]);
    n[channel] -= take;
    remove -= take;
    if (remove > 0.0) {
      double rest = 0.0;
      for (std::size_t k = 0; k < fractions.size(); ++k) rest += n[k];
      const double keep = rest > 0.0 ? std::max(0.0, 1.0 - remove / rest) : 0.0;
      for (std::size_t k = 0; k < fractions.size(); ++k) n[k] *= keep;
    }
  }
  double sum = 0.0;
  for (std::size_t k = 0; k < fractions.size(); ++k) sum += n[k];
  if (sum > 0.0 && std::isfinite(sum)) {
    for (std::size_t k = 0; k < fractions.size(); ++k) fractions[k] = n[k] / sum;
  }
  // target == 0: distribution is kept as the (empty) reservoir's history.
  nuclei = target;
  return std::abs(delta);
}

}  // namespace

ThermochemistryOperator::CoupledRowOutcome ThermochemistryOperator::integrateCoupled(
    const ThermochemistrySolver& light, const MetalIonLadderSolver* ion_solver, bool ions_enabled,
    const ThermochemistryCellInput& start_input, const MetalIonRowState& start_ions, double dt_seconds) const {
  CoupledRowOutcome out;
  out.ions = start_ions;
  ThermochemistryCellInput input = start_input;
  const MetalIonElementSet ion_set = m_config.solver.metal_ion_elements;
  const double rho = input.mass_density_cgs;
  const auto temperature_from = [&](const ThermochemistryCellInput& in, const NetworkState& n, double u) {
    const double ne = light.totalElectronDensity(in, n);
    const double particles = light.particleDensity(in, n, ne);
    return particles > 0.0 ? (m_config.solver.adiabatic_index - 1.0) * u * rho / (k_boltzmann_erg_per_k * particles)
                           : m_config.solver.temperature_floor_k;
  };
  // Operator-split half of the ladder at the frozen (T, n_e, radiation)
  // context of `target` (P1-49-02: one authoritative context).
  const auto advance_ions = [&](ThermochemistryCellInput& target, double temperature_k) -> bool {
    const MetalIonEnvironment env = light.metalIonEnvironment(target, target.species_cm3, temperature_k);
    for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
      const auto element = static_cast<MetalIonElement>(e);
      if (!ion_set.contains(element)) continue;
      const std::span<double> fractions(out.ions.fractions.data() + ion_set.stageOffset(element),
                                        metalIonStageCount(element));
      if (ion_solver->advance(element, fractions, env, 0.5 * dt_seconds, out.ion_report) !=
          SourceAdvanceStatus::kCompleted) {
        return false;
      }
    }
    refreshMetalIonInput(target, out.ions);
    return true;
  };
  const auto fail_ions = [&](const char* phase) {
    out.status = SourceAdvanceStatus::kRetryShorterInterval;
    out.suggested_retry_seconds = 0.5 * dt_seconds;
    out.diagnostic = std::string("metal-ion ladder ") + phase + " half did not converge within the substep cap";
    return out;
  };
  const double potential_begin = ions_enabled ? ionPotentialDensity(out.ions) : 0.0;
  try {
    if (ions_enabled) {
      const double t_start = temperature_from(input, input.species_cm3, input.specific_internal_energy_erg_g);
      if (!advance_ions(input, std::max(t_start, m_config.solver.temperature_floor_k))) return fail_ions("first");
    }
    light.advanceCell(input, dt_seconds, out.output, out.ledger);
    if (out.output.status != SourceAdvanceStatus::kCompleted) {
      out.status = out.output.status;
      out.suggested_retry_seconds = out.output.suggested_retry_seconds;
      out.diagnostic = out.output.diagnostic;
      return out;
    }
    if (ions_enabled) {
      ThermochemistryCellInput end_input = input;
      end_input.species_cm3 = out.output.species_cm3;
      end_input.specific_internal_energy_erg_g = out.output.specific_internal_energy_erg_g;
      // P1-49-03: reconcile C/O ladder nuclei with the END CO abundance of
      // the light network before the closing half (conservative transfer).
      if (light.speciesCount() > static_cast<std::size_t>(ChemicalSpecies::kCO)) {
        const double n_co_end = out.output.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kCO)];
        for (const auto element : {MetalIonElement::kC, MetalIonElement::kO}) {
          if (!ion_set.contains(element)) continue;
          const std::size_t e = static_cast<std::size_t>(element);
          const ChemicalElement chemical = chemicalElementOf(element);
          const double element_nuclei =
              end_input.gas_element_mass_fraction[static_cast<std::size_t>(chemical)] * rho /
              (chemicalElementInfo(chemical).atomic_mass_amu * k_amu_g);
          out.ledger.metal_ion_co_transfer_nuclei += transferLadderNuclei(
              std::span<double>(out.ions.fractions.data() + ion_set.stageOffset(element),
                                metalIonStageCount(element)),
              out.ions.nuclei_cm3[e], element_nuclei - n_co_end,
              element == MetalIonElement::kC ? 1U : 0U);
        }
      }
      refreshMetalIonInput(end_input, out.ions);
      const double t_end = temperature_from(end_input, end_input.species_cm3,
                                            end_input.specific_internal_energy_erg_g);
      if (!advance_ions(end_input, std::max(t_end, m_config.solver.temperature_floor_k))) return fail_ions("closing");
      // P1-49-05: EOS refresh at the AUTHORITATIVE end thermal energy u_end
      // (unchanged: no energy is manufactured) with the final ion
      // composition. Conserved: u (thermal), element nuclei, CO/atomic
      // partition. Derived: n_e, mu, T (P and c_s follow from u and are
      // refreshed by the caller for every committed chemistry row).
      const double ne_final = light.totalElectronDensity(end_input, out.output.species_cm3);
      const double particles = light.particleDensity(end_input, out.output.species_cm3, ne_final);
      out.output.mean_molecular_weight = ThermochemistrySolver::meanMolecularWeight(rho, particles);
      out.output.temperature_k = temperature_from(end_input, out.output.species_cm3,
                                                  out.output.specific_internal_energy_erg_g);
      const double n_h = end_input.gas_element_mass_fraction[0] * rho /
          (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
      out.output.electron_per_hydrogen = n_h > 0.0 ? ne_final / n_h : 0.0;
      out.ion_electron_cm3 = end_input.metal_ion_electron_cm3;
      // P1-49-04: closing charge check on the identical inventory.
      double charge = end_input.metal_ion_electron_cm3;
      const std::size_t stored = std::max<std::size_t>(light.speciesCount(), k_hybrid_species_count);
      for (std::size_t sp = 0; sp < stored; ++sp) {
        charge += chemicalSpeciesInfo(static_cast<ChemicalSpecies>(sp)).charge * out.output.species_cm3[sp];
      }
      if (n_h > 0.0) {
        out.ledger.charge_residual =
            std::max(out.ledger.charge_residual, std::abs(std::max(charge, 0.0) - ne_final) / n_h);
      }
      out.ledger.accepted_substeps += out.ion_report.accepted_substeps;
      out.ledger.rejected_substeps += out.ion_report.rejected_substeps;
      out.ledger.max_error_estimate = std::max(out.ledger.max_error_estimate, out.ion_report.max_error_estimate);
      out.ledger.metal_ion_normalization_residual =
          std::max(out.ledger.metal_ion_normalization_residual, out.ion_report.normalization_residual);
      out.ledger.metal_table_clamps += out.ion_report.range_clamps;
      // P1-49-06: stored-potential partition (see SourceLedger).
      out.ion_potential_erg_cm3 = ionPotentialDensity(out.ions);
      out.ledger.metal_ion_potential_delta_erg += out.ion_potential_erg_cm3 - potential_begin;
    }
  } catch (const std::range_error& error) {
    out.status = SourceAdvanceStatus::kRejectedInvalid;
    out.diagnostic = std::string("metal-ion dataset range: ") + error.what();
  }
  return out;
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
  // ions(dt/2) at the step-start (T, n_e, radiation) -> light network +
  // energy over dt with the ion distribution frozen (ion losses evaluated
  // implicitly at the Newton temperature) -> C/O <-> CO reservoir transfer
  // -> ions(dt/2) at the end (T, n_e, radiation) -> EOS refresh.
  const bool ions = !demoted && rowCarriesMetalIons(chemistry_tag);
  MetalIonRowState start_ions;
  MetalIonAdvanceReport attach_report;
  const MetalIonElementSet ion_set = m_config.solver.metal_ion_elements;
  if (ions) {
    try {
      attachMetalIonState(state, row, rowTemperatureK(state, row, input), input, start_ions, attach_report);
    } catch (const std::range_error& error) {
      result.status = SourceAdvanceStatus::kRejectedInvalid;
      result.diagnostic = "row " + std::to_string(row) + ": metal-ion dataset range: " + error.what();
      return result;
    }
  }
  // `input` (with the attached step-start ladder) is the shared start state
  // of production and shadow reference.
  CoupledRowOutcome production =
      integrateCoupled(solver, m_solver.metalIonSolver(), ions, input, start_ions, interval.physical_seconds);
  if (production.status != SourceAdvanceStatus::kCompleted) {
    result.status = production.status;
    result.suggested_retry_seconds = production.suggested_retry_seconds;
    result.diagnostic = "row " + std::to_string(row) + " (gas_cell_id " +
        std::to_string(state.gas_cells.gas_cell_id[row]) + "): " + production.diagnostic;
    return result;
  }
  production.ledger.metal_table_clamps += attach_report.range_clamps;
  const ThermochemistryCellOutput& output = production.output;
  SourceLedger& per_volume = production.ledger;
  const MetalIonRowState& ion_state = production.ions;
  if (m_reference_solver != nullptr && shadow.samples < m_config.shadow.max_cells) {
    // Deterministic sample keyed by stable id and accepted step.
    const std::uint64_t key = core::fnv1a64(std::to_string(state.gas_cells.gas_cell_id[row]) + ":" +
                                            std::to_string(interval.accepted_step_index));
    const double u01 = static_cast<double>(key >> 11U) * (1.0 / 9007199254740992.0);
    if (u01 < m_config.shadow.sample_fraction) {
      // P1-49-07 / patch 51: the reference ALWAYS runs the full detailed
      // coupled chronology (its own higher-accuracy light network AND, when
      // ladders are configured, its own higher-accuracy ladder:
      // m_reference_solver->metalIonSolver(), 10x tighter tolerance, 4x
      // substep cap), from the identical step-start state, on private stack
      // copies only (no production state is touched; no heap allocation, so
      // no governor admission is required; cost bounded by shadow.max_cells).
      // For a demoted (hybrid) row the production branch stays genuinely
      // reduced; the reference reconstructs the row's detailed ladder
      // privately from the authoritative element lanes minus CO and the
      // frozen stored stage distribution (or, for a row that never carried
      // ladder state, the documented ionization-equilibrium projection at
      // the step-start (T, n_e, radiation) - the same initializer a promoted
      // row would use). This measures exactly the error of demotion.
      const bool reference_ions = m_reference_solver->metalIonSolver() != nullptr;
      ThermochemistryCellInput reference_input = input;
      MetalIonRowState reference_start = start_ions;
      bool reference_ready = true;
      if (reference_ions && !ions) {
        MetalIonAdvanceReport reference_attach;
        try {
          attachMetalIonState(state, row, rowTemperatureK(state, row, reference_input), reference_input,
                              reference_start, reference_attach);
        } catch (const std::range_error&) {
          reference_ready = false;  // reference outside its dataset: no sample (production unaffected)
        }
      }
      const CoupledRowOutcome reference = reference_ready
          ? integrateCoupled(*m_reference_solver, reference_ions ? m_reference_solver->metalIonSolver() : nullptr,
                             reference_ions, reference_input, reference_start, interval.physical_seconds)
          : CoupledRowOutcome{.status = SourceAdvanceStatus::kRejectedInvalid};
      if (reference.status == SourceAdvanceStatus::kCompleted) {
        ++shadow.samples;
        NetworkState prod_species = output.species_cm3;
        if (demoted) {
          for (std::size_t sp = solver.speciesCount(); sp < m_solver.speciesCount(); ++sp) {
            prod_species[sp] = input.species_cm3[sp];  // frozen reservoir
          }
        }
        double species = 0.0;
        for (std::size_t sp = 0; sp < m_solver.speciesCount(); ++sp) {
          const double scale = std::max({std::abs(reference.output.species_cm3[sp]), std::abs(prod_species[sp]),
                                         1.0e-10 * std::max(input.species_cm3[0], 1.0e-300)});
          species = std::max(species, std::abs(reference.output.species_cm3[sp] - prod_species[sp]) / scale);
        }
        shadow.max_relative_species_difference = std::max(shadow.max_relative_species_difference, species);
        shadow.max_relative_electron_difference = std::max(
            shadow.max_relative_electron_difference,
            std::abs(reference.output.electron_per_hydrogen - output.electron_per_hydrogen) /
                std::max(reference.output.electron_per_hydrogen, 1.0e-12));
        shadow.max_relative_energy_difference = std::max(
            shadow.max_relative_energy_difference,
            std::abs(reference.output.specific_internal_energy_erg_g - output.specific_internal_energy_erg_g) /
                reference.output.specific_internal_energy_erg_g);
        if (reference_ions) {
          ++shadow.ion_samples;
          // Production ladder: evolved (detailed row) or the frozen /
          // projected distribution a reduced row carries (hybrid row).
          const MetalIonRowState& production_ladder = ions ? ion_state : reference_start;
          const auto metal_free_electrons = [&](const MetalIonRowState& ladder) {
            double electrons = 0.0;
            for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
              const auto element = static_cast<MetalIonElement>(e);
              if (!ion_set.contains(element)) continue;
              electrons += ladder.nuclei_cm3[e] *
                  metalIonElectronsPerNucleus(std::span<const double>(
                      ladder.fractions.data() + ion_set.stageOffset(element), metalIonStageCount(element)));
            }
            return electrons;
          };
          double stage = 0.0;
          for (std::size_t k = 0; k < ion_set.stageCount(); ++k) {
            stage = std::max(stage, std::abs(reference.ions.fractions[k] - production_ladder.fractions[k]));
          }
          shadow.max_ion_stage_fraction_difference = std::max(shadow.max_ion_stage_fraction_difference, stage);
          // Metal free electrons the production model actually counts: the
          // full ladder for a detailed row; none for a reduced row (the
          // hybrid network carries no metal-ion electrons). Normalized by the
          // reference TOTAL electron density, i.e. the n_e error of the model.
          const double n_h = input.gas_element_mass_fraction[0] * input.mass_density_cgs /
              (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
          const double production_metal_e = ions ? metal_free_electrons(ion_state) : 0.0;
          const double reference_metal_e = metal_free_electrons(reference.ions);
          const double ne_reference_total = std::max(reference.output.electron_per_hydrogen * n_h, 1.0e-300);
          shadow.max_relative_ion_electron_difference = std::max(
              shadow.max_relative_ion_electron_difference,
              std::abs(reference_metal_e - production_metal_e) / ne_reference_total);
          shadow.max_relative_ion_temperature_difference = std::max(
              shadow.max_relative_ion_temperature_difference,
              std::abs(reference.output.temperature_k - output.temperature_k) /
                  std::max(reference.output.temperature_k, 1.0e-300));
          const double production_potential = ions ? production.ion_potential_erg_cm3
                                                   : ionPotentialDensity(reference_start);
          const double e_scale = std::max({std::abs(reference.ion_potential_erg_cm3),
                                           std::abs(production_potential), 1.0e-300});
          shadow.max_relative_ion_potential_difference = std::max(
              shadow.max_relative_ion_potential_difference,
              std::abs(reference.ion_potential_erg_cm3 - production_potential) / e_scale);
        }
        const double ratio =
            static_cast<double>(reference.ledger.accepted_substeps + reference.ledger.rejected_substeps) /
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
  result.ledger.metal_ion_potential_delta_erg *= volume_cm3;
  result.ledger.metal_ion_co_transfer_nuclei *= volume_cm3;

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
    // the END authoritative nuclei (element lane minus END CO, P1-49-03).
    staged.ion_stage_count = static_cast<std::uint32_t>(ion_set.stageCount());
    staged.ion_state_initialized_from_equilibrium = ion_state.initialized_from_equilibrium;
    const double n_co_end = m_solver.speciesCount() > static_cast<std::size_t>(ChemicalSpecies::kCO)
        ? output.species_cm3[static_cast<std::size_t>(ChemicalSpecies::kCO)]
        : 0.0;
    for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
      const auto element = static_cast<MetalIonElement>(e);
      if (!ion_set.contains(element)) continue;
      const ChemicalElement chemical = chemicalElementOf(element);
      const double amu = chemicalElementInfo(chemical).atomic_mass_amu;
      const std::size_t offset = ion_set.stageOffset(element);
      double staged_nuclei = 0.0;
      for (std::size_t k = 0; k < metalIonStageCount(element); ++k) {
        staged.ion_stage_mass_code[offset + k] = ion_state.nuclei_cm3[e] * ion_state.fractions[offset + k] *
            amu * k_amu_g / input.mass_density_cgs * mass;
        staged_nuclei += ion_state.nuclei_cm3[e] * ion_state.fractions[offset + k];
      }
      double expected = input.gas_element_mass_fraction[static_cast<std::size_t>(chemical)] *
          input.mass_density_cgs / (amu * k_amu_g);
      if (chemical == ChemicalElement::kC || chemical == ChemicalElement::kO) expected -= n_co_end;
      expected = std::max(expected, 0.0);
      if (expected > 0.0) {
        result.ledger.metal_ion_inventory_residual =
            std::max(result.ledger.metal_ion_inventory_residual, std::abs(staged_nuclei - expected) / expected);
      }
    }
  }
  if (!std::isfinite(staged.internal_energy_code) || !(staged.internal_energy_code > 0.0)) {
    result.status = SourceAdvanceStatus::kRejectedInvalid;
    result.diagnostic = "thermochemistry produced a non-positive internal energy at row " + std::to_string(row);
  }
  return result;
}

double ThermochemistryOperator::reconcileMetalIonStageLanes(core::GasCompositionSidecar& sidecar,
                                                            std::size_t row) const {
  if (m_layout.ionStageCount() == 0U) return 0.0;
  const MetalIonElementSet set = m_config.solver.metal_ion_elements;
  const bool has_co = m_solver.speciesCount() > static_cast<std::size_t>(ChemicalSpecies::kCO);
  const double co_mass = has_co ? std::max(sidecar.value(m_layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCO)), row), 0.0)
                                : 0.0;
  const double co_amu = chemicalSpeciesInfo(ChemicalSpecies::kCO).mass_amu;
  double max_correction = 0.0;
  for (std::size_t e = 0; e < k_metal_ion_element_count; ++e) {
    const auto element = static_cast<MetalIonElement>(e);
    if (!set.contains(element)) continue;
    const ChemicalElement chemical = chemicalElementOf(element);
    double target = std::max(sidecar.value(m_layout.elementLane(chemical), row), 0.0);
    if (chemical == ChemicalElement::kC || chemical == ChemicalElement::kO) {
      target -= co_mass * chemicalElementInfo(chemical).atomic_mass_amu / co_amu;
    }
    target = std::max(target, 0.0);
    const std::size_t offset = set.stageOffset(element);
    double sum = 0.0;
    for (std::size_t k = 0; k < metalIonStageCount(element); ++k) {
      sum += std::max(sidecar.value(m_layout.ionStageLane(offset + k), row), 0.0);
    }
    if (!(sum > 0.0) || !std::isfinite(sum)) continue;
    const double scale = target / sum;
    max_correction = std::max(max_correction, std::abs(scale - 1.0));
    for (std::size_t k = 0; k < metalIonStageCount(element); ++k) {
      const std::uint32_t lane = m_layout.ionStageLane(offset + k);
      sidecar.setValue(lane, row, std::max(sidecar.value(lane, row), 0.0) * scale);
    }
  }
  return max_correction;
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
  // Patch 51 (A03/C02): the cooling estimate uses the SAME model that
  // advances this row (fidelityModelForRow): a reduced row is costed by the
  // hybrid solver, whose element table / fine-structure terms still include
  // C/O/Si/Fe (the hybrid solver owns no ladder, so nothing is excluded); a
  // detailed row by the detailed solver with its attached ladder (ladder
  // losses replace the excluded table elements). Never both, never neither.
  const ThermochemistryCellInput input = cellInput(state, row, interval);
  const ThermochemistrySolver& row_solver = fidelityModelForRow(state, row);
  double mu = state.gas_composition.mean_molecular_weight[row];
  if (!(mu > 0.0) || !std::isfinite(mu)) mu = 1.22;
  const double temperature_k = (m_config.solver.adiabatic_index - 1.0) * mu * k_amu_g *
      input.specific_internal_energy_erg_g / k_boltzmann_erg_per_k;
  if (!(temperature_k > 0.0) || !(input.mass_density_cgs > 0.0) ||
      !(input.specific_internal_energy_erg_g > 0.0)) {
    return std::numeric_limits<double>::infinity();
  }
  const CoolingHeatingTerms terms =
      row_solver.coolingHeating(input, input.species_cm3, temperature_k);
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
    const CoolingHeatingTerms terms = fidelityModelForRow(state, row).coolingHeating(input, input.species_cm3, t);
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
