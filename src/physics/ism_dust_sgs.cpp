#include "cosmosim/physics/ism_dust_sgs.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {
namespace {
constexpr double k_newton_g_cgs = 6.67430e-8;
constexpr double k_pi = 3.14159265358979323846;
}  // namespace

// --------------------------------------------------------------------------
// ISM regimes
// --------------------------------------------------------------------------
double jeansLengthCm(double temperature_k, double mean_molecular_weight, double rho_cgs, double adiabatic_index,
                     double extra_velocity_dispersion_sq_cm2_s2) {
  if (!(rho_cgs > 0.0) || !(mean_molecular_weight > 0.0)) return std::numeric_limits<double>::infinity();
  const double cs2 = adiabatic_index * k_boltzmann_erg_per_k * std::max(temperature_k, 0.0) /
          (mean_molecular_weight * k_amu_g) + std::max(extra_velocity_dispersion_sq_cm2_s2, 0.0);
  return std::sqrt(k_pi * cs2 / (k_newton_g_cgs * rho_cgs));
}

IsmRegimeDecision decideIsmRegime(const core::SimulationState& state, const IsmRegimeConfig& config,
                                  const ThermochemistryUnitBridge& units, const SourceInterval& interval,
                                  std::uint32_t row) {
  const auto& sidecar = state.gas_composition;
  IsmRegimeDecision decision;
  decision.tag = sidecar.ism_regime_tag[row];
  decision.dwell_years = sidecar.ism_regime_dwell_years[row] + interval.physical_seconds / k_seconds_per_year;
  const double a = interval.evaluation_scale_factor;
  const double frame_density = units.comoving ? 1.0 / (a * a * a) : 1.0;
  const double frame_length = units.comoving ? a : 1.0;
  const double h = config.hysteresis_fraction;
  const auto& descriptor = sidecar.descriptor();
  const double rho_code = state.gas_cells.density_code[row];
  const double rho_cgs = rho_code * frame_density * units.density_code_to_cgs;
  double mass = 0.0;
  for (std::uint32_t l = 0; l < descriptor.element_lane_count; ++l) mass += sidecar.value(descriptor.element_lane_begin + l, row);
  for (std::uint32_t l = 0; l < descriptor.dust_lane_count; ++l) mass += sidecar.value(descriptor.dust_lane_begin + l, row);
  const double x_h = mass > 0.0 ? sidecar.value(descriptor.element_lane_begin, row) / mass : 0.75;
  const double n_h = x_h * rho_cgs / (chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g);
  double mu = sidecar.mean_molecular_weight[row];
  if (!(mu > 0.0) || !std::isfinite(mu)) mu = 1.22;  // neutral reference until first chemistry commit
  const double u_cgs = state.gas_cells.internal_energy_code[row] * units.velocity_code_to_cm_s * units.velocity_code_to_cm_s;
  const double t_k = (config.adiabatic_index - 1.0) * mu * k_amu_g * u_cgs / k_boltzmann_erg_per_k;
  const double volume_code = state.cells.mass_code[row] / std::max(rho_code, 1.0e-300);
  const double dx_cm = std::cbrt(std::max(volume_code, 0.0)) * frame_length * units.length_code_to_cm;
  const double lambda = jeansLengthCm(t_k, mu, rho_cgs, config.adiabatic_index);
  const double required = config.jeans_resolution_cells * dx_cm;
  std::uint8_t desired = decision.tag;
  if (decision.tag == core::k_gas_ism_regime_effective) {
    if (lambda > required * (1.0 + h) || n_h < config.effective_threshold_n_h_cm3 * (1.0 - h)) {
      desired = core::k_gas_ism_regime_resolved;
    }
  } else if (lambda < required * (1.0 - h) && n_h >= config.effective_threshold_n_h_cm3 * (1.0 + h)) {
    desired = core::k_gas_ism_regime_effective;
  }
  if (desired != decision.tag) {
    if (decision.dwell_years < config.minimum_dwell_years) {
      decision.held_by_dwell = true;
    } else {
      decision.tag = desired;
      decision.dwell_years = 0.0;
      decision.changed = true;
    }
  }
  return decision;
}

IsmRegimeReport updateIsmRegimes(core::SimulationState& state, const IsmRegimeConfig& config,
                                 const ThermochemistryUnitBridge& units, const SourceInterval& interval,
                                 std::span<const std::uint32_t> rows) {
  IsmRegimeReport report;
  auto& sidecar = state.gas_composition;
  if (!sidecar.enabled()) {
    throw std::runtime_error("ISM regime manager requires the gas composition sidecar");
  }
  sidecar.requireAlignedWith(state.gas_cells.gas_cell_id, "updateIsmRegimes");
  const std::size_t count = rows.empty() ? state.cells.size() : rows.size();
  // Decide for all rows from the unmodified state first, then apply.
  std::vector<std::pair<std::uint32_t, IsmRegimeDecision>> decisions;
  decisions.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const std::uint32_t row = rows.empty() ? static_cast<std::uint32_t>(i) : rows[i];
    decisions.emplace_back(row, decideIsmRegime(state, config, units, interval, row));
  }
  for (const auto& [row, d] : decisions) {
    sidecar.ism_regime_tag[row] = d.tag;
    sidecar.ism_regime_dwell_years[row] = d.dwell_years;
    if (d.held_by_dwell) ++report.held_by_dwell;
    if (d.changed) {
      report.changed_rows.push_back(row);
      if (d.tag == core::k_gas_ism_regime_effective) ++report.to_effective; else ++report.to_resolved;
    }
    if (d.tag == core::k_gas_ism_regime_effective) ++report.effective_cells; else ++report.resolved_cells;
  }
  return report;
}

// --------------------------------------------------------------------------
// Dust
// --------------------------------------------------------------------------
namespace {
// Largest dust mass of reservoir d that the gas element lanes can still
// supply, limited by the key element of the stoichiometry.
double condensableDustMass(const core::GasCompositionSidecar& sidecar, const CompositionLayout& layout,
                           std::size_t row, std::size_t d) {
  const auto& fractions = dustElementMassFractions(static_cast<DustReservoir>(d));
  const bool detailed = sidecar.descriptor().species_lane_count >= k_detailed_species_count;
  double co_mass = 0.0, cii_mass = 0.0;
  if (detailed) {
    co_mass = sidecar.value(layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCO)), row);
    cii_mass = sidecar.value(layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCII)), row);
  }
  const double a_c = chemicalElementInfo(ChemicalElement::kC).atomic_mass_amu;
  const double a_o = chemicalElementInfo(ChemicalElement::kO).atomic_mass_amu;
  double limit = std::numeric_limits<double>::infinity();
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    if (fractions[e] <= 0.0) continue;
    double available = sidecar.value(layout.elementLane(static_cast<ChemicalElement>(e)), row);
    // Molecule-bound nuclei cannot condense: CO keeps its C and O.
    if (detailed && e == static_cast<std::size_t>(ChemicalElement::kC)) available = cii_mass;
    if (detailed && e == static_cast<std::size_t>(ChemicalElement::kO)) available -= co_mass * a_o / (a_c + a_o);
    limit = std::min(limit, std::max(available, 0.0) / fractions[e]);
  }
  return std::isfinite(limit) ? std::max(limit, 0.0) : 0.0;
}

void transferGasToDust(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout, std::size_t row,
                       std::size_t d, double dust_mass_delta) {
  // dust_mass_delta > 0: gas -> dust; < 0: dust -> gas. Stoichiometric.
  const auto& fractions = dustElementMassFractions(static_cast<DustReservoir>(d));
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    if (fractions[e] <= 0.0) continue;
    const std::uint32_t lane = layout.elementLane(static_cast<ChemicalElement>(e));
    sidecar.setValue(lane, row, std::max(sidecar.value(lane, row) - fractions[e] * dust_mass_delta, 0.0));
  }
  const std::uint32_t dust_lane = layout.dustLane(d);
  sidecar.setValue(dust_lane, row, std::max(sidecar.value(dust_lane, row) + dust_mass_delta, 0.0));
  // Species lanes of the networks (H/He/C/O) are untouched by refractory
  // transfers except carbon: C+ (and CO) track gas-phase carbon in the
  // detailed network, so the carbon removed/returned is applied to C+.
  if (fractions[static_cast<std::size_t>(ChemicalElement::kC)] > 0.0 &&
      sidecar.descriptor().species_lane_count >= k_detailed_species_count) {
    const std::uint32_t cii = layout.speciesLane(static_cast<std::size_t>(ChemicalSpecies::kCII));
    const double delta_c = fractions[static_cast<std::size_t>(ChemicalElement::kC)] * dust_mass_delta;
    sidecar.setValue(cii, row, std::max(sidecar.value(cii, row) - delta_c, 0.0));
  }
}
}  // namespace

void advanceDustCell(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout, std::size_t row,
                     double n_h_cm3, double temperature_k, double dt_seconds, const DustConfig& config,
                     DustUpdateLedger& ledger) {
  if (sidecar.descriptor().dust_lane_count == 0U || !(dt_seconds > 0.0)) return;
  double metals = 0.0, mass = 0.0;
  for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
    const double m = sidecar.value(layout.elementLane(static_cast<ChemicalElement>(e)), row);
    mass += m;
    if (e >= 2) metals += m;
  }
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    const double dust = sidecar.value(layout.dustLane(d), row);
    mass += dust;
    metals += dust;
  }
  if (!(mass > 0.0)) return;
  const double z_rel = (metals / mass) / solarMetalMassFraction();
  const double t = std::max(temperature_k, 1.0);
  const double n_h = std::max(n_h_cm3, 0.0);
  const double tau_sp = config.sputtering_timescale_ref_s / std::max(n_h / 1.0e-3, 1.0e-300) /
      (std::pow(2.0e6 / t, 2.5) + 1.0);
  const double tau_acc = (z_rel > 0.0 && n_h > 0.0)
      ? config.accretion_timescale_ref_s / ((n_h / 100.0) * std::sqrt(t / 50.0) * z_rel)
      : std::numeric_limits<double>::infinity();
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    const double d0 = sidecar.value(layout.dustLane(d), row);
    if (!(d0 > 0.0)) continue;  // growth needs seed grains; sputtering needs dust
    // Operator split within the cell: sputtering (exact) then accretion (exact).
    const double after_sputter = d0 * std::exp(-3.0 * dt_seconds / tau_sp);
    transferGasToDust(sidecar, layout, row, d, after_sputter - d0);
    ledger.dust_to_gas_code += d0 - after_sputter;
    if (std::isfinite(tau_acc) && after_sputter > 0.0) {
      const double d_max = after_sputter + condensableDustMass(sidecar, layout, row, d);
      const double grown = d_max / (1.0 + (d_max / after_sputter - 1.0) * std::exp(-dt_seconds / tau_acc));
      const double delta = std::clamp(grown - after_sputter, 0.0, d_max - after_sputter);
      transferGasToDust(sidecar, layout, row, d, delta);
      ledger.gas_to_dust_code += delta;
    }
  }
  ++ledger.cells;
}

double destroyDustBySupernova(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout,
                              std::size_t row, double swept_gas_mass_code, double cell_gas_mass_code) {
  if (sidecar.descriptor().dust_lane_count == 0U || !(cell_gas_mass_code > 0.0)) return 0.0;
  const double fraction = std::clamp(swept_gas_mass_code / cell_gas_mass_code, 0.0, 1.0);
  double destroyed_total = 0.0;
  for (std::size_t d = 0; d < sidecar.descriptor().dust_lane_count; ++d) {
    const double destroyed = fraction * sidecar.value(layout.dustLane(d), row);
    transferGasToDust(sidecar, layout, row, d, -destroyed);
    destroyed_total += destroyed;
  }
  return destroyed_total;
}

EjectaPartition partitionEjectaIntoDust(const std::array<double, k_chemical_element_count>& ejecta,
                                        double condensation_efficiency) {
  EjectaPartition result;
  result.gas_elements = ejecta;
  const double delta = std::clamp(condensation_efficiency, 0.0, 1.0);
  if (delta == 0.0) return result;
  for (std::size_t d = 0; d < k_dust_reservoir_count; ++d) {
    const auto& fractions = dustElementMassFractions(static_cast<DustReservoir>(d));
    double limit = std::numeric_limits<double>::infinity();
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
      if (fractions[e] > 0.0) limit = std::min(limit, result.gas_elements[e] / fractions[e]);
    }
    const double dust = std::isfinite(limit) ? delta * std::max(limit, 0.0) : 0.0;
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
      result.gas_elements[e] = std::max(result.gas_elements[e] - fractions[e] * dust, 0.0);
    }
    result.dust[d] = dust;
  }
  return result;
}

// --------------------------------------------------------------------------
// SGS
// --------------------------------------------------------------------------
double sgsSpecificEnergyAfter(double k0, double dt, double dx, double c_eps) noexcept {
  if (!(k0 > 0.0) || !(dt > 0.0) || !(dx > 0.0)) return std::max(k0, 0.0);
  const double denom = 1.0 + c_eps * std::sqrt(k0) * dt / (2.0 * dx);
  return k0 / (denom * denom);
}

SgsDissipation sgsDissipate(double sgs_energy, double mass, double dt_code_time, double cell_size_code,
                            const SgsConfig& config) noexcept {
  SgsDissipation out;
  out.sgs_energy = std::max(sgs_energy, 0.0);
  if (!(mass > 0.0)) return out;
  const double k0 = out.sgs_energy / mass;
  const double k1 = sgsSpecificEnergyAfter(k0, dt_code_time, cell_size_code, config.dissipation_coefficient);
  out.sgs_energy = mass * k1;
  out.dissipated = mass * (k0 - k1);
  return out;
}

double advanceSgsCell(core::SimulationState& state, const CompositionLayout& layout, std::size_t row,
                      double dt_code_time, double cell_size_code, const SgsConfig& config, SgsLedger& ledger) {
  auto& sidecar = state.gas_composition;
  if (!layout.hasSgs()) return 0.0;
  const std::uint32_t lane = sidecar.descriptor().sgs_energy_lane;
  const double mass = state.cells.mass_code[row];
  if (!(mass > 0.0)) return 0.0;
  const double e_sgs = sidecar.value(lane, row);
  const double k0 = e_sgs / mass;
  const double k1 = sgsSpecificEnergyAfter(k0, dt_code_time, cell_size_code, config.dissipation_coefficient);
  const double dissipated = mass * (k0 - k1);
  sidecar.setValue(lane, row, mass * k1);
  state.gas_cells.internal_energy_code[row] += dissipated / mass;  // heat, exact transfer
  ledger.dissipated_to_heat_code += dissipated;
  return dissipated;
}

}  // namespace cosmosim::physics
