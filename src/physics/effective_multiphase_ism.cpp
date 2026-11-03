#include "cosmosim/physics/effective_multiphase_ism.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cosmosim::physics {
namespace {

constexpr double k_proton_mass_kg = 1.67262192369e-27;
constexpr double k_boltzmann_j_per_k = 1.380649e-23;
constexpr double k_boltzmann_erg_per_k = 1.380649e-16;
constexpr double k_small = 1.0e-30;

[[nodiscard]] double stableColdFraction(double y) {
  if (!std::isfinite(y) || y <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Algebraically equivalent to 1 + 1/(2y) - sqrt(1/y + 1/(4y^2)),
  // but avoids catastrophic cancellation for y << 1.
  const double root = std::sqrt(1.0 + 4.0 * y);
  return std::clamp((root - 1.0) / (root + 1.0), 0.0, 1.0);
}

[[nodiscard]] std::uint64_t fnvMix(std::uint64_t hash, std::uint64_t value) {
  for (int byte_index = 0; byte_index < 8; ++byte_index) {
    hash ^= (value >> (8 * byte_index)) & 0xffULL;
    hash *= 1099511628211ULL;
  }
  return hash;
}

[[nodiscard]] std::uint64_t hashDouble(std::uint64_t hash, double value) {
  return fnvMix(hash, std::bit_cast<std::uint64_t>(value));
}

}  // namespace

EffectiveMultiphaseEosTable::EffectiveMultiphaseEosTable(
    EffectiveMultiphaseEosConfig config,
    core::UnitSystem units,
    CoolingRateProvider cooling_provider)
    : m_config(std::move(config)),
      m_units(std::move(units)),
      m_cooling_provider(std::move(cooling_provider)) {
  if (m_config.n_h_threshold_cgs <= 0.0 || m_config.hydrogen_mass_fraction <= 0.0 ||
      m_config.hydrogen_mass_fraction > 1.0 || m_config.t_star_at_threshold_code <= 0.0 ||
      m_config.evaporation_factor_at_threshold <= 0.0 || m_config.q_eos < 0.0 ||
      m_config.q_eos > 1.0 || m_config.massive_star_fraction < 0.0 ||
      m_config.massive_star_fraction >= 1.0 || m_config.table_bins < 8U ||
      m_config.max_density_ratio <= 1.0 || m_config.adiabatic_index <= 1.0) {
    throw std::invalid_argument("EffectiveMultiphaseEosTable received invalid configuration");
  }
  m_threshold_density_cgs =
      m_config.n_h_threshold_cgs * (k_proton_mass_kg * 1.0e3) /
      m_config.hydrogen_mass_fraction;
  m_threshold_density_phys_code = m_units.densitySiToCode(m_threshold_density_cgs * 1.0e3);
  if (!std::isfinite(m_threshold_density_phys_code) || m_threshold_density_phys_code <= 0.0) {
    throw std::invalid_argument("EffectiveMultiphaseEosTable threshold conversion is invalid");
  }
  build();
}

const EffectiveMultiphaseEosConfig& EffectiveMultiphaseEosTable::config() const noexcept {
  return m_config;
}

const core::UnitSystem& EffectiveMultiphaseEosTable::units() const noexcept { return m_units; }

double EffectiveMultiphaseEosTable::thresholdDensityPhysCode() const noexcept {
  return m_threshold_density_phys_code;
}

double EffectiveMultiphaseEosTable::thresholdDensityCgs() const noexcept {
  return m_threshold_density_cgs;
}

double EffectiveMultiphaseEosTable::specificInternalEnergyFromTemperatureCode(
    double temperature_k,
    double mean_molecular_weight) const {
  if (temperature_k <= 0.0 || mean_molecular_weight <= 0.0) {
    throw std::invalid_argument("effective EOS temperature conversion requires positive inputs");
  }
  const double u_si = k_boltzmann_j_per_k * temperature_k /
      ((m_config.adiabatic_index - 1.0) * mean_molecular_weight * k_proton_mass_kg);
  return u_si / (m_units.velocity_si_per_code * m_units.velocity_si_per_code);
}

EffectiveMultiphaseEosEntry EffectiveMultiphaseEosTable::evaluateDirect(
    double density_phys_code) const {
  EffectiveMultiphaseEosEntry entry;
  entry.density_phys_code = density_phys_code;
  entry.density_ratio = density_phys_code / m_threshold_density_phys_code;
  if (!std::isfinite(entry.density_ratio) || entry.density_ratio < 1.0) {
    return entry;
  }

  const double u_cold = specificInternalEnergyFromTemperatureCode(
      m_config.cold_phase_temperature_k, m_config.mean_molecular_weight_cold);
  const double u_iso = specificInternalEnergyFromTemperatureCode(
      m_config.isothermal_temperature_k, m_config.mean_molecular_weight_ionized);
  entry.star_formation_timescale_code = m_config.t_star_at_threshold_code /
      std::sqrt(entry.density_ratio);
  const double evaporation = m_config.evaporation_factor_at_threshold *
      std::pow(entry.density_ratio, m_config.evaporation_density_exponent);
  entry.specific_internal_energy_hot_code = u_cold +
      m_config.supernova_specific_energy_code / (1.0 + evaporation);

  const double hot_temperature_k = std::max(
      m_config.cold_phase_temperature_k,
      entry.specific_internal_energy_hot_code *
          m_units.velocity_si_per_code * m_units.velocity_si_per_code *
          (m_config.adiabatic_index - 1.0) *
          m_config.mean_molecular_weight_ionized * k_proton_mass_kg /
          k_boltzmann_j_per_k);
  const double n_h = m_config.n_h_threshold_cgs * entry.density_ratio;
  const CoolingRateResult rates = m_cooling_provider.lookupRates(CoolingRateQuery{
      .temperature_k = hot_temperature_k,
      .hydrogen_number_density_cgs = n_h,
      .metallicity_mass_fraction = 0.0,
      .redshift = 0.0,
  });
  const double net_cooling = rates.netCoolingRateErgCm3S();
  if (!std::isfinite(net_cooling) || net_cooling <= 0.0) {
    return entry;
  }
  const double thermal_energy_density_erg_cm3 =
      1.5 * n_h * k_boltzmann_erg_per_k * hot_temperature_k;
  const double cooling_time_s = thermal_energy_density_erg_cm3 / net_cooling;
  entry.cooling_timescale_code = cooling_time_s / m_units.timeSiPerCode();

  const double beta = m_config.massive_star_fraction;
  const double denominator = beta * m_config.supernova_specific_energy_code -
      (1.0 - beta) * u_cold;
  if (!std::isfinite(denominator) || denominator <= 0.0 ||
      !std::isfinite(entry.cooling_timescale_code) || entry.cooling_timescale_code <= 0.0) {
    return entry;
  }
  const double y = (entry.star_formation_timescale_code / entry.cooling_timescale_code) *
      entry.specific_internal_energy_hot_code / denominator;
  entry.cold_mass_fraction = stableColdFraction(y);
  if (!std::isfinite(entry.cold_mass_fraction)) {
    return entry;
  }
  const double u_full = (1.0 - entry.cold_mass_fraction) *
      entry.specific_internal_energy_hot_code + entry.cold_mass_fraction * u_cold;
  entry.specific_internal_energy_eff_code =
      m_config.q_eos * u_full + (1.0 - m_config.q_eos) * u_iso;
  entry.pressure_phys_code = (m_config.adiabatic_index - 1.0) *
      entry.density_phys_code * entry.specific_internal_energy_eff_code;
  entry.valid = std::isfinite(entry.specific_internal_energy_eff_code) &&
      entry.specific_internal_energy_eff_code > 0.0 &&
      std::isfinite(entry.pressure_phys_code) && entry.pressure_phys_code > 0.0;
  return entry;
}

void EffectiveMultiphaseEosTable::build() {
  m_entries.clear();
  m_entries.reserve(m_config.table_bins);
  const double log_max = std::log(m_config.max_density_ratio);
  for (std::uint32_t index = 0; index < m_config.table_bins; ++index) {
    const double fraction = static_cast<double>(index) /
        static_cast<double>(m_config.table_bins - 1U);
    const double ratio = std::exp(fraction * log_max);
    EffectiveMultiphaseEosEntry entry = evaluateDirect(m_threshold_density_phys_code * ratio);
    if (!entry.valid) {
      throw std::runtime_error("effective multiphase EOS table construction produced invalid equilibrium");
    }
    m_entries.push_back(entry);
  }
  normalizeThresholdContinuity();
  computeSignalSpeedDerivative();
  computeHash();
}

void EffectiveMultiphaseEosTable::normalizeThresholdContinuity() {
  const double u_iso = specificInternalEnergyFromTemperatureCode(
      m_config.isothermal_temperature_k, m_config.mean_molecular_weight_ionized);
  const double offset = u_iso - m_entries.front().specific_internal_energy_eff_code;
  const double u_cold = specificInternalEnergyFromTemperatureCode(
      m_config.cold_phase_temperature_k, m_config.mean_molecular_weight_cold);
  for (EffectiveMultiphaseEosEntry& entry : m_entries) {
    entry.specific_internal_energy_eff_code = std::max(
        entry.specific_internal_energy_eff_code + offset, u_cold);
    entry.pressure_phys_code = (m_config.adiabatic_index - 1.0) *
        entry.density_phys_code * entry.specific_internal_energy_eff_code;
  }
}

void EffectiveMultiphaseEosTable::computeSignalSpeedDerivative() {
  for (std::size_t index = 0; index < m_entries.size(); ++index) {
    const std::size_t left = index == 0 ? 0 : index - 1;
    const std::size_t right = index + 1 >= m_entries.size() ? m_entries.size() - 1 : index + 1;
    const double density_delta = m_entries[right].density_phys_code - m_entries[left].density_phys_code;
    double derivative = density_delta > 0.0
        ? (m_entries[right].pressure_phys_code - m_entries[left].pressure_phys_code) / density_delta
        : 0.0;
    if (!std::isfinite(derivative) || derivative <= 0.0) {
      derivative = m_config.adiabatic_index * m_entries[index].pressure_phys_code /
          std::max(m_entries[index].density_phys_code, k_small);
    }
    m_entries[index].signal_speed_squared_code = std::max(derivative, k_small);
  }
}

EffectiveMultiphaseEosLookup EffectiveMultiphaseEosTable::lookup(double density_phys_code) const {
  EffectiveMultiphaseEosLookup result;
  if (!std::isfinite(density_phys_code) || density_phys_code < m_threshold_density_phys_code) {
    return result;
  }
  result.above_threshold = true;
  const double ratio = std::clamp(
      density_phys_code / m_threshold_density_phys_code, 1.0, m_config.max_density_ratio);
  const double position = std::log(ratio) / std::log(m_config.max_density_ratio) *
      static_cast<double>(m_entries.size() - 1U);
  const std::size_t left = std::min<std::size_t>(
      static_cast<std::size_t>(std::floor(position)), m_entries.size() - 1U);
  const std::size_t right = std::min(left + 1U, m_entries.size() - 1U);
  const double fraction = std::clamp(position - static_cast<double>(left), 0.0, 1.0);
  const auto lerp = [&](double a, double b) { return a + fraction * (b - a); };
  const EffectiveMultiphaseEosEntry& a = m_entries[left];
  const EffectiveMultiphaseEosEntry& b = m_entries[right];
  result.entry = EffectiveMultiphaseEosEntry{
      .density_phys_code = density_phys_code,
      .density_ratio = ratio,
      .cold_mass_fraction = std::clamp(lerp(a.cold_mass_fraction, b.cold_mass_fraction), 0.0, 1.0),
      .specific_internal_energy_hot_code = lerp(a.specific_internal_energy_hot_code, b.specific_internal_energy_hot_code),
      .specific_internal_energy_eff_code = lerp(a.specific_internal_energy_eff_code, b.specific_internal_energy_eff_code),
      .pressure_phys_code = lerp(a.pressure_phys_code, b.pressure_phys_code),
      .signal_speed_squared_code = std::max(lerp(a.signal_speed_squared_code, b.signal_speed_squared_code), k_small),
      .star_formation_timescale_code = lerp(a.star_formation_timescale_code, b.star_formation_timescale_code),
      .cooling_timescale_code = lerp(a.cooling_timescale_code, b.cooling_timescale_code),
      .valid = a.valid && b.valid,
  };
  // Pressure is reconstructed at the requested density using interpolated u,
  // avoiding a density mismatch between table nodes.
  result.entry.pressure_phys_code = (m_config.adiabatic_index - 1.0) *
      density_phys_code * result.entry.specific_internal_energy_eff_code;
  result.valid = result.entry.valid && result.entry.specific_internal_energy_eff_code > 0.0 &&
      result.entry.signal_speed_squared_code > 0.0;
  return result;
}

std::span<const EffectiveMultiphaseEosEntry> EffectiveMultiphaseEosTable::entries() const noexcept {
  return m_entries;
}

void EffectiveMultiphaseEosTable::computeHash() {
  std::uint64_t hash = 1469598103934665603ULL;
  hash = fnvMix(hash, kEffectiveMultiphaseEosSchemaVersion);
  for (const char c : m_config.parameter_set) hash = fnvMix(hash, static_cast<unsigned char>(c));
  hash = hashDouble(hash, m_threshold_density_phys_code);
  hash = hashDouble(hash, m_config.q_eos);
  hash = hashDouble(hash, m_config.massive_star_fraction);
  hash = fnvMix(hash, m_entries.size());
  for (const EffectiveMultiphaseEosEntry& entry : m_entries) {
    hash = hashDouble(hash, entry.density_phys_code);
    hash = hashDouble(hash, entry.cold_mass_fraction);
    hash = hashDouble(hash, entry.specific_internal_energy_eff_code);
    hash = hashDouble(hash, entry.pressure_phys_code);
    hash = hashDouble(hash, entry.signal_speed_squared_code);
    hash = hashDouble(hash, entry.star_formation_timescale_code);
  }
  m_table_hash = hash == 0U ? 1U : hash;
}

std::uint64_t EffectiveMultiphaseEosTable::tableHash() const noexcept { return m_table_hash; }

std::string EffectiveMultiphaseEosTable::tableHashHex() const {
  std::ostringstream stream;
  stream << std::hex << std::setw(16) << std::setfill('0') << m_table_hash;
  return stream.str();
}

std::string EffectiveMultiphaseEosTable::coolingReferenceDescription() const {
  std::ostringstream stream;
  stream << "primordial_reference_z0_metallicity0_uv="
         << core::uvBackgroundModelToString(m_cooling_provider.config().uv_background_model)
         << "_shielding="
         << core::selfShieldingModelToString(m_cooling_provider.config().self_shielding_model);
  return stream.str();
}

EffectiveIsmThermodynamicClosure::EffectiveIsmThermodynamicClosure(
    EffectiveMultiphaseEosTable table)
    : m_table(std::make_shared<const EffectiveMultiphaseEosTable>(std::move(table))) {}

EffectiveIsmThermodynamicClosure::EffectiveIsmThermodynamicClosure(
    std::shared_ptr<const EffectiveMultiphaseEosTable> table)
    : m_table(std::move(table)) {
  if (m_table == nullptr) {
    throw std::invalid_argument("effective ISM thermodynamic closure requires an EOS table");
  }
}

hydro::HydroThermodynamicClosureResult EffectiveIsmThermodynamicClosure::evaluate(
    std::size_t cell_index,
    const hydro::HydroConservedState& conserved,
    const hydro::HydroPrimitiveState& ideal_primitive,
    double scale_factor,
    double redshift) const {
  (void)cell_index;
  (void)redshift;
  hydro::HydroThermodynamicClosureResult result{
      .pressure_comoving = ideal_primitive.pressure_comoving,
      .signal_speed_squared_code = ideal_primitive.signal_speed_squared_code,
      .target_specific_internal_energy_code = ideal_primitive.specific_internal_energy_code,
      .uses_effective_ism = false,
      .valid = true,
  };
  const double a = std::max(scale_factor, 1.0e-12);
  const double density_phys_code = conserved.mass_density_comoving / (a * a * a);
  const EffectiveMultiphaseEosLookup equilibrium = m_table->lookup(density_phys_code);
  if (!equilibrium.above_threshold || !equilibrium.valid) {
    return result;
  }
  result.target_specific_internal_energy_code = equilibrium.entry.specific_internal_energy_eff_code;
  const bool hot_above_eos = ideal_primitive.specific_internal_energy_code >
      equilibrium.entry.specific_internal_energy_eff_code *
          (1.0 + m_table->config().hot_excess_tolerance);
  if (hot_above_eos) {
    return result;
  }
  result.pressure_comoving = (m_table->config().adiabatic_index - 1.0) *
      conserved.mass_density_comoving * equilibrium.entry.specific_internal_energy_eff_code;
  result.signal_speed_squared_code = equilibrium.entry.signal_speed_squared_code;
  result.uses_effective_ism = true;
  result.valid = std::isfinite(result.pressure_comoving) && result.pressure_comoving > 0.0 &&
      std::isfinite(result.signal_speed_squared_code) && result.signal_speed_squared_code > 0.0;
  return result;
}

const EffectiveMultiphaseEosTable& EffectiveIsmThermodynamicClosure::table() const noexcept {
  return *m_table;
}

EffectiveIsmEnergyRelaxationSource::EffectiveIsmEnergyRelaxationSource(
    const EffectiveIsmThermodynamicClosure& closure)
    : m_closure(closure) {}

hydro::HydroConservedState EffectiveIsmEnergyRelaxationSource::sourceForCell(
    std::size_t cell_index,
    const hydro::HydroConservedState& conserved,
    const hydro::HydroPrimitiveState& primitive,
    const hydro::HydroSourceContext& context) const {
  hydro::HydroConservedState source;
  const hydro::HydroThermodynamicClosureResult closure = m_closure.evaluate(
      cell_index,
      conserved,
      primitive,
      context.update.scale_factor,
      context.redshift);
  const double target_u = closure.target_specific_internal_energy_code;
  const double actual_u = primitive.specific_internal_energy_code;
  if (!std::isfinite(target_u) || target_u <= 0.0 || !std::isfinite(actual_u) || actual_u <= 0.0 ||
      context.update.dt_code <= 0.0) {
    return source;
  }

  double adjustment_fraction = 1.0;
  if (m_closure.table().config().relaxation == core::EffectiveIsmEosRelaxation::kFiniteTimescale) {
    const double timescale = m_closure.table().config().relaxation_timescale_code;
    if (timescale <= 0.0) return source;
    adjustment_fraction = -std::expm1(-context.update.dt_code / timescale);
  }

  double delta_u = 0.0;
  if (actual_u < target_u) {
    delta_u = adjustment_fraction * (target_u - actual_u);
  } else if (m_closure.table().config().relaxation == core::EffectiveIsmEosRelaxation::kFiniteTimescale &&
      actual_u > target_u * (1.0 + m_closure.table().config().hot_excess_tolerance)) {
    delta_u = adjustment_fraction * (target_u - actual_u);
  }
  if (delta_u == 0.0) return source;

  const double delta_energy_density = conserved.mass_density_comoving * delta_u;
  source.total_energy_density_comoving = delta_energy_density / context.update.dt_code;
  // Audit F7: the ledger is an integrated code energy,
  //   Delta E_code = rho_comoving * Delta u * V_comoving  (= m_gas * Delta u),
  // so it is invariant under refinement (one parent vs eight children) and
  // correct for unequal-volume AMR cells. The cell volume is the authoritative
  // patch geometry volume injected by HydroCoreSolver.
  if (!(context.cell_volume_comoving > 0.0) || !std::isfinite(context.cell_volume_comoving)) {
    throw std::invalid_argument(
        "effective-ISM energy ledger requires HydroSourceContext::cell_volume_comoving > 0 "
        "(integrated energy = rho_comoving * delta_u * V_comoving)");
  }
  const double delta_energy_code = delta_energy_density * context.cell_volume_comoving;
  if (delta_energy_code > 0.0) {
    m_ledger.energy_added_code += delta_energy_code;
  } else {
    m_ledger.energy_removed_code += -delta_energy_code;
  }
  m_ledger.net_energy_adjustment_code += delta_energy_code;
  ++m_ledger.adjusted_cell_count;
  return source;
}

void EffectiveIsmEnergyRelaxationSource::resetLedger() const noexcept { m_ledger = {}; }
EffectiveIsmEnergyLedger EffectiveIsmEnergyRelaxationSource::ledger() const noexcept { return m_ledger; }

EffectiveMultiphaseEosConfig makeEffectiveMultiphaseEosConfig(
    const core::PhysicsConfig& physics_config) {
  return EffectiveMultiphaseEosConfig{
      .parameter_set = physics_config.sf_effective_parameter_set,
      .hydrogen_mass_fraction = 0.76,
      .n_h_threshold_cgs = physics_config.sf_effective_n_h_threshold_cgs,
      .min_baryon_overdensity = physics_config.sf_effective_min_baryon_overdensity,
      .t_star_at_threshold_code = physics_config.sf_effective_t_star_at_threshold_code,
      .evaporation_factor_at_threshold = physics_config.sf_effective_evaporation_factor_at_threshold,
      .evaporation_density_exponent = physics_config.sf_effective_evaporation_density_exponent,
      .cold_phase_temperature_k = physics_config.sf_effective_cold_phase_temperature_k,
      .supernova_specific_energy_code = physics_config.sf_effective_supernova_specific_energy_code,
      .massive_star_fraction = physics_config.sf_effective_massive_star_fraction,
      .q_eos = physics_config.sf_effective_q_eos,
      .isothermal_temperature_k = physics_config.sf_effective_isothermal_temperature_k,
      .hot_excess_tolerance = physics_config.sf_effective_hot_excess_tolerance,
      .relaxation = physics_config.sf_effective_eos_relaxation,
      .relaxation_timescale_code = physics_config.sf_effective_eos_relaxation_timescale_code,
      .table_bins = physics_config.sf_effective_eos_table_bins,
      .max_density_ratio = physics_config.sf_effective_eos_max_density_ratio,
      .birth_mass_convention = physics_config.sf_effective_birth_mass_convention,
      .feedback_coupling = physics_config.sf_effective_feedback_coupling,
  };
}

CoolingRateProvider makeEffectiveIsmReferenceCoolingProvider(
    const core::PhysicsConfig& physics_config) {
  CoolingModelConfig config;
  config.uv_background_model = physics_config.uv_background_model;
  config.self_shielding_model = physics_config.self_shielding_model;
  config.temperature_floor_k = physics_config.temperature_floor_k;
  config.enable_metal_line_cooling = false;
  return CoolingRateProvider(config);
}

}  // namespace cosmosim::physics
