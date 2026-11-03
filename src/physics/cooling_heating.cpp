#include "cosmosim/physics/cooling_heating.hpp"

#include "cosmosim/physics/source_interval.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cosmosim::physics {
namespace {

constexpr double k_small = 1.0e-30;

[[nodiscard]] double clampTemperature(double temperature_k) {
  return std::max(temperature_k, 1.0);
}

[[nodiscard]] double safeLog10(double value) {
  return std::log10(std::max(value, k_small));
}

[[nodiscard]] double linearInterpolate(double x0, double x1, double y0, double y1, double x) {
  if (std::abs(x1 - x0) < 1.0e-16) {
    return y0;
  }
  const double t = (x - x0) / (x1 - x0);
  return (1.0 - t) * y0 + t * y1;
}

[[nodiscard]] double uvHeatingScale(core::UvBackgroundModel uv_background_model, double redshift) {
  switch (uv_background_model) {
    case core::UvBackgroundModel::kNone:
      return 0.0;
    case core::UvBackgroundModel::kHm12: {
      const double z_peak = 2.0;
      const double width = 2.5;
      const double dz = (redshift - z_peak) / width;
      return std::exp(-dz * dz);
    }
    case core::UvBackgroundModel::kFg20:
      return 1.0 / (1.0 + std::pow((1.0 + redshift) / 4.0, 2.0));
  }
  throw std::invalid_argument("CoolingRateProvider: unhandled UvBackgroundModel enum value");
}

[[nodiscard]] double shieldingFactor(
    core::SelfShieldingModel self_shielding_model,
    double hydrogen_number_density_cgs) {
  switch (self_shielding_model) {
    case core::SelfShieldingModel::kNone:
      return 1.0;
    case core::SelfShieldingModel::kRahmati13Like: {
      const double threshold = 1.0e-2;
      return 1.0 / (1.0 + std::pow(hydrogen_number_density_cgs / threshold, 1.6));
    }
  }
  throw std::invalid_argument("CoolingRateProvider: unhandled SelfShieldingModel enum value");
}

[[nodiscard]] double effectiveHydrogenNumberDensity(
    double mass_density_comoving,
    double default_hydrogen_number_density_cgs) {
  if (default_hydrogen_number_density_cgs > 0.0) {
    return default_hydrogen_number_density_cgs;
  }
  return std::max(mass_density_comoving, 1.0e-10);
}

}  // namespace

double CoolingRateResult::netCoolingRateErgCm3S() const {
  return cooling_rate_erg_cm3_s - heating_rate_erg_cm3_s;
}

MetalLineCoolingTable MetalLineCoolingTable::loadFromTextFile(
    const std::string& path,
    const std::string& table_label) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("MetalLineCoolingTable: failed to open file '" + path + "'");
  }

  MetalLineCoolingTable table;
  table.m_provenance.source_path = path;
  table.m_provenance.source_format = "text_two_column";
  table.m_provenance.table_label = table_label;

  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') {
      continue;
    }
    std::istringstream stream(line);
    double log10_temperature_k = 0.0;
    double log10_cooling_rate_erg_cm3_s = 0.0;
    stream >> log10_temperature_k >> log10_cooling_rate_erg_cm3_s;
    if (stream.fail()) {
      throw std::runtime_error("MetalLineCoolingTable: malformed row in '" + path + "'");
    }
    table.m_log10_temperature_k.push_back(log10_temperature_k);
    table.m_log10_cooling_rate_erg_cm3_s.push_back(log10_cooling_rate_erg_cm3_s);
  }

  if (table.m_log10_temperature_k.empty()) {
    throw std::runtime_error("MetalLineCoolingTable: no data rows in '" + path + "'");
  }

  for (std::size_t i = 1; i < table.m_log10_temperature_k.size(); ++i) {
    if (table.m_log10_temperature_k[i] <= table.m_log10_temperature_k[i - 1]) {
      throw std::runtime_error("MetalLineCoolingTable: log10(T) must be strictly increasing");
    }
  }

  table.m_provenance.row_count = table.m_log10_temperature_k.size();
  return table;
}

bool MetalLineCoolingTable::empty() const noexcept {
  return m_log10_temperature_k.empty();
}

double MetalLineCoolingTable::lookupCoolingRateErgCm3S(double temperature_k) const {
  if (m_log10_temperature_k.empty()) {
    return 0.0;
  }

  const double target = safeLog10(clampTemperature(temperature_k));
  const auto upper = std::upper_bound(m_log10_temperature_k.begin(), m_log10_temperature_k.end(), target);
  if (upper == m_log10_temperature_k.begin()) {
    return std::pow(10.0, m_log10_cooling_rate_erg_cm3_s.front());
  }
  if (upper == m_log10_temperature_k.end()) {
    return std::pow(10.0, m_log10_cooling_rate_erg_cm3_s.back());
  }

  const std::size_t right_index = static_cast<std::size_t>(std::distance(m_log10_temperature_k.begin(), upper));
  const std::size_t left_index = right_index - 1;
  const double interpolated_log10_rate = linearInterpolate(
      m_log10_temperature_k[left_index],
      m_log10_temperature_k[right_index],
      m_log10_cooling_rate_erg_cm3_s[left_index],
      m_log10_cooling_rate_erg_cm3_s[right_index],
      target);
  return std::pow(10.0, interpolated_log10_rate);
}

const CoolingTableProvenance& MetalLineCoolingTable::provenance() const noexcept {
  return m_provenance;
}

CoolingRateProvider::CoolingRateProvider(CoolingModelConfig config) : m_config(std::move(config)) {
  if (m_config.temperature_floor_k <= 0.0) {
    throw std::invalid_argument("CoolingRateProvider: temperature_floor_k must be > 0");
  }
  if (m_config.max_fractional_energy_change_per_substep <= 0.0 ||
      m_config.max_fractional_energy_change_per_substep >= 1.0) {
    throw std::invalid_argument(
        "CoolingRateProvider: max_fractional_energy_change_per_substep must be in (0, 1)");
  }
  if (m_config.max_subcycles <= 0) {
    throw std::invalid_argument("CoolingRateProvider: max_subcycles must be > 0");
  }
  if (m_config.temperature_per_internal_energy_k <= 0.0) {
    throw std::invalid_argument("CoolingRateProvider: temperature_per_internal_energy_k must be > 0");
  }
}

void CoolingRateProvider::setMetalLineTable(MetalLineCoolingTable table) {
  m_metal_line_table = std::move(table);
}

bool CoolingRateProvider::hasMetalLineTable() const noexcept {
  return m_metal_line_table.has_value();
}

CoolingRateResult CoolingRateProvider::lookupRates(const CoolingRateQuery& query) const {
  const double temperature_k = std::max(query.temperature_k, m_config.temperature_floor_k);
  const double n_h = std::max(query.hydrogen_number_density_cgs, 1.0e-12);

  // Primordial baseline: compact analytic fit with explicit UV/shielding controls.
  const double primordial_cooling_coeff = 1.0e-22 * std::sqrt(temperature_k / 1.0e4);
  const double primordial_cooling = primordial_cooling_coeff * n_h * n_h;

  const double uv_scale = uvHeatingScale(m_config.uv_background_model, query.redshift);
  const double shielding = shieldingFactor(m_config.self_shielding_model, n_h);
  const double primordial_heating = 2.0e-24 * uv_scale * shielding * n_h;

  double metal_line_cooling = 0.0;
  if (m_config.enable_metal_line_cooling && m_metal_line_table.has_value()) {
    metal_line_cooling = m_metal_line_table->lookupCoolingRateErgCm3S(temperature_k) *
        std::max(query.metallicity_mass_fraction, 0.0);
  }

  CoolingRateResult result;
  result.cooling_rate_erg_cm3_s = primordial_cooling + metal_line_cooling;
  result.heating_rate_erg_cm3_s = primordial_heating;
  return result;
}

const CoolingModelConfig& CoolingRateProvider::config() const noexcept {
  return m_config;
}

std::optional<CoolingTableProvenance> CoolingRateProvider::metalLineTableProvenance() const {
  if (!m_metal_line_table.has_value()) {
    return std::nullopt;
  }
  return m_metal_line_table->provenance();
}

CoolingSourceIntegrator::CoolingSourceIntegrator(
    double energy_floor_code,
    double specific_energy_rate_cgs_to_code)
    : m_energy_floor_code(energy_floor_code),
      m_specific_energy_rate_cgs_to_code(specific_energy_rate_cgs_to_code) {
  if (m_energy_floor_code <= 0.0) {
    throw std::invalid_argument("CoolingSourceIntegrator: energy_floor_code must be > 0");
  }
  if (!std::isfinite(m_specific_energy_rate_cgs_to_code) ||
      m_specific_energy_rate_cgs_to_code <= 0.0) {
    throw std::invalid_argument(
        "CoolingSourceIntegrator: specific-energy rate conversion must be finite and > 0");
  }
}

CoolingIntegrationResult CoolingSourceIntegrator::integrateSpecificInternalEnergy(
    double specific_internal_energy_code,
    double mass_density_physical_cgs,
    double dt_code,
    const CoolingRateQuery& query,
    const CoolingRateProvider& rate_provider) const {
  if (dt_code <= 0.0 || !std::isfinite(dt_code)) {
    throw std::invalid_argument("CoolingSourceIntegrator: dt_code must be finite and > 0");
  }
  // CHUI Part 2 audit F2 repair. The former explicit loop exited after
  // max_subcycles with time still remaining and the caller accepted the
  // partial result. The interval is now integrated with an implicit
  // (backward-Euler) scalar step, Newton-solved, with step-doubling local
  // error control and Richardson correction; max_subcycles bounds the number
  // of step ATTEMPTS. The result is complete only when the long-double clock
  // closes the requested interval exactly; otherwise `completed` is false,
  // the energy returned is the UNCHANGED input energy and the caller must
  // reject the stage (CoolingHeatingSource throws SourceIntervalIncompleteError).
  const CoolingModelConfig& config = rate_provider.config();
  const double rho = std::max(mass_density_physical_cgs, 1.0e-40);
  const auto rate = [&](double energy) {
    CoolingRateQuery local_query = query;
    local_query.temperature_k = std::max(energy * config.temperature_per_internal_energy_k,
                                         config.temperature_floor_k);
    const CoolingRateResult rates = rate_provider.lookupRates(local_query);
    return (rates.heating_rate_erg_cm3_s - rates.cooling_rate_erg_cm3_s) / rho *
        m_specific_energy_rate_cgs_to_code;
  };
  const double energy_begin = std::max(specific_internal_energy_code, m_energy_floor_code);
  // Implicit step e = e0 + h f(e): damped Newton, numerical derivative.
  const auto implicit_step = [&](double e0, double h, double& e) -> bool {
    e = e0;
    for (int iteration = 0; iteration < 30; ++iteration) {
      const double f = rate(e);
      const double delta = 1.0e-7 * std::max(e, m_energy_floor_code);
      const double dfde = (rate(e + delta) - f) / delta;
      const double g = e - e0 - h * f;
      const double dg = 1.0 - h * dfde;
      if (!std::isfinite(f) || !std::isfinite(dg)) return false;
      double step = (dg != 0.0) ? -g / dg : -g;
      // Keep the iterate above the physical floor (damping, not clipping of
      // the converged answer).
      if (e + step < m_energy_floor_code) step = 0.5 * (m_energy_floor_code - e) + 0.5 * step * 0.5;
      e = std::max(e + step, m_energy_floor_code);
      if (std::abs(step) <= 1.0e-12 * std::max(e, m_energy_floor_code)) return true;
    }
    return false;
  };

  CoolingIntegrationResult result;
  result.requested_dt_code = dt_code;
  const double rtol = std::min(config.max_fractional_energy_change_per_substep, 0.1) * 1.0e-2;
  const double f0 = rate(energy_begin);
  double h = std::abs(f0) > 0.0
      ? std::min(dt_code, config.max_fractional_energy_change_per_substep * energy_begin / std::abs(f0))
      : dt_code;
  const double h_min = 1.0e-14 * dt_code;
  long double t_done = 0.0L;
  const long double t_request = static_cast<long double>(dt_code);
  double energy = energy_begin;
  int attempts = 0;
  while (!physics::sourceIntervalClosed(t_done, t_request)) {
    if (attempts >= config.max_subcycles) {
      result.completed = false;
      result.diagnostics.hit_subcycle_limit = true;
      result.diagnostics.subcycles_used = attempts;
      result.completed_dt_code = static_cast<double>(t_done);
      result.specific_internal_energy_code = specific_internal_energy_code;  // unchanged
      result.diagnostics.suggested_next_dt_code = std::max(static_cast<double>(t_done), h_min);
      return result;
    }
    ++attempts;
    const double remaining = static_cast<double>(t_request - t_done);
    bool final_step = false;
    if (h >= remaining * (1.0 - 1.0e-12)) { h = remaining; final_step = true; }
    double e_full = 0.0, e_half = 0.0, e_half2 = 0.0;
    if (!implicit_step(energy, h, e_full) || !implicit_step(energy, 0.5 * h, e_half) ||
        !implicit_step(e_half, 0.5 * h, e_half2)) {
      ++result.diagnostics.rejected_attempts;
      h *= 0.25;
      if (h < h_min) { attempts = config.max_subcycles; }
      continue;
    }
    const double err = std::abs(e_full - e_half2) /
        (rtol * std::max(e_half2, energy) + 1.0e-6 * m_energy_floor_code);
    if (err > 1.0) {
      ++result.diagnostics.rejected_attempts;
      h *= std::max(0.2, 0.9 / std::sqrt(err));
      if (h < h_min) { attempts = config.max_subcycles; }
      continue;
    }
    result.diagnostics.max_error_estimate = std::max(result.diagnostics.max_error_estimate, err);
    const double richardson = 2.0 * e_half2 - e_full;
    double e_new = richardson > m_energy_floor_code ? richardson : e_half2;
    if (e_new <= m_energy_floor_code) {
      // Physical floor reached by a converged step: counted, not hidden.
      e_new = m_energy_floor_code;
      ++result.diagnostics.floor_activations;
    }
    energy = e_new;
    ++result.diagnostics.accepted_substeps;
    t_done = final_step ? t_request : t_done + static_cast<long double>(h);
    h *= std::min(4.0, 0.9 / std::sqrt(std::max(err, 1.0e-10)));
  }
  result.completed = true;
  result.completed_dt_code = dt_code;
  result.diagnostics.subcycles_used = attempts;
  result.diagnostics.hit_subcycle_limit = false;
  const double final_rate = rate(energy);
  result.diagnostics.suggested_next_dt_code = config.max_fractional_energy_change_per_substep *
      std::max(energy, m_energy_floor_code) / std::max(std::abs(final_rate), 1.0e-30);
  result.specific_internal_energy_code = energy;
  return result;
}

CoolingHeatingSource::CoolingHeatingSource(CoolingRateProvider provider, CoolingSourceIntegrator integrator)
    : m_rate_provider(std::move(provider)), m_integrator(std::move(integrator)) {}

hydro::HydroConservedState CoolingHeatingSource::sourceForCell(
    std::size_t cell_index,
    const hydro::HydroConservedState& conserved,
    const hydro::HydroPrimitiveState& primitive,
    const hydro::HydroSourceContext& context) const {
  const double kinetic_density = 0.5 *
      (conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
       conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
       conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving) /
      std::max(conserved.mass_density_comoving, 1.0e-16);
  const double internal_density = std::max(conserved.total_energy_density_comoving - kinetic_density, 0.0);
  const double specific_internal_energy = internal_density / std::max(conserved.mass_density_comoving, 1.0e-16);

  hydro::HydroCellSourceProperties properties;
  if (context.source_property_provider != nullptr) {
    properties = context.source_property_provider->sourcePropertiesForCell(cell_index);
  } else {
    if (cell_index >= context.hydrogen_number_density_cgs.size() ||
        cell_index >= context.temperature_k.size()) {
      throw std::runtime_error(
          "CoolingHeatingSource requires explicit physical hydrogen density and temperature for every cooled cell");
    }
    properties.hydrogen_number_density_cgs = context.hydrogen_number_density_cgs[cell_index];
    properties.temperature_k = context.temperature_k[cell_index];
    properties.metallicity_mass_fraction =
        (cell_index < context.metallicity_mass_fraction.size())
            ? context.metallicity_mass_fraction[cell_index]
            : 0.0;
    properties.mass_density_physical_cgs =
        (cell_index < context.mass_density_physical_cgs.size())
            ? context.mass_density_physical_cgs[cell_index]
            : 0.0;
  }
  const double hydrogen_density = properties.hydrogen_number_density_cgs;
  const double temperature_k = properties.temperature_k;
  if (!(hydrogen_density > 0.0) || !std::isfinite(hydrogen_density) ||
      !(temperature_k > 0.0) || !std::isfinite(temperature_k)) {
    throw std::runtime_error(
        "CoolingHeatingSource requires finite positive physical hydrogen density and temperature");
  }
  const double metallicity = properties.metallicity_mass_fraction;

  const CoolingRateQuery query{
      .temperature_k = temperature_k,
      .hydrogen_number_density_cgs = hydrogen_density,
      .metallicity_mass_fraction = metallicity,
      .redshift = context.redshift};

  if (primitive.uses_effective_ism) {
    // The effective multiphase closure owns its unresolved thermal equilibrium;
    // applying the ordinary cooling source as well would double-count it.
    return {};
  }
  const double physical_mass_density_cgs = properties.mass_density_physical_cgs;
  if (!(physical_mass_density_cgs > 0.0) || !std::isfinite(physical_mass_density_cgs)) {
    throw std::runtime_error(
        "CoolingHeatingSource requires finite positive physical CGS mass density");
  }

  const CoolingIntegrationResult integrated = m_integrator.integrateSpecificInternalEnergy(
      specific_internal_energy,
      physical_mass_density_cgs,
      context.update.dt_code,
      query,
      m_rate_provider);
  if (!integrated.completed) {
    // Audit F2: never accept a partial interval. The hydro stage converts
    // this into a StageTransaction rejection (kNumericalRetry) and nothing
    // is committed.
    throw SourceIntervalIncompleteError(
        "cooling source could not integrate the accepted interval at solver row " +
            std::to_string(cell_index) + " (completed " + std::to_string(integrated.completed_dt_code) +
            " of " + std::to_string(context.update.dt_code) + " code time after " +
            std::to_string(integrated.diagnostics.subcycles_used) + " step attempts)",
        integrated.diagnostics.suggested_next_dt_code);
  }

  hydro::HydroConservedState source{};
  source.total_energy_density_comoving =
      (integrated.specific_internal_energy_code - specific_internal_energy) *
      conserved.mass_density_comoving / std::max(context.update.dt_code, 1.0e-30);
  return source;
}

}  // namespace cosmosim::physics
