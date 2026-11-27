#include "cosmosim/core/cosmology.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::core {

namespace {

[[nodiscard]] double adaptiveSimpsonRecursive(
    const std::function<double(double)>& integrand,
    double left,
    double right,
    double f_left,
    double f_mid,
    double f_right,
    double whole,
    double tolerance,
    int depth_remaining) {
  const double mid = 0.5 * (left + right);
  const double left_mid = 0.5 * (left + mid);
  const double right_mid = 0.5 * (mid + right);
  const double f_left_mid = integrand(left_mid);
  const double f_right_mid = integrand(right_mid);
  const double left_part = (mid - left) * (f_left + 4.0 * f_left_mid + f_mid) / 6.0;
  const double right_part = (right - mid) * (f_mid + 4.0 * f_right_mid + f_right) / 6.0;
  const double refined = left_part + right_part;
  const double error = refined - whole;
  if (!std::isfinite(refined) || !std::isfinite(error)) {
    throw std::runtime_error("adaptive cosmology integration produced a non-finite estimate");
  }
  if (std::abs(error) <= 15.0 * tolerance) {
    return refined + error / 15.0;
  }
  if (depth_remaining <= 0) {
    throw std::runtime_error("adaptive cosmology integration failed to converge within recursion limit");
  }
  return adaptiveSimpsonRecursive(
             integrand, left, mid, f_left, f_left_mid, f_mid, left_part,
             0.5 * tolerance, depth_remaining - 1) +
         adaptiveSimpsonRecursive(
             integrand, mid, right, f_mid, f_right_mid, f_right, right_part,
             0.5 * tolerance, depth_remaining - 1);
}

[[nodiscard]] double adaptiveSimpson(
    const std::function<double(double)>& integrand,
    double left,
    double right,
    double tolerance) {
  const double mid = 0.5 * (left + right);
  const double f_left = integrand(left);
  const double f_mid = integrand(mid);
  const double f_right = integrand(right);
  const double whole = (right - left) * (f_left + 4.0 * f_mid + f_right) / 6.0;
  return adaptiveSimpsonRecursive(
      integrand, left, right, f_left, f_mid, f_right, whole, tolerance, 24);
}

// --------------------------------------------------------------------------
// Fermi-Dirac moment integrals (S01/N01). Gauss-Legendre nodes are computed
// by Newton iteration on P_n (no transcribed node tables).
constexpr int k_fd_gauss_order = 16;
constexpr int k_fd_panel_count = 40;
constexpr double k_fd_x_max = 80.0;  // e^-80 x^6 < 1e-22 relative contribution
constexpr double k_fd_table_y_min = 1.0e-5;
constexpr double k_fd_table_y_max = 1.0e6;
constexpr std::size_t k_fd_table_points = 4001U;

struct GaussLegendreRule {
  std::array<double, k_fd_gauss_order> nodes{};
  std::array<double, k_fd_gauss_order> weights{};
};

[[nodiscard]] GaussLegendreRule buildGaussLegendreRule() {
  GaussLegendreRule rule;
  constexpr int n = k_fd_gauss_order;
  for (int i = 0; i < n; ++i) {
    // Chebyshev-like initial guess for the i-th root of P_n on (-1, 1).
    double x = std::cos(constants::k_pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
    double derivative = 0.0;
    for (int iteration = 0; iteration < 100; ++iteration) {
      double p0 = 1.0;
      double p1 = x;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / static_cast<double>(k);
        p0 = p1;
        p1 = p2;
      }
      derivative = static_cast<double>(n) * (x * p1 - p0) / (x * x - 1.0);
      const double step = p1 / derivative;
      x -= step;
      if (std::abs(step) < 1.0e-16) break;
    }
    rule.nodes[static_cast<std::size_t>(i)] = x;
    rule.weights[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - x * x) * derivative * derivative);
  }
  return rule;
}

[[nodiscard]] const GaussLegendreRule& gaussLegendreRule() {
  static const GaussLegendreRule rule = buildGaussLegendreRule();
  return rule;
}

[[nodiscard]] double fermiDiracOccupation(double x) {
  // 1/(e^x + 1) evaluated without overflow for large x.
  return x > 40.0 ? std::exp(-x) / (1.0 + std::exp(-x)) : 1.0 / (std::exp(x) + 1.0);
}

template <typename Integrand>
[[nodiscard]] double integrateFermiDirac(const Integrand& integrand) {
  const GaussLegendreRule& rule = gaussLegendreRule();
  const double width = k_fd_x_max / static_cast<double>(k_fd_panel_count);
  double sum = 0.0;
  for (int panel = 0; panel < k_fd_panel_count; ++panel) {
    const double left = width * static_cast<double>(panel);
    const double half = 0.5 * width;
    const double mid = left + half;
    double panel_sum = 0.0;
    for (std::size_t k = 0; k < rule.nodes.size(); ++k) {
      const double x = mid + half * rule.nodes[k];
      panel_sum += rule.weights[k] * integrand(x) * fermiDiracOccupation(x);
    }
    sum += half * panel_sum;
  }
  return sum;
}

void requireAdmissibleMassRatio(double y) {
  if (!std::isfinite(y) || y < 0.0) {
    throw std::invalid_argument("Fermi-Dirac moment requires finite y = m c^2 / (k T) >= 0");
  }
}

[[nodiscard]] double energyQuadrature(double y) {
  return integrateFermiDirac([y](double x) { return x * x * std::sqrt(x * x + y * y); });
}

[[nodiscard]] double pressureQuadrature(double y) {
  return integrateFermiDirac([y](double x) {
    const double energy = std::sqrt(x * x + y * y);
    return energy > 0.0 ? (x * x * x * x) / (3.0 * energy) : 0.0;
  });
}

struct FermiDiracTable {
  double ln_y_min = 0.0;
  double step = 0.0;
  std::vector<double> ln_energy;
  std::vector<double> ln_pressure;
};

[[nodiscard]] FermiDiracTable buildFermiDiracTable() {
  FermiDiracTable table;
  table.ln_y_min = std::log(k_fd_table_y_min);
  table.step = (std::log(k_fd_table_y_max) - table.ln_y_min) / static_cast<double>(k_fd_table_points - 1U);
  table.ln_energy.resize(k_fd_table_points);
  table.ln_pressure.resize(k_fd_table_points);
  for (std::size_t i = 0; i < k_fd_table_points; ++i) {
    const double y = std::exp(table.ln_y_min + table.step * static_cast<double>(i));
    table.ln_energy[i] = std::log(energyQuadrature(y));
    table.ln_pressure[i] = std::log(pressureQuadrature(y));
  }
  return table;
}

[[nodiscard]] const FermiDiracTable& fermiDiracTable() {
  static const FermiDiracTable table = buildFermiDiracTable();
  return table;
}

[[nodiscard]] double interpolateLnTable(const FermiDiracTable& table, const std::vector<double>& values, double y) {
  const double u = (std::log(y) - table.ln_y_min) / table.step;
  const auto last = static_cast<std::ptrdiff_t>(values.size()) - 1;
  auto j = static_cast<std::ptrdiff_t>(std::floor(u));
  j = std::clamp<std::ptrdiff_t>(j, 1, last - 2);
  const double t = u - static_cast<double>(j);  // position relative to node j
  // 4-point Lagrange on nodes j-1, j, j+1, j+2 (offsets -1, 0, 1, 2).
  const double lm1 = -t * (t - 1.0) * (t - 2.0) / 6.0;
  const double l0 = (t + 1.0) * (t - 1.0) * (t - 2.0) / 2.0;
  const double l1 = -(t + 1.0) * t * (t - 2.0) / 2.0;
  const double l2 = (t + 1.0) * t * (t - 1.0) / 6.0;
  const auto at = [&values](std::ptrdiff_t index) { return values[static_cast<std::size_t>(index)]; };
  return std::exp(lm1 * at(j - 1) + l0 * at(j) + l1 * at(j + 1) + l2 * at(j + 2));
}

// Massless and nonrelativistic limits (expansion coefficients from the
// closed-form moments int x^n/(e^x+1) = (1 - 2^-n) n! zeta(n+1)).
[[nodiscard]] double energyMasslessMoment() { return 7.0 * std::pow(constants::k_pi, 4) / 120.0; }
[[nodiscard]] double firstMoment() { return constants::k_pi * constants::k_pi / 12.0; }       // int x f
[[nodiscard]] double secondMoment() { return 1.5 * constants::k_zeta_3; }                     // int x^2 f
[[nodiscard]] double fourthMoment() { return 22.5 * constants::k_zeta_5; }                    // int x^4 f

}  // namespace

double instantaneousDecouplingNeutrinoTemperatureRatio() { return std::cbrt(4.0 / 11.0); }

double fermiDiracEnergyIntegralQuadrature(double y) {
  requireAdmissibleMassRatio(y);
  return energyQuadrature(y);
}

double fermiDiracPressureIntegralQuadrature(double y) {
  requireAdmissibleMassRatio(y);
  return pressureQuadrature(y);
}

double fermiDiracNumberIntegral() { return secondMoment(); }

double fermiDiracMeanMomentumOverTemperature() { return energyMasslessMoment() / secondMoment(); }

double fermiDiracEnergyIntegral(double y) {
  requireAdmissibleMassRatio(y);
  if (y < k_fd_table_y_min) {
    // sqrt(x^2+y^2) = x + y^2/(2x) + O(y^4/x^3): remainder O(y^4 ln y) < 1e-19.
    return energyMasslessMoment() + 0.5 * y * y * firstMoment();
  }
  if (y > k_fd_table_y_max) {
    // sqrt(x^2+y^2) = y + x^2/(2y) + O(x^4/y^3): relative remainder < 1e-16.
    return y * secondMoment() + fourthMoment() / (2.0 * y);
  }
  return interpolateLnTable(fermiDiracTable(), fermiDiracTable().ln_energy, y);
}

double fermiDiracPressureIntegral(double y) {
  requireAdmissibleMassRatio(y);
  if (y < k_fd_table_y_min) {
    return energyMasslessMoment() / 3.0 - y * y * firstMoment() / 6.0;
  }
  if (y > k_fd_table_y_max) {
    // x^4/(3 sqrt(x^2+y^2)) = x^4/(3y) - x^6/(6 y^3) + ...; second term relative < 1e-11.
    const double sixth_moment = (63.0 / 64.0) * 720.0 * 1.008349277381922827;  // (1-2^-6) 6! zeta(7)
    return fourthMoment() / (3.0 * y) - sixth_moment / (6.0 * y * y * y);
  }
  return interpolateLnTable(fermiDiracTable(), fermiDiracTable().ln_pressure, y);
}



LambdaCdmBackground::LambdaCdmBackground(CosmologyBackgroundConfig config) : m_config(config) {
  if (!std::isfinite(m_config.hubble_param) || !(m_config.hubble_param > 0.0)) {
    throw std::invalid_argument("hubble_param must be finite and positive");
  }
  if (!std::isfinite(m_config.omega_matter) || !std::isfinite(m_config.omega_lambda) ||
      !std::isfinite(m_config.omega_radiation) || !std::isfinite(m_config.omega_curvature)) {
    throw std::invalid_argument("cosmology density parameters must be finite");
  }
  if (m_config.omega_matter < 0.0 || m_config.omega_lambda < 0.0 || m_config.omega_radiation < 0.0) {
    throw std::invalid_argument("matter, lambda, and radiation density parameters must be non-negative");
  }
  constexpr double k_density_closure_tolerance = 1.0e-10;
  const RadiationSpeciesConfig& species = m_config.radiation_species;
  if (!species.enabled) {
    if (!(species == RadiationSpeciesConfig{})) {
      throw std::invalid_argument(
          "radiation_species parameters are set but radiation_species.enabled == false (inert values)");
    }
    const double density_sum = m_config.omega_radiation + m_config.omega_matter +
        m_config.omega_curvature + m_config.omega_lambda;
    if (!std::isfinite(density_sum) ||
        std::abs(density_sum - 1.0) > k_density_closure_tolerance) {
      throw std::invalid_argument(
          "cosmology density parameters must sum to 1 within 1e-10 so H(a=1) == H0");
    }
    m_omega_cb0 = m_config.omega_matter;
    return;
  }
  // S01: relativistic species derived from T_cmb0, n_ur and eigenstate masses.
  if (m_config.omega_radiation != 0.0) {
    throw std::invalid_argument(
        "omega_radiation must be 0 when radiation species are derived from T_cmb0 (no double counting)");
  }
  if (!std::isfinite(species.t_cmb0_k) || !(species.t_cmb0_k > 0.0)) {
    throw std::invalid_argument("radiation species require a finite T_cmb0 > 0 K");
  }
  if (!std::isfinite(species.n_ur) || species.n_ur < 0.0) {
    throw std::invalid_argument("n_ur (massless effective species) must be finite and >= 0");
  }
  if (!species.ncdm_masses_ev.empty() &&
      (!std::isfinite(species.ncdm_temperature_ratio) || !(species.ncdm_temperature_ratio > 0.0))) {
    throw std::invalid_argument("massive neutrino eigenstates require a finite T_ncdm/T_cmb ratio > 0");
  }
  if (species.ncdm_masses_ev.empty() && species.ncdm_temperature_ratio != 0.0) {
    throw std::invalid_argument("ncdm_temperature_ratio is set but no massive eigenstate is configured");
  }
  for (const double mass_ev : species.ncdm_masses_ev) {
    if (!std::isfinite(mass_ev) || !(mass_ev > 0.0)) {
      throw std::invalid_argument(
          "massive neutrino masses must be finite and > 0 eV (massless species belong in n_ur)");
    }
  }
  const double rho_crit0 = criticalDensity0Si();
  const double c2 = constants::k_speed_of_light_si * constants::k_speed_of_light_si;
  const double t4 = std::pow(species.t_cmb0_k, 4);
  m_omega_gamma0 = constants::k_radiation_constant_si * t4 / (c2 * rho_crit0);
  const double ur_ratio4 = std::pow(instantaneousDecouplingNeutrinoTemperatureRatio(), 4);
  m_omega_ur0 = species.n_ur * (7.0 / 8.0) * ur_ratio4 * m_omega_gamma0;
  if (!species.ncdm_masses_ev.empty()) {
    m_omega_ncdm_unit0 = (7.0 / 8.0) * std::pow(species.ncdm_temperature_ratio, 4) * m_omega_gamma0;
    const double kt_ncdm0_j = constants::k_boltzmann_si * species.ncdm_temperature_ratio * species.t_cmb0_k;
    m_ncdm_y0.reserve(species.ncdm_masses_ev.size());
    for (const double mass_ev : species.ncdm_masses_ev) {
      m_ncdm_y0.push_back(mass_ev * constants::k_electron_volt_si / kt_ncdm0_j);
    }
    for (std::size_t i = 0; i < m_ncdm_y0.size(); ++i) {
      m_omega_ncdm0 += massiveNeutrinoDensityParameter(i, 1.0);
    }
  }
  m_omega_cb0 = m_config.omega_matter - m_omega_ncdm0;
  if (!(m_omega_cb0 > 0.0)) {
    throw std::invalid_argument(
        "omega_matter must exceed the present massive-neutrino density Omega_ncdm,0 = " +
        std::to_string(m_omega_ncdm0) + " (omega_matter includes massive neutrinos; Omega_cb must be > 0)");
  }
  const double density_sum = m_omega_gamma0 + m_omega_ur0 + m_config.omega_matter +
      m_config.omega_curvature + m_config.omega_lambda;
  if (!std::isfinite(density_sum) || std::abs(density_sum - 1.0) > k_density_closure_tolerance) {
    const double required_lambda = 1.0 - (m_omega_gamma0 + m_omega_ur0 + m_config.omega_matter +
                                          m_config.omega_curvature);
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.17g", required_lambda);
    throw std::invalid_argument(
        std::string("with photons/neutrinos derived from T_cmb0 the density parameters must close to 1 within "
                    "1e-10; strict closure requires omega_lambda = ") + buffer +
        " (Omega_gamma + Omega_ur + Omega_matter + Omega_k + Omega_lambda = 1; nothing is renormalized)");
  }
}

const CosmologyBackgroundConfig& LambdaCdmBackground::config() const { return m_config; }

double LambdaCdmBackground::hubble0Si() const {
  return m_config.hubble_param * constants::k_hubble_100_km_s_mpc_si;
}

double LambdaCdmBackground::eFactor(double scale_factor) const {
  if (!std::isfinite(scale_factor) || !(scale_factor > 0.0)) {
    throw std::invalid_argument("scale_factor must be finite and positive");
  }
  const double a2 = scale_factor * scale_factor;
  const double a3 = a2 * scale_factor;
  const double a4 = a2 * a2;

  double density_sum = 0.0;
  if (!m_config.radiation_species.enabled) {
    // Legacy: E(a)^2 = omega_r a^-4 + omega_m a^-3 + omega_k a^-2 + omega_lambda.
    density_sum = m_config.omega_radiation / a4 + m_config.omega_matter / a3 +
                  m_config.omega_curvature / a2 + m_config.omega_lambda;
  } else {
    // S01: E(a)^2 = (Omega_gamma + Omega_ur) a^-4 + Omega_cb a^-3 + sum_i Omega_ncdm,i(a)
    //               + Omega_k a^-2 + Omega_lambda.
    density_sum = (m_omega_gamma0 + m_omega_ur0) / a4 + m_omega_cb0 / a3 +
                  totalMassiveNeutrinoDensityParameter(scale_factor) +
                  m_config.omega_curvature / a2 + m_config.omega_lambda;
  }
  if (!std::isfinite(density_sum) || !(density_sum > 0.0)) {
    throw std::invalid_argument("density sum for H(a) must be finite and positive");
  }
  return std::sqrt(density_sum);
}

double LambdaCdmBackground::hubbleSi(double scale_factor) const {
  const double hubble = hubble0Si() * eFactor(scale_factor);
  if (!std::isfinite(hubble) || !(hubble > 0.0)) {
    throw std::runtime_error("H(a) is non-finite or non-positive");
  }
  return hubble;
}

double LambdaCdmBackground::criticalDensitySi(double scale_factor) const {
  const double h = hubbleSi(scale_factor);
  // rho_crit(a) = 3 H(a)^2 / (8 pi G).
  return (3.0 * h * h) / (8.0 * constants::k_pi * constants::k_newton_g_si);
}


double LambdaCdmBackground::cosmicTimeIntervalSi(
    double scale_factor_begin,
    double scale_factor_end) const {
  if (!std::isfinite(scale_factor_begin) || !std::isfinite(scale_factor_end) ||
      scale_factor_begin <= 0.0 || scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale factors must be finite and positive");
  }
  if (scale_factor_begin == scale_factor_end) {
    return 0.0;
  }
  const double sign = scale_factor_end > scale_factor_begin ? 1.0 : -1.0;
  const double left = std::min(scale_factor_begin, scale_factor_end);
  const double right = std::max(scale_factor_begin, scale_factor_end);
  const auto integrand = [this](double scale_factor) {
    return 1.0 / (scale_factor * hubbleSi(scale_factor));
  };
  // Relative target is tightened below double-precision needs of the scheduler.
  const double estimate = (right - left) * integrand(0.5 * (left + right));
  const double tolerance = std::max(std::abs(estimate) * 1.0e-11, 1.0e-6);
  return sign * adaptiveSimpson(integrand, left, right, tolerance);
}

double LambdaCdmBackground::cosmicTimeSinceScaleFactorSi(
    double scale_factor,
    double lower_scale_factor) const {
  if (scale_factor < lower_scale_factor) {
    throw std::invalid_argument("scale_factor must be >= lower_scale_factor");
  }
  return cosmicTimeIntervalSi(lower_scale_factor, scale_factor);
}

// ------------------------------------------------------------ S01 accounting

bool LambdaCdmBackground::radiationSpeciesEnabled() const noexcept {
  return m_config.radiation_species.enabled;
}

double LambdaCdmBackground::criticalDensity0Si() const {
  const double h0 = hubble0Si();
  return (3.0 * h0 * h0) / (8.0 * constants::k_pi * constants::k_newton_g_si);
}

double LambdaCdmBackground::photonDensityParameter0() const noexcept { return m_omega_gamma0; }

double LambdaCdmBackground::masslessNeutrinoDensityParameter0() const noexcept { return m_omega_ur0; }

std::size_t LambdaCdmBackground::massiveNeutrinoCount() const noexcept { return m_ncdm_y0.size(); }

namespace {
void requireEigenstate(std::size_t eigenstate, std::size_t count) {
  if (eigenstate >= count) {
    throw std::out_of_range("massive neutrino eigenstate index out of range");
  }
}
void requirePositiveScaleFactor(double scale_factor) {
  if (!std::isfinite(scale_factor) || !(scale_factor > 0.0)) {
    throw std::invalid_argument("scale_factor must be finite and positive");
  }
}
}  // namespace

double LambdaCdmBackground::massiveNeutrinoMassOverTemperature(std::size_t eigenstate, double scale_factor) const {
  requireEigenstate(eigenstate, m_ncdm_y0.size());
  requirePositiveScaleFactor(scale_factor);
  return m_ncdm_y0[eigenstate] * scale_factor;
}

double LambdaCdmBackground::massiveNeutrinoDensityParameter(std::size_t eigenstate, double scale_factor) const {
  const double y = massiveNeutrinoMassOverTemperature(eigenstate, scale_factor);
  const double a4 = std::pow(scale_factor, 4);
  return m_omega_ncdm_unit0 / a4 * fermiDiracEnergyIntegral(y) / energyMasslessMoment();
}

double LambdaCdmBackground::massiveNeutrinoPressureParameter(std::size_t eigenstate, double scale_factor) const {
  const double y = massiveNeutrinoMassOverTemperature(eigenstate, scale_factor);
  const double a4 = std::pow(scale_factor, 4);
  return m_omega_ncdm_unit0 / a4 * fermiDiracPressureIntegral(y) / energyMasslessMoment();
}

double LambdaCdmBackground::totalMassiveNeutrinoDensityParameter(double scale_factor) const {
  double total = 0.0;
  for (std::size_t i = 0; i < m_ncdm_y0.size(); ++i) total += massiveNeutrinoDensityParameter(i, scale_factor);
  return total;
}

double LambdaCdmBackground::totalMassiveNeutrinoPressureParameter(double scale_factor) const {
  double total = 0.0;
  for (std::size_t i = 0; i < m_ncdm_y0.size(); ++i) total += massiveNeutrinoPressureParameter(i, scale_factor);
  return total;
}

double LambdaCdmBackground::coldMatterDensityParameter0() const noexcept { return m_omega_cb0; }

double LambdaCdmBackground::massiveNeutrinoMeanMomentumOverMass(std::size_t eigenstate, double scale_factor) const {
  const double y = massiveNeutrinoMassOverTemperature(eigenstate, scale_factor);
  return fermiDiracMeanMomentumOverTemperature() / y;
}

double LambdaCdmBackground::matterRadiationEqualityScaleFactor() const {
  if (!radiationSpeciesEnabled()) {
    throw std::logic_error("matter-radiation equality requires radiation species to be enabled");
  }
  // Radiation-like part: photons + massless + 3 P_ncdm; matter-like part:
  // Omega_cb a^-3 + (rho_ncdm - 3 P_ncdm). Monotone difference in ln a.
  const auto excess = [this](double a) {
    const double a4 = std::pow(a, 4);
    const double p3 = 3.0 * totalMassiveNeutrinoPressureParameter(a);
    const double radiation = (m_omega_gamma0 + m_omega_ur0) / a4 + p3;
    const double matter = m_omega_cb0 / (a4 / a) + totalMassiveNeutrinoDensityParameter(a) - p3;
    return radiation - matter;
  };
  double ln_lo = std::log(1.0e-12);
  double ln_hi = 0.0;
  if (!(excess(std::exp(ln_lo)) > 0.0) || !(excess(1.0) < 0.0)) {
    throw std::runtime_error("matter-radiation equality is not bracketed in a in [1e-12, 1]");
  }
  for (int iteration = 0; iteration < 200 && (ln_hi - ln_lo) > 1.0e-14; ++iteration) {
    const double ln_mid = 0.5 * (ln_lo + ln_hi);
    if (excess(std::exp(ln_mid)) > 0.0) {
      ln_lo = ln_mid;
    } else {
      ln_hi = ln_mid;
    }
  }
  return std::exp(0.5 * (ln_lo + ln_hi));
}

// ------------------------------------------------------------ S01 CMB

double LambdaCdmBackground::cmbTemperatureK(double scale_factor) const {
  if (!radiationSpeciesEnabled()) {
    throw std::logic_error(
        "T_cmb(a) requested but the background has no radiation species (cosmology.radiation_background=none); "
        "no independent CMB temperature authority exists");
  }
  requirePositiveScaleFactor(scale_factor);
  return m_config.radiation_species.t_cmb0_k / scale_factor;
}

double LambdaCdmBackground::cmbEnergyDensitySi(double scale_factor) const {
  const double t = cmbTemperatureK(scale_factor);
  return constants::k_radiation_constant_si * std::pow(t, 4);
}

double LambdaCdmBackground::cmbPhotonNumberDensitySi(double scale_factor) const {
  const double t = cmbTemperatureK(scale_factor);
  const double k_over_hbar_c = constants::k_boltzmann_si * t /
                               (constants::k_reduced_planck_si * constants::k_speed_of_light_si);
  return 2.0 * constants::k_zeta_3 / (constants::k_pi * constants::k_pi) * k_over_hbar_c * k_over_hbar_c *
         k_over_hbar_c;
}

// ------------------------------------------------------------ S01 distances

double LambdaCdmBackground::conformalTimeIntervalSi(double scale_factor_begin, double scale_factor_end) const {
  requirePositiveScaleFactor(scale_factor_begin);
  requirePositiveScaleFactor(scale_factor_end);
  if (scale_factor_begin == scale_factor_end) return 0.0;
  const double sign = scale_factor_end > scale_factor_begin ? 1.0 : -1.0;
  const double left = std::log(std::min(scale_factor_begin, scale_factor_end));
  const double right = std::log(std::max(scale_factor_begin, scale_factor_end));
  // d eta = da / (a^2 H) = d ln a / (a H).
  const auto integrand = [this](double ln_a) {
    const double a = std::exp(ln_a);
    return 1.0 / (a * hubbleSi(a));
  };
  const double estimate = (right - left) * integrand(0.5 * (left + right));
  const double tolerance = std::max(std::abs(estimate) * 1.0e-12, 1.0e-9);
  return sign * adaptiveSimpson(integrand, left, right, tolerance);
}

double LambdaCdmBackground::comovingDistanceSi(double scale_factor) const {
  requirePositiveScaleFactor(scale_factor);
  if (scale_factor > 1.0) {
    throw std::invalid_argument("comoving distance is defined for past light cones (a <= 1)");
  }
  return constants::k_speed_of_light_si * conformalTimeIntervalSi(scale_factor, 1.0);
}

double LambdaCdmBackground::transverseComovingDistanceSi(double scale_factor) const {
  const double chi = comovingDistanceSi(scale_factor);
  const double omega_k = m_config.omega_curvature;
  if (omega_k == 0.0) return chi;
  const double hubble_distance = constants::k_speed_of_light_si / hubble0Si();
  const double root = std::sqrt(std::abs(omega_k));
  return omega_k > 0.0 ? hubble_distance / root * std::sinh(root * chi / hubble_distance)
                       : hubble_distance / root * std::sin(root * chi / hubble_distance);
}

double LambdaCdmBackground::angularDiameterDistanceSi(double scale_factor) const {
  return scale_factor * transverseComovingDistanceSi(scale_factor);
}

double LambdaCdmBackground::luminosityDistanceSi(double scale_factor) const {
  return transverseComovingDistanceSi(scale_factor) / scale_factor;
}

double LambdaCdmBackground::scaleFactorAtComovingDistanceSi(double comoving_distance_si, double a_min) const {
  if (!std::isfinite(comoving_distance_si) || comoving_distance_si < 0.0) {
    throw std::invalid_argument("comoving distance must be finite and >= 0");
  }
  requirePositiveScaleFactor(a_min);
  if (comoving_distance_si == 0.0) return 1.0;
  if (comoving_distance_si > comovingDistanceSi(a_min)) {
    throw std::out_of_range("requested comoving distance lies beyond a_min (no extrapolation)");
  }
  double ln_lo = std::log(a_min);  // chi(lo) >= target
  double ln_hi = 0.0;              // chi(hi) = 0 < target
  for (int iteration = 0; iteration < 200 && (ln_hi - ln_lo) > 1.0e-15; ++iteration) {
    const double ln_mid = 0.5 * (ln_lo + ln_hi);
    if (comovingDistanceSi(std::exp(ln_mid)) > comoving_distance_si) {
      ln_lo = ln_mid;
    } else {
      ln_hi = ln_mid;
    }
  }
  return std::exp(0.5 * (ln_lo + ln_hi));
}

}  // namespace cosmosim::core
