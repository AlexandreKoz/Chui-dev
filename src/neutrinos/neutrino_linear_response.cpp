// CHUI-EXT-PHYS-REF-v1 N03 massive-neutrino PM source models (source only,
// NOT compiled or executed in pass 01).

#include "cosmosim/neutrinos/neutrino_linear_response.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::neutrinos {
namespace {

namespace k = core::constants;

constexpr double k_char_table_y_max = 50.0;
constexpr double k_char_table_step = 0.01;
constexpr double k_char_x_max = 40.0;  // e^-40 x^4 below 1e-11 of the normalization

[[nodiscard]] double fdOccupation(double x) { return x > 40.0 ? std::exp(-x) / (1.0 + std::exp(-x)) : 1.0 / (std::exp(x) + 1.0); }

struct Gauss16 {
  std::array<double, 16> x{};
  std::array<double, 16> w{};
};

[[nodiscard]] Gauss16 buildGauss16() {
  Gauss16 rule;
  constexpr int n = 16;
  for (int i = 0; i < n; ++i) {
    double x = std::cos(k::k_pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
    double dp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0;
      double p1 = x;
      for (int kk = 2; kk <= n; ++kk) {
        const double p2 = ((2.0 * kk - 1.0) * x * p1 - (kk - 1.0) * p0) / static_cast<double>(kk);
        p0 = p1;
        p1 = p2;
      }
      dp = static_cast<double>(n) * (x * p1 - p0) / (x * x - 1.0);
      const double step = p1 / dp;
      x -= step;
      if (std::abs(step) < 1.0e-16) break;
    }
    rule.x[static_cast<std::size_t>(i)] = x;
    rule.w[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - x * x) * dp * dp);
  }
  return rule;
}

[[nodiscard]] const Gauss16& gauss16() {
  static const Gauss16 rule = buildGauss16();
  return rule;
}

[[nodiscard]] double sinc(double t) {
  if (std::abs(t) < 1.0e-4) {
    const double t2 = t * t;
    return 1.0 - t2 / 6.0 + t2 * t2 / 120.0;
  }
  return std::sin(t) / t;
}

[[nodiscard]] double characteristicQuadrature(double y) {
  // Panel width <= 1/(1 + y) keeps >= 16 nodes per sin(x y) period.
  const auto panels = static_cast<std::size_t>(std::ceil(k_char_x_max * (1.0 + y)));
  const double width = k_char_x_max / static_cast<double>(panels);
  const Gauss16& rule = gauss16();
  double sum = 0.0;
  for (std::size_t p = 0; p < panels; ++p) {
    const double mid = width * (static_cast<double>(p) + 0.5);
    double panel = 0.0;
    for (std::size_t i = 0; i < 16U; ++i) {
      const double x = mid + 0.5 * width * rule.x[i];
      panel += rule.w[i] * x * x * fdOccupation(x) * sinc(x * y);
    }
    sum += 0.5 * width * panel;
  }
  return sum / (1.5 * k::k_zeta_3);
}

[[nodiscard]] const std::vector<double>& characteristicTable() {
  static const std::vector<double> table = [] {
    const auto n = static_cast<std::size_t>(std::llround(k_char_table_y_max / k_char_table_step)) + 1U;
    std::vector<double> values(n);
    for (std::size_t i = 0; i < n; ++i) values[i] = characteristicQuadrature(k_char_table_step * static_cast<double>(i));
    return values;
  }();
  return table;
}

void appendBytes(std::vector<std::byte>& out, const void* data, std::size_t bytes) {
  const auto* begin = static_cast<const std::byte*>(data);
  out.insert(out.end(), begin, begin + bytes);
}

void appendU64(std::vector<std::byte>& out, std::uint64_t value) { appendBytes(out, &value, sizeof(value)); }
void appendF64(std::vector<std::byte>& out, double value) { appendBytes(out, &value, sizeof(value)); }
void appendVector(std::vector<std::byte>& out, const std::vector<double>& values) {
  appendU64(out, values.size());
  if (!values.empty()) appendBytes(out, values.data(), values.size() * sizeof(double));
}

struct Reader {
  std::span<const std::byte> bytes;
  std::size_t offset = 0;
  void read(void* out, std::size_t n) {
    if (offset + n > bytes.size()) throw std::runtime_error(std::string(k_neutrino_response_module) + " sidecar truncated");
    std::memcpy(out, bytes.data() + offset, n);
    offset += n;
  }
  std::uint64_t u64() {
    std::uint64_t v = 0;
    read(&v, sizeof(v));
    return v;
  }
  double f64() {
    double v = 0.0;
    read(&v, sizeof(v));
    return v;
  }
  std::vector<double> vec(std::size_t max_size) {
    const std::uint64_t n = u64();
    if (n > max_size) throw std::runtime_error(std::string(k_neutrino_response_module) + " sidecar array too large");
    std::vector<double> v(static_cast<std::size_t>(n));
    if (n > 0) read(v.data(), static_cast<std::size_t>(n) * sizeof(double));
    return v;
  }
};

[[nodiscard]] std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
  return h;
}

}  // namespace

double fermiDiracCharacteristicFunctionQuadrature(double y) {
  if (!std::isfinite(y) || y < 0.0) throw std::invalid_argument("characteristic function needs y >= 0");
  return characteristicQuadrature(y);
}

double fermiDiracCharacteristicFunction(double y) {
  if (!std::isfinite(y) || y < 0.0) throw std::invalid_argument("characteristic function needs y >= 0");
  if (y > k_char_table_y_max) {
    // int x f sin(xy) dx = 1/(2y^3) + 1/(2y^5) + 3/(2y^7) + O(y^-9) from the
    // Taylor series of f = 1/(e^x+1) at 0; F = that / (y (3/2) zeta(3)).
    const double y2 = y * y;
    const double y4 = y2 * y2;
    return (0.5 / y4 + 0.5 / (y4 * y2) + 1.5 / (y4 * y4)) / (1.5 * k::k_zeta_3);
  }
  const std::vector<double>& table = characteristicTable();
  const double u = y / k_char_table_step;
  const auto last = static_cast<std::ptrdiff_t>(table.size()) - 1;
  auto j = static_cast<std::ptrdiff_t>(std::floor(u));
  j = std::clamp<std::ptrdiff_t>(j, 1, last - 2);
  const double t = u - static_cast<double>(j);
  const double lm1 = -t * (t - 1.0) * (t - 2.0) / 6.0;
  const double l0 = (t + 1.0) * (t - 1.0) * (t - 2.0) / 2.0;
  const double l1 = -(t + 1.0) * t * (t - 2.0) / 2.0;
  const double l2 = (t + 1.0) * t * (t - 1.0) / 6.0;
  const auto at = [&table](std::ptrdiff_t i) { return table[static_cast<std::size_t>(i)]; };
  return lm1 * at(j - 1) + l0 * at(j) + l1 * at(j + 1) + l2 * at(j + 2);
}

double superconformalTimeIntervalSi(const core::LambdaCdmBackground& background, double a0, double a1) {
  if (!(a0 > 0.0) || !(a1 > 0.0)) throw std::invalid_argument("superconformal time needs a > 0");
  if (a0 == a1) return 0.0;
  // s = int d ln a / (a^2 H), composite 16-point Gauss-Legendre in ln a with
  // panels of width <= 0.01 (smooth integrand; relative error < 1e-14).
  const double l0 = std::log(a0);
  const double l1 = std::log(a1);
  const auto panels = static_cast<std::size_t>(std::max(1.0, std::ceil(std::abs(l1 - l0) / 0.01)));
  const double width = (l1 - l0) / static_cast<double>(panels);
  const Gauss16& rule = gauss16();
  double sum = 0.0;
  for (std::size_t p = 0; p < panels; ++p) {
    const double mid = l0 + width * (static_cast<double>(p) + 0.5);
    double panel = 0.0;
    for (std::size_t i = 0; i < 16U; ++i) {
      const double a = std::exp(mid + 0.5 * width * rule.x[i]);
      panel += rule.w[i] / (a * a * background.hubbleSi(a));
    }
    sum += 0.5 * width * panel;
  }
  return sum;
}

NeutrinoSpectralSource::NeutrinoSpectralSource(const core::LambdaCdmBackground& background,
                                               spectral::BoltzmannTransferSet tables, NeutrinoSourceSettings settings,
                                               const parallel::MpiContext* mpi, core::MemoryGovernor* governor)
    : m_background(background), m_tables(std::move(tables)), m_settings(settings), m_mpi(mpi) {
  using Method = core::NeutrinoMethodConfig;
  if (m_settings.method != Method::kTransferRatio && m_settings.method != Method::kLinearResponse) {
    throw std::invalid_argument("NeutrinoSpectralSource supports transfer_ratio and linear_response only");
  }
  if (m_background.massiveNeutrinoCount() == 0U) {
    throw std::invalid_argument("NeutrinoSpectralSource needs at least one massive eigenstate");
  }
  if (m_settings.k_bins < 8U || !(m_settings.k_min_h_mpc > 0.0) || !(m_settings.k_max_h_mpc > m_settings.k_min_h_mpc) ||
      !(m_settings.hubble_param > 0.0)) {
    throw std::invalid_argument("NeutrinoSpectralSource: invalid k binning");
  }
  if (m_tables.commonKMin() > m_settings.k_min_h_mpc || m_tables.commonKMax() < m_settings.k_max_h_mpc) {
    throw std::invalid_argument("Boltzmann tables do not cover the configured response k range (no extrapolation)");
  }
  const auto& masses = m_background.config().radiation_species.ncdm_masses_ev;
  if (m_settings.method == Method::kLinearResponse) {
    for (const double m : masses) {
      if (std::abs(m - masses.front()) > 1.0e-12 * masses.front()) {
        throw std::invalid_argument(
            "linear_response requires degenerate eigenstate masses (one Fermi-Dirac velocity scale); "
            "non-degenerate response kernels are not implemented");
      }
    }
    if (m_settings.history_capacity < 16U) throw std::invalid_argument("history capacity must be >= 16");
    const double ratio = m_background.config().radiation_species.ncdm_temperature_ratio;
    const double t_nu0 = ratio * m_background.config().radiation_species.t_cmb0_k;
    m_velocity_scale_m_s = k::k_boltzmann_si * t_nu0 * k::k_speed_of_light_si / (masses.front() * k::k_electron_volt_si);
  }
  m_pm_k_to_si = 1.0 / k::k_megaparsec_si;
  const std::size_t bins = m_settings.k_bins;
  const double ln_lo = std::log(m_settings.k_min_h_mpc);
  const double ln_hi = std::log(m_settings.k_max_h_mpc);
  m_bin_ln_edges.resize(bins + 1U);
  m_bin_k_h_mpc.resize(bins);
  for (std::size_t i = 0; i <= bins; ++i) {
    m_bin_ln_edges[i] = ln_lo + (ln_hi - ln_lo) * static_cast<double>(i) / static_cast<double>(bins);
  }
  for (std::size_t i = 0; i < bins; ++i) m_bin_k_h_mpc[i] = std::exp(0.5 * (m_bin_ln_edges[i] + m_bin_ln_edges[i + 1U]));
  m_current_ratio.assign(bins, 0.0);
  m_current_multiplier.assign(bins, 1.0);
  if (m_settings.method == Method::kLinearResponse) {
    // Bounded authoritative history: reserved once (capacity x (s, a, A, R)).
    const std::uint64_t entry_bytes = (2U + 2U * bins) * sizeof(double);
    const std::uint64_t bytes = entry_bytes * m_settings.history_capacity;
    if (governor != nullptr) {
      m_reservation = governor->reserve(core::MemoryClass::kCanonicalPersistent, bytes, k_neutrino_response_module);
    }
    m_history_s.reserve(m_settings.history_capacity);
    m_history_a.reserve(m_settings.history_capacity);
    m_history_amplitude.reserve(m_settings.history_capacity * bins);
    m_history_ratio.reserve(m_settings.history_capacity * bins);
    if (m_reservation.pending()) m_reservation.commit();
  }
}

std::unique_ptr<NeutrinoSpectralSource> NeutrinoSpectralSource::fromConfig(const core::SimulationConfig& config,
                                                                           const parallel::MpiContext* mpi,
                                                                           core::MemoryGovernor* governor) {
  const core::NeutrinoPerturbationConfig& nu = config.extended.neutrinos;
  if (nu.method != core::NeutrinoMethodConfig::kTransferRatio && nu.method != core::NeutrinoMethodConfig::kLinearResponse) {
    return nullptr;
  }
  const core::LambdaCdmBackground background(core::cosmologyBackgroundConfigFromSimulation(config));
  std::vector<spectral::BoltzmannTransferEpoch> epochs;
  for (std::size_t i = 0; i < nu.table_paths.size(); ++i) {
    spectral::BoltzmannLoadRequest request;
    request.path = nu.table_paths[i];
    request.expected_sha256 = nu.table_sha256.at(i);
    request.format = nu.table_format;
    request.declared_gauge = nu.table_gauge;
    request.declared_redshift = nu.table_redshifts.at(i);
    request.omega_baryon = config.cosmology.omega_baryon;
    epochs.push_back(spectral::loadBoltzmannTransferEpoch(request, background));
  }
  if (nu.method == core::NeutrinoMethodConfig::kLinearResponse) {
    for (const auto& epoch : epochs) {
      if (!epoch.hasVelocities()) {
        throw std::runtime_error("neutrinos.method=linear_response: table '" + epoch.path +
                                 "' has no Newtonian-gauge theta_cb/theta_ncdm (initial velocity term)");
      }
    }
  }
  NeutrinoSourceSettings settings;
  settings.method = nu.method;
  settings.k_bins = nu.response_k_bins;
  settings.history_capacity = nu.response_history_capacity;
  settings.k_min_h_mpc = nu.response_k_min_h_mpc;
  settings.k_max_h_mpc = nu.response_k_max_h_mpc;
  settings.hubble_param = config.cosmology.hubble_param;
  settings.omega_baryon = config.cosmology.omega_baryon;
  return std::make_unique<NeutrinoSpectralSource>(background, spectral::BoltzmannTransferSet(std::move(epochs)),
                                                  settings, mpi, governor);
}

double NeutrinoSpectralSource::densityFraction(double scale_factor) const {
  const double omega_nu = m_background.totalMassiveNeutrinoDensityParameter(scale_factor);
  const double omega_cb = m_background.coldMatterDensityParameter0() / (scale_factor * scale_factor * scale_factor);
  return omega_nu / omega_cb;
}

double NeutrinoSpectralSource::superconformalTimeSi(double scale_factor) const {
  // s measured from the first response epoch (only differences enter).
  return superconformalTimeIntervalSi(m_background, m_initial_scale_factor, scale_factor);
}

std::size_t NeutrinoSpectralSource::binIndex(double k_pm) const noexcept {
  if (!(k_pm > 0.0)) return k_no_bin;
  const double ln_k = std::log(k_pm / m_settings.hubble_param);
  if (ln_k < m_bin_ln_edges.front() || ln_k >= m_bin_ln_edges.back()) return k_no_bin;
  const double width = (m_bin_ln_edges.back() - m_bin_ln_edges.front()) / static_cast<double>(m_bin_k_h_mpc.size());
  const auto bin = static_cast<std::size_t>((ln_k - m_bin_ln_edges.front()) / width);
  return std::min(bin, m_bin_k_h_mpc.size() - 1U);
}

void NeutrinoSpectralSource::requireWavenumberDomain(double k_min_pm, double k_max_pm) const {
  const double lo = k_min_pm / m_settings.hubble_param;
  const double hi = k_max_pm / m_settings.hubble_param;
  if (lo < m_settings.k_min_h_mpc || hi > m_settings.k_max_h_mpc) {
    throw std::runtime_error("PM mesh wavenumbers [" + std::to_string(lo) + ", " + std::to_string(hi) +
                             "] h/Mpc exceed neutrinos.response_k_min/max_h_mpc [" +
                             std::to_string(m_settings.k_min_h_mpc) + ", " + std::to_string(m_settings.k_max_h_mpc) +
                             "] (no extrapolation of the neutrino source)");
  }
}

bool NeutrinoSpectralSource::needsObservation() const noexcept {
  return m_settings.method == core::NeutrinoMethodConfig::kLinearResponse;
}

void NeutrinoSpectralSource::beginSolve(double scale_factor) {
  if (!std::isfinite(scale_factor) || !(scale_factor > 0.0)) throw std::invalid_argument("invalid PM scale factor");
  m_solve_scale_factor = scale_factor;
  m_solve_open = true;
  if (m_settings.method != core::NeutrinoMethodConfig::kLinearResponse) return;
  const std::size_t bins = m_bin_k_h_mpc.size();
  // Deterministic replay after a rejected step: entries at a >= current are discarded.
  while (!m_history_a.empty() && m_history_a.back() >= scale_factor) {
    m_history_a.pop_back();
    m_history_s.pop_back();
    m_history_amplitude.resize(m_history_amplitude.size() - bins);
    m_history_ratio.resize(m_history_ratio.size() - bins);
  }
  if (m_history_a.empty()) {
    m_initial_scale_factor = scale_factor;
    m_initial_delta_ratio.assign(bins, 0.0);
    m_initial_theta_ratio_s.assign(bins, 0.0);
    // theta_s = a theta_conformal; CLASS theta in 1/Mpc with c = 1 -> c / Mpc in SI.
    const double theta_unit = scale_factor * k::k_speed_of_light_si / k::k_megaparsec_si;
    for (std::size_t b = 0; b < bins; ++b) {
      m_initial_delta_ratio[b] = m_tables.ncdmToColdRatio(m_bin_k_h_mpc[b], scale_factor);
      m_initial_theta_ratio_s[b] = theta_unit * m_tables.thetaNcdmOverDeltaCold(m_bin_k_h_mpc[b], scale_factor);
    }
  }
}

void NeutrinoSpectralSource::finalizeObservation(std::span<const double> local_bin_power,
                                                 std::span<const double> local_bin_modes) {
  if (!m_solve_open) throw std::logic_error("finalizeObservation without beginSolve");
  m_solve_open = false;
  const double a = m_solve_scale_factor;
  if (m_settings.method == core::NeutrinoMethodConfig::kTransferRatio) {
    computeTransferRatio(a);
    return;
  }
  const std::size_t bins = m_bin_k_h_mpc.size();
  if (local_bin_power.size() != bins || local_bin_modes.size() != bins) {
    throw std::logic_error("neutrino response observation has the wrong bin count");
  }
  std::vector<double> power(local_bin_power.begin(), local_bin_power.end());
  std::vector<double> modes(local_bin_modes.begin(), local_bin_modes.end());
  if (m_mpi != nullptr && m_mpi->isEnabled() && m_mpi->worldSize() > 1) {
    m_mpi->allreduceSumDoublesInPlace(power);
    m_mpi->allreduceSumDoublesInPlace(modes);
  }
  // rms amplitude per bin; bins without mesh modes interpolate ln A in ln k
  // between observed bins (growth ratios are smooth in k) and take the nearest
  // observed bin beyond the observed range. At least one bin must be observed.
  std::vector<double> amplitude(bins, 0.0);
  std::vector<std::size_t> observed;
  for (std::size_t b = 0; b < bins; ++b) {
    if (modes[b] > 0.0 && power[b] > 0.0) {
      amplitude[b] = std::sqrt(power[b] / modes[b]);
      observed.push_back(b);
    }
  }
  if (observed.empty()) {
    throw std::runtime_error("neutrino linear response: no PM mode fell into the response k bins");
  }
  for (std::size_t b = 0; b < bins; ++b) {
    if (amplitude[b] > 0.0) continue;
    const auto upper = std::lower_bound(observed.begin(), observed.end(), b);
    if (upper == observed.begin()) {
      amplitude[b] = amplitude[observed.front()];
    } else if (upper == observed.end()) {
      amplitude[b] = amplitude[observed.back()];
    } else {
      const std::size_t hi = *upper;
      const std::size_t lo = *(upper - 1);
      const double t = (std::log(m_bin_k_h_mpc[b]) - std::log(m_bin_k_h_mpc[lo])) /
                       (std::log(m_bin_k_h_mpc[hi]) - std::log(m_bin_k_h_mpc[lo]));
      amplitude[b] = std::exp((1.0 - t) * std::log(amplitude[lo]) + t * std::log(amplitude[hi]));
    }
  }
  computeLinearResponse(a, amplitude);
}

void NeutrinoSpectralSource::computeTransferRatio(double a) {
  const double f = densityFraction(a);
  for (std::size_t b = 0; b < m_bin_k_h_mpc.size(); ++b) {
    m_current_ratio[b] = m_tables.ncdmToColdRatio(m_bin_k_h_mpc[b], a);
    m_current_multiplier[b] = 1.0 + f * m_current_ratio[b];
  }
}

double NeutrinoSpectralSource::evaluateLinearResponse(std::size_t bin, double s_now, double a_now,
                                                      std::span<const double> s_nodes,
                                                      std::span<const double> a_nodes,
                                                      std::span<const double> amplitude_over_now,
                                                      std::span<const double> ratio_nodes) const {
  if (s_nodes.empty() || s_nodes.size() != a_nodes.size() || s_nodes.size() != amplitude_over_now.size() ||
      s_nodes.size() != ratio_nodes.size()) {
    throw std::invalid_argument("linear response history arrays are inconsistent");
  }
  static_cast<void>(a_now);
  const double k_si = m_bin_k_h_mpc.at(bin) * m_settings.hubble_param * m_pm_k_to_si;
  const double ku0 = k_si * m_velocity_scale_m_s;
  const double rho_crit0 = m_background.criticalDensity0Si();
  const double omega_cb0 = m_background.coldMatterDensityParameter0();
  const double four_pi_g = 4.0 * k::k_pi * k::k_newton_g_si;
  // Initial free-streaming term (node 0 is the initial epoch s_i).
  const double tau0 = s_now - s_nodes.front();
  double ratio = (m_initial_delta_ratio.at(bin) - m_initial_theta_ratio_s.at(bin) * tau0) *
                 fermiDiracCharacteristicFunction(ku0 * tau0) * amplitude_over_now.front();
  // Gravitational response: trapezoid over the history nodes plus the current
  // endpoint, where the (s - s') weight vanishes.
  const auto integrand = [&](std::size_t j) {
    const double a = a_nodes[j];
    const double lag = s_now - s_nodes[j];
    const double rho_cb = omega_cb0 * rho_crit0 / (a * a * a);
    const double rho_nu = m_background.totalMassiveNeutrinoDensityParameter(a) * rho_crit0;
    const double a4 = a * a * a * a;
    return lag * fermiDiracCharacteristicFunction(ku0 * lag) * four_pi_g * a4 * (rho_cb + rho_nu * ratio_nodes[j]) *
           amplitude_over_now[j];
  };
  double integral = 0.0;
  for (std::size_t j = 0; j + 1U < s_nodes.size(); ++j) {
    integral += 0.5 * (s_nodes[j + 1U] - s_nodes[j]) * (integrand(j) + integrand(j + 1U));
  }
  integral += 0.5 * (s_now - s_nodes.back()) * integrand(s_nodes.size() - 1U);  // endpoint value 0
  ratio += integral;
  return ratio;
}

void NeutrinoSpectralSource::computeLinearResponse(double a, std::span<const double> amplitude) {
  const std::size_t bins = m_bin_k_h_mpc.size();
  if (m_history_s.size() >= m_settings.history_capacity) {
    throw std::runtime_error("neutrino linear-response history reached neutrinos.response_history_capacity = " +
                             std::to_string(m_settings.history_capacity) +
                             " entries (refused; entries are never dropped). Raise the capacity.");
  }
  const double s_now = superconformalTimeSi(a);
  const std::size_t n = m_history_s.size();
  std::vector<double> amp_ratio(n);
  std::vector<double> ratios(n);
  for (std::size_t b = 0; b < bins; ++b) {
    double r = m_initial_delta_ratio[b];
    if (n > 0U) {
      for (std::size_t j = 0; j < n; ++j) {
        amp_ratio[j] = m_history_amplitude[j * bins + b] / amplitude[b];
        ratios[j] = m_history_ratio[j * bins + b];
      }
      r = evaluateLinearResponse(b, s_now, a, m_history_s, m_history_a, amp_ratio, ratios);
    }
    m_current_ratio[b] = r;
  }
  m_history_s.push_back(s_now);
  m_history_a.push_back(a);
  m_history_amplitude.insert(m_history_amplitude.end(), amplitude.begin(), amplitude.end());
  m_history_ratio.insert(m_history_ratio.end(), m_current_ratio.begin(), m_current_ratio.end());
  const double f = densityFraction(a);
  for (std::size_t b = 0; b < bins; ++b) m_current_multiplier[b] = 1.0 + f * m_current_ratio[b];
}

double NeutrinoSpectralSource::multiplier(double k_pm) const noexcept {
  const double ln_k = std::log(k_pm / m_settings.hubble_param);
  const std::size_t bins = m_bin_k_h_mpc.size();
  const double ln_first = std::log(m_bin_k_h_mpc.front());
  const double ln_last = std::log(m_bin_k_h_mpc.back());
  if (!(ln_k > ln_first)) return m_current_multiplier.front();
  if (!(ln_k < ln_last)) return m_current_multiplier.back();
  const double width = (ln_last - ln_first) / static_cast<double>(bins - 1U);
  const double u = (ln_k - ln_first) / width;
  const auto lo = std::min(static_cast<std::size_t>(u), bins - 2U);
  const double t = u - static_cast<double>(lo);
  return (1.0 - t) * m_current_multiplier[lo] + t * m_current_multiplier[lo + 1U];
}

std::uint64_t NeutrinoSpectralSource::fingerprint() const noexcept {
  std::uint64_t h = 0x6e75726573706f6eULL;
  h = mix(h, static_cast<std::uint64_t>(m_settings.method));
  h = mix(h, m_bin_k_h_mpc.size());
  h = mix(h, m_history_s.size());
  if (!m_history_s.empty()) h = mix(h, std::bit_cast<std::uint64_t>(m_history_s.back()));
  h = mix(h, std::bit_cast<std::uint64_t>(m_initial_scale_factor));
  return h;
}

core::ModuleSidecarBlock NeutrinoSpectralSource::serialize() const {
  core::ModuleSidecarBlock block;
  block.module_name = k_neutrino_response_module;
  block.schema_version = k_neutrino_response_schema;
  std::vector<std::byte>& out = block.payload;
  appendU64(out, static_cast<std::uint64_t>(m_settings.method));
  appendU64(out, m_bin_k_h_mpc.size());
  appendF64(out, m_settings.k_min_h_mpc);
  appendF64(out, m_settings.k_max_h_mpc);
  appendF64(out, m_initial_scale_factor);
  appendVector(out, m_initial_delta_ratio);
  appendVector(out, m_initial_theta_ratio_s);
  appendVector(out, m_history_s);
  appendVector(out, m_history_a);
  appendVector(out, m_history_amplitude);
  appendVector(out, m_history_ratio);
  appendVector(out, m_current_ratio);
  appendVector(out, m_current_multiplier);
  return block;
}

void NeutrinoSpectralSource::restore(const core::ModuleSidecarBlock& block) {
  if (block.module_name != k_neutrino_response_module || block.schema_version != k_neutrino_response_schema ||
      block.particle_indexed) {
    throw std::runtime_error(std::string(k_neutrino_response_module) + " sidecar has an incompatible schema");
  }
  Reader reader{block.payload, 0U};
  const std::size_t bins = m_bin_k_h_mpc.size();
  if (reader.u64() != static_cast<std::uint64_t>(m_settings.method) || reader.u64() != bins ||
      reader.f64() != m_settings.k_min_h_mpc || reader.f64() != m_settings.k_max_h_mpc) {
    throw std::runtime_error(std::string(k_neutrino_response_module) +
                             " sidecar was written with a different neutrino method or k binning (restart refused)");
  }
  const std::size_t cap = m_settings.history_capacity;
  m_initial_scale_factor = reader.f64();
  m_initial_delta_ratio = reader.vec(bins);
  m_initial_theta_ratio_s = reader.vec(bins);
  std::vector<double> s = reader.vec(cap);
  std::vector<double> a = reader.vec(cap);
  std::vector<double> amp = reader.vec(cap * bins);
  std::vector<double> ratio = reader.vec(cap * bins);
  if (a.size() != s.size() || amp.size() != s.size() * bins || ratio.size() != s.size() * bins) {
    throw std::runtime_error(std::string(k_neutrino_response_module) + " sidecar history is inconsistent");
  }
  // Keep the reserved capacity (assign into the pre-reserved vectors).
  m_history_s.assign(s.begin(), s.end());
  m_history_a.assign(a.begin(), a.end());
  m_history_amplitude.assign(amp.begin(), amp.end());
  m_history_ratio.assign(ratio.begin(), ratio.end());
  m_current_ratio = reader.vec(bins);
  m_current_multiplier = reader.vec(bins);
  if (m_current_ratio.size() != bins || m_current_multiplier.size() != bins || reader.offset != block.payload.size()) {
    throw std::runtime_error(std::string(k_neutrino_response_module) + " sidecar has trailing or missing data");
  }
}

bool NeutrinoSpectralSource::restoreFrom(const core::ModuleSidecarRegistry& registry) {
  const core::ModuleSidecarBlock* block = registry.find(k_neutrino_response_module);
  if (block == nullptr) return false;
  restore(*block);
  return true;
}

}  // namespace cosmosim::neutrinos
