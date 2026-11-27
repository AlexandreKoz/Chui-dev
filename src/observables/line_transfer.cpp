// CHUI-EXT-PHYS-REF-v1 L01 / RA01 line transfer (source only, NOT compiled or
// executed in pass 01).

#include "cosmosim/observables/line_transfer.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <sstream>
#include <stdexcept>

#include "core/internal/sha256.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/spectral/reference_fft.hpp"

namespace cosmosim::observables {
namespace {

namespace k = core::constants;

[[nodiscard]] double wrapDelta(double d, const VelocityAxis& axis) {
  if (!axis.periodic) return d;
  const double period = axis.u_max - axis.u_min;
  d -= period * std::floor(d / period + 0.5);
  return d;
}

[[nodiscard]] double profileAt(double b, double a, double dv) {
  const double x = dv / b;
  return spectral::voigtH(a, x) / (std::sqrt(k::k_pi) * b);
}

[[nodiscard]] std::string trim(const std::string& s) {
  std::size_t lo = 0;
  std::size_t hi = s.size();
  while (lo < hi && std::isspace(static_cast<unsigned char>(s[lo]))) ++lo;
  while (hi > lo && std::isspace(static_cast<unsigned char>(s[hi - 1]))) --hi;
  return s.substr(lo, hi - lo);
}

}  // namespace

double segmentAveragedProfile(const spectral::AtomicLine& line, LineProfileModel model, double temperature_k,
                              double u0, double u1, double u_p, const VelocityAxis& axis) {
  if (!(temperature_k > 0.0)) throw std::invalid_argument("line profile needs T > 0");
  const double b = std::sqrt(2.0 * k::k_boltzmann_si * temperature_k / line.atom_mass_kg);
  const double a = model == LineProfileModel::kVoigt ? line.einstein_a_s * line.rest_wavelength_m / (4.0 * k::k_pi * b)
                                                     : 0.0;
  const double span = u1 - u0;
  if (std::abs(span) < 1.0e-12 * b) return profileAt(b, a, wrapDelta(u_p - u0, axis));
  const auto subintervals = static_cast<std::size_t>(std::ceil(std::abs(span) / (0.25 * b)));
  const double h = span / static_cast<double>(subintervals);
  // 3-point Gauss-Legendre (nodes 0, +-sqrt(3/5); weights 8/9, 5/9).
  const double node = std::sqrt(3.0 / 5.0);
  double sum = 0.0;
  for (std::size_t s = 0; s < subintervals; ++s) {
    const double mid = u0 + h * (static_cast<double>(s) + 0.5);
    sum += (5.0 / 9.0) * profileAt(b, a, wrapDelta(u_p - (mid - 0.5 * h * node), axis)) +
           (8.0 / 9.0) * profileAt(b, a, wrapDelta(u_p - mid, axis)) +
           (5.0 / 9.0) * profileAt(b, a, wrapDelta(u_p - (mid + 0.5 * h * node), axis));
  }
  return 0.5 * sum / static_cast<double>(subintervals);  // (1/span) * sum * h/2
}

std::vector<double> lymanAlphaOpticalDepth(const std::vector<LineSegment>& segments, const VelocityAxis& axis,
                                           LineProfileModel model) {
  if (axis.pixels == 0U || !(axis.u_max > axis.u_min)) throw std::invalid_argument("invalid velocity axis");
  const spectral::AtomicLine& line = spectral::hydrogenLymanAlpha();
  const double sigma_lambda = spectral::lineIntegratedCrossSectionSi(line) * line.rest_wavelength_m;
  std::vector<double> tau(axis.pixels, 0.0);
  for (const LineSegment& seg : segments) {
    if (!(seg.absorber_density_m3 > 0.0) || !(seg.path_proper_m > 0.0)) continue;
    const double column = seg.absorber_density_m3 * seg.path_proper_m;
    for (std::size_t p = 0; p < axis.pixels; ++p) {
      tau[p] += column * sigma_lambda *
                segmentAveragedProfile(line, model, seg.temperature_k, seg.velocity_begin_m_s, seg.velocity_end_m_s,
                                       axis.center(p), axis);
    }
  }
  return tau;
}

ForestStatistics forestStatistics(const std::vector<double>& tau, const VelocityAxis& axis, std::size_t pdf_bins,
                                  double redshift, const std::vector<LineSegment>& segments) {
  if (tau.size() != axis.pixels || pdf_bins == 0U) throw std::invalid_argument("forest statistics input mismatch");
  ForestStatistics out;
  std::vector<double> flux(tau.size());
  for (std::size_t p = 0; p < tau.size(); ++p) {
    flux[p] = std::exp(-tau[p]);
    out.mean_flux += flux[p];
  }
  out.mean_flux /= static_cast<double>(tau.size());
  out.flux_pdf.assign(pdf_bins, 0.0);
  for (const double f : flux) {
    const auto bin = std::min(static_cast<std::size_t>(f * static_cast<double>(pdf_bins)), pdf_bins - 1U);
    out.flux_pdf[bin] += 1.0 / (static_cast<double>(flux.size()) / static_cast<double>(pdf_bins));
  }
  const double lambda_obs = spectral::hydrogenLymanAlpha().rest_wavelength_m * (1.0 + redshift);
  for (const double f : flux) out.equivalent_width_m += (1.0 - f) * lambda_obs * axis.width() / k::k_speed_of_light_si;
  for (const LineSegment& seg : segments) out.neutral_column_m2 += seg.absorber_density_m3 * seg.path_proper_m;
  if (spectral::isPowerOfTwo(axis.pixels) && out.mean_flux > 0.0) {
    std::vector<std::complex<double>> delta(axis.pixels);
    for (std::size_t p = 0; p < axis.pixels; ++p) delta[p] = flux[p] / out.mean_flux - 1.0;
    spectral::fft1d(delta, false);
    const double length = axis.u_max - axis.u_min;
    const double n = static_cast<double>(axis.pixels);
    for (std::size_t m = 1; m <= axis.pixels / 2U; ++m) {
      out.power_k_s_per_m.push_back(2.0 * k::k_pi * static_cast<double>(m) / length);
      out.power_m_per_s.push_back(length * std::norm(delta[m]) / (n * n));
    }
  }
  return out;
}

double spinTemperature(const SpinTemperatureInputs& in) {
  if (!(in.kinetic_temperature_k > 0.0) || !(in.radiation_temperature_k > 0.0) || in.collisional_coupling < 0.0 ||
      in.lyman_alpha_coupling < 0.0) {
    throw std::invalid_argument("spin temperature needs T_K, T_gamma > 0 and couplings >= 0");
  }
  const double t_c = in.color_temperature_k > 0.0 ? in.color_temperature_k : in.kinetic_temperature_k;
  const double numerator = 1.0 + in.collisional_coupling + in.lyman_alpha_coupling;
  const double denominator = 1.0 / in.radiation_temperature_k + in.collisional_coupling / in.kinetic_temperature_k +
                             in.lyman_alpha_coupling / t_c;
  return numerator / denominator;
}

CollisionalCouplingTable parseCollisionalCouplingTable(const std::string& bytes, const std::string& origin) {
  CollisionalCouplingTable table;
  std::string schema;
  std::istringstream stream(bytes);
  std::string line;
  while (std::getline(stream, line)) {
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t.front() == '#') {
      const auto colon = t.find(':');
      if (colon == std::string::npos) continue;
      const std::string key = trim(t.substr(1, colon - 1));
      if (key == "schema") schema = trim(t.substr(colon + 1));
      if (key == "reference") table.reference = trim(t.substr(colon + 1));
      continue;
    }
    std::istringstream row(t);
    double tk = 0.0, hh = 0.0, eh = 0.0, ph = 0.0;
    row >> tk >> hh >> eh >> ph;
    if (row.fail() || !(tk > 0.0) || hh < 0.0 || eh < 0.0 || ph < 0.0) {
      throw std::runtime_error("21 cm collision table '" + origin + "': malformed row '" + t + "'");
    }
    if (!table.temperature_k.empty() && !(tk > table.temperature_k.back())) {
      throw std::runtime_error("21 cm collision table '" + origin + "': temperatures must increase");
    }
    table.temperature_k.push_back(tk);
    table.kappa_hh.push_back(hh);
    table.kappa_eh.push_back(eh);
    table.kappa_ph.push_back(ph);
  }
  if (schema != "chui_21cm_collision_rates_v1" || table.reference.empty() || table.temperature_k.size() < 2U) {
    throw std::runtime_error("21 cm collision table '" + origin + "': schema/reference/rows missing");
  }
  table.sha256 = core::internal::sha256Hex(bytes);
  return table;
}

double CollisionalCouplingTable::coupling(double t_k, double t_gamma, double n_h, double n_e, double n_p) const {
  if (!(t_k >= temperature_k.front()) || !(t_k <= temperature_k.back())) {
    throw std::out_of_range("kinetic temperature outside the 21 cm collision table (no extrapolation)");
  }
  const auto it = std::upper_bound(temperature_k.begin(), temperature_k.end(), t_k);
  std::size_t hi = static_cast<std::size_t>(std::distance(temperature_k.begin(), it));
  hi = std::clamp<std::size_t>(hi, 1U, temperature_k.size() - 1U);
  const std::size_t lo = hi - 1U;
  const double t = std::log(t_k / temperature_k[lo]) / std::log(temperature_k[hi] / temperature_k[lo]);
  const auto interp = [&](const std::vector<double>& v) {
    if (v[lo] > 0.0 && v[hi] > 0.0) return std::exp((1.0 - t) * std::log(v[lo]) + t * std::log(v[hi]));
    return (1.0 - t) * v[lo] + t * v[hi];
  };
  const spectral::AtomicLine& line = spectral::hydrogen21cm();
  const double t_star = k::k_planck_si * line.rest_frequency_hz / k::k_boltzmann_si;
  const double rate = 1.0e-6 * (interp(kappa_hh) * n_h + interp(kappa_eh) * n_e + interp(kappa_ph) * n_p);
  return t_star / (line.einstein_a_s * t_gamma) * rate;
}

double lymanAlphaCoupling(double p_alpha, double t_gamma) {
  if (!(p_alpha >= 0.0) || !(t_gamma > 0.0)) throw std::invalid_argument("x_alpha needs P_alpha >= 0, T_gamma > 0");
  const spectral::AtomicLine& line = spectral::hydrogen21cm();
  const double t_star = k::k_planck_si * line.rest_frequency_hz / k::k_boltzmann_si;
  return 4.0 * p_alpha * t_star / (27.0 * line.einstein_a_s * t_gamma);
}

std::vector<double> hiBrightnessTemperature(const std::vector<LineSegment>& segments, const VelocityAxis& axis,
                                            HiTransferModel model, double background_k) {
  if (axis.pixels == 0U) throw std::invalid_argument("invalid velocity axis");
  const spectral::AtomicLine& line = spectral::hydrogen21cm();
  const double nu0 = line.rest_frequency_hz;
  const double coefficient = 3.0 * std::pow(k::k_speed_of_light_si, 3) * line.einstein_a_s * k::k_planck_si /
                             (32.0 * k::k_pi * nu0 * nu0 * k::k_boltzmann_si);
  std::vector<double> t_b(axis.pixels, background_k);
  for (std::size_t p = 0; p < axis.pixels; ++p) {
    for (auto it = segments.rbegin(); it != segments.rend(); ++it) {
      const LineSegment& seg = *it;
      if (!(seg.absorber_density_m3 > 0.0) || !(seg.path_proper_m > 0.0)) continue;
      if (!(seg.spin_temperature_k > 0.0)) throw std::invalid_argument("21 cm transfer needs T_s > 0 per segment");
      const double phi = segmentAveragedProfile(line, LineProfileModel::kDopplerOnly, seg.temperature_k,
                                                seg.velocity_begin_m_s, seg.velocity_end_m_s, axis.center(p), axis);
      const double dtau = coefficient * seg.absorber_density_m3 * seg.path_proper_m * phi / seg.spin_temperature_k;
      if (model == HiTransferModel::kFullOpticalDepth) {
        const double absorbed = -std::expm1(-dtau);
        t_b[p] = t_b[p] * (1.0 - absorbed) + seg.spin_temperature_k * absorbed;
      } else {
        t_b[p] += dtau * (seg.spin_temperature_k - background_k);
      }
    }
  }
  return t_b;
}

}  // namespace cosmosim::observables
