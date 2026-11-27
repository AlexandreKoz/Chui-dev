// CHUI-EXT-PHYS-REF-v1 S02: atomic line data and Voigt profiles (source only,
// NOT compiled or executed in pass 01).

#include "cosmosim/spectral/atomic_data.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::spectral {
namespace {

namespace k = core::constants;

constexpr int k_weideman_terms = 32;

struct WeidemanCoefficients {
  double l = 0.0;
  std::array<double, k_weideman_terms> a{};  // a_1 .. a_N
};

[[nodiscard]] WeidemanCoefficients buildWeideman() {
  // a_n = (1/(2M)) sum_{k=-M+1}^{M-1} f(t_k) cos(n theta_k), M = 2N,
  // theta_k = k pi / M, t_k = L tan(theta_k / 2), f(t) = exp(-t^2) (L^2 + t^2);
  // the k = +-M node has t = infinity and f = 0. f is even, so a_n is real.
  WeidemanCoefficients c;
  const int n_terms = k_weideman_terms;
  const int m = 2 * n_terms;
  c.l = std::sqrt(static_cast<double>(n_terms) / std::sqrt(2.0));
  for (int n = 1; n <= n_terms; ++n) {
    double sum = 0.0;
    for (int kk = -m + 1; kk <= m - 1; ++kk) {
      const double theta = static_cast<double>(kk) * k::k_pi / static_cast<double>(m);
      const double t = c.l * std::tan(0.5 * theta);
      const double f = std::exp(-t * t) * (c.l * c.l + t * t);
      sum += f * std::cos(static_cast<double>(n) * theta);
    }
    c.a[static_cast<std::size_t>(n - 1)] = sum / (2.0 * static_cast<double>(m));
  }
  return c;
}

[[nodiscard]] const WeidemanCoefficients& weideman() {
  static const WeidemanCoefficients coefficients = buildWeideman();
  return coefficients;
}

struct Gauss8 {
  std::array<double, 8> x{};
  std::array<double, 8> w{};
};

[[nodiscard]] Gauss8 buildGauss8() {
  Gauss8 rule;
  constexpr int n = 8;
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

[[nodiscard]] double gauss8(const std::function<double(double)>& f, double l, double r) {
  static const Gauss8 rule = buildGauss8();
  const double half = 0.5 * (r - l);
  const double mid = 0.5 * (r + l);
  double s = 0.0;
  for (std::size_t i = 0; i < 8U; ++i) s += rule.w[i] * f(mid + half * rule.x[i]);
  return half * s;
}

[[nodiscard]] double adaptiveGauss(const std::function<double(double)>& f, double l, double r, double whole,
                                   double tol, int depth) {
  const double mid = 0.5 * (l + r);
  const double left = gauss8(f, l, mid);
  const double right = gauss8(f, mid, r);
  if (std::abs(left + right - whole) <= tol || depth <= 0) {
    if (depth <= 0 && std::abs(left + right - whole) > tol) {
      throw std::runtime_error("Voigt reference quadrature did not converge");
    }
    return left + right;
  }
  return adaptiveGauss(f, l, mid, left, 0.5 * tol, depth - 1) + adaptiveGauss(f, mid, r, right, 0.5 * tol, depth - 1);
}

[[nodiscard]] AtomicLine makeLymanAlpha() {
  AtomicLine line;
  line.name = "HI_Lyman_alpha_1215";
  line.rest_wavelength_m = 1215.67e-10;
  line.rest_frequency_hz = k::k_speed_of_light_si / line.rest_wavelength_m;
  line.oscillator_strength = 0.4164;
  line.einstein_a_s = 6.265e8;
  line.lower_weight = 2.0;
  line.upper_weight = 6.0;
  line.atom_mass_kg = 1.00782503223 * k::k_atomic_mass_unit_si;
  line.provenance = "recalled_lead: NIST ASD H I 1s-2p (lambda_vac 1215.67 A, f 0.4164, A 6.265e8 s^-1); verify";
  return line;
}

[[nodiscard]] AtomicLine makeHydrogen21cm() {
  AtomicLine line;
  line.name = "HI_21cm_hyperfine";
  line.rest_frequency_hz = 1420.405751768e6;
  line.rest_wavelength_m = k::k_speed_of_light_si / line.rest_frequency_hz;
  line.einstein_a_s = 2.85e-15;
  line.lower_weight = 1.0;
  line.upper_weight = 3.0;
  // f_lu from A_ul: A = (2 pi e^2 nu^2 / (eps0 m_e c^3)) (g_l/g_u) f.
  const double nu = line.rest_frequency_hz;
  line.oscillator_strength = line.einstein_a_s * k::k_vacuum_permittivity_si * k::k_electron_mass_si *
                             std::pow(k::k_speed_of_light_si, 3) * line.upper_weight /
                             (2.0 * k::k_pi * k::k_elementary_charge_si * k::k_elementary_charge_si * nu * nu *
                              line.lower_weight);
  line.atom_mass_kg = 1.00782503223 * k::k_atomic_mass_unit_si;
  line.provenance = "recalled_lead: nu_10 = 1420.405751768 MHz, A_10 = 2.85e-15 s^-1 (1% literature spread); verify";
  return line;
}

}  // namespace

const AtomicLine& hydrogenLymanAlpha() {
  static const AtomicLine line = makeLymanAlpha();
  return line;
}

const AtomicLine& hydrogen21cm() {
  static const AtomicLine line = makeHydrogen21cm();
  return line;
}

double lineIntegratedCrossSectionSi(const AtomicLine& line) {
  return k::k_elementary_charge_si * k::k_elementary_charge_si * line.oscillator_strength /
         (4.0 * k::k_vacuum_permittivity_si * k::k_electron_mass_si * k::k_speed_of_light_si);
}

double dopplerWidthHz(const AtomicLine& line, double temperature_k, double turbulent_b_m_s) {
  if (!std::isfinite(temperature_k) || temperature_k < 0.0 || !std::isfinite(turbulent_b_m_s) || turbulent_b_m_s < 0.0) {
    throw std::invalid_argument("Doppler width requires finite T >= 0 and b_turb >= 0");
  }
  const double b2 = 2.0 * k::k_boltzmann_si * temperature_k / line.atom_mass_kg + turbulent_b_m_s * turbulent_b_m_s;
  if (!(b2 > 0.0)) throw std::invalid_argument("Doppler width vanishes (T = 0 and b_turb = 0)");
  return line.rest_frequency_hz / k::k_speed_of_light_si * std::sqrt(b2);
}

double naturalDampingParameter(const AtomicLine& line, double doppler_width_hz) {
  if (!(doppler_width_hz > 0.0)) throw std::invalid_argument("damping parameter needs Delta nu_D > 0");
  return line.einstein_a_s / (4.0 * k::k_pi * doppler_width_hz);
}

std::complex<double> faddeeva(std::complex<double> z) {
  if (z.imag() < 0.0) throw std::domain_error("faddeeva implemented for Im z >= 0 only");
  const WeidemanCoefficients& c = weideman();
  const std::complex<double> iz(-z.imag(), z.real());  // i z
  const std::complex<double> denominator = c.l - iz;
  const std::complex<double> big_z = (c.l + iz) / denominator;
  // Horner evaluation of sum_{n=0}^{N-1} a_{n+1} Z^n.
  std::complex<double> polynomial(0.0, 0.0);
  for (int n = k_weideman_terms - 1; n >= 0; --n) {
    polynomial = polynomial * big_z + c.a[static_cast<std::size_t>(n)];
  }
  return 2.0 * polynomial / (denominator * denominator) + 1.0 / (std::sqrt(k::k_pi) * denominator);
}

double voigtH(double a, double x) {
  if (!std::isfinite(a) || a < 0.0 || !std::isfinite(x)) {
    throw std::invalid_argument("voigtH requires finite a >= 0 and finite x");
  }
  if (a == 0.0) return std::exp(-x * x);
  return faddeeva(std::complex<double>(x, a)).real();
}

double voigtHQuadrature(double a, double x) {
  if (!std::isfinite(a) || !(a > 1.0e-6) || !std::isfinite(x)) {
    throw std::invalid_argument("voigtHQuadrature requires a > 1e-6 and finite x");
  }
  // Substituting y = x + a tan(phi): H = (1/pi) int_{-pi/2}^{pi/2} exp(-(x + a tan phi)^2) dphi.
  const std::function<double(double)> integrand = [a, x](double phi) {
    const double y = x + a * std::tan(phi);
    return std::exp(-y * y);
  };
  const double peak = std::atan(-x / a);  // integrand maximum (= 1)
  const double lo = -0.5 * k::k_pi;
  const double hi = 0.5 * k::k_pi;
  const double scale = std::max(std::exp(-x * x), a / (std::sqrt(k::k_pi) * (x * x + a * a)));
  const double tol = 1.0e-12 * scale;
  const double left = adaptiveGauss(integrand, lo, peak, gauss8(integrand, lo, peak), tol, 40);
  const double right = adaptiveGauss(integrand, peak, hi, gauss8(integrand, peak, hi), tol, 40);
  return (left + right) / k::k_pi;
}

double voigtProfilePerHz(double a, double x, double doppler_width_hz) {
  if (!(doppler_width_hz > 0.0)) throw std::invalid_argument("Voigt profile needs Delta nu_D > 0");
  return voigtH(a, x) / (std::sqrt(k::k_pi) * doppler_width_hz);
}

}  // namespace cosmosim::spectral
