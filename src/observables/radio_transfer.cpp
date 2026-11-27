// CHUI-EXT-PHYS-REF-v1 RA03-RA05 radio continuum, Faraday and Stokes transfer
// (source only, NOT compiled or executed in pass 01).

#include "cosmosim/observables/radio_transfer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::observables {
namespace {

namespace k = core::constants;

// Gaussian-CGS constants derived from the SI-defining values.
constexpr double k_c_cgs = k::k_speed_of_light_si * 1.0e2;
constexpr double k_q_esu = k::k_elementary_charge_si * k_c_cgs / 10.0;
constexpr double k_me_cgs = k::k_electron_mass_si * 1.0e3;
constexpr double k_kb_cgs = k::k_boltzmann_si * 1.0e7;
constexpr double k_h_cgs = k::k_planck_si * 1.0e7;
constexpr double k_cgs_emissivity_to_si = 0.1;   // erg s^-1 cm^-3 -> W m^-3
constexpr double k_cgs_absorption_to_si = 100.0;  // cm^-1 -> m^-1

struct KernelTable {
  double ln_x_min = 0.0;
  double step = 0.0;
  std::vector<double> ln_f;
  std::vector<double> ln_g;
};

[[nodiscard]] double kernelIntegral(double x, double order, bool divide_cosh) {
  // int_0^tmax e^{-x cosh t} cosh(order t) [/ cosh t] dt, 8-point GL panels of width <= 0.25.
  const double t_max = std::max(2.0, std::log(2.0 * 60.0 / x) + 3.0);
  const auto panels = static_cast<std::size_t>(std::ceil(t_max / 0.25));
  const double width = t_max / static_cast<double>(panels);
  // Gauss-Legendre 8-point nodes/weights by Newton iteration (computed once).
  static const std::pair<std::array<double, 8>, std::array<double, 8>> rule = [] {
    std::array<double, 8> xs{};
    std::array<double, 8> ws{};
    constexpr int n = 8;
    for (int i = 0; i < n; ++i) {
      double z = std::cos(k::k_pi * (i + 0.75) / (n + 0.5));
      double dp = 0.0;
      for (int it = 0; it < 100; ++it) {
        double p0 = 1.0;
        double p1 = z;
        for (int j = 2; j <= n; ++j) {
          const double p2 = ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p0) / j;
          p0 = p1;
          p1 = p2;
        }
        dp = n * (z * p1 - p0) / (z * z - 1.0);
        const double step = p1 / dp;
        z -= step;
        if (std::abs(step) < 1.0e-16) break;
      }
      xs[static_cast<std::size_t>(i)] = z;
      ws[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
    return std::make_pair(xs, ws);
  }();
  double sum = 0.0;
  for (std::size_t p = 0; p < panels; ++p) {
    const double mid = width * (static_cast<double>(p) + 0.5);
    for (std::size_t i = 0; i < 8U; ++i) {
      const double t = mid + 0.5 * width * rule.first[i];
      const double c = std::cosh(t);
      double v = std::exp(-x * c) * std::cosh(order * t);
      if (divide_cosh) v /= c;
      sum += 0.5 * width * rule.second[i] * v;
    }
  }
  return sum;
}

[[nodiscard]] const KernelTable& kernelTable() {
  static const KernelTable table = [] {
    KernelTable t;
    constexpr std::size_t n = 1201;
    t.ln_x_min = std::log(1.0e-8);
    t.step = (std::log(60.0) - t.ln_x_min) / static_cast<double>(n - 1U);
    t.ln_f.resize(n);
    t.ln_g.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      const double x = std::exp(t.ln_x_min + t.step * static_cast<double>(i));
      t.ln_f[i] = std::log(x * kernelIntegral(x, 5.0 / 3.0, true));
      t.ln_g[i] = std::log(x * kernelIntegral(x, 2.0 / 3.0, false));
    }
    return t;
  }();
  return table;
}

[[nodiscard]] double interpolateKernel(const std::vector<double>& lane, double x) {
  const KernelTable& t = kernelTable();
  const double u = (std::log(x) - t.ln_x_min) / t.step;
  const auto last = static_cast<std::ptrdiff_t>(lane.size()) - 1;
  auto j = static_cast<std::ptrdiff_t>(std::floor(u));
  j = std::clamp<std::ptrdiff_t>(j, 1, last - 2);
  const double s = u - static_cast<double>(j);
  const double lm1 = -s * (s - 1.0) * (s - 2.0) / 6.0;
  const double l0 = (s + 1.0) * (s - 1.0) * (s - 2.0) / 2.0;
  const double l1 = -(s + 1.0) * s * (s - 2.0) / 2.0;
  const double l2 = (s + 1.0) * s * (s - 1.0) / 6.0;
  const auto at = [&](std::ptrdiff_t i) { return lane[static_cast<std::size_t>(i)]; };
  return std::exp(lm1 * at(j - 1) + l0 * at(j) + l1 * at(j + 1) + l2 * at(j + 2));
}

using Matrix5 = std::array<double, 25>;

[[nodiscard]] Matrix5 multiply(const Matrix5& a, const Matrix5& b) {
  Matrix5 c{};
  for (std::size_t i = 0; i < 5U; ++i) {
    for (std::size_t kk = 0; kk < 5U; ++kk) {
      const double aik = a[i * 5U + kk];
      for (std::size_t j = 0; j < 5U; ++j) c[i * 5U + j] += aik * b[kk * 5U + j];
    }
  }
  return c;
}

// Solves D X = N (5x5) by Gaussian elimination with partial pivoting.
[[nodiscard]] Matrix5 solve(Matrix5 d, Matrix5 n) {
  for (std::size_t col = 0; col < 5U; ++col) {
    std::size_t pivot = col;
    for (std::size_t r = col + 1U; r < 5U; ++r) {
      if (std::abs(d[r * 5U + col]) > std::abs(d[pivot * 5U + col])) pivot = r;
    }
    if (std::abs(d[pivot * 5U + col]) < 1.0e-300) throw std::runtime_error("Pade denominator is singular");
    if (pivot != col) {
      for (std::size_t j = 0; j < 5U; ++j) {
        std::swap(d[col * 5U + j], d[pivot * 5U + j]);
        std::swap(n[col * 5U + j], n[pivot * 5U + j]);
      }
    }
    for (std::size_t r = 0; r < 5U; ++r) {
      if (r == col) continue;
      const double f = d[r * 5U + col] / d[col * 5U + col];
      if (f == 0.0) continue;
      for (std::size_t j = 0; j < 5U; ++j) {
        d[r * 5U + j] -= f * d[col * 5U + j];
        n[r * 5U + j] -= f * n[col * 5U + j];
      }
    }
  }
  for (std::size_t r = 0; r < 5U; ++r) {
    const double inv = 1.0 / d[r * 5U + r];
    for (std::size_t j = 0; j < 5U; ++j) n[r * 5U + j] *= inv;
  }
  return n;
}

[[nodiscard]] Matrix5 expm(const Matrix5& m) {
  double norm = 0.0;
  for (std::size_t j = 0; j < 5U; ++j) {
    double col = 0.0;
    for (std::size_t i = 0; i < 5U; ++i) col += std::abs(m[i * 5U + j]);
    norm = std::max(norm, col);
  }
  int squarings = 0;
  if (norm > 0.5) squarings = static_cast<int>(std::ceil(std::log2(norm / 0.5)));
  const double scale = std::ldexp(1.0, -squarings);
  Matrix5 x{};
  for (std::size_t i = 0; i < 25U; ++i) x[i] = m[i] * scale;
  // Pade(6,6): c_k = (2q - k)! q! / ((2q)! k! (q - k)!), q = 6 (computed).
  constexpr int q = 6;
  std::array<double, q + 1> c{};
  for (int kk = 0; kk <= q; ++kk) {
    c[static_cast<std::size_t>(kk)] = std::tgamma(2.0 * q - kk + 1.0) * std::tgamma(q + 1.0) /
                                      (std::tgamma(2.0 * q + 1.0) * std::tgamma(kk + 1.0) * std::tgamma(q - kk + 1.0));
  }
  Matrix5 power{};
  for (std::size_t i = 0; i < 5U; ++i) power[i * 5U + i] = 1.0;
  Matrix5 num{};
  Matrix5 den{};
  for (int kk = 0; kk <= q; ++kk) {
    const double sign = (kk % 2 == 0) ? 1.0 : -1.0;
    for (std::size_t i = 0; i < 25U; ++i) {
      num[i] += c[static_cast<std::size_t>(kk)] * power[i];
      den[i] += sign * c[static_cast<std::size_t>(kk)] * power[i];
    }
    power = multiply(power, x);
  }
  Matrix5 result = solve(den, num);
  for (int s = 0; s < squarings; ++s) result = multiply(result, result);
  return result;
}

}  // namespace

double planckBnu(double nu, double t) {
  if (!(nu > 0.0) || !(t > 0.0)) throw std::invalid_argument("Planck function needs nu, T > 0");
  const double x = k::k_planck_si * nu / (k::k_boltzmann_si * t);
  return 2.0 * k::k_planck_si * nu * nu * nu / (k::k_speed_of_light_si * k::k_speed_of_light_si) / std::expm1(x);
}

double rayleighJeansTemperature(double intensity, double nu) {
  return intensity * k::k_speed_of_light_si * k::k_speed_of_light_si / (2.0 * k::k_boltzmann_si * nu * nu);
}

double freeFreeGaunt(double nu, double t, double z) {
  const double arg = z * (nu / 1.0e9) * std::pow(t / 1.0e4, -1.5);
  return std::log(std::exp(5.960 - std::sqrt(3.0) / k::k_pi * std::log(arg)) + std::exp(1.0));
}

double freeFreeEmissivity(double nu, double t, double n_e, double z2_n_i) {
  if (!(t > 0.0) || n_e < 0.0 || z2_n_i < 0.0) throw std::invalid_argument("free-free needs T > 0, n >= 0");
  const double prefactor = std::pow(2.0, 5) * k::k_pi * std::pow(k_q_esu, 6) / (3.0 * k_me_cgs * std::pow(k_c_cgs, 3)) *
                           std::sqrt(2.0 * k::k_pi / (3.0 * k_kb_cgs * k_me_cgs));
  const double eps = prefactor * (n_e * 1.0e-6) * (z2_n_i * 1.0e-6) / std::sqrt(t) *
                     std::exp(-k_h_cgs * nu / (k_kb_cgs * t)) * freeFreeGaunt(nu, t, 1.0);
  return eps / (4.0 * k::k_pi) * k_cgs_emissivity_to_si;
}

double freeFreeAbsorption(double nu, double t, double n_e, double z2_n_i) {
  return freeFreeEmissivity(nu, t, n_e, z2_n_i) / planckBnu(nu, t);
}

double synchrotronKernelF(double x) {
  if (!(x > 0.0)) throw std::invalid_argument("synchrotron kernel needs x > 0");
  if (x < 1.0e-8) return 4.0 * k::k_pi / (std::sqrt(3.0) * std::tgamma(1.0 / 3.0)) * std::cbrt(0.5 * x);
  if (x > 60.0) return std::sqrt(0.5 * k::k_pi * x) * std::exp(-x);
  return interpolateKernel(kernelTable().ln_f, x);
}

double synchrotronKernelG(double x) {
  if (!(x > 0.0)) throw std::invalid_argument("synchrotron kernel needs x > 0");
  if (x < 1.0e-8) return x * 0.5 * std::tgamma(2.0 / 3.0) * std::pow(0.5 * x, -2.0 / 3.0);
  if (x > 60.0) return std::sqrt(0.5 * k::k_pi * x) * std::exp(-x);
  return interpolateKernel(kernelTable().ln_g, x);
}

double synchrotronCriticalFrequency(double gamma, double b_tesla, double sin_alpha) {
  return 3.0 / (4.0 * k::k_pi) * gamma * gamma * k_q_esu * (b_tesla * 1.0e4) * sin_alpha / (k_me_cgs * k_c_cgs);
}

SynchrotronCoefficients synchrotronFromSpectrum(const ElectronSpectrumSample& e, double nu, double b, double sin_alpha) {
  const std::size_t n = e.gamma.size();
  if (n < 3U || e.number_per_m3.size() != n) throw std::invalid_argument("electron spectrum needs >= 3 samples");
  if (!(sin_alpha > 0.0) || !(b > 0.0)) return {};  // no emission perpendicular to nothing: physically zero
  const double b_g = b * 1.0e4;
  const double p_prefactor = std::sqrt(3.0) * std::pow(k_q_esu, 3) * b_g * sin_alpha / (k_me_cgs * k_c_cgs * k_c_cgs);
  double j_f = 0.0;
  double j_g = 0.0;
  double absorption = 0.0;
  for (std::size_t i = 0; i + 1U < n; ++i) {
    // Trapezoid in ln gamma of N P gamma; derivative of N / gamma^2 by central differences in ln gamma.
    const auto integrand = [&](std::size_t idx, double& f_part, double& g_part, double& a_part) {
      const double g = e.gamma[idx];
      const double x = nu / synchrotronCriticalFrequency(g, b, sin_alpha);
      const double n_cgs = e.number_per_m3[idx] * 1.0e-6;
      const double pf = p_prefactor * synchrotronKernelF(x);
      f_part = n_cgs * pf * g;
      g_part = n_cgs * p_prefactor * synchrotronKernelG(x) * g;
      const std::size_t lo = idx == 0 ? 0 : idx - 1U;
      const std::size_t hi = idx + 1U == n ? idx : idx + 1U;
      const double y_lo = e.number_per_m3[lo] * 1.0e-6 / (e.gamma[lo] * e.gamma[lo]);
      const double y_hi = e.number_per_m3[hi] * 1.0e-6 / (e.gamma[hi] * e.gamma[hi]);
      const double dy_dg = (y_hi - y_lo) / (e.gamma[hi] - e.gamma[lo]);
      a_part = pf * g * g * dy_dg * g;  // P gamma^2 d(N/gamma^2)/dgamma, times gamma for d ln gamma
    };
    double f0 = 0.0, g0 = 0.0, a0 = 0.0, f1 = 0.0, g1 = 0.0, a1 = 0.0;
    integrand(i, f0, g0, a0);
    integrand(i + 1U, f1, g1, a1);
    const double dl = std::log(e.gamma[i + 1U] / e.gamma[i]);
    j_f += 0.5 * dl * (f0 + f1);
    j_g += 0.5 * dl * (g0 + g1);
    absorption += 0.5 * dl * (a0 + a1);
  }
  SynchrotronCoefficients out;
  out.j_i = j_f / (4.0 * k::k_pi) * k_cgs_emissivity_to_si;
  out.linear_fraction = j_f > 0.0 ? j_g / j_f : 0.0;
  out.alpha_i = -absorption / (8.0 * k::k_pi * k_me_cgs * nu * nu) * k_cgs_absorption_to_si;
  return out;
}

double synchrotronPowerLawEmissivity(double c_m3, double p, double nu, double b, double sin_alpha) {
  const double c_cgs = c_m3 * 1.0e-6;
  const double bs = b * 1.0e4 * sin_alpha;
  const double omega = 2.0 * k::k_pi * nu;
  const double p_omega = std::sqrt(3.0) * std::pow(k_q_esu, 3) * c_cgs * bs / (2.0 * k::k_pi * k_me_cgs * k_c_cgs * k_c_cgs * (p + 1.0)) *
                         std::tgamma(p / 4.0 + 19.0 / 12.0) * std::tgamma(p / 4.0 - 1.0 / 12.0) *
                         std::pow(k_me_cgs * k_c_cgs * omega / (3.0 * k_q_esu * bs), -(p - 1.0) / 2.0);
  return 2.0 * k::k_pi * p_omega / (4.0 * k::k_pi) * k_cgs_emissivity_to_si;
}

double synchrotronPowerLawAbsorption(double c_m3, double p, double nu, double b, double sin_alpha) {
  const double c_cgs = c_m3 * 1.0e-6;
  const double bs = b * 1.0e4 * sin_alpha;
  const double alpha = std::sqrt(3.0) * std::pow(k_q_esu, 3) / (8.0 * k::k_pi * k_me_cgs) *
                       std::pow(3.0 * k_q_esu / (2.0 * k::k_pi * std::pow(k_me_cgs, 3) * std::pow(k_c_cgs, 5)), p / 2.0) *
                       c_cgs * std::pow(bs, (p + 2.0) / 2.0) * std::tgamma((3.0 * p + 2.0) / 12.0) *
                       std::tgamma((3.0 * p + 22.0) / 12.0) * std::pow(nu, -(p + 4.0) / 2.0);
  return alpha * k_cgs_absorption_to_si;
}

double synchrotronPowerLawPolarization(double p) { return (p + 1.0) / (p + 7.0 / 3.0); }

double rotationMeasureCoefficient() {
  const double e = k::k_elementary_charge_si;
  return e * e * e / (8.0 * k::k_pi * k::k_pi * k::k_vacuum_permittivity_si * k::k_electron_mass_si *
                      k::k_electron_mass_si * std::pow(k::k_speed_of_light_si, 3));
}

double faradayRotationRate(double wavelength, double n_e, double b_par) {
  return 2.0 * rotationMeasureCoefficient() * wavelength * wavelength * n_e * b_par;
}

StokesMatrix propagationMatrix(double ai, double aq, double au, double av, double rq, double ru, double rv) {
  return {ai, aq, au, av, aq, ai, rv, ru, au, -rv, ai, rq, av, -ru, -rq, ai};
}

Stokes stokesStep(const Stokes& s_in, const Stokes& j, const StokesMatrix& kmat, double ds) {
  if (!(ds >= 0.0)) throw std::invalid_argument("Stokes step needs ds >= 0");
  Matrix5 m{};
  for (std::size_t r = 0; r < 4U; ++r) {
    for (std::size_t c = 0; c < 4U; ++c) m[r * 5U + c] = -kmat[r * 4U + c] * ds;
    m[r * 5U + 4U] = j[r] * ds;
  }
  const Matrix5 e = expm(m);
  Stokes out{};
  for (std::size_t r = 0; r < 4U; ++r) {
    double v = e[r * 5U + 4U];
    for (std::size_t c = 0; c < 4U; ++c) v += e[r * 5U + c] * s_in[c];
    out[r] = v;
  }
  return out;
}

DispersionMeasures integrateDispersion(const std::vector<DispersionSegment>& segments) {
  DispersionMeasures out;
  const double c = rotationMeasureCoefficient();
  for (const DispersionSegment& s : segments) {
    const double zp1 = 1.0 + s.redshift;
    out.dm_m2 += s.electron_density_m3 * s.path_proper_m / zp1;
    out.rm_rad_m2 += c * s.electron_density_m3 * s.b_parallel_tesla * s.path_proper_m / (zp1 * zp1);
  }
  return out;
}

}  // namespace cosmosim::observables
