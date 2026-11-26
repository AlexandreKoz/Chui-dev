#include "cosmosim/radiation/radiation_closure.hpp"

#include <algorithm>
#include <cmath>

namespace cosmosim::radiation {
namespace {

constexpr double k_isotropic_f = 1.0e-8;

[[nodiscard]] double norm3(const std::array<double, 3>& v) noexcept {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Real roots of x^3 + a x^2 + b x + c = 0 (all real for hyperbolic M1;
// round-off with a slightly negative discriminant takes the real parts).
void cubicRealRoots(double a, double b, double c, double& r_min, double& r_max) noexcept {
  const double a3 = a / 3.0;
  const double p = b - a * a / 3.0;
  const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
  if (std::abs(p) < 1.0e-300) {
    const double root = std::cbrt(-q) - a3;
    r_min = root;
    r_max = root;
    return;
  }
  if (p < 0.0) {
    const double m = 2.0 * std::sqrt(-p / 3.0);
    double arg = 3.0 * q / (p * m);
    arg = std::clamp(arg, -1.0, 1.0);
    const double theta = std::acos(arg) / 3.0;
    const double k_two_pi_3 = 2.0943951023931954923;
    const double x0 = m * std::cos(theta) - a3;
    const double x1 = m * std::cos(theta - k_two_pi_3) - a3;
    const double x2 = m * std::cos(theta - 2.0 * k_two_pi_3) - a3;
    r_min = std::min({x0, x1, x2});
    r_max = std::max({x0, x1, x2});
    return;
  }
  // p > 0: one real root (should not occur for admissible M1 states).
  const double disc = q * q / 4.0 + p * p * p / 27.0;
  const double s = std::sqrt(std::max(disc, 0.0));
  const double root = std::cbrt(-q / 2.0 + s) + std::cbrt(-q / 2.0 - s) - a3;
  const double real_pair = -0.5 * (root + a3) - a3;  // real part of the complex pair
  r_min = std::min(root, real_pair);
  r_max = std::max(root, real_pair);
}

}  // namespace

double m1EddingtonFactor(double f) noexcept {
  f = std::clamp(f, 0.0, 1.0);
  return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * std::sqrt(4.0 - 3.0 * f * f));
}

double m1EddingtonFactorDerivative(double f) noexcept {
  f = std::clamp(f, 0.0, 1.0);
  const double s = std::sqrt(4.0 - 3.0 * f * f);
  const double d = 5.0 + 2.0 * s;
  return (8.0 * f * d + 6.0 * f * (3.0 + 4.0 * f * f) / s) / (d * d);
}

std::array<double, 9> m1PressureTensor(const M1State& st, double c) noexcept {
  std::array<double, 9> p{};
  if (!(st.energy > 0.0)) return p;
  const double fn = norm3(st.flux);
  const double f = fn / (c * st.energy);
  if (f < k_isotropic_f) {
    p[0] = p[4] = p[8] = st.energy / 3.0;
    return p;
  }
  const double chi = m1EddingtonFactor(f);
  const double a = 0.5 * (1.0 - chi);
  const double b = 0.5 * (3.0 * chi - 1.0);
  const std::array<double, 3> n{st.flux[0] / fn, st.flux[1] / fn, st.flux[2] / fn};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) p[3 * i + j] = st.energy * ((i == j ? a : 0.0) + b * n[i] * n[j]);
  }
  return p;
}

bool m1Admissible(const M1State& st, double c, double tol) noexcept {
  if (!std::isfinite(st.energy) || !std::isfinite(st.flux[0]) || !std::isfinite(st.flux[1]) ||
      !std::isfinite(st.flux[2]) || st.energy < 0.0) {
    return false;
  }
  return norm3(st.flux) <= c * st.energy * (1.0 + tol) + 1.0e-300;
}

std::array<double, 4> m1PhysicalFlux(const M1State& st, int axis, double c) noexcept {
  const std::array<double, 9> p = m1PressureTensor(st, c);
  const double c2 = c * c;
  return {st.flux[axis], c2 * p[3 * axis + 0], c2 * p[3 * axis + 1], c2 * p[3 * axis + 2]};
}

M1WaveSpeeds m1WaveSpeeds(const M1State& st, int axis, double c) noexcept {
  M1WaveSpeeds out{-c, c};
  if (!(st.energy > 0.0)) return out;
  const double r = norm3(st.flux);
  const double f = std::min(r / (c * st.energy), 1.0);
  if (f < k_isotropic_f) {
    const double s = c / std::sqrt(3.0);
    return {-s, s};
  }
  const double fn = st.flux[axis];
  const double ft = std::sqrt(std::max(r * r - fn * fn, 0.0));
  const double nx = fn / r;
  const double ny = ft / r;
  const double chi = m1EddingtonFactor(f);
  const double dchi = m1EddingtonFactorDerivative(f);
  const double a = 0.5 * (1.0 - chi);
  const double b = 0.5 * (3.0 * chi - 1.0);
  const double da = -0.5 * dchi;
  const double db = 1.5 * dchi;
  const double b_over_cf = b / (c * f);  // E b / |F|
  // Dimensionless Jacobian J = A / c in the (normal, in-plane) frame, using
  // variables (E, F/c): rows: d(F_x/c)/dU, d(c P_xx)/dU, d(c P_xy)/dU.
  // dPxx/dE = (a + b nx^2) - f (da + db nx^2)
  const double pxx_e = (a + b * nx * nx) - f * (da + db * nx * nx);
  const double pxx_fx = (da + db * nx * nx) * nx + 2.0 * b_over_cf * c * nx * ny * ny;   // x c (per F/c)
  const double pxx_fy = (da + db * nx * nx) * ny - 2.0 * b_over_cf * c * nx * nx * ny;
  const double pxy_e = b * nx * ny - f * db * nx * ny;
  const double pxy_fx = db * nx * nx * ny + b_over_cf * c * ny * (ny * ny - nx * nx);
  const double pxy_fy = db * nx * ny * ny + b_over_cf * c * nx * (nx * nx - ny * ny);
  // Matrix M (eigenvalues lambda/c):
  //  [ 0      1      0     ]
  //  [ pxx_e  pxx_fx pxx_fy]
  //  [ pxy_e  pxy_fx pxy_fy]
  const double m10 = pxx_e, m11 = pxx_fx, m12 = pxx_fy;
  const double m20 = pxy_e, m21 = pxy_fx, m22 = pxy_fy;
  const double trace = m11 + m22;
  const double minors = (0.0 * m11 - 1.0 * m10) + (0.0 * m22 - 0.0 * m20) + (m11 * m22 - m12 * m21);
  const double det = 0.0 * (m11 * m22 - m12 * m21) - 1.0 * (m10 * m22 - m12 * m20) + 0.0;
  double lo = 0.0;
  double hi = 0.0;
  cubicRealRoots(-trace, minors, -det, lo, hi);
  const double lambda_z = (3.0 * chi - 1.0) / (2.0 * f) * nx;
  lo = std::min(lo, lambda_z);
  hi = std::max(hi, lambda_z);
  out.min = std::clamp(lo, -1.0, 1.0) * c;
  out.max = std::clamp(hi, -1.0, 1.0) * c;
  if (!std::isfinite(out.min) || !std::isfinite(out.max) || out.min > out.max) return {-c, c};
  return out;
}

double m1AdmissibleFraction(const M1State& base, const M1State& delta, double c) noexcept {
  auto at = [&](double t) {
    M1State s;
    s.energy = base.energy + t * delta.energy;
    for (int k = 0; k < 3; ++k) s.flux[k] = base.flux[k] + t * delta.flux[k];
    return s;
  };
  if (m1Admissible(at(1.0), c, 0.0)) return 1.0;  // convex cone: whole segment admissible
  if (!m1Admissible(base, c, 0.0)) return 0.0;
  double lo = 0.0;
  double hi = 1.0;
  for (int it = 0; it < 60; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (m1Admissible(at(mid), c, 0.0)) lo = mid; else hi = mid;
  }
  return lo;
}

}  // namespace cosmosim::radiation
