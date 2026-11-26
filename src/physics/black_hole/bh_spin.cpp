#include "cosmosim/physics/black_hole/bh_spin.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::physics {
namespace {

[[nodiscard]] double norm(const std::array<double, 3>& a) noexcept {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}
void requireSpin(double a) {
  if (!std::isfinite(a) || std::abs(a) > 1.0) {
    throw std::domain_error("Kerr spin parameter must satisfy |a*| <= 1");
  }
}

}  // namespace

double kerrIscoRadius(double a) {
  requireSpin(a);
  const double z1 = 1.0 + std::cbrt(1.0 - a * a) * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
  const double z2 = std::sqrt(3.0 * a * a + z1 * z1);
  const double sgn = (a > 0.0) ? 1.0 : (a < 0.0 ? -1.0 : 0.0);
  // (3 - Z1) -> 0 at a -> 1: clamp the product's round-off at zero.
  const double prod = std::max((3.0 - z1) * (3.0 + z1 + 2.0 * z2), 0.0);
  return 3.0 + z2 - sgn * std::sqrt(prod);
}

double kerrIscoEnergy(double a) {
  const double r = kerrIscoRadius(a);
  return std::sqrt(1.0 - 2.0 / (3.0 * r));
}

double kerrRadiativeEfficiency(double a) { return 1.0 - kerrIscoEnergy(a); }

double kerrIscoSpecificAngularMomentum(double a) {
  const double r = kerrIscoRadius(a);
  const double sr = std::sqrt(r);
  const double denom_inner = std::max(r * sr - 3.0 * sr + 2.0 * a, 0.0);
  if (!(denom_inner > 0.0)) {
    // Near-extremal prograde conditioning: r -> 1, numerator and denominator
    // both -> 0; use the analytic extremal limit l -> 2/sqrt(3).
    return 2.0 / std::sqrt(3.0);
  }
  return (r * r - 2.0 * a * sr + a * a) / (std::pow(r, 0.75) * std::sqrt(denom_inner));
}

double kerrIrreducibleMass(double mass, double a) noexcept {
  const double aa = std::min(std::abs(a), 1.0);
  return mass * std::sqrt(0.5 * (1.0 + std::sqrt(1.0 - aa * aa)));
}

double kerrHorizonAngularVelocity(double a, double r_g, double c) noexcept {
  const double aa = std::min(std::abs(a), 1.0);
  const double r_h = r_g * (1.0 + std::sqrt(1.0 - aa * aa));
  return r_h > 0.0 ? aa * c / (2.0 * r_h) : 0.0;
}

KerrAccretionResult accreteKerrThinDisk(const KerrAccretionInputs& in) {
  if (!(in.mass_code > 0.0) || !(in.inner_mass_code >= 0.0) || !(in.newton_g_code > 0.0) ||
      !(in.speed_of_light_code > 0.0) || !(in.max_relative_mass_step > 0.0)) {
    throw std::invalid_argument("accreteKerrThinDisk: invalid inputs");
  }
  const double axis_norm = norm(in.disk_axis);
  if (!(axis_norm > 0.0)) throw std::invalid_argument("accreteKerrThinDisk: undefined disk axis");
  const std::array<double, 3> n{in.disk_axis[0] / axis_norm, in.disk_axis[1] / axis_norm,
                                in.disk_axis[2] / axis_norm};
  KerrAccretionResult out;
  out.mass_code = in.mass_code;
  out.spin_angular_momentum_code = in.spin_angular_momentum_code;
  const double g = in.newton_g_code;
  const double c = in.speed_of_light_code;
  const double j0 = norm(in.spin_angular_momentum_code);
  if (j0 > 0.0) {
    const double cosine = (in.spin_angular_momentum_code[0] * n[0] + in.spin_angular_momentum_code[1] * n[1] +
                           in.spin_angular_momentum_code[2] * n[2]) / j0;
    if (std::abs(cosine) < 1.0 - 1.0e-9) {
      out.misaligned_rejected = true;
      return out;
    }
  }
  double remaining = in.inner_mass_code;
  while (remaining > 0.0) {
    const double m = out.mass_code;
    const double dm0 = std::min(remaining, in.max_relative_mass_step * m);
    const double j_par = out.spin_angular_momentum_code[0] * n[0] + out.spin_angular_momentum_code[1] * n[1] +
        out.spin_angular_momentum_code[2] * n[2];
    const double a = std::clamp(c * j_par / (g * m * m), -1.0, 1.0);
    const double e = kerrIscoEnergy(a);
    const double l = kerrIscoSpecificAngularMomentum(a);
    const double m_new = m + e * dm0;
    const double dj = l * g * m / c * dm0;
    for (std::size_t k = 0; k < 3; ++k) out.spin_angular_momentum_code[k] += dj * n[k];
    out.radiated_energy_code += (1.0 - e) * dm0 * c * c;
    out.mass_code = m_new;
    const double j_max = g * m_new * m_new / c;
    const double j_now = norm(out.spin_angular_momentum_code);
    if (j_now > j_max) {
      const double excess = j_now - j_max;
      for (std::size_t k = 0; k < 3; ++k) {
        const double unit = out.spin_angular_momentum_code[k] / j_now;
        out.saturation_angular_momentum_code[k] += excess * unit;
        out.spin_angular_momentum_code[k] = j_max * unit;
      }
    }
    remaining -= dm0;
    ++out.substeps;
  }
  return out;
}

BlandfordZnajekResult extractBlandfordZnajek(const BlandfordZnajekInputs& in) {
  if (!(in.mass_code > 0.0) || !(in.newton_g_code > 0.0) || !(in.speed_of_light_code > 0.0) || !(in.dt_code >= 0.0) ||
      !(in.field_angular_velocity_ratio > 0.0) || !(in.field_angular_velocity_ratio < 1.0) || !(in.kappa > 0.0)) {
    throw std::invalid_argument("extractBlandfordZnajek: invalid inputs (requires 0 < Omega_F/Omega_H < 1)");
  }
  BlandfordZnajekResult out;
  out.mass_code = in.mass_code;
  out.spin_angular_momentum_code = in.spin_angular_momentum_code;
  const double g = in.newton_g_code;
  const double c = in.speed_of_light_code;
  auto spin_of = [&](double m, const std::array<double, 3>& j) { return std::min(c * norm(j) / (g * m * m), 1.0); };
  out.irreducible_mass_before_code = kerrIrreducibleMass(in.mass_code, spin_of(in.mass_code, in.spin_angular_momentum_code));
  double t = 0.0;
  constexpr int k_max_substeps = 100000;
  for (int step = 0; step < k_max_substeps && t < in.dt_code; ++step) {
    const double m = out.mass_code;
    const double j = norm(out.spin_angular_momentum_code);
    const double a = spin_of(m, out.spin_angular_momentum_code);
    if (!(a > 1.0e-12) || !(j > 0.0)) break;  // zero-spin limit: no extraction
    const double r_g = g * m / (c * c);
    const double omega_h = kerrHorizonAngularVelocity(a, r_g, c);
    const double x = omega_h * r_g / c;
    const double f = 1.0 + 1.38 * x * x - 9.2 * x * x * x * x;
    const double power = in.kappa / (4.0 * core::constants::k_pi * c) * in.magnetic_flux_gaussian_code *
        in.magnetic_flux_gaussian_code * omega_h * omega_h * std::max(f, 0.0);
    if (!(power > 0.0)) break;
    const double omega_f = in.field_angular_velocity_ratio * omega_h;
    // Substep limited to 1% of the current angular momentum.
    const double h = std::min(in.dt_code - t, 0.01 * j * omega_f / power);
    const double de = power * h;
    const double dj = de / omega_f;
    for (std::size_t k = 0; k < 3; ++k) {
      const double unit = out.spin_angular_momentum_code[k] / j;
      out.spin_angular_momentum_code[k] -= dj * unit;
      out.extracted_angular_momentum_code[k] += dj * unit;
    }
    out.mass_code -= de / (c * c);
    out.extracted_energy_code += de;
    t += h;
  }
  out.irreducible_mass_after_code =
      kerrIrreducibleMass(out.mass_code, spin_of(out.mass_code, out.spin_angular_momentum_code));
  out.area_theorem_satisfied =
      out.irreducible_mass_after_code >= out.irreducible_mass_before_code * (1.0 - 1.0e-12);
  return out;
}

}  // namespace cosmosim::physics
