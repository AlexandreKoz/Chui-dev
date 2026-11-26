#include "cosmosim/physics/black_hole/bh_dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::physics {
namespace {

constexpr double k_sqrt_pi = 1.7724538509055160273;

[[nodiscard]] double norm(const std::array<double, 3>& a) noexcept {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

[[nodiscard]] std::array<double, 3> totalAcceleration(const std::array<double, 3>& v_bh,
                                                      const BhFrictionBackground& dm,
                                                      const BhFrictionBackground& st,
                                                      const BhFrictionInputs& in) noexcept {
  std::array<double, 3> a{};
  for (const BhFrictionBackground* bg : {&dm, &st}) {
    if (bg->samples == 0U || !(bg->density_proper > 0.0)) continue;
    const std::array<double, 3> v_rel{v_bh[0] - bg->bulk_velocity[0], v_bh[1] - bg->bulk_velocity[1],
                                      v_bh[2] - bg->bulk_velocity[2]};
    const std::array<double, 3> ak = chandrasekharAcceleration(v_rel, *bg, in);
    a[0] += ak[0];
    a[1] += ak[1];
    a[2] += ak[2];
  }
  return a;
}

}  // namespace

double chandrasekharBracket(double x) noexcept {
  if (!(x > 0.0)) return 0.0;
  if (x < 0.05) {
    const double x2 = x * x;
    return 4.0 * x2 * x / (3.0 * k_sqrt_pi) * (1.0 - 0.6 * x2 + (3.0 / 14.0) * x2 * x2 - x2 * x2 * x2 / 18.0);
  }
  return std::erf(x) - 2.0 * x * std::exp(-x * x) / k_sqrt_pi;
}

double tremmelCoulombLogarithm(double softening, double g, double m, double speed) noexcept {
  if (!(softening > 0.0) || !(speed > 0.0) || !(m > 0.0) || !(g > 0.0)) return 0.0;
  const double b_min = g * m / (speed * speed);
  if (!(b_min > 0.0) || softening <= b_min) return 0.0;
  return std::log(softening / b_min);
}

std::array<double, 3> chandrasekharAcceleration(const std::array<double, 3>& v_rel, const BhFrictionBackground& bg,
                                                const BhFrictionInputs& in) noexcept {
  const double v = norm(v_rel);
  if (!(v > 0.0) || !(bg.density_proper > 0.0) || bg.samples == 0U) return {0.0, 0.0, 0.0};
  const double ln_lambda = tremmelCoulombLogarithm(in.softening_proper, in.newton_g_code, in.bh_mass_code, v);
  if (!(ln_lambda > 0.0)) return {0.0, 0.0, 0.0};
  const double x = bg.sigma_1d > 0.0 ? v / (std::sqrt(2.0) * bg.sigma_1d) : std::numeric_limits<double>::infinity();
  const double bracket = std::isfinite(x) ? chandrasekharBracket(x) : 1.0;  // cold background limit
  const double coef = -4.0 * core::constants::k_pi * in.newton_g_code * in.newton_g_code * in.bh_mass_code *
      bg.density_proper * ln_lambda * bracket / (v * v * v);
  return {coef * v_rel[0], coef * v_rel[1], coef * v_rel[2]};
}

BhFrictionResult integrateDynamicalFriction(const std::array<double, 3>& v0, const BhFrictionBackground& dm,
                                            const BhFrictionBackground& st, const BhFrictionInputs& in,
                                            double dt) {
  if (!(dt >= 0.0) || !std::isfinite(dt) || !(in.timestep_eta > 0.0)) {
    throw std::invalid_argument("integrateDynamicalFriction: invalid interval or eta");
  }
  BhFrictionResult out;
  std::array<double, 3> v = v0;
  const std::array<double, 3> a0 = totalAcceleration(v, dm, st, in);
  const double a0n = norm(a0);
  if (!(a0n > 0.0) || dt == 0.0) return out;
  // Reference speed for the limiter: relative to the mass-weighted background.
  auto relative_speed = [&](const std::array<double, 3>& vv) {
    double best = 0.0;
    for (const BhFrictionBackground* bg : {&dm, &st}) {
      if (bg->samples == 0U) continue;
      best = std::max(best, std::hypot(vv[0] - bg->bulk_velocity[0], vv[1] - bg->bulk_velocity[1],
                                       vv[2] - bg->bulk_velocity[2]));
    }
    return best;
  };
  out.timestep_limit_code = in.timestep_eta * relative_speed(v) / a0n;
  double t = 0.0;
  while (t < dt) {
    if (out.substeps >= in.max_substeps) {
      throw std::runtime_error("integrateDynamicalFriction: substep budget exhausted (stiff friction; reduce dt)");
    }
    const std::array<double, 3> a = totalAcceleration(v, dm, st, in);
    const double an = norm(a);
    if (!(an > 0.0)) break;
    const double h = std::min(dt - t, std::max(in.timestep_eta * relative_speed(v) / an, 1.0e-12 * dt));
    std::array<double, 3> vm{v[0] + 0.5 * h * a[0], v[1] + 0.5 * h * a[1], v[2] + 0.5 * h * a[2]};
    const std::array<double, 3> am = totalAcceleration(vm, dm, st, in);
    v = {v[0] + h * am[0], v[1] + h * am[1], v[2] + h * am[2]};
    t += h;
    ++out.substeps;
  }
  out.velocity_change = {v[0] - v0[0], v[1] - v0[1], v[2] - v0[2]};
  out.applied = true;
  return out;
}

}  // namespace cosmosim::physics
