#include "cosmosim/radiation/radiation_matter_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cosmosim/radiation/radiation_closure.hpp"

namespace cosmosim::radiation {
namespace {

using Vec3 = std::array<double, 3>;
[[nodiscard]] double dot(const Vec3& a, const Vec3& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

struct GreyEval {
  double r_e = 0.0;
  double r_t = 0.0;
  Vec3 flux{};
  Vec3 momentum{};
  double total_energy = 0.0;
  double specific_energy = 0.0;
  bool valid = true;
};

}  // namespace

RadiationLocalResult solveGreyThermalImplicit(RadiationMatterCell& cell, double& e_rad, Vec3& g_rad, double alpha,
                                              const RadiationCouplingConstants& k) {
  RadiationLocalResult out;
  const double c = k.c_code;
  const double ch = k.c_hat_code;
  const double rho = cell.density;
  if (!(rho > 0.0) || !(cell.temperature_per_specific_energy > 0.0) || !(alpha >= 0.0) || !(e_rad >= 0.0)) {
    throw std::invalid_argument("solveGreyThermalImplicit: invalid cell state");
  }
  const double kp = rho * k.kappa_absorption_code;
  const double kf = rho * (k.kappa_absorption_code + k.kappa_scattering_code);
  const double a_rad = k.radiation_constant_code;
  const Vec3 v = cell.velocity;
  const double e0 = e_rad;
  const Vec3 f0{g_rad[0] * c / ch, g_rad[1] * c / ch, g_rad[2] * c / ch};
  const Vec3 p0{rho * v[0], rho * v[1], rho * v[2]};
  const double etot0 = rho * (cell.specific_internal_energy + 0.5 * dot(v, v));
  const double tpu = cell.temperature_per_specific_energy;
  // Frozen Eddington tensor D = P/E of the right-hand-side state.
  std::array<double, 9> d{};
  {
    const M1State st{e0, f0};
    const std::array<double, 9> p = m1PressureTensor(st, c);
    for (int i = 0; i < 9; ++i) d[i] = e0 > 0.0 ? p[i] / e0 : (i % 4 == 0 ? 1.0 / 3.0 : 0.0);
  }
  auto evaluate = [&](double t, double e) {
    GreyEval ev;
    const double at4 = a_rad * t * t * t * t;
    const Vec3 dv{d[0] * v[0] + d[1] * v[1] + d[2] * v[2], d[3] * v[0] + d[4] * v[1] + d[5] * v[2],
                  d[6] * v[0] + d[7] * v[1] + d[8] * v[2]};
    const double denom = 1.0 + alpha * ch * kf;
    for (int i = 0; i < 3; ++i) {
      ev.flux[i] = (f0[i] + alpha * ch * (kp * at4 * v[i] + kf * (v[i] * e + e * dv[i]))) / denom;
    }
    const double g0 = kp * (e - at4) + (kf - 2.0 * kp) * dot(v, ev.flux) / (c * c);
    Vec3 gm{};
    for (int i = 0; i < 3; ++i) gm[i] = kf * ev.flux[i] / c - kp * at4 * v[i] / c - kf * (v[i] * e + e * dv[i]) / c;
    for (int i = 0; i < 3; ++i) ev.momentum[i] = p0[i] + alpha * gm[i];
    ev.total_energy = etot0 + alpha * c * g0;
    ev.specific_energy = (ev.total_energy - dot(ev.momentum, ev.momentum) / (2.0 * rho)) / rho;
    ev.r_e = e - e0 + alpha * ch * g0;
    ev.r_t = t - ev.specific_energy * tpu;
    ev.valid = std::isfinite(ev.r_e) && std::isfinite(ev.r_t) && std::isfinite(ev.total_energy);
    return ev;
  };
  double t = cell.specific_internal_energy * tpu;
  double e = e0;
  const double t_scale = std::max(t, 1.0e-300);
  const double e_scale = std::max({e0, a_rad * t * t * t * t, 1.0e-300});
  auto norm = [&](const GreyEval& ev) { return std::max(std::abs(ev.r_e) / e_scale, std::abs(ev.r_t) / t_scale); };
  GreyEval cur = evaluate(t, e);
  for (out.iterations = 0; out.iterations < k.max_iterations; ++out.iterations) {
    out.residual = norm(cur);
    if (cur.valid && out.residual <= k.relative_tolerance) break;
    // Forward-difference Jacobian (deterministic relative steps).
    const double ht = 1.0e-7 * std::max(t, t_scale * 1.0e-6);
    const double he = 1.0e-7 * std::max(e, e_scale * 1.0e-6);
    const GreyEval et = evaluate(t + ht, e);
    const GreyEval ee = evaluate(t, e + he);
    const double j11 = (et.r_e - cur.r_e) / ht, j12 = (ee.r_e - cur.r_e) / he;
    const double j21 = (et.r_t - cur.r_t) / ht, j22 = (ee.r_t - cur.r_t) / he;
    const double det = j11 * j22 - j12 * j21;
    if (!std::isfinite(det) || det == 0.0) {
      out.converged = false;
      return out;
    }
    const double dt_step = (-cur.r_e * j22 + cur.r_t * j12) / det;
    const double de_step = (-cur.r_t * j11 + cur.r_e * j21) / det;
    double lambda = 1.0;
    bool accepted = false;
    for (int ls = 0; ls < 40; ++ls) {
      const double t_try = t + lambda * dt_step;
      const double e_try = e + lambda * de_step;
      if (t_try > 0.0 && e_try >= 0.0) {
        const GreyEval trial = evaluate(t_try, e_try);
        if (trial.valid && (norm(trial) < out.residual || norm(trial) <= k.relative_tolerance)) {
          t = t_try;
          e = e_try;
          cur = trial;
          accepted = true;
          break;
        }
      }
      lambda *= 0.5;
    }
    if (!accepted) {
      out.converged = false;
      return out;
    }
  }
  out.residual = norm(cur);
  if (!(out.residual <= k.relative_tolerance) || !(cur.specific_energy > 0.0)) {
    out.converged = false;
    return out;
  }
  out.converged = true;
  out.gas_energy_change = cur.total_energy - etot0;
  for (int i = 0; i < 3; ++i) out.gas_momentum_change[i] = cur.momentum[i] - p0[i];
  e_rad = e;
  for (int i = 0; i < 3; ++i) g_rad[i] = cur.flux[i] * ch / c;
  cell.specific_internal_energy = cur.specific_energy;
  for (int i = 0; i < 3; ++i) cell.velocity[i] = cur.momentum[i] / rho;
  return out;
}

RadiationLocalResult solvePhotoionizationImplicit(const RadiationMatterCell& cell, const RadiationSpectralDefinition& sp,
                                                  std::span<double> e_rad, std::span<std::array<double, 3>> g_rad,
                                                  double alpha, const RadiationCouplingConstants& k) {
  RadiationLocalResult out;
  if (e_rad.size() != sp.groupCount() || g_rad.size() != sp.groupCount()) {
    throw std::invalid_argument("solvePhotoionizationImplicit: group count mismatch");
  }
  for (std::size_t g = 0; g < sp.groupCount(); ++g) {
    // Opacity per proper code length.
    double chi = 0.0;
    std::array<double, k_radiation_absorber_count> chi_s{};
    for (std::size_t s = 0; s < k_radiation_absorber_count; ++s) {
      chi_s[s] = cell.absorber_number_density_cgs[s] * sp.groups[g].sigma_energy_cm2[s] * k.length_cm_per_code;
      chi += chi_s[s];
    }
    const double damp = 1.0 / (1.0 + alpha * k.c_hat_code * chi);
    e_rad[g] *= damp;
    for (int i = 0; i < 3; ++i) g_rad[g][i] *= damp;
    // Matter-side: energy alpha c chi E (ledgered for chemistry), momentum
    // alpha chi F/c with F = G c/c_hat.
    for (std::size_t s = 0; s < k_radiation_absorber_count; ++s) {
      out.absorbed_energy_density[s] += alpha * k.c_code * chi_s[s] * e_rad[g];
    }
    for (int i = 0; i < 3; ++i) out.gas_momentum_change[i] += alpha * chi * g_rad[g][i] / k.c_hat_code;
  }
  out.iterations = 1U;
  return out;
}

physics::TransportedPhotoRates transportedPhotoRates(const RadiationSpectralDefinition& sp,
                                                     std::span<const double> e_cgs) {
  physics::TransportedPhotoRates r;
  r.present = true;
  constexpr double c_cgs = 2.99792458e10;
  for (std::size_t g = 0; g < sp.groupCount() && g < e_cgs.size(); ++g) {
    const RadiationGroup& grp = sp.groups[g];
    if (!(grp.mean_photon_energy_ev > 0.0) || !(e_cgs[g] > 0.0)) continue;
    const double photons = e_cgs[g] / (grp.mean_photon_energy_ev * k_ev_erg);
    const double gamma[3] = {c_cgs * grp.sigma_number_cm2[0] * photons, c_cgs * grp.sigma_number_cm2[1] * photons,
                             c_cgs * grp.sigma_number_cm2[2] * photons};
    const double excess[3] = {
        (grp.absorbed_energy_per_ionization_ev[0] - k_ionization_potential_hi_ev) * k_ev_erg,
        (grp.absorbed_energy_per_ionization_ev[1] - k_ionization_potential_hei_ev) * k_ev_erg,
        (grp.absorbed_energy_per_ionization_ev[2] - k_ionization_potential_heii_ev) * k_ev_erg};
    r.gamma_hi_s += gamma[0];
    r.gamma_hei_s += gamma[1];
    r.gamma_heii_s += gamma[2];
    r.heat_hi_erg_s += gamma[0] * std::max(excess[0], 0.0);
    r.heat_hei_erg_s += gamma[1] * std::max(excess[1], 0.0);
    r.heat_heii_erg_s += gamma[2] * std::max(excess[2], 0.0);
  }
  return r;
}

physics::TransportedPhotoRates RadiationPhotoRateCache::ratesForRow(const core::SimulationState&,
                                                                    std::uint32_t row) const {
  if (row >= m_rates.size()) return physics::TransportedPhotoRates{};
  return m_rates[row];
}

}  // namespace cosmosim::radiation
