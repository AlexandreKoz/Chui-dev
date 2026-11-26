#include "cosmosim/physics/black_hole/bh_supply.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cosmosim::physics {
namespace {

[[nodiscard]] std::array<double, 3> cross(const std::array<double, 3>& a, const std::array<double, 3>& b) noexcept {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
[[nodiscard]] double norm(const std::array<double, 3>& a) noexcept {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}
void add3(std::array<double, 3>& a, const std::array<double, 3>& b) noexcept {
  a[0] += b[0];
  a[1] += b[1];
  a[2] += b[2];
}

}  // namespace

double bhKernelShape(double q) noexcept {
  if (!(q >= 0.0) || q >= 1.0) return 0.0;
  const double one_minus = 1.0 - q;
  const double t = one_minus * one_minus;
  return t * t * (1.0 + 4.0 * q);
}

double bhCellKernelWeight(const std::array<double, 3>& offset, const std::array<double, 3>& half_width,
                          double h) noexcept {
  if (!(h > 0.0)) return 0.0;
  constexpr double k_node = 0.57735026918962576451;  // 1/sqrt(3)
  double sum = 0.0;
  for (int ix = -1; ix <= 1; ix += 2) {
    for (int iy = -1; iy <= 1; iy += 2) {
      for (int iz = -1; iz <= 1; iz += 2) {
        const double x = offset[0] + ix * k_node * half_width[0];
        const double y = offset[1] + iy * k_node * half_width[1];
        const double z = offset[2] + iz * k_node * half_width[2];
        sum += bhKernelShape(std::sqrt(x * x + y * y + z * z) / h);
      }
    }
  }
  return sum / 8.0;
}

void BhPhaseSums::merge(const BhPhaseSums& rhs) {
  weight_mass += rhs.weight_mass;
  weight_volume += rhs.weight_volume;
  weight_mass_cs2 += rhs.weight_mass_cs2;
  add3(weight_momentum, rhs.weight_momentum);
  add3(weight_angular_momentum, rhs.weight_angular_momentum);
  weight_mass_radius += rhs.weight_mass_radius;
  mass += rhs.mass;
  cells += rhs.cells;
}

void BhKernelStatistics::merge(const BhKernelStatistics& rhs) {
  hot.merge(rhs.hot);
  cold.merge(rhs.cold);
  star_mass += rhs.star_mass;
  add3(star_mass_position, rhs.star_mass_position);
  add3(star_momentum, rhs.star_momentum);
  add3(star_r_cross_p, rhs.star_r_cross_p);
  stars += rhs.stars;
  add3(cold_mass_position, rhs.cold_mass_position);
  add3(cold_momentum, rhs.cold_momentum);
  add3(cold_r_cross_p, rhs.cold_r_cross_p);
  cold_counter_rotating_mass += rhs.cold_counter_rotating_mass;
  star_counter_rotating_mass += rhs.star_counter_rotating_mass;
  if (rhs.min_cell_width_proper > 0.0) {
    min_cell_width_proper = min_cell_width_proper > 0.0
        ? std::min(min_cell_width_proper, rhs.min_cell_width_proper)
        : rhs.min_cell_width_proper;
  }
}

void BhRadiusBinSums::merge(const BhRadiusBinSums& rhs) {
  for (std::size_t i = 0; i < k_bh_kernel_radius_bins; ++i) {
    enclosed_gas_mass[i] += rhs.enclosed_gas_mass[i];
    enclosed_cells[i] += rhs.enclosed_cells[i];
  }
}

std::array<double, k_bh_kernel_radius_bins> bhKernelRadiusCandidates(double r_min, double r_max) {
  if (!(r_min > 0.0) || !(r_max >= r_min) || !std::isfinite(r_max)) {
    throw std::invalid_argument("bhKernelRadiusCandidates: requires 0 < r_min <= r_max < inf");
  }
  std::array<double, k_bh_kernel_radius_bins> radii{};
  const double ratio = r_max / r_min;
  for (std::size_t i = 0; i < k_bh_kernel_radius_bins; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(k_bh_kernel_radius_bins - 1U);
    radii[i] = r_min * std::pow(ratio, t);
  }
  radii.back() = r_max;
  return radii;
}

std::size_t bhSelectKernelRadius(const BhRadiusBinSums& sums, double target_mass_code) {
  for (std::size_t i = 0; i < k_bh_kernel_radius_bins; ++i) {
    if (sums.enclosed_gas_mass[i] >= target_mass_code) return i;
  }
  return k_bh_kernel_radius_bins - 1U;
}

std::array<double, 3> bhDiskAxis(const BhKernelStatistics& s) {
  const double m = s.star_mass + s.cold.mass;
  if (!(m > 0.0)) return {0.0, 0.0, 0.0};
  std::array<double, 3> v_bulk{};
  for (std::size_t a = 0; a < 3; ++a) v_bulk[a] = (s.star_momentum[a] + s.cold_momentum[a]) / m;
  std::array<double, 3> mr = s.star_mass_position;
  add3(mr, s.cold_mass_position);
  std::array<double, 3> j = s.star_r_cross_p;
  add3(j, s.cold_r_cross_p);
  const std::array<double, 3> correction = cross(mr, v_bulk);
  for (std::size_t a = 0; a < 3; ++a) j[a] -= correction[a];
  const double n = norm(j);
  if (!(n > 0.0) || !std::isfinite(n)) return {0.0, 0.0, 0.0};
  return {j[0] / n, j[1] / n, j[2] / n};
}

double bondiHoyleRate(double alpha, double newton_g, double bh_mass, double density, double cs, double v) {
  const double denom2 = cs * cs + v * v;
  if (!(denom2 > 0.0) || !std::isfinite(denom2)) {
    throw std::domain_error("bondiHoyleRate: degenerate c_s^2 + v_rel^2 (no physical EOS floor)");
  }
  return 4.0 * core::constants::k_pi * alpha * newton_g * newton_g * bh_mass * bh_mass * density /
      (denom2 * std::sqrt(denom2));
}

double rosasGuevaraSuppression(double cs, double v_phi, double c_visc) noexcept {
  if (!(c_visc > 0.0)) return 1.0;
  if (!(v_phi > 0.0)) return 1.0;
  const double ratio = cs / v_phi;
  return std::min(1.0, ratio * ratio * ratio / c_visc);
}

double hopkinsQuataertTorqueRateMsunPerYear(double eps_t, double f_d, double m_bh_msun, double m_enc_msun,
                                            double r0_pc, double m_d_msun, double f_gas) noexcept {
  if (!(eps_t > 0.0) || !(f_d > 0.0) || !(m_d_msun > 0.0) || !(f_gas > 0.0) || !(r0_pc > 0.0) ||
      !(m_bh_msun > 0.0) || !(m_enc_msun > 0.0)) {
    return 0.0;
  }
  const double f0 = 0.31 * f_d * f_d * std::pow(m_d_msun / 1.0e9, -1.0 / 3.0);
  return eps_t * std::pow(f_d, 2.5) * std::pow(m_bh_msun / 1.0e8, 1.0 / 6.0) * (m_enc_msun / 1.0e9) *
      std::pow(r0_pc / 100.0, -1.5) / (1.0 + f0 / f_gas);
}

BhSupplyEstimate estimateBhSupply(const BhKernelStatistics& s, const BhSupplyInputs& in,
                                  const ProperFrameContext& frame) {
  BhSupplyEstimate out;
  const double g = frame.newtonGCode();
  out.kernel_adequate = (s.hot.cells + s.cold.cells) >= in.kernel_min_cells;
  out.enclosed_baryon_mass_code = s.hot.mass + s.cold.mass + s.star_mass;

  if (in.use_bondi && s.hot.weight_mass > 0.0 && s.hot.weight_volume > 0.0) {
    out.hot_density_proper = s.hot.weight_mass / s.hot.weight_volume;
    out.hot_sound_speed = std::sqrt(std::max(s.hot.weight_mass_cs2 / s.hot.weight_mass, 0.0));
    std::array<double, 3> dv{};
    for (std::size_t a = 0; a < 3; ++a) {
      dv[a] = in.bh_velocity_peculiar[a] - s.hot.weight_momentum[a] / s.hot.weight_mass;
    }
    out.hot_relative_speed = norm(dv);
    const double denom2 = out.hot_sound_speed * out.hot_sound_speed + out.hot_relative_speed * out.hot_relative_speed;
    if (!(denom2 > 0.0)) {
      out.bondi_estimator_rejected = true;
    } else {
      out.bondi_unsuppressed_rate_code = bondiHoyleRate(in.alpha, g, in.bh_mass_code, out.hot_density_proper,
                                                        out.hot_sound_speed, out.hot_relative_speed);
      out.bondi_radius_proper = g * in.bh_mass_code / denom2;
      if (s.min_cell_width_proper > 0.0) out.bondi_radius_over_dx = out.bondi_radius_proper / s.min_cell_width_proper;
      if (in.softening_proper > 0.0) out.bondi_radius_over_softening = out.bondi_radius_proper / in.softening_proper;
      if (in.c_visc > 0.0 && s.hot.weight_mass_radius > 0.0) {
        out.hot_rotation_speed = norm(s.hot.weight_angular_momentum) / s.hot.weight_mass_radius;
        out.angular_suppression_factor =
            rosasGuevaraSuppression(out.hot_sound_speed, out.hot_rotation_speed, in.c_visc);
      }
      out.hot_rate_code = out.bondi_unsuppressed_rate_code * out.angular_suppression_factor;
    }
  }

  if (in.use_torque && in.torque_epsilon > 0.0) {
    const double m_disk_gas = std::max(0.0, s.cold.mass - 2.0 * s.cold_counter_rotating_mass);
    const double m_disk_star = std::max(0.0, s.star_mass - 2.0 * s.star_counter_rotating_mass);
    out.disk_mass_code = m_disk_gas + m_disk_star;
    out.disk_fraction = out.enclosed_baryon_mass_code > 0.0 ? out.disk_mass_code / out.enclosed_baryon_mass_code : 0.0;
    out.disk_gas_fraction = out.disk_mass_code > 0.0 ? m_disk_gas / out.disk_mass_code : 0.0;
    const std::array<double, 3> axis = bhDiskAxis(s);
    out.torque_adequate = (s.cold.cells + s.stars) >= in.torque_min_samples && norm(axis) > 0.0;
    if (out.torque_adequate) {
      const double msun_per_code = frame.mass_si_per_code / core::constants::k_solar_mass_si;
      const double r0_pc = in.kernel_radius_proper * frame.length_si_per_code / core::constants::k_parsec_si;
      const double rate_msun_yr = hopkinsQuataertTorqueRateMsunPerYear(
          in.torque_epsilon, out.disk_fraction, in.bh_mass_code * msun_per_code,
          out.enclosed_baryon_mass_code * msun_per_code, r0_pc, out.disk_mass_code * msun_per_code,
          out.disk_gas_fraction);
      const double code_mass_per_msun = 1.0 / msun_per_code;
      const double code_time_per_year = k_year_si / frame.timeSiPerCode();
      out.cold_rate_code = rate_msun_yr * code_mass_per_msun / code_time_per_year;
      // The torque supply can only draw on cold gas actually present.
      if (!(s.cold.mass > 0.0)) out.cold_rate_code = 0.0;
    }
  }
  return out;
}

BhResolvedFlux resolvedControlSurfaceFlux(const BhResolvedFluxSample* samples, std::size_t count, double r0,
                                          double dr) noexcept {
  BhResolvedFlux out;
  if (!(dr > 0.0) || samples == nullptr) return out;
  for (std::size_t i = 0; i < count; ++i) {
    const BhResolvedFluxSample& s = samples[i];
    if (std::abs(s.radius_proper - r0) > 0.5 * dr) continue;
    out.inward_rate_code += s.mass_code * std::max(-s.radial_velocity_relative, 0.0) / dr;
    out.net_rate_code += -s.mass_code * s.radial_velocity_relative / dr;
    ++out.samples;
  }
  return out;
}

}  // namespace cosmosim::physics
