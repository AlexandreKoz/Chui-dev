// CHUI-EXT-PHYS-REF-v1 S03 lightcone planning (source only, NOT compiled in pass 01).

#include "cosmosim/lightcone/lightcone_plan.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::lightcone {

ComovingDistanceTable::ComovingDistanceTable(const core::LambdaCdmBackground& background, double a_min,
                                             std::size_t points)
    : m_a_min(a_min) {
  if (!(a_min > 0.0) || !(a_min < 1.0) || points < 16U) throw std::invalid_argument("invalid comoving distance table");
  m_ln_a.resize(points);
  m_chi.resize(points);
  const double l0 = std::log(a_min);
  for (std::size_t i = 0; i < points; ++i) {
    m_ln_a[i] = l0 * (1.0 - static_cast<double>(i) / static_cast<double>(points - 1U));
  }
  m_ln_a.back() = 0.0;
  // Cumulative from a = 1 inward, segment by segment (each segment integrated
  // by the background authority's adaptive conformal-time integral).
  m_chi.back() = 0.0;
  for (std::size_t i = points - 1U; i > 0; --i) {
    const double a_lo = std::exp(m_ln_a[i - 1U]);
    const double a_hi = std::exp(m_ln_a[i]);
    m_chi[i - 1U] = m_chi[i] + core::constants::k_speed_of_light_si * background.conformalTimeIntervalSi(a_lo, a_hi);
  }
}

double ComovingDistanceTable::comovingDistanceM(double a) const {
  if (!(a >= m_a_min) || !(a <= 1.0)) throw std::out_of_range("scale factor outside the comoving distance table");
  const double x = std::log(a);
  const auto it = std::upper_bound(m_ln_a.begin(), m_ln_a.end(), x);
  auto hi = static_cast<std::size_t>(std::distance(m_ln_a.begin(), it));
  hi = std::clamp<std::size_t>(hi, 1U, m_ln_a.size() - 1U);
  const std::size_t lo = hi - 1U;
  const double t = (x - m_ln_a[lo]) / (m_ln_a[hi] - m_ln_a[lo]);
  // chi is smooth in ln a; cubic Hermite with finite-difference slopes.
  const auto slope = [&](std::size_t i) {
    const std::size_t l = i == 0U ? 0U : i - 1U;
    const std::size_t r = std::min(i + 1U, m_ln_a.size() - 1U);
    return (m_chi[r] - m_chi[l]) / (m_ln_a[r] - m_ln_a[l]);
  };
  const double h = m_ln_a[hi] - m_ln_a[lo];
  const double t2 = t * t;
  const double t3 = t2 * t;
  return (2.0 * t3 - 3.0 * t2 + 1.0) * m_chi[lo] + (t3 - 2.0 * t2 + t) * h * slope(lo) +
         (-2.0 * t3 + 3.0 * t2) * m_chi[hi] + (t3 - t2) * h * slope(hi);
}

double ComovingDistanceTable::scaleFactorAt(double chi) const {
  if (!(chi >= 0.0) || !(chi <= m_chi.front())) throw std::out_of_range("comoving distance outside the table");
  // chi decreases monotonically with ln a: bisection on the interpolant.
  double lo = std::log(m_a_min);
  double hi = 0.0;
  for (int it = 0; it < 200 && hi - lo > 1.0e-15; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (comovingDistanceM(std::exp(mid)) > chi) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return std::exp(0.5 * (lo + hi));
}

std::vector<LightconeShell> planLightconeShells(const ComovingDistanceTable& table,
                                                const std::vector<double>& snapshot_scale_factors, double chi_min_m) {
  if (snapshot_scale_factors.empty()) throw std::invalid_argument("lightcone needs at least one snapshot");
  for (std::size_t i = 1; i < snapshot_scale_factors.size(); ++i) {
    if (!(snapshot_scale_factors[i] > snapshot_scale_factors[i - 1U])) {
      throw std::invalid_argument("snapshot scale factors must be strictly ascending");
    }
  }
  const std::size_t n = snapshot_scale_factors.size();
  std::vector<double> chi(n);
  for (std::size_t i = 0; i < n; ++i) chi[i] = table.comovingDistanceM(snapshot_scale_factors[i]);
  std::vector<LightconeShell> shells(n);
  for (std::size_t i = 0; i < n; ++i) {
    LightconeShell& shell = shells[i];
    shell.snapshot_index = i;
    shell.snapshot_scale_factor = snapshot_scale_factors[i];
    shell.chi_outer_m = i == 0U ? chi[0] : 0.5 * (chi[i] + chi[i - 1U]);
    shell.chi_inner_m = i + 1U == n ? chi_min_m : 0.5 * (chi[i] + chi[i + 1U]);
    if (!(shell.chi_outer_m > shell.chi_inner_m)) throw std::invalid_argument("degenerate lightcone shell");
    shell.a_outer = table.scaleFactorAt(shell.chi_outer_m);
    shell.a_inner = shell.chi_inner_m == 0.0 ? 1.0 : table.scaleFactorAt(shell.chi_inner_m);
  }
  return shells;
}

ReplicationMetadata replicationMetadata(double box_size_m, double chi_max_m) {
  if (!(box_size_m > 0.0) || !(chi_max_m >= 0.0)) throw std::invalid_argument("invalid replication request");
  ReplicationMetadata meta;
  meta.box_size_m = box_size_m;
  meta.chi_max_m = chi_max_m;
  meta.replicas_along_ray = static_cast<std::uint64_t>(std::ceil(chi_max_m / box_size_m));
  meta.repeated_structures = meta.replicas_along_ray > 1U;
  meta.missing_mode_k_max_per_m = 2.0 * core::constants::k_pi / box_size_m;
  return meta;
}

TemporalWeights temporalWeights(const std::vector<double>& a, double a_emit) {
  if (a.size() < 2U) throw std::invalid_argument("temporal interpolation needs >= 2 snapshots");
  if (!(a_emit >= a.front()) || !(a_emit <= a.back())) {
    throw std::out_of_range("emission scale factor outside the snapshot range (no temporal extrapolation)");
  }
  const auto it = std::upper_bound(a.begin(), a.end(), a_emit);
  std::size_t hi = static_cast<std::size_t>(std::distance(a.begin(), it));
  hi = std::clamp<std::size_t>(hi, 1U, a.size() - 1U);
  TemporalWeights w;
  w.lower_snapshot = hi - 1U;
  w.upper_snapshot = hi;
  w.upper_weight = (a_emit - a[hi - 1U]) / (a[hi] - a[hi - 1U]);
  return w;
}

double observedRedshift(double z_cos, double v_los) {
  return (1.0 + z_cos) * (1.0 + v_los / core::constants::k_speed_of_light_si) - 1.0;
}

double hubbleVelocityMS(const core::LambdaCdmBackground& background, double a, double dx) {
  return a * background.hubbleSi(a) * dx;
}

}  // namespace cosmosim::lightcone
