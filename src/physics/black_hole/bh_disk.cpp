#include "cosmosim/physics/black_hole/bh_disk.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cosmosim::physics {

double oneMinusExpNeg(double x) noexcept { return -std::expm1(-x); }

BhDiskStepResult integrateScalarDisk(const BhDiskStepInputs& in) {
  if (!(in.viscous_time_code > 0.0) || !std::isfinite(in.viscous_time_code) || !(in.wind_loading >= 0.0) ||
      !std::isfinite(in.wind_loading) || !(in.supply_rate_code >= 0.0) || !std::isfinite(in.supply_rate_code) ||
      !(in.disk_mass_code >= 0.0) || !std::isfinite(in.disk_mass_code) || !(in.dt_code >= 0.0) ||
      !std::isfinite(in.dt_code) || !(in.inner_rate_cap_code > 0.0)) {
    throw std::invalid_argument("integrateScalarDisk: invalid disk inputs (requires M>=0, S>=0, t_visc>0, "
                                "eta_w>=0, dt>=0, cap>0)");
  }
  BhDiskStepResult out;
  const double one_plus_eta = 1.0 + in.wind_loading;
  const double tau = in.viscous_time_code / one_plus_eta;
  const double s = in.supply_rate_code;
  const double cap = in.inner_rate_cap_code;
  const double m_star = std::isfinite(cap) ? cap * in.viscous_time_code : std::numeric_limits<double>::infinity();
  const double k_capped = std::isfinite(cap) ? s - one_plus_eta * cap : 0.0;
  double m = in.disk_mass_code;
  double remaining = in.dt_code;
  double inner = 0.0;
  double wind = 0.0;
  // At most three analytic pieces: U -> C, C -> U, or a single regime.
  for (std::uint32_t piece = 0; piece < 4U && remaining > 0.0; ++piece) {
    ++out.pieces;
    const bool capped = std::isfinite(cap) && (m > m_star || (m == m_star && k_capped > 0.0));
    if (capped) {
      double t = remaining;
      if (k_capped < 0.0) t = std::min(remaining, (m - m_star) / (-k_capped));
      if (piece == 3U) t = remaining;
      const double m_next = (t == remaining || k_capped >= 0.0) ? m + k_capped * t : m_star;
      inner += cap * t;
      wind += in.wind_loading * cap * t;
      out.capped_time_code += t;
      m = std::max(m_next, 0.0);
      remaining -= t;
    } else {
      const double equilibrium = s * tau;
      double t = remaining;
      if (std::isfinite(cap) && equilibrium > m_star && m < m_star) {
        const double crossing = tau * std::log((equilibrium - m) / (equilibrium - m_star));
        if (crossing < remaining) t = crossing;
      }
      if (piece == 3U) t = remaining;
      const double decay = oneMinusExpNeg(t / tau);  // 1 - e^{-t/tau}
      double m_next = m * (1.0 - decay) + equilibrium * decay;
      if (t < remaining) m_next = m_star;  // exact crossing state (removes log round-off)
      const double drained = m + s * t - m_next;
      inner += drained / one_plus_eta;
      wind += in.wind_loading * drained / one_plus_eta;
      m = std::max(m_next, 0.0);
      remaining -= t;
    }
  }
  out.disk_mass_end_code = m;
  out.drained_mass_code = in.disk_mass_code + s * in.dt_code - m;
  // Assign the exact ledger: inner + wind == drained to round-off.
  const double piece_total = inner + wind;
  if (piece_total > 0.0) {
    const double scale = out.drained_mass_code / piece_total;
    inner *= scale;
    wind *= scale;
  }
  out.inner_mass_code = std::max(inner, 0.0);
  out.wind_mass_code = std::max(wind, 0.0);
  return out;
}

}  // namespace cosmosim::physics
