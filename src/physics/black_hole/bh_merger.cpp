#include "cosmosim/physics/black_hole/bh_merger.hpp"

#include <algorithm>
#include <cmath>
#include <tuple>

namespace cosmosim::physics {

bool bhPairMergeEligible(const BhMergerCandidate& a, const BhMergerCandidate& b, const BhMergerCriteria& c,
                         const ProperFrameContext& frame, const PeriodicBox& box, double* energy_out,
                         double* separation_out) noexcept {
  if (a.particle_id == b.particle_id || !(a.dynamical_mass_code > 0.0) || !(b.dynamical_mass_code > 0.0)) {
    return false;
  }
  const std::array<double, 3> d_stored = box.delta(a.position_stored, b.position_stored);
  const std::array<double, 3> dr{frame.properLength(d_stored[0]), frame.properLength(d_stored[1]),
                                 frame.properLength(d_stored[2])};
  // Peculiar velocity difference plus Hubble flow across the proper separation.
  const double h = frame.hubbleCode();
  const std::array<double, 3> dv{b.velocity_peculiar[0] - a.velocity_peculiar[0] + h * dr[0],
                                 b.velocity_peculiar[1] - a.velocity_peculiar[1] + h * dr[1],
                                 b.velocity_peculiar[2] - a.velocity_peculiar[2] + h * dr[2]};
  const double r = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
  const double eps = std::max(a.softening_proper, b.softening_proper);
  if (!(eps > 0.0) || !(r < c.radius_softening_factor * eps)) return false;
  const double m1 = a.dynamical_mass_code;
  const double m2 = b.dynamical_mass_code;
  const double mu = m1 * m2 / (m1 + m2);
  const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
  const double energy = 0.5 * mu * v2 - c.newton_g_code * m1 * m2 / std::sqrt(r * r + eps * eps);
  if (!(energy < 0.0)) return false;
  const double radial = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
  if (!(radial < 0.0) && !(r < eps)) return false;
  if (energy_out != nullptr) *energy_out = energy;
  if (separation_out != nullptr) *separation_out = r;
  return true;
}

std::vector<BhMergerPair> selectDisjointMergerPairs(std::span<const BhMergerCandidate> cand,
                                                    const BhMergerCriteria& c, const ProperFrameContext& frame,
                                                    const PeriodicBox& box) {
  std::vector<BhMergerPair> all;
  for (std::size_t i = 0; i < cand.size(); ++i) {
    for (std::size_t j = i + 1; j < cand.size(); ++j) {
      double e = 0.0;
      double r = 0.0;
      if (!bhPairMergeEligible(cand[i], cand[j], c, frame, box, &e, &r)) continue;
      const bool i_survives = cand[i].dynamical_mass_code > cand[j].dynamical_mass_code ||
          (cand[i].dynamical_mass_code == cand[j].dynamical_mass_code && cand[i].particle_id < cand[j].particle_id);
      BhMergerPair pair;
      pair.survivor_index = static_cast<std::uint32_t>(i_survives ? i : j);
      pair.absorbed_index = static_cast<std::uint32_t>(i_survives ? j : i);
      pair.survivor_id = cand[pair.survivor_index].particle_id;
      pair.absorbed_id = cand[pair.absorbed_index].particle_id;
      pair.binding_energy_code = e;
      pair.separation_proper = r;
      all.push_back(pair);
    }
  }
  std::sort(all.begin(), all.end(), [](const BhMergerPair& x, const BhMergerPair& y) {
    const auto kx = std::make_tuple(x.binding_energy_code, std::min(x.survivor_id, x.absorbed_id),
                                    std::max(x.survivor_id, x.absorbed_id));
    const auto ky = std::make_tuple(y.binding_energy_code, std::min(y.survivor_id, y.absorbed_id),
                                    std::max(y.survivor_id, y.absorbed_id));
    return kx < ky;
  });
  std::vector<std::uint8_t> used(cand.size(), 0U);
  std::vector<BhMergerPair> chosen;
  for (const BhMergerPair& p : all) {
    if (used[p.survivor_index] != 0U || used[p.absorbed_index] != 0U) continue;
    used[p.survivor_index] = 1U;
    used[p.absorbed_index] = 1U;
    chosen.push_back(p);
  }
  return chosen;
}

BlackHoleEngineRecord mergeBlackHoleEngineRecords(const BlackHoleEngineRecord& s, const BlackHoleEngineRecord& a,
                                                  std::uint64_t absorbed_id) {
  BlackHoleEngineRecord m = s;
  m.horizon_mass_code += a.horizon_mass_code;
  m.disk_mass_code += a.disk_mass_code;
  m.disk_metal_mass_code += a.disk_metal_mass_code;
  m.pending_wind_metal_mass_code += a.pending_wind_metal_mass_code;
  for (std::size_t e = 0; e < m.disk_element_mass_code.size(); ++e) {
    m.disk_element_mass_code[e] += a.disk_element_mass_code[e];
    m.pending_wind_element_mass_code[e] += a.pending_wind_element_mass_code[e];
  }
  // Internal angular momenta of the components; the ORBITAL angular momentum
  // of the pair is NOT added to the horizon spin (that needs an NR remnant
  // fit; P1). It is accounted in the merger diagnostic by the caller.
  for (std::size_t k = 0; k < 3; ++k) {
    m.disk_angular_momentum_code[k] += a.disk_angular_momentum_code[k];
    m.spin_angular_momentum_code[k] += a.spin_angular_momentum_code[k];
  }
  m.pending_thermal_energy_code += a.pending_thermal_energy_code;
  m.pending_kick_energy_code += a.pending_kick_energy_code;
  m.pending_wind_energy_code += a.pending_wind_energy_code;
  m.pending_wind_mass_code += a.pending_wind_mass_code;
  m.pending_jet_energy_code += a.pending_jet_energy_code;
  m.pending_radiation_energy_code += a.pending_radiation_energy_code;
  m.pending_magnetic_energy_code += a.pending_magnetic_energy_code;
  m.cumulative_supply_requested_mass_code += a.cumulative_supply_requested_mass_code;
  m.cumulative_withdrawn_mass_code += a.cumulative_withdrawn_mass_code;
  m.cumulative_unmet_supply_mass_code += a.cumulative_unmet_supply_mass_code;
  m.cumulative_inner_inflow_mass_code += a.cumulative_inner_inflow_mass_code;
  m.cumulative_disk_wind_mass_code += a.cumulative_disk_wind_mass_code;
  m.cumulative_mass_defect_code += a.cumulative_mass_defect_code;
  m.cumulative_radiated_energy_code += a.cumulative_radiated_energy_code;
  m.cumulative_mechanical_energy_code += a.cumulative_mechanical_energy_code;
  m.cumulative_deposited_energy_code += a.cumulative_deposited_energy_code;
  m.cumulative_escaped_radiation_energy_code += a.cumulative_escaped_radiation_energy_code;
  m.cumulative_captured_gas_energy_code += a.cumulative_captured_gas_energy_code;
  m.cumulative_capture_dissipation_code += a.cumulative_capture_dissipation_code;
  m.cumulative_df_reservoir_energy_code += a.cumulative_df_reservoir_energy_code;
  m.cumulative_spin_extracted_energy_code += a.cumulative_spin_extracted_energy_code;
  m.merger_count = s.merger_count + a.merger_count + 1U;
  m.last_merged_partner_id = absorbed_id;
  m.flags |= kBhFlagMergedRemnant;
  return m;
}

double petersCircularInspiralTime(double a, double m1, double m2, double g, double c) noexcept {
  if (!(a > 0.0) || !(m1 > 0.0) || !(m2 > 0.0) || !(g > 0.0) || !(c > 0.0)) return 0.0;
  const double c5 = c * c * c * c * c;
  const double a4 = a * a * a * a;
  return 5.0 / 256.0 * c5 * a4 / (g * g * g * m1 * m2 * (m1 + m2));
}

}  // namespace cosmosim::physics
