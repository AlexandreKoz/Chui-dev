#include "cosmosim/physics/agn/agn_injection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <optional>
#include <stdexcept>
#include <utility>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/black_hole/bh_engine_state.hpp"
#include "cosmosim/physics/black_hole/bh_supply.hpp"
#include "cosmosim/physics/source_transaction.hpp"

namespace cosmosim::physics {
namespace {

using Vec3 = std::array<double, 3>;
[[nodiscard]] double dot(const Vec3& a, const Vec3& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
[[nodiscard]] double norm(const Vec3& a) noexcept { return std::sqrt(dot(a, a)); }

// Statistics accumulator slots (round kAgnStatistics). "max" slots merge by max.
enum Acc : std::size_t {
  kSumOmega = 0,
  kSumWeightMass = 1,
  kSumWeightMassCv = 2,
  kKickNet = 3,  // 3..5
  kKickB = 6,
  kKickC = 7,
  kLobePlusOmega = 8,
  kLobePlusB = 9,
  kLobePlusC = 10,
  kLobeMinusOmega = 11,
  kLobeMinusB = 12,
  kLobeMinusC = 13,
  kWindDirection = 14,  // 14..16: sum w m r_hat
  kCells = 17,
  kMaxOmegaOverMass = 18,         // max slot
  kMaxLobeOmegaOverMass = 19,     // max slot
  kMinCellWidth = 20,             // min slot (stored as -width for max merge)
  // Commit-round replies.
  kGasEnergy = 24,
  kGasMomentum = 25,  // 25..27
  kGasMass = 28,
  kThermalShare = 29,
  kKickShare = 30,
  kJetShare = 31,
  kWindShare = 32,
};
// Commit parameters.
enum Par : std::size_t {
  kCosHalfAngle = 0,
  kThermalQ = 1,
  kThermalReservoir = 2,
  kThermalNorm = 3,
  kKickA = 4,
  kJetA = 5,
  kLobePlusNorm = 6,
  kLobeMinusNorm = 7,
  kWindMass = 8,
  kWindVelocity = 9,
  kWindSpecificThermal = 10,
  kWindNorm = 11,
  kWindElementFraction = 12,  // 12..22
  kWindMetalFraction = 23,
};

struct Pending {
  std::uint32_t bh_row = 0U;
  std::uint32_t particle_row = 0U;
  std::uint64_t id = 0U;
  Vec3 position{};
  Vec3 velocity{};
  BlackHoleEngineRecord record{};
  std::array<double, k_bh_reply_accumulator_count> stats{};
  double thermal_q = 0.0;
  double kick_q = 0.0;
  double jet_q = 0.0;
  double wind_q = 0.0;
  double wind_mass = 0.0;
};

void mergeStats(std::array<double, k_bh_reply_accumulator_count>& into,
                const std::array<double, k_bh_reply_accumulator_count>& from, bool first) {
  for (std::size_t k = 0; k < k_bh_reply_accumulator_count; ++k) {
    if (k == kMaxOmegaOverMass || k == kMaxLobeOmegaOverMass || k == kMinCellWidth) {
      into[k] = first ? from[k] : std::max(into[k], from[k]);
    } else {
      into[k] += from[k];
    }
  }
}

}  // namespace

bool solveEnergyConsistentAmplitude(double c, double b, double q, double& a) noexcept {
  a = 0.0;
  if (!(q >= 0.0) || !(c >= 0.0) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(q)) return false;
  if (q == 0.0) return true;
  if (c == 0.0) {
    if (b > 0.0) {
      a = q / b;
      return true;
    }
    return false;
  }
  const double disc = b * b + 4.0 * c * q;  // >= b^2 >= 0
  const double root = std::sqrt(disc);
  // Cancellation-free form of the nonnegative root.
  a = b >= 0.0 ? 2.0 * q / (b + root) : (root - b) / (2.0 * c);
  return std::isfinite(a) && a >= 0.0;
}

double cellKineticIncrement(double m, const std::array<double, 3>& p, double dm, const std::array<double, 3>& dp) noexcept {
  const Vec3 p2{p[0] + dp[0], p[1] + dp[1], p[2] + dp[2]};
  return dot(p2, p2) / (2.0 * (m + dm)) - dot(p, p) / (2.0 * m);
}

AgnFeedbackInjector::AgnFeedbackInjector(core::AgnFeedbackConfig config) : m_config(std::move(config)) {
  if (m_config.model != core::AgnFeedbackModelConfig::kChannelsV1) {
    throw std::invalid_argument("AgnFeedbackInjector requires physics.agn_feedback_model=channels_v1");
  }
}

AgnInjectionReport AgnFeedbackInjector::apply(core::SimulationState& state, const AgnInjectionContext& ctx,
                                              BhSupportExchange& exchange) const {
  ctx.frame.validate();
  AgnInjectionReport report;
  BlackHoleEngineStore store(state);
  const std::uint32_t rank = exchange.worldRank();
  const double c_light = ctx.frame.speedOfLightCode();
  const double v_max = m_config.max_velocity_c_fraction * c_light;
  const double cos_half = std::cos(m_config.jet_half_opening_deg * core::constants::k_pi / 180.0);
  const double km_s_code = 1.0e3 / ctx.frame.velocity_si_per_code;

  std::vector<Pending> sources;
  for (std::size_t b = 0; b < state.black_holes.size(); ++b) {
    Pending s;
    s.bh_row = static_cast<std::uint32_t>(b);
    s.particle_row = state.black_holes.particle_index[b];
    s.id = state.particle_sidecar.particle_id[s.particle_row];
    if (!store.has(s.id)) continue;
    s.record = store.get(s.id);
    s.position = {state.particles.position_x_comoving[s.particle_row], state.particles.position_y_comoving[s.particle_row],
                  state.particles.position_z_comoving[s.particle_row]};
    s.velocity = {state.particles.velocity_x_peculiar[s.particle_row], state.particles.velocity_y_peculiar[s.particle_row],
                  state.particles.velocity_z_peculiar[s.particle_row]};
    // Radiation channel disposition (exactly once).
    BlackHoleEngineRecord& r = s.record;
    switch (m_config.radiation_coupling) {
      case core::AgnRadiationCouplingConfig::kNone:
        r.cumulative_escaped_radiation_energy_code += r.pending_radiation_energy_code;
        report.radiation_escaped_code += r.pending_radiation_energy_code;
        r.pending_radiation_energy_code = 0.0;
        break;
      case core::AgnRadiationCouplingConfig::kLocalEmpirical: {
        const double coupled = m_config.local_empirical_fraction * r.pending_radiation_energy_code;
        const double escaped = r.pending_radiation_energy_code - coupled;
        r.pending_thermal_energy_code += coupled;
        r.cumulative_escaped_radiation_energy_code += escaped;
        report.radiation_locally_coupled_code += coupled;
        report.radiation_escaped_code += escaped;
        r.pending_radiation_energy_code = 0.0;
        break;
      }
      case core::AgnRadiationCouplingConfig::kRtTransport:
        break;  // consumed by the RT AGN source provider
    }
    const bool any = r.pending_thermal_energy_code > 0.0 || r.pending_kick_energy_code > 0.0 ||
        r.pending_jet_energy_code > 0.0 || r.pending_wind_energy_code > 0.0 || r.pending_wind_mass_code > 0.0;
    if (!any) {
      store.put(s.id, r);
      continue;
    }
    sources.push_back(s);
  }
  std::sort(sources.begin(), sources.end(), [](const Pending& a, const Pending& b) { return a.id < b.id; });
  report.sources = sources.size();

  // Round A: support statistics.
  std::vector<BhSupportQuery> queries;
  for (std::size_t i = 0; i < sources.size(); ++i) {
    const Pending& s = sources[i];
    BhSupportQuery q;
    q.bh_id = s.id;
    q.query_index = i;
    q.origin_rank = rank;
    q.round = BhSupportRound::kAgnStatistics;
    q.position_stored = s.position;
    q.velocity_peculiar = s.velocity;
    q.kernel_radius_proper = s.record.last_kernel_radius_proper_code;
    q.support_radius_stored = ctx.frame.storedLength(q.kernel_radius_proper);
    const double jn = norm(s.record.disk_angular_momentum_code);
    const double sn = norm(s.record.spin_angular_momentum_code);
    const Vec3& j = sn > 0.0 ? s.record.spin_angular_momentum_code : s.record.disk_angular_momentum_code;
    const double n = sn > 0.0 ? sn : jn;
    q.axis = n > 0.0 ? Vec3{j[0] / n, j[1] / n, j[2] / n} : Vec3{0.0, 0.0, 0.0};
    q.bh_mass_code = state.particles.mass_code[s.particle_row];
    q.parameters[kCosHalfAngle] = cos_half;
    queries.push_back(q);
  }
  auto evaluate = [&](std::span<const BhSupportQuery> routed, bool commit_round,
                      SourceTransactionBatch* batch) -> std::vector<BhSupportReply> {
    std::vector<BhSupportReply> replies(routed.size());
    std::vector<SourceTransactionPlan> plans(routed.size());
    for (std::size_t i = 0; i < routed.size(); ++i) {
      replies[i].bh_id = routed[i].bh_id;
      replies[i].query_index = routed[i].query_index;
      replies[i].round = routed[i].round;
      replies[i].responder_rank = rank;
      replies[i].accumulators[kMinCellWidth] = -std::numeric_limits<double>::max();
      plans[i].id = SourceTransactionId{routed[i].bh_id, ctx.interval_index, SourceChannel::kAgnThermal, rank};
      plans[i].topology_epoch = state.cellIndexGeneration();
      plans[i].owner_rank = routed[i].origin_rank;
      plans[i].intended_time_code = ctx.time_end_code;
    }
    if (!routed.empty()) {
      const BhProximityFilter filter(routed, ctx.box);
      std::vector<std::uint32_t> hits;
      for (std::size_t row = 0; row < state.cells.size(); ++row) {
        if (!ctx.leaves.owned(row)) continue;
        const Vec3 x{state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                     state.cells.center_z_comoving[row]};
        filter.matches(x, hits);
        if (hits.empty()) continue;
        std::array<double, 3> hw_s{};
        if (!bhCellHalfWidthStored(state, static_cast<std::uint32_t>(row), hw_s)) continue;
        const Vec3 hw{ctx.frame.properLength(hw_s[0]), ctx.frame.properLength(hw_s[1]), ctx.frame.properLength(hw_s[2])};
        const double volume = 8.0 * hw[0] * hw[1] * hw[2];
        const double m = state.cells.mass_code[row];
        const Vec3 v{state.gas_cells.velocity_x_peculiar[row], state.gas_cells.velocity_y_peculiar[row],
                     state.gas_cells.velocity_z_peculiar[row]};
        const double u = state.gas_cells.internal_energy_code[row];
        const double t = state.gas_cells.temperature_code[row];
        const double cv = t > 0.0 ? u / t : 0.0;
        for (const std::uint32_t qi : hits) {
          const BhSupportQuery& q = routed[qi];
          const Vec3 d = ctx.box.delta(q.position_stored, x);
          const Vec3 r{ctx.frame.properLength(d[0]), ctx.frame.properLength(d[1]), ctx.frame.properLength(d[2])};
          const double rn = norm(r);
          if (rn > q.kernel_radius_proper || !(m > 0.0)) continue;
          const double w = bhCellKernelWeight(r, hw, q.kernel_radius_proper);
          if (!(w > 0.0)) continue;
          const double omega = w * volume;
          const Vec3 rhat = rn > 0.0 ? Vec3{r[0] / rn, r[1] / rn, r[2] / rn} : Vec3{0.0, 0.0, 0.0};
          const double mu = norm(q.axis) > 0.0 && rn > 0.0 ? dot(rhat, q.axis) : 0.0;
          const bool lobe_plus = norm(q.axis) > 0.0 && mu >= q.parameters[kCosHalfAngle];
          const bool lobe_minus = norm(q.axis) > 0.0 && -mu >= q.parameters[kCosHalfAngle];
          auto& acc = replies[qi].accumulators;
          if (!commit_round) {
            acc[kSumOmega] += omega;
            acc[kSumWeightMass] += w * m;
            acc[kSumWeightMassCv] += w * m * cv;
            for (std::size_t a = 0; a < 3; ++a) {
              acc[kKickNet + a] += omega * rhat[a];
              acc[kWindDirection + a] += w * m * rhat[a];
            }
            acc[kKickB] += omega * dot(v, rhat);
            acc[kKickC] += omega * omega / (2.0 * m);
            if (lobe_plus) {
              acc[kLobePlusOmega] += omega;
              acc[kLobePlusB] += omega * dot(v, q.axis);
              acc[kLobePlusC] += omega * omega / (2.0 * m);
            }
            if (lobe_minus) {
              acc[kLobeMinusOmega] += omega;
              acc[kLobeMinusB] -= omega * dot(v, q.axis);
              acc[kLobeMinusC] += omega * omega / (2.0 * m);
            }
            acc[kCells] += 1.0;
            acc[kMaxOmegaOverMass] = std::max(acc[kMaxOmegaOverMass], omega / m);
            if (lobe_plus || lobe_minus) acc[kMaxLobeOmegaOverMass] = std::max(acc[kMaxLobeOmegaOverMass], omega / m);
            acc[kMinCellWidth] = std::max(acc[kMinCellWidth], -2.0 * std::min({hw[0], hw[1], hw[2]}));
            continue;
          }
          // Commit round: build the combined per-cell increment of all channels.
          const auto& par = q.parameters;
          GasCellIncrement inc;
          inc.gas_cell_id = state.gas_cells.gas_cell_id[row];
          inc.local_row = static_cast<std::uint32_t>(row);
          const Vec3 p{m * v[0], m * v[1], m * v[2]};
          double de = 0.0;
          Vec3 dp{};
          double dm = 0.0;
          // Thermal.
          if (par[kThermalQ] > 0.0 && par[kThermalNorm] > 0.0) {
            const double share = par[kThermalReservoir] > 0.5 ? w * m * cv / par[kThermalNorm] : omega / par[kThermalNorm];
            de += share * par[kThermalQ];
            acc[kThermalShare] += share * par[kThermalQ];
          }
          // Kicks (radial pattern): momentum-only, energy = actual kinetic change.
          if (par[kKickA] > 0.0) {
            const Vec3 dpk{par[kKickA] * omega * rhat[0], par[kKickA] * omega * rhat[1], par[kKickA] * omega * rhat[2]};
            const double dk = cellKineticIncrement(m, p, 0.0, dpk);
            for (std::size_t a = 0; a < 3; ++a) dp[a] += dpk[a];
            de += dk;
            acc[kKickShare] += dk;
          }
          // Jets: balanced lobes.
          if (par[kJetA] > 0.0 && (lobe_plus || lobe_minus)) {
            const double norm_lobe = lobe_plus ? par[kLobePlusNorm] : par[kLobeMinusNorm];
            const double sign = lobe_plus ? 1.0 : -1.0;
            if (norm_lobe > 0.0) {
              const double s = par[kJetA] * omega / norm_lobe * sign;
              const Vec3 dpj{s * q.axis[0], s * q.axis[1], s * q.axis[2]};
              const Vec3 p_now{p[0] + dp[0], p[1] + dp[1], p[2] + dp[2]};
              const double dk = cellKineticIncrement(m, p_now, 0.0, dpj);
              for (std::size_t a = 0; a < 3; ++a) dp[a] += dpj[a];
              de += dk;
              acc[kJetShare] += dk;
            }
          }
          // Mass-loaded wind (BH rest frame radial at v_w, Galilean boost).
          if (par[kWindMass] > 0.0 && par[kWindNorm] > 0.0) {
            const double dmw = par[kWindMass] * w * m / par[kWindNorm];
            const double vw = par[kWindVelocity];
            const Vec3 dp_rest{dmw * vw * rhat[0], dmw * vw * rhat[1], dmw * vw * rhat[2]};
            const Vec3& vb = q.velocity_peculiar;
            const double e_rest = dmw * (par[kWindSpecificThermal] + 0.5 * vw * vw);
            const double e = e_rest + dot(vb, dp_rest) + 0.5 * dmw * dot(vb, vb);
            for (std::size_t a = 0; a < 3; ++a) dp[a] += dp_rest[a] + dmw * vb[a];
            dm += dmw;
            de += e;
            acc[kWindShare] += e;
            for (std::size_t el = 0; el < k_chemical_element_count; ++el) {
              inc.delta.element_mass_code[el] += dmw * par[kWindElementFraction + el];
            }
            inc.delta.metal_mass_code += dmw * par[kWindMetalFraction];
          }
          if (de == 0.0 && dm == 0.0 && dp[0] == 0.0 && dp[1] == 0.0 && dp[2] == 0.0) continue;
          inc.delta.mass_code = dm;
          inc.delta.momentum_code = dp;
          inc.delta.total_energy_code = de;
          acc[kGasEnergy] += de;
          for (std::size_t a = 0; a < 3; ++a) acc[kGasMomentum + a] += dp[a];
          acc[kGasMass] += dm;
          ++replies[qi].rows;
          plans[qi].cells.push_back(inc);
        }
      }
    }
    if (commit_round) {
      for (std::size_t qi = 0; qi < routed.size(); ++qi) {
        if (plans[qi].cells.empty()) continue;
        ConservedIncrement gas{};
        for (const GasCellIncrement& c : plans[qi].cells) gas += c.delta;
        plans[qi].source_increment = gas.negated();  // the source (BH engine budgets) supplies the increment
        replies[qi].transferred = gas;
        batch->add(std::move(plans[qi]));
      }
    }
    return replies;
  };

  {
    const auto routed = exchange.routeQueries(queries);
    const auto replies = exchange.returnReplies(evaluate(routed, false, nullptr));
    std::vector<std::uint8_t> seen(sources.size(), 0U);
    for (const BhSupportReply& r : replies) {
      if (r.query_index >= sources.size() || sources[r.query_index].id != r.bh_id) {
        throw std::runtime_error("AgnFeedbackInjector: statistics reply does not match its query");
      }
      mergeStats(sources[r.query_index].stats, r.accumulators, seen[r.query_index] == 0U);
      seen[r.query_index] = 1U;
    }
  }

  // Owner: channel amplitudes.
  for (std::size_t i = 0; i < sources.size(); ++i) {
    Pending& s = sources[i];
    BlackHoleEngineRecord& r = s.record;
    const auto& st = s.stats;
    BhSupportQuery& q = queries[i];
    q.round = BhSupportRound::kAgnCommit;
    q.parameters.fill(0.0);
    q.parameters[kCosHalfAngle] = cos_half;
    r.flags &= ~kBhFlagPendingFeedbackHeld;
    if (!(st[kCells] > 0.0) || !(st[kSumOmega] > 0.0)) {
      r.flags |= kBhFlagPendingFeedbackHeld;
      ++report.held_no_support;
      continue;
    }
    const double m_bh = q.bh_mass_code;
    // Refinement adequacy (AGN05): injection diameter / N_j vs current width.
    const double width = -st[kMinCellWidth];
    if (m_config.refinement_enabled && width > 0.0) {
      const double required = 2.0 * q.kernel_radius_proper / m_config.refinement_injection_cells;
      if (width > required) {
        report.refinement_requests.push_back(
            AgnRefinementRequest{s.id, s.position, q.support_radius_stored, required, width});
        ++report.unresolved_injection;
      }
    }
    // Thermal.
    if (r.pending_thermal_energy_code > 0.0) {
      if (m_config.thermal_mode == core::AgnThermalModeConfig::kContinuous) {
        q.parameters[kThermalQ] = r.pending_thermal_energy_code;
        q.parameters[kThermalNorm] = st[kSumOmega];
        s.thermal_q = r.pending_thermal_energy_code;
      } else {
        const double capacity = m_config.thermal_delta_t_k * st[kSumWeightMassCv];
        if (capacity > 0.0 && r.pending_thermal_energy_code >= capacity) {
          q.parameters[kThermalQ] = capacity;
          q.parameters[kThermalReservoir] = 1.0;
          q.parameters[kThermalNorm] = st[kSumWeightMassCv];
          s.thermal_q = capacity;
        } else {
          ++report.reservoir_waiting;
        }
      }
    }
    // Kicks with recoil in the same budget.
    if (r.pending_kick_energy_code > 0.0) {
      const Vec3 net{st[kKickNet], st[kKickNet + 1], st[kKickNet + 2]};
      const double b = st[kKickB] - dot(net, s.velocity);
      const double c = st[kKickC] + dot(net, net) / (2.0 * m_bh);
      double a = 0.0;
      if (solveEnergyConsistentAmplitude(c, b, r.pending_kick_energy_code, a) &&
          a * st[kMaxOmegaOverMass] <= v_max) {
        q.parameters[kKickA] = a;
        s.kick_q = r.pending_kick_energy_code;
      } else {
        r.flags |= kBhFlagPendingFeedbackHeld;
        ++report.held_velocity_domain;
      }
    }
    // Jets: independent lobe normalization (equal and opposite momenta).
    if (r.pending_jet_energy_code > 0.0) {
      const double wp = st[kLobePlusOmega];
      const double wm = st[kLobeMinusOmega];
      if (norm(q.axis) == 0.0) {
        r.flags |= kBhFlagPendingFeedbackHeld;
        ++report.held_axis_undefined;
      } else if (!(wp > 0.0) || !(wm > 0.0)) {
        r.flags |= kBhFlagPendingFeedbackHeld;
        ++report.held_lobe_imbalance;
      } else {
        const double b = st[kLobePlusB] / wp + st[kLobeMinusB] / wm;
        const double c = st[kLobePlusC] / (wp * wp) + st[kLobeMinusC] / (wm * wm);
        double a = 0.0;
        if (solveEnergyConsistentAmplitude(c, b, r.pending_jet_energy_code, a) &&
            a * st[kMaxLobeOmegaOverMass] / std::min(wp, wm) <= v_max) {
          q.parameters[kJetA] = a;
          q.parameters[kLobePlusNorm] = wp;
          q.parameters[kLobeMinusNorm] = wm;
          s.jet_q = r.pending_jet_energy_code;
          report.max_lobe_imbalance = std::max(report.max_lobe_imbalance, std::abs(wp - wm) / (wp + wm));
        } else {
          r.flags |= kBhFlagPendingFeedbackHeld;
          ++report.held_velocity_domain;
        }
      }
    }
    // Mass-loaded wind.
    if (r.pending_wind_mass_code > 0.0 && st[kSumWeightMass] > 0.0) {
      const double vw = m_config.wind_velocity_km_s * km_s_code;
      double m_launch = r.pending_wind_mass_code;
      const double e_avail = r.pending_wind_energy_code;
      const Vec3 dir{st[kWindDirection] / st[kSumWeightMass], st[kWindDirection + 1] / st[kSumWeightMass],
                     st[kWindDirection + 2] / st[kSumWeightMass]};
      double m_after = m_bh - m_launch;
      double dk_bh = 0.0;
      // Launchable mass at v_w with the recoil kinetic change charged to the
      // same budget (fixed-point on the second-order recoil term).
      for (int iteration = 0; iteration < 4; ++iteration) {
        const Vec3 recoil{m_launch * vw * dir[0], m_launch * vw * dir[1], m_launch * vw * dir[2]};
        m_after = m_bh - m_launch;
        dk_bh = m_after > 0.0 ? -dot(recoil, s.velocity) + dot(recoil, recoil) / (2.0 * m_after) : 0.0;
        if (!(vw > 0.0) || 0.5 * m_launch * vw * vw + dk_bh <= e_avail) break;
        m_launch = std::max(0.0, 2.0 * (e_avail - std::max(dk_bh, 0.0)) / (vw * vw)) * (1.0 - 1.0e-12);
      }
      const double e_used = e_avail;
      const double thermal_total = e_used - 0.5 * m_launch * vw * vw - dk_bh;
      if (m_launch > 0.0 && m_after > 0.0 && thermal_total >= 0.0) {
        q.parameters[kWindMass] = m_launch;
        q.parameters[kWindVelocity] = vw;
        q.parameters[kWindSpecificThermal] = thermal_total / m_launch;
        q.parameters[kWindNorm] = st[kSumWeightMass];
        const double inv = 1.0 / r.pending_wind_mass_code;
        for (std::size_t el = 0; el < k_chemical_element_count; ++el) {
          q.parameters[kWindElementFraction + el] = r.pending_wind_element_mass_code[el] * inv;
        }
        q.parameters[kWindMetalFraction] = r.pending_wind_metal_mass_code * inv;
        s.wind_q = e_used;
        s.wind_mass = m_launch;
        if (m_launch < r.pending_wind_mass_code) r.flags |= kBhFlagPendingFeedbackHeld;
      } else {
        r.flags |= kBhFlagPendingFeedbackHeld;
      }
    }
  }

  // Round B: commit (prepare -> vote -> commit), then BH-side updates.
  SourceTransactionBatch batch(TransactionFloors{ctx.mass_floor_code, ctx.specific_internal_energy_floor_code});
  std::vector<BhSupportReply> commit_replies;
  {
    const auto routed = exchange.routeQueries(queries);
    std::vector<BhSupportReply> local = evaluate(routed, true, &batch);
    bool ok = true;
    std::string failure;
    if (batch.size() > 0U) {
      ok = batch.prepare(state) && batch.validate();
      if (!ok) failure = batch.report().failure;
    }
    if (!exchange.allAgree(ok)) {
      throw std::runtime_error("AGN injection transaction rejected collectively (no state mutated)" +
                               (failure.empty() ? std::string() : ": " + failure));
    }
    if (batch.size() > 0U) {
      if (!batch.authorize(true)) throw std::logic_error("AGN injection authorization failed after unanimous vote");
      report.touched_gas_rows = batch.commit(state);
      report.gas_mutated = !report.touched_gas_rows.empty();
    }
    commit_replies = exchange.returnReplies(local);
  }
  std::vector<ConservedIncrement> gas_by_source(sources.size());
  std::vector<std::array<double, 4>> shares(sources.size(), {0.0, 0.0, 0.0, 0.0});
  for (const BhSupportReply& r : commit_replies) {
    if (r.query_index >= sources.size() || sources[r.query_index].id != r.bh_id) {
      throw std::runtime_error("AgnFeedbackInjector: commit reply does not match its query");
    }
    gas_by_source[r.query_index] += r.transferred;
    shares[r.query_index][0] += r.accumulators[kThermalShare];
    shares[r.query_index][1] += r.accumulators[kKickShare];
    shares[r.query_index][2] += r.accumulators[kJetShare];
    shares[r.query_index][3] += r.accumulators[kWindShare];
  }
  for (std::size_t i = 0; i < sources.size(); ++i) {
    Pending& s = sources[i];
    BlackHoleEngineRecord& r = s.record;
    const ConservedIncrement& gas = gas_by_source[i];
    const std::size_t p = s.particle_row;
    const double m_old = state.particles.mass_code[p];
    const Vec3 mom_old{m_old * s.velocity[0], m_old * s.velocity[1], m_old * s.velocity[2]};
    const double m_new = m_old - gas.mass_code;
    if (!(m_new > 0.0)) throw std::runtime_error("AgnFeedbackInjector: BH dynamical mass became non-positive");
    const Vec3 v_new{(mom_old[0] - gas.momentum_code[0]) / m_new, (mom_old[1] - gas.momentum_code[1]) / m_new,
                     (mom_old[2] - gas.momentum_code[2]) / m_new};
    const double dk_bh = 0.5 * m_new * dot(v_new, v_new) - 0.5 * m_old * dot(s.velocity, s.velocity);
    state.particles.velocity_x_peculiar[p] = v_new[0];
    state.particles.velocity_y_peculiar[p] = v_new[1];
    state.particles.velocity_z_peculiar[p] = v_new[2];
    // Budget debits use the released amounts.
    const double released = s.thermal_q + s.kick_q + s.jet_q + s.wind_q;
    if (s.thermal_q > 0.0) r.pending_thermal_energy_code = std::max(r.pending_thermal_energy_code - s.thermal_q, 0.0);
    if (s.kick_q > 0.0) r.pending_kick_energy_code = 0.0;
    if (s.jet_q > 0.0) r.pending_jet_energy_code = 0.0;
    if (s.wind_mass > 0.0) {
      const double frac = s.wind_mass / r.pending_wind_mass_code;
      for (std::size_t el = 0; el < k_chemical_element_count; ++el) r.pending_wind_element_mass_code[el] *= (1.0 - frac);
      r.pending_wind_metal_mass_code *= (1.0 - frac);
      r.pending_wind_mass_code -= s.wind_mass;
      r.pending_wind_energy_code = std::max(r.pending_wind_energy_code - s.wind_q, 0.0);
    }
    if (released > 0.0) {
      const double residual = std::abs(gas.total_energy_code + dk_bh - released) / released;
      report.max_energy_closure_residual = std::max(report.max_energy_closure_residual, residual);
    }
    r.cumulative_deposited_energy_code += gas.total_energy_code + dk_bh;
    report.thermal_injected_code += shares[i][0];
    report.kinetic_budget_injected_code += s.kick_q;
    report.jet_budget_injected_code += s.jet_q;
    report.wind_budget_injected_code += s.wind_q;
    report.wind_mass_injected_code += s.wind_mass;
    state.particles.mass_code[p] = r.dynamicalMassCode();
    store.put(s.id, r);
  }
  store.store(state);
  if (!sources.empty()) state.notifyParticleMassMutation();
  return report;
}

}  // namespace cosmosim::physics
