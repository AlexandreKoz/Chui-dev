#include "cosmosim/physics/black_hole/bh_engine.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/physics/black_hole/bh_disk.hpp"
#include "cosmosim/physics/black_hole/bh_dynamics.hpp"
#include "cosmosim/physics/black_hole/bh_merger.hpp"
#include "cosmosim/physics/black_hole/bh_supply.hpp"

namespace cosmosim::physics {
namespace {

using Vec3 = std::array<double, 3>;

[[nodiscard]] Vec3 cross(const Vec3& a, const Vec3& b) noexcept {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
[[nodiscard]] double dot(const Vec3& a, const Vec3& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
[[nodiscard]] double norm(const Vec3& a) noexcept { return std::sqrt(dot(a, a)); }

// Query parameter slots (documented wire contract).
constexpr std::size_t k_param_temperature_split = 0U;
constexpr std::size_t k_param_scale_hot = 1U;
constexpr std::size_t k_param_scale_cold = 2U;
constexpr std::size_t k_param_max_fraction = 3U;
constexpr std::size_t k_param_mass_floor = 4U;
// Reply accumulator slots (withdrawal round).
constexpr std::size_t k_acc_angular_momentum = 0U;  // [0..2] J about the BH (proper r x code momentum)
constexpr std::size_t k_acc_thermal_energy = 3U;
constexpr std::size_t k_acc_kinetic_energy = 4U;
constexpr std::size_t k_acc_requested_mass = 5U;

struct LocalBh {
  std::uint32_t bh_row = 0U;
  std::uint32_t particle_row = 0U;
  std::uint64_t id = 0U;
  Vec3 position{};
  Vec3 velocity{};
  BlackHoleEngineRecord record{};
  // Per-interval working state.
  double kernel_radius_proper = 0.0;
  BhKernelStatistics stats{};
  BhSupplyEstimate supply{};
  double request_hot = 0.0;
  double request_cold = 0.0;
  ConservedIncrement withdrawn{};
  Vec3 withdrawn_angular_momentum{};
  double withdrawn_thermal = 0.0;
  double withdrawn_kinetic = 0.0;
};

[[nodiscard]] Vec3 particlePosition(const core::SimulationState& s, std::size_t p) {
  return {s.particles.position_x_comoving[p], s.particles.position_y_comoving[p], s.particles.position_z_comoving[p]};
}
[[nodiscard]] Vec3 particleVelocity(const core::SimulationState& s, std::size_t p) {
  return {s.particles.velocity_x_peculiar[p], s.particles.velocity_y_peculiar[p], s.particles.velocity_z_peculiar[p]};
}
[[nodiscard]] Vec3 cellPosition(const core::SimulationState& s, std::size_t r) {
  return {s.cells.center_x_comoving[r], s.cells.center_y_comoving[r], s.cells.center_z_comoving[r]};
}
[[nodiscard]] Vec3 cellVelocity(const core::SimulationState& s, std::size_t r) {
  return {s.gas_cells.velocity_x_peculiar[r], s.gas_cells.velocity_y_peculiar[r], s.gas_cells.velocity_z_peculiar[r]};
}
[[nodiscard]] Vec3 properOffset(const ProperFrameContext& f, const Vec3& d) noexcept {
  return {f.properLength(d[0]), f.properLength(d[1]), f.properLength(d[2])};
}

BhSupportQuery makeQuery(const LocalBh& bh, std::size_t slot, std::uint32_t rank, BhSupportRound round) {
  BhSupportQuery q;
  q.bh_id = bh.id;
  q.query_index = slot;
  q.origin_rank = rank;
  q.round = round;
  q.position_stored = bh.position;
  q.velocity_peculiar = bh.velocity;
  q.bh_mass_code = bh.record.horizon_mass_code;
  return q;
}

// --- Round evaluation over this rank's owned leaf rows -------------------

std::vector<BhSupportReply> evaluateRadiusBins(const core::SimulationState& state, const BhIntervalContext& ctx,
                                               const PeriodicBox& box, std::span<const BhSupportQuery> routed,
                                               std::uint32_t rank) {
  std::vector<BhSupportReply> replies(routed.size());
  for (std::size_t i = 0; i < routed.size(); ++i) {
    replies[i].bh_id = routed[i].bh_id;
    replies[i].query_index = routed[i].query_index;
    replies[i].round = routed[i].round;
    replies[i].responder_rank = rank;
  }
  if (routed.empty()) return replies;
  const BhProximityFilter filter(routed, box);
  std::vector<std::uint32_t> hits;
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    if (!ctx.leaves.owned(row)) continue;
    const Vec3 x = cellPosition(state, row);
    filter.matches(x, hits);
    for (const std::uint32_t qi : hits) {
      const Vec3 r = properOffset(ctx.frame, box.delta(routed[qi].position_stored, x));
      const double dist = norm(r);
      for (std::size_t b = 0; b < k_bh_kernel_radius_bins; ++b) {
        if (dist <= routed[qi].radius_bins_proper[b]) {
          replies[qi].bins.enclosed_gas_mass[b] += state.cells.mass_code[row];
          ++replies[qi].bins.enclosed_cells[b];
        }
      }
    }
  }
  return replies;
}

std::vector<BhSupportReply> evaluateKernelStatistics(const core::SimulationState& state, const BhIntervalContext& ctx,
                                                     const PeriodicBox& box, std::span<const BhSupportQuery> routed,
                                                     std::uint32_t rank, bool counter_rotation) {
  std::vector<BhSupportReply> replies(routed.size());
  for (std::size_t i = 0; i < routed.size(); ++i) {
    replies[i].bh_id = routed[i].bh_id;
    replies[i].query_index = routed[i].query_index;
    replies[i].round = routed[i].round;
    replies[i].responder_rank = rank;
  }
  if (routed.empty()) return replies;
  const BhProximityFilter filter(routed, box);
  std::vector<std::uint32_t> hits;
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    if (!ctx.leaves.owned(row)) continue;
    const Vec3 x = cellPosition(state, row);
    filter.matches(x, hits);
    if (hits.empty()) continue;
    std::array<double, 3> hw_stored{};
    if (!bhCellHalfWidthStored(state, static_cast<std::uint32_t>(row), hw_stored)) continue;
    const Vec3 hw = properOffset(ctx.frame, hw_stored);
    const double m = state.cells.mass_code[row];
    const Vec3 v = cellVelocity(state, row);
    const double cs = state.gas_cells.sound_speed_code[row];
    const double temperature = state.gas_cells.temperature_code[row];
    const double volume = 8.0 * hw[0] * hw[1] * hw[2];
    for (const std::uint32_t qi : hits) {
      const BhSupportQuery& q = routed[qi];
      const Vec3 r = properOffset(ctx.frame, box.delta(q.position_stored, x));
      if (norm(r) > q.kernel_radius_proper) continue;
      const bool hot = temperature >= q.parameters[k_param_temperature_split];
      BhKernelStatistics& s = replies[qi].stats;
      if (counter_rotation) {
        if (!hot) {
          const Vec3 dv{v[0] - q.bulk_velocity[0], v[1] - q.bulk_velocity[1], v[2] - q.bulk_velocity[2]};
          if (dot(cross(r, dv), q.axis) < 0.0) s.cold_counter_rotating_mass += m;
        }
        continue;
      }
      const double w = bhCellKernelWeight(r, hw, q.kernel_radius_proper);
      BhPhaseSums& phase = hot ? s.hot : s.cold;
      const Vec3 vrel{v[0] - q.velocity_peculiar[0], v[1] - q.velocity_peculiar[1], v[2] - q.velocity_peculiar[2]};
      const Vec3 j = cross(r, vrel);
      phase.weight_mass += w * m;
      phase.weight_volume += w * volume;
      phase.weight_mass_cs2 += w * m * cs * cs;
      for (std::size_t a = 0; a < 3; ++a) {
        phase.weight_momentum[a] += w * m * v[a];
        phase.weight_angular_momentum[a] += w * m * j[a];
      }
      phase.weight_mass_radius += w * m * norm(r);
      phase.mass += m;
      ++phase.cells;
      if (!hot) {
        const Vec3 rxp = cross(r, v);
        for (std::size_t a = 0; a < 3; ++a) {
          s.cold_mass_position[a] += m * r[a];
          s.cold_momentum[a] += m * v[a];
          s.cold_r_cross_p[a] += m * rxp[a];
        }
      }
      const double width = 2.0 * std::min({hw[0], hw[1], hw[2]});
      s.min_cell_width_proper = s.min_cell_width_proper > 0.0 ? std::min(s.min_cell_width_proper, width) : width;
    }
  }
  // Stars (owned particles only; stars are particles, not mesh).
  for (std::size_t si = 0; si < state.star_particles.size(); ++si) {
    const std::uint32_t p = state.star_particles.particle_index[si];
    if (p >= state.particles.size()) continue;
    const Vec3 x = particlePosition(state, p);
    filter.matches(x, hits);
    if (hits.empty()) continue;
    const double m = state.particles.mass_code[p];
    const Vec3 v = particleVelocity(state, p);
    for (const std::uint32_t qi : hits) {
      const BhSupportQuery& q = routed[qi];
      const Vec3 r = properOffset(ctx.frame, box.delta(q.position_stored, x));
      if (norm(r) > q.kernel_radius_proper) continue;
      BhKernelStatistics& s = replies[qi].stats;
      if (counter_rotation) {
        const Vec3 dv{v[0] - q.bulk_velocity[0], v[1] - q.bulk_velocity[1], v[2] - q.bulk_velocity[2]};
        if (dot(cross(r, dv), q.axis) < 0.0) s.star_counter_rotating_mass += m;
        continue;
      }
      const Vec3 rxp = cross(r, v);
      s.star_mass += m;
      for (std::size_t a = 0; a < 3; ++a) {
        s.star_mass_position[a] += m * r[a];
        s.star_momentum[a] += m * v[a];
        s.star_r_cross_p[a] += m * rxp[a];
      }
      ++s.stars;
    }
  }
  return replies;
}

// Withdrawal: joint limitation of overlapping requests per row, one plan per
// routed query, local prepare/validate, collective vote, commit.
std::vector<BhSupportReply> evaluateWithdrawal(core::SimulationState& state, const BhIntervalContext& ctx,
                                               const PeriodicBox& box, std::span<const BhSupportQuery> routed,
                                               BhSupportExchange& exchange, const TransactionFloors& floors,
                                               BhEngineIntervalReport& report) {
  const std::uint32_t rank = exchange.worldRank();
  std::vector<BhSupportReply> replies(routed.size());
  for (std::size_t i = 0; i < routed.size(); ++i) {
    replies[i].bh_id = routed[i].bh_id;
    replies[i].query_index = routed[i].query_index;
    replies[i].round = routed[i].round;
    replies[i].responder_rank = rank;
  }
  struct Request {
    std::uint32_t row;
    std::uint32_t query;
    double mass;
  };
  std::vector<Request> requests;
  if (!routed.empty()) {
    const BhProximityFilter filter(routed, box);
    std::vector<std::uint32_t> hits;
    for (std::size_t row = 0; row < state.cells.size(); ++row) {
      if (!ctx.leaves.owned(row)) continue;
      const Vec3 x = cellPosition(state, row);
      filter.matches(x, hits);
      if (hits.empty()) continue;
      std::array<double, 3> hw_stored{};
      if (!bhCellHalfWidthStored(state, static_cast<std::uint32_t>(row), hw_stored)) continue;
      const Vec3 hw = properOffset(ctx.frame, hw_stored);
      const double m = state.cells.mass_code[row];
      const double temperature = state.gas_cells.temperature_code[row];
      for (const std::uint32_t qi : hits) {
        const BhSupportQuery& q = routed[qi];
        const Vec3 r = properOffset(ctx.frame, box.delta(q.position_stored, x));
        if (norm(r) > q.kernel_radius_proper) continue;
        const bool hot = temperature >= q.parameters[k_param_temperature_split];
        const double scale = hot ? q.parameters[k_param_scale_hot] : q.parameters[k_param_scale_cold];
        if (!(scale > 0.0)) continue;
        const double w = bhCellKernelWeight(r, hw, q.kernel_radius_proper);
        const double req = scale * w * m;
        if (req > 0.0) requests.push_back(Request{static_cast<std::uint32_t>(row), qi, req});
      }
    }
  }
  // Joint limitation per row (deterministic: rows sorted, then query order).
  std::sort(requests.begin(), requests.end(), [](const Request& a, const Request& b) {
    return a.row != b.row ? a.row < b.row : a.query < b.query;
  });
  std::unordered_map<std::uint32_t, double> row_scale;
  for (std::size_t i = 0; i < requests.size();) {
    std::size_t j = i;
    double total = 0.0;
    while (j < requests.size() && requests[j].row == requests[i].row) total += requests[j++].mass;
    const std::uint32_t row = requests[i].row;
    const double m = state.cells.mass_code[row];
    const double max_fraction = routed[requests[i].query].parameters[k_param_max_fraction];
    const double available = std::max(0.0, std::min(max_fraction * m, m - floors.mass_floor_code));
    row_scale[row] = total > available ? (total > 0.0 ? available / total : 0.0) : 1.0;
    i = j;
  }
  // Plans (one per routed query).
  SourceTransactionBatch batch(floors);
  std::vector<SourceTransactionPlan> plans(routed.size());
  const bool composition = state.gas_composition.enabled();
  std::optional<CompositionLayout> layout;
  if (composition) layout = CompositionLayout::fromDescriptor(state.gas_composition.descriptor());
  for (std::size_t qi = 0; qi < routed.size(); ++qi) {
    plans[qi].id = SourceTransactionId{routed[qi].bh_id, ctx.interval_index, SourceChannel::kBhAccretionWithdrawal,
                                       static_cast<std::uint32_t>(rank)};
    plans[qi].topology_epoch = state.cellIndexGeneration();
    plans[qi].owner_rank = routed[qi].origin_rank;
    plans[qi].evaluation_scale_factor = ctx.frame.scale_factor;
    plans[qi].intended_time_code = ctx.time_end_code;
  }
  for (const Request& rq : requests) {
    const std::uint32_t row = rq.row;
    const double m = state.cells.mass_code[row];
    const double actual = rq.mass * row_scale[row];
    if (!(actual > 0.0)) {
      replies[rq.query].accumulators[k_acc_requested_mass] += rq.mass;
      continue;
    }
    const double f = actual / m;
    const Vec3 v = cellVelocity(state, row);
    const double u = state.gas_cells.internal_energy_code[row];
    GasCellIncrement inc;
    inc.gas_cell_id = state.gas_cells.gas_cell_id[row];
    inc.local_row = row;
    inc.withdrawal_fraction = f;
    inc.delta.mass_code = -actual;
    for (std::size_t a = 0; a < 3; ++a) inc.delta.momentum_code[a] = -actual * v[a];
    const double kinetic = 0.5 * actual * dot(v, v);
    const double thermal = actual * u;
    inc.delta.total_energy_code = -(thermal + kinetic);
    if (composition) {
      for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
        inc.delta.element_mass_code[e] =
            -f * state.gas_composition.value(layout->elementLane(static_cast<ChemicalElement>(e)), row);
      }
    }
    inc.delta.metal_mass_code = -f * state.gas_cells.metal_mass_code[row];
    const Vec3 r = properOffset(ctx.frame, box.delta(routed[rq.query].position_stored, cellPosition(state, row)));
    const Vec3 vrel{v[0] - routed[rq.query].velocity_peculiar[0], v[1] - routed[rq.query].velocity_peculiar[1],
                    v[2] - routed[rq.query].velocity_peculiar[2]};
    const Vec3 j = cross(r, vrel);
    for (std::size_t a = 0; a < 3; ++a) {
      inc.delta.angular_momentum_code[a] = -actual * j[a];
      replies[rq.query].accumulators[k_acc_angular_momentum + a] += actual * j[a];
    }
    replies[rq.query].accumulators[k_acc_thermal_energy] += thermal;
    replies[rq.query].accumulators[k_acc_kinetic_energy] += kinetic;
    replies[rq.query].accumulators[k_acc_requested_mass] += rq.mass;
    ++replies[rq.query].rows;
    plans[rq.query].cells.push_back(inc);
  }
  for (std::size_t qi = 0; qi < routed.size(); ++qi) {
    ConservedIncrement gas{};
    for (const GasCellIncrement& c : plans[qi].cells) gas += c.delta;
    plans[qi].source_increment = gas.negated();
    replies[qi].transferred = plans[qi].source_increment;  // positive: received by the BH
    if (!plans[qi].cells.empty()) batch.add(std::move(plans[qi]));
  }
  bool ok = true;
  std::string failure;
  if (batch.size() > 0U) {
    ok = batch.prepare(state) && batch.validate();
    if (!ok) failure = batch.report().failure;
  }
  const bool agreed = exchange.allAgree(ok);
  if (!agreed) {
    throw std::runtime_error("BH withdrawal transaction rejected collectively (no state mutated)" +
                             (failure.empty() ? std::string() : ": " + failure));
  }
  if (batch.size() > 0U) {
    if (!batch.authorize(true)) throw std::logic_error("BH withdrawal authorization failed after a unanimous vote");
    const std::vector<std::uint32_t> touched = batch.commit(state);
    report.touched_gas_rows.insert(report.touched_gas_rows.end(), touched.begin(), touched.end());
    report.gas_mutated = report.gas_mutated || !touched.empty();
    report.max_transaction_residual =
        std::max(report.max_transaction_residual, batch.report().max_relative_closure_residual);
  }
  return replies;
}

[[nodiscard]] std::vector<LocalBh> collectLocalBlackHoles(const core::SimulationState& state,
                                                          const BlackHoleEngineStore& store) {
  std::vector<LocalBh> out;
  out.reserve(state.black_holes.size());
  for (std::size_t b = 0; b < state.black_holes.size(); ++b) {
    LocalBh bh;
    bh.bh_row = static_cast<std::uint32_t>(b);
    bh.particle_row = state.black_holes.particle_index[b];
    if (bh.particle_row >= state.particles.size()) {
      throw std::runtime_error("BlackHoleEngine: BH sidecar row references a missing particle");
    }
    bh.id = state.particle_sidecar.particle_id[bh.particle_row];
    bh.position = particlePosition(state, bh.particle_row);
    bh.velocity = particleVelocity(state, bh.particle_row);
    bh.record = store.get(bh.id);
    out.push_back(bh);
  }
  std::sort(out.begin(), out.end(), [](const LocalBh& a, const LocalBh& b) { return a.id < b.id; });
  return out;
}

void mergeReplies(std::vector<LocalBh>& bhs, std::span<const BhSupportReply> replies, BhSupportRound round,
                  std::vector<BhRadiusBinSums>* bins) {
  for (const BhSupportReply& r : replies) {
    if (r.query_index >= bhs.size() || bhs[r.query_index].id != r.bh_id || r.round != round) {
      throw std::runtime_error("BlackHoleEngine: support reply does not match its query (protocol violation)");
    }
    LocalBh& bh = bhs[r.query_index];
    switch (round) {
      case BhSupportRound::kRadiusBins:
        (*bins)[r.query_index].merge(r.bins);
        break;
      case BhSupportRound::kKernelStatistics:
      case BhSupportRound::kCounterRotation:
        bh.stats.merge(r.stats);
        break;
      case BhSupportRound::kWithdrawal:
        bh.withdrawn += r.transferred;
        for (std::size_t a = 0; a < 3; ++a) bh.withdrawn_angular_momentum[a] += r.accumulators[k_acc_angular_momentum + a];
        bh.withdrawn_thermal += r.accumulators[k_acc_thermal_energy];
        bh.withdrawn_kinetic += r.accumulators[k_acc_kinetic_energy];
        break;
      default:
        throw std::logic_error("BlackHoleEngine: unexpected support round");
    }
  }
}

}  // namespace

BlackHoleEngine::BlackHoleEngine(BlackHoleEngineParameters p) : m_parameters(std::move(p)) {
  if (m_parameters.model.model != core::BhEngineModelConfig::kReferenceV1) {
    throw std::invalid_argument("BlackHoleEngine requires physics.bh_engine_model=reference_v1");
  }
  if (!(m_parameters.kernel_radius_min_proper_code > 0.0) ||
      !(m_parameters.kernel_radius_max_proper_code >= m_parameters.kernel_radius_min_proper_code)) {
    throw std::invalid_argument("BlackHoleEngine: invalid kernel radius range");
  }
  if (m_parameters.model.disk_model == core::BhDiskModelConfig::kScalarReservoir &&
      !(m_parameters.viscous_time_code > 0.0)) {
    throw std::invalid_argument("BlackHoleEngine: scalar reservoir requires a positive viscous time");
  }
  if (!(m_parameters.efficiencies.total() < 1.0)) {
    throw std::invalid_argument("BlackHoleEngine: channel efficiencies must sum below 1");
  }
}

BhEngineIntervalReport BlackHoleEngine::advanceInterval(core::SimulationState& state, const BhIntervalContext& ctx,
                                                        BhSupportExchange& exchange) const {
  ctx.frame.validate();
  if (!(ctx.dt_proper_code > 0.0) || !std::isfinite(ctx.dt_proper_code)) {
    throw std::invalid_argument("BlackHoleEngine: interval requires a finite positive proper dt");
  }
  const BlackHoleEngineParameters& P = m_parameters;
  const core::BlackHoleEngineConfig& M = P.model;
  BhEngineIntervalReport report;
  BlackHoleEngineStore store(state);
  (void)store.adoptLegacyBlackHoles(state);
  std::vector<LocalBh> bhs = collectLocalBlackHoles(state, store);
  report.local_black_holes = bhs.size();
  const std::uint32_t rank = exchange.worldRank();
  const double dt = ctx.dt_proper_code;
  const EddingtonConstantsCode edd =
      EddingtonConstantsCode::fromFrame(ctx.frame, P.proton_mass_si, P.thomson_cross_section_si);
  const auto radii = bhKernelRadiusCandidates(P.kernel_radius_min_proper_code, P.kernel_radius_max_proper_code);
  const bool use_torque = M.supply_model == core::BhSupplyModelConfig::kBondiTorque ||
      M.supply_model == core::BhSupplyModelConfig::kTorque;
  const bool use_bondi = M.supply_model != core::BhSupplyModelConfig::kTorque;
  const double support_stored = ctx.frame.storedLength(P.kernel_radius_max_proper_code);

  // Round 1: kernel radius by enclosed gas mass.
  std::vector<BhSupportQuery> queries;
  queries.reserve(bhs.size());
  for (std::size_t i = 0; i < bhs.size(); ++i) {
    BhSupportQuery q = makeQuery(bhs[i], i, rank, BhSupportRound::kRadiusBins);
    q.support_radius_stored = support_stored;
    q.radius_bins_proper = radii;
    queries.push_back(q);
  }
  std::vector<BhRadiusBinSums> bins(bhs.size());
  {
    const auto routed = exchange.routeQueries(queries);
    const auto replies = exchange.returnReplies(evaluateRadiusBins(state, ctx, P.box, routed, rank));
    mergeReplies(bhs, replies, BhSupportRound::kRadiusBins, &bins);
  }
  // Round 2: kernel statistics at the selected radius.
  for (std::size_t i = 0; i < bhs.size(); ++i) {
    const std::size_t b = bhSelectKernelRadius(bins[i], M.kernel_target_gas_mass_code);
    bhs[i].kernel_radius_proper = radii[b];
    BhSupportQuery& q = queries[i];
    q.round = BhSupportRound::kKernelStatistics;
    q.kernel_radius_proper = radii[b];
    q.support_radius_stored = ctx.frame.storedLength(radii[b]);
    q.parameters[k_param_temperature_split] = M.phase_temperature_split_k;
  }
  {
    const auto routed = exchange.routeQueries(queries);
    const auto replies = exchange.returnReplies(evaluateKernelStatistics(state, ctx, P.box, routed, rank, false));
    mergeReplies(bhs, replies, BhSupportRound::kKernelStatistics, nullptr);
  }
  // Round 3: counter-rotating masses about the measured baryonic axis.
  if (use_torque) {
    for (std::size_t i = 0; i < bhs.size(); ++i) {
      BhSupportQuery& q = queries[i];
      q.round = BhSupportRound::kCounterRotation;
      q.axis = bhDiskAxis(bhs[i].stats);
      const double m = bhs[i].stats.star_mass + bhs[i].stats.cold.mass;
      for (std::size_t a = 0; a < 3; ++a) {
        q.bulk_velocity[a] = m > 0.0 ? (bhs[i].stats.star_momentum[a] + bhs[i].stats.cold_momentum[a]) / m : 0.0;
      }
    }
    const auto routed = exchange.routeQueries(queries);
    const auto replies = exchange.returnReplies(evaluateKernelStatistics(state, ctx, P.box, routed, rank, true));
    mergeReplies(bhs, replies, BhSupportRound::kCounterRotation, nullptr);
  }
  // Supply estimate and requested masses.
  for (std::size_t i = 0; i < bhs.size(); ++i) {
    LocalBh& bh = bhs[i];
    BhSupplyInputs in;
    in.bh_mass_code = bh.record.horizon_mass_code;
    in.bh_velocity_peculiar = bh.velocity;
    in.kernel_radius_proper = bh.kernel_radius_proper;
    in.alpha = P.alpha_bondi;
    in.c_visc = M.supply_model == core::BhSupplyModelConfig::kBondiAngularSuppressed ? M.angular_suppression_c_visc : 0.0;
    in.use_bondi = use_bondi;
    in.use_torque = use_torque;
    in.torque_epsilon = M.torque_epsilon;
    in.kernel_min_cells = M.kernel_min_cells;
    in.torque_min_samples = M.torque_min_samples;
    in.softening_proper = ctx.frame.properLength(P.bh_softening_stored_code);
    bh.supply = estimateBhSupply(bh.stats, in, ctx.frame);
    bh.record.flags &= ~(kBhFlagKernelInadequate | kBhFlagTorqueInadequate | kBhFlagEddingtonCapped |
                         kBhFlagSupportTruncated | kBhFlagBondiUnresolved);
    if (!bh.supply.kernel_adequate) {
      bh.record.flags |= kBhFlagKernelInadequate;
      ++report.kernel_inadequate;
    }
    if (use_torque && !bh.supply.torque_adequate) {
      bh.record.flags |= kBhFlagTorqueInadequate;
      ++report.torque_inadequate;
    }
    if (bh.supply.bondi_radius_over_dx < 1.0) bh.record.flags |= kBhFlagBondiUnresolved;
    bh.request_hot = bh.supply.hot_rate_code * dt;
    bh.request_cold = bh.supply.cold_rate_code * dt;
    if (M.disk_model == core::BhDiskModelConfig::kInstantaneous &&
        M.eddington_cap == core::BhEddingtonCapConfig::kReferenceRate) {
      // Legacy-equivalent comparator: cap the REQUEST (no disk to hold excess).
      const double cap =
          eddingtonReferenceRateCode(bh.record.horizon_mass_code, M.eddington_reference_efficiency, edd) * dt;
      const double total = bh.request_hot + bh.request_cold;
      if (total > cap && total > 0.0) {
        bh.request_hot *= cap / total;
        bh.request_cold *= cap / total;
        bh.record.flags |= kBhFlagEddingtonCapped;
        ++report.eddington_capped;
      }
    }
    BhSupportQuery& q = queries[i];
    q.round = BhSupportRound::kWithdrawal;
    q.parameters[k_param_scale_hot] =
        bh.stats.hot.weight_mass > 0.0 ? bh.request_hot / bh.stats.hot.weight_mass : 0.0;
    q.parameters[k_param_scale_cold] =
        bh.stats.cold.weight_mass > 0.0 ? bh.request_cold / bh.stats.cold.weight_mass : 0.0;
    q.parameters[k_param_max_fraction] = P.max_withdrawal_fraction;
    q.parameters[k_param_mass_floor] = P.mass_floor_code;
  }
  // Round 4: withdrawal transaction (prepare -> vote -> commit).
  {
    const TransactionFloors floors{P.mass_floor_code, P.specific_internal_energy_floor_code};
    const auto routed = exchange.routeQueries(queries);
    const auto replies = exchange.returnReplies(evaluateWithdrawal(state, ctx, P.box, routed, exchange, floors, report));
    mergeReplies(bhs, replies, BhSupportRound::kWithdrawal, nullptr);
  }
  // Owner update: disk, horizon ledger, budgets, particle mass/momentum.
  const BhMassLedger ledger = M.mass_ledger == core::BhMassLedgerConfig::kRelativisticDefect
      ? BhMassLedger::kRelativisticDefect
      : BhMassLedger::kLegacyNewtonian;
  for (LocalBh& bh : bhs) {
    BlackHoleEngineRecord& rec = bh.record;
    const double requested = bh.request_hot + bh.request_cold;
    const double withdrawn = bh.withdrawn.mass_code;
    rec.cumulative_supply_requested_mass_code += requested;
    rec.cumulative_withdrawn_mass_code += withdrawn;
    rec.cumulative_unmet_supply_mass_code += std::max(requested - withdrawn, 0.0);
    if (requested > 0.0 && withdrawn < requested * (1.0 - 1.0e-12)) rec.flags |= kBhFlagSupportTruncated;
    rec.cumulative_captured_gas_energy_code += bh.withdrawn.total_energy_code;
    // Disk composition and captured orbital angular momentum.
    rec.disk_mass_code += withdrawn;
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) rec.disk_element_mass_code[e] += bh.withdrawn.element_mass_code[e];
    rec.disk_metal_mass_code += bh.withdrawn.metal_mass_code;
    for (std::size_t a = 0; a < 3; ++a) rec.disk_angular_momentum_code[a] += bh.withdrawn_angular_momentum[a];
    const double disk_before = rec.disk_mass_code;  // after supply credited
    // Integrate disk over dt (substeps bound the relative BH-mass change).
    double inner = 0.0;
    double wind = 0.0;
    double drained = 0.0;
    double horizon = rec.horizon_mass_code;
    double defect = 0.0;
    BhEmissionBudget total{};
    total.bh_particle_id = bh.id;
    total.interval_index = ctx.interval_index;
    total.interval_begin_code = ctx.time_begin_code;
    total.interval_end_code = ctx.time_end_code;
    if (M.disk_model == core::BhDiskModelConfig::kScalarReservoir) {
      // Supply S enters uniformly over the interval: undo the instantaneous
      // credit above and integrate with S = withdrawn/dt.
      double m_disk = disk_before - withdrawn;
      const double s_rate = withdrawn / dt;
      const double estimate = std::max(m_disk + withdrawn, 0.0);
      std::uint32_t n = 1U;
      if (horizon > 0.0 && estimate > 0.0) {
        const double ratio = estimate / (M.disk_max_relative_bh_mass_change * horizon);
        n = static_cast<std::uint32_t>(std::clamp(std::ceil(ratio), 1.0, static_cast<double>(P.max_disk_substeps)));
      }
      const double h = dt / n;
      for (std::uint32_t k = 0; k < n; ++k) {
        BhDiskStepInputs din;
        din.disk_mass_code = m_disk;
        din.supply_rate_code = s_rate;
        din.viscous_time_code = P.viscous_time_code;
        din.wind_loading = M.disk_wind_loading;
        din.dt_code = h;
        if (M.eddington_cap == core::BhEddingtonCapConfig::kReferenceRate) {
          din.inner_rate_cap_code = eddingtonReferenceRateCode(horizon, M.eddington_reference_efficiency, edd);
        }
        const BhDiskStepResult dres = integrateScalarDisk(din);
        if (dres.capped_time_code > 0.0) rec.flags |= kBhFlagEddingtonCapped;
        const BhEmissionBudget sub = computeBhEmissionBudget(ledger, dres.inner_mass_code, h, horizon, P.efficiencies,
                                                             M.eddington_reference_efficiency, edd);
        horizon += sub.horizon_mass_increment_code;
        defect += sub.mass_defect_code;
        total.radiative_energy_code += sub.radiative_energy_code;
        total.thermal_energy_code += sub.thermal_energy_code;
        total.kick_energy_code += sub.kick_energy_code;
        total.wind_energy_code += sub.wind_energy_code;
        total.jet_energy_code += sub.jet_energy_code;
        inner += dres.inner_mass_code;
        wind += dres.wind_mass_code;
        drained += dres.drained_mass_code;
        m_disk = dres.disk_mass_end_code;
      }
      report.disk_substeps += n;
      if ((rec.flags & kBhFlagEddingtonCapped) != 0U) ++report.eddington_capped;
      rec.disk_mass_code = m_disk;
    } else {
      inner = rec.disk_mass_code;  // instantaneous: everything withdrawn reaches the engine now
      drained = inner;
      const BhEmissionBudget sub = computeBhEmissionBudget(ledger, inner, dt, horizon, P.efficiencies,
                                                           M.eddington_reference_efficiency, edd);
      horizon += sub.horizon_mass_increment_code;
      defect += sub.mass_defect_code;
      total.radiative_energy_code = sub.radiative_energy_code;
      total.thermal_energy_code = sub.thermal_energy_code;
      total.kick_energy_code = sub.kick_energy_code;
      total.wind_energy_code = sub.wind_energy_code;
      total.jet_energy_code = sub.jet_energy_code;
      rec.disk_mass_code = 0.0;
    }
    // Disk composition drains proportionally (inner engine + wind).
    const double drain_fraction = disk_before > 0.0 ? std::clamp(drained / disk_before, 0.0, 1.0) : 0.0;
    const double wind_fraction = drained > 0.0 ? wind / drained : 0.0;
    for (std::size_t e = 0; e < k_chemical_element_count; ++e) {
      const double out = rec.disk_element_mass_code[e] * drain_fraction;
      rec.disk_element_mass_code[e] -= out;
      rec.pending_wind_element_mass_code[e] += out * wind_fraction;
    }
    {
      const double out = rec.disk_metal_mass_code * drain_fraction;
      rec.disk_metal_mass_code -= out;
      rec.pending_wind_metal_mass_code += out * wind_fraction;
    }
    // Angular momentum of drained material leaves the disk vector (P0 keeps
    // no spin: the inner part is recorded with the horizon mass, the wind part
    // leaves with the ejecta; vector J evolution is P1).
    for (std::size_t a = 0; a < 3; ++a) rec.disk_angular_momentum_code[a] *= (1.0 - drain_fraction);
    rec.pending_wind_mass_code += wind;
    rec.cumulative_inner_inflow_mass_code += inner;
    rec.cumulative_disk_wind_mass_code += wind;
    rec.cumulative_mass_defect_code += defect;
    rec.horizon_mass_code = horizon;
    rec.pending_radiation_energy_code += total.radiative_energy_code;
    rec.pending_thermal_energy_code += total.thermal_energy_code;
    rec.pending_kick_energy_code += total.kick_energy_code;
    rec.pending_wind_energy_code += total.wind_energy_code;
    rec.pending_jet_energy_code += total.jet_energy_code;
    rec.cumulative_radiated_energy_code += total.radiative_energy_code;
    const double mechanical = total.thermal_energy_code + total.kick_energy_code + total.wind_energy_code +
        total.jet_energy_code;
    rec.cumulative_mechanical_energy_code += mechanical;
    total.inner_mass_code = inner;
    total.horizon_mass_increment_code = horizon - bh.record.horizon_mass_code;
    total.mass_defect_code = defect;
    total.luminosity_code = total.radiative_energy_code / dt;
    const double l_edd = eddingtonLuminosityCode(horizon, edd);
    total.luminosity_eddington_ratio = l_edd > 0.0 ? total.luminosity_code / l_edd : 0.0;
    const double mdot_ref = eddingtonReferenceRateCode(horizon, M.eddington_reference_efficiency, edd);
    total.rate_eddington_reference_ratio = mdot_ref > 0.0 ? (inner / dt) / mdot_ref : 0.0;
    rec.last_supply_rate_code = withdrawn / dt;
    rec.last_inner_rate_code = inner / dt;
    rec.last_eddington_rate_ratio = total.rate_eddington_reference_ratio;
    rec.last_luminosity_eddington_ratio = total.luminosity_eddington_ratio;
    rec.last_kernel_radius_proper_code = bh.kernel_radius_proper;
    rec.last_bondi_radius_over_dx = bh.supply.bondi_radius_over_dx;
    rec.last_hot_gas_mass_code = bh.stats.hot.mass;
    rec.last_cold_gas_mass_code = bh.stats.cold.mass;
    rec.source_time_code = ctx.time_end_code;
    rec.interval_index = ctx.interval_index;
    ++rec.event_sequence;
    // Particle: momentum conservation of the capture; the mass defect leaves
    // isotropically in the BH frame (velocity unchanged by the defect).
    const std::size_t p = bh.particle_row;
    const double m_old = state.particles.mass_code[p];
    const double m_capture = m_old + withdrawn;
    Vec3 v_new = bh.velocity;
    if (m_capture > 0.0) {
      for (std::size_t a = 0; a < 3; ++a) v_new[a] = (m_old * bh.velocity[a] + bh.withdrawn.momentum_code[a]) / m_capture;
    }
    const double ke_before = 0.5 * m_old * dot(bh.velocity, bh.velocity) + bh.withdrawn_kinetic;
    const double ke_after = 0.5 * m_capture * dot(v_new, v_new);
    rec.cumulative_capture_dissipation_code += ke_before - ke_after;  // >= 0: inelastic capture
    state.particles.velocity_x_peculiar[p] = v_new[0];
    state.particles.velocity_y_peculiar[p] = v_new[1];
    state.particles.velocity_z_peculiar[p] = v_new[2];
    state.particles.mass_code[p] = rec.dynamicalMassCode();
    state.black_holes.subgrid_mass_code[bh.bh_row] = rec.horizon_mass_code;
    state.black_holes.accretion_rate_code[bh.bh_row] = inner / dt;
    state.black_holes.eddington_ratio[bh.bh_row] = total.rate_eddington_reference_ratio;
    state.black_holes.feedback_energy_code[bh.bh_row] = mechanical + total.radiative_energy_code;
    state.black_holes.cumulative_accreted_mass_code[bh.bh_row] += inner;
    state.black_holes.cumulative_feedback_energy_code[bh.bh_row] += mechanical + total.radiative_energy_code;
    state.black_holes.duty_cycle_total_time_code[bh.bh_row] += ctx.time_end_code - ctx.time_begin_code;
    report.particle_mass_mutated = true;
    report.requested_mass_code += requested;
    report.withdrawn_mass_code += withdrawn;
    report.unmet_mass_code += std::max(requested - withdrawn, 0.0);
    report.inner_mass_code += inner;
    report.disk_wind_mass_code += wind;
    report.mass_defect_code += defect;
    report.radiated_energy_code += total.radiative_energy_code;
    report.mechanical_energy_code += mechanical;
    report.captured_gas_energy_code += bh.withdrawn.total_energy_code;
    report.budgets.push_back(total);
    bh.velocity = v_new;
  }

  // Dynamical friction (local DM/star background; reaction on the sampled
  // particles conserves momentum exactly; the energy residual is a named
  // unresolved reservoir).
  if (M.dynamics_model == core::BhDynamicsModelConfig::kUnresolvedChandrasekhar && !bhs.empty()) {
    const double eps_proper = ctx.frame.properLength(P.bh_softening_stored_code);
    struct Sample {
      std::uint32_t row;
      double r;
    };
    // One bounded pass over local particles: candidate samples per BH from
    // the proximity filter (no N_BH x N_particle scan).
    std::vector<BhSupportQuery> df_queries;
    df_queries.reserve(bhs.size());
    for (std::size_t i = 0; i < bhs.size(); ++i) {
      BhSupportQuery q = makeQuery(bhs[i], i, rank, BhSupportRound::kKernelStatistics);
      q.support_radius_stored = 16.0 * P.bh_softening_stored_code;
      df_queries.push_back(q);
    }
    std::vector<std::vector<Sample>> dm_by_bh(bhs.size());
    std::vector<std::vector<Sample>> star_by_bh(bhs.size());
    {
      const BhProximityFilter df_filter(df_queries, P.box);
      std::vector<std::uint32_t> hits;
      for (std::size_t p = 0; p < state.particles.size(); ++p) {
        const auto species = static_cast<core::ParticleSpecies>(state.particleSpeciesTag(p));
        if (species != core::ParticleSpecies::kDarkMatter && species != core::ParticleSpecies::kStar) continue;
        const Vec3 x = particlePosition(state, p);
        df_filter.matches(x, hits);
        for (const std::uint32_t qi : hits) {
          const double r = norm(properOffset(ctx.frame, P.box.delta(bhs[qi].position, x)));
          (species == core::ParticleSpecies::kDarkMatter ? dm_by_bh[qi] : star_by_bh[qi])
              .push_back(Sample{static_cast<std::uint32_t>(p), r});
        }
      }
    }
    for (std::size_t bi = 0; bi < bhs.size(); ++bi) {
      LocalBh& bh = bhs[bi];
      std::vector<Sample>& dm_samples = dm_by_bh[bi];
      std::vector<Sample>& star_samples = star_by_bh[bi];
      auto background = [&](std::vector<Sample>& samples, std::vector<std::uint32_t>& used) {
        BhFrictionBackground bg;
        std::sort(samples.begin(), samples.end(), [](const Sample& a, const Sample& b) {
          return a.r != b.r ? a.r < b.r : a.row < b.row;
        });
        const std::size_t k = std::min<std::size_t>(samples.size(), M.df_neighbor_count);
        if (k < M.df_neighbor_count || k == 0U) return bg;  // inadequate sampling: no correction
        double mass = 0.0;
        Vec3 mom{};
        for (std::size_t i = 0; i < k; ++i) {
          const std::uint32_t p = samples[i].row;
          const double m = state.particles.mass_code[p];
          const Vec3 v = particleVelocity(state, p);
          mass += m;
          for (std::size_t a = 0; a < 3; ++a) mom[a] += m * v[a];
          used.push_back(p);
        }
        const double r_k = samples[k - 1U].r;
        if (!(mass > 0.0) || !(r_k > 0.0)) return bg;
        for (std::size_t a = 0; a < 3; ++a) bg.bulk_velocity[a] = mom[a] / mass;
        double s2 = 0.0;
        for (std::size_t i = 0; i < k; ++i) {
          const std::uint32_t p = samples[i].row;
          const Vec3 v = particleVelocity(state, p);
          const Vec3 dv{v[0] - bg.bulk_velocity[0], v[1] - bg.bulk_velocity[1], v[2] - bg.bulk_velocity[2]};
          s2 += state.particles.mass_code[p] * dot(dv, dv);
        }
        bg.sigma_1d = std::sqrt(s2 / (3.0 * mass));
        bg.density_proper = mass / (4.0 / 3.0 * core::constants::k_pi * r_k * r_k * r_k);
        bg.samples = k;
        return bg;
      };
      std::vector<std::uint32_t> used;
      const BhFrictionBackground dm = background(dm_samples, used);
      const BhFrictionBackground st = background(star_samples, used);
      BhFrictionInputs fin;
      fin.bh_mass_code = bh.record.horizon_mass_code;
      fin.newton_g_code = ctx.frame.newtonGCode();
      fin.softening_proper = eps_proper;
      fin.timestep_eta = M.df_timestep_eta;
      const BhFrictionResult fr = integrateDynamicalFriction(bh.velocity, dm, st, fin, dt);
      if (!fr.applied || used.empty()) continue;
      const std::size_t p = bh.particle_row;
      const double m_bh = state.particles.mass_code[p];
      const Vec3 dp{m_bh * fr.velocity_change[0], m_bh * fr.velocity_change[1], m_bh * fr.velocity_change[2]};
      double m_bg = 0.0;
      for (const std::uint32_t q : used) m_bg += state.particles.mass_code[q];
      if (!(m_bg > 0.0)) continue;
      double dke = 0.5 * m_bh * (dot(Vec3{bh.velocity[0] + fr.velocity_change[0], bh.velocity[1] + fr.velocity_change[1],
                                          bh.velocity[2] + fr.velocity_change[2]},
                                     Vec3{bh.velocity[0] + fr.velocity_change[0], bh.velocity[1] + fr.velocity_change[1],
                                          bh.velocity[2] + fr.velocity_change[2]}) -
                                 dot(bh.velocity, bh.velocity));
      const Vec3 dv_bg{-dp[0] / m_bg, -dp[1] / m_bg, -dp[2] / m_bg};
      for (const std::uint32_t q : used) {
        const Vec3 v = particleVelocity(state, q);
        const Vec3 v2{v[0] + dv_bg[0], v[1] + dv_bg[1], v[2] + dv_bg[2]};
        dke += 0.5 * state.particles.mass_code[q] * (dot(v2, v2) - dot(v, v));
        state.particles.velocity_x_peculiar[q] = v2[0];
        state.particles.velocity_y_peculiar[q] = v2[1];
        state.particles.velocity_z_peculiar[q] = v2[2];
      }
      for (std::size_t a = 0; a < 3; ++a) bh.velocity[a] += fr.velocity_change[a];
      state.particles.velocity_x_peculiar[p] = bh.velocity[0];
      state.particles.velocity_y_peculiar[p] = bh.velocity[1];
      state.particles.velocity_z_peculiar[p] = bh.velocity[2];
      bh.record.cumulative_df_reservoir_energy_code += -dke;  // dissipated into the unresolved background
      report.df_reservoir_energy_code += -dke;
      ++report.df_applied;
    }
  }

  for (const LocalBh& bh : bhs) store.put(bh.id, bh.record);

  // Bound-pair mergers among local BHs (cross-rank pairs are not evaluated
  // in pass 01; recorded gap).
  if (M.merger_model == core::BhMergerModelConfig::kBoundPair && bhs.size() > 1U) {
    std::vector<BhMergerCandidate> cand;
    cand.reserve(bhs.size());
    for (const LocalBh& bh : bhs) {
      BhMergerCandidate c;
      c.particle_id = bh.id;
      c.local_bh_row = bh.bh_row;
      c.owner_rank = rank;
      c.position_stored = bh.position;
      c.velocity_peculiar = bh.velocity;
      c.dynamical_mass_code = state.particles.mass_code[bh.particle_row];
      c.softening_proper = ctx.frame.properLength(P.bh_softening_stored_code);
      cand.push_back(c);
    }
    const BhMergerCriteria crit{M.merger_radius_softening_factor, ctx.frame.newtonGCode()};
    const std::vector<BhMergerPair> pairs = selectDisjointMergerPairs(cand, crit, ctx.frame, P.box);
    if (!pairs.empty()) {
      std::vector<std::uint32_t> outbound;
      for (const BhMergerPair& pair : pairs) {
        const LocalBh& s = bhs[pair.survivor_index];
        const LocalBh& a = bhs[pair.absorbed_index];
        const double ms = state.particles.mass_code[s.particle_row];
        const double ma = state.particles.mass_code[a.particle_row];
        const double mt = ms + ma;
        const Vec3 d = P.box.delta(s.position, a.position);
        Vec3 pos{s.position[0] + d[0] * ma / mt, s.position[1] + d[1] * ma / mt, s.position[2] + d[2] * ma / mt};
        pos = P.box.wrap(pos);
        const Vec3 vel{(ms * s.velocity[0] + ma * a.velocity[0]) / mt, (ms * s.velocity[1] + ma * a.velocity[1]) / mt,
                       (ms * s.velocity[2] + ma * a.velocity[2]) / mt};
        const BlackHoleEngineRecord merged = mergeBlackHoleEngineRecords(store.get(s.id), store.get(a.id), a.id);
        store.put(s.id, merged);
        store.erase(a.id);
        const std::size_t ps = s.particle_row;
        state.particles.position_x_comoving[ps] = pos[0];
        state.particles.position_y_comoving[ps] = pos[1];
        state.particles.position_z_comoving[ps] = pos[2];
        state.particles.velocity_x_peculiar[ps] = vel[0];
        state.particles.velocity_y_peculiar[ps] = vel[1];
        state.particles.velocity_z_peculiar[ps] = vel[2];
        state.particles.mass_code[ps] = merged.dynamicalMassCode();
        state.black_holes.subgrid_mass_code[s.bh_row] = merged.horizon_mass_code;
        state.black_holes.cumulative_accreted_mass_code[s.bh_row] +=
            state.black_holes.cumulative_accreted_mass_code[a.bh_row];
        outbound.push_back(a.particle_row);
        report.removed_particle_ids.push_back(a.id);
        ++report.mergers;
      }
      // Remove absorbed particles through the existing migration commit (a
      // local outbound with no destination), then remap scheduler lanes by ID.
      std::vector<std::uint64_t> ids_before(state.particle_sidecar.particle_id.begin(),
                                            state.particle_sidecar.particle_id.end());
      std::sort(outbound.begin(), outbound.end());
      core::ParticleMigrationCommit commit;
      commit.world_rank = static_cast<int>(rank);
      commit.outbound_local_indices = outbound;
      commit.preserve_gas_cell_state = true;
      state.commitParticleMigration(commit);
      if (ctx.particle_scheduler != nullptr) {
        std::vector<std::uint64_t> ids_after(state.particle_sidecar.particle_id.begin(),
                                             state.particle_sidecar.particle_id.end());
        core::remapSchedulerByParticleId(*ctx.particle_scheduler, ids_before, ids_after);
      }
    }
  }
  store.store(state);
  if (report.particle_mass_mutated || report.mergers > 0U) state.notifyParticleMassMutation();
  std::sort(report.touched_gas_rows.begin(), report.touched_gas_rows.end());
  report.touched_gas_rows.erase(std::unique(report.touched_gas_rows.begin(), report.touched_gas_rows.end()),
                                report.touched_gas_rows.end());
  return report;
}

BlackHoleEngineParameters makeBlackHoleEngineParameters(const core::PhysicsConfig& physics,
                                                        const core::UnitSystem& units, double softening_stored,
                                                        const PeriodicBox& box, double mass_floor_code) {
  BlackHoleEngineParameters p;
  p.model = physics.bh_engine;
  p.efficiencies.radiative = physics.bh_epsilon_r;
  p.efficiencies.thermal = physics.agn.thermal_efficiency;
  p.efficiencies.kick = physics.agn.kick_efficiency;
  p.efficiencies.wind = physics.agn.wind_efficiency;
  p.efficiencies.jet = physics.agn.jet_efficiency;
  p.alpha_bondi = physics.bh_alpha_bondi;
  p.proton_mass_si = physics.bh_proton_mass_si;
  p.thomson_cross_section_si = physics.bh_thomson_cross_section_si;
  const double time_si = units.timeSiPerCode();
  p.viscous_time_code = physics.bh_engine.disk_viscous_time_years * k_year_si / time_si;
  const double kpc_code = core::constants::k_kiloparsec_si / units.length_si_per_code;
  p.kernel_radius_min_proper_code = physics.bh_engine.kernel_radius_min_kpc_phys * kpc_code;
  p.kernel_radius_max_proper_code = physics.bh_engine.kernel_radius_max_kpc_phys * kpc_code;
  p.bh_softening_stored_code = softening_stored;
  p.box = box;
  p.mass_floor_code = mass_floor_code;
  return p;
}

}  // namespace cosmosim::physics
