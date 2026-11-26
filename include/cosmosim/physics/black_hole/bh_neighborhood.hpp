#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH03/AGN02 support infrastructure: owned-leaf
// geometry, a bounded proximity filter and the owner-target support exchange
// shared by the BH supply, withdrawal and AGN injection rounds.
//
// Protocol (memory/performance companion section 5; AGN companion section 3):
// the BH owner issues a compact POD query per BH and round; the exchange
// delivers it ONLY to ranks whose owned-leaf bounding box intersects the
// query sphere (the local rank included); each participant evaluates the
// query over its OWNED LEAF cells (never ghosts or covered parents) and
// returns a fixed-size sufficient-statistic reply to the origin rank. No
// global recipient list, no root-global BH or cell catalog, no all-rank gas
// copy. Rounds are collective: every rank calls every round (empty ranks
// participate with empty sets). The serial implementation is the identity.
//
// Cost: one O(N_owned_leaf) pass per round with O(1) bucket lookups; work and
// scratch scale with the number of (cell, query) matches, not N_BH x N_cell.

#include <array>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <vector>

#include "cosmosim/physics/black_hole/bh_supply.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"
#include "cosmosim/physics/source_transaction.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// Half widths (stored coordinates) of a gas row from its patch geometry;
// falls back to a cube of volume m/rho for non-AMR descriptors. Returns
// false when neither is available.
[[nodiscard]] bool bhCellHalfWidthStored(const core::SimulationState& state, std::uint32_t row,
                                         std::array<double, 3>& half_width_stored) noexcept;

enum class BhSupportRound : std::uint32_t {
  kRadiusBins = 0,
  kKernelStatistics = 1,
  kCounterRotation = 2,
  kWithdrawal = 3,
  kAgnStatistics = 4,
  kAgnCommit = 5,
};

inline constexpr std::size_t k_bh_query_parameter_count = 24U;
inline constexpr std::size_t k_bh_reply_accumulator_count = 48U;

// Trivially copyable wire record (serialized with memcpy by the MPI exchange).
struct BhSupportQuery {
  std::uint64_t bh_id = 0U;
  std::uint64_t query_index = 0U;   // owner-local slot
  std::uint32_t origin_rank = 0U;
  BhSupportRound round = BhSupportRound::kRadiusBins;
  std::array<double, 3> position_stored{};
  std::array<double, 3> velocity_peculiar{};
  double support_radius_stored = 0.0;   // query sphere (stored coordinates)
  double kernel_radius_proper = 0.0;    // kernel h (proper)
  std::array<double, k_bh_kernel_radius_bins> radius_bins_proper{};
  std::array<double, 3> axis{};
  std::array<double, 3> bulk_velocity{};
  double bh_mass_code = 0.0;
  std::array<double, k_bh_query_parameter_count> parameters{};  // round-specific (documented per round)
};

struct BhSupportReply {
  std::uint64_t bh_id = 0U;
  std::uint64_t query_index = 0U;
  BhSupportRound round = BhSupportRound::kRadiusBins;
  std::uint32_t responder_rank = 0U;
  BhRadiusBinSums bins{};
  BhKernelStatistics stats{};
  ConservedIncrement transferred{};  // actual committed gas increment (withdrawal: gas -> BH sign convention below)
  std::array<double, k_bh_reply_accumulator_count> accumulators{};  // round-specific sums
  std::uint64_t rows = 0U;
};

// Collective owner-target exchange.
class BhSupportExchange {
 public:
  virtual ~BhSupportExchange() = default;
  // Deliver every query to each rank whose owned domain intersects its
  // sphere; returns the queries this rank must evaluate (deterministic order:
  // origin rank, then query_index).
  [[nodiscard]] virtual std::vector<BhSupportQuery> routeQueries(std::span<const BhSupportQuery> local_queries) = 0;
  // Return replies to their origin ranks; returns the replies addressed to
  // this rank (deterministic order: query_index, then responder rank).
  [[nodiscard]] virtual std::vector<BhSupportReply> returnReplies(std::span<const BhSupportReply> replies) = 0;
  // Unanimous vote (logical AND over all ranks).
  [[nodiscard]] virtual bool allAgree(bool local_ok) = 0;
  [[nodiscard]] virtual std::uint32_t worldRank() const noexcept = 0;
};

// Serial identity exchange (single rank).
class LocalBhSupportExchange final : public BhSupportExchange {
 public:
  [[nodiscard]] std::vector<BhSupportQuery> routeQueries(std::span<const BhSupportQuery> local_queries) override;
  [[nodiscard]] std::vector<BhSupportReply> returnReplies(std::span<const BhSupportReply> replies) override;
  [[nodiscard]] bool allAgree(bool local_ok) override { return local_ok; }
  [[nodiscard]] std::uint32_t worldRank() const noexcept override { return 0U; }
};

// Bounded spatial bucket lookup of query spheres (stored coordinates).
class BhProximityFilter {
 public:
  BhProximityFilter(std::span<const BhSupportQuery> queries, const PeriodicBox& box);
  // Indices of queries whose sphere contains `position`; `out` is cleared.
  void matches(const std::array<double, 3>& position, std::vector<std::uint32_t>& out) const;
  [[nodiscard]] bool empty() const noexcept { return m_queries.empty(); }

 private:
  [[nodiscard]] std::int64_t bucketKey(std::int64_t ix, std::int64_t iy, std::int64_t iz) const noexcept;
  [[nodiscard]] std::array<std::int64_t, 3> bucketOf(const std::array<double, 3>& p) const noexcept;

  std::span<const BhSupportQuery> m_queries;
  PeriodicBox m_box;
  double m_bucket = 1.0;
  std::array<std::int64_t, 3> m_bucket_count{1, 1, 1};
  std::unordered_map<std::int64_t, std::vector<std::uint32_t>> m_buckets;
};

// Owned-leaf row predicate: empty mask means every row is an owned leaf.
struct BhOwnedLeafView {
  std::span<const std::uint8_t> owned_leaf_mask;
  [[nodiscard]] bool owned(std::size_t row) const noexcept {
    return owned_leaf_mask.empty() || (row < owned_leaf_mask.size() && owned_leaf_mask[row] != 0U);
  }
};

}  // namespace cosmosim::physics
