// C04-04 sampling identity, global K and deterministic aggregation fixtures.
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// The same physical regions must be selected for 1, 2 and 4 rank
// decompositions (including uneven and empty ranks) with one global K, no
// rank-local ID-zero duplicates, and published aggregates that do not depend
// on rank order.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <tuple>
#include <vector>

#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_shadow_reference.hpp"

namespace {
namespace h = cosmosim::hydro;

struct Region {
  int level = 0;
  std::array<std::int64_t, 3> origin{};
  std::array<std::uint32_t, 3> extent{8U, 8U, 8U};
};

[[nodiscard]] std::vector<Region> physicalRegions() {
  std::vector<Region> regions;
  for (std::int64_t k = 0; k < 4; ++k) {
    for (std::int64_t j = 0; j < 4; ++j) {
      for (std::int64_t i = 0; i < 4; ++i) {
        regions.push_back(Region{.level = 0, .origin = {8 * i, 8 * j, 8 * k}});
      }
    }
  }
  for (std::int64_t i = 0; i < 6; ++i) {
    regions.push_back(Region{.level = 1, .origin = {16 + 8 * i, 16, 0}});
  }
  return regions;
}

[[nodiscard]] h::HydroShadowCandidate candidate(const Region& region, std::uint64_t rank_local_patch_id) {
  return h::HydroShadowCandidate{
      .patch_id = rank_local_patch_id,
      .level = region.level,
      .region_key = h::hydroShadowPhysicalRegionKey(region.level, region.origin, region.extent)};
}

// Simulates the runtime protocol: every rank submits its local top K, the
// bounded allgather concatenates them in rank order, every rank ranks again.
[[nodiscard]] std::vector<h::HydroShadowCandidate> globalSelection(
    const std::vector<std::vector<h::HydroShadowCandidate>>& per_rank, std::size_t k,
    std::uint64_t step, std::uint64_t seed) {
  std::vector<h::HydroShadowCandidate> wire;
  for (const auto& local : per_rank) {
    const auto top = h::rankHydroShadowCandidates(local, k, step, seed);
    wire.insert(wire.end(), top.begin(), top.end());
  }
  // The wire carries only (region_key, level): rank-local patch IDs are lost.
  for (auto& entry : wire) entry.patch_id = 0U;
  return h::rankHydroShadowCandidates(wire, k, step, seed);
}

void checkSelectionIndependentOfDecomposition() {
  const auto regions = physicalRegions();
  constexpr std::size_t k = 5U;
  constexpr std::uint64_t step = 48U;
  constexpr std::uint64_t seed = 20261117U;
  // Reference: one rank, patch IDs 1..N.
  std::vector<std::vector<h::HydroShadowCandidate>> single(1);
  for (std::size_t r = 0; r < regions.size(); ++r) single[0].push_back(candidate(regions[r], r + 1U));
  const auto expected = globalSelection(single, k, step, seed);
  assert(expected.size() == k);
  for (const unsigned ranks : {2U, 4U}) {
    for (const bool uneven : {false, true}) {
      std::vector<std::vector<h::HydroShadowCandidate>> per_rank(ranks);
      for (std::size_t r = 0; r < regions.size(); ++r) {
        // Different owner mapping, different (rank-local, even colliding)
        // patch IDs, reversed submission order. Uneven: rank 0 owns 3/4, the
        // last rank owns nothing (zero-cell rank still takes part).
        unsigned owner = static_cast<unsigned>((r * 7U + ranks) % ranks);
        if (uneven) owner = (r % 4U == 0U && ranks > 2U) ? 1U : 0U;
        per_rank[owner].insert(per_rank[owner].begin(), candidate(regions[r], 0U));
      }
      if (uneven) assert(per_rank.back().empty());
      const auto selected = globalSelection(per_rank, k, step, seed);
      assert(selected.size() == expected.size());
      for (std::size_t i = 0; i < expected.size(); ++i) {
        assert(selected[i].region_key == expected[i].region_key);
        assert(selected[i].level == expected[i].level);
      }
    }
  }
}

// The former fixed-grid adapter submitted {patch_id 0, level 0} from EVERY
// rank. Physical keys never collide that way: K=1 selects exactly one region.
void checkNoRankLocalIdZeroDuplicates() {
  std::vector<std::vector<h::HydroShadowCandidate>> per_rank(4);
  for (std::size_t rank = 0; rank < per_rank.size(); ++rank) {
    per_rank[rank].push_back(candidate(Region{.level = 0, .origin = {8 * static_cast<std::int64_t>(rank), 0, 0}}, 0U));
  }
  const auto selected = globalSelection(per_rank, 1U, 0U, 0U);
  assert(selected.size() == 1U);
  // A duplicated physical region (two ranks claiming the same key) collapses.
  std::vector<std::vector<h::HydroShadowCandidate>> duplicated(2);
  duplicated[0].push_back(candidate(Region{.level = 0, .origin = {0, 0, 0}}, 5U));
  duplicated[1].push_back(candidate(Region{.level = 0, .origin = {0, 0, 0}}, 9U));
  assert(globalSelection(duplicated, 4U, 0U, 0U).size() == 1U);
  // Global K is a single global bound (not per rank).
  std::vector<std::vector<h::HydroShadowCandidate>> many(4);
  const auto regions = physicalRegions();
  for (std::size_t r = 0; r < regions.size(); ++r) many[r % 4U].push_back(candidate(regions[r], 0U));
  assert(globalSelection(many, 3U, 7U, 1U).size() == 3U);
}

void checkKeysAndLattice() {
  const auto a = h::hydroShadowPhysicalRegionKey(0, {0, 0, 0}, {8U, 8U, 8U});
  assert(a != 0U);
  assert(a == h::hydroShadowPhysicalRegionKey(0, {0, 0, 0}, {8U, 8U, 8U}));
  assert(a != h::hydroShadowPhysicalRegionKey(1, {0, 0, 0}, {8U, 8U, 8U}));
  assert(a != h::hydroShadowPhysicalRegionKey(0, {8, 0, 0}, {8U, 8U, 8U}));
  assert(a != h::hydroShadowPhysicalRegionKey(0, {0, 0, 0}, {8U, 8U, 1U}));
  const auto origin = h::hydroShadowLatticeOrigin({0.5, 0.25, 0.0}, {0.0, 0.0, 0.0}, {0.0625, 0.0625, 0.0625});
  assert(origin.has_value() && (*origin)[0] == 8 && (*origin)[1] == 4 && (*origin)[2] == 0);
  assert(!h::hydroShadowLatticeOrigin({0.51, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0625, 0.0625, 0.0625}).has_value());
  assert(!h::hydroShadowLatticeOrigin({0.5, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0625, 0.0625}).has_value());
}

// Region results are reduced after sorting by physical key: the published
// floating sums are independent of the rank order in which they arrive.
void checkDeterministicRegionAggregation() {
  std::vector<h::HydroShadowRegionResult> results;
  for (std::uint64_t r = 0; r < 9U; ++r) {
    h::HydroShadowRegionResult result;
    result.region_key = 1000U + 37U * ((r * 5U) % 9U);
    result.level = static_cast<int>(r % 2U);
    result.comparison_class = static_cast<h::HydroShadowComparisonClass>(r % 4U);
    result.compared_cells = 100U + r;
    result.compared_volume = 0.1 * static_cast<double>(r + 1U);
    for (std::size_t lane = 0; lane < 7U; ++lane) {
      result.components[lane].weighted_abs_difference = 1.0e-16 * static_cast<double>((r + 3U) * (lane + 1U)) + 1.0;
      result.components[lane].linf_difference = static_cast<double>(r) * 1.0e-3;
    }
    results.push_back(result);
  }
  const auto reduce = [](std::vector<h::HydroShadowRegionResult> in) {
    std::vector<std::uint8_t> wire;
    for (const auto& result : in) h::appendHydroShadowRegionResultWire(result, wire);
    assert(wire.size() == in.size() * h::k_hydro_shadow_region_result_wire_bytes);
    std::vector<h::HydroShadowRegionResult> decoded;
    for (std::size_t offset = 0; offset < wire.size(); offset += h::k_hydro_shadow_region_result_wire_bytes) {
      decoded.push_back(h::readHydroShadowRegionResultWire(
          std::span<const std::uint8_t>(wire.data() + offset, h::k_hydro_shadow_region_result_wire_bytes)));
    }
    std::stable_sort(decoded.begin(), decoded.end(), [](const auto& a, const auto& b) {
      return std::tuple{a.region_key, a.level, static_cast<int>(a.comparison_class), a.compared_cells,
                        a.compared_composition_cells} <
             std::tuple{b.region_key, b.level, static_cast<int>(b.comparison_class), b.compared_cells,
                        b.compared_composition_cells};
    });
    double sum = 0.0;
    for (const auto& result : decoded) sum += result.components[4].weighted_abs_difference;
    return sum;
  };
  const double forward = reduce(results);
  std::reverse(results.begin(), results.end());
  const double backward = reduce(results);
  std::rotate(results.begin(), results.begin() + 4, results.end());
  const double rotated = reduce(results);
  assert(std::memcmp(&forward, &backward, sizeof(double)) == 0);
  assert(std::memcmp(&forward, &rotated, sizeof(double)) == 0);
  // Wire round trip preserves every field exactly.
  std::vector<std::uint8_t> one;
  h::appendHydroShadowRegionResultWire(results.front(), one);
  const auto decoded = h::readHydroShadowRegionResultWire(one);
  assert(decoded.region_key == results.front().region_key);
  assert(decoded.comparison_class == results.front().comparison_class);
  assert(decoded.compared_volume == results.front().compared_volume);
}

// Dependency cones on an 8^3 patch: one stale face slot excludes only the
// cells within the halfstep stencil radius of its owner; suppressed axes are
// ignored (2D patch).
void checkDependencyConeOn8Cubed() {
  auto geometry = h::makeCartesianPatchGeometry({.nx = 8U, .ny = 8U, .nz = 8U,
                                                 .cell_width_x_comoving = 1.0,
                                                 .cell_width_y_comoving = 1.0,
                                                 .cell_width_z_comoving = 1.0});
  h::appendCartesianBoundaryGhostFaces(geometry, h::HydroBoundaryKind::kReflective);
  std::vector<std::uint8_t> stale(geometry.ghost_cells.size(), 0U);
  // Mark every x-lower face slot stale.
  std::size_t stale_slots = 0U;
  for (std::size_t slot = 0; slot < geometry.ghost_cells.size(); ++slot) {
    const auto& ghost = geometry.ghost_cells[slot];
    if (ghost.axis == h::HydroFaceAxis::kX && ghost.side == h::HydroFaceSide::kLower) {
      stale[slot] = 1U;
      ++stale_slots;
    }
  }
  assert(stale_slots == 64U);
  const auto cone = h::hydroShadowDependencyConeMask(geometry, stale, h::k_hydro_shadow_halfstep_dependency_radius);
  std::size_t excluded = 0U;
  for (std::size_t row = 0; row < geometry.cellCount(); ++row) {
    const auto ijk = geometry.cellIjk(row);
    const bool expect = ijk[0] <= static_cast<std::size_t>(h::k_hydro_shadow_halfstep_dependency_radius);
    assert((cone[row] != 0U) == expect);
    excluded += cone[row];
  }
  assert(excluded == 5U * 64U);
  // 2D (nz = 1): the suppressed axis never widens the cone.
  auto flat = h::makeCartesianPatchGeometry({.nx = 8U, .ny = 8U, .nz = 1U, .cell_width_x_comoving = 1.0,
                                             .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  h::appendCartesianBoundaryGhostFaces(flat, h::HydroBoundaryKind::kOpen);
  std::vector<std::uint8_t> one(flat.ghost_cells.size(), 0U);
  one.front() = 1U;
  const auto flat_cone = h::hydroShadowDependencyConeMask(flat, one, 1);
  std::size_t flat_excluded = 0U;
  for (const auto flag : flat_cone) flat_excluded += flag;
  assert(flat_excluded <= 9U && flat_excluded >= 4U);
}

}  // namespace

int main() {
  checkKeysAndLattice();
  checkSelectionIndependentOfDecomposition();
  checkNoRankLocalIdZeroDuplicates();
  checkDeterministicRegionAggregation();
  checkDependencyConeOn8Cubed();
  return 0;
}
