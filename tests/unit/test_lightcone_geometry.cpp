// CHUI-EXT-PHYS-REF-v1 gates VI/VC/VL (S03): leaf ray traversal (ordering,
// path-length closure, AMR leaves, covered-parent exclusion, ties, periodic
// and open boundaries, hole and capacity refusal), rank segment merging,
// comoving-distance tables, lightcone shells, replication and temporal
// interpolation. Expected values are analytic. NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/lightcone/lightcone_plan.hpp"
#include "cosmosim/lightcone/ray_geometry.hpp"

namespace lc = cosmosim::lightcone;
namespace c = cosmosim::core;

namespace {

bool closeAbs(double a, double b, double tol) { return std::abs(a - b) <= tol; }

// EXT-S03-T01: axis-aligned periodic ray through a uniform 4^3 grid.
void testUniformAxisRay() {
  const auto geometry = lc::LeafGeometry::uniformGrid({4, 4, 4}, {0, 0, 0}, {1, 1, 1}, lc::RayBoundary::kPeriodic);
  const auto segments = geometry.traceStraightRay({0.05, 0.3, 0.6}, {1, 0, 0}, 1.0, 64);
  CHUI_REQUIRE(segments.size() == 5U);
  const double expected_length[5] = {0.2, 0.25, 0.25, 0.25, 0.05};
  const std::uint32_t expected_ix[5] = {0, 1, 2, 3, 0};
  double total = 0.0;
  for (std::size_t i = 0; i < 5U; ++i) {
    CHUI_REQUIRE(closeAbs(segments[i].length(), expected_length[i], 1.0e-12));
    CHUI_REQUIRE(segments[i].cell_row == expected_ix[i] + 4U * (1U + 4U * 2U));
    if (i > 0) CHUI_REQUIRE(closeAbs(segments[i].s_begin, segments[i - 1].s_end, 0.0));
    total += segments[i].length();
  }
  CHUI_REQUIRE(closeAbs(total, 1.0, 1.0e-14));
}

// EXT-S03-T02: oblique ray, midpoints lie in their cells, lengths close.
void testObliqueRay() {
  const auto geometry = lc::LeafGeometry::uniformGrid({8, 8, 8}, {0, 0, 0}, {2, 2, 2}, lc::RayBoundary::kPeriodic);
  const double n = std::sqrt(1.0 + 4.0 + 9.0);
  const std::array<double, 3> d{1.0 / n, 2.0 / n, 3.0 / n};
  const auto segments = geometry.traceStraightRay({0.1, 0.2, 0.3}, d, 7.5, 4096);
  double total = 0.0;
  for (const auto& seg : segments) {
    lc::LeafCellRef ref;
    CHUI_REQUIRE(geometry.locate(seg.midpoint, ref));
    CHUI_REQUIRE(ref.cell_row == seg.cell_row);
    CHUI_REQUIRE(seg.length() > 0.0);
    total += seg.length();
  }
  CHUI_REQUIRE(closeAbs(total, 7.5, 1.0e-12));
  // A ray exactly along a grid edge resolves ties deterministically (two runs identical).
  const auto a = geometry.traceStraightRay({0.25, 0.25, 0.1}, {0, 0, 1}, 1.0, 64);
  const auto b = geometry.traceStraightRay({0.25, 0.25, 0.1}, {0, 0, 1}, 1.0, 64);
  CHUI_REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) CHUI_REQUIRE(a[i].cell_row == b[i].cell_row);
}

// EXT-S03-T03: two-level leaf tiling through PatchSoa with a metadata-only
// parent that must never be traversed; open boundary; hole and capacity.
void testAmrLeavesAndRefusals() {
  c::PatchSoa patches;
  patches.resize(3);
  // Parent covering the whole domain (cell_count 0: metadata only).
  const auto set = [&](std::size_t p, std::uint64_t id, int level, std::array<double, 3> o, std::array<double, 3> e,
                       std::array<std::uint16_t, 3> dims, std::uint32_t first, std::uint32_t count) {
    patches.patch_id[p] = id;
    patches.level[p] = level;
    patches.origin_x_comoving[p] = o[0];
    patches.origin_y_comoving[p] = o[1];
    patches.origin_z_comoving[p] = o[2];
    patches.extent_x_comoving[p] = e[0];
    patches.extent_y_comoving[p] = e[1];
    patches.extent_z_comoving[p] = e[2];
    patches.cell_dim_x[p] = dims[0];
    patches.cell_dim_y[p] = dims[1];
    patches.cell_dim_z[p] = dims[2];
    patches.first_cell[p] = first;
    patches.cell_count[p] = count;
  };
  set(0, 100, 0, {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, 0, 0);          // parent (not a leaf)
  set(1, 101, 1, {0, 0, 0}, {0.5, 1, 1}, {2, 4, 4}, 0, 32);       // fine left half
  set(2, 102, 0, {0.5, 0, 0}, {0.5, 1, 1}, {1, 2, 2}, 32, 4);     // coarse right half
  const auto geometry = lc::LeafGeometry::fromPatches(patches, {0, 0, 0}, {1, 1, 1}, lc::RayBoundary::kOpen);
  CHUI_REQUIRE(geometry.leafCount() == 2U && geometry.cellCount() == 36U);
  const auto segments = geometry.traceStraightRay({0.01, 0.1, 0.1}, {1, 0, 0}, 10.0, 64);
  CHUI_REQUIRE(segments.size() == 3U);  // two fine cells, one coarse cell, then the open face
  CHUI_REQUIRE(segments[0].patch_id == 101U && segments[2].patch_id == 102U);
  CHUI_REQUIRE(closeAbs(segments[2].s_end, 0.99, 1.0e-12));
  CHUI_REQUIRE(segments[2].cell_row == 32U);
  // Hole: remove the right half -> probe in an uncovered region throws.
  const auto holey = lc::LeafGeometry({geometry.leaf(0)}, {0, 0, 0}, {1, 1, 1}, lc::RayBoundary::kOpen);
  bool threw = false;
  try {
    (void)holey.traceStraightRay({0.01, 0.1, 0.1}, {1, 0, 0}, 10.0, 64);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  threw = false;
  try {
    (void)geometry.traceStraightRay({0.01, 0.1, 0.1}, {1, 0, 0}, 10.0, 2);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

// EXT-S03-T04: merging per-rank segment lists and overlap refusal.
void testMerge() {
  const auto geometry = lc::LeafGeometry::uniformGrid({4, 4, 4}, {0, 0, 0}, {1, 1, 1}, lc::RayBoundary::kPeriodic);
  const auto all = geometry.traceStraightRay({0.05, 0.3, 0.6}, {1, 0, 0}, 1.0, 64);
  std::vector<lc::RaySegment> rank0{all[1], all[3]};
  std::vector<lc::RaySegment> rank1{all[0], all[2], all[4]};
  const auto merged = lc::mergeRankSegments({rank1, rank0});
  CHUI_REQUIRE(merged.size() == all.size());
  for (std::size_t i = 0; i < all.size(); ++i) CHUI_REQUIRE(merged[i].s_begin == all[i].s_begin);
  bool threw = false;
  try {
    (void)lc::mergeRankSegments({std::vector<lc::RaySegment>{all[1]}, std::vector<lc::RaySegment>{all[1]}});
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

// EXT-S03-T05: comoving distances (EdS analytic), shells, replication, time weights.
void testLightconePlan() {
  c::CosmologyBackgroundConfig eds;
  eds.hubble_param = 0.7;
  eds.omega_matter = 1.0;
  eds.omega_lambda = 0.0;
  const c::LambdaCdmBackground bg(eds);
  const lc::ComovingDistanceTable table(bg, 0.05, 4096);
  const double dh = cosmosim::core::constants::k_speed_of_light_si / bg.hubble0Si();
  for (double a : {0.06, 0.1, 0.33, 0.77, 1.0}) {
    const double expected = 2.0 * dh * (1.0 - std::sqrt(a));
    CHUI_REQUIRE(std::abs(table.comovingDistanceM(a) - expected) <= 1.0e-8 * dh);
    if (a < 1.0) CHUI_REQUIRE(std::abs(table.scaleFactorAt(expected) - a) <= 1.0e-8 * a);
  }
  const std::vector<double> snaps{0.2, 0.25, 0.5, 1.0};
  const auto shells = lc::planLightconeShells(table, snaps);
  CHUI_REQUIRE(shells.size() == 4U);
  for (std::size_t i = 1; i < shells.size(); ++i) CHUI_REQUIRE(shells[i].chi_outer_m == shells[i - 1].chi_inner_m);
  CHUI_REQUIRE(shells.back().chi_inner_m == 0.0 && shells.back().a_inner == 1.0);
  const auto meta = lc::replicationMetadata(3.0e24, 1.0e26);
  CHUI_REQUIRE(meta.replicas_along_ray == 34U && meta.repeated_structures);
  const auto w = lc::temporalWeights(snaps, 0.375);
  CHUI_REQUIRE(w.lower_snapshot == 1U && w.upper_snapshot == 2U && std::abs(w.upper_weight - 0.5) < 1.0e-14);
  bool threw = false;
  try {
    (void)lc::temporalWeights(snaps, 0.1);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  CHUI_REQUIRE(std::abs(lc::observedRedshift(1.0, 0.0) - 1.0) < 1.0e-15);
}

}  // namespace

int main() {
  testUniformAxisRay();
  testObliqueRay();
  testAmrLeavesAndRefusals();
  testMerge();
  testLightconePlan();
  return 0;
}
