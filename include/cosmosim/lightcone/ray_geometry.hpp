#pragma once

// CHUI-EXT-PHYS-REF-v1 S03: read-only leaf ray service shared by every
// observable (Ly-a forest, HI 21 cm, continuum/Stokes, SZ, DM/RM). It never
// mutates simulation state; it maps rays onto the AMR LEAF tiling.
//
// Geometry contract
//   * Leaves are the populated patches (cell_count > 0) of core::PatchSoa;
//     metadata-only parents own no cells, so covered parents can never be
//     traversed (no double counting). Cell rows follow the canonical
//     amr::linearPatchCellIndex order i + nx (j + ny k) from first_cell.
//   * Point location uses a uniform bucket index over the domain (each leaf
//     listed in every bucket it overlaps): O(1) expected per lookup instead
//     of a scan over all patches. Leaves are half-open boxes [lo, hi).
//   * A ray is traced in COMOVING length s from its start point along a unit
//     direction. At every face crossing the next cell is located at the probe
//     point p + nudge d (nudge = 1e-6 of the finest cell width), so ties on
//     faces/edges/corners resolve deterministically toward the direction of
//     travel. Segments are emitted in strictly increasing s.
//   * Boundaries: kPeriodic wraps the position into the domain (replicated
//     boxes along long rays); kOpen terminates the ray at the domain face.
//   * A hole in the tiling (no leaf contains a probe point inside the domain)
//     throws: the geometry contract (exact tiling, amr_geometry_coverage) is
//     violated and no observable may silently skip volume.
//   * Capacity: at most `max_segments` segments per ray; reaching the limit
//     throws (bounded memory; never truncates silently).
// MPI: each rank builds the geometry of its locally owned leaves; additive
// line-of-sight integrals (optical depth, column density, DM, y) are reduced
// across ranks per ray; ORDERED radiative transfer needs the merged segment
// list (mergeRankSegments) on the ray owner (exchange edge documented in the
// ledger, single-rank path implemented).

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::lightcone {

enum class RayBoundary : std::uint8_t { kPeriodic = 0, kOpen = 1 };

struct LeafPatch {
  std::uint64_t patch_id = 0;
  std::int32_t level = 0;
  std::array<double, 3> origin{};
  std::array<double, 3> extent{};
  std::array<std::uint32_t, 3> dims{};
  std::uint32_t first_cell = 0;
};

struct LeafCellRef {
  std::uint32_t leaf = 0;
  std::uint32_t cell_row = 0;
  std::array<std::uint32_t, 3> ijk{};
};

struct RaySegment {
  std::uint64_t patch_id = 0;  // global leaf identity (rank-independent)
  std::uint32_t cell_row = 0;  // rank-local row of the cell in CellSoa
  std::uint32_t leaf = 0;      // rank-local leaf index
  double s_begin = 0.0;  // comoving length along the ray from its start
  double s_end = 0.0;
  std::array<double, 3> midpoint{};  // wrapped comoving position
  [[nodiscard]] double length() const noexcept { return s_end - s_begin; }
};

class LeafGeometry {
 public:
  LeafGeometry(std::vector<LeafPatch> leaves, std::array<double, 3> domain_origin, std::array<double, 3> domain_extent,
               RayBoundary boundary, std::size_t buckets_per_axis = 0);
  // Populated patches of a SimulationState (legacy undeclared geometry is refused).
  [[nodiscard]] static LeafGeometry fromPatches(const core::PatchSoa& patches, std::array<double, 3> domain_origin,
                                                std::array<double, 3> domain_extent, RayBoundary boundary);
  // A single uniform n_x x n_y x n_z grid filling the domain (rows 0..N-1).
  [[nodiscard]] static LeafGeometry uniformGrid(std::array<std::uint32_t, 3> dims, std::array<double, 3> domain_origin,
                                                std::array<double, 3> domain_extent, RayBoundary boundary);

  [[nodiscard]] bool locate(const std::array<double, 3>& point, LeafCellRef& out) const;
  [[nodiscard]] std::vector<RaySegment> traceStraightRay(const std::array<double, 3>& start,
                                                         const std::array<double, 3>& unit_direction,
                                                         double max_length, std::size_t max_segments) const;
  [[nodiscard]] std::size_t leafCount() const noexcept { return m_leaves.size(); }
  [[nodiscard]] std::size_t cellCount() const noexcept { return m_cell_count; }
  [[nodiscard]] const LeafPatch& leaf(std::size_t i) const { return m_leaves.at(i); }
  [[nodiscard]] double finestCellWidth() const noexcept { return m_finest_width; }
  [[nodiscard]] std::array<double, 3> domainOrigin() const noexcept { return m_origin; }
  [[nodiscard]] std::array<double, 3> domainExtent() const noexcept { return m_extent; }
  [[nodiscard]] RayBoundary boundary() const noexcept { return m_boundary; }
  [[nodiscard]] std::uint64_t indexBytes() const noexcept;
  [[nodiscard]] std::array<double, 3> wrap(const std::array<double, 3>& p) const noexcept;
  // Cell center of a located cell (comoving).
  [[nodiscard]] std::array<double, 3> cellCenter(const LeafCellRef& ref) const;
  [[nodiscard]] double cellWidth(std::uint32_t leaf, std::size_t axis) const;

 private:
  [[nodiscard]] std::size_t bucketOf(const std::array<double, 3>& p) const noexcept;
  std::vector<LeafPatch> m_leaves;
  std::array<double, 3> m_origin{};
  std::array<double, 3> m_extent{};
  RayBoundary m_boundary = RayBoundary::kPeriodic;
  std::size_t m_buckets = 1;
  std::vector<std::uint32_t> m_bucket_begin;  // CSR offsets, size buckets^3 + 1
  std::vector<std::uint32_t> m_bucket_leaves;
  std::size_t m_cell_count = 0;
  double m_finest_width = 0.0;
};

// Merges per-rank segment lists of one ray into increasing s (ties broken by
// patch id so the order is rank-count independent); throws on overlap (two
// owners of one path interval, relative tolerance 1e-12 of the ray length).
[[nodiscard]] std::vector<RaySegment> mergeRankSegments(std::vector<std::vector<RaySegment>> per_rank);

// Unit vector from spherical angles (theta from +z, phi from +x toward +y).
[[nodiscard]] std::array<double, 3> unitVectorFromAngles(double theta, double phi);

}  // namespace cosmosim::lightcone
