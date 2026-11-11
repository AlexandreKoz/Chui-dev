#pragma once

// CHUI patch 58 (P1-GEOMETRY): exact AMR leaf-tiling proof.
//
// Problem. The production AMR hydro domain is the union of the populated
// (leaf) patches of every rank. Ghost filling treats the faces of that union's
// bounding box as the only legal non-periodic walls, so the populated leaves
// must tile the box exactly: no hole (would be a fictitious internal wall) and
// no overlap (two owners of one volume). A summed-volume comparison does not
// prove this — an overlap can cancel a hole of equal volume.
//
// Method (exact, any number of levels/ranks/patch shapes). Every populated
// patch is mapped onto one integer lattice whose unit is the finest populated
// cell width per axis; the mapping is accepted only when every patch face lies
// on a lattice plane (within a documented roundoff tolerance) and its cell
// width is an integer number of lattice units. For a half-open lattice box
// P = [a, b) the indicator is, by inclusion-exclusion,
//
//     1_P(x) = sum_{c in corners(P)} s_P(c) H(x - c),
//     H(y)   = prod_axis [y_axis >= 0],
//     s_P(c) = (-1)^(number of axes on which c takes the upper bound b).
//
// Summing over all patches, the coverage multiplicity m(x) = sum_P 1_P(x) is
// the signed point measure W = sum_P sum_c s_P(c) delta_c pushed through H.
// The orthant functions H(x - c) of distinct lattice points are linearly
// independent (evaluate at a support point that is minimal in the product
// order: only its own term is non-zero), so
//
//     m == 1_B  <=>  W == W_B  (the 8 signed corners of the bounding box B).
//
// m == 1_B means every lattice cell of B is covered exactly once and nothing
// outside B is covered: non-overlap and complete coverage, simultaneously and
// exactly, in integer arithmetic. Corners of neighbouring patches cancel, so
// after aggregation a valid tiling leaves exactly 8 points, and a rank's
// partial aggregate only retains the corners of its own owned region
// (surface-proportional), which bounds the distributed reduction.
//
// This header is MPI-free. The collective orchestration (vote, reductions,
// bounded gather, memory admission, cache) lives in
// src/workflows/internal/amr_geometry_assurance.{hpp,cpp}.

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

namespace cosmosim::amr {

// One populated patch as seen by the coverage proof. `cell_count` is the
// number of real cells the patch owns (0 = metadata-only/inactive parent,
// which is never part of the tiling).
struct AmrPatchGeometryRecord {
  std::uint64_t patch_id = 0;
  std::int32_t level = 0;
  std::uint32_t owning_rank = 0;
  std::uint32_t cell_count = 0;
  std::array<double, 3> origin_comoving{0.0, 0.0, 0.0};
  std::array<double, 3> extent_comoving{0.0, 0.0, 0.0};
  std::array<std::uint32_t, 3> cell_dims{0U, 0U, 0U};
};

// All six geometry lanes are zero: the documented legacy non-AMR descriptor
// (PatchSoa comment). It declares no geometry and can never be a leaf of a
// production AMR tiling.
[[nodiscard]] bool isLegacyUndeclaredAmrPatchGeometry(const AmrPatchGeometryRecord& record) noexcept;

// Empty string when the record is a well-formed geometric leaf; otherwise a
// diagnostic naming patch ID, rank, level, axis and offending value. Checks:
// level >= 0; cell_count > 0 and equal to the product of the cell dimensions
// (overflow-checked); per axis finite origin, finite positive extent, finite
// upper bound strictly above the origin (catches absorption overflow), finite
// positive cell width.
[[nodiscard]] std::string diagnoseAmrPatchGeometry(const AmrPatchGeometryRecord& record);

// Largest lattice span per axis admitted by the proof (2^40 finest cells). It
// keeps every lattice coordinate exactly representable in double and int64
// and keeps the snap tolerance below 1e-3 lattice units.
inline constexpr std::int64_t k_amr_coverage_max_lattice_cells_per_axis = std::int64_t{1} << 40;

struct AmrCoverageLattice {
  std::array<double, 3> lower_comoving{0.0, 0.0, 0.0};
  std::array<double, 3> unit_comoving{0.0, 0.0, 0.0};  // finest populated cell width
  std::array<std::int64_t, 3> cells{0, 0, 0};          // bounding box in lattice units
};

// Builds the lattice from the (globally reduced) bounding box and finest cell
// width. Throws std::runtime_error when the bounding box is not an integer
// number of lattice units, is degenerate, or exceeds the admitted span.
[[nodiscard]] AmrCoverageLattice makeAmrCoverageLattice(
    const std::array<double, 3>& lower_comoving,
    const std::array<double, 3>& upper_comoving,
    const std::array<double, 3>& finest_cell_width_comoving);

// Snap tolerance in lattice units for a coordinate of magnitude
// `span_lattice_units`: 1e-6 + 64 * eps * span (absorbs the roundoff of
// origin + k * width arithmetic; always << 0.5).
[[nodiscard]] double amrCoverageSnapTolerance(double span_lattice_units) noexcept;

struct AmrSignedLatticeCorner {
  std::array<std::int64_t, 3> index{0, 0, 0};
  std::int64_t weight = 0;
};
static_assert(std::is_trivially_copyable_v<AmrSignedLatticeCorner>);
static_assert(sizeof(AmrSignedLatticeCorner) == 32U);

// Appends the 8 signed corners of `record` (already accepted by
// diagnoseAmrPatchGeometry). Throws std::runtime_error naming the patch when a
// face is not on a lattice plane, a cell width is not an integer number of
// lattice units, or the patch leaves the lattice box.
void appendAmrPatchSignedCorners(
    const AmrCoverageLattice& lattice,
    const AmrPatchGeometryRecord& record,
    std::vector<AmrSignedLatticeCorner>& corners);

// Sorts by lattice index, merges equal points (checked int64 addition) and
// removes zero-weight points. O(n log n), in place.
void aggregateAmrSignedCorners(std::vector<AmrSignedLatticeCorner>& corners);

struct AmrTilingVerdict {
  bool exact_tiling = false;
  // Aggregated points that differ from the bounding-box measure (0 when exact).
  std::uint64_t residual_points = 0;
  std::string diagnostic;
};

// `aggregated` must be the aggregate of ALL populated patches (all ranks).
// Exact verdict; the diagnostic lists up to four residual points in comoving
// coordinates and classifies the defect from the summed lattice volume
// (diagnostic only: overlap, hole, or overlap compensating a hole).
[[nodiscard]] AmrTilingVerdict verifyAmrExactTiling(
    const AmrCoverageLattice& lattice,
    std::span<const AmrSignedLatticeCorner> aggregated,
    long double summed_lattice_volume);

// Lattice volume of one accepted record (diagnostic classification only).
[[nodiscard]] long double amrPatchLatticeVolume(
    const AmrCoverageLattice& lattice, const AmrPatchGeometryRecord& record);

// Serial convenience (one address space; tests and single-rank callers):
// diagnoses every record, builds the lattice from their bounding box and
// finest width, and returns the exact verdict. Never throws for malformed
// geometry: the verdict carries the diagnostic instead.
[[nodiscard]] AmrTilingVerdict verifyAmrPatchSetTiling(
    std::span<const AmrPatchGeometryRecord> records);

}  // namespace cosmosim::amr
