#include "cosmosim/amr/amr_ghost_fill.hpp"
#include "amr/internal/amr_retained_state_admission.hpp"
#include "amr/internal/amr_periodic_topology.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory_resource>

#include "cosmosim/core/governed_scratch_arena.hpp"
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"

namespace cosmosim::amr {
namespace {

constexpr double k_geometry_tol = 1.0e-10;
constexpr double k_time_tol = 1.0e-12;

struct SelectedSourceCell {
  const AmrHydroGhostFillPatch* patch = nullptr;
  const AmrHydroSparseGhostSource* sparse_remote = nullptr;
  std::size_t cell = 0U;
  std::size_t fine_stencil_ordinal = 0U;
};

struct CellSelection {
  std::pmr::vector<SelectedSourceCell> cells;
  explicit CellSelection(std::pmr::memory_resource* resource) : cells(resource) {}
};

[[nodiscard]] const PatchDescriptor& sourcePatch(const SelectedSourceCell& source) {
  if (source.patch != nullptr) {
    return source.patch->geometry->patch;
  }
  if (source.sparse_remote != nullptr && source.sparse_remote->patch != nullptr) {
    return *source.sparse_remote->patch;
  }
  throw std::logic_error("fillAmrHydroGhostCells: empty selected source cell");
}

[[nodiscard]] const PatchDescriptor& selectionPatch(const CellSelection& selection) {
  if (selection.cells.empty()) {
    throw std::logic_error("fillAmrHydroGhostCells: empty source selection");
  }
  return sourcePatch(selection.cells.front());
}

[[nodiscard]] const AmrHydroSparseRemoteCell* findSparseRemoteCell(
    const AmrHydroSparseGhostSource& source,
    std::size_t patch_local_cell) {
  const auto it = std::lower_bound(
      source.cells.begin(), source.cells.end(), patch_local_cell,
      [](const AmrHydroSparseRemoteCell& cell, std::size_t offset) {
        return static_cast<std::size_t>(cell.patch_local_cell) < offset;
      });
  if (it == source.cells.end() ||
      static_cast<std::size_t>(it->patch_local_cell) != patch_local_cell) {
    return nullptr;
  }
  return &*it;
}

[[nodiscard]] double cellWidth(
    const hydro::HydroPatchGeometry& geometry,
    hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return geometry.cell_width_x_comoving;
    case hydro::HydroFaceAxis::kY:
      return geometry.cell_width_y_comoving;
    case hydro::HydroFaceAxis::kZ:
      return geometry.cell_width_z_comoving;
  }
  return 0.0;
}

[[nodiscard]] int axisIndex(hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return 0;
    case hydro::HydroFaceAxis::kY:
      return 1;
    case hydro::HydroFaceAxis::kZ:
      return 2;
  }
  return 0;
}

[[nodiscard]] double sideSign(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? -1.0 : 1.0;
}

[[nodiscard]] std::array<double, 3> patchCellWidths(const PatchDescriptor& patch) {
  return {
      patch.extent_comov[0] / static_cast<double>(patch.cell_dims[0]),
      patch.extent_comov[1] / static_cast<double>(patch.cell_dims[1]),
      patch.extent_comov[2] / static_cast<double>(patch.cell_dims[2]),
  };
}

[[nodiscard]] std::array<double, 3> cellCenter(
    const PatchDescriptor& patch,
    std::size_t cell) {
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t i = cell % nx;
  const std::size_t j = (cell / nx) % ny;
  const std::size_t k = cell / (nx * ny);
  const auto widths = patchCellWidths(patch);
  return {
      patch.origin_comov[0] + (static_cast<double>(i) + 0.5) * widths[0],
      patch.origin_comov[1] + (static_cast<double>(j) + 0.5) * widths[1],
      patch.origin_comov[2] + (static_cast<double>(k) + 0.5) * widths[2],
  };
}

[[nodiscard]] bool containsPoint(
    const PatchDescriptor& patch,
    const std::array<double, 3>& point) {
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const double lower = patch.origin_comov[axis] - k_geometry_tol;
    const double upper = patch.origin_comov[axis] + patch.extent_comov[axis] + k_geometry_tol;
    if (point[axis] < lower || point[axis] >= upper) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] std::size_t cellContainingPoint(
    const PatchDescriptor& patch,
    const std::array<double, 3>& point) {
  const auto widths = patchCellWidths(patch);
  std::array<std::size_t, 3> ijk{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    double coordinate = (point[axis] - patch.origin_comov[axis]) / widths[axis];
    coordinate = std::clamp(
        coordinate,
        0.0,
        static_cast<double>(patch.cell_dims[axis]) - k_geometry_tol);
    ijk[axis] = static_cast<std::size_t>(coordinate);
    if (ijk[axis] >= patch.cell_dims[axis]) {
      ijk[axis] = patch.cell_dims[axis] - 1U;
    }
  }
  return ijk[0] + static_cast<std::size_t>(patch.cell_dims[0]) *
                      (ijk[1] + static_cast<std::size_t>(patch.cell_dims[1]) * ijk[2]);
}

[[nodiscard]] internal::AmrPeriodicFrame frameForTarget(const AmrHydroGhostFillPatch& target) {
  const auto& extent = target.periodic_box_extent_comoving;
  return internal::makeAmrPeriodicFrame(
      extent[0] > 0.0 || extent[1] > 0.0 || extent[2] > 0.0, extent);
}

// Probe point of a ghost at distance `layer` cell widths from its owner along
// the outward normal, canonicalized into the periodic box (Part 3: wrapped
// neighbours across the box boundary resolve to their real patch).
[[nodiscard]] std::array<double, 3> ghostProbePoint(
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroGhostCell& ghost,
    std::uint8_t layer,
    const internal::AmrPeriodicFrame& frame) {
  std::array<double, 3> point = cellCenter(patch_geometry.patch, ghost.owner_real_cell);
  point[axisIndex(ghost.axis)] += sideSign(ghost.side) * static_cast<double>(std::max<std::uint8_t>(layer, 1U)) *
      cellWidth(patch_geometry.geometry, ghost.axis);
  if (ghost.transverse_axis < 3U) {
    // P1-HYDRO-02 edge stencil ghost: also displaced along the transverse axis.
    const auto widths = patchCellWidths(patch_geometry.patch);
    point[ghost.transverse_axis] += static_cast<double>(ghost.transverse_offset) * widths[ghost.transverse_axis];
  }
  return internal::wrapPoint(frame, point);
}


[[nodiscard]] std::size_t alignedGhostBoundaryIndex(
    double coordinate,
    double source_origin,
    double source_width,
    std::size_t source_dim,
    bool upper_boundary) {
  if (!(source_width > 0.0) || !std::isfinite(coordinate) ||
      !std::isfinite(source_origin) || !std::isfinite(source_width)) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: fine-to-coarse source geometry is non-finite or degenerate");
  }
  const double scaled = (coordinate - source_origin) / source_width;
  const double rounded = std::round(scaled);
  const double tolerance = 1.0e-9 * std::max(1.0, std::abs(scaled));
  if (std::abs(scaled - rounded) > tolerance) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: fine-to-coarse ghost volume is not aligned to source-cell boundaries");
  }
  if (rounded < 0.0 || rounded > static_cast<double>(source_dim)) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: fine source patch does not cover the complete coarse ghost volume");
  }
  const std::size_t index = static_cast<std::size_t>(rounded);
  if ((!upper_boundary && index >= source_dim) || (upper_boundary && index == 0U)) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: fine source range for coarse ghost is empty");
  }
  return index;
}

[[nodiscard]] std::array<std::array<double, 3>, 2> ghostCellBounds(
    const PatchDescriptor& target_patch,
    const hydro::HydroGhostCell& ghost,
    std::uint8_t layer,
    const internal::AmrPeriodicFrame& frame) {
  auto center = cellCenter(target_patch, ghost.owner_real_cell);
  const auto target_widths = patchCellWidths(target_patch);
  center[axisIndex(ghost.axis)] += sideSign(ghost.side) *
      static_cast<double>(std::max<std::uint8_t>(layer, 1U)) * target_widths[axisIndex(ghost.axis)];
  if (ghost.transverse_axis < 3U) {
    center[ghost.transverse_axis] +=
        static_cast<double>(ghost.transverse_offset) * target_widths[ghost.transverse_axis];
  }
  std::array<std::array<double, 3>, 2> bounds{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    bounds[0][axis] = center[axis] - 0.5 * target_widths[axis];
    bounds[1][axis] = center[axis] + 0.5 * target_widths[axis];
    // A ghost volume beyond a periodic box face is the wrapped volume of the
    // neighbour on the opposite side (one rigid shift by +-L).
    const double shift = internal::intervalWrapShift(frame, axis, bounds[0][axis], bounds[1][axis]);
    bounds[0][axis] += shift;
    bounds[1][axis] += shift;
  }
  return bounds;
}

template <class TVisit>
void forEachFineSourceCellIntersectingGhost(
    const PatchDescriptor& target_patch,
    const hydro::HydroGhostCell& ghost,
    std::uint8_t layer,
    const internal::AmrPeriodicFrame& frame,
    const PatchDescriptor& source_patch, TVisit&& visit) {
  if (source_patch.level <= target_patch.level) {
    return;
  }
  const auto bounds = ghostCellBounds(target_patch, ghost, layer, frame);
  const auto source_widths = patchCellWidths(source_patch);
  std::array<std::size_t, 3> begin{};
  std::array<std::size_t, 3> end{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const double source_lower = source_patch.origin_comov[axis];
    const double source_upper = source_lower + source_patch.extent_comov[axis];
    const double lower = std::max(bounds[0][axis], source_lower);
    const double upper = std::min(bounds[1][axis], source_upper);
    const double scale = std::max({1.0, std::abs(lower), std::abs(upper)});
    if (upper - lower <= k_geometry_tol * scale) {
      return;
    }
    const std::size_t dim = source_patch.cell_dims[axis];
    begin[axis] = alignedGhostBoundaryIndex(
        lower, source_lower, source_widths[axis], dim, false);
    end[axis] = alignedGhostBoundaryIndex(
        upper, source_lower, source_widths[axis], dim, true);
    if (end[axis] <= begin[axis]) {
      return;
    }
  }

  const std::size_t nx = source_patch.cell_dims[0];
  const std::size_t ny = source_patch.cell_dims[1];
  (void)core::checkedSizeProduct3(
      end[0] - begin[0], end[1] - begin[1], end[2] - begin[2],
      "fillAmrHydroGhostCells fine-source patch intersection count");
  for (std::size_t k = begin[2]; k < end[2]; ++k) {
    for (std::size_t j = begin[1]; j < end[1]; ++j) {
      for (std::size_t i = begin[0]; i < end[0]; ++i) {
        visit(i + nx * (j + ny * k));
      }
    }
  }
}

void validateAndOrderFineSourceCoverage(
    const PatchDescriptor& target_patch,
    const hydro::HydroGhostCell& ghost,
    std::uint8_t layer,
    const internal::AmrPeriodicFrame& frame,
    CellSelection& selection) {
  if (selection.cells.empty()) {
    return;
  }
  const PatchDescriptor& first_source = sourcePatch(selection.cells.front());
  if (first_source.level <= target_patch.level) {
    throw std::logic_error(
        "fillAmrHydroGhostCells: fine-source coverage contains a non-fine patch");
  }
  const auto target_widths = patchCellWidths(target_patch);
  const auto source_widths = patchCellWidths(first_source);
  std::array<std::size_t, 3> refinement_ratio{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (!(target_widths[axis] > 0.0) || !(source_widths[axis] > 0.0)) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: fine-to-coarse cell width is non-positive");
    }
    const double ratio = target_widths[axis] / source_widths[axis];
    const double rounded = std::round(ratio);
    const double tolerance = 1.0e-9 * std::max(1.0, std::abs(ratio));
    if (rounded < 1.0 || std::abs(ratio - rounded) > tolerance) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: coarse ghost width is not an integer fine-cell multiple");
    }
    refinement_ratio[axis] = core::checkedIntegralNarrow<std::size_t>(
        static_cast<std::uint64_t>(rounded),
        "fillAmrHydroGhostCells fine-to-coarse refinement ratio");
  }
  const std::size_t expected_count = core::checkedSizeProduct3(
      refinement_ratio[0], refinement_ratio[1], refinement_ratio[2],
      "fillAmrHydroGhostCells complete fine-source stencil size");
  const auto bounds = ghostCellBounds(target_patch, ghost, layer, frame);

  for (SelectedSourceCell& selected : selection.cells) {
    const PatchDescriptor& source = sourcePatch(selected);
    if (source.level != first_source.level) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: coarse ghost is covered by mixed fine AMR levels");
    }
    const auto widths = patchCellWidths(source);
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const double width_scale = std::max({1.0, std::abs(widths[axis]), std::abs(source_widths[axis])});
      if (std::abs(widths[axis] - source_widths[axis]) > k_geometry_tol * width_scale) {
        throw std::runtime_error(
            "fillAmrHydroGhostCells: coarse ghost is covered by inconsistent fine-cell widths");
      }
    }

    const auto center = cellCenter(source, selected.cell);
    std::array<std::size_t, 3> ordinal{};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const double first_center = bounds[0][axis] + 0.5 * source_widths[axis];
      const double scaled = (center[axis] - first_center) / source_widths[axis];
      const double rounded = std::round(scaled);
      const double tolerance = 1.0e-9 * std::max(1.0, std::abs(scaled));
      if (std::abs(scaled - rounded) > tolerance || rounded < 0.0 ||
          rounded >= static_cast<double>(refinement_ratio[axis])) {
        throw std::runtime_error(
            "fillAmrHydroGhostCells: fine source cell is outside the aligned coarse ghost stencil");
      }
      ordinal[axis] = static_cast<std::size_t>(rounded);
    }
    const std::size_t linear = ordinal[0] + refinement_ratio[0] *
        (ordinal[1] + refinement_ratio[1] * ordinal[2]);
    if (linear >= expected_count) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: fine source cell is outside the aligned coarse ghost stencil");
    }
    selected.fine_stencil_ordinal = linear;
  }

  // Sorting the already-required selection replaces a second coverage bitmap.
  // Contiguous ordinals prove both completeness and non-overlap, including
  // multi-patch fine coverage, without allocating the theoretical stencil.
  std::sort(
      selection.cells.begin(), selection.cells.end(),
      [](const SelectedSourceCell& lhs, const SelectedSourceCell& rhs) {
        return lhs.fine_stencil_ordinal < rhs.fine_stencil_ordinal;
      });
  for (std::size_t i = 1U; i < selection.cells.size(); ++i) {
    if (selection.cells[i - 1U].fine_stencil_ordinal ==
        selection.cells[i].fine_stencil_ordinal) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: fine source patches overlap within a coarse ghost stencil");
    }
  }
  if (selection.cells.size() != expected_count) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: fine source patches do not cover the complete coarse ghost volume");
  }
  for (std::size_t i = 0U; i < selection.cells.size(); ++i) {
    if (selection.cells[i].fine_stencil_ordinal != i) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: fine source patches do not cover the complete coarse ghost volume");
    }
  }
}

void validatePatchView(const AmrHydroGhostFillPatch& patch) {
  if (patch.geometry == nullptr || patch.conserved == nullptr) {
    throw std::invalid_argument("fillAmrHydroGhostCells: patch view must provide geometry and conserved storage");
  }
  if (patch.conserved->size() < patch.geometry->geometry.totalCellStorageCount()) {
    throw std::invalid_argument("fillAmrHydroGhostCells: conserved storage is missing AMR ghost rows");
  }
  if (!patch.available_real_cells.empty() &&
      patch.available_real_cells.size() != patch.geometry->geometry.cellCount()) {
    throw std::invalid_argument("fillAmrHydroGhostCells: sparse source availability does not match real-cell count");
  }
  if (!std::isfinite(patch.target_state_time_code) || !std::isfinite(patch.ghost_fill_time_code) ||
      !std::isfinite(patch.source_current_state_time_code)) {
    throw std::invalid_argument("fillAmrHydroGhostCells: ghost-fill times must be finite");
  }
}

[[nodiscard]] bool staleRemotePatch(const AmrHydroGhostFillPatch& patch) {
  return patch.residency == AmrGhostSourceResidency::kRemoteReadOnly &&
      patch.ghost_hydro_epoch != patch.expected_ghost_hydro_epoch;
}

[[nodiscard]] bool timesMatch(double lhs, double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= k_time_tol * scale;
}

[[nodiscard]] bool sourceCellAvailable(
    const AmrHydroGhostFillPatch& patch,
    std::size_t cell) {
  return patch.available_real_cells.empty() ||
      (cell < patch.available_real_cells.size() && patch.available_real_cells[cell] != 0U);
}

void requireSourceCellAvailable(
    const AmrHydroGhostFillPatch& patch,
    std::size_t cell) {
  if (!sourceCellAvailable(patch, cell)) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: sparse remote interface payload is missing a required source cell");
  }
}

// Fills `selection` in place (scratch-arena storage reserved by the caller).
// M01: it intentionally returns nothing; returning the pmr selection by value
// copied it through the default (ungoverned) memory resource on every ghost.
void selectSourceCells(
    const AmrHydroGhostFillPatch& target,
    const hydro::HydroGhostCell& ghost,
    std::uint8_t layer,
    AmrHydroBoundaryClass boundary_class,
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    CellSelection& selection,
    bool force_any_relation = false) {
  const internal::AmrPeriodicFrame frame = frameForTarget(target);
  const auto probe = ghostProbePoint(*target.geometry, ghost, layer, frame);
  // Deep stencil-only layers (layer > 1) and edge stencil ghosts may cross
  // into a patch of another level than the face-adjacent layer (2:1 balance
  // is a FACE property); any legal level relation is accepted.
  const bool any_relation = force_any_relation || layer > 1U || ghost.transverse_axis < 3U;
  const PatchDescriptor& target_patch = target.geometry->patch;
  selection.cells.clear();

  // Same-level and coarse-to-fine ghost fills resolve one source cell by the
  // ghost probe point. Fine-to-coarse is handled separately below because one
  // coarse ghost volume may span several fine child patches.
  for (const AmrHydroGhostFillPatch& candidate : patches) {
    if (candidate.geometry == nullptr || candidate.conserved == nullptr) {
      continue;
    }
    const PatchDescriptor& source_patch = candidate.geometry->patch;
    if (source_patch.patch_id == target_patch.patch_id ||
        !containsPoint(source_patch, probe)) {
      continue;
    }
    const bool same_level =
        source_patch.level == target_patch.level &&
        (any_relation || boundary_class == AmrHydroBoundaryClass::kSameLevel);
    const bool coarse_to_fine =
        source_patch.level < target_patch.level &&
        (any_relation || boundary_class == AmrHydroBoundaryClass::kCoarseFine);
    if (!same_level && !coarse_to_fine) {
      continue;
    }
    const std::size_t source_cell = cellContainingPoint(source_patch, probe);
    requireSourceCellAvailable(candidate, source_cell);
    selection.cells.push_back(SelectedSourceCell{
        .patch = &candidate,
        .cell = source_cell});
    return;
  }

  for (const AmrHydroSparseGhostSource& candidate : remote_sources) {
    if (candidate.patch == nullptr) {
      continue;
    }
    const PatchDescriptor& source_patch = *candidate.patch;
    if (source_patch.patch_id == target_patch.patch_id ||
        !containsPoint(source_patch, probe)) {
      continue;
    }
    const bool same_level =
        source_patch.level == target_patch.level &&
        (any_relation || boundary_class == AmrHydroBoundaryClass::kSameLevel);
    const bool coarse_to_fine =
        source_patch.level < target_patch.level &&
        (any_relation || boundary_class == AmrHydroBoundaryClass::kCoarseFine);
    if (!same_level && !coarse_to_fine) {
      continue;
    }
    const std::size_t source_cell = cellContainingPoint(source_patch, probe);
    if (findSparseRemoteCell(candidate, source_cell) == nullptr) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: sparse remote interface payload is missing a required source cell");
    }
    selection.cells.push_back(SelectedSourceCell{
        .sparse_remote = &candidate,
        .cell = source_cell});
    return;
  }

  if (!any_relation && boundary_class != AmrHydroBoundaryClass::kCoarseFine) {
    return;
  }

  for (const AmrHydroGhostFillPatch& candidate : patches) {
    if (candidate.geometry == nullptr || candidate.conserved == nullptr) {
      continue;
    }
    const PatchDescriptor& source_patch = candidate.geometry->patch;
    if (source_patch.patch_id == target_patch.patch_id ||
        source_patch.level <= target_patch.level) {
      continue;
    }
    forEachFineSourceCellIntersectingGhost(
        target_patch, ghost, layer, frame, source_patch, [&](std::size_t cell) {
      requireSourceCellAvailable(candidate, cell);
      selection.cells.push_back(SelectedSourceCell{
          .patch = &candidate,
          .cell = cell});
    });
  }

  for (const AmrHydroSparseGhostSource& candidate : remote_sources) {
    if (candidate.patch == nullptr) {
      continue;
    }
    const PatchDescriptor& source_patch = *candidate.patch;
    if (source_patch.patch_id == target_patch.patch_id ||
        source_patch.level <= target_patch.level) {
      continue;
    }
    forEachFineSourceCellIntersectingGhost(
        target_patch, ghost, layer, frame, source_patch, [&](std::size_t cell) {
      if (findSparseRemoteCell(candidate, cell) == nullptr) {
        throw std::runtime_error(
            "fillAmrHydroGhostCells: sparse remote fine-to-coarse payload is missing a required source cell");
      }
      selection.cells.push_back(SelectedSourceCell{
          .sparse_remote = &candidate,
          .cell = cell});
    });
  }

  validateAndOrderFineSourceCoverage(target_patch, ghost, layer, frame, selection);
  return;
}

[[nodiscard]] hydro::HydroConservedState averageSourceState(const CellSelection& selection) {
  hydro::HydroConservedState state;
  for (const SelectedSourceCell& source : selection.cells) {
    if (source.patch != nullptr) {
      state += source.patch->conserved->loadCell(source.cell);
    } else {
      const AmrHydroSparseRemoteCell* remote_cell =
          findSparseRemoteCell(*source.sparse_remote, source.cell);
      if (remote_cell == nullptr) {
        throw std::runtime_error(
            "fillAmrHydroGhostCells: sparse remote source disappeared during ghost averaging");
      }
      state += remote_cell->conserved;
    }
  }
  const double inv_count = 1.0 / static_cast<double>(selection.cells.size());
  return inv_count * state;
}

// ---------------------------------------------------------------------------
// P1-55-01 / H02: dual-energy auxiliary thermal lane of AMR ghosts.
//
// Contract. The auxiliary lane is the comoving internal-energy DENSITY rho*e
// (same units and volume convention as total_energy_density_comoving; e is
// the canonical SPECIFIC internal energy GasCellSidecar::internal_energy_code).
// Ghost cells are temporary derived views, so the auxiliary value of a ghost
// is formed by EXACTLY the operator that forms its conserved state:
//   * same level / coarse->fine injection: the value of the one source cell;
//   * fine->coarse restriction: the equal-volume mean of rho*e over the
//     complete fine stencil (the same mean as the conserved densities, i.e.
//     the cell-integrated thermal energy divided by the coarse volume);
//   * temporal coarse->fine: (1-alpha) (rho e)_start + alpha (rho e)_end from
//     the persistent temporal history, alpha identical to the conserved state;
//   * sparse MPI sources: rho_remote * e_remote of the owner's canonical state
//     (a remote recursive interval is interpolated in rho*e before division,
//     so this product equals the interpolated density);
//   * physical / deep-layer / edge closures: the scalar of the cell the
//     closure copies or mirrors (momentum parity does not act on rho*e).
// A ghost whose operator touches any source with unknown (NaN / non-positive)
// auxiliary energy is stored as NaN: explicitly invalid, never a finite but
// unrelated value. The solver then uses E - K for that row and counts it.
// The lane is written only when the TARGET SoA has it enabled; real cells of
// the target and of every local source must already hold rho*e (or NaN).
// ---------------------------------------------------------------------------
[[nodiscard]] bool validAuxiliary(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] double sourceAuxiliary(const SelectedSourceCell& source) {
  constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();
  if (source.patch != nullptr) {
    const hydro::HydroConservedStateSoa& soa = *source.patch->conserved;
    if (!soa.auxiliaryInternalEnergyEnabled()) {
      return k_nan;
    }
    const auto auxiliary = soa.auxiliaryInternalEnergyDensityComoving();
    return source.cell < auxiliary.size() && validAuxiliary(auxiliary[source.cell])
        ? auxiliary[source.cell] : k_nan;
  }
  const AmrHydroSparseRemoteCell* remote_cell =
      source.sparse_remote != nullptr ? findSparseRemoteCell(*source.sparse_remote, source.cell) : nullptr;
  if (remote_cell == nullptr) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: sparse remote source disappeared during auxiliary averaging");
  }
  const double rho_e = remote_cell->conserved.mass_density_comoving * remote_cell->internal_energy_code;
  return validAuxiliary(remote_cell->internal_energy_code) && validAuxiliary(rho_e) ? rho_e : k_nan;
}

// Equal-volume mean over exactly the cells averageSourceState() used.
[[nodiscard]] double averageSourceAuxiliary(const CellSelection& selection) {
  double sum = 0.0;
  for (const SelectedSourceCell& source : selection.cells) {
    const double value = sourceAuxiliary(source);
    if (!validAuxiliary(value)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    sum += value;
  }
  return sum / static_cast<double>(selection.cells.size());
}

[[nodiscard]] double cellAuxiliary(const hydro::HydroConservedStateSoa& soa, std::size_t cell) {
  if (!soa.auxiliaryInternalEnergyEnabled()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const auto auxiliary = soa.auxiliaryInternalEnergyDensityComoving();
  return cell < auxiliary.size() ? auxiliary[cell] : std::numeric_limits<double>::quiet_NaN();
}

void storeGhostAuxiliary(
    AmrHydroGhostFillPatch& target,
    std::size_t ghost_cell,
    double value,
    AmrHydroGhostFillDiagnostics& diagnostics) {
  if (!target.conserved->auxiliaryInternalEnergyEnabled()) {
    return;
  }
  std::span<double> auxiliary = target.conserved->auxiliaryInternalEnergyDensityComoving();
  if (ghost_cell >= auxiliary.size()) {
    throw std::logic_error("fillAmrHydroGhostCells: auxiliary ghost row is outside target storage");
  }
  if (validAuxiliary(value)) {
    auxiliary[ghost_cell] = value;
    ++diagnostics.auxiliary_internal_energy_ghosts_transferred;
  } else {
    auxiliary[ghost_cell] = std::numeric_limits<double>::quiet_NaN();
    ++diagnostics.auxiliary_internal_energy_ghosts_invalid;
  }
}

[[nodiscard]] bool selectionIsStale(const CellSelection& selection) {
  return std::any_of(
      selection.cells.begin(), selection.cells.end(),
      [](const SelectedSourceCell& source) {
        if (source.patch != nullptr) {
          return staleRemotePatch(*source.patch);
        }
        return source.sparse_remote != nullptr &&
            source.sparse_remote->ghost_hydro_epoch !=
                source.sparse_remote->expected_ghost_hydro_epoch;
      });
}

[[nodiscard]] double sourceTime(const SelectedSourceCell& source) {
  if (source.patch != nullptr) {
    return source.patch->source_current_state_time_code;
  }
  if (source.sparse_remote != nullptr) {
    return source.sparse_remote->source_current_state_time_code;
  }
  throw std::logic_error("fillAmrHydroGhostCells: selected source cell has no owner");
}

[[nodiscard]] bool selectionTimesMatch(
    const CellSelection& selection,
    double target_time_code) {
  return std::all_of(
      selection.cells.begin(), selection.cells.end(),
      [target_time_code](const SelectedSourceCell& source) {
        return timesMatch(sourceTime(source), target_time_code);
      });
}

[[nodiscard]] bool selectionIsSparseRemote(const CellSelection& selection) noexcept {
  return !selection.cells.empty() && selection.cells.front().sparse_remote != nullptr;
}

[[nodiscard]] const AmrHydroGhostFillPatch* selectionLocalPatch(
    const CellSelection& selection) noexcept {
  return selection.cells.empty() ? nullptr : selection.cells.front().patch;
}

[[nodiscard]] bool finiteConserved(const hydro::HydroConservedState& state) {
  return std::isfinite(state.mass_density_comoving) &&
      std::isfinite(state.momentum_density_x_comoving) &&
      std::isfinite(state.momentum_density_y_comoving) &&
      std::isfinite(state.momentum_density_z_comoving) &&
      std::isfinite(state.total_energy_density_comoving) &&
      std::isfinite(state.metal_mass_density_comoving);
}

[[nodiscard]] hydro::HydroConservedState temporalStateForCell(
    const AmrHydroGhostFillPatch& source,
    std::size_t source_cell,
    double fill_time_code,
    double adiabatic_index,
    bool& endpoint_fill,
    double& auxiliary,
    AmrHydroGhostFillDiagnostics& diagnostics) {
  auxiliary = std::numeric_limits<double>::quiet_NaN();
  if (source.temporal_boundary_history == nullptr) {
    ++diagnostics.temporal_history_missing_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse-to-fine temporal history is required but unavailable");
  }
  const auto* record = source.temporal_boundary_history->findByPatchId(source.geometry->patch.patch_id);
  if (record == nullptr) {
    ++diagnostics.temporal_history_missing_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: no temporal history record for coarse source patch");
  }
  if (!record->end_state_valid || record->patch_level != source.geometry->patch.level) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse temporal history is incomplete or level-invalid");
  }
  if (record->patch_geometry_fingerprint != amrPatchGeometryFingerprint(source.geometry->patch)) {
    ++diagnostics.temporal_geometry_mismatch_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse temporal history geometry fingerprint is stale");
  }
  if (record->gas_cell_identity_generation != source.geometry->source_gas_cell_identity_generation) {
    ++diagnostics.temporal_identity_mismatch_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse temporal history gas identity generation is stale");
  }
  const double duration = record->interval_end_code - record->interval_start_code;
  if (!std::isfinite(duration) || duration <= 0.0 || !std::isfinite(fill_time_code)) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporal history interval is invalid");
  }
  double alpha = (fill_time_code - record->interval_start_code) / duration;
  const double endpoint_tol = k_time_tol * std::max({1.0, std::abs(record->interval_start_code),
      std::abs(record->interval_end_code), std::abs(fill_time_code)});
  if (alpha < -endpoint_tol || alpha > 1.0 + endpoint_tol) {
    ++diagnostics.temporal_time_out_of_range_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: requested coarse temporal ghost lies outside validated interval");
  }
  if (alpha <= endpoint_tol) {
    alpha = 0.0;
    endpoint_fill = true;
  } else if (alpha >= 1.0 - endpoint_tol) {
    alpha = 1.0;
    endpoint_fill = true;
  }
  if (source_cell >= source.geometry->gas_cell_ids.size()) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporal source cell is outside source geometry");
  }
  const std::uint64_t gas_cell_id = source.geometry->gas_cell_ids[source_cell];
  const auto cell_it = std::find_if(
      record->cells.begin(), record->cells.end(),
      [gas_cell_id, source_cell](const core::AmrTemporalBoundaryHistoryCellRecord& cell) {
        return cell.gas_cell_id == gas_cell_id && cell.patch_local_cell == source_cell;
      });
  if (cell_it == record->cells.end()) {
    ++diagnostics.temporal_identity_mismatch_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporal history cannot resolve coarse source by stable gas_cell_id");
  }
  const hydro::HydroConservedState start{
      .mass_density_comoving = cell_it->start_mass_density_comoving,
      .momentum_density_x_comoving = cell_it->start_momentum_density_x_comoving,
      .momentum_density_y_comoving = cell_it->start_momentum_density_y_comoving,
      .momentum_density_z_comoving = cell_it->start_momentum_density_z_comoving,
      .total_energy_density_comoving = cell_it->start_total_energy_density_comoving,
      .metal_mass_density_comoving = cell_it->start_metal_mass_density_comoving};
  const hydro::HydroConservedState end{
      .mass_density_comoving = cell_it->end_mass_density_comoving,
      .momentum_density_x_comoving = cell_it->end_momentum_density_x_comoving,
      .momentum_density_y_comoving = cell_it->end_momentum_density_y_comoving,
      .momentum_density_z_comoving = cell_it->end_momentum_density_z_comoving,
      .total_energy_density_comoving = cell_it->end_total_energy_density_comoving,
      .metal_mass_density_comoving = cell_it->end_metal_mass_density_comoving};
  const hydro::HydroConservedState result = (1.0 - alpha) * start + alpha * end;
  if (!finiteConserved(result) || result.mass_density_comoving <= 0.0) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporally interpolated conserved state is non-finite or non-positive");
  }
  const hydro::HydroPrimitiveState primitive = hydro::HydroCoreSolver::primitiveFromConserved(result, adiabatic_index);
  if (!std::isfinite(primitive.rho_comoving) || !std::isfinite(primitive.pressure_comoving) ||
      !std::isfinite(primitive.vel_x_peculiar) || !std::isfinite(primitive.vel_y_peculiar) ||
      !std::isfinite(primitive.vel_z_peculiar) || primitive.rho_comoving <= 0.0 ||
      primitive.pressure_comoving <= 0.0) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: primitive recovery of temporal conserved state is inadmissible");
  }
  // H02: the auxiliary thermal lane uses the same endpoints and alpha. A
  // history captured without the lane (NaN; e.g. a pre-v25 checkpoint) keeps
  // the ghost explicitly invalid instead of borrowing a current-time value.
  if (validAuxiliary(cell_it->start_internal_energy_density_comoving) &&
      validAuxiliary(cell_it->end_internal_energy_density_comoving)) {
    auxiliary = (1.0 - alpha) * cell_it->start_internal_energy_density_comoving +
        alpha * cell_it->end_internal_energy_density_comoving;
  }
  return result;
}

[[noreturn]] void rejectTemporal(
    AmrHydroGhostDescriptor& descriptor,
    AmrHydroGhostFillStatus status,
    std::size_t& counter,
    AmrHydroGhostFillDiagnostics& diagnostics,
    std::string_view message) {
  descriptor.fill_status = status;
  ++counter;
  ++diagnostics.unresolved_ghosts;
  throw std::runtime_error(std::string("fillAmrHydroGhostCells: ") + std::string(message));
}

// ---------------------------------------------------------------------------
// P1-55-02 / P2-55-03: domain-versus-missing-data decision for stencil-only
// ghosts (deep normal layers, edge/corner stencil ghosts, and deep physical
// reflective layers of patches thinner than the stencil).
//
// Every such ghost is described by the target-relative integer cell index of
// the stencil cell it stands for (owner index + layer along the normal, +
// transverse offset for edge ghosts). For each axis on which that index
// leaves the TARGET patch the global boundary is classified:
//   * periodic axis: no closure; the probe wraps through the periodic frame
//     and must resolve to a real (local, coarse/fine or MPI) source;
//   * non-periodic axis with a declared global domain
//     (AmrHydroGhostFillPatch::domain_extent_comoving): the domain walls are
//     the aligned domain faces in target-cell units;
//   * non-periodic axis without a declared domain: only a target face that is
//     itself classified kPhysical (no neighbour patch) is a wall. An
//     interface face never becomes a wall by absence of data.
// Index components beyond a wall are mapped with the physical closure:
// reflective -> mirror i' = 2 w - 1 - i (normal momentum of THAT axis
// reversed; repeated for a domain narrower than the stencil), open ->
// zero-gradient onto the last in-domain row. The mapped in-domain stencil
// cell is then resolved like any other stencil cell: a real target cell
// (including an in-patch periodic wrap), or a probe through the ordinary
// local / coarse->fine (temporal) / fine->coarse / sparse MPI selection with
// any legal level relation. Missing data at the mapped location fails closed.
// A deep normal layer whose stencil cell is provably (or not provably out of)
// inside the domain and has no source is an in-domain missing source and is
// rejected deterministically: it is never continued from a shallower ghost.
// Closures never create faces or flux registers (ghost-only semantics).
// ---------------------------------------------------------------------------
constexpr int k_max_closure_reflections_per_axis = 8;

[[nodiscard]] std::size_t standardShellGhostSlot(
    const hydro::HydroPatchGeometry& geometry, std::size_t normal_axis, std::size_t side,
    const std::array<std::size_t, 3>& ijk) {
  const std::size_t nx = geometry.nx;
  const std::size_t ny = geometry.ny;
  const std::size_t nz = geometry.nz;
  switch (normal_axis) {
    case 0U: return 2U * (ijk[2] * ny + ijk[1]) + side;
    case 1U: return 2U * ny * nz + 2U * (ijk[2] * nx + ijk[0]) + side;
    default: return 2U * (ny * nz + nx * nz) + 2U * (ijk[1] * nx + ijk[0]) + side;
  }
}

void reverseMomentumComponent(hydro::HydroConservedState& state, std::size_t axis) {
  switch (axis) {
    case 0U: state.momentum_density_x_comoving = -state.momentum_density_x_comoving; break;
    case 1U: state.momentum_density_y_comoving = -state.momentum_density_y_comoving; break;
    default: state.momentum_density_z_comoving = -state.momentum_density_z_comoving; break;
  }
}

[[nodiscard]] hydro::HydroFaceAxis faceAxisFromIndex(std::size_t axis) {
  switch (axis) {
    case 0U: return hydro::HydroFaceAxis::kX;
    case 1U: return hydro::HydroFaceAxis::kY;
    default: return hydro::HydroFaceAxis::kZ;
  }
}

[[nodiscard]] std::string ghostLocationText(
    const AmrHydroGhostFillPatch& target,
    const AmrHydroGhostDescriptor& descriptor) {
  const hydro::HydroGhostCell& ghost = target.geometry->geometry.ghost_cells.at(descriptor.ghost_slot);
  std::string text = "patch " + std::to_string(target.geometry->patch.patch_id) +
      " level " + std::to_string(static_cast<unsigned>(target.geometry->patch.level)) +
      " axis " + std::to_string(axisIndex(ghost.axis)) +
      (ghost.side == hydro::HydroFaceSide::kLower ? " side lower" : " side upper") +
      " layer " + std::to_string(static_cast<unsigned>(descriptor.ghost_layer)) +
      " epoch " + std::to_string(target.ghost_hydro_epoch);
  if (ghost.transverse_axis < 3U) {
    text += " transverse axis " + std::to_string(static_cast<unsigned>(ghost.transverse_axis)) +
        " offset " + std::to_string(static_cast<int>(ghost.transverse_offset));
  }
  return text;
}

// The face-adjacent standard-shell ghost of the target face (axis, side) is
// a non-periodic physical ghost: the target owns no neighbour on that face.
[[nodiscard]] bool targetFaceIsPhysicalWall(
    const AmrHydroGhostFillPatch& target, std::size_t axis, std::size_t side) {
  const hydro::HydroPatchGeometry& geometry = target.geometry->geometry;
  const std::size_t slot = standardShellGhostSlot(geometry, axis, side, {0U, 0U, 0U});
  const AmrHydroGhostDescriptor& face = target.geometry->ghosts.at(slot);
  const hydro::HydroGhostCell& face_ghost = geometry.ghost_cells.at(slot);
  if (face.edge_stencil || face.ghost_layer != 1U ||
      static_cast<std::size_t>(axisIndex(face_ghost.axis)) != axis ||
      (face_ghost.side == hydro::HydroFaceSide::kLower ? 0U : 1U) != side) {
    throw std::logic_error("fillAmrHydroGhostCells: standard ghost shell layout mismatch");
  }
  return face.boundary_class == AmrHydroBoundaryClass::kPhysical &&
      face_ghost.boundary_kind != hydro::HydroBoundaryKind::kPeriodic;
}

[[nodiscard]] bool domainDeclared(const AmrHydroGhostFillPatch& target, std::size_t axis) noexcept {
  const double extent = target.domain_extent_comoving[axis];
  return extent > 0.0 && std::isfinite(extent) && std::isfinite(target.domain_origin_comoving[axis]);
}

[[nodiscard]] std::ptrdiff_t alignedDomainIndex(double coordinate, double patch_origin, double width) {
  const double scaled = (coordinate - patch_origin) / width;
  const double rounded = std::round(scaled);
  if (!std::isfinite(scaled) || !(width > 0.0) ||
      std::abs(scaled - rounded) > 1.0e-6 * std::max(1.0, std::abs(scaled))) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: patch cells are not aligned to the declared global domain faces");
  }
  return static_cast<std::ptrdiff_t>(rounded);
}

struct ClosureAxisBounds {
  bool periodic = false;
  bool has_lower = false;
  bool has_upper = false;
  std::ptrdiff_t lower = 0;  // first in-domain target-relative cell index
  std::ptrdiff_t upper = 0;  // one past the last in-domain index
};

[[nodiscard]] ClosureAxisBounds closureBounds(const AmrHydroGhostFillPatch& target, std::size_t axis) {
  ClosureAxisBounds bounds;
  if (frameForTarget(target).axisPeriodic(axis)) {
    bounds.periodic = true;
    return bounds;
  }
  const PatchDescriptor& patch = target.geometry->patch;
  const auto n = static_cast<std::ptrdiff_t>(patch.cell_dims[axis]);
  if (domainDeclared(target, axis)) {
    const double width = patch.extent_comov[axis] / static_cast<double>(patch.cell_dims[axis]);
    bounds.has_lower = true;
    bounds.has_upper = true;
    bounds.lower = alignedDomainIndex(target.domain_origin_comoving[axis], patch.origin_comov[axis], width);
    bounds.upper = alignedDomainIndex(
        target.domain_origin_comoving[axis] + target.domain_extent_comoving[axis], patch.origin_comov[axis], width);
    if (bounds.lower > 0 || bounds.upper < n) {
      throw std::runtime_error(
          "fillAmrHydroGhostCells: patch " + std::to_string(patch.patch_id) +
          " extends outside the declared global domain");
    }
    return bounds;
  }
  if (targetFaceIsPhysicalWall(target, axis, 0U)) {
    bounds.has_lower = true;
    bounds.lower = 0;
  }
  if (targetFaceIsPhysicalWall(target, axis, 1U)) {
    bounds.has_upper = true;
    bounds.upper = n;
  }
  return bounds;
}

// A face classified physical (no neighbour in the local + remote descriptor
// set) must lie on a declared non-periodic domain wall. Otherwise neighbour
// metadata is missing (e.g. an underdelivered MPI interface) and the face
// would silently receive an outflow/reflective closure inside the domain.
void validateTargetPhysicalFacesOnDomainWalls(
    const AmrHydroGhostFillPatch& target, AmrHydroGhostFillDiagnostics& diagnostics) {
  const internal::AmrPeriodicFrame frame = frameForTarget(target);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (frame.axisPeriodic(axis) || !domainDeclared(target, axis)) continue;
    const ClosureAxisBounds bounds = closureBounds(target, axis);
    const auto n = static_cast<std::ptrdiff_t>(target.geometry->patch.cell_dims[axis]);
    for (std::size_t side = 0; side < 2U; ++side) {
      if (!targetFaceIsPhysicalWall(target, axis, side)) continue;
      const bool on_wall = side == 0U ? bounds.lower == 0 : bounds.upper == n;
      if (!on_wall) {
        ++diagnostics.physical_face_inside_domain_rejections;
        ++diagnostics.unresolved_ghosts;
        throw std::runtime_error(
            "fillAmrHydroGhostCells: patch " + std::to_string(target.geometry->patch.patch_id) +
            " level " + std::to_string(static_cast<unsigned>(target.geometry->patch.level)) +
            " axis " + std::to_string(axis) + (side == 0U ? " lower" : " upper") +
            " face has no neighbour but lies inside the global domain (missing interface source)");
      }
    }
  }
}

struct ClosureOutcome {
  bool closed = false;     // at least one axis crossed a physical wall
  bool reflected = false;  // at least one reflective mirror was applied
};

void fillAmrGhost(
    AmrHydroGhostFillPatch& target,
    AmrHydroGhostDescriptor& descriptor,
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index,
    AmrHydroGhostFillDiagnostics& diagnostics,
    CellSelection& selection,
    const hydro::HydroGhostCell* probe_ghost = nullptr,
    std::uint8_t probe_layer = 0U,
    bool mapped_probe = false);

// require_wall: the stencil cell is known to have no in-domain source at its
// unmapped location (deep normal layer with an empty selection, or a deep
// physical layer beyond the target wall). If no wall is crossed it is an
// in-domain missing source and is rejected.
[[nodiscard]] ClosureOutcome fillClosureMappedGhost(
    AmrHydroGhostFillPatch& target,
    AmrHydroGhostDescriptor& descriptor,
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index,
    AmrHydroGhostFillDiagnostics& diagnostics,
    CellSelection& selection,
    std::array<std::ptrdiff_t, 3> index,
    bool require_wall) {
  const hydro::HydroPatchGeometry& geometry = target.geometry->geometry;
  const hydro::HydroGhostCell& ghost = geometry.ghost_cells.at(descriptor.ghost_slot);
  const std::array<std::ptrdiff_t, 3> dims{
      static_cast<std::ptrdiff_t>(geometry.nx), static_cast<std::ptrdiff_t>(geometry.ny),
      static_cast<std::ptrdiff_t>(geometry.nz)};
  const auto rejectUnsupported = [&](std::string_view why) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedUnsupportedClosureGeometry;
    ++diagnostics.unsupported_closure_geometry_rejections;
    ++diagnostics.unresolved_ghosts;
    throw std::runtime_error("fillAmrHydroGhostCells: unsupported physical closure geometry (" +
                             std::string(why) + ") for " + ghostLocationText(target, descriptor));
  };
  ClosureOutcome outcome;
  std::array<bool, 3> reverse{false, false, false};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (index[axis] >= 0 && index[axis] < dims[axis]) continue;  // no target face crossed
    const ClosureAxisBounds bounds = closureBounds(target, axis);
    if (bounds.periodic) continue;
    for (int pass = 0;; ++pass) {
      const bool below = bounds.has_lower && index[axis] < bounds.lower;
      const bool above = bounds.has_upper && index[axis] >= bounds.upper;
      if (!below && !above) break;
      if (pass >= k_max_closure_reflections_per_axis) {
        rejectUnsupported("global domain narrower than the reconstruction stencil");
      }
      outcome.closed = true;
      if (target.physical_boundary_kind == hydro::HydroBoundaryKind::kReflective) {
        index[axis] = below ? 2 * bounds.lower - 1 - index[axis] : 2 * bounds.upper - 1 - index[axis];
        reverse[axis] = !reverse[axis];
        outcome.reflected = true;
      } else if (target.physical_boundary_kind == hydro::HydroBoundaryKind::kOpen) {
        index[axis] = below ? bounds.lower : bounds.upper - 1;
      } else {
        rejectUnsupported("hierarchy boundary kind has no stencil closure");
      }
    }
  }
  if (require_wall && !outcome.closed) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kMissingSource;
    ++diagnostics.missing_source_records;
    ++diagnostics.in_domain_missing_source_rejections;
    ++diagnostics.unresolved_ghosts;
    throw std::runtime_error(
        "fillAmrHydroGhostCells: in-domain stencil ghost has no local, coarse/fine or MPI source (" +
        ghostLocationText(target, descriptor) + "); refusing to extrapolate missing exchange data");
  }

  // Resolve the mapped in-domain stencil cell.
  std::array<std::size_t, 3> clamped{};
  std::array<std::size_t, 2> out_axis{};
  std::array<std::ptrdiff_t, 2> out_offset{};  // signed distance beyond the target face
  std::size_t out_count = 0U;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (index[axis] < 0) {
      clamped[axis] = 0U;
      out_axis.at(out_count) = axis;
      out_offset.at(out_count++) = index[axis];
    } else if (index[axis] >= dims[axis]) {
      clamped[axis] = static_cast<std::size_t>(dims[axis] - 1);
      out_axis.at(out_count) = axis;
      out_offset.at(out_count++) = index[axis] - (dims[axis] - 1);
    } else {
      clamped[axis] = static_cast<std::size_t>(index[axis]);
    }
  }
  const auto applyParity = [&](hydro::HydroConservedState state) {
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (reverse[axis]) reverseMomentumComponent(state, axis);
    }
    return state;
  };
  std::size_t target_cell = hydro::k_invalid_cell_index;
  hydro::HydroGhostCell probe = ghost;
  std::uint8_t probe_layer = 1U;
  if (out_count == 0U) {
    target_cell = geometry.linearCellIndex(clamped[0], clamped[1], clamped[2]);
  } else {
    const std::ptrdiff_t layer_distance = out_offset[0] < 0 ? -out_offset[0] : out_offset[0];
    if (layer_distance > std::numeric_limits<std::uint8_t>::max() ||
        (out_count == 2U && (out_offset[1] > std::numeric_limits<std::int8_t>::max() ||
                             out_offset[1] < -std::numeric_limits<std::int8_t>::max()))) {
      rejectUnsupported("mapped stencil distance exceeds the ghost probe range");
    }
    probe.owner_real_cell = geometry.linearCellIndex(clamped[0], clamped[1], clamped[2]);
    probe.axis = faceAxisFromIndex(out_axis[0]);
    probe.side = out_offset[0] < 0 ? hydro::HydroFaceSide::kLower : hydro::HydroFaceSide::kUpper;
    probe_layer = static_cast<std::uint8_t>(layer_distance);
    probe.transverse_axis = out_count == 2U ? static_cast<std::uint8_t>(out_axis[1]) : std::uint8_t{3U};
    probe.transverse_offset = out_count == 2U ? static_cast<std::int8_t>(out_offset[1]) : std::int8_t{0};
    const auto point = ghostProbePoint(*target.geometry, probe, probe_layer, frameForTarget(target));
    if (containsPoint(target.geometry->patch, point)) {
      // In-patch periodic wrap: the stencil cell is a real cell of this patch.
      target_cell = cellContainingPoint(target.geometry->patch, point);
    }
  }
  if (target_cell != hydro::k_invalid_cell_index) {
    target.conserved->storeCell(descriptor.ghost_cell, applyParity(target.conserved->loadCell(target_cell)));
    storeGhostAuxiliary(target, descriptor.ghost_cell, cellAuxiliary(*target.conserved, target_cell), diagnostics);
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledSameLevel;
    selection.cells.clear();
    if (descriptor.edge_stencil && descriptor.source_gas_cell_id == 0U &&
        target_cell < target.geometry->gas_cell_ids.size()) {
      descriptor.source_gas_cell_id = target.geometry->gas_cell_ids[target_cell];
    }
  } else {
    fillAmrGhost(target, descriptor, patches, remote_sources, adiabatic_index, diagnostics, selection,
                 &probe, probe_layer, /*mapped_probe=*/true);
    if (reverse[0] || reverse[1] || reverse[2]) {
      target.conserved->storeCell(descriptor.ghost_cell,
                                  applyParity(target.conserved->loadCell(descriptor.ghost_cell)));
    }
    // Stable source identity only for a single-cell source (composition
    // origin / register bookkeeping). The auxiliary lane never uses it.
    if (descriptor.edge_stencil && descriptor.source_gas_cell_id == 0U && selection.cells.size() == 1U) {
      const SelectedSourceCell& source = selection.cells.front();
      if (source.patch != nullptr && source.cell < source.patch->geometry->gas_cell_ids.size()) {
        descriptor.source_gas_cell_id = source.patch->geometry->gas_cell_ids[source.cell];
      } else if (source.sparse_remote != nullptr) {
        const AmrHydroSparseRemoteCell* remote = findSparseRemoteCell(*source.sparse_remote, source.cell);
        if (remote != nullptr) descriptor.source_gas_cell_id = remote->gas_cell_id;
      }
    }
  }
  if (!descriptor.edge_stencil && outcome.closed) {
    descriptor.fill_status = outcome.reflected ? AmrHydroGhostFillStatus::kFilledDeepLayerReflected
                                               : AmrHydroGhostFillStatus::kFilledDeepLayerExtrapolated;
    if (outcome.reflected) {
      ++diagnostics.deep_layer_reflected_ghosts;
    } else {
      ++diagnostics.deep_layer_extrapolated_ghosts;
    }
  }
  return outcome;
}

[[nodiscard]] std::array<std::ptrdiff_t, 3> deepLayerStencilIndex(
    const hydro::HydroPatchGeometry& geometry,
    const hydro::HydroGhostCell& ghost,
    std::uint8_t layer) {
  const std::array<std::size_t, 3> owner = geometry.cellIjk(ghost.owner_real_cell);
  std::array<std::ptrdiff_t, 3> index{
      static_cast<std::ptrdiff_t>(owner[0]), static_cast<std::ptrdiff_t>(owner[1]),
      static_cast<std::ptrdiff_t>(owner[2])};
  index[static_cast<std::size_t>(axisIndex(ghost.axis))] +=
      (ghost.side == hydro::HydroFaceSide::kLower ? -1 : 1) * static_cast<std::ptrdiff_t>(layer);
  if (ghost.transverse_axis < 3U) {
    index[ghost.transverse_axis] += static_cast<std::ptrdiff_t>(ghost.transverse_offset);
  }
  return index;
}

void fillAmrGhost(
    AmrHydroGhostFillPatch& target,
    AmrHydroGhostDescriptor& descriptor,
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index,
    AmrHydroGhostFillDiagnostics& diagnostics,
    CellSelection& selection,
    const hydro::HydroGhostCell* probe_ghost,
    std::uint8_t probe_layer,
    bool mapped_probe) {
  const hydro::HydroGhostCell& ghost = probe_ghost != nullptr
      ? *probe_ghost
      : target.geometry->geometry.ghost_cells.at(descriptor.ghost_slot);
  const std::uint8_t layer = probe_ghost != nullptr ? probe_layer : descriptor.ghost_layer;
  selectSourceCells(target, ghost, layer, descriptor.boundary_class,
                    patches, remote_sources, selection, descriptor.edge_stencil || mapped_probe);
  if (selection.cells.empty() && probe_ghost == nullptr && descriptor.ghost_layer > 1U &&
      !descriptor.edge_stencil) {
    // P1-55-02: a deep stencil-only layer without a source is legal ONLY
    // beyond a proven non-periodic physical wall; otherwise it is rejected.
    (void)fillClosureMappedGhost(
        target, descriptor, patches, remote_sources, adiabatic_index, diagnostics, selection,
        deepLayerStencilIndex(target.geometry->geometry, ghost, descriptor.ghost_layer),
        /*require_wall=*/true);
    return;
  }
  if (selection.cells.empty()) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kMissingSource;
    ++diagnostics.missing_source_records;
    ++diagnostics.unresolved_ghosts;
    if (mapped_probe) ++diagnostics.in_domain_missing_source_rejections;
    throw std::runtime_error(
        std::string("fillAmrHydroGhostCells: unresolved AMR hydro ghost (") +
        (mapped_probe ? "closure-mapped in-domain stencil cell, " : "") +
        ghostLocationText(target, descriptor) + ")");
  }
  if (selectionIsStale(selection)) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedStaleRemote;
    ++diagnostics.stale_epoch_rejections;
    ++diagnostics.unresolved_ghosts;
    throw std::runtime_error("fillAmrHydroGhostCells: stale source AMR hydro ghost epoch");
  }

  const PatchDescriptor& source_patch = selectionPatch(selection);
  const bool same_level = source_patch.level == target.geometry->patch.level;
  const bool coarse_to_fine = source_patch.level < target.geometry->patch.level;
  const bool fine_to_coarse = source_patch.level > target.geometry->patch.level;
  if (same_level) {
    if (!selectionTimesMatch(selection, target.ghost_fill_time_code)) {
      rejectTemporal(
          descriptor,
          AmrHydroGhostFillStatus::kRejectedTemporalSameLevelMismatch,
          diagnostics.temporal_same_level_mismatch_rejections,
          diagnostics,
          "same-level source state time does not match requested ghost-fill time");
    }
    target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
    storeGhostAuxiliary(target, descriptor.ghost_cell, averageSourceAuxiliary(selection), diagnostics);
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledSameLevel;
    ++diagnostics.same_level_ghosts_filled;
    return;
  }
  if (fine_to_coarse) {
    if (!selectionTimesMatch(selection, target.ghost_fill_time_code)) {
      rejectTemporal(
          descriptor,
          AmrHydroGhostFillStatus::kRejectedTemporalFineToCoarseMismatch,
          diagnostics.temporal_fine_to_coarse_misalignment_rejections,
          diagnostics,
          "fine-to-coarse ghost access is not at a synchronization time");
    }
    target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
    storeGhostAuxiliary(target, descriptor.ghost_cell, averageSourceAuxiliary(selection), diagnostics);
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledFineToCoarse;
    ++diagnostics.fine_to_coarse_ghosts_filled;
    return;
  }

  if (coarse_to_fine && target.enable_temporal_coarse_to_fine && selectionIsSparseRemote(selection)) {
    // G1.1: a sparse remote coarse source is accepted ONLY when its state is
    // labelled with exactly the requested fill time, i.e. the distributed
    // recursive driver already interpolated it linearly in conserved
    // variables between the exchanged start/end states of the coarse
    // interval (the same law as temporalStateForCell). A current (end-of-
    // interval) remote coarse state carries its own clock and is rejected.
    if (!selectionTimesMatch(selection, target.ghost_fill_time_code)) {
      rejectTemporal(
          descriptor,
          AmrHydroGhostFillStatus::kRejectedTemporalHistoryMissing,
          diagnostics.temporal_history_missing_rejections,
          diagnostics,
          "remote coarse ghost is not temporally interpolated to the fill time");
    }
    target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
    storeGhostAuxiliary(target, descriptor.ghost_cell, averageSourceAuxiliary(selection), diagnostics);
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledCoarseToFineTemporal;
    ++diagnostics.coarse_to_fine_ghosts_filled;
    ++diagnostics.temporal_coarse_to_fine_ghosts_filled;
    return;
  }
  if (coarse_to_fine && target.enable_temporal_coarse_to_fine) {
    const AmrHydroGhostFillPatch* temporal_source = selectionLocalPatch(selection);
    if (selectionIsSparseRemote(selection) || temporal_source == nullptr ||
        (temporal_source->residency == AmrGhostSourceResidency::kRemoteReadOnly &&
         temporal_source->temporal_boundary_history == nullptr)) {
      rejectTemporal(
          descriptor,
          AmrHydroGhostFillStatus::kRejectedTemporalHistoryMissing,
          diagnostics.temporal_history_missing_rejections,
          diagnostics,
          "remote temporal interpolation requires an explicitly cached local temporal history");
    }
    hydro::HydroConservedState accumulated;
    double accumulated_auxiliary = 0.0;
    bool auxiliary_valid = true;
    bool endpoint_fill = false;
    const std::size_t missing_before = diagnostics.temporal_history_missing_rejections;
    const std::size_t invalid_before = diagnostics.temporal_history_invalid_rejections;
    const std::size_t out_of_range_before = diagnostics.temporal_time_out_of_range_rejections;
    const std::size_t geometry_before = diagnostics.temporal_geometry_mismatch_rejections;
    const std::size_t identity_before = diagnostics.temporal_identity_mismatch_rejections;
    try {
      for (const SelectedSourceCell& source : selection.cells) {
        if (source.patch == nullptr) {
          throw std::logic_error(
              "fillAmrHydroGhostCells: temporal coarse source unexpectedly uses sparse remote state");
        }
        bool cell_endpoint = false;
        double cell_auxiliary = std::numeric_limits<double>::quiet_NaN();
        accumulated += temporalStateForCell(
            *source.patch,
            source.cell,
            target.ghost_fill_time_code,
            adiabatic_index,
            cell_endpoint,
            cell_auxiliary,
            diagnostics);
        endpoint_fill = endpoint_fill || cell_endpoint;
        auxiliary_valid = auxiliary_valid && validAuxiliary(cell_auxiliary);
        if (auxiliary_valid) accumulated_auxiliary += cell_auxiliary;
      }
    } catch (const std::runtime_error&) {
      if (diagnostics.temporal_time_out_of_range_rejections > out_of_range_before) {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalTimeOutOfRange;
      } else if (diagnostics.temporal_history_missing_rejections > missing_before) {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalHistoryMissing;
      } else if (diagnostics.temporal_history_invalid_rejections > invalid_before ||
                 diagnostics.temporal_geometry_mismatch_rejections > geometry_before ||
                 diagnostics.temporal_identity_mismatch_rejections > identity_before) {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalHistoryInvalid;
      } else {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalHistoryInvalid;
      }
      ++diagnostics.unresolved_ghosts;
      throw;
    }
    const double inv_count = 1.0 / static_cast<double>(selection.cells.size());
    target.conserved->storeCell(descriptor.ghost_cell, inv_count * accumulated);
    storeGhostAuxiliary(target, descriptor.ghost_cell,
                        auxiliary_valid ? inv_count * accumulated_auxiliary
                                        : std::numeric_limits<double>::quiet_NaN(),
                        diagnostics);
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledCoarseToFineTemporal;
    ++diagnostics.coarse_to_fine_ghosts_filled;
    ++diagnostics.temporal_coarse_to_fine_ghosts_filled;
    if (endpoint_fill) {
      ++diagnostics.temporal_endpoint_ghosts_filled;
    }
    return;
  }

  if (coarse_to_fine && !selectionTimesMatch(selection, target.ghost_fill_time_code)) {
    rejectTemporal(
        descriptor,
        AmrHydroGhostFillStatus::kRejectedTemporalHistoryMissing,
        diagnostics.temporal_history_missing_rejections,
        diagnostics,
        "coarse-to-fine ghost fill is time-misaligned and lacks temporal interpolation");
  }
  target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
  storeGhostAuxiliary(target, descriptor.ghost_cell, averageSourceAuxiliary(selection), diagnostics);
  descriptor.fill_status = AmrHydroGhostFillStatus::kFilledCoarseToFine;
  ++diagnostics.coarse_to_fine_ghosts_filled;
}

// P1-HYDRO-02 edge stencil ghost: the stencil cell one layer beyond the
// owner's normal face (axis n) and |offset| cells beyond its transverse face
// (axis t). It is resolved through the closure mapper above: crossed
// non-periodic physical walls apply the physical closure per axis (open copy,
// reflective mirror with that axis' momentum reversed, including mirrors that
// land beyond the opposite face of a thin patch, which are then resolved from
// the neighbour patch there); the remaining in-domain location is a target
// real cell or a probe resolved with any legal level relation. Missing data
// fails closed; nothing is extrapolated and no degeneracy is substituted.
void fillEdgeStencilGhost(
    AmrHydroGhostFillPatch& target,
    AmrHydroGhostDescriptor& descriptor,
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index,
    AmrHydroGhostFillDiagnostics& diagnostics,
    CellSelection& selection) {
  const hydro::HydroPatchGeometry& geometry = target.geometry->geometry;
  const hydro::HydroGhostCell& ghost = geometry.ghost_cells.at(descriptor.ghost_slot);
  const std::size_t n = static_cast<std::size_t>(axisIndex(ghost.axis));
  if (ghost.transverse_axis > 2U || ghost.transverse_axis == n || ghost.transverse_offset == 0) {
    throw std::logic_error("fillAmrHydroGhostCells: malformed edge stencil ghost metadata");
  }
  const ClosureOutcome outcome = fillClosureMappedGhost(
      target, descriptor, patches, remote_sources, adiabatic_index, diagnostics, selection,
      deepLayerStencilIndex(geometry, ghost, 1U), /*require_wall=*/false);
  ++diagnostics.edge_stencil_ghosts_filled;
  if (outcome.closed) ++diagnostics.edge_stencil_physical_closures;
}

void markStaleRemoteTarget(
    AmrHydroGhostFillPatch& patch,
    AmrHydroGhostFillDiagnostics& diagnostics) {
  for (AmrHydroGhostDescriptor& ghost : patch.geometry->ghosts) {
    ghost.fill_status = AmrHydroGhostFillStatus::kRejectedStaleRemote;
    ++diagnostics.stale_epoch_rejections;
    ++diagnostics.unresolved_ghosts;
  }
}

// Temporal-history endpoint of one coarse source row. The thermal authority
// is the canonical specific internal energy, exactly as in
// loadAmrHydroConservedStateInto (same positive density, same pressure
// floor), so an alpha = 0 / 1 temporal ghost equals the state a synchronized
// fill reads at that time; the pressure cache is only the legacy fallback for
// rows that never set e. The auxiliary lane is rho * e of the same row (NaN
// when e is unknown), i.e. the value the solver's dual-energy lane holds.
struct HistoryEndpointState {
  hydro::HydroConservedState conserved{};
  double internal_energy_density_comoving = std::numeric_limits<double>::quiet_NaN();
};

[[nodiscard]] HistoryEndpointState historyEndpointForStateRow(
    const core::SimulationState& state,
    std::uint32_t row,
    double adiabatic_index) {
  if (row >= state.cells.size() || row >= state.gas_cells.size()) {
    throw std::runtime_error("AMR temporal history source row is outside SimulationState");
  }
  const double stored_density = state.gas_cells.density_code[row];
  const double rho = stored_density > 0.0 ? stored_density : std::max(state.cells.mass_code[row], 1.0e-14);
  const double canonical_e = state.gas_cells.internal_energy_code[row];
  const double pressure = canonical_e > 0.0
      ? std::max((adiabatic_index - 1.0) * rho * canonical_e, 1.0e-14)
      : std::max(state.gas_cells.pressure_code[row], 1.0e-14);
  const hydro::HydroPrimitiveState primitive{
      .rho_comoving = rho,
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = pressure,
      .metallicity_mass_fraction = std::clamp(
          state.gas_cells.metal_mass_code[row] /
              std::max(state.cells.mass_code[row], 1.0e-30),
          0.0,
          1.0)};
  HistoryEndpointState endpoint;
  endpoint.conserved = hydro::HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index);
  if (!finiteConserved(endpoint.conserved) || endpoint.conserved.mass_density_comoving <= 0.0) {
    throw std::runtime_error("AMR temporal history source state is non-finite or non-positive");
  }
  if (std::isfinite(canonical_e) && canonical_e > 0.0) {
    endpoint.internal_energy_density_comoving = endpoint.conserved.mass_density_comoving * canonical_e;
  }
  return endpoint;
}

void assignStart(core::AmrTemporalBoundaryHistoryCellRecord& record, const HistoryEndpointState& endpoint) {
  const hydro::HydroConservedState& state = endpoint.conserved;
  record.start_mass_density_comoving = state.mass_density_comoving;
  record.start_momentum_density_x_comoving = state.momentum_density_x_comoving;
  record.start_momentum_density_y_comoving = state.momentum_density_y_comoving;
  record.start_momentum_density_z_comoving = state.momentum_density_z_comoving;
  record.start_total_energy_density_comoving = state.total_energy_density_comoving;
  record.start_metal_mass_density_comoving = state.metal_mass_density_comoving;
  record.start_internal_energy_density_comoving = endpoint.internal_energy_density_comoving;
}

void assignEnd(core::AmrTemporalBoundaryHistoryCellRecord& record, const HistoryEndpointState& endpoint) {
  const hydro::HydroConservedState& state = endpoint.conserved;
  record.end_mass_density_comoving = state.mass_density_comoving;
  record.end_momentum_density_x_comoving = state.momentum_density_x_comoving;
  record.end_momentum_density_y_comoving = state.momentum_density_y_comoving;
  record.end_momentum_density_z_comoving = state.momentum_density_z_comoving;
  record.end_total_energy_density_comoving = state.total_energy_density_comoving;
  record.end_metal_mass_density_comoving = state.metal_mass_density_comoving;
  record.end_internal_energy_density_comoving = endpoint.internal_energy_density_comoving;
}

}  // namespace

std::uint64_t amrPatchGeometryFingerprint(const PatchDescriptor& patch) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&hash](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(patch.patch_id);
  mix(patch.parent_patch_id);
  mix(static_cast<std::uint64_t>(patch.level));
  mix(patch.morton_key);
  for (const double value : patch.origin_comov) {
    mix(std::bit_cast<std::uint64_t>(value));
  }
  for (const double value : patch.extent_comov) {
    mix(std::bit_cast<std::uint64_t>(value));
  }
  for (const std::uint16_t value : patch.cell_dims) {
    mix(static_cast<std::uint64_t>(value));
  }
  return hash == 0U ? 1U : hash;
}

namespace {

using HistoryRecord = core::AmrTemporalBoundaryHistoryRecord;
using HistoryCell = core::AmrTemporalBoundaryHistoryCellRecord;

std::uint64_t historyCandidateBound(std::size_t records, std::size_t cells) {
  auto bytes = static_cast<std::uint64_t>(core::checkedSizeAdd(
      core::checkedSizeMultiply(records, sizeof(HistoryRecord), "AMR history outer capacity"),
      core::checkedSizeMultiply(cells, sizeof(HistoryCell), "AMR history nested capacity"),
      "AMR history candidate capacity"));
  // One outer allocation and at most one nested allocation per record. This
  // allowance covers construction; retained allocator uncertainty remains in
  // the process policy's opaque allocator reserve, as for other state stores.
  return core::checkedMemoryBytesAdd(bytes,
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          core::checkedSizeAdd(records, 1U, "AMR history allocation count"), 64U,
          "AMR history allocation overhead")), "AMR history admission");
}

std::size_t historyCellCount(const core::SimulationState& state) {
  std::size_t count = 0U;
  for (const auto& record : state.amr_temporal_boundary_history.records()) {
    count = core::checkedSizeAdd(count, record.cells.size(), "AMR history cell count");
  }
  return count;
}

// Exact-capacity copying avoids implementation-dependent nested vector growth.
// The caller has admitted the complete candidate while canonical history lives.
void appendHistoryCopy(std::vector<HistoryRecord>& candidate, const HistoryRecord& source) {
  HistoryRecord record{
      .patch_id = source.patch_id,
      .patch_level = source.patch_level,
      .patch_geometry_fingerprint = source.patch_geometry_fingerprint,
      .gas_cell_identity_generation = source.gas_cell_identity_generation,
      .interval_start_code = source.interval_start_code,
      .interval_end_code = source.interval_end_code,
      .end_state_valid = source.end_state_valid,
      .cells = {}};
  record.cells.reserve(source.cells.size());
  if (record.cells.capacity() != source.cells.size()) {
    throw std::length_error("AMR history nested capacity exceeds admission");
  }
  record.cells.insert(record.cells.end(), source.cells.begin(), source.cells.end());
  candidate.push_back(std::move(record));
}

void publishHistory(core::SimulationState& state,
                    internal::AmrRetainedStateAdmission& admission,
                    std::vector<HistoryRecord>&& records) {
  core::AmrTemporalBoundaryHistoryStore candidate;
  candidate.assign(std::move(records));
  const auto old_owned = state.amr_temporal_boundary_history.ownedCapacityBytes();
  const auto new_owned = candidate.ownedCapacityBytes();
  static_assert(std::is_nothrow_move_assignable_v<core::AmrTemporalBoundaryHistoryStore>);
  admission.commit(old_owned, new_owned, [&]() noexcept {
    state.amr_temporal_boundary_history = std::move(candidate);
  });
}

void captureHistoryStart(core::SimulationState& state,
                         std::span<const PatchDescriptor> patches, int level,
                         double interval_start_code, double adiabatic_index,
                         core::MemoryGovernor* governor) {
  if (!std::isfinite(interval_start_code)) {
    throw std::invalid_argument("captureAmrTemporalBoundaryHistoryStart requires finite interval_start_code");
  }
  for (const auto& record : state.amr_temporal_boundary_history.records()) {
    if (level < 0 || static_cast<int>(record.patch_level) == level) {
      throw std::runtime_error("cannot begin an AMR temporal interval while that level's interval is active");
    }
  }
  const auto isCoarseSource = [&](const PatchDescriptor& patch) {
    if (level >= 0 && static_cast<int>(patch.level) != level) return false;
    return std::any_of(patches.begin(), patches.end(), [&](const PatchDescriptor& other) {
      return level < 0 ? other.level > patch.level : static_cast<int>(other.level) == level + 1;
    });
  };
  std::size_t coarse_count = 0U;
  std::size_t cells = historyCellCount(state);
  std::uint64_t geometry_bound = 0U;
  std::uint64_t geometry_scratch_bound = 0U;
  for (const auto& patch : patches) {
    if (!isCoarseSource(patch)) continue;
    if (state.amr_temporal_boundary_history.findByPatchId(patch.patch_id) != nullptr) {
      throw std::runtime_error("AMR history candidate duplicates a retained patch identity");
    }
    for (const auto& other : patches) {
      if (&other != &patch && other.patch_id == patch.patch_id) {
        throw std::runtime_error("AMR history candidate has duplicate patch descriptors");
      }
    }
    ++coarse_count;
    const auto capacity = amrHydroGeometryCapacity(patch);
    cells = core::checkedSizeAdd(cells, capacity.real_cells, "AMR history cells");
    geometry_bound = std::max(geometry_bound, capacity.retained_bytes);
    geometry_scratch_bound = std::max(geometry_scratch_bound, capacity.construction_scratch_bytes);
  }
  if (coarse_count == 0U) return;
  const auto count = core::checkedSizeAdd(coarse_count, state.amr_temporal_boundary_history.size(),
                                         "AMR history merged count");
  internal::AmrRetainedStateAdmission admission(
      state, governor, historyCandidateBound(count, cells), "amr.temporal_history.replacement");
  core::MemoryReservation geometry_reservation;
  if (governor != nullptr) {
    geometry_reservation = governor->reserve(core::MemoryClass::kPhaseResident,
        core::checkedMemoryBytesAdd(geometry_bound, geometry_scratch_bound,
                                    "AMR temporal history geometry and scratch"),
        "amr.temporal_history.geometry");
    geometry_reservation.commit();
  }
  std::vector<HistoryRecord> records;
  records.reserve(count);
  if (records.capacity() != count) throw std::length_error("AMR history outer capacity exceeds admission");
  for (const auto& patch : patches) {
    if (!isCoarseSource(patch)) continue;
    const AmrHydroPatchGeometry geometry = buildAmrHydroPatchGeometry(state, patch);
    HistoryRecord record;
    record.patch_id = patch.patch_id;
    record.patch_level = patch.level;
    record.patch_geometry_fingerprint = amrPatchGeometryFingerprint(patch);
    record.gas_cell_identity_generation = state.gasCellIdentityGeneration();
    record.interval_start_code = interval_start_code;
    record.interval_end_code = interval_start_code;
    record.cells.reserve(geometry.real_cells.size());
    if (record.cells.capacity() != geometry.real_cells.size()) {
      throw std::length_error("AMR history cell capacity exceeds admission");
    }
    for (const auto& descriptor : geometry.real_cells) {
      HistoryCell cell;
      cell.gas_cell_id = descriptor.gas_cell_id;
      cell.patch_local_cell = descriptor.patch_local_cell;
      assignStart(cell, historyEndpointForStateRow(state, descriptor.local_cell_row, adiabatic_index));
      record.cells.push_back(cell);
    }
    records.push_back(std::move(record));
  }
  for (const auto& record : state.amr_temporal_boundary_history.records()) appendHistoryCopy(records, record);
  publishHistory(state, admission, std::move(records));
}

}  // namespace

void captureAmrTemporalBoundaryHistoryStart(
    core::SimulationState& state, std::span<const PatchDescriptor> patches,
    double interval_start_code, double adiabatic_index, core::MemoryGovernor* governor) {
  captureHistoryStart(state, patches, -1, interval_start_code, adiabatic_index, governor);
}

void captureAmrTemporalBoundaryHistoryEnd(
    core::SimulationState& state, std::span<const PatchDescriptor> patches,
    double interval_end_code, double adiabatic_index, core::MemoryGovernor* governor) {
  captureAmrTemporalBoundaryHistoryEndForLevel(state, patches, -1, interval_end_code, adiabatic_index, governor);
}

void captureAmrTemporalBoundaryHistoryStartForLevel(
    core::SimulationState& state, std::span<const PatchDescriptor> patches, int level,
    double interval_start_code, double adiabatic_index, core::MemoryGovernor* governor) {
  if (level < 0) throw std::invalid_argument("AMR history start requires level >= 0");
  captureHistoryStart(state, patches, level, interval_start_code, adiabatic_index, governor);
}

void retireAmrTemporalBoundaryHistoryLevel(core::SimulationState& state, int level,
                                          core::MemoryGovernor* governor) {
  // Erase only moves/destroys vector owners; it cannot allocate or fail halfway.
  const auto old_owned = state.amr_temporal_boundary_history.ownedCapacityBytes();
  auto new_owned = old_owned;
  for (const auto& record : state.amr_temporal_boundary_history.records()) {
    if (level < 0 || static_cast<int>(record.patch_level) == level) {
      new_owned -= static_cast<std::uint64_t>(core::checkedSizeMultiply(
          record.cells.capacity(), sizeof(HistoryCell), "AMR retired history cells"));
    }
  }
  internal::AmrRetainedStateAdmission admission(state, governor, 0U, "amr.temporal_history.retirement");
  admission.commitShrink(old_owned, new_owned, [&]() noexcept {
    state.amr_temporal_boundary_history.eraseLevel(level);
  });
}

void captureAmrTemporalBoundaryHistoryEndForLevel(
    core::SimulationState& state, std::span<const PatchDescriptor> patches, int level,
    double interval_end_code, double adiabatic_index, core::MemoryGovernor* governor) {
  if (!std::isfinite(interval_end_code)) {
    throw std::invalid_argument("captureAmrTemporalBoundaryHistoryEnd requires finite interval_end_code");
  }
  const auto existing = state.amr_temporal_boundary_history.records();
  if (std::none_of(existing.begin(), existing.end(), [&](const auto& record) {
        return level < 0 || static_cast<int>(record.patch_level) == level;
      })) return;
  internal::AmrRetainedStateAdmission admission(state, governor,
      historyCandidateBound(existing.size(), historyCellCount(state)), "amr.temporal_history.end_candidate");
  std::vector<HistoryRecord> records;
  records.reserve(existing.size());
  if (records.capacity() != existing.size()) throw std::length_error("AMR history outer capacity exceeds admission");
  for (const auto& source : existing) appendHistoryCopy(records, source);
  for (auto& record : records) {
    if (level >= 0 && static_cast<int>(record.patch_level) != level) continue;
    const auto patch_it = std::find_if(patches.begin(), patches.end(),
        [&](const auto& patch) { return patch.patch_id == record.patch_id; });
    if (patch_it == patches.end() || record.patch_geometry_fingerprint != amrPatchGeometryFingerprint(*patch_it)) {
      throw std::runtime_error("cannot complete AMR temporal history after patch geometry changed");
    }
    if (record.gas_cell_identity_generation != state.gasCellIdentityGeneration()) {
      throw std::runtime_error("cannot complete AMR temporal history after gas-cell identity generation changed");
    }
    if (interval_end_code <= record.interval_start_code) {
      throw std::invalid_argument("AMR temporal history interval must have positive duration");
    }
    for (auto& cell : record.cells) {
      const auto row = state.rowForGasCellId(cell.gas_cell_id);
      if (!row.has_value()) throw std::runtime_error("cannot complete AMR history: missing coarse gas_cell_id");
      if (state.cells.patch_index[*row] >= state.patches.size() ||
          state.patches.patch_id[state.cells.patch_index[*row]] != record.patch_id) {
        throw std::runtime_error("cannot complete AMR temporal history because coarse gas cell changed patch ownership");
      }
      assignEnd(cell, historyEndpointForStateRow(state, *row, adiabatic_index));
    }
    record.interval_end_code = interval_end_code;
    record.end_state_valid = true;
  }
  publishHistory(state, admission, std::move(records));
}

void retireAmrTemporalBoundaryHistory(core::SimulationState& state, core::MemoryGovernor* governor) {
  retireAmrTemporalBoundaryHistoryLevel(state, -1, governor);
}

std::uint64_t amrGhostFillScratchBytes(
    std::span<const AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources) {
  std::size_t source_cells = 0U;
  for (const auto& patch : patches) {
    if (patch.geometry != nullptr) {
      source_cells = core::checkedSizeAdd(source_cells,
          patch.geometry->geometry.cellCount(), "AMR ghost local source extent");
    }
  }
  for (const auto& source : remote_sources) {
    source_cells = core::checkedSizeAdd(source_cells, source.cells.size(),
                                      "AMR ghost remote source extent");
  }
  return core::checkedMemoryBytesAdd(
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          source_cells, sizeof(SelectedSourceCell), "AMR ghost selection capacity")),
      256U, "AMR ghost selection alignment");
}

AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index, core::MemoryGovernor* governor) {
  const std::uint64_t scratch_bytes = amrGhostFillScratchBytes(patches, remote_sources);
  core::GovernedScratchArena arena(governor, core::MemoryClass::kScratchArena,
                                   scratch_bytes, "amr.ghost_fill.selection");
  CellSelection selection(arena.resource());
  const std::size_t source_capacity = static_cast<std::size_t>(
      (scratch_bytes - 256U) / sizeof(SelectedSourceCell));
  selection.cells.reserve(source_capacity);

  AmrHydroGhostFillDiagnostics diagnostics;
  for (const AmrHydroGhostFillPatch& patch : patches) {
    validatePatchView(patch);
  }
  for (const AmrHydroSparseGhostSource& source : remote_sources) {
    if (source.patch == nullptr) {
      throw std::invalid_argument(
          "fillAmrHydroGhostCells: sparse remote source is missing its patch descriptor");
    }
    if (!std::isfinite(source.source_current_state_time_code)) {
      throw std::invalid_argument(
          "fillAmrHydroGhostCells: sparse remote source time must be finite");
    }
    const std::size_t patch_cell_count =
        static_cast<std::size_t>(source.patch->cell_dims[0]) *
        static_cast<std::size_t>(source.patch->cell_dims[1]) *
        static_cast<std::size_t>(source.patch->cell_dims[2]);
    std::size_t previous_offset = 0U;
    bool have_previous = false;
    for (const AmrHydroSparseRemoteCell& cell : source.cells) {
      const std::size_t offset = cell.patch_local_cell;
      if (offset >= patch_cell_count || cell.gas_cell_id == 0U) {
        throw std::invalid_argument(
            "fillAmrHydroGhostCells: sparse remote cell is outside its patch or has zero identity");
      }
      if (have_previous && offset <= previous_offset) {
        throw std::invalid_argument(
            "fillAmrHydroGhostCells: sparse remote cells must be strictly sorted by patch-local offset");
      }
      previous_offset = offset;
      have_previous = true;
    }
  }

  for (AmrHydroGhostFillPatch& patch : patches) {
    if (!patch.requires_ghost_fill) {
      continue;
    }
    if (staleRemotePatch(patch)) {
      markStaleRemoteTarget(patch, diagnostics);
      continue;
    }
    validateTargetPhysicalFacesOnDomainWalls(patch, diagnostics);
    hydro::fillHydroBoundaryGhostCells(*patch.conserved, patch.geometry->geometry, adiabatic_index);
    // H02: physical ghosts carry the auxiliary lane of the cell they copy or
    // mirror (no-op when the target SoA has no auxiliary lane).
    hydro::fillHydroBoundaryGhostAuxiliaryInternalEnergy(*patch.conserved, patch.geometry->geometry);
    for (AmrHydroGhostDescriptor& ghost : patch.geometry->ghosts) {
      if (ghost.edge_stencil) {
        fillEdgeStencilGhost(
            patch, ghost, patches, remote_sources, adiabatic_index, diagnostics, selection);
        continue;
      }
      if (ghost.boundary_class == AmrHydroBoundaryClass::kPhysical) {
        if (ghost.physical_closure_mapped) {
          // P2-55-03: deep reflective layer of a patch thinner than the
          // mirror depth; the mirror lies beyond the opposite face.
          const hydro::HydroGhostCell& cell = patch.geometry->geometry.ghost_cells.at(ghost.ghost_slot);
          (void)fillClosureMappedGhost(
              patch, ghost, patches, remote_sources, adiabatic_index, diagnostics, selection,
              deepLayerStencilIndex(patch.geometry->geometry, cell, ghost.ghost_layer),
              /*require_wall=*/true);
          continue;
        }
        ghost.fill_status = AmrHydroGhostFillStatus::kFilledPhysicalBoundary;
        ++diagnostics.physical_ghosts_filled;
        continue;
      }
      fillAmrGhost(
          patch, ghost, patches, remote_sources, adiabatic_index, diagnostics, selection);
    }
  }
  return diagnostics;
}

AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index) {
  return fillAmrHydroGhostCells(patches, remote_sources, adiabatic_index, nullptr);
}

AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    double adiabatic_index) {
  return fillAmrHydroGhostCells(patches, {}, adiabatic_index);
}

}  // namespace cosmosim::amr
