#include "cosmosim/hydro/hydro_cartesian_patch.hpp"

#include "cosmosim/core/checked_arithmetic.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmosim::hydro {
namespace {

void validateSpec(const HydroCartesianPatchSpec& spec) {
  if (spec.nx == 0 || spec.ny == 0 || spec.nz == 0) {
    throw std::invalid_argument("Cartesian hydro patch dimensions must be positive");
  }
  if (!std::isfinite(spec.origin_x_comoving) || !std::isfinite(spec.origin_y_comoving) ||
      !std::isfinite(spec.origin_z_comoving)) {
    throw std::invalid_argument("Cartesian hydro patch origins must be finite");
  }
  if (!std::isfinite(spec.cell_width_x_comoving) || spec.cell_width_x_comoving <= 0.0 ||
      !std::isfinite(spec.cell_width_y_comoving) || spec.cell_width_y_comoving <= 0.0 ||
      !std::isfinite(spec.cell_width_z_comoving) || spec.cell_width_z_comoving <= 0.0) {
    throw std::invalid_argument("Cartesian hydro patch cell widths must be finite and positive");
  }
}


}  // namespace

std::size_t HydroPatchGeometry::cellCount() const {
  return core::checkedSizeProduct3(nx, ny, nz, "HydroPatchGeometry::cellCount");
}

std::size_t hydroEdgeGhostCount(std::size_t nx, std::size_t ny, std::size_t nz, std::size_t depth) {
  if (depth == 0U) {
    return 0U;
  }
  const std::array<std::size_t, 3> extent{nx, ny, nz};
  std::size_t count = 0U;
  for (std::size_t n = 0U; n < 3U; ++n) {
    for (std::size_t t = 0U; t < 3U; ++t) {
      if (t == n || extent[t] <= 1U) continue;
      const std::size_t u = 3U - n - t;
      // 2 normal sides x 2 transverse sides x depth x extent[u].
      count = core::checkedSizeAdd(
          count,
          core::checkedSizeMultiply(core::checkedSizeMultiply(4U, depth, "edge ghost depth"),
                                    extent[u], "edge ghost block"),
          "edge ghost count");
    }
  }
  return count;
}

std::size_t HydroPatchGeometry::edgeGhostSlot(
    std::size_t normal_axis, std::size_t normal_side,
    std::size_t transverse_axis, std::size_t transverse_side,
    std::size_t depth, std::size_t third_index) const {
  if (edge_ghost_depth == 0U || edge_ghost_slot_begin == k_invalid_ghost_cell_slot ||
      normal_axis > 2U || transverse_axis > 2U || transverse_axis == normal_axis ||
      normal_side > 1U || transverse_side > 1U || depth == 0U || depth > edge_ghost_depth) {
    return k_invalid_ghost_cell_slot;
  }
  const std::array<std::size_t, 3> extent{nx, ny, nz};
  const std::size_t d = edge_ghost_depth;
  std::size_t offset = 0U;
  for (std::size_t n = 0U; n < 3U; ++n) {
    for (std::size_t side = 0U; side < 2U; ++side) {
      for (std::size_t t = 0U; t < 3U; ++t) {
        if (t == n || extent[t] <= 1U) continue;
        const std::size_t u = 3U - n - t;
        if (n == normal_axis && side == normal_side && t == transverse_axis) {
          if (third_index >= extent[u]) return k_invalid_ghost_cell_slot;
          const std::size_t slot = edge_ghost_slot_begin + offset +
              (transverse_side * d + (depth - 1U)) * extent[u] + third_index;
          return slot < ghost_cells.size() ? slot : k_invalid_ghost_cell_slot;
        }
        offset += 2U * d * extent[u];
      }
    }
  }
  return k_invalid_ghost_cell_slot;
}

std::size_t HydroPatchGeometry::totalCellStorageCount() const {
  return core::checkedSizeAdd(cellCount(), ghost_cells.size(), "HydroPatchGeometry::totalCellStorageCount");
}

std::size_t HydroPatchGeometry::linearCellIndex(std::size_t i, std::size_t j, std::size_t k) const {
  if (i >= nx || j >= ny || k >= nz) {
    throw std::out_of_range("HydroPatchGeometry::linearCellIndex index out of range");
  }
  return (k * ny + j) * nx + i;
}

std::array<std::size_t, 3> HydroPatchGeometry::cellIjk(std::size_t row) const {
  if (row >= cellCount() || nx == 0 || ny == 0) {
    throw std::out_of_range("HydroPatchGeometry::cellIjk row out of range");
  }
  const std::size_t i = row % nx;
  const std::size_t j = (row / nx) % ny;
  const std::size_t k = row / (nx * ny);
  return {i, j, k};
}

std::size_t HydroPatchGeometry::neighborCell(
    std::size_t row,
    int di,
    int dj,
    int dk) const {
  if (nx == 0 || ny == 0 || nz == 0 || row >= cellCount()) {
    return k_invalid_cell_index;
  }
  const std::size_t i = row % nx;
  const std::size_t j = (row / nx) % ny;
  const std::size_t k = row / core::checkedSizeMultiply(nx, ny, "HydroPatchGeometry::neighborCell");

  const auto apply_delta = [](std::size_t index, int delta, std::size_t extent) -> std::size_t {
    if (delta < 0) {
      const std::size_t magnitude = static_cast<std::size_t>(-static_cast<long long>(delta));
      return magnitude > index ? k_invalid_cell_index : index - magnitude;
    }
    const std::size_t magnitude = static_cast<std::size_t>(delta);
    if (magnitude >= extent || index > extent - 1U - magnitude) {
      return k_invalid_cell_index;
    }
    return index + magnitude;
  };

  const std::size_t ni = apply_delta(i, di, nx);
  const std::size_t nj = apply_delta(j, dj, ny);
  const std::size_t nk = apply_delta(k, dk, nz);
  if (ni == k_invalid_cell_index || nj == k_invalid_cell_index || nk == k_invalid_cell_index) {
    return k_invalid_cell_index;
  }
  return (nk * ny + nj) * nx + ni;
}

HydroPatchGeometry makeCartesianPatchGeometry(const HydroCartesianPatchSpec& spec) {
  validateSpec(spec);
  const std::size_t cell_count = core::checkedSizeProduct3(
      spec.nx, spec.ny, spec.nz, "makeCartesianPatchGeometry");

  HydroPatchGeometry geometry;
  geometry.nx = spec.nx;
  geometry.ny = spec.ny;
  geometry.nz = spec.nz;
  geometry.origin_x_comoving = spec.origin_x_comoving;
  geometry.origin_y_comoving = spec.origin_y_comoving;
  geometry.origin_z_comoving = spec.origin_z_comoving;
  geometry.cell_width_x_comoving = spec.cell_width_x_comoving;
  geometry.cell_width_y_comoving = spec.cell_width_y_comoving;
  geometry.cell_width_z_comoving = spec.cell_width_z_comoving;
  geometry.cell_volume_comoving =
      spec.cell_width_x_comoving * spec.cell_width_y_comoving * spec.cell_width_z_comoving;
  const double x_face_area = spec.cell_width_y_comoving * spec.cell_width_z_comoving;
  const double y_face_area = spec.cell_width_x_comoving * spec.cell_width_z_comoving;
  const double z_face_area = spec.cell_width_x_comoving * spec.cell_width_y_comoving;
  if (!std::isfinite(geometry.cell_volume_comoving) || geometry.cell_volume_comoving <= 0.0 ||
      !std::isfinite(x_face_area) || !std::isfinite(y_face_area) || !std::isfinite(z_face_area)) {
    throw std::overflow_error("Cartesian hydro patch derived volume/face areas must be finite and positive");
  }

  const std::size_t x_faces = (spec.nx > 1U)
      ? core::checkedSizeProduct3(spec.nx - 1U, spec.ny, spec.nz, "Cartesian hydro x-face count")
      : 0U;
  const std::size_t y_faces = (spec.ny > 1U)
      ? core::checkedSizeProduct3(spec.nx, spec.ny - 1U, spec.nz, "Cartesian hydro y-face count")
      : 0U;
  const std::size_t z_faces = (spec.nz > 1U)
      ? core::checkedSizeProduct3(spec.nx, spec.ny, spec.nz - 1U, "Cartesian hydro z-face count")
      : 0U;
  const std::size_t face_count = core::checkedSizeAdd(
      core::checkedSizeAdd(x_faces, y_faces, "Cartesian hydro face count"),
      z_faces,
      "Cartesian hydro face count");
  geometry.faces.reserve(face_count);

  for (std::size_t k = 0; k < spec.nz; ++k) {
    for (std::size_t j = 0; j < spec.ny; ++j) {
      for (std::size_t i = 0; i + 1U < spec.nx; ++i) {
        const std::size_t owner = (k * spec.ny + j) * spec.nx + i;
        const std::size_t neighbor = owner + 1U;
        geometry.faces.push_back(HydroFace{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = geometry.neighborCell(owner, -1, 0, 0),
            .neighbor_plus_cell = geometry.neighborCell(neighbor, 1, 0, 0),
            .area_comoving = x_face_area,
            .normal_x = 1.0,
            .normal_y = 0.0,
            .normal_z = 0.0,
            .axis = HydroFaceAxis::kX,
            .owner_minus2_cell = geometry.neighborCell(owner, -2, 0, 0),
            .neighbor_plus2_cell = geometry.neighborCell(neighbor, 2, 0, 0)});
      }
    }
  }

  for (std::size_t k = 0; k < spec.nz; ++k) {
    for (std::size_t j = 0; j + 1U < spec.ny; ++j) {
      for (std::size_t i = 0; i < spec.nx; ++i) {
        const std::size_t owner = (k * spec.ny + j) * spec.nx + i;
        const std::size_t neighbor = owner + spec.nx;
        geometry.faces.push_back(HydroFace{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = geometry.neighborCell(owner, 0, -1, 0),
            .neighbor_plus_cell = geometry.neighborCell(neighbor, 0, 1, 0),
            .area_comoving = y_face_area,
            .normal_x = 0.0,
            .normal_y = 1.0,
            .normal_z = 0.0,
            .axis = HydroFaceAxis::kY,
            .owner_minus2_cell = geometry.neighborCell(owner, 0, -2, 0),
            .neighbor_plus2_cell = geometry.neighborCell(neighbor, 0, 2, 0)});
      }
    }
  }

  for (std::size_t k = 0; k + 1U < spec.nz; ++k) {
    for (std::size_t j = 0; j < spec.ny; ++j) {
      for (std::size_t i = 0; i < spec.nx; ++i) {
        const std::size_t owner = (k * spec.ny + j) * spec.nx + i;
        const std::size_t neighbor = owner + spec.nx * spec.ny;
        geometry.faces.push_back(HydroFace{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = geometry.neighborCell(owner, 0, 0, -1),
            .neighbor_plus_cell = geometry.neighborCell(neighbor, 0, 0, 1),
            .area_comoving = z_face_area,
            .normal_x = 0.0,
            .normal_y = 0.0,
            .normal_z = 1.0,
            .axis = HydroFaceAxis::kZ,
            .owner_minus2_cell = geometry.neighborCell(owner, 0, 0, -2),
            .neighbor_plus2_cell = geometry.neighborCell(neighbor, 0, 0, 2)});
      }
    }
  }

  if (geometry.cellCount() != cell_count) {
    throw std::logic_error("Cartesian hydro patch cell count mismatch");
  }
  return geometry;
}

std::array<std::size_t, 3> chooseNearCubicCartesianFactors(std::size_t cell_count) {
  if (cell_count == 0) {
    return {0, 0, 0};
  }

  std::array<std::size_t, 3> best{cell_count, 1U, 1U};
  double best_score = std::numeric_limits<double>::infinity();
  for (std::size_t nz = 1; nz * nz * nz <= cell_count; ++nz) {
    if (cell_count % nz != 0) {
      continue;
    }
    const std::size_t rem = cell_count / nz;
    for (std::size_t ny = nz; ny * ny <= rem; ++ny) {
      if (rem % ny != 0) {
        continue;
      }
      const std::size_t nx = rem / ny;
      if (nx < ny) {
        continue;
      }
      const double score =
          static_cast<double>(nx - ny) * static_cast<double>(nx - ny) +
          static_cast<double>(ny - nz) * static_cast<double>(ny - nz) +
          static_cast<double>(nx - nz) * static_cast<double>(nx - nz);
      if (score < best_score) {
        best = {nx, ny, nz};
        best_score = score;
      }
    }
  }
  return best;
}

}  // namespace cosmosim::hydro
