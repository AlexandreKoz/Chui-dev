#include "cosmosim/hydro/hydro_boundary_conditions.hpp"

#include <array>
#include <stdexcept>

namespace cosmosim::hydro {
namespace {

[[nodiscard]] double faceArea(const HydroPatchGeometry& geometry, HydroFaceAxis axis) {
  switch (axis) {
    case HydroFaceAxis::kX:
      return geometry.cell_width_y_comoving * geometry.cell_width_z_comoving;
    case HydroFaceAxis::kY:
      return geometry.cell_width_x_comoving * geometry.cell_width_z_comoving;
    case HydroFaceAxis::kZ:
      return geometry.cell_width_x_comoving * geometry.cell_width_y_comoving;
  }
  return 0.0;
}

void setNormal(HydroFace& face, HydroFaceAxis axis, HydroFaceSide side) {
  const double sign = side == HydroFaceSide::kLower ? -1.0 : 1.0;
  face.normal_x = axis == HydroFaceAxis::kX ? sign : 0.0;
  face.normal_y = axis == HydroFaceAxis::kY ? sign : 0.0;
  face.normal_z = axis == HydroFaceAxis::kZ ? sign : 0.0;
}

[[nodiscard]] std::size_t boundaryOwner(
    const HydroPatchGeometry& geometry,
    HydroFaceAxis axis,
    HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  switch (axis) {
    case HydroFaceAxis::kX: {
      const std::size_t i = side == HydroFaceSide::kLower ? 0U : geometry.nx - 1U;
      return geometry.linearCellIndex(i, a, b);
    }
    case HydroFaceAxis::kY: {
      const std::size_t j = side == HydroFaceSide::kLower ? 0U : geometry.ny - 1U;
      return geometry.linearCellIndex(a, j, b);
    }
    case HydroFaceAxis::kZ: {
      const std::size_t k = side == HydroFaceSide::kLower ? 0U : geometry.nz - 1U;
      return geometry.linearCellIndex(a, b, k);
    }
  }
  return k_invalid_cell_index;
}

[[nodiscard]] std::size_t boundarySource(
    const HydroPatchGeometry& geometry,
    HydroBoundaryKind boundary_kind,
    HydroFaceAxis axis,
    HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  if (boundary_kind != HydroBoundaryKind::kPeriodic) {
    return boundaryOwner(geometry, axis, side, a, b);
  }

  switch (axis) {
    case HydroFaceAxis::kX: {
      const std::size_t i = side == HydroFaceSide::kLower ? geometry.nx - 1U : 0U;
      return geometry.linearCellIndex(i, a, b);
    }
    case HydroFaceAxis::kY: {
      const std::size_t j = side == HydroFaceSide::kLower ? geometry.ny - 1U : 0U;
      return geometry.linearCellIndex(a, j, b);
    }
    case HydroFaceAxis::kZ: {
      const std::size_t k = side == HydroFaceSide::kLower ? geometry.nz - 1U : 0U;
      return geometry.linearCellIndex(a, b, k);
    }
  }
  return k_invalid_cell_index;
}

void appendBoundaryGhost(
    HydroPatchGeometry& geometry,
    HydroBoundaryKind boundary_kind,
    HydroFaceAxis axis,
    HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  const std::size_t ghost_slot = geometry.ghost_cells.size();
  const std::size_t ghost_cell = geometry.cellCount() + ghost_slot;
  const std::size_t owner = boundaryOwner(geometry, axis, side, a, b);
  const std::size_t source = boundarySource(geometry, boundary_kind, axis, side, a, b);

  geometry.ghost_cells.push_back(HydroGhostCell{
      .owner_real_cell = owner,
      .source_real_cell = source,
      .ghost_cell = ghost_cell,
      .ghost_slot = ghost_slot,
      .boundary_kind = boundary_kind,
      .axis = axis,
      .side = side,
      .mutation_rights = HydroGhostMutationRights::kWritablePhysicalBoundaryScratch});

  HydroFace face{
      .owner_cell = owner,
      .neighbor_cell = ghost_cell,
      .owner_minus_cell = k_invalid_cell_index,
      .neighbor_plus_cell = k_invalid_cell_index,
      .ghost_cell_slot = ghost_slot,
      .area_comoving = faceArea(geometry, axis),
      .axis = axis};
  setNormal(face, axis, side);
  geometry.faces.push_back(face);
}

void appendAxisBoundaryGhosts(HydroPatchGeometry& geometry, HydroBoundaryKind boundary_kind, HydroFaceAxis axis) {
  switch (axis) {
    case HydroFaceAxis::kX:
      for (std::size_t k = 0; k < geometry.nz; ++k) {
        for (std::size_t j = 0; j < geometry.ny; ++j) {
          appendBoundaryGhost(geometry, boundary_kind, axis, HydroFaceSide::kLower, j, k);
          appendBoundaryGhost(geometry, boundary_kind, axis, HydroFaceSide::kUpper, j, k);
        }
      }
      break;
    case HydroFaceAxis::kY:
      for (std::size_t k = 0; k < geometry.nz; ++k) {
        for (std::size_t i = 0; i < geometry.nx; ++i) {
          appendBoundaryGhost(geometry, boundary_kind, axis, HydroFaceSide::kLower, i, k);
          appendBoundaryGhost(geometry, boundary_kind, axis, HydroFaceSide::kUpper, i, k);
        }
      }
      break;
    case HydroFaceAxis::kZ:
      for (std::size_t j = 0; j < geometry.ny; ++j) {
        for (std::size_t i = 0; i < geometry.nx; ++i) {
          appendBoundaryGhost(geometry, boundary_kind, axis, HydroFaceSide::kLower, i, j);
          appendBoundaryGhost(geometry, boundary_kind, axis, HydroFaceSide::kUpper, i, j);
        }
      }
      break;
  }
}

void reflectNormalVelocity(HydroPrimitiveState& primitive, HydroFaceAxis axis) {
  switch (axis) {
    case HydroFaceAxis::kX:
      primitive.vel_x_peculiar = -primitive.vel_x_peculiar;
      break;
    case HydroFaceAxis::kY:
      primitive.vel_y_peculiar = -primitive.vel_y_peculiar;
      break;
    case HydroFaceAxis::kZ:
      primitive.vel_z_peculiar = -primitive.vel_z_peculiar;
      break;
  }
}

void validateGhostMetadata(const HydroPatchGeometry& geometry, const HydroGhostCell& ghost) {
  const std::size_t real_cell_count = geometry.cellCount();
  if (ghost.owner_real_cell >= real_cell_count || ghost.source_real_cell >= real_cell_count) {
    throw std::invalid_argument("Hydro ghost metadata must reference real owner/source cells");
  }
  if (ghost.ghost_slot >= geometry.ghost_cells.size()) {
    throw std::invalid_argument("Hydro ghost metadata has an invalid ghost slot");
  }
  if (ghost.ghost_cell != real_cell_count + ghost.ghost_slot) {
    throw std::invalid_argument("Hydro ghost metadata row does not match real-cell offset");
  }
}

[[nodiscard]] std::size_t periodicWrap(
    const HydroPatchGeometry& geometry, std::size_t row, HydroFaceAxis axis, long long delta) {
  const std::array<std::size_t, 3> ijk = geometry.cellIjk(row);
  std::array<std::size_t, 3> wrapped = ijk;
  const std::size_t axis_index = axis == HydroFaceAxis::kX ? 0U : (axis == HydroFaceAxis::kY ? 1U : 2U);
  const std::array<std::size_t, 3> extent{geometry.nx, geometry.ny, geometry.nz};
  const long long n = static_cast<long long>(extent[axis_index]);
  long long index = (static_cast<long long>(ijk[axis_index]) + delta) % n;
  if (index < 0) {
    index += n;
  }
  wrapped[axis_index] = static_cast<std::size_t>(index);
  return geometry.linearCellIndex(wrapped[0], wrapped[1], wrapped[2]);
}

[[nodiscard]] double normalSignAlongAxis(const HydroFace& face) {
  switch (face.axis) {
    case HydroFaceAxis::kX:
      return face.normal_x >= 0.0 ? 1.0 : -1.0;
    case HydroFaceAxis::kY:
      return face.normal_y >= 0.0 ? 1.0 : -1.0;
    case HydroFaceAxis::kZ:
      return face.normal_z >= 0.0 ? 1.0 : -1.0;
  }
  return 1.0;
}

}  // namespace

std::size_t completePeriodicReconstructionStencils(HydroPatchGeometry& geometry) {
  const std::size_t real_cell_count = geometry.cellCount();
  for (const HydroGhostCell& ghost : geometry.ghost_cells) {
    if (ghost.boundary_kind != HydroBoundaryKind::kPeriodic) {
      throw std::invalid_argument(
          "periodic reconstruction stencils require every ghost cell to be periodic");
    }
  }
  std::size_t filled = 0U;
  for (HydroFace& face : geometry.faces) {
    if (face.owner_cell >= real_cell_count) {
      throw std::invalid_argument("periodic stencil completion requires real face owners");
    }
    std::size_t neighbor_real = face.neighbor_cell;
    if (face.ghost_cell_slot != k_invalid_ghost_cell_slot) {
      neighbor_real = geometry.ghost_cells.at(face.ghost_cell_slot).source_real_cell;
    }
    if (neighbor_real == k_invalid_cell_index || neighbor_real >= real_cell_count) {
      continue;
    }
    const long long step = normalSignAlongAxis(face) > 0.0 ? 1 : -1;
    const auto fill = [&](std::size_t& slot, std::size_t value) {
      if (slot == k_invalid_cell_index) {
        slot = value;
        ++filled;
      }
    };
    fill(face.owner_minus_cell, periodicWrap(geometry, face.owner_cell, face.axis, -step));
    fill(face.owner_minus2_cell, periodicWrap(geometry, face.owner_cell, face.axis, -2 * step));
    fill(face.neighbor_plus_cell, periodicWrap(geometry, neighbor_real, face.axis, step));
    fill(face.neighbor_plus2_cell, periodicWrap(geometry, neighbor_real, face.axis, 2 * step));
  }
  return filled;
}

void appendCartesianBoundaryGhostFaces(
    HydroPatchGeometry& geometry,
    HydroBoundaryKind boundary_kind) {
  if (boundary_kind == HydroBoundaryKind::kInterior || boundary_kind == HydroBoundaryKind::kImportedMpi) {
    throw std::invalid_argument("Physical Cartesian boundary ghosts require periodic, open, or reflective boundary kind");
  }
  if (geometry.nx == 0 || geometry.ny == 0 || geometry.nz == 0 || geometry.cellCount() == 0) {
    throw std::invalid_argument("Physical Cartesian boundary ghosts require a non-empty patch");
  }
  if (!geometry.ghost_cells.empty()) {
    throw std::invalid_argument("Cartesian boundary ghosts have already been appended to this geometry");
  }

  appendAxisBoundaryGhosts(geometry, boundary_kind, HydroFaceAxis::kX);
  appendAxisBoundaryGhosts(geometry, boundary_kind, HydroFaceAxis::kY);
  appendAxisBoundaryGhosts(geometry, boundary_kind, HydroFaceAxis::kZ);
}

void fillHydroBoundaryGhostCells(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    double adiabatic_index) {
  if (conserved.size() < geometry.totalCellStorageCount()) {
    throw std::invalid_argument("Hydro conserved state does not include storage for boundary ghosts");
  }

  for (const HydroGhostCell& ghost : geometry.ghost_cells) {
    validateGhostMetadata(geometry, ghost);
    if (ghost.mutation_rights == HydroGhostMutationRights::kReadOnlyImported ||
        ghost.boundary_kind == HydroBoundaryKind::kImportedMpi) {
      continue;
    }
    HydroPrimitiveState primitive =
        HydroCoreSolver::primitiveFromConserved(conserved.loadCell(ghost.source_real_cell), adiabatic_index);
    if (ghost.boundary_kind == HydroBoundaryKind::kReflective) {
      reflectNormalVelocity(primitive, ghost.axis);
    }
    conserved.storeCell(
        ghost.ghost_cell,
        HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index));
  }
}

}  // namespace cosmosim::hydro
