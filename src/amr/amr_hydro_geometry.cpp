#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"

#include "cosmosim/amr/amr_patch_indexing.hpp"
#include "cosmosim/core/memory_governor.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace cosmosim::amr {
namespace {

enum BoundarySlot : std::size_t {
  kLowerX = 0,
  kUpperX = 1,
  kLowerY = 2,
  kUpperY = 3,
  kLowerZ = 4,
  kUpperZ = 5,
};

[[nodiscard]] double checkedCellWidth(double extent_comov, std::uint16_t cell_dim, const char* axis) {
  if (extent_comov <= 0.0) {
    throw std::invalid_argument(std::string("AMR hydro patch extent must be positive on ") + axis);
  }
  return extent_comov / static_cast<double>(cell_dim);
}

[[nodiscard]] hydro::HydroGhostMutationRights mutationRightsForBoundary(
    AmrHydroBoundaryClass boundary_class) {
  return boundary_class == AmrHydroBoundaryClass::kPhysical
      ? hydro::HydroGhostMutationRights::kWritablePhysicalBoundaryScratch
      : hydro::HydroGhostMutationRights::kReadOnlyImported;
}

[[nodiscard]] hydro::HydroBoundaryKind hydroBoundaryKindForBoundary(
    AmrHydroBoundaryClass boundary_class,
    hydro::HydroBoundaryKind physical_boundary_kind) {
  return boundary_class == AmrHydroBoundaryClass::kPhysical
      ? physical_boundary_kind
      : hydro::HydroBoundaryKind::kImportedMpi;
}

[[nodiscard]] AmrHydroGhostFillStatus fillStatusForBoundary(AmrHydroBoundaryClass boundary_class) {
  switch (boundary_class) {
    case AmrHydroBoundaryClass::kPhysical:
      return AmrHydroGhostFillStatus::kUnfilledPhysicalBoundary;
    case AmrHydroBoundaryClass::kSameLevel:
      return AmrHydroGhostFillStatus::kUnfilledSameLevel;
    case AmrHydroBoundaryClass::kCoarseFine:
      return AmrHydroGhostFillStatus::kUnfilledCoarseFine;
  }
  return AmrHydroGhostFillStatus::kUnfilledPhysicalBoundary;
}

void setNormal(hydro::HydroFace& face, hydro::HydroFaceAxis axis, hydro::HydroFaceSide side) {
  const double sign = side == hydro::HydroFaceSide::kLower ? -1.0 : 1.0;
  face.normal_x = axis == hydro::HydroFaceAxis::kX ? sign : 0.0;
  face.normal_y = axis == hydro::HydroFaceAxis::kY ? sign : 0.0;
  face.normal_z = axis == hydro::HydroFaceAxis::kZ ? sign : 0.0;
}

[[nodiscard]] double faceArea(const hydro::HydroPatchGeometry& geometry, hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return geometry.cell_width_y_comoving * geometry.cell_width_z_comoving;
    case hydro::HydroFaceAxis::kY:
      return geometry.cell_width_x_comoving * geometry.cell_width_z_comoving;
    case hydro::HydroFaceAxis::kZ:
      return geometry.cell_width_x_comoving * geometry.cell_width_y_comoving;
  }
  return 0.0;
}

[[nodiscard]] std::size_t boundaryOwner(
    const hydro::HydroPatchGeometry& geometry,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX: {
      const std::size_t i = side == hydro::HydroFaceSide::kLower ? 0U : geometry.nx - 1U;
      return geometry.linearCellIndex(i, a, b);
    }
    case hydro::HydroFaceAxis::kY: {
      const std::size_t j = side == hydro::HydroFaceSide::kLower ? 0U : geometry.ny - 1U;
      return geometry.linearCellIndex(a, j, b);
    }
    case hydro::HydroFaceAxis::kZ: {
      const std::size_t k = side == hydro::HydroFaceSide::kLower ? 0U : geometry.nz - 1U;
      return geometry.linearCellIndex(a, b, k);
    }
  }
  return hydro::k_invalid_cell_index;
}

[[nodiscard]] BoundarySlot boundarySlot(hydro::HydroFaceAxis axis, hydro::HydroFaceSide side) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return side == hydro::HydroFaceSide::kLower ? kLowerX : kUpperX;
    case hydro::HydroFaceAxis::kY:
      return side == hydro::HydroFaceSide::kLower ? kLowerY : kUpperY;
    case hydro::HydroFaceAxis::kZ:
      return side == hydro::HydroFaceSide::kLower ? kLowerZ : kUpperZ;
  }
  return kLowerX;
}

[[nodiscard]] std::size_t wrappedPeriodicSource(
    const hydro::HydroPatchGeometry& geometry,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return geometry.linearCellIndex(side == hydro::HydroFaceSide::kLower ? geometry.nx - 1U : 0U, a, b);
    case hydro::HydroFaceAxis::kY:
      return geometry.linearCellIndex(a, side == hydro::HydroFaceSide::kLower ? geometry.ny - 1U : 0U, b);
    case hydro::HydroFaceAxis::kZ:
      return geometry.linearCellIndex(a, b, side == hydro::HydroFaceSide::kLower ? geometry.nz - 1U : 0U);
  }
  return hydro::k_invalid_cell_index;
}

[[nodiscard]] std::uint64_t sourceGasCellIdForGhost(
    const AmrHydroPatchGeometry& patch_geometry,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t a,
    std::size_t b,
    AmrHydroBoundaryClass boundary_class) {
  if (boundary_class != AmrHydroBoundaryClass::kPhysical) {
    return 0U;
  }

  const auto& geometry = patch_geometry.geometry;
  std::size_t source = hydro::k_invalid_cell_index;
  switch (axis) {
    case hydro::HydroFaceAxis::kX: {
      const std::size_t i = side == hydro::HydroFaceSide::kLower ? geometry.nx - 1U : 0U;
      source = geometry.linearCellIndex(i, a, b);
      break;
    }
    case hydro::HydroFaceAxis::kY: {
      const std::size_t j = side == hydro::HydroFaceSide::kLower ? geometry.ny - 1U : 0U;
      source = geometry.linearCellIndex(a, j, b);
      break;
    }
    case hydro::HydroFaceAxis::kZ: {
      const std::size_t k = side == hydro::HydroFaceSide::kLower ? geometry.nz - 1U : 0U;
      source = geometry.linearCellIndex(a, b, k);
      break;
    }
  }
  return patch_geometry.gas_cell_ids.at(source);
}

void appendBoundaryGhost(
    AmrHydroPatchGeometry& patch_geometry,
    const AmrHydroGeometryOptions& options,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  auto& geometry = patch_geometry.geometry;
  const AmrHydroBoundaryClass boundary_class =
      options.boundary_classes[static_cast<std::size_t>(boundarySlot(axis, side))];
  const std::size_t ghost_slot = geometry.ghost_cells.size();
  const std::size_t ghost_cell = geometry.cellCount() + ghost_slot;
  const std::size_t owner = boundaryOwner(geometry, axis, side, a, b);
  const std::uint64_t owner_gas_cell_id = patch_geometry.gas_cell_ids.at(owner);
  const std::uint64_t source_gas_cell_id =
      sourceGasCellIdForGhost(patch_geometry, axis, side, a, b, boundary_class);

  // Periodic physical boundary: the ghost is the wrapped cell of the SAME
  // patch (valid only when the patch spans the periodic box along this axis;
  // the orchestrator refuses other periodic topologies). Its flux is
  // deposited on that wrapped cell by HydroCoreSolver, so mass, momentum,
  // energy and metals cross the box boundary. Before this repair the source
  // was the owner itself, which made every periodic box face a zero-flux wall.
  const bool periodic_wrap = boundary_class == AmrHydroBoundaryClass::kPhysical &&
      options.physical_boundary_kind == hydro::HydroBoundaryKind::kPeriodic;
  const std::size_t wrapped_source = periodic_wrap
      ? wrappedPeriodicSource(geometry, axis, side, a, b)
      : owner;
  geometry.ghost_cells.push_back(hydro::HydroGhostCell{
      .owner_real_cell = owner,
      .source_real_cell = wrapped_source,
      .ghost_cell = ghost_cell,
      .ghost_slot = ghost_slot,
      .boundary_kind = hydroBoundaryKindForBoundary(boundary_class, options.physical_boundary_kind),
      .axis = axis,
      .side = side,
      .mutation_rights = mutationRightsForBoundary(boundary_class)});

  hydro::HydroFace face{
      .owner_cell = owner,
      .neighbor_cell = ghost_cell,
      .owner_minus_cell = hydro::k_invalid_cell_index,
      .neighbor_plus_cell = hydro::k_invalid_cell_index,
      .ghost_cell_slot = ghost_slot,
      .area_comoving = faceArea(geometry, axis),
      .axis = axis};
  setNormal(face, axis, side);
  geometry.faces.push_back(face);

  patch_geometry.ghosts.push_back(AmrHydroGhostDescriptor{
      .patch_id = patch_geometry.patch.patch_id,
      .source_patch_id = 0U,
      .source_gas_cell_id = source_gas_cell_id,
      .owner_gas_cell_id = owner_gas_cell_id,
      .owner_local_cell_row = patch_geometry.local_cell_rows.at(owner),
      .ghost_slot = ghost_slot,
      .ghost_cell = ghost_cell,
      .boundary_class = boundary_class,
      .fill_status = fillStatusForBoundary(boundary_class)});
}

void appendAxisBoundaryGhosts(
    AmrHydroPatchGeometry& patch_geometry,
    const AmrHydroGeometryOptions& options,
    hydro::HydroFaceAxis axis) {
  const auto& geometry = patch_geometry.geometry;
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      for (std::size_t k = 0; k < geometry.nz; ++k) {
        for (std::size_t j = 0; j < geometry.ny; ++j) {
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kLower, j, k);
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kUpper, j, k);
        }
      }
      break;
    case hydro::HydroFaceAxis::kY:
      for (std::size_t k = 0; k < geometry.nz; ++k) {
        for (std::size_t i = 0; i < geometry.nx; ++i) {
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kLower, i, k);
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kUpper, i, k);
        }
      }
      break;
    case hydro::HydroFaceAxis::kZ:
      for (std::size_t j = 0; j < geometry.ny; ++j) {
        for (std::size_t i = 0; i < geometry.nx; ++i) {
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kLower, i, j);
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kUpper, i, j);
        }
      }
      break;
  }
}

[[nodiscard]] AmrHydroFaceDescriptor makeFaceDescriptor(
    const AmrHydroPatchGeometry& patch_geometry,
    std::size_t face_index) {
  const hydro::HydroFace& face = patch_geometry.geometry.faces.at(face_index);
  AmrHydroBoundaryClass boundary_class = AmrHydroBoundaryClass::kPhysical;
  std::uint64_t ghost_source_gas_cell_id = 0U;
  if (face.ghost_cell_slot != hydro::k_invalid_ghost_cell_slot) {
    const auto ghost_it = std::find_if(
        patch_geometry.ghosts.begin(),
        patch_geometry.ghosts.end(),
        [&face](const AmrHydroGhostDescriptor& ghost) {
          return ghost.ghost_slot == face.ghost_cell_slot;
        });
    if (ghost_it != patch_geometry.ghosts.end()) {
      boundary_class = ghost_it->boundary_class;
      ghost_source_gas_cell_id = ghost_it->source_gas_cell_id;
    }
  }

  return AmrHydroFaceDescriptor{
      .patch_id = patch_geometry.patch.patch_id,
      .face_id = (patch_geometry.patch.patch_id << 32U) | static_cast<std::uint64_t>(face_index),
      .owner_gas_cell_id = patch_geometry.gas_cell_ids.at(face.owner_cell),
      .neighbor_gas_cell_id = face.neighbor_cell < patch_geometry.gas_cell_ids.size()
          ? patch_geometry.gas_cell_ids[face.neighbor_cell]
          : 0U,
      .ghost_source_gas_cell_id = ghost_source_gas_cell_id,
      .owner_patch_cell = face.owner_cell,
      .neighbor_patch_cell = face.neighbor_cell,
      .ghost_slot = face.ghost_cell_slot,
      .area_comoving = face.area_comoving,
      .normal_x = face.normal_x,
      .normal_y = face.normal_y,
      .normal_z = face.normal_z,
      .axis = face.axis,
      .boundary_class = boundary_class};
}

[[nodiscard]] double positiveDensityForCell(const core::SimulationState& state, std::uint32_t row) {
  const double density = state.gas_cells.density_code.at(row);
  if (density > 0.0) {
    return density;
  }
  const double volume = 1.0;
  return std::max(state.cells.mass_code.at(row) / volume, 1.0e-14);
}


// ---------------------------------------------------------------------------
// CHUI Part 3 Checkpoint A: deep (stencil-only) ghost layers.
//
// The one-layer shell produced by appendAxisBoundaryGhosts gives every boundary
// face its neighbour ghost (layer 1). Layers 2..N are appended AFTER the
// standard shell, so the standard slot formula used by the hydro stencil
// helpers is untouched. Deep ghosts own no face: they are read by
// reconstruction only, never advanced and never deposit flux.
//
// Source cell of a deep PHYSICAL ghost (filled by
// hydro::fillHydroBoundaryGhostCells):
//   open        -> the boundary owner cell (zero-gradient extrapolation),
//   reflective  -> mirror cell at depth d-1 inside the patch (normal velocity
//                  reflected by the filler). When the patch is narrower than
//                  d cells along the axis the mirror lies beyond the opposite
//                  face: the layer is created with a provisional owner source
//                  and physical_closure_mapped = true, and the AMR ghost
//                  filler resolves the exact mirror (P2-55-03),
//   periodic    -> never created (the stencil already wraps to real cells).
// Deep NON-PHYSICAL ghosts (same-level, coarse/fine, remote) are filled by the
// AMR ghost filler from the source patch cell containing the probe at
// owner_center + d * width * outward_normal (wrapped into the periodic box).
// ---------------------------------------------------------------------------
struct DeepGhostLayerInfo {
  std::size_t standard_slot = hydro::k_invalid_ghost_cell_slot;
  std::uint8_t layer = 1;
};

[[nodiscard]] std::size_t axisDimension(const hydro::HydroPatchGeometry& geometry, hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX: return geometry.nx;
    case hydro::HydroFaceAxis::kY: return geometry.ny;
    case hydro::HydroFaceAxis::kZ: return geometry.nz;
  }
  return 0U;
}

[[nodiscard]] std::size_t mirrorCellAtDepth(
    const hydro::HydroPatchGeometry& geometry,
    std::size_t owner,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t depth) {
  const std::array<std::size_t, 3> ijk = geometry.cellIjk(owner);
  std::array<std::size_t, 3> mirrored = ijk;
  const std::size_t axis_index = axis == hydro::HydroFaceAxis::kX ? 0U : (axis == hydro::HydroFaceAxis::kY ? 1U : 2U);
  const std::size_t n = axisDimension(geometry, axis);
  if (depth == 0U || depth > n) {
    return hydro::k_invalid_cell_index;
  }
  mirrored[axis_index] = side == hydro::HydroFaceSide::kLower ? depth - 1U : n - depth;
  return geometry.linearCellIndex(mirrored[0], mirrored[1], mirrored[2]);
}

void appendDeepGhostLayers(
    AmrHydroPatchGeometry& patch_geometry,
    const AmrHydroGeometryOptions& options,
    std::vector<DeepGhostLayerInfo>& layer_info) {
  auto& geometry = patch_geometry.geometry;
  const std::size_t standard_count = geometry.ghost_cells.size();
  layer_info.assign(standard_count, DeepGhostLayerInfo{});
  for (std::size_t slot = 0; slot < standard_count; ++slot) {
    layer_info[slot] = DeepGhostLayerInfo{.standard_slot = slot, .layer = 1U};
  }
  if (options.ghost_layers <= 1U) {
    return;
  }
  for (std::size_t standard_slot = 0; standard_slot < standard_count; ++standard_slot) {
    const hydro::HydroGhostCell standard = geometry.ghost_cells[standard_slot];
    const AmrHydroGhostDescriptor standard_descriptor = patch_geometry.ghosts.at(standard_slot);
    const bool physical = standard_descriptor.boundary_class == AmrHydroBoundaryClass::kPhysical;
    if (physical && options.physical_boundary_kind == hydro::HydroBoundaryKind::kPeriodic) {
      continue;
    }
    for (std::uint8_t layer = 2U; layer <= options.ghost_layers; ++layer) {
      std::size_t source = standard.owner_real_cell;
      bool closure_mapped = false;
      if (physical && options.physical_boundary_kind == hydro::HydroBoundaryKind::kReflective) {
        source = mirrorCellAtDepth(geometry, standard.owner_real_cell, standard.axis, standard.side, layer);
        if (source == hydro::k_invalid_cell_index) {
          // Provisional, always-valid source for the boundary filler; the
          // AMR ghost filler overwrites it with the exact mirror state.
          source = standard.owner_real_cell;
          closure_mapped = true;
        }
      }
      const std::size_t ghost_slot = geometry.ghost_cells.size();
      const std::size_t ghost_cell = geometry.cellCount() + ghost_slot;
      hydro::HydroGhostCell deep = standard;
      deep.source_real_cell = source;
      deep.ghost_cell = ghost_cell;
      deep.ghost_slot = ghost_slot;
      geometry.ghost_cells.push_back(deep);
      AmrHydroGhostDescriptor descriptor = standard_descriptor;
      descriptor.ghost_slot = ghost_slot;
      descriptor.ghost_cell = ghost_cell;
      descriptor.ghost_layer = layer;
      descriptor.source_patch_id = 0U;
      descriptor.source_gas_cell_id = physical && !closure_mapped ? patch_geometry.gas_cell_ids.at(source) : 0U;
      descriptor.physical_closure_mapped = closure_mapped;
      patch_geometry.ghosts.push_back(descriptor);
      layer_info.push_back(DeepGhostLayerInfo{.standard_slot = standard_slot, .layer = layer});
    }
  }
}

// Fill still-unavailable stencil entries that step outward from a ghost onto
// the next deeper layer of the same standard ghost column.
void linkDeepGhostStencils(
    AmrHydroPatchGeometry& patch_geometry,
    const std::vector<DeepGhostLayerInfo>& layer_info) {
  auto& geometry = patch_geometry.geometry;
  const std::size_t real_cells = geometry.cellCount();
  if (std::all_of(layer_info.begin(), layer_info.end(),
                  [](const DeepGhostLayerInfo& info) { return info.layer == 1U; })) {
    return;
  }
  std::unordered_map<std::uint64_t, std::size_t> deeper;
  deeper.reserve(layer_info.size());
  for (std::size_t slot = 0; slot < layer_info.size(); ++slot) {
    deeper.emplace((static_cast<std::uint64_t>(layer_info[slot].standard_slot) << 8U) |
                       static_cast<std::uint64_t>(layer_info[slot].layer),
                   real_cells + slot);
  }
  const auto outward = [&](std::size_t cell, hydro::HydroFaceAxis axis, int dir) -> std::size_t {
    if (cell == hydro::k_invalid_cell_index || cell < real_cells) {
      return hydro::k_invalid_cell_index;
    }
    const std::size_t slot = cell - real_cells;
    if (slot >= geometry.ghost_cells.size()) {
      return hydro::k_invalid_cell_index;
    }
    const hydro::HydroGhostCell& ghost = geometry.ghost_cells[slot];
    if (ghost.axis != axis) {
      return hydro::k_invalid_cell_index;
    }
    const bool moving_away = (ghost.side == hydro::HydroFaceSide::kUpper) == (dir > 0);
    if (!moving_away) {
      return hydro::k_invalid_cell_index;
    }
    const DeepGhostLayerInfo& info = layer_info.at(slot);
    const auto it = deeper.find((static_cast<std::uint64_t>(info.standard_slot) << 8U) |
                                (static_cast<std::uint64_t>(info.layer) + 1U));
    return it == deeper.end() ? hydro::k_invalid_cell_index : it->second;
  };
  const auto normalSign = [](const hydro::HydroFace& face) {
    switch (face.axis) {
      case hydro::HydroFaceAxis::kX: return face.normal_x >= 0.0 ? 1 : -1;
      case hydro::HydroFaceAxis::kY: return face.normal_y >= 0.0 ? 1 : -1;
      case hydro::HydroFaceAxis::kZ: return face.normal_z >= 0.0 ? 1 : -1;
    }
    return 1;
  };
  for (hydro::HydroFace& face : geometry.faces) {
    const int along = normalSign(face);
    if (face.owner_minus2_cell == hydro::k_invalid_cell_index) {
      face.owner_minus2_cell = outward(face.owner_minus_cell, face.axis, -along);
    }
    if (face.neighbor_plus_cell == hydro::k_invalid_cell_index) {
      face.neighbor_plus_cell = outward(face.neighbor_cell, face.axis, along);
    }
    if (face.neighbor_plus2_cell == hydro::k_invalid_cell_index) {
      face.neighbor_plus2_cell = outward(face.neighbor_plus_cell, face.axis, along);
    }
  }
}

// ---------------------------------------------------------------------------
// P1-HYDRO-02: edge stencil ghost block.
//
// The CTU transverse predictor of a face ghost g (normal axis n) needs the
// tangential slab g + o e_t, |o| <= stencil half-width, for both transverse
// axes t. Where g + o e_t leaves the patch along t, the slab cell lies beyond
// the patch along BOTH n and t (a face-edge region; with two such axes the
// patch corner region is covered by the two one-dimensional slabs, so no
// three-axis corner cell is ever read). Those cells are keyed by position:
// (n, normal side, t, transverse side, depth along t, index along the third
// axis), appended after the deep layers in the order of
// hydro::HydroPatchGeometry::edgeGhostSlot. Each is stencil-only and
// read-only (kReadOnlyImported, kImportedMpi), so the physical boundary
// filler never writes it and no face or flux register refers to it.
// ---------------------------------------------------------------------------
[[nodiscard]] std::size_t standardShellSlot(
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

void appendEdgeStencilGhosts(
    AmrHydroPatchGeometry& patch_geometry,
    const AmrHydroGeometryOptions& options) {
  auto& geometry = patch_geometry.geometry;
  const std::size_t depth = options.edge_stencil_depth;
  if (depth == 0U) {
    return;
  }
  if (depth > static_cast<std::size_t>(std::numeric_limits<std::int8_t>::max())) {
    throw std::invalid_argument("buildAmrHydroPatchGeometry: edge stencil depth exceeds the int8 offset range");
  }
  const std::array<std::size_t, 3> extent{geometry.nx, geometry.ny, geometry.nz};
  const std::size_t real_cells = geometry.cellCount();
  const std::size_t expected = hydro::hydroEdgeGhostCount(geometry.nx, geometry.ny, geometry.nz, depth);
  if (core::checkedSizeAdd(geometry.ghost_cells.size(), expected, "edge ghost capacity") >
          geometry.ghost_cells.capacity() ||
      patch_geometry.ghosts.size() + expected > patch_geometry.ghosts.capacity()) {
    throw std::length_error("buildAmrHydroPatchGeometry: edge stencil ghosts exceed the admitted capacity");
  }
  geometry.edge_ghost_slot_begin = geometry.ghost_cells.size();
  geometry.edge_ghost_depth = static_cast<std::uint8_t>(depth);
  constexpr std::array<hydro::HydroFaceAxis, 3> k_axes{
      hydro::HydroFaceAxis::kX, hydro::HydroFaceAxis::kY, hydro::HydroFaceAxis::kZ};
  for (std::size_t n = 0U; n < 3U; ++n) {
    for (std::size_t side = 0U; side < 2U; ++side) {
      for (std::size_t t = 0U; t < 3U; ++t) {
        if (t == n || extent[t] <= 1U) continue;
        const std::size_t u = 3U - n - t;
        for (std::size_t transverse_side = 0U; transverse_side < 2U; ++transverse_side) {
          for (std::size_t d = 1U; d <= depth; ++d) {
            for (std::size_t third = 0U; third < extent[u]; ++third) {
              std::array<std::size_t, 3> ijk{};
              ijk[n] = side == 0U ? 0U : extent[n] - 1U;
              ijk[t] = transverse_side == 0U ? 0U : extent[t] - 1U;
              ijk[u] = third;
              const std::size_t owner = geometry.linearCellIndex(ijk[0], ijk[1], ijk[2]);
              const std::size_t standard_slot = standardShellSlot(geometry, n, side, ijk);
              const hydro::HydroGhostCell standard = geometry.ghost_cells.at(standard_slot);
              if (standard.owner_real_cell != owner || standard.axis != k_axes[n] ||
                  standard.side != (side == 0U ? hydro::HydroFaceSide::kLower : hydro::HydroFaceSide::kUpper)) {
                throw std::logic_error("buildAmrHydroPatchGeometry: standard ghost shell layout mismatch");
              }
              const AmrHydroGhostDescriptor standard_descriptor = patch_geometry.ghosts.at(standard_slot);
              const std::size_t slot = geometry.ghost_cells.size();
              hydro::HydroGhostCell edge = standard;
              edge.source_real_cell = owner;  // metadata only; never read (read-only imported)
              edge.ghost_slot = slot;
              edge.ghost_cell = real_cells + slot;
              edge.boundary_kind = hydro::HydroBoundaryKind::kImportedMpi;
              edge.mutation_rights = hydro::HydroGhostMutationRights::kReadOnlyImported;
              edge.origin_gas_cell_id = 0U;
              edge.origin_rank = -1;
              edge.transverse_axis = static_cast<std::uint8_t>(t);
              edge.transverse_offset = static_cast<std::int8_t>(
                  transverse_side == 0U ? -static_cast<int>(d) : static_cast<int>(d));
              geometry.ghost_cells.push_back(edge);
              AmrHydroGhostDescriptor descriptor = standard_descriptor;
              descriptor.ghost_slot = slot;
              descriptor.ghost_cell = edge.ghost_cell;
              descriptor.ghost_layer = 1U;
              descriptor.source_patch_id = 0U;
              descriptor.source_gas_cell_id = 0U;
              descriptor.edge_stencil = true;
              patch_geometry.ghosts.push_back(descriptor);
              if (geometry.edgeGhostSlot(n, side, t, transverse_side, d, third) != slot) {
                throw std::logic_error("buildAmrHydroPatchGeometry: edge ghost slot formula mismatch");
              }
            }
          }
        }
      }
    }
  }
  if (geometry.ghost_cells.size() - geometry.edge_ghost_slot_begin != expected) {
    throw std::logic_error("buildAmrHydroPatchGeometry: edge ghost count mismatch");
  }
}

}  // namespace

AmrHydroGeometryCapacity amrHydroGeometryCapacity(
    const PatchDescriptor& patch, std::uint8_t ghost_layers, std::uint8_t edge_stencil_depth) {
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t nz = patch.cell_dims[2];
  const std::size_t n = checkedPatchCellCount(patch, "AMR geometry capacity");
  const auto mul = [](std::size_t a, std::size_t b) {
    return core::checkedSizeMultiply(a, b, "AMR geometry topology");
  };
  const auto add = [](std::size_t a, std::size_t b) {
    return core::checkedSizeAdd(a, b, "AMR geometry topology");
  };
  const std::size_t surface = add(add(mul(nx, ny), mul(nx, nz)), mul(ny, nz));
  const std::size_t face_ghosts = mul(2U, surface);
  const std::size_t ghosts = add(mul(face_ghosts, std::max<std::size_t>(1U, ghost_layers)),
      hydro::hydroEdgeGhostCount(nx, ny, nz, edge_stencil_depth));
  const std::size_t internal = add(add(mul(nx - 1U, mul(ny, nz)),
      mul(ny - 1U, mul(nx, nz))), mul(nz - 1U, mul(nx, ny)));
  const std::size_t faces = add(internal, face_ghosts);
  const auto bytes = [](std::size_t count, std::size_t width) {
    return static_cast<std::uint64_t>(core::checkedSizeMultiply(
        count, width, "AMR geometry storage"));
  };
  std::uint64_t retained = 0U;
  for (const std::size_t width : {sizeof(AmrHydroCellDescriptor),
       sizeof(std::uint64_t), sizeof(std::uint32_t)}) {
    retained = core::checkedMemoryBytesAdd(retained, bytes(n, width), "AMR geometry storage");
  }
  for (const std::size_t width : {sizeof(hydro::HydroFace),
       sizeof(hydro::HydroFluxRegisterFace), sizeof(AmrHydroFaceDescriptor)}) {
    retained = core::checkedMemoryBytesAdd(retained, bytes(faces, width), "AMR geometry storage");
  }
  for (const std::size_t width : {sizeof(hydro::HydroGhostCell),
       sizeof(AmrHydroGhostDescriptor)}) {
    retained = core::checkedMemoryBytesAdd(retained, bytes(ghosts, width), "AMR geometry storage");
  }
  return {n, ghosts, faces, retained,
      bytes(n, 2U * sizeof(std::uint32_t))};
}

std::uint64_t AmrHydroPatchGeometry::ownedCapacityBytes() const {
  std::uint64_t bytes = 0U;
  const auto add = [&bytes](std::size_t capacity, std::size_t width) {
    bytes = core::checkedMemoryBytesAdd(bytes,
        static_cast<std::uint64_t>(core::checkedSizeMultiply(
            capacity, width, "AMR geometry physical capacity")),
        "AMR geometry physical capacity");
  };
  add(real_cells.capacity(), sizeof(AmrHydroCellDescriptor));
  add(gas_cell_ids.capacity(), sizeof(std::uint64_t));
  add(local_cell_rows.capacity(), sizeof(std::uint32_t));
  add(ghosts.capacity(), sizeof(AmrHydroGhostDescriptor));
  add(faces.capacity(), sizeof(AmrHydroFaceDescriptor));
  add(geometry.faces.capacity(), sizeof(hydro::HydroFace));
  add(geometry.ghost_cells.capacity(), sizeof(hydro::HydroGhostCell));
  add(geometry.flux_register_faces.capacity(), sizeof(hydro::HydroFluxRegisterFace));
  return bytes;
}

std::span<const std::uint64_t> AmrHydroPatchGeometry::gasCellIds() const noexcept {
  return gas_cell_ids;
}

std::vector<std::size_t> AmrHydroPatchGeometry::internalFaceIndices() const {
  std::vector<std::size_t> indices;
  indices.reserve(faces.size());
  for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
    if (faces[face_index].ghost_slot == hydro::k_invalid_ghost_cell_slot) {
      indices.push_back(face_index);
    }
  }
  return indices;
}

AmrHydroPatchGeometry buildAmrHydroPatchGeometry(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    const AmrHydroGeometryOptions& options) {
  const std::size_t expected_cells = checkedPatchCellCount(patch, "buildAmrHydroPatchGeometry");
  state.requireGasCellIdentityMapCoversDenseRows("buildAmrHydroPatchGeometry");
  if (!state.gasCellIdentityMapMatchesSidecarLanes()) {
    throw std::runtime_error("buildAmrHydroPatchGeometry: gas-cell identity sidecar mirrors are stale");
  }

  // The descriptor gives the exact local row count. Avoid a second global
  // scan merely to discover capacity, and never grow a patch row list through
  // an unadmitted sequence of geometric reallocations.
  std::vector<std::uint32_t> rows;
  rows.reserve(expected_cells);
  for (const core::GasCellIdentityRecord& record : state.gas_cell_identity.records()) {
    if (record.owning_patch_id == patch.patch_id) rows.push_back(record.local_cell_row);
  }
  std::sort(rows.begin(), rows.end());
  if (rows.size() != expected_cells) {
    throw std::runtime_error("buildAmrHydroPatchGeometry: patch gas-cell coverage does not match PatchDescriptor cell_dims");
  }
  const std::vector<std::uint32_t> row_by_patch_cell = buildPatchLocalRowMap(
      state, patch, std::span<const std::uint32_t>(rows.data(), rows.size()), "buildAmrHydroPatchGeometry");

  AmrHydroPatchGeometry result;
  result.patch = patch;
  result.source_gas_cell_identity_generation = state.gasCellIdentityGeneration();
  result.geometry.nx = patch.cell_dims[0];
  result.geometry.ny = patch.cell_dims[1];
  result.geometry.nz = patch.cell_dims[2];
  result.geometry.origin_x_comoving = patch.origin_comov[0];
  result.geometry.origin_y_comoving = patch.origin_comov[1];
  result.geometry.origin_z_comoving = patch.origin_comov[2];
  result.geometry.cell_width_x_comoving = checkedCellWidth(patch.extent_comov[0], patch.cell_dims[0], "x");
  result.geometry.cell_width_y_comoving = checkedCellWidth(patch.extent_comov[1], patch.cell_dims[1], "y");
  result.geometry.cell_width_z_comoving = checkedCellWidth(patch.extent_comov[2], patch.cell_dims[2], "z");
  result.geometry.cell_volume_comoving =
      result.geometry.cell_width_x_comoving *
      result.geometry.cell_width_y_comoving *
      result.geometry.cell_width_z_comoving;

  const auto capacity = amrHydroGeometryCapacity(patch, options.ghost_layers, options.edge_stencil_depth);
  result.geometry.faces.reserve(capacity.faces);
  result.geometry.ghost_cells.reserve(capacity.ghost_cells);
  result.geometry.flux_register_faces.reserve(capacity.faces);
  result.ghosts.reserve(capacity.ghost_cells);
  result.faces.reserve(capacity.faces);
  result.real_cells.reserve(row_by_patch_cell.size());
  result.gas_cell_ids.reserve(row_by_patch_cell.size());
  result.local_cell_rows.reserve(row_by_patch_cell.size());
  for (std::size_t patch_cell = 0; patch_cell < row_by_patch_cell.size(); ++patch_cell) {
    const std::uint32_t row = row_by_patch_cell[patch_cell];
    const auto* record = state.gas_cell_identity.findByLocalRow(row);
    if (record == nullptr || record->owning_patch_id != patch.patch_id) {
      throw std::runtime_error("buildAmrHydroPatchGeometry: identity row does not belong to requested patch");
    }
    if (row >= state.cells.size() || row >= state.gas_cells.size()) {
      throw std::runtime_error("buildAmrHydroPatchGeometry: identity row is outside gas-cell storage");
    }
    result.real_cells.push_back(AmrHydroCellDescriptor{
        .patch_id = patch.patch_id,
        .gas_cell_id = record->gas_cell_id,
        .local_cell_row = row,
        .patch_local_cell = patch_cell});
    result.gas_cell_ids.push_back(record->gas_cell_id);
    result.local_cell_rows.push_back(row);
  }

  const auto append_internal_face =
      [&result](std::size_t owner, std::size_t neighbor, hydro::HydroFaceAxis axis) {
        hydro::HydroFace face{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = hydro::k_invalid_cell_index,
            .neighbor_plus_cell = hydro::k_invalid_cell_index,
            .area_comoving = faceArea(result.geometry, axis),
            .axis = axis};
        switch (axis) {
          case hydro::HydroFaceAxis::kX:
            face.owner_minus_cell = result.geometry.neighborCell(owner, -1, 0, 0);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 1, 0, 0);
            face.owner_minus2_cell = result.geometry.neighborCell(owner, -2, 0, 0);
            face.neighbor_plus2_cell = result.geometry.neighborCell(neighbor, 2, 0, 0);
            face.normal_x = 1.0;
            break;
          case hydro::HydroFaceAxis::kY:
            face.owner_minus_cell = result.geometry.neighborCell(owner, 0, -1, 0);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 0, 1, 0);
            face.owner_minus2_cell = result.geometry.neighborCell(owner, 0, -2, 0);
            face.neighbor_plus2_cell = result.geometry.neighborCell(neighbor, 0, 2, 0);
            face.normal_y = 1.0;
            break;
          case hydro::HydroFaceAxis::kZ:
            face.owner_minus_cell = result.geometry.neighborCell(owner, 0, 0, -1);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 0, 0, 1);
            face.owner_minus2_cell = result.geometry.neighborCell(owner, 0, 0, -2);
            face.neighbor_plus2_cell = result.geometry.neighborCell(neighbor, 0, 0, 2);
            face.normal_z = 1.0;
            break;
        }
        result.geometry.faces.push_back(face);
      };

  for (std::size_t k = 0; k < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i + 1U < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + 1U, hydro::HydroFaceAxis::kX);
      }
    }
  }
  for (std::size_t k = 0; k < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j + 1U < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + result.geometry.nx, hydro::HydroFaceAxis::kY);
      }
    }
  }
  for (std::size_t k = 0; k + 1U < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + result.geometry.nx * result.geometry.ny, hydro::HydroFaceAxis::kZ);
      }
    }
  }

  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kX);
  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kY);
  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kZ);
  // P2-55-04: the face-adjacent shell is the only face-carrying ghost prefix;
  // deep layers and the edge block below are stencil-only.
  result.geometry.face_ghost_slot_count = result.geometry.ghost_cells.size();
  // Reconstruction stencils through the standard one-layer shell (periodic
  // wraps stay in-patch), then the deep stencil-only layers (Part 3) give
  // full PLM/PPM order across patch/rank/coarse-fine/open/reflective faces.
  (void)hydro::completeReconstructionStencils(result.geometry);
  {
    std::vector<DeepGhostLayerInfo> layer_info;
    appendDeepGhostLayers(result, options, layer_info);
    linkDeepGhostStencils(result, layer_info);
  }
  appendEdgeStencilGhosts(result, options);

  result.faces.reserve(result.geometry.faces.size());
  for (std::size_t face_index = 0; face_index < result.geometry.faces.size(); ++face_index) {
    result.faces.push_back(makeFaceDescriptor(result, face_index));
  }
  return result;
}

AmrHydroPatchGeometry buildRemoteAmrHydroPatchGeometry(
    const PatchDescriptor& patch,
    std::span<const std::uint64_t> gas_cell_ids,
    std::uint64_t source_gas_cell_identity_generation,
    const AmrHydroGeometryOptions& options,
    std::span<const std::uint8_t> available_real_cells) {
  const std::size_t expected_cells = checkedPatchCellCount(patch, "buildRemoteAmrHydroPatchGeometry");
  if (gas_cell_ids.size() != expected_cells) {
    throw std::runtime_error("buildRemoteAmrHydroPatchGeometry: remote gas-cell coverage does not match patch geometry");
  }
  if (!available_real_cells.empty() && available_real_cells.size() != expected_cells) {
    throw std::runtime_error("buildRemoteAmrHydroPatchGeometry: sparse remote availability mask does not match patch geometry");
  }
  AmrHydroPatchGeometry result;
  result.patch = patch;
  result.source_gas_cell_identity_generation = source_gas_cell_identity_generation;
  result.geometry.nx = patch.cell_dims[0];
  result.geometry.ny = patch.cell_dims[1];
  result.geometry.nz = patch.cell_dims[2];
  result.geometry.origin_x_comoving = patch.origin_comov[0];
  result.geometry.origin_y_comoving = patch.origin_comov[1];
  result.geometry.origin_z_comoving = patch.origin_comov[2];
  result.geometry.cell_width_x_comoving = checkedCellWidth(patch.extent_comov[0], patch.cell_dims[0], "x");
  result.geometry.cell_width_y_comoving = checkedCellWidth(patch.extent_comov[1], patch.cell_dims[1], "y");
  result.geometry.cell_width_z_comoving = checkedCellWidth(patch.extent_comov[2], patch.cell_dims[2], "z");
  result.geometry.cell_volume_comoving =
      result.geometry.cell_width_x_comoving *
      result.geometry.cell_width_y_comoving *
      result.geometry.cell_width_z_comoving;
  const auto capacity = amrHydroGeometryCapacity(patch, options.ghost_layers);
  result.geometry.faces.reserve(capacity.faces);
  result.geometry.ghost_cells.reserve(capacity.ghost_cells);
  result.geometry.flux_register_faces.reserve(capacity.faces);
  result.ghosts.reserve(capacity.ghost_cells);
  result.faces.reserve(capacity.faces);
  result.real_cells.reserve(expected_cells);
  result.gas_cell_ids.assign(gas_cell_ids.begin(), gas_cell_ids.end());
  result.local_cell_rows.assign(expected_cells, core::kInvalidGasCellRow);
  for (std::size_t patch_cell = 0; patch_cell < expected_cells; ++patch_cell) {
    const bool available = available_real_cells.empty() || available_real_cells[patch_cell] != 0U;
    if (available && result.gas_cell_ids[patch_cell] == 0U) {
      throw std::runtime_error("buildRemoteAmrHydroPatchGeometry: available remote gas-cell IDs must be non-zero");
    }
    if (!available && result.gas_cell_ids[patch_cell] != 0U) {
      throw std::runtime_error("buildRemoteAmrHydroPatchGeometry: unavailable sparse remote cell unexpectedly carries an ID");
    }
    result.real_cells.push_back(AmrHydroCellDescriptor{
        .patch_id = patch.patch_id,
        .gas_cell_id = result.gas_cell_ids[patch_cell],
        .local_cell_row = core::kInvalidGasCellRow,
        .patch_local_cell = patch_cell});
  }

  const auto append_internal_face =
      [&result](std::size_t owner, std::size_t neighbor, hydro::HydroFaceAxis axis) {
        hydro::HydroFace face{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = hydro::k_invalid_cell_index,
            .neighbor_plus_cell = hydro::k_invalid_cell_index,
            .area_comoving = faceArea(result.geometry, axis),
            .axis = axis};
        switch (axis) {
          case hydro::HydroFaceAxis::kX:
            face.owner_minus_cell = result.geometry.neighborCell(owner, -1, 0, 0);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 1, 0, 0);
            face.owner_minus2_cell = result.geometry.neighborCell(owner, -2, 0, 0);
            face.neighbor_plus2_cell = result.geometry.neighborCell(neighbor, 2, 0, 0);
            face.normal_x = 1.0;
            break;
          case hydro::HydroFaceAxis::kY:
            face.owner_minus_cell = result.geometry.neighborCell(owner, 0, -1, 0);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 0, 1, 0);
            face.owner_minus2_cell = result.geometry.neighborCell(owner, 0, -2, 0);
            face.neighbor_plus2_cell = result.geometry.neighborCell(neighbor, 0, 2, 0);
            face.normal_y = 1.0;
            break;
          case hydro::HydroFaceAxis::kZ:
            face.owner_minus_cell = result.geometry.neighborCell(owner, 0, 0, -1);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 0, 0, 1);
            face.owner_minus2_cell = result.geometry.neighborCell(owner, 0, 0, -2);
            face.neighbor_plus2_cell = result.geometry.neighborCell(neighbor, 0, 0, 2);
            face.normal_z = 1.0;
            break;
        }
        result.geometry.faces.push_back(face);
      };

  for (std::size_t k = 0; k < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i + 1U < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + 1U, hydro::HydroFaceAxis::kX);
      }
    }
  }
  for (std::size_t k = 0; k < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j + 1U < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + result.geometry.nx, hydro::HydroFaceAxis::kY);
      }
    }
  }
  for (std::size_t k = 0; k + 1U < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + result.geometry.nx * result.geometry.ny, hydro::HydroFaceAxis::kZ);
      }
    }
  }
  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kX);
  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kY);
  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kZ);
  // P2-55-04: the face-adjacent shell is the only face-carrying ghost prefix;
  // deep layers and the edge block below are stencil-only.
  result.geometry.face_ghost_slot_count = result.geometry.ghost_cells.size();
  // Reconstruction stencils through the standard one-layer shell (periodic
  // wraps stay in-patch), then the deep stencil-only layers (Part 3) give
  // full PLM/PPM order across patch/rank/coarse-fine/open/reflective faces.
  (void)hydro::completeReconstructionStencils(result.geometry);
  {
    std::vector<DeepGhostLayerInfo> layer_info;
    appendDeepGhostLayers(result, options, layer_info);
    linkDeepGhostStencils(result, layer_info);
  }
  result.faces.reserve(result.geometry.faces.size());
  for (std::size_t face_index = 0; face_index < result.geometry.faces.size(); ++face_index) {
    result.faces.push_back(makeFaceDescriptor(result, face_index));
  }
  return result;
}

void loadAmrHydroConservedStateInto(
    const core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    double adiabatic_index,
    hydro::HydroConservedStateSoa& conserved) {
  state.requireGasCellIdentityMapFresh(
      patch_geometry.source_gas_cell_identity_generation,
      "loadAmrHydroConservedState");
  conserved.resize(patch_geometry.geometry.totalCellStorageCount());
  for (std::size_t patch_cell = 0; patch_cell < patch_geometry.real_cells.size(); ++patch_cell) {
    const std::uint32_t row = patch_geometry.real_cells[patch_cell].local_cell_row;
    // Thermal authority is the canonical specific internal energy (same
    // contract as the fixed-grid path and the refresh helper). The derived
    // pressure cache is only a legacy fallback for states that never set e;
    // using it as authority would mis-build E for effective-EOS cells.
    const double rho = positiveDensityForCell(state, row);
    const double canonical_e = state.gas_cells.internal_energy_code.at(row);
    const double pressure = canonical_e > 0.0
        ? std::max((adiabatic_index - 1.0) * rho * canonical_e, 1.0e-14)
        : std::max(state.gas_cells.pressure_code.at(row), 1.0e-14);
    const hydro::HydroPrimitiveState primitive{
        .rho_comoving = rho,
        .vel_x_peculiar = state.gas_cells.velocity_x_peculiar.at(row),
        .vel_y_peculiar = state.gas_cells.velocity_y_peculiar.at(row),
        .vel_z_peculiar = state.gas_cells.velocity_z_peculiar.at(row),
        .pressure_comoving = pressure,
        .metallicity_mass_fraction = std::clamp(
            state.gas_cells.metal_mass_code.at(row) /
                std::max(state.cells.mass_code.at(row), 1.0e-30),
            0.0,
            1.0)};
    conserved.storeCell(
        patch_cell,
        hydro::HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index));
  }
}

hydro::HydroConservedStateSoa loadAmrHydroConservedState(
    const core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    double adiabatic_index) {
  hydro::HydroConservedStateSoa conserved;
  loadAmrHydroConservedStateInto(state, patch_geometry, adiabatic_index, conserved);
  return conserved;
}

}  // namespace cosmosim::amr
