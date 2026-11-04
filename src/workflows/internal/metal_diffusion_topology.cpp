#include "workflows/internal/metal_diffusion_topology.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace cosmosim::workflows::internal {
namespace {

// Per-node allowances for node-based std::map/std::set scratch (libstdc++/
// libc++ red-black node: 3 pointers + colour + payload, rounded to the
// allocator's 16-byte granularity). A model, not measured allocator usage.
constexpr std::size_t k_map_node_bytes = 96U;   // key (16) + vector header (24) + node (48) + pad
constexpr std::size_t k_set_node_bytes = 64U;   // tuple key (12) + node (48) + pad

struct DirectedFaceRef {
  std::uint32_t cell = 0U;
  std::uint32_t face_index = 0U;
};

struct PatchBoundaryFace {
  std::uint32_t cell = 0U;
  std::uint8_t axis = 0U;
  std::int8_t side = 0;
  double plane = 0.0;
  double transverse_u_min = 0.0;
  double transverse_u_max = 0.0;
  double transverse_v_min = 0.0;
  double transverse_v_max = 0.0;
  double normal_width = 0.0;
};

struct NeighborAccumulator {
  double velocity_area_sum[3]{};
  double distance_area_sum = 0.0;
  double area_sum = 0.0;
};

[[nodiscard]] bool patchHasOwnedLeafCells(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t patch_index,
    std::uint32_t world_rank) {
  if (patch_index >= state.patches.size() ||
      state.patches.owning_rank[patch_index] != world_rank) {
    return false;
  }
  const std::uint32_t first = state.patches.first_cell[patch_index];
  const std::uint32_t count = state.patches.cell_count[patch_index];
  if (count == 0U || first >= owned_leaf_mask.size() ||
      static_cast<std::size_t>(first) + static_cast<std::size_t>(count) >
          owned_leaf_mask.size()) {
    return false;
  }
  return std::any_of(
      owned_leaf_mask.begin() + static_cast<std::ptrdiff_t>(first),
      owned_leaf_mask.begin() + static_cast<std::ptrdiff_t>(first + count),
      [](std::uint8_t value) { return value != 0U; });
}

}  // namespace

std::uint64_t metalDiffusionTopologyConstructionBoundBytes(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t world_rank,
    std::size_t remote_cell_count) {
  // Exact counts from patch dimensions (owned patches with owned leaves):
  //  I = interior faces, B = local patch-boundary cell faces,
  //  R = remote boundary cell faces (<= 6 per remote cell).
  // Cross/periodic/remote faces <= B + R (each fine face matches at most one
  // coarse face under 2:1 balance and each match consumes a distinct pair).
  constexpr std::uint64_t k_max = std::numeric_limits<std::uint64_t>::max();
  const auto add = [](std::uint64_t a, std::uint64_t b) { return a > k_max - b ? k_max : a + b; };
  const auto mul = [](std::uint64_t a, std::uint64_t b) { return (a != 0U && b > k_max / a) ? k_max : a * b; };
  std::uint64_t interior = 0U;
  std::uint64_t boundary = 0U;
  for (std::uint32_t p = 0; p < state.patches.size(); ++p) {
    if (!patchHasOwnedLeafCells(state, owned_leaf_mask, p, world_rank)) continue;
    const std::uint64_t nx = state.patches.cell_dim_x[p];
    const std::uint64_t ny = state.patches.cell_dim_y[p];
    const std::uint64_t nz = state.patches.cell_dim_z[p];
    if (nx == 0U || ny == 0U || nz == 0U) continue;
    interior = add(interior, add(add(mul(nx - 1U, mul(ny, nz)), mul(nx, mul(ny - 1U, nz))), mul(nx, mul(ny, nz - 1U))));
    boundary = add(boundary, mul(2U, add(add(mul(ny, nz), mul(nx, nz)), mul(nx, ny))));
  }
  const std::uint64_t remote_boundary = mul(6U, static_cast<std::uint64_t>(remote_cell_count));
  const std::uint64_t matched = add(boundary, remote_boundary);
  std::uint64_t bytes = 0U;
  bytes = add(bytes, mul(add(interior, matched), sizeof(physics::MetalDiffusionFace)));        // faces
  bytes = add(bytes, mul(mul(2U, matched), sizeof(RemoteDiffusionFace)));                       // remote faces (+growth)
  bytes = add(bytes, mul(state.cells.size(), sizeof(double)));                                  // strain
  bytes = add(bytes, mul(matched, sizeof(PatchBoundaryFace)));                                  // boundary list
  bytes = add(bytes, mul(matched, k_map_node_bytes + 2U * sizeof(std::size_t)));                // plane map
  bytes = add(bytes, mul(matched, k_set_node_bytes));                                           // unique sets
  bytes = add(bytes, mul(add(mul(2U, add(interior, matched)), matched), sizeof(DirectedFaceRef)));  // refs
  // std::vector geometric growth can transiently hold old+new buffers of
  // the push_back-built lists; allow one doubling of the largest list.
  bytes = add(bytes, mul(add(interior, matched), sizeof(physics::MetalDiffusionFace)));
  return bytes;
}

MetalDiffusionTopologyResult buildMetalDiffusionTopology(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t world_rank,
    double scale_factor,
    core::BoundaryCondition hydro_boundary,
    std::span<const RemoteDiffusionCell> remote_cells,
    const std::optional<std::array<double, 6>>& global_envelope_physical) {
  if (!std::isfinite(scale_factor) || !(scale_factor > 0.0)) {
    throw std::invalid_argument(
        "metal diffusion topology requires finite positive scale factor");
  }
  if (owned_leaf_mask.size() != state.cells.size()) {
    throw std::invalid_argument(
        "metal diffusion topology cell/mask extent mismatch");
  }

  MetalDiffusionTopologyResult result;
  result.strain_magnitude_code.assign(state.cells.size(), 0.0);
  std::vector<PatchBoundaryFace> patch_boundary_faces;
  {
    // Exact capacity (matches metalDiffusionTopologyConstructionBoundBytes):
    // no geometric-growth slack in the boundary list.
    std::size_t boundary_capacity = 6U * remote_cells.size();
    for (std::uint32_t p = 0; p < state.patches.size(); ++p) {
      if (!patchHasOwnedLeafCells(state, owned_leaf_mask, p, world_rank)) continue;
      const std::size_t nx = state.patches.cell_dim_x[p];
      const std::size_t ny = state.patches.cell_dim_y[p];
      const std::size_t nz = state.patches.cell_dim_z[p];
      boundary_capacity += 2U * (ny * nz + nx * nz + nx * ny);
    }
    patch_boundary_faces.reserve(boundary_capacity);
  }
  double topology_scale = 1.0;
  std::array<double, 3> domain_min{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()};
  std::array<double, 3> domain_max{
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};

  for (std::uint32_t patch_index = 0;
       patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.owning_rank[patch_index] != world_rank) {
      continue;
    }
    const bool has_owned_leaf_cells = patchHasOwnedLeafCells(
        state, owned_leaf_mask, patch_index, world_rank);
    const std::uint32_t first = state.patches.first_cell[patch_index];
    const std::uint32_t nx = state.patches.cell_dim_x[patch_index];
    const std::uint32_t ny = state.patches.cell_dim_y[patch_index];
    const std::uint32_t nz = state.patches.cell_dim_z[patch_index];
    if (nx == 0U || ny == 0U || nz == 0U) {
      throw std::runtime_error(
          "owned leaf diffusion patch has zero Cartesian dimension");
    }
    const std::uint64_t expected_count =
        static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny) *
        static_cast<std::uint64_t>(nz);
    if (expected_count != state.patches.cell_count[patch_index]) {
      throw std::runtime_error(
          "owned leaf diffusion patch cell_count disagrees with Cartesian dimensions");
    }
    if (static_cast<std::uint64_t>(first) + expected_count > state.cells.size()) {
      throw std::runtime_error(
          "owned leaf diffusion patch cell range exceeds dense gas storage");
    }

    const double dx = state.patches.extent_x_comoving[patch_index] /
        static_cast<double>(nx) * scale_factor;
    const double dy = state.patches.extent_y_comoving[patch_index] /
        static_cast<double>(ny) * scale_factor;
    const double dz = state.patches.extent_z_comoving[patch_index] /
        static_cast<double>(nz) * scale_factor;
    const double ox = state.patches.origin_x_comoving[patch_index] * scale_factor;
    const double oy = state.patches.origin_y_comoving[patch_index] * scale_factor;
    const double oz = state.patches.origin_z_comoving[patch_index] * scale_factor;
    if (!(dx > 0.0) || !(dy > 0.0) || !(dz > 0.0) ||
        !std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz) ||
        !std::isfinite(ox) || !std::isfinite(oy) || !std::isfinite(oz)) {
      throw std::runtime_error(
          "owned leaf diffusion patch has non-finite or non-positive geometry");
    }

    domain_min[0] = std::min(domain_min[0], ox);
    domain_min[1] = std::min(domain_min[1], oy);
    domain_min[2] = std::min(domain_min[2], oz);
    domain_max[0] = std::max(
        domain_max[0], ox + dx * static_cast<double>(nx));
    domain_max[1] = std::max(
        domain_max[1], oy + dy * static_cast<double>(ny));
    domain_max[2] = std::max(
        domain_max[2], oz + dz * static_cast<double>(nz));
    topology_scale = std::max({
        topology_scale, std::abs(ox), std::abs(oy), std::abs(oz),
        std::abs(dx * static_cast<double>(nx)),
        std::abs(dy * static_cast<double>(ny)),
        std::abs(dz * static_cast<double>(nz))});

    // Covered/coarse patches still define the physical domain envelope used
    // for periodic wrap detection, but only authoritative leaf cells may
    // contribute diffusion faces.
    if (!has_owned_leaf_cells) {
      continue;
    }

    const auto row = [first, nx, ny](
                         std::uint32_t i, std::uint32_t j,
                         std::uint32_t k) -> std::uint32_t {
      return first + i + nx * (j + ny * k);
    };
    const auto addTopologyFace = [&](std::uint32_t left,
                                     std::uint32_t right,
                                     std::uint8_t axis,
                                     double area,
                                     double distance,
                                     physics::MetalDiffusionBoundaryKind kind) {
      if (left >= owned_leaf_mask.size() || right >= owned_leaf_mask.size() ||
          owned_leaf_mask[left] == 0U || owned_leaf_mask[right] == 0U) {
        return;
      }
      const physics::MetalDiffusionFace face{
          .left_cell = left,
          .right_cell = right,
          .area_code = area,
          .center_distance_code = distance,
          .boundary_kind = kind,
          .axis = axis};
      result.faces.push_back(face);
    };
    const auto addBoundaryFace = [&](std::uint32_t cell,
                                     std::uint8_t axis,
                                     std::int8_t side,
                                     double plane,
                                     double u0,
                                     double u1,
                                     double v0,
                                     double v1,
                                     double normal_width) {
      if (cell >= owned_leaf_mask.size() || owned_leaf_mask[cell] == 0U) {
        return;
      }
      patch_boundary_faces.push_back(PatchBoundaryFace{
          .cell = cell,
          .axis = axis,
          .side = side,
          .plane = plane,
          .transverse_u_min = u0,
          .transverse_u_max = u1,
          .transverse_v_min = v0,
          .transverse_v_max = v1,
          .normal_width = normal_width});
    };

    for (std::uint32_t k = 0; k < nz; ++k) {
      for (std::uint32_t j = 0; j < ny; ++j) {
        for (std::uint32_t i = 0; i < nx; ++i) {
          const std::uint32_t left = row(i, j, k);
          if (i + 1U < nx) {
            addTopologyFace(
                left, row(i + 1U, j, k), 0U, dy * dz, dx,
                physics::MetalDiffusionBoundaryKind::kInternal);
            ++result.interior_patch_face_count;
          }
          if (j + 1U < ny) {
            addTopologyFace(
                left, row(i, j + 1U, k), 1U, dx * dz, dy,
                physics::MetalDiffusionBoundaryKind::kInternal);
            ++result.interior_patch_face_count;
          }
          if (k + 1U < nz) {
            addTopologyFace(
                left, row(i, j, k + 1U), 2U, dx * dy, dz,
                physics::MetalDiffusionBoundaryKind::kInternal);
            ++result.interior_patch_face_count;
          }

          const double x0 = ox + static_cast<double>(i) * dx;
          const double x1 = x0 + dx;
          const double y0 = oy + static_cast<double>(j) * dy;
          const double y1 = y0 + dy;
          const double z0 = oz + static_cast<double>(k) * dz;
          const double z1 = z0 + dz;
          if (i == 0U) {
            addBoundaryFace(left, 0U, -1, x0, y0, y1, z0, z1, dx);
          }
          if (i + 1U == nx) {
            addBoundaryFace(left, 0U, +1, x1, y0, y1, z0, z1, dx);
          }
          if (j == 0U) {
            addBoundaryFace(left, 1U, -1, y0, x0, x1, z0, z1, dy);
          }
          if (j + 1U == ny) {
            addBoundaryFace(left, 1U, +1, y1, x0, x1, z0, z1, dy);
          }
          if (k == 0U) {
            addBoundaryFace(left, 2U, -1, z0, x0, x1, y0, y1, dz);
          }
          if (k + 1U == nz) {
            addBoundaryFace(left, 2U, +1, z1, x0, x1, y0, y1, dz);
          }
        }
      }
    }
  }

  // G2.1: remote owned leaf cells contribute their patch-boundary faces as
  // virtual rows [local_count, local_count + remote_count). Geometry uses the
  // owner's exact arithmetic so both ranks see identical planes.
  const std::size_t local_count = state.cells.size();
  if (remote_cells.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) - local_count) {
    throw std::overflow_error("metal diffusion topology remote cell count exceeds uint32 row space");
  }
  for (std::size_t r = 0; r < remote_cells.size(); ++r) {
    const RemoteDiffusionCell& cell = remote_cells[r];
    const auto virtual_row = static_cast<std::uint32_t>(local_count + r);
    std::array<double, 3> lo{}, hi{}, width{};
    for (int c = 0; c < 3; ++c) {
      if (cell.patch_dims[c] == 0U || cell.index[c] >= cell.patch_dims[c]) {
        throw std::runtime_error("remote diffusion cell has an invalid patch index");
      }
      width[c] = cell.cell_width_comoving[c] * scale_factor;
      const double origin = cell.patch_origin_comoving[c] * scale_factor;
      lo[c] = origin + static_cast<double>(cell.index[c]) * width[c];
      hi[c] = lo[c] + width[c];
      if (!(width[c] > 0.0) || !std::isfinite(lo[c]) || !std::isfinite(hi[c])) {
        throw std::runtime_error("remote diffusion cell has non-finite or non-positive geometry");
      }
    }
    const auto push = [&](std::uint8_t axis, std::int8_t side) {
      const int u = axis == 0U ? 1 : 0;
      const int v = axis == 2U ? 1 : 2;
      patch_boundary_faces.push_back(PatchBoundaryFace{
          .cell = virtual_row, .axis = axis, .side = side,
          .plane = side < 0 ? lo[axis] : hi[axis],
          .transverse_u_min = lo[u], .transverse_u_max = hi[u],
          .transverse_v_min = lo[v], .transverse_v_max = hi[v],
          .normal_width = width[axis]});
    };
    for (std::uint8_t axis = 0U; axis < 3U; ++axis) {
      if (cell.index[axis] == 0U) push(axis, -1);
      if (cell.index[axis] + 1U == cell.patch_dims[axis]) push(axis, +1);
    }
  }
  if (global_envelope_physical.has_value()) {
    for (int c = 0; c < 3; ++c) {
      domain_min[c] = (*global_envelope_physical)[c];
      domain_max[c] = (*global_envelope_physical)[3 + c];
      topology_scale = std::max({topology_scale, std::abs(domain_min[c]), std::abs(domain_max[c])});
    }
  }

  if (result.faces.empty() && patch_boundary_faces.empty()) {
    return result;
  }

  const double face_tolerance = 1.0e-10 * topology_scale;
  const double quantization = std::max(face_tolerance, 1.0e-14);
  std::map<std::pair<std::uint8_t, std::int64_t>, std::vector<std::size_t>>
      lower_faces_by_plane;
  for (std::size_t index = 0; index < patch_boundary_faces.size(); ++index) {
    const PatchBoundaryFace& candidate = patch_boundary_faces[index];
    if (candidate.side >= 0) {
      continue;
    }
    const auto bin = static_cast<std::int64_t>(
        std::llround(candidate.plane / quantization));
    lower_faces_by_plane[{candidate.axis, bin}].push_back(index);
  }

  std::set<std::tuple<std::uint32_t, std::uint32_t, std::uint8_t>>
      unique_cross_patch_faces;
  std::set<std::tuple<std::uint32_t, std::uint32_t, std::uint8_t>>
      unique_remote_faces;
  const auto addOverlappingInterface = [&](const PatchBoundaryFace& upper,
                                           const PatchBoundaryFace& lower,
                                           physics::MetalDiffusionBoundaryKind kind,
                                           std::uint64_t& counter) {
    if (upper.cell == lower.cell || upper.axis != lower.axis) {
      return;
    }
    const double overlap_u =
        std::min(upper.transverse_u_max, lower.transverse_u_max) -
        std::max(upper.transverse_u_min, lower.transverse_u_min);
    const double overlap_v =
        std::min(upper.transverse_v_max, lower.transverse_v_max) -
        std::max(upper.transverse_v_min, lower.transverse_v_min);
    if (!(overlap_u > face_tolerance) || !(overlap_v > face_tolerance)) {
      return;
    }
    const bool upper_remote = upper.cell >= local_count;
    const bool lower_remote = lower.cell >= local_count;
    if (upper_remote && lower_remote) {
      return;  // remote-remote interfaces belong to their owners
    }
    if (upper_remote || lower_remote) {
      // Local "upper" face (side +1) touches the remote cell's lower face:
      // the remote cell lies on the +axis side, and vice versa.
      const PatchBoundaryFace& local_face = upper_remote ? lower : upper;
      const PatchBoundaryFace& remote_face = upper_remote ? upper : lower;
      const std::uint32_t remote_index = remote_face.cell - static_cast<std::uint32_t>(local_count);
      const std::int8_t local_side = local_face.side;
      // A direct and a periodic adjacency between the same pair are distinct
      // faces (e.g. two ranks splitting a periodic axis), so the key carries
      // the side and the boundary kind.
      if (!unique_remote_faces.emplace(local_face.cell, remote_index,
                                       static_cast<std::uint8_t>(upper.axis * 4U + (local_side > 0 ? 2U : 0U) +
                                           (kind == physics::MetalDiffusionBoundaryKind::kPeriodic ? 1U : 0U)))
               .second) {
        return;
      }
      const RemoteDiffusionCell& remote = remote_cells[remote_index];
      result.remote_faces.push_back(RemoteDiffusionFace{
          .local_row = local_face.cell,
          .remote_index = remote_index,
          .local_gas_cell_id = state.gas_cells.gas_cell_id[local_face.cell],
          .remote_gas_cell_id = remote.gas_cell_id,
          .remote_rank = remote.owner_rank,
          .axis = upper.axis,
          .local_side = local_side,
          .periodic = kind == physics::MetalDiffusionBoundaryKind::kPeriodic,
          .area_code = overlap_u * overlap_v,
          .center_distance_code = 0.5 * (upper.normal_width + lower.normal_width)});
      ++result.remote_face_count;
      return;
    }
    const std::uint32_t lo = std::min(upper.cell, lower.cell);
    const std::uint32_t hi = std::max(upper.cell, lower.cell);
    if (!unique_cross_patch_faces.emplace(lo, hi, upper.axis).second) {
      return;
    }
    const physics::MetalDiffusionFace face{
        .left_cell = upper.cell,
        .right_cell = lower.cell,
        .area_code = overlap_u * overlap_v,
        .center_distance_code = 0.5 *
            (upper.normal_width + lower.normal_width),
        .boundary_kind = kind,
        .axis = upper.axis};
    result.faces.push_back(face);
    ++counter;
  };

  for (const PatchBoundaryFace& upper : patch_boundary_faces) {
    if (upper.side <= 0) {
      continue;
    }
    const auto upper_bin = static_cast<std::int64_t>(
        std::llround(upper.plane / quantization));
    for (std::int64_t delta_bin = -1; delta_bin <= 1; ++delta_bin) {
      const auto found = lower_faces_by_plane.find(
          {upper.axis, upper_bin + delta_bin});
      if (found == lower_faces_by_plane.end()) {
        continue;
      }
      for (const std::size_t lower_index : found->second) {
        const PatchBoundaryFace& lower = patch_boundary_faces[lower_index];
        if (std::abs(upper.plane - lower.plane) > face_tolerance) {
          continue;
        }
        addOverlappingInterface(
            upper, lower, physics::MetalDiffusionBoundaryKind::kInternal,
            result.cross_patch_face_count);
      }
    }
  }

  if (hydro_boundary == core::BoundaryCondition::kPeriodic) {
    for (std::uint8_t axis = 0U; axis < 3U; ++axis) {
      if (!std::isfinite(domain_min[axis]) ||
          !std::isfinite(domain_max[axis]) ||
          !(domain_max[axis] > domain_min[axis])) {
        throw std::runtime_error(
            "periodic metal diffusion requires finite positive domain extent");
      }
      for (const PatchBoundaryFace& upper : patch_boundary_faces) {
        if (upper.axis != axis || upper.side <= 0 ||
            std::abs(upper.plane - domain_max[axis]) > face_tolerance) {
          continue;
        }
        for (const PatchBoundaryFace& lower : patch_boundary_faces) {
          if (lower.axis != axis || lower.side >= 0 ||
              std::abs(lower.plane - domain_min[axis]) > face_tolerance) {
            continue;
          }
          addOverlappingInterface(
              upper, lower, physics::MetalDiffusionBoundaryKind::kPeriodic,
              result.periodic_wrap_face_count);
        }
      }
    }
  }

  // Build a compact directed face-reference list instead of retaining six
  // NeighborAccumulator objects (~240 B) per gas cell. The refs are sorted by
  // cell and face index so each cell reconstructs the exact same face-order
  // weighted derivative with one stack-local accumulator.
  // Remote faces are referenced from their local side only, at face index
  // faces.size() + remote face index.
  if (result.faces.size() + result.remote_faces.size() > static_cast<std::size_t>(
          std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error(
        "metal diffusion topology face count exceeds uint32 reference capacity");
  }
  std::vector<DirectedFaceRef> directed_face_refs;
  if (result.faces.size() >
      std::numeric_limits<std::size_t>::max() / 2U) {
    throw std::overflow_error(
        "metal diffusion topology directed reference count overflow");
  }
  directed_face_refs.reserve(result.faces.size() * 2U + result.remote_faces.size());
  for (std::size_t remote_index = 0U; remote_index < result.remote_faces.size(); ++remote_index) {
    directed_face_refs.push_back(DirectedFaceRef{
        .cell = result.remote_faces[remote_index].local_row,
        .face_index = static_cast<std::uint32_t>(result.faces.size() + remote_index)});
  }
  for (std::size_t face_index = 0U; face_index < result.faces.size(); ++face_index) {
    const auto face_row = static_cast<std::uint32_t>(face_index);
    const physics::MetalDiffusionFace& face = result.faces[face_index];
    directed_face_refs.push_back(DirectedFaceRef{.cell = face.left_cell,
                                                 .face_index = face_row});
    directed_face_refs.push_back(DirectedFaceRef{.cell = face.right_cell,
                                                 .face_index = face_row});
  }
  std::sort(
      directed_face_refs.begin(), directed_face_refs.end(),
      [](const DirectedFaceRef& lhs, const DirectedFaceRef& rhs) {
        if (lhs.cell != rhs.cell) {
          return lhs.cell < rhs.cell;
        }
        return lhs.face_index < rhs.face_index;
      });

  const auto checkedCapacityBytes = [](std::size_t capacity, std::size_t width) {
    if (capacity != 0U &&
        width > std::numeric_limits<std::uint64_t>::max() / capacity) {
      throw std::overflow_error(
          "metal diffusion topology scratch byte count overflow");
    }
    return static_cast<std::uint64_t>(capacity) *
        static_cast<std::uint64_t>(width);
  };
  const std::uint64_t patch_boundary_scratch_bytes = checkedCapacityBytes(
      patch_boundary_faces.capacity(), sizeof(PatchBoundaryFace));
  const std::uint64_t directed_face_ref_scratch_bytes = checkedCapacityBytes(
      directed_face_refs.capacity(), sizeof(DirectedFaceRef));
  if (directed_face_ref_scratch_bytes >
      std::numeric_limits<std::uint64_t>::max() - patch_boundary_scratch_bytes) {
    throw std::overflow_error(
        "metal diffusion topology aggregate scratch byte count overflow");
  }
  // G0.5: node-based containers are included with the documented per-node
  // allowance (k_map_node_bytes / k_set_node_bytes), matching the admission
  // bound metalDiffusionTopologyConstructionBoundBytes.
  std::uint64_t lower_index_capacity = 0U;
  for (const auto& entry : lower_faces_by_plane) lower_index_capacity += entry.second.capacity();
  const std::uint64_t map_bytes =
      checkedCapacityBytes(lower_faces_by_plane.size(), k_map_node_bytes) +
      checkedCapacityBytes(static_cast<std::size_t>(lower_index_capacity), sizeof(std::size_t));
  const std::uint64_t set_bytes = checkedCapacityBytes(
      unique_cross_patch_faces.size() + unique_remote_faces.size(), k_set_node_bytes);
  result.construction_scratch_high_water_bytes =
      patch_boundary_scratch_bytes + directed_face_ref_scratch_bytes + map_bytes + set_bytes;

  const std::array<std::span<const double>, 3> velocity_fields{
      state.gas_cells.velocity_x_peculiar,
      state.gas_cells.velocity_y_peculiar,
      state.gas_cells.velocity_z_peculiar};
  std::size_t ref_begin = 0U;
  while (ref_begin < directed_face_refs.size()) {
    const std::uint32_t cell_index = directed_face_refs[ref_begin].cell;
    std::size_t ref_end = ref_begin + 1U;
    while (ref_end < directed_face_refs.size() &&
           directed_face_refs[ref_end].cell == cell_index) {
      ++ref_end;
    }
    if (cell_index >= owned_leaf_mask.size() ||
        owned_leaf_mask[cell_index] == 0U) {
      ref_begin = ref_end;
      continue;
    }

    std::array<std::array<NeighborAccumulator, 2>, 3> accumulators{};
    for (std::size_t ref_index = ref_begin; ref_index < ref_end; ++ref_index) {
      const std::uint32_t face_ref = directed_face_refs[ref_index].face_index;
      if (face_ref >= result.faces.size()) {
        // Remote neighbour: its step-start velocity from the exchange, so the
        // strain (and hence kappa) does not depend on the decomposition.
        const RemoteDiffusionFace& rf = result.remote_faces[face_ref - result.faces.size()];
        if (rf.axis >= 3U) {
          throw std::runtime_error("metal diffusion topology remote face has invalid axis");
        }
        NeighborAccumulator& accumulator = accumulators[rf.axis][rf.local_side > 0 ? 1U : 0U];
        const RemoteDiffusionCell& remote = remote_cells[rf.remote_index];
        for (std::size_t component = 0U; component < 3U; ++component) {
          accumulator.velocity_area_sum[component] += rf.area_code * remote.velocity_peculiar[component];
        }
        accumulator.distance_area_sum += rf.area_code * rf.center_distance_code;
        accumulator.area_sum += rf.area_code;
        continue;
      }
      const physics::MetalDiffusionFace& face = result.faces[face_ref];
      if (face.axis >= 3U) {
        throw std::runtime_error("metal diffusion topology face has invalid axis");
      }
      const bool is_left = cell_index == face.left_cell;
      const std::uint32_t neighbor = is_left ? face.right_cell : face.left_cell;
      const std::size_t direction = is_left ? 1U : 0U;
      NeighborAccumulator& accumulator = accumulators[face.axis][direction];
      for (std::size_t component = 0U; component < 3U; ++component) {
        accumulator.velocity_area_sum[component] +=
            face.area_code * velocity_fields[component][neighbor];
      }
      accumulator.distance_area_sum +=
          face.area_code * face.center_distance_code;
      accumulator.area_sum += face.area_code;
    }

    physics::MetalDiffusionVelocityGradient gradient{};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const NeighborAccumulator& lower = accumulators[axis][0];
      const NeighborAccumulator& upper = accumulators[axis][1];
      for (std::size_t component = 0U; component < 3U; ++component) {
        const double self_velocity = velocity_fields[component][cell_index];
        double derivative = 0.0;
        if (lower.area_sum > 0.0 && upper.area_sum > 0.0) {
          const double lower_velocity =
              lower.velocity_area_sum[component] / lower.area_sum;
          const double upper_velocity =
              upper.velocity_area_sum[component] / upper.area_sum;
          const double lower_distance =
              lower.distance_area_sum / lower.area_sum;
          const double upper_distance =
              upper.distance_area_sum / upper.area_sum;
          derivative = (upper_velocity - lower_velocity) /
              std::max(lower_distance + upper_distance, 1.0e-30);
        } else if (upper.area_sum > 0.0) {
          derivative =
              (upper.velocity_area_sum[component] / upper.area_sum -
               self_velocity) /
              std::max(upper.distance_area_sum / upper.area_sum, 1.0e-30);
        } else if (lower.area_sum > 0.0) {
          derivative =
              (self_velocity -
               lower.velocity_area_sum[component] / lower.area_sum) /
              std::max(lower.distance_area_sum / lower.area_sum, 1.0e-30);
        }
        gradient.grad[component][axis] = derivative;
      }
    }
    result.strain_magnitude_code[cell_index] =
        physics::traceFreeStrainMagnitude(gradient);
    ref_begin = ref_end;
  }

  return result;
}

}  // namespace cosmosim::workflows::internal
