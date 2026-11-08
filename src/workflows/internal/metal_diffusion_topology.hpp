#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/ism_dust_sgs.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"

namespace cosmosim::workflows::internal {

// G2.1: a remote owned leaf cell adjacent to a local patch (one ghost layer
// delivered by the directed AMR exchange). Geometry is reproduced with the
// owner's exact arithmetic (patch origin + index * width, then * a) so both
// ranks compute bitwise-identical face planes.
struct RemoteDiffusionCell {
  std::uint64_t gas_cell_id = 0U;
  int owner_rank = 0;
  std::array<double, 3> patch_origin_comoving{};
  std::array<double, 3> cell_width_comoving{};
  std::array<std::uint32_t, 3> index{};
  std::array<std::uint32_t, 3> patch_dims{};
  std::array<double, 3> velocity_peculiar{};
  double density_physical_code = 0.0;
  // T02: owner's AMR level of this leaf (face clock = max of both sides).
  std::int32_t level = 0;
};

// Interface between a local owned leaf row and a remote owned leaf cell.
// local_side = +1 when the remote cell lies on the +axis side of the local
// cell. The unique flux authority is the side with the smaller gas_cell_id.
using RemoteDiffusionFace = physics::SgsRemoteShearFace;

struct MetalDiffusionTopologyResult {
  std::vector<physics::MetalDiffusionFace> faces;
  std::vector<RemoteDiffusionFace> remote_faces;
  std::vector<double> strain_magnitude_code;
  std::uint64_t construction_scratch_high_water_bytes = 0U;
  std::uint64_t interior_patch_face_count = 0U;
  std::uint64_t cross_patch_face_count = 0U;
  std::uint64_t periodic_wrap_face_count = 0U;
  std::uint64_t remote_face_count = 0U;
};

// G0.5: exact-count admission bound for buildMetalDiffusionTopology's
// simultaneously live allocations (result faces/remote faces/strain plus the
// boundary-face list, plane map, unique-face set and directed references),
// computed from patch dimensions BEFORE construction. Node-based containers
// use the documented per-node allowances in the implementation.
[[nodiscard]] std::uint64_t metalDiffusionTopologyConstructionBoundBytes(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t world_rank,
    std::size_t remote_cell_count);

// Construct the authoritative local diffusion interface graph from AMR patch
// geometry. Faces are oriented from the negative-axis cell to the positive-axis
// cell, including periodic wraps. The same graph is used to reconstruct the
// trace-free strain magnitude that sets the Smagorinsky diffusivity. Gradient
// accumulation is cell-local scratch; no 3x3 gradient is retained per gas cell.
[[nodiscard]] MetalDiffusionTopologyResult buildMetalDiffusionTopology(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t world_rank,
    double scale_factor,
    core::BoundaryCondition hydro_boundary,
    std::span<const RemoteDiffusionCell> remote_cells = {},
    // Global physical domain envelope {min_x,min_y,min_z,max_x,max_y,max_z};
    // REQUIRED on more than one rank for periodic wraps (the local patch
    // envelope is not the box). Absent: the local envelope (single rank).
    const std::optional<std::array<double, 6>>& global_envelope_physical = std::nullopt);

}  // namespace cosmosim::workflows::internal
