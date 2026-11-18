#include "workflows/internal/amr_migration_payload.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace cosmosim::workflows::internal {
namespace {

[[nodiscard]] parallel::AmrPatchCellPayloadRecord makePatchCellPayloadRecord(
    const core::SimulationState& state,
    std::size_t patch_index,
    std::uint32_t offset,
    int world_rank) {
  const std::uint32_t first_cell = state.patches.first_cell[patch_index];
  const std::uint32_t cell_index = first_cell + offset;
  parallel::AmrPatchCellPayloadRecord record;
  record.patch_id = state.patches.patch_id[patch_index];
  record.owner_rank = world_rank;
  record.local_cell_offset = offset;
  record.patch_index = static_cast<std::uint32_t>(patch_index);
  record.center_x_comoving = state.cells.center_x_comoving[cell_index];
  record.center_y_comoving = state.cells.center_y_comoving[cell_index];
  record.center_z_comoving = state.cells.center_z_comoving[cell_index];
  record.mass_code = state.cells.mass_code[cell_index];
  record.time_bin = state.cells.time_bin[cell_index];
  const core::GasCellIdentityRecord& identity = gasCellIdentityRecordForLocalRow(
      state, cell_index, "makePatchCellPayloadRecord");
  record.gas_cell_id = identity.gas_cell_id;
  record.parent_particle_id = identity.parent_particle_id.value_or(0U);
  record.velocity_x_peculiar =
      cell_index < state.gas_cells.velocity_x_peculiar.size()
          ? state.gas_cells.velocity_x_peculiar[cell_index]
          : 0.0;
  record.velocity_y_peculiar =
      cell_index < state.gas_cells.velocity_y_peculiar.size()
          ? state.gas_cells.velocity_y_peculiar[cell_index]
          : 0.0;
  record.velocity_z_peculiar =
      cell_index < state.gas_cells.velocity_z_peculiar.size()
          ? state.gas_cells.velocity_z_peculiar[cell_index]
          : 0.0;
  record.density_code = cell_index < state.gas_cells.density_code.size()
      ? state.gas_cells.density_code[cell_index]
      : 0.0;
  record.pressure_code = cell_index < state.gas_cells.pressure_code.size()
      ? state.gas_cells.pressure_code[cell_index]
      : 0.0;
  record.internal_energy_code =
      cell_index < state.gas_cells.internal_energy_code.size()
          ? state.gas_cells.internal_energy_code[cell_index]
          : 0.0;
  record.temperature_code = cell_index < state.gas_cells.temperature_code.size()
      ? state.gas_cells.temperature_code[cell_index]
      : 0.0;
  record.sound_speed_code = cell_index < state.gas_cells.sound_speed_code.size()
      ? state.gas_cells.sound_speed_code[cell_index]
      : 0.0;
  record.metal_mass_code = cell_index < state.gas_cells.metal_mass_code.size()
      ? state.gas_cells.metal_mass_code[cell_index]
      : 0.0;
  // I01: ship the owner's regime authority with the hydro state.
  record.ism_regime_tag = state.gas_composition.enabled() &&
          cell_index < state.gas_composition.ism_regime_tag.size()
      ? state.gas_composition.ism_regime_tag[cell_index]
      : parallel::k_amr_payload_ism_regime_tag_unknown;
  parallel::validateAmrPatchCellPayloadRecord(record);
  return record;
}

[[nodiscard]] std::uint8_t faceBit(parallel::AmrPatchBoundaryFace face) noexcept {
  return static_cast<std::uint8_t>(face);
}

[[nodiscard]] std::size_t faceIndex(parallel::AmrPatchBoundaryFace face) noexcept {
  switch (face) {
    case parallel::AmrPatchBoundaryFace::kXLower: return 0U;
    case parallel::AmrPatchBoundaryFace::kXUpper: return 1U;
    case parallel::AmrPatchBoundaryFace::kYLower: return 2U;
    case parallel::AmrPatchBoundaryFace::kYUpper: return 3U;
    case parallel::AmrPatchBoundaryFace::kZLower: return 4U;
    case parallel::AmrPatchBoundaryFace::kZUpper: return 5U;
  }
  return 0U;
}

[[nodiscard]] std::uint16_t faceDepth(
    const parallel::AmrPatchBoundaryCellRequest& request,
    parallel::AmrPatchBoundaryFace face) {
  const bool selected = (request.boundary_face_mask & faceBit(face)) != 0U;
  const std::uint16_t depth = request.boundary_face_depths[faceIndex(face)];
  if (selected && depth == 0U) {
    throw std::invalid_argument("AMR boundary-cell request selected a face with zero depth");
  }
  if (!selected && depth != 0U) {
    throw std::invalid_argument("AMR boundary-cell request carries depth for an unselected face");
  }
  return selected ? depth : 0U;
}

void mergeRequest(
    parallel::AmrPatchBoundaryCellRequest& destination,
    const parallel::AmrPatchBoundaryCellRequest& source) {
  if (destination.patch_id == 0U) {
    destination.patch_id = source.patch_id;
  }
  if (destination.patch_id != source.patch_id) {
    throw std::invalid_argument("AMR boundary-cell request merge crossed patch identities");
  }
  destination.boundary_face_mask |= source.boundary_face_mask;
  for (std::size_t face = 0; face < destination.boundary_face_depths.size(); ++face) {
    destination.boundary_face_depths[face] = std::max(
        destination.boundary_face_depths[face], source.boundary_face_depths[face]);
  }
}

}  // namespace

std::vector<parallel::AmrPatchPayloadRecord>
buildMigrationAmrPatchPayloadRecords(
    const core::SimulationState& state,
    int world_rank) {
  std::vector<parallel::AmrPatchPayloadRecord> records;
  if (world_rank < 0 || state.patches.size() == 0) {
    return records;
  }
  records.reserve(state.patches.size());
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.owning_rank[patch_index] != static_cast<std::uint32_t>(world_rank)) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    if (cell_count == 0U) {
      continue;
    }
    if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) >
        state.cells.size()) {
      throw std::runtime_error("AMR patch payload build found a patch range outside CellSoa");
    }
    parallel::AmrPatchPayloadRecord record;
    record.patch_id = state.patches.patch_id[patch_index];
    record.parent_patch_id = state.patches.parent_patch_id[patch_index];
    record.morton_key = state.patches.morton_key[patch_index];
    record.owner_rank = world_rank;
    record.level = static_cast<std::uint32_t>(
        std::max<std::int32_t>(state.patches.level[patch_index], 0));
    record.first_cell = first_cell;
    record.cell_count = cell_count;
    record.origin_x_comoving = state.patches.origin_x_comoving[patch_index];
    record.origin_y_comoving = state.patches.origin_y_comoving[patch_index];
    record.origin_z_comoving = state.patches.origin_z_comoving[patch_index];
    record.extent_x_comoving = state.patches.extent_x_comoving[patch_index];
    record.extent_y_comoving = state.patches.extent_y_comoving[patch_index];
    record.extent_z_comoving = state.patches.extent_z_comoving[patch_index];
    record.cell_dim_x = state.patches.cell_dim_x[patch_index];
    record.cell_dim_y = state.patches.cell_dim_y[patch_index];
    record.cell_dim_z = state.patches.cell_dim_z[patch_index];
    record.decomposition_epoch = state.gasCellIdentityGeneration();
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::uint32_t cell_index = first_cell + offset;
      record.cell_mass_sum_code += state.cells.mass_code[cell_index];
      if (cell_index < state.gas_cells.internal_energy_code.size()) {
        record.gas_internal_energy_sum_code += state.gas_cells.internal_energy_code[cell_index];
      }
    }
    parallel::validateAmrPatchPayloadRecord(record);
    records.push_back(record);
  }
  return records;
}

void fillMigrationAmrPatchBoundaryCellPayloadChunk(
    const core::SimulationState& state,
    int world_rank,
    std::span<const parallel::AmrPatchBoundaryCellRequest> requests,
    std::uint64_t first_record,
    std::size_t max_records,
    std::vector<parallel::AmrPatchCellPayloadRecord>& output) {
  output.clear();
  if (world_rank < 0 || requests.empty() || max_records == 0U) {
    return;
  }
  output.reserve(max_records);

  std::unordered_map<std::uint64_t, parallel::AmrPatchBoundaryCellRequest> request_by_patch_id;
  request_by_patch_id.reserve(requests.size());
  for (const parallel::AmrPatchBoundaryCellRequest& request : requests) {
    if (request.patch_id == 0U || request.boundary_face_mask == 0U) {
      throw std::invalid_argument("AMR boundary-cell request must name a patch and at least one face");
    }
    auto [it, inserted] = request_by_patch_id.emplace(request.patch_id, request);
    if (!inserted) {
      mergeRequest(it->second, request);
    }
  }

  std::uint64_t logical_index = 0U;
  const auto emit = [&](std::size_t patch_index, std::uint32_t offset) -> bool {
    if (logical_index++ < first_record) {
      return false;
    }
    output.push_back(makePatchCellPayloadRecord(
        state, patch_index, offset, world_rank));
    return output.size() == max_records;
  };

  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    const std::uint64_t patch_id = state.patches.patch_id[patch_index];
    const auto request_it = request_by_patch_id.find(patch_id);
    if (request_it == request_by_patch_id.end()) {
      continue;
    }
    const auto& request = request_it->second;
    if (state.patches.owning_rank[patch_index] != static_cast<std::uint32_t>(world_rank)) {
      throw std::invalid_argument("AMR boundary-cell request selected a patch not owned by the local rank");
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) >
        state.cells.size()) {
      throw std::runtime_error("AMR boundary-cell payload build found a patch range outside CellSoa");
    }
    const std::uint16_t nx = state.patches.cell_dim_x[patch_index];
    const std::uint16_t ny = state.patches.cell_dim_y[patch_index];
    const std::uint16_t nz = state.patches.cell_dim_z[patch_index];
    const std::uint64_t geometry_cell_count = static_cast<std::uint64_t>(nx) * ny * nz;
    if (geometry_cell_count != cell_count) {
      throw std::runtime_error("AMR boundary-cell payload build found cell_dims/cell_count mismatch");
    }

    const std::size_t x_lower = faceDepth(request, parallel::AmrPatchBoundaryFace::kXLower);
    const std::size_t x_upper = faceDepth(request, parallel::AmrPatchBoundaryFace::kXUpper);
    const std::size_t y_lower = faceDepth(request, parallel::AmrPatchBoundaryFace::kYLower);
    const std::size_t y_upper = faceDepth(request, parallel::AmrPatchBoundaryFace::kYUpper);
    const std::size_t z_lower = faceDepth(request, parallel::AmrPatchBoundaryFace::kZLower);
    const std::size_t z_upper = faceDepth(request, parallel::AmrPatchBoundaryFace::kZUpper);
    if (x_lower > nx || x_upper > nx || y_lower > ny || y_upper > ny ||
        z_lower > nz || z_upper > nz) {
      throw std::invalid_argument("AMR boundary-cell request depth exceeds patch dimension");
    }

    for (std::uint32_t k = 0; k < nz; ++k) {
      const bool selected_z =
          (z_lower != 0U && k < z_lower) ||
          (z_upper != 0U && static_cast<std::size_t>(k) >= static_cast<std::size_t>(nz) - z_upper);
      for (std::uint32_t j = 0; j < ny; ++j) {
        const bool selected_y =
            (y_lower != 0U && j < y_lower) ||
            (y_upper != 0U && static_cast<std::size_t>(j) >= static_cast<std::size_t>(ny) - y_upper);
        const std::uint32_t row_base = nx * (j + static_cast<std::uint32_t>(ny) * k);
        if (selected_z || selected_y) {
          for (std::uint32_t i = 0; i < nx; ++i) {
            if (emit(patch_index, row_base + i)) {
              return;
            }
          }
          continue;
        }

        const std::size_t lower_count = std::min<std::size_t>(x_lower, nx);
        for (std::size_t i = 0; i < lower_count; ++i) {
          if (emit(patch_index, row_base + static_cast<std::uint32_t>(i))) {
            return;
          }
        }
        const std::size_t upper_begin = std::max<std::size_t>(
            lower_count, static_cast<std::size_t>(nx) - std::min<std::size_t>(x_upper, nx));
        for (std::size_t i = upper_begin; i < nx; ++i) {
          if (emit(patch_index, row_base + static_cast<std::uint32_t>(i))) {
            return;
          }
        }
      }
    }
    request_by_patch_id.erase(request_it);
  }
  if (!request_by_patch_id.empty()) {
    throw std::invalid_argument("AMR boundary-cell request referenced an unknown local patch_id");
  }
}

}  // namespace cosmosim::workflows::internal
