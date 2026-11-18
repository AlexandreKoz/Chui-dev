#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "parallel/internal/distributed_memory_support.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
PmMeshOwnershipDescriptor PmSlabLayout::ownershipDescriptor(
    std::uint64_t decomposition_epoch,
    std::string decomposition_mode) const {
  PmMeshOwnershipDescriptor descriptor{
      .decomposition_mode = std::move(decomposition_mode),
      .owner_rank = world_rank,
      .decomposition_epoch = decomposition_epoch,
      .global_nx = global_nx,
      .global_ny = global_ny,
      .global_nz = global_nz,
      .begin_x = owned_x.begin_x,
      .end_x = owned_x.end_x,
  };
  validatePmMeshOwnershipDescriptor(descriptor);
  return descriptor;
}

void validatePmMeshOwnershipDescriptor(const PmMeshOwnershipDescriptor& descriptor) {
  if (descriptor.decomposition_mode != "slab" && descriptor.decomposition_mode != "pencil") {
    throw std::invalid_argument("PM mesh ownership descriptor has unsupported decomposition mode");
  }
  if (descriptor.owner_rank < 0) {
    throw std::invalid_argument("PM mesh owner_rank must be non-negative");
  }
  if (descriptor.global_nx == 0 || descriptor.global_ny == 0 || descriptor.global_nz == 0) {
    throw std::invalid_argument("PM mesh descriptor global dimensions must be positive");
  }
  if (descriptor.begin_x > descriptor.end_x || descriptor.end_x > descriptor.global_nx) {
    throw std::invalid_argument("PM mesh x ownership range is invalid");
  }
}

void validateTreePseudoParticleDescriptor(const TreePseudoParticleDescriptor& descriptor) {
  if (descriptor.wire_version != 1U) {
    throw std::invalid_argument("tree pseudo-particle wire version is unsupported");
  }
  if (descriptor.source_rank < 0) {
    throw std::invalid_argument("tree pseudo-particle source_rank must be non-negative");
  }
}

void validateTreePseudoParticlePacket(const TreePseudoParticlePacket& packet) {
  validateTreePseudoParticleDescriptor(packet.descriptor);
  if (packet.source_count == 0 && packet.mass_code != 0.0) {
    throw std::invalid_argument("empty tree pseudo-particle packet cannot carry non-zero mass");
  }
  if (packet.mass_code < 0.0 || !std::isfinite(packet.mass_code)) {
    throw std::invalid_argument("tree pseudo-particle mass must be finite and non-negative");
  }
  const std::array values{packet.center_x_comoving, packet.center_y_comoving, packet.center_z_comoving,
                          packet.min_x_comoving, packet.max_x_comoving, packet.min_y_comoving,
                          packet.max_y_comoving, packet.min_z_comoving, packet.max_z_comoving};
  for (const double value : values) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("tree pseudo-particle packet contains non-finite geometry");
    }
  }
  if (packet.min_x_comoving > packet.max_x_comoving || packet.min_y_comoving > packet.max_y_comoving ||
      packet.min_z_comoving > packet.max_z_comoving) {
    throw std::invalid_argument("tree pseudo-particle packet bounds are invalid");
  }
  if (packet.child_count > 8U) {
    throw std::invalid_argument("tree pseudo-particle packet child count is invalid");
  }
  if (packet.is_leaf > 1U) {
    throw std::invalid_argument("tree pseudo-particle packet leaf flag is invalid");
  }
  if (packet.geometry_frame > 1U) {
    throw std::invalid_argument("tree pseudo-particle packet geometry frame is invalid");
  }
}

void validateHydroGhostCellDescriptor(const HydroGhostCellDescriptor& descriptor) {
  if (descriptor.gas_cell_id == 0) {
    throw std::invalid_argument("hydro ghost cell descriptor requires a non-zero gas_cell_id");
  }
  if (descriptor.owner_rank < 0 || descriptor.consumer_rank < 0) {
    throw std::invalid_argument("hydro ghost cell ranks must be non-negative");
  }
  if (descriptor.owner_rank == descriptor.consumer_rank) {
    throw std::invalid_argument("hydro ghost cell must be consumed on a non-owner rank");
  }
  if (!descriptor.boundary_state_only) {
    throw std::invalid_argument("hydro ghost cells are boundary exchange state, not authoritative conserved truth");
  }
}

void validateHydroGhostCellRequest(const HydroGhostCellRequest& request) {
  validateHydroGhostCellDescriptor(request.descriptor);
  if (request.face_key == 0) {
    throw std::invalid_argument("hydro ghost cell request requires a non-zero face key");
  }
  if (request.axis > 2U || request.side > 1U) {
    throw std::invalid_argument("hydro ghost cell request carries invalid face orientation metadata");
  }
}

void validateHydroGhostCellPayloadRecord(const HydroGhostCellPayloadRecord& record) {
  validateHydroGhostCellDescriptor(record.descriptor);
  if (record.face_key == 0) {
    throw std::invalid_argument("hydro ghost cell payload requires a non-zero face key");
  }
  if (!std::isfinite(record.mass_density_comoving) ||
      !std::isfinite(record.momentum_density_x_comoving) ||
      !std::isfinite(record.momentum_density_y_comoving) ||
      !std::isfinite(record.momentum_density_z_comoving) ||
      !std::isfinite(record.total_energy_density_comoving) ||
      !std::isfinite(record.metal_mass_density_comoving)) {
    throw std::invalid_argument("hydro ghost cell payload contains non-finite conserved state");
  }
  if (record.mass_density_comoving <= 0.0) {
    throw std::invalid_argument("hydro ghost cell payload requires positive mass density");
  }
  if (record.metal_mass_density_comoving < 0.0 ||
      record.metal_mass_density_comoving > record.mass_density_comoving * (1.0 + 1.0e-12)) {
    throw std::invalid_argument("hydro ghost cell payload contains invalid metal mass density");
  }
}

void validateAmrPatchExchangeDescriptor(const AmrPatchExchangeDescriptor& descriptor) {
  if (descriptor.owner_rank < 0 || descriptor.peer_rank < 0) {
    throw std::invalid_argument("AMR patch exchange ranks must be non-negative");
  }
  if (!descriptor.metadata_only && descriptor.owner_rank != descriptor.peer_rank) {
    throw std::invalid_argument("remote AMR patch exchange cannot mutate authoritative patch metadata");
  }
}

void validateAmrPatchPayloadRecord(const AmrPatchPayloadRecord& record) {
  if (record.owner_rank < 0) {
    throw std::invalid_argument("AMR patch payload owner_rank must be non-negative");
  }
  if (record.patch_id == 0) {
    throw std::invalid_argument("AMR patch payload patch_id must be non-zero");
  }
  if (record.cell_count == 0) {
    throw std::invalid_argument("AMR patch payload cannot describe an empty patch");
  }
  if (record.extent_x_comoving <= 0.0 || record.extent_y_comoving <= 0.0 || record.extent_z_comoving <= 0.0 ||
      record.cell_dim_x == 0U || record.cell_dim_y == 0U || record.cell_dim_z == 0U) {
    throw std::invalid_argument("AMR patch payload requires explicit positive patch geometry");
  }
  const std::uint64_t geometry_cells =
      static_cast<std::uint64_t>(record.cell_dim_x) *
      static_cast<std::uint64_t>(record.cell_dim_y) *
      static_cast<std::uint64_t>(record.cell_dim_z);
  if (geometry_cells != record.cell_count) {
    throw std::invalid_argument("AMR patch payload cell_count does not match explicit patch geometry");
  }
  if (!std::isfinite(record.origin_x_comoving) || !std::isfinite(record.origin_y_comoving) ||
      !std::isfinite(record.origin_z_comoving) || !std::isfinite(record.extent_x_comoving) ||
      !std::isfinite(record.extent_y_comoving) || !std::isfinite(record.extent_z_comoving)) {
    throw std::invalid_argument("AMR patch payload contains non-finite patch geometry");
  }
  if (!std::isfinite(record.cell_mass_sum_code) || !std::isfinite(record.gas_internal_energy_sum_code)) {
    throw std::invalid_argument("AMR patch payload contains non-finite cell sums");
  }
}

void validateAmrPatchCellPayloadRecord(const AmrPatchCellPayloadRecord& record) {
  if (record.owner_rank < 0) {
    throw std::invalid_argument("AMR patch cell payload owner_rank must be non-negative");
  }
  if (record.patch_id == 0) {
    throw std::invalid_argument("AMR patch cell payload patch_id must be non-zero");
  }
  if (record.gas_cell_id == 0) {
    throw std::invalid_argument("AMR patch cell payload must carry stable gas-cell identity");
  }
  if (!std::isfinite(record.center_x_comoving) || !std::isfinite(record.center_y_comoving) ||
      !std::isfinite(record.center_z_comoving) || !std::isfinite(record.mass_code) ||
      !std::isfinite(record.velocity_x_peculiar) || !std::isfinite(record.velocity_y_peculiar) ||
      !std::isfinite(record.velocity_z_peculiar) || !std::isfinite(record.density_code) ||
      !std::isfinite(record.pressure_code) ||
      !std::isfinite(record.internal_energy_code) || !std::isfinite(record.temperature_code) ||
      !std::isfinite(record.sound_speed_code) || !std::isfinite(record.metal_mass_code)) {
    throw std::invalid_argument("AMR patch cell payload contains non-finite state");
  }
  if (record.density_code <= 0.0 || record.pressure_code <= 0.0) {
    throw std::invalid_argument("AMR patch cell payload requires positive thermodynamic state");
  }
  if (record.metal_mass_code < 0.0 ||
      record.metal_mass_code > record.mass_code + 1.0e-12 * std::max(1.0, std::abs(record.mass_code))) {
    throw std::invalid_argument("AMR patch cell payload metal mass must lie within the gas-cell mass");
  }
  if (record.ism_regime_tag > 1U && record.ism_regime_tag != k_amr_payload_ism_regime_tag_unknown) {
    throw std::invalid_argument("AMR patch cell payload carries an invalid ISM regime tag");
  }
}

void validateAmrFluxRegisterPayloadRecord(const AmrFluxRegisterPayloadRecord& record) {
  if (record.register_key == 0U || record.coarse_patch_id == 0U || record.coarse_gas_cell_id == 0U) {
    throw std::invalid_argument("AMR flux-register payload requires non-zero stable identity fields");
  }
  if (record.source_rank < 0 || record.owner_rank < 0) {
    throw std::invalid_argument("AMR flux-register payload ranks must be non-negative");
  }
  if (record.axis > 2U || record.orientation > 1U) {
    throw std::invalid_argument("AMR flux-register payload carries invalid face identity");
  }
  if (record.face_area_comov <= 0.0 || record.coarse_area_comov <= 0.0 ||
      record.fine_area_comov <= 0.0 || record.dt_code <= 0.0) {
    throw std::invalid_argument("AMR flux-register payload requires positive area and timestep metadata");
  }
  if (record.coarse_face_count == 0U && record.fine_face_count == 0U) {
    throw std::invalid_argument("AMR flux-register payload must carry at least one face contribution");
  }
  const std::array values{
      record.coarse_mass_flux_code,
      record.coarse_momentum_x_flux_code,
      record.coarse_momentum_y_flux_code,
      record.coarse_momentum_z_flux_code,
      record.coarse_total_energy_flux_code,
      record.coarse_metal_mass_flux_code,
      record.fine_mass_flux_code,
      record.fine_momentum_x_flux_code,
      record.fine_momentum_y_flux_code,
      record.fine_momentum_z_flux_code,
      record.fine_total_energy_flux_code,
      record.fine_metal_mass_flux_code,
      record.face_area_comov,
      record.coarse_area_comov,
      record.fine_area_comov,
      record.dt_code};
  for (const double value : values) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("AMR flux-register payload contains non-finite values");
    }
  }
}

void validateHydroConservativeFluxCorrectionRecord(const HydroConservativeFluxCorrectionRecord& record) {
  if (record.gas_cell_id == 0) {
    throw std::invalid_argument("hydro conservative flux correction requires a non-zero gas_cell_id");
  }
  if (record.source_rank < 0 || record.owner_rank < 0) {
    throw std::invalid_argument("hydro conservative flux correction ranks must be non-negative");
  }
  if (!std::isfinite(record.delta_mass_density_comoving) ||
      !std::isfinite(record.delta_momentum_density_x_comoving) ||
      !std::isfinite(record.delta_momentum_density_y_comoving) ||
      !std::isfinite(record.delta_momentum_density_z_comoving) ||
      !std::isfinite(record.delta_total_energy_density_comoving) ||
      !std::isfinite(record.delta_metal_mass_density_comoving)) {
    throw std::invalid_argument("hydro conservative flux correction contains non-finite state");
  }
}

void recordDistributedProfiling(
    core::ProfilerSession* profiler,
    const LoadBalanceMetrics& metrics,
    std::uint64_t ghost_exchange_send_bytes,
    std::uint64_t ghost_exchange_recv_bytes) {
  if (profiler == nullptr) {
    return;
  }

  profiler->counters().setCount("parallel.ghost_exchange_send_bytes", ghost_exchange_send_bytes);
  profiler->counters().setCount("parallel.ghost_exchange_recv_bytes", ghost_exchange_recv_bytes);
  profiler->counters().setCount("parallel.total_memory_bytes", metrics.total_memory_bytes);

  const std::uint64_t imbalance_ppm = (metrics.weighted_imbalance_ratio <= 0.0)
                                          ? 0ULL
                                          : static_cast<std::uint64_t>(std::llround(metrics.weighted_imbalance_ratio * 1.0e6));
  profiler->counters().setCount("parallel.weighted_imbalance_ratio_ppm", imbalance_ppm);
  auto record_component_total = [&](std::string_view name, const std::vector<double>& values) {
    const double total = std::accumulate(values.begin(), values.end(), 0.0);
    profiler->counters().setCount(std::string(name), static_cast<std::uint64_t>(std::llround(std::max(0.0, total))));
  };
  record_component_total("parallel.weight_component_particle_count", metrics.particle_count_cost_by_rank);
  record_component_total("parallel.weight_component_gas_cell", metrics.gas_cell_cost_by_rank);
  record_component_total("parallel.weight_component_tree_interaction", metrics.tree_interaction_cost_by_rank);
  record_component_total("parallel.weight_component_pm_mesh", metrics.pm_mesh_cost_by_rank);
  record_component_total("parallel.weight_component_amr_patch", metrics.amr_patch_cost_by_rank);
  record_component_total("parallel.weight_component_active_fraction", metrics.active_fraction_cost_by_rank);
  record_component_total("parallel.weight_component_memory_pressure", metrics.memory_pressure_cost_by_rank);
  record_component_total("parallel.weight_component_transient_memory", metrics.transient_memory_cost_by_rank);
  record_component_total("parallel.weight_component_source_event", metrics.source_event_cost_by_rank);
  record_component_total("parallel.weight_component_communication", metrics.communication_cost_by_rank);
  record_component_total("parallel.weight_component_gpu_occupancy", metrics.gpu_occupancy_cost_by_rank);

  if (!metrics.weighted_load_by_rank.empty()) {
    const auto max_rank_it = std::max_element(metrics.weighted_load_by_rank.begin(), metrics.weighted_load_by_rank.end());
    const std::size_t max_rank = static_cast<std::size_t>(std::distance(metrics.weighted_load_by_rank.begin(), max_rank_it));
    profiler->counters().setCount("parallel.max_weighted_load_rank", static_cast<std::uint64_t>(max_rank));
    profiler->counters().setCount(
        "parallel.max_weighted_load",
        static_cast<std::uint64_t>(std::llround(std::max(0.0, *max_rank_it))));

    struct ComponentView {
      std::string_view name;
      const std::vector<double>* values;
    };
    const std::array<ComponentView, 12> components{{
        {"particle_count", &metrics.particle_count_cost_by_rank},
        {"gas_cell", &metrics.gas_cell_cost_by_rank},
        {"tree_interaction", &metrics.tree_interaction_cost_by_rank},
        {"pm_mesh", &metrics.pm_mesh_cost_by_rank},
        {"amr_patch", &metrics.amr_patch_cost_by_rank},
        {"active_fraction", &metrics.active_fraction_cost_by_rank},
        {"memory_pressure", &metrics.memory_pressure_cost_by_rank},
        {"transient_memory", &metrics.transient_memory_cost_by_rank},
        {"source_event", &metrics.source_event_cost_by_rank},
        {"communication", &metrics.communication_cost_by_rank},
        {"gpu_occupancy", &metrics.gpu_occupancy_cost_by_rank},
        {"generic_work", &metrics.generic_work_cost_by_rank},
    }};
    std::string_view dominant_component = "none";
    double dominant_value = 0.0;
    for (const ComponentView& component : components) {
      if (component.values->size() <= max_rank) {
        continue;
      }
      const double value = (*component.values)[max_rank];
      if (value > dominant_value) {
        dominant_value = value;
        dominant_component = component.name;
      }
    }
    profiler->counters().setCount(
        "parallel.max_rank_dominant_component_cost",
        static_cast<std::uint64_t>(std::llround(std::max(0.0, dominant_value))));
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.hotspot",
        .severity = metrics.weighted_imbalance_ratio > 1.25 ? core::RuntimeEventSeverity::kWarning
                                                            : core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.domain_decomposition",
        .step_index = std::nullopt,
        .simulation_time_code = std::nullopt,
        .scale_factor = std::nullopt,
        .message = "domain decomposition load hotspot attribution",
        .payload = {{"max_rank", std::to_string(max_rank)},
                    {"max_rank_load", std::to_string(*max_rank_it)},
                    {"mean_load", std::to_string(metrics.mean_weighted_load)},
                    {"weighted_imbalance_ratio", std::to_string(metrics.weighted_imbalance_ratio)},
                    {"memory_imbalance_ratio", std::to_string(metrics.memory_imbalance_ratio)},
                    {"dominant_component", std::string(dominant_component)},
                    {"dominant_component_cost", std::to_string(dominant_value)}}});
  }

  const std::uint64_t bytes_moved = ghost_exchange_send_bytes + ghost_exchange_recv_bytes;
  profiler->addBytesMoved(bytes_moved);
}


}  // namespace cosmosim::parallel
