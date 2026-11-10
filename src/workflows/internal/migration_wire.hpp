#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal::migration_wire {

inline constexpr std::uint32_t k_particle_record_wire_version = 3U;  // v3: CHUI Part 2 gas composition payload
inline constexpr std::uint32_t k_dmo_particle_record_wire_version = 1U;
// v4: CHUI Part 2 gas composition payload. v5 (P1-55-01): temporal-history
// cells carry the dual-energy auxiliary rho*e start/end endpoints.
inline constexpr std::uint32_t k_amr_patch_record_wire_version = 5U;
inline constexpr std::uint32_t k_fragment_wire_version = 1U;
inline constexpr std::size_t k_fragment_header_bytes =
    sizeof(std::uint32_t) + sizeof(std::uint64_t) * 3U + sizeof(std::uint32_t);


struct DmoParticleMigrationRecord {
  std::uint64_t particle_id = 0U;
  std::uint64_t sfc_key = 0U;
  double position_x_comoving = 0.0;
  double position_y_comoving = 0.0;
  double position_z_comoving = 0.0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double mass_code = 0.0;
};

struct FragmentView {
  std::uint64_t record_sequence = 0U;
  std::uint64_t record_total_bytes = 0U;
  std::uint64_t fragment_offset = 0U;
  std::span<const std::uint8_t> payload{};
};

struct ParticleWireReferenceWidths {
  std::size_t dark_matter_bytes = 0U;
  std::size_t gas_extra_bytes = 0U;
  std::size_t star_extra_bytes = 0U;
  std::size_t black_hole_extra_bytes = 0U;
  std::size_t tracer_extra_bytes = 0U;
};

struct PacketCapacityPlan {
  std::size_t per_peer_packet_bytes = 0U;
  std::size_t fragment_payload_bytes = 0U;
  std::size_t aggregate_packet_bytes = 0U;
};

[[nodiscard]] PacketCapacityPlan planPacketCapacity(
    std::size_t transport_round_limit_bytes,
    std::size_t rank_count);

[[nodiscard]] std::vector<std::uint8_t> encodeDmoParticleMigrationRecord(
    const DmoParticleMigrationRecord& record);
[[nodiscard]] DmoParticleMigrationRecord decodeDmoParticleMigrationRecord(
    std::span<const std::uint8_t> bytes);
[[nodiscard]] constexpr std::size_t dmoParticleMigrationWireBytes() noexcept {
  return sizeof(std::uint32_t) + sizeof(std::uint64_t) * 2U + sizeof(double) * 7U;
}

[[nodiscard]] std::vector<std::uint8_t> encodeParticleMigrationRecord(
    const core::ParticleMigrationRecord& record);
[[nodiscard]] core::ParticleMigrationRecord decodeParticleMigrationRecord(
    std::span<const std::uint8_t> bytes);

[[nodiscard]] std::vector<std::uint8_t> encodeAmrPatchMigrationRecord(
    const core::AmrPatchMigrationRecord& record);
[[nodiscard]] core::AmrPatchMigrationRecord decodeAmrPatchMigrationRecord(
    std::span<const std::uint8_t> bytes);

[[nodiscard]] std::vector<std::uint8_t> encodeFragment(
    std::uint64_t record_sequence,
    std::uint64_t record_total_bytes,
    std::uint64_t fragment_offset,
    std::span<const std::uint8_t> payload);
[[nodiscard]] FragmentView decodeFragment(std::span<const std::uint8_t> bytes);

// Source-grounded upper bounds that require no population-scale temporary
// encoding. Runtime migration admission uses these before record/wire staging.
[[nodiscard]] std::size_t estimateParticleMigrationWireUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_particle_index);
[[nodiscard]] std::size_t estimateAmrPatchMigrationWireUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_patch_index);
[[nodiscard]] std::size_t estimateParticleMigrationDynamicHeapUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_particle_index);
[[nodiscard]] std::size_t estimateAmrPatchMigrationDynamicHeapUpperBoundBytes(
    const core::SimulationState& state,
    std::uint32_t local_patch_index);

[[nodiscard]] ParticleWireReferenceWidths particleWireReferenceWidths();

}  // namespace cosmosim::workflows::internal::migration_wire
