#pragma once

// CHUI distributed composition transport payloads (wire schema 2, patch 50).
//
// Owner-targeted record streams carry the composition sidecar on the
// distributed production AMR path:
//
//  * donor fractions: for every remote ghost source cell that is a donor of a
//    composition face transfer, its owner sends the lane mass fractions
//    (lane mass / gas mass, gas mass = elements + dust) at the step-start
//    state; the receiver advects remote inflow with exactly the owner's
//    composition (request/response by stable gas_cell_id);
//  * corrections: a composition lane delta computed on a rank for a cell
//    owned by another rank (fine-side mirrored transfer onto a remote coarse
//    cell at a coarse/fine interface) is delivered to that owner and applied
//    inside its composition transaction;
//  * element-diffusion interface state and equal/opposite interface deltas
//    (workflows/internal/distributed_element_diffusion.cpp).
//
// P0-49-01 (schema 2): the former fixed 32-lane POD could not represent the
// selectable C02 ion-rich layouts (e.g. C+O = 37 lanes, C/O/Si/Fe + dust +
// SGS = 82 lanes). A layout of L lanes is now split into
// ceil(L / k_composition_lane_chunk_width) fixed-size chunk records. Each
// chunk carries the wire schema version, the stream kind, the authoritative
// layout signature (GasCompositionLayoutDescriptor::layout_signature), the
// total lane count, chunk index/count/lane range, the exchange sequence and
// the stable gas_cell_id plus source/destination ranks. The receiver
// reassembles by (gas_cell_id, source_rank) and refuses truncation,
// duplicated or missing chunks, a foreign layout, a stale sequence or a
// misaddressed record (exactly-once application per (cell, source)).
//
// Bounds: L <= k_max_composition_layout_lanes (static upper bound of every
// layout CompositionLayout::make can construct; checked there). The byte
// exchange runs in bounded rounds of at most
// k_max_destination_exchange_round_bytes_per_peer bytes per peer, so no MPI
// int count/displacement overflows; the receive capacity of every round is
// admitted through the optional MemoryGovernor before allocation, and a local
// refusal (or any local narrowing failure) is turned into a collective error
// BEFORE any rank enters the payload collective (no deadlock, no partial
// application). All ranks always participate (empty streams allowed); receive
// order is deterministic (ascending source rank, then send order).
// Without MPI (or one rank) records addressed to the local rank are returned
// and any other destination is refused.

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::parallel {

inline constexpr std::uint32_t k_composition_wire_schema_version = 2U;
// Payload doubles per chunk record (fixed so every record is a POD).
inline constexpr std::size_t k_composition_lane_chunk_width = 32U;
// Static upper bound of every constructible composition layout:
// 11 elements + 10 detailed light species + 58 C/O/Si/Fe ion stages + 2 dust
// + 1 SGS = 82 <= 128. physics::CompositionLayout::make refuses a layout
// above this bound (no silent truncation can occur at the wire).
inline constexpr std::size_t k_max_composition_layout_lanes = 128U;
inline constexpr std::size_t k_max_composition_lane_chunks =
    (k_max_composition_layout_lanes + k_composition_lane_chunk_width - 1U) / k_composition_lane_chunk_width;
// Per-peer payload cap of one Alltoallv round (keeps every MPI int count and
// displacement far below INT_MAX for any rank count used in production).
inline constexpr std::uint64_t k_max_destination_exchange_round_bytes_per_peer = 64ULL * 1024ULL * 1024ULL;

enum class CompositionLaneStream : std::uint32_t {
  kDonorFraction = 1U,
  kCorrection = 2U,
  kDiffusionState = 3U,
  kDiffusionDelta = 4U,
};

struct CompositionDonorRequestRecord {
  std::uint64_t gas_cell_id = 0U;
  int requester_rank = 0;
  int owner_rank = 0;
};

// One chunk of one cell's composition lane vector.
struct CompositionLaneChunkRecord {
  std::uint64_t gas_cell_id = 0U;
  std::uint64_t layout_signature = 0U;
  std::uint64_t exchange_sequence = 0U;
  std::int32_t source_rank = 0;
  std::int32_t destination_rank = 0;
  std::uint32_t schema_version = 0U;
  std::uint32_t stream_kind = 0U;
  std::uint32_t lane_count = 0U;
  std::uint16_t chunk_index = 0U;
  std::uint16_t chunk_count = 0U;
  std::uint32_t chunk_lane_begin = 0U;
  std::uint32_t chunk_lane_count = 0U;
  // Stream-defined per-cell scalar (diffusion: composition gas mass); the
  // same value is repeated in every chunk and checked on reassembly.
  double scalar = 0.0;
  std::array<double, k_composition_lane_chunk_width> lanes{};
};

// Identity every chunk of one exchange must carry.
struct CompositionLaneStreamIdentity {
  CompositionLaneStream stream = CompositionLaneStream::kDonorFraction;
  std::uint64_t layout_signature = 0U;
  std::uint32_t lane_count = 0U;
  std::uint64_t exchange_sequence = 0U;
};

// Reassembled per-(cell, source) lane vectors, sorted by (gas_cell_id,
// source_rank); lanes are stored flat (entry i at [i*lane_count, ...)).
struct AssembledCompositionLanes {
  std::size_t lane_count = 0U;
  std::vector<std::uint64_t> gas_cell_id;
  std::vector<int> source_rank;
  std::vector<double> scalar;
  std::vector<double> lanes;

  [[nodiscard]] std::size_t size() const noexcept { return gas_cell_id.size(); }
  [[nodiscard]] std::span<const double> lanesOf(std::size_t entry) const {
    return std::span<const double>(lanes.data() + entry * lane_count, lane_count);
  }
};

// ceil(lane_count / chunk width); throws std::length_error above the static
// layout bound. 0 lanes -> 0 chunks.
[[nodiscard]] std::size_t compositionLaneChunkCount(std::size_t lane_count);

// Bytes that coexist for `entries` lane vectors of `lane_count` lanes on one
// side of a chunked exchange (send records + bucketed send copy + receive
// records + reassembled flat lanes). Saturating.
[[nodiscard]] std::uint64_t compositionLaneExchangeAdmissionBytes(std::uint64_t entries, std::size_t lane_count);

// Splits one cell's lane vector into chunk records addressed to
// destination_rank (records and destinations grow in lockstep).
void appendCompositionLaneChunks(
    const CompositionLaneStreamIdentity& identity,
    std::uint64_t gas_cell_id,
    int source_rank,
    int destination_rank,
    double scalar,
    std::span<const double> lanes,
    std::vector<CompositionLaneChunkRecord>& records,
    std::vector<int>& destinations);

// Validates and reassembles received chunks (see header comment). Throws on
// any schema/identity/addressing/coverage violation.
[[nodiscard]] AssembledCompositionLanes assembleCompositionLaneChunks(
    std::span<const CompositionLaneChunkRecord> records,
    const CompositionLaneStreamIdentity& expected,
    int receiver_rank);

// Byte-level primitive: record i (record_bytes wide) goes to destination[i].
// Collective. When governor != nullptr every round's receive capacity is
// reserved (kCommunication, owner label) before allocation; refusal is
// collective before the payload Alltoallv.
[[nodiscard]] std::vector<std::byte> executeBlockingDestinationRecordExchange(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor = nullptr,
    std::string_view owner_label = "parallel.destination_record_exchange");

[[nodiscard]] std::vector<CompositionDonorRequestRecord> exchangeCompositionDonorRequests(
    const MpiContext& mpi_context,
    std::span<const CompositionDonorRequestRecord> requests,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor = nullptr);

// Chunked lane exchange (collective). Validates every outbound record
// (schema, source == this rank, destination matches) before the collective.
[[nodiscard]] std::vector<CompositionLaneChunkRecord> exchangeCompositionLaneChunks(
    const MpiContext& mpi_context,
    std::span<const CompositionLaneChunkRecord> records,
    std::span<const int> destinations,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor = nullptr,
    std::string_view owner_label = "parallel.composition_lane_exchange");

}  // namespace cosmosim::parallel
