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
// refusal, allocation failure or narrowing failure is turned into a
// collective error BEFORE any rank enters the next collective (patch 52:
// this now also covers the control-plane arrays allocated before the count
// exchange; see exchangeDestinationRecordsInto). All ranks always participate (empty streams allowed); receive
// order is deterministic (ascending source rank, then send order).
// Without MPI (or one rank) records addressed to the local rank are returned
// and any other destination is refused.

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
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

// Caller-owned bytes of one side of a chunked lane exchange for `entries`
// lane vectors of `lane_count` lanes: the outbound chunk records and their
// destination ranks, plus the reassembled flat lanes/ids/scalars and the
// reassembly sort order. Saturating.
//
// Patch 52 (P1-51-02): the exchange primitive itself now admits and owns
// every allocation it makes (control plane, bucket/cursor, round staging and
// the final receive buffer, see exchangeDestinationRecordsInto), so this
// caller-side figure deliberately EXCLUDES them; adding them here would
// double-charge the same bytes against the node budget.
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

// Received payload whose memory admission lives exactly as long as the
// payload (patch 51, M01/D01; lifetime order fixed in patch 52, P1-51-02):
// `lease` (kCommunication) covers the owned capacity of `records`/`bytes`.
// `lease` is declared FIRST so that implicit destruction (reverse
// declaration order) frees the storage BEFORE the admission is released;
// move assignment likewise frees the old storage before releasing the old
// lease. releaseStorage() frees early in the same order.
struct GovernedRecordBuffer {
  core::MemoryReservation lease;
  std::vector<std::byte> bytes;

  GovernedRecordBuffer() = default;
  GovernedRecordBuffer(const GovernedRecordBuffer&) = delete;
  GovernedRecordBuffer& operator=(const GovernedRecordBuffer&) = delete;
  GovernedRecordBuffer(GovernedRecordBuffer&&) noexcept = default;
  GovernedRecordBuffer& operator=(GovernedRecordBuffer&& other) noexcept {
    if (this != &other) {
      bytes = std::move(other.bytes);
      lease = std::move(other.lease);
    }
    return *this;
  }
  ~GovernedRecordBuffer() = default;
  void releaseStorage() noexcept {
    std::vector<std::byte>().swap(bytes);
    lease.release();
  }
};
template <class Record>
struct GovernedRecords {
  core::MemoryReservation lease;
  std::vector<Record> records;

  GovernedRecords() = default;
  GovernedRecords(const GovernedRecords&) = delete;
  GovernedRecords& operator=(const GovernedRecords&) = delete;
  GovernedRecords(GovernedRecords&&) noexcept = default;
  GovernedRecords& operator=(GovernedRecords&& other) noexcept {
    if (this != &other) {
      records = std::move(other.records);
      lease = std::move(other.lease);
    }
    return *this;
  }
  ~GovernedRecords() = default;
  void releaseStorage() noexcept {
    std::vector<Record>().swap(records);
    lease.release();
  }
};

// Receive-buffer allocator of the governed primitive: called at most once,
// inside the guarded admission phase (step 4 below, i.e. BEFORE the
// collective readiness vote), with the total received byte count; returns
// the base of a caller-owned buffer of at least that many bytes. It may
// throw: the failure is made collective before any payload round.
using DestinationReceiveAllocator = std::byte* (*)(void* owner, std::size_t bytes);

// Collective-safe destination record exchange (patch 52, P1-51-01). Every
// rank calls it in the same order, possibly with an empty stream. Protocol:
//   1. [local, guarded] argument/destination validation, control-plane
//      admission (count/offset/displacement/cursor arrays, kCommunication)
//      and allocation, per-destination tally;
//   2. [collective] failure vote "destination record exchange preparation";
//   3. [collective] MPI_Alltoall of per-peer record counts (no allocation
//      or throwing operation between 2 and 3);
//   4. [local, noexcept arithmetic] round plan; [collective] Allreduce(max)
//      of the round count;
//   5. [local, guarded] MPI int-range checks of every round count and
//      displacement, admission of the final receive buffer (`lease`,
//      returned to the caller) and of the scratch (bucket index + maximal
//      round staging), allocation of the receive buffer (via `allocate`),
//      bucket and staging capacity;
//   6. [collective] failure vote "destination record exchange admission";
//   7. [collective x rounds] MPI_Alltoallv payload rounds. Between 6 and the
//      last round nothing allocates: staging is resized within reserved
//      capacity, counts are already range-checked.
// An MPI error return inside steps 3/4/7 throws on that rank and remains
// communicator-fatal (default MPI error handler); ordinary local failures
// (bad_alloc, admission refusal, narrowing) are always reported through the
// votes in 2/6 so no peer is stranded inside a collective. On any throw the
// scratch/control reservations are released by RAII and the caller-owned
// receive buffer is destroyed by the caller's unwinding.
// Without MPI (or one rank) records addressed to the local rank are copied
// into the receive buffer (admitted when a governor is given) and any other
// destination is refused.
[[nodiscard]] core::MemoryReservation exchangeDestinationRecordsInto(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    core::MemoryGovernor* governor,
    std::string_view owner_label,
    DestinationReceiveAllocator allocate,
    void* allocate_owner);

// Typed governed front end: received records land directly in the returned
// typed vector (no byte buffer, no second typed copy) and the lease covering
// them is returned with them.
template <class Record>
[[nodiscard]] GovernedRecords<Record> exchangeGovernedRecords(
    const MpiContext& mpi_context,
    std::span<const Record> records,
    std::span<const int> destinations,
    core::MemoryGovernor* governor,
    std::string_view owner_label) {
  static_assert(std::is_trivially_copyable_v<Record>);
  static_assert(std::is_default_constructible_v<Record>);
  GovernedRecords<Record> out;
  const std::span<const std::byte> bytes(
      reinterpret_cast<const std::byte*>(records.data()), records.size() * sizeof(Record));
  out.lease = exchangeDestinationRecordsInto(
      mpi_context, bytes, sizeof(Record), destinations, governor, owner_label,
      [](void* owner, std::size_t total) -> std::byte* {
        auto& target = *static_cast<std::vector<Record>*>(owner);
        if (total % sizeof(Record) != 0U) {
          throw std::runtime_error("destination record exchange: received byte count is not a whole record stream");
        }
        target.resize(total / sizeof(Record));
        return reinterpret_cast<std::byte*>(target.data());
      },
      &out.records);
  return out;
}

// Ungoverned byte-level primitive (same collective-safe protocol, no
// admission). Retained for API compatibility; production streams use the
// governed front ends.
[[nodiscard]] std::vector<std::byte> executeBlockingDestinationRecordExchange(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    std::uint64_t exchange_sequence);

// Governed byte-level primitive (protocol of exchangeDestinationRecordsInto).
[[nodiscard]] GovernedRecordBuffer executeGovernedDestinationRecordExchange(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor,
    std::string_view owner_label);

[[nodiscard]] GovernedRecords<CompositionDonorRequestRecord> exchangeCompositionDonorRequests(
    const MpiContext& mpi_context,
    std::span<const CompositionDonorRequestRecord> requests,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor = nullptr);

// Chunked lane exchange (collective). Validates every outbound record
// (schema, source == this rank, destination matches) before the collective.
// Records are received directly into the returned typed buffer.
[[nodiscard]] GovernedRecords<CompositionLaneChunkRecord> exchangeCompositionLaneChunks(
    const MpiContext& mpi_context,
    std::span<const CompositionLaneChunkRecord> records,
    std::span<const int> destinations,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor = nullptr,
    std::string_view owner_label = "parallel.composition_lane_exchange");

}  // namespace cosmosim::parallel
