#include "cosmosim/parallel/composition_exchange.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
namespace {

static_assert(std::is_trivially_copyable_v<CompositionLaneChunkRecord>);
static_assert(std::is_trivially_copyable_v<CompositionDonorRequestRecord>);
static_assert(k_max_composition_lane_chunks <= std::numeric_limits<std::uint16_t>::max());

[[nodiscard]] std::uint64_t saturatingMul(std::uint64_t a, std::uint64_t b) noexcept {
  return (a != 0U && b > std::numeric_limits<std::uint64_t>::max() / a)
      ? std::numeric_limits<std::uint64_t>::max() : a * b;
}

[[nodiscard]] std::uint64_t saturatingAdd(std::uint64_t a, std::uint64_t b) noexcept {
  return a > std::numeric_limits<std::uint64_t>::max() - b ? std::numeric_limits<std::uint64_t>::max() : a + b;
}

// Destination allocator: called once (inside the admission phase, before
// the collective failure vote) with the total received byte count; returns
// the base of a buffer of at least that many bytes that outlives the call.
using ReceiveAllocator = std::byte* (*)(void* owner, std::size_t bytes);

// Core collective implementation shared by the byte and typed front ends.
// Returns the lease covering the final receive buffer (or an invalid
// reservation without a governor).
core::MemoryReservation exchangeIntoDestination(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    core::MemoryGovernor* governor,
    std::string_view owner_label,
    ReceiveAllocator allocate,
    void* allocate_owner);

template <class Record>
[[nodiscard]] GovernedRecords<Record> exchangeGovernedTyped(
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
  out.lease = exchangeIntoDestination(
      mpi_context, bytes, sizeof(Record), destinations, governor, owner_label,
      [](void* owner, std::size_t total) -> std::byte* {
        auto& target = *static_cast<std::vector<Record>*>(owner);
        if (total % sizeof(Record) != 0U) {
          throw std::runtime_error("composition exchange: received byte count is not a whole record stream");
        }
        target.resize(total / sizeof(Record));
        return reinterpret_cast<std::byte*>(target.data());
      },
      &out.records);
  return out;
}

}  // namespace

std::size_t compositionLaneChunkCount(std::size_t lane_count) {
  if (lane_count > k_max_composition_layout_lanes) {
    throw std::length_error(
        "composition lane count " + std::to_string(lane_count) +
        " exceeds the static composition layout bound " + std::to_string(k_max_composition_layout_lanes));
  }
  return (lane_count + k_composition_lane_chunk_width - 1U) / k_composition_lane_chunk_width;
}

std::uint64_t compositionLaneExchangeAdmissionBytes(std::uint64_t entries, std::size_t lane_count) {
  const std::uint64_t chunks = saturatingMul(entries, compositionLaneChunkCount(lane_count));
  // outbound records + destination ranks + bucket index (sender side) and
  // received records + reassembled flat lanes/ids/scalars (receiver side).
  std::uint64_t bytes = saturatingMul(chunks, 2U * sizeof(CompositionLaneChunkRecord) + sizeof(int) +
                                                  sizeof(std::size_t));
  bytes = saturatingAdd(bytes, saturatingMul(entries, saturatingAdd(saturatingMul(lane_count, sizeof(double)),
                                                                    sizeof(std::uint64_t) + sizeof(int) +
                                                                        sizeof(double))));
  return bytes;
}

void appendCompositionLaneChunks(
    const CompositionLaneStreamIdentity& identity,
    std::uint64_t gas_cell_id,
    int source_rank,
    int destination_rank,
    double scalar,
    std::span<const double> lanes,
    std::vector<CompositionLaneChunkRecord>& records,
    std::vector<int>& destinations) {
  if (gas_cell_id == 0U) {
    throw std::invalid_argument("composition lane chunk: gas_cell_id 0 is not a stable cell identifier");
  }
  if (lanes.size() != identity.lane_count) {
    throw std::invalid_argument("composition lane chunk: lane vector does not match the stream lane count");
  }
  const std::size_t chunk_count = compositionLaneChunkCount(lanes.size());
  for (std::size_t c = 0; c < chunk_count; ++c) {
    CompositionLaneChunkRecord record{};
    record.gas_cell_id = gas_cell_id;
    record.layout_signature = identity.layout_signature;
    record.exchange_sequence = identity.exchange_sequence;
    record.source_rank = source_rank;
    record.destination_rank = destination_rank;
    record.schema_version = k_composition_wire_schema_version;
    record.stream_kind = static_cast<std::uint32_t>(identity.stream);
    record.lane_count = identity.lane_count;
    record.chunk_index = static_cast<std::uint16_t>(c);
    record.chunk_count = static_cast<std::uint16_t>(chunk_count);
    const std::size_t begin = c * k_composition_lane_chunk_width;
    const std::size_t count = std::min(k_composition_lane_chunk_width, lanes.size() - begin);
    record.chunk_lane_begin = static_cast<std::uint32_t>(begin);
    record.chunk_lane_count = static_cast<std::uint32_t>(count);
    record.scalar = scalar;
    std::copy_n(lanes.begin() + static_cast<std::ptrdiff_t>(begin), count, record.lanes.begin());
    records.push_back(record);
    destinations.push_back(destination_rank);
  }
}

AssembledCompositionLanes assembleCompositionLaneChunks(
    std::span<const CompositionLaneChunkRecord> records,
    const CompositionLaneStreamIdentity& expected,
    int receiver_rank) {
  AssembledCompositionLanes out;
  out.lane_count = expected.lane_count;
  const std::size_t chunk_count = compositionLaneChunkCount(expected.lane_count);
  if (records.empty()) return out;
  if (chunk_count == 0U) {
    throw std::runtime_error("composition lane chunks received for a zero-lane stream");
  }
  for (const auto& r : records) {
    if (r.schema_version != k_composition_wire_schema_version) {
      throw std::runtime_error("composition lane chunk: wire schema version mismatch (got " +
                               std::to_string(r.schema_version) + ")");
    }
    if (r.stream_kind != static_cast<std::uint32_t>(expected.stream)) {
      throw std::runtime_error("composition lane chunk: record belongs to a different stream");
    }
    if (r.layout_signature != expected.layout_signature || r.lane_count != expected.lane_count) {
      throw std::runtime_error("composition lane chunk: sender composition layout differs from the receiver");
    }
    if (r.exchange_sequence != expected.exchange_sequence) {
      throw std::runtime_error("composition lane chunk: stale or foreign exchange sequence");
    }
    if (r.destination_rank != receiver_rank) {
      throw std::runtime_error("composition lane chunk: record delivered to the wrong rank");
    }
    if (r.gas_cell_id == 0U || r.chunk_count != chunk_count || r.chunk_index >= chunk_count) {
      throw std::runtime_error("composition lane chunk: malformed chunk identity");
    }
    const std::size_t begin = static_cast<std::size_t>(r.chunk_index) * k_composition_lane_chunk_width;
    const std::size_t count = std::min(k_composition_lane_chunk_width, expected.lane_count - begin);
    if (r.chunk_lane_begin != begin || r.chunk_lane_count != count) {
      throw std::runtime_error("composition lane chunk: chunk lane range is inconsistent with the layout");
    }
  }
  std::vector<std::size_t> order(records.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&records](std::size_t a, std::size_t b) {
    const auto& l = records[a];
    const auto& r = records[b];
    if (l.gas_cell_id != r.gas_cell_id) return l.gas_cell_id < r.gas_cell_id;
    if (l.source_rank != r.source_rank) return l.source_rank < r.source_rank;
    return l.chunk_index < r.chunk_index;
  });
  if (order.size() % chunk_count != 0U) {
    throw std::runtime_error("composition lane chunk: incomplete chunk set (lost or truncated lanes)");
  }
  const std::size_t entries = order.size() / chunk_count;
  out.gas_cell_id.reserve(entries);
  out.source_rank.reserve(entries);
  out.scalar.reserve(entries);
  out.lanes.resize(core::checkedSizeMultiply(entries, expected.lane_count, "assembled composition lanes"));
  for (std::size_t e = 0; e < entries; ++e) {
    const auto& first = records[order[e * chunk_count]];
    for (std::size_t c = 0; c < chunk_count; ++c) {
      const auto& r = records[order[e * chunk_count + c]];
      if (r.gas_cell_id != first.gas_cell_id || r.source_rank != first.source_rank || r.chunk_index != c) {
        throw std::runtime_error("composition lane chunk: duplicated or missing chunk for cell " +
                                 std::to_string(first.gas_cell_id));
      }
      if (!(r.scalar == first.scalar) && !(std::isnan(r.scalar) && std::isnan(first.scalar))) {
        throw std::runtime_error("composition lane chunk: per-cell scalar differs between chunks");
      }
      std::copy_n(r.lanes.begin(), r.chunk_lane_count,
                  out.lanes.begin() + static_cast<std::ptrdiff_t>(e * expected.lane_count + r.chunk_lane_begin));
    }
    if (!out.gas_cell_id.empty() && out.gas_cell_id.back() == first.gas_cell_id &&
        out.source_rank.back() == first.source_rank) {
      throw std::runtime_error("composition lane chunk: duplicate lane vector for one (cell, source)");
    }
    out.gas_cell_id.push_back(first.gas_cell_id);
    out.source_rank.push_back(first.source_rank);
    out.scalar.push_back(first.scalar);
  }
  return out;
}

namespace {

core::MemoryReservation exchangeIntoDestination(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    core::MemoryGovernor* governor,
    std::string_view owner_label,
    ReceiveAllocator allocate,
    void* allocate_owner) {
  const int world_rank = mpi_context.worldRank();
  const int world_size = std::max(1, mpi_context.worldSize());
  // Every local argument check is made collective before any payload
  // collective so one malformed rank cannot strand its peers.
  std::exception_ptr failure;
  try {
    if (record_bytes == 0U || records.size() != destination_rank.size() * record_bytes) {
      throw std::invalid_argument("destination record exchange: record stream/destination size mismatch");
    }
    for (const int destination : destination_rank) {
      if (destination < 0 || destination >= world_size) {
        throw std::invalid_argument("destination record exchange: destination rank outside MPI world");
      }
    }
  } catch (...) {
    failure = std::current_exception();
  }
  core::MemoryReservation lease;
  if (!mpi_context.isEnabled() || world_size <= 1) {
    if (failure != nullptr) std::rethrow_exception(failure);
    for (const int destination : destination_rank) {
      if (destination != world_rank) {
        throw std::runtime_error("destination record exchange: remote destination without MPI");
      }
    }
    if (governor != nullptr && !records.empty()) {
      lease = governor->reserve(core::MemoryClass::kCommunication, records.size(), owner_label);
      lease.commit();
    }
    std::byte* target = allocate(allocate_owner, records.size());
    if (!records.empty()) std::memcpy(target, records.data(), records.size());
    return lease;
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  mpi_context.rethrowCollectivePreparationFailure(failure, "destination record exchange arguments");
  const std::size_t ranks = static_cast<std::size_t>(world_size);
  std::vector<std::uint64_t> send_counts(ranks, 0U);
  for (const int destination : destination_rank) {
    ++send_counts[static_cast<std::size_t>(destination)];
  }
  std::vector<std::uint64_t> recv_counts(ranks, 0U);
  if (MPI_Alltoall(send_counts.data(), 1, MPI_UINT64_T, recv_counts.data(), 1, MPI_UINT64_T,
                   MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("destination record exchange: count exchange failed");
  }
  // Bounded rounds: per-peer cap so that one round's total staging
  // (world_size * cap) stays below INT_MAX/2 bytes.
  const std::uint64_t int_safe_total = static_cast<std::uint64_t>(std::numeric_limits<int>::max() / 2);
  const std::uint64_t per_peer_bytes_cap = std::max<std::uint64_t>(
      record_bytes, std::min<std::uint64_t>(k_max_destination_exchange_round_bytes_per_peer,
                                            int_safe_total / static_cast<std::uint64_t>(ranks)));
  const std::uint64_t per_peer_records = std::max<std::uint64_t>(1U, per_peer_bytes_cap / record_bytes);
  std::uint64_t local_rounds = 0U;
  std::uint64_t recv_total_records = 0U;
  std::uint64_t max_send_round_bytes = 0U;  // round 0 is the largest round
  std::uint64_t max_recv_round_bytes = 0U;
  for (std::size_t r = 0; r < ranks; ++r) {
    local_rounds = std::max(local_rounds, (send_counts[r] + per_peer_records - 1U) / per_peer_records);
    local_rounds = std::max(local_rounds, (recv_counts[r] + per_peer_records - 1U) / per_peer_records);
    recv_total_records = saturatingAdd(recv_total_records, recv_counts[r]);
    max_send_round_bytes = saturatingAdd(
        max_send_round_bytes, saturatingMul(std::min(per_peer_records, send_counts[r]), record_bytes));
    max_recv_round_bytes = saturatingAdd(
        max_recv_round_bytes, saturatingMul(std::min(per_peer_records, recv_counts[r]), record_bytes));
  }
  const std::uint64_t rounds = mpi_context.allreduceMaxUint64(local_rounds);

  // Admission + EVERY allocation of the exchange happens here, before the
  // collective vote: final receive buffer (lease, returned), bucket index and
  // round staging at maximum round size (scratch reservation, local).
  core::MemoryReservation scratch;
  std::byte* recv_base = nullptr;
  std::vector<std::size_t> bucket;                 // send order bucketed by destination
  std::vector<std::size_t> send_offsets(ranks + 1U, 0U);
  std::vector<std::size_t> recv_offsets(ranks + 1U, 0U);
  std::vector<std::byte> send_stage;
  std::vector<std::byte> recv_stage;
  std::vector<int> send_bytes(ranks);
  std::vector<int> send_displ(ranks);
  std::vector<int> recv_bytes(ranks);
  std::vector<int> recv_displ(ranks);
  failure = nullptr;
  try {
    const std::uint64_t recv_bytes_total = saturatingMul(recv_total_records, record_bytes);
    const std::uint64_t scratch_bytes = saturatingAdd(
        saturatingAdd(max_send_round_bytes, max_recv_round_bytes),
        saturatingMul(destination_rank.size(), sizeof(std::size_t)));
    if (governor != nullptr) {
      if (recv_bytes_total != 0U) {
        lease = governor->reserve(core::MemoryClass::kCommunication, recv_bytes_total, owner_label);
        lease.commit();
      }
      if (scratch_bytes != 0U) {
        scratch = governor->reserve(core::MemoryClass::kCommunication, scratch_bytes, owner_label);
        scratch.commit();
      }
    }
    for (std::size_t r = 0; r < ranks; ++r) {
      send_offsets[r + 1U] = core::checkedSizeAdd(send_offsets[r],
          core::checkedIntegralNarrow<std::size_t>(send_counts[r], "send count"), "send offsets");
      recv_offsets[r + 1U] = core::checkedSizeAdd(recv_offsets[r],
          core::checkedIntegralNarrow<std::size_t>(recv_counts[r], "recv count"), "recv offsets");
    }
    recv_base = allocate(allocate_owner,
                         core::checkedSizeMultiply(recv_offsets[ranks], record_bytes, "recv buffer bytes"));
    bucket.resize(destination_rank.size());
    std::vector<std::size_t> cursor(send_offsets.begin(), send_offsets.end() - 1);
    for (std::size_t i = 0; i < destination_rank.size(); ++i) {
      bucket[cursor[static_cast<std::size_t>(destination_rank[i])]++] = i;
    }
    send_stage.reserve(core::checkedIntegralNarrow<std::size_t>(max_send_round_bytes, "send stage"));
    recv_stage.reserve(core::checkedIntegralNarrow<std::size_t>(max_recv_round_bytes, "recv stage"));
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, "destination record exchange admission");

  for (std::uint64_t round = 0; round < rounds; ++round) {
    const std::uint64_t first = round * per_peer_records;
    std::size_t send_round_total = 0U;
    std::size_t recv_round_total = 0U;
    for (std::size_t r = 0; r < ranks; ++r) {
      const std::uint64_t sc = send_counts[r] > first ? std::min(per_peer_records, send_counts[r] - first) : 0U;
      const std::uint64_t qc = recv_counts[r] > first ? std::min(per_peer_records, recv_counts[r] - first) : 0U;
      send_bytes[r] = static_cast<int>(sc * record_bytes);
      send_displ[r] = static_cast<int>(send_round_total);
      recv_bytes[r] = static_cast<int>(qc * record_bytes);
      recv_displ[r] = static_cast<int>(recv_round_total);
      send_round_total += static_cast<std::size_t>(sc * record_bytes);
      recv_round_total += static_cast<std::size_t>(qc * record_bytes);
    }
    // Within reserved capacity (round 0 is maximal): no reallocation, no
    // throw between collectives.
    send_stage.resize(send_round_total);
    recv_stage.resize(recv_round_total);
    for (std::size_t r = 0; r < ranks; ++r) {
      const std::size_t count = static_cast<std::size_t>(send_bytes[r]) / record_bytes;
      for (std::size_t k = 0; k < count; ++k) {
        const std::size_t source_index = bucket[send_offsets[r] + static_cast<std::size_t>(first) + k];
        std::memcpy(send_stage.data() + static_cast<std::size_t>(send_displ[r]) + k * record_bytes,
                    records.data() + source_index * record_bytes, record_bytes);
      }
    }
    if (MPI_Alltoallv(send_stage.data(), send_bytes.data(), send_displ.data(), MPI_BYTE,
                      recv_stage.data(), recv_bytes.data(), recv_displ.data(), MPI_BYTE,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("destination record exchange: payload exchange failed");
    }
    for (std::size_t r = 0; r < ranks; ++r) {
      if (recv_bytes[r] == 0) continue;
      std::memcpy(recv_base + (recv_offsets[r] + static_cast<std::size_t>(first)) * record_bytes,
                  recv_stage.data() + static_cast<std::size_t>(recv_displ[r]),
                  static_cast<std::size_t>(recv_bytes[r]));
    }
  }
  (void)world_rank;
  return lease;  // staging + scratch reservation released here
#else
  (void)allocate;
  (void)allocate_owner;
  (void)owner_label;
  throw std::runtime_error("destination record exchange requires MPI support when MPI context is enabled");
#endif
}

}  // namespace

std::vector<std::byte> executeBlockingDestinationRecordExchange(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  std::vector<std::byte> out;
  (void)exchangeIntoDestination(
      mpi_context, records, record_bytes, destination_rank, nullptr, "parallel.destination_record_exchange",
      [](void* owner, std::size_t total) -> std::byte* {
        auto& target = *static_cast<std::vector<std::byte>*>(owner);
        target.resize(total);
        return target.data();
      },
      &out);
  return out;
}

GovernedRecordBuffer executeGovernedDestinationRecordExchange(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor,
    std::string_view owner_label) {
  (void)exchange_sequence;
  GovernedRecordBuffer out;
  out.lease = exchangeIntoDestination(
      mpi_context, records, record_bytes, destination_rank, governor, owner_label,
      [](void* owner, std::size_t total) -> std::byte* {
        auto& target = *static_cast<std::vector<std::byte>*>(owner);
        target.resize(total);
        return target.data();
      },
      &out.bytes);
  return out;
}

GovernedRecords<CompositionDonorRequestRecord> exchangeCompositionDonorRequests(
    const MpiContext& mpi_context,
    std::span<const CompositionDonorRequestRecord> requests,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor) {
  std::vector<int> destinations;
  destinations.reserve(requests.size());
  std::exception_ptr failure;
  try {
    for (const auto& request : requests) {
      if (request.gas_cell_id == 0U || request.requester_rank != mpi_context.worldRank()) {
        throw std::invalid_argument("composition donor request is malformed or not from this rank");
      }
      destinations.push_back(request.owner_rank);
    }
  } catch (...) {
    failure = std::current_exception();
    destinations.clear();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, "composition donor request validation");
  (void)exchange_sequence;
  return exchangeGovernedTyped(mpi_context, requests, std::span<const int>(destinations), governor,
                               "parallel.composition_donor_requests");
}

GovernedRecords<CompositionLaneChunkRecord> exchangeCompositionLaneChunks(
    const MpiContext& mpi_context,
    std::span<const CompositionLaneChunkRecord> records,
    std::span<const int> destinations,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* governor,
    std::string_view owner_label) {
  std::exception_ptr failure;
  try {
    if (records.size() != destinations.size()) {
      throw std::invalid_argument("composition lane exchange: record/destination count mismatch");
    }
    for (std::size_t i = 0; i < records.size(); ++i) {
      const auto& r = records[i];
      if (r.schema_version != k_composition_wire_schema_version || r.gas_cell_id == 0U ||
          r.source_rank != mpi_context.worldRank() || r.destination_rank != destinations[i] ||
          r.exchange_sequence != exchange_sequence || r.lane_count > k_max_composition_layout_lanes ||
          r.chunk_lane_count > k_composition_lane_chunk_width) {
        throw std::invalid_argument("composition lane chunk record is malformed or not from this rank");
      }
    }
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, "composition lane chunk validation");
  return exchangeGovernedTyped(mpi_context, records, destinations, governor, owner_label);
}

}  // namespace cosmosim::parallel
