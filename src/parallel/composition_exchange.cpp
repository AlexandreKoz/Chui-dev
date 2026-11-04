#include "cosmosim/parallel/composition_exchange.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
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

template <class Record>
[[nodiscard]] std::vector<Record> exchangeTyped(
    const MpiContext& mpi_context,
    std::span<const Record> records,
    std::vector<int> destinations,
    std::uint64_t exchange_sequence) {
  static_assert(std::is_trivially_copyable_v<Record>);
  const std::span<const std::byte> bytes(
      reinterpret_cast<const std::byte*>(records.data()), records.size() * sizeof(Record));
  const std::vector<std::byte> received = executeBlockingDestinationRecordExchange(
      mpi_context, bytes, sizeof(Record), destinations, exchange_sequence);
  if (received.size() % sizeof(Record) != 0U) {
    throw std::runtime_error("composition exchange received a truncated record stream");
  }
  std::vector<Record> result(received.size() / sizeof(Record));
  if (!result.empty()) {
    std::memcpy(result.data(), received.data(), received.size());
  }
  return result;
}

}  // namespace

std::vector<std::byte> executeBlockingDestinationRecordExchange(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  if (record_bytes == 0U || records.size() != destination_rank.size() * record_bytes) {
    throw std::invalid_argument("destination record exchange: record stream/destination size mismatch");
  }
  const int world_rank = mpi_context.worldRank();
  const int world_size = std::max(1, mpi_context.worldSize());
  for (const int destination : destination_rank) {
    if (destination < 0 || destination >= world_size) {
      throw std::invalid_argument("destination record exchange: destination rank outside MPI world");
    }
  }
  if (!mpi_context.isEnabled() || world_size <= 1) {
    for (const int destination : destination_rank) {
      if (destination != world_rank) {
        throw std::runtime_error("destination record exchange: remote destination without MPI");
      }
    }
    return std::vector<std::byte>(records.begin(), records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
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
  // Stable bucketing by destination preserves the caller's send order.
  std::vector<std::size_t> send_offsets(ranks + 1U, 0U);
  for (std::size_t r = 0; r < ranks; ++r) {
    send_offsets[r + 1U] = core::checkedSizeAdd(send_offsets[r],
        core::checkedIntegralNarrow<std::size_t>(send_counts[r], "send count"), "send offsets");
  }
  std::vector<std::byte> send_buffer(records.size());
  std::vector<std::size_t> cursor(send_offsets.begin(), send_offsets.end() - 1);
  for (std::size_t i = 0; i < destination_rank.size(); ++i) {
    const std::size_t slot = cursor[static_cast<std::size_t>(destination_rank[i])]++;
    std::memcpy(send_buffer.data() + slot * record_bytes, records.data() + i * record_bytes, record_bytes);
  }
  std::vector<int> send_bytes(ranks);
  std::vector<int> send_displ(ranks);
  std::vector<int> recv_bytes(ranks);
  std::vector<int> recv_displ(ranks);
  std::size_t total_recv = 0U;
  for (std::size_t r = 0; r < ranks; ++r) {
    send_bytes[r] = core::checkedIntegralNarrow<int>(
        core::checkedSizeMultiply(static_cast<std::size_t>(send_counts[r]), record_bytes, "send bytes"), "send bytes");
    send_displ[r] = core::checkedIntegralNarrow<int>(
        core::checkedSizeMultiply(send_offsets[r], record_bytes, "send displ"), "send displ");
    recv_bytes[r] = core::checkedIntegralNarrow<int>(
        core::checkedSizeMultiply(static_cast<std::size_t>(recv_counts[r]), record_bytes, "recv bytes"), "recv bytes");
    recv_displ[r] = core::checkedIntegralNarrow<int>(total_recv, "recv displ");
    total_recv = core::checkedSizeAdd(total_recv, static_cast<std::size_t>(recv_bytes[r]), "recv total");
  }
  std::vector<std::byte> recv_buffer(total_recv);
  if (MPI_Alltoallv(send_buffer.data(), send_bytes.data(), send_displ.data(), MPI_BYTE,
                    recv_buffer.data(), recv_bytes.data(), recv_displ.data(), MPI_BYTE,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("destination record exchange: payload exchange failed");
  }
  (void)world_rank;
  return recv_buffer;
#else
  throw std::runtime_error("destination record exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<CompositionDonorRequestRecord> exchangeCompositionDonorRequests(
    const MpiContext& mpi_context,
    std::span<const CompositionDonorRequestRecord> requests,
    std::uint64_t exchange_sequence) {
  std::vector<int> destinations;
  destinations.reserve(requests.size());
  for (const auto& request : requests) {
    if (request.gas_cell_id == 0U || request.requester_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("composition donor request is malformed or not from this rank");
    }
    destinations.push_back(request.owner_rank);
  }
  return exchangeTyped(mpi_context, requests, std::move(destinations), exchange_sequence);
}

std::vector<CompositionDonorFractionRecord> exchangeCompositionDonorFractions(
    const MpiContext& mpi_context,
    std::span<const CompositionDonorFractionRecord> responses,
    std::uint64_t exchange_sequence) {
  std::vector<int> destinations;
  destinations.reserve(responses.size());
  for (const auto& response : responses) {
    if (response.gas_cell_id == 0U || response.owner_rank != mpi_context.worldRank() ||
        response.lane_count > k_max_composition_lanes) {
      throw std::invalid_argument("composition donor fraction record is malformed or not owned here");
    }
    destinations.push_back(response.requester_rank);
  }
  return exchangeTyped(mpi_context, responses, std::move(destinations), exchange_sequence);
}

std::vector<CompositionCorrectionRecord> exchangeCompositionCorrections(
    const MpiContext& mpi_context,
    std::span<const CompositionCorrectionRecord> corrections,
    std::uint64_t exchange_sequence) {
  std::vector<int> destinations;
  destinations.reserve(corrections.size());
  for (const auto& correction : corrections) {
    if (correction.gas_cell_id == 0U || correction.source_rank != mpi_context.worldRank() ||
        correction.lane_count > k_max_composition_lanes) {
      throw std::invalid_argument("composition correction record is malformed or not from this rank");
    }
    destinations.push_back(correction.owner_rank);
  }
  return exchangeTyped(mpi_context, corrections, std::move(destinations), exchange_sequence);
}

}  // namespace cosmosim::parallel
