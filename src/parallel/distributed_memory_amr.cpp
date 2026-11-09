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
DirectedAmrPatchCellTransferPlan planDirectedAmrPatchCellTransfer(
    std::uint64_t logical_record_count,
    std::size_t transport_round_limit_bytes) {
  const std::size_t requested_limit = transport_round_limit_bytes == 0U
      ? mpiTransportRoundLimitBytes()
      : transport_round_limit_bytes;
  const std::size_t bounded_limit = std::min(
      requested_limit,
      static_cast<std::size_t>(std::numeric_limits<int>::max()));
  const std::size_t records_per_round = bounded_limit / sizeof(AmrPatchCellPayloadRecord);
  if (records_per_round == 0U) {
    throw std::invalid_argument(
        "directed AMR patch-cell transport round cannot hold one complete record");
  }
  const std::uint64_t round_count_u64 = logical_record_count / records_per_round +
      (logical_record_count % records_per_round == 0U ? 0U : 1U);
  return DirectedAmrPatchCellTransferPlan{
      .logical_record_count = logical_record_count,
      .records_per_round = records_per_round,
      .round_count = round_count_u64,
  };
}

std::vector<std::size_t> planDirectedAmrPatchCellTransferRounds(
    std::uint64_t logical_record_count,
    std::size_t transport_round_limit_bytes) {
  constexpr std::uint64_t k_max_materialized_round_count = 65'536U;
  const DirectedAmrPatchCellTransferPlan plan = planDirectedAmrPatchCellTransfer(
      logical_record_count, transport_round_limit_bytes);
  if (plan.round_count > k_max_materialized_round_count) {
    throw std::length_error(
        "directed AMR diagnostic round materialization exceeds bounded metadata cap; "
        "use the constant-space transfer plan");
  }
  const std::size_t round_count = core::checkedIntegralNarrow<std::size_t>(
      plan.round_count, "directed AMR patch-cell materialized round count");
  std::vector<std::size_t> rounds;
  rounds.reserve(round_count);
  std::uint64_t remaining = plan.logical_record_count;
  while (remaining != 0U) {
    const std::size_t records = static_cast<std::size_t>(std::min<std::uint64_t>(
        remaining, static_cast<std::uint64_t>(plan.records_per_round)));
    rounds.push_back(records);
    remaining -= records;
  }
  return rounds;
}


namespace {

// Conservative per-entry overhead estimate of one std::unordered_map node plus
// its bucket slot (pointer-linked node + cached hash + bucket pointer); used
// only for governor admission of transient indices, not as a measurement.
inline constexpr std::size_t k_directed_amr_hash_node_overhead_bytes = 48U;

struct AmrRankEnvelopeRecord {
  int rank = -1;
  std::uint8_t has_patches = 0;
  std::uint64_t decomposition_epoch = 0;
  // Patch 53: number of local authoritative patches; bounds the descriptors a
  // peer can send, so receivers admit their staged result before streaming.
  std::uint64_t patch_count = 0;
  double min_x_comoving = 0.0;
  double max_x_comoving = 0.0;
  double min_y_comoving = 0.0;
  double max_y_comoving = 0.0;
  double min_z_comoving = 0.0;
  double max_z_comoving = 0.0;
  double max_cell_width_comoving = 0.0;
};

[[nodiscard]] double amrPatchMaxX(const AmrPatchPayloadRecord& record) noexcept {
  return record.origin_x_comoving + record.extent_x_comoving;
}

[[nodiscard]] double amrPatchMaxY(const AmrPatchPayloadRecord& record) noexcept {
  return record.origin_y_comoving + record.extent_y_comoving;
}

[[nodiscard]] double amrPatchMaxZ(const AmrPatchPayloadRecord& record) noexcept {
  return record.origin_z_comoving + record.extent_z_comoving;
}

[[nodiscard]] double amrPatchMaxCellWidth(const AmrPatchPayloadRecord& record) noexcept {
  const double dx = record.extent_x_comoving / static_cast<double>(std::max<std::uint16_t>(record.cell_dim_x, 1U));
  const double dy = record.extent_y_comoving / static_cast<double>(std::max<std::uint16_t>(record.cell_dim_y, 1U));
  const double dz = record.extent_z_comoving / static_cast<double>(std::max<std::uint16_t>(record.cell_dim_z, 1U));
  return std::max({dx, dy, dz});
}

[[nodiscard, maybe_unused]] std::uint64_t checkedU64AddLocalEarly(std::uint64_t lhs, std::uint64_t rhs, const char* label) {
  if (lhs > std::numeric_limits<std::uint64_t>::max() - rhs) {
    throw std::overflow_error(std::string(label) + " overflows uint64_t");
  }
  return lhs + rhs;
}

[[nodiscard, maybe_unused]] AmrRankEnvelopeRecord buildLocalAmrEnvelope(
    std::span<const AmrPatchPayloadRecord> records,
    int world_rank) {
  AmrRankEnvelopeRecord envelope;
  envelope.rank = world_rank;
  if (records.empty()) {
    return envelope;
  }
  envelope.has_patches = 1U;
  envelope.patch_count = static_cast<std::uint64_t>(records.size());
  envelope.decomposition_epoch = records.front().decomposition_epoch;
  envelope.min_x_comoving = std::numeric_limits<double>::infinity();
  envelope.min_y_comoving = std::numeric_limits<double>::infinity();
  envelope.min_z_comoving = std::numeric_limits<double>::infinity();
  envelope.max_x_comoving = -std::numeric_limits<double>::infinity();
  envelope.max_y_comoving = -std::numeric_limits<double>::infinity();
  envelope.max_z_comoving = -std::numeric_limits<double>::infinity();
  for (const AmrPatchPayloadRecord& record : records) {
    validateAmrPatchPayloadRecord(record);
    if (record.owner_rank != world_rank) {
      throw std::invalid_argument("directed AMR envelope can only be built from local authoritative patches");
    }
    if (record.decomposition_epoch != envelope.decomposition_epoch) {
      throw std::invalid_argument("directed AMR envelope found mixed decomposition epochs in one exchange");
    }
    envelope.min_x_comoving = std::min(envelope.min_x_comoving, record.origin_x_comoving);
    envelope.min_y_comoving = std::min(envelope.min_y_comoving, record.origin_y_comoving);
    envelope.min_z_comoving = std::min(envelope.min_z_comoving, record.origin_z_comoving);
    envelope.max_x_comoving = std::max(envelope.max_x_comoving, amrPatchMaxX(record));
    envelope.max_y_comoving = std::max(envelope.max_y_comoving, amrPatchMaxY(record));
    envelope.max_z_comoving = std::max(envelope.max_z_comoving, amrPatchMaxZ(record));
    envelope.max_cell_width_comoving = std::max(envelope.max_cell_width_comoving, amrPatchMaxCellWidth(record));
  }
  return envelope;
}

[[nodiscard]] bool amrIntervalsOverlap(double amin, double amax, double bmin, double bmax, double reach) noexcept {
  return (amin - reach) <= bmax && (bmin - reach) <= amax;
}

// CHUI Part 3: periodic image overlap. On a periodic axis the peer interval
// is also tested shifted by +-L, so ranks that meet only through the box wrap
// are discovered as peers.
[[nodiscard]] bool amrIntervalsOverlapPeriodic(
    double amin, double amax, double bmin, double bmax, double reach, double length) noexcept {
  if (amrIntervalsOverlap(amin, amax, bmin, bmax, reach)) {
    return true;
  }
  if (!(length > 0.0)) {
    return false;
  }
  return amrIntervalsOverlap(amin, amax, bmin + length, bmax + length, reach) ||
      amrIntervalsOverlap(amin, amax, bmin - length, bmax - length, reach);
}

[[nodiscard]] double topologyAxisLength(const AmrInterfaceTopologyOptions& topology, std::size_t axis) noexcept {
  const double length = topology.periodic_box_extent_comoving[axis];
  return (length > 0.0 && std::isfinite(length)) ? length : 0.0;
}

// Reach in cells: deep ghost layers need up to ghost_layers source widths.
[[nodiscard]] double topologyReachScale(const AmrInterfaceTopologyOptions& topology) noexcept {
  return static_cast<double>(std::max<std::uint16_t>(topology.ghost_layers, 1U));
}

[[nodiscard, maybe_unused]] bool amrEnvelopeMayNeedPeer(
    const AmrRankEnvelopeRecord& local,
    const AmrRankEnvelopeRecord& peer,
    const AmrInterfaceTopologyOptions& topology = {}) noexcept {
  if (local.has_patches == 0U || peer.has_patches == 0U || local.rank == peer.rank) {
    return false;
  }
  const double reach = topologyReachScale(topology) *
      std::max(local.max_cell_width_comoving, peer.max_cell_width_comoving);
  return amrIntervalsOverlapPeriodic(local.min_x_comoving, local.max_x_comoving, peer.min_x_comoving, peer.max_x_comoving, reach, topologyAxisLength(topology, 0U)) &&
      amrIntervalsOverlapPeriodic(local.min_y_comoving, local.max_y_comoving, peer.min_y_comoving, peer.max_y_comoving, reach, topologyAxisLength(topology, 1U)) &&
      amrIntervalsOverlapPeriodic(local.min_z_comoving, local.max_z_comoving, peer.min_z_comoving, peer.max_z_comoving, reach, topologyAxisLength(topology, 2U));
}

[[nodiscard]] bool amrPatchesMayShareInterface(
    const AmrPatchPayloadRecord& lhs,
    const AmrPatchPayloadRecord& rhs) noexcept {
  const double reach = std::max(amrPatchMaxCellWidth(lhs), amrPatchMaxCellWidth(rhs));
  return amrIntervalsOverlap(lhs.origin_x_comoving, amrPatchMaxX(lhs), rhs.origin_x_comoving, amrPatchMaxX(rhs), reach) &&
      amrIntervalsOverlap(lhs.origin_y_comoving, amrPatchMaxY(lhs), rhs.origin_y_comoving, amrPatchMaxY(rhs), reach) &&
      amrIntervalsOverlap(lhs.origin_z_comoving, amrPatchMaxZ(lhs), rhs.origin_z_comoving, amrPatchMaxZ(rhs), reach);
}

// Patch 53 (P0-52-01): collective-safe candidate-peer discovery.
struct AmrCandidatePeerSet {
  std::vector<int> peers;                          // ascending
  std::uint64_t remote_patch_upper_bound = 0U;     // sum of candidate peers' patch counts
  std::uint64_t max_peer_patch_count = 0U;
};

[[nodiscard, maybe_unused]] AmrCandidatePeerSet discoverCandidateAmrPeers(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_patch_records,
    DirectedAmrExchangeDiagnostics* diagnostics,
    const AmrInterfaceTopologyOptions& topology,
    core::MemoryGovernor* memory_governor,
    core::MemoryReservation& control_lease,
    std::exception_ptr prior_local_failure = nullptr) {
  // prior_local_failure: the caller's own E0 preparation failure, folded into
  // the first vote so it is agreed before the Allgather (one collective).
  const int world_size = mpi_context.worldSize();
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)local_patch_records;
  (void)topology;
  (void)memory_governor;
  (void)control_lease;
#endif
  if (!mpi_context.isEnabled() || world_size <= 1) {
    if (prior_local_failure != nullptr) std::rethrow_exception(prior_local_failure);
    if (diagnostics != nullptr) {
      diagnostics->control_plane_bytes += sizeof(AmrRankEnvelopeRecord);
    }
    return {};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_rank = mpi_context.worldRank();
  static_assert(std::is_trivially_copyable_v<AmrRankEnvelopeRecord>);
  // E0: envelope (validates every local record), control-plane admission and
  // allocation. Nothing that can throw remains between the vote and the
  // Allgather (the int narrowing of a compile-time record size cannot fail).
  AmrRankEnvelopeRecord local_envelope;
  std::vector<AmrRankEnvelopeRecord> envelopes;
  AmrCandidatePeerSet result;
  std::exception_ptr failure = prior_local_failure;
  try {
    if (failure != nullptr) std::rethrow_exception(failure);
    local_envelope = buildLocalAmrEnvelope(local_patch_records, world_rank);
    const std::size_t ranks = static_cast<std::size_t>(world_size);
    const std::size_t control_bytes = core::checkedSizeAdd(
        core::checkedSizeMultiply(ranks, sizeof(AmrRankEnvelopeRecord), "directed AMR envelope bytes"),
        core::checkedSizeMultiply(ranks, sizeof(int), "directed AMR candidate peer bytes"),
        "directed AMR control-plane bytes");
    if (memory_governor != nullptr) {
      control_lease = memory_governor->reserve(
          core::MemoryClass::kCommunication, control_bytes, "parallel.directed_amr.control_plane");
      control_lease.commit();
    }
    envelopes.resize(ranks);
    result.peers.reserve(ranks);
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, "directed AMR envelope preflight");
  if (MPI_Allgather(
          &local_envelope,
          static_cast<int>(sizeof(AmrRankEnvelopeRecord)),
          MPI_BYTE,
          envelopes.data(),
          static_cast<int>(sizeof(AmrRankEnvelopeRecord)),
          MPI_BYTE,
          MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR control-plane envelope Allgather failed");
  }
  if (diagnostics != nullptr) {
    diagnostics->control_plane_bytes += static_cast<std::uint64_t>(world_size) * sizeof(AmrRankEnvelopeRecord);
  }
  // E2 (first half): symmetric peer graph. Every pair is decided on the
  // rank-ordered envelopes (lower rank first) so both members evaluate the
  // bitwise-identical floating-point expression and agree on the edge.
  try {
    for (std::size_t r = 0; r < envelopes.size(); ++r) {
      if (envelopes[r].rank != static_cast<int>(r)) {
        throw std::runtime_error("directed AMR control-plane envelope returned invalid rank metadata");
      }
    }
    const AmrRankEnvelopeRecord& self = envelopes[static_cast<std::size_t>(world_rank)];
    for (int peer = 0; peer < world_size; ++peer) {
      if (peer == world_rank) continue;
      const AmrRankEnvelopeRecord& other = envelopes[static_cast<std::size_t>(peer)];
      const bool edge = peer < world_rank ? amrEnvelopeMayNeedPeer(other, self, topology)
                                          : amrEnvelopeMayNeedPeer(self, other, topology);
      if (!edge) continue;
      result.peers.push_back(peer);  // within reserved capacity
      result.remote_patch_upper_bound = checkedU64AddLocalEarly(
          result.remote_patch_upper_bound, other.patch_count, "directed AMR remote patch upper bound");
      result.max_peer_patch_count = std::max(result.max_peer_patch_count, other.patch_count);
    }
  } catch (...) {
    failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(failure, "directed AMR candidate peer graph");
  if (diagnostics != nullptr) {
    diagnostics->candidate_peer_count = static_cast<std::uint64_t>(result.peers.size());
  }
  return result;
#else
  throw std::runtime_error("directed AMR peer discovery requires MPI support when MPI context is enabled");
#endif
}

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
// Patch 53 pair-epoch health handshake. Returns the peer's health. An MPI
// error return is fatal for the communicator and is rethrown as such.
[[nodiscard]] bool exchangePairHealth(
    const MpiContext& mpi_context,
    int peer_rank,
    bool local_healthy,
    int tag_base,
    std::uint64_t exchange_sequence) {
  int local_flag = local_healthy ? 1 : 0;
  int peer_flag = 0;
  const int tag = ghostExchangeSequencedTag(tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence);
  if (MPI_Sendrecv(&local_flag, 1, MPI_INT, peer_rank, tag, &peer_flag, 1, MPI_INT, peer_rank, tag,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR pair-health Sendrecv failed (communicator-fatal)");
  }
  return peer_flag != 0;
}
#endif

template <typename T>
[[nodiscard]] std::vector<T> exchangePodRecordsWithPeer(
    const MpiContext& mpi_context,
    int peer_rank,
    std::span<const T> local_records,
    int count_tag_base,
    int payload_tag_base,
    std::uint64_t exchange_sequence,
    const char* caller,
    std::uint64_t* transport_round_count = nullptr,
    core::MemoryGovernor* memory_governor = nullptr,
    core::MemoryReservation* receive_lease = nullptr) {
  // Pair-epoch contract (patch 53): both ranks traverse the same message
  // sequence; every fallible local step precedes the readiness message, and
  // after it only bounded, non-throwing arithmetic and payload messages run.
  static_assert(std::is_trivially_copyable_v<T>);
  if (peer_rank < 0 || peer_rank >= mpi_context.worldSize() || peer_rank == mpi_context.worldRank()) {
    throw std::invalid_argument(std::string(caller) + ": invalid peer rank");
  }
  if (!mpi_context.isEnabled()) {
    if (!local_records.empty()) {
      throw std::runtime_error(std::string(caller) + ": non-empty peer exchange requires MPI");
    }
    return {};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const std::uint64_t send_count = static_cast<std::uint64_t>(local_records.size());
  std::uint64_t recv_count = 0;
  if (MPI_Sendrecv(
          const_cast<std::uint64_t*>(&send_count),
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &recv_count,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error(std::string(caller) + ": directed AMR record-count Sendrecv failed");
  }

  std::vector<T> received;
  std::size_t receive_count = 0U;
  std::size_t records_per_round = 0U;
  std::exception_ptr local_preparation_failure;
  try {
    receive_count = core::checkedIntegralNarrow<std::size_t>(
        recv_count, std::string(caller) + " receive record count");
    if (memory_governor != nullptr && receive_lease != nullptr && receive_count != 0U) {
      *receive_lease = memory_governor->reserve(
          core::MemoryClass::kCommunication,
          core::checkedSizeMultiply(receive_count, sizeof(T), std::string(caller) + " receive bytes"),
          caller);
      receive_lease->commit();
    }
    received.resize(receive_count);
    const std::size_t round_limit_bytes = std::min(
        mpiTransportRoundLimitBytes(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()));
    records_per_round = round_limit_bytes / sizeof(T);
    if (records_per_round == 0U) {
      throw std::overflow_error(
          std::string(caller) + ": one record exceeds the bounded MPI transport round");
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }

  int local_ready = local_preparation_failure == nullptr ? 1 : 0;
  int peer_ready = 0;
  if (MPI_Sendrecv(
          &local_ready,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(
              count_tag_base + 100, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &peer_ready,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(
              count_tag_base + 100, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error(std::string(caller) + ": directed AMR preparation-readiness Sendrecv failed");
  }
  if (local_ready == 0) {
    std::rethrow_exception(local_preparation_failure);
  }
  if (peer_ready == 0) {
    throw std::runtime_error(std::string(caller) + ": peer rejected directed AMR receive preparation");
  }

  std::size_t send_offset = 0U;
  std::size_t recv_offset = 0U;
  std::uint64_t rounds = 0U;
  while (send_offset < local_records.size() || recv_offset < received.size()) {
    const std::size_t send_records = std::min(records_per_round, local_records.size() - send_offset);
    const std::size_t recv_records = std::min(records_per_round, received.size() - recv_offset);
    const std::size_t send_bytes = core::checkedSizeMultiply(
        send_records, sizeof(T), std::string(caller) + " bounded send bytes");
    const std::size_t recv_bytes = core::checkedSizeMultiply(
        recv_records, sizeof(T), std::string(caller) + " bounded receive bytes");
    const int send_count_bytes = core::checkedIntegralNarrow<int>(
        send_bytes, std::string(caller) + " bounded send count");
    const int recv_count_bytes = core::checkedIntegralNarrow<int>(
        recv_bytes, std::string(caller) + " bounded receive count");
    T* recv_pointer = recv_records == 0U ? nullptr : received.data() + recv_offset;
    const T* send_pointer = send_records == 0U ? nullptr : local_records.data() + send_offset;
    if (MPI_Sendrecv(
            const_cast<T*>(send_pointer),
            send_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            recv_pointer,
            recv_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(std::string(caller) + ": bounded directed AMR payload Sendrecv failed");
    }
    send_offset += send_records;
    recv_offset += recv_records;
    ++rounds;
  }
  if (transport_round_count != nullptr) {
    *transport_round_count += rounds;
  }
  return received;
#else
  throw std::runtime_error(std::string(caller) + ": MPI support is not compiled in");
#endif
}

void exchangeAmrPatchCellRecordStreamWithPeer(
    const MpiContext& mpi_context,
    int peer_rank,
    std::uint64_t logical_send_count,
    std::span<const AmrPatchBoundaryCellRequest> local_requests,
    std::span<const AmrPatchPayloadRecord> remote_patches,
    std::span<const AmrPatchBoundaryCellRequest> remote_requests,
    const DirectedAmrPatchCellPayloadProvider& provider,
    const DirectedAmrPatchCellPayloadConsumer& consumer,
    const DirectedAmrPatchCellAdmission& admission,
    int count_tag_base,
    int payload_tag_base,
    std::uint64_t exchange_sequence,
    std::size_t local_transport_round_limit_bytes,
    DirectedAmrExchangeDiagnostics& diagnostics,
    core::MemoryGovernor* chunk_governor = nullptr) {
  // chunk_governor != nullptr: this call admits its two streaming chunks
  // (the caller did not size/admit the transport window itself).
  if (!provider || !consumer) {
    throw std::invalid_argument(
        "directed AMR patch-cell streaming requires both producer and consumer callbacks");
  }
  if (!mpi_context.isEnabled()) {
    if (logical_send_count != 0U) {
      throw std::runtime_error(
          "directed AMR patch-cell streaming requires MPI for non-empty payloads");
    }
    return;
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::uint64_t logical_receive_count = 0U;
  if (MPI_Sendrecv(
          &logical_send_count,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &logical_receive_count,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR patch-cell record-count Sendrecv failed");
  }

  const std::size_t local_limit = std::min(
      local_transport_round_limit_bytes == 0U
          ? mpiTransportRoundLimitBytes()
          : local_transport_round_limit_bytes,
      static_cast<std::size_t>(std::numeric_limits<int>::max()));
  std::uint64_t local_limit_u64 = static_cast<std::uint64_t>(local_limit);
  std::uint64_t peer_limit_u64 = 0U;
  if (MPI_Sendrecv(
          &local_limit_u64,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 1, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &peer_limit_u64,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 1, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR patch-cell transport-limit Sendrecv failed");
  }
  DirectedAmrPatchCellTransferPlan send_plan;
  DirectedAmrPatchCellTransferPlan receive_plan;
  // Patch 53: EVERY fallible step before the first payload round is inside
  // this guarded phase (limit narrowing, planning, caller admission, chunk
  // admission and capacity, diagnostics arithmetic). After the pairwise
  // readiness message the round loop only resizes within reserved capacity.
  core::MemoryReservation chunk_lease;
  std::vector<AmrPatchCellPayloadRecord> send_chunk;
  std::vector<AmrPatchCellPayloadRecord> receive_chunk;
  std::exception_ptr local_admission_failure;
  try {
    const std::size_t agreed_limit = std::min(
        local_limit,
        core::checkedIntegralNarrow<std::size_t>(
            peer_limit_u64, "directed AMR peer transport limit"));
    send_plan = planDirectedAmrPatchCellTransfer(
        logical_send_count, agreed_limit);
    receive_plan = planDirectedAmrPatchCellTransfer(
        logical_receive_count, agreed_limit);
    if (send_plan.records_per_round != receive_plan.records_per_round) {
      throw std::logic_error("directed AMR send/receive round plans disagree on the agreed limit");
    }
    if (admission) {
      admission(peer_rank, remote_patches, remote_requests, logical_receive_count);
    }
    const std::size_t max_records = send_plan.records_per_round;
    const std::size_t chunk_bytes = core::checkedSizeMultiply(
        max_records, sizeof(AmrPatchCellPayloadRecord), "directed AMR streamed chunk bytes");
    if (chunk_governor != nullptr && chunk_bytes != 0U) {
      chunk_lease = chunk_governor->reserve(
          core::MemoryClass::kCommunication,
          core::checkedSizeMultiply(chunk_bytes, 2U, "directed AMR streamed chunk pair bytes"),
          "parallel.directed_amr.cell_stream_chunks");
      chunk_lease.commit();
    }
    send_chunk.reserve(max_records);
    receive_chunk.reserve(max_records);
    diagnostics.patch_cell_send_capacity_high_water_bytes = std::max(
        diagnostics.patch_cell_send_capacity_high_water_bytes,
        static_cast<std::uint64_t>(core::checkedSizeMultiply(
            send_chunk.capacity(), sizeof(AmrPatchCellPayloadRecord),
            "directed AMR streamed send capacity")));
    diagnostics.patch_cell_receive_capacity_high_water_bytes = std::max(
        diagnostics.patch_cell_receive_capacity_high_water_bytes,
        static_cast<std::uint64_t>(core::checkedSizeMultiply(
            receive_chunk.capacity(), sizeof(AmrPatchCellPayloadRecord),
            "directed AMR streamed receive capacity")));
    diagnostics.communication_workspace_high_water_bytes = std::max(
        diagnostics.communication_workspace_high_water_bytes,
        core::checkedMemoryBytesAdd(
            diagnostics.patch_cell_send_capacity_high_water_bytes,
            diagnostics.patch_cell_receive_capacity_high_water_bytes,
            "directed AMR streamed simultaneous transport capacity"));
  } catch (...) {
    local_admission_failure = std::current_exception();
  }
  int local_admitted = local_admission_failure == nullptr ? 1 : 0;
  int peer_admitted = 0;
  if (MPI_Sendrecv(
          &local_admitted,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 2, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &peer_admitted,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 2, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR patch-cell admission-readiness Sendrecv failed");
  }
  if (local_admitted == 0) {
    std::rethrow_exception(local_admission_failure);
  }
  if (peer_admitted == 0) {
    throw std::runtime_error("directed AMR peer rejected patch-cell memory admission");
  }

  const std::uint64_t round_count = std::max(send_plan.round_count, receive_plan.round_count);

  std::uint64_t send_offset = 0U;
  std::uint64_t receive_offset = 0U;
  for (std::uint64_t round = 0U; round < round_count; ++round) {
    const std::size_t send_records = round < send_plan.round_count
        ? static_cast<std::size_t>(std::min<std::uint64_t>(
              send_plan.logical_record_count - send_offset,
              static_cast<std::uint64_t>(send_plan.records_per_round)))
        : 0U;
    const std::size_t receive_records = round < receive_plan.round_count
        ? static_cast<std::size_t>(std::min<std::uint64_t>(
              receive_plan.logical_record_count - receive_offset,
              static_cast<std::uint64_t>(receive_plan.records_per_round)))
        : 0U;

    std::exception_ptr local_producer_failure;
    try {
      send_chunk.clear();
      if (send_records != 0U) {
        provider(local_requests, send_offset, send_records, send_chunk);
      }
      if (send_chunk.size() != send_records) {
        throw std::runtime_error(
            "directed AMR patch-cell producer did not provide the requested streamed record count");
      }
    } catch (...) {
      local_producer_failure = std::current_exception();
    }
    int local_ready = local_producer_failure == nullptr ? 1 : 0;
    int peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 3, mpi_context.worldRank(), peer_rank, exchange_sequence),
            &peer_ready,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 3, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error("directed AMR patch-cell producer-readiness Sendrecv failed");
    }
    if (local_ready == 0) {
      std::rethrow_exception(local_producer_failure);
    }
    if (peer_ready == 0) {
      throw std::runtime_error("directed AMR peer rejected patch-cell producer preparation");
    }

    receive_chunk.resize(receive_records);
    const std::size_t send_bytes = core::checkedSizeMultiply(
        send_records, sizeof(AmrPatchCellPayloadRecord),
        "directed AMR streamed send bytes");
    const std::size_t receive_bytes = core::checkedSizeMultiply(
        receive_records, sizeof(AmrPatchCellPayloadRecord),
        "directed AMR streamed receive bytes");
    const int send_count_bytes = core::checkedIntegralNarrow<int>(
        send_bytes, "directed AMR streamed send count");
    const int receive_count_bytes = core::checkedIntegralNarrow<int>(
        receive_bytes, "directed AMR streamed receive count");
    if (MPI_Sendrecv(
            send_records == 0U ? nullptr : send_chunk.data(),
            send_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            receive_records == 0U ? nullptr : receive_chunk.data(),
            receive_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error("directed AMR bounded streamed patch-cell Sendrecv failed");
    }

    std::exception_ptr local_consumer_failure;
    try {
      if (!receive_chunk.empty()) {
        consumer(peer_rank, receive_chunk);
      }
    } catch (...) {
      local_consumer_failure = std::current_exception();
    }
    int local_consumed = local_consumer_failure == nullptr ? 1 : 0;
    int peer_consumed = 0;
    if (MPI_Sendrecv(
            &local_consumed,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 4, mpi_context.worldRank(), peer_rank, exchange_sequence),
            &peer_consumed,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 4, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error("directed AMR patch-cell consumer-readiness Sendrecv failed");
    }
    if (local_consumed == 0) {
      std::rethrow_exception(local_consumer_failure);
    }
    if (peer_consumed == 0) {
      throw std::runtime_error("directed AMR peer rejected streamed patch-cell consumption");
    }

    send_offset += static_cast<std::uint64_t>(send_records);
    receive_offset += static_cast<std::uint64_t>(receive_records);
    ++diagnostics.patch_cell_transport_round_count;
  }
  if (send_offset != logical_send_count) {
    throw std::logic_error("directed AMR streamed patch-cell sender did not cover logical payload");
  }
  if (receive_offset != logical_receive_count) {
    throw std::logic_error("directed AMR streamed patch-cell receiver did not cover logical payload");
  }
#else
  (void)peer_rank;
  (void)logical_send_count;
  (void)local_requests;
  (void)remote_patches;
  (void)remote_requests;
  (void)provider;
  (void)consumer;
  (void)admission;
  (void)count_tag_base;
  (void)payload_tag_base;
  (void)exchange_sequence;
  (void)local_transport_round_limit_bytes;
  (void)diagnostics;
  (void)chunk_governor;
  throw std::runtime_error("directed AMR patch-cell streaming requires MPI support");
#endif
}

[[nodiscard]] std::uint8_t amrBoundaryFaceBit(AmrPatchBoundaryFace face) noexcept {
  return static_cast<std::uint8_t>(face);
}

[[nodiscard]] double amrInterfaceTolerance(
    const AmrPatchPayloadRecord& lhs,
    const AmrPatchPayloadRecord& rhs) noexcept {
  const double scale = std::max({
      1.0,
      std::abs(lhs.origin_x_comoving), std::abs(amrPatchMaxX(lhs)),
      std::abs(lhs.origin_y_comoving), std::abs(amrPatchMaxY(lhs)),
      std::abs(lhs.origin_z_comoving), std::abs(amrPatchMaxZ(lhs)),
      std::abs(rhs.origin_x_comoving), std::abs(amrPatchMaxX(rhs)),
      std::abs(rhs.origin_y_comoving), std::abs(amrPatchMaxY(rhs)),
      std::abs(rhs.origin_z_comoving), std::abs(amrPatchMaxZ(rhs))});
  return 1.0e-10 * scale;
}

[[nodiscard]] bool amrIntervalsShareArea(
    double lhs_min,
    double lhs_max,
    double rhs_min,
    double rhs_max,
    double tolerance) noexcept {
  return std::min(lhs_max, rhs_max) - std::max(lhs_min, rhs_min) > tolerance;
}

[[nodiscard]] std::uint8_t amrPatchBoundaryMaskForPeer(
    const AmrPatchPayloadRecord& local,
    const AmrPatchPayloadRecord& remote,
    const AmrInterfaceTopologyOptions& topology = {}) noexcept {
  if (local.owner_rank == remote.owner_rank) {
    return 0U;
  }
  const std::array<double, 3> local_min{
      local.origin_x_comoving, local.origin_y_comoving, local.origin_z_comoving};
  const std::array<double, 3> local_max{
      amrPatchMaxX(local), amrPatchMaxY(local), amrPatchMaxZ(local)};
  const std::array<double, 3> remote_min{
      remote.origin_x_comoving, remote.origin_y_comoving, remote.origin_z_comoving};
  const std::array<double, 3> remote_max{
      amrPatchMaxX(remote), amrPatchMaxY(remote), amrPatchMaxZ(remote)};
  const double tolerance = amrInterfaceTolerance(local, remote);
  std::uint8_t mask = 0U;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t transverse_a = (axis + 1U) % 3U;
    const std::size_t transverse_b = (axis + 2U) % 3U;
    if (!amrIntervalsShareArea(
            local_min[transverse_a], local_max[transverse_a],
            remote_min[transverse_a], remote_max[transverse_a], tolerance) ||
        !amrIntervalsShareArea(
            local_min[transverse_b], local_max[transverse_b],
            remote_min[transverse_b], remote_max[transverse_b], tolerance)) {
      continue;
    }
    // Minimum-image normal gap on periodic axes (CHUI Part 3): a lower box
    // face at 0 meets a remote upper face at L across the wrap.
    const double length = topologyAxisLength(topology, axis);
    const auto gap = [length](double a, double b) {
      double delta = a - b;
      if (length > 0.0) {
        delta -= length * std::round(delta / length);
      }
      return std::abs(delta);
    };
    if (gap(local_min[axis], remote_max[axis]) <= tolerance) {
      mask |= amrBoundaryFaceBit(axis == 0U ? AmrPatchBoundaryFace::kXLower :
          (axis == 1U ? AmrPatchBoundaryFace::kYLower : AmrPatchBoundaryFace::kZLower));
    }
    if (gap(local_max[axis], remote_min[axis]) <= tolerance) {
      mask |= amrBoundaryFaceBit(axis == 0U ? AmrPatchBoundaryFace::kXUpper :
          (axis == 1U ? AmrPatchBoundaryFace::kYUpper : AmrPatchBoundaryFace::kZUpper));
    }
  }
  return mask;
}

[[nodiscard]] std::size_t amrBoundaryFaceIndex(AmrPatchBoundaryFace face) noexcept {
  switch (face) {
    case AmrPatchBoundaryFace::kXLower: return 0U;
    case AmrPatchBoundaryFace::kXUpper: return 1U;
    case AmrPatchBoundaryFace::kYLower: return 2U;
    case AmrPatchBoundaryFace::kYUpper: return 3U;
    case AmrPatchBoundaryFace::kZLower: return 4U;
    case AmrPatchBoundaryFace::kZUpper: return 5U;
  }
  return 0U;
}

[[nodiscard]] std::size_t amrBoundaryFaceAxis(AmrPatchBoundaryFace face) noexcept {
  switch (face) {
    case AmrPatchBoundaryFace::kXLower:
    case AmrPatchBoundaryFace::kXUpper:
      return 0U;
    case AmrPatchBoundaryFace::kYLower:
    case AmrPatchBoundaryFace::kYUpper:
      return 1U;
    case AmrPatchBoundaryFace::kZLower:
    case AmrPatchBoundaryFace::kZUpper:
      return 2U;
  }
  return 0U;
}

[[nodiscard]] std::uint16_t amrPatchDimension(
    const AmrPatchPayloadRecord& patch,
    std::size_t axis) {
  switch (axis) {
    case 0U: return patch.cell_dim_x;
    case 1U: return patch.cell_dim_y;
    case 2U: return patch.cell_dim_z;
    default: throw std::out_of_range("directed AMR boundary axis out of range");
  }
}

[[nodiscard]] double amrPatchExtent(
    const AmrPatchPayloadRecord& patch,
    std::size_t axis) {
  switch (axis) {
    case 0U: return patch.extent_x_comoving;
    case 1U: return patch.extent_y_comoving;
    case 2U: return patch.extent_z_comoving;
    default: throw std::out_of_range("directed AMR boundary axis out of range");
  }
}

[[nodiscard]] std::uint16_t requiredSourceBoundaryDepth(
    const AmrPatchPayloadRecord& source,
    const AmrPatchPayloadRecord& target,
    std::size_t axis,
    std::uint16_t ghost_layers) {
  const std::uint16_t source_dim = amrPatchDimension(source, axis);
  const std::uint16_t target_dim = amrPatchDimension(target, axis);
  if (source_dim == 0U || target_dim == 0U) {
    throw std::invalid_argument("directed AMR boundary depth requires positive patch dimensions");
  }
  const double source_width = amrPatchExtent(source, axis) / static_cast<double>(source_dim);
  const double target_width = amrPatchExtent(target, axis) / static_cast<double>(target_dim);
  if (!(source_width > 0.0) || !(target_width > 0.0) ||
      !std::isfinite(source_width) || !std::isfinite(target_width)) {
    throw std::invalid_argument("directed AMR boundary depth requires finite positive cell widths");
  }

  const double ratio = target_width / source_width;
  const double ratio_tolerance = 1.0e-10 * std::max(1.0, std::abs(ratio));
  const std::size_t layers = std::max<std::uint16_t>(ghost_layers, 1U);
  // Deep stencil-only ghost layers (CHUI Part 3) sit up to `layers` target
  // widths from the face; deeper source cells beyond this patch belong to the
  // next patch along the normal, which sends its own request.
  const auto clampDeep = [source_dim](double cells) {
    const double rounded_up = std::ceil(cells - 1.0e-9);
    return static_cast<std::uint16_t>(std::clamp<double>(rounded_up, 1.0, static_cast<double>(source_dim)));
  };
  if (ratio <= 1.0 + ratio_tolerance) {
    return layers == 1U ? 1U : clampDeep(static_cast<double>(layers) * ratio);
  }
  const double rounded = std::round(ratio);
  if (std::abs(ratio - rounded) > ratio_tolerance || rounded < 1.0 ||
      rounded > static_cast<double>(std::numeric_limits<std::uint16_t>::max())) {
    throw std::runtime_error(
        "directed AMR fine-to-coarse boundary depth is not an aligned integer cell-width ratio");
  }
  const auto depth = static_cast<std::uint16_t>(rounded);
  if (depth > source_dim) {
    throw std::runtime_error(
        "directed AMR fine-to-coarse boundary depth exceeds source patch dimension");
  }
  return layers == 1U ? depth : clampDeep(static_cast<double>(layers) * rounded);
}

[[nodiscard]] AmrPatchBoundaryCellRequest amrPatchBoundaryRequestForPeer(
    const AmrPatchPayloadRecord& source,
    const AmrPatchPayloadRecord& target,
    const AmrInterfaceTopologyOptions& topology) {
  AmrPatchBoundaryCellRequest request;
  request.patch_id = source.patch_id;
  request.boundary_face_mask = amrPatchBoundaryMaskForPeer(source, target, topology);
  constexpr std::array<AmrPatchBoundaryFace, 6> k_faces{
      AmrPatchBoundaryFace::kXLower,
      AmrPatchBoundaryFace::kXUpper,
      AmrPatchBoundaryFace::kYLower,
      AmrPatchBoundaryFace::kYUpper,
      AmrPatchBoundaryFace::kZLower,
      AmrPatchBoundaryFace::kZUpper};
  for (const AmrPatchBoundaryFace face : k_faces) {
    if ((request.boundary_face_mask & amrBoundaryFaceBit(face)) == 0U) {
      continue;
    }
    request.boundary_face_depths[amrBoundaryFaceIndex(face)] =
        requiredSourceBoundaryDepth(source, target, amrBoundaryFaceAxis(face), topology.ghost_layers);
  }
  return request;
}

void mergeBoundaryRequest(
    AmrPatchBoundaryCellRequest& destination,
    const AmrPatchBoundaryCellRequest& source) {
  if (destination.patch_id == 0U) {
    destination.patch_id = source.patch_id;
  }
  if (destination.patch_id != source.patch_id) {
    throw std::invalid_argument("cannot merge directed AMR boundary requests for different patches");
  }
  destination.boundary_face_mask |= source.boundary_face_mask;
  for (std::size_t face = 0; face < destination.boundary_face_depths.size(); ++face) {
    destination.boundary_face_depths[face] = std::max(
        destination.boundary_face_depths[face], source.boundary_face_depths[face]);
  }
}

[[nodiscard]] std::uint16_t requestFaceDepth(
    const AmrPatchBoundaryCellRequest& request,
    AmrPatchBoundaryFace face) {
  const bool selected = (request.boundary_face_mask & amrBoundaryFaceBit(face)) != 0U;
  const std::uint16_t depth = request.boundary_face_depths[amrBoundaryFaceIndex(face)];
  if (selected && depth == 0U) {
    throw std::invalid_argument("directed AMR selected boundary face has zero source depth");
  }
  if (!selected && depth != 0U) {
    throw std::invalid_argument("directed AMR unselected boundary face has non-zero source depth");
  }
  return selected ? depth : 0U;
}

[[nodiscard]] std::size_t requestedBoundaryCellCount(
    const AmrPatchPayloadRecord& patch,
    const AmrPatchBoundaryCellRequest& request) {
  if (request.patch_id != patch.patch_id || request.boundary_face_mask == 0U) {
    throw std::invalid_argument("directed AMR boundary count request does not match patch metadata");
  }
  const auto selected_positions = [&request](
      std::uint16_t dim,
      AmrPatchBoundaryFace lower,
      AmrPatchBoundaryFace upper) -> std::size_t {
    const std::size_t lower_depth = requestFaceDepth(request, lower);
    const std::size_t upper_depth = requestFaceDepth(request, upper);
    if (lower_depth > dim || upper_depth > dim) {
      throw std::invalid_argument("directed AMR boundary depth exceeds patch dimension");
    }
    return std::min<std::size_t>(
        dim, core::checkedSizeAdd(lower_depth, upper_depth,
                                  "directed AMR selected boundary depth"));
  };
  const std::size_t nx = patch.cell_dim_x;
  const std::size_t ny = patch.cell_dim_y;
  const std::size_t nz = patch.cell_dim_z;
  const std::size_t total = core::checkedSizeProduct3(
      nx, ny, nz, "directed AMR boundary cell count");
  const std::size_t interior_x = nx - selected_positions(
      patch.cell_dim_x, AmrPatchBoundaryFace::kXLower, AmrPatchBoundaryFace::kXUpper);
  const std::size_t interior_y = ny - selected_positions(
      patch.cell_dim_y, AmrPatchBoundaryFace::kYLower, AmrPatchBoundaryFace::kYUpper);
  const std::size_t interior_z = nz - selected_positions(
      patch.cell_dim_z, AmrPatchBoundaryFace::kZLower, AmrPatchBoundaryFace::kZUpper);
  const std::size_t unselected = core::checkedSizeProduct3(
      interior_x, interior_y, interior_z, "directed AMR unselected interior cell count");
  return total - unselected;
}

[[nodiscard]] AmrPatchBoundaryCellRequest oneLayerBoundaryRequest(
    const AmrPatchPayloadRecord& patch,
    std::uint8_t mask) {
  AmrPatchBoundaryCellRequest request;
  request.patch_id = patch.patch_id;
  request.boundary_face_mask = mask;
  constexpr std::array<AmrPatchBoundaryFace, 6> k_faces{
      AmrPatchBoundaryFace::kXLower,
      AmrPatchBoundaryFace::kXUpper,
      AmrPatchBoundaryFace::kYLower,
      AmrPatchBoundaryFace::kYUpper,
      AmrPatchBoundaryFace::kZLower,
      AmrPatchBoundaryFace::kZUpper};
  for (const AmrPatchBoundaryFace face : k_faces) {
    if ((mask & amrBoundaryFaceBit(face)) != 0U) {
      request.boundary_face_depths[amrBoundaryFaceIndex(face)] = 1U;
    }
  }
  return request;
}

[[nodiscard]] bool patchCellOffsetMatchesBoundaryRequest(
    const AmrPatchPayloadRecord& patch,
    std::uint32_t offset,
    const AmrPatchBoundaryCellRequest& request) {
  if (request.patch_id != patch.patch_id || patch.cell_dim_x == 0U ||
      patch.cell_dim_y == 0U || patch.cell_dim_z == 0U || offset >= patch.cell_count) {
    return false;
  }
  const std::size_t nx = patch.cell_dim_x;
  const std::size_t ny = patch.cell_dim_y;
  const std::size_t nz = patch.cell_dim_z;
  const std::size_t plane = core::checkedSizeMultiply(nx, ny, "directed AMR patch plane");
  const std::size_t i = offset % nx;
  const std::size_t j = (offset / nx) % ny;
  const std::size_t k = offset / plane;
  const std::size_t x_lower = requestFaceDepth(request, AmrPatchBoundaryFace::kXLower);
  const std::size_t x_upper = requestFaceDepth(request, AmrPatchBoundaryFace::kXUpper);
  const std::size_t y_lower = requestFaceDepth(request, AmrPatchBoundaryFace::kYLower);
  const std::size_t y_upper = requestFaceDepth(request, AmrPatchBoundaryFace::kYUpper);
  const std::size_t z_lower = requestFaceDepth(request, AmrPatchBoundaryFace::kZLower);
  const std::size_t z_upper = requestFaceDepth(request, AmrPatchBoundaryFace::kZUpper);
  if (x_lower > nx || x_upper > nx || y_lower > ny || y_upper > ny ||
      z_lower > nz || z_upper > nz) {
    throw std::invalid_argument("directed AMR boundary depth exceeds patch dimension");
  }
  return (x_lower != 0U && i < x_lower) ||
      (x_upper != 0U && i >= nx - x_upper) ||
      (y_lower != 0U && j < y_lower) ||
      (y_upper != 0U && j >= ny - y_upper) ||
      (z_lower != 0U && k < z_lower) ||
      (z_upper != 0U && k >= nz - z_upper);
}

[[nodiscard]] std::uint64_t checkedU64AddLocal(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 addition overflow");
  }
  return lhs + rhs;
}

[[nodiscard]] std::uint64_t checkedRecordTrafficBytes(
    std::size_t sent_records,
    std::size_t received_records,
    std::size_t record_bytes,
    std::string_view context) {
  const std::size_t record_count = core::checkedSizeAdd(
      sent_records, received_records, std::string(context) + " record count");
  const std::size_t bytes = core::checkedSizeMultiply(
      record_count, record_bytes, std::string(context) + " byte count");
  return core::checkedIntegralNarrow<std::uint64_t>(
      bytes, std::string(context) + " uint64 byte count");
}

[[nodiscard]] std::unordered_map<std::uint64_t, AmrPatchBoundaryCellRequest>
requestByPatchId(std::span<const AmrPatchBoundaryCellRequest> requests) {
  std::unordered_map<std::uint64_t, AmrPatchBoundaryCellRequest> result;
  result.reserve(requests.size());
  for (const AmrPatchBoundaryCellRequest& request : requests) {
    if (request.patch_id == 0U || request.boundary_face_mask == 0U) {
      throw std::invalid_argument("directed AMR boundary request is empty or malformed");
    }
    auto [it, inserted] = result.emplace(request.patch_id, request);
    if (!inserted) {
      mergeBoundaryRequest(it->second, request);
    }
  }
  return result;
}

}  // namespace

std::vector<AmrPatchBoundaryCellRequest> planDirectedAmrPatchBoundaryCellRequests(
    std::span<const AmrPatchPayloadRecord> local_patch_records,
    std::span<const AmrPatchPayloadRecord> remote_patch_records,
    const AmrInterfaceTopologyOptions& topology) {
  std::vector<AmrPatchBoundaryCellRequest> requests;
  requests.reserve(local_patch_records.size());
  for (const AmrPatchPayloadRecord& local_record : local_patch_records) {
    validateAmrPatchPayloadRecord(local_record);
    AmrPatchBoundaryCellRequest combined;
    combined.patch_id = local_record.patch_id;
    for (const AmrPatchPayloadRecord& remote_record : remote_patch_records) {
      validateAmrPatchPayloadRecord(remote_record);
      const AmrPatchBoundaryCellRequest one =
          amrPatchBoundaryRequestForPeer(local_record, remote_record, topology);
      if (one.boundary_face_mask != 0U) {
        mergeBoundaryRequest(combined, one);
      }
    }
    if (combined.boundary_face_mask != 0U) {
      requests.push_back(combined);
    }
  }
  return requests;
}

std::size_t directedAmrPatchBoundaryCellCount(
    const AmrPatchPayloadRecord& patch,
    const AmrPatchBoundaryCellRequest& request) {
  validateAmrPatchPayloadRecord(patch);
  if (request.boundary_face_mask == 0U) {
    return 0U;
  }
  return requestedBoundaryCellCount(patch, request);
}

std::size_t directedAmrPatchBoundaryCellCount(
    const AmrPatchPayloadRecord& patch,
    std::uint8_t boundary_face_mask) {
  validateAmrPatchPayloadRecord(patch);
  if (boundary_face_mask == 0U) {
    return 0U;
  }
  return requestedBoundaryCellCount(
      patch, oneLayerBoundaryRequest(patch, boundary_face_mask));
}

DirectedAmrPatchPayloadExchange executeBlockingDirectedAmrPatchPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_patch_records,
    const DirectedAmrPatchCellPayloadProvider& cell_payload_provider,
    const DirectedAmrPatchCellPayloadConsumer& cell_payload_consumer,
    const DirectedAmrPatchCellAdmission& cell_payload_admission,
    std::size_t transport_round_limit_bytes,
    std::uint64_t exchange_sequence,
    const AmrInterfaceTopologyOptions& topology,
    core::MemoryGovernor* memory_governor) {
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)exchange_sequence;
  (void)topology;
  (void)transport_round_limit_bytes;
  (void)cell_payload_admission;
#endif
  const int world_rank = mpi_context.worldRank();
  const bool distributed = mpi_context.isEnabled() && mpi_context.worldSize() > 1;
  // E0 [local, guarded]: local patch index (admitted), record validation,
  // callback presence. Its failure is folded into the first collective vote.
  core::MemoryReservation local_index_lease;
  std::unordered_map<std::uint64_t, AmrPatchPayloadRecord> local_patch_by_id;
  std::exception_ptr local_failure;
  try {
    if (memory_governor != nullptr && !local_patch_records.empty()) {
      local_index_lease = memory_governor->reserve(
          core::MemoryClass::kCommunication,
          core::checkedSizeMultiply(local_patch_records.size(),
                                    sizeof(AmrPatchPayloadRecord) + k_directed_amr_hash_node_overhead_bytes,
                                    "directed AMR local patch index bytes"),
          "parallel.directed_amr.local_patch_index");
      local_index_lease.commit();
    }
    local_patch_by_id.reserve(local_patch_records.size());
    for (const AmrPatchPayloadRecord& record : local_patch_records) {
      validateAmrPatchPayloadRecord(record);
      if (record.owner_rank != world_rank) {
        throw std::invalid_argument("directed AMR patch exchange received non-local authoritative patch metadata");
      }
      const auto [it, inserted] = local_patch_by_id.emplace(record.patch_id, record);
      if (!inserted) {
        throw std::invalid_argument("directed AMR patch exchange found duplicate local patch metadata");
      }
    }
    if ((!cell_payload_provider || !cell_payload_consumer) && !local_patch_records.empty()) {
      throw std::invalid_argument(
          "directed AMR patch exchange requires boundary-cell producer and consumer callbacks");
    }
  } catch (...) {
    local_failure = std::current_exception();
  }

  DirectedAmrPatchPayloadExchange staged;
  if (!distributed) {
    if (local_failure != nullptr) std::rethrow_exception(local_failure);
    return staged;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_patch_count_tag_base = 8910;
  constexpr int k_patch_payload_tag_base = 9910;
  constexpr int k_cell_count_tag_base = 10910;
  constexpr int k_cell_payload_tag_base = 11910;
  constexpr int k_pair_health_tag_base = 7910;
  constexpr int k_pair_plan_ready_tag_base = 7920;
  // E0 vote + E1 Allgather + E2 peer graph (two votes) inside discovery.
  core::MemoryReservation control_lease;
  const AmrCandidatePeerSet candidates = discoverCandidateAmrPeers(
      mpi_context, local_patch_records, &staged.diagnostics, topology, memory_governor, control_lease,
      local_failure);
  // E2 (second half) [local, guarded]: admit and reserve the staged result at
  // its upper bound (sum of candidate peers' patch counts) so pair epochs
  // only append within capacity.
  local_failure = nullptr;
  try {
    const std::size_t bound = core::checkedIntegralNarrow<std::size_t>(
        candidates.remote_patch_upper_bound, "directed AMR staged descriptor bound");
    if (memory_governor != nullptr && bound != 0U) {
      staged.patch_payload_lease = memory_governor->reserve(
          core::MemoryClass::kCommunication,
          core::checkedSizeMultiply(bound, sizeof(AmrPatchPayloadRecord), "directed AMR staged descriptor bytes"),
          "parallel.directed_amr.remote_patch_descriptors");
      staged.patch_payload_lease.commit();
    }
    staged.patch_payloads_received.reserve(bound);
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(local_failure, "directed AMR staged descriptor admission");

  // E3: pair epochs. `exchange_failure` is the first local failure of this
  // rank; once set, every remaining pair epoch is opened with health 0 and
  // skipped symmetrically by both ranks.
  std::exception_ptr exchange_failure;
  for (const int peer_rank : candidates.peers) {
    const bool peer_healthy = exchangePairHealth(
        mpi_context, peer_rank, exchange_failure == nullptr, k_pair_health_tag_base, exchange_sequence);
    if (exchange_failure != nullptr) continue;
    if (!peer_healthy) {
      exchange_failure = std::make_exception_ptr(std::runtime_error(
          "directed AMR peer rank " + std::to_string(peer_rank) + " reported an earlier local failure"));
      continue;
    }
    try {
      // (a) descriptors: internal guarded receive admission + readiness.
      core::MemoryReservation remote_descriptor_lease;
      std::vector<AmrPatchPayloadRecord> remote_patches = exchangePodRecordsWithPeer(
          mpi_context,
          peer_rank,
          local_patch_records,
          k_patch_count_tag_base,
          k_patch_payload_tag_base,
          exchange_sequence,
          "executeBlockingDirectedAmrPatchPayloadExchange/patch",
          nullptr,
          memory_governor,
          &remote_descriptor_lease);
      // (b) [local, guarded] validation, planning, indices, staging; then a
      // pairwise readiness message before any further message of the pair.
      core::MemoryReservation plan_lease;
      std::vector<AmrPatchBoundaryCellRequest> local_requests;
      std::vector<AmrPatchBoundaryCellRequest> remote_requests;
      std::unordered_map<std::uint64_t, AmrPatchBoundaryCellRequest> local_request_by_id;
      std::unordered_map<std::uint64_t, AmrPatchBoundaryCellRequest> remote_request_by_id;
      std::unordered_map<std::uint64_t, AmrPatchPayloadRecord> remote_patch_by_id;
      std::unordered_map<std::uint64_t, std::size_t> observed_local_count;
      std::unordered_map<std::uint64_t, std::size_t> observed_remote_count;
      std::uint64_t logical_send_count = 0U;
      std::exception_ptr plan_failure;
      try {
        if (remote_patches.size() > candidates.max_peer_patch_count) {
          throw std::runtime_error("directed AMR peer sent more descriptors than its gathered patch count");
        }
        for (const AmrPatchPayloadRecord& record : remote_patches) {
          validateAmrPatchPayloadRecord(record);
          if (record.owner_rank != peer_rank) {
            throw std::runtime_error(
                "directed AMR patch exchange received metadata from a rank that is not the owner");
          }
        }
        const std::size_t n_local = local_patch_records.size();
        const std::size_t n_remote = remote_patches.size();
        if (memory_governor != nullptr && (n_local + n_remote) != 0U) {
          // Requests (<= n_local, <= n_remote), three request/patch indices
          // and two coverage counters, each with a conservative hash-node
          // overhead estimate; all live for this pair epoch only.
          const std::size_t entry_bytes = sizeof(AmrPatchBoundaryCellRequest) + sizeof(AmrPatchPayloadRecord) +
              3U * k_directed_amr_hash_node_overhead_bytes + 2U * sizeof(std::size_t);
          plan_lease = memory_governor->reserve(
              core::MemoryClass::kCommunication,
              core::checkedSizeMultiply(core::checkedSizeAdd(n_local, n_remote, "directed AMR pair entries"),
                                        entry_bytes, "directed AMR pair planning bytes"),
              "parallel.directed_amr.pair_planning");
          plan_lease.commit();
        }
        local_requests = planDirectedAmrPatchBoundaryCellRequests(local_patch_records, remote_patches, topology);
        remote_requests = planDirectedAmrPatchBoundaryCellRequests(remote_patches, local_patch_records, topology);
        if (!local_requests.empty() || !remote_requests.empty()) {
          local_request_by_id = requestByPatchId(local_requests);
          remote_request_by_id = requestByPatchId(remote_requests);
          remote_patch_by_id.reserve(n_remote);
          for (const AmrPatchPayloadRecord& patch : remote_patches) {
            remote_patch_by_id.emplace(patch.patch_id, patch);
          }
          observed_local_count.reserve(local_requests.size());
          observed_remote_count.reserve(remote_requests.size());
          for (const AmrPatchBoundaryCellRequest& request : local_requests) {
            const auto patch_it = local_patch_by_id.find(request.patch_id);
            if (patch_it == local_patch_by_id.end()) {
              throw std::logic_error("directed AMR local request names an unknown local patch");
            }
            logical_send_count = checkedU64AddLocal(
                logical_send_count,
                static_cast<std::uint64_t>(requestedBoundaryCellCount(patch_it->second, request)),
                "directed AMR logical send boundary-cell count");
          }
          for (const AmrPatchBoundaryCellRequest& request : remote_requests) {
            if (!remote_patch_by_id.contains(request.patch_id)) {
              throw std::logic_error("directed AMR remote request names an unknown remote patch");
            }
          }
          std::uint64_t peer_interface_count = 0U;
          for (const AmrPatchPayloadRecord& local_record : local_patch_records) {
            for (const AmrPatchPayloadRecord& remote_record : remote_patches) {
              if (amrPatchBoundaryMaskForPeer(local_record, remote_record, topology) != 0U) {
                peer_interface_count = checkedU64AddLocal(
                    peer_interface_count, 1U, "directed AMR peer interface count");
              }
            }
          }
          staged.diagnostics.remote_interface_count = checkedU64AddLocal(
              staged.diagnostics.remote_interface_count, peer_interface_count,
              "directed AMR accumulated remote interface count");
          for (const AmrPatchPayloadRecord& remote_patch : remote_patches) {
            if (remote_request_by_id.contains(remote_patch.patch_id)) {
              if (staged.patch_payloads_received.size() == staged.patch_payloads_received.capacity()) {
                throw std::logic_error("directed AMR staged descriptors exceed their admitted capacity");
              }
              staged.patch_payloads_received.push_back(remote_patch);  // within reserved capacity
            }
          }
        }
        staged.diagnostics.directed_patch_descriptor_records_sent = checkedU64AddLocal(
            staged.diagnostics.directed_patch_descriptor_records_sent,
            static_cast<std::uint64_t>(n_local), "directed AMR descriptor records sent");
        staged.diagnostics.directed_patch_descriptor_records_received = checkedU64AddLocal(
            staged.diagnostics.directed_patch_descriptor_records_received,
            static_cast<std::uint64_t>(n_remote), "directed AMR descriptor records received");
        staged.diagnostics.patch_descriptor_bytes = checkedU64AddLocal(
            staged.diagnostics.patch_descriptor_bytes,
            checkedRecordTrafficBytes(n_local, n_remote, sizeof(AmrPatchPayloadRecord),
                                      "directed AMR patch descriptor traffic"),
            "directed AMR accumulated patch descriptor traffic");
      } catch (...) {
        plan_failure = std::current_exception();
      }
      const bool peer_planned = exchangePairHealth(
          mpi_context, peer_rank, plan_failure == nullptr, k_pair_plan_ready_tag_base, exchange_sequence);
      if (plan_failure != nullptr) std::rethrow_exception(plan_failure);
      if (!peer_planned) {
        throw std::runtime_error("directed AMR peer rank " + std::to_string(peer_rank) +
                                 " rejected its descriptor validation/planning");
      }
      // Symmetric decision: this rank's local requests are the peer's remote
      // requests (identical function, identical argument order and bytes).
      if (local_requests.empty() && remote_requests.empty()) {
        continue;
      }
      ++staged.diagnostics.neighbor_peer_count;

      // (c) streamed boundary cells: callbacks run inside the stream's own
      // guarded phases, each followed by a pairwise readiness message.
      const DirectedAmrPatchCellPayloadProvider validated_provider =
          [&](std::span<const AmrPatchBoundaryCellRequest> requests,
              std::uint64_t first_record,
              std::size_t max_records,
              std::vector<AmrPatchCellPayloadRecord>& output) {
            cell_payload_provider(requests, first_record, max_records, output);
            for (const AmrPatchCellPayloadRecord& record : output) {
              validateAmrPatchCellPayloadRecord(record);
              if (record.owner_rank != world_rank) {
                throw std::invalid_argument(
                    "directed AMR boundary payload producer returned a non-local cell");
              }
              const auto patch_it = local_patch_by_id.find(record.patch_id);
              const auto request_it = local_request_by_id.find(record.patch_id);
              if (patch_it == local_patch_by_id.end() || request_it == local_request_by_id.end() ||
                  !patchCellOffsetMatchesBoundaryRequest(
                      patch_it->second, record.local_cell_offset, request_it->second)) {
                throw std::invalid_argument(
                    "directed AMR boundary payload producer returned a non-interface cell");
              }
              ++observed_local_count[record.patch_id];
            }
          };
      const DirectedAmrPatchCellPayloadConsumer validated_consumer =
          [&](int source_rank, std::span<const AmrPatchCellPayloadRecord> records) {
            if (source_rank != peer_rank) {
              throw std::runtime_error("directed AMR patch-cell consumer received wrong peer rank");
            }
            for (const AmrPatchCellPayloadRecord& record : records) {
              validateAmrPatchCellPayloadRecord(record);
              if (record.owner_rank != peer_rank) {
                throw std::runtime_error(
                    "directed AMR patch-cell stream received payload from a rank that is not the owner");
              }
              const auto patch_it = remote_patch_by_id.find(record.patch_id);
              const auto request_it = remote_request_by_id.find(record.patch_id);
              if (patch_it == remote_patch_by_id.end() || request_it == remote_request_by_id.end() ||
                  !patchCellOffsetMatchesBoundaryRequest(
                      patch_it->second, record.local_cell_offset, request_it->second)) {
                throw std::runtime_error(
                    "directed AMR patch-cell stream received a non-interface cell payload");
              }
              ++observed_remote_count[record.patch_id];
            }
            cell_payload_consumer(source_rank, records);
          };

      exchangeAmrPatchCellRecordStreamWithPeer(
          mpi_context,
          peer_rank,
          logical_send_count,
          local_requests,
          remote_patches,
          remote_requests,
          validated_provider,
          validated_consumer,
          cell_payload_admission,
          k_cell_count_tag_base,
          k_cell_payload_tag_base,
          exchange_sequence,
          transport_round_limit_bytes,
          staged.diagnostics,
          transport_round_limit_bytes == 0U ? memory_governor : nullptr);

      // (d) [local] coverage checks after the last message of this pair
      // epoch: a failure only marks this rank failed (announced to its
      // remaining peers, agreed at E4).
      for (const AmrPatchBoundaryCellRequest& request : local_requests) {
        const AmrPatchPayloadRecord& patch = local_patch_by_id.at(request.patch_id);
        if (observed_local_count[request.patch_id] != requestedBoundaryCellCount(patch, request)) {
          throw std::runtime_error(
              "directed AMR streamed producer did not cover requested local boundary cells");
        }
      }
      std::uint64_t logical_receive_count = 0U;
      for (const AmrPatchBoundaryCellRequest& request : remote_requests) {
        const AmrPatchPayloadRecord& patch = remote_patch_by_id.at(request.patch_id);
        const std::size_t expected = requestedBoundaryCellCount(patch, request);
        if (observed_remote_count[request.patch_id] != expected) {
          throw std::runtime_error("directed AMR streamed remote boundary coverage mismatch");
        }
        logical_receive_count = checkedU64AddLocal(
            logical_receive_count, static_cast<std::uint64_t>(expected),
            "directed AMR logical receive boundary-cell count");
      }
      staged.diagnostics.directed_patch_cell_records_sent = checkedU64AddLocal(
          staged.diagnostics.directed_patch_cell_records_sent, logical_send_count,
          "directed AMR accumulated sent boundary-cell records");
      staged.diagnostics.directed_patch_cell_records_received = checkedU64AddLocal(
          staged.diagnostics.directed_patch_cell_records_received, logical_receive_count,
          "directed AMR accumulated received boundary-cell records");
      const std::uint64_t traffic_records = checkedU64AddLocal(
          logical_send_count, logical_receive_count, "directed AMR streamed traffic record count");
      staged.diagnostics.patch_cell_payload_bytes = checkedU64AddLocal(
          staged.diagnostics.patch_cell_payload_bytes,
          checkedRecordTrafficBytes(
              core::checkedIntegralNarrow<std::size_t>(traffic_records, "directed AMR streamed traffic records"),
              0U, sizeof(AmrPatchCellPayloadRecord), "directed AMR streamed traffic bytes"),
          "directed AMR accumulated streamed patch-cell traffic");
    } catch (...) {
      // Every throw inside the epoch either happened symmetrically at a
      // pairwise readiness message (the peer throws at the same message) or
      // after the epoch's last message. Either way the pair is finished.
      exchange_failure = std::current_exception();
    }
  }
  // E4 [collective]: publish only on unanimous success.
  mpi_context.rethrowCollectivePreparationFailure(exchange_failure, "directed AMR patch exchange");
  staged.diagnostics.remote_patch_ghost_count =
      static_cast<std::uint64_t>(staged.patch_payloads_received.size());
  return staged;
#else
  throw std::runtime_error("directed AMR patch payload exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<AmrPatchPayloadRecord> executeBlockingAmrPatchPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  for (const AmrPatchPayloadRecord& record : local_records) {
    validateAmrPatchPayloadRecord(record);
    if (record.owner_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("AMR patch payload compatibility exchange received non-local authoritative patch metadata");
    }
  }
  return std::vector<AmrPatchPayloadRecord>(local_records.begin(), local_records.end());
}

std::vector<AmrPatchCellPayloadRecord> executeBlockingAmrPatchCellPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchCellPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
  for (const AmrPatchCellPayloadRecord& record : local_records) {
    validateAmrPatchCellPayloadRecord(record);
    if (record.owner_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("directed AMR patch-cell compatibility exchange received non-local authoritative cell payload");
    }
  }
  (void)exchange_sequence;
  return std::vector<AmrPatchCellPayloadRecord>(local_records.begin(), local_records.end());
}

AmrFluxExchangeStagingPlan planAmrFluxExchangeStaging(
    std::size_t total_inbound_capacity,
    std::size_t max_peer_send_count,
    std::size_t max_peer_receive_count,
    std::size_t world_size) {
  // Peak physical staging for the sequential blocking protocol:
  //   retained inbound result
  //   + one active peer outbound staging buffer (max over peers)
  //   + one active peer receive staging buffer (max over peers)
  //   + O(world_size) count metadata (send/recv u64 + active-peer int).
  // Peer buffers are mutually exclusive across loop iterations, so only the
  // single largest buffer in each direction is charged.
  const std::size_t inbound_bytes = core::checkedSizeMultiply(
      total_inbound_capacity, sizeof(AmrFluxRegisterPayloadRecord),
      "AMR flux-register inbound staging bytes");
  const std::size_t outbound_bytes = core::checkedSizeMultiply(
      max_peer_send_count, sizeof(AmrFluxRegisterPayloadRecord),
      "AMR flux-register peak outbound staging bytes");
  const std::size_t receive_bytes = core::checkedSizeMultiply(
      max_peer_receive_count, sizeof(AmrFluxRegisterPayloadRecord),
      "AMR flux-register peak receive staging bytes");
  const std::size_t send_counts_bytes = core::checkedSizeMultiply(
      world_size, sizeof(std::uint64_t), "AMR flux-register send count metadata");
  const std::size_t recv_counts_bytes = core::checkedSizeMultiply(
      world_size, sizeof(std::uint64_t), "AMR flux-register recv count metadata");
  const std::size_t peer_list_bytes = core::checkedSizeMultiply(
      world_size, sizeof(int), "AMR flux-register active-peer metadata");
  std::size_t peak = core::checkedSizeAdd(
      inbound_bytes, outbound_bytes, "AMR flux-register staging peak");
  peak = core::checkedSizeAdd(peak, receive_bytes, "AMR flux-register staging peak");
  peak = core::checkedSizeAdd(peak, send_counts_bytes, "AMR flux-register staging peak");
  peak = core::checkedSizeAdd(peak, recv_counts_bytes, "AMR flux-register staging peak");
  peak = core::checkedSizeAdd(peak, peer_list_bytes, "AMR flux-register staging peak");
  AmrFluxExchangeStagingPlan plan;
  plan.total_inbound_capacity = total_inbound_capacity;
  plan.max_peer_send_count = max_peer_send_count;
  plan.max_peer_receive_count = max_peer_receive_count;
  plan.peak_reservation_bytes =
      core::checkedIntegralNarrow<std::uint64_t>(peak, "AMR flux-register staging peak width");
  return plan;
}

std::uint64_t amrFluxExchangeStagingPeakBytes(const AmrFluxExchangeStagingPlan& plan) {
  return plan.peak_reservation_bytes;
}

AmrFluxRegisterExchangeResult executeGovernedAmrFluxRegisterPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrFluxRegisterPayloadRecord> local_records,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* memory_governor) {
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)exchange_sequence;
#endif
  const int world_rank = mpi_context.worldRank();
  const int world_size = mpi_context.worldSize();
  const bool distributed = mpi_context.isEnabled() && world_size > 1;
  AmrFluxRegisterExchangeResult result;
  // F0 [local, guarded]: record validation + count control plane.
  std::vector<std::uint64_t> send_record_counts;
  std::vector<std::uint64_t> recv_record_counts;
  std::exception_ptr local_preparation_failure;
  try {
    for (const AmrFluxRegisterPayloadRecord& record : local_records) {
      validateAmrFluxRegisterPayloadRecord(record);
      if (record.source_rank != world_rank) {
        throw std::invalid_argument("AMR flux-register payload source rank does not match MPI context");
      }
      if (record.owner_rank < 0 || record.owner_rank >= world_size) {
        throw std::invalid_argument("AMR flux-register payload owner rank is outside MPI world");
      }
    }
    if (distributed) {
      send_record_counts.resize(static_cast<std::size_t>(world_size), 0U);
      recv_record_counts.resize(static_cast<std::size_t>(world_size), 0U);
      for (const AmrFluxRegisterPayloadRecord& record : local_records) {
        if (record.owner_rank == world_rank) {
          continue;
        }
        auto& count = send_record_counts[static_cast<std::size_t>(record.owner_rank)];
        if (count == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error("AMR flux-register owner-target count overflows uint64_t");
        }
        ++count;
      }
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  if (!distributed) {
    if (local_preparation_failure != nullptr) std::rethrow_exception(local_preparation_failure);
    if (memory_governor != nullptr && !local_records.empty()) {
      result.lease = memory_governor->reserve(
          core::MemoryClass::kCommunication,
          core::checkedSizeMultiply(local_records.size(), sizeof(AmrFluxRegisterPayloadRecord),
                                    "AMR flux-register serial result bytes"),
          "amr.flux_exchange.inbound_records");
      result.lease.commit();
    }
    result.records.assign(local_records.begin(), local_records.end());
    return result;
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_flux_count_tag_base = 12910;
  constexpr int k_flux_payload_tag_base = 13910;
  constexpr int k_flux_pair_health_tag_base = 7930;
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure, "AMR flux-register owner-target count preparation");

  if (MPI_Alltoall(
          send_record_counts.data(), 1, MPI_UINT64_T,
          recv_record_counts.data(), 1, MPI_UINT64_T,
          MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("AMR flux-register owner-target count exchange failed");
  }

  // F2 [local, guarded]: active peers (symmetric: p is active for q iff q
  // sends to or receives from p), inbound admission and capacity.
  std::vector<int> active_peers;
  std::size_t max_peer_send_count = 0U;
  local_preparation_failure = nullptr;
  try {
    std::size_t total_inbound_records = 0U;
    active_peers.reserve(static_cast<std::size_t>(world_size));
    for (int peer_rank = 0; peer_rank < world_size; ++peer_rank) {
      if (peer_rank == world_rank) {
        continue;
      }
      const std::size_t peer = static_cast<std::size_t>(peer_rank);
      if (send_record_counts[peer] != 0U || recv_record_counts[peer] != 0U) {
        active_peers.push_back(peer_rank);
      }
      total_inbound_records = core::checkedSizeAdd(
          total_inbound_records,
          core::checkedIntegralNarrow<std::size_t>(
              recv_record_counts[peer], "AMR flux-register inbound owner count"),
          "AMR flux-register inbound record count");
      max_peer_send_count = std::max(
          max_peer_send_count,
          core::checkedIntegralNarrow<std::size_t>(send_record_counts[peer], "AMR flux-register peer send count"));
    }
    const std::size_t local_owner_count = static_cast<std::size_t>(std::count_if(
        local_records.begin(), local_records.end(),
        [world_rank](const AmrFluxRegisterPayloadRecord& record) {
          return record.owner_rank == world_rank;
        }));
    const std::size_t total_inbound_capacity = core::checkedSizeAdd(
        total_inbound_records, local_owner_count,
        "AMR flux-register inbound reserve");
    // Retained result + one active peer's outbound staging (the per-peer
    // receive buffer is admitted by exchangePodRecordsWithPeer).
    const std::uint64_t inbound_bytes = core::checkedSizeMultiply(
        core::checkedSizeAdd(total_inbound_capacity, max_peer_send_count, "AMR flux-register staged records"),
        sizeof(AmrFluxRegisterPayloadRecord),
        "AMR flux-register inbound bytes");
    if (memory_governor != nullptr && inbound_bytes > 0U) {
      result.lease = memory_governor->reserve(
          core::MemoryClass::kCommunication, inbound_bytes,
          "amr.flux_exchange.inbound_records");
      result.lease.commit();
    }
    result.records.reserve(total_inbound_capacity);
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure, "AMR flux-register active-peer preparation");

  // F3: pair epochs in ascending peer rank with the pair-health handshake;
  // F4: one collective vote before the result is published.
  std::exception_ptr exchange_failure;
  for (const int peer_rank : active_peers) {
    const bool peer_healthy = exchangePairHealth(
        mpi_context, peer_rank, exchange_failure == nullptr, k_flux_pair_health_tag_base, exchange_sequence);
    if (exchange_failure != nullptr) continue;
    if (!peer_healthy) {
      exchange_failure = std::make_exception_ptr(std::runtime_error(
          "AMR flux-register peer rank " + std::to_string(peer_rank) + " reported an earlier local failure"));
      continue;
    }
    try {
      const std::size_t peer = static_cast<std::size_t>(peer_rank);
      // Outbound staging precedes the pair's first count message, so its
      // failure is agreed by a dedicated pairwise readiness message.
      std::vector<AmrFluxRegisterPayloadRecord> outbound_to_peer;
      std::exception_ptr stage_failure;
      try {
        outbound_to_peer.reserve(static_cast<std::size_t>(send_record_counts[peer]));
        for (const AmrFluxRegisterPayloadRecord& record : local_records) {
          if (record.owner_rank == peer_rank) {
            outbound_to_peer.push_back(record);
          }
        }
        if (outbound_to_peer.size() != send_record_counts[peer]) {
          throw std::runtime_error("AMR flux-register owner-target payload count changed after count exchange");
        }
      } catch (...) {
        stage_failure = std::current_exception();
      }
      const bool peer_staged = exchangePairHealth(
          mpi_context, peer_rank, stage_failure == nullptr, k_flux_pair_health_tag_base + 1, exchange_sequence);
      if (stage_failure != nullptr) std::rethrow_exception(stage_failure);
      if (!peer_staged) {
        throw std::runtime_error("AMR flux-register peer rank " + std::to_string(peer_rank) +
                                 " failed to stage its outbound records");
      }
      core::MemoryReservation receive_lease;
      std::vector<AmrFluxRegisterPayloadRecord> received_from_peer = exchangePodRecordsWithPeer(
          mpi_context,
          peer_rank,
          std::span<const AmrFluxRegisterPayloadRecord>(outbound_to_peer),
          k_flux_count_tag_base,
          k_flux_payload_tag_base,
          exchange_sequence,
          "executeBlockingAmrFluxRegisterPayloadExchange",
          nullptr,
          memory_governor,
          &receive_lease);
      // After the epoch's last message: local checks only mark this rank.
      if (received_from_peer.size() != recv_record_counts[peer]) {
        throw std::runtime_error("AMR flux-register owner-target receive count changed after count exchange");
      }
      for (const AmrFluxRegisterPayloadRecord& record : received_from_peer) {
        validateAmrFluxRegisterPayloadRecord(record);
        if (record.owner_rank != world_rank || record.source_rank != peer_rank) {
          throw std::runtime_error("AMR flux-register directed exchange returned stale source/owner rank metadata");
        }
      }
      // Within the reserved capacity (sum of the agreed receive counts).
      result.records.insert(result.records.end(), received_from_peer.begin(), received_from_peer.end());
    } catch (...) {
      exchange_failure = std::current_exception();
    }
  }
  try {
    if (exchange_failure == nullptr) {
      for (const AmrFluxRegisterPayloadRecord& record : local_records) {
        if (record.owner_rank == world_rank) {
          result.records.push_back(record);  // within reserved capacity
        }
      }
    }
  } catch (...) {
    exchange_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(exchange_failure, "AMR flux-register exchange");
  return result;
#else
  (void)memory_governor;
  throw std::runtime_error("AMR flux-register payload exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<AmrFluxRegisterPayloadRecord> executeBlockingAmrFluxRegisterPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrFluxRegisterPayloadRecord> local_records,
    std::uint64_t exchange_sequence,
    core::MemoryGovernor* memory_governor) {
  AmrFluxRegisterExchangeResult governed = executeGovernedAmrFluxRegisterPayloadExchange(
      mpi_context, local_records, exchange_sequence, memory_governor);
  return std::move(governed.records);  // storage moves out; lease released here (compatibility only)
}

}  // namespace cosmosim::parallel
