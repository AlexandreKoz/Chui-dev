#pragma once

// CHUI Part 3 Checkpoint B — distributed composition transport payloads.
//
// Two owner-targeted record streams complete the composition sidecar on the
// distributed production AMR path (Part 2 refused it):
//
//  * donor fractions: for every remote ghost source cell that is a donor of a
//    composition face transfer, its owner sends the lane mass fractions
//    (lane mass / gas mass, gas mass = elements + dust) at the step-start
//    state; the receiver advects remote inflow with exactly the owner's
//    composition (request/response by stable gas_cell_id);
//  * corrections: a composition lane delta computed on a rank for a cell
//    owned by another rank (fine-side mirrored transfer onto a remote coarse
//    cell at a coarse/fine interface) is delivered to that owner and applied
//    inside its composition transaction.
//
// Records are fixed-size PODs (k_max_composition_lanes) so the exchange is a
// single byte-level MPI_Alltoallv; counts are exchanged first (Alltoall), all
// ranks always participate (no deadlock on empty sets), and the receive
// order is deterministic (ascending source rank, then send order).
// Without MPI (or one rank) records addressed to the local rank are returned
// and any other destination is refused.

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::parallel {

inline constexpr std::size_t k_max_composition_lanes = 32U;

struct CompositionDonorRequestRecord {
  std::uint64_t gas_cell_id = 0U;
  int requester_rank = 0;
  int owner_rank = 0;
};

struct CompositionDonorFractionRecord {
  std::uint64_t gas_cell_id = 0U;
  int owner_rank = 0;
  int requester_rank = 0;
  std::uint32_t lane_count = 0U;
  std::uint32_t reserved = 0U;
  std::array<double, k_max_composition_lanes> fraction{};
};

struct CompositionCorrectionRecord {
  std::uint64_t gas_cell_id = 0U;
  int source_rank = 0;
  int owner_rank = 0;
  std::uint32_t lane_count = 0U;
  std::uint32_t reserved = 0U;
  std::array<double, k_max_composition_lanes> lane_delta_code{};
};

// Byte-level primitive: record i (record_bytes wide) goes to destination[i].
[[nodiscard]] std::vector<std::byte> executeBlockingDestinationRecordExchange(
    const MpiContext& mpi_context,
    std::span<const std::byte> records,
    std::size_t record_bytes,
    std::span<const int> destination_rank,
    std::uint64_t exchange_sequence);

[[nodiscard]] std::vector<CompositionDonorRequestRecord> exchangeCompositionDonorRequests(
    const MpiContext& mpi_context,
    std::span<const CompositionDonorRequestRecord> requests,
    std::uint64_t exchange_sequence);

[[nodiscard]] std::vector<CompositionDonorFractionRecord> exchangeCompositionDonorFractions(
    const MpiContext& mpi_context,
    std::span<const CompositionDonorFractionRecord> responses,
    std::uint64_t exchange_sequence);

[[nodiscard]] std::vector<CompositionCorrectionRecord> exchangeCompositionCorrections(
    const MpiContext& mpi_context,
    std::span<const CompositionCorrectionRecord> corrections,
    std::uint64_t exchange_sequence);

}  // namespace cosmosim::parallel
