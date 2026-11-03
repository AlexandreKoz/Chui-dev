#pragma once

#include <cstdint>
#include <string_view>

#include "cosmosim/core/stage_transaction.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

// Communicator-wide StageTransaction vote: one MPI_Allreduce(MIN) of a 0/1
// flag. Serial and single-rank executions return the local vote without any
// MPI call, so non-MPI builds never touch MPI. Must be invoked by every rank
// at the same transaction boundary, outside OpenMP parallel regions.
class MpiStageCollectiveDecision final : public core::StageCollectiveDecision {
 public:
  explicit MpiStageCollectiveDecision(const parallel::MpiContext& mpi_context) noexcept
      : m_mpi_context(mpi_context) {}

  [[nodiscard]] bool allAccept(bool local_accept, std::string_view phase) const override {
    (void)phase;
    if (!m_mpi_context.isEnabled() || m_mpi_context.worldSize() <= 1) {
      return local_accept;
    }
    return m_mpi_context.allreduceMinUint64(local_accept ? 1U : 0U) == 1U;
  }

 private:
  const parallel::MpiContext& m_mpi_context;
};

}  // namespace cosmosim::workflows::internal
