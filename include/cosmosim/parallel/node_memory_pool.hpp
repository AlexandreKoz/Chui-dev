#pragma once

#include <cstdint>
#include <memory>

#include "cosmosim/core/node_memory_admission.hpp"

namespace cosmosim::parallel {

class MpiContext;

// Node-shared accounting counter for parallel.node_memory_admission =
// shared_pool.
//
// MPI builds with an active world: the ranks of each shared-memory node
// (MPI_COMM_TYPE_SHARED) collectively allocate one 8-byte MPI-3 shared window
// segment on node-local rank 0; every rank then updates that counter with
// lock-free std::atomic_ref operations (no MPI call per reservation, so the
// MPI_THREAD_FUNNELED model is respected and no rank can block another).
// Construction and destruction are COLLECTIVE over MPI_COMM_WORLD (node
// split + window allocation/free) and must be executed by every rank at the
// same workflow boundary, on the main thread.
//
// Serial / non-MPI builds or an inactive world: a process-local pool.
class NodeSharedMemoryPool final {
 public:
  NodeSharedMemoryPool(const MpiContext& mpi_context, std::uint64_t node_capacity_bytes);
  ~NodeSharedMemoryPool();

  NodeSharedMemoryPool(const NodeSharedMemoryPool&) = delete;
  NodeSharedMemoryPool& operator=(const NodeSharedMemoryPool&) = delete;
  NodeSharedMemoryPool(NodeSharedMemoryPool&&) = delete;
  NodeSharedMemoryPool& operator=(NodeSharedMemoryPool&&) = delete;

  [[nodiscard]] core::NodeMemoryPool& pool() noexcept;
  [[nodiscard]] bool usesMpiSharedWindow() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace cosmosim::parallel
