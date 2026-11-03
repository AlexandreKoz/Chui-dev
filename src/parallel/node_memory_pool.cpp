#include "cosmosim/parallel/node_memory_pool.hpp"

#include <stdexcept>
#include <string>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {

struct NodeSharedMemoryPool::Impl {
  std::unique_ptr<core::ProcessLocalNodeMemoryPool> local_pool;
  std::unique_ptr<core::AtomicCounterNodeMemoryPool> shared_pool;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  MPI_Comm node_comm = MPI_COMM_NULL;
  MPI_Win window = MPI_WIN_NULL;
#endif
};

NodeSharedMemoryPool::NodeSharedMemoryPool(const MpiContext& mpi_context, std::uint64_t node_capacity_bytes)
    : m_impl(std::make_unique<Impl>()) {
  if (node_capacity_bytes == 0U) {
    throw std::invalid_argument("node-shared memory pool requires a positive node capacity");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (mpi_context.isEnabled() && mpi_context.worldSize() > 1) {
    if (MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &m_impl->node_comm) !=
            MPI_SUCCESS ||
        m_impl->node_comm == MPI_COMM_NULL) {
      throw std::runtime_error("node-shared memory pool: MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) failed");
    }
    int node_rank = 0;
    MPI_Comm_rank(m_impl->node_comm, &node_rank);
    void* base = nullptr;
    const MPI_Aint local_bytes = node_rank == 0 ? static_cast<MPI_Aint>(sizeof(std::uint64_t)) : 0;
    if (MPI_Win_allocate_shared(local_bytes, static_cast<int>(sizeof(std::uint64_t)), MPI_INFO_NULL,
                                m_impl->node_comm, &base, &m_impl->window) != MPI_SUCCESS) {
      MPI_Comm_free(&m_impl->node_comm);
      throw std::runtime_error("node-shared memory pool: MPI_Win_allocate_shared failed");
    }
    MPI_Aint segment_bytes = 0;
    int displacement_unit = 0;
    void* counter_address = nullptr;
    MPI_Win_shared_query(m_impl->window, 0, &segment_bytes, &displacement_unit, &counter_address);
    if (counter_address == nullptr || segment_bytes < static_cast<MPI_Aint>(sizeof(std::uint64_t))) {
      MPI_Win_free(&m_impl->window);
      MPI_Comm_free(&m_impl->node_comm);
      throw std::runtime_error("node-shared memory pool: shared counter segment unavailable");
    }
    auto* counter = static_cast<std::uint64_t*>(counter_address);
    if (node_rank == 0) {
      *counter = 0U;
    }
    // Publish the zeroed counter before any rank reserves.
    MPI_Win_sync(m_impl->window);
    MPI_Barrier(m_impl->node_comm);
    m_impl->shared_pool = std::make_unique<core::AtomicCounterNodeMemoryPool>(counter, node_capacity_bytes);
    return;
  }
#else
  (void)mpi_context;
#endif
  m_impl->local_pool = std::make_unique<core::ProcessLocalNodeMemoryPool>(node_capacity_bytes);
}

NodeSharedMemoryPool::~NodeSharedMemoryPool() {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_impl != nullptr && m_impl->window != MPI_WIN_NULL) {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized == 0) {
      MPI_Win_free(&m_impl->window);
      MPI_Comm_free(&m_impl->node_comm);
    }
  }
#endif
}

core::NodeMemoryPool& NodeSharedMemoryPool::pool() noexcept {
  if (m_impl->shared_pool != nullptr) {
    return *m_impl->shared_pool;
  }
  return *m_impl->local_pool;
}

bool NodeSharedMemoryPool::usesMpiSharedWindow() const noexcept {
  return m_impl->shared_pool != nullptr;
}

}  // namespace cosmosim::parallel
