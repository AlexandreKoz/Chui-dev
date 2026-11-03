#include "cosmosim/parallel/node_memory_pool.hpp"

#include <atomic>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {

bool nodeSharedAtomicsPlatformSupported() noexcept {
#if (defined(__GNUC__) || defined(__clang__)) && \
    (defined(__x86_64__) || defined(__aarch64__) || defined(__powerpc64__))
  return std::atomic_ref<std::uint64_t>::is_always_lock_free;
#else
  return false;
#endif
}

struct NodeSharedMemoryPool::Impl {
  std::unique_ptr<core::ProcessLocalNodeMemoryPool> local_pool;
  std::unique_ptr<core::AtomicCounterNodeMemoryPool> shared_pool;
  std::string fallback_reason;
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
    // Every step below is executed by every rank; the local verdict is only
    // acted upon after the world-wide agreement, so ranks never diverge.
    std::string local_reason;
    if (!nodeSharedAtomicsPlatformSupported()) {
      local_reason = "toolchain/ISA not on the address-free lock-free 64-bit atomic allowlist";
    }
    if (MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &m_impl->node_comm) !=
            MPI_SUCCESS ||
        m_impl->node_comm == MPI_COMM_NULL) {
      throw std::runtime_error("node-shared memory pool: MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) failed");
    }
    int node_rank = 0;
    MPI_Comm_rank(m_impl->node_comm, &node_rank);
    void* base = nullptr;
    const MPI_Aint local_bytes = node_rank == 0 ? static_cast<MPI_Aint>(sizeof(std::uint64_t)) : 0;
    const int allocate_status = MPI_Win_allocate_shared(
        local_bytes, static_cast<int>(sizeof(std::uint64_t)), MPI_INFO_NULL, m_impl->node_comm, &base,
        &m_impl->window);
    void* counter_address = nullptr;
    if (allocate_status != MPI_SUCCESS || m_impl->window == MPI_WIN_NULL) {
      local_reason = "MPI_Win_allocate_shared failed";
    } else {
      int* memory_model = nullptr;
      int has_model = 0;
      MPI_Win_get_attr(m_impl->window, MPI_WIN_MODEL, &memory_model, &has_model);
      if (has_model == 0 || memory_model == nullptr || *memory_model != MPI_WIN_UNIFIED) {
        local_reason = "MPI window memory model is not MPI_WIN_UNIFIED";
      }
      MPI_Aint segment_bytes = 0;
      int displacement_unit = 0;
      MPI_Win_shared_query(m_impl->window, 0, &segment_bytes, &displacement_unit, &counter_address);
      if (counter_address == nullptr || segment_bytes < static_cast<MPI_Aint>(sizeof(std::uint64_t))) {
        local_reason = "shared counter segment unavailable";
      } else if (reinterpret_cast<std::uintptr_t>(counter_address) % alignof(std::uint64_t) != 0U) {
        local_reason = "shared counter is not 8-byte aligned";
      }
    }
    int local_ok = local_reason.empty() ? 1 : 0;
    int all_ok = 0;
    MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (all_ok == 1) {
      auto* counter = static_cast<std::uint64_t*>(counter_address);
      if (node_rank == 0) {
        std::atomic_ref<std::uint64_t>(*counter).store(0U, std::memory_order_release);
      }
      // Publish the zeroed counter before any rank reserves.
      MPI_Win_sync(m_impl->window);
      MPI_Barrier(m_impl->node_comm);
      m_impl->shared_pool = std::make_unique<core::AtomicCounterNodeMemoryPool>(counter, node_capacity_bytes);
      return;
    }
    m_impl->fallback_reason = local_reason.empty()
        ? std::string("another rank could not establish the node-shared counter")
        : local_reason;
    if (m_impl->window != MPI_WIN_NULL) {
      MPI_Win_free(&m_impl->window);
    }
    MPI_Comm_free(&m_impl->node_comm);
    return;
  }
#else
  (void)mpi_context;
#endif
  // Single process: a process-local counter is exact.
  m_impl->local_pool = std::make_unique<core::ProcessLocalNodeMemoryPool>(node_capacity_bytes);
}

NodeSharedMemoryPool::~NodeSharedMemoryPool() {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_impl != nullptr && m_impl->window != MPI_WIN_NULL) {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized == 0 && std::uncaught_exceptions() == 0) {
      MPI_Win_free(&m_impl->window);
      MPI_Comm_free(&m_impl->node_comm);
    }
  }
#endif
}

bool NodeSharedMemoryPool::established() const noexcept {
  return m_impl->shared_pool != nullptr || m_impl->local_pool != nullptr;
}

const std::string& NodeSharedMemoryPool::fallbackReason() const noexcept {
  return m_impl->fallback_reason;
}

core::NodeMemoryPool& NodeSharedMemoryPool::pool() {
  if (m_impl->shared_pool != nullptr) {
    return *m_impl->shared_pool;
  }
  if (m_impl->local_pool != nullptr) {
    return *m_impl->local_pool;
  }
  throw std::logic_error("node-shared memory pool was not established: " + m_impl->fallback_reason);
}

bool NodeSharedMemoryPool::usesMpiSharedWindow() const noexcept {
  return m_impl->shared_pool != nullptr;
}

}  // namespace cosmosim::parallel
