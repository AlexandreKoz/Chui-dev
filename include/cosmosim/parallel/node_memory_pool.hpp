#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "cosmosim/core/node_memory_admission.hpp"

namespace cosmosim::parallel {

class MpiContext;

// Node-shared accounting counter for parallel.node_memory_admission =
// shared_pool.
//
// Mechanism. MPI builds with more than one world rank: the ranks of each
// shared-memory node (MPI_COMM_TYPE_SHARED) collectively allocate one 8-byte
// MPI-3 shared window segment on node-local rank 0 and update that counter
// with lock-free std::atomic_ref compare-and-swap (no MPI call per
// reservation, MPI_THREAD_FUNNELED is respected, nobody waits, no deadlock).
//
// The C++ standard does not by itself guarantee that atomics on memory mapped
// into several processes are atomic across those processes; it only
// recommends that lock-free atomics be address-free. The guarantee is
// therefore CLAIMED ONLY when every condition below holds on every rank:
//   1. std::atomic_ref<uint64_t> is always lock-free (compile time) and the
//      toolchain/ISA pair is one whose lock-free atomics are address-free
//      hardware instructions on cache-coherent memory: GCC/Clang on x86-64,
//      AArch64 or POWER64 (compile-time allowlist);
//   2. the window reports MPI_WIN_MODEL == MPI_WIN_UNIFIED (public and
//      private copies are the same memory, so direct load/store is coherent);
//   3. the counter address returned by MPI_Win_shared_query is 8-byte aligned;
//   4. all of the above are agreed by an MPI_Allreduce(MIN) over the world at
//      construction, so either every rank uses the pool or none does.
// If any condition fails the window is released and established() is false:
// the caller must fall back to the deterministic shared_equal ceilings (the
// reference workflow does, and reports the reason). Serial / non-MPI or
// single-rank runs use a process-local counter, which is trivially correct.
//
// Lifetime. Construction and normal destruction are COLLECTIVE over
// MPI_COMM_WORLD and must run on the main thread at the same workflow
// boundary on all ranks. Governors must be detached (or destroyed) before the
// pool. During exception unwinding the destructor does not call the collective
// MPI_Win_free (a failing rank must not hang waiting for peers); the window is
// left to MPI_Finalize/MPI_Abort.
class NodeSharedMemoryPool final {
 public:
  NodeSharedMemoryPool(const MpiContext& mpi_context, std::uint64_t node_capacity_bytes);
  ~NodeSharedMemoryPool();

  NodeSharedMemoryPool(const NodeSharedMemoryPool&) = delete;
  NodeSharedMemoryPool& operator=(const NodeSharedMemoryPool&) = delete;
  NodeSharedMemoryPool(NodeSharedMemoryPool&&) = delete;
  NodeSharedMemoryPool& operator=(NodeSharedMemoryPool&&) = delete;

  // True when the node-wide sum guarantee holds (see conditions above).
  [[nodiscard]] bool established() const noexcept;
  // Reason when not established (empty otherwise).
  [[nodiscard]] const std::string& fallbackReason() const noexcept;
  // Requires established().
  [[nodiscard]] core::NodeMemoryPool& pool();
  [[nodiscard]] bool usesMpiSharedWindow() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

// Compile-time part of condition 1 (exposed for tests and diagnostics).
[[nodiscard]] bool nodeSharedAtomicsPlatformSupported() noexcept;

}  // namespace cosmosim::parallel
