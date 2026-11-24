#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"

namespace cosmosim::mhd {

// Exact retained-capacity model of the MHD allocations made by
// MhdLevelStepper (CHUI-MHD-PERF v1.0 sections 4-6, 9). Every term mirrors one
// MhdDoubleArray lane allocated by the stepper/workspaces, rounded to the
// 64-byte allocator request (mhdAlignedLaneBytes), so MP01 can compare it
// against ownedCapacityBytes() after allocation. Formulas, per patch of
// (nx, ny, nz) cells, ghost g (2 PLM / 3 PPM), staggered tangential ghost 2:
//   primitives      9 x (nx+2g)(ny+2g)(nz+2g)
//   staggered       sum_a (n_a+1)(n_t1+4)(n_t2+4)
//   face signals    2 x sum_a (n_a+1)(n_t1+2)(n_t2+2)
//   stage EMF       N_e = nx(ny+1)(nz+1) + (nx+1)ny(nz+1) + (nx+1)(ny+1)nz
//   boundary flux   6 x 2 x sum_a n_t1 n_t2
//   rollback image  6 x N_c (gas) + N_f (faces)          [kPhaseResident]
//   stage-0 image   6 x N_c (gas) + N_f (faces)          [kPhaseResident]
//   non-ideal       3 x (nx+2)(ny+2)(nz+2) coefficients + N_c stage rate
//                   + 2 x sum_a n_t1 n_t2 boundary energy
//     super-step    N_f (Y_{j-2} faces) + 2 N_c + N_e    (RKL2/explicit lanes)
//   edge integral   3 x N_e (when produce_edge_integrals)
// Canonical face flux 8 N_f is charged separately as kCanonicalPersistent.
struct MhdPatchMemoryEstimate {
  std::uint64_t canonical_face_flux_bytes = 0U;
  std::uint64_t primitive_cache_bytes = 0U;
  std::uint64_t staggered_cache_bytes = 0U;
  std::uint64_t face_signal_bytes = 0U;
  std::uint64_t stage_emf_bytes = 0U;
  std::uint64_t boundary_flux_bytes = 0U;
  std::uint64_t rollback_image_bytes = 0U;
  std::uint64_t stage0_image_bytes = 0U;
  std::uint64_t nonideal_bytes = 0U;
  std::uint64_t edge_integral_bytes = 0U;

  [[nodiscard]] std::uint64_t phaseResidentBytes() const noexcept {
    return rollback_image_bytes + stage0_image_bytes;
  }
  [[nodiscard]] std::uint64_t scratchArenaBytes() const noexcept {
    return primitive_cache_bytes + staggered_cache_bytes + face_signal_bytes + stage_emf_bytes +
        boundary_flux_bytes + nonideal_bytes + edge_integral_bytes;
  }
  [[nodiscard]] std::uint64_t phasePlusScratchBytes() const noexcept {
    return phaseResidentBytes() + scratchArenaBytes();
  }
};

struct MhdMemoryRequest {
  MhdPatchGeometry geometry{};
  MhdReconstructionMode reconstruction = MhdReconstructionMode::kPpm;
  bool nonideal = false;
  bool super_stepping = false;
  bool edge_integrals = false;
};

[[nodiscard]] MhdPatchMemoryEstimate estimatePatchMemory(const MhdMemoryRequest& request);

// Exact geometric counts (PERF section 4).
[[nodiscard]] std::uint64_t mhdFaceCount(const MhdPatchGeometry& geometry);
[[nodiscard]] std::uint64_t mhdEdgeCount(const MhdPatchGeometry& geometry);

// Resource limits (PERF sections 3, 6, 8, 10, 12, 17). Defaults are the
// workstation_48gib_v01 profile.
struct MhdResourceLimits {
  std::uint64_t scratch_bytes_per_active_cell = k_mhd_max_scratch_bytes_per_active_cell;
  std::uint64_t comm_budget_bytes_per_rank = k_mhd_default_comm_budget_bytes_per_rank;
  std::uint64_t comm_max_bytes_per_rank = k_mhd_max_comm_budget_bytes_per_rank;
  std::uint64_t io_chunk_bytes_per_rank = k_mhd_io_chunk_target_bytes_per_rank;
  std::uint64_t host_safe_bytes = 42ULL * 1024ULL * 1024ULL * 1024ULL;
  std::uint64_t host_nominal_bytes = 36ULL * 1024ULL * 1024ULL * 1024ULL;
  std::uint64_t gpu_device_target_bytes = 6ULL * 1024ULL * 1024ULL * 1024ULL;
  std::uint32_t active_batch_target_patches = 8U;
};

// MP02 check: MHD-owned phase-resident + scratch bytes per active real cell
// for the co-resident set of `patch_count` patches of the request shape.
[[nodiscard]] double scratchBytesPerActiveCell(const MhdMemoryRequest& request,
                                               std::size_t patch_count);

// Deterministic batch size: the largest number of patches (<= target, >= 1)
// whose phase+scratch fits both the per-cell cap and `available_bytes`.
// Returns 0 when even one patch does not fit (caller must fail preflight).
[[nodiscard]] std::size_t selectMhdBatchPatches(const MhdMemoryRequest& request,
                                                std::uint64_t available_bytes,
                                                const MhdResourceLimits& limits);

// Whole-host aggregate preflight (PERF section 3.1): the sum over ranks on one
// host of their governor hard limits must not exceed host_safe_bytes; a rank
// hard limit of 0 (unlimited) is not admissible for a performance-qualified
// MHD run (MP11). Throws MHD_MEMORY_BUDGET_EXCEEDED.
void requireHostEnvelope(std::span<const std::uint64_t> rank_hard_limits_on_host,
                         const MhdResourceLimits& limits,
                         bool performance_qualified);

// Governed admission of an MHD workspace set BEFORE any of it is allocated.
// Holds one kPhaseResident and one kScratchArena reservation for its lifetime
// (retained capacity), plus the canonical face flux charge when requested.
// A refusal is mapped to MHD_MEMORY_BUDGET_EXCEEDED with requested/safe bytes
// and never mutates state.
class MhdMemoryAdmission {
 public:
  MhdMemoryAdmission() = default;
  MhdMemoryAdmission(core::MemoryGovernor& governor,
                     const MhdPatchMemoryEstimate& per_patch,
                     std::size_t patch_count,
                     bool charge_canonical_faces,
                     std::string owner);

  [[nodiscard]] std::uint64_t phaseResidentBytes() const noexcept { return m_phase_bytes; }
  [[nodiscard]] std::uint64_t scratchBytes() const noexcept { return m_scratch_bytes; }
  [[nodiscard]] std::uint64_t canonicalBytes() const noexcept { return m_canonical_bytes; }
  [[nodiscard]] bool admitted() const noexcept { return m_admitted; }
  // Commit after the allocations exist; verifies the actual retained capacity
  // does not exceed the admitted bytes (unexpected growth is a failure).
  void commit(std::uint64_t actual_phase_plus_scratch_bytes);
  void release() noexcept;

 private:
  core::MemoryReservation m_phase{};
  core::MemoryReservation m_scratch{};
  core::MemoryReservation m_canonical{};
  std::uint64_t m_phase_bytes = 0U;
  std::uint64_t m_scratch_bytes = 0U;
  std::uint64_t m_canonical_bytes = 0U;
  bool m_admitted = false;
  std::string m_owner;
};

}  // namespace cosmosim::mhd
