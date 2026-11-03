#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace cosmosim::core {

// Node-aware CPU memory admission.
//
// Problem. MemoryGovernor enforces a RANK-local hard limit. Co-located MPI
// ranks share one node's physical memory (and usually one job cgroup), so a
// per-rank limit derived from "node RAM" without dividing by the number of
// co-located ranks over-commits the node by that factor.
//
// Policy (provably conservative, no runtime communication):
//   node_limit   = min(explicit node budget, cgroup limit, scheduler job
//                      limit, physical RAM) over the values that are known
//   node_usable  = node_limit - node_reserve            (subtracted ONCE)
//   rank_share   = floor(node_usable / ranks_on_node)
//   rank_limit   = process_budget == 0 ? rank_share : min(process_budget, rank_share)
// Because sum_r rank_limit_r <= ranks_on_node * rank_share <= node_usable,
// the governed demand of all co-located ranks can never exceed the node
// limit minus the reserve, regardless of how unevenly the ranks allocate.
// Rank-local external-library reserves (MPI/FFTW/HDF5/allocator) remain in
// MemoryGovernorPolicy::external_runtime_reserve_bytes and are charged INSIDE
// each rank's share, so nothing is subtracted twice.
//
// Assumptions / limits (documented, not hidden):
//  * A discovered cgroup or scheduler limit is assumed to be shared by all
//    co-located ranks of this job (true for job-level cgroups; conservative
//    when each task has its own cgroup).
//  * Physical RAM alone is optimistic on nodes shared with other jobs; set an
//    explicit node budget on such systems.
//  * Governed accounting is an estimate of CHUI-owned allocations, not OS
//    resident memory; ProcessMemoryObservation remains diagnostic evidence.
//  * If no node limit can be determined and no explicit budget is given, the
//    admission refuses (std::runtime_error) instead of guessing.
enum class NodeMemoryLimitSource : std::uint8_t {
  kNone = 0,
  kExplicitNodeBudget = 1,
  kCgroupV2 = 2,
  kCgroupV1 = 3,
  kSchedulerJob = 4,
  kPhysicalRam = 5,
};

[[nodiscard]] std::string_view nodeMemoryLimitSourceLabel(NodeMemoryLimitSource source) noexcept;

struct NodeMemoryDiscovery {
  std::optional<std::uint64_t> physical_ram_bytes;
  std::optional<std::uint64_t> cgroup_limit_bytes;
  NodeMemoryLimitSource cgroup_source = NodeMemoryLimitSource::kNone;
  // SLURM_MEM_PER_NODE (MiB) converted to bytes when present.
  std::optional<std::uint64_t> scheduler_job_limit_bytes;
};

// Pure parsers (unit-testable without touching the host).
[[nodiscard]] std::optional<std::uint64_t> parseProcMeminfoTotalBytes(std::string_view meminfo_text) noexcept;
// cgroup v2 memory.max / v1 memory.limit_in_bytes. "max" and v1 "unlimited"
// sentinels (>= 2^62) return std::nullopt.
[[nodiscard]] std::optional<std::uint64_t> parseCgroupMemoryLimitBytes(std::string_view text) noexcept;
[[nodiscard]] std::optional<std::uint64_t> parseSchedulerMemPerNodeMiB(std::string_view text) noexcept;

// Host discovery (Linux /proc, /sys/fs/cgroup, sysconf, SLURM environment).
// Never throws; unavailable values stay std::nullopt.
[[nodiscard]] NodeMemoryDiscovery discoverNodeMemory() noexcept;

struct NodeMemoryAdmissionInput {
  std::uint64_t explicit_node_budget_bytes = 0U;  // 0 = not set
  std::uint64_t node_reserve_bytes = 0U;
  std::uint32_t ranks_on_node = 1U;
  std::uint64_t process_budget_bytes = 0U;  // 0 = unlimited (historical semantics)
  NodeMemoryDiscovery discovery{};
};

struct NodeMemoryAdmissionResult {
  NodeMemoryLimitSource node_limit_source = NodeMemoryLimitSource::kNone;
  std::uint64_t node_limit_bytes = 0U;
  std::uint64_t node_usable_bytes = 0U;
  std::uint64_t rank_share_bytes = 0U;
  // Value for MemoryGovernorPolicy::hard_limit_bytes (never 0 here: a node
  // limited admission always yields a finite rank limit).
  std::uint64_t rank_hard_limit_bytes = 0U;
  std::uint32_t ranks_on_node = 1U;
  bool process_budget_was_binding = false;
  std::string explanation;
};

// Throws std::invalid_argument for ranks_on_node == 0 or a reserve that
// consumes the whole node limit, std::runtime_error when no node limit is
// known, and std::overflow_error never (all arithmetic is checked/divisive).
[[nodiscard]] NodeMemoryAdmissionResult computeNodeMemoryAdmission(const NodeMemoryAdmissionInput& input);

}  // namespace cosmosim::core
