#include "cosmosim/core/node_memory_admission.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace cosmosim::core {
namespace {

constexpr std::uint64_t k_unlimited_sentinel_floor = 1ULL << 62U;

[[nodiscard]] std::string_view trimView(std::string_view text) noexcept {
  while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front())) != 0) {
    text.remove_prefix(1U);
  }
  while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back())) != 0) {
    text.remove_suffix(1U);
  }
  return text;
}

[[nodiscard]] std::optional<std::uint64_t> parseUnsigned(std::string_view text) noexcept {
  text = trimView(text);
  if (text.empty()) {
    return std::nullopt;
  }
  std::uint64_t value = 0U;
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
  if (error != std::errc{} || end != text.data() + text.size()) {
    return std::nullopt;
  }
  return value;
}

[[nodiscard]] std::optional<std::string> readSmallFile(const std::string& path) noexcept {
  try {
    std::ifstream stream(path);
    if (!stream) {
      return std::nullopt;
    }
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
  } catch (...) {
    return std::nullopt;
  }
}

[[nodiscard]] std::optional<std::uint64_t> checkedMiBToBytes(std::uint64_t mib) noexcept {
  constexpr std::uint64_t k_mib = 1024ULL * 1024ULL;
  if (mib > std::numeric_limits<std::uint64_t>::max() / k_mib) {
    return std::nullopt;
  }
  return mib * k_mib;
}

}  // namespace

std::string_view nodeMemoryLimitSourceLabel(NodeMemoryLimitSource source) noexcept {
  switch (source) {
    case NodeMemoryLimitSource::kNone:
      return "none";
    case NodeMemoryLimitSource::kExplicitNodeBudget:
      return "explicit_node_budget";
    case NodeMemoryLimitSource::kCgroupV2:
      return "cgroup_v2";
    case NodeMemoryLimitSource::kCgroupV1:
      return "cgroup_v1";
    case NodeMemoryLimitSource::kSchedulerJob:
      return "scheduler_job";
    case NodeMemoryLimitSource::kPhysicalRam:
      return "physical_ram";
  }
  return "unknown";
}

std::optional<std::uint64_t> parseProcMeminfoTotalBytes(std::string_view meminfo_text) noexcept {
  constexpr std::string_view k_key = "MemTotal:";
  const std::size_t position = meminfo_text.find(k_key);
  if (position == std::string_view::npos) {
    return std::nullopt;
  }
  std::string_view rest = meminfo_text.substr(position + k_key.size());
  const std::size_t line_end = rest.find('\n');
  if (line_end != std::string_view::npos) {
    rest = rest.substr(0U, line_end);
  }
  rest = trimView(rest);
  const std::size_t unit_position = rest.find(' ');
  std::string_view number = rest.substr(0U, unit_position);
  std::string_view unit = unit_position == std::string_view::npos
      ? std::string_view{} : trimView(rest.substr(unit_position));
  const std::optional<std::uint64_t> value = parseUnsigned(number);
  if (!value.has_value()) {
    return std::nullopt;
  }
  if (unit.empty()) {
    return value;
  }
  if (unit == "kB" || unit == "KB" || unit == "kb") {
    if (*value > std::numeric_limits<std::uint64_t>::max() / 1024ULL) {
      return std::nullopt;
    }
    return *value * 1024ULL;
  }
  return std::nullopt;
}

std::optional<std::uint64_t> parseCgroupMemoryLimitBytes(std::string_view text) noexcept {
  text = trimView(text);
  if (text.empty() || text == "max") {
    return std::nullopt;
  }
  const std::optional<std::uint64_t> value = parseUnsigned(text);
  if (!value.has_value() || *value == 0U || *value >= k_unlimited_sentinel_floor) {
    return std::nullopt;
  }
  return value;
}

std::optional<std::uint64_t> parseSchedulerMemPerNodeMiB(std::string_view text) noexcept {
  text = trimView(text);
  // SLURM may append a unit suffix (M/G); accept plain MiB and M/G suffixes.
  std::uint64_t multiplier_mib = 1U;
  if (!text.empty() && (text.back() == 'M' || text.back() == 'm')) {
    text.remove_suffix(1U);
  } else if (!text.empty() && (text.back() == 'G' || text.back() == 'g')) {
    text.remove_suffix(1U);
    multiplier_mib = 1024U;
  }
  const std::optional<std::uint64_t> value = parseUnsigned(text);
  if (!value.has_value() || *value == 0U) {
    return std::nullopt;
  }
  if (*value > std::numeric_limits<std::uint64_t>::max() / multiplier_mib) {
    return std::nullopt;
  }
  return checkedMiBToBytes(*value * multiplier_mib);
}

NodeMemoryDiscovery discoverNodeMemory() noexcept {
  NodeMemoryDiscovery discovery;
  try {
    if (const auto meminfo = readSmallFile("/proc/meminfo"); meminfo.has_value()) {
      discovery.physical_ram_bytes = parseProcMeminfoTotalBytes(*meminfo);
    }
#if defined(__unix__) || defined(__APPLE__)
    if (!discovery.physical_ram_bytes.has_value()) {
      const long pages = ::sysconf(_SC_PHYS_PAGES);
      const long page_size = ::sysconf(_SC_PAGESIZE);
      if (pages > 0 && page_size > 0 &&
          static_cast<std::uint64_t>(pages) <=
              std::numeric_limits<std::uint64_t>::max() / static_cast<std::uint64_t>(page_size)) {
        discovery.physical_ram_bytes =
            static_cast<std::uint64_t>(pages) * static_cast<std::uint64_t>(page_size);
      }
    }
#endif
    // cgroup membership: "0::/path" (v2) or "N:memory:/path" (v1).
    if (const auto cgroup = readSmallFile("/proc/self/cgroup"); cgroup.has_value()) {
      std::istringstream lines(*cgroup);
      std::string line;
      while (std::getline(lines, line)) {
        if (line.rfind("0::", 0U) == 0U) {
          const std::string relative = line.substr(3U);
          for (const std::string& path : {"/sys/fs/cgroup" + relative + "/memory.max",
                                          std::string("/sys/fs/cgroup/memory.max")}) {
            if (const auto text = readSmallFile(path); text.has_value()) {
              if (const auto limit = parseCgroupMemoryLimitBytes(*text); limit.has_value()) {
                discovery.cgroup_limit_bytes = limit;
                discovery.cgroup_source = NodeMemoryLimitSource::kCgroupV2;
                break;
              }
            }
          }
        } else if (const std::size_t marker = line.find(":memory:"); marker != std::string::npos) {
          const std::string relative = line.substr(marker + 8U);
          for (const std::string& path :
               {"/sys/fs/cgroup/memory" + relative + "/memory.limit_in_bytes",
                std::string("/sys/fs/cgroup/memory/memory.limit_in_bytes")}) {
            if (const auto text = readSmallFile(path); text.has_value()) {
              if (const auto limit = parseCgroupMemoryLimitBytes(*text); limit.has_value()) {
                if (!discovery.cgroup_limit_bytes.has_value() || *limit < *discovery.cgroup_limit_bytes) {
                  discovery.cgroup_limit_bytes = limit;
                  discovery.cgroup_source = NodeMemoryLimitSource::kCgroupV1;
                }
                break;
              }
            }
          }
        }
      }
    }
    if (const char* slurm = std::getenv("SLURM_MEM_PER_NODE"); slurm != nullptr) {
      discovery.scheduler_job_limit_bytes = parseSchedulerMemPerNodeMiB(slurm);
    }
  } catch (...) {
    // Discovery is best effort by contract; missing values stay nullopt.
  }
  return discovery;
}

NodeMemoryAdmissionResult computeNodeMemoryAdmission(const NodeMemoryAdmissionInput& input) {
  if (input.ranks_on_node == 0U) {
    throw std::invalid_argument("node memory admission requires ranks_on_node > 0");
  }
  NodeMemoryAdmissionResult result;
  result.ranks_on_node = input.ranks_on_node;
  const auto consider = [&result](std::optional<std::uint64_t> value, NodeMemoryLimitSource source) {
    if (!value.has_value() || *value == 0U) {
      return;
    }
    if (result.node_limit_source == NodeMemoryLimitSource::kNone || *value < result.node_limit_bytes) {
      result.node_limit_bytes = *value;
      result.node_limit_source = source;
    }
  };
  consider(input.explicit_node_budget_bytes != 0U
               ? std::optional<std::uint64_t>(input.explicit_node_budget_bytes)
               : std::nullopt,
           NodeMemoryLimitSource::kExplicitNodeBudget);
  consider(input.discovery.cgroup_limit_bytes, input.discovery.cgroup_source == NodeMemoryLimitSource::kNone
                                                   ? NodeMemoryLimitSource::kCgroupV2
                                                   : input.discovery.cgroup_source);
  consider(input.discovery.scheduler_job_limit_bytes, NodeMemoryLimitSource::kSchedulerJob);
  consider(input.discovery.physical_ram_bytes, NodeMemoryLimitSource::kPhysicalRam);
  if (result.node_limit_source == NodeMemoryLimitSource::kNone) {
    throw std::runtime_error(
        "node memory admission could not determine any node memory limit (no explicit "
        "parallel.node_memory_budget_bytes, cgroup, scheduler, or physical-RAM value); "
        "set parallel.node_memory_budget_bytes explicitly");
  }
  if (input.node_reserve_bytes >= result.node_limit_bytes) {
    throw std::invalid_argument(
        "parallel.node_memory_reserve_bytes consumes the entire node memory limit (" +
        std::to_string(result.node_limit_bytes) + " bytes from " +
        std::string(nodeMemoryLimitSourceLabel(result.node_limit_source)) + ")");
  }
  result.node_usable_bytes = result.node_limit_bytes - input.node_reserve_bytes;
  result.rank_share_bytes = result.node_usable_bytes / input.ranks_on_node;
  if (result.rank_share_bytes == 0U) {
    throw std::invalid_argument("node memory admission produced a zero-byte rank share");
  }
  result.rank_hard_limit_bytes = result.rank_share_bytes;
  if (input.process_budget_bytes != 0U && input.process_budget_bytes < result.rank_share_bytes) {
    result.rank_hard_limit_bytes = input.process_budget_bytes;
    result.process_budget_was_binding = true;
  }
  std::ostringstream explanation;
  explanation << "node_limit=" << result.node_limit_bytes
              << " source=" << nodeMemoryLimitSourceLabel(result.node_limit_source)
              << " node_reserve=" << input.node_reserve_bytes
              << " ranks_on_node=" << input.ranks_on_node
              << " rank_share=" << result.rank_share_bytes
              << " rank_hard_limit=" << result.rank_hard_limit_bytes
              << (result.process_budget_was_binding ? " (process budget binding)" : "");
  result.explanation = explanation.str();
  return result;
}

static_assert(std::atomic_ref<std::uint64_t>::is_always_lock_free,
              "node-shared memory accounting requires lock-free 64-bit atomics");

AtomicCounterNodeMemoryPool::AtomicCounterNodeMemoryPool(std::uint64_t* counter, std::uint64_t capacity_bytes)
    : m_counter(counter), m_capacity_bytes(capacity_bytes) {
  if (m_counter == nullptr) {
    throw std::invalid_argument("node memory pool requires a counter");
  }
  if (reinterpret_cast<std::uintptr_t>(m_counter) % alignof(std::uint64_t) != 0U) {
    throw std::invalid_argument("node memory pool counter must be 8-byte aligned");
  }
  if (m_capacity_bytes == 0U) {
    throw std::invalid_argument("node memory pool requires a positive capacity");
  }
}

bool AtomicCounterNodeMemoryPool::tryAcquire(std::uint64_t bytes) noexcept {
  if (bytes == 0U) {
    return true;
  }
  std::atomic_ref<std::uint64_t> counter(*m_counter);
  std::uint64_t used = counter.load(std::memory_order_relaxed);
  while (true) {
    if (bytes > m_capacity_bytes || used > m_capacity_bytes - bytes) {
      return false;
    }
    if (counter.compare_exchange_weak(used, used + bytes, std::memory_order_acq_rel,
                                      std::memory_order_relaxed)) {
      return true;
    }
  }
}

void AtomicCounterNodeMemoryPool::forceAcquire(std::uint64_t bytes) noexcept {
  std::atomic_ref<std::uint64_t> counter(*m_counter);
  std::uint64_t used = counter.load(std::memory_order_relaxed);
  while (true) {
    const std::uint64_t next = used > std::numeric_limits<std::uint64_t>::max() - bytes
        ? std::numeric_limits<std::uint64_t>::max() : used + bytes;
    if (counter.compare_exchange_weak(used, next, std::memory_order_acq_rel, std::memory_order_relaxed)) {
      return;
    }
  }
}

void AtomicCounterNodeMemoryPool::release(std::uint64_t bytes) noexcept {
  if (bytes == 0U) {
    return;
  }
  std::atomic_ref<std::uint64_t> counter(*m_counter);
  std::uint64_t used = counter.load(std::memory_order_relaxed);
  while (true) {
    const std::uint64_t next = used >= bytes ? used - bytes : 0U;
    if (counter.compare_exchange_weak(used, next, std::memory_order_acq_rel, std::memory_order_relaxed)) {
      return;
    }
  }
}

std::uint64_t AtomicCounterNodeMemoryPool::capacityBytes() const noexcept { return m_capacity_bytes; }

std::uint64_t AtomicCounterNodeMemoryPool::usedBytes() const noexcept {
  return std::atomic_ref<std::uint64_t>(*m_counter).load(std::memory_order_acquire);
}

ProcessLocalNodeMemoryPool::ProcessLocalNodeMemoryPool(std::uint64_t capacity_bytes)
    : m_pool(&m_counter, capacity_bytes) {}

bool ProcessLocalNodeMemoryPool::tryAcquire(std::uint64_t bytes) noexcept { return m_pool.tryAcquire(bytes); }
void ProcessLocalNodeMemoryPool::forceAcquire(std::uint64_t bytes) noexcept { m_pool.forceAcquire(bytes); }
void ProcessLocalNodeMemoryPool::release(std::uint64_t bytes) noexcept { m_pool.release(bytes); }
std::uint64_t ProcessLocalNodeMemoryPool::capacityBytes() const noexcept { return m_pool.capacityBytes(); }
std::uint64_t ProcessLocalNodeMemoryPool::usedBytes() const noexcept { return m_pool.usedBytes(); }

}  // namespace cosmosim::core
