#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/node_memory_admission.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <utility>

namespace cosmosim::core {
namespace {

constexpr std::uint64_t k_basis_point_denominator = 10000U;

[[nodiscard]] bool validMemoryClass(MemoryClass memory_class) noexcept {
  return memoryClassIndex(memory_class) <
      static_cast<std::size_t>(MemoryClass::kCount);
}

[[nodiscard]] std::uint64_t scaledBasisPointThreshold(
    std::uint64_t hard_limit_bytes,
    std::uint16_t basis_points) noexcept {
  const std::uint64_t whole = hard_limit_bytes / k_basis_point_denominator;
  const std::uint64_t remainder = hard_limit_bytes % k_basis_point_denominator;
  return whole * static_cast<std::uint64_t>(basis_points) +
      (remainder * static_cast<std::uint64_t>(basis_points)) /
          k_basis_point_denominator;
}

[[nodiscard]] std::pair<std::uint64_t, std::uint64_t> pressureThresholds(
    std::uint64_t hard_limit_bytes,
    std::uint16_t amber_basis_points,
    std::uint16_t red_basis_points) noexcept {
  if (hard_limit_bytes <= 1U) {
    return {1U, 1U};
  }
  if (hard_limit_bytes == 2U) {
    return {1U, 2U};
  }
  const std::uint64_t amber_raw = scaledBasisPointThreshold(
      hard_limit_bytes, amber_basis_points);
  const std::uint64_t amber = std::clamp<std::uint64_t>(
      amber_raw, 1U, hard_limit_bytes - 2U);
  const std::uint64_t red_raw = scaledBasisPointThreshold(
      hard_limit_bytes, red_basis_points);
  const std::uint64_t red = std::clamp<std::uint64_t>(
      red_raw, amber + 1U, hard_limit_bytes - 1U);
  return {amber, red};
}

[[nodiscard]] std::string rejectionDiagnostic(
    MemoryClass memory_class,
    std::string_view owner,
    std::uint64_t requested_bytes,
    const MemoryGovernorSnapshot& snapshot,
    std::string_view reason) {
  std::ostringstream out;
  out << "memory reservation rejected"
      << " owner=" << owner
      << " class=" << memoryClassLabel(memory_class)
      << " requested_bytes=" << requested_bytes
      << " committed_bytes=" << snapshot.committed_bytes
      << " reserved_bytes=" << snapshot.reserved_bytes
      << " baseline_owned_bytes=" << snapshot.baseline_owned_bytes
      << " external_runtime_reserve_bytes=" << snapshot.external_runtime_reserve_bytes
      << " planned_overlap_reserve_bytes=" << snapshot.planned_overlap_reserve_bytes
      << " accounted_bytes=" << snapshot.accounted_bytes
      << " hard_limit_bytes=" << snapshot.hard_limit_bytes
      << " headroom_bytes=" << snapshot.headroom_bytes
      << " pressure=" << memoryPressureLabel(snapshot.pressure)
      << " reason=" << reason;
  return out.str();
}

}  // namespace


DeterministicBatchSizingResult selectDeterministicBatchSize(
    const MemoryGovernorSnapshot& snapshot,
    const DeterministicBatchSizingPolicy& policy) {
  if (policy.requested_max_items == 0U || policy.bytes_per_item == 0U ||
      policy.minimum_items == 0U || policy.alignment_items == 0U ||
      policy.headroom_use_basis_points == 0U ||
      policy.headroom_use_basis_points > k_basis_point_denominator) {
    throw std::invalid_argument("deterministic batch sizing policy is invalid");
  }
  if (policy.minimum_items > policy.requested_max_items) {
    throw std::invalid_argument(
        "deterministic batch minimum_items exceeds requested_max_items");
  }

  DeterministicBatchSizingResult result;
  result.selected_items = policy.requested_max_items;
  if (snapshot.headroom_bytes == std::numeric_limits<std::uint64_t>::max()) {
    result.usable_headroom_bytes = snapshot.headroom_bytes;
  }

  if (snapshot.headroom_bytes == std::numeric_limits<std::uint64_t>::max()) {
    if (policy.requested_max_items >
        std::numeric_limits<std::uint64_t>::max() / policy.bytes_per_item) {
      throw std::overflow_error("deterministic batch selected byte count overflow");
    }
    result.selected_bytes = policy.requested_max_items * policy.bytes_per_item;
    return result;
  }

  const std::uint64_t after_fixed =
      snapshot.headroom_bytes > policy.fixed_reserve_bytes
          ? snapshot.headroom_bytes - policy.fixed_reserve_bytes
          : 0U;
  const std::uint64_t whole = after_fixed / k_basis_point_denominator;
  const std::uint64_t remainder = after_fixed % k_basis_point_denominator;
  result.usable_headroom_bytes =
      whole * static_cast<std::uint64_t>(policy.headroom_use_basis_points) +
      (remainder * static_cast<std::uint64_t>(policy.headroom_use_basis_points)) /
          k_basis_point_denominator;

  std::uint64_t by_headroom = result.usable_headroom_bytes / policy.bytes_per_item;
  if (by_headroom >= policy.alignment_items) {
    by_headroom -= by_headroom % policy.alignment_items;
  }
  if (by_headroom < policy.minimum_items) {
    result.selected_items = 0U;
    result.selected_bytes = 0U;
    result.constrained_by_headroom = true;
    return result;
  }
  result.selected_items = std::min(policy.requested_max_items, by_headroom);
  if (result.selected_items >
      std::numeric_limits<std::uint64_t>::max() / policy.bytes_per_item) {
    throw std::overflow_error("deterministic batch selected byte count overflow");
  }
  result.selected_bytes = result.selected_items * policy.bytes_per_item;
  result.constrained_by_headroom = result.selected_items < policy.requested_max_items;
  return result;
}

DeterministicBatchSizingResult selectDeterministicBatchSizeForWorkspaces(
    const MemoryGovernorSnapshot& snapshot,
    const DeterministicBatchSizingPolicy& policy,
    std::span<const DeterministicBatchWorkspace> workspaces) {
  if (policy.requested_max_items == 0U || policy.minimum_items == 0U ||
      policy.minimum_items > policy.requested_max_items ||
      policy.alignment_items == 0U || policy.headroom_use_basis_points == 0U ||
      policy.headroom_use_basis_points > k_basis_point_denominator || workspaces.empty()) {
    throw std::invalid_argument("deterministic workspace batch policy is invalid");
  }
  // Validate representability independently of available headroom. This
  // prevents a huge requested count from wrapping into an apparently cheap
  // allocation and makes unlimited mode use the same checked contract.
  for (const auto& workspace : workspaces) {
    if (workspace.bytes_per_item != 0U && policy.requested_max_items >
        std::numeric_limits<std::uint64_t>::max() / workspace.bytes_per_item) {
      throw std::overflow_error("deterministic workspace batch byte count overflows uint64");
    }
  }
  const auto requiredBytes = [&](std::uint64_t count) {
    std::uint64_t bytes = 0U;
    for (const auto& workspace : workspaces) {
      if (count <= workspace.retained_capacity_items) { continue; }
      const std::uint64_t additional = count * workspace.bytes_per_item;
      if (additional > std::numeric_limits<std::uint64_t>::max() - bytes) {
        return std::optional<std::uint64_t>{};
      }
      bytes += additional;
    }
    return std::optional<std::uint64_t>{bytes};
  };
  // Reuse the established headroom and basis-point policy, including its
  // fixed-reserve subtraction. The dummy unit width only obtains the usable
  // byte budget; actual workspace costs are evaluated below.
  const auto budget = selectDeterministicBatchSize(
      snapshot, DeterministicBatchSizingPolicy{
          .requested_max_items = 1U,
          .bytes_per_item = 1U,
          .fixed_reserve_bytes = policy.fixed_reserve_bytes,
          .minimum_items = 1U,
          .alignment_items = 1U,
          .headroom_use_basis_points = policy.headroom_use_basis_points,
      });
  const std::uint64_t usable = budget.usable_headroom_bytes;
  const auto fits = [&](std::uint64_t count) {
    const auto bytes = requiredBytes(count);
    return bytes.has_value() && *bytes <= usable;
  };
  if (snapshot.headroom_bytes == std::numeric_limits<std::uint64_t>::max() &&
      !requiredBytes(policy.requested_max_items).has_value()) {
    throw std::overflow_error("deterministic workspace aggregate byte count overflows uint64");
  }
  DeterministicBatchSizingResult result;
  result.usable_headroom_bytes = usable;
  if (!fits(policy.minimum_items)) {
    result.constrained_by_headroom = true;
    return result;
  }
  std::uint64_t low = policy.minimum_items;
  std::uint64_t high = policy.requested_max_items;
  while (low < high) {
    const std::uint64_t mid = low + (high - low) / 2U + (high - low) % 2U;
    if (fits(mid)) { low = mid; } else { high = mid - 1U; }
  }
  if (low >= policy.alignment_items) {
    low -= low % policy.alignment_items;
  }
  if (low < policy.minimum_items) {
    result.constrained_by_headroom = true;
    return result;
  }
  result.selected_items = low;
  result.selected_bytes = *requiredBytes(low);
  result.constrained_by_headroom = low < policy.requested_max_items;
  return result;
}

std::size_t memoryClassIndex(MemoryClass memory_class) noexcept {
  return static_cast<std::size_t>(memory_class);
}

std::string_view memoryClassLabel(MemoryClass memory_class) noexcept {
  switch (memory_class) {
    case MemoryClass::kCanonicalPersistent:
      return "canonical_persistent";
    case MemoryClass::kPersistentCache:
      return "persistent_cache";
    case MemoryClass::kPhaseResident:
      return "phase_resident";
    case MemoryClass::kScratchArena:
      return "scratch_arena";
    case MemoryClass::kCommunication:
      return "communication";
    case MemoryClass::kDiagnostic:
      return "diagnostic";
    case MemoryClass::kExternalRuntime:
      return "external_runtime";
    case MemoryClass::kCount:
      return "invalid";
  }
  return "invalid";
}

std::string_view memoryPressureLabel(MemoryPressure pressure) noexcept {
  switch (pressure) {
    case MemoryPressure::kGreen:
      return "green";
    case MemoryPressure::kAmber:
      return "amber";
    case MemoryPressure::kRed:
      return "red";
    case MemoryPressure::kTrip:
      return "trip";
  }
  return "trip";
}

std::uint64_t checkedMemoryBytesAdd(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 byte addition overflow");
  }
  return lhs + rhs;
}

MemoryReservation::MemoryReservation(
    MemoryGovernor* governor,
    MemoryClass memory_class,
    std::uint64_t bytes,
    std::uint16_t owner_slot) noexcept
    : m_governor(governor),
      m_memory_class(memory_class),
      m_bytes(bytes),
      m_state(State::kPending),
      m_owner_slot(owner_slot) {}

MemoryReservation::~MemoryReservation() noexcept {
  release();
}

MemoryReservation::MemoryReservation(MemoryReservation&& other) noexcept
    : m_governor(std::exchange(other.m_governor, nullptr)),
      m_memory_class(other.m_memory_class),
      m_bytes(std::exchange(other.m_bytes, 0U)),
      m_state(std::exchange(other.m_state, State::kInert)),
      m_owner_slot(std::exchange(other.m_owner_slot, std::uint16_t{0U})) {}

MemoryReservation& MemoryReservation::operator=(MemoryReservation&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  release();
  m_governor = std::exchange(other.m_governor, nullptr);
  m_memory_class = other.m_memory_class;
  m_bytes = std::exchange(other.m_bytes, 0U);
  m_state = std::exchange(other.m_state, State::kInert);
  m_owner_slot = std::exchange(other.m_owner_slot, std::uint16_t{0U});
  return *this;
}

bool MemoryReservation::valid() const noexcept {
  return m_state == State::kPending || m_state == State::kCommitted;
}

bool MemoryReservation::pending() const noexcept {
  return m_state == State::kPending;
}

bool MemoryReservation::committed() const noexcept {
  return m_state == State::kCommitted;
}

std::uint64_t MemoryReservation::bytes() const noexcept {
  return m_bytes;
}

MemoryClass MemoryReservation::memoryClass() const noexcept {
  return m_memory_class;
}

void MemoryReservation::commit() {
  if (m_state != State::kPending || m_governor == nullptr) {
    throw std::logic_error("MemoryReservation.commit requires one live pending reservation");
  }
  m_governor->commitReservation(m_memory_class, m_bytes);
  m_state = State::kCommitted;
}

void MemoryReservation::reconcileBaselineOwnedAndRelease(
    std::uint64_t baseline_owned_bytes) {
  if (m_state != State::kCommitted || m_governor == nullptr) {
    throw std::logic_error(
        "MemoryReservation.reconcileBaselineOwnedAndRelease requires one live committed reservation");
  }
  m_governor->reconcileBaselineAndReleaseCommittedReservation(
      m_memory_class, m_bytes, baseline_owned_bytes, m_owner_slot);
  m_state = State::kReleased;
  m_governor = nullptr;
  m_bytes = 0U;
}

void MemoryReservation::release() noexcept {
  if (m_governor == nullptr || m_state == State::kInert || m_state == State::kReleased) {
    m_state = State::kReleased;
    m_governor = nullptr;
    m_bytes = 0U;
    return;
  }
  if (m_state == State::kPending) {
    m_governor->releasePendingReservation(m_memory_class, m_bytes, m_owner_slot);
  } else if (m_state == State::kCommitted) {
    m_governor->releaseCommittedReservation(m_memory_class, m_bytes, m_owner_slot);
  }
  m_state = State::kReleased;
  m_governor = nullptr;
  m_bytes = 0U;
}

MemoryGovernor::MemoryGovernor(MemoryGovernorPolicy policy)
    : m_policy(policy) {
  if (!std::isfinite(m_policy.safety_margin_fraction) ||
      m_policy.safety_margin_fraction < 0.0 ||
      m_policy.safety_margin_fraction > 1.0) {
    throw std::invalid_argument(
        "MemoryGovernor safety margin must be finite and within [0,1]");
  }
  if (m_policy.amber_basis_points == 0U ||
      m_policy.amber_basis_points >= m_policy.red_basis_points ||
      m_policy.red_basis_points >= k_basis_point_denominator) {
    throw std::invalid_argument(
        "MemoryGovernor pressure thresholds require 0 < amber < red < 10000 basis points");
  }
  updatePeakAccountedLocked();
}

MemoryGovernor::~MemoryGovernor() {
  detachNodePool();
}

MemoryReservation MemoryGovernor::reserve(
    MemoryClass memory_class,
    std::uint64_t requested_bytes,
    std::string_view owner) {
  if (!validMemoryClass(memory_class)) {
    throw std::invalid_argument("MemoryGovernor.reserve received invalid memory class");
  }
  if (owner.empty()) {
    throw std::invalid_argument("MemoryGovernor.reserve requires a non-empty owner label");
  }

  std::lock_guard lock(m_mutex);
  const std::uint64_t current_raw = rawAccountedBytesLocked();
  std::uint64_t prospective_raw = current_raw;
  try {
    prospective_raw = checkedMemoryBytesAdd(
        current_raw, requested_bytes, "MemoryGovernor.reserve prospective demand");
  } catch (const std::overflow_error&) {
    ++m_rejection_count;
    ++m_owner_slots[ownerSlotLocked(owner)].rejection_count;
    MemoryGovernorSnapshot current;
    current.hard_limit_bytes = m_policy.hard_limit_bytes;
    current.baseline_owned_bytes = m_baseline_owned_bytes;
    current.external_runtime_reserve_bytes = m_policy.external_runtime_reserve_bytes;
    current.planned_overlap_reserve_bytes = m_policy.planned_overlap_reserve_bytes;
    current.committed_bytes = m_committed_bytes;
    current.reserved_bytes = m_reserved_bytes;
    current.accounted_bytes = current_raw;
    current.safety_adjusted_accounted_bytes = safetyAdjustedBytes(current_raw);
    current.headroom_bytes = headroomBytesLocked(current_raw);
    current.peak_committed_bytes = m_peak_committed_bytes;
    current.peak_reserved_bytes = m_peak_reserved_bytes;
    current.peak_accounted_bytes = m_peak_accounted_bytes;
    current.rejection_count = m_rejection_count;
    current.pressure = MemoryPressure::kTrip;
    current.committed_by_class = m_committed_by_class;
    current.reserved_by_class = m_reserved_by_class;
    throw std::overflow_error(rejectionDiagnostic(
        memory_class,
        owner,
        requested_bytes,
        current,
        "uint64_accounting_overflow"));
  }

  const std::uint64_t prospective_adjusted = safetyAdjustedBytes(prospective_raw);
  if (m_policy.hard_limit_bytes != 0U &&
      prospective_adjusted > m_policy.hard_limit_bytes) {
    ++m_rejection_count;
    ++m_owner_slots[ownerSlotLocked(owner)].rejection_count;
    MemoryGovernorSnapshot current;
    current.hard_limit_bytes = m_policy.hard_limit_bytes;
    current.baseline_owned_bytes = m_baseline_owned_bytes;
    current.external_runtime_reserve_bytes = m_policy.external_runtime_reserve_bytes;
    current.planned_overlap_reserve_bytes = m_policy.planned_overlap_reserve_bytes;
    current.committed_bytes = m_committed_bytes;
    current.reserved_bytes = m_reserved_bytes;
    current.accounted_bytes = current_raw;
    current.safety_adjusted_accounted_bytes = safetyAdjustedBytes(current_raw);
    current.headroom_bytes = headroomBytesLocked(current_raw);
    current.peak_committed_bytes = m_peak_committed_bytes;
    current.peak_reserved_bytes = m_peak_reserved_bytes;
    current.peak_accounted_bytes = m_peak_accounted_bytes;
    current.rejection_count = m_rejection_count;
    current.pressure = MemoryPressure::kTrip;
    current.committed_by_class = m_committed_by_class;
    current.reserved_by_class = m_reserved_by_class;
    throw MemoryAdmissionError(rejectionDiagnostic(
        memory_class,
        owner,
        requested_bytes,
        current,
        "hard_limit_would_be_exceeded"));
  }

  if (m_node_pool != nullptr && !m_node_pool->tryAcquire(requested_bytes)) {
    ++m_rejection_count;
    ++m_owner_slots[ownerSlotLocked(owner)].rejection_count;
    MemoryGovernorSnapshot current;
    current.hard_limit_bytes = m_policy.hard_limit_bytes;
    current.committed_bytes = m_committed_bytes;
    current.reserved_bytes = m_reserved_bytes;
    current.accounted_bytes = current_raw;
    current.headroom_bytes = headroomBytesLocked(current_raw);
    current.rejection_count = m_rejection_count;
    current.pressure = MemoryPressure::kTrip;
    throw MemoryAdmissionError(rejectionDiagnostic(
        memory_class, owner, requested_bytes, current,
        "node_pool_exhausted used=" + std::to_string(m_node_pool->usedBytes()) +
            " capacity=" + std::to_string(m_node_pool->capacityBytes())));
  }
  m_node_pool_held_bytes += m_node_pool != nullptr ? requested_bytes : 0U;
  const std::size_t class_index = memoryClassIndex(memory_class);
  m_reserved_bytes = checkedMemoryBytesAdd(
      m_reserved_bytes, requested_bytes, "MemoryGovernor.reserve reserved total");
  m_reserved_by_class[class_index] = checkedMemoryBytesAdd(
      m_reserved_by_class[class_index], requested_bytes,
      "MemoryGovernor.reserve class reserved total");
  m_peak_reserved_bytes = std::max(m_peak_reserved_bytes, m_reserved_bytes);
  updatePeakAccountedLocked();
  const std::uint16_t owner_slot = ownerSlotLocked(owner);
  OwnerSlot& slot = m_owner_slots[owner_slot];
  slot.live_bytes = slot.live_bytes > std::numeric_limits<std::uint64_t>::max() - requested_bytes
      ? std::numeric_limits<std::uint64_t>::max()
      : slot.live_bytes + requested_bytes;
  slot.peak_live_bytes = std::max(slot.peak_live_bytes, slot.live_bytes);
  ++slot.reservation_count;
  return MemoryReservation(this, memory_class, requested_bytes, owner_slot);
}

std::uint16_t MemoryGovernor::ownerSlotLocked(std::string_view owner) noexcept {
  const std::size_t label_length = std::min(owner.size(), k_memory_owner_label_capacity);
  const std::string_view key = owner.substr(0U, label_length);
  constexpr std::size_t k_other_slot = k_memory_owner_telemetry_slots - 1U;
  for (std::size_t index = 0; index < k_other_slot; ++index) {
    OwnerSlot& slot = m_owner_slots[index];
    if (!slot.used) {
      slot.used = true;
      std::copy(key.begin(), key.end(), slot.label.begin());
      slot.label[label_length] = '\0';
      return static_cast<std::uint16_t>(index);
    }
    if (std::string_view(slot.label.data()) == key) {
      return static_cast<std::uint16_t>(index);
    }
  }
  OwnerSlot& other = m_owner_slots[k_other_slot];
  if (!other.used) {
    other.used = true;
    constexpr std::string_view k_other_label = "<other>";
    std::copy(k_other_label.begin(), k_other_label.end(), other.label.begin());
    other.label[k_other_label.size()] = '\0';
  }
  return static_cast<std::uint16_t>(k_other_slot);
}

void MemoryGovernor::attachNodePool(NodeMemoryPool* pool) {
  std::lock_guard lock(m_mutex);
  if (m_node_pool != nullptr) {
    throw std::logic_error("MemoryGovernor already has a node pool attached");
  }
  if (pool == nullptr) {
    return;
  }
  const std::uint64_t accounted = rawAccountedBytesLocked();
  if (!pool->tryAcquire(accounted)) {
    ++m_rejection_count;
    throw MemoryAdmissionError(
        "memory node pool cannot hold this rank's accounted bytes at attach: accounted=" +
        std::to_string(accounted) + " used=" + std::to_string(pool->usedBytes()) +
        " capacity=" + std::to_string(pool->capacityBytes()));
  }
  m_node_pool = pool;
  m_node_pool_held_bytes = accounted;
}

void MemoryGovernor::detachNodePool() noexcept {
  std::lock_guard lock(m_mutex);
  if (m_node_pool != nullptr) {
    m_node_pool->release(m_node_pool_held_bytes);
  }
  m_node_pool = nullptr;
  m_node_pool_held_bytes = 0U;
}

std::uint64_t MemoryGovernor::nodePoolHeldBytes() const {
  std::lock_guard lock(m_mutex);
  return m_node_pool_held_bytes;
}

void MemoryGovernor::releaseNodePoolBytesLocked(std::uint64_t bytes) noexcept {
  if (m_node_pool == nullptr || bytes == 0U) {
    return;
  }
  const std::uint64_t released = std::min(bytes, m_node_pool_held_bytes);
  m_node_pool->release(released);
  m_node_pool_held_bytes -= released;
}

void MemoryGovernor::releaseOwnerBytesLocked(std::uint16_t owner_slot, std::uint64_t bytes) noexcept {
  if (owner_slot >= m_owner_slots.size()) {
    return;
  }
  OwnerSlot& slot = m_owner_slots[owner_slot];
  slot.live_bytes = slot.live_bytes >= bytes ? slot.live_bytes - bytes : 0U;
}

std::vector<MemoryOwnerTelemetry> MemoryGovernor::ownerTelemetry() const {
  std::lock_guard lock(m_mutex);
  std::vector<MemoryOwnerTelemetry> result;
  for (const OwnerSlot& slot : m_owner_slots) {
    if (!slot.used) {
      continue;
    }
    result.push_back(MemoryOwnerTelemetry{
        .label = std::string(slot.label.data()),
        .live_bytes = slot.live_bytes,
        .peak_live_bytes = slot.peak_live_bytes,
        .reservation_count = slot.reservation_count,
        .rejection_count = slot.rejection_count});
  }
  return result;
}

void MemoryGovernor::setBaselineOwnedBytes(std::uint64_t baseline_owned_bytes) {
  std::lock_guard lock(m_mutex);
  // Validate the replacement baseline before mutating state so a reconciliation
  // overflow cannot leave the governor with an unusable partially-applied value.
  std::uint64_t prospective = baseline_owned_bytes;
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.external_runtime_reserve_bytes,
      "MemoryGovernor baseline external runtime reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.planned_overlap_reserve_bytes,
      "MemoryGovernor baseline planned overlap reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, m_committed_bytes, "MemoryGovernor baseline committed bytes");
  prospective = checkedMemoryBytesAdd(
      prospective, m_reserved_bytes, "MemoryGovernor baseline reserved bytes");
  adjustNodePoolForBaselineLocked(m_baseline_owned_bytes, baseline_owned_bytes);
  m_baseline_owned_bytes = baseline_owned_bytes;
  m_peak_accounted_bytes = std::max(m_peak_accounted_bytes, prospective);
}

void MemoryGovernor::adjustNodePoolForBaselineLocked(
    std::uint64_t old_baseline, std::uint64_t new_baseline) noexcept {
  if (m_node_pool == nullptr) {
    return;
  }
  if (new_baseline > old_baseline) {
    // Baseline is already-resident memory: record it truthfully.
    const std::uint64_t delta = new_baseline - old_baseline;
    m_node_pool->forceAcquire(delta);
    m_node_pool_held_bytes = m_node_pool_held_bytes > std::numeric_limits<std::uint64_t>::max() - delta
        ? std::numeric_limits<std::uint64_t>::max() : m_node_pool_held_bytes + delta;
  } else {
    releaseNodePoolBytesLocked(old_baseline - new_baseline);
  }
}

MemoryGovernorSnapshot MemoryGovernor::snapshot() const {
  std::lock_guard lock(m_mutex);
  const std::uint64_t raw = rawAccountedBytesLocked();
  const std::uint64_t adjusted = safetyAdjustedBytes(raw);
  MemoryGovernorSnapshot result;
  result.hard_limit_bytes = m_policy.hard_limit_bytes;
  result.baseline_owned_bytes = m_baseline_owned_bytes;
  result.external_runtime_reserve_bytes = m_policy.external_runtime_reserve_bytes;
  result.planned_overlap_reserve_bytes = m_policy.planned_overlap_reserve_bytes;
  result.committed_bytes = m_committed_bytes;
  result.reserved_bytes = m_reserved_bytes;
  result.accounted_bytes = raw;
  result.safety_adjusted_accounted_bytes = adjusted;
  result.headroom_bytes = headroomBytesLocked(raw);
  result.peak_committed_bytes = m_peak_committed_bytes;
  result.peak_reserved_bytes = m_peak_reserved_bytes;
  result.peak_accounted_bytes = m_peak_accounted_bytes;
  result.rejection_count = m_rejection_count;
  result.pressure = pressureForAdjustedBytes(adjusted);
  result.committed_by_class = m_committed_by_class;
  result.reserved_by_class = m_reserved_by_class;
  return result;
}

const MemoryGovernorPolicy& MemoryGovernor::policy() const noexcept {
  return m_policy;
}

void MemoryGovernor::commitReservation(MemoryClass memory_class, std::uint64_t bytes) {
  std::lock_guard lock(m_mutex);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_reserved_bytes < bytes || m_reserved_by_class[class_index] < bytes) {
    throw std::logic_error("MemoryGovernor reservation accounting underflow during commit");
  }
  m_reserved_bytes -= bytes;
  m_reserved_by_class[class_index] -= bytes;
  m_committed_bytes = checkedMemoryBytesAdd(
      m_committed_bytes, bytes, "MemoryGovernor.commit committed total");
  m_committed_by_class[class_index] = checkedMemoryBytesAdd(
      m_committed_by_class[class_index], bytes,
      "MemoryGovernor.commit class committed total");
  m_peak_committed_bytes = std::max(m_peak_committed_bytes, m_committed_bytes);
  updatePeakAccountedLocked();
}

void MemoryGovernor::releasePendingReservation(
    MemoryClass memory_class,
    std::uint64_t bytes,
    std::uint16_t owner_slot) noexcept {
  std::lock_guard lock(m_mutex);
  releaseOwnerBytesLocked(owner_slot, bytes);
  releaseNodePoolBytesLocked(bytes);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_reserved_bytes < bytes || m_reserved_by_class[class_index] < bytes) {
    std::terminate();
  }
  m_reserved_bytes -= bytes;
  m_reserved_by_class[class_index] -= bytes;
}

void MemoryGovernor::releaseCommittedReservation(
    MemoryClass memory_class,
    std::uint64_t bytes,
    std::uint16_t owner_slot) noexcept {
  std::lock_guard lock(m_mutex);
  releaseOwnerBytesLocked(owner_slot, bytes);
  releaseNodePoolBytesLocked(bytes);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_committed_bytes < bytes || m_committed_by_class[class_index] < bytes) {
    std::terminate();
  }
  m_committed_bytes -= bytes;
  m_committed_by_class[class_index] -= bytes;
}

void MemoryGovernor::reconcileBaselineAndReleaseCommittedReservation(
    MemoryClass memory_class,
    std::uint64_t bytes,
    std::uint64_t baseline_owned_bytes,
    std::uint16_t owner_slot) {
  std::lock_guard lock(m_mutex);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_committed_bytes < bytes || m_committed_by_class[class_index] < bytes) {
    throw std::logic_error(
        "MemoryGovernor committed accounting underflow during baseline reconciliation");
  }

  const std::uint64_t remaining_committed = m_committed_bytes - bytes;
  std::uint64_t prospective = baseline_owned_bytes;
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.external_runtime_reserve_bytes,
      "MemoryGovernor reconciled baseline external runtime reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.planned_overlap_reserve_bytes,
      "MemoryGovernor reconciled baseline planned overlap reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, remaining_committed,
      "MemoryGovernor reconciled baseline remaining committed bytes");
  prospective = checkedMemoryBytesAdd(
      prospective, m_reserved_bytes,
      "MemoryGovernor reconciled baseline reserved bytes");

  m_committed_bytes = remaining_committed;
  m_committed_by_class[class_index] -= bytes;
  releaseNodePoolBytesLocked(bytes);
  adjustNodePoolForBaselineLocked(m_baseline_owned_bytes, baseline_owned_bytes);
  m_baseline_owned_bytes = baseline_owned_bytes;
  releaseOwnerBytesLocked(owner_slot, bytes);
  m_peak_accounted_bytes = std::max(m_peak_accounted_bytes, prospective);
}

std::uint64_t MemoryGovernor::rawAccountedBytesLocked() const {
  std::uint64_t total = m_baseline_owned_bytes;
  total = checkedMemoryBytesAdd(
      total, m_policy.external_runtime_reserve_bytes,
      "MemoryGovernor raw external runtime reserve");
  total = checkedMemoryBytesAdd(
      total, m_policy.planned_overlap_reserve_bytes,
      "MemoryGovernor raw planned overlap reserve");
  total = checkedMemoryBytesAdd(total, m_committed_bytes, "MemoryGovernor raw committed bytes");
  total = checkedMemoryBytesAdd(total, m_reserved_bytes, "MemoryGovernor raw reserved bytes");
  return total;
}

std::uint64_t MemoryGovernor::safetyAdjustedBytes(std::uint64_t raw_bytes) const {
  if (raw_bytes == 0U || m_policy.safety_margin_fraction == 0.0) {
    return raw_bytes;
  }
  const long double scaled = static_cast<long double>(raw_bytes) *
      (1.0L + static_cast<long double>(m_policy.safety_margin_fraction));
  if (!std::isfinite(scaled) ||
      scaled > static_cast<long double>(std::numeric_limits<std::uint64_t>::max())) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  return static_cast<std::uint64_t>(std::ceil(scaled));
}

std::uint64_t MemoryGovernor::headroomBytesLocked(std::uint64_t raw_bytes) const {
  if (m_policy.hard_limit_bytes == 0U) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  if (safetyAdjustedBytes(raw_bytes) > m_policy.hard_limit_bytes) {
    return 0U;
  }
  if (m_policy.safety_margin_fraction == 0.0) {
    return m_policy.hard_limit_bytes - raw_bytes;
  }

  // Find the greatest raw demand whose safety-adjusted demand still fits the
  // configured process ceiling. This uses exactly the same rounding rule as
  // admission rather than relying on a lossy algebraic inversion.
  std::uint64_t low = raw_bytes;
  std::uint64_t high = m_policy.hard_limit_bytes;
  while (low < high) {
    const std::uint64_t delta = high - low;
    const std::uint64_t mid = low + delta / 2U + delta % 2U;
    if (safetyAdjustedBytes(mid) <= m_policy.hard_limit_bytes) {
      low = mid;
    } else {
      high = mid - 1U;
    }
  }
  return low - raw_bytes;
}

MemoryPressure MemoryGovernor::pressureForAdjustedBytes(
    std::uint64_t adjusted_bytes) const noexcept {
  if (m_policy.hard_limit_bytes == 0U) {
    return MemoryPressure::kGreen;
  }
  if (adjusted_bytes > m_policy.hard_limit_bytes) {
    return MemoryPressure::kTrip;
  }
  if (adjusted_bytes == 0U) {
    return MemoryPressure::kGreen;
  }
  const auto [amber, red] = pressureThresholds(
      m_policy.hard_limit_bytes,
      m_policy.amber_basis_points,
      m_policy.red_basis_points);
  if (adjusted_bytes >= red) {
    return MemoryPressure::kRed;
  }
  if (adjusted_bytes >= amber) {
    return MemoryPressure::kAmber;
  }
  return MemoryPressure::kGreen;
}

void MemoryGovernor::updatePeakAccountedLocked() {
  const std::uint64_t raw = rawAccountedBytesLocked();
  m_peak_accounted_bytes = std::max(m_peak_accounted_bytes, raw);
}

}  // namespace cosmosim::core
