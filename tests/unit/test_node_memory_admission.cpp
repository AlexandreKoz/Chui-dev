// Node-aware CPU memory admission and bounded owner telemetry.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants (synthetic ceilings only; no physical machine is assumed):
//  * parsers: /proc/meminfo kB, cgroup "max"/v1 sentinel -> unlimited,
//    SLURM MiB/G suffixes;
//  * node limit = min of known sources; reserve subtracted exactly once;
//    rank share = floor(usable / ranks); sum of rank shares <= usable for
//    8/16/48/128 GiB nodes and 1..64 co-located ranks;
//  * a smaller process budget binds; 0 keeps the node share;
//  * invalid inputs (0 ranks, reserve >= limit, no limit) are refused;
//  * MemoryGovernor: reserve -> commit -> release returns accounting to the
//    baseline exactly once, uint64 overflow is an explicit error, refusal does
//    not change accounting, per-owner peaks are bounded and recorded.
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/node_memory_admission.hpp"
#include "../support/test_require.hpp"

namespace {

namespace co = cosmosim::core;
constexpr std::uint64_t k_gib = 1024ULL * 1024ULL * 1024ULL;

void testParsers() {
  CHUI_REQUIRE(co::parseProcMeminfoTotalBytes("MemTotal:       49152000 kB\nMemFree: 1 kB\n").value() ==
               49152000ULL * 1024ULL);
  CHUI_REQUIRE(!co::parseProcMeminfoTotalBytes("MemFree: 1 kB\n").has_value());
  CHUI_REQUIRE(!co::parseCgroupMemoryLimitBytes("max\n").has_value());
  CHUI_REQUIRE(!co::parseCgroupMemoryLimitBytes("9223372036854771712").has_value());
  CHUI_REQUIRE(co::parseCgroupMemoryLimitBytes("17179869184\n").value() == 16ULL * k_gib);
  CHUI_REQUIRE(co::parseSchedulerMemPerNodeMiB("4096").value() == 4ULL * k_gib);
  CHUI_REQUIRE(co::parseSchedulerMemPerNodeMiB("8G").value() == 8ULL * k_gib);
  CHUI_REQUIRE(!co::parseSchedulerMemPerNodeMiB("abc").has_value());
}

void testSyntheticNodesNeverOvercommit() {
  for (std::uint64_t node_gib : {8ULL, 16ULL, 48ULL, 128ULL}) {
    for (std::uint32_t ranks : {1U, 2U, 3U, 7U, 16U, 64U}) {
      co::NodeMemoryAdmissionInput input;
      input.discovery.physical_ram_bytes = node_gib * k_gib;
      input.node_reserve_bytes = k_gib / 2U;
      input.ranks_on_node = ranks;
      const auto result = co::computeNodeMemoryAdmission(input);
      CHUI_REQUIRE(result.node_limit_source == co::NodeMemoryLimitSource::kPhysicalRam);
      CHUI_REQUIRE(result.node_usable_bytes == node_gib * k_gib - k_gib / 2U);
      // floor division => ranks * share <= usable (no co-located overcommit).
      CHUI_REQUIRE(result.rank_share_bytes <= result.node_usable_bytes / ranks);
      CHUI_REQUIRE(result.node_usable_bytes - result.rank_share_bytes * ranks < ranks);
      CHUI_REQUIRE(result.rank_hard_limit_bytes == result.rank_share_bytes);
    }
  }
}

void testMinimumOfKnownSourcesAndReserveOnce() {
  co::NodeMemoryAdmissionInput input;
  input.discovery.physical_ram_bytes = 48ULL * k_gib;
  input.discovery.cgroup_limit_bytes = 42ULL * k_gib;
  input.discovery.cgroup_source = co::NodeMemoryLimitSource::kCgroupV2;
  input.explicit_node_budget_bytes = 36ULL * k_gib;  // illustrative historical nominal target
  input.node_reserve_bytes = 4ULL * k_gib;
  input.ranks_on_node = 4U;
  const auto result = co::computeNodeMemoryAdmission(input);
  CHUI_REQUIRE(result.node_limit_source == co::NodeMemoryLimitSource::kExplicitNodeBudget);
  CHUI_REQUIRE(result.node_limit_bytes == 36ULL * k_gib);
  CHUI_REQUIRE(result.node_usable_bytes == 32ULL * k_gib);
  CHUI_REQUIRE(result.rank_share_bytes == 8ULL * k_gib);
  input.explicit_node_budget_bytes = 0U;
  const auto cgroup = co::computeNodeMemoryAdmission(input);
  CHUI_REQUIRE(cgroup.node_limit_source == co::NodeMemoryLimitSource::kCgroupV2);
  CHUI_REQUIRE(cgroup.rank_share_bytes == (42ULL - 4ULL) * k_gib / 4ULL);
}

void testProcessBudgetBinding() {
  co::NodeMemoryAdmissionInput input;
  input.discovery.physical_ram_bytes = 128ULL * k_gib;
  input.ranks_on_node = 2U;
  input.process_budget_bytes = 10ULL * k_gib;
  const auto result = co::computeNodeMemoryAdmission(input);
  CHUI_REQUIRE(result.process_budget_was_binding);
  CHUI_REQUIRE(result.rank_hard_limit_bytes == 10ULL * k_gib);
  input.process_budget_bytes = 0U;
  const auto unlimited_process = co::computeNodeMemoryAdmission(input);
  CHUI_REQUIRE(!unlimited_process.process_budget_was_binding);
  CHUI_REQUIRE(unlimited_process.rank_hard_limit_bytes == 64ULL * k_gib);
}

template <class TException, class TCallable>
bool throws(TCallable&& callable) {
  try {
    callable();
  } catch (const TException&) {
    return true;
  }
  return false;
}

void testInvalidInputsRefused() {
  co::NodeMemoryAdmissionInput input;
  input.discovery.physical_ram_bytes = 8ULL * k_gib;
  input.ranks_on_node = 0U;
  CHUI_REQUIRE(throws<std::invalid_argument>([&] { (void)co::computeNodeMemoryAdmission(input); }));
  input.ranks_on_node = 1U;
  input.node_reserve_bytes = 8ULL * k_gib;
  CHUI_REQUIRE(throws<std::invalid_argument>([&] { (void)co::computeNodeMemoryAdmission(input); }));
  co::NodeMemoryAdmissionInput unknown;
  CHUI_REQUIRE(throws<std::runtime_error>([&] { (void)co::computeNodeMemoryAdmission(unknown); }));
}

void testGovernorLifecycleOverflowAndOwnerTelemetry() {
  co::MemoryGovernor governor(co::MemoryGovernorPolicy{.hard_limit_bytes = 1000U});
  {
    auto first = governor.reserve(co::MemoryClass::kScratchArena, 400U, "hydro.ctu_transverse");
    first.commit();
    auto second = governor.reserve(co::MemoryClass::kScratchArena, 500U, "hydro.ctu_transverse");
    CHUI_REQUIRE(governor.snapshot().accounted_bytes == 900U);
    // Refusal leaves accounting unchanged and counts the owner rejection.
    CHUI_REQUIRE(throws<co::MemoryAdmissionError>(
        [&] { (void)governor.reserve(co::MemoryClass::kScratchArena, 200U, "sources.feedback"); }));
    CHUI_REQUIRE(governor.snapshot().accounted_bytes == 900U);
    second.release();
    second.release();  // idempotent: released exactly once
    CHUI_REQUIRE(governor.snapshot().accounted_bytes == 400U);
  }
  CHUI_REQUIRE(governor.snapshot().accounted_bytes == 0U);
  bool found_ctu = false;
  bool found_feedback = false;
  for (const auto& owner : governor.ownerTelemetry()) {
    if (owner.label == "hydro.ctu_transverse") {
      found_ctu = true;
      CHUI_REQUIRE(owner.peak_live_bytes == 900U);
      CHUI_REQUIRE(owner.live_bytes == 0U);
      CHUI_REQUIRE(owner.reservation_count == 2U);
    }
    if (owner.label == "sources.feedback") {
      found_feedback = true;
      CHUI_REQUIRE(owner.rejection_count == 1U);
    }
  }
  CHUI_REQUIRE(found_ctu && found_feedback);

  co::MemoryGovernor unlimited;
  auto big = unlimited.reserve(co::MemoryClass::kScratchArena, std::numeric_limits<std::uint64_t>::max() - 1U, "big");
  CHUI_REQUIRE(throws<std::overflow_error>(
      [&] { (void)unlimited.reserve(co::MemoryClass::kScratchArena, 2U, "overflow"); }));
  big.release();
}

// shared_pool: two co-located "ranks" (governors) share one node pool.
// Uneven use is admitted, the node-wide sum is bounded, releases return
// capacity, baseline reconciliation records resident memory truthfully, and
// detach/destruction return every held byte.
void testNodeSharedPoolBoundsTheSumNotTheShare() {
  co::ProcessLocalNodeMemoryPool pool(1000U);
  {
    co::MemoryGovernor heavy;
    co::MemoryGovernor light;
    heavy.attachNodePool(&pool);
    light.attachNodePool(&pool);
    auto big = heavy.reserve(co::MemoryClass::kScratchArena, 800U, "rank0.amr");  // > equal share
    auto small = light.reserve(co::MemoryClass::kScratchArena, 150U, "rank1.hydro");
    CHUI_REQUIRE(pool.usedBytes() == 950U);
    CHUI_REQUIRE(throws<co::MemoryAdmissionError>(
        [&] { (void)light.reserve(co::MemoryClass::kScratchArena, 100U, "rank1.hydro"); }));
    CHUI_REQUIRE(pool.usedBytes() == 950U);  // refusal holds nothing
    big.release();
    CHUI_REQUIRE(pool.usedBytes() == 150U);
    auto now_fits = light.reserve(co::MemoryClass::kScratchArena, 100U, "rank1.hydro");
    CHUI_REQUIRE(pool.usedBytes() == 250U);
    light.setBaselineOwnedBytes(900U);  // resident memory is recorded even above capacity
    CHUI_REQUIRE(pool.usedBytes() == 1150U);
    CHUI_REQUIRE(throws<co::MemoryAdmissionError>(
        [&] { (void)heavy.reserve(co::MemoryClass::kScratchArena, 1U, "rank0.amr"); }));
    light.setBaselineOwnedBytes(0U);
    CHUI_REQUIRE(pool.usedBytes() == 250U);
    now_fits.release();
    small.release();
    CHUI_REQUIRE(pool.usedBytes() == 0U);
  }
  CHUI_REQUIRE(pool.usedBytes() == 0U);
}

void testOwnerTableIsBounded() {
  co::MemoryGovernor governor;
  for (std::size_t i = 0; i < 2U * co::k_memory_owner_telemetry_slots; ++i) {
    auto reservation = governor.reserve(co::MemoryClass::kDiagnostic, 1U, "owner." + std::to_string(i));
    (void)reservation;
  }
  const auto owners = governor.ownerTelemetry();
  CHUI_REQUIRE(owners.size() == co::k_memory_owner_telemetry_slots);
  CHUI_REQUIRE(owners.back().label == "<other>");
  CHUI_REQUIRE(owners.back().reservation_count == co::k_memory_owner_telemetry_slots + 1U);
}

}  // namespace

int main() {
  testParsers();
  testSyntheticNodesNeverOvercommit();
  testMinimumOfKnownSourcesAndReserveOnce();
  testProcessBudgetBinding();
  testInvalidInputsRefused();
  testGovernorLifecycleOverflowAndOwnerTelemetry();
  testOwnerTableIsBounded();
  testNodeSharedPoolBoundsTheSumNotTheShare();
  return 0;
}
