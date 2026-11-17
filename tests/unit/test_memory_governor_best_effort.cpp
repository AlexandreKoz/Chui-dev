// M01 / C04-05 best-effort admission and optional-work reclaim fixture.
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// v10 §5.1/§5.2: optional diagnostics cannot reserve mandatory-physics RAM
// and are deferred first on refusal. reserveBestEffort() is refused outside
// the GREEN band of a finite hard limit; a refused mandatory reserve()
// invokes the registered reclaimer once (registering thread only, never for
// kDiagnostic, never re-entrantly) and retries.
#include <cassert>
#include <cstdint>
#include <optional>
#include <thread>

#include "cosmosim/core/memory_governor.hpp"

namespace {

using cosmosim::core::MemoryAdmissionError;
using cosmosim::core::MemoryClass;
using cosmosim::core::MemoryGovernor;
using cosmosim::core::MemoryGovernorPolicy;
using cosmosim::core::MemoryReservation;

constexpr std::uint64_t k_limit = 1000U;  // amber at 850 bytes (default 8500 bp)

[[nodiscard]] MemoryGovernor makeGovernor() {
  return MemoryGovernor(MemoryGovernorPolicy{.hard_limit_bytes = k_limit});
}

void checkBestEffortBand() {
  MemoryGovernor governor = makeGovernor();
  // Inside the green band: admitted.
  MemoryReservation small = governor.reserveBestEffort(MemoryClass::kDiagnostic, 400U, "test.optional.small");
  // Would reach amber (400 + 500 = 900 >= 850): refused although below the hard limit.
  bool refused = false;
  try {
    (void)governor.reserveBestEffort(MemoryClass::kDiagnostic, 500U, "test.optional.large");
  } catch (const MemoryAdmissionError&) {
    refused = true;
  }
  assert(refused);
  // A mandatory reservation of the same size is still admitted (below the hard limit).
  MemoryReservation mandatory = governor.reserve(MemoryClass::kScratchArena, 500U, "test.mandatory");
  (void)small;
  (void)mandatory;
  // Unlimited policy: best effort behaves like reserve().
  MemoryGovernor unlimited;
  MemoryReservation any = unlimited.reserveBestEffort(MemoryClass::kDiagnostic, 1U << 30U, "test.unlimited");
  (void)any;
}

void checkReclaimBeforeMandatoryRefusal() {
  MemoryGovernor governor = makeGovernor();
  std::optional<MemoryReservation> optional_work;
  optional_work.emplace(governor.reserveBestEffort(MemoryClass::kDiagnostic, 800U, "test.optional"));
  std::uint64_t calls = 0U;
  governor.setBestEffortReclaimer(
      [&] {
        ++calls;
        const bool had = optional_work.has_value();
        optional_work.reset();  // releases the 800 bytes
        return had ? std::uint64_t{800U} : std::uint64_t{0U};
      },
      std::this_thread::get_id());
  // 800 + 600 > 1000: refused, reclaimer releases the optional bytes, retry admitted.
  MemoryReservation mandatory = governor.reserve(MemoryClass::kPhaseResident, 600U, "test.mandatory");
  (void)mandatory;
  assert(calls == 1U);
  assert(!optional_work.has_value());
  assert(governor.bestEffortReclaimCount() == 1U);

  // Nothing left to reclaim: the refusal propagates.
  bool refused = false;
  try {
    (void)governor.reserve(MemoryClass::kPhaseResident, 600U, "test.mandatory.second");
  } catch (const MemoryAdmissionError&) {
    refused = true;
  }
  assert(refused);
  assert(calls == 2U);

  // Diagnostics never trigger the reclaimer.
  refused = false;
  try {
    (void)governor.reserve(MemoryClass::kDiagnostic, 600U, "test.diagnostic");
  } catch (const MemoryAdmissionError&) {
    refused = true;
  }
  assert(refused);
  assert(calls == 2U);

  // Another thread never runs the reclaimer (its owner is single-threaded).
  bool other_refused = false;
  std::thread other([&] {
    try {
      (void)governor.reserve(MemoryClass::kPhaseResident, 600U, "test.other_thread");
    } catch (const MemoryAdmissionError&) {
      other_refused = true;
    }
  });
  other.join();
  assert(other_refused);
  assert(calls == 2U);

  governor.clearBestEffortReclaimer();
}

}  // namespace

int main() {
  checkBestEffortBand();
  checkReclaimBeforeMandatoryRefusal();
  return 0;
}
