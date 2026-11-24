// CHUI-MHD-PERF-v1 resource gates in source: MP02 cap, MP06 rejection before
// mutation, MP07 no MHD buffers when disabled, MP11 finite governor, batch
// sizing and whole-host envelope. Source written, NOT EXECUTED.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "mhd_test_fixtures.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/mhd/mhd_memory.hpp"

namespace {

using namespace cosmosim::mhd;
using cosmosim::mhd::test::GasPrimitive;
using cosmosim::mhd::test::TestLevel;

void expectMemoryError(auto&& action) {
  bool rejected = false;
  try {
    action();
  } catch (const MhdError& error) {
    rejected = error.code() == MhdErrorCode::kMemoryBudgetExceeded;
  }
  CHUI_REQUIRE(rejected);
}

// MP06: an insufficient governor refuses the MHD workspace before any stage
// allocation or state mutation; the governor holds nothing afterwards.
void testAdmissionRejectsBeforeMutation() {
  TestLevel level({16, 16, 16}, {16, 16, 16}, {1.0 / 16, 1.0 / 16, 1.0 / 16},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {0.1, 0.0, 0.0}, {}, 1.0,
                   k_mhd_gamma_default);
  const double rho = level.cellValue(0, {3, 3, 3});
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  cosmosim::core::MemoryGovernorPolicy policy{};
  policy.hard_limit_bytes = 1024U * 1024U;  // far below one 16^3 workspace set
  cosmosim::core::MemoryGovernor governor(policy);
  expectMemoryError([&]() { stepper.admitMemory(governor, {}); });
  CHUI_REQUIRE(!stepper.memoryAdmitted());
  CHUI_REQUIRE(governor.snapshot().reserved_bytes == 0U);
  CHUI_REQUIRE(stepper.ownedCapacityBytes() == 0U);  // nothing was allocated
  CHUI_REQUIRE(level.cellValue(0, {3, 3, 3}) == rho);
}

// Admitted capacity bounds the actual retained capacity; a mode that needs
// unadmitted workspaces is refused before work starts.
void testAdmittedCapacityAndModeGuard() {
  TestLevel level({16, 16, 16}, {16, 16, 16}, {1.0 / 16, 1.0 / 16, 1.0 / 16},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {0.1, 0.0, 0.0}, {}, 1.0,
                   k_mhd_gamma_default);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  cosmosim::core::MemoryGovernorPolicy policy{};
  policy.hard_limit_bytes = 512ULL * 1024ULL * 1024ULL;
  cosmosim::core::MemoryGovernor governor(policy);
  MhdResourceLimits limits{};
  limits.scratch_bytes_per_active_cell = 600U;  // 16^3 PPM is ~484 B/cell
  stepper.admitMemory(governor, {}, limits);
  CHUI_REQUIRE(stepper.memoryAdmitted());
  const MhdLevelAdvanceResult result =
      stepper.advance(cosmosim::mhd::test::isolatedInterval(0.0, 1.0e-3), {});
  CHUI_REQUIRE(result.accepted);
  CHUI_REQUIRE(governor.snapshot().reserved_bytes + governor.snapshot().committed_bytes >=
               result.scratch_capacity_bytes);
  // Plateau (MP05 in miniature): a second identical interval does not grow
  // the retained capacity.
  const std::uint64_t capacity = stepper.ownedCapacityBytes();
  CHUI_REQUIRE(stepper.advance(cosmosim::mhd::test::isolatedInterval(1.0e-3, 1.0e-3), {}).accepted);
  CHUI_REQUIRE(stepper.ownedCapacityBytes() == capacity);
  // Edge-integral lanes were not admitted.
  MhdLevelAdvanceOptions integrals{};
  integrals.produce_edge_integrals = true;
  expectMemoryError([&]() {
    (void)stepper.advance(cosmosim::mhd::test::isolatedInterval(2.0e-3, 1.0e-3), integrals);
  });
}

// MP02: the default 512 B/cell cap rejects a configuration that exceeds it
// (16^3 PPM + RKL2 lanes) rather than silently allocating.
void testScratchCap() {
  MhdMemoryRequest request{};
  request.geometry.nx = request.geometry.ny = request.geometry.nz = 16;
  request.geometry.dx_comoving = request.geometry.dy_comoving = request.geometry.dz_comoving = 1.0;
  request.nonideal = true;
  request.super_stepping = true;
  CHUI_REQUIRE(scratchBytesPerActiveCell(request, 1U) > 512.0);
  const MhdResourceLimits limits{};
  CHUI_REQUIRE(selectMhdBatchPatches(request, 1ULL << 40U, limits) == 0U);
  request.geometry.nx = request.geometry.ny = request.geometry.nz = 32;
  CHUI_REQUIRE(selectMhdBatchPatches(request, 1ULL << 40U, limits) == 8U);
  const std::uint64_t per_patch = estimatePatchMemory(request).phasePlusScratchBytes();
  CHUI_REQUIRE(selectMhdBatchPatches(request, 3U * per_patch + 1U, limits) == 3U);
}

// MP04/MP11: whole-host envelope over ranks and the finite-limit rule.
void testHostEnvelope() {
  const MhdResourceLimits limits{};
  const std::uint64_t gib = 1024ULL * 1024ULL * 1024ULL;
  const std::array<std::uint64_t, 4> fits{10 * gib, 10 * gib, 10 * gib, 10 * gib};
  requireHostEnvelope(fits, limits, true);
  const std::array<std::uint64_t, 4> over{11 * gib, 11 * gib, 11 * gib, 11 * gib};
  expectMemoryError([&]() { requireHostEnvelope(over, limits, true); });
  const std::array<std::uint64_t, 2> unlimited{0U, 10 * gib};
  expectMemoryError([&]() { requireHostEnvelope(unlimited, limits, true); });
  requireHostEnvelope(unlimited, limits, false);  // legacy accounting allowed unqualified
}

// MP07: an MHD-disabled SimulationState carries no magnetic lanes and the
// memory report has no MHD entries.
void testDisabledHasNoMhdMemory() {
  const cosmosim::core::SimulationState state{};
  CHUI_REQUIRE(state.magnetic_faces.empty());
  const cosmosim::core::MemoryReport report = cosmosim::core::collectSimulationMemoryReport(state);
  for (const cosmosim::core::MemoryEntry& entry : report.entries) {
    CHUI_REQUIRE(entry.label.rfind("mhd.", 0) != 0U);
  }
}

}  // namespace

int main() {
  testAdmissionRejectsBeforeMutation();
  testAdmittedCapacityAndModeGuard();
  testScratchCap();
  testHostEnvelope();
  testDisabledHasNoMhdMemory();
  return 0;
}
