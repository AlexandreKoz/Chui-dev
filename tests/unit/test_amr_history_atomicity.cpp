// SOURCE ONLY: NOT RUN — USER-RESERVED.
// Exact history/identity checks target transaction semantics, not convergence.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_ghost_fill.hpp"
#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "../support/test_require.hpp"
#include "../support/scoped_allocation_failure.hpp"

namespace {
namespace amr = cosmosim::amr;
namespace core = cosmosim::core;
constexpr double k_gamma = 1.4;

core::SimulationState makeState() {
  core::SimulationState state;
  state.resizeCells(8U);
  state.resizePatches(4U);
  std::vector<core::GasCellIdentityRecord> identities;
  constexpr std::array<std::uint8_t, 4> levels{0U, 1U, 1U, 2U};
  for (std::uint32_t p = 0U; p < 4U; ++p) {
    const double extent = 1.0 / static_cast<double>(1U << levels[p]);
    amr::writePatchDescriptorToStateRow(state, p, amr::PatchDescriptor{
        .patch_id = p + 1U, .parent_patch_id = 0U, .level = levels[p], .morton_key = p + 1U,
        .origin_comov = {2.0 * p, 0.0, 0.0}, .extent_comov = {extent, extent, extent},
        .cell_dims = {2U, 1U, 1U}});
    state.patches.first_cell[p] = 2U * p;
    state.patches.cell_count[p] = 2U;
    state.patches.owning_rank[p] = 0U;
    for (std::uint32_t c = 0U; c < 2U; ++c) {
      const auto row = 2U * p + c;
      state.cells.patch_index[row] = p;
      state.cells.center_x_comoving[row] = 2.0 * p + (0.25 + 0.5 * c) * extent;
      state.cells.center_y_comoving[row] = 0.5 * extent;
      state.cells.center_z_comoving[row] = 0.5 * extent;
      state.cells.mass_code[row] = 0.5 * extent * extent * extent;
      state.gas_cells.gas_cell_id[row] = 100U + row;
      state.gas_cells.density_code[row] = 1.0;
      state.gas_cells.pressure_code[row] = 1.0;
      state.gas_cells.internal_energy_code[row] = 1.0 / (k_gamma - 1.0);
      identities.push_back({.gas_cell_id = 100U + row, .parent_particle_id = std::nullopt,
                            .owning_patch_id = p + 1U, .local_cell_row = row});
    }
  }
  state.gas_cell_identity.assign(std::move(identities));
  return state;
}

std::uint64_t baseline(const core::SimulationState& state) {
  return core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state));
}

// Encode fields explicitly: structure padding is not a scientific bit target.
std::vector<std::uint64_t> historyBits(const core::SimulationState& state) {
  std::vector<std::uint64_t> result;
  const auto scalar = [&](double v) { result.push_back(std::bit_cast<std::uint64_t>(v)); };
  for (const auto& record : state.amr_temporal_boundary_history.records()) {
    result.insert(result.end(), {record.patch_id, record.patch_level, record.patch_geometry_fingerprint,
        record.gas_cell_identity_generation, record.end_state_valid ? 1U : 0U,
        static_cast<std::uint64_t>(record.cells.size())});
    scalar(record.interval_start_code);
    scalar(record.interval_end_code);
    for (const auto& cell : record.cells) {
      result.push_back(cell.gas_cell_id);
      result.push_back(cell.patch_local_cell);
      for (const double value : {cell.start_mass_density_comoving, cell.start_momentum_density_x_comoving,
          cell.start_momentum_density_y_comoving, cell.start_momentum_density_z_comoving,
          cell.start_total_energy_density_comoving, cell.start_metal_mass_density_comoving,
          cell.end_mass_density_comoving, cell.end_momentum_density_x_comoving,
          cell.end_momentum_density_y_comoving, cell.end_momentum_density_z_comoving,
          cell.end_total_energy_density_comoving, cell.end_metal_mass_density_comoving}) scalar(value);
    }
  }
  return result;
}

struct HistorySnapshot {
  std::vector<std::uint64_t> bits;
  const core::AmrTemporalBoundaryHistoryRecord* records = nullptr;
  const core::AmrTemporalBoundaryHistoryCellRecord* cells = nullptr;
  std::uint64_t bytes = 0U;
  explicit HistorySnapshot(const core::SimulationState& state)
      : bits(historyBits(state)), records(state.amr_temporal_boundary_history.records().data()),
        cells(state.amr_temporal_boundary_history.records().front().cells.data()),
        bytes(state.amr_temporal_boundary_history.ownedCapacityBytes()) {}
  void requireUnchanged(const core::SimulationState& state) const {
    CHUI_REQUIRE(historyBits(state) == bits);
    CHUI_REQUIRE(state.amr_temporal_boundary_history.records().data() == records);
    CHUI_REQUIRE(state.amr_temporal_boundary_history.records().front().cells.data() == cells);
    CHUI_REQUIRE(state.amr_temporal_boundary_history.ownedCapacityBytes() == bytes);
  }
};

void requireNoLeases(const core::MemoryGovernor& governor, std::uint64_t expected_baseline) {
  const auto snapshot = governor.snapshot();
  CHUI_REQUIRE(snapshot.committed_bytes == 0U);
  CHUI_REQUIRE(snapshot.reserved_bytes == 0U);
  CHUI_REQUIRE(snapshot.baseline_owned_bytes == expected_baseline);
}

void testStartRefusalAndLateValidationPreserveParent() {
  auto state = makeState();
  const auto patches = amr::buildProductionAmrPatchDescriptors(state);
  amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 0, 0.0, k_gamma);
  const HistorySnapshot before(state);
  const auto bytes_before = baseline(state);
  core::MemoryGovernor tight({.hard_limit_bytes = bytes_before});
  tight.setBaselineOwnedBytes(bytes_before);
  bool refused = false;
  try {
    amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 1, 0.0, k_gamma, &tight);
  } catch (const core::MemoryAdmissionError&) { refused = true; }
  CHUI_REQUIRE(refused);
  before.requireUnchanged(state);
  requireNoLeases(tight, bytes_before);

  // First new level-1 patch succeeds; second patch's last cell is invalid.
  state.gas_cells.density_code[5] = -1.0;
  core::MemoryGovernor governor;
  governor.setBaselineOwnedBytes(baseline(state));
  bool invalid = false;
  try {
    amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 1, 0.0, k_gamma, &governor);
  } catch (const std::invalid_argument&) { invalid = true; }
  CHUI_REQUIRE(invalid);
  before.requireUnchanged(state);
  requireNoLeases(governor, baseline(state));
  state.gas_cells.density_code[5] = 1.0;
  auto duplicate = patches;
  duplicate.push_back(patches[1]);
  bool duplicate_rejected = false;
  try {
    amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, duplicate, 1, 0.0, k_gamma, &governor);
  } catch (const std::runtime_error&) { duplicate_rejected = true; }
  CHUI_REQUIRE(duplicate_rejected);
  before.requireUnchanged(state);
  requireNoLeases(governor, baseline(state));
}

void testEndFailurePreservesEveryLevel() {
  auto state = makeState();
  const auto patches = amr::buildProductionAmrPatchDescriptors(state);
  amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 0, 0.0, k_gamma);
  amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 1, 0.0, k_gamma);
  const HistorySnapshot before(state);
  core::MemoryGovernor governor;
  governor.setBaselineOwnedBytes(baseline(state));
  state.gas_cells.density_code[5] = -1.0;
  bool failed = false;
  try {
    amr::captureAmrTemporalBoundaryHistoryEndForLevel(state, patches, 1, 0.5, k_gamma, &governor);
  } catch (const std::invalid_argument&) { failed = true; }
  CHUI_REQUIRE(failed);
  before.requireUnchanged(state);
  requireNoLeases(governor, baseline(state));
  state.gas_cells.density_code[5] = 1.0;
  auto changed_geometry = patches;
  changed_geometry[2].origin_comov[0] += 0.1;
  failed = false;
  try {
    amr::captureAmrTemporalBoundaryHistoryEndForLevel(state, changed_geometry, 1, 0.5, k_gamma, &governor);
  } catch (const std::runtime_error&) { failed = true; }
  CHUI_REQUIRE(failed);
  before.requireUnchanged(state);
  core::MemoryGovernor tight({.hard_limit_bytes = baseline(state)});
  tight.setBaselineOwnedBytes(baseline(state));
  bool refused = false;
  try {
    amr::captureAmrTemporalBoundaryHistoryEndForLevel(state, patches, 1, 0.5, k_gamma, &tight);
  } catch (const core::MemoryAdmissionError&) { refused = true; }
  CHUI_REQUIRE(refused);
  before.requireUnchanged(state);
  requireNoLeases(tight, baseline(state));
}

void testEveryStartAndEndAllocationFailure() {
  for (const bool end_capture : {false, true}) {
    bool reached_success = false;
    std::size_t failure_count = 0U;
    for (std::ptrdiff_t allocation = 0; allocation < 2048; ++allocation) {
      auto state = makeState();
      const auto patches = amr::buildProductionAmrPatchDescriptors(state);
      amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 0, 0.0, k_gamma);
      if (end_capture) amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 1, 0.0, k_gamma);
      const HistorySnapshot before(state);
      const auto bytes_before = baseline(state);
      core::MemoryGovernor governor;
      governor.setBaselineOwnedBytes(bytes_before);
      bool failed = false;
      {
        cosmosim::test_support::ScopedAllocationFailure inject(allocation);
        try {
          if (end_capture) {
            amr::captureAmrTemporalBoundaryHistoryEndForLevel(state, patches, 1, 0.5, k_gamma, &governor);
          } else {
            amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 1, 0.0, k_gamma, &governor);
          }
        } catch (const std::bad_alloc&) { failed = true; }
      }
      if (failed) {
        ++failure_count;
        before.requireUnchanged(state);
        requireNoLeases(governor, bytes_before);
      } else {
        reached_success = true;
        CHUI_REQUIRE(state.amr_temporal_boundary_history.size() == 3U);
        CHUI_REQUIRE(!state.amr_temporal_boundary_history.findByPatchId(1U)->end_state_valid);
        CHUI_REQUIRE(state.amr_temporal_boundary_history.findByPatchId(2U)->end_state_valid == end_capture);
        CHUI_REQUIRE(state.amr_temporal_boundary_history.findByPatchId(3U)->end_state_valid == end_capture);
        requireNoLeases(governor, baseline(state));
        break;
      }
    }
    CHUI_REQUIRE(reached_success);
    CHUI_REQUIRE(failure_count > 1U);  // Exercise late nested allocations too.
  }
}

void testRetirementNeedsNoAllocationAndReconciles() {
  auto state = makeState();
  const auto patches = amr::buildProductionAmrPatchDescriptors(state);
  core::MemoryGovernor governor;
  governor.setBaselineOwnedBytes(baseline(state));
  amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 0, 0.0, k_gamma, &governor);
  amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 1, 0.0, k_gamma, &governor);
  // Un-governed API itself requires no heap storage, even when every new
  // allocation would fail. Governed retirement additionally inspects reports.
  {
    cosmosim::test_support::ScopedAllocationFailure inject(0);
    amr::retireAmrTemporalBoundaryHistoryLevel(state, 1);
  }
  CHUI_REQUIRE(state.amr_temporal_boundary_history.size() == 1U);
  governor.setBaselineOwnedBytes(baseline(state));
  amr::captureAmrTemporalBoundaryHistoryStartForLevel(state, patches, 1, 0.0, k_gamma, &governor);
  amr::retireAmrTemporalBoundaryHistoryLevel(state, 1, &governor);
  CHUI_REQUIRE(state.amr_temporal_boundary_history.size() == 1U);
  requireNoLeases(governor, baseline(state));
  amr::retireAmrTemporalBoundaryHistory(state, &governor);
  CHUI_REQUIRE(state.amr_temporal_boundary_history.empty());
  requireNoLeases(governor, baseline(state));
}
}  // namespace

int main() {
  testStartRefusalAndLateValidationPreserveParent();
  testEndFailurePreservesEveryLevel();
  testEveryStartAndEndAllocationFailure();
  testRetirementNeedsNoAllocationAndReconciles();
}
