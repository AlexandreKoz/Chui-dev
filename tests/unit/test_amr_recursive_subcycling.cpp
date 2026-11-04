// CHUI Part 3 Milestone I — recursive (>= 3 level) AMR subcycling.
// STATUS: written as source only. NOT RUN — USER-RESERVED.
//
// Hierarchy (leaf-only, 2:1 balanced INCLUDING the periodic wrap, box
// [0,3) x [0,1) x [0,1)): roots A, B, C (2^3 cells each) at x = 0, 1, 2;
// B and C refined into level-1 octets; the B child at origin (1.5, 0, 0)
// refined again into a level-2 octet (its neighbours are all level 1).
// A further test checks that a non-2:1 hierarchy is refused.
// Invariants exercised through the production driver
// advanceProductionAmrHydroRecursive():
//  1. Clock structure: substeps_by_level == {1, 2, 4} for one coarse step
//     (level L advances dt / 2^L, i.e. different physical step sizes).
//  2. Uniform state (rho, v, P constant) stays uniform to 1e-12 relative at
//     every level (temporal ghosts, deep ghosts and reflux introduce no
//     spurious gradient).
//  3. Conservation: with a density/pressure perturbation, total mass and
//     x-momentum over all leaves are conserved to 1e-12 relative (periodic
//     box; coarse/fine reflux at both nested interfaces).
//  4. Synchronization: after the step no pending flux register and no open
//     temporal interval remain (legal regrid/output/checkpoint boundary).
//  5. Start-boundary refusal: a recursive step attempted with an open
//     temporal interval is refused before mutation (illegal mid-cycle call).
//  6. Single-level equivalence: without refinement the recursive driver
//     equals the synchronized reference sweep to 1e-14 relative.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "../support/test_require.hpp"
#include "../support/scoped_allocation_failure.hpp"

namespace {

namespace amr = cosmosim::amr;
namespace core = cosmosim::core;
namespace hy = cosmosim::hydro;

constexpr double k_gamma = 5.0 / 3.0;

core::SimulationState makeRoots(bool perturbed) {
  core::SimulationState state;
  state.resizeCells(24U);
  state.resizePatches(3U);
  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t p = 0; p < 3U; ++p) {
    amr::writePatchDescriptorToStateRow(state, p, amr::PatchDescriptor{
        .patch_id = 1U + p, .parent_patch_id = 0U, .level = 0U, .morton_key = 1U + p,
        .origin_comov = {static_cast<double>(p), 0.0, 0.0}, .extent_comov = {1.0, 1.0, 1.0},
        .cell_dims = {2U, 2U, 2U}});
    state.patches.first_cell[p] = 8U * p;
    state.patches.cell_count[p] = 8U;
    state.patches.owning_rank[p] = 0;
    for (std::uint32_t c = 0; c < 8U; ++c) {
      const std::uint32_t row = 8U * p + c;
      const double x = static_cast<double>(p) + 0.25 + 0.5 * static_cast<double>(c % 2U);
      const double rho = perturbed ? 1.0 + 0.3 * std::sin(3.0 * x) : 1.0;
      const double pressure = perturbed ? 1.0 + 0.2 * std::cos(2.0 * x) : 1.0;
      state.cells.center_x_comoving[row] = x;
      state.cells.center_y_comoving[row] = 0.25 + 0.5 * static_cast<double>((c / 2U) % 2U);
      state.cells.center_z_comoving[row] = 0.25 + 0.5 * static_cast<double>(c / 4U);
      state.cells.patch_index[row] = p;
      state.cells.mass_code[row] = rho * 0.125;
      state.gas_cells.gas_cell_id[row] = 100U + row;
      state.gas_cells.density_code[row] = rho;
      state.gas_cells.pressure_code[row] = pressure;
      state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
      state.gas_cells.velocity_x_peculiar[row] = 0.2;
      state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
      state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
      records.push_back(core::GasCellIdentityRecord{.gas_cell_id = 100U + row, .parent_particle_id = std::nullopt,
                                                    .owning_patch_id = 1U + p, .local_cell_row = row});
    }
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

amr::ProductionAmrHydroOptions baseOptions(const core::SimulationState& state, core::MemoryGovernor& governor) {
  governor.setBaselineOwnedBytes(core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state)));
  amr::ProductionAmrHydroOptions options;
  options.physical_boundary_kind = hy::HydroBoundaryKind::kPeriodic;
  options.periodic_box_extent_comoving = {3.0, 1.0, 1.0};
  options.adiabatic_index = k_gamma;
  options.density_floor = 1.0e-12;
  options.pressure_floor = 1.0e-12;
  options.refinement_ratio = 2U;
  options.regrid_memory_governor = &governor;
  return options;
}

core::SimulationState makeThreeLevels(bool perturbed) {
  core::SimulationState state = makeRoots(perturbed);
  core::MemoryGovernor governor;
  const auto options = baseOptions(state, governor);
  const auto roots = amr::buildProductionAmrPatchDescriptors(state);
  const auto b = std::find_if(roots.begin(), roots.end(), [](const amr::PatchDescriptor& d) { return d.patch_id == 2U; });
  (void)amr::refineProductionPatchInSimulationState(state, *b, 10U, 1000U, options);
  {
    core::MemoryGovernor governor_c;
    const auto options_c = baseOptions(state, governor_c);
    const auto after_b = amr::buildProductionAmrPatchDescriptors(state);
    const auto c = std::find_if(after_b.begin(), after_b.end(), [](const amr::PatchDescriptor& d) { return d.patch_id == 3U; });
    CHUI_REQUIRE(c != after_b.end());
    (void)amr::refineProductionPatchInSimulationState(state, *c, 20U, 2000U, options_c);
  }
  core::MemoryGovernor governor2;
  const auto options2 = baseOptions(state, governor2);
  const auto level1 = amr::buildProductionAmrPatchDescriptors(state);
  const auto corner = std::find_if(level1.begin(), level1.end(), [](const amr::PatchDescriptor& d) {
    return d.level == 1U && std::abs(d.origin_comov[0] - 1.5) < 1.0e-12 && d.origin_comov[1] == 0.0 &&
        d.origin_comov[2] == 0.0;
  });
  CHUI_REQUIRE(corner != level1.end());
  (void)amr::refineProductionPatchInSimulationState(state, *corner, 30U, 3000U, options2);
  return state;
}

struct Totals {
  double mass = 0.0;
  double momentum_x = 0.0;
};

Totals totals(const core::SimulationState& state) {
  Totals t;
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    const auto patch = amr::patchDescriptorFromStateRow(state, state.cells.patch_index[row]);
    const double volume = patch.extent_comov[0] * patch.extent_comov[1] * patch.extent_comov[2] /
        static_cast<double>(patch.cell_dims[0] * patch.cell_dims[1] * patch.cell_dims[2]);
    const double mass = state.gas_cells.density_code[row] * volume;
    t.mass += mass;
    t.momentum_x += mass * state.gas_cells.velocity_x_peculiar[row];
  }
  return t;
}

amr::ProductionAmrHydroDiagnostics step(core::SimulationState& state, bool recursive, double dt) {
  core::MemoryGovernor governor;
  auto options = baseOptions(state, governor);
  const hy::HydroCoreSolver solver(k_gamma);
  hy::HllcRiemannSolver riemann;
  const hy::HydroUpdateContext update{.dt_code = dt, .time_begin_code = 0.0, .time_end_code = dt,
                                      .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const hy::HydroSourceContext source_context{.update = update};
  return recursive
      ? amr::advanceProductionAmrHydroRecursive(state, {}, update, source_context, solver, riemann, {}, options)
      : amr::advanceProductionAmrHydro(state, {}, update, source_context, solver, riemann, {}, options);
}

void testClockStructureAndUniformState() {
  core::SimulationState state = makeThreeLevels(false);
  const auto diagnostics = step(state, true, 0.02);
  CHUI_REQUIRE(diagnostics.substeps_by_level.size() == 3U);
  CHUI_REQUIRE(diagnostics.substeps_by_level[0] == 1U);
  CHUI_REQUIRE(diagnostics.substeps_by_level[1] == 2U);
  CHUI_REQUIRE(diagnostics.substeps_by_level[2] == 4U);
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    CHUI_REQUIRE(std::abs(state.gas_cells.density_code[row] - 1.0) < 1.0e-12);
    CHUI_REQUIRE(std::abs(state.gas_cells.velocity_x_peculiar[row] - 0.2) < 1.0e-12);
  }
  CHUI_REQUIRE(state.pending_flux_registers.empty());
  CHUI_REQUIRE(state.amr_temporal_boundary_history.empty());
}

void testConservationAcrossNestedInterfaces() {
  core::SimulationState state = makeThreeLevels(true);
  const Totals before = totals(state);
  (void)step(state, true, 0.02);
  const Totals after = totals(state);
  CHUI_REQUIRE(std::abs(after.mass - before.mass) <= 1.0e-12 * before.mass);
  CHUI_REQUIRE(std::abs(after.momentum_x - before.momentum_x) <= 1.0e-12 * std::max(1.0, std::abs(before.momentum_x)));
  CHUI_REQUIRE(state.pending_flux_registers.empty());
  CHUI_REQUIRE(state.amr_temporal_boundary_history.empty());
}

void testIllegalMidCycleStartRefused() {
  core::SimulationState state = makeThreeLevels(true);
  amr::captureAmrTemporalBoundaryHistoryStartForLevel(
      state, amr::buildProductionAmrPatchDescriptors(state), 0, 0.0, k_gamma);
  const std::vector<double> rho_before(state.gas_cells.density_code.begin(), state.gas_cells.density_code.end());
  bool refused = false;
  try {
    (void)step(state, true, 0.02);
  } catch (const std::runtime_error&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
  CHUI_REQUIRE(std::equal(rho_before.begin(), rho_before.end(), state.gas_cells.density_code.begin()));
}

void testUnbalancedHierarchyRefused() {
  // Refine only B, then its child at x = 1.0 (touching A at level 0): a
  // level-2 patch touches level 0 -> refused before any mutation.
  core::SimulationState state = makeRoots(true);
  core::MemoryGovernor governor;
  const auto options = baseOptions(state, governor);
  const auto roots = amr::buildProductionAmrPatchDescriptors(state);
  const auto b = std::find_if(roots.begin(), roots.end(), [](const amr::PatchDescriptor& d) { return d.patch_id == 2U; });
  (void)amr::refineProductionPatchInSimulationState(state, *b, 10U, 1000U, options);
  core::MemoryGovernor governor2;
  const auto options2 = baseOptions(state, governor2);
  const auto level1 = amr::buildProductionAmrPatchDescriptors(state);
  const auto touching = std::find_if(level1.begin(), level1.end(), [](const amr::PatchDescriptor& d) {
    return d.level == 1U && std::abs(d.origin_comov[0] - 1.0) < 1.0e-12 && d.origin_comov[1] == 0.0 &&
        d.origin_comov[2] == 0.0;
  });
  CHUI_REQUIRE(touching != level1.end());
  (void)amr::refineProductionPatchInSimulationState(state, *touching, 30U, 3000U, options2);
  const std::vector<double> rho_before(state.gas_cells.density_code.begin(), state.gas_cells.density_code.end());
  bool refused = false;
  try {
    (void)step(state, true, 0.02);
  } catch (const std::runtime_error&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
  CHUI_REQUIRE(std::equal(rho_before.begin(), rho_before.end(), state.gas_cells.density_code.begin()));
}

void testSingleLevelEqualsSynchronizedReference() {
  core::SimulationState recursive = makeRoots(true);
  core::SimulationState reference = makeRoots(true);
  (void)step(recursive, true, 0.02);
  (void)step(reference, false, 0.02);
  for (std::size_t row = 0; row < recursive.cells.size(); ++row) {
    const double a = recursive.gas_cells.density_code[row];
    const double b = reference.gas_cells.density_code[row];
    CHUI_REQUIRE(std::abs(a - b) <= 1.0e-14 * std::max(1.0, std::abs(b)));
  }
}

void testRecursiveWorkingSetRefusalAndAllocationFailure() {
  for (const bool refined : {false, true}) {
    auto state = refined ? makeThreeLevels(true) : makeRoots(true);
    const auto before = state;
    const auto baseline = core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state));
    const hy::HydroCoreSolver solver(k_gamma);
    hy::HllcRiemannSolver riemann;
    const hy::HydroUpdateContext update{.dt_code = 0.02, .time_begin_code = 0.0, .time_end_code = 0.02,
                                        .scale_factor = 1.0, .hubble_rate_code = 0.0};
    const hy::HydroSourceContext source{.update = update};
    const auto unchanged = [&] {
      CHUI_REQUIRE(std::memcmp(before.cells.mass_code.data(), state.cells.mass_code.data(),
          state.cells.size() * sizeof(double)) == 0);
      for (const auto lanes : {
          std::pair{&before.gas_cells.density_code, &state.gas_cells.density_code},
          std::pair{&before.gas_cells.internal_energy_code, &state.gas_cells.internal_energy_code},
          std::pair{&before.gas_cells.pressure_code, &state.gas_cells.pressure_code},
          std::pair{&before.gas_cells.velocity_x_peculiar, &state.gas_cells.velocity_x_peculiar},
          std::pair{&before.gas_cells.metal_mass_code, &state.gas_cells.metal_mass_code}}) {
        CHUI_REQUIRE(std::memcmp(lanes.first->data(), lanes.second->data(),
            lanes.first->size() * sizeof(double)) == 0);
      }
      CHUI_REQUIRE(state.amr_temporal_boundary_history.empty());
      CHUI_REQUIRE(state.pending_flux_registers.empty());
      CHUI_REQUIRE(state.gasCellIdentityGeneration() == before.gasCellIdentityGeneration());
      CHUI_REQUIRE(state.gravitySourceGeneration() == before.gravitySourceGeneration());
    };
    core::MemoryGovernor tight({.hard_limit_bytes = baseline});
    auto options = baseOptions(state, tight);
    bool refused = false;
    try {
      (void)amr::advanceProductionAmrHydroRecursive(state, {}, update, source, solver, riemann, {}, options);
    } catch (const core::MemoryAdmissionError&) { refused = true; }
    CHUI_REQUIRE(refused);
    unchanged();
    CHUI_REQUIRE(tight.snapshot().committed_bytes == 0U);
    CHUI_REQUIRE(tight.snapshot().reserved_bytes == 0U);
    CHUI_REQUIRE(tight.snapshot().baseline_owned_bytes == baseline);

    core::MemoryGovernor governor;
    options = baseOptions(state, governor);
    bool failed = false;
    std::size_t selected_rows = 0U;
    for (std::size_t p = 0U; p < state.patches.size(); ++p) {
      if (state.patches.level[p] == (refined ? 1 : 0)) selected_rows += state.patches.cell_count[p];
    }
    {
      // Select one exact level-row allocation, after driver admission.
      cosmosim::test_support::ScopedAllocationFailure inject(0, selected_rows * sizeof(std::uint32_t));
      try {
        (void)amr::advanceProductionAmrHydroRecursive(state, {}, update, source, solver, riemann, {}, options);
      } catch (const std::bad_alloc&) { failed = true; }
    }
    CHUI_REQUIRE(failed);
    unchanged();
    CHUI_REQUIRE(governor.snapshot().peak_committed_bytes > 0U);
    CHUI_REQUIRE(governor.snapshot().committed_bytes == 0U);
    CHUI_REQUIRE(governor.snapshot().reserved_bytes == 0U);
    CHUI_REQUIRE(governor.snapshot().baseline_owned_bytes == baseline);
  }
}

void testCoverageRejectsOverlapAndMissingRows() {
  auto state = makeRoots(false);
  CHUI_REQUIRE(amr::hasProductionAmrHydroCoverage(state));
  state.patches.first_cell[1] = 0U;
  CHUI_REQUIRE(!amr::hasProductionAmrHydroCoverage(state));
  state = makeRoots(false);
  state.patches.cell_count[1] = 0U;
  CHUI_REQUIRE(!amr::hasProductionAmrHydroCoverage(state));
}

}  // namespace

int main() {
  testRecursiveWorkingSetRefusalAndAllocationFailure();
  testCoverageRejectsOverlapAndMissingRows();
  testClockStructureAndUniformState();
  testConservationAcrossNestedInterfaces();
  testIllegalMidCycleStartRefused();
  testUnbalancedHierarchyRefused();
  testSingleLevelEqualsSynchronizedReference();
  return 0;
}
