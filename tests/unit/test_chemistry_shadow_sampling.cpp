// C03/C04 chemistry shadow sampling fixture.
// NOT COMPILED OR RUN in this handoff (the owner reserves build/test stages).
//
// The sampled shadow disagreement steers the adaptive chemistry controller,
// so its sample set must not depend on row order or MPI decomposition:
//   * selection is decided by the (gas_cell_id, accepted step) key alone;
//     forward and reverse row order give identical samples and maxima;
//   * the per-stage cost cap never silently reshapes the set: rows cut by it
//     are counted (cap_truncated);
//   * the controller acts only when the interval-GLOBAL key-selected count
//     fits the budget (decomposition independent), checked on a 1-way vs
//     2-way split where the cap binds in one layout only;
//   * failed references are counted (reference_failures), never dropped.
// Exit code 77 (registered as SKIP_RETURN_CODE) when the detailed reference
// backend cannot be constructed in this build (external data required).
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <exception>
#include <iostream>
#include <numeric>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"

namespace {

using namespace cosmosim;
namespace h = cosmosim::hydro;
namespace p = cosmosim::physics;

constexpr std::size_t k_rows = 64U;

[[nodiscard]] p::ThermochemistryUnitBridge units() {
  p::ThermochemistryUnitBridge u;
  u.density_code_to_cgs = 1.6726e-24;
  u.velocity_code_to_cm_s = 1.0e5;
  u.length_code_to_cm = 3.0857e21;
  u.mass_code_to_g = 1.989e33;
  return u;
}

[[nodiscard]] core::SimulationState gas(const p::CompositionLayout& layout) {
  core::SimulationState state;
  state.resizeCells(k_rows);
  for (std::size_t r = 0; r < k_rows; ++r) {
    const double rho = 1.0 + 0.01 * static_cast<double>(r);
    state.gas_cells.gas_cell_id[r] = 9000U + 7U * r;
    state.gas_cells.density_code[r] = rho;
    state.cells.mass_code[r] = rho;
    state.gas_cells.internal_energy_code[r] = 50.0 + static_cast<double>(r);
    state.gas_cells.metal_mass_code[r] = 0.0134 * rho;
    state.gas_cells.temperature_code[r] =
        h::temperatureKelvinPerInternalEnergyCode(5.0 / 3.0, 1.0e3) * state.gas_cells.internal_energy_code[r];
  }
  p::enableGasComposition(state, layout, {});
  return state;
}

[[nodiscard]] p::ChemistryShadowReport runRowList(const p::ThermochemistryOperator& chemistry,
                                                  const core::SimulationState& state,
                                                  const std::vector<std::uint32_t>& order);

[[nodiscard]] p::ChemistryShadowReport runRows(const p::ThermochemistryOperator& chemistry,
                                               const core::SimulationState& state, bool reverse) {
  std::vector<std::uint32_t> order(k_rows);
  std::iota(order.begin(), order.end(), 0U);
  if (reverse) std::reverse(order.begin(), order.end());
  return runRowList(chemistry, state, order);
}

// The runtime's controller rule (source_runtime updateGlobalChemistryControl):
// act only when the interval-global key-selected count fits the budget and no
// reference failed.
[[nodiscard]] bool controllerActs(std::uint64_t global_selected, std::uint64_t global_failures,
                                  std::uint32_t budget) {
  return global_failures == 0U && global_selected <= budget;
}

p::ChemistryShadowReport runRowList(const p::ThermochemistryOperator& chemistry,
                                    const core::SimulationState& state,
                                    const std::vector<std::uint32_t>& order) {
  const p::SourceInterval interval{.physical_seconds = 3.15576e12, .evaluation_scale_factor = 1.0,
                                   .accepted_step_index = 11U};
  p::ChemistryShadowReport shadow;
  for (const std::uint32_t row : order) {
    (void)chemistry.advanceRow(state, interval, row, p::k_chemistry_tag_detailed, shadow);
  }
  return shadow;
}

[[nodiscard]] p::ThermochemistryOperatorConfig config(std::uint32_t max_cells) {
  p::ThermochemistryOperatorConfig c;
  c.solver.backend = p::ThermochemistryBackend::kFixedDetailed;
  c.solver.allow_unqualified_rates = true;
  c.units = units();
  c.shadow.sample_fraction = 0.5;
  c.shadow.max_cells = max_cells;
  return c;
}

}  // namespace

int main() {
  const p::CompositionLayout layout = p::CompositionLayout::make(
      {.network = p::CompositionNetworkKind::kDetailedCarbonOxygen});
  try {
    const p::ThermochemistryOperator probe(config(1U << 20U), layout);
    (void)probe;
  } catch (const std::exception& error) {
    std::cout << "SKIP: detailed reference backend unavailable: " << error.what() << "\n";
    return 77;
  }
  const core::SimulationState state = gas(layout);

  // Uncapped: identical sample set and maxima in either row order.
  const p::ThermochemistryOperator uncapped(config(1U << 20U), layout);
  const p::ChemistryShadowReport forward = runRows(uncapped, state, false);
  const p::ChemistryShadowReport backward = runRows(uncapped, state, true);
  assert(forward.cap_truncated == 0U && backward.cap_truncated == 0U);
  assert(forward.samples == backward.samples);
  assert(forward.reference_failures == backward.reference_failures);
  assert(forward.max_relative_species_difference == backward.max_relative_species_difference);
  assert(forward.max_relative_energy_difference == backward.max_relative_energy_difference);
  const std::uint64_t selected = forward.samples + forward.reference_failures;
  assert(selected > 0U && selected < k_rows);  // key-based ~half of the rows

  assert(forward.selected == selected && backward.selected == selected);

  // Capped: the cut is counted, never silent; the selected total is unchanged.
  if (selected > 3U) {
    const p::ThermochemistryOperator capped(config(3U), layout);
    const p::ChemistryShadowReport cut = runRows(capped, state, false);
    assert(cut.samples + cut.reference_failures == 3U);
    assert(cut.cap_truncated == selected - 3U);
    assert(cut.selected == selected);
  }

  // Decomposition independence of the controller decision (rev-8 audit
  // example): budget below the one-rank selected count but above each half.
  // One "rank" truncates, two "ranks" do not; the decision must agree.
  if (selected >= 4U) {
    const auto budget = static_cast<std::uint32_t>(selected - 1U);
    const p::ThermochemistryOperator op(config(budget), layout);
    std::vector<std::uint32_t> all(k_rows);
    std::iota(all.begin(), all.end(), 0U);
    const std::vector<std::uint32_t> lower(all.begin(), all.begin() + k_rows / 2U);
    const std::vector<std::uint32_t> upper(all.begin() + k_rows / 2U, all.end());
    const p::ChemistryShadowReport one = runRowList(op, state, all);
    const p::ChemistryShadowReport a = runRowList(op, state, lower);
    const p::ChemistryShadowReport b = runRowList(op, state, upper);
    assert(one.selected == a.selected + b.selected);
    assert(controllerActs(one.selected, one.reference_failures, budget) ==
           controllerActs(a.selected + b.selected, a.reference_failures + b.reference_failures, budget));
    assert(!controllerActs(one.selected, 0U, budget));  // over budget: held on every layout
  }
  return 0;
}
