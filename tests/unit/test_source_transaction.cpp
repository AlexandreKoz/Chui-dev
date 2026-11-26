// CHUI-RT-BH-AGN-REF-v1 B02 transaction tests: exact-once commit, closure,
// floors, replay and stale-epoch refusal with NO mutation. NOT RUN.
#include "../support/test_require.hpp"

#include <cmath>
#include <stdexcept>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/source_transaction.hpp"

namespace p = cosmosim::physics;

namespace {

cosmosim::core::SimulationState makeGas() {
  cosmosim::core::SimulationState s;
  s.resizeCells(2);
  for (std::size_t i = 0; i < 2; ++i) {
    s.cells.mass_code[i] = 2.0;
    s.gas_cells.density_code[i] = 1.0;
    s.gas_cells.internal_energy_code[i] = 3.0;
    s.gas_cells.velocity_x_peculiar[i] = 1.0;
    s.gas_cells.metal_mass_code[i] = 0.02;
    s.gas_cells.gas_cell_id[i] = 100U + i;
  }
  return s;
}

p::SourceTransactionPlan withdrawal(const cosmosim::core::SimulationState& s, double f) {
  p::SourceTransactionPlan plan;
  plan.id = p::SourceTransactionId{7U, 1U, p::SourceChannel::kBhAccretionWithdrawal, 0U};
  plan.topology_epoch = s.cellIndexGeneration();
  p::GasCellIncrement inc;
  inc.gas_cell_id = 100U;
  inc.local_row = 0U;
  inc.withdrawal_fraction = f;
  inc.delta.mass_code = -f * 2.0;
  inc.delta.momentum_code = {-f * 2.0, 0.0, 0.0};
  inc.delta.total_energy_code = -f * 2.0 * (3.0 + 0.5);
  plan.cells.push_back(inc);
  plan.source_increment = inc.delta.negated();
  return plan;
}

void testWithdrawalLeavesSpecificStateUnchanged() {
  auto s = makeGas();
  p::SourceTransactionBatch batch;
  batch.add(withdrawal(s, 0.25));
  CHUI_REQUIRE(batch.prepare(s) && batch.validate() && batch.authorize(true));
  const auto rows = batch.commit(s);
  CHUI_REQUIRE(rows.size() == 1U);
  CHUI_REQUIRE(std::abs(s.cells.mass_code[0] - 1.5) < 1.0e-15);
  CHUI_REQUIRE(std::abs(s.gas_cells.internal_energy_code[0] - 3.0) < 1.0e-14);
  CHUI_REQUIRE(std::abs(s.gas_cells.velocity_x_peculiar[0] - 1.0) < 1.0e-15);
  CHUI_REQUIRE(std::abs(s.gas_cells.metal_mass_code[0] - 0.015) < 1.0e-15);
  CHUI_REQUIRE(std::abs(s.gas_cells.density_code[0] - 0.75) < 1.0e-15);
  CHUI_REQUIRE(batch.report().max_relative_closure_residual <= 1.0e-12);
  bool threw = false;
  try {
    (void)batch.commit(s);  // exactly once
  } catch (const std::logic_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testFloorAndReplayRejectionDoNotMutate() {
  auto s = makeGas();
  p::SourceTransactionBatch over;
  over.add(withdrawal(s, 0.6));
  auto second = withdrawal(s, 0.6);
  second.id.source_id = 8U;
  over.add(second);  // overlapping withdrawals exceed the cell mass
  CHUI_REQUIRE(!over.prepare(s));
  CHUI_REQUIRE(s.cells.mass_code[0] == 2.0);
  p::SourceTransactionBatch replay;
  replay.add(withdrawal(s, 0.1));
  replay.add(withdrawal(s, 0.1));  // identical identity
  CHUI_REQUIRE(!replay.prepare(s));
  p::SourceTransactionBatch stale;
  auto plan = withdrawal(s, 0.1);
  plan.topology_epoch += 1U;
  stale.add(plan);
  CHUI_REQUIRE(!stale.prepare(s));
  p::SourceTransactionBatch vote;
  vote.add(withdrawal(s, 0.1));
  CHUI_REQUIRE(vote.prepare(s) && vote.validate());
  CHUI_REQUIRE(!vote.authorize(false));  // one participant refused: nothing commits
  CHUI_REQUIRE(s.cells.mass_code[0] == 2.0);
}

void testInjectionThermalIncrementFromEnergyMinusKinetic() {
  auto s = makeGas();
  p::SourceTransactionBatch batch;
  p::SourceTransactionPlan plan;
  plan.id = p::SourceTransactionId{9U, 1U, p::SourceChannel::kAgnWind, 0U};
  plan.topology_epoch = s.cellIndexGeneration();
  p::GasCellIncrement inc;
  inc.gas_cell_id = 101U;
  inc.local_row = 1U;
  inc.delta.mass_code = 1.0;
  inc.delta.momentum_code = {0.0, 0.0, 0.0};  // injected at rest
  inc.delta.total_energy_code = 0.0;
  plan.cells.push_back(inc);
  plan.source_increment = inc.delta.negated();
  batch.add(plan);
  CHUI_REQUIRE(batch.prepare(s) && batch.validate() && batch.authorize(true));
  (void)batch.commit(s);
  // Inelastic mixing: momentum conserved, the kinetic loss heats the gas.
  const double ke_before = 0.5 * 2.0 * 1.0;
  const double ke_after = 0.5 * (2.0 * 2.0) / 3.0;
  const double u_expected = (2.0 * 3.0 + ke_before - ke_after) / 3.0;
  CHUI_REQUIRE(std::abs(s.gas_cells.velocity_x_peculiar[1] - 2.0 / 3.0) < 1.0e-15);
  CHUI_REQUIRE(std::abs(s.gas_cells.internal_energy_code[1] - u_expected) < 1.0e-14);
}

}  // namespace

int main() {
  testWithdrawalLeavesSpecificStateUnchanged();
  testFloorAndReplayRejectionDoNotMutate();
  testInjectionThermalIncrementFromEnergyMinusKinetic();
  return 0;
}
