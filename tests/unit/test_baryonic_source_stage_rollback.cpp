// CHUI Part 2 corrective gate — local source-stage atomicity (NOT RUN —
// USER-RESERVED).
//
// Failure injection: one row cannot integrate the interval (substep budget of
// one on a stiff, hot, dense cell) while other rows would change their ISM
// regime, chemistry-fidelity tag, dwell counters, dust and SGS lanes. The
// stage is driven through a StageTransaction exactly as the SourceRuntime
// does. After the rejection every canonical gas lane, every composition lane,
// all tags/dwell/mu metadata and the published (committed) diagnostics must
// be bit-identical to the step-start values. A successful stage commits all
// staged changes together and publishes the diagnostics.
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/stage_transaction.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/baryonic_source_stage.hpp"

namespace {

using namespace cosmosim;

struct Snapshot {
  std::vector<double> gas;
  std::vector<double> lanes;
  std::vector<std::uint8_t> tags;
  std::vector<double> meta;
};

void appendBits(std::vector<double>& out, const double* begin, std::size_t n) { out.insert(out.end(), begin, begin + n); }

Snapshot snapshot(const core::SimulationState& s) {
  Snapshot snap;
  const std::size_t n = s.cells.size();
  for (const auto* lane : {&s.gas_cells.density_code, &s.gas_cells.internal_energy_code, &s.gas_cells.pressure_code,
                           &s.gas_cells.temperature_code, &s.gas_cells.sound_speed_code, &s.gas_cells.metal_mass_code,
                           &s.gas_cells.velocity_x_peculiar}) {
    appendBits(snap.gas, lane->data(), n);
  }
  appendBits(snap.gas, s.cells.mass_code.data(), n);
  for (std::size_t l = 0; l < s.gas_composition.laneCount(); ++l) {
    const auto v = s.gas_composition.lane(l);
    appendBits(snap.lanes, v.data(), v.size());
  }
  snap.tags.insert(snap.tags.end(), s.gas_composition.ism_regime_tag.begin(), s.gas_composition.ism_regime_tag.end());
  snap.tags.insert(snap.tags.end(), s.gas_composition.chemistry_model_tag.begin(),
                   s.gas_composition.chemistry_model_tag.end());
  appendBits(snap.meta, s.gas_composition.model_dwell_years.data(), n);
  appendBits(snap.meta, s.gas_composition.ism_regime_dwell_years.data(), n);
  appendBits(snap.meta, s.gas_composition.mean_molecular_weight.data(), n);
  return snap;
}

bool bitIdentical(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

physics::ThermochemistryUnitBridge units() {
  physics::ThermochemistryUnitBridge u;
  u.density_code_to_cgs = 1.6726e-24;  // code density 1 ~ n_H 0.76
  u.velocity_code_to_cm_s = 1.0e5;
  u.length_code_to_cm = 3.0857e21;     // kpc
  u.mass_code_to_g = 1.989e33;
  return u;
}

core::SimulationState makeState(const physics::CompositionLayout& layout) {
  core::SimulationState state;
  state.resizeCells(3);
  for (std::size_t r = 0; r < 3; ++r) {
    state.gas_cells.gas_cell_id[r] = 77U + r;
    state.cells.mass_code[r] = 1.0;
    state.gas_cells.density_code[r] = 50.0;              // dense: regime change candidate
    state.gas_cells.internal_energy_code[r] = 0.01;      // cold
    state.gas_cells.temperature_code[r] = 50.0;
    state.gas_cells.metal_mass_code[r] = 0.0134;
  }
  state.gas_cells.internal_energy_code[2] = 1.0e4;        // hot dense: stiff, fails with one substep
  physics::enableGasComposition(state, layout, {});
  for (std::size_t r = 0; r < 3; ++r) {
    state.gas_composition.setValue(layout.dustLane(0), r, 1.0e-6);
    state.gas_composition.setValue(state.gas_composition.descriptor().sgs_energy_lane, r, 0.5);
    state.gas_composition.model_dwell_years[r] = 5.0e6;
    state.gas_composition.ism_regime_dwell_years[r] = 5.0e6;
  }
  return state;
}

physics::BaryonicSourceStageConfig stageConfig() {
  physics::BaryonicSourceStageConfig c;
  c.regime = physics::IsmRegimeConfig{.effective_threshold_n_h_cm3 = 0.13, .minimum_dwell_years = 0.0};
  c.dust = physics::DustConfig{};
  c.sgs = physics::SgsConfig{};
  return c;
}

void testRejectionPreservesEverything() {
  const auto layout = physics::CompositionLayout::make({.network = physics::CompositionNetworkKind::kHybridLight,
                                                        .dust = physics::CompositionDustKind::kTwoReservoirV1,
                                                        .sgs_energy = true});
  auto state = makeState(layout);
  physics::ThermochemistryOperatorConfig op_config;
  op_config.solver.backend = physics::ThermochemistryBackend::kFixedHybrid;
  op_config.solver.max_substeps = 1U;  // failure injection
  op_config.units = units();
  op_config.solver.allow_unqualified_rates = true;
  const physics::ThermochemistryOperator chemistry(op_config, layout);
  physics::BaryonicSourceStage stage(chemistry, stageConfig());
  const Snapshot before = snapshot(state);
  const physics::SourceInterval interval{.physical_seconds = 3.15576e15, .evaluation_scale_factor = 1.0,
                                         .accepted_step_index = 9U};
  core::StageTransaction transaction("test.baryonic_stage", nullptr);
  const bool prepared = transaction.prepare([&] { stage.prepare(state, interval, {}); });
  assert(!prepared);
  assert(!transaction.validate(true, core::StageRejectionKind::kNone, {}));
  assert(transaction.outcome().rejection == core::StageRejectionKind::kNumericalRetry);
  stage.discard();
  const Snapshot after = snapshot(state);
  assert(bitIdentical(before.gas, after.gas));
  assert(bitIdentical(before.lanes, after.lanes));
  assert(before.tags == after.tags);
  // NaN mu (unknown) compares bitwise identical as well.
  assert(bitIdentical(before.meta, after.meta));
  assert(stage.lastCommittedReport().chemistry.cells == 0U);
  assert(stage.lastCommittedReport().regime.changed_rows.empty());
  assert(!stage.prepared());
}

void testSuccessCommitsAllTogether() {
  const auto layout = physics::CompositionLayout::make({.network = physics::CompositionNetworkKind::kHybridLight,
                                                        .dust = physics::CompositionDustKind::kTwoReservoirV1,
                                                        .sgs_energy = true});
  auto state = makeState(layout);
  state.gas_cells.internal_energy_code[2] = 0.01;  // no stiff row
  physics::ThermochemistryOperatorConfig op_config;
  op_config.solver.backend = physics::ThermochemistryBackend::kFixedHybrid;
  op_config.units = units();
  op_config.solver.allow_unqualified_rates = true;
  const physics::ThermochemistryOperator chemistry(op_config, layout);
  physics::BaryonicSourceStage stage(chemistry, stageConfig());
  double element_before = 0.0, energy_before = 0.0;
  for (std::size_t r = 0; r < 3; ++r) {
    element_before += physics::compositionGasMass(state.gas_composition, r);
    energy_before += state.gas_composition.value(state.gas_composition.descriptor().sgs_energy_lane, r);
  }
  const physics::SourceInterval interval{.physical_seconds = 3.15576e12, .evaluation_scale_factor = 1.0};
  core::StageTransaction transaction("test.baryonic_stage", nullptr);
  assert(transaction.prepare([&] { stage.prepare(state, interval, {}); }));
  assert(transaction.validate(true, core::StageRejectionKind::kNone, {}));
  std::vector<std::uint32_t> refresh;
  (void)transaction.commit([&] { refresh = stage.commit(state); });
  assert(!refresh.empty());
  double element_after = 0.0;
  for (std::size_t r = 0; r < 3; ++r) element_after += physics::compositionGasMass(state.gas_composition, r);
  assert(std::abs(element_after - element_before) < 1.0e-14 * element_before);  // gas + dust conserved
  // Regime decisions were applied together with the chemistry (dense cold
  // 1-kpc cells are Jeans-unresolved -> effective) and dwell restarted.
  for (std::size_t r = 0; r < 3; ++r) {
    assert(state.gas_composition.ism_regime_tag[r] == core::k_gas_ism_regime_effective);
    assert(state.gas_composition.ism_regime_dwell_years[r] == 0.0);
  }
  assert(stage.lastCommittedReport().regime.to_effective == 3U);
  assert(stage.lastCommittedReport().sgs.dissipated_to_heat_code > 0.0);
  (void)energy_before;
}

}  // namespace

int main() {
  testRejectionPreservesEverything();
  testSuccessCommitsAllTogether();
  return 0;
}
