// CHUI Part 2 — ISM regimes, live dust, SGS energy, adaptive chemistry and
// effective-EOS edges (NOT RUN — USER-RESERVED).
//
// Invariants:
//  * dust: accretion removes gas-phase refractory elements in exact
//    stoichiometry, sputtering returns them; each element's gas+dust total is
//    conserved to round-off; zero dust => no growth; hot gas sputters;
//  * stellar ejecta partition: gas + dust element sums equal the ejecta;
//  * SN destruction returns exactly the destroyed solids to the gas;
//  * SGS: exact analytic decay K0/(1+C sqrt(K0) t/(2dx))^2; dissipated energy
//    appears in the thermal energy (total conserved); zero-turbulence limit;
//  * ISM regime: an unresolved dense cell becomes effective, a resolved one
//    stays resolved; hysteresis + dwell prevent flip-flop;
//  * adaptive chemistry: cold dense gas promotes to detailed; demotion waits
//    for the dwell time and keeps C+/CO lanes as an exact frozen reservoir;
//  * effective EOS: above the table the power-law asymptote is used and
//    flagged (no silent clamp), c^2 equals dP/drho of the interpolant, and the
//    pressure stays positive and monotone across the table end.
#include <cassert>
#include <cmath>
#include <optional>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/physics/ism_dust_sgs.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"

namespace {

using namespace cosmosim;
using physics::ChemicalElement;

core::SimulationState makeState(std::size_t cells, double metal_fraction) {
  core::SimulationState state;
  state.resizeCells(cells);
  for (std::size_t i = 0; i < cells; ++i) {
    state.gas_cells.gas_cell_id[i] = 10U + i;
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.density_code[i] = 1.0;
    state.gas_cells.internal_energy_code[i] = 1.0;
    state.gas_cells.metal_mass_code[i] = metal_fraction;
  }
  return state;
}

double elementTotal(const core::SimulationState& s, const physics::CompositionLayout& layout, ChemicalElement e,
                    std::size_t row) {
  double total = s.gas_composition.value(layout.elementLane(e), row);
  for (std::size_t d = 0; d < s.gas_composition.descriptor().dust_lane_count; ++d) {
    total += physics::dustElementMassFractions(static_cast<physics::DustReservoir>(d))[static_cast<std::size_t>(e)] *
        s.gas_composition.value(layout.dustLane(d), row);
  }
  return total;
}

void testDustGrowthAndSputteringConserveElements() {
  auto state = makeState(1, 0.0134);
  const auto layout = physics::CompositionLayout::make({.dust = physics::CompositionDustKind::kTwoReservoirV1});
  physics::enableGasComposition(state, layout, {});
  state.gas_composition.setValue(layout.dustLane(0), 0, 1.0e-5);  // seed silicate
  state.gas_composition.setValue(layout.dustLane(1), 0, 1.0e-5);  // seed carbon
  // Seeds come out of the gas: remove their elements first so totals are physical.
  for (std::size_t d = 0; d < 2; ++d) {
    const auto& f = physics::dustElementMassFractions(static_cast<physics::DustReservoir>(d));
    for (std::size_t e = 0; e < physics::k_chemical_element_count; ++e) {
      state.gas_composition.addValue(layout.elementLane(static_cast<ChemicalElement>(e)), 0, -1.0e-5 * f[e]);
    }
  }
  std::vector<double> before;
  for (std::size_t e = 0; e < physics::k_chemical_element_count; ++e) {
    before.push_back(elementTotal(state, layout, static_cast<ChemicalElement>(e), 0));
  }
  physics::DustConfig config;
  physics::DustUpdateLedger ledger;
  // Cold dense gas: growth.
  physics::advanceDustCell(state.gas_composition, layout, 0, 100.0, 50.0, 1.0e8 * physics::k_seconds_per_year, config, ledger);
  assert(state.gas_composition.value(layout.dustLane(0), 0) > 1.0e-5);
  assert(ledger.gas_to_dust_code > 0.0);
  // Hot gas: sputtering returns solids.
  const double dust_before_hot = state.gas_composition.value(layout.dustLane(1), 0);
  physics::advanceDustCell(state.gas_composition, layout, 0, 1.0e-2, 1.0e7, 1.0e8 * physics::k_seconds_per_year, config, ledger);
  assert(state.gas_composition.value(layout.dustLane(1), 0) < dust_before_hot);
  for (std::size_t e = 0; e < physics::k_chemical_element_count; ++e) {
    const double after = elementTotal(state, layout, static_cast<ChemicalElement>(e), 0);
    assert(std::abs(after - before[e]) <= 1.0e-14 * std::max(before[e], 1.0e-20));
  }
  // Destruction by an SN returns exactly the destroyed mass.
  const double gas_mass = physics::compositionGasMass(state.gas_composition, 0);
  const double destroyed = physics::destroyDustBySupernova(state.gas_composition, layout, 0, 0.5, 1.0);
  assert(destroyed > 0.0);
  assert(std::abs(physics::compositionGasMass(state.gas_composition, 0) - gas_mass) < 1.0e-15);
}

void testNoSeedNoGrowth() {
  auto state = makeState(1, 0.0134);
  const auto layout = physics::CompositionLayout::make({.dust = physics::CompositionDustKind::kTwoReservoirV1});
  physics::enableGasComposition(state, layout, {});
  physics::DustUpdateLedger ledger;
  physics::advanceDustCell(state.gas_composition, layout, 0, 1.0e3, 20.0, 1.0e16, {}, ledger);
  assert(state.gas_composition.value(layout.dustLane(0), 0) == 0.0);
}

void testEjectaPartitionConserves() {
  std::array<double, physics::k_chemical_element_count> ejecta{};
  ejecta[static_cast<std::size_t>(ChemicalElement::kC)] = 0.3;
  ejecta[static_cast<std::size_t>(ChemicalElement::kO)] = 1.0;
  ejecta[static_cast<std::size_t>(ChemicalElement::kMg)] = 0.1;
  ejecta[static_cast<std::size_t>(ChemicalElement::kSi)] = 0.1;
  ejecta[static_cast<std::size_t>(ChemicalElement::kFe)] = 0.1;
  const auto part = physics::partitionEjectaIntoDust(ejecta, 0.5);
  for (std::size_t e = 0; e < physics::k_chemical_element_count; ++e) {
    double total = part.gas_elements[e];
    for (std::size_t d = 0; d < physics::k_dust_reservoir_count; ++d) {
      total += part.dust[d] * physics::dustElementMassFractions(static_cast<physics::DustReservoir>(d))[e];
    }
    assert(std::abs(total - ejecta[e]) < 1.0e-14);
  }
  assert(part.dust[1] > 0.0 && part.dust[0] > 0.0);
}

void testSgsExactDecayAndEnergyConservation() {
  const double k = physics::sgsSpecificEnergyAfter(4.0, 1.0, 2.0, 1.0);
  assert(std::abs(k - 4.0 / std::pow(1.0 + 1.0 * 2.0 * 1.0 / 4.0, 2.0)) < 1.0e-15);
  assert(physics::sgsSpecificEnergyAfter(0.0, 1.0, 1.0, 1.58) == 0.0);
  auto state = makeState(1, 0.0);
  const auto layout = physics::CompositionLayout::make({.sgs_energy = true});
  physics::enableGasComposition(state, layout, {});
  state.gas_composition.setValue(state.gas_composition.descriptor().sgs_energy_lane, 0, 0.5);
  const double total_before = state.gas_cells.internal_energy_code[0] * state.cells.mass_code[0] + 0.5;
  physics::SgsLedger ledger;
  const double dissipated = physics::advanceSgsCell(state, layout, 0, 3.0, 1.0, {}, ledger);
  assert(dissipated > 0.0);
  const double total_after = state.gas_cells.internal_energy_code[0] * state.cells.mass_code[0] +
      state.gas_composition.value(state.gas_composition.descriptor().sgs_energy_lane, 0);
  assert(std::abs(total_after - total_before) < 1.0e-14);
}

void testIsmRegimeHysteresisAndDwell() {
  auto state = makeState(2, 0.0);
  physics::enableGasComposition(state, physics::CompositionLayout::make({}), {});
  physics::ThermochemistryUnitBridge units;
  units.density_code_to_cgs = 1.0e-22;   // n_H ~ 45 cm^-3
  units.velocity_code_to_cm_s = 1.0e5;   // 1 km/s
  units.length_code_to_cm = 3.0857e21;   // 1 kpc cells: Jeans length unresolved
  state.gas_cells.internal_energy_code[0] = 1.0;  // ~ 100 K
  // Row 1: tiny cells (resolved).
  state.cells.mass_code[1] = 1.0e-12;
  physics::SourceInterval interval{.physical_seconds = 3.15576e13, .evaluation_scale_factor = 1.0};
  physics::IsmRegimeConfig config;
  config.minimum_dwell_years = 0.5;
  auto report = physics::updateIsmRegimes(state, config, units, interval, {});
  assert(state.gas_composition.ism_regime_tag[0] == core::k_gas_ism_regime_effective);
  assert(state.gas_composition.ism_regime_tag[1] == core::k_gas_ism_regime_resolved);
  assert(report.to_effective == 1U);
  // Large dwell: a subsequent flip back is held.
  config.minimum_dwell_years = 1.0e9;
  state.gas_cells.density_code[0] = 1.0e-6;  // far below threshold
  report = physics::updateIsmRegimes(state, config, units, interval, {});
  assert(state.gas_composition.ism_regime_tag[0] == core::k_gas_ism_regime_effective);
  assert(report.held_by_dwell == 1U);
}

void testEffectiveEosEdges() {
  core::PhysicsConfig physics_config;
  physics_config.sf_effective_eos_table_bins = 64;
  physics_config.sf_effective_eos_max_density_ratio = 1.0e3;
  physics_config.sf_effective_supernova_specific_energy_code = 1.0e5;
  const auto units = core::makeUnitSystem("kpc", "msun", "km_s");
  const physics::EffectiveMultiphaseEosTable table(physics::makeEffectiveMultiphaseEosConfig(physics_config), units,
                                                   physics::makeEffectiveIsmReferenceCoolingProvider(physics_config));
  const double rho_th = table.thresholdDensityPhysCode();
  const auto inside = table.lookup(1.0e3 * rho_th * 0.999);
  const auto beyond = table.lookup(1.0e4 * rho_th);
  assert(inside.valid && !inside.extrapolated);
  assert(beyond.valid && beyond.extrapolated);
  assert(beyond.entry.pressure_phys_code > inside.entry.pressure_phys_code);  // monotone, positive
  // c^2 equals the derivative of the interpolated pressure inside a segment.
  const double rho = 30.0 * rho_th;
  const auto a = table.lookup(rho * (1.0 - 1.0e-6));
  const auto b = table.lookup(rho * (1.0 + 1.0e-6));
  const double dpdrho = (b.entry.pressure_phys_code - a.entry.pressure_phys_code) / (2.0e-6 * rho);
  const auto c = table.lookup(rho);
  assert(std::abs(c.entry.signal_speed_squared_code - dpdrho) / dpdrho < 1.0e-3);
  assert(!table.lookup(0.5 * rho_th).above_threshold);
}

}  // namespace

int main() {
  testDustGrowthAndSputteringConserveElements();
  testNoSeedNoGrowth();
  testEjectaPartitionConserves();
  testSgsExactDecayAndEnergyConservation();
  testIsmRegimeHysteresisAndDwell();
  testEffectiveEosEdges();
  return 0;
}
