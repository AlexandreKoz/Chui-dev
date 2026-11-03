// CHUI Part 2 — composition substrate (NOT RUN — USER-RESERVED).
//
// Invariants exercised:
//  * element catalog: A09 solar pattern gives X~0.738, Y~0.249, Z~0.0134
//    (Asplund et al. 2009 Table 4 present-day photosphere values).
//  * layout manifest: stable signature, disjoint lane ranges, zero lanes when
//    disabled (DMO / plain hydro carry no composition storage).
//  * initialization: element masses sum exactly to gas mass, metals follow
//    solar pattern scaled to the cell metal mass.
//  * prolong/restrict: per-lane conservation to round-off.
//  * migration payload: pack -> wire -> unpack reproduces lanes and history;
//    width mismatch is refused (no silently lost lanes).
//  * row identity guard: a permuted composition row is detected (audit risk
//    "cooling/composition row permutation").
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/gas_composition_sidecar.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "../../src/workflows/internal/migration_wire.hpp"

namespace {

using namespace cosmosim;
using physics::ChemicalElement;

bool near(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max({std::abs(a), std::abs(b), 1.0e-300});
}

void testSolarPattern() {
  const auto& f = physics::solarElementMassFractions();
  double sum = 0.0;
  for (const double v : f) sum += v;
  assert(near(sum, 1.0, 1.0e-14));
  assert(std::abs(f[0] - 0.7381) < 2.0e-3);
  assert(std::abs(f[1] - 0.2485) < 2.0e-3);
  assert(std::abs(physics::solarMetalMassFraction() - 0.0134) < 1.0e-3);
  // Oxygen is the most abundant metal by mass in A09.
  for (std::size_t e = 2; e < physics::k_chemical_element_count; ++e) {
    assert(f[static_cast<std::size_t>(ChemicalElement::kO)] >= f[e]);
  }
}

void testLayoutManifest() {
  const auto elements_only = physics::CompositionLayout::make({});
  assert(elements_only.descriptor().lane_count == 11U);
  assert(elements_only.descriptor().species_lane_count == 0U);
  const auto hybrid = physics::CompositionLayout::make(
      {.network = physics::CompositionNetworkKind::kHybridLight,
       .dust = physics::CompositionDustKind::kTwoReservoirV1,
       .sgs_energy = true});
  assert(hybrid.descriptor().lane_count == 11U + 8U + 2U + 1U);
  assert(hybrid.hasSgs());
  assert(hybrid.descriptor().layout_signature != elements_only.descriptor().layout_signature);
  // Signature is deterministic.
  const auto again = physics::CompositionLayout::make(hybrid.request());
  assert(again.descriptor() == hybrid.descriptor());
  bool threw = false;
  try { (void)elements_only.speciesLane(0); } catch (const std::out_of_range&) { threw = true; }
  assert(threw);
}

void testInitialElementMasses() {
  const physics::PrimordialCompositionPolicy policy{};
  const auto primordial = physics::initialElementMasses(2.0, 0.0, policy);
  assert(near(primordial[1], 2.0 * 0.2449, 1.0e-14));
  assert(near(primordial[0] + primordial[1], 2.0, 1.0e-14));
  const double z = physics::solarMetalMassFraction();
  const auto solar = physics::initialElementMasses(1.0, z, policy);
  double sum = 0.0, metals = 0.0;
  for (std::size_t e = 0; e < solar.size(); ++e) {
    sum += solar[e];
    if (e >= 2) metals += solar[e];
  }
  assert(near(sum, 1.0, 1.0e-13));
  assert(near(metals, z, 1.0e-13));
  // At solar metallicity helium reaches the solar value.
  assert(std::abs(solar[1] - physics::solarElementMassFractions()[1]) < 1.0e-12);
}

void testProlongRestrictConservation() {
  const std::vector<double> parent{3.0, 1.5, 0.25, 1.0e-9};
  const std::vector<double> share{0.1, 0.2, 0.3, 0.15, 0.05, 0.05, 0.1, 0.05};
  std::vector<double> children(parent.size() * share.size());
  physics::prolongComposition(parent, share, children);
  std::vector<double> back(parent.size());
  physics::restrictComposition(children, parent.size(), back);
  for (std::size_t l = 0; l < parent.size(); ++l) assert(near(back[l], parent[l], 1.0e-14));
}

core::SimulationState makeState(std::size_t cells) {
  core::SimulationState state;
  state.resizeCells(cells);
  for (std::size_t i = 0; i < cells; ++i) {
    state.gas_cells.gas_cell_id[i] = 100U + i;
    state.cells.mass_code[i] = 1.0 + static_cast<double>(i);
    state.gas_cells.metal_mass_code[i] = 0.01 * state.cells.mass_code[i];
  }
  return state;
}

void testEnableAndMirror() {
  auto state = makeState(4);
  const auto layout = physics::CompositionLayout::make(
      {.network = physics::CompositionNetworkKind::kHybridLight});
  physics::enableGasComposition(state, layout, {});
  assert(state.gas_composition.enabled());
  for (std::size_t row = 0; row < 4; ++row) {
    assert(near(physics::compositionGasMass(state.gas_composition, row),
                state.cells.mass_code[row], 1.0e-13));
    // Species lanes of H sum to the H element lane.
    const double h = state.gas_composition.value(layout.elementLane(ChemicalElement::kH), row);
    const double hi = state.gas_composition.value(
        layout.speciesLane(static_cast<std::size_t>(physics::ChemicalSpecies::kHI)), row);
    assert(near(h, hi, 1.0e-15));
  }
  const double drift = physics::synchronizeMetalMirrorFromComposition(state, layout, {});
  assert(drift < 1.0e-12);
}

void testPackUnpackAndWire() {
  auto state = makeState(2);
  const auto layout = physics::CompositionLayout::make(
      {.network = physics::CompositionNetworkKind::kDetailedCarbonOxygen,
       .dust = physics::CompositionDustKind::kTwoReservoirV1});
  physics::enableGasComposition(state, layout, {});
  state.gas_composition.chemistry_model_tag[1] = 2U;
  state.gas_composition.model_dwell_years[1] = 3.5e6;
  const std::vector<double> payload = state.gas_composition.packRow(1);
  assert(payload.size() == state.gas_composition.packedRowWidth());

  core::ParticleMigrationRecord record;
  record.particle_id = 7U;
  record.species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
  record.has_gas_cell_fields = true;
  record.gas_cell_fields.gas_cell_id = 101U;
  record.gas_cell_fields.has_parent_particle = 1U;
  record.gas_cell_fields.parent_particle_id = 7U;
  record.gas_cell_fields.composition_payload = payload;
  const auto encoded = workflows::internal::migration_wire::encodeParticleMigrationRecord(record);
  const auto decoded = workflows::internal::migration_wire::decodeParticleMigrationRecord(encoded);
  assert(decoded.gas_cell_fields.composition_payload == payload);

  auto receiver = state.gas_composition.emptyLike(1);
  receiver.unpackRow(0, decoded.gas_cell_fields.composition_payload, 101U);
  assert(receiver.chemistry_model_tag[0] == 2U);
  assert(receiver.model_dwell_years[0] == 3.5e6);
  for (std::size_t l = 0; l < receiver.laneCount(); ++l) {
    assert(receiver.value(l, 0) == state.gas_composition.value(l, 1));
  }
  std::vector<double> truncated = payload;
  truncated.pop_back();
  bool threw = false;
  try { receiver.unpackRow(0, truncated, 101U); } catch (const std::runtime_error&) { threw = true; }
  assert(threw);
}

void testRowPermutationGuard() {
  auto state = makeState(3);
  physics::enableGasComposition(state, physics::CompositionLayout::make({}), {});
  state.gas_composition.requireAlignedWith(state.gas_cells.gas_cell_id, "test");
  std::swap(state.gas_cells.gas_cell_id[0], state.gas_cells.gas_cell_id[2]);
  bool threw = false;
  try {
    state.gas_composition.requireAlignedWith(state.gas_cells.gas_cell_id, "test");
  } catch (const std::runtime_error&) { threw = true; }
  assert(threw);
}

void testDisabledHasNoStorage() {
  core::SimulationState state;
  state.resizeCells(1000);
  assert(!state.gas_composition.enabled());
  assert(state.gas_composition.ownedCapacityBytes() == 0U);
  assert(state.gas_composition.packRow(0).empty());
}

void testRemoveFractionTransfersElements() {
  auto state = makeState(1);
  const auto layout = physics::CompositionLayout::make(
      {.dust = physics::CompositionDustKind::kTwoReservoirV1});
  physics::enableGasComposition(state, layout, {});
  state.gas_composition.setValue(layout.dustLane(0), 0, 0.002);
  const double before = physics::compositionGasMass(state.gas_composition, 0);
  const auto removed = physics::removeGasFractionFromComposition(
      state.gas_composition, layout, 0, 0.25);
  double removed_sum = 0.0;
  for (const double m : removed) removed_sum += m;
  const double after = physics::compositionGasMass(state.gas_composition, 0);
  assert(near(before, after + removed_sum, 1.0e-14));
  assert(near(removed_sum, 0.25 * before, 1.0e-13));
}

}  // namespace

int main() {
  testSolarPattern();
  testLayoutManifest();
  testInitialElementMasses();
  testProlongRestrictConservation();
  testEnableAndMirror();
  testPackUnpackAndWire();
  testRowPermutationGuard();
  testDisabledHasNoStorage();
  testRemoveFractionTransfersElements();
  return 0;
}
