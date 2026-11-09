// CHUI patch 53 — C02 metal-ion dataset integration (NOT RUN — USER-RESERVED).
//
// Fixture cases (always): the loader grammar refuses an incomplete stage, a
// duplicate rate line, an ionizing bare nucleus, a recombining neutral atom,
// a duplicate/too-high ionization potential and an element absent from the
// file; the solver refuses an unverified table without the research opt-in
// and accepts it with the opt-in; adaptive fidelity with ladders and no
// element table is refused (demoted rows would lose the ladder elements'
// cooling).
//
// Real-dataset cases (when CHUI_C02_TABLE names the user-local
// chui_metal_ion_rates_v1_C_O_Si_Fe_unverified.txt, sha256 633a66ee...):
//   * every subset of {C, O, Si, Fe} loads; the full set has 58 stages with
//     the three temperature processes and the energy contract (54 edges);
//   * the declared verification is unverified (never upgraded);
//   * -inf entries evaluate to exact zero rates at the boundaries (bare
//     nucleus ionization, neutral recombination) and inside the grid;
//   * loadBaryonicDataProducts refuses a wrong SHA-256;
//   * T outside [1e4, 1e8] K throws under the reject policy and is counted
//     under clamp.
// The fixture numbers are arbitrary test values, not atomic data.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/physics/baryonic_runtime_config.hpp"
#include "cosmosim/physics/metal_ion_network.hpp"
#include "cosmosim/physics/thermochemistry.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"
#include "../support/test_require.hpp"

namespace {
namespace ph = cosmosim::physics;

template <class F>
bool throws(F&& f) {
  try {
    f();
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

std::string carbonFixture(bool with_contract) {
  std::string text =
      "schema chui_metal_ion_rates_v1\n"
      "provenance synthetic unit-test fixture (not physical data)\n"
      "license test-only\n"
      "verification unverified\n"
      "temperature_grid_log10_k 2 4 5\n";
  for (int k = 0; k < 7; ++k) {
    const std::string stage = std::to_string(k);
    text += "rate C " + stage + " collisional_ionization " + (k == 6 ? "-inf -inf" : "-10 -9") + "\n";
    text += "rate C " + stage + " recombination " + (k == 0 ? "-inf -inf" : "-12 -12") + "\n";
    text += "rate C " + stage + " cooling -23 -22\n";
  }
  if (with_contract) {
    text += "energy_convention collisional_ionization_work_in_cooling\n";
    for (int k = 0; k < 6; ++k) text += "ionization_potential C " + std::to_string(k) + " " + std::to_string(k + 1) + "\n";
  }
  return text;
}

void testGrammarRefusals() {
  const auto c = ph::MetalIonElementSet::parse("C");
  CHUI_REQUIRE(!throws([&] { (void)ph::MetalIonRateTable::parse(carbonFixture(true), "ok", c); }));
  std::string missing = carbonFixture(true);
  missing.erase(missing.find("rate C 3 cooling"), std::string("rate C 3 cooling -23 -22\n").size());
  CHUI_REQUIRE(throws([&] { (void)ph::MetalIonRateTable::parse(missing, "missing", c); }));
  CHUI_REQUIRE(throws([&] {
    (void)ph::MetalIonRateTable::parse(carbonFixture(true) + "rate C 2 cooling -23 -22\n", "dup", c);
  }));
  std::string bare = carbonFixture(true);
  bare.replace(bare.find("rate C 6 collisional_ionization -inf -inf"),
               std::string("rate C 6 collisional_ionization -inf -inf").size(), "rate C 6 collisional_ionization -10 -9");
  CHUI_REQUIRE(throws([&] { (void)ph::MetalIonRateTable::parse(bare, "bare", c); }));
  std::string neutral = carbonFixture(true);
  neutral.replace(neutral.find("rate C 0 recombination -inf -inf"),
                  std::string("rate C 0 recombination -inf -inf").size(), "rate C 0 recombination -12 -12");
  CHUI_REQUIRE(throws([&] { (void)ph::MetalIonRateTable::parse(neutral, "neutral", c); }));
  CHUI_REQUIRE(throws([&] {
    (void)ph::MetalIonRateTable::parse(carbonFixture(true) + "ionization_potential C 0 1\n", "dup_ip", c);
  }));
  CHUI_REQUIRE(throws([&] {
    (void)ph::MetalIonRateTable::parse(carbonFixture(false) + "ionization_potential C 6 1\n", "bare_ip", c);
  }));
  CHUI_REQUIRE(throws([&] {
    (void)ph::MetalIonRateTable::parse(carbonFixture(true), "absent", ph::MetalIonElementSet::parse("C,Fe"));
  }));
}

ph::ThermochemistryConfig detailedConfig(const ph::MetalIonRateTable* table, const char* elements, bool opt_in) {
  ph::ThermochemistryConfig config;
  config.backend = ph::ThermochemistryBackend::kFixedDetailed;
  config.metal_ion_table = table;
  config.metal_ion_elements = ph::MetalIonElementSet::parse(elements);
  config.allow_unqualified_rates = opt_in;
  return config;
}

void testOptInGate() {
  const auto table = ph::MetalIonRateTable::parse(carbonFixture(true), "fixture", ph::MetalIonElementSet::parse("C"));
  CHUI_REQUIRE(!table.declaredVerified());
  CHUI_REQUIRE(throws([&] { ph::ThermochemistrySolver solver(detailedConfig(&table, "C", false)); }));
  CHUI_REQUIRE(!throws([&] { ph::ThermochemistrySolver solver(detailedConfig(&table, "C", true)); }));
  // A table without the energy contract is refused even with the opt-in.
  const auto no_contract =
      ph::MetalIonRateTable::parse(carbonFixture(false), "fixture", ph::MetalIonElementSet::parse("C"));
  CHUI_REQUIRE(throws([&] { ph::ThermochemistrySolver solver(detailedConfig(&no_contract, "C", true)); }));
}

void testAdaptiveLadderNeedsElementTable() {
  const auto table = ph::MetalIonRateTable::parse(carbonFixture(true), "fixture", ph::MetalIonElementSet::parse("C"));
  ph::ThermochemistryOperatorConfig op;
  op.solver = detailedConfig(&table, "C", true);
  op.adaptive.enabled = true;
  op.units.velocity_code_to_cm_s = 1.0e5;
  op.units.length_code_to_cm = 3.0856775814913673e21;
  op.units.mass_code_to_g = 1.98847e43;
  op.units.density_code_to_cgs = op.units.mass_code_to_g / std::pow(op.units.length_code_to_cm, 3.0);
  const auto layout = ph::CompositionLayout::make({.network = ph::CompositionNetworkKind::kDetailedCarbonOxygen,
                                                   .metal_ion_element_mask = ph::MetalIonElementSet::parse("C").mask()});
  CHUI_REQUIRE(throws([&] { ph::ThermochemistryOperator chemistry(op, layout); }));
}

void testRealDataset(const std::string& path) {
  for (const char* subset : {"C", "O", "Si", "Fe", "C,O", "Si,Fe", "C,O,Si,Fe"}) {
    const auto set = ph::MetalIonElementSet::parse(subset);
    const auto table = ph::MetalIonRateTable::loadFromFile(path, set);
    CHUI_REQUIRE(!table.declaredVerified());
    CHUI_REQUIRE(table.energyConvention() == "collisional_ionization_work_in_cooling");
    for (std::size_t e = 0; e < ph::k_metal_ion_element_count; ++e) {
      const auto element = static_cast<ph::MetalIonElement>(e);
      if (!set.contains(element)) continue;
      CHUI_REQUIRE(table.hasIonizationEnergyContract(element));
    }
  }
  const auto all = ph::MetalIonElementSet::parse("C,O,Si,Fe");
  CHUI_REQUIRE(all.stageCount() == 58U);
  const auto table = ph::MetalIonRateTable::loadFromFile(path, all);
  CHUI_REQUIRE(std::abs(table.minimumTemperatureK() - 1.0e4) <= 1.0e-6 * 1.0e4);
  CHUI_REQUIRE(std::abs(table.maximumTemperatureK() - 1.0e8) <= 1.0e-6 * 1.0e8);
  bool clamped = false;
  for (std::size_t e = 0; e < ph::k_metal_ion_element_count; ++e) {
    const auto element = static_cast<ph::MetalIonElement>(e);
    const std::size_t stages = ph::metalIonStageCount(element);
    const auto bare = table.evaluate(element, stages - 1U, 1.0e6, 0.0, ph::MetalIonRangePolicy::kReject, clamped);
    CHUI_REQUIRE(bare.collisional_ionization_cm3_s == 0.0);
    const auto neutral = table.evaluate(element, 0U, 1.0e6, 0.0, ph::MetalIonRangePolicy::kReject, clamped);
    CHUI_REQUIRE(neutral.recombination_cm3_s == 0.0);
    CHUI_REQUIRE(neutral.cooling_erg_cm3_s > 0.0);
  }
  CHUI_REQUIRE(!clamped);
  CHUI_REQUIRE(throws([&] {
    (void)table.evaluate(ph::MetalIonElement::kC, 1U, 5.0e3, 0.0, ph::MetalIonRangePolicy::kReject, clamped);
  }));
  (void)table.evaluate(ph::MetalIonElement::kC, 1U, 5.0e3, 0.0, ph::MetalIonRangePolicy::kClampAndCount, clamped);
  CHUI_REQUIRE(clamped);
  cosmosim::core::PhysicsConfig physics;
  physics.chemistry_metal_ion_elements = "C,O,Si,Fe";
  physics.metal_ion_rate_table_path = path;
  physics.metal_ion_rate_table_sha256 = std::string(64, '0');
  CHUI_REQUIRE(throws([&] { (void)ph::loadBaryonicDataProducts(physics); }));
  physics.metal_ion_rate_table_sha256 = "633a66ee793e686ff932dd891fb7ac8f1e198924582a7fcf3b57ce5afa6c2d00";
  CHUI_REQUIRE(!throws([&] { (void)ph::loadBaryonicDataProducts(physics); }));
}

}  // namespace

int main() {
  testGrammarRefusals();
  testOptInGate();
  testAdaptiveLadderNeedsElementTable();
  if (const char* path = std::getenv("CHUI_C02_TABLE"); path != nullptr && *path != '\0') {
    testRealDataset(path);
  }
  return 0;
}
