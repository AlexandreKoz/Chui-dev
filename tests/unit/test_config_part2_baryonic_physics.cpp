// CHUI Part 2 .param.txt keys: defaults keep Part 1 behaviour, normalized
// dump round trip, fail-closed validation of incoherent model combinations.
// STATUS: written as source only, NOT RUN (user-reserved).
#include <stdexcept>
#include <string>

#include "cosmosim/core/config.hpp"
#include "../support/test_require.hpp"

namespace {

namespace co = cosmosim::core;

const std::string k_base = "[mode]\nmode = isolated_galaxy\n[physics]\nuv_background_model = none\n";

bool rejects(const std::string& physics_lines) {
  try {
    (void)co::loadFrozenConfigFromString(k_base + physics_lines, "part2_reject");
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

void testDefaultsAreLegacy() {
  const auto frozen = co::loadFrozenConfigFromString("[mode]\nmode = zoom_in\n", "part2_defaults");
  CHUI_REQUIRE(frozen.config.physics.thermochemistry_backend == co::ThermochemistryBackendMode::kLegacyCooling);
  CHUI_REQUIRE(frozen.config.physics.stellar_population_model == co::StellarPopulationModelConfig::kLegacyTable);
  CHUI_REQUIRE(frozen.config.physics.fb_supernova_model == co::SupernovaFeedbackModelConfig::kLegacy);
  CHUI_REQUIRE(frozen.normalized_text.find("thermochemistry_backend = legacy_cooling") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("dust_model = none") != std::string::npos);
}

void testRoundTripDetailedFullPhysics() {
  const std::string text = k_base +
      "metal_species_mode = core_elements\nthermochemistry_backend = fixed_detailed\n"
      "metal_cooling_mode = fine_structure_cold\ndust_model = two_reservoir_v1\n"
      "sgs_turbulence_model = budgeted_v1\nchemistry_adaptive_selection = local_v1\n"
      "stellar_population_model = chabrier_reference_v1\nfb_supernova_model = mechanical_resolved_unresolved\n"
      "isrf_g0_habing = 1.5\nthermochemistry_max_substeps = 256\n";
  const auto frozen = co::loadFrozenConfigFromString(text, "part2_roundtrip");
  const auto reparsed = co::loadFrozenConfigFromString(frozen.normalized_text, "part2_reparse");
  for (const auto* c : {&frozen.config, &reparsed.config}) {
    CHUI_REQUIRE(c->physics.thermochemistry_backend == co::ThermochemistryBackendMode::kFixedDetailed);
    CHUI_REQUIRE(c->physics.dust_model == co::DustModelConfig::kTwoReservoirV1);
    CHUI_REQUIRE(c->physics.sgs_turbulence_model == co::SgsTurbulenceModelConfig::kBudgetedV1);
    CHUI_REQUIRE(c->physics.chemistry_adaptive_selection == co::ChemistryAdaptiveSelectionConfig::kLocalV1);
    CHUI_REQUIRE(c->physics.isrf_g0_habing == 1.5);
    CHUI_REQUIRE(c->physics.thermochemistry_max_substeps == 256U);
  }
  CHUI_REQUIRE(frozen.normalized_text == reparsed.normalized_text);
}

void testFailClosedCombinations() {
  // core_elements without a Part 2 backend and vice versa.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\n"));
  CHUI_REQUIRE(rejects("thermochemistry_backend = fixed_hybrid\n"));
  // Legacy analytic UV proxies refused by Part 2 backends.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "uv_background_model = hm12\n"));
  // Treecool without a path.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "uv_background_model = treecool_table\n"));
  // Table metal cooling without a file.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "metal_cooling_mode = element_table\n"));
  // Dust needs a network backend.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = equilibrium_primordial\n"
                       "dust_model = two_reservoir_v1\n"));
  // Adaptive selection needs the detailed layout.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "chemistry_adaptive_selection = local_v1\n"));
  // Effective ISM + Part 2 needs the regime manager (no double cooling).
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "star_formation_model = effective_multiphase_tng_like\n"));
  // Yield-table population needs a table; mechanical feedback needs a Part 2 population.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "stellar_population_model = chabrier_yield_table_v1\n"
                       "fb_supernova_model = mechanical_resolved_unresolved\n"));
  CHUI_REQUIRE(rejects("fb_supernova_model = mechanical_resolved_unresolved\n"));
  // Part 2 population with feedback on requires the mechanical deposition.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "stellar_population_model = chabrier_reference_v1\n"));
  // Numerical ranges.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_relative_tolerance = 0.5\n"));
  CHUI_REQUIRE(rejects("ssp_cc_mass_min_msun = 50\nssp_cc_mass_max_msun = 40\n"));
}

}  // namespace

int main() {
  testDefaultsAreLegacy();
  testRoundTripDetailedFullPhysics();
  testFailClosedCombinations();
  return 0;
}
