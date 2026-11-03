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
      "thermochemistry_allow_unqualified_rates = true\n"
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

void testUvBackgroundMismatchKeyRoundTrip() {
  const std::string text = k_base +
      "metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
      "thermochemistry_allow_unqualified_rates = true\nmetal_cooling_mode = element_table\n"
      "metal_cooling_table_path = /data/wss09_chui_v1.txt\n"
      "metal_cooling_table_sha256 = " + std::string(64, 'a') + "\n"
      "metal_cooling_accept_uv_background_mismatch = true\n";
  const auto frozen = co::loadFrozenConfigFromString(text, "part25_uvb");
  CHUI_REQUIRE(frozen.config.physics.metal_cooling_accept_uv_background_mismatch);
  CHUI_REQUIRE(frozen.normalized_text.find("metal_cooling_accept_uv_background_mismatch = true") != std::string::npos);
  const auto reparsed = co::loadFrozenConfigFromString(frozen.normalized_text, "part25_uvb_reparse");
  CHUI_REQUIRE(reparsed.normalized_text == frozen.normalized_text);
  const auto defaults = co::loadFrozenConfigFromString("[mode]\nmode = zoom_in\n", "part25_defaults");
  CHUI_REQUIRE(!defaults.config.physics.metal_cooling_accept_uv_background_mismatch);
  CHUI_REQUIRE(defaults.normalized_text.find("metal_cooling_accept_uv_background_mismatch = false") !=
               std::string::npos);
}

void testFailClosedCombinations() {
  // core_elements without a Part 2 backend and vice versa.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\n"));
  CHUI_REQUIRE(rejects("thermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"));
  // Legacy analytic UV proxies refused by Part 2 backends.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "uv_background_model = hm12\n"));
  // Treecool without a path.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "uv_background_model = treecool_table\n"));
  // Table metal cooling without a file.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "metal_cooling_mode = element_table\n"));
  // Dust needs a network backend.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = equilibrium_primordial\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "dust_model = two_reservoir_v1\n"));
  // Adaptive selection needs the detailed layout.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "chemistry_adaptive_selection = local_v1\n"));
  // Effective ISM + Part 2 needs the regime manager (no double cooling).
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "star_formation_model = effective_multiphase_tng_like\n"));
  // Yield-table population needs a table; mechanical feedback needs a Part 2 population.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "stellar_population_model = chabrier_yield_table_v1\n"
                       "fb_supernova_model = mechanical_resolved_unresolved\n"));
  CHUI_REQUIRE(rejects("fb_supernova_model = mechanical_resolved_unresolved\n"));
  // Part 2 population with feedback on requires the mechanical deposition.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "stellar_population_model = chabrier_reference_v1\n"));
  // Unqualified rates require the explicit opt-in; clustered associations
  // are not implemented and are refused.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = equilibrium_primordial\n"));
  CHUI_REQUIRE(rejects("ssp_clustered_associations = true\n"));
  // Part 2.5: a data path without its SHA-256 (or with a malformed one) is
  // refused, so a replaced file cannot silently change a restarted run.
  for (const std::string sha : {std::string(), std::string("ABC"), std::string(64, 'G')}) {
    CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                         "thermochemistry_allow_unqualified_rates = true\nmetal_cooling_mode = element_table\n"
                         "metal_cooling_table_path = /data/t.txt\nmetal_cooling_table_sha256 = " + sha + "\n"));
  }
  {
    // TREECOOL path without checksum (own base: k_base already sets the UVB model).
    const std::string treecool = "[mode]\nmode = isolated_galaxy\n[physics]\nmetal_species_mode = core_elements\n"
                                 "thermochemistry_backend = fixed_hybrid\nthermochemistry_allow_unqualified_rates = true\n"
                                 "uv_background_model = treecool_table\nuv_background_table_path = /data/treecool\n";
    bool refused = false;
    try { (void)co::loadFrozenConfigFromString(treecool, "part25_treecool_nosha"); }
    catch (const std::exception&) { refused = true; }
    CHUI_REQUIRE(refused);
    const auto ok = co::loadFrozenConfigFromString(
        treecool + "uv_background_table_sha256 = " + std::string(64, '0') + "\n", "part25_treecool_sha");
    CHUI_REQUIRE(ok.config.physics.uv_background_table_sha256 == std::string(64, '0'));
  }
  // Part 2.5: UVB-mismatch acceptance is a no-op without an element table.
  CHUI_REQUIRE(rejects("metal_species_mode = core_elements\nthermochemistry_backend = fixed_hybrid\n"
                       "thermochemistry_allow_unqualified_rates = true\n"
                       "metal_cooling_accept_uv_background_mismatch = true\n"));
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
  testUvBackgroundMismatchKeyRoundTrip();
  return 0;
}
