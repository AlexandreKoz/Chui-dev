// CHUI Part 3 Checkpoints F/I — global adaptive-chemistry control law and
// four-mode Part 3 configuration decks.
// STATUS: written as source only. NOT RUN — USER-RESERVED.
//
// Invariants:
//  * nextChemistryPromotionBias is deterministic, bounded in [0, 4], raises
//    the detailed-region bias on large GLOBAL shadow disagreement, relaxes it
//    when the shadow agrees, and is unchanged without samples;
//  * every configs/chui_part3/<mode>_part3.param.txt deck parses, validates,
//    selects automatic regrid and its normalized dump re-parses identically
//    (the decks are meaningful inputs, not samples that fail at load).
#include <filesystem>
#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"
#include "../support/test_require.hpp"

#ifndef COSMOSIM_SOURCE_DIR
#define COSMOSIM_SOURCE_DIR "."
#endif

int main() {
  namespace ph = cosmosim::physics;
  namespace co = cosmosim::core;
  CHUI_REQUIRE(ph::nextChemistryPromotionBias(0.0, 0.5, 10U) == 0.5);
  CHUI_REQUIRE(ph::nextChemistryPromotionBias(3.9, 0.5, 10U) == 4.0);
  CHUI_REQUIRE(ph::nextChemistryPromotionBias(1.0, 0.001, 10U) == 0.75);
  CHUI_REQUIRE(ph::nextChemistryPromotionBias(0.1, 0.001, 10U) == 0.0);
  CHUI_REQUIRE(ph::nextChemistryPromotionBias(1.0, 0.05, 10U) == 1.0);
  CHUI_REQUIRE(ph::nextChemistryPromotionBias(1.0, 0.5, 0U) == 1.0);
  for (const char* mode : {"cosmo_cube", "zoom_in", "isolated_galaxy", "isolated_cluster"}) {
    const std::filesystem::path deck =
        std::filesystem::path(COSMOSIM_SOURCE_DIR) / "configs" / "chui_part3" / (std::string(mode) + "_part3.param.txt");
    const co::FrozenConfig frozen = co::loadFrozenConfigFromFile(deck);
    CHUI_REQUIRE(frozen.config.numerics.amr_regrid_mode == co::AmrRegridConfigMode::kAutomatic);
    CHUI_REQUIRE(frozen.config.numerics.global_step_max_retries == 2);
    const co::FrozenConfig reparsed = co::loadFrozenConfigFromString(frozen.normalized_text, "part3_deck_reparse");
    CHUI_REQUIRE(reparsed.normalized_text == frozen.normalized_text);
  }
  // Final-design decks (2026-11-04): recursive subcycling, data checksums
  // present, and the WSS09+FG20 deck is the ONLY one that acknowledges the
  // UVB mismatch. Parsing only; data files are not read by the config layer.
  for (const char* deck_name : {"cosmo_cube_final", "zoom_in_final", "isolated_galaxy_final", "isolated_cluster_final",
                                "cosmo_cube_wss09_experimental_final"}) {
    const std::filesystem::path deck =
        std::filesystem::path(COSMOSIM_SOURCE_DIR) / "configs" / "chui_final" / (std::string(deck_name) + ".param.txt");
    const co::FrozenConfig frozen = co::loadFrozenConfigFromFile(deck);
    CHUI_REQUIRE(frozen.config.numerics.amr_time_stepping == co::AmrTimeSteppingConfigMode::kRecursiveSubcycling);
    const auto& physics = frozen.config.physics;
    if (physics.stellar_population_model == co::StellarPopulationModelConfig::kChabrierYieldTableV1) {
      CHUI_REQUIRE(physics.ssp_yield_table_sha256.size() == 64U && physics.ssp_wind_table_sha256.size() == 64U);
    }
    CHUI_REQUIRE(physics.metal_cooling_accept_uv_background_mismatch ==
                 (std::string(deck_name) == "cosmo_cube_wss09_experimental_final"));
    const co::FrozenConfig reparsed = co::loadFrozenConfigFromString(frozen.normalized_text, "final_deck_reparse");
    CHUI_REQUIRE(reparsed.normalized_text == frozen.normalized_text);
  }
  return 0;
}
