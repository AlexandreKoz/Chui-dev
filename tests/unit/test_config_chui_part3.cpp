// CHUI Part 3 configuration keys: typed validation, normalized dump round
// trip and refusals.
// STATUS: written as source only. NOT RUN — USER-RESERVED.
//
// Keys: numerics.amr_regrid_mode / _interval_steps / _max_operations /
// amr_refine_cell_mass_code / amr_refine_jump_threshold /
// amr_refine_jeans_cells / amr_derefine_hysteresis,
// numerics.global_step_max_retries, physics.ssp_wind_model /
// ssp_wind_table_path / ssp_wind_table_sha256, physics.ssp_association_*.
#include <stdexcept>
#include <string>

#include "cosmosim/core/config.hpp"
#include "../support/test_require.hpp"

namespace {

namespace co = cosmosim::core;

bool rejects(const std::string& text) {
  try {
    (void)co::loadFrozenConfigFromString(text, "part3_reject");
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

void testDefaultsKeepHistoricalBehaviour() {
  const auto frozen = co::loadFrozenConfigFromString("[mode]\nmode = cosmo_cube\n", "part3_defaults");
  CHUI_REQUIRE(frozen.config.numerics.amr_regrid_mode == co::AmrRegridConfigMode::kOff);
  CHUI_REQUIRE(frozen.config.numerics.global_step_max_retries == 0);
  CHUI_REQUIRE(frozen.config.physics.ssp_wind_model == co::SspWindModelConfig::kNone);
  for (const char* line : {"amr_regrid_mode = off", "global_step_max_retries = 0", "amr_refine_jeans_cells = 4",
                           "ssp_wind_model = none", "ssp_association_radius_code = 0"}) {
    CHUI_REQUIRE(frozen.normalized_text.find(line) != std::string::npos);
  }
}

void testRegridRoundTrip() {
  const std::string text =
      "[mode]\nmode = cosmo_cube\n[numerics]\namr_regrid_mode = automatic\namr_regrid_interval_steps = 4\n"
      "amr_regrid_max_operations = 16\namr_refine_jump_threshold = 0.3\namr_refine_jeans_cells = 8\n"
      "amr_derefine_hysteresis = 0.25\namr_max_level = 6\nglobal_step_max_retries = 3\n";
  const auto frozen = co::loadFrozenConfigFromString(text, "part3_regrid");
  const auto reparsed = co::loadFrozenConfigFromString(frozen.normalized_text, "part3_regrid_reparse");
  for (const auto* c : {&frozen.config, &reparsed.config}) {
    CHUI_REQUIRE(c->numerics.amr_regrid_mode == co::AmrRegridConfigMode::kAutomatic);
    CHUI_REQUIRE(c->numerics.amr_regrid_interval_steps == 4);
    CHUI_REQUIRE(c->numerics.amr_regrid_max_operations == 16);
    CHUI_REQUIRE(c->numerics.amr_refine_jump_threshold == 0.3);
    CHUI_REQUIRE(c->numerics.amr_refine_jeans_cells == 8.0);
    CHUI_REQUIRE(c->numerics.amr_derefine_hysteresis == 0.25);
    CHUI_REQUIRE(c->numerics.global_step_max_retries == 3);
  }
  CHUI_REQUIRE(frozen.normalized_text == reparsed.normalized_text);
}

void testRefusals() {
  const std::string base = "[mode]\nmode = cosmo_cube\n[numerics]\n";
  // automatic regrid without any criterion
  CHUI_REQUIRE(rejects(base + "amr_regrid_mode = automatic\namr_refine_jeans_cells = 0\n"));
  CHUI_REQUIRE(rejects(base + "amr_regrid_interval_steps = 0\n"));
  CHUI_REQUIRE(rejects(base + "amr_derefine_hysteresis = 1.5\n"));
  CHUI_REQUIRE(rejects(base + "amr_refine_jump_threshold = -1\n"));
  CHUI_REQUIRE(rejects(base + "amr_regrid_mode = sometimes\n"));
  CHUI_REQUIRE(rejects(base + "global_step_max_retries = 99\n"));
  CHUI_REQUIRE(rejects(base + "global_step_max_retries = 2\nhierarchical_max_rung = 2\n"));
  const std::string galaxy = "[mode]\nmode = isolated_galaxy\n[physics]\nuv_background_model = none\n";
  // associations need a radius and a Part 2 population model
  CHUI_REQUIRE(rejects(galaxy + "ssp_clustered_associations = true\n"));
  // wind table needs a path + checksum and a Part 2 population model
  CHUI_REQUIRE(rejects(galaxy + "ssp_wind_model = table_v1\n"));
  CHUI_REQUIRE(rejects(galaxy + "ssp_wind_table_path = /tmp/wind.txt\n"));
  CHUI_REQUIRE(rejects(galaxy + "ssp_wind_model = nonsense\n"));
}

// R1 (A03): physics-aware AMR indicator keys round-trip, default to 0 and
// refuse unsupported combinations explicitly.
void testPhysicsIndicatorKeys() {
  const std::string base = "[mode]\nmode = cosmo_cube\n[numerics]\namr_regrid_mode = automatic\namr_max_level = 6\n";
  const auto frozen = co::loadFrozenConfigFromString(
      base + "amr_refine_metallicity_jump = 0.5\namr_refine_vorticity = 0.25\namr_refine_jeans_cells = 0\n",
      "r1_keys");
  const auto reparsed = co::loadFrozenConfigFromString(frozen.normalized_text, "r1_keys_reparse");
  for (const auto* c : {&frozen.config, &reparsed.config}) {
    CHUI_REQUIRE(c->numerics.amr_refine_metallicity_jump == 0.5);
    CHUI_REQUIRE(c->numerics.amr_refine_vorticity == 0.25);
    CHUI_REQUIRE(c->numerics.amr_refine_cooling_length_cells == 0.0);
    CHUI_REQUIRE(c->numerics.amr_refine_species_jump == 0.0);
  }
  CHUI_REQUIRE(frozen.normalized_text == reparsed.normalized_text);
  CHUI_REQUIRE(rejects(base + "amr_refine_vorticity = -1\n"));
  // Cooling length needs the legacy in-solver cooling provider.
  CHUI_REQUIRE(rejects(base + "amr_refine_cooling_length_cells = 4\n[physics]\nthermochemistry_backend = fixed_hybrid\n"
                              "thermochemistry_allow_unqualified_rates = true\nuv_background_model = none\n"));
  CHUI_REQUIRE(rejects(base + "amr_refine_cooling_length_cells = 4\n[physics]\nenable_cooling = false\n"));
  // Species interfaces need a composition network backend.
  CHUI_REQUIRE(rejects(base + "amr_refine_species_jump = 0.1\n"));
  // R2: event lookahead needs automatic regrid and a Part 2 population model.
  CHUI_REQUIRE(rejects(base + "amr_predictive_event_lookahead_years = 1.0e7\n"));
  CHUI_REQUIRE(rejects("[mode]\nmode = cosmo_cube\n[numerics]\namr_predictive_event_lookahead_years = 1.0e7\n"));
}

}  // namespace

int main() {
  testPhysicsIndicatorKeys();
  testDefaultsKeepHistoricalBehaviour();
  testRegridRoundTrip();
  testRefusals();
  return 0;
}
