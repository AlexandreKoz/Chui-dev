// CHUI-EXT-PHYS-REF-v1 gates V00/VI (I02): the example decks of
// configs/extended_physics load, normalize and round-trip; disabled sections
// emit no text; every refused capability names its missing edge; inert keys
// are rejected. NOT RUN in pass 01 (source-only campaign).
#include "../support/test_require.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "cosmosim/core/config.hpp"

namespace c = cosmosim::core;

namespace {

const std::filesystem::path k_decks = std::filesystem::path(COSMOSIM_SOURCE_DIR) / "configs" / "extended_physics";

std::string readFile(const std::filesystem::path& path) {
  std::ifstream in(path);
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

bool rejects(const std::string& text, const std::string& needle) {
  try {
    (void)c::loadFrozenConfigFromString(text, "ext_phys_case");
  } catch (const c::ConfigError& error) {
    return std::string(error.what()).find(needle) != std::string::npos;
  }
  return false;
}

// EXT-I02-T01: decks admitted (the linear-response deck is admitted at config
// level; its tables are verified when the gravity runtime is built).
void testDecks() {
  for (const char* name : {"neutrino_background_dmo.param.txt", "neutrino_linear_response_dmo.param.txt",
                           "isolated_galaxy_hi_cube_rime.param.txt"}) {
    const auto frozen = c::loadFrozenConfigFromFile(k_decks / name);
    const auto again = c::loadFrozenConfigFromString(frozen.normalized_text, "round_trip");
    CHUI_REQUIRE(again.normalized_text == frozen.normalized_text);
  }
  const auto hi = c::loadFrozenConfigFromFile(k_decks / "isolated_galaxy_hi_cube_rime.param.txt");
  CHUI_REQUIRE(hi.config.extended.observables.hi_cube && hi.config.extended.instruments.enabled);
  CHUI_REQUIRE(hi.normalized_text.find("[instruments]") != std::string::npos);
}

// EXT-I02-T02 / V00: disabled sections and refusals.
void testRefusals() {
  const std::string base = "mode.mode = isolated_galaxy\n";
  const auto legacy = c::loadFrozenConfigFromString(base, "legacy");
  for (const char* section : {"[neutrinos]", "[cosmic_rays]", "[observables]", "[instruments]"}) {
    CHUI_REQUIRE(legacy.normalized_text.find(section) == std::string::npos);
  }
  CHUI_REQUIRE(rejects(base + "observables.forest_pixels = 64\n", "would have no effect"));
  CHUI_REQUIRE(rejects(base + "instruments.sefd_jy = 10\n", "would have no effect"));
  CHUI_REQUIRE(rejects(base + "cosmic_rays.model = spectral\n", "mhd.enabled=true"));
  CHUI_REQUIRE(rejects(base + "observables.enabled = true\nobservables.interval_steps = 1\n"
                              "observables.sightlines_per_side = 4\nobservables.sz_maps = true\n", "legacy cooling"));
  // Part 2 base (species lanes) from the repository deck.
  const std::string obs = readFile(std::filesystem::path(COSMOSIM_SOURCE_DIR) / "configs" / "chui_part2" /
                                   "isolated_galaxy_part2.param.txt") +
                          "\nobservables.enabled = true\nobservables.interval_steps = 1\nobservables.sightlines_per_side = 4\n";
  CHUI_REQUIRE(rejects(obs, "at least one product"));
  CHUI_REQUIRE(rejects(obs + "observables.lyman_alpha_forest = true\nobservables.forest_pixels = 64\n", "cosmological"));
  CHUI_REQUIRE(rejects(obs + "observables.hi_cube = true\nobservables.hi_channels = 8\n"
                             "observables.hi_velocity_half_width_km_s = 100\ninstruments.enabled = true\n"
                             "instruments.backend = accelerated\n", "FAIL_CLOSED"));
}

}  // namespace

int main() {
  testDecks();
  testRefusals();
  return 0;
}
