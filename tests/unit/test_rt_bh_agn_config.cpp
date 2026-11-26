// CHUI-RT-BH-AGN-REF-v1 B04 admission tests (V04 disabled equivalence, V05/V07
// admission). NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <string>

#include "cosmosim/core/config.hpp"

namespace {

using cosmosim::core::ConfigError;
using cosmosim::core::loadFrozenConfigFromString;

const std::string k_base = "mode.mode = isolated_galaxy\n";

bool rejects(const std::string& text, const std::string& needle) {
  try {
    (void)loadFrozenConfigFromString(text, "rt_bh_agn_case");
  } catch (const ConfigError& error) {
    return std::string(error.what()).find(needle) != std::string::npos;
  }
  return false;
}

void testDisabledDefaultsKeepNormalizedTextUnchanged() {
  const auto frozen = loadFrozenConfigFromString(k_base, "defaults");
  CHUI_REQUIRE(frozen.normalized_text.find("rt_enabled") == std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("bh_engine_model") == std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("agn_feedback_model") == std::string::npos);
  CHUI_REQUIRE(!frozen.config.physics.rt.enabled);
  CHUI_REQUIRE(frozen.config.physics.bh_engine.model == cosmosim::core::BhEngineModelConfig::kLegacyHostCell);
}

void testInertKeysRejected() {
  CHUI_REQUIRE(rejects(k_base + "physics.rt_cfl = 0.2\n", "physics.rt_enabled=false"));
  CHUI_REQUIRE(rejects(k_base + "physics.bh_torque_epsilon = 0.1\n", "legacy_host_cell"));
  CHUI_REQUIRE(rejects(k_base + "physics.agn_thermal_efficiency = 0.01\n", "legacy_thermal_host"));
}

std::string referenceEngine() {
  return k_base +
      "physics.enable_black_hole_agn = true\n"
      "physics.bh_engine_model = reference_v1\n"
      "physics.bh_supply_model = bondi\n"
      "physics.bh_kernel_radius_min_kpc_phys = 0.1\n"
      "physics.bh_kernel_radius_max_kpc_phys = 1.0\n"
      "physics.bh_kernel_target_gas_mass_code = 1.0e-3\n"
      "physics.bh_disk_viscous_time_years = 1.0e6\n"
      "physics.agn_feedback_model = channels_v1\n"
      "physics.agn_thermal_efficiency = 0.015\n";
}

void testReferenceEngineAdmittedAndRoundTrips() {
  const auto frozen = loadFrozenConfigFromString(referenceEngine(), "reference_engine");
  CHUI_REQUIRE(frozen.config.physics.bh_engine.model == cosmosim::core::BhEngineModelConfig::kReferenceV1);
  CHUI_REQUIRE(frozen.normalized_text.find("bh_engine_model = reference_v1") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("agn_feedback_model = channels_v1") != std::string::npos);
  const auto again = loadFrozenConfigFromString(frozen.normalized_text, "round_trip");
  CHUI_REQUIRE(again.normalized_text == frozen.normalized_text);
}

void testReferenceEngineRefusals() {
  // Torque without explicit epsilon (no default calibration).
  std::string torque = referenceEngine();
  torque.replace(torque.find("bh_supply_model = bondi"), std::string("bh_supply_model = bondi").size(),
                 "bh_supply_model = bondi_torque");
  CHUI_REQUIRE(rejects(torque, "bh_torque_epsilon"));
  // Legacy-only feedback keys with the reference engine.
  CHUI_REQUIRE(rejects(referenceEngine() + "physics.bh_epsilon_f = 0.05\n", "legacy_host_cell-only"));
  // P1 models refused.
  CHUI_REQUIRE(rejects(referenceEngine() + "physics.bh_spin_model = kerr_thin_disk\n", "P1"));
  CHUI_REQUIRE(rejects(referenceEngine() + "physics.agn_magnetic_injection = true\n", "P1"));
  // Wind without a mass origin.
  CHUI_REQUIRE(rejects(referenceEngine() + "physics.agn_wind_efficiency = 0.01\nphysics.agn_wind_velocity_km_s = 1000\n",
                       "bh_disk_wind_loading"));
  // Superluminal-domain jet rejected, not capped.
  CHUI_REQUIRE(rejects(referenceEngine() + "physics.agn_jet_efficiency = 0.01\nphysics.agn_jet_velocity_km_s = 2.0e5\n",
                       "nonrelativistic"));
}

void testRadiationAdmission() {
  const std::string rt = k_base + "physics.rt_enabled = true\n";
  CHUI_REQUIRE(rejects(rt + "physics.rt_transport = vet\nphysics.rt_spectrum_upper_edge_ev = 100\n", "P1"));
  CHUI_REQUIRE(rejects(rt, "54.4"));  // ionizing3 needs an explicit upper edge
  CHUI_REQUIRE(rejects(rt + "physics.rt_spectrum_upper_edge_ev = 100\n", "FAIL_CLOSED"));  // photoionization
  const std::string grey = k_base +
      "physics.enable_cooling = false\n"
      "physics.rt_enabled = true\n"
      "physics.rt_spectral_profile = grey\n"
      "physics.rt_opacity_model = constant_reference\n"
      "physics.rt_constant_absorption_cm2_g = 0.4\n"
      "physics.rt_allow_reference_models = true\n"
      "physics.rt_boundary = outflow\n";
  const auto frozen = loadFrozenConfigFromString(grey, "grey_reference");
  CHUI_REQUIRE(frozen.config.physics.rt.enabled);
  CHUI_REQUIRE(frozen.normalized_text.find("rt_spectral_profile = grey") != std::string::npos);
  CHUI_REQUIRE(rejects(grey + "physics.rt_c_fraction = 1.5\n", "(0, 1]"));
}

}  // namespace

int main() {
  testDisabledDefaultsKeepNormalizedTextUnchanged();
  testInertKeysRejected();
  testReferenceEngineAdmittedAndRoundTrips();
  testReferenceEngineRefusals();
  testRadiationAdmission();
  return 0;
}
