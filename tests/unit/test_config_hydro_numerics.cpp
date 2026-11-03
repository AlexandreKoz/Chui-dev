// Typed .param.txt keys introduced in Part 1 (hydro method selectors, dual
// energy, node memory admission): defaults, normalized-dump round trip and
// fail-closed validation.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
#include <stdexcept>
#include <string>

#include "cosmosim/core/config.hpp"
#include "../support/test_require.hpp"

namespace {

namespace co = cosmosim::core;

void testDefaultsPreserveHistoricalPath() {
  const auto frozen = co::loadFrozenConfigFromString("[mode]\nmode = zoom_in\n", "hydro_numerics_defaults");
  CHUI_REQUIRE(frozen.config.numerics.hydro_reconstruction == co::HydroReconstructionMode::kPiecewiseLinear);
  CHUI_REQUIRE(frozen.config.numerics.hydro_slope_limiter == co::HydroSlopeLimiterMode::kMonotonizedCentral);
  CHUI_REQUIRE(frozen.config.numerics.hydro_integration_scheme == co::HydroIntegrationMode::kMusclHancock);
  CHUI_REQUIRE(frozen.config.numerics.hydro_riemann_solver == co::HydroRiemannSolverMode::kHllc);
  CHUI_REQUIRE(frozen.config.numerics.hydro_dual_energy == co::HydroDualEnergyConfigMode::kOff);
  CHUI_REQUIRE(frozen.config.parallel.node_memory_admission == co::NodeMemoryAdmissionMode::kOff);
  CHUI_REQUIRE(frozen.normalized_text.find("hydro_reconstruction = piecewise_linear") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("node_memory_admission = off") != std::string::npos);
}

void testRoundTrip() {
  const auto frozen = co::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\nhydro_reconstruction = piecewise_parabolic\n"
      "hydro_slope_limiter = van_leer\nhydro_integration_scheme = ctu\nhydro_riemann_solver = hlle\n"
      "hydro_dual_energy = thermal_fraction\nhydro_dual_energy_eta = 0.01\n"
      "[parallel]\nnode_memory_admission = shared_equal\nnode_memory_budget_bytes = 34359738368\n"
      "node_memory_reserve_bytes = 1073741824\n",
      "hydro_numerics_roundtrip");
  const auto reparsed = co::loadFrozenConfigFromString(frozen.normalized_text, "hydro_numerics_reparse");
  for (const auto* config : {&frozen.config, &reparsed.config}) {
    CHUI_REQUIRE(config->numerics.hydro_reconstruction == co::HydroReconstructionMode::kPiecewiseParabolic);
    CHUI_REQUIRE(config->numerics.hydro_slope_limiter == co::HydroSlopeLimiterMode::kVanLeer);
    CHUI_REQUIRE(config->numerics.hydro_integration_scheme == co::HydroIntegrationMode::kCornerTransportUpwind);
    CHUI_REQUIRE(config->numerics.hydro_riemann_solver == co::HydroRiemannSolverMode::kHlle);
    CHUI_REQUIRE(config->numerics.hydro_dual_energy == co::HydroDualEnergyConfigMode::kThermalFraction);
    CHUI_REQUIRE(config->numerics.hydro_dual_energy_eta == 0.01);
    CHUI_REQUIRE(config->parallel.node_memory_admission == co::NodeMemoryAdmissionMode::kSharedEqual);
    CHUI_REQUIRE(config->parallel.node_memory_budget_bytes == 34359738368ULL);
    CHUI_REQUIRE(config->parallel.node_memory_reserve_bytes == 1073741824ULL);
  }
}

bool rejects(const std::string& text) {
  try {
    (void)co::loadFrozenConfigFromString(text, "hydro_numerics_invalid");
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

void testFailClosedValidation() {
  CHUI_REQUIRE(rejects("[mode]\nmode = zoom_in\n[numerics]\nhydro_reconstruction = weno5\n"));
  CHUI_REQUIRE(rejects("[mode]\nmode = zoom_in\n[parallel]\nnode_memory_admission = dynamic\n"));
  const auto pool = co::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[parallel]\nnode_memory_admission = shared_pool\n", "node_pool_mode");
  CHUI_REQUIRE(pool.config.parallel.node_memory_admission == co::NodeMemoryAdmissionMode::kSharedPool);
  CHUI_REQUIRE(pool.normalized_text.find("node_memory_admission = shared_pool") != std::string::npos);
  CHUI_REQUIRE(rejects("[mode]\nmode = zoom_in\n[numerics]\nhydro_integration_scheme = rk3\n"));
  CHUI_REQUIRE(rejects("[mode]\nmode = zoom_in\n[numerics]\nhydro_riemann_solver = roe\n"));
  CHUI_REQUIRE(rejects("[mode]\nmode = zoom_in\n[numerics]\nhydro_dual_energy_eta = 1.5\n"));
  // A node budget without the admission mode would be silently ignored: refused.
  CHUI_REQUIRE(rejects("[mode]\nmode = zoom_in\n[parallel]\nnode_memory_budget_bytes = 1024\n"));
  CHUI_REQUIRE(rejects(
      "[mode]\nmode = zoom_in\n[parallel]\nnode_memory_admission = shared_equal\n"
      "node_memory_budget_bytes = 1024\nnode_memory_reserve_bytes = 1024\n"));
}

}  // namespace

int main() {
  testDefaultsPreserveHistoricalPath();
  testRoundTrip();
  testFailClosedValidation();
  return 0;
}
