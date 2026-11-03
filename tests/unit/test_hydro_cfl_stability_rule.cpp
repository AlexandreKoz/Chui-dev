// Scheme-consistent hydro CFL rule (Part 1 stabilization).
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Discrete stability: the directionally unsplit MUSCL-Hancock update
//   U^{n+1} = U^n - sum_d dt/(a dx_d) (F_{d,+} - F_{d,-})
// with first-order states reduces, for linear advection with speeds v_d >= 0,
// to unsplit donor cell whose amplification factor is
//   g(theta) = 1 - sum_d nu_d (1 - exp(-i theta_d)),  nu_d = v_d dt/(a dx_d).
// At theta = (pi, pi, pi), g = 1 - 2 sum_d nu_d, so |g| <= 1 requires
// sum_d nu_d <= 1. The previous per-axis rule (C = 0.4 per axis) admits
// sum_d nu_d = 1.2 for 3D diagonal flow (|g| = 1.4, unstable); the
// dimension-summed rule keeps sum_d nu_d <= 0.4. CTU keeps the per-axis rule.
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "../support/test_require.hpp"

namespace {

namespace co = cosmosim::core;

co::DirectionalCflTimeStepInput diagonal(co::HydroCflStabilityRule rule) {
  return co::DirectionalCflTimeStepInput{
      .cell_width_axis_code = {1.0, 1.0, 1.0},
      .velocity_axis_code = {1.0, 1.0, 1.0},
      .sound_speed_code = 0.0,
      .coordinate_frame = co::CoordinateFrame::kPhysical,
      .scale_factor = 1.0,
      .stability_rule = rule};
}

double worstAmplification(const std::array<double, 3>& nu) {
  // Maximum |g| over a coarse theta grid including the (pi,pi,pi) corner.
  double worst = 0.0;
  constexpr int k_samples = 8;
  for (int a = 0; a <= k_samples; ++a) {
    for (int b = 0; b <= k_samples; ++b) {
      for (int c = 0; c <= k_samples; ++c) {
        const std::array<double, 3> theta{
            3.141592653589793 * a / k_samples, 3.141592653589793 * b / k_samples,
            3.141592653589793 * c / k_samples};
        std::complex<double> g = 1.0;
        for (int d = 0; d < 3; ++d) {
          g -= nu[d] * (1.0 - std::exp(std::complex<double>(0.0, -theta[d])));
        }
        worst = std::max(worst, std::abs(g));
      }
    }
  }
  return worst;
}

void testDiagonalFlowPreviousRuleWasUnstable() {
  const double dt_axis = co::computeDirectionalCflTimeStep(diagonal(co::HydroCflStabilityRule::kPerAxisMaximum), 0.4);
  const double dt_sum = co::computeDirectionalCflTimeStep(diagonal(co::HydroCflStabilityRule::kDimensionSummed), 0.4);
  CHUI_REQUIRE(std::abs(dt_axis - 0.4) < 1.0e-15);
  CHUI_REQUIRE(std::abs(dt_sum - 0.4 / 3.0) < 1.0e-15);
  CHUI_REQUIRE(worstAmplification({dt_axis, dt_axis, dt_axis}) > 1.0 + 1.0e-9);   // old rule: unstable
  CHUI_REQUIRE(worstAmplification({dt_sum, dt_sum, dt_sum}) <= 1.0 + 1.0e-12);    // new rule: stable
}

void testOneDimensionalStepUnchanged() {
  // Only x active (y, z cells span the domain): both rules give 0.4 dx/(|v|+c).
  co::DirectionalCflTimeStepInput input = diagonal(co::HydroCflStabilityRule::kDimensionSummed);
  input.velocity_axis_code = {2.0, 0.0, 0.0};
  input.sound_speed_code = 1.0;
  input.active_axis = co::activeHydroAxes({1.0, 10.0, 10.0}, {32.0, 10.0, 10.0});
  CHUI_REQUIRE(input.active_axis[0] && !input.active_axis[1] && !input.active_axis[2]);
  CHUI_REQUIRE(std::abs(co::computeDirectionalCflTimeStep(input, 0.4) - 0.4 / 3.0) < 1.0e-15);
  input.stability_rule = co::HydroCflStabilityRule::kPerAxisMaximum;
  CHUI_REQUIRE(std::abs(co::computeDirectionalCflTimeStep(input, 0.4) - 0.4 / 3.0) < 1.0e-15);
}

void testAnisotropicCellsHighSoundSpeedAndZeroVelocity() {
  co::DirectionalCflTimeStepInput input = diagonal(co::HydroCflStabilityRule::kDimensionSummed);
  input.cell_width_axis_code = {1.0, 0.01, 1.0};
  input.velocity_axis_code = {0.0, 0.0, 0.0};
  input.sound_speed_code = 100.0;
  // rate = 100/1 + 100/0.01 + 100/1 = 10200
  CHUI_REQUIRE(std::abs(co::computeDirectionalCflTimeStep(input, 0.4) - 0.4 / 10200.0) < 1.0e-18);
  input.sound_speed_code = 0.0;  // static, cold: no signal, no limit
  CHUI_REQUIRE(std::isinf(co::computeDirectionalCflTimeStep(input, 0.4)));
}

void testComovingCrossingLength() {
  // Comoving widths, peculiar speeds: crossing length a dx, so dt scales with a.
  co::DirectionalCflTimeStepInput input = diagonal(co::HydroCflStabilityRule::kDimensionSummed);
  input.coordinate_frame = co::CoordinateFrame::kComoving;
  input.scale_factor = 0.25;
  CHUI_REQUIRE(std::abs(co::computeDirectionalCflTimeStep(input, 0.4) - 0.25 * 0.4 / 3.0) < 1.0e-15);
}

void testGuardUsesTheSameRule() {
  const auto input = diagonal(co::HydroCflStabilityRule::kDimensionSummed);
  // dt accepted under the old per-axis rule must now be rejected by the guard.
  const auto diagnostics = co::makeHydroCflDiagnostics(0U, input, 0.4, 0.4);
  CHUI_REQUIRE(diagnostics.stability_rule == co::HydroCflStabilityRule::kDimensionSummed);
  CHUI_REQUIRE(diagnostics.safety_factor < 1.0);
  bool rejected = false;
  try {
    co::assertHydroCflStable(diagnostics);
  } catch (const std::exception&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);
}

void testConfigMapping() {
  co::SimulationConfig config;
  CHUI_REQUIRE(co::hydroCflStabilityRuleFor(config) == co::HydroCflStabilityRule::kDimensionSummed);
  config.numerics.hydro_integration_scheme = co::HydroIntegrationMode::kCornerTransportUpwind;
  CHUI_REQUIRE(co::hydroCflStabilityRuleFor(config) == co::HydroCflStabilityRule::kPerAxisMaximum);
}

}  // namespace

int main() {
  testDiagonalFlowPreviousRuleWasUnstable();
  testOneDimensionalStepUnchanged();
  testAnisotropicCellsHighSoundSpeedAndZeroVelocity();
  testComovingCrossingLength();
  testGuardUsesTheSameRule();
  testConfigMapping();
  return 0;
}
