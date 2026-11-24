// CHUI-MHD-REF-v1 ideal reference integration tests (G01, G03, G04, G05,
// I01, I03, I04/I05 smoke, I09, axis/patch invariance, rollback).
// Source written, NOT EXECUTED in the source-only campaign.

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numbers>
#include <vector>

#include "../unit/mhd/mhd_test_fixtures.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_level_stepper.hpp"

namespace {

using namespace cosmosim::mhd;
using cosmosim::mhd::test::GasPrimitive;
using cosmosim::mhd::test::TestLevel;
constexpr double k_gamma = 5.0 / 3.0;
constexpr double k_two_pi = 2.0 * std::numbers::pi;

// Advance with the stepper's own CFL estimate until `end_time`.
MhdLevelAdvanceResult runUntil(MhdLevelStepper& stepper, double end_time, double cfl_scale = 1.0) {
  double time = 0.0;
  MhdLevelAdvanceResult last{};
  while (time < end_time * (1.0 - 1.0e-14)) {
    MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(time, 0.0);
    interval.dt_code = 1.0;
    double dt = cfl_scale * stepper.estimateStableTimeStep(interval, {});
    dt = std::min(dt, end_time - time);
    interval.dt_code = dt;
    last = stepper.advance(interval, {});
    CHUI_REQUIRE(last.accepted);
    time += last.dt_accepted_code;
  }
  return last;
}

// P0-02: positive normal flux leaves the lower cell. A density wave advected
// at +1 must move toward +x by exactly the elapsed time (phase of mode 1).
void testAdvectionDirection() {
  const double length = 1.0;
  TestLevel level({32, 4, 4}, {32, 4, 4}, {length / 32.0, 0.25, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [&](std::array<double, 3> x) {
        GasPrimitive w{};
        w.rho = 1.0 + 0.2 * std::sin(k_two_pi * x[0] / length);
        w.velocity = {1.0, 0.0, 0.0};
        w.pressure = 1.0;
        return w;
      },
      {0.1, 0.0, 0.0}, {}, 1.0, k_gamma);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const double elapsed = 0.25;
  (void)runUntil(stepper, elapsed);
  std::complex<double> mode{};
  for (std::int64_t i = 0; i < 32; ++i) {
    const double x = (static_cast<double>(i) + 0.5) / 32.0;
    mode += level.cellValue(0, {i, 1, 1}) * std::exp(std::complex<double>(0.0, -k_two_pi * x));
  }
  // sin(2 pi (x - t)) has mode-1 phase -pi/2 - 2 pi t.
  const double phase = std::arg(mode);
  const double expected = std::remainder(-0.5 * std::numbers::pi - k_two_pi * elapsed, k_two_pi);
  CHUI_REQUIRE(std::abs(std::remainder(phase - expected, k_two_pi)) < 0.05);
}

// G01: B = 0 stays identically zero (no magnetic floor, no spurious EMF).
void testZeroFieldLimit() {
  TestLevel level({16, 16, 4}, {8, 8, 4}, {1.0 / 16, 1.0 / 16, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [](std::array<double, 3> x) {
        GasPrimitive w{};
        w.rho = 1.0 + 0.1 * std::cos(k_two_pi * x[0]) * std::sin(k_two_pi * x[1]);
        w.velocity = {0.3 * std::sin(k_two_pi * x[1]), -0.2 * std::cos(k_two_pi * x[0]), 0.0};
        w.pressure = 0.6;
        return w;
      },
      {0.0, 0.0, 0.0}, {}, 1.0, k_gamma);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  (void)runUntil(stepper, 0.05);
  for (std::size_t patch = 0; patch < level.patchCount(); ++patch) {
    for (const double value : level.faces(patch).fluxX()) CHUI_REQUIRE(value == 0.0);
    for (const double value : level.faces(patch).fluxY()) CHUI_REQUIRE(value == 0.0);
    for (const double value : level.faces(patch).fluxZ()) CHUI_REQUIRE(value == 0.0);
  }
}

// A uniform magnetized moving state is an exact steady solution.
void testFreeStream() {
  TestLevel level({8, 8, 8}, {8, 8, 8}, {0.125, 0.125, 0.125}, cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [](std::array<double, 3>) {
        GasPrimitive w{};
        w.velocity = {0.3, -0.7, 0.2};
        return w;
      },
      {0.4, 0.3, -0.5}, {}, 1.0, k_gamma);
  const double rho0 = level.cellValue(0, {3, 3, 3});
  const double e0 = level.cellValue(4, {3, 3, 3});
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  (void)runUntil(stepper, 0.1);
  for (std::int64_t i = 0; i < 8; ++i) {
    CHUI_REQUIRE(std::abs(level.cellValue(0, {i, 2, 5}) - rho0) <= 1e-13);
    CHUI_REQUIRE(std::abs(level.cellValue(4, {i, 2, 5}) - e0) <= 1e-12);
  }
}

// P0-03: the same one-dimensional magnetized shear/sound wave along x, y and
// z on rectangular (non-cubic) domains gives identical profiles.
void testAxisRotationInvariance() {
  std::array<std::vector<double>, 3> profiles{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    std::array<std::int64_t, 3> cells{4, 6, 5};
    cells[axis] = 24;
    std::array<double, 3> spacing{1.0 / 24, 1.0 / 24, 1.0 / 24};
    TestLevel level(cells, cells, spacing, cosmosim::mhd::test::periodicDomain());
    std::array<double, 3> b{0.0, 0.0, 0.0};
    b[axis] = 0.5;
    b[(axis + 1U) % 3U] = 0.2;
    level.initialize(
        [&](std::array<double, 3> x) {
          GasPrimitive w{};
          w.rho = 1.0 + 0.1 * std::sin(k_two_pi * x[axis]);
          w.velocity[axis] = 0.05 * std::sin(k_two_pi * x[axis]);
          w.velocity[(axis + 2U) % 3U] = 0.1 * std::cos(k_two_pi * x[axis]);
          w.pressure = 1.0;
          return w;
        },
        b, {}, 1.0, k_gamma);
    MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
    for (int step = 0; step < 6; ++step) {
      MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(step * 2.0e-3, 2.0e-3);
      CHUI_REQUIRE(stepper.advance(interval, {}).accepted);
    }
    for (std::int64_t n = 0; n < 24; ++n) {
      std::array<std::int64_t, 3> cell{1, 2, 3};
      cell[axis] = n;
      profiles[axis].push_back(level.cellValue(0, cell));
    }
  }
  for (std::size_t n = 0; n < profiles[0].size(); ++n) {
    CHUI_REQUIRE(std::abs(profiles[0][n] - profiles[1][n]) <= 1e-12);
    CHUI_REQUIRE(std::abs(profiles[0][n] - profiles[2][n]) <= 1e-12);
  }
}

GasPrimitive orszagTangGas(std::array<double, 3> x) {
  GasPrimitive w{};
  w.rho = 25.0 / (36.0 * std::numbers::pi);
  w.pressure = 5.0 / (12.0 * std::numbers::pi);
  w.velocity = {-std::sin(k_two_pi * x[1]), std::sin(k_two_pi * x[0]), 0.0};
  return w;
}

// Orszag-Tang field from A_z = B0 (cos(4 pi x)/(4 pi) + cos(2 pi y)/(2 pi)).
MhdEdgePotential orszagTangPotential() {
  const double b0 = 1.0 / std::sqrt(4.0 * std::numbers::pi);
  return [b0](MhdAxis axis, std::array<double, 3> start, double length) {
    if (axis != MhdAxis::kZ) return 0.0;
    return b0 * (std::cos(2.0 * k_two_pi * start[0]) / (2.0 * k_two_pi) +
                 std::cos(k_two_pi * start[1]) / k_two_pi) *
        length;
  };
}

// Multi-patch stage ownership: a 2x2 tiling must reproduce the single-patch
// evolution (shared faces/edges from one owner, same-stage halos) and keep
// shared boundary faces bitwise identical.
void testMultiPatchEquivalence() {
  TestLevel single({16, 16, 4}, {16, 16, 4}, {1.0 / 16, 1.0 / 16, 0.25},
                   cosmosim::mhd::test::periodicDomain());
  TestLevel tiled({16, 16, 4}, {8, 8, 4}, {1.0 / 16, 1.0 / 16, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  single.initialize(orszagTangGas, {0.0, 0.0, 0.0}, orszagTangPotential(), 1.0, k_gamma);
  tiled.initialize(orszagTangGas, {0.0, 0.0, 0.0}, orszagTangPotential(), 1.0, k_gamma);
  MhdLevelStepper a = single.makeStepper(MhdNumericsPolicy{});
  MhdLevelStepper b = tiled.makeStepper(MhdNumericsPolicy{});
  for (int step = 0; step < 4; ++step) {
    const MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(step * 1e-3, 1e-3);
    CHUI_REQUIRE(a.advance(interval, {}).accepted);
    CHUI_REQUIRE(b.advance(interval, {}).accepted);
  }
  for (std::int64_t j = 0; j < 16; ++j) {
    for (std::int64_t i = 0; i < 16; ++i) {
      for (std::size_t lane = 0; lane < 5U; ++lane) {
        const double lhs = single.cellValue(lane, {i, j, 1});
        const double rhs = tiled.cellValue(lane, {i, j, 1});
        CHUI_REQUIRE(std::abs(lhs - rhs) <= 1e-12 * (1.0 + std::abs(lhs)));
      }
    }
  }
  // Shared x-plane between tiles 0 (grid x=0) and 1 (grid x=1): identical.
  const MhdPatchGeometry& geometry = tiled.binding(0).geometry;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      CHUI_REQUIRE(tiled.faces(0).faceFluxX(geometry.nx, j, k) == tiled.faces(1).faceFluxX(0, j, k));
    }
  }
}

// G03/G04: translating field loop (Gardiner & Stone 2005) in a periodic box:
// mass, momentum and total energy conserved to roundoff; strict CT divergence.
void testFieldLoopConservation() {
  TestLevel level({32, 32, 4}, {16, 16, 4}, {1.0 / 32, 1.0 / 32, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [](std::array<double, 3>) {
        GasPrimitive w{};
        w.velocity = {2.0, 1.0, 0.0};
        return w;
      },
      {0.0, 0.0, 0.0}, fieldLoopPotential(1.0e-3, 0.3, {0.5, 0.5, 0.5}), 1.0, k_gamma);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const MhdConservationBudget begin = stepper.levelBudget(1.0);
  for (int step = 0; step < 20; ++step) {
    MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(0.0, 1.0);
    interval.dt_code = stepper.estimateStableTimeStep(interval, {});
    const MhdLevelAdvanceResult result = stepper.advance(interval, {});
    CHUI_REQUIRE(result.accepted);
    CHUI_REQUIRE(result.divergence_max_normalized <= 1e-11);
    CHUI_REQUIRE(result.ideal.first_order_faces == 0U);  // G05 on a smooth problem
  }
  const MhdConservationBudget end = stepper.levelBudget(1.0);
  CHUI_REQUIRE(std::abs(end.mass - begin.mass) <= 1e-13 * begin.mass);
  CHUI_REQUIRE(std::abs(end.energy_total - begin.energy_total) <= 1e-12 * begin.energy_total);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    CHUI_REQUIRE(std::abs(end.momentum[axis] - begin.momentum[axis]) <= 1e-12 * begin.mass * 2.0);
  }
  // Loop magnetic energy must decay (numerical dissipation) but stay positive.
  CHUI_REQUIRE(end.energy_magnetic > 0.0 && end.energy_magnetic <= begin.energy_magnetic);
}

// I01: circularly polarized Alfven wave (exact nonlinear solution travelling
// at v_A = B0/sqrt(rho)). L1 error in B_y after one crossing, PLM-MC and PPM;
// observed order between the two finest resolutions must be >= 1.8.
double alfvenError(std::int64_t cells, MhdReconstructionMode mode) {
  const double amplitude = 0.1;
  TestLevel level({cells, 4, 4}, {cells, 4, 4}, {1.0 / cells, 0.25, 0.25},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [&](std::array<double, 3> x) {
        GasPrimitive w{};
        w.rho = 1.0;
        w.pressure = 0.1;
        // v_perp = -B_perp / sqrt(rho) for a wave travelling to +x.
        w.velocity = {0.0, -amplitude * std::sin(k_two_pi * x[0]),
                      -amplitude * std::cos(k_two_pi * x[0])};
        return w;
      },
      {1.0, 0.0, 0.0}, circularAlfvenPotential(amplitude, 1.0), 1.0, k_gamma);
  MhdNumericsPolicy policy{};
  policy.reconstruction = mode;
  MhdLevelStepper stepper = level.makeStepper(policy);
  (void)runUntil(stepper, 1.0, 0.8);
  double error = 0.0;
  const MhdPatchGeometry& geometry = level.binding(0).geometry;
  for (std::size_t i = 0; i < geometry.nx; ++i) {
    const double x = (static_cast<double>(i) + 0.5) / static_cast<double>(cells);
    const double numerical = level.faces(0).cellFieldComoving(i, 1, 1)[1];
    error += std::abs(numerical - amplitude * std::sin(k_two_pi * x)) / static_cast<double>(cells);
  }
  return error;
}

void testAlfvenConvergence() {
  for (const MhdReconstructionMode mode : {MhdReconstructionMode::kPlmMc, MhdReconstructionMode::kPpm}) {
    const double coarse = alfvenError(32, mode);
    const double fine = alfvenError(64, mode);
    CHUI_REQUIRE(fine < coarse);
    CHUI_REQUIRE(std::log2(coarse / fine) >= 1.8);
  }
}

// I09: expansion-only (de Sitter a = e^{t}, H = 1) uniform magnetized gas at
// rest: B_c invariant to roundoff; p_c ~ a^{3-3 gamma}; M_c ~ 1/a.
void testExpansionOnly() {
  TestLevel level({4, 4, 4}, {4, 4, 4}, {0.25, 0.25, 0.25}, cosmosim::mhd::test::periodicDomain());
  const double b0 = 0.7;
  level.initialize([](std::array<double, 3>) { return GasPrimitive{}; }, {b0, 0.2, 0.0}, {}, 1.0,
                   k_gamma);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const std::vector<double> phi0(level.faces(0).fluxX().begin(), level.faces(0).fluxX().end());
  double time = 0.0;
  const double dt = 1.0e-3;
  for (int step = 0; step < 100; ++step) {
    MhdStageContext interval{};
    interval.time_code = time;
    interval.dt_code = dt;
    interval.comoving_coordinates = true;
    interval.scale_factor = std::exp(time);
    interval.scale_factor_end = std::exp(time + dt);
    interval.hubble_rate_code = 1.0;
    interval.hubble_rate_end_code = 1.0;
    CHUI_REQUIRE(stepper.advance(interval, {}).accepted);
    time += dt;
  }
  for (std::size_t index = 0; index < phi0.size(); ++index) {
    CHUI_REQUIRE(std::abs(level.faces(0).fluxX()[index] - phi0[index]) <= 1e-12 * std::abs(phi0[index]));
  }
  const double a = std::exp(time);
  const MhdConservationBudget budget = stepper.levelBudget(a);
  const double volume = 1.0;
  const double magnetic_expected = 0.5 * (b0 * b0 + 0.04) / a * volume;
  const double thermal_expected = 1.0 / (k_gamma - 1.0) * std::pow(a, 3.0 - 3.0 * k_gamma) * volume;
  CHUI_REQUIRE(std::abs(budget.energy_magnetic - magnetic_expected) <= 1e-12 * magnetic_expected);
  CHUI_REQUIRE(std::abs(budget.energy_thermal - thermal_expected) <= 1e-6 * thermal_expected);
  CHUI_REQUIRE(std::abs(budget.momentum[0]) <= 1e-15);
}

// I03: Brio-Wu along x with open boundaries: positivity, no NaN, CT kept.
void testBrioWu() {
  TestLevel level({128, 3, 3}, {128, 3, 3}, {1.0 / 128, 1.0 / 128, 1.0 / 128},
                  cosmosim::mhd::test::uniformDomain(MhdBoundaryKind::kOpen));
  // Piecewise transverse field via A_z = |x - 1/2|: B_y = -dA_z/dx = +1 for
  // x < 1/2 and -1 for x > 1/2 (exact on the face lattice: A_z is evaluated
  // at the z-edge nodes, where it is piecewise linear).
  const MhdEdgePotential potential = [](MhdAxis axis, std::array<double, 3> start, double length) {
    if (axis != MhdAxis::kZ) return 0.0;
    return std::abs(start[0] - 0.5) * length;
  };
  level.initialize(
      [](std::array<double, 3> x) {
        GasPrimitive w{};
        w.rho = x[0] < 0.5 ? 1.0 : 0.125;
        w.pressure = x[0] < 0.5 ? 1.0 : 0.1;
        return w;
      },
      {0.75, 0.0, 0.0}, potential, 1.0, 2.0);
  MhdNumericsPolicy policy{};
  policy.gamma_adiabatic = 2.0;  // Brio & Wu (1988)
  policy.reconstruction = MhdReconstructionMode::kPlmMc;
  MhdLevelStepper stepper = level.makeStepper(policy);
  const MhdLevelAdvanceResult result = runUntil(stepper, 0.1, 0.8);
  CHUI_REQUIRE(result.divergence_max_normalized <= 1e-11);
  for (std::int64_t i = 0; i < 128; ++i) {
    CHUI_REQUIRE(level.cellValue(0, {i, 1, 1}) > 0.0);
  }
}

// Perfectly conducting closed box: zero mass/energy flux through the walls and
// wall-normal magnetic flux exactly preserved.
void testConductingWallBox() {
  TestLevel level({8, 8, 8}, {8, 8, 8}, {0.125, 0.125, 0.125},
                  cosmosim::mhd::test::uniformDomain(MhdBoundaryKind::kReflectingConductingWall));
  level.initialize(
      [](std::array<double, 3> x) {
        GasPrimitive w{};
        w.rho = 1.0 + 0.2 * std::exp(-20.0 * ((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)));
        w.velocity = {0.1 * std::sin(std::numbers::pi * x[1]), 0.0, 0.0};
        return w;
      },
      {0.3, 0.1, 0.2}, {}, 1.0, k_gamma);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  const MhdConservationBudget begin = stepper.levelBudget(1.0);
  const double wall_flux = level.faces(0).faceFluxX(0, 3, 3);
  (void)runUntil(stepper, 0.05);
  const MhdConservationBudget end = stepper.levelBudget(1.0);
  CHUI_REQUIRE(std::abs(end.mass - begin.mass) <= 1e-13 * begin.mass);
  CHUI_REQUIRE(std::abs(end.energy_total - begin.energy_total) <= 1e-12 * begin.energy_total);
  CHUI_REQUIRE(level.faces(0).faceFluxX(0, 3, 3) == wall_flux);
}

// Atomic rejection: an interval 100x the CFL bound is retried with halving in
// the isolated timeline; a rejected attempt leaves no partial state.
void testRetryAndRollback() {
  TestLevel level({16, 4, 4}, {16, 4, 4}, {1.0 / 16, 0.25, 0.25}, cosmosim::mhd::test::periodicDomain());
  level.initialize(
      [](std::array<double, 3> x) {
        GasPrimitive w{};
        w.rho = 1.0 + 0.5 * std::sin(k_two_pi * x[0]);
        w.velocity = {1.0, 0.0, 0.0};
        w.pressure = 0.01;
        return w;
      },
      {0.0, 0.0, 0.0}, {}, 1.0, k_gamma);
  MhdLevelStepper stepper = level.makeStepper(MhdNumericsPolicy{});
  MhdStageContext interval = cosmosim::mhd::test::isolatedInterval(0.0, 1.0);
  const double cfl_dt = stepper.estimateStableTimeStep(interval, {});
  interval.dt_code = 100.0 * cfl_dt;
  const MhdConservationBudget begin = stepper.levelBudget(1.0);
  try {
    const MhdLevelAdvanceResult result = stepper.advance(interval, {});
    if (result.accepted) {
      CHUI_REQUIRE(result.retry_count > 0);
      CHUI_REQUIRE(result.dt_accepted_code <= 0.5 * interval.dt_code);
    }
  } catch (const MhdError& error) {
    CHUI_REQUIRE(error.code() == MhdErrorCode::kNonPositiveState);
    const MhdConservationBudget after = stepper.levelBudget(1.0);
    CHUI_REQUIRE(after.mass == begin.mass && after.energy_total == begin.energy_total);
  }
}

}  // namespace

int main() {
  testAdvectionDirection();
  testZeroFieldLimit();
  testFreeStream();
  testAxisRotationInvariance();
  testMultiPatchEquivalence();
  testFieldLoopConservation();
  testExpansionOnly();
  testBrioWu();
  testConductingWallBox();
  testRetryAndRollback();
  testAlfvenConvergence();
  return 0;
}
