// Audit F1 regression and Hancock-predictor contract tests.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants exercised (all with analytic references):
//  * signed forward half-step advection: ramp [1,2,3,4], u_n=+1, dt/(a dx)=0.2
//    gives the owner face density 2.5 - 0.1 = 2.4 (the defective code gave 2.6);
//  * reversed velocity mirrors the sign (owner face 2.6 for u_n=-1);
//  * identical results along x, y and z with matching per-axis Courant factors;
//  * the same physical setup expressed with a reversed face normal (stencil
//    ordered along -x) yields the mirrored interface states;
//  * the comoving factor: dt/(a dx) at a=2 equals dt'/dx with dt'=dt/2;
//  * pressure-gradient acceleration has the physical sign du/dt = -dP/dx / rho;
//  * passive metallicity is advected by -u dZ/dx only;
//  * the limited slope never creates new extrema (monotone face values).
#include <array>
#include <cmath>
#include <cstddef>

#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "../support/test_require.hpp"

namespace {

using cosmosim::hydro::HydroFace;
using cosmosim::hydro::HydroFaceAxis;
using cosmosim::hydro::HydroPrimitiveState;
using cosmosim::hydro::HydroReconstructionPolicy;
using cosmosim::hydro::MusclHancockReconstruction;

constexpr double k_gamma = 5.0 / 3.0;
constexpr double k_tol = 1.0e-12;

HydroPrimitiveState rampState(double rho, double u_axis, HydroFaceAxis axis, double pressure = 1.0) {
  HydroPrimitiveState state;
  state.rho_comoving = rho;
  state.pressure_comoving = pressure;
  state.specific_internal_energy_code = pressure / ((k_gamma - 1.0) * rho);
  state.signal_speed_squared_code = k_gamma * pressure / rho;
  switch (axis) {
    case HydroFaceAxis::kX: state.vel_x_peculiar = u_axis; break;
    case HydroFaceAxis::kY: state.vel_y_peculiar = u_axis; break;
    case HydroFaceAxis::kZ: state.vel_z_peculiar = u_axis; break;
  }
  return state;
}

HydroFace axisFace(HydroFaceAxis axis, double sign = 1.0) {
  HydroFace face{.owner_cell = 1, .neighbor_cell = 2, .owner_minus_cell = 0, .neighbor_plus_cell = 3,
                 .area_comoving = 1.0, .axis = axis};
  switch (axis) {
    case HydroFaceAxis::kX: face.normal_x = sign; break;
    case HydroFaceAxis::kY: face.normal_y = sign; break;
    case HydroFaceAxis::kZ: face.normal_z = sign; break;
  }
  return face;
}

MusclHancockReconstruction makeReconstruction(std::array<double, 3> courant) {
  return MusclHancockReconstruction(HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = 0.0,
      .dt_over_cell_width_code = courant,
      .rho_floor = 1.0e-12,
      .pressure_floor = 1.0e-12,
      .enable_muscl_hancock_predictor = true,
      .adiabatic_index = k_gamma});
}

void testAuditF1SignedRampOnEveryAxis() {
  for (HydroFaceAxis axis : {HydroFaceAxis::kX, HydroFaceAxis::kY, HydroFaceAxis::kZ}) {
    std::array<double, 3> courant{0.0, 0.0, 0.0};
    courant[axis == HydroFaceAxis::kX ? 0 : (axis == HydroFaceAxis::kY ? 1 : 2)] = 0.2;
    const auto reconstruction = makeReconstruction(courant);
    const auto q0 = rampState(1.0, 1.0, axis);
    const auto q1 = rampState(2.0, 1.0, axis);
    const auto q2 = rampState(3.0, 1.0, axis);
    const auto q3 = rampState(4.0, 1.0, axis);
    HydroPrimitiveState left;
    HydroPrimitiveState right;
    CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(q0, q1, q2, q3, axisFace(axis), left, right));
    // Spatial face value 2.5; forward half step: -(dt/2)(u d rho/dx) = -0.1.
    CHUI_REQUIRE(std::abs(left.rho_comoving - 2.4) < k_tol);
    // Neighbor: 3 - 0.5 - 0.1 = 2.4 as well (constant u, linear profile).
    CHUI_REQUIRE(std::abs(right.rho_comoving - 2.4) < k_tol);
    // Uniform velocity and pressure are preserved exactly.
    CHUI_REQUIRE(std::abs(left.pressure_comoving - 1.0) < k_tol);
    CHUI_REQUIRE(std::abs(right.pressure_comoving - 1.0) < k_tol);
  }
}

void testReversedVelocityMirrorsTheSign() {
  const auto reconstruction = makeReconstruction({0.2, 0.0, 0.0});
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(
      rampState(1.0, -1.0, HydroFaceAxis::kX), rampState(2.0, -1.0, HydroFaceAxis::kX),
      rampState(3.0, -1.0, HydroFaceAxis::kX), rampState(4.0, -1.0, HydroFaceAxis::kX),
      axisFace(HydroFaceAxis::kX), left, right));
  CHUI_REQUIRE(std::abs(left.rho_comoving - 2.6) < k_tol);
  CHUI_REQUIRE(std::abs(right.rho_comoving - 2.6) < k_tol);
}

void testReversedNormalOrientationIsConsistent() {
  // Physical setup: rho increases with x, u_x=+1. Describe the same face
  // with normal -x: owner is the cell at larger x and the stencil runs along -x.
  const auto reconstruction = makeReconstruction({0.2, 0.0, 0.0});
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(
      rampState(4.0, 1.0, HydroFaceAxis::kX), rampState(3.0, 1.0, HydroFaceAxis::kX),
      rampState(2.0, 1.0, HydroFaceAxis::kX), rampState(1.0, 1.0, HydroFaceAxis::kX),
      axisFace(HydroFaceAxis::kX, -1.0), left, right));
  // The physical interface state is unique: 2.4 on both sides.
  CHUI_REQUIRE(std::abs(left.rho_comoving - 2.4) < k_tol);
  CHUI_REQUIRE(std::abs(right.rho_comoving - 2.4) < k_tol);
}

void testComovingScaleFactorEnteringTheCourantFactor() {
  // dt=0.4, a=2, dx=1 must equal dt=0.2, a=1, dx=1.
  const auto comoving = cosmosim::hydro::hydroTransportCourantFactors(0.4, 2.0, 1.0, 1.0, 1.0);
  const auto static_frame = cosmosim::hydro::hydroTransportCourantFactors(0.2, 1.0, 1.0, 1.0, 1.0);
  CHUI_REQUIRE(std::abs(comoving[0] - static_frame[0]) < k_tol);
  CHUI_REQUIRE(std::abs(comoving[0] - 0.2) < k_tol);
  // Anisotropic widths: lambda_axis = dt/(a dx_axis).
  const auto anisotropic = cosmosim::hydro::hydroTransportCourantFactors(0.4, 2.0, 0.5, 1.0, 2.0);
  CHUI_REQUIRE(std::abs(anisotropic[0] - 0.4) < k_tol);
  CHUI_REQUIRE(std::abs(anisotropic[1] - 0.2) < k_tol);
  CHUI_REQUIRE(std::abs(anisotropic[2] - 0.1) < k_tol);
  const auto reconstruction = makeReconstruction(comoving);
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(
      rampState(1.0, 1.0, HydroFaceAxis::kX), rampState(2.0, 1.0, HydroFaceAxis::kX),
      rampState(3.0, 1.0, HydroFaceAxis::kX), rampState(4.0, 1.0, HydroFaceAxis::kX),
      axisFace(HydroFaceAxis::kX), left, right));
  CHUI_REQUIRE(std::abs(left.rho_comoving - 2.4) < k_tol);
}

void testPressureGradientAccelerationSign() {
  // Static gas, uniform density 1, pressure ramp 1,2,3,4 along +x.
  // du/dt = -(1/rho) dP/dx = -1 per cell; owner face u = 0 + 0.5*0.2*(-1) = -0.1.
  const auto reconstruction = makeReconstruction({0.2, 0.0, 0.0});
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(
      rampState(1.0, 0.0, HydroFaceAxis::kX, 1.0), rampState(1.0, 0.0, HydroFaceAxis::kX, 2.0),
      rampState(1.0, 0.0, HydroFaceAxis::kX, 3.0), rampState(1.0, 0.0, HydroFaceAxis::kX, 4.0),
      axisFace(HydroFaceAxis::kX), left, right));
  CHUI_REQUIRE(std::abs(left.vel_x_peculiar + 0.1) < k_tol);
  CHUI_REQUIRE(std::abs(right.vel_x_peculiar + 0.1) < k_tol);
  // Transverse velocity is not accelerated by a normal pressure gradient.
  CHUI_REQUIRE(std::abs(left.vel_y_peculiar) < k_tol);
  // Static gas: pressure itself is not advected (u=0, du/dx=0).
  CHUI_REQUIRE(std::abs(left.pressure_comoving - 2.5) < k_tol);
}

void testPassiveMetallicityAdvection() {
  auto make = [](double z) {
    auto s = rampState(1.0, 1.0, HydroFaceAxis::kX);
    s.metallicity_mass_fraction = z;
    return s;
  };
  const auto reconstruction = makeReconstruction({0.2, 0.0, 0.0});
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(
      make(0.01), make(0.02), make(0.03), make(0.04), axisFace(HydroFaceAxis::kX), left, right));
  // 0.025 - 0.5*0.2*1*0.01 = 0.024.
  CHUI_REQUIRE(std::abs(left.metallicity_mass_fraction - 0.024) < k_tol);
}

void testLimitedFaceValuesStayWithinNeighbourBounds() {
  // Step profile: the limited, unpredicted face value must lie between the
  // two adjacent cell averages for every limiter (no new extrema).
  for (auto limiter : {cosmosim::hydro::HydroSlopeLimiter::kMinmod,
                       cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
                       cosmosim::hydro::HydroSlopeLimiter::kVanLeer}) {
    MusclHancockReconstruction reconstruction(HydroReconstructionPolicy{
        .limiter = limiter, .rho_floor = 1.0e-12, .pressure_floor = 1.0e-12,
        .enable_muscl_hancock_predictor = false, .adiabatic_index = k_gamma});
    HydroPrimitiveState left;
    HydroPrimitiveState right;
    CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(
        rampState(1.0, 0.0, HydroFaceAxis::kX), rampState(1.0, 0.0, HydroFaceAxis::kX),
        rampState(10.0, 0.0, HydroFaceAxis::kX), rampState(10.0, 0.0, HydroFaceAxis::kX),
        axisFace(HydroFaceAxis::kX), left, right));
    CHUI_REQUIRE(left.rho_comoving >= 1.0 - k_tol && left.rho_comoving <= 10.0 + k_tol);
    CHUI_REQUIRE(right.rho_comoving >= 1.0 - k_tol && right.rho_comoving <= 10.0 + k_tol);
  }
}

void testIdealGasFaceStateIsThermodynamicallyClosed() {
  const auto reconstruction = makeReconstruction({0.2, 0.0, 0.0});
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  CHUI_REQUIRE(reconstruction.reconstructFaceFromPrimitiveStates(
      rampState(1.0, 0.3, HydroFaceAxis::kX, 1.0), rampState(1.2, 0.4, HydroFaceAxis::kX, 1.3),
      rampState(1.4, 0.5, HydroFaceAxis::kX, 1.6), rampState(1.6, 0.6, HydroFaceAxis::kX, 1.9),
      axisFace(HydroFaceAxis::kX), left, right));
  for (const auto& state : {left, right}) {
    CHUI_REQUIRE(std::abs(state.specific_internal_energy_code -
                          state.pressure_comoving / ((k_gamma - 1.0) * state.rho_comoving)) < 1.0e-12);
    CHUI_REQUIRE(std::abs(state.signal_speed_squared_code -
                          k_gamma * state.pressure_comoving / state.rho_comoving) < 1.0e-12);
  }
}

}  // namespace

int main() {
  testAuditF1SignedRampOnEveryAxis();
  testReversedVelocityMirrorsTheSign();
  testReversedNormalOrientationIsConsistent();
  testComovingScaleFactorEnteringTheCourantFactor();
  testPressureGradientAccelerationSign();
  testPassiveMetallicityAdvection();
  testLimitedFaceValuesStayWithinNeighbourBounds();
  testIdealGasFaceStateIsThermodynamicallyClosed();
  return 0;
}
