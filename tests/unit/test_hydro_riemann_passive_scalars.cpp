// Audit F3 regression: passive conserved scalars through every Riemann branch.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants:
//  * HLLE intermediate state, identical states rho=2, v_n=0.5, P=1, Z=0.2:
//    mass flux 1 and metal flux exactly 0.2 (the defective HLL branch gave 0);
//  * HLLC forced fallback (diverging flow from the audit) keeps
//    F(rho Z) = Z F(rho) for uniform Z;
//  * both face orientations (normal +x / -x) and both flow directions;
//  * HLLC contact preservation: stationary contact gives zero mass and metal
//    flux; a moving contact advects the UPWIND metallicity;
//  * near-degenerate wave speeds (zero velocity, equal states) give finite,
//    consistent fluxes.
#include <cmath>
#include <memory>

#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "../support/test_require.hpp"

namespace {

using cosmosim::hydro::HllcRiemannSolver;
using cosmosim::hydro::HlleRiemannSolver;
using cosmosim::hydro::HydroFace;
using cosmosim::hydro::HydroPrimitiveState;

constexpr double k_gamma = 1.4;

HydroFace xFace(double sign) {
  return HydroFace{.owner_cell = 0, .neighbor_cell = 1, .area_comoving = 1.0, .normal_x = sign};
}

HydroPrimitiveState state(double rho, double vx, double pressure, double z) {
  HydroPrimitiveState s;
  s.rho_comoving = rho;
  s.vel_x_peculiar = vx;
  s.pressure_comoving = pressure;
  s.metallicity_mass_fraction = z;
  return s;
}

bool near(double a, double b, double tol = 1.0e-12) {
  return std::abs(a - b) <= tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

void testAuditF3HlleUniformState() {
  const HydroPrimitiveState uniform = state(2.0, 0.5, 1.0, 0.2);
  HlleRiemannSolver hlle;
  const auto flux = hlle.computeFlux(uniform, uniform, xFace(1.0), k_gamma);
  CHUI_REQUIRE(near(flux.mass_density_comoving, 1.0));
  CHUI_REQUIRE(near(flux.metal_mass_density_comoving, 0.2));
}

void testAuditF3HllcForcedFallbackKeepsUniformMetallicity() {
  // Strong diverging flow from the audit probe: HLLC p* <= 0 -> HLLE.
  const HydroPrimitiveState left = state(2.0, -1.0, 0.1, 0.2);
  const HydroPrimitiveState right = state(2.0, 2.0, 0.1, 0.2);
  HllcRiemannSolver hllc;
  const auto flux = hllc.computeFlux(left, right, xFace(1.0), k_gamma);
  CHUI_REQUIRE(hllc.fallbackCount() == 1U);
  CHUI_REQUIRE(std::isfinite(flux.mass_density_comoving));
  CHUI_REQUIRE(near(flux.metal_mass_density_comoving, 0.2 * flux.mass_density_comoving, 1.0e-12));
}

void testBothOrientationsAndFlowDirections() {
  HlleRiemannSolver hlle;
  HllcRiemannSolver hllc;
  for (double sign : {1.0, -1.0}) {
    for (double velocity : {0.7, -0.7}) {
      const HydroPrimitiveState left = state(1.0, velocity, 1.0, 0.05);
      const HydroPrimitiveState right = state(1.3, velocity, 1.2, 0.05);
      for (const cosmosim::hydro::HydroRiemannSolver* solver :
           {static_cast<const cosmosim::hydro::HydroRiemannSolver*>(&hlle),
            static_cast<const cosmosim::hydro::HydroRiemannSolver*>(&hllc)}) {
        const auto flux = solver->computeFlux(left, right, xFace(sign), k_gamma);
        // Uniform Z: the metal flux is exactly Z times the mass flux in
        // every branch (upwind, intermediate, star states).
        CHUI_REQUIRE(near(flux.metal_mass_density_comoving, 0.05 * flux.mass_density_comoving, 1.0e-11));
      }
    }
  }
}

void testHllcContactPreservationAndUpwindMetallicity() {
  HllcRiemannSolver hllc;
  // Stationary contact: equal P, u=0, density and metallicity jump.
  const auto stationary = hllc.computeFlux(
      state(1.0, 0.0, 1.0, 0.01), state(4.0, 0.0, 1.0, 0.03), xFace(1.0), k_gamma);
  CHUI_REQUIRE(std::abs(stationary.mass_density_comoving) < 1.0e-13);
  CHUI_REQUIRE(std::abs(stationary.metal_mass_density_comoving) < 1.0e-13);
  // Moving contact (u=+0.4): HLLC is exact for a contact: F_rho = rho_L u,
  // F_Z = rho_L Z_L u.
  const auto moving = hllc.computeFlux(
      state(1.0, 0.4, 1.0, 0.01), state(4.0, 0.4, 1.0, 0.03), xFace(1.0), k_gamma);
  CHUI_REQUIRE(near(moving.mass_density_comoving, 0.4, 1.0e-12));
  CHUI_REQUIRE(near(moving.metal_mass_density_comoving, 0.004, 1.0e-12));
  // Reverse flow: upwind is the right state.
  const auto reverse = hllc.computeFlux(
      state(1.0, -0.4, 1.0, 0.01), state(4.0, -0.4, 1.0, 0.03), xFace(1.0), k_gamma);
  CHUI_REQUIRE(near(reverse.mass_density_comoving, -1.6, 1.0e-12));
  CHUI_REQUIRE(near(reverse.metal_mass_density_comoving, -0.048, 1.0e-12));
  CHUI_REQUIRE(hllc.fallbackCount() == 0U);
}

void testNearDegenerateWaveSpeeds() {
  HlleRiemannSolver hlle;
  const HydroPrimitiveState rest = state(1.0, 0.0, 1.0, 0.3);
  const auto flux = hlle.computeFlux(rest, rest, xFace(1.0), k_gamma);
  CHUI_REQUIRE(std::abs(flux.mass_density_comoving) < 1.0e-14);
  CHUI_REQUIRE(std::abs(flux.metal_mass_density_comoving) < 1.0e-14);
  CHUI_REQUIRE(near(flux.momentum_density_x_comoving, 1.0));
}

void testFactoryAndNames() {
  using cosmosim::hydro::HydroRiemannSolverKind;
  CHUI_REQUIRE(cosmosim::hydro::hydroRiemannSolverKindFromString("hllc") == HydroRiemannSolverKind::kHllc);
  CHUI_REQUIRE(cosmosim::hydro::hydroRiemannSolverKindFromString("hlle") == HydroRiemannSolverKind::kHlle);
  CHUI_REQUIRE(cosmosim::hydro::hydroRiemannSolverKindToString(HydroRiemannSolverKind::kHlle) == "hlle");
  const auto solver = cosmosim::hydro::makeHydroRiemannSolver(HydroRiemannSolverKind::kHlle);
  const auto flux = solver->computeFlux(state(2.0, 0.5, 1.0, 0.2), state(2.0, 0.5, 1.0, 0.2), xFace(1.0), k_gamma);
  CHUI_REQUIRE(near(flux.metal_mass_density_comoving, 0.2));
}

}  // namespace

int main() {
  testAuditF3HlleUniformState();
  testAuditF3HllcForcedFallbackKeepsUniformMetallicity();
  testBothOrientationsAndFlowDirections();
  testHllcContactPreservationAndUpwindMetallicity();
  testNearDegenerateWaveSpeeds();
  testFactoryAndNames();
  return 0;
}
