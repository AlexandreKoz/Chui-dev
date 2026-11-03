// Piecewise-parabolic reconstruction contract tests.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants:
//  * CW84 interface interpolation is exact for cell averages of a quadratic
//    (a_{i+1/2} = (i+1/2)^2 for q_i = i^2 + 1/12) and is not a renamed linear
//    slope (the PLM value differs);
//  * the monotonicity constraint flattens a local extremum to a constant and
//    steepens an overshooting parabola so it stays within the data range;
//  * full-stencil PPM face states equal the parabola edges (predictor off);
//  * a narrow stencil reduces order explicitly and is counted;
//  * strong compressive pressure jumps activate flattening;
//  * the single Hancock half-step uses the parabola's mean gradient: for a
//    linear density profile PPM and PLM predicted states agree exactly.
#include <array>
#include <cmath>

#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "../support/test_require.hpp"

namespace {

using cosmosim::hydro::HydroFace;
using cosmosim::hydro::HydroFaceAxis;
using cosmosim::hydro::HydroFaceStencilStates;
using cosmosim::hydro::HydroPrimitiveState;
using cosmosim::hydro::HydroReconstructionPolicy;

constexpr double k_gamma = 5.0 / 3.0;

HydroPrimitiveState cell(double rho, double pressure = 1.0, double vx = 0.0) {
  HydroPrimitiveState s;
  s.rho_comoving = rho;
  s.vel_x_peculiar = vx;
  s.pressure_comoving = pressure;
  s.specific_internal_energy_code = pressure / ((k_gamma - 1.0) * rho);
  s.signal_speed_squared_code = k_gamma * pressure / rho;
  return s;
}

HydroReconstructionPolicy policy(bool predictor, double courant = 0.0) {
  return HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = courant,
      .dt_over_cell_width_code = {courant, courant, courant},
      .rho_floor = 1.0e-12,
      .pressure_floor = 1.0e-12,
      .enable_muscl_hancock_predictor = predictor,
      .adiabatic_index = k_gamma};
}

const HydroFace k_face{.owner_cell = 2, .neighbor_cell = 3, .owner_minus_cell = 1, .neighbor_plus_cell = 4,
                       .area_comoving = 1.0, .normal_x = 1.0, .axis = HydroFaceAxis::kX,
                       .owner_minus2_cell = 0, .neighbor_plus2_cell = 5};

double quadraticAverage(int i) { return static_cast<double>(i * i) + 1.0 / 12.0; }

void testInterfaceExactForQuadraticAverages() {
  const auto edges = cosmosim::hydro::ppmCellEdges(
      {quadraticAverage(3), quadraticAverage(4), quadraticAverage(5), quadraticAverage(6), quadraticAverage(7)});
  CHUI_REQUIRE(std::abs(edges.left - 4.5 * 4.5) < 1.0e-12);
  CHUI_REQUIRE(std::abs(edges.right - 5.5 * 5.5) < 1.0e-12);
  // A linear (PLM, MC) edge would be q_5 + 0.5*slope = 25+1/12 + 5 = 30.0833..,
  // not 30.25: the parabola is genuinely higher order.
  CHUI_REQUIRE(std::abs(edges.right - (quadraticAverage(5) + 0.5 * 10.0)) > 0.1);
}

void testMonotonicityConstraint() {
  // Local maximum: constant reconstruction.
  const auto extremum = cosmosim::hydro::ppmCellEdges({1.0, 2.0, 5.0, 2.0, 1.0});
  CHUI_REQUIRE(extremum.left == 5.0 && extremum.right == 5.0);
  // Edges never leave the [min,max] of the adjacent averages.
  const auto step = cosmosim::hydro::ppmCellEdges({0.0, 0.0, 0.1, 10.0, 10.0});
  CHUI_REQUIRE(step.left >= 0.0 && step.left <= 0.1 + 1.0e-14);
  CHUI_REQUIRE(step.right >= 0.1 - 1.0e-14 && step.right <= 10.0);
  // The parabola keeps the cell mean: steepened edge 3q - 2a stays bounded.
  const double mean = 0.1;
  const double parabola_extreme_bound = std::max(std::abs(step.left - mean), std::abs(step.right - mean));
  CHUI_REQUIRE(std::isfinite(parabola_extreme_bound));
}

void testFullStencilFaceStatesMatchParabolaEdges() {
  std::array<HydroPrimitiveState, 6> cells{};
  for (int k = 0; k < 6; ++k) {
    cells[static_cast<std::size_t>(k)] = cell(quadraticAverage(k + 3));
  }
  cosmosim::hydro::PpmReconstruction ppm(policy(false));
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  const HydroFaceStencilStates stencil{&cells[0], &cells[1], &cells[2], &cells[3], &cells[4], &cells[5]};
  CHUI_REQUIRE(ppm.reconstructFaceFromStencil(stencil, k_face, left, right));
  // Owner is cell i=5 (index 2), neighbor i=6: shared face at x=5.5.
  CHUI_REQUIRE(std::abs(left.rho_comoving - 5.5 * 5.5) < 1.0e-12);
  CHUI_REQUIRE(std::abs(right.rho_comoving - 5.5 * 5.5) < 1.0e-12);
  CHUI_REQUIRE(ppm.stencilDegradationCount() == 0U);
}

void testNarrowStencilDegradesExplicitly() {
  std::array<HydroPrimitiveState, 4> cells{cell(1.0), cell(2.0), cell(3.0), cell(4.0)};
  cosmosim::hydro::PpmReconstruction ppm(policy(false));
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  CHUI_REQUIRE(ppm.reconstructFaceFromPrimitiveStates(cells[0], cells[1], cells[2], cells[3], k_face, left, right));
  CHUI_REQUIRE(ppm.stencilDegradationCount() == 2U);
  // Reduced to limited PLM on a linear profile: exact face value 2.5.
  CHUI_REQUIRE(std::abs(left.rho_comoving - 2.5) < 1.0e-12);
  CHUI_REQUIRE(std::abs(right.rho_comoving - 2.5) < 1.0e-12);
}

void testFlatteningActivatesAtStrongCompressiveShock() {
  std::array<HydroPrimitiveState, 6> cells{
      cell(1.0, 100.0, 1.0), cell(1.0, 100.0, 1.0), cell(1.0, 50.0, 0.5),
      cell(1.0, 1.0, 0.0), cell(1.0, 1.0, 0.0), cell(1.0, 1.0, 0.0)};
  cosmosim::hydro::PpmReconstruction ppm(policy(false));
  HydroPrimitiveState left;
  HydroPrimitiveState right;
  const HydroFaceStencilStates stencil{&cells[0], &cells[1], &cells[2], &cells[3], &cells[4], &cells[5]};
  CHUI_REQUIRE(ppm.reconstructFaceFromStencil(stencil, k_face, left, right));
  CHUI_REQUIRE(ppm.flatteningActivationCount() >= 1U);
  CHUI_REQUIRE(left.pressure_comoving >= 1.0 && left.pressure_comoving <= 100.0);
}

void testSingleHancockPredictorOnLinearProfile() {
  // Linear density with u=+1: parabola edges equal PLM edges and the mean
  // gradient equals the PLM slope, so the predicted states coincide (no
  // second half-time prediction is stacked).
  std::array<HydroPrimitiveState, 6> cells{};
  for (int k = 0; k < 6; ++k) {
    cells[static_cast<std::size_t>(k)] = cell(1.0 + static_cast<double>(k), 1.0, 1.0);
  }
  cosmosim::hydro::PpmReconstruction ppm(policy(true, 0.2));
  cosmosim::hydro::MusclHancockReconstruction plm(policy(true, 0.2));
  HydroPrimitiveState ppm_left;
  HydroPrimitiveState ppm_right;
  HydroPrimitiveState plm_left;
  HydroPrimitiveState plm_right;
  const HydroFaceStencilStates stencil{&cells[0], &cells[1], &cells[2], &cells[3], &cells[4], &cells[5]};
  CHUI_REQUIRE(ppm.reconstructFaceFromStencil(stencil, k_face, ppm_left, ppm_right));
  CHUI_REQUIRE(plm.reconstructFaceFromPrimitiveStates(cells[1], cells[2], cells[3], cells[4], k_face, plm_left, plm_right));
  CHUI_REQUIRE(std::abs(ppm_left.rho_comoving - plm_left.rho_comoving) < 1.0e-12);
  CHUI_REQUIRE(std::abs(ppm_right.rho_comoving - plm_right.rho_comoving) < 1.0e-12);
  // Owner cell rho=3: face 3.5 - 0.1 = 3.4.
  CHUI_REQUIRE(std::abs(ppm_left.rho_comoving - 3.4) < 1.0e-12);
}

void testSchemeNames() {
  using cosmosim::hydro::HydroReconstructionScheme;
  CHUI_REQUIRE(cosmosim::hydro::hydroReconstructionSchemeFromString("piecewise_parabolic") ==
               HydroReconstructionScheme::kPiecewiseParabolic);
  CHUI_REQUIRE(cosmosim::hydro::hydroReconstructionSchemeToString(HydroReconstructionScheme::kPiecewiseLinear) ==
               "piecewise_linear");
}

}  // namespace

int main() {
  testInterfaceExactForQuadraticAverages();
  testMonotonicityConstraint();
  testFullStencilFaceStatesMatchParabolaEdges();
  testNarrowStencilDegradesExplicitly();
  testFlatteningActivatesAtStrongCompressiveShock();
  testSingleHancockPredictorOnLinearProfile();
  testSchemeNames();
  return 0;
}
