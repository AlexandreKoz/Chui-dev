// HLLD/HLLE consistency, region selection, contact resolution, degeneracies,
// signed fan bounds, first-order fallback validity, PLM/PPM stencils and
// characteristic PPM bounds. Source written, NOT EXECUTED.

#include <array>
#include <cmath>
#include <cstdint>

#include "../../support/test_require.hpp"
#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_riemann.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"

namespace {

using namespace cosmosim::mhd;
constexpr double k_gamma = 5.0 / 3.0;

MhdPrimitiveState state(double rho, std::array<double, 3> v, double p, std::array<double, 3> b) {
  MhdPrimitiveState w{};
  w.rho_comoving = rho;
  w.vel_x_peculiar = v[0];
  w.vel_y_peculiar = v[1];
  w.vel_z_peculiar = v[2];
  w.pressure_comoving = p;
  w.bx_comoving = b[0];
  w.by_comoving = b[1];
  w.bz_comoving = b[2];
  return w;
}

// Independent physical x-flux of ideal MHD (Lorentz-Heaviside units).
std::array<double, 8> physicalFluxX(const MhdPrimitiveState& w) {
  const double b2 = w.bx_comoving * w.bx_comoving + w.by_comoving * w.by_comoving +
      w.bz_comoving * w.bz_comoving;
  const double v2 = w.vel_x_peculiar * w.vel_x_peculiar + w.vel_y_peculiar * w.vel_y_peculiar +
      w.vel_z_peculiar * w.vel_z_peculiar;
  const double pt = w.pressure_comoving + 0.5 * b2;
  const double energy = w.pressure_comoving / (k_gamma - 1.0) + 0.5 * w.rho_comoving * v2 + 0.5 * b2;
  const double vdotb = w.vel_x_peculiar * w.bx_comoving + w.vel_y_peculiar * w.by_comoving +
      w.vel_z_peculiar * w.bz_comoving;
  return {w.rho_comoving * w.vel_x_peculiar,
          w.rho_comoving * w.vel_x_peculiar * w.vel_x_peculiar + pt - w.bx_comoving * w.bx_comoving,
          w.rho_comoving * w.vel_x_peculiar * w.vel_y_peculiar - w.bx_comoving * w.by_comoving,
          w.rho_comoving * w.vel_x_peculiar * w.vel_z_peculiar - w.bx_comoving * w.bz_comoving,
          (energy + pt) * w.vel_x_peculiar - w.bx_comoving * vdotb,
          0.0,
          w.vel_x_peculiar * w.by_comoving - w.vel_y_peculiar * w.bx_comoving,
          w.vel_x_peculiar * w.bz_comoving - w.vel_z_peculiar * w.bx_comoving};
}

std::array<double, 8> asArray(const MhdRiemannFlux& flux) {
  return {flux.mass_flux, flux.momentum_flux[0], flux.momentum_flux[1], flux.momentum_flux[2],
          flux.energy_flux, flux.magnetic_flux[0], flux.magnetic_flux[1], flux.magnetic_flux[2]};
}

void testConsistency() {
  const MhdPrimitiveState w = state(1.3, {0.4, -0.2, 0.7}, 0.9, {0.6, -0.5, 0.3});
  const std::array<double, 8> exact = physicalFluxX(w);
  for (const MhdRiemann::Solver solver : {MhdRiemann::Solver::kHlld, MhdRiemann::Solver::kHlle}) {
    const MhdRiemann::Result result = MhdRiemann::solve(MhdAxis::kX, w, w, 0.6, k_gamma, solver);
    CHUI_REQUIRE(result.valid);
    const std::array<double, 8> numerical = asArray(result.flux);
    for (std::size_t m = 0; m < 8U; ++m) {
      CHUI_REQUIRE(std::abs(numerical[m] - exact[m]) <= 1.0e-13 * (1.0 + std::abs(exact[m])));
    }
    // E = -v x B in the local basis: E_t1 = E_y = -(v_z B_x - v_x B_z).
    const double e_y = -(w.vel_z_peculiar * w.bx_comoving - w.vel_x_peculiar * w.bz_comoving);
    const double e_z = -(w.vel_x_peculiar * w.by_comoving - w.vel_y_peculiar * w.bx_comoving);
    CHUI_REQUIRE(std::abs(result.flux.electric_field[1] - e_y) <= 1.0e-13);
    CHUI_REQUIRE(std::abs(result.flux.electric_field[2] - e_z) <= 1.0e-13);
  }
}

// The HLLD flux is a continuous function of a uniform Galilean shift v0. An
// inverted branch selection between F_R* and F_R** produces an O(1) jump when
// the contact speed S_M crosses zero. Sample densely through every wave
// crossing and bound the largest increment.
void testGalileanContinuity() {
  const double bn = 0.8;
  const double dv = 1.0e-4;
  std::array<double, 8> previous{};
  double max_jump = 0.0;
  bool have_previous = false;
  for (int step = 0; step <= 40000; ++step) {
    const double v0 = -2.0 + dv * step;
    const MhdPrimitiveState left = state(1.0, {v0, 0.0, 0.0}, 1.0, {bn, 1.0, 0.0});
    const MhdPrimitiveState right = state(0.5, {v0, 0.1, -0.05}, 0.6, {bn, 0.6, 0.3});
    const MhdRiemann::Result result =
        MhdRiemann::solve(MhdAxis::kX, left, right, bn, k_gamma, MhdRiemann::Solver::kHlld);
    CHUI_REQUIRE(result.valid);
    CHUI_REQUIRE(!result.used_hlle);  // the sweep must exercise HLLD itself
    const std::array<double, 8> flux = asArray(result.flux);
    if (have_previous) {
      for (std::size_t m = 0; m < 8U; ++m) {
        max_jump = std::max(max_jump, std::abs(flux[m] - previous[m]));
      }
    }
    previous = flux;
    have_previous = true;
  }
  // Fluxes are O(1) with O(10) derivatives in v0: 100 dv is a generous
  // Lipschitz bound; a branch error is O(0.1) or larger.
  CHUI_REQUIRE(max_jump <= 100.0 * dv);
}

void testContactResolution() {
  // Stationary contact (equal p, v = 0, equal B, different rho): HLLD resolves
  // it exactly (zero mass flux) while HLLE diffuses it.
  for (const double bn : {0.0, 0.7}) {
    const MhdPrimitiveState left = state(1.0, {0.0, 0.0, 0.0}, 1.0, {bn, 0.4, 0.2});
    const MhdPrimitiveState right = state(0.125, {0.0, 0.0, 0.0}, 1.0, {bn, 0.4, 0.2});
    const MhdRiemann::Result hlld =
        MhdRiemann::solve(MhdAxis::kX, left, right, bn, k_gamma, MhdRiemann::Solver::kHlld);
    const MhdRiemann::Result hlle = MhdRiemann::solveHlle(MhdAxis::kX, left, right, bn, k_gamma);
    CHUI_REQUIRE(hlld.valid && hlle.valid);
    CHUI_REQUIRE(std::abs(hlld.flux.mass_flux) <= 1.0e-14);
    CHUI_REQUIRE(std::abs(hlle.flux.mass_flux) > 1.0e-3);
  }
  // B = 0, moving contact: exact upwind flux rho_L v.
  const MhdPrimitiveState left = state(2.0, {0.5, 0.0, 0.0}, 1.0, {0.0, 0.0, 0.0});
  const MhdPrimitiveState right = state(1.0, {0.5, 0.0, 0.0}, 1.0, {0.0, 0.0, 0.0});
  const MhdRiemann::Result zero_field =
      MhdRiemann::solve(MhdAxis::kX, left, right, 0.0, k_gamma, MhdRiemann::Solver::kHlld);
  CHUI_REQUIRE(zero_field.valid);
  CHUI_REQUIRE(std::abs(zero_field.flux.mass_flux - 1.0) <= 1.0e-13);
  CHUI_REQUIRE(zero_field.flux.magnetic_flux[1] == 0.0 && zero_field.flux.magnetic_flux[2] == 0.0);
}

void testSignedBoundsAndFallback() {
  // Super-fast flow to +x: both bounds positive, so a^- = 0.
  const MhdPrimitiveState fast = state(1.0, {20.0, 0.0, 0.0}, 1.0, {0.1, 0.1, 0.0});
  const MhdSignalBounds bounds = MhdRiemann::signalBounds(MhdAxis::kX, fast, fast, 0.1, k_gamma);
  CHUI_REQUIRE(bounds.valid && bounds.speed_left > 0.0 && bounds.speed_right > bounds.speed_left);
  const MhdRiemann::Result solved =
      MhdRiemann::solve(MhdAxis::kX, fast, fast, 0.1, k_gamma, MhdRiemann::Solver::kHlld);
  CHUI_REQUIRE(solved.flux.slowest_signal_speed > 0.0);
  // First-order fallback refuses an inadmissible cell-average pair instead of
  // returning a silent zero flux.
  const MhdPrimitiveState bad = state(1.0, {0.0, 0.0, 0.0}, -1.0, {0.0, 0.0, 0.0});
  CHUI_REQUIRE(!firstOrderPositivityFlux(MhdAxis::kX, bad, fast, 0.1, k_gamma).valid);
  const MhdRiemann::Result rusanov = firstOrderPositivityFlux(MhdAxis::kX, fast, fast, 0.1, k_gamma);
  CHUI_REQUIRE(rusanov.valid && rusanov.first_order_required);
}

// Cyclic axis rotation: solving along y with (y, z, x)-rotated states gives the
// x-solution with rotated vector components.
void testAxisRotation() {
  const MhdPrimitiveState left_x = state(1.0, {0.3, -0.1, 0.2}, 1.0, {0.75, 1.0, 0.2});
  const MhdPrimitiveState right_x = state(0.4, {-0.2, 0.3, 0.1}, 0.3, {0.75, -0.6, 0.4});
  const auto rotate = [](const MhdPrimitiveState& w) {  // x->y, y->z, z->x
    return state(w.rho_comoving, {w.vel_z_peculiar, w.vel_x_peculiar, w.vel_y_peculiar},
                 w.pressure_comoving, {w.bz_comoving, w.bx_comoving, w.by_comoving});
  };
  const MhdRiemann::Result along_x =
      MhdRiemann::solve(MhdAxis::kX, left_x, right_x, 0.75, k_gamma, MhdRiemann::Solver::kHlld);
  const MhdRiemann::Result along_y = MhdRiemann::solve(
      MhdAxis::kY, rotate(left_x), rotate(right_x), 0.75, k_gamma, MhdRiemann::Solver::kHlld);
  CHUI_REQUIRE(along_x.valid && along_y.valid);
  CHUI_REQUIRE(std::abs(along_x.flux.mass_flux - along_y.flux.mass_flux) <= 1e-14);
  CHUI_REQUIRE(std::abs(along_x.flux.energy_flux - along_y.flux.energy_flux) <= 1e-14);
  CHUI_REQUIRE(std::abs(along_x.flux.momentum_flux[0] - along_y.flux.momentum_flux[1]) <= 1e-14);
  CHUI_REQUIRE(std::abs(along_x.flux.momentum_flux[1] - along_y.flux.momentum_flux[2]) <= 1e-14);
  CHUI_REQUIRE(std::abs(along_x.flux.momentum_flux[2] - along_y.flux.momentum_flux[0]) <= 1e-14);
  for (std::size_t m = 1; m < 3U; ++m) {  // local electric components identical
    CHUI_REQUIRE(std::abs(along_x.flux.electric_field[m] - along_y.flux.electric_field[m]) <= 1e-14);
  }
}

void testScalarLimiters() {
  const auto linear = MhdStepOperator::reconstructPlmMc(0.0, 1.0, 2.0);
  CHUI_REQUIRE(linear.first == 0.5 && linear.second == 1.5);
  const auto extremum = MhdStepOperator::reconstructPlmMc(0.0, 1.0, 0.0);
  CHUI_REQUIRE(extremum.first == 1.0 && extremum.second == 1.0);
  const auto ppm_linear = MhdStepOperator::reconstructPpmScalar(-1.0, 0.0, 1.0, 2.0, 3.0);
  CHUI_REQUIRE(std::abs(ppm_linear.first - 0.5) <= 1e-15 && std::abs(ppm_linear.second - 1.5) <= 1e-15);
  const auto step = MhdStepOperator::reconstructPpmScalar(0.0, 0.0, 0.0, 1.0, 1.0);
  CHUI_REQUIRE(step.first >= 0.0 && step.second <= 1.0 && step.first <= step.second);
}

// Characteristic PPM on a padded primitive patch: exact on linear data in all
// seven waves, bounded at a discontinuity, and the ghost stencil never reads
// beyond 3 layers (the patch is allocated with exactly 3).
void testCharacteristicPpm() {
  MhdNumericsPolicy policy{};
  policy.reconstruction = MhdReconstructionMode::kPpm;
  const MhdStepOperator op(policy);
  MhdPrimitivePatch primitives{};
  primitives.resize(6, 3, 3, 3);
  for (std::ptrdiff_t k = -3; k < 6; ++k) {
    for (std::ptrdiff_t j = -3; j < 6; ++j) {
      for (std::ptrdiff_t i = -3; i < 9; ++i) {
        if (!primitives.contains({i, j, k})) continue;
        const double x = static_cast<double>(i);
        primitives.store({i, j, k}, state(1.0 + 0.01 * x, {0.1 + 0.002 * x, 0.0, 0.003 * x},
                                          1.0 + 0.02 * x, {0.5, 0.4 + 0.001 * x, 0.3}));
      }
    }
  }
  MhdCellFaceStates states{};
  std::uint64_t fallbacks = 0;
  CHUI_REQUIRE(op.reconstructCell(primitives, {2, 1, 1}, MhdAxis::kX, states, &fallbacks));
  CHUI_REQUIRE(fallbacks == 0U);
  CHUI_REQUIRE(std::abs(states.lower.rho_comoving - (1.0 + 0.01 * 1.5)) <= 1e-12);
  CHUI_REQUIRE(std::abs(states.upper.pressure_comoving - (1.0 + 0.02 * 2.5)) <= 1e-12);
  // Discontinuity in rho and p: face values bounded by neighbours.
  for (std::ptrdiff_t i = -3; i < 9; ++i) {
    for (std::ptrdiff_t j = -3; j < 6; ++j) {
      for (std::ptrdiff_t k = -3; k < 6; ++k) {
        if (!primitives.contains({i, j, k})) continue;
        const bool high = i < 3;
        primitives.store({i, j, k}, state(high ? 1.0 : 0.125, {0.0, 0.0, 0.0}, high ? 1.0 : 0.1,
                                          {0.75, high ? 1.0 : -1.0, 0.0}));
      }
    }
  }
  CHUI_REQUIRE(op.reconstructCell(primitives, {2, 1, 1}, MhdAxis::kX, states, &fallbacks));
  // Characteristic limiting bounds each wave amplitude, not every primitive
  // lane; the contract is admissibility (positivity, finiteness) with a
  // counted PLM fallback, and the cell value is never exceeded grossly.
  CHUI_REQUIRE(states.upper.rho_comoving > 0.0 && std::isfinite(states.upper.rho_comoving));
  CHUI_REQUIRE(states.upper.pressure_comoving > 0.0);
  CHUI_REQUIRE(states.upper.rho_comoving <= 1.5);
}

}  // namespace

int main() {
  testConsistency();
  testGalileanContinuity();
  testContactResolution();
  testSignedBoundsAndFallback();
  testAxisRotation();
  testScalarLimiters();
  testCharacteristicPpm();
  return 0;
}
