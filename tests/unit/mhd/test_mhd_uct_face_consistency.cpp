// CHUI-MHD pass 04 (M03-08, reconciliation R-UCT): the edge EMF of the ideal
// stage must reduce EXACTLY (to roundoff) to the transverse induction flux of
// the face solver actually used for the gas update, for one-dimensional
// solutions aligned with every axis, for HLLD+PPM, HLLD+PLM, HLLE+PPM and
// HLLE+PLM, on anisotropic cells and across a contact/pressure jump:
//   edge along t2 (a = d, b = t1):  E_t2 = -F_d(B_t1) * len_t2 / sqrt(a)
//   edge along t1 (a = t2, b = d):  E_t1 = +F_d(B_t2) * len_t1 / sqrt(a)
// An oblique (two-dimensional) wave checks that the corrected EMF keeps the
// discrete divergence at roundoff after the CT update.
// NOT_RUN: source-only campaign (no configure/build/test in pass 04).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>

#include "../../support/test_require.hpp"
#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_riemann.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"

namespace {

using namespace cosmosim::mhd;

constexpr double k_gamma = 5.0 / 3.0;

MhdPatchBoundary periodicBoundary() {
  MhdPatchBoundary boundary{};
  boundary.lower = {MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic,
                    MhdBoundaryKind::kPeriodic};
  boundary.upper = boundary.lower;
  return boundary;
}

struct Profile {
  std::size_t direction = 0U;  // variation axis d
  double length = 1.0;         // period along d
  bool oblique = false;        // add variation along t1 as well
};

// Edge vector potential of B_t1 = amp cos(k x_d), B_t2 = amp sin(k x_d)
// (+ optional oblique part depending on x_t1): exact edge integrals.
MhdEdgePotential profilePotential(const Profile& profile, double amplitude) {
  const double k = 2.0 * std::numbers::pi / profile.length;
  const std::size_t d = profile.direction;
  const std::size_t t1 = (d + 1U) % 3U;
  const std::size_t t2 = (d + 2U) % 3U;
  return [=](MhdAxis axis, std::array<double, 3> start, double length) {
    const std::size_t along = mhdAxisIndex(axis);
    double value = 0.0;
    if (along == t1) {  // A_t1 = -amp cos(k x_d)/k (constant along the edge)
      value += -amplitude * std::cos(k * start[d]) / k * length;
    } else if (along == t2) {  // A_t2 = -amp sin(k x_d)/k
      value += -amplitude * std::sin(k * start[d]) / k * length;
      if (profile.oblique) {  // + 0.3 amp sin(k x_t1)/k: B_d gains a t1-dependent part
        value += 0.3 * amplitude * std::sin(k * start[t1]) / k * length;
      }
    }
    return value;
  };
}

void fillGas(const MhdPatchGeometry& geometry, const MhdFaceFluxStorage& faces,
             const Profile& profile, MhdGasStateSoa& gas) {
  const std::size_t d = profile.direction;
  const std::size_t t1 = (d + 1U) % 3U;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::array<double, 3> x{(static_cast<double>(i) + 0.5) * geometry.dx_comoving,
                                      (static_cast<double>(j) + 0.5) * geometry.dy_comoving,
                                      (static_cast<double>(k) + 0.5) * geometry.dz_comoving};
        const double s = x[d] / profile.length;
        // Smooth background with a contact/pressure jump in the middle.
        double rho = (s < 0.5 ? 1.0 : 0.4) + 0.05 * std::sin(2.0 * std::numbers::pi * s);
        double pressure = (s < 0.5 ? 1.0 : 0.6);
        if (profile.oblique) {
          rho += 0.05 * std::cos(2.0 * std::numbers::pi * x[t1] / profile.length);
        }
        std::array<double, 3> velocity{};
        velocity[d] = 0.3;
        velocity[t1] = 0.1 * std::sin(2.0 * std::numbers::pi * s);
        const std::array<double, 3> b = faces.cellFieldComoving(i, j, k);
        MhdConservedState conserved{};
        conserved.rho_comoving = rho;
        conserved.momentum_density_x_comoving = rho * velocity[0];
        conserved.momentum_density_y_comoving = rho * velocity[1];
        conserved.momentum_density_z_comoving = rho * velocity[2];
        conserved.total_energy_density_comoving = pressure / (k_gamma - 1.0) +
            0.5 * rho * (velocity[0] * velocity[0] + velocity[1] * velocity[1] +
                         velocity[2] * velocity[2]) +
            0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
        gas.storeGas(geometry.cellIndex(i, j, k), conserved);
      }
    }
  }
}

std::size_t edgeIndexAlong(const MhdPatchGeometry& g, std::size_t axis, const std::array<std::size_t, 3>& e) {
  switch (axis) {
    case 0U: return g.edgeIndexX(e[0], e[1], e[2]);
    case 1U: return g.edgeIndexY(e[0], e[1], e[2]);
    default: return g.edgeIndexZ(e[0], e[1], e[2]);
  }
}

std::span<const double> emfAlong(const MhdEdgeEmfStorage& emf, std::size_t axis) {
  switch (axis) {
    case 0U: return emf.emfX();
    case 1U: return emf.emfY();
    default: return emf.emfZ();
  }
}

double spacing(const MhdPatchGeometry& g, std::size_t axis) {
  return axis == 0U ? g.dx_comoving : (axis == 1U ? g.dy_comoving : g.dz_comoving);
}

void checkAlignedReduction(std::size_t direction, MhdReconstructionMode reconstruction,
                           MhdRiemann::Solver solver) {
  MhdNumericsPolicy policy{};
  policy.reconstruction = reconstruction;
  policy.riemann = solver;
  policy.gamma_adiabatic = k_gamma;
  const MhdStepOperator op(policy);
  std::array<std::size_t, 3> n{6U, 6U, 6U};
  n[direction] = 24U;
  // Anisotropic cells (metric check): different widths on every axis.
  const std::array<double, 3> widths{0.7 / static_cast<double>(n[0]),
                                     1.3 / static_cast<double>(n[1]),
                                     0.9 / static_cast<double>(n[2])};
  MhdPatchGeometry geometry{};
  geometry.nx = n[0];
  geometry.ny = n[1];
  geometry.nz = n[2];
  geometry.dx_comoving = widths[0];
  geometry.dy_comoving = widths[1];
  geometry.dz_comoving = widths[2];
  Profile profile{};
  profile.direction = direction;
  profile.length = widths[direction] * static_cast<double>(n[direction]);
  MhdFaceFluxStorage faces(geometry);
  std::array<double, 3> uniform{};
  uniform[direction] = 1.0;  // guide field
  (void)seedFaceFlux(geometry, uniform, profilePotential(profile, 0.4), faces.fluxX(),
                     faces.fluxY(), faces.fluxZ());
  MhdGasStateSoa gas(geometry.cellCount());
  fillGas(geometry, faces, profile, gas);
  const MhdPatchBoundary boundary = periodicBoundary();
  MhdStageWorkspace workspace{};
  workspace.resize(geometry, reconstruction);
  const MhdFaceFluxView view = makeFaceFluxView(faces, 1U, 1U, 0);
  CHUI_REQUIRE(op.recoverPrimitives(geometry, gas, view, 1.0, workspace) == 0U);
  MhdStepOperator::fillPhysicalGhosts(boundary, workspace);
  MhdGasStateSoa target = gas;
  MhdStageCounters counters{};
  op.evaluateIdealStage(geometry, boundary, 1.0, 0.0, false, 1.0e-3, workspace, target, counters);
  CHUI_REQUIRE(workspace.signals.induction_correction_valid);

  const std::size_t d = direction;
  const std::size_t t1 = (d + 1U) % 3U;
  const std::size_t t2 = (d + 2U) % 3U;
  double scale = 0.0;
  double worst = 0.0;
  for (std::size_t c2 = 0; c2 < n[t2]; ++c2) {
    for (std::size_t c1 = 0; c1 < n[t1]; ++c1) {
      for (std::size_t f = 0; f <= n[d]; ++f) {
        MhdIndex3 left_cell{};
        left_cell[d] = static_cast<std::ptrdiff_t>(f) - 1;
        left_cell[t1] = static_cast<std::ptrdiff_t>(c1);
        left_cell[t2] = static_cast<std::ptrdiff_t>(c2);
        MhdIndex3 right_cell = left_cell;
        right_cell[d] = static_cast<std::ptrdiff_t>(f);
        MhdCellFaceStates left_states{};
        MhdCellFaceStates right_states{};
        CHUI_REQUIRE(op.reconstructCell(workspace.primitives, left_cell, static_cast<MhdAxis>(d),
                                        left_states, nullptr));
        CHUI_REQUIRE(op.reconstructCell(workspace.primitives, right_cell, static_cast<MhdAxis>(d),
                                        right_states, nullptr));
        const double bn = workspace.staggered.at(static_cast<MhdAxis>(d), right_cell);
        const MhdRiemann::Result solved = MhdRiemann::solve(
            static_cast<MhdAxis>(d), left_states.upper, right_states.lower, bn, k_gamma, solver);
        CHUI_REQUIRE(solved.valid);
        // Edge along t2 at (node f along d, node c1 along t1, cell c2 along t2).
        std::array<std::size_t, 3> edge_t2{};
        edge_t2[d] = f;
        edge_t2[t1] = c1;
        edge_t2[t2] = c2;
        const double expected_t2 = -solved.flux.magnetic_flux[t1] * spacing(geometry, t2);
        const double actual_t2 = emfAlong(workspace.emf, t2)[edgeIndexAlong(geometry, t2, edge_t2)];
        // Edge along t1 at (node f along d, cell c1 along t1, node c2 along t2).
        const double expected_t1 = solved.flux.magnetic_flux[t2] * spacing(geometry, t1);
        const double actual_t1 = emfAlong(workspace.emf, t1)[edgeIndexAlong(geometry, t1, edge_t2)];
        scale = std::max({scale, std::abs(expected_t1), std::abs(expected_t2)});
        worst = std::max({worst, std::abs(actual_t2 - expected_t2), std::abs(actual_t1 - expected_t1)});
      }
    }
  }
  CHUI_REQUIRE(scale > 0.0);
  CHUI_REQUIRE(worst <= 1.0e-12 * scale);
}

void checkObliqueDivergence(MhdReconstructionMode reconstruction) {
  MhdNumericsPolicy policy{};
  policy.reconstruction = reconstruction;
  policy.gamma_adiabatic = k_gamma;
  const MhdStepOperator op(policy);
  MhdPatchGeometry geometry{};
  geometry.nx = 16U;
  geometry.ny = 16U;
  geometry.nz = 6U;
  geometry.dx_comoving = 1.0 / 16.0;
  geometry.dy_comoving = 1.0 / 16.0;
  geometry.dz_comoving = 0.5 / 6.0;
  Profile profile{};
  profile.direction = 0U;
  profile.length = 1.0;
  profile.oblique = true;
  MhdFaceFluxStorage faces(geometry);
  (void)seedFaceFlux(geometry, {1.0, 0.5, 0.0}, profilePotential(profile, 0.3), faces.fluxX(),
                     faces.fluxY(), faces.fluxZ());
  MhdGasStateSoa gas(geometry.cellCount());
  fillGas(geometry, faces, profile, gas);
  const MhdPatchBoundary boundary = periodicBoundary();
  MhdStageWorkspace workspace{};
  workspace.resize(geometry, reconstruction);
  const MhdFaceFluxView view = makeFaceFluxView(faces, 1U, 1U, 0);
  CHUI_REQUIRE(op.recoverPrimitives(geometry, gas, view, 1.0, workspace) == 0U);
  MhdStepOperator::fillPhysicalGhosts(boundary, workspace);
  MhdGasStateSoa target = gas;
  MhdStageCounters counters{};
  op.evaluateIdealStage(geometry, boundary, 1.0, 0.0, false, 1.0e-3, workspace, target, counters);
  applyInductionUpdate(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ(), workspace.emf, 1.0e-3);
  const double normalizer =
      magneticDivergenceNormalizer(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ());
  CHUI_REQUIRE(maxNormalizedFaceDivergence(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ(),
                                           normalizer) <= 1.0e-13);
}

}  // namespace

int main() {
  for (std::size_t direction = 0; direction < 3U; ++direction) {
    checkAlignedReduction(direction, MhdReconstructionMode::kPpm, MhdRiemann::Solver::kHlld);
    checkAlignedReduction(direction, MhdReconstructionMode::kPlmMc, MhdRiemann::Solver::kHlld);
    checkAlignedReduction(direction, MhdReconstructionMode::kPpm, MhdRiemann::Solver::kHlle);
    checkAlignedReduction(direction, MhdReconstructionMode::kPlmMc, MhdRiemann::Solver::kHlle);
  }
  checkObliqueDivergence(MhdReconstructionMode::kPpm);
  checkObliqueDivergence(MhdReconstructionMode::kPlmMc);
  return 0;
}
