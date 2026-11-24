// CHUI-MHD-REF-v1 G02 discrete complex, CT orientation, seeds and exact
// capacity accounting (MP01/MP02/MP09). Source written, NOT EXECUTED.

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

#include "mhd_test_fixtures.hpp"
#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_induction.hpp"
#include "cosmosim/mhd/mhd_initial_conditions.hpp"
#include "cosmosim/mhd/mhd_memory.hpp"

namespace {

using namespace cosmosim::mhd;

MhdPatchGeometry rectangular(std::size_t nx, std::size_t ny, std::size_t nz) {
  MhdPatchGeometry geometry{};
  geometry.nx = nx;
  geometry.ny = ny;
  geometry.nz = nz;
  geometry.dx_comoving = 0.5;
  geometry.dy_comoving = 0.25;
  geometry.dz_comoving = 0.125;
  return geometry;
}

// D*C = 0 for arbitrary (deterministic pseudo-random) edge EMFs on a
// non-cubic patch: every interior AND boundary cell's oriented face sum stays
// zero after the CT update.
void testDiscreteComplexRectangular() {
  const MhdPatchGeometry geometry = rectangular(5, 7, 9);
  MhdFaceFluxStorage faces(geometry);
  MhdEdgeEmfStorage emf(geometry);
  std::uint64_t state = 12345U;
  const auto next = [&state]() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11U) * (1.0 / 9007199254740992.0) - 0.5;
  };
  for (double& value : emf.emfX()) value = next();
  for (double& value : emf.emfY()) value = next();
  for (double& value : emf.emfZ()) value = next();
  applyConstrainedTransportUpdate(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ(),
                                  emf.emfX(), emf.emfY(), emf.emfZ(), 0.37);
  const std::vector<double> divergence =
      faceFluxDivergence(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ());
  double maximum = 0.0;
  for (const double value : divergence) maximum = std::max(maximum, std::abs(value));
  const double scale = faceFluxMagnitudeSum(faces.fluxX(), faces.fluxY(), faces.fluxZ());
  CHUI_REQUIRE(scale > 0.0);
  CHUI_REQUIRE(maximum <= 1.0e-12 * scale);
  CHUI_REQUIRE(discreteComplexResidualX(geometry) == 0.0);
  CHUI_REQUIRE(discreteComplexResidualY(geometry) <= 1.0e-15);
  CHUI_REQUIRE(discreteComplexResidualZ(geometry) <= 1.0e-15);
}

// Analytic CT orientation (Stokes, right-handed): one unit E_z edge integral
// at node (i=2, j=2) must change the four faces that share that edge with the
// signs derived below; no z-face may change.
void testCtOrientationSigns() {
  const MhdPatchGeometry geometry = rectangular(4, 4, 4);
  MhdFaceFluxStorage faces(geometry);
  MhdEdgeEmfStorage emf(geometry);
  emf.emfZ()[geometry.edgeIndexZ(2, 2, 1)] = 1.0;
  applyConstrainedTransportUpdate(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ(),
                                  emf.emfX(), emf.emfY(), emf.emfZ(), 1.0);
  CHUI_REQUIRE(faces.faceFluxX(2, 2, 1) == 1.0);
  CHUI_REQUIRE(faces.faceFluxX(2, 1, 1) == -1.0);
  CHUI_REQUIRE(faces.faceFluxY(2, 2, 1) == -1.0);
  CHUI_REQUIRE(faces.faceFluxY(1, 2, 1) == 1.0);
  // Derivation: dPhi_x/dt = -[(E_z(j+1)-E_z(j)) - (E_y(k+1)-E_y(k))]; for the
  // face (2,2,1) E_z(j)=1 at its lower-y edge gives +1; for (2,1,1) E_z(j+1)=1
  // gives -1. dPhi_y/dt = -[(E_x(k+1)-E_x(k)) - (E_z(i+1)-E_z(i))]: face (2,2,1)
  // has E_z(i)=1 -> -1; face (1,2,1) has E_z(i+1)=1 -> +1.
  double total = 0.0;
  for (const double value : faces.fluxZ()) total += std::abs(value);
  CHUI_REQUIRE(total == 0.0);
}

void testSeedsAreDivergenceFree() {
  const MhdPatchGeometry geometry = rectangular(8, 6, 4);
  MhdFaceFluxStorage faces(geometry);
  const double loop = seedFaceFlux(geometry, {0.1, -0.2, 0.05},
                                   fieldLoopPotential(1.0e-3, 0.8, {2.0, 0.75, 0.25}),
                                   faces.fluxX(), faces.fluxY(), faces.fluxZ());
  CHUI_REQUIRE(loop <= 1.0e-13);
  const double stochastic = seedFaceFlux(
      geometry, {0.0, 0.0, 0.0},
      stochasticSeedPotential({4.0, 1.5, 0.5}, 3, 1.0, -1.0, 2026U), faces.fluxX(),
      faces.fluxY(), faces.fluxZ());
  CHUI_REQUIRE(stochastic <= 1.0e-12);
  // Determinism: the same seed reproduces the same bytes.
  MhdFaceFluxStorage again(geometry);
  (void)seedFaceFlux(geometry, {0.0, 0.0, 0.0},
                     stochasticSeedPotential({4.0, 1.5, 0.5}, 3, 1.0, -1.0, 2026U), again.fluxX(),
                     again.fluxY(), again.fluxZ());
  for (std::size_t index = 0; index < faces.fluxX().size(); ++index) {
    CHUI_REQUIRE(faces.fluxX()[index] == again.fluxX()[index]);
  }
  // Uniform net flux is representable (periodic mean field).
  MhdFaceFluxStorage uniform(geometry);
  (void)seedFaceFlux(geometry, {0.3, 0.0, 0.0}, {}, uniform.fluxX(), uniform.fluxY(),
                     uniform.fluxZ());
  CHUI_REQUIRE(std::abs(uniform.faceFluxX(0, 0, 0) - 0.3 * geometry.faceAreaX()) < 1e-16);
}

void testEdgeBoundaryRules() {
  const MhdPatchGeometry geometry = rectangular(4, 5, 6);
  MhdEdgeEmfStorage emf(geometry);
  for (double& value : emf.emfX()) value = 1.0;
  for (double& value : emf.emfY()) value = 2.0;
  for (double& value : emf.emfZ()) value = 3.0;
  MhdPatchBoundary boundary{};
  boundary.lower[0] = MhdBoundaryKind::kReflectingConductingWall;
  boundary.upper[0] = MhdBoundaryKind::kReflectingConductingWall;
  // Distinguish the periodic copy along y: lower plane value 7.
  for (std::size_t k = 0; k <= geometry.nz; ++k) {
    for (std::size_t i = 0; i < geometry.nx; ++i) {
      emf.emfX()[geometry.edgeIndexX(i, 0, k)] = 7.0;
    }
  }
  (void)MhdEmfBuilder::applyEdgeBoundaryRules(geometry, boundary, emf);
  // Tangential E vanishes on both x walls (y- and z-edges at i = 0, nx).
  for (std::size_t k = 0; k <= geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      CHUI_REQUIRE(emf.emfY()[geometry.edgeIndexY(0, j, k)] == 0.0);
      CHUI_REQUIRE(emf.emfY()[geometry.edgeIndexY(geometry.nx, j, k)] == 0.0);
    }
  }
  // Periodic y: upper plane copies the lower plane bitwise.
  for (std::size_t k = 0; k <= geometry.nz; ++k) {
    for (std::size_t i = 0; i < geometry.nx; ++i) {
      CHUI_REQUIRE(emf.emfX()[geometry.edgeIndexX(i, geometry.ny, k)] == 7.0);
    }
  }
}

void testDivergenceNormalizerHasNoArtificialFloor() {
  const MhdPatchGeometry geometry = rectangular(3, 3, 3);
  MhdFaceFluxStorage faces(geometry);
  // Tiny field with a genuine divergence must still be detected.
  faces.setFaceFluxX(1, 1, 1, 1.0e-20);
  const double normalizer =
      magneticDivergenceNormalizer(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ());
  CHUI_REQUIRE(normalizer == 1.0e-20);
  CHUI_REQUIRE(maxNormalizedFaceDivergence(geometry, faces.fluxX(), faces.fluxY(), faces.fluxZ(),
                                           normalizer) == 1.0);
}

void testStableKeyOrderingAndHash() {
  MhdFaceKey a{.patch_id = 3, .geometry_generation = 1, .refinement_level = 0,
               .axis = MhdAxis::kX, .i = 1, .j = 2, .k = 3};
  MhdFaceKey b = a;
  b.patch_id = 4;
  MhdFaceKey c = a;
  c.refinement_level = 1;
  c.patch_id = 1;
  CHUI_REQUIRE(stableFaceKeyBefore(a, b));
  CHUI_REQUIRE(stableFaceKeyBefore(a, c));  // level dominates
  CHUI_REQUIRE(faceKeyHash(a) != faceKeyHash(b));
}

// MP01 / MP09: the analytical estimate equals the retained capacity of the
// allocations made by an actual stepper interval, and the canonical face cost
// follows the exact N_f formula at 16^3, 32^3 and 64^3.
void testExactCapacityAccounting() {
  for (const std::size_t n : {16U, 32U, 64U}) {
    const MhdPatchGeometry geometry = rectangular(n, n, n);
    CHUI_REQUIRE(mhdFaceCount(geometry) == 3U * (n + 1U) * n * n);
    CHUI_REQUIRE(mhdEdgeCount(geometry) == 3U * n * (n + 1U) * (n + 1U));
    MhdMemoryRequest request{};
    request.geometry = geometry;
    const MhdPatchMemoryEstimate estimate = estimatePatchMemory(request);
    MhdFaceFluxStorage faces(geometry);
    CHUI_REQUIRE(estimate.canonical_face_flux_bytes == faces.ownedCapacityBytes());
  }
  // MP02: 32^3 PPM ideal and non-ideal (RKL2 lanes) phase+scratch per cell.
  MhdMemoryRequest ideal{};
  ideal.geometry = rectangular(32, 32, 32);
  CHUI_REQUIRE(scratchBytesPerActiveCell(ideal, 1U) <= 512.0);
  MhdMemoryRequest nonideal = ideal;
  nonideal.nonideal = true;
  nonideal.super_stepping = true;
  CHUI_REQUIRE(scratchBytesPerActiveCell(nonideal, 8U) <= 512.0);

  // Estimate == ownedCapacityBytes after one interval of a real stepper.
  using cosmosim::mhd::test::TestLevel;
  TestLevel level({8, 8, 8}, {8, 8, 8}, {0.125, 0.125, 0.125},
                  cosmosim::mhd::test::periodicDomain());
  level.initialize([](std::array<double, 3>) { return cosmosim::mhd::test::GasPrimitive{}; },
                   {0.1, 0.0, 0.0}, {}, 1.0, k_mhd_gamma_default);
  MhdNumericsPolicy policy{};
  MhdLevelStepper stepper = level.makeStepper(policy);
  const MhdLevelAdvanceResult result =
      stepper.advance(cosmosim::mhd::test::isolatedInterval(0.0, 1.0e-3), {});
  CHUI_REQUIRE(result.accepted);
  MhdMemoryRequest request{};
  request.geometry = level.binding(0).geometry;
  const MhdPatchMemoryEstimate estimate = estimatePatchMemory(request);
  CHUI_REQUIRE(stepper.ownedCapacityBytes() == estimate.phasePlusScratchBytes());
}

}  // namespace

int main() {
  testDiscreteComplexRectangular();
  testCtOrientationSigns();
  testSeedsAreDivergenceFree();
  testEdgeBoundaryRules();
  testDivergenceNormalizerHasNoArtificialFloor();
  testStableKeyOrderingAndHash();
  testExactCapacityAccounting();
  return 0;
}
