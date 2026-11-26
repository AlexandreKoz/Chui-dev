#pragma once

// CHUI-RT-BH-AGN-REF-v1, AGN companion section 6 (P1-AGN01): CT-compatible
// magnetic injection. STATUS: P1 source implementation, NOT qualified and NOT
// wired (physics.agn_magnetic_injection is refused at admission).
//
// Geometry: the repository's staggered face-integrated comoving flux layout
// (mhd_state.hpp: faceIndexX = i + (nx+1)(j + ny k), etc.; edge layouts
// edgeIndexX/Y/Z). The injection is an EDGE vector potential increment
// A_e = n_e f(r_e) (n = jet axis, f compact (1 - (r/h)^2)^3), whose discrete
// Stokes circulation gives face-flux increments
//   dPhi_x = dz [A_z(j+1) - A_z(j)] - dy [A_y(k+1) - A_y(k)]   (cyclic),
// so the discrete divergence of every cell is unchanged to round-off (div of
// a discrete curl vanishes identically) and the compact support carries no net
// boundary flux or monopole. The field is toroidal about n (regular on axis).
// Energy (Lorentz-Heaviside, e_B = B^2/2, cell-centred average of the two
// face fluxes per axis divided by the face area):
//   dE_B(s) = s sum_i V_i B_i . b_i + s^2/2 sum_i V_i |b_i|^2,
// including the cross term with the pre-existing field; the amplitude s >= 0
// solves dE_B(s) = Q with the stable quadratic root (no feasible root ->
// reject/defer). The CT update must still be applied by the MHD authority at
// a synchronized stage (with AMR reflux-curl); that wiring is the open gap.

#include <array>
#include <cstdint>
#include <vector>

namespace cosmosim::physics {

struct MagneticPatchGeometry {
  std::uint32_t nx = 0U;
  std::uint32_t ny = 0U;
  std::uint32_t nz = 0U;
  std::array<double, 3> origin{};
  std::array<double, 3> width{};
};

struct MagneticFaceFluxes {
  std::vector<double> x;  // (nx+1) ny nz
  std::vector<double> y;  // nx (ny+1) nz
  std::vector<double> z;  // nx ny (nz+1)
};

[[nodiscard]] MagneticFaceFluxes buildToroidalMagneticInjection(const MagneticPatchGeometry& geometry,
                                                               const std::array<double, 3>& centre,
                                                               const std::array<double, 3>& axis, double radius);
// Max |div Phi| over cells (sum of signed face fluxes per cell).
[[nodiscard]] double magneticDiscreteDivergenceMax(const MagneticPatchGeometry& geometry,
                                                   const MagneticFaceFluxes& faces);
struct MagneticEnergyQuadratic {
  double linear = 0.0;     // sum V B . b
  double quadratic = 0.0;  // 0.5 sum V |b|^2
};
[[nodiscard]] MagneticEnergyQuadratic magneticInjectionEnergyTerms(const MagneticPatchGeometry& geometry,
                                                                   const MagneticFaceFluxes& existing,
                                                                   const MagneticFaceFluxes& unit_injection);
// Amplitude s >= 0 with dE_B(s) = Q; false when infeasible.
[[nodiscard]] bool solveMagneticInjectionAmplitude(const MagneticEnergyQuadratic& terms, double energy,
                                                   double& amplitude) noexcept;

}  // namespace cosmosim::physics
