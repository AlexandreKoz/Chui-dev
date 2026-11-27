#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_riemann.hpp"

namespace cosmosim::mhd {

// Integer lattice coordinate triple (x, y, z). Ghost coordinates are negative
// or >= extent. Used by every padded MHD stage buffer so that axis handling is
// written once and permuted cyclically (normal a, tangential a+1, a+2).
using MhdIndex3 = std::array<std::ptrdiff_t, 3>;

[[nodiscard]] constexpr std::size_t mhdAxisIndex(MhdAxis axis) noexcept {
  return static_cast<std::size_t>(axis);
}
[[nodiscard]] constexpr MhdAxis mhdNextAxis(MhdAxis axis, int shift = 1) noexcept {
  return static_cast<MhdAxis>((static_cast<int>(axis) + shift) % 3);
}

// Monotonized-central limited slope (van Leer 1977): returns the limited
// centered difference used by PLM-MC, the corner extrapolation and the
// transverse staggered-field reconstruction.
[[nodiscard]] inline double mcLimitedSlope(double minus, double center, double plus) noexcept {
  const double left = center - minus;
  const double right = plus - center;
  if (left * right <= 0.0) {
    return 0.0;
  }
  const double central = 0.5 * (plus - minus);
  const double bound = 2.0 * (std::abs(left) < std::abs(right) ? std::abs(left)
                                                               : std::abs(right));
  const double magnitude = std::abs(central) < bound ? std::abs(central) : bound;
  return central > 0.0 ? magnitude : -magnitude;
}

// Padded cell-centred primitive cache for one patch and one stage. Ghost
// layers are filled from physical boundary kinds or (kNeighborPatch) by the
// level stepper. Phase-local scratch, never an authority.
//
// Magnetic lanes hold B_eff = B_c / sqrt(a) (the Riemann-normalized field),
// pressure holds p_c = a^3 p_phys, rho holds rho_c = a^3 rho_phys, velocity
// is physical peculiar velocity; metal holds the metal mass fraction Z.
struct MhdPrimitivePatch {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  std::size_t ghost = 0;
  std::size_t px = 0;
  std::size_t py = 0;
  std::size_t pz = 0;

  MhdDoubleArray rho;
  MhdDoubleArray vx;
  MhdDoubleArray vy;
  MhdDoubleArray vz;
  MhdDoubleArray pressure;
  MhdDoubleArray bx;
  MhdDoubleArray by;
  MhdDoubleArray bz;
  MhdDoubleArray metal;

  void resize(std::size_t patch_nx, std::size_t patch_ny, std::size_t patch_nz,
              std::size_t ghost_cells);
  [[nodiscard]] std::size_t index(std::ptrdiff_t i, std::ptrdiff_t j,
                                  std::ptrdiff_t k) const noexcept {
    const auto pi = static_cast<std::size_t>(i + static_cast<std::ptrdiff_t>(ghost));
    const auto pj = static_cast<std::size_t>(j + static_cast<std::ptrdiff_t>(ghost));
    const auto pk = static_cast<std::size_t>(k + static_cast<std::ptrdiff_t>(ghost));
    return pi + px * (pj + py * pk);
  }
  [[nodiscard]] std::size_t index(const MhdIndex3& c) const noexcept {
    return index(c[0], c[1], c[2]);
  }
  [[nodiscard]] bool contains(const MhdIndex3& c) const noexcept;
  [[nodiscard]] const MhdDoubleArray& velocity(MhdAxis axis) const noexcept;
  [[nodiscard]] const MhdDoubleArray& field(MhdAxis axis) const noexcept;
  [[nodiscard]] MhdPrimitiveState load(const MhdIndex3& c) const noexcept;
  void store(const MhdIndex3& c, const MhdPrimitiveState& state) noexcept;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Padded staggered normal field B_eff on the three face families. A face of
// normal axis a has NO normal ghosts (index 0..n_a) and `tangential_ghost`
// ghost rows in each tangential direction. Filled from the authoritative face
// flux: B_eff = Phi / (A_face * sqrt(a)).
struct MhdStaggeredFieldPatch {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  std::size_t tangential_ghost = 0;
  std::array<MhdDoubleArray, 3> normal_field{};

  void resize(std::size_t patch_nx, std::size_t patch_ny, std::size_t patch_nz,
              std::size_t ghost);
  [[nodiscard]] std::array<std::size_t, 3> extents(MhdAxis normal) const noexcept;
  [[nodiscard]] bool contains(MhdAxis normal, const MhdIndex3& c) const noexcept;
  [[nodiscard]] std::size_t index(MhdAxis normal, const MhdIndex3& c) const noexcept;
  [[nodiscard]] double at(MhdAxis normal, const MhdIndex3& c) const noexcept {
    return normal_field[mhdAxisIndex(normal)][index(normal, c)];
  }
  [[nodiscard]] double& at(MhdAxis normal, const MhdIndex3& c) noexcept {
    return normal_field[mhdAxisIndex(normal)][index(normal, c)];
  }
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Per-face signed HLL fan data on a ghost-extended face lattice: faces of
// normal a span 0..n_a along a and -1..n_t along each tangential direction, so
// every edge of the patch (including boundary edges) has its four adjacent
// faces available without clamping. Values are the non-negative upwind
// weights a^+ = max(0, S_R), a^- = max(0, -S_L) of the face fan.
struct MhdFaceSignalPatch {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  std::array<MhdDoubleArray, 3> speed_plus{};
  std::array<MhdDoubleArray, 3> speed_minus{};
  // Pass 04 (M03-08, reconciliation R-UCT): face-consistency correction of
  // the transverse induction flux, same index space as the speeds.
  //   induction_correction[n][s] = F_face(B_t) - F_HLL^PLM(B_t),  t = (n+1+s)%3,
  // where F_face is the flux of the ACTUAL face solver used for the gas update
  // (HLLD, its HLLE fallback or the first-order positivity flux, on the
  // primary PLM/PPM states) and F_HLL^PLM is the HLL flux on exactly the
  // PLM-MC states and signed fan bounds that UCT-HLL uses at that face in its
  // one-dimensional reduction. Valid only when `induction_correction_valid`.
  std::array<std::array<MhdDoubleArray, 2>, 3> induction_correction{};
  bool induction_correction_valid = false;

  void resize(const MhdPatchGeometry& patch_geometry);
  [[nodiscard]] std::array<std::size_t, 3> extents(MhdAxis normal) const noexcept;
  [[nodiscard]] bool contains(MhdAxis normal, const MhdIndex3& c) const noexcept;
  [[nodiscard]] std::size_t index(MhdAxis normal, const MhdIndex3& c) const noexcept;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
};

// Fill the tangential ghost rows of the staggered normal field for one patch
// from its interior faces using the physical boundary kinds:
//   periodic: exact cyclic wrap (a periodic patch must carry identical
//             boundary-face values on both ends; this is validated);
//   open:     zero-gradient copy of the nearest interior row;
//   wall:     field components tangential to a conducting wall are mirrored
//             with a sign change (odd), consistent with fillGhostPrimitives.
// Sides of kind kNeighborPatch are left untouched for the level stepper.
void fillStaggeredGhosts(const MhdPatchBoundary& boundary, MhdStaggeredFieldPatch& staggered);

// Upwind constrained-transport EMF builder (CHUI-MHD-REF v1.0 sections 3,
// 6.1, 7). The ideal edge EMF is the genuine two-dimensional HLL (UCT-HLL)
// electric field of Londrillo & Del Zanna (2004) / Del Zanna et al. (2007),
// Mignone & Del Zanna (2021) eq. form, for an edge along axis c with
// tangential axes (a, b) = (c+1, c+2) cyclic and E_c = -(v_a B_b - v_b B_a):
//
//   E_c = [a_a+ a_b+ E(--) + a_a+ a_b- E(-+) + a_a- a_b+ E(+-) + a_a- a_b- E(++)]
//         / [(a_a+ + a_a-)(a_b+ + a_b-)]
//       + a_a+ a_a- / (a_a+ + a_a-) * (B_b(+a) - B_b(-a))
//       - a_b+ a_b- / (a_b+ + a_b-) * (B_a(+b) - B_a(-b))
//
// where a_a+- are the maxima of the a-normal face fan weights over the two
// a-faces adjacent to the edge (in b), a_b+- likewise for the b-faces,
// E(sa,sb) uses the corner velocity of the cell in quadrant (sa,sb) and the
// staggered fields B_b(sa) (from the b-face on side sa, reconstructed along a
// to the edge) and B_a(sb) (from the a-face on side sb, reconstructed along b).
// Corner velocities are v_cell + sa*0.5*slope_a + sb*0.5*slope_b with MC
// slopes, and staggered fields use MC slopes along the transverse direction.
//
// Pass 04 face-consistent correction (reconciliation R-UCT, M03-08). The
// primary REF profile solves HLLD on PPM (or PLM) states, so plain UCT-HLL on
// PLM-MC corner states does NOT reduce to the face flux actually used for the
// gas update. When MhdFaceSignalPatch::induction_correction_valid, the edge
// EMF is
//   E_c = E_UCT-HLL
//         - [a_b+ D_a(B_b)|_(b-) + a_b- D_a(B_b)|_(b+)] / (a_b+ + a_b-)
//         + [a_a+ D_b(B_a)|_(a-) + a_a- D_b(B_a)|_(a+)] / (a_a+ + a_a-),
// D_n(B_t) = F_face(B_t) - F_HLL^PLM(B_t) on the two adjacent faces of each
// family, upwind-weighted with the same signed fan weights UCT-HLL applies to
// the corresponding states. Derivation of the 1-D reduction: with no
// variation along b, both a-faces carry the same D_a and the same fan, the
// b-face L/R states coincide so D_b = F(U) - F(U) = 0, and E_UCT-HLL =
// -F_HLL^PLM(B_b); hence E_c = -F_HLL^PLM - D_a = -F_face(B_b) exactly (to
// roundoff) for HLLD+PPM, HLLD+PLM and every fallback flux. Symmetrically
// E_c = +F_face(B_a) for variation along b only. The
// correction is an edge-local function of face data, so D*C = 0 (divergence)
// is untouched; wall and periodic edge rules are applied afterwards as
// before. This is UCT-HLL plus a face-consistency correction, NOT the
// Mignone & Del Zanna HLLD-coefficient UCT; it is labelled R-UCT.
//
// One-dimensional reduction: with no variation along b the expression reduces
// exactly to -F_HLL(B_b) on the a-face (and with no variation along a, to
// +F_HLL(B_a) on the b-face), evaluated with PLM-MC edge states and the
// face's signed fan bounds. All four weights are non-negative and
// a_+ + a_- > 0 for any admissible state (fast speed > 0); a zero denominator
// is rejected as MHD_INVALID_FIELD rather than masked.
//
// Output EMFs are oriented edge line integrals of E_c = E_eff / sqrt(a) (the
// comoving B_c-CT EMF) times the comoving edge length, i.e. rates per unit
// time; the caller multiplies by the stage dt. Conducting-wall edges have
// their wall-tangential EMF set to zero; periodic self-wrapped patches copy
// the lower boundary edge plane to the upper one so the two storage copies of
// a logical edge are bitwise identical.
class MhdEmfBuilder {
 public:
  struct Diagnostics {
    std::uint64_t edges_built = 0U;
    std::uint64_t wall_edges_zeroed = 0U;
  };

  static void buildUctHll(const MhdPatchGeometry& geometry,
                          const MhdPatchBoundary& boundary,
                          const MhdPrimitivePatch& primitives,
                          const MhdStaggeredFieldPatch& staggered,
                          const MhdFaceSignalPatch& signals,
                          double scale_factor,
                          MhdEdgeEmfStorage& emf,
                          Diagnostics* diagnostics = nullptr);

  // Enforce the conducting-wall (E_t = 0) and periodic self-wrap edge rules on
  // an arbitrary edge-EMF field (shared by the ideal, Ohmic/ambipolar and Hall
  // operators so every EMF family satisfies the same boundary contract).
  static std::uint64_t applyEdgeBoundaryRules(const MhdPatchGeometry& geometry,
                                              const MhdPatchBoundary& boundary,
                                              MhdEdgeEmfStorage& emf);
};

// Inductive (curl) update of the magnetic face fluxes:
//   Phi <- Phi - weight * C(emf_rate).
void applyInductionUpdate(const MhdPatchGeometry& geometry,
                          std::span<double> phi_x,
                          std::span<double> phi_y,
                          std::span<double> phi_z,
                          const MhdEdgeEmfStorage& emf_rate,
                          double weight);

// Accumulate weight * interval EMF into a running edge-integral authority
// (used for AMR/MPI time-integrated EMF registers). No per-substep history.
void accumulateEdgeEmfIntegral(const MhdPatchGeometry& geometry,
                               MhdEdgeEmfStorage& accumulator,
                               const MhdEdgeEmfStorage& interval_emf,
                               double weight);

// Divergence normalizer (CHUI-MHD-REF v1.0 section 16.1, G03):
//   max(B_ref * A_ref, eps_scale), B_ref*A_ref = max |Phi|,
//   eps_scale = 1e-30 * B_unit_code * A_ref (B_unit_code = 1 in code units).
// A zero field therefore never reports an artificially excellent divergence
// through division by a large floor.
[[nodiscard]] double magneticDivergenceNormalizer(
    const MhdPatchGeometry& geometry,
    std::span<const double> phi_x,
    std::span<const double> phi_y,
    std::span<const double> phi_z);

}  // namespace cosmosim::mhd
