#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <span>
#include <string>

#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_state.hpp"

namespace cosmosim::mhd {

// Initial magnetic face flux construction (CHUI-MHD-REF v1.0 section 10):
//   Phi = Phi_uniform + C(A_edge),
// a uniform net-flux part plus the discrete curl of an edge-integrated vector
// potential. Because D*C = 0 the result is exactly (roundoff) divergence free
// on every patch, and periodic domains may carry a nonzero mean flux. All
// values are comoving: B_c = a_initial^2 B_phys. No cell-centred field is ever
// projected into face flux.

// Provenance of a seeded field (recorded into run provenance by the caller).
struct MhdSeedProvenance {
  std::string seed_model;
  std::array<double, 3> uniform_b_comoving{};
  double amplitude_comoving = 0.0;
  std::uint64_t random_seed = 0U;
  double scale_factor_initial = 1.0;
  double divergence_max_normalized = 0.0;
  std::uint32_t scheme_version = k_mhd_scheme_version;
};

// Line integral of the vector potential A_axis along one edge of the global
// comoving lattice, evaluated by the caller (e.g. analytic A or a deterministic
// spectral seed). Arguments: edge axis, comoving start point of the edge, edge
// length. Returns integral of A . dl in comoving code units.
using MhdEdgePotential =
    std::function<double(MhdAxis axis, std::array<double, 3> start_comoving, double length)>;

// Overwrite the patch faces with Phi_uniform + C(A). `potential` may be empty
// (uniform field only). The edge potential is evaluated from GLOBAL comoving
// positions so neighbouring patches produce identical shared faces.
// Returns the normalized divergence (G03 normalizer) of the result.
double seedFaceFlux(const MhdPatchGeometry& geometry,
                    std::array<double, 3> uniform_b_comoving,
                    const MhdEdgePotential& potential,
                    std::span<double> phi_x,
                    std::span<double> phi_y,
                    std::span<double> phi_z);

// Pass 04 (M03-10): composite-hierarchy edge authority. Wraps `inner` so that
// the line integral of any lattice edge (start on the lattice of a level l,
// length 2^(finest_depth - (l - base)) finest cells) is the PAIRWISE sum of
// `inner` evaluated on unit sub-edges of the finest lattice, whose start
// points are recomputed from integer finest-lattice coordinates
// origin + n * (base_spacing / 2^finest_depth). Consequently
//   I(coarse edge) == I(left fine half) + I(right fine half)
// with identical floating-point association, i.e. coarse and fine curl seeds
// are additive on the composite complex. `finest_depth` is the number of
// refinement levels below the base level present anywhere in the GLOBAL
// hierarchy (0 for a single level). Edges whose start/length are not on the
// finest lattice or whose length is not a power-of-two multiple of the finest
// spacing throw MHD_INVALID_FIELD (never silently resampled). Cost per edge
// is 2^(finest_depth - depth(edge)) inner evaluations; depth is limited to 20.
[[nodiscard]] MhdEdgePotential compositeEdgePotential(MhdEdgePotential inner,
                                                      std::array<double, 3> origin_comoving,
                                                      std::array<double, 3> base_spacing_comoving,
                                                      std::int32_t finest_depth);

// Gauss-quadrature (3-point) edge integral of an analytic potential A(x).
[[nodiscard]] double integrateEdgePotential(
    const std::function<std::array<double, 3>(std::array<double, 3>)>& vector_potential,
    MhdAxis axis,
    std::array<double, 3> start_comoving,
    double length);

// Field loop (Gardiner & Stone 2005): A_z = A0 (R - r) for r < R in the x-y
// plane centred at `center`; B is confined to the loop with |B| = A0.
[[nodiscard]] MhdEdgePotential fieldLoopPotential(double amplitude, double radius,
                                                  std::array<double, 3> center);

// Transverse part of a circularly polarized Alfven wave along x with
// k = 2 pi / wavelength: A = (0, amplitude sin(k x)/k, amplitude cos(k x)/k)
// gives B_y = -dA_z/dx = amplitude sin(k x), B_z = dA_y/dx = amplitude cos(k x).
// The guide field B_x is supplied through the uniform part of seedFaceFlux.
[[nodiscard]] MhdEdgePotential circularAlfvenPotential(double amplitude, double wavelength);

// Deterministic divergence-free stochastic seed: a finite sum of Fourier modes
// with Hermitian symmetry on the periodic box (box_comoving), modes |k| <=
// k_max (in units of 2 pi / L), random amplitudes from a SplitMix64 stream
// seeded by `seed`, each vector potential mode perpendicular to k. The field
// B = curl A is therefore solenoidal and the discrete faces divergence-free.
[[nodiscard]] MhdEdgePotential stochasticSeedPotential(std::array<double, 3> box_comoving,
                                                       std::int32_t k_max,
                                                       double amplitude,
                                                       double spectral_index,
                                                       std::uint64_t seed);

}  // namespace cosmosim::mhd
