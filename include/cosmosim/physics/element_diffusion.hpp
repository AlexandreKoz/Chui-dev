#pragma once

// CHUI Part 3 Checkpoint C — element-resolved conservative turbulent mixing.
//
// Physics: unresolved turbulence mixes gas-phase material and the dust that
// is dynamically coupled to it. Every MASS lane of the composition sidecar
// (eleven element lanes, every tracked species lane and every dust lane) is a
// passive mass fraction x_l = m_l / M with M = sum(elements) + sum(dust) (the
// composition gas mass), evolved with the SAME Smagorinsky diffusivity as the
// historical total-Z operator:
//
//   dm_{i,l}/dt = - sum_f C_f (x_{i,l} - x_{j,l}),   C_f = A_f <rho kappa>_H / d_f
//
// (<.>_H: harmonic mean of the two cells; the same conductance as
// MetalDiffusionModel). The SGS energy lane is NOT a mass and is not mixed by
// this operator (its own model owns its transport).
//
// Properties that hold by construction (not by tolerance):
//  * pairwise equal and opposite face transfers: every lane total is
//    conserved exactly up to roundoff;
//  * one linear operator with identical weights for all lanes: every linear
//    constraint between lanes is preserved, in particular
//    sum(species of element e) == element lane e and the composition gas mass
//    of each cell (sum_l x_l = 1 on both sides of each face);
//  * explicit positivity: with dt_sub <= cfl * min_i M_i / sum_{f in i} C_f
//    (cfl <= 1) each new lane value is a convex combination, so no lane can
//    become negative and no limiter (which would break the linear
//    constraints) is needed. A step needing more than max_subcycles fails
//    closed (the caller rejects the source interval collectively).
//  * open and reflective boundary faces carry zero diffusive flux (no-flux
//    mixing boundary); periodic faces are ordinary internal pairs.
//
// The operator works on a scratch copy and writes the sidecar only after the
// whole interval succeeded (transactional; a refusal leaves it unchanged).

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/gas_composition_sidecar.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"

namespace cosmosim::physics {

struct ElementDiffusionReport {
  std::uint32_t subcycles = 0U;
  std::uint64_t faces_evaluated = 0U;
  double stable_dt_code = 0.0;
  std::size_t mixed_lane_count = 0U;
  double max_relative_lane_residual = 0.0;
  double max_relative_cell_mass_change = 0.0;
};

[[nodiscard]] std::uint64_t elementDiffusionWorkspaceBytes(std::size_t lanes, std::size_t rows) noexcept;

[[nodiscard]] ElementDiffusionReport advanceElementResolvedDiffusion(
    core::GasCompositionSidecar& sidecar,
    std::span<const MetalDiffusionFace> faces,
    std::span<const double> rho_kappa_code,
    double dt_code,
    double parabolic_cfl,
    std::uint32_t max_subcycles,
    std::vector<double>& workspace);

}  // namespace cosmosim::physics
