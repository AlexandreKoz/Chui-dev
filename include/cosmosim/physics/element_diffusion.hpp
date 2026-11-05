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
#include <exception>
#include <span>
#include <vector>

#include "cosmosim/core/gas_composition_sidecar.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"

namespace cosmosim::core {
class MemoryGovernor;
}

namespace cosmosim::physics {

struct ElementDiffusionReport {
  std::uint32_t subcycles = 0U;
  std::uint64_t faces_evaluated = 0U;
  double stable_dt_code = 0.0;
  std::size_t mixed_lane_count = 0U;
  double max_relative_lane_residual = 0.0;
  double max_relative_cell_mass_change = 0.0;
};

// G0.5: the ONE admission estimate for advanceElementResolvedDiffusion's
// simultaneously live scratch: values (lanes x rows), per-row gas mass and
// conductance sum, the per-subcycle delta (lanes x rows), the lane-index,
// before and scale vectors (O(lanes)). `lanes` may be the full sidecar lane
// count (an upper bound of the mixed lanes). Callers must NOT add the delta
// separately. Remote-coupling buffers are admitted by the coupling owner.
// Saturates at UINT64_MAX instead of wrapping.
[[nodiscard]] std::uint64_t elementDiffusionWorkspaceBytes(std::size_t lanes, std::size_t rows) noexcept;

// G2.1 distributed coupling. The operator is MPI-free; a distributed caller
// supplies this interface. Contract (every method is COLLECTIVE: it is
// entered by every rank the same number of times in the same order, also by
// ranks without rows or faces):
//  * reduceFailure(local): if any rank holds a failure, rethrows on every
//    rank (the local exception where present, a peer-failure error
//    otherwise). Called before every other collective so a local error never
//    leaves peers blocked in an exchange.
//  * addRemoteConductance(sum): adds, for every local row, the conductance
//    of its remote interface faces (identical values on both sides).
//  * globalMin / globalSum: deterministic reductions.
//  * addRemoteFaceTransfers(values, gas_mass, lanes, dt_sub, delta): exchanges
//    FRESH local boundary-row lane values (this subcycle), lets the unique
//    face authority compute every remote face transfer exactly once, and adds
//    the resulting equal-and-opposite contributions to the local rows of
//    `delta` (layout lane-major: delta[k * rows + row]).
class ElementDiffusionRemoteCoupling {
 public:
  virtual ~ElementDiffusionRemoteCoupling() = default;
  virtual void reduceFailure(std::exception_ptr local_failure, const char* stage) = 0;
  virtual void addRemoteConductance(std::span<double> conductance_sum) = 0;
  [[nodiscard]] virtual double globalMin(double value) = 0;
  virtual void globalSum(std::span<double> values) = 0;
  virtual void addRemoteFaceTransfers(std::span<const double> values, std::span<const double> gas_mass,
                                      std::size_t lane_count, double dt_sub, std::span<double> delta) = 0;
};

// Local (coupling == nullptr) or distributed element-resolved diffusion.
// With a coupling, stable dt and the subcycle count are global (identical on
// every rank), positivity and lane conservation are verified globally, and
// every rank commits or every rank throws.
[[nodiscard]] ElementDiffusionReport advanceElementResolvedDiffusion(
    core::GasCompositionSidecar& sidecar,
    std::span<const MetalDiffusionFace> faces,
    std::span<const double> rho_kappa_code,
    double dt_code,
    double parabolic_cfl,
    std::uint32_t max_subcycles,
    std::vector<double>& workspace,
    ElementDiffusionRemoteCoupling* coupling = nullptr,
    core::MemoryGovernor* memory_governor = nullptr);

}  // namespace cosmosim::physics
