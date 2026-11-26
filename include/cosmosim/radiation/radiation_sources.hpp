#pragma once

// CHUI-RT-BH-AGN-REF-v1, RT companion section 7 / AGN companion section 5
// (RT07, AGN01): time-integrated radiation sources.
//
// AGN: the radiative budget accumulated by the BH engine in
// chui_bh_engine_v1::pending_radiation_energy_code (accretion powered, never
// recomputed from BH mass) is emitted EXACTLY ONCE over the accepted
// interval: luminosity L = E_pending / dt, group split by the normalized AGN
// spectral shape (L_g = L int_g phi), deposited over the leaf cells of the
// radiation records inside the BH kernel with cell-volume kernel weights
// (sum of weights = 1), isotropic in the source rest frame (no net emitted
// momentum; moving-source recoil O(v/c) is a named neglected term). The
// pending budget is debited at the same time (restart truth); if no record
// cell lies inside the kernel the budget is retained and flagged.
// Stellar sources: FAIL_CLOSED (no age/metallicity luminosity provider).

#include <cstdint>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"
#include "cosmosim/radiation/radiation_runtime.hpp"

namespace cosmosim::radiation {

struct AgnRadiationSourceReport {
  std::uint64_t sources = 0U;
  std::uint64_t held_no_support = 0U;
  double emitted_energy_code = 0.0;  // proper code energy scheduled for the interval
};

// Builds deposits for every local BH with a pending radiative budget and
// debits the engine store. `group_fractions` = radiationGroupLuminosityFractions.
[[nodiscard]] AgnRadiationSourceReport buildAgnRadiationSources(core::SimulationState& state,
                                                                const std::vector<double>& group_fractions,
                                                                const physics::ProperFrameContext& frame,
                                                                const physics::PeriodicBox& box, double dt_proper_code,
                                                                std::vector<RadiationSourceDeposit>& deposits);

}  // namespace cosmosim::radiation
