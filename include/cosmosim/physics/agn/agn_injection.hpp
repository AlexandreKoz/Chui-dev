#pragma once

// CHUI-RT-BH-AGN-REF-v1, AGN companion (AGN01-AGN05): conservative injection
// of the BH engine's pending channel budgets into owned leaf gas.
//
// Budgets come only from the chui_bh_engine_v1 pending lanes (the injector
// never recomputes a luminosity). Channels (P0):
//  * thermal   continuous: dE_i = omega_i/sum(omega) Q_th (omega_i = integral_cell K dV);
//              reservoir : released when Q >= dT sum_i w_i m_i c_v,i with c_v,i = u_i/T_i
//              (actual cell EOS ratio; ideal-gas constant c_v assumption) and distributed
//              proportional to w_i m_i c_v,i (uniform dT on the kernel-weighted mass).
//  * kicks     radial pattern q_i = omega_i r_hat_i; amplitude from the energy-consistent
//              quadratic C A^2 + B A = Q including pre-existing gas velocities and the BH
//              recoil (B -= Q_net . v_bh, C += |Q_net|^2/(2 M_bh)); stable root.
//  * jets      bipolar cones (half opening angle) about the disk angular-momentum axis;
//              each lobe normalized independently so the lobe momenta are exactly equal
//              and opposite (zero recoil); amplitude from the same quadratic. Mass is the
//              ambient gas (entrainment transfers momentum, never creates mass).
//  * winds     mass from the disk-wind reservoir (pending_wind_mass), radial in the BH
//              rest frame at v_w; dp = dm v_w r_hat + dm v_bh, dE = dm (u_w + |v_w|^2/2)
//              + v_bh . dp_rest + dm |v_bh|^2/2; u_w from the energy remainder; if the
//              budget cannot launch the mass at v_w only the launchable mass is ejected
//              (rest stays pending, flagged). Recoil R = sum dm v_w r_hat goes to the BH
//              and its kinetic change is charged to the wind budget.
//  * radiation none -> recorded as escaping; local_empirical -> declared fraction to the
//              thermal channel (no RT); rt_transport -> left pending for the RT source
//              provider (transported exactly once).
// The per-cell update is the conservative transaction convention: the thermal
// increment is dE - dK with dK = |p+dp|^2/(2(m+dm)) - |p|^2/(2m) (never the
// requested kinetic energy deposited again as heat). Magnetic faces are never
// touched (MHD authority). Nonrelativistic domain: any cell velocity increment
// above max_velocity_c_fraction c holds the channel (rejected, never capped).
// Missing support (no owned leaf cells, empty lobe, undefined axis) retains
// the budget with a diagnostic flag.

#include <array>
#include <cstdint>
#include <vector>

#include "cosmosim/core/config_rt_bh_agn.hpp"
#include "cosmosim/physics/black_hole/bh_neighborhood.hpp"
#include "cosmosim/physics/rt_bh_agn_frame.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// Nonnegative root of C A^2 + B A = Q (Q >= 0, C >= 0). Returns false when
// no nonnegative root exists (C = 0 and B <= 0 with Q > 0).
[[nodiscard]] bool solveEnergyConsistentAmplitude(double c, double b, double q, double& amplitude) noexcept;
[[nodiscard]] double cellKineticIncrement(double m, const std::array<double, 3>& p, double dm,
                                          const std::array<double, 3>& dp) noexcept;

struct AgnRefinementRequest {
  std::uint64_t bh_id = 0U;
  std::array<double, 3> position_stored{};
  double support_radius_stored = 0.0;
  double required_cell_width_proper = 0.0;  // injection diameter / N_j
  double current_cell_width_proper = 0.0;
};

struct AgnInjectionReport {
  std::uint64_t sources = 0U;
  std::uint64_t held_no_support = 0U;
  std::uint64_t held_lobe_imbalance = 0U;
  std::uint64_t held_velocity_domain = 0U;
  std::uint64_t held_axis_undefined = 0U;
  std::uint64_t reservoir_waiting = 0U;
  std::uint64_t unresolved_injection = 0U;
  double thermal_injected_code = 0.0;
  double kinetic_budget_injected_code = 0.0;
  double jet_budget_injected_code = 0.0;
  double wind_budget_injected_code = 0.0;
  double wind_mass_injected_code = 0.0;
  double radiation_escaped_code = 0.0;
  double radiation_locally_coupled_code = 0.0;
  double max_energy_closure_residual = 0.0;
  double max_lobe_imbalance = 0.0;
  std::vector<std::uint32_t> touched_gas_rows;
  std::vector<AgnRefinementRequest> refinement_requests;
  bool gas_mutated = false;
};

struct AgnInjectionContext {
  ProperFrameContext frame{};
  double time_end_code = 0.0;
  std::uint64_t interval_index = 0U;
  BhOwnedLeafView leaves{};
  PeriodicBox box{};
  double mass_floor_code = 1.0e-20;
  double specific_internal_energy_floor_code = 0.0;
};

class AgnFeedbackInjector {
 public:
  explicit AgnFeedbackInjector(core::AgnFeedbackConfig config);
  // Collective. Reads/writes pending budgets in the chui_bh_engine_v1 store.
  [[nodiscard]] AgnInjectionReport apply(core::SimulationState& state, const AgnInjectionContext& context,
                                         BhSupportExchange& exchange) const;
  [[nodiscard]] const core::AgnFeedbackConfig& config() const noexcept { return m_config; }

 private:
  core::AgnFeedbackConfig m_config;
};

}  // namespace cosmosim::physics
