#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "cosmosim/physics/source_interval.hpp"
#include "cosmosim/physics/thermochemistry_operator.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// ===========================================================================
// 1) ISM regime manager (CHUI Part 2, ism_regime_model=effective_or_resolved)
//
// Regimes: EFFECTIVE (unresolved multiphase; the SH03-type effective EOS and
// its relaxation own the thermal state; resolved thermochemistry skipped;
// effective SF prescription) versus RESOLVED (physical thermochemistry,
// resolved-collapse SF). Criterion (Truelove et al. 1997 Jeans condition):
//   lambda_J = sqrt(pi c_s^2 / (G rho_phys)),  resolved iff lambda_J >= N_J dx
// with c_s from the thermal state and the composition mu, dx = V^(1/3)
// proper. A cell is EFFECTIVE iff n_H >= n_th (effective threshold) AND the
// Jeans length is unresolved. Hysteresis: enter effective below
// N_J (1 - h) dx / leave above N_J (1 + h) dx (density gate likewise), and a
// minimum dwell time. Decisions are local and deterministic, made only at
// the accepted source synchronization point (never mid-solve).
// Transition consistency: mass, elements, species, SGS energy and the
// canonical thermal energy e are unchanged by a switch; only ownership
// changes (P_eff = (gamma-1) rho u_eff equals the ideal pressure at the same
// u, so pressure is continuous across the switch).
// ===========================================================================
struct IsmRegimeConfig {
  double effective_threshold_n_h_cm3 = 0.13;
  double jeans_resolution_cells = 4.0;
  double hysteresis_fraction = 0.2;
  double minimum_dwell_years = 1.0e6;
  double adiabatic_index = 5.0 / 3.0;
};

struct IsmRegimeReport {
  std::uint64_t effective_cells = 0U;
  std::uint64_t resolved_cells = 0U;
  std::uint64_t to_effective = 0U;
  std::uint64_t to_resolved = 0U;
  std::uint64_t held_by_dwell = 0U;
  std::vector<std::uint32_t> changed_rows;
};

// Provisional decision for one row from its step-start state (PURE).
struct IsmRegimeDecision {
  std::uint8_t tag = 0U;
  double dwell_years = 0.0;
  bool changed = false;
  bool held_by_dwell = false;
};
[[nodiscard]] IsmRegimeDecision decideIsmRegime(const core::SimulationState& state, const IsmRegimeConfig& config,
                                                const ThermochemistryUnitBridge& units,
                                                const SourceInterval& interval, std::uint32_t row);

// Evaluate and APPLY regime tags for `rows` (all rows when empty). Direct
// utility for tests/tools; the production source stage stages decisions
// through BaryonicSourceStage and applies them only on commit.
IsmRegimeReport updateIsmRegimes(core::SimulationState& state, const IsmRegimeConfig& config,
                                 const ThermochemistryUnitBridge& units, const SourceInterval& interval,
                                 std::span<const std::uint32_t> rows);
[[nodiscard]] double jeansLengthCm(double temperature_k, double mean_molecular_weight, double rho_cgs,
                                   double adiabatic_index, double extra_velocity_dispersion_sq_cm2_s2 = 0.0);

// ===========================================================================
// 2) Live dust: two compositions (silicate MgFeSiO4, carbonaceous C).
//
// Processes (all exact, element-conserving transfers between element lanes
// and dust lanes in the reservoir stoichiometry; nuclei are never destroyed):
//  * accretion (grain growth), Hirashita (2000, PASJ 52, 585) / Asano et al.
//    (2013, EP&S 65, 213) scaling
//      tau_acc = tau_0 (n_H/100 cm^-3)^-1 (T/50 K)^-1/2 (Z/Z_sun)^-1
//      dD/dt = D (1 - D/D_max)/tau_acc  (logistic; exact solution used)
//    D_max = D + condensable gas mass limited by the key element;
//  * thermal sputtering, Tsai & Mathews (1995, ApJ 448, 84)
//      tau_sp = tau_sp,0 (n_H/1e-3 cm^-3)^-1 [(2e6 K / T)^2.5 + 1]
//      dD/dt = -3 D / tau_sp  (exact exponential);
//  * SN destruction, McKee (1989): destroyed = eps M_s(E/1e51 erg) D/M_gas,
//    applied by feedback deposition (once per committed SN);
//  * stellar condensation: a fraction delta of returned refractory elements
//    (AGB / core-collapse efficiencies, Dwek 1998, ApJ 501, 643 parameter
//    form) enters the dust lanes at injection — the same mass is never also
//    injected into gas.
// Default efficiencies/timescales are literature-scale RESEARCH defaults and
// are NOT calibrated (EXPERIMENTAL/UNQUALIFIED).
// Not modelled (named): grain size distribution, shattering, coagulation,
// non-thermal sputtering in shocks, grain charging.
// ===========================================================================
struct DustConfig {
  double accretion_timescale_ref_s = 2.0e7 * k_seconds_per_year;
  double sputtering_timescale_ref_s = 1.7e8 * k_seconds_per_year;
  double sn_swept_mass_per_1e51_erg_g = 6800.0 * k_solar_mass_g;
  double sn_destruction_efficiency = 0.3;
  double agb_condensation_efficiency = 0.2;
  double cc_condensation_efficiency = 0.1;
};

struct DustUpdateLedger {
  double gas_to_dust_code = 0.0;   // accretion
  double dust_to_gas_code = 0.0;   // sputtering
  std::uint64_t cells = 0U;
};

// One cell, isochoric, over dt seconds. Lanes are modified in place through
// the sidecar (caller stages by working on a copy when transactional).
void advanceDustCell(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout, std::size_t row,
                     double n_h_cm3, double temperature_k, double dt_seconds, const DustConfig& config,
                     DustUpdateLedger& ledger);
// Returns destroyed dust mass (code) returned to the gas element lanes.
double destroyDustBySupernova(core::GasCompositionSidecar& sidecar, const CompositionLayout& layout,
                              std::size_t row, double swept_gas_mass_code, double cell_gas_mass_code);
// Split stellar ejecta element masses into gas elements and dust lanes.
struct EjectaPartition {
  std::array<double, k_chemical_element_count> gas_elements{};
  std::array<double, k_dust_reservoir_count> dust{};
};
[[nodiscard]] EjectaPartition partitionEjectaIntoDust(const std::array<double, k_chemical_element_count>& ejecta,
                                                      double condensation_efficiency);

// ===========================================================================
// 3) SGS turbulent energy, budgeted_v1 (RESEARCH CLOSURE, uncalibrated).
//   E_sgs = M K (extensive lane). Sources: an explicit fraction of feedback
//   energy (taken FROM the SN budget, never added on top), advection/mixing
//   with the composition transport. Sink: dissipation into heat,
//     dK/dt = -C_eps K^{3/2} / dx   =>  K(t) = K0 / (1 + C_eps sqrt(K0) t/(2 dx))^2
//   with the dissipated energy M (K0 - K) added to the canonical thermal e
//   (exact transfer, ledgered). Resolved shear uses a traceless
//   Smagorinsky stress nu_t = C_v Delta^2 |S*| on the authoritative AMR face
//   graph. Its pairwise impulse and stress work update resolved momentum and
//   energy; the positive production rho nu_t |S*|^2 is removed from resolved
//   total energy and added to E_sgs in the same staged commit. Pressure
//   support: K enters ONLY the
//   star-formation Jeans/virial support as sigma^2 = (2/3) K; it is NOT added
//   to the hydro EOS pressure (no double counting). The trace is excluded:
//   isotropic SGS pressure remains deliberately outside this research model.
// ===========================================================================
struct SgsConfig {
  double dissipation_coefficient = 1.58;  // Schmidt, Niemeyer & Hillebrandt 2006
  double eddy_viscosity_coefficient = 0.05;  // uncalibrated C_v (Smagorinsky form)
};
struct SgsLedger {
  double dissipated_to_heat_code = 0.0;
  double injected_from_feedback_code = 0.0;
  double produced_from_resolved_code = 0.0;
  double resolved_stress_work_code = 0.0;
  double momentum_impulse_l1_code = 0.0;
  std::uint64_t shear_faces = 0U;
  std::uint64_t shear_cells = 0U;
};

// Local half of a face whose other authoritative cell resides on another
// rank. Geometry is expressed in the same physical-code length convention as
// MetalDiffusionFace. local_side is +1 when the remote cell is on the +axis
// side. The stable ids/rank make the distributed topology agreement
// auditable; the remote state is the owner's step-start payload.
struct SgsRemoteShearFace {
  std::uint32_t local_row = 0U;
  std::uint32_t remote_index = 0U;
  std::uint64_t local_gas_cell_id = 0U;
  std::uint64_t remote_gas_cell_id = 0U;
  int remote_rank = 0;
  std::uint8_t axis = 0U;
  std::int8_t local_side = 0;
  bool periodic = false;
  double area_code = 0.0;
  double center_distance_code = 0.0;
  double remote_density_physical_code = 0.0;
  std::array<double, 3> remote_velocity_peculiar{};
};

// Five full-row staging lanes. prepareSgsShearTransfer never mutates state;
// commitSgsShearTransfer is allocation-free and cannot fail after validation.
// This split lets SourceRuntime coordinate the MPI accept/reject vote before
// any rank publishes momentum, thermal energy, or E_sgs.
struct SgsShearTransferWorkspace {
  std::vector<double> velocity_x;
  std::vector<double> velocity_y;
  std::vector<double> velocity_z;
  std::vector<double> internal_energy;
  std::vector<double> sgs_energy;
  std::vector<std::uint8_t> touched;
  bool prepared = false;

  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
  [[nodiscard]] static std::uint64_t requiredBytes(std::size_t rows) noexcept;
  void clear() noexcept;
};

[[nodiscard]] SgsLedger prepareSgsShearTransfer(
    const core::SimulationState& state,
    std::span<const std::uint8_t> active_rows,
    std::span<const double> trace_free_strain_magnitude_code,
    std::span<const MetalDiffusionFace> local_faces,
    std::span<const SgsRemoteShearFace> remote_faces,
    double physical_length_per_stored_length,
    double dt_code_time,
    const SgsConfig& config,
    SgsShearTransferWorkspace& workspace);

void commitSgsShearTransfer(
    core::SimulationState& state,
    SgsShearTransferWorkspace& workspace) noexcept;
// Returns dissipated energy (code units) and updates e and the SGS lane.
double advanceSgsCell(core::SimulationState& state, const CompositionLayout& layout, std::size_t row,
                      double dt_code_time, double cell_size_code, const SgsConfig& config, SgsLedger& ledger);

// Pure lane-level dissipation: returns {new extensive SGS energy, dissipated}.
struct SgsDissipation {
  double sgs_energy = 0.0;
  double dissipated = 0.0;
};
[[nodiscard]] SgsDissipation sgsDissipate(double sgs_energy, double mass, double dt_code_time,
                                          double cell_size_code, const SgsConfig& config) noexcept;

// Exact analytic SGS decay (exposed for tests).
[[nodiscard]] double sgsSpecificEnergyAfter(double k0, double dt, double dx, double c_eps) noexcept;

}  // namespace cosmosim::physics
