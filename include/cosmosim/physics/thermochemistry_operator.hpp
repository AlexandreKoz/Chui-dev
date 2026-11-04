#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/gas_composition_sidecar.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/source_interval.hpp"
#include "cosmosim/physics/thermochemistry.hpp"

namespace cosmosim::core {
class SimulationState;
class MemoryGovernor;
}

namespace cosmosim::physics {

// ISM regime tags stored in GasCompositionSidecar::ism_regime_tag.
inline constexpr std::uint8_t k_ism_regime_resolved = core::k_gas_ism_regime_resolved;
inline constexpr std::uint8_t k_ism_regime_effective = core::k_gas_ism_regime_effective;

// Unit bridge between canonical code state and proper cgs physics.
struct ThermochemistryUnitBridge {
  double density_code_to_cgs = 1.0;          // UnitSystem::densityCodeToCgs(1)
  double velocity_code_to_cm_s = 1.0;        // velocity_si_per_code * 100
  double length_code_to_cm = 1.0;            // length_si_per_code * 100
  double mass_code_to_g = 1.0;
  bool comoving = false;                     // rho_phys = rho_code / a^3
};

// ---------------------------------------------------------------------------
// Local adaptive chemistry fidelity (CHUI Part 2 foundation; global error/
// resource-aware scheduling and statistical shadow qualification are Part 3).
// Requires the detailed layout. chemistry_model_tag: 0 = DETAILED (hybrid
// light network + C+/CO reactions), 1 = HYBRID (light network only, C+/CO
// lanes FROZEN). Demotion keeps the C+/CO lanes as an exact frozen
// reservoir (no species vanishes, nuclei and chemical energy unchanged), so
// promotion resumes from the stored molecular state (reversible mapping).
// Selection inputs: T, n_H, molecular fraction f_H2, CO fraction of carbon;
// promote if any cold/dense/molecular indicator fires, demote only when all
// are clear by the hysteresis margin and the dwell time has elapsed.
// Decisions are deterministic per row and taken only at the accepted source
// synchronization point (before the solve, never mid-solve).
// ---------------------------------------------------------------------------
inline constexpr std::uint8_t k_chemistry_tag_detailed = 0U;
inline constexpr std::uint8_t k_chemistry_tag_hybrid = 1U;

struct ChemistryAdaptiveConfig {
  bool enabled = false;
  double cold_temperature_k = 3.0e3;
  double dense_n_h_cm3 = 1.0;
  double molecular_fraction_threshold = 1.0e-3;
  double co_carbon_fraction_threshold = 1.0e-6;
  double hysteresis_fraction = 0.3;
  double minimum_dwell_years = 1.0e6;
  // R3 (C03): timescale/stiffness criterion. Promote when the network's own
  // net cooling time t_cool = u / (Lambda + chem - Gamma) is shorter than
  // cooling_time_ratio x the source interval (stiff, rapidly cooling gas in
  // which the composition-dependent cooling of the detailed network must be
  // followed); demotion additionally requires t_cool > ratio dt (1 + h).
  // 0 disables. EXPERIMENTAL/UNQUALIFIED (uncalibrated).
  double cooling_time_ratio = 0.0;
  // Previous accepted, normalized integrator-error and conservative residual
  // history. Zero disables. These use persisted per-cell history and thus do
  // not depend on row order, rank, retry, or restart.
  double error_estimate_threshold = 0.0;
  double residual_fraction_threshold = 0.0;
};

struct ChemistryAdaptiveReport {
  std::uint64_t promoted = 0U;
  std::uint64_t demoted = 0U;
  std::uint64_t held_by_dwell = 0U;
  std::uint64_t promoted_by_cooling_time = 0U;  // R3: promotions where the timescale criterion fired
  std::uint64_t promoted_by_error_history = 0U;
  std::uint64_t detailed_cells = 0U;
  std::uint64_t hybrid_cells = 0U;
};

// Bounded shadow reference: deterministic sample keyed by (gas_cell_id,
// accepted step); runs the DETAILED solver at 0.1x tolerance on a private
// copy of the cell input and compares with the production result. Never
// modifies state; capped by max_cells per stage; no domain-sized arrays.
struct ChemistryShadowConfig {
  double sample_fraction = 0.0;
  std::uint32_t max_cells = 0U;
};
struct ChemistryShadowReport {
  std::uint64_t samples = 0U;
  double max_relative_species_difference = 0.0;
  double max_relative_electron_difference = 0.0;
  double max_relative_energy_difference = 0.0;
  double mean_substep_cost_ratio = 0.0;  // reference / production substeps
};

struct ThermochemistryOperatorConfig {
  ThermochemistryConfig solver{};
  ThermochemistryUnitBridge units{};
  // Live dust lanes present => D/D_MW from the solid mass; otherwise, when
  // true, D/D_MW = Z/Z_sun (explicit metallicity-scaled dust proxy; never
  // both, so dust catalysis/heating cannot be double counted).
  bool metallicity_scaled_dust_proxy = false;
  // Milky Way dust-to-gas mass ratio used to normalize D (Draine 2011).
  double milky_way_dust_to_gas = 0.01;
  ChemistryAdaptiveConfig adaptive{};
  ChemistryShadowConfig shadow{};
};

// Staged per-row result (bounded: one record per processed row).
struct ThermochemistryStagedRow {
  std::uint32_t row = 0;
  double internal_energy_code = 0.0;
  double mean_molecular_weight = 0.0;
  double temperature_k = 0.0;
  NetworkState species_mass_code{};  // species lane masses (network backends)
};

// Result of advancing one row (pure: no state mutation).
struct ThermochemistryRowResult {
  SourceAdvanceStatus status = SourceAdvanceStatus::kCompleted;
  std::string diagnostic;
  double suggested_retry_seconds = 0.0;
  ThermochemistryStagedRow staged{};
  SourceLedger ledger{};  // volume-integrated [erg]
};

// Provisional chemistry-fidelity decision of one row (pure).
struct ChemistryFidelityDecision {
  std::uint8_t tag = k_chemistry_tag_detailed;
  double dwell_years = 0.0;
  bool promoted = false;
  bool demoted = false;
  bool held_by_dwell = false;
  bool cooling_time_trigger = false;  // R3: the timescale criterion requested detailed chemistry
  bool error_history_trigger = false;
  double cooling_time_seconds = 0.0;  // diagnostic (+inf when net heating)
};

// ---------------------------------------------------------------------------
// Thermochemistry kernels for the SourceRuntime (CHUI Part 2). When a Part 2
// backend is selected, the hydro in-solver CoolingHeatingSource is NOT
// registered, so cooling is advanced exactly once per accepted interval.
// All member functions are PURE with respect to SimulationState: the
// transactional owner (BaryonicSourceStage) stages every result and writes
// state only after the collective accept (local source-stage atomicity).
// Effective-regime rows are never passed here by the stage (the effective
// multiphase closure owns their thermal state — no double cooling).
// ---------------------------------------------------------------------------
class ThermochemistryOperator {
 public:
  ThermochemistryOperator(ThermochemistryOperatorConfig config, CompositionLayout layout);

  // Advance one row with the given (possibly provisional) fidelity tag.
  [[nodiscard]] ThermochemistryRowResult advanceRow(const core::SimulationState& state,
                                                    const SourceInterval& interval, std::uint32_t row,
                                                    std::uint8_t chemistry_tag,
                                                    ChemistryShadowReport& shadow) const;
  // Local fidelity decision from the row's step-start state (no-op decision
  // returning the current tag unless adaptive.enabled).
  [[nodiscard]] ChemistryFidelityDecision decideChemistryFidelity(const core::SimulationState& state,
                                                                  const SourceInterval& interval,
                                                                  std::uint32_t row) const;

  // CHUI Part 3 (C03/C04 RESEARCH, EXPERIMENTAL/UNQUALIFIED): globally
  // synchronized, error-aware widening of the detailed-chemistry region.
  // promotion_bias b >= 0 scales the promotion thresholds by s = 1 + b
  // (cold: T < T_cold s, dense: n_H > n_dense / s, molecular / CO fractions
  // divided by s) and the demotion thresholds identically, so hysteresis is
  // preserved. The runtime derives b from the GLOBAL (MPI-max) sampled
  // shadow disagreement of the previous accepted step, identically on every
  // rank, and persists it (restart-reproducible). b = 0 is the Part 2 model.
  void setGlobalPromotionBias(double bias);
  [[nodiscard]] double globalPromotionBias() const noexcept { return m_promotion_bias; }

  [[nodiscard]] const ThermochemistrySolver& solver() const noexcept { return m_solver; }
  [[nodiscard]] const CompositionLayout& layout() const noexcept { return m_layout; }
  [[nodiscard]] const ThermochemistryOperatorConfig& config() const noexcept { return m_config; }

  // Build the solver input of one row (also used by the shadow sampler).
  [[nodiscard]] ThermochemistryCellInput cellInput(const core::SimulationState& state, std::uint32_t row,
                                                   const SourceInterval& interval) const;

 private:
  ThermochemistryOperatorConfig m_config;
  CompositionLayout m_layout;
  ThermochemistrySolver m_solver;
  std::unique_ptr<ThermochemistrySolver> m_hybrid_solver;     // adaptive demoted rows
  std::unique_ptr<ThermochemistrySolver> m_reference_solver;  // shadow reference
  double m_promotion_bias = 0.0;
};

// Deterministic control law of the global promotion bias from the global
// maximum relative shadow species difference (exposed for tests):
//   d > 0.1  -> b + 0.5 (capped at 4),  d < 0.01 -> b - 0.25 (floored at 0),
//   otherwise unchanged; no samples -> unchanged.
[[nodiscard]] double nextChemistryPromotionBias(double bias, double global_max_shadow_difference,
                                                std::uint64_t global_samples) noexcept;

}  // namespace cosmosim::physics
