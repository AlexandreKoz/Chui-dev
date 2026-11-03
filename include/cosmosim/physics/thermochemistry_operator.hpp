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
};

struct ChemistryAdaptiveReport {
  std::uint64_t promoted = 0U;
  std::uint64_t demoted = 0U;
  std::uint64_t held_by_dwell = 0U;
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
};

}  // namespace cosmosim::physics
