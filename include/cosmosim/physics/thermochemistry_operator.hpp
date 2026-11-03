#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/source_interval.hpp"
#include "cosmosim/physics/thermochemistry.hpp"

namespace cosmosim::core {
class SimulationState;
class MemoryGovernor;
}

namespace cosmosim::physics {

// ISM regime tags stored in GasCompositionSidecar::ism_regime_tag.
inline constexpr std::uint8_t k_ism_regime_resolved = 0U;
inline constexpr std::uint8_t k_ism_regime_effective = 1U;

// Unit bridge between canonical code state and proper cgs physics.
struct ThermochemistryUnitBridge {
  double density_code_to_cgs = 1.0;          // UnitSystem::densityCodeToCgs(1)
  double velocity_code_to_cm_s = 1.0;        // velocity_si_per_code * 100
  double length_code_to_cm = 1.0;            // length_si_per_code * 100
  double mass_code_to_g = 1.0;
  bool comoving = false;                     // rho_phys = rho_code / a^3
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
};

// Staged per-row result (bounded: one record per processed row).
struct ThermochemistryStagedRow {
  std::uint32_t row = 0;
  double internal_energy_code = 0.0;
  double mean_molecular_weight = 0.0;
  double temperature_k = 0.0;
  NetworkState species_mass_code{};  // species lane masses (network backends)
};

// ---------------------------------------------------------------------------
// One owner for thermochemistry in the SourceRuntime (CHUI Part 2). When a
// Part 2 backend is selected, the hydro in-solver CoolingHeatingSource is NOT
// registered, so cooling is advanced exactly once per accepted interval.
// prepare(): reads canonical state (e, rho, composition), integrates every
//   resolved-regime row over the physical interval, stages the results; any
//   row that cannot complete => SourceIntervalIncompleteError (whole stage
//   rejected; nothing staged is committed). Effective-regime rows are skipped
//   entirely (the effective multiphase closure owns their thermal state — the
//   no-double-cooling rule).
// commit(): writes e, species lanes and the mu cache for staged rows, then
//   returns the touched rows for the thermodynamic refresh (F4 contract).
// ---------------------------------------------------------------------------
class ThermochemistryOperator {
 public:
  ThermochemistryOperator(ThermochemistryOperatorConfig config, CompositionLayout layout);

  [[nodiscard]] static std::uint64_t stagingBytes(std::size_t rows) noexcept {
    return static_cast<std::uint64_t>(rows) * sizeof(ThermochemistryStagedRow);
  }

  // rows empty => all gas rows.
  void prepare(const core::SimulationState& state, const SourceInterval& interval,
               std::span<const std::uint32_t> rows);
  std::vector<std::uint32_t> commit(core::SimulationState& state);

  [[nodiscard]] const SourceLedger& ledger() const noexcept { return m_ledger; }
  [[nodiscard]] const ThermochemistrySolver& solver() const noexcept { return m_solver; }
  [[nodiscard]] std::uint64_t skippedEffectiveRows() const noexcept { return m_skipped_effective; }

  // Build the solver input of one row (also used by the shadow sampler).
  [[nodiscard]] ThermochemistryCellInput cellInput(const core::SimulationState& state, std::uint32_t row,
                                                   const SourceInterval& interval) const;

 private:
  ThermochemistryOperatorConfig m_config;
  CompositionLayout m_layout;
  ThermochemistrySolver m_solver;
  std::vector<ThermochemistryStagedRow> m_staged;
  SourceLedger m_ledger;
  std::uint64_t m_skipped_effective = 0U;
  bool m_prepared = false;
};

}  // namespace cosmosim::physics
