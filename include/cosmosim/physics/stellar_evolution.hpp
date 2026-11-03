#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

enum class StellarYieldChannel : std::uint8_t {
  kWindsAgb = 0,
  kCoreCollapseSn = 1,
  kTypeIaSn = 2,
  // Source compatibility alias.
  kWinds = kWindsAgb,
};

constexpr std::size_t k_stellar_yield_channel_count = 3;

struct StellarEvolutionCumulative {
  double return_fraction_total = 0.0;
  // Total ejected metal mass, including returned birth-composition metals.
  double total_ejected_metal_fraction_total = 0.0;
  // Net newly synthesized metal mass. This is not the total ejecta.
  double newly_synthesized_metal_fraction_total = 0.0;
  double event_count_per_initial_mass_code_total = 0.0;
  double energy_erg_per_initial_mass_code = 0.0;
  std::array<double, k_stellar_yield_channel_count> return_fraction_channel{};
  std::array<double, k_stellar_yield_channel_count> total_ejected_metal_fraction_channel{};
  std::array<double, k_stellar_yield_channel_count> newly_synthesized_metal_fraction_channel{};
  std::array<double, k_stellar_yield_channel_count> event_count_per_initial_mass_code_channel{};
  std::array<double, k_stellar_yield_channel_count> energy_erg_per_initial_mass_channel{};

  // Compatibility lane for older callers. It mirrors total_ejected_metal_fraction_total.
  double metal_yield_fraction_total = 0.0;
  std::array<double, k_stellar_yield_channel_count> metal_yield_fraction_channel{};
};

struct StellarEvolutionIntervalBudget {
  double returned_mass_code = 0.0;
  double returned_metals_code = 0.0;
  double newly_synthesized_metals_code = 0.0;
  double event_count = 0.0;
  double feedback_energy_erg = 0.0;
  double remnant_change_code = 0.0;
  std::array<double, k_stellar_yield_channel_count> returned_mass_channel_code{};
  std::array<double, k_stellar_yield_channel_count> returned_metals_channel_code{};
  std::array<double, k_stellar_yield_channel_count> newly_synthesized_metals_channel_code{};
  std::array<double, k_stellar_yield_channel_count> event_count_channel{};
  std::array<double, k_stellar_yield_channel_count> feedback_energy_channel_erg{};
};

struct StellarEvolutionTable {
  std::string table_id = "builtin_zero_yield";
  std::string table_version = "v2";
  std::string source_path = "builtin";
  std::string source_papers;
  std::string source_repository;
  std::string redistribution_license;
  std::string sha256;
  std::string imf;
  std::string stellar_mass_range;
  std::string solar_abundance_reference;
  bool production_calibrated = false;

  // Tensor-product grid. Flattened fields use z-major ordering:
  // flat_index = metallicity_index * age_yr.size() + age_index.
  std::vector<double> age_yr;
  std::vector<double> birth_metallicity_mass_fraction;
  std::vector<double> return_fraction_total;
  std::vector<double> total_ejected_metal_fraction_total;
  std::vector<double> newly_synthesized_metal_fraction_total;
  std::vector<double> event_count_per_initial_mass_code_total;
  std::vector<double> energy_erg_per_initial_mass_code;
  std::array<std::vector<double>, k_stellar_yield_channel_count> return_fraction_channel;
  std::array<std::vector<double>, k_stellar_yield_channel_count> total_ejected_metal_fraction_channel;
  std::array<std::vector<double>, k_stellar_yield_channel_count> newly_synthesized_metal_fraction_channel;
  std::array<std::vector<double>, k_stellar_yield_channel_count> event_count_per_initial_mass_code_channel;
  std::array<std::vector<double>, k_stellar_yield_channel_count> energy_erg_per_initial_mass_channel;

  // Compatibility mirrors accepted by older one-dimensional fixtures. If the
  // v2 total-ejecta vectors are empty, these are interpreted as a single-Z table.
  std::vector<double> metal_yield_fraction_total;
  std::array<std::vector<double>, k_stellar_yield_channel_count> metal_yield_fraction_channel;

  [[nodiscard]] std::size_t flatIndex(std::size_t metallicity_index, std::size_t age_index) const;
  [[nodiscard]] bool isConsistent() const noexcept;
  void requireConsistent() const;
  [[nodiscard]] StellarEvolutionCumulative evaluate(
      double age_yr_value,
      double birth_metallicity_mass_fraction_value) const;
  [[nodiscard]] StellarEvolutionCumulative evaluateAtAgeYears(double age_yr_value) const;
  [[nodiscard]] StellarEvolutionIntervalBudget integrateInterval(
      double age_begin_yr,
      double age_end_yr,
      double initial_mass_code,
      double birth_metallicity_mass_fraction_value) const;
  [[nodiscard]] StellarEvolutionIntervalBudget integrateInterval(
      double age_begin_yr,
      double age_end_yr,
      double initial_mass_code) const;

  [[nodiscard]] static StellarEvolutionTable loadFromTextFile(
      const std::string& path,
      const std::string& source_tag = "stellar_evolution_table");
  [[nodiscard]] static StellarEvolutionTable makeBuiltinReference();
};

struct StellarEvolutionConfig {
  bool enabled = true;
  std::string table_path;
  std::uint32_t metadata_schema_version = 2;
  // Retained only for the source-compatible legacy API. Production runtime
  // passes a physical elapsed-time interval explicitly.
  double hubble_time_years = 1.44e10;
  bool require_production_calibrated_table = false;
};

struct StellarEvolutionStarBudget {
  std::uint32_t star_index = 0;
  std::uint32_t particle_index = 0;
  double star_age_begin_years = 0.0;
  double star_age_end_years = 0.0;
  double mass_old_code = 0.0;
  double mass_new_code = 0.0;
  StellarEvolutionIntervalBudget interval;
};

struct StellarEvolutionStepCounters {
  std::uint64_t scanned_stars = 0;
  std::uint64_t evolved_stars = 0;
  // Stars skipped because they were born at the end of this interval (F5).
  std::uint64_t scanned_newborn_stars = 0;
  double returned_mass_code = 0.0;
  double returned_metals_code = 0.0;
  double newly_synthesized_metals_code = 0.0;
  double event_count = 0.0;
  double feedback_energy_erg = 0.0;
};

struct StellarEvolutionStepReport {
  StellarEvolutionStepCounters counters;
  std::vector<StellarEvolutionStarBudget> budgets;
};

// Gate-4 physical ownership: production fills a caller-owned reusable governed
// budget buffer. The buffer capacity is admitted via the process
// MemoryGovernor BEFORE growth; the callee never grows beyond admitted
// capacity without a new admission. Actual retained capacity() is reconciled
// so the buffer is reused across batches.
struct StellarEvolutionBatchWorkspace {
  std::vector<StellarEvolutionStarBudget> budgets;
  core::MemoryReservation reservation;
  std::uint64_t admitted_bytes = 0U;
};

[[nodiscard]] std::uint64_t stellarEvolutionBatchStagingBytes(std::size_t batch_capacity);

struct StellarEvolutionRuntimeView {
  std::span<const std::uint32_t> active_star_indices;
  std::span<const std::uint32_t> particle_index;
  std::span<const double> birth_mass_code;
  std::span<const double> formation_scale_factor;
  std::span<const double> birth_metallicity_mass_fraction;
  std::span<double> stellar_age_years_last;
  std::span<double> stellar_returned_mass_cumulative_code;
  std::span<double> stellar_returned_metals_cumulative_code;
  std::span<double> stellar_newly_synthesized_metals_cumulative_code;
  std::span<double> stellar_feedback_energy_cumulative_erg;
  std::array<std::span<double>, k_stellar_yield_channel_count> returned_mass_channel_cumulative_code;
  std::array<std::span<double>, k_stellar_yield_channel_count> returned_metals_channel_cumulative_code;
  std::array<std::span<double>, k_stellar_yield_channel_count> feedback_energy_channel_cumulative_erg;
  std::span<double> particle_mass_code;
  // CHUI Part 2 audit F5: birth chronology. Stars are created by the source
  // stage at the END of the accepted interval [t_begin, t_end] (t_birth =
  // t_end). A star whose birth_tick equals newborn_birth_tick was born at
  // t_end of the interval being evolved and therefore has ZERO elapsed age
  // in it: no mass return, yields, winds or SN energy for time before birth.
  // UINT64_MAX disables the rule (legacy diagnostic callers).
  std::span<const std::uint64_t> birth_tick;
  std::uint64_t newborn_birth_tick = std::numeric_limits<std::uint64_t>::max();
};

// Elapsed stellar age (years) of one star inside an accepted interval of
// interval_years, given its birth tick (F5 chronology rule above).
[[nodiscard]] double stellarElapsedYearsInInterval(std::uint64_t birth_tick, std::uint64_t newborn_birth_tick,
                                                   double interval_years) noexcept;

class StellarEvolutionBookkeeper {
 public:
  StellarEvolutionBookkeeper(StellarEvolutionConfig config, StellarEvolutionTable table);

  [[nodiscard]] const StellarEvolutionConfig& config() const noexcept;
  [[nodiscard]] const StellarEvolutionTable& table() const noexcept;

  [[nodiscard]] double evaluateStarAgeYears(
      double formation_scale_factor,
      double current_scale_factor,
      const core::LambdaCdmBackground& background) const;
  // Compatibility overload. It intentionally uses the same persisted age
  // authority and is not called by production SourceRuntime.
  [[nodiscard]] double evaluateStarAgeYears(
      double formation_scale_factor,
      double current_scale_factor) const;

  [[nodiscard]] StellarEvolutionStepReport evaluateElapsedYears(
      const core::SimulationState& state,
      std::span<const std::uint32_t> active_star_indices,
      double elapsed_years) const;
  [[nodiscard]] StellarEvolutionStepReport evaluateElapsedYearsFromView(
      StellarEvolutionRuntimeView view,
      double elapsed_years) const;
  // Gate-4 production path: fills the caller-owned workspace budgets buffer
  // (admitted via governor before growth) and returns only counters. The
  // detailed budgets remain available in workspace.budgets for event
  // construction within the same admitted capacity; no hidden heap vector is
  // grown. Diagnostic callers may keep using evaluateElapsedYears.
  // newborn_birth_tick: F5 chronology (stars with that birth tick were born
  // at the end of this interval and are not evolved). Production passes the
  // current step's birth tick; UINT64_MAX keeps the legacy behaviour.
  [[nodiscard]] StellarEvolutionStepCounters evaluateElapsedYearsGoverned(
      const core::SimulationState& state,
      std::span<const std::uint32_t> active_star_indices,
      double elapsed_years,
      core::MemoryGovernor* governor,
      StellarEvolutionBatchWorkspace& workspace,
      std::uint64_t newborn_birth_tick = std::numeric_limits<std::uint64_t>::max()) const;
  void commitBudgets(
      core::SimulationState& state,
      const StellarEvolutionStepReport& report) const;
  void commitBudgetsFromView(
      StellarEvolutionRuntimeView view,
      const StellarEvolutionStepReport& report) const;

  [[nodiscard]] StellarEvolutionStepReport applyElapsedYears(
      core::SimulationState& state,
      std::span<const std::uint32_t> active_star_indices,
      double elapsed_years) const;

  // Source-compatible legacy entry points. They no longer derive age from
  // scale factor; dt_code is mapped to elapsed years only for isolated tests.
  [[nodiscard]] StellarEvolutionStepReport apply(
      core::SimulationState& state,
      std::span<const std::uint32_t> active_star_indices,
      double current_scale_factor,
      double dt_code) const;
  [[nodiscard]] StellarEvolutionStepReport applyFromView(
      StellarEvolutionRuntimeView view,
      double current_scale_factor,
      double dt_code) const;

  [[nodiscard]] core::ModuleSidecarBlock buildMetadataSidecar(
      const StellarEvolutionStepReport& report) const;

 private:
  [[nodiscard]] StellarEvolutionRuntimeView makeRuntimeView(
      core::SimulationState& state,
      std::span<const std::uint32_t> active_star_indices) const;

  StellarEvolutionConfig m_config;
  StellarEvolutionTable m_table;
};

[[nodiscard]] StellarEvolutionConfig makeStellarEvolutionConfig(
    const core::PhysicsConfig& physics_config);
[[nodiscard]] StellarEvolutionTable loadStellarEvolutionTable(
    const core::PhysicsConfig& physics_config);

}  // namespace cosmosim::physics
