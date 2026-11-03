#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

enum class StellarFeedbackMode : std::uint8_t {
  kThermal = 0,
  kKinetic = 1,
  kMomentum = 2,
  kThermalKineticMomentum = 3,
};

enum class StellarFeedbackVariant : std::uint8_t {
  kNone = 0,
  kDelayedCooling = 1,
  kStochastic = 2,
};

struct StellarFeedbackConfig {
  bool enabled = true;
  StellarFeedbackMode mode = StellarFeedbackMode::kThermalKineticMomentum;
  StellarFeedbackVariant variant = StellarFeedbackVariant::kNone;
  bool use_returned_mass_budget = true;
  double epsilon_thermal = 0.6;
  double epsilon_kinetic = 0.3;
  double epsilon_momentum = 0.1;
  double sn_energy_erg_per_mass_code = 1.0e49;
  double momentum_code_per_mass_code = 3.0e3;
  // Converts physical total energy in erg to code total energy (M_code V_code^2).
  double total_energy_code_per_erg = 1.0;
  std::uint32_t neighbor_count = 8;
  double delayed_cooling_time_code = 0.0;
  double stochastic_event_probability = 0.25;
  std::uint64_t random_seed = 42424242ull;
  std::uint32_t metadata_schema_version = 2;
};

struct StellarFeedbackGeometryView {
  std::span<const double> particle_position_x_comoving;
  std::span<const double> particle_position_y_comoving;
  std::span<const double> particle_position_z_comoving;
  std::span<const double> cell_center_x_comoving;
  std::span<const double> cell_center_y_comoving;
  std::span<const double> cell_center_z_comoving;
  // Optional filters. Empty spans retain legacy local-state behavior.
  std::span<const std::uint64_t> gas_cell_id;
  std::span<const std::uint8_t> is_owned_leaf;
  std::span<const std::uint32_t> candidate_cell_indices;
};

// Phase-resident, locality-preserving index for exact nearest-cell feedback
// queries.  The index stores only cell rows sorted by x coordinate; per-event
// target scratch remains O(neighbor_count), not O(N_gas).  It is rebuilt from
// current geometry and never becomes restart authority.
class StellarFeedbackSpatialIndex {
 public:
  void rebuild(const StellarFeedbackGeometryView& geometry_view);
  void rebuildFromCellIndices(
      const StellarFeedbackGeometryView& geometry_view,
      std::vector<std::uint32_t> cell_indices);
  void clear() noexcept;

  [[nodiscard]] std::span<const std::uint32_t> sortedCellIndices() const noexcept;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;
  [[nodiscard]] std::uint64_t highWaterBytes() const noexcept;

 private:
  std::vector<std::uint32_t> m_sorted_cell_indices;
  std::uint64_t m_high_water_bytes = 0U;
};

class StellarFeedbackCellVolumeProvider {
 public:
  virtual ~StellarFeedbackCellVolumeProvider() = default;
  [[nodiscard]] virtual double cellVolumeCode(std::uint32_t cell_index) const = 0;
};

// Receives every gas cell whose mass, density, metals or internal energy is
// about to be changed by feedback deposition (audit F4). Implementations must
// not throw from capacity growth they did not admit beforehand; the source
// runtime's implementation reserves through the MemoryGovernor.
class StellarFeedbackThermodynamicTouchSink {
 public:
  virtual ~StellarFeedbackThermodynamicTouchSink() = default;
  virtual void recordGasCellMutation(std::uint32_t cell_index, double internal_energy_before_code) = 0;
};

struct StellarFeedbackDepositionView {
  std::span<double> cell_mass_code;
  std::span<double> gas_density_code;
  std::span<double> gas_internal_energy_code;
  std::span<double> gas_metal_mass_code;
  // Optional. If absent, volume is derived once from old mass/old density.
  std::span<const double> cell_volume_code;
  // Optional phase-local provider used by AMR workflows to derive exact cell
  // geometry without retaining a full-population volume vector.
  const StellarFeedbackCellVolumeProvider* cell_volume_provider = nullptr;
  // Optional (audit F4): notified before each gas-cell thermodynamic mutation.
  StellarFeedbackThermodynamicTouchSink* thermodynamic_touch_sink = nullptr;
};

struct StellarFeedbackEvent {
  std::uint32_t star_index = 0;
  double returned_mass_code = 0.0;
  double returned_metals_code = 0.0;
  double feedback_energy_erg = 0.0;
};

struct StellarFeedbackBudget {
  double source_mass_code = 0.0;
  double returned_mass_code = 0.0;
  double returned_metals_code = 0.0;
  double total_energy_erg = 0.0;
  double thermal_energy_erg = 0.0;
  double kinetic_energy_erg = 0.0;
  double momentum_budget_code = 0.0;
};

struct StellarFeedbackTarget {
  std::uint32_t cell_index = 0;
  double weight = 0.0;
  double radial_dx_comoving = 0.0;
  double radial_dy_comoving = 0.0;
  double radial_dz_comoving = 0.0;
};

struct StellarFeedbackStarReport {
  std::uint32_t star_index = 0;
  std::uint32_t particle_index = 0;
  StellarFeedbackBudget budget;
  std::size_t target_count = 0;
  bool stochastic_event_fired = true;
  bool delayed_cooling_applied = false;
  double deposited_mass_code = 0.0;
  double deposited_metals_code = 0.0;
  double deposited_thermal_energy_erg = 0.0;
  double deposited_kinetic_energy_erg = 0.0;
  double deposited_momentum_code = 0.0;
  double unresolved_mass_code = 0.0;
  double unresolved_metals_code = 0.0;
  double unresolved_thermal_energy_erg = 0.0;
  double unresolved_kinetic_energy_erg = 0.0;
  double unresolved_momentum_code = 0.0;
};

struct StellarFeedbackStepCounters {
  std::uint64_t scanned_stars = 0;
  std::uint64_t feedback_stars = 0;
  std::uint64_t target_cells_visited = 0;
  double source_mass_code = 0.0;
  double deposited_mass_code = 0.0;
  double deposited_metals_code = 0.0;
  double deposited_thermal_energy_erg = 0.0;
  double deposited_kinetic_energy_erg = 0.0;
  double deposited_momentum_code = 0.0;
  double unresolved_mass_code = 0.0;
  double unresolved_metals_code = 0.0;
  double unresolved_thermal_energy_erg = 0.0;
  double unresolved_kinetic_energy_erg = 0.0;
  double unresolved_momentum_code = 0.0;
};

struct StellarFeedbackStepReport {
  StellarFeedbackStepCounters counters;
  std::vector<StellarFeedbackStarReport> star_reports;
};

// Compatibility/runtime token only. Persistent carry authority lives in
// StarParticleSidecar. M2C intentionally keeps no population-scale duplicate
// carry vectors here: unresolved budgets must survive restart/migration in one
// authority, not in a second runtime mirror.
struct StellarFeedbackModuleState {
  void ensureStarCapacity(std::size_t star_count) noexcept;
  [[nodiscard]] constexpr std::uint64_t ownedCapacityBytes() const noexcept {
    return 0U;
  }
};

class StellarFeedbackModel {
 public:
  explicit StellarFeedbackModel(StellarFeedbackConfig config);

  [[nodiscard]] const StellarFeedbackConfig& config() const noexcept;
  [[nodiscard]] StellarFeedbackBudget computeBudget(
      double source_mass_code,
      double returned_mass_code,
      double returned_metals_code) const;
  [[nodiscard]] StellarFeedbackBudget computeBudgetFromEnergy(
      double source_mass_code,
      double returned_mass_code,
      double returned_metals_code,
      double feedback_energy_erg) const;

  [[nodiscard]] std::vector<StellarFeedbackTarget> selectTargets(
      const StellarFeedbackGeometryView& geometry_view,
      std::uint32_t particle_index) const;
  [[nodiscard]] std::vector<StellarFeedbackTarget> selectTargets(
      const StellarFeedbackGeometryView& geometry_view,
      const StellarFeedbackSpatialIndex& spatial_index,
      std::uint32_t particle_index) const;
  [[nodiscard]] std::vector<StellarFeedbackTarget> selectTargets(
      const core::SimulationState& state,
      std::uint32_t particle_index) const;

  [[nodiscard]] StellarFeedbackStepReport applyWithViews(
      core::SimulationState& state,
      StellarFeedbackModuleState& module_state,
      const StellarFeedbackGeometryView& geometry_view,
      StellarFeedbackDepositionView deposition_view,
      std::span<const std::uint32_t> active_star_indices,
      std::span<const double> returned_mass_delta_code,
      std::span<const double> returned_metals_delta_code,
      double dt_code,
      std::span<const double> feedback_energy_delta_erg = {}) const;

  [[nodiscard]] StellarFeedbackStepReport applyEventsWithViews(
      core::SimulationState& state,
      StellarFeedbackModuleState& module_state,
      const StellarFeedbackGeometryView& geometry_view,
      const StellarFeedbackSpatialIndex* spatial_index,
      StellarFeedbackDepositionView deposition_view,
      std::span<const StellarFeedbackEvent> events,
      double dt_code) const;

  // Production path: performs identical deposition and accumulates counters
  // without materializing a per-star report vector. The caller discards the
  // star_reports anyway, so this avoids a population-scale allocation.
  void applyEventsCountersOnly(
      core::SimulationState& state,
      StellarFeedbackModuleState& module_state,
      const StellarFeedbackGeometryView& geometry_view,
      const StellarFeedbackSpatialIndex* spatial_index,
      StellarFeedbackDepositionView deposition_view,
      std::span<const StellarFeedbackEvent> events,
      double dt_code,
      StellarFeedbackStepReport& report_out) const;

  [[nodiscard]] StellarFeedbackStepReport apply(
      core::SimulationState& state,
      StellarFeedbackModuleState& module_state,
      std::span<const std::uint32_t> active_star_indices,
      std::span<const double> returned_mass_delta_code,
      std::span<const double> returned_metals_delta_code,
      double dt_code,
      std::span<const double> feedback_energy_delta_erg = {}) const;

  [[nodiscard]] core::ModuleSidecarBlock buildMetadataSidecar(
      const StellarFeedbackStepReport& report) const;

 private:
  [[nodiscard]] static std::string modeToString(StellarFeedbackMode mode);
  [[nodiscard]] static std::string variantToString(StellarFeedbackVariant variant);
  [[nodiscard]] bool stochasticEventFires(
      std::uint32_t star_index, std::uint64_t step_seed) const;

  StellarFeedbackConfig m_config;
};

[[nodiscard]] StellarFeedbackConfig makeStellarFeedbackConfig(
    const core::PhysicsConfig& physics_config);

}  // namespace cosmosim::physics
