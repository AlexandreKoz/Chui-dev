#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <memory>
#include <memory_resource>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"

namespace cosmosim::physics {

inline constexpr std::uint32_t kStarFormationModelSchemaVersion = 3;
inline constexpr std::uint32_t kStarFormationRngKeySchemaVersion = 1;

struct StarFormationConfig {
  bool enabled = true;
  core::StarFormationModelKind model = core::StarFormationModelKind::kLegacySchmidtThreshold;

  // Legacy Schmidt-threshold compatibility model.
  double density_threshold_code = 10.0;
  double temperature_threshold_k = 1.0e4;
  double min_converging_flow_rate_code = 0.0;

  // Shared and adaptive_bound_jeans parameters.
  double epsilon_ff = 0.01;
  double bound_alpha_vir_max = 1.0;
  bool require_converging_flow = true;
  core::StarFormationCollapseTimescale collapse_timescale =
      core::StarFormationCollapseTimescale::kFreeFall;
  double jeans_mass_floor_code = 100.0;
  core::StarParticleMassPolicy star_particle_mass_policy =
      core::StarParticleMassPolicy::kFixed;
  double target_star_particle_mass_code = 0.1;
  double target_star_particle_mass_fraction = 0.25;
  double min_star_particle_mass_code = 0.1;
  double max_star_particle_mass_code = 1.0e30;
  std::uint32_t max_spawn_particles_per_cell_step = 8;
  double max_fractional_mass_conversion = 0.25;
  double min_remaining_gas_fraction = 0.01;
  double min_remaining_gas_mass_code = 1.0e-12;
  double temperature_safety_ceiling_k = 0.0;
  bool stochastic_spawning = true;
  std::uint64_t random_seed = 123456789ull;

  double effective_min_baryon_overdensity = 0.0;
  double effective_hot_excess_tolerance = 1.0;
  double effective_massive_star_fraction = 0.1;
  core::EffectiveIsmBirthMassConvention effective_birth_mass_convention =
      core::EffectiveIsmBirthMassConvention::kInitialSspMass;

  double newton_g_code = core::constants::k_newton_g_si;
  bool density_is_comoving = false;
  bool geometry_is_comoving = true;
  std::uint32_t metadata_schema_version = kStarFormationModelSchemaVersion;
};

enum class StarFormationRejectionReason : std::uint8_t {
  kNone,
  kInactive,
  kNotOwned,
  kNonLeaf,
  kGhost,
  kInvalidState,
  kNonPositiveMass,
  kNonPositiveVolume,
  kNotConverging,
  kUnbound,
  kJeansStable,
  kLegacyDensity,
  kLegacyTemperature,
  kAdaptiveMassFloor,
  kEffectiveBelowDensityThreshold,
  kEffectiveBelowOverdensityThreshold,
  kEffectiveHotAboveEos,
  kEffectiveInvalidEquilibrium,
};

struct StarFormationCellInput {
  // local_cell_row is transient and is used only to apply a validated plan.
  std::uint32_t cell_index = 0;
  std::uint64_t gas_cell_id = 0;
  std::uint32_t owning_rank = 0;
  bool is_active = true;
  bool is_owned = true;
  bool is_leaf = true;
  bool is_ghost = false;

  double gas_mass_code = 0.0;
  double gas_density_code = 0.0;
  double cell_volume_code = 0.0;
  double gas_temperature_k = 0.0;
  double gas_sound_speed_code = 0.0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double velocity_divergence_code = 0.0;
  double velocity_gradient_frobenius_sq_code = 0.0;
  double gas_metal_mass_code = 0.0;
  double metallicity_mass_fraction = 0.0;  // legacy/injected compatibility lane
  double gas_specific_internal_energy_code = 0.0;
  double baryon_overdensity = 0.0;
  bool is_cosmological = false;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
};

struct StarFormationRuntimeView {
  std::span<const std::uint32_t> active_cell_indices;
  std::span<const double> center_x_comoving;
  std::span<const double> center_y_comoving;
  std::span<const double> center_z_comoving;
  std::span<double> gas_mass_code;
  std::span<double> gas_density_code;
  std::span<const double> gas_temperature_k;
  std::span<const double> velocity_divergence_code;
  std::span<const double> metallicity_mass_fraction;
};

struct StarFormationCounters {
  std::uint64_t scanned_cells = 0;
  std::uint64_t rejected_inactive = 0;
  std::uint64_t rejected_not_owned = 0;
  std::uint64_t rejected_non_leaf = 0;
  std::uint64_t rejected_ghost = 0;
  std::uint64_t rejected_invalid_state = 0;
  std::uint64_t rejected_non_positive_mass = 0;
  std::uint64_t rejected_non_positive_volume = 0;
  std::uint64_t rejected_not_converging = 0;
  std::uint64_t rejected_unbound = 0;
  std::uint64_t rejected_jeans_stable = 0;
  std::uint64_t rejected_legacy_density_threshold = 0;
  std::uint64_t rejected_legacy_temperature_threshold = 0;
  std::uint64_t rejected_adaptive_mass_floor = 0;
  std::uint64_t rejected_effective_below_density_threshold = 0;
  std::uint64_t rejected_effective_below_overdensity_threshold = 0;
  std::uint64_t rejected_effective_hot_above_eos = 0;
  std::uint64_t rejected_effective_invalid_equilibrium = 0;
  std::uint64_t eligible_cells = 0;
  std::uint64_t effective_cells_scanned = 0;
  std::uint64_t effective_cells_above_threshold = 0;
  std::uint64_t effective_cells_on_eos = 0;
  std::uint64_t effective_cells_hot_above_eos = 0;
  std::uint64_t effective_cells_invalid_equilibrium = 0;
  double effective_cold_gas_mass_code = 0.0;
  double effective_instantaneous_sfr_code = 0.0;
  double effective_min_cold_fraction = 1.0;
  double effective_max_cold_fraction = 0.0;
  double effective_min_pressure_code = 0.0;
  double effective_max_pressure_code = 0.0;
  std::uint64_t spawn_events = 0;
  std::uint64_t spawned_particles = 0;
  double expected_spawn_mass_code = 0.0;
  double spawned_mass_code = 0.0;
  double minimum_free_fall_time_code = 0.0;
  double minimum_collapse_time_code = 0.0;
  double maximum_spawn_probability = 0.0;
  double maximum_fractional_mass_conversion = 0.0;

  double gas_mass_removed_code = 0.0;
  double star_mass_created_code = 0.0;
  double mass_residual_code = 0.0;
  double gas_momentum_removed_x_code = 0.0;
  double gas_momentum_removed_y_code = 0.0;
  double gas_momentum_removed_z_code = 0.0;
  double star_momentum_created_x_code = 0.0;
  double star_momentum_created_y_code = 0.0;
  double star_momentum_created_z_code = 0.0;
  double momentum_residual_norm_code = 0.0;
  double gas_metal_mass_removed_code = 0.0;
  double star_metal_mass_created_code = 0.0;
  double metal_mass_residual_code = 0.0;
  double gas_internal_energy_removed_code = 0.0;
  double star_kinetic_energy_created_code = 0.0;
  double star_formation_internal_energy_sink_code = 0.0;
};

struct StarFormationCellOutcome {
  bool eligible = false;
  StarFormationRejectionReason rejection_reason = StarFormationRejectionReason::kNone;
  double physical_density_code = 0.0;
  double physical_cell_scale_code = 0.0;
  double virial_parameter = 0.0;
  double jeans_mass_code = 0.0;
  double free_fall_time_code = 0.0;
  double compression_time_code = 0.0;
  double collapse_time_code = 0.0;
  double cold_mass_fraction = 0.0;
  double effective_specific_internal_energy_code = 0.0;
  double effective_pressure_code = 0.0;
  double effective_signal_speed_squared_code = 0.0;
  double sfr_density_rate_code = 0.0;
  double expected_spawn_mass_code = 0.0;
  std::uint32_t spawned_particle_count = 0;
  double spawned_mass_code = 0.0;
  double target_particle_mass_code = 0.0;
  double spawn_probability = 0.0;
  double random_u01 = 0.0;
};

struct StarFormationStepReport {
  StarFormationCounters counters;
  std::vector<std::uint32_t> spawned_from_cells;
  std::vector<std::uint64_t> birth_keys;
};

// Exact internal birth-plan element width for owner-local staging admission.
[[nodiscard]] std::size_t starFormationBirthPlanBytes() noexcept;

[[nodiscard]] std::uint64_t starFormationBirthKey(
    std::uint64_t gas_cell_id,
    std::uint64_t global_integration_tick,
    std::uint32_t birth_attempt_ordinal,
    std::uint32_t model_schema_version = kStarFormationModelSchemaVersion);

[[nodiscard]] double starFormationUniform01(
    std::uint64_t global_seed,
    std::uint64_t gas_cell_id,
    std::uint64_t global_integration_tick,
    std::uint32_t birth_attempt_ordinal,
    std::uint32_t rng_key_schema_version = kStarFormationRngKeySchemaVersion);

[[nodiscard]] std::uint64_t starFormationParticleIdFromBirthKey(
    std::uint64_t birth_key,
    std::uint32_t collision_ordinal = 0U);

// Exact, allocation-bounded batch precommit used by local and distributed registries.
// Existing IDs and immutable birth keys are validated before any state mutation.
[[nodiscard]] std::vector<std::uint64_t> precommitStarParticleIdsExact(
    std::span<const std::uint64_t> existing_particle_ids,
    std::span<const std::uint64_t> birth_keys);

class ParticleIdPrecommit {
 public:
  virtual ~ParticleIdPrecommit() = default;
  // Releases a result lease after the caller finishes its birth transaction.
  // Legacy/standalone registries have no retained result lease.
  virtual void finishPrecommit() noexcept {}
  virtual void preparePopulationGrowth(core::SimulationState&, std::size_t,
                                       core::ParticleSpecies) {}
  [[nodiscard]] virtual std::vector<std::uint64_t> precommit(
      const core::SimulationState& state,
      std::span<const std::uint64_t> birth_keys) = 0;
};

class LocalParticleIdRegistry final : public ParticleIdPrecommit {
 public:
  [[nodiscard]] std::vector<std::uint64_t> precommit(
      const core::SimulationState& state,
      std::span<const std::uint64_t> birth_keys) override;
};

class StarFormationModel {
 public:
  explicit StarFormationModel(
      StarFormationConfig config,
      std::shared_ptr<const EffectiveMultiphaseEosTable> effective_eos_table = nullptr);

  [[nodiscard]] const StarFormationConfig& config() const noexcept;
  [[nodiscard]] bool isEligible(const StarFormationCellInput& cell) const;
  [[nodiscard]] double physicalDensityCode(double stored_density_code, double scale_factor) const;
  [[nodiscard]] double physicalCellScaleCode(double stored_volume_code, double scale_factor) const;
  [[nodiscard]] double virialParameter(
      double physical_density_code,
      double physical_cell_scale_code,
      double sound_speed_code,
      double velocity_gradient_frobenius_sq_code) const;
  [[nodiscard]] double jeansMassCode(double physical_density_code, double sound_speed_code) const;
  [[nodiscard]] double freeFallTimeCode(double gas_density_code) const;
  [[nodiscard]] double compressionTimeCode(double velocity_divergence_code) const;
  [[nodiscard]] double sourceTimeStepLimitCode(double collapse_time_code) const;
  [[nodiscard]] double sfrDensityRateCode(double gas_density_code) const;
  [[nodiscard]] double expectedSpawnMassCode(const StarFormationCellInput& cell, double dt_code) const;
  [[nodiscard]] StarFormationCellOutcome sampleCellOutcome(
      const StarFormationCellInput& cell,
      double dt_code,
      std::uint64_t step_index,
      std::uint32_t ignored_rank_local_seed_offset = 0,
      double scale_factor = 1.0) const;

  [[nodiscard]] StarFormationStepReport applyFromInputs(
      core::SimulationState& state,
      std::span<const StarFormationCellInput> cell_inputs,
      double dt_code,
      double scale_factor,
      std::uint64_t global_integration_tick,
      ParticleIdPrecommit* id_precommit = nullptr,
      std::pmr::memory_resource* metadata_scratch = nullptr) const;

  [[nodiscard]] StarFormationStepReport applyFromView(
      core::SimulationState& state,
      StarFormationRuntimeView view,
      double dt_code,
      double scale_factor,
      std::uint64_t step_index,
      std::uint32_t owning_rank = 0) const;

  [[nodiscard]] StarFormationStepReport apply(
      core::SimulationState& state,
      std::span<const std::uint32_t> active_cell_indices,
      std::span<const double> velocity_divergence_code,
      std::span<const double> metallicity_mass_fraction,
      double dt_code,
      double scale_factor,
      std::uint64_t step_index,
      std::uint32_t owning_rank = 0) const;

  [[nodiscard]] core::ModuleSidecarBlock buildMetadataSidecar(
      const StarFormationCounters& counters,
      std::string_view configuration_hash = {},
      const StarFormationCounters* cumulative_counters = nullptr) const;

  // CHUI Part 2 birth chronology: code time of the END of the accepted
  // interval in which the source stage creates stars (t_birth = t_end).
  // Recorded in the stellar population record (F5).
  void setBirthTimeCode(double time_end_code) noexcept { m_birth_time_code = time_end_code; }

 private:
  [[nodiscard]] StarFormationCellOutcome evaluateCell(
      const StarFormationCellInput& cell,
      double dt_code,
      double scale_factor,
      std::uint64_t global_integration_tick) const;

  StarFormationConfig m_config;
  std::shared_ptr<const EffectiveMultiphaseEosTable> m_effective_eos_table;
  double m_birth_time_code = 0.0;
};

[[nodiscard]] StarFormationConfig makeStarFormationConfig(const core::PhysicsConfig& physics_config);

class StarFormationCallback final : public core::IntegrationCallback {
 public:
  explicit StarFormationCallback(StarFormationModel model, std::uint32_t owning_rank = 0);

  [[nodiscard]] std::string_view callbackName() const override;
  [[nodiscard]] std::span<const core::IntegrationStage> integrationStages() const override;
  [[nodiscard]] std::span<const core::StageContract> stageContracts() const override;
  void onStage(core::StepContext& context) override;

  // Test/compatibility injection seam. Production SourceRuntime builds these values
  // from authoritative hydro and gas-metal state instead.
  void setVelocityDivergenceCode(std::span<const double> velocity_divergence_code);
  void setMetallicityMassFraction(std::span<const double> metallicity_mass_fraction);
  void setRankLocalSeedOffset(std::uint32_t owning_rank);

  [[nodiscard]] const StarFormationStepReport& lastStepReport() const noexcept;

 private:
  void ensureFieldSizes(std::size_t cell_count);

  StarFormationModel m_model;
  std::uint32_t m_owning_rank = 0;
  std::vector<std::uint32_t> m_full_cell_indices;
  std::vector<double> m_velocity_divergence_code;
  std::vector<double> m_metallicity_mass_fraction;
  StarFormationStepReport m_last_step_report;
};

}  // namespace cosmosim::physics
