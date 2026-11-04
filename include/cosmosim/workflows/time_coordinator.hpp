#pragma once

#include <cstdint>
#include <optional>
#include <memory>
#include <span>
#include <vector>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/output_restart_runtime.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::core {
struct UnitSystem;
}

namespace cosmosim::io {
struct RestartReadResult;
}

namespace cosmosim::physics {
class EffectiveMultiphaseEosTable;
}

namespace cosmosim::workflows {

class GravityRuntime;
class HydroAmrRuntime;
class SourceRuntime;
struct RuntimeServices;
namespace internal {
class MigrationBalanceRuntime;
struct PendingOutputBoundary;
}

// Owns the two scheduler authorities, integrator truth, and ordered output
// cadence for one rung-zero runtime. Construction imports restart state or
// initializes fresh state through the single factory below.
class RungZeroTimeState {
 public:
  [[nodiscard]] core::HierarchicalTimeBinScheduler& particleScheduler() noexcept;
  [[nodiscard]] const core::HierarchicalTimeBinScheduler& particleScheduler() const noexcept;
  [[nodiscard]] core::HierarchicalTimeBinScheduler& gasCellScheduler() noexcept;
  [[nodiscard]] const core::HierarchicalTimeBinScheduler& gasCellScheduler() const noexcept;
  [[nodiscard]] core::IntegratorState& integratorState() noexcept;
  [[nodiscard]] const core::IntegratorState& integratorState() const noexcept;
  [[nodiscard]] internal::PendingOutputBoundary& pendingOutput() noexcept;
  [[nodiscard]] const internal::PendingOutputBoundary& pendingOutput() const noexcept;

 private:
  friend class TimeCoordinator;
  friend RungZeroTimeState initializeRungZeroTimeState(
      const core::SimulationConfig&,
      const ReferenceWorkflowOptions&,
      core::SimulationState&,
      const core::UnitSystem&,
      const core::LambdaCdmBackground*,
      const io::RestartReadResult*);

  explicit RungZeroTimeState(std::uint8_t max_bin);

  core::HierarchicalTimeBinScheduler m_particle_scheduler;
  core::HierarchicalTimeBinScheduler m_gas_cell_scheduler;
  core::IntegratorState m_integrator_state;
  internal::PendingOutputBoundary m_pending_output;
};

[[nodiscard]] RungZeroTimeState initializeRungZeroTimeState(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    core::SimulationState& state,
    const core::UnitSystem& units,
    const core::LambdaCdmBackground* cosmology_background,
    const io::RestartReadResult* restart_state);

// Rung-zero KDK lifecycle and typed stage-dispatch owner. The lower-level core
// orchestrator retains timeline/PM invariants but receives this coordinator's
// typed dispatcher; no production owner receives StepContext through its public
// interface.
// G1.2: the criterion that limited this rank's timestep candidate at the
// last adaptive-bin refresh (the rank whose candidate equals the global
// minimum reports the step's actual limiter).
struct TimestepLimiterRecord {
  double dt_time_code = 0.0;
  const char* criterion = "none";
  bool gas_cell = false;
  std::uint32_t element = 0U;
  int amr_level = -1;
  double hydro_cfl_level_relaxation = 1.0;
};

class TimeCoordinator {
 public:
  TimeCoordinator(
      const RuntimeServices& services,
      RungZeroTimeState& time_state,
      GravityRuntime& gravity,
      HydroAmrRuntime& hydro_amr,
      SourceRuntime& source,
      RuntimeExecutionPlan execution_plan,
      const internal::MigrationBalanceRuntime& migration_balance,
      std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table = nullptr) noexcept;

  // Execute one bounded rung-zero segment. This owns active-set formation,
  // ordered endpoint/output clipping, adaptive-bin refresh, rebalance, and
  // output-boundary dispatch while preserving the existing KDK stage order.
  void runRungZeroSegment(
      const core::SimulationConfig& config,
      const ReferenceWorkflowOptions& options,
      core::SimulationState& state,
      const core::LambdaCdmBackground* cosmology_background,
      std::vector<std::uint64_t>& expected_global_particle_ids,
      ReferenceWorkflowReport& report,
      core::ProfilerSession& profiler,
      const core::ModePolicy& mode_policy,
      bool restoring_from_restart);

  void executeSingleStep(
      core::SimulationState& state,
      core::IntegratorState& integrator_state,
      core::ActiveSetDescriptor active_set,
      const core::LambdaCdmBackground* cosmology_background,
      core::TransientStepWorkspace* workspace,
      const core::ModePolicy* mode_policy,
      core::ProfilerSession* profiler_session,
      std::optional<std::uint64_t> expected_scheduler_tick,
      core::StepBoundaryKind requested_boundary_kind);

  void executeOutputBoundary(
      core::SimulationState& state,
      core::IntegratorState& integrator_state,
      core::ProfilerSession* profiler_session,
      core::StepBoundaryKind requested_boundary_kind);

  [[nodiscard]] double updateAdaptiveTimeBins(
      core::SimulationState& state,
      core::HierarchicalTimeBinScheduler& particle_scheduler,
      core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
      const core::IntegratorState& integrator_state,
      const core::SimulationConfig& config,
      std::span<const double> particle_accel_x,
      std::span<const double> particle_accel_y,
      std::span<const double> particle_accel_z,
      std::span<const double> cell_accel_x,
      std::span<const double> cell_accel_y,
      std::span<const double> cell_accel_z,
      std::span<const std::uint32_t> active_particle_indices,
      std::span<const std::uint32_t> active_cell_indices,
      bool update_all_elements);
  [[nodiscard]] const TimestepLimiterRecord& lastTimestepLimiter() const noexcept { return m_last_timestep_limiter; }

 private:
  void runHierarchicalSegment(
      const core::SimulationConfig&, const ReferenceWorkflowOptions&,
      core::SimulationState&, const core::LambdaCdmBackground*,
      std::vector<std::uint64_t>&, ReferenceWorkflowReport&, core::ProfilerSession&,
      const core::ModePolicy&);
  void coordinatePmRefreshDirective(core::StepContext& context);
  void dispatchStage(
      core::StepContext& context,
      bool require_output_safe_boundary);

  const RuntimeServices& m_services;
  RungZeroTimeState& m_time_state;
  GravityRuntime& m_gravity;
  HydroAmrRuntime& m_hydro_amr;
  SourceRuntime& m_source;
  RuntimeExecutionPlan m_execution_plan;
  const internal::MigrationBalanceRuntime& m_migration_balance;
  // Cleared for each KDK step and populated only by legal source-stage births.
  // Keeping this distinct from the ownership ledger prevents migrations from
  // being mistaken for particle creation.
  std::vector<std::uint64_t> m_newly_created_particle_ids;
  // Immutable configuration-derived EOS table; build once per coordinator,
  // not once per adaptive-bin refresh.
  std::shared_ptr<const physics::EffectiveMultiphaseEosTable> m_effective_eos_table;
  core::StepOrchestrator m_lifecycle;
  TimestepLimiterRecord m_last_timestep_limiter{};
};

}  // namespace cosmosim::workflows
