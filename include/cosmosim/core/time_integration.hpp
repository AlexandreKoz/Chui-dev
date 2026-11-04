#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::core {

class HierarchicalTimeBinScheduler;
struct IntegratorState;

// Explicit step stage contract shared by gravity, hydro, physics, analysis, and I/O modules.
enum class IntegrationStage : std::uint8_t {
  kGravityKickPre = 0,
  kDrift = 1,
  kForceRefresh = 2,
  kHydroUpdate = 3,
  kSourceTerms = 4,
  kGravityKickPost = 5,
  kAnalysisHooks = 6,
  kOutputCheck = 7,
};
enum class StageDataDomain : std::uint16_t {
  kNone = 0,
  kParticles = 1U << 0U,
  kGasCells = 1U << 1U,
  kMeshCells = 1U << 2U,
  kGhostCells = 1U << 3U,
  kPmField = 1U << 4U,
  kTreeState = 1U << 5U,
  kOutputState = 1U << 6U,
  kRestartState = 1U << 7U,
  kDiagnostics = 1U << 8U,
};
[[nodiscard]] constexpr StageDataDomain operator|(StageDataDomain lhs, StageDataDomain rhs) noexcept {
  return static_cast<StageDataDomain>(static_cast<std::uint16_t>(lhs) | static_cast<std::uint16_t>(rhs));
}
[[nodiscard]] constexpr bool hasStageDataDomain(StageDataDomain value, StageDataDomain mask) noexcept {
  return (static_cast<std::uint16_t>(value) & static_cast<std::uint16_t>(mask)) != 0U;
}
enum class StageSyncRequirement : std::uint8_t {
  kNone = 0,
  kLocalOnly = 1,
  kPmRefreshBoundary = 2,
  kGlobal = 3,
  // Force/acceleration evaluation may run on local active sets.  A PM long-range
  // refresh is a stricter sub-event and is authorized separately through
  // PmRefreshDirective only at legal global PM-refresh boundaries.
  kForceEvaluation = 4,
};
enum class StageActiveSetFamily : std::uint8_t {
  kNone = 0,
  kAllParticles = 1,
  kActiveParticles = 2,
  kGasCells = 3,
  kMeshCells = 4,
  kGhostCells = 5,
  kOutputState = 6,
  kRestartState = 7,
};
enum class StageSafety : std::uint8_t { kUnsafe = 0, kSafe = 1 };
enum class StageSubsystem : std::uint8_t { kCore = 0, kGravity = 1, kHydro = 2, kSources = 3, kAnalysis = 4, kOutput = 5 };
struct StageContract {
  IntegrationStage stage = IntegrationStage::kGravityKickPre;
  StageDataDomain required_inputs = StageDataDomain::kNone;
  StageDataDomain mutated_state = StageDataDomain::kNone;
  StageDataDomain produced_outputs = StageDataDomain::kNone;
  StageDataDomain allowed_side_effects = StageDataDomain::kNone;
  StageSyncRequirement sync_requirements = StageSyncRequirement::kNone;
  StageActiveSetFamily active_set_family = StageActiveSetFamily::kNone;
  StageSafety restart_safety = StageSafety::kUnsafe;
  StageSafety output_safety = StageSafety::kUnsafe;
  StageSubsystem owner = StageSubsystem::kCore;
};

[[nodiscard]] std::string_view integrationStageName(IntegrationStage stage);
[[nodiscard]] constexpr std::size_t integrationStageCount() noexcept { return 8; }
[[nodiscard]] constexpr std::size_t integrationStageIndex(IntegrationStage stage) noexcept {
  return static_cast<std::size_t>(stage);
}

// Baseline stepping family; hierarchical bins can reuse the same stage contract.
enum class TimeStepScheme : std::uint8_t {
  kKickDriftKick = 0,
};

// Deterministic source coupling around the hydro update.  Discrete source
// events (star birth, stellar returns/feedback and BH transactions) execute
// only in kStrangPostHalf; smooth deterministic operators execute in both
// half steps.  kUnsplitCompatibility is retained for direct low-level callers
// that dispatch a source stage outside the canonical scheduler.
enum class SourceSplitPhase : std::uint8_t {
  kUnsplitCompatibility = 0,
  kStrangPreHalf = 1,
  kStrangPostHalf = 2,
};

// Explicit integration-boundary classes distinguish local active-bin work from
// globally coherent restart/output points and legal PM-refresh surfaces.
enum class StepBoundaryKind : std::uint8_t {
  kLocalActiveBinStep = 0,
  kGlobalSynchronizationPoint = 1,
  kPmRefreshPoint = 2,
  kSnapshotPoint = 3,
  kCheckpointPoint = 4,
};

[[nodiscard]] std::string_view stepBoundaryKindName(StepBoundaryKind kind);
[[nodiscard]] bool isRestartSafeBoundary(StepBoundaryKind kind) noexcept;
[[nodiscard]] bool isOutputSafeBoundary(StepBoundaryKind kind) noexcept;

struct StepBoundaryState {
  StepBoundaryKind kind = StepBoundaryKind::kGlobalSynchronizationPoint;
  bool restart_safe = true;
  bool output_safe = true;
  bool pm_refresh_allowed = true;
  bool local_substep = false;
  // Synchronized coarse-block boundary strictly inside a multirate PM cycle
  // (TreePmPmReusePolicy::kMultirateOuterKick). Positions and rungs are
  // synchronized but momenta carry an open outer PM half-kick, so the
  // boundary is neither output- nor restart-safe; pending snapshot and
  // checkpoint requests are held until the cycle closes (never dropped).
  bool multirate_pm_cycle_interior = false;
};

// Outer PM level of the multirate hierarchical TreePM integrator (2026-10-26;
// stabilized round 3/4). A PM cycle of length T (fixed at its opening by the
// time coordinator) is tiled exactly by m >= 1 coarse blocks of possibly
// different durations Δc_j (re-planned per block from the reference criteria).
// With
//   H = H_kin + H_tree,short + H_PM,long
// the cycle update is the r-RESPA composition
//   K_PM(t0 → t0+T/2) ∘ B_1 ∘ … ∘ B_m ∘ K_PM(t0+T/2 → t0+T)
// where each B_j is the existing tree-only power-of-two block (fine-rung
// tree kicks around all-source drifts) and K_PM(·) is a long-range kick
// with the field evaluated at the cycle endpoint positions. Cosmological
// kicks use the same convention as the block driver: Δu = A ∫dt/a / a_now.
// Disabled plans (enabled = false) reproduce the reference scheme: PM half
// kicks at both endpoints of every coarse block.
struct HierarchicalPmCyclePlan {
  bool enabled = false;
  bool opens_cycle = true;
  bool closes_cycle = true;
  // Opening kick factor ∫_{t0}^{t0+T/2} dt/a / a(t0); used when opens_cycle.
  double opening_pm_kick_factor_code = 0.0;
  // Closing kick lower limit a(t0+T/2) (cosmological runs); the closing
  // integral runs to the actual closing scale factor.
  double cycle_mid_scale_factor = 0.0;
  // Closing kick factor for non-cosmological runs (T/2 in code time).
  double closing_half_interval_time_code = 0.0;
};

// Deterministic restart-boundary verdict. The diagnostic string is part of the
// runtime contract for failed restart requests and must include enough phase,
// scheduler, and PM-cadence context to explain why a checkpoint would not be
// reloadable as persistent truth.
struct RestartBoundaryDecision {
  bool restart_safe = false;
  StepBoundaryKind current_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint;
  StepBoundaryKind last_completed_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint;
  bool inside_kdk_step = false;
  bool last_completed_restart_safe = false;
  bool local_substep_active = false;
  bool pm_refresh_legal = true;
  bool pm_refresh_commit_pending = false;
  std::uint64_t step_index = 0;
  std::optional<std::uint64_t> scheduler_tick;
  std::string diagnostic;
};

enum class PmSyncStage : std::uint8_t {
  kNone = 0,
  kInitialLongRangeBootstrap = 1,
  kScheduledLongRangeRefresh = 2,
};

// Integrator-issued TreePM refresh directive.  Gravity callbacks may own solver
// buffers, but they must consume this placement contract instead of inventing
// legal refresh surfaces internally.  Local-bin force-refresh surfaces require a
// coherent predicted source view because inactive particles are not physically
// drifted in the persistent state during the local substep.
struct PmRefreshDirective {
  enum class Reason : std::uint8_t {
    kNone = 0,
    kInitialForceBootstrap = 1,
    kScheduledForceRefreshStage = 2,
    kSourceMutationForceRefresh = 3,
  };
  bool force_refresh_surface = false;
  bool cadence_opportunity_allowed = false;
  bool initial_cache_bootstrap_allowed = false;
  bool requires_predicted_inactive_sources = false;
  // Rank-local evidence that this stage wants a cadence event.  In distributed
  // execution the workflow globalizes this request before materializing the
  // authoritative event; it is not itself committed cadence truth.
  bool sync_event_requested = false;
  bool has_sync_event = false;
  bool refresh_long_range_field = false;
  bool solver_executed = false;
  PmSyncStage sync_stage = PmSyncStage::kNone;
  std::uint64_t gravity_kick_opportunity = 0;
  std::uint64_t field_version = 0;
  std::uint64_t last_refresh_opportunity = 0;
  std::uint64_t field_built_step_index = 0;
  double field_built_scale_factor = 1.0;
  double force_evaluation_scale_factor = 1.0;
  Reason reason = Reason::kNone;
};

struct CosmologicalStepFactors {
  bool cosmological = false;
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  double dt_time_code = 0.0;
  double time_si_per_code = 1.0;
  double dt_time_si = 0.0;
  double scale_factor_begin = 1.0;
  double scale_factor_midpoint = 1.0;
  double scale_factor_end = 1.0;
  double redshift_begin = 0.0;
  double redshift_end = 0.0;
  double hubble_begin_code = 0.0;
  double hubble_end_code = 0.0;
  double drift_factor_code = 0.0;
  // Particle velocity lanes named velocity_*_peculiar store physical peculiar
  // velocities u = a dx/dt.  The KDK kick sub-operators therefore apply the
  // exact homogeneous Hubble drag u <- (a_begin/a_end) u and a distinct
  // force kick for comoving acceleration kernels, not the drift integral.
  double first_kick_factor_code = 0.0;
  double second_kick_factor_code = 0.0;
  double first_hubble_drag_factor = 1.0;
  double second_hubble_drag_factor = 1.0;
  double hubble_drag_factor = 1.0;
};

class CosmologicalTimeline {
 public:
  explicit CosmologicalTimeline(const LambdaCdmBackground* background = nullptr, double time_si_per_code = 1.0);

  [[nodiscard]] CosmologicalStepFactors prepareStep(
      double current_time_code,
      double current_scale_factor,
      double dt_time_code) const;

  void commitStep(IntegratorState& integrator_state, const CosmologicalStepFactors& step) const;

 private:
  const LambdaCdmBackground* m_background = nullptr;
  double m_time_si_per_code = 1.0;
};

// Hierarchical stepping metadata kept outside particle arrays for auditable ownership.
struct TimeBinContext {
  bool hierarchical_enabled = false;
  std::uint8_t active_bin = 0;
  std::uint8_t max_bin = 0;
};

struct PmSyncEvent {
  std::uint64_t gravity_kick_opportunity = 0;
  bool refresh_long_range_field = false;
  std::uint64_t field_version = 0;
  std::uint64_t last_refresh_opportunity = 0;
  std::uint64_t field_built_step_index = 0;
  double field_built_scale_factor = 1.0;
};

// Restartable PM-cadence authority. The TreePM solver may own field buffers, but
// this state owns the legality of long-range refresh opportunities across restarts.
struct PmSynchronizationPersistentState {
  std::uint64_t cadence_steps = 1;
  std::uint64_t gravity_kick_opportunity = 0;
  std::uint64_t last_refresh_opportunity = 0;
  std::uint64_t field_version = 0;
  std::uint64_t last_refresh_step_index = 0;
  double last_refresh_scale_factor = 1.0;
  bool refresh_commit_pending = false;
  std::uint64_t pending_refresh_opportunity = 0;
  std::uint64_t pending_refresh_field_version = 0;
};

// Scheduler-owned PM cadence state. Gravity callbacks execute the solver, but cadence
// legality and refresh boundaries are represented here to avoid free-floating PM truth.
class PmSynchronizationState {
 public:
  void reset(std::uint64_t cadence_steps = 1);
  // Preview is side-effect free so distributed runtime code can globalize
  // rank-local evidence before any cadence epoch becomes authoritative.
  [[nodiscard]] PmSyncEvent previewKickOpportunity(
      std::uint64_t step_index,
      double scale_factor,
      bool has_long_range_field) const;
  // Commit the accepted opportunity only after the corresponding solver
  // callback completed successfully. Refresh events remain pending until
  // commitRefresh() publishes the new PM field version.
  void commitKickOpportunity(const PmSyncEvent& event);
  [[nodiscard]] PmSyncEvent registerKickOpportunity(
      std::uint64_t step_index,
      double scale_factor,
      bool has_long_range_field);
  void commitRefresh(const PmSyncEvent& event);

  [[nodiscard]] std::uint64_t cadenceSteps() const noexcept { return m_cadence_steps; }
  [[nodiscard]] std::uint64_t gravityKickOpportunity() const noexcept { return m_gravity_kick_opportunity; }
  [[nodiscard]] std::uint64_t fieldVersion() const noexcept { return m_field_version; }
  [[nodiscard]] std::uint64_t lastRefreshOpportunity() const noexcept { return m_last_refresh_opportunity; }
  [[nodiscard]] std::uint64_t lastRefreshStepIndex() const noexcept { return m_last_refresh_step_index; }
  [[nodiscard]] double lastRefreshScaleFactor() const noexcept { return m_last_refresh_scale_factor; }
  [[nodiscard]] bool refreshCommitPending() const noexcept { return m_refresh_commit_pending; }

  [[nodiscard]] PmSynchronizationPersistentState exportPersistentState() const;
  void importPersistentState(const PmSynchronizationPersistentState& persistent_state);

 private:
  std::uint64_t m_cadence_steps = 1;
  std::uint64_t m_gravity_kick_opportunity = 0;
  std::uint64_t m_last_refresh_opportunity = 0;
  std::uint64_t m_field_version = 0;
  std::uint64_t m_last_refresh_step_index = 0;
  double m_last_refresh_scale_factor = 1.0;
  bool m_refresh_commit_pending = false;
  std::uint64_t m_pending_refresh_opportunity = 0;
  std::uint64_t m_pending_refresh_field_version = 0;
};


// Persistent integrator state tracked by the orchestrator.
struct IntegratorState {
  double current_time_code = 0.0;
  double current_scale_factor = 1.0;
  double current_redshift = 0.0;
  double current_hubble_rate_code = 0.0;
  double time_si_per_code = 1.0;
  double dt_time_code = 0.0;
  double last_drift_factor_code = 0.0;
  double last_first_kick_factor_code = 0.0;
  double last_second_kick_factor_code = 0.0;
  double last_first_hubble_drag_factor = 1.0;
  double last_second_hubble_drag_factor = 1.0;
  StepBoundaryKind current_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint;
  StepBoundaryKind last_completed_boundary_kind = StepBoundaryKind::kCheckpointPoint;
  bool inside_kdk_step = false;
  bool last_completed_restart_safe = true;
  std::uint64_t step_index = 0;
  TimeStepScheme scheme = TimeStepScheme::kKickDriftKick;
  TimeBinContext time_bins;
  PmSynchronizationState pm_sync_state;
  bool pm_refresh_enabled = false;
  bool pm_long_range_field_valid = false;
  // Authoritative SimulationState gravity generation represented by the
  // committed long-range PM field. Zero means unknown/unbound after restart.
  std::uint64_t pm_source_generation = 0;
};

// Explicit compact active-set descriptor with optional subset spans.
struct ActiveSetDescriptor {
  std::span<const std::uint32_t> particle_indices;
  std::span<const std::uint32_t> cell_indices;
  bool particles_are_subset = false;
  bool cells_are_subset = false;
  // Runtime-truth metadata: active descriptors are derived products of the
  // scheduler/state generation that created them, never independent active-set
  // authority. Stale descriptors must fail before solver kernels consume them.
  bool particles_from_scheduler = false;
  bool cells_from_scheduler = false;
  bool has_generation_metadata = false;
  std::uint64_t source_particle_index_generation = 0;
  std::uint64_t source_cell_index_generation = 0;
  std::uint64_t source_scheduler_tick = 0;
  // Optional distributed consensus supplied by the workflow owner.  Local
  // subset sizes alone cannot decide whether a rank is at the same global
  // synchronization boundary as its peers (including zero-work ranks).
  bool has_global_synchronization_metadata = false;
  bool globally_complete_active_set = false;

  [[nodiscard]] bool hasParticleSubset(std::size_t total_particle_count) const noexcept;
  [[nodiscard]] bool hasCellSubset(std::size_t total_cell_count) const noexcept;
};

// Per-stage execution context passed to callbacks without mutating interface contracts.
struct StepContext {
  SimulationState& state;
  IntegratorState& integrator_state;
  ActiveSetDescriptor active_set;
  GravityParticleKernelView active_gravity_particles;
  bool has_active_gravity_particles = false;
  TransientStepWorkspace* workspace = nullptr;
  const LambdaCdmBackground* cosmology_background = nullptr;
  const ModePolicy* mode_policy = nullptr;
  ProfilerSession* profiler_session = nullptr;
  // Non-owning scheduler seams used by source modules that append new particles at
  // the legal source mutation boundary. Unit tests may leave these null.
  HierarchicalTimeBinScheduler* particle_scheduler = nullptr;
  HierarchicalTimeBinScheduler* gas_cell_scheduler = nullptr;
  // Per-step append-only sink for authoritative numeric IDs created by legal
  // source mutations. The workflow ownership ledger consumes these IDs before
  // any subsequent decomposition or migration validation.
  std::vector<std::uint64_t>* newly_created_particle_ids = nullptr;
  CosmologicalStepFactors timeline_step;
  StepBoundaryState boundary;
  PmRefreshDirective pm_refresh_directive;
  IntegrationStage stage = IntegrationStage::kGravityKickPre;
  SourceSplitPhase source_split_phase = SourceSplitPhase::kUnsplitCompatibility;
  // Transient integrator-issued block KDK directive. Rung/bin authority stays
  // in HierarchicalTimeBinScheduler; these factors describe one kick surface.
  // Forces use scale-free A, velocities store p/a at the common position epoch.
  struct HierarchicalKdkDirective {
    bool enabled = false;
    bool include_long_range_force = false;
    bool synchronization_end = false;
    bool force_only_synchronization = false;
    std::uint64_t coarse_source_generation = 0U;
    std::uint64_t evaluation_tick = 0U;
    std::array<double, 13U> tree_kick_factor_code{};
    double pm_kick_factor_code = 0.0;
  } hierarchical_kdk;
};

inline constexpr std::size_t k_hierarchical_timeline_capacity = 8193U;

// Materialize a previously requested PM cadence event without committing it.
// The caller supplies the authoritative long-range-validity verdict: serial
// dispatch uses the local value, while distributed workflow dispatch supplies
// a rank-global verdict.
void materializePmRefreshDirective(
    StepContext& context,
    bool has_long_range_field);

// Dependency-safe dispatch seam used by workflow-level typed resource views.
// Core retains KDK/timeline invariants while the workflow owns task selection;
// no workflow type appears in this interface.
using StageDispatchFunction =
    std::function<void(StepContext& context, bool require_output_safe_boundary)>;

// Stage-bound handler interface implemented by gravity, hydro, source, analysis, and output modules.
// Ownership stays with the caller; StepOrchestrator stores non-owning pointers in
// deterministic registration order.  Each handler declares the exact stage set it
// may receive, and production dispatch only iterates the current stage bucket.
class IntegrationCallback {
 public:
  virtual ~IntegrationCallback() = default;

  [[nodiscard]] virtual std::string_view callbackName() const = 0;
  [[nodiscard]] virtual std::span<const IntegrationStage> integrationStages() const = 0;
  [[nodiscard]] virtual std::span<const StageContract> stageContracts() const = 0;
  virtual void onStage(StepContext& context) = 0;
};

// Stage scheduler isolates ordering from solver implementation details.
class StageScheduler {
 public:
  [[nodiscard]] std::vector<IntegrationStage> schedule(
      const IntegratorState& integrator_state,
      const ActiveSetDescriptor& active_set) const;

  [[nodiscard]] static std::span<const IntegrationStage> kickDriftKickOrder();
};

[[nodiscard]] bool isCanonicalIntegrationStageOrder(std::span<const IntegrationStage> ordered_stages);

// Single authoritative step orchestrator for current baseline stepping.
class StepOrchestrator {
 public:
  explicit StepOrchestrator(StageScheduler scheduler = {});
  // Collisionless synchronized-block KDK. All positions share each fine
  // epoch; target kicks/forces use scheduler subsets. Bins remain fixed within
  // a coarse block, and PM kicks occur only at its two endpoint surfaces.
  // Distributed callers supply their existing collective failure authority
  // for local preparation between stage dispatches; empty is serial-only.
  void executeHierarchicalBlockWithDispatcher(
      SimulationState& state, IntegratorState& integrator_state,
      HierarchicalTimeBinScheduler& particle_scheduler,
      HierarchicalTimeBinScheduler& gas_cell_scheduler,
      const StageDispatchFunction& dispatcher,
      const LambdaCdmBackground* cosmology_background,
      TransientStepWorkspace& workspace, const ModePolicy* mode_policy,
      ProfilerSession* profiler_session = nullptr,
      const std::function<void(std::exception_ptr, std::string_view)>& preparation_failure_sync = {},
      const HierarchicalPmCyclePlan& pm_cycle = {}) const;

  void registerCallback(IntegrationCallback& callback);
  [[nodiscard]] std::size_t callbackCount() const noexcept;
  [[nodiscard]] std::span<IntegrationCallback* const> handlersFor(IntegrationStage stage) const noexcept;
  [[nodiscard]] std::span<const StageContract> contractsFor(IntegrationStage stage) const noexcept;
  [[nodiscard]] std::optional<StageContract> contractForHandlerStage(
      const IntegrationCallback& callback,
      IntegrationStage stage) const noexcept;

  void executeOutputBoundary(
      SimulationState& state,
      IntegratorState& integrator_state,
      ProfilerSession* profiler_session = nullptr,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

  void executeOutputBoundaryWithDispatcher(
      SimulationState& state,
      IntegratorState& integrator_state,
      const StageDispatchFunction& dispatcher,
      ProfilerSession* profiler_session = nullptr,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

  void executeSingleStep(
      SimulationState& state,
      IntegratorState& integrator_state,
      ActiveSetDescriptor active_set,
      const LambdaCdmBackground* cosmology_background,
      TransientStepWorkspace* workspace = nullptr,
      const ModePolicy* mode_policy = nullptr,
      ProfilerSession* profiler_session = nullptr,
      std::optional<std::uint64_t> expected_scheduler_tick = std::nullopt,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

  void executeSingleStepWithDispatcher(
      SimulationState& state,
      IntegratorState& integrator_state,
      ActiveSetDescriptor active_set,
      const StageDispatchFunction& dispatcher,
      const LambdaCdmBackground* cosmology_background,
      TransientStepWorkspace* workspace = nullptr,
      const ModePolicy* mode_policy = nullptr,
      ProfilerSession* profiler_session = nullptr,
      std::optional<std::uint64_t> expected_scheduler_tick = std::nullopt,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

  void executeSchedulerSubstep(
      SimulationState& state,
      IntegratorState& integrator_state,
      const HierarchicalTimeBinScheduler& scheduler,
      std::span<const std::uint32_t> active_particle_indices,
      std::span<const std::uint32_t> active_cell_indices,
      const LambdaCdmBackground* cosmology_background,
      TransientStepWorkspace* workspace = nullptr,
      const ModePolicy* mode_policy = nullptr,
      ProfilerSession* profiler_session = nullptr,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

 private:
  StageScheduler m_scheduler;
  void dispatchStageHandlers(StepContext& context, bool require_output_safe_boundary) const;

  std::array<std::vector<IntegrationCallback*>, integrationStageCount()> m_handlers_by_stage;
  std::array<std::vector<StageContract>, integrationStageCount()> m_contracts_by_stage;
  std::size_t m_callback_count = 0;
};

// da/dt = a H(a) for standard FLRW backgrounds.
[[nodiscard]] double computeScaleFactorRate(const LambdaCdmBackground& background, double scale_factor);

// Production scale-factor evolution: invert the FLRW time integral dt = integral da/(a H(a)).
[[nodiscard]] double advanceScaleFactorByCosmicTime(
    const LambdaCdmBackground& background,
    double scale_factor,
    double dt_time_si,
    std::uint32_t midpoint_samples = 64);

// dt estimate for an intended delta-a increment around the current scale factor.
[[nodiscard]] double estimateDeltaTimeFromScaleFactorStep(
    const LambdaCdmBackground& background,
    double scale_factor,
    double delta_scale_factor);

// Drift prefactor integral: integral_{a0}^{a1} da / (a^2 H(a)).
[[nodiscard]] double computeComovingDriftFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples = 16);

// Kick prefactor for comoving acceleration terms proportional to 1/a.
[[nodiscard]] double computeComovingKickFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples = 16);

// Hubble drag factor for dv/dt = -H(a) v over [a0, a1].
[[nodiscard]] double computeHubbleDragFactor(double scale_factor_begin, double scale_factor_end);

[[nodiscard]] StepBoundaryState classifyStepBoundary(
    const SimulationState& state,
    const ActiveSetDescriptor& active_set,
    bool scheduler_owned_substep,
    StepBoundaryKind requested_kind = StepBoundaryKind::kGlobalSynchronizationPoint);

[[nodiscard]] RestartBoundaryDecision evaluateRestartBoundary(
    const IntegratorState& integrator_state,
    std::optional<std::uint64_t> scheduler_tick = std::nullopt);
[[nodiscard]] bool canWriteRestart(
    const IntegratorState& integrator_state,
    std::optional<std::uint64_t> scheduler_tick = std::nullopt);

void assertCanWriteSnapshotAtBoundary(const IntegratorState& integrator_state);
void assertCanWriteCheckpointAtBoundary(
    const IntegratorState& integrator_state,
    std::optional<std::uint64_t> scheduler_tick = std::nullopt);

}  // namespace cosmosim::core

#include "cosmosim/core/time_scheduler.hpp"

