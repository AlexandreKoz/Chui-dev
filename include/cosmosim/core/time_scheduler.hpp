#pragma once

#include "cosmosim/core/time_integration.hpp"

namespace cosmosim::core {

// Hot metadata sidecar for element-local time-bin ownership and scheduling state.
struct TimeBinHotMetadata {
  std::vector<std::uint8_t> bin_index;
  std::vector<std::uint64_t> next_activation_tick;
  std::vector<std::uint8_t> active_flag;
  std::vector<std::uint8_t> pending_bin_index;

  [[nodiscard]] std::size_t size() const noexcept { return bin_index.size(); }
};

// Diagnostics emitted by hierarchical scheduler for auditing and pathological collapse detection.
struct TimeBinDiagnostics {
  std::vector<std::uint32_t> occupancy_by_bin;
  std::vector<std::uint32_t> active_count_by_bin;
  std::uint32_t active_elements = 0;
  std::uint32_t promoted_elements = 0;
  std::uint32_t demoted_elements = 0;
  std::uint32_t clipped_to_min_dt = 0;
  std::uint32_t clipped_to_max_dt = 0;
  std::uint32_t illegal_transition_attempts = 0;
  std::uint32_t deferred_coarsening_events = 0;
  std::uint32_t collapse_candidates = 0;
  double active_fraction = 0.0;
  std::uint8_t most_active_bin = 0;
};

struct TimeBinMappingResult {
  std::uint8_t bin_index = 0;
  bool clipped_to_min = false;
  bool clipped_to_max = false;
};

// Auditable labels attached to solver timestep criteria before scheduler reconciliation.
enum class TimeStepCandidateSource : std::uint8_t {
  kHydroCfl = 0,
  kGravityAcceleration = 1,
  kCosmologyExpansion = 2,
  kSourceTerm = 3,
  kUserClamp = 4,
};

[[nodiscard]] std::string_view timeStepCandidateSourceName(TimeStepCandidateSource source);

struct TimeStepCandidateSubmission {
  std::uint32_t element_index = 0;
  double dt_time_code = 0.0;
  TimeStepCandidateSource source = TimeStepCandidateSource::kUserClamp;
  std::string label;
};

struct TimeStepReconciliationResult {
  std::uint32_t submitted_candidates = 0;
  std::uint32_t elements_with_candidates = 0;
  std::uint32_t committed_transition_requests = 0;
  std::uint32_t clipped_to_min_dt = 0;
  std::uint32_t clipped_to_max_dt = 0;
  std::array<std::uint32_t, 5> limiting_candidates_by_source{};
  TimeStepCandidateSource dominant_limiting_source = TimeStepCandidateSource::kUserClamp;
};

// Persisted scheduler state required for exact restart continuation.
struct TimeBinPersistentState {
  std::uint64_t current_tick = 0;
  std::uint8_t max_bin = 0;
  std::vector<std::uint8_t> bin_index;
  std::vector<std::uint64_t> next_activation_tick;
  std::vector<std::uint8_t> active_flag;
  std::vector<std::uint8_t> pending_bin_index;
};

// Typed limits that normalize physical timestep proposals into the discrete bin hierarchy.
struct TimeStepLimits {
  double min_dt_time_code = 0.0;
  double max_dt_time_code = 0.0;
  std::uint8_t max_bin = 0;
};

// Compact inputs for conservative CFL and gravity criteria hooks.
struct CflTimeStepInput {
  double cell_width_code = 0.0;
  double flow_speed_code = 0.0;
  double sound_speed_code = 0.0;
};

// Discrete stability rule of the hydro time-integration family.
//  * kPerAxisMaximum: dt = C min_d L dx_d / (|v_d| + c). Appropriate for
//    corner-transport-upwind (CTU) integration, whose transverse corrections
//    make the stability region per-axis (simplified 3D variant: nu_d <= 1/2,
//    so C = 0.4 is admissible). Historical default rule.
//  * kDimensionSummed: dt = C / sum_{active d} (|v_d| + c) / (L dx_d), i.e.
//    sum_d nu_d <= C. Required by the directionally UNSPLIT MUSCL-Hancock
//    update U^{n+1} = U^n - sum_d dt/(a dx_d) (F_{d,+} - F_{d,-}): for linear
//    advection with first-order states its amplification factor is that of
//    unsplit donor cell, stable iff sum_d nu_d <= 1 (the Hancock half step
//    does not enlarge the region). With C = 0.4 the 1D step is unchanged and
//    3D diagonal flow is kept at sum nu <= 0.4 instead of up to 1.2.
// L = a in comoving coordinates (peculiar speeds, comoving widths), 1 otherwise.
enum class HydroCflStabilityRule : std::uint8_t {
  kPerAxisMaximum = 0,
  kDimensionSummed = 1,
};

struct DirectionalCflTimeStepInput {
  // Ordered x, y, z stored-coordinate widths. velocity_axis_code and
  // sound_speed_code are physical peculiar/signal speeds. For a comoving
  // coordinate frame the crossing length is a * dx_comoving, so the CFL
  // time scales linearly with scale_factor. Proper-coordinate runs use the
  // stored width directly.
  std::array<double, 3> cell_width_axis_code{0.0, 0.0, 0.0};
  std::array<double, 3> velocity_axis_code{0.0, 0.0, 0.0};
  double sound_speed_code = 0.0;
  CoordinateFrame coordinate_frame = CoordinateFrame::kPhysical;
  double scale_factor = 1.0;
  HydroCflStabilityRule stability_rule = HydroCflStabilityRule::kPerAxisMaximum;
  // Axes along which the hydro solver actually transports (an axis whose
  // single cell spans the whole domain carries no flux differences). An
  // inactive axis is excluded from both rules.
  std::array<bool, 3> active_axis{true, true, true};
};

// Active hydro axes for a cell: an axis is inactive only when one cell spans
// the whole domain extent along it (degenerate 1D/2D setups). This criterion
// is per cell and rank-independent, so all ranks agree.
[[nodiscard]] std::array<bool, 3> activeHydroAxes(
    const std::array<double, 3>& cell_width_axis_code,
    const std::array<double, 3>& domain_extent_axis_code);

// numerics.hydro_integration_scheme -> stability rule (ctu: per-axis,
// muscl_hancock: dimension-summed). Single mapping used by TimeCoordinator
// and the hydro accepted-step guard.
[[nodiscard]] HydroCflStabilityRule hydroCflStabilityRuleFor(const SimulationConfig& config) noexcept;
// Configured domain extent in the stored coordinate frame
// (cosmology.box_size_{x,y,z}_mpc_comoving, the convention of the gravity
// domain and the fixed-grid layout).
[[nodiscard]] std::array<double, 3> hydroDomainExtentCode(const SimulationConfig& config) noexcept;

struct HydroCflDiagnostics {
  std::uint32_t local_row = 0;
  std::uint64_t gas_cell_id = 0;
  bool has_gas_cell_id = false;
  std::uint64_t patch_id = 0;
  bool has_patch_id = false;
  std::uint32_t patch_row = 0;
  bool has_patch_row = false;
  double proposed_dt_time_code = std::numeric_limits<double>::infinity();
  double accepted_dt_time_code = 0.0;
  double cfl_number = 0.0;
  double safety_factor = std::numeric_limits<double>::infinity();
  std::array<double, 3> cell_width_axis_code{0.0, 0.0, 0.0};
  std::array<double, 3> velocity_axis_code{0.0, 0.0, 0.0};
  double sound_speed_code = 0.0;
  CoordinateFrame coordinate_frame = CoordinateFrame::kPhysical;
  double scale_factor = 1.0;
  std::uint8_t limiting_axis = 0;
  HydroCflStabilityRule stability_rule = HydroCflStabilityRule::kPerAxisMaximum;
};

struct GravityTimeStepInput {
  // Length and acceleration must use the same coordinate system. For
  // cosmological comoving softening, convert TreePM's scale-free A to the
  // coordinate acceleration A/a^3 before calling computeGravityTimeStep.
  double softening_length_code = 0.0;
  double acceleration_magnitude_code = 0.0;
};

struct ComovingGravityTimeStepInput {
  double softening_length_comoving_code = 0.0;
  // TreePM's scale-free kernel A = G sum_j m_j (x_j-x_i)/|x_j-x_i|^3.
  double scale_free_acceleration_magnitude_code = 0.0;
  double scale_factor = 1.0;
};

// Independent coarse PM displacement criterion. The force/softening criterion
// still selects particle rungs; this limits the interval between PM kicks.
struct ComovingDisplacementTimeStepInput {
  double mesh_or_split_length_comoving_code = 0.0;
  double velocity_magnitude_peculiar_code = 0.0;
  double scale_free_acceleration_magnitude_code = 0.0;
  double scale_factor = 1.0;
};

[[nodiscard]] double computeComovingDisplacementTimeStep(
    const ComovingDisplacementTimeStepInput& input);

// Narrow runtime view consumed by adaptive timestep criteria. It is deliberately
// smaller than SimulationState: criteria loops read only velocities, sound speed,
// softening, source-module rates, and force-output spans needed to propose bins.
struct TimeStepParticleCriteriaView {
  std::span<const double> velocity_x_peculiar;
  std::span<const double> velocity_y_peculiar;
  std::span<const double> velocity_z_peculiar;
  std::span<const std::uint32_t> species_tag;
  bool homogeneous_dmo_species = false;
  std::span<const double> gravity_softening_comoving;
  std::span<const double> accel_x_comoving;
  std::span<const double> accel_y_comoving;
  std::span<const double> accel_z_comoving;
  std::span<const std::uint32_t> black_hole_particle_index;
  std::span<const double> black_hole_subgrid_mass_code;
  std::span<const double> black_hole_accretion_rate_code;
};

struct TimeStepGasCellCriteriaView {
  std::span<const std::uint32_t> gas_particle_index_by_cell;
  std::span<const std::uint64_t> gas_cell_id_by_cell;
  std::span<const std::uint64_t> patch_id_by_cell;
  std::span<const std::uint32_t> patch_row_by_cell;
  std::span<const double> cell_width_x_code;
  std::span<const double> cell_width_y_code;
  std::span<const double> cell_width_z_code;
  std::span<const double> cell_mass_code;
  std::span<const double> velocity_x_peculiar;
  std::span<const double> velocity_y_peculiar;
  std::span<const double> velocity_z_peculiar;
  std::span<const double> density_code;
  std::span<const double> temperature_code;
  std::span<const double> sound_speed_code;
  // Peculiar-velocity divergence evaluated with physical cell spacing, in
  // inverse code-time units. Invalid/non-patch cells carry NaN and are ignored.
  std::span<const double> velocity_divergence_code;
  // Local conservative parabolic limit for the configured metal-diffusion
  // operator. Infinity means diffusion imposes no limit for this cell.
  std::span<const double> metal_diffusion_dt_code;
  std::span<const double> accel_x_comoving;
  std::span<const double> accel_y_comoving;
  std::span<const double> accel_z_comoving;
  // G1.2 (recursive AMR subcycling): per-cell factor r^(L - L_min,global)
  // by which the COARSE interval may exceed the cell's own hydro CFL bound,
  // because a level-L cell advances with dt / r^(L - L_min). Applies to the
  // hydro CFL candidate ONLY; gravity, expansion, source, diffusion and
  // chemistry candidates are never relaxed (those operators run at the coarse
  // cadence). Empty = 1 (synchronized AMR, fixed grid).
  std::span<const double> hydro_cfl_level_relaxation;
};

struct AdaptiveTimeStepCriteriaView {
  TimeStepParticleCriteriaView particles;
  TimeStepGasCellCriteriaView gas_cells;
};

// Criterion hooks return a finite positive timestep in code units. Positive
// infinity is the sole explicit sentinel for "this criterion imposes no limit";
// NaN, -infinity, zero, and negative values are contract violations and fail closed.
using CriteriaHook = std::function<double(std::uint32_t)>;

struct TimeStepCriteriaHooks {
  CriteriaHook cfl_hook;
  CriteriaHook gravity_hook;
  CriteriaHook source_hook;
  CriteriaHook user_clamp_hook;
};

class TimeStepCriteriaRegistry {
 public:
  void registerCflHook(CriteriaHook hook);
  void registerGravityHook(CriteriaHook hook);
  void registerSourceHook(CriteriaHook hook);
  void registerUserClampHook(CriteriaHook hook);

  [[nodiscard]] const TimeStepCriteriaHooks& hooks() const noexcept;

 private:
  TimeStepCriteriaHooks m_hooks;
};

// Physical representation policy. The logical scheduler authority is unchanged:
// max_bin==0 may store the homogeneous rung-zero state implicitly while the
// generic representation retains the full hierarchical machinery.
enum class SchedulerRepresentation : std::uint8_t {
  kGenericHierarchical = 0,
  kUniformRungZero = 1,
};

// Integer timeline scheduler with power-of-two bins and compact active set extraction.
class RetainedCapacityTransaction;

class HierarchicalTimeBinScheduler {
 public:
  static constexpr std::uint8_t k_unset_pending_bin = 0xFF;
  static constexpr std::uint8_t k_max_representable_bin = 62;

  explicit HierarchicalTimeBinScheduler(std::uint8_t max_bin = 0);

  void reset(std::uint32_t element_count, std::uint8_t initial_bin, std::uint64_t start_tick = 0);
  void appendElements(std::uint32_t new_element_count, std::uint8_t initial_bin, std::uint64_t first_activation_tick);
  void planAppendCapacity(RetainedCapacityTransaction& plan,
                          std::uint32_t new_element_count, std::uint8_t initial_bin);
  void setElementBin(std::uint32_t element_index, std::uint8_t bin_index, std::uint64_t current_tick);
  void submitCandidateTimeStep(
      std::uint32_t element_index,
      double dt_time_code,
      const TimeStepLimits& limits,
      TimeStepCandidateSource source,
      std::string_view label = {});
  void submitCandidateBin(
      std::uint32_t element_index,
      std::uint8_t target_bin,
      TimeStepCandidateSource source,
      std::string_view label = {});
  [[nodiscard]] TimeStepReconciliationResult reconcileCandidateTransitions();

  [[nodiscard]] std::span<const std::uint32_t> activeElements() const noexcept;

  [[nodiscard]] std::uint64_t currentTick() const noexcept;
  [[nodiscard]] std::uint8_t maxBin() const noexcept;
  [[nodiscard]] std::uint32_t elementCount() const noexcept;
  [[nodiscard]] SchedulerRepresentation representation() const noexcept { return m_representation; }
  [[nodiscard]] std::uint8_t binIndex(std::uint32_t element_index) const;
  [[nodiscard]] std::uint64_t nextActivationTick(std::uint32_t element_index) const;
  [[nodiscard]] bool isElementActive(std::uint32_t element_index) const;
  [[nodiscard]] std::uint8_t pendingBinIndex(std::uint32_t element_index) const;

  [[nodiscard]] bool isBinActiveAtTick(std::uint8_t bin_index, std::uint64_t tick) const;
  [[nodiscard]] std::uint64_t binPeriodTicks(std::uint8_t bin_index) const;

  std::span<const std::uint32_t> beginSubstep();
  void endSubstep();
  // Close an all-active aligned block endpoint without advancing its tick.
  // Every interval is closed; next block may change bins and reopen here.
  void closeSynchronizationPoint();
  // Reconcile an entire closed aligned block using already-submitted physical
  // candidates. Reserve exact bin occupancies before changing membership;
  // reuse the existing candidate lane, without another N-element bin array.
  void planSynchronizedCandidateCapacity(RetainedCapacityTransaction& plan);
  void commitSynchronizedCandidates();

  [[nodiscard]] const TimeBinHotMetadata& hotMetadata() const noexcept;
  [[nodiscard]] const TimeBinDiagnostics& diagnostics() const noexcept;
  // Population-scale scheduler accounting distinguishes live logical elements,
  // retained vector capacity, and the historical retained-capacity high-water.
  [[nodiscard]] std::uint64_t logicalSizeBytes() const;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const;
  [[nodiscard]] std::uint64_t ownedCapacityHighWaterBytes() const noexcept;

  [[nodiscard]] TimeBinPersistentState exportPersistentState() const;
  // Same content as exportPersistentState, written into caller-owned storage
  // whose lane capacity is reused (no allocation when capacity suffices).
  void exportPersistentStateInto(TimeBinPersistentState& persistent_state) const;
  void importPersistentState(const TimeBinPersistentState& persistent_state);

 private:
  std::uint8_t clampBin(std::uint8_t requested) const noexcept;
  void eraseFromBin(std::uint32_t element_index, std::uint8_t bin_index);
  void insertIntoBin(std::uint32_t element_index, std::uint8_t bin_index);
  void rebuildActiveSet();
  void requestBinTransition(std::uint32_t element_index, std::uint8_t target_bin);
  void applyPendingTransitions();
  void validateInternalState(std::string_view source_label) const;
  void validateTransitionRequest(
      std::uint32_t element_index,
      std::uint8_t target_bin,
      std::string_view source_label) const;
  void refreshOwnedCapacityHighWater();
  void materializeGenericRepresentationForAppend();

  std::uint64_t m_current_tick = 0;
  std::uint32_t m_element_count = 0;
  std::uint8_t m_max_bin = 0;
  SchedulerRepresentation m_representation = SchedulerRepresentation::kGenericHierarchical;
  std::uint8_t m_uniform_pending_bin_index = k_unset_pending_bin;
  TimeBinHotMetadata m_hot;
  std::vector<std::vector<std::uint32_t>> m_elements_by_bin;
  std::vector<std::size_t> m_position_in_bin;
  std::vector<std::uint32_t> m_active_elements;
  // Reusable radix-sort staging for deterministic active-index order without O(A log A) comparison sorting.
  std::vector<std::uint32_t> m_active_sort_scratch;
  TimeBinDiagnostics m_diagnostics;
  std::vector<std::uint8_t> m_candidate_bin_index;
  std::vector<TimeStepCandidateSource> m_candidate_source;
  TimeStepReconciliationResult m_last_reconciliation;
  std::uint64_t m_owned_capacity_high_water_bytes = 0U;
  bool m_substep_open = false;
};

[[nodiscard]] ActiveSetDescriptor makeSchedulerActiveSetDescriptor(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices = {});

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state);

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    std::uint64_t expected_scheduler_tick);

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    const HierarchicalTimeBinScheduler& scheduler);

// Scheduler-backed reorder by rung. This is the only production-safe time-bin
// ordering path: it consumes scheduler authority directly instead of public
// particle/cell mirrors.
[[nodiscard]] ParticleReorderMap buildParticleReorderMapByScheduler(
    const SimulationState& state,
    const HierarchicalTimeBinScheduler& scheduler);


// Stable scheduler identity namespace. Current production scheduling is particle-row
// backed, with gas cells treated as particle-bound finite-volume carriers. Future AMR
// patches or mesh cells must enter through an explicit kind/stable_id pair rather
// than by overloading local row indices.
enum class ScheduledElementKind : std::uint8_t {
  kParticle = 0,
  kParticleBoundGasCell = 1,
  kAmrPatch = 2,
};

struct ScheduledElementKey {
  ScheduledElementKind kind = ScheduledElementKind::kParticle;
  std::uint64_t stable_id = 0;

  [[nodiscard]] friend bool operator==(const ScheduledElementKey&, const ScheduledElementKey&) = default;
};

struct TimeBinSchedulerIdentityRecord {
  // Full scheduler authority for one physical element identity. Migration packets
  // may carry time_bin as a diagnostic mirror, but exact continuation requires
  // these scheduler-owned lanes matched by stable particle_id/gas_cell_id.
  std::uint64_t element_id = 0;
  std::uint8_t bin_index = 0;
  std::uint64_t next_activation_tick = 0;
  std::uint8_t pending_bin_index = HierarchicalTimeBinScheduler::k_unset_pending_bin;
};

[[nodiscard]] std::vector<TimeBinSchedulerIdentityRecord> exportParticleSchedulerIdentityRecords(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> particle_indices);

void attachSchedulerFieldsToParticleMigrationRecords(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> particle_indices,
    std::span<ParticleMigrationRecord> records);

[[nodiscard]] TimeBinPersistentState rebuildSchedulerPersistentStateFromMigrationRecords(
    std::uint64_t current_tick,
    std::uint8_t max_bin,
    std::span<const ParticleMigrationRecord> records,
    std::span<const std::uint64_t> destination_element_ids);

void rebuildSchedulerFromParticleMigrationRecords(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const ParticleMigrationRecord> records,
    std::span<const std::uint64_t> destination_particle_ids);

[[nodiscard]] TimeBinPersistentState remapSchedulerPersistentStateByParticleId(
    const TimeBinPersistentState& source_state,
    std::span<const std::uint64_t> source_particle_ids,
    std::span<const std::uint64_t> destination_particle_ids);

void remapSchedulerByParticleId(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const std::uint64_t> source_particle_ids,
    std::span<const std::uint64_t> destination_particle_ids);

void remapSchedulerByParticleReorderMap(
    HierarchicalTimeBinScheduler& scheduler,
    const ParticleReorderMap& reorder_map);

[[nodiscard]] TimeBinPersistentState rebuildSchedulerPersistentStateFromIdentityRecords(
    std::uint64_t current_tick,
    std::uint8_t max_bin,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    std::span<const std::uint64_t> destination_element_ids);

void rebuildSchedulerFromParticleIdentityRecords(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    std::span<const std::uint64_t> destination_particle_ids);

// Gas-cell scheduler authority is keyed by stable gas_cell_id, never by a
// parent particle row.  These helpers are the remap/restart seam used by
// cell migration and restart reconstruction.
[[nodiscard]] std::vector<TimeBinSchedulerIdentityRecord> exportGasCellSchedulerIdentityRecords(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> cell_indices);

void rebuildSchedulerFromGasCellIdentityRecords(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    const SimulationState& state);

void syncGasCellTimeBinMirrorsFromGasCellScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state);

// Legacy compatibility helper only.  New workflow paths must use the
// gas-cell scheduler helper above so parentless and split cells remain
// independently schedulable.
void syncGasCellTimeBinMirrorsFromParticleScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state);

[[nodiscard]] TimeBinMappingResult mapDtToTimeBin(double dt_time_code, const TimeStepLimits& limits);
[[nodiscard]] double binIndexToDt(std::uint8_t bin_index, const TimeStepLimits& limits);
enum class TimeBinMirrorDomain : std::uint8_t {
  kParticles = 0,
  kCells = 1,
  kParticlesAndCells = 2,
};

void syncTimeBinMirrorsFromScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state,
    TimeBinMirrorDomain domain = TimeBinMirrorDomain::kParticles);
[[nodiscard]] bool timeBinMirrorsMatchScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain = TimeBinMirrorDomain::kParticles);
void debugAssertTimeBinMirrorAuthorityInvariant(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain = TimeBinMirrorDomain::kParticles);

[[nodiscard]] double computeCflTimeStep(const CflTimeStepInput& input, double c_cfl);
[[nodiscard]] double computeDirectionalCflTimeStep(const DirectionalCflTimeStepInput& input, double c_cfl);
[[nodiscard]] HydroCflDiagnostics makeHydroCflDiagnostics(
    std::uint32_t local_row,
    const DirectionalCflTimeStepInput& input,
    double c_cfl,
    double accepted_dt_time_code,
    std::uint64_t gas_cell_id = 0,
    std::optional<std::uint64_t> patch_id = std::nullopt,
    std::optional<std::uint32_t> patch_row = std::nullopt);
void assertHydroCflStable(
    const HydroCflDiagnostics& diagnostics,
    double relative_tolerance = 1.0e-12);
[[nodiscard]] double computeGravityTimeStep(const GravityTimeStepInput& input, double eta);
[[nodiscard]] double computeComovingGravityTimeStep(
    const ComovingGravityTimeStepInput& input,
    double eta);

[[nodiscard]] double combineTimeStepCriteria(
    std::uint32_t element_index,
    const TimeStepCriteriaHooks& hooks,
    double fallback_dt_time_code);

[[nodiscard]] double computeCosmologyExpansionTimeStep(
    const LambdaCdmBackground& background,
    double scale_factor,
    double max_delta_ln_a,
    double max_hubble_time_fraction,
    double time_si_per_code = 1.0);

}  // namespace cosmosim::core
