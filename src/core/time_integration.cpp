#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"

#include "cosmosim/core/checked_arithmetic.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace cosmosim::core {
namespace {

// KDK remains the particle/gravity authority.  The two source dispatches are
// the deterministic Strang reference S(dt/2) H(dt) S(dt/2).  SourceRuntime
// keeps stochastic/discrete events out of the pre half-step and applies them
// exactly once at their absolute end-of-interval epoch in the post half-step.
constexpr std::array<IntegrationStage, 9> k_kick_drift_kick_order = {
    IntegrationStage::kGravityKickPre,
    IntegrationStage::kDrift,
    IntegrationStage::kForceRefresh,
    IntegrationStage::kSourceTerms,
    IntegrationStage::kHydroUpdate,
    IntegrationStage::kSourceTerms,
    IntegrationStage::kGravityKickPost,
    IntegrationStage::kAnalysisHooks,
    IntegrationStage::kOutputCheck,
};

// Deterministic LSD radix ordering for unique 32-bit local indices.  Four
// contiguous byte passes avoid comparison sorting while preserving the existing
// O(1) swap-with-last bin membership updates.  Scratch capacity is retained
// across substeps so steady-state scheduling does not allocate.
void radixSortLocalIndices(
    std::vector<std::uint32_t>& values,
    std::vector<std::uint32_t>& scratch) {
  if (values.size() < 2U) {
    return;
  }
  scratch.resize(values.size());
  constexpr std::size_t k_radix = 256U;
  constexpr std::uint32_t k_mask = 0xffU;
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    std::array<std::size_t, k_radix> offsets{};
    for (const std::uint32_t value : values) {
      ++offsets[(value >> shift) & k_mask];
    }
    std::size_t prefix = 0U;
    for (std::size_t bucket = 0; bucket < k_radix; ++bucket) {
      const std::size_t count = offsets[bucket];
      offsets[bucket] = prefix;
      prefix += count;
    }
    for (const std::uint32_t value : values) {
      scratch[offsets[(value >> shift) & k_mask]++] = value;
    }
    values.swap(scratch);
  }
}


[[nodiscard]] std::string stageContractErrorPrefix(std::string_view callback_name, IntegrationStage stage) {
  return "IntegrationCallback '" + std::string(callback_name) + "' stage contract for '" +
      std::string(integrationStageName(stage)) + "' ";
}

void validateStageContractShape(std::string_view callback_name, const StageContract& contract) {
  const auto prefix = stageContractErrorPrefix(callback_name, contract.stage);

  if (contract.owner != StageSubsystem::kAnalysis && contract.required_inputs == StageDataDomain::kNone) {
    throw std::invalid_argument(prefix + "must declare required_inputs");
  }
  if (contract.owner != StageSubsystem::kAnalysis && contract.produced_outputs == StageDataDomain::kNone) {
    throw std::invalid_argument(prefix + "must declare produced_outputs");
  }
  if (contract.owner != StageSubsystem::kAnalysis && contract.stage != IntegrationStage::kOutputCheck &&
      contract.mutated_state == StageDataDomain::kNone) {
    throw std::invalid_argument(prefix + "must declare mutated_state");
  }
  if (contract.owner != StageSubsystem::kAnalysis && contract.stage != IntegrationStage::kOutputCheck &&
      contract.active_set_family == StageActiveSetFamily::kNone) {
    throw std::invalid_argument(prefix + "must declare active_set_family");
  }
  if (contract.stage == IntegrationStage::kOutputCheck) {
    if (contract.sync_requirements != StageSyncRequirement::kGlobal) {
      throw std::invalid_argument(prefix + "must require a global output boundary");
    }
    if (contract.active_set_family != StageActiveSetFamily::kOutputState) {
      throw std::invalid_argument(prefix + "must use the output-state active-set family");
    }
    if (contract.restart_safety != StageSafety::kSafe || contract.output_safety != StageSafety::kSafe) {
      throw std::invalid_argument(prefix + "must be restart-safe and output-safe");
    }
  }
  if (contract.owner == StageSubsystem::kGravity && contract.stage == IntegrationStage::kForceRefresh &&
      contract.sync_requirements != StageSyncRequirement::kForceEvaluation) {
    throw std::invalid_argument(prefix + "must declare force-evaluation synchronization");
  }
}

void validateRuntimeStageContract(
    const StageContract& contract,
    const StepContext& context,
    bool require_output_safe_boundary) {
  if (contract.stage != context.stage) {
    throw std::runtime_error("stage contract metadata mismatch for dispatched callback");
  }
  if (require_output_safe_boundary &&
      (contract.output_safety != StageSafety::kSafe || contract.restart_safety != StageSafety::kSafe)) {
    throw std::runtime_error("output boundary callback must declare output-safe and restart-safe contract");
  }
  if (contract.stage == IntegrationStage::kOutputCheck) {
    if (!context.boundary.output_safe || !context.boundary.restart_safe) {
      throw std::runtime_error("output boundary stage reached an unsafe integration boundary");
    }
    if (context.integrator_state.inside_kdk_step) {
      throw std::runtime_error("output boundary stage reached an open KDK step");
    }
  }
  if (contract.sync_requirements == StageSyncRequirement::kPmRefreshBoundary &&
      (!context.boundary.pm_refresh_allowed || context.boundary.local_substep)) {
    throw std::runtime_error("PM-refresh-boundary contract reached an illegal local integration step");
  }
}

void materializeLocalPmRefreshDirective(StepContext& context) {
  auto& directive = context.pm_refresh_directive;
  if (!directive.sync_event_requested) {
    return;
  }
  if (!directive.cadence_opportunity_allowed) {
    if (directive.reason == PmRefreshDirective::Reason::kSourceMutationForceRefresh) {
      throw std::runtime_error(
          "gravity source state changed after the scheduled PM refresh on a boundary where long-range repair is illegal");
    }
    throw std::runtime_error(
        "PM cadence event requested outside an integrator-authorized refresh opportunity");
  }
  materializePmRefreshDirective(
      context, context.integrator_state.pm_long_range_field_valid);
}

void validatePmRefreshDirectiveLegality(const StepContext& context, std::string_view phase) {
  const auto& directive = context.pm_refresh_directive;
  const bool initial_surface_exception =
      directive.reason == PmRefreshDirective::Reason::kInitialForceBootstrap &&
      directive.sync_stage == PmSyncStage::kInitialLongRangeBootstrap &&
      context.integrator_state.step_index == 0U &&
      directive.initial_cache_bootstrap_allowed;
  const bool coordinated_local_force_refresh =
      context.boundary.local_substep &&
      context.active_set.has_global_synchronization_metadata &&
      directive.reason == PmRefreshDirective::Reason::kScheduledForceRefreshStage &&
      directive.sync_stage == PmSyncStage::kScheduledLongRangeRefresh;
  if (directive.refresh_long_range_field && !directive.has_sync_event) {
    throw std::runtime_error("PM refresh directive requested long-range refresh without a sync event at " + std::string(phase));
  }
  if (directive.has_sync_event) {
    if (directive.sync_stage == PmSyncStage::kNone) {
      throw std::runtime_error("PM refresh directive sync event has no typed PM sync stage at " + std::string(phase));
    }
    if (!directive.cadence_opportunity_allowed) {
      throw std::runtime_error("PM sync event reached a non-cadence opportunity at " + std::string(phase));
    }
    if (!context.boundary.pm_refresh_allowed ||
        (context.boundary.local_substep && !initial_surface_exception &&
         !coordinated_local_force_refresh)) {
      throw std::runtime_error("PM sync event reached an illegal local integration boundary at " + std::string(phase));
    }
    if (!std::isfinite(directive.force_evaluation_scale_factor) || directive.force_evaluation_scale_factor <= 0.0) {
      throw std::runtime_error("PM sync event has invalid force-evaluation scale factor at " + std::string(phase));
    }
  } else if (directive.sync_stage != PmSyncStage::kNone) {
    throw std::runtime_error("typed PM sync stage was set without a sync event at " + std::string(phase));
  }
  if (context.boundary.local_substep && !initial_surface_exception &&
      !coordinated_local_force_refresh &&
      (directive.has_sync_event || directive.refresh_long_range_field || directive.cadence_opportunity_allowed)) {
    throw std::runtime_error("local active-bin step attempted to authorize PM long-range synchronization at " + std::string(phase));
  }
}

[[nodiscard]] double redshiftFromScaleFactor(double scale_factor) {
  if (scale_factor <= 0.0) {
    throw std::invalid_argument("scale_factor must be positive");
  }
  return 1.0 / scale_factor - 1.0;
}

[[nodiscard]] double midpointIntegrateCosmicTime(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples) {
  if (scale_factor_begin <= 0.0 || scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale factors must be positive for cosmic-time integral");
  }
  if (scale_factor_end < scale_factor_begin) {
    throw std::invalid_argument("scale_factor_end must be >= scale_factor_begin");
  }
  const std::uint32_t samples = std::max<std::uint32_t>(midpoint_samples, 1U);
  const double delta_a = (scale_factor_end - scale_factor_begin) / static_cast<double>(samples);
  if (delta_a == 0.0) {
    return 0.0;
  }

  double accum = 0.0;
  for (std::uint32_t i = 0; i < samples; ++i) {
    const double a_mid = scale_factor_begin + (static_cast<double>(i) + 0.5) * delta_a;
    accum += 1.0 / (a_mid * background.hubbleSi(a_mid));
  }
  return accum * delta_a;
}

[[nodiscard]] std::size_t timestepSourceIndex(TimeStepCandidateSource source) noexcept {
  return static_cast<std::size_t>(source);
}

[[nodiscard]] double midpointIntegrateDriftLike(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples) {
  if (scale_factor_begin <= 0.0 || scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale factors must be positive for drift-like integrals");
  }
  if (scale_factor_end < scale_factor_begin) {
    throw std::invalid_argument("scale_factor_end must be >= scale_factor_begin");
  }

  const std::uint32_t samples = std::max<std::uint32_t>(midpoint_samples, 1U);
  const double delta_a = (scale_factor_end - scale_factor_begin) / static_cast<double>(samples);
  if (delta_a == 0.0) {
    return 0.0;
  }

  double accum = 0.0;
  for (std::uint32_t i = 0; i < samples; ++i) {
    const double a_mid = scale_factor_begin + (static_cast<double>(i) + 0.5) * delta_a;
    accum += 1.0 / (a_mid * a_mid * background.hubbleSi(a_mid));
  }
  return accum * delta_a;
}

[[nodiscard]] std::uint64_t powerOfTwo(std::uint8_t exponent) {
  if (exponent >= 63) {
    throw std::invalid_argument("bin exponent too large for 64-bit tick period");
  }
  return 1ULL << exponent;
}

[[nodiscard]] double requireValidTimeStepCriterion(
    double value,
    std::string_view criterion_name,
    std::uint32_t element_index) {
  if (std::isinf(value) && value > 0.0) {
    return value;
  }
  if (!std::isfinite(value) || value <= 0.0) {
    std::ostringstream out;
    out << "invalid timestep criterion"
        << "; criterion=" << criterion_name
        << "; element=" << element_index
        << "; value=" << value
        << "; expected finite positive value or +infinity for explicit no-limit";
    throw std::runtime_error(out.str());
  }
  return value;
}

[[nodiscard]] std::string timeBinContextMessage(
    std::string_view invariant,
    std::string_view source_label,
    std::uint64_t current_tick,
    std::uint32_t element_index,
    std::uint8_t bin,
    std::uint64_t next_activation_tick,
    std::uint64_t integrator_step = 0,
    double integrator_time = 0.0,
    std::string_view pm_sync_state = "n/a") {
  std::ostringstream out;
  out << "time-bin invariant failure: " << invariant
      << "; source=" << source_label
      << "; element=" << element_index
      << "; current_tick=" << current_tick
      << "; bin=" << static_cast<unsigned>(bin)
      << "; next_activation_tick=" << next_activation_tick
      << "; integrator_step=" << integrator_step
      << "; integrator_time=" << integrator_time
      << "; pm_sync=" << pm_sync_state;
  return out.str();
}

[[nodiscard]] std::string schedulerContextMessage(
    std::string_view invariant,
    std::string_view source_label,
    std::uint64_t current_tick,
    std::uint64_t integrator_step = 0,
    double integrator_time = 0.0,
    std::string_view pm_sync_state = "n/a") {
  std::ostringstream out;
  out << "scheduler invariant failure: " << invariant
      << "; source=" << source_label
      << "; current_tick=" << current_tick
      << "; integrator_step=" << integrator_step
      << "; integrator_time=" << integrator_time
      << "; pm_sync=" << pm_sync_state;
  return out.str();
}

[[nodiscard]] std::string pmContextMessage(
    std::string_view invariant,
    std::string_view source_label,
    const PmSyncEvent& event,
    std::uint64_t cadence_steps,
    std::uint64_t current_field_version) {
  std::ostringstream out;
  out << "PM synchronization invariant failure: " << invariant
      << "; source=" << source_label
      << "; opportunity=" << event.gravity_kick_opportunity
      << "; cadence_steps=" << cadence_steps
      << "; refresh=" << (event.refresh_long_range_field ? "true" : "false")
      << "; event_field_version=" << event.field_version
      << "; current_field_version=" << current_field_version
      << "; last_refresh_opportunity=" << event.last_refresh_opportunity
      << "; field_built_step=" << event.field_built_step_index
      << "; field_built_scale_factor=" << event.field_built_scale_factor;
  return out.str();
}

}  // namespace

std::string_view timeStepCandidateSourceName(TimeStepCandidateSource source) {
  switch (source) {
    case TimeStepCandidateSource::kHydroCfl:
      return "hydro_cfl";
    case TimeStepCandidateSource::kGravityAcceleration:
      return "gravity_acceleration";
    case TimeStepCandidateSource::kCosmologyExpansion:
      return "cosmology_expansion";
    case TimeStepCandidateSource::kSourceTerm:
      return "source_term";
    case TimeStepCandidateSource::kUserClamp:
      return "user_clamp";
  }
  return "unknown";
}

std::string_view stepBoundaryKindName(StepBoundaryKind kind) {
  switch (kind) {
    case StepBoundaryKind::kLocalActiveBinStep:
      return "local_active_bin_step";
    case StepBoundaryKind::kGlobalSynchronizationPoint:
      return "global_synchronization_point";
    case StepBoundaryKind::kPmRefreshPoint:
      return "pm_refresh_point";
    case StepBoundaryKind::kSnapshotPoint:
      return "snapshot_point";
    case StepBoundaryKind::kCheckpointPoint:
      return "checkpoint_point";
  }
  return "unknown";
}

bool isRestartSafeBoundary(StepBoundaryKind kind) noexcept {
  return kind == StepBoundaryKind::kGlobalSynchronizationPoint ||
         kind == StepBoundaryKind::kSnapshotPoint ||
         kind == StepBoundaryKind::kCheckpointPoint;
}

bool isOutputSafeBoundary(StepBoundaryKind kind) noexcept {
  return kind == StepBoundaryKind::kGlobalSynchronizationPoint ||
         kind == StepBoundaryKind::kSnapshotPoint ||
         kind == StepBoundaryKind::kCheckpointPoint;
}

CosmologicalTimeline::CosmologicalTimeline(const LambdaCdmBackground* background, double time_si_per_code)
    : m_background(background), m_time_si_per_code(time_si_per_code) {
  if (m_time_si_per_code <= 0.0 || !std::isfinite(m_time_si_per_code)) {
    throw std::invalid_argument("CosmologicalTimeline: time_si_per_code must be finite and positive");
  }
}

CosmologicalStepFactors CosmologicalTimeline::prepareStep(
    double current_time_code,
    double current_scale_factor,
    double dt_time_code) const {
  if (dt_time_code <= 0.0) {
    throw std::invalid_argument("dt_time_code must be positive");
  }
  if (current_scale_factor <= 0.0) {
    throw std::invalid_argument("current_scale_factor must be positive");
  }

  CosmologicalStepFactors step{};
  step.cosmological = m_background != nullptr;
  step.time_begin_code = current_time_code;
  step.time_end_code = current_time_code + dt_time_code;
  step.dt_time_code = dt_time_code;
  step.time_si_per_code = m_time_si_per_code;
  step.dt_time_si = dt_time_code * m_time_si_per_code;
  step.scale_factor_begin = current_scale_factor;
  step.redshift_begin = redshiftFromScaleFactor(current_scale_factor);

  if (m_background == nullptr) {
    step.scale_factor_midpoint = current_scale_factor;
    step.scale_factor_end = current_scale_factor;
    step.redshift_end = step.redshift_begin;
    step.drift_factor_code = dt_time_code;
    step.first_kick_factor_code = 0.5 * dt_time_code;
    step.second_kick_factor_code = 0.5 * dt_time_code;
    step.first_hubble_drag_factor = 1.0;
    step.second_hubble_drag_factor = 1.0;
    step.hubble_drag_factor = 1.0;
    return step;
  }

  step.hubble_begin_code = m_background->hubbleSi(current_scale_factor) * m_time_si_per_code;
  step.scale_factor_midpoint = advanceScaleFactorByCosmicTime(*m_background, current_scale_factor, 0.5 * step.dt_time_si);
  step.scale_factor_end = advanceScaleFactorByCosmicTime(*m_background, current_scale_factor, step.dt_time_si);
  step.redshift_end = redshiftFromScaleFactor(step.scale_factor_end);
  step.hubble_end_code = m_background->hubbleSi(step.scale_factor_end) * m_time_si_per_code;
  step.drift_factor_code = computeComovingDriftFactor(*m_background, step.scale_factor_begin, step.scale_factor_end, 64) / m_time_si_per_code;
  step.first_kick_factor_code = computeComovingKickFactor(*m_background, step.scale_factor_begin, step.scale_factor_midpoint, 64) / m_time_si_per_code;
  step.second_kick_factor_code = computeComovingKickFactor(*m_background, step.scale_factor_midpoint, step.scale_factor_end, 64) / m_time_si_per_code;
  step.first_hubble_drag_factor = computeHubbleDragFactor(step.scale_factor_begin, step.scale_factor_midpoint);
  step.second_hubble_drag_factor = computeHubbleDragFactor(step.scale_factor_midpoint, step.scale_factor_end);
  step.hubble_drag_factor = computeHubbleDragFactor(step.scale_factor_begin, step.scale_factor_end);
  return step;
}

void CosmologicalTimeline::commitStep(IntegratorState& integrator_state, const CosmologicalStepFactors& step) const {
  if (integrator_state.inside_kdk_step) {
    throw std::runtime_error("CosmologicalTimeline::commitStep called while integrator is marked inside a KDK step");
  }
  integrator_state.current_time_code = step.time_end_code;
  integrator_state.current_scale_factor = step.scale_factor_end;
  integrator_state.current_redshift = step.redshift_end;
  integrator_state.current_hubble_rate_code = step.hubble_end_code;
  integrator_state.time_si_per_code = step.time_si_per_code;
  integrator_state.last_drift_factor_code = step.drift_factor_code;
  integrator_state.last_first_kick_factor_code = step.first_kick_factor_code;
  integrator_state.last_second_kick_factor_code = step.second_kick_factor_code;
  integrator_state.last_first_hubble_drag_factor = step.first_hubble_drag_factor;
  integrator_state.last_second_hubble_drag_factor = step.second_hubble_drag_factor;
  ++integrator_state.step_index;
}

void PmSynchronizationState::reset(std::uint64_t cadence_steps) {
  if (cadence_steps == 0) {
    throw std::invalid_argument("PM synchronization cadence must be >= 1");
  }
  m_cadence_steps = cadence_steps;
  m_gravity_kick_opportunity = 0;
  m_last_refresh_opportunity = 0;
  m_field_version = 0;
  m_last_refresh_step_index = 0;
  m_last_refresh_scale_factor = 1.0;
  m_refresh_commit_pending = false;
  m_pending_refresh_opportunity = 0;
  m_pending_refresh_field_version = 0;
}

PmSyncEvent PmSynchronizationState::previewKickOpportunity(
    std::uint64_t step_index,
    double scale_factor,
    bool has_long_range_field) const {
  if (m_refresh_commit_pending) {
    const PmSyncEvent pending{
        .gravity_kick_opportunity = m_pending_refresh_opportunity,
        .refresh_long_range_field = true,
        .field_version = m_pending_refresh_field_version,
        .last_refresh_opportunity = m_pending_refresh_opportunity,
        .field_built_step_index = m_last_refresh_step_index,
        .field_built_scale_factor = m_last_refresh_scale_factor,
    };
    throw std::runtime_error(pmContextMessage(
        "previous refresh event was not committed before next kick opportunity",
        "PmSynchronizationState::previewKickOpportunity",
        pending,
        m_cadence_steps,
        m_field_version));
  }
  if (m_gravity_kick_opportunity == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("PM synchronization gravity kick opportunity overflow");
  }
  const std::uint64_t next_opportunity = m_gravity_kick_opportunity + 1U;
  const bool refresh = !has_long_range_field ||
      ((next_opportunity - m_last_refresh_opportunity) >= m_cadence_steps);
  if (refresh && m_field_version == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("PM synchronization field version overflow");
  }
  return PmSyncEvent{
      .gravity_kick_opportunity = next_opportunity,
      .refresh_long_range_field = refresh,
      .field_version = refresh ? (m_field_version + 1U) : m_field_version,
      .last_refresh_opportunity = refresh ? next_opportunity : m_last_refresh_opportunity,
      .field_built_step_index = refresh ? step_index : m_last_refresh_step_index,
      .field_built_scale_factor = refresh ? scale_factor : m_last_refresh_scale_factor,
  };
}

void PmSynchronizationState::commitKickOpportunity(const PmSyncEvent& event) {
  if (m_refresh_commit_pending) {
    throw std::runtime_error(
        "PmSynchronizationState::commitKickOpportunity: previous refresh event is still pending");
  }
  if (m_gravity_kick_opportunity == std::numeric_limits<std::uint64_t>::max() ||
      event.gravity_kick_opportunity != m_gravity_kick_opportunity + 1U) {
    throw std::runtime_error(pmContextMessage(
        "kick opportunity does not immediately follow committed cadence state",
        "PmSynchronizationState::commitKickOpportunity",
        event,
        m_cadence_steps,
        m_field_version));
  }
  if (event.refresh_long_range_field) {
    if (m_field_version == std::numeric_limits<std::uint64_t>::max() ||
        event.field_version != m_field_version + 1U ||
        event.last_refresh_opportunity != event.gravity_kick_opportunity ||
        !std::isfinite(event.field_built_scale_factor) ||
        event.field_built_scale_factor <= 0.0) {
      throw std::runtime_error(pmContextMessage(
          "refresh opportunity metadata is inconsistent with committed cadence state",
          "PmSynchronizationState::commitKickOpportunity",
          event,
          m_cadence_steps,
          m_field_version));
    }
  } else if (event.field_version != m_field_version ||
             event.last_refresh_opportunity != m_last_refresh_opportunity ||
             event.field_built_step_index != m_last_refresh_step_index ||
             event.field_built_scale_factor != m_last_refresh_scale_factor) {
    throw std::runtime_error(pmContextMessage(
        "reuse opportunity metadata is inconsistent with committed cadence state",
        "PmSynchronizationState::commitKickOpportunity",
        event,
        m_cadence_steps,
        m_field_version));
  }

  m_gravity_kick_opportunity = event.gravity_kick_opportunity;
  if (event.refresh_long_range_field) {
    m_refresh_commit_pending = true;
    m_pending_refresh_opportunity = event.gravity_kick_opportunity;
    m_pending_refresh_field_version = event.field_version;
  }
}

PmSyncEvent PmSynchronizationState::registerKickOpportunity(
    std::uint64_t step_index,
    double scale_factor,
    bool has_long_range_field) {
  const PmSyncEvent event = previewKickOpportunity(
      step_index, scale_factor, has_long_range_field);
  commitKickOpportunity(event);
  return event;
}

void PmSynchronizationState::commitRefresh(const PmSyncEvent& event) {
  if (!event.refresh_long_range_field) {
    return;
  }
  if (!m_refresh_commit_pending || event.gravity_kick_opportunity != m_pending_refresh_opportunity) {
    throw std::runtime_error(pmContextMessage(
        "refresh event does not match pending kick opportunity",
        "PmSynchronizationState::commitRefresh",
        event,
        m_cadence_steps,
        m_field_version));
  }
  if (event.gravity_kick_opportunity != m_gravity_kick_opportunity) {
    throw std::runtime_error(pmContextMessage(
        "event does not match current kick opportunity",
        "PmSynchronizationState::commitRefresh",
        event,
        m_cadence_steps,
        m_field_version));
  }
  if (event.field_version != m_field_version + 1U || event.field_version != m_pending_refresh_field_version) {
    throw std::runtime_error(pmContextMessage(
        "field version transition is illegal",
        "PmSynchronizationState::commitRefresh",
        event,
        m_cadence_steps,
        m_field_version));
  }
  m_field_version = event.field_version;
  m_last_refresh_opportunity = event.last_refresh_opportunity;
  m_last_refresh_step_index = event.field_built_step_index;
  m_last_refresh_scale_factor = event.field_built_scale_factor;
  m_refresh_commit_pending = false;
  m_pending_refresh_opportunity = 0;
  m_pending_refresh_field_version = 0;
}


PmSynchronizationPersistentState PmSynchronizationState::exportPersistentState() const {
  return PmSynchronizationPersistentState{
      .cadence_steps = m_cadence_steps,
      .gravity_kick_opportunity = m_gravity_kick_opportunity,
      .last_refresh_opportunity = m_last_refresh_opportunity,
      .field_version = m_field_version,
      .last_refresh_step_index = m_last_refresh_step_index,
      .last_refresh_scale_factor = m_last_refresh_scale_factor,
      .refresh_commit_pending = m_refresh_commit_pending,
      .pending_refresh_opportunity = m_pending_refresh_opportunity,
      .pending_refresh_field_version = m_pending_refresh_field_version,
  };
}

void PmSynchronizationState::importPersistentState(const PmSynchronizationPersistentState& persistent_state) {
  if (persistent_state.cadence_steps == 0) {
    throw std::invalid_argument("PmSynchronizationState::importPersistentState: cadence_steps must be >= 1");
  }
  if (persistent_state.last_refresh_opportunity > persistent_state.gravity_kick_opportunity) {
    throw std::invalid_argument(
        "PmSynchronizationState::importPersistentState: last refresh opportunity is ahead of current kick opportunity");
  }
  if (persistent_state.refresh_commit_pending) {
    if (persistent_state.pending_refresh_opportunity == 0 ||
        persistent_state.pending_refresh_opportunity != persistent_state.gravity_kick_opportunity ||
        persistent_state.pending_refresh_field_version != persistent_state.field_version + 1U) {
      throw std::invalid_argument(
          "PmSynchronizationState::importPersistentState: pending refresh metadata is inconsistent");
    }
  } else if (persistent_state.pending_refresh_opportunity != 0 ||
             persistent_state.pending_refresh_field_version != 0) {
    throw std::invalid_argument(
        "PmSynchronizationState::importPersistentState: nonzero pending refresh metadata without pending refresh");
  }

  m_cadence_steps = persistent_state.cadence_steps;
  m_gravity_kick_opportunity = persistent_state.gravity_kick_opportunity;
  m_last_refresh_opportunity = persistent_state.last_refresh_opportunity;
  m_field_version = persistent_state.field_version;
  m_last_refresh_step_index = persistent_state.last_refresh_step_index;
  m_last_refresh_scale_factor = persistent_state.last_refresh_scale_factor;
  m_refresh_commit_pending = persistent_state.refresh_commit_pending;
  m_pending_refresh_opportunity = persistent_state.pending_refresh_opportunity;
  m_pending_refresh_field_version = persistent_state.pending_refresh_field_version;
}


void materializePmRefreshDirective(
    StepContext& context,
    bool has_long_range_field) {
  auto& directive = context.pm_refresh_directive;
  if (!directive.sync_event_requested) {
    throw std::logic_error(
        "PM refresh directive materialization requires a requested sync event");
  }
  if (directive.has_sync_event) {
    throw std::logic_error(
        "PM refresh directive sync event was materialized more than once");
  }
  if (!directive.cadence_opportunity_allowed) {
    throw std::runtime_error(
        "PM refresh directive materialization reached a non-cadence opportunity");
  }
  if (!std::isfinite(directive.force_evaluation_scale_factor) ||
      directive.force_evaluation_scale_factor <= 0.0) {
    throw std::runtime_error(
        "PM refresh directive materialization requires a finite positive force-evaluation scale factor");
  }

  bool field_available_for_cadence = has_long_range_field;
  switch (directive.reason) {
    case PmRefreshDirective::Reason::kInitialForceBootstrap:
      directive.sync_stage = PmSyncStage::kInitialLongRangeBootstrap;
      field_available_for_cadence = false;
      break;
    case PmRefreshDirective::Reason::kScheduledForceRefreshStage:
      directive.sync_stage = PmSyncStage::kScheduledLongRangeRefresh;
      break;
    case PmRefreshDirective::Reason::kSourceMutationForceRefresh:
      directive.sync_stage = PmSyncStage::kScheduledLongRangeRefresh;
      field_available_for_cadence = false;
      break;
    case PmRefreshDirective::Reason::kNone:
      throw std::logic_error(
          "PM refresh directive materialization requires a typed refresh reason");
  }

  const PmSyncEvent event = context.integrator_state.pm_sync_state.previewKickOpportunity(
      context.integrator_state.step_index,
      directive.force_evaluation_scale_factor,
      field_available_for_cadence);
  if ((directive.reason == PmRefreshDirective::Reason::kInitialForceBootstrap ||
       directive.reason == PmRefreshDirective::Reason::kSourceMutationForceRefresh) &&
      !event.refresh_long_range_field) {
    throw std::logic_error(
        "mandatory PM invalidation materialized as a cadence reuse event");
  }

  directive.has_sync_event = true;
  directive.refresh_long_range_field = event.refresh_long_range_field;
  directive.gravity_kick_opportunity = event.gravity_kick_opportunity;
  directive.field_version = event.field_version;
  directive.last_refresh_opportunity = event.last_refresh_opportunity;
  directive.field_built_step_index = event.field_built_step_index;
  directive.field_built_scale_factor = event.field_built_scale_factor;
}

std::string_view integrationStageName(IntegrationStage stage) {
  switch (stage) {
    case IntegrationStage::kGravityKickPre:
      return "gravity_kick_pre";
    case IntegrationStage::kDrift:
      return "drift";
    case IntegrationStage::kForceRefresh:
      return "force_refresh";
    case IntegrationStage::kHydroUpdate:
      return "hydro_update";
    case IntegrationStage::kSourceTerms:
      return "source_terms";
    case IntegrationStage::kGravityKickPost:
      return "gravity_kick_post";
    case IntegrationStage::kAnalysisHooks:
      return "analysis_hooks";
    case IntegrationStage::kOutputCheck:
      return "output_check";
  }
  return "unknown";
}

bool ActiveSetDescriptor::hasParticleSubset(std::size_t total_particle_count) const noexcept {
  return particles_are_subset && particle_indices.size() < total_particle_count;
}

bool ActiveSetDescriptor::hasCellSubset(std::size_t total_cell_count) const noexcept {
  return cells_are_subset && cell_indices.size() < total_cell_count;
}

ActiveSetDescriptor makeSchedulerActiveSetDescriptor(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices) {
  ActiveSetDescriptor descriptor{
      .particle_indices = active_particle_indices,
      .cell_indices = active_cell_indices,
      .particles_are_subset = active_particle_indices.size() != state.particles.size(),
      .cells_are_subset = active_cell_indices.size() != state.cells.size(),
      .particles_from_scheduler = true,
      .cells_from_scheduler = !active_cell_indices.empty(),
      .has_generation_metadata = true,
      .source_particle_index_generation = state.particleIndexGeneration(),
      .source_cell_index_generation = state.cellIndexGeneration(),
      .source_scheduler_tick = scheduler.currentTick(),
  };
  debugAssertActiveSetDescriptorFresh(descriptor, state, scheduler);
  return descriptor;
}

ParticleReorderMap buildParticleReorderMapByScheduler(
    const SimulationState& state,
    const HierarchicalTimeBinScheduler& scheduler) {
  if (scheduler.elementCount() != state.particles.size()) {
    throw std::invalid_argument(
        "buildParticleReorderMapByScheduler: scheduler element count must match particle count");
  }
  ParticleReorderMap reorder_map;
  reorder_map.new_to_old_index.resize(state.particles.size());
  std::iota(reorder_map.new_to_old_index.begin(), reorder_map.new_to_old_index.end(), 0U);
  std::stable_sort(
      reorder_map.new_to_old_index.begin(),
      reorder_map.new_to_old_index.end(),
      [&](std::uint32_t lhs, std::uint32_t rhs) {
        return std::tuple{scheduler.binIndex(lhs), lhs} < std::tuple{scheduler.binIndex(rhs), rhs};
      });
  reorder_map.old_to_new_index.resize(state.particles.size());
  for (std::size_t new_index = 0; new_index < reorder_map.new_to_old_index.size(); ++new_index) {
    reorder_map.old_to_new_index[reorder_map.new_to_old_index[new_index]] =
        checkedLocalParticleRow(new_index, "buildParticleReorderMapByScheduler new row");
  }
  return reorder_map;
}

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state) {
  if (active_set.has_generation_metadata &&
      active_set.source_particle_index_generation != state.particleIndexGeneration()) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor particle generation is stale",
        "debugAssertActiveSetDescriptorFresh",
        active_set.source_scheduler_tick));
  }
  if (active_set.has_generation_metadata &&
      active_set.source_cell_index_generation != state.cellIndexGeneration()) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor cell generation is stale",
        "debugAssertActiveSetDescriptorFresh",
        active_set.source_scheduler_tick));
  }
  for (const std::uint32_t particle_index : active_set.particle_indices) {
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range(timeBinContextMessage(
          "ActiveSetDescriptor contains stale particle index",
          "debugAssertActiveSetDescriptorFresh",
          active_set.source_scheduler_tick,
          particle_index,
          0,
          active_set.source_scheduler_tick));
    }
  }
  for (const std::uint32_t cell_index : active_set.cell_indices) {
    if (cell_index >= state.cells.size()) {
      throw std::out_of_range(timeBinContextMessage(
          "ActiveSetDescriptor contains stale cell index",
          "debugAssertActiveSetDescriptorFresh",
          active_set.source_scheduler_tick,
          cell_index,
          0,
          active_set.source_scheduler_tick));
    }
  }
}

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    std::uint64_t expected_scheduler_tick) {
  debugAssertActiveSetDescriptorFresh(active_set, state);
  if (!active_set.has_generation_metadata || !active_set.particles_from_scheduler) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor must carry scheduler provenance before solver callbacks",
        "debugAssertActiveSetDescriptorFresh",
        expected_scheduler_tick));
  }
  if ((active_set.particles_from_scheduler || active_set.cells_from_scheduler) &&
      active_set.source_scheduler_tick != expected_scheduler_tick) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor scheduler tick is stale",
        "debugAssertActiveSetDescriptorFresh",
        expected_scheduler_tick));
  }
}

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    const HierarchicalTimeBinScheduler& scheduler) {
  debugAssertActiveSetDescriptorFresh(active_set, state, scheduler.currentTick());
  if (active_set.particles_from_scheduler) {
    const auto active = scheduler.activeElements();
    if (active.size() != active_set.particle_indices.size() ||
        !std::equal(active.begin(), active.end(), active_set.particle_indices.begin())) {
      throw std::runtime_error(schedulerContextMessage(
          "ActiveSetDescriptor particle indices do not match scheduler active set",
          "debugAssertActiveSetDescriptorFresh",
          scheduler.currentTick()));
    }
  }
}

std::vector<IntegrationStage> StageScheduler::schedule(
    const IntegratorState& integrator_state,
    const ActiveSetDescriptor& /*active_set*/) const {
  if (integrator_state.scheme != TimeStepScheme::kKickDriftKick) {
    throw std::invalid_argument("unsupported timestep scheme");
  }

  return std::vector<IntegrationStage>(k_kick_drift_kick_order.begin(), k_kick_drift_kick_order.end());
}

std::span<const IntegrationStage> StageScheduler::kickDriftKickOrder() { return k_kick_drift_kick_order; }

bool isCanonicalIntegrationStageOrder(std::span<const IntegrationStage> ordered_stages) {
  return ordered_stages.size() == k_kick_drift_kick_order.size() &&
         std::equal(
             ordered_stages.begin(),
             ordered_stages.end(),
             k_kick_drift_kick_order.begin(),
             k_kick_drift_kick_order.end());
}

StepOrchestrator::StepOrchestrator(StageScheduler scheduler) : m_scheduler(std::move(scheduler)) {}

void StepOrchestrator::registerCallback(IntegrationCallback& callback) {
  const auto stages = callback.integrationStages();
  const auto contracts = callback.stageContracts();
  if (stages.empty()) {
    throw std::invalid_argument("IntegrationCallback '" + std::string(callback.callbackName()) + "' declares no stages");
  }
  if (contracts.empty()) {
    throw std::invalid_argument("IntegrationCallback '" + std::string(callback.callbackName()) + "' declares no stage contracts");
  }
  if (contracts.size() != stages.size()) {
    throw std::invalid_argument(
        "IntegrationCallback '" + std::string(callback.callbackName()) +
        "' must declare exactly one executable contract for each registered stage");
  }
  for (std::size_t i = 0; i < contracts.size(); ++i) {
    const auto stage = contracts[i].stage;
    const bool stage_registered = std::find(stages.begin(), stages.end(), stage) != stages.end();
    if (!stage_registered) {
      throw std::invalid_argument(
          "IntegrationCallback '" + std::string(callback.callbackName()) +
          "' declares a stage contract for an unregistered stage");
    }
    if (std::find_if(contracts.begin(), contracts.begin() + static_cast<std::ptrdiff_t>(i),
                     [stage](const StageContract& existing) { return existing.stage == stage; }) !=
        contracts.begin() + static_cast<std::ptrdiff_t>(i)) {
      throw std::invalid_argument(
          "IntegrationCallback '" + std::string(callback.callbackName()) +
          "' declares duplicate executable contracts for the same stage");
    }
    validateStageContractShape(callback.callbackName(), contracts[i]);
  }
  for (const IntegrationStage stage : stages) {
    const std::size_t index = integrationStageIndex(stage);
    if (index >= m_handlers_by_stage.size()) {
      throw std::out_of_range("IntegrationCallback '" + std::string(callback.callbackName()) + "' declares an invalid stage");
    }
    auto& handlers = m_handlers_by_stage[index];
    auto& stage_contracts = m_contracts_by_stage[index];
    if (std::find(handlers.begin(), handlers.end(), &callback) != handlers.end()) {
      throw std::invalid_argument(
          "IntegrationCallback '" + std::string(callback.callbackName()) + "' declares a duplicate stage");
    }
    handlers.push_back(&callback);
    const auto contract_it = std::find_if(contracts.begin(), contracts.end(), [stage](const StageContract& contract) {
      return contract.stage == stage;
    });
    if (contract_it == contracts.end()) {
      throw std::invalid_argument(
          "IntegrationCallback '" + std::string(callback.callbackName()) + "' is missing a contract for stage");
    }
    stage_contracts.push_back(*contract_it);
  }
  ++m_callback_count;
}

std::size_t StepOrchestrator::callbackCount() const noexcept { return m_callback_count; }

std::span<IntegrationCallback* const> StepOrchestrator::handlersFor(IntegrationStage stage) const noexcept {
  const std::size_t index = integrationStageIndex(stage);
  if (index >= m_handlers_by_stage.size()) {
    return {};
  }
  const auto& handlers = m_handlers_by_stage[index];
  return {handlers.data(), handlers.size()};
}
std::span<const StageContract> StepOrchestrator::contractsFor(IntegrationStage stage) const noexcept {
  const std::size_t index = integrationStageIndex(stage);
  if (index >= m_contracts_by_stage.size()) {
    return {};
  }
  const auto& contracts = m_contracts_by_stage[index];
  return {contracts.data(), contracts.size()};
}
std::optional<StageContract> StepOrchestrator::contractForHandlerStage(
    const IntegrationCallback& callback,
    IntegrationStage stage) const noexcept {
  const auto handlers = handlersFor(stage);
  const auto contracts = contractsFor(stage);
  for (std::size_t i = 0; i < handlers.size() && i < contracts.size(); ++i) {
    if (handlers[i] == &callback) {
      return contracts[i];
    }
  }
  return std::nullopt;
}

void StepOrchestrator::dispatchStageHandlers(
    StepContext& context,
    bool require_output_safe_boundary) const {
  materializeLocalPmRefreshDirective(context);
  const std::string stage_name = "stage." + std::string(integrationStageName(context.stage));
  COSMOSIM_PROFILE_SCOPE(context.profiler_session, stage_name);
  if (context.profiler_session != nullptr) {
    context.profiler_session->counters().addCount(stage_name + ".invocations", 1);
  }

  for (auto* callback : handlersFor(context.stage)) {
    const auto contract = contractForHandlerStage(*callback, context.stage);
    if (!contract.has_value()) {
      throw std::runtime_error("registered callback missing executable stage contract");
    }
    validateRuntimeStageContract(*contract, context, require_output_safe_boundary);
    validatePmRefreshDirectiveLegality(context, std::string("before callback ") + std::string(callback->callbackName()));
    const std::string callback_phase = "callback." + std::string(callback->callbackName());
    COSMOSIM_PROFILE_SCOPE(context.profiler_session, callback_phase);
    callback->onStage(context);
    validatePmRefreshDirectiveLegality(context, std::string("after callback ") + std::string(callback->callbackName()));
    if (context.profiler_session != nullptr) {
      context.profiler_session->counters().addCount(callback_phase + ".invocations", 1);
    }
  }
}

void StepOrchestrator::executeOutputBoundary(
    SimulationState& state,
    IntegratorState& integrator_state,
    ProfilerSession* profiler_session,
    StepBoundaryKind requested_boundary_kind) const {
  executeOutputBoundaryWithDispatcher(
      state,
      integrator_state,
      [this](StepContext& context, bool require_output_safe_boundary) {
        dispatchStageHandlers(context, require_output_safe_boundary);
      },
      profiler_session,
      requested_boundary_kind);
}

void StepOrchestrator::executeOutputBoundaryWithDispatcher(
    SimulationState& state,
    IntegratorState& integrator_state,
    const StageDispatchFunction& dispatcher,
    ProfilerSession* profiler_session,
    StepBoundaryKind requested_boundary_kind) const {
  if (!dispatcher) {
    throw std::invalid_argument("output boundary dispatcher must be callable");
  }
  if (integrator_state.inside_kdk_step) {
    throw std::runtime_error("output boundary dispatch requested inside an open KDK step");
  }

  ActiveSetDescriptor output_active_set{};
  const StepBoundaryState boundary = classifyStepBoundary(
      state,
      output_active_set,
      false,
      requested_boundary_kind);

  StepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .active_set = output_active_set,
      .active_gravity_particles = {},
      .has_active_gravity_particles = false,
      .workspace = nullptr,
      .cosmology_background = nullptr,
      .mode_policy = nullptr,
      .profiler_session = profiler_session,
      .particle_scheduler = nullptr,
      .gas_cell_scheduler = nullptr,
      .newly_created_particle_ids = nullptr,
      .timeline_step = {},
      .boundary = boundary,
      .pm_refresh_directive = {},
      .stage = IntegrationStage::kOutputCheck,
  };

  dispatcher(context, true);
}

void StepOrchestrator::executeSingleStep(
    SimulationState& state,
    IntegratorState& integrator_state,
    ActiveSetDescriptor active_set,
    const LambdaCdmBackground* cosmology_background,
    TransientStepWorkspace* workspace,
    const ModePolicy* mode_policy,
    ProfilerSession* profiler_session,
    std::optional<std::uint64_t> expected_scheduler_tick,
    StepBoundaryKind requested_boundary_kind) const {
  executeSingleStepWithDispatcher(
      state,
      integrator_state,
      active_set,
      [this](StepContext& context, bool require_output_safe_boundary) {
        dispatchStageHandlers(context, require_output_safe_boundary);
      },
      cosmology_background,
      workspace,
      mode_policy,
      profiler_session,
      expected_scheduler_tick,
      requested_boundary_kind);
}

void StepOrchestrator::executeSingleStepWithDispatcher(
    SimulationState& state,
    IntegratorState& integrator_state,
    ActiveSetDescriptor active_set,
    const StageDispatchFunction& dispatcher,
    const LambdaCdmBackground* cosmology_background,
    TransientStepWorkspace* workspace,
    const ModePolicy* mode_policy,
    ProfilerSession* profiler_session,
    std::optional<std::uint64_t> expected_scheduler_tick,
    StepBoundaryKind requested_boundary_kind) const {
  if (!dispatcher) {
    throw std::invalid_argument("single-step dispatcher must be callable");
  }
  if (integrator_state.dt_time_code <= 0.0) {
    throw std::invalid_argument("dt_time_code must be positive");
  }
  if (expected_scheduler_tick.has_value()) {
    debugAssertActiveSetDescriptorFresh(active_set, state, *expected_scheduler_tick);
  } else {
    if (active_set.particles_from_scheduler || active_set.cells_from_scheduler) {
      throw std::invalid_argument(
          "StepOrchestrator::executeSingleStep: scheduler-derived active sets require explicit scheduler tick; "
          "use executeSchedulerSubstep or pass expected_scheduler_tick");
    }
    debugAssertActiveSetDescriptorFresh(active_set, state);
  }

  COSMOSIM_PROFILE_SCOPE(profiler_session, "step_orchestrator.execute_single_step");

  if (profiler_session != nullptr) {
    profiler_session->counters().setCount(
        "active_particles",
        static_cast<std::uint64_t>(active_set.particle_indices.size()));
    profiler_session->counters().setCount(
        "active_cells",
        static_cast<std::uint64_t>(active_set.cell_indices.size()));
    profiler_session->counters().addCount("step_invocations", 1);
  }

  StepBoundaryState boundary = classifyStepBoundary(
      state,
      active_set,
      expected_scheduler_tick.has_value(),
      requested_boundary_kind);
  if (boundary.local_substep && integrator_state.step_index == 0U &&
      integrator_state.pm_refresh_enabled && !integrator_state.pm_long_range_field_valid) {
    boundary.pm_refresh_allowed = true;
  }
  CosmologicalTimeline timeline(cosmology_background, integrator_state.time_si_per_code);
  const CosmologicalStepFactors timeline_step = timeline.prepareStep(
      integrator_state.current_time_code,
      integrator_state.current_scale_factor,
      integrator_state.dt_time_code);
  integrator_state.current_boundary_kind = boundary.kind;
  integrator_state.inside_kdk_step = true;

  try {
  StepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .active_set = active_set,
      .active_gravity_particles = {},
      .has_active_gravity_particles = false,
      .workspace = workspace,
      .cosmology_background = cosmology_background,
      .mode_policy = mode_policy,
      .profiler_session = profiler_session,
      .particle_scheduler = nullptr,
      .gas_cell_scheduler = nullptr,
      .newly_created_particle_ids = nullptr,
      .timeline_step = timeline_step,
      .boundary = boundary,
      .pm_refresh_directive = {},
      .stage = IntegrationStage::kGravityKickPre,
  };

  const auto ordered_stages = m_scheduler.schedule(integrator_state, active_set);
  if (!isCanonicalIntegrationStageOrder(ordered_stages)) {
    throw std::runtime_error("stage scheduler order deviates from canonical kick-drift-kick contract");
  }

  std::uint32_t source_dispatch_index = 0U;
  for (const auto stage : ordered_stages) {
    if (stage == IntegrationStage::kOutputCheck) {
      continue;
    }
    context.stage = stage;
    context.source_split_phase = SourceSplitPhase::kUnsplitCompatibility;
    if (stage == IntegrationStage::kSourceTerms) {
      context.source_split_phase = source_dispatch_index++ == 0U
          ? SourceSplitPhase::kStrangPreHalf
          : SourceSplitPhase::kStrangPostHalf;
    }
    context.boundary = boundary;
    context.pm_refresh_directive = {};
    context.has_active_gravity_particles = false;
    context.active_gravity_particles = {};
    if (stage == IntegrationStage::kDrift) {
      if (workspace == nullptr) {
        throw std::runtime_error(
            "drift stage requires TransientStepWorkspace for the active gravity view");
      }
      if (active_set.hasParticleSubset(state.particles.size())) {
        context.active_gravity_particles = buildGravityParticleKernelView(
            state,
            active_set.particle_indices,
            *workspace);
      } else {
        context.active_gravity_particles = buildGravityParticleKernelViewAllParticlesDirect(state, *workspace);
      }
      context.has_active_gravity_particles = true;
    }
    if (stage == IntegrationStage::kGravityKickPost &&
        integrator_state.pm_refresh_enabled &&
        integrator_state.pm_source_generation != state.gravitySourceGeneration()) {
      // Hydro/source work may mutate authoritative gas mass, source membership,
      // or ownership after the scheduled post-drift PM refresh. Record local
      // evidence here, but do not advance cadence state: the workflow must first
      // globalize the request so every rank enters the same TreePM path.
      context.pm_refresh_directive.force_refresh_surface = true;
      context.pm_refresh_directive.reason = PmRefreshDirective::Reason::kSourceMutationForceRefresh;
      context.pm_refresh_directive.force_evaluation_scale_factor = timeline_step.scale_factor_end;
      context.pm_refresh_directive.cadence_opportunity_allowed = boundary.pm_refresh_allowed;
      context.pm_refresh_directive.sync_event_requested = true;
    }
    if (stage == IntegrationStage::kGravityKickPre) {
      // The initial-condition surface is globally synchronized even when the
      // first hierarchical scheduler tick activates only a subset.  This is
      // the sole local-substep exception: subsequent PM refreshes remain
      // restricted to consensus global boundaries.
      context.pm_refresh_directive.initial_cache_bootstrap_allowed =
          !boundary.local_substep ||
          (integrator_state.step_index == 0U && !integrator_state.pm_long_range_field_valid);
      if (integrator_state.pm_refresh_enabled && context.pm_refresh_directive.initial_cache_bootstrap_allowed &&
          !integrator_state.pm_long_range_field_valid) {
        context.pm_refresh_directive.reason = PmRefreshDirective::Reason::kInitialForceBootstrap;
        context.pm_refresh_directive.force_evaluation_scale_factor = timeline_step.scale_factor_begin;
        context.pm_refresh_directive.cadence_opportunity_allowed = true;
        context.pm_refresh_directive.sync_event_requested = true;
      }
    } else if (stage == IntegrationStage::kForceRefresh) {
      context.pm_refresh_directive.force_refresh_surface = true;
      context.pm_refresh_directive.requires_predicted_inactive_sources = boundary.local_substep;
      context.pm_refresh_directive.reason = PmRefreshDirective::Reason::kScheduledForceRefreshStage;
      context.pm_refresh_directive.force_evaluation_scale_factor = timeline_step.scale_factor_end;
      // A scheduler-local active set is still a collective force-evaluation
      // surface when the workflow supplied global synchronization metadata:
      // every rank enters this stage, including ranks with zero local targets.
      // Inactive sources are predicted to the evaluation epoch above, so
      // cadence one refreshes the PM field here instead of silently freezing it
      // until a fully active particle tick.
      const bool legal_pm_refresh_boundary = boundary.pm_refresh_allowed;
      context.pm_refresh_directive.cadence_opportunity_allowed = legal_pm_refresh_boundary;
      if (integrator_state.pm_refresh_enabled && legal_pm_refresh_boundary) {
        context.pm_refresh_directive.sync_event_requested = true;
      }
    }
    const std::size_t particle_count_before_stage = state.particles.size();
    dispatcher(context, false);
    // Position drift is an authoritative gravity-source mutation. Hydro and
    // source stages only force a new source epoch when gas exists; DMO source
    // callbacks therefore do not manufacture a second PM solve. Lower-level
    // membership/refinement/migration mutations also advance the same epoch
    // through SimulationState generation APIs.
    if (stage == IntegrationStage::kDrift) {
      // Drift is a globally ordered source-position epoch even on locally
      // empty ranks. Keeping the generation sequence rank-consistent lets
      // TreePM reject stale distributed state without synthesizing O(N) hashes.
      state.bumpGravitySourceGeneration();
    }
    if ((stage == IntegrationStage::kForceRefresh ||
         stage == IntegrationStage::kGravityKickPost ||
         (stage == IntegrationStage::kGravityKickPre && context.pm_refresh_directive.has_sync_event)) &&
        context.pm_refresh_directive.has_sync_event) {
      if (!context.pm_refresh_directive.solver_executed) {
        throw std::runtime_error("integrator-issued PM refresh directive was not consumed by a TreePM callback");
      }
      const PmSyncEvent event{
          .gravity_kick_opportunity = context.pm_refresh_directive.gravity_kick_opportunity,
          .refresh_long_range_field = context.pm_refresh_directive.refresh_long_range_field,
          .field_version = context.pm_refresh_directive.field_version,
          .last_refresh_opportunity = context.pm_refresh_directive.last_refresh_opportunity,
          .field_built_step_index = context.pm_refresh_directive.field_built_step_index,
          .field_built_scale_factor = context.pm_refresh_directive.field_built_scale_factor,
      };
      // Cadence truth advances only after the solver consumed the globally
      // authorized directive. A failing callback therefore cannot leave one
      // rank with an accepted opportunity while peers remain on the old epoch.
      if (!context.pm_refresh_directive.refresh_long_range_field &&
          !integrator_state.pm_long_range_field_valid) {
        throw std::runtime_error(
            "PM sync reused a long-range field before the integrator marked one valid");
      }
      integrator_state.pm_sync_state.commitKickOpportunity(event);
      if (context.pm_refresh_directive.refresh_long_range_field) {
        integrator_state.pm_sync_state.commitRefresh(event);
        integrator_state.pm_long_range_field_valid = true;
        integrator_state.pm_source_generation = state.gravitySourceGeneration();
      }
    }
    if (stage == IntegrationStage::kDrift && state.hasHomogeneousDmoMetadata()) {
      if (active_set.particle_indices.size() != state.particles.size()) {
        throw std::runtime_error("compact homogeneous DMO drift requires the all-active rung-zero population");
      }
      state.updateAllParticleDriftEpoch(timeline_step.time_end_code, timeline_step.scale_factor_end);
    } else if (stage == IntegrationStage::kDrift &&
        state.particle_sidecar.last_drift_time_code.size() == state.particles.size() &&
        state.particle_sidecar.last_drift_scale_factor.size() == state.particles.size()) {
      for (const std::uint32_t particle_index : active_set.particle_indices) {
        if (particle_index >= state.particle_sidecar.last_drift_time_code.size()) {
          throw std::out_of_range("drift epoch update particle index out of range");
        }
        state.particle_sidecar.last_drift_time_code[particle_index] = timeline_step.time_end_code;
        state.particle_sidecar.last_drift_scale_factor[particle_index] = timeline_step.scale_factor_end;
      }
    } else if (stage == IntegrationStage::kDrift && !state.particle_sidecar.last_drift_time_code.empty()) {
      throw std::runtime_error("particle drift-time sidecar is partially sized before drift epoch update");
    }
    if (stage == IntegrationStage::kSourceTerms && state.particles.size() > particle_count_before_stage) {
      if (state.hasHomogeneousDmoMetadata()) {
        throw std::runtime_error("particle creation committed while metadata remained homogeneous DMO; materialize before species mutation");
      }
      if (state.particle_sidecar.last_drift_time_code.size() != state.particles.size() ||
          state.particle_sidecar.last_drift_scale_factor.size() != state.particles.size()) {
        throw std::runtime_error("particle drift-time sidecar is not sized after source-term mutation");
      }
      for (std::size_t particle_index = particle_count_before_stage; particle_index < state.particles.size(); ++particle_index) {
        state.particle_sidecar.last_drift_time_code[particle_index] = timeline_step.time_end_code;
        state.particle_sidecar.last_drift_scale_factor[particle_index] = timeline_step.scale_factor_end;
      }
    }
  }

  integrator_state.inside_kdk_step = false;
  timeline.commitStep(integrator_state, timeline_step);
  integrator_state.last_completed_boundary_kind = boundary.kind;
  integrator_state.last_completed_restart_safe = boundary.restart_safe;
  } catch (...) {
    // Solver callbacks may have already mutated physical state, so pretending to
    // roll back the complete KDK step would be unsafe. Mark bookkeeping as
    // incomplete/restart-unsafe and rethrow; callers must abort or rebuild from
    // a previously committed checkpoint rather than continue from partial truth.
    integrator_state.inside_kdk_step = false;
    integrator_state.current_boundary_kind = boundary.kind;
    integrator_state.last_completed_restart_safe = false;
    throw;
  }
}

void StepOrchestrator::executeHierarchicalBlockWithDispatcher(
    SimulationState& state, IntegratorState& integrator,
    HierarchicalTimeBinScheduler& particles, HierarchicalTimeBinScheduler& cells,
    const StageDispatchFunction& dispatcher, const LambdaCdmBackground* background,
    TransientStepWorkspace& workspace, const ModePolicy* mode_policy,
    ProfilerSession* profiler,
    const std::function<void(std::exception_ptr, std::string_view)>& preparation_failure_sync,
    const HierarchicalPmCyclePlan& pm_cycle) const {
  // Core owns the operator; the workflow supplies its existing collective
  // failure authority. No rank may enter the next solver collective after a
  // peer failed local scheduler/timeline preparation.
  const auto prepare = [&](std::string_view phase, auto&& operation) {
    std::exception_ptr failure;
    try { operation(); } catch (...) { failure = std::current_exception(); }
    if (preparation_failure_sync) preparation_failure_sync(failure, phase);
    else if (failure) std::rethrow_exception(failure);
  };
  std::uint64_t periods = 0U;
  std::uint64_t start_tick = 0U;
  std::uint64_t source_begin = 0U;
  CosmologicalStepFactors coarse{};
  const CosmologicalTimeline timeline(background, integrator.time_si_per_code);
  std::array<double, k_hierarchical_timeline_capacity> scale{};
  const auto all_rows = std::span<const std::uint32_t>(workspace.gravity_particle_index_scratch);
  std::span<const std::uint32_t> active;
  prepare("hierarchical block timeline and scheduler preparation", [&] {
    if (particles.maxBin() == 0U || particles.maxBin() > 12U ||
        cells.maxBin() != particles.maxBin() || cells.elementCount() != 0U ||
        state.cells.size() != 0U || !state.hasHomogeneousDmoMetadata() ||
        particles.elementCount() != state.particles.size() ||
        particles.currentTick() != cells.currentTick() || integrator.inside_kdk_step) {
      throw std::logic_error("hierarchical KDK requires synchronized homogeneous DMO scheduler/state");
    }
    periods = particles.binPeriodTicks(particles.maxBin());
    start_tick = particles.currentTick();
    if (start_tick % periods != 0U || start_tick > std::numeric_limits<std::uint64_t>::max() - periods ||
        state.gravitySourceGeneration() >= std::numeric_limits<std::uint64_t>::max() - periods) {
      throw std::overflow_error("hierarchical KDK block is unaligned or its epoch would overflow");
    }
    coarse = timeline.prepareStep(integrator.current_time_code,
        integrator.current_scale_factor, integrator.dt_time_code * static_cast<double>(periods));
    source_begin = state.gravitySourceGeneration();
    // Bounded by maxBin<=12, independent of N. The workflow reserves this stack
    // live set through the process MemoryGovernor before entering this method.
    scale[0] = coarse.scale_factor_begin;
    for (std::size_t half_tick = 1U; half_tick < 2U * periods; ++half_tick) {
      scale[half_tick] = background != nullptr
          ? advanceScaleFactorByCosmicTime(*background, coarse.scale_factor_begin,
              0.5 * static_cast<double>(half_tick) * integrator.dt_time_code * integrator.time_si_per_code)
          : coarse.scale_factor_begin;
    }
    scale[2U * periods] = coarse.scale_factor_end;
    for (std::size_t i = 1U; i <= 2U * periods; ++i) {
      if (!std::isfinite(scale[i]) || scale[i] <= 0.0 || scale[i] < scale[i - 1U]) {
        throw std::runtime_error("hierarchical cosmological timeline is not finite positive monotonic");
      }
    }
    if (all_rows.size() != state.particles.size()) throw std::logic_error("hierarchical drift identity scratch is not admitted");
    for (std::size_t row = 0U; row < all_rows.size(); ++row) {
      if (all_rows[row] != row) throw std::logic_error("hierarchical all-source drift requires the exact identity row set");
    }
    if (!std::isfinite(integrator.dt_time_code) || integrator.dt_time_code <= 0.0 ||
        state.homogeneousDmoLastDriftTimeCode() != integrator.current_time_code ||
        state.homogeneousDmoLastDriftScaleFactor() != integrator.current_scale_factor) {
      throw std::logic_error("hierarchical block requires positive quantum and synchronized source position epoch");
    }
    active = particles.beginSubstep();
    static_cast<void>(cells.beginSubstep());
    if (active.size() != all_rows.size()) throw std::logic_error("block start is not all-active");
    integrator.inside_kdk_step = true;
    integrator.last_completed_restart_safe = false;
  });
  const auto kick_integral = [&](std::size_t begin, std::size_t end, double a_now) {
    // dp/dt=A/a. Store u=p/a at the COMMON current position epoch, so kick
    // increments are integral(dt/a)/a_now, without a separate Hubble drag.
    return background != nullptr
        ? computeComovingDriftFactor(*background, scale[begin], scale[end], 64U) /
            (integrator.time_si_per_code * a_now)
        : 0.5 * static_cast<double>(end - begin) * integrator.dt_time_code;
  };
  // PM surfaces of this block. Reference (plan disabled): both coarse
  // endpoints. Multirate: the opening surface only when the block opens a PM
  // cycle, the closing surface (refresh + kick) only when it closes one.
  const bool pm_opens = !pm_cycle.enabled || pm_cycle.opens_cycle;
  const bool pm_closes = !pm_cycle.enabled || pm_cycle.closes_cycle;
  if (pm_cycle.enabled) {
    prepare("hierarchical multirate PM cycle plan validation", [&] {
      if (pm_cycle.opens_cycle && !(std::isfinite(pm_cycle.opening_pm_kick_factor_code) &&
                                    pm_cycle.opening_pm_kick_factor_code > 0.0)) {
        throw std::logic_error("multirate PM cycle opening kick factor must be finite positive");
      }
      if (pm_cycle.closes_cycle &&
          (background != nullptr
               ? !(std::isfinite(pm_cycle.cycle_mid_scale_factor) && pm_cycle.cycle_mid_scale_factor > 0.0)
               : !(std::isfinite(pm_cycle.closing_half_interval_time_code) &&
                   pm_cycle.closing_half_interval_time_code > 0.0))) {
        throw std::logic_error("multirate PM cycle closing kick parameters must be finite positive");
      }
    });
  }
  const auto closing_pm_kick = [&](double a_end) {
    if (!pm_cycle.enabled) return kick_integral(periods, 2U * periods, a_end);
    return background != nullptr
        ? computeComovingDriftFactor(*background, pm_cycle.cycle_mid_scale_factor, a_end, 64U) /
            (integrator.time_si_per_code * a_end)
        : pm_cycle.closing_half_interval_time_code;
  };
  try {
    for (std::uint64_t offset = 0U; offset < periods; ++offset) {
      StepContext context{.state = state, .integrator_state = integrator};
      prepare("hierarchical opening kick preparation", [&] {
        const double a_begin = scale[2U * offset];
        const double a_end = scale[2U * (offset + 1U)];
        CosmologicalStepFactors fine = coarse;
        fine.time_begin_code = coarse.time_begin_code + static_cast<double>(offset) * integrator.dt_time_code;
        fine.time_end_code = offset + 1U == periods ? coarse.time_end_code
            : coarse.time_begin_code + static_cast<double>(offset + 1U) * integrator.dt_time_code;
        fine.dt_time_code = fine.time_end_code - fine.time_begin_code;
        if (!(fine.dt_time_code > 0.0) || !std::isfinite(fine.dt_time_code)) {
          throw std::runtime_error("hierarchical fine interval is not finite positive representable");
        }
        fine.dt_time_si = fine.dt_time_code * integrator.time_si_per_code;
        fine.scale_factor_begin = a_begin;
        fine.scale_factor_midpoint = scale[2U * offset + 1U];
        fine.scale_factor_end = a_end;
        fine.redshift_begin = 1.0 / a_begin - 1.0;
        fine.redshift_end = 1.0 / a_end - 1.0;
        fine.hubble_begin_code = background != nullptr ? background->hubbleSi(a_begin) * integrator.time_si_per_code : 0.0;
        fine.hubble_end_code = background != nullptr ? background->hubbleSi(a_end) * integrator.time_si_per_code : 0.0;
        fine.hubble_drag_factor = a_begin / a_end;
        // D: x += p integral(dt/a^2). Since u=p/a_begin, multiply that
        // integral by a_begin. This retains exact free Hubble drag in D.
        if (background != nullptr) {
          constexpr std::uint32_t k_samples = 64U;
          const double da = (a_end - a_begin) / k_samples;
          double drift = 0.0;
          for (std::uint32_t sample = 0U; sample < k_samples; ++sample) {
            const double a = a_begin + (static_cast<double>(sample) + 0.5) * da;
            drift += 1.0 / (a * a * a * background->hubbleSi(a));
          }
          fine.drift_factor_code = a_begin * drift * da / integrator.time_si_per_code;
        } else fine.drift_factor_code = fine.dt_time_code;
        context.active_set = makeSchedulerActiveSetDescriptor(particles, state, active, {});
        context.active_gravity_particles = {};
        context.has_active_gravity_particles = false;
        context.workspace = &workspace;
        context.particle_scheduler = &particles;
        context.gas_cell_scheduler = &cells;
        context.cosmology_background = background;
        context.mode_policy = mode_policy;
        context.profiler_session = profiler;
        context.timeline_step = fine;

        context.hierarchical_kdk.enabled = true;
        context.hierarchical_kdk.coarse_source_generation = source_begin;
        context.hierarchical_kdk.evaluation_tick = particles.currentTick();
        context.active_set.has_global_synchronization_metadata = true;
        context.active_set.globally_complete_active_set = offset == 0U;
        context.boundary = {.kind = StepBoundaryKind::kLocalActiveBinStep,
            .restart_safe = false, .output_safe = false, .pm_refresh_allowed = false, .local_substep = true};
        integrator.current_boundary_kind = context.boundary.kind;
        for (std::uint8_t bin = 0U; bin <= particles.maxBin(); ++bin) {
          const auto period = particles.binPeriodTicks(bin);
          if (offset % period == 0U) context.hierarchical_kdk.tree_kick_factor_code[bin] =
              kick_integral(2U * offset, 2U * offset + period, a_begin);
        }
        if (offset == 0U && pm_opens) {
          context.hierarchical_kdk.pm_kick_factor_code = pm_cycle.enabled
              ? pm_cycle.opening_pm_kick_factor_code : kick_integral(0U, periods, a_begin);
        }
        context.stage = IntegrationStage::kGravityKickPre;
      });
      dispatcher(context, false);
      prepare("hierarchical all-source drift preparation", [&] {
        particles.endSubstep();
        cells.endSubstep();
        context.stage = IntegrationStage::kDrift;
        context.active_set = {.particle_indices = all_rows, .cell_indices = {}};
        context.active_gravity_particles = buildGravityParticleKernelViewAllParticlesDirect(state, workspace);
        context.has_active_gravity_particles = true;
      });
      dispatcher(context, false);
      prepare("hierarchical closing force preparation", [&] {
        state.bumpGravitySourceGeneration();
        state.updateAllParticleDriftEpoch(context.timeline_step.time_end_code, context.timeline_step.scale_factor_end);
        if (state.gravitySourceGeneration() != source_begin + offset + 1U) {
          throw std::logic_error("hierarchical source epoch changed outside the authorized all-source drift");
        }
        active = particles.beginSubstep();
        static_cast<void>(cells.beginSubstep());
        context.active_set = makeSchedulerActiveSetDescriptor(particles, state, active, {});
        context.active_set.has_global_synchronization_metadata = true;
        context.active_set.globally_complete_active_set = offset + 1U == periods;
        context.has_active_gravity_particles = false;
        context.hierarchical_kdk.evaluation_tick = particles.currentTick();
        context.hierarchical_kdk.synchronization_end = offset + 1U == periods;
        const bool closing_pm_surface = offset + 1U == periods && pm_closes;
        context.hierarchical_kdk.include_long_range_force = closing_pm_surface;
        context.hierarchical_kdk.tree_kick_factor_code.fill(0.0);
        for (std::uint8_t bin = 0U; bin <= particles.maxBin(); ++bin) {
          const auto period = particles.binPeriodTicks(bin);
          if ((offset + 1U) % period == 0U) context.hierarchical_kdk.tree_kick_factor_code[bin] =
              kick_integral(2U * (offset + 1U) - period, 2U * (offset + 1U), context.timeline_step.scale_factor_end);
        }
        context.hierarchical_kdk.pm_kick_factor_code = closing_pm_surface
            ? closing_pm_kick(context.timeline_step.scale_factor_end) : 0.0;
        context.boundary.pm_refresh_allowed = closing_pm_surface;
        context.pm_refresh_directive = {.force_refresh_surface = true,
            .cadence_opportunity_allowed = closing_pm_surface,
            .sync_event_requested = closing_pm_surface,
            .force_evaluation_scale_factor = context.timeline_step.scale_factor_end,
            .reason = PmRefreshDirective::Reason::kScheduledForceRefreshStage};
        context.stage = IntegrationStage::kForceRefresh;
      });
      dispatcher(context, false);
      prepare("hierarchical PM endpoint commit", [&] {
        if (context.hierarchical_kdk.synchronization_end && pm_closes && integrator.pm_refresh_enabled &&
            !context.pm_refresh_directive.has_sync_event) {
          throw std::logic_error("coarse PM endpoint omitted synchronized refresh event");
        }
        if (context.pm_refresh_directive.has_sync_event) {
          const auto& d = context.pm_refresh_directive;
          if (!d.solver_executed || !d.refresh_long_range_field) throw std::logic_error("coarse PM endpoint did not refresh");
          if (d.field_built_step_index != integrator.step_index ||
              d.field_built_scale_factor != context.timeline_step.scale_factor_end) {
            throw std::logic_error("coarse PM endpoint carries a stale field epoch");
          }
          const PmSyncEvent event{d.gravity_kick_opportunity, d.refresh_long_range_field,
              d.field_version, d.last_refresh_opportunity, d.field_built_step_index, d.field_built_scale_factor};
          integrator.pm_sync_state.commitKickOpportunity(event);
          integrator.pm_sync_state.commitRefresh(event);
          integrator.pm_long_range_field_valid = true;
          integrator.pm_source_generation = state.gravitySourceGeneration();
        }
        context.pm_refresh_directive = {};
        context.stage = IntegrationStage::kSourceTerms;
        context.source_split_phase = SourceSplitPhase::kStrangPreHalf;
      });
      dispatcher(context, false);
      context.stage = IntegrationStage::kHydroUpdate;
      context.source_split_phase = SourceSplitPhase::kUnsplitCompatibility;
      dispatcher(context, false);
      context.stage = IntegrationStage::kSourceTerms;
      context.source_split_phase = SourceSplitPhase::kStrangPostHalf;
      dispatcher(context, false);
      context.stage = IntegrationStage::kGravityKickPost;
      context.source_split_phase = SourceSplitPhase::kUnsplitCompatibility;
      dispatcher(context, false);
      prepare("hierarchical synchronization commit", [&] {
        context.stage = IntegrationStage::kAnalysisHooks;
        if (offset + 1U == periods) {
          particles.closeSynchronizationPoint();
          cells.closeSynchronizationPoint();
          integrator.inside_kdk_step = false;
          timeline.commitStep(integrator, coarse);
          integrator.current_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint;
          integrator.last_completed_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint;
          // A multirate block that does not close its PM cycle leaves momenta
          // with an open outer PM half-kick: synchronized but not restart/output safe.
          integrator.last_completed_restart_safe = pm_closes;
          state.metadata.step_index = integrator.step_index;
          state.metadata.scale_factor = integrator.current_scale_factor;
          context.boundary = {.kind = StepBoundaryKind::kGlobalSynchronizationPoint,
              .restart_safe = pm_closes, .output_safe = pm_closes,
              .multirate_pm_cycle_interior = !pm_closes};
        }
      });
      dispatcher(context, false);
      context.stage = IntegrationStage::kOutputCheck;
      dispatcher(context, offset + 1U == periods);
    }
  } catch (...) {
    integrator.inside_kdk_step = false;
    integrator.last_completed_restart_safe = false;
    integrator.current_boundary_kind = StepBoundaryKind::kLocalActiveBinStep;
    throw;
  }
}

void TimeStepCriteriaRegistry::registerCflHook(CriteriaHook hook) { m_hooks.cfl_hook = std::move(hook); }

void TimeStepCriteriaRegistry::registerGravityHook(CriteriaHook hook) { m_hooks.gravity_hook = std::move(hook); }

void TimeStepCriteriaRegistry::registerSourceHook(CriteriaHook hook) { m_hooks.source_hook = std::move(hook); }

void TimeStepCriteriaRegistry::registerUserClampHook(CriteriaHook hook) {
  m_hooks.user_clamp_hook = std::move(hook);
}

const TimeStepCriteriaHooks& TimeStepCriteriaRegistry::hooks() const noexcept { return m_hooks; }


void StepOrchestrator::executeSchedulerSubstep(
    SimulationState& state,
    IntegratorState& integrator_state,
    const HierarchicalTimeBinScheduler& scheduler,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices,
    const LambdaCdmBackground* cosmology_background,
    TransientStepWorkspace* workspace,
    const ModePolicy* mode_policy,
    ProfilerSession* profiler_session,
    StepBoundaryKind requested_boundary_kind) const {
  const ActiveSetDescriptor active_set = makeSchedulerActiveSetDescriptor(
      scheduler,
      state,
      active_particle_indices,
      active_cell_indices);
  executeSingleStep(
      state,
      integrator_state,
      active_set,
      cosmology_background,
      workspace,
      mode_policy,
      profiler_session,
      scheduler.currentTick(),
      requested_boundary_kind);
}

HierarchicalTimeBinScheduler::HierarchicalTimeBinScheduler(std::uint8_t max_bin) : m_max_bin(max_bin) {
  m_representation = max_bin == 0U ? SchedulerRepresentation::kUniformRungZero : SchedulerRepresentation::kGenericHierarchical;
  if (max_bin > k_max_representable_bin) {
    throw std::invalid_argument(
        "HierarchicalTimeBinScheduler: max_bin exceeds the 64-bit tick representation limit");
  }
}

void HierarchicalTimeBinScheduler::reset(
    std::uint32_t element_count,
    std::uint8_t initial_bin,
    std::uint64_t start_tick) {
  m_current_tick = start_tick;
  m_element_count = element_count;
  m_representation = m_max_bin == 0U
      ? SchedulerRepresentation::kUniformRungZero
      : SchedulerRepresentation::kGenericHierarchical;
  m_uniform_pending_bin_index = k_unset_pending_bin;
  m_last_reconciliation = {};
  m_substep_open = false;
  m_active_elements.clear();
  m_active_sort_scratch.clear();
  m_hot = {};
  m_elements_by_bin.clear();
  m_position_in_bin.clear();
  m_candidate_bin_index.clear();
  m_candidate_source.clear();
  m_diagnostics = {};
  m_diagnostics.occupancy_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  m_diagnostics.active_count_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);

  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    if (initial_bin != 0U) {
      throw std::invalid_argument("uniform rung-zero scheduler requires initial_bin == 0");
    }
    m_active_elements.resize(element_count);
    std::iota(m_active_elements.begin(), m_active_elements.end(), 0U);
    m_diagnostics.occupancy_by_bin[0] = element_count;
    // Permanently release any generic capacity when transitioning to uniform mode.
    TimeBinHotMetadata{}.bin_index.swap(m_hot.bin_index);
    TimeBinHotMetadata{}.next_activation_tick.swap(m_hot.next_activation_tick);
    TimeBinHotMetadata{}.active_flag.swap(m_hot.active_flag);
    TimeBinHotMetadata{}.pending_bin_index.swap(m_hot.pending_bin_index);
    std::vector<std::vector<std::uint32_t>>{}.swap(m_elements_by_bin);
    std::vector<std::size_t>{}.swap(m_position_in_bin);
    std::vector<std::uint32_t>{}.swap(m_active_sort_scratch);
    std::vector<std::uint8_t>{}.swap(m_candidate_bin_index);
    std::vector<TimeStepCandidateSource>{}.swap(m_candidate_source);
  } else {
    const std::uint8_t clamped_bin = clampBin(initial_bin);
    m_hot.bin_index.assign(element_count, clamped_bin);
    m_hot.next_activation_tick.assign(element_count, m_current_tick);
    m_hot.active_flag.assign(element_count, 0);
    m_hot.pending_bin_index.assign(element_count, k_unset_pending_bin);
    m_candidate_bin_index.assign(element_count, k_unset_pending_bin);
    m_candidate_source.assign(element_count, TimeStepCandidateSource::kUserClamp);
    m_position_in_bin.resize(element_count, 0);
    m_elements_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, {});
    m_elements_by_bin[clamped_bin].reserve(element_count);
    for (std::uint32_t element = 0; element < element_count; ++element) {
      m_position_in_bin[element] = m_elements_by_bin[clamped_bin].size();
      m_elements_by_bin[clamped_bin].push_back(element);
    }
    m_diagnostics.occupancy_by_bin[clamped_bin] = element_count;
  }
  refreshOwnedCapacityHighWater();
}


void HierarchicalTimeBinScheduler::planAppendCapacity(
    RetainedCapacityTransaction& plan,
    std::uint32_t new_element_count, std::uint8_t initial_bin) {
  if (new_element_count == 0U) return;
  const std::size_t target = checkedSizeAdd(m_element_count, new_element_count,
      "scheduler population capacity target");
  (void)checkedIntegralNarrow<std::uint32_t>(target,
      "scheduler population capacity index range");
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    if (initial_bin != 0U) {
      throw std::invalid_argument("uniform rung-zero append requires initial_bin == 0");
    }
    plan.add(m_active_sort_scratch, target);
    plan.add(m_elements_by_bin, 1U);
    // Population creation is the boundary where a homogeneous-DMO scheduler
    // may need a future activation tick. Budget the one-time generic
    // materialization here so appendElements() never creates an ungoverned
    // O(N) compatibility lane when that happens. Bin-0 membership is staged
    // through the already-governed sort scratch so any open active-set span
    // remains valid until the normal substep invalidation boundary.
    plan.add(m_hot.bin_index, target);
    plan.add(m_hot.next_activation_tick, target);
    plan.add(m_hot.active_flag, target);
    plan.add(m_hot.pending_bin_index, target);
    plan.add(m_position_in_bin, target);
    plan.add(m_candidate_bin_index, target);
    plan.add(m_candidate_source, target);
    return;
  }
  plan.add(m_hot.bin_index, target);
  plan.add(m_hot.next_activation_tick, target);
  plan.add(m_hot.active_flag, target);
  plan.add(m_hot.pending_bin_index, target);
  plan.add(m_position_in_bin, target);
  plan.add(m_candidate_bin_index, target);
  plan.add(m_candidate_source, target);
  if (m_elements_by_bin.empty()) {
    throw std::logic_error("scheduler population append requires initialized bins");
  }
  auto& members = m_elements_by_bin[clampBin(initial_bin)];
  plan.add(members, checkedSizeAdd(members.size(), new_element_count,
      "scheduler bin population capacity target"));
}


void HierarchicalTimeBinScheduler::materializeGenericRepresentationForAppend() {
  if (m_representation != SchedulerRepresentation::kUniformRungZero) {
    return;
  }
  const std::size_t count = static_cast<std::size_t>(m_element_count);
  m_hot.bin_index.assign(count, 0U);
  m_hot.next_activation_tick.assign(count, m_current_tick);
  m_hot.active_flag.assign(count, m_substep_open ? 1U : 0U);
  m_hot.pending_bin_index.assign(count, m_uniform_pending_bin_index);
  m_position_in_bin.resize(count);
  std::iota(m_position_in_bin.begin(), m_position_in_bin.end(), std::size_t{0});
  m_candidate_bin_index.assign(count, k_unset_pending_bin);
  m_candidate_source.assign(count, TimeStepCandidateSource::kUserClamp);

  m_active_sort_scratch.resize(count);
  std::iota(m_active_sort_scratch.begin(), m_active_sort_scratch.end(), 0U);
  m_elements_by_bin.clear();
  m_elements_by_bin.resize(1U);
  m_elements_by_bin[0] = std::move(m_active_sort_scratch);
  // Keep m_active_elements untouched while a substep is open: StepContext may
  // still hold a span into the current active set. The normal endSubstep()/
  // next beginSubstep rebuild is the safe invalidation boundary.
  m_uniform_pending_bin_index = k_unset_pending_bin;
  m_representation = SchedulerRepresentation::kGenericHierarchical;
  refreshOwnedCapacityHighWater();
  validateInternalState("HierarchicalTimeBinScheduler::materializeGenericRepresentationForAppend");
}

void HierarchicalTimeBinScheduler::appendElements(
    std::uint32_t new_element_count,
    std::uint8_t initial_bin,
    std::uint64_t first_activation_tick) {
  if (new_element_count == 0U) return;
  const std::uint32_t old_count = m_element_count;
  if (new_element_count > std::numeric_limits<std::uint32_t>::max() - old_count) {
    throw std::overflow_error("HierarchicalTimeBinScheduler::appendElements element count overflow");
  }
  const std::uint32_t total_count = old_count + new_element_count;
  if (m_representation == SchedulerRepresentation::kUniformRungZero &&
      initial_bin == 0U && first_activation_tick == m_current_tick) {
    m_active_elements.resize(total_count);
    for (std::uint32_t element = old_count; element < total_count; ++element) {
      m_active_elements[element] = element;
    }
    m_element_count = total_count;
    m_diagnostics.occupancy_by_bin[0] = total_count;
    refreshOwnedCapacityHighWater();
    return;
  }
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    if (initial_bin != 0U) {
      throw std::invalid_argument("uniform rung-zero scheduler append requires initial_bin == 0");
    }
    materializeGenericRepresentationForAppend();
  }
  const std::uint8_t clamped_bin = clampBin(initial_bin);
  m_hot.bin_index.resize(total_count, clamped_bin);
  m_hot.next_activation_tick.resize(total_count, first_activation_tick);
  m_hot.active_flag.resize(total_count, 0);
  m_hot.pending_bin_index.resize(total_count, k_unset_pending_bin);
  m_position_in_bin.resize(total_count, 0);
  m_candidate_bin_index.resize(total_count, k_unset_pending_bin);
  m_candidate_source.resize(total_count, TimeStepCandidateSource::kUserClamp);
  for (std::uint32_t element = old_count; element < total_count; ++element) {
    auto& members = m_elements_by_bin[clamped_bin];
    m_position_in_bin[element] = members.size();
    members.push_back(element);
  }
  m_element_count = total_count;
  m_diagnostics.occupancy_by_bin[clamped_bin] += new_element_count;
  refreshOwnedCapacityHighWater();
  validateInternalState("HierarchicalTimeBinScheduler::appendElements");
}


void HierarchicalTimeBinScheduler::setElementBin(
    std::uint32_t element_index,
    std::uint8_t bin_index,
    std::uint64_t current_tick) {
  if (element_index >= m_element_count) throw std::out_of_range("element_index out of range");
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    if (bin_index != 0U || current_tick != m_current_tick) {
      throw std::invalid_argument("uniform rung-zero scheduler cannot commit a nonzero/asynchronous bin");
    }
    return;
  }
  const std::uint8_t clamped_bin = clampBin(bin_index);
  const std::uint8_t old_bin = m_hot.bin_index[element_index];
  if (old_bin == clamped_bin) {
    const std::uint64_t period_ticks = binPeriodTicks(clamped_bin);
    m_hot.next_activation_tick[element_index] =
        (current_tick % period_ticks == 0) ? current_tick : ((current_tick / period_ticks) * period_ticks + period_ticks);
    return;
  }
  eraseFromBin(element_index, old_bin);
  insertIntoBin(element_index, clamped_bin);
  m_hot.bin_index[element_index] = clamped_bin;
  const std::uint64_t period_ticks = binPeriodTicks(clamped_bin);
  m_hot.next_activation_tick[element_index] =
      (current_tick % period_ticks == 0) ? current_tick : ((current_tick / period_ticks) * period_ticks + period_ticks);
  refreshOwnedCapacityHighWater();
}


void HierarchicalTimeBinScheduler::submitCandidateTimeStep(
    std::uint32_t element_index,
    double dt_time_code,
    const TimeStepLimits& limits,
    TimeStepCandidateSource source,
    std::string_view label) {
  const auto mapped = mapDtToTimeBin(dt_time_code, limits);
  submitCandidateBin(element_index, mapped.bin_index, source, label);
  if (mapped.clipped_to_min) {
    ++m_last_reconciliation.clipped_to_min_dt;
  }
  if (mapped.clipped_to_max) {
    ++m_last_reconciliation.clipped_to_max_dt;
  }
}

void HierarchicalTimeBinScheduler::submitCandidateBin(
    std::uint32_t element_index,
    std::uint8_t target_bin,
    TimeStepCandidateSource source,
    std::string_view label) {
  if (element_index >= m_element_count) {
    throw std::out_of_range("element_index out of range");
  }
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    (void)label;
    ++m_last_reconciliation.submitted_candidates;
    ++m_last_reconciliation.elements_with_candidates;
    const auto source_index = timestepSourceIndex(source);
    if (source_index < m_last_reconciliation.limiting_candidates_by_source.size()) {
      ++m_last_reconciliation.limiting_candidates_by_source[source_index];
      m_last_reconciliation.dominant_limiting_source = source;
    }
    return;
  }
  if (m_candidate_bin_index.size() != m_hot.size()) {
    m_candidate_bin_index.assign(m_hot.size(), k_unset_pending_bin);
    m_candidate_source.assign(m_hot.size(), TimeStepCandidateSource::kUserClamp);
  }
  // The compatibility label is intentionally transient and is not retained or
  // copied into scheduler state. Typed candidate provenance drives production
  // arbitration and later human-readable diagnostics.
  (void)label;
  const std::uint8_t clamped = clampBin(target_bin);
  ++m_last_reconciliation.submitted_candidates;
  if (m_candidate_bin_index[element_index] == k_unset_pending_bin || clamped < m_candidate_bin_index[element_index]) {
    m_candidate_bin_index[element_index] = clamped;
    m_candidate_source[element_index] = source;
  }
  refreshOwnedCapacityHighWater();
}

TimeStepReconciliationResult HierarchicalTimeBinScheduler::reconcileCandidateTransitions() {
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    TimeStepReconciliationResult result = m_last_reconciliation;
    // Rung zero has no representable transition. Candidate diagnostics remain
    // aggregate-only; no transition request is committed or stored per element.
    result.committed_transition_requests = 0U;
    m_last_reconciliation = {};
    return result;
  }
  if (m_candidate_bin_index.size() != m_hot.size()) {
    m_candidate_bin_index.assign(m_hot.size(), k_unset_pending_bin);
    m_candidate_source.assign(m_hot.size(), TimeStepCandidateSource::kUserClamp);
  }
  TimeStepReconciliationResult result = m_last_reconciliation;
  for (std::uint32_t element = 0; element < m_candidate_bin_index.size(); ++element) {
    const std::uint8_t candidate = m_candidate_bin_index[element];
    if (candidate == k_unset_pending_bin) {
      continue;
    }
    ++result.elements_with_candidates;
    const std::size_t source_index = timestepSourceIndex(m_candidate_source[element]);
    if (source_index < result.limiting_candidates_by_source.size()) {
      ++result.limiting_candidates_by_source[source_index];
      if (result.limiting_candidates_by_source[source_index] >
          result.limiting_candidates_by_source[timestepSourceIndex(result.dominant_limiting_source)]) {
        result.dominant_limiting_source = m_candidate_source[element];
      }
    }
    validateTransitionRequest(
        element,
        candidate,
        timeStepCandidateSourceName(m_candidate_source[element]));
    requestBinTransition(element, candidate);
    ++result.committed_transition_requests;
    m_candidate_bin_index[element] = k_unset_pending_bin;
    m_candidate_source[element] = TimeStepCandidateSource::kUserClamp;
  }
  m_last_reconciliation = {};
  return result;
}

void HierarchicalTimeBinScheduler::requestBinTransition(
    std::uint32_t element_index,
    std::uint8_t target_bin) {
  if (element_index >= m_element_count) {
    throw std::out_of_range("element_index out of range");
  }
  validateTransitionRequest(element_index, target_bin, "HierarchicalTimeBinScheduler::requestBinTransition");
  if (m_representation == SchedulerRepresentation::kUniformRungZero) return;
  m_hot.pending_bin_index[element_index] = clampBin(target_bin);
}

std::span<const std::uint32_t> HierarchicalTimeBinScheduler::activeElements() const noexcept {
  return m_active_elements;
}

std::uint64_t HierarchicalTimeBinScheduler::currentTick() const noexcept { return m_current_tick; }

std::uint8_t HierarchicalTimeBinScheduler::maxBin() const noexcept { return m_max_bin; }

std::uint32_t HierarchicalTimeBinScheduler::elementCount() const noexcept { return m_element_count; }

std::uint8_t HierarchicalTimeBinScheduler::binIndex(std::uint32_t element_index) const {
  if (element_index >= m_element_count) throw std::out_of_range("scheduler element index out of range");
  return m_representation == SchedulerRepresentation::kUniformRungZero ? 0U : m_hot.bin_index[element_index];
}

std::uint64_t HierarchicalTimeBinScheduler::nextActivationTick(std::uint32_t element_index) const {
  if (element_index >= m_element_count) throw std::out_of_range("scheduler element index out of range");
  return m_representation == SchedulerRepresentation::kUniformRungZero ? m_current_tick : m_hot.next_activation_tick[element_index];
}

bool HierarchicalTimeBinScheduler::isElementActive(std::uint32_t element_index) const {
  if (element_index >= m_element_count) throw std::out_of_range("scheduler element index out of range");
  return m_representation == SchedulerRepresentation::kUniformRungZero ? m_substep_open : m_hot.active_flag[element_index] != 0U;
}

std::uint8_t HierarchicalTimeBinScheduler::pendingBinIndex(std::uint32_t element_index) const {
  if (element_index >= m_element_count) throw std::out_of_range("scheduler element index out of range");
  return m_representation == SchedulerRepresentation::kUniformRungZero ? m_uniform_pending_bin_index : m_hot.pending_bin_index[element_index];
}

bool HierarchicalTimeBinScheduler::isBinActiveAtTick(std::uint8_t bin_index, std::uint64_t tick) const {
  const std::uint8_t clamped_bin = clampBin(bin_index);
  return tick % binPeriodTicks(clamped_bin) == 0;
}

std::uint64_t HierarchicalTimeBinScheduler::binPeriodTicks(std::uint8_t bin_index) const {
  return powerOfTwo(clampBin(bin_index));
}

std::span<const std::uint32_t> HierarchicalTimeBinScheduler::beginSubstep() {
  if (m_substep_open) {
    throw std::runtime_error(schedulerContextMessage(
        "beginSubstep called while a substep is already open",
        "HierarchicalTimeBinScheduler::beginSubstep",
        m_current_tick));
  }
#ifndef NDEBUG
  validateInternalState("HierarchicalTimeBinScheduler::beginSubstep");
#endif
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    m_diagnostics.active_elements = m_element_count;
    m_diagnostics.active_fraction = m_element_count == 0U ? 0.0 : 1.0;
    m_diagnostics.active_count_by_bin[0] = m_element_count;
    m_diagnostics.occupancy_by_bin[0] = m_element_count;
  } else {
    rebuildActiveSet();
  }
  refreshOwnedCapacityHighWater();
  m_substep_open = true;
#ifndef NDEBUG
  validateInternalState("HierarchicalTimeBinScheduler::beginSubstep.active_set_created");
#endif
  return m_active_elements;
}

void HierarchicalTimeBinScheduler::endSubstep() {
  if (!m_substep_open) {
    throw std::runtime_error(schedulerContextMessage(
        "endSubstep called without an open substep",
        "HierarchicalTimeBinScheduler::endSubstep",
        m_current_tick));
  }
  if (m_current_tick == std::numeric_limits<std::uint64_t>::max()) {
    throw std::runtime_error(schedulerContextMessage(
        "global integer time cannot advance monotonically without overflowing",
        "HierarchicalTimeBinScheduler::endSubstep",
        m_current_tick));
  }
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    (void)reconcileCandidateTransitions();
    ++m_current_tick;
    m_substep_open = false;
    refreshOwnedCapacityHighWater();
    return;
  }
  for (const std::uint32_t element : m_active_elements) {
    const std::uint64_t period_ticks = binPeriodTicks(m_hot.bin_index[element]);
    if (period_ticks > std::numeric_limits<std::uint64_t>::max() - m_current_tick) {
      throw std::runtime_error(timeBinContextMessage(
          "next activation tick would overflow before substep commit",
          "HierarchicalTimeBinScheduler::endSubstep",
          m_current_tick,
          element,
          m_hot.bin_index[element],
          m_hot.next_activation_tick[element]));
    }
  }

  (void)reconcileCandidateTransitions();
  applyPendingTransitions();
  for (const std::uint32_t element : m_active_elements) {
    const std::uint8_t bin = m_hot.bin_index[element];
    const std::uint64_t period_ticks = binPeriodTicks(bin);
    m_hot.next_activation_tick[element] = m_current_tick + period_ticks;
    m_hot.active_flag[element] = 0;
  }
  ++m_current_tick;
  m_substep_open = false;
  refreshOwnedCapacityHighWater();
#ifndef NDEBUG
  validateInternalState("HierarchicalTimeBinScheduler::endSubstep.committed");
#endif
}

void HierarchicalTimeBinScheduler::planSynchronizedCandidateCapacity(RetainedCapacityTransaction& plan) {
  if (m_substep_open || m_current_tick % binPeriodTicks(m_max_bin) != 0U || m_max_bin > 12U) {
    throw std::logic_error("block bin assignment requires closed aligned scheduler authority");
  }
  std::array<std::size_t, 13U> counts{};
  for (std::uint32_t row = 0U; row < m_element_count; ++row) {
    const auto bin = m_candidate_bin_index[row];
    if (bin == k_unset_pending_bin || bin > m_max_bin ||
        m_hot.next_activation_tick[row] != m_current_tick || m_hot.active_flag[row] != 0U ||
        m_hot.pending_bin_index[row] != k_unset_pending_bin) {
      throw std::logic_error("synchronized bin assignment has incomplete candidates or an open particle interval");
    }
    ++counts[bin];
  }
  for (std::uint8_t bin = 0U; bin <= m_max_bin; ++bin) plan.add(m_elements_by_bin[bin], counts[bin]);
  plan.add(m_active_elements, m_element_count);
  plan.add(m_active_sort_scratch, m_element_count);
}

void HierarchicalTimeBinScheduler::commitSynchronizedCandidates() {
  if (m_substep_open || m_current_tick % binPeriodTicks(m_max_bin) != 0U) {
    throw std::logic_error("synchronized bin commit is outside its legal boundary");
  }
  std::array<std::size_t, 13U> counts{};
  if (m_max_bin > 12U) throw std::logic_error("block bin commit exceeds bounded hierarchy");
  for (std::uint32_t row = 0U; row < m_element_count; ++row) {
    const auto bin = m_candidate_bin_index[row];
    if (bin == k_unset_pending_bin || bin > m_max_bin ||
        m_hot.active_flag[row] != 0U || m_hot.next_activation_tick[row] != m_current_tick ||
        m_hot.pending_bin_index[row] != k_unset_pending_bin) {
      throw std::logic_error("block bin candidate is missing or its interval is open");
    }
    ++counts[bin];
  }
  if (m_active_elements.capacity() < m_element_count || m_active_sort_scratch.capacity() < m_element_count) {
    throw std::logic_error("block active-set storage was not admitted");
  }
  for (std::uint8_t bin = 0U; bin <= m_max_bin; ++bin) {
    if (counts[bin] > m_elements_by_bin[bin].capacity()) throw std::logic_error("block bin storage was not admitted");
    m_elements_by_bin[bin].clear();
  }
  for (std::uint32_t row = 0U; row < m_element_count; ++row) {
    const auto bin = m_candidate_bin_index[row];
    auto& members = m_elements_by_bin[bin];
    m_position_in_bin[row] = members.size();
    members.push_back(row);
    m_hot.bin_index[row] = bin;
    m_hot.next_activation_tick[row] = m_current_tick;
    m_hot.pending_bin_index[row] = k_unset_pending_bin;
    m_candidate_bin_index[row] = k_unset_pending_bin;
    m_candidate_source[row] = TimeStepCandidateSource::kUserClamp;
  }
  m_last_reconciliation = {};
  refreshOwnedCapacityHighWater();
  validateInternalState("HierarchicalTimeBinScheduler::commitSynchronizedCandidates");
}

void HierarchicalTimeBinScheduler::closeSynchronizationPoint() {
  if (!m_substep_open || m_current_tick % binPeriodTicks(m_max_bin) != 0U ||
      m_active_elements.size() != m_element_count) {
    throw std::logic_error("scheduler block closure requires an open aligned all-active synchronization point");
  }
  if (m_representation != SchedulerRepresentation::kUniformRungZero) {
    for (const auto row : m_active_elements) {
      if (m_hot.pending_bin_index[row] != k_unset_pending_bin ||
          m_candidate_bin_index[row] != k_unset_pending_bin) {
        throw std::logic_error("block closure cannot discard an uncommitted bin transition");
      }
      m_hot.active_flag[row] = 0U;
      m_hot.next_activation_tick[row] = m_current_tick;
    }
  }
  m_substep_open = false;
  validateInternalState("HierarchicalTimeBinScheduler::closeSynchronizationPoint");
}

const TimeBinHotMetadata& HierarchicalTimeBinScheduler::hotMetadata() const noexcept { return m_hot; }

const TimeBinDiagnostics& HierarchicalTimeBinScheduler::diagnostics() const noexcept { return m_diagnostics; }

std::uint64_t HierarchicalTimeBinScheduler::logicalSizeBytes() const {
  std::uint64_t total = 0U;
  const auto add_bytes = [&total](std::uint64_t bytes, std::string_view label) {
    if (bytes > std::numeric_limits<std::uint64_t>::max() - total) {
      throw std::overflow_error(
          "HierarchicalTimeBinScheduler logical memory accounting overflow at " +
          std::string(label));
    }
    total += bytes;
  };
  const auto add_vector = [&add_bytes](const auto& values, std::string_view label) {
    using Value = typename std::decay_t<decltype(values)>::value_type;
    const std::uint64_t size = static_cast<std::uint64_t>(values.size());
    constexpr std::uint64_t element_bytes = static_cast<std::uint64_t>(sizeof(Value));
    if (size != 0U &&
        element_bytes > std::numeric_limits<std::uint64_t>::max() / size) {
      throw std::overflow_error(
          "HierarchicalTimeBinScheduler vector logical memory accounting overflow at " +
          std::string(label));
    }
    add_bytes(size * element_bytes, label);
  };

  add_vector(m_hot.bin_index, "hot.bin_index");
  add_vector(m_hot.next_activation_tick, "hot.next_activation_tick");
  add_vector(m_hot.active_flag, "hot.active_flag");
  add_vector(m_hot.pending_bin_index, "hot.pending_bin_index");
  add_vector(m_elements_by_bin, "elements_by_bin.outer");
  for (const auto& bin : m_elements_by_bin) {
    add_vector(bin, "elements_by_bin.inner");
  }
  add_vector(m_position_in_bin, "position_in_bin");
  add_vector(m_active_elements, "active_elements");
  add_vector(m_active_sort_scratch, "active_sort_scratch");
  add_vector(m_diagnostics.occupancy_by_bin, "diagnostics.occupancy_by_bin");
  add_vector(m_diagnostics.active_count_by_bin, "diagnostics.active_count_by_bin");
  add_vector(m_candidate_bin_index, "candidate_bin_index");
  add_vector(m_candidate_source, "candidate_source");
  return total;
}

std::uint64_t HierarchicalTimeBinScheduler::ownedCapacityBytes() const {
  std::uint64_t total = 0U;
  const auto add_bytes = [&total](std::uint64_t bytes, std::string_view label) {
    if (bytes > std::numeric_limits<std::uint64_t>::max() - total) {
      throw std::overflow_error(
          "HierarchicalTimeBinScheduler memory accounting overflow at " +
          std::string(label));
    }
    total += bytes;
  };
  const auto add_vector = [&add_bytes](const auto& values, std::string_view label) {
    using Value = typename std::decay_t<decltype(values)>::value_type;
    const std::uint64_t capacity = static_cast<std::uint64_t>(values.capacity());
    constexpr std::uint64_t element_bytes = static_cast<std::uint64_t>(sizeof(Value));
    if (capacity != 0U &&
        element_bytes > std::numeric_limits<std::uint64_t>::max() / capacity) {
      throw std::overflow_error(
          "HierarchicalTimeBinScheduler vector memory accounting overflow at " +
          std::string(label));
    }
    add_bytes(capacity * element_bytes, label);
  };

  add_vector(m_hot.bin_index, "hot.bin_index");
  add_vector(m_hot.next_activation_tick, "hot.next_activation_tick");
  add_vector(m_hot.active_flag, "hot.active_flag");
  add_vector(m_hot.pending_bin_index, "hot.pending_bin_index");
  add_vector(m_elements_by_bin, "elements_by_bin.outer");
  for (const auto& bin : m_elements_by_bin) {
    add_vector(bin, "elements_by_bin.inner");
  }
  add_vector(m_position_in_bin, "position_in_bin");
  add_vector(m_active_elements, "active_elements");
  add_vector(m_active_sort_scratch, "active_sort_scratch");
  add_vector(m_diagnostics.occupancy_by_bin, "diagnostics.occupancy_by_bin");
  add_vector(m_diagnostics.active_count_by_bin, "diagnostics.active_count_by_bin");
  add_vector(m_candidate_bin_index, "candidate_bin_index");
  add_vector(m_candidate_source, "candidate_source");
  return total;
}

std::uint64_t HierarchicalTimeBinScheduler::ownedCapacityHighWaterBytes() const noexcept {
  return m_owned_capacity_high_water_bytes;
}

void HierarchicalTimeBinScheduler::refreshOwnedCapacityHighWater() {
  m_owned_capacity_high_water_bytes = std::max(
      m_owned_capacity_high_water_bytes, ownedCapacityBytes());
}

TimeBinPersistentState HierarchicalTimeBinScheduler::exportPersistentState() const {
  TimeBinPersistentState persistent_state;
  exportPersistentStateInto(persistent_state);
  return persistent_state;
}

void HierarchicalTimeBinScheduler::exportPersistentStateInto(TimeBinPersistentState& persistent_state) const {
  validateInternalState("HierarchicalTimeBinScheduler::exportPersistentState");
  persistent_state.current_tick = m_current_tick;
  persistent_state.max_bin = m_max_bin;
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    persistent_state.bin_index.assign(m_element_count, 0U);
    persistent_state.next_activation_tick.assign(m_element_count, m_current_tick);
    persistent_state.active_flag.assign(m_element_count, m_substep_open ? 1U : 0U);
    persistent_state.pending_bin_index.assign(m_element_count, m_uniform_pending_bin_index);
  } else {
    persistent_state.bin_index.assign(m_hot.bin_index.begin(), m_hot.bin_index.end());
    persistent_state.next_activation_tick.assign(m_hot.next_activation_tick.begin(), m_hot.next_activation_tick.end());
    persistent_state.active_flag.assign(m_hot.active_flag.begin(), m_hot.active_flag.end());
    persistent_state.pending_bin_index.assign(m_hot.pending_bin_index.begin(), m_hot.pending_bin_index.end());
  }
}

void HierarchicalTimeBinScheduler::importPersistentState(const TimeBinPersistentState& persistent_state) {
  if (persistent_state.max_bin > k_max_representable_bin) {
    throw std::invalid_argument(
        "TimeBinPersistentState max_bin exceeds the 64-bit tick representation limit");
  }
  if (persistent_state.bin_index.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::overflow_error(
        "TimeBinPersistentState element count exceeds the scheduler uint32 index representation");
  }
  if (persistent_state.bin_index.size() != persistent_state.next_activation_tick.size() ||
      persistent_state.bin_index.size() != persistent_state.active_flag.size() ||
      persistent_state.bin_index.size() != persistent_state.pending_bin_index.size()) {
    throw std::invalid_argument("TimeBinPersistentState arrays must have matching sizes");
  }

  for (std::size_t element_index = 0; element_index < persistent_state.bin_index.size(); ++element_index) {
    if (persistent_state.bin_index[element_index] > persistent_state.max_bin) {
      throw std::invalid_argument("TimeBinPersistentState bin_index exceeds max_bin");
    }
    if (persistent_state.active_flag[element_index] > 1U) {
      throw std::invalid_argument("TimeBinPersistentState active_flag must be 0 or 1");
    }
    const bool pending_is_unset = persistent_state.pending_bin_index[element_index] == k_unset_pending_bin;
    if (!pending_is_unset && persistent_state.pending_bin_index[element_index] > persistent_state.max_bin) {
      throw std::invalid_argument("TimeBinPersistentState pending_bin_index exceeds max_bin");
    }
  }

  const std::uint8_t uniform_pending_bin = persistent_state.pending_bin_index.empty()
      ? k_unset_pending_bin
      : persistent_state.pending_bin_index.front();
  const bool uniform_rung_zero = persistent_state.max_bin == 0U &&
      (uniform_pending_bin == k_unset_pending_bin || uniform_pending_bin == 0U) &&
      std::all_of(persistent_state.bin_index.begin(), persistent_state.bin_index.end(), [](std::uint8_t v) { return v == 0U; }) &&
      std::all_of(persistent_state.next_activation_tick.begin(), persistent_state.next_activation_tick.end(),
                  [&](std::uint64_t tick) { return tick == persistent_state.current_tick; }) &&
      std::all_of(persistent_state.active_flag.begin(), persistent_state.active_flag.end(), [](std::uint8_t v) { return v == 0U; }) &&
      std::all_of(persistent_state.pending_bin_index.begin(), persistent_state.pending_bin_index.end(),
                  [&](std::uint8_t v) { return v == uniform_pending_bin; });
  if (uniform_rung_zero) {
    m_max_bin = 0U;
    m_representation = SchedulerRepresentation::kUniformRungZero;
    m_uniform_pending_bin_index = uniform_pending_bin;
    m_current_tick = persistent_state.current_tick;
    m_element_count = checkedIntegralNarrow<std::uint32_t>(persistent_state.bin_index.size(), "uniform scheduler import element count");
    m_hot = {};
    m_elements_by_bin.clear(); std::vector<std::vector<std::uint32_t>>{}.swap(m_elements_by_bin);
    m_position_in_bin.clear(); std::vector<std::size_t>{}.swap(m_position_in_bin);
    m_active_sort_scratch.clear(); std::vector<std::uint32_t>{}.swap(m_active_sort_scratch);
    m_candidate_bin_index.clear(); std::vector<std::uint8_t>{}.swap(m_candidate_bin_index);
    m_candidate_source.clear(); std::vector<TimeStepCandidateSource>{}.swap(m_candidate_source);
    m_active_elements.resize(m_element_count);
    std::iota(m_active_elements.begin(), m_active_elements.end(), 0U);
    m_substep_open = false;
    m_last_reconciliation = {};
    m_diagnostics = {};
    m_diagnostics.occupancy_by_bin.assign(1U, m_element_count);
    m_diagnostics.active_count_by_bin.assign(1U, 0U);
    refreshOwnedCapacityHighWater();
    validateInternalState("HierarchicalTimeBinScheduler::importPersistentState.uniform");
    return;
  }

  m_current_tick = persistent_state.current_tick;
  m_max_bin = persistent_state.max_bin;
  m_representation = SchedulerRepresentation::kGenericHierarchical;
  m_uniform_pending_bin_index = k_unset_pending_bin;
  m_element_count = checkedIntegralNarrow<std::uint32_t>(persistent_state.bin_index.size(), "scheduler import element count");
  m_hot.bin_index = persistent_state.bin_index;
  m_hot.next_activation_tick = persistent_state.next_activation_tick;
  m_hot.active_flag = persistent_state.active_flag;
  m_hot.pending_bin_index = persistent_state.pending_bin_index;
  m_candidate_bin_index.assign(m_hot.bin_index.size(), k_unset_pending_bin);
  m_candidate_source.assign(m_hot.bin_index.size(), TimeStepCandidateSource::kUserClamp);
  m_last_reconciliation = {};

  m_elements_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, {});
  m_position_in_bin.assign(m_hot.bin_index.size(), 0);

  std::array<std::size_t, 256U> import_occupancy{};
  for (const auto bin : m_hot.bin_index) ++import_occupancy[clampBin(bin)];
  for (std::size_t bin = 0U; bin < m_elements_by_bin.size(); ++bin) {
    m_elements_by_bin[bin].reserve(import_occupancy[bin]);
  }
  for (std::size_t element_index = 0; element_index < m_hot.bin_index.size(); ++element_index) {
    const std::uint8_t bin = clampBin(m_hot.bin_index[element_index]);
    m_hot.bin_index[element_index] = bin;
    auto& members = m_elements_by_bin[bin];
    m_position_in_bin[element_index] = members.size();
    members.push_back(checkedIntegralNarrow<std::uint32_t>(
        element_index, "HierarchicalTimeBinScheduler::importPersistentState element index"));
  }

  m_active_elements.clear();
  m_substep_open = false;
  m_diagnostics = {};
  m_diagnostics.occupancy_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  m_diagnostics.active_count_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  for (std::size_t bin = 0; bin < m_elements_by_bin.size(); ++bin) {
    m_diagnostics.occupancy_by_bin[bin] = checkedIntegralNarrow<std::uint32_t>(
        m_elements_by_bin[bin].size(), "scheduler bin occupancy");
  }
  refreshOwnedCapacityHighWater();
  validateInternalState("HierarchicalTimeBinScheduler::importPersistentState");
}

void HierarchicalTimeBinScheduler::validateInternalState(std::string_view source_label) const {
  if (m_representation == SchedulerRepresentation::kUniformRungZero) {
    if (m_max_bin != 0U || m_active_elements.size() != m_element_count ||
        !m_hot.bin_index.empty() || !m_hot.next_activation_tick.empty() ||
        !m_hot.active_flag.empty() || !m_hot.pending_bin_index.empty() ||
        !m_elements_by_bin.empty() || !m_position_in_bin.empty() ||
        !m_active_sort_scratch.empty() || !m_candidate_bin_index.empty() || !m_candidate_source.empty()) {
      throw std::runtime_error(schedulerContextMessage(
          "uniform rung-zero physical representation invariant violated", source_label, m_current_tick));
    }
    for (std::uint32_t i = 0; i < m_element_count; ++i) {
      if (m_active_elements[i] != i) {
        throw std::runtime_error(schedulerContextMessage(
            "uniform rung-zero active identity lane is not canonical", source_label, m_current_tick));
      }
    }
    return;
  }
  if (m_hot.size() != m_element_count) {
    throw std::runtime_error(schedulerContextMessage("generic scheduler element count mismatch", source_label, m_current_tick));
  }
  if (m_hot.bin_index.size() != m_hot.next_activation_tick.size() ||
      m_hot.bin_index.size() != m_hot.active_flag.size() ||
      m_hot.bin_index.size() != m_hot.pending_bin_index.size()) {
    throw std::runtime_error(schedulerContextMessage(
        "hot metadata arrays have mismatched sizes",
        source_label,
        m_current_tick));
  }
  if (m_position_in_bin.size() != m_hot.size()) {
    throw std::runtime_error(schedulerContextMessage(
        "bin-position mirror size does not match hot metadata",
        source_label,
        m_current_tick));
  }
  for (std::uint32_t element = 0; element < m_hot.size(); ++element) {
    const std::uint8_t bin = m_hot.bin_index[element];
    const std::uint64_t next_tick = m_hot.next_activation_tick[element];
    if (bin > m_max_bin) {
      throw std::runtime_error(timeBinContextMessage(
          "bin exceeds scheduler max_bin",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    if (m_hot.active_flag[element] > 1U) {
      throw std::runtime_error(timeBinContextMessage(
          "active flag is not boolean",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    if (!m_substep_open && m_hot.active_flag[element] != 0U) {
      throw std::runtime_error(timeBinContextMessage(
          "stale active flag outside an open substep",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    const std::uint64_t period = binPeriodTicks(bin);
    if (next_tick < m_current_tick) {
      throw std::runtime_error(timeBinContextMessage(
          "next activation tick is behind the global scheduler tick",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    if (next_tick > m_current_tick && next_tick % period != 0U) {
      throw std::runtime_error(timeBinContextMessage(
          "future next activation tick is not aligned to bin period",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    const std::uint8_t pending = m_hot.pending_bin_index[element];
    if (pending != k_unset_pending_bin && pending > m_max_bin) {
      throw std::runtime_error(timeBinContextMessage(
          "pending bin exceeds scheduler max_bin",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
  }
}

void HierarchicalTimeBinScheduler::validateTransitionRequest(
    std::uint32_t element_index,
    std::uint8_t target_bin,
    std::string_view source_label) const {
  if (element_index >= m_element_count) {
    throw std::out_of_range("element_index out of range");
  }
  if (m_representation == SchedulerRepresentation::kUniformRungZero && target_bin != 0U) {
    throw std::invalid_argument("uniform rung-zero scheduler rejects nonzero transition");
  }
  (void)target_bin;
  (void)source_label;
}

std::uint8_t HierarchicalTimeBinScheduler::clampBin(std::uint8_t requested) const noexcept {
  return std::min(requested, m_max_bin);
}

void HierarchicalTimeBinScheduler::eraseFromBin(std::uint32_t element_index, std::uint8_t bin_index) {
  auto& members = m_elements_by_bin[bin_index];
  const std::size_t remove_pos = m_position_in_bin[element_index];

  if (remove_pos >= members.size()) {
    throw std::runtime_error("bin position metadata corrupted");
  }

  const std::uint32_t last_element = members.back();
  members[remove_pos] = last_element;
  m_position_in_bin[last_element] = remove_pos;
  members.pop_back();
}

void HierarchicalTimeBinScheduler::insertIntoBin(std::uint32_t element_index, std::uint8_t bin_index) {
  auto& members = m_elements_by_bin[bin_index];
  m_position_in_bin[element_index] = members.size();
  members.push_back(element_index);
}

void HierarchicalTimeBinScheduler::rebuildActiveSet() {
  m_active_elements.clear();

  if (m_diagnostics.occupancy_by_bin.size() != m_elements_by_bin.size()) {
    m_diagnostics.occupancy_by_bin.assign(m_elements_by_bin.size(), 0U);
  }
  if (m_diagnostics.active_count_by_bin.size() != m_elements_by_bin.size()) {
    m_diagnostics.active_count_by_bin.assign(m_elements_by_bin.size(), 0U);
  }

  std::fill(m_diagnostics.active_count_by_bin.begin(), m_diagnostics.active_count_by_bin.end(), 0U);
  std::fill(m_diagnostics.occupancy_by_bin.begin(), m_diagnostics.occupancy_by_bin.end(), 0U);

  for (std::size_t bin = 0; bin < m_elements_by_bin.size(); ++bin) {
    const auto& members = m_elements_by_bin[bin];
    m_diagnostics.occupancy_by_bin[bin] = checkedIntegralNarrow<std::uint32_t>(
        members.size(), "scheduler active-set bin occupancy");

    if (!isBinActiveAtTick(static_cast<std::uint8_t>(bin), m_current_tick)) {
      continue;
    }

    for (const std::uint32_t element : members) {
      if (m_hot.next_activation_tick[element] != m_current_tick) {
        continue;
      }
      m_active_elements.push_back(element);
      m_hot.active_flag[element] = 1U;
      ++m_diagnostics.active_count_by_bin[bin];
    }
  }

  radixSortLocalIndices(m_active_elements, m_active_sort_scratch);

  m_diagnostics.active_elements = checkedIntegralNarrow<std::uint32_t>(
      m_active_elements.size(), "scheduler active element count");
  const auto total_elements = static_cast<double>(m_hot.size());
  m_diagnostics.active_fraction = total_elements > 0.0
      ? static_cast<double>(m_diagnostics.active_elements) / total_elements
      : 0.0;

  m_diagnostics.most_active_bin = 0;
  std::uint32_t best = 0;
  for (std::uint8_t bin = 0; bin <= m_max_bin; ++bin) {
    if (m_diagnostics.active_count_by_bin[bin] > best) {
      best = m_diagnostics.active_count_by_bin[bin];
      m_diagnostics.most_active_bin = bin;
    }
  }

  m_diagnostics.collapse_candidates = 0;
  if (m_hot.size() > 0) {
    const auto finest_occupancy = m_diagnostics.occupancy_by_bin[0];
    if (finest_occupancy * 4U >= m_hot.size() * 3U) {
      m_diagnostics.collapse_candidates = finest_occupancy;
    }
  }
}

void HierarchicalTimeBinScheduler::applyPendingTransitions() {
  for (const std::uint32_t element : m_active_elements) {
    const std::uint8_t pending = m_hot.pending_bin_index[element];
    if (pending == k_unset_pending_bin) {
      continue;
    }

    const std::uint8_t old_bin = m_hot.bin_index[element];
    const std::uint8_t new_bin = clampBin(pending);
    if (new_bin == old_bin) {
      m_hot.pending_bin_index[element] = k_unset_pending_bin;
      continue;
    }

    const std::uint64_t new_period = binPeriodTicks(new_bin);
    if (m_current_tick % new_period != 0) {
      // Coarsening to a larger period is legal only at a shared boundary. Keep
      // the request authoritative and pending while the element continues on
      // its current bin; the next aligned activation commits it.
      if (new_bin > old_bin) {
        ++m_diagnostics.deferred_coarsening_events;
        continue;
      }
      ++m_diagnostics.illegal_transition_attempts;
      throw std::runtime_error(timeBinContextMessage(
          "pending refinement is unexpectedly unaligned",
          "HierarchicalTimeBinScheduler::applyPendingTransitions",
          m_current_tick,
          element,
          old_bin,
          m_hot.next_activation_tick[element]));
    }

    eraseFromBin(element, old_bin);
    insertIntoBin(element, new_bin);
    m_hot.bin_index[element] = new_bin;

    if (new_bin < old_bin) {
      ++m_diagnostics.promoted_elements;
    } else {
      ++m_diagnostics.demoted_elements;
    }

    m_hot.pending_bin_index[element] = k_unset_pending_bin;
  }
}

TimeBinMappingResult mapDtToTimeBin(double dt_time_code, const TimeStepLimits& limits) {
  if (limits.min_dt_time_code <= 0.0 || limits.max_dt_time_code <= 0.0) {
    throw std::invalid_argument("timestep limits must be positive");
  }
  if (limits.max_dt_time_code < limits.min_dt_time_code) {
    throw std::invalid_argument("max_dt_time_code must be >= min_dt_time_code");
  }

  TimeBinMappingResult result{};
  const double clamped_dt = std::clamp(dt_time_code, limits.min_dt_time_code, limits.max_dt_time_code);
  result.clipped_to_min = clamped_dt == limits.min_dt_time_code && dt_time_code < limits.min_dt_time_code;
  result.clipped_to_max = clamped_dt == limits.max_dt_time_code && dt_time_code > limits.max_dt_time_code;

  std::uint8_t mapped_bin = 0;
  for (std::uint8_t bin = 0; bin <= limits.max_bin; ++bin) {
    const double dt_for_bin = binIndexToDt(bin, limits);
    if (dt_for_bin <= clamped_dt) {
      mapped_bin = bin;
    } else {
      break;
    }
  }
  result.bin_index = mapped_bin;
  return result;
}

double binIndexToDt(std::uint8_t bin_index, const TimeStepLimits& limits) {
  return limits.min_dt_time_code * static_cast<double>(powerOfTwo(std::min(bin_index, limits.max_bin)));
}


namespace {

void validatePersistentStateShapeForRemap(
    const TimeBinPersistentState& state,
    std::size_t expected_size,
    std::string_view caller) {
  if (state.bin_index.size() != expected_size || state.next_activation_tick.size() != expected_size ||
      state.active_flag.size() != expected_size || state.pending_bin_index.size() != expected_size) {
    throw std::invalid_argument(std::string(caller) + ": scheduler persistent state size does not match identity lane");
  }
}

[[nodiscard]] std::unordered_map<std::uint64_t, std::size_t> buildUniqueIdentityIndex(
    std::span<const std::uint64_t> ids,
    std::string_view caller) {
  std::unordered_map<std::uint64_t, std::size_t> index_by_id;
  index_by_id.reserve(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] == 0U) {
      throw std::invalid_argument(std::string(caller) + ": zero identity is not a valid scheduler remap key");
    }
    if (!index_by_id.emplace(ids[i], i).second) {
      throw std::invalid_argument(std::string(caller) + ": duplicate stable identity in scheduler remap keys");
    }
  }
  return index_by_id;
}

}  // namespace

std::vector<TimeBinSchedulerIdentityRecord> exportParticleSchedulerIdentityRecords(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> particle_indices) {
  if (scheduler.elementCount() != state.particles.size()) {
    throw std::invalid_argument("exportParticleSchedulerIdentityRecords: scheduler size mismatch");
  }
  std::vector<TimeBinSchedulerIdentityRecord> records;
  records.reserve(particle_indices.size());
  for (const std::uint32_t particle_index : particle_indices) {
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range("exportParticleSchedulerIdentityRecords: particle index out of range");
    }
    records.push_back(TimeBinSchedulerIdentityRecord{
        .element_id = state.particle_sidecar.particle_id[particle_index],
        .bin_index = scheduler.binIndex(particle_index),
        .next_activation_tick = scheduler.nextActivationTick(particle_index),
        .pending_bin_index = scheduler.pendingBinIndex(particle_index),
    });
  }
  return records;
}

void attachSchedulerFieldsToParticleMigrationRecords(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> particle_indices,
    std::span<ParticleMigrationRecord> records) {
  if (particle_indices.size() != records.size()) {
    throw std::invalid_argument(
        "attachSchedulerFieldsToParticleMigrationRecords: particle_indices and records sizes must match");
  }
  if (scheduler.elementCount() != state.particles.size()) {
    throw std::invalid_argument("attachSchedulerFieldsToParticleMigrationRecords: scheduler size mismatch");
  }
  for (std::size_t i = 0; i < particle_indices.size(); ++i) {
    const std::uint32_t particle_index = particle_indices[i];
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range("attachSchedulerFieldsToParticleMigrationRecords: particle index out of range");
    }
    if (records[i].particle_id != state.particle_sidecar.particle_id[particle_index]) {
      throw std::invalid_argument(
          "attachSchedulerFieldsToParticleMigrationRecords: record particle_id does not match source particle index");
    }
    records[i].has_scheduler_fields = true;
    records[i].scheduler_fields.bin_index = scheduler.binIndex(particle_index);
    records[i].scheduler_fields.next_activation_tick = scheduler.nextActivationTick(particle_index);
    records[i].scheduler_fields.pending_bin_index = scheduler.pendingBinIndex(particle_index);
    records[i].time_bin = scheduler.binIndex(particle_index);
  }
}


std::vector<ParticleMigrationRecord> SimulationState::packParticleMigrationRecords(
    std::span<const std::uint32_t> local_indices,
    const HierarchicalTimeBinScheduler& scheduler) const {
  auto records = packParticleMigrationRecordsCore(local_indices);
  attachSchedulerFieldsToParticleMigrationRecords(scheduler, *this, local_indices, records);
  return records;
}

TimeBinPersistentState rebuildSchedulerPersistentStateFromMigrationRecords(
    std::uint64_t current_tick,
    std::uint8_t max_bin,
    std::span<const ParticleMigrationRecord> records,
    std::span<const std::uint64_t> destination_element_ids) {
  std::vector<TimeBinSchedulerIdentityRecord> identity_records;
  identity_records.reserve(records.size());
  for (const ParticleMigrationRecord& record : records) {
    if (!record.has_scheduler_fields) {
      throw std::invalid_argument(
          "rebuildSchedulerPersistentStateFromMigrationRecords: migration record is missing scheduler authority fields");
    }
    if (record.scheduler_fields.bin_index != record.time_bin) {
      throw std::invalid_argument(
          "rebuildSchedulerPersistentStateFromMigrationRecords: time_bin mirror is stale relative to scheduler fields");
    }
    identity_records.push_back(TimeBinSchedulerIdentityRecord{
        .element_id = record.particle_id,
        .bin_index = record.scheduler_fields.bin_index,
        .next_activation_tick = record.scheduler_fields.next_activation_tick,
        .pending_bin_index = record.scheduler_fields.pending_bin_index,
    });
  }
  return rebuildSchedulerPersistentStateFromIdentityRecords(
      current_tick, max_bin, identity_records, destination_element_ids);
}

void rebuildSchedulerFromParticleMigrationRecords(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const ParticleMigrationRecord> records,
    std::span<const std::uint64_t> destination_particle_ids) {
  scheduler.importPersistentState(rebuildSchedulerPersistentStateFromMigrationRecords(
      scheduler.currentTick(),
      scheduler.maxBin(),
      records,
      destination_particle_ids));
}

TimeBinPersistentState remapSchedulerPersistentStateByParticleId(
    const TimeBinPersistentState& source_state,
    std::span<const std::uint64_t> source_particle_ids,
    std::span<const std::uint64_t> destination_particle_ids) {
  validatePersistentStateShapeForRemap(
      source_state,
      source_particle_ids.size(),
      "remapSchedulerPersistentStateByParticleId");
  const auto source_index_by_id =
      buildUniqueIdentityIndex(source_particle_ids, "remapSchedulerPersistentStateByParticleId");
  (void)buildUniqueIdentityIndex(
      destination_particle_ids,
      "remapSchedulerPersistentStateByParticleId");

  TimeBinPersistentState remapped;
  remapped.current_tick = source_state.current_tick;
  remapped.max_bin = source_state.max_bin;
  remapped.bin_index.resize(destination_particle_ids.size());
  remapped.next_activation_tick.resize(destination_particle_ids.size());
  remapped.active_flag.assign(destination_particle_ids.size(), 0U);
  remapped.pending_bin_index.resize(destination_particle_ids.size());

  for (std::size_t new_index = 0; new_index < destination_particle_ids.size(); ++new_index) {
    const auto found = source_index_by_id.find(destination_particle_ids[new_index]);
    if (found == source_index_by_id.end()) {
      throw std::invalid_argument(
          "remapSchedulerPersistentStateByParticleId: missing scheduler authority for destination particle_id");
    }
    const std::size_t old_index = found->second;
    remapped.bin_index[new_index] = source_state.bin_index[old_index];
    remapped.next_activation_tick[new_index] = source_state.next_activation_tick[old_index];
    remapped.pending_bin_index[new_index] = source_state.pending_bin_index[old_index];
  }
  return remapped;
}

void remapSchedulerByParticleId(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const std::uint64_t> source_particle_ids,
    std::span<const std::uint64_t> destination_particle_ids) {
  scheduler.importPersistentState(remapSchedulerPersistentStateByParticleId(
      scheduler.exportPersistentState(),
      source_particle_ids,
      destination_particle_ids));
}

void remapSchedulerByParticleReorderMap(
    HierarchicalTimeBinScheduler& scheduler,
    const ParticleReorderMap& reorder_map) {
  const auto source_state = scheduler.exportPersistentState();
  if (!reorder_map.isConsistent(source_state.bin_index.size())) {
    throw std::invalid_argument("remapSchedulerByParticleReorderMap: inconsistent reorder map for scheduler");
  }
  TimeBinPersistentState remapped;
  remapped.current_tick = source_state.current_tick;
  remapped.max_bin = source_state.max_bin;
  remapped.bin_index.resize(source_state.bin_index.size());
  remapped.next_activation_tick.resize(source_state.bin_index.size());
  remapped.active_flag.assign(source_state.bin_index.size(), 0U);
  remapped.pending_bin_index.resize(source_state.bin_index.size());
  for (std::size_t new_index = 0; new_index < reorder_map.new_to_old_index.size(); ++new_index) {
    const std::uint32_t old_index = reorder_map.new_to_old_index[new_index];
    remapped.bin_index[new_index] = source_state.bin_index[old_index];
    remapped.next_activation_tick[new_index] = source_state.next_activation_tick[old_index];
    remapped.pending_bin_index[new_index] = source_state.pending_bin_index[old_index];
  }
  scheduler.importPersistentState(remapped);
}

TimeBinPersistentState rebuildSchedulerPersistentStateFromIdentityRecords(
    std::uint64_t current_tick,
    std::uint8_t max_bin,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    std::span<const std::uint64_t> destination_element_ids) {
  std::unordered_map<std::uint64_t, TimeBinSchedulerIdentityRecord> record_by_id;
  record_by_id.reserve(records.size());
  for (const auto& record : records) {
    if (record.element_id == 0U) {
      throw std::invalid_argument("rebuildSchedulerPersistentStateFromIdentityRecords: zero element_id is invalid");
    }
    if (record.bin_index > max_bin) {
      throw std::invalid_argument("rebuildSchedulerPersistentStateFromIdentityRecords: bin_index exceeds max_bin");
    }
    if (record.pending_bin_index != HierarchicalTimeBinScheduler::k_unset_pending_bin &&
        record.pending_bin_index > max_bin) {
      throw std::invalid_argument("rebuildSchedulerPersistentStateFromIdentityRecords: pending_bin_index exceeds max_bin");
    }
    if (!record_by_id.emplace(record.element_id, record).second) {
      throw std::invalid_argument("rebuildSchedulerPersistentStateFromIdentityRecords: duplicate scheduler identity record");
    }
  }
  (void)buildUniqueIdentityIndex(destination_element_ids, "rebuildSchedulerPersistentStateFromIdentityRecords");

  TimeBinPersistentState rebuilt;
  rebuilt.current_tick = current_tick;
  rebuilt.max_bin = max_bin;
  rebuilt.bin_index.resize(destination_element_ids.size());
  rebuilt.next_activation_tick.resize(destination_element_ids.size());
  rebuilt.active_flag.assign(destination_element_ids.size(), 0U);
  rebuilt.pending_bin_index.resize(destination_element_ids.size());
  for (std::size_t row = 0; row < destination_element_ids.size(); ++row) {
    const auto found = record_by_id.find(destination_element_ids[row]);
    if (found == record_by_id.end()) {
      throw std::invalid_argument(
          "rebuildSchedulerPersistentStateFromIdentityRecords: ParticleMigrationRecord::time_bin mirror is not scheduler authority");
    }
    rebuilt.bin_index[row] = found->second.bin_index;
    rebuilt.next_activation_tick[row] = found->second.next_activation_tick;
    rebuilt.pending_bin_index[row] = found->second.pending_bin_index;
  }
  return rebuilt;
}

void rebuildSchedulerFromParticleIdentityRecords(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    std::span<const std::uint64_t> destination_particle_ids) {
  scheduler.importPersistentState(rebuildSchedulerPersistentStateFromIdentityRecords(
      scheduler.currentTick(),
      scheduler.maxBin(),
      records,
      destination_particle_ids));
}

std::vector<TimeBinSchedulerIdentityRecord> exportGasCellSchedulerIdentityRecords(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> cell_indices) {
  state.requireGasCellIdentityMapCoversDenseRows("exportGasCellSchedulerIdentityRecords");
  if (scheduler.elementCount() != state.cells.size()) {
    throw std::invalid_argument("exportGasCellSchedulerIdentityRecords: scheduler size mismatch");
  }
  std::vector<TimeBinSchedulerIdentityRecord> records;
  records.reserve(cell_indices.size());
  for (const std::uint32_t cell_row : cell_indices) {
    if (cell_row >= state.cells.size()) {
      throw std::out_of_range("exportGasCellSchedulerIdentityRecords: cell row out of range");
    }
    const auto gas_cell_id = state.gas_cell_identity.gasCellIdForLocalRow(cell_row);
    if (!gas_cell_id.has_value()) {
      throw std::runtime_error("exportGasCellSchedulerIdentityRecords: identity map is missing cell row");
    }
    records.push_back(TimeBinSchedulerIdentityRecord{
        .element_id = *gas_cell_id,
        .bin_index = scheduler.binIndex(cell_row),
        .next_activation_tick = scheduler.nextActivationTick(cell_row),
        .pending_bin_index = scheduler.pendingBinIndex(cell_row),
    });
  }
  return records;
}

void rebuildSchedulerFromGasCellIdentityRecords(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    const SimulationState& state) {
  state.requireGasCellIdentityMapCoversDenseRows("rebuildSchedulerFromGasCellIdentityRecords");
  std::vector<std::uint64_t> destination_gas_cell_ids(state.cells.size());
  for (std::uint32_t cell_row = 0; cell_row < state.cells.size(); ++cell_row) {
    const auto gas_cell_id = state.gas_cell_identity.gasCellIdForLocalRow(cell_row);
    if (!gas_cell_id.has_value()) {
      throw std::runtime_error("rebuildSchedulerFromGasCellIdentityRecords: identity map is missing cell row");
    }
    destination_gas_cell_ids[cell_row] = *gas_cell_id;
  }
  scheduler.importPersistentState(rebuildSchedulerPersistentStateFromIdentityRecords(
      scheduler.currentTick(), scheduler.maxBin(), records, destination_gas_cell_ids));
}

void syncGasCellTimeBinMirrorsFromGasCellScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state) {
  if (state.cells.size() == 0U) {
    return;
  }
  state.requireGasCellIdentityMapCoversDenseRows("syncGasCellTimeBinMirrorsFromGasCellScheduler");
  if (scheduler.elementCount() != state.cells.size()) {
    throw std::invalid_argument("syncGasCellTimeBinMirrorsFromGasCellScheduler: scheduler size mismatch");
  }
  for (std::size_t cell_row = 0; cell_row < state.cells.size(); ++cell_row) {
    state.cells.time_bin[cell_row] = scheduler.binIndex(checkedLocalCellRow(cell_row, "sync gas scheduler row"));
  }
}

void syncGasCellTimeBinMirrorsFromParticleScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state) {
  if (state.cells.size() == 0) {
    return;
  }
  state.requireGasCellIdentityMapCoversDenseRows("syncGasCellTimeBinMirrorsFromParticleScheduler");
  if (scheduler.elementCount() < state.particles.size()) {
    throw std::invalid_argument("syncGasCellTimeBinMirrorsFromParticleScheduler: scheduler lacks particle bin entries");
  }
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const auto parent_id = parentParticleIdForGasCellRow(state, cell_index);
    if (!parent_id.has_value()) {
      continue;
    }
    const auto particle_it =
        std::find(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end(), *parent_id);
    if (particle_it == state.particle_sidecar.particle_id.end()) {
      throw std::invalid_argument("syncGasCellTimeBinMirrorsFromParticleScheduler: parent particle is not local");
    }
    const auto particle_distance =
        std::distance(state.particle_sidecar.particle_id.begin(), particle_it);
    const std::uint32_t particle_index = checkedLocalParticleRow(
        static_cast<std::size_t>(particle_distance),
        "syncGasCellTimeBinMirrorsFromParticleScheduler parent row");
    state.cells.time_bin[cell_index] = scheduler.binIndex(particle_index);
  }
}

void syncTimeBinMirrorsFromScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state,
    TimeBinMirrorDomain domain) {
  auto sync_particles = [&]() {
    if (state.hasHomogeneousDmoMetadata() && scheduler.maxBin() == 0U && state.particles.time_bin.empty()) return;
    // A heterogeneous scheduler needs only the existing one-byte bin mirror;
    // uniform species/ownership/drift metadata remains physically compact.
    if (state.hasHomogeneousDmoMetadata() && state.particles.time_bin.empty()) {
      state.particles.time_bin.resize(state.particles.size(), 0U);
    }
    if (scheduler.elementCount() < state.particles.size()) {
      throw std::invalid_argument("syncTimeBinMirrorsFromScheduler: scheduler lacks particle bin entries");
    }
    for (std::size_t i = 0; i < state.particles.size(); ++i) {
      state.particles.time_bin[i] = scheduler.binIndex(checkedLocalParticleRow(i, "sync particle scheduler row"));
    }
  };
  auto sync_cells = [&]() {
    syncGasCellTimeBinMirrorsFromParticleScheduler(scheduler, state);
  };

  switch (domain) {
    case TimeBinMirrorDomain::kParticles:
      sync_particles();
      break;
    case TimeBinMirrorDomain::kCells:
      sync_cells();
      break;
    case TimeBinMirrorDomain::kParticlesAndCells:
      sync_particles();
      sync_cells();
      break;
  }
}

bool timeBinMirrorsMatchScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain) {
  auto particles_match = [&]() {
    if (scheduler.elementCount() < state.particles.size()) {
      return false;
    }
    for (std::size_t i = 0; i < state.particles.size(); ++i) {
      if (state.particleTimeBin(i) != scheduler.binIndex(checkedLocalParticleRow(i, "timeBinMirrorsMatchScheduler particle row"))) {
        return false;
      }
    }
    return true;
  };
  auto cells_match = [&]() {
    if (state.cells.size() == 0) {
      return true;
    }
    if (!state.gas_cell_identity.isConsistent() || !state.gas_cell_identity.coversDenseLocalRows(state.cells.size())) {
      return false;
    }
    if (scheduler.elementCount() < state.particles.size()) {
      return false;
    }
    for (std::uint32_t i = 0; i < state.cells.size(); ++i) {
      const auto parent_id = parentParticleIdForGasCellRow(state, i);
      if (!parent_id.has_value()) {
        continue;
      }
      const auto particle_it =
          std::find(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end(), *parent_id);
      if (particle_it == state.particle_sidecar.particle_id.end()) {
        return false;
      }
      const auto particle_distance =
          std::distance(state.particle_sidecar.particle_id.begin(), particle_it);
      const std::uint32_t particle_index = checkedLocalParticleRow(
          static_cast<std::size_t>(particle_distance),
          "timeBinMirrorsMatchScheduler parent row");
      if (state.cells.time_bin[i] != scheduler.binIndex(particle_index)) {
        return false;
      }
    }
    return true;
  };

  switch (domain) {
    case TimeBinMirrorDomain::kParticles:
      return particles_match();
    case TimeBinMirrorDomain::kCells:
      return cells_match();
    case TimeBinMirrorDomain::kParticlesAndCells:
      return particles_match() && cells_match();
  }
  return false;
}

void debugAssertTimeBinMirrorAuthorityInvariant(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain) {
  const auto throw_particle = [&](std::uint32_t index) {
    throw std::runtime_error(timeBinContextMessage(
        "time_bin mirror authority invariant violated for particle mirror",
        "debugAssertTimeBinMirrorAuthorityInvariant",
        scheduler.currentTick(),
        index,
        index < scheduler.elementCount() ? scheduler.binIndex(index) : 0,
        index < scheduler.elementCount() ? scheduler.nextActivationTick(index) : scheduler.currentTick()));
  };
  const auto throw_cell = [&](std::uint32_t index) {
    throw std::runtime_error(timeBinContextMessage(
        "time_bin mirror authority invariant violated for cell mirror",
        "debugAssertTimeBinMirrorAuthorityInvariant",
        scheduler.currentTick(),
        index,
        index < scheduler.elementCount() ? scheduler.binIndex(index) : 0,
        index < scheduler.elementCount() ? scheduler.nextActivationTick(index) : scheduler.currentTick()));
  };
  if ((domain == TimeBinMirrorDomain::kParticles || domain == TimeBinMirrorDomain::kParticlesAndCells) &&
      scheduler.elementCount() < state.particles.size()) {
    throw_particle(checkedLocalParticleRow(
        scheduler.elementCount(), "debugAssertTimeBinMirrorAuthorityInvariant particle scheduler size"));
  }
  if ((domain == TimeBinMirrorDomain::kCells || domain == TimeBinMirrorDomain::kParticlesAndCells) &&
      scheduler.elementCount() < state.cells.size()) {
    throw_cell(checkedLocalCellRow(
        scheduler.elementCount(), "debugAssertTimeBinMirrorAuthorityInvariant cell scheduler size"));
  }
  if (domain == TimeBinMirrorDomain::kParticles || domain == TimeBinMirrorDomain::kParticlesAndCells) {
    for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
      if (state.particleTimeBin(i) != scheduler.binIndex(checkedLocalParticleRow(i, "timeBinMirrorsMatchScheduler particle row"))) {
        throw_particle(i);
      }
    }
  }
  if (domain == TimeBinMirrorDomain::kCells || domain == TimeBinMirrorDomain::kParticlesAndCells) {
    if (state.cells.size() > 0) {
      requireGasCellIdentityMapCoversDenseRows(state, "debugAssertTimeBinMirrorAuthorityInvariant");
    }
    for (std::uint32_t i = 0; i < state.cells.size(); ++i) {
      const std::uint32_t scheduler_index = checkedLocalCellRow(
          i, "debugAssertTimeBinMirrorAuthorityInvariant cell row");
      if (state.cells.time_bin[i] != scheduler.binIndex(scheduler_index)) {
        throw_cell(i);
      }
    }
  }
}

double computeCflTimeStep(const CflTimeStepInput& input, double c_cfl) {
  if (input.cell_width_code <= 0.0) {
    throw std::invalid_argument("cell_width_code must be positive");
  }
  if (c_cfl <= 0.0) {
    throw std::invalid_argument("c_cfl must be positive");
  }

  const double denom = std::abs(input.flow_speed_code) + std::max(input.sound_speed_code, 0.0);
  if (denom <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  return c_cfl * (input.cell_width_code / denom);
}

std::array<bool, 3> activeHydroAxes(
    const std::array<double, 3>& cell_width_axis_code,
    const std::array<double, 3>& domain_extent_axis_code) {
  std::array<bool, 3> active{true, true, true};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const double extent = domain_extent_axis_code[axis];
    if (std::isfinite(extent) && extent > 0.0 &&
        cell_width_axis_code[axis] >= extent * (1.0 - 1.0e-12)) {
      active[axis] = false;
    }
  }
  return active;
}

HydroCflStabilityRule hydroCflStabilityRuleFor(const SimulationConfig& config) noexcept {
  return config.numerics.hydro_integration_scheme == HydroIntegrationMode::kCornerTransportUpwind
      ? HydroCflStabilityRule::kPerAxisMaximum
      : HydroCflStabilityRule::kDimensionSummed;
}

std::array<double, 3> hydroDomainExtentCode(const SimulationConfig& config) noexcept {
  return {config.cosmology.box_size_x_mpc_comoving,
          config.cosmology.box_size_y_mpc_comoving,
          config.cosmology.box_size_z_mpc_comoving};
}

double computeDirectionalCflTimeStep(const DirectionalCflTimeStepInput& input, double c_cfl) {
  if (c_cfl <= 0.0 || !std::isfinite(c_cfl)) {
    throw std::invalid_argument("c_cfl must be finite and positive");
  }
  if (!std::isfinite(input.sound_speed_code)) {
    throw std::invalid_argument("sound_speed_code must be finite");
  }
  if (!std::isfinite(input.scale_factor) || input.scale_factor <= 0.0) {
    throw std::invalid_argument("directional CFL scale_factor must be finite and positive");
  }
  const double length_scale = input.coordinate_frame == CoordinateFrame::kComoving
      ? input.scale_factor
      : 1.0;
  const double sound_speed = std::max(input.sound_speed_code, 0.0);
  double dt = std::numeric_limits<double>::infinity();
  double summed_rate = 0.0;
  for (std::size_t axis = 0; axis < input.cell_width_axis_code.size(); ++axis) {
    const double cell_width_stored = input.cell_width_axis_code[axis];
    if (cell_width_stored <= 0.0 || !std::isfinite(cell_width_stored)) {
      throw std::invalid_argument("cell_width_axis_code entries must be finite and positive");
    }
    if (!std::isfinite(input.velocity_axis_code[axis])) {
      throw std::invalid_argument("velocity_axis_code entries must be finite");
    }
    if (!input.active_axis[axis]) {
      continue;
    }
    const double denom = std::abs(input.velocity_axis_code[axis]) + sound_speed;
    if (denom > 0.0) {
      const double crossing_length_code = length_scale * cell_width_stored;
      dt = std::min(dt, c_cfl * crossing_length_code / denom);
      summed_rate += denom / crossing_length_code;
    }
  }
  if (input.stability_rule == HydroCflStabilityRule::kDimensionSummed) {
    return summed_rate > 0.0 ? c_cfl / summed_rate : std::numeric_limits<double>::infinity();
  }
  return dt;
}

HydroCflDiagnostics makeHydroCflDiagnostics(
    std::uint32_t local_row,
    const DirectionalCflTimeStepInput& input,
    double c_cfl,
    double accepted_dt_time_code,
    std::uint64_t gas_cell_id,
    std::optional<std::uint64_t> patch_id,
    std::optional<std::uint32_t> patch_row) {
  HydroCflDiagnostics diagnostics;
  diagnostics.local_row = local_row;
  diagnostics.gas_cell_id = gas_cell_id;
  diagnostics.has_gas_cell_id = gas_cell_id != 0U;
  diagnostics.patch_id = patch_id.value_or(0U);
  diagnostics.has_patch_id = patch_id.has_value();
  diagnostics.patch_row = patch_row.value_or(0U);
  diagnostics.has_patch_row = patch_row.has_value();
  diagnostics.accepted_dt_time_code = accepted_dt_time_code;
  diagnostics.cfl_number = c_cfl;
  diagnostics.cell_width_axis_code = input.cell_width_axis_code;
  diagnostics.velocity_axis_code = input.velocity_axis_code;
  diagnostics.sound_speed_code = input.sound_speed_code;
  diagnostics.coordinate_frame = input.coordinate_frame;
  diagnostics.scale_factor = input.scale_factor;
  diagnostics.proposed_dt_time_code = computeDirectionalCflTimeStep(input, c_cfl);
  diagnostics.safety_factor = diagnostics.accepted_dt_time_code > 0.0
      ? diagnostics.proposed_dt_time_code / diagnostics.accepted_dt_time_code
      : std::numeric_limits<double>::infinity();

  double limiting_dt = std::numeric_limits<double>::infinity();
  const double sound_speed = std::max(input.sound_speed_code, 0.0);
  const double length_scale = input.coordinate_frame == CoordinateFrame::kComoving
      ? input.scale_factor
      : 1.0;
  diagnostics.stability_rule = input.stability_rule;
  for (std::size_t axis = 0; axis < input.cell_width_axis_code.size(); ++axis) {
    if (!input.active_axis[axis]) {
      continue;
    }
    const double denom = std::abs(input.velocity_axis_code[axis]) + sound_speed;
    const double axis_dt = denom > 0.0
        ? c_cfl * length_scale * input.cell_width_axis_code[axis] / denom
        : std::numeric_limits<double>::infinity();
    if (axis_dt < limiting_dt) {
      limiting_dt = axis_dt;
      diagnostics.limiting_axis = static_cast<std::uint8_t>(axis);
    }
  }
  return diagnostics;
}

void assertHydroCflStable(
    const HydroCflDiagnostics& diagnostics,
    double relative_tolerance) {
  if (relative_tolerance < 0.0 || !std::isfinite(relative_tolerance)) {
    throw std::invalid_argument("relative_tolerance must be finite and non-negative");
  }
  if (diagnostics.accepted_dt_time_code <= 0.0 || !std::isfinite(diagnostics.accepted_dt_time_code)) {
    throw std::invalid_argument("accepted hydro dt must be finite and positive");
  }
  if (!std::isfinite(diagnostics.proposed_dt_time_code)) {
    return;
  }
  const double allowed = diagnostics.proposed_dt_time_code *
      (1.0 + relative_tolerance);
  if (diagnostics.accepted_dt_time_code <= allowed) {
    return;
  }

  std::ostringstream out;
  out << "hydro CFL violation before state mutation"
      << "; local_row=" << diagnostics.local_row
      << "; gas_cell_id=" << (diagnostics.has_gas_cell_id ? std::to_string(diagnostics.gas_cell_id) : "n/a")
      << "; patch_id=" << (diagnostics.has_patch_id ? std::to_string(diagnostics.patch_id) : "n/a")
      << "; patch_row=" << (diagnostics.has_patch_row ? std::to_string(diagnostics.patch_row) : "n/a")
      << "; proposed_dt=" << diagnostics.proposed_dt_time_code
      << "; accepted_dt=" << diagnostics.accepted_dt_time_code
      << "; cfl_number=" << diagnostics.cfl_number
      << "; safety_factor=" << diagnostics.safety_factor
      << "; limiting_axis=" << static_cast<unsigned>(diagnostics.limiting_axis)
      << "; velocity=(" << diagnostics.velocity_axis_code[0] << ","
      << diagnostics.velocity_axis_code[1] << ","
      << diagnostics.velocity_axis_code[2] << ")"
      << "; sound_speed=" << diagnostics.sound_speed_code
      << "; coordinate_frame="
      << (diagnostics.coordinate_frame == CoordinateFrame::kComoving ? "comoving" : "physical")
      << "; scale_factor=" << diagnostics.scale_factor
      << "; cell_width=(" << diagnostics.cell_width_axis_code[0] << ","
      << diagnostics.cell_width_axis_code[1] << ","
      << diagnostics.cell_width_axis_code[2] << ")";
  throw std::runtime_error(out.str());
}

double computeGravityTimeStep(const GravityTimeStepInput& input, double eta) {
  if (!std::isfinite(input.softening_length_code) || input.softening_length_code <= 0.0) {
    throw std::invalid_argument("softening_length_code must be finite and positive");
  }
  if (!std::isfinite(input.acceleration_magnitude_code)) {
    throw std::invalid_argument("acceleration_magnitude_code must be finite");
  }
  if (!std::isfinite(eta) || eta <= 0.0) {
    throw std::invalid_argument("eta must be finite and positive");
  }

  const double accel = std::abs(input.acceleration_magnitude_code);
  if (accel == 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  return eta * std::sqrt(input.softening_length_code / accel);
}

double computeComovingGravityTimeStep(
    const ComovingGravityTimeStepInput& input,
    double eta) {
  if (!std::isfinite(input.scale_factor) || input.scale_factor <= 0.0) {
    throw std::invalid_argument("comoving gravity timestep scale_factor must be finite and positive");
  }
  if (!std::isfinite(input.scale_free_acceleration_magnitude_code)) {
    throw std::invalid_argument(
        "comoving gravity timestep scale-free acceleration must be finite");
  }
  const double scale_factor_cubed =
      input.scale_factor * input.scale_factor * input.scale_factor;
  return computeGravityTimeStep(
      {.softening_length_code = input.softening_length_comoving_code,
       .acceleration_magnitude_code =
           input.scale_free_acceleration_magnitude_code / scale_factor_cubed},
      eta);
}

double computeComovingDisplacementTimeStep(const ComovingDisplacementTimeStepInput& input) {
  if (!std::isfinite(input.mesh_or_split_length_comoving_code) ||
      input.mesh_or_split_length_comoving_code <= 0.0 ||
      !std::isfinite(input.velocity_magnitude_peculiar_code) ||
      input.velocity_magnitude_peculiar_code < 0.0 ||
      !std::isfinite(input.scale_free_acceleration_magnitude_code) ||
      input.scale_free_acceleration_magnitude_code < 0.0 ||
      !std::isfinite(input.scale_factor) || input.scale_factor <= 0.0) {
    throw std::invalid_argument("comoving displacement timestep requires finite positive length/scale and nonnegative speed/acceleration");
  }
  // At the synchronized epoch v_com=u/a, g_com=A/a^3. In an expanding
  // background these frozen-epoch denominators are conservative. Solve
  // v_com dt + g_com dt^2/2 <= ell without cancellation. This is a local
  // resolution-crossing estimate, not a bound on subsequent force growth.
  const long double a = input.scale_factor;
  const long double v = input.velocity_magnitude_peculiar_code/a;
  const long double g = ((input.scale_free_acceleration_magnitude_code/a)/a)/a;
  const long double ell = input.mesh_or_split_length_comoving_code;
  if (v == 0.0L && g == 0.0L) return std::numeric_limits<double>::infinity();
  return static_cast<double>(2.0L*ell/(v + std::hypot(v,std::sqrt(2.0L*g*ell))));
}

double combineTimeStepCriteria(
    std::uint32_t element_index,
    const TimeStepCriteriaHooks& hooks,
    double fallback_dt_time_code) {
  if (!std::isfinite(fallback_dt_time_code) || fallback_dt_time_code <= 0.0) {
    throw std::invalid_argument("fallback_dt_time_code must be finite and positive");
  }

  double dt = fallback_dt_time_code;
  if (hooks.cfl_hook) {
    dt = std::min(dt, requireValidTimeStepCriterion(hooks.cfl_hook(element_index), "cfl", element_index));
  }
  if (hooks.gravity_hook) {
    dt = std::min(dt, requireValidTimeStepCriterion(hooks.gravity_hook(element_index), "gravity", element_index));
  }
  if (hooks.source_hook) {
    dt = std::min(dt, requireValidTimeStepCriterion(hooks.source_hook(element_index), "source", element_index));
  }
  if (hooks.user_clamp_hook) {
    dt = std::min(dt, requireValidTimeStepCriterion(hooks.user_clamp_hook(element_index), "user_clamp", element_index));
  }

  if (!std::isfinite(dt) || dt <= 0.0) {
    throw std::runtime_error("combineTimeStepCriteria produced a non-finite or non-positive accepted timestep");
  }
  return dt;
}

double computeScaleFactorRate(const LambdaCdmBackground& background, double scale_factor) {
  if (!std::isfinite(scale_factor) || scale_factor <= 0.0) {
    throw std::invalid_argument("scale_factor must be finite and positive");
  }
  const double rate = scale_factor * background.hubbleSi(scale_factor);
  if (!std::isfinite(rate) || rate <= 0.0) {
    throw std::runtime_error("computed scale-factor rate must be finite and positive");
  }
  return rate;
}

double computeCosmologyExpansionTimeStep(
    const LambdaCdmBackground& background,
    double scale_factor,
    double max_delta_ln_a,
    double max_hubble_time_fraction,
    double time_si_per_code) {
  if (!std::isfinite(scale_factor) || scale_factor <= 0.0) {
    throw std::invalid_argument("scale_factor must be finite and positive");
  }
  if (!std::isfinite(max_delta_ln_a) || max_delta_ln_a <= 0.0 ||
      !std::isfinite(max_hubble_time_fraction) || max_hubble_time_fraction <= 0.0 ||
      !std::isfinite(time_si_per_code) || time_si_per_code <= 0.0) {
    throw std::invalid_argument("cosmology timestep limits and time_si_per_code must be finite and positive");
  }
  const double hubble = background.hubbleSi(scale_factor);
  if (!std::isfinite(hubble) || hubble <= 0.0) {
    throw std::invalid_argument("H(a) must be finite and positive");
  }
  const double dt_si = std::min(max_delta_ln_a / hubble, max_hubble_time_fraction / hubble);
  return dt_si / time_si_per_code;
}

double advanceScaleFactorByCosmicTime(
    const LambdaCdmBackground& background,
    double scale_factor,
    double dt_time_si,
    std::uint32_t midpoint_samples) {
  if (!std::isfinite(scale_factor) || scale_factor <= 0.0) {
    throw std::invalid_argument("scale_factor must be finite and positive");
  }
  if (!std::isfinite(dt_time_si) || dt_time_si < 0.0) {
    throw std::invalid_argument("dt_time_si must be finite and non-negative");
  }
  if (dt_time_si == 0.0) {
    return scale_factor;
  }

  double lo = scale_factor;
  double hi = scale_factor * 1.001;
  if (hi <= lo) {
    hi = lo + 1.0e-12;
  }

  const std::uint32_t samples = std::max<std::uint32_t>(midpoint_samples, 8U);
  for (std::uint32_t iter = 0; iter < 256U; ++iter) {
    const double covered = midpointIntegrateCosmicTime(background, scale_factor, hi, samples);
    if (covered >= dt_time_si) {
      break;
    }
    lo = hi;
    hi *= 2.0;
    if (!std::isfinite(hi) || hi > 1.0e12) {
      throw std::runtime_error("failed to bracket FLRW scale-factor advance");
    }
  }

  if (midpointIntegrateCosmicTime(background, scale_factor, hi, samples) < dt_time_si) {
    throw std::runtime_error("failed to bracket FLRW scale-factor advance");
  }

  lo = scale_factor;
  for (std::uint32_t iter = 0; iter < 96U; ++iter) {
    const double mid = 0.5 * (lo + hi);
    const double covered = midpointIntegrateCosmicTime(background, scale_factor, mid, samples);
    if (covered < dt_time_si) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

double estimateDeltaTimeFromScaleFactorStep(
    const LambdaCdmBackground& background,
    double scale_factor,
    double delta_scale_factor) {
  if (delta_scale_factor < 0.0) {
    throw std::invalid_argument("delta_scale_factor must be non-negative");
  }
  return midpointIntegrateCosmicTime(background, scale_factor, scale_factor + delta_scale_factor, 64);
}

double computeComovingDriftFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples) {
  return midpointIntegrateDriftLike(background, scale_factor_begin, scale_factor_end, midpoint_samples);
}

double computeComovingKickFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples) {
  if (scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale_factor_end must be positive for peculiar-velocity kick integral");
  }

  // The runtime velocity arrays store physical peculiar velocities u = a dx/dt,
  // while TreePM returns the scale-free comoving Newtonian acceleration kernel A.
  // Holding A fixed over a kick substep gives
  //   du/dt + H u = A/a^2,
  // hence u_1 = (a_0/a_1) u_0 + A/a_1 * integral(dt/a).
  // The Hubble-drag factor is tracked separately; this returns only the A
  // multiplier.  For short steps it reduces to dt/a^2 in code units.
  return midpointIntegrateDriftLike(background, scale_factor_begin, scale_factor_end, midpoint_samples) /
      scale_factor_end;
}

double computeHubbleDragFactor(double scale_factor_begin, double scale_factor_end) {
  if (scale_factor_begin <= 0.0 || scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale factors must be positive");
  }
  if (scale_factor_end < scale_factor_begin) {
    throw std::invalid_argument("scale_factor_end must be >= scale_factor_begin");
  }

  return scale_factor_begin / scale_factor_end;
}

StepBoundaryState classifyStepBoundary(
    const SimulationState& state,
    const ActiveSetDescriptor& active_set,
    bool scheduler_owned_substep,
    StepBoundaryKind requested_kind) {
  const bool particle_subset = active_set.hasParticleSubset(state.particles.size());
  const bool cell_subset = active_set.hasCellSubset(state.cells.size());
  const bool local = scheduler_owned_substep &&
      (active_set.has_global_synchronization_metadata
           ? !active_set.globally_complete_active_set
           : (particle_subset || cell_subset));
  StepBoundaryState boundary{};
  boundary.local_substep = local;
  if (local) {
    boundary.kind = StepBoundaryKind::kLocalActiveBinStep;
    boundary.restart_safe = false;
    boundary.output_safe = false;
    // Local active-bin ownership does not by itself make an MPI collective
    // unsafe.  A workflow may explicitly certify that every rank enters the
    // force-evaluation surface via global synchronization metadata.
    boundary.pm_refresh_allowed = active_set.has_global_synchronization_metadata;
    return boundary;
  }
  boundary.kind = requested_kind;
  boundary.restart_safe = isRestartSafeBoundary(requested_kind);
  boundary.output_safe = isOutputSafeBoundary(requested_kind);
  boundary.pm_refresh_allowed = requested_kind == StepBoundaryKind::kGlobalSynchronizationPoint ||
      requested_kind == StepBoundaryKind::kPmRefreshPoint ||
      requested_kind == StepBoundaryKind::kCheckpointPoint ||
      requested_kind == StepBoundaryKind::kSnapshotPoint;
  return boundary;
}

RestartBoundaryDecision evaluateRestartBoundary(
    const IntegratorState& integrator_state,
    std::optional<std::uint64_t> scheduler_tick) {
  const PmSynchronizationPersistentState pm_sync_state =
      integrator_state.pm_sync_state.exportPersistentState();
  const bool current_local_substep =
      integrator_state.current_boundary_kind == StepBoundaryKind::kLocalActiveBinStep;
  const bool last_completed_local_substep =
      integrator_state.last_completed_boundary_kind == StepBoundaryKind::kLocalActiveBinStep;
  const bool local_substep_active =
      current_local_substep || last_completed_local_substep;
  const bool pm_refresh_legal = !pm_sync_state.refresh_commit_pending;

  RestartBoundaryDecision decision{
      .current_boundary_kind = integrator_state.current_boundary_kind,
      .last_completed_boundary_kind = integrator_state.last_completed_boundary_kind,
      .inside_kdk_step = integrator_state.inside_kdk_step,
      .last_completed_restart_safe = integrator_state.last_completed_restart_safe,
      .local_substep_active = local_substep_active,
      .pm_refresh_legal = pm_refresh_legal,
      .pm_refresh_commit_pending = pm_sync_state.refresh_commit_pending,
      .step_index = integrator_state.step_index,
      .scheduler_tick = scheduler_tick,
      .diagnostic = {},
  };
  decision.restart_safe = !decision.inside_kdk_step &&
      decision.last_completed_restart_safe &&
      !decision.local_substep_active &&
      decision.pm_refresh_legal &&
      isRestartSafeBoundary(decision.last_completed_boundary_kind);

  std::ostringstream out;
  out << "restart boundary " << (decision.restart_safe ? "safe" : "unsafe")
      << "; current_boundary_kind=" << stepBoundaryKindName(decision.current_boundary_kind)
      << "; last_completed_boundary_kind=" << stepBoundaryKindName(decision.last_completed_boundary_kind)
      << "; inside_kdk_step=" << (decision.inside_kdk_step ? "true" : "false")
      << "; last_completed_restart_safe=" << (decision.last_completed_restart_safe ? "true" : "false")
      << "; local_substep_active=" << (decision.local_substep_active ? "true" : "false")
      << "; pm_refresh_legal=" << (decision.pm_refresh_legal ? "true" : "false")
      << "; pm_refresh_commit_pending=" << (decision.pm_refresh_commit_pending ? "true" : "false")
      << "; step_index=" << decision.step_index;
  if (decision.scheduler_tick.has_value()) {
    out << "; scheduler_tick=" << *decision.scheduler_tick;
  } else {
    out << "; scheduler_tick=n/a";
  }
  if (decision.inside_kdk_step) {
    out << "; reason=half-step restart is not represented or reloadable";
  } else if (decision.local_substep_active) {
    out << "; reason=local active-bin substep restart is not represented or reloadable";
  } else if (!decision.pm_refresh_legal) {
    out << "; reason=PM refresh transition has an uncommitted long-range field update";
  } else if (!decision.last_completed_restart_safe ||
             !isRestartSafeBoundary(decision.last_completed_boundary_kind)) {
    out << "; reason=last completed integration boundary is not a restart-safe global checkpoint boundary";
  }
  decision.diagnostic = out.str();
  return decision;
}

bool canWriteRestart(
    const IntegratorState& integrator_state,
    std::optional<std::uint64_t> scheduler_tick) {
  return evaluateRestartBoundary(integrator_state, scheduler_tick).restart_safe;
}

void assertCanWriteSnapshotAtBoundary(const IntegratorState& integrator_state) {
  if (integrator_state.inside_kdk_step || !integrator_state.last_completed_restart_safe ||
      !isOutputSafeBoundary(integrator_state.last_completed_boundary_kind)) {
    throw std::runtime_error(
        "snapshot output requested at unsafe integration boundary: " +
        std::string(stepBoundaryKindName(integrator_state.last_completed_boundary_kind)));
  }
}

void assertCanWriteCheckpointAtBoundary(
    const IntegratorState& integrator_state,
    std::optional<std::uint64_t> scheduler_tick) {
  const RestartBoundaryDecision decision = evaluateRestartBoundary(integrator_state, scheduler_tick);
  if (!decision.restart_safe) {
    throw std::runtime_error("restart checkpoint requested at unsafe integration boundary: " + decision.diagnostic);
  }
}

}  // namespace cosmosim::core
