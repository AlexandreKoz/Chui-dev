#include "cosmosim/workflows/runtime_module_registry.hpp"

#include <algorithm>
#include <exception>
#include "cosmosim/workflows/runtime_services.hpp"
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace cosmosim::workflows {

namespace internal {

class RuntimeTaskGrantBinder final {
 public:
  template <class View>
  static void bind(
      View& view,
      std::span<const RuntimeResourceAccess> resources) noexcept {
    view.bindDeclaredResources(resources);
  }

  template <class View>
  static void clear(View& view) noexcept {
    view.clearDeclaredResources();
  }
};

}  // namespace internal

namespace {

[[nodiscard]] bool viewSupportsStage(
    RuntimeStageViewKind view_kind,
    core::IntegrationStage stage) noexcept {
  switch (view_kind) {
    case RuntimeStageViewKind::kStageAudit:
      return true;
    case RuntimeStageViewKind::kDriftParticles:
      return stage == core::IntegrationStage::kDrift;
    case RuntimeStageViewKind::kGravity:
      return stage == core::IntegrationStage::kGravityKickPre ||
             stage == core::IntegrationStage::kForceRefresh ||
             stage == core::IntegrationStage::kGravityKickPost;
    case RuntimeStageViewKind::kHydroAmr:
      return stage == core::IntegrationStage::kHydroUpdate;
    case RuntimeStageViewKind::kSourceMutation:
      return stage == core::IntegrationStage::kSourceTerms;
    case RuntimeStageViewKind::kAnalysis:
      return stage == core::IntegrationStage::kAnalysisHooks;
    case RuntimeStageViewKind::kOutputRestart:
      return stage == core::IntegrationStage::kOutputCheck;
  }
  return false;
}

[[nodiscard]] std::optional<RuntimeResourceAccessMode> allowedModeForView(
    RuntimeStageViewKind view_kind,
    RuntimeResourceKey resource) noexcept {
  using Mode = RuntimeResourceAccessMode;
  using Key = RuntimeResourceKey;
  switch (view_kind) {
    case RuntimeStageViewKind::kStageAudit:
      return resource == Key::kDiagnostics
          ? std::optional<Mode>(Mode::kWrite)
          : std::nullopt;
    case RuntimeStageViewKind::kDriftParticles:
      switch (resource) {
        case Key::kParticlePosition: return Mode::kReadWrite;
        case Key::kParticleVelocity: return Mode::kReadWrite;
        case Key::kMigrationOwnership: return Mode::kRead;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kGravity:
      switch (resource) {
        case Key::kParticlePosition: return Mode::kRead;
        case Key::kParticleVelocity: return Mode::kReadWrite;
        case Key::kParticleGravitySource: return Mode::kRead;
        case Key::kHydroConservedState: return Mode::kReadWrite;
        case Key::kHydroPrimitiveState: return Mode::kReadWrite;
        case Key::kMigrationOwnership: return Mode::kRead;
        case Key::kSchedulerTruth: return Mode::kRead;
        case Key::kGravityAcceleration: return Mode::kWrite;
        case Key::kIntegratorTruth: return Mode::kReadWrite;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kHydroAmr:
      switch (resource) {
        case Key::kParticleVelocity: return Mode::kReadWrite;
        case Key::kSourceMutationState: return Mode::kReadWrite;
        case Key::kHydroConservedState: return Mode::kReadWrite;
        case Key::kHydroPrimitiveState: return Mode::kReadWrite;
        case Key::kAmrPatchState: return Mode::kReadWrite;
        case Key::kEffectiveIsmThermodynamics: return Mode::kReadWrite;
        case Key::kGravityAcceleration: return Mode::kRead;
        case Key::kMigrationOwnership: return Mode::kRead;
        case Key::kIntegratorTruth: return Mode::kRead;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kSourceMutation:
      switch (resource) {
        case Key::kSourceMutationState: return Mode::kReadWrite;
        case Key::kParticlePosition: return Mode::kReadWrite;
        case Key::kParticleVelocity: return Mode::kReadWrite;
        case Key::kParticleIdentity: return Mode::kReadWrite;
        case Key::kParticleSpeciesIndex: return Mode::kReadWrite;
        case Key::kHydroConservedState: return Mode::kReadWrite;
        case Key::kHydroPrimitiveState: return Mode::kRead;
        case Key::kAmrPatchState: return Mode::kRead;
        case Key::kEffectiveIsmThermodynamics: return Mode::kRead;
        case Key::kMigrationOwnership: return Mode::kReadWrite;
        case Key::kIntegratorTruth: return Mode::kRead;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kAnalysis:
      switch (resource) {
        case Key::kParticlePosition:
        case Key::kParticleVelocity:
        case Key::kParticleGravitySource:
        case Key::kHydroPrimitiveState:
        case Key::kSourceMutationState:
        case Key::kMigrationOwnership:
        case Key::kIntegratorTruth:
          return Mode::kRead;
        case Key::kDiagnostics:
          return Mode::kWrite;
        default:
          return std::nullopt;
      }
    case RuntimeStageViewKind::kOutputRestart:
      switch (resource) {
        case Key::kParticlePosition:
        case Key::kParticleVelocity:
        case Key::kHydroConservedState:
        case Key::kSourceMutationState:
        case Key::kMigrationOwnership:
        case Key::kSchedulerTruth:
        case Key::kIntegratorTruth:
          return Mode::kRead;
        case Key::kOutputRestartState:
          return Mode::kReadWrite;
        case Key::kDiagnostics:
          return Mode::kWrite;
        default:
          return std::nullopt;
      }
  }
  return std::nullopt;
}

[[nodiscard]] bool viewAllowsAccess(
    RuntimeStageViewKind view_kind,
    RuntimeResourceAccess requested) noexcept {
  const auto allowed_mode = allowedModeForView(view_kind, requested.resource);
  return allowed_mode.has_value() && runtimeResourceAccessSatisfies(
      RuntimeResourceAccess{requested.resource, *allowed_mode}, requested);
}

[[nodiscard]] bool accessModesConflict(
    RuntimeResourceAccessMode lhs,
    RuntimeResourceAccessMode rhs) noexcept {
  return lhs != RuntimeResourceAccessMode::kRead || rhs != RuntimeResourceAccessMode::kRead;
}

[[nodiscard]] bool hasResourceConflict(
    const RuntimeTaskDeclaration& lhs,
    const RuntimeTaskDeclaration& rhs) noexcept {
  for (const RuntimeResourceAccess lhs_access : lhs.resources) {
    for (const RuntimeResourceAccess rhs_access : rhs.resources) {
      if (lhs_access.resource == rhs_access.resource &&
          accessModesConflict(lhs_access.mode, rhs_access.mode)) {
        return true;
      }
    }
  }
  return false;
}

void validateDescriptor(const RuntimeModuleDescriptor& descriptor) {
  if (descriptor.module_id.empty()) {
    throw std::invalid_argument("runtime module descriptor has an empty module_id");
  }
  if (!descriptor.factory) {
    throw std::invalid_argument(
        "runtime module '" + descriptor.module_id + "' has no factory");
  }
  std::set<std::string> task_ids;
  for (const RuntimeTaskDeclaration& task : descriptor.stage_tasks) {
    if (task.task_id.empty()) {
      throw std::invalid_argument(
          "runtime module '" + descriptor.module_id + "' has an empty task_id");
    }
    if (!task_ids.insert(task.task_id).second) {
      throw std::invalid_argument(
          "runtime module '" + descriptor.module_id +
          "' declares duplicate task '" + task.task_id + "'");
    }
    if (!viewSupportsStage(task.view_kind, task.stage)) {
      throw std::invalid_argument(
          "runtime task '" + task.task_id +
          "' declares a view that cannot execute its stage");
    }
    if (task.scheduling.memory_ownership == RuntimeTaskMemoryOwnership::kDispatcherOwned &&
        !task.scheduling.peak_is_known) {
      throw std::invalid_argument(
          "dispatcher-owned runtime task requires a complete incremental peak contract");
    }
    if (task.scheduling.estimated_peak_bytes != 0U &&
        task.scheduling.memory_class == core::MemoryClass::kExternalRuntime) {
      throw std::invalid_argument(
          "runtime task '" + task.task_id +
          "' cannot reserve opaque external-runtime memory as a task peak");
    }
    for (const std::string& dependency : task.dependencies) {
      if (dependency.empty()) {
        throw std::invalid_argument(
            "runtime task '" + task.task_id + "' declares an empty dependency");
      }
    }
    if (task.resources.empty()) {
      throw std::invalid_argument(
          "runtime task '" + task.task_id + "' declares no resource access");
    }
    std::set<RuntimeResourceKey> resources;
    for (const RuntimeResourceAccess access : task.resources) {
      if (!resources.insert(access.resource).second) {
        throw std::invalid_argument(
            "runtime task '" + task.task_id +
            "' declares duplicate access for resource '" +
            std::string(runtimeResourceKeyName(access.resource)) + "'");
      }
      if (!viewAllowsAccess(task.view_kind, access)) {
        throw std::invalid_argument(
            "runtime task '" + task.task_id +
            "' declares access outside its typed view for resource '" +
            std::string(runtimeResourceKeyName(access.resource)) + "'");
      }
    }
  }
}

[[nodiscard]] std::vector<std::size_t> resolveModuleOrder(
    const std::vector<RuntimeModuleDescriptor>& descriptors) {
  std::unordered_map<std::string, std::size_t> by_id;
  for (std::size_t index = 0; index < descriptors.size(); ++index) {
    if (!by_id.emplace(descriptors[index].module_id, index).second) {
      throw std::invalid_argument(
          "duplicate runtime module id '" + descriptors[index].module_id + "'");
    }
  }

  for (const RuntimeModuleDescriptor& descriptor : descriptors) {
    for (const std::string& incompatible : descriptor.incompatibilities) {
      if (by_id.contains(incompatible)) {
        throw std::invalid_argument(
            "runtime module '" + descriptor.module_id +
            "' is incompatible with registered module '" + incompatible + "'");
      }
    }
    for (const std::string& prerequisite : descriptor.prerequisites) {
      if (!by_id.contains(prerequisite)) {
        throw std::invalid_argument(
            "runtime module '" + descriptor.module_id +
            "' requires missing module '" + prerequisite + "'");
      }
      if (prerequisite == descriptor.module_id) {
        throw std::invalid_argument(
            "runtime module '" + descriptor.module_id + "' requires itself");
      }
    }
  }

  std::vector<std::size_t> ordered;
  std::set<std::size_t> emitted;
  while (ordered.size() != descriptors.size()) {
    std::vector<std::size_t> ready;
    for (std::size_t index = 0; index < descriptors.size(); ++index) {
      if (emitted.contains(index)) {
        continue;
      }
      const bool prerequisites_ready = std::all_of(
          descriptors[index].prerequisites.begin(),
          descriptors[index].prerequisites.end(),
          [&](const std::string& prerequisite) {
            return emitted.contains(by_id.at(prerequisite));
          });
      if (prerequisites_ready) {
        ready.push_back(index);
      }
    }
    if (ready.empty()) {
      throw std::invalid_argument("runtime module prerequisite graph contains a cycle");
    }
    std::sort(
        ready.begin(),
        ready.end(),
        [&](std::size_t lhs, std::size_t rhs) {
          return std::tuple{
                     descriptors[lhs].construction_ordinal,
                     descriptors[lhs].module_id} <
                 std::tuple{
                     descriptors[rhs].construction_ordinal,
                     descriptors[rhs].module_id};
        });
    for (const std::size_t index : ready) {
      emitted.insert(index);
      ordered.push_back(index);
    }
  }
  return ordered;
}

}  // namespace

std::uint64_t incrementalMemoryBeyondRetained(
    std::uint64_t physical_peak_bytes, std::uint64_t retained_bytes) noexcept {
  return physical_peak_bytes > retained_bytes
      ? physical_peak_bytes - retained_bytes : 0U;
}

RuntimeTaskMemoryEstimate evaluateRuntimeTaskMemoryEstimate(
    const RuntimeExecutionPlan::PlannedTask& task) {
  if (task.estimate_memory) {
    return task.estimate_memory();
  }
  const std::uint64_t bytes = task.estimate_incremental_bytes
      ? task.estimate_incremental_bytes()
      : task.declaration.scheduling.estimated_peak_bytes;
  return {bytes, task.declaration.scheduling.peak_is_known,
          task.declaration.scheduling.peak_is_known
              ? std::string_view{} : std::string_view{"owner peak is not complete"}};
}

namespace {

template <class View>
class TaskGrantScope final {
 public:
  TaskGrantScope(
      View& view,
      std::span<const RuntimeResourceAccess> resources) noexcept
      : m_view(view) {
    internal::RuntimeTaskGrantBinder::bind(m_view, resources);
  }

  ~TaskGrantScope() {
    internal::RuntimeTaskGrantBinder::clear(m_view);
  }

  TaskGrantScope(const TaskGrantScope&) = delete;
  TaskGrantScope& operator=(const TaskGrantScope&) = delete;

 private:
  View& m_view;
};

template <class View, class Task>
void executeTypedStage(
    std::span<const RuntimeExecutionPlan::PlannedTask> tasks,
    core::IntegrationStage stage,
    RuntimeStageViewKind expected_view,
    View& view,
    const RuntimeServices* services) {
  view.requireFresh();
  for (const RuntimeExecutionPlan::PlannedTask& task : tasks) {
    if (task.declaration.stage != stage) {
      continue;
    }
    if (task.declaration.view_kind != expected_view) {
      continue;
    }
    if (!std::holds_alternative<Task>(task.task)) {
      throw std::logic_error(
          "frozen runtime task '" + task.declaration.task_id +
          "' has an inconsistent typed stage view");
    }
    core::MemoryReservation task_reservation;
    if (services != nullptr && services->memory_governor != nullptr) {
      const auto& profile = task.declaration.scheduling;
      std::exception_ptr admission_failure;
      try {
        const RuntimeTaskMemoryEstimate estimate =
            evaluateRuntimeTaskMemoryEstimate(task);
        if (profile.memory_ownership == RuntimeTaskMemoryOwnership::kDispatcherOwned) {
          if (!estimate.complete) {
            throw std::logic_error(
                "dispatcher-owned task requires a complete memory peak contract");
          }
          task_reservation = services->memory_governor->reserve(
              profile.memory_class, estimate.incremental_bytes,
              task.declaration.task_id);
          task_reservation.commit();
        } else if (estimate.complete &&
                   (task.estimate_memory || task.estimate_incremental_bytes)) {
          // Preflight only: the owner already reserves its physical workspace.
          // Holding this duplicate reservation through execution would charge
          // the same byte range twice and spuriously reject feasible work.
          auto preflight = services->memory_governor->reserve(
              profile.memory_class, estimate.incremental_bytes,
              task.declaration.task_id);
          preflight.release();
        }
        // Incomplete owner models remain allocation-time governed. A partial
        // estimate must never be promoted to a certificate of concurrency.
      } catch (...) {
        admission_failure = std::current_exception();
      }
      FailureCoordinator(*services).rethrowCollectiveFailure(
          admission_failure, "runtime task memory admission");
    }
    TaskGrantScope<View> grant_scope(view, task.declaration.resources);
    std::get<Task>(task.task)(view);
  }
}

}  // namespace

RuntimeStageViewKind taskViewKind(const RuntimeStageTaskFunction& task) noexcept {
  return std::visit(
      [](const auto& typed_task) {
        using Task = std::decay_t<decltype(typed_task)>;
        if constexpr (std::is_same_v<Task, StageAuditTask>) {
          return RuntimeStageViewKind::kStageAudit;
        } else if constexpr (std::is_same_v<Task, DriftStageTask>) {
          return RuntimeStageViewKind::kDriftParticles;
        } else if constexpr (std::is_same_v<Task, GravityStageTask>) {
          return RuntimeStageViewKind::kGravity;
        } else if constexpr (std::is_same_v<Task, HydroAmrStageTask>) {
          return RuntimeStageViewKind::kHydroAmr;
        } else if constexpr (std::is_same_v<Task, SourceMutationStageTask>) {
          return RuntimeStageViewKind::kSourceMutation;
        } else if constexpr (std::is_same_v<Task, AnalysisStageTask>) {
          return RuntimeStageViewKind::kAnalysis;
        } else {
          return RuntimeStageViewKind::kOutputRestart;
        }
      },
      task);
}

std::string_view runtimeResourceKeyName(RuntimeResourceKey resource) noexcept {
  switch (resource) {
    case RuntimeResourceKey::kParticlePosition: return "particle_position";
    case RuntimeResourceKey::kParticleVelocity: return "particle_velocity";
    case RuntimeResourceKey::kParticleIdentity: return "particle_identity";
    case RuntimeResourceKey::kParticleSpeciesIndex: return "particle_species_index";
    case RuntimeResourceKey::kParticleGravitySource: return "particle_gravity_source";
    case RuntimeResourceKey::kGravityAcceleration: return "gravity_acceleration";
    case RuntimeResourceKey::kHydroConservedState: return "hydro_conserved_state";
    case RuntimeResourceKey::kHydroPrimitiveState: return "hydro_primitive_state";
    case RuntimeResourceKey::kAmrPatchState: return "amr_patch_state";
    case RuntimeResourceKey::kEffectiveIsmThermodynamics: return "effective_ism_thermodynamics";
    case RuntimeResourceKey::kSourceMutationState: return "source_mutation_state";
    case RuntimeResourceKey::kMigrationOwnership: return "migration_ownership";
    case RuntimeResourceKey::kSchedulerTruth: return "scheduler_truth";
    case RuntimeResourceKey::kIntegratorTruth: return "integrator_truth";
    case RuntimeResourceKey::kOutputRestartState: return "output_restart_state";
    case RuntimeResourceKey::kDiagnostics: return "diagnostics";
  }
  return "unknown";
}

bool runtimeTasksMayOverlap(
    const RuntimeTaskDeclaration& lhs,
    const RuntimeTaskDeclaration& rhs,
    const core::MemoryGovernorSnapshot& memory_snapshot) noexcept {
  if (!lhs.scheduling.peak_is_known || !rhs.scheduling.peak_is_known ||
      !lhs.dependencies.empty() || !rhs.dependencies.empty() ||
      lhs.scheduling.lifetime_boundary != RuntimeTaskLifetimeBoundary::kTaskEnd ||
      rhs.scheduling.lifetime_boundary != RuntimeTaskLifetimeBoundary::kTaskEnd) {
    return false;
  }
  if (hasResourceConflict(lhs, rhs)) {
    return false;
  }
  const auto high = RuntimeTaskPressureClass::kHigh;
  if ((lhs.scheduling.memory_bandwidth_pressure == high &&
       rhs.scheduling.memory_bandwidth_pressure == high) ||
      (lhs.scheduling.communication_pressure == high &&
       rhs.scheduling.communication_pressure == high) ||
      (lhs.scheduling.compute_pressure == high &&
       rhs.scheduling.compute_pressure == high)) {
    return false;
  }
  const std::uint64_t lhs_peak = lhs.scheduling.estimated_peak_bytes;
  const std::uint64_t rhs_peak = rhs.scheduling.estimated_peak_bytes;
  {
    if (rhs_peak > std::numeric_limits<std::uint64_t>::max() - lhs_peak) {
      return false;
    }
    const std::uint64_t simultaneous_peak = lhs_peak + rhs_peak;
    if (memory_snapshot.headroom_bytes != std::numeric_limits<std::uint64_t>::max() &&
        simultaneous_peak > memory_snapshot.headroom_bytes) {
      return false;
    }
  }
  return true;
}

RuntimeTaskMemoryEstimate RuntimeExecutionPlan::taskMemoryEstimate(
    std::string_view task_id) const {
  const std::size_t separator = task_id.find("::");
  if (separator == std::string_view::npos) {
    throw std::invalid_argument("runtime task memory estimate requires module::task_id");
  }
  for (const PlannedTask& task : m_tasks) {
    if (task.module_id == task_id.substr(0U, separator) &&
        task.declaration.task_id == task_id.substr(separator + 2U)) {
      return evaluateRuntimeTaskMemoryEstimate(task);
    }
  }
  throw std::out_of_range("runtime task memory estimate requested for unknown task");
}

std::size_t RuntimeExecutionPlan::moduleCount() const noexcept {
  return m_module_instances.size();
}

std::size_t RuntimeExecutionPlan::taskCount() const noexcept {
  return m_tasks.size();
}

std::span<const std::string> RuntimeExecutionPlan::orderedModuleIds() const noexcept {
  return m_ordered_module_ids;
}

void RuntimeExecutionPlan::executeAuditStage(
    core::IntegrationStage stage,
    AnalysisStageView& view) const {
  executeTypedStage<AnalysisStageView, StageAuditTask>(
      m_tasks, stage, RuntimeStageViewKind::kStageAudit, view, m_services);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    DriftParticleStageView& view) const {
  executeTypedStage<DriftParticleStageView, DriftStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kDriftParticles, view, m_services);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    GravityStageView& view) const {
  executeTypedStage<GravityStageView, GravityStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kGravity, view, m_services);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    HydroAmrStageView& view) const {
  executeTypedStage<HydroAmrStageView, HydroAmrStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kHydroAmr, view, m_services);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    SourceMutationStageView& view) const {
  executeTypedStage<SourceMutationStageView, SourceMutationStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kSourceMutation, view, m_services);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    AnalysisStageView& view) const {
  executeTypedStage<AnalysisStageView, AnalysisStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kAnalysis, view, m_services);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    OutputRestartStageView& view) const {
  executeTypedStage<OutputRestartStageView, OutputRestartStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kOutputRestart, view, m_services);
}

void RuntimeModuleRegistry::registerModule(RuntimeModuleDescriptor descriptor) {
  if (m_frozen) {
    throw std::logic_error("runtime module registry is already frozen");
  }
  validateDescriptor(descriptor);
  if (std::any_of(
          m_descriptors.begin(),
          m_descriptors.end(),
          [&](const RuntimeModuleDescriptor& existing) {
            return existing.module_id == descriptor.module_id;
          })) {
    throw std::invalid_argument(
        "duplicate runtime module id '" + descriptor.module_id + "'");
  }
  m_descriptors.push_back(std::move(descriptor));
}

bool RuntimeModuleRegistry::frozen() const noexcept { return m_frozen; }

std::size_t RuntimeModuleRegistry::descriptorCount() const noexcept {
  return m_descriptors.size();
}

std::span<const RuntimeModuleDescriptor> RuntimeModuleRegistry::descriptors() const noexcept {
  return m_descriptors;
}

RuntimeExecutionPlan RuntimeModuleRegistry::freezeAndInstantiate(
    const RuntimeModuleFactoryContext& context) {
  if (m_frozen) {
    throw std::logic_error("runtime module registry may be instantiated only once");
  }
  m_frozen = true;

  const std::vector<std::size_t> order = resolveModuleOrder(m_descriptors);
  RuntimeExecutionPlan plan;
  plan.m_services = &context.services;
  for (const std::size_t descriptor_index : order) {
    const RuntimeModuleDescriptor& descriptor = m_descriptors[descriptor_index];
    plan.m_ordered_module_ids.push_back(descriptor.module_id);
    RuntimeModuleInstance instance = descriptor.factory(context);

    std::unordered_map<std::string, RuntimeStageTaskContribution*> contribution_by_id;
    for (RuntimeStageTaskContribution& contribution : instance.stage_tasks) {
      if (contribution.task_id.empty() ||
          !contribution_by_id.emplace(contribution.task_id, &contribution).second) {
        throw std::invalid_argument(
            "runtime module factory '" + descriptor.module_id +
            "' returned an empty or duplicate task id");
      }
    }
    if (contribution_by_id.size() != descriptor.stage_tasks.size()) {
      throw std::invalid_argument(
          "runtime module factory '" + descriptor.module_id +
          "' did not return exactly its declared stage tasks");
    }
    for (const RuntimeTaskDeclaration& declaration : descriptor.stage_tasks) {
      const auto contribution_it = contribution_by_id.find(declaration.task_id);
      if (contribution_it == contribution_by_id.end()) {
        throw std::invalid_argument(
            "runtime module factory '" + descriptor.module_id +
            "' omitted declared task '" + declaration.task_id + "'");
      }
      if (taskViewKind(contribution_it->second->task) != declaration.view_kind) {
        throw std::invalid_argument(
            "runtime module factory task '" + declaration.task_id +
            "' does not match its declared typed view");
      }
      plan.m_tasks.push_back(RuntimeExecutionPlan::PlannedTask{
          .module_id = descriptor.module_id,
          .declaration = declaration,
          .task = std::move(contribution_it->second->task),
          .estimate_incremental_bytes = std::move(contribution_it->second->estimate_incremental_bytes),
          .estimate_memory = std::move(contribution_it->second->estimate_memory),
      });
    }
    plan.m_module_instances.push_back(std::move(instance));
  }

  std::stable_sort(
      plan.m_tasks.begin(),
      plan.m_tasks.end(),
      [](const RuntimeExecutionPlan::PlannedTask& lhs,
         const RuntimeExecutionPlan::PlannedTask& rhs) {
        return std::tuple{
                   static_cast<std::uint8_t>(lhs.declaration.stage),
                   lhs.declaration.ordinal,
                   lhs.module_id,
                   lhs.declaration.task_id} <
               std::tuple{
                   static_cast<std::uint8_t>(rhs.declaration.stage),
                   rhs.declaration.ordinal,
                   rhs.module_id,
                   rhs.declaration.task_id};
      });

  std::unordered_map<std::string, std::size_t> task_position_by_key;
  task_position_by_key.reserve(plan.m_tasks.size());
  for (std::size_t index = 0; index < plan.m_tasks.size(); ++index) {
    const std::string key = plan.m_tasks[index].module_id + "::" +
        plan.m_tasks[index].declaration.task_id;
    if (!task_position_by_key.emplace(key, index).second) {
      throw std::invalid_argument("duplicate fully-qualified runtime task key '" + key + "'");
    }
  }
  for (std::size_t index = 0; index < plan.m_tasks.size(); ++index) {
    const RuntimeExecutionPlan::PlannedTask& task = plan.m_tasks[index];
    for (const std::string& dependency : task.declaration.dependencies) {
      const auto dependency_it = task_position_by_key.find(dependency);
      if (dependency_it == task_position_by_key.end()) {
        throw std::invalid_argument(
            "runtime task '" + task.module_id + "::" + task.declaration.task_id +
            "' requires missing task dependency '" + dependency + "'");
      }
      if (dependency_it->second >= index) {
        throw std::invalid_argument(
            "runtime task dependency graph contradicts deterministic stage/ordinal order: '" +
            dependency + "' must precede '" + task.module_id + "::" +
            task.declaration.task_id + "'");
      }
    }
  }
  return plan;
}

}  // namespace cosmosim::workflows
