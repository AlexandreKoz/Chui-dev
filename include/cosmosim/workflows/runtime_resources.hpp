#pragma once

#include <cstdint>
#include <span>
#include <string_view>

namespace cosmosim::core {
struct IntegratorState;
class SimulationState;
class HierarchicalTimeBinScheduler;
}

namespace cosmosim::workflows {

class RuntimeExecutionPlan;
class TimeCoordinator;
namespace internal {
class RuntimeStageAccess;
class RuntimeTaskGrantBinder;
struct RuntimeStageResourceBundle;
}

// Resource declarations are part of the executable task contract. They live
// beside the stage views so the dispatcher can bind one task-scoped grant
// without introducing a dependency cycle through the module registry.
enum class RuntimeResourceAccessMode : std::uint8_t {
  kRead,
  kWrite,
  kReadWrite,
};

enum class RuntimeResourceKey : std::uint8_t {
  kParticlePosition,
  kParticleVelocity,
  kParticleIdentity,
  kParticleSpeciesIndex,
  kParticleGravitySource,
  kGravityAcceleration,
  kHydroConservedState,
  kHydroPrimitiveState,
  kAmrPatchState,
  kEffectiveIsmThermodynamics,
  kSourceMutationState,
  kMigrationOwnership,
  kSchedulerTruth,
  kIntegratorTruth,
  kOutputRestartState,
  kDiagnostics,
  // Pass 04 (M03-11): the authoritative oriented magnetic face flux
  // (SimulationState::magnetic_faces). Only the kHydroUpdate hydro/MHD owner
  // holds write access; output/restart reads it. Appended last so existing
  // ordinals are unchanged.
  kMagneticFaceState,
};

struct RuntimeResourceAccess {
  RuntimeResourceKey resource = RuntimeResourceKey::kParticlePosition;
  RuntimeResourceAccessMode mode = RuntimeResourceAccessMode::kRead;
};

[[nodiscard]] bool runtimeResourceAccessSatisfies(
    RuntimeResourceAccess granted,
    RuntimeResourceAccess required) noexcept;

// A transient stage view is derived from these runtime authorities. It is
// never restart truth: owners rebuild views after every invalidating event.
struct RuntimeResourceEpoch {
  std::uint64_t particle_index_generation = 0;
  std::uint64_t cell_index_generation = 0;
  std::uint64_t gas_identity_generation = 0;
  std::uint64_t scheduler_tick = 0;
  std::uint64_t step_index = 0;
};

enum class RuntimeEpochField : std::uint32_t {
  kNone = 0,
  kParticleIndex = 1U << 0U,
  kCellIndex = 1U << 1U,
  kGasIdentity = 1U << 2U,
  kSchedulerTick = 1U << 3U,
  kStepIndex = 1U << 4U,
};

[[nodiscard]] constexpr RuntimeEpochField operator|(
    RuntimeEpochField lhs,
    RuntimeEpochField rhs) noexcept {
  return static_cast<RuntimeEpochField>(
      static_cast<std::uint32_t>(lhs) | static_cast<std::uint32_t>(rhs));
}

[[nodiscard]] constexpr RuntimeEpochField operator&(
    RuntimeEpochField lhs,
    RuntimeEpochField rhs) noexcept {
  return static_cast<RuntimeEpochField>(
      static_cast<std::uint32_t>(lhs) & static_cast<std::uint32_t>(rhs));
}

class RuntimeResourceEpochSource {
 public:
  virtual ~RuntimeResourceEpochSource() = default;
  [[nodiscard]] virtual RuntimeResourceEpoch currentRuntimeEpoch() const noexcept = 0;
};

// Adapter over existing runtime authorities. It borrows them without exposing
// any of their mutation APIs through the stage-view surface.
class SimulationRuntimeEpochSource final : public RuntimeResourceEpochSource {
 public:
  SimulationRuntimeEpochSource(
      const core::SimulationState& state,
      const core::HierarchicalTimeBinScheduler& scheduler,
      const core::IntegratorState& integrator_state) noexcept;

  [[nodiscard]] RuntimeResourceEpoch currentRuntimeEpoch() const noexcept override;

 private:
  const core::SimulationState* m_state = nullptr;
  const core::HierarchicalTimeBinScheduler* m_scheduler = nullptr;
  const core::IntegratorState* m_integrator_state = nullptr;
};

// Capability carried by every production stage view. The capability exposes
// freshness validation only; it cannot mutate an epoch or recover the state
// object from which the view was derived.
class RuntimeResourceLease {
 public:
  RuntimeResourceLease(
      const RuntimeResourceEpochSource& source,
      RuntimeEpochField guarded_fields) noexcept;

  [[nodiscard]] const RuntimeResourceEpoch& capturedEpoch() const noexcept;
  [[nodiscard]] RuntimeEpochField guardedFields() const noexcept;
  [[nodiscard]] bool isFresh() const noexcept;
  void requireFresh(std::string_view view_name) const;

 private:
  const RuntimeResourceEpochSource* m_source = nullptr;
  RuntimeResourceEpoch m_captured_epoch{};
  RuntimeEpochField m_guarded_fields = RuntimeEpochField::kNone;
};

// Public stage views deliberately expose freshness and their task-scoped
// declaration only. They do not carry StepContext, SimulationState, or a public
// owner escape. Built-in owners use the source-private RuntimeStageAccess bridge,
// which verifies this task's descriptor grant before returning internal state.
#define COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW(Type)                         \
  class Type {                                                            \
   public:                                                                \
    explicit Type(RuntimeResourceLease lease) noexcept;                   \
    void requireFresh() const;                                            \
    [[nodiscard]] std::span<const RuntimeResourceAccess>                  \
    declaredResources() const noexcept;                                   \
                                                                          \
   private:                                                               \
    friend class RuntimeExecutionPlan;                                    \
    friend class TimeCoordinator;                                         \
    friend class internal::RuntimeStageAccess;                            \
    friend class internal::RuntimeTaskGrantBinder;                       \
    Type(                                                                 \
        RuntimeResourceLease lease,                                       \
        internal::RuntimeStageResourceBundle& resources) noexcept;        \
    void bindDeclaredResources(                                           \
        std::span<const RuntimeResourceAccess> resources) noexcept;       \
    void clearDeclaredResources() noexcept;                               \
    void requireDeclaredResources(                                        \
        std::span<const RuntimeResourceAccess> required,                  \
        std::string_view caller) const;                                   \
                                                                          \
    RuntimeResourceLease m_lease;                                         \
    internal::RuntimeStageResourceBundle* m_resources = nullptr;          \
    std::span<const RuntimeResourceAccess> m_declared_resources;          \
  };

COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW(DriftParticleStageView)
COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW(GravityStageView)
COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW(HydroAmrStageView)
COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW(SourceMutationStageView)
COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW(AnalysisStageView)
COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW(OutputRestartStageView)

#undef COSMOSIM_DECLARE_RUNTIME_STAGE_VIEW

class MigrationOwnershipView {
 public:
  explicit MigrationOwnershipView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  RuntimeResourceLease m_lease;
};

class TimeCriteriaStageView {
 public:
  explicit TimeCriteriaStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  RuntimeResourceLease m_lease;
};

}  // namespace cosmosim::workflows
