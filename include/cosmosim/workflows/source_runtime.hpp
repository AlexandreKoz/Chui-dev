#pragma once

#include <cstdint>
#include <memory>
#include <span>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::physics {
class EffectiveMultiphaseEosTable;
}  // namespace cosmosim::physics

namespace cosmosim::core {
struct MemoryReport;
}

namespace cosmosim::workflows {

struct RuntimeServices;

// Owns construction and deterministic registration order for the existing
// rung-zero source callbacks. No source model or update law is changed here.
class SourceRuntime {
 public:
  virtual ~SourceRuntime() = default;

  virtual void execute(SourceMutationStageView& view) = 0;
  // T02: smooth thermochemistry/dust/SGS on one recursive AMR level
  // subinterval. The caller is the AMR clock owner and invokes both phases on
  // every rank, including empty local spans. Discrete sources and diffusion
  // retain their separate global chronology until their dedicated contracts
  // are connected; this method must never execute them implicitly.
  virtual void executeAmrLevelSmoothSources(
      core::StepContext& parent_context,
      std::span<const std::uint32_t> level_rows,
      const hydro::HydroUpdateContext& level_update,
      core::SourceSplitPhase phase) = 0;
  virtual void setAmrLevelSourceCouplingEnabled(bool enabled) noexcept = 0;
  [[nodiscard]] virtual RuntimeTaskMemoryEstimate estimateMemory(
      const core::SimulationState&) const {
    return {0U, false, "source owner does not publish a complete memory model"};
  }
  [[nodiscard]] virtual core::MemoryReport memoryReport() const = 0;

  // G1.6 whole-step transaction: persistent source-controller state that
  // lives OUTSIDE SimulationState (model clamp/diagnostic counters, the
  // particle-ID occupancy cache, last published ledgers) is checkpointed at
  // the coordinator's step rollback point and restored when the step is
  // rejected, so a retried step is bit-identical to a step that never
  // failed. Everything authoritative (chemistry controller, populations,
  // event ordinals, composition) is in SimulationState and restored with it.
  // All three calls are rank-uniform (the coordinator calls them on every
  // rank); restore performs no collective (collective caches are invalidated
  // and rebuilt at their next collective use).
  virtual void captureRollbackCheckpoint() {}
  virtual void restoreRollbackCheckpoint() {}
  virtual void discardRollbackCheckpoint() noexcept {}
};

[[nodiscard]] std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const core::UnitSystem& units,
    std::uint32_t world_rank,
    const parallel::MpiContext& mpi_context,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table = nullptr,
    const RuntimeServices* runtime_services = nullptr);

}  // namespace cosmosim::workflows
