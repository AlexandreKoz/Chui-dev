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
}

namespace cosmosim::hydro {
struct HydroUpdateContext;
}  // namespace cosmosim::hydro

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

  // T01 deterministic Strang reference S(dt/2) H(dt) S(dt/2). The time
  // integrator dispatches this stage twice per (block) interval with
  // StepContext::source_split_phase = kStrangPreHalf, then kStrangPostHalf.
  //   * Smooth operators (thermochemistry/ISM/dust/SGS, element diffusion
  //     and SGS shear transfer) advance the physical half interval
  //     [t_begin, t_mid] or [t_mid, t_end] at that half's physical midpoint.
  //   * Discrete/stochastic events (star birth, SSP/stellar returns, winds,
  //     SN/feedback deposition, BH seeding/accretion/feedback) execute ONLY in
  //     kStrangPostHalf, once, at the absolute end-of-interval epoch with the
  //     full interval duration (never twice per coarse step).
  //   * The adaptive chemistry controller accumulates shadow observations from
  //     every half/level of one interval and publishes exactly once when the
  //     post half closes the interval. A rejected trial is rolled back through
  //     capture/restoreRollbackCheckpoint() and SimulationState restoration.
  // kUnsplitCompatibility (direct low-level dispatch) runs smooth sources over
  // the full interval, then events, then publishes; it opens and closes the
  // interval in one call.
  virtual void execute(SourceMutationStageView& view) = 0;
  // T02: smooth thermochemistry/ISM/dust/SGS on one recursive AMR level
  // subinterval [level t_begin, t_end], for the pre-hydro or post-hydro half
  // of that level's own Strang split. The AMR clock owner calls both phases on
  // every rank for every level clock, including empty local spans (the stage
  // vote is collective). Observations accumulate into the open coarse source
  // interval; nothing is published here. Under recursive AMR the level
  // halves also integrate element diffusion and SGS shear on the faces owned
  // by this level clock (same-level faces in both halves, coarse/fine faces
  // whose coarser side is this level in the pre half only; their closing half
  // is executeAmrLevelInterfaceSynchronization). Discrete events keep their
  // synchronized coarse chronology in execute(). Throws
  // std::logic_error when composition did not enable the coupling or when no
  // coarse source interval is open.
  virtual void executeAmrLevelSmoothSources(
      core::StepContext& parent_context,
      std::span<const std::uint32_t> level_rows,
      const hydro::HydroUpdateContext& level_update,
      int level,
      core::SourceSplitPhase phase) = 0;
  // Composition-root switch, set by HydroAmrRuntime::bindSourceRuntime() when
  // the configured AMR stepping is recursive subcycling. While enabled and the
  // hierarchy has production AMR coverage, execute() skips its global smooth
  // thermochemistry halves (the AMR levels own them) and verifies at interval
  // close that every level pre half had a matching post half.
  virtual void setAmrLevelSourceCouplingEnabled(bool enabled) noexcept = 0;
  // T02 (patch 51): closing half of the coarse/fine interface split of level
  // `level`, called by the recursive AMR driver after all children of the
  // level interval, its reflux and history retirement (every rank, same
  // order). Integrates the diffusion/SGS faces whose coarser side is on
  // `level` and whose finer side is deeper, for dt_level/2, with both sides
  // at the level end time. Must follow the matching level pre half.
  virtual void executeAmrLevelInterfaceSynchronization(
      core::StepContext& parent_context,
      const hydro::HydroUpdateContext& level_update,
      int level) = 0;
  [[nodiscard]] virtual RuntimeTaskMemoryEstimate estimateMemory(
      const core::SimulationState&) const {
    return {0U, false, "source owner does not publish a complete memory model"};
  }
  [[nodiscard]] virtual core::MemoryReport memoryReport() const = 0;

  // A03 bridge: the production regrid planner may ask the source owner for
  // the same thermochemistry model's cooling time (code units, +inf where no
  // cooling) on every dense gas row. Returns false when the owner has no
  // Part-2 thermochemistry; the caller then has no cooling-length criterion.
  [[nodiscard]] virtual bool fillAmrCoolingTimes(
      core::StepContext&, std::span<double>) { return false; }

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
