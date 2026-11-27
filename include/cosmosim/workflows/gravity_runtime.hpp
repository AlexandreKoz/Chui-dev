#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <vector>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::workflows {

struct RuntimeServices;

// Gravity-owned startup normalization and diagnostic policy description.
// These are called by the composition lifecycle before the first gravity task.
void initializeGravityState(
    core::SimulationState& state,
    const core::SimulationConfig& config);
[[nodiscard]] std::string describeGravitySofteningPolicy(
    const core::SimulationConfig& config,
    const core::SimulationState& state);

// Row-aligned synchronized hierarchical force endpoint captured at a
// multirate PM cycle opening (transactional rollback only; never persisted).
// Rows follow SimulationState row order at capture. Unlike the ordinary
// restart force cache (io::GravityForceCachePersistentState, remapped by
// stable particle ID because a restart may reorder particles), this endpoint
// is restored by row with no ID lookup and no allocation; its validity is
// tied to the frozen layout of the cycle: the particle count, the
// SimulationState particle-index generation (advanced by every reorder,
// resize, compaction or migration) and the caller's identity/layout
// certificate (fingerprint of the ID lane and every frozen lane) must all
// match at restore, otherwise the restore is refused before any write.
struct TransactionalForceEndpoint {
  std::vector<double> total_x_comoving;
  std::vector<double> total_y_comoving;
  std::vector<double> total_z_comoving;
  std::vector<double> tree_x_comoving;
  std::vector<double> tree_y_comoving;
  std::vector<double> tree_z_comoving;
  std::vector<double> pm_x_comoving;
  std::vector<double> pm_y_comoving;
  std::vector<double> pm_z_comoving;
  std::vector<double> nearest_pair_r2_comoving;
  std::uint64_t particle_count = 0U;
  std::uint64_t particle_index_generation = 0U;
  std::uint64_t layout_certificate = 0U;
  std::uint64_t pm_field_version = 0U;
  bool valid = false;

  static constexpr std::size_t k_lane_count = 10U;
  // Bytes of lane capacity currently retained (what a governor reservation
  // for this object must cover).
  [[nodiscard]] std::uint64_t retainedCapacityBytes() const noexcept;
  // Smallest retained capacity over all lanes (capture of N rows allocates
  // nothing when this is >= N).
  [[nodiscard]] std::size_t minimumLaneCapacity() const noexcept;
  // Reserves every lane to `rows` (call only after the governor admitted the
  // bytes; strong guarantee per lane: a failure leaves earlier lanes reserved).
  void reserveRows(std::size_t rows);
  // Frees all lane capacity and invalidates the endpoint.
  void release() noexcept;
};

// Narrow read surface consumed by hydro/source orchestration. Dense-row
// acceleration mirrors are invalidated by particle migration through
// commitParticleDecompositionChange().
class GravityAccelerationProvider {
 public:
  virtual ~GravityAccelerationProvider() = default;

  [[nodiscard]] virtual std::span<const double> cellAccelX() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> cellAccelY() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> cellAccelZ() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> particleAccelX() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> particleAccelY() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> particleAccelZ() const noexcept = 0;
};

// Persistence-facing gravity surface. Output/restart depends on this contract,
// never on the concrete solver implementation.
class GravityRestartStateProvider {
 public:
  virtual ~GravityRestartStateProvider() = default;

  [[nodiscard]] virtual std::uint64_t decompositionEpoch() const noexcept = 0;
  [[nodiscard]] virtual const gravity::PmGridShape& pmGridShape() const noexcept = 0;
  [[nodiscard]] virtual const parallel::DistributedExecutionTopology&
  runtimeTopology() const noexcept = 0;
  [[nodiscard]] virtual io::GravityForceCachePersistentView
  restartForceCacheView(const core::SimulationState& state) const = 0;
};

// Owner interface for rung-zero TreePM/PM cadence and force-cache lifecycle.
// Stage dispatch enters only through the typed GravityStageView contract.
class GravityRuntime : public GravityAccelerationProvider,
                       public GravityRestartStateProvider {
 public:
  ~GravityRuntime() override = default;

  virtual void execute(GravityStageView& view) = 0;
  // Current retained-capacity-aware phase model. The owner continues to hold
  // its physical reservation; the dispatcher only performs a preflight.
  // At the serial task boundary this may evict derived optional caches to make
  // headroom for mandatory work; it never changes source/force/scheduler truth.
  [[nodiscard]] virtual RuntimeTaskMemoryEstimate estimateMemory(
      const core::SimulationState&) const {
    return {0U, false, "gravity owner does not publish a complete memory model"};
  }

  [[nodiscard]] virtual std::size_t pmGridSize() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t longRangeRefreshCount() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t longRangeReuseCount() const noexcept = 0;
  [[nodiscard]] virtual int pmCadenceSteps() const noexcept = 0;
  [[nodiscard]] virtual double gravitationalConstantCode() const noexcept = 0;
  [[nodiscard]] virtual std::span<const ReferenceWorkflowReport::TreePmCadenceRecord>
  cadenceRecords() const noexcept = 0;
  // Multirate outer-PM cycle transaction (TreePmPmReusePolicy::kMultirateOuterKick).
  // Every rank calls these in the same order (the time coordinator decides
  // from globally reduced criteria).
  // * beginCycleTransaction: opens a bounded diagnostics journal. Accepted-
  //   trajectory diagnostics produced until commit/rollback (PM refresh/reuse
  //   counts, TreePM cadence records) are held there, so a rejected attempt
  //   can never evict accepted history from the bounded recent window.
  // * commitCycleTransaction: publishes the journal into accepted history.
  // * rollbackCycleTransaction: moves the journal into the rejected-attempt
  //   (execution-cost) totals, drops every source-history-dependent derived
  //   gravity state (tree build/refit candidates and topology maps, particle
  //   ghost layer, borrowed-source certificates, per-row force validity) and
  //   restores the cycle-opening endpoint force ROW-ALIGNED (no ID remap,
  //   no allocation; see TransactionalForceEndpoint) for the restored state.
  //   All restore preconditions are validated before the first write; on
  //   refusal nothing in the runtime changes and the transaction stays open
  //   (the caller's collective failure handling ends the run on every rank).
  //   Content-keyed caches with provable validity are retained (see
  //   TreePmCoordinator::invalidateSourceHistoryForRollback).
  // * captureTransactionalForceEndpoint: row-aligned copy of the complete
  //   synchronized endpoint at the current positions into caller-owned
  //   (governed, capacity-reused) storage; refuses an incomplete endpoint.
  virtual void beginCycleTransaction() = 0;
  virtual void commitCycleTransaction() noexcept = 0;
  virtual void captureTransactionalForceEndpoint(
      const core::SimulationState& state, std::uint64_t layout_certificate,
      TransactionalForceEndpoint& endpoint) const = 0;
  // Read-only check of every rollbackCycleTransaction restore precondition;
  // lets the caller refuse a rollback before mutating any authoritative state.
  virtual void validateTransactionalForceEndpoint(
      const TransactionalForceEndpoint& opening_endpoint,
      const core::SimulationState& restored_state,
      std::uint64_t layout_certificate) const = 0;
  virtual void rollbackCycleTransaction(
      const TransactionalForceEndpoint& opening_endpoint,
      const core::SimulationState& restored_state,
      std::uint64_t layout_certificate) = 0;
  [[nodiscard]] virtual bool cycleTransactionOpen() const noexcept = 0;
  // Pass 04 (M04-01): after a rollback whose caller restored only part of the
  // state (the multirate PM-cycle rollback restores particle lanes, not the
  // module sidecar registry), re-publish runtime-owned module sidecars from
  // the journal-restored live owners so sidecar and live state agree. The
  // global-step rollback restores the whole SimulationState copy and needs no
  // call. Default: no runtime-owned module sidecar.
  virtual void synchronizeModuleSidecarsAfterRollback(core::SimulationState& state) {
    static_cast<void>(state);
  }
  // Execution-cost totals of rejected attempts (not part of the trajectory).
  [[nodiscard]] virtual std::uint64_t rejectedLongRangeRefreshCount() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t rejectedLongRangeReuseCount() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t rejectedCadenceRecordCount() const noexcept = 0;
  [[nodiscard]] virtual core::MemoryReport memoryReport() const = 0;
  [[nodiscard]] virtual const parallel::DecompositionRuntimeMeasurements&
  lastRuntimeDecompositionMeasurements() const noexcept = 0;
  virtual void shutdownMpiResources() = 0;

  virtual void restoreDecompositionEpoch(std::uint64_t decomposition_epoch) = 0;
  // Advances the decomposition epoch and explicitly invalidates authoritative
  // top-domain geometry freshness: after this call and before a reinstall,
  // authoritativeDomainGeometryMatches(...) is false. Ownership data,
  // migration state, and the epoch itself are not cleared.
  virtual void commitParticleDecompositionChange() = 0;
  // Install authoritative top-domain leaves at segment/restart/ownership-
  // change boundaries together with the physical source generation they
  // cover. Replaces both the stable seed leaf set (decomposition-local
  // identities/SFC intervals used as the refit seed) and the current
  // published leaf set. Geometry freshness is source_generation equality,
  // independent of ownership decomposition_epoch.
  virtual void installAuthoritativeTopDomainLeaves(
      std::span<const parallel::TopDomainLeaf> leaves,
      std::uint64_t source_generation) = 0;
  [[nodiscard]] virtual bool authoritativeDomainGeometryMatches(
      std::uint64_t source_generation) const noexcept = 0;
  virtual void importRestartForceCache(
      const io::GravityForceCachePersistentState& cache,
      const core::SimulationState& state) = 0;
  // Hierarchical KDK endpoint reuse (rank-local; callers agree collectively).
  // True when the last all-active synchronization endpoint left Tree and PM
  // split components for every row at the current source positions, with PM
  // from the integrator's committed field version: the next block may open
  // with these forces instead of re-evaluating at the same positions.
  [[nodiscard]] virtual bool hierarchicalEndpointForceReusable(
      const core::SimulationState& state, const core::IntegratorState& integrator) const = 0;
  // Multirate outer-PM cycle interior (TreePmPmReusePolicy::kMultirateOuterKick):
  // true when the Tree component (and total-force history) of every row is
  // current at the present source positions. The PM component is not
  // required: its field belongs to the open cycle's opening endpoint.
  [[nodiscard]] virtual bool hierarchicalShortRangeEndpointReusable(
      const core::SimulationState& state) const = 0;
  // Per-row squared nearest direct-pair separation from the latest Tree
  // solve (+infinity: none inside r_cut). Empty outside hierarchical KDK.
  [[nodiscard]] virtual std::span<const double> hierarchicalNearestPairR2() const noexcept = 0;

};

[[nodiscard]] std::unique_ptr<GravityRuntime> makeGravityRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const RuntimeServices& services,
    std::filesystem::path zoom_region_path = {});

}  // namespace cosmosim::workflows
