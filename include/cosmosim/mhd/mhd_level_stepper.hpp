#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_memory.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"

namespace cosmosim::mhd {

// One authoritative patch of a synchronized level. Gas and face storage are
// owned by the caller (ultimately core::SimulationState); the stepper only
// holds views plus its own governed scratch and rollback images.
struct MhdLevelPatchBinding {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t level = 0;
  // Lower corner in the level's integer cell lattice.
  std::array<std::int64_t, 3> origin_cell{};
  MhdPatchGeometry geometry{};
  MhdGasStateSoa* gas = nullptr;
  MhdFaceFluxMutableView faces{};
};

// Rectangular level domain with physical boundary kinds per axis. Allowed:
// kPeriodic (both sides), kOpen, kReflectingConductingWall.
struct MhdLevelDomain {
  std::array<std::int64_t, 3> cells{};
  std::array<MhdBoundaryKind, 3> lower{MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic,
                                       MhdBoundaryKind::kPeriodic};
  std::array<MhdBoundaryKind, 3> upper{MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic,
                                       MhdBoundaryKind::kPeriodic};
};

// Timeline sampler supplied by the authoritative cosmological timeline owner:
// returns (a, H) at a time inside the requested interval. Required for an
// in-stepper dt retry in comoving runs (the MHD module never recomputes a(t));
// without it a rejected comoving interval is returned to the scheduler.
using MhdTimelineSampler = std::function<std::pair<double, double>(double time_code)>;

struct MhdLevelAdvanceOptions {
  // nullptr selects the ideal profile (no non-ideal allocation or work).
  const MhdNonIdealRateOperator* nonideal = nullptr;
  MhdTimelineSampler timeline{};
  // When true, the time-integrated edge EMF of the whole accepted interval
  // (all D/H/I/H/D stages, SSPRK2 weights and RKL2 recurrence) is produced per
  // patch: Phi_end = Phi_begin - C(I). This is the AMR/MPI register input.
  bool produce_edge_integrals = false;
};

struct MhdLevelAdvanceResult {
  bool accepted = false;
  double dt_requested_code = 0.0;
  double dt_accepted_code = 0.0;
  // When !accepted: the scheduler must re-issue at most this interval.
  double dt_suggested_code = 0.0;
  std::int32_t retry_count = 0;
  std::string last_rejection;
  MhdStageCounters ideal{};
  MhdTransportStageDiagnostics transport{};
  std::uint64_t ideal_stage_evaluations_attempted = 0U;
  std::uint64_t ideal_stage_evaluations_accepted = 0U;
  double fallback_fraction = 0.0;
  bool fallback_fraction_gate_exceeded = false;
  MhdConservationBudget budget_begin{};
  MhdConservationBudget budget_end{};
  double divergence_max_normalized = 0.0;
  std::uint64_t scratch_capacity_bytes = 0U;
};

// Synchronized single-level MHD interval owner for a set of same-rank patches
// that tile the level domain with a uniform patch shape.
//
// Ownership contract (CHUI-MHD-REF v1.0 sections 5.2, 6.1, 7.5, 8.5, 15):
//  * One interval = D(dt/2) H(dt/2) I(dt) H(dt/2) D(dt/2), all stages of all
//    patches at a common stage time; the stepper is the single committer.
//  * Every stage halo is filled from the neighbours' SAME-stage state; shared
//    faces, boundary-face gas/Poynting fluxes and edge EMFs are reconciled from
//    the deterministic owner (lowest patch_id, then lowest local coordinates),
//    so each logical face/edge has one value and conservation is exact.
//  * Any failure restores every patch to the interval-start image (gas and
//    faces). Retryable failures (non-positive state, stiffness growth) halve
//    dt at most max_step_retries times; non-retryable failures (invalid
//    field, divergence contract, bad coefficient, stale geometry, unsupported
//    physics, memory) propagate after rollback.
//  * Multi-level, multi-rank and recursive operation are NOT implemented here:
//    mixed levels or kImportedMpi sides are rejected at construction.
class MhdLevelStepper {
 public:
  MhdLevelStepper(MhdNumericsPolicy policy, MhdLevelDomain domain,
                  std::vector<MhdLevelPatchBinding> patches);

  [[nodiscard]] std::size_t patchCount() const noexcept { return m_patches.size(); }
  [[nodiscard]] const MhdLevelPatchBinding& binding(std::size_t index) const {
    return m_patches.at(index).binding;
  }
  [[nodiscard]] const MhdPatchBoundary& patchBoundary(std::size_t index) const {
    return m_patches.at(index).boundary;
  }
  // Time-integrated edge EMF of the last accepted interval (only when the
  // interval ran with produce_edge_integrals).
  [[nodiscard]] const MhdEdgeEmfStorage& edgeIntegral(std::size_t index) const;

  [[nodiscard]] MhdLevelAdvanceResult advance(const MhdStageContext& interval,
                                              const MhdLevelAdvanceOptions& options);

  // Scheduler candidate: min over patches of the ideal fast CFL bound and, when
  // non-ideal is active, the largest interval whose Hall subcycle and RKL2
  // stage forecasts respect the configured per-interval cost caps.
  [[nodiscard]] double estimateStableTimeStep(const MhdStageContext& interval,
                                              const MhdLevelAdvanceOptions& options);

  [[nodiscard]] MhdConservationBudget levelBudget(double scale_factor) const;
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;

  // Governed admission of every stepper allocation BEFORE it exists (PERF
  // sections 5-6): the exact estimate for this patch set and the requested
  // mode is reserved as kPhaseResident (rollback + stage-0 images) and
  // kScratchArena (stage caches), checked against the 512 B/active-cell cap.
  // Refusal throws MHD_MEMORY_BUDGET_EXCEEDED with no state mutated. Once
  // admitted, advance() rejects any request needing more than was admitted
  // and verifies retained capacity after the first allocation.
  void admitMemory(core::MemoryGovernor& governor, const MhdLevelAdvanceOptions& mode,
                   const MhdResourceLimits& limits = {});
  [[nodiscard]] bool memoryAdmitted() const noexcept { return m_admission.admitted(); }

 private:
  struct PatchRuntime {
    MhdLevelPatchBinding binding{};
    MhdPatchBoundary boundary{};
    std::array<std::int64_t, 3> grid{};
    MhdStageWorkspace workspace{};
    MhdNonIdealWorkspace nonideal{};
    MhdGasStateSoa rollback_gas{};
    std::array<MhdDoubleArray, 3> rollback_faces{};
    MhdGasStateSoa stage0_gas{};
    std::array<MhdDoubleArray, 3> stage0_faces{};
    // Time-integrated EMF tracked as part of the magnetic state.
    MhdEdgeEmfStorage integral{};
    MhdEdgeEmfStorage stage0_integral{};
    MhdEdgeEmfStorage yjm2_integral{};
  };

  enum class NonIdealPlanKind : std::uint8_t { kNone, kExplicit, kRkl2 };

  void validateTopology();
  [[nodiscard]] std::size_t patchAtGrid(const std::array<std::int64_t, 3>& grid) const;
  [[nodiscard]] std::span<double> faceLane(PatchRuntime& patch, std::size_t axis) noexcept;

  void canonicalizeMirrors();
  void snapshotRollback();
  void restoreRollback();
  void snapshotStage0(bool gas_all_lanes);
  void prepareStage(double scale_factor);
  void fillNeighbourGhosts();
  void reconcileEdges();
  // Copies every shared face from its owner; returns the largest |difference|
  // that was overwritten (0 when the mirrors were already identical).
  double reconcileFaces();
  void reconcileBoundaryGasFluxes();
  void reconcileBoundaryEnergyFluxes();
  void applyCt(PatchRuntime& patch, const MhdEdgeEmfStorage& emf, double weight);
  void combineMagnetic(PatchRuntime& patch, double alpha, double beta);

  void idealSsprk2(const MhdStageContext& interval, MhdLevelAdvanceResult& result);
  [[nodiscard]] double evaluateNonIdealRate(MhdNonIdealOperatorKind kind, double scale_factor,
                                            MhdLevelAdvanceResult& result);
  void explicitNonIdealSubstep(MhdNonIdealOperatorKind kind, const MhdStageContext& sub,
                               MhdLevelAdvanceResult& result);
  void rkl2Cycle(const MhdStageContext& sub, std::int32_t stages, MhdLevelAdvanceResult& result);
  void diffusionHalf(const MhdStageContext& sub, MhdLevelAdvanceResult& result);
  void hallHalf(const MhdStageContext& sub, MhdLevelAdvanceResult& result);
  void forecastNonIdeal(const MhdStageContext& interval, MhdLevelAdvanceResult& result);
  void validateAccepted(double scale_factor, MhdLevelAdvanceResult& result);

  MhdNumericsPolicy m_policy;
  MhdStepOperator m_operator;
  MhdLevelDomain m_domain;
  std::array<std::int64_t, 3> m_patch_cells{};
  std::array<std::int64_t, 3> m_grid_count{};
  std::vector<PatchRuntime> m_patches;
  std::vector<std::size_t> m_grid_lookup;
  // Interval-scoped state.
  const MhdNonIdealRateOperator* m_nonideal = nullptr;
  bool m_track_integral = false;
  double m_forecast_dt_diffusion = 0.0;
  double m_forecast_dt_hall = 0.0;
  std::int64_t m_interval_rkl_stages = 0;
  std::int64_t m_interval_hall_subcycles = 0;
  MhdMemoryAdmission m_admission{};
  MhdMemoryRequest m_admitted_request{};
  bool m_admission_committed = false;
};

}  // namespace cosmosim::mhd
