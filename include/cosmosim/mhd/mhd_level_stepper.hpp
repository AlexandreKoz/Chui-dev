#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/mhd/mhd_geometry.hpp"
#include "cosmosim/mhd/mhd_memory.hpp"
#include "cosmosim/mhd/mhd_state.hpp"
#include "cosmosim/mhd/mhd_step_operator.hpp"
#include "cosmosim/mhd/mhd_transport.hpp"

namespace cosmosim::mhd {

// One authoritative LOCAL patch of a synchronized leaf hierarchy. Gas and face
// storage are owned by the caller (ultimately core::SimulationState); the
// stepper only holds views plus its own governed scratch and rollback images.
struct MhdLevelPatchBinding {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t level = 0;
  // Lower corner in the integer cell lattice of `level` (level L lattice has
  // 2^(L - base_level) times the base-level cell count per axis).
  std::array<std::int64_t, 3> origin_cell{};
  MhdPatchGeometry geometry{};
  MhdGasStateSoa* gas = nullptr;
  MhdFaceFluxMutableView faces{};
  // Owning MPI rank; must equal the transport rank for a local binding.
  std::int32_t owner_rank = 0;
};

// Descriptor of a patch owned by ANOTHER rank (global topology directory
// entry; no storage). Every rank must pass the same global directory
// (local bindings + remote descriptors) so exchange plans agree exactly.
struct MhdRemotePatchDescriptor {
  std::uint64_t patch_id = 0U;
  std::uint64_t geometry_generation = 0U;
  std::int32_t level = 0;
  std::array<std::int64_t, 3> origin_cell{};
  std::array<std::int64_t, 3> cells{};
  std::int32_t owner_rank = 0;
};

// Rectangular domain of the BASE level with physical boundary kinds per axis.
// Allowed: kPeriodic (both sides), kOpen, kReflectingConductingWall. Finer
// levels cover the same physical domain with ratio-2 lattices.
struct MhdLevelDomain {
  std::array<std::int64_t, 3> cells{};
  std::array<MhdBoundaryKind, 3> lower{MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic,
                                       MhdBoundaryKind::kPeriodic};
  std::array<MhdBoundaryKind, 3> upper{MhdBoundaryKind::kPeriodic, MhdBoundaryKind::kPeriodic,
                                       MhdBoundaryKind::kPeriodic};
  std::int32_t base_level = 0;
  // Comoving position of the lattice origin (patch comoving origins must equal
  // origin_comoving + origin_cell * spacing).
  std::array<double, 3> origin_comoving{};
};

// Collective transport used for remote halo, owner reconciliation and
// coarse/fine override payloads, and for collective error agreement. The MHD
// module never calls MPI directly (core-only dependency); the workflow layer
// supplies an implementation on top of parallel::exchangeBoundedAlltoallBytes.
// Every method is COLLECTIVE: all ranks of the hierarchy call it in the same
// order (zero-patch ranks included).
class MhdHaloTransport {
 public:
  virtual ~MhdHaloTransport() = default;
  [[nodiscard]] virtual int rank() const noexcept = 0;
  [[nodiscard]] virtual int size() const noexcept = 0;
  // send[r] = payload for rank r (send[rank()] is ignored). Returns recv[r].
  [[nodiscard]] virtual std::vector<std::vector<std::uint8_t>> exchange(
      const std::vector<std::vector<std::uint8_t>>& send) = 0;
  [[nodiscard]] virtual std::uint64_t allreduceMaxUint64(std::uint64_t value) = 0;
  [[nodiscard]] virtual double allreduceMinDouble(double value) = 0;
  [[nodiscard]] virtual double allreduceMaxDouble(double value) = 0;
  virtual void allreduceSumDoubles(std::span<double> values) = 0;
};

// Single-rank transport (no communication).
class MhdSerialHaloTransport final : public MhdHaloTransport {
 public:
  [[nodiscard]] int rank() const noexcept override { return 0; }
  [[nodiscard]] int size() const noexcept override { return 1; }
  [[nodiscard]] std::vector<std::vector<std::uint8_t>> exchange(
      const std::vector<std::vector<std::uint8_t>>& send) override;
  [[nodiscard]] std::uint64_t allreduceMaxUint64(std::uint64_t value) override { return value; }
  [[nodiscard]] double allreduceMinDouble(double value) override { return value; }
  [[nodiscard]] double allreduceMaxDouble(double value) override { return value; }
  void allreduceSumDoubles(std::span<double>) override {}
};

// Timeline sampler supplied by the authoritative cosmological timeline owner:
// returns (a, H) at a time inside the requested interval. Used only by the
// standalone kInternalTimeline retry mode.
using MhdTimelineSampler = std::function<std::pair<double, double>(double time_code)>;

// Who owns the retry of a rejected interval.
//  * kExternal (default, mandatory in production): a retryable failure rolls
//    back every patch and returns accepted=false with dt_suggested_code; the
//    outer scheduler owns the timeline and retries the WHOLE interval
//    (including KDK/source state). The stepper never commits a shortened
//    interval on its own.
//  * kInternalTimeline: standalone/reference use only (tests, benchmarks):
//    the stepper halves dt and retries internally, using the supplied
//    timeline sampler or a constant timeline.
enum class MhdRetryMode : std::uint8_t { kExternal = 0, kInternalTimeline = 1 };

struct MhdLevelStepperConfig {
  // Upper bound on simultaneously active patch workspaces (PERF section 6:
  // at most eight 32^3-equivalent patches by default).
  std::uint32_t batch_patches = 8U;
  MhdRetryMode retry_mode = MhdRetryMode::kExternal;
  // Non-owning; nullptr selects MhdSerialHaloTransport. Must outlive the stepper.
  MhdHaloTransport* transport = nullptr;
  // Production: advance()/estimateStableTimeStep() refuse to allocate before
  // admitMemory() succeeded.
  bool require_memory_admission = false;
  // Hierarchy epoch written into every exchange payload header and checked by
  // the receiver (MHD_MPI_EPOCH_MISMATCH on disagreement).
  std::uint64_t hierarchy_epoch = 0U;
  // Communication budget for one exchange phase (send + receive payloads),
  // PERF section 8.1 (64 MiB default, 128 MiB explicit maximum).
  std::uint64_t comm_budget_bytes = k_mhd_default_comm_budget_bytes_per_rank;
};

struct MhdLevelAdvanceOptions {
  // nullptr selects the ideal profile (no non-ideal allocation or work).
  const MhdNonIdealRateOperator* nonideal = nullptr;
  MhdTimelineSampler timeline{};
  // When true, the time-integrated edge EMF of the whole accepted interval
  // (all D/H/I/H/D stages, SSPRK2 weights and RKL2 recurrence) is produced per
  // patch: Phi_end = Phi_begin - C(I).
  bool produce_edge_integrals = false;
  // When true, the accepted-interval boundary-plane flux integrals are
  // produced per patch: gas I_gas = sum_stage w F/a (six lanes; the cell
  // update is -/+ I/dx) and non-ideal Poynting energy I_E = sum_stage w F_E,
  // combined with exactly the SSPRK2 / explicit / RKL2 / Hall recurrences
  // that act on the state. External register input and boundary-energy audit.
  bool produce_boundary_flux_integrals = false;
};

struct MhdLevelAdvanceResult {
  bool accepted = false;
  double dt_requested_code = 0.0;
  double dt_accepted_code = 0.0;
  // When !accepted: the owner must re-issue at most this interval.
  double dt_suggested_code = 0.0;
  std::int32_t retry_count = 0;
  std::uint32_t rejected_attempts = 0U;
  std::string last_rejection;
  // Counters of the ACCEPTED attempt only (physics budgets).
  MhdStageCounters ideal{};
  MhdTransportStageDiagnostics transport{};
  // Counters over every attempt, accepted or rolled back (cost accounting).
  MhdStageCounters ideal_attempted{};
  MhdTransportStageDiagnostics transport_attempted{};
  std::uint64_t ideal_stage_evaluations_attempted = 0U;
  std::uint64_t ideal_stage_evaluations_accepted = 0U;
  std::uint64_t ghost_prolongation_slope_fallbacks = 0U;
  double fallback_fraction = 0.0;
  bool fallback_fraction_gate_exceeded = false;
  MhdConservationBudget budget_begin{};
  MhdConservationBudget budget_end{};
  double divergence_max_normalized = 0.0;
  std::uint64_t scratch_capacity_bytes = 0U;
  std::uint64_t mpi_bytes_sent = 0U;
  std::uint64_t mpi_bytes_received = 0U;
  std::uint64_t exchange_phases = 0U;
};

// Honest resource plan of a stepper (PERF sections 4-9). Level-resident bytes
// are unavoidable per local patch for a synchronized stage (rollback image,
// stage-0 image, stage increment, edge EMF, boundary records, non-ideal and
// super-stepping lanes, optional integrals, remote halo imports); batch
// scratch is the reusable tile-local workspace of `batch_patches` slots.
struct MhdLevelMemoryPlan {
  std::uint64_t level_resident_bytes = 0U;
  std::uint64_t batch_scratch_bytes = 0U;
  std::uint64_t topology_bytes = 0U;
  std::uint64_t communication_bytes = 0U;
  std::uint64_t canonical_face_flux_bytes = 0U;
  std::uint64_t active_cells = 0U;
  std::uint32_t batch_patches = 0U;
  [[nodiscard]] std::uint64_t phasePlusScratchBytes() const noexcept {
    return level_resident_bytes + batch_scratch_bytes;
  }
  [[nodiscard]] double bytesPerActiveCell() const noexcept {
    return active_cells == 0U ? 0.0
                              : static_cast<double>(phasePlusScratchBytes()) /
            static_cast<double>(active_cells);
  }
};

// Synchronized MHD interval owner for a LEAF-TILED patch hierarchy (one or
// more ratio-2 levels, arbitrary rectangular patch extents, nonzero origins,
// sparse fine coverage, local and remote patches, zero-patch ranks).
//
// Ownership contract (CHUI-MHD-REF v1.0 sections 5.2, 6.1, 7.3-7.5, 8.5, 15):
//  * All levels advance with ONE common dt and common stage times
//    (synchronized AMR sweeps, amr.mhd_subcycling=false). One interval =
//    D(dt/2) H(dt/2) I(dt) H(dt/2) D(dt/2); the stepper is the single committer.
//  * Every stage reads only the immutable stage state (canonical gas + faces);
//    patch kernels run in bounded batches whose results go to per-patch stage
//    increments, so a batch can never read an already-updated neighbour.
//  * Ghost data: same-level leaves copy, finer leaves restrict (conserved
//    average + summed face flux), coarser leaves prolong (minmod-limited
//    primitives; divergence-irrelevant ghost faces by face sampling). Remote
//    ghost data are exchanged once per stage through the transport.
//  * Shared same-level faces, edges, gas and Poynting boundary records have
//    one deterministic owner (lowest (patch_id, local coordinates)); coarse
//    interface records, edges and faces take the fine-side values (mean of the
//    four fine records, sum of the two fine edges, sum of the four fine faces).
//    Gas mass/momentum/energy and magnetic flux are therefore conserved and
//    D*C = 0 holds on the composite leaf complex at every stage.
//  * Any failure restores every local patch from preallocated rollback images
//    with no-throw copies (RAII); collective error agreement precedes every
//    exchange so all ranks roll back together. Retryable failures are handed
//    to the owner (kExternal) or retried internally (kInternalTimeline).
//  * Distributed recursive subcycling is NOT implemented here (fail closed at
//    the workflow; REF 7.4 tier mhd_ref_amr_subcycled).
class MhdLevelStepper {
 public:
  MhdLevelStepper(MhdNumericsPolicy policy, MhdLevelDomain domain,
                  std::vector<MhdLevelPatchBinding> patches,
                  std::vector<MhdRemotePatchDescriptor> remote_patches = {},
                  MhdLevelStepperConfig config = {});
  // Not copyable or movable: the stepper holds views into caller storage and
  // a pointer to its own serial transport. Hold it by value at a fixed
  // address or through std::unique_ptr (prvalue initialization is fine).
  MhdLevelStepper(const MhdLevelStepper&) = delete;
  MhdLevelStepper& operator=(const MhdLevelStepper&) = delete;
  MhdLevelStepper(MhdLevelStepper&&) = delete;
  MhdLevelStepper& operator=(MhdLevelStepper&&) = delete;
  ~MhdLevelStepper() = default;

  // Local patches, in the order supplied.
  [[nodiscard]] std::size_t patchCount() const noexcept { return m_local.size(); }
  [[nodiscard]] std::size_t globalPatchCount() const noexcept { return m_topo.size(); }
  [[nodiscard]] std::int32_t levelCount() const noexcept { return m_level_count; }
  [[nodiscard]] const MhdLevelPatchBinding& binding(std::size_t index) const {
    return m_local.at(index).binding;
  }
  [[nodiscard]] const MhdPatchBoundary& patchBoundary(std::size_t index) const {
    return m_local.at(index).boundary;
  }
  // Time-integrated edge EMF of the last accepted interval (only when the
  // interval ran with produce_edge_integrals).
  [[nodiscard]] const MhdEdgeEmfStorage& edgeIntegral(std::size_t index) const;
  // Accepted-interval boundary flux integrals (produce_boundary_flux_integrals).
  [[nodiscard]] const MhdBoundaryFaceFluxes& boundaryGasFluxIntegral(std::size_t index) const;
  [[nodiscard]] const MhdBoundaryFaceEnergyFluxes& boundaryEnergyFluxIntegral(
      std::size_t index) const;
  // Number of coarse/fine interface faces whose coarse records are overridden.
  [[nodiscard]] std::uint64_t coarseFineInterfaceFaceCount() const noexcept {
    return m_interface_face_count;
  }

  // COLLECTIVE.
  [[nodiscard]] MhdLevelAdvanceResult advance(const MhdStageContext& interval,
                                              const MhdLevelAdvanceOptions& options);

  // COLLECTIVE scheduler candidate (rank-global minimum): ideal fast CFL and,
  // with non-ideal active, the largest interval whose Hall subcycle and RKL2
  // stage forecasts respect the per-interval cost caps. Never mutates
  // canonical gas or face state (scratch only; admitted scratch when the
  // stepper requires admission).
  [[nodiscard]] double estimateStableTimeStep(const MhdStageContext& interval,
                                              const MhdLevelAdvanceOptions& options);

  // COLLECTIVE read-only validation of the canonical composite state: shared
  // same-level face mirrors and coarse interface faces versus the sum of fine
  // faces. Returns the rank-global maximum normalized mismatch; throws
  // MHD_INVALID_FIELD when it exceeds roundoff. Nothing is written.
  double validateCanonicalConsistency();

  // COLLECTIVE (rank-global sums; divergence is the rank-global maximum).
  [[nodiscard]] MhdConservationBudget levelBudget(double scale_factor);
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const noexcept;

  [[nodiscard]] MhdLevelMemoryPlan planMemory(const MhdLevelAdvanceOptions& mode,
                                              const MhdResourceLimits& limits = {}) const;
  // Governed admission of every stepper allocation BEFORE it exists (PERF
  // sections 5-9): level-resident state as kPhaseResident, batch scratch as
  // kScratchArena, topology metadata as kPhaseResident, exchange buffers as
  // kCommunication, all checked against the 512 B/active-cell cap and the
  // communication budget. Refusal throws MHD_MEMORY_BUDGET_EXCEEDED with no
  // state mutated. After admission, advance() rejects any mode needing more
  // than was admitted and commits measured capacity after first allocation.
  void admitMemory(core::MemoryGovernor& governor, const MhdLevelAdvanceOptions& mode,
                   const MhdResourceLimits& limits = {});
  [[nodiscard]] bool memoryAdmitted() const noexcept { return m_admission.admitted(); }

 private:
  static constexpr std::size_t k_npos = static_cast<std::size_t>(-1);

  // ---- Topology --------------------------------------------------------
  struct Box {
    std::array<std::int64_t, 3> lo{};
    std::array<std::int64_t, 3> hi{};
    [[nodiscard]] bool empty() const noexcept {
      return lo[0] >= hi[0] || lo[1] >= hi[1] || lo[2] >= hi[2];
    }
    [[nodiscard]] std::int64_t volume() const noexcept {
      return empty() ? 0 : (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
    }
  };

  struct TopoPatch {
    std::uint64_t patch_id = 0U;
    std::uint64_t geometry_generation = 0U;
    std::int32_t level = 0;
    std::int32_t owner_rank = 0;
    Box box{};
    std::size_t local = k_npos;
  };

  enum class ItemKind : std::uint8_t {
    // Halo (stage state -> ghost cells / ghost staggered faces / samples).
    kGhostCellDirect = 0,
    kGhostCellRestricted,
    kGhostFaceDirect,
    kGhostFaceRestricted,
    kSampleCellDirect,
    kSampleCellRestricted,
    kSampleFaceDirect,
    // Same-level owner reconciliation (src = owner candidate).
    kFaceOwner,
    kEdgeOwner,
    kGasRecordOwner,
    kEnergyRecordOwner,
    // Coarse <- fine interface overrides (src = fine leaf).
    kFaceFine,
    kEdgeFine,
    kGasRecordFine,
    kEnergyRecordFine,
  };

  // A box-to-box exchange item. Index boxes are in the DESTINATION patch's
  // local coordinates of the item family (cells; faces with the node along
  // `axis`; edges with the cell along `axis` and nodes along the others).
  // Source local coordinates: src = scale * dst + src_offset.
  struct Item {
    ItemKind kind = ItemKind::kGhostCellDirect;
    std::uint8_t axis = 0U;
    std::uint8_t scale = 1U;
    std::uint32_t src = 0U;
    std::uint32_t dst = 0U;
    Box dst_box{};
    std::array<std::int64_t, 3> src_offset{};
    // Halo items with a remote source: offset (doubles) in the destination's
    // import buffer. Sample items: slot index in the destination sample table.
    std::uint64_t import_offset = 0U;
    std::uint32_t slot = 0U;
  };

  // Prolonged ghost cell of a fine patch from coarse samples.
  struct ProlongCell {
    MhdIndex3 dst{};
    std::int32_t center = -1;
    std::array<std::int32_t, 6> neighbour{-1, -1, -1, -1, -1, -1};  // (axis, -/+)
    std::array<std::int8_t, 3> half{};  // -1 lower fine half, +1 upper
  };
  // Prolonged ghost staggered face of a fine patch from coarse face samples.
  struct ProlongFace {
    std::uint8_t axis = 0U;
    MhdIndex3 dst{};
    std::int32_t first = -1;
    std::int32_t second = -1;  // >= 0: odd node, mean of first/second
  };

  struct HaloPlan {
    std::vector<Item> items;            // ghost boxes + samples, deterministic order
    std::vector<ProlongCell> prolong_cells;
    std::vector<ProlongFace> prolong_faces;
    std::uint32_t cell_sample_count = 0U;
    std::uint32_t face_sample_count = 0U;
    std::uint64_t import_doubles = 0U;  // remote-source values
  };

  // ---- Per-local-patch runtime -----------------------------------------
  struct PatchRuntime {
    MhdLevelPatchBinding binding{};
    MhdPatchBoundary boundary{};
    std::uint32_t topo = 0U;
    HaloPlan halo{};
    MhdDoubleArray halo_import;
    // Level-resident synchronized-stage state.
    MhdGasStateSoa rollback_gas{};
    std::array<MhdDoubleArray, 3> rollback_faces{};
    MhdGasStateSoa stage0_gas{};
    std::array<MhdDoubleArray, 3> stage0_faces{};
    MhdGasStateSoa increment{};
    MhdEdgeEmfStorage emf{};
    MhdBoundaryFaceFluxes records{};
    // Non-ideal level-resident lanes.
    MhdDoubleArray energy_rate;
    MhdBoundaryFaceEnergyFluxes energy_records{};
    std::array<MhdDoubleArray, 3> yjm2_faces{};
    MhdDoubleArray yjm2_energy;
    MhdDoubleArray energy_rate_y0;
    MhdEdgeEmfStorage emf_y0{};
    // Optional accepted-interval integrals (edge EMF and boundary fluxes).
    MhdEdgeEmfStorage integral{};
    MhdEdgeEmfStorage stage0_integral{};
    MhdEdgeEmfStorage yjm2_integral{};
    MhdBoundaryFaceFluxes gas_integral{};
    MhdBoundaryFaceFluxes stage0_gas_integral{};
    MhdBoundaryFaceEnergyFluxes energy_integral{};
    MhdBoundaryFaceEnergyFluxes stage0_energy_integral{};
    MhdBoundaryFaceEnergyFluxes yjm2_energy_integral{};
    MhdBoundaryFaceEnergyFluxes energy_records_y0{};
    MhdBoundaryFaceEnergyFluxes rollback_energy_integral{};
    MhdBoundaryFaceFluxes rollback_gas_integral{};
  };

  // Reusable batch workspace (one per active slot).
  struct BatchSlot {
    MhdStageWorkspace workspace{};
    MhdNonIdealWorkspace nonideal{};
    MhdDoubleArray cell_samples;  // 9 lanes per sample
    MhdDoubleArray face_samples;
    std::vector<double> values;   // local halo item values (capacity retained)
    MhdStageCounters counters{};
    MhdTransportStageDiagnostics transport{};
    std::uint64_t slope_fallbacks = 0U;
    double stable_dt = 0.0;
  };

  enum class NonIdealPhase : std::uint8_t { kEmf = 0, kEnergy = 1 };

  // Topology construction (mhd_level_topology.cpp).
  void buildTopology(std::vector<MhdRemotePatchDescriptor> remote);
  void buildSpatialIndex();
  [[nodiscard]] std::array<std::int64_t, 3> levelCells(std::int32_t level) const noexcept;
  [[nodiscard]] bool periodicAxis(std::size_t axis) const noexcept;
  // Patches of `level` whose (shifted) box intersects `query`, with the
  // periodic image shift; deterministic order (patch_id ascending, shift).
  void queryLevel(std::int32_t level, const Box& query,
                  std::vector<std::pair<std::uint32_t, std::array<std::int64_t, 3>>>& out) const;
  [[nodiscard]] std::int64_t findCellHolder(std::int32_t level,
                                            const std::array<std::int64_t, 3>& cell,
                                            std::array<std::int64_t, 3>& shift) const;
  void planHalo(std::uint32_t topo_index, HaloPlan& plan) const;
  void planReconciliation();
  void planCrossLevel();
  void finalizeItemLists();
  [[nodiscard]] MhdPatchBoundary boundaryForPatch(const TopoPatch& patch) const;

  // Exchange framework.
  template <class Packer, class Applier>
  void runItems(std::uint64_t phase_id, std::span<const Item> items, Packer&& pack,
                Applier&& apply, bool remote_only);
  void agree(const std::exception_ptr& local_failure, std::string_view phase);
  template <class Fn>
  void collective(std::string_view phase, Fn&& fn);
  [[nodiscard]] std::uint64_t topologyCapacityBytes() const noexcept;
  void exchangeHalo(double scale_factor);
  void computeItemValues(const Item& item, double scale_factor, std::vector<double>& out) const;
  [[nodiscard]] std::size_t itemValueCount(const Item& item) const noexcept;

  // Batch processing.
  template <class Fn>
  void forEachLocalPatchBatched(Fn&& fn);
  void prepareSlot(std::size_t patch, BatchSlot& slot, double scale_factor);
  void fillExternalGhosts(std::size_t patch, BatchSlot& slot, double scale_factor);

  // Reconciliation and overrides (collective).
  void reconcileGasRecords();
  void reconcileEnergyRecords();
  void reconcileEdges();
  void reconcileFaces(bool compare_only, double* max_mismatch);
  void applyItemsWrite(std::span<const Item> items, std::uint64_t phase_id);

  // Transaction.
  void ensureAllocated(const MhdLevelAdvanceOptions& options);
  void snapshotRollback();
  void restoreRollback() noexcept;
  void snapshotStage0(bool gas_all_lanes);
  void combineMagnetic(PatchRuntime& patch, double alpha, double beta);
  void applyCt(PatchRuntime& patch, const MhdEdgeEmfStorage& emf, double weight);

  // Operators.
  void idealStage(double alpha, double beta, double weight, double scale_factor,
                  double hubble, bool expansion, MhdStageCounters& counters,
                  MhdLevelAdvanceResult& result);
  void idealSsprk2(const MhdStageContext& interval, MhdStageCounters& counters,
                   MhdLevelAdvanceResult& result);
  [[nodiscard]] double evaluateNonIdealRate(MhdNonIdealOperatorKind kind, double scale_factor,
                                            MhdTransportStageDiagnostics& diagnostics,
                                            MhdLevelAdvanceResult& result);
  void explicitNonIdealSubstep(MhdNonIdealOperatorKind kind, const MhdStageContext& sub,
                               MhdTransportStageDiagnostics& diagnostics,
                               MhdLevelAdvanceResult& result);
  void rkl2Cycle(const MhdStageContext& sub, std::int32_t stages,
                 MhdTransportStageDiagnostics& diagnostics, MhdLevelAdvanceResult& result);
  void diffusionHalf(const MhdStageContext& sub, MhdTransportStageDiagnostics& diagnostics,
                     MhdLevelAdvanceResult& result);
  void hallHalf(const MhdStageContext& sub, MhdTransportStageDiagnostics& diagnostics,
                MhdLevelAdvanceResult& result);
  void forecastNonIdeal(const MhdStageContext& interval,
                        MhdTransportStageDiagnostics& diagnostics, MhdLevelAdvanceResult& result,
                        bool enforce_caps);
  void validateAccepted(double scale_factor, MhdLevelAdvanceResult& result);
  void requireModeAdmitted(const MhdLevelAdvanceOptions& options) const;
  [[nodiscard]] std::span<double> faceLane(PatchRuntime& patch, std::size_t axis) noexcept;
  [[nodiscard]] std::span<const double> faceLane(const PatchRuntime& patch,
                                                 std::size_t axis) const noexcept;

  MhdNumericsPolicy m_policy;
  MhdStepOperator m_operator;
  MhdLevelDomain m_domain;
  MhdLevelStepperConfig m_config;
  MhdSerialHaloTransport m_serial_transport{};
  MhdHaloTransport* m_transport = nullptr;
  std::int32_t m_level_count = 1;
  std::int32_t m_max_level = 0;

  std::vector<TopoPatch> m_topo;                       // sorted (level, patch_id)
  std::vector<PatchRuntime> m_local;                   // caller order
  // Sparse spatial index per level: sorted (bin key, topo index).
  std::vector<std::vector<std::pair<std::array<std::int64_t, 3>, std::uint32_t>>> m_bins;
  std::vector<std::array<std::int64_t, 3>> m_bin_size;
  // Phase item lists (only items with a local source or a local destination).
  std::vector<Item> m_halo_items;          // remote-source halo values (exchange only)
  std::vector<Item> m_face_self_items;
  std::vector<Item> m_face_cross_items;
  std::vector<Item> m_edge_self_items;
  std::vector<Item> m_edge_cross_items;
  std::vector<Item> m_record_self_items;    // gas + energy share geometry
  std::vector<Item> m_record_cross_items;
  std::vector<Item> m_face_fine_items;
  std::vector<std::vector<Item>> m_edge_fine_items;   // by coarse level, descending use
  std::vector<Item> m_record_fine_items;
  // Energy-record variants (same geometry, one value per face).
  std::vector<Item> m_energy_self_items;
  std::vector<Item> m_energy_cross_items;
  std::vector<Item> m_energy_fine_items;
  std::uint64_t m_interface_face_count = 0U;

  std::vector<BatchSlot> m_slots;
  std::vector<double> m_exchange_scratch;
  // Exchange staging (kCommunication; capacity retained across phases).
  std::vector<std::vector<double>> m_send_values;
  std::vector<std::vector<std::uint8_t>> m_send_bytes;
  std::vector<std::vector<double>> m_recv_values;
  std::uint64_t m_phase_counter = 0U;
  std::uint64_t m_bytes_sent = 0U;
  std::uint64_t m_bytes_received = 0U;
  std::uint64_t m_exchange_phases = 0U;
  double m_compare_max = 0.0;

  // Interval-scoped state.
  const MhdNonIdealRateOperator* m_nonideal = nullptr;
  bool m_track_integral = false;
  bool m_track_flux_integrals = false;
  bool m_allocated_nonideal = false;
  bool m_allocated_super = false;
  bool m_allocated_integral = false;
  bool m_allocated_flux_integrals = false;
  double m_forecast_dt_diffusion = 0.0;
  double m_forecast_dt_hall = 0.0;
  std::int64_t m_interval_rkl_stages = 0;
  std::int64_t m_interval_hall_subcycles = 0;
  MhdMemoryAdmission m_admission{};
  core::MemoryReservation m_comm_reservation{};
  core::MemoryReservation m_topology_reservation{};
  MhdLevelAdvanceOptions m_admitted_mode{};
  bool m_admitted_nonideal = false;
  bool m_admitted_super = false;
  MhdLevelMemoryPlan m_admitted_plan{};
  bool m_admission_committed = false;
};

}  // namespace cosmosim::mhd
