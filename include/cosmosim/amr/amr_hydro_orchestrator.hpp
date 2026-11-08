#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <span>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/amr/amr_ghost_fill.hpp"
#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace cosmosim::amr {

enum class ProductionAmrHydroSweepMode : std::uint8_t {
  kSynchronized = 0,
  kLocalSubcycling = 1,
};

// Derived thermodynamic caches written by AMR operations that mutate gas
// (scatter after a patch advance, reflux, refine/derefine writes). With a
// zero temperature factor and no closure the legacy behaviour is kept
// (ideal P/c_s, temperature_code = internal_energy_code) for focused callers;
// production passes the runtime-wide temperature factor and the active
// effective-EOS closure so P, c_s and T are consistent with the selected
// model immediately after the mutation (no stale cache can reach the CFL).
struct AmrThermodynamicCachePolicy {
  double temperature_k_per_internal_energy_code = 0.0;
  const hydro::HydroThermodynamicClosure* closure = nullptr;
  double scale_factor = 1.0;
  double redshift = 0.0;
  // Dual energy at reflux: when the corrected (E - K) is below eta E the
  // reflux keeps the pre-reflux specific internal energy instead of the
  // cancellation-dominated E - K, and ledgers the difference in
  // RefluxDiagnostics::dual_energy_reflux_thermal_correction_code.
  bool dual_energy = false;
  double dual_energy_eta = 1.0e-3;
};

// CHUI Part 3 Milestone I: reflux interval of one coarse/fine level pair
// (coarse level = index into ProductionAmrHydroOptions::
// reflux_intervals_by_coarse_level). A level-L sweep records registers of two
// pairs: (L-1, L) as the fine side (substep k of the level-(L-1) interval)
// and (L, L+1) as the coarse side (its own interval, substep 0).
struct AmrLevelRefluxInterval {
  double interval_start_code = 0.0;
  double interval_end_code = 0.0;
  double coarse_dt_code = 0.0;
  std::uint32_t substeps_per_interval = 1;
  std::uint32_t substep_index = 0;
  bool valid = false;
};

// Workflow-owned smooth baryonic source seam for recursive AMR. The AMR
// library owns level clocks and calls this operator on every rank before and
// after each level transport interval. The implementation may perform MPI
// collectives; empty local row spans are therefore meaningful participation,
// not a reason to skip the call.
// Phase identity of one level interval [t_b, t_e] of level L (patch 51, T02):
//  * kPreHydroHalf: level L and every finer level are synchronized at t_b
//    (children have not started); called BEFORE the level-L temporal history
//    start is captured.
//  * kPostHydroHalf: level L at t_e, finer levels still at t_b; called before
//    the level-L history end is captured.
//  * kPostChildrenSynchronization: called only for a level that HAS finer
//    levels, after all children substeps, the complete reflux and the
//    retirement of the level-L temporal history; level L and every finer
//    level are synchronized at t_e. Operators that couple level L to finer
//    cells (coarse/fine interface diffusion and SGS faces) must use only
//    kPreHydroHalf and this phase, where both sides are contemporaneous.
enum class RecursiveAmrSourcePhase : std::uint8_t {
  kPreHydroHalf = 0,
  kPostHydroHalf = 1,
  kPostChildrenSynchronization = 2,
};

class RecursiveAmrLevelSourceOperator {
 public:
  virtual ~RecursiveAmrLevelSourceOperator() = default;
  virtual void apply(
      core::SimulationState& state,
      std::span<const std::uint32_t> level_rows,
      const hydro::HydroUpdateContext& level_update,
      int level,
      RecursiveAmrSourcePhase phase) = 0;
};

struct ProductionAmrHydroOptions {
  hydro::HydroBoundaryKind physical_boundary_kind = hydro::HydroBoundaryKind::kOpen;
  double adiabatic_index = 5.0 / 3.0;
  double density_floor = 1.0e-10;
  double pressure_floor = 1.0e-10;
  ProductionAmrHydroSweepMode sweep_mode = ProductionAmrHydroSweepMode::kSynchronized;
  std::uint32_t refinement_ratio = 2;
  bool persist_incomplete_flux_registers = false;
  double reflux_interval_start_code = 0.0;
  double reflux_interval_end_code = 0.0;
  double reflux_coarse_dt_code = 0.0;
  std::uint32_t expected_fine_substeps = 1;
  std::uint32_t fine_substep_index = 0;
  // When non-empty, register intervals are taken per entry from this table
  // (indexed by the register's coarse level) instead of the scalar fields
  // above; an entry whose level has no valid interval fails closed.
  std::span<const AmrLevelRefluxInterval> reflux_intervals_by_coarse_level{};
  // Explicit AMR hydro state/ghost time.  Ghosts are consumed at update
  // step-start before the in-solver MUSCL-Hancock predictor evolves them.
  double state_time_code = 0.0;
  double ghost_fill_time_code = 0.0;
  bool enable_temporal_coarse_to_fine = false;
  hydro::HydroActiveBatchPolicy active_batch_policy{};
  // Optional production conversion factor used to derive cooling/source
  // inputs patch-locally from canonical SimulationState instead of retaining
  // full-population thermodynamic helper vectors. A zero value preserves the
  // legacy span-provided source-context path used by focused tests.
  double physical_density_cgs_per_density_code = 0.0;
  // Internal subcycling conversion, applied only to patch-local source lanes.
  // The global source context remains a borrowed immutable view.
  double source_density_rescale = 1.0;
  // Production regrid is an explicit memory transaction. Regrid entry points
  // fail closed when this authority is absent; ordinary hydro stepping does
  // not require it. Runtime composition must pass the process MemoryGovernor.
  core::MemoryGovernor* regrid_memory_governor = nullptr;
  // Non-owning synchronous hook. Null retains transport-only recursive AMR.
  RecursiveAmrLevelSourceOperator* recursive_level_source_operator = nullptr;
  // Spatial reconstruction and limiter used for every patch advance. The
  // time-integration family (MUSCL-Hancock / CTU) and dual-energy policy are
  // owned by the HydroCoreSolver passed to the advance entry points. PPM uses
  // the multi-layer normal shell. CTU reconstructs imported/physical ghost
  // transverse increments from complete tangential face-slab stencils; only
  // the edge/corner region where the selected reconstruction's full stencil
  // crosses a second interface uses the explicit counted normal-only fallback.
  hydro::HydroReconstructionScheme reconstruction_scheme = hydro::HydroReconstructionScheme::kPiecewiseLinear;
  hydro::HydroSlopeLimiter slope_limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral;
  // Periodic box extent in the patch coordinate frame (comoving), required
  // when physical_boundary_kind == kPeriodic (box [0, L) per axis). CHUI
  // Part 3: a patch spanning the box along an axis wraps in-patch; any other
  // patch face on a periodic box face is an ordinary same-level or
  // coarse/fine interface with its wrapped neighbour (local or remote; the
  // minimum-image face test and wrapped ghost probes/volumes/register keys).
  // Only a hierarchy that leaves a periodic box face without any wrapped
  // neighbour is refused (std::runtime_error).
  std::array<double, 3> periodic_box_extent_comoving{0.0, 0.0, 0.0};
  // Ghost layers per boundary face: 0 = derived from reconstruction_scheme
  // (PCM 1, PLM 2, PPM 3) so reconstruction keeps full order across patch,
  // coarse/fine and rank interfaces. A non-zero value overrides (tests).
  std::uint8_t ghost_layers = 0;
  AmrThermodynamicCachePolicy thermodynamic_cache{};
};

// Complete physical workspace for immutable step-start ghost snapshots.
// The patch count includes inactive patches; only active shells consume state.
[[nodiscard]] std::uint64_t amrPreparedGhostWorkspaceBytes(
    std::size_t patch_count, std::size_t prepared_ghost_count);

struct ProductionAmrHydroDiagnostics {
  // Numerical-method telemetry accumulated from per-patch HydroProfileEvent.
  std::uint64_t riemann_solve_count = 0;
  std::uint64_t positivity_fallback_count = 0;
  std::uint64_t stencil_degradation_count = 0;
  std::uint64_t ctu_transverse_corrected_states = 0;
  std::uint64_t ctu_transverse_degraded_states = 0;
  std::uint64_t ctu_transverse_boundary_proxy_states = 0;
  // Ghost-side CTU states corrected from a complete tangential ghost-slab
  // stencil, and the count of such complete ghost-slab axes (H01 telemetry).
  std::uint64_t ctu_transverse_ghost_corrected_states = 0;
  std::uint64_t ctu_transverse_ghost_axes = 0;
  std::uint64_t ctu_transverse_inadmissible_states = 0;
  std::uint64_t dual_energy_auxiliary_cells = 0;
  // Ghost rows whose auxiliary energy could not be resolved from a local
  // source cell (remote sources): dual-energy protection degraded there.
  std::uint64_t dual_energy_unresolved_ghosts = 0;
  // Sum over patches of HydroConservationReport::dual_energy_sync_delta
  // (non-conservative thermal correction applied by the scatter).
  double dual_energy_thermal_correction_code = 0.0;
  std::size_t patch_count = 0;
  std::size_t advanced_patch_count = 0;
  std::size_t active_cell_count = 0;
  std::size_t active_face_count = 0;
  std::size_t flux_register_entry_count = 0;
  std::size_t levels_advanced = 0;
  std::size_t max_level_advanced = 0;
  std::size_t subcycled_level_step_count = 0;
  std::size_t pending_register_created_count = 0;
  std::size_t pending_register_completed_count = 0;
  std::size_t pending_register_deferred_count = 0;
  std::size_t pending_register_applied_count = 0;
  std::size_t pending_register_rejected_count = 0;
  std::vector<std::uint32_t> substeps_by_level;
  std::uint64_t scratch_high_water_bytes = 0;
  std::uint64_t prepared_ghost_capacity_bytes = 0;
  // Actual nested geometry capacity, excluding the parent vector header.
  std::uint64_t geometry_capacity_bytes = 0;
  std::uint64_t max_patch_conserved_bytes = 0;
  std::size_t active_batch_capacity_cells = 0;
  std::size_t face_batch_capacity = 0;
  AmrHydroGhostFillDiagnostics ghost_fill;
  RefluxDiagnostics reflux;
};

struct DistributedAmrRemotePatch {
  PatchDescriptor patch;
  int owner_rank = 0;
  std::uint64_t ghost_hydro_epoch = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  std::vector<AmrHydroSparseRemoteCell> boundary_cells;
  // G1.1: physical-time label of boundary_cells. NaN (default) = the sweep's
  // state_time_code (synchronized legacy behaviour). The distributed
  // recursive driver labels current states with their level clock and
  // coarse sources temporally interpolated to the fill time with that time,
  // so a stale coarse state can never satisfy a time-consistency check.
  double state_time_code = std::numeric_limits<double>::quiet_NaN();
};

struct DistributedAmrHydroExchange {
  int local_rank = 0;
  std::uint64_t ghost_hydro_epoch = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  std::span<const DistributedAmrRemotePatch> remote_patches;
  std::vector<FluxRegisterEntry>* outbound_remote_flux_registers = nullptr;
  // G2.3: when set (synchronized distributed sweep), local register entries
  // that are INCOMPLETE because their other side lives on another rank are
  // returned here instead of being skipped; the caller merges them with the
  // inbound remote entries (mergeFluxRegisterEntriesByKey) and applies once.
  std::vector<FluxRegisterEntry>* deferred_local_flux_registers = nullptr;
};

// G2.3: combine register halves with the same register_key (local coarse
// side + any number of remote fine sides, or vice versa). Area-weighted flux
// averages, areas and face counts are summed; input order (caller-sorted:
// local first, then inbound by source rank) fixes the floating-point
// summation order, so the result is deterministic for a given decomposition.
// Throws on incompatible metadata (coarse cell, level, axis, orientation, dt).
[[nodiscard]] std::vector<FluxRegisterEntry> mergeFluxRegisterEntriesByKey(
    std::span<const FluxRegisterEntry> entries);

struct ProductionAmrRegridDiagnostics {
  std::size_t refined_patch_count = 0;
  std::size_t derefined_patch_count = 0;
  std::size_t created_gas_cell_count = 0;
  std::size_t retired_gas_cell_count = 0;
  double conserved_mass_before = 0.0;
  double conserved_mass_after = 0.0;
  double conserved_momentum_x_before = 0.0;
  double conserved_momentum_x_after = 0.0;
  double conserved_momentum_y_before = 0.0;
  double conserved_momentum_y_after = 0.0;
  double conserved_momentum_z_before = 0.0;
  double conserved_momentum_z_after = 0.0;
  double conserved_total_energy_before = 0.0;
  double conserved_total_energy_after = 0.0;
};

[[nodiscard]] bool hasProductionAmrHydroCoverage(const core::SimulationState& state);

// Ghost layers the production sweep builds per boundary face (see
// ProductionAmrHydroOptions::ghost_layers). The distributed exchange planner
// must request the same depth (parallel::AmrInterfaceTopologyOptions).
[[nodiscard]] std::uint8_t productionAmrHydroGhostLayers(const ProductionAmrHydroOptions& options);

[[nodiscard]] bool patchStateRowHasExplicitGeometry(
    const core::SimulationState& state,
    std::size_t patch_index);

[[nodiscard]] PatchDescriptor patchDescriptorFromStateRow(
    const core::SimulationState& state,
    std::size_t patch_index);

void writePatchDescriptorToStateRow(
    core::SimulationState& state,
    std::size_t patch_index,
    const PatchDescriptor& patch);

[[nodiscard]] std::vector<PatchDescriptor> buildProductionAmrPatchDescriptors(
    const core::SimulationState& state);

void populateAmrHydroFluxRegisterFaces(
    AmrHydroPatchGeometry& patch_geometry,
    std::span<const PatchDescriptor> all_patches);

// Periodic-aware overload (CHUI Part 3): a non-zero extent enables wrapped
// source-patch resolution for coarse/fine register keys across the box face.
void populateAmrHydroFluxRegisterFaces(
    AmrHydroPatchGeometry& patch_geometry,
    std::span<const PatchDescriptor> all_patches,
    const std::array<double, 3>& periodic_box_extent_comoving);

void scatterAmrHydroConservedState(
    core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroConservedStateSoa& conserved,
    double adiabatic_index);

[[nodiscard]] RefluxDiagnostics applyFluxRegistersToSimulationState(
    core::SimulationState& state,
    std::span<const FluxRegisterEntry> entries,
    std::span<const PatchDescriptor> all_patches,
    double adiabatic_index,
    const AmrThermodynamicCachePolicy& thermodynamic_cache = {});

[[nodiscard]] std::size_t mergeFluxRegistersIntoPendingStore(
    core::SimulationState& state,
    std::span<const FluxRegisterEntry> entries,
    const ProductionAmrHydroOptions& options);

[[nodiscard]] RefluxDiagnostics applyCompletePendingFluxRegistersToSimulationState(
    core::SimulationState& state,
    std::span<const PatchDescriptor> all_patches,
    double adiabatic_index,
    const AmrThermodynamicCachePolicy& thermodynamic_cache = {},
    core::MemoryGovernor* governor = nullptr);

// Gate-2 active-level row storage contract. coarse_rows and fine_rows coexist,
// so admission must include BOTH physical capacities. The workspace owns both
// buffers; capacity is admitted via the existing regrid governor before growth
// and reconciled to real capacity() when retained across calls.
struct AmrActiveLevelRowWorkspace {
  std::vector<std::uint32_t> coarse_rows;
  std::vector<std::uint32_t> fine_rows;
  core::MemoryReservation reservation;
  std::uint64_t admitted_bytes = 0U;
};

[[nodiscard]] std::uint64_t activeLevelRowStagingBytes(
    std::size_t coarse_capacity,
    std::size_t fine_capacity);
void fillActiveRowsForLevelInto(
    const core::SimulationState& state,
    int level,
    std::span<const std::uint32_t> requested_rows,
    std::vector<std::uint32_t>& out_rows);
void prepareAmrActiveLevelRows(
    const core::SimulationState& state,
    std::span<const std::uint32_t> requested_rows,
    int min_level,
    int max_level,
    core::MemoryGovernor* governor,
    AmrActiveLevelRowWorkspace& workspace);

[[nodiscard]] ProductionAmrHydroDiagnostics advanceProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrHydroDiagnostics advanceDistributedProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options,
    const DistributedAmrHydroExchange& exchange);

// G1.1 — collective services for the DISTRIBUTED recursive driver. Every
// method is collective over the world communicator and is entered by every
// rank the same number of times in the same order, including ranks that hold
// no cell on the level being advanced. Implementations live in the workflow
// layer (MPI); the AMR library stays MPI-free.
class RecursiveAmrDistributedHooks {
 public:
  virtual ~RecursiveAmrDistributedHooks() = default;
  // Rethrows on every rank if any rank holds a failure.
  virtual void reduceFailure(std::exception_ptr local_failure, const char* stage) = 0;
  // Bitwise OR over ranks (populated-level set, 256 levels).
  virtual void globalLevelUnion(std::array<std::uint64_t, 4>& level_bits) = 0;
  // Bitwise OR / max of a scalar flag over ranks.
  [[nodiscard]] virtual std::uint64_t globalMax(std::uint64_t value) = 0;
  // Synchronized start: one all-level boundary exchange; returns the remote
  // neighbour patch descriptors (cross-rank 2:1 validation, regrid planner).
  [[nodiscard]] virtual std::span<const PatchDescriptor> prepare(double t_begin) = 0;
  // Before the level-`level` sweep at t_begin: exchanges the CURRENT boundary
  // states of levels {level, finer_level} (finer_level < 0: none). With
  // capture_start the remote level-`level` cells are retained as the start
  // of the level interval [t_begin, t_end]. Returns the remote sources for
  // the sweep: current same/finer-level patches (state_time_code = t_begin)
  // plus coarser_level patches linearly interpolated in conserved variables
  // to t_begin from their retained interval (coarser_level < 0: none).
  [[nodiscard]] virtual std::span<const DistributedAmrRemotePatch> beginLevelSweep(
      int level, int coarser_level, int finer_level, double t_begin, double t_end, bool capture_start) = 0;
  // After the sweep: routes outbound remote register halves to the coarse
  // owners and merges inbound halves into the canonical pending store with
  // the SAME level options (identical interval tags on every rank).
  virtual void exchangeLevelFluxRegisters(core::SimulationState& state,
                                          std::span<const FluxRegisterEntry> outbound,
                                          const ProductionAmrHydroOptions& level_options) = 0;
  // After the level-`level` sweep (finer level exists): retains the END
  // states of remote level-`level` cells (exchange of {level, finer_level}).
  virtual void captureLevelEnd(int level, int finer_level, double t_end) = 0;
  // Drops the retained remote interval of `level` (local operation).
  virtual void retireLevel(int level) noexcept = 0;
};

// G1.1: distributed recursive Berger-Oliger subcycling. Same time-ownership
// and reflux contract as advanceProductionAmrHydroRecursive, with:
//  * the GLOBAL populated-level set (ranks with no cells on a level still
//    traverse it and enter every collective);
//  * 2:1 face balance validated between local patches and the remote
//    neighbour descriptors (every cross-rank pair is checked by both sides);
//  * same-level/finer ghosts exchanged at each level clock, coarser remote
//    ghosts temporally interpolated from exchanged start/end states;
//  * remote fine-side flux registers merged into the coarse owner's pending
//    store after every sweep and applied once the level interval completes;
//  * every local failure reduced before the next collective.
[[nodiscard]] ProductionAmrHydroDiagnostics advanceDistributedProductionAmrHydroRecursive(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& coarse_update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options,
    int local_rank,
    std::uint64_t ghost_hydro_epoch,
    RecursiveAmrDistributedHooks& hooks);

// CHUI Part 3 Milestone I: recursive Berger-Oliger subcycling over every
// level present (>= 2, no upper limit). Level L advances r^(L - L_min) times
// per coarse step with dt_L = dt / r^(L - L_min); coarse-to-fine ghosts use
// the parent level's nested temporal history; fine-to-coarse ghosts use the
// synchronized fine state; the (L, L+1) reflux is applied once level L+1 has
// completed all substeps of the level-L interval. Single-rank (local
// hierarchy); the distributed form is advanceDistributedProductionAmrHydroRecursive.
[[nodiscard]] ProductionAmrHydroDiagnostics advanceProductionAmrHydroRecursive(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& coarse_update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrHydroDiagnostics advanceProductionAmrHydroSubcycled(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& coarse_update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t first_child_patch_id,
    std::uint64_t first_child_gas_cell_id,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t replacement_gas_cell_id,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options = {});

}  // namespace cosmosim::amr
