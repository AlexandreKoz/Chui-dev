#pragma once

#include <cstddef>
#include <cstdint>
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
  // Spatial reconstruction and limiter used for every patch advance. The
  // time-integration family (MUSCL-Hancock / CTU) and dual-energy policy are
  // owned by the HydroCoreSolver passed to the advance entry points. PPM and
  // CTU use the patch-interior stencils; within two cells of a patch boundary
  // PPM reduces to PLM and CTU transverse terms from non-periodic ghosts are
  // omitted (both counted in HydroProfileEvent / diagnostics). Multi-layer
  // AMR ghost shells are Part 3 scope.
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
  std::uint64_t ctu_transverse_degraded_states = 0;
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
};

struct DistributedAmrHydroExchange {
  int local_rank = 0;
  std::uint64_t ghost_hydro_epoch = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  std::span<const DistributedAmrRemotePatch> remote_patches;
  std::vector<FluxRegisterEntry>* outbound_remote_flux_registers = nullptr;
};

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
    const AmrThermodynamicCachePolicy& thermodynamic_cache = {});

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
