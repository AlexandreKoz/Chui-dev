#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::amr {

enum class AmrGhostSourceResidency : std::uint8_t {
  kLocal,
  kRemoteReadOnly,
};

struct AmrHydroGhostFillPatch {
  AmrHydroPatchGeometry* geometry = nullptr;
  hydro::HydroConservedStateSoa* conserved = nullptr;
  AmrGhostSourceResidency residency = AmrGhostSourceResidency::kLocal;
  std::uint64_t ghost_hydro_epoch = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  // The finite-volume solver consumes ghost states before MUSCL-Hancock
  // reconstruction.  ghost_fill_time_code is therefore the step-start time
  // for the patch update, not an implicit predictor half-time.
  double target_state_time_code = 0.0;
  double ghost_fill_time_code = 0.0;
  double source_current_state_time_code = 0.0;
  const core::AmrTemporalBoundaryHistoryStore* temporal_boundary_history = nullptr;
  // Legacy dense-remote coverage mask retained only for compatibility with
  // validation callers that still provide a fully materialized remote patch.
  // Production distributed hydro uses AmrHydroSparseGhostSource instead.
  std::span<const std::uint8_t> available_real_cells;
  bool enable_temporal_coarse_to_fine = false;
  bool requires_ghost_fill = true;
  // CHUI Part 3: periodic box extent of the patch coordinate frame. All zero
  // (default) = non-periodic. When set, ghost probe points and fine-to-coarse
  // ghost volumes of THIS target are canonicalized into [0, L) so neighbours
  // across the periodic box boundary resolve to their real (local or remote)
  // patch instead of being treated as missing.
  std::array<double, 3> periodic_box_extent_comoving{0.0, 0.0, 0.0};
};

struct AmrHydroSparseRemoteCell {
  std::uint32_t patch_local_cell = 0U;
  std::uint64_t gas_cell_id = 0U;
  hydro::HydroConservedState conserved{};
};

struct AmrHydroSparseGhostSource {
  const PatchDescriptor* patch = nullptr;
  int owner_rank = 0;
  std::uint64_t ghost_hydro_epoch = 0U;
  std::uint64_t expected_ghost_hydro_epoch = 0U;
  double source_current_state_time_code = 0.0;
  std::span<const AmrHydroSparseRemoteCell> cells;
};

struct AmrHydroGhostFillDiagnostics {
  std::size_t physical_ghosts_filled = 0;
  std::size_t same_level_ghosts_filled = 0;
  std::size_t coarse_to_fine_ghosts_filled = 0;
  std::size_t fine_to_coarse_ghosts_filled = 0;
  std::size_t temporal_coarse_to_fine_ghosts_filled = 0;
  std::size_t temporal_endpoint_ghosts_filled = 0;
  std::size_t skipped_remote_ghosts = 0;
  std::size_t stale_epoch_rejections = 0;
  std::size_t temporal_same_level_mismatch_rejections = 0;
  std::size_t temporal_history_missing_rejections = 0;
  std::size_t temporal_history_invalid_rejections = 0;
  std::size_t temporal_time_out_of_range_rejections = 0;
  std::size_t temporal_fine_to_coarse_misalignment_rejections = 0;
  std::size_t temporal_geometry_mismatch_rejections = 0;
  std::size_t temporal_identity_mismatch_rejections = 0;
  std::size_t missing_source_records = 0;
  std::size_t unresolved_ghosts = 0;
  // Deep stencil-only layers continued by zero gradient (domain edge behind a
  // thin neighbour patch): counted reconstruction-order reductions.
  std::size_t deep_layer_extrapolated_ghosts = 0;
};

[[nodiscard]] std::uint64_t amrPatchGeometryFingerprint(const PatchDescriptor& patch);

// Capture the start/end conserved states for coarse patches participating in a
// local AMR subcycle.  The history is persistent state: callers may checkpoint
// after the coarse advance and resume fine substeps using the same records.
void captureAmrTemporalBoundaryHistoryStart(
    core::SimulationState& state,
    std::span<const PatchDescriptor> patches,
    double interval_start_code,
    double adiabatic_index, core::MemoryGovernor* governor = nullptr);
void captureAmrTemporalBoundaryHistoryEnd(
    core::SimulationState& state,
    std::span<const PatchDescriptor> patches,
    double interval_end_code,
    double adiabatic_index, core::MemoryGovernor* governor = nullptr);
void retireAmrTemporalBoundaryHistory(core::SimulationState& state, core::MemoryGovernor* governor = nullptr);

// CHUI Part 3 Milestone I: nested per-level intervals for recursive
// subcycling (>= 3 levels). StartForLevel captures only level-`level` coarse
// sources and refuses if that level already has an open interval; records of
// other levels are preserved. EndForLevel completes only that level (level
// -1 = all, the historical behaviour); RetireLevel removes only that level.
// Start/end build admitted candidates and publish once. Any allocation,
// identity, geometry, state or admission failure preserves every old record.
// Retirement allocates no record storage; pass the same governor used at
// capture to reconcile released nested capacities. Defaults preserve callers.
void captureAmrTemporalBoundaryHistoryStartForLevel(
    core::SimulationState& state, std::span<const PatchDescriptor> patches, int level,
    double interval_start_code, double adiabatic_index, core::MemoryGovernor* governor = nullptr);
void captureAmrTemporalBoundaryHistoryEndForLevel(
    core::SimulationState& state, std::span<const PatchDescriptor> patches, int level,
    double interval_end_code, double adiabatic_index, core::MemoryGovernor* governor = nullptr);
void retireAmrTemporalBoundaryHistoryLevel(
    core::SimulationState& state, int level, core::MemoryGovernor* governor = nullptr);


// A single reusable selection is sufficient: at most one record per source
// cell can contribute to any ghost. No theoretical refinement-volume bitmap
// or per-candidate cell-index vector is materialized.
[[nodiscard]] std::uint64_t amrGhostFillScratchBytes(
    std::span<const AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources);

[[nodiscard]] AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index, core::MemoryGovernor* governor);

[[nodiscard]] AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    double adiabatic_index);

[[nodiscard]] AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index);

}  // namespace cosmosim::amr
