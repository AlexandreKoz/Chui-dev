#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::amr {

enum class AmrHydroBoundaryClass : std::uint8_t {
  kPhysical,
  kSameLevel,
  kCoarseFine,
};

enum class AmrHydroGhostFillStatus : std::uint8_t {
  kUnfilledPhysicalBoundary,
  kUnfilledSameLevel,
  kUnfilledCoarseFine,
  kFilledPhysicalBoundary,
  kFilledSameLevel,
  kFilledCoarseToFine,
  kFilledCoarseToFineTemporal,
  kFilledFineToCoarse,
  kSkippedRemote,
  kRejectedStaleRemote,
  kRejectedTemporalSameLevelMismatch,
  kRejectedTemporalHistoryMissing,
  kRejectedTemporalHistoryInvalid,
  kRejectedTemporalTimeOutOfRange,
  kRejectedTemporalFineToCoarseMismatch,
  kMissingSource,
  kFilledDeepLayerExtrapolated,
  kFilledDeepLayerReflected,
};

struct AmrHydroGeometryOptions {
  hydro::HydroBoundaryKind physical_boundary_kind = hydro::HydroBoundaryKind::kOpen;
  std::array<AmrHydroBoundaryClass, 6> boundary_classes = {
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
  };
  // CHUI Part 3: number of ghost layers per boundary face. Layer 1 is the
  // face-adjacent ghost that carries the boundary face; layers 2..N are
  // stencil-only ghosts (no face, never advanced, never receive flux) that
  // give PLM (N = 2) and PPM (N = 3) full-order reconstruction across patch,
  // coarse/fine, rank and open/reflective boundaries. In-patch periodic wraps
  // never need deep layers (the stencil resolves to real wrapped cells).
  // A deep layer that cannot be formed (patch narrower than the mirror depth
  // for reflective boundaries) is omitted and the resulting order reduction
  // is counted by the solver's stencil-degradation counter.
  std::uint8_t ghost_layers = 1;
  // P1-HYDRO-02: transverse depth of the edge stencil ghost block (0 = none).
  // CTU sets it to the reconstruction stencil half-width so the tangential
  // ghost slab of every face ghost is complete across patch edges and
  // corners (same-level, coarse/fine, periodic and MPI). Edge ghosts are
  // stencil-only: no face, no flux register, never advanced.
  std::uint8_t edge_stencil_depth = 0;
};

struct AmrHydroCellDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t gas_cell_id = 0;
  std::uint32_t local_cell_row = core::kInvalidGasCellRow;
  std::size_t patch_local_cell = hydro::k_invalid_cell_index;
};

struct AmrHydroGhostDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t source_patch_id = 0;
  std::uint64_t source_gas_cell_id = 0;
  std::uint64_t owner_gas_cell_id = 0;
  std::uint32_t owner_local_cell_row = core::kInvalidGasCellRow;
  std::size_t ghost_slot = hydro::k_invalid_ghost_cell_slot;
  std::size_t ghost_cell = hydro::k_invalid_cell_index;
  AmrHydroBoundaryClass boundary_class = AmrHydroBoundaryClass::kPhysical;
  AmrHydroGhostFillStatus fill_status = AmrHydroGhostFillStatus::kUnfilledPhysicalBoundary;
  // 1 = face-adjacent ghost; >1 = stencil-only deep layer at distance
  // ghost_layer cell widths from the owner cell along the outward normal.
  std::uint8_t ghost_layer = 1;
  // P1-HYDRO-02: edge stencil ghost (see hydro::HydroGhostCell::
  // transverse_axis). boundary_class is copied from the owner's face ghost;
  // the filler resolves edge ghosts with any legal level relation and applies
  // the physical domain closure on every axis that leaves the domain.
  bool edge_stencil = false;
};

struct AmrHydroFaceDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t face_id = 0;
  std::uint64_t owner_gas_cell_id = 0;
  std::uint64_t neighbor_gas_cell_id = 0;
  std::uint64_t ghost_source_gas_cell_id = 0;
  std::size_t owner_patch_cell = hydro::k_invalid_cell_index;
  std::size_t neighbor_patch_cell = hydro::k_invalid_cell_index;
  std::size_t ghost_slot = hydro::k_invalid_ghost_cell_slot;
  double area_comoving = 0.0;
  double normal_x = 0.0;
  double normal_y = 0.0;
  double normal_z = 0.0;
  hydro::HydroFaceAxis axis = hydro::HydroFaceAxis::kX;
  AmrHydroBoundaryClass boundary_class = AmrHydroBoundaryClass::kPhysical;
};

struct AmrHydroPatchGeometry {
  PatchDescriptor patch;
  hydro::HydroPatchGeometry geometry;
  std::vector<AmrHydroCellDescriptor> real_cells;
  std::vector<std::uint64_t> gas_cell_ids;
  std::vector<std::uint32_t> local_cell_rows;
  std::vector<AmrHydroGhostDescriptor> ghosts;
  std::vector<AmrHydroFaceDescriptor> faces;
  std::uint64_t source_gas_cell_identity_generation = 0;

  // Physical capacities of the eight nested geometry arrays, excluding this
  // object's own storage in its parent vector. No logical-size substitution.
  [[nodiscard]] std::uint64_t ownedCapacityBytes() const;
  [[nodiscard]] std::span<const std::uint64_t> gasCellIds() const noexcept;
  [[nodiscard]] std::vector<std::size_t> internalFaceIndices() const;
};

struct AmrHydroGeometryCapacity {
  std::size_t real_cells = 0;
  std::size_t ghost_cells = 0;
  std::size_t faces = 0;
  std::uint64_t retained_bytes = 0;
  std::uint64_t construction_scratch_bytes = 0;
};

// Checked, exact topology counts for the regular patch representation.
// retained_bytes excludes the AmrHydroPatchGeometry object itself.
// ghost_layers >= 1 bounds the deep stencil-only layers (see
// AmrHydroGeometryOptions::ghost_layers); the default keeps the historical
// one-layer bound for existing callers.
[[nodiscard]] AmrHydroGeometryCapacity amrHydroGeometryCapacity(
    const PatchDescriptor& patch, std::uint8_t ghost_layers = 1, std::uint8_t edge_stencil_depth = 0);

[[nodiscard]] AmrHydroPatchGeometry buildAmrHydroPatchGeometry(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    const AmrHydroGeometryOptions& options = {});

[[nodiscard]] AmrHydroPatchGeometry buildRemoteAmrHydroPatchGeometry(
    const PatchDescriptor& patch,
    std::span<const std::uint64_t> gas_cell_ids,
    std::uint64_t source_gas_cell_identity_generation,
    const AmrHydroGeometryOptions& options = {},
    std::span<const std::uint8_t> available_real_cells = {});

// Refill an already admitted SoA. Production reuses one allocation across
// patches; the original returning overload remains for standalone callers.
void loadAmrHydroConservedStateInto(
    const core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    double adiabatic_index,
    hydro::HydroConservedStateSoa& conserved);

[[nodiscard]] hydro::HydroConservedStateSoa loadAmrHydroConservedState(
    const core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    double adiabatic_index);

}  // namespace cosmosim::amr
