#pragma once

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {

// Add one ghost row per physical patch boundary face. Ghost rows are transient
// hydro scratch rows addressed after the real cell range and are not
// SimulationState cells.
void appendCartesianBoundaryGhostFaces(
    HydroPatchGeometry& geometry,
    HydroBoundaryKind boundary_kind);

// Complete the along-normal reconstruction stencils (owner_minus2,
// owner_minus, neighbor_plus, neighbor_plus2) of every face of a fully
// periodic single-patch geometry by periodic wrap of real cells. Periodic
// ghost neighbors are resolved through HydroGhostCell::source_real_cell.
// Without this, faces within two cells of a periodic boundary would be
// reconstructed at reduced order although the wrapped data are available.
// Only unset (k_invalid_cell_index) stencil entries are written. Requires
// every ghost cell of the geometry to be periodic; returns the number of
// stencil entries filled.
std::size_t completePeriodicReconstructionStencils(HydroPatchGeometry& geometry);

// Fill physical-boundary ghost rows from authoritative real cells. Imported MPI
// ghosts are read-only and intentionally skipped by this helper.
void fillHydroBoundaryGhostCells(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    double adiabatic_index);

}  // namespace cosmosim::hydro
