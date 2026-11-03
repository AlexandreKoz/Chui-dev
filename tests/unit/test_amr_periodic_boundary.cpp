// AMR periodic physical boundary and ghost-layer stencil regression tests.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Defect repaired: an AMR patch spanning a periodic box built its physical
// ghosts with source_real_cell = owner, so each ghost copied the boundary
// cell itself and the solver deposited the box-face flux back onto that same
// cell: every periodic box face was a zero-flux wall (no mass, momentum,
// energy or metal crossed it). Note that zero NET global boundary
// contribution is physically correct for a periodic box; the defect was the
// zero LOCAL flux through each box face.
//
// Invariants:
//  * periodic physical ghosts copy the wrapped cell of the same patch;
//  * with the box face evaluated once (upper copy), one donor-cell step of a
//    uniform-velocity contact with metallicity moves mass and metals from the
//    last cell into the first exactly like an interior face, and globally
//    conserves mass, momentum, total energy and metal mass;
//  * stencil completion: PLM/PPM stencils use the existing ghost rows (open
//    boundary) and wrap in-patch for periodic boundaries.
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "../support/test_require.hpp"

namespace {

namespace hy = cosmosim::hydro;
constexpr double k_gamma = 5.0 / 3.0;
constexpr std::uint32_t k_nx = 4U;

// One root patch spanning the unit periodic box along x (4 x 1 x 1 cells).
cosmosim::core::SimulationState makeSpanningPatchState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(k_nx);
  state.resizePatches(1);
  state.patches.patch_id[0] = 7;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = k_nx;
  state.patches.owning_rank[0] = 0;
  state.patches.parent_patch_id[0] = 0;
  state.patches.morton_key[0] = 7;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 1.0;
  state.patches.extent_y_comoving[0] = 0.25;
  state.patches.extent_z_comoving[0] = 0.25;
  state.patches.cell_dim_x[0] = k_nx;
  state.patches.cell_dim_y[0] = 1;
  state.patches.cell_dim_z[0] = 1;
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  for (std::uint32_t row = 0; row < k_nx; ++row) {
    const double rho = 1.0 + 0.25 * static_cast<double>(row);
    state.cells.patch_index[row] = 0;
    state.cells.center_x_comoving[row] = (static_cast<double>(row) + 0.5) * 0.25;
    state.cells.center_y_comoving[row] = 0.125;
    state.cells.center_z_comoving[row] = 0.125;
    state.cells.mass_code[row] = rho * 0.25 * 0.25 * 0.25;
    state.gas_cells.gas_cell_id[row] = 500U + row;
    state.gas_cells.density_code[row] = rho;
    state.gas_cells.internal_energy_code[row] = 1.0 / ((k_gamma - 1.0) * rho);  // P = 1
    state.gas_cells.pressure_code[row] = 1.0;
    state.gas_cells.velocity_x_peculiar[row] = 0.5;
    state.gas_cells.metal_mass_code[row] = (0.01 + 0.01 * row) * state.cells.mass_code[row];
    records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = 500U + row, .parent_particle_id = std::nullopt, .owning_patch_id = 7, .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

cosmosim::amr::PatchDescriptor spanningDescriptor() {
  cosmosim::amr::PatchDescriptor patch;
  patch.patch_id = 7;
  patch.origin_comov = {0.0, 0.0, 0.0};
  patch.extent_comov = {1.0, 0.25, 0.25};
  patch.cell_dims = {static_cast<std::uint16_t>(k_nx), 1, 1};
  return patch;
}

void testPeriodicGhostsWrapAndFluxCrossesTheBoxFace() {
  const auto state = makeSpanningPatchState();
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.physical_boundary_kind = hy::HydroBoundaryKind::kPeriodic;
  const auto patch_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(state, spanningDescriptor(), options);
  const hy::HydroPatchGeometry& geometry = patch_geometry.geometry;
  for (const hy::HydroGhostCell& ghost : geometry.ghost_cells) {
    if (ghost.axis != hy::HydroFaceAxis::kX) {
      continue;
    }
    const std::size_t expected = ghost.side == hy::HydroFaceSide::kLower ? k_nx - 1U : 0U;
    CHUI_REQUIRE(ghost.source_real_cell == expected);
    CHUI_REQUIRE(ghost.boundary_kind == hy::HydroBoundaryKind::kPeriodic);
  }
  hy::HydroConservedStateSoa conserved = cosmosim::amr::loadAmrHydroConservedState(state, patch_geometry, k_gamma);
  hy::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);
  std::vector<std::size_t> active_cells;
  for (std::size_t c = 0; c < geometry.cellCount(); ++c) active_cells.push_back(c);
  std::vector<std::size_t> active_faces;
  for (std::size_t f = 0; f < geometry.faces.size(); ++f) {
    const hy::HydroFace& face = geometry.faces[f];
    if (face.ghost_cell_slot != hy::k_invalid_ghost_cell_slot &&
        geometry.ghost_cells[face.ghost_cell_slot].side == hy::HydroFaceSide::kLower) {
      continue;
    }
    active_faces.push_back(f);
  }
  const std::vector<std::size_t> all_cells = active_cells;
  const auto before = hy::HydroCoreSolver::conservationTotals(conserved, geometry, all_cells);
  const hy::HydroConservedState first_before = conserved.loadCell(0U);
  const hy::HydroConservedState last_before = conserved.loadCell(k_nx - 1U);
  const hy::HydroConservedState second_before = conserved.loadCell(1U);
  const double dt = 0.1;
  const hy::HydroUpdateContext update{.dt_code = dt, .scale_factor = 1.0};
  const hy::HydroSourceContext source_context{.update = update};
  const hy::HydroCoreSolver solver(k_gamma);
  const hy::PiecewiseConstantReconstruction pcm;
  hy::HllcRiemannSolver riemann;
  hy::HydroScratchBuffers scratch;
  hy::HydroProfileEvent profile;
  solver.advancePatchActiveSetWithScratch(
      conserved, geometry, hy::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update, pcm, riemann, {}, source_context, scratch, nullptr, &profile);
  const auto after = hy::HydroCoreSolver::conservationTotals(conserved, geometry, all_cells);
  const double nu = 0.5 * dt / 0.25;
  // Uniform u, P: HLLC is exact for the contact -> donor cell. Cell 0
  // receives from cell n-1 across the periodic box face.
  CHUI_REQUIRE(std::abs(conserved.loadCell(0U).mass_density_comoving -
                        (first_before.mass_density_comoving -
                         nu * (first_before.mass_density_comoving - last_before.mass_density_comoving))) < 1.0e-12);
  CHUI_REQUIRE(std::abs(conserved.loadCell(0U).metal_mass_density_comoving -
                        (first_before.metal_mass_density_comoving -
                         nu * (first_before.metal_mass_density_comoving - last_before.metal_mass_density_comoving))) <
               1.0e-12);
  CHUI_REQUIRE(std::abs(conserved.loadCell(1U).mass_density_comoving -
                        (second_before.mass_density_comoving -
                         nu * (second_before.mass_density_comoving - first_before.mass_density_comoving))) < 1.0e-12);
  CHUI_REQUIRE(std::abs(after.mass - before.mass) < 1.0e-14);
  CHUI_REQUIRE(std::abs(after.momentum_x - before.momentum_x) < 1.0e-14);
  CHUI_REQUIRE(std::abs(after.total_energy - before.total_energy) < 1.0e-13);
  CHUI_REQUIRE(std::abs(after.metal_mass - before.metal_mass) < 1.0e-15);
}

void testGhostLayerStencilCompletionOpenBoundary() {
  hy::HydroPatchGeometry geometry = hy::makeCartesianPatchGeometry(hy::HydroCartesianPatchSpec{
      .nx = 4, .ny = 1, .nz = 1,
      .cell_width_x_comoving = 1.0, .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  hy::appendCartesianBoundaryGhostFaces(geometry, hy::HydroBoundaryKind::kOpen);
  CHUI_REQUIRE(hy::completeReconstructionStencils(geometry) > 0U);
  const std::size_t real = geometry.cellCount();
  for (const hy::HydroFace& face : geometry.faces) {
    if (face.axis != hy::HydroFaceAxis::kX) continue;
    if (face.ghost_cell_slot == hy::k_invalid_ghost_cell_slot && face.owner_cell == 0U) {
      // Interior face 0|1: owner_minus is the lower ghost row of cell 0.
      CHUI_REQUIRE(face.owner_minus_cell >= real);
      CHUI_REQUIRE(geometry.ghost_cells[face.owner_minus_cell - real].owner_real_cell == 0U);
      CHUI_REQUIRE(face.owner_minus2_cell == hy::k_invalid_cell_index);  // second layer unavailable
    }
    if (face.ghost_cell_slot != hy::k_invalid_ghost_cell_slot &&
        geometry.ghost_cells[face.ghost_cell_slot].side == hy::HydroFaceSide::kUpper) {
      // Upper boundary face 3|ghost: interior stencil on the owner side.
      CHUI_REQUIRE(face.owner_minus_cell == 2U && face.owner_minus2_cell == 1U);
      CHUI_REQUIRE(face.neighbor_plus_cell == hy::k_invalid_cell_index);
    }
  }
}

}  // namespace

int main() {
  testPeriodicGhostsWrapAndFluxCrossesTheBoxFace();
  testGhostLayerStencilCompletionOpenBoundary();
  return 0;
}
