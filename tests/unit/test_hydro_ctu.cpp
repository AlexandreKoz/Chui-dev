// Corner-transport-upwind (CTU) integration and periodic-boundary tests.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Reference: for linear advection with constant (u, v) > 0 and first-order
// (piecewise-constant) states, Colella's (1990) donor-cell CTU update is
//   q^{n+1}_{ij} = (1-nx)(1-ny) q_ij + nx(1-ny) q_{i-1,j}
//                  + ny(1-nx) q_{i,j-1} + nx ny q_{i-1,j-1},
// with nx = u dt/(a dx), ny = v dt/(a dy). A density contact advected by the
// Euler equations with uniform velocity and pressure is exactly this problem,
// and HLLC resolves isolated contacts exactly, so the CHUI CTU path must
// reproduce the formula to round-off. The directionally unsplit donor-cell
// (MUSCL-Hancock family, no transverse terms) lacks the corner term and must
// differ. The test is repeated at a=2 with dt doubled (same nx, ny) to check
// that the transverse increments and the flux update share the 1/a factor.
//
// Also covered: Riemann evaluation count (2 per face for CTU), conservation
// residual, free-stream preservation in 3D with PLM+CTU, and the periodic
// boundary rule that each box face is evaluated exactly once.
#include <cmath>
#include <cstddef>
#include <vector>

#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "../support/test_require.hpp"

namespace {

namespace hy = cosmosim::hydro;

constexpr double k_gamma = 5.0 / 3.0;
constexpr std::size_t k_n = 4U;
constexpr double k_velocity = 0.3;

double densityAt(std::size_t i, std::size_t j) {
  return 1.0 + 0.3 * static_cast<double>((i * 7U + j * 3U) % 5U) / 4.0;
}

struct PeriodicBox {
  hy::HydroPatchGeometry geometry;
  hy::HydroConservedStateSoa conserved;
  std::vector<std::size_t> active_cells;
  std::vector<std::size_t> active_faces;
};

PeriodicBox makeBox(std::size_t nz, bool include_lower_periodic_duplicates) {
  PeriodicBox box;
  box.geometry = hy::makeCartesianPatchGeometry(hy::HydroCartesianPatchSpec{
      .nx = k_n, .ny = k_n, .nz = nz,
      .cell_width_x_comoving = 1.0, .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  hy::appendCartesianBoundaryGhostFaces(box.geometry, hy::HydroBoundaryKind::kPeriodic);
  (void)hy::completePeriodicReconstructionStencils(box.geometry);
  box.conserved.resize(box.geometry.totalCellStorageCount());
  for (std::size_t k = 0; k < nz; ++k) {
    for (std::size_t j = 0; j < k_n; ++j) {
      for (std::size_t i = 0; i < k_n; ++i) {
        hy::HydroPrimitiveState primitive;
        primitive.rho_comoving = densityAt(i, j);
        primitive.vel_x_peculiar = k_velocity;
        primitive.vel_y_peculiar = k_velocity;
        primitive.pressure_comoving = 1.0;
        box.conserved.storeCell(box.geometry.linearCellIndex(i, j, k),
            hy::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
      }
    }
  }
  hy::fillHydroBoundaryGhostCells(box.conserved, box.geometry, k_gamma);
  for (std::size_t cell = 0; cell < box.geometry.cellCount(); ++cell) {
    box.active_cells.push_back(cell);
  }
  for (std::size_t face = 0; face < box.geometry.faces.size(); ++face) {
    const hy::HydroFace& f = box.geometry.faces[face];
    if (f.ghost_cell_slot != hy::k_invalid_ghost_cell_slot && !include_lower_periodic_duplicates &&
        box.geometry.ghost_cells[f.ghost_cell_slot].side == hy::HydroFaceSide::kLower) {
      continue;
    }
    box.active_faces.push_back(face);
  }
  return box;
}

hy::HydroProfileEvent advance(PeriodicBox& box, const hy::HydroCoreSolver& solver,
                              const hy::HydroReconstruction& reconstruction, double dt, double a) {
  const hy::HydroUpdateContext update{.dt_code = dt, .scale_factor = a};
  const hy::HydroSourceContext source_context{.update = update};
  hy::HllcRiemannSolver riemann;
  hy::HydroScratchBuffers scratch;
  hy::HydroProfileEvent profile;
  solver.advancePatchActiveSetWithScratch(
      box.conserved, box.geometry,
      hy::HydroActiveSetView{.active_cells = box.active_cells, .active_faces = box.active_faces},
      update, reconstruction, riemann, {}, source_context, scratch, nullptr, &profile);
  return profile;
}

double ctuReference(std::size_t i, std::size_t j, double nu) {
  const std::size_t im = (i + k_n - 1U) % k_n;
  const std::size_t jm = (j + k_n - 1U) % k_n;
  return (1.0 - nu) * (1.0 - nu) * densityAt(i, j) + nu * (1.0 - nu) * densityAt(im, j) +
      nu * (1.0 - nu) * densityAt(i, jm) + nu * nu * densityAt(im, jm);
}

double unsplitDonorReference(std::size_t i, std::size_t j, double nu) {
  const std::size_t im = (i + k_n - 1U) % k_n;
  const std::size_t jm = (j + k_n - 1U) % k_n;
  return densityAt(i, j) - nu * (densityAt(i, j) - densityAt(im, j)) - nu * (densityAt(i, j) - densityAt(i, jm));
}

void testDonorCellCtuMatchesExactCornerTransport(double dt, double a) {
  const double nu = k_velocity * dt / a;
  PeriodicBox box = makeBox(1U, false);
  const hy::HydroCoreSolver ctu_solver(
      k_gamma, hy::HydroNumericsPolicy{.integration_scheme = hy::HydroIntegrationScheme::kCornerTransportUpwind});
  const hy::PiecewiseConstantReconstruction pcm;
  const hy::HydroProfileEvent profile = advance(box, ctu_solver, pcm, dt, a);
  CHUI_REQUIRE(profile.riemann_solve_count == 2U * box.active_faces.size());
  CHUI_REQUIRE(profile.ctu_transverse_degraded_states == 0U);
  double max_unsplit_difference = 0.0;
  for (std::size_t j = 0; j < k_n; ++j) {
    for (std::size_t i = 0; i < k_n; ++i) {
      const double rho = box.conserved.loadCell(box.geometry.linearCellIndex(i, j, 0U)).mass_density_comoving;
      CHUI_REQUIRE(std::abs(rho - ctuReference(i, j, nu)) < 1.0e-11);
      max_unsplit_difference = std::max(max_unsplit_difference, std::abs(rho - unsplitDonorReference(i, j, nu)));
    }
  }
  // The corner-coupling term nu^2 (q_ij - q_{i-1,j} - q_{i,j-1} + q_{i-1,j-1})
  // is resolved: CTU is not the directionally unsplit donor-cell scheme.
  CHUI_REQUIRE(max_unsplit_difference > 1.0e-4);
  // Conservative bookkeeping is exact.
  CHUI_REQUIRE(std::abs(profile.conservation.residual.mass) < 1.0e-12);
  CHUI_REQUIRE(std::abs(profile.conservation.residual.total_energy) < 1.0e-11);
}

void testMusclHancockFamilyIsUnsplitDonorCell() {
  const double dt = 0.5;
  PeriodicBox box = makeBox(1U, false);
  const hy::HydroCoreSolver solver(k_gamma);
  const hy::PiecewiseConstantReconstruction pcm;
  const hy::HydroProfileEvent profile = advance(box, solver, pcm, dt, 1.0);
  CHUI_REQUIRE(profile.riemann_solve_count == box.active_faces.size());
  for (std::size_t j = 0; j < k_n; ++j) {
    for (std::size_t i = 0; i < k_n; ++i) {
      const double rho = box.conserved.loadCell(box.geometry.linearCellIndex(i, j, 0U)).mass_density_comoving;
      CHUI_REQUIRE(std::abs(rho - unsplitDonorReference(i, j, k_velocity * dt)) < 1.0e-11);
    }
  }
}

void testPeriodicBoxFaceEvaluatedOnce() {
  // Including the lower-side periodic ghost faces as well as the upper-side
  // ones applies every box-boundary flux twice; the corrected rule (upper
  // copies only, as the fixed-grid runtime now selects) reproduces the
  // interior donor-cell update in the boundary cells too.
  const double dt = 0.5;
  PeriodicBox duplicated = makeBox(1U, true);
  const hy::HydroCoreSolver solver(k_gamma);
  const hy::PiecewiseConstantReconstruction pcm;
  (void)advance(duplicated, solver, pcm, dt, 1.0);
  const double rho_boundary_duplicated =
      duplicated.conserved.loadCell(duplicated.geometry.linearCellIndex(0U, 1U, 0U)).mass_density_comoving;
  CHUI_REQUIRE(std::abs(rho_boundary_duplicated - unsplitDonorReference(0U, 1U, k_velocity * dt)) > 1.0e-6);
}

void testFreeStreamPreservationWithPlmCtuIn3d() {
  PeriodicBox box = makeBox(k_n, false);
  // Overwrite with a uniform state.
  for (std::size_t cell = 0; cell < box.geometry.cellCount(); ++cell) {
    hy::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.3;
    primitive.vel_x_peculiar = 0.2;
    primitive.vel_y_peculiar = -0.1;
    primitive.vel_z_peculiar = 0.05;
    primitive.pressure_comoving = 0.7;
    box.conserved.storeCell(cell, hy::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }
  hy::fillHydroBoundaryGhostCells(box.conserved, box.geometry, k_gamma);
  const hy::HydroConservedState reference = box.conserved.loadCell(0U);
  const hy::HydroCoreSolver solver(
      k_gamma, hy::HydroNumericsPolicy{.integration_scheme = hy::HydroIntegrationScheme::kCornerTransportUpwind});
  const auto courant = hy::hydroTransportCourantFactors(0.3, 1.0, 1.0, 1.0, 1.0);
  const hy::MusclHancockReconstruction plm(hy::HydroReconstructionPolicy{
      .dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
      .rho_floor = 1.0e-12, .pressure_floor = 1.0e-12, .adiabatic_index = k_gamma});
  const hy::HydroProfileEvent profile = advance(box, solver, plm, 0.3, 1.0);
  for (std::size_t cell = 0; cell < box.geometry.cellCount(); ++cell) {
    const hy::HydroConservedState updated = box.conserved.loadCell(cell);
    CHUI_REQUIRE(std::abs(updated.mass_density_comoving - reference.mass_density_comoving) < 1.0e-13);
    CHUI_REQUIRE(std::abs(updated.total_energy_density_comoving - reference.total_energy_density_comoving) < 1.0e-13);
  }
  CHUI_REQUIRE(profile.ctu_transverse_corrected_states == 2U * box.active_faces.size());
}

void testImportedGhostUsesCountedOneSidedCtuClosure() {
  hy::HydroPatchGeometry geometry = hy::makeCartesianPatchGeometry(
      hy::HydroCartesianPatchSpec{
          .nx = 3U, .ny = 3U, .nz = 1U,
          .cell_width_x_comoving = 1.0, .cell_width_y_comoving = 1.0,
          .cell_width_z_comoving = 1.0});
  hy::appendCartesianBoundaryGhostFaces(geometry, hy::HydroBoundaryKind::kOpen);
  (void)hy::completeReconstructionStencils(geometry);
  hy::HydroConservedStateSoa conserved;
  conserved.resize(geometry.totalCellStorageCount());
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    hy::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.0 + 0.01 * static_cast<double>(cell);
    primitive.vel_x_peculiar = 0.15;
    primitive.vel_y_peculiar = -0.08;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(cell, hy::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }
  for (hy::HydroGhostCell& ghost : geometry.ghost_cells) {
    ghost.boundary_kind = hy::HydroBoundaryKind::kImportedMpi;
    ghost.mutation_rights = hy::HydroGhostMutationRights::kReadOnlyImported;
    CHUI_REQUIRE(ghost.owner_real_cell < geometry.cellCount());
    conserved.storeCell(ghost.ghost_cell, conserved.loadCell(ghost.owner_real_cell));
  }
  std::vector<std::size_t> active_cells(geometry.cellCount());
  for (std::size_t i = 0; i < active_cells.size(); ++i) active_cells[i] = i;
  std::vector<std::size_t> active_faces(geometry.faces.size());
  for (std::size_t i = 0; i < active_faces.size(); ++i) active_faces[i] = i;
  const hy::HydroCoreSolver solver(
      k_gamma, hy::HydroNumericsPolicy{
          .integration_scheme = hy::HydroIntegrationScheme::kCornerTransportUpwind});
  const auto courant = hy::hydroTransportCourantFactors(0.05, 1.0, 1.0, 1.0, 1.0);
  const hy::MusclHancockReconstruction reconstruction(hy::HydroReconstructionPolicy{
      .dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
      .rho_floor = 1.0e-12, .pressure_floor = 1.0e-12,
      .adiabatic_index = k_gamma});
  hy::HllcRiemannSolver riemann;
  hy::HydroScratchBuffers scratch;
  hy::HydroProfileEvent profile;
  const hy::HydroUpdateContext update{.dt_code = 0.05, .scale_factor = 1.0};
  const hy::HydroSourceContext source{.update = update};
  solver.advancePatchActiveSetWithScratch(
      conserved, geometry,
      hy::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update, reconstruction, riemann, {}, source, scratch, nullptr, &profile);
  CHUI_REQUIRE(profile.ctu_transverse_boundary_proxy_states > 0U);
}

void testPeriodicStencilCompletion() {
  hy::HydroPatchGeometry geometry = hy::makeCartesianPatchGeometry(hy::HydroCartesianPatchSpec{
      .nx = 5, .ny = 1, .nz = 1,
      .cell_width_x_comoving = 1.0, .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  hy::appendCartesianBoundaryGhostFaces(geometry, hy::HydroBoundaryKind::kPeriodic);
  CHUI_REQUIRE(hy::completePeriodicReconstructionStencils(geometry) > 0U);
  for (const hy::HydroFace& face : geometry.faces) {
    if (face.axis != hy::HydroFaceAxis::kX) {
      continue;
    }
    if (face.ghost_cell_slot == hy::k_invalid_ghost_cell_slot && face.owner_cell == 0U) {
      // Interior face 0|1: owner_minus wraps to cell 4, owner_minus2 to 3.
      CHUI_REQUIRE(face.owner_minus_cell == 4U && face.owner_minus2_cell == 3U);
      CHUI_REQUIRE(face.neighbor_plus_cell == 2U && face.neighbor_plus2_cell == 3U);
    }
    if (face.ghost_cell_slot != hy::k_invalid_ghost_cell_slot &&
        geometry.ghost_cells[face.ghost_cell_slot].side == hy::HydroFaceSide::kUpper) {
      // Upper box face 4|0 (normal +x): stencil 2,3 | 4 | 0(ghost) -> 1,2.
      CHUI_REQUIRE(face.owner_cell == 4U);
      CHUI_REQUIRE(face.owner_minus_cell == 3U && face.owner_minus2_cell == 2U);
      CHUI_REQUIRE(face.neighbor_plus_cell == 1U && face.neighbor_plus2_cell == 2U);
    }
  }
}

}  // namespace

int main() {
  testDonorCellCtuMatchesExactCornerTransport(0.5, 1.0);
  testDonorCellCtuMatchesExactCornerTransport(1.0, 2.0);
  testMusclHancockFamilyIsUnsplitDonorCell();
  testPeriodicBoxFaceEvaluatedOnce();
  testFreeStreamPreservationWithPlmCtuIn3d();
  testImportedGhostUsesCountedOneSidedCtuClosure();
  testPeriodicStencilCompletion();
  return 0;
}
