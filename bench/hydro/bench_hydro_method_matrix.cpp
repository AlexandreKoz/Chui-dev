// Hydro method-matrix benchmark (CPU): reconstruction {PCM, PLM, PPM} x
// integration {MUSCL-Hancock, CTU} x Riemann {HLLC, HLLE}.
//
// STATUS: source only. Authored without execution; NO timing or accuracy
// claim is made by this file. Run manually (see docs/profiling.md):
//   ./build/<preset>/bench_hydro_method_matrix [cells_per_axis=32] [steps=8]
//
// Problem: oblique advection of a smooth density wave
//   rho = 1 + 0.2 sin(2 pi (x + y + z)/L), u = (u0, u0, u0), P = 1
// on a fully periodic N^3 patch. The exact solution is the translated wave,
// so the reported L1 density error at a fixed number of steps, together with
// wall time, supports "error at fixed wall time" comparisons instead of
// assuming the most expensive method is best.
//
// Reported per configuration (one line, key=value):
//   cells_per_s, faces_per_s, riemann_solves, riemann_fallbacks,
//   positivity_fallbacks, stencil_degradations, ctu_degraded_states,
//   scratch_high_water_bytes, moved_bytes, l1_density_error, wall_ms.
// Tile/batch sizes (HydroActiveBatchPolicy::max_active_cells) are candidates
// passed as the optional third argument; none is assumed universally best.
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

namespace hy = cosmosim::hydro;

constexpr double k_gamma = 5.0 / 3.0;
constexpr double k_two_pi = 6.283185307179586;
constexpr double k_u0 = 0.5;

double exactDensity(double x, double y, double z, double time, double box) {
  return 1.0 + 0.2 * std::sin(k_two_pi * ((x - k_u0 * time) + (y - k_u0 * time) + (z - k_u0 * time)) / box);
}

struct Configuration {
  hy::HydroReconstructionScheme reconstruction;
  hy::HydroIntegrationScheme integration;
  hy::HydroRiemannSolverKind riemann;
};

std::unique_ptr<hy::HydroReconstruction> makeReconstruction(
    hy::HydroReconstructionScheme scheme, const std::array<double, 3>& courant) {
  const hy::HydroReconstructionPolicy policy{
      .dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
      .rho_floor = 1.0e-12, .pressure_floor = 1.0e-12, .adiabatic_index = k_gamma};
  switch (scheme) {
    case hy::HydroReconstructionScheme::kPiecewiseConstant:
      return std::make_unique<hy::PiecewiseConstantReconstruction>();
    case hy::HydroReconstructionScheme::kPiecewiseLinear:
      return std::make_unique<hy::MusclHancockReconstruction>(policy);
    case hy::HydroReconstructionScheme::kPiecewiseParabolic:
      return std::make_unique<hy::PpmReconstruction>(policy);
  }
  return nullptr;
}

void runConfiguration(const Configuration& configuration, std::size_t n, std::size_t steps,
                      std::size_t batch_cells) {
  const double box = 1.0;
  const double dx = box / static_cast<double>(n);
  hy::HydroPatchGeometry geometry = hy::makeCartesianPatchGeometry(hy::HydroCartesianPatchSpec{
      .nx = n, .ny = n, .nz = n,
      .cell_width_x_comoving = dx, .cell_width_y_comoving = dx, .cell_width_z_comoving = dx});
  hy::appendCartesianBoundaryGhostFaces(geometry, hy::HydroBoundaryKind::kPeriodic);
  (void)hy::completePeriodicReconstructionStencils(geometry);
  hy::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t i = 0; i < n; ++i) {
        hy::HydroPrimitiveState primitive;
        primitive.rho_comoving = exactDensity((i + 0.5) * dx, (j + 0.5) * dx, (k + 0.5) * dx, 0.0, box);
        primitive.vel_x_peculiar = k_u0;
        primitive.vel_y_peculiar = k_u0;
        primitive.vel_z_peculiar = k_u0;
        primitive.pressure_comoving = 1.0;
        conserved.storeCell(geometry.linearCellIndex(i, j, k),
                            hy::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
      }
    }
  }
  std::vector<std::size_t> active_cells(geometry.cellCount());
  for (std::size_t c = 0; c < active_cells.size(); ++c) active_cells[c] = c;
  std::vector<std::size_t> active_faces;
  for (std::size_t f = 0; f < geometry.faces.size(); ++f) {
    const hy::HydroFace& face = geometry.faces[f];
    if (face.ghost_cell_slot != hy::k_invalid_ghost_cell_slot &&
        geometry.ghost_cells[face.ghost_cell_slot].side == hy::HydroFaceSide::kLower) {
      continue;  // periodic box face evaluated once (upper copy)
    }
    active_faces.push_back(f);
  }
  // Signal speed ~ sqrt(3) u0 + c; Courant 0.4 (<= 1/2 for 3D CTU).
  const double c_sound = std::sqrt(k_gamma);
  const double dt = 0.4 * dx / (std::sqrt(3.0) * k_u0 + c_sound);
  const hy::HydroUpdateContext update{.dt_code = dt, .scale_factor = 1.0};
  const hy::HydroSourceContext source_context{.update = update};
  const auto courant = hy::hydroTransportCourantFactors(dt, 1.0, dx, dx, dx);
  const auto reconstruction = makeReconstruction(configuration.reconstruction, courant);
  const auto riemann = hy::makeHydroRiemannSolver(configuration.riemann);
  const hy::HydroCoreSolver solver(k_gamma, hy::HydroNumericsPolicy{.integration_scheme = configuration.integration});
  const hy::HydroActiveBatchPolicy batch =
      solver.requiredWorkspaceFeatures(hy::HydroActiveBatchPolicy{.max_active_cells = batch_cells});
  hy::HydroScratchBuffers scratch;
  (void)scratch.prepareForActiveSet(active_cells.size(), active_faces.size(), geometry.cellCount(),
                                    geometry.totalCellStorageCount(), batch);
  hy::HydroProfileEvent profile;
  const auto start = std::chrono::steady_clock::now();
  for (std::size_t step = 0; step < steps; ++step) {
    hy::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);
    solver.advancePatchActiveSetBatchedWithScratch(
        conserved, geometry, hy::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
        update, *reconstruction, *riemann, {}, source_context, scratch, nullptr, batch, &profile);
  }
  const double wall_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
  double l1 = 0.0;
  const double time = dt * static_cast<double>(steps);
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t i = 0; i < n; ++i) {
        const double rho = conserved.loadCell(geometry.linearCellIndex(i, j, k)).mass_density_comoving;
        l1 += std::abs(rho - exactDensity((i + 0.5) * dx, (j + 0.5) * dx, (k + 0.5) * dx, time, box));
      }
    }
  }
  l1 /= static_cast<double>(geometry.cellCount());
  const double seconds = wall_ms * 1.0e-3;
  std::cout << "reconstruction=" << hy::hydroReconstructionSchemeToString(configuration.reconstruction)
            << " integration=" << hy::hydroIntegrationSchemeToString(configuration.integration)
            << " riemann=" << hy::hydroRiemannSolverKindToString(configuration.riemann)
            << " cells_per_axis=" << n << " steps=" << steps << " batch_cells=" << batch_cells
            << " wall_ms=" << wall_ms
            << " cells_per_s=" << (seconds > 0.0 ? geometry.cellCount() * steps / seconds : 0.0)
            << " faces_per_s=" << (seconds > 0.0 ? active_faces.size() * steps / seconds : 0.0)
            << " riemann_solves=" << profile.riemann_solve_count
            << " riemann_fallbacks=" << profile.riemann_fallback_count
            << " positivity_fallbacks=" << profile.positivity_fallback_count
            << " stencil_degradations=" << profile.stencil_degradation_count
            << " ctu_degraded_states=" << profile.ctu_transverse_degraded_states
            << " scratch_high_water_bytes=" << scratch.high_water.capacity_bytes
            << " moved_bytes=" << profile.bytes_moved
            << " l1_density_error=" << l1 << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  const std::size_t n = argc > 1 ? static_cast<std::size_t>(std::strtoull(argv[1], nullptr, 10)) : 32U;
  const std::size_t steps = argc > 2 ? static_cast<std::size_t>(std::strtoull(argv[2], nullptr, 10)) : 8U;
  const std::size_t batch = argc > 3 ? static_cast<std::size_t>(std::strtoull(argv[3], nullptr, 10)) : 0U;
  for (auto reconstruction : {hy::HydroReconstructionScheme::kPiecewiseConstant,
                              hy::HydroReconstructionScheme::kPiecewiseLinear,
                              hy::HydroReconstructionScheme::kPiecewiseParabolic}) {
    for (auto integration : {hy::HydroIntegrationScheme::kMusclHancock,
                             hy::HydroIntegrationScheme::kCornerTransportUpwind}) {
      for (auto riemann : {hy::HydroRiemannSolverKind::kHllc, hy::HydroRiemannSolverKind::kHlle}) {
        runConfiguration(Configuration{reconstruction, integration, riemann}, n, steps, batch);
      }
    }
  }
  return 0;
}
