// Dual-energy (auxiliary internal energy + total-energy sync) tests.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Problem 1: a cold contact advected at high Mach (rho varies, u = 10,
// P = 1e-4 uniform, K/(rho e) ~ 3e5, thermal fraction 3e-6 < eta). The
// exact solution keeps P uniform. The auxiliary internal-energy update
// advects rho e with the upwinded mass flux (F_rho e_upwind = u P/(gamma-1),
// uniform) and has zero pdV work for uniform u, so after the E sync P is
// exact to round-off; the switch must fire.
// Problem 2: Galilean-boosted smooth sound wave (U0 = 1e4, du = 1e-2,
// P = 1): the auxiliary pressure must stay within 1e-3 of the rest-frame
// pressure after one step (the physical change is ~1e-6).
//
// Invariants:
//  * relative pressure error < 1e-8 (problem 1) and < 1e-3 (problem 2);
//  * the switch fired (dual_energy_sync_count > 0) and the measured total
//    energy change is reported in dual_energy_sync_delta, so the conservation
//    residual (which excludes it) stays at round-off;
//  * dual energy requires the planned workspace lane (fail-closed contract).
//
// Evidence limit (documented in docs/hydro_numerics.md): a one-step unit
// test cannot demonstrate the multi-step high-Mach benefit; that is a
// human-run validation target. The synchronized scheme also cannot recover
// thermal energy below ulp(E) (K/(rho e) >~ 1e13), which would require a
// persistent entropy lane.
#include <cmath>
#include <cstddef>
#include <stdexcept>
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
constexpr std::size_t k_n = 16U;

struct Result {
  double max_relative_pressure_error = 0.0;
  hy::HydroProfileEvent profile;
};

Result runCold(const hy::HydroNumericsPolicy& numerics, double velocity, double pressure,
               double density_amplitude, double velocity_amplitude) {
  hy::HydroPatchGeometry geometry = hy::makeCartesianPatchGeometry(hy::HydroCartesianPatchSpec{
      .nx = k_n, .ny = 1, .nz = 1,
      .cell_width_x_comoving = 1.0, .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  hy::appendCartesianBoundaryGhostFaces(geometry, hy::HydroBoundaryKind::kPeriodic);
  (void)hy::completePeriodicReconstructionStencils(geometry);
  hy::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  for (std::size_t i = 0; i < k_n; ++i) {
    hy::HydroPrimitiveState primitive;
    const double phase = 2.0 * 3.141592653589793 * static_cast<double>(i) / k_n;
    primitive.rho_comoving = 1.0 + density_amplitude * std::sin(phase);
    primitive.vel_x_peculiar = velocity + velocity_amplitude * std::sin(phase);
    primitive.pressure_comoving = pressure;
    conserved.storeCell(i, hy::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }
  hy::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);
  std::vector<std::size_t> active_cells;
  for (std::size_t i = 0; i < k_n; ++i) active_cells.push_back(i);
  std::vector<std::size_t> active_faces;
  for (std::size_t face = 0; face < geometry.faces.size(); ++face) {
    const hy::HydroFace& f = geometry.faces[face];
    if (f.ghost_cell_slot != hy::k_invalid_ghost_cell_slot &&
        geometry.ghost_cells[f.ghost_cell_slot].side == hy::HydroFaceSide::kLower) {
      continue;
    }
    active_faces.push_back(face);
  }
  // Courant number u dt/dx = 0.4.
  const double dt = 0.4 / velocity;
  const hy::HydroUpdateContext update{.dt_code = dt, .scale_factor = 1.0};
  const hy::HydroSourceContext source_context{.update = update};
  const auto courant = hy::hydroTransportCourantFactors(dt, 1.0, 1.0, 1.0, 1.0);
  const hy::MusclHancockReconstruction plm(hy::HydroReconstructionPolicy{
      .dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
      .rho_floor = 1.0e-30, .pressure_floor = 1.0e-30, .adiabatic_index = k_gamma});
  const hy::HydroCoreSolver solver(k_gamma, numerics);
  hy::HllcRiemannSolver riemann;
  hy::HydroScratchBuffers scratch;
  Result result;
  solver.advancePatchActiveSetWithScratch(
      conserved, geometry,
      hy::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update, plm, riemann, {}, source_context, scratch, nullptr, &result.profile);
  for (std::size_t i = 0; i < k_n; ++i) {
    const hy::HydroConservedState cell = conserved.loadCell(i);
    const double kinetic = 0.5 * cell.momentum_density_x_comoving * cell.momentum_density_x_comoving /
        cell.mass_density_comoving;
    const double cell_pressure = (k_gamma - 1.0) * (cell.total_energy_density_comoving - kinetic);
    result.max_relative_pressure_error = std::max(
        result.max_relative_pressure_error, std::abs(cell_pressure - pressure) / pressure);
  }
  return result;
}

void testHighMachContactKeepsThermalPressure() {
  hy::HydroNumericsPolicy on;
  on.dual_energy.mode = hy::HydroDualEnergyMode::kThermalFraction;
  on.dual_energy.thermal_fraction_switch = 1.0e-3;
  const Result with_aux = runCold(on, 10.0, 1.0e-4, 0.5, 0.0);
  CHUI_REQUIRE(with_aux.max_relative_pressure_error < 1.0e-8);
  CHUI_REQUIRE(with_aux.profile.conservation.dual_energy_sync_count > 0U);
  CHUI_REQUIRE(with_aux.profile.dual_energy_auxiliary_cells > 0U);
  // Mass and momentum are untouched by the sync; the residual excludes the
  // separately reported sync delta.
  CHUI_REQUIRE(std::abs(with_aux.profile.conservation.residual.mass) < 1.0e-12);
  CHUI_REQUIRE(std::abs(with_aux.profile.conservation.residual.total_energy) <
               1.0e-9 * std::abs(with_aux.profile.conservation.before.total_energy));
}

void testBoostedSmoothWaveThermalPressure() {
  hy::HydroNumericsPolicy on;
  on.dual_energy.mode = hy::HydroDualEnergyMode::kThermalFraction;
  const Result boosted = runCold(on, 1.0e4, 1.0, 0.0, 1.0e-2);
  CHUI_REQUIRE(boosted.max_relative_pressure_error < 1.0e-3);
  CHUI_REQUIRE(boosted.profile.conservation.dual_energy_sync_count > 0U);
}

void testDualEnergyRequiresPlannedWorkspace() {
  hy::HydroNumericsPolicy on;
  on.dual_energy.mode = hy::HydroDualEnergyMode::kThermalFraction;
  const hy::HydroCoreSolver solver(k_gamma, on);
  const hy::HydroActiveBatchPolicy planned = solver.requiredWorkspaceFeatures(hy::HydroActiveBatchPolicy{});
  CHUI_REQUIRE(planned.dual_energy_workspace);
  CHUI_REQUIRE(!planned.ctu_transverse_workspace);
  // Invalid switch values are rejected at construction.
  bool rejected = false;
  try {
    hy::HydroNumericsPolicy bad = on;
    bad.dual_energy.thermal_fraction_switch = 1.5;
    const hy::HydroCoreSolver invalid(k_gamma, bad);
    (void)invalid;
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);
}

}  // namespace

int main() {
  testHighMachContactKeepsThermalPressure();
  testBoostedSmoothWaveThermalPressure();
  testDualEnergyRequiresPlannedWorkspace();
  return 0;
}
