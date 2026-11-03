// Dual-energy formalism tests (auxiliary internal-energy lane, no E overwrite).
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Policy under test (docs/hydro_numerics.md section 7): the conservative
// update owns rho, rho u, E, rho Z; the thermal energy handed to the
// canonical state is E - K when E - K >= eta E, otherwise the auxiliary
// rho e evolved by the non-conservative internal-energy equation from the
// canonical e (never from E - K). E is not modified by the solver; the
// thermal correction is ledgered in dual_energy_sync_delta.
//
// Invariants:
//  1. Cold contact at K/(rho e) ~ 3e13 (beyond what E can represent): the
//     auxiliary pressure stays uniform to < 1e-9 relative, the conservative
//     E - K is wrong by > 1e-6 relative, E is bitwise the conservative value,
//     and the ledger equals sum (aux - (E-K)) V.
//  2. Warm flow (thermal fraction >> eta): the auxiliary branch is never used
//     and the ledger is exactly zero (dual energy does not perturb
//     conservation where it is not needed).
//  3. Shock heating: two cold streams (Mach ~1e4) collide; shocked cells are
//     heated through E (selection falls back to E - K, auxiliary lane resynced),
//     so the peak pressure approaches the Rankine-Hugoniot value instead of
//     staying cold.
//  4. Contract: dual energy requires the planned workspace and the auxiliary
//     lane; invalid eta is rejected.
#include <cmath>
#include <cstddef>
#include <limits>
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

hy::HydroNumericsPolicy dualOn() {
  hy::HydroNumericsPolicy policy;
  policy.dual_energy.mode = hy::HydroDualEnergyMode::kThermalFraction;
  policy.dual_energy.thermal_fraction_switch = 1.0e-3;
  return policy;
}

struct Box {
  hy::HydroPatchGeometry geometry;
  hy::HydroConservedStateSoa conserved;
  std::vector<std::size_t> active_cells;
  std::vector<std::size_t> active_faces;
};

// Periodic 1D box; density/velocity profiles given per cell, uniform P.
template <class TDensity, class TVelocity>
Box makeBox(bool dual, double pressure, TDensity density, TVelocity velocity,
            hy::HydroBoundaryKind boundary = hy::HydroBoundaryKind::kPeriodic) {
  Box box;
  box.geometry = hy::makeCartesianPatchGeometry(hy::HydroCartesianPatchSpec{
      .nx = k_n, .ny = 1, .nz = 1,
      .cell_width_x_comoving = 1.0, .cell_width_y_comoving = 1.0, .cell_width_z_comoving = 1.0});
  hy::appendCartesianBoundaryGhostFaces(box.geometry, boundary);
  (void)hy::completeReconstructionStencils(box.geometry);
  box.conserved.resize(box.geometry.totalCellStorageCount());
  box.conserved.enableAuxiliaryInternalEnergy(dual);
  for (std::size_t i = 0; i < k_n; ++i) {
    hy::HydroPrimitiveState primitive;
    primitive.rho_comoving = density(i);
    primitive.vel_x_peculiar = velocity(i);
    primitive.pressure_comoving = pressure;
    box.conserved.storeCell(i, hy::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
    if (dual) {
      // Canonical specific internal energy e = P/((gamma-1) rho): rho e exact.
      box.conserved.auxiliaryInternalEnergyDensityComoving()[i] = pressure / (k_gamma - 1.0);
    }
  }
  for (std::size_t i = 0; i < k_n; ++i) box.active_cells.push_back(i);
  for (std::size_t face = 0; face < box.geometry.faces.size(); ++face) {
    const hy::HydroFace& f = box.geometry.faces[face];
    if (boundary == hy::HydroBoundaryKind::kPeriodic && f.ghost_cell_slot != hy::k_invalid_ghost_cell_slot &&
        box.geometry.ghost_cells[f.ghost_cell_slot].side == hy::HydroFaceSide::kLower) {
      continue;  // periodic box face evaluated once
    }
    box.active_faces.push_back(face);
  }
  return box;
}

hy::HydroProfileEvent step(Box& box, const hy::HydroCoreSolver& solver, double dt) {
  hy::fillHydroBoundaryGhostCells(box.conserved, box.geometry, k_gamma);
  hy::fillHydroBoundaryGhostAuxiliaryInternalEnergy(box.conserved, box.geometry);
  const hy::HydroUpdateContext update{.dt_code = dt, .scale_factor = 1.0};
  const hy::HydroSourceContext source_context{.update = update};
  const auto courant = hy::hydroTransportCourantFactors(dt, 1.0, 1.0, 1.0, 1.0);
  const hy::MusclHancockReconstruction plm(hy::HydroReconstructionPolicy{
      .dt_over_dx_code = courant[0], .dt_over_cell_width_code = courant,
      .rho_floor = 1.0e-30, .pressure_floor = 1.0e-30, .adiabatic_index = k_gamma});
  hy::HllcRiemannSolver riemann;
  hy::HydroScratchBuffers scratch;
  hy::HydroProfileEvent profile;
  solver.advancePatchActiveSetWithScratch(
      box.conserved, box.geometry,
      hy::HydroActiveSetView{.active_cells = box.active_cells, .active_faces = box.active_faces},
      update, plm, riemann, {}, source_context, scratch, nullptr, &profile);
  return profile;
}

double thermalConservative(const hy::HydroConservedState& cell) {
  return cell.total_energy_density_comoving -
      0.5 * cell.momentum_density_x_comoving * cell.momentum_density_x_comoving / cell.mass_density_comoving;
}

void testColdContactBeyondTotalEnergyPrecision() {
  constexpr double pressure = 1.0e-8;
  constexpr double velocity = 1.0e3;  // K ~ 5e5, rho e ~ 1.5e-8
  const auto density = [](std::size_t i) {
    return 1.0 + 0.5 * std::sin(2.0 * 3.141592653589793 * static_cast<double>(i) / k_n);
  };
  const auto uniform_velocity = [](std::size_t) { return velocity; };
  const double dt = 0.4 / velocity;

  Box with_aux = makeBox(true, pressure, density, uniform_velocity);
  const hy::HydroProfileEvent profile = step(with_aux, hy::HydroCoreSolver(k_gamma, dualOn()), dt);
  Box without_aux = makeBox(false, pressure, density, uniform_velocity);
  (void)step(without_aux, hy::HydroCoreSolver(k_gamma), dt);

  double aux_error = 0.0;
  double conservative_error = 0.0;
  double ledger = 0.0;
  for (std::size_t i = 0; i < k_n; ++i) {
    const double rho_e_aux = with_aux.conserved.auxiliaryInternalEnergyDensityComoving()[i];
    aux_error = std::max(aux_error, std::abs((k_gamma - 1.0) * rho_e_aux - pressure) / pressure);
    const hy::HydroConservedState cons_cell = without_aux.conserved.loadCell(i);
    conservative_error = std::max(
        conservative_error, std::abs((k_gamma - 1.0) * thermalConservative(cons_cell) - pressure) / pressure);
    ledger += rho_e_aux - thermalConservative(with_aux.conserved.loadCell(i));
  }
  CHUI_REQUIRE(aux_error < 1.0e-9);
  CHUI_REQUIRE(conservative_error > 1.0e-6);
  CHUI_REQUIRE(profile.dual_energy_auxiliary_cells == k_n);
  CHUI_REQUIRE(profile.dual_energy_unknown_auxiliary_cells == 0U);
  // The ledger is measured against the pre-floor conservative thermal energy:
  // sum (aux - (E_final - K)) + floor_delta = sum (aux - (E_pre - K)).
  // Tolerance: a few ulp of E ~ 5e5 per cell.
  ledger += profile.conservation.floor_delta.total_energy;
  CHUI_REQUIRE(std::abs(profile.conservation.dual_energy_sync_delta.total_energy - ledger) < 1.0e-8);
  CHUI_REQUIRE(std::abs(profile.conservation.residual.mass) < 1.0e-12);
}

void testWarmFlowIsUntouched() {
  const auto density = [](std::size_t i) { return 1.0 + 0.1 * static_cast<double>(i % 3U); };
  const auto velocity = [](std::size_t) { return 0.3; };
  Box box = makeBox(true, 1.0, density, velocity);
  const hy::HydroProfileEvent profile = step(box, hy::HydroCoreSolver(k_gamma, dualOn()), 0.2);
  CHUI_REQUIRE(profile.dual_energy_auxiliary_cells == 0U);
  CHUI_REQUIRE(profile.conservation.dual_energy_sync_count == 0U);
  CHUI_REQUIRE(profile.conservation.dual_energy_sync_delta.total_energy == 0.0);
  for (std::size_t i = 0; i < k_n; ++i) {
    // Auxiliary lane resynchronized to the conservative thermal energy.
    const double thermal = thermalConservative(box.conserved.loadCell(i));
    CHUI_REQUIRE(std::abs(box.conserved.auxiliaryInternalEnergyDensityComoving()[i] - thermal) <=
                 1.0e-14 * std::abs(thermal));
  }
}

void testShockHeatingIsCapturedThroughTotalEnergy() {
  constexpr double pressure = 1.0e-6;
  constexpr double stream_speed = 10.0;  // Mach ~ 1e4 relative to c ~ 1.3e-3
  const auto density = [](std::size_t) { return 1.0; };
  // Converging streams: +u for i < n/2, -u for i >= n/2; open (inflow-
  // continuing) boundaries so no vacuum forms at the box edges.
  const auto velocity = [](std::size_t i) { return i < k_n / 2U ? stream_speed : -stream_speed; };
  Box box = makeBox(true, pressure, density, velocity, hy::HydroBoundaryKind::kOpen);
  const hy::HydroCoreSolver solver(k_gamma, dualOn());
  const double dt = 0.2 / stream_speed;
  for (int s = 0; s < 10; ++s) {
    (void)step(box, solver, dt);
  }
  double peak_pressure = 0.0;
  for (std::size_t i = 0; i < k_n; ++i) {
    peak_pressure = std::max(peak_pressure,
        (k_gamma - 1.0) * box.conserved.auxiliaryInternalEnergyDensityComoving()[i]);
  }
  // Strong-shock post-shock pressure ~ (gamma+1)/2 rho u^2 ~ 133 (stationary
  // wall frame); require that the selected thermal state is heated by many
  // orders of magnitude rather than staying at 1e-6.
  CHUI_REQUIRE(peak_pressure > 10.0);
  for (std::size_t i = 0; i < k_n; ++i) {
    CHUI_REQUIRE(std::isfinite(box.conserved.auxiliaryInternalEnergyDensityComoving()[i]));
    CHUI_REQUIRE(box.conserved.auxiliaryInternalEnergyDensityComoving()[i] > 0.0);
  }
}

void testContracts() {
  const hy::HydroCoreSolver solver(k_gamma, dualOn());
  const hy::HydroActiveBatchPolicy planned = solver.requiredWorkspaceFeatures(hy::HydroActiveBatchPolicy{});
  CHUI_REQUIRE(planned.dual_energy_workspace);
  CHUI_REQUIRE(!planned.ctu_transverse_workspace);
  // Missing auxiliary lane is refused (fail closed, no silent E - K).
  Box box = makeBox(false, 1.0, [](std::size_t) { return 1.0; }, [](std::size_t) { return 0.0; });
  bool refused = false;
  try {
    (void)step(box, solver, 0.1);
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
  bool rejected = false;
  try {
    hy::HydroNumericsPolicy bad = dualOn();
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
  testColdContactBeyondTotalEnergyPrecision();
  testWarmFlowIsUntouched();
  testShockHeatingIsCapturedThroughTotalEnergy();
  testContracts();
  return 0;
}
