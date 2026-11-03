// Audit F7: effective-ISM energy ledger is an integrated energy.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants:
//  * ledger entry = rho_comoving * delta_u * V_comoving (code energy);
//  * one parent cell of volume V and eight equivalent children of volume V/8
//    (same density and thermodynamic state) record the same total energy
//    (the defective density-sum ledger grew eightfold);
//  * unequal-volume cells weight by their own volume;
//  * a missing cell volume is refused instead of silently recording an
//    energy density;
//  * HydroCoreSolver injects the patch volume, so the ledger matches the
//    solver's own volume-integrated source_delta.total_energy.
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "../support/test_require.hpp"

namespace {

namespace hy = cosmosim::hydro;
namespace ph = cosmosim::physics;

ph::EffectiveMultiphaseEosTable makeTable() {
  cosmosim::core::PhysicsConfig physics;
  physics.sf_effective_q_eos = 0.3;
  physics.sf_effective_eos_table_bins = 128;
  physics.sf_effective_eos_max_density_ratio = 1.0e5;
  physics.sf_effective_supernova_specific_energy_code = 1.0e5;
  const auto units = cosmosim::core::makeUnitSystem("kpc", "msun", "km_s");
  return ph::EffectiveMultiphaseEosTable(
      ph::makeEffectiveMultiphaseEosConfig(physics), units,
      ph::makeEffectiveIsmReferenceCoolingProvider(physics));
}

struct ColdCell {
  hy::HydroPrimitiveState primitive;
  hy::HydroConservedState conserved;
  double target_u = 0.0;
};

ColdCell coldEffectiveCell(const ph::EffectiveMultiphaseEosTable& table) {
  const double rho = 10.0 * table.thresholdDensityPhysCode();
  const auto equilibrium = table.lookup(rho);
  ColdCell cell;
  cell.target_u = equilibrium.entry.specific_internal_energy_eff_code;
  cell.primitive.rho_comoving = rho;
  cell.primitive.specific_internal_energy_code = 0.5 * cell.target_u;
  cell.primitive.pressure_comoving = (5.0 / 3.0 - 1.0) * rho * cell.primitive.specific_internal_energy_code;
  cell.primitive.signal_speed_squared_code = (5.0 / 3.0) * cell.primitive.pressure_comoving / rho;
  cell.conserved = hy::HydroCoreSolver::conservedFromPrimitive(cell.primitive, 5.0 / 3.0);
  return cell;
}

bool near(double a, double b, double tol) {
  return std::abs(a - b) <= tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

void testParentVersusEightChildren() {
  ph::EffectiveIsmThermodynamicClosure closure(makeTable());
  const ColdCell cell = coldEffectiveCell(closure.table());
  ph::EffectiveIsmEnergyRelaxationSource parent_source(closure);
  ph::EffectiveIsmEnergyRelaxationSource child_source(closure);
  hy::HydroSourceContext context{.update = hy::HydroUpdateContext{.dt_code = 0.25, .scale_factor = 1.0}};
  context.cell_volume_comoving = 8.0;
  (void)parent_source.sourceForCell(0U, cell.conserved, cell.primitive, context);
  context.cell_volume_comoving = 1.0;
  for (std::size_t child = 0; child < 8U; ++child) {
    (void)child_source.sourceForCell(child, cell.conserved, cell.primitive, context);
  }
  const double expected = cell.primitive.rho_comoving * (cell.target_u - cell.primitive.specific_internal_energy_code) * 8.0;
  CHUI_REQUIRE(near(parent_source.ledger().energy_added_code, expected, 5.0e-3));
  CHUI_REQUIRE(near(child_source.ledger().energy_added_code, parent_source.ledger().energy_added_code, 1.0e-12));
  CHUI_REQUIRE(child_source.ledger().adjusted_cell_count == 8U);
}

void testUnequalVolumesWeightByOwnVolume() {
  ph::EffectiveIsmThermodynamicClosure closure(makeTable());
  const ColdCell cell = coldEffectiveCell(closure.table());
  ph::EffectiveIsmEnergyRelaxationSource unit_source(closure);
  ph::EffectiveIsmEnergyRelaxationSource mixed_source(closure);
  hy::HydroSourceContext context{.update = hy::HydroUpdateContext{.dt_code = 0.25, .scale_factor = 1.0}};
  context.cell_volume_comoving = 1.0;
  (void)unit_source.sourceForCell(0U, cell.conserved, cell.primitive, context);
  (void)mixed_source.sourceForCell(0U, cell.conserved, cell.primitive, context);
  context.cell_volume_comoving = 3.0;
  (void)mixed_source.sourceForCell(1U, cell.conserved, cell.primitive, context);
  CHUI_REQUIRE(near(mixed_source.ledger().net_energy_adjustment_code,
                    4.0 * unit_source.ledger().net_energy_adjustment_code, 1.0e-12));
}

// Coordinate convention: the same physical cell described in an isolated
// (a = 1) frame and in a comoving frame at a = 0.5 has rho_c = a^3 rho_phys
// and V_c = V_phys / a^3, so rho_c V_c = m and the integrated ledger
// m * delta_u must agree. The closure forms the physical density rho_c/a^3,
// so the equilibrium target is the same in both descriptions.
void testComovingAndIsolatedConventionsAgree() {
  ph::EffectiveIsmThermodynamicClosure closure(makeTable());
  const ColdCell physical = coldEffectiveCell(closure.table());
  const double a = 0.5;
  const double volume_physical = 2.0;
  ColdCell comoving = physical;
  comoving.primitive.rho_comoving = physical.primitive.rho_comoving * a * a * a;
  comoving.primitive.pressure_comoving = physical.primitive.pressure_comoving * a * a * a;
  comoving.conserved = hy::HydroCoreSolver::conservedFromPrimitive(comoving.primitive, 5.0 / 3.0);
  ph::EffectiveIsmEnergyRelaxationSource isolated_source(closure);
  ph::EffectiveIsmEnergyRelaxationSource comoving_source(closure);
  hy::HydroSourceContext isolated{.update = hy::HydroUpdateContext{.dt_code = 0.25, .scale_factor = 1.0}};
  isolated.cell_volume_comoving = volume_physical;
  hy::HydroSourceContext cosmological{
      .update = hy::HydroUpdateContext{.dt_code = 0.25, .scale_factor = a}, .redshift = 1.0 / a - 1.0};
  cosmological.cell_volume_comoving = volume_physical / (a * a * a);
  (void)isolated_source.sourceForCell(0U, physical.conserved, physical.primitive, isolated);
  (void)comoving_source.sourceForCell(0U, comoving.conserved, comoving.primitive, cosmological);
  CHUI_REQUIRE(isolated_source.ledger().adjusted_cell_count == 1U);
  CHUI_REQUIRE(near(comoving_source.ledger().net_energy_adjustment_code,
                    isolated_source.ledger().net_energy_adjustment_code, 1.0e-10));
}

void testMissingVolumeIsRefused() {
  ph::EffectiveIsmThermodynamicClosure closure(makeTable());
  const ColdCell cell = coldEffectiveCell(closure.table());
  ph::EffectiveIsmEnergyRelaxationSource source(closure);
  const hy::HydroSourceContext context{.update = hy::HydroUpdateContext{.dt_code = 0.25, .scale_factor = 1.0}};
  bool refused = false;
  try {
    (void)source.sourceForCell(0U, cell.conserved, cell.primitive, context);
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
}

void testSolverInjectsPatchVolumeAndLedgerMatchesSourceDelta() {
  ph::EffectiveIsmThermodynamicClosure closure(makeTable());
  const ColdCell cell = coldEffectiveCell(closure.table());
  ph::EffectiveIsmEnergyRelaxationSource source(closure);
  // 2x1x1 patch with dx = 0.5: V = 0.125. No faces between identical cells
  // produce flux; only the source acts.
  hy::HydroPatchGeometry geometry = hy::makeCartesianPatchGeometry(hy::HydroCartesianPatchSpec{
      .nx = 2, .ny = 1, .nz = 1,
      .cell_width_x_comoving = 0.5, .cell_width_y_comoving = 0.5, .cell_width_z_comoving = 0.5});
  hy::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  conserved.storeCell(0U, cell.conserved);
  conserved.storeCell(1U, cell.conserved);
  const std::vector<std::size_t> active_cells{0U, 1U};
  std::vector<std::size_t> active_faces;
  for (std::size_t face = 0; face < geometry.faces.size(); ++face) active_faces.push_back(face);
  const hy::HydroUpdateContext update{.dt_code = 0.25, .scale_factor = 1.0};
  const hy::HydroSourceContext source_context{.update = update};
  const hy::HydroCoreSolver solver(5.0 / 3.0);
  const hy::PiecewiseConstantReconstruction pcm;
  hy::HllcRiemannSolver riemann;
  hy::HydroScratchBuffers scratch;
  hy::HydroProfileEvent profile;
  const hy::HydroSourceTerm* const sources[] = {&source};
  solver.advancePatchActiveSetWithScratch(
      conserved, geometry, hy::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update, pcm, riemann, sources, source_context, scratch, nullptr, &profile);
  CHUI_REQUIRE(source.ledger().adjusted_cell_count == 2U);
  CHUI_REQUIRE(near(source.ledger().net_energy_adjustment_code,
                    profile.conservation.source_delta.total_energy, 1.0e-10));
}

}  // namespace

int main() {
  testParentVersusEightChildren();
  testUnequalVolumesWeightByOwnVolume();
  testComovingAndIsolatedConventionsAgree();
  testMissingVolumeIsRefused();
  testSolverInjectsPatchVolumeAndLedgerMatchesSourceDelta();
  return 0;
}
