#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"

namespace {

bool nearlyEqual(double a, double b, double tolerance = 1.0e-8) {
  return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

cosmosim::physics::EffectiveMultiphaseEosTable makeTable(double q_eos = 0.3) {
  cosmosim::core::PhysicsConfig physics;
  physics.sf_effective_q_eos = q_eos;
  physics.sf_effective_eos_table_bins = 128;
  physics.sf_effective_eos_max_density_ratio = 1.0e5;
  physics.sf_effective_supernova_specific_energy_code = 1.0e5;
  const auto units = cosmosim::core::makeUnitSystem("kpc", "msun", "km_s");
  return cosmosim::physics::EffectiveMultiphaseEosTable(
      cosmosim::physics::makeEffectiveMultiphaseEosConfig(physics),
      units,
      cosmosim::physics::makeEffectiveIsmReferenceCoolingProvider(physics));
}

void testThresholdAndDensityScaling() {
  const auto table = makeTable();
  assert(table.thresholdDensityPhysCode() > 0.0);
  assert(table.thresholdDensityCgs() > 0.0);
  const auto below = table.lookup(0.99 * table.thresholdDensityPhysCode());
  assert(!below.above_threshold);
  const auto threshold = table.lookup(table.thresholdDensityPhysCode());
  const auto dense = table.lookup(100.0 * table.thresholdDensityPhysCode());
  assert(threshold.above_threshold && threshold.valid);
  assert(dense.above_threshold && dense.valid);
  assert(nearlyEqual(
      dense.entry.star_formation_timescale_code,
      threshold.entry.star_formation_timescale_code / 10.0,
      2.0e-3));
  assert(dense.entry.cold_mass_fraction >= 0.0 && dense.entry.cold_mass_fraction <= 1.0);
  assert(dense.entry.pressure_phys_code > 0.0);
  assert(dense.entry.signal_speed_squared_code > 0.0);
}

void testTableDeterminismMonotonicityAndDerivative() {
  const auto first = makeTable();
  const auto second = makeTable();
  assert(first.tableHash() == second.tableHash());
  assert(!first.tableHashHex().empty());
  const auto entries = first.entries();
  assert(entries.size() == 128U);
  for (std::size_t i = 0; i < entries.size(); ++i) {
    assert(entries[i].valid);
    assert(entries[i].cold_mass_fraction >= 0.0 && entries[i].cold_mass_fraction <= 1.0);
    assert(entries[i].pressure_phys_code > 0.0);
    assert(entries[i].signal_speed_squared_code > 0.0);
    if (i > 0U) {
      assert(entries[i].density_phys_code > entries[i - 1U].density_phys_code);
      assert(entries[i].pressure_phys_code >= entries[i - 1U].pressure_phys_code);
    }
  }
  const std::size_t mid = entries.size() / 2U;
  const double finite_difference =
      (entries[mid + 1U].pressure_phys_code - entries[mid - 1U].pressure_phys_code) /
      (entries[mid + 1U].density_phys_code - entries[mid - 1U].density_phys_code);
  assert(nearlyEqual(entries[mid].signal_speed_squared_code, finite_difference, 0.08));
}

void testQeosLimitsAndThresholdContinuity() {
  const auto iso = makeTable(0.0);
  const auto softened = makeTable(0.3);
  const auto full = makeTable(1.0);
  const double density = 10.0 * iso.thresholdDensityPhysCode();
  const auto a = iso.lookup(density);
  const auto b = softened.lookup(density);
  const auto c = full.lookup(density);
  assert(a.valid && b.valid && c.valid);
  const double lo = std::min(a.entry.specific_internal_energy_eff_code,
                             c.entry.specific_internal_energy_eff_code);
  const double hi = std::max(a.entry.specific_internal_energy_eff_code,
                             c.entry.specific_internal_energy_eff_code);
  assert(b.entry.specific_internal_energy_eff_code >= lo);
  assert(b.entry.specific_internal_energy_eff_code <= hi);
  const auto just_above = softened.lookup(
      softened.thresholdDensityPhysCode() * (1.0 + 1.0e-8));
  assert(just_above.valid);
  const double u_iso = softened.specificInternalEnergyFromTemperatureCode(
      softened.config().isothermal_temperature_k,
      softened.config().mean_molecular_weight_ionized);
  assert(nearlyEqual(just_above.entry.specific_internal_energy_eff_code, u_iso, 2.0e-3));
}

void testHydroClosureAndEnergyLedger() {
  auto table = makeTable();
  const double rho_phys = 10.0 * table.thresholdDensityPhysCode();
  const auto equilibrium = table.lookup(rho_phys);
  assert(equilibrium.valid);
  cosmosim::physics::EffectiveIsmThermodynamicClosure closure(std::move(table));

  cosmosim::hydro::HydroPrimitiveState primitive;
  primitive.rho_comoving = rho_phys;
  primitive.specific_internal_energy_code = 0.5 * equilibrium.entry.specific_internal_energy_eff_code;
  primitive.pressure_comoving = (5.0 / 3.0 - 1.0) * primitive.rho_comoving *
      primitive.specific_internal_energy_code;
  primitive.signal_speed_squared_code = (5.0 / 3.0) * primitive.pressure_comoving /
      primitive.rho_comoving;
  const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, 5.0 / 3.0);
  const auto result = closure.evaluate(0U, conserved, primitive, 1.0, 0.0);
  assert(result.valid && result.uses_effective_ism);
  assert(nearlyEqual(result.target_specific_internal_energy_code,
                     equilibrium.entry.specific_internal_energy_eff_code, 2.0e-3));
  assert(result.signal_speed_squared_code > 0.0);

  cosmosim::physics::EffectiveIsmEnergyRelaxationSource source(closure);
  cosmosim::hydro::HydroSourceContext context{
      .update = cosmosim::hydro::HydroUpdateContext{.dt_code = 0.25, .scale_factor = 1.0},
      .redshift = 0.0};
  // Audit F7 contract: the ledger integrates over the authoritative volume.
  context.cell_volume_comoving = 1.0;
  const auto delta = source.sourceForCell(0U, conserved, primitive, context);
  assert(delta.total_energy_density_comoving > 0.0);
  const auto ledger = source.ledger();
  assert(ledger.adjusted_cell_count == 1U);
  assert(ledger.energy_added_code > 0.0);
  assert(ledger.energy_removed_code == 0.0);
  assert(nearlyEqual(ledger.net_energy_adjustment_code, ledger.energy_added_code));

  auto hot = primitive;
  hot.specific_internal_energy_code = equilibrium.entry.specific_internal_energy_eff_code * 3.0;
  hot.pressure_comoving = (5.0 / 3.0 - 1.0) * hot.rho_comoving * hot.specific_internal_energy_code;
  hot.signal_speed_squared_code = (5.0 / 3.0) * hot.pressure_comoving / hot.rho_comoving;
  const auto hot_conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(hot, 5.0 / 3.0);
  const auto hot_result = closure.evaluate(0U, hot_conserved, hot, 1.0, 0.0);
  assert(!hot_result.uses_effective_ism);
  assert(nearlyEqual(hot_result.pressure_comoving, hot.pressure_comoving));
}

void testClosureSharesReadOnlyTable() {
  auto table = std::make_shared<const cosmosim::physics::EffectiveMultiphaseEosTable>(
      makeTable());
  cosmosim::physics::EffectiveIsmThermodynamicClosure closure(table);
  assert(&closure.table() == table.get());
  assert(closure.table().tableHash() == table->tableHash());
}

}  // namespace

int main() {
  testThresholdAndDensityScaling();
  testTableDeterminismMonotonicityAndDerivative();
  testQeosLimitsAndThresholdContinuity();
  testHydroClosureAndEnergyLedger();
  testClosureSharesReadOnlyTable();
  return 0;
}
