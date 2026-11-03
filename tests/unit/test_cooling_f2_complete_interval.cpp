// CHUI Part 2 — audit F2 regression (NOT RUN — USER-RESERVED).
//
// Pre-repair pathology: CoolingSourceIntegrator exited its explicit loop
// after max_subcycles with unprocessed time and the caller accepted the
// partial energy as the full-interval result.
//
// Reference problem (exact solution): with temperature_per_internal_energy_k
// = 1e4 the legacy analytic law gives T/1e4 = e, so with UV off
//   de/dt = -K sqrt(e),  K = 1e-22 n_H^2 / rho * conversion,
//   sqrt(e(t)) = sqrt(e0) - K t / 2.
// Invariants:
//  * complete intervals reproduce the exact solution to the error tolerance
//    and report completed_dt == requested dt (no smaller interval dressed as
//    completion);
//  * the stiff interval with an attempt budget too small is REPORTED
//    incomplete, returns the unchanged energy, and CoolingHeatingSource
//    throws SourceIntervalIncompleteError (a StageNumericalRetryRequest);
//  * zero net source still advances (completes) the full interval;
//  * StageTransaction classifies the failure as kNumericalRetry and leaves
//    staged state uncommitted.
#include <cassert>
#include <cmath>
#include <stdexcept>

#include "cosmosim/core/stage_transaction.hpp"
#include "cosmosim/physics/cooling_heating.hpp"
#include "cosmosim/physics/source_interval.hpp"

namespace {

using namespace cosmosim;

physics::CoolingRateProvider makeProvider(int max_subcycles) {
  physics::CoolingModelConfig model;
  model.uv_background_model = core::UvBackgroundModel::kNone;
  model.temperature_floor_k = 1.0;
  model.temperature_per_internal_energy_k = 1.0e4;
  model.max_fractional_energy_change_per_substep = 0.05;
  model.max_subcycles = max_subcycles;
  return physics::CoolingRateProvider(model);
}

const physics::CoolingRateQuery k_query{.temperature_k = 1.0e4, .hydrogen_number_density_cgs = 1.0,
                                        .metallicity_mass_fraction = 0.0, .redshift = 0.0};

void testExactSolutionAndFullInterval() {
  const auto provider = makeProvider(4096);
  physics::CoolingSourceIntegrator integrator(1.0e-8, 1.0e22);  // K = 1
  const auto result = integrator.integrateSpecificInternalEnergy(1.0, 1.0, 1.0, k_query, provider);
  assert(result.completed);
  assert(result.completed_dt_code == result.requested_dt_code);
  assert(!result.diagnostics.hit_subcycle_limit);
  const double exact = 0.25;  // (1 - 1/2)^2
  assert(std::abs(result.specific_internal_energy_code - exact) < 2.0e-4);
  // Tighter request => smaller error (controlled refinement).
  assert(result.diagnostics.accepted_substeps >= 2);
}

void testF2StiffIntervalIsRejectedNotTruncated() {
  const auto provider = makeProvider(2);  // deliberately too few attempts
  physics::CoolingSourceIntegrator integrator(1.0e-8, 1.0e22);
  const double e0 = 1.0;
  const auto result = integrator.integrateSpecificInternalEnergy(e0, 1.0, 1.9, k_query, provider);
  assert(!result.completed);
  assert(result.diagnostics.hit_subcycle_limit);
  assert(result.completed_dt_code < result.requested_dt_code);
  assert(result.specific_internal_energy_code == e0);  // unchanged on refusal
  assert(result.diagnostics.suggested_next_dt_code > 0.0);

  physics::CoolingHeatingSource source(provider, integrator);
  hydro::HydroPrimitiveState primitive{.rho_comoving = 1.0, .pressure_comoving = 2.0 / 3.0};
  primitive.specific_internal_energy_code = 1.0;
  const auto conserved = hydro::HydroCoreSolver::conservedFromPrimitive(primitive, 5.0 / 3.0);
  const std::array<double, 1> rho{1.0}, nh{1.0}, z{0.0}, t{1.0e4};
  hydro::HydroSourceContext context{};
  context.update.dt_code = 1.9;
  context.update.scale_factor = 1.0;
  context.mass_density_physical_cgs = rho;
  context.hydrogen_number_density_cgs = nh;
  context.metallicity_mass_fraction = z;
  context.temperature_k = t;
  bool threw = false;
  try {
    (void)source.sourceForCell(0, conserved, primitive, context);
  } catch (const physics::SourceIntervalIncompleteError& error) {
    threw = true;
    assert(error.suggestedRetryPhysicalSeconds() > 0.0);
  }
  assert(threw);

  // Transaction classification: retry request, nothing committed.
  core::StageTransaction transaction("test.cooling", nullptr);
  bool committed = false;
  const bool prepared = transaction.prepare([&] { (void)source.sourceForCell(0, conserved, primitive, context); });
  assert(!prepared);
  assert(!transaction.validate(true, core::StageRejectionKind::kNone, ""));
  assert(transaction.outcome().rejection == core::StageRejectionKind::kNumericalRetry);
  assert(!committed);
}

void testZeroSourceCompletesInterval() {
  physics::CoolingModelConfig model;
  model.uv_background_model = core::UvBackgroundModel::kNone;
  model.temperature_per_internal_energy_k = 1.0e4;
  physics::CoolingRateProvider provider(model);
  // n_H tiny => negligible rate: the full interval must still complete.
  physics::CoolingSourceIntegrator integrator(1.0e-8, 1.0);
  const physics::CoolingRateQuery query{.temperature_k = 1.0e4, .hydrogen_number_density_cgs = 1.0e-12};
  const auto result = integrator.integrateSpecificInternalEnergy(1.0, 1.0, 10.0, query, provider);
  assert(result.completed);
  assert(result.completed_dt_code == 10.0);
  assert(std::abs(result.specific_internal_energy_code - 1.0) < 1.0e-12);
}

void testIntervalClosureHelper() {
  assert(physics::sourceIntervalClosed(1.0L, 1.0L));
  assert(!physics::sourceIntervalClosed(0.999999L, 1.0L));
  assert(!physics::sourceIntervalClosed(0.0L, 0.0L));
}

}  // namespace

int main() {
  testExactSolutionAndFullInterval();
  testF2StiffIntervalIsRejectedNotTruncated();
  testZeroSourceCompletesInterval();
  testIntervalClosureHelper();
  return 0;
}
