// CHUI-RT-BH-AGN-REF-v1 V05/V06/V07 pure-numerics tests: disk analytic
// solutions, ledgers, supply closures, friction limits, Kerr limits, BZ area
// theorem, merger criteria, energy-consistent kick root, CT magnetic
// injection. Tests call the production functions (no reimplementation of the
// equation under test beyond independent closed forms). NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "cosmosim/core/units.hpp"
#include "cosmosim/physics/agn/agn_injection.hpp"
#include "cosmosim/physics/agn/agn_magnetic.hpp"
#include "cosmosim/physics/black_hole/bh_budget.hpp"
#include "cosmosim/physics/black_hole/bh_disk.hpp"
#include "cosmosim/physics/black_hole/bh_dynamics.hpp"
#include "cosmosim/physics/black_hole/bh_merger.hpp"
#include "cosmosim/physics/black_hole/bh_spin.hpp"
#include "cosmosim/physics/black_hole/bh_supply.hpp"

namespace p = cosmosim::physics;

namespace {

bool close(double a, double b, double rel, double abs_floor = 1.0e-300) {
  return std::abs(a - b) <= std::max(rel * std::max(std::abs(a), std::abs(b)), abs_floor);
}

void testDiskUncappedMatchesExactSolution() {
  p::BhDiskStepInputs in;
  in.disk_mass_code = 2.0;
  in.supply_rate_code = 0.5;
  in.viscous_time_code = 3.0;
  in.wind_loading = 0.5;
  in.dt_code = 1.7;
  const auto r = p::integrateScalarDisk(in);
  const double tau = 3.0 / 1.5;
  const double expected = 2.0 * std::exp(-1.7 / tau) + 0.5 * tau * (1.0 - std::exp(-1.7 / tau));
  CHUI_REQUIRE(close(r.disk_mass_end_code, expected, 1.0e-13));
  CHUI_REQUIRE(close(r.drained_mass_code, 2.0 + 0.5 * 1.7 - expected, 1.0e-12));
  CHUI_REQUIRE(close(r.inner_mass_code + r.wind_mass_code, r.drained_mass_code, 1.0e-14));
  CHUI_REQUIRE(close(r.wind_mass_code, 0.5 * r.inner_mass_code, 1.0e-13));
  // Tiny dt: expm1 path keeps the ledger exact (no cancellation).
  in.dt_code = 1.0e-14;
  const auto t = p::integrateScalarDisk(in);
  CHUI_REQUIRE(close(t.drained_mass_code, 2.0 + 0.5e-14 - t.disk_mass_end_code, 1.0e-6, 1.0e-28));
  CHUI_REQUIRE(t.disk_mass_end_code > 0.0);
}

void testDiskCappedRetainsExcess() {
  p::BhDiskStepInputs in;
  in.disk_mass_code = 10.0;
  in.supply_rate_code = 2.0;
  in.viscous_time_code = 1.0;
  in.wind_loading = 0.0;
  in.inner_rate_cap_code = 1.0;  // M* = 1
  in.dt_code = 2.0;
  const auto r = p::integrateScalarDisk(in);
  // Capped the whole interval: inner = cap dt, the excess stays in the disk.
  CHUI_REQUIRE(close(r.inner_mass_code, 2.0, 1.0e-13));
  CHUI_REQUIRE(close(r.disk_mass_end_code, 10.0 + (2.0 - 1.0) * 2.0, 1.0e-13));
  CHUI_REQUIRE(close(r.capped_time_code, 2.0, 1.0e-13));
  // Draining through the threshold: piecewise analytic, exact ledger.
  in.supply_rate_code = 0.0;
  in.disk_mass_code = 1.5;
  in.dt_code = 5.0;
  const auto d = p::integrateScalarDisk(in);
  CHUI_REQUIRE(d.pieces >= 2U);
  CHUI_REQUIRE(close(d.inner_mass_code + d.wind_mass_code, 1.5 - d.disk_mass_end_code, 1.0e-13));
  CHUI_REQUIRE(close(d.capped_time_code, 0.5, 1.0e-12));
}

void testLedgersAndEddington() {
  p::EddingtonConstantsCode c{1.0, 1.0, 1.0, 10.0};
  p::BhChannelEfficiencies eff;
  eff.radiative = 0.1;
  eff.thermal = 0.01;
  const auto legacy = p::computeBhEmissionBudget(p::BhMassLedger::kLegacyNewtonian, 2.0, 1.0, 5.0, eff, 0.1, c);
  CHUI_REQUIRE(legacy.horizon_mass_increment_code == 2.0);
  CHUI_REQUIRE(legacy.mass_defect_code == 0.0);
  const auto rel = p::computeBhEmissionBudget(p::BhMassLedger::kRelativisticDefect, 2.0, 1.0, 5.0, eff, 0.1, c);
  CHUI_REQUIRE(close(rel.horizon_mass_increment_code, (1.0 - 0.11) * 2.0, 1.0e-15));
  CHUI_REQUIRE(close(rel.mass_defect_code * 100.0, rel.radiative_energy_code + rel.thermal_energy_code, 1.0e-14));
  const double l_edd = p::eddingtonLuminosityCode(5.0, c);
  CHUI_REQUIRE(close(l_edd, 4.0 * cosmosim::core::constants::k_pi * 5.0 * 10.0, 1.0e-15));
  CHUI_REQUIRE(close(p::eddingtonReferenceRateCode(5.0, 0.1, c), l_edd / (0.1 * 100.0), 1.0e-15));
}

void testSupplyClosures() {
  const double bondi = p::bondiHoyleRate(1.0, 2.0, 3.0, 4.0, 3.0, 4.0);
  CHUI_REQUIRE(close(bondi, 4.0 * cosmosim::core::constants::k_pi * 4.0 * 9.0 * 4.0 / 125.0, 1.0e-14));
  CHUI_REQUIRE(p::rosasGuevaraSuppression(1.0, 0.5, 2.0) == 1.0);
  CHUI_REQUIRE(close(p::rosasGuevaraSuppression(1.0, 2.0, 1.0), 0.125, 1.0e-15));
  // Torque reference evaluation at the pivot values (f_d = 1, f_gas = 1).
  const double f0 = 0.31;
  const double torque = p::hopkinsQuataertTorqueRateMsunPerYear(0.1, 1.0, 1.0e8, 1.0e9, 100.0, 1.0e9, 1.0);
  CHUI_REQUIRE(close(torque, 0.1 / (1.0 + f0), 1.0e-14));
  CHUI_REQUIRE(p::hopkinsQuataertTorqueRateMsunPerYear(0.1, 0.0, 1.0e8, 1.0e9, 100.0, 0.0, 1.0) == 0.0);
  // Kernel shape: compact, normalized at the centre, zero at the edge.
  CHUI_REQUIRE(p::bhKernelShape(0.0) == 1.0);
  CHUI_REQUIRE(p::bhKernelShape(1.0) == 0.0);
}

void testFrictionLimits() {
  // Small-X series continuity with the closed form.
  const double a = p::chandrasekharBracket(0.0499999);
  const double b = p::chandrasekharBracket(0.0500001);
  CHUI_REQUIRE(close(a, b, 1.0e-5));
  CHUI_REQUIRE(p::chandrasekharBracket(0.0) == 0.0);
  CHUI_REQUIRE(close(p::chandrasekharBracket(30.0), 1.0, 1.0e-12));
  // Resolved strong-scattering radius: zero correction.
  CHUI_REQUIRE(p::tremmelCoulombLogarithm(1.0, 1.0, 1.0, 0.5) == 0.0);
  p::BhFrictionBackground bg;
  bg.density_proper = 1.0;
  bg.sigma_1d = 1.0;
  bg.samples = 64U;
  p::BhFrictionInputs in;
  in.bh_mass_code = 1.0e-3;
  in.newton_g_code = 1.0;
  in.softening_proper = 1.0;
  const auto res = p::integrateDynamicalFriction({1.0, 0.0, 0.0}, bg, p::BhFrictionBackground{}, in, 1.0);
  CHUI_REQUIRE(res.applied);
  CHUI_REQUIRE(res.velocity_change[0] < 0.0);           // decelerates
  CHUI_REQUIRE(1.0 + res.velocity_change[0] > 0.0);     // never reverses
}

void testKerrLimits() {
  CHUI_REQUIRE(close(p::kerrIscoRadius(0.0), 6.0, 1.0e-14));
  CHUI_REQUIRE(close(p::kerrRadiativeEfficiency(0.0), 1.0 - std::sqrt(8.0 / 9.0), 1.0e-14));  // 0.05719
  CHUI_REQUIRE(close(p::kerrIscoRadius(1.0), 1.0, 1.0e-6));
  CHUI_REQUIRE(close(p::kerrIscoRadius(-1.0), 9.0, 1.0e-12));
  CHUI_REQUIRE(close(p::kerrIscoSpecificAngularMomentum(0.0), 2.0 * std::sqrt(3.0), 1.0e-12));
  CHUI_REQUIRE(close(p::kerrIscoSpecificAngularMomentum(1.0), 2.0 / std::sqrt(3.0), 1.0e-6));
  // Spin-up from zero along the disk axis stays physical.
  p::KerrAccretionInputs in;
  in.mass_code = 1.0;
  in.disk_axis = {0.0, 0.0, 1.0};
  in.inner_mass_code = 3.0;
  in.newton_g_code = 1.0;
  in.speed_of_light_code = 1.0;
  const auto r = p::accreteKerrThinDisk(in);
  const double spin = r.spin_angular_momentum_code[2] / (r.mass_code * r.mass_code);
  CHUI_REQUIRE(spin > 0.0 && spin <= 1.0 + 1.0e-12);
  CHUI_REQUIRE(r.mass_code > 1.0 && r.mass_code < 4.0);
  // Misaligned flows need a warp model: rejected.
  in.spin_angular_momentum_code = {0.5, 0.0, 0.0};
  CHUI_REQUIRE(p::accreteKerrThinDisk(in).misaligned_rejected);
}

void testBlandfordZnajekAreaTheorem() {
  p::BlandfordZnajekInputs in;
  in.mass_code = 1.0;
  in.spin_angular_momentum_code = {0.0, 0.0, 0.9};
  in.magnetic_flux_gaussian_code = 5.0;
  in.newton_g_code = 1.0;
  in.speed_of_light_code = 1.0;
  in.dt_code = 50.0;
  const auto r = p::extractBlandfordZnajek(in);
  CHUI_REQUIRE(r.extracted_energy_code > 0.0);
  CHUI_REQUIRE(r.area_theorem_satisfied);
  CHUI_REQUIRE(close(r.mass_code, 1.0 - r.extracted_energy_code, 1.0e-12));
  in.spin_angular_momentum_code = {0.0, 0.0, 0.0};
  CHUI_REQUIRE(p::extractBlandfordZnajek(in).extracted_energy_code == 0.0);  // zero-spin limit
}

void testMergerCriteria() {
  const auto units = cosmosim::core::makeUnitSystem("kpc", "msun", "km_s");
  const auto frame = p::ProperFrameContext::isolated(units);
  p::BhMergerCriteria crit{2.0, 1.0};
  p::BhMergerCandidate a;
  a.particle_id = 1;
  a.dynamical_mass_code = 1.0;
  a.softening_proper = 1.0;
  p::BhMergerCandidate b = a;
  b.particle_id = 2;
  b.position_stored = {0.5, 0.0, 0.0};
  b.velocity_peculiar = {-0.1, 0.0, 0.0};
  const p::PeriodicBox box{};
  CHUI_REQUIRE(p::bhPairMergeEligible(a, b, crit, frame, box));
  b.velocity_peculiar = {10.0, 0.0, 0.0};  // unbound flyby
  CHUI_REQUIRE(!p::bhPairMergeEligible(a, b, crit, frame, box));
  std::vector<p::BhMergerCandidate> three{a, b, b};
  three[1].velocity_peculiar = {-0.1, 0.0, 0.0};
  three[2].particle_id = 3;
  three[2].velocity_peculiar = {-0.1, 0.0, 0.0};
  const auto pairs = p::selectDisjointMergerPairs(three, crit, frame, box);
  CHUI_REQUIRE(pairs.size() == 1U);  // disjoint: a BH merges at most once per event
}

void testKickQuadraticAndMagneticInjection() {
  double amp = 0.0;
  CHUI_REQUIRE(p::solveEnergyConsistentAmplitude(2.0, -1.0, 3.0, amp));
  CHUI_REQUIRE(close(2.0 * amp * amp - amp, 3.0, 1.0e-13));
  CHUI_REQUIRE(p::solveEnergyConsistentAmplitude(0.0, 2.0, 3.0, amp) && close(amp, 1.5, 1.0e-15));
  CHUI_REQUIRE(!p::solveEnergyConsistentAmplitude(0.0, -1.0, 3.0, amp));
  // Kinetic increment of a momentum kick at fixed mass: v.dp + |dp|^2/2m.
  CHUI_REQUIRE(close(p::cellKineticIncrement(2.0, {2.0, 0.0, 0.0}, 0.0, {1.0, 0.0, 0.0}), 1.0 + 0.25, 1.0e-15));
  p::MagneticPatchGeometry g;
  g.nx = g.ny = g.nz = 8U;
  g.width = {1.0, 1.0, 1.0};
  const auto inj = p::buildToroidalMagneticInjection(g, {4.0, 4.0, 4.0}, {0.0, 0.0, 1.0}, 3.0);
  CHUI_REQUIRE(p::magneticDiscreteDivergenceMax(g, inj) < 1.0e-14);
  const auto terms = p::magneticInjectionEnergyTerms(g, inj, inj);  // existing = injection pattern
  CHUI_REQUIRE(terms.quadratic > 0.0 && close(terms.linear, 2.0 * terms.quadratic, 1.0e-12));
  double s = 0.0;
  CHUI_REQUIRE(p::solveMagneticInjectionAmplitude(terms, 1.0, s));
  CHUI_REQUIRE(close(s * terms.linear + s * s * terms.quadratic, 1.0, 1.0e-12));
}

}  // namespace

int main() {
  testDiskUncappedMatchesExactSolution();
  testDiskCappedRetainsExcess();
  testLedgersAndEddington();
  testSupplyClosures();
  testFrictionLimits();
  testKerrLimits();
  testBlandfordZnajekAreaTheorem();
  testMergerCriteria();
  testKickQuadraticAndMagneticInjection();
  return 0;
}
