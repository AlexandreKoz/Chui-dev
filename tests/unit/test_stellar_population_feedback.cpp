// CHUI Part 2 — F5 chronology, SSP model, event replay, mechanical feedback
// (NOT RUN — USER-RESERVED).
//
// Invariants:
//  * F5: a star born at t_end (birth tick == interval tick) has zero elapsed
//    age in that interval; any other star gets the full interval;
//  * IMF normalized to unit mass; cumulative CC count per Msun for a Chabrier
//    IMF with 8-40 Msun is ~0.01 (literature range 0.009-0.012);
//  * returned mass never exceeds formed mass; AGB+CC return over a Hubble
//    time is O(0.3-0.45) of the formed mass for this IMF/IFMR;
//  * Ia cumulative number reaches N_Ia at t_max and is zero before t_min;
//  * step subdivision (Part 2.5 review): the additivity check is only
//    meaningful if NO non-monotone clamp fired (massiveReturn clamps a
//    negative quadrature difference element-wise); the test asserts the
//    clamp counter stayed zero, that every sub-interval keeps element
//    closure, and that the cumulative return is non-decreasing for the
//    total and every element on a fine age grid;
//  * singleton yield-table axes report clamping when queried off-grid;
//  * event chronology is invariant under step subdivision and restart
//    (same thresholds => same realized ordinals);
//  * patch 52 (T02): the physical age of each realized event (threshold
//    crossing of Lambda_c) is identical under step subdivision, lies in its
//    own interval, and nextStellarEventAgeYears predicts exactly the first
//    event the stage then realizes (so the time authority's publication-lag
//    clip is consistent with the chronology);
//  * mechanical feedback: mass, every element, and vector momentum (lab
//    frame: m_ej v_star) conserved exactly; injected kinetic + thermal energy
//    equals the budget; asymmetric stencil still has zero net star-frame
//    momentum; unresolved cells use the KO15 terminal momentum.
#include <algorithm>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/source_interval.hpp"
#include "cosmosim/physics/stellar_population.hpp"
#include "cosmosim/physics/supernova_feedback.hpp"

namespace {

using namespace cosmosim;

void testF5Chronology() {
  assert(physics::stellarElapsedYearsInInterval(42U, 42U, 1.0e6) == 0.0);
  assert(physics::stellarElapsedYearsInInterval(41U, 42U, 1.0e6) == 1.0e6);
  assert(physics::stellarElapsedYearsInInterval(42U, std::numeric_limits<std::uint64_t>::max(), 1.0e6) == 1.0e6);
}

void testSspBudgets() {
  const physics::SspModel ssp(physics::SspConfig{});
  assert(std::abs(ssp.massBetween(0.1, 100.0) - 1.0) < 1.0e-6);
  const double n_cc = ssp.numberBetween(8.0, 40.0);
  assert(n_cc > 0.008 && n_cc < 0.013);
  std::array<double, physics::k_chemical_element_count> x{};
  x[0] = 0.74; x[1] = 0.25; x[4] = 0.01;
  const auto agb = ssp.massiveReturn(0, 0.0, 1.37e10, x);
  const auto cc = ssp.massiveReturn(1, 0.0, 1.37e10, x);
  const double total = agb.mass_fraction + cc.mass_fraction;
  assert(total > 0.25 && total < 0.5 && total < 1.0);
  double el = 0.0;
  for (std::size_t e = 0; e < x.size(); ++e) el += agb.element_fraction[e] + cc.element_fraction[e];
  assert(std::abs(el - total) < 1.0e-10);  // element closure
  // A young population (< shortest lifetime) returns nothing.
  assert(ssp.massiveReturn(1, 0.0, 1.0e5, x).mass_fraction == 0.0);
  assert(ssp.iaCumulativeNumber(1.0e7) == 0.0);
  assert(std::abs(ssp.iaCumulativeNumber(2.0e10) - 1.3e-3) < 1.0e-12);
  // Subdivision: integral over [0,t1] + [t1,t2] equals [0,t2].
  const std::uint64_t clamps_before = ssp.nonmonotone_return_clamps;
  const auto a = ssp.massiveReturn(1, 0.0, 1.0e7, x);
  const auto b = ssp.massiveReturn(1, 1.0e7, 4.0e7, x);
  const auto c = ssp.massiveReturn(1, 0.0, 4.0e7, x);
  // Telescoping holds only when no difference was clamped; otherwise the
  // equality below could pass or fail for the wrong reason.
  assert(ssp.nonmonotone_return_clamps == clamps_before);
  assert(c.mass_fraction > 0.0 && a.mass_fraction >= 0.0 && b.mass_fraction > 0.0);
  // Exactly additive (cumulative-difference construction), up to round-off.
  assert(std::abs(a.mass_fraction + b.mass_fraction - c.mass_fraction) < 1.0e-12 * c.mass_fraction);
  for (const auto* part : {&a, &b, &c}) {
    double sum = 0.0;
    for (std::size_t e = 0; e < x.size(); ++e) sum += part->element_fraction[e];
    assert(std::abs(sum - part->mass_fraction) <= 1.0e-12 * std::max(part->mass_fraction, 1.0e-300));
  }
  for (std::size_t e = 0; e < x.size(); ++e) {
    assert(std::abs(a.element_fraction[e] + b.element_fraction[e] - c.element_fraction[e]) <=
           1.0e-12 * std::max(c.mass_fraction, 1.0e-300));
  }
  // Monotonicity of the cumulative return on a fine log-age grid (both
  // channels, total and per element): the exact function is
  // non-decreasing; allow only double round-off of the running sum.
  for (std::size_t channel = 0; channel < 2; ++channel) {
    physics::SspModel::ChannelReturn previous{};
    for (int k = 0; k <= 400; ++k) {
      const double age = std::pow(10.0, 6.0 + 4.15 * k / 400.0);
      const auto cum = ssp.cumulativeMassiveReturn(channel, age, x);
      assert(cum.mass_fraction >= previous.mass_fraction - 1.0e-14);
      for (std::size_t e = 0; e < x.size(); ++e) {
        assert(cum.element_fraction[e] >= previous.element_fraction[e] - 1.0e-14);
      }
      previous = cum;
    }
  }
}

void testEventReplayIsSubdivisionAndRestartInvariant() {
  std::array<double, physics::k_chemical_element_count> x{};
  x[0] = 0.75; x[1] = 0.25;
  auto one = physics::makeStellarPopulationRecord(x, 0.0, 0xabcdefULL);
  auto many = one;
  const auto big = physics::realizeEvents(one, 0U, 37.5);
  std::uint32_t total = 0;
  for (double lam = 0.5; lam <= 37.5; lam += 0.5) total += physics::realizeEvents(many, 0U, lam).count;
  assert(big.count == total);
  assert(one.cc_realized == many.cc_realized && one.cc_next_threshold == many.cc_next_threshold);
  // Restart: a copied record continues identically; a rejected stage that
  // never stores its copy consumes nothing.
  auto restored = one;
  assert(physics::realizeEvents(restored, 0U, 60.0).count == physics::realizeEvents(one, 0U, 60.0).count);
  assert(big.count > 20U && big.count < 60U);  // Poisson(37.5)
}

core::SimulationState makeGas() {
  core::SimulationState state;
  state.resizeCells(27);
  std::size_t row = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k, ++row) {
        state.gas_cells.gas_cell_id[row] = 1000U + row;
        state.cells.center_x_comoving[row] = i + 0.5;
        state.cells.center_y_comoving[row] = j + 0.5;
        state.cells.center_z_comoving[row] = k + 0.5;
        state.cells.mass_code[row] = 1.0;
        state.gas_cells.density_code[row] = 1.0;
        state.gas_cells.internal_energy_code[row] = 1.0;
        state.gas_cells.velocity_x_peculiar[row] = 0.1 * i;
      }
  return state;
}

void testMechanicalFeedbackConservation() {
  auto state = makeGas();
  physics::FeedbackUnits units;  // code == cgs for the test
  units.mass_code_to_g = physics::k_solar_mass_g;
  units.velocity_code_to_cm_s = 1.0e5;
  units.length_code_to_cm = 3.0857e18;
  units.density_code_to_cgs = 1.0e-24;
  physics::MechanicalFeedbackDeposition dep({.neighbor_count = 8}, units);
  std::vector<std::uint32_t> rows(27);
  for (std::uint32_t r = 0; r < 27; ++r) rows[r] = r;
  dep.setCandidates(state, rows);
  physics::FeedbackSourceBudget src;
  src.source_id = 7;
  src.position_comoving = {1.3, 1.6, 1.4};  // asymmetric stencil
  src.velocity_peculiar = {2.0, -1.0, 0.5};
  src.ejecta_mass_code = 10.0;
  src.ejecta_elements_code[0] = 7.0; src.ejecta_elements_code[1] = 2.5; src.ejecta_elements_code[4] = 0.5;
  src.supernova_count = 1.0;
  src.supernova_energy_erg = 1.0e51;
  double mass0 = 0.0, e0 = 0.0;
  std::array<double, 3> p0{};
  for (std::size_t r = 0; r < 27; ++r) {
    const double m = state.cells.mass_code[r];
    const double vx = state.gas_cells.velocity_x_peculiar[r];
    mass0 += m;
    p0[0] += m * vx;
    e0 += m * state.gas_cells.internal_energy_code[r] + 0.5 * m * vx * vx;
  }
  dep.stage(state, src);
  (void)dep.commit(state);
  double mass1 = 0.0, e1 = 0.0;
  std::array<double, 3> p1{};
  for (std::size_t r = 0; r < 27; ++r) {
    const double m = state.cells.mass_code[r];
    const double vx = state.gas_cells.velocity_x_peculiar[r], vy = state.gas_cells.velocity_y_peculiar[r],
                 vz = state.gas_cells.velocity_z_peculiar[r];
    mass1 += m;
    p1[0] += m * vx; p1[1] += m * vy; p1[2] += m * vz;
    e1 += m * state.gas_cells.internal_energy_code[r] + 0.5 * m * (vx * vx + vy * vy + vz * vz);
  }
  const auto& l = dep.ledger();
  assert(std::abs(mass1 - mass0 - 10.0) < 1.0e-10);
  for (int c = 0; c < 3; ++c) assert(std::abs(p1[c] - p0[c] - 10.0 * src.velocity_peculiar[c]) < 1.0e-8);
  const double e_code = 1.0e51 / (units.mass_code_to_g * 1.0e10);
  const double budget = e_code + 0.5 * 10.0 * (4.0 + 1.0 + 0.25);
  assert(std::abs((e1 - e0) - budget) < 1.0e-8 * budget);
  assert(std::abs(l.kinetic_energy_injected_code + l.thermal_energy_injected_code - budget) < 1.0e-8 * budget);
  assert(l.resolved_cells + l.unresolved_cells == 8U);
  assert(l.radial_momentum_delivered_code <= l.radial_momentum_requested_code + 1.0e-9);
}

std::string yieldTable(double h_net) {
  // Synthetic UNIT-TEST fixture (not astrophysical data): one Z, two masses.
  std::string t = "# chui_ssp_yield_table_v1\n# citation = unit-test fixture (not astrophysical data)\n"
                  "# license = CC0 test fixture\n# net_yield_sum_tolerance = 1e-3\n"
                  "elements H He C N O Ne Mg Si S Ca Fe\n";
  const std::string row_ok = " 2.0 " + std::to_string(h_net) + " 0.01 0.02 0 0.05 0 0 0 0 0 0\n";
  t += "channel agb\nmetallicity 0.02\nmass 1 5\nrow 0 0" + row_ok + "row 0 1" + row_ok;
  t += "channel cc\nmetallicity 0.02\nmass 10 30\nrow 0 0" + row_ok + "row 0 1" + row_ok;
  t += "channel ia\nmetallicity 0.02\nrow 0 1.4 0 0 0.05 0 0.15 0 0.01 0.15 0.08 0.01 0.95\n";
  return t;
}

void testYieldTableIsNeverRebalanced() {
  // Balanced: -0.08 H balances +0.08 metals -> accepted, values untouched.
  const auto table = physics::SspYieldTable::loadFromString(yieldTable(-0.08), "fixture");
  bool clamped = false;
  const auto row = table.massive(0, 3.0, 0.02, clamped);
  assert(std::abs(row.net_yield_msun[0] + 0.08) < 1.0e-15);
  assert(table.max_net_yield_residual < 1.0e-12);
  // Part 2.5: the fixture has ONE metallicity (0.02). On-grid queries are
  // not clamps; off-grid metallicities use that row and are REPORTED.
  assert(!clamped);
  bool off_grid = false;
  (void)table.massive(0, 3.0, 0.001, off_grid);
  assert(off_grid);
  bool ia_on = false, ia_off = false;
  (void)table.snia(0.02, ia_on);
  (void)table.snia(0.001, ia_off);
  assert(!ia_on && ia_off);
  // Unbalanced by 0.08 / 2.0 = 4% > tolerance: refused, never moved into H.
  bool threw = false;
  try { (void)physics::SspYieldTable::loadFromString(yieldTable(0.0), "fixture"); }
  catch (const physics::SourceModelUnavailableError&) { threw = true; }
  assert(threw);
}

void testTerminalMomentumScaling() {
  const double p1 = physics::kimOstriker2015TerminalMomentumMsunKms(1.0e51, 1.0);
  assert(std::abs(p1 - 2.17e5) < 1.0);
  assert(physics::kimOstriker2015TerminalMomentumMsunKms(1.0e51, 100.0) < p1);
}

// G3.2: fixed-node cumulative return. A 200-way subdivision of [1 Myr, 13.7
// Gyr] (as a population crossing many steps / AMR substeps) sums to the
// single-interval return to round-off, with NO non-monotone clamp in either
// channel; every sub-interval return is non-negative per element.
void testSspFixedNodeSubdivisionIsMonotone() {
  const physics::SspModel ssp(physics::SspConfig{});
  std::array<double, physics::k_chemical_element_count> x{};
  x[0] = 0.72; x[1] = 0.26; x[4] = 0.012; x[2] = 0.003; x[3] = 0.001; x[6] = 0.004;
  const std::uint64_t clamps_before = ssp.nonmonotone_return_clamps;
  for (std::size_t channel = 0; channel < 2; ++channel) {
    const double t0 = 1.0e6, t1 = 1.37e10;
    const auto whole = ssp.massiveReturn(channel, t0, t1, x);
    double sum = 0.0;
    std::array<double, physics::k_chemical_element_count> sum_e{};
    double previous_age = t0;
    for (int k = 1; k <= 200; ++k) {
      const double age = k == 200 ? t1 : t0 * std::pow(t1 / t0, static_cast<double>(k) / 200.0);
      const auto part = ssp.massiveReturn(channel, previous_age, age, x);
      assert(part.mass_fraction >= 0.0);
      for (std::size_t e = 0; e < x.size(); ++e) {
        assert(part.element_fraction[e] >= 0.0);
        sum_e[e] += part.element_fraction[e];
      }
      sum += part.mass_fraction;
      previous_age = age;
    }
    assert(std::abs(sum - whole.mass_fraction) <= 1.0e-11 * std::max(whole.mass_fraction, 1.0e-300));
    for (std::size_t e = 0; e < x.size(); ++e) {
      assert(std::abs(sum_e[e] - whole.element_fraction[e]) <= 1.0e-11 * std::max(whole.mass_fraction, 1.0e-300));
    }
  }
  assert(ssp.nonmonotone_return_clamps == clamps_before);
}

// R2 (A04): the event forecast is the actual chronology. Scanning a 1000 Msun
// population in 1 Myr steps, forecastStellarEvent(age, 1 Myr) is true exactly
// in the steps where realizeEvents (the source-stage law) then realizes >= 1
// event, and the forecast never mutates the record.
void testEventForecastMatchesChronology() {
  const physics::SspModel ssp(physics::SspConfig{});
  std::array<double, physics::k_chemical_element_count> x{};
  x[0] = 0.74; x[1] = 0.25; x[4] = 0.01;
  physics::StellarPopulationRecord record = physics::makeStellarPopulationRecord(x, 0.0, 0x5eedULL);
  const double birth = 1000.0;
  double z = 0.0;
  for (std::size_t e = 2; e < x.size(); ++e) z += x[e];
  int agreements = 0;
  for (int step = 0; step < 60; ++step) {
    const double age = 1.0e6 * step;
    const physics::StellarPopulationRecord before = record;
    const bool forecast = physics::forecastStellarEvent(ssp, record, age, 1.0e6, birth, true);
    assert(std::memcmp(&before, &record, sizeof(record)) == 0);  // nothing consumed
    const auto cc = physics::realizeEvents(record, 0U, ssp.ccCumulativeNumber(age + 1.0e6, z) * birth);
    const auto ia = physics::realizeEvents(record, 1U, ssp.iaCumulativeNumber(age + 1.0e6) * birth);
    assert(forecast == (cc.count + ia.count > 0U));
    if (forecast) ++agreements;
  }
  assert(agreements > 0);  // the 3-40 Myr CC window of a 1000 Msun population has events
  // Expectation mode: >= 0.5 expected events in the window (heuristic rule).
  physics::StellarPopulationRecord fresh = physics::makeStellarPopulationRecord(x, 0.0, 1U);
  assert(!physics::forecastStellarEvent(ssp, fresh, 0.0, 1.0e5, birth, false));      // before the first SN
  assert(physics::forecastStellarEvent(ssp, fresh, 5.0e6, 1.0e7, birth, false));
  fresh.initialized = 0U;
  assert(!physics::forecastStellarEvent(ssp, fresh, 5.0e6, 1.0e7, birth, true));     // no record => no forecast
}

}  // namespace

void testEventPhysicalTimeIsSubdivisionInvariant() {
  const physics::SspModel ssp(physics::SspConfig{});
  std::array<double, physics::k_chemical_element_count> x{};
  x[0] = 0.74; x[1] = 0.25; x[4] = 0.01;
  double z = 0.0;
  for (std::size_t e = 2; e < x.size(); ++e) z += x[e];
  const double birth = 1000.0;
  const auto collect = [&](double step_years) {
    physics::StellarPopulationRecord record = physics::makeStellarPopulationRecord(x, 0.0, 0xabcdULL);
    std::vector<double> ages;
    for (double age_b = 0.0; age_b < 4.0e7; age_b += step_years) {
      const double age_e = age_b + step_years;
      const double predicted = physics::nextStellarEventAgeYears(ssp, record, age_b, birth);
      const physics::StellarPopulationRecord before = record;
      const auto cc = physics::realizeEvents(record, 0U, ssp.ccCumulativeNumber(age_e, z) * birth);
      const auto chron = physics::realizedEventChronology(ssp, before, 0U, cc.count, z, birth, age_b, age_e);
      assert(chron.catch_up_count == 0U);
      if (cc.count > 0U) {
        assert(chron.earliest_age_years >= age_b && chron.latest_age_years <= age_e);
        // Ia cannot fire before snia_delay_min, so the next event is CC here.
        if (age_e < physics::SspConfig{}.snia_delay_min_years) {
          assert(std::abs(predicted - chron.earliest_age_years) <= 1.0e-9 * age_e);
        }
        ages.push_back(chron.sum_age_years / static_cast<double>(cc.count));
      }
      (void)physics::realizeEvents(record, 1U, ssp.iaCumulativeNumber(age_e) * birth);
    }
    return std::make_pair(record.cc_realized, ages.empty() ? 0.0 : ages.front());
  };
  const auto coarse = collect(1.0e6);
  const auto fine = collect(2.5e5);
  assert(coarse.first == fine.first && coarse.first > 0U);  // same realized ordinals
  // The first event's physical age is the same (it is a single event in
  // its interval for both cadences here or the mean of identical ages).
  assert(std::abs(coarse.second - fine.second) <= 1.0e6);
  // Direct check: the same threshold has one age whatever the bracket.
  const double threshold = physics::eventExponentialDraw(0xabcdULL, 0U, 0U);
  const double a1 = physics::stellarEventAgeYears(ssp, 0U, threshold, z, birth, 0.0, 4.0e7);
  const double a2 = physics::stellarEventAgeYears(ssp, 0U, threshold, z, birth, 1.0e6, 1.0e8);
  if (a1 > 1.0e6) assert(std::abs(a1 - a2) <= 1.0e-9 * a1);
  assert(!std::isfinite(physics::stellarEventAgeYears(ssp, 0U, 1.0e30, z, birth, 0.0, 1.0e10)));
}

int main() {
  testF5Chronology();
  testSspBudgets();
  testSspFixedNodeSubdivisionIsMonotone();
  testEventForecastMatchesChronology();
  testEventPhysicalTimeIsSubdivisionInvariant();
  testEventReplayIsSubdivisionAndRestartInvariant();
  testMechanicalFeedbackConservation();
  testTerminalMomentumScaling();
  testYieldTableIsNeverRebalanced();
  return 0;
}
