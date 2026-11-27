// CHUI-EXT-PHYS-REF-v1 gate VCR (CR01-CR08): advection conservation and gas
// ledger, field-aligned diffusion (row conservation, variance growth 2 kappa t),
// Sharma-Hammett monotonicity at oblique field, RKL2 vs explicit, streaming
// heating balance, two-moment diffusive limit vs one-moment and exact energy
// exchange, power-law reconstruction, adiabatic momentum-space work
// (momentum-space part dE/dt = -(1/3) E div u, number conserved),
// synchrotron cooling vs the analytic solution, injection normalization,
// loss-provider refusal. Expected values are analytic. NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/cosmic_rays/cosmic_ray_sources.hpp"
#include "cosmosim/cosmic_rays/cosmic_ray_spectrum.hpp"
#include "cosmosim/cosmic_rays/cosmic_ray_transport.hpp"

namespace cr = cosmosim::cosmic_rays;

namespace {

bool closeRel(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

struct Medium {
  std::vector<double> ux, uy, uz, bx, by, bz, vs;
  explicit Medium(std::size_t n) : ux(n, 0.0), uy(n, 0.0), uz(n, 0.0), bx(n, 1.0), by(n, 0.0), bz(n, 0.0), vs() {}
  cr::CrMediumView view() const { return {ux, uy, uz, bx, by, bz, vs}; }
};

double sum(const std::vector<double>& v) {
  double s = 0.0;
  for (const double x : v) s += x;
  return s;
}

// EXT-CR01-T01: periodic advection conserves energy; ledger residual vanishes for uniform u.
void testAdvection() {
  cr::CrGrid grid;
  grid.n = {64, 1, 1};
  grid.dx = {1.0 / 64.0, 1.0, 1.0};
  Medium m(grid.cells());
  std::fill(m.ux.begin(), m.ux.end(), 1.0);
  std::vector<double> e(grid.cells(), 1.0);
  for (std::size_t i = 16; i < 32; ++i) e[i] = 2.0;
  const double total0 = sum(e);
  cr::CrGasExchange gas;
  cr::CrStepLedger ledger;
  const double dt = cr::advectionTimestep(grid, m.view(), 0.5);
  for (int s = 0; s < 20; ++s) cr::advectGrey(grid, e, m.view(), 4.0 / 3.0, dt, gas, ledger);
  CHUI_REQUIRE(closeRel(sum(e), total0, 1.0e-14));
  CHUI_REQUIRE(std::abs(ledger.residual()) < 1.0e-12 * total0 * grid.cellVolume());
  CHUI_REQUIRE(*std::min_element(e.begin(), e.end()) >= 1.0 - 1.0e-14);  // donor cell: no undershoot
}

// EXT-CR02-T01: field along x, kappa_perp = 0: each row conserves its energy,
// and the x variance of a Gaussian grows by 2 kappa t.
void testFieldAlignedDiffusion() {
  cr::CrGrid grid;
  grid.n = {128, 8, 1};
  grid.dx = {1.0 / 128.0, 1.0 / 8.0, 1.0};
  Medium m(grid.cells());
  cr::CrDiffusionModel model;
  model.kappa_parallel = 1.0e-3;
  std::vector<double> e(grid.cells());
  const double sigma0 = 0.05;
  for (std::size_t j = 0; j < 8; ++j) {
    for (std::size_t i = 0; i < 128; ++i) {
      const double x = (i + 0.5) / 128.0 - 0.5;
      e[grid.index(i, j, 0)] = (1.0 + j) * std::exp(-x * x / (2.0 * sigma0 * sigma0));
    }
  }
  std::vector<double> row0(8, 0.0);
  for (std::size_t j = 0; j < 8; ++j) {
    for (std::size_t i = 0; i < 128; ++i) row0[j] += e[grid.index(i, j, 0)];
  }
  const double dt = cr::diffusionExplicitTimestep(grid, model, 0.4);
  const double t_end = 0.5;
  const int steps = static_cast<int>(std::ceil(t_end / dt));
  cr::CrStepLedger ledger;
  for (int s = 0; s < steps; ++s) cr::diffuseExplicit(grid, e, m.view(), model, t_end / steps, ledger);
  for (std::size_t j = 0; j < 8; ++j) {
    double row = 0.0;
    double var = 0.0;
    for (std::size_t i = 0; i < 128; ++i) {
      const double x = (i + 0.5) / 128.0 - 0.5;
      row += e[grid.index(i, j, 0)];
      var += x * x * e[grid.index(i, j, 0)];
    }
    CHUI_REQUIRE(closeRel(row, row0[j], 1.0e-12));
    CHUI_REQUIRE(std::abs(var / row - (sigma0 * sigma0 + 2.0 * model.kappa_parallel * t_end)) < 2.0e-5);
  }
}

// EXT-CR02-T02: oblique field, limited transverse gradients: no new extrema.
void testMonotonicity() {
  cr::CrGrid grid;
  grid.n = {48, 48, 1};
  grid.dx = {1.0 / 48.0, 1.0 / 48.0, 1.0};
  Medium m(grid.cells());
  for (std::size_t c = 0; c < grid.cells(); ++c) {
    m.bx[c] = std::cos(0.5);
    m.by[c] = std::sin(0.5);
  }
  cr::CrDiffusionModel model;
  model.kappa_parallel = 1.0;
  std::vector<double> e(grid.cells(), 1.0);
  for (std::size_t j = 20; j < 28; ++j) {
    for (std::size_t i = 20; i < 28; ++i) e[grid.index(i, j, 0)] = 10.0;
  }
  const double dt = cr::diffusionExplicitTimestep(grid, model, 0.4);
  cr::CrStepLedger ledger;
  const double total0 = sum(e);
  for (int s = 0; s < 100; ++s) cr::diffuseExplicit(grid, e, m.view(), model, dt, ledger);
  CHUI_REQUIRE(*std::min_element(e.begin(), e.end()) >= 1.0 - 1.0e-12);
  CHUI_REQUIRE(*std::max_element(e.begin(), e.end()) <= 10.0 + 1.0e-12);
  CHUI_REQUIRE(closeRel(sum(e), total0, 1.0e-12));
  // EXT-CR02-T03: one RKL2 super-step of 20 explicit steps vs 20 explicit steps.
  std::vector<double> a = e;
  std::vector<double> b = e;
  for (int s = 0; s < 20; ++s) cr::diffuseExplicit(grid, a, m.view(), model, dt, ledger);
  const std::int32_t stages = cr::diffuseRkl2(grid, b, m.view(), model, 20.0 * dt, 0.4, ledger);
  CHUI_REQUIRE(stages >= 2);
  double diff = 0.0;
  double norm = 0.0;
  for (std::size_t c = 0; c < a.size(); ++c) {
    diff += (a[c] - b[c]) * (a[c] - b[c]);
    norm += a[c] * a[c];
  }
  CHUI_REQUIRE(std::sqrt(diff / norm) < 1.0e-3);
  CHUI_REQUIRE(closeRel(sum(b), sum(a), 1.0e-12));
}

// EXT-CR02-T04: streaming heats the gas by exactly what the CR lose (periodic).
void testStreaming() {
  cr::CrGrid grid;
  grid.n = {64, 1, 1};
  grid.dx = {1.0 / 64.0, 1.0, 1.0};
  Medium m(grid.cells());
  m.vs.assign(grid.cells(), 0.1);
  std::vector<double> e(grid.cells());
  for (std::size_t i = 0; i < 64; ++i) e[i] = 2.0 + std::sin(2.0 * 3.141592653589793 * (i + 0.5) / 64.0);
  cr::CrGasExchange gas;
  cr::CrStepLedger ledger;
  const double dt = cr::streamingTimestep(grid, m.view(), 0.02, 0.4);
  for (int s = 0; s < 50; ++s) cr::streamGrey(grid, e, m.view(), 4.0 / 3.0, 0.02, dt, gas, ledger);
  CHUI_REQUIRE(ledger.streaming_heating > 0.0);
  CHUI_REQUIRE(std::abs(ledger.cr_energy_change + ledger.streaming_heating) < 1.0e-12 * ledger.streaming_heating);
  for (const double h : gas.energy) CHUI_REQUIRE(h >= 0.0);
}

// EXT-CR03-T01: two-moment in the stiff (diffusive) limit reproduces one-moment
// diffusion; energy exchange with the gas is exact.
void testTwoMoment() {
  cr::CrGrid grid;
  grid.n = {128, 1, 1};
  grid.dx = {1.0 / 128.0, 1.0, 1.0};
  Medium m(grid.cells());
  cr::CrTwoMomentSettings settings;
  settings.reduced_speed = 50.0;  // c_r dx / kappa = 39 >> 1: stiff regime
  settings.diffusion.kappa_parallel = 1.0e-2;
  cr::CrTwoMomentState state;
  state.e.resize(grid.cells());
  state.f.assign(grid.cells(), {0.0, 0.0, 0.0});
  std::vector<double> reference(grid.cells());
  for (std::size_t i = 0; i < 128; ++i) {
    const double x = (i + 0.5) / 128.0 - 0.5;
    state.e[i] = 1.0 + std::exp(-x * x / (2.0 * 0.05 * 0.05));
    reference[i] = state.e[i];
  }
  cr::CrTwoMomentReport report;
  const double t_end = 0.05;
  CHUI_REQUIRE(cr::advanceTwoMoment(grid, state, m.view(), settings, t_end, report));
  cr::CrStepLedger ledger;
  cr::CrDiffusionModel model = settings.diffusion;
  const double dt = cr::diffusionExplicitTimestep(grid, model, 0.4);
  const int steps = static_cast<int>(std::ceil(t_end / dt));
  for (int s = 0; s < steps; ++s) cr::diffuseExplicit(grid, reference, m.view(), model, t_end / steps, ledger);
  double diff = 0.0;
  double norm = 0.0;
  for (std::size_t i = 0; i < 128; ++i) {
    diff += (state.e[i] - reference[i]) * (state.e[i] - reference[i]);
    norm += (reference[i] - 1.0) * (reference[i] - 1.0);
  }
  CHUI_REQUIRE(std::sqrt(diff / norm) < 0.03);  // relaxation-time lag O(tau / t)
  CHUI_REQUIRE(std::abs(report.cr_energy_change + report.gas_energy_gain) < 1.0e-12);
}

// EXT-CR04-T01: power-law reconstruction recovers the slope exactly.
void testReconstruction() {
  const double rest = cr::crSpeciesRestEnergyJ(cr::CrSpecies::kProton);
  for (const double q : {-2.0, 0.0, 2.2, 4.5, 8.0}) {
    const double n = cr::powerLawNumber(1.0, 3.0, q, 2.5);
    const double e = cr::powerLawEnergy(1.0, 3.0, q, 2.5, rest);
    const auto bin = cr::reconstructPowerLaw(1.0, 3.0, n, e, rest, -10.0, 15.0);
    CHUI_REQUIRE(std::abs(bin.q - q) < 1.0e-9 && closeRel(bin.f0, 2.5, 1.0e-9) && !bin.clamped);
  }
  // q = 3 limit of the number integral: 4 pi f0 lo^3 ln(hi/lo).
  CHUI_REQUIRE(closeRel(cr::powerLawNumber(1.0, 2.0, 3.0, 1.0), 4.0 * 3.141592653589793 * std::log(2.0), 1.0e-12));
}

// EXT-CR04-T02: adiabatic compression of an ultrarelativistic spectrum well inside
// the grid. The MOMENTUM-SPACE part alone (per fixed volume) conserves number
// and gives dE/dt = A E with A = -div u / 3; the spatial number compression
// (3 A n) is separate, so the total is -(4/3) E div u as for a gamma = 4/3 gas.
void testAdiabaticMomentumSpace() {
  cr::CrSpectralSettings settings;
  settings.species = cr::CrSpecies::kElectron;
  cr::CrSpectralState state(cr::CrMomentumGrid::logarithmic(1.0e1, 1.0e7, 24), settings, 1, nullptr);
  cr::CrInjectionLedger injection;
  cr::injectCosmicRays(state, 0, 1.0, {2.2, 1.0e3, 1.0e5}, injection);
  CHUI_REQUIRE(closeRel(state.totalEnergy(0), 1.0, 1.0e-12));
  const double n0 = state.totalNumber(0);
  const double e0 = state.totalEnergy(0);
  cr::CrMomentumRate rate;
  rate.adiabatic_coefficient = 1.0e-3;  // compression: -div u / 3 > 0
  const double dt = cr::momentumTimestep(state.grid(), rate, 0.005);  // forward-Euler error << tolerance
  const int steps = static_cast<int>(std::ceil(0.9 / (rate.adiabatic_coefficient * dt)));
  cr::CrSpectralLedger ledger;
  for (int s = 0; s < steps; ++s) CHUI_REQUIRE(cr::advanceMomentumSpace(state, 0, rate, dt, ledger));
  CHUI_REQUIRE(closeRel(state.totalNumber(0), n0, 1.0e-12));  // no edge losses: spectrum interior
  const double t = steps * dt;
  // T ~ p m c^2 for p~ >= 1e3 (relative 1e-3 rest-mass correction).
  CHUI_REQUIRE(closeRel(state.totalEnergy(0), e0 * std::exp(rate.adiabatic_coefficient * t), 2.0e-3));
  CHUI_REQUIRE(closeRel(state.totalEnergy(0) - e0, ledger.work, 1.0e-10));
}

// EXT-CR06-T01: synchrotron cooling of a narrow spectrum vs p(t) = p0 / (1 + C p0 t).
void testSynchrotronCooling() {
  const double coefficient = cr::electronRadiativeLossCoefficient(1.0e-9, 0.0);  // 10 microgauss
  CHUI_REQUIRE(coefficient > 0.0);
  cr::CrSpectralSettings settings;
  settings.species = cr::CrSpecies::kElectron;
  cr::CrSpectralState state(cr::CrMomentumGrid::logarithmic(1.0e2, 1.0e6, 80), settings, 1, nullptr);
  cr::CrInjectionLedger injection;
  cr::injectCosmicRays(state, 0, 1.0, {2.0, 4.0e4, 5.0e4}, injection);
  cr::CrMomentumRate rate;
  rate.radiative_coefficient = coefficient;
  const double p0 = 4.5e4;
  const double t_cool = 1.0 / (coefficient * p0);
  const double dt = cr::momentumTimestep(state.grid(), rate, 0.2);
  const int steps = static_cast<int>(std::ceil(t_cool / dt));
  cr::CrSpectralLedger ledger;
  for (int s = 0; s < steps; ++s) CHUI_REQUIRE(cr::advanceMomentumSpace(state, 0, rate, t_cool / steps, ledger));
  // Mean momentum weighted by number ~ p0 / 2 after one cooling time (first-order scheme: 10% band).
  double num = 0.0;
  double pn = 0.0;
  for (std::size_t b = 0; b < state.grid().bins(); ++b) {
    const double pc = std::sqrt(state.grid().edges[b] * state.grid().edges[b + 1]);
    num += state.number(b)[0];
    pn += pc * state.number(b)[0];
  }
  CHUI_REQUIRE(std::abs(pn / num - p0 / 2.0) < 0.1 * p0 / 2.0);
  CHUI_REQUIRE(ledger.work < 0.0);
}

// EXT-CR05-T01: grey proton losses refuse without a table; with a table the
// split between gas heating and escape closes exactly.
void testLosses() {
  std::vector<double> e(4, 1.0);
  std::vector<double> n_e(4, 1.0);
  std::vector<double> heat(4, 0.0);
  cr::CrGreyLossLedger ledger;
  bool threw = false;
  try {
    cr::applyGreyProtonLosses(e, n_e, 1.0, cr::CrLossProvider{}, heat, 1.0, ledger);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw && e[0] == 1.0);
  const auto table = cr::parseCrLossTable(
      "# schema: chui_cr_loss_rates_v1\n# reference: synthetic test values (not physical)\n"
      "proton_coulomb_cm3_s 2.0e-3\nproton_hadronic_cm3_s 1.0e-3\nhadronic_thermalized_fraction 0.25\n",
      "synthetic");
  cr::applyGreyProtonLosses(e, n_e, 100.0, cr::CrLossProvider(table), heat, 1.0, ledger);
  CHUI_REQUIRE(closeRel(e[0], std::exp(-3.0e-3 * 100.0), 1.0e-14));
  CHUI_REQUIRE(closeRel(ledger.cr_energy_lost, ledger.gas_heating + ledger.escaping, 1.0e-14));
  CHUI_REQUIRE(closeRel(ledger.escaping / ledger.cr_energy_lost, 0.75 / 3.0, 1.0e-12));
  CHUI_REQUIRE(!cr::CrLossProvider(table).secondaryElectronSource().known());
}

}  // namespace

int main() {
  testAdvection();
  testFieldAlignedDiffusion();
  testMonotonicity();
  testStreaming();
  testTwoMoment();
  testReconstruction();
  testAdiabaticMomentumSpace();
  testSynchrotronCooling();
  testLosses();
  return 0;
}
