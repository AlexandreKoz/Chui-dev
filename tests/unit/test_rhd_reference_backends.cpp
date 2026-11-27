// CHUI-EXT-PHYS-REF-v1 gate VR (R03-R06, F01, F02): FLD limiter limits,
// thick-limit Gaussian diffusion against the analytic spreading law, energy
// conservation and LTE relaxation; discrete-ordinates exact uniform solution
// and isotropic Eddington tensor; VET energy budget closure and thick-limit
// closure; causal-ray light-travel, exact attenuation and exactly-once ledger,
// queue refusal and restart round trip; HG Monte Carlo asymmetry moment, slab
// transmission against exp(-tau), energy closure; stellar luminosity table
// parsing, exact age integral and domain refusals; F01 workspace formula.
// Expected values are analytic. NOT RUN in pass 01.
#include "../support/test_require.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/lightcone/ray_geometry.hpp"
#include "cosmosim/radiation/radiation_causal_transport.hpp"
#include "cosmosim/radiation/radiation_reference_transport.hpp"
#include "cosmosim/radiation/radiation_runtime.hpp"
#include "cosmosim/radiation/radiation_stellar_sources.hpp"

namespace r = cosmosim::radiation;
namespace lc = cosmosim::lightcone;

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

bool closeRel(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

// EXT-R03-T01: Levermore-Pomraning limits and thick-limit Gaussian diffusion.
void testFld() {
  CHUI_REQUIRE(closeRel(r::levermorePomraningLimiter(0.0), 1.0 / 3.0, 1.0e-15));
  CHUI_REQUIRE(closeRel(r::levermorePomraningLimiter(1.0e8) * 1.0e8, 1.0, 1.0e-6));
  r::UniformPatch patch;
  patch.n = {128, 1, 1};
  patch.dx = {1.0 / 128.0, 1.0, 1.0};
  patch.boundary = r::ReferenceBoundary::kPeriodic;
  const std::size_t n = patch.cells();
  const double chi = 1.0e3;  // optically thick: R << 1
  std::vector<double> e(n);
  const double sigma0 = 0.05;
  for (std::size_t i = 0; i < n; ++i) {
    const double x = (static_cast<double>(i) + 0.5) / 128.0 - 0.5;
    e[i] = 1.0 + std::exp(-x * x / (2.0 * sigma0 * sigma0));  // background keeps R small
  }
  const std::vector<double> zero(n, 0.0);
  const std::vector<double> chi_s(n, chi);
  r::FldSettings settings;
  double total0 = 0.0;
  for (const double v : e) total0 += v;
  const double d = 1.0 / (3.0 * chi);  // c = 1
  const double t_end = 0.5 * sigma0 * sigma0 / d;  // sigma^2 doubles: sigma^2 = sigma0^2 + 2 D t
  const int steps = 200;
  for (int s = 0; s < steps; ++s) {
    (void)r::advanceFldImplicit(patch, e, zero, chi_s, zero, zero, t_end / steps, settings);
  }
  double total1 = 0.0;
  for (const double v : e) total1 += v;
  CHUI_REQUIRE(closeRel(total0, total1, 1.0e-10));  // conservative fluxes, periodic, no absorption
  // Peak amplitude of a 1-D Gaussian scales as sigma0 / sigma = 1/sqrt(2).
  const double peak = e[n / 2U] - 1.0;
  CHUI_REQUIRE(std::abs(peak - 1.0 / std::sqrt(2.0)) < 0.02);
  // LTE relaxation: absorbing medium relaxes to a T^4.
  std::vector<double> e2(n, 0.0);
  const std::vector<double> chi_a(n, 10.0);
  const std::vector<double> temperature(n, 1.0);
  for (int s = 0; s < 50; ++s) (void)r::advanceFldImplicit(patch, e2, chi_a, zero, temperature, zero, 1.0, settings);
  CHUI_REQUIRE(closeRel(e2[7], 1.0, 1.0e-10));
}

// EXT-R03-T02: discrete ordinates in a uniform periodic emitting medium is exact.
void testOrdinates() {
  r::UniformPatch patch;
  patch.n = {8, 8, 8};
  patch.dx = {0.125, 0.125, 0.125};
  const std::size_t n = patch.cells();
  const std::vector<double> eta(n, 2.0);
  const std::vector<double> chi_a(n, 4.0);
  const std::vector<double> chi_s(n, 0.0);
  r::OrdinatesSettings settings;
  settings.directions = 128;
  const auto sol = r::solveDiscreteOrdinates(patch, eta, chi_a, chi_s, settings);
  CHUI_REQUIRE(sol.converged);
  for (const double i : sol.intensity) CHUI_REQUIRE(closeRel(i, 0.5, 1.0e-8));  // I = eta / chi
  const auto f = r::eddingtonTensors(patch, sol);
  CHUI_REQUIRE(closeRel(f[0][0] + f[0][4] + f[0][8], 1.0, 1.0e-12));
  CHUI_REQUIRE(std::abs(f[0][0] - 1.0 / 3.0) < 0.02 && std::abs(f[0][1]) < 0.02);
  // Pure scatterer with absorption: J converges to eta / chi_a in the uniform periodic case.
  const std::vector<double> chi_s2(n, 3.0);
  const std::vector<double> chi_a2(n, 1.0);
  r::OrdinatesSettings slow = settings;
  slow.max_sweeps = 4000;  // source iteration contracts by ~(chi_s + w)/(chi_t + w) per sweep
  const auto scat = r::solveDiscreteOrdinates(patch, eta, chi_a2, chi_s2, slow);
  CHUI_REQUIRE(scat.converged && closeRel(scat.mean_intensity[5], 2.0, 1.0e-7));
}

// EXT-R03-T03: VET energy budget and thick-limit closure.
void testVet() {
  r::UniformPatch patch;
  patch.n = {8, 8, 8};
  patch.dx = {0.125, 0.125, 0.125};
  patch.boundary = r::ReferenceBoundary::kPeriodic;
  const std::size_t n = patch.cells();
  r::VetSettings settings;
  settings.ordinates.directions = 96;
  r::VetTransport vet(patch, settings);
  r::VetState state;
  state.energy.assign(n, 0.0);
  state.flux.assign(n, {0.0, 0.0, 0.0});
  const std::vector<double> chi_a(n, 2.0);
  const std::vector<double> chi_s(n, 0.0);
  const std::vector<double> temperature(n, 1.0);
  double e0 = 0.0;
  const double volume = 0.125 * 0.125 * 0.125;
  const auto report = vet.advance(state, chi_a, chi_s, temperature, 0.3);
  double e1 = 0.0;
  for (const double v : state.energy) e1 += v * volume;
  CHUI_REQUIRE(closeRel(e1 - e0, report.emitted_energy - report.absorbed_energy - report.boundary_outflow_energy, 1.0e-10));
  CHUI_REQUIRE(report.closure_updates >= 1U && report.substeps >= 1U);
  CHUI_REQUIRE(std::abs(state.eddington[3][0] - 1.0 / 3.0) < 0.02);
  // Relaxation toward a T^4 = 1 with rate c chi_a: E(t) = 1 - exp(-c chi_a t) (uniform medium).
  CHUI_REQUIRE(std::abs(state.energy[17] - (1.0 - std::exp(-2.0 * 0.3))) < 0.02);
}

// EXT-R04-T01: causal rays: light travel distance, exact attenuation, ledger closure.
void testCausalRays() {
  const auto geometry = lc::LeafGeometry::uniformGrid({16, 16, 16}, {0, 0, 0}, {1, 1, 1}, lc::RayBoundary::kOpen);
  r::CausalRayTransport transport(geometry, 1.0, 1024, 1.0e-14, 4096);
  transport.emitIsotropic({0.5, 0.5, 0.5}, 1.0, 0.0, 64, 1U);
  const std::vector<double> chi(geometry.cellCount(), 2.0);
  std::vector<double> absorbed(geometry.cellCount(), 0.0);
  transport.advance(0.2, chi, absorbed);
  for (const auto& p : transport.packets()) {
    CHUI_REQUIRE(closeRel(p.travelled, 0.2, 1.0e-12));                       // c t
    CHUI_REQUIRE(closeRel(p.energy, p.emitted_energy * std::exp(-2.0 * 0.2), 1.0e-12));  // exp(-chi c t)
  }
  CHUI_REQUIRE(std::abs(transport.ledger().closureResidual(transport.inFlightEnergy())) < 1.0e-14);
  // Restart round trip.
  const auto bytes = transport.serialize();
  r::CausalRayTransport restored(geometry, 1.0, 1024, 1.0e-14, 4096);
  restored.restore(bytes);
  CHUI_REQUIRE(restored.packets().size() == transport.packets().size() && restored.currentTime() == 0.2);
  // All packets leave the unit box by t = 1 (max distance to a face is sqrt(3)/2 < 1).
  transport.advance(1.0, chi, absorbed);
  CHUI_REQUIRE(transport.packets().empty());
  CHUI_REQUIRE(std::abs(transport.ledger().closureResidual(0.0)) < 1.0e-14);
  double deposited = 0.0;
  for (const double v : absorbed) deposited += v;
  CHUI_REQUIRE(closeRel(deposited, transport.ledger().absorbed, 1.0e-12));
  bool threw = false;
  try {
    transport.emitIsotropic({0.5, 0.5, 0.5}, 1.0, 1.0, 2048, 2U);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

// EXT-R05-T01: HG sampling and slab transmission.
void testHenyeyGreenstein() {
  for (const double g : {0.0, 0.5, -0.3}) {
    double sum = 0.0;
    const int n = 200000;
    for (int i = 0; i < n; ++i) sum += r::sampleHenyeyGreensteinCosine(g, (i + 0.5) / n);
    CHUI_REQUIRE(std::abs(sum / n - g) < 1.0e-3);  // <mu> = g (stratified quadrature of the inverse CDF)
  }
  // Pure absorber thin slab: escape fraction of an isotropic central source in
  // a sphere-like cube is not analytic, so use a 1-D slab along x with
  // transverse periodicity and a source on the -x face directed... instead:
  // uniform absorbing cube, compare the escaped fraction with the exact
  // angular average of exp(-chi l(direction)) from the cube centre.
  const auto geometry = lc::LeafGeometry::uniformGrid({8, 8, 8}, {0, 0, 0}, {1, 1, 1}, lc::RayBoundary::kOpen);
  const std::vector<double> chi_a(geometry.cellCount(), 1.5);
  const std::vector<double> chi_s(geometry.cellCount(), 0.0);
  r::HgMonteCarloSettings settings;
  settings.packets = 200000;
  settings.batches = 20;
  const auto mc = r::runHenyeyGreensteinMonteCarlo(geometry, chi_a, chi_s, {0.5, 0.5, 0.5}, 1.0, settings);
  CHUI_REQUIRE(std::abs(mc.emitted - mc.absorbed - mc.escaped - mc.truncated) < 1.0e-12);
  // Exact: <exp(-chi l(n))> over directions, l = distance from the centre to the cube face.
  double exact = 0.0;
  const auto dirs = r::fibonacciDirections(20000);
  for (const auto& d : dirs) {
    double l = 1e300;
    for (std::size_t a = 0; a < 3U; ++a) {
      if (d[a] != 0.0) l = std::min(l, 0.5 / std::abs(d[a]));
    }
    exact += std::exp(-1.5 * l);
  }
  exact /= static_cast<double>(dirs.size());
  CHUI_REQUIRE(std::abs(mc.escaped_fraction_mean - exact) < 5.0 * mc.escaped_fraction_standard_error + 1.0e-3);
}

std::string stellarTable() {
  std::string text =
      "# schema: chui_ssp_group_luminosity_v1\n# source: synthetic_test_table\n# imf: test\n"
      "# group_edges_ev: 13.598434, 24.587389, 54.417765, 100\n# metallicities: 0.001, 0.02\n"
      "# columns: age_yr metallicity group L_erg_s_per_msun Q_photons_s_per_msun\n";
  for (double age : {0.0, 1.0e6, 1.0e7}) {
    for (double z : {0.001, 0.02}) {
      for (int g = 0; g < 3; ++g) {
        const double l = (age < 2.0e6 ? 1.0e34 : 0.0) * (z < 0.01 ? 2.0 : 1.0) / (g + 1);
        text += std::to_string(age) + " " + std::to_string(z) + " " + std::to_string(g) + " " + std::to_string(l) +
                " " + std::to_string(l / 3.0e-11) + "\n";
      }
    }
  }
  return text;
}

// EXT-F02-T01: stellar luminosity table: parsing, exact piecewise-linear age
// integral, metallicity interpolation in log Z, domain refusal.
void testStellarTable() {
  const auto table = r::parseStellarLuminosityTable(stellarTable(), "synthetic");
  CHUI_REQUIRE(table.groupCount() == 3U && table.age_yr.size() == 3U);
  // L(age) = 1e34 on [0, 1e6], linear to 0 at 1e7: int_0^5e6 = 1e34 (1e6 + 0.5 (4e6)(1 + 5/9)) yr.
  const double year_s = 3.15576e7;
  const double l_mid = 1.0e34 * (1.0 - 4.0 / 9.0);
  const double expected = (1.0e34 * 1.0e6 + 0.5 * 4.0e6 * (1.0e34 + l_mid)) * year_s;
  CHUI_REQUIRE(closeRel(table.integratedEnergy(0, 0.0, 5.0e6, 0.02), expected, 1.0e-12));
  // log-Z midpoint of 0.001 and 0.02: L = 1.5e34 at age 0 for group 0.
  CHUI_REQUIRE(closeRel(table.luminosity(0, 0.0, std::sqrt(0.001 * 0.02)), 1.5e34, 1.0e-12));
  bool threw = false;
  try {
    (void)table.luminosity(0, 2.0e7, 0.02);  // beyond the last age
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  threw = false;
  try {
    (void)table.luminosity(0, 1.0e6, 0.03);  // beyond the last metallicity
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
  // F01 workspace formula: 25 lanes of nvalues plus 3 gas copies.
  CHUI_REQUIRE(r::radiationIntervalWorkspaceBytes(100, 10, true) == 25U * 100U * 8U + 3U * 10U * (4U + 48U));
  CHUI_REQUIRE(r::radiationIntervalWorkspaceBytes(100, 10, false) == 25U * 100U * 8U);
}

}  // namespace

int main() {
  testFld();
  testOrdinates();
  testVet();
  testCausalRays();
  testHenyeyGreenstein();
  testStellarTable();
  return 0;
}
