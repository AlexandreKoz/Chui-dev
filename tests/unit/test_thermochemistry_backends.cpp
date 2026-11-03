// CHUI Part 2 — thermochemistry networks and backends (NOT RUN — USER-RESERVED).
//
// Declared invariants and reference data:
//  * every manifest reaction balances H/He/C/O nuclei and charge (integer
//    stoichiometry); a deliberately unbalanced reaction is detected;
//  * case-A H recombination at 1e4 K agrees with Osterbrock & Ferland (2006)
//    alpha_A = 4.18e-13 cm^3 s^-1 to the known ~25% accuracy of the KWH96 fit;
//  * Compton coefficient = 4 sigma_T a k_B T_cmb^4/(m_e c) = 1.017e-37 T^4;
//  * analytic reference backend reproduces u(t) = u_eq + (u0 - u_eq) e^{-t/tau}
//    and the error decreases when the tolerance is tightened;
//  * equilibrium backend: fully ionized at 1e7 K (mu -> ~0.59), neutral at
//    3e3 K without radiation (mu -> ~1.22);
//  * hybrid network: element totals and charge conserved to round-off over
//    a strongly ionizing/cooling trajectory; species stay non-negative;
//  * cold dusty neutral gas forms H2 (dust catalysis channel);
//  * detailed network: C nuclei conserved while CO forms in shielded H2 gas;
//  * F2 for networks: too small a substep budget returns
//    kRetryShorterInterval with completed < requested, never success;
//  * UVB table: log interpolation, zero above z_max; Rahmati factor limits;
//  * metal table: provenance mandatory, H/He columns refused, asinh-space
//    interpolation is sign safe, out-of-range temperature refused.
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/physics/chemical_network.hpp"
#include "cosmosim/physics/metal_cooling_table.hpp"
#include "cosmosim/physics/radiation_background.hpp"
#include "cosmosim/physics/thermochemistry.hpp"

namespace {

using namespace cosmosim::physics;
using S = ChemicalSpecies;
constexpr std::size_t idx(S s) { return static_cast<std::size_t>(s); }

ThermochemistryCellInput primordialInput(double n_h, double u_erg_g) {
  ThermochemistryCellInput input;
  const double x_h = 1.0 - 0.2449;
  input.gas_element_mass_fraction[0] = x_h;
  input.gas_element_mass_fraction[1] = 0.2449;
  input.mass_density_cgs = n_h * chemicalElementInfo(ChemicalElement::kH).atomic_mass_amu * k_amu_g / x_h;
  input.specific_internal_energy_erg_g = u_erg_g;
  const double n_he = 0.2449 * input.mass_density_cgs / (chemicalElementInfo(ChemicalElement::kHe).atomic_mass_amu * k_amu_g);
  input.species_cm3[idx(S::kHI)] = n_h;
  input.species_cm3[idx(S::kHeI)] = n_he;
  return input;
}

double uFromTemperature(double t, double mu) { return 1.5 * k_boltzmann_erg_per_k * t / (mu * k_amu_g); }

void testManifestBalance() {
  assert(validateNetworkManifest(hybridLightNetworkManifest()).ok);
  assert(validateNetworkManifest(detailedCarbonOxygenNetworkManifest()).ok);
  NetworkManifest broken = hybridLightNetworkManifest();
  broken.reactions[0].net_change[k_slot_electron] = 0;  // H + e -> H+ without electron release
  assert(!validateNetworkManifest(broken).ok);
  for (const auto& r : detailedCarbonOxygenNetworkManifest().reactions) assert(!r.provenance.empty());
}

void testRateReferences() {
  const double alpha = rates::recombinationHII(1.0e4);
  assert(std::abs(alpha - 4.18e-13) / 4.18e-13 < 0.25);
  const double expected = 1.0178e-37 * std::pow(2.7255, 4.0);
  assert(std::abs(cooling::comptonCoolingCoefficient(0.0) - expected) / expected < 2.0e-3);
  // Collisional ionization vanishes in cold gas, recombination dominates.
  assert(rates::collisionalIonizationHI(3.0e3) < 1.0e-30);
  assert(rates::h2SelfShielding(0.0, 100.0) > 0.99);
  assert(rates::h2SelfShielding(1.0e21, 100.0) < 1.0e-4);
}

void testAnalyticReferenceExactSolution() {
  for (const double rtol : {1.0e-3, 1.0e-5}) {
    ThermochemistryConfig config;
    config.backend = ThermochemistryBackend::kAnalyticReference;
    config.relative_tolerance = rtol;
    config.analytic_relaxation_time_s = 1.0e13;
    config.analytic_equilibrium_temperature_k = 1.0e4;
    const ThermochemistrySolver solver(config);
    ThermochemistryCellInput input = primordialInput(1.0, 1.0e14);
    const double u_eq = k_boltzmann_erg_per_k * 1.0e4 / ((2.0 / 3.0) * 0.59 * k_amu_g);
    const double dt = 2.5e13;
    ThermochemistryCellOutput out;
    SourceLedger ledger;
    solver.advanceCell(input, dt, out, ledger);
    assert(out.status == SourceAdvanceStatus::kCompleted);
    assert(out.completed_seconds == dt);
    const double exact = u_eq + (1.0e14 - u_eq) * std::exp(-dt / 1.0e13);
    assert(std::abs(out.specific_internal_energy_erg_g - exact) / exact < 30.0 * rtol);
  }
}

void testEquilibriumLimits() {
  ThermochemistryConfig config;
  config.backend = ThermochemistryBackend::kEquilibriumPrimordial;
  const ThermochemistrySolver solver(config);
  const RateEnvironment env{};  // no radiation
  const NetworkState hot = solver.equilibriumPrimordialState(1.0e7, 1.0, 0.08, env);
  assert(hot[idx(S::kHII)] > 0.999);
  assert(hot[idx(S::kHeIII)] > 0.99 * 0.08);
  const NetworkState cold = solver.equilibriumPrimordialState(3.0e3, 1.0, 0.08, env);
  assert(cold[idx(S::kHI)] > 0.999);
  ThermochemistryCellInput input = primordialInput(1.0, 1.0);
  double ne_hot = hot[idx(S::kHII)] + hot[idx(S::kHeII)] + 2.0 * hot[idx(S::kHeIII)];
  const double mu_hot = ThermochemistrySolver::meanMolecularWeight(
      input.mass_density_cgs, solver.particleDensity(input, hot, ne_hot));
  assert(std::abs(mu_hot - 0.59) < 0.02);
  const double mu_cold = ThermochemistrySolver::meanMolecularWeight(
      input.mass_density_cgs, solver.particleDensity(input, cold, 0.0));
  assert(std::abs(mu_cold - 1.22) < 0.03);
}

void testHybridConservationOnIonizingTrajectory() {
  ThermochemistryConfig config;
  config.backend = ThermochemistryBackend::kFixedHybrid;
  const ThermochemistrySolver solver(config);
  // Neutral gas suddenly at ~2e5 K: collisional ionization + cooling.
  ThermochemistryCellInput input = primordialInput(0.1, uFromTemperature(2.0e5, 1.22));
  const double n_h = 0.1;
  const double n_he = input.species_cm3[idx(S::kHeI)];
  ThermochemistryCellOutput out;
  SourceLedger ledger;
  solver.advanceCell(input, 3.15576e13, out, ledger);
  assert(out.status == SourceAdvanceStatus::kCompleted);
  const auto& n = out.species_cm3;
  for (std::size_t s = 0; s < k_hybrid_species_count; ++s) assert(n[s] >= 0.0);
  const double h = n[idx(S::kHI)] + n[idx(S::kHII)] + n[idx(S::kHm)] + 2.0 * n[idx(S::kH2)] + 2.0 * n[idx(S::kH2p)];
  assert(std::abs(h - n_h) / n_h < 1.0e-12);
  const double he = n[idx(S::kHeI)] + n[idx(S::kHeII)] + n[idx(S::kHeIII)];
  assert(std::abs(he - n_he) / n_he < 1.0e-12);
  assert(ledger.charge_residual < 1.0e-12);
  assert(n[idx(S::kHII)] > 0.5 * n_h);         // collisionally ionized
  assert(ledger.radiative_loss_erg > 0.0);
  assert(ledger.chemical_binding_delta_erg > 0.0);  // ionization stored energy
  assert(out.mean_molecular_weight < 1.0);
}

void testDustCatalysedH2Formation() {
  ThermochemistryConfig config;
  config.backend = ThermochemistryBackend::kFixedHybrid;
  const ThermochemistrySolver solver(config);
  ThermochemistryCellInput input = primordialInput(100.0, uFromTemperature(100.0, 1.22));
  input.dust_to_gas_relative_mw = 1.0;
  ThermochemistryCellOutput out;
  SourceLedger ledger;
  solver.advanceCell(input, 3.15576e14, out, ledger);  // 10 Myr
  assert(out.status == SourceAdvanceStatus::kCompleted);
  // t_form ~ 1/(3e-17 * 100) ~ 10 Myr: an O(1) molecular fraction.
  assert(2.0 * out.species_cm3[idx(S::kH2)] / 100.0 > 0.1);
  // No dust and no electrons => negligible H2.
  input.dust_to_gas_relative_mw = 0.0;
  ThermochemistryCellOutput dry;
  solver.advanceCell(input, 3.15576e14, dry, ledger);
  assert(dry.species_cm3[idx(S::kH2)] < 1.0e-3 * out.species_cm3[idx(S::kH2)]);
}

void testDetailedCarbonConservation() {
  ThermochemistryConfig config;
  config.backend = ThermochemistryBackend::kFixedDetailed;
  config.isrf_g0_habing = 1.0;
  config.metal_mode = MetalCoolingMode::kFineStructureCold;
  const ThermochemistrySolver solver(config);
  ThermochemistryCellInput input = primordialInput(1.0e4, uFromTemperature(30.0, 2.3));
  input.gas_element_mass_fraction[2] = 2.4e-3;  // C
  input.gas_element_mass_fraction[4] = 5.8e-3;  // O
  input.dust_to_gas_relative_mw = 1.0;
  const double n_c = 2.4e-3 * input.mass_density_cgs / (chemicalElementInfo(ChemicalElement::kC).atomic_mass_amu * k_amu_g);
  input.species_cm3[idx(S::kHI)] = 0.0;
  input.species_cm3[idx(S::kH2)] = 0.5e4;
  input.species_cm3[idx(S::kCII)] = n_c;
  ThermochemistryCellOutput out;
  SourceLedger ledger;
  solver.advanceCell(input, 3.15576e13, out, ledger);
  assert(out.status == SourceAdvanceStatus::kCompleted);
  const double c_total = out.species_cm3[idx(S::kCII)] + out.species_cm3[idx(S::kCO)];
  assert(std::abs(c_total - n_c) / n_c < 1.0e-12);
  assert(out.species_cm3[idx(S::kCO)] > 0.0);
}

void testNetworkF2Rejection() {
  ThermochemistryConfig config;
  config.backend = ThermochemistryBackend::kFixedHybrid;
  config.max_substeps = 1U;
  const ThermochemistrySolver solver(config);
  ThermochemistryCellInput input = primordialInput(10.0, uFromTemperature(1.0e6, 0.6));
  ThermochemistryCellOutput out;
  SourceLedger ledger;
  solver.advanceCell(input, 3.15576e15, out, ledger);
  assert(out.status == SourceAdvanceStatus::kRetryShorterInterval);
  assert(out.completed_seconds < 3.15576e15);
  assert(out.suggested_retry_seconds > 0.0);
}

void testUvBackgroundTableAndShielding() {
  std::vector<double> lz{0.0, std::log10(3.0), std::log10(7.0)};
  std::vector<UvBackgroundRates> rows{{1e-13, 1e-13, 1e-15, 1e-24, 1e-24, 1e-25},
                                      {1e-12, 1e-12, 1e-14, 1e-23, 1e-23, 1e-24},
                                      {1e-14, 1e-14, 0.0, 1e-25, 1e-25, 0.0}};
  const auto table = UvBackgroundTable::fromRows(lz, rows);
  const auto mid = table.ratesAtRedshift(std::sqrt(3.0) - 1.0);  // halfway in log(1+z)
  assert(std::abs(mid.gamma_hi_s - std::sqrt(1e-13 * 1e-12)) / mid.gamma_hi_s < 1.0e-12);
  assert(table.ratesAtRedshift(10.0).gamma_hi_s == 0.0);  // above z_max: no background
  assert(rahmati2013ShieldingFactor(1.0e-6, 1.0e4, 1.0e-12) > 0.999);
  const double dense = rahmati2013ShieldingFactor(1.0e3, 1.0e4, 1.0e-12);
  assert(dense < 1.0e-2 && dense > 0.0);
  bool threw = false;
  try { (void)UvBackgroundTable::loadTreecool(""); } catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
}

std::string tinyMetalTable(bool with_h) {
  std::string t = "# chui_metal_cooling_table_v1\n# citation = unit-test fixture (not astrophysical data)\n"
                  "# license = CC0 test fixture\n# solar_reference = Asplund et al. 2009\n";
  t += with_h ? "elements H O\n" : "elements C O\n";
  t += "log10_T 4 5\nlog10_nH -2 0\nredshift 0\nreference_mass_fraction 2.4e-3 5.8e-3\n";
  t += with_h ? "lambda H 1e-22 1e-22 1e-22 1e-22\n" : "lambda C 1e-23 1e-23 -1e-24 -1e-24\n";
  t += "lambda O 1e-22 1e-22 1e-21 1e-21\n";
  return t;
}

void testMetalTable() {
  const auto table = MetalCoolingTable::loadFromString(tinyMetalTable(false), "fixture");
  MetalCoolingTable::Query q;
  q.temperature_k = std::pow(10.0, 4.5);
  q.hydrogen_number_density_cm3 = 0.1;
  q.gas_mass_fraction[2] = 2.4e-3;
  q.gas_mass_fraction[4] = 5.8e-3;
  const auto r = table.evaluate(q, MetalTableRangePolicy::kReject);
  assert(std::isfinite(r.net_cooling_erg_cm3_s));
  // Depletion halves O and C => cooling halves (linear abundance scaling).
  q.gas_mass_fraction[2] *= 0.5;
  q.gas_mass_fraction[4] *= 0.5;
  const auto half = table.evaluate(q, MetalTableRangePolicy::kReject);
  assert(std::abs(half.net_cooling_erg_cm3_s - 0.5 * r.net_cooling_erg_cm3_s) <=
         1.0e-12 * std::abs(r.net_cooling_erg_cm3_s));
  q.temperature_k = 1.0e6;
  bool threw = false;
  try { (void)table.evaluate(q, MetalTableRangePolicy::kClampDensityAndRedshift); }
  catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
  threw = false;
  try { (void)MetalCoolingTable::loadFromString(tinyMetalTable(true), "fixture"); }
  catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
  std::string no_citation = tinyMetalTable(false);
  no_citation.erase(no_citation.find("# citation"), no_citation.find("# license") - no_citation.find("# citation"));
  threw = false;
  try { (void)MetalCoolingTable::loadFromString(no_citation, "fixture"); }
  catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
}

void testMissingTableRefusesModel() {
  ThermochemistryConfig config;
  config.backend = ThermochemistryBackend::kFixedHybrid;
  config.metal_mode = MetalCoolingMode::kElementTable;
  bool threw = false;
  try { const ThermochemistrySolver solver(config); } catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
}

}  // namespace

int main() {
  testManifestBalance();
  testRateReferences();
  testAnalyticReferenceExactSolution();
  testEquilibriumLimits();
  testHybridConservationOnIonizingTrajectory();
  testDustCatalysedH2Formation();
  testDetailedCarbonConservation();
  testNetworkF2Rejection();
  testUvBackgroundTableAndShielding();
  testMetalTable();
  testMissingTableRefusesModel();
  return 0;
}
