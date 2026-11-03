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
//    interpolation is sign safe, out-of-range temperature refused;
//  * Part 2.5: WSS09 eq. 4 abundance scaling relative to HYDROGEN
//    ((X_e/X_H)/(X_e,ref/X_H,ref)); uv_background and
//    reference_hydrogen_mass_fraction mandatory; UVB mismatch refused unless
//    explicitly accepted;
//  * Part 2.5 interpolation safety: a SINGLETON axis (T, n_H or z) covers
//    only its own coordinate (tolerance 1e-8 x max(1,|value|)). A
//    one-redshift (z = 0) table evaluated at z = 3 is REFUSED under kReject
//    and clamped + reported under kClampDensityAndRedshift; temperature is
//    never clamped.
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
  config.allow_unqualified_rates = true;
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
  config.allow_unqualified_rates = true;
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
  config.allow_unqualified_rates = true;
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
  config.allow_unqualified_rates = true;
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
  config.allow_unqualified_rates = true;
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
  config.allow_unqualified_rates = true;
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
                  "# license = CC0 test fixture\n# solar_reference = Asplund et al. 2009\n"
                  "# uv_background = none\n";
  t += with_h ? "elements H O\n" : "elements C O\n";
  t += "log10_T 4 5\nlog10_nH -2 0\nredshift 0\nreference_mass_fraction 2.4e-3 5.8e-3\n"
       "reference_hydrogen_mass_fraction 0.7065\n";
  t += with_h ? "lambda H 1e-22 1e-22 1e-22 1e-22\n" : "lambda C 1e-23 1e-23 -1e-24 -1e-24\n";
  t += "lambda O 1e-22 1e-22 1e-21 1e-21\n";
  return t;
}

void testMetalTable() {
  const auto table = MetalCoolingTable::loadFromString(tinyMetalTable(false), "fixture");
  MetalCoolingTable::Query q;
  q.temperature_k = std::pow(10.0, 4.5);
  q.hydrogen_number_density_cm3 = 0.1;
  q.gas_mass_fraction[0] = 0.7065;
  q.gas_mass_fraction[2] = 2.4e-3;
  q.gas_mass_fraction[4] = 5.8e-3;
  const auto r = table.evaluate(q, MetalTableRangePolicy::kReject);
  assert(std::isfinite(r.net_cooling_erg_cm3_s));
  {
    // At exactly the reference abundances the result is the table value
    // times n_H^2. At log T = 4.5, mid n_H column interpolation in asinh
    // space: evaluate the reference from the same table points.
    MetalCoolingTable::Query q_ref = q;
    q_ref.temperature_k = 1.0e4;
    q_ref.hydrogen_number_density_cm3 = 1.0e-2;  // grid corner: C = 1e-23, O = 1e-22
    const auto corner = table.evaluate(q_ref, MetalTableRangePolicy::kReject);
    assert(std::abs(corner.net_cooling_erg_cm3_s - (1.0e-23 + 1.0e-22) * 1.0e-4) <= 1.0e-12 * 1.1e-26);
    // WSS09 eq. 4: halving X_H at fixed X_e and n_H doubles [e/H] => doubles
    // the metal term. (The pre-Part-2.5 scaling X_e/X_e,ref left it fixed.)
    q_ref.gas_mass_fraction[0] = 0.5 * 0.7065;
    const auto doubled = table.evaluate(q_ref, MetalTableRangePolicy::kReject);
    assert(std::abs(doubled.net_cooling_erg_cm3_s - 2.0 * corner.net_cooling_erg_cm3_s) <=
           1.0e-12 * std::abs(corner.net_cooling_erg_cm3_s));
  }
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
  for (const std::string& drop : {std::string("# uv_background = none\n"),
                                  std::string("reference_hydrogen_mass_fraction 0.7065\n")}) {
    std::string missing = tinyMetalTable(false);
    missing.erase(missing.find(drop), drop.size());
    threw = false;
    try { (void)MetalCoolingTable::loadFromString(missing, "fixture"); }
    catch (const SourceModelUnavailableError&) { threw = true; }
    assert(threw);
  }
  std::string no_citation = tinyMetalTable(false);
  no_citation.erase(no_citation.find("# citation"), no_citation.find("# license") - no_citation.find("# citation"));
  threw = false;
  try { (void)MetalCoolingTable::loadFromString(no_citation, "fixture"); }
  catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
}

// Fixture with configurable axes; every lambda value is 1e-22 (C only).
std::string axisMetalTable(const std::string& log_t, std::size_t nt, const std::string& log_n, std::size_t nn,
                           const std::string& z, std::size_t nz) {
  std::string t = "# chui_metal_cooling_table_v1\n# citation = unit-test fixture (not astrophysical data)\n"
                  "# license = CC0 test fixture\n# solar_reference = fixture\n# uv_background = none\n"
                  "elements C\nlog10_T " + log_t + "\nlog10_nH " + log_n + "\nredshift " + z +
                  "\nreference_mass_fraction 2.4e-3\nreference_hydrogen_mass_fraction 0.7065\nlambda C";
  for (std::size_t i = 0; i < nt * nn * nz; ++i) t += " 1e-22";
  return t + "\n";
}

bool metalRejects(const MetalCoolingTable& table, const MetalCoolingTable::Query& q, MetalTableRangePolicy policy) {
  try { (void)table.evaluate(q, policy); } catch (const SourceModelUnavailableError&) { return true; }
  return false;
}

void testMetalTableSingletonAxes() {
  MetalCoolingTable::Query q;
  q.temperature_k = std::pow(10.0, 4.5);
  q.hydrogen_number_density_cm3 = 0.1;
  q.gas_mass_fraction[0] = 0.7065;
  q.gas_mass_fraction[2] = 2.4e-3;
  const double expected = 1.0e-22 * 0.1 * 0.1;
  {
    // One-redshift table (as converted from WSS09 z_0.000.hdf5 alone).
    const auto table = MetalCoolingTable::loadFromString(axisMetalTable("4 5", 2, "-2 0", 2, "0", 1), "z0");
    q.redshift = 0.0;
    auto r = table.evaluate(q, MetalTableRangePolicy::kReject);
    assert(!r.clamped && std::abs(r.net_cooling_erg_cm3_s - expected) <= 1.0e-12 * expected);
    q.redshift = 1.0e-12;  // round-off of 1/a - 1 at a = 1: inside tolerance, not a clamp
    r = table.evaluate(q, MetalTableRangePolicy::kReject);
    assert(!r.clamped);
    q.redshift = 3.0;  // the regression: z = 0 rates must not be used at z = 3
    assert(metalRejects(table, q, MetalTableRangePolicy::kReject));
    r = table.evaluate(q, MetalTableRangePolicy::kClampDensityAndRedshift);
    assert(r.clamped);  // explicit clamping is allowed but always reported
    assert(std::abs(r.net_cooling_erg_cm3_s - expected) <= 1.0e-12 * expected);  // values unchanged
    q.redshift = 1.0e-6;  // beyond tolerance: a genuine (small) offset is still out of range
    assert(metalRejects(table, q, MetalTableRangePolicy::kReject));
    q.redshift = 0.0;
  }
  {
    // One-density table.
    const auto table = MetalCoolingTable::loadFromString(axisMetalTable("4 5", 2, "-1", 1, "0 1", 2), "n1");
    assert(!table.evaluate(q, MetalTableRangePolicy::kReject).clamped);  // n_H = 0.1 matches
    q.hydrogen_number_density_cm3 = 1.0;
    assert(metalRejects(table, q, MetalTableRangePolicy::kReject));
    assert(table.evaluate(q, MetalTableRangePolicy::kClampDensityAndRedshift).clamped);
    q.hydrogen_number_density_cm3 = 0.1;
  }
  {
    // One-temperature table: temperature is NEVER clamped, under any policy.
    const auto table = MetalCoolingTable::loadFromString(axisMetalTable("4", 1, "-2 0", 2, "0 1", 2), "t1");
    q.temperature_k = 1.0e4;
    assert(!table.evaluate(q, MetalTableRangePolicy::kReject).clamped);
    q.temperature_k = 2.0e4;
    assert(metalRejects(table, q, MetalTableRangePolicy::kReject));
    assert(metalRejects(table, q, MetalTableRangePolicy::kClampDensityAndRedshift));
  }
  {
    // Multi-point axes: an edge query within round-off is on the grid; a
    // real overshoot is refused (reject) or clamped and reported.
    const auto table = MetalCoolingTable::loadFromString(axisMetalTable("4 5", 2, "-2 0", 2, "0 1", 2), "multi");
    q.temperature_k = 1.0e4;
    q.hydrogen_number_density_cm3 = 1.0 * (1.0 + 1.0e-14);
    q.redshift = 1.0 + 1.0e-12;
    assert(!table.evaluate(q, MetalTableRangePolicy::kReject).clamped);
    q.redshift = 1.5;
    assert(metalRejects(table, q, MetalTableRangePolicy::kReject));
    assert(table.evaluate(q, MetalTableRangePolicy::kClampDensityAndRedshift).clamped);
    q.temperature_k = 1.0e5 * 1.001;
    assert(metalRejects(table, q, MetalTableRangePolicy::kClampDensityAndRedshift));
  }
}

void testMetalTableUvBackgroundPolicy() {
  // Collisional table ("none") is compatible only with no CHUI UVB.
  const auto collisional = MetalCoolingTable::loadFromString(tinyMetalTable(false), "fixture");
  assert(collisional.checkUvBackgroundCompatibility("", false).compatible);
  bool threw = false;
  try { (void)collisional.checkUvBackgroundCompatibility("abc123", false); }
  catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
  // Photoionized table (WSS09-like, HM01) against a FG20 TREECOOL: refused
  // by default, accepted only explicitly (and flagged).
  std::string hm01 = tinyMetalTable(false);
  hm01.replace(hm01.find("# uv_background = none"), std::string("# uv_background = none").size(),
               "# uv_background = HM01+CMB (Haardt & Madau 2001), CLOUDY 07.02");
  const auto photo = MetalCoolingTable::loadFromString(hm01, "fixture");
  threw = false;
  try { (void)photo.checkUvBackgroundCompatibility("fg20sha", false); }
  catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
  threw = false;
  try { (void)photo.checkUvBackgroundCompatibility("", false); }  // photoionized table, CHUI without UVB
  catch (const SourceModelUnavailableError&) { threw = true; }
  assert(threw);
  const auto accepted = photo.checkUvBackgroundCompatibility("fg20sha", true);
  assert(!accepted.compatible && accepted.accepted_mismatch && !accepted.diagnostic.empty());
  // Exact TREECOOL identity declared by the table => compatible.
  std::string matched = hm01;
  matched.insert(matched.find("elements"), "# uv_background_treecool_sha256 = fg20sha\n");
  const auto same = MetalCoolingTable::loadFromString(matched, "fixture");
  assert(same.checkUvBackgroundCompatibility("fg20sha", false).compatible);
}

void testUnqualifiedRatesRefusedByDefault() {
  // Corrective gate: unverified coefficients must not run silently.
  for (const auto backend : {ThermochemistryBackend::kEquilibriumPrimordial, ThermochemistryBackend::kFixedHybrid,
                             ThermochemistryBackend::kFixedDetailed}) {
    ThermochemistryConfig config;
    config.backend = backend;
    bool threw = false;
    try { const ThermochemistrySolver solver(config); } catch (const SourceModelUnavailableError&) { threw = true; }
    assert(threw);
  }
  assert(!hybridLightNetworkManifest().qualified());
  ThermochemistryConfig analytic;
  analytic.backend = ThermochemistryBackend::kAnalyticReference;
  const ThermochemistrySolver ok(analytic);  // reference model needs no rate data
  (void)ok;
}

void testMissingTableRefusesModel() {
  ThermochemistryConfig config;
  config.allow_unqualified_rates = true;
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
  testMetalTableUvBackgroundPolicy();
  testMetalTableSingletonAxes();
  testMissingTableRefusesModel();
  testUnqualifiedRatesRefusedByDefault();
  return 0;
}
