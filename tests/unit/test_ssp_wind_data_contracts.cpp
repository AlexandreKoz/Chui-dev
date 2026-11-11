// CHUI patch 58 (S02/F02): explicit stellar-yield coverage policy and the
// machine-readable wind-table IMF/metallicity compatibility contract.
// Authored without execution: NOT RUN — USER-RESERVED.
//
// Every table below is a SYNTHETIC unit-test fixture (not astrophysical data);
// the numbers only exercise the gates and carry no scientific meaning.

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/physics/source_interval.hpp"
#include "cosmosim/physics/stellar_population.hpp"
#include "cosmosim/physics/stellar_wind.hpp"
#include "../support/test_require.hpp"

namespace {

namespace ph = cosmosim::physics;

const std::string k_row_ok = " 2.0 -0.08 0.01 0.02 0 0.05 0 0 0 0 0 0\n";  // net yields sum to 0
const std::string k_ia_row = " 1.4 0 0 0.05 0 0.15 0 0.01 0.15 0.08 0.01 0.95\n";

// Two metallicities on every channel; AGB mass axis `agb_lo..8`, CC 8..40.
std::string coveringYieldTable(const std::string& agb_lo) {
  std::string t = "# chui_ssp_yield_table_v1\n# citation = unit-test fixture (not astrophysical data)\n"
                  "# license = CC0 test fixture\nelements H He C N O Ne Mg Si S Ca Fe\n";
  t += "channel agb\nmetallicity 0.001 0.02\nmass " + agb_lo + " 8\n";
  for (int z = 0; z < 2; ++z) {
    for (int m = 0; m < 2; ++m) t += "row " + std::to_string(z) + " " + std::to_string(m) + k_row_ok;
  }
  t += "channel cc\nmetallicity 0.001 0.02\nmass 8 40\n";
  for (int z = 0; z < 2; ++z) {
    for (int m = 0; m < 2; ++m) t += "row " + std::to_string(z) + " " + std::to_string(m) + k_row_ok;
  }
  t += "channel ia\nmetallicity 0.001 0.02\nrow 0" + k_ia_row + "row 1" + k_ia_row;
  return t;
}

// Staged-data shape: CC grid 10..30 Msun (< configured 8..40), single Z.
std::string narrowYieldTable() {
  std::string t = "# chui_ssp_yield_table_v1\n# citation = unit-test fixture (not astrophysical data)\n"
                  "# license = CC0 test fixture\nelements H He C N O Ne Mg Si S Ca Fe\n";
  t += "channel agb\nmetallicity 0.02\nmass 1 5\nrow 0 0" + k_row_ok + "row 0 1" + k_row_ok;
  t += "channel cc\nmetallicity 0.02\nmass 10 30\nrow 0 0" + k_row_ok + "row 0 1" + k_row_ok;
  t += "channel ia\nmetallicity 0.02\nrow 0" + k_ia_row;
  return t;
}

std::array<double, ph::k_chemical_element_count> birthMixture(double metals) {
  std::array<double, ph::k_chemical_element_count> x{};
  x[4] = metals;  // all metals as O (fixture)
  x[1] = 0.26;
  x[0] = 1.0 - 0.26 - metals;
  return x;
}

template <class F>
bool throwsUnavailable(F&& f) {
  try {
    f();
  } catch (const ph::SourceModelUnavailableError&) {
    return true;
  }
  return false;
}

ph::SspConfig tableConfig(ph::SspYieldCoveragePolicy policy) {
  ph::SspConfig config;  // IMF 0.1..100, CC 8..40 (defaults)
  config.yield_mode = ph::SspYieldMode::kTable;
  config.yield_coverage = policy;
  return config;
}

void testYieldStrictConstructionRefusal() {
  const auto table = ph::SspYieldTable::loadFromString(narrowYieldTable(), "narrow");
  const auto assessment = ph::assessSspYieldCoverage(tableConfig(ph::SspYieldCoveragePolicy::kReject), table);
  CHUI_REQUIRE(assessment.channels[1].configured_mass_range_uncovered);  // CC 8..40 vs 10..30
  CHUI_REQUIRE(assessment.channels[2].single_metallicity);               // Ia single Z
  CHUI_REQUIRE(!assessment.channels[2].has_mass_axis);                   // no Ia mass-domain check
  CHUI_REQUIRE(!assessment.strict_refusal.empty());
  CHUI_REQUIRE(throwsUnavailable([&] { (void)ph::SspModel(tableConfig(ph::SspYieldCoveragePolicy::kReject), &table); }));
  // Research policy: admitted; out-of-grid use is the counted edge row.
  const ph::SspModel research(tableConfig(ph::SspYieldCoveragePolicy::kEdgeClampResearch), &table);
  const auto cc = research.massiveReturn(1, 1.0e6, 5.0e7, birthMixture(0.004));
  CHUI_REQUIRE(cc.mass_fraction > 0.0);
  CHUI_REQUIRE(research.table_clamps > 0U);
  CHUI_REQUIRE(research.table_mass_clamps > 0U);
  CHUI_REQUIRE(research.table_metallicity_clamps > 0U);
}

void testYieldStrictAtUseRefusal() {
  // AGB axis starts at 1 Msun: admitted at construction (CC fully covered,
  // AGB upper end covered), lower AGB masses enforced at use.
  const auto table = ph::SspYieldTable::loadFromString(coveringYieldTable("1"), "covering");
  const ph::SspModel strict(tableConfig(ph::SspYieldCoveragePolicy::kReject), &table);
  const auto assessment = ph::assessSspYieldCoverage(strict.config(), table);
  CHUI_REQUIRE(assessment.strict_refusal.empty());
  CHUI_REQUIRE(assessment.channels[0].lower_mass_enforced_at_use);
  // In-grid use: no clamp, no throw.
  const auto young = strict.massiveReturn(1, 1.0e6, 5.0e7, birthMixture(0.01));
  CHUI_REQUIRE(young.mass_fraction > 0.0 && strict.table_clamps == 0U);
  // Turnoff below 1 Msun after ~1e10 yr: refused, never edge-clamped.
  CHUI_REQUIRE(throwsUnavailable([&] { (void)strict.massiveReturn(0, 1.0e6, 1.37e10, birthMixture(0.01)); }));
  // Birth metallicity outside 0.001..0.02: refused for CC and for Ia.
  CHUI_REQUIRE(throwsUnavailable([&] { (void)strict.massiveReturn(1, 1.0e6, 5.0e7, birthMixture(0.05)); }));
  CHUI_REQUIRE(throwsUnavailable([&] { (void)strict.iaReturnPerEvent(birthMixture(0.05)); }));
  CHUI_REQUIRE(strict.table_clamps == 0U);  // refusals are not counted as accepted clamps
  // Separate-axis interpolation flags.
  bool mass_clamped = false, z_clamped = false;
  (void)table.massive(1, 50.0, 0.01, mass_clamped, z_clamped);
  CHUI_REQUIRE(mass_clamped && !z_clamped);
}

const char* k_v2_head =
    "# chui_ssp_wind_table_v2\n# source = synthetic analytic test fixture (not a physical table)\n"
    "# license = CC0 test fixture\n# normalization = per_msun_initial_stellar_mass\n"
    "# energy_unit = erg_per_msun\n# momentum_unit = msun_km_s_per_msun\n"
    "# metallicity_definition = metal_mass_fraction\n";
const char* k_two_z_body =
    "metallicities 2 0.001 0.02\nages 2\n6.0 1.0e47 10.0 2.0e47 20.0\n7.0 4.0e47 40.0 8.0e47 80.0\n";
const char* k_one_z_body = "metallicities 1 0.014\nages 2\n6.0 1.0 1.0\n7.0 2.0 2.0\n";

std::string v2Table(const std::string& imf, const std::string& lo, const std::string& hi, const char* body) {
  return std::string(k_v2_head) + "# imf = " + imf + "\n# imf_mass_min_msun = " + lo + "\n# imf_mass_max_msun = " +
      hi + "\n" + body;
}

bool loadRefused(const std::string& text) {
  try {
    (void)ph::StellarWindTable::loadFromString(text, "fixture");
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

ph::StellarWindCompatibilityReport assess(const ph::StellarWindTable& table, ph::StellarWindCompatibilityPolicy policy) {
  return ph::assessStellarWindCompatibility(table, ph::k_chui_ssp_imf_identity, 0.1, 100.0, policy);
}

void testWindV2Compatible() {
  const auto table = ph::StellarWindTable::loadFromString(
      v2Table(ph::k_chui_ssp_imf_identity, "0.1", "100", k_two_z_body), "v2");
  CHUI_REQUIRE(table.metadata().schema == ph::StellarWindTableSchema::kV2);
  CHUI_REQUIRE(table.metadata().imf_identity == ph::k_chui_ssp_imf_identity);
  const auto strict = assess(table, ph::StellarWindCompatibilityPolicy::kStrict);
  CHUI_REQUIRE(strict.admitted && strict.imf_identity_matches && strict.imf_mass_limits_match);
  CHUI_REQUIRE(strict.approximations == "none");
  CHUI_REQUIRE(table.interval(1.0e6, 2.0e6, 0.5).metallicity_clamped);  // at-use refusal is the runtime's job
}

void testWindStrictRefusalsAndResearchAcknowledgement() {
  const auto kroupa = ph::StellarWindTable::loadFromString(v2Table("kroupa2001", "0.1", "100", k_two_z_body), "k");
  CHUI_REQUIRE(!assess(kroupa, ph::StellarWindCompatibilityPolicy::kStrict).admitted);
  const auto kroupa_research = assess(kroupa, ph::StellarWindCompatibilityPolicy::kApproximateResearch);
  CHUI_REQUIRE(kroupa_research.admitted && kroupa_research.approximations.find("kroupa2001") != std::string::npos);

  const auto limits = ph::StellarWindTable::loadFromString(
      v2Table(ph::k_chui_ssp_imf_identity, "0.08", "120", k_two_z_body), "limits");
  const auto limits_strict = assess(limits, ph::StellarWindCompatibilityPolicy::kStrict);
  CHUI_REQUIRE(!limits_strict.admitted && limits_strict.imf_identity_matches && !limits_strict.imf_mass_limits_match);

  const auto unknown = ph::StellarWindTable::loadFromString(v2Table("unknown", "0.1", "100", k_two_z_body), "u");
  CHUI_REQUIRE(!assess(unknown, ph::StellarWindCompatibilityPolicy::kStrict).admitted);

  const auto single = ph::StellarWindTable::loadFromString(
      v2Table(ph::k_chui_ssp_imf_identity, "0.1", "100", k_one_z_body), "single");
  const auto single_strict = assess(single, ph::StellarWindCompatibilityPolicy::kStrict);
  CHUI_REQUIRE(!single_strict.admitted && single_strict.single_metallicity);

  // v1 (SB99-converter-shaped header): loads unchanged, IMF undeclared, the
  // free-text limitation is NOT metadata; strict refuses, research names it.
  const std::string v1 =
      "# chui_ssp_wind_table_v1\n# source = synthetic fixture\n# license = CC0 test fixture\n"
      "# converter = test\n# quantities = cumulative wind energy and momentum\n"
      "# limitation = fixture IMF differs from the CHUI Chabrier SSP IMF; single metallicity\n" +
      std::string(k_one_z_body);
  const auto legacy = ph::StellarWindTable::loadFromString(v1, "v1");
  CHUI_REQUIRE(legacy.metadata().schema == ph::StellarWindTableSchema::kV1);
  CHUI_REQUIRE(legacy.metadata().imf_identity.empty());
  CHUI_REQUIRE(!legacy.limitations().empty());
  const auto legacy_strict = assess(legacy, ph::StellarWindCompatibilityPolicy::kStrict);
  CHUI_REQUIRE(!legacy_strict.admitted && !legacy_strict.imf_declared);
  const auto legacy_research = assess(legacy, ph::StellarWindCompatibilityPolicy::kApproximateResearch);
  CHUI_REQUIRE(legacy_research.admitted && legacy_research.approximations.find("undeclared") != std::string::npos);
}

void testWindMalformedDeclarations() {
  const std::string good = v2Table(ph::k_chui_ssp_imf_identity, "0.1", "100", k_two_z_body);
  CHUI_REQUIRE(!loadRefused(good));
  // Missing IMF.
  CHUI_REQUIRE(loadRefused(std::string(k_v2_head) + "# imf_mass_min_msun = 0.1\n# imf_mass_max_msun = 100\n" +
                           k_two_z_body));
  // Duplicate (contradictory) IMF declarations.
  CHUI_REQUIRE(loadRefused(good + "# imf = kroupa2001\n"));
  // Malformed values.
  CHUI_REQUIRE(loadRefused(v2Table("Chabrier 2003", "0.1", "100", k_two_z_body)));
  CHUI_REQUIRE(loadRefused(v2Table(ph::k_chui_ssp_imf_identity, "0.1", "0.05", k_two_z_body)));
  CHUI_REQUIRE(loadRefused(v2Table(ph::k_chui_ssp_imf_identity, "nan", "100", k_two_z_body)));
  CHUI_REQUIRE(loadRefused(v2Table(ph::k_chui_ssp_imf_identity, "0.1x", "100", k_two_z_body)));
  std::string bad_norm = good;
  bad_norm.replace(bad_norm.find("per_msun_initial_stellar_mass"), 29, "per_star");
  CHUI_REQUIRE(loadRefused(bad_norm));
  // v2 key inside a v1 table; two schema headers; unknown schema version.
  CHUI_REQUIRE(loadRefused("# chui_ssp_wind_table_v1\n# source = x\n# imf = kroupa2001\n" + std::string(k_one_z_body)));
  CHUI_REQUIRE(loadRefused("# chui_ssp_wind_table_v1\n" + good));
  CHUI_REQUIRE(loadRefused("# chui_ssp_wind_table_v3\n# source = x\n" + std::string(k_one_z_body)));
}

}  // namespace

int main() {
  testYieldStrictConstructionRefusal();
  testYieldStrictAtUseRefusal();
  testWindV2Compatible();
  testWindStrictRefusalsAndResearchAcknowledgement();
  testWindMalformedDeclarations();
  return 0;
}
