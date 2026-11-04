// CHUI Part 2.5 — independent rate and cooling reference cases
// (NOT RUN — USER-RESERVED).
//
// Every expected number below was computed OUTSIDE the C++ code base, by a
// separate script that implements the PUBLISHED equations as printed in the
// cited papers (docs/chemistry_rate_provenance.md, "Independent reference
// values"). A mismatch therefore means the implementation no longer matches
// the publication, not merely that it changed. Tolerances: 1e-9 relative for
// closed-form transcriptions (double rounding only); physical cross-checks
// state their own tolerance.
//
// Sources:
//  KWH96  Katz, Weinberg & Hernquist 1996, ApJS 105, 19, Tables 1-2, eq. 23
//  GP98   Galli & Palla 1998, A&A 335, 403, Table 1 (H3, H5, H7, H8)
//  A97    Abel et al. 1997, New Astron. 2, 181 (k10, k12, k18)
//  GA08   Glover & Abel 2008, MNRAS 388, 1627, Table A1 reaction 9, eqs. 13-17
//  F92    Ferland et al. 1992 case A/B fits, as printed in GA08 Table A1
//  DB96   Draine & Bertoldi 1996, ApJ 468, 269, eq. 37
//  W03    Wolfire et al. 2003, ApJ 587, 278, eqs. 23-24 (phi_PAH = 0.5)
//  HM89   Hollenbach & McKee 1989 fine-structure data as tabulated by KROME
//         (Grassi et al. 2014, Table B)
#include <cmath>
#include <string>

#include "cosmosim/physics/chemical_network.hpp"
#include "../support/test_require.hpp"

namespace {

using namespace cosmosim::physics;
using S = ChemicalSpecies;
constexpr std::size_t idx(S s) { return static_cast<std::size_t>(s); }

bool near(double value, double expected, double rtol) {
  return std::abs(value - expected) <= rtol * std::abs(expected);
}
constexpr double k_exact = 1.0e-9;

void testKwh96Table2() {
  CHUI_REQUIRE(near(rates::recombinationHII(1.0e4), 5.0971216268e-13, k_exact));
  CHUI_REQUIRE(near(rates::recombinationHeIII(1.0e4), 2.0388486507e-12, k_exact));
  CHUI_REQUIRE(near(rates::recombinationHeII(1.0e4), 4.3141104433e-13, k_exact));
  CHUI_REQUIRE(near(rates::dielectronicHeII(1.0e5), 6.1051443167e-13, k_exact));
  CHUI_REQUIRE(near(rates::collisionalIonizationHI(2.0e4), 2.1397210736e-12, k_exact));
  CHUI_REQUIRE(near(rates::collisionalIonizationHeI(2.0e4), 1.4810844953e-15, k_exact));
  CHUI_REQUIRE(near(rates::collisionalIonizationHeII(2.0e4), 1.0743770256e-23, k_exact));
  CHUI_REQUIRE(near(rates::collisionalIonizationHI(1.0e5), 1.9088405716e-09, k_exact));
  CHUI_REQUIRE(near(rates::collisionalIonizationHeI(1.0e5), 2.1694570885e-10, k_exact));
  CHUI_REQUIRE(near(rates::collisionalIonizationHeII(1.0e5), 1.6243646986e-12, k_exact));
}

void testCaseAversusCaseB() {
  // F92 fits (independent of KWH96): case A 4.297e-13, case B 2.592e-13.
  const double a = rates::recombinationHIICaseAFerland1992(1.0e4);
  const double b = rates::recombinationHIICaseBFerland1992(1.0e4);
  CHUI_REQUIRE(near(a, 4.2969509300e-13, k_exact));
  CHUI_REQUIRE(near(b, 2.5918095507e-13, k_exact));
  // The KWH96 fit used by CHUI is case-A-like: within 20% of case A and far
  // (> 90%) above case B. If someone silently swapped to case B this fails.
  const double kwh = rates::recombinationHII(1.0e4);
  CHUI_REQUIRE(std::abs(kwh / a - 1.0) < 0.20);
  CHUI_REQUIRE(kwh / b > 1.9);
}

void testGalliPalla1998() {
  CHUI_REQUIRE(near(rates::formationHminus(1.0e3), 8.0042332219e-16, k_exact));      // H3
  CHUI_REQUIRE(near(rates::associativeDetachment(100.0), 1.5e-9, k_exact));         // H5, T <= 300 K
  CHUI_REQUIRE(near(rates::associativeDetachment(1.0e3), 1.2361181730e-09, k_exact));
  CHUI_REQUIRE(near(rates::formationH2plus(1.0e3), 4.8607132352e-18, k_exact));      // H8
  CHUI_REQUIRE(near(rates::mutualNeutralization(1.0e3), 2.4078053118e-07, k_exact)); // H7
  // H8 published validity 1-32000 K: the value at 3.2e4 K is held beyond.
  CHUI_REQUIRE(near(rates::formationH2plus(3.2e4), 7.0225570932e-16, k_exact));
  CHUI_REQUIRE(rates::formationH2plus(1.0e6) == rates::formationH2plus(3.2e4));
}

void testAbel1997() {
  CHUI_REQUIRE(near(rates::chargeTransferH2plus(500.0), 6.4e-10, k_exact));            // k10
  CHUI_REQUIRE(near(rates::electronDissociationH2(3.0e3), 5.0441202456e-24, k_exact));  // k12
  CHUI_REQUIRE(near(rates::dissociativeRecombinationH2plus(100.0), 1.0e-8, k_exact));   // k18, T < 617 K
  CHUI_REQUIRE(near(rates::dissociativeRecombinationH2plus(1.0e3), 6.9274584753e-09, k_exact));
}

void testGloverAbel2008Reaction9() {
  // Limits (GA08 Table A1 reaction 9).
  CHUI_REQUIRE(near(rates::collisionalDissociationH2ByHLowDensity(3.0e3), 8.3630042370e-20, k_exact));
  CHUI_REQUIRE(near(rates::collisionalDissociationH2ByHLte(3.0e3), 1.5536984947e-15, k_exact));
  // Critical densities (eqs. 14-17) for x_H = 0.9, x_H2 = 0.1, x_He = 0.08.
  CHUI_REQUIRE(near(rates::h2DissociationCriticalDensity(3.0e3, 1.0, 0.0, 0.0), 1.3431207062e+03, k_exact));
  CHUI_REQUIRE(near(rates::h2DissociationCriticalDensity(3.0e3, 0.0, 1.0, 0.0), 9.2822637427e+05, k_exact));
  CHUI_REQUIRE(near(rates::h2DissociationCriticalDensity(3.0e3, 0.0, 0.0, 1.0), 1.0392638956e+05, k_exact));
  CHUI_REQUIRE(near(rates::h2DissociationCriticalDensity(3.0e3, 0.9, 0.1, 0.08), 1.4904045741e+03, k_exact));
  // Log-space blend (eq. 13) at three densities.
  CHUI_REQUIRE(near(rates::collisionalDissociationH2ByH(3.0e3, 1.0e-2, 0.9, 0.1, 0.08), 8.3635558213e-20, 1.0e-8));
  CHUI_REQUIRE(near(rates::collisionalDissociationH2ByH(3.0e3, 1.0e3, 0.9, 0.1, 0.08), 4.3305534887e-18, 1.0e-8));
  CHUI_REQUIRE(near(rates::collisionalDissociationH2ByH(3.0e3, 1.0e8, 0.9, 0.1, 0.08), 1.5534708933e-15, 1.0e-8));
  // Physical behaviour that the removed A97 k13 transcription violated:
  // monotonically increasing with T in both limits up to 1e5 K.
  double previous = 0.0;
  for (double t = 1.0e3; t <= 1.0e5; t *= 1.5) {
    const double k = rates::collisionalDissociationH2ByH(t, 1.0, 1.0, 0.0, 0.0);
    CHUI_REQUIRE(k > previous);
    previous = k;
  }
  // Cold gas: no underflow to NaN.
  CHUI_REQUIRE(std::isfinite(rates::collisionalDissociationH2ByH(10.0, 1.0, 1.0, 0.0, 0.0)));
  // R16 is present in the manifest and balanced.
  bool found = false;
  for (const auto& r : hybridLightNetworkManifest().reactions) found = found || r.id == "R16";
  CHUI_REQUIRE(found);
  CHUI_REQUIRE(validateNetworkManifest(hybridLightNetworkManifest()).ok);
}

void testH2SelfShieldingThermalDoppler() {
  // b5 = sqrt(k T / m_H) / 1 km/s = 0.9083 at 100 K (the pre-Part-2.5 code
  // used 0.0917, ten times too small).
  CHUI_REQUIRE(near(rates::h2ThermalDopplerB5(100.0), 9.0828328526e-01, 1.0e-9));
  CHUI_REQUIRE(near(rates::h2SelfShielding(1.0e14, 100.0), 6.8006009867e-01, 1.0e-9));
  CHUI_REQUIRE(near(rates::h2SelfShielding(1.0e18, 100.0), 7.5343554122e-04, 1.0e-9));
}

void testWolfirePhotoelectric() {
  // T = 100 K, n_H = 30, n_e = 0.03 cm^-3, G0 = 1.7, D/D_MW = 1, phi_PAH = 0.5.
  CHUI_REQUIRE(k_wolfire2003_phi_pah == 0.5);
  CHUI_REQUIRE(near(cooling::photoelectricHeating(100.0, 30.0, 0.03, 1.7, 1.0), 2.0147474478e-24, 1.0e-9));
}

void testFineStructure() {
  // Two-level populations with HM89/KROME data; O I Delta E = 227.7 K.
  CHUI_REQUIRE(near(cooling::fineStructureOI(100.0, 1.0e-3, 100.0), 1.7801196802e-26, 1.0e-9));
  CHUI_REQUIRE(near(cooling::fineStructureCII(100.0, 1.0e-4, 100.0, 0.01), 7.8858398320e-26, 1.0e-9));
}

void testKwh96CoolingWithoutCenTriplet() {
  // HeII-only gas at 3e4 K, n_e = 1, n_He+ = 0.1, z = 0: KWH96 Table 1 HeII
  // recombination + dielectronic + excitation + free-free (eq. 23 Gaunt
  // factor) + Compton. The Cen (1992) He 2^3S term (8.99e-30 here, 8e-5
  // relative) is NOT part of KWH96 and must not be present.
  RateEnvironment env;
  env.temperature_k = 3.0e4;
  env.hydrogen_nuclei_cm3 = 1.0;
  env.redshift = 0.0;
  NetworkState n{};
  n[idx(S::kHeII)] = 0.1;
  const auto terms = cooling::primordial(env, n, k_hybrid_species_count, 1.0);
  CHUI_REQUIRE(near(terms.radiative_cooling, 1.0819664145e-25, 1.0e-9));
}

void testProvenanceStatuses() {
  // Gate policy: secondary verification and unverified records never qualify.
  CHUI_REQUIRE(!rateVerificationQualifies(RateVerification::kVerifiedAgainstSecondarySource));
  CHUI_REQUIRE(!rateVerificationQualifies(RateVerification::kUnverifiedTranscription));
  CHUI_REQUIRE(rateVerificationQualifies(RateVerification::kVerifiedAgainstSource));
  CHUI_REQUIRE(!hybridLightNetworkManifest().qualified());     // R18, R13/R15/R17 secondary only (D4)
  CHUI_REQUIRE(!detailedCarbonOxygenNetworkManifest().qualified());
  // The equilibrium rate set is primary-verified (or derived) end to end;
  // the backend nevertheless stays gated until USER-RESERVED validation.
  for (const auto& c : equilibriumPrimordialManifest().cooling_records) {
    const std::string id(c.id);
    if (id == "C01" || id == "C02" || id == "C03") CHUI_REQUIRE(rateVerificationQualifies(c.verification));
    if (id == "C10") CHUI_REQUIRE(c.verification == RateVerification::kVerifiedAgainstSecondarySource);
  }
  for (const auto& r : equilibriumPrimordialManifest().reactions) {
    CHUI_REQUIRE(r.verification == RateVerification::kVerifiedAgainstSource);
  }
  for (const auto& r : hybridLightNetworkManifest().reactions) {
    // D4 (2026-11-04): R18 checked against secondary compilations; the gate stays.
    if (r.id == "R18") CHUI_REQUIRE(r.verification == RateVerification::kVerifiedAgainstSecondarySource);
    CHUI_REQUIRE(!r.provenance.empty());
  }
  NetworkManifest empty;
  CHUI_REQUIRE(!empty.qualified());  // no provenance records => not qualified
}

}  // namespace

int main() {
  testKwh96Table2();
  testCaseAversusCaseB();
  testGalliPalla1998();
  testAbel1997();
  testGloverAbel2008Reaction9();
  testH2SelfShieldingThermalDoppler();
  testWolfirePhotoelectric();
  testFineStructure();
  testKwh96CoolingWithoutCenTriplet();
  testProvenanceStatuses();
  return 0;
}
