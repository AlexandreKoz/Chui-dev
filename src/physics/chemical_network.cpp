#include "cosmosim/physics/chemical_network.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cosmosim::physics {
namespace {

using S = ChemicalSpecies;
constexpr std::size_t idx(S s) { return static_cast<std::size_t>(s); }

constexpr double k_ev = 1.602176634e-12;
constexpr double k_h2_binding_erg = 4.48 * k_ev;
// Black & Dalgarno (1977, ApJS 34, 405): ~0.4 eV heat per photodissociation.
constexpr double k_h2_photodissociation_heat_erg = 0.4 * k_ev;
// Ionization potentials [erg] for the chemical binding ledger.
constexpr double k_chi_hi = 13.598 * k_ev;
constexpr double k_chi_hei = 24.587 * k_ev;
constexpr double k_chi_heii = 54.418 * k_ev;

[[nodiscard]] double tClamp(double t) { return std::max(t, 1.0); }
[[nodiscard]] double t5Factor(double t) { return 1.0 / (1.0 + std::sqrt(t / 1.0e5)); }

std::array<std::int8_t, k_network_slot_count> change(
    std::initializer_list<std::pair<std::size_t, int>> entries) {
  std::array<std::int8_t, k_network_slot_count> c{};
  for (const auto& [slot, delta] : entries) c[slot] = static_cast<std::int8_t>(c[slot] + delta);
  return c;
}
constexpr std::size_t E = k_slot_electron;
constexpr std::size_t O = k_slot_atomic_oxygen;

// --- rate functions R [cm^-3 s^-1] ---------------------------------------
double rCollIonHI(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::collisionalIonizationHI(e.temperature_k) * n[idx(S::kHI)] * ne;
}
double rRecHII(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::recombinationHII(e.temperature_k) * n[idx(S::kHII)] * ne;
}
double rCollIonHeI(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::collisionalIonizationHeI(e.temperature_k) * n[idx(S::kHeI)] * ne;
}
double rRecHeII(const RateEnvironment& e, const NetworkState& n, double ne) {
  return (rates::recombinationHeII(e.temperature_k) + rates::dielectronicHeII(e.temperature_k)) *
      n[idx(S::kHeII)] * ne;
}
double rCollIonHeII(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::collisionalIonizationHeII(e.temperature_k) * n[idx(S::kHeII)] * ne;
}
double rRecHeIII(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::recombinationHeIII(e.temperature_k) * n[idx(S::kHeIII)] * ne;
}
double rPhotoHI(const RateEnvironment& e, const NetworkState& n, double) {
  return e.gamma_hi_s * n[idx(S::kHI)];
}
double rPhotoHeI(const RateEnvironment& e, const NetworkState& n, double) {
  return e.gamma_hei_s * n[idx(S::kHeI)];
}
double rPhotoHeII(const RateEnvironment& e, const NetworkState& n, double) {
  return e.gamma_heii_s * n[idx(S::kHeII)];
}
double rHminusForm(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::formationHminus(e.temperature_k) * n[idx(S::kHI)] * ne;
}
double rAssocDetach(const RateEnvironment& e, const NetworkState& n, double) {
  return rates::associativeDetachment(e.temperature_k) * n[idx(S::kHm)] * n[idx(S::kHI)];
}
double rH2plusForm(const RateEnvironment& e, const NetworkState& n, double) {
  return rates::formationH2plus(e.temperature_k) * n[idx(S::kHI)] * n[idx(S::kHII)];
}
double rH2plusChargeTransfer(const RateEnvironment& e, const NetworkState& n, double) {
  return rates::chargeTransferH2plus(e.temperature_k) * n[idx(S::kH2p)] * n[idx(S::kHI)];
}
double rMutualNeutralization(const RateEnvironment& e, const NetworkState& n, double) {
  return rates::mutualNeutralization(e.temperature_k) * n[idx(S::kHm)] * n[idx(S::kHII)];
}
double rH2ElectronDissociation(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::electronDissociationH2(e.temperature_k) * n[idx(S::kH2)] * ne;
}
double rH2CollisionalDissociationH(const RateEnvironment& e, const NetworkState& n, double) {
  // GA08 eq. 14: x_H = n_HI / n, x_H2 = 2 n_H2 / n, x_He = n_He / n with n the
  // hydrogen-nuclei number density.
  const double nh = e.hydrogen_nuclei_cm3;
  if (!(nh > 0.0)) return 0.0;
  const double x_h = n[idx(S::kHI)] / nh;
  const double x_h2 = 2.0 * n[idx(S::kH2)] / nh;
  const double x_he = std::max(e.helium_nuclei_cm3, 0.0) / nh;
  return rates::collisionalDissociationH2ByH(e.temperature_k, nh, x_h, x_h2, x_he) * n[idx(S::kH2)] *
      n[idx(S::kHI)];
}
double rH2plusDissRec(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::dissociativeRecombinationH2plus(e.temperature_k) * n[idx(S::kH2p)] * ne;
}
double rDustH2(const RateEnvironment& e, const NetworkState& n, double) {
  return rates::dustH2FormationCoefficient(e.dust_to_gas_relative_mw) * e.hydrogen_nuclei_cm3 *
      n[idx(S::kHI)];
}
double rH2Photodissociation(const RateEnvironment& e, const NetworkState& n, double) {
  // Glover & Mac Low (2007, ApJS 169, 239) eq. 32 normalization 3.3e-11 G0
  // (Habing), DB96 eq. 37 self-shielding with the thermal Doppler parameter,
  // and dust attenuation exp(-tau_d,1000) with tau_d,1000 = 2e-21 N_H
  // (GML07 eq. 34), i.e. exp(-3.74 A_V) for N_H/A_V = 1.87e21 cm^-2.
  const double shield = rates::h2SelfShielding(e.h2_column_cm2, e.temperature_k) *
      std::exp(-3.74 * e.visual_extinction_mag);
  return 3.3e-11 * e.g0_habing * shield * n[idx(S::kH2)];
}
// Nelson & Langer (1997, ApJ 482, 796) reduced CO chemistry, CH_x in steady
// state: R = k0 n(C+) n(H2) beta, beta = k1 x_O / (k1 x_O + Gamma_CHx / n),
// k0 = 5e-16, k1 = 5e-10 cm^3 s^-1, Gamma_CHx = 5e-10 G0 exp(-2.5 A_V) s^-1
// (verified against the Glover & Clark 2012, MNRAS 421, 9, eqs. 2-5
// restatement, which also corrects typos in NL97 eqs. 18/20; the NL97 text
// itself was not accessible). NL97 assumes fully molecular hydrogen.
double rCOFormation(const RateEnvironment& e, const NetworkState& n, double) {
  const double n_o = std::max(e.oxygen_gas_nuclei_cm3 - n[idx(S::kCO)], 0.0);
  const double k1no = 5.0e-10 * n_o;
  const double gamma_chx = 5.0e-10 * e.g0_habing * std::exp(-2.5 * e.visual_extinction_mag);
  const double beta = (k1no + gamma_chx) > 0.0 ? k1no / (k1no + gamma_chx) : 0.0;
  return 5.0e-16 * n[idx(S::kCII)] * n[idx(S::kH2)] * beta;
}
// NL97: Gamma_CO = 1e-10 G0 exp(-2.5 A_V) s^-1 (dust shielding only; CO and
// H2 line self-shielding are NOT modelled — named in omitted_channels).
double rCOPhotodissociation(const RateEnvironment& e, const NetworkState& n, double) {
  return 1.0e-10 * e.g0_habing * std::exp(-2.5 * e.visual_extinction_mag) * n[idx(S::kCO)];
}

ReactionDefinition withStatus(ReactionDefinition r, RateVerification v, double t_min = 1.0, double t_max = 1.0e9) {
  r.verification = v;
  r.valid_temperature_min_k = t_min;
  r.valid_temperature_max_k = t_max;
  return r;
}
constexpr RateVerification k_primary = RateVerification::kVerifiedAgainstSource;
constexpr RateVerification k_secondary = RateVerification::kVerifiedAgainstSecondarySource;
[[maybe_unused]] constexpr RateVerification k_unverified = RateVerification::kUnverifiedTranscription;

std::vector<ReactionDefinition> primordialIonizationReactions() {
  using std::size_t;
  const size_t HI = idx(S::kHI), HII = idx(S::kHII), HeI = idx(S::kHeI), HeII = idx(S::kHeII),
               HeIII = idx(S::kHeIII);
  return {
      withStatus({"R01", "H + e -> H+ + 2e", "Katz, Weinberg & Hernquist 1996 ApJS 105 19, Table 2 (Black 1981; Cen 1992)",
                  change({{HI, -1}, {HII, +1}, {E, +1}}), rCollIonHI}, k_primary),
      withStatus({"R02", "H+ + e -> H + gamma (case-A-like fit)", "KWH96 Table 2 (Black 1981; Cen 1992)",
                  change({{HII, -1}, {HI, +1}, {E, -1}}), rRecHII}, k_primary),
      withStatus({"R03", "He + e -> He+ + 2e", "KWH96 Table 2 (Black 1981; Cen 1992)",
                  change({{HeI, -1}, {HeII, +1}, {E, +1}}), rCollIonHeI}, k_primary),
      withStatus({"R04", "He+ + e -> He (radiative + dielectronic)",
                  "KWH96 Table 2 (radiative + dielectronic; Black 1981; Cen 1992)",
                  change({{HeII, -1}, {HeI, +1}, {E, -1}}), rRecHeII}, k_primary),
      withStatus({"R05", "He+ + e -> He++ + 2e", "KWH96 Table 2 (Black 1981; Cen 1992)",
                  change({{HeII, -1}, {HeIII, +1}, {E, +1}}), rCollIonHeII}, k_primary),
      withStatus({"R06", "He++ + e -> He+ (case-A-like fit)", "KWH96 Table 2 (Black 1981; Cen 1992)",
                  change({{HeIII, -1}, {HeII, +1}, {E, -1}}), rRecHeIII}, k_primary),
      // Photo-rates are DATA (the user-supplied TREECOOL table, its column
      // contract verified against the Faucher-Giguere 2020 README) times the
      // Rahmati et al. 2013 eqs. 13-14 shielding factor (verified).
      withStatus({"R07", "H + gamma -> H+ + e (UVB, shielded)",
                  "tabulated UVB (TREECOOL) x Rahmati et al. 2013 MNRAS 430 2427 eqs. 13-14",
                  change({{HI, -1}, {HII, +1}, {E, +1}}), rPhotoHI}, k_primary),
      withStatus({"R08", "He + gamma -> He+ + e", "tabulated UVB (TREECOOL) x Rahmati 2013 factor",
                  change({{HeI, -1}, {HeII, +1}, {E, +1}}), rPhotoHeI}, k_primary),
      withStatus({"R09", "He+ + gamma -> He++ + e", "tabulated UVB (TREECOOL) x Rahmati 2013 factor",
                  change({{HeII, -1}, {HeIII, +1}, {E, +1}}), rPhotoHeII}, k_primary),
  };
}

std::vector<ReactionDefinition> hybridReactions() {
  using std::size_t;
  const size_t HI = idx(S::kHI), HII = idx(S::kHII), Hm = idx(S::kHm), H2 = idx(S::kH2), H2p = idx(S::kH2p);
  std::vector<ReactionDefinition> r = primordialIonizationReactions();
  const std::vector<ReactionDefinition> molecular{
      withStatus({"R10", "H + e -> H- + gamma", "Galli & Palla 1998 A&A 335 403, Table 1 H3",
                  change({{HI, -1}, {Hm, +1}, {E, -1}}), rHminusForm}, k_primary),
      withStatus({"R11", "H- + H -> H2 + e", "Galli & Palla 1998 Table 1 H5",
                  change({{Hm, -1}, {HI, -1}, {H2, +1}, {E, +1}}), rAssocDetach}, k_primary),
      // GP98 H8 states validity 1 <= T <= 32000 K; the value is held beyond.
      withStatus({"R12", "H + H+ -> H2+ + gamma", "Galli & Palla 1998 Table 1 H8 (valid 1-32000 K)",
                  change({{HI, -1}, {HII, -1}, {H2p, +1}}), rH2plusForm}, k_primary, 1.0, 3.2e4),
      withStatus({"R13", "H2+ + H -> H2 + H+",
                  "Karpas, Anicich & Huntress 1979 JChPh 70 2877, as listed by Abel et al. 1997 NewA 2 181 k10",
                  change({{H2p, -1}, {HI, -1}, {H2, +1}, {HII, +1}}), rH2plusChargeTransfer}, k_secondary),
      // GP98 H7 ("fit by PAMS"): no validity range stated; the 4.4e-13 T term
      // makes it grow above ~1e5 K where H- is negligible.
      withStatus({"R14", "H- + H+ -> 2H", "Galli & Palla 1998 Table 1 H7",
                  change({{Hm, -1}, {HII, -1}, {HI, +2}}), rMutualNeutralization}, k_primary),
      withStatus({"R15", "H2 + e -> 2H + e", "Donahue & Shull 1991 ApJ 383 511, as listed by Abel et al. 1997 k12",
                  change({{H2, -1}, {HI, +2}}), rH2ElectronDissociation}, k_secondary),
      // R16 restored in Part 2.5 from Glover & Abel 2008 (the printed
      // Abel et al. 1997 k13 form decreases above ~1 eV and is unusable).
      withStatus({"R16", "H2 + H -> 3H",
                  "Glover & Abel 2008 MNRAS 388 1627, Table A1 reaction 9 (v=0: Mac Low & Shull 1986; LTE: Lepp & "
                  "Shull 1983) with the critical-density blend of eqs. 13-17",
                  change({{H2, -1}, {HI, +2}}), rH2CollisionalDissociationH}, k_primary),
      withStatus({"R17", "H2+ + e -> 2H", "Schneider et al. 1994 fit, as listed by Abel et al. 1997 k18",
                  change({{H2p, -1}, {E, -1}, {HI, +2}}), rH2plusDissRec}, k_secondary),
      withStatus({"R18", "H + H (grain) -> H2",
                  "Jura 1975 ApJ 197 575, R = 3e-17 cm^3 s^-1, scaled by D/D_MW; secondary: Le Petit et al. 2006 "
                  "(Meudon PDR, astro-ph/0602150), Gry et al. 2002 FUSE 3.1-4.5e-17; primary not read",
                  change({{HI, -2}, {H2, +1}}), rDustH2}, k_secondary),
      withStatus({"R19", "H2 + gamma(FUV) -> 2H",
                  "Glover & Mac Low 2007 ApJS 169 239 eqs. 32 and 34; Draine & Bertoldi 1996 ApJ 468 269 eq. 37",
                  change({{H2, -1}, {HI, +2}}), rH2Photodissociation}, k_primary),
  };
  r.insert(r.end(), molecular.begin(), molecular.end());
  return r;
}

std::vector<CoolingChannelRecord> primordialCoolingRecords() {
  return {
      {"C01", "collisional excitation HI, HeII; collisional ionization HI, HeI, HeII; recombination HII, HeII, "
              "HeIII; dielectronic HeII; free-free with the Spitzer & Hart 1979 Gaunt factor",
       "Katz, Weinberg & Hernquist 1996 ApJS 105 19, Table 1 and eq. 23", k_primary},
      {"C02", "Compton exchange with the CMB, 4 sigma_T a k T_gamma^4 / (m_e c) (T - T_gamma) n_e, "
              "T_gamma = 2.7255 (1+z) K (cosmological modes only)",
       "CODATA constants; Fixsen 2009 ApJ 707 916 (equivalent to KWH96 5.41e-36 n_e T (1+z)^4 with T0 = 2.7 K)",
       RateVerification::kDerivedFromConstants},
      {"C03", "UVB photoheating q_HI, q_HeI, q_HeII times the Rahmati 2013 shielding factor",
       "tabulated TREECOOL data; Rahmati et al. 2013 eqs. 13-14", k_primary},
  };
}

}  // namespace

std::string_view rateVerificationLabel(RateVerification v) noexcept {
  switch (v) {
    case RateVerification::kUnverifiedTranscription: return "unverified";
    case RateVerification::kDerivedFromConstants: return "derived_from_constants";
    case RateVerification::kVerifiedAgainstSource: return "verified_primary";
    case RateVerification::kVerifiedAgainstSecondarySource: return "verified_secondary";
  }
  return "unknown";
}

std::string_view thermochemistryBackendLabel(ThermochemistryBackend backend) noexcept {
  switch (backend) {
    case ThermochemistryBackend::kLegacyCooling: return "legacy_cooling";
    case ThermochemistryBackend::kAnalyticReference: return "analytic_reference";
    case ThermochemistryBackend::kEquilibriumPrimordial: return "equilibrium_primordial";
    case ThermochemistryBackend::kFixedHybrid: return "fixed_hybrid";
    case ThermochemistryBackend::kFixedDetailed: return "fixed_detailed";
  }
  return "unknown";
}

ThermochemistryBackend thermochemistryBackendFromLabel(std::string_view label) {
  for (const auto backend : {ThermochemistryBackend::kLegacyCooling, ThermochemistryBackend::kAnalyticReference,
                             ThermochemistryBackend::kEquilibriumPrimordial, ThermochemistryBackend::kFixedHybrid,
                             ThermochemistryBackend::kFixedDetailed}) {
    if (thermochemistryBackendLabel(backend) == label) return backend;
  }
  throw std::invalid_argument("unknown thermochemistry backend '" + std::string(label) + "'");
}

namespace {
// Channels shared by every backend that evaluates cooling::primordial plus
// the optional metal / photoelectric terms of ThermochemistrySolver.
void appendOptionalCoolingRecords(std::vector<CoolingChannelRecord>& records) {
  records.push_back({"C10", "C+ 158 um and O I 63 um fine structure (two-level), metal_cooling_mode with fine structure",
                     "Hollenbach & McKee 1989 ApJ 342 306 Table 8, verified against Grassi et al. 2014 (KROME) Table B",
                     k_secondary});
  records.push_back({"C11", "photoelectric heating (needs isrf_g0_habing > 0 and dust)",
                     "Wolfire et al. 2003 ApJ 587 278 eqs. 23-24 with phi_PAH = 0.5", k_primary});
  // The element-resolved metal table is DATA: its provenance (citation,
  // license, solar reference, UV background) is carried and checked by the
  // table header (MetalCoolingTable), not by this rate manifest.
}
}  // namespace

const NetworkManifest& equilibriumPrimordialManifest() {
  static const NetworkManifest manifest = [] {
    NetworkManifest m;
    m.backend = ThermochemistryBackend::kEquilibriumPrimordial;
    m.name = "chui_equilibrium_primordial_kwh96";
    m.version = "1";
    m.species_count = 0;  // ionization-equilibrium closure, no stored species
    m.reactions = primordialIonizationReactions();
    m.cooling_records = primordialCoolingRecords();
    appendOptionalCoolingRecords(m.cooling_records);
    m.cooling_channels = {"KWH96 Table 1 primordial cooling", "Compton (cosmological modes only)",
                          "UVB photoheating", "optional metal / fine-structure / photoelectric terms"};
    m.omitted_channels = {
        "He 2^3S collisional excitation n_e^2 n_He+ term of Cen 1992 (not in KWH96 Table 1; primary not "
        "accessible; removed in Part 2.5)",
        "case-B recombination for self-shielded gas (case-A-like KWH96 fits used everywhere)",
        "molecules (H2, H-, H2+): equilibrium backend is atomic only",
    };
    m.regime = "atomic H/He in ionization equilibrium (KWH96 assumptions), optically thin except the local "
               "Rahmati et al. 2013 approximation; 1e4 K <~ T <= 1e9 K for the cooling fits";
    return m;
  }();
  return manifest;
}

const NetworkManifest& hybridLightNetworkManifest() {
  static const NetworkManifest manifest = [] {
    NetworkManifest m;
    m.backend = ThermochemistryBackend::kFixedHybrid;
    m.name = "chui_hybrid_light";
    m.version = "2";  // Part 2.5: R16 restored, b5 and phi_PAH corrected, Cen 2^3S term removed
    m.species_count = k_hybrid_species_count;
    m.reactions = hybridReactions();
    m.cooling_records = primordialCoolingRecords();
    m.cooling_records.push_back(
        {"C04", "H2 rovibrational cooling: LTE (Hollenbach & McKee 1979) with the n -> 0 limit and critical-density "
                "interpolation of Galli & Palla 1998 eqs. A5-A7",
         "LTE term: Hollenbach & McKee 1979 ApJS 41 555 verified against Grassi et al. 2014 eq. 10; low-density "
         "limit and interpolation: Galli & Palla 1998 Appendix A (primary)",
         k_secondary});
    m.cooling_records.push_back({"C05", "H2 formation heating 4.48 eV x n/(n + n_cr) and dissociation cooling (R15, R16)",
                                 "binding energy 4.48 eV (constants); thermalized fraction follows C04", k_secondary});
    m.cooling_records.push_back({"C06", "H2 photodissociation heating 0.4 eV per event",
                                 "Black & Dalgarno 1977 ApJS 34 405; secondary: Glover & Mac Low 2007 (astro-ph/0605120) adopt "
                                 "0.4 eV per dissociation (some sources attribute it to Stephens & Dalgarno 1973); "
                                 "primary not read", k_secondary});
    appendOptionalCoolingRecords(m.cooling_records);
    m.cooling_channels = {
        "KWH96 collisional excitation HI, HeII; collisional ionization HI, HeI, HeII; recombination "
        "HII, HeII, HeIII; dielectronic HeII; bremsstrahlung (Gaunt factor of Spitzer & Hart 1979 fit)",
        "Compton cooling/heating off the CMB (cosmological modes only)",
        "H2 rovibrational: Hollenbach & McKee 1979 LTE + Galli & Palla 1998 low-density limit",
        "metal cooling: element table (BLOCKED BY DATA unless supplied) and/or C+ 158um + OI 63um fine structure",
        "photoheating (UVB table), photoelectric heating (Wolfire et al. 2003, needs G0 > 0 and dust)",
        "H2 formation heating 4.48 eV x n/(n+n_cr) (dust and gas-phase) and dissociation cooling (R15, R16)",
    };
    m.omitted_channels = {
        "H- photodetachment and H2+ photodissociation (no Lyman-Werner/IR spectrum beyond G0)",
        "H- + e -> H + 2e and H2 + H+ -> H2+ + H (relevant only at T > ~1e4 K where H2 is negligible)",
        "three-body H2 formation (n_H > 1e8 cm^-3 regime, outside target)",
        "H2 collisional dissociation by H2 and He partners (GA08 reactions 10-11); only H2 + H (R16) and "
        "H2 + e (R15) destroy H2 collisionally",
        "He 2^3S collisional excitation n_e^2 n_He+ term of Cen 1992 (not in KWH96 Table 1; removed in Part 2.5)",
        "turbulent Doppler broadening in H2 self-shielding (thermal b only)",
        "case-B recombination for self-shielded gas",
        "HD chemistry and HD cooling",
        "metal electron donors other than C+ (Mg, Si, S, Fe photoionization in neutral gas)",
        "grain recombination cooling (Wolfire et al. 2003 eq. 25) and grain-assisted recombination",
    };
    m.regime = "10 K <= T <= 1e9 K, n_H <= 1e6 cm^-3, optically thin to H/He ionizing photons except "
               "the local Rahmati et al. 2013 self-shielding approximation";
    return m;
  }();
  return manifest;
}

const NetworkManifest& detailedCarbonOxygenNetworkManifest() {
  static const NetworkManifest manifest = [] {
    NetworkManifest m = hybridLightNetworkManifest();
    m.backend = ThermochemistryBackend::kFixedDetailed;
    m.name = "chui_detailed_hybrid_co_nl97";
    m.version = "2";
    m.species_count = k_detailed_species_count;
    const std::size_t CII = idx(S::kCII), CO = idx(S::kCO), H2 = idx(S::kH2), HI = idx(S::kHI);
    m.reactions.push_back(withStatus({"R20", "C+ + H2 + O + e -> CO + 2H (CHx steady state)",
                                      "Nelson & Langer 1997 ApJ 482 796, verified against Glover & Clark 2012 MNRAS "
                                      "421 9 eqs. 2-5",
                                      change({{CII, -1}, {H2, -1}, {O, -1}, {E, -1}, {CO, +1}, {HI, +2}}),
                                      rCOFormation},
                                     k_secondary, 10.0, 1.0e4));
    m.reactions.push_back(withStatus({"R21", "CO + gamma -> C+ + O + e",
                                      "Nelson & Langer 1997, verified against Glover & Clark 2012 eq. 5",
                                      change({{CO, -1}, {CII, +1}, {O, +1}, {E, +1}}), rCOPhotodissociation},
                                     k_secondary));
    m.cooling_channels.push_back("C+ 158 um fine structure from the evolved C+ abundance");
    m.omitted_channels.push_back("neutral atomic carbon (NL97 convention: non-CO carbon is C+)");
    m.omitted_channels.push_back("CO rotational cooling (no validated escape-probability table bundled)");
    m.omitted_channels.push_back("CO and H2 line self-shielding of CO (dust attenuation only)");
    m.regime = "cold/molecular enriched gas, 10 K <= T <= 1e4 K for the C/O channels; hybrid regime otherwise";
    return m;
  }();
  return manifest;
}

bool NetworkManifest::qualified() const noexcept {
  if (cooling_records.empty()) return false;  // no provenance => not qualified
  for (const CoolingChannelRecord& record : cooling_records) {
    if (!rateVerificationQualifies(record.verification)) return false;
  }
  for (const ReactionDefinition& reaction : reactions) {
    if (!rateVerificationQualifies(reaction.verification)) return false;
  }
  return true;
}

NetworkValidationReport validateNetworkManifest(const NetworkManifest& manifest) {
  NetworkValidationReport report;
  for (const ReactionDefinition& reaction : manifest.reactions) {
    int dh = 0, dhe = 0, dc = 0, d_o = 0, dq = 0;
    for (std::size_t s = 0; s < k_chemical_species_count; ++s) {
      const int delta = reaction.net_change[s];
      if (delta == 0) continue;
      if (s >= manifest.species_count) {
        report.ok = false;
        report.diagnostic += std::string(reaction.id) + " touches species outside the network; ";
      }
      const auto& info = chemicalSpeciesInfo(static_cast<ChemicalSpecies>(s));
      dh += delta * info.n_h;
      dhe += delta * info.n_he;
      dc += delta * info.n_c;
      d_o += delta * info.n_o;
      dq += delta * info.charge;
    }
    dq -= reaction.net_change[k_slot_electron];
    d_o += reaction.net_change[k_slot_atomic_oxygen];
    if (dh != 0 || dhe != 0 || dc != 0 || d_o != 0 || dq != 0) {
      report.ok = false;
      report.diagnostic += std::string(reaction.id) + " (" + std::string(reaction.equation) +
          ") violates element/charge balance; ";
    }
    if (reaction.rate == nullptr || reaction.provenance.empty()) {
      report.ok = false;
      report.diagnostic += std::string(reaction.id) + " lacks a rate law or provenance; ";
    }
  }
  return report;
}

double electronDensity(const NetworkState& n, std::size_t species_count) noexcept {
  double ne = 0.0;
  for (std::size_t s = 0; s < species_count; ++s) {
    ne += chemicalSpeciesInfo(static_cast<ChemicalSpecies>(s)).charge * n[s];
  }
  return std::max(ne, 0.0);
}

namespace rates {

double collisionalIonizationHI(double t) {
  t = tClamp(t);
  return 5.85e-11 * std::sqrt(t) * std::exp(-157809.1 / t) * t5Factor(t);
}
double collisionalIonizationHeI(double t) {
  t = tClamp(t);
  return 2.38e-11 * std::sqrt(t) * std::exp(-285335.4 / t) * t5Factor(t);
}
double collisionalIonizationHeII(double t) {
  t = tClamp(t);
  return 5.68e-12 * std::sqrt(t) * std::exp(-631515.0 / t) * t5Factor(t);
}
double recombinationHII(double t) {
  t = tClamp(t);
  return 8.4e-11 / std::sqrt(t) * std::pow(t / 1.0e3, -0.2) / (1.0 + std::pow(t / 1.0e6, 0.7));
}
double recombinationHeII(double t) {
  t = tClamp(t);
  return 1.5e-10 * std::pow(t, -0.6353);
}
double dielectronicHeII(double t) {
  t = tClamp(t);
  return 1.9e-3 * std::pow(t, -1.5) * std::exp(-470000.0 / t) * (1.0 + 0.3 * std::exp(-94000.0 / t));
}
double recombinationHeIII(double t) {
  t = tClamp(t);
  return 3.36e-10 / std::sqrt(t) * std::pow(t / 1.0e3, -0.2) / (1.0 + std::pow(t / 1.0e6, 0.7));
}
double formationHminus(double t) {
  t = tClamp(t);
  return 1.4e-18 * std::pow(t, 0.928) * std::exp(-t / 16200.0);
}
double associativeDetachment(double t) {
  t = tClamp(t);
  return t < 300.0 ? 1.5e-9 : 4.0e-9 * std::pow(t, -0.17);
}
double recombinationHIICaseAFerland1992(double t) {
  t = tClamp(t);
  return 1.269e-13 * std::pow(315614.0 / t, 1.503) * std::pow(1.0 + std::pow(604625.0 / t, 0.470), -1.923);
}
double recombinationHIICaseBFerland1992(double t) {
  t = tClamp(t);
  return 2.753e-14 * std::pow(315614.0 / t, 1.5) * std::pow(1.0 + std::pow(115188.0 / t, 0.407), -2.242);
}
double formationH2plus(double t) {
  // GP98 Table 1 H8 is published as valid for 1 <= T <= 32000 K (verified,
  // Part 2.5). Beyond the upper limit the 3.2e4 K value is held (the cubic
  // in log T would otherwise grow without bound); documented in
  // docs/chemistry_rate_provenance.md, never extrapolated.
  t = std::clamp(t, 1.0, 3.2e4);
  const double lt = std::log10(t);
  return std::pow(10.0, -19.38 - 1.523 * lt + 1.118 * lt * lt - 0.1269 * lt * lt * lt);
}
double mutualNeutralization(double t) {
  t = tClamp(t);
  return std::max(5.7e-6 / std::sqrt(t) + 6.3e-8 - 9.2e-11 * std::sqrt(t) + 4.4e-13 * t, 0.0);
}
double chargeTransferH2plus(double) { return 6.4e-10; }
double electronDissociationH2(double t) {
  t = tClamp(t);
  return 5.6e-11 * std::sqrt(t) * std::exp(-102124.0 / t);
}
double collisionalDissociationH2ByHLowDensity(double t) {
  t = tClamp(t);
  // GA08 Table A1 reaction 9, v=0 limit: 6.67e-12 T^0.5 exp[-(1 + 63593/T)].
  return std::pow(10.0, std::log10(6.67e-12) + 0.5 * std::log10(t) - (1.0 + 63593.0 / t) / std::log(10.0));
}
double collisionalDissociationH2ByHLte(double t) {
  t = tClamp(t);
  // GA08 Table A1 reaction 9, LTE limit: 3.52e-9 exp(-43900/T).
  return std::pow(10.0, std::log10(3.52e-9) - (43900.0 / t) / std::log(10.0));
}
double h2DissociationCriticalDensity(double t, double x_h, double x_h2, double x_he) {
  t = tClamp(t);
  const double lt4 = std::log10(t / 1.0e4);
  const double ncr_h = std::pow(10.0, 3.0 - 0.416 * lt4 - 0.327 * lt4 * lt4);            // GA08 eq. 15
  const double ncr_h2 = std::pow(10.0, 4.845 - 1.3 * lt4 + 1.62 * lt4 * lt4);            // GA08 eq. 16
  const double ncr_he = std::pow(10.0, 5.0792 * (1.0 - 1.23e-5 * (t - 2000.0)));         // GA08 eq. 17
  const double inverse = std::max(x_h, 0.0) / ncr_h + std::max(x_h2, 0.0) / ncr_h2 +      // GA08 eq. 14
      std::max(x_he, 0.0) / ncr_he;
  return inverse > 0.0 ? 1.0 / inverse : std::numeric_limits<double>::infinity();
}
double collisionalDissociationH2ByH(double t, double n_h_nuclei, double x_h, double x_h2, double x_he) {
  t = tClamp(t);
  const double log_v0 = std::log10(6.67e-12) + 0.5 * std::log10(t) - (1.0 + 63593.0 / t) / std::log(10.0);
  const double log_lte = std::log10(3.52e-9) - (43900.0 / t) / std::log(10.0);
  const double ncr = h2DissociationCriticalDensity(t, x_h, x_h2, x_he);
  // GA08 eq. 13 in log space: weight w = (n/ncr)/(1 + n/ncr).
  const double ratio = std::isfinite(ncr) ? std::max(n_h_nuclei, 0.0) / ncr : 0.0;
  const double w = ratio / (1.0 + ratio);
  return std::pow(10.0, w * log_lte + (1.0 - w) * log_v0);
}
double dissociativeRecombinationH2plus(double t) {
  t = tClamp(t);
  return t < 617.0 ? 1.0e-8 : 1.32e-6 * std::pow(t, -0.76);
}
double dustH2FormationCoefficient(double dust_to_gas_relative_mw) {
  return 3.0e-17 * std::max(dust_to_gas_relative_mw, 0.0);
}
double h2ThermalDopplerB5(double temperature_k) {
  // b = sqrt(2 k T / m_H2) = sqrt(k T / m_H); m_H = 1.6735575e-24 g.
  constexpr double k_b = 1.380649e-16;
  constexpr double m_h = 1.6735575e-24;
  return std::sqrt(k_b * std::max(temperature_k, 1.0) / m_h) / 1.0e5;
}
double h2SelfShielding(double h2_column_cm2, double temperature_k) {
  // DB96 eq. 37: f = 0.965/(1 + x/b5)^2 + 0.035/sqrt(1+x) exp[-8.5e-4 sqrt(1+x)],
  // x = N_H2 / 5e14 cm^-2, b5 = b / 1 km/s. Part 2.5 correction: the former
  // b5 = 0.29 sqrt(T/1e3) was 10x smaller than the thermal value 0.0908
  // sqrt(T) (over-shielding by up to two orders of magnitude in the
  // power-law regime). Thermal broadening only.
  const double b5 = h2ThermalDopplerB5(temperature_k);
  const double x = std::max(h2_column_cm2, 0.0) / 5.0e14;
  return 0.965 / std::pow(1.0 + x / b5, 2.0) +
      0.035 / std::sqrt(1.0 + x) * std::exp(-8.5e-4 * std::sqrt(1.0 + x));
}

}  // namespace rates

namespace cooling {

double comptonCoolingCoefficient(double redshift) {
  // 4 sigma_T a k_B T_gamma^4 / (m_e c), T_gamma = 2.7255 (1+z) K.
  constexpr double sigma_t = 6.6524587321e-25;
  constexpr double a_rad = 7.565733250e-15;
  constexpr double k_b = 1.380649e-16;
  constexpr double m_e = 9.1093837015e-28;
  constexpr double c = 2.99792458e10;
  const double t_gamma = 2.7255 * (1.0 + std::max(redshift, 0.0));
  return 4.0 * sigma_t * a_rad * k_b * std::pow(t_gamma, 4.0) / (m_e * c);
}

double fineStructureCII(double t, double n_cii, double n_h, double n_e) {
  t = tClamp(t);
  constexpr double a_ul = 2.4e-6;
  constexpr double delta_e_k = 91.2;
  constexpr double g_ratio = 2.0;  // g_u / g_l = 4 / 2
  const double c_ul = 8.0e-10 * std::pow(t / 100.0, 0.07) * n_h + 2.8e-7 * std::pow(t / 100.0, -0.5) * n_e;
  const double c_lu = c_ul * g_ratio * std::exp(-delta_e_k / t);
  const double upper_fraction = c_lu / (c_lu + c_ul + a_ul);
  return n_cii * upper_fraction * a_ul * k_boltzmann_erg_per_k * delta_e_k;
}

double fineStructureOI(double t, double n_oi, double n_h) {
  t = tClamp(t);
  constexpr double a_ul = 8.9e-5;
  constexpr double delta_e_k = 227.7;  // hc / (k 63.18 um)
  constexpr double g_ratio = 3.0 / 5.0;
  const double c_ul = 9.2e-11 * std::pow(t / 100.0, 0.67) * n_h;
  const double c_lu = c_ul * g_ratio * std::exp(-delta_e_k / t);
  const double upper_fraction = c_lu / (c_lu + c_ul + a_ul);
  return n_oi * upper_fraction * a_ul * k_boltzmann_erg_per_k * delta_e_k;
}

double photoelectricHeating(double t, double n_h, double n_e, double g0, double d_rel) {
  if (!(g0 > 0.0) || !(d_rel > 0.0) || !(n_h > 0.0)) return 0.0;
  t = tClamp(t);
  // Wolfire et al. 2003 eqs. 23-24: psi = G0 T^1/2 / (n_e phi_PAH).
  const double psi = g0 * std::sqrt(t) / (std::max(n_e, 1.0e-30) * k_wolfire2003_phi_pah);
  const double eps = 4.9e-2 / (1.0 + 4.0e-3 * std::pow(psi, 0.73)) +
      3.7e-2 * std::pow(t / 1.0e4, 0.7) / (1.0 + 2.0e-4 * psi);
  return 1.3e-24 * eps * g0 * n_h * d_rel;
}

CoolingHeatingTerms primordial(const RateEnvironment& env, const NetworkState& n,
                               std::size_t species_count, double n_e) {
  CoolingHeatingTerms terms;
  const double t = tClamp(env.temperature_k);
  const double f5 = t5Factor(t);
  const double sqrt_t = std::sqrt(t);
  const double hi = n[idx(S::kHI)], hii = n[idx(S::kHII)];
  const double hei = n[idx(S::kHeI)], heii = n[idx(S::kHeII)], heiii = n[idx(S::kHeIII)];
  // Collisional ionization: energy goes into the ionization (chemical) reservoir.
  terms.chemical_cooling += k_chi_hi * rates::collisionalIonizationHI(t) * hi * n_e;
  terms.chemical_cooling += k_chi_hei * rates::collisionalIonizationHeI(t) * hei * n_e;
  terms.chemical_cooling += k_chi_heii * rates::collisionalIonizationHeII(t) * heii * n_e;
  // KWH96 radiative terms.
  double rad = 0.0;
  // (The Cen 1992 He 2^3S n_e^2 n_He+ term is NOT part of KWH96 Table 1 and
  // was removed in Part 2.5; see omitted_channels.)
  const double t3 = std::pow(t / 1.0e3, -0.2) / (1.0 + std::pow(t / 1.0e6, 0.7));
  rad += 8.70e-27 * sqrt_t * t3 * n_e * hii;
  rad += 1.55e-26 * std::pow(t, 0.3647) * n_e * heii;
  rad += 3.48e-26 * sqrt_t * t3 * n_e * heiii;
  rad += 1.24e-13 * std::pow(t, -1.5) * std::exp(-470000.0 / t) * (1.0 + 0.3 * std::exp(-94000.0 / t)) *
      n_e * heii;
  rad += 7.50e-19 * f5 * std::exp(-118348.0 / t) * n_e * hi;
  rad += 5.54e-17 * std::pow(t, -0.397) * f5 * std::exp(-473638.0 / t) * n_e * heii;
  const double lt = std::log10(t);
  const double gaunt = 1.1 + 0.34 * std::exp(-(5.5 - lt) * (5.5 - lt) / 3.0);
  rad += 1.42e-27 * gaunt * sqrt_t * (hii + heii + 4.0 * heiii) * n_e;
  // Compton exchange with the CMB (signed: heats gas below T_CMB).
  const double t_gamma = 2.7255 * (1.0 + std::max(env.redshift, 0.0));
  const double compton = comptonCoolingCoefficient(env.redshift) * (t - t_gamma) * n_e;
  if (compton >= 0.0) rad += compton; else terms.heating -= compton;
  // Photoheating from the (shielded) UVB.
  terms.heating += hi * env.heat_hi_erg_s + hei * env.heat_hei_erg_s + heii * env.heat_heii_erg_s;
  // H2 rovibrational cooling, critical-density interpolation.
  if (species_count > idx(S::kH2) && n[idx(S::kH2)] > 0.0) {
    const double t3k = t / 1.0e3;
    const double l_rot = 9.5e-22 * std::pow(t3k, 3.76) / (1.0 + 0.12 * std::pow(t3k, 2.1)) *
            std::exp(-std::pow(0.13 / t3k, 3.0)) + 3.0e-24 * std::exp(-0.51 / t3k);
    const double l_vib = 6.7e-19 * std::exp(-5.86 / t3k) + 1.6e-18 * std::exp(-11.7 / t3k);
    const double l_lte = l_rot + l_vib;
    const double lt_c = std::log10(std::clamp(t, 10.0, 1.0e4));
    const double l0 = std::pow(10.0, -103.0 + 97.59 * lt_c - 48.05 * lt_c * lt_c +
                                         10.80 * lt_c * lt_c * lt_c - 0.9032 * lt_c * lt_c * lt_c * lt_c);
    const double nh = std::max(env.hydrogen_nuclei_cm3, 1.0e-30);
    const double per_molecule = l_lte / (1.0 + l_lte / (nh * l0));
    rad += per_molecule * n[idx(S::kH2)];
    // Formation heating / dissociation cooling with the same critical
    // density n_cr = L_LTE / L_0 (fraction of 4.48 eV thermalized).
    const double thermalized = nh / (nh + l_lte / l0);
    const double r_form = rates::dustH2FormationCoefficient(env.dust_to_gas_relative_mw) * nh * hi +
        rates::associativeDetachment(t) * n[idx(S::kHm)] * hi;
    terms.heating += k_h2_binding_erg * thermalized * r_form;
    terms.chemical_cooling += k_h2_binding_erg * rates::electronDissociationH2(t) * n_e * n[idx(S::kH2)];
    {
      // R16 (GA08 reaction 9) removes the same 4.48 eV binding energy.
      const double x_h = hi / nh, x_h2 = 2.0 * n[idx(S::kH2)] / nh;
      const double x_he = std::max(env.helium_nuclei_cm3, 0.0) / nh;
      terms.chemical_cooling += k_h2_binding_erg *
          rates::collisionalDissociationH2ByH(t, nh, x_h, x_h2, x_he) * hi * n[idx(S::kH2)];
    }
    const double shield = rates::h2SelfShielding(env.h2_column_cm2, t) * std::exp(-3.74 * env.visual_extinction_mag);
    terms.heating += k_h2_photodissociation_heat_erg * 3.3e-11 * env.g0_habing * shield * n[idx(S::kH2)];
  }
  terms.radiative_cooling += rad;
  return terms;
}

}  // namespace cooling

}  // namespace cosmosim::physics
