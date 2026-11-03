#include "cosmosim/physics/chemical_network.hpp"

#include <algorithm>
#include <cmath>
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
double rH2plusDissRec(const RateEnvironment& e, const NetworkState& n, double ne) {
  return rates::dissociativeRecombinationH2plus(e.temperature_k) * n[idx(S::kH2p)] * ne;
}
double rDustH2(const RateEnvironment& e, const NetworkState& n, double) {
  return rates::dustH2FormationCoefficient(e.dust_to_gas_relative_mw) * e.hydrogen_nuclei_cm3 *
      n[idx(S::kHI)];
}
double rH2Photodissociation(const RateEnvironment& e, const NetworkState& n, double) {
  // Glover & Mac Low (2007, ApJS 169, 239) eq. 32 normalization, Habing G0,
  // DB96 self-shielding and dust attenuation exp(-3.74 A_V) (DB96).
  const double shield = rates::h2SelfShielding(e.h2_column_cm2, e.temperature_k) *
      std::exp(-3.74 * e.visual_extinction_mag);
  return 3.3e-11 * e.g0_habing * shield * n[idx(S::kH2)];
}
// Nelson & Langer (1997, ApJ 482, 796) reduced CO chemistry, CH_x in steady
// state: R = k0 n(C+) n(H2) beta, beta = k1 n(O) / (k1 n(O) + Gamma_CHx),
// k0 = 5e-16, k1 = 5e-10 cm^3 s^-1, Gamma_CHx = 5e-10 G0 s^-1.
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

std::vector<ReactionDefinition> hybridReactions() {
  using std::size_t;
  const size_t HI = idx(S::kHI), HII = idx(S::kHII), Hm = idx(S::kHm), H2 = idx(S::kH2),
               H2p = idx(S::kH2p), HeI = idx(S::kHeI), HeII = idx(S::kHeII), HeIII = idx(S::kHeIII);
  const char* kwh = "Katz, Weinberg & Hernquist 1996 ApJS 105 19, Table 2 (Cen 1992)";
  const char* gp98 = "Galli & Palla 1998 A&A 335 403";
  return {
      {"R01", "H + e -> H+ + 2e", kwh, change({{HI, -1}, {HII, +1}, {E, +1}}), rCollIonHI},
      {"R02", "H+ + e -> H + gamma (case A)", kwh, change({{HII, -1}, {HI, +1}, {E, -1}}), rRecHII},
      {"R03", "He + e -> He+ + 2e", kwh, change({{HeI, -1}, {HeII, +1}, {E, +1}}), rCollIonHeI},
      {"R04", "He+ + e -> He (radiative + dielectronic)", kwh, change({{HeII, -1}, {HeI, +1}, {E, -1}}), rRecHeII},
      {"R05", "He+ + e -> He++ + 2e", kwh, change({{HeII, -1}, {HeIII, +1}, {E, +1}}), rCollIonHeII},
      {"R06", "He++ + e -> He+ (case A)", kwh, change({{HeIII, -1}, {HeII, +1}, {E, -1}}), rRecHeIII},
      {"R07", "H + gamma -> H+ + e (UVB, shielded)", "tabulated UVB (TREECOOL) x Rahmati et al. 2013 MNRAS 430 2427",
       change({{HI, -1}, {HII, +1}, {E, +1}}), rPhotoHI},
      {"R08", "He + gamma -> He+ + e", "tabulated UVB (TREECOOL)", change({{HeI, -1}, {HeII, +1}, {E, +1}}), rPhotoHeI},
      {"R09", "He+ + gamma -> He++ + e", "tabulated UVB (TREECOOL)", change({{HeII, -1}, {HeIII, +1}, {E, +1}}), rPhotoHeII},
      {"R10", "H + e -> H- + gamma", gp98, change({{HI, -1}, {Hm, +1}, {E, -1}}), rHminusForm},
      {"R11", "H- + H -> H2 + e", gp98, change({{Hm, -1}, {HI, -1}, {H2, +1}, {E, +1}}), rAssocDetach},
      {"R12", "H + H+ -> H2+ + gamma", gp98, change({{HI, -1}, {HII, -1}, {H2p, +1}}), rH2plusForm},
      {"R13", "H2+ + H -> H2 + H+", "Karpas, Anicich & Huntress 1979 JChPh 70 2877",
       change({{H2p, -1}, {HI, -1}, {H2, +1}, {HII, +1}}), rH2plusChargeTransfer},
      {"R14", "H- + H+ -> 2H", gp98, change({{Hm, -1}, {HII, -1}, {HI, +2}}), rMutualNeutralization},
      {"R15", "H2 + e -> 2H + e", "Donahue & Shull 1991 ApJ 383 511 (via Abel et al. 1997 k12)",
       change({{H2, -1}, {HI, +2}}), rH2ElectronDissociation},
      // R16 (H2 + H -> 3H) REMOVED by the corrective gate: the transcribed
      // Dove & Mandy / Abel et al. k13 form decreased with temperature above
      // 1e4 K and vanished by 1e6 K (physically implausible; transcription
      // defect). It is listed in omitted_channels until the coefficient is
      // re-derived from the original publication.
      {"R17", "H2+ + e -> 2H", "Schneider et al. 1994 fit (Abel et al. 1997 k18)",
       change({{H2p, -1}, {E, -1}, {HI, +2}}), rH2plusDissRec},
      {"R18", "H + H (grain) -> H2", "Jura 1975 ApJ 197 575, scaled by D/D_MW",
       change({{HI, -2}, {H2, +1}}), rDustH2},
      {"R19", "H2 + gamma(FUV) -> 2H", "Glover & Mac Low 2007 ApJS 169 239; Draine & Bertoldi 1996 ApJ 468 269",
       change({{H2, -1}, {HI, +2}}), rH2Photodissociation},
  };
}

}  // namespace

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

const NetworkManifest& hybridLightNetworkManifest() {
  static const NetworkManifest manifest = [] {
    NetworkManifest m;
    m.backend = ThermochemistryBackend::kFixedHybrid;
    m.name = "chui_hybrid_light";
    m.version = "1";
    m.species_count = k_hybrid_species_count;
    m.reactions = hybridReactions();
    m.cooling_channels = {
        "KWH96 collisional excitation HI, HeII; collisional ionization HI, HeI, HeII; recombination "
        "HII, HeII, HeIII; dielectronic HeII; bremsstrahlung (Gaunt factor of Spitzer & Hart 1979 fit)",
        "Compton cooling/heating off the CMB (cosmological modes only)",
        "H2 rovibrational: Hollenbach & McKee 1979 LTE + Galli & Palla 1998 low-density limit",
        "metal cooling: element table (BLOCKED BY DATA unless supplied) and/or C+ 158um + OI 63um fine structure",
        "photoheating (UVB table), photoelectric heating (Bakes & Tielens 1994, needs G0 > 0 and dust)",
        "H2 formation heating 4.48 eV x n/(n+n_cr) (dust and gas-phase) and dissociation cooling",
    };
    m.omitted_channels = {
        "H- photodetachment and H2+ photodissociation (no Lyman-Werner/IR spectrum beyond G0)",
        "H- + e -> H + 2e and H2 + H+ -> H2+ + H (relevant only at T > ~1e4 K where H2 is negligible)",
        "three-body H2 formation (n_H > 1e8 cm^-3 regime, outside target)",
        "collisional dissociation H2 + H -> 3H (removed: transcription defect; H2 destruction in warm/hot "
        "atomic gas is UNDERESTIMATED, mitigated only by H2 + e dissociation)",
        "HD chemistry and HD cooling",
        "metal electron donors other than C+ (Mg, Si, S, Fe photoionization in neutral gas)",
        "grain recombination cooling and grain-assisted recombination",
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
    m.version = "1";
    m.species_count = k_detailed_species_count;
    const std::size_t CII = idx(S::kCII), CO = idx(S::kCO), H2 = idx(S::kH2), HI = idx(S::kHI);
    m.reactions.push_back({"R20", "C+ + H2 + O + e -> CO + 2H (CHx steady state)",
                           "Nelson & Langer 1997 ApJ 482 796",
                           change({{CII, -1}, {H2, -1}, {O, -1}, {E, -1}, {CO, +1}, {HI, +2}}), rCOFormation,
                           10.0, 1.0e4});
    m.reactions.push_back({"R21", "CO + gamma -> C+ + O + e", "Nelson & Langer 1997 ApJ 482 796",
                           change({{CO, -1}, {CII, +1}, {O, +1}, {E, +1}}), rCOPhotodissociation});
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
  if (cooling_verification == RateVerification::kUnverifiedTranscription) return false;
  for (const ReactionDefinition& reaction : reactions) {
    if (reaction.verification == RateVerification::kUnverifiedTranscription) return false;
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
double formationH2plus(double t) {
  // UNQUALIFIED: the 3.2e4 K upper clamp is an implementation choice (not a
  // verified validity limit of the published fit) that stops the log-cubic
  // from diverging; documented in docs/chemistry_rate_provenance.md.
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
double dissociativeRecombinationH2plus(double t) {
  t = tClamp(t);
  return t < 617.0 ? 1.0e-8 : 1.32e-6 * std::pow(t, -0.76);
}
double dustH2FormationCoefficient(double dust_to_gas_relative_mw) {
  return 3.0e-17 * std::max(dust_to_gas_relative_mw, 0.0);
}
double h2SelfShielding(double h2_column_cm2, double temperature_k) {
  // DB96 eq. 37 with b = 3 km/s thermal+turbulent default scaled by
  // sqrt(T/1e3) floor: b5 = b / 1e5 cm/s.
  const double b5 = std::max(std::sqrt(std::max(temperature_k, 1.0) / 1.0e3) * 0.29, 1.0e-2);
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
  constexpr double delta_e_k = 228.0;
  constexpr double g_ratio = 3.0 / 5.0;
  const double c_ul = 9.2e-11 * std::pow(t / 100.0, 0.67) * n_h;
  const double c_lu = c_ul * g_ratio * std::exp(-delta_e_k / t);
  const double upper_fraction = c_lu / (c_lu + c_ul + a_ul);
  return n_oi * upper_fraction * a_ul * k_boltzmann_erg_per_k * delta_e_k;
}

double photoelectricHeating(double t, double n_h, double n_e, double g0, double d_rel) {
  if (!(g0 > 0.0) || !(d_rel > 0.0) || !(n_h > 0.0)) return 0.0;
  t = tClamp(t);
  const double psi = g0 * std::sqrt(t) / std::max(n_e, 1.0e-30);
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
  rad += 5.01e-27 * std::pow(t, -0.1687) * f5 * std::exp(-55338.0 / t) * n_e * n_e * heii;
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
    const double shield = rates::h2SelfShielding(env.h2_column_cm2, t) * std::exp(-3.74 * env.visual_extinction_mag);
    terms.heating += k_h2_photodissociation_heat_erg * 3.3e-11 * env.g0_habing * shield * n[idx(S::kH2)];
  }
  terms.radiative_cooling += rad;
  return terms;
}

}  // namespace cooling

}  // namespace cosmosim::physics
