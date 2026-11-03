#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// CHUI Part 2 thermochemical network manifest.
//
// Representation
// --------------
// A network state is a vector of NUMBER DENSITIES n_i [cm^-3, proper] over
// the stored species (ChemicalSpecies order). Free electrons are derived by
// charge neutrality, n_e = sum_i q_i n_i (>= 0 enforced by the solver), and
// neutral atomic oxygen in the detailed network is derived as
// n_O = n_O,gas - n_CO. Every reaction carries an explicit net stoichiometric
// change over 12 slots: the stored species, the electron and implicit atomic
// oxygen. validateNetworkManifest() proves for every reaction that the net
// change conserves H, He, C, O nuclei and electric charge exactly (integer
// arithmetic), so the stiff solver conserves them to round-off (backward
// Euler preserves linear invariants of the RHS).
//
// Rate provenance labels
// ----------------------
// Every rate cites its publication. Coefficients were TRANSCRIBED from the
// cited papers for this source pass and must be independently re-checked by
// the owner against the originals (NOT RUN — USER-RESERVED). Where a channel
// is not implemented it is NAMED in `omitted_channels` — the network does not
// claim it.
// ---------------------------------------------------------------------------

inline constexpr std::size_t k_network_slot_count = k_chemical_species_count + 2U;
inline constexpr std::size_t k_slot_electron = k_chemical_species_count;      // 10
inline constexpr std::size_t k_slot_atomic_oxygen = k_chemical_species_count + 1U;  // 11

enum class ThermochemistryBackend : std::uint8_t {
  // Pre-Part-2 in-solver cooling source (legacy analytic reference law),
  // kept for configuration compatibility; F2-safe integrator.
  kLegacyCooling = 0,
  // REFERENCE ONLY: Newtonian relaxation du/dt = -(u - u_eq)/t_relax + H_0,
  // closed-form solution, for numerical verification.
  kAnalyticReference = 1,
  // Ionization-equilibrium H/He closure (KWH96 method), equilibrium cooling.
  kEquilibriumPrimordial = 2,
  // Non-equilibrium H, H+, H-, H2, H2+, He, He+, He++ network + metal cooling.
  kFixedHybrid = 3,
  // Hybrid + reduced cold C+/CO network (Nelson & Langer 1997 type).
  kFixedDetailed = 4,
};

[[nodiscard]] std::string_view thermochemistryBackendLabel(ThermochemistryBackend backend) noexcept;
[[nodiscard]] ThermochemistryBackend thermochemistryBackendFromLabel(std::string_view label);

// Physical environment of one rate evaluation (proper cgs units).
struct RateEnvironment {
  double temperature_k = 1.0e4;
  double hydrogen_nuclei_cm3 = 1.0;       // n_H (all chemical forms)
  double oxygen_gas_nuclei_cm3 = 0.0;     // gas-phase O nuclei (after depletion)
  double carbon_gas_nuclei_cm3 = 0.0;     // gas-phase C nuclei (after depletion)
  // Shielded photo-rates [s^-1] and photoheating per photoionization [erg].
  double gamma_hi_s = 0.0;
  double gamma_hei_s = 0.0;
  double gamma_heii_s = 0.0;
  double heat_hi_erg_s = 0.0;   // eps_HI = q_HI (erg s^-1 per HI atom)
  double heat_hei_erg_s = 0.0;
  double heat_heii_erg_s = 0.0;
  // FUV interstellar radiation field, Habing units (G0 = 1: 1.6e-3 erg cm^-2 s^-1).
  double g0_habing = 0.0;
  double visual_extinction_mag = 0.0;
  double h2_column_cm2 = 0.0;
  // Dust-to-gas ratio relative to the Milky Way value (0.01).
  double dust_to_gas_relative_mw = 0.0;
  double redshift = 0.0;
};

using NetworkState = std::array<double, k_chemical_species_count>;

// Rate function: returns the volumetric reaction rate R [cm^-3 s^-1].
using ReactionRateFunction = double (*)(const RateEnvironment&, const NetworkState&, double n_e);

struct ReactionDefinition {
  std::string_view id;
  std::string_view equation;
  std::string_view provenance;
  // Net change per reaction event over k_network_slot_count slots.
  std::array<std::int8_t, k_network_slot_count> net_change{};
  ReactionRateFunction rate = nullptr;
  double valid_temperature_min_k = 1.0;
  double valid_temperature_max_k = 1.0e9;
};

struct NetworkManifest {
  ThermochemistryBackend backend = ThermochemistryBackend::kFixedHybrid;
  std::string name;
  std::string version;
  std::size_t species_count = 0;  // stored species slots [0, species_count)
  std::vector<ReactionDefinition> reactions;
  std::vector<std::string> cooling_channels;
  std::vector<std::string> omitted_channels;
  std::string regime;  // validity statement
};

// Manifests are immutable singletons: one catalog per model, never per cell.
[[nodiscard]] const NetworkManifest& hybridLightNetworkManifest();
[[nodiscard]] const NetworkManifest& detailedCarbonOxygenNetworkManifest();

struct NetworkValidationReport {
  bool ok = true;
  std::string diagnostic;
};
[[nodiscard]] NetworkValidationReport validateNetworkManifest(const NetworkManifest& manifest);

[[nodiscard]] double electronDensity(const NetworkState& n, std::size_t species_count) noexcept;

// ---------------------------------------------------------------------------
// Rate coefficients (cm^3 s^-1 unless noted). Exposed for unit tests.
// ---------------------------------------------------------------------------
namespace rates {
// Katz, Weinberg & Hernquist (1996, ApJS 105, 19), Table 2 (after Cen 1992).
double collisionalIonizationHI(double t);
double collisionalIonizationHeI(double t);
double collisionalIonizationHeII(double t);
double recombinationHII(double t);        // case A
double recombinationHeII(double t);       // radiative, case A
double dielectronicHeII(double t);
double recombinationHeIII(double t);      // case A
// Galli & Palla (1998, A&A 335, 403), Table 1/3.
double formationHminus(double t);         // H + e -> H- + gamma
double associativeDetachment(double t);   // H- + H -> H2 + e
double formationH2plus(double t);         // H + H+ -> H2+ + gamma
double mutualNeutralization(double t);    // H- + H+ -> 2H
// Karpas, Anicich & Huntress (1979, J. Chem. Phys. 70, 2877).
double chargeTransferH2plus(double t);    // H2+ + H -> H2 + H+
// Donahue & Shull (1991, ApJ 383, 511) as listed by Abel et al. (1997).
double electronDissociationH2(double t);  // H2 + e -> 2H + e
// Dove & Mandy (1986, ApJ 311, L93) low-density limit as listed by Abel et
// al. (1997, NewA 2, 181), k13.
double hydrogenDissociationH2(double t);  // H2 + H -> 3H
// Abel et al. (1997) k18 (Schneider et al. 1994 fit).
double dissociativeRecombinationH2plus(double t);  // H2+ + e -> 2H
// Jura (1975, ApJ 197, 575): R_dust = 3e-17 (D/D_MW) n_H n_HI cm^-3 s^-1.
double dustH2FormationCoefficient(double dust_to_gas_relative_mw);
// Draine & Bertoldi (1996, ApJ 468, 269) H2 self-shielding function.
double h2SelfShielding(double h2_column_cm2, double temperature_k);
}  // namespace rates

// ---------------------------------------------------------------------------
// Volumetric cooling/heating terms [erg cm^-3 s^-1] for a network state.
// Radiative losses leave the gas; chemical terms change the chemical
// potential reservoir (ionization / dissociation energy) and are ledgered
// separately so no reaction heat is counted twice.
// ---------------------------------------------------------------------------
struct CoolingHeatingTerms {
  double radiative_cooling = 0.0;  // excitation, recombination, ff, Compton, H2, metals
  double chemical_cooling = 0.0;   // collisional ionization + dissociation (to binding reservoir)
  double heating = 0.0;            // photoheating, photoelectric, H2 formation
  double metal_cooling = 0.0;      // subset of radiative_cooling (diagnostic)
  [[nodiscard]] double net() const noexcept { return heating - radiative_cooling - chemical_cooling; }
};

namespace cooling {
// KWH96 Table 1 primordial terms, Hollenbach & McKee (1979) LTE + Galli &
// Palla (1998) low-density H2 cooling, Compton cooling off the CMB
// (T_CMB = 2.7255 (1+z) K, Fixsen 2009).
CoolingHeatingTerms primordial(const RateEnvironment& env, const NetworkState& n,
                               std::size_t species_count, double n_e);
// Fine-structure cooling of C+ (158 um) and O I (63 um), two-level atoms
// with collision rates of Hollenbach & McKee (1989, ApJ 342, 306) Table 8.
double fineStructureCII(double t, double n_cii, double n_h, double n_e);
double fineStructureOI(double t, double n_oi, double n_h);
// Photoelectric heating, Bakes & Tielens (1994, ApJ 427, 822) as in
// Wolfire et al. (2003, ApJ 587, 278) eq. 19-20, scaled by D/D_MW.
double photoelectricHeating(double t, double n_h, double n_e, double g0, double d_rel);
double comptonCoolingCoefficient(double redshift);  // erg s^-1 K^-1 per electron
}  // namespace cooling

}  // namespace cosmosim::physics
