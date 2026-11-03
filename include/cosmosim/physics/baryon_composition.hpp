#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/gas_composition_sidecar.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// Element vector (CHUI Part 2 baseline candidate, not a frozen inventory).
// Conserved owner: GasCompositionSidecar element lanes hold the TOTAL gas
// nuclei of each element in all gas-phase chemical forms (atoms, ions,
// molecules), in code mass units. Dust solids are held in separate dust lanes
// (their element content is fixed by reservoir stoichiometry); stars hold
// their own birth composition. Z (GasCellSidecar::metal_mass_code) is a
// DERIVED compatibility mirror equal to gas-phase + dust metal nuclei
// (everything except H and He) whenever element-resolved chemistry is on.
// ---------------------------------------------------------------------------
enum class ChemicalElement : std::uint8_t {
  kH = 0, kHe, kC, kN, kO, kNe, kMg, kSi, kS, kCa, kFe,
};
inline constexpr std::size_t k_chemical_element_count = 11U;

struct ChemicalElementInfo {
  std::string_view symbol;
  // Standard atomic weight (IUPAC 2009 conventional values), amu.
  double atomic_mass_amu = 0.0;
  // Asplund, Grevesse, Sauval & Scott (2009, ARA&A 47, 481) present-day solar
  // photospheric abundance A(X) = log10(n_X/n_H) + 12 (Table 1).
  double solar_log_epsilon = 0.0;
};

[[nodiscard]] const std::array<ChemicalElementInfo, k_chemical_element_count>& chemicalElementTable();
[[nodiscard]] const ChemicalElementInfo& chemicalElementInfo(ChemicalElement element);
// Solar mass fractions X_e derived from the A09 abundances (sum = 1).
[[nodiscard]] const std::array<double, k_chemical_element_count>& solarElementMassFractions();
[[nodiscard]] double solarMetalMassFraction();

// Physical constants used across the baryonic modules (CODATA 2018).
inline constexpr double k_amu_g = 1.66053906660e-24;
inline constexpr double k_proton_mass_g = 1.67262192369e-24;
inline constexpr double k_electron_mass_g = 9.1093837015e-28;
inline constexpr double k_boltzmann_erg_per_k = 1.380649e-16;
inline constexpr double k_seconds_per_year = 31557600.0;  // Julian year
inline constexpr double k_solar_mass_g = 1.98847e33;
inline constexpr double k_erg_per_ev = 1.602176634e-12;

// ---------------------------------------------------------------------------
// Chemical species catalog shared by every network manifest. Order is the
// persistent species-slot order (layout v1); never reorder, only append.
// Electrons are DERIVED from charge neutrality and never stored or advected.
// In the detailed C/O network all non-CO gas-phase carbon is carried as C+
// (Nelson & Langer 1997 convention: C is assumed photoionized) and neutral
// atomic oxygen is the derived remainder O_total - O(CO).
// ---------------------------------------------------------------------------
enum class ChemicalSpecies : std::uint8_t {
  kHI = 0, kHII, kHm, kH2, kH2p, kHeI, kHeII, kHeIII, kCII, kCO,
};
inline constexpr std::size_t k_chemical_species_count = 10U;
inline constexpr std::size_t k_hybrid_species_count = 8U;    // kHI..kHeIII
inline constexpr std::size_t k_detailed_species_count = 10U; // + kCII, kCO

struct ChemicalSpeciesInfo {
  std::string_view name;
  // Nuclei per particle for the elements that appear in any network species.
  std::uint8_t n_h = 0;
  std::uint8_t n_he = 0;
  std::uint8_t n_c = 0;
  std::uint8_t n_o = 0;
  std::int8_t charge = 0;
  // Particle mass in amu (nuclear content minus electron deficit/excess).
  double mass_amu = 0.0;
};
[[nodiscard]] const ChemicalSpeciesInfo& chemicalSpeciesInfo(ChemicalSpecies species);
[[nodiscard]] inline constexpr double k_electron_mass_amu() { return 5.48579909065e-4; }

// ---------------------------------------------------------------------------
// Composition layout: the physics manifest that gives meaning to the generic
// core::GasCompositionSidecar lanes.
// ---------------------------------------------------------------------------
enum class CompositionNetworkKind : std::uint8_t {
  // Element lanes only (equilibrium/analytic chemistry derive species).
  kElementsOnly = 0,
  // Non-equilibrium light network (see chemical_network.hpp manifest).
  kHybridLight = 1,
  // Hybrid + reduced cold C/O network.
  kDetailedCarbonOxygen = 2,
};

enum class CompositionDustKind : std::uint8_t {
  kNone = 0,
  // Two-composition live dust: silicate (MgFeSiO4) and carbonaceous (C).
  kTwoReservoirV1 = 1,
};

struct CompositionLayoutRequest {
  CompositionNetworkKind network = CompositionNetworkKind::kElementsOnly;
  CompositionDustKind dust = CompositionDustKind::kNone;
  bool sgs_energy = false;
};

struct CompositionLane {
  std::string name;
  std::string unit;  // "mass_code" or "energy_code"
};

class CompositionLayout {
 public:
  static constexpr std::uint32_t k_layout_version = 1U;

  [[nodiscard]] static CompositionLayout make(const CompositionLayoutRequest& request);

  [[nodiscard]] const core::GasCompositionLayoutDescriptor& descriptor() const noexcept {
    return m_descriptor;
  }
  [[nodiscard]] const std::vector<CompositionLane>& lanes() const noexcept { return m_lanes; }
  [[nodiscard]] const CompositionLayoutRequest& request() const noexcept { return m_request; }
  [[nodiscard]] std::uint32_t elementLane(ChemicalElement element) const noexcept {
    return m_descriptor.element_lane_begin + static_cast<std::uint32_t>(element);
  }
  // Species lane index for network species slot (see chemical_network.hpp).
  [[nodiscard]] std::uint32_t speciesLane(std::size_t species_slot) const;
  [[nodiscard]] std::uint32_t dustLane(std::size_t dust_slot) const;
  [[nodiscard]] bool hasSgs() const noexcept {
    return m_descriptor.sgs_energy_lane != core::GasCompositionLayoutDescriptor::k_no_lane;
  }

 private:
  CompositionLayoutRequest m_request{};
  core::GasCompositionLayoutDescriptor m_descriptor{};
  std::vector<CompositionLane> m_lanes;
};

// Dust reservoir stoichiometry (mass fraction of each element in the solid).
// Silicate: olivine MgFeSiO4 (Draine 2003 ARA&A 41, 241, "astrosilicate"
// composition approximated by MgFeSiO4); carbonaceous: pure C (graphite/
// amorphous carbon). These are compositional definitions, not fitted rates.
inline constexpr std::size_t k_dust_reservoir_count = 2U;
enum class DustReservoir : std::uint8_t { kSilicate = 0, kCarbonaceous = 1 };
[[nodiscard]] const std::array<double, k_chemical_element_count>& dustElementMassFractions(
    DustReservoir reservoir);

// ---------------------------------------------------------------------------
// Composition initialization policies (explicit; never silent).
// ---------------------------------------------------------------------------
struct PrimordialCompositionPolicy {
  // Primordial helium mass fraction Y_p (Aver, Olive & Skillman 2015, JCAP
  // 07, 011: Y_p = 0.2449 +- 0.0040).
  double helium_mass_fraction = 0.2449;
};

// Element masses for a cell of gas mass m with total metal mass m_Z: metals
// follow the A09 solar pattern; H/He split the remainder with Y_p scaled so
// that X+Y+Z=1 (Y = Y_p + (Y_sun - Y_p) Z/Z_sun, linear helium enrichment,
// capped so H stays positive). Documented initialization, not a measurement.
[[nodiscard]] std::array<double, k_chemical_element_count> initialElementMasses(
    double gas_mass_code, double metal_mass_code, const PrimordialCompositionPolicy& policy);

// ---------------------------------------------------------------------------
// Conservative AMR transfer kernels (interfaces for Part 3 regrid wiring).
// prolong: children receive the parent's mass fractions times their own
//   mass share (exact per-lane conservation: sum children = parent).
// restrict: parent lane = sum of child lanes.
// ---------------------------------------------------------------------------
void prolongComposition(std::span<const double> parent_lanes,
                        std::span<const double> child_mass_share,
                        std::span<double> child_lanes_row_major);
void restrictComposition(std::span<const double> child_lanes_row_major,
                         std::size_t lane_count,
                         std::span<double> parent_lanes);

// Recompute GasCellSidecar::metal_mass_code from element (+dust) lanes for
// the listed rows (all rows when empty). Returns max relative drift seen
// between the old mirror and the recomputed value (diagnostic).
double synchronizeMetalMirrorFromComposition(core::SimulationState& state,
                                             const CompositionLayout& layout,
                                             std::span<const std::uint32_t> rows);
// Same, using the state's own descriptor (element lanes in catalog order).
double synchronizeMetalMirrorFromComposition(core::SimulationState& state,
                                             std::span<const std::uint32_t> rows);

// Enable the composition sidecar on a state that has gas rows, initializing
// from the existing total-metal lane with the explicit policy. Throws if the
// sidecar is already enabled with a different layout.
void enableGasComposition(core::SimulationState& state, const CompositionLayout& layout,
                          const PrimordialCompositionPolicy& policy);

// Total extensive mass held by a row's composition lanes that are part of the
// gas mass (elements + dust); must equal CellSoa::mass_code to tolerance.
[[nodiscard]] double compositionGasMass(const core::GasCompositionSidecar& sidecar,
                                        std::size_t row);

// Uniform removal of a fraction of a cell's gas (star formation / accretion):
// scales every extensive lane; returns removed element masses (+ dust folded
// into elements by stoichiometry) for the receiving reservoir.
std::array<double, k_chemical_element_count> removeGasFractionFromComposition(
    core::GasCompositionSidecar& sidecar, const CompositionLayout& layout, std::size_t row,
    double removed_fraction);

}  // namespace cosmosim::physics
