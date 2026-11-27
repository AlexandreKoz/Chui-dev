#pragma once

// CHUI-EXT-PHYS-REF-v1 S03/I03: read-only science accessor shared by every
// observable. One canonical authority per quantity: the SAME hydrogen,
// electron, temperature, velocity and magnetic samples feed Ly-a, HI, SZ, DM,
// free-free, synchrotron and Faraday products. Values are converted on the
// fly to proper SI (no copy of the gas state; nothing is ever written back).
//
// Sources (all from core::SimulationState, rows = CellSoa rows):
//   n_X     composition species lanes (extensive code mass) / (proper cell
//           volume * particle mass); species slot order is the persistent
//           physics::ChemicalSpecies order (HI, HII, H-, H2, H2+, HeI, HeII,
//           HeIII, ...). Requires a composition layout with >= 8 species
//           lanes; otherwise every species density is UNAVAILABLE.
//   n_e     charge neutrality over the stored ionic species (H+, H2+, He+,
//           2 He++, minus H-); metals' electrons are not included (stated).
//   T       GasCellSidecar::temperature_code (kelvin cache of the thermal
//           refresh).
//   v       GasCellSidecar::velocity_*_peculiar (code velocity -> m/s, peculiar).
//   B       supplied by the caller as cell-centred SI samples of the CT
//           authority (MHD runs); absent -> UNAVAILABLE (polarized and
//           synchrotron products refuse).
//   V_cell  from the AMR patch geometry of the row (comoving code length^3),
//           proper volume a^3 V; legacy undeclared patch geometry refuses.

#include <array>
#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/spectral/provider_status.hpp"

namespace cosmosim::observables {

class GasScienceAccessor {
 public:
  GasScienceAccessor(const core::SimulationState& state, const core::UnitSystem& units, double scale_factor);

  [[nodiscard]] std::size_t rows() const noexcept { return m_rows; }
  [[nodiscard]] double scaleFactor() const noexcept { return m_a; }
  [[nodiscard]] bool hasSpecies() const noexcept { return m_species_begin != k_no_species; }

  [[nodiscard]] spectral::ProviderValue neutralHydrogenDensity(std::size_t row) const;  // m^-3
  [[nodiscard]] spectral::ProviderValue protonDensity(std::size_t row) const;           // n_HII
  [[nodiscard]] spectral::ProviderValue electronDensity(std::size_t row) const;
  [[nodiscard]] spectral::ProviderValue ionZ2Density(std::size_t row) const;  // sum Z^2 n_i (free-free)
  [[nodiscard]] spectral::ProviderValue temperature(std::size_t row) const;   // K
  [[nodiscard]] std::array<double, 3> peculiarVelocity(std::size_t row) const;  // m/s
  [[nodiscard]] double properVolume(std::size_t row) const;                     // m^3

  // Optional cell-centred magnetic field (Tesla), sized rows() or empty.
  void attachMagneticField(std::span<const double> bx, std::span<const double> by, std::span<const double> bz);
  [[nodiscard]] bool hasMagneticField() const noexcept { return !m_bx.empty(); }
  [[nodiscard]] std::array<double, 3> magneticField(std::size_t row) const;  // throws when unavailable

  [[nodiscard]] std::string provenance() const;

 private:
  static constexpr std::size_t k_no_species = static_cast<std::size_t>(-1);
  [[nodiscard]] spectral::ProviderValue speciesDensity(std::size_t row, std::size_t slot) const;

  const core::SimulationState& m_state;
  core::UnitSystem m_units;
  double m_a = 1.0;
  std::size_t m_rows = 0;
  std::size_t m_species_begin = k_no_species;
  std::vector<double> m_row_volume_comoving_code;  // per row, from patch geometry
  std::span<const double> m_bx, m_by, m_bz;
};

}  // namespace cosmosim::observables
