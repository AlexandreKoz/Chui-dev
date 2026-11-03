#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// Element-resolved metal cooling table (CHUI ASCII format v1).
//
// Intended content: per-element net metal-line cooling of photo-/collisionally
// ionized gas computed with an equilibrium-ionization code (e.g. the Wiersma,
// Schaye & Smith 2009 CLOUDY tables, or Ploeckinger & Schaye 2020). The data
// are NOT bundled; the owner converts a licensed source table into this
// format and records its citation/license/checksum in the header. A missing
// file stops the metal-table model selection with SourceModelUnavailableError.
//
// Normalization: column = Lambda_e / n_H^2 [erg cm^3 s^-1] for element e at
// the table's REFERENCE abundance (header solar_reference). The runtime
// scales by the GAS-PHASE (post-depletion) element abundance relative to that
// reference: Lambda_e = table_e(T, n_H, z) n_H^2 (n_e/n_H,gas / n_e/n_H,eq)
// (x_e,gas / x_e,ref), where the optional electron correction uses the
// table's equilibrium n_e/n_H column (Wiersma et al. 2009 eq. 4). H and He
// are NEVER taken from this table (the non-equilibrium primordial network owns
// them): a header declaring H or He columns is refused.
//
// Format (whitespace separated, '#' comments ignored except header keys):
//   # chui_metal_cooling_table_v1
//   # citation = <text>          (required)
//   # license = <text>           (required)
//   # solar_reference = <text>   (required)
//   elements C N O Ne Mg Si S Ca Fe
//   log10_T <N_T values, strictly increasing>
//   log10_nH <N_n values, strictly increasing>
//   redshift <N_z values, strictly increasing>        (>= 1 value)
//   reference_mass_fraction <one value per element>
//   ne_over_nH <N_z * N_T * N_n values>                (optional block)
//   lambda <element> <N_z * N_T * N_n values>          (one block per element)
// Index order: z slowest, then T, then n_H fastest.
//
// Interpolation: trilinear in (z, log T, log n_H) of the transformed value
// asinh(Lambda / 1e-30): smooth, sign-safe (net metal heating can be
// negative under photoionization) and continuous across grid cells — never a
// blind log of a signed quantity.
// Out of range: T outside the table => kReject throws; kClampDensityOnly
// clamps n_H and z to the table edges (counted) but still rejects T.
// ---------------------------------------------------------------------------
enum class MetalTableRangePolicy : std::uint8_t {
  kReject = 0,
  kClampDensityAndRedshift = 1,
};

struct MetalCoolingTableProvenance {
  std::string source_path;
  std::string sha256_hex;
  std::string citation;
  std::string license;
  std::string solar_reference;
  std::size_t grid_points = 0;
};

class MetalCoolingTable {
 public:
  [[nodiscard]] static MetalCoolingTable loadFromFile(const std::string& path,
                                                      const std::string& expected_sha256_hex = {});
  [[nodiscard]] static MetalCoolingTable loadFromString(const std::string& text, const std::string& label);

  struct Query {
    double temperature_k = 1.0e4;
    double hydrogen_number_density_cm3 = 1.0;
    double redshift = 0.0;
    double electron_per_hydrogen = -1.0;  // < 0: no electron correction
    // Gas-phase element mass fractions (post depletion), catalog order.
    std::array<double, k_chemical_element_count> gas_mass_fraction{};
    double hydrogen_mass_fraction = 0.75;
  };
  struct Result {
    double net_cooling_erg_cm3_s = 0.0;  // volumetric, signed
    bool clamped = false;
  };

  [[nodiscard]] Result evaluate(const Query& query, MetalTableRangePolicy policy) const;
  [[nodiscard]] const MetalCoolingTableProvenance& provenance() const noexcept { return m_provenance; }
  [[nodiscard]] double temperatureMinK() const noexcept;
  [[nodiscard]] double temperatureMaxK() const noexcept;

 private:
  [[nodiscard]] double interpolate(const std::vector<double>& block, double z, double lt, double ln,
                                   bool& clamped, MetalTableRangePolicy policy) const;

  MetalCoolingTableProvenance m_provenance;
  std::vector<ChemicalElement> m_elements;
  std::vector<double> m_reference_mass_fraction;
  std::vector<double> m_log_t;
  std::vector<double> m_log_n;
  std::vector<double> m_redshift;
  std::vector<double> m_ne_over_nh;  // optional, stored as asinh(n_e/n_H)
  std::vector<std::vector<double>> m_lambda_transformed;  // per element
};

}  // namespace cosmosim::physics
