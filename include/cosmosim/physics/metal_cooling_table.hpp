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
// the table's REFERENCE abundance (header solar_reference). Sign convention:
// POSITIVE = net cooling, negative = net (photo)heating (WSS09 Net_Cooling).
// The runtime applies Wiersma, Schaye & Smith (2009, MNRAS 393, 99) eq. 4:
//   Lambda_e = table_e(T, n_H, z) n_H^2 [(n_e/n_H)_gas / (n_e/n_H)_ref]
//              x 10^[e/H],   10^[e/H] = (n_e/n_H) / (n_e/n_H)_ref
//            = (X_e / X_H) / (X_e,ref / X_H,ref)
// with GAS-PHASE (post-depletion) mass fractions X and the table's
// reference_hydrogen_mass_fraction X_H,ref (Part 2.5 correction: the former
// scaling X_e / X_e,ref omitted the hydrogen normalization, a 7.6% error for
// X_H = 0.76 against the CLOUDY solar X_H = 0.7065). The electron correction
// uses the table's reference-abundance n_e/n_H block (WSS09 Solar/
// Electron_density_over_n_h) when present. H and He are NEVER taken from this
// table (the non-equilibrium primordial network owns them, including Compton
// cooling): a header declaring H or He columns is refused.
//
// UV background: tables computed under photoionization embed a specific UVB
// (WSS09: Haardt & Madau 2001 + CMB, CLOUDY 07.02). The header MUST declare
// `uv_background` (`none` for a collisional/CIE table). Optionally
// `uv_background_treecool_sha256` names the TREECOOL file the table was
// computed with. The runtime compatibility policy is
// checkUvBackgroundCompatibility() (see below).
//
// Format (whitespace separated, '#' comments ignored except header keys):
//   # chui_metal_cooling_table_v1
//   # citation = <text>          (required)
//   # license = <text>           (required)
//   # solar_reference = <text>   (required)
//   # uv_background = <none | text identity>        (required, Part 2.5)
//   # uv_background_treecool_sha256 = <hex>         (optional)
//   elements C N O Ne Mg Si S Ca Fe
//   log10_T <N_T values, strictly increasing>
//   log10_nH <N_n values, strictly increasing>
//   redshift <N_z values, strictly increasing>        (>= 1 value)
//   reference_mass_fraction <one value per element>
//   reference_hydrogen_mass_fraction <one value>       (required, Part 2.5)
//   ne_over_nH <N_z * N_T * N_n values>                (optional block)
//   lambda <element> <N_z * N_T * N_n values>          (one block per element)
// Index order: z slowest, then T, then n_H fastest.
//
// Interpolation: trilinear in (z, log T, log n_H) of the transformed value
// asinh(Lambda / 1e-30): smooth, sign-safe (net metal heating can be
// negative under photoionization) and continuous across grid cells — never a
// blind log of a signed quantity.
// Out of range: T outside the table => always throws; kReject also throws
// for n_H or z outside the table; kClampDensityAndRedshift clamps n_H and z
// to the table edges and sets Result::clamped (counted by the thermo-
// chemistry ledger as metal_table_clamps). A SINGLETON axis covers exactly
// its one coordinate (match tolerance 1e-8 x max(1, |value|)): e.g. a
// one-redshift table is refused at any other redshift under kReject, and is
// clamped AND reported under kClampDensityAndRedshift — never silently
// treated as valid at all redshifts (Part 2.5 correction).
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
  std::string uv_background;                 // "none" or the declared UVB identity
  std::string uv_background_treecool_sha256;  // optional
  std::size_t grid_points = 0;
};

// UV-background compatibility between a metal table and the configured CHUI
// background (Part 2.5 policy). Compatible iff
//   table uv_background == "none"  and CHUI has no UVB, or
//   table uv_background_treecool_sha256 equals the SHA-256 of the configured
//   TREECOOL file.
// Any other combination (e.g. a WSS09/HM01 table with a FG20 TREECOOL, or a
// photoionized table with no CHUI UVB) is a MISMATCH: the metal photoheating
// and ionization balance embedded in the table do not match the H/He
// photo-rates. A mismatch is refused unless accept_mismatch is true
// (physics.metal_cooling_accept_uv_background_mismatch), in which case the
// returned diagnostic is logged and recorded with the run provenance.
struct UvBackgroundCompatibility {
  bool compatible = false;
  bool accepted_mismatch = false;
  std::string diagnostic;
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
  // chui_uvb_sha256: empty when CHUI runs without a tabulated UVB.
  // Throws SourceModelUnavailableError on a refused mismatch.
  [[nodiscard]] UvBackgroundCompatibility checkUvBackgroundCompatibility(const std::string& chui_uvb_sha256,
                                                                         bool accept_mismatch) const;
  [[nodiscard]] double referenceHydrogenMassFraction() const noexcept { return m_reference_hydrogen_mass_fraction; }
  [[nodiscard]] const MetalCoolingTableProvenance& provenance() const noexcept { return m_provenance; }
  [[nodiscard]] double temperatureMinK() const noexcept;
  [[nodiscard]] double temperatureMaxK() const noexcept;

 private:
  [[nodiscard]] double interpolate(const std::vector<double>& block, double z, double lt, double ln,
                                   bool& clamped, MetalTableRangePolicy policy) const;

  MetalCoolingTableProvenance m_provenance;
  std::vector<ChemicalElement> m_elements;
  std::vector<double> m_reference_mass_fraction;
  double m_reference_hydrogen_mass_fraction = 0.0;
  std::vector<double> m_log_t;
  std::vector<double> m_log_n;
  std::vector<double> m_redshift;
  std::vector<double> m_ne_over_nh;  // optional, stored as asinh(n_e/n_H)
  std::vector<std::vector<double>> m_lambda_transformed;  // per element
};

}  // namespace cosmosim::physics
