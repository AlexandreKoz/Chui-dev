#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// Uniform ionizing background (CHUI Part 2).
//
// Production source: a user-supplied table in the TREECOOL format used by
// GADGET/AREPO (columns: log10(1+z), Gamma_HI, Gamma_HeI, Gamma_HeII [s^-1],
// q_HI, q_HeI, q_HeII [erg s^-1 per atom]). Openly published examples are the
// Faucher-Giguere (2020, MNRAS 493, 1614) and Haardt & Madau (2012, ApJ 746,
// 125) TREECOOL files; they are NOT bundled here (licensing/redistribution is
// the owner's decision). Missing file => SourceModelUnavailableError for the
// production model, with acquisition instructions. Above the table's maximum
// redshift the background is zero (pre-reionization), below its minimum the
// last row is held (documented clamp).
//
// The pre-Part-2 kHm12 / kFg20 enumerators implement ANALYTIC PROXIES of an
// overall heating scale only; they are labelled legacy_analytic_proxy and are
// refused by every Part 2 thermochemistry backend.
// ---------------------------------------------------------------------------
struct UvBackgroundRates {
  double gamma_hi_s = 0.0;
  double gamma_hei_s = 0.0;
  double gamma_heii_s = 0.0;
  double heat_hi_erg_s = 0.0;
  double heat_hei_erg_s = 0.0;
  double heat_heii_erg_s = 0.0;
};

struct UvBackgroundTableProvenance {
  std::string source_path;
  std::string sha256_hex;
  std::size_t row_count = 0;
  double redshift_min = 0.0;
  double redshift_max = 0.0;
};

class UvBackgroundTable {
 public:
  // Throws SourceModelUnavailableError on missing/malformed files or when an
  // expected SHA-256 is given and does not match.
  [[nodiscard]] static UvBackgroundTable loadTreecool(const std::string& path,
                                                      const std::string& expected_sha256_hex = {});
  // Construct from explicit rows (tests / programmatic use).
  [[nodiscard]] static UvBackgroundTable fromRows(std::vector<double> log10_one_plus_z,
                                                  std::vector<UvBackgroundRates> rows);

  [[nodiscard]] UvBackgroundRates ratesAtRedshift(double redshift) const;
  [[nodiscard]] const UvBackgroundTableProvenance& provenance() const noexcept { return m_provenance; }
  [[nodiscard]] bool empty() const noexcept { return m_log10_one_plus_z.empty(); }

 private:
  UvBackgroundTableProvenance m_provenance;
  std::vector<double> m_log10_one_plus_z;
  std::vector<UvBackgroundRates> m_rows;
};

// Rahmati et al. (2013, MNRAS 430, 2427) eq. 14 photoionization suppression
// fitted to RT simulations:
//  Gamma/Gamma_UVB = 0.98 [1 + (n_H/n_0)^1.64]^-2.28 + 0.02 [1 + n_H/n_0]^-0.84
//  n_0 = 6.73e-3 (sigma/2.49e-18)^(-2/3) (T/1e4)^0.17 (Gamma_UVB/1e-12)^(2/3)
// with f_g = 0.17 and the z-independent effective cross-section default; the
// same factor is applied to HeI/HeII rates and photoheating (approximation,
// documented). This is a calibrated LOCAL approximation, not radiative transfer.
[[nodiscard]] double rahmati2013ShieldingFactor(double hydrogen_number_density_cm3, double temperature_k,
                                                double gamma_hi_uvb_s,
                                                double effective_cross_section_cm2 = 2.49e-18);

// Shared hex-encoded SHA-256 of a file's bytes (FIPS 180-4 implementation);
// used for every external data file manifest.
[[nodiscard]] std::string sha256FileHex(const std::string& path);

}  // namespace cosmosim::physics
