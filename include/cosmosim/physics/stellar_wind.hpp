#pragma once

// CHUI Part 3 Checkpoint C — time-dependent stellar-wind / pre-SN mechanical
// feedback channel (F02).
//
// Channel separation (no double counting by construction):
//  * MASS and ELEMENT return of massive-star winds is already part of the
//    SSP core-collapse continuous return: yield tables give the total ejecta
//    M_ini - M_final of a star, which includes its wind mass loss. This
//    channel therefore returns NO mass, metals or dust.
//  * This channel adds only the mechanical ENERGY and radial MOMENTUM of
//    winds (and, if the supplied table includes them, other pre-SN mechanical
//    sources) as functions of SSP age and metallicity, per unit birth mass.
//  * CC and Ia supernova energy/momentum stay in the discrete event channel;
//    AGB return stays in the AGB continuous channel. The radiation-pressure
//    prescription `fb_presn_wind_momentum_lsun_per_c` (EXPERIMENTAL user
//    closure) is mutually exclusive with this channel (config refusal).
//
// Data: `chui_ssp_wind_table_v1` text (owner-converted from a licensed SSP
// synthesis output, e.g. Starburst99 mechanical luminosity / momentum tables;
// nothing is bundled). Format:
//
//   # chui_ssp_wind_table_v1
//   # source = <free text provenance; required>
//   metallicities <N> Z_1 ... Z_N              (strictly increasing, >= 0)
//   ages <M>                                    (number of age rows, >= 2)
//   <log10 age/yr> <E_1> <p_1> ... <E_N> <p_N>   (M rows, strictly increasing age)
//
// E = cumulative mechanical energy [erg / Msun of birth mass], p = cumulative
// radial momentum [Msun km/s / Msun of birth mass], both non-decreasing in
// age (refused otherwise: a cumulative budget can never decrease). Interval
// budgets are differences of cumulative values, so they are exactly additive
// under step subdivision and restart. Before the first tabulated age the
// cumulative value is linear in age from zero; after the last it is constant.
// A metallicity outside the table is clamped and COUNTED (reported), never
// silently extrapolated.
//
// Patch 58 (F02): versioned, machine-readable compatibility metadata.
// `chui_ssp_wind_table_v2` has the v1 data layout plus REQUIRED header
// declarations, each exactly once (a duplicate or a malformed value refuses
// the table):
//
//   # chui_ssp_wind_table_v2
//   # source = <provenance; required>
//   # license = <licence statement; required>
//   # imf = <identity token [a-z0-9_]+, e.g. chabrier2003_individual,
//           kroupa2001; "unknown" allowed and never compatible>
//   # imf_mass_min_msun = <finite > 0>
//   # imf_mass_max_msun = <finite > imf_mass_min_msun>
//   # normalization = per_msun_initial_stellar_mass
//   # energy_unit = erg_per_msun
//   # momentum_unit = msun_km_s_per_msun
//   # metallicity_definition = metal_mass_fraction
//
// The metallicity grid remains the structural `metallicities` line. A v1
// table declares none of this: its IMF is UNDECLARED, and its free-text
// `# limitation =` lines are surfaced as text only — never parsed as
// scientific metadata. A v1 table carrying any v2-reserved key is refused
// (ambiguous schema). SHA-256 binding is unchanged.
//
// Compatibility with the CHUI SSP (assessStellarWindCompatibility):
//   strict (default): v2, imf == k_chui_ssp_imf_identity, IMF mass limits
//     equal to the configured SSP limits (relative 1e-9), more than one
//     metallicity; out-of-grid metallicity at use refuses the stage.
//   approximate_research: any declared or undeclared IMF / single-Z table is
//     admitted as a NAMED approximation (recorded; no IMF reweighting, no
//     multi-metallicity synthesis); metallicity clamps counted as before.

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace cosmosim::physics {

enum class StellarWindTableSchema : std::uint8_t { kV1 = 1, kV2 = 2 };

struct StellarWindTableMetadata {
  StellarWindTableSchema schema = StellarWindTableSchema::kV1;
  std::string license;                 // v2 required; v1 optional (informational)
  std::string imf_identity;            // v2 only; empty = undeclared (v1)
  double imf_mass_min_msun = 0.0;      // v2 only
  double imf_mass_max_msun = 0.0;      // v2 only
  std::string normalization;           // v2 only
  std::string energy_unit;             // v2 only
  std::string momentum_unit;           // v2 only
  std::string metallicity_definition;  // v2 only
};

enum class StellarWindCompatibilityPolicy : std::uint8_t { kStrict = 0, kApproximateResearch = 1 };

struct StellarWindCompatibilityReport {
  bool imf_declared = false;
  bool imf_identity_matches = false;
  bool imf_mass_limits_match = false;
  bool single_metallicity = false;
  bool admitted = false;
  // Non-empty when the strict policy refuses, or the named approximations
  // the research policy accepted ("none" when exact-compatible).
  std::string refusal;
  std::string approximations;
};

[[nodiscard]] const char* stellarWindCompatibilityPolicyName(StellarWindCompatibilityPolicy policy) noexcept;

struct StellarWindIntervalBudget {
  double energy_erg_per_msun = 0.0;
  double momentum_msun_kms_per_msun = 0.0;
  bool metallicity_clamped = false;
};

class StellarWindTable {
 public:
  [[nodiscard]] static StellarWindTable loadFromFile(const std::string& path, const std::string& expected_sha256);
  [[nodiscard]] static StellarWindTable loadFromString(const std::string& text, const std::string& label);

  [[nodiscard]] StellarWindIntervalBudget interval(double age_begin_years, double age_end_years,
                                                   double metallicity) const;
  [[nodiscard]] const std::string& source() const noexcept { return m_source; }
  // Patch 53: declared model limitations ('# limitation =' header lines, e.g.
  // the SB99 Kroupa IMF vs the CHUI Chabrier SSP IMF, single metallicity),
  // joined with '; ' and surfaced in the run's data-assumption record.
  [[nodiscard]] const std::string& limitations() const noexcept { return m_limitations; }
  [[nodiscard]] std::size_t metallicityCount() const noexcept { return m_z.size(); }
  [[nodiscard]] std::size_t ageCount() const noexcept { return m_log_age.size(); }
  [[nodiscard]] double metallicityMin() const noexcept { return m_z.empty() ? 0.0 : m_z.front(); }
  [[nodiscard]] double metallicityMax() const noexcept { return m_z.empty() ? 0.0 : m_z.back(); }
  [[nodiscard]] const StellarWindTableMetadata& metadata() const noexcept { return m_metadata; }
  [[nodiscard]] const std::string& sha256() const noexcept { return m_sha256; }

 private:
  [[nodiscard]] double cumulative(std::size_t quantity, double age_years, double metallicity, bool& clamped) const;
  std::string m_source;
  std::string m_limitations;
  std::string m_sha256;
  StellarWindTableMetadata m_metadata;
  std::vector<double> m_z;
  std::vector<double> m_log_age;
  // [age][z][quantity] flattened, quantity 0 = energy, 1 = momentum.
  std::vector<double> m_values;
};

// Compares the table's declared assumptions with the CHUI SSP actually
// selected (IMF identity token + configured IMF mass limits) under `policy`.
[[nodiscard]] StellarWindCompatibilityReport assessStellarWindCompatibility(
    const StellarWindTable& table, std::string_view ssp_imf_identity, double ssp_imf_mass_min_msun,
    double ssp_imf_mass_max_msun, StellarWindCompatibilityPolicy policy);

}  // namespace cosmosim::physics
