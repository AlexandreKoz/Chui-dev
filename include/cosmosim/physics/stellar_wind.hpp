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

#include <cstdint>
#include <string>
#include <vector>

namespace cosmosim::physics {

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

 private:
  [[nodiscard]] double cumulative(std::size_t quantity, double age_years, double metallicity, bool& clamped) const;
  std::string m_source;
  std::string m_limitations;
  std::vector<double> m_z;
  std::vector<double> m_log_age;
  // [age][z][quantity] flattened, quantity 0 = energy, 1 = momentum.
  std::vector<double> m_values;
};

}  // namespace cosmosim::physics
