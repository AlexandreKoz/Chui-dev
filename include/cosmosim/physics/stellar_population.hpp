#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

// ===========================================================================
// CHUI Part 2 single stellar population (SSP) model.
//
// IMF: Chabrier (2003, PASP 115, 763) individual-star form, Table 1:
//   xi(log m) = 0.158 exp[-(log m - log 0.079)^2 / (2 * 0.69^2)],  m <= 1
//   xi(log m) = 0.0443 m^-1.3,                                      m > 1
// normalized to unit mass over [m_min, m_max] (Msun).
// Lifetimes: Raiteri, Villata & Navarro (1996, A&A 315, 105) eq. 3,
//   log t[yr] = a0(Z) + a1(Z) log M + a2(Z) (log M)^2, Z clamped to
//   [7e-5, 0.03] (fit range).
// Remnants: Kalirai et al. (2008, ApJ 676, 594) IFMR m_f = 0.109 m + 0.394
//   below m_cc,min; 1.4 Msun neutron stars in [m_cc,min, m_cc,max]; above
//   m_cc,max direct collapse with no ejecta (configurable mass limit).
// Type Ia: power-law delay-time distribution R(t) proportional to t^-1 on
//   [t_min, t_max] normalized to N_Ia per formed Msun (Maoz & Mannucci 2012,
//   PASA 29, 447: N_Ia ~ 1.3e-3 Msun^-1), Chandrasekhar ejecta 1.38 Msun
//   taken from the remnant reservoir.
// Element yields: REFERENCE mode returns the BIRTH element mixture with an
//   optional configured net metal yield in the solar pattern (labelled
//   REFERENCE ONLY, not astrophysical yields). PRODUCTION mode requires a
//   user-supplied chui_ssp_yield_table_v1 (net yields per element vs
//   progenitor mass and metallicity for AGB/CC, per-event ejecta for Ia) with
//   mandatory citation/license; nothing is bundled or synthesized.
// ===========================================================================

enum class SspYieldMode : std::uint8_t { kReference = 0, kTable = 1 };

struct SspConfig {
  double imf_mass_min_msun = 0.1;
  double imf_mass_max_msun = 100.0;
  double cc_mass_min_msun = 8.0;
  double cc_mass_max_msun = 40.0;
  double snia_per_msun = 1.3e-3;
  double snia_delay_min_years = 4.0e7;
  double snia_delay_max_years = 1.37e10;
  double energy_per_sn_erg = 1.0e51;
  double reference_net_metal_yield = 0.0;
  bool stochastic_events = true;
  SspYieldMode yield_mode = SspYieldMode::kReference;
};

// Yield table, chui_ssp_yield_table_v1 (see docs/stellar_evolution.md).
class SspYieldTable {
 public:
  [[nodiscard]] static SspYieldTable loadFromFile(const std::string& path, const std::string& expected_sha256 = {});
  [[nodiscard]] static SspYieldTable loadFromString(const std::string& text, const std::string& label);
  // Net yield (Msun per star of progenitor mass m) for channel 0=AGB, 1=CC,
  // and ejecta mass; linear in mass and metallicity, clamped to the grid
  // (clamping counted by the caller via the returned flag).
  struct Row {
    double ejecta_msun = 0.0;
    std::array<double, k_chemical_element_count> net_yield_msun{};
  };
  [[nodiscard]] Row massive(std::size_t channel, double mass_msun, double metallicity, bool& clamped) const;
  [[nodiscard]] Row snia(double metallicity, bool& clamped) const;
  [[nodiscard]] const std::string& citation() const noexcept { return m_citation; }
  [[nodiscard]] const std::string& license() const noexcept { return m_license; }
  [[nodiscard]] const std::string& sha256() const noexcept { return m_sha256; }
  std::uint64_t balanced_rows = 0;  // rows whose net-yield imbalance was moved to H

 private:
  struct Channel {
    std::vector<double> metallicity;
    std::vector<double> mass;  // empty for Ia
    std::vector<Row> rows;     // [z][m] or [z]
  };
  std::array<Channel, 3> m_channels{};  // agb, cc, ia
  std::string m_citation, m_license, m_sha256;
};

class SspModel {
 public:
  SspModel(SspConfig config, const SspYieldTable* table = nullptr);

  [[nodiscard]] const SspConfig& config() const noexcept { return m_config; }
  [[nodiscard]] double lifetimeYears(double mass_msun, double metallicity) const;
  [[nodiscard]] double turnoffMassMsun(double age_years, double metallicity) const;
  [[nodiscard]] double remnantMassMsun(double mass_msun) const;
  // Per unit formed stellar mass [Msun^-1 or Msun/Msun]:
  [[nodiscard]] double numberBetween(double m_lo, double m_hi) const;
  [[nodiscard]] double massBetween(double m_lo, double m_hi) const;
  [[nodiscard]] double ccCumulativeNumber(double age_years, double metallicity) const;
  [[nodiscard]] double iaCumulativeNumber(double age_years) const;

  struct ChannelReturn {
    double mass_fraction = 0.0;  // returned mass / formed mass
    std::array<double, k_chemical_element_count> element_fraction{};  // per formed mass
  };
  // Continuous return of AGB (channel 0) and CC (channel 1) ejecta for ages
  // [age_begin, age_end] given the birth element mass fractions; Ia handled
  // per event.
  [[nodiscard]] ChannelReturn massiveReturn(std::size_t channel, double age_begin, double age_end,
                                            const std::array<double, k_chemical_element_count>& birth_x) const;
  [[nodiscard]] ChannelReturn iaReturnPerEvent(const std::array<double, k_chemical_element_count>& birth_x) const;
  std::uint64_t mutable table_clamps = 0;

 private:
  [[nodiscard]] double phi(double m) const;  // dN/dm normalized
  [[nodiscard]] std::array<double, k_chemical_element_count> ejectaElements(
      std::size_t channel, double m, double metallicity,
      const std::array<double, k_chemical_element_count>& birth_x, double& ejecta) const;

  SspConfig m_config;
  const SspYieldTable* m_table = nullptr;
  double m_normalization = 1.0;
};

// ---------------------------------------------------------------------------
// Persistent per-star population record (particle-indexed module sidecar
// "chui_stellar_population_v1", migrates with the star and is checkpointed).
// Event chronology: discrete SNe of channel c are the jumps of a Poisson
// process in the star's cumulative expected count Lambda_c(age). The record
// stores only the NEXT threshold Lambda_next (O(1) memory per star, no
// materialized lifetime queue) and the realized ordinal. Thresholds are
// Lambda_{k+1} = Lambda_k + E_k with E_k = -ln(1 - u_k), u_k a counter-based
// hash of (birth_key, channel, k): realizations are independent of step
// subdivision, row order, rank and restart (identity = birth_key, channel,
// ordinal). A rejected stage never advances the stored record.
// ---------------------------------------------------------------------------
struct StellarPopulationRecord {
  std::array<double, k_chemical_element_count> birth_element_mass_fraction{};
  double birth_time_code = 0.0;
  double cc_next_threshold = 0.0;
  double ia_next_threshold = 0.0;
  std::uint32_t cc_realized = 0;
  std::uint32_t ia_realized = 0;
  std::uint64_t event_stream_key = 0;  // birth_key or association key
  double returned_mass_code = 0.0;     // cumulative, all channels
  std::uint32_t initialized = 0;
  std::uint32_t reserved = 0;
};
inline constexpr const char* k_stellar_population_module = "chui_stellar_population_v1";
inline constexpr std::uint32_t k_stellar_population_schema = 1U;

[[nodiscard]] double eventExponentialDraw(std::uint64_t stream_key, std::uint32_t channel, std::uint32_t ordinal);
// Initialize thresholds for a new star.
[[nodiscard]] StellarPopulationRecord makeStellarPopulationRecord(
    const std::array<double, k_chemical_element_count>& birth_x, double birth_time_code, std::uint64_t stream_key);
// Realize events of one channel up to cumulative expected count `lambda_end`
// (record updated only through the returned copy).
struct RealizedEvents {
  std::uint32_t count = 0;
  std::uint32_t first_ordinal = 0;
};
[[nodiscard]] RealizedEvents realizeEvents(StellarPopulationRecord& record, std::uint32_t channel, double lambda_end);

class StellarPopulationStore {
 public:
  explicit StellarPopulationStore(const core::SimulationState& state);
  [[nodiscard]] bool has(std::uint64_t particle_id) const;
  [[nodiscard]] StellarPopulationRecord get(std::uint64_t particle_id) const;
  void put(std::uint64_t particle_id, const StellarPopulationRecord& record);
  // Write the block back (replaces the module sidecar). Rows of particles
  // that no longer exist are dropped; missing rows are refused.
  void store(core::SimulationState& state) const;
  [[nodiscard]] std::size_t size() const noexcept { return m_ids.size(); }

 private:
  std::vector<std::uint64_t> m_ids;  // sorted
  std::vector<StellarPopulationRecord> m_records;
};

}  // namespace cosmosim::physics
