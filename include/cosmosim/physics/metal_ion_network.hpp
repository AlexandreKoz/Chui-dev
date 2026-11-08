#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/physics/source_interval.hpp"

namespace cosmosim::physics {

// ---------------------------------------------------------------------------
// C02 fixed detailed metal-ion reference (CHUI): non-equilibrium ionization
// ladders of selected metals carried by the fixed_detailed backend.
//
// Species inventory (explicit, physically fixed by atomic number Z):
//   C  : C I  .. C VII   (7 stages, charge 0..6)
//   O  : O I  .. O IX    (9 stages, charge 0..8)
//   Si : Si I .. Si XV   (15 stages, charge 0..14)
//   Fe : Fe I .. Fe XXVII(27 stages, charge 0..26)
// Any subset may be selected (physics.chemistry_metal_ion_elements).
//
// Equations (proper cgs, per element X with gas-phase nuclei n_X, stage
// fractions x_k, sum_k x_k = 1):
//   dx_k/dt = I_{k-1} x_{k-1} + R_{k+1} x_{k+1} - (I_k + R_k) x_k,
//   I_k = n_e C_k(T) + P_k(z)        (collisional + photo ionization k->k+1)
//   R_k = n_e alpha_k(T)             (radiative + dielectronic k->k-1)
// with C_Z = 0 and alpha_0 = 0. The matrix columns sum to zero, so the
// element inventory is a linear invariant; backward Euler preserves it and
// produces non-negative fractions (M-matrix). Free electrons contributed by
// the ladder are n_e,X = n_X sum_k k x_k; radiative losses are
// Lambda_X = n_e n_X sum_k x_k L_k(T) [erg cm^-3 s^-1]; metal photoheating
// (when the dataset carries it) is n_X sum_k x_k H_k(z).
//
// Rates are NEVER compiled into the code: every coefficient comes from an
// external dataset (schema chui_metal_ion_rates_v1, see
// docs/chemistry_rate_provenance.md) that names its provenance and licence
// and declares its own verification status. Charge exchange with H/He is not
// represented by schema v1 (named omitted channel); a dataset line with an
// unknown process is refused rather than silently ignored.
// ---------------------------------------------------------------------------

enum class MetalIonElement : std::uint8_t { kC = 0, kO = 1, kSi = 2, kFe = 3 };
inline constexpr std::size_t k_metal_ion_element_count = 4U;
inline constexpr std::array<std::uint8_t, k_metal_ion_element_count> k_metal_ion_atomic_number{6U, 8U, 14U, 26U};
// Sum over all selectable elements of (Z + 1) stages.
inline constexpr std::size_t k_metal_ion_stage_capacity = 7U + 9U + 15U + 27U;

[[nodiscard]] ChemicalElement chemicalElementOf(MetalIonElement element) noexcept;
[[nodiscard]] std::string_view metalIonElementSymbol(MetalIonElement element) noexcept;
[[nodiscard]] constexpr std::size_t metalIonStageCount(MetalIonElement element) noexcept {
  return static_cast<std::size_t>(k_metal_ion_atomic_number[static_cast<std::size_t>(element)]) + 1U;
}
// Spectroscopic stage label, e.g. (kO, 5) -> "O_VI".
[[nodiscard]] std::string metalIonStageName(MetalIonElement element, std::size_t stage);

// Ordered element subset; lane/stage offsets follow element order C, O, Si, Fe.
class MetalIonElementSet {
 public:
  constexpr MetalIonElementSet() = default;
  [[nodiscard]] static MetalIonElementSet fromMask(std::uint8_t mask);
  // "none" or a comma-separated list of C, O, Si, Fe (case-sensitive
  // symbols, no duplicates). Throws std::invalid_argument otherwise.
  [[nodiscard]] static MetalIonElementSet parse(std::string_view text);
  // Inverse of stageCount() for layout reconstruction (subset sums of
  // {7, 9, 15, 27} are unique). Throws when no subset matches.
  [[nodiscard]] static MetalIonElementSet fromStageCount(std::size_t stage_count);

  [[nodiscard]] std::uint8_t mask() const noexcept { return m_mask; }
  [[nodiscard]] bool empty() const noexcept { return m_mask == 0U; }
  [[nodiscard]] bool contains(MetalIonElement element) const noexcept {
    return (m_mask & (1U << static_cast<unsigned>(element))) != 0U;
  }
  [[nodiscard]] bool containsChemicalElement(ChemicalElement element) const noexcept;
  [[nodiscard]] std::size_t stageCount() const noexcept;
  // Offset of the element's first stage within the packed stage vector.
  [[nodiscard]] std::size_t stageOffset(MetalIonElement element) const;
  [[nodiscard]] std::string label() const;  // "none" or "C,O,Si,Fe"
  [[nodiscard]] bool operator==(const MetalIonElementSet& rhs) const noexcept { return m_mask == rhs.m_mask; }

 private:
  std::uint8_t m_mask = 0U;
};

enum class MetalIonRangePolicy : std::uint8_t {
  // Temperature outside the dataset grid rejects the source interval.
  kReject = 0,
  // Hold the boundary value; every clamp is counted and reported.
  kClampAndCount = 1,
};

struct MetalIonStageRates {
  double collisional_ionization_cm3_s = 0.0;  // k -> k+1, per free electron
  double recombination_cm3_s = 0.0;           // k -> k-1 (radiative + dielectronic)
  double cooling_erg_cm3_s = 0.0;             // L_k: loss per n_e n_{X^k}
  double photoionization_s = 0.0;             // k -> k+1 (dataset UVB)
  double photoheating_erg_s = 0.0;            // per X^k ion
};

// Immutable rate dataset. Format (text, '#' comments, whitespace separated):
//   schema chui_metal_ion_rates_v1
//   provenance <rest of line>               (required, non-empty)
//   license <rest of line>                  (required, non-empty)
//   verification verified_against_source | unverified
//   temperature_grid_log10_k <N> t_0 .. t_{N-1}    (strictly increasing, N>=2)
//   [redshift_grid <M> z_0 .. z_{M-1}]             (strictly increasing)
//   rate <El> <stage> collisional_ionization <N log10 values>
//   rate <El> <stage> recombination          <N log10 values>
//   rate <El> <stage> cooling                <N log10 values>
//   [photo <El> <stage> ionization <M log10 values>]
//   [photo <El> <stage> heating    <M log10 values>]
// "-inf" denotes an exactly zero coefficient. Every selected element needs
// all three temperature processes for every stage (ionization of the bare
// stage and recombination of the neutral stage must be -inf). Photo lines,
// when present for one stage, are required for all ionizable stages of that
// element and need the redshift grid.
class MetalIonRateTable {
 public:
  [[nodiscard]] static MetalIonRateTable loadFromFile(const std::string& path, MetalIonElementSet required);
  [[nodiscard]] static MetalIonRateTable parse(std::string_view text, std::string_view source_label,
                                               MetalIonElementSet required);

  [[nodiscard]] MetalIonElementSet elements() const noexcept { return m_elements; }
  [[nodiscard]] const std::string& provenance() const noexcept { return m_provenance; }
  [[nodiscard]] const std::string& license() const noexcept { return m_license; }
  [[nodiscard]] const std::string& sourceLabel() const noexcept { return m_source_label; }
  [[nodiscard]] bool declaredVerified() const noexcept { return m_declared_verified; }
  [[nodiscard]] bool hasPhotoRates(MetalIonElement element) const noexcept {
    return m_has_photo[static_cast<std::size_t>(element)];
  }
  // True when the dataset declares a supported energy convention and every
  // ionization potential (stages 0..Z-1) of `element`.
  [[nodiscard]] bool hasIonizationEnergyContract(MetalIonElement element) const noexcept;
  [[nodiscard]] const std::string& energyConvention() const noexcept { return m_energy_convention; }
  // chi_k (k -> k+1) in erg; throws std::out_of_range when absent.
  [[nodiscard]] double ionizationPotentialErg(MetalIonElement element, std::size_t stage) const;
  // Stored ionization potential per nucleus for a stage distribution,
  // sum_k x_k sum_{j<k} chi_j [erg]. Requires hasIonizationEnergyContract.
  [[nodiscard]] double storedIonizationEnergyPerNucleusErg(MetalIonElement element,
                                                           std::span<const double> fractions) const;
  [[nodiscard]] double minimumTemperatureK() const noexcept;
  [[nodiscard]] double maximumTemperatureK() const noexcept;
  // Rates of one stage at (T, z). `clamped` is set when T (or z) was outside
  // the grid under kClampAndCount; kReject throws std::range_error.
  [[nodiscard]] MetalIonStageRates evaluate(MetalIonElement element, std::size_t stage, double temperature_k,
                                            double redshift, MetalIonRangePolicy policy, bool& clamped) const;

 private:
  struct Series {
    std::vector<double> log10_values;  // -inf allowed (zero)
  };
  [[nodiscard]] std::size_t seriesIndex(MetalIonElement element, std::size_t stage, std::size_t process) const;
  [[nodiscard]] static double interpolate(const std::vector<double>& grid, const Series& series,
                                          double coordinate, MetalIonRangePolicy policy, bool& clamped,
                                          const char* what);

  MetalIonElementSet m_elements{};
  std::string m_provenance;
  std::string m_license;
  std::string m_source_label;
  bool m_declared_verified = false;
  std::vector<double> m_log_t;
  std::vector<double> m_redshift;
  std::array<bool, k_metal_ion_element_count> m_has_photo{};
  std::string m_energy_convention;
  // [element][stage] chi_k in erg; NaN = not supplied.
  std::vector<double> m_ionization_potential_erg;
  // [element][stage][process]; process 0..2 temperature, 3..4 photo.
  std::vector<Series> m_series;
};

// Local radiation context of the metal-ion ladder (P1-49-02). The single
// authority is ThermochemistrySolver::metalIonEnvironment(), used by the
// equilibrium initializer, both operator-split ion half-steps, and the
// ion cooling/photoheating terms inside the Newton energy solve.
enum class MetalIonPhotoAttenuationModel : std::uint8_t {
  // No attenuation (self_shielding off): dataset UVB rates as tabulated.
  kNone = 0,
  // Grey approximation: the local Rahmati (2013) HI self-shielding factor
  // Gamma_HI,shielded / Gamma_HI,UVB is applied to every metal-ion photo
  // rate and photoheating rate. This is a stated frequency-independent
  // approximation (NOT frequency-resolved transfer): ions with thresholds
  // above 13.6 eV are in reality shielded differently. It is the same factor
  // the light network applies to its He photo rates.
  kHiSelfShieldingGrey = 1,
};

struct MetalIonEnvironment {
  double temperature_k = 1.0e4;
  double electron_cm3 = 0.0;
  double redshift = 0.0;
  // Multiplicative attenuation of the dataset photo rates (1 = none). Set
  // only through ThermochemistrySolver::metalIonEnvironment().
  double photo_attenuation = 1.0;
  MetalIonPhotoAttenuationModel attenuation_model = MetalIonPhotoAttenuationModel::kNone;
};

struct MetalIonAdvanceReport {
  std::uint64_t accepted_substeps = 0U;
  std::uint64_t rejected_substeps = 0U;
  std::uint64_t range_clamps = 0U;
  double max_error_estimate = 0.0;
  // max |sum_k x_k - 1| removed by the final round-off projection.
  double normalization_residual = 0.0;
};

// Pure numerical kernel; no state, no allocation.
class MetalIonLadderSolver {
 public:
  MetalIonLadderSolver(const MetalIonRateTable& table, MetalIonRangePolicy policy, double relative_tolerance,
                       double absolute_tolerance, std::uint32_t max_substeps);

  // Advance the stage fractions of one element over exactly dt_seconds at a
  // frozen environment (operator-split half of the thermochemistry
  // interval). Returns kCompleted or kRetryShorterInterval (never partial).
  [[nodiscard]] SourceAdvanceStatus advance(MetalIonElement element, std::span<double> fractions,
                                            const MetalIonEnvironment& env, double dt_seconds,
                                            MetalIonAdvanceReport& report) const;
  // Ionization equilibrium (collisional + photo) at env; sum = 1.
  void equilibrium(MetalIonElement element, std::span<double> fractions, const MetalIonEnvironment& env,
                   MetalIonAdvanceReport& report) const;
  // Volumetric radiative loss of one element [erg cm^-3 s^-1].
  [[nodiscard]] double cooling(MetalIonElement element, std::span<const double> stage_cm3,
                               const MetalIonEnvironment& env, MetalIonAdvanceReport& report) const;
  // Volumetric photoheating of one element [erg cm^-3 s^-1].
  [[nodiscard]] double photoheating(MetalIonElement element, std::span<const double> stage_cm3,
                                    const MetalIonEnvironment& env, MetalIonAdvanceReport& report) const;
  [[nodiscard]] const MetalIonRateTable& table() const noexcept { return *m_table; }

 private:
  void rates(MetalIonElement element, const MetalIonEnvironment& env, std::span<double> ionize,
             std::span<double> recombine, MetalIonAdvanceReport& report) const;
  static void backwardEulerStep(std::span<const double> ionize, std::span<const double> recombine, double dt,
                                std::span<double> x);

  const MetalIonRateTable* m_table = nullptr;
  MetalIonRangePolicy m_policy = MetalIonRangePolicy::kReject;
  double m_relative_tolerance = 1.0e-4;
  double m_absolute_tolerance = 1.0e-12;
  std::uint32_t m_max_substeps = 4096U;
};

// Free electrons per nucleus, sum_k k x_k.
[[nodiscard]] double metalIonElectronsPerNucleus(std::span<const double> fractions) noexcept;

}  // namespace cosmosim::physics
