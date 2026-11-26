#pragma once

// CHUI-RT-BH-AGN-REF-v1, RT companion section 3 (RT01): spectral contract.
//
// Groups are defined by PHYSICAL edges (eV). Profiles (no hardcoded count):
//   grey      : one bolometric thermal group (RHD reference; not an ionization model)
//   ionizing3 : [13.598434, 24.587389), [24.587389, 54.417765), [54.417765, upper)
//               (ionization potentials of H I, He I, He II; the rounded spec
//               table values 13.6/24.6/54.4 are replaced by the thresholds
//               themselves so no group straddles a threshold)
//   broad_uv_ir: 0.01,1,5.8,11.2,13.6,24.6,54.4,100 (defined; admission FAIL_CLOSED)
// Within-group shape s(nu) is DECLARED (SED table, or an analytic reference
// shape) and every derived group quantity comes from the same Gauss-Legendre
// quadrature in log(nu) (no magic constants):
//   <h nu>_g, sigma_E,g,s = int sigma s / int s, sigma_N,g,s = int sigma s/(h nu) / int s/(h nu),
//   eps_g,s = int sigma s / int sigma s/(h nu) (energy absorbed per ionization),
//   phi_low,g = s(nu_l)/int_g s dnu (lower-edge spectral density for the
//   conservative frequency-redshift flux).
// Identity used by the coupling (exact by construction):
//   sigma_N (E/<h nu>) eps = sigma_E E.
// The discrete photon diagnostic N_g = E_g/<h nu>_g is closure dependent and
// is NOT an independently conserved moment.
//
// Cross sections (cm^2):
//   H I, He II: exact nonrelativistic hydrogenic ground-state form
//     (Osterbrock & Ferland 2006 eq. 2.4), sigma_0 = 6.30e-18 / Z^2 cm^2;
//     verification: kDerivedFromConstants-equivalent (closed form).
//   He I: Verner et al. (1996, ApJ 465, 487) analytic fit, TRANSCRIBED
//     parameters (E0 = 13.61 eV, sigma0 = 949.2 Mb, ya = 1.469, P = 3.188,
//     yw = 2.039, y0 = 0.4434, y1 = 2.136); status kUnverifiedTranscription,
//     admitted only together with physics.thermochemistry_allow_unqualified_rates.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace cosmosim::radiation {

enum class RadiationAbsorber : std::uint8_t { kHI = 0, kHeI = 1, kHeII = 2 };
inline constexpr std::size_t k_radiation_absorber_count = 3U;

inline constexpr double k_ionization_potential_hi_ev = 13.598434;
inline constexpr double k_ionization_potential_hei_ev = 24.587389;
inline constexpr double k_ionization_potential_heii_ev = 54.417765;
inline constexpr double k_ev_erg = 1.602176634e-12;
inline constexpr double k_planck_erg_s = 6.62607015e-27;
inline constexpr double k_boltzmann_ev_per_k = 8.617333262e-5;

[[nodiscard]] double photoionizationCrossSectionCm2(RadiationAbsorber absorber, double photon_energy_ev) noexcept;
[[nodiscard]] double ionizationPotentialEv(RadiationAbsorber absorber) noexcept;

enum class SpectralShapeKind : std::uint8_t { kTable = 0, kBlackbody = 1, kPowerLaw = 2, kFlat = 3 };

// Tabulated L_nu shape (arbitrary normalization) on a strictly increasing
// photon-energy grid; log-log interpolation; NO extrapolation outside the
// table (queries outside the domain throw).
struct TabulatedSpectrum {
  std::string schema;             // "chui_ssp_sed_table_v1" / "chui_agn_sed_table_v1"
  std::string sha256;
  std::vector<double> energy_ev;  // strictly increasing, > 0
  std::vector<double> l_nu;       // >= 0
  [[nodiscard]] double evaluate(double energy_ev) const;
};

// Load a two-column (energy_eV, L_nu) text table with a '# schema: <name>'
// header; verifies the SHA-256 of the file bytes. Throws on any defect.
[[nodiscard]] TabulatedSpectrum loadTabulatedSpectrum(const std::string& path, const std::string& expected_sha256,
                                                      const std::string& expected_schema);

struct SpectralShape {
  SpectralShapeKind kind = SpectralShapeKind::kFlat;
  double blackbody_temperature_k = 0.0;
  double power_law_index = 0.0;  // L_nu ~ nu^-index
  const TabulatedSpectrum* table = nullptr;
  [[nodiscard]] double evaluate(double energy_ev) const;  // relative L_nu
  [[nodiscard]] std::string describe() const;
};

struct RadiationGroup {
  double lower_ev = 0.0;
  double upper_ev = 0.0;
  double mean_photon_energy_ev = 0.0;
  double lower_edge_density_per_ev = 0.0;  // phi_low (1/eV)
  std::array<double, k_radiation_absorber_count> sigma_energy_cm2{};
  std::array<double, k_radiation_absorber_count> sigma_number_cm2{};
  std::array<double, k_radiation_absorber_count> absorbed_energy_per_ionization_ev{};
  double shape_energy_integral = 0.0;  // int_g s dE (relative units) for source normalization
};

struct RadiationSpectralDefinition {
  std::string profile;                  // grey | ionizing3 | broad_uv_ir
  bool grey = false;
  std::vector<RadiationGroup> groups;
  std::vector<double> edges_ev;         // empty for grey
  std::uint32_t quadrature_points_per_group = 32U;
  std::string shape_description;
  std::string sha256;                   // canonical hash of edges + shape + quadrature + cross-section ids
  [[nodiscard]] std::size_t groupCount() const noexcept { return groups.size(); }
};

[[nodiscard]] std::vector<double> radiationProfileEdgesEv(const std::string& profile, double upper_edge_ev);
// Build groups for a profile and declared within-group shape.
[[nodiscard]] RadiationSpectralDefinition buildRadiationSpectralDefinition(const std::string& profile,
                                                                          double upper_edge_ev,
                                                                          const SpectralShape& shape,
                                                                          std::uint32_t quadrature_points = 32U);
// Fraction of a source's luminosity in each group for the given shape:
// int_g s dE / int_{all groups} s dE (spectral normalization of sources).
[[nodiscard]] std::vector<double> radiationGroupLuminosityFractions(const RadiationSpectralDefinition& definition,
                                                                    const SpectralShape& shape);

}  // namespace cosmosim::radiation
