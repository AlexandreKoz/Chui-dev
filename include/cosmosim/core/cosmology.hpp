#pragma once

#include <cstddef>
#include <vector>

namespace cosmosim::core {

// Relativistic species of the homogeneous background (CHUI-EXT-PHYS-REF-v1
// S01/N01). Disabled by default: the legacy LambdaCDM formula (including the
// legacy `omega_radiation` lane, which no config key feeds) is reproduced
// bit-for-bit when `enabled == false`.
//
// When enabled:
//   * photons:      Omega_gamma,0 = a_rad T_cmb0^4 / (c^2 rho_crit,0), rho ~ a^-4
//   * massless nu:  Omega_ur,0 = n_ur (7/8) (4/11)^(4/3) Omega_gamma,0 (N_eff
//                   convention: each unit is one fermionic species with
//                   g = 2 (nu + nubar) at the instantaneous-decoupling
//                   temperature (4/11)^(1/3) T_cmb).
//   * massive nu:   each eigenstate i (mass m_i, g = 2) follows the exact
//                   relativistic Fermi-Dirac moments with
//                   T_nu(a) = ratio T_cmb0 / a  (no chemical potential):
//                     rho_i c^2 = g/(2 pi^2) (k T_nu)^4/(hbar c)^3 I_rho(y),
//                     P_i       = g/(2 pi^2) (k T_nu)^4/(hbar c)^3 I_P(y),
//                     y = m_i c^2 / (k T_nu) = y_i,0 a,
//                     I_rho(y) = int x^2 sqrt(x^2+y^2)/(e^x+1) dx,
//                     I_P(y)   = int x^4 / (3 sqrt(x^2+y^2)) /(e^x+1) dx.
//   * Omega_matter (the existing `cosmology.omega_matter` / GADGET Omega0)
//     is the TOTAL present-day clustering density INCLUDING the present
//     energy density of massive eigenstates (CLASS convention
//     Omega_m = Omega_b + Omega_cdm + Omega_ncdm). The cold part that the
//     N-body particles represent is Omega_cb = Omega_matter - Omega_ncdm,0,
//     so massive neutrinos are never double counted:
//       E^2(a) = (Omega_gamma + Omega_ur) a^-4 + Omega_cb a^-3
//              + sum_i Omega_ncdm,i(a) + Omega_k a^-2 + Omega_lambda.
//   * Closure is strict: Omega_gamma + Omega_ur + Omega_matter + Omega_k +
//     Omega_lambda = 1 within 1e-10; the error message names the required
//     Omega_lambda. Nothing is silently renormalized.
struct RadiationSpeciesConfig {
  bool enabled = false;
  double t_cmb0_k = 0.0;            // present CMB temperature [K], required > 0
  double n_ur = 0.0;                // massless effective species (N_eff convention)
  std::vector<double> ncdm_masses_ev;  // one entry per massive eigenstate (g = 2 each)
  double ncdm_temperature_ratio = 0.0;  // T_ncdm,0 / T_cmb,0 (> 0 when masses present)
  bool operator==(const RadiationSpeciesConfig&) const = default;
};

// LambdaCDM background parameters.
// hubble_param is the dimensionless little-h: H0 = h * 100 km/s/Mpc.
struct CosmologyBackgroundConfig {
  double hubble_param = 0.674;  // h = H0 / (100 km/s/Mpc)
  double omega_matter = 0.315;
  double omega_lambda = 0.685;
  double omega_radiation = 0.0;  // legacy lane; must stay 0 when radiation_species.enabled
  double omega_curvature = 0.0;
  RadiationSpeciesConfig radiation_species;
};

// Instantaneous-decoupling neutrino/photon temperature ratio (4/11)^(1/3)
// from entropy conservation through e+e- annihilation (derived, not fitted).
[[nodiscard]] double instantaneousDecouplingNeutrinoTemperatureRatio();

// Dimensionless Fermi-Dirac moment integrals (see RadiationSpeciesConfig).
// The *Quadrature variants evaluate the integral directly (composite
// Gauss-Legendre on [0, 80], 40 panels x 16 nodes); the plain variants use a
// shared read-only log-spaced table (y in [1e-5, 1e6]) with 4-point Lagrange
// interpolation in (ln y, ln I) and verified two-term asymptotic expansions
// outside the table. Both throw for y < 0 or non-finite y.
[[nodiscard]] double fermiDiracEnergyIntegral(double y);
[[nodiscard]] double fermiDiracPressureIntegral(double y);
[[nodiscard]] double fermiDiracNumberIntegral();  // int x^2/(e^x+1) = (3/2) zeta(3)
[[nodiscard]] double fermiDiracEnergyIntegralQuadrature(double y);
[[nodiscard]] double fermiDiracPressureIntegralQuadrature(double y);
// Mean |p| c / (k T) of the relativistic Fermi-Dirac distribution,
// <x> = int x^3 f / int x^2 f = (7 pi^4/120) / ((3/2) zeta(3)).
[[nodiscard]] double fermiDiracMeanMomentumOverTemperature();

class LambdaCdmBackground {
 public:
  explicit LambdaCdmBackground(CosmologyBackgroundConfig config = {});

  // Hubble constant at z=0 in SI units [s^-1].
  [[nodiscard]] const CosmologyBackgroundConfig& config() const;
  [[nodiscard]] double hubble0Si() const;
  // E(a) = H(a) / H0.
  [[nodiscard]] double eFactor(double scale_factor) const;
  // H(a) = H0 * E(a); legacy formula when radiation species are disabled.
  [[nodiscard]] double hubbleSi(double scale_factor) const;
  // rho_crit(a) = 3 H(a)^2 / (8 pi G), SI units [kg m^-3].
  [[nodiscard]] double criticalDensitySi(double scale_factor) const;
  // Exact FLRW elapsed proper time: integral da / (a H(a)), SI seconds.
  [[nodiscard]] double cosmicTimeIntervalSi(
      double scale_factor_begin,
      double scale_factor_end) const;
  // Proper cosmic time since a configurable positive lower bound.
  [[nodiscard]] double cosmicTimeSinceScaleFactorSi(
      double scale_factor,
      double lower_scale_factor = 1.0e-8) const;

  // ---------------------------------------------------- S01 species accounting
  [[nodiscard]] bool radiationSpeciesEnabled() const noexcept;
  [[nodiscard]] double criticalDensity0Si() const;
  [[nodiscard]] double photonDensityParameter0() const noexcept;            // Omega_gamma,0
  [[nodiscard]] double masslessNeutrinoDensityParameter0() const noexcept;  // Omega_ur,0
  [[nodiscard]] std::size_t massiveNeutrinoCount() const noexcept;
  // rho_ncdm,i(a) / rho_crit,0 (energy density, includes rest mass).
  [[nodiscard]] double massiveNeutrinoDensityParameter(std::size_t eigenstate, double scale_factor) const;
  // P_ncdm,i(a) / (c^2 rho_crit,0).
  [[nodiscard]] double massiveNeutrinoPressureParameter(std::size_t eigenstate, double scale_factor) const;
  [[nodiscard]] double totalMassiveNeutrinoDensityParameter(double scale_factor) const;
  [[nodiscard]] double totalMassiveNeutrinoPressureParameter(double scale_factor) const;
  // Omega_cb = Omega_matter - sum_i Omega_ncdm,i(1) (cold matter carried by particles).
  [[nodiscard]] double coldMatterDensityParameter0() const noexcept;
  // y_i(a) = m_i c^2 / (k T_ncdm(a)).
  [[nodiscard]] double massiveNeutrinoMassOverTemperature(std::size_t eigenstate, double scale_factor) const;
  // <|p|> c / (m_i c^2) = <x> / y_i(a): mean relativistic factor for validity checks.
  [[nodiscard]] double massiveNeutrinoMeanMomentumOverMass(std::size_t eigenstate, double scale_factor) const;
  // Scale factor of matter-radiation equality: (Omega_gamma + Omega_ur) a^-4 +
  // relativistic part of ncdm equals matter; bisection on the same background.
  [[nodiscard]] double matterRadiationEqualityScaleFactor() const;

  // ---------------------------------------------------- S01 CMB background
  // Standard adiabatic evolution T_cmb(a) = T_cmb0 / a; throws when disabled.
  [[nodiscard]] double cmbTemperatureK(double scale_factor) const;
  // u_cmb = a_rad T^4 [J m^-3] (proper).
  [[nodiscard]] double cmbEnergyDensitySi(double scale_factor) const;
  // n_gamma = (2 zeta(3) / pi^2) (k T / (hbar c))^3 [m^-3] (proper).
  [[nodiscard]] double cmbPhotonNumberDensitySi(double scale_factor) const;

  // ---------------------------------------------------- S01 distances (SI m)
  // Line-of-sight comoving distance chi(a) = c int_a^1 da' / (a'^2 H(a')).
  [[nodiscard]] double comovingDistanceSi(double scale_factor) const;
  // Transverse comoving distance (curvature aware).
  [[nodiscard]] double transverseComovingDistanceSi(double scale_factor) const;
  [[nodiscard]] double angularDiameterDistanceSi(double scale_factor) const;
  [[nodiscard]] double luminosityDistanceSi(double scale_factor) const;
  // Conformal-time interval eta_end - eta_begin = int dt / a = int da/(a^2 H), SI s.
  [[nodiscard]] double conformalTimeIntervalSi(double scale_factor_begin, double scale_factor_end) const;
  // Inverse of comovingDistanceSi on [a_min, 1] (bisection on the monotone map).
  [[nodiscard]] double scaleFactorAtComovingDistanceSi(double comoving_distance_si, double a_min = 1.0e-3) const;

 private:
  CosmologyBackgroundConfig m_config;
  double m_omega_gamma0 = 0.0;
  double m_omega_ur0 = 0.0;
  double m_omega_ncdm_unit0 = 0.0;    // massless-equivalent Omega of one eigenstate at a = 1
  std::vector<double> m_ncdm_y0;      // y_i at a = 1
  double m_omega_ncdm0 = 0.0;         // sum_i Omega_ncdm,i(a = 1)
  double m_omega_cb0 = 0.0;
};

}  // namespace cosmosim::core
