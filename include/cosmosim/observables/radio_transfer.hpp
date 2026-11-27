#pragma once

// CHUI-EXT-PHYS-REF-v1 RA03-RA05: radio continuum emission/absorption,
// Faraday propagation, full-Stokes formal transfer and DM/RM.
// All coefficients are evaluated in Gaussian-CGS internally (constants
// DERIVED from the SI-defining ones: e_esu = e c[cm/s]/10) and returned in SI:
// j [W m^-3 Hz^-1 sr^-1], alpha [m^-1], rotation rates [rad m^-1].
//
// Free-free (Rybicki & Lightman eq. 5.14; recalled lead):
//   eps_ff = (2^5 pi e^6 / (3 m c^3)) (2 pi / (3 k m))^(1/2) Z^2 n_e n_i T^-1/2 e^(-h nu/kT) g_ff,
//   j = eps_ff / (4 pi); alpha = j / B_nu(T) (Kirchhoff, LTE);
//   g_ff from Draine (2011) eq. 10.9: g = ln{exp[5.960 - (sqrt3/pi) ln(Z nu_9 T_4^-3/2)] + e}
//   (recalled lead; flagged for verification).
// Synchrotron (RL eqs. 6.17c, 6.18, 6.31, 6.50; recalled leads):
//   single particle P(nu) = sqrt(3) q^3 B sin(alpha) / (m c^2) F(nu / nu_c),
//   nu_c = (3 / 4 pi) gamma^2 q B sin(alpha) / (m c),
//   F(x) = x int_x^inf K_5/3 = x int_0^inf cosh(5t/3) e^(-x cosh t) / cosh t dt (derived identity),
//   G(x) = x K_2/3(x) = x int_0^inf e^(-x cosh t) cosh(2t/3) dt,
//   j_I = (1/4 pi) int N(gamma) P dgamma, linear fraction Pi = int N P_G / int N P_F,
//   alpha_I = -(c^2 / 8 pi nu^2) int dE P(nu,E) E^2 d/dE[N(E)/E^2].
//   alpha is the angle between B and the line of sight (isotropic electrons).
//   Power-law closed forms (RL 6.36, 6.53, 6.38) are provided as an
//   INDEPENDENT comparator of the kernel integration. Circular emission,
//   Faraday conversion and polarized absorption (alpha_Q, alpha_U) are named
//   omissions in pass 01.
// Faraday rotation (cold plasma, SI): d(2 chi)/ds = rho_V,
//   rho_V = e^3 lambda^2 n_e B_par / (4 pi^2 eps0 m_e^2 c^3),
//   B_par = B . k (k = propagation direction toward the observer): positive
//   RM for B pointing toward the observer (IAU convention).
// Stokes transfer dS/ds = j - K S, S = (I, Q, U, V), with
//   K = [[aI, aQ, aU, aV], [aQ, aI, rV, rU], [aU, -rV, aI, rQ], [aV, -rU, -rQ, aI]];
//   constant-coefficient segments are solved EXACTLY by the 5x5 augmented
//   exponential exp([[-K, j], [0, 0]] ds) (Pade(6,6) with scaling and
//   squaring, coefficients computed), valid for singular K (pure rotation).
// DM/RM (observer frame): DM = int n_e dl / (1 + z), RM = C int n_e B_par dl / (1 + z)^2,
//   C = e^3 / (8 pi^2 eps0 m_e^2 c^3).

#include <array>
#include <cstddef>
#include <vector>

namespace cosmosim::observables {

using Stokes = std::array<double, 4>;
using StokesMatrix = std::array<double, 16>;  // row-major 4 x 4

[[nodiscard]] double planckBnu(double nu_hz, double t_k);                 // W m^-2 Hz^-1 sr^-1
[[nodiscard]] double rayleighJeansTemperature(double intensity_si, double nu_hz);
[[nodiscard]] double freeFreeGaunt(double nu_hz, double t_k, double z_charge);
[[nodiscard]] double freeFreeEmissivity(double nu_hz, double t_k, double n_e_m3, double z2_n_i_m3);
[[nodiscard]] double freeFreeAbsorption(double nu_hz, double t_k, double n_e_m3, double z2_n_i_m3);

[[nodiscard]] double synchrotronKernelF(double x);
[[nodiscard]] double synchrotronKernelG(double x);
[[nodiscard]] double synchrotronCriticalFrequency(double gamma, double b_tesla, double sin_alpha);

struct ElectronSpectrumSample {
  std::vector<double> gamma;          // strictly increasing Lorentz factors
  std::vector<double> number_per_m3;  // N(gamma) [m^-3 per unit gamma]
};

struct SynchrotronCoefficients {
  double j_i = 0.0;
  double linear_fraction = 0.0;
  double alpha_i = 0.0;
};
[[nodiscard]] SynchrotronCoefficients synchrotronFromSpectrum(const ElectronSpectrumSample& electrons, double nu_hz,
                                                              double b_tesla, double sin_alpha);
// Power-law N(gamma) = C gamma^-p (C in m^-3): closed forms (comparator).
[[nodiscard]] double synchrotronPowerLawEmissivity(double c_m3, double p, double nu_hz, double b_tesla, double sin_alpha);
[[nodiscard]] double synchrotronPowerLawAbsorption(double c_m3, double p, double nu_hz, double b_tesla, double sin_alpha);
[[nodiscard]] double synchrotronPowerLawPolarization(double p);

[[nodiscard]] double faradayRotationRate(double wavelength_m, double n_e_m3, double b_parallel_tesla);
[[nodiscard]] double rotationMeasureCoefficient();  // C [rad m^-2 per (m^-3 T m)]

[[nodiscard]] StokesMatrix propagationMatrix(double alpha_i, double alpha_q, double alpha_u, double alpha_v,
                                             double rho_q, double rho_u, double rho_v);
// Exact constant-coefficient step over length ds.
[[nodiscard]] Stokes stokesStep(const Stokes& s_in, const Stokes& emission, const StokesMatrix& k, double ds);

struct DispersionSegment {
  double path_proper_m = 0.0;
  double electron_density_m3 = 0.0;
  double b_parallel_tesla = 0.0;
  double redshift = 0.0;
};
struct DispersionMeasures {
  double dm_m2 = 0.0;   // [m^-2] (pc cm^-3 = 3.0857e22 m^-2)
  double rm_rad_m2 = 0.0;
};
[[nodiscard]] DispersionMeasures integrateDispersion(const std::vector<DispersionSegment>& segments);

}  // namespace cosmosim::observables
