#pragma once

// CHUI-EXT-PHYS-REF-v1 S02: shared read-only atomic line data and line
// profiles for Lyman-alpha (L01-L03) and the HI 21 cm line (RA01-RA02).
//
// Line constants carry an explicit provenance tag. Values marked
// `recalled_lead` are standard textbook/NIST values entered from the
// implementing agent's knowledge WITHOUT external verification in pass 01;
// they must be checked against the NIST Atomic Spectra Database before
// qualification (ledger item). No observable silently substitutes a
// different constant: every consumer reads this one table.
//
// Voigt profile: H(a, x) = Re w(x + i a), w the Faddeeva function, evaluated
// by Weideman's rational expansion (SIAM J. Numer. Anal. 31, 1497, 1994;
// recalled lead) with N = 32 terms, coefficients computed at start-up by a
// direct DFT (no transcribed coefficient table). H is normalized such that
// int H(a,x) dx = sqrt(pi). A direct-quadrature reference
// (voigtHQuadrature) is provided as an independent comparator.

#include <complex>
#include <string>

namespace cosmosim::spectral {

struct AtomicLine {
  std::string name;
  double rest_wavelength_m = 0.0;   // vacuum
  double rest_frequency_hz = 0.0;
  double oscillator_strength = 0.0;  // absorption f_lu
  double einstein_a_s = 0.0;         // spontaneous decay A_ul
  double upper_weight = 0.0;         // g_u
  double lower_weight = 0.0;         // g_l
  double atom_mass_kg = 0.0;
  std::string provenance;            // e.g. "recalled_lead: NIST ASD"
};

// H I Lyman-alpha 1s-2p (fine-structure doublet treated as one line with the
// summed oscillator strength): lambda = 1215.67 A, f = 0.4164,
// A = 6.265e8 s^-1, g_l = 2, g_u = 6 (recalled_lead, NIST ASD).
[[nodiscard]] const AtomicLine& hydrogenLymanAlpha();

// H I 21 cm hyperfine: nu = 1420.405751768 MHz, A_10 = 2.85e-15 s^-1,
// g_0 = 1, g_1 = 3 (recalled_lead; A_10 literature values range
// 2.85e-15 .. 2.88e-15 s^-1, a 1% systematic recorded in the ledger).
[[nodiscard]] const AtomicLine& hydrogen21cm();

// Classical absorption cross-section normalization
//   sigma(nu) = (pi e^2 / (m_e c)) f phi(nu)    [Gaussian electrostatic units]
//             = (e^2 / (4 eps0 m_e c)) f phi(nu) [SI], phi in s (per Hz).
[[nodiscard]] double lineIntegratedCrossSectionSi(const AtomicLine& line);  // m^2 Hz

// Doppler width Delta nu_D = (nu_0 / c) sqrt(2 k T / m) and damping
// parameter a = A / (4 pi Delta nu_D) (natural broadening only).
[[nodiscard]] double dopplerWidthHz(const AtomicLine& line, double temperature_k, double turbulent_b_m_s = 0.0);
[[nodiscard]] double naturalDampingParameter(const AtomicLine& line, double doppler_width_hz);

// Faddeeva w(z) for Im z >= 0 (Weideman N = 32).
[[nodiscard]] std::complex<double> faddeeva(std::complex<double> z);
// Voigt-Hjerting H(a, x) = Re w(x + i a), a >= 0; a = 0 returns exp(-x^2).
[[nodiscard]] double voigtH(double a, double x);
// Independent reference: H(a,x) = (a/pi) int exp(-y^2)/((x-y)^2 + a^2) dy by
// adaptive Gauss-Kronrod-free composite Gauss-Legendre over y in [-12, 12]
// with panel refinement around y = x; valid for a > 1e-6.
[[nodiscard]] double voigtHQuadrature(double a, double x);

// Normalized line profile phi(nu) [1/Hz] with int phi dnu = 1:
//   phi = H(a, x) / (sqrt(pi) Delta nu_D),  x = (nu - nu_0) / Delta nu_D.
[[nodiscard]] double voigtProfilePerHz(double a, double x, double doppler_width_hz);

}  // namespace cosmosim::spectral
