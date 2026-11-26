#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 7 (P1-BH02/P1-BH03): Kerr
// thin-disk accretion, joint mass/angular-momentum evolution and
// Blandford-Znajek-informed spin extraction.
//
// STATUS: P1 source implementation. Production admission refuses
// physics.bh_spin_model != none (no qualification; published warp/alignment
// coefficients not pinned). Available to the explicitly named reference API
// and tests only.
//
// Kerr ISCO (Bardeen, Press & Teukolsky 1972), signed a relative to the disk
// angular momentum (a < 0 counter-aligned), r in GM/c^2:
//   Z1 = 1 + (1-a^2)^{1/3}[(1+a)^{1/3} + (1-a)^{1/3}],  Z2 = sqrt(3a^2 + Z1^2)
//   r_isco = 3 + Z2 - sgn(a) sqrt((3-Z1)(3+Z1+2Z2))
//   E_isco = sqrt(1 - 2/(3 r_isco)),  eps = 1 - E_isco
//   l_isco = (r^2 - 2a sqrt(r) + a^2)/(r^{3/4} sqrt(r^{3/2} - 3 sqrt(r) + 2a))  [GM/c]
// Ideal thin disk: dM = E_isco dm0, dJ = l_isco (G M / c) dm0 along the disk
// axis. Integrated in dm0 substeps bounded by a relative mass increment; the
// equilibrium a* -> 1 is approached without overshoot; an update that would
// exceed |a*| = 1 is NOT clipped silently: the excess angular momentum is
// returned as a named saturation term (photon capture is not modelled; the
// Thorne a* = 0.998 cap is a separate, not-implemented model).
//
// BZ jet (Tchekhovskoy, Narayan & McKinney 2010 high-spin form, Gaussian):
//   P = (kappa / (4 pi c)) Phi^2 Omega_H^2 f(Omega_H),  f = 1 + 1.38 x^2 - 9.2 x^4,
//   x = Omega_H r_g / c,  r_H = r_g (1 + sqrt(1 - a^2)),  Omega_H = a c / (2 r_H)
// (coefficients 1.38, -9.2 and kappa are TRANSCRIBED, verification status
// kUnverifiedTranscription). Extraction removes energy AND angular momentum,
// dE = P dt, dJ = P/Omega_F dt with Omega_F = omega_F Omega_H, 0 < omega_F < 1,
// so the horizon area (irreducible mass) never decreases (checked).

#include <array>
#include <cstdint>

namespace cosmosim::physics {

[[nodiscard]] double kerrIscoRadius(double signed_spin);
[[nodiscard]] double kerrIscoEnergy(double signed_spin);
[[nodiscard]] double kerrRadiativeEfficiency(double signed_spin);
[[nodiscard]] double kerrIscoSpecificAngularMomentum(double signed_spin);
[[nodiscard]] double kerrIrreducibleMass(double mass, double spin_magnitude) noexcept;
[[nodiscard]] double kerrHorizonAngularVelocity(double spin_magnitude, double gravitational_radius,
                                                double speed_of_light) noexcept;

struct KerrAccretionInputs {
  double mass_code = 0.0;
  std::array<double, 3> spin_angular_momentum_code{};
  std::array<double, 3> disk_axis{};  // unit vector of the disk angular momentum
  double inner_mass_code = 0.0;       // dm0 reaching the ISCO
  double newton_g_code = 0.0;
  double speed_of_light_code = 0.0;
  double max_relative_mass_step = 1.0e-3;
};

struct KerrAccretionResult {
  double mass_code = 0.0;
  std::array<double, 3> spin_angular_momentum_code{};
  double radiated_energy_code = 0.0;            // (1 - E_isco) dm0 c^2 integrated
  std::array<double, 3> saturation_angular_momentum_code{};  // named excess at |a*| = 1
  std::uint32_t substeps = 0U;
  bool misaligned_rejected = false;            // warp model required (not implemented)
};

// Aligned/counter-aligned ideal thin-disk accretion. Misaligned flows
// (|cos theta| < 1 - 1e-9 with nonzero spin) are rejected: a warp model is
// required and none is pinned.
[[nodiscard]] KerrAccretionResult accreteKerrThinDisk(const KerrAccretionInputs& inputs);

struct BlandfordZnajekInputs {
  double mass_code = 0.0;
  std::array<double, 3> spin_angular_momentum_code{};
  double magnetic_flux_gaussian_code = 0.0;  // Phi_BH in the declared (Gaussian) flux convention
  double kappa = 0.053;                      // geometry factor (split monopole ~ 1/(6 pi)); transcription
  double field_angular_velocity_ratio = 0.5; // omega_F = Omega_F / Omega_H in (0, 1)
  double newton_g_code = 0.0;
  double speed_of_light_code = 0.0;
  double dt_code = 0.0;
};

struct BlandfordZnajekResult {
  double extracted_energy_code = 0.0;
  std::array<double, 3> extracted_angular_momentum_code{};  // removed from the BH (along the spin axis)
  double mass_code = 0.0;                                   // after extraction
  std::array<double, 3> spin_angular_momentum_code{};
  double irreducible_mass_before_code = 0.0;
  double irreducible_mass_after_code = 0.0;
  bool area_theorem_satisfied = true;
};

[[nodiscard]] BlandfordZnajekResult extractBlandfordZnajek(const BlandfordZnajekInputs& inputs);

}  // namespace cosmosim::physics
