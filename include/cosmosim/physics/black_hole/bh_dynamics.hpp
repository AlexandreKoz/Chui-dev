#pragma once

// CHUI-RT-BH-AGN-REF-v1, BH companion section 3.2 (BH02): unresolved
// collisionless dynamical friction on top of the TreePM orbit.
//
// Chandrasekhar (1943) for an isotropic Maxwellian background:
//   a_df = -4 pi G^2 M rho lnL [erf(X) - 2X e^{-X^2}/sqrt(pi)] v / |v|^3,
//   X = |v| / (sqrt(2) sigma),
// small-X series bracket = 4X^3/(3 sqrt pi) (1 - 3X^2/5 + 3X^4/14 - X^6/18).
// Resolution-aware unresolved range (Tremmel et al. 2015, MNRAS 451, 1868):
// the TreePM force resolves the wake beyond the gravitational softening, so
// the correction uses b_max = softening (proper) and b_min = G M / v^2;
// lnL = ln(b_max/b_min) when positive, else ZERO (no correction when the
// strong-scattering radius is resolved). Dark matter and stars are separate
// backgrounds (rho_k, sigma_k, bulk v_k); gas drag needs a Mach-dependent
// closure and is NOT this formula (not implemented: gas is excluded).
// This is a proposed resolution-aware policy, not a calibrated friction law.

#include <array>
#include <cstdint>

namespace cosmosim::physics {

[[nodiscard]] double chandrasekharBracket(double x) noexcept;

struct BhFrictionBackground {
  double density_proper = 0.0;      // rho of the species (proper code density)
  double sigma_1d = 0.0;            // 1-D velocity dispersion
  std::array<double, 3> bulk_velocity{};
  std::uint64_t samples = 0U;
};

struct BhFrictionInputs {
  double bh_mass_code = 0.0;
  double newton_g_code = 0.0;
  double softening_proper = 0.0;    // b_max
  double timestep_eta = 0.1;        // |dv| <= eta |v| per substep
  std::uint32_t max_substeps = 100000U;
};

// Acceleration (code units) of a BH moving with velocity v relative to one
// background. Returns zero when lnL <= 0 or the background is empty.
[[nodiscard]] std::array<double, 3> chandrasekharAcceleration(const std::array<double, 3>& v_rel,
                                                              const BhFrictionBackground& background,
                                                              const BhFrictionInputs& inputs) noexcept;
[[nodiscard]] double tremmelCoulombLogarithm(double softening_proper, double newton_g, double bh_mass,
                                             double speed) noexcept;

struct BhFrictionResult {
  std::array<double, 3> velocity_change{};  // applied to the BH (code velocity)
  double timestep_limit_code = 0.0;         // eta |v| / |a| at entry (diagnostic / scheduler hint)
  std::uint32_t substeps = 0U;
  bool applied = false;
};

// Integrate dv/dt = sum_k a_k(v - v_k) over dt with frozen backgrounds (RK2
// midpoint substeps bounded by eta |v_rel|/|a|; the relative speed never
// reverses sign because friction only decelerates).
[[nodiscard]] BhFrictionResult integrateDynamicalFriction(const std::array<double, 3>& bh_velocity,
                                                          const BhFrictionBackground& dark_matter,
                                                          const BhFrictionBackground& stars,
                                                          const BhFrictionInputs& inputs, double dt_code);

}  // namespace cosmosim::physics
