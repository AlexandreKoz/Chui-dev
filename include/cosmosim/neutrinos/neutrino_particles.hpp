#pragma once

// CHUI-EXT-PHYS-REF-v1 N04: collisionless massive-neutrino particle kernels
// (ordinary sampled particles and weighted delta-f markers). REFERENCE
// KERNELS ONLY in pass 01: `neutrinos.method = particles | delta_f` is
// refused at admission because the canonical particle model has no neutrino
// species (see config_extended_physics.cpp). Nothing here is reachable from
// the production runtime.
//
// Phase-space variables: comoving position x [m], canonical comoving momentum
// q = a p_phys [kg m/s] (conserved in the unperturbed background).
// Weak-field metric ds^2 = a^2[-(1+2Psi)c^2 dtau^2 + (1-2Phi)dx^2] with
// Phi = Psi (no anisotropic stress). From the point-particle Lagrangian
// L = -m c a sqrt((1+2Psi)c^2 - (1-2Phi)v^2) (derived in the numerical
// contracts, N04), to first order in the potential:
//   drift:  dx/dt = q c^2 / (a eps c^2) = q / (a eps),
//           eps c^2 = sqrt(q^2 c^2 + m^2 c^4 a^2)   (eps = a gamma m)
//   kick:   dq/dt = -(eps / a) (1 + v^2/c^2) grad_x Psi,   v/c = q c / (eps c^2)
// The d Phi / d t energy term and O(Psi) corrections to q are neglected
// (quasi-static weak field). Nonrelativistic limit: dq/dt = -m grad Psi.
//
// Delta-f (weights): f = f_ref + delta f with the analytic homogeneous
// reference f_ref(q) = 1/(exp(|q| c / (k T_nu,0)) + 1). Liouville: f is
// conserved along a trajectory, so a marker sampled from f_ref at q_i carries
// f_i = f_ref(q_i) and its weight is w = (f - f_ref(q)) / f = 1 - f_ref(q)/f_i
// (Elbers et al. 2021 formulation; recalled lead). The marker's
// gravitational contribution is w times its energy mass; the homogeneous
// background enters H(a) only (no extra mean density). Weights are SIGNED and
// are never clipped. Positive-mass tree multipoles/opening criteria are
// unsafe near cancellation, so the only admissible force model for signed
// weights here is the PM-only, explicitly scale-limited one
// (requireNonNegativeTreeMasses refuses signed masses on the tree path).

#include <array>
#include <cstdint>
#include <span>

#include "cosmosim/core/cosmology.hpp"

namespace cosmosim::neutrinos {

// Counter-based SplitMix64 uniform in (0, 1): reproducible across ranks and
// restarts for a given (seed, particle id, draw counter).
[[nodiscard]] double neutrinoCounterUniform(std::uint64_t seed, std::uint64_t particle_id, std::uint64_t counter) noexcept;

// Inverse CDF of the momentum magnitude x = |q| c / (k T) of the relativistic
// FD distribution x^2/(e^x + 1): bisection on the tabulated CDF (2048 nodes,
// x in [0, 40]) with linear refinement. u in (0, 1).
[[nodiscard]] double sampleFermiDiracMomentumMagnitude(double u);

struct NeutrinoMarker {
  std::array<double, 3> position_m{};
  std::array<double, 3> momentum_si{};  // canonical comoving q [kg m/s]
  double rest_mass_particle_kg = 0.0;   // ordinary marker rest mass (number of neutrinos x m_nu)
  double initial_reference_density = 0.0;  // f_i = f_ref(q_i) (delta-f marker)
  std::uint64_t particle_id = 0;
};

struct NeutrinoSpeciesScales {
  double mass_ev = 0.0;
  double temperature_k = 0.0;  // T_nu,0
  [[nodiscard]] double momentumScaleSi() const;  // k T_nu,0 / c
  [[nodiscard]] double massKg() const;
};

// Isotropic thermal sample (no bulk perturbation; IC displacement/velocity
// fields from transfer tables are an N02 integration edge, FAIL_CLOSED).
[[nodiscard]] NeutrinoMarker sampleThermalMarker(const NeutrinoSpeciesScales& scales, std::uint64_t seed,
                                                 std::uint64_t particle_id, const std::array<double, 3>& position_m,
                                                 double rest_mass_particle_kg);

// eps(q, a) = sqrt(q^2 c^2 + m^2 c^4 a^2) / c^2  [kg] (comoving energy / c^2).
[[nodiscard]] double comovingEnergyMass(const NeutrinoSpeciesScales& scales, double q_magnitude_si, double a);

// Exact drift over [a0, a1] for constant q: Delta x = q int dt / (a eps)
// = q int da / (a^2 H eps(q, a)); 16-point Gauss-Legendre in ln a.
void driftMarker(NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales,
                 const core::LambdaCdmBackground& background, double a0, double a1);

// Kick with a frozen gradient of the peculiar potential grad_x Psi [m s^-2]
// (comoving gradient of the potential in m^2/s^2): dq = -(eps/a)(1+v^2/c^2)
// grad Psi dt evaluated with q at the interval start (first-order, explicit;
// the KDK caller supplies the half-interval).
void kickMarker(NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales, const std::array<double, 3>& grad_psi,
                const core::LambdaCdmBackground& background, double a0, double a1);

// Delta-f weight w = 1 - f_ref(q) / f_i (signed, never clipped).
[[nodiscard]] double deltaFWeight(const NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales);
// Energy mass of a marker: rest mass x gamma(q, a) (eps / (m a)).
[[nodiscard]] double markerEnergyMass(const NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales, double a);
// Gravitational (signed) estimator mass of a delta-f marker: w x energy mass.
[[nodiscard]] double deltaFGravitationalMass(const NeutrinoMarker& marker, const NeutrinoSpeciesScales& scales, double a);

// Throws std::invalid_argument when any mass is negative: signed delta-f
// estimator masses must never enter the positive-mass tree multipoles.
void requireNonNegativeTreeMasses(std::span<const double> masses);

}  // namespace cosmosim::neutrinos
