#pragma once

// CHUI-RT-BH-AGN-REF-v1, RT companion section 2 (RT02): M1 closure algebra.
//
//   f = |F| / (c E),  chi(f) = (3 + 4 f^2) / (5 + 2 sqrt(4 - 3 f^2))   (Levermore 1984)
//   P = E [ (1 - chi)/2 I + (3 chi - 1)/2 n n ],  n = F/|F|
// Degeneracies: F = 0 -> P = E I/3 exactly; E = 0 -> F must be 0 (P = 0).
// Admissible set: E >= 0, |F| <= c E (convex cone).
//
// Characteristic speed bounds (HLL): the exact eigenvalues of the M1 flux
// Jacobian along a face normal are computed from the analytic 3x3 Jacobian in
// the frame (normal, in-plane transverse) plus the decoupled out-of-plane
// eigenvalue c (3 chi - 1)/(2 f) mu (mu = n . normal); the cubic is solved in
// closed (trigonometric) form. Bounds are clamped to [-c, c] (exact for
// admissible states up to round-off). No ad hoc opacity wave-speed correction.

#include <array>

namespace cosmosim::radiation {

struct M1State {
  double energy = 0.0;
  std::array<double, 3> flux{};
};

[[nodiscard]] double m1EddingtonFactor(double reduced_flux) noexcept;
[[nodiscard]] double m1EddingtonFactorDerivative(double reduced_flux) noexcept;
// Radiation pressure tensor (row-major 3x3) for speed c.
[[nodiscard]] std::array<double, 9> m1PressureTensor(const M1State& state, double c) noexcept;
// |F| <= c E (1 + tol) and E >= 0, finite.
[[nodiscard]] bool m1Admissible(const M1State& state, double c, double tolerance = 1.0e-12) noexcept;

struct M1WaveSpeeds {
  double min = 0.0;
  double max = 0.0;
};
// Minimal/maximal characteristic speeds of the M1 system along unit normal
// `axis` (0,1,2) for speed c.
[[nodiscard]] M1WaveSpeeds m1WaveSpeeds(const M1State& state, int axis, double c) noexcept;

// Physical flux of (E, F) along axis: (F_axis, c^2 P_axis,*).
[[nodiscard]] std::array<double, 4> m1PhysicalFlux(const M1State& state, int axis, double c) noexcept;

// Largest theta in [0,1] such that base + theta * delta is admissible
// (base admissible). Closed form from the quadratic of the cone boundary.
[[nodiscard]] double m1AdmissibleFraction(const M1State& base, const M1State& delta, double c) noexcept;

}  // namespace cosmosim::radiation
