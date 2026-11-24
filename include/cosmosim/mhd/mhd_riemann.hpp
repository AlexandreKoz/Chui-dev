#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string_view>

#include "cosmosim/mhd/mhd_geometry.hpp"

namespace cosmosim::mhd {

// Cell primitive state in comoving code units. bx/by/bz are the derived cell
// fields B_c (normal components reconstructed from face flux by the caller);
// the face-shared normal component passed to the Riemann solver always wins.
struct MhdPrimitiveState {
  double rho_comoving = 0.0;
  double vel_x_peculiar = 0.0;
  double vel_y_peculiar = 0.0;
  double vel_z_peculiar = 0.0;
  double pressure_comoving = 0.0;
  double bx_comoving = 0.0;
  double by_comoving = 0.0;
  double bz_comoving = 0.0;
  double metallicity_mass_fraction = 0.0;

  [[nodiscard]] double magneticEnergyDensity() const noexcept {
    return 0.5 * (bx_comoving * bx_comoving + by_comoving * by_comoving +
                  bz_comoving * bz_comoving);
  }
};

// Conserved state used by the ideal operator. The energy lane is the single
// authoritative total energy: thermal + kinetic + magnetic.
struct MhdConservedState {
  double rho_comoving = 0.0;
  double momentum_density_x_comoving = 0.0;
  double momentum_density_y_comoving = 0.0;
  double momentum_density_z_comoving = 0.0;
  double total_energy_density_comoving = 0.0;
  double bx_comoving = 0.0;
  double by_comoving = 0.0;
  double bz_comoving = 0.0;
  double metal_mass_density_comoving = 0.0;

  MhdConservedState& operator+=(const MhdConservedState& rhs) noexcept;
  MhdConservedState& operator-=(const MhdConservedState& rhs) noexcept;
};

[[nodiscard]] MhdConservedState operator+(MhdConservedState lhs,
                                          const MhdConservedState& rhs) noexcept;
[[nodiscard]] MhdConservedState operator-(MhdConservedState lhs,
                                          const MhdConservedState& rhs) noexcept;
[[nodiscard]] MhdConservedState operator*(double scalar, MhdConservedState state) noexcept;

// Numerical flux at one face. momentum_flux/magnetic_flux are absolute
// components of the flux vector along the face normal; electric_field holds the
// two tangential ideal electric-field components (E = -u x B); the normal
// component is identically zero and unused.
struct MhdRiemannFlux {
  double mass_flux = 0.0;
  std::array<double, 3> momentum_flux{};
  double energy_flux = 0.0;
  std::array<double, 3> magnetic_flux{};
  std::array<double, 3> electric_field{};
  double max_fast_speed = 0.0;
  double slowest_signal_speed = 0.0;
  double fastest_signal_speed = 0.0;
  bool used_hlle_fallback = false;
};

// Recover primitives from conserved state. Returns false for a non-positive or
// non-finite state so the caller can fall back before any floor is applied.
[[nodiscard]] bool recoverPrimitive(const MhdConservedState& conserved,
                                    double gamma_adiabatic,
                                    MhdPrimitiveState& primitive) noexcept;

[[nodiscard]] MhdConservedState primitiveToConserved(
    const MhdPrimitiveState& primitive,
    double gamma_adiabatic) noexcept;

// Sound/fast/slow signal speeds for a state with a prescribed normal field.
// `normal_sign` selects the direction (used only for reporting).
[[nodiscard]] double fastMagnetosonicSpeed(const MhdPrimitiveState& primitive,
                                           double gamma_adiabatic,
                                           MhdAxis axis) noexcept;
[[nodiscard]] double alfvenSpeedAlong(const MhdPrimitiveState& primitive,
                                      MhdAxis axis) noexcept;

// Ideal-MHD Riemann solvers. The normal magnetic field must be the single
// face-shared value on both sides (enforced by the caller).
class MhdRiemann {
 public:
  enum class Solver : std::uint8_t { kHlld = 0, kHlle = 1 };

  struct Result {
    MhdRiemannFlux flux{};
    bool valid = false;
    bool used_hlle = false;
    bool first_order_required = false;
  };

  [[nodiscard]] static Result solve(MhdAxis axis,
                                    const MhdPrimitiveState& left,
                                    const MhdPrimitiveState& right,
                                    double shared_normal_b_comoving,
                                    double gamma_adiabatic,
                                    Solver preferred = Solver::kHlld) noexcept;

  // Explicit HLLE entry point (also used as the HLLD degeneracy fallback).
  [[nodiscard]] static Result solveHlle(MhdAxis axis,
                                        const MhdPrimitiveState& left,
                                        const MhdPrimitiveState& right,
                                        double shared_normal_b_comoving,
                                        double gamma_adiabatic) noexcept;
};

// First-order positivity-preserving local fallback used only after HLLE has
// also failed. It is a documented Rusanov/Lax-Friedrichs flux over the
// supplied states; it is never a silent zero flux.
[[nodiscard]] MhdRiemannFlux firstOrderPositivityFlux(
    MhdAxis axis,
    const MhdPrimitiveState& left,
    const MhdPrimitiveState& right,
    double shared_normal_b_comoving,
    double gamma_adiabatic) noexcept;

// Primitive/axis permutation helpers exposed for tests.
[[nodiscard]] std::array<double, 3> permuteToLocal(MhdAxis axis,
                                                   const std::array<double, 3>& global) noexcept;
[[nodiscard]] std::array<double, 3> permuteToGlobal(MhdAxis axis,
                                                    const std::array<double, 3>& local) noexcept;

}  // namespace cosmosim::mhd
