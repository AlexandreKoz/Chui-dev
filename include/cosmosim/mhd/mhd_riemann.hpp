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

// Numerical flux at one face. momentum_flux/magnetic_flux are GLOBAL (x,y,z)
// components of the flux vector along the face normal. electric_field holds the
// ideal electric field in the LOCAL right-handed basis {n, t1, t2} with
// t1/t2 = cyclic successors of the normal axis (x:(y,z), y:(z,x), z:(x,y));
// E_n = 0, E_t1 = F_B,t2, E_t2 = -F_B,t1. slowest/fastest_signal_speed are the
// SIGNED outer fan bounds S_L <= S_R used by HLL upwinding.
//
// Cosmological runs pass B_eff = B_c/sqrt(a) in the primitive states and the
// shared normal field; the returned fluxes are then the comoving-volume
// gas+magnetic fluxes and electric_field is E_eff = -(u x B_eff). The CT EMF
// for B_c is E_eff/sqrt(a) (CHUI-MHD-REF v1.0 section 4.3).
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

// Signed HLL fan bounds at one face (Davis/Einfeldt estimates).
struct MhdSignalBounds {
  double speed_left = 0.0;
  double speed_right = 0.0;
  bool valid = false;
};

// Ideal-MHD Riemann solvers. The normal magnetic field must be the single
// face-shared value on both sides (enforced by the solver, which overrides the
// normal component of both input states).
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

  // Signed outer fan bounds only (used for EMF upwinding on ghost-row faces
  // whose gas flux is not applied).
  [[nodiscard]] static MhdSignalBounds signalBounds(MhdAxis axis,
                                                    const MhdPrimitiveState& left,
                                                    const MhdPrimitiveState& right,
                                                    double shared_normal_b_comoving,
                                                    double gamma_adiabatic) noexcept;

  // Explicit HLLE entry point (also used as the HLLD degeneracy fallback).
  [[nodiscard]] static Result solveHlle(MhdAxis axis,
                                        const MhdPrimitiveState& left,
                                        const MhdPrimitiveState& right,
                                        double shared_normal_b_comoving,
                                        double gamma_adiabatic) noexcept;
};

// First-order positivity-preserving local fallback used only after HLLE has
// also failed on the reconstructed states. It is a Rusanov/local Lax-Friedrichs
// flux between the two adjacent CELL-AVERAGE states (not the rejected
// reconstructed states) with the shared face normal field. Returns
// valid=false for an inadmissible pair; it never returns a silent zero flux.
[[nodiscard]] MhdRiemann::Result firstOrderPositivityFlux(
    MhdAxis axis,
    const MhdPrimitiveState& left_cell_average,
    const MhdPrimitiveState& right_cell_average,
    double shared_normal_b_comoving,
    double gamma_adiabatic) noexcept;

// Primitive/axis permutation helpers exposed for tests.
[[nodiscard]] std::array<double, 3> permuteToLocal(MhdAxis axis,
                                                   const std::array<double, 3>& global) noexcept;
[[nodiscard]] std::array<double, 3> permuteToGlobal(MhdAxis axis,
                                                    const std::array<double, 3>& local) noexcept;

}  // namespace cosmosim::mhd
