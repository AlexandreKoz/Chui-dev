#pragma once

#include <cstdint>
#include <memory>
#include <string_view>

namespace cosmosim::hydro {

struct HydroConservedState;
struct HydroFace;
struct HydroPrimitiveState;

class HydroRiemannSolver {
 public:
  virtual ~HydroRiemannSolver() = default;

  [[nodiscard]] virtual HydroConservedState computeFlux(
      const HydroPrimitiveState& left_state,
      const HydroPrimitiveState& right_state,
      const HydroFace& face,
      double adiabatic_index) const = 0;

  [[nodiscard]] virtual std::uint64_t fallbackCount() const { return 0; }
};

// Runtime-selectable approximate Riemann solver family.
enum class HydroRiemannSolverKind : std::uint8_t {
  // HLLC (Toro, Spruce & Speares 1994) with Davis wave-speed bounds; falls
  // back to HLLE when the contact/star states are degenerate or inadmissible.
  kHllc = 0,
  // HLLE (Harten-Lax-van Leer with Einfeldt 1988 wave-speed bounds).
  kHlle = 1,
};

[[nodiscard]] std::string_view hydroRiemannSolverKindToString(HydroRiemannSolverKind kind);
[[nodiscard]] HydroRiemannSolverKind hydroRiemannSolverKindFromString(std::string_view name);
[[nodiscard]] std::unique_ptr<HydroRiemannSolver> makeHydroRiemannSolver(HydroRiemannSolverKind kind);

// Passive-scalar contract (audit F3 repair): every lane of
// HydroConservedState, including metal_mass_density_comoving and any lane
// added later through hydro_field_schema.hpp, is evolved through EVERY flux
// branch (supersonic upwind, HLL intermediate state, HLLC star states and the
// HLLC->HLLE fallback). The HLL intermediate flux is evaluated with the
// lane-generic HydroConservedState arithmetic, so a uniform mass fraction X
// always gives F(rho X) = X F(rho) exactly.
//
// HLLE wave speeds (Einfeldt 1988, EOS-agnostic form):
//   d^2 = (sqrt(rho_L) c_L^2 + sqrt(rho_R) c_R^2)/(sqrt(rho_L)+sqrt(rho_R))
//         + eta2 (u_R - u_L)^2,   eta2 = sqrt(rho_L rho_R)/(2 (sqrt(rho_L)+sqrt(rho_R))^2)
//   S_L = min(u_L - c_L, u_R - c_R, u~ - d),  S_R = max(u_L + c_L, u_R + c_R, u~ + d)
// with u~ the sqrt(rho)-weighted normal velocity. Taking the min/max over
// the Davis and Einfeldt estimates keeps the fallback at least as
// dissipative as either bound.
class HlleRiemannSolver final : public HydroRiemannSolver {
 public:
  [[nodiscard]] HydroConservedState computeFlux(
      const HydroPrimitiveState& left_state,
      const HydroPrimitiveState& right_state,
      const HydroFace& face,
      double adiabatic_index) const override;
};

// HLLC is the default approximate Riemann solver with HLLE fallback for positivity/speed degeneracy.
class HllcRiemannSolver final : public HydroRiemannSolver {
 public:
  [[nodiscard]] HydroConservedState computeFlux(
      const HydroPrimitiveState& left_state,
      const HydroPrimitiveState& right_state,
      const HydroFace& face,
      double adiabatic_index) const override;

  [[nodiscard]] std::uint64_t fallbackCount() const override;

 private:
  mutable std::uint64_t m_fallback_count = 0;
};

}  // namespace cosmosim::hydro
