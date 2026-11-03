#pragma once

#include <cstdint>
#include <array>
#include <string_view>

namespace cosmosim::hydro {

struct HydroFace;
struct HydroPrimitiveState;
class HydroConservedStateSoa;
class HydroPrimitiveCacheSoa;

enum class HydroSlopeLimiter {
  kMinmod,
  kMonotonizedCentral,
  kVanLeer,
};

// Spatial reconstruction family. This is orthogonal to the time-integration
// family (HydroIntegrationScheme in hydro_core_solver.hpp) and to the Riemann
// solver. Every reconstruction supplies spatial left/right interface values;
// the single forward half-step (Hancock) predictor is applied once, inside the
// reconstruction, using the cell-centred state and the reconstructed mean
// gradient. No reconstruction stacks a second half-time prediction.
enum class HydroReconstructionScheme : std::uint8_t {
  // First-order Godunov reference: interface values are the cell averages.
  kPiecewiseConstant = 0,
  // Limited piecewise-linear (MUSCL) reconstruction, stencil half-width 2.
  kPiecewiseLinear = 1,
  // Colella & Woodward (1984) piecewise-parabolic reconstruction,
  // stencil half-width 3 (cells i-2..i+3 around face i+1/2).
  kPiecewiseParabolic = 2,
};

[[nodiscard]] std::string_view hydroSlopeLimiterToString(HydroSlopeLimiter limiter);
[[nodiscard]] HydroSlopeLimiter hydroSlopeLimiterFromString(std::string_view name);
[[nodiscard]] double applyHydroSlopeLimiter(HydroSlopeLimiter limiter, double delta_minus, double delta_plus);

[[nodiscard]] std::string_view hydroReconstructionSchemeToString(HydroReconstructionScheme scheme);
[[nodiscard]] HydroReconstructionScheme hydroReconstructionSchemeFromString(std::string_view name);

struct HydroReconstructionPolicy {
  HydroSlopeLimiter limiter = HydroSlopeLimiter::kMonotonizedCentral;
  // Backward-compatible scalar used when dt_over_cell_width_code is left at zero.
  //
  // UNIT CONTRACT (audit F1 follow-up): this is the transport Courant factor
  // of the conservative update, i.e. dt_code / (a * dx_comoving) in comoving
  // cosmological runs (a = HydroUpdateContext::scale_factor, the same epoch
  // used by the flux update dt/(a V)), and dt_code / dx in non-expanding
  // runs. Use hydroTransportCourantFactors() to construct it consistently.
  double dt_over_dx_code = 0.0;
  // Ordered as x, y, z. Values are dt_code / (a * cell_width_axis_comoving).
  std::array<double, 3> dt_over_cell_width_code{0.0, 0.0, 0.0};
  double rho_floor = 1.0e-12;
  double pressure_floor = 1.0e-12;
  bool enable_muscl_hancock_predictor = true;
  double adiabatic_index = 5.0 / 3.0;
};

// Transport Courant factors dt/(a dx_axis) shared by the Hancock predictor,
// the CTU transverse correction and the conservative flux update.
// scale_factor must be the flux-update epoch (normally the interval midpoint);
// non-expanding runs pass 1.
[[nodiscard]] std::array<double, 3> hydroTransportCourantFactors(
    double dt_code,
    double scale_factor,
    double cell_width_x_comoving,
    double cell_width_y_comoving,
    double cell_width_z_comoving);

// Batch-local view of the per-face reconstruction stencil. Cells are ordered
// along the face normal:
//   owner_minus2 -> owner_minus -> owner | face | neighbor -> neighbor_plus -> neighbor_plus2
// A null pointer means that stencil cell is unavailable (patch boundary, AMR
// ghost shell narrower than the scheme, imported ghost). Reconstructions must
// treat missing cells explicitly (documented local order reduction, counted in
// stencilDegradationCount()) instead of fabricating values.
struct HydroFaceStencilStates {
  const HydroPrimitiveState* owner_minus2 = nullptr;
  const HydroPrimitiveState* owner_minus = nullptr;
  const HydroPrimitiveState* owner = nullptr;
  const HydroPrimitiveState* neighbor = nullptr;
  const HydroPrimitiveState* neighbor_plus = nullptr;
  const HydroPrimitiveState* neighbor_plus2 = nullptr;
};

class HydroReconstruction {
 public:
  virtual ~HydroReconstruction() = default;

  // Batch-local primitive reconstruction entry point. The caller owns the
  // four primitive stencil states, so reconstruction never requires a
  // population-sized primitive cache merely to preserve MUSCL-Hancock
  // semantics. Stencil order follows HydroFaceStencilStates.
  [[nodiscard]] virtual bool reconstructFaceFromPrimitiveStates(
      const HydroPrimitiveState& left_minus,
      const HydroPrimitiveState& left_cell,
      const HydroPrimitiveState& right_cell,
      const HydroPrimitiveState& right_plus,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const;

  // Wide-stencil entry point used by the production solver. The default
  // implementation substitutes missing near cells with the adjacent cell
  // (zero one-sided slope, i.e. the historical 4-point behavior) and forwards
  // to reconstructFaceFromPrimitiveStates.
  [[nodiscard]] virtual bool reconstructFaceFromStencil(
      const HydroFaceStencilStates& stencil,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const;

  [[nodiscard]] virtual bool reconstructFaceFromCache(
      const HydroPrimitiveCacheSoa& primitive_cache,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const;

  virtual void reconstructFace(
      const HydroConservedStateSoa& conserved,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state,
      double adiabatic_index) const = 0;

  // Number of cells required on each side of a face (1 = PCM, 2 = PLM, 3 = PPM).
  [[nodiscard]] virtual int requiredStencilHalfWidth() const { return 1; }
  [[nodiscard]] virtual std::uint64_t limiterClipCount() const { return 0; }
  // Interface states whose predicted value was inadmissible (rho or P below
  // floor, non-finite) and were replaced by an admissible lower-order state.
  [[nodiscard]] virtual std::uint64_t positivityFallbackCount() const { return 0; }
  // Face sides reconstructed at reduced order because the stencil was too
  // narrow (patch boundary / ghost shell narrower than the scheme).
  [[nodiscard]] virtual std::uint64_t stencilDegradationCount() const { return 0; }
};

class PiecewiseConstantReconstruction final : public HydroReconstruction {
 public:
  [[nodiscard]] bool reconstructFaceFromPrimitiveStates(
      const HydroPrimitiveState& left_minus,
      const HydroPrimitiveState& left_cell,
      const HydroPrimitiveState& right_cell,
      const HydroPrimitiveState& right_plus,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  [[nodiscard]] bool reconstructFaceFromCache(
      const HydroPrimitiveCacheSoa& primitive_cache,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  void reconstructFace(
      const HydroConservedStateSoa& conserved,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state,
      double adiabatic_index) const override;
};

// MUSCL-Hancock face reconstruction with runtime-selectable slope limiter.
// Slopes are taken from explicit per-face Cartesian stencil metadata.
//
// Predictor (audit F1 repair). With s the coordinate along the face normal
// n (owner -> neighbor), u_n = u . n, Delta q the limited cell slope along s
// and lambda = dt/(a dx) the transport Courant factor, each side's interface
// state is
//   W_face = W_i + (+/-1/2) Delta W + (lambda/2) * D(W_i, Delta W),
// where '+' is the owner (left) side, '-' the neighbor (right) side, and the
// FORWARD primitive time derivative (already carrying the physical sign) is
//   D rho   = -(u_n Delta rho + rho Delta u_n)
//   D u     = -u_n Delta u - n (Delta P / rho)
//   D P     = -(u_n Delta P + Gamma1 P Delta u_n)
//   D e     = -(u_n Delta e + (P/rho) Delta u_n)
//   D Z     = -u_n Delta Z
// evaluated at the cell centre W_i (standard Hancock, Toro 2009 sect. 14.4).
// Gamma1 = rho c^2 / P is the local effective adiabatic exponent supplied by
// the thermodynamic closure (ideal gas: gamma). The signal speed is rebuilt
// as c^2 = Gamma1 P / rho at the predicted state and ideal-gas states keep
// e = P/((gamma-1) rho). For the audit ramp [1,2,3,4], u_n=+1, lambda=0.2,
// the owner face value is 2.5 - 0.1 = 2.4 (the defective version gave 2.6).
// An inadmissible predicted state falls back to the unpredicted limited
// state, then to the first-order cell state, then (only if the cell average
// itself is below the floors) to floor clamping; each fallback is counted.
class MusclHancockReconstruction final : public HydroReconstruction {
 public:
  explicit MusclHancockReconstruction(HydroReconstructionPolicy policy = {});

  [[nodiscard]] bool reconstructFaceFromPrimitiveStates(
      const HydroPrimitiveState& left_minus,
      const HydroPrimitiveState& left_cell,
      const HydroPrimitiveState& right_cell,
      const HydroPrimitiveState& right_plus,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  [[nodiscard]] bool reconstructFaceFromCache(
      const HydroPrimitiveCacheSoa& primitive_cache,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  void reconstructFace(
      const HydroConservedStateSoa& conserved,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state,
      double adiabatic_index) const override;

  [[nodiscard]] HydroReconstructionPolicy policy() const;
  [[nodiscard]] int requiredStencilHalfWidth() const override { return 2; }
  [[nodiscard]] std::uint64_t limiterClipCount() const override;
  [[nodiscard]] std::uint64_t positivityFallbackCount() const override;

 private:
  HydroReconstructionPolicy m_policy;
  mutable std::uint64_t m_limiter_clip_count = 0;
  mutable std::uint64_t m_positivity_fallback_count = 0;
};

// Piecewise-parabolic reconstruction (Colella & Woodward 1984, J. Comput.
// Phys. 54, 174), applied component-wise to primitive variables.
//
//  1. Interface interpolation (CW84 eq. 1.6-1.8) with MC-limited slopes:
//       a_{i+1/2} = q_i + (q_{i+1}-q_i)/2 - (dq_{i+1} - dq_i)/6
//     which reduces to the fourth-order value 7/12(q_i+q_{i+1}) -
//     1/12(q_{i-1}+q_{i+2}) for smooth unlimited data.
//  2. Monotonicity constraint (CW84 eq. 1.10): extrema reset the cell to a
//     constant; overshooting parabolas are steepened toward the far edge.
//  3. Shock flattening (CW84 appendix, simplified to the per-cell
//     coefficient): chi = clamp(10 (|dp1/dp2| - 0.75), 0, 1) when
//     |p_{i+1}-p_{i-1}|/min(p) > 0.33 and u_{i-1} > u_{i+1} (compression);
//     edges are blended toward the limited piecewise-linear edges with
//     weight chi. The upwind-neighbour max(chi_i, chi_{i+-1}) of CW84 needs a
//     seven-cell stencil and is intentionally omitted.
//  4. Contact steepening is NOT applied (documented choice: known to create
//     staircase artifacts in multidimensional flows).
//  5. Time evolution: exactly one forward Hancock half step (same equations
//     as MusclHancockReconstruction) using the parabola's mean gradient
//     Delta q = a_R - a_L at the cell centre. No CW84 characteristic-tracing
//     integral is applied in addition, so PPM never stacks two predictors.
// Stencil: cells i-2..i+3 for face i+1/2. If owner_minus2 / neighbor_plus2
// is unavailable the affected side reduces to limited piecewise-linear
// (counted in stencilDegradationCount()); if owner_minus / neighbor_plus is
// unavailable it reduces to first order on that side.
class PpmReconstruction final : public HydroReconstruction {
 public:
  explicit PpmReconstruction(HydroReconstructionPolicy policy = {});

  [[nodiscard]] bool reconstructFaceFromPrimitiveStates(
      const HydroPrimitiveState& left_minus,
      const HydroPrimitiveState& left_cell,
      const HydroPrimitiveState& right_cell,
      const HydroPrimitiveState& right_plus,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  [[nodiscard]] bool reconstructFaceFromStencil(
      const HydroFaceStencilStates& stencil,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  [[nodiscard]] bool reconstructFaceFromCache(
      const HydroPrimitiveCacheSoa& primitive_cache,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  void reconstructFace(
      const HydroConservedStateSoa& conserved,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state,
      double adiabatic_index) const override;

  [[nodiscard]] HydroReconstructionPolicy policy() const;
  [[nodiscard]] int requiredStencilHalfWidth() const override { return 3; }
  [[nodiscard]] std::uint64_t limiterClipCount() const override;
  [[nodiscard]] std::uint64_t positivityFallbackCount() const override;
  [[nodiscard]] std::uint64_t stencilDegradationCount() const override;
  // Cells for which the flattening coefficient was > 0.
  [[nodiscard]] std::uint64_t flatteningActivationCount() const;

 private:
  HydroReconstructionPolicy m_policy;
  mutable std::uint64_t m_limiter_clip_count = 0;
  mutable std::uint64_t m_positivity_fallback_count = 0;
  mutable std::uint64_t m_stencil_degradation_count = 0;
  mutable std::uint64_t m_flattening_activation_count = 0;
};

// Per-cell parabola edge values for one scalar (exposed for unit testing of
// the CW84 interpolation and monotonicity constraint).
struct HydroParabolaEdges {
  double left = 0.0;
  double right = 0.0;
};
// q[0..4] = q_{i-2}..q_{i+2}; returns the monotonized edge values of cell i
// (before flattening).
[[nodiscard]] HydroParabolaEdges ppmCellEdges(const std::array<double, 5>& q);

}  // namespace cosmosim::hydro
