#include "cosmosim/hydro/hydro_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {
namespace {

constexpr double k_small = 1.0e-14;

// Component layout of the batch-local primitive vector used by the
// reconstruction kernels. Kept as a fixed-size array so the per-face kernels
// never allocate and the compiler can unroll the component loops.
constexpr std::size_t k_rho = 0U;
constexpr std::size_t k_vx = 1U;
constexpr std::size_t k_vy = 2U;
constexpr std::size_t k_vz = 3U;
constexpr std::size_t k_pressure = 4U;
constexpr std::size_t k_internal_energy = 5U;
constexpr std::size_t k_signal_speed_squared = 6U;
constexpr std::size_t k_metallicity = 7U;
constexpr std::size_t k_component_count = 8U;
using PrimitiveVector = std::array<double, k_component_count>;

[[nodiscard]] PrimitiveVector toVector(const HydroPrimitiveState& state) {
  return PrimitiveVector{
      state.rho_comoving,
      state.vel_x_peculiar,
      state.vel_y_peculiar,
      state.vel_z_peculiar,
      state.pressure_comoving,
      state.specific_internal_energy_code,
      state.signal_speed_squared_code,
      state.metallicity_mass_fraction};
}

[[nodiscard]] HydroPrimitiveState fromVector(const PrimitiveVector& v, bool uses_effective_ism) {
  HydroPrimitiveState state;
  state.rho_comoving = v[k_rho];
  state.vel_x_peculiar = v[k_vx];
  state.vel_y_peculiar = v[k_vy];
  state.vel_z_peculiar = v[k_vz];
  state.pressure_comoving = v[k_pressure];
  state.specific_internal_energy_code = v[k_internal_energy];
  state.signal_speed_squared_code = v[k_signal_speed_squared];
  state.uses_effective_ism = uses_effective_ism;
  state.metallicity_mass_fraction = v[k_metallicity];
  return state;
}

[[nodiscard]] bool isAdmissible(const HydroPrimitiveState& state, double rho_floor, double pressure_floor) {
  return std::isfinite(state.rho_comoving) && std::isfinite(state.pressure_comoving) &&
      std::isfinite(state.vel_x_peculiar) && std::isfinite(state.vel_y_peculiar) &&
      std::isfinite(state.vel_z_peculiar) && std::isfinite(state.specific_internal_energy_code) &&
      std::isfinite(state.signal_speed_squared_code) &&
      state.rho_comoving >= rho_floor && state.pressure_comoving >= pressure_floor &&
      state.specific_internal_energy_code > 0.0 && state.signal_speed_squared_code > 0.0;
}

// Last-resort clamp, used only when the cell average itself violates the
// floors (the face states must be admissible inputs to the Riemann solver).
// Returns true when a value changed.
[[nodiscard]] bool enforcePrimitiveFloors(HydroPrimitiveState& state, double rho_floor, double pressure_floor) {
  bool changed = false;
  if (!(state.rho_comoving >= rho_floor)) {
    state.rho_comoving = rho_floor;
    changed = true;
  }
  if (!(state.pressure_comoving >= pressure_floor)) {
    state.pressure_comoving = pressure_floor;
    changed = true;
  }
  if (!(state.specific_internal_energy_code > 0.0)) {
    state.specific_internal_energy_code = pressure_floor /
        std::max(state.rho_comoving, rho_floor);
    changed = true;
  }
  if (!(state.signal_speed_squared_code > 0.0)) {
    state.signal_speed_squared_code = state.pressure_comoving /
        std::max(state.rho_comoving, rho_floor);
    changed = true;
  }
  const double metallicity_before = state.metallicity_mass_fraction;
  state.metallicity_mass_fraction = std::clamp(state.metallicity_mass_fraction, 0.0, 1.0);
  changed = changed || state.metallicity_mass_fraction != metallicity_before;
  return changed;
}

[[nodiscard]] double minmodLimiter(double a, double b) {
  if (a * b <= 0.0) {
    return 0.0;
  }
  return (std::abs(a) < std::abs(b)) ? a : b;
}

[[nodiscard]] std::size_t axisIndex(HydroFaceAxis axis) {
  switch (axis) {
    case HydroFaceAxis::kX:
      return 0U;
    case HydroFaceAxis::kY:
      return 1U;
    case HydroFaceAxis::kZ:
      return 2U;
  }
  return 0U;
}

[[nodiscard]] double dtOverCellWidth(const HydroReconstructionPolicy& policy, HydroFaceAxis axis) {
  const double axis_value = policy.dt_over_cell_width_code[axisIndex(axis)];
  return axis_value > 0.0 ? axis_value : policy.dt_over_dx_code;
}

void validatePolicy(const HydroReconstructionPolicy& policy, const char* who) {
  if (!(policy.rho_floor > 0.0) || !(policy.pressure_floor > 0.0)) {
    throw std::invalid_argument(std::string(who) + " reconstruction requires positive floors");
  }
  if (!(policy.dt_over_dx_code >= 0.0) ||
      !(policy.dt_over_cell_width_code[0] >= 0.0) ||
      !(policy.dt_over_cell_width_code[1] >= 0.0) ||
      !(policy.dt_over_cell_width_code[2] >= 0.0)) {
    throw std::invalid_argument(std::string(who) + " reconstruction requires non-negative CFL policy values");
  }
  if (!(policy.adiabatic_index > 1.0)) {
    throw std::invalid_argument(std::string(who) + " reconstruction requires adiabatic_index > 1");
  }
}

// Frame of the face normal. Stencil order follows the normal, so the normal
// velocity is u . n and the momentum response to a normal pressure gradient
// is -n (dP/ds)/rho.
struct NormalFrame {
  double nx = 1.0;
  double ny = 0.0;
  double nz = 0.0;
};

[[nodiscard]] NormalFrame normalFrame(const HydroFace& face) {
  return NormalFrame{face.normal_x, face.normal_y, face.normal_z};
}

[[nodiscard]] double dotNormal(const PrimitiveVector& v, const NormalFrame& n) {
  return v[k_vx] * n.nx + v[k_vy] * n.ny + v[k_vz] * n.nz;
}

// Forward primitive time derivative times (lambda/2), evaluated at the cell
// centre W_c with the along-normal slope dW (see header for the equations).
[[nodiscard]] PrimitiveVector hancockHalfStepIncrement(
    const PrimitiveVector& center,
    const PrimitiveVector& slope,
    const NormalFrame& n,
    double half_lambda,
    double fallback_gamma) {
  const double rho = std::max(center[k_rho], k_small);
  const double pressure = std::max(center[k_pressure], k_small);
  const double u_n = dotNormal(center, n);
  const double du_n = dotNormal(slope, n);
  const double gamma1 = center[k_signal_speed_squared] > 0.0
      ? std::max(rho * center[k_signal_speed_squared] / pressure, 1.0)
      : fallback_gamma;
  const double inv_rho = 1.0 / rho;

  PrimitiveVector d{};
  d[k_rho] = -(u_n * slope[k_rho] + rho * du_n);
  d[k_vx] = -u_n * slope[k_vx] - n.nx * slope[k_pressure] * inv_rho;
  d[k_vy] = -u_n * slope[k_vy] - n.ny * slope[k_pressure] * inv_rho;
  d[k_vz] = -u_n * slope[k_vz] - n.nz * slope[k_pressure] * inv_rho;
  d[k_pressure] = -(u_n * slope[k_pressure] + gamma1 * pressure * du_n);
  d[k_internal_energy] = -(u_n * slope[k_internal_energy] + pressure * inv_rho * du_n);
  d[k_signal_speed_squared] = 0.0;  // rebuilt from Gamma1 P / rho below
  d[k_metallicity] = -u_n * slope[k_metallicity];
  for (double& value : d) {
    value *= half_lambda;
  }
  return d;
}

// Rebuild the thermodynamically dependent lanes of a face state from the
// reconstructed (rho, P): c^2 = Gamma1 P / rho with Gamma1 frozen at the
// cell centre, and e = P/((gamma-1) rho) for ideal-gas states so the face
// state handed to the Riemann solver satisfies its own EOS.
void closeFaceThermodynamics(
    PrimitiveVector& face,
    const PrimitiveVector& center,
    bool uses_effective_ism,
    double gamma) {
  const double rho_c = std::max(center[k_rho], k_small);
  const double p_c = std::max(center[k_pressure], k_small);
  const double gamma1 = center[k_signal_speed_squared] > 0.0
      ? std::max(rho_c * center[k_signal_speed_squared] / p_c, 1.0)
      : gamma;
  if (face[k_rho] > 0.0 && face[k_pressure] > 0.0) {
    face[k_signal_speed_squared] = gamma1 * face[k_pressure] / face[k_rho];
    if (!uses_effective_ism) {
      face[k_internal_energy] = face[k_pressure] / ((gamma - 1.0) * face[k_rho]);
    }
  }
  face[k_metallicity] = std::clamp(face[k_metallicity], 0.0, 1.0);
}

struct SideReconstruction {
  PrimitiveVector center{};
  PrimitiveVector face_spatial{};  // spatial interface value (no prediction)
  PrimitiveVector slope{};         // mean along-normal gradient used by Hancock
  bool uses_effective_ism = false;
};

// Select the admissible face state: predicted -> spatial -> first order ->
// clamp. Returns the number of fallback steps taken (0 = predicted state).
[[nodiscard]] int finalizeSide(
    const SideReconstruction& side,
    const NormalFrame& n,
    double half_lambda,
    bool predict,
    const HydroReconstructionPolicy& policy,
    HydroPrimitiveState& out) {
  if (predict && half_lambda > 0.0) {
    PrimitiveVector predicted = side.face_spatial;
    const PrimitiveVector d =
        hancockHalfStepIncrement(side.center, side.slope, n, half_lambda, policy.adiabatic_index);
    for (std::size_t c = 0; c < k_component_count; ++c) {
      predicted[c] += d[c];
    }
    closeFaceThermodynamics(predicted, side.center, side.uses_effective_ism, policy.adiabatic_index);
    out = fromVector(predicted, side.uses_effective_ism);
    if (isAdmissible(out, policy.rho_floor, policy.pressure_floor)) {
      return 0;
    }
  }
  PrimitiveVector spatial = side.face_spatial;
  closeFaceThermodynamics(spatial, side.center, side.uses_effective_ism, policy.adiabatic_index);
  out = fromVector(spatial, side.uses_effective_ism);
  if (isAdmissible(out, policy.rho_floor, policy.pressure_floor)) {
    return (predict && half_lambda > 0.0) ? 1 : 0;
  }
  out = fromVector(side.center, side.uses_effective_ism);
  if (isAdmissible(out, policy.rho_floor, policy.pressure_floor)) {
    return 2;
  }
  (void)enforcePrimitiveFloors(out, policy.rho_floor, policy.pressure_floor);
  return 3;
}

[[nodiscard]] double limitedSlope(
    HydroSlopeLimiter limiter, double qm, double q, double qp, std::uint64_t& clip_counter) {
  const double dm = q - qm;
  const double dp = qp - q;
  const double slope = applyHydroSlopeLimiter(limiter, dm, dp);
  if (std::abs(slope) < std::abs(0.5 * (dm + dp))) {
    ++clip_counter;
  }
  return slope;
}

// Limited piecewise-linear side. `upstream` is the stencil cell on the far
// side of `center` from the face, `downstream` the cell across the face.
// sign = +1 for the owner (face at s=+1/2), -1 for the neighbor (s=-1/2).
[[nodiscard]] SideReconstruction plmSide(
    const HydroPrimitiveState& upstream_state,
    const HydroPrimitiveState& center_state,
    const HydroPrimitiveState& downstream_state,
    double sign,
    HydroSlopeLimiter limiter,
    std::uint64_t& clip_counter) {
  // Along-normal ordering: for the owner, (upstream, center, downstream) =
  // (owner_minus, owner, neighbor); for the neighbor they are passed as
  // (neighbor_plus, neighbor, owner) and the slope sign is restored below.
  const PrimitiveVector qm = toVector(upstream_state);
  const PrimitiveVector q = toVector(center_state);
  const PrimitiveVector qp = toVector(downstream_state);
  SideReconstruction side;
  side.center = q;
  side.uses_effective_ism = center_state.uses_effective_ism;
  for (std::size_t c = 0; c < k_component_count; ++c) {
    // slope_toward_face is d q per cell in the direction of the face.
    const double slope_toward_face = limitedSlope(limiter, qm[c], q[c], qp[c], clip_counter);
    side.face_spatial[c] = q[c] + 0.5 * slope_toward_face;
    // Convert to the along-normal slope used by the predictor.
    side.slope[c] = sign * slope_toward_face;
  }
  return side;
}

[[nodiscard]] SideReconstruction constantSide(const HydroPrimitiveState& center_state) {
  SideReconstruction side;
  side.center = toVector(center_state);
  side.face_spatial = side.center;
  side.slope.fill(0.0);
  side.uses_effective_ism = center_state.uses_effective_ism;
  return side;
}

}  // namespace

std::string_view hydroSlopeLimiterToString(HydroSlopeLimiter limiter) {
  switch (limiter) {
    case HydroSlopeLimiter::kMinmod:
      return "minmod";
    case HydroSlopeLimiter::kMonotonizedCentral:
      return "mc";
    case HydroSlopeLimiter::kVanLeer:
      return "van_leer";
  }
  return "unknown";
}

HydroSlopeLimiter hydroSlopeLimiterFromString(std::string_view name) {
  if (name == "minmod") {
    return HydroSlopeLimiter::kMinmod;
  }
  if (name == "mc" || name == "monotonized_central") {
    return HydroSlopeLimiter::kMonotonizedCentral;
  }
  if (name == "van_leer") {
    return HydroSlopeLimiter::kVanLeer;
  }
  throw std::invalid_argument("Unsupported hydro slope limiter: " + std::string(name));
}

std::string_view hydroReconstructionSchemeToString(HydroReconstructionScheme scheme) {
  switch (scheme) {
    case HydroReconstructionScheme::kPiecewiseConstant:
      return "piecewise_constant";
    case HydroReconstructionScheme::kPiecewiseLinear:
      return "piecewise_linear";
    case HydroReconstructionScheme::kPiecewiseParabolic:
      return "piecewise_parabolic";
  }
  return "unknown";
}

HydroReconstructionScheme hydroReconstructionSchemeFromString(std::string_view name) {
  if (name == "piecewise_constant" || name == "pcm") {
    return HydroReconstructionScheme::kPiecewiseConstant;
  }
  if (name == "piecewise_linear" || name == "plm") {
    return HydroReconstructionScheme::kPiecewiseLinear;
  }
  if (name == "piecewise_parabolic" || name == "ppm") {
    return HydroReconstructionScheme::kPiecewiseParabolic;
  }
  throw std::invalid_argument("Unsupported hydro reconstruction scheme: " + std::string(name));
}

double applyHydroSlopeLimiter(HydroSlopeLimiter limiter, double delta_minus, double delta_plus) {
  switch (limiter) {
    case HydroSlopeLimiter::kMinmod:
      return minmodLimiter(delta_minus, delta_plus);
    case HydroSlopeLimiter::kMonotonizedCentral: {
      const double mm = minmodLimiter(delta_minus, delta_plus);
      return minmodLimiter(0.5 * (delta_minus + delta_plus), 2.0 * mm);
    }
    case HydroSlopeLimiter::kVanLeer: {
      const double product = delta_minus * delta_plus;
      if (product <= 0.0) {
        return 0.0;
      }
      return 2.0 * product / (delta_minus + delta_plus);
    }
  }
  return 0.0;
}

std::array<double, 3> hydroTransportCourantFactors(
    double dt_code,
    double scale_factor,
    double cell_width_x_comoving,
    double cell_width_y_comoving,
    double cell_width_z_comoving) {
  if (!(dt_code >= 0.0) || !std::isfinite(dt_code)) {
    throw std::invalid_argument("hydroTransportCourantFactors requires finite dt_code >= 0");
  }
  if (!(scale_factor > 0.0) || !std::isfinite(scale_factor)) {
    throw std::invalid_argument("hydroTransportCourantFactors requires finite scale_factor > 0");
  }
  if (!(cell_width_x_comoving > 0.0) || !(cell_width_y_comoving > 0.0) ||
      !(cell_width_z_comoving > 0.0)) {
    throw std::invalid_argument("hydroTransportCourantFactors requires positive cell widths");
  }
  const double dt_over_a = dt_code / scale_factor;
  return {dt_over_a / cell_width_x_comoving,
          dt_over_a / cell_width_y_comoving,
          dt_over_a / cell_width_z_comoving};
}

bool HydroReconstruction::reconstructFaceFromCache(
    const HydroPrimitiveCacheSoa& primitive_cache,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  (void)primitive_cache;
  (void)face;
  (void)left_state;
  (void)right_state;
  return false;
}

bool HydroReconstruction::reconstructFaceFromPrimitiveStates(
    const HydroPrimitiveState& left_minus,
    const HydroPrimitiveState& left_cell,
    const HydroPrimitiveState& right_cell,
    const HydroPrimitiveState& right_plus,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  (void)left_minus;
  (void)left_cell;
  (void)right_cell;
  (void)right_plus;
  (void)face;
  (void)left_state;
  (void)right_state;
  return false;
}

bool HydroReconstruction::reconstructFaceFromStencil(
    const HydroFaceStencilStates& stencil,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  if (stencil.owner == nullptr) {
    throw std::invalid_argument("hydro reconstruction stencil requires the owner cell");
  }
  const HydroPrimitiveState& owner = *stencil.owner;
  const HydroPrimitiveState& neighbor = stencil.neighbor != nullptr ? *stencil.neighbor : owner;
  const HydroPrimitiveState& owner_minus = stencil.owner_minus != nullptr ? *stencil.owner_minus : owner;
  const HydroPrimitiveState& neighbor_plus =
      stencil.neighbor_plus != nullptr ? *stencil.neighbor_plus : neighbor;
  return reconstructFaceFromPrimitiveStates(
      owner_minus, owner, neighbor, neighbor_plus, face, left_state, right_state);
}

void PiecewiseConstantReconstruction::reconstructFace(
    const HydroConservedStateSoa& conserved,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state,
    double adiabatic_index) const {
  left_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.owner_cell), adiabatic_index);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return;
  }
  right_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.neighbor_cell), adiabatic_index);
}

bool PiecewiseConstantReconstruction::reconstructFaceFromCache(
    const HydroPrimitiveCacheSoa& primitive_cache,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  left_state = primitive_cache.loadCell(face.owner_cell);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return true;
  }
  right_state = primitive_cache.loadCell(face.neighbor_cell);
  return true;
}

bool PiecewiseConstantReconstruction::reconstructFaceFromPrimitiveStates(
    const HydroPrimitiveState& left_minus,
    const HydroPrimitiveState& left_cell,
    const HydroPrimitiveState& right_cell,
    const HydroPrimitiveState& right_plus,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  (void)left_minus;
  (void)right_plus;
  (void)face;
  left_state = left_cell;
  right_state = right_cell;
  return true;
}

MusclHancockReconstruction::MusclHancockReconstruction(HydroReconstructionPolicy policy) : m_policy(policy) {
  validatePolicy(m_policy, "MUSCL-Hancock");
}

bool MusclHancockReconstruction::reconstructFaceFromCache(
    const HydroPrimitiveCacheSoa& primitive_cache,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  const HydroPrimitiveState left_cell = primitive_cache.loadCell(face.owner_cell);
  if (face.neighbor_cell == k_invalid_cell_index) {
    left_state = left_cell;
    right_state = left_cell;
    return true;
  }
  const HydroPrimitiveState right_cell = primitive_cache.loadCell(face.neighbor_cell);
  const HydroPrimitiveState left_minus = face.owner_minus_cell == k_invalid_cell_index
      ? left_cell
      : primitive_cache.loadCell(face.owner_minus_cell);
  const HydroPrimitiveState right_plus = face.neighbor_plus_cell == k_invalid_cell_index
      ? right_cell
      : primitive_cache.loadCell(face.neighbor_plus_cell);

  return reconstructFaceFromPrimitiveStates(
      left_minus,
      left_cell,
      right_cell,
      right_plus,
      face,
      left_state,
      right_state);
}

bool MusclHancockReconstruction::reconstructFaceFromPrimitiveStates(
    const HydroPrimitiveState& left_minus,
    const HydroPrimitiveState& left_cell,
    const HydroPrimitiveState& right_cell,
    const HydroPrimitiveState& right_plus,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  const NormalFrame n = normalFrame(face);
  // Owner: slope along +s from (owner_minus, owner, neighbor).
  const SideReconstruction owner_side =
      plmSide(left_minus, left_cell, right_cell, +1.0, m_policy.limiter, m_limiter_clip_count);
  // Neighbor: slope toward the face is along -s, from (neighbor_plus,
  // neighbor, owner); plmSide restores the along-normal sign.
  const SideReconstruction neighbor_side =
      plmSide(right_plus, right_cell, left_cell, -1.0, m_policy.limiter, m_limiter_clip_count);

  const double half_lambda = 0.5 * dtOverCellWidth(m_policy, face.axis);
  const bool predict = m_policy.enable_muscl_hancock_predictor;
  if (finalizeSide(owner_side, n, half_lambda, predict, m_policy, left_state) != 0) {
    ++m_positivity_fallback_count;
  }
  if (finalizeSide(neighbor_side, n, half_lambda, predict, m_policy, right_state) != 0) {
    ++m_positivity_fallback_count;
  }
  return true;
}

void MusclHancockReconstruction::reconstructFace(
    const HydroConservedStateSoa& conserved,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state,
    double adiabatic_index) const {
  left_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.owner_cell), adiabatic_index);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return;
  }
  right_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.neighbor_cell), adiabatic_index);
}

HydroReconstructionPolicy MusclHancockReconstruction::policy() const { return m_policy; }
std::uint64_t MusclHancockReconstruction::limiterClipCount() const { return m_limiter_clip_count; }
std::uint64_t MusclHancockReconstruction::positivityFallbackCount() const { return m_positivity_fallback_count; }

// ----------------------------------------------------------------------------
// PPM
// ----------------------------------------------------------------------------

HydroParabolaEdges ppmCellEdges(const std::array<double, 5>& q) {
  // q = {q_{i-2}, q_{i-1}, q_i, q_{i+1}, q_{i+2}}
  // MC-limited cell slopes (CW84 eq. 1.7-1.8).
  const auto mc_slope = [](double qm, double q0, double qp) {
    const double dm = q0 - qm;
    const double dp = qp - q0;
    if (dm * dp <= 0.0) {
      return 0.0;
    }
    const double central = 0.5 * (qp - qm);
    const double bound = 2.0 * std::min(std::abs(dm), std::abs(dp));
    return std::copysign(std::min(std::abs(central), bound), central);
  };
  const double dq_im1 = mc_slope(q[0], q[1], q[2]);
  const double dq_i = mc_slope(q[1], q[2], q[3]);
  const double dq_ip1 = mc_slope(q[2], q[3], q[4]);
  // CW84 eq. 1.6: interface values.
  double left = q[1] + 0.5 * (q[2] - q[1]) - (dq_i - dq_im1) / 6.0;
  double right = q[2] + 0.5 * (q[3] - q[2]) - (dq_ip1 - dq_i) / 6.0;
  // Keep each interface value between its neighbouring averages.
  left = std::clamp(left, std::min(q[1], q[2]), std::max(q[1], q[2]));
  right = std::clamp(right, std::min(q[2], q[3]), std::max(q[2], q[3]));

  // CW84 eq. 1.10 monotonicity constraint.
  const double q0 = q[2];
  if ((right - q0) * (q0 - left) <= 0.0) {
    return HydroParabolaEdges{q0, q0};
  }
  const double delta = right - left;
  const double six_mean_offset = 6.0 * delta * (q0 - 0.5 * (left + right));
  const double delta_squared = delta * delta;
  if (six_mean_offset > delta_squared) {
    left = 3.0 * q0 - 2.0 * right;
  } else if (-delta_squared > six_mean_offset) {
    right = 3.0 * q0 - 2.0 * left;
  }
  return HydroParabolaEdges{left, right};
}

namespace {

// CW84 flattening coefficient for cell i along the face normal.
// p/u_n arrays hold cells i-2..i+2.
[[nodiscard]] double ppmFlatteningCoefficient(
    const std::array<double, 5>& pressure, const std::array<double, 5>& normal_velocity) {
  constexpr double k_epsilon = 0.33;
  constexpr double k_omega1 = 0.75;
  constexpr double k_omega2 = 10.0;
  const double dp1 = pressure[3] - pressure[1];
  const double dp2 = pressure[4] - pressure[0];
  const double p_min = std::max(std::min(pressure[3], pressure[1]), k_small);
  const bool shock_candidate = std::abs(dp1) / p_min > k_epsilon &&
      normal_velocity[1] - normal_velocity[3] > 0.0;
  if (!shock_candidate) {
    return 0.0;
  }
  const double ratio = std::abs(dp2) > k_small ? std::abs(dp1 / dp2) : 1.0;
  return std::clamp(k_omega2 * (ratio - k_omega1), 0.0, 1.0);
}

// PPM side reconstruction for the cell `cells[2]` with stencil ordered along
// the face normal (cells[0..4] = i-2..i+2 along +s). face_on_plus_side is true
// for the owner (face at s=+1/2) and false for the neighbor (s=-1/2).
[[nodiscard]] SideReconstruction ppmSide(
    const std::array<const HydroPrimitiveState*, 5>& cells,
    bool face_on_plus_side,
    const NormalFrame& n,
    HydroSlopeLimiter limiter,
    std::uint64_t& clip_counter,
    std::uint64_t& flattening_counter) {
  std::array<PrimitiveVector, 5> v{};
  for (std::size_t k = 0; k < 5U; ++k) {
    v[k] = toVector(*cells[k]);
  }
  std::array<double, 5> pressure{};
  std::array<double, 5> normal_velocity{};
  for (std::size_t k = 0; k < 5U; ++k) {
    pressure[k] = v[k][k_pressure];
    normal_velocity[k] = dotNormal(v[k], n);
  }
  const double chi = ppmFlatteningCoefficient(pressure, normal_velocity);
  if (chi > 0.0) {
    ++flattening_counter;
  }

  SideReconstruction side;
  side.center = v[2];
  side.uses_effective_ism = cells[2]->uses_effective_ism;
  for (std::size_t c = 0; c < k_component_count; ++c) {
    const HydroParabolaEdges edges = ppmCellEdges({v[0][c], v[1][c], v[2][c], v[3][c], v[4][c]});
    const double plm_slope = limitedSlope(limiter, v[1][c], v[2][c], v[3][c], clip_counter);
    const double left = (1.0 - chi) * edges.left + chi * (v[2][c] - 0.5 * plm_slope);
    const double right = (1.0 - chi) * edges.right + chi * (v[2][c] + 0.5 * plm_slope);
    side.face_spatial[c] = face_on_plus_side ? right : left;
    // Mean along-normal gradient of the (flattened) parabola.
    side.slope[c] = right - left;
  }
  return side;
}

}  // namespace

PpmReconstruction::PpmReconstruction(HydroReconstructionPolicy policy) : m_policy(policy) {
  validatePolicy(m_policy, "PPM");
}

bool PpmReconstruction::reconstructFaceFromStencil(
    const HydroFaceStencilStates& stencil,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  if (stencil.owner == nullptr) {
    throw std::invalid_argument("PPM reconstruction stencil requires the owner cell");
  }
  const HydroPrimitiveState& owner = *stencil.owner;
  const HydroPrimitiveState& neighbor = stencil.neighbor != nullptr ? *stencil.neighbor : owner;
  const NormalFrame n = normalFrame(face);

  // Owner cell i parabola needs i-2..i+2 = owner_minus2, owner_minus, owner,
  // neighbor, neighbor_plus.
  SideReconstruction owner_side;
  if (stencil.owner_minus2 != nullptr && stencil.owner_minus != nullptr &&
      stencil.neighbor != nullptr && stencil.neighbor_plus != nullptr) {
    owner_side = ppmSide(
        {stencil.owner_minus2, stencil.owner_minus, &owner, &neighbor, stencil.neighbor_plus},
        true, n, m_policy.limiter, m_limiter_clip_count, m_flattening_activation_count);
  } else {
    ++m_stencil_degradation_count;
    if (stencil.owner_minus != nullptr) {
      owner_side = plmSide(*stencil.owner_minus, owner, neighbor, +1.0,
          m_policy.limiter, m_limiter_clip_count);
    } else {
      owner_side = constantSide(owner);
    }
  }

  // Neighbor cell i+1 parabola needs i-1..i+3 = owner_minus, owner, neighbor,
  // neighbor_plus, neighbor_plus2.
  SideReconstruction neighbor_side;
  if (stencil.neighbor != nullptr && stencil.owner_minus != nullptr &&
      stencil.neighbor_plus != nullptr && stencil.neighbor_plus2 != nullptr) {
    neighbor_side = ppmSide(
        {stencil.owner_minus, &owner, &neighbor, stencil.neighbor_plus, stencil.neighbor_plus2},
        false, n, m_policy.limiter, m_limiter_clip_count, m_flattening_activation_count);
  } else {
    ++m_stencil_degradation_count;
    if (stencil.neighbor_plus != nullptr) {
      neighbor_side = plmSide(*stencil.neighbor_plus, neighbor, owner, -1.0,
          m_policy.limiter, m_limiter_clip_count);
    } else {
      neighbor_side = constantSide(neighbor);
    }
  }

  const double half_lambda = 0.5 * dtOverCellWidth(m_policy, face.axis);
  const bool predict = m_policy.enable_muscl_hancock_predictor;
  if (finalizeSide(owner_side, n, half_lambda, predict, m_policy, left_state) != 0) {
    ++m_positivity_fallback_count;
  }
  if (finalizeSide(neighbor_side, n, half_lambda, predict, m_policy, right_state) != 0) {
    ++m_positivity_fallback_count;
  }
  return true;
}

bool PpmReconstruction::reconstructFaceFromPrimitiveStates(
    const HydroPrimitiveState& left_minus,
    const HydroPrimitiveState& left_cell,
    const HydroPrimitiveState& right_cell,
    const HydroPrimitiveState& right_plus,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  // Four-point callers cannot provide the PPM stencil: both sides reduce to
  // limited piecewise-linear and are counted as degraded.
  return reconstructFaceFromStencil(
      HydroFaceStencilStates{
          .owner_minus2 = nullptr,
          .owner_minus = &left_minus,
          .owner = &left_cell,
          .neighbor = &right_cell,
          .neighbor_plus = &right_plus,
          .neighbor_plus2 = nullptr},
      face, left_state, right_state);
}

bool PpmReconstruction::reconstructFaceFromCache(
    const HydroPrimitiveCacheSoa& primitive_cache,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  const HydroPrimitiveState owner = primitive_cache.loadCell(face.owner_cell);
  if (face.neighbor_cell == k_invalid_cell_index) {
    left_state = owner;
    right_state = owner;
    return true;
  }
  const HydroPrimitiveState neighbor = primitive_cache.loadCell(face.neighbor_cell);
  HydroPrimitiveState owner_minus;
  HydroPrimitiveState owner_minus2;
  HydroPrimitiveState neighbor_plus;
  HydroPrimitiveState neighbor_plus2;
  HydroFaceStencilStates stencil{.owner = &owner, .neighbor = &neighbor};
  if (face.owner_minus_cell != k_invalid_cell_index) {
    owner_minus = primitive_cache.loadCell(face.owner_minus_cell);
    stencil.owner_minus = &owner_minus;
  }
  if (face.owner_minus2_cell != k_invalid_cell_index) {
    owner_minus2 = primitive_cache.loadCell(face.owner_minus2_cell);
    stencil.owner_minus2 = &owner_minus2;
  }
  if (face.neighbor_plus_cell != k_invalid_cell_index) {
    neighbor_plus = primitive_cache.loadCell(face.neighbor_plus_cell);
    stencil.neighbor_plus = &neighbor_plus;
  }
  if (face.neighbor_plus2_cell != k_invalid_cell_index) {
    neighbor_plus2 = primitive_cache.loadCell(face.neighbor_plus2_cell);
    stencil.neighbor_plus2 = &neighbor_plus2;
  }
  return reconstructFaceFromStencil(stencil, face, left_state, right_state);
}

void PpmReconstruction::reconstructFace(
    const HydroConservedStateSoa& conserved,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state,
    double adiabatic_index) const {
  left_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.owner_cell), adiabatic_index);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return;
  }
  right_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.neighbor_cell), adiabatic_index);
}

HydroReconstructionPolicy PpmReconstruction::policy() const { return m_policy; }
std::uint64_t PpmReconstruction::limiterClipCount() const { return m_limiter_clip_count; }
std::uint64_t PpmReconstruction::positivityFallbackCount() const { return m_positivity_fallback_count; }
std::uint64_t PpmReconstruction::stencilDegradationCount() const { return m_stencil_degradation_count; }
std::uint64_t PpmReconstruction::flatteningActivationCount() const { return m_flattening_activation_count; }

}  // namespace cosmosim::hydro
