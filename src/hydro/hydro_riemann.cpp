#include "cosmosim/hydro/hydro_riemann.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {
namespace {

constexpr double k_small = 1.0e-14;

[[nodiscard]] double dot3(double ax, double ay, double az, double bx, double by, double bz) {
  return ax * bx + ay * by + az * bz;
}

[[nodiscard]] double soundSpeed(const HydroPrimitiveState& primitive, double adiabatic_index) {
  const double c_squared = primitive.signal_speed_squared_code > 0.0
      ? primitive.signal_speed_squared_code
      : adiabatic_index * std::max(primitive.pressure_comoving, k_small) /
          std::max(primitive.rho_comoving, k_small);
  return std::sqrt(std::max(c_squared, k_small));
}

[[nodiscard]] HydroConservedState eulerPhysicalFlux(
    const HydroPrimitiveState& primitive,
    const HydroFace& face,
    double adiabatic_index) {
  const double mass_density = primitive.rho_comoving;
  const double vel_n = dot3(
      primitive.vel_x_peculiar,
      primitive.vel_y_peculiar,
      primitive.vel_z_peculiar,
      face.normal_x,
      face.normal_y,
      face.normal_z);

  HydroConservedState flux;
  flux.mass_density_comoving = mass_density * vel_n;
  flux.momentum_density_x_comoving = mass_density * primitive.vel_x_peculiar * vel_n + primitive.pressure_comoving * face.normal_x;
  flux.momentum_density_y_comoving = mass_density * primitive.vel_y_peculiar * vel_n + primitive.pressure_comoving * face.normal_y;
  flux.momentum_density_z_comoving = mass_density * primitive.vel_z_peculiar * vel_n + primitive.pressure_comoving * face.normal_z;

  const HydroConservedState conserved = HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index);
  flux.total_energy_density_comoving = (conserved.total_energy_density_comoving + primitive.pressure_comoving) * vel_n;
  flux.metal_mass_density_comoving =
      primitive.rho_comoving * std::clamp(primitive.metallicity_mass_fraction, 0.0, 1.0) * vel_n;
  return flux;
}

struct HllWaveSpeeds {
  double s_left = 0.0;
  double s_right = 0.0;
};

[[nodiscard]] double normalVelocity(const HydroPrimitiveState& state, const HydroFace& face) {
  return dot3(state.vel_x_peculiar, state.vel_y_peculiar, state.vel_z_peculiar,
      face.normal_x, face.normal_y, face.normal_z);
}

// Einfeldt (1988) bounds combined with the Davis direct estimates. The
// averaged signal speed uses only rho, u_n and c^2, so it is valid for the
// effective-ISM closure as well as for an ideal gas.
[[nodiscard]] HllWaveSpeeds einfeldtWaveSpeeds(
    const HydroPrimitiveState& left_state,
    const HydroPrimitiveState& right_state,
    const HydroFace& face,
    double adiabatic_index) {
  const double u_left = normalVelocity(left_state, face);
  const double u_right = normalVelocity(right_state, face);
  const double c_left = soundSpeed(left_state, adiabatic_index);
  const double c_right = soundSpeed(right_state, adiabatic_index);
  const double sqrt_rho_left = std::sqrt(std::max(left_state.rho_comoving, k_small));
  const double sqrt_rho_right = std::sqrt(std::max(right_state.rho_comoving, k_small));
  const double weight_sum = sqrt_rho_left + sqrt_rho_right;
  const double u_average = (sqrt_rho_left * u_left + sqrt_rho_right * u_right) / weight_sum;
  const double eta2 = 0.5 * sqrt_rho_left * sqrt_rho_right / (weight_sum * weight_sum);
  const double du = u_right - u_left;
  const double d_squared =
      (sqrt_rho_left * c_left * c_left + sqrt_rho_right * c_right * c_right) / weight_sum +
      eta2 * du * du;
  const double d = std::sqrt(std::max(d_squared, k_small));
  return HllWaveSpeeds{
      .s_left = std::min({u_left - c_left, u_right - c_right, u_average - d}),
      .s_right = std::max({u_left + c_left, u_right + c_right, u_average + d})};
}

// Lane-generic HLL flux. All conserved lanes, including passive scalars,
// pass through the same formula (audit F3 repair).
[[nodiscard]] HydroConservedState hllFlux(
    const HydroPrimitiveState& left_state,
    const HydroPrimitiveState& right_state,
    const HydroFace& face,
    double adiabatic_index) {
  const HllWaveSpeeds speeds = einfeldtWaveSpeeds(left_state, right_state, face, adiabatic_index);
  const HydroConservedState f_left = eulerPhysicalFlux(left_state, face, adiabatic_index);
  if (speeds.s_left >= 0.0) {
    return f_left;
  }
  const HydroConservedState f_right = eulerPhysicalFlux(right_state, face, adiabatic_index);
  if (speeds.s_right <= 0.0) {
    return f_right;
  }
  const HydroConservedState u_left = HydroCoreSolver::conservedFromPrimitive(left_state, adiabatic_index);
  const HydroConservedState u_right = HydroCoreSolver::conservedFromPrimitive(right_state, adiabatic_index);
  const double inv_speed_delta = 1.0 / std::max(speeds.s_right - speeds.s_left, k_small);
  return inv_speed_delta *
      (speeds.s_right * f_left - speeds.s_left * f_right +
       (speeds.s_left * speeds.s_right) * (u_right - u_left));
}

}  // namespace

HydroConservedState HlleRiemannSolver::computeFlux(
    const HydroPrimitiveState& left_state,
    const HydroPrimitiveState& right_state,
    const HydroFace& face,
    double adiabatic_index) const {
  return hllFlux(left_state, right_state, face, adiabatic_index);
}

HydroConservedState HllcRiemannSolver::computeFlux(
    const HydroPrimitiveState& left_state,
    const HydroPrimitiveState& right_state,
    const HydroFace& face,
    double adiabatic_index) const {
  const double un_left = dot3(left_state.vel_x_peculiar, left_state.vel_y_peculiar, left_state.vel_z_peculiar,
      face.normal_x, face.normal_y, face.normal_z);
  const double un_right = dot3(right_state.vel_x_peculiar, right_state.vel_y_peculiar, right_state.vel_z_peculiar,
      face.normal_x, face.normal_y, face.normal_z);

  const double c_left = soundSpeed(left_state, adiabatic_index);
  const double c_right = soundSpeed(right_state, adiabatic_index);

  const double s_left = std::min(un_left - c_left, un_right - c_right);
  const double s_right = std::max(un_left + c_left, un_right + c_right);

  if (s_left >= 0.0) {
    return eulerPhysicalFlux(left_state, face, adiabatic_index);
  }
  if (s_right <= 0.0) {
    return eulerPhysicalFlux(right_state, face, adiabatic_index);
  }

  const HydroConservedState u_left = HydroCoreSolver::conservedFromPrimitive(left_state, adiabatic_index);
  const HydroConservedState u_right = HydroCoreSolver::conservedFromPrimitive(right_state, adiabatic_index);
  const HydroConservedState f_left = eulerPhysicalFlux(left_state, face, adiabatic_index);
  const HydroConservedState f_right = eulerPhysicalFlux(right_state, face, adiabatic_index);

  const double numerator = right_state.pressure_comoving - left_state.pressure_comoving +
      left_state.rho_comoving * un_left * (s_left - un_left) -
      right_state.rho_comoving * un_right * (s_right - un_right);
  const double denominator = left_state.rho_comoving * (s_left - un_left) -
      right_state.rho_comoving * (s_right - un_right);
  if (std::abs(denominator) < k_small) {
    ++m_fallback_count;
    return hllFlux(left_state, right_state, face, adiabatic_index);
  }
  const double s_star = numerator / denominator;

  const double p_star_left = left_state.pressure_comoving + left_state.rho_comoving * (s_left - un_left) * (s_star - un_left);
  const double p_star_right = right_state.pressure_comoving + right_state.rho_comoving * (s_right - un_right) * (s_star - un_right);
  const double p_star = 0.5 * (p_star_left + p_star_right);
  if (p_star <= 0.0) {
    ++m_fallback_count;
    return hllFlux(left_state, right_state, face, adiabatic_index);
  }

  auto starState = [&](const HydroPrimitiveState& prim, const HydroConservedState& cons, double s_k, double u_k_n) {
    const double denom = s_k - s_star;
    if (std::abs(denom) < k_small) {
      return HydroConservedState{};
    }
    const double rho_star = prim.rho_comoving * (s_k - u_k_n) / denom;
    if (!std::isfinite(rho_star) || rho_star <= 0.0) {
      return HydroConservedState{};
    }
    HydroConservedState star;
    star.mass_density_comoving = rho_star;

    const double velocity_projection = dot3(prim.vel_x_peculiar, prim.vel_y_peculiar, prim.vel_z_peculiar,
        face.normal_x, face.normal_y, face.normal_z);
    const double delta_normal = s_star - velocity_projection;
    star.momentum_density_x_comoving = rho_star * (prim.vel_x_peculiar + delta_normal * face.normal_x);
    star.momentum_density_y_comoving = rho_star * (prim.vel_y_peculiar + delta_normal * face.normal_y);
    star.momentum_density_z_comoving = rho_star * (prim.vel_z_peculiar + delta_normal * face.normal_z);

    const double denom_energy = prim.rho_comoving * (s_k - u_k_n);
    if (std::abs(denom_energy) < k_small) {
      return HydroConservedState{};
    }
    const double energy_num =
        (s_k - u_k_n) * cons.total_energy_density_comoving - prim.pressure_comoving * u_k_n + p_star * s_star;
    star.total_energy_density_comoving = rho_star * energy_num / denom_energy;
    if (!std::isfinite(star.total_energy_density_comoving) || star.total_energy_density_comoving <= 0.0) {
      return HydroConservedState{};
    }
    star.metal_mass_density_comoving =
        rho_star * std::clamp(prim.metallicity_mass_fraction, 0.0, 1.0);
    return star;
  };

  const HydroConservedState u_star_left = starState(left_state, u_left, s_left, un_left);
  const HydroConservedState u_star_right = starState(right_state, u_right, s_right, un_right);

  if (u_star_left.mass_density_comoving <= 0.0 || u_star_right.mass_density_comoving <= 0.0 ||
      !std::isfinite(u_star_left.total_energy_density_comoving) || !std::isfinite(u_star_right.total_energy_density_comoving)) {
    ++m_fallback_count;
    return hllFlux(left_state, right_state, face, adiabatic_index);
  }

  if (s_star >= 0.0) {
    return f_left + s_left * (u_star_left - u_left);
  }
  return f_right + s_right * (u_star_right - u_right);
}

std::string_view hydroRiemannSolverKindToString(HydroRiemannSolverKind kind) {
  switch (kind) {
    case HydroRiemannSolverKind::kHllc:
      return "hllc";
    case HydroRiemannSolverKind::kHlle:
      return "hlle";
  }
  return "unknown";
}

HydroRiemannSolverKind hydroRiemannSolverKindFromString(std::string_view name) {
  if (name == "hllc") {
    return HydroRiemannSolverKind::kHllc;
  }
  if (name == "hlle") {
    return HydroRiemannSolverKind::kHlle;
  }
  throw std::invalid_argument("Unsupported hydro Riemann solver: " + std::string(name));
}

std::unique_ptr<HydroRiemannSolver> makeHydroRiemannSolver(HydroRiemannSolverKind kind) {
  switch (kind) {
    case HydroRiemannSolverKind::kHllc:
      return std::make_unique<HllcRiemannSolver>();
    case HydroRiemannSolverKind::kHlle:
      return std::make_unique<HlleRiemannSolver>();
  }
  throw std::invalid_argument("unhandled HydroRiemannSolverKind");
}

std::uint64_t HllcRiemannSolver::fallbackCount() const {
  return m_fallback_count;
}

}  // namespace cosmosim::hydro
