#include "cosmosim/mhd/mhd_riemann.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace cosmosim::mhd {
namespace {

struct LocalPrimitive {
  double rho = 0.0;
  double vn = 0.0;
  double vt1 = 0.0;
  double vt2 = 0.0;
  double p = 0.0;
  double bn = 0.0;
  double bt1 = 0.0;
  double bt2 = 0.0;
};

// Eight-component local conserved vector in the order
// (rho, m_n, m_t1, m_t2, E, B_n, B_t1, B_t2).
using LocalVector = std::array<double, 8>;

enum LocalIndex : std::size_t {
  kRho = 0,
  kMn = 1,
  kMt1 = 2,
  kMt2 = 3,
  kEnergy = 4,
  kBn = 5,
  kBt1 = 6,
  kBt2 = 7,
};

LocalPrimitive toLocal(const MhdPrimitiveState& primitive, MhdAxis axis) noexcept {
  const std::array<double, 3> velocity{
      primitive.vel_x_peculiar, primitive.vel_y_peculiar, primitive.vel_z_peculiar};
  const std::array<double, 3> field{
      primitive.bx_comoving, primitive.by_comoving, primitive.bz_comoving};
  const std::array<double, 3> local_velocity = permuteToLocal(axis, velocity);
  const std::array<double, 3> local_field = permuteToLocal(axis, field);
  LocalPrimitive local{};
  local.rho = primitive.rho_comoving;
  local.vn = local_velocity[0];
  local.vt1 = local_velocity[1];
  local.vt2 = local_velocity[2];
  local.p = primitive.pressure_comoving;
  local.bn = local_field[0];
  local.bt1 = local_field[1];
  local.bt2 = local_field[2];
  return local;
}

double magneticPressure(const LocalPrimitive& state) noexcept {
  return 0.5 * (state.bn * state.bn + state.bt1 * state.bt1 + state.bt2 * state.bt2);
}

double totalPressure(const LocalPrimitive& state) noexcept {
  return state.p + magneticPressure(state);
}

LocalVector conservedFromPrimitive(const LocalPrimitive& state, double gamma) noexcept {
  LocalVector u{};
  u[kRho] = state.rho;
  u[kMn] = state.rho * state.vn;
  u[kMt1] = state.rho * state.vt1;
  u[kMt2] = state.rho * state.vt2;
  const double kinetic = 0.5 * state.rho *
      (state.vn * state.vn + state.vt1 * state.vt1 + state.vt2 * state.vt2);
  const double magnetic = magneticPressure(state);
  u[kEnergy] = state.p / (gamma - 1.0) + kinetic + magnetic;
  u[kBn] = state.bn;
  u[kBt1] = state.bt1;
  u[kBt2] = state.bt2;
  return u;
}

// Ideal-MHD flux along the local normal.
LocalVector fluxFromPrimitive(const LocalPrimitive& state, double gamma) noexcept {
  const LocalVector u = conservedFromPrimitive(state, gamma);
  const double pt = totalPressure(state);
  const double v_dot_b = state.vn * state.bn + state.vt1 * state.bt1 + state.vt2 * state.bt2;
  LocalVector flux{};
  flux[kRho] = state.rho * state.vn;
  flux[kMn] = state.rho * state.vn * state.vn + pt - state.bn * state.bn;
  flux[kMt1] = state.rho * state.vn * state.vt1 - state.bn * state.bt1;
  flux[kMt2] = state.rho * state.vn * state.vt2 - state.bn * state.bt2;
  flux[kEnergy] = (u[kEnergy] + pt) * state.vn - state.bn * v_dot_b;
  flux[kBn] = 0.0;
  flux[kBt1] = state.vn * state.bt1 - state.vt1 * state.bn;
  flux[kBt2] = state.vn * state.bt2 - state.vt2 * state.bn;
  return flux;
}

double fastSpeed(const LocalPrimitive& state, double gamma) noexcept {
  if (!(state.rho > 0.0) || !(gamma > 1.0)) {
    return 0.0;
  }
  const double sound_squared = gamma * state.p / state.rho;
  const double alfven_squared =
      (state.bn * state.bn + state.bt1 * state.bt1 + state.bt2 * state.bt2) / state.rho;
  const double alfven_normal_squared = state.bn * state.bn / state.rho;
  const double sum = sound_squared + alfven_squared;
  const double discriminant = std::max(
      0.0, sum * sum - 4.0 * sound_squared * alfven_normal_squared);
  return std::sqrt(std::max(0.0, 0.5 * (sum + std::sqrt(discriminant))));
}

LocalVector globalFluxFromLocal(const LocalVector& local_flux, MhdAxis axis) noexcept {
  const std::array<double, 3> local_momentum{local_flux[kMn], local_flux[kMt1], local_flux[kMt2]};
  const std::array<double, 3> local_field{local_flux[kBn], local_flux[kBt1], local_flux[kBt2]};
  const std::array<double, 3> global_momentum = permuteToGlobal(axis, local_momentum);
  const std::array<double, 3> global_field = permuteToGlobal(axis, local_field);
  LocalVector global{};
  global[kRho] = local_flux[kRho];
  global[kMn] = global_momentum[0];
  global[kMt1] = global_momentum[1];
  global[kMt2] = global_momentum[2];
  global[kEnergy] = local_flux[kEnergy];
  global[kBn] = global_field[0];
  global[kBt1] = global_field[1];
  global[kBt2] = global_field[2];
  return global;
}

MhdRiemannFlux buildFluxResult(const LocalVector& local_flux,
                               MhdAxis axis,
                               double max_fast,
                               bool used_hlle) noexcept {
  const LocalVector global = globalFluxFromLocal(local_flux, axis);
  MhdRiemannFlux result{};
  result.mass_flux = global[kRho];
  result.momentum_flux = {global[kMn], global[kMt1], global[kMt2]};
  result.energy_flux = global[kEnergy];
  result.magnetic_flux = {global[kBn], global[kBt1], global[kBt2]};
  // The electric-field components are recovered from the magnetic flux lanes:
  //   E_t1 = F_B,t2,  E_t2 = -F_B,t1,  E_n = 0.
  // This relation holds identically for both HLLD and HLLE because the
  // magnetic flux is F_B = v_n B - B_n v and E = -v x B.
  result.electric_field = {0.0, local_flux[kBt2], -local_flux[kBt1]};
  result.max_fast_speed = max_fast;
  result.slowest_signal_speed = 0.0;
  result.fastest_signal_speed = max_fast;
  result.used_hlle_fallback = used_hlle;
  return result;
}

bool finiteState(const LocalPrimitive& state) noexcept {
  return std::isfinite(state.rho) && std::isfinite(state.p) &&
      std::isfinite(state.vn) && std::isfinite(state.vt1) && std::isfinite(state.vt2) &&
      std::isfinite(state.bn) && std::isfinite(state.bt1) && std::isfinite(state.bt2) &&
      state.rho > 0.0 && state.p > 0.0;
}

}  // namespace

std::array<double, 3> permuteToLocal(MhdAxis axis,
                                     const std::array<double, 3>& global) noexcept {
  switch (axis) {
    case MhdAxis::kX:
      return {global[0], global[1], global[2]};
    case MhdAxis::kY:
      return {global[1], global[2], global[0]};
    case MhdAxis::kZ:
      return {global[2], global[0], global[1]};
  }
  return {global[0], global[1], global[2]};
}

std::array<double, 3> permuteToGlobal(MhdAxis axis,
                                      const std::array<double, 3>& local) noexcept {
  switch (axis) {
    case MhdAxis::kX:
      return {local[0], local[1], local[2]};
    case MhdAxis::kY:
      return {local[2], local[0], local[1]};
    case MhdAxis::kZ:
      return {local[1], local[2], local[0]};
  }
  return {local[0], local[1], local[2]};
}

MhdConservedState& MhdConservedState::operator+=(const MhdConservedState& rhs) noexcept {
  rho_comoving += rhs.rho_comoving;
  momentum_density_x_comoving += rhs.momentum_density_x_comoving;
  momentum_density_y_comoving += rhs.momentum_density_y_comoving;
  momentum_density_z_comoving += rhs.momentum_density_z_comoving;
  total_energy_density_comoving += rhs.total_energy_density_comoving;
  bx_comoving += rhs.bx_comoving;
  by_comoving += rhs.by_comoving;
  bz_comoving += rhs.bz_comoving;
  metal_mass_density_comoving += rhs.metal_mass_density_comoving;
  return *this;
}

MhdConservedState& MhdConservedState::operator-=(const MhdConservedState& rhs) noexcept {
  rho_comoving -= rhs.rho_comoving;
  momentum_density_x_comoving -= rhs.momentum_density_x_comoving;
  momentum_density_y_comoving -= rhs.momentum_density_y_comoving;
  momentum_density_z_comoving -= rhs.momentum_density_z_comoving;
  total_energy_density_comoving -= rhs.total_energy_density_comoving;
  bx_comoving -= rhs.bx_comoving;
  by_comoving -= rhs.by_comoving;
  bz_comoving -= rhs.bz_comoving;
  metal_mass_density_comoving -= rhs.metal_mass_density_comoving;
  return *this;
}

MhdConservedState operator+(MhdConservedState lhs, const MhdConservedState& rhs) noexcept {
  lhs += rhs;
  return lhs;
}

MhdConservedState operator-(MhdConservedState lhs, const MhdConservedState& rhs) noexcept {
  lhs -= rhs;
  return lhs;
}

MhdConservedState operator*(double scalar, MhdConservedState state) noexcept {
  state.rho_comoving *= scalar;
  state.momentum_density_x_comoving *= scalar;
  state.momentum_density_y_comoving *= scalar;
  state.momentum_density_z_comoving *= scalar;
  state.total_energy_density_comoving *= scalar;
  state.bx_comoving *= scalar;
  state.by_comoving *= scalar;
  state.bz_comoving *= scalar;
  state.metal_mass_density_comoving *= scalar;
  return state;
}

bool recoverPrimitive(const MhdConservedState& conserved,
                      double gamma_adiabatic,
                      MhdPrimitiveState& primitive) noexcept {
  if (!(gamma_adiabatic > 1.0) || !(conserved.rho_comoving > 0.0) ||
      !std::isfinite(conserved.rho_comoving)) {
    return false;
  }
  const double rho = conserved.rho_comoving;
  const double vx = conserved.momentum_density_x_comoving / rho;
  const double vy = conserved.momentum_density_y_comoving / rho;
  const double vz = conserved.momentum_density_z_comoving / rho;
  const double kinetic =
      0.5 * rho * (vx * vx + vy * vy + vz * vz);
  const double magnetic = 0.5 * (conserved.bx_comoving * conserved.bx_comoving +
                                 conserved.by_comoving * conserved.by_comoving +
                                 conserved.bz_comoving * conserved.bz_comoving);
  const double pressure = (conserved.total_energy_density_comoving - kinetic - magnetic) *
      (gamma_adiabatic - 1.0);
  if (!std::isfinite(pressure) || !(pressure > 0.0) ||
      !std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz)) {
    return false;
  }
  primitive.rho_comoving = rho;
  primitive.vel_x_peculiar = vx;
  primitive.vel_y_peculiar = vy;
  primitive.vel_z_peculiar = vz;
  primitive.pressure_comoving = pressure;
  primitive.bx_comoving = conserved.bx_comoving;
  primitive.by_comoving = conserved.by_comoving;
  primitive.bz_comoving = conserved.bz_comoving;
  primitive.metallicity_mass_fraction = conserved.metal_mass_density_comoving / rho;
  return true;
}

MhdConservedState primitiveToConserved(const MhdPrimitiveState& primitive,
                                       double gamma_adiabatic) noexcept {
  MhdConservedState conserved{};
  conserved.rho_comoving = primitive.rho_comoving;
  conserved.momentum_density_x_comoving = primitive.rho_comoving * primitive.vel_x_peculiar;
  conserved.momentum_density_y_comoving = primitive.rho_comoving * primitive.vel_y_peculiar;
  conserved.momentum_density_z_comoving = primitive.rho_comoving * primitive.vel_z_peculiar;
  const double kinetic = 0.5 * primitive.rho_comoving *
      (primitive.vel_x_peculiar * primitive.vel_x_peculiar +
       primitive.vel_y_peculiar * primitive.vel_y_peculiar +
       primitive.vel_z_peculiar * primitive.vel_z_peculiar);
  const double magnetic = primitive.magneticEnergyDensity();
  conserved.total_energy_density_comoving =
      primitive.pressure_comoving / (gamma_adiabatic - 1.0) + kinetic + magnetic;
  conserved.bx_comoving = primitive.bx_comoving;
  conserved.by_comoving = primitive.by_comoving;
  conserved.bz_comoving = primitive.bz_comoving;
  conserved.metal_mass_density_comoving =
      primitive.rho_comoving * primitive.metallicity_mass_fraction;
  return conserved;
}

double fastMagnetosonicSpeed(const MhdPrimitiveState& primitive,
                             double gamma_adiabatic,
                             MhdAxis axis) noexcept {
  return fastSpeed(toLocal(primitive, axis), gamma_adiabatic);
}

double alfvenSpeedAlong(const MhdPrimitiveState& primitive, MhdAxis axis) noexcept {
  if (!(primitive.rho_comoving > 0.0)) {
    return 0.0;
  }
  const LocalPrimitive local = toLocal(primitive, axis);
  return std::abs(local.bn) / std::sqrt(local.rho);
}

MhdRiemann::Result MhdRiemann::solveHlle(MhdAxis axis,
                                         const MhdPrimitiveState& left,
                                         const MhdPrimitiveState& right,
                                         double shared_normal_b_comoving,
                                         double gamma_adiabatic) noexcept {
  Result result{};
  LocalPrimitive local_left = toLocal(left, axis);
  LocalPrimitive local_right = toLocal(right, axis);
  local_left.bn = shared_normal_b_comoving;
  local_right.bn = shared_normal_b_comoving;
  if (!finiteState(local_left) || !finiteState(local_right)) {
    return result;
  }
  const double fast_left = fastSpeed(local_left, gamma_adiabatic);
  const double fast_right = fastSpeed(local_right, gamma_adiabatic);
  const double speed_left = std::min(local_left.vn - fast_left, local_right.vn - fast_right);
  const double speed_right = std::max(local_left.vn + fast_left, local_right.vn + fast_right);
  const LocalVector u_left = conservedFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector u_right = conservedFromPrimitive(local_right, gamma_adiabatic);
  const LocalVector f_left = fluxFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector f_right = fluxFromPrimitive(local_right, gamma_adiabatic);
  LocalVector flux{};
  if (speed_left >= 0.0) {
    flux = f_left;
  } else if (speed_right <= 0.0) {
    flux = f_right;
  } else {
    const double denominator = speed_right - speed_left;
    if (!(denominator > 0.0) || !std::isfinite(denominator)) {
      return result;
    }
    for (std::size_t index = 0; index < flux.size(); ++index) {
      flux[index] = (speed_right * f_left[index] - speed_left * f_right[index] +
                     speed_left * speed_right * (u_right[index] - u_left[index])) /
          denominator;
    }
  }
  result.flux = buildFluxResult(flux, axis, std::max(fast_left, fast_right), true);
  result.valid = true;
  result.used_hlle = true;
  return result;
}

MhdRiemann::Result MhdRiemann::solve(MhdAxis axis,
                                     const MhdPrimitiveState& left,
                                     const MhdPrimitiveState& right,
                                     double shared_normal_b_comoving,
                                     double gamma_adiabatic,
                                     Solver preferred) noexcept {
  Result result{};
  LocalPrimitive local_left = toLocal(left, axis);
  LocalPrimitive local_right = toLocal(right, axis);
  // The single face-shared normal component overrides both sides.
  local_left.bn = shared_normal_b_comoving;
  local_right.bn = shared_normal_b_comoving;
  if (!finiteState(local_left) || !finiteState(local_right)) {
    return result;
  }
  if (preferred == Solver::kHlle) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }

  const double fast_left = fastSpeed(local_left, gamma_adiabatic);
  const double fast_right = fastSpeed(local_right, gamma_adiabatic);
  const double speed_left = std::min(local_left.vn - fast_left, local_right.vn - fast_right);
  const double speed_right = std::max(local_left.vn + fast_left, local_right.vn + fast_right);
  const LocalVector u_left = conservedFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector u_right = conservedFromPrimitive(local_right, gamma_adiabatic);
  const LocalVector f_left = fluxFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector f_right = fluxFromPrimitive(local_right, gamma_adiabatic);
  const double pt_left = totalPressure(local_left);
  const double pt_right = totalPressure(local_right);

  if (!(speed_right > speed_left)) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }
  if (speed_left >= 0.0) {
    result.flux = buildFluxResult(f_left, axis, fast_left, false);
    result.valid = true;
    return result;
  }
  if (speed_right <= 0.0) {
    result.flux = buildFluxResult(f_right, axis, fast_right, false);
    result.valid = true;
    return result;
  }

  const double denom_contact =
      (speed_right - local_right.vn) * local_right.rho -
      (speed_left - local_left.vn) * local_left.rho;
  if (!(std::abs(denom_contact) > std::numeric_limits<double>::min())) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }
  const double speed_contact =
      ((speed_right - local_right.vn) * local_right.rho * local_right.vn -
       (speed_left - local_left.vn) * local_left.rho * local_left.vn -
       pt_right + pt_left) /
      denom_contact;
  if (!std::isfinite(speed_contact) || speed_contact <= speed_left ||
      speed_contact >= speed_right) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }

  const double pt_star = pt_left +
      local_left.rho * (speed_left - local_left.vn) * (speed_contact - local_left.vn);
  if (!std::isfinite(pt_star) || pt_star <= 0.0) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }

  struct StarSide {
    double rho = 0.0;
    double vt1 = 0.0;
    double vt2 = 0.0;
    double bt1 = 0.0;
    double bt2 = 0.0;
    double energy = 0.0;
  };
  const auto computeStar = [&](const LocalPrimitive& state, const LocalVector& u,
                               double speed) -> StarSide {
    StarSide star{};
    const double denom =
        state.rho * (speed - state.vn) * (speed - speed_contact) -
        state.bn * state.bn;
    if (!(std::abs(denom) > std::numeric_limits<double>::min()) || !std::isfinite(denom)) {
      return star;
    }
    star.rho = state.rho * (speed - state.vn) / (speed - speed_contact);
    star.vt1 = state.vt1 -
        state.bn * state.bt1 * (speed_contact - state.vn) / denom;
    star.vt2 = state.vt2 -
        state.bn * state.bt2 * (speed_contact - state.vn) / denom;
    star.bt1 = state.bt1 *
        (state.rho * (speed - state.vn) * (speed - state.vn) - state.bn * state.bn) / denom;
    star.bt2 = state.bt2 *
        (state.rho * (speed - state.vn) * (speed - state.vn) - state.bn * state.bn) / denom;
    const double v_dot_b =
        state.vn * state.bn + state.vt1 * state.bt1 + state.vt2 * state.bt2;
    const double v_star_dot_b_star =
        speed_contact * state.bn + star.vt1 * star.bt1 + star.vt2 * star.bt2;
    const double numerator =
        (speed - state.vn) * u[kEnergy] - totalPressure(state) * state.vn +
        pt_star * speed_contact + state.bn * (v_dot_b - v_star_dot_b_star);
    star.energy = numerator / (speed - speed_contact);
    return star;
  };

  const StarSide star_left = computeStar(local_left, u_left, speed_left);
  const StarSide star_right = computeStar(local_right, u_right, speed_right);
  const auto starValid = [](const StarSide& star) {
    return std::isfinite(star.rho) && star.rho > 0.0 && std::isfinite(star.energy) &&
        std::isfinite(star.vt1) && std::isfinite(star.vt2) && std::isfinite(star.bt1) &&
        std::isfinite(star.bt2);
  };
  if (!starValid(star_left) || !starValid(star_right)) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }

  const auto starVector = [&](const LocalPrimitive& state, const StarSide& star) {
    LocalVector u{};
    u[kRho] = star.rho;
    u[kMn] = star.rho * speed_contact;
    u[kMt1] = star.rho * star.vt1;
    u[kMt2] = star.rho * star.vt2;
    u[kEnergy] = star.energy;
    u[kBn] = state.bn;
    u[kBt1] = star.bt1;
    u[kBt2] = star.bt2;
    return u;
  };

  const LocalVector u_star_left = starVector(local_left, star_left);
  const LocalVector u_star_right = starVector(local_right, star_right);
  LocalVector f_star_left{};
  LocalVector f_star_right{};
  for (std::size_t index = 0; index < f_star_left.size(); ++index) {
    f_star_left[index] = f_left[index] + speed_left * (u_star_left[index] - u_left[index]);
    f_star_right[index] = f_right[index] + speed_right * (u_star_right[index] - u_right[index]);
  }

  const double alfven_left = std::abs(shared_normal_b_comoving) / std::sqrt(star_left.rho);
  const double alfven_right = std::abs(shared_normal_b_comoving) / std::sqrt(star_right.rho);
  const double speed_star_left = speed_contact - alfven_left;
  const double speed_star_right = speed_contact + alfven_right;

  // Double-star (tangential discontinuity) states.
  const double sqrt_rho_left = std::sqrt(star_left.rho);
  const double sqrt_rho_right = std::sqrt(star_right.rho);
  const double sqrt_sum = sqrt_rho_left + sqrt_rho_right;
  const double sign_bn = shared_normal_b_comoving > 0.0 ? 1.0
      : (shared_normal_b_comoving < 0.0 ? -1.0 : 0.0);
  StarSide double_left = star_left;
  StarSide double_right = star_right;

  const auto buildDouble = [&](const StarSide& star, double side_sign, StarSide& out) {
    out.rho = star.rho;
    out.vt1 = (sqrt_rho_left * star_left.vt1 + sqrt_rho_right * star_right.vt1 +
               sign_bn * (star_right.bt1 - star_left.bt1)) /
        sqrt_sum;
    out.vt2 = (sqrt_rho_left * star_left.vt2 + sqrt_rho_right * star_right.vt2 +
               sign_bn * (star_right.bt2 - star_left.bt2)) /
        sqrt_sum;
    out.bt1 = (sqrt_rho_left * star_right.bt1 + sqrt_rho_right * star_left.bt1 +
               sign_bn * sqrt_rho_left * sqrt_rho_right *
                   (star_right.vt1 - star_left.vt1)) /
        sqrt_sum;
    out.bt2 = (sqrt_rho_left * star_right.bt2 + sqrt_rho_right * star_left.bt2 +
               sign_bn * sqrt_rho_left * sqrt_rho_right *
                   (star_right.vt2 - star_left.vt2)) /
        sqrt_sum;
    const double v_dot_b_star =
        speed_contact * shared_normal_b_comoving + star.vt1 * star.bt1 +
        star.vt2 * star.bt2;
    const double v_dot_b_double =
        speed_contact * shared_normal_b_comoving + out.vt1 * out.bt1 +
        out.vt2 * out.bt2;
    out.energy = star.energy -
        side_sign * std::sqrt(star.rho) * (v_dot_b_star - v_dot_b_double) * sign_bn;
    return std::isfinite(out.energy);
  };

  const bool double_ok_left = buildDouble(star_left, 1.0, double_left);
  const bool double_ok_right = buildDouble(star_right, -1.0, double_right);

  const auto doubleVector = [&](const LocalPrimitive& state, const StarSide& star) {
    LocalVector u{};
    u[kRho] = star.rho;
    u[kMn] = star.rho * speed_contact;
    u[kMt1] = star.rho * star.vt1;
    u[kMt2] = star.rho * star.vt2;
    u[kEnergy] = star.energy;
    u[kBn] = state.bn;
    u[kBt1] = star.bt1;
    u[kBt2] = star.bt2;
    return u;
  };

  LocalVector flux = f_left;
  if (speed_star_left >= 0.0) {
    flux = f_star_left;
  } else if (speed_contact >= 0.0) {
    if (!double_ok_left) {
      return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
    }
    const LocalVector u_double_left = doubleVector(local_left, double_left);
    for (std::size_t index = 0; index < flux.size(); ++index) {
      flux[index] = f_star_left[index] +
          speed_star_left * (u_double_left[index] - u_star_left[index]);
    }
  } else if (speed_star_right <= 0.0) {
    if (!double_ok_right) {
      return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
    }
    const LocalVector u_double_right = doubleVector(local_right, double_right);
    for (std::size_t index = 0; index < flux.size(); ++index) {
      flux[index] = f_star_right[index] +
          speed_star_right * (u_double_right[index] - u_star_right[index]);
    }
  } else {
    flux = f_star_right;
  }

  if (!std::isfinite(flux[kRho]) || !std::isfinite(flux[kEnergy])) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }
  result.flux = buildFluxResult(flux, axis, std::max(fast_left, fast_right), false);
  result.valid = true;
  return result;
}

MhdRiemannFlux firstOrderPositivityFlux(MhdAxis axis,
                                        const MhdPrimitiveState& left,
                                        const MhdPrimitiveState& right,
                                        double shared_normal_b_comoving,
                                        double gamma_adiabatic) noexcept {
  MhdRiemannFlux result{};
  LocalPrimitive local_left = toLocal(left, axis);
  LocalPrimitive local_right = toLocal(right, axis);
  local_left.bn = shared_normal_b_comoving;
  local_right.bn = shared_normal_b_comoving;
  if (!finiteState(local_left) || !finiteState(local_right)) {
    return result;
  }
  const double fast_left = fastSpeed(local_left, gamma_adiabatic);
  const double fast_right = fastSpeed(local_right, gamma_adiabatic);
  const double speed = std::max(std::abs(local_left.vn) + fast_left,
                                std::abs(local_right.vn) + fast_right);
  const LocalVector u_left = conservedFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector u_right = conservedFromPrimitive(local_right, gamma_adiabatic);
  const LocalVector f_left = fluxFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector f_right = fluxFromPrimitive(local_right, gamma_adiabatic);
  LocalVector flux{};
  for (std::size_t index = 0; index < flux.size(); ++index) {
    flux[index] = 0.5 * (f_left[index] + f_right[index]) -
        0.5 * speed * (u_right[index] - u_left[index]);
  }
  return buildFluxResult(flux, axis, std::max(fast_left, fast_right), true);
}

}  // namespace cosmosim::mhd
