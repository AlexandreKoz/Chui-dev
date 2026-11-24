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
                               double speed_left,
                               double speed_right,
                               bool used_hlle) noexcept {
  const LocalVector global = globalFluxFromLocal(local_flux, axis);
  MhdRiemannFlux result{};
  result.mass_flux = global[kRho];
  result.momentum_flux = {global[kMn], global[kMt1], global[kMt2]};
  result.energy_flux = global[kEnergy];
  result.magnetic_flux = {global[kBn], global[kBt1], global[kBt2]};
  // Ideal electric field in the local right-handed basis (n, t1, t2):
  //   F_B,t1 = v_n B_t1 - v_t1 B_n = -E_t2,  F_B,t2 = v_n B_t2 - v_t2 B_n = E_t1,
  // so E_t1 = F_B,t2 and E_t2 = -F_B,t1 (E = -v x B). Stored in LOCAL order
  // {E_n=0, E_t1, E_t2}; callers map to global components with permuteToGlobal.
  result.electric_field = {0.0, local_flux[kBt2], -local_flux[kBt1]};
  result.max_fast_speed = max_fast;
  // Signed outermost fan bounds (S_L <= S_R) in the face-normal direction.
  result.slowest_signal_speed = speed_left;
  result.fastest_signal_speed = speed_right;
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
  if (!finiteState(local_left) || !finiteState(local_right) || !(gamma_adiabatic > 1.0)) {
    return result;
  }
  const double fast_left = fastSpeed(local_left, gamma_adiabatic);
  const double fast_right = fastSpeed(local_right, gamma_adiabatic);
  // Davis/Einfeldt outer bounds.
  const double speed_left = std::min(local_left.vn - fast_left, local_right.vn - fast_right);
  const double speed_right = std::max(local_left.vn + fast_left, local_right.vn + fast_right);
  if (!std::isfinite(speed_left) || !std::isfinite(speed_right) ||
      !(speed_right > speed_left)) {
    return result;
  }
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
    for (std::size_t index = 0; index < flux.size(); ++index) {
      flux[index] = (speed_right * f_left[index] - speed_left * f_right[index] +
                     speed_left * speed_right * (u_right[index] - u_left[index])) /
          denominator;
    }
  }
  flux[kBn] = 0.0;
  for (const double value : flux) {
    if (!std::isfinite(value)) {
      return result;
    }
  }
  result.flux = buildFluxResult(flux, axis, std::max(fast_left, fast_right), speed_left,
                                speed_right, true);
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
  if (preferred == Solver::kHlle) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }
  LocalPrimitive local_left = toLocal(left, axis);
  LocalPrimitive local_right = toLocal(right, axis);
  // The single face-shared normal component overrides both sides.
  local_left.bn = shared_normal_b_comoving;
  local_right.bn = shared_normal_b_comoving;
  if (!finiteState(local_left) || !finiteState(local_right) || !(gamma_adiabatic > 1.0)) {
    return result;
  }
  const double bn = shared_normal_b_comoving;
  const double fast_left = fastSpeed(local_left, gamma_adiabatic);
  const double fast_right = fastSpeed(local_right, gamma_adiabatic);
  const double speed_left = std::min(local_left.vn - fast_left, local_right.vn - fast_right);
  const double speed_right = std::max(local_left.vn + fast_left, local_right.vn + fast_right);
  const double max_fast = std::max(fast_left, fast_right);
  if (!std::isfinite(speed_left) || !std::isfinite(speed_right) ||
      !(speed_right > speed_left)) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }
  const LocalVector u_left = conservedFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector u_right = conservedFromPrimitive(local_right, gamma_adiabatic);
  const LocalVector f_left = fluxFromPrimitive(local_left, gamma_adiabatic);
  const LocalVector f_right = fluxFromPrimitive(local_right, gamma_adiabatic);
  if (speed_left >= 0.0) {
    result.flux = buildFluxResult(f_left, axis, max_fast, speed_left, speed_right, false);
    result.valid = true;
    return result;
  }
  if (speed_right <= 0.0) {
    result.flux = buildFluxResult(f_right, axis, max_fast, speed_left, speed_right, false);
    result.valid = true;
    return result;
  }
  const double pt_left = totalPressure(local_left);
  const double pt_right = totalPressure(local_right);

  // Miyoshi & Kusano (2005) eq. (38): contact (entropy) speed.
  const double sdl = speed_left - local_left.vn;    // < 0
  const double sdr = speed_right - local_right.vn;  // > 0
  const double denom_contact = sdr * local_right.rho - sdl * local_left.rho;
  if (!(denom_contact > 0.0) || !std::isfinite(denom_contact)) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }
  const double speed_contact =
      (sdr * local_right.rho * local_right.vn - sdl * local_left.rho * local_left.vn -
       pt_right + pt_left) /
      denom_contact;
  if (!std::isfinite(speed_contact) || !(speed_contact > speed_left) ||
      !(speed_contact < speed_right)) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }
  // Eq. (41), symmetric form of the star total pressure.
  const double pt_star =
      (sdr * local_right.rho * pt_left - sdl * local_left.rho * pt_right +
       local_left.rho * local_right.rho * sdr * sdl * (local_right.vn - local_left.vn)) /
      denom_contact;
  if (!std::isfinite(pt_star) || !(pt_star > 0.0)) {
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
  // Relative threshold for the degenerate (S_L* -> S_L, B_t -> 0 and
  // B_n^2 >= gamma p) configuration of eqs. (44)-(47): the tangential
  // components are then continuous across the fast wave (MK05 section 4).
  constexpr double k_degenerate_relative = 1.0e-12;
  const auto computeStar = [&](const LocalPrimitive& state, const LocalVector& u,
                               double speed) -> StarSide {
    StarSide star{};
    const double sd = speed - state.vn;
    const double sm_minus_u = speed_contact - state.vn;
    star.rho = state.rho * sd / (speed - speed_contact);
    const double denom = state.rho * sd * (speed - speed_contact) - bn * bn;
    if (std::abs(denom) <= k_degenerate_relative * pt_star) {
      star.vt1 = state.vt1;
      star.vt2 = state.vt2;
      star.bt1 = state.bt1;
      star.bt2 = state.bt2;
    } else {
      const double velocity_factor = bn * sm_minus_u / denom;
      const double field_factor = (state.rho * sd * sd - bn * bn) / denom;
      star.vt1 = state.vt1 - state.bt1 * velocity_factor;
      star.vt2 = state.vt2 - state.bt2 * velocity_factor;
      star.bt1 = state.bt1 * field_factor;
      star.bt2 = state.bt2 * field_factor;
    }
    const double v_dot_b =
        state.vn * bn + state.vt1 * state.bt1 + state.vt2 * state.bt2;
    const double v_star_dot_b_star =
        speed_contact * bn + star.vt1 * star.bt1 + star.vt2 * star.bt2;
    // Eq. (48).
    star.energy = (sd * u[kEnergy] - totalPressure(state) * state.vn +
                   pt_star * speed_contact + bn * (v_dot_b - v_star_dot_b_star)) /
        (speed - speed_contact);
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
  // Star-state thermal admissibility: E* - rho* |v*|^2/2 - |B*|^2/2 > 0.
  const auto starThermalPositive = [&](const StarSide& star) {
    const double kinetic = 0.5 * star.rho *
        (speed_contact * speed_contact + star.vt1 * star.vt1 + star.vt2 * star.vt2);
    const double magnetic = 0.5 * (bn * bn + star.bt1 * star.bt1 + star.bt2 * star.bt2);
    return star.energy - kinetic - magnetic > 0.0;
  };
  if (!starThermalPositive(star_left) || !starThermalPositive(star_right)) {
    return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
  }

  const auto stateVector = [&](const StarSide& star) {
    LocalVector u{};
    u[kRho] = star.rho;
    u[kMn] = star.rho * speed_contact;
    u[kMt1] = star.rho * star.vt1;
    u[kMt2] = star.rho * star.vt2;
    u[kEnergy] = star.energy;
    u[kBn] = bn;
    u[kBt1] = star.bt1;
    u[kBt2] = star.bt2;
    return u;
  };

  const LocalVector u_star_left = stateVector(star_left);
  const LocalVector u_star_right = stateVector(star_right);
  LocalVector f_star_left{};
  LocalVector f_star_right{};
  for (std::size_t index = 0; index < f_star_left.size(); ++index) {
    f_star_left[index] = f_left[index] + speed_left * (u_star_left[index] - u_left[index]);
    f_star_right[index] = f_right[index] + speed_right * (u_star_right[index] - u_right[index]);
  }

  LocalVector flux{};
  const double sqrt_rho_left = std::sqrt(star_left.rho);
  const double sqrt_rho_right = std::sqrt(star_right.rho);
  // Eq. (51): rotational (Alfven) discontinuity speeds.
  const double speed_star_left = speed_contact - std::abs(bn) / sqrt_rho_left;
  const double speed_star_right = speed_contact + std::abs(bn) / sqrt_rho_right;
  // B_n -> 0 collapses both rotational waves onto the contact; the double-star
  // region has zero width and the flux is F_L* or F_R* (MK05 section 4).
  const bool vanishing_normal_field = bn * bn <= k_degenerate_relative * pt_star;
  if (vanishing_normal_field) {
    flux = speed_contact >= 0.0 ? f_star_left : f_star_right;
  } else if (speed_star_left >= 0.0) {
    // S_L <= 0 <= S_L*.
    flux = f_star_left;
  } else if (speed_star_right <= 0.0) {
    // S_R* <= 0 <= S_R.
    flux = f_star_right;
  } else {
    // Double-star states, eqs. (59)-(63).
    const double sqrt_sum = sqrt_rho_left + sqrt_rho_right;
    const double sign_bn = bn > 0.0 ? 1.0 : -1.0;
    StarSide double_state{};
    double_state.vt1 = (sqrt_rho_left * star_left.vt1 + sqrt_rho_right * star_right.vt1 +
                        sign_bn * (star_right.bt1 - star_left.bt1)) /
        sqrt_sum;
    double_state.vt2 = (sqrt_rho_left * star_left.vt2 + sqrt_rho_right * star_right.vt2 +
                        sign_bn * (star_right.bt2 - star_left.bt2)) /
        sqrt_sum;
    double_state.bt1 = (sqrt_rho_left * star_right.bt1 + sqrt_rho_right * star_left.bt1 +
                        sign_bn * sqrt_rho_left * sqrt_rho_right *
                            (star_right.vt1 - star_left.vt1)) /
        sqrt_sum;
    double_state.bt2 = (sqrt_rho_left * star_right.bt2 + sqrt_rho_right * star_left.bt2 +
                        sign_bn * sqrt_rho_left * sqrt_rho_right *
                            (star_right.vt2 - star_left.vt2)) /
        sqrt_sum;
    const double v_dot_b_double = speed_contact * bn + double_state.vt1 * double_state.bt1 +
        double_state.vt2 * double_state.bt2;
    const auto vDotBStar = [&](const StarSide& star) {
      return speed_contact * bn + star.vt1 * star.bt1 + star.vt2 * star.bt2;
    };
    if (speed_contact >= 0.0) {
      // S_L* <= 0 <= S_M: F_L** = F_L* + S_L*(U_L** - U_L*).
      StarSide left_double = double_state;
      left_double.rho = star_left.rho;
      left_double.energy = star_left.energy -
          sqrt_rho_left * (vDotBStar(star_left) - v_dot_b_double) * sign_bn;
      const LocalVector u_double = stateVector(left_double);
      for (std::size_t index = 0; index < flux.size(); ++index) {
        flux[index] = f_star_left[index] +
            speed_star_left * (u_double[index] - u_star_left[index]);
      }
    } else {
      // S_M <= 0 <= S_R*: F_R** = F_R* + S_R*(U_R** - U_R*).
      StarSide right_double = double_state;
      right_double.rho = star_right.rho;
      right_double.energy = star_right.energy +
          sqrt_rho_right * (vDotBStar(star_right) - v_dot_b_double) * sign_bn;
      const LocalVector u_double = stateVector(right_double);
      for (std::size_t index = 0; index < flux.size(); ++index) {
        flux[index] = f_star_right[index] +
            speed_star_right * (u_double[index] - u_star_right[index]);
      }
    }
  }
  flux[kBn] = 0.0;
  for (const double value : flux) {
    if (!std::isfinite(value)) {
      return solveHlle(axis, left, right, shared_normal_b_comoving, gamma_adiabatic);
    }
  }
  result.flux = buildFluxResult(flux, axis, max_fast, speed_left, speed_right, false);
  result.valid = true;
  return result;
}

MhdSignalBounds MhdRiemann::signalBounds(MhdAxis axis,
                                         const MhdPrimitiveState& left,
                                         const MhdPrimitiveState& right,
                                         double shared_normal_b_comoving,
                                         double gamma_adiabatic) noexcept {
  MhdSignalBounds bounds{};
  LocalPrimitive local_left = toLocal(left, axis);
  LocalPrimitive local_right = toLocal(right, axis);
  local_left.bn = shared_normal_b_comoving;
  local_right.bn = shared_normal_b_comoving;
  if (!finiteState(local_left) || !finiteState(local_right) || !(gamma_adiabatic > 1.0)) {
    return bounds;
  }
  const double fast_left = fastSpeed(local_left, gamma_adiabatic);
  const double fast_right = fastSpeed(local_right, gamma_adiabatic);
  bounds.speed_left = std::min(local_left.vn - fast_left, local_right.vn - fast_right);
  bounds.speed_right = std::max(local_left.vn + fast_left, local_right.vn + fast_right);
  bounds.valid = std::isfinite(bounds.speed_left) && std::isfinite(bounds.speed_right) &&
      bounds.speed_right > bounds.speed_left;
  return bounds;
}

MhdRiemann::Result firstOrderPositivityFlux(MhdAxis axis,
                                            const MhdPrimitiveState& left_cell_average,
                                            const MhdPrimitiveState& right_cell_average,
                                            double shared_normal_b_comoving,
                                            double gamma_adiabatic) noexcept {
  // Local Lax-Friedrichs (Rusanov) flux between the two adjacent CELL-AVERAGE
  // states with the shared face normal field. With an admissible pair of
  // cell averages this flux is positivity preserving under the ideal CFL
  // bound. An inadmissible pair returns valid=false: there is no silent zero
  // flux; the caller must reject the stage and retry the interval.
  MhdRiemann::Result result{};
  LocalPrimitive local_left = toLocal(left_cell_average, axis);
  LocalPrimitive local_right = toLocal(right_cell_average, axis);
  local_left.bn = shared_normal_b_comoving;
  local_right.bn = shared_normal_b_comoving;
  if (!finiteState(local_left) || !finiteState(local_right) || !(gamma_adiabatic > 1.0)) {
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
  flux[kBn] = 0.0;
  for (const double value : flux) {
    if (!std::isfinite(value)) {
      return result;
    }
  }
  result.flux = buildFluxResult(flux, axis, std::max(fast_left, fast_right), -speed, speed,
                                true);
  result.valid = true;
  result.used_hlle = true;
  result.first_order_required = true;
  return result;
}

}  // namespace cosmosim::mhd
