#include "cosmosim/mhd/mhd_step_operator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace cosmosim::mhd {
namespace {

MhdFailureContext makeContext(std::string_view owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::string(owner);
  context.detail = std::move(detail);
  return context;
}

[[nodiscard]] double spacingAlong(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.dx_comoving
                    : (axis == 1U ? geometry.dy_comoving : geometry.dz_comoving);
}

[[nodiscard]] std::size_t faceIndex(const MhdPatchGeometry& geometry, std::size_t axis,
                                    const MhdIndex3& c) noexcept {
  const auto i = static_cast<std::size_t>(c[0]);
  const auto j = static_cast<std::size_t>(c[1]);
  const auto k = static_cast<std::size_t>(c[2]);
  return axis == 0U ? geometry.faceIndexX(i, j, k)
                    : (axis == 1U ? geometry.faceIndexY(i, j, k) : geometry.faceIndexZ(i, j, k));
}

[[nodiscard]] std::span<const double> faceLane(const MhdFaceFluxView& faces,
                                               std::size_t axis) noexcept {
  return axis == 0U ? faces.phi_x : (axis == 1U ? faces.phi_y : faces.phi_z);
}

[[nodiscard]] double faceArea(const MhdPatchGeometry& geometry, std::size_t axis) noexcept {
  return axis == 0U ? geometry.faceAreaX()
                    : (axis == 1U ? geometry.faceAreaY() : geometry.faceAreaZ());
}

// Cell-centred comoving field B_c: arithmetic mean of the two bordering face
// fields (CHUI-MHD-REF v1.0 section 5.1).
[[nodiscard]] std::array<double, 3> cellFieldComoving(const MhdPatchGeometry& geometry,
                                                      const MhdFaceFluxView& faces,
                                                      std::size_t i, std::size_t j,
                                                      std::size_t k) noexcept {
  return {0.5 * (faces.phi_x[geometry.faceIndexX(i, j, k)] +
                 faces.phi_x[geometry.faceIndexX(i + 1U, j, k)]) /
              geometry.faceAreaX(),
          0.5 * (faces.phi_y[geometry.faceIndexY(i, j, k)] +
                 faces.phi_y[geometry.faceIndexY(i, j + 1U, k)]) /
              geometry.faceAreaY(),
          0.5 * (faces.phi_z[geometry.faceIndexZ(i, j, k)] +
                 faces.phi_z[geometry.faceIndexZ(i, j, k + 1U)]) /
              geometry.faceAreaZ()};
}

// Recover a primitive state from comoving conserved gas lanes plus the
// face-derived B_c. Magnetic energy is B_c^2/(2a); the returned magnetic
// lanes are B_eff = B_c/sqrt(a). Returns false for any inadmissible state.
[[nodiscard]] bool recoverComoving(const MhdConservedState& gas,
                                   const std::array<double, 3>& b_comoving,
                                   double scale_factor,
                                   double gamma,
                                   MhdPrimitiveState& out) noexcept {
  const double rho = gas.rho_comoving;
  if (!(rho > 0.0) || !std::isfinite(rho) || !(scale_factor > 0.0)) {
    return false;
  }
  const double vx = gas.momentum_density_x_comoving / rho;
  const double vy = gas.momentum_density_y_comoving / rho;
  const double vz = gas.momentum_density_z_comoving / rho;
  const double kinetic = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
  const double b2 = b_comoving[0] * b_comoving[0] + b_comoving[1] * b_comoving[1] +
      b_comoving[2] * b_comoving[2];
  const double magnetic = 0.5 * b2 / scale_factor;
  const double thermal = gas.total_energy_density_comoving - kinetic - magnetic;
  const double pressure = (gamma - 1.0) * thermal;
  const double metal_fraction = gas.metal_mass_density_comoving / rho;
  if (!std::isfinite(pressure) || !(pressure > 0.0) || !std::isfinite(vx) ||
      !std::isfinite(vy) || !std::isfinite(vz) || !std::isfinite(b2) ||
      !std::isfinite(metal_fraction) || metal_fraction < 0.0 || metal_fraction > 1.0) {
    return false;
  }
  const double inverse_sqrt_a = 1.0 / std::sqrt(scale_factor);
  out.rho_comoving = rho;
  out.vel_x_peculiar = vx;
  out.vel_y_peculiar = vy;
  out.vel_z_peculiar = vz;
  out.pressure_comoving = pressure;
  out.bx_comoving = b_comoving[0] * inverse_sqrt_a;
  out.by_comoving = b_comoving[1] * inverse_sqrt_a;
  out.bz_comoving = b_comoving[2] * inverse_sqrt_a;
  out.metallicity_mass_fraction = metal_fraction;
  return true;
}

// ---------------------------------------------------------------------------
// MHD primitive eigensystem (local frame) for characteristic PPM.
//
// Primitive vector W = (rho, v_n, v_t1, v_t2, p, B_t1, B_t2); the normal field
// is a parameter. Right eigenvectors are derived from the quasi-linear
// primitive system W_t + A W_n = 0 (B_n constant) with the Roe/Ryu-Jones
// normalization (alpha_f, alpha_s, beta_t1, beta_t2) that stays regular at
// the triple-umbilic and B_t = 0 degeneracies. Column order:
//   0: v-c_f, 1: v-c_a, 2: v-c_s, 3: v (entropy), 4: v+c_s, 5: v+c_a, 6: v+c_f.
// Left eigenvectors are obtained by Gauss-Jordan inversion with partial
// pivoting; a singular basis is reported and the caller falls back to PLM.
// ---------------------------------------------------------------------------
using Vector7 = std::array<double, 7>;
using Matrix7 = std::array<Vector7, 7>;  // [row][column]

[[nodiscard]] bool buildRightEigenvectors(const Vector7& w, double bn, double gamma,
                                          Matrix7& right) noexcept {
  const double rho = w[0];
  const double pressure = w[4];
  if (!(rho > 0.0) || !(pressure > 0.0)) {
    return false;
  }
  const double sqrt_rho = std::sqrt(rho);
  const double a2 = gamma * pressure / rho;
  const double a = std::sqrt(a2);
  const double bt1 = w[5];
  const double bt2 = w[6];
  const double bx2 = bn * bn / rho;
  const double bt_mag2 = (bt1 * bt1 + bt2 * bt2) / rho;
  const double total2 = a2 + bx2 + bt_mag2;
  const double discriminant =
      std::sqrt(std::max(0.0, total2 * total2 - 4.0 * a2 * bx2));
  const double cf2 = 0.5 * (total2 + discriminant);
  const double cs2 = std::max(0.0, 0.5 * (total2 - discriminant));
  const double cf = std::sqrt(cf2);
  const double cs = std::sqrt(cs2);
  double alpha_f = 1.0;
  double alpha_s = 0.0;
  const double gap = cf2 - cs2;
  if (gap > 1.0e-12 * total2) {
    alpha_f = std::sqrt(std::clamp((a2 - cs2) / gap, 0.0, 1.0));
    alpha_s = std::sqrt(std::clamp((cf2 - a2) / gap, 0.0, 1.0));
  }
  const double bt = std::sqrt(bt1 * bt1 + bt2 * bt2);
  double beta1 = 1.0 / std::sqrt(2.0);
  double beta2 = 1.0 / std::sqrt(2.0);
  if (bt > 1.0e-12 * std::sqrt(rho * total2)) {
    beta1 = bt1 / bt;
    beta2 = bt2 / bt;
  }
  const double s = bn >= 0.0 ? 1.0 : -1.0;
  for (Vector7& row : right) {
    row.fill(0.0);
  }
  const auto setColumn = [&right](std::size_t column, const Vector7& values) {
    for (std::size_t row = 0; row < 7U; ++row) {
      right[row][column] = values[row];
    }
  };
  // Fast waves.
  setColumn(0, {rho * alpha_f, -alpha_f * cf, alpha_s * cs * beta1 * s,
                alpha_s * cs * beta2 * s, rho * a2 * alpha_f, alpha_s * sqrt_rho * a * beta1,
                alpha_s * sqrt_rho * a * beta2});
  setColumn(6, {rho * alpha_f, alpha_f * cf, -alpha_s * cs * beta1 * s,
                -alpha_s * cs * beta2 * s, rho * a2 * alpha_f, alpha_s * sqrt_rho * a * beta1,
                alpha_s * sqrt_rho * a * beta2});
  // Slow waves.
  setColumn(2, {rho * alpha_s, -alpha_s * cs, -alpha_f * cf * beta1 * s,
                -alpha_f * cf * beta2 * s, rho * a2 * alpha_s, -alpha_f * sqrt_rho * a * beta1,
                -alpha_f * sqrt_rho * a * beta2});
  setColumn(4, {rho * alpha_s, alpha_s * cs, alpha_f * cf * beta1 * s,
                alpha_f * cf * beta2 * s, rho * a2 * alpha_s, -alpha_f * sqrt_rho * a * beta1,
                -alpha_f * sqrt_rho * a * beta2});
  // Alfven waves (velocity and field perturbations perpendicular to B_t).
  setColumn(1, {0.0, 0.0, -beta2, beta1, 0.0, -s * sqrt_rho * beta2, s * sqrt_rho * beta1});
  setColumn(5, {0.0, 0.0, -beta2, beta1, 0.0, s * sqrt_rho * beta2, -s * sqrt_rho * beta1});
  // Entropy wave.
  setColumn(3, {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  return true;
}

[[nodiscard]] bool invertMatrix7(const Matrix7& input, Matrix7& inverse) noexcept {
  Matrix7 work = input;
  for (std::size_t row = 0; row < 7U; ++row) {
    inverse[row].fill(0.0);
    inverse[row][row] = 1.0;
  }
  double scale = 0.0;
  for (const Vector7& row : input) {
    for (const double value : row) {
      scale = std::max(scale, std::abs(value));
    }
  }
  if (!(scale > 0.0) || !std::isfinite(scale)) {
    return false;
  }
  for (std::size_t column = 0; column < 7U; ++column) {
    std::size_t pivot = column;
    for (std::size_t row = column + 1U; row < 7U; ++row) {
      if (std::abs(work[row][column]) > std::abs(work[pivot][column])) {
        pivot = row;
      }
    }
    if (!(std::abs(work[pivot][column]) > 1.0e-13 * scale)) {
      return false;
    }
    std::swap(work[pivot], work[column]);
    std::swap(inverse[pivot], inverse[column]);
    const double inverse_pivot = 1.0 / work[column][column];
    for (std::size_t k = 0; k < 7U; ++k) {
      work[column][k] *= inverse_pivot;
      inverse[column][k] *= inverse_pivot;
    }
    for (std::size_t row = 0; row < 7U; ++row) {
      if (row == column) {
        continue;
      }
      const double factor = work[row][column];
      if (factor == 0.0) {
        continue;
      }
      for (std::size_t k = 0; k < 7U; ++k) {
        work[row][k] -= factor * work[column][k];
        inverse[row][k] -= factor * inverse[column][k];
      }
    }
  }
  return true;
}

[[nodiscard]] Vector7 multiply(const Matrix7& matrix, const Vector7& vector) noexcept {
  Vector7 result{};
  for (std::size_t row = 0; row < 7U; ++row) {
    double sum = 0.0;
    for (std::size_t column = 0; column < 7U; ++column) {
      sum += matrix[row][column] * vector[column];
    }
    result[row] = sum;
  }
  return result;
}

// Colella & Woodward (1984) monotonicity constraints on one parabola given
// the cell value and limited interface values (all in one characteristic or
// scalar variable). Interface values are first confined between the adjacent
// cell means (CW84 eq. 1.10), then the parabola is made monotone.
void limitParabola(double c_minus, double c, double c_plus, double& lower,
                   double& upper) noexcept {
  lower = std::clamp(lower, std::min(c_minus, c), std::max(c_minus, c));
  upper = std::clamp(upper, std::min(c, c_plus), std::max(c, c_plus));
  if ((upper - c) * (c - lower) <= 0.0) {
    lower = c;
    upper = c;
    return;
  }
  const double range = upper - lower;
  const double center = c - 0.5 * (lower + upper);
  if (range * center > range * range / 6.0) {
    lower = 3.0 * c - 2.0 * upper;
  } else if (-range * range / 6.0 > range * center) {
    upper = 3.0 * c - 2.0 * lower;
  }
}

[[nodiscard]] Vector7 toLocal7(const MhdPrimitiveState& state, MhdAxis axis) noexcept {
  const std::array<double, 3> v = permuteToLocal(
      axis, {state.vel_x_peculiar, state.vel_y_peculiar, state.vel_z_peculiar});
  const std::array<double, 3> b =
      permuteToLocal(axis, {state.bx_comoving, state.by_comoving, state.bz_comoving});
  return {state.rho_comoving, v[0], v[1], v[2], state.pressure_comoving, b[1], b[2]};
}

[[nodiscard]] double normalField(const MhdPrimitiveState& state, MhdAxis axis) noexcept {
  return permuteToLocal(axis, {state.bx_comoving, state.by_comoving, state.bz_comoving})[0];
}

[[nodiscard]] MhdPrimitiveState fromLocal7(const Vector7& w, double bn, double metal,
                                           MhdAxis axis) noexcept {
  const std::array<double, 3> v = permuteToGlobal(axis, {w[1], w[2], w[3]});
  const std::array<double, 3> b = permuteToGlobal(axis, {bn, w[5], w[6]});
  MhdPrimitiveState state{};
  state.rho_comoving = w[0];
  state.vel_x_peculiar = v[0];
  state.vel_y_peculiar = v[1];
  state.vel_z_peculiar = v[2];
  state.pressure_comoving = w[4];
  state.bx_comoving = b[0];
  state.by_comoving = b[1];
  state.bz_comoving = b[2];
  state.metallicity_mass_fraction = metal;
  return state;
}

[[nodiscard]] MhdPrimitiveState plmCell(const MhdPrimitiveState& minus,
                                        const MhdPrimitiveState& center,
                                        const MhdPrimitiveState& plus,
                                        bool upper) noexcept {
  const auto lane = [upper](double m, double c, double p) {
    const double slope = mcLimitedSlope(m, c, p);
    return c + (upper ? 0.5 : -0.5) * slope;
  };
  MhdPrimitiveState out{};
  out.rho_comoving = lane(minus.rho_comoving, center.rho_comoving, plus.rho_comoving);
  out.vel_x_peculiar = lane(minus.vel_x_peculiar, center.vel_x_peculiar, plus.vel_x_peculiar);
  out.vel_y_peculiar = lane(minus.vel_y_peculiar, center.vel_y_peculiar, plus.vel_y_peculiar);
  out.vel_z_peculiar = lane(minus.vel_z_peculiar, center.vel_z_peculiar, plus.vel_z_peculiar);
  out.pressure_comoving =
      lane(minus.pressure_comoving, center.pressure_comoving, plus.pressure_comoving);
  out.bx_comoving = lane(minus.bx_comoving, center.bx_comoving, plus.bx_comoving);
  out.by_comoving = lane(minus.by_comoving, center.by_comoving, plus.by_comoving);
  out.bz_comoving = lane(minus.bz_comoving, center.bz_comoving, plus.bz_comoving);
  out.metallicity_mass_fraction =
      lane(minus.metallicity_mass_fraction, center.metallicity_mass_fraction,
           plus.metallicity_mass_fraction);
  return out;
}

[[nodiscard]] bool admissiblePrimitive(const MhdPrimitiveState& state) noexcept {
  return std::isfinite(state.rho_comoving) && state.rho_comoving > 0.0 &&
      std::isfinite(state.pressure_comoving) && state.pressure_comoving > 0.0 &&
      std::isfinite(state.vel_x_peculiar) && std::isfinite(state.vel_y_peculiar) &&
      std::isfinite(state.vel_z_peculiar) && std::isfinite(state.bx_comoving) &&
      std::isfinite(state.by_comoving) && std::isfinite(state.bz_comoving) &&
      std::isfinite(state.metallicity_mass_fraction) &&
      state.metallicity_mass_fraction >= 0.0 && state.metallicity_mass_fraction <= 1.0;
}

}  // namespace

std::string_view mhdReconstructionName(MhdReconstructionMode mode) noexcept {
  switch (mode) {
    case MhdReconstructionMode::kPlmMc: return "plm_mc";
    case MhdReconstructionMode::kPpm: return "ppm";
  }
  return "unknown";
}

// ---------------------------------------------------------------------------
// MhdGasStateSoa
// ---------------------------------------------------------------------------

MhdGasStateSoa::MhdGasStateSoa(std::size_t cell_count) {
  resize(cell_count);
}

void MhdGasStateSoa::resize(std::size_t cell_count) {
  for (MhdDoubleArray* lane : {&m_rho, &m_mx, &m_my, &m_mz, &m_energy, &m_metal}) {
    if (lane->size() != cell_count) {
      lane->resize(cell_count, 0.0);
    }
  }
}

std::uint64_t MhdGasStateSoa::ownedCapacityBytes() const noexcept {
  return m_rho.ownedCapacityBytes() + m_mx.ownedCapacityBytes() +
      m_my.ownedCapacityBytes() + m_mz.ownedCapacityBytes() +
      m_energy.ownedCapacityBytes() + m_metal.ownedCapacityBytes();
}

std::span<double> MhdGasStateSoa::lane(std::size_t index) noexcept {
  switch (index) {
    case 0: return m_rho.span();
    case 1: return m_mx.span();
    case 2: return m_my.span();
    case 3: return m_mz.span();
    case 4: return m_energy.span();
    default: return m_metal.span();
  }
}

std::span<const double> MhdGasStateSoa::lane(std::size_t index) const noexcept {
  switch (index) {
    case 0: return m_rho.span();
    case 1: return m_mx.span();
    case 2: return m_my.span();
    case 3: return m_mz.span();
    case 4: return m_energy.span();
    default: return m_metal.span();
  }
}

MhdConservedState MhdGasStateSoa::loadGas(std::size_t cell) const noexcept {
  MhdConservedState state{};
  state.rho_comoving = m_rho[cell];
  state.momentum_density_x_comoving = m_mx[cell];
  state.momentum_density_y_comoving = m_my[cell];
  state.momentum_density_z_comoving = m_mz[cell];
  state.total_energy_density_comoving = m_energy[cell];
  state.metal_mass_density_comoving = m_metal[cell];
  return state;
}

void MhdGasStateSoa::storeGas(std::size_t cell, const MhdConservedState& state) noexcept {
  m_rho[cell] = state.rho_comoving;
  m_mx[cell] = state.momentum_density_x_comoving;
  m_my[cell] = state.momentum_density_y_comoving;
  m_mz[cell] = state.momentum_density_z_comoving;
  m_energy[cell] = state.total_energy_density_comoving;
  m_metal[cell] = state.metal_mass_density_comoving;
}

void MhdGasStateSoa::copyFrom(const MhdGasStateSoa& other) {
  resize(other.size());
  for (std::size_t lane_index = 0; lane_index < k_lane_count; ++lane_index) {
    const std::span<const double> source = other.lane(lane_index);
    std::copy(source.begin(), source.end(), lane(lane_index).begin());
  }
}

void MhdGasStateSoa::assignLinearCombination(double alpha, const MhdGasStateSoa& lhs,
                                             double beta, const MhdGasStateSoa& rhs) {
  if (lhs.size() != rhs.size()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_gas_state", "linear combination of mismatched states"));
  }
  resize(lhs.size());
  for (std::size_t lane_index = 0; lane_index < k_lane_count; ++lane_index) {
    const std::span<const double> a = lhs.lane(lane_index);
    const std::span<const double> b = rhs.lane(lane_index);
    std::span<double> out = lane(lane_index);
    for (std::size_t cell = 0; cell < out.size(); ++cell) {
      out[cell] = alpha * a[cell] + beta * b[cell];
    }
  }
}

// ---------------------------------------------------------------------------
// MhdBoundaryFaceFluxes / counters / workspace
// ---------------------------------------------------------------------------

void MhdBoundaryFaceFluxes::resize(const MhdPatchGeometry& geometry) {
  nx = geometry.nx;
  ny = geometry.ny;
  nz = geometry.nz;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t count = faceCount(axis) * k_values_per_face;
    for (std::size_t side = 0; side < 2U; ++side) {
      if (values[axis][side].size() != count) {
        values[axis][side].resize(count, 0.0);
      }
    }
  }
}

std::size_t MhdBoundaryFaceFluxes::faceCount(std::size_t axis) const noexcept {
  const std::array<std::size_t, 3> n{nx, ny, nz};
  return n[(axis + 1U) % 3U] * n[(axis + 2U) % 3U];
}

double* MhdBoundaryFaceFluxes::at(std::size_t axis, std::size_t side, std::size_t t1,
                                  std::size_t t2) noexcept {
  const std::array<std::size_t, 3> n{nx, ny, nz};
  const std::size_t n_t1 = n[(axis + 1U) % 3U];
  return values[axis][side].data() + (t1 + n_t1 * t2) * k_values_per_face;
}

const double* MhdBoundaryFaceFluxes::at(std::size_t axis, std::size_t side, std::size_t t1,
                                        std::size_t t2) const noexcept {
  const std::array<std::size_t, 3> n{nx, ny, nz};
  const std::size_t n_t1 = n[(axis + 1U) % 3U];
  return values[axis][side].data() + (t1 + n_t1 * t2) * k_values_per_face;
}

std::uint64_t MhdBoundaryFaceFluxes::ownedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const auto& sides : values) {
    for (const MhdDoubleArray& lane : sides) {
      bytes += lane.ownedCapacityBytes();
    }
  }
  return bytes;
}

MhdStageCounters& MhdStageCounters::operator+=(const MhdStageCounters& rhs) noexcept {
  faces_solved += rhs.faces_solved;
  hlle_fallback_faces += rhs.hlle_fallback_faces;
  first_order_faces += rhs.first_order_faces;
  ppm_characteristic_fallback_cells += rhs.ppm_characteristic_fallback_cells;
  edges_built += rhs.edges_built;
  return *this;
}

void MhdStageWorkspace::resize(const MhdPatchGeometry& geometry,
                               MhdReconstructionMode reconstruction) {
  const std::size_t ghost = mhdGhostCells(reconstruction);
  if (primitives.nx != geometry.nx || primitives.ny != geometry.ny ||
      primitives.nz != geometry.nz || primitives.ghost != ghost) {
    primitives.resize(geometry.nx, geometry.ny, geometry.nz, ghost);
  }
  if (staggered.nx != geometry.nx || staggered.ny != geometry.ny ||
      staggered.nz != geometry.nz ||
      staggered.tangential_ghost != k_mhd_staggered_tangential_ghost) {
    staggered.resize(geometry.nx, geometry.ny, geometry.nz, k_mhd_staggered_tangential_ghost);
  }
  if (signals.nx != geometry.nx || signals.ny != geometry.ny || signals.nz != geometry.nz) {
    signals.resize(geometry);
  }
  if (emf.geometry().nx != geometry.nx || emf.geometry().ny != geometry.ny ||
      emf.geometry().nz != geometry.nz) {
    emf.resize(geometry);
  }
  if (boundary_fluxes.nx != geometry.nx || boundary_fluxes.ny != geometry.ny ||
      boundary_fluxes.nz != geometry.nz) {
    boundary_fluxes.resize(geometry);
  }
}

std::uint64_t MhdStageWorkspace::ownedCapacityBytes() const noexcept {
  return primitives.ownedCapacityBytes() + staggered.ownedCapacityBytes() +
      signals.ownedCapacityBytes() + emf.ownedCapacityBytes() +
      boundary_fluxes.ownedCapacityBytes();
}

// ---------------------------------------------------------------------------
// MhdStepOperator
// ---------------------------------------------------------------------------

MhdStepOperator::MhdStepOperator(MhdNumericsPolicy policy) : m_policy(policy) {
  if (!(m_policy.gamma_adiabatic > 1.0) || !std::isfinite(m_policy.gamma_adiabatic)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_step_operator", "gamma_adiabatic must exceed 1"));
  }
  if (!(m_policy.cfl_ideal > 0.0) || !(m_policy.cfl_ideal < 1.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_step_operator", "ideal CFL must lie in (0,1)"));
  }
  if (m_policy.max_step_retries < 0 || m_policy.max_step_retries > 8) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_step_operator", "max_step_retries must lie in [0,8]"));
  }
}

std::pair<double, double> MhdStepOperator::reconstructPlmMc(double c_minus, double c,
                                                            double c_plus) noexcept {
  const double slope = mcLimitedSlope(c_minus, c, c_plus);
  return {c - 0.5 * slope, c + 0.5 * slope};
}

std::pair<double, double> MhdStepOperator::reconstructPpmScalar(
    double c_minus2, double c_minus, double c, double c_plus, double c_plus2) noexcept {
  double lower = (7.0 / 12.0) * (c_minus + c) - (1.0 / 12.0) * (c_minus2 + c_plus);
  double upper = (7.0 / 12.0) * (c + c_plus) - (1.0 / 12.0) * (c_minus + c_plus2);
  limitParabola(c_minus, c, c_plus, lower, upper);
  return {lower, upper};
}

void MhdStepOperator::fillGhostPrimitives(const MhdPatchBoundary& boundary,
                                          MhdPrimitivePatch& primitives) {
  const std::array<std::ptrdiff_t, 3> extents{static_cast<std::ptrdiff_t>(primitives.nx),
                                              static_cast<std::ptrdiff_t>(primitives.ny),
                                              static_cast<std::ptrdiff_t>(primitives.nz)};
  const auto g = static_cast<std::ptrdiff_t>(primitives.ghost);
  std::array<MhdDoubleArray*, 3> velocity{&primitives.vx, &primitives.vy, &primitives.vz};
  std::array<MhdDoubleArray*, 3> field{&primitives.bx, &primitives.by, &primitives.bz};
  // Sequential axis passes over the full padded range: a corner ghost is
  // finalized by the last axis pass, which reads values already produced by the
  // earlier passes (or by the level stepper for neighbour sides).
  for (std::size_t t = 0; t < 3U; ++t) {
    const std::size_t o1 = (t + 1U) % 3U;
    const std::size_t o2 = (t + 2U) % 3U;
    for (std::ptrdiff_t c2 = -g; c2 < extents[o2] + g; ++c2) {
      for (std::ptrdiff_t c1 = -g; c1 < extents[o1] + g; ++c1) {
        for (std::ptrdiff_t ct = -g; ct < extents[t] + g; ++ct) {
          if (ct >= 0 && ct < extents[t]) {
            continue;
          }
          const bool lower_side = ct < 0;
          const MhdBoundaryKind kind = lower_side ? boundary.lower[t] : boundary.upper[t];
          if (kind == MhdBoundaryKind::kNeighborPatch || kind == MhdBoundaryKind::kImportedMpi) {
            continue;
          }
          std::ptrdiff_t source = ct;
          bool wall = false;
          switch (kind) {
            case MhdBoundaryKind::kPeriodic:
              source = ((ct % extents[t]) + extents[t]) % extents[t];
              break;
            case MhdBoundaryKind::kOpen:
              source = lower_side ? 0 : extents[t] - 1;
              break;
            case MhdBoundaryKind::kReflectingConductingWall:
              source = lower_side ? (-1 - ct) : (2 * extents[t] - 1 - ct);
              wall = true;
              break;
            default:
              break;
          }
          MhdIndex3 target{};
          target[t] = ct;
          target[o1] = c1;
          target[o2] = c2;
          MhdIndex3 origin = target;
          origin[t] = source;
          const std::size_t to = primitives.index(target);
          const std::size_t from = primitives.index(origin);
          for (MhdDoubleArray* lane : {&primitives.rho, &primitives.vx, &primitives.vy,
                                       &primitives.vz, &primitives.pressure, &primitives.bx,
                                       &primitives.by, &primitives.bz, &primitives.metal}) {
            (*lane)[to] = (*lane)[from];
          }
          if (wall) {
            // Wall-normal velocity odd (zero mass/energy flux), wall-normal
            // field even (wall flux preserved), wall-tangential field odd
            // (zero Poynting/tangential stress flux). Tangential E at the wall
            // is enforced exactly by MhdEmfBuilder::applyEdgeBoundaryRules.
            (*velocity[t])[to] = -(*velocity[t])[to];
            (*field[o1])[to] = -(*field[o1])[to];
            (*field[o2])[to] = -(*field[o2])[to];
          }
        }
      }
    }
  }
}

bool MhdStepOperator::recoverCellPrimitive(const MhdConservedState& gas,
                                           const std::array<double, 3>& b_comoving,
                                           double scale_factor,
                                           MhdPrimitiveState& out) const noexcept {
  return recoverComoving(gas, b_comoving, scale_factor, m_policy.gamma_adiabatic, out);
}

std::array<double, 3> MhdStepOperator::cellFieldFromFaces(const MhdPatchGeometry& geometry,
                                                          std::span<const double> phi_x,
                                                          std::span<const double> phi_y,
                                                          std::span<const double> phi_z,
                                                          std::size_t i, std::size_t j,
                                                          std::size_t k) noexcept {
  return {0.5 * (phi_x[geometry.faceIndexX(i, j, k)] + phi_x[geometry.faceIndexX(i + 1U, j, k)]) /
              geometry.faceAreaX(),
          0.5 * (phi_y[geometry.faceIndexY(i, j, k)] + phi_y[geometry.faceIndexY(i, j + 1U, k)]) /
              geometry.faceAreaY(),
          0.5 * (phi_z[geometry.faceIndexZ(i, j, k)] + phi_z[geometry.faceIndexZ(i, j, k + 1U)]) /
              geometry.faceAreaZ()};
}

void MhdStepOperator::fillPhysicalGhosts(const MhdPatchBoundary& boundary,
                                         MhdStageWorkspace& workspace) {
  fillGhostPrimitives(boundary, workspace.primitives);
  fillStaggeredGhosts(boundary, workspace.staggered);
}

std::size_t MhdStepOperator::recoverPrimitives(const MhdPatchGeometry& geometry,
                                               const MhdGasStateSoa& gas,
                                               const MhdFaceFluxView& faces,
                                               double scale_factor,
                                               MhdStageWorkspace& workspace) const {
  if (!faces.valid() || gas.size() != geometry.cellCount()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_primitive_recovery",
                               "gas or face-flux lanes do not match the patch geometry"));
  }
  if (!(scale_factor > 0.0) || !std::isfinite(scale_factor)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_primitive_recovery", "non-positive scale factor"));
  }
  workspace.resize(geometry, m_policy.reconstruction);
  std::size_t bad = 0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t cell = geometry.cellIndex(i, j, k);
        MhdPrimitiveState primitive{};
        if (!recoverComoving(gas.loadGas(cell), cellFieldComoving(geometry, faces, i, j, k),
                             scale_factor, m_policy.gamma_adiabatic, primitive)) {
          ++bad;
          continue;
        }
        workspace.primitives.store({static_cast<std::ptrdiff_t>(i),
                                    static_cast<std::ptrdiff_t>(j),
                                    static_cast<std::ptrdiff_t>(k)},
                                   primitive);
      }
    }
  }
  const double inverse_sqrt_a = 1.0 / std::sqrt(scale_factor);
  const std::array<std::size_t, 3> n{geometry.nx, geometry.ny, geometry.nz};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::span<const double> lane = faceLane(faces, axis);
    const double inverse_area = 1.0 / faceArea(geometry, axis);
    const std::size_t t1 = (axis + 1U) % 3U;
    const std::size_t t2 = (axis + 2U) % 3U;
    for (std::size_t b = 0; b < n[t2]; ++b) {
      for (std::size_t a = 0; a < n[t1]; ++a) {
        for (std::size_t f = 0; f <= n[axis]; ++f) {
          MhdIndex3 c{};
          c[axis] = static_cast<std::ptrdiff_t>(f);
          c[t1] = static_cast<std::ptrdiff_t>(a);
          c[t2] = static_cast<std::ptrdiff_t>(b);
          const double value = lane[faceIndex(geometry, axis, c)];
          if (!std::isfinite(value)) {
            throw MhdError(MhdErrorCode::kInvalidField,
                           makeContext("mhd_primitive_recovery", "non-finite face flux"));
          }
          workspace.staggered.at(static_cast<MhdAxis>(axis), c) =
              value * inverse_area * inverse_sqrt_a;
        }
      }
    }
  }
  return bad;
}

bool MhdStepOperator::reconstructCell(const MhdPrimitivePatch& primitives,
                                      const MhdIndex3& cell,
                                      MhdAxis axis,
                                      MhdCellFaceStates& states,
                                      std::uint64_t* characteristic_fallbacks) const {
  const std::size_t along = mhdAxisIndex(axis);
  const auto neighbour = [&](std::ptrdiff_t offset) {
    MhdIndex3 c = cell;
    c[along] += offset;
    return primitives.load(c);
  };
  const MhdPrimitiveState center = neighbour(0);
  if (!admissiblePrimitive(center)) {
    return false;
  }
  const MhdPrimitiveState minus = neighbour(-1);
  const MhdPrimitiveState plus = neighbour(+1);
  const auto plm = [&]() {
    states.lower = plmCell(minus, center, plus, false);
    states.upper = plmCell(minus, center, plus, true);
    return admissiblePrimitive(states.lower) && admissiblePrimitive(states.upper);
  };
  if (m_policy.reconstruction == MhdReconstructionMode::kPlmMc) {
    return plm();
  }

  // Bounded characteristic PPM.
  const MhdPrimitiveState minus2 = neighbour(-2);
  const MhdPrimitiveState plus2 = neighbour(+2);
  const std::array<Vector7, 5> w{toLocal7(minus2, axis), toLocal7(minus, axis),
                                 toLocal7(center, axis), toLocal7(plus, axis),
                                 toLocal7(plus2, axis)};
  const double bn = normalField(center, axis);
  Matrix7 right{};
  Matrix7 left{};
  const bool basis_ok = buildRightEigenvectors(w[2], bn, m_policy.gamma_adiabatic, right) &&
      invertMatrix7(right, left);
  if (!basis_ok) {
    if (characteristic_fallbacks != nullptr) {
      ++(*characteristic_fallbacks);
    }
    return plm();
  }
  Vector7 lower_face{};
  Vector7 upper_face{};
  for (std::size_t m = 0; m < 7U; ++m) {
    lower_face[m] = (7.0 / 12.0) * (w[1][m] + w[2][m]) - (1.0 / 12.0) * (w[0][m] + w[3][m]);
    upper_face[m] = (7.0 / 12.0) * (w[2][m] + w[3][m]) - (1.0 / 12.0) * (w[1][m] + w[4][m]);
  }
  const Vector7 char_minus = multiply(left, w[1]);
  const Vector7 char_center = multiply(left, w[2]);
  const Vector7 char_plus = multiply(left, w[3]);
  Vector7 char_lower = multiply(left, lower_face);
  Vector7 char_upper = multiply(left, upper_face);
  Vector7 delta_lower{};
  Vector7 delta_upper{};
  for (std::size_t m = 0; m < 7U; ++m) {
    limitParabola(char_minus[m], char_center[m], char_plus[m], char_lower[m], char_upper[m]);
    delta_lower[m] = char_lower[m] - char_center[m];
    delta_upper[m] = char_upper[m] - char_center[m];
  }
  const Vector7 back_lower = multiply(right, delta_lower);
  const Vector7 back_upper = multiply(right, delta_upper);
  Vector7 primitive_lower{};
  Vector7 primitive_upper{};
  for (std::size_t m = 0; m < 7U; ++m) {
    primitive_lower[m] = w[2][m] + back_lower[m];
    primitive_upper[m] = w[2][m] + back_upper[m];
  }
  // Passive metal fraction: scalar bounded PPM (advected with the entropy
  // wave); limitParabola keeps it within the neighbouring cell means.
  const std::pair<double, double> metal = reconstructPpmScalar(
      minus2.metallicity_mass_fraction, minus.metallicity_mass_fraction,
      center.metallicity_mass_fraction, plus.metallicity_mass_fraction,
      plus2.metallicity_mass_fraction);
  states.lower = fromLocal7(primitive_lower, bn, metal.first, axis);
  states.upper = fromLocal7(primitive_upper, bn, metal.second, axis);
  if (!admissiblePrimitive(states.lower) || !admissiblePrimitive(states.upper)) {
    // Positivity fallback for this cell only; counted, never silent.
    if (characteristic_fallbacks != nullptr) {
      ++(*characteristic_fallbacks);
    }
    return plm();
  }
  return true;
}

void MhdStepOperator::evaluateIdealStage(const MhdPatchGeometry& geometry,
                                         const MhdPatchBoundary& boundary,
                                         double a_stage,
                                         double hubble_stage,
                                         bool apply_expansion_source,
                                         double flux_weight,
                                         MhdStageWorkspace& workspace,
                                         MhdGasStateSoa& gas_target,
                                         MhdStageCounters& counters) const {
  if (!(a_stage > 0.0) || !std::isfinite(a_stage) || !std::isfinite(hubble_stage) ||
      !std::isfinite(flux_weight)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_ideal", "invalid stage scale factor, H or weight"));
  }
  if (gas_target.size() != geometry.cellCount()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_ideal", "gas target does not match the patch"));
  }
  const MhdPrimitivePatch& primitives = workspace.primitives;
  const std::array<std::ptrdiff_t, 3> n{static_cast<std::ptrdiff_t>(geometry.nx),
                                        static_cast<std::ptrdiff_t>(geometry.ny),
                                        static_cast<std::ptrdiff_t>(geometry.nz)};
  const double gamma = m_policy.gamma_adiabatic;
  std::array<std::span<double>, MhdGasStateSoa::k_lane_count> lanes{};
  for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
    lanes[lane] = gas_target.lane(lane);
  }
  const auto cellIndex = [&](const MhdIndex3& c) {
    return geometry.cellIndex(static_cast<std::size_t>(c[0]), static_cast<std::size_t>(c[1]),
                              static_cast<std::size_t>(c[2]));
  };

  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const MhdAxis normal = static_cast<MhdAxis>(axis);
    const std::size_t t1 = (axis + 1U) % 3U;
    const std::size_t t2 = (axis + 2U) % 3U;
    const double weight_per_length = flux_weight / (a_stage * spacingAlong(geometry, axis));
    std::span<double> speed_plus = workspace.signals.speed_plus[axis].span();
    std::span<double> speed_minus = workspace.signals.speed_minus[axis].span();
    for (std::ptrdiff_t b = -1; b <= n[t2]; ++b) {
      for (std::ptrdiff_t a = -1; a <= n[t1]; ++a) {
        const bool interior_line = a >= 0 && a < n[t1] && b >= 0 && b < n[t2];
        MhdIndex3 cell{};
        cell[t1] = a;
        cell[t2] = b;
        cell[axis] = -1;
        MhdCellFaceStates previous{};
        if (!reconstructCell(primitives, cell, normal, previous,
                             &counters.ppm_characteristic_fallback_cells)) {
          throw MhdError(MhdErrorCode::kNonPositiveState,
                         makeContext("mhd_ideal", "inadmissible ghost cell before the sweep"));
        }
        for (std::ptrdiff_t f = 0; f <= n[axis]; ++f) {
          cell[axis] = f;
          // Cell f (a ghost cell when f == n) supplies the right state of face f.
          MhdCellFaceStates current{};
          if (!reconstructCell(primitives, cell, normal, current,
                               &counters.ppm_characteristic_fallback_cells)) {
            throw MhdError(MhdErrorCode::kNonPositiveState,
                           makeContext("mhd_ideal", "inadmissible cell in the face sweep"));
          }
          MhdIndex3 face = cell;  // face f lies between cells f-1 and f
          const double bn = workspace.staggered.at(normal, face);
          const std::size_t signal_index = workspace.signals.index(normal, face);
          const MhdPrimitiveState& left = previous.upper;
          const MhdPrimitiveState& right = current.lower;
          if (!interior_line) {
            // Ghost-row face: only the signed fan is needed for UCT-HLL.
            const MhdSignalBounds bounds =
                MhdRiemann::signalBounds(normal, left, right, bn, gamma);
            if (!bounds.valid) {
              throw MhdError(MhdErrorCode::kNonPositiveState,
                             makeContext("mhd_ideal", "invalid signal fan on a ghost-row face"));
            }
            speed_plus[signal_index] = std::max(0.0, bounds.speed_right);
            speed_minus[signal_index] = std::max(0.0, -bounds.speed_left);
            previous = current;
            continue;
          }
          MhdRiemann::Result solved =
              MhdRiemann::solve(normal, left, right, bn, gamma, m_policy.riemann);
          double metal_upwind = 0.0;
          if (!solved.valid) {
            MhdIndex3 left_cell = cell;
            left_cell[axis] = f - 1;
            const MhdPrimitiveState left_average = primitives.load(left_cell);
            const MhdPrimitiveState right_average = primitives.load(cell);
            solved = firstOrderPositivityFlux(normal, left_average, right_average, bn, gamma);
            if (!solved.valid) {
              std::ostringstream stream;
              stream << "no admissible flux (HLLD, HLLE and first-order cell-average "
                        "Rusanov all failed) on axis " << axis << " face (" << face[0] << ','
                     << face[1] << ',' << face[2] << ')';
              throw MhdError(MhdErrorCode::kNonPositiveState,
                             makeContext("mhd_ideal", stream.str()));
            }
            metal_upwind = solved.flux.mass_flux >= 0.0
                ? left_average.metallicity_mass_fraction
                : right_average.metallicity_mass_fraction;
            ++counters.first_order_faces;
          } else {
            metal_upwind = solved.flux.mass_flux >= 0.0 ? left.metallicity_mass_fraction
                                                        : right.metallicity_mass_fraction;
            if (solved.used_hlle && m_policy.riemann == MhdRiemann::Solver::kHlld) {
              ++counters.hlle_fallback_faces;
            }
          }
          ++counters.faces_solved;
          const MhdRiemannFlux& flux = solved.flux;
          speed_plus[signal_index] = std::max(0.0, flux.fastest_signal_speed);
          speed_minus[signal_index] = std::max(0.0, -flux.slowest_signal_speed);
          const std::array<double, MhdGasStateSoa::k_lane_count> values{
              flux.mass_flux, flux.momentum_flux[0], flux.momentum_flux[1],
              flux.momentum_flux[2], flux.energy_flux, flux.mass_flux * metal_upwind};
          if (f == 0 || f == n[axis]) {
            double* record = workspace.boundary_fluxes.at(
                axis, f == 0 ? 0U : 1U, static_cast<std::size_t>(a), static_cast<std::size_t>(b));
            std::copy(values.begin(), values.end(), record);
          } else {
            // Positive normal flux leaves the lower cell and enters the upper.
            MhdIndex3 lower_cell = cell;
            lower_cell[axis] = f - 1;
            const std::size_t lower_index = cellIndex(lower_cell);
            const std::size_t upper_index = cellIndex(cell);
            for (std::size_t lane = 0; lane < values.size(); ++lane) {
              lanes[lane][lower_index] -= weight_per_length * values[lane];
              lanes[lane][upper_index] += weight_per_length * values[lane];
            }
          }
          previous = current;
        }
      }
    }
  }

  // Comoving expansion source at the stage state (CHUI-MHD-REF v1.0 4.3):
  //   dm_c/dt = -H m_c,  dE_c/dt = -H (2 K_c + 3 p_c + M_c),  M_c = B_eff^2/2.
  if (apply_expansion_source && hubble_stage != 0.0) {
    for (std::ptrdiff_t k = 0; k < n[2]; ++k) {
      for (std::ptrdiff_t j = 0; j < n[1]; ++j) {
        for (std::ptrdiff_t i = 0; i < n[0]; ++i) {
          const MhdIndex3 c{i, j, k};
          const MhdPrimitiveState state = primitives.load(c);
          const std::size_t index = cellIndex(c);
          const double rho = state.rho_comoving;
          const double kinetic = 0.5 * rho *
              (state.vel_x_peculiar * state.vel_x_peculiar +
               state.vel_y_peculiar * state.vel_y_peculiar +
               state.vel_z_peculiar * state.vel_z_peculiar);
          const double magnetic = state.magneticEnergyDensity();
          const double rate = -hubble_stage * flux_weight;
          lanes[1][index] += rate * rho * state.vel_x_peculiar;
          lanes[2][index] += rate * rho * state.vel_y_peculiar;
          lanes[3][index] += rate * rho * state.vel_z_peculiar;
          lanes[4][index] +=
              rate * (2.0 * kinetic + 3.0 * state.pressure_comoving + magnetic);
        }
      }
    }
  }

  MhdEmfBuilder::Diagnostics emf_diagnostics{};
  MhdEmfBuilder::buildUctHll(geometry, boundary, primitives, workspace.staggered,
                             workspace.signals, a_stage, workspace.emf, &emf_diagnostics);
  counters.edges_built += emf_diagnostics.edges_built;
}

void MhdStepOperator::applyBoundaryFaceFluxes(const MhdPatchGeometry& geometry,
                                              const MhdBoundaryFaceFluxes& fluxes,
                                              double a_stage,
                                              double flux_weight,
                                              MhdGasStateSoa& gas_target) {
  const std::array<std::size_t, 3> n{geometry.nx, geometry.ny, geometry.nz};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t t1 = (axis + 1U) % 3U;
    const std::size_t t2 = (axis + 2U) % 3U;
    const double weight_per_length = flux_weight / (a_stage * spacingAlong(geometry, axis));
    for (std::size_t side = 0; side < 2U; ++side) {
      const double sign = side == 0U ? 1.0 : -1.0;
      for (std::size_t b = 0; b < n[t2]; ++b) {
        for (std::size_t a = 0; a < n[t1]; ++a) {
          std::array<std::size_t, 3> c{};
          c[axis] = side == 0U ? 0U : n[axis] - 1U;
          c[t1] = a;
          c[t2] = b;
          const std::size_t cell = geometry.cellIndex(c[0], c[1], c[2]);
          const double* record = fluxes.at(axis, side, a, b);
          for (std::size_t lane = 0; lane < MhdGasStateSoa::k_lane_count; ++lane) {
            gas_target.lane(lane)[cell] += sign * weight_per_length * record[lane];
          }
        }
      }
    }
  }
}

double MhdStepOperator::estimateIdealTimeStep(const MhdPatchGeometry& geometry,
                                              const MhdGasStateSoa& gas,
                                              const MhdFaceFluxView& faces,
                                              double scale_factor,
                                              double cfl) const {
  if (!faces.valid() || gas.size() != geometry.cellCount() || !(scale_factor > 0.0) ||
      !(cfl > 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_ideal_timestep", "invalid timestep inputs"));
  }
  const std::array<double, 3> inverse_spacing{
      1.0 / (scale_factor * geometry.dx_comoving), 1.0 / (scale_factor * geometry.dy_comoving),
      1.0 / (scale_factor * geometry.dz_comoving)};
  double maximum_rate = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t cell = geometry.cellIndex(i, j, k);
        MhdPrimitiveState primitive{};
        if (!recoverComoving(gas.loadGas(cell), cellFieldComoving(geometry, faces, i, j, k),
                             scale_factor, m_policy.gamma_adiabatic, primitive)) {
          MhdFailureContext failure =
              makeContext("mhd_ideal_timestep", "inadmissible cell in the CFL estimate");
          failure.patch_id = faces.patch_id;
          failure.geometry_generation = faces.geometry_generation;
          failure.scale_factor = scale_factor;
          throw MhdError(MhdErrorCode::kNonPositiveState, std::move(failure));
        }
        const std::array<double, 3> velocity{primitive.vel_x_peculiar, primitive.vel_y_peculiar,
                                             primitive.vel_z_peculiar};
        double rate = 0.0;
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          rate += (std::abs(velocity[axis]) +
                   fastMagnetosonicSpeed(primitive, m_policy.gamma_adiabatic,
                                         static_cast<MhdAxis>(axis))) *
              inverse_spacing[axis];
        }
        maximum_rate = std::max(maximum_rate, rate);
      }
    }
  }
  if (!(maximum_rate > 0.0) || !std::isfinite(maximum_rate)) {
    throw MhdError(MhdErrorCode::kNonPositiveState,
                   makeContext("mhd_ideal_timestep", "non-positive or non-finite signal rate"));
  }
  return cfl / maximum_rate;
}

MhdConservationBudget MhdStepOperator::patchBudget(const MhdPatchGeometry& geometry,
                                                   const MhdGasStateSoa& gas,
                                                   const MhdFaceFluxView& faces,
                                                   double scale_factor) const {
  MhdConservationBudget budget{};
  if (!faces.valid() || gas.size() != geometry.cellCount() || !(scale_factor > 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_budget", "invalid budget inputs"));
  }
  const double volume = geometry.cellVolume();
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t cell = geometry.cellIndex(i, j, k);
        const MhdConservedState conserved = gas.loadGas(cell);
        const std::array<double, 3> b = cellFieldComoving(geometry, faces, i, j, k);
        const double rho = conserved.rho_comoving;
        const double momentum2 =
            conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
            conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
            conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving;
        const double kinetic = rho > 0.0 ? 0.5 * momentum2 / rho : 0.0;
        const double magnetic = 0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor;
        budget.mass += rho * volume;
        budget.momentum[0] += conserved.momentum_density_x_comoving * volume;
        budget.momentum[1] += conserved.momentum_density_y_comoving * volume;
        budget.momentum[2] += conserved.momentum_density_z_comoving * volume;
        budget.energy_kinetic += kinetic * volume;
        budget.energy_magnetic += magnetic * volume;
        budget.energy_total += conserved.total_energy_density_comoving * volume;
      }
    }
  }
  budget.energy_thermal = budget.energy_total - budget.energy_kinetic - budget.energy_magnetic;
  const double normalizer =
      magneticDivergenceNormalizer(geometry, faces.phi_x, faces.phi_y, faces.phi_z);
  budget.divergence_normalizer = normalizer;
  budget.divergence_max_normalized =
      maxNormalizedFaceDivergence(geometry, faces.phi_x, faces.phi_y, faces.phi_z, normalizer);
  return budget;
}

}  // namespace cosmosim::mhd
