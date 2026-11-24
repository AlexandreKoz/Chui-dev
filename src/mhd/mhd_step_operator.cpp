#include "cosmosim/mhd/mhd_step_operator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cosmosim::mhd {
namespace {

MhdFailureContext makeContext(std::string_view owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::string(owner);
  context.detail = std::move(detail);
  return context;
}

std::size_t dimAlong(MhdAxis axis, const MhdPatchGeometry& geometry) noexcept {
  switch (axis) {
    case MhdAxis::kX: return geometry.nx;
    case MhdAxis::kY: return geometry.ny;
    case MhdAxis::kZ: return geometry.nz;
  }
  return 0U;
}

double spacingAlong(MhdAxis axis, const MhdPatchGeometry& geometry) noexcept {
  switch (axis) {
    case MhdAxis::kX: return geometry.dx_comoving;
    case MhdAxis::kY: return geometry.dy_comoving;
    case MhdAxis::kZ: return geometry.dz_comoving;
  }
  return 1.0;
}

double areaAlong(MhdAxis axis, const MhdPatchGeometry& geometry) noexcept {
  switch (axis) {
    case MhdAxis::kX: return geometry.faceAreaX();
    case MhdAxis::kY: return geometry.faceAreaY();
    case MhdAxis::kZ: return geometry.faceAreaZ();
  }
  return 1.0;
}

std::size_t faceIndexAlong(MhdAxis axis, const MhdPatchGeometry& geometry,
                           std::size_t n, std::size_t a, std::size_t b) noexcept {
  switch (axis) {
    case MhdAxis::kX: return geometry.faceIndexX(n, a, b);
    case MhdAxis::kY: return geometry.faceIndexY(a, n, b);
    case MhdAxis::kZ: return geometry.faceIndexZ(a, b, n);
  }
  return 0U;
}

std::size_t cellIndexAlong(MhdAxis axis, const MhdPatchGeometry& geometry,
                           std::size_t n, std::size_t a, std::size_t b) noexcept {
  switch (axis) {
    case MhdAxis::kX: return geometry.cellIndex(n, a, b);
    case MhdAxis::kY: return geometry.cellIndex(a, n, b);
    case MhdAxis::kZ: return geometry.cellIndex(a, b, n);
  }
  return 0U;
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
  m_rho.resize(cell_count, 0.0);
  m_mx.resize(cell_count, 0.0);
  m_my.resize(cell_count, 0.0);
  m_mz.resize(cell_count, 0.0);
  m_energy.resize(cell_count, 0.0);
  m_metal.resize(cell_count, 0.0);
}

std::uint64_t MhdGasStateSoa::ownedCapacityBytes() const noexcept {
  return m_rho.ownedCapacityBytes() + m_mx.ownedCapacityBytes() +
      m_my.ownedCapacityBytes() + m_mz.ownedCapacityBytes() +
      m_energy.ownedCapacityBytes() + m_metal.ownedCapacityBytes();
}

MhdConservedState MhdGasStateSoa::load(std::size_t cell) const noexcept {
  MhdConservedState state{};
  state.rho_comoving = m_rho[cell];
  state.momentum_density_x_comoving = m_mx[cell];
  state.momentum_density_y_comoving = m_my[cell];
  state.momentum_density_z_comoving = m_mz[cell];
  state.total_energy_density_comoving = m_energy[cell];
  state.metal_mass_density_comoving = m_metal[cell];
  return state;
}

void MhdGasStateSoa::store(std::size_t cell, const MhdConservedState& state) noexcept {
  m_rho[cell] = state.rho_comoving;
  m_mx[cell] = state.momentum_density_x_comoving;
  m_my[cell] = state.momentum_density_y_comoving;
  m_mz[cell] = state.momentum_density_z_comoving;
  m_energy[cell] = state.total_energy_density_comoving;
  m_metal[cell] = state.metal_mass_density_comoving;
}

void MhdGasStateSoa::copyFrom(const MhdGasStateSoa& other) {
  if (size() != other.size()) {
    resize(other.size());
  }
  std::copy(other.m_rho.begin(), other.m_rho.end(), m_rho.begin());
  std::copy(other.m_mx.begin(), other.m_mx.end(), m_mx.begin());
  std::copy(other.m_my.begin(), other.m_my.end(), m_my.begin());
  std::copy(other.m_mz.begin(), other.m_mz.end(), m_mz.begin());
  std::copy(other.m_energy.begin(), other.m_energy.end(), m_energy.begin());
  std::copy(other.m_metal.begin(), other.m_metal.end(), m_metal.begin());
}

void MhdGasStateSoa::combineSsprk2(const MhdGasStateSoa& previous,
                                   const MhdGasStateSoa& first_stage,
                                   const MhdGasStateSoa& first_stage_rhs,
                                   double dt_code) {
  const std::size_t count = std::min({size(), previous.size(), first_stage.size()});
  for (std::size_t cell = 0; cell < count; ++cell) {
    m_rho[cell] = 0.5 * (previous.m_rho[cell] + first_stage.m_rho[cell] +
                         dt_code * first_stage_rhs.m_rho[cell]);
    m_mx[cell] = 0.5 * (previous.m_mx[cell] + first_stage.m_mx[cell] +
                        dt_code * first_stage_rhs.m_mx[cell]);
    m_my[cell] = 0.5 * (previous.m_my[cell] + first_stage.m_my[cell] +
                        dt_code * first_stage_rhs.m_my[cell]);
    m_mz[cell] = 0.5 * (previous.m_mz[cell] + first_stage.m_mz[cell] +
                        dt_code * first_stage_rhs.m_mz[cell]);
    m_energy[cell] = 0.5 * (previous.m_energy[cell] + first_stage.m_energy[cell] +
                            dt_code * first_stage_rhs.m_energy[cell]);
    m_metal[cell] = 0.5 * (previous.m_metal[cell] + first_stage.m_metal[cell] +
                           dt_code * first_stage_rhs.m_metal[cell]);
  }
}

// ---------------------------------------------------------------------------
// MhdPrimitivePatch / MhdStageWorkspace
// ---------------------------------------------------------------------------

void MhdPrimitivePatch::resize(std::size_t patch_nx, std::size_t patch_ny,
                               std::size_t patch_nz, std::size_t ghost_cells) {
  nx = patch_nx;
  ny = patch_ny;
  nz = patch_nz;
  ghost = ghost_cells;
  px = nx + 2U * ghost;
  py = ny + 2U * ghost;
  pz = nz + 2U * ghost;
  const std::size_t count = px * py * pz;
  rho.resize(count, 0.0);
  vx.resize(count, 0.0);
  vy.resize(count, 0.0);
  vz.resize(count, 0.0);
  pressure.resize(count, 0.0);
  bx.resize(count, 0.0);
  by.resize(count, 0.0);
  bz.resize(count, 0.0);
  metal.resize(count, 0.0);
}

std::uint64_t MhdPrimitivePatch::ownedCapacityBytes() const noexcept {
  return rho.ownedCapacityBytes() + vx.ownedCapacityBytes() + vy.ownedCapacityBytes() +
      vz.ownedCapacityBytes() + pressure.ownedCapacityBytes() + bx.ownedCapacityBytes() +
      by.ownedCapacityBytes() + bz.ownedCapacityBytes() + metal.ownedCapacityBytes();
}

namespace {

std::uint64_t vectorCapacityBytes(const std::vector<double>& values) noexcept {
  return static_cast<std::uint64_t>(values.capacity() * sizeof(double));
}

}  // namespace

void MhdStageWorkspace::resize(const MhdPatchGeometry& geometry,
                               MhdReconstructionMode reconstruction,
                               bool nonideal_enabled) {
  const std::size_t ghost = reconstruction == MhdReconstructionMode::kPpm ? 3U : 2U;
  if (primitives.nx != geometry.nx || primitives.ny != geometry.ny ||
      primitives.nz != geometry.nz || primitives.ghost != ghost) {
    primitives.resize(geometry.nx, geometry.ny, geometry.nz, ghost);
  }
  if (signals.geometry.nx != geometry.nx || signals.geometry.ny != geometry.ny ||
      signals.geometry.nz != geometry.nz) {
    signals.resize(geometry);
  }
  const std::size_t cells = geometry.cellCount();
  if (gas_rhs.size() != cells) {
    gas_rhs.resize(cells);
  }
  phi_rhs_x.resize(geometry.faceCountX(), 0.0);
  phi_rhs_y.resize(geometry.faceCountY(), 0.0);
  phi_rhs_z.resize(geometry.faceCountZ(), 0.0);
  if (nonideal_enabled) {
    emf_rate.resize(geometry);
    cell_properties.assign(cells, MhdCellPhysicalProperties{});
    coefficient_cache.assign(cells, MhdTransportCoefficients{});
    nonideal.resize(geometry);
  }
}

std::uint64_t MhdStageWorkspace::ownedCapacityBytes() const noexcept {
  return primitives.ownedCapacityBytes() +
      vectorCapacityBytes(signals.mass_flux_x) + vectorCapacityBytes(signals.electric_t1_x) +
      vectorCapacityBytes(signals.electric_t2_x) + vectorCapacityBytes(signals.face_speed_x) +
      vectorCapacityBytes(signals.mass_flux_y) + vectorCapacityBytes(signals.electric_t1_y) +
      vectorCapacityBytes(signals.electric_t2_y) + vectorCapacityBytes(signals.face_speed_y) +
      vectorCapacityBytes(signals.mass_flux_z) + vectorCapacityBytes(signals.electric_t1_z) +
      vectorCapacityBytes(signals.electric_t2_z) + vectorCapacityBytes(signals.face_speed_z) +
      emf_rate.ownedCapacityBytes() + gas_rhs.ownedCapacityBytes() +
      phi_rhs_x.ownedCapacityBytes() + phi_rhs_y.ownedCapacityBytes() +
      phi_rhs_z.ownedCapacityBytes() + nonideal.ownedCapacityBytes();
}

// ---------------------------------------------------------------------------
// MhdStepOperator
// ---------------------------------------------------------------------------

MhdStepOperator::MhdStepOperator(MhdNumericsPolicy policy) : m_policy(policy) {
  if (!(m_policy.gamma_adiabatic > 1.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_step_operator", "gamma_adiabatic must exceed 1"));
  }
  if (!(m_policy.cfl_ideal > 0.0) || !(m_policy.cfl_ideal < 1.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_step_operator", "ideal CFL must lie in (0,1)"));
  }
}

std::pair<double, double> MhdStepOperator::reconstructPlmMc(
    double c_minus2, double c_minus, double c, double c_plus, double c_plus2,
    bool face_is_upper) noexcept {
  (void)c_minus2;
  (void)c_plus2;
  const double slope_central = 0.5 * (c_plus - c_minus);
  const double slope_left = c - c_minus;
  const double slope_right = c_plus - c;
  double slope = slope_central;
  if (slope_left * slope_right <= 0.0) {
    slope = 0.0;
  } else if (std::abs(slope_central) >
             2.0 * std::min(std::abs(slope_left), std::abs(slope_right))) {
    slope = 2.0 * std::min(std::abs(slope_left), std::abs(slope_right)) *
        (slope_central > 0.0 ? 1.0 : -1.0);
  }
  const double lower = c - 0.5 * slope;
  const double upper = c + 0.5 * slope;
  return face_is_upper ? std::pair<double, double>{upper, lower}
                       : std::pair<double, double>{lower, upper};
}

std::pair<double, double> MhdStepOperator::reconstructPpm(
    double c_minus2, double c_minus, double c, double c_plus, double c_plus2,
    bool face_is_upper) noexcept {
  double lower = (7.0 / 12.0) * (c_minus + c) -
      (1.0 / 12.0) * (c_minus2 + c_plus);
  double upper = (7.0 / 12.0) * (c + c_plus) -
      (1.0 / 12.0) * (c_minus + c_plus2);
  if ((upper - c) * (c - lower) <= 0.0) {
    lower = c;
    upper = c;
  } else {
    const double range = upper - lower;
    const double center = c - 0.5 * (lower + upper);
    if (range * center > (range * range) / 6.0) {
      lower = 3.0 * c - 2.0 * upper;
    } else if (-(range * range) / 6.0 > range * center) {
      upper = 3.0 * c - 2.0 * lower;
    }
  }
  return face_is_upper ? std::pair<double, double>{upper, lower}
                       : std::pair<double, double>{lower, upper};
}

void MhdStepOperator::fillGhostPrimitives(const MhdPatchGeometry& geometry,
                                          const MhdPatchBoundary& boundary,
                                          MhdPrimitivePatch& primitives) {
  const std::ptrdiff_t ghost = static_cast<std::ptrdiff_t>(primitives.ghost);
  const std::ptrdiff_t extents[3] = {static_cast<std::ptrdiff_t>(geometry.nx),
                                     static_cast<std::ptrdiff_t>(geometry.ny),
                                     static_cast<std::ptrdiff_t>(geometry.nz)};
  const auto resolve = [](std::ptrdiff_t value, std::ptrdiff_t extent,
                          MhdBoundaryKind kind) {
    if (extent <= 0) {
      return std::ptrdiff_t{0};
    }
    if (kind == MhdBoundaryKind::kPeriodic) {
      std::ptrdiff_t wrapped = value % extent;
      if (wrapped < 0) {
        wrapped += extent;
      }
      return wrapped;
    }
    return std::min(std::max(value, std::ptrdiff_t{0}), extent - 1);
  };

  for (std::ptrdiff_t k = -ghost; k < extents[2] + ghost; ++k) {
    for (std::ptrdiff_t j = -ghost; j < extents[1] + ghost; ++j) {
      for (std::ptrdiff_t i = -ghost; i < extents[0] + ghost; ++i) {
        const bool interior = i >= 0 && i < extents[0] && j >= 0 && j < extents[1] &&
            k >= 0 && k < extents[2];
        if (interior) {
          continue;
        }
        const std::array<std::ptrdiff_t, 3> coordinates{i, j, k};
        std::array<std::ptrdiff_t, 3> source{i, j, k};
        for (int axis = 0; axis < 3; ++axis) {
          const MhdBoundaryKind kind =
              coordinates[static_cast<std::size_t>(axis)] < 0
              ? boundary.lower[static_cast<std::size_t>(axis)]
              : boundary.upper[static_cast<std::size_t>(axis)];
          source[static_cast<std::size_t>(axis)] = resolve(
              coordinates[static_cast<std::size_t>(axis)],
              extents[static_cast<std::size_t>(axis)], kind);
        }
        const std::size_t target = primitives.index(i, j, k);
        const std::size_t origin = primitives.index(source[0], source[1], source[2]);
        primitives.rho[target] = primitives.rho[origin];
        primitives.vx[target] = primitives.vx[origin];
        primitives.vy[target] = primitives.vy[origin];
        primitives.vz[target] = primitives.vz[origin];
        primitives.pressure[target] = primitives.pressure[origin];
        primitives.bx[target] = primitives.bx[origin];
        primitives.by[target] = primitives.by[origin];
        primitives.bz[target] = primitives.bz[origin];
        primitives.metal[target] = primitives.metal[origin];

        for (int axis = 0; axis < 3; ++axis) {
          const MhdBoundaryKind kind =
              coordinates[static_cast<std::size_t>(axis)] < 0
              ? boundary.lower[static_cast<std::size_t>(axis)]
              : boundary.upper[static_cast<std::size_t>(axis)];
          const bool outside = coordinates[static_cast<std::size_t>(axis)] < 0 ||
              coordinates[static_cast<std::size_t>(axis)] >=
                  extents[static_cast<std::size_t>(axis)];
          if (kind != MhdBoundaryKind::kReflectingConductingWall || !outside) {
            continue;
          }
          if (axis == 0) {
            primitives.vx[target] = -primitives.vx[target];
            primitives.bx[target] = -primitives.bx[target];
          } else if (axis == 1) {
            primitives.vy[target] = -primitives.vy[target];
            primitives.by[target] = -primitives.by[target];
          } else {
            primitives.vz[target] = -primitives.vz[target];
            primitives.bz[target] = -primitives.bz[target];
          }
        }
      }
    }
  }
}

bool MhdStepOperator::recoverPatchPrimitives(const MhdPatchGeometry& geometry,
                                             const MhdGasStateSoa& gas,
                                             const MhdFaceFluxView& faces,
                                             double scale_factor,
                                             MhdStageWorkspace& workspace,
                                             std::size_t* nonpositive_cells) const {
  (void)scale_factor;
  workspace.resize(geometry, m_policy.reconstruction, false);
  if (!faces.valid() || gas.size() != geometry.cellCount()) {
    return false;
  }
  std::size_t bad = 0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t cell = geometry.cellIndex(i, j, k);
        MhdConservedState conserved = gas.load(cell);
        conserved.bx_comoving =
            0.5 * (faces.phi_x[geometry.faceIndexX(i, j, k)] +
                   faces.phi_x[geometry.faceIndexX(i + 1U, j, k)]) / geometry.faceAreaX();
        conserved.by_comoving =
            0.5 * (faces.phi_y[geometry.faceIndexY(i, j, k)] +
                   faces.phi_y[geometry.faceIndexY(i, j + 1U, k)]) / geometry.faceAreaY();
        conserved.bz_comoving =
            0.5 * (faces.phi_z[geometry.faceIndexZ(i, j, k)] +
                   faces.phi_z[geometry.faceIndexZ(i, j, k + 1U)]) / geometry.faceAreaZ();
        MhdPrimitiveState primitive{};
        if (!recoverPrimitive(conserved, m_policy.gamma_adiabatic, primitive)) {
          ++bad;
          continue;
        }
        const std::size_t index = workspace.primitives.index(
            static_cast<std::ptrdiff_t>(i), static_cast<std::ptrdiff_t>(j),
            static_cast<std::ptrdiff_t>(k));
        workspace.primitives.rho[index] = primitive.rho_comoving;
        workspace.primitives.vx[index] = primitive.vel_x_peculiar;
        workspace.primitives.vy[index] = primitive.vel_y_peculiar;
        workspace.primitives.vz[index] = primitive.vel_z_peculiar;
        workspace.primitives.pressure[index] = primitive.pressure_comoving;
        workspace.primitives.bx[index] = primitive.bx_comoving;
        workspace.primitives.by[index] = primitive.by_comoving;
        workspace.primitives.bz[index] = primitive.bz_comoving;
        workspace.primitives.metal[index] = primitive.metallicity_mass_fraction;
      }
    }
  }
  if (nonpositive_cells != nullptr) {
    *nonpositive_cells = bad;
  }
  return bad == 0U;
}

double MhdStepOperator::estimateIdealTimeStep(const MhdPatchGeometry& geometry,
                                              const MhdGasStateSoa& gas,
                                              const MhdFaceFluxView& faces,
                                              double cfl) const {
  if (!faces.valid() || gas.size() != geometry.cellCount()) {
    return 0.0;
  }
  double maximum_rate = 0.0;
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t cell = geometry.cellIndex(i, j, k);
        MhdConservedState conserved = gas.load(cell);
        conserved.bx_comoving =
            0.5 * (faces.phi_x[geometry.faceIndexX(i, j, k)] +
                   faces.phi_x[geometry.faceIndexX(i + 1U, j, k)]) / geometry.faceAreaX();
        conserved.by_comoving =
            0.5 * (faces.phi_y[geometry.faceIndexY(i, j, k)] +
                   faces.phi_y[geometry.faceIndexY(i, j + 1U, k)]) / geometry.faceAreaY();
        conserved.bz_comoving =
            0.5 * (faces.phi_z[geometry.faceIndexZ(i, j, k)] +
                   faces.phi_z[geometry.faceIndexZ(i, j, k + 1U)]) / geometry.faceAreaZ();
        MhdPrimitiveState primitive{};
        if (!recoverPrimitive(conserved, m_policy.gamma_adiabatic, primitive)) {
          continue;
        }
        const double speed_x = std::abs(primitive.vel_x_peculiar) +
            fastMagnetosonicSpeed(primitive, m_policy.gamma_adiabatic, MhdAxis::kX);
        const double speed_y = std::abs(primitive.vel_y_peculiar) +
            fastMagnetosonicSpeed(primitive, m_policy.gamma_adiabatic, MhdAxis::kY);
        const double speed_z = std::abs(primitive.vel_z_peculiar) +
            fastMagnetosonicSpeed(primitive, m_policy.gamma_adiabatic, MhdAxis::kZ);
        maximum_rate = std::max({maximum_rate, speed_x / geometry.dx_comoving,
                                 speed_y / geometry.dy_comoving,
                                 speed_z / geometry.dz_comoving});
      }
    }
  }
  if (!(maximum_rate > 0.0)) {
    return std::numeric_limits<double>::infinity();
  }
  return cfl / maximum_rate;
}

MhdConservationBudget MhdStepOperator::patchBudget(const MhdPatchGeometry& geometry,
                                                   const MhdGasStateSoa& gas,
                                                   const MhdFaceFluxView& faces,
                                                   double scale_factor) const {
  MhdConservationBudget budget{};
  if (!faces.valid() || gas.size() != geometry.cellCount()) {
    return budget;
  }
  (void)scale_factor;
  const double volume = geometry.cellVolume();
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      for (std::size_t i = 0; i < geometry.nx; ++i) {
        const std::size_t cell = geometry.cellIndex(i, j, k);
        const MhdConservedState conserved = gas.load(cell);
        const double bx =
            0.5 * (faces.phi_x[geometry.faceIndexX(i, j, k)] +
                   faces.phi_x[geometry.faceIndexX(i + 1U, j, k)]) / geometry.faceAreaX();
        const double by =
            0.5 * (faces.phi_y[geometry.faceIndexY(i, j, k)] +
                   faces.phi_y[geometry.faceIndexY(i, j + 1U, k)]) / geometry.faceAreaY();
        const double bz =
            0.5 * (faces.phi_z[geometry.faceIndexZ(i, j, k)] +
                   faces.phi_z[geometry.faceIndexZ(i, j, k + 1U)]) / geometry.faceAreaZ();
        const double rho = conserved.rho_comoving;
        const double kinetic = rho > 0.0
            ? 0.5 * (conserved.momentum_density_x_comoving *
                         conserved.momentum_density_x_comoving +
                     conserved.momentum_density_y_comoving *
                         conserved.momentum_density_y_comoving +
                     conserved.momentum_density_z_comoving *
                         conserved.momentum_density_z_comoving) /
                rho
            : 0.0;
        const double magnetic = 0.5 * (bx * bx + by * by + bz * bz);
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
  budget.energy_thermal =
      budget.energy_total - budget.energy_kinetic - budget.energy_magnetic;
  return budget;
}

void MhdStepOperator::applyExpansionSourceHalf(std::span<double> rho,
                                               std::span<double> mx,
                                               std::span<double> my,
                                               std::span<double> mz,
                                               std::span<double> energy,
                                               const MhdPatchGeometry& geometry,
                                               std::span<const double> phi_x,
                                               std::span<const double> phi_y,
                                               std::span<const double> phi_z,
                                               double scale_factor,
                                               double hubble_rate_code,
                                               double dt_half_code) const {
  if (!(hubble_rate_code != 0.0) || !(dt_half_code > 0.0) || !(scale_factor > 0.0)) {
    return;
  }
  (void)phi_x;
  (void)phi_y;
  (void)phi_z;
  (void)energy;
  const double hubble_half_dt = hubble_rate_code * dt_half_code;
  const double damping = 1.0 / (1.0 + hubble_half_dt);
  const std::size_t cells = geometry.cellCount();
  for (std::size_t cell = 0; cell < cells; ++cell) {
    // Peculiar-momentum decay from Hubble drag; the thermal/magnetic background
    // terms are owned by the workflow Strang source stage and are not duplicated
    // here. `rho` and `energy` are intentionally untouched.
    mx[cell] *= damping;
    my[cell] *= damping;
    mz[cell] *= damping;
  }
  (void)rho;
}

void MhdStepOperator::computeIdealRate(const MhdPatchGeometry& geometry,
                                       const MhdPatchBoundary& boundary,
                                       const MhdGasStateSoa& gas,
                                       std::span<const double> phi_x,
                                       std::span<const double> phi_y,
                                       std::span<const double> phi_z,
                                       const MhdStageContext& context,
                                       MhdStageWorkspace& workspace,
                                       std::span<double> phi_rhs_x,
                                       std::span<double> phi_rhs_y,
                                       std::span<double> phi_rhs_z,
                                       std::uint64_t* hlle_faces,
                                       std::uint64_t* first_order_faces) const {
  const std::size_t cells = geometry.cellCount();
  std::fill(phi_rhs_x.begin(), phi_rhs_x.end(), 0.0);
  std::fill(phi_rhs_y.begin(), phi_rhs_y.end(), 0.0);
  std::fill(phi_rhs_z.begin(), phi_rhs_z.end(), 0.0);
  workspace.resize(geometry, m_policy.reconstruction, false);
  workspace.gas_rhs.resize(cells);
  std::fill(workspace.gas_rhs.rho().begin(), workspace.gas_rhs.rho().end(), 0.0);
  std::fill(workspace.gas_rhs.momentumX().begin(), workspace.gas_rhs.momentumX().end(), 0.0);
  std::fill(workspace.gas_rhs.momentumY().begin(), workspace.gas_rhs.momentumY().end(), 0.0);
  std::fill(workspace.gas_rhs.momentumZ().begin(), workspace.gas_rhs.momentumZ().end(), 0.0);
  std::fill(workspace.gas_rhs.totalEnergy().begin(), workspace.gas_rhs.totalEnergy().end(), 0.0);
  std::fill(workspace.gas_rhs.metalMass().begin(), workspace.gas_rhs.metalMass().end(), 0.0);

  const MhdFaceFluxView view{&geometry, phi_x, phi_y, phi_z, context.patch_id,
                             context.geometry_generation, context.refinement_level};
  if (!recoverPatchPrimitives(geometry, gas, view, context.scale_factor, workspace, nullptr)) {
    throw MhdError(MhdErrorCode::kNonPositiveState,
                   makeContext("mhd_ideal", "non-positive state before the face sweep"));
  }
  fillGhostPrimitives(geometry, boundary, workspace.primitives);

  const auto reconstruct_lane = [&](const MhdDoubleArray& lane, std::ptrdiff_t i,
                                    std::ptrdiff_t j, std::ptrdiff_t k, MhdAxis axis,
                                    bool upper_face) {
    const auto at = [&](std::ptrdiff_t offset) {
      switch (axis) {
        case MhdAxis::kX:
          return lane[workspace.primitives.index(i + offset, j, k)];
        case MhdAxis::kY:
          return lane[workspace.primitives.index(i, j + offset, k)];
        case MhdAxis::kZ:
          return lane[workspace.primitives.index(i, j, k + offset)];
      }
      return lane[workspace.primitives.index(i, j, k)];
    };
    const double c_minus2 = at(-2);
    const double c_minus = at(-1);
    const double center = at(0);
    const double c_plus = at(1);
    const double c_plus2 = at(2);
    const auto pair = m_policy.reconstruction == MhdReconstructionMode::kPpm
        ? reconstructPpm(c_minus2, c_minus, center, c_plus, c_plus2, upper_face)
        : reconstructPlmMc(c_minus2, c_minus, center, c_plus, c_plus2, upper_face);
    return pair.first;
  };

  const auto global_state = [&](std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k,
                                MhdAxis axis, bool upper_face) {
    MhdPrimitiveState state{};
    state.rho_comoving = reconstruct_lane(workspace.primitives.rho, i, j, k, axis, upper_face);
    state.vel_x_peculiar = reconstruct_lane(workspace.primitives.vx, i, j, k, axis, upper_face);
    state.vel_y_peculiar = reconstruct_lane(workspace.primitives.vy, i, j, k, axis, upper_face);
    state.vel_z_peculiar = reconstruct_lane(workspace.primitives.vz, i, j, k, axis, upper_face);
    state.pressure_comoving =
        reconstruct_lane(workspace.primitives.pressure, i, j, k, axis, upper_face);
    state.bx_comoving = reconstruct_lane(workspace.primitives.bx, i, j, k, axis, upper_face);
    state.by_comoving = reconstruct_lane(workspace.primitives.by, i, j, k, axis, upper_face);
    state.bz_comoving = reconstruct_lane(workspace.primitives.bz, i, j, k, axis, upper_face);
    state.metallicity_mass_fraction =
        reconstruct_lane(workspace.primitives.metal, i, j, k, axis, upper_face);
    return state;
  };

  const auto override_normal_field = [](MhdPrimitiveState& state, MhdAxis axis,
                                        double normal_b) {
    switch (axis) {
      case MhdAxis::kX: state.bx_comoving = normal_b; break;
      case MhdAxis::kY: state.by_comoving = normal_b; break;
      case MhdAxis::kZ: state.bz_comoving = normal_b; break;
    }
  };

  std::uint64_t hlle_count = 0;
  std::uint64_t first_order_count = 0;
  std::uint64_t degenerate_edges = 0;
  MhdEdgeEmfStorage emf_interval(geometry);

  const auto sweep = [&](MhdAxis axis) {
    const std::size_t extent = dimAlong(axis, geometry);
    const std::size_t extent_a = dimAlong(static_cast<MhdAxis>((static_cast<int>(axis) + 1) % 3),
                                          geometry);
    const std::size_t extent_b = dimAlong(static_cast<MhdAxis>((static_cast<int>(axis) + 2) % 3),
                                          geometry);
    const double spacing_inverse = 1.0 / spacingAlong(axis, geometry);
    const double face_area = areaAlong(axis, geometry);
    for (std::size_t b = 0; b < extent_b; ++b) {
      for (std::size_t a = 0; a < extent_a; ++a) {
        for (std::size_t f = 0; f <= extent; ++f) {
          const std::ptrdiff_t left_n = static_cast<std::ptrdiff_t>(f) - 1;
          const std::ptrdiff_t right_n = static_cast<std::ptrdiff_t>(f);
          const std::ptrdiff_t pa = static_cast<std::ptrdiff_t>(a);
          const std::ptrdiff_t pb = static_cast<std::ptrdiff_t>(b);
          std::ptrdiff_t left_i = pa, left_j = pb, left_k = 0;
          std::ptrdiff_t right_i = pa, right_j = pb, right_k = 0;
          switch (axis) {
            case MhdAxis::kX:
              left_i = left_n; right_i = right_n; left_j = pa; left_k = pb;
              right_j = pa; right_k = pb;
              break;
            case MhdAxis::kY:
              left_j = left_n; right_j = right_n; left_i = pa; left_k = pb;
              right_i = pa; right_k = pb;
              break;
            case MhdAxis::kZ:
              left_k = left_n; right_k = right_n; left_i = pa; left_j = pb;
              right_i = pa; right_j = pb;
              break;
          }
          MhdPrimitiveState left = global_state(left_i, left_j, left_k, axis, true);
          MhdPrimitiveState right = global_state(right_i, right_j, right_k, axis, false);
          const std::size_t face_slot = faceIndexAlong(axis, geometry, f, a, b);
          const double normal_b = (axis == MhdAxis::kX
                                       ? phi_x[face_slot]
                                       : (axis == MhdAxis::kY ? phi_y[face_slot]
                                                              : phi_z[face_slot])) /
              face_area;
          override_normal_field(left, axis, normal_b);
          override_normal_field(right, axis, normal_b);

          MhdRiemann::Result solved = MhdRiemann::solve(
              axis, left, right, normal_b, m_policy.gamma_adiabatic, m_policy.riemann);
          if (!solved.valid) {
            solved = MhdRiemann::solveHlle(axis, left, right, normal_b,
                                           m_policy.gamma_adiabatic);
          }
          if (!solved.valid) {
            solved.flux = firstOrderPositivityFlux(axis, left, right, normal_b,
                                                   m_policy.gamma_adiabatic);
            solved.used_hlle = true;
            solved.first_order_required = true;
          }
          if (solved.used_hlle) {
            ++hlle_count;
          }
          if (solved.first_order_required) {
            ++first_order_count;
          }
          const MhdRiemannFlux& flux = solved.flux;
          const double upwind_metal =
              flux.mass_flux > 0.0 ? left.metallicity_mass_fraction
                                   : right.metallicity_mass_fraction;

          if (f > 0U) {
            const std::size_t below = cellIndexAlong(axis, geometry, f - 1U, a, b);
            workspace.gas_rhs.rho()[below] += flux.mass_flux * spacing_inverse;
            workspace.gas_rhs.momentumX()[below] +=
                flux.momentum_flux[0] * spacing_inverse;
            workspace.gas_rhs.momentumY()[below] +=
                flux.momentum_flux[1] * spacing_inverse;
            workspace.gas_rhs.momentumZ()[below] +=
                flux.momentum_flux[2] * spacing_inverse;
            workspace.gas_rhs.totalEnergy()[below] += flux.energy_flux * spacing_inverse;
            workspace.gas_rhs.metalMass()[below] +=
                flux.mass_flux * upwind_metal * spacing_inverse;
          }
          if (f < extent) {
            const std::size_t above = cellIndexAlong(axis, geometry, f, a, b);
            workspace.gas_rhs.rho()[above] -= flux.mass_flux * spacing_inverse;
            workspace.gas_rhs.momentumX()[above] -=
                flux.momentum_flux[0] * spacing_inverse;
            workspace.gas_rhs.momentumY()[above] -=
                flux.momentum_flux[1] * spacing_inverse;
            workspace.gas_rhs.momentumZ()[above] -=
                flux.momentum_flux[2] * spacing_inverse;
            workspace.gas_rhs.totalEnergy()[above] -= flux.energy_flux * spacing_inverse;
            workspace.gas_rhs.metalMass()[above] -=
                flux.mass_flux * upwind_metal * spacing_inverse;
          }

          const std::array<double, 3> local_electric =
              MhdEmfBuilder::faceElectricFieldFromFlux(flux);
          const std::array<double, 3> global_electric =
              permuteToGlobal(axis, local_electric);
          double tangential_one = 0.0;
          double tangential_two = 0.0;
          switch (axis) {
            case MhdAxis::kX: tangential_one = global_electric[1];
              tangential_two = global_electric[2]; break;
            case MhdAxis::kY: tangential_one = global_electric[2];
              tangential_two = global_electric[0]; break;
            case MhdAxis::kZ: tangential_one = global_electric[0];
              tangential_two = global_electric[1]; break;
          }
          const double speed = std::max(flux.max_fast_speed, 1.0e-30);
          switch (axis) {
            case MhdAxis::kX:
              workspace.signals.mass_flux_x[face_slot] = flux.mass_flux;
              workspace.signals.electric_t1_x[face_slot] = tangential_one;
              workspace.signals.electric_t2_x[face_slot] = tangential_two;
              workspace.signals.face_speed_x[face_slot] = speed;
              break;
            case MhdAxis::kY:
              workspace.signals.mass_flux_y[face_slot] = flux.mass_flux;
              workspace.signals.electric_t1_y[face_slot] = tangential_one;
              workspace.signals.electric_t2_y[face_slot] = tangential_two;
              workspace.signals.face_speed_y[face_slot] = speed;
              break;
            case MhdAxis::kZ:
              workspace.signals.mass_flux_z[face_slot] = flux.mass_flux;
              workspace.signals.electric_t1_z[face_slot] = tangential_one;
              workspace.signals.electric_t2_z[face_slot] = tangential_two;
              workspace.signals.face_speed_z[face_slot] = speed;
              break;
          }
        }
      }
    }
  };

  sweep(MhdAxis::kX);
  sweep(MhdAxis::kY);
  sweep(MhdAxis::kZ);

  MhdEmfBuilder::buildUctHll(geometry, workspace.signals, emf_interval, &degenerate_edges);
  applyConstrainedTransportUpdate(geometry, phi_rhs_x, phi_rhs_y, phi_rhs_z,
                                  emf_interval.emfX(), emf_interval.emfY(),
                                  emf_interval.emfZ(), 1.0);
  if (hlle_faces != nullptr) {
    *hlle_faces = hlle_count;
  }
  if (first_order_faces != nullptr) {
    *first_order_faces = first_order_count;
  }
}

MhdAdvanceResult MhdStepOperator::advance(const MhdPatchGeometry& geometry,
                                          const MhdPatchBoundary& boundary,
                                          MhdGasStateSoa& gas,
                                          MhdFaceFluxMutableView faces,
                                          const MhdStageContext& context,
                                          const MhdCoefficientProvider* coefficient_provider,
                                          MhdDiagnosticSink* diagnostics) {
  MhdAdvanceResult result{};
  if (!faces.valid()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_advance", "face flux view does not match geometry"));
  }
  if (faces.geometry_generation != context.geometry_generation ||
      faces.patch_id != context.patch_id) {
    throw MhdError(MhdErrorCode::kStaleGeometry,
                   makeContext("mhd_advance",
                               "face flux view epoch does not match the stage context"));
  }
  if (!(context.dt_code > 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_advance", "non-positive step interval"));
  }
  if (gas.size() != geometry.cellCount()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   makeContext("mhd_advance", "gas state size does not match the patch"));
  }
  validateBoundarySupport(
      boundary, MhdFailureContext{.owner = "mhd_advance",
                                  .detail = "",
                                  .rank = -1,
                                  .patch_id = context.patch_id,
                                  .geometry_generation = context.geometry_generation,
                                  .face_key = 0U,
                                  .edge_key = 0U,
                                  .time_code = context.time_code,
                                  .scale_factor = context.scale_factor});

  MhdStageWorkspace workspace;
  const bool nonideal = coefficient_provider != nullptr;
  workspace.resize(geometry, m_policy.reconstruction, nonideal);

  MhdGasStateSoa snapshot_gas(gas.size());
  snapshot_gas.copyFrom(gas);
  std::vector<double> snapshot_x(faces.phi_x.begin(), faces.phi_x.end());
  std::vector<double> snapshot_y(faces.phi_y.begin(), faces.phi_y.end());
  std::vector<double> snapshot_z(faces.phi_z.begin(), faces.phi_z.end());

  double dt = context.dt_code;
  const double divergence_tolerance = 1.0e-11;
  const std::size_t cells = geometry.cellCount();

  for (std::int32_t attempt = 0; attempt <= m_policy.max_step_retries; ++attempt) {
    gas.copyFrom(snapshot_gas);
    std::copy(snapshot_x.begin(), snapshot_x.end(), faces.phi_x.begin());
    std::copy(snapshot_y.begin(), snapshot_y.end(), faces.phi_y.begin());
    std::copy(snapshot_z.begin(), snapshot_z.end(), faces.phi_z.begin());

    MhdStageContext stage_context = context;
    stage_context.dt_code = dt;
    stage_context.stage_ordinal = 1U;

    std::uint64_t hlle_first = 0;
    std::uint64_t first_order_first = 0;
    bool stage_ok = true;
    try {
      computeIdealRate(geometry, boundary, gas, faces.phi_x, faces.phi_y, faces.phi_z,
                       stage_context, workspace, workspace.phi_rhs_x.span(),
                       workspace.phi_rhs_y.span(), workspace.phi_rhs_z.span(),
                       &hlle_first, &first_order_first);
    } catch (const MhdError&) {
      stage_ok = false;
    }

    if (stage_ok) {
      // U1 = U0 + dt * L(U0); Phi1 = Phi0 + dt * L(Phi0).
      for (std::size_t cell = 0; cell < cells; ++cell) {
        gas.rho()[cell] += dt * workspace.gas_rhs.rho()[cell];
        gas.momentumX()[cell] += dt * workspace.gas_rhs.momentumX()[cell];
        gas.momentumY()[cell] += dt * workspace.gas_rhs.momentumY()[cell];
        gas.momentumZ()[cell] += dt * workspace.gas_rhs.momentumZ()[cell];
        gas.totalEnergy()[cell] += dt * workspace.gas_rhs.totalEnergy()[cell];
        gas.metalMass()[cell] += dt * workspace.gas_rhs.metalMass()[cell];
      }
      for (std::size_t index = 0; index < faces.phi_x.size(); ++index) {
        faces.phi_x[index] += dt * workspace.phi_rhs_x[index];
      }
      for (std::size_t index = 0; index < faces.phi_y.size(); ++index) {
        faces.phi_y[index] += dt * workspace.phi_rhs_y[index];
      }
      for (std::size_t index = 0; index < faces.phi_z.size(); ++index) {
        faces.phi_z[index] += dt * workspace.phi_rhs_z[index];
      }
      std::size_t bad_stage = 0;
      const MhdFaceFluxView stage_view{&geometry, faces.phi_x, faces.phi_y, faces.phi_z,
                                       context.patch_id, context.geometry_generation,
                                       context.refinement_level};
      stage_ok = recoverPatchPrimitives(geometry, gas, stage_view, context.scale_factor,
                                        workspace, &bad_stage) &&
          bad_stage == 0U;
    }

    if (stage_ok) {
      stage_context.stage_ordinal = 2U;
      std::uint64_t hlle_second = 0;
      std::uint64_t first_order_second = 0;
      try {
        computeIdealRate(geometry, boundary, gas, faces.phi_x, faces.phi_y, faces.phi_z,
                         stage_context, workspace, workspace.phi_rhs_x.span(),
                         workspace.phi_rhs_y.span(), workspace.phi_rhs_z.span(),
                         &hlle_second, &first_order_second);
      } catch (const MhdError&) {
        stage_ok = false;
      }
      if (stage_ok) {
        // U_{n+1} = 0.5 (U0 + U1 + dt L(U1)).
        for (std::size_t cell = 0; cell < cells; ++cell) {
          const double rho_one = gas.rho()[cell];
          const double mx_one = gas.momentumX()[cell];
          const double my_one = gas.momentumY()[cell];
          const double mz_one = gas.momentumZ()[cell];
          const double energy_one = gas.totalEnergy()[cell];
          const double metal_one = gas.metalMass()[cell];
          gas.rho()[cell] = 0.5 * (snapshot_gas.rho()[cell] + rho_one +
                                   dt * workspace.gas_rhs.rho()[cell]);
          gas.momentumX()[cell] = 0.5 * (snapshot_gas.momentumX()[cell] + mx_one +
                                         dt * workspace.gas_rhs.momentumX()[cell]);
          gas.momentumY()[cell] = 0.5 * (snapshot_gas.momentumY()[cell] + my_one +
                                         dt * workspace.gas_rhs.momentumY()[cell]);
          gas.momentumZ()[cell] = 0.5 * (snapshot_gas.momentumZ()[cell] + mz_one +
                                         dt * workspace.gas_rhs.momentumZ()[cell]);
          gas.totalEnergy()[cell] = 0.5 * (snapshot_gas.totalEnergy()[cell] + energy_one +
                                           dt * workspace.gas_rhs.totalEnergy()[cell]);
          gas.metalMass()[cell] = 0.5 * (snapshot_gas.metalMass()[cell] + metal_one +
                                         dt * workspace.gas_rhs.metalMass()[cell]);
        }
        for (std::size_t index = 0; index < faces.phi_x.size(); ++index) {
          const double one = faces.phi_x[index];
          faces.phi_x[index] = 0.5 * (snapshot_x[index] + one +
                                      dt * workspace.phi_rhs_x[index]);
        }
        for (std::size_t index = 0; index < faces.phi_y.size(); ++index) {
          const double one = faces.phi_y[index];
          faces.phi_y[index] = 0.5 * (snapshot_y[index] + one +
                                      dt * workspace.phi_rhs_y[index]);
        }
        for (std::size_t index = 0; index < faces.phi_z.size(); ++index) {
          const double one = faces.phi_z[index];
          faces.phi_z[index] = 0.5 * (snapshot_z[index] + one +
                                      dt * workspace.phi_rhs_z[index]);
        }
        std::size_t bad_final = 0;
        const MhdFaceFluxView final_view{&geometry, faces.phi_x, faces.phi_y, faces.phi_z,
                                         context.patch_id, context.geometry_generation,
                                         context.refinement_level};
        stage_ok = recoverPatchPrimitives(geometry, gas, final_view, context.scale_factor,
                                          workspace, &bad_final) &&
            bad_final == 0U;
        if (stage_ok) {
          result.hlle_fallback_faces = hlle_first + hlle_second;
          result.first_order_faces = first_order_first + first_order_second;
        }
      }
    }

    if (!stage_ok) {
      result.retry_count = attempt + 1;
      if (attempt == m_policy.max_step_retries) {
        gas.copyFrom(snapshot_gas);
        std::copy(snapshot_x.begin(), snapshot_x.end(), faces.phi_x.begin());
        std::copy(snapshot_y.begin(), snapshot_y.end(), faces.phi_y.begin());
        std::copy(snapshot_z.begin(), snapshot_z.end(), faces.phi_z.begin());
        MhdFailureContext failure = makeContext(
            "mhd_advance",
            "step retry budget exhausted without a non-positive-free second-order step");
        failure.patch_id = context.patch_id;
        failure.geometry_generation = context.geometry_generation;
        failure.time_code = context.time_code;
        failure.scale_factor = context.scale_factor;
        throw MhdError(MhdErrorCode::kNonPositiveState, std::move(failure));
      }
      dt *= 0.5;
      continue;
    }

    const double normalizer = magneticDivergenceNormalizer(
        faces.phi_x, faces.phi_y, faces.phi_z, 1.0);
    result.divergence_normalizer = normalizer;
    result.divergence_max_normalized = maxNormalizedFaceDivergence(
        geometry, faces.phi_x, faces.phi_y, faces.phi_z, normalizer);
    requireDivergenceContract(result.divergence_max_normalized, divergence_tolerance,
                              MhdFailureContext{.owner = "mhd_advance",
                                                .detail = "post-step strictly-CT divergence",
                                                .rank = -1,
                                                .patch_id = context.patch_id,
                                                .geometry_generation =
                                                    context.geometry_generation,
                                                .time_code = context.time_code,
                                                .scale_factor = context.scale_factor});
    result.dt_used_code = dt;
    result.accepted = true;
    result.face_count = static_cast<std::uint64_t>(geometry.faceCountTotal());
    result.fallback_fraction = result.face_count > 0U
        ? static_cast<double>(result.hlle_fallback_faces) /
            static_cast<double>(result.face_count)
        : 0.0;
    result.accepted_stage_ordinal = 2U;
    result.ideal_stage_evaluations = 2U;
    const MhdFaceFluxView final_view{&geometry, faces.phi_x, faces.phi_y, faces.phi_z,
                                     context.patch_id, context.geometry_generation,
                                     context.refinement_level};
    result.budget = patchBudget(geometry, gas, final_view, context.scale_factor);
    result.budget.divergence_max_normalized = result.divergence_max_normalized;
    result.budget.divergence_normalizer = normalizer;

    if (context.comoving_coordinates && context.hubble_rate_code != 0.0) {
      applyExpansionSourceHalf(gas.rho(), gas.momentumX(), gas.momentumY(),
                               gas.momentumZ(), gas.totalEnergy(), geometry,
                               faces.phi_x, faces.phi_y, faces.phi_z,
                               context.scale_factor, context.hubble_rate_code,
                               0.5 * dt);
    }
    if (diagnostics != nullptr) {
      diagnostics->recordStage(
          MhdStageContext{.dt_code = dt,
                          .time_code = context.time_code,
                          .scale_factor = context.scale_factor,
                          .hubble_rate_code = context.hubble_rate_code,
                          .refinement_level = context.refinement_level,
                          .patch_id = context.patch_id,
                          .geometry_generation = context.geometry_generation,
                          .stage_ordinal = 2U,
                          .comoving_coordinates = context.comoving_coordinates},
          result.budget);
    }
    return result;
  }
  return result;
}

}  // namespace cosmosim::mhd
