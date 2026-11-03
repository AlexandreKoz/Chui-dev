#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/memory_governor.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cosmosim::hydro {
namespace {

constexpr double k_small = 1.0e-14;
constexpr std::size_t k_max_faces_per_active_cell = 6U;

[[nodiscard]] double dot3(double ax, double ay, double az, double bx, double by, double bz) {
  return ax * bx + ay * by + az * bz;
}

[[nodiscard]] double norm3(double x, double y, double z) {
  return std::sqrt(dot3(x, y, z, x, y, z));
}

[[nodiscard]] std::array<double, 3> gravityAtCell(
    std::size_t cell_index,
    std::span<const double> ax,
    std::span<const double> ay,
    std::span<const double> az) {
  const double gx = (cell_index < ax.size()) ? ax[cell_index] : 0.0;
  const double gy = (cell_index < ay.size()) ? ay[cell_index] : 0.0;
  const double gz = (cell_index < az.size()) ? az[cell_index] : 0.0;
  return {gx, gy, gz};
}

void validateAdvanceInputs(
    const HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroUpdateContext& update,
    const HydroSourceContext& source_context,
    double adiabatic_index) {
  if (conserved.size() == 0) {
    throw std::invalid_argument("Hydro advance requires at least one cell");
  }
  if (geometry.cell_volume_comoving <= 0.0) {
    throw std::invalid_argument("Hydro advance requires positive cell_volume_comoving");
  }
  if (update.dt_code <= 0.0) {
    throw std::invalid_argument("Hydro advance requires dt_code > 0");
  }
  if (update.scale_factor <= 0.0) {
    throw std::invalid_argument("Hydro advance requires scale_factor > 0");
  }
  if (adiabatic_index <= 1.0) {
    throw std::invalid_argument("Hydro solver requires adiabatic_index > 1");
  }
  if (std::abs(source_context.update.dt_code - update.dt_code) > 1.0e-16 ||
      std::abs(source_context.update.scale_factor - update.scale_factor) > 1.0e-16 ||
      std::abs(source_context.update.hubble_rate_code - update.hubble_rate_code) > 1.0e-16) {
    throw std::invalid_argument("Hydro source_context.update must match update passed to advancePatch");
  }

  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  for (const HydroFace& face : geometry.faces) {
    if (face.owner_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face owner index is out of range");
    }
    if (face.owner_cell >= real_cell_count) {
      throw std::invalid_argument("Hydro face owner must be an authoritative real cell");
    }
    if (face.neighbor_cell != k_invalid_cell_index && face.neighbor_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face neighbor index is out of range");
    }
    if (face.owner_minus_cell != k_invalid_cell_index && face.owner_minus_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face owner-minus stencil index is out of range");
    }
    if (face.neighbor_plus_cell != k_invalid_cell_index && face.neighbor_plus_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face neighbor-plus stencil index is out of range");
    }
    if (face.ghost_cell_slot != k_invalid_ghost_cell_slot) {
      if (face.ghost_cell_slot >= geometry.ghost_cells.size()) {
        throw std::invalid_argument("Hydro face references an invalid ghost cell slot");
      }
      const HydroGhostCell& ghost = geometry.ghost_cells[face.ghost_cell_slot];
      if (face.neighbor_cell != ghost.ghost_cell || face.owner_cell != ghost.owner_real_cell) {
        throw std::invalid_argument("Hydro face ghost metadata does not match face endpoints");
      }
    }
    if (face.area_comoving <= 0.0) {
      throw std::invalid_argument("Hydro face area_comoving must be positive");
    }
    const double normal_norm = norm3(face.normal_x, face.normal_y, face.normal_z);
    if (std::abs(normal_norm - 1.0) > 1.0e-10) {
      throw std::invalid_argument("Hydro face normal must be unit length");
    }
  }
  if (!geometry.flux_register_faces.empty() &&
      geometry.flux_register_faces.size() != geometry.faces.size()) {
    throw std::invalid_argument("Hydro flux-register metadata must be empty or match face count");
  }
  for (const HydroFluxRegisterFace& register_face : geometry.flux_register_faces) {
    if (register_face.role == HydroFluxRegisterFaceRole::kNone) {
      continue;
    }
    if (register_face.register_key == 0U) {
      throw std::invalid_argument("Hydro flux-register face requires nonzero register_key");
    }
    if (register_face.coarse_patch_id == 0U) {
      throw std::invalid_argument("Hydro flux-register face requires nonzero coarse_patch_id");
    }
    if (register_face.coarse_cell_index == k_invalid_cell_index) {
      throw std::invalid_argument("Hydro flux-register face requires valid coarse_cell_index");
    }
    if (register_face.coarse_orientation_sign != 1.0 &&
        register_face.coarse_orientation_sign != -1.0) {
      throw std::invalid_argument("Hydro flux-register coarse_orientation_sign must be +1 or -1");
    }
  }
}

void validateActiveSet(const HydroActiveSetView& active_set, const HydroConservedStateSoa& conserved, const HydroPatchGeometry& geometry) {
  if (active_set.active_cells.empty()) {
    throw std::invalid_argument("Hydro active set requires at least one active cell");
  }

  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  for (std::size_t cell_index : active_set.active_cells) {
    if (cell_index >= conserved.size()) {
      throw std::invalid_argument("Hydro active cell index is out of range");
    }
    if (cell_index >= real_cell_count) {
      throw std::invalid_argument("Hydro active cells must be authoritative real cells");
    }
  }
  for (std::size_t face_index : active_set.active_faces) {
    if (face_index >= geometry.faces.size()) {
      throw std::invalid_argument("Hydro active face index is out of range");
    }
  }
}

// Thermal-state selection of the dual-energy formalism for one cell
// (Bryan et al. 1995). Returns true when the auxiliary internal energy was
// used. The conservative thermal energy E - K is authoritative unless it is
// non-positive or a fraction smaller than eta of E, where its cancellation
// error dominates; then rho e_aux supplies e, P = (gamma-1) rho e and c^2.
[[nodiscard]] bool selectDualEnergyThermalState(
    const HydroConservedStateSoa& conserved,
    std::size_t cell_index,
    const HydroConservedState& cell,
    const HydroDualEnergyPolicy& dual_energy,
    double adiabatic_index,
    HydroPrimitiveState& primitive) {
  if (dual_energy.mode == HydroDualEnergyMode::kOff || !conserved.auxiliaryInternalEnergyEnabled()) {
    return false;
  }
  const std::span<const double> auxiliary = conserved.auxiliaryInternalEnergyDensityComoving();
  if (cell_index >= auxiliary.size()) {
    return false;
  }
  const double rho_e_aux = auxiliary[cell_index];
  if (!std::isfinite(rho_e_aux) || !(rho_e_aux > 0.0)) {
    return false;
  }
  const double rho = std::max(cell.mass_density_comoving, k_small);
  const double kinetic = 0.5 *
      (cell.momentum_density_x_comoving * cell.momentum_density_x_comoving +
       cell.momentum_density_y_comoving * cell.momentum_density_y_comoving +
       cell.momentum_density_z_comoving * cell.momentum_density_z_comoving) / rho;
  const double thermal_conservative = cell.total_energy_density_comoving - kinetic;
  if (thermal_conservative > 0.0 &&
      thermal_conservative >= dual_energy.thermal_fraction_switch * cell.total_energy_density_comoving) {
    return false;
  }
  primitive.specific_internal_energy_code = rho_e_aux / rho;
  primitive.pressure_comoving = (adiabatic_index - 1.0) * rho_e_aux;
  primitive.signal_speed_squared_code = adiabatic_index * primitive.pressure_comoving / rho;
  return true;
}

[[nodiscard]] HydroPrimitiveState closedPrimitiveForCell(
    const HydroConservedStateSoa& conserved,
    std::size_t cell_index,
    double adiabatic_index,
    const HydroSourceContext& source_context,
    const HydroDualEnergyPolicy& dual_energy) {
  const HydroConservedState cell = conserved.loadCell(cell_index);
  HydroPrimitiveState primitive = HydroCoreSolver::primitiveFromConserved(cell, adiabatic_index);
  (void)selectDualEnergyThermalState(conserved, cell_index, cell, dual_energy, adiabatic_index, primitive);
  if (source_context.thermodynamic_closure != nullptr) {
    const HydroThermodynamicClosureResult closure =
        source_context.thermodynamic_closure->evaluate(
            cell_index,
            cell,
            primitive,
            source_context.update.scale_factor,
            source_context.redshift);
    if (closure.valid && closure.pressure_comoving > 0.0 &&
        closure.signal_speed_squared_code > 0.0) {
      primitive.pressure_comoving = closure.pressure_comoving;
      primitive.signal_speed_squared_code = closure.signal_speed_squared_code;
      primitive.uses_effective_ism = closure.uses_effective_ism;
    }
  }
  return primitive;
}

void fillPrimitiveCache(
    const HydroConservedStateSoa& conserved,
    const HydroActiveSetView& active_set,
    const HydroPatchGeometry& geometry,
    double adiabatic_index,
    const HydroSourceContext& source_context,
    const HydroDualEnergyPolicy& dual_energy,
    HydroPrimitiveCacheSoa& primitive_cache) {
  if (primitive_cache.size() != conserved.size()) {
    primitive_cache.resize(conserved.size());
  }
  const auto store = [&](std::size_t cell_index) {
    if (cell_index != k_invalid_cell_index) {
      primitive_cache.storeCell(
          cell_index,
          closedPrimitiveForCell(conserved, cell_index, adiabatic_index, source_context, dual_energy));
    }
  };
  for (std::size_t cell_index : active_set.active_cells) {
    store(cell_index);
  }
  for (std::size_t face_index : active_set.active_faces) {
    const HydroFace& face = geometry.faces[face_index];
    store(face.owner_cell);
    store(face.neighbor_cell);
    store(face.owner_minus_cell);
    store(face.neighbor_plus_cell);
    store(face.owner_minus2_cell);
    store(face.neighbor_plus2_cell);
  }
}

[[nodiscard]] std::size_t hydroFaceBatchLimit(
    const HydroActiveBatchPolicy& policy,
    std::size_t active_face_count) {
  const std::size_t max_active_cells = policy.max_active_cells == 0U
      ? k_hydro_automatic_active_batch_max_cells
      : policy.max_active_cells;
  if (max_active_cells > std::numeric_limits<std::size_t>::max() / k_max_faces_per_active_cell) {
    throw std::overflow_error("hydro active-batch face bound overflow");
  }
  const std::size_t face_limit = max_active_cells * k_max_faces_per_active_cell;
  return std::max<std::size_t>(1U, std::min(active_face_count, face_limit));
}

template <typename T>
[[nodiscard]] std::uint64_t vectorCapacityBytes(const std::vector<T>& values) {
  const std::uint64_t capacity = static_cast<std::uint64_t>(values.capacity());
  constexpr std::uint64_t element_bytes = static_cast<std::uint64_t>(sizeof(T));
  if (capacity > std::numeric_limits<std::uint64_t>::max() / element_bytes) {
    throw std::overflow_error("hydro scratch capacity byte count overflow");
  }
  return capacity * element_bytes;
}

[[nodiscard]] HydroConservationTotals totalsFromDelta(
    const HydroConservedState& delta,
    const HydroConservedState& reference_state,
    double cell_volume_comoving) {
  const double rho = std::max(reference_state.mass_density_comoving, k_small);
  const double old_kinetic_density = 0.5 *
      (reference_state.momentum_density_x_comoving * reference_state.momentum_density_x_comoving +
       reference_state.momentum_density_y_comoving * reference_state.momentum_density_y_comoving +
       reference_state.momentum_density_z_comoving * reference_state.momentum_density_z_comoving) /
      rho;
  const HydroConservedState updated = reference_state + delta;
  const double updated_rho = std::max(updated.mass_density_comoving, k_small);
  const double new_kinetic_density = 0.5 *
      (updated.momentum_density_x_comoving * updated.momentum_density_x_comoving +
       updated.momentum_density_y_comoving * updated.momentum_density_y_comoving +
       updated.momentum_density_z_comoving * updated.momentum_density_z_comoving) /
      updated_rho;

  HydroConservationTotals totals;
  totals.mass = delta.mass_density_comoving * cell_volume_comoving;
  totals.momentum_x = delta.momentum_density_x_comoving * cell_volume_comoving;
  totals.momentum_y = delta.momentum_density_y_comoving * cell_volume_comoving;
  totals.momentum_z = delta.momentum_density_z_comoving * cell_volume_comoving;
  totals.total_energy = delta.total_energy_density_comoving * cell_volume_comoving;
  totals.internal_energy =
      ((updated.total_energy_density_comoving - new_kinetic_density) -
       (reference_state.total_energy_density_comoving - old_kinetic_density)) *
      cell_volume_comoving;
  return totals;
}

[[nodiscard]] std::size_t fluxNeighborTargetCell(
    const HydroPatchGeometry& geometry,
    const HydroFace& face) {
  if (face.ghost_cell_slot == k_invalid_ghost_cell_slot) {
    return face.neighbor_cell;
  }
  if (face.ghost_cell_slot >= geometry.ghost_cells.size()) {
    throw std::invalid_argument("Hydro face references an invalid ghost cell slot");
  }

  const HydroGhostCell& ghost = geometry.ghost_cells[face.ghost_cell_slot];
  if (ghost.boundary_kind == HydroBoundaryKind::kPeriodic) {
    // Periodic ghosts are only reconstruction scratch.  The conservative
    // update must land on the wrapped authoritative real cell, otherwise a
    // production periodic patch leaks flux into non-authoritative ghost storage.
    return ghost.source_real_cell;
  }

  // Open and reflective physical-boundary ghosts are not authoritative cells.
  // Imported MPI ghosts are scratch rows that intentionally accumulate the
  // conservative opposite-side face delta; the workflow restores the row and
  // routes that delta back to the authoritative owner by stable gas_cell_id.
  if (ghost.boundary_kind == HydroBoundaryKind::kImportedMpi) {
    return ghost.ghost_cell;
  }
  return k_invalid_cell_index;
}

void maybeRecordFluxRegister(
    const HydroPatchGeometry& geometry,
    std::size_t face_index,
    const HydroUpdateContext& update,
    const HydroConservedState& flux_code,
    HydroFluxRegisterSink* flux_register_sink) {
  if (flux_register_sink == nullptr || geometry.flux_register_faces.empty()) {
    return;
  }
  const HydroFluxRegisterFace& register_face = geometry.flux_register_faces[face_index];
  if (register_face.role == HydroFluxRegisterFaceRole::kNone) {
    return;
  }
  const HydroConservedState oriented_flux = register_face.coarse_orientation_sign * flux_code;
  flux_register_sink->recordFaceFlux(HydroFluxRegisterRecord{
      .role = register_face.role,
      .register_key = register_face.register_key,
      .coarse_patch_id = register_face.coarse_patch_id,
      .coarse_gas_cell_id = register_face.coarse_gas_cell_id,
      .coarse_cell_index = register_face.coarse_cell_index,
      .level = register_face.level,
      .axis = register_face.axis,
      .orientation = register_face.orientation,
      .face_area_comoving = geometry.faces[face_index].area_comoving,
      .dt_code = update.dt_code,
      .flux_code = oriented_flux});
}

}  // namespace

HydroConservedState& HydroConservedState::operator+=(const HydroConservedState& rhs) {
  mass_density_comoving += rhs.mass_density_comoving;
  momentum_density_x_comoving += rhs.momentum_density_x_comoving;
  momentum_density_y_comoving += rhs.momentum_density_y_comoving;
  momentum_density_z_comoving += rhs.momentum_density_z_comoving;
  total_energy_density_comoving += rhs.total_energy_density_comoving;
  metal_mass_density_comoving += rhs.metal_mass_density_comoving;
  return *this;
}

HydroConservedState& HydroConservedState::operator-=(const HydroConservedState& rhs) {
  mass_density_comoving -= rhs.mass_density_comoving;
  momentum_density_x_comoving -= rhs.momentum_density_x_comoving;
  momentum_density_y_comoving -= rhs.momentum_density_y_comoving;
  momentum_density_z_comoving -= rhs.momentum_density_z_comoving;
  total_energy_density_comoving -= rhs.total_energy_density_comoving;
  metal_mass_density_comoving -= rhs.metal_mass_density_comoving;
  return *this;
}

HydroConservedState operator+(HydroConservedState lhs, const HydroConservedState& rhs) {
  lhs += rhs;
  return lhs;
}

HydroConservedState operator-(HydroConservedState lhs, const HydroConservedState& rhs) {
  lhs -= rhs;
  return lhs;
}

HydroConservedState operator*(double scalar, HydroConservedState state) {
  state.mass_density_comoving *= scalar;
  state.momentum_density_x_comoving *= scalar;
  state.momentum_density_y_comoving *= scalar;
  state.momentum_density_z_comoving *= scalar;
  state.total_energy_density_comoving *= scalar;
  state.metal_mass_density_comoving *= scalar;
  return state;
}

HydroConservationTotals& HydroConservationTotals::operator+=(const HydroConservationTotals& rhs) {
  mass += rhs.mass;
  momentum_x += rhs.momentum_x;
  momentum_y += rhs.momentum_y;
  momentum_z += rhs.momentum_z;
  total_energy += rhs.total_energy;
  internal_energy += rhs.internal_energy;
  metal_mass += rhs.metal_mass;
  return *this;
}

HydroConservationTotals& HydroConservationTotals::operator-=(const HydroConservationTotals& rhs) {
  mass -= rhs.mass;
  momentum_x -= rhs.momentum_x;
  momentum_y -= rhs.momentum_y;
  momentum_z -= rhs.momentum_z;
  total_energy -= rhs.total_energy;
  internal_energy -= rhs.internal_energy;
  metal_mass -= rhs.metal_mass;
  return *this;
}

HydroConservationTotals operator+(HydroConservationTotals lhs, const HydroConservationTotals& rhs) {
  lhs += rhs;
  return lhs;
}

HydroConservationTotals operator-(HydroConservationTotals lhs, const HydroConservationTotals& rhs) {
  lhs -= rhs;
  return lhs;
}

HydroConservedStateSoa::HydroConservedStateSoa(std::size_t cell_count)
    : m_mass_density_comoving(cell_count, 0.0),
      m_momentum_density_x_comoving(cell_count, 0.0),
      m_momentum_density_y_comoving(cell_count, 0.0),
      m_momentum_density_z_comoving(cell_count, 0.0),
      m_total_energy_density_comoving(cell_count, 0.0),
      m_metal_mass_density_comoving(cell_count, 0.0) {}

void HydroConservedStateSoa::resize(std::size_t cell_count) {
  if (m_auxiliary_internal_energy_enabled) {
    m_auxiliary_internal_energy_density_comoving.resize(
        cell_count, std::numeric_limits<double>::quiet_NaN());
  }
  m_mass_density_comoving.resize(cell_count, 0.0);
  m_momentum_density_x_comoving.resize(cell_count, 0.0);
  m_momentum_density_y_comoving.resize(cell_count, 0.0);
  m_momentum_density_z_comoving.resize(cell_count, 0.0);
  m_total_energy_density_comoving.resize(cell_count, 0.0);
  m_metal_mass_density_comoving.resize(cell_count, 0.0);
}

std::size_t HydroConservedStateSoa::size() const { return m_mass_density_comoving.size(); }

std::uint64_t HydroConservedStateSoa::currentSizeBytes() const {
  const std::uint64_t k_lane_count = m_auxiliary_internal_energy_enabled ? 7U : 6U;
  constexpr std::uint64_t k_lane_bytes = sizeof(double);
  const std::uint64_t cell_count = static_cast<std::uint64_t>(size());
  if (cell_count > std::numeric_limits<std::uint64_t>::max() / k_lane_count / k_lane_bytes) {
    throw std::overflow_error("hydro conserved SoA logical byte count overflow");
  }
  return cell_count * k_lane_count * k_lane_bytes;
}

std::uint64_t HydroConservedStateSoa::ownedCapacityBytes() const {
  std::uint64_t total = 0U;
  const auto add = [&total](const std::vector<double>& values) {
    const std::uint64_t bytes = vectorCapacityBytes(values);
    if (total > std::numeric_limits<std::uint64_t>::max() - bytes) {
      throw std::overflow_error("hydro conserved SoA aggregate capacity byte count overflow");
    }
    total += bytes;
  };
  add(m_mass_density_comoving);
  add(m_momentum_density_x_comoving);
  add(m_momentum_density_y_comoving);
  add(m_momentum_density_z_comoving);
  add(m_total_energy_density_comoving);
  add(m_metal_mass_density_comoving);
  add(m_auxiliary_internal_energy_density_comoving);
  return total;
}

void HydroConservedStateSoa::enableAuxiliaryInternalEnergy(bool enabled) {
  m_auxiliary_internal_energy_enabled = enabled;
  if (enabled) {
    m_auxiliary_internal_energy_density_comoving.assign(
        size(), std::numeric_limits<double>::quiet_NaN());
  } else {
    m_auxiliary_internal_energy_density_comoving.clear();
    m_auxiliary_internal_energy_density_comoving.shrink_to_fit();
  }
}

bool HydroConservedStateSoa::auxiliaryInternalEnergyEnabled() const noexcept {
  return m_auxiliary_internal_energy_enabled;
}

std::span<double> HydroConservedStateSoa::auxiliaryInternalEnergyDensityComoving() {
  return m_auxiliary_internal_energy_density_comoving;
}

std::span<const double> HydroConservedStateSoa::auxiliaryInternalEnergyDensityComoving() const {
  return m_auxiliary_internal_energy_density_comoving;
}

HydroConservedState HydroConservedStateSoa::loadCell(std::size_t cell_index) const {
  return HydroConservedState{
      .mass_density_comoving = m_mass_density_comoving.at(cell_index),
      .momentum_density_x_comoving = m_momentum_density_x_comoving.at(cell_index),
      .momentum_density_y_comoving = m_momentum_density_y_comoving.at(cell_index),
      .momentum_density_z_comoving = m_momentum_density_z_comoving.at(cell_index),
      .total_energy_density_comoving = m_total_energy_density_comoving.at(cell_index),
      .metal_mass_density_comoving = m_metal_mass_density_comoving.at(cell_index)};
}

void HydroConservedStateSoa::storeCell(std::size_t cell_index, const HydroConservedState& cell_state) {
  m_mass_density_comoving.at(cell_index) = cell_state.mass_density_comoving;
  m_momentum_density_x_comoving.at(cell_index) = cell_state.momentum_density_x_comoving;
  m_momentum_density_y_comoving.at(cell_index) = cell_state.momentum_density_y_comoving;
  m_momentum_density_z_comoving.at(cell_index) = cell_state.momentum_density_z_comoving;
  m_total_energy_density_comoving.at(cell_index) = cell_state.total_energy_density_comoving;
  m_metal_mass_density_comoving.at(cell_index) = cell_state.metal_mass_density_comoving;
}

std::span<double> HydroConservedStateSoa::massDensityComoving() { return m_mass_density_comoving; }
std::span<const double> HydroConservedStateSoa::massDensityComoving() const { return m_mass_density_comoving; }
std::span<double> HydroConservedStateSoa::momentumDensityXComoving() { return m_momentum_density_x_comoving; }
std::span<const double> HydroConservedStateSoa::momentumDensityXComoving() const { return m_momentum_density_x_comoving; }
std::span<double> HydroConservedStateSoa::momentumDensityYComoving() { return m_momentum_density_y_comoving; }
std::span<const double> HydroConservedStateSoa::momentumDensityYComoving() const { return m_momentum_density_y_comoving; }
std::span<double> HydroConservedStateSoa::momentumDensityZComoving() { return m_momentum_density_z_comoving; }
std::span<const double> HydroConservedStateSoa::momentumDensityZComoving() const { return m_momentum_density_z_comoving; }
std::span<double> HydroConservedStateSoa::totalEnergyDensityComoving() { return m_total_energy_density_comoving; }
std::span<const double> HydroConservedStateSoa::totalEnergyDensityComoving() const { return m_total_energy_density_comoving; }
std::span<double> HydroConservedStateSoa::metalMassDensityComoving() { return m_metal_mass_density_comoving; }
std::span<const double> HydroConservedStateSoa::metalMassDensityComoving() const { return m_metal_mass_density_comoving; }

HydroPrimitiveCacheSoa::HydroPrimitiveCacheSoa(std::size_t cell_count)
    : m_rho_comoving(cell_count, 0.0),
      m_vel_x_peculiar(cell_count, 0.0),
      m_vel_y_peculiar(cell_count, 0.0),
      m_vel_z_peculiar(cell_count, 0.0),
      m_pressure_comoving(cell_count, 0.0),
      m_specific_internal_energy_code(cell_count, 0.0),
      m_signal_speed_squared_code(cell_count, 0.0),
      m_uses_effective_ism(cell_count, 0U),
      m_metallicity_mass_fraction(cell_count, 0.0) {}

void HydroPrimitiveCacheSoa::resize(std::size_t cell_count) {
  m_rho_comoving.resize(cell_count, 0.0);
  m_vel_x_peculiar.resize(cell_count, 0.0);
  m_vel_y_peculiar.resize(cell_count, 0.0);
  m_vel_z_peculiar.resize(cell_count, 0.0);
  m_pressure_comoving.resize(cell_count, 0.0);
  m_specific_internal_energy_code.resize(cell_count, 0.0);
  m_signal_speed_squared_code.resize(cell_count, 0.0);
  m_uses_effective_ism.resize(cell_count, 0U);
  m_metallicity_mass_fraction.resize(cell_count, 0.0);
}

std::size_t HydroPrimitiveCacheSoa::size() const { return m_rho_comoving.size(); }

HydroPrimitiveState HydroPrimitiveCacheSoa::loadCell(std::size_t cell_index) const {
  return HydroPrimitiveState{
      .rho_comoving = m_rho_comoving.at(cell_index),
      .vel_x_peculiar = m_vel_x_peculiar.at(cell_index),
      .vel_y_peculiar = m_vel_y_peculiar.at(cell_index),
      .vel_z_peculiar = m_vel_z_peculiar.at(cell_index),
      .pressure_comoving = m_pressure_comoving.at(cell_index),
      .specific_internal_energy_code = m_specific_internal_energy_code.at(cell_index),
      .signal_speed_squared_code = m_signal_speed_squared_code.at(cell_index),
      .uses_effective_ism = m_uses_effective_ism.at(cell_index) != 0U,
      .metallicity_mass_fraction = m_metallicity_mass_fraction.at(cell_index)};
}

void HydroPrimitiveCacheSoa::storeCell(std::size_t cell_index, const HydroPrimitiveState& primitive_state) {
  m_rho_comoving.at(cell_index) = primitive_state.rho_comoving;
  m_vel_x_peculiar.at(cell_index) = primitive_state.vel_x_peculiar;
  m_vel_y_peculiar.at(cell_index) = primitive_state.vel_y_peculiar;
  m_vel_z_peculiar.at(cell_index) = primitive_state.vel_z_peculiar;
  m_pressure_comoving.at(cell_index) = primitive_state.pressure_comoving;
  m_specific_internal_energy_code.at(cell_index) = primitive_state.specific_internal_energy_code;
  m_signal_speed_squared_code.at(cell_index) = primitive_state.signal_speed_squared_code;
  m_uses_effective_ism.at(cell_index) = primitive_state.uses_effective_ism ? 1U : 0U;
  m_metallicity_mass_fraction.at(cell_index) = primitive_state.metallicity_mass_fraction;
}

HydroConservedState ComovingGravityExpansionSource::sourceForCell(
    std::size_t cell_index,
    const HydroConservedState& conserved,
    const HydroPrimitiveState& primitive,
    const HydroSourceContext& context) const {
  const std::array<double, 3> g = gravityAtCell(
      cell_index,
      context.gravity_accel_x_peculiar,
      context.gravity_accel_y_peculiar,
      context.gravity_accel_z_peculiar);

  if (!std::isfinite(context.update.scale_factor) || context.update.scale_factor <= 0.0) {
    throw std::invalid_argument("comoving gravity source requires finite scale_factor > 0");
  }
  if (!std::isfinite(context.update.hubble_rate_code) || context.update.hubble_rate_code < 0.0) {
    throw std::invalid_argument("comoving gravity source requires finite hubble_rate_code >= 0");
  }
  // TreePM returns the scale-free comoving particle kernel A. The conserved
  // gas momentum stores rho_c u with physical peculiar velocity u, whose
  // equation is du/dt + H u = A/a^2. Keep the scale-factor response here so
  // particles and gas share one TreePM normalization without double counting.
  const double inverse_scale_factor_squared =
      1.0 / (context.update.scale_factor * context.update.scale_factor);
  const std::array<double, 3> peculiar_acceleration{
      g[0] * inverse_scale_factor_squared,
      g[1] * inverse_scale_factor_squared,
      g[2] * inverse_scale_factor_squared,
  };
  const double hubble_rate = context.update.hubble_rate_code;
  const double internal_energy_density = std::max(
      conserved.total_energy_density_comoving -
          0.5 * (conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
                 conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
                 conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving) /
              std::max(conserved.mass_density_comoving, k_small),
      0.0);
  const double kinetic_energy_density = std::max(
      conserved.total_energy_density_comoving - internal_energy_density,
      0.0);

  HydroConservedState source;
  source.mass_density_comoving = 0.0;
  source.momentum_density_x_comoving =
      primitive.rho_comoving * peculiar_acceleration[0] -
      hubble_rate * conserved.momentum_density_x_comoving;
  source.momentum_density_y_comoving =
      primitive.rho_comoving * peculiar_acceleration[1] -
      hubble_rate * conserved.momentum_density_y_comoving;
  source.momentum_density_z_comoving =
      primitive.rho_comoving * peculiar_acceleration[2] -
      hubble_rate * conserved.momentum_density_z_comoving;

  const double work_gravity = primitive.rho_comoving *
      (primitive.vel_x_peculiar * peculiar_acceleration[0] +
       primitive.vel_y_peculiar * peculiar_acceleration[1] +
       primitive.vel_z_peculiar * peculiar_acceleration[2]);
  // For rho_c, rho_c u, and rho_c(e + u^2/2), homogeneous expansion gives
  // -H rho_c u^2 - 3 H P. This is -2 H E for gamma=5/3, while the pressure
  // form remains correct for every gamma supported by HydroCoreSolver.
  source.total_energy_density_comoving =
      work_gravity - hubble_rate *
          (2.0 * kinetic_energy_density + 3.0 * primitive.pressure_comoving);
  source.metal_mass_density_comoving = 0.0;
  return source;
}

HydroCoreSolver::HydroCoreSolver(double adiabatic_index, HydroNumericsPolicy numerics)
    : m_adiabatic_index(adiabatic_index), m_numerics(numerics) {
  if (m_adiabatic_index <= 1.0) {
    throw std::invalid_argument("HydroCoreSolver requires adiabatic_index > 1");
  }
  if (m_numerics.dual_energy.mode != HydroDualEnergyMode::kOff &&
      (!std::isfinite(m_numerics.dual_energy.thermal_fraction_switch) ||
       m_numerics.dual_energy.thermal_fraction_switch <= 0.0 ||
       m_numerics.dual_energy.thermal_fraction_switch >= 1.0)) {
    throw std::invalid_argument("HydroCoreSolver dual-energy thermal_fraction_switch must be in (0, 1)");
  }
}

double HydroCoreSolver::adiabaticIndex() const { return m_adiabatic_index; }

const HydroNumericsPolicy& HydroCoreSolver::numericsPolicy() const noexcept { return m_numerics; }

HydroActiveBatchPolicy HydroCoreSolver::requiredWorkspaceFeatures(
    HydroActiveBatchPolicy batch_policy) const noexcept {
  batch_policy.ctu_transverse_workspace = batch_policy.ctu_transverse_workspace ||
      m_numerics.integration_scheme == HydroIntegrationScheme::kCornerTransportUpwind;
  batch_policy.dual_energy_workspace = batch_policy.dual_energy_workspace ||
      m_numerics.dual_energy.mode != HydroDualEnergyMode::kOff;
  return batch_policy;
}

std::string_view hydroIntegrationSchemeToString(HydroIntegrationScheme scheme) {
  switch (scheme) {
    case HydroIntegrationScheme::kMusclHancock:
      return "muscl_hancock";
    case HydroIntegrationScheme::kCornerTransportUpwind:
      return "ctu";
  }
  return "unknown";
}

HydroIntegrationScheme hydroIntegrationSchemeFromString(std::string_view name) {
  if (name == "muscl_hancock") {
    return HydroIntegrationScheme::kMusclHancock;
  }
  if (name == "ctu" || name == "corner_transport_upwind") {
    return HydroIntegrationScheme::kCornerTransportUpwind;
  }
  throw std::invalid_argument("Unsupported hydro integration scheme: " + std::string(name));
}

std::string_view hydroDualEnergyModeToString(HydroDualEnergyMode mode) {
  switch (mode) {
    case HydroDualEnergyMode::kOff:
      return "off";
    case HydroDualEnergyMode::kThermalFraction:
      return "thermal_fraction";
  }
  return "unknown";
}

HydroDualEnergyMode hydroDualEnergyModeFromString(std::string_view name) {
  if (name == "off") {
    return HydroDualEnergyMode::kOff;
  }
  if (name == "thermal_fraction") {
    return HydroDualEnergyMode::kThermalFraction;
  }
  throw std::invalid_argument("Unsupported hydro dual-energy mode: " + std::string(name));
}

namespace {
// Maximum staged primitive stencil cells per face: owner_minus2, owner_minus,
// owner, neighbor, neighbor_plus, neighbor_plus2.
constexpr std::size_t k_max_stencil_cells_per_face = 6U;

struct HydroScratchCapacityPlan {
  std::size_t touched_staging = 0;
  std::size_t touched_unique = 0;
  std::size_t source_active = 0;
  std::size_t full_active = 0;
  std::size_t face_batch = 0;
  std::size_t stencil_staging = 0;
  std::size_t stencil_unique = 0;
  std::size_t ctu_cells = 0;
  std::size_t dual_energy_cells = 0;
};

[[nodiscard]] HydroScratchCapacityPlan planHydroScratch(
    std::size_t active_cells, std::size_t active_faces,
    std::size_t real_cells, std::size_t total_cells,
    const HydroActiveBatchPolicy& policy) {
  if (real_cells > total_cells || active_cells > real_cells) {
    throw std::invalid_argument("hydro scratch extent is inconsistent with geometry");
  }
  // The solver appends raw face incidences before sorting/deduplicating.
  // Reserve their raw cardinality, not the smaller unique-cell cardinality.
  const std::size_t touched = core::checkedSizeAdd(
      active_cells, core::checkedSizeMultiply(active_faces, 2U,
          "hydro touched-cell bound"), "hydro touched-cell bound");
  const std::size_t face_batch = hydroFaceBatchLimit(policy, active_faces);
  const std::size_t stencil = core::checkedSizeMultiply(face_batch, k_max_stencil_cells_per_face,
      "hydro primitive stencil bound");
  const std::size_t touched_unique = std::min(total_cells, touched);
  return {touched, touched_unique, active_cells,
      std::min(total_cells, touched), face_batch, stencil,
      std::min(total_cells, stencil),
      policy.ctu_transverse_workspace ? touched_unique : 0U,
      policy.dual_energy_workspace ? touched_unique : 0U};
}

template <typename T>
void reserveHydroVector(std::vector<T>& values, std::size_t count) {
  if (count > values.max_size()) throw std::length_error("hydro scratch capacity exceeds max_size");
  if (count > values.capacity()) {
    values.reserve(count);
    if (values.capacity() > count) {
      throw std::length_error("hydro scratch allocator exceeded admitted capacity");
    }
  }
}

[[nodiscard]] std::uint64_t scratchPlanBytes(const HydroScratchCapacityPlan& plan) {
  std::uint64_t bytes = 0U;
  const auto add = [&bytes](std::size_t count, std::size_t width) {
    bytes = core::checkedMemoryBytesAdd(bytes,
        static_cast<std::uint64_t>(core::checkedSizeMultiply(count, width,
            "hydro scratch physical bound")), "hydro scratch physical bound");
  };
  add(plan.touched_unique, sizeof(HydroConservedState));
  add(plan.face_batch, 2U * sizeof(HydroPrimitiveState) + sizeof(HydroConservedState));
  add(plan.touched_staging, sizeof(std::size_t));
  add(plan.source_active, sizeof(std::size_t));
  add(plan.stencil_staging, sizeof(std::size_t));
  add(plan.stencil_staging, sizeof(HydroPrimitiveState));
  add(plan.full_active, sizeof(std::size_t));
  add(plan.ctu_cells, sizeof(std::array<HydroConservedState, 3>));
  add(plan.dual_energy_cells, 2U * sizeof(double));
  return bytes;
}
}  // namespace

std::uint64_t HydroScratchBuffers::workspaceBytes(
    std::size_t active_cell_count, std::size_t active_face_count,
    std::size_t real_cell_count, std::size_t total_cell_count,
    const HydroActiveBatchPolicy& batch_policy) {
  return scratchPlanBytes(planHydroScratch(active_cell_count, active_face_count,
      real_cell_count, total_cell_count, batch_policy));
}

std::uint64_t HydroScratchBuffers::prepareForActiveSet(
    std::size_t active_cell_count, std::size_t active_face_count,
    std::size_t real_cell_count, std::size_t total_cell_count,
    const HydroActiveBatchPolicy& batch_policy) {
  const auto plan = planHydroScratch(active_cell_count, active_face_count,
      real_cell_count, total_cell_count, batch_policy);
  reserveHydroVector(cell_delta, plan.touched_unique);
  reserveHydroVector(left_states, plan.face_batch);
  reserveHydroVector(right_states, plan.face_batch);
  reserveHydroVector(fluxes, plan.face_batch);
  reserveHydroVector(touched_cells, plan.touched_staging);
  reserveHydroVector(source_active_cells, plan.source_active);
  reserveHydroVector(primitive_stencil_cells, plan.stencil_staging);
  reserveHydroVector(primitive_stencil_states, plan.stencil_staging);
  reserveHydroVector(full_active_cells, plan.full_active);
  reserveHydroVector(ctu_transverse_delta, plan.ctu_cells);
  reserveHydroVector(aux_internal_energy_advection, plan.dual_energy_cells);
  reserveHydroVector(aux_velocity_divergence, plan.dual_energy_cells);
  // full_active_faces is a compatibility lane unused by the sparse solver.
  // Its pre-existing capacity is still included in ownedCapacityBytes().
  const std::uint64_t actual = ownedCapacityBytes();
  high_water.capacity_bytes = std::max(high_water.capacity_bytes, actual);
  return actual;
}

void HydroScratchBuffers::resizeFaceBatch(
    std::size_t active_cell_batch_count,
    std::size_t active_face_batch_count) {
  if (active_face_batch_count > std::numeric_limits<std::size_t>::max() / k_max_stencil_cells_per_face) {
    throw std::overflow_error("hydro primitive stencil batch size overflow");
  }
  const std::size_t primitive_stencil_capacity = active_face_batch_count * k_max_stencil_cells_per_face;
  left_states.resize(active_face_batch_count);
  right_states.resize(active_face_batch_count);
  fluxes.resize(active_face_batch_count);
  primitive_stencil_cells.reserve(primitive_stencil_capacity);
  primitive_stencil_states.reserve(primitive_stencil_capacity);
  high_water.active_cell_batch_capacity = std::max(
      high_water.active_cell_batch_capacity, active_cell_batch_count);
  high_water.face_batch_capacity = std::max(
      high_water.face_batch_capacity, left_states.capacity());
  high_water.primitive_stencil_capacity = std::max(
      high_water.primitive_stencil_capacity, primitive_stencil_states.capacity());
  high_water.touched_cell_capacity = std::max(
      high_water.touched_cell_capacity, touched_cells.capacity());
  high_water.capacity_bytes = std::max(high_water.capacity_bytes, ownedCapacityBytes());
}

std::uint64_t HydroScratchBuffers::ownedCapacityBytes() const {
  std::uint64_t total = 0U;
  const auto add = [&total](std::uint64_t value) {
    if (total > std::numeric_limits<std::uint64_t>::max() - value) {
      throw std::overflow_error("hydro scratch aggregate capacity byte count overflow");
    }
    total += value;
  };
  add(vectorCapacityBytes(cell_delta));
  add(vectorCapacityBytes(left_states));
  add(vectorCapacityBytes(right_states));
  add(vectorCapacityBytes(fluxes));
  add(vectorCapacityBytes(touched_cells));
  add(vectorCapacityBytes(source_active_cells));
  add(vectorCapacityBytes(primitive_stencil_cells));
  add(vectorCapacityBytes(primitive_stencil_states));
  add(vectorCapacityBytes(full_active_cells));
  add(vectorCapacityBytes(full_active_faces));
  add(vectorCapacityBytes(ctu_transverse_delta));
  add(vectorCapacityBytes(aux_internal_energy_advection));
  add(vectorCapacityBytes(aux_velocity_divergence));
  return total;
}

HydroConservedState HydroCoreSolver::conservedFromPrimitive(
    const HydroPrimitiveState& primitive,
    double adiabatic_index) {
  if (primitive.rho_comoving <= 0.0) {
    throw std::invalid_argument("Primitive state requires rho_comoving > 0");
  }
  if (primitive.pressure_comoving <= 0.0) {
    throw std::invalid_argument("Primitive state requires pressure_comoving > 0");
  }

  const double velocity_squared =
      primitive.vel_x_peculiar * primitive.vel_x_peculiar +
      primitive.vel_y_peculiar * primitive.vel_y_peculiar +
      primitive.vel_z_peculiar * primitive.vel_z_peculiar;

  HydroConservedState conserved;
  conserved.mass_density_comoving = primitive.rho_comoving;
  conserved.momentum_density_x_comoving = primitive.rho_comoving * primitive.vel_x_peculiar;
  conserved.momentum_density_y_comoving = primitive.rho_comoving * primitive.vel_y_peculiar;
  conserved.momentum_density_z_comoving = primitive.rho_comoving * primitive.vel_z_peculiar;
  const double specific_internal_energy = primitive.specific_internal_energy_code > 0.0
      ? primitive.specific_internal_energy_code
      : primitive.pressure_comoving / ((adiabatic_index - 1.0) * primitive.rho_comoving);
  conserved.total_energy_density_comoving =
      primitive.rho_comoving * specific_internal_energy + 0.5 * primitive.rho_comoving * velocity_squared;
  conserved.metal_mass_density_comoving = primitive.rho_comoving *
      std::clamp(primitive.metallicity_mass_fraction, 0.0, 1.0);
  return conserved;
}

HydroPrimitiveState HydroCoreSolver::primitiveFromConserved(
    const HydroConservedState& conserved,
    double adiabatic_index) {
  if (conserved.mass_density_comoving <= 0.0) {
    throw std::invalid_argument("Conserved state requires mass_density_comoving > 0");
  }

  const double inv_rho = 1.0 / conserved.mass_density_comoving;
  HydroPrimitiveState primitive;
  primitive.rho_comoving = conserved.mass_density_comoving;
  primitive.vel_x_peculiar = conserved.momentum_density_x_comoving * inv_rho;
  primitive.vel_y_peculiar = conserved.momentum_density_y_comoving * inv_rho;
  primitive.vel_z_peculiar = conserved.momentum_density_z_comoving * inv_rho;

  const double kinetic_density = 0.5 *
      (conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
       conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
       conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving) *
      inv_rho;
  const double internal_density = conserved.total_energy_density_comoving - kinetic_density;
  primitive.specific_internal_energy_code = std::max(internal_density * inv_rho, k_small);
  primitive.pressure_comoving = std::max((adiabatic_index - 1.0) * internal_density, k_small);
  primitive.signal_speed_squared_code = std::max(
      adiabatic_index * primitive.pressure_comoving * inv_rho, k_small);
  primitive.uses_effective_ism = false;
  primitive.metallicity_mass_fraction = std::clamp(
      conserved.metal_mass_density_comoving * inv_rho, 0.0, 1.0);
  return primitive;
}

HydroConservationTotals HydroCoreSolver::conservationTotals(
    const HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    std::span<const std::size_t> cells) {
  if (geometry.cell_volume_comoving <= 0.0) {
    throw std::invalid_argument("Hydro conservation totals require positive cell_volume_comoving");
  }

  HydroConservationTotals totals;
  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  for (std::size_t cell_index : cells) {
    if (cell_index >= conserved.size()) {
      throw std::out_of_range("Hydro conservation totals cell index is out of range");
    }
    if (cell_index >= real_cell_count) {
      throw std::invalid_argument("Hydro conservation totals require authoritative real cells");
    }
    const HydroConservedState state = conserved.loadCell(cell_index);
    const double rho = std::max(state.mass_density_comoving, k_small);
    const double kinetic_density = 0.5 *
        (state.momentum_density_x_comoving * state.momentum_density_x_comoving +
         state.momentum_density_y_comoving * state.momentum_density_y_comoving +
         state.momentum_density_z_comoving * state.momentum_density_z_comoving) /
        rho;
    totals.mass += state.mass_density_comoving * geometry.cell_volume_comoving;
    totals.momentum_x += state.momentum_density_x_comoving * geometry.cell_volume_comoving;
    totals.momentum_y += state.momentum_density_y_comoving * geometry.cell_volume_comoving;
    totals.momentum_z += state.momentum_density_z_comoving * geometry.cell_volume_comoving;
    totals.total_energy += state.total_energy_density_comoving * geometry.cell_volume_comoving;
    totals.internal_energy +=
        (state.total_energy_density_comoving - kinetic_density) * geometry.cell_volume_comoving;
    totals.metal_mass += state.metal_mass_density_comoving * geometry.cell_volume_comoving;
  }
  return totals;
}

void HydroCoreSolver::advancePatch(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  HydroScratchBuffers scratch;
  advancePatchWithScratch(
      conserved,
      geometry,
      update,
      reconstruction,
      riemann_solver,
      source_terms,
      source_context,
      scratch,
      nullptr,
      profile,
      flux_register_sink);
}

void HydroCoreSolver::advancePatchWithScratch(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context,
    HydroScratchBuffers& scratch,
    HydroPrimitiveCacheSoa* primitive_cache,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  scratch.full_active_cells.resize(conserved.size());
  for (std::size_t i = 0; i < scratch.full_active_cells.size(); ++i) {
    scratch.full_active_cells[i] = i;
  }
  scratch.full_active_faces.resize(geometry.faces.size());
  for (std::size_t i = 0; i < scratch.full_active_faces.size(); ++i) {
    scratch.full_active_faces[i] = i;
  }
  advancePatchActiveSetWithScratch(
      conserved,
      geometry,
      HydroActiveSetView{.active_cells = scratch.full_active_cells, .active_faces = scratch.full_active_faces},
      update,
      reconstruction,
      riemann_solver,
      source_terms,
      source_context,
      scratch,
      primitive_cache,
      profile,
      flux_register_sink);
}

void HydroCoreSolver::advancePatchActiveSet(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroActiveSetView& active_set,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  HydroScratchBuffers scratch;
  advancePatchActiveSetWithScratch(
      conserved,
      geometry,
      active_set,
      update,
      reconstruction,
      riemann_solver,
      source_terms,
      source_context,
      scratch,
      nullptr,
      profile,
      flux_register_sink);
}

void HydroCoreSolver::advancePatchActiveSetWithScratch(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroActiveSetView& active_set,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context,
    HydroScratchBuffers& scratch,
    HydroPrimitiveCacheSoa* primitive_cache,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  advancePatchActiveSetBatchedWithScratch(
      conserved,
      geometry,
      active_set,
      update,
      reconstruction,
      riemann_solver,
      source_terms,
      source_context,
      scratch,
      primitive_cache,
      requiredWorkspaceFeatures(HydroActiveBatchPolicy{}),
      profile,
      flux_register_sink);
}

namespace {

constexpr std::size_t k_invalid_touched_slot = static_cast<std::size_t>(-1);

[[nodiscard]] std::size_t findTouchedSlot(const HydroScratchBuffers& scratch, std::size_t cell_index) {
  const auto it = std::lower_bound(
      scratch.touched_cells.begin(), scratch.touched_cells.end(), cell_index);
  if (it == scratch.touched_cells.end() || *it != cell_index) {
    return k_invalid_touched_slot;
  }
  return static_cast<std::size_t>(std::distance(scratch.touched_cells.begin(), it));
}

[[nodiscard]] std::size_t requireTouchedSlot(const HydroScratchBuffers& scratch, std::size_t cell_index) {
  const std::size_t slot = findTouchedSlot(scratch, cell_index);
  if (slot == k_invalid_touched_slot) {
    throw std::logic_error("hydro sparse touched-cell lookup missed a touched cell");
  }
  return slot;
}

[[nodiscard]] double faceNormalVelocity(const HydroPrimitiveState& state, const HydroFace& face) {
  return dot3(state.vel_x_peculiar, state.vel_y_peculiar, state.vel_z_peculiar,
      face.normal_x, face.normal_y, face.normal_z);
}

[[nodiscard]] std::size_t faceAxisIndex(HydroFaceAxis axis) {
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

// Real cell whose CTU transverse increment describes the given face side.
// Periodic ghosts map to their wrapped source cell; open/reflective/imported
// ghosts carry no transverse information (returns k_invalid_cell_index).
[[nodiscard]] std::size_t transverseSourceCell(
    const HydroPatchGeometry& geometry,
    std::size_t cell_index,
    std::size_t real_cell_count) {
  if (cell_index < real_cell_count) {
    return cell_index;
  }
  if (cell_index == k_invalid_cell_index || cell_index - real_cell_count >= geometry.ghost_cells.size()) {
    return k_invalid_cell_index;
  }
  const HydroGhostCell& ghost = geometry.ghost_cells[cell_index - real_cell_count];
  if (ghost.boundary_kind == HydroBoundaryKind::kPeriodic) {
    return ghost.source_real_cell;
  }
  return k_invalid_cell_index;
}

[[nodiscard]] bool admissiblePrimitive(const HydroPrimitiveState& state, double rho_floor, double pressure_floor) {
  return std::isfinite(state.rho_comoving) && std::isfinite(state.pressure_comoving) &&
      std::isfinite(state.vel_x_peculiar) && std::isfinite(state.vel_y_peculiar) &&
      std::isfinite(state.vel_z_peculiar) && state.rho_comoving >= rho_floor &&
      state.pressure_comoving >= pressure_floor && state.specific_internal_energy_code > 0.0 &&
      state.signal_speed_squared_code > 0.0;
}

// Apply half of the transverse conservative increments to a face state.
// Ideal-gas states are reconverted with the solver EOS; effective-ISM states
// keep their barotropic closure linearized at the predicted state
// (P' = P + c^2 (rho' - rho), c'^2 = c^2).
[[nodiscard]] bool applyTransverseCorrection(
    HydroPrimitiveState& state,
    const HydroConservedState& half_transverse,
    double adiabatic_index,
    double rho_floor,
    double pressure_floor) {
  const HydroConservedState corrected =
      HydroCoreSolver::conservedFromPrimitive(state, adiabatic_index) + half_transverse;
  if (!std::isfinite(corrected.mass_density_comoving) || corrected.mass_density_comoving < rho_floor) {
    return false;
  }
  HydroPrimitiveState candidate = HydroCoreSolver::primitiveFromConserved(corrected, adiabatic_index);
  if (state.uses_effective_ism) {
    candidate.pressure_comoving = state.pressure_comoving +
        state.signal_speed_squared_code * (candidate.rho_comoving - state.rho_comoving);
    candidate.signal_speed_squared_code = state.signal_speed_squared_code;
    candidate.uses_effective_ism = true;
  }
  const double internal_density = corrected.total_energy_density_comoving -
      0.5 * (corrected.momentum_density_x_comoving * corrected.momentum_density_x_comoving +
             corrected.momentum_density_y_comoving * corrected.momentum_density_y_comoving +
             corrected.momentum_density_z_comoving * corrected.momentum_density_z_comoving) /
          corrected.mass_density_comoving;
  if (!(internal_density > 0.0) || !admissiblePrimitive(candidate, rho_floor, pressure_floor)) {
    return false;
  }
  state = candidate;
  return true;
}

}  // namespace

void HydroCoreSolver::advancePatchActiveSetBatchedWithScratch(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroActiveSetView& active_set,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context_in,
    HydroScratchBuffers& scratch,
    HydroPrimitiveCacheSoa* primitive_cache,
    const HydroActiveBatchPolicy& batch_policy,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  const std::uint64_t limiter_before = reconstruction.limiterClipCount();
  const std::uint64_t positivity_before = reconstruction.positivityFallbackCount();
  const std::uint64_t degradation_before = reconstruction.stencilDegradationCount();
  const std::uint64_t riemann_fallback_before = riemann_solver.fallbackCount();

  const bool use_ctu = m_numerics.integration_scheme == HydroIntegrationScheme::kCornerTransportUpwind;
  const bool use_dual_energy = m_numerics.dual_energy.mode != HydroDualEnergyMode::kOff;
  if (use_ctu && !batch_policy.ctu_transverse_workspace) {
    throw std::invalid_argument(
        "CTU hydro advance requires HydroActiveBatchPolicy::ctu_transverse_workspace; "
        "plan the batch with HydroCoreSolver::requiredWorkspaceFeatures()");
  }
  if (use_dual_energy && !batch_policy.dual_energy_workspace) {
    throw std::invalid_argument(
        "dual-energy hydro advance requires HydroActiveBatchPolicy::dual_energy_workspace; "
        "plan the batch with HydroCoreSolver::requiredWorkspaceFeatures()");
  }

  validateAdvanceInputs(conserved, geometry, update, source_context_in, m_adiabatic_index);
  validateActiveSet(active_set, conserved, geometry);
  if (use_dual_energy && !conserved.auxiliaryInternalEnergyEnabled()) {
    throw std::invalid_argument(
        "dual-energy hydro advance requires the auxiliary internal-energy lane "
        "(HydroConservedStateSoa::enableAuxiliaryInternalEnergy) initialized from the canonical "
        "specific internal energy");
  }
  // Source terms integrate densities over the authoritative cell volume
  // (audit F7); the copy is a handful of spans/pointers, not field data.
  HydroSourceContext source_context = source_context_in;
  source_context.cell_volume_comoving = geometry.cell_volume_comoving;
  HydroConservationReport conservation_report;

  scratch.source_active_cells.assign(active_set.active_cells.begin(), active_set.active_cells.end());
  std::sort(scratch.source_active_cells.begin(), scratch.source_active_cells.end());
  scratch.source_active_cells.erase(
      std::unique(scratch.source_active_cells.begin(), scratch.source_active_cells.end()),
      scratch.source_active_cells.end());

  scratch.touched_cells.clear();
  if (active_set.active_faces.size() > (std::numeric_limits<std::size_t>::max() -
      active_set.active_cells.size()) / 2U) {
    throw std::overflow_error("hydro touched-cell staging count overflow");
  }
  scratch.touched_cells.reserve(active_set.active_faces.size() * 2U + active_set.active_cells.size());
  scratch.touched_cells.insert(
      scratch.touched_cells.end(), active_set.active_cells.begin(), active_set.active_cells.end());
  for (std::size_t face_index : active_set.active_faces) {
    const HydroFace& face = geometry.faces[face_index];
    scratch.touched_cells.push_back(face.owner_cell);
    const std::size_t neighbor_target = fluxNeighborTargetCell(geometry, face);
    if (neighbor_target != k_invalid_cell_index) {
      scratch.touched_cells.push_back(neighbor_target);
    }
  }
  std::sort(scratch.touched_cells.begin(), scratch.touched_cells.end());
  scratch.touched_cells.erase(
      std::unique(scratch.touched_cells.begin(), scratch.touched_cells.end()),
      scratch.touched_cells.end());
  scratch.cell_delta.assign(scratch.touched_cells.size(), HydroConservedState{});
  if (use_ctu) {
    scratch.ctu_transverse_delta.assign(scratch.touched_cells.size(), std::array<HydroConservedState, 3>{});
  }
  if (use_dual_energy) {
    scratch.aux_internal_energy_advection.assign(scratch.touched_cells.size(), 0.0);
    scratch.aux_velocity_divergence.assign(scratch.touched_cells.size(), 0.0);
  }

  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  scratch.full_active_cells.clear();
  scratch.full_active_cells.reserve(scratch.touched_cells.size());
  for (const std::size_t cell_index : scratch.touched_cells) {
    if (cell_index < real_cell_count) {
      scratch.full_active_cells.push_back(cell_index);
    }
  }
  conservation_report.before = conservationTotals(conserved, geometry, scratch.full_active_cells);
  conservation_report.cell_count = static_cast<std::uint64_t>(scratch.full_active_cells.size());

  // Compatibility callers may still provide a full cache. Production M2A
  // paths pass nullptr and construct only batch-local primitive stencils.
  if (primitive_cache != nullptr) {
    fillPrimitiveCache(
        conserved, active_set, geometry, m_adiabatic_index, source_context, m_numerics.dual_energy,
        *primitive_cache);
  }

  const std::size_t face_batch_limit = hydroFaceBatchLimit(batch_policy, active_set.active_faces.size());
  scratch.resizeFaceBatch(
      batch_policy.max_active_cells == 0U ? k_hydro_automatic_active_batch_max_cells : batch_policy.max_active_cells,
      face_batch_limit);

  const auto total_start = std::chrono::steady_clock::now();
  double reconstruct_ms = 0.0;
  double riemann_ms = 0.0;
  double accumulate_ms = 0.0;
  std::uint64_t riemann_solve_count = 0U;
  std::uint64_t ctu_corrected = 0U;
  std::uint64_t ctu_degraded = 0U;
  std::uint64_t ctu_inadmissible = 0U;
  const double flux_scale = update.dt_code / (update.scale_factor * geometry.cell_volume_comoving);
  // Floors used only to judge admissibility of CTU-corrected face states.
  constexpr double k_state_floor = 1.0e-30;

  // Stage primitives for one face batch and fill left/right interface states
  // with the selected reconstruction (normal-only Hancock prediction).
  const auto reconstructBatch = [&](std::size_t batch_begin, std::size_t batch_count) {
    if (primitive_cache == nullptr) {
      scratch.primitive_stencil_cells.clear();
      for (std::size_t slot = 0; slot < batch_count; ++slot) {
        const HydroFace& face = geometry.faces[active_set.active_faces[batch_begin + slot]];
        const std::size_t candidates[k_max_stencil_cells_per_face] = {
            face.owner_minus2_cell, face.owner_minus_cell, face.owner_cell,
            face.neighbor_cell, face.neighbor_plus_cell, face.neighbor_plus2_cell};
        for (const std::size_t cell : candidates) {
          if (cell != k_invalid_cell_index) {
            scratch.primitive_stencil_cells.push_back(cell);
          }
        }
      }
      std::sort(scratch.primitive_stencil_cells.begin(), scratch.primitive_stencil_cells.end());
      scratch.primitive_stencil_cells.erase(
          std::unique(scratch.primitive_stencil_cells.begin(), scratch.primitive_stencil_cells.end()),
          scratch.primitive_stencil_cells.end());
      scratch.primitive_stencil_states.resize(scratch.primitive_stencil_cells.size());
      for (std::size_t slot = 0; slot < scratch.primitive_stencil_cells.size(); ++slot) {
        scratch.primitive_stencil_states[slot] = closedPrimitiveForCell(
            conserved,
            scratch.primitive_stencil_cells[slot],
            m_adiabatic_index,
            source_context,
            m_numerics.dual_energy);
      }
      scratch.high_water.primitive_stencil_capacity = std::max(
          scratch.high_water.primitive_stencil_capacity,
          scratch.primitive_stencil_states.capacity());
    }

    const auto primitive_for = [&](std::size_t cell_index) -> const HydroPrimitiveState* {
      if (cell_index == k_invalid_cell_index) {
        return nullptr;
      }
      const auto it = std::lower_bound(
          scratch.primitive_stencil_cells.begin(),
          scratch.primitive_stencil_cells.end(),
          cell_index);
      if (it == scratch.primitive_stencil_cells.end() || *it != cell_index) {
        throw std::logic_error("hydro primitive stencil lookup missed a staged cell");
      }
      return &scratch.primitive_stencil_states[static_cast<std::size_t>(
          std::distance(scratch.primitive_stencil_cells.begin(), it))];
    };

    for (std::size_t slot = 0; slot < batch_count; ++slot) {
      const std::size_t face_index = active_set.active_faces[batch_begin + slot];
      const HydroFace& face = geometry.faces[face_index];
      bool consumed = false;
      if (primitive_cache != nullptr) {
        consumed = reconstruction.reconstructFaceFromCache(
            *primitive_cache, face, scratch.left_states[slot], scratch.right_states[slot]);
      } else {
        const HydroFaceStencilStates stencil{
            .owner_minus2 = primitive_for(face.owner_minus2_cell),
            .owner_minus = primitive_for(face.owner_minus_cell),
            .owner = primitive_for(face.owner_cell),
            .neighbor = primitive_for(face.neighbor_cell),
            .neighbor_plus = primitive_for(face.neighbor_plus_cell),
            .neighbor_plus2 = primitive_for(face.neighbor_plus2_cell)};
        consumed = reconstruction.reconstructFaceFromStencil(
            stencil, face, scratch.left_states[slot], scratch.right_states[slot]);
      }
      if (!consumed) {
        reconstruction.reconstructFace(
            conserved,
            face,
            scratch.left_states[slot],
            scratch.right_states[slot],
            m_adiabatic_index);
      }
    }
  };

  // Half of the CTU transverse increment of `cell` excluding `face_axis`.
  const auto halfTransverse = [&](std::size_t cell, std::size_t face_axis,
                                  HydroConservedState& out) -> bool {
    const std::size_t source_cell = transverseSourceCell(geometry, cell, real_cell_count);
    if (source_cell == k_invalid_cell_index) {
      return false;
    }
    const std::size_t touched = findTouchedSlot(scratch, source_cell);
    if (touched == k_invalid_touched_slot) {
      return false;
    }
    out = HydroConservedState{};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (axis != face_axis) {
        out += scratch.ctu_transverse_delta[touched][axis];
      }
    }
    out = 0.5 * out;
    return true;
  };

  // CTU pass 1: first-order-in-time transverse fluxes from normal-predicted
  // states. These fluxes never touch authoritative state or flux registers.
  if (use_ctu) {
    for (std::size_t batch_begin = 0; batch_begin < active_set.active_faces.size();
         batch_begin += face_batch_limit) {
      const std::size_t batch_count = std::min(
          face_batch_limit, active_set.active_faces.size() - batch_begin);
      scratch.left_states.resize(batch_count);
      scratch.right_states.resize(batch_count);
      scratch.fluxes.resize(batch_count);
      const auto reconstruct_start = std::chrono::steady_clock::now();
      reconstructBatch(batch_begin, batch_count);
      const auto riemann_start = std::chrono::steady_clock::now();
      for (std::size_t slot = 0; slot < batch_count; ++slot) {
        const std::size_t face_index = active_set.active_faces[batch_begin + slot];
        scratch.fluxes[slot] = riemann_solver.computeFlux(
            scratch.left_states[slot], scratch.right_states[slot],
            geometry.faces[face_index], m_adiabatic_index);
      }
      riemann_solve_count += static_cast<std::uint64_t>(batch_count);
      const auto accumulate_start = std::chrono::steady_clock::now();
      for (std::size_t slot = 0; slot < batch_count; ++slot) {
        const HydroFace& face = geometry.faces[active_set.active_faces[batch_begin + slot]];
        const std::size_t axis = faceAxisIndex(face.axis);
        const HydroConservedState face_delta = (flux_scale * face.area_comoving) * scratch.fluxes[slot];
        scratch.ctu_transverse_delta[requireTouchedSlot(scratch, face.owner_cell)][axis] -= face_delta;
        const std::size_t neighbor_target = fluxNeighborTargetCell(geometry, face);
        if (neighbor_target != k_invalid_cell_index) {
          scratch.ctu_transverse_delta[requireTouchedSlot(scratch, neighbor_target)][axis] += face_delta;
        }
      }
      const auto accumulate_stop = std::chrono::steady_clock::now();
      reconstruct_ms += std::chrono::duration<double, std::milli>(riemann_start - reconstruct_start).count();
      riemann_ms += std::chrono::duration<double, std::milli>(accumulate_start - riemann_start).count();
      accumulate_ms += std::chrono::duration<double, std::milli>(accumulate_stop - accumulate_start).count();
    }
  }

  // Final pass: the only fluxes that update conserved state and registers.
  for (std::size_t batch_begin = 0; batch_begin < active_set.active_faces.size(); batch_begin += face_batch_limit) {
    const std::size_t batch_count = std::min(
        face_batch_limit, active_set.active_faces.size() - batch_begin);
    scratch.left_states.resize(batch_count);
    scratch.right_states.resize(batch_count);
    scratch.fluxes.resize(batch_count);

    const auto reconstruct_start = std::chrono::steady_clock::now();
    reconstructBatch(batch_begin, batch_count);
    if (use_ctu) {
      for (std::size_t slot = 0; slot < batch_count; ++slot) {
        const HydroFace& face = geometry.faces[active_set.active_faces[batch_begin + slot]];
        const std::size_t axis = faceAxisIndex(face.axis);
        HydroConservedState half_transverse;
        if (halfTransverse(face.owner_cell, axis, half_transverse)) {
          if (applyTransverseCorrection(scratch.left_states[slot], half_transverse,
                  m_adiabatic_index, k_state_floor, k_state_floor)) {
            ++ctu_corrected;
          } else {
            ++ctu_inadmissible;
          }
        } else {
          ++ctu_degraded;
        }
        if (face.neighbor_cell != k_invalid_cell_index &&
            halfTransverse(face.neighbor_cell, axis, half_transverse)) {
          if (applyTransverseCorrection(scratch.right_states[slot], half_transverse,
                  m_adiabatic_index, k_state_floor, k_state_floor)) {
            ++ctu_corrected;
          } else {
            ++ctu_inadmissible;
          }
        } else {
          ++ctu_degraded;
        }
      }
    }
    const auto reconstruct_stop = std::chrono::steady_clock::now();

    const auto riemann_start = std::chrono::steady_clock::now();
    for (std::size_t slot = 0; slot < batch_count; ++slot) {
      const std::size_t face_index = active_set.active_faces[batch_begin + slot];
      scratch.fluxes[slot] = riemann_solver.computeFlux(
          scratch.left_states[slot],
          scratch.right_states[slot],
          geometry.faces[face_index],
          m_adiabatic_index);
      maybeRecordFluxRegister(
          geometry,
          face_index,
          update,
          scratch.fluxes[slot],
          flux_register_sink);
    }
    riemann_solve_count += static_cast<std::uint64_t>(batch_count);
    const auto riemann_stop = std::chrono::steady_clock::now();

    const auto accumulate_start = std::chrono::steady_clock::now();
    for (std::size_t slot = 0; slot < batch_count; ++slot) {
      const std::size_t face_index = active_set.active_faces[batch_begin + slot];
      const HydroFace& face = geometry.faces[face_index];
      const double face_scale = flux_scale * face.area_comoving;
      const HydroConservedState face_delta = face_scale * scratch.fluxes[slot];
      const std::size_t owner_slot = requireTouchedSlot(scratch, face.owner_cell);
      scratch.cell_delta[owner_slot] -= face_delta;
      const std::size_t neighbor_target = fluxNeighborTargetCell(geometry, face);
      const std::size_t neighbor_slot = neighbor_target != k_invalid_cell_index
          ? requireTouchedSlot(scratch, neighbor_target)
          : k_invalid_touched_slot;
      if (neighbor_slot != k_invalid_touched_slot) {
        scratch.cell_delta[neighbor_slot] += face_delta;
      }
      if (use_dual_energy) {
        // Upwinded advective flux of rho*e and the face normal velocity for
        // the pdV term of the auxiliary internal-energy equation.
        const double mass_flux = scratch.fluxes[slot].mass_density_comoving;
        const double upwind_e = mass_flux >= 0.0
            ? scratch.left_states[slot].specific_internal_energy_code
            : scratch.right_states[slot].specific_internal_energy_code;
        const double advected = face_scale * mass_flux * upwind_e;
        const double face_velocity = face_scale * 0.5 *
            (faceNormalVelocity(scratch.left_states[slot], face) +
             faceNormalVelocity(scratch.right_states[slot], face));
        scratch.aux_internal_energy_advection[owner_slot] -= advected;
        scratch.aux_velocity_divergence[owner_slot] += face_velocity;
        if (neighbor_slot != k_invalid_touched_slot) {
          scratch.aux_internal_energy_advection[neighbor_slot] += advected;
          scratch.aux_velocity_divergence[neighbor_slot] -= face_velocity;
        }
      }
    }
    const auto accumulate_stop = std::chrono::steady_clock::now();

    reconstruct_ms += std::chrono::duration<double, std::milli>(
        reconstruct_stop - reconstruct_start).count();
    riemann_ms += std::chrono::duration<double, std::milli>(
        riemann_stop - riemann_start).count();
    accumulate_ms += std::chrono::duration<double, std::milli>(
        accumulate_stop - accumulate_start).count();
  }
  scratch.high_water.capacity_bytes = std::max(
      scratch.high_water.capacity_bytes, scratch.ownedCapacityBytes());

  const auto source_start = std::chrono::steady_clock::now();
  std::uint64_t dual_energy_auxiliary_cells = 0U;
  std::uint64_t dual_energy_unknown_auxiliary_cells = 0U;
  std::uint64_t dual_energy_floor_cells = 0U;
  for (std::size_t touched_slot = 0; touched_slot < scratch.touched_cells.size(); ++touched_slot) {
    const std::size_t cell_index = scratch.touched_cells[touched_slot];
    const HydroConservedState old_cell = conserved.loadCell(cell_index);
    HydroPrimitiveState primitive = primitive_cache != nullptr
        ? primitive_cache->loadCell(cell_index)
        : closedPrimitiveForCell(
            conserved, cell_index, m_adiabatic_index, source_context, m_numerics.dual_energy);

    HydroConservedState source_total;
    if (std::binary_search(
            scratch.source_active_cells.begin(),
            scratch.source_active_cells.end(),
            cell_index)) {
      for (const HydroSourceTerm* source : source_terms) {
        if (source == nullptr) {
          continue;
        }
        source_total += source->sourceForCell(cell_index, old_cell, primitive, source_context);
      }
    }

    const HydroConservedState flux_delta = scratch.cell_delta[touched_slot];
    const HydroConservedState source_delta = update.dt_code * source_total;
    conservation_report.flux_delta += totalsFromDelta(
        flux_delta,
        old_cell,
        geometry.cell_volume_comoving);
    conservation_report.source_delta += totalsFromDelta(
        source_delta,
        old_cell + flux_delta,
        geometry.cell_volume_comoving);

    HydroConservedState updated = old_cell + flux_delta + source_delta;

    // Dual energy (Bryan et al. 1995 policy, no total-energy overwrite).
    // Authority: the conservative update owns rho, rho u, E and rho Z; the
    // thermal energy handed to the canonical state is either E - K or the
    // auxiliary rho e_aux evolved by the non-conservative equation
    //   (rho e)^{n+1} = (rho e)^n - dt/(aV) sum_f A F_rho e_upwind
    //                   - P^n dt/(aV) sum_f A u_face.n + dt (S_E - u.S_m),
    // with (rho e)^n taken from the auxiliary lane (initialized from the
    // canonical specific internal energy, so it never suffers the E - K
    // cancellation) and S_E - u.S_m the thermal share of the hydro sources
    // (cooling, effective-ISM relaxation; -3 H P for expansion).
    // Selection after the update: E - K if it is positive and >= eta E
    // (shocks heat through E and are captured here; the auxiliary lane is
    // resynchronized to E - K), otherwise rho e_aux if positive (cold
    // hypersonic flow, strong rarefactions). E itself is never modified here:
    // the difference (selected - (E - K)) V is a non-conservative thermal
    // correction that enters the canonical state when the runtime stores e
    // from the auxiliary lane, and it is ledgered in dual_energy_sync_delta.
    if (use_dual_energy && cell_index < real_cell_count && updated.mass_density_comoving > k_small &&
        conserved.auxiliaryInternalEnergyEnabled()) {
      std::span<double> auxiliary = conserved.auxiliaryInternalEnergyDensityComoving();
      const double old_rho = std::max(old_cell.mass_density_comoving, k_small);
      const double old_kinetic = 0.5 *
          (old_cell.momentum_density_x_comoving * old_cell.momentum_density_x_comoving +
           old_cell.momentum_density_y_comoving * old_cell.momentum_density_y_comoving +
           old_cell.momentum_density_z_comoving * old_cell.momentum_density_z_comoving) / old_rho;
      double old_internal = auxiliary[cell_index];
      if (!std::isfinite(old_internal) || !(old_internal > 0.0)) {
        old_internal = std::max(old_cell.total_energy_density_comoving - old_kinetic, 0.0);
        ++dual_energy_unknown_auxiliary_cells;
      }
      const double source_internal = source_total.total_energy_density_comoving -
          (primitive.vel_x_peculiar * source_total.momentum_density_x_comoving +
           primitive.vel_y_peculiar * source_total.momentum_density_y_comoving +
           primitive.vel_z_peculiar * source_total.momentum_density_z_comoving);
      const double auxiliary_internal = old_internal +
          scratch.aux_internal_energy_advection[touched_slot] -
          primitive.pressure_comoving * scratch.aux_velocity_divergence[touched_slot] +
          update.dt_code * source_internal;
      const double new_kinetic = 0.5 *
          (updated.momentum_density_x_comoving * updated.momentum_density_x_comoving +
           updated.momentum_density_y_comoving * updated.momentum_density_y_comoving +
           updated.momentum_density_z_comoving * updated.momentum_density_z_comoving) /
          updated.mass_density_comoving;
      const double conservative_internal = updated.total_energy_density_comoving - new_kinetic;
      const bool conservative_reliable = conservative_internal > 0.0 &&
          conservative_internal >=
              m_numerics.dual_energy.thermal_fraction_switch * updated.total_energy_density_comoving;
      double selected_internal = conservative_internal;
      if (!conservative_reliable && std::isfinite(auxiliary_internal) && auxiliary_internal > 0.0) {
        selected_internal = auxiliary_internal;
        ++dual_energy_auxiliary_cells;
      } else if (!(conservative_internal > 0.0)) {
        // Neither branch admissible: thermal floor (counted, ledgered).
        selected_internal = k_small;
        ++dual_energy_floor_cells;
      }
      auxiliary[cell_index] = selected_internal;
      if (selected_internal != conservative_internal) {
        HydroConservedState thermal_correction;
        thermal_correction.total_energy_density_comoving = selected_internal - conservative_internal;
        conservation_report.dual_energy_sync_delta += totalsFromDelta(
            thermal_correction, updated, geometry.cell_volume_comoving);
        ++conservation_report.dual_energy_sync_count;
      }
    }

    const HydroConservedState before_floor = updated;
    updated.mass_density_comoving = std::max(updated.mass_density_comoving, k_small);
    const double kinetic_density = 0.5 *
        (updated.momentum_density_x_comoving * updated.momentum_density_x_comoving +
         updated.momentum_density_y_comoving * updated.momentum_density_y_comoving +
         updated.momentum_density_z_comoving * updated.momentum_density_z_comoving) /
        std::max(updated.mass_density_comoving, k_small);
    const double minimum_total_energy_density = kinetic_density + k_small;
    const bool applied_internal_energy_floor =
        updated.total_energy_density_comoving < minimum_total_energy_density;
    if (applied_internal_energy_floor) {
      updated.total_energy_density_comoving = minimum_total_energy_density;
      ++conservation_report.internal_energy_floor_count;
    }
    if (updated.mass_density_comoving != before_floor.mass_density_comoving ||
        updated.total_energy_density_comoving != before_floor.total_energy_density_comoving) {
      conservation_report.floor_delta += totalsFromDelta(
          updated - before_floor,
          before_floor,
          geometry.cell_volume_comoving);
    }
    conserved.storeCell(cell_index, updated);
    if (primitive_cache != nullptr) {
      primitive_cache->storeCell(
          cell_index,
          closedPrimitiveForCell(
              conserved, cell_index, m_adiabatic_index, source_context, m_numerics.dual_energy));
    }
  }
  const auto source_stop = std::chrono::steady_clock::now();
  const auto total_stop = std::chrono::steady_clock::now();
  conservation_report.after = conservationTotals(conserved, geometry, scratch.full_active_cells);
  // The dual-energy thermal correction is applied to the canonical state by
  // the runtime store, not to E here, so it is not part of this residual.
  conservation_report.residual = conservation_report.after - conservation_report.before -
      conservation_report.flux_delta - conservation_report.source_delta - conservation_report.floor_delta;

  scratch.high_water.touched_cell_capacity = std::max(
      scratch.high_water.touched_cell_capacity, scratch.touched_cells.capacity());
  scratch.high_water.capacity_bytes = std::max(
      scratch.high_water.capacity_bytes, scratch.ownedCapacityBytes());

  if (profile != nullptr) {
    profile->face_count += static_cast<std::uint64_t>(active_set.active_faces.size());
    const std::uint64_t limiter_after = reconstruction.limiterClipCount();
    const std::uint64_t positivity_after = reconstruction.positivityFallbackCount();
    const std::uint64_t degradation_after = reconstruction.stencilDegradationCount();
    const std::uint64_t riemann_fallback_after = riemann_solver.fallbackCount();
    profile->limiter_clip_count += (limiter_after >= limiter_before) ? (limiter_after - limiter_before) : 0U;
    profile->positivity_fallback_count +=
        (positivity_after >= positivity_before) ? (positivity_after - positivity_before) : 0U;
    profile->stencil_degradation_count +=
        (degradation_after >= degradation_before) ? (degradation_after - degradation_before) : 0U;
    profile->riemann_fallback_count +=
        (riemann_fallback_after >= riemann_fallback_before) ? (riemann_fallback_after - riemann_fallback_before) : 0U;
    profile->riemann_solve_count += riemann_solve_count;
    profile->ctu_transverse_corrected_states += ctu_corrected;
    profile->ctu_transverse_degraded_states += ctu_degraded;
    profile->ctu_transverse_inadmissible_states += ctu_inadmissible;
    profile->dual_energy_auxiliary_cells += dual_energy_auxiliary_cells;
    profile->dual_energy_unknown_auxiliary_cells += dual_energy_unknown_auxiliary_cells;
    profile->dual_energy_floor_cells += dual_energy_floor_cells;
    profile->reconstruct_ms += reconstruct_ms;
    profile->riemann_ms += riemann_ms;
    profile->accumulate_ms += accumulate_ms;
    profile->source_ms += std::chrono::duration<double, std::milli>(source_stop - source_start).count();
    profile->total_ms += std::chrono::duration<double, std::milli>(total_stop - total_start).count();
    const std::uint64_t touched_count = static_cast<std::uint64_t>(scratch.touched_cells.size());
    // Face traffic: two reconstructed primitive states (9 lanes each) plus
    // one flux (6 lanes) per Riemann solve; cell traffic: load+store of the
    // conserved state (6 lanes each way is counted once, historical rule).
    const std::uint64_t solves = riemann_solve_count;
    if (touched_count > std::numeric_limits<std::uint64_t>::max() / 6U ||
        solves > std::numeric_limits<std::uint64_t>::max() / 10U) {
      throw std::overflow_error("hydro profile moved-element count overflow");
    }
    const std::uint64_t touched_terms = touched_count * 6U;
    const std::uint64_t face_terms = solves * 10U;
    if (touched_terms > std::numeric_limits<std::uint64_t>::max() - face_terms) {
      throw std::overflow_error("hydro profile moved-element aggregate overflow");
    }
    const std::uint64_t moved_values = touched_terms + face_terms;
    if (moved_values > std::numeric_limits<std::uint64_t>::max() / sizeof(double)) {
      throw std::overflow_error("hydro profile byte count overflow");
    }
    const std::uint64_t moved_bytes = moved_values * sizeof(double);
    if (profile->bytes_moved > std::numeric_limits<std::uint64_t>::max() - moved_bytes) {
      throw std::overflow_error("hydro profile cumulative byte count overflow");
    }
    profile->bytes_moved += moved_bytes;
    profile->conservation = conservation_report;
  }
}

}  // namespace cosmosim::hydro
