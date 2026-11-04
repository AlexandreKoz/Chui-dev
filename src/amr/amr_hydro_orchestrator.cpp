#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/gas_thermodynamic_refresh.hpp"

#include "cosmosim/amr/amr_patch_indexing.hpp"
#include "amr/internal/amr_retained_state_admission.hpp"
#include "amr/internal/amr_periodic_topology.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/governed_scratch_arena.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <memory_resource>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace cosmosim::amr {

std::uint64_t amrPreparedGhostWorkspaceBytes(
    std::size_t patch_count, std::size_t prepared_ghost_count) {
  const std::size_t offset_count = core::checkedSizeAdd(
      patch_count, 1U, "AMR prepared ghost offsets");
  return core::checkedMemoryBytesAdd(
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          prepared_ghost_count, sizeof(hydro::HydroConservedState),
          "AMR prepared ghost states")),
      core::checkedMemoryBytesAdd(
          static_cast<std::uint64_t>(core::checkedSizeMultiply(
              offset_count, sizeof(std::size_t), "AMR prepared ghost offsets")),
          256U, "AMR prepared ghost alignment"),
      "AMR prepared ghost coexistence");
}

namespace {

constexpr double k_geometry_tol = 1.0e-10;

// Ghost layers needed for full-order reconstruction across patch interfaces:
// the reconstruction stencil half-width (PCM 1, PLM 2, PPM 3). An explicit
// non-zero ProductionAmrHydroOptions::ghost_layers overrides (tests).
[[nodiscard]] std::uint8_t amrHydroGhostLayersForOptions(const ProductionAmrHydroOptions& options) {
  if (options.ghost_layers != 0U) {
    return options.ghost_layers;
  }
  switch (options.reconstruction_scheme) {
    case hydro::HydroReconstructionScheme::kPiecewiseConstant: return 1U;
    case hydro::HydroReconstructionScheme::kPiecewiseLinear: return 2U;
    case hydro::HydroReconstructionScheme::kPiecewiseParabolic: return 3U;
  }
  return 1U;
}

[[nodiscard]] std::size_t product(std::array<std::uint16_t, 3> dims) {
  return core::checkedSizeProduct3(
      static_cast<std::size_t>(dims[0]),
      static_cast<std::size_t>(dims[1]),
      static_cast<std::size_t>(dims[2]),
      "AMR patch cell-count product");
}

[[nodiscard]] std::size_t linearIndex(
    std::array<std::uint16_t, 3> dims,
    std::size_t i,
    std::size_t j,
    std::size_t k) {
  const std::size_t yz = core::checkedSizeAdd(
      j,
      core::checkedSizeMultiply(static_cast<std::size_t>(dims[1]), k, "AMR linear index yz"),
      "AMR linear index yz sum");
  return core::checkedSizeAdd(
      i,
      core::checkedSizeMultiply(static_cast<std::size_t>(dims[0]), yz, "AMR linear index x"),
      "AMR linear index total");
}

[[nodiscard]] std::vector<double> sortedUnique(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  values.erase(
      std::unique(
          values.begin(),
          values.end(),
          [](double lhs, double rhs) {
            const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
            return std::abs(lhs - rhs) <= k_geometry_tol * scale;
          }),
      values.end());
  return values;
}

[[nodiscard]] double spacingForUnique(const std::vector<double>& values, double fallback) {
  if (values.size() <= 1U) {
    return std::max(1.0e-12, fallback);
  }
  return std::max(1.0e-12, (values.back() - values.front()) / static_cast<double>(values.size() - 1U));
}

[[nodiscard]] std::array<std::uint16_t, 3> nearCubicFactors(std::size_t count) {
  if (count == 0U || count > std::numeric_limits<std::uint16_t>::max()) {
    throw std::invalid_argument("AMR production patch descriptor requires positive representable cell count");
  }
  std::array<std::uint16_t, 3> best{1, 1, static_cast<std::uint16_t>(count)};
  std::size_t best_surface = std::numeric_limits<std::size_t>::max();
  for (std::size_t nx = 1; nx <= count; ++nx) {
    if (count % nx != 0U) {
      continue;
    }
    const std::size_t yz = count / nx;
    for (std::size_t ny = 1; ny <= yz; ++ny) {
      if (yz % ny != 0U) {
        continue;
      }
      const std::size_t nz = yz / ny;
      const std::size_t surface = nx * ny + nx * nz + ny * nz;
      if (surface < best_surface && nx <= std::numeric_limits<std::uint16_t>::max() &&
          ny <= std::numeric_limits<std::uint16_t>::max() && nz <= std::numeric_limits<std::uint16_t>::max()) {
        best = {static_cast<std::uint16_t>(nx), static_cast<std::uint16_t>(ny), static_cast<std::uint16_t>(nz)};
        best_surface = surface;
      }
    }
  }
  return best;
}

[[nodiscard]] std::array<double, 3> cellWidths(const PatchDescriptor& patch) {
  return {
      patch.extent_comov[0] / static_cast<double>(patch.cell_dims[0]),
      patch.extent_comov[1] / static_cast<double>(patch.cell_dims[1]),
      patch.extent_comov[2] / static_cast<double>(patch.cell_dims[2]),
  };
}

[[nodiscard]] std::array<double, 3> cellCenter(const PatchDescriptor& patch, std::size_t cell) {
  const auto widths = cellWidths(patch);
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t i = cell % nx;
  const std::size_t j = (cell / nx) % ny;
  const std::size_t k = cell / (nx * ny);
  return {
      patch.origin_comov[0] + (static_cast<double>(i) + 0.5) * widths[0],
      patch.origin_comov[1] + (static_cast<double>(j) + 0.5) * widths[1],
      patch.origin_comov[2] + (static_cast<double>(k) + 0.5) * widths[2],
  };
}

[[nodiscard]] bool containsPoint(const PatchDescriptor& patch, const std::array<double, 3>& point) {
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const double lower = patch.origin_comov[axis] - k_geometry_tol;
    const double upper = patch.origin_comov[axis] + patch.extent_comov[axis] + k_geometry_tol;
    if (point[axis] < lower || point[axis] >= upper) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] std::size_t cellContainingPoint(
    const PatchDescriptor& patch,
    const std::array<double, 3>& point) {
  const auto widths = cellWidths(patch);
  std::array<std::size_t, 3> ijk{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    double coordinate = (point[axis] - patch.origin_comov[axis]) / widths[axis];
    coordinate = std::clamp(coordinate, 0.0, static_cast<double>(patch.cell_dims[axis]) - k_geometry_tol);
    ijk[axis] = static_cast<std::size_t>(coordinate);
    if (ijk[axis] >= patch.cell_dims[axis]) {
      ijk[axis] = patch.cell_dims[axis] - 1U;
    }
  }
  return linearIndex(patch.cell_dims, ijk[0], ijk[1], ijk[2]);
}

[[nodiscard]] int axisIndex(hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return 0;
    case hydro::HydroFaceAxis::kY:
      return 1;
    case hydro::HydroFaceAxis::kZ:
      return 2;
  }
  return 0;
}

[[nodiscard]] double sideSign(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? -1.0 : 1.0;
}

[[nodiscard]] hydro::HydroFaceSide oppositeSide(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? hydro::HydroFaceSide::kUpper : hydro::HydroFaceSide::kLower;
}

[[nodiscard]] std::array<double, 3> ghostProbePoint(
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroGhostCell& ghost,
    std::uint8_t layer = 1U,
    const internal::AmrPeriodicFrame& frame = {}) {
  std::array<double, 3> point = cellCenter(patch_geometry.patch, ghost.owner_real_cell);
  const auto widths = cellWidths(patch_geometry.patch);
  point[axisIndex(ghost.axis)] += sideSign(ghost.side) *
      static_cast<double>(std::max<std::uint8_t>(layer, 1U)) * widths[axisIndex(ghost.axis)];
  return internal::wrapPoint(frame, point);
}

[[nodiscard]] std::optional<PatchDescriptor> sourcePatchForGhost(
    const AmrHydroPatchGeometry& target,
    const hydro::HydroGhostCell& ghost,
    std::span<const PatchDescriptor> all_patches,
    const internal::AmrPeriodicFrame& frame = {},
    std::uint8_t layer = 1U) {
  const std::array<double, 3> probe = ghostProbePoint(target, ghost, layer, frame);
  for (const PatchDescriptor& candidate : all_patches) {
    if (candidate.patch_id == target.patch.patch_id) {
      continue;
    }
    if (containsPoint(candidate, probe)) {
      return candidate;
    }
  }
  return std::nullopt;
}

[[nodiscard]] bool intervalsOverlap(double a0, double a1, double b0, double b1) {
  return std::max(a0, b0) < std::min(a1, b1) + k_geometry_tol;
}

[[nodiscard]] bool touchesOnSide(
    const PatchDescriptor& target,
    const PatchDescriptor& candidate,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    const internal::AmrPeriodicFrame& frame = {}) {
  const int normal_axis = axisIndex(axis);
  const double target_face = side == hydro::HydroFaceSide::kLower
      ? target.origin_comov[normal_axis]
      : target.origin_comov[normal_axis] + target.extent_comov[normal_axis];
  const double candidate_face = side == hydro::HydroFaceSide::kLower
      ? candidate.origin_comov[normal_axis] + candidate.extent_comov[normal_axis]
      : candidate.origin_comov[normal_axis];
  const double scale = std::max({1.0, std::abs(target_face), std::abs(candidate_face)});
  // Periodic axes compare face positions with the minimum image, so the
  // lower box face (0) touches an upper face at L. Non-periodic axes keep the
  // exact historical comparison.
  const double face_gap = internal::minImageDifference(
      frame, static_cast<std::size_t>(normal_axis), target_face, candidate_face);
  if (std::abs(face_gap) > k_geometry_tol * scale) {
    return false;
  }
  for (int other_axis = 0; other_axis < 3; ++other_axis) {
    if (other_axis == normal_axis) {
      continue;
    }
    const double a0 = target.origin_comov[other_axis];
    const double a1 = a0 + target.extent_comov[other_axis];
    const double b0 = candidate.origin_comov[other_axis];
    const double b1 = b0 + candidate.extent_comov[other_axis];
    if (!intervalsOverlap(a0, a1, b0, b1)) {
      return false;
    }
  }
  return true;
}

// Periodic AMR topology (CHUI Part 3 Checkpoint A; box convention [0, L)
// per axis, the same as core::canonicalPeriodicPosition and the TreePM
// domain). A patch that spans the box along an axis wraps in-patch (the
// geometry builder sources the ghost from the opposite face of the SAME
// patch). Any other patch face lying on a periodic box face must touch at
// least one wrapped neighbour patch (same level or coarse/fine, local or
// remote) through the minimum-image face test: that face is then an ordinary
// same-level / coarse-fine interface with wrapped ghost sampling, wrapped
// flux-register keys and the ordinary distributed exchange. Only a hierarchy
// that does not tile the periodic box (no wrapped neighbour at all) is still
// refused, because closing the face or inventing a ghost would be wrong.
void validatePeriodicAmrPatchTopology(
    const PatchDescriptor& patch,
    std::span<const PatchDescriptor> patches,
    const std::array<double, 3>& box_extent,
    const std::array<AmrHydroBoundaryClass, 6>& boundary_classes) {
  (void)patches;
  for (int axis = 0; axis < 3; ++axis) {
    const double length = box_extent[static_cast<std::size_t>(axis)];
    if (!(length > 0.0) || !std::isfinite(length)) {
      throw std::invalid_argument(
          "periodic AMR hydro requires ProductionAmrHydroOptions::periodic_box_extent_comoving > 0");
    }
    const double lower = patch.origin_comov[static_cast<std::size_t>(axis)];
    const double upper = lower + patch.extent_comov[static_cast<std::size_t>(axis)];
    const double tolerance = k_geometry_tol * std::max(1.0, length);
    if (lower < -tolerance || upper > length + tolerance) {
      throw std::runtime_error(
          "periodic AMR hydro: patch " + std::to_string(patch.patch_id) +
          " lies outside the canonical periodic box [0, L) on axis " + std::to_string(axis));
    }
    const bool spans = std::abs(patch.extent_comov[static_cast<std::size_t>(axis)] - length) <= tolerance;
    if (spans) {
      continue;
    }
    // A non-spanning patch on a periodic axis has no physical face at all:
    // both sides (interior or across the wrap) must touch a neighbour patch.
    const std::size_t lower_slot = static_cast<std::size_t>(2 * axis);
    if (boundary_classes[lower_slot] == AmrHydroBoundaryClass::kPhysical ||
        boundary_classes[lower_slot + 1U] == AmrHydroBoundaryClass::kPhysical) {
      throw std::runtime_error(
          "periodic AMR hydro: patch " + std::to_string(patch.patch_id) +
          " does not span periodic axis " + std::to_string(axis) +
          " and a face has no (local or remote, direct or wrapped) neighbour patch; the AMR "
          "hierarchy does not tile the periodic box (refusing instead of closing the face)");
    }
  }
}

[[nodiscard]] AmrHydroBoundaryClass boundaryClassForSide(
    const PatchDescriptor& target,
    std::span<const PatchDescriptor> patches,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    const internal::AmrPeriodicFrame& frame = {}) {
  for (const PatchDescriptor& candidate : patches) {
    if (candidate.patch_id == target.patch_id) {
      continue;
    }
    if (!touchesOnSide(target, candidate, axis, side, frame)) {
      continue;
    }
    return candidate.level == target.level
        ? AmrHydroBoundaryClass::kSameLevel
        : AmrHydroBoundaryClass::kCoarseFine;
  }
  return AmrHydroBoundaryClass::kPhysical;
}

[[nodiscard]] std::uint64_t makeRegisterKey(
    std::uint64_t coarse_patch_id,
    std::size_t coarse_cell_index,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide orientation) {
  std::uint64_t key = 1469598103934665603ULL;
  const auto mix = [&key](std::uint64_t value) {
    key ^= value;
    key *= 1099511628211ULL;
  };
  mix(coarse_patch_id);
  mix(static_cast<std::uint64_t>(coarse_cell_index) + 1U);
  mix(static_cast<std::uint64_t>(axisIndex(axis)) + 11U);
  mix(orientation == hydro::HydroFaceSide::kLower ? 23U : 29U);
  return key == 0U ? 1U : key;
}

[[nodiscard]] double safeDensity(const core::SimulationState& state, std::uint32_t row, double floor) {
  return std::max(state.gas_cells.density_code[row], floor);
}

[[nodiscard]] hydro::HydroPrimitiveState primitiveForRow(
    const core::SimulationState& state,
    std::uint32_t row,
    double density_floor,
    double pressure_floor,
    double adiabatic_index) {
  const double rho = safeDensity(state, row, density_floor);
  const double internal_energy_floor = pressure_floor > 0.0
      ? pressure_floor / ((adiabatic_index - 1.0) * std::max(rho, 1.0e-30))
      : 0.0;
  // Canonical thermal authority: specific internal energy. The pressure
  // cache is used only when e was never set (legacy/test states); taking
  // max(P_cache, (gamma-1) rho e) silently injected energy whenever the cache
  // was stale or held an effective-EOS pressure.
  const double canonical_e = state.gas_cells.internal_energy_code[row];
  const double internal_energy = std::max(canonical_e, internal_energy_floor);
  const double pressure = canonical_e > 0.0
      ? std::max((adiabatic_index - 1.0) * rho * internal_energy, pressure_floor)
      : std::max(state.gas_cells.pressure_code[row], pressure_floor);
  return hydro::HydroPrimitiveState{
      .rho_comoving = rho,
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = pressure,
      .metallicity_mass_fraction = std::clamp(
          state.gas_cells.metal_mass_code[row] /
              std::max(state.cells.mass_code[row], 1.0e-30),
          0.0,
          1.0)};
}

[[nodiscard]] ConservedState volumeIntegratedForRow(
    const core::SimulationState& state,
    std::uint32_t row,
    double volume,
    const ProductionAmrHydroOptions& options) {
  const hydro::HydroConservedState density = hydro::HydroCoreSolver::conservedFromPrimitive(
      primitiveForRow(state, row, options.density_floor, options.pressure_floor, options.adiabatic_index),
      options.adiabatic_index);
  return ConservedState{
      .mass_code = density.mass_density_comoving * volume,
      .momentum_x_code = density.momentum_density_x_comoving * volume,
      .momentum_y_code = density.momentum_density_y_comoving * volume,
      .momentum_z_code = density.momentum_density_z_comoving * volume,
      .total_energy_code = density.total_energy_density_comoving * volume,
      .metal_mass_code = density.metal_mass_density_comoving * volume};
}


void refreshTraceableGhostSourceIds(
    std::vector<AmrHydroPatchGeometry>& geometries,
    const internal::AmrPeriodicFrame& frame = {}) {
  std::unordered_map<std::uint64_t, std::size_t> patch_slot_by_id;
  patch_slot_by_id.reserve(geometries.size());
  for (std::size_t slot = 0; slot < geometries.size(); ++slot) {
    patch_slot_by_id.emplace(geometries[slot].patch.patch_id, slot);
  }
  for (AmrHydroPatchGeometry& target : geometries) {
    for (AmrHydroGhostDescriptor& ghost_descriptor : target.ghosts) {
      if (ghost_descriptor.source_patch_id == 0U) {
        continue;
      }
      const auto source_slot_it = patch_slot_by_id.find(ghost_descriptor.source_patch_id);
      if (source_slot_it == patch_slot_by_id.end()) {
        continue;
      }
      const AmrHydroPatchGeometry& source = geometries[source_slot_it->second];
      const hydro::HydroGhostCell& ghost = target.geometry.ghost_cells.at(ghost_descriptor.ghost_slot);
      const std::array<double, 3> probe =
          ghostProbePoint(target, ghost, ghost_descriptor.ghost_layer, frame);
      if (!containsPoint(source.patch, probe)) {
        continue;
      }
      const std::size_t source_cell = cellContainingPoint(source.patch, probe);
      if (source_cell < source.gas_cell_ids.size()) {
        ghost_descriptor.source_gas_cell_id = source.gas_cell_ids[source_cell];
        for (std::size_t face_index = 0; face_index < target.geometry.faces.size(); ++face_index) {
          const hydro::HydroFace& face = target.geometry.faces[face_index];
          if (face.ghost_cell_slot == ghost_descriptor.ghost_slot &&
              face_index < target.geometry.flux_register_faces.size() &&
              target.geometry.flux_register_faces[face_index].role == hydro::HydroFluxRegisterFaceRole::kFine &&
              target.geometry.flux_register_faces[face_index].coarse_gas_cell_id == 0U) {
            target.geometry.flux_register_faces[face_index].coarse_gas_cell_id = ghost_descriptor.source_gas_cell_id;
          }
        }
      }
    }
  }
}

[[nodiscard]] const AmrHydroSparseRemoteCell* findRemoteBoundaryCell(
    const DistributedAmrRemotePatch& remote,
    std::size_t patch_local_cell) {
  const auto it = std::lower_bound(
      remote.boundary_cells.begin(),
      remote.boundary_cells.end(),
      patch_local_cell,
      [](const AmrHydroSparseRemoteCell& cell, std::size_t offset) {
        return static_cast<std::size_t>(cell.patch_local_cell) < offset;
      });
  if (it == remote.boundary_cells.end() ||
      static_cast<std::size_t>(it->patch_local_cell) != patch_local_cell) {
    return nullptr;
  }
  return &*it;
}

void refreshRemoteTraceableGhostSourceIds(
    std::vector<AmrHydroPatchGeometry>& geometries,
    std::span<const DistributedAmrRemotePatch> remote_patches,
    const internal::AmrPeriodicFrame& frame = {}) {
  std::unordered_map<std::uint64_t, const DistributedAmrRemotePatch*> remote_by_id;
  remote_by_id.reserve(remote_patches.size());
  for (const DistributedAmrRemotePatch& remote : remote_patches) {
    remote_by_id.emplace(remote.patch.patch_id, &remote);
  }
  for (AmrHydroPatchGeometry& target : geometries) {
    for (AmrHydroGhostDescriptor& ghost_descriptor : target.ghosts) {
      const auto remote_it = remote_by_id.find(ghost_descriptor.source_patch_id);
      if (remote_it == remote_by_id.end()) {
        continue;
      }
      const DistributedAmrRemotePatch& remote = *remote_it->second;
      const hydro::HydroGhostCell& ghost =
          target.geometry.ghost_cells.at(ghost_descriptor.ghost_slot);
      const std::array<double, 3> probe =
          ghostProbePoint(target, ghost, ghost_descriptor.ghost_layer, frame);
      if (!containsPoint(remote.patch, probe)) {
        continue;
      }
      const std::size_t source_cell = cellContainingPoint(remote.patch, probe);
      const AmrHydroSparseRemoteCell* remote_cell =
          findRemoteBoundaryCell(remote, source_cell);
      if (remote_cell == nullptr) {
        continue;
      }
      ghost_descriptor.source_gas_cell_id = remote_cell->gas_cell_id;
      for (std::size_t face_index = 0;
           face_index < target.geometry.faces.size(); ++face_index) {
        const hydro::HydroFace& face = target.geometry.faces[face_index];
        if (face.ghost_cell_slot == ghost_descriptor.ghost_slot &&
            face_index < target.geometry.flux_register_faces.size() &&
            target.geometry.flux_register_faces[face_index].role ==
                hydro::HydroFluxRegisterFaceRole::kFine &&
            target.geometry.flux_register_faces[face_index].coarse_gas_cell_id == 0U) {
          target.geometry.flux_register_faces[face_index].coarse_gas_cell_id =
              remote_cell->gas_cell_id;
        }
      }
    }
  }
}

// Rewrite P, c_s and T of one gas row after an AMR mutation, according to
// the cache policy (see AmrThermodynamicCachePolicy).
void storeDerivedGasThermodynamics(
    core::SimulationState& state,
    std::uint32_t row,
    double adiabatic_index,
    const AmrThermodynamicCachePolicy& policy) {
  if (policy.temperature_k_per_internal_energy_code <= 0.0 && policy.closure == nullptr) {
    state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];  // legacy
    return;
  }
  if (policy.temperature_k_per_internal_energy_code <= 0.0) {
    state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];  // legacy T
  }
  (void)hydro::refreshGasThermodynamicCell(
      state, row,
      hydro::GasThermodynamicRefreshPolicy{
          .adiabatic_index = adiabatic_index,
          .density_floor_code = 1.0e-14,
          .pressure_floor_code = 1.0e-14,
          .scale_factor = policy.scale_factor,
          .redshift = policy.redshift,
          .closure = policy.closure,
          .temperature_k_per_internal_energy_code = policy.temperature_k_per_internal_energy_code});
}

void writeVolumeIntegratedToRow(
    core::SimulationState& state,
    std::uint32_t row,
    const ConservedState& volume_state,
    double volume,
    const ProductionAmrHydroOptions& options) {
  if (volume <= 0.0) {
    throw std::invalid_argument("writeVolumeIntegratedToRow requires positive volume");
  }
  hydro::HydroConservedState density_state{
      .mass_density_comoving = volume_state.mass_code / volume,
      .momentum_density_x_comoving = volume_state.momentum_x_code / volume,
      .momentum_density_y_comoving = volume_state.momentum_y_code / volume,
      .momentum_density_z_comoving = volume_state.momentum_z_code / volume,
      .total_energy_density_comoving = volume_state.total_energy_code / volume,
      .metal_mass_density_comoving = volume_state.metal_mass_code / volume};
  hydro::HydroPrimitiveState primitive =
      hydro::HydroCoreSolver::primitiveFromConserved(density_state, options.adiabatic_index);
  primitive.rho_comoving = std::max(primitive.rho_comoving, options.density_floor);
  primitive.pressure_comoving = std::max(primitive.pressure_comoving, options.pressure_floor);
  state.cells.mass_code[row] = primitive.rho_comoving * volume;
  state.gas_cells.density_code[row] = primitive.rho_comoving;
  state.gas_cells.pressure_code[row] = primitive.pressure_comoving;
  state.gas_cells.internal_energy_code[row] =
      primitive.pressure_comoving / ((options.adiabatic_index - 1.0) * primitive.rho_comoving);
  state.gas_cells.velocity_x_peculiar[row] = primitive.vel_x_peculiar;
  state.gas_cells.velocity_y_peculiar[row] = primitive.vel_y_peculiar;
  state.gas_cells.velocity_z_peculiar[row] = primitive.vel_z_peculiar;
  state.gas_cells.metal_mass_code[row] =
      std::clamp(primitive.metallicity_mass_fraction, 0.0, 1.0) * state.cells.mass_code[row];
  state.gas_cells.sound_speed_code[row] =
      std::sqrt(std::max(0.0, options.adiabatic_index * primitive.pressure_comoving / primitive.rho_comoving));
  storeDerivedGasThermodynamics(state, row, options.adiabatic_index, options.thermodynamic_cache);
}

[[nodiscard]] std::vector<std::uint32_t> patchLocalRowsForPatch(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::string_view context) {
  return rowsForPatchLocalCellOrder(state, patch, context);
}

[[nodiscard]] std::array<double, 5> totalsForRows(
    const core::SimulationState& state,
    std::span<const std::uint32_t> rows,
    double volume,
    const ProductionAmrHydroOptions& options) {
  std::array<double, 5> totals{};
  for (const std::uint32_t row : rows) {
    const ConservedState conserved = volumeIntegratedForRow(state, row, volume, options);
    totals[0] += conserved.mass_code;
    totals[1] += conserved.momentum_x_code;
    totals[2] += conserved.momentum_y_code;
    totals[3] += conserved.momentum_z_code;
    totals[4] += conserved.total_energy_code;
  }
  return totals;
}

void copyCellRow(core::SimulationState& state, std::uint32_t dst, const core::SimulationState& old, std::uint32_t src) {
  state.cells.center_x_comoving[dst] = old.cells.center_x_comoving[src];
  state.cells.center_y_comoving[dst] = old.cells.center_y_comoving[src];
  state.cells.center_z_comoving[dst] = old.cells.center_z_comoving[src];
  state.cells.mass_code[dst] = old.cells.mass_code[src];
  state.cells.time_bin[dst] = old.cells.time_bin[src];
  state.gas_cells.gas_cell_id[dst] = old.gas_cells.gas_cell_id[src];
  state.gas_cells.parent_particle_id[dst] = old.gas_cells.parent_particle_id[src];
  state.gas_cells.velocity_x_peculiar[dst] = old.gas_cells.velocity_x_peculiar[src];
  state.gas_cells.velocity_y_peculiar[dst] = old.gas_cells.velocity_y_peculiar[src];
  state.gas_cells.velocity_z_peculiar[dst] = old.gas_cells.velocity_z_peculiar[src];
  state.gas_cells.density_code[dst] = old.gas_cells.density_code[src];
  state.gas_cells.pressure_code[dst] = old.gas_cells.pressure_code[src];
  state.gas_cells.internal_energy_code[dst] = old.gas_cells.internal_energy_code[src];
  state.gas_cells.metal_mass_code[dst] = old.gas_cells.metal_mass_code[src];
  state.gas_cells.temperature_code[dst] = old.gas_cells.temperature_code[src];
  state.gas_cells.sound_speed_code[dst] = old.gas_cells.sound_speed_code[src];
}

void rebuildIdentityFromSidecars(
    core::SimulationState& state,
    std::uint64_t identity_generation) {
  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    if (state.gas_cells.gas_cell_id[row] == 0U) {
      throw std::runtime_error("production AMR hydro generated a zero gas_cell_id");
    }
    const std::uint32_t patch_index = state.cells.patch_index[row];
    if (patch_index >= state.patches.size()) {
      throw std::runtime_error("production AMR hydro generated an invalid patch index");
    }
    const std::uint64_t parent_id = state.gas_cells.parent_particle_id[row];
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = state.gas_cells.gas_cell_id[row],
        .parent_particle_id = parent_id == 0U ? std::optional<std::uint64_t>{} : std::optional<std::uint64_t>{parent_id},
        .owning_patch_id = state.patches.patch_id[patch_index],
        .local_cell_row = core::checkedIntegralNarrow<std::uint32_t>(
            row, "production AMR identity local row")});
  }
  state.restoreGasCellIdentityRecords(std::move(records), identity_generation);
}

[[nodiscard]] bool idRangeOverflows(std::uint64_t first, std::size_t count) {
  return count == 0U || first == 0U || first > std::numeric_limits<std::uint64_t>::max() - static_cast<std::uint64_t>(count - 1U);
}

void validatePatchIdRangeAvailable(
    const core::SimulationState& state,
    std::uint64_t first_patch_id,
    std::size_t count,
    std::uint64_t allowed_existing_patch_id = 0U) {
  if (idRangeOverflows(first_patch_id, count)) {
    throw std::invalid_argument("production AMR regrid patch-id range is zero or overflows");
  }
  // Gate-3: the requested IDs are a single contiguous interval
  // [first, first+count). No O(N) hash set is needed: an existing ID collides
  // iff it lies in that half-open interval. Deterministic, allocation-free.
  const std::uint64_t last_exclusive = first_patch_id + static_cast<std::uint64_t>(count);
  for (const std::uint64_t patch_id : state.patches.patch_id) {
    if (patch_id != allowed_existing_patch_id && patch_id >= first_patch_id &&
        patch_id < last_exclusive) {
      throw std::invalid_argument("production AMR regrid patch-id range collides with existing state");
    }
  }
}

void validateGasCellIdRangeAvailable(
    const core::SimulationState& state,
    std::uint64_t first_gas_cell_id,
    std::size_t count) {
  if (idRangeOverflows(first_gas_cell_id, count)) {
    throw std::invalid_argument("production AMR regrid gas-cell-id range is zero or overflows");
  }
  // Gate-3: same contiguous-interval check as patch IDs. Eliminates the
  // population-scale unordered_set without changing collision semantics.
  const std::uint64_t last_exclusive = first_gas_cell_id + static_cast<std::uint64_t>(count);
  for (const std::uint64_t gas_cell_id : state.gas_cells.gas_cell_id) {
    if (gas_cell_id >= first_gas_cell_id && gas_cell_id < last_exclusive) {
      throw std::invalid_argument("production AMR regrid gas-cell-id range collides with existing state");
    }
  }
  for (const core::GasCellIdentityRecord& record : state.gas_cell_identity.records()) {
    if (record.gas_cell_id >= first_gas_cell_id && record.gas_cell_id < last_exclusive) {
      throw std::invalid_argument("production AMR regrid gas-cell-id range collides with identity map");
    }
  }
}

[[nodiscard]] std::uint64_t nextFreePatchId(const core::SimulationState& state, std::size_t count) {
  std::uint64_t candidate = 1U;
  for (const std::uint64_t patch_id : state.patches.patch_id) {
    if (patch_id == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("production AMR next patch-id seed overflows uint64");
    }
    candidate = std::max(candidate, patch_id + 1U);
  }
  validatePatchIdRangeAvailable(state, candidate, count);
  return candidate;
}

[[nodiscard]] std::uint64_t nextFreeGasCellId(const core::SimulationState& state, std::size_t count) {
  std::uint64_t candidate = 1U;
  for (const std::uint64_t gas_cell_id : state.gas_cells.gas_cell_id) {
    if (gas_cell_id == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("production AMR next gas-cell-id seed overflows uint64");
    }
    candidate = std::max(candidate, gas_cell_id + 1U);
  }
  validateGasCellIdRangeAvailable(state, candidate, count);
  return candidate;
}

[[nodiscard]] std::optional<std::uint32_t> patchIndexById(
    const core::SimulationState& state,
    std::uint64_t patch_id) {
  for (std::size_t index = 0; index < state.patches.size(); ++index) {
    if (state.patches.patch_id[index] == patch_id) {
      return core::checkedIntegralNarrow<std::uint32_t>(
          index, "production AMR patch index");
    }
  }
  return std::nullopt;
}

void requireProductionRegridSynchronizationQuiescent(
    const core::SimulationState& state,
    std::string_view operation) {
  if (!state.amr_temporal_boundary_history.empty()) {
    throw std::runtime_error(
        std::string(operation) +
        " is prohibited while restart-resumable temporal boundary history is active");
  }
  // PendingFluxRegisterRecord identifies the coarse target but does not carry
  // every fine-patch contributor identity. Therefore intersection-local
  // invalidation cannot be proven for an arbitrary regrid. Fail closed until
  // all pending reflux authority has drained.
  if (!state.pending_flux_registers.empty()) {
    throw std::runtime_error(
        std::string(operation) +
        " is prohibited while pending AMR flux-register synchronization is active");
  }
}

[[nodiscard]] std::uint64_t nextIdentityGeneration(
    const core::SimulationState& state) {
  const std::uint64_t generation = state.gasCellIdentityGeneration();
  if (generation == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("production AMR gas-cell identity generation overflows uint64");
  }
  return generation + 1U;
}

[[nodiscard]] std::uint64_t productionAmrCandidateReservationBytes(
    std::size_t cell_count,
    std::size_t patch_count,
    std::size_t extra_scratch_bytes) {
  constexpr std::size_t k_cell_soa_bytes_per_cell =
      4U * sizeof(double) + sizeof(std::uint8_t) + sizeof(std::uint32_t);
  constexpr std::size_t k_gas_sidecar_bytes_per_cell =
      2U * sizeof(std::uint64_t) + 9U * sizeof(double);
  constexpr std::size_t k_patch_soa_bytes_per_patch =
      3U * sizeof(std::uint64_t) + sizeof(std::int32_t) +
      3U * sizeof(std::uint32_t) + 6U * sizeof(double) +
      3U * sizeof(std::uint16_t);
  using GasIdLookupValue = std::pair<const std::uint64_t, std::size_t>;
  using LocalRowLookupValue = std::pair<const std::uint32_t, std::size_t>;
  // GasCellIdentityMap::assign builds two unordered maps. Reserve a
  // conservative per-record allowance for both map values, node links and
  // bucket pointers so the governor sees the hidden lookup allocation too.
  constexpr std::size_t k_identity_lookup_bytes_per_cell =
      sizeof(GasIdLookupValue) + sizeof(LocalRowLookupValue) +
      8U * sizeof(void*);
  const std::size_t per_cell = core::checkedSizeAdd(
      core::checkedSizeAdd(
          k_cell_soa_bytes_per_cell, k_gas_sidecar_bytes_per_cell,
          "production AMR candidate per-cell state bytes"),
      core::checkedSizeAdd(
          sizeof(core::GasCellIdentityRecord), k_identity_lookup_bytes_per_cell,
          "production AMR candidate identity bytes"),
      "production AMR candidate per-cell total bytes");
  std::size_t bytes = core::checkedSizeMultiply(
      cell_count, per_cell, "production AMR candidate cell bytes");
  bytes = core::checkedSizeAdd(
      bytes,
      core::checkedSizeMultiply(
          patch_count, k_patch_soa_bytes_per_patch,
          "production AMR candidate patch bytes"),
      "production AMR candidate state bytes");
  bytes = core::checkedSizeAdd(
      bytes, extra_scratch_bytes, "production AMR candidate transaction bytes");
  return core::checkedIntegralNarrow<std::uint64_t>(
      bytes, "production AMR candidate reservation width");
}

struct ProductionAmrRegridAdmission {
  core::MemoryReservation reservation;
  std::uint64_t process_baseline_before = 0U;
  std::uint64_t state_baseline_before = 0U;
};

[[nodiscard]] ProductionAmrRegridAdmission admitProductionAmrRegrid(
    const core::SimulationState& state,
    const ProductionAmrHydroOptions& options,
    std::uint64_t requested_bytes,
    std::string_view owner) {
  if (options.regrid_memory_governor == nullptr) {
    throw std::logic_error(
        "production AMR regrid requires the process MemoryGovernor in ProductionAmrHydroOptions");
  }
  core::MemoryGovernor& governor = *options.regrid_memory_governor;
  const std::uint64_t state_baseline = core::memoryReportBaselineOwnedBytes(
      core::collectSimulationMemoryReport(state));
  const core::MemoryGovernorSnapshot snapshot = governor.snapshot();
  if (snapshot.baseline_owned_bytes < state_baseline) {
    throw std::logic_error(
        "production AMR regrid governor baseline is stale: SimulationState owned capacity exceeds process baseline");
  }
  return ProductionAmrRegridAdmission{
      .reservation = governor.reserve(core::MemoryClass::kPhaseResident, requested_bytes, owner),
      .process_baseline_before = snapshot.baseline_owned_bytes,
      .state_baseline_before = state_baseline};
}

void commitProductionAmrCandidate(
    core::SimulationState& state,
    core::SimulationState&& candidate,
    ProductionAmrRegridAdmission& admission) {
  static_assert(std::is_nothrow_move_assignable_v<core::CellSoa>);
  static_assert(std::is_nothrow_move_assignable_v<core::GasCellSidecar>);
  static_assert(std::is_nothrow_move_assignable_v<core::PatchSoa>);
  static_assert(std::is_nothrow_move_assignable_v<core::GasCellIdentityMap>);

  if (state.gas_composition.enabled()) {
    // CHUI Part 2 fail-closed boundary: the conservative composition
    // prolongation/restriction kernels exist (physics::prolongComposition /
    // physics::restrictComposition) but regrid row construction does not yet
    // carry them; wiring them into runtime refinement is Part 3 (audit F6).
    // Refuse before any mutation instead of erasing species/dust lanes.
    throw std::runtime_error(
        "production AMR regrid with an enabled gas composition sidecar is not supported until Part 3 "
        "wires conservative composition prolongation/restriction into regrid; refusing before mutation");
  }
  admission.reservation.commit();
  state.cells = std::move(candidate.cells);
  state.gas_cells = std::move(candidate.gas_cells);
  state.patches = std::move(candidate.patches);
  state.gas_cell_identity = std::move(candidate.gas_cell_identity);
  state.bumpCellIndexGeneration();

  const std::uint64_t state_baseline_after = core::memoryReportBaselineOwnedBytes(
      core::collectSimulationMemoryReport(state));
  std::uint64_t process_baseline_after = admission.process_baseline_before;
  if (state_baseline_after >= admission.state_baseline_before) {
    process_baseline_after = core::checkedMemoryBytesAdd(
        process_baseline_after,
        state_baseline_after - admission.state_baseline_before,
        "production AMR regrid retained baseline growth");
  } else {
    const std::uint64_t released = admission.state_baseline_before - state_baseline_after;
    if (released > process_baseline_after) {
      throw std::logic_error(
          "production AMR regrid retained baseline release exceeds process baseline");
    }
    process_baseline_after -= released;
  }
  admission.reservation.reconcileBaselineOwnedAndRelease(process_baseline_after);
}

}  // namespace

bool patchStateRowHasExplicitGeometry(
    const core::SimulationState& state,
    std::size_t patch_index) {
  if (patch_index >= state.patches.size()) {
    return false;
  }
  return state.patches.extent_x_comoving[patch_index] > 0.0 &&
      state.patches.extent_y_comoving[patch_index] > 0.0 &&
      state.patches.extent_z_comoving[patch_index] > 0.0 &&
      state.patches.cell_dim_x[patch_index] != 0U &&
      state.patches.cell_dim_y[patch_index] != 0U &&
      state.patches.cell_dim_z[patch_index] != 0U;
}

PatchDescriptor patchDescriptorFromStateRow(
    const core::SimulationState& state,
    std::size_t patch_index) {
  if (patch_index >= state.patches.size()) {
    throw std::out_of_range("patchDescriptorFromStateRow: patch_index out of range");
  }
  if (!patchStateRowHasExplicitGeometry(state, patch_index)) {
    throw std::runtime_error("patchDescriptorFromStateRow: patch row lacks explicit restart-authoritative geometry");
  }
  PatchDescriptor descriptor;
  descriptor.patch_id = state.patches.patch_id[patch_index];
  descriptor.parent_patch_id = state.patches.parent_patch_id[patch_index];
  descriptor.level = static_cast<std::uint8_t>(std::max<std::int32_t>(state.patches.level[patch_index], 0));
  descriptor.morton_key = state.patches.morton_key[patch_index];
  descriptor.origin_comov = {
      state.patches.origin_x_comoving[patch_index],
      state.patches.origin_y_comoving[patch_index],
      state.patches.origin_z_comoving[patch_index]};
  descriptor.extent_comov = {
      state.patches.extent_x_comoving[patch_index],
      state.patches.extent_y_comoving[patch_index],
      state.patches.extent_z_comoving[patch_index]};
  descriptor.cell_dims = {
      state.patches.cell_dim_x[patch_index],
      state.patches.cell_dim_y[patch_index],
      state.patches.cell_dim_z[patch_index]};
  if (descriptor.patch_id == 0U || product(descriptor.cell_dims) != state.patches.cell_count[patch_index]) {
    throw std::runtime_error("patchDescriptorFromStateRow: invalid explicit patch descriptor lanes");
  }
  return descriptor;
}

void writePatchDescriptorToStateRow(
    core::SimulationState& state,
    std::size_t patch_index,
    const PatchDescriptor& patch) {
  if (patch_index >= state.patches.size()) {
    throw std::out_of_range("writePatchDescriptorToStateRow: patch_index out of range");
  }
  if (patch.patch_id == 0U || patch.extent_comov[0] <= 0.0 || patch.extent_comov[1] <= 0.0 ||
      patch.extent_comov[2] <= 0.0 || patch.cell_dims[0] == 0U || patch.cell_dims[1] == 0U ||
      patch.cell_dims[2] == 0U) {
    throw std::invalid_argument("writePatchDescriptorToStateRow: invalid explicit patch geometry");
  }
  state.patches.patch_id[patch_index] = patch.patch_id;
  state.patches.parent_patch_id[patch_index] = patch.parent_patch_id;
  state.patches.level[patch_index] = static_cast<std::int32_t>(patch.level);
  state.patches.morton_key[patch_index] = patch.morton_key;
  state.patches.origin_x_comoving[patch_index] = patch.origin_comov[0];
  state.patches.origin_y_comoving[patch_index] = patch.origin_comov[1];
  state.patches.origin_z_comoving[patch_index] = patch.origin_comov[2];
  state.patches.extent_x_comoving[patch_index] = patch.extent_comov[0];
  state.patches.extent_y_comoving[patch_index] = patch.extent_comov[1];
  state.patches.extent_z_comoving[patch_index] = patch.extent_comov[2];
  state.patches.cell_dim_x[patch_index] = patch.cell_dims[0];
  state.patches.cell_dim_y[patch_index] = patch.cell_dims[1];
  state.patches.cell_dim_z[patch_index] = patch.cell_dims[2];
}

std::uint8_t productionAmrHydroGhostLayers(const ProductionAmrHydroOptions& options) {
  return amrHydroGhostLayersForOptions(options);
}

bool hasProductionAmrHydroCoverage(const core::SimulationState& state) {
  if (state.cells.size() == 0U || state.patches.size() == 0U) {
    return false;
  }
  if (!state.gas_cell_identity.isConsistent() ||
      !state.gas_cell_identity.coversDenseLocalRows(state.cells.size()) ||
      !state.gasCellIdentityMapMatchesSidecarLanes()) {
    return false;
  }
  std::vector<std::uint8_t> covered(state.cells.size(), 0U);
  for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
    const std::uint32_t first = state.patches.first_cell[patch];
    const std::uint32_t count = state.patches.cell_count[patch];
    if (count == 0U) {
      continue;
    }
    if (!patchStateRowHasExplicitGeometry(state, patch)) {
      return false;
    }
    if (static_cast<std::size_t>(state.patches.cell_dim_x[patch]) *
            static_cast<std::size_t>(state.patches.cell_dim_y[patch]) *
            static_cast<std::size_t>(state.patches.cell_dim_z[patch]) != count) {
      return false;
    }
    if (first + count > state.cells.size()) {
      return false;
    }
    for (std::uint32_t offset = 0; offset < count; ++offset) {
      const std::uint32_t row = first + offset;
      if (covered[row] != 0U || state.cells.patch_index[row] != patch) {
        return false;
      }
      const auto* record = state.gas_cell_identity.findByLocalRow(row);
      if (record == nullptr || record->owning_patch_id != state.patches.patch_id[patch]) {
        return false;
      }
      covered[row] = 1U;
    }
  }
  return std::all_of(covered.begin(), covered.end(), [](std::uint8_t flag) { return flag != 0U; });
}

std::vector<PatchDescriptor> buildProductionAmrPatchDescriptors(const core::SimulationState& state) {
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error(
        "buildProductionAmrPatchDescriptors requires complete SimulationState AMR patch coverage with explicit patch geometry");
  }
  std::vector<PatchDescriptor> descriptors;
  descriptors.reserve(state.patches.size());
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.cell_count[patch_index] == 0U) {
      continue;
    }
    descriptors.push_back(patchDescriptorFromStateRow(state, patch_index));
  }
  return descriptors;
}

void populateAmrHydroFluxRegisterFaces(
    AmrHydroPatchGeometry& patch_geometry,
    std::span<const PatchDescriptor> all_patches) {
  populateAmrHydroFluxRegisterFaces(patch_geometry, all_patches, {0.0, 0.0, 0.0});
}

void populateAmrHydroFluxRegisterFaces(
    AmrHydroPatchGeometry& patch_geometry,
    std::span<const PatchDescriptor> all_patches,
    const std::array<double, 3>& periodic_box_extent_comoving) {
  const internal::AmrPeriodicFrame frame = internal::makeAmrPeriodicFrame(
      periodic_box_extent_comoving[0] > 0.0 || periodic_box_extent_comoving[1] > 0.0 ||
          periodic_box_extent_comoving[2] > 0.0,
      periodic_box_extent_comoving);
  patch_geometry.geometry.flux_register_faces.assign(
      patch_geometry.geometry.faces.size(),
      hydro::HydroFluxRegisterFace{});
  // Deep stencil-only ghosts own no face and no register; they only need the
  // (possibly different, possibly wrapped) source patch so that their stable
  // source gas-cell identity and dual-energy source can be traced.
  for (AmrHydroGhostDescriptor& deep : patch_geometry.ghosts) {
    if (deep.ghost_layer <= 1U || deep.boundary_class == AmrHydroBoundaryClass::kPhysical) {
      continue;
    }
    const hydro::HydroGhostCell& ghost = patch_geometry.geometry.ghost_cells.at(deep.ghost_slot);
    const std::optional<PatchDescriptor> source_patch =
        sourcePatchForGhost(patch_geometry, ghost, all_patches, frame, deep.ghost_layer);
    deep.source_patch_id = source_patch.has_value() ? source_patch->patch_id : 0U;
  }
  for (std::size_t face_index = 0; face_index < patch_geometry.geometry.faces.size(); ++face_index) {
    const hydro::HydroFace& face = patch_geometry.geometry.faces[face_index];
    if (face.ghost_cell_slot == hydro::k_invalid_ghost_cell_slot) {
      continue;
    }
    const auto ghost_it = std::find_if(
        patch_geometry.ghosts.begin(),
        patch_geometry.ghosts.end(),
        [&face](const AmrHydroGhostDescriptor& ghost) { return ghost.ghost_slot == face.ghost_cell_slot; });
    if (ghost_it == patch_geometry.ghosts.end() ||
        ghost_it->boundary_class == AmrHydroBoundaryClass::kPhysical) {
      continue;
    }
    const hydro::HydroGhostCell& ghost = patch_geometry.geometry.ghost_cells.at(face.ghost_cell_slot);
    const std::optional<PatchDescriptor> source_patch =
        sourcePatchForGhost(patch_geometry, ghost, all_patches, frame);
    if (!source_patch.has_value()) {
      continue;
    }
    AmrHydroGhostDescriptor& ghost_descriptor = *ghost_it;
    ghost_descriptor.source_patch_id = source_patch->patch_id;
    const std::array<double, 3> probe = ghostProbePoint(patch_geometry, ghost, 1U, frame);
    [[maybe_unused]] const std::size_t source_cell = cellContainingPoint(*source_patch, probe);
    // Source gas-cell IDs are patched after all local AMR hydro geometries exist, because this
    // metadata belongs to the source patch's gas_cell_id vector, not the target patch.
    if (ghost_descriptor.boundary_class != AmrHydroBoundaryClass::kCoarseFine) {
      continue;
    }

    const bool target_is_coarse = patch_geometry.patch.level < source_patch->level;
    const PatchDescriptor& coarse_patch = target_is_coarse ? patch_geometry.patch : *source_patch;
    const std::size_t coarse_cell = target_is_coarse ? face.owner_cell : cellContainingPoint(coarse_patch, probe);
    const hydro::HydroFaceSide coarse_orientation = target_is_coarse ? ghost.side : oppositeSide(ghost.side);
    patch_geometry.geometry.flux_register_faces[face_index] = hydro::HydroFluxRegisterFace{
        .role = target_is_coarse ? hydro::HydroFluxRegisterFaceRole::kCoarse : hydro::HydroFluxRegisterFaceRole::kFine,
        .register_key = makeRegisterKey(coarse_patch.patch_id, coarse_cell, ghost.axis, coarse_orientation),
        .coarse_patch_id = coarse_patch.patch_id,
        .coarse_gas_cell_id = target_is_coarse ? patch_geometry.gas_cell_ids.at(coarse_cell) : ghost_descriptor.source_gas_cell_id,
        .coarse_cell_index = coarse_cell,
        .level = static_cast<int>(coarse_patch.level),
        .axis = ghost.axis,
        .orientation = coarse_orientation,
        .coarse_orientation_sign = target_is_coarse ? 1.0 : -1.0};
  }
}

namespace {
// A single sorted mirror index replaces the two population-scale hash maps
// previously rebuilt for every patch scatter. Duplicate-parent semantics are
// preserved: only a parent referenced by exactly one gas cell is mirrored.
struct AmrParentMirrorIndex {
  using Row = std::pair<std::uint64_t, std::uint32_t>;
  std::pmr::vector<Row> rows;
  std::pmr::vector<std::uint64_t> parent_ids;

  explicit AmrParentMirrorIndex(std::pmr::memory_resource* resource)
      : rows(resource), parent_ids(resource) {}

  [[nodiscard]] static std::uint64_t workspaceBytes(const core::SimulationState& state) {
    return core::checkedMemoryBytesAdd(
        static_cast<std::uint64_t>(core::checkedSizeMultiply(
            state.particles.size(), sizeof(Row), "AMR parent mirror rows")),
        core::checkedMemoryBytesAdd(
            static_cast<std::uint64_t>(core::checkedSizeMultiply(
                state.gas_cell_identity.records().size(), sizeof(std::uint64_t),
                "AMR parent mirror references")),
            256U, "AMR parent mirror alignment"), "AMR parent mirror workspace");
  }

  void build(const core::SimulationState& state) {
    rows.reserve(state.particles.size());
    parent_ids.reserve(state.gas_cell_identity.records().size());
    for (std::uint32_t row = 0; row < state.particles.size(); ++row) {
      rows.emplace_back(state.particle_sidecar.particle_id[row], row);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& record : state.gas_cell_identity.records()) {
      if (record.parent_particle_id.has_value()) parent_ids.push_back(*record.parent_particle_id);
    }
    std::sort(parent_ids.begin(), parent_ids.end());
  }

  [[nodiscard]] std::optional<std::uint32_t> uniqueParentRow(std::uint64_t id) const {
    const auto refs = std::equal_range(parent_ids.begin(), parent_ids.end(), id);
    if (refs.second - refs.first != 1) return std::nullopt;
    const auto it = std::lower_bound(rows.begin(), rows.end(), Row{id, 0U});
    if (it == rows.end() || it->first != id) return std::nullopt;
    return it->second;
  }
};

void scatterAmrHydroConservedStateWithIndex(
    core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroConservedStateSoa& conserved,
    double adiabatic_index,
    const AmrParentMirrorIndex& parent_index,
    const AmrThermodynamicCachePolicy& thermodynamic_cache = {}) {
  state.requireGasCellIdentityMapFresh(
      patch_geometry.source_gas_cell_identity_generation,
      "scatterAmrHydroConservedState");
  if (conserved.size() < patch_geometry.geometry.cellCount()) {
    throw std::invalid_argument("scatterAmrHydroConservedState: conserved storage lacks real AMR cells");
  }
  for (std::size_t patch_cell = 0; patch_cell < patch_geometry.real_cells.size(); ++patch_cell) {
    const AmrHydroCellDescriptor& cell = patch_geometry.real_cells[patch_cell];
    const auto row = state.rowForGasCellId(cell.gas_cell_id);
    if (!row.has_value() || *row != cell.local_cell_row) {
      throw std::runtime_error("scatterAmrHydroConservedState: gas_cell_id resolved to a different row");
    }
    const hydro::HydroPrimitiveState primitive =
        hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(patch_cell), adiabatic_index);
    const double rho = std::max(primitive.rho_comoving, 1.0e-14);
    double pressure = std::max(primitive.pressure_comoving, 1.0e-14);
    if (conserved.auxiliaryInternalEnergyEnabled()) {
      // Dual energy: store the solver-selected thermal energy.
      const double rho_e = conserved.auxiliaryInternalEnergyDensityComoving()[patch_cell];
      if (std::isfinite(rho_e) && rho_e > 0.0) {
        pressure = std::max((adiabatic_index - 1.0) * rho_e, 1.0e-14);
      }
    }
    state.gas_cells.density_code[*row] = rho;
    state.gas_cells.pressure_code[*row] = pressure;
    state.gas_cells.internal_energy_code[*row] = pressure / ((adiabatic_index - 1.0) * rho);
    state.gas_cells.velocity_x_peculiar[*row] = primitive.vel_x_peculiar;
    state.gas_cells.velocity_y_peculiar[*row] = primitive.vel_y_peculiar;
    state.gas_cells.velocity_z_peculiar[*row] = primitive.vel_z_peculiar;
    state.gas_cells.sound_speed_code[*row] = std::sqrt(std::max(0.0, adiabatic_index * pressure / rho));
    storeDerivedGasThermodynamics(state, *row, adiabatic_index, thermodynamic_cache);
    state.cells.mass_code[*row] = rho * patch_geometry.geometry.cell_volume_comoving;
    state.gas_cells.metal_mass_code[*row] =
        std::clamp(primitive.metallicity_mass_fraction, 0.0, 1.0) * state.cells.mass_code[*row];

    const auto parent_id = state.parentParticleIdForGasCellId(cell.gas_cell_id);
    if (parent_id.has_value()) {
      if (const auto parent_row = parent_index.uniqueParentRow(*parent_id);
          parent_row.has_value()) {
        state.particles.mass_code[*parent_row] = state.cells.mass_code[*row];
        state.particles.velocity_x_peculiar[*parent_row] = primitive.vel_x_peculiar;
        state.particles.velocity_y_peculiar[*parent_row] = primitive.vel_y_peculiar;
        state.particles.velocity_z_peculiar[*parent_row] = primitive.vel_z_peculiar;
      }
    }
  }
}

}  // namespace

void scatterAmrHydroConservedState(
    core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroConservedStateSoa& conserved,
    double adiabatic_index) {
  core::GovernedScratchArena arena(nullptr, core::MemoryClass::kScratchArena,
      AmrParentMirrorIndex::workspaceBytes(state), "amr.parent_mirror.compatibility");
  AmrParentMirrorIndex index(arena.resource());
  index.build(state);
  scatterAmrHydroConservedStateWithIndex(state, patch_geometry, conserved,
      adiabatic_index, index);
}

RefluxDiagnostics applyFluxRegistersToSimulationState(
    core::SimulationState& state,
    std::span<const FluxRegisterEntry> entries,
    std::span<const PatchDescriptor> all_patches,
    double adiabatic_index,
    const AmrThermodynamicCachePolicy& thermodynamic_cache) {
  // Gate-2 failure atomicity: phase 1 fully validates every record that would
  // be committed before ANY conserved-state mutation. The stale-mapping throw
  // below is the only throwing rejection; all other rejections are counted as
  // skips. Validating first guarantees a later record cannot throw after
  // earlier records already mutated canonical gas state.
  for (const FluxRegisterEntry& entry : entries) {
    if (!entry.isComplete()) {
      continue;
    }
    const double area_scale = std::max({1.0, std::abs(entry.coarse_area_comov), std::abs(entry.fine_area_comov)});
    if (entry.coarse_area_comov <= 0.0 || entry.fine_area_comov <= 0.0 ||
        std::abs(entry.coarse_area_comov - entry.fine_area_comov) > 1.0e-10 * area_scale) {
      continue;
    }
    if (entry.coarse_gas_cell_id == 0U) {
      continue;
    }
    const auto patch_it = std::find_if(
        all_patches.begin(),
        all_patches.end(),
        [&entry](const PatchDescriptor& patch) { return patch.patch_id == entry.coarse_patch_id; });
    if (patch_it == all_patches.end()) {
      continue;
    }
    const auto row_opt = state.gas_cell_identity.rowForGasCellId(entry.coarse_gas_cell_id);
    if (!row_opt.has_value()) {
      continue;
    }
    const std::uint32_t row = *row_opt;
    const auto* record = state.gas_cell_identity.findByGasCellId(entry.coarse_gas_cell_id);
    if (record == nullptr || record->owning_patch_id != entry.coarse_patch_id) {
      continue;
    }
    const auto patch_index = patchIndexById(state, entry.coarse_patch_id);
    if (!patch_index.has_value() || row >= state.cells.patch_index.size() ||
        state.cells.patch_index[row] != *patch_index ||
        state.gas_cells.gas_cell_id[row] != entry.coarse_gas_cell_id ||
        entry.coarse_cell_index >= product(patch_it->cell_dims) ||
        patchLocalCellForRow(state, *patch_it, row, "applyFluxRegistersToSimulationState").linear_index != entry.coarse_cell_index) {
      throw std::runtime_error("applyFluxRegistersToSimulationState rejected a stale AMR reflux target mapping");
    }
  }
  RefluxDiagnostics diagnostics;
  ProductionAmrHydroOptions options;
  options.adiabatic_index = adiabatic_index;
  for (const FluxRegisterEntry& entry : entries) {
    if (!entry.isComplete()) {
      ++diagnostics.skipped_incomplete_register_count;
      continue;
    }
    const double area_scale = std::max({1.0, std::abs(entry.coarse_area_comov), std::abs(entry.fine_area_comov)});
    if (entry.coarse_area_comov <= 0.0 || entry.fine_area_comov <= 0.0 ||
        std::abs(entry.coarse_area_comov - entry.fine_area_comov) > 1.0e-10 * area_scale) {
      ++diagnostics.skipped_area_mismatch_count;
      continue;
    }
    if (entry.coarse_gas_cell_id == 0U) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const auto patch_it = std::find_if(
        all_patches.begin(),
        all_patches.end(),
        [&entry](const PatchDescriptor& patch) { return patch.patch_id == entry.coarse_patch_id; });
    if (patch_it == all_patches.end()) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const auto row_opt = state.gas_cell_identity.rowForGasCellId(entry.coarse_gas_cell_id);
    if (!row_opt.has_value()) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const std::uint32_t row = *row_opt;
    const auto* record = state.gas_cell_identity.findByGasCellId(entry.coarse_gas_cell_id);
    if (record == nullptr || record->owning_patch_id != entry.coarse_patch_id) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const auto patch_index = patchIndexById(state, entry.coarse_patch_id);
    if (!patch_index.has_value() || row >= state.cells.patch_index.size() ||
        state.cells.patch_index[row] != *patch_index ||
        state.gas_cells.gas_cell_id[row] != entry.coarse_gas_cell_id ||
        entry.coarse_cell_index >= product(patch_it->cell_dims) ||
        patchLocalCellForRow(state, *patch_it, row, "applyFluxRegistersToSimulationState").linear_index != entry.coarse_cell_index) {
      throw std::runtime_error("applyFluxRegistersToSimulationState rejected a stale AMR reflux target mapping");
    }
    const double volume = patch_it->extent_comov[0] * patch_it->extent_comov[1] * patch_it->extent_comov[2] /
        static_cast<double>(std::max<std::size_t>(product(patch_it->cell_dims), 1U));
    ++diagnostics.complete_register_count;
    const hydro::HydroPrimitiveState old_primitive = primitiveForRow(
        state, row, 1.0e-14, 1.0e-14, adiabatic_index);
    hydro::HydroConservedState conserved =
        hydro::HydroCoreSolver::conservedFromPrimitive(old_primitive, adiabatic_index);
    ConservedState delta_flux = entry.fine_face_flux_code - entry.coarse_face_flux_code;
    delta_flux *= (entry.coarse_area_comov * entry.dt_code / volume);
    conserved.mass_density_comoving -= delta_flux.mass_code;
    conserved.momentum_density_x_comoving -= delta_flux.momentum_x_code;
    conserved.momentum_density_y_comoving -= delta_flux.momentum_y_code;
    conserved.momentum_density_z_comoving -= delta_flux.momentum_z_code;
    conserved.total_energy_density_comoving -= delta_flux.total_energy_code;
    conserved.metal_mass_density_comoving -= delta_flux.metal_mass_code;
    const double old_internal_density = old_primitive.pressure_comoving / (adiabatic_index - 1.0);
    hydro::HydroPrimitiveState primitive =
        hydro::HydroCoreSolver::primitiveFromConserved(conserved, adiabatic_index);
    if (thermodynamic_cache.dual_energy) {
      const double kinetic = 0.5 *
          (conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
           conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
           conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving) /
          std::max(conserved.mass_density_comoving, 1.0e-14);
      const double thermal_conservative = conserved.total_energy_density_comoving - kinetic;
      if (!(thermal_conservative > 0.0) ||
          thermal_conservative < thermodynamic_cache.dual_energy_eta * conserved.total_energy_density_comoving) {
        // Keep the pre-reflux specific internal energy (the corrected E - K
        // is cancellation-dominated); mass, momentum, E and metals are still
        // corrected conservatively. Ledger the thermal correction.
        const double kept_e = old_primitive.pressure_comoving /
            ((adiabatic_index - 1.0) * std::max(old_primitive.rho_comoving, 1.0e-14));
        const double kept_rho_e = std::max(primitive.rho_comoving, 1.0e-14) * kept_e;
        primitive.pressure_comoving = std::max((adiabatic_index - 1.0) * kept_rho_e, 1.0e-14);
        ++diagnostics.dual_energy_reflux_cells;
        diagnostics.dual_energy_reflux_thermal_correction_code += (kept_rho_e - thermal_conservative) * volume;
      }
    }
    state.gas_cells.density_code[row] = primitive.rho_comoving;
    state.gas_cells.pressure_code[row] = primitive.pressure_comoving;
    state.gas_cells.internal_energy_code[row] =
        primitive.pressure_comoving / ((adiabatic_index - 1.0) * std::max(primitive.rho_comoving, 1.0e-14));
    state.gas_cells.velocity_x_peculiar[row] = primitive.vel_x_peculiar;
    state.gas_cells.velocity_y_peculiar[row] = primitive.vel_y_peculiar;
    state.gas_cells.velocity_z_peculiar[row] = primitive.vel_z_peculiar;
    state.gas_cells.sound_speed_code[row] = std::sqrt(std::max(0.0, adiabatic_index * primitive.pressure_comoving /
        std::max(primitive.rho_comoving, 1.0e-14)));
    storeDerivedGasThermodynamics(state, row, adiabatic_index, thermodynamic_cache);
    state.cells.mass_code[row] = primitive.rho_comoving * volume;
    state.gas_cells.metal_mass_code[row] =
        std::clamp(primitive.metallicity_mass_fraction, 0.0, 1.0) * state.cells.mass_code[row];
    const double new_internal_density = primitive.pressure_comoving / (adiabatic_index - 1.0);

    diagnostics.corrected_cells += 1U;
    diagnostics.corrected_mass_code += std::abs(delta_flux.mass_code * volume);
    diagnostics.corrected_momentum_x_code += std::abs(delta_flux.momentum_x_code * volume);
    diagnostics.corrected_momentum_y_code += std::abs(delta_flux.momentum_y_code * volume);
    diagnostics.corrected_momentum_z_code += std::abs(delta_flux.momentum_z_code * volume);
    diagnostics.corrected_total_energy_code += std::abs(delta_flux.total_energy_code * volume);
    diagnostics.corrected_metal_mass_code += std::abs(delta_flux.metal_mass_code * volume);
    diagnostics.corrected_energy_code += std::abs(delta_flux.total_energy_code * volume);
    diagnostics.corrected_internal_energy_code += std::abs((new_internal_density - old_internal_density) * volume);
  }

  return diagnostics;
}

namespace {

[[nodiscard]] std::uint8_t axisToStorage(hydro::HydroFaceAxis axis) {
  return static_cast<std::uint8_t>(axisIndex(axis));
}

[[nodiscard]] std::uint8_t sideToStorage(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? 0U : 1U;
}

[[nodiscard]] hydro::HydroFaceAxis axisFromStorage(std::uint8_t axis) {
  switch (axis) {
    case 0U:
      return hydro::HydroFaceAxis::kX;
    case 1U:
      return hydro::HydroFaceAxis::kY;
    case 2U:
      return hydro::HydroFaceAxis::kZ;
    default:
      throw std::runtime_error("pending AMR flux register has invalid face axis metadata");
  }
}

[[nodiscard]] hydro::HydroFaceSide sideFromStorage(std::uint8_t side) {
  if (side > 1U) {
    throw std::runtime_error("pending AMR flux register has invalid face-side metadata");
  }
  return side == 0U ? hydro::HydroFaceSide::kLower : hydro::HydroFaceSide::kUpper;
}

void addFluxIntegral(
    double& mass,
    double& momentum_x,
    double& momentum_y,
    double& momentum_z,
    double& total_energy,
    double& metal_mass,
    const ConservedState& flux,
    double area_comov,
    double dt_code) {
  const double scale = area_comov * dt_code;
  mass += flux.mass_code * scale;
  momentum_x += flux.momentum_x_code * scale;
  momentum_y += flux.momentum_y_code * scale;
  momentum_z += flux.momentum_z_code * scale;
  total_energy += flux.total_energy_code * scale;
  metal_mass += flux.metal_mass_code * scale;
}

[[nodiscard]] ConservedState averageFluxFromIntegral(
    double mass,
    double momentum_x,
    double momentum_y,
    double momentum_z,
    double total_energy,
    double metal_mass,
    double area_comov,
    double dt_code) {
  const double denom = area_comov * dt_code;
  if (denom <= 0.0) {
    return ConservedState{};
  }
  return ConservedState{
      .mass_code = mass / denom,
      .momentum_x_code = momentum_x / denom,
      .momentum_y_code = momentum_y / denom,
      .momentum_z_code = momentum_z / denom,
      .total_energy_code = total_energy / denom,
      .metal_mass_code = metal_mass / denom};
}

void validatePendingCompatible(
    const core::PendingFluxRegisterRecord& pending,
    const FluxRegisterEntry& entry) {
  if (pending.coarse_patch_id != entry.coarse_patch_id ||
      pending.coarse_gas_cell_id != entry.coarse_gas_cell_id ||
      pending.coarse_cell_index != entry.coarse_cell_index ||
      pending.level != entry.level ||
      pending.axis != axisToStorage(entry.axis) ||
      pending.orientation != sideToStorage(entry.orientation)) {
    throw std::runtime_error("pending AMR flux register received incompatible stable identity metadata for an existing key");
  }
}

[[nodiscard]] core::PendingFluxRegisterRecord makePendingRecordSeed(
    const core::SimulationState& state,
    const FluxRegisterEntry& entry,
    const ProductionAmrHydroOptions& options) {
  const double coarse_dt = options.reflux_coarse_dt_code > 0.0 ? options.reflux_coarse_dt_code : entry.dt_code;
  const double interval_start = options.reflux_interval_start_code;
  const double interval_end = options.reflux_interval_end_code > interval_start
      ? options.reflux_interval_end_code
      : interval_start + coarse_dt;
  const double expected_area = std::max({entry.face_area_comov, entry.coarse_area_comov, entry.fine_area_comov, 0.0});
  return core::PendingFluxRegisterRecord{
      .register_key = entry.register_key,
      .coarse_patch_id = entry.coarse_patch_id,
      .coarse_gas_cell_id = entry.coarse_gas_cell_id,
      .coarse_cell_index = entry.coarse_cell_index,
      .level = entry.level,
      .axis = axisToStorage(entry.axis),
      .orientation = sideToStorage(entry.orientation),
      .expected_area_comov = expected_area,
      .interval_start_code = interval_start,
      .interval_end_code = interval_end,
      .coarse_dt_code = coarse_dt,
      .expected_fine_substeps = std::max<std::uint32_t>(1U, options.expected_fine_substeps),
      .gas_cell_identity_generation = state.gasCellIdentityGeneration(),
      .patch_geometry_generation = state.cellIndexGeneration()};
}

[[nodiscard]] FluxRegisterEntry fluxEntryFromPendingRecord(
    const core::PendingFluxRegisterRecord& pending) {
  const double area = pending.expected_area_comov;
  const double dt = pending.coarse_dt_code > 0.0
      ? pending.coarse_dt_code
      : std::max(0.0, pending.interval_end_code - pending.interval_start_code);
  return FluxRegisterEntry{
      .register_key = pending.register_key,
      .coarse_patch_id = pending.coarse_patch_id,
      .coarse_gas_cell_id = pending.coarse_gas_cell_id,
      .coarse_cell_index = pending.coarse_cell_index,
      .level = pending.level,
      .axis = axisFromStorage(pending.axis),
      .orientation = sideFromStorage(pending.orientation),
      .coarse_face_flux_code = averageFluxFromIntegral(
          pending.coarse_mass_flux_integral_code,
          pending.coarse_momentum_x_flux_integral_code,
          pending.coarse_momentum_y_flux_integral_code,
          pending.coarse_momentum_z_flux_integral_code,
          pending.coarse_total_energy_flux_integral_code,
          pending.coarse_metal_mass_flux_integral_code,
          area,
          dt),
      .fine_face_flux_code = averageFluxFromIntegral(
          pending.fine_mass_flux_integral_code,
          pending.fine_momentum_x_flux_integral_code,
          pending.fine_momentum_y_flux_integral_code,
          pending.fine_momentum_z_flux_integral_code,
          pending.fine_total_energy_flux_integral_code,
          pending.fine_metal_mass_flux_integral_code,
          area,
          dt),
      .face_area_comov = area,
      .coarse_area_comov = pending.coarse_area_accumulated_comov,
      .fine_area_comov = pending.fine_area_accumulated_comov,
      .dt_code = dt,
      .coarse_face_count = pending.coarse_face_count,
      .fine_face_count = pending.fine_face_count};
}

void accumulateDiagnostics(ProductionAmrHydroDiagnostics& lhs, const ProductionAmrHydroDiagnostics& rhs) {
  lhs.patch_count = std::max(lhs.patch_count, rhs.patch_count);
  lhs.advanced_patch_count += rhs.advanced_patch_count;
  lhs.riemann_solve_count += rhs.riemann_solve_count;
  lhs.positivity_fallback_count += rhs.positivity_fallback_count;
  lhs.stencil_degradation_count += rhs.stencil_degradation_count;
  lhs.ctu_transverse_degraded_states += rhs.ctu_transverse_degraded_states;
  lhs.dual_energy_auxiliary_cells += rhs.dual_energy_auxiliary_cells;
  lhs.dual_energy_unresolved_ghosts += rhs.dual_energy_unresolved_ghosts;
  lhs.dual_energy_thermal_correction_code += rhs.dual_energy_thermal_correction_code;
  lhs.active_cell_count += rhs.active_cell_count;
  lhs.active_face_count += rhs.active_face_count;
  lhs.flux_register_entry_count += rhs.flux_register_entry_count;
  lhs.pending_register_created_count += rhs.pending_register_created_count;
  lhs.pending_register_completed_count += rhs.pending_register_completed_count;
  lhs.pending_register_deferred_count += rhs.pending_register_deferred_count;
  lhs.pending_register_applied_count += rhs.pending_register_applied_count;
  lhs.pending_register_rejected_count += rhs.pending_register_rejected_count;
  lhs.ghost_fill.physical_ghosts_filled += rhs.ghost_fill.physical_ghosts_filled;
  lhs.ghost_fill.same_level_ghosts_filled += rhs.ghost_fill.same_level_ghosts_filled;
  lhs.ghost_fill.coarse_to_fine_ghosts_filled += rhs.ghost_fill.coarse_to_fine_ghosts_filled;
  lhs.ghost_fill.fine_to_coarse_ghosts_filled += rhs.ghost_fill.fine_to_coarse_ghosts_filled;
  lhs.ghost_fill.temporal_coarse_to_fine_ghosts_filled += rhs.ghost_fill.temporal_coarse_to_fine_ghosts_filled;
  lhs.ghost_fill.temporal_endpoint_ghosts_filled += rhs.ghost_fill.temporal_endpoint_ghosts_filled;
  lhs.ghost_fill.skipped_remote_ghosts += rhs.ghost_fill.skipped_remote_ghosts;
  lhs.ghost_fill.stale_epoch_rejections += rhs.ghost_fill.stale_epoch_rejections;
  lhs.ghost_fill.temporal_same_level_mismatch_rejections += rhs.ghost_fill.temporal_same_level_mismatch_rejections;
  lhs.ghost_fill.temporal_history_missing_rejections += rhs.ghost_fill.temporal_history_missing_rejections;
  lhs.ghost_fill.temporal_history_invalid_rejections += rhs.ghost_fill.temporal_history_invalid_rejections;
  lhs.ghost_fill.temporal_time_out_of_range_rejections += rhs.ghost_fill.temporal_time_out_of_range_rejections;
  lhs.ghost_fill.temporal_fine_to_coarse_misalignment_rejections += rhs.ghost_fill.temporal_fine_to_coarse_misalignment_rejections;
  lhs.ghost_fill.temporal_geometry_mismatch_rejections += rhs.ghost_fill.temporal_geometry_mismatch_rejections;
  lhs.ghost_fill.temporal_identity_mismatch_rejections += rhs.ghost_fill.temporal_identity_mismatch_rejections;
  lhs.ghost_fill.missing_source_records += rhs.ghost_fill.missing_source_records;
  lhs.ghost_fill.unresolved_ghosts += rhs.ghost_fill.unresolved_ghosts;
  lhs.reflux.complete_register_count += rhs.reflux.complete_register_count;
  lhs.reflux.skipped_incomplete_register_count += rhs.reflux.skipped_incomplete_register_count;
  lhs.reflux.skipped_area_mismatch_count += rhs.reflux.skipped_area_mismatch_count;
  lhs.reflux.skipped_missing_target_count += rhs.reflux.skipped_missing_target_count;
  lhs.reflux.corrected_cells += rhs.reflux.corrected_cells;
  lhs.reflux.corrected_mass_code += rhs.reflux.corrected_mass_code;
  lhs.reflux.corrected_momentum_x_code += rhs.reflux.corrected_momentum_x_code;
  lhs.reflux.corrected_momentum_y_code += rhs.reflux.corrected_momentum_y_code;
  lhs.reflux.corrected_momentum_z_code += rhs.reflux.corrected_momentum_z_code;
  lhs.reflux.corrected_total_energy_code += rhs.reflux.corrected_total_energy_code;
  lhs.reflux.corrected_metal_mass_code += rhs.reflux.corrected_metal_mass_code;
  lhs.reflux.corrected_energy_code += rhs.reflux.corrected_energy_code;
  lhs.reflux.corrected_internal_energy_code += rhs.reflux.corrected_internal_energy_code;
  lhs.reflux.dual_energy_reflux_cells += rhs.reflux.dual_energy_reflux_cells;
  lhs.reflux.dual_energy_reflux_thermal_correction_code += rhs.reflux.dual_energy_reflux_thermal_correction_code;
}

[[nodiscard]] std::vector<std::uint32_t> activeRowsForLevel(
    const core::SimulationState& state,
    int level,
    std::span<const std::uint32_t> requested_rows) {
  const bool all_requested = requested_rows.empty();
  std::vector<std::uint32_t> rows;
  if (!all_requested) {
    // When a specific active set is requested, iterate the requested rows
    // directly and filter by patch level. This avoids constructing a
    // population-scale hash set; the scan cost is O(|requested_rows|).
    rows.reserve(requested_rows.size());
    for (const std::uint32_t row : requested_rows) {
      if (row >= state.cells.patch_index.size()) {
        continue;
      }
      const std::uint32_t patch_index = state.cells.patch_index[row];
      if (patch_index < state.patches.size() && state.patches.level[patch_index] == level) {
        rows.push_back(row);
      }
    }
  } else {
    // When all cells are active, scan once without a hash set.
    rows.reserve(state.cells.size());
    for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
      if (row >= state.cells.patch_index.size()) {
        continue;
      }
      const std::uint32_t patch_index = state.cells.patch_index[row];
      if (patch_index < state.patches.size() && state.patches.level[patch_index] == level) {
        rows.push_back(row);
      }
    }
  }
  return rows;
}

}  // namespace

std::uint64_t activeLevelRowStagingBytes(
    std::size_t coarse_capacity,
    std::size_t fine_capacity) {
  const std::size_t coarse_bytes = core::checkedSizeMultiply(
      coarse_capacity, sizeof(std::uint32_t), "AMR active coarse row staging bytes");
  const std::size_t fine_bytes = core::checkedSizeMultiply(
      fine_capacity, sizeof(std::uint32_t), "AMR active fine row staging bytes");
  return core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeAdd(coarse_bytes, fine_bytes, "AMR active level row staging peak"),
      "AMR active level row staging width");
}

void fillActiveRowsForLevelInto(
    const core::SimulationState& state,
    int level,
    std::span<const std::uint32_t> requested_rows,
    std::vector<std::uint32_t>& out_rows) {
  out_rows.clear();
  const bool all_requested = requested_rows.empty();
  if (!all_requested) {
    for (const std::uint32_t row : requested_rows) {
      if (row >= state.cells.patch_index.size()) {
        continue;
      }
      const std::uint32_t patch_index = state.cells.patch_index[row];
      if (patch_index < state.patches.size() && state.patches.level[patch_index] == level) {
        out_rows.push_back(row);
      }
    }
  } else {
    for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
      if (row >= state.cells.patch_index.size()) {
        continue;
      }
      const std::uint32_t patch_index = state.cells.patch_index[row];
      if (patch_index < state.patches.size() && state.patches.level[patch_index] == level) {
        out_rows.push_back(row);
      }
    }
  }
}

void prepareAmrActiveLevelRows(
    const core::SimulationState& state,
    std::span<const std::uint32_t> requested_rows,
    int min_level,
    int max_level,
    core::MemoryGovernor* governor,
    AmrActiveLevelRowWorkspace& workspace) {
  // Exact-count admission: scan (allocation-free) to determine per-level
  // needs, then admit EXACTLY coarse_need+fine_need. Both lists coexist, so
  // admission includes BOTH physical capacities. This is tighter than the
  // 2x worst-case bound and preserves headroom for the sweep workspace below.
  std::size_t coarse_need = 0U;
  std::size_t fine_need = 0U;
  if (requested_rows.empty()) {
    for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
      if (row >= state.cells.patch_index.size()) {
        continue;
      }
      const std::uint32_t patch_index = state.cells.patch_index[row];
      if (patch_index >= state.patches.size()) {
        continue;
      }
      const int row_level = static_cast<int>(state.patches.level[patch_index]);
      if (row_level == min_level) {
        ++coarse_need;
      }
      if (row_level == max_level) {
        ++fine_need;
      }
    }
  } else {
    for (const std::uint32_t row : requested_rows) {
      if (row >= state.cells.patch_index.size()) {
        continue;
      }
      const std::uint32_t patch_index = state.cells.patch_index[row];
      if (patch_index >= state.patches.size()) {
        continue;
      }
      const int row_level = static_cast<int>(state.patches.level[patch_index]);
      if (row_level == min_level) {
        ++coarse_need;
      }
      if (row_level == max_level) {
        ++fine_need;
      }
    }
  }
  const std::uint64_t required_bytes =
      activeLevelRowStagingBytes(coarse_need, fine_need);
  const std::uint64_t retained_bytes = core::checkedMemoryBytesAdd(
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          workspace.coarse_rows.capacity(), sizeof(std::uint32_t),
          "AMR active coarse retained bytes")),
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          workspace.fine_rows.capacity(), sizeof(std::uint32_t),
          "AMR active fine retained bytes")),
      "AMR active level retained bytes");
  // Reuse retained capacity when sufficient: no new admission, no growth.
  if (retained_bytes >= required_bytes && required_bytes > 0U) {
    fillActiveRowsForLevelInto(state, min_level, requested_rows, workspace.coarse_rows);
    fillActiveRowsForLevelInto(state, max_level, requested_rows, workspace.fine_rows);
    return;
  }
  // Otherwise admit growth BEFORE vector growth. The reservation covers both
  // physical capacities simultaneously.
  core::MemoryReservation growth;
  if (governor != nullptr && required_bytes > 0U) {
    growth = governor->reserve(
        core::MemoryClass::kPhaseResident, required_bytes,
        "amr.active_level_rows.coarse_fine");
    growth.commit();
  }
  if (workspace.coarse_rows.capacity() < coarse_need) {
    workspace.coarse_rows.reserve(coarse_need);
  }
  if (workspace.fine_rows.capacity() < fine_need) {
    workspace.fine_rows.reserve(fine_need);
  }
  fillActiveRowsForLevelInto(state, min_level, requested_rows, workspace.coarse_rows);
  fillActiveRowsForLevelInto(state, max_level, requested_rows, workspace.fine_rows);
  if (governor != nullptr) {
    // Retain admission for the actual capacity so repeated subcycles reuse it.
    // Release the transient growth lease and transfer to the workspace lease.
    growth.release();
    const std::uint64_t actual_bytes = core::checkedMemoryBytesAdd(
        static_cast<std::uint64_t>(core::checkedSizeMultiply(
            workspace.coarse_rows.capacity(), sizeof(std::uint32_t),
            "AMR active coarse actual bytes")),
        static_cast<std::uint64_t>(core::checkedSizeMultiply(
            workspace.fine_rows.capacity(), sizeof(std::uint32_t),
            "AMR active fine actual bytes")),
        "AMR active level actual bytes");
    workspace.reservation.release();
    if (actual_bytes > 0U) {
      workspace.reservation =
          governor->reserve(core::MemoryClass::kPhaseResident, actual_bytes,
                            "amr.active_level_rows.retained");
      workspace.reservation.commit();
    }
    workspace.admitted_bytes = actual_bytes;
  }
}

std::size_t mergeFluxRegistersIntoPendingStore(
    core::SimulationState& state,
    std::span<const FluxRegisterEntry> entries,
    const ProductionAmrHydroOptions& options) {
  const std::size_t max_records = core::checkedSizeAdd(
      state.pending_flux_registers.size(), entries.size(),
      "AMR pending register replacement count");
  const std::uint64_t replacement_bound = static_cast<std::uint64_t>(
      core::checkedSizeMultiply(max_records, sizeof(core::PendingFluxRegisterRecord),
                                "AMR pending register replacement bytes"));
  internal::AmrRetainedStateAdmission admission(
      state, options.regrid_memory_governor, replacement_bound,
      "amr.pending_flux_registers.replacement");
  core::PendingFluxRegisterStore candidate_store;
  std::vector<core::PendingFluxRegisterRecord> records;
  records.reserve(max_records);
  if (records.capacity() != max_records) {
    throw std::length_error("AMR pending register capacity exceeds admission");
  }
  records.insert(records.end(), state.pending_flux_registers.records().begin(),
                 state.pending_flux_registers.records().end());
  candidate_store.assign(std::move(records));
  std::size_t created_count = 0;
  auto& pending_store = candidate_store;
  for (const FluxRegisterEntry& entry : entries) {
    if (entry.register_key == 0U) {
      continue;
    }
    auto* pending = pending_store.findByRegisterKey(entry.register_key);
    if (pending == nullptr) {
      core::PendingFluxRegisterRecord seed = makePendingRecordSeed(state, entry, options);
      pending = &pending_store.upsertByRegisterKey(seed);
      ++created_count;
    } else {
      validatePendingCompatible(*pending, entry);
      pending->expected_fine_substeps = std::max(pending->expected_fine_substeps, options.expected_fine_substeps);
      pending->expected_area_comov = std::max({pending->expected_area_comov, entry.face_area_comov, entry.coarse_area_comov, entry.fine_area_comov});
    }

    if (entry.coarse_face_count > 0U && entry.coarse_area_comov > 0.0) {
      pending->coarse_face_count += entry.coarse_face_count;
      pending->coarse_area_accumulated_comov = std::max(pending->coarse_area_accumulated_comov, entry.coarse_area_comov);
      addFluxIntegral(
          pending->coarse_mass_flux_integral_code,
          pending->coarse_momentum_x_flux_integral_code,
          pending->coarse_momentum_y_flux_integral_code,
          pending->coarse_momentum_z_flux_integral_code,
          pending->coarse_total_energy_flux_integral_code,
          pending->coarse_metal_mass_flux_integral_code,
          entry.coarse_face_flux_code,
          entry.coarse_area_comov,
          entry.dt_code);
    }
    if (entry.fine_face_count > 0U && entry.fine_area_comov > 0.0) {
      pending->fine_face_count += entry.fine_face_count;
      pending->fine_area_accumulated_comov = std::max(pending->fine_area_accumulated_comov, entry.fine_area_comov);
      addFluxIntegral(
          pending->fine_mass_flux_integral_code,
          pending->fine_momentum_x_flux_integral_code,
          pending->fine_momentum_y_flux_integral_code,
          pending->fine_momentum_z_flux_integral_code,
          pending->fine_total_energy_flux_integral_code,
          pending->fine_metal_mass_flux_integral_code,
          entry.fine_face_flux_code,
          entry.fine_area_comov,
          entry.dt_code);
      const std::uint32_t expected = std::max<std::uint32_t>(1U, pending->expected_fine_substeps);
      const std::uint32_t substep = expected == 1U ? 0U : std::min(options.fine_substep_index, expected - 1U);
      if (substep < 64U) {
        const std::uint64_t bit = 1ULL << substep;
        if ((pending->fine_substep_coverage_mask & bit) == 0U) {
          pending->fine_substep_coverage_mask |= bit;
          pending->completed_fine_substeps += 1U;
        }
      } else {
        pending->completed_fine_substeps += 1U;
      }
    }
  }
  const std::uint64_t old_owned = state.pending_flux_registers.ownedCapacityBytes();
  const std::uint64_t new_owned = candidate_store.ownedCapacityBytes();
  static_assert(std::is_nothrow_move_assignable_v<core::PendingFluxRegisterStore>);
  admission.commit(old_owned, new_owned, [&]() noexcept {
    state.pending_flux_registers = std::move(candidate_store);
  });
  return created_count;
}

RefluxDiagnostics applyCompletePendingFluxRegistersToSimulationState(
    core::SimulationState& state,
    std::span<const PatchDescriptor> all_patches,
    double adiabatic_index,
    const AmrThermodynamicCachePolicy& thermodynamic_cache) {
  RefluxDiagnostics diagnostics;
  // Gate-2 atomicity: phase 1 determines complete records, phase 2 fully
  // validates every convertible target (including stale mappings and invalid
  // axis/side metadata that fluxEntryFromPendingRecord can throw on) BEFORE
  // any conserved-state mutation, phase 3 applies, phase 4 erases. Validation
  // failure leaves canonical gas state AND pending register state unchanged so
  // retry remains possible. No second SimulationState is created.
  struct PendingCandidate {
    std::size_t pending_index = 0U;
    FluxRegisterEntry entry{};
  };
  std::vector<PendingCandidate> candidates;
  candidates.reserve(state.pending_flux_registers.size());
  auto mutable_records = state.pending_flux_registers.mutableRecords();
  for (std::size_t index = 0; index < mutable_records.size(); ++index) {
    const core::PendingFluxRegisterRecord& pending = mutable_records[index];
    if (!pending.isComplete()) {
      ++diagnostics.skipped_incomplete_register_count;
      continue;
    }
    if (pending.gas_cell_identity_generation != state.gasCellIdentityGeneration()) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    // Conversion can throw on invalid axis/side metadata. Doing it here,
    // before any mutation, keeps the operation atomic.
    FluxRegisterEntry entry = fluxEntryFromPendingRecord(pending);
    candidates.push_back(PendingCandidate{.pending_index = index, .entry = entry});
  }
  // Phase 2: validate all candidate targets before mutating any canonical row.
  // Pure validation pass (no mutation): replicate the throwing stale-mapping
  // check for every candidate. This is the atomicity gate.
  for (const auto& candidate : candidates) {
    const FluxRegisterEntry& entry = candidate.entry;
    if (!entry.isComplete()) {
      continue;
    }
    const double area_scale = std::max({1.0, std::abs(entry.coarse_area_comov), std::abs(entry.fine_area_comov)});
    if (entry.coarse_area_comov <= 0.0 || entry.fine_area_comov <= 0.0 ||
        std::abs(entry.coarse_area_comov - entry.fine_area_comov) > 1.0e-10 * area_scale) {
      continue;
    }
    if (entry.coarse_gas_cell_id == 0U) {
      continue;
    }
    const auto patch_it = std::find_if(
        all_patches.begin(), all_patches.end(),
        [&entry](const PatchDescriptor& patch) { return patch.patch_id == entry.coarse_patch_id; });
    if (patch_it == all_patches.end()) {
      continue;
    }
    const auto row_opt = state.gas_cell_identity.rowForGasCellId(entry.coarse_gas_cell_id);
    if (!row_opt.has_value()) {
      continue;
    }
    const std::uint32_t row = *row_opt;
    const auto* record = state.gas_cell_identity.findByGasCellId(entry.coarse_gas_cell_id);
    if (record == nullptr || record->owning_patch_id != entry.coarse_patch_id) {
      continue;
    }
    const auto patch_index = patchIndexById(state, entry.coarse_patch_id);
    if (!patch_index.has_value() || row >= state.cells.patch_index.size() ||
        state.cells.patch_index[row] != *patch_index ||
        state.gas_cells.gas_cell_id[row] != entry.coarse_gas_cell_id ||
        entry.coarse_cell_index >= product(patch_it->cell_dims) ||
        patchLocalCellForRow(state, *patch_it, row, "applyFluxRegistersToSimulationState").linear_index != entry.coarse_cell_index) {
      throw std::runtime_error("applyFluxRegistersToSimulationState rejected a stale AMR reflux target mapping");
    }
  }
  // Phase 3: validation succeeded, so apply all corrections. Each single-entry
  // apply re-validates but cannot hit the stale-mapping throw because phase 2
  // already proved every target. Phase 4 marks/erases only truly applied rows.
  for (const auto& candidate : candidates) {
    const FluxRegisterEntry& entry = candidate.entry;
    RefluxDiagnostics one = applyFluxRegistersToSimulationState(
        state, std::span<const FluxRegisterEntry>(&entry, 1), all_patches, adiabatic_index, thermodynamic_cache);
    diagnostics.complete_register_count += one.complete_register_count;
    diagnostics.skipped_incomplete_register_count += one.skipped_incomplete_register_count;
    diagnostics.skipped_area_mismatch_count += one.skipped_area_mismatch_count;
    diagnostics.skipped_missing_target_count += one.skipped_missing_target_count;
    diagnostics.corrected_cells += one.corrected_cells;
    diagnostics.corrected_mass_code += one.corrected_mass_code;
    diagnostics.corrected_momentum_x_code += one.corrected_momentum_x_code;
    diagnostics.corrected_momentum_y_code += one.corrected_momentum_y_code;
    diagnostics.corrected_momentum_z_code += one.corrected_momentum_z_code;
    diagnostics.corrected_total_energy_code += one.corrected_total_energy_code;
    diagnostics.corrected_metal_mass_code += one.corrected_metal_mass_code;
    diagnostics.corrected_energy_code += one.corrected_energy_code;
    diagnostics.corrected_internal_energy_code += one.corrected_internal_energy_code;
    diagnostics.dual_energy_reflux_cells += one.dual_energy_reflux_cells;
    diagnostics.dual_energy_reflux_thermal_correction_code += one.dual_energy_reflux_thermal_correction_code;
    if (one.complete_register_count == 1U) {
      // Mark as consumed by zeroing the key so eraseCompleted is a no-op
      // for already-applied records. The register_key==0 sentinel is already
      // skipped by mergeFluxRegistersIntoPendingStore.
      mutable_records[candidate.pending_index].register_key = 0U;
    }
  }
  // Erase records whose key was zeroed (successfully applied).
  state.pending_flux_registers.eraseCompleted();
  return diagnostics;
}

namespace {

// The synchronized sweep owns one reusable target SoA, seven local-source
// lanes, active indexes and the complete sparse solver scratch. Admission and
// physical preparation finish before the first canonical cell update. This
// owner excludes geometry, prepared ghosts and flux-register storage, which
// have independent leases. No governor operations occur inside numerical loops.
struct AmrSynchronizedSweepWorkspace {
  core::MemoryReservation reservation;
  hydro::HydroConservedStateSoa conserved;
  std::array<std::vector<double>, 7> source_fields;
  std::vector<std::size_t> active_cells;
  std::vector<std::size_t> active_faces;
  hydro::HydroScratchBuffers scratch;
  std::uint64_t actual_capacity_bytes = 0U;

  explicit AmrSynchronizedSweepWorkspace(
      std::span<const AmrHydroPatchGeometry> geometries,
      const ProductionAmrHydroOptions& options) {
    std::size_t max_real = 0U;
    std::size_t max_total = 0U;
    std::size_t max_faces = 0U;
    for (const auto& geometry : geometries) {
      max_real = std::max(max_real, geometry.geometry.cellCount());
      max_total = std::max(max_total, geometry.geometry.totalCellStorageCount());
      max_faces = std::max(max_faces, geometry.geometry.faces.size());
    }
    const auto bytes = [](std::size_t count, std::size_t width) {
      return static_cast<std::uint64_t>(core::checkedSizeMultiply(
          count, width, "AMR synchronized sweep workspace"));
    };
    const bool dual_energy = options.active_batch_policy.dual_energy_workspace;
    std::uint64_t bound = bytes(max_total, (dual_energy ? 14U : 13U) * sizeof(double));
    bound = core::checkedMemoryBytesAdd(bound,
        bytes(max_real, sizeof(std::size_t)), "AMR active cells");
    bound = core::checkedMemoryBytesAdd(bound,
        bytes(max_faces, sizeof(std::size_t)), "AMR active faces");
    bound = core::checkedMemoryBytesAdd(bound,
        hydro::HydroScratchBuffers::workspaceBytes(max_real, max_faces,
            max_real, max_total, options.active_batch_policy), "AMR solver scratch");
    if (options.regrid_memory_governor != nullptr) {
      reservation = options.regrid_memory_governor->reserve(
          core::MemoryClass::kPhaseResident, bound, "amr.synchronized_sweep.workspace");
    }
    conserved.enableAuxiliaryInternalEnergy(dual_energy);
    conserved.resize(max_total);
    for (auto& field : source_fields) {
      field.reserve(max_total);
      if (field.capacity() != max_total) {
        throw std::length_error("AMR source field exceeded admitted capacity");
      }
    }
    active_cells.reserve(max_real);
    active_faces.reserve(max_faces);
    if (active_cells.capacity() != max_real || active_faces.capacity() != max_faces) {
      throw std::length_error("AMR active index exceeded admitted capacity");
    }
    (void)scratch.prepareForActiveSet(max_real, max_faces, max_real, max_total,
                                     options.active_batch_policy);
    actual_capacity_bytes = conserved.ownedCapacityBytes();
    for (const auto& field : source_fields) {
      actual_capacity_bytes = core::checkedMemoryBytesAdd(actual_capacity_bytes,
          bytes(field.capacity(), sizeof(double)), "AMR source physical capacity");
    }
    actual_capacity_bytes = core::checkedMemoryBytesAdd(actual_capacity_bytes,
        bytes(active_cells.capacity(), sizeof(std::size_t)), "AMR active physical capacity");
    actual_capacity_bytes = core::checkedMemoryBytesAdd(actual_capacity_bytes,
        bytes(active_faces.capacity(), sizeof(std::size_t)), "AMR active physical capacity");
    actual_capacity_bytes = core::checkedMemoryBytesAdd(actual_capacity_bytes,
        scratch.ownedCapacityBytes(), "AMR solver physical capacity");
    if (actual_capacity_bytes > bound) {
      throw std::length_error("AMR synchronized workspace exceeds admitted physical bound");
    }
    if (reservation.pending()) reservation.commit();
  }
};

[[nodiscard]] ProductionAmrHydroDiagnostics advanceProductionAmrHydroSynchronizedImpl(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options,
    const DistributedAmrHydroExchange* distributed_exchange) {
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("advanceProductionAmrHydro requires complete AMR patch coverage in SimulationState");
  }
  if (distributed_exchange != nullptr && options.persist_incomplete_flux_registers) {
    throw std::runtime_error(
        "distributed AMR pending flux-register accumulation is gated off; synchronized remote reflux must complete in-step");
  }
  ProductionAmrHydroDiagnostics diagnostics;
  const std::vector<PatchDescriptor> descriptors = buildProductionAmrPatchDescriptors(state);
  std::vector<PatchDescriptor> all_descriptors = descriptors;
  if (distributed_exchange != nullptr) {
    all_descriptors.reserve(descriptors.size() + distributed_exchange->remote_patches.size());
    for (const DistributedAmrRemotePatch& remote : distributed_exchange->remote_patches) {
      all_descriptors.push_back(remote.patch);
    }
  }
  diagnostics.patch_count = all_descriptors.size();
  const bool periodic_box = options.physical_boundary_kind == hydro::HydroBoundaryKind::kPeriodic;
  const internal::AmrPeriodicFrame periodic_frame =
      internal::makeAmrPeriodicFrame(periodic_box, options.periodic_box_extent_comoving);
  const std::uint8_t ghost_layers = amrHydroGhostLayersForOptions(options);

  const bool all_cells_active = active_cell_rows.empty();
  // Sparse active membership is a sorted, phase-local index rather than a
  // population-scale collection of hash nodes. Keep the original input and
  // patch iteration order unchanged; only membership lookup is transformed.
  const std::uint64_t active_lookup_bytes = all_cells_active ? 0U :
      core::checkedMemoryBytesAdd(
          core::checkedIntegralNarrow<std::uint64_t>(
              core::checkedSizeMultiply(active_cell_rows.size(), sizeof(std::uint32_t),
                                        "AMR active lookup bytes"),
              "AMR active lookup byte width"),
          256U, "AMR active lookup alignment");
  core::GovernedScratchArena active_lookup_arena(
      options.regrid_memory_governor, core::MemoryClass::kScratchArena,
      active_lookup_bytes, "amr.active_cell_lookup");
  std::pmr::vector<std::uint32_t> active_lookup(active_lookup_arena.resource());
  if (!all_cells_active) {
    active_lookup.assign(active_cell_rows.begin(), active_cell_rows.end());
    std::sort(active_lookup.begin(), active_lookup.end());
    active_lookup.erase(std::unique(active_lookup.begin(), active_lookup.end()),
                        active_lookup.end());
    if (!active_lookup.empty() && active_lookup.back() >= state.cells.size()) {
      throw std::out_of_range("AMR active cell row exceeds canonical cell extent");
    }
  }
  const auto isActiveCell = [&active_lookup](std::uint32_t row) {
    return std::binary_search(active_lookup.begin(), active_lookup.end(), row);
  };

  // Admit the entire retained geometry set before construction. The row-map
  // scratch has a separate lease and can be released after the last builder.
  // Other stage buffers retain their own owners; this is not a whole-stage
  // certificate and must not be added to an existing owner-held lease.
  std::uint64_t geometry_bytes = static_cast<std::uint64_t>(core::checkedSizeMultiply(
      descriptors.size(), sizeof(AmrHydroPatchGeometry), "AMR geometry objects"));
  std::uint64_t row_scratch_bytes = 0U;
  for (const PatchDescriptor& patch : descriptors) {
    const auto capacity = amrHydroGeometryCapacity(patch, ghost_layers);
    geometry_bytes = core::checkedMemoryBytesAdd(
        geometry_bytes, capacity.retained_bytes, "AMR geometry live set");
    row_scratch_bytes = std::max(row_scratch_bytes, capacity.construction_scratch_bytes);
  }
  core::MemoryReservation geometry_reservation;
  if (options.regrid_memory_governor != nullptr) {
    geometry_reservation = options.regrid_memory_governor->reserve(
        core::MemoryClass::kPhaseResident, geometry_bytes, "amr.geometry.retained");
  }
  std::vector<AmrHydroPatchGeometry> geometries;
  geometries.reserve(descriptors.size());
  if (geometries.capacity() != descriptors.size()) {
    throw std::runtime_error("AMR geometry object capacity exceeds admitted bound");
  }
  if (geometry_reservation.pending()) geometry_reservation.commit();
  core::MemoryReservation row_scratch_reservation;
  if (options.regrid_memory_governor != nullptr && row_scratch_bytes != 0U) {
    row_scratch_reservation = options.regrid_memory_governor->reserve(
        core::MemoryClass::kScratchArena, row_scratch_bytes, "amr.geometry.row_map");
    row_scratch_reservation.commit();
  }
  for (const PatchDescriptor& patch : descriptors) {
    AmrHydroGeometryOptions geometry_options;
    geometry_options.physical_boundary_kind = options.physical_boundary_kind;
    geometry_options.ghost_layers = ghost_layers;
    geometry_options.boundary_classes = {
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kLower, periodic_frame),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kUpper, periodic_frame),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kLower, periodic_frame),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kUpper, periodic_frame),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kLower, periodic_frame),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kUpper, periodic_frame),
    };
    if (periodic_box) {
      validatePeriodicAmrPatchTopology(
          patch, all_descriptors, options.periodic_box_extent_comoving, geometry_options.boundary_classes);
    }
    geometries.push_back(buildAmrHydroPatchGeometry(state, patch, geometry_options));
    populateAmrHydroFluxRegisterFaces(geometries.back(), all_descriptors,
        periodic_box ? options.periodic_box_extent_comoving : std::array<double, 3>{0.0, 0.0, 0.0});
    const std::uint64_t actual = geometries.back().ownedCapacityBytes();
    if (actual > amrHydroGeometryCapacity(patch, ghost_layers).retained_bytes) {
      throw std::runtime_error("AMR geometry physical capacity exceeds admitted topology bound");
    }
    diagnostics.geometry_capacity_bytes = core::checkedMemoryBytesAdd(
        diagnostics.geometry_capacity_bytes, actual, "AMR geometry actual capacity");
  }
  row_scratch_reservation.release();

  std::vector<AmrHydroSparseGhostSource> remote_sources;
  if (distributed_exchange != nullptr) {
    remote_sources.reserve(distributed_exchange->remote_patches.size());
    for (const DistributedAmrRemotePatch& remote : distributed_exchange->remote_patches) {
      if (remote.owner_rank == distributed_exchange->local_rank) {
        throw std::runtime_error("distributed AMR remote patch source cannot be owned by the local rank");
      }
      std::size_t previous_offset = 0U;
      bool have_previous = false;
      const std::size_t patch_cell_count =
          static_cast<std::size_t>(remote.patch.cell_dims[0]) *
          static_cast<std::size_t>(remote.patch.cell_dims[1]) *
          static_cast<std::size_t>(remote.patch.cell_dims[2]);
      for (const AmrHydroSparseRemoteCell& cell : remote.boundary_cells) {
        if (cell.patch_local_cell >= patch_cell_count || cell.gas_cell_id == 0U ||
            (have_previous && cell.patch_local_cell <= previous_offset)) {
          throw std::runtime_error(
              "distributed AMR sparse remote boundary cells are unsorted, duplicated, or out of range");
        }
        previous_offset = cell.patch_local_cell;
        have_previous = true;
      }
      remote_sources.push_back(AmrHydroSparseGhostSource{
          .patch = &remote.patch,
          .owner_rank = remote.owner_rank,
          .ghost_hydro_epoch = remote.ghost_hydro_epoch,
          .expected_ghost_hydro_epoch = remote.expected_ghost_hydro_epoch,
          .source_current_state_time_code = options.state_time_code,
          .cells = remote.boundary_cells});
    }
  }

  refreshTraceableGhostSourceIds(geometries, periodic_frame);
  if (distributed_exchange != nullptr) {
    refreshRemoteTraceableGhostSourceIds(geometries, distributed_exchange->remote_patches, periodic_frame);
  }

  const auto patchRequiresGhostFill = [&](std::size_t patch_index) {
    return all_cells_active || std::any_of(
        geometries[patch_index].real_cells.begin(),
        geometries[patch_index].real_cells.end(),
        [&isActiveCell](const AmrHydroCellDescriptor& cell) {
          return isActiveCell(cell.local_cell_row);
        });
  };
  const auto addGhostDiagnostics = [](AmrHydroGhostFillDiagnostics& target,
                                      const AmrHydroGhostFillDiagnostics& source) {
    target.physical_ghosts_filled += source.physical_ghosts_filled;
    target.same_level_ghosts_filled += source.same_level_ghosts_filled;
    target.coarse_to_fine_ghosts_filled += source.coarse_to_fine_ghosts_filled;
    target.fine_to_coarse_ghosts_filled += source.fine_to_coarse_ghosts_filled;
    target.temporal_coarse_to_fine_ghosts_filled += source.temporal_coarse_to_fine_ghosts_filled;
    target.temporal_endpoint_ghosts_filled += source.temporal_endpoint_ghosts_filled;
    target.skipped_remote_ghosts += source.skipped_remote_ghosts;
    target.stale_epoch_rejections += source.stale_epoch_rejections;
    target.temporal_same_level_mismatch_rejections += source.temporal_same_level_mismatch_rejections;
    target.temporal_history_missing_rejections += source.temporal_history_missing_rejections;
    target.temporal_history_invalid_rejections += source.temporal_history_invalid_rejections;
    target.temporal_time_out_of_range_rejections += source.temporal_time_out_of_range_rejections;
    target.temporal_fine_to_coarse_misalignment_rejections += source.temporal_fine_to_coarse_misalignment_rejections;
    target.temporal_geometry_mismatch_rejections += source.temporal_geometry_mismatch_rejections;
    target.temporal_identity_mismatch_rejections += source.temporal_identity_mismatch_rejections;
    target.missing_source_records += source.missing_source_records;
    target.unresolved_ghosts += source.unresolved_ghosts;
    target.deep_layer_extrapolated_ghosts += source.deep_layer_extrapolated_ghosts;
  };

  // Snapshot only the ghost shell needed by each active local patch while the
  // entire canonical SimulationState is still at the common step-start epoch.
  // This preserves the old synchronized-sweep semantics without retaining a
  // full conserved SoA for every local patch simultaneously.
  // One bounded, contiguous step-start snapshot replaces nested vectors.
  // Offsets are immutable after construction, so no snapshot address can be
  // invalidated by a later patch. Inactive patches occupy no shell storage.
  std::size_t prepared_ghost_count = 0U;
  for (std::size_t patch_index = 0; patch_index < geometries.size(); ++patch_index) {
    if (patchRequiresGhostFill(patch_index)) {
      prepared_ghost_count = core::checkedSizeAdd(
          prepared_ghost_count, geometries[patch_index].geometry.ghost_cells.size(),
          "AMR prepared ghost total");
    }
  }
  const std::size_t ghost_offset_count = core::checkedSizeAdd(
      geometries.size(), 1U, "AMR prepared ghost offsets");
  const std::uint64_t prepared_ghost_bytes = amrPreparedGhostWorkspaceBytes(
      geometries.size(), prepared_ghost_count);
  core::GovernedScratchArena prepared_ghost_arena(
      options.regrid_memory_governor, core::MemoryClass::kPhaseResident,
      prepared_ghost_bytes, "amr.ghosts.step_start_snapshots");
  std::pmr::vector<std::size_t> prepared_ghost_offsets(prepared_ghost_arena.resource());
  std::pmr::vector<hydro::HydroConservedState> prepared_ghost_states(prepared_ghost_arena.resource());
  prepared_ghost_offsets.reserve(ghost_offset_count);
  prepared_ghost_offsets.push_back(0U);
  for (std::size_t patch_index = 0; patch_index < geometries.size(); ++patch_index) {
    const std::size_t count = patchRequiresGhostFill(patch_index)
        ? geometries[patch_index].geometry.ghost_cells.size() : 0U;
    prepared_ghost_offsets.push_back(core::checkedSizeAdd(
        prepared_ghost_offsets.back(), count, "AMR prepared ghost offset"));
  }
  if (prepared_ghost_offsets.back() != prepared_ghost_count) {
    throw std::logic_error("AMR prepared ghost membership changed during preparation");
  }
  prepared_ghost_states.resize(prepared_ghost_count);
  diagnostics.prepared_ghost_capacity_bytes = static_cast<std::uint64_t>(
      core::checkedSizeMultiply(prepared_ghost_states.capacity(),
                                sizeof(hydro::HydroConservedState),
                                "AMR prepared ghost actual capacity"));
  for (std::size_t patch_index = 0; patch_index < geometries.size(); ++patch_index) {
    if (!patchRequiresGhostFill(patch_index)) {
      continue;
    }
    const auto touchesTarget = [&](std::size_t source_index) {
      if (source_index == patch_index) return false;
      bool touches = false;
      for (const hydro::HydroFaceAxis axis : {
               hydro::HydroFaceAxis::kX,
               hydro::HydroFaceAxis::kY,
               hydro::HydroFaceAxis::kZ}) {
        touches = touches || touchesOnSide(
            geometries[patch_index].patch, geometries[source_index].patch,
            axis, hydro::HydroFaceSide::kLower, periodic_frame);
        touches = touches || touchesOnSide(
            geometries[patch_index].patch, geometries[source_index].patch,
            axis, hydro::HydroFaceSide::kUpper, periodic_frame);
      }
      return touches;
    };
    std::size_t source_count = 0U;
    std::size_t source_cells = 0U;
    for (std::size_t source_index = 0; source_index < geometries.size(); ++source_index) {
      if (touchesTarget(source_index)) {
        source_count = core::checkedSizeAdd(source_count, 1U, "AMR ghost source count");
        source_cells = core::checkedSizeAdd(source_cells,
            geometries[source_index].geometry.totalCellStorageCount(),
            "AMR ghost source storage");
      }
    }
    const std::size_t target_cells =
        geometries[patch_index].geometry.totalCellStorageCount();
    const auto ghostBytes = [](std::size_t count, std::size_t width) {
      return static_cast<std::uint64_t>(core::checkedSizeMultiply(
          count, width, "AMR ghost preparation physical capacity"));
    };
    std::uint64_t ghost_preparation_bound = ghostBytes(
        core::checkedSizeAdd(target_cells, source_cells, "AMR ghost preparation cells"),
        6U * sizeof(double));
    ghost_preparation_bound = core::checkedMemoryBytesAdd(
        ghost_preparation_bound, ghostBytes(source_count, sizeof(std::size_t)),
        "AMR ghost source indices");
    ghost_preparation_bound = core::checkedMemoryBytesAdd(
        ghost_preparation_bound, ghostBytes(source_count, sizeof(hydro::HydroConservedStateSoa)),
        "AMR ghost source objects");
    ghost_preparation_bound = core::checkedMemoryBytesAdd(
        ghost_preparation_bound, ghostBytes(core::checkedSizeAdd(source_count, 1U,
            "AMR ghost view count"), sizeof(AmrHydroGhostFillPatch)),
        "AMR ghost preparation views");
    // This lease includes the target, all simultaneously live local source
    // SoAs and their descriptor/index capacities. It excludes the already
    // governed immutable shell and the ghost-selection scratch's own lease.
    // The lease is declared before every owned container, so failure releases
    // it only after the partially constructed physical storage is destroyed.
    core::MemoryReservation ghost_preparation_reservation;
    if (options.regrid_memory_governor != nullptr) {
      ghost_preparation_reservation = options.regrid_memory_governor->reserve(
          core::MemoryClass::kPhaseResident, ghost_preparation_bound,
          "amr.ghosts.local_source_preparation");
      ghost_preparation_reservation.commit();
    }
    hydro::HydroConservedStateSoa target_state = loadAmrHydroConservedState(
        state, geometries[patch_index], options.adiabatic_index);
    const std::uint64_t target_bytes = target_state.ownedCapacityBytes();
    diagnostics.max_patch_conserved_bytes = std::max(
        diagnostics.max_patch_conserved_bytes, target_bytes);
    if (target_bytes > ghostBytes(target_cells, 6U * sizeof(double))) {
      throw std::length_error("AMR ghost target exceeds admitted physical capacity");
    }
    std::vector<std::size_t> local_source_indices;
    local_source_indices.reserve(source_count);
    for (std::size_t source_index = 0; source_index < geometries.size(); ++source_index) {
      if (touchesTarget(source_index)) local_source_indices.push_back(source_index);
    }
    if (local_source_indices.capacity() != source_count ||
        local_source_indices.size() != source_count) {
      throw std::length_error("AMR ghost source index exceeds admitted capacity");
    }

    std::vector<hydro::HydroConservedStateSoa> local_source_states;
    local_source_states.reserve(source_count);
    if (local_source_states.capacity() != source_count) {
      throw std::length_error("AMR ghost source object capacity exceeds admission");
    }
    for (const std::size_t source_index : local_source_indices) {
      local_source_states.push_back(loadAmrHydroConservedState(
          state, geometries[source_index], options.adiabatic_index));
      const std::uint64_t source_bytes = local_source_states.back().ownedCapacityBytes();
      if (source_bytes > ghostBytes(
              geometries[source_index].geometry.totalCellStorageCount(), 6U * sizeof(double))) {
        throw std::length_error("AMR ghost source exceeds admitted physical capacity");
      }
      diagnostics.max_patch_conserved_bytes = std::max(
          diagnostics.max_patch_conserved_bytes, source_bytes);
    }

    std::vector<AmrHydroGhostFillPatch> ghost_views;
    ghost_views.reserve(core::checkedSizeAdd(source_count, 1U, "AMR ghost view capacity"));
    if (ghost_views.capacity() != source_count + 1U) {
      throw std::length_error("AMR ghost view capacity exceeds admission");
    }
    ghost_views.push_back(AmrHydroGhostFillPatch{
        .geometry = &geometries[patch_index],
        .conserved = &target_state,
        .residency = AmrGhostSourceResidency::kLocal,
        .ghost_hydro_epoch = state.gasCellIdentityGeneration(),
        .expected_ghost_hydro_epoch = state.gasCellIdentityGeneration(),
        .target_state_time_code = options.state_time_code,
        .ghost_fill_time_code = options.ghost_fill_time_code,
        .source_current_state_time_code = options.state_time_code,
        .temporal_boundary_history = options.enable_temporal_coarse_to_fine
            ? &state.amr_temporal_boundary_history
            : nullptr,
        .enable_temporal_coarse_to_fine = options.enable_temporal_coarse_to_fine,
        .requires_ghost_fill = true,
        .periodic_box_extent_comoving = periodic_box
            ? options.periodic_box_extent_comoving : std::array<double, 3>{0.0, 0.0, 0.0}});
    for (std::size_t slot = 0; slot < local_source_indices.size(); ++slot) {
      ghost_views.push_back(AmrHydroGhostFillPatch{
          .geometry = &geometries[local_source_indices[slot]],
          .conserved = &local_source_states[slot],
          .residency = AmrGhostSourceResidency::kLocal,
          .ghost_hydro_epoch = state.gasCellIdentityGeneration(),
          .expected_ghost_hydro_epoch = state.gasCellIdentityGeneration(),
          .target_state_time_code = options.state_time_code,
          .ghost_fill_time_code = options.ghost_fill_time_code,
          .source_current_state_time_code = options.state_time_code,
          .temporal_boundary_history = options.enable_temporal_coarse_to_fine
              ? &state.amr_temporal_boundary_history
              : nullptr,
          .enable_temporal_coarse_to_fine = options.enable_temporal_coarse_to_fine,
          .requires_ghost_fill = false});
    }
    addGhostDiagnostics(
        diagnostics.ghost_fill,
        fillAmrHydroGhostCells(ghost_views, remote_sources, options.adiabatic_index,
                               options.regrid_memory_governor));

    const std::size_t prepared_begin = prepared_ghost_offsets[patch_index];
    std::size_t ghost_slot = 0U;
    for (const hydro::HydroGhostCell& ghost : geometries[patch_index].geometry.ghost_cells) {
      prepared_ghost_states[prepared_begin + ghost_slot++] =
          target_state.loadCell(ghost.ghost_cell);
    }
    if (prepared_begin + ghost_slot != prepared_ghost_offsets[patch_index + 1U]) {
      throw std::logic_error("AMR prepared ghost extent differs from admitted geometry");
    }
  }

  // All legal reflux keys are already present in the immutable geometry.
  // Admit a bounded flat lookup and accumulator before the first cell update.
  // The collection scratch is released after construction, while the fixed
  // accumulator remains live through reflux. No hash-node growth occurs in
  // the numerical sweep.
  std::size_t flux_face_count = 0U;
  for (const auto& geometry : geometries) {
    for (const auto& face : geometry.geometry.flux_register_faces) {
      if (face.role != hydro::HydroFluxRegisterFaceRole::kNone) {
        flux_face_count = core::checkedSizeAdd(flux_face_count, 1U,
                                               "AMR flux-register face count");
      }
    }
  }
  const std::uint64_t key_scratch_bytes = static_cast<std::uint64_t>(
      core::checkedSizeMultiply(flux_face_count, sizeof(std::uint64_t),
                                "AMR flux-register key collection"));
  core::GovernedScratchArena flux_arena(
      options.regrid_memory_governor, core::MemoryClass::kPhaseResident,
      core::checkedMemoryBytesAdd(
          FluxRegisterAccumulator::fixedWorkspaceBytes(flux_face_count),
          key_scratch_bytes, "AMR flux-register preparation coexistence"),
      "amr.flux_registers.fixed_workspace");
  FluxRegisterAccumulator flux_registers(flux_arena.resource());
  {
    std::pmr::vector<std::uint64_t> keys(flux_arena.resource());
    keys.reserve(flux_face_count);
    for (const auto& geometry : geometries) {
      for (const auto& face : geometry.geometry.flux_register_faces) {
        if (face.role != hydro::HydroFluxRegisterFaceRole::kNone) {
          keys.push_back(face.register_key);
        }
      }
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    flux_registers.prepareFixedKeys(keys);
  }
  AmrSynchronizedSweepWorkspace sweep_workspace(geometries, options);
  core::GovernedScratchArena parent_mirror_arena(
      options.regrid_memory_governor, core::MemoryClass::kScratchArena,
      AmrParentMirrorIndex::workspaceBytes(state), "amr.parent_mirror.index");
  AmrParentMirrorIndex parent_mirror_index(parent_mirror_arena.resource());
  parent_mirror_index.build(state);
  for (std::size_t patch_index = 0; patch_index < geometries.size(); ++patch_index) {
    const AmrHydroPatchGeometry& patch_geometry = geometries[patch_index];
    auto& active_patch_cells = sweep_workspace.active_cells;
    active_patch_cells.clear();
    active_patch_cells.reserve(patch_geometry.real_cells.size());
    for (const AmrHydroCellDescriptor& cell : patch_geometry.real_cells) {
      if (all_cells_active || isActiveCell(cell.local_cell_row)) {
        active_patch_cells.push_back(cell.patch_local_cell);
      }
    }
    if (active_patch_cells.empty()) {
      continue;
    }
    auto& conserved_state = sweep_workspace.conserved;
    loadAmrHydroConservedStateInto(state, patch_geometry,
                                  options.adiabatic_index, conserved_state);
    const std::size_t prepared_begin = prepared_ghost_offsets[patch_index];
    if (prepared_ghost_offsets[patch_index + 1U] - prepared_begin !=
        patch_geometry.geometry.ghost_cells.size()) {
      throw std::logic_error("AMR hydro active patch is missing its step-start ghost snapshot");
    }
    for (std::size_t ghost_slot = 0; ghost_slot < patch_geometry.geometry.ghost_cells.size(); ++ghost_slot) {
      conserved_state.storeCell(
          patch_geometry.geometry.ghost_cells[ghost_slot].ghost_cell,
          prepared_ghost_states[prepared_begin + ghost_slot]);
    }
    if (conserved_state.auxiliaryInternalEnergyEnabled()) {
      // Dual energy: auxiliary thermal lane from the canonical e of real
      // cells; physical ghosts copy their source; same-level / coarse-fine
      // ghosts keep NaN and fall back to E - K (counted by the solver).
      std::span<double> auxiliary = conserved_state.auxiliaryInternalEnergyDensityComoving();
      std::fill(auxiliary.begin(), auxiliary.end(), std::numeric_limits<double>::quiet_NaN());
      for (std::size_t patch_cell = 0; patch_cell < patch_geometry.real_cells.size(); ++patch_cell) {
        const std::uint32_t row = patch_geometry.real_cells[patch_cell].local_cell_row;
        const double canonical_e = state.gas_cells.internal_energy_code[row];
        if (canonical_e > 0.0) {
          auxiliary[patch_cell] = conserved_state.loadCell(patch_cell).mass_density_comoving * canonical_e;
        }
      }
      hydro::fillHydroBoundaryGhostAuxiliaryInternalEnergy(conserved_state, patch_geometry.geometry);
      // Same-level / coarse-fine ghosts: rho_ghost * e_source from the
      // canonical e of the local source cell. Remote or unknown sources stay
      // NaN; their face states then use E - K and are counted, so degraded
      // dual-energy protection at such boundaries is reported, not silent.
      for (const AmrHydroGhostDescriptor& ghost : patch_geometry.ghosts) {
        if (ghost.ghost_cell >= auxiliary.size() || std::isfinite(auxiliary[ghost.ghost_cell])) {
          continue;
        }
        const auto source_row = ghost.source_gas_cell_id != 0U
            ? state.rowForGasCellId(ghost.source_gas_cell_id)
            : std::optional<std::uint32_t>{};
        const double source_e = source_row.has_value() ? state.gas_cells.internal_energy_code[*source_row] : 0.0;
        if (source_e > 0.0) {
          auxiliary[ghost.ghost_cell] = conserved_state.loadCell(ghost.ghost_cell).mass_density_comoving * source_e;
        } else {
          ++diagnostics.dual_energy_unresolved_ghosts;
        }
      }
    }
    // buildAmrHydroPatchGeometry emits real cells in increasing patch-local
    // order. This checked invariant permits allocation-free membership queries
    // without reordering the active set passed to the numerical solver.
    if (!std::is_sorted(active_patch_cells.begin(), active_patch_cells.end())) {
      throw std::logic_error("AMR patch-local active cells are not ordered");
    }
    const auto patchCellActive = [&active_patch_cells](std::size_t row) {
      return std::binary_search(active_patch_cells.begin(), active_patch_cells.end(), row);
    };
    auto& active_faces = sweep_workspace.active_faces;
    active_faces.clear();
    for (std::size_t face_index = 0; face_index < patch_geometry.geometry.faces.size(); ++face_index) {
      const hydro::HydroFace& face = patch_geometry.geometry.faces[face_index];
      std::size_t neighbor_for_activity = face.neighbor_cell;
      if (face.ghost_cell_slot != hydro::k_invalid_ghost_cell_slot) {
        const hydro::HydroGhostCell& ghost = patch_geometry.geometry.ghost_cells[face.ghost_cell_slot];
        if (ghost.boundary_kind == hydro::HydroBoundaryKind::kPeriodic) {
          // In-patch periodic wrap: the lower and upper ghost faces describe
          // the same physical box face. Evaluate it once (upper copy) and
          // let the wrapped source cell drive activity.
          if (ghost.side == hydro::HydroFaceSide::kLower) {
            continue;
          }
          neighbor_for_activity = ghost.source_real_cell;
        }
      }
      if (patchCellActive(face.owner_cell) ||
          (neighbor_for_activity < patch_geometry.geometry.cellCount() && patchCellActive(neighbor_for_activity))) {
        active_faces.push_back(face_index);
      }
    }

    auto& gravity_x = sweep_workspace.source_fields[0];
    gravity_x.assign(conserved_state.size(), 0.0);
    auto& gravity_y = sweep_workspace.source_fields[1];
    gravity_y.assign(conserved_state.size(), 0.0);
    auto& gravity_z = sweep_workspace.source_fields[2];
    gravity_z.assign(conserved_state.size(), 0.0);
    auto& mass_density_physical_cgs = sweep_workspace.source_fields[3];
    mass_density_physical_cgs.assign(conserved_state.size(), 0.0);
    auto& hydrogen = sweep_workspace.source_fields[4];
    hydrogen.assign(conserved_state.size(), 0.0);
    auto& metallicity = sweep_workspace.source_fields[5];
    metallicity.assign(conserved_state.size(), 0.0);
    auto& temperature = sweep_workspace.source_fields[6];
    temperature.assign(conserved_state.size(), 0.0);
    constexpr double k_hydrogen_mass_fraction = 0.76;
    constexpr double k_proton_mass_g = 1.67262192369e-24;
    for (std::size_t local_cell = 0; local_cell < patch_geometry.real_cells.size(); ++local_cell) {
      const std::uint32_t row = patch_geometry.real_cells[local_cell].local_cell_row;
      if (row < global_source_context.gravity_accel_x_peculiar.size()) {
        gravity_x[local_cell] = global_source_context.gravity_accel_x_peculiar[row];
      }
      if (row < global_source_context.gravity_accel_y_peculiar.size()) {
        gravity_y[local_cell] = global_source_context.gravity_accel_y_peculiar[row];
      }
      if (row < global_source_context.gravity_accel_z_peculiar.size()) {
        gravity_z[local_cell] = global_source_context.gravity_accel_z_peculiar[row];
      }
      if (row < global_source_context.mass_density_physical_cgs.size()) {
        mass_density_physical_cgs[local_cell] =
            global_source_context.mass_density_physical_cgs[row];
      } else if (options.physical_density_cgs_per_density_code > 0.0) {
        mass_density_physical_cgs[local_cell] =
            std::max(state.gas_cells.density_code[row], options.density_floor) *
            options.physical_density_cgs_per_density_code;
      }
      if (row < global_source_context.hydrogen_number_density_cgs.size()) {
        hydrogen[local_cell] = global_source_context.hydrogen_number_density_cgs[row];
      } else if (mass_density_physical_cgs[local_cell] > 0.0) {
        hydrogen[local_cell] = k_hydrogen_mass_fraction *
            mass_density_physical_cgs[local_cell] / k_proton_mass_g;
      }
      if (options.source_density_rescale != 1.0) {
        mass_density_physical_cgs[local_cell] *= options.source_density_rescale;
        hydrogen[local_cell] *= options.source_density_rescale;
      }
      if (row < global_source_context.metallicity_mass_fraction.size()) {
        metallicity[local_cell] = global_source_context.metallicity_mass_fraction[row];
      } else {
        const double gas_mass_code = state.cells.mass_code[row];
        metallicity[local_cell] = gas_mass_code > 0.0
            ? std::clamp(state.gas_cells.metal_mass_code[row] / gas_mass_code, 0.0, 1.0)
            : 0.0;
      }
      if (row < global_source_context.temperature_k.size()) {
        temperature[local_cell] = global_source_context.temperature_k[row];
      } else {
        temperature[local_cell] = state.gas_cells.temperature_code[row];
      }
    }
    const hydro::HydroSourceContext patch_source_context{
        .update = update,
        .gravity_accel_x_peculiar = gravity_x,
        .gravity_accel_y_peculiar = gravity_y,
        .gravity_accel_z_peculiar = gravity_z,
        .mass_density_physical_cgs = mass_density_physical_cgs,
        .hydrogen_number_density_cgs = hydrogen,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = temperature,
        .thermodynamic_closure = global_source_context.thermodynamic_closure,
        .redshift = global_source_context.redshift,
        .face_mass_flux_sink = global_source_context.face_mass_flux_sink};

    const auto widths = cellWidths(patch_geometry.patch);
    // Audit F1 follow-up: predictor Courant factors are dt/(a dx), matching
    // the conservative update dt/(a V) sum A F at the same epoch.
    const std::array<double, 3> courant_factors = hydro::hydroTransportCourantFactors(
        update.dt_code,
        update.scale_factor,
        std::max(widths[0], 1.0e-12),
        std::max(widths[1], 1.0e-12),
        std::max(widths[2], 1.0e-12));
    const hydro::HydroReconstructionPolicy reconstruction_policy{
        .limiter = options.slope_limiter,
        .dt_over_dx_code = courant_factors[0],
        .dt_over_cell_width_code = courant_factors,
        .rho_floor = options.density_floor > 0.0 ? options.density_floor : 1.0e-30,
        .pressure_floor = options.pressure_floor > 0.0 ? options.pressure_floor : 1.0e-30,
        .enable_muscl_hancock_predictor = true,
        .adiabatic_index = options.adiabatic_index};
    hydro::PiecewiseConstantReconstruction pcm_reconstruction;
    hydro::MusclHancockReconstruction plm_reconstruction(reconstruction_policy);
    hydro::PpmReconstruction ppm_reconstruction(reconstruction_policy);
    const hydro::HydroReconstruction& reconstruction =
        options.reconstruction_scheme == hydro::HydroReconstructionScheme::kPiecewiseConstant
            ? static_cast<const hydro::HydroReconstruction&>(pcm_reconstruction)
            : (options.reconstruction_scheme == hydro::HydroReconstructionScheme::kPiecewiseParabolic
                   ? static_cast<const hydro::HydroReconstruction&>(ppm_reconstruction)
                   : static_cast<const hydro::HydroReconstruction&>(plm_reconstruction));
    auto& scratch = sweep_workspace.scratch;
    const std::uint64_t scratch_capacity_before = scratch.ownedCapacityBytes();
    hydro::HydroProfileEvent profile;
    // CHUI Part 2 composition transport: resolve this patch's real rows and
    // ghost sources into stable local identities for the face mass sink.
    std::vector<std::uint32_t> composition_dense_rows;
    std::vector<std::uint64_t> composition_ghost_origins;
    if (global_source_context.face_mass_flux_sink != nullptr) {
      composition_dense_rows.resize(patch_geometry.geometry.cellCount(), core::kInvalidGasCellRow);
      for (const AmrHydroCellDescriptor& cell : patch_geometry.real_cells) {
        if (cell.patch_local_cell < composition_dense_rows.size()) {
          composition_dense_rows[cell.patch_local_cell] = cell.local_cell_row;
        }
      }
      composition_ghost_origins.assign(patch_geometry.geometry.ghost_cells.size(), 0U);
      for (const AmrHydroGhostDescriptor& ghost : patch_geometry.ghosts) {
        if (ghost.ghost_slot < composition_ghost_origins.size()) {
          composition_ghost_origins[ghost.ghost_slot] = ghost.source_gas_cell_id;
        }
      }
      global_source_context.face_mass_flux_sink->beginAdvance(
          patch_geometry.geometry, composition_dense_rows, composition_ghost_origins);
    }
    struct CompositionSinkScope {
      hydro::HydroFaceMassFluxSink* sink;
      ~CompositionSinkScope() { if (sink != nullptr) sink->endAdvance(); }
    } composition_sink_scope{global_source_context.face_mass_flux_sink};
    solver.advancePatchActiveSetBatchedWithScratch(
        conserved_state,
        patch_geometry.geometry,
        hydro::HydroActiveSetView{.active_cells = active_patch_cells, .active_faces = active_faces},
        update,
        reconstruction,
        riemann_solver,
        source_terms,
        patch_source_context,
        scratch,
        nullptr,
        options.active_batch_policy,
        &profile,
        &flux_registers);
    if (scratch.ownedCapacityBytes() != scratch_capacity_before) {
      throw std::logic_error("AMR solver scratch grew after physical admission: before=" +
          std::to_string(scratch_capacity_before) + " after=" +
          std::to_string(scratch.ownedCapacityBytes()) + " touched=" +
          std::to_string(scratch.touched_cells.capacity()) + " delta=" +
          std::to_string(scratch.cell_delta.capacity()) + " full=" +
          std::to_string(scratch.full_active_cells.capacity()) + " stencil=" +
          std::to_string(scratch.primitive_stencil_cells.capacity()) + " faces=" +
          std::to_string(active_faces.size()));
    }
    scatterAmrHydroConservedStateWithIndex(
        state, patch_geometry, conserved_state, options.adiabatic_index,
        parent_mirror_index, options.thermodynamic_cache);
    diagnostics.advanced_patch_count += 1U;
    diagnostics.riemann_solve_count += profile.riemann_solve_count;
    diagnostics.positivity_fallback_count += profile.positivity_fallback_count;
    diagnostics.stencil_degradation_count += profile.stencil_degradation_count;
    diagnostics.ctu_transverse_degraded_states += profile.ctu_transverse_degraded_states;
    diagnostics.dual_energy_auxiliary_cells += profile.dual_energy_auxiliary_cells;
    diagnostics.dual_energy_thermal_correction_code += profile.conservation.dual_energy_sync_delta.total_energy;
    diagnostics.active_cell_count += active_patch_cells.size();
    diagnostics.active_face_count += active_faces.size();
    diagnostics.scratch_high_water_bytes = std::max(
        diagnostics.scratch_high_water_bytes, scratch.high_water.capacity_bytes);
    diagnostics.active_batch_capacity_cells = std::max(
        diagnostics.active_batch_capacity_cells, scratch.high_water.active_cell_batch_capacity);
    diagnostics.face_batch_capacity = std::max(
        diagnostics.face_batch_capacity, scratch.high_water.face_batch_capacity);
  }
  const std::vector<FluxRegisterEntry> entries = flux_registers.entries();
  diagnostics.flux_register_entry_count = entries.size();
  std::vector<FluxRegisterEntry> local_entries;
  std::vector<FluxRegisterEntry> remote_entries;
  if (distributed_exchange != nullptr) {
    local_entries.reserve(entries.size());
    remote_entries.reserve(entries.size());
    for (const FluxRegisterEntry& entry : entries) {
      if (entry.coarse_gas_cell_id != 0U && state.rowForGasCellId(entry.coarse_gas_cell_id).has_value()) {
        local_entries.push_back(entry);
      } else {
        remote_entries.push_back(entry);
      }
    }
    if (distributed_exchange->outbound_remote_flux_registers != nullptr) {
      distributed_exchange->outbound_remote_flux_registers->insert(
          distributed_exchange->outbound_remote_flux_registers->end(),
          remote_entries.begin(),
          remote_entries.end());
    } else if (!remote_entries.empty()) {
      throw std::runtime_error("distributed AMR generated remote reflux entries without an exchange sink");
    }
  } else {
    local_entries = entries;
  }
  if (options.persist_incomplete_flux_registers) {
    diagnostics.pending_register_created_count = mergeFluxRegistersIntoPendingStore(state, local_entries, options);
    std::size_t complete_pending = 0;
    for (const core::PendingFluxRegisterRecord& pending : state.pending_flux_registers.records()) {
      if (pending.isComplete()) {
        ++complete_pending;
      }
    }
    diagnostics.pending_register_completed_count = complete_pending;
    diagnostics.pending_register_deferred_count = state.pending_flux_registers.size() - complete_pending;
    diagnostics.reflux = applyCompletePendingFluxRegistersToSimulationState(
        state, descriptors, options.adiabatic_index, options.thermodynamic_cache);
    diagnostics.pending_register_applied_count = diagnostics.reflux.complete_register_count;
    diagnostics.pending_register_rejected_count = diagnostics.reflux.skipped_area_mismatch_count +
        diagnostics.reflux.skipped_missing_target_count;
  } else {
    diagnostics.reflux = applyFluxRegistersToSimulationState(
        state, local_entries, descriptors, options.adiabatic_index, options.thermodynamic_cache);
  }
  return diagnostics;
}

}  // namespace

ProductionAmrHydroDiagnostics advanceProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options) {
  if (options.sweep_mode == ProductionAmrHydroSweepMode::kLocalSubcycling) {
    return advanceProductionAmrHydroSubcycled(
        state,
        active_cell_rows,
        update,
        global_source_context,
        solver,
        riemann_solver,
        source_terms,
        options);
  }
  return advanceProductionAmrHydroSynchronizedImpl(
      state,
      active_cell_rows,
      update,
      global_source_context,
      solver,
      riemann_solver,
      source_terms,
      options,
      nullptr);
}

ProductionAmrHydroDiagnostics advanceDistributedProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options,
    const DistributedAmrHydroExchange& exchange) {
  if (options.sweep_mode == ProductionAmrHydroSweepMode::kLocalSubcycling ||
      options.enable_temporal_coarse_to_fine) {
    throw std::runtime_error(
        "distributed AMR hydro currently supports synchronized patch execution only; "
        "MPI temporal subcycling is not implemented");
  }
  return advanceProductionAmrHydroSynchronizedImpl(
      state,
      active_cell_rows,
      update,
      global_source_context,
      solver,
      riemann_solver,
      source_terms,
      options,
      &exchange);
}

ProductionAmrHydroDiagnostics advanceProductionAmrHydroSubcycled(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& coarse_update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options) {
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("advanceProductionAmrHydroSubcycled requires complete AMR patch coverage in SimulationState");
  }
  if (options.refinement_ratio < 2U) {
    throw std::invalid_argument("AMR hydro subcycling requires refinement_ratio >= 2");
  }
  if (!std::isfinite(options.state_time_code)) {
    throw std::invalid_argument("AMR hydro subcycling requires finite state_time_code");
  }
  const std::vector<PatchDescriptor> descriptors = buildProductionAmrPatchDescriptors(state);
  int min_level = std::numeric_limits<int>::max();
  int max_level = std::numeric_limits<int>::min();
  for (const PatchDescriptor& patch : descriptors) {
    const int patch_level = static_cast<int>(patch.level);
    min_level = std::min(min_level, patch_level);
    max_level = std::max(max_level, patch_level);
  }
  if (min_level == std::numeric_limits<int>::max()) {
    return ProductionAmrHydroDiagnostics{};
  }
  // The temporal-history representation below records one coarse interval and
  // is intentionally scoped to a single local coarse/fine pair.  Rejecting a
  // deeper hierarchy is safer than silently reusing an interval from the wrong
  // parent fine substep; recursive multi-level temporal histories remain work
  // for a later scheduler-integrated implementation.
  if (max_level - min_level > 1) {
    throw std::runtime_error(
        "local temporal AMR subcycling currently supports one coarse/fine level pair; "
        "a deeper hierarchy requires nested temporal boundary histories");
  }
  if (!state.amr_temporal_boundary_history.empty()) {
    throw std::runtime_error(
        "cannot start a new local AMR subcycle while restart-resumable temporal boundary history is active");
  }

  ProductionAmrHydroDiagnostics diagnostics;
  diagnostics.patch_count = descriptors.size();
  diagnostics.levels_advanced = static_cast<std::size_t>(max_level - min_level + 1);
  diagnostics.max_level_advanced = static_cast<std::size_t>(std::max(0, max_level));
  diagnostics.substeps_by_level.assign(static_cast<std::size_t>(max_level + 1), 0U);

  // Gate-2B: coarse_rows and fine_rows coexist. Admit BOTH capacities via the
  // existing regrid governor BEFORE growth, reuse capacity across the fine
  // substep loop (allocated once here), and reconcile actual capacity().
  // Tight-headroom rejection happens here, before any solver mutation below.
  AmrActiveLevelRowWorkspace level_rows;
  prepareAmrActiveLevelRows(
      state, active_cell_rows, min_level, max_level,
      options.regrid_memory_governor, level_rows);
  const std::vector<std::uint32_t>& coarse_rows = level_rows.coarse_rows;
  const std::vector<std::uint32_t>& fine_rows = level_rows.fine_rows;
  const double coarse_start_code = options.state_time_code;
  const double coarse_end_code = coarse_start_code + coarse_update.dt_code;
  captureAmrTemporalBoundaryHistoryStart(state, descriptors, coarse_start_code, options.adiabatic_index, options.regrid_memory_governor);

  if (!coarse_rows.empty()) {
    ProductionAmrHydroOptions coarse_options = options;
    coarse_options.sweep_mode = ProductionAmrHydroSweepMode::kSynchronized;
    coarse_options.persist_incomplete_flux_registers = true;
    coarse_options.reflux_interval_start_code = coarse_start_code;
    coarse_options.reflux_interval_end_code = coarse_end_code;
    coarse_options.reflux_coarse_dt_code = coarse_update.dt_code;
    coarse_options.expected_fine_substeps = options.refinement_ratio;
    coarse_options.fine_substep_index = 0U;
    coarse_options.state_time_code = coarse_start_code;
    coarse_options.ghost_fill_time_code = coarse_start_code;
    coarse_options.enable_temporal_coarse_to_fine = false;
    const ProductionAmrHydroDiagnostics coarse_diag = advanceProductionAmrHydro(
        state, coarse_rows, coarse_update, global_source_context, solver, riemann_solver, source_terms, coarse_options);
    accumulateDiagnostics(diagnostics, coarse_diag);
    diagnostics.substeps_by_level[static_cast<std::size_t>(min_level)] = 1U;
    diagnostics.subcycled_level_step_count += 1U;
  }
  captureAmrTemporalBoundaryHistoryEnd(state, descriptors, coarse_end_code, options.adiabatic_index, options.regrid_memory_governor);

  const double fine_dt_code = coarse_update.dt_code / static_cast<double>(options.refinement_ratio);
  for (std::uint32_t substep = 0; substep < options.refinement_ratio; ++substep) {
    if (fine_rows.empty()) {
      continue;
    }
    const double fine_start_code = coarse_start_code + static_cast<double>(substep) * fine_dt_code;
    hydro::HydroUpdateContext fine_update = coarse_update;
    fine_update.dt_code = fine_dt_code;
    fine_update.time_begin_code = fine_start_code;
    fine_update.time_end_code = fine_start_code + fine_dt_code;
    const double fraction_begin = static_cast<double>(substep) /
        static_cast<double>(options.refinement_ratio);
    const double fraction_end = static_cast<double>(substep + 1U) /
        static_cast<double>(options.refinement_ratio);
    const double fraction_mid = 0.5 * (fraction_begin + fraction_end);
    const auto interpolate = [](double begin, double end, double fraction) {
      return begin + fraction * (end - begin);
    };
    fine_update.scale_factor_begin = interpolate(
        coarse_update.scale_factor_begin, coarse_update.scale_factor_end, fraction_begin);
    fine_update.scale_factor_end = interpolate(
        coarse_update.scale_factor_begin, coarse_update.scale_factor_end, fraction_end);
    fine_update.scale_factor = interpolate(
        coarse_update.scale_factor_begin, coarse_update.scale_factor_end, fraction_mid);
    fine_update.hubble_rate_begin_code = interpolate(
        coarse_update.hubble_rate_begin_code, coarse_update.hubble_rate_end_code, fraction_begin);
    fine_update.hubble_rate_end_code = interpolate(
        coarse_update.hubble_rate_begin_code, coarse_update.hubble_rate_end_code, fraction_end);
    fine_update.hubble_rate_code = interpolate(
        coarse_update.hubble_rate_begin_code, coarse_update.hubble_rate_end_code, fraction_mid);

    hydro::HydroSourceContext fine_source_context = global_source_context;
    fine_source_context.update = fine_update;
    fine_source_context.redshift = fine_update.scale_factor > 0.0
        ? std::max(0.0, 1.0 / fine_update.scale_factor - 1.0)
        : global_source_context.redshift;
    double density_rescale = 1.0;
    if (fine_update.comoving_coordinates &&
        coarse_update.scale_factor > 0.0 && fine_update.scale_factor > 0.0) {
      density_rescale = std::pow(
          coarse_update.scale_factor / fine_update.scale_factor, 3.0);
    }
    ProductionAmrHydroOptions fine_options = options;
    fine_options.source_density_rescale = options.source_density_rescale * density_rescale;
    fine_options.sweep_mode = ProductionAmrHydroSweepMode::kSynchronized;
    fine_options.persist_incomplete_flux_registers = true;
    fine_options.reflux_interval_start_code = coarse_start_code;
    fine_options.reflux_interval_end_code = coarse_end_code;
    fine_options.reflux_coarse_dt_code = coarse_update.dt_code;
    fine_options.expected_fine_substeps = options.refinement_ratio;
    fine_options.fine_substep_index = substep;
    fine_options.state_time_code = fine_start_code;
    // Ghosts are consumed before the predictor inside HydroCoreSolver, hence
    // this is the physical beginning of the fine substep.
    fine_options.ghost_fill_time_code = fine_start_code;
    fine_options.enable_temporal_coarse_to_fine = true;
    const ProductionAmrHydroDiagnostics fine_diag = advanceProductionAmrHydro(
        state, fine_rows, fine_update, fine_source_context, solver, riemann_solver, source_terms, fine_options);
    accumulateDiagnostics(diagnostics, fine_diag);
    diagnostics.substeps_by_level[static_cast<std::size_t>(max_level)] += 1U;
    diagnostics.subcycled_level_step_count += 1U;
  }

  const RefluxDiagnostics post_child_reflux = applyCompletePendingFluxRegistersToSimulationState(
      state, buildProductionAmrPatchDescriptors(state), options.adiabatic_index, options.thermodynamic_cache);
  diagnostics.reflux.complete_register_count += post_child_reflux.complete_register_count;
  diagnostics.reflux.skipped_incomplete_register_count += post_child_reflux.skipped_incomplete_register_count;
  diagnostics.reflux.skipped_area_mismatch_count += post_child_reflux.skipped_area_mismatch_count;
  diagnostics.reflux.skipped_missing_target_count += post_child_reflux.skipped_missing_target_count;
  diagnostics.reflux.corrected_cells += post_child_reflux.corrected_cells;
  diagnostics.reflux.corrected_mass_code += post_child_reflux.corrected_mass_code;
  diagnostics.reflux.corrected_momentum_x_code += post_child_reflux.corrected_momentum_x_code;
  diagnostics.reflux.corrected_momentum_y_code += post_child_reflux.corrected_momentum_y_code;
  diagnostics.reflux.corrected_momentum_z_code += post_child_reflux.corrected_momentum_z_code;
  diagnostics.reflux.corrected_total_energy_code += post_child_reflux.corrected_total_energy_code;
  diagnostics.reflux.corrected_metal_mass_code += post_child_reflux.corrected_metal_mass_code;
  diagnostics.reflux.corrected_energy_code += post_child_reflux.corrected_energy_code;
  diagnostics.reflux.corrected_internal_energy_code += post_child_reflux.corrected_internal_energy_code;
  diagnostics.reflux.dual_energy_reflux_cells += post_child_reflux.dual_energy_reflux_cells;
  diagnostics.reflux.dual_energy_reflux_thermal_correction_code +=
      post_child_reflux.dual_energy_reflux_thermal_correction_code;
  diagnostics.pending_register_applied_count += post_child_reflux.complete_register_count;
  diagnostics.pending_register_rejected_count += post_child_reflux.skipped_area_mismatch_count +
      post_child_reflux.skipped_missing_target_count;

  std::size_t complete_pending = 0;
  for (const core::PendingFluxRegisterRecord& pending : state.pending_flux_registers.records()) {
    if (pending.isComplete()) {
      ++complete_pending;
    }
  }
  diagnostics.pending_register_completed_count = std::max(diagnostics.pending_register_completed_count, complete_pending);
  diagnostics.pending_register_deferred_count = state.pending_flux_registers.size() - complete_pending;
  if (state.pending_flux_registers.empty()) {
    retireAmrTemporalBoundaryHistory(state);
  }
  return diagnostics;
}

ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t first_child_patch_id,
    std::uint64_t first_child_gas_cell_id,
    const ProductionAmrHydroOptions& options) {
  requireProductionRegridSynchronizationQuiescent(state, "production AMR refine");
  if (first_child_patch_id == 0U || first_child_gas_cell_id == 0U) {
    throw std::invalid_argument("production AMR refine requires nonzero child patch and gas-cell ID seeds");
  }
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("production AMR refine requires complete SimulationState patch coverage");
  }
  const auto parent_patch_index = patchIndexById(state, parent_patch.patch_id);
  if (!parent_patch_index.has_value()) {
    throw std::runtime_error("production AMR refine target patch is absent from SimulationState");
  }
  if (parent_patch.level == std::numeric_limits<std::uint8_t>::max()) {
    throw std::overflow_error("production AMR refine child level overflows uint8");
  }
  if (parent_patch.morton_key >
      (std::numeric_limits<std::uint64_t>::max() - 8U) / 8U) {
    throw std::overflow_error("production AMR refine child Morton key overflows uint64");
  }

  const std::size_t parent_cell_count = product(parent_patch.cell_dims);
  if (parent_cell_count == 0U || state.cells.size() < parent_cell_count || state.patches.size() == 0U) {
    throw std::runtime_error("production AMR refine has invalid parent-cell/patch population");
  }
  const std::size_t child_cell_count = core::checkedSizeMultiply(
      parent_cell_count, 8U, "production AMR refine child cell count");
  const std::size_t new_cell_count = core::checkedSizeAdd(
      state.cells.size() - parent_cell_count, child_cell_count,
      "production AMR refine new cell count");
  const std::size_t new_patch_count = core::checkedSizeAdd(
      state.patches.size() - 1U, 8U, "production AMR refine new patch count");
  // Gate-3 live-set: parent_rows (parent_cell_count u32) and new_child_rows
  // (child_cell_count u32) coexist with the candidate hierarchy (new_cell /
  // new_patch state) before commit. Admission includes BOTH row buffers plus
  // candidate state; unrelated phase maxima are NOT summed.
  const std::size_t refine_row_scratch_bytes = core::checkedSizeAdd(
      core::checkedSizeMultiply(
          parent_cell_count, sizeof(std::uint32_t),
          "production AMR refine parent-row scratch bytes"),
      core::checkedSizeMultiply(
          child_cell_count, sizeof(std::uint32_t),
          "production AMR refine diagnostic row scratch bytes"),
      "production AMR refine row scratch bytes");
  ProductionAmrRegridAdmission admission = admitProductionAmrRegrid(
      state, options,
      productionAmrCandidateReservationBytes(new_cell_count, new_patch_count, refine_row_scratch_bytes),
      "amr.production_refine.transaction");

  const std::vector<std::uint32_t> parent_rows =
      patchLocalRowsForPatch(state, parent_patch, "production AMR refine");
  if (parent_rows.size() != parent_cell_count) {
    throw std::runtime_error("production AMR refine parent row coverage does not match descriptor");
  }
  validatePatchIdRangeAvailable(state, first_child_patch_id, 8U);
  validateGasCellIdRangeAvailable(state, first_child_gas_cell_id, child_cell_count);
  const double parent_volume =
      parent_patch.extent_comov[0] * parent_patch.extent_comov[1] *
      parent_patch.extent_comov[2] / static_cast<double>(parent_cell_count);
  const auto before = totalsForRows(state, parent_rows, parent_volume, options);
  const std::uint64_t identity_generation = nextIdentityGeneration(state);

  core::SimulationState candidate;
  candidate.resizeCells(new_cell_count);
  candidate.resizePatches(new_patch_count);
  const std::array<double, 3> child_extent = {
      parent_patch.extent_comov[0] * 0.5,
      parent_patch.extent_comov[1] * 0.5,
      parent_patch.extent_comov[2] * 0.5};
  const double child_volume = parent_volume / 8.0;

  std::size_t write_row = 0U;
  std::size_t write_patch = 0U;
  std::vector<std::uint32_t> new_child_rows;
  new_child_rows.reserve(child_cell_count);
  std::uint64_t next_child_gas_cell_id = first_child_gas_cell_id;
  for (std::size_t old_patch_index = 0; old_patch_index < state.patches.size(); ++old_patch_index) {
    if (state.patches.patch_id[old_patch_index] != parent_patch.patch_id) {
      const std::uint32_t first = state.patches.first_cell[old_patch_index];
      const std::uint32_t count = state.patches.cell_count[old_patch_index];
      writePatchDescriptorToStateRow(candidate, write_patch, patchDescriptorFromStateRow(state, old_patch_index));
      candidate.patches.first_cell[write_patch] = core::checkedIntegralNarrow<std::uint32_t>(
          write_row, "production AMR refine kept patch first row");
      candidate.patches.cell_count[write_patch] = count;
      candidate.patches.owning_rank[write_patch] = state.patches.owning_rank[old_patch_index];
      for (std::uint32_t offset = 0; offset < count; ++offset) {
        const std::uint32_t dst = core::checkedIntegralNarrow<std::uint32_t>(
            write_row, "production AMR refine kept cell destination row");
        copyCellRow(candidate, dst, state, first + offset);
        candidate.cells.patch_index[dst] = core::checkedIntegralNarrow<std::uint32_t>(
            write_patch, "production AMR refine kept cell patch index");
        ++write_row;
      }
      ++write_patch;
      continue;
    }

    for (std::uint8_t octant = 0; octant < 8U; ++octant) {
      const std::array<double, 3> child_origin = {
          parent_patch.origin_comov[0] + (((octant & 1U) != 0U) ? child_extent[0] : 0.0),
          parent_patch.origin_comov[1] + (((octant & 2U) != 0U) ? child_extent[1] : 0.0),
          parent_patch.origin_comov[2] + (((octant & 4U) != 0U) ? child_extent[2] : 0.0)};
      const std::array<double, 3> child_widths = {
          child_extent[0] / static_cast<double>(parent_patch.cell_dims[0]),
          child_extent[1] / static_cast<double>(parent_patch.cell_dims[1]),
          child_extent[2] / static_cast<double>(parent_patch.cell_dims[2])};
      PatchDescriptor child_patch = parent_patch;
      child_patch.patch_id = first_child_patch_id + static_cast<std::uint64_t>(octant);
      child_patch.parent_patch_id = parent_patch.patch_id;
      child_patch.level = static_cast<std::uint8_t>(parent_patch.level + 1U);
      child_patch.morton_key = parent_patch.morton_key * 8U +
          static_cast<std::uint64_t>(octant) + 1U;
      child_patch.origin_comov = child_origin;
      child_patch.extent_comov = child_extent;
      writePatchDescriptorToStateRow(candidate, write_patch, child_patch);
      candidate.patches.first_cell[write_patch] = core::checkedIntegralNarrow<std::uint32_t>(
          write_row, "production AMR refine child patch first row");
      candidate.patches.cell_count[write_patch] = core::checkedIntegralNarrow<std::uint32_t>(
          parent_cell_count, "production AMR refine child patch cell count");
      candidate.patches.owning_rank[write_patch] = state.patches.owning_rank[*parent_patch_index];
      for (std::size_t child_cell = 0; child_cell < parent_cell_count; ++child_cell) {
        const std::size_t i = child_cell % parent_patch.cell_dims[0];
        const std::size_t j =
            (child_cell / parent_patch.cell_dims[0]) % parent_patch.cell_dims[1];
        const std::size_t xy = core::checkedSizeMultiply(
            static_cast<std::size_t>(parent_patch.cell_dims[0]),
            static_cast<std::size_t>(parent_patch.cell_dims[1]),
            "production AMR refine parent xy cell count");
        const std::size_t k = child_cell / xy;
        const std::array<double, 3> center = {
            child_origin[0] + (static_cast<double>(i) + 0.5) * child_widths[0],
            child_origin[1] + (static_cast<double>(j) + 0.5) * child_widths[1],
            child_origin[2] + (static_cast<double>(k) + 0.5) * child_widths[2]};
        const std::size_t parent_cell = cellContainingPoint(parent_patch, center);
        const std::uint32_t parent_row = parent_rows[parent_cell];
        const ConservedState parent_conserved =
            volumeIntegratedForRow(state, parent_row, parent_volume, options);
        const ConservedState child_conserved = parent_conserved * (1.0 / 8.0);
        const std::uint32_t dst = core::checkedIntegralNarrow<std::uint32_t>(
            write_row, "production AMR refine child destination row");
        candidate.cells.center_x_comoving[dst] = center[0];
        candidate.cells.center_y_comoving[dst] = center[1];
        candidate.cells.center_z_comoving[dst] = center[2];
        candidate.cells.time_bin[dst] = state.cells.time_bin[parent_row];
        candidate.cells.patch_index[dst] = core::checkedIntegralNarrow<std::uint32_t>(
            write_patch, "production AMR refine child patch index");
        candidate.gas_cells.gas_cell_id[dst] = next_child_gas_cell_id++;
        candidate.gas_cells.parent_particle_id[dst] =
            state.gas_cells.parent_particle_id[parent_row];
        writeVolumeIntegratedToRow(candidate, dst, child_conserved, child_volume, options);
        new_child_rows.push_back(dst);
        ++write_row;
      }
      ++write_patch;
    }
  }
  if (write_row != new_cell_count || write_patch != new_patch_count) {
    throw std::runtime_error("production AMR refine internal row accounting failed");
  }
  rebuildIdentityFromSidecars(candidate, identity_generation);
  if (!hasProductionAmrHydroCoverage(candidate) ||
      !candidate.gas_cell_identity.coversDenseLocalRows(candidate.cells.size()) ||
      !candidate.gasCellIdentityMapMatchesSidecarLanes()) {
    throw std::runtime_error("production AMR refine candidate failed authority validation");
  }
  const auto after = totalsForRows(candidate, new_child_rows, child_volume, options);
  const ProductionAmrRegridDiagnostics diagnostics{
      .refined_patch_count = 1,
      .created_gas_cell_count = new_child_rows.size(),
      .retired_gas_cell_count = parent_rows.size(),
      .conserved_mass_before = before[0],
      .conserved_mass_after = after[0],
      .conserved_momentum_x_before = before[1],
      .conserved_momentum_x_after = after[1],
      .conserved_momentum_y_before = before[2],
      .conserved_momentum_y_after = after[2],
      .conserved_momentum_z_before = before[3],
      .conserved_momentum_z_after = after[3],
      .conserved_total_energy_before = before[4],
      .conserved_total_energy_after = after[4]};
  commitProductionAmrCandidate(state, std::move(candidate), admission);
  return diagnostics;
}

ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options) {
  const std::size_t child_cell_count = product(parent_patch.cell_dims) * 8U;
  return refineProductionPatchInSimulationState(
      state, parent_patch, nextFreePatchId(state, 8U), nextFreeGasCellId(state, child_cell_count), options);
}

ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t replacement_gas_cell_id,
    const ProductionAmrHydroOptions& options) {
  requireProductionRegridSynchronizationQuiescent(state, "production AMR derefine");
  if (replacement_gas_cell_id == 0U) {
    throw std::invalid_argument("production AMR derefine requires a nonzero replacement gas_cell_id seed");
  }
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("production AMR derefine requires complete SimulationState patch coverage");
  }
  if (parent_patch.level == std::numeric_limits<std::uint8_t>::max()) {
    throw std::overflow_error("production AMR derefine child-level query overflows uint8");
  }
  const std::size_t parent_cell_count = product(parent_patch.cell_dims);
  const std::size_t child_cell_count = core::checkedSizeMultiply(
      parent_cell_count, 8U, "production AMR derefine child cell count");
  if (state.cells.size() < child_cell_count || state.patches.size() < 8U) {
    throw std::runtime_error("production AMR derefine state is too small for an exact child octet");
  }
  const std::size_t kept_cell_count = core::checkedSizeAdd(
      state.cells.size() - child_cell_count, parent_cell_count,
      "production AMR derefine kept cell count");
  const std::size_t kept_patch_count = core::checkedSizeAdd(
      state.patches.size() - 8U, 1U,
      "production AMR derefine kept patch count");
  // Gate-3 live-set: child_rows (child u32) + restricted triple (parent
  // Conserved/u64/u8) + one transient per-child row buffer (parent u32, reused
  // sequentially) coexist with the kept-cell/kept-patch candidate before
  // commit. Admission includes all simultaneously-live scratch.
  std::size_t scratch_bytes = core::checkedSizeMultiply(
      parent_cell_count,
      sizeof(ConservedState) + sizeof(std::uint64_t) + sizeof(std::uint8_t),
      "production AMR derefine restriction scratch bytes");
  scratch_bytes = core::checkedSizeAdd(
      scratch_bytes,
      core::checkedSizeMultiply(
          child_cell_count, sizeof(std::uint32_t),
          "production AMR derefine child-row scratch bytes"),
      "production AMR derefine scratch bytes");
  scratch_bytes = core::checkedSizeAdd(
      scratch_bytes,
      core::checkedSizeMultiply(
          parent_cell_count, sizeof(std::uint32_t),
          "production AMR derefine transient per-child row scratch bytes"),
      "production AMR derefine scratch bytes");
  ProductionAmrRegridAdmission admission = admitProductionAmrRegrid(
      state, options,
      productionAmrCandidateReservationBytes(kept_cell_count, kept_patch_count, scratch_bytes),
      "amr.production_derefine.transaction");

  std::vector<PatchDescriptor> children;
  children.reserve(8U);
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    const PatchDescriptor patch = patchDescriptorFromStateRow(state, patch_index);
    const bool parent_link_ok =
        patch.parent_patch_id == parent_patch.patch_id || patch.parent_patch_id == 0U;
    const bool same_parent_extent =
        parent_link_ok && patch.level == static_cast<std::uint8_t>(parent_patch.level + 1U) &&
        patch.origin_comov[0] >= parent_patch.origin_comov[0] - k_geometry_tol &&
        patch.origin_comov[1] >= parent_patch.origin_comov[1] - k_geometry_tol &&
        patch.origin_comov[2] >= parent_patch.origin_comov[2] - k_geometry_tol &&
        patch.origin_comov[0] + patch.extent_comov[0] <=
            parent_patch.origin_comov[0] + parent_patch.extent_comov[0] + k_geometry_tol &&
        patch.origin_comov[1] + patch.extent_comov[1] <=
            parent_patch.origin_comov[1] + parent_patch.extent_comov[1] + k_geometry_tol &&
        patch.origin_comov[2] + patch.extent_comov[2] <=
            parent_patch.origin_comov[2] + parent_patch.extent_comov[2] + k_geometry_tol;
    if (same_parent_extent) {
      children.push_back(patch);
    }
  }
  const std::array<double, 3> expected_child_extent{
      parent_patch.extent_comov[0] * 0.5,
      parent_patch.extent_comov[1] * 0.5,
      parent_patch.extent_comov[2] * 0.5};
  std::array<bool, 8> octant_seen{};
  for (const PatchDescriptor& child : children) {
    if (child.level != static_cast<std::uint8_t>(parent_patch.level + 1U) ||
        child.cell_dims != parent_patch.cell_dims) {
      throw std::runtime_error(
          "production AMR derefine child level/cell dimensions do not match exact octant contract");
    }
    std::uint8_t octant = 0U;
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double scale = std::max(
          {1.0, std::abs(child.extent_comov[axis]), std::abs(expected_child_extent[axis])});
      if (std::abs(child.extent_comov[axis] - expected_child_extent[axis]) >
          k_geometry_tol * scale) {
        throw std::runtime_error(
            "production AMR derefine child extent does not match exact octant contract");
      }
      const double lower = parent_patch.origin_comov[axis];
      const double upper = parent_patch.origin_comov[axis] + expected_child_extent[axis];
      const double child_origin = child.origin_comov[axis];
      const double lower_scale = std::max({1.0, std::abs(child_origin), std::abs(lower)});
      const double upper_scale = std::max({1.0, std::abs(child_origin), std::abs(upper)});
      if (std::abs(child_origin - lower) <= k_geometry_tol * lower_scale) {
        continue;
      }
      if (std::abs(child_origin - upper) <= k_geometry_tol * upper_scale) {
        octant |= static_cast<std::uint8_t>(1U << axis);
        continue;
      }
      throw std::runtime_error(
          "production AMR derefine child origin does not match an exact parent octant");
    }
    if (octant_seen[octant]) {
      throw std::runtime_error("production AMR derefine detected duplicate child octant");
    }
    octant_seen[octant] = true;
  }
  if (children.size() != 8U ||
      !std::all_of(octant_seen.begin(), octant_seen.end(), [](bool seen) { return seen; })) {
    throw std::runtime_error(
        "production AMR derefine requires exactly eight child patches covering the parent descriptor");
  }

  const double child_volume =
      children.front().extent_comov[0] * children.front().extent_comov[1] *
      children.front().extent_comov[2] /
      static_cast<double>(product(children.front().cell_dims));
  std::vector<std::uint32_t> child_rows;
  child_rows.reserve(child_cell_count);
  for (const PatchDescriptor& child : children) {
    const std::vector<std::uint32_t> rows =
        patchLocalRowsForPatch(state, child, "production AMR derefine child coverage");
    child_rows.insert(child_rows.end(), rows.begin(), rows.end());
  }
  if (child_rows.size() != child_cell_count) {
    throw std::runtime_error("production AMR derefine child row coverage does not match exact octet");
  }
  const auto before = totalsForRows(state, child_rows, child_volume, options);
  validatePatchIdRangeAvailable(state, parent_patch.patch_id, 1U);
  validateGasCellIdRangeAvailable(state, replacement_gas_cell_id, parent_cell_count);

  const double parent_volume =
      parent_patch.extent_comov[0] * parent_patch.extent_comov[1] *
      parent_patch.extent_comov[2] / static_cast<double>(parent_cell_count);
  std::vector<ConservedState> restricted(parent_cell_count);
  std::vector<std::uint64_t> restricted_parent_particle_id(
      parent_cell_count, std::numeric_limits<std::uint64_t>::max());
  std::vector<std::uint8_t> restricted_time_bin(
      parent_cell_count, std::numeric_limits<std::uint8_t>::max());
  for (const PatchDescriptor& child : children) {
    const std::vector<std::uint32_t> rows =
        patchLocalRowsForPatch(state, child, "production AMR derefine child restriction");
    for (std::size_t child_cell = 0; child_cell < rows.size(); ++child_cell) {
      const std::array<double, 3> center = cellCenter(child, child_cell);
      const std::size_t parent_cell = cellContainingPoint(parent_patch, center);
      const std::uint32_t child_row = rows[child_cell];
      restricted[parent_cell] +=
          volumeIntegratedForRow(state, child_row, child_volume, options);
      const std::uint64_t child_parent = state.gas_cells.parent_particle_id[child_row];
      std::uint64_t& parent_slot = restricted_parent_particle_id[parent_cell];
      if (parent_slot == std::numeric_limits<std::uint64_t>::max()) {
        parent_slot = child_parent;
      } else if (parent_slot != child_parent) {
        parent_slot = 0U;
      }
      restricted_time_bin[parent_cell] =
          std::min(restricted_time_bin[parent_cell], state.cells.time_bin[child_row]);
    }
  }

  // Gate-3: exactly eight child IDs; sorted adjacent comparison replaces the
  // small hash set with deterministic allocation-free duplicate detection.
  std::array<std::uint64_t, 8U> child_patch_ids{};
  for (std::size_t i = 0; i < children.size() && i < child_patch_ids.size(); ++i) {
    child_patch_ids[i] = children[i].patch_id;
  }
  {
    std::array<std::uint64_t, 8U> sorted_ids = child_patch_ids;
    std::sort(sorted_ids.begin(), sorted_ids.end());
    for (std::size_t i = 1; i < sorted_ids.size(); ++i) {
      if (sorted_ids[i - 1U] == sorted_ids[i]) {
        throw std::runtime_error("production AMR derefine detected duplicate child patch id");
      }
    }
  }
  auto isChildPatch = [&child_patch_ids](std::uint64_t patch_id) {
    return std::find(child_patch_ids.begin(), child_patch_ids.end(), patch_id) !=
        child_patch_ids.end();
  };
  const std::uint64_t identity_generation = nextIdentityGeneration(state);
  core::SimulationState candidate;
  candidate.resizeCells(kept_cell_count);
  candidate.resizePatches(kept_patch_count);
  std::size_t write_row = 0U;
  std::size_t write_patch = 0U;
  std::vector<std::uint32_t> parent_rows;
  parent_rows.reserve(parent_cell_count);
  bool inserted_parent = false;
  std::uint64_t next_parent_gas_cell_id = replacement_gas_cell_id;
  for (std::size_t old_patch_index = 0; old_patch_index < state.patches.size(); ++old_patch_index) {
    const std::uint64_t old_patch_id = state.patches.patch_id[old_patch_index];
    if (isChildPatch(old_patch_id)) {
      if (!inserted_parent) {
        writePatchDescriptorToStateRow(candidate, write_patch, parent_patch);
        candidate.patches.first_cell[write_patch] = core::checkedIntegralNarrow<std::uint32_t>(
            write_row, "production AMR derefine parent first row");
        candidate.patches.cell_count[write_patch] = core::checkedIntegralNarrow<std::uint32_t>(
            parent_cell_count, "production AMR derefine parent cell count");
        candidate.patches.owning_rank[write_patch] = state.patches.owning_rank[old_patch_index];
        for (std::size_t cell = 0; cell < parent_cell_count; ++cell) {
          const std::uint32_t dst = core::checkedIntegralNarrow<std::uint32_t>(
              write_row, "production AMR derefine parent destination row");
          const std::array<double, 3> center = cellCenter(parent_patch, cell);
          candidate.cells.center_x_comoving[dst] = center[0];
          candidate.cells.center_y_comoving[dst] = center[1];
          candidate.cells.center_z_comoving[dst] = center[2];
          candidate.cells.time_bin[dst] =
              restricted_time_bin[cell] == std::numeric_limits<std::uint8_t>::max()
              ? 0U
              : restricted_time_bin[cell];
          candidate.cells.patch_index[dst] = core::checkedIntegralNarrow<std::uint32_t>(
              write_patch, "production AMR derefine parent patch index");
          candidate.gas_cells.gas_cell_id[dst] = next_parent_gas_cell_id++;
          candidate.gas_cells.parent_particle_id[dst] =
              restricted_parent_particle_id[cell] == std::numeric_limits<std::uint64_t>::max()
              ? 0U
              : restricted_parent_particle_id[cell];
          writeVolumeIntegratedToRow(candidate, dst, restricted[cell], parent_volume, options);
          parent_rows.push_back(dst);
          ++write_row;
        }
        ++write_patch;
        inserted_parent = true;
      }
      continue;
    }
    const std::uint32_t first = state.patches.first_cell[old_patch_index];
    const std::uint32_t count = state.patches.cell_count[old_patch_index];
    writePatchDescriptorToStateRow(candidate, write_patch, patchDescriptorFromStateRow(state, old_patch_index));
    candidate.patches.first_cell[write_patch] = core::checkedIntegralNarrow<std::uint32_t>(
        write_row, "production AMR derefine kept patch first row");
    candidate.patches.cell_count[write_patch] = count;
    candidate.patches.owning_rank[write_patch] = state.patches.owning_rank[old_patch_index];
    for (std::uint32_t offset = 0; offset < count; ++offset) {
      const std::uint32_t dst = core::checkedIntegralNarrow<std::uint32_t>(
          write_row, "production AMR derefine kept cell destination row");
      copyCellRow(candidate, dst, state, first + offset);
      candidate.cells.patch_index[dst] = core::checkedIntegralNarrow<std::uint32_t>(
          write_patch, "production AMR derefine kept cell patch index");
      ++write_row;
    }
    ++write_patch;
  }
  if (!inserted_parent || write_row != kept_cell_count || write_patch != kept_patch_count) {
    throw std::runtime_error("production AMR derefine internal row accounting failed");
  }
  rebuildIdentityFromSidecars(candidate, identity_generation);
  if (!hasProductionAmrHydroCoverage(candidate) ||
      !candidate.gas_cell_identity.coversDenseLocalRows(candidate.cells.size()) ||
      !candidate.gasCellIdentityMapMatchesSidecarLanes()) {
    throw std::runtime_error("production AMR derefine candidate failed authority validation");
  }
  const auto after = totalsForRows(candidate, parent_rows, parent_volume, options);
  const ProductionAmrRegridDiagnostics diagnostics{
      .derefined_patch_count = 1,
      .created_gas_cell_count = parent_rows.size(),
      .retired_gas_cell_count = child_rows.size(),
      .conserved_mass_before = before[0],
      .conserved_mass_after = after[0],
      .conserved_momentum_x_before = before[1],
      .conserved_momentum_x_after = after[1],
      .conserved_momentum_y_before = before[2],
      .conserved_momentum_y_after = after[2],
      .conserved_momentum_z_before = before[3],
      .conserved_momentum_z_after = after[3],
      .conserved_total_energy_before = before[4],
      .conserved_total_energy_after = after[4]};
  commitProductionAmrCandidate(state, std::move(candidate), admission);
  return diagnostics;
}

ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options) {
  return derefineProductionPatchInSimulationState(
      state, parent_patch, nextFreeGasCellId(state, product(parent_patch.cell_dims)), options);
}

}  // namespace cosmosim::amr
