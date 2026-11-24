#include "cosmosim/mhd/mhd_runtime_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "cosmosim/mhd/mhd_validation.hpp"

namespace cosmosim::mhd {
namespace {

[[nodiscard]] MhdFailureContext adapterContext(std::string owner, std::string detail) {
  MhdFailureContext context{};
  context.owner = std::move(owner);
  context.detail = std::move(detail);
  return context;
}

}  // namespace

MhdNumericsPolicy mhdNumericsPolicyFromConfig(const core::MhdConfig& config) {
  MhdNumericsPolicy policy{};
  policy.reconstruction = config.reconstruction == core::MhdReconstructionConfigMode::kPlmMc
      ? MhdReconstructionMode::kPlmMc
      : MhdReconstructionMode::kPpm;
  policy.riemann = config.riemann_solver == core::MhdRiemannSolverConfigMode::kHlle
      ? MhdRiemann::Solver::kHlle
      : MhdRiemann::Solver::kHlld;
  policy.emf_scheme = MhdEmfScheme::kUctHll;
  policy.time_integrator = MhdTimeIntegrator::kSsprk2;
  policy.divergence_policy = MhdDivergencePolicy::kStrictCt;
  policy.cfl_ideal = config.cfl_ideal;
  policy.max_step_retries = config.max_step_retries;
  policy.fallback_fraction_limit = config.fallback_fraction_limit;
  policy.floor_energy_fraction_limit = config.floor_energy_fraction_limit;
  policy.deterministic = config.deterministic;
  policy.stiffness.cfl_diffusion = config.cfl_diffusion;
  policy.stiffness.cfl_hall = config.cfl_hall;
  policy.stiffness.rkl_switch_ratio = config.nonideal.rkl_switch_ratio;
  policy.stiffness.rkl_max_stages = config.nonideal.rkl_max_stages;
  policy.stiffness.hall_subcycle_cost_cap =
      config.performance.hall_subcycle_cost_cap_per_level_interval;
  policy.stiffness.rkl_stage_cost_cap =
      config.performance.rkl_stage_cost_cap_per_level_interval;
  if (config.precision != core::MhdPrecisionConfigMode::kFloat64) {
    throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                   adapterContext("mhd_config", "only float64 MHD precision is admitted"));
  }
  if (config.magnetic_representation != core::MhdMagneticRepresentation::kFaceFluxCt) {
    throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                   adapterContext("mhd_config",
                                  "only face-flux constrained transport is admitted"));
  }
  if (config.spatial_method != core::MhdSpatialMethod::kFiniteVolumeUnsplit) {
    throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                   adapterContext("mhd_config",
                                  "only finite_volume_unsplit MHD is admitted"));
  }
  if (config.riemann_fallback != core::MhdRiemannFallbackConfigMode::kHlle) {
    throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                   adapterContext("mhd_config", "only HLLE fallback is admitted"));
  }
  return policy;
}

MhdPatchBoundary mhdBoundaryFromCondition(core::BoundaryCondition condition) {
  MhdBoundaryKind kind = MhdBoundaryKind::kOpen;
  switch (condition) {
    case core::BoundaryCondition::kPeriodic:
      kind = MhdBoundaryKind::kPeriodic;
      break;
    case core::BoundaryCondition::kOpen:
      kind = MhdBoundaryKind::kOpen;
      break;
    case core::BoundaryCondition::kReflective:
      kind = MhdBoundaryKind::kReflectingConductingWall;
      break;
  }
  MhdPatchBoundary boundary{};
  boundary.lower = {kind, kind, kind};
  boundary.upper = {kind, kind, kind};
  return boundary;
}

std::size_t findMagneticFacePatch(const core::MagneticFaceSoa& storage,
                                  std::uint64_t patch_id,
                                  std::uint64_t geometry_generation) {
  for (std::size_t index = 0; index < storage.patches.size(); ++index) {
    const core::MagneticFacePatchRecord& record = storage.patches[index];
    if (record.patch_id == patch_id && record.geometry_generation == geometry_generation) {
      return index;
    }
  }
  return storage.patches.size();
}

std::size_t ensureMagneticFacePatch(core::MagneticFaceSoa& storage,
                                    const MhdAdapterPatchBinding& binding) {
  if (!binding.geometry.valid()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "invalid patch geometry"));
  }
  const auto matches = [&](const core::MagneticFacePatchRecord& record) {
    return record.patch_id == binding.patch_id &&
        record.geometry_generation == binding.geometry_generation;
  };
  for (std::size_t index = 0; index < storage.patches.size(); ++index) {
    const core::MagneticFacePatchRecord& record = storage.patches[index];
    if (!matches(record)) {
      continue;
    }
    if (record.nx != binding.geometry.nx || record.ny != binding.geometry.ny ||
        record.nz != binding.geometry.nz) {
      throw MhdError(
          MhdErrorCode::kStaleGeometry,
          adapterContext("mhd_adapter",
                         "patch dimensions changed without an explicit AMR regrid transaction"));
    }
    return index;
  }

  // Reject a stale generation for the same patch id: a regrid must commit a
  // fresh set of records atomically, never leave two generations resident.
  for (const core::MagneticFacePatchRecord& record : storage.patches) {
    if (record.patch_id == binding.patch_id) {
      throw MhdError(MhdErrorCode::kStaleGeometry,
                     adapterContext("mhd_adapter",
                                    "multiple geometry generations for one patch id"));
    }
  }

  core::MagneticFacePatchRecord record{};
  record.patch_id = binding.patch_id;
  record.geometry_generation = binding.geometry_generation;
  record.level = binding.level;
  record.nx = static_cast<std::uint32_t>(binding.geometry.nx);
  record.ny = static_cast<std::uint32_t>(binding.geometry.ny);
  record.nz = static_cast<std::uint32_t>(binding.geometry.nz);
  record.origin_x_comoving = binding.geometry.origin_x_comoving;
  record.origin_y_comoving = binding.geometry.origin_y_comoving;
  record.origin_z_comoving = binding.geometry.origin_z_comoving;
  record.cell_width_x_comoving = binding.geometry.dx_comoving;
  record.cell_width_y_comoving = binding.geometry.dy_comoving;
  record.cell_width_z_comoving = binding.geometry.dz_comoving;
  record.flux_x_offset = storage.face_flux_x_comoving.size();
  record.flux_y_offset = storage.face_flux_y_comoving.size();
  record.flux_z_offset = storage.face_flux_z_comoving.size();
  const std::size_t base_x = storage.face_flux_x_comoving.size();
  const std::size_t base_y = storage.face_flux_y_comoving.size();
  const std::size_t base_z = storage.face_flux_z_comoving.size();
  storage.patches.push_back(record);
  storage.face_flux_x_comoving.resize(base_x + binding.geometry.faceCountX(), 0.0);
  storage.face_flux_y_comoving.resize(base_y + binding.geometry.faceCountY(), 0.0);
  storage.face_flux_z_comoving.resize(base_z + binding.geometry.faceCountZ(), 0.0);
  return storage.patches.size() - 1U;
}

MhdFaceFluxView makeMagneticFaceView(const core::MagneticFaceSoa& storage,
                                     const MhdAdapterPatchBinding& binding,
                                     std::size_t patch_index) {
  if (patch_index >= storage.patches.size()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "face-flux patch index out of range"));
  }
  const core::MagneticFacePatchRecord& record = storage.patches[patch_index];
  MhdFaceFluxView view{};
  view.geometry = &binding.geometry;
  view.phi_x = std::span<const double>(storage.face_flux_x_comoving.data() + record.flux_x_offset,
                                       binding.geometry.faceCountX());
  view.phi_y = std::span<const double>(storage.face_flux_y_comoving.data() + record.flux_y_offset,
                                       binding.geometry.faceCountY());
  view.phi_z = std::span<const double>(storage.face_flux_z_comoving.data() + record.flux_z_offset,
                                       binding.geometry.faceCountZ());
  view.patch_id = record.patch_id;
  view.geometry_generation = record.geometry_generation;
  view.refinement_level = record.level;
  return view;
}

MhdFaceFluxMutableView makeMutableMagneticFaceView(core::MagneticFaceSoa& storage,
                                                   const MhdAdapterPatchBinding& binding,
                                                   std::size_t patch_index) {
  if (patch_index >= storage.patches.size()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "face-flux patch index out of range"));
  }
  const core::MagneticFacePatchRecord& record = storage.patches[patch_index];
  MhdFaceFluxMutableView view{};
  view.geometry = &binding.geometry;
  view.phi_x = std::span<double>(storage.face_flux_x_comoving.data() + record.flux_x_offset,
                                 binding.geometry.faceCountX());
  view.phi_y = std::span<double>(storage.face_flux_y_comoving.data() + record.flux_y_offset,
                                 binding.geometry.faceCountY());
  view.phi_z = std::span<double>(storage.face_flux_z_comoving.data() + record.flux_z_offset,
                                 binding.geometry.faceCountZ());
  view.patch_id = record.patch_id;
  view.geometry_generation = record.geometry_generation;
  view.refinement_level = record.level;
  return view;
}

void loadMhdGasFromState(const core::SimulationState& state,
                         const MhdPatchGeometry& geometry,
                         std::size_t first_geometry_row,
                         std::span<const std::uint32_t> dense_row_by_geometry_row,
                         MhdGasStateSoa& gas) {
  const std::size_t cell_count = geometry.cellCount();
  gas.resize(cell_count);
  if (first_geometry_row > dense_row_by_geometry_row.size() ||
      cell_count > dense_row_by_geometry_row.size() - first_geometry_row) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "patch cell range exceeds the dense row map"));
  }
  const double volume = geometry.cellVolume();
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    const std::uint32_t dense_row = dense_row_by_geometry_row[first_geometry_row + cell];
    const double rho = state.gas_cells.density_code[dense_row];
    const double vx = state.gas_cells.velocity_x_peculiar[dense_row];
    const double vy = state.gas_cells.velocity_y_peculiar[dense_row];
    const double vz = state.gas_cells.velocity_z_peculiar[dense_row];
    const double internal = state.gas_cells.internal_energy_code[dense_row];
    const double metal_mass = state.gas_cells.metal_mass_code[dense_row];
    MhdConservedState conserved{};
    conserved.rho_comoving = rho;
    conserved.momentum_density_x_comoving = rho * vx;
    conserved.momentum_density_y_comoving = rho * vy;
    conserved.momentum_density_z_comoving = rho * vz;
    // Magnetic energy is added by the caller after the field view is bound;
    // recoverPatchPrimitives always re-derives B from the authoritative faces,
    // so the magnetic term here is only the initial total-energy bookkeeping
    // contribution. The operator subtracts the same face-derived term.
    conserved.total_energy_density_comoving = rho * internal +
        0.5 * rho * (vx * vx + vy * vy + vz * vz);
    conserved.metal_mass_density_comoving = volume > 0.0 ? metal_mass / volume : 0.0;
    gas.store(cell, conserved);
  }
}

void storeMhdGasToState(core::SimulationState& state,
                        const MhdPatchGeometry& geometry,
                        std::size_t first_geometry_row,
                        std::span<const std::uint32_t> geometry_row_by_dense_row,
                        const MhdGasStateSoa& gas,
                        double density_floor_code,
                        double pressure_floor_code) {
  const std::size_t cell_count = geometry.cellCount();
  if (gas.size() != cell_count) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "gas state size does not match the patch"));
  }
  const double volume = geometry.cellVolume();
  const double gamma = k_mhd_gamma_default;
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    const std::size_t geometry_row = first_geometry_row + cell;
    if (geometry_row >= geometry_row_by_dense_row.size()) {
      throw MhdError(MhdErrorCode::kInvalidField,
                     adapterContext("mhd_adapter", "patch cell range exceeds the dense row map"));
    }
    const std::uint32_t dense_row = geometry_row_by_dense_row[geometry_row];
    MhdConservedState conserved = gas.load(cell);
    const double rho = std::max(conserved.rho_comoving, density_floor_code);
    const double vx = conserved.momentum_density_x_comoving / std::max(rho, 1.0e-30);
    const double vy = conserved.momentum_density_y_comoving / std::max(rho, 1.0e-30);
    const double vz = conserved.momentum_density_z_comoving / std::max(rho, 1.0e-30);
    const double kinetic =
        0.5 * rho * (vx * vx + vy * vy + vz * vz);
    const double magnetic = 0.5 * (conserved.bx_comoving * conserved.bx_comoving +
                                   conserved.by_comoving * conserved.by_comoving +
                                   conserved.bz_comoving * conserved.bz_comoving);
    double internal_density =
        conserved.total_energy_density_comoving - kinetic - magnetic;
    const double pressure_floor = std::max(pressure_floor_code, 0.0);
    if (!(internal_density > 0.0)) {
      internal_density = pressure_floor / (gamma - 1.0);
    }
    const double internal_specific = internal_density / std::max(rho, 1.0e-30);
    state.gas_cells.density_code[dense_row] = rho;
    state.gas_cells.velocity_x_peculiar[dense_row] = vx;
    state.gas_cells.velocity_y_peculiar[dense_row] = vy;
    state.gas_cells.velocity_z_peculiar[dense_row] = vz;
    state.gas_cells.internal_energy_code[dense_row] = internal_specific;
    state.gas_cells.pressure_code[dense_row] =
        std::max((gamma - 1.0) * rho * internal_specific, pressure_floor);
    state.gas_cells.metal_mass_code[dense_row] =
        conserved.metal_mass_density_comoving * volume;
    state.cells.mass_code[dense_row] = rho * volume;
  }
}

bool seedUniformComovingField(core::MagneticFaceSoa& storage,
                              std::size_t patch_index,
                              double bx_comoving,
                              double by_comoving,
                              double bz_comoving) {
  if (patch_index >= storage.patches.size()) {
    return false;
  }
  const core::MagneticFacePatchRecord& record = storage.patches[patch_index];
  if (record.nx == 0U || record.ny == 0U || record.nz == 0U) {
    return false;
  }
  const std::uint64_t nx = record.nx;
  const std::uint64_t ny = record.ny;
  const std::uint64_t nz = record.nz;
  const double area_x = record.cell_width_y_comoving * record.cell_width_z_comoving;
  const double area_y = record.cell_width_x_comoving * record.cell_width_z_comoving;
  const double area_z = record.cell_width_x_comoving * record.cell_width_y_comoving;
  std::fill_n(storage.face_flux_x_comoving.begin() + record.flux_x_offset,
              (nx + 1U) * ny * nz, bx_comoving * area_x);
  std::fill_n(storage.face_flux_y_comoving.begin() + record.flux_y_offset,
              nx * (ny + 1U) * nz, by_comoving * area_y);
  std::fill_n(storage.face_flux_z_comoving.begin() + record.flux_z_offset,
              nx * ny * (nz + 1U), bz_comoving * area_z);
  return true;
}

double patchMagneticEnergy(const core::MagneticFaceSoa& storage, std::size_t patch_index) {
  if (patch_index >= storage.patches.size()) {
    return 0.0;
  }
  const core::MagneticFacePatchRecord& record = storage.patches[patch_index];
  const std::uint64_t nx = record.nx;
  const std::uint64_t ny = record.ny;
  const std::uint64_t nz = record.nz;
  const double volume =
      record.cell_width_x_comoving * record.cell_width_y_comoving *
      record.cell_width_z_comoving;
  const double area_x = record.cell_width_y_comoving * record.cell_width_z_comoving;
  const double area_y = record.cell_width_x_comoving * record.cell_width_z_comoving;
  const double area_z = record.cell_width_x_comoving * record.cell_width_y_comoving;
  if (volume <= 0.0 || area_x <= 0.0 || area_y <= 0.0 || area_z <= 0.0) {
    return 0.0;
  }
  const double* phi_x = storage.face_flux_x_comoving.data() + record.flux_x_offset;
  const double* phi_y = storage.face_flux_y_comoving.data() + record.flux_y_offset;
  const double* phi_z = storage.face_flux_z_comoving.data() + record.flux_z_offset;
  const auto index_x = [&](std::uint64_t i, std::uint64_t j, std::uint64_t k) {
    return i + (nx + 1U) * (j + ny * k);
  };
  const auto index_y = [&](std::uint64_t i, std::uint64_t j, std::uint64_t k) {
    return i + nx * (j + (ny + 1U) * k);
  };
  const auto index_z = [&](std::uint64_t i, std::uint64_t j, std::uint64_t k) {
    return i + nx * (j + ny * k);
  };
  double energy = 0.0;
  for (std::uint64_t k = 0; k < nz; ++k) {
    for (std::uint64_t j = 0; j < ny; ++j) {
      for (std::uint64_t i = 0; i < nx; ++i) {
        const double bx = 0.5 * (phi_x[index_x(i, j, k)] + phi_x[index_x(i + 1U, j, k)]) / area_x;
        const double by = 0.5 * (phi_y[index_y(i, j, k)] + phi_y[index_y(i, j + 1U, k)]) / area_y;
        const double bz = 0.5 * (phi_z[index_z(i, j, k)] + phi_z[index_z(i, j, k + 1U)]) / area_z;
        energy += 0.5 * (bx * bx + by * by + bz * bz) * volume;
      }
    }
  }
  return energy;
}

MhdConfiguredCoefficientProvider::MhdConfiguredCoefficientProvider(
    const core::MhdConfig& config)
    : m_constant(config.nonideal.eta_ohmic_phys_cgs,
                 config.nonideal.eta_hall_phys_cgs,
                 config.nonideal.eta_ambipolar_phys_cgs) {
  if (config.nonideal.coefficient_provider != core::MhdCoefficientProviderConfigMode::kConstant) {
    throw MhdError(
        MhdErrorCode::kBadCoefficient,
        adapterContext("mhd_coefficient_provider",
                       "only the constant cgs provider is instantiated in this source pass; "
                       "table and chemistry providers fail closed"));
  }
}

std::string MhdConfiguredCoefficientProvider::provenanceId() const {
  return m_constant.provenanceId();
}

void MhdConfiguredCoefficientProvider::validateAdmission() const {
  m_constant.validateAdmission();
}

MhdTransportCoefficients MhdConfiguredCoefficientProvider::evaluate(
    std::size_t cell_index,
    const MhdCellPhysicalProperties& properties) const {
  return m_constant.evaluate(cell_index, properties);
}

}  // namespace cosmosim::mhd
