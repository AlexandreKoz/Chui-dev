#include "cosmosim/mhd/mhd_runtime_adapter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
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

[[nodiscard]] std::uint64_t checkedMultiply(std::uint64_t lhs, std::uint64_t rhs) {
  if (lhs != 0U && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "face count multiplication overflow"));
  }
  return lhs * rhs;
}

[[nodiscard]] std::uint64_t checkedAdd(std::uint64_t lhs, std::uint64_t rhs) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "face offset addition overflow"));
  }
  return lhs + rhs;
}

struct FaceCounts {
  std::uint64_t x = 0U;
  std::uint64_t y = 0U;
  std::uint64_t z = 0U;
};

[[nodiscard]] FaceCounts checkedFaceCounts(std::uint64_t nx, std::uint64_t ny, std::uint64_t nz) {
  return {checkedMultiply(checkedMultiply(checkedAdd(nx, 1U), ny), nz),
          checkedMultiply(checkedMultiply(nx, checkedAdd(ny, 1U)), nz),
          checkedMultiply(checkedMultiply(nx, ny), checkedAdd(nz, 1U))};
}

[[nodiscard]] std::array<double, 3> cellFieldFromView(const MhdPatchGeometry& geometry,
                                                      const MhdFaceFluxView& faces,
                                                      std::size_t cell) {
  const std::size_t i = cell % geometry.nx;
  const std::size_t j = (cell / geometry.nx) % geometry.ny;
  const std::size_t k = cell / (geometry.nx * geometry.ny);
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

void requireRowMap(const core::SimulationState& state, const MhdPatchGeometry& geometry,
                   std::span<const std::uint32_t> rows) {
  if (rows.size() != geometry.cellCount()) {
    throw MhdError(MhdErrorCode::kStaleGeometry,
                   adapterContext("mhd_adapter", "patch cell -> gas row map has the wrong size"));
  }
  const std::size_t gas_rows = state.gas_cells.density_code.size();
  std::set<std::uint32_t> seen;
  for (const std::uint32_t row : rows) {
    if (row >= gas_rows || !seen.insert(row).second) {
      throw MhdError(MhdErrorCode::kStaleGeometry,
                     adapterContext("mhd_adapter",
                                    "patch cell -> gas row map has an out-of-range or duplicate "
                                    "row (stale identity binding)"));
    }
  }
}

}  // namespace

MhdNumericsPolicy mhdNumericsPolicyFromConfig(const core::MhdConfig& config) {
  if (config.precision != core::MhdPrecisionConfigMode::kFloat64) {
    throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                   adapterContext("mhd_config", "only float64 MHD precision is admitted"));
  }
  if (config.magnetic_representation != core::MhdMagneticRepresentation::kFaceFluxCt ||
      config.spatial_method != core::MhdSpatialMethod::kFiniteVolumeUnsplit ||
      config.riemann_fallback != core::MhdRiemannFallbackConfigMode::kHlle ||
      config.emf_scheme != core::MhdEmfSchemeConfigMode::kUctHll ||
      config.time_integrator != core::MhdTimeIntegratorConfigMode::kSsprk2 ||
      config.divergence_policy != core::MhdDivergencePolicyConfigMode::kStrictCt) {
    throw MhdError(MhdErrorCode::kUnsupportedPhysics,
                   adapterContext("mhd_config",
                                  "only the CHUI-MHD-REF-v1 face-flux CT / FV unsplit / SSPRK2 / "
                                  "UCT-HLL / strict-CT / HLLE-fallback contract is admitted"));
  }
  MhdNumericsPolicy policy{};
  policy.reconstruction = config.reconstruction == core::MhdReconstructionConfigMode::kPlmMc
      ? MhdReconstructionMode::kPlmMc
      : MhdReconstructionMode::kPpm;
  policy.riemann = config.riemann_solver == core::MhdRiemannSolverConfigMode::kHlle
      ? MhdRiemann::Solver::kHlle
      : MhdRiemann::Solver::kHlld;
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
  policy.stiffness.rkl_stage_cost_cap = config.performance.rkl_stage_cost_cap_per_level_interval;
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

core::MagneticFaceSoa buildMagneticFaceLayout(std::span<const MhdAdapterPatchBinding> bindings,
                                              std::string coefficient_table_sha256) {
  core::MagneticFaceSoa storage{};
  storage.coefficient_table_sha256 = std::move(coefficient_table_sha256);
  storage.patches.reserve(bindings.size());
  std::set<std::uint64_t> ids;
  FaceCounts total{};
  for (const MhdAdapterPatchBinding& binding : bindings) {
    const MhdPatchGeometry& geometry = binding.geometry;
    if (!geometry.valid() || geometry.nx > std::numeric_limits<std::uint32_t>::max() ||
        geometry.ny > std::numeric_limits<std::uint32_t>::max() ||
        geometry.nz > std::numeric_limits<std::uint32_t>::max() ||
        !std::isfinite(geometry.origin_x_comoving) || !std::isfinite(geometry.origin_y_comoving) ||
        !std::isfinite(geometry.origin_z_comoving)) {
      throw MhdError(MhdErrorCode::kInvalidField,
                     adapterContext("mhd_adapter", "invalid patch geometry in face layout"));
    }
    if (!ids.insert(binding.patch_id).second) {
      throw MhdError(MhdErrorCode::kStaleGeometry,
                     adapterContext("mhd_adapter",
                                    "duplicate patch_id in face layout (two generations of one "
                                    "patch cannot be resident)"));
    }
    const FaceCounts counts = checkedFaceCounts(geometry.nx, geometry.ny, geometry.nz);
    core::MagneticFacePatchRecord record{};
    record.patch_id = binding.patch_id;
    record.geometry_generation = binding.geometry_generation;
    record.level = binding.level;
    record.nx = static_cast<std::uint32_t>(geometry.nx);
    record.ny = static_cast<std::uint32_t>(geometry.ny);
    record.nz = static_cast<std::uint32_t>(geometry.nz);
    record.origin_x_comoving = geometry.origin_x_comoving;
    record.origin_y_comoving = geometry.origin_y_comoving;
    record.origin_z_comoving = geometry.origin_z_comoving;
    record.cell_width_x_comoving = geometry.dx_comoving;
    record.cell_width_y_comoving = geometry.dy_comoving;
    record.cell_width_z_comoving = geometry.dz_comoving;
    record.flux_x_offset = total.x;
    record.flux_y_offset = total.y;
    record.flux_z_offset = total.z;
    total.x = checkedAdd(total.x, counts.x);
    total.y = checkedAdd(total.y, counts.y);
    total.z = checkedAdd(total.z, counts.z);
    storage.patches.push_back(record);
  }
  if (total.x > std::numeric_limits<std::size_t>::max() ||
      total.y > std::numeric_limits<std::size_t>::max() ||
      total.z > std::numeric_limits<std::size_t>::max()) {
    throw MhdError(MhdErrorCode::kMemoryBudgetExceeded,
                   adapterContext("mhd_adapter", "face layout exceeds the addressable size"));
  }
  // One sizing per lane (no incremental growth).
  storage.face_flux_x_comoving.assign(static_cast<std::size_t>(total.x), 0.0);
  storage.face_flux_y_comoving.assign(static_cast<std::size_t>(total.y), 0.0);
  storage.face_flux_z_comoving.assign(static_cast<std::size_t>(total.z), 0.0);
  return storage;
}

void validateMagneticFaceSoa(const core::MagneticFaceSoa& storage) {
  if (!storage.isConsistent()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter",
                                  "MagneticFaceSoa records/offsets do not tile the face lanes"));
  }
  std::set<std::uint64_t> ids;
  for (const core::MagneticFacePatchRecord& record : storage.patches) {
    if (!ids.insert(record.patch_id).second) {
      throw MhdError(MhdErrorCode::kStaleGeometry,
                     adapterContext("mhd_adapter", "duplicate patch_id in MagneticFaceSoa"));
    }
    (void)checkedFaceCounts(record.nx, record.ny, record.nz);
  }
  const auto finite = [](const auto& lane) {
    return std::all_of(lane.begin(), lane.end(), [](double v) { return std::isfinite(v); });
  };
  if (!finite(storage.face_flux_x_comoving) || !finite(storage.face_flux_y_comoving) ||
      !finite(storage.face_flux_z_comoving)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "non-finite authoritative face flux"));
  }
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

MhdPatchGeometry magneticPatchGeometry(const core::MagneticFacePatchRecord& record) {
  MhdPatchGeometry geometry{};
  geometry.nx = record.nx;
  geometry.ny = record.ny;
  geometry.nz = record.nz;
  geometry.origin_x_comoving = record.origin_x_comoving;
  geometry.origin_y_comoving = record.origin_y_comoving;
  geometry.origin_z_comoving = record.origin_z_comoving;
  geometry.dx_comoving = record.cell_width_x_comoving;
  geometry.dy_comoving = record.cell_width_y_comoving;
  geometry.dz_comoving = record.cell_width_z_comoving;
  return geometry;
}

namespace {

void requireBindingMatchesRecord(const core::MagneticFaceSoa& storage,
                                 const MhdAdapterPatchBinding& binding,
                                 std::size_t patch_index) {
  if (patch_index >= storage.patches.size()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "face-flux patch index out of range"));
  }
  const core::MagneticFacePatchRecord& record = storage.patches[patch_index];
  if (record.patch_id != binding.patch_id ||
      record.geometry_generation != binding.geometry_generation ||
      record.nx != binding.geometry.nx || record.ny != binding.geometry.ny ||
      record.nz != binding.geometry.nz) {
    throw MhdError(MhdErrorCode::kStaleGeometry,
                   adapterContext("mhd_adapter",
                                  "binding (patch_id, generation, shape) does not match the face "
                                  "record; rebind after regrid/migration"));
  }
  const FaceCounts counts = checkedFaceCounts(record.nx, record.ny, record.nz);
  if (checkedAdd(record.flux_x_offset, counts.x) > storage.face_flux_x_comoving.size() ||
      checkedAdd(record.flux_y_offset, counts.y) > storage.face_flux_y_comoving.size() ||
      checkedAdd(record.flux_z_offset, counts.z) > storage.face_flux_z_comoving.size()) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "face record offsets exceed the face lanes"));
  }
}

}  // namespace

MhdFaceFluxView makeMagneticFaceView(const core::MagneticFaceSoa& storage,
                                     const MhdAdapterPatchBinding& binding,
                                     std::size_t patch_index) {
  requireBindingMatchesRecord(storage, binding, patch_index);
  const core::MagneticFacePatchRecord& record = storage.patches[patch_index];
  MhdFaceFluxView view{};
  view.geometry = &binding.geometry;
  view.phi_x = std::span<const double>(
      storage.face_flux_x_comoving.data() + record.flux_x_offset, binding.geometry.faceCountX());
  view.phi_y = std::span<const double>(
      storage.face_flux_y_comoving.data() + record.flux_y_offset, binding.geometry.faceCountY());
  view.phi_z = std::span<const double>(
      storage.face_flux_z_comoving.data() + record.flux_z_offset, binding.geometry.faceCountZ());
  view.patch_id = record.patch_id;
  view.geometry_generation = record.geometry_generation;
  view.refinement_level = record.level;
  return view;
}

MhdFaceFluxMutableView makeMutableMagneticFaceView(core::MagneticFaceSoa& storage,
                                                   const MhdAdapterPatchBinding& binding,
                                                   std::size_t patch_index) {
  requireBindingMatchesRecord(storage, binding, patch_index);
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
                         std::span<const std::uint32_t> dense_row_by_patch_cell,
                         const MhdFaceFluxView& faces,
                         double scale_factor,
                         MhdGasStateSoa& gas) {
  requireRowMap(state, geometry, dense_row_by_patch_cell);
  if (!faces.valid() || !(scale_factor > 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "invalid face view or scale factor on load"));
  }
  const std::size_t cell_count = geometry.cellCount();
  gas.resize(cell_count);
  const double volume = geometry.cellVolume();
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    const std::uint32_t row = dense_row_by_patch_cell[cell];
    const double rho = state.gas_cells.density_code[row];
    const double vx = state.gas_cells.velocity_x_peculiar[row];
    const double vy = state.gas_cells.velocity_y_peculiar[row];
    const double vz = state.gas_cells.velocity_z_peculiar[row];
    const double internal = state.gas_cells.internal_energy_code[row];
    const std::array<double, 3> b = cellFieldFromView(geometry, faces, cell);
    MhdConservedState conserved{};
    conserved.rho_comoving = rho;
    conserved.momentum_density_x_comoving = rho * vx;
    conserved.momentum_density_y_comoving = rho * vy;
    conserved.momentum_density_z_comoving = rho * vz;
    conserved.total_energy_density_comoving = rho * internal +
        0.5 * rho * (vx * vx + vy * vy + vz * vz) +
        0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor;
    conserved.metal_mass_density_comoving = state.gas_cells.metal_mass_code[row] / volume;
    if (!(rho > 0.0) || !(internal > 0.0) || !std::isfinite(conserved.total_energy_density_comoving)) {
      throw MhdError(MhdErrorCode::kNonPositiveState,
                     adapterContext("mhd_adapter", "canonical gas row is inadmissible on load"));
    }
    gas.storeGas(cell, conserved);
  }
}

void storeMhdGasToState(core::SimulationState& state,
                        const MhdPatchGeometry& geometry,
                        std::span<const std::uint32_t> dense_row_by_patch_cell,
                        const MhdGasStateSoa& gas,
                        const MhdFaceFluxView& faces,
                        double scale_factor,
                        double gamma_adiabatic) {
  requireRowMap(state, geometry, dense_row_by_patch_cell);
  if (gas.size() != geometry.cellCount() || !faces.valid() || !(scale_factor > 0.0) ||
      !(gamma_adiabatic > 1.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "invalid store inputs"));
  }
  const std::size_t cell_count = geometry.cellCount();
  // Validate every cell first so no row is written from an inadmissible state.
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    const MhdConservedState conserved = gas.loadGas(cell);
    const std::array<double, 3> b = cellFieldFromView(geometry, faces, cell);
    const double rho = conserved.rho_comoving;
    const double momentum2 =
        conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
        conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
        conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving;
    const double thermal = conserved.total_energy_density_comoving -
        (rho > 0.0 ? 0.5 * momentum2 / rho : 0.0) -
        0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor;
    if (!(rho > 0.0) || !(thermal > 0.0) || !std::isfinite(thermal)) {
      MhdFailureContext failure =
          adapterContext("mhd_adapter", "accepted MHD state has non-positive rho or thermal "
                                        "energy; no floor is applied");
      failure.scale_factor = scale_factor;
      throw MhdError(MhdErrorCode::kNonPositiveState, std::move(failure));
    }
  }
  const double volume = geometry.cellVolume();
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    const std::uint32_t row = dense_row_by_patch_cell[cell];
    const MhdConservedState conserved = gas.loadGas(cell);
    const std::array<double, 3> b = cellFieldFromView(geometry, faces, cell);
    const double rho = conserved.rho_comoving;
    const double vx = conserved.momentum_density_x_comoving / rho;
    const double vy = conserved.momentum_density_y_comoving / rho;
    const double vz = conserved.momentum_density_z_comoving / rho;
    const double thermal = conserved.total_energy_density_comoving -
        0.5 * rho * (vx * vx + vy * vy + vz * vz) -
        0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor;
    state.gas_cells.density_code[row] = rho;
    state.gas_cells.velocity_x_peculiar[row] = vx;
    state.gas_cells.velocity_y_peculiar[row] = vy;
    state.gas_cells.velocity_z_peculiar[row] = vz;
    state.gas_cells.internal_energy_code[row] = thermal / rho;
    state.gas_cells.pressure_code[row] = (gamma_adiabatic - 1.0) * thermal;
    state.gas_cells.metal_mass_code[row] = conserved.metal_mass_density_comoving * volume;
    state.cells.mass_code[row] = rho * volume;
  }
}

double patchMagneticEnergy(const core::MagneticFaceSoa& storage, std::size_t patch_index,
                           double scale_factor) {
  if (patch_index >= storage.patches.size() || !(scale_factor > 0.0)) {
    throw MhdError(MhdErrorCode::kInvalidField,
                   adapterContext("mhd_adapter", "invalid magnetic energy request"));
  }
  const core::MagneticFacePatchRecord& record = storage.patches[patch_index];
  MhdAdapterPatchBinding binding{};
  binding.patch_id = record.patch_id;
  binding.geometry_generation = record.geometry_generation;
  binding.level = record.level;
  binding.geometry = magneticPatchGeometry(record);
  const MhdFaceFluxView view = makeMagneticFaceView(storage, binding, patch_index);
  const MhdPatchGeometry& geometry = binding.geometry;
  double energy = 0.0;
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const std::array<double, 3> b = cellFieldFromView(geometry, view, cell);
    energy += 0.5 * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / scale_factor *
        geometry.cellVolume();
  }
  return energy;
}

std::unique_ptr<MhdCoefficientProvider> makeConfiguredCoefficientProvider(
    const core::MhdConfig& config) {
  std::unique_ptr<MhdCoefficientProvider> provider;
  switch (config.nonideal.coefficient_provider) {
    case core::MhdCoefficientProviderConfigMode::kConstant:
      if (!config.nonideal.coefficient_table_path.empty()) {
        throw MhdError(MhdErrorCode::kBadCoefficient,
                       adapterContext("mhd_coefficient_provider",
                                      "constant provider refuses a coefficient table path"));
      }
      provider = std::make_unique<ConstantMhdCoefficientProvider>(
          config.nonideal.eta_ohmic_phys_cgs, config.nonideal.eta_hall_phys_cgs,
          config.nonideal.eta_ambipolar_phys_cgs);
      break;
    case core::MhdCoefficientProviderConfigMode::kTable: {
      if (config.nonideal.eta_ohmic_phys_cgs != 0.0 || config.nonideal.eta_hall_phys_cgs != 0.0 ||
          config.nonideal.eta_ambipolar_phys_cgs != 0.0) {
        throw MhdError(MhdErrorCode::kBadCoefficient,
                       adapterContext("mhd_coefficient_provider",
                                      "table provider refuses constant eta_*_phys_cgs keys "
                                      "(mixed coefficient authority)"));
      }
      TableMhdCoefficientProvider::Options options{};
      options.composition_id = config.nonideal.composition_id;
      options.grain_model_id = config.nonideal.grain_model_id;
      options.expected_sha256 = config.nonideal.expected_table_sha256;
      provider = std::make_unique<TableMhdCoefficientProvider>(
          std::filesystem::path(config.nonideal.coefficient_table_path), std::move(options));
      break;
    }
    case core::MhdCoefficientProviderConfigMode::kChemistry:
      provider = std::make_unique<ChemistryMhdCoefficientProvider>();
      break;  // validateAdmission below fails closed (MHD_UNSUPPORTED_PHYSICS)
    case core::MhdCoefficientProviderConfigMode::kNone:
      throw MhdError(MhdErrorCode::kBadCoefficient,
                     adapterContext("mhd_coefficient_provider",
                                    "non-ideal transport requires a real coefficient provider"));
  }
  provider->validateAdmission();
  return provider;
}

}  // namespace cosmosim::mhd
