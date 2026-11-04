#include "cosmosim/workflows/time_coordinator.hpp"
#include "cosmosim/core/stage_transaction.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <cstring>
#include <limits>
#include <new>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/workflows/source_runtime.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/runtime_console_reporter.hpp"
#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/star_formation_geometry.hpp"
#include "cosmosim/workflows/migration_balance_runtime.hpp"
#include "cosmosim/workflows/output_restart_runtime.hpp"

namespace cosmosim::workflows {
namespace {

// Rollback-store accounting constants belong at namespace scope: C++20 local
// classes cannot declare static data members.
constexpr std::uint64_t k_multirate_scheduler_row_bytes =
    sizeof(std::uint8_t) + sizeof(std::uint64_t) + sizeof(std::uint8_t) + sizeof(std::uint8_t);
constexpr std::uint64_t k_multirate_rollback_vector_count =
    7U + TransactionalForceEndpoint::k_lane_count + 8U;
constexpr std::uint64_t k_multirate_rollback_allocation_overhead_bytes = 64U;

constexpr RuntimeEpochField k_particle_stage_epochs =
    RuntimeEpochField::kParticleIndex |
    RuntimeEpochField::kSchedulerTick |
    RuntimeEpochField::kStepIndex;
constexpr RuntimeEpochField k_state_stage_epochs =
    RuntimeEpochField::kParticleIndex |
    RuntimeEpochField::kCellIndex |
    RuntimeEpochField::kGasIdentity |
    RuntimeEpochField::kSchedulerTick |
    RuntimeEpochField::kStepIndex;

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}

// Rank-local rung histogram "n0,n1,...,nmax" of the committed assignment.
[[nodiscard]] std::string formatBinOccupancy(
    const core::HierarchicalTimeBinScheduler& scheduler, std::span<const std::uint32_t> rows) {
  std::vector<std::uint64_t> counts(static_cast<std::size_t>(scheduler.maxBin()) + 1U, 0U);
  for (const auto row : rows) {
    const auto bin = static_cast<std::size_t>(scheduler.binIndex(row));
    if (bin < counts.size()) ++counts[bin];
  }
  std::string text;
  for (std::size_t bin = 0; bin < counts.size(); ++bin) {
    if (bin != 0U) text += ',';
    text += std::to_string(counts[bin]);
  }
  return text;
}

void ensureSchedulerCoversState(
    std::size_t required_size,
    core::HierarchicalTimeBinScheduler& scheduler,
    std::string_view scheduler_name) {
  if (required_size > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error(std::string(scheduler_name) + " element count exceeds uint32 range");
  }
  const std::uint32_t required_count = static_cast<std::uint32_t>(required_size);
  if (required_count <= scheduler.elementCount()) {
    return;
  }
  const std::uint8_t new_element_bin = scheduler.maxBin() > 0 ? 1U : 0U;
  const std::uint64_t bin_period = 1ULL << new_element_bin;
  if (scheduler.currentTick() > std::numeric_limits<std::uint64_t>::max() - bin_period) {
    throw std::overflow_error(std::string(scheduler_name) + " next activation tick overflows uint64");
  }
  const std::uint64_t first_activation_tick =
      ((scheduler.currentTick() / bin_period) + 1ULL) * bin_period;
  scheduler.appendElements(
      required_count - scheduler.elementCount(), new_element_bin, first_activation_tick);
}

void initializeSchedulerBins(
    const core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler) {
  if (state.particles.size() > static_cast<std::size_t>(
          std::numeric_limits<std::uint32_t>::max()) ||
      state.cells.size() > static_cast<std::size_t>(
          std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error(
        "reference workflow scheduler element count exceeds uint32 range");
  }
  const std::uint32_t particle_count =
      static_cast<std::uint32_t>(state.particles.size());
  const std::uint32_t cell_count =
      static_cast<std::uint32_t>(state.cells.size());
  const std::uint8_t particle_default_bin =
      0U;
  const std::uint8_t gas_default_bin =
      gas_cell_scheduler.maxBin() > 0 ? 1U : 0U;
  particle_scheduler.reset(particle_count, particle_default_bin, 0U);
  gas_cell_scheduler.reset(cell_count, gas_default_bin, 0U);
  if (cell_count > 0U) {
    state.requireGasCellIdentityMapCoversDenseRows(
        "initialize independent gas-cell scheduler");
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      gas_cell_scheduler.setElementBin(
          cell_index, 0U, gas_cell_scheduler.currentTick());
    }
  }
}

void ensureSchedulersCoverState(
    const core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler) {
  ensureSchedulerCoversState(state.particles.size(), particle_scheduler, "particle scheduler");
  ensureSchedulerCoversState(state.cells.size(), gas_cell_scheduler, "gas-cell scheduler");
}

void syncTimeBinsFromSchedulers(
    const core::HierarchicalTimeBinScheduler& particle_scheduler,
    const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    core::SimulationState& state) {
  core::syncTimeBinMirrorsFromScheduler(
      particle_scheduler, state, core::TimeBinMirrorDomain::kParticles);
  core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, state);
}

[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
}

[[nodiscard]] gravity::TreeSofteningSpeciesPolicy speciesSofteningByTag(
    const core::SimulationConfig& config) {
  gravity::TreeSofteningSpeciesPolicy values{};
  values.enabled = config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
  values.epsilon_comoving_by_species.fill(config.numerics.gravity_softening_kpc_comoving * 1.0e-3);
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] =
      (config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_dark_matter_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] =
      (config.numerics.gravity_softening_gas_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_gas_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)] =
      (config.numerics.gravity_softening_star_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_star_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)] =
      (config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_black_hole_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)] =
      (config.numerics.gravity_softening_tracer_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_tracer_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)];
  return values;
}

struct AdaptiveTimeStepCriteriaStorage {
  std::vector<std::uint32_t> gas_particle_index_by_cell;
  std::vector<std::uint64_t> gas_cell_id_by_cell;
  std::vector<std::uint64_t> patch_id_by_cell;
  std::vector<std::uint32_t> patch_row_by_cell;
  std::vector<double> cell_width_x_code;
  std::vector<double> cell_width_y_code;
  std::vector<double> cell_width_z_code;
  std::vector<double> velocity_divergence_code;
  std::vector<double> metal_diffusion_dt_code;
  std::vector<double> hydro_cfl_level_relaxation;
};

struct LocalGasCellCflMetadata {
  std::vector<double> cell_width_x_code;
  std::vector<double> cell_width_y_code;
  std::vector<double> cell_width_z_code;
  std::vector<std::uint64_t> patch_id_by_cell;
  std::vector<std::uint32_t> patch_row_by_cell;
};

[[nodiscard]] internal::CartesianGasCellRowLayout requireCartesianGasCellRowLayout(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::string_view caller) {
  internal::CartesianGasCellLayoutBuildResult result =
      internal::buildCartesianGasCellRowLayout(state, config);
  if (!result.ok()) {
    throw std::runtime_error(
        std::string(caller) + ": fixed Cartesian hydro geometry rejected: " + result.diagnostic);
  }
  return std::move(result.layout);
}

[[nodiscard]] LocalGasCellCflMetadata buildLocalGasCellCflMetadata(
    const core::SimulationState& state,
    const core::SimulationConfig& config) {
  LocalGasCellCflMetadata metadata;
  const std::size_t cell_count = state.cells.size();
  metadata.cell_width_x_code.resize(cell_count);
  metadata.cell_width_y_code.resize(cell_count);
  metadata.cell_width_z_code.resize(cell_count);
  metadata.patch_id_by_cell.assign(cell_count, 0U);
  metadata.patch_row_by_cell.assign(cell_count, std::numeric_limits<std::uint32_t>::max());
  if (cell_count == 0U) {
    return metadata;
  }
  state.requireGasCellIdentityMapCoversDenseRows("hydro CFL metadata construction");

  if (amr::hasProductionAmrHydroCoverage(state)) {
    const std::vector<amr::PatchDescriptor> descriptors = amr::buildProductionAmrPatchDescriptors(state);
    for (const amr::PatchDescriptor& descriptor : descriptors) {
      const amr::AmrHydroPatchGeometry patch_geometry = amr::buildAmrHydroPatchGeometry(state, descriptor);
      for (const amr::AmrHydroCellDescriptor& cell : patch_geometry.real_cells) {
        if (cell.local_cell_row >= cell_count || cell.patch_local_cell == hydro::k_invalid_cell_index) {
          throw std::runtime_error("AMR CFL metadata received an invalid physical cell descriptor");
        }
        const auto* identity = state.gas_cell_identity.findByLocalRow(cell.local_cell_row);
        if (identity == nullptr || identity->gas_cell_id != cell.gas_cell_id ||
            identity->owning_patch_id != descriptor.patch_id) {
          throw std::runtime_error("AMR CFL metadata rejected stale gas-cell identity or patch ownership");
        }
        const std::uint32_t row = cell.local_cell_row;
        metadata.cell_width_x_code[row] = patch_geometry.geometry.cell_width_x_comoving;
        metadata.cell_width_y_code[row] = patch_geometry.geometry.cell_width_y_comoving;
        metadata.cell_width_z_code[row] = patch_geometry.geometry.cell_width_z_comoving;
        metadata.patch_id_by_cell[row] = descriptor.patch_id;
        metadata.patch_row_by_cell[row] = static_cast<std::uint32_t>(cell.patch_local_cell);
      }
    }
    for (std::uint32_t row = 0; row < cell_count; ++row) {
      if (!std::isfinite(metadata.cell_width_x_code[row]) || metadata.cell_width_x_code[row] <= 0.0 ||
          !std::isfinite(metadata.cell_width_y_code[row]) || metadata.cell_width_y_code[row] <= 0.0 ||
          !std::isfinite(metadata.cell_width_z_code[row]) || metadata.cell_width_z_code[row] <= 0.0 ||
          metadata.patch_id_by_cell[row] == 0U ||
          metadata.patch_row_by_cell[row] == std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("AMR CFL metadata did not cover every authoritative gas-cell row");
      }
    }
    return metadata;
  }

  const internal::CartesianGasCellRowLayout layout = requireCartesianGasCellRowLayout(
      state, config, "hydro CFL metadata construction");
  for (std::uint32_t row = 0; row < cell_count; ++row) {
    metadata.cell_width_x_code[row] = layout.spec.cell_width_x_comoving;
    metadata.cell_width_y_code[row] = layout.spec.cell_width_y_comoving;
    metadata.cell_width_z_code[row] = layout.spec.cell_width_z_comoving;
    const auto* identity = state.gas_cell_identity.findByLocalRow(row);
    if (identity == nullptr || identity->gas_cell_id == 0U) {
      throw std::runtime_error("fixed-patch CFL metadata rejected incomplete gas-cell identity coverage");
    }
    if (state.patches.size() != 0U) {
      if (identity->owning_patch_id == 0U) {
        throw std::runtime_error("fixed-patch CFL metadata rejected a gas cell without explicit patch ownership");
      }
      bool patch_found = false;
      for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
        if (state.patches.patch_id[patch] == identity->owning_patch_id) {
          if (state.cells.patch_index[row] != patch) {
            throw std::runtime_error("fixed-patch CFL metadata rejected stale dense patch-index mirror");
          }
          patch_found = true;
          break;
        }
      }
      if (!patch_found) {
        throw std::runtime_error("fixed-patch CFL metadata rejected an identity patch absent from PatchSoa");
      }
      metadata.patch_id_by_cell[row] = identity->owning_patch_id;
      metadata.patch_row_by_cell[row] = layout.geometry_row_by_dense_row.at(row);
    }
  }
  return metadata;
}

[[nodiscard]] core::AdaptiveTimeStepCriteriaView buildAdaptiveTimeStepCriteriaView(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z,
    double scale_factor,
    AdaptiveTimeStepCriteriaStorage& storage,
    int recursive_amr_global_min_level = -1,
    std::uint32_t amr_refinement_ratio = 2U) {
  state.requireGasCellIdentityMapCoversDenseRows("adaptive time-bin view construction");
  // G1.2: level relaxation of the hydro CFL candidate under recursive
  // subcycling (global minimum level agreed by the caller on every rank).
  storage.hydro_cfl_level_relaxation.clear();
  if (recursive_amr_global_min_level >= 0 && state.cells.size() != 0U) {
    storage.hydro_cfl_level_relaxation.assign(state.cells.size(), 1.0);
    for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
      const std::uint32_t patch = row < state.cells.patch_index.size() ? state.cells.patch_index[row]
                                                                       : std::numeric_limits<std::uint32_t>::max();
      if (patch >= state.patches.size()) {
        throw std::runtime_error("recursive AMR CFL relaxation: gas cell without patch ownership");
      }
      const int level = state.patches.level[patch];
      if (level < recursive_amr_global_min_level) {
        throw std::runtime_error("recursive AMR CFL relaxation: cell level below the global minimum level");
      }
      storage.hydro_cfl_level_relaxation[row] =
          std::pow(static_cast<double>(amr_refinement_ratio), static_cast<double>(level - recursive_amr_global_min_level));
    }
  }
  // Audit F4: the gas CFL candidate reads GasCellSidecar::sound_speed_code.
  // A source stage that changed density/internal energy must have run the
  // thermodynamic refresh contract before the next timestep is selected.
  if (!state.gasThermodynamicCachesFresh()) {
    throw std::logic_error(
        "timestep selection refused stale gas thermodynamic caches: a source stage mutated gas "
        "thermodynamics without hydro::refreshGasThermodynamics() (mutation epoch " +
        std::to_string(state.gasThermodynamicMutationEpoch()) + ", refresh epoch " +
        std::to_string(state.gasThermodynamicRefreshEpoch()) + ")");
  }
  const auto particle_row_by_id = internal::buildParticleRowById(state);
  constexpr std::uint32_t k_no_parent_particle = std::numeric_limits<std::uint32_t>::max();
  storage.gas_particle_index_by_cell.clear();
  storage.gas_particle_index_by_cell.reserve(state.cells.size());
  storage.gas_cell_id_by_cell.assign(state.gas_cells.gas_cell_id.begin(), state.gas_cells.gas_cell_id.end());
  LocalGasCellCflMetadata cfl_metadata = buildLocalGasCellCflMetadata(state, config);
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const auto parent_row = internal::parentParticleRowForGasCellRow(
        state, cell_index, particle_row_by_id, "adaptive time-bin view construction");
    storage.gas_particle_index_by_cell.push_back(parent_row.value_or(k_no_parent_particle));
  }
  storage.patch_id_by_cell = std::move(cfl_metadata.patch_id_by_cell);
  storage.patch_row_by_cell = std::move(cfl_metadata.patch_row_by_cell);
  storage.cell_width_x_code = std::move(cfl_metadata.cell_width_x_code);
  storage.cell_width_y_code = std::move(cfl_metadata.cell_width_y_code);
  storage.cell_width_z_code = std::move(cfl_metadata.cell_width_z_code);
  storage.velocity_divergence_code.assign(
      state.cells.size(), std::numeric_limits<double>::quiet_NaN());
  storage.metal_diffusion_dt_code.assign(
      state.cells.size(), std::numeric_limits<double>::infinity());
  const double length_to_physical = config.units.coordinate_frame ==
          core::CoordinateFrame::kComoving
      ? std::max(scale_factor, 1.0e-12)
      : 1.0;
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const internal::StarFormationPatchCellGeometry geometry =
        internal::starFormationPatchCellGeometry(state, cell_index);
    if (!geometry.valid) {
      continue;
    }
    const std::array<std::span<const double>, 3> velocity_fields{
        state.gas_cells.velocity_x_peculiar,
        state.gas_cells.velocity_y_peculiar,
        state.gas_cells.velocity_z_peculiar};
    const std::array<double, 3> spacing_phys_code{
        geometry.dx_stored * length_to_physical,
        geometry.dy_stored * length_to_physical,
        geometry.dz_stored * length_to_physical};
    physics::MetalDiffusionVelocityGradient velocity_gradient;
    bool gradient_valid = true;
    for (std::size_t component = 0; component < velocity_fields.size(); ++component) {
      for (std::size_t axis = 0; axis < spacing_phys_code.size(); ++axis) {
        velocity_gradient.grad[component][axis] =
            internal::starFormationDerivativeAtCell(
                velocity_fields[component], geometry, static_cast<int>(axis),
                spacing_phys_code[axis], cell_index);
        gradient_valid = gradient_valid &&
            std::isfinite(velocity_gradient.grad[component][axis]);
      }
    }
    if (gradient_valid) {
      storage.velocity_divergence_code[cell_index] =
          velocity_gradient.grad[0][0] + velocity_gradient.grad[1][1] +
          velocity_gradient.grad[2][2];
      if (config.physics.enable_metal_diffusion &&
          config.physics.metal_diffusion_model ==
              core::MetalDiffusionModel::kSmagorinsky) {
        const double filter_length = std::cbrt(
            spacing_phys_code[0] * spacing_phys_code[1] * spacing_phys_code[2]);
        const double strain = physics::traceFreeStrainMagnitude(velocity_gradient);
        const double kappa = std::clamp(
            config.physics.metal_diffusion_coefficient * filter_length *
                filter_length * strain,
            config.physics.metal_diffusion_coefficient_floor_code,
            config.physics.metal_diffusion_coefficient_ceiling_code);
        if (kappa > 0.0) {
          const double minimum_spacing_squared = std::min({
              spacing_phys_code[0] * spacing_phys_code[0],
              spacing_phys_code[1] * spacing_phys_code[1],
              spacing_phys_code[2] * spacing_phys_code[2]});
          storage.metal_diffusion_dt_code[cell_index] =
              config.physics.metal_diffusion_cfl * minimum_spacing_squared /
              (6.0 * kappa);
        }
      }
    }
  }
  return core::AdaptiveTimeStepCriteriaView{
      .particles = core::TimeStepParticleCriteriaView{
          .velocity_x_peculiar = state.particles.velocity_x_peculiar,
          .velocity_y_peculiar = state.particles.velocity_y_peculiar,
          .velocity_z_peculiar = state.particles.velocity_z_peculiar,
          .species_tag = state.particle_sidecar.species_tag,
          .homogeneous_dmo_species = state.hasHomogeneousDmoMetadata(),
          .gravity_softening_comoving = state.particle_sidecar.gravity_softening_comoving,
          .accel_x_comoving = particle_accel_x,
          .accel_y_comoving = particle_accel_y,
          .accel_z_comoving = particle_accel_z,
          .black_hole_particle_index = state.black_holes.particle_index,
          .black_hole_subgrid_mass_code = state.black_holes.subgrid_mass_code,
          .black_hole_accretion_rate_code = state.black_holes.accretion_rate_code,
      },
      .gas_cells = core::TimeStepGasCellCriteriaView{
          .gas_particle_index_by_cell = storage.gas_particle_index_by_cell,
          .gas_cell_id_by_cell = storage.gas_cell_id_by_cell,
          .patch_id_by_cell = storage.patch_id_by_cell,
          .patch_row_by_cell = storage.patch_row_by_cell,
          .cell_width_x_code = storage.cell_width_x_code,
          .cell_width_y_code = storage.cell_width_y_code,
          .cell_width_z_code = storage.cell_width_z_code,
          .cell_mass_code = state.cells.mass_code,
          .velocity_x_peculiar = state.gas_cells.velocity_x_peculiar,
          .velocity_y_peculiar = state.gas_cells.velocity_y_peculiar,
          .velocity_z_peculiar = state.gas_cells.velocity_z_peculiar,
          .density_code = state.gas_cells.density_code,
          .temperature_code = state.gas_cells.temperature_code,
          .sound_speed_code = state.gas_cells.sound_speed_code,
          .velocity_divergence_code = storage.velocity_divergence_code,
          .metal_diffusion_dt_code = storage.metal_diffusion_dt_code,
          .accel_x_comoving = cell_accel_x,
          .accel_y_comoving = cell_accel_y,
          .accel_z_comoving = cell_accel_z,
          .hydro_cfl_level_relaxation = storage.hydro_cfl_level_relaxation,
      },
  };
}

enum class SchedulerElementFamily {
  kParticles,
  kGasCells,
};

[[nodiscard]] double updateAdaptiveTimeBinsFromView(
    const core::AdaptiveTimeStepCriteriaView& view,
    core::HierarchicalTimeBinScheduler& scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    const physics::EffectiveMultiphaseEosTable* effective_eos_table,
    SchedulerElementFamily element_family,
    std::span<const std::uint32_t> requested_elements,
    bool update_all_elements,
    TimestepLimiterRecord* limiter = nullptr) {
  if (integrator_state.dt_time_code <= 0.0) {
    throw std::invalid_argument("adaptive time-bin update requires dt_time_code > 0");
  }
  const core::TimeStepLimits limits{
      .min_dt_time_code = integrator_state.dt_time_code,
      .max_dt_time_code = integrator_state.dt_time_code * static_cast<double>(1ULL << scheduler.maxBin()),
      .max_bin = scheduler.maxBin(),
  };
  const std::size_t particle_count = view.particles.velocity_x_peculiar.size();
  const std::size_t cell_count = view.gas_cells.cell_mass_code.size();
  if (view.particles.velocity_y_peculiar.size() != particle_count ||
      view.particles.velocity_z_peculiar.size() != particle_count ||
      (!view.particles.homogeneous_dmo_species && view.particles.species_tag.size() != particle_count) ||
      (view.particles.homogeneous_dmo_species && !view.particles.species_tag.empty())) {
    throw std::invalid_argument("particle timestep criteria view has mismatched extents");
  }
  if (view.gas_cells.density_code.size() != cell_count ||
      view.gas_cells.velocity_x_peculiar.size() != cell_count ||
      view.gas_cells.velocity_y_peculiar.size() != cell_count ||
      view.gas_cells.velocity_z_peculiar.size() != cell_count ||
      view.gas_cells.temperature_code.size() != cell_count ||
      view.gas_cells.sound_speed_code.size() != cell_count ||
      view.gas_cells.velocity_divergence_code.size() != cell_count ||
      view.gas_cells.metal_diffusion_dt_code.size() != cell_count ||
      view.gas_cells.gas_particle_index_by_cell.size() != cell_count ||
      view.gas_cells.gas_cell_id_by_cell.size() != cell_count ||
      view.gas_cells.patch_id_by_cell.size() != cell_count ||
      view.gas_cells.patch_row_by_cell.size() != cell_count ||
      view.gas_cells.cell_width_x_code.size() != cell_count ||
      view.gas_cells.cell_width_y_code.size() != cell_count ||
      view.gas_cells.cell_width_z_code.size() != cell_count) {
    throw std::invalid_argument("gas-cell timestep criteria view has mismatched extents");
  }
  const auto species_softening = speciesSofteningByTag(config);
  const double global_softening = config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
  const core::UnitSystem runtime_units = core::makeUnitSystem(
      config.units.length_unit,
      config.units.mass_unit,
      config.units.velocity_unit);
  const double newton_g_code = newtonGCodeFromUnits(runtime_units);
  const double gravity_scale_factor = std::max(integrator_state.current_scale_factor, 1.0e-12);
  std::optional<double> cosmology_dt;
  const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
  if (mode_policy.cosmological_comoving_frame &&
      config.units.coordinate_frame == core::CoordinateFrame::kComoving &&
      integrator_state.current_scale_factor > 0.0 &&
      config.cosmology.hubble_param > 0.0) {
    core::CosmologyBackgroundConfig background_config;
    background_config.hubble_param = config.cosmology.hubble_param;
    background_config.omega_matter = config.cosmology.omega_matter;
    background_config.omega_lambda = config.cosmology.omega_lambda;
    const core::LambdaCdmBackground background(background_config);
    cosmology_dt = core::computeCosmologyExpansionTimeStep(
        background,
        integrator_state.current_scale_factor,
        config.numerics.cosmology_max_delta_ln_a,
        config.numerics.cosmology_max_hubble_time_fraction,
        integrator_state.time_si_per_code);
  }
  const auto star_formation_dt_for_cell = [&](std::uint32_t cell_index) -> std::optional<double> {
    if (!config.physics.enable_star_formation || cell_index >= cell_count ||
        (config.physics.star_formation_model !=
             core::StarFormationModelKind::kEffectiveMultiphaseTngLike &&
         config.physics.sf_epsilon_ff <= 0.0)) {
      return std::nullopt;
    }
    const double gas_mass = view.gas_cells.cell_mass_code[cell_index];
    const double stored_density = view.gas_cells.density_code[cell_index];
    const double temperature = view.gas_cells.temperature_code[cell_index];
    if (!(gas_mass > 0.0) || !(stored_density > 0.0) ||
        !std::isfinite(gas_mass) || !std::isfinite(stored_density) || !std::isfinite(temperature)) {
      return std::nullopt;
    }
    if (config.physics.star_formation_model ==
        core::StarFormationModelKind::kLegacySchmidtThreshold) {
      const double scale_factor = std::max(integrator_state.current_scale_factor, 1.0e-12);
      const double physical_density = config.units.coordinate_frame == core::CoordinateFrame::kComoving
          ? stored_density / (scale_factor * scale_factor * scale_factor)
          : stored_density;
      if (physical_density < config.physics.sf_density_threshold_code ||
          temperature > config.physics.sf_temperature_threshold_k) {
        return std::nullopt;
      }
      const double t_ff_code = std::sqrt(3.0 * core::constants::k_pi /
          (32.0 * newton_g_code * std::max(physical_density, 1.0e-30)));
      return std::max(
          1.0e-12,
          config.numerics.source_max_fractional_change * t_ff_code /
              config.physics.sf_epsilon_ff);
    }

    const double scale_factor = std::max(integrator_state.current_scale_factor, 1.0e-12);
    const double physical_density = config.units.coordinate_frame == core::CoordinateFrame::kComoving
        ? stored_density / (scale_factor * scale_factor * scale_factor)
        : stored_density;
    if (config.physics.star_formation_model ==
        core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
      if (effective_eos_table == nullptr) return std::nullopt;
      const auto equilibrium = effective_eos_table->lookup(physical_density);
      if (!equilibrium.above_threshold || !equilibrium.valid ||
          equilibrium.entry.cold_mass_fraction <= 0.0) {
        return std::nullopt;
      }
      const double long_lived_factor =
          config.physics.sf_effective_birth_mass_convention ==
              core::EffectiveIsmBirthMassConvention::kLongLivedMass
          ? (1.0 - config.physics.sf_effective_massive_star_fraction)
          : 1.0;
      const double rate_per_mass = long_lived_factor *
          equilibrium.entry.cold_mass_fraction /
          std::max(equilibrium.entry.star_formation_timescale_code, 1.0e-30);
      const double maximum_fraction = std::clamp(
          std::min(config.physics.sf_max_fractional_mass_conversion,
                   config.numerics.source_max_fractional_change),
          1.0e-12, 1.0 - 1.0e-12);
      return std::max(1.0e-12, -std::log1p(-maximum_fraction) / rate_per_mass);
    }
    if (config.physics.sf_temperature_safety_ceiling_k > 0.0 &&
        temperature > config.physics.sf_temperature_safety_ceiling_k) {
      return std::nullopt;
    }
    const double t_ff_code = std::sqrt(3.0 * core::constants::k_pi /
        (32.0 * newton_g_code * std::max(physical_density, 1.0e-30)));
    double collapse_time_code = t_ff_code;
    if (config.physics.sf_collapse_timescale ==
        core::StarFormationCollapseTimescale::kMinimumFreeFallOrCompression) {
      const double divergence_code = view.gas_cells.velocity_divergence_code[cell_index];
      if (std::isfinite(divergence_code) && divergence_code < 0.0) {
        collapse_time_code = std::min(t_ff_code, -1.0 / divergence_code);
      }
    }
    const double maximum_fraction = std::clamp(
        std::min(
            config.physics.sf_max_fractional_mass_conversion,
            config.numerics.source_max_fractional_change),
        1.0e-12,
        1.0 - 1.0e-12);
    const double dt_limit = -collapse_time_code * std::log1p(-maximum_fraction) /
        config.physics.sf_epsilon_ff;
    return std::max(1.0e-12, dt_limit);
  };
  const std::size_t bh_count = view.particles.black_hole_particle_index.size();
  if (view.particles.black_hole_subgrid_mass_code.size() != bh_count ||
      view.particles.black_hole_accretion_rate_code.size() != bh_count) {
    throw std::invalid_argument("black-hole timestep criteria view has mismatched extents");
  }
  std::vector<std::size_t> bh_row_by_particle;
  if (config.physics.enable_black_hole_agn && bh_count != 0U) {
    bh_row_by_particle.assign(
        particle_count, std::numeric_limits<std::size_t>::max());
    for (std::size_t bh_index = 0; bh_index < bh_count; ++bh_index) {
      const std::uint32_t particle_index =
          view.particles.black_hole_particle_index[bh_index];
      if (particle_index >= particle_count) {
        throw std::invalid_argument(
            "black-hole sidecar references a particle row outside the particle extent");
      }
      if (bh_row_by_particle[particle_index] !=
          std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument(
            "black-hole sidecar contains duplicate particle-row ownership");
      }
      bh_row_by_particle[particle_index] = bh_index;
    }
  }
  const auto black_hole_dt_for_particle = [&](std::uint32_t particle_index) -> std::optional<double> {
    if (!config.physics.enable_black_hole_agn || view.particles.homogeneous_dmo_species ||
        particle_index >= view.particles.species_tag.size() ||
        view.particles.species_tag[particle_index] != static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
      return std::nullopt;
    }
    if (bh_row_by_particle.empty()) {
      return std::nullopt;
    }
    const std::size_t bh_index = bh_row_by_particle[particle_index];
    if (bh_index == std::numeric_limits<std::size_t>::max()) {
      return std::nullopt;
    }
    const double mass = std::max(
        view.particles.black_hole_subgrid_mass_code[bh_index], 1.0e-30);
    const double mdot = std::max(
        view.particles.black_hole_accretion_rate_code[bh_index], 0.0);
    if (mdot <= 0.0) {
      return std::nullopt;
    }
    return std::max(
        1.0e-12, config.numerics.source_max_fractional_change * mass / mdot);
  };
  double local_min_dt_time_code = std::numeric_limits<double>::infinity();
  TimestepLimiterRecord local_limiter{};
  local_limiter.dt_time_code = std::numeric_limits<double>::infinity();
  const auto consider_dt = [&](double value, const char* criterion = "unspecified", std::uint32_t element = 0U,
                               int amr_level = -1, double relaxation = 1.0) {
    if (std::isinf(value) && value > 0.0) {
      return;  // Explicit no-limit sentinel; another criterion or endpoint must bound the step.
    }
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::runtime_error(
          "adaptive timestep criterion produced an invalid value; expected finite positive or +infinity");
    }
    if (value < local_min_dt_time_code) {
      local_limiter = TimestepLimiterRecord{
          .dt_time_code = value, .criterion = criterion,
          .gas_cell = element_family == SchedulerElementFamily::kGasCells, .element = element,
          .amr_level = amr_level, .hydro_cfl_level_relaxation = relaxation};
    }
    local_min_dt_time_code = std::min(local_min_dt_time_code, value);
  };
  struct LimiterPublisher {
    TimestepLimiterRecord* out;
    const TimestepLimiterRecord& record;
    ~LimiterPublisher() {
      if (out != nullptr && record.dt_time_code < out->dt_time_code) *out = record;
    }
  } limiter_publisher{limiter, local_limiter};
  if (cosmology_dt.has_value()) consider_dt(*cosmology_dt, "cosmology_expansion");
  const std::size_t expected_element_count =
      element_family == SchedulerElementFamily::kParticles ? particle_count : cell_count;
  if (scheduler.elementCount() != expected_element_count) {
    throw std::invalid_argument(
        "adaptive time-bin scheduler extent does not match its declared element family");
  }
  const auto for_each_requested_element = [&](auto&& callback) {
    if (update_all_elements) {
      for (std::uint32_t element = 0; element < expected_element_count;
           ++element) {
        callback(element);
      }
      return;
    }
    for (const std::uint32_t element : requested_elements) {
      if (element >= expected_element_count) {
        throw std::out_of_range(
            "active timestep-criteria element is outside scheduler extent");
      }
      callback(element);
    }
  };
  if (element_family == SchedulerElementFamily::kGasCells) {
    for_each_requested_element([&](std::uint32_t cell_index) {
      const std::uint32_t gas_index = view.gas_cells.gas_particle_index_by_cell[cell_index];
      if (gas_index != std::numeric_limits<std::uint32_t>::max() && gas_index >= particle_count) {
        throw std::out_of_range("gas-particle index in timestep criteria view is out of range");
      }
      const double vx = view.gas_cells.velocity_x_peculiar[cell_index];
      const double vy = view.gas_cells.velocity_y_peculiar[cell_index];
      const double vz = view.gas_cells.velocity_z_peculiar[cell_index];
      const std::array<double, 3> cell_widths{
          view.gas_cells.cell_width_x_code[cell_index],
          view.gas_cells.cell_width_y_code[cell_index],
          view.gas_cells.cell_width_z_code[cell_index]};
      // Same rule as the accepted-step guard in HydroAmrRuntime
      // (core::hydroCflStabilityRuleFor / activeHydroAxes).
      const core::DirectionalCflTimeStepInput hydro_cfl_input{
          .cell_width_axis_code = cell_widths,
          .velocity_axis_code = {vx, vy, vz},
          .sound_speed_code = std::max(view.gas_cells.sound_speed_code[cell_index], 0.0),
          .coordinate_frame = config.units.coordinate_frame,
          .scale_factor = gravity_scale_factor,
          .stability_rule = core::hydroCflStabilityRuleFor(config),
          .active_axis = core::activeHydroAxes(cell_widths, core::hydroDomainExtentCode(config)),
      };
      // G1.2: under recursive subcycling the coarse interval is bounded by
      // r^(L - L_min) times the cell's own CFL step (it is subcycled).
      const double cfl_relaxation = view.gas_cells.hydro_cfl_level_relaxation.empty()
          ? 1.0 : view.gas_cells.hydro_cfl_level_relaxation[cell_index];
      if (!(cfl_relaxation >= 1.0) || !std::isfinite(cfl_relaxation)) {
        throw std::runtime_error("hydro CFL level relaxation must be finite and >= 1");
      }
      const double cfl_dt =
          core::computeDirectionalCflTimeStep(hydro_cfl_input, 0.4) * cfl_relaxation;
      const int cfl_level = view.gas_cells.hydro_cfl_level_relaxation.empty()
          ? -1 : static_cast<int>(std::lround(std::log2(cfl_relaxation)));
      consider_dt(cfl_dt, "gas_cell_hydro_cfl", cell_index, cfl_level, cfl_relaxation);
      const double ax = (cell_index < view.gas_cells.accel_x_comoving.size()) ? view.gas_cells.accel_x_comoving[cell_index] : 0.0;
      const double ay = (cell_index < view.gas_cells.accel_y_comoving.size()) ? view.gas_cells.accel_y_comoving[cell_index] : 0.0;
      const double az = (cell_index < view.gas_cells.accel_z_comoving.size()) ? view.gas_cells.accel_z_comoving[cell_index] : 0.0;
      const double amag = std::sqrt(ax * ax + ay * ay + az * az);
      const double eps = gas_index != std::numeric_limits<std::uint32_t>::max() &&
              !view.particles.gravity_softening_comoving.empty()
          ? view.particles.gravity_softening_comoving[gas_index]
          : species_softening.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)];
      const double gravity_dt = core::computeComovingGravityTimeStep(
          {.softening_length_comoving_code = std::max(eps, 1.0e-12),
           .scale_free_acceleration_magnitude_code = amag,
           .scale_factor = gravity_scale_factor},
          config.numerics.gravity_timestep_eta);
      consider_dt(gravity_dt, "gas_cell_gravity_acceleration", cell_index);
      scheduler.submitCandidateTimeStep(
          cell_index, cfl_dt, limits, core::TimeStepCandidateSource::kHydroCfl, "gas_cell_hydro_cfl");
      scheduler.submitCandidateTimeStep(
          cell_index, gravity_dt, limits, core::TimeStepCandidateSource::kGravityAcceleration, "gas_cell_gravity_acceleration");
      if (cosmology_dt.has_value()) {
        scheduler.submitCandidateTimeStep(
            cell_index, *cosmology_dt, limits, core::TimeStepCandidateSource::kCosmologyExpansion, "gas_cell_cosmology_expansion");
      }
      if (const auto source_dt = star_formation_dt_for_cell(cell_index); source_dt.has_value()) {
        consider_dt(*source_dt, "gas_cell_star_formation_source", cell_index);
        scheduler.submitCandidateTimeStep(
            cell_index, *source_dt, limits, core::TimeStepCandidateSource::kSourceTerm, "gas_cell_star_formation_source");
      }
      const double diffusion_dt = view.gas_cells.metal_diffusion_dt_code[cell_index];
      if (config.physics.enable_metal_diffusion && std::isfinite(diffusion_dt) &&
          diffusion_dt > 0.0) {
        consider_dt(diffusion_dt, "gas_cell_metal_diffusion_parabolic", cell_index);
        scheduler.submitCandidateTimeStep(
            cell_index, diffusion_dt, limits, core::TimeStepCandidateSource::kSourceTerm,
            "gas_cell_metal_diffusion_parabolic");
      }
    });
    return local_min_dt_time_code;
  }

  for_each_requested_element([&](std::uint32_t particle_index) {
    const double ax = (particle_index < view.particles.accel_x_comoving.size()) ? view.particles.accel_x_comoving[particle_index] : 0.0;
    const double ay = (particle_index < view.particles.accel_y_comoving.size()) ? view.particles.accel_y_comoving[particle_index] : 0.0;
    const double az = (particle_index < view.particles.accel_z_comoving.size()) ? view.particles.accel_z_comoving[particle_index] : 0.0;
    const double amag = std::sqrt(ax * ax + ay * ay + az * az);
    const std::size_t particle_species = view.particles.homogeneous_dmo_species
        ? static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)
        : static_cast<std::size_t>(view.particles.species_tag[particle_index]);
    const double eps = !view.particles.gravity_softening_comoving.empty()
        ? view.particles.gravity_softening_comoving[particle_index]
        : ((particle_species < species_softening.epsilon_comoving_by_species.size())
              ? species_softening.epsilon_comoving_by_species[particle_species]
              : global_softening);
    const double gravity_dt = core::computeComovingGravityTimeStep(
        {.softening_length_comoving_code = std::max(eps, 1.0e-12),
         .scale_free_acceleration_magnitude_code = amag,
         .scale_factor = gravity_scale_factor},
        config.numerics.gravity_timestep_eta);
    consider_dt(gravity_dt, "particle_gravity_acceleration", particle_index);
    scheduler.submitCandidateTimeStep(
        particle_index,
        gravity_dt,
        limits,
        core::TimeStepCandidateSource::kGravityAcceleration,
        "particle_gravity_acceleration");
    if (cosmology_dt.has_value()) {
      scheduler.submitCandidateTimeStep(
          particle_index,
          *cosmology_dt,
          limits,
          core::TimeStepCandidateSource::kCosmologyExpansion,
          "particle_cosmology_expansion");
    }
    if (const auto source_dt = black_hole_dt_for_particle(particle_index); source_dt.has_value()) {
      consider_dt(*source_dt, "particle_black_hole_source", particle_index);
      scheduler.submitCandidateTimeStep(
          particle_index,
          *source_dt,
          limits,
          core::TimeStepCandidateSource::kSourceTerm,
          "particle_black_hole_source");
    }
  });
  return local_min_dt_time_code;
}

[[nodiscard]] double updateAdaptiveTimeBinFamilies(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    const physics::EffectiveMultiphaseEosTable* effective_eos_table,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices,
    bool update_all_elements,
    int recursive_amr_global_min_level = -1,
    TimestepLimiterRecord* limiter = nullptr) {
  AdaptiveTimeStepCriteriaStorage storage;
  const core::AdaptiveTimeStepCriteriaView view = buildAdaptiveTimeStepCriteriaView(
      state,
      config,
      particle_accel_x,
      particle_accel_y,
      particle_accel_z,
      cell_accel_x,
      cell_accel_y,
      cell_accel_z,
      integrator_state.current_scale_factor,
      storage,
      recursive_amr_global_min_level,
      amr::ProductionAmrHydroOptions{}.refinement_ratio);
  const double particle_min_dt = updateAdaptiveTimeBinsFromView(
      view,
      particle_scheduler,
      integrator_state,
      config,
      effective_eos_table,
      SchedulerElementFamily::kParticles,
      active_particle_indices,
      update_all_elements,
      limiter);
  const double gas_min_dt = updateAdaptiveTimeBinsFromView(
      view,
      gas_cell_scheduler,
      integrator_state,
      config,
      effective_eos_table,
      SchedulerElementFamily::kGasCells,
      active_cell_indices,
      update_all_elements,
      limiter);
  return std::min(particle_min_dt, gas_min_dt);
}



}  // namespace

RungZeroTimeState::RungZeroTimeState(std::uint8_t max_bin)
    : m_particle_scheduler(max_bin),
      m_gas_cell_scheduler(max_bin) {}

core::HierarchicalTimeBinScheduler&
RungZeroTimeState::particleScheduler() noexcept {
  return m_particle_scheduler;
}

const core::HierarchicalTimeBinScheduler&
RungZeroTimeState::particleScheduler() const noexcept {
  return m_particle_scheduler;
}

core::HierarchicalTimeBinScheduler&
RungZeroTimeState::gasCellScheduler() noexcept {
  return m_gas_cell_scheduler;
}

const core::HierarchicalTimeBinScheduler&
RungZeroTimeState::gasCellScheduler() const noexcept {
  return m_gas_cell_scheduler;
}

core::IntegratorState& RungZeroTimeState::integratorState() noexcept {
  return m_integrator_state;
}

const core::IntegratorState& RungZeroTimeState::integratorState() const noexcept {
  return m_integrator_state;
}

internal::PendingOutputBoundary& RungZeroTimeState::pendingOutput() noexcept {
  return m_pending_output;
}

const internal::PendingOutputBoundary&
RungZeroTimeState::pendingOutput() const noexcept {
  return m_pending_output;
}

RungZeroTimeState initializeRungZeroTimeState(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    core::SimulationState& state,
    const core::UnitSystem& units,
    const core::LambdaCdmBackground* cosmology_background,
    const io::RestartReadResult* restart_state) {
  if (config.numerics.hierarchical_max_rung > 0 &&
      (state.cells.size() != 0U || state.star_particles.size() != 0U || state.black_holes.size() != 0U || state.tracers.size() != 0U)) {
    throw std::logic_error("hierarchical KDK currently requires a zero-gas DMO population");
  }
  const std::uint8_t max_bin = static_cast<std::uint8_t>(
      std::max(0, std::min(config.numerics.hierarchical_max_rung, 12)));
  RungZeroTimeState time_state(max_bin);
  if (restart_state != nullptr) {
    time_state.m_particle_scheduler.importPersistentState(
        restart_state->scheduler_state);
    time_state.m_gas_cell_scheduler.importPersistentState(
        restart_state->gas_cell_scheduler_state);
    time_state.m_integrator_state = restart_state->integrator_state;
    if (restart_state->diagnostics.restart_schema_version <
        io::restartSchema().version) {
      time_state.m_integrator_state.pm_refresh_enabled = true;
    }
    if (time_state.m_particle_scheduler.elementCount() !=
        state.particles.size()) {
      throw std::runtime_error(
          "ReferenceWorkflow restart payload particle scheduler coverage does not match SimulationState");
    }
    if (time_state.m_gas_cell_scheduler.elementCount() != state.cells.size()) {
      throw std::runtime_error(
          "ReferenceWorkflow restart payload gas-cell scheduler coverage does not match SimulationState");
    }
    state.requireGasCellIdentityMapCoversDenseRows(
        "ReferenceWorkflow restart resume");
  } else {
    initializeSchedulerBins(
        state,
        time_state.m_particle_scheduler,
        time_state.m_gas_cell_scheduler);
    if (!options.initial_particle_scheduler_identity_records.empty()) {
      core::rebuildSchedulerFromParticleIdentityRecords(
          time_state.m_particle_scheduler,
          options.initial_particle_scheduler_identity_records,
          state.particle_sidecar.particle_id);
    }
    core::IntegratorState& integrator_state = time_state.m_integrator_state;
    integrator_state.step_index = options.step_index;
    integrator_state.current_time_code = config.numerics.t_code_begin;
    integrator_state.time_si_per_code = units.timeSiPerCode();
    integrator_state.current_scale_factor = cosmology_background != nullptr
        ? config.numerics.a_begin
        : 1.0;
    integrator_state.current_redshift = integrator_state.current_scale_factor > 0.0
        ? 1.0 / integrator_state.current_scale_factor - 1.0
        : 0.0;
    integrator_state.current_hubble_rate_code = cosmology_background != nullptr
        ? cosmology_background->hubbleSi(
              integrator_state.current_scale_factor) *
              integrator_state.time_si_per_code
        : 0.0;
    integrator_state.dt_time_code = options.dt_time_code > 0.0
        ? options.dt_time_code
        : (config.numerics.t_code_end - config.numerics.t_code_begin);
    integrator_state.time_bins.hierarchical_enabled = true;
    integrator_state.time_bins.max_bin =
        time_state.m_particle_scheduler.maxBin();
    integrator_state.pm_refresh_enabled = true;
    // Multirate outer-PM cycles refresh at every cycle-closing opportunity:
    // the cycle length lives in the hierarchical driver, not in the sync state.
    integrator_state.pm_sync_state.reset(
        config.numerics.treepm_pm_reuse_policy == core::TreePmPmReusePolicy::kMultirateOuterKick ? 1U :
        static_cast<std::uint64_t>(
        std::max(config.numerics.treepm_update_cadence_steps, 1)));
    state.updateAllParticleDriftEpoch(integrator_state.current_time_code,
        integrator_state.current_scale_factor);
  }
  syncTimeBinsFromSchedulers(
      time_state.m_particle_scheduler,
      time_state.m_gas_cell_scheduler,
      state);
  const io::OutputCadencePersistentState* restored_output =
      restart_state != nullptr ? &restart_state->output_cadence_state : nullptr;
  time_state.m_pending_output = internal::initializeOutputCadence(
      config,
      options,
      time_state.m_integrator_state,
      restored_output);
  return time_state;
}

TimeCoordinator::TimeCoordinator(
    const RuntimeServices& services,
    RungZeroTimeState& time_state,
    GravityRuntime& gravity,
    HydroAmrRuntime& hydro_amr,
    SourceRuntime& source,
    RuntimeExecutionPlan execution_plan,
    const internal::MigrationBalanceRuntime& migration_balance,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table) noexcept
    : m_services(services),
      m_time_state(time_state),
      m_gravity(gravity),
      m_hydro_amr(hydro_amr),
      m_source(source),
      m_execution_plan(std::move(execution_plan)),
      m_migration_balance(migration_balance),
      m_effective_eos_table(std::move(effective_eos_table)) {}

void TimeCoordinator::runRungZeroSegment(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    core::SimulationState& state,
    const core::LambdaCdmBackground* cosmology_background,
    std::vector<std::uint64_t>& expected_global_particle_ids,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    const core::ModePolicy& mode_policy,
    bool restoring_from_restart) {
  if (config.numerics.hierarchical_max_rung > 0) {
    runHierarchicalSegment(config, options, state, cosmology_background,
        expected_global_particle_ids, report, profiler, mode_policy);
    return;
  }
  core::HierarchicalTimeBinScheduler& particle_scheduler =
      m_time_state.m_particle_scheduler;
  core::HierarchicalTimeBinScheduler& gas_cell_scheduler =
      m_time_state.m_gas_cell_scheduler;
  core::IntegratorState& integrator_state = m_time_state.m_integrator_state;
  internal::PendingOutputBoundary& pending_output = m_time_state.m_pending_output;
  if (!restoring_from_restart) {
    static_cast<void>(updateAdaptiveTimeBins(
        state,
        particle_scheduler,
        gas_cell_scheduler,
        integrator_state,
        config,
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        true));
  }
  ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);

  const auto install_authoritative_domain_geometry = [&]() {
    std::exception_ptr install_failure;
    try {
      core::MemoryReservation geometry_reservation;
      {
        const auto leaves = m_migration_balance.authoritativeTopDomainLeaves(
            state, m_gravity.decompositionEpoch(), &geometry_reservation);
        m_gravity.installAuthoritativeTopDomainLeaves(
            leaves, state.gravitySourceGeneration());
      }
      if (geometry_reservation.valid()) {
        const std::array geometry_reports{
            core::collectSimulationMemoryReport(state),
            core::collectSchedulerMemoryReport(
                particle_scheduler, gas_cell_scheduler),
            m_gravity.memoryReport(),
            m_hydro_amr.memoryReport(),
            m_source.memoryReport(),
            m_services.profiler.retainedEventMemoryReport()};
        geometry_reservation.reconcileBaselineOwnedAndRelease(
            core::memoryReportBaselineOwnedBytes(
                core::mergeMemoryReports(geometry_reports)));
      }
    } catch (...) {
      install_failure = std::current_exception();
    }
    FailureCoordinator(m_services).rethrowCollectiveFailure(
        install_failure, "authoritative top-domain installation");
  };
  // Domain geometry is a decomposition-owned contract. Install at the segment
  // boundary (including restart) and again only after ownership decomposition
  // changes; each install replaces both the stable seed leaf set and the
  // current published set. commitParticleDecompositionChange() explicitly
  // invalidates geometry freshness first, so stale routing geometry is never
  // consumable between commit and reinstall. Between those events
  // GravityRuntime refits published leaf bounds from the unchanged seed set
  // against the current source generation before each force solve; TreePM
  // validates freshness (generation equality) and the existing conservative
  // coverage check before routing.
  install_authoritative_domain_geometry();

  const std::uint64_t run_start_step_index = integrator_state.step_index;
  const std::uint64_t configured_segment_steps = options.max_steps_override > 0
      ? options.max_steps_override
      : static_cast<std::uint64_t>(std::max(config.numerics.max_global_steps, 0));
  const std::uint64_t target_step_index = integrator_state.step_index + configured_segment_steps;
  core::TransientStepWorkspace workspace(m_services.memory_governor);
  // CHUI Part 3 Checkpoint D: collective global-step rollback (rung zero).
  // With numerics.global_step_max_retries > 0, the complete authoritative
  // step-start state (SimulationState incl. composition/star/module
  // sidecars and AMR hierarchy, integrator, both schedulers, pending output,
  // gravity endpoint force) is captured under MemoryGovernor admission. A
  // stage that rejects with kNumericalRetry on any rank (collective verdict)
  // restores everything and retries the SAME physical model with half the
  // interval (legal retry). Any other rejection kind (model downgrade,
  // resource refusal, invalid state) is an illegal downgrade request and is
  // rethrown on every rank. Event chronology is restored with the population
  // sidecar, so retried events are never duplicated. Output is written only
  // after an accepted step (the pending-output latch is restored too).
  // Cost: one full copy of the canonical state per step while enabled (an
  // explicit, governed reference mode; bounded delta journals are future
  // optimization work). Default 0 keeps the historical behaviour.
  struct GlobalStepRetryPoint {
    std::optional<core::SimulationState> state;
    core::IntegratorState integrator{};
    std::optional<core::HierarchicalTimeBinScheduler> particles;
    std::optional<core::HierarchicalTimeBinScheduler> cells;
    internal::PendingOutputBoundary pending{};
    TransactionalForceEndpoint force_endpoint{};
    std::uint64_t layout_certificate = 0U;
    core::MemoryReservation reservation;
  };
  const std::uint32_t global_step_max_retries =
      static_cast<std::uint32_t>(std::max(config.numerics.global_step_max_retries, 0));
  std::uint32_t global_step_retries = 0U;
  double global_step_retry_dt_limit = 0.0;
  const auto layoutCertificate = [&state]() {
    std::uint64_t hash = 1469598103934665603ULL;
    const auto mix = [&hash](std::uint64_t v) { hash ^= v; hash *= 1099511628211ULL; };
    mix(state.particles.size());
    for (const std::uint64_t id : state.particle_sidecar.particle_id) mix(id);
    mix(state.cells.size());
    return hash;
  };
  while (integrator_state.step_index < target_step_index &&
         integrator_state.current_time_code < config.numerics.t_code_end) {
    const RuntimeConsoleReporter::Clock::time_point console_step_begin =
        RuntimeConsoleReporter::Clock::now();
    const std::uint64_t pm_refresh_count_before = m_gravity.longRangeRefreshCount();
    const std::uint64_t pm_reuse_count_before = m_gravity.longRangeReuseCount();
    GlobalStepRetryPoint retry_point;
    if (global_step_max_retries > 0U) {
      std::exception_ptr capture_failure;
      try {
        if (m_services.memory_governor != nullptr) {
          retry_point.reservation = m_services.memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state)),
              "time.global_step_rollback.state_copy");
          retry_point.reservation.commit();
        }
        retry_point.state.emplace(state);
        retry_point.integrator = integrator_state;
        retry_point.particles.emplace(particle_scheduler);
        retry_point.cells.emplace(gas_cell_scheduler);
        retry_point.pending = pending_output;
        retry_point.layout_certificate = layoutCertificate();
        m_gravity.captureTransactionalForceEndpoint(state, retry_point.layout_certificate, retry_point.force_endpoint);
        m_gravity.beginCycleTransaction();
      } catch (...) {
        capture_failure = std::current_exception();
      }
      FailureCoordinator(m_services).rethrowCollectiveFailure(capture_failure, "global step rollback capture");
    }
    // Rung zero has one physical timestep authority. Re-evaluate all local
    // criteria before constructing KDK stage times, then agree the minimum
    // collectively so every rank advances the identical interval.
    double local_physical_dt = std::numeric_limits<double>::infinity();
    std::exception_ptr local_timestep_failure;
    try {
      local_physical_dt = updateAdaptiveTimeBins(
          state, particle_scheduler, gas_cell_scheduler, integrator_state, config,
          m_gravity.particleAccelX(), m_gravity.particleAccelY(), m_gravity.particleAccelZ(),
          m_gravity.cellAccelX(), m_gravity.cellAccelY(), m_gravity.cellAccelZ(),
          {}, {}, true);
      const bool explicit_no_limit = std::isinf(local_physical_dt) && local_physical_dt > 0.0;
      if ((!std::isfinite(local_physical_dt) && !explicit_no_limit) ||
          local_physical_dt <= 0.0) {
        throw std::runtime_error(
            "local global-timestep candidate is invalid; expected finite positive or +infinity");
      }
    } catch (...) {
      local_timestep_failure = std::current_exception();
    }
    m_services.mpi_context.rethrowCollectivePreparationFailure(
        local_timestep_failure, "global timestep candidate selection");
    double accepted_dt = m_services.mpi_context.allreduceMinDouble(local_physical_dt);
    {
      // G1.2: record the limiting criterion and AMR level of this decision
      // (the rank whose candidate equals the global minimum is the limiter).
      const TimestepLimiterRecord& limiter = m_last_timestep_limiter;
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "time.timestep_limiter",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.time",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "local timestep limiter for the agreed global candidate",
          .payload = {{"criterion", limiter.criterion},
                      {"gas_cell", limiter.gas_cell ? "true" : "false"},
                      {"element_row", std::to_string(limiter.element)},
                      {"amr_level_offset", std::to_string(limiter.amr_level)},
                      {"hydro_cfl_level_relaxation", formatRuntimeDouble(limiter.hydro_cfl_level_relaxation)},
                      {"local_candidate_dt_time_code", formatRuntimeDouble(local_physical_dt)},
                      {"global_candidate_dt_time_code", formatRuntimeDouble(accepted_dt)},
                      {"is_global_limiter", local_physical_dt == accepted_dt ? "true" : "false"}},
      });
    }
    if (options.dt_time_code > 0.0) accepted_dt = std::min(accepted_dt, options.dt_time_code);
    const double step_begin_scale_factor = integrator_state.current_scale_factor;
    const double remaining_time_code =
        config.numerics.t_code_end - integrator_state.current_time_code;
    if (!std::isfinite(remaining_time_code) || remaining_time_code <= 0.0) {
      throw std::runtime_error(
          "ReferenceWorkflow computed an invalid remaining endpoint interval");
    }
    double ordered_step_limit_time_code = config.numerics.t_code_end;
    bool limited_by_output_event = false;
    if (pending_output.snapshot_interval_time_code > 0.0 &&
        pending_output.next_snapshot_time_code < ordered_step_limit_time_code &&
        pending_output.next_snapshot_time_code > integrator_state.current_time_code) {
      ordered_step_limit_time_code = pending_output.next_snapshot_time_code;
      limited_by_output_event = true;
    }
    const double ordered_remaining_time_code =
        ordered_step_limit_time_code - integrator_state.current_time_code;
    if (!std::isfinite(ordered_remaining_time_code) || ordered_remaining_time_code <= 0.0) {
      throw std::runtime_error("ReferenceWorkflow computed an invalid ordered timeline interval");
    }
    if (accepted_dt > ordered_remaining_time_code) {
      const double unclipped_dt_time_code = accepted_dt;
      accepted_dt = ordered_remaining_time_code;
      if (limited_by_output_event && std::isfinite(unclipped_dt_time_code)) {
        pending_output.restart_resume_dt_time_code = unclipped_dt_time_code;
      }
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = limited_by_output_event ? "time.output_event_clip" : "time.endpoint_clip",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.time",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = limited_by_output_event
              ? "integration interval clipped to an ordered output event"
              : "integration interval clipped to the configured endpoint",
          .payload = {{"unclipped_dt_time_code", formatRuntimeDouble(unclipped_dt_time_code)},
                      {"clipped_dt_time_code", formatRuntimeDouble(accepted_dt)},
                      {limited_by_output_event ? "output_event_time_code" : "endpoint_time_code",
                       formatRuntimeDouble(ordered_step_limit_time_code)}},
      });
    }
    if (!std::isfinite(accepted_dt) || accepted_dt <= 0.0) {
      throw std::runtime_error("ReferenceWorkflow global timestep selection produced an invalid dt");
    }
    if (global_step_retry_dt_limit > 0.0 && accepted_dt > global_step_retry_dt_limit) {
      // Legal retry: same physics, shorter interval (identical on all ranks).
      accepted_dt = global_step_retry_dt_limit;
      profiler.counters().addCount("time.global_step_retry_shortened_intervals");
    }
    integrator_state.dt_time_code = accepted_dt;
    const double accepted_dt_time_code_for_console = accepted_dt;

    const std::span<const std::uint32_t> active_particles =
        particle_scheduler.beginSubstep();
    const std::span<const std::uint32_t> active_cells =
        gas_cell_scheduler.beginSubstep();
    if (particle_scheduler.currentTick() != gas_cell_scheduler.currentTick()) {
      throw std::runtime_error("particle and gas-cell schedulers lost their shared integer timeline");
    }
    const bool local_has_active_work = !active_particles.empty() || !active_cells.empty();
    const std::uint64_t active_work_rank_count =
        m_services.mpi_context.allreduceSumUint64(local_has_active_work ? 1ULL : 0ULL);
    if (active_work_rank_count == 0ULL) {
      particle_scheduler.endSubstep();
      gas_cell_scheduler.endSubstep();
      const bool particle_decomposition_changed = m_migration_balance.rebalance(
          state,
          particle_scheduler,
          gas_cell_scheduler,
          parallel::DecompositionRuntimeMeasurements{},
          active_particles,
          expected_global_particle_ids,
          integrator_state.step_index);
      if (particle_decomposition_changed) {
        m_gravity.commitParticleDecompositionChange();
        install_authoritative_domain_geometry();
      }
      ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
      syncTimeBinsFromSchedulers(particle_scheduler, gas_cell_scheduler, state);
      continue;
    }

    workspace.clear();
    profiler.counters().addCount("workflow_workspace_reuses", 1U);
    profiler.counters().setCount("scheduler_active_index_copy_bytes", 0U);
    internal::latchOutputRequestForCompletedStep(
        config,
        options,
        integrator_state.step_index + 1U,
        integrator_state.current_time_code + integrator_state.dt_time_code,
        pending_output);
    const double resume_dt_after_step = pending_output.restart_resume_dt_time_code;
    const core::StepBoundaryKind requested_boundary =
        internal::requestedBoundaryForPendingOutput(pending_output);
    const bool console_snapshot_requested = pending_output.snapshot_due;
    const bool console_restart_requested = pending_output.checkpoint_due;
    // One fused reduction for the four population counts (previously four
    // sequential collectives per step); values and semantics are unchanged.
    std::array<std::uint64_t, 4U> global_population_counts{
        static_cast<std::uint64_t>(active_particles.size()),
        static_cast<std::uint64_t>(state.particles.size()),
        static_cast<std::uint64_t>(active_cells.size()),
        static_cast<std::uint64_t>(state.cells.size())};
    m_services.mpi_context.allreduceSumUint64sInPlace(global_population_counts);
    const std::uint64_t global_active_particle_count = global_population_counts[0];
    const std::uint64_t global_particle_count = global_population_counts[1];
    const std::uint64_t global_active_cell_count = global_population_counts[2];
    const std::uint64_t global_cell_count = global_population_counts[3];
    core::ActiveSetDescriptor active_set = core::makeSchedulerActiveSetDescriptor(
        particle_scheduler, state, active_particles, active_cells);
    active_set.has_global_synchronization_metadata = true;
    active_set.globally_complete_active_set =
        global_active_particle_count == global_particle_count &&
        global_active_cell_count == global_cell_count;
    m_newly_created_particle_ids.clear();
    if (global_step_max_retries == 0U) {
      executeSingleStep(
          state, integrator_state, active_set, cosmology_background, &workspace, &mode_policy, &profiler,
          particle_scheduler.currentTick(), requested_boundary);
    } else {
      bool rejected = false;
      std::uint64_t illegal = 0U;
      std::string rejection_reason;
      try {
        executeSingleStep(
            state, integrator_state, active_set, cosmology_background, &workspace, &mode_policy, &profiler,
            particle_scheduler.currentTick(), requested_boundary);
      } catch (const core::StageTransactionRejected& rejection) {
        rejected = true;
        rejection_reason = rejection.what();
        illegal = (rejection.kind() == core::StageRejectionKind::kNumericalRetry ||
                   rejection.kind() == core::StageRejectionKind::kPeerRejected) ? 0U : 1U;
      }
      // Stage rejections are collective by contract (StageTransaction with
      // the MPI decision); the verdict on legality is reduced so every rank
      // takes the same branch.
      const std::uint64_t any_rejected = m_services.mpi_context.allreduceMaxUint64(rejected ? 1U : 0U);
      const std::uint64_t any_illegal = m_services.mpi_context.allreduceMaxUint64(illegal);
      if (any_rejected != 0U) {
        if (any_illegal != 0U || global_step_retries >= global_step_max_retries) {
          throw std::runtime_error(
              "global step rejected and not retryable (" +
              std::string(any_illegal != 0U ? "illegal model downgrade / resource / invalid-state rejection"
                                             : "retry limit exhausted") +
              "): " + rejection_reason);
        }
        std::exception_ptr restore_failure;
        try {
          const std::uint64_t generation_reached = state.gravitySourceGeneration();
          state = std::move(*retry_point.state);
          while (state.gravitySourceGeneration() <= generation_reached) state.bumpGravitySourceGeneration();
          state.bumpCellIndexGeneration();
          integrator_state = retry_point.integrator;
          particle_scheduler = std::move(*retry_point.particles);
          gas_cell_scheduler = std::move(*retry_point.cells);
          pending_output = retry_point.pending;
          m_gravity.rollbackCycleTransaction(retry_point.force_endpoint, state, retry_point.layout_certificate);
          integrator_state.pm_source_generation = state.gravitySourceGeneration();
          syncTimeBinsFromSchedulers(particle_scheduler, gas_cell_scheduler, state);
        } catch (...) {
          restore_failure = std::current_exception();
        }
        FailureCoordinator(m_services).rethrowCollectiveFailure(restore_failure, "global step rollback restore");
        ++global_step_retries;
        global_step_retry_dt_limit = 0.5 * accepted_dt;
        profiler.counters().addCount("time.global_step_rollbacks");
        continue;
      }
      if (m_gravity.cycleTransactionOpen()) m_gravity.commitCycleTransaction();
      global_step_retries = 0U;
      global_step_retry_dt_limit = 0.0;
    }
    if (cosmology_background != nullptr) {
      const double a_next = integrator_state.current_scale_factor;
      if (!std::isfinite(a_next) || a_next <= 0.0 || a_next + 1.0e-14 < step_begin_scale_factor) {
        throw std::runtime_error("cosmological step violated finite positive monotonic scale-factor invariant");
      }
      const double delta_ln_a = std::log(a_next / step_begin_scale_factor);
      const double allowed_delta = config.numerics.cosmology_max_delta_ln_a * (1.0 + 1.0e-10);
      if (!std::isfinite(delta_ln_a) || delta_ln_a > allowed_delta) {
        throw std::runtime_error("cosmological step exceeded numerics.cosmology_max_delta_ln_a");
      }
      if ((config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kScaleFactor ||
           config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kLogScaleFactor) &&
          a_next > config.numerics.a_end + 1.0e-12 * std::max(1.0, config.numerics.a_end)) {
        throw std::runtime_error("cosmological step crossed the authoritative scale-factor endpoint");
      }
    }
    expected_global_particle_ids.insert(
        expected_global_particle_ids.end(),
        m_newly_created_particle_ids.begin(),
        m_newly_created_particle_ids.end());
    if (resume_dt_after_step > 0.0) {
      integrator_state.dt_time_code = resume_dt_after_step;
    }

    const std::array runtime_reports{
        core::collectSimulationMemoryReport(state, &workspace),
        core::collectSchedulerMemoryReport(
            particle_scheduler, gas_cell_scheduler),
        m_gravity.memoryReport(),
        m_hydro_amr.memoryReport(),
        m_source.memoryReport(),
        profiler.retainedEventMemoryReport()};
    core::MemoryReport merged_runtime_memory_report =
        core::mergeMemoryReports(runtime_reports);
    if (m_services.memory_governor != nullptr) {
      const std::uint64_t governor_baseline_bytes =
          core::memoryReportBaselineOwnedBytes(merged_runtime_memory_report);
      m_services.memory_governor->setBaselineOwnedBytes(governor_baseline_bytes);
      core::attachMemoryGovernorSnapshot(
          merged_runtime_memory_report, *m_services.memory_governor);
    }
    core::attachProcessMemoryObservation(
        merged_runtime_memory_report, core::observeProcessMemory());
    attachDistributedMemoryTelemetry(
        merged_runtime_memory_report, m_services);
    profiler.setMemoryReport(std::move(merged_runtime_memory_report));
    state.metadata.step_index = integrator_state.step_index;
    state.metadata.scale_factor = integrator_state.current_scale_factor;
    ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
    static_cast<void>(updateAdaptiveTimeBins(
        state,
        particle_scheduler,
        gas_cell_scheduler,
        integrator_state,
        config,
        m_gravity.particleAccelX(),
        m_gravity.particleAccelY(),
        m_gravity.particleAccelZ(),
        m_gravity.cellAccelX(),
        m_gravity.cellAccelY(),
        m_gravity.cellAccelZ(),
        active_particles,
        active_cells,
        false));
    profiler.counters().addCount(
        "timestep_particle_criteria_evaluations",
        static_cast<std::uint64_t>(active_particles.size()));
    profiler.counters().addCount(
        "timestep_gas_cell_criteria_evaluations",
        static_cast<std::uint64_t>(active_cells.size()));
    ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
    particle_scheduler.endSubstep();
    gas_cell_scheduler.endSubstep();

    parallel::DecompositionRuntimeMeasurements rebalance_measurements =
        m_gravity.lastRuntimeDecompositionMeasurements();
    const hydro::HydroProfileEvent& hydro_profile = m_hydro_amr.lastHydroProfile();
    const std::uint64_t hydro_ghost_exchange_bytes =
        m_hydro_amr.ghostExchangeBytesRecent();
    rebalance_measurements.hydro_face_fluxes_recent = hydro_profile.face_count;
    rebalance_measurements.hydro_wall_ms_recent = hydro_profile.total_ms;
    rebalance_measurements.ghost_exchange_bytes_recent += hydro_ghost_exchange_bytes;
    rebalance_measurements.has_measurements = rebalance_measurements.has_measurements ||
        hydro_profile.face_count > 0 || hydro_profile.total_ms > 0.0 ||
        hydro_ghost_exchange_bytes > 0;
    const bool particle_decomposition_changed = m_migration_balance.rebalance(
        state,
        particle_scheduler,
        gas_cell_scheduler,
        rebalance_measurements,
        active_particles,
        expected_global_particle_ids,
        integrator_state.step_index);
    if (particle_decomposition_changed) {
      m_gravity.commitParticleDecompositionChange();
      install_authoritative_domain_geometry();
      if (m_services.console_reporter != nullptr) {
        m_services.console_reporter->emitDecomposition(
            integrator_state.step_index, m_gravity.decompositionEpoch());
      }
    }
    ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
    syncTimeBinsFromSchedulers(particle_scheduler, gas_cell_scheduler, state);
    if (integrator_state.last_completed_restart_safe &&
        integrator_state.last_completed_boundary_kind !=
            core::StepBoundaryKind::kLocalActiveBinStep) {
      executeOutputBoundary(
          state,
          integrator_state,
          &profiler,
          requested_boundary);
    }

    if (m_services.console_reporter != nullptr) {
      RuntimeConsoleStepStatus console_status;
      console_status.step_index = integrator_state.step_index;
      console_status.t_code = integrator_state.current_time_code;
      console_status.dt_time_code = accepted_dt_time_code_for_console;
      if (mode_policy.cosmological_comoving_frame) {
        console_status.a_scale = integrator_state.current_scale_factor;
        if (std::isfinite(integrator_state.current_scale_factor) &&
            integrator_state.current_scale_factor > 0.0) {
          console_status.redshift =
              1.0 / integrator_state.current_scale_factor - 1.0;
        }
      }
      console_status.active_particle_count = global_active_particle_count;
      console_status.total_particle_count = global_particle_count;
      console_status.active_cell_count = global_active_cell_count;
      console_status.total_cell_count = global_cell_count;
      console_status.wall_step_seconds = std::chrono::duration<double>(
          RuntimeConsoleReporter::Clock::now() - console_step_begin).count();
      const bool pm_refreshed =
          m_gravity.longRangeRefreshCount() > pm_refresh_count_before;
      const bool pm_reused =
          m_gravity.longRangeReuseCount() > pm_reuse_count_before;
      console_status.pm_activity = pm_refreshed
          ? "refresh"
          : (pm_reused ? "reuse" : "none");
      const bool snapshot_flushed =
          console_snapshot_requested && !pending_output.snapshot_due;
      const bool restart_flushed =
          console_restart_requested && !pending_output.checkpoint_due;
      if (snapshot_flushed && restart_flushed) {
        console_status.output_activity = "snapshot+restart";
      } else if (snapshot_flushed) {
        console_status.output_activity = "snapshot";
      } else if (restart_flushed) {
        console_status.output_activity = "restart";
      } else {
        console_status.output_activity = "none";
      }
      if (const core::MemoryReport* memory_report = profiler.memoryReport();
          memory_report != nullptr && memory_report->governor_snapshot.has_value()) {
        console_status.memory_headroom_bytes =
            memory_report->governor_snapshot->headroom_bytes;
        console_status.memory_pressure = std::string(core::memoryPressureLabel(
            memory_report->governor_snapshot->pressure));
      }
      m_services.console_reporter->emitStep(console_status);
    }
  }

  report.completed_steps = integrator_state.step_index - run_start_step_index;
  report.final_time_code = integrator_state.current_time_code;
  report.final_scale_factor = integrator_state.current_scale_factor;
  report.final_hydro_cfl_diagnostics = m_hydro_amr.lastHydroCflDiagnostics();
  report.final_hydro_imported_mpi_ghosts =
      static_cast<std::uint64_t>(m_hydro_amr.remoteImportedGhostCount());
  report.final_hydro_remote_interface_faces =
      static_cast<std::uint64_t>(m_hydro_amr.remoteInterfaceFaceCount());
  report.final_hydro_remote_stale_payloads =
      static_cast<std::uint64_t>(m_hydro_amr.remoteStaleInvalidPayloadCount());
  report.final_hydro_directed_amr_remote_interfaces =
      m_hydro_amr.directedAmrRemoteInterfaceCount();
  report.final_hydro_directed_amr_patch_cell_records_received =
      m_hydro_amr.directedAmrPatchCellRecordsReceived();
}

void TimeCoordinator::runHierarchicalSegment(
    const core::SimulationConfig& config, const ReferenceWorkflowOptions& options,
    core::SimulationState& state, const core::LambdaCdmBackground* background,
    std::vector<std::uint64_t>& expected_ids, ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler, const core::ModePolicy& mode_policy) {
  auto& particles = m_time_state.m_particle_scheduler;
  auto& cells = m_time_state.m_gas_cell_scheduler;
  auto& integrator = m_time_state.m_integrator_state;
  auto& pending = m_time_state.m_pending_output;
  const auto& mpi = m_services.mpi_context;
  std::exception_ptr eligibility_failure;
  try {
    if (state.cells.size() != 0U || state.star_particles.size() != 0U ||
        state.black_holes.size() != 0U || state.tracers.size() != 0U ||
        !state.hasHomogeneousDmoMetadata() || particles.maxBin() == 0U ||
        particles.maxBin() != cells.maxBin()) {
      throw std::logic_error("hierarchical workflow requires compact collisionless DMO on one scheduler hierarchy");
    }
  } catch (...) { eligibility_failure = std::current_exception(); }
  FailureCoordinator(m_services).rethrowCollectiveFailure(eligibility_failure, "hierarchical DMO eligibility");
  core::TransientStepWorkspace workspace(m_services.memory_governor);
  const auto owned_runtime_bytes = [&]() {
    const std::array reports{core::collectSimulationMemoryReport(state, &workspace),
        core::collectSchedulerMemoryReport(particles, cells), m_gravity.memoryReport(),
        m_hydro_amr.memoryReport(), m_source.memoryReport(), profiler.retainedEventMemoryReport()};
    return core::memoryReportBaselineOwnedBytes(core::mergeMemoryReports(reports));
  };
  const auto reconcile_runtime = [&]() {
    core::MemoryReport memory;
    std::exception_ptr memory_failure;
    try {
      const std::array reports{core::collectSimulationMemoryReport(state, &workspace),
          core::collectSchedulerMemoryReport(particles, cells), m_gravity.memoryReport(),
          m_hydro_amr.memoryReport(), m_source.memoryReport(), profiler.retainedEventMemoryReport()};
      memory = core::mergeMemoryReports(reports);
      if (m_services.memory_governor != nullptr) {
        m_services.memory_governor->setBaselineOwnedBytes(core::memoryReportBaselineOwnedBytes(memory));
        core::attachMemoryGovernorSnapshot(memory, *m_services.memory_governor);
      }
      core::attachProcessMemoryObservation(memory, core::observeProcessMemory());
    } catch (...) { memory_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(memory_failure, "hierarchical local memory report preparation");
    attachDistributedMemoryTelemetry(memory, m_services);
    profiler.setMemoryReport(std::move(memory));
  };
  const auto install_geometry = [&]() {
    std::exception_ptr failure;
    try {
      core::MemoryReservation reservation;
      {
        const auto leaves = m_migration_balance.authoritativeTopDomainLeaves(state, m_gravity.decompositionEpoch(), &reservation);
        m_gravity.installAuthoritativeTopDomainLeaves(leaves, state.gravitySourceGeneration());
      }
      if (reservation.valid()) reservation.reconcileBaselineOwnedAndRelease(owned_runtime_bytes());
    } catch (...) { failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(failure, "hierarchical domain geometry installation");
  };
  install_geometry();
  const std::uint64_t requested_steps = options.max_steps_override > 0U
      ? options.max_steps_override : static_cast<std::uint64_t>(config.numerics.max_global_steps);
  if (integrator.step_index > std::numeric_limits<std::uint64_t>::max() - requested_steps) {
    throw std::overflow_error("hierarchical coarse segment step count overflow");
  }
  const std::uint64_t run_start_step = integrator.step_index;
  const std::uint64_t final_step = integrator.step_index + requested_steps;
  // Multirate outer PM cycle state (TreePmPmReusePolicy::kMultirateOuterKick;
  // core::HierarchicalPmCyclePlan; contract in docs/time_integration.md).
  //
  // * The OUTER interval T (PM half-kick lengths) is fixed when the cycle
  //   opens, from physical criteria only: target blocks x the opening block
  //   interval, the PM displacement criterion, and output/end boundaries.
  //   The external step budget never shortens a cycle (segmentation
  //   invariance); a budget reached inside a cycle is extended to the cycle
  //   end (soft stop) and the executed step count is reported.
  // * INNER blocks are re-planned at every block from the reference criteria
  //   (particle gravity/encounter/growth, expansion, explicit dt) and clipped
  //   to the cycle end, so no block ever uses a quantum coarser than the
  //   current global particle criterion.
  // * A PM displacement requirement discovered inside an open cycle rejects
  //   the cycle: the authoritative state captured at the cycle opening is
  //   restored transactionally on every rank and the cycle is retried with a
  //   shorter T. A single-block cycle cannot be rejected, so retries
  //   terminate; exceeding the retry bound fails closed.
  // * Checkpoints/output occur only at cycle ends; a segment always starts
  //   and ends at a closed cycle, so no cycle state is persisted.
  const bool multirate_pm =
      config.numerics.treepm_pm_reuse_policy == core::TreePmPmReusePolicy::kMultirateOuterKick;
  const std::uint64_t multirate_blocks_per_cycle = multirate_pm
      ? static_cast<std::uint64_t>(std::max(config.numerics.treepm_update_cadence_steps, 1)) : 1U;
  constexpr unsigned k_multirate_max_retries_per_cycle = 64U;
  // Rollback store (multirate PM cycle transactions). Hierarchical KDK runs
  // homogeneous DMO with frozen ownership inside a cycle (no rebalance or
  // migration at interior boundaries), so a block can only mutate: the
  // phase-space lanes, the SFC key lane, the common drift epoch, the
  // metadata, the gravity source generation, the integrator, both
  // schedulers, the pending output boundary and gravity-derived state.
  // Authoritative state among these is journaled here; gravity-derived state
  // is restored (endpoint force) or invalidated by
  // GravityRuntime::rollbackCycleTransaction. Everything else (identity,
  // mass, species/flags/rank, softening, index generation, cell/gas
  // containers) must be unchanged at rollback: a fingerprint taken at the
  // opening verifies it and the rollback fails closed otherwise.
  // The store lives for the whole segment: its lanes keep their capacity and
  // its governed reservation is grown (never re-reserved per cycle) so
  // successive cycles reuse one allocation; both are released at segment end.
  struct MultirateRollbackStore {
    // Memory model (docs/time_integration.md §5): every buffer below is owned
    // here and covered by `reservation`, sized from the actual retained
    // capacities (not a fixed per-particle constant).
    //   per particle row: 6 phase lanes (8 B) + 10 endpoint lanes (8 B)
    //                     + scheduler bin/next/active/pending (1+8+1+1 B)
    //   per SFC-key row:  8 B; per gas-cell scheduler row: 11 B (0 for DMO)
    //   fixed: StateMetadata string capacities + a conservative allocator
    //          overhead of k_multirate_rollback_allocation_overhead_bytes per vector.
    std::vector<double> position_x_comoving;
    std::vector<double> position_y_comoving;
    std::vector<double> position_z_comoving;
    std::vector<double> velocity_x_peculiar;
    std::vector<double> velocity_y_peculiar;
    std::vector<double> velocity_z_peculiar;
    std::vector<std::uint64_t> sfc_key;
    core::StateMetadata metadata{};
    double drift_time_code = 0.0;
    double drift_scale_factor = 0.0;
    std::uint64_t frozen_fingerprint = 0U;
    core::IntegratorState integrator{};
    core::TimeBinPersistentState particle_scheduler{};
    core::TimeBinPersistentState cell_scheduler{};
    // Row-aligned cycle-opening force endpoint (restored without ID remap).
    TransactionalForceEndpoint force_endpoint{};
    internal::PendingOutputBoundary pending{};
    bool armed = false;  // holds the rollback point of the open cycle
    core::MemoryReservation reservation;
    std::uint64_t reserved_bytes = 0U;
    std::uint64_t peak_reserved_bytes = 0U;
    std::uint64_t peak_retained_capacity_bytes = 0U;
    std::uint64_t growths = 0U;

    [[nodiscard]] static std::uint64_t schedulerCapacityBytes(const core::TimeBinPersistentState& scheduler) noexcept {
      return static_cast<std::uint64_t>(scheduler.bin_index.capacity()) +
          static_cast<std::uint64_t>(scheduler.next_activation_tick.capacity()) * sizeof(std::uint64_t) +
          static_cast<std::uint64_t>(scheduler.active_flag.capacity()) +
          static_cast<std::uint64_t>(scheduler.pending_bin_index.capacity());
    }
    [[nodiscard]] static std::uint64_t metadataBytes(const core::StateMetadata& value) noexcept {
      std::uint64_t bytes = 0U;
      for (const auto* text : {&value.run_name, &value.normalized_config_hash_hex, &value.snapshot_stem,
                               &value.restart_stem}) {
        bytes += static_cast<std::uint64_t>(text->capacity()) + 1U;
      }
      return bytes;
    }
    // Bytes actually retained by every owned buffer (+ overhead allowance).
    [[nodiscard]] std::uint64_t retainedCapacityBytes() const noexcept {
      std::uint64_t bytes = k_multirate_rollback_vector_count * k_multirate_rollback_allocation_overhead_bytes + metadataBytes(metadata);
      for (const auto* lane : {&position_x_comoving, &position_y_comoving, &position_z_comoving,
                               &velocity_x_peculiar, &velocity_y_peculiar, &velocity_z_peculiar}) {
        bytes += static_cast<std::uint64_t>(lane->capacity()) * sizeof(double);
      }
      bytes += static_cast<std::uint64_t>(sfc_key.capacity()) * sizeof(std::uint64_t);
      bytes += force_endpoint.retainedCapacityBytes();
      bytes += schedulerCapacityBytes(particle_scheduler) + schedulerCapacityBytes(cell_scheduler);
      return bytes;
    }
    // Governed bytes needed for fresh storage of the given extents (checked).
    [[nodiscard]] static std::uint64_t requiredBytes(
        std::size_t rows, std::size_t sfc_rows, std::size_t cell_rows, const core::StateMetadata& next_metadata) {
      const std::uint64_t per_row = 6U * sizeof(double) +
          TransactionalForceEndpoint::k_lane_count * sizeof(double) + k_multirate_scheduler_row_bytes;
      std::uint64_t bytes = core::checkedSizeMultiply(
          static_cast<std::uint64_t>(rows), per_row, "multirate rollback snapshot particle bytes");
      bytes = core::checkedMemoryBytesAdd(bytes, core::checkedSizeMultiply(
          static_cast<std::uint64_t>(sfc_rows), sizeof(std::uint64_t), "multirate rollback snapshot SFC bytes"),
          "multirate rollback snapshot bytes");
      bytes = core::checkedMemoryBytesAdd(bytes, core::checkedSizeMultiply(
          static_cast<std::uint64_t>(cell_rows), k_multirate_scheduler_row_bytes, "multirate rollback snapshot cell bytes"),
          "multirate rollback snapshot bytes");
      return core::checkedMemoryBytesAdd(bytes,
          k_multirate_rollback_vector_count * k_multirate_rollback_allocation_overhead_bytes + metadataBytes(next_metadata),
          "multirate rollback snapshot bytes");
    }
    // True when a capture of these extents allocates nothing.
    [[nodiscard]] bool coversRows(std::size_t rows, std::size_t sfc_rows, std::size_t cell_rows) const noexcept {
      for (const auto* lane : {&position_x_comoving, &position_y_comoving, &position_z_comoving,
                               &velocity_x_peculiar, &velocity_y_peculiar, &velocity_z_peculiar}) {
        if (lane->capacity() < rows) return false;
      }
      const auto scheduler_covers = [](const core::TimeBinPersistentState& scheduler, std::size_t count) {
        return scheduler.bin_index.capacity() >= count && scheduler.next_activation_tick.capacity() >= count &&
            scheduler.active_flag.capacity() >= count && scheduler.pending_bin_index.capacity() >= count;
      };
      return sfc_key.capacity() >= sfc_rows && force_endpoint.minimumLaneCapacity() >= rows &&
          scheduler_covers(particle_scheduler, rows) && scheduler_covers(cell_scheduler, cell_rows);
    }
    // Allocates every lane to the given extents (after governor admission).
    void reserveRows(std::size_t rows, std::size_t sfc_rows, std::size_t cell_rows) {
      for (auto* lane : {&position_x_comoving, &position_y_comoving, &position_z_comoving,
                         &velocity_x_peculiar, &velocity_y_peculiar, &velocity_z_peculiar}) {
        lane->reserve(rows);
      }
      sfc_key.reserve(sfc_rows);
      force_endpoint.reserveRows(rows);
      for (auto [scheduler, count] : {std::pair{&particle_scheduler, rows}, std::pair{&cell_scheduler, cell_rows}}) {
        scheduler->bin_index.reserve(count);
        scheduler->next_activation_tick.reserve(count);
        scheduler->active_flag.reserve(count);
        scheduler->pending_bin_index.reserve(count);
      }
    }
    void releaseStorage() noexcept {
      armed = false;
      for (auto* lane : {&position_x_comoving, &position_y_comoving, &position_z_comoving,
                         &velocity_x_peculiar, &velocity_y_peculiar, &velocity_z_peculiar}) {
        std::vector<double>().swap(*lane);
      }
      std::vector<std::uint64_t>().swap(sfc_key);
      force_endpoint.release();
      particle_scheduler = core::TimeBinPersistentState{};
      cell_scheduler = core::TimeBinPersistentState{};
      metadata = core::StateMetadata{};
      reservation.release();
      reserved_bytes = 0U;
    }
  };
  struct MultirateCycle {
    bool open = false;
    std::uint64_t index = 0U;            // cycles opened in this segment (diagnostic)
    std::uint64_t completed_blocks = 0U;
    double opening_time_code = 0.0;
    double end_time_code = 0.0;
    double length_time_code = 0.0;
    const char* length_limiter = "target_blocks";
    unsigned retries = 0U;               // rejections of the current opening
    double retry_length_limit = 0.0;     // > 0 after a rejection: bound on T at reopening
    std::uint64_t forced_rejections = 0U;  // test seam bookkeeping (see ReferenceWorkflowOptions)
    // Block evidence of the open attempt; published when the cycle closes,
    // discarded by a rollback (rejected attempts never reach the report).
    std::vector<ReferenceWorkflowReport::HierarchicalBlockRecord> staged_records;
  } pm_cycle;
  MultirateRollbackStore rollback_store;
  // Released on every exit path (normal end, collective failure).
  struct RollbackStoreRelease {
    MultirateRollbackStore& store;
    ~RollbackStoreRelease() { store.releaseStorage(); }
  } rollback_store_release{rollback_store};
  std::uint64_t forced_rollbacks_remaining =
      multirate_pm ? options.multirate_test_forced_cycle_rollbacks : 0U;
  const std::uint64_t forced_rejections_per_cycle =
      std::max<std::uint64_t>(options.multirate_test_forced_rejections_per_cycle, 1U);
  std::uint64_t snapshot_refusals_remaining =
      multirate_pm ? options.multirate_test_snapshot_refusals : 0U;
  const bool snapshot_refusal_this_rank = options.multirate_test_snapshot_refusal_ranks.empty() ||
      std::find(options.multirate_test_snapshot_refusal_ranks.begin(),
                options.multirate_test_snapshot_refusal_ranks.end(),
                m_services.mpi_context.worldRank()) != options.multirate_test_snapshot_refusal_ranks.end();
  core::HierarchicalPmCyclePlan pm_cycle_plan{};
  // Hierarchical scheduler-preparation caches (2026-10-27). Both hold only
  // derived values and are rebuilt from authoritative state on any miss; a
  // hit performs O(1) key comparisons, never an O(N) validation scan.
  // * Identity row list 0..N-1 (the all-active source/target index lane) in
  //   the segment workspace's scratch arena. Key: row count,
  //   SimulationState::particleIndexGeneration() (advanced by every reorder,
  //   resize, compaction, migration commit and population replacement) and
  //   the scratch span itself (data pointer and extent, so a foreign
  //   re-preparation is detected). Force/source generations are not part of
  //   the key: drifting, kicking or rolling back positions never changes row
  //   identity. The arena is reset before a rebuild: the identity list is its
  //   only allocation in this segment, so it stays bounded by one N-row list
  //   (previously every block allocated a fresh list without a reset, so the
  //   arena grew by 4N bytes per block).
  // * Local maximum particle mass (encounter criterion). Key: row count,
  //   SimulationState::particleMassGeneration() (advanced by every dense-row
  //   change and every in-place mass write) and the mass lane address. The
  //   global maximum is still reduced every block (one collective, as before),
  //   so ranks with and without a local hit stay in lock step.
  struct IdentityRowCache {
    bool valid = false;
    std::size_t rows = 0U;
    std::uint64_t particle_index_generation = 0U;
    const std::uint32_t* data = nullptr;
  } identity_rows;
  struct LocalMassMaxCache {
    bool valid = false;
    std::size_t rows = 0U;
    std::uint64_t particle_mass_generation = 0U;
    const double* data = nullptr;
    double value = 0.0;
  } local_mass_max_cache;
  if (multirate_pm && !integrator.last_completed_restart_safe && integrator.step_index > 0U) {
    throw std::logic_error("multirate PM cycle must open at a restart-safe synchronization boundary");
  }
  if (multirate_pm && m_gravity.cycleTransactionOpen()) {
    throw std::logic_error("multirate PM segment started with an open gravity cycle transaction");
  }
  // FNV-1a over everything a cycle must leave unchanged (see above).
  const auto frozen_state_fingerprint = [&]() {
    std::uint64_t hash = 14695981039346656037ULL;
    // Word-wise FNV-1a (8-byte words, then the byte tail): one pass over the
    // frozen lanes at each multi-block opening and each rollback.
    const auto mix_bytes = [&hash](const void* data, std::size_t bytes) {
      const auto* p = static_cast<const unsigned char*>(data);
      std::size_t i = 0U;
      for (; i + sizeof(std::uint64_t) <= bytes; i += sizeof(std::uint64_t)) {
        std::uint64_t word = 0U;
        std::memcpy(&word, p + i, sizeof(word));
        hash ^= word;
        hash *= 1099511628211ULL;
      }
      for (; i < bytes; ++i) { hash ^= p[i]; hash *= 1099511628211ULL; }
    };
    const auto mix_value = [&](std::uint64_t value) { mix_bytes(&value, sizeof(value)); };
    const auto mix_lane = [&](const auto& lane) {
      mix_value(static_cast<std::uint64_t>(lane.size()));
      if (!lane.empty()) mix_bytes(lane.data(), lane.size() * sizeof(lane[0]));
    };
    mix_value(static_cast<std::uint64_t>(state.particles.size()));
    mix_lane(state.particles.mass_code);
    mix_lane(state.particle_sidecar.particle_id);
    mix_lane(state.particle_sidecar.species_tag);
    mix_lane(state.particle_sidecar.particle_flags);
    mix_lane(state.particle_sidecar.owning_rank);
    mix_lane(state.particle_sidecar.gravity_softening_comoving);
    mix_lane(state.particle_sidecar.has_gravity_softening_override);
    mix_value(static_cast<std::uint64_t>(state.cells.size()));
    mix_value(static_cast<std::uint64_t>(state.gas_cells.size()));
    mix_value(state.particleIndexGeneration());
    mix_value(static_cast<std::uint64_t>(state.particleMetadataRepresentation()));
    return hash;
  };
  // Captures the rollback point and opens the gravity cycle transaction.
  // Returns false on EVERY rank when admission failed on ANY rank (governor
  // refusal or allocation failure, both reduced collectively before any rank
  // advances): the caller then applies the configured memory policy. Any
  // other failure is collective and fatal.
  //
  // Admission and allocation order (no undercounted peak, no old/new overlap):
  // 1. Steady state (existing capacity covers the extents and a reservation
  //    exists): nothing is allocated; the capture writes into retained
  //    capacity that the reservation already covers.
  // 2. Growth: the old storage and its reservation are released FIRST, the
  //    governor is asked for the exact fresh-storage bytes, and only then are
  //    the lanes allocated; if the allocator retains more than requested, the
  //    reservation is grown with the new grant taken before the old one is
  //    dropped (transient over-count, never under-count).
  // 3. Any refusal (governor or std::bad_alloc) releases the storage on that
  //    rank; the collective vote then makes every rank fall back or fail.
  const auto take_cycle_snapshot = [&]() -> bool {
    std::exception_ptr failure;
    bool admitted_local = true;
    auto& store = rollback_store;
    try {
      if (store.armed || m_gravity.cycleTransactionOpen()) {
        throw std::logic_error("multirate PM cycle snapshot requested while a rollback point is armed");
      }
      if (!state.hasHomogeneousDmoMetadata()) {
        throw std::logic_error("multirate PM cycle rollback requires homogeneous DMO particle metadata");
      }
      const bool simulated_refusal = snapshot_refusals_remaining > 0U && snapshot_refusal_this_rank;
      const std::size_t rows = state.particles.size();
      const std::size_t sfc_rows = state.particle_sidecar.sfc_key.size();
      const std::size_t cell_rows = static_cast<std::size_t>(cells.elementCount());
      const bool governed = m_services.memory_governor != nullptr;
      if (simulated_refusal && !options.multirate_test_snapshot_allocation_failure) {
        admitted_local = false;
      } else if (simulated_refusal || !store.coversRows(rows, sfc_rows, cell_rows) ||
                 (governed && store.reserved_bytes == 0U)) {
        store.releaseStorage();
        const std::uint64_t bytes =
            MultirateRollbackStore::requiredBytes(rows, sfc_rows, cell_rows, state.metadata);
        if (governed) {
          try {
            store.reservation = m_services.memory_governor->reserve(
                core::MemoryClass::kPhaseResident, bytes, "time.multirate_cycle_rollback_snapshot");
            store.reservation.commit();
            store.reserved_bytes = bytes;
          } catch (const core::MemoryAdmissionError&) {
            admitted_local = false;
          }
        }
        if (admitted_local) {
          // Test seam: an allocation failure after a granted reservation.
          if (simulated_refusal) throw std::bad_alloc();
          store.reserveRows(rows, sfc_rows, cell_rows);
          ++store.growths;
        }
      }
      if (admitted_local) {
        const auto copy = [](const auto& lane, auto& out) { out.assign(lane.begin(), lane.end()); };
        copy(state.particles.position_x_comoving, store.position_x_comoving);
        copy(state.particles.position_y_comoving, store.position_y_comoving);
        copy(state.particles.position_z_comoving, store.position_z_comoving);
        copy(state.particles.velocity_x_peculiar, store.velocity_x_peculiar);
        copy(state.particles.velocity_y_peculiar, store.velocity_y_peculiar);
        copy(state.particles.velocity_z_peculiar, store.velocity_z_peculiar);
        copy(state.particle_sidecar.sfc_key, store.sfc_key);
        store.metadata = state.metadata;
        store.drift_time_code = state.homogeneousDmoLastDriftTimeCode();
        store.drift_scale_factor = state.homogeneousDmoLastDriftScaleFactor();
        store.frozen_fingerprint = frozen_state_fingerprint();
        store.integrator = integrator;
        particles.exportPersistentStateInto(store.particle_scheduler);
        cells.exportPersistentStateInto(store.cell_scheduler);
        store.pending = pending;
        m_gravity.captureTransactionalForceEndpoint(state, store.frozen_fingerprint, store.force_endpoint);
        // Reconcile the reservation with what is actually retained (allocator
        // rounding, metadata strings): grow, never shrink below capacity.
        const std::uint64_t retained = store.retainedCapacityBytes();
        if (governed && retained > store.reserved_bytes) {
          try {
            auto grown = m_services.memory_governor->reserve(
                core::MemoryClass::kPhaseResident, retained, "time.multirate_cycle_rollback_snapshot");
            grown.commit();
            store.reservation = std::move(grown);  // new grant held before the old is released
            store.reserved_bytes = retained;
          } catch (const core::MemoryAdmissionError&) {
            admitted_local = false;
          }
        }
        store.peak_reserved_bytes = std::max(store.peak_reserved_bytes, store.reserved_bytes);
        store.peak_retained_capacity_bytes = std::max(store.peak_retained_capacity_bytes, retained);
      }
    } catch (const std::bad_alloc&) {
      // Allocation failure of the snapshot lanes is an admission refusal,
      // not a fatal error: drop the partial copy and its capacity.
      admitted_local = false;
    } catch (...) { failure = std::current_exception(); }
    // A locally refused rank holds no rollback storage at all (reservation
    // and retained capacity go together, so the governor never under-counts).
    if (!admitted_local) store.releaseStorage();
    FailureCoordinator(m_services).rethrowCollectiveFailure(failure, "multirate PM cycle snapshot");
    // Collective admission: every rank opens the same cycle shape.
    const bool admitted = mpi.allreduceMinUint64(admitted_local ? 1U : 0U) == 1U;
    if (snapshot_refusals_remaining > 0U) --snapshot_refusals_remaining;
    if (!admitted) {
      store.armed = false;
      return false;
    }
    std::exception_ptr begin_failure;
    try {
      m_gravity.beginCycleTransaction();
      store.armed = true;
    } catch (...) { begin_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(begin_failure, "multirate PM cycle transaction begin");
    return true;
  };
  // Transactional rollback to the cycle opening. Every rank executes it (the
  // decision derives from globally reduced criteria). Order:
  // 1. Verify the frozen identity/layout certificate (fingerprint of the ID
  //    lane and every frozen lane) and the row count: row identity is proven,
  //    not assumed from a matching count.
  // 2. Restore the authoritative lanes and records (no allocation: equal
  //    extents), advance the gravity source generation strictly beyond every
  //    generation the rejected attempt used.
  // 3. GravityRuntime::rollbackCycleTransaction restores the opening endpoint
  //    ROW-ALIGNED (no ID hash map, no allocation; refusal is checked before
  //    any write), moves the attempt's diagnostics to execution cost and drops
  //    every source-history-dependent derived gravity state.
  // A failure at any step on any rank fails every rank (FailureCoordinator).
  const auto rollback_cycle = [&]() {
    auto& store = rollback_store;
    // Phase 1 (read-only, every rank): prove row identity and every restore
    // precondition, then agree collectively. A refusal on ANY rank fails the
    // run on EVERY rank before any rank has restored anything.
    std::exception_ptr validation_failure;
    std::uint64_t layout_certificate = 0U;
    // Prepare all potentially allocating scheduler reconstruction BEFORE any
    // authoritative lane is restored. The prepared objects are installed by
    // noexcept swap only after EVERY MPI rank has validated and prepared.
    std::optional<core::HierarchicalTimeBinScheduler> prepared_particles;
    std::optional<core::HierarchicalTimeBinScheduler> prepared_cells;
    core::StateMetadata prepared_metadata;
    // Preflight the transient reconstruction rather than allocating a second
    // O(N) scheduler without MemoryGovernor admission. Each generic imported
    // row needs at most 25 B of represented lanes; 48 B/row plus bin/list
    // overhead is deliberately conservative for allocator capacity rounding.
    // This is temporary admission for the rollback preparation, distinct from
    // the retained cycle-opening snapshot reservation.
    core::MemoryReservation prepared_scheduler_reservation;
    static_assert(std::is_nothrow_swappable_v<core::HierarchicalTimeBinScheduler>);
    static_assert(std::is_nothrow_swappable_v<core::StateMetadata>);
    try {
      if (options.multirate_test_rollback_failure_rank == m_services.mpi_context.worldRank()) {
        // Test seam: break the row-identity proof on this rank only, so the
        // production row-aligned validation refuses the restore here while
        // every peer would accept it.
        store.force_endpoint.layout_certificate ^= 1U;
      }
      if (!store.armed || !m_gravity.cycleTransactionOpen()) {
        throw std::logic_error("multirate PM cycle rollback requested without an armed rollback point");
      }
      layout_certificate = frozen_state_fingerprint();
      if (state.particles.size() != store.position_x_comoving.size() ||
          state.particle_sidecar.sfc_key.size() != store.sfc_key.size() ||
          layout_certificate != store.frozen_fingerprint) {
        throw std::logic_error(
            "multirate PM cycle rollback: state outside the rollback lanes changed inside the cycle");
      }
      m_gravity.validateTransactionalForceEndpoint(store.force_endpoint, state, layout_certificate);
      if (!std::isfinite(store.drift_time_code) || !std::isfinite(store.drift_scale_factor) ||
          store.drift_scale_factor <= 0.0) {
        throw std::logic_error("multirate PM cycle rollback: invalid saved drift epoch");
      }
      if (state.particles.time_bin.size() != state.particles.size()) {
        throw std::logic_error("multirate PM cycle rollback: time-bin mirror not materialized");
      }
      if (m_services.memory_governor != nullptr) {
        const auto rows = core::checkedSizeAdd(store.particle_scheduler.bin_index.size(),
            store.cell_scheduler.bin_index.size(), "multirate rollback scheduler prep rows");
        const auto bin_count = core::checkedSizeAdd(
            static_cast<std::size_t>(store.particle_scheduler.max_bin) + 1U,
            static_cast<std::size_t>(store.cell_scheduler.max_bin) + 1U,
            "multirate rollback scheduler prep bins");
        const auto bytes = core::checkedSizeAdd(
            core::checkedSizeMultiply(rows, 48U, "multirate rollback scheduler prep rows"),
            core::checkedSizeAdd(
                core::checkedSizeMultiply(bin_count, 256U, "multirate rollback scheduler prep bins"),
                4096U, "multirate rollback scheduler prep overhead"),
            "multirate rollback scheduler prep reservation");
        prepared_scheduler_reservation = m_services.memory_governor->reserve(
            core::MemoryClass::kPhaseResident, static_cast<std::uint64_t>(bytes),
            "time.multirate_rollback_scheduler_preparation");
        prepared_scheduler_reservation.commit();
      }
      // This import reconstructs per-bin vectors and may allocate or throw.
      // It must not run after positions/velocities have been restored.
      prepared_particles.emplace(store.particle_scheduler.max_bin);
      prepared_particles->importPersistentState(store.particle_scheduler);
      prepared_cells.emplace(store.cell_scheduler.max_bin);
      prepared_cells->importPersistentState(store.cell_scheduler);
      prepared_metadata = store.metadata;
      // Catch unexpectedly large retained capacities before the first
      // authoritative restore. An undersized model is a fail-closed error.
      if (prepared_scheduler_reservation.valid()) {
        const auto actual = core::checkedSizeAdd(
            static_cast<std::size_t>(prepared_particles->ownedCapacityBytes()),
            static_cast<std::size_t>(prepared_cells->ownedCapacityBytes()),
            "multirate rollback actual prepared scheduler bytes");
        const auto estimated = core::checkedSizeAdd(
            core::checkedSizeMultiply(
                core::checkedSizeAdd(store.particle_scheduler.bin_index.size(),
                    store.cell_scheduler.bin_index.size(), "rollback prepared count"),
                48U, "rollback prepared row bytes"),
            core::checkedSizeAdd(
                core::checkedSizeMultiply(
                    static_cast<std::size_t>(store.particle_scheduler.max_bin) +
                        static_cast<std::size_t>(store.cell_scheduler.max_bin) + 2U,
                    256U, "rollback prepared bin bytes"),
                4096U, "rollback prepared overhead"),
            "rollback prepared admitted bytes");
        if (actual > estimated) {
          throw std::logic_error("multirate PM cycle rollback: prepared schedulers exceed governed admission");
        }
      }
    } catch (...) { validation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(validation_failure, "multirate PM cycle rollback validation/preparation");
    // Phase 2 (all ranks prepared): install the preconstructed schedulers
    // without allocation. A refusal during scheduler import therefore cannot
    // leave particle lanes restored with partially rebuilt time-bin state.
    std::exception_ptr failure;
    try {
      const std::uint64_t generation_reached = state.gravitySourceGeneration();
      const auto restore = [](const auto& saved, auto& lane) { std::copy(saved.begin(), saved.end(), lane.begin()); };
      restore(store.position_x_comoving, state.particles.position_x_comoving);
      restore(store.position_y_comoving, state.particles.position_y_comoving);
      restore(store.position_z_comoving, state.particles.position_z_comoving);
      restore(store.velocity_x_peculiar, state.particles.velocity_x_peculiar);
      restore(store.velocity_y_peculiar, state.particles.velocity_y_peculiar);
      restore(store.velocity_z_peculiar, state.particles.velocity_z_peculiar);
      restore(store.sfc_key, state.particle_sidecar.sfc_key);
      using std::swap;
      swap(state.metadata, prepared_metadata);
      state.updateAllParticleDriftEpoch(store.drift_time_code, store.drift_scale_factor);
      integrator = store.integrator;
      swap(particles, *prepared_particles);
      swap(cells, *prepared_cells);
      pending = store.pending;
      while (state.gravitySourceGeneration() <= generation_reached) state.bumpGravitySourceGeneration();
      m_gravity.rollbackCycleTransaction(store.force_endpoint, state, layout_certificate);
      integrator.pm_source_generation = state.gravitySourceGeneration();
      syncTimeBinsFromSchedulers(particles, cells, state);
      report.multirate_rejected_blocks += pm_cycle.completed_blocks;
      store.armed = false;
    } catch (...) { failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(failure, "multirate PM cycle rollback");
    install_geometry();
    std::exception_ptr reconcile_failure;
    try { reconcile_runtime(); } catch (...) { reconcile_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(reconcile_failure, "multirate PM cycle rollback reconciliation");
  };
  while ((integrator.step_index < final_step || pm_cycle.open) &&
         integrator.current_time_code < config.numerics.t_code_end) {
    const bool pm_cycle_interior_start = multirate_pm && pm_cycle.open;
    if (pm_cycle_interior_start && integrator.step_index >= final_step) {
      // Soft stop: the step budget was reached inside an open cycle. Finish
      // the cycle (never shorten it for an execution budget).
      profiler.counters().addCount("time.hierarchical_multirate_soft_stop_extension_blocks");
    }
    const auto console_begin = RuntimeConsoleReporter::Clock::now();
    const double block_a_begin = integrator.current_scale_factor;
    std::exception_ptr preparation_failure;
    try {
      const std::size_t particle_rows = state.particles.size();
      const bool identity_hit = !options.hierarchical_test_disable_preparation_caches &&
          identity_rows.valid && identity_rows.rows == particle_rows &&
          identity_rows.particle_index_generation == state.particleIndexGeneration() &&
          workspace.gravity_particle_index_scratch.size() == particle_rows &&
          workspace.gravity_particle_index_scratch.data() == identity_rows.data;
      if (identity_hit) {
        ++report.hierarchical_identity_row_reuses;
      } else {
        identity_rows.valid = false;
        workspace.scratch.reset();
        workspace.prepareGravityParticleIndexScratch(particle_rows);
        for (std::size_t row = 0U; row < particle_rows; ++row) {
          workspace.gravity_particle_index_scratch[row] = core::checkedIntegralNarrow<std::uint32_t>(row, "hierarchical identity row");
        }
        identity_rows = IdentityRowCache{.valid = true, .rows = particle_rows,
            .particle_index_generation = state.particleIndexGeneration(),
            .data = workspace.gravity_particle_index_scratch.data()};
        ++report.hierarchical_identity_row_rebuilds;
      }
      core::RetainedCapacityTransaction mirror_plan(owned_runtime_bytes);
      mirror_plan.add(state.particles.time_bin, state.particles.size());
      mirror_plan.execute(m_services.memory_governor, core::MemoryClass::kPhaseResident, "time.hierarchical_bin_mirror");
      state.particles.time_bin.resize(state.particles.size(), 0U);
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical workspace admission");
    const auto all_rows = std::span<const std::uint32_t>(workspace.gravity_particle_index_scratch);
    // Endpoint force reuse. The previous block closed with an all-active
    // Tree+PM force at exactly these positions (no drift since); when every
    // rank still holds it complete (no migration, same committed PM field),
    // the block opens with it. Otherwise (first block, after a migration, or
    // a pre-v24 restart without the persisted endpoint) a synchronized
    // bootstrap evaluates it. Restart schema v24 persists the endpoint, so
    // uninterrupted and resumed runs take the same branch with the same
    // forces. The decision is collective.
    bool endpoint_reusable_local = false;
    preparation_failure = {};
    try {
      // Inside a multirate PM cycle the next block needs only the Tree
      // component at the current positions; the PM field stays the one built
      // at the cycle opening (its kick was already applied).
      endpoint_reusable_local = pm_cycle_interior_start
          ? m_gravity.hierarchicalShortRangeEndpointReusable(state)
          : m_gravity.hierarchicalEndpointForceReusable(state, integrator);
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical endpoint reuse check");
    const bool reuse_endpoint = mpi.allreduceMinDouble(endpoint_reusable_local ? 1.0 : 0.0) > 0.5;
    if (reuse_endpoint) {
      profiler.counters().addCount("time.hierarchical_endpoint_force_reuse");
      // The reused PM component was built at these source positions (the
      // reuse check proves it); restate that for a resumed run, whose
      // integrator does not persist this derived generation. Inside a
      // multirate cycle the PM field belongs to the cycle opening instead.
      if (!pm_cycle_interior_start) integrator.pm_source_generation = state.gravitySourceGeneration();
    } else if (pm_cycle_interior_start) {
      // Multirate cycle interior (e.g. after a migration): a synchronized
      // Tree-only force at the current positions. No PM refresh, no PM kick.
      profiler.counters().addCount("time.hierarchical_multirate_tree_bootstrap");
      core::StepContext sync{.state = state, .integrator_state = integrator,
          .active_set = {.particle_indices = all_rows, .cell_indices = {}},
          .active_gravity_particles = {}, .has_active_gravity_particles = false,
          .workspace = &workspace, .cosmology_background = background,
          .mode_policy = &mode_policy, .profiler_session = &profiler,
          .timeline_step = {.time_begin_code = integrator.current_time_code,
              .time_end_code = integrator.current_time_code,
              .scale_factor_begin = integrator.current_scale_factor,
              .scale_factor_end = integrator.current_scale_factor},
          .boundary = {.kind = core::StepBoundaryKind::kGlobalSynchronizationPoint,
              .restart_safe = false, .output_safe = false, .pm_refresh_allowed = false,
              .local_substep = false, .multirate_pm_cycle_interior = true},
          .pm_refresh_directive = {.force_refresh_surface = true, .cadence_opportunity_allowed = false,
              .sync_event_requested = false, .force_evaluation_scale_factor = integrator.current_scale_factor,
              .reason = core::PmRefreshDirective::Reason::kScheduledForceRefreshStage},
          .stage = core::IntegrationStage::kForceRefresh};
      sync.hierarchical_kdk = {.enabled = true, .include_long_range_force = false,
          .force_only_synchronization = true, .coarse_source_generation = state.gravitySourceGeneration(),
          .evaluation_tick = particles.currentTick()};
      dispatchStage(sync, false);
      if (sync.pm_refresh_directive.has_sync_event || !sync.pm_refresh_directive.solver_executed) {
        throw std::logic_error("multirate cycle-interior bootstrap must be a Tree-only synchronized solve");
      }
    } else {
      profiler.counters().addCount("time.hierarchical_endpoint_force_bootstrap");
      core::StepContext sync{.state = state, .integrator_state = integrator,
          .active_set = {.particle_indices = all_rows, .cell_indices = {}},
          .active_gravity_particles = {}, .has_active_gravity_particles = false,
          .workspace = &workspace, .cosmology_background = background,
          .mode_policy = &mode_policy, .profiler_session = &profiler,
          .timeline_step = {.time_begin_code = integrator.current_time_code,
              .time_end_code = integrator.current_time_code,
              .scale_factor_begin = integrator.current_scale_factor,
              .scale_factor_end = integrator.current_scale_factor},
          .boundary = {.kind = core::StepBoundaryKind::kGlobalSynchronizationPoint},
          .pm_refresh_directive = {.force_refresh_surface = true, .cadence_opportunity_allowed = true,
              .sync_event_requested = true, .force_evaluation_scale_factor = integrator.current_scale_factor,
              .reason = core::PmRefreshDirective::Reason::kScheduledForceRefreshStage},
          .stage = core::IntegrationStage::kForceRefresh};
      sync.hierarchical_kdk = {.enabled = true, .include_long_range_force = true,
          .force_only_synchronization = true, .coarse_source_generation = state.gravitySourceGeneration(),
          .evaluation_tick = particles.currentTick()};
      dispatchStage(sync, false);
      const auto& d = sync.pm_refresh_directive;
      if (!d.has_sync_event || !d.solver_executed || !d.refresh_long_range_field) {
        throw std::logic_error("hierarchical synchronized bootstrap did not commit a fresh split force");
      }
      const core::PmSyncEvent event{d.gravity_kick_opportunity, d.refresh_long_range_field,
          d.field_version, d.last_refresh_opportunity, d.field_built_step_index, d.field_built_scale_factor};
      integrator.pm_sync_state.commitKickOpportunity(event);
      integrator.pm_sync_state.commitRefresh(event);
      integrator.pm_long_range_field_valid = true;
      integrator.pm_source_generation = state.gravitySourceGeneration();
    }
    // Phase boundary: the synchronized bootstrap may grow gravity-owned
    // retained capacity (PM plans, tree and residual workspaces). Reconcile
    // the governor baseline before the rung-membership transaction measures
    // retained bytes; otherwise that transaction correctly reports a stale
    // baseline (repair register B-10). reconcile_runtime is collective.
    {
      std::exception_ptr reconcile_failure;
      try { reconcile_runtime(); } catch (...) { reconcile_failure = std::current_exception(); }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          reconcile_failure, "hierarchical post-bootstrap memory reconciliation");
    }
    const auto ax = m_gravity.particleAccelX();
    const auto ay = m_gravity.particleAccelY();
    const auto az = m_gravity.particleAccelZ();
    const auto nearest_pair_r2 = m_gravity.hierarchicalNearestPairR2();
    const auto softening = speciesSofteningByTag(config);
    const double eps = softening.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)];
    const double gravity_timestep_eta = config.numerics.gravity_timestep_eta;
    // Softened-encounter criterion. The acceleration criterion
    // eta sqrt(eps a^3/|A|) is blind to a softened pair passing through
    // coincidence (|A| -> 0 at r -> 0). The pair's tidal rate
    // T = G m/(r^2+eps^2)^{3/2} (comoving, scale-free; physical T/a^3) stays
    // finite and maximal there, so dt <= eta sqrt(a^3/T) bounds the step by
    // the softened two-body dynamical time eta sqrt(eps^3 a^3/(G m)) at
    // coincidence while being weaker than the acceleration criterion by
    // sqrt(r/eps) beyond r ~ eps. r is the nearest direct pair from the
    // synchronized Tree solve (restart-persisted with the v24 endpoint);
    // m_max, the global maximum particle mass, makes it conservative for
    // unequal masses. Nodes accepted as multipoles are not neighbours here.
    const bool encounter_criterion = config.numerics.hierarchical_encounter_timestep_enabled;
    const std::size_t mass_rows = state.particles.mass_code.size();
    const bool mass_max_hit = !options.hierarchical_test_disable_preparation_caches &&
        local_mass_max_cache.valid && local_mass_max_cache.rows == mass_rows &&
        local_mass_max_cache.particle_mass_generation == state.particleMassGeneration() &&
        local_mass_max_cache.data == state.particles.mass_code.data();
    if (!mass_max_hit) {
      // all_rows is the identity list over every row: scan the whole lane.
      double scanned = 0.0;
      for (const double mass : state.particles.mass_code) scanned = std::max(scanned, mass);
      local_mass_max_cache = LocalMassMaxCache{.valid = true, .rows = mass_rows,
          .particle_mass_generation = state.particleMassGeneration(),
          .data = state.particles.mass_code.data(), .value = scanned};
      ++report.hierarchical_mass_max_scans;
    } else {
      ++report.hierarchical_mass_max_reuses;
    }
    const double local_mass_max = local_mass_max_cache.value;
    const double mass_max = -mpi.allreduceMinDouble(-local_mass_max);
    const double encounter_g_mass = m_gravity.gravitationalConstantCode() * mass_max;
    const double a_cubed = integrator.current_scale_factor * integrator.current_scale_factor *
        integrator.current_scale_factor;
    // Timestep growth limiter: a particle's interval may at most double from
    // one synchronization point to the next (refinement is unrestricted).
    // The acceleration-only criterion sqrt(2 eta eps/|a|) is blind to the
    // instant a softened pair passes through coincidence, where |a| -> 0;
    // without this limit such a pair jumped to the coarsest rung for a whole
    // block and crossed its encounter in one step. The previous interval is
    // reconstructed from restart-persisted state only (the previous quantum,
    // integrator.dt_time_code, times the committed rung period), so resumed
    // and uninterrupted runs apply the same limit. No history before the
    // first completed step.
    const double previous_quantum = integrator.dt_time_code;
    const bool growth_limited = integrator.step_index > 0U &&
        std::isfinite(previous_quantum) && previous_quantum > 0.0;
    const auto criterion = [&](std::uint32_t row) {
      double dt = core::computeComovingGravityTimeStep({.softening_length_comoving_code = std::max(eps, 1.0e-12),
          .scale_free_acceleration_magnitude_code = std::sqrt(ax[row] * ax[row] + ay[row] * ay[row] + az[row] * az[row]),
          .scale_factor = integrator.current_scale_factor}, gravity_timestep_eta);
      if (encounter_criterion && std::isfinite(nearest_pair_r2[row]) && encounter_g_mass > 0.0) {
        const double d = nearest_pair_r2[row] + eps * eps;
        dt = std::min(dt, gravity_timestep_eta * std::sqrt(a_cubed * d * std::sqrt(d) / encounter_g_mass));
      }
      if (growth_limited) {
        dt = std::min(dt, 2.0 * previous_quantum *
            static_cast<double>(particles.binPeriodTicks(particles.binIndex(row))));
      }
      return dt;
    };
    double local_min = std::numeric_limits<double>::infinity();
    double local_speed_max = 0.0;
    double local_acceleration_max = 0.0;
    double local_mass_sum = 0.0;
    double local_mass_speed2_sum = 0.0;
    preparation_failure = {};
    try {
      if (ax.size() != state.particles.size() || ay.size() != ax.size() || az.size() != ax.size()) {
        throw std::logic_error("hierarchical timestep criteria lack synchronized total force rows");
      }
      if (encounter_criterion && nearest_pair_r2.size() != state.particles.size()) {
        throw std::logic_error("hierarchical encounter criterion lacks synchronized nearest-pair rows");
      }
      for (const auto row : all_rows) {
        local_min = std::min(local_min, criterion(row));
        const double speed = std::hypot(state.particles.velocity_x_peculiar[row],
            state.particles.velocity_y_peculiar[row], state.particles.velocity_z_peculiar[row]);
        const double acceleration = std::hypot(ax[row], ay[row], az[row]);
        if (!std::isfinite(speed) || !std::isfinite(acceleration)) {
          throw std::runtime_error("hierarchical coarse criterion received nonfinite speed or force");
        }
        local_speed_max = std::max(local_speed_max, speed);
        local_acceleration_max = std::max(local_acceleration_max, acceleration);
        const double mass = state.particles.mass_code[row];
        local_mass_sum += mass;
        local_mass_speed2_sum += mass * speed * speed;
      }
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical timestep criteria");
    const auto periods = particles.binPeriodTicks(particles.maxBin());
    double coarse_dt = mpi.allreduceMinDouble(local_min) * static_cast<double>(periods);
    const double speed_max = -mpi.allreduceMinDouble(-local_speed_max);
    const double acceleration_max = -mpi.allreduceMinDouble(-local_acceleration_max);
    const double dx = config.cosmology.box_size_x_mpc_comoving / config.numerics.treepm_pm_grid_nx;
    const double dy = config.cosmology.box_size_y_mpc_comoving / config.numerics.treepm_pm_grid_ny;
    const double dz = config.cosmology.box_size_z_mpc_comoving / config.numerics.treepm_pm_grid_nz;
    const double split_scale = config.numerics.treepm_asmth_cells*std::cbrt(dx*dy*dz);
    // Collective in both modes so every rank issues identical reductions.
    const double global_mass_sum = mpi.allreduceSumDouble(local_mass_sum);
    const double global_mass_speed2_sum = mpi.allreduceSumDouble(local_mass_speed2_sum);
    const double global_particle_count = mpi.allreduceSumDouble(static_cast<double>(all_rows.size()));
    const double displacement_scale_factor = background != nullptr ? integrator.current_scale_factor : 1.0;
    double displacement_dt = std::numeric_limits<double>::infinity();
    if (config.numerics.hierarchical_coarse_displacement_criterion ==
        core::HierarchicalCoarseDisplacementCriterion::kRmsSplitScale) {
      // The long-range field is Gaussian-smoothed on r_s, so a coarse PM
      // interval only needs the typical (mass-weighted RMS) particle to move
      // a fraction of min(r_s, mean interparticle spacing). Fast particles in
      // dense regions keep their own fine rungs for the short-range force.
      const double box_volume = config.cosmology.box_size_x_mpc_comoving *
          config.cosmology.box_size_y_mpc_comoving * config.cosmology.box_size_z_mpc_comoving;
      const double mean_spacing = global_particle_count > 0.0
          ? std::cbrt(box_volume / global_particle_count) : split_scale;
      const double speed_rms_comoving = global_mass_sum > 0.0
          ? std::sqrt(global_mass_speed2_sum / global_mass_sum) / displacement_scale_factor : 0.0;
      if (!std::isfinite(speed_rms_comoving) || !std::isfinite(mean_spacing) || mean_spacing <= 0.0) {
        throw std::runtime_error("hierarchical RMS displacement criterion received invalid inputs");
      }
      if (speed_rms_comoving > 0.0) {
        displacement_dt = config.numerics.hierarchical_coarse_rms_displacement_factor *
            std::min(split_scale, mean_spacing) / speed_rms_comoving;
      }
    } else {
      displacement_dt = core::computeComovingDisplacementTimeStep({
          .mesh_or_split_length_comoving_code = std::min({dx,dy,dz,split_scale}),
          .velocity_magnitude_peculiar_code = speed_max,
          .scale_free_acceleration_magnitude_code = acceleration_max,
          .scale_factor = displacement_scale_factor});
    }
    // Diagnostic only: name the criterion that bounds this coarse interval.
    const double particle_coarse_dt = coarse_dt;
    const char* coarse_limiter = "particle_gravity";
    const auto limit_coarse = [&](double candidate, const char* name) {
      if (candidate < coarse_dt) {
        coarse_dt = candidate;
        coarse_limiter = name;
      }
    };
    const double expansion_dt = background != nullptr
        ? core::computeCosmologyExpansionTimeStep(*background, integrator.current_scale_factor,
              config.numerics.cosmology_max_delta_ln_a, config.numerics.cosmology_max_hubble_time_fraction,
              integrator.time_si_per_code)
        : std::numeric_limits<double>::infinity();
    double boundary_time = config.numerics.t_code_end;
    if (pending.snapshot_interval_time_code > 0.0 && pending.next_snapshot_time_code > integrator.current_time_code) {
      boundary_time = std::min(boundary_time, pending.next_snapshot_time_code);
    }
    if (!multirate_pm) {
      // Reference hierarchical scheme: one coarse block = one PM interval.
      limit_coarse(displacement_dt, "pm_displacement");
      limit_coarse(expansion_dt, "cosmology_expansion");
      if (options.dt_time_code > 0.0) limit_coarse(options.dt_time_code, "explicit_dt");
      limit_coarse(boundary_time - integrator.current_time_code, "output_or_end_boundary");
    } else {
      // Inner block bound: the per-block reference criteria only. The PM
      // displacement criterion and output/end boundaries bound the outer
      // cycle instead (below).
      limit_coarse(expansion_dt, "cosmology_expansion");
      if (options.dt_time_code > 0.0) limit_coarse(options.dt_time_code, "explicit_dt");
      const double t_now = integrator.current_time_code;
      if (!pm_cycle.open) {
        double length = coarse_dt * static_cast<double>(multirate_blocks_per_cycle);
        const char* limiter = "target_blocks";
        if (displacement_dt < length) { length = displacement_dt; limiter = "pm_displacement"; }
        if (pm_cycle.retry_length_limit > 0.0 && pm_cycle.retry_length_limit < length) {
          length = pm_cycle.retry_length_limit;
          limiter = "rollback_retry";
        }
        bool boundary_clipped = false;
        if (boundary_time - t_now <= length) {
          length = boundary_time - t_now;
          limiter = "output_or_end_boundary";
          boundary_clipped = true;
        }
        if (!std::isfinite(length) || !(length > 0.0) || t_now + length <= t_now) {
          profiler.counters().addCount("time.hierarchical_multirate_recovery_failures");
          throw std::runtime_error("multirate PM cycle length is not finite positive representable");
        }
        // Rollback point, before rung assignment and the opening PM kick.
        // A cycle whose first block reaches its end has no interior boundary
        // and cannot be rejected: no snapshot is needed. When snapshot
        // admission fails on any rank (collectively agreed), the configured
        // memory policy decides: a one-block cycle on every rank (its single
        // block still honours the particle criteria and, since
        // coarse_dt < length <= Δt_PMdisp, the PM displacement bound), or a
        // collective failure (reproducibility mode).
        if (coarse_dt < length * (1.0 - 1.0e-12) && !take_cycle_snapshot()) {
          if (config.numerics.treepm_multirate_rollback_memory_policy ==
              core::MultirateRollbackMemoryPolicy::kFail) {
            profiler.counters().addCount("time.hierarchical_multirate_memory_policy_failures");
            throw std::runtime_error(
                "multirate PM cycle rollback snapshot was refused by the memory governor on at least one rank "
                "and numerics.treepm_multirate_rollback_memory_policy = fail");
          }
          length = coarse_dt;
          limiter = "rollback_memory";
          boundary_clipped = false;
          ++report.multirate_memory_limited_cycles;
          profiler.counters().addCount("time.hierarchical_multirate_cycle_memory_limited");
        }
        pm_cycle.open = true;
        ++pm_cycle.index;
        pm_cycle.completed_blocks = 0U;
        pm_cycle.opening_time_code = t_now;
        pm_cycle.end_time_code = boundary_clipped ? boundary_time : t_now + length;
        pm_cycle.length_time_code = pm_cycle.end_time_code - t_now;
        pm_cycle.length_limiter = limiter;
        // Reset per cycle, not per attempt: a retried opening keeps its
        // rejection history (retries, forced-rollback bookkeeping).
        if (pm_cycle.retries == 0U) pm_cycle.forced_rejections = 0U;
        profiler.counters().addCount(std::string("time.hierarchical_multirate_cycle_limited_by.") + limiter);
        if (boundary_clipped) profiler.counters().addCount("time.hierarchical_multirate_cycle_shortened_by_boundary");
        const double half_cycle_time_code = 0.5 * pm_cycle.length_time_code;
        pm_cycle_plan = {};
        pm_cycle_plan.enabled = true;
        if (background != nullptr) {
          const double a0 = integrator.current_scale_factor;
          pm_cycle_plan.cycle_mid_scale_factor = core::advanceScaleFactorByCosmicTime(
              *background, a0, half_cycle_time_code * integrator.time_si_per_code);
          pm_cycle_plan.opening_pm_kick_factor_code = core::computeComovingDriftFactor(
              *background, a0, pm_cycle_plan.cycle_mid_scale_factor, 64U) / (integrator.time_si_per_code * a0);
        } else {
          pm_cycle_plan.opening_pm_kick_factor_code = half_cycle_time_code;
          pm_cycle_plan.closing_half_interval_time_code = half_cycle_time_code;
        }
      } else {
        // Interior boundary: the outer interval fixed at the opening must still
        // satisfy the PM displacement criterion with the current velocities.
        const bool displacement_violated = displacement_dt < pm_cycle.length_time_code * (1.0 - 1.0e-12);
        const bool forced = !displacement_violated && forced_rollbacks_remaining > 0U &&
            pm_cycle.forced_rejections < forced_rejections_per_cycle;
        if (displacement_violated || forced) {
          ++pm_cycle.retries;
          if (pm_cycle.retries > k_multirate_max_retries_per_cycle) {
            profiler.counters().addCount("time.hierarchical_multirate_recovery_failures");
            throw std::runtime_error(
                "multirate PM cycle could not satisfy the PM displacement criterion after " +
                std::to_string(k_multirate_max_retries_per_cycle) + " retries");
          }
          // Genuine rejection: shrink T to the newly required bound (at least
          // halving, so a single-block cycle is reached in O(log M) retries).
          // Test seam: retry with the identical T (rollback completeness).
          const double retry_limit = forced
              ? pm_cycle.length_time_code
              : std::min(0.5 * pm_cycle.length_time_code, displacement_dt);
          const unsigned retries = pm_cycle.retries;
          const std::uint64_t forced_rejections = pm_cycle.forced_rejections + (forced ? 1U : 0U);
          if (forced) --forced_rollbacks_remaining;
          profiler.counters().addCount("time.hierarchical_multirate_cycle_rejected_retried");
          rollback_cycle();
          pm_cycle.staged_records.clear();
          pm_cycle.open = false;
          pm_cycle.retries = retries;
          pm_cycle.retry_length_limit = retry_limit;
          pm_cycle.forced_rejections = forced_rejections;
          ++report.multirate_pm_cycle_rollbacks;
          if (forced) ++report.multirate_forced_test_rollbacks;
          continue;
        }
      }
      // Block interval inside the open cycle: the reference per-block bound,
      // clipped so the blocks tile [opening, end] exactly; a sliver shorter
      // than a quarter block is avoided by splitting the remainder in two.
      const double remaining = pm_cycle.end_time_code - t_now;
      bool closes = false;
      if (coarse_dt >= remaining * (1.0 - 1.0e-12)) {
        coarse_dt = remaining;
        // Endpoint rounding: the driver ends the block at t + quantum*periods
        // = t + coarse_dt (periods is a power of two, so the split is exact),
        // but fl(t + fl(end - t)) can miss the cycle end by an ulp, leaving an
        // unrepresentable sliver before t_end. Nudge the interval by at most a
        // few ulps so the closing tick lands exactly on the planned end.
        for (int nudge = 0; nudge < 4 && t_now + coarse_dt != pm_cycle.end_time_code; ++nudge) {
          coarse_dt = std::nextafter(coarse_dt, t_now + coarse_dt < pm_cycle.end_time_code
              ? std::numeric_limits<double>::infinity() : 0.0);
        }
        if (t_now + coarse_dt != pm_cycle.end_time_code) {
          profiler.counters().addCount("time.hierarchical_multirate_cycle_end_inexact");
        }
        coarse_limiter = "multirate_cycle_end";
        closes = true;
      } else if (remaining - coarse_dt < 0.25 * coarse_dt) {
        coarse_dt = 0.5 * remaining;
        coarse_limiter = "multirate_cycle_end_split";
      }
      pm_cycle_plan.opens_cycle = pm_cycle.completed_blocks == 0U;
      pm_cycle_plan.closes_cycle = closes;
    }
    profiler.counters().addCount(std::string("time.hierarchical_coarse_limited_by.") + coarse_limiter);
    const double quantum = coarse_dt / static_cast<double>(periods);
    if (!std::isfinite(quantum) || quantum <= 0.0 || integrator.current_time_code + quantum <= integrator.current_time_code) {
      // Fail closed before any state is committed: no representable rung can
      // honour the required interval. Identical on every rank (global inputs).
      profiler.counters().addCount("time.hierarchical_quantum_unrepresentable");
      throw std::runtime_error("hierarchical quantum is not finite positive representable");
    }
    integrator.dt_time_code = quantum;
    ReferenceWorkflowReport::HierarchicalBlockRecord block_record{
        .step_index = integrator.step_index,
        .time_begin_code = integrator.current_time_code,
        .time_end_code = integrator.current_time_code,
        .coarse_interval_time_code = coarse_dt,
        .quantum_time_code = quantum,
        .global_min_particle_criterion_time_code = particle_coarse_dt / static_cast<double>(periods),
        .coarse_limiter = coarse_limiter,
        .multirate = multirate_pm,
        .pm_cycle_index = multirate_pm ? pm_cycle.index : 0U,
        .pm_cycle_block_index = multirate_pm ? pm_cycle.completed_blocks : 0U,
        .opens_pm_cycle = !multirate_pm || pm_cycle_plan.opens_cycle,
        .closes_pm_cycle = !multirate_pm || pm_cycle_plan.closes_cycle,
        .pm_cycle_opening_time_code = multirate_pm ? pm_cycle.opening_time_code : integrator.current_time_code,
        .pm_cycle_end_time_code = multirate_pm ? pm_cycle.end_time_code : integrator.current_time_code + coarse_dt,
        .pm_cycle_retries = multirate_pm ? pm_cycle.retries : 0U,
        .pm_displacement_limit_time_code = displacement_dt,
        .pm_cycle_length_limiter = multirate_pm ? std::string(pm_cycle.length_limiter) : std::string("coarse_block"),
    };
    const core::TimeStepLimits limits{.min_dt_time_code = quantum,
        .max_dt_time_code = coarse_dt, .max_bin = particles.maxBin()};
    preparation_failure = {};
    try {
      for (const auto row : all_rows) particles.submitCandidateTimeStep(row, criterion(row), limits,
          core::TimeStepCandidateSource::kGravityAcceleration, "hierarchical_synchronized_gravity");
      core::RetainedCapacityTransaction bins_plan(owned_runtime_bytes);
      particles.planSynchronizedCandidateCapacity(bins_plan);
      bins_plan.execute(m_services.memory_governor, core::MemoryClass::kPhaseResident, "time.hierarchical_bin_membership");
      particles.commitSynchronizedCandidates();
      syncTimeBinsFromSchedulers(particles, cells, state);
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical synchronized rung assignment");
    preparation_failure = {};
    try {
      internal::latchOutputRequestForCompletedStep(config, options, integrator.step_index + 1U,
        integrator.current_time_code + coarse_dt, pending);
      profiler.recordEvent(core::RuntimeEvent{.event_kind = "time.hierarchical_block", .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "core.time", .step_index = integrator.step_index,
        .simulation_time_code = integrator.current_time_code, .scale_factor = integrator.current_scale_factor,
        .message = "source-implemented synchronized power-of-two DMO KDK block; qualification pending",
        .payload = {{"quantum_time_code", formatRuntimeDouble(quantum)},
            {"coarse_interval_time_code", formatRuntimeDouble(coarse_dt)},
            {"pm_displacement_limit_time_code", formatRuntimeDouble(displacement_dt)}, {"fine_ticks", std::to_string(periods)},
            {"particle_coarse_limit_time_code", formatRuntimeDouble(particle_coarse_dt)},
            {"expansion_limit_time_code", formatRuntimeDouble(expansion_dt)},
            {"coarse_limiter", coarse_limiter},
            {"coarse_displacement_criterion",
             config.numerics.hierarchical_coarse_displacement_criterion ==
                     core::HierarchicalCoarseDisplacementCriterion::kRmsSplitScale
                 ? "rms_split_scale" : "max_speed_mesh"},
            {"rank_local_bin_occupancy", formatBinOccupancy(particles, all_rows)},
            {"pm_policy", multirate_pm ? "multirate_outer_kick" : "coarse_endpoint_half_kicks"},
            {"pm_cycle_index", std::to_string(multirate_pm ? pm_cycle.index : 0U)},
            {"pm_cycle_length_time_code", formatRuntimeDouble(multirate_pm ? pm_cycle.length_time_code : coarse_dt)},
            {"pm_cycle_length_limiter", multirate_pm ? pm_cycle.length_limiter : "coarse_block"},
            {"pm_cycle_block_index", std::to_string(multirate_pm ? pm_cycle.completed_blocks : 0U)},
            {"rung_assignment", "fixed_within_block"}}});
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical output and event preparation");
    core::MemoryReservation timeline_reservation;
    preparation_failure = {};
    try {
      if (m_services.memory_governor != nullptr) {
        timeline_reservation = m_services.memory_governor->reserve(core::MemoryClass::kPhaseResident,
            core::checkedSizeMultiply(core::k_hierarchical_timeline_capacity, sizeof(double), "hierarchical timeline bytes"),
            "time.hierarchical_timeline");
        timeline_reservation.commit();
      }
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical timeline admission");
    parallel::DecompositionRuntimeMeasurements block_work{};
    m_lifecycle.executeHierarchicalBlockWithDispatcher(state, integrator, particles, cells,
        [&](core::StepContext& context, bool safe_output) {
          if (context.stage == core::IntegrationStage::kAnalysisHooks && context.hierarchical_kdk.synchronization_end) {
            preparation_failure = {};
            try {
              if (background != nullptr &&
                  (std::log(integrator.current_scale_factor / block_a_begin) >
                      config.numerics.cosmology_max_delta_ln_a * (1.0 + 1.0e-10) ||
                   ((config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kScaleFactor ||
                     config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kLogScaleFactor) &&
                    integrator.current_scale_factor > config.numerics.a_end + 1.0e-12 * std::max(1.0, config.numerics.a_end)))) {
                throw std::runtime_error("hierarchical block crossed the configured cosmological interval or endpoint");
              }
            } catch (...) { preparation_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical cosmological endpoint");
            std::exception_ptr boundary_failure;
            try { reconcile_runtime(); } catch (...) { boundary_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(boundary_failure, "hierarchical pre-migration memory reconciliation");
            // Multirate: ownership is frozen inside an open PM cycle so the
            // opening snapshot stays a complete rollback point (no decomposition
            // epoch, geometry or row identity changes to undo). Rebalancing
            // happens at the cycle end, a closed synchronization point.
            const bool migrated = !context.boundary.multirate_pm_cycle_interior &&
                m_migration_balance.rebalance(state, particles, cells,
                block_work, context.active_set.particle_indices,
                expected_ids, integrator.step_index);
            if (migrated) {
              m_gravity.commitParticleDecompositionChange();
              install_geometry();
              // Ownership commit invalidates the former target-row span. Analysis
              // and output consume coherent state, never the departed active view.
              context.active_set = {};
            }
            boundary_failure = {};
            try {
              syncTimeBinsFromSchedulers(particles, cells, state);
              reconcile_runtime();
            } catch (...) { boundary_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(boundary_failure, "hierarchical post-migration mirrors");
          }
          dispatchStage(context, safe_output);
          if (context.stage == core::IntegrationStage::kForceRefresh) {
            std::exception_ptr work_failure;
            try {
              const auto current = m_gravity.lastRuntimeDecompositionMeasurements();
              const auto add = [](std::uint64_t a, std::uint64_t b) {
                return core::checkedMemoryBytesAdd(a, b, "hierarchical block work counter");
              };
              block_work.tree_pair_evaluations_recent = add(block_work.tree_pair_evaluations_recent, current.tree_pair_evaluations_recent);
              block_work.incoming_tree_pair_evaluations_recent = add(block_work.incoming_tree_pair_evaluations_recent, current.incoming_tree_pair_evaluations_recent);
              block_work.tree_remote_request_bytes_recent = add(block_work.tree_remote_request_bytes_recent, current.tree_remote_request_bytes_recent);
              block_work.pm_mesh_cells_touched_recent = add(block_work.pm_mesh_cells_touched_recent, current.pm_mesh_cells_touched_recent);
              block_work.pm_fft_transpose_bytes_recent = add(block_work.pm_fft_transpose_bytes_recent, current.pm_fft_transpose_bytes_recent);
              block_work.ghost_exchange_bytes_recent = add(block_work.ghost_exchange_bytes_recent, current.ghost_exchange_bytes_recent);
              block_work.tree_wall_ms_recent += current.tree_wall_ms_recent;
              block_work.pm_wall_ms_recent += current.pm_wall_ms_recent;
              block_work.has_measurements = true;
              block_work.has_spatial_tree_work = current.has_spatial_tree_work;
              block_work.spatial_tree_work_per_target = current.spatial_tree_work_per_target;
            } catch (...) { work_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(work_failure, "hierarchical block work aggregation");
          }
        }, background, workspace, &mode_policy, &profiler,
        [&](std::exception_ptr failure, std::string_view phase) {
          FailureCoordinator(m_services).rethrowCollectiveFailure(failure, phase);
        }, pm_cycle_plan);
    if (multirate_pm) {
      // Advance or close the outer PM cycle after the block completed.
      if (pm_cycle_plan.closes_cycle) {
        pm_cycle.open = false;
        if (m_gravity.cycleTransactionOpen()) m_gravity.commitCycleTransaction();
        rollback_store.armed = false;
        pm_cycle.retries = 0U;
        pm_cycle.retry_length_limit = 0.0;
        ++report.multirate_pm_cycles_completed;
        profiler.counters().addCount("time.hierarchical_multirate_cycles_completed");
      } else {
        ++pm_cycle.completed_blocks;
      }
    }
    {
      // Bounded per-block evidence for regression tests (recent 512 blocks).
      constexpr std::size_t k_block_record_limit = 512U;
      const auto publish = [&](const ReferenceWorkflowReport::HierarchicalBlockRecord& record) {
        if (report.hierarchical_block_records.size() == k_block_record_limit) {
          report.hierarchical_block_records.erase(report.hierarchical_block_records.begin());
        }
        report.hierarchical_block_records.push_back(record);
      };
      block_record.time_end_code = integrator.current_time_code;
      if (!multirate_pm) {
        publish(block_record);
      } else {
        pm_cycle.staged_records.push_back(block_record);
        if (pm_cycle_plan.closes_cycle) {
          for (const auto& record : pm_cycle.staged_records) publish(record);
          pm_cycle.staged_records.clear();
        }
      }
    }
    if (timeline_reservation.valid()) timeline_reservation.release();
    preparation_failure = {};
    try { reconcile_runtime(); } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical block memory report");
    profiler.counters().addCount("time.hierarchical_coarse_blocks");
    profiler.counters().addCount("time.hierarchical_fine_drifts", periods);
    const auto global_particles = mpi.allreduceSumUint64(static_cast<std::uint64_t>(state.particles.size()));
    if (m_services.console_reporter != nullptr) {
      RuntimeConsoleStepStatus status;
      status.step_index = integrator.step_index;
      status.t_code = integrator.current_time_code;
      status.dt_time_code = coarse_dt;
      status.a_scale = integrator.current_scale_factor;
      status.redshift = integrator.current_redshift;
      status.active_particle_count = global_particles;
      status.total_particle_count = global_particles;
      status.wall_step_seconds = std::chrono::duration<double>(RuntimeConsoleReporter::Clock::now() - console_begin).count();
      status.pm_activity = "coarse_endpoints";
      status.output_activity = "coarse_sync";
      m_services.console_reporter->emitStep(status);
    }
  }
  report.completed_steps = integrator.step_index - run_start_step;
  report.requested_steps = requested_steps;
  // Accepted-trajectory value: blocks beyond the budget that closed an open
  // cycle (blocks of rejected attempts were rolled back with step_index; the
  // profiler counter time.hierarchical_multirate_soft_stop_extension_blocks
  // counts executed work, rejected attempts included).
  report.multirate_soft_stop_extension_blocks = multirate_pm && report.completed_steps > requested_steps
      ? report.completed_steps - requested_steps : 0U;
  if (multirate_pm && (pm_cycle.open || rollback_store.armed || m_gravity.cycleTransactionOpen())) {
    throw std::logic_error("hierarchical segment ended with an open multirate PM cycle");
  }
  // Rollback storage: report the governed high-water and release it (the
  // RAII guard also releases it on failure paths).
  report.multirate_rollback_snapshot_peak_bytes = rollback_store.peak_reserved_bytes;
  report.multirate_rollback_snapshot_peak_capacity_bytes = rollback_store.peak_retained_capacity_bytes;
  report.multirate_rollback_snapshot_allocations = rollback_store.growths;
  rollback_store.releaseStorage();
  report.multirate_rollback_snapshot_retained_bytes = rollback_store.reserved_bytes;
  report.final_time_code = integrator.current_time_code;
  report.final_scale_factor = integrator.current_scale_factor;
}

void TimeCoordinator::executeSingleStep(
    core::SimulationState& state,
    core::IntegratorState& integrator_state,
    core::ActiveSetDescriptor active_set,
    const core::LambdaCdmBackground* cosmology_background,
    core::TransientStepWorkspace* workspace,
    const core::ModePolicy* mode_policy,
    core::ProfilerSession* profiler_session,
    std::optional<std::uint64_t> expected_scheduler_tick,
    core::StepBoundaryKind requested_boundary_kind) {
  m_lifecycle.executeSingleStepWithDispatcher(
      state,
      integrator_state,
      active_set,
      [this](core::StepContext& context, bool require_output_safe_boundary) {
        dispatchStage(context, require_output_safe_boundary);
      },
      cosmology_background,
      workspace,
      mode_policy,
      profiler_session,
      expected_scheduler_tick,
      requested_boundary_kind);
}

void TimeCoordinator::executeOutputBoundary(
    core::SimulationState& state,
    core::IntegratorState& integrator_state,
    core::ProfilerSession* profiler_session,
    core::StepBoundaryKind requested_boundary_kind) {
  m_lifecycle.executeOutputBoundaryWithDispatcher(
      state,
      integrator_state,
      [this](core::StepContext& context, bool require_output_safe_boundary) {
        dispatchStage(context, require_output_safe_boundary);
      },
      profiler_session,
      requested_boundary_kind);
}

double TimeCoordinator::updateAdaptiveTimeBins(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices,
    bool update_all_elements) {
  if (config.physics.enable_star_formation &&
      config.physics.star_formation_model ==
          core::StarFormationModelKind::kEffectiveMultiphaseTngLike &&
      m_effective_eos_table == nullptr) {
    const core::UnitSystem runtime_units = core::makeUnitSystem(
        config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);
    m_effective_eos_table = std::make_shared<physics::EffectiveMultiphaseEosTable>(
        physics::makeEffectiveMultiphaseEosConfig(config.physics),
        runtime_units,
        physics::makeEffectiveIsmReferenceCoolingProvider(config.physics));
  }
  // G1.2: recursive subcycling relaxes the hydro CFL candidate by level.
  // The GLOBAL minimum populated level is agreed here (every caller is a
  // rank-uniform collective point of the rung-zero loop) so the relaxation,
  // and hence the accepted interval, does not depend on the decomposition.
  int recursive_amr_global_min_level = -1;
  if (config.numerics.amr_time_stepping == core::AmrTimeSteppingConfigMode::kRecursiveSubcycling) {
    std::uint64_t local_min_level = std::numeric_limits<std::uint64_t>::max();
    for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
      if (state.patches.cell_count[patch] == 0U || state.patches.level[patch] < 0) continue;
      local_min_level = std::min<std::uint64_t>(local_min_level, static_cast<std::uint64_t>(state.patches.level[patch]));
    }
    const std::uint64_t global_min_level = m_services.mpi_context.allreduceMinUint64(local_min_level);
    if (global_min_level != std::numeric_limits<std::uint64_t>::max()) {
      recursive_amr_global_min_level = static_cast<int>(global_min_level);
    }
  }
  m_last_timestep_limiter = TimestepLimiterRecord{};
  m_last_timestep_limiter.dt_time_code = std::numeric_limits<double>::infinity();
  return updateAdaptiveTimeBinFamilies(
      state,
      particle_scheduler,
      gas_cell_scheduler,
      integrator_state,
      config,
      m_effective_eos_table.get(),
      particle_accel_x,
      particle_accel_y,
      particle_accel_z,
      cell_accel_x,
      cell_accel_y,
      cell_accel_z,
      active_particle_indices,
      active_cell_indices,
      update_all_elements,
      recursive_amr_global_min_level,
      &m_last_timestep_limiter);
}

void TimeCoordinator::coordinatePmRefreshDirective(core::StepContext& context) {
  if (context.stage != core::IntegrationStage::kGravityKickPre &&
      context.stage != core::IntegrationStage::kForceRefresh &&
      context.stage != core::IntegrationStage::kGravityKickPost) {
    return;
  }

  if (context.hierarchical_kdk.enabled) {
    auto& d = context.pm_refresh_directive;
    const auto& mpi = m_services.mpi_context;
    const auto votes = mpi.allreduceSumUint64(d.sync_event_requested ? 1ULL : 0ULL);
    if (votes != 0U && votes != static_cast<std::uint64_t>(mpi.worldSize())) {
      throw std::logic_error("hierarchical PM synchronization directive diverged across ranks");
    }
    if (d.sync_event_requested) {
      std::exception_ptr directive_failure;
      try {
        if (!d.cadence_opportunity_allowed || !context.boundary.pm_refresh_allowed ||
            !context.hierarchical_kdk.include_long_range_force) {
          throw std::logic_error("hierarchical PM refresh requested outside a synchronized split endpoint");
        }
      } catch (...) { directive_failure = std::current_exception(); }
      FailureCoordinator(m_services).rethrowCollectiveFailure(directive_failure, "hierarchical PM endpoint eligibility");
      const bool all_valid = mpi.allreduceSumUint64(context.integrator_state.pm_long_range_field_valid ? 1ULL : 0ULL) ==
          static_cast<std::uint64_t>(mpi.worldSize());
      directive_failure = {};
      try {
        core::materializePmRefreshDirective(context, all_valid);
        if (!d.refresh_long_range_field) throw std::logic_error("hierarchical PM endpoint cadence must refresh");
      } catch (...) { directive_failure = std::current_exception(); }
      FailureCoordinator(m_services).rethrowCollectiveFailure(directive_failure, "hierarchical PM endpoint materialization");
    }
    return;
  }

  auto& directive = context.pm_refresh_directive;
  const auto& mpi_context = m_services.mpi_context;
  const std::uint64_t world_size = static_cast<std::uint64_t>(
      std::max(mpi_context.worldSize(), 1));

  const auto materialize_single_rank = [&]() {
    if (!directive.sync_event_requested) {
      return;
    }
    if (!directive.cadence_opportunity_allowed) {
      if (directive.reason == core::PmRefreshDirective::Reason::kSourceMutationForceRefresh) {
        throw std::runtime_error(
            "gravity source state changed after the scheduled PM refresh on a boundary where long-range repair is illegal");
      }
      throw std::runtime_error(
          "PM cadence event requested outside an integrator-authorized refresh opportunity");
    }
    core::materializePmRefreshDirective(
        context, context.integrator_state.pm_long_range_field_valid);
  };

  if (!mpi_context.isEnabled() || world_size <= 1U) {
    materialize_single_rank();
    return;
  }

  // All ranks reach this fixed-size reduction before any branch-dependent
  // TreePM work. Local predicates are evidence only; the values below define
  // the one distributed refresh decision for this stage.
  std::array<std::uint64_t, 5> votes{
      context.integrator_state.pm_refresh_enabled ? 1ULL : 0ULL,
      context.boundary.pm_refresh_allowed ? 1ULL : 0ULL,
      context.boundary.local_substep ? 1ULL : 0ULL,
      context.integrator_state.pm_long_range_field_valid ? 1ULL : 0ULL,
      context.integrator_state.pm_source_generation != context.state.gravitySourceGeneration() ? 1ULL : 0ULL,
  };
  mpi_context.allreduceSumUint64sInPlace(votes);

  const bool pm_enablement_consistent = votes[0] == 0U || votes[0] == world_size;
  const bool global_pm_enabled = votes[0] == world_size;
  const bool global_refresh_allowed = votes[1] == world_size;
  const bool any_local_substep = votes[2] != 0U;
  const bool global_field_valid = votes[3] == world_size;
  const bool field_validity_mixed = votes[3] != 0U && votes[3] != world_size;
  const bool global_source_changed = votes[4] != 0U;

  std::exception_ptr coordination_failure;
  try {
    if (!pm_enablement_consistent) {
      throw std::runtime_error(
          "TreePM PM-refresh enablement diverged across ranks before cadence coordination");
    }
    if (!global_pm_enabled && !global_field_valid) {
      throw std::runtime_error(
          "TreePM long-range PM field is invalid while PM refresh is disabled");
    }

    // Rebuild the directive from global semantics rather than preserving a
    // rank-local candidate. The event itself remains uncommitted until the
    // gravity callback succeeds and the core integrator accepts it.
    directive.has_sync_event = false;
    directive.refresh_long_range_field = false;
    directive.solver_executed = false;
    directive.sync_stage = core::PmSyncStage::kNone;
    directive.gravity_kick_opportunity = 0U;
    directive.field_version = 0U;
    directive.last_refresh_opportunity = 0U;
    directive.field_built_step_index = 0U;
    directive.field_built_scale_factor = 1.0;

    switch (context.stage) {
      case core::IntegrationStage::kGravityKickPre: {
        const bool global_initial_bootstrap_needed =
            global_pm_enabled && !global_field_valid;
        const bool global_initial_bootstrap_allowed =
            !any_local_substep ||
            (context.integrator_state.step_index == 0U &&
             global_initial_bootstrap_needed);
        directive.initial_cache_bootstrap_allowed =
            global_initial_bootstrap_allowed;
        directive.sync_event_requested = global_initial_bootstrap_needed;
        directive.cadence_opportunity_allowed =
            global_initial_bootstrap_needed && global_initial_bootstrap_allowed;
        if (global_initial_bootstrap_needed) {
          if (!global_initial_bootstrap_allowed) {
            throw std::runtime_error(
                "TreePM distributed initial PM bootstrap is not legal at this integration boundary");
          }
          directive.reason =
              core::PmRefreshDirective::Reason::kInitialForceBootstrap;
          directive.force_evaluation_scale_factor =
              context.timeline_step.scale_factor_begin;
        } else {
          directive.reason = core::PmRefreshDirective::Reason::kNone;
        }
        break;
      }
      case core::IntegrationStage::kForceRefresh: {
        directive.force_refresh_surface = true;
        directive.requires_predicted_inactive_sources = any_local_substep;
        directive.reason =
            core::PmRefreshDirective::Reason::kScheduledForceRefreshStage;
        directive.force_evaluation_scale_factor =
            context.timeline_step.scale_factor_end;
        directive.cadence_opportunity_allowed =
            global_pm_enabled && global_refresh_allowed;
        directive.sync_event_requested = directive.cadence_opportunity_allowed;
        if (global_pm_enabled && !global_field_valid &&
            !global_refresh_allowed) {
          throw std::runtime_error(
              "TreePM distributed PM field is invalid but no rank-global refresh boundary is available");
        }
        break;
      }
      case core::IntegrationStage::kGravityKickPost: {
        directive.force_refresh_surface = global_source_changed;
        directive.cadence_opportunity_allowed = false;
        directive.sync_event_requested = false;
        directive.reason = core::PmRefreshDirective::Reason::kNone;
        if (global_source_changed) {
          if (!global_pm_enabled) {
            throw std::runtime_error(
                "TreePM gravity sources changed while PM refresh is disabled");
          }
          if (!global_refresh_allowed) {
            throw std::runtime_error(
                "gravity source state changed after the scheduled PM refresh on a boundary where long-range repair is illegal");
          }
          directive.reason =
              core::PmRefreshDirective::Reason::kSourceMutationForceRefresh;
          directive.force_evaluation_scale_factor =
              context.timeline_step.scale_factor_end;
          directive.cadence_opportunity_allowed = true;
          directive.sync_event_requested = true;
        } else if (global_pm_enabled && !global_field_valid) {
          if (field_validity_mixed) {
            throw std::runtime_error(
                "TreePM long-range PM validity diverged across ranks without an authorized refresh request");
          }
          throw std::runtime_error(
              "TreePM post-kick reached an invalid long-range PM field without a source-mutation refresh request");
        }
        break;
      }
      default:
        break;
    }

    if (directive.sync_event_requested) {
      core::materializePmRefreshDirective(context, global_field_valid);
    }
  } catch (...) {
    coordination_failure = std::current_exception();
  }

  FailureCoordinator(m_services).rethrowCollectiveFailure(
      coordination_failure, "TreePM PM-refresh rank-global directive coordination");
}

void TimeCoordinator::dispatchStage(
    core::StepContext& context,
    bool require_output_safe_boundary) {
  (void)require_output_safe_boundary;
  context.particle_scheduler = &m_time_state.m_particle_scheduler;
  context.gas_cell_scheduler = &m_time_state.m_gas_cell_scheduler;
  context.newly_created_particle_ids = &m_newly_created_particle_ids;
  const std::string stage_name =
      "stage." + std::string(core::integrationStageName(context.stage));
  COSMOSIM_PROFILE_SCOPE(context.profiler_session, stage_name);
  if (context.profiler_session != nullptr) {
    context.profiler_session->counters().addCount(
        stage_name + ".invocations", 1U);
  }

  coordinatePmRefreshDirective(context);

  SimulationRuntimeEpochSource epoch_source(
      context.state, m_time_state.m_particle_scheduler, context.integrator_state);
  internal::RuntimeStageResourceBundle stage_resources(context);
  AnalysisStageView audit_view(
      RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
  if (!context.hierarchical_kdk.force_only_synchronization) {
    m_execution_plan.executeAuditStage(context.stage, audit_view);
  }
  if (context.hierarchical_kdk.enabled &&
      (context.stage == core::IntegrationStage::kHydroUpdate ||
       context.stage == core::IntegrationStage::kSourceTerms ||
       ((context.stage == core::IntegrationStage::kAnalysisHooks || context.stage == core::IntegrationStage::kOutputCheck) &&
        !context.hierarchical_kdk.synchronization_end))) return;
  // Interior multirate PM cycle boundaries carry an open outer PM half-kick:
  // analysis would observe intermediate momenta and a rollback would make it
  // non-reproducible. Output is held by OutputRestartRuntime; analysis waits
  // for the cycle end.
  if (context.hierarchical_kdk.enabled && context.boundary.multirate_pm_cycle_interior &&
      context.stage == core::IntegrationStage::kAnalysisHooks) return;

  switch (context.stage) {
    case core::IntegrationStage::kGravityKickPre:
    case core::IntegrationStage::kForceRefresh:
    case core::IntegrationStage::kGravityKickPost: {
      GravityStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kDrift: {
      DriftParticleStageView view(
          RuntimeResourceLease(epoch_source, k_particle_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kHydroUpdate: {
      HydroAmrStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kSourceTerms: {
      SourceMutationStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kAnalysisHooks: {
      AnalysisStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kOutputCheck: {
      OutputRestartStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
  }
}

}  // namespace cosmosim::workflows
