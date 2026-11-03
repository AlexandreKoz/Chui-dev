#include "cosmosim/workflows/migration_balance_runtime.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "workflows/internal/amr_migration_payload.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/migration_wire.hpp"
#include "workflows/internal/runtime_decomposition_source_storage.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows::internal {
namespace {

void syncTimeBinsFromScheduler(
    const core::HierarchicalTimeBinScheduler& scheduler,
    core::SimulationState& state) {
  core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, core::TimeBinMirrorDomain::kParticles);
}

struct GasCellMigrationRecord {
  std::uint64_t particle_id = 0;
  std::uint64_t patch_id = 0;
  std::uint32_t patch_local_cell_offset = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double mass_code = 0.0;
  std::uint8_t time_bin = 0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double density_code = 0.0;
  double pressure_code = 0.0;
  double internal_energy_code = 0.0;
  double temperature_code = 0.0;
  double sound_speed_code = 0.0;
};

struct LocalPatchMetadata {
  std::uint64_t patch_id = 0;
  std::int32_t level = 0;
  std::uint32_t original_cell_count = 0;
  std::uint64_t parent_patch_id = 0;
  std::uint64_t morton_key = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double extent_x_comoving = 0.0;
  double extent_y_comoving = 0.0;
  double extent_z_comoving = 0.0;
  std::uint16_t cell_dim_x = 0;
  std::uint16_t cell_dim_y = 0;
  std::uint16_t cell_dim_z = 0;
  std::uint32_t owning_rank = 0;

  [[nodiscard]] bool hasExplicitCartesianGeometry() const noexcept {
    return cell_dim_x != 0U || cell_dim_y != 0U || cell_dim_z != 0U ||
        extent_x_comoving != 0.0 || extent_y_comoving != 0.0 || extent_z_comoving != 0.0;
  }
};

struct LocalGasCellMigrationState {
  std::unordered_map<std::uint64_t, GasCellMigrationRecord> gas_records_by_particle_id;
  std::unordered_map<std::uint64_t, LocalPatchMetadata> patch_metadata_by_patch_id;
};

static_assert(std::is_trivially_copyable_v<GasCellMigrationRecord>);
static_assert(std::is_trivially_copyable_v<LocalPatchMetadata>);



[[nodiscard]] LocalGasCellMigrationState collectLocalGasCellRecords(
    const core::SimulationState& state,
    std::span<const std::uint32_t> gas_particle_indices) {
  core::legacyRequireParticleBoundGasCellContract(state, "collectLocalGasCellRecords legacy import path");
  if (!state.patches.isConsistent()) {
    throw std::runtime_error("collectLocalGasCellRecords requires consistent retained AMR patch metadata");
  }

  LocalGasCellMigrationState retained{};
  retained.gas_records_by_particle_id.reserve(gas_particle_indices.size());
  retained.patch_metadata_by_patch_id.reserve(state.patches.size());
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  std::vector<std::uint8_t> keep_mask(state.particles.size(), 0U);
  for (const std::uint32_t particle_index : gas_particle_indices) {
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range("gas particle index out of range while collecting gas migration records");
    }
    keep_mask[particle_index] = 1U;
  }

  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint32_t particle_index = gas_globals[cell_index];
    if (keep_mask[particle_index] == 0U) {
      continue;
    }
    const std::uint32_t old_patch_index = state.cells.patch_index[cell_index];
    if (old_patch_index >= state.patches.size()) {
      throw std::runtime_error("retained gas cell refers to missing AMR patch metadata during ownership compaction");
    }
    const std::uint32_t first_cell = state.patches.first_cell[old_patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[old_patch_index];
    if (cell_index < first_cell || cell_index >= static_cast<std::size_t>(first_cell) + cell_count) {
      throw std::runtime_error("retained gas cell lies outside its AMR patch cell range during ownership compaction");
    }
    const std::uint64_t patch_id = state.patches.patch_id[old_patch_index];
    if (patch_id == 0ULL) {
      throw std::runtime_error("retained gas cell refers to an AMR patch with an invalid stable patch ID");
    }

    const LocalPatchMetadata metadata{
        .patch_id = patch_id,
        .level = state.patches.level[old_patch_index],
        .original_cell_count = cell_count,
        .parent_patch_id = state.patches.parent_patch_id[old_patch_index],
        .morton_key = state.patches.morton_key[old_patch_index],
        .origin_x_comoving = state.patches.origin_x_comoving[old_patch_index],
        .origin_y_comoving = state.patches.origin_y_comoving[old_patch_index],
        .origin_z_comoving = state.patches.origin_z_comoving[old_patch_index],
        .extent_x_comoving = state.patches.extent_x_comoving[old_patch_index],
        .extent_y_comoving = state.patches.extent_y_comoving[old_patch_index],
        .extent_z_comoving = state.patches.extent_z_comoving[old_patch_index],
        .cell_dim_x = state.patches.cell_dim_x[old_patch_index],
        .cell_dim_y = state.patches.cell_dim_y[old_patch_index],
        .cell_dim_z = state.patches.cell_dim_z[old_patch_index],
        .owning_rank = state.patches.owning_rank[old_patch_index],
    };
    const auto [metadata_it, metadata_inserted] =
        retained.patch_metadata_by_patch_id.emplace(patch_id, metadata);
    if (!metadata_inserted &&
        (metadata_it->second.level != metadata.level ||
         metadata_it->second.original_cell_count != metadata.original_cell_count ||
         metadata_it->second.parent_patch_id != metadata.parent_patch_id ||
         metadata_it->second.morton_key != metadata.morton_key ||
         metadata_it->second.owning_rank != metadata.owning_rank)) {
      throw std::runtime_error("duplicate stable AMR patch ID has conflicting metadata during ownership compaction");
    }

    GasCellMigrationRecord record;
    record.particle_id = state.particle_sidecar.particle_id[particle_index];
    record.patch_id = patch_id;
    record.patch_local_cell_offset = static_cast<std::uint32_t>(cell_index - first_cell);
    record.center_x_comoving = state.cells.center_x_comoving[cell_index];
    record.center_y_comoving = state.cells.center_y_comoving[cell_index];
    record.center_z_comoving = state.cells.center_z_comoving[cell_index];
    record.mass_code = state.cells.mass_code[cell_index];
    record.time_bin = state.cells.time_bin[cell_index];
    record.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell_index];
    record.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell_index];
    record.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell_index];
    record.density_code = state.gas_cells.density_code[cell_index];
    record.pressure_code = state.gas_cells.pressure_code[cell_index];
    record.internal_energy_code = state.gas_cells.internal_energy_code[cell_index];
    record.temperature_code = state.gas_cells.temperature_code[cell_index];
    record.sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
    const auto [record_it, record_inserted] =
        retained.gas_records_by_particle_id.emplace(record.particle_id, record);
    if (!record_inserted) {
      throw std::runtime_error("duplicate retained gas particle ID during ownership compaction");
    }
  }
  return retained;
}

[[nodiscard]] std::vector<std::uint64_t> gasParticleIdByOldCellIndex(const core::SimulationState& state) {
  core::legacyRequireParticleBoundGasCellContract(state, "gasParticleIdByOldCellIndex legacy import path");
  std::vector<std::uint64_t> cell_particle_ids(state.cells.size(), 0ULL);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    cell_particle_ids[cell_index] =
        core::parentParticleIdForGasCellRow(state, static_cast<std::uint32_t>(cell_index)).value();
  }
  return cell_particle_ids;
}

void rebuildLocalGasStateFromParticleIds(
    core::SimulationState& state,
    const LocalGasCellMigrationState& retained,
    std::span<const std::uint64_t> old_cell_particle_id) {
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  if (gas_globals.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error("rebuildLocalGasStateFromParticleIds exceeds uint32_t cell-index capacity");
  }

  core::CellSoa rebuilt_cells;
  core::GasCellSidecar rebuilt_gas;
  core::PatchSoa rebuilt_patches;
  rebuilt_cells.resize(gas_globals.size());
  rebuilt_gas.resize(gas_globals.size());

  struct RebuiltPatchRange {
    const LocalPatchMetadata* metadata = nullptr;
    std::uint32_t first_cell = 0;
    std::uint32_t cell_count = 0;
    std::uint32_t previous_patch_local_offset = 0;
    bool has_previous_patch_local_offset = false;
  };
  std::vector<RebuiltPatchRange> rebuilt_patch_ranges;
  rebuilt_patch_ranges.reserve(retained.patch_metadata_by_patch_id.size());
  std::unordered_map<std::uint64_t, std::uint32_t> rebuilt_patch_index_by_patch_id;
  rebuilt_patch_index_by_patch_id.reserve(retained.patch_metadata_by_patch_id.size());

  std::unordered_map<std::uint64_t, std::uint32_t> new_cell_index_by_particle_id;
  new_cell_index_by_particle_id.reserve(gas_globals.size());
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint64_t particle_id = state.particle_sidecar.particle_id[gas_globals[cell_index]];
    const auto record_it = retained.gas_records_by_particle_id.find(particle_id);
    if (record_it == retained.gas_records_by_particle_id.end()) {
      throw std::runtime_error("missing gas-cell migration record for local gas particle after ownership compaction");
    }
    const GasCellMigrationRecord& record = record_it->second;
    const auto metadata_it = retained.patch_metadata_by_patch_id.find(record.patch_id);
    if (metadata_it == retained.patch_metadata_by_patch_id.end()) {
      throw std::runtime_error("retained gas cell lost its stable AMR patch metadata during ownership compaction");
    }

    std::uint32_t rebuilt_patch_index = 0;
    const auto patch_index_it = rebuilt_patch_index_by_patch_id.find(record.patch_id);
    if (patch_index_it == rebuilt_patch_index_by_patch_id.end()) {
      if (rebuilt_patch_ranges.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
        throw std::overflow_error("rebuildLocalGasStateFromParticleIds exceeds uint32_t patch-index capacity");
      }
      rebuilt_patch_index = static_cast<std::uint32_t>(rebuilt_patch_ranges.size());
      rebuilt_patch_index_by_patch_id.emplace(record.patch_id, rebuilt_patch_index);
      rebuilt_patch_ranges.push_back(RebuiltPatchRange{
          .metadata = &metadata_it->second,
          .first_cell = static_cast<std::uint32_t>(cell_index),
      });
    } else {
      rebuilt_patch_index = patch_index_it->second;
      if (rebuilt_patch_index + 1U != rebuilt_patch_ranges.size()) {
        throw std::runtime_error(
            "retained gas cells from one AMR patch are non-contiguous after particle compaction; "
            "the legacy particle-bound gas-cell contract cannot preserve this topology");
      }
    }
    RebuiltPatchRange& patch_range = rebuilt_patch_ranges[rebuilt_patch_index];
    if (patch_range.has_previous_patch_local_offset &&
        record.patch_local_cell_offset <= patch_range.previous_patch_local_offset) {
      throw std::runtime_error("retained AMR patch-local cell ordering is not strictly increasing after compaction");
    }
    patch_range.previous_patch_local_offset = record.patch_local_cell_offset;
    patch_range.has_previous_patch_local_offset = true;
    ++patch_range.cell_count;

    rebuilt_gas.gas_cell_id[cell_index] = particle_id;
    rebuilt_gas.parent_particle_id[cell_index] = particle_id;
    rebuilt_cells.center_x_comoving[cell_index] = record.center_x_comoving;
    rebuilt_cells.center_y_comoving[cell_index] = record.center_y_comoving;
    rebuilt_cells.center_z_comoving[cell_index] = record.center_z_comoving;
    rebuilt_cells.mass_code[cell_index] = record.mass_code;
    rebuilt_cells.time_bin[cell_index] = record.time_bin;
    rebuilt_cells.patch_index[cell_index] = rebuilt_patch_index;
    rebuilt_gas.velocity_x_peculiar[cell_index] = record.velocity_x_peculiar;
    rebuilt_gas.velocity_y_peculiar[cell_index] = record.velocity_y_peculiar;
    rebuilt_gas.velocity_z_peculiar[cell_index] = record.velocity_z_peculiar;
    rebuilt_gas.density_code[cell_index] = record.density_code;
    rebuilt_gas.pressure_code[cell_index] = record.pressure_code;
    rebuilt_gas.internal_energy_code[cell_index] = record.internal_energy_code;
    rebuilt_gas.temperature_code[cell_index] = record.temperature_code;
    rebuilt_gas.sound_speed_code[cell_index] = record.sound_speed_code;
    const auto [index_it, inserted] =
        new_cell_index_by_particle_id.emplace(particle_id, static_cast<std::uint32_t>(cell_index));
    if (!inserted) {
      throw std::runtime_error("duplicate local gas particle ID after ownership compaction");
    }
  }

  rebuilt_patches.resize(rebuilt_patch_ranges.size());
  for (std::size_t patch_index = 0; patch_index < rebuilt_patch_ranges.size(); ++patch_index) {
    const RebuiltPatchRange& range = rebuilt_patch_ranges[patch_index];
    if (range.metadata == nullptr || range.metadata->patch_id == 0ULL || range.cell_count == 0U) {
      throw std::runtime_error("invalid retained AMR patch range while rebuilding local gas state");
    }
    if (range.metadata->hasExplicitCartesianGeometry() &&
        range.cell_count != range.metadata->original_cell_count) {
      throw std::runtime_error(
          "ownership compaction would split a geometry-bearing AMR patch; migrate complete patches or preserve "
          "explicit sparse-cell topology before rebuilding local gas state");
    }
    const LocalPatchMetadata& metadata = *range.metadata;
    rebuilt_patches.patch_id[patch_index] = metadata.patch_id;
    rebuilt_patches.level[patch_index] = metadata.level;
    rebuilt_patches.first_cell[patch_index] = range.first_cell;
    rebuilt_patches.cell_count[patch_index] = range.cell_count;
    rebuilt_patches.parent_patch_id[patch_index] = metadata.parent_patch_id;
    rebuilt_patches.morton_key[patch_index] = metadata.morton_key;
    rebuilt_patches.origin_x_comoving[patch_index] = metadata.origin_x_comoving;
    rebuilt_patches.origin_y_comoving[patch_index] = metadata.origin_y_comoving;
    rebuilt_patches.origin_z_comoving[patch_index] = metadata.origin_z_comoving;
    rebuilt_patches.extent_x_comoving[patch_index] = metadata.extent_x_comoving;
    rebuilt_patches.extent_y_comoving[patch_index] = metadata.extent_y_comoving;
    rebuilt_patches.extent_z_comoving[patch_index] = metadata.extent_z_comoving;
    rebuilt_patches.cell_dim_x[patch_index] = metadata.cell_dim_x;
    rebuilt_patches.cell_dim_y[patch_index] = metadata.cell_dim_y;
    rebuilt_patches.cell_dim_z[patch_index] = metadata.cell_dim_z;
    rebuilt_patches.owning_rank[patch_index] = metadata.owning_rank;
  }

  if (state.gas_composition.enabled()) {
    throw std::runtime_error(
        "rebuildLocalGasStateFromParticleIds (legacy particle-bound gas import) cannot carry the CHUI "
        "gas composition sidecar; use gas-cell-ID migration (commitGasCellMigration/commitAmrPatchMigration)");
  }
  state.cells = std::move(rebuilt_cells);
  state.gas_cells = std::move(rebuilt_gas);
  state.patches = std::move(rebuilt_patches);
  state.bumpCellIndexGeneration();
  if (!state.patches.isConsistent()) {
    throw std::runtime_error("rebuildLocalGasStateFromParticleIds produced inconsistent AMR patch metadata");
  }
  state.refreshGasCellIdentityMapFromSidecarLanes();
  core::legacyRequireParticleBoundGasCellContract(state, "rebuildLocalGasStateFromParticleIds legacy import path");

  const auto remap_host_cell = [&](std::uint32_t old_cell_index) {
    if (old_cell_index >= old_cell_particle_id.size()) {
      throw std::runtime_error(
          "rebuildLocalGasStateFromParticleIds: sidecar host_cell_index does not refer to an old gas cell");
    }
    const std::uint64_t host_particle_id = old_cell_particle_id[old_cell_index];
    const auto found = new_cell_index_by_particle_id.find(host_particle_id);
    if (found == new_cell_index_by_particle_id.end()) {
      throw std::runtime_error(
          "rebuildLocalGasStateFromParticleIds: sidecar host gas cell was removed during ownership compaction");
    }
    return found->second;
  };

  for (std::size_t row = 0; row < state.black_holes.size(); ++row) {
    state.black_holes.host_cell_index[row] = remap_host_cell(state.black_holes.host_cell_index[row]);
  }
  for (std::size_t row = 0; row < state.tracers.size(); ++row) {
    state.tracers.host_cell_index[row] = remap_host_cell(state.tracers.host_cell_index[row]);
  }
}

void compactStateToCurrentOwner(
    core::SimulationState& state,
    int world_rank) {
  if (world_rank < 0) {
    throw std::invalid_argument("compactStateToCurrentOwner requires non-negative world rank");
  }
  std::vector<std::uint32_t> stale_ghost_indices;
  stale_ghost_indices.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const bool owned_here = state.particleOwningRank(particle_index) == static_cast<std::uint32_t>(world_rank);
    if (!owned_here) {
      stale_ghost_indices.push_back(static_cast<std::uint32_t>(particle_index));
    }
  }
  core::ParticleMigrationCommit commit;
  commit.world_rank = world_rank;
  commit.stale_local_ghost_indices = std::move(stale_ghost_indices);
  commit.preserve_gas_cell_state = true;
  state.commitParticleMigration(commit);

  if (state.patches.size() != 0 || state.cells.size() != 0) {
    if (!state.patches.isConsistent()) {
      throw std::runtime_error("compactStateToCurrentOwner requires consistent AMR patch metadata");
    }
    std::vector<std::uint32_t> stale_patch_indices;
    stale_patch_indices.reserve(state.patches.size());
    for (std::uint32_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      if (state.patches.owning_rank[patch_index] != static_cast<std::uint32_t>(world_rank)) {
        stale_patch_indices.push_back(patch_index);
      }
    }
    if (!stale_patch_indices.empty()) {
      core::AmrPatchMigrationCommit patch_commit;
      patch_commit.world_rank = world_rank;
      patch_commit.stale_local_ghost_patch_indices = std::move(stale_patch_indices);
      state.commitAmrPatchMigration(patch_commit);
    }
  }
}


[[nodiscard]] parallel::DecompositionConfig makeWorkflowDecompositionConfig(
    const core::SimulationConfig& config,
    int world_size) {
  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = world_size;
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = config.cosmology.box_size_x_mpc_comoving;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = config.cosmology.box_size_y_mpc_comoving;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = config.cosmology.box_size_z_mpc_comoving;
  decomposition_config.owned_particle_weight = 0.0;
  decomposition_config.active_target_weight = 0.0;
  decomposition_config.remote_tree_interaction_weight = 0.0;
  decomposition_config.work_weight = 0.0;
  decomposition_config.memory_weight = 0.0;
  decomposition_config.component_weights = parallel::DecompositionWeightCoefficients{
      .particle_count = config.parallel.decomposition_particle_count_weight,
      .gas_cell = config.parallel.decomposition_gas_cell_weight,
      .tree_interaction = config.parallel.decomposition_tree_interaction_weight,
      .pm_mesh = config.parallel.decomposition_pm_mesh_weight,
      .amr_patch = config.parallel.decomposition_amr_patch_weight,
      .active_fraction = config.parallel.decomposition_active_fraction_weight,
      .memory_pressure = config.parallel.decomposition_memory_pressure_weight,
      .gpu_occupancy = config.parallel.decomposition_gpu_occupancy_weight,
      .generic_work = config.parallel.decomposition_generic_work_weight,
  };
  return decomposition_config;
}

[[nodiscard]] std::uint64_t declaredMajorTaskTransientReserveBytes(
    const core::SimulationConfig& config) {
  const std::uint64_t hydro_cells = config.numerics.hydro_active_batch_max_cells == 0U
      ? static_cast<std::uint64_t>(hydro::k_hydro_automatic_active_batch_max_cells)
      : config.numerics.hydro_active_batch_max_cells;
  const std::size_t hydro_bytes = core::checkedSizeMultiply(
      core::checkedIntegralNarrow<std::size_t>(hydro_cells, "hydro batch cell count"),
      core::checkedIntegralNarrow<std::size_t>(hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell,
                                               "hydro batch scratch bytes per cell"),
      "hydro major-task transient reserve");
  const std::size_t communication_bytes = core::checkedSizeMultiply(
      parallel::mpiTransportRoundLimitBytes(), 2U,
      "bidirectional communication major-task transient reserve");
  constexpr std::size_t k_feedback_batch_max = 4096U;
  const std::size_t feedback_bytes = core::checkedSizeMultiply(
      k_feedback_batch_max,
      sizeof(physics::StellarEvolutionStarBudget) + sizeof(physics::StellarFeedbackEvent) + sizeof(std::uint32_t),
      "feedback major-task transient reserve");
  return core::checkedIntegralNarrow<std::uint64_t>(
      std::max({hydro_bytes, communication_bytes, feedback_bytes}),
      "major-task transient reserve byte width");
}

void applyMemoryAwareDecompositionEnvelope(
    parallel::DecompositionConfig& decomposition_config,
    std::uint64_t local_persistent_bytes,
    const core::SimulationConfig& config,
    const RuntimeServices& services) {
  if (services.memory_governor == nullptr) {
    return;
  }
  const core::MemoryGovernorSnapshot snapshot = services.memory_governor->snapshot();
  if (snapshot.hard_limit_bytes == 0U) {
    return;
  }
  if (local_persistent_bytes > std::numeric_limits<std::uint64_t>::max() - snapshot.headroom_bytes) {
    throw std::overflow_error("decomposition headroom plus persistent estimate overflows uint64");
  }
  const std::uint64_t local_ceiling = local_persistent_bytes + snapshot.headroom_bytes;
  const std::uint64_t global_ceiling = services.mpi_context.allreduceMinUint64(local_ceiling);
  const std::uint64_t transient_reserve = declaredMajorTaskTransientReserveBytes(config);
  if (global_ceiling <= transient_reserve) {
    throw std::runtime_error(
        "memory-aware decomposition has no rank-local headroom for the declared major-task transient reserve");
  }
  decomposition_config.max_rank_memory_bytes = global_ceiling;
  decomposition_config.rank_transient_reserve_bytes = transient_reserve;
}

[[nodiscard]] parallel::DecompositionFeedbackCoefficients makeWorkflowFeedbackCoefficients(
    const core::SimulationConfig& config) {
  return parallel::DecompositionFeedbackCoefficients{
      .measured_tree_pair = config.parallel.decomposition_measured_tree_pair_weight,
      .measured_pm_cell = config.parallel.decomposition_measured_pm_cell_weight,
      .measured_amr_cell = config.parallel.decomposition_measured_amr_cell_weight,
      .measured_hydro_face = config.parallel.decomposition_measured_hydro_face_weight,
      .measured_wall_ms = config.parallel.decomposition_measured_wall_ms_weight,
  };
}

[[maybe_unused]] [[nodiscard]] std::vector<parallel::DecompositionItem> buildRuntimeDecompositionItems(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_rank,
    std::span<const std::uint32_t> active_particle_indices) {
  std::vector<std::uint8_t> active_mask(state.particles.size(), 0U);
  for (const std::uint32_t pidx : active_particle_indices) {
    if (pidx < active_mask.size()) {
      active_mask[pidx] = 1U;
    }
  }
  std::vector<std::uint32_t> patch_cell_count(state.patches.size(), 0U);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_index];
    if (patch_index < patch_cell_count.size()) {
      ++patch_cell_count[patch_index];
    }
  }

  std::vector<parallel::DecompositionItem> items;
  items.reserve(state.particles.size() + state.patches.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const std::uint32_t species_tag = state.particleSpeciesTag(particle_index);
    parallel::DecompositionItem item;
    item.entity_id = state.particle_sidecar.particle_id[particle_index];
    item.kind = parallel::DecompositionEntityKind::kParticle;
    item.current_owner_rank = state.particle_sidecar.owning_rank.empty()
        ? world_rank
        : static_cast<int>(state.particleOwningRank(particle_index));
    item.x_comov = state.particles.position_x_comoving[particle_index];
    item.y_comov = state.particles.position_y_comoving[particle_index];
    item.z_comov = state.particles.position_z_comoving[particle_index];
    item.has_spatial_bounds = true;
    item.min_x_comov = item.max_x_comov = item.x_comov;
    item.min_y_comov = item.max_y_comov = item.y_comov;
    item.min_z_comov = item.max_z_comov = item.z_comov;
    item.active_target_count_recent = active_mask[particle_index] != 0U ? 1U : 0U;
    item.remote_tree_interactions_recent = 1U;
    item.memory_bytes = estimateParticleMemoryBytesForDecomposition(state, species_tag);
    double amr_patch_cost = 0.0;
    if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) && !state.cells.patch_index.empty()) {
      std::vector<std::uint32_t> seen_patch_indices;
      for (const std::uint32_t cell_row : state.gas_cell_identity.rowsForParentParticleId(item.entity_id)) {
        if (cell_row >= state.cells.patch_index.size()) {
          continue;
        }
        const std::uint32_t patch_index = state.cells.patch_index[cell_row];
        if (patch_index >= patch_cell_count.size() ||
            std::find(seen_patch_indices.begin(), seen_patch_indices.end(), patch_index) != seen_patch_indices.end()) {
          continue;
        }
        seen_patch_indices.push_back(patch_index);
        amr_patch_cost += static_cast<double>(patch_cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max<std::int32_t>(state.patches.level[patch_index], 0)));
      }
    }
    item.work_components = parallel::DecompositionWorkComponents{
        .particle_count_cost = 1.0,
        .gas_cell_cost = species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) ? 1.0 : 0.0,
        .tree_interaction_cost = 1.0,
        .pm_mesh_cost = 1.0,
        .amr_patch_cost = amr_patch_cost,
        .active_fraction_cost = static_cast<double>(item.active_target_count_recent),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .transient_memory_cost = (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) &&
                                  item.active_target_count_recent != 0U)
            ? static_cast<double>(hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell)
            : 0.0,
        .source_event_cost = (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar) ||
                              species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole))
            ? 1.0
            : 0.0,
        .communication_cost = static_cast<double>(item.remote_tree_interactions_recent) + 1.0,
        .gpu_occupancy_cost = 0.0,
        .generic_work_cost = 1.0,
        .has_explicit_components = true,
    };
    items.push_back(item);
  }

  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.cell_count[patch_index] == 0U) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    if (first_cell >= state.cells.size()) {
      continue;
    }
    parallel::DecompositionItem item;
    item.entity_id = state.patches.patch_id[patch_index];
    item.kind = parallel::DecompositionEntityKind::kAmrPatch;
    item.current_owner_rank = state.patches.owning_rank.empty()
        ? world_rank
        : static_cast<int>(state.patches.owning_rank[patch_index]);
    item.x_comov = state.cells.center_x_comoving[first_cell];
    item.y_comov = state.cells.center_y_comoving[first_cell];
    item.z_comov = state.cells.center_z_comoving[first_cell];
    item.has_spatial_bounds = true;
    item.min_x_comov = item.max_x_comov = item.x_comov;
    item.min_y_comov = item.max_y_comov = item.y_comov;
    item.min_z_comov = item.max_z_comov = item.z_comov;
    const std::uint64_t patch_end =
        static_cast<std::uint64_t>(first_cell) +
        static_cast<std::uint64_t>(state.patches.cell_count[patch_index]);
    if (patch_end > state.cells.size()) {
      throw std::runtime_error(
          "runtime decomposition patch cell range is outside authoritative CellSoa");
    }
    for (std::uint64_t cell = first_cell; cell < patch_end; ++cell) {
      const std::size_t cell_index = static_cast<std::size_t>(cell);
      item.min_x_comov = std::min(item.min_x_comov, state.cells.center_x_comoving[cell_index]);
      item.max_x_comov = std::max(item.max_x_comov, state.cells.center_x_comoving[cell_index]);
      item.min_y_comov = std::min(item.min_y_comov, state.cells.center_y_comoving[cell_index]);
      item.max_y_comov = std::max(item.max_y_comov, state.cells.center_y_comoving[cell_index]);
      item.min_z_comov = std::min(item.min_z_comov, state.cells.center_z_comoving[cell_index]);
      item.max_z_comov = std::max(item.max_z_comov, state.cells.center_z_comoving[cell_index]);
    }
    item.memory_bytes = static_cast<std::uint64_t>(state.patches.cell_count[patch_index]) *
        static_cast<std::uint64_t>(sizeof(double) * 8U + sizeof(std::uint32_t) * 2U);
    item.work_components = parallel::DecompositionWorkComponents{
        .amr_patch_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max<std::int32_t>(state.patches.level[patch_index], 0))),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .transient_memory_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
            static_cast<double>(hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell),
        .communication_cost = 2.0 * (
            static_cast<double>(state.patches.cell_dim_x[patch_index]) * state.patches.cell_dim_y[patch_index] +
            static_cast<double>(state.patches.cell_dim_x[patch_index]) * state.patches.cell_dim_z[patch_index] +
            static_cast<double>(state.patches.cell_dim_y[patch_index]) * state.patches.cell_dim_z[patch_index]),
        .generic_work_cost = static_cast<double>(state.patches.cell_count[patch_index]),
        .has_explicit_components = true,
    };
    items.push_back(item);
  }
  parallel::applyRuntimeDecompositionFeedback(items, parallel::DecompositionRuntimeMeasurements{}, makeWorkflowFeedbackCoefficients(config));
  return items;
}


void syncTimeBinsFromScheduler(
    const core::HierarchicalTimeBinScheduler& scheduler,
    core::SimulationState& state);




struct MigrationExchangeStats {
  std::uint64_t packet_send_capacity_bytes = 0U;
  std::uint64_t packet_receive_capacity_bytes = 0U;
  std::uint64_t record_assembly_capacity_bytes = 0U;
  std::uint64_t communication_high_water_bytes = 0U;
  std::uint64_t wire_bytes_sent = 0U;
  std::uint64_t wire_bytes_received = 0U;
  std::uint64_t physical_exchange_rounds = 0U;
};

template <typename TRecord>
struct MigrationExchangeResult {
  std::vector<TRecord> records;
  MigrationExchangeStats stats{};
};

template <typename TRecord, typename TEncode, typename TDecode>
[[nodiscard]] MigrationExchangeResult<TRecord> exchangeRuntimeMigrationPayloads(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<TRecord>>& records_by_rank,
    std::string_view phase_name,
    TEncode&& encode_record,
    TDecode&& decode_record) {
  const int world_size = mpi_context.worldSize();
  if (world_size <= 0 || records_by_rank.size() != static_cast<std::size_t>(world_size)) {
    throw std::invalid_argument("migration packet exchange rank extent must match positive MPI world size");
  }
  if (world_size == 1) {
    return MigrationExchangeResult<TRecord>{
        .records = records_by_rank.empty() ? std::vector<TRecord>{} : records_by_rank.front(),
        .stats = {},
    };
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("runtime migration packet exchange requires MPI for world_size > 1");
  }
#if COSMOSIM_ENABLE_MPI
  struct OutboundCursor {
    std::size_t record_index = 0U;
    std::size_t fragment_offset = 0U;
    std::vector<std::uint8_t> encoded_record;
  };
  struct InboundAssembly {
    std::uint64_t record_sequence = 0U;
    std::uint64_t total_bytes = 0U;
    std::size_t received_bytes = 0U;
    std::vector<std::uint8_t> encoded_record;
    bool active = false;
  };

  const std::size_t round_limit = parallel::mpiTransportRoundLimitBytes();
  const std::size_t rank_count = static_cast<std::size_t>(world_size);
  const migration_wire::PacketCapacityPlan packet_plan =
      migration_wire::planPacketCapacity(round_limit, rank_count);
  const std::size_t fragment_payload_limit = packet_plan.fragment_payload_bytes;

  std::vector<OutboundCursor> outbound(rank_count);
  std::vector<InboundAssembly> inbound_assembly(rank_count);
  MigrationExchangeResult<TRecord> result;
  std::vector<std::vector<std::uint8_t>> send_payloads(rank_count);

  const auto localHasPending = [&]() {
    for (std::size_t peer = 0; peer < rank_count; ++peer) {
      if (outbound[peer].record_index < records_by_rank[peer].size()) return true;
    }
    return false;
  };
  const auto assemblyCapacityBytes = [&]() {
    std::uint64_t total = 0U;
    for (const auto& assembly : inbound_assembly) {
      total = core::checkedMemoryBytesAdd(
          total,
          core::checkedIntegralNarrow<std::uint64_t>(
              assembly.encoded_record.capacity(), "migration record assembly capacity"),
          "migration record assembly aggregate capacity");
    }
    return total;
  };
  const auto payloadCapacityBytes = [](const std::vector<std::vector<std::uint8_t>>& payloads) {
    std::uint64_t total = 0U;
    for (const auto& payload : payloads) {
      total = core::checkedMemoryBytesAdd(
          total,
          core::checkedIntegralNarrow<std::uint64_t>(payload.capacity(), "migration packet capacity"),
          "migration packet aggregate capacity");
    }
    return total;
  };

  while (mpi_context.allreduceMaxUint64(localHasPending() ? 1U : 0U) != 0U) {
    std::exception_ptr local_failure;
    try {
      for (auto& payload : send_payloads) payload.clear();
      for (std::size_t peer = 0; peer < rank_count; ++peer) {
        auto& cursor = outbound[peer];
        if (cursor.record_index >= records_by_rank[peer].size()) continue;
        if (cursor.encoded_record.empty()) {
          cursor.encoded_record = encode_record(records_by_rank[peer][cursor.record_index]);
          cursor.fragment_offset = 0U;
          if (cursor.encoded_record.empty()) {
            throw std::runtime_error("migration record encoder returned an empty record");
          }
        }
        const std::size_t remaining = cursor.encoded_record.size() - cursor.fragment_offset;
        const std::size_t fragment_bytes = std::min(fragment_payload_limit, remaining);
        send_payloads[peer] = migration_wire::encodeFragment(
            static_cast<std::uint64_t>(cursor.record_index),
            static_cast<std::uint64_t>(cursor.encoded_record.size()),
            static_cast<std::uint64_t>(cursor.fragment_offset),
            std::span<const std::uint8_t>(
                cursor.encoded_record.data() + cursor.fragment_offset,
                fragment_bytes));
        result.stats.wire_bytes_sent = core::checkedMemoryBytesAdd(
            result.stats.wire_bytes_sent,
            core::checkedIntegralNarrow<std::uint64_t>(
                send_payloads[peer].size(), "migration packet sent bytes"),
            "migration cumulative wire bytes sent");
        cursor.fragment_offset = core::checkedSizeAdd(
            cursor.fragment_offset, fragment_bytes, "migration outbound fragment progress");
        if (cursor.fragment_offset == cursor.encoded_record.size()) {
          cursor.encoded_record.clear();
          ++cursor.record_index;
          cursor.fragment_offset = 0U;
        }
      }
    } catch (...) {
      local_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(local_failure, phase_name);

    result.stats.packet_send_capacity_bytes = std::max(
        result.stats.packet_send_capacity_bytes, payloadCapacityBytes(send_payloads));
    const auto recv_payloads = parallel::exchangeBoundedAlltoallBytes(mpi_context, send_payloads);
    result.stats.packet_receive_capacity_bytes = std::max(
        result.stats.packet_receive_capacity_bytes, payloadCapacityBytes(recv_payloads));
    ++result.stats.physical_exchange_rounds;

    local_failure = nullptr;
    try {
      for (std::size_t peer = 0; peer < rank_count; ++peer) {
        const auto& payload = recv_payloads[peer];
        if (payload.empty()) continue;
        result.stats.wire_bytes_received = core::checkedMemoryBytesAdd(
            result.stats.wire_bytes_received,
            core::checkedIntegralNarrow<std::uint64_t>(payload.size(), "migration packet received bytes"),
            "migration cumulative wire bytes received");
        const migration_wire::FragmentView fragment = migration_wire::decodeFragment(payload);
        auto& assembly = inbound_assembly[peer];
        if (!assembly.active) {
          if (fragment.fragment_offset != 0U) {
            throw std::runtime_error("migration record fragment stream did not begin at offset zero");
          }
          assembly.active = true;
          assembly.record_sequence = fragment.record_sequence;
          assembly.total_bytes = fragment.record_total_bytes;
          assembly.received_bytes = 0U;
          assembly.encoded_record.resize(core::checkedIntegralNarrow<std::size_t>(
              fragment.record_total_bytes, "migration inbound record bytes"));
        }
        if (fragment.record_sequence != assembly.record_sequence ||
            fragment.record_total_bytes != assembly.total_bytes ||
            fragment.fragment_offset != assembly.received_bytes) {
          throw std::runtime_error("migration fragment sequence/offset contract mismatch");
        }
        if (fragment.payload.size() > assembly.encoded_record.size() - assembly.received_bytes) {
          throw std::runtime_error("migration fragment exceeds inbound record assembly capacity");
        }
        std::copy(
            fragment.payload.begin(), fragment.payload.end(),
            assembly.encoded_record.begin() + static_cast<std::ptrdiff_t>(assembly.received_bytes));
        assembly.received_bytes = core::checkedSizeAdd(
            assembly.received_bytes, fragment.payload.size(), "migration inbound fragment progress");
        if (assembly.received_bytes == assembly.encoded_record.size()) {
          result.records.push_back(decode_record(std::span<const std::uint8_t>(assembly.encoded_record)));
          assembly.encoded_record.clear();
          assembly.active = false;
          ++assembly.record_sequence;
          assembly.total_bytes = 0U;
          assembly.received_bytes = 0U;
        }
      }
      result.stats.record_assembly_capacity_bytes = std::max(
          result.stats.record_assembly_capacity_bytes, assemblyCapacityBytes());
      const std::uint64_t packet_capacity = core::checkedMemoryBytesAdd(
          result.stats.packet_send_capacity_bytes,
          result.stats.packet_receive_capacity_bytes,
          "migration packet send/receive capacity");
      result.stats.communication_high_water_bytes = std::max(
          result.stats.communication_high_water_bytes,
          core::checkedMemoryBytesAdd(
              packet_capacity,
              result.stats.record_assembly_capacity_bytes,
              "migration packet/assembly high water"));
    } catch (...) {
      local_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(local_failure, phase_name);
  }

  std::exception_ptr completion_failure;
  try {
    for (const auto& assembly : inbound_assembly) {
      if (assembly.active || !assembly.encoded_record.empty()) {
        throw std::runtime_error("migration packet exchange ended with an incomplete record assembly");
      }
    }
  } catch (...) {
    completion_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(completion_failure, phase_name);
  return result;
#else
  throw std::runtime_error("runtime migration packet exchange requires an MPI-enabled build");
#endif
}


struct CompactDmoExchangeResult {
  MigrationExchangeStats stats{};
  std::vector<std::uint64_t> inbound_particle_ids;
};

[[nodiscard]] CompactDmoExchangeResult exchangeRuntimeCompactDmoParticles(
    const parallel::MpiContext& mpi_context,
    const core::SimulationState& state,
    const std::vector<std::vector<std::uint32_t>>& local_indices_by_rank,
    std::span<const std::uint64_t> inbound_count_by_rank,
    core::ParticleSoa& candidate_particles,
    core::ParticleSidecar& candidate_sidecar,
    std::size_t inbound_begin,
    bool materialize_candidate_sfc) {
  const int world_size = mpi_context.worldSize();
  const std::size_t rank_count = static_cast<std::size_t>(world_size);
  if (world_size <= 0 || local_indices_by_rank.size() != rank_count ||
      inbound_count_by_rank.size() != rank_count) {
    throw std::invalid_argument("compact DMO migration exchange rank extent mismatch");
  }
  CompactDmoExchangeResult result;
  std::size_t total_inbound = 0U;
  std::vector<std::size_t> inbound_peer_begin(rank_count, inbound_begin);
  for (std::size_t peer = 0; peer < rank_count; ++peer) {
    inbound_peer_begin[peer] = core::checkedSizeAdd(inbound_begin, total_inbound,
        "compact DMO inbound peer begin");
    total_inbound = core::checkedSizeAdd(total_inbound,
        core::checkedIntegralNarrow<std::size_t>(inbound_count_by_rank[peer],
            "compact DMO inbound peer count"),
        "compact DMO total inbound count");
  }
  if (core::checkedSizeAdd(inbound_begin, total_inbound, "compact DMO candidate inbound extent") !=
      candidate_particles.size()) {
    throw std::logic_error("compact DMO migration candidate inbound extent mismatch");
  }
  result.inbound_particle_ids.reserve(total_inbound);
  if (world_size == 1) {
    if (total_inbound != 0U || !local_indices_by_rank.front().empty()) {
      throw std::logic_error("single-rank compact DMO migration cannot contain remote records");
    }
    return result;
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("compact DMO migration exchange requires MPI for world_size > 1");
  }
#if COSMOSIM_ENABLE_MPI
  struct OutboundCursor {
    std::size_t record_index = 0U;
    std::size_t fragment_offset = 0U;
    std::vector<std::uint8_t> encoded_record;
  };
  struct InboundAssembly {
    std::uint64_t record_sequence = 0U;
    std::uint64_t total_bytes = 0U;
    std::size_t received_bytes = 0U;
    std::vector<std::uint8_t> encoded_record;
    bool active = false;
  };
  const std::size_t round_limit = parallel::mpiTransportRoundLimitBytes();
  const migration_wire::PacketCapacityPlan packet_plan =
      migration_wire::planPacketCapacity(round_limit, rank_count);
  const std::size_t fragment_payload_limit = packet_plan.fragment_payload_bytes;
  std::vector<OutboundCursor> outbound(rank_count);
  std::vector<InboundAssembly> inbound_assembly(rank_count);
  std::vector<std::vector<std::uint8_t>> send_payloads(rank_count);
  const auto localHasPending = [&]() {
    for (std::size_t peer = 0; peer < rank_count; ++peer) {
      if (outbound[peer].record_index < local_indices_by_rank[peer].size()) return true;
    }
    return false;
  };
  const auto payloadCapacityBytes = [](const std::vector<std::vector<std::uint8_t>>& payloads) {
    std::uint64_t total = 0U;
    for (const auto& payload : payloads) {
      total = core::checkedMemoryBytesAdd(total,
          core::checkedIntegralNarrow<std::uint64_t>(payload.capacity(), "compact DMO packet capacity"),
          "compact DMO packet aggregate capacity");
    }
    return total;
  };
  const auto assemblyCapacityBytes = [&]() {
    std::uint64_t total = 0U;
    for (const auto& assembly : inbound_assembly) {
      total = core::checkedMemoryBytesAdd(total,
          core::checkedIntegralNarrow<std::uint64_t>(assembly.encoded_record.capacity(),
              "compact DMO record assembly capacity"),
          "compact DMO record assembly aggregate capacity");
    }
    return total;
  };
  while (mpi_context.allreduceMaxUint64(localHasPending() ? 1U : 0U) != 0U) {
    std::exception_ptr local_failure;
    try {
      for (auto& payload : send_payloads) payload.clear();
      for (std::size_t peer = 0; peer < rank_count; ++peer) {
        auto& cursor = outbound[peer];
        if (cursor.record_index >= local_indices_by_rank[peer].size()) continue;
        if (cursor.encoded_record.empty()) {
          const std::uint32_t row = local_indices_by_rank[peer][cursor.record_index];
          if (row >= state.particles.size()) throw std::out_of_range("compact DMO outbound row out of range");
          cursor.encoded_record = migration_wire::encodeDmoParticleMigrationRecord(
              migration_wire::DmoParticleMigrationRecord{
                  .particle_id = state.particle_sidecar.particle_id[row],
                  .sfc_key = state.particleSfcKey(row),
                  .position_x_comoving = state.particles.position_x_comoving[row],
                  .position_y_comoving = state.particles.position_y_comoving[row],
                  .position_z_comoving = state.particles.position_z_comoving[row],
                  .velocity_x_peculiar = state.particles.velocity_x_peculiar[row],
                  .velocity_y_peculiar = state.particles.velocity_y_peculiar[row],
                  .velocity_z_peculiar = state.particles.velocity_z_peculiar[row],
                  .mass_code = state.particles.mass_code[row],
              });
          cursor.fragment_offset = 0U;
        }
        const std::size_t remaining = cursor.encoded_record.size() - cursor.fragment_offset;
        const std::size_t fragment_bytes = std::min(fragment_payload_limit, remaining);
        send_payloads[peer] = migration_wire::encodeFragment(
            static_cast<std::uint64_t>(cursor.record_index),
            static_cast<std::uint64_t>(cursor.encoded_record.size()),
            static_cast<std::uint64_t>(cursor.fragment_offset),
            std::span<const std::uint8_t>(cursor.encoded_record.data() + cursor.fragment_offset,
                                          fragment_bytes));
        result.stats.wire_bytes_sent = core::checkedMemoryBytesAdd(
            result.stats.wire_bytes_sent,
            core::checkedIntegralNarrow<std::uint64_t>(send_payloads[peer].size(),
                "compact DMO sent packet bytes"),
            "compact DMO cumulative wire bytes sent");
        cursor.fragment_offset = core::checkedSizeAdd(cursor.fragment_offset, fragment_bytes,
                                                       "compact DMO outbound fragment progress");
        if (cursor.fragment_offset == cursor.encoded_record.size()) {
          cursor.encoded_record.clear();
          ++cursor.record_index;
          cursor.fragment_offset = 0U;
        }
      }
    } catch (...) { local_failure = std::current_exception(); }
    mpi_context.rethrowCollectivePreparationFailure(local_failure,
        "compact DMO migration bounded packet preparation");
    result.stats.packet_send_capacity_bytes = std::max(result.stats.packet_send_capacity_bytes,
                                                        payloadCapacityBytes(send_payloads));
    const auto recv_payloads = parallel::exchangeBoundedAlltoallBytes(mpi_context, send_payloads);
    result.stats.packet_receive_capacity_bytes = std::max(result.stats.packet_receive_capacity_bytes,
                                                           payloadCapacityBytes(recv_payloads));
    ++result.stats.physical_exchange_rounds;
    local_failure = nullptr;
    try {
      for (std::size_t peer = 0; peer < rank_count; ++peer) {
        const auto& payload = recv_payloads[peer];
        if (payload.empty()) continue;
        result.stats.wire_bytes_received = core::checkedMemoryBytesAdd(
            result.stats.wire_bytes_received,
            core::checkedIntegralNarrow<std::uint64_t>(payload.size(), "compact DMO received packet bytes"),
            "compact DMO cumulative wire bytes received");
        const migration_wire::FragmentView fragment = migration_wire::decodeFragment(payload);
        auto& assembly = inbound_assembly[peer];
        if (!assembly.active) {
          if (fragment.fragment_offset != 0U) throw std::runtime_error("compact DMO fragment stream did not start at zero");
          if (fragment.record_sequence >= inbound_count_by_rank[peer]) throw std::runtime_error("compact DMO inbound sequence exceeds admitted count");
          assembly.active = true;
          assembly.record_sequence = fragment.record_sequence;
          assembly.total_bytes = fragment.record_total_bytes;
          assembly.received_bytes = 0U;
          assembly.encoded_record.resize(core::checkedIntegralNarrow<std::size_t>(
              fragment.record_total_bytes, "compact DMO inbound record bytes"));
        }
        if (fragment.record_sequence != assembly.record_sequence ||
            fragment.record_total_bytes != assembly.total_bytes ||
            fragment.fragment_offset != assembly.received_bytes) {
          throw std::runtime_error("compact DMO fragment sequence/offset contract mismatch");
        }
        if (fragment.payload.size() > assembly.encoded_record.size() - assembly.received_bytes) {
          throw std::runtime_error("compact DMO fragment exceeds record assembly capacity");
        }
        std::copy(fragment.payload.begin(), fragment.payload.end(),
                  assembly.encoded_record.begin() + static_cast<std::ptrdiff_t>(assembly.received_bytes));
        assembly.received_bytes = core::checkedSizeAdd(assembly.received_bytes, fragment.payload.size(),
                                                        "compact DMO inbound fragment progress");
        if (assembly.received_bytes == assembly.encoded_record.size()) {
          const auto record = migration_wire::decodeDmoParticleMigrationRecord(assembly.encoded_record);
          if (!std::isfinite(record.position_x_comoving) || !std::isfinite(record.position_y_comoving) ||
              !std::isfinite(record.position_z_comoving) || !std::isfinite(record.velocity_x_peculiar) ||
              !std::isfinite(record.velocity_y_peculiar) || !std::isfinite(record.velocity_z_peculiar) ||
              !std::isfinite(record.mass_code) || record.mass_code < 0.0) {
            throw std::runtime_error("compact DMO migration received non-finite/negative particle state");
          }
          const std::size_t destination = core::checkedSizeAdd(
              inbound_peer_begin[peer],
              core::checkedIntegralNarrow<std::size_t>(assembly.record_sequence,
                  "compact DMO inbound sequence width"),
              "compact DMO inbound destination row");
          if (destination >= candidate_particles.size()) throw std::runtime_error("compact DMO inbound destination outside candidate");
          candidate_sidecar.particle_id[destination] = record.particle_id;
          if (materialize_candidate_sfc) candidate_sidecar.sfc_key[destination] = record.sfc_key;
          candidate_particles.position_x_comoving[destination] = record.position_x_comoving;
          candidate_particles.position_y_comoving[destination] = record.position_y_comoving;
          candidate_particles.position_z_comoving[destination] = record.position_z_comoving;
          candidate_particles.velocity_x_peculiar[destination] = record.velocity_x_peculiar;
          candidate_particles.velocity_y_peculiar[destination] = record.velocity_y_peculiar;
          candidate_particles.velocity_z_peculiar[destination] = record.velocity_z_peculiar;
          candidate_particles.mass_code[destination] = record.mass_code;
          result.inbound_particle_ids.push_back(record.particle_id);
          assembly.encoded_record.clear();
          assembly.active = false;
          assembly.total_bytes = 0U;
          assembly.received_bytes = 0U;
        }
      }
      result.stats.record_assembly_capacity_bytes = std::max(result.stats.record_assembly_capacity_bytes,
                                                               assemblyCapacityBytes());
      result.stats.communication_high_water_bytes = std::max(
          result.stats.communication_high_water_bytes,
          core::checkedMemoryBytesAdd(
              core::checkedMemoryBytesAdd(result.stats.packet_send_capacity_bytes,
                                          result.stats.packet_receive_capacity_bytes,
                                          "compact DMO packet capacities"),
              result.stats.record_assembly_capacity_bytes,
              "compact DMO communication high water"));
    } catch (...) { local_failure = std::current_exception(); }
    mpi_context.rethrowCollectivePreparationFailure(local_failure,
        "compact DMO migration bounded packet receive");
  }
  for (const auto& assembly : inbound_assembly) {
    if (assembly.active || !assembly.encoded_record.empty()) {
      throw std::runtime_error("compact DMO migration ended with incomplete inbound record");
    }
  }
  if (result.inbound_particle_ids.size() != total_inbound) {
    throw std::runtime_error("compact DMO migration received count differs from admitted count");
  }
  return result;
#else
  throw std::runtime_error("compact DMO migration exchange requires an MPI-enabled build");
#endif
}

MigrationExchangeResult<core::ParticleMigrationRecord> exchangeRuntimeParticleMigrationRecords(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<core::ParticleMigrationRecord>>& records_by_rank) {
  return exchangeRuntimeMigrationPayloads<core::ParticleMigrationRecord>(
      mpi_context,
      records_by_rank,
      "particle migration bounded packet exchange",
      [](const core::ParticleMigrationRecord& record) {
        return migration_wire::encodeParticleMigrationRecord(record);
      },
      [](std::span<const std::uint8_t> bytes) {
        return migration_wire::decodeParticleMigrationRecord(bytes);
      });
}

MigrationExchangeResult<core::AmrPatchMigrationRecord> exchangeRuntimeAmrPatchMigrationRecords(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<core::AmrPatchMigrationRecord>>& records_by_rank) {
  return exchangeRuntimeMigrationPayloads<core::AmrPatchMigrationRecord>(
      mpi_context,
      records_by_rank,
      "AMR patch migration bounded packet exchange",
      [](const core::AmrPatchMigrationRecord& record) {
        return migration_wire::encodeAmrPatchMigrationRecord(record);
      },
      [](std::span<const std::uint8_t> bytes) {
        return migration_wire::decodeAmrPatchMigrationRecord(bytes);
      });
}


enum class MigrationTransactionRepresentation : std::uint8_t {
  kCompactHomogeneousDmo = 0,
  kGenericFullPhysics = 1,
};

struct MigrationAdmissionPlan {
  MigrationTransactionRepresentation representation = MigrationTransactionRepresentation::kGenericFullPhysics;
  std::uint64_t outbound_particle_count = 0U;
  std::uint64_t inbound_particle_count = 0U;
  std::uint64_t outbound_particle_record_capacity_count = 0U;
  std::uint64_t inbound_particle_record_capacity_count = 0U;
  std::uint64_t outbound_patch_count = 0U;
  std::uint64_t inbound_patch_count = 0U;
  std::uint64_t outbound_wire_bytes = 0U;
  std::uint64_t inbound_wire_bytes = 0U;
  std::uint64_t outbound_dynamic_heap_bytes = 0U;
  std::uint64_t inbound_dynamic_heap_bytes = 0U;
  std::uint64_t packet_staging_bytes = 0U;
  std::uint64_t record_capacity_bytes = 0U;
  std::uint64_t scheduler_remap_bytes = 0U;
  std::uint64_t index_map_bytes = 0U;
  std::uint64_t commit_coexistence_bytes = 0U;
  std::uint64_t metadata_materialization_bytes = 0U;
  std::uint64_t candidate_state_capacity_bytes = 0U;
  std::uint64_t selection_capacity_bytes = 0U;
  std::uint64_t record_assembly_bound_bytes = 0U;
  std::uint64_t exact_ownership_workspace_bound_bytes = 0U;
  std::uint64_t uniform_scheduler_capacity_bytes = 0U;
  bool candidate_requires_sfc_lane = false;
  std::uint64_t requested_extra_bytes = 0U;
  std::vector<std::uint64_t> outbound_particle_count_by_rank;
  std::vector<std::uint64_t> inbound_particle_count_by_rank;
  std::vector<std::uint64_t> outbound_particle_record_capacity_count_by_rank;
  std::vector<std::uint64_t> outbound_patch_count_by_rank;
};

[[nodiscard]] std::vector<std::uint64_t> exchangeMigrationControlU64(
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> send_values,
    std::string_view phase_name) {
  const int world_size = mpi_context.worldSize();
  if (world_size <= 0 || send_values.size() != static_cast<std::size_t>(world_size)) {
    throw std::invalid_argument("migration control exchange rank extent mismatch");
  }
  if (world_size == 1) return std::vector<std::uint64_t>(send_values.begin(), send_values.end());
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("migration control exchange requires MPI when world size exceeds one");
  }
#if COSMOSIM_ENABLE_MPI
  std::vector<std::uint64_t> recv_values(static_cast<std::size_t>(world_size), 0U);
  std::exception_ptr local_failure;
  try {
    if (MPI_Alltoall(
            send_values.data(), 1, MPI_UINT64_T,
            recv_values.data(), 1, MPI_UINT64_T,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(std::string(phase_name) + " MPI_Alltoall failed");
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(local_failure, phase_name);
  return recv_values;
#else
  throw std::runtime_error("migration control exchange requires an MPI-enabled build");
#endif
}

[[nodiscard]] std::uint64_t sumU64(
    std::span<const std::uint64_t> values,
    std::string_view context) {
  std::uint64_t total = 0U;
  for (const std::uint64_t value : values) {
    total = core::checkedMemoryBytesAdd(total, value, context);
  }
  return total;
}

[[nodiscard]] std::uint32_t findLocalPatchIndexById(
    const core::SimulationState& state,
    std::uint64_t patch_id) {
  const auto it = std::find(state.patches.patch_id.begin(), state.patches.patch_id.end(), patch_id);
  if (it == state.patches.patch_id.end()) {
    throw std::runtime_error("runtime rebalance references an AMR patch ID absent from local state");
  }
  return core::checkedIntegralNarrow<std::uint32_t>(
      static_cast<std::size_t>(std::distance(state.patches.patch_id.begin(), it)),
      "runtime AMR patch local index");
}


[[nodiscard]] bool compactDmoMigrationEligible(
    const core::SimulationState& state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const parallel::RuntimeRebalancePlan& rebalance) {
  return state.hasHomogeneousDmoMetadata() &&
      scheduler.representation() == core::SchedulerRepresentation::kUniformRungZero &&
      scheduler.maxBin() == 0U &&
      state.cells.size() == 0U && state.patches.size() == 0U &&
      state.star_particles.size() == 0U && state.black_holes.size() == 0U && state.tracers.size() == 0U &&
      state.sidecars.size() == 0U &&
      state.particle_sidecar.gravity_softening_comoving.empty() &&
      state.particle_sidecar.has_gravity_softening_override.empty() &&
      rebalance.amr_patch_ownership_updates.empty();
}

[[nodiscard]] bool compactDmoMetadataConsensus(
    const core::SimulationState& state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const parallel::MpiContext& mpi_context) {
  if (!state.hasHomogeneousDmoMetadata() ||
      scheduler.representation() != core::SchedulerRepresentation::kUniformRungZero) return false;
  const auto agrees_u64 = [&](std::uint64_t value) {
    return mpi_context.allreduceMinUint64(value) == mpi_context.allreduceMaxUint64(value);
  };
  return agrees_u64(state.homogeneousDmoParticleFlags()) &&
      agrees_u64(std::bit_cast<std::uint64_t>(state.homogeneousDmoLastDriftTimeCode())) &&
      agrees_u64(std::bit_cast<std::uint64_t>(state.homogeneousDmoLastDriftScaleFactor())) &&
      agrees_u64(scheduler.currentTick()) && agrees_u64(scheduler.maxBin());
}

[[nodiscard]] MigrationAdmissionPlan buildCompactDmoMigrationAdmissionPlan(
    const core::SimulationState& state,
    const parallel::RuntimeRebalancePlan& rebalance,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    bool exact_ownership_audit_enabled) {
  const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
  MigrationAdmissionPlan plan;
  plan.representation = MigrationTransactionRepresentation::kCompactHomogeneousDmo;
  plan.outbound_particle_count_by_rank.assign(rank_count, 0U);
  plan.outbound_particle_record_capacity_count_by_rank.assign(rank_count, 0U);
  plan.outbound_patch_count_by_rank.assign(rank_count, 0U);
  for (const auto& intent : rebalance.particle_migrations) {
    if (intent.old_owner_rank != world_rank || intent.new_owner_rank == world_rank) continue;
    if (intent.new_owner_rank < 0 || intent.new_owner_rank >= mpi_context.worldSize() ||
        intent.item_index >= state.particles.size() ||
        state.particle_sidecar.particle_id[intent.item_index] != intent.particle_id) {
      throw std::runtime_error("compact DMO migration intent failed authoritative local-row validation");
    }
    auto& count = plan.outbound_particle_count_by_rank[static_cast<std::size_t>(intent.new_owner_rank)];
    count = core::checkedMemoryBytesAdd(count, 1U, "compact DMO outbound count");
  }
  plan.inbound_particle_count_by_rank = exchangeMigrationControlU64(
      mpi_context, plan.outbound_particle_count_by_rank,
      "compact DMO particle-count admission exchange");
  plan.outbound_particle_count = sumU64(plan.outbound_particle_count_by_rank,
                                         "compact DMO outbound particle count");
  plan.inbound_particle_count = sumU64(plan.inbound_particle_count_by_rank,
                                        "compact DMO inbound particle count");
  const std::uint64_t local_rows = core::checkedIntegralNarrow<std::uint64_t>(
      state.particles.size(), "compact DMO local particle rows");
  if (plan.outbound_particle_count > local_rows) {
    throw std::runtime_error("compact DMO outbound migration count exceeds local particle population");
  }
  const std::uint64_t final_rows = core::checkedMemoryBytesAdd(
      local_rows - plan.outbound_particle_count, plan.inbound_particle_count,
      "compact DMO final particle rows");

  const std::uint64_t local_nonuniform_sfc = state.homogeneousDmoSfcKeyIsUniform() ? 0U : 1U;
  const bool any_nonuniform_sfc = mpi_context.allreduceSumUint64(local_nonuniform_sfc) != 0U;
  const std::uint64_t local_uniform_sfc = state.homogeneousDmoSfcKeyIsUniform()
      ? state.homogeneousDmoUniformSfcKey() : 0U;
  const bool uniform_key_disagreement = !any_nonuniform_sfc &&
      mpi_context.allreduceMinUint64(local_uniform_sfc) != mpi_context.allreduceMaxUint64(local_uniform_sfc);
  plan.candidate_requires_sfc_lane = any_nonuniform_sfc || uniform_key_disagreement;

  const std::uint64_t bytes_per_candidate = sizeof(double) * 7U + sizeof(std::uint64_t) +
      (plan.candidate_requires_sfc_lane ? sizeof(std::uint64_t) : 0U);
  plan.candidate_state_capacity_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(final_rows, "compact DMO candidate row count"),
          core::checkedIntegralNarrow<std::size_t>(bytes_per_candidate, "compact DMO candidate row bytes"),
          "compact DMO candidate state capacity"),
      "compact DMO candidate state byte width");
  plan.commit_coexistence_bytes = plan.candidate_state_capacity_bytes;
  plan.selection_capacity_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeAdd(
          core::checkedSizeMultiply(
              core::checkedIntegralNarrow<std::size_t>(plan.outbound_particle_count,
                  "compact DMO selection count"),
              sizeof(std::uint32_t) * 2U, "compact DMO selection capacities"),
          core::checkedSizeMultiply(rank_count, sizeof(std::vector<std::uint32_t>),
                                    "compact DMO peer selection headers"),
          "compact DMO selection total capacity"),
      "compact DMO selection byte width");
  plan.uniform_scheduler_capacity_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(final_rows,
              "compact DMO scheduler final rows"),
          sizeof(std::uint32_t), "compact DMO uniform scheduler identity lane"),
      "compact DMO uniform scheduler byte width");
  const std::uint64_t round_limit = core::checkedIntegralNarrow<std::uint64_t>(
      parallel::mpiTransportRoundLimitBytes(), "compact DMO MPI round limit");
  plan.record_assembly_bound_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(rank_count, migration_wire::dmoParticleMigrationWireBytes(),
                                "compact DMO record assembly bound"),
      "compact DMO record assembly bound width");
  plan.packet_staging_bytes = core::checkedMemoryBytesAdd(
      core::checkedMemoryBytesAdd(round_limit, round_limit, "compact DMO MPI internal round buffers"),
      core::checkedMemoryBytesAdd(round_limit, round_limit, "compact DMO packet send/receive vectors"),
      "compact DMO packet staging");
  plan.packet_staging_bytes = core::checkedMemoryBytesAdd(
      plan.packet_staging_bytes, plan.record_assembly_bound_bytes,
      "compact DMO packet plus assembly staging");
  plan.outbound_wire_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(plan.outbound_particle_count,
              "compact DMO outbound wire count"),
          migration_wire::dmoParticleMigrationWireBytes() + migration_wire::k_fragment_header_bytes,
          "compact DMO outbound wire upper bound"),
      "compact DMO outbound wire byte width");
  plan.inbound_wire_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(plan.inbound_particle_count,
              "compact DMO inbound wire count"),
          migration_wire::dmoParticleMigrationWireBytes() + migration_wire::k_fragment_header_bytes,
          "compact DMO inbound wire upper bound"),
      "compact DMO inbound wire byte width");
  const std::uint64_t inbound_duplicate_staging = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(plan.inbound_particle_count,
              "compact DMO inbound duplicate count"),
          sizeof(std::uint64_t), "compact DMO inbound duplicate staging"),
      "compact DMO inbound duplicate byte width");
  plan.exact_ownership_workspace_bound_bytes = exact_ownership_audit_enabled
      ? parallel::k_exact_ownership_validation_workspace_limit_bytes : 0U;
  plan.requested_extra_bytes = plan.candidate_state_capacity_bytes;
  for (const auto [bytes, label] : std::array{
           std::pair{plan.selection_capacity_bytes, std::string_view{"compact DMO outbound selection"}},
           std::pair{plan.packet_staging_bytes, std::string_view{"compact DMO bounded packets"}},
           std::pair{inbound_duplicate_staging, std::string_view{"compact DMO inbound duplicate staging"}},
           std::pair{plan.uniform_scheduler_capacity_bytes, std::string_view{"compact DMO uniform scheduler candidate"}},
           std::pair{plan.exact_ownership_workspace_bound_bytes, std::string_view{"compact DMO exact ownership workspace"}},
       }) {
    plan.requested_extra_bytes = core::checkedMemoryBytesAdd(plan.requested_extra_bytes, bytes, label);
  }
  return plan;
}

[[nodiscard]] MigrationAdmissionPlan buildMigrationAdmissionPlan(
    const core::SimulationState& state,
    const parallel::RuntimeRebalancePlan& rebalance,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    std::uint64_t transaction_baseline_before) {
  const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
  MigrationAdmissionPlan plan;
  plan.representation = MigrationTransactionRepresentation::kGenericFullPhysics;
  plan.outbound_particle_count_by_rank.assign(rank_count, 0U);
  plan.outbound_particle_record_capacity_count_by_rank.assign(rank_count, 0U);
  plan.outbound_patch_count_by_rank.assign(rank_count, 0U);
  std::vector<std::uint64_t> send_wire(rank_count, 0U);
  std::vector<std::uint64_t> send_dynamic_heap(rank_count, 0U);
  std::vector<std::uint64_t> send_max_record_wire(rank_count, 0U);
  std::vector<std::uint64_t> send_state_bytes(rank_count, 0U);
  std::vector<std::uint64_t> send_patch_particle_upper_bound(rank_count, 0U);

  std::exception_ptr local_sizing_failure;
  try {
    std::uint64_t max_particle_wire = 0U;
    std::uint64_t max_particle_dynamic_heap = 0U;
    std::uint64_t max_particle_state_bytes = 0U;
    for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const auto local_index = core::checkedIntegralNarrow<std::uint32_t>(
        particle_index, "migration admission particle scan index");
    max_particle_wire = std::max(
        max_particle_wire,
        core::checkedIntegralNarrow<std::uint64_t>(
            migration_wire::estimateParticleMigrationWireUpperBoundBytes(state, local_index),
            "particle migration maximum wire upper bound"));
    max_particle_dynamic_heap = std::max(
        max_particle_dynamic_heap,
        core::checkedIntegralNarrow<std::uint64_t>(
            migration_wire::estimateParticleMigrationDynamicHeapUpperBoundBytes(state, local_index),
            "particle migration maximum dynamic heap upper bound"));
    max_particle_state_bytes = std::max(
        max_particle_state_bytes,
        estimateParticleMemoryBytesForDecomposition(
            state, state.particleSpeciesTag(particle_index)));
    }

    for (const parallel::ParticleMigrationIntent& intent : rebalance.particle_migrations) {
    if (intent.old_owner_rank != world_rank) continue;
    if (intent.new_owner_rank < 0 || intent.new_owner_rank >= mpi_context.worldSize()) {
      throw std::invalid_argument("migration admission encountered out-of-range particle target rank");
    }
    if (intent.item_index >= state.particles.size()) {
      throw std::runtime_error("migration admission particle item index is not a local particle row");
    }
    const std::size_t target = static_cast<std::size_t>(intent.new_owner_rank);
    const auto local_index = core::checkedIntegralNarrow<std::uint32_t>(
        intent.item_index, "migration admission particle local index");
    const std::uint64_t wire = core::checkedIntegralNarrow<std::uint64_t>(
        migration_wire::estimateParticleMigrationWireUpperBoundBytes(state, local_index),
        "particle migration wire upper bound");
    const std::uint64_t dynamic_heap = core::checkedIntegralNarrow<std::uint64_t>(
        migration_wire::estimateParticleMigrationDynamicHeapUpperBoundBytes(state, local_index),
        "particle migration dynamic heap upper bound");
    const std::uint64_t state_bytes = estimateParticleMemoryBytesForDecomposition(
        state, state.particleSpeciesTag(local_index));
    plan.outbound_particle_count_by_rank[target] = core::checkedMemoryBytesAdd(
        plan.outbound_particle_count_by_rank[target], 1U,
        "particle migration outbound count");
    send_wire[target] = core::checkedMemoryBytesAdd(
        send_wire[target], wire, "particle migration outbound wire bytes");
    send_dynamic_heap[target] = core::checkedMemoryBytesAdd(
        send_dynamic_heap[target], dynamic_heap,
        "particle migration outbound dynamic heap bytes");
    send_max_record_wire[target] = std::max(send_max_record_wire[target], wire);
    send_state_bytes[target] = core::checkedMemoryBytesAdd(
        send_state_bytes[target], state_bytes,
        "particle migration outbound state bytes");
    }

    for (const parallel::AmrPatchOwnershipUpdate& update : rebalance.amr_patch_ownership_updates) {
    if (update.old_owner_rank != world_rank) continue;
    if (update.new_owner_rank < 0 || update.new_owner_rank >= mpi_context.worldSize()) {
      throw std::invalid_argument("migration admission encountered out-of-range AMR patch target rank");
    }
    const std::uint32_t patch_index = findLocalPatchIndexById(state, update.patch_id);
    const std::size_t target = static_cast<std::size_t>(update.new_owner_rank);
    const std::uint64_t patch_cell_count = core::checkedIntegralNarrow<std::uint64_t>(
        state.patches.cell_count[patch_index], "AMR migration cell count");
    const std::uint64_t wire = core::checkedIntegralNarrow<std::uint64_t>(
        migration_wire::estimateAmrPatchMigrationWireUpperBoundBytes(state, patch_index),
        "AMR migration wire upper bound");
    const std::uint64_t dynamic_heap = core::checkedIntegralNarrow<std::uint64_t>(
        migration_wire::estimateAmrPatchMigrationDynamicHeapUpperBoundBytes(state, patch_index),
        "AMR migration dynamic heap upper bound");
    const std::uint64_t state_bytes = core::checkedIntegralNarrow<std::uint64_t>(
        core::checkedSizeMultiply(
            state.patches.cell_count[patch_index],
            sizeof(double) * 8U + sizeof(std::uint32_t) * 2U,
            "AMR migration state byte upper bound"),
        "AMR migration state byte width");
    plan.outbound_patch_count_by_rank[target] = core::checkedMemoryBytesAdd(
        plan.outbound_patch_count_by_rank[target], 1U,
        "AMR migration outbound patch count");
    send_wire[target] = core::checkedMemoryBytesAdd(
        send_wire[target], wire, "AMR migration outbound wire bytes");
    send_dynamic_heap[target] = core::checkedMemoryBytesAdd(
        send_dynamic_heap[target], dynamic_heap,
        "AMR migration outbound dynamic heap bytes");
    send_max_record_wire[target] = std::max(send_max_record_wire[target], wire);
    send_state_bytes[target] = core::checkedMemoryBytesAdd(
        send_state_bytes[target], state_bytes,
        "AMR migration outbound state bytes");

    // Patch reassignment can induce parent-gas-particle migration in addition
    // to explicit particle intents.  Account for at most one particle per
    // patch cell here without materializing the parent-index map before
    // governor admission.  Double-counting an explicit overlapping particle
    // intent is intentionally conservative; under-counting this path is not.
    send_patch_particle_upper_bound[target] = core::checkedMemoryBytesAdd(
        send_patch_particle_upper_bound[target], patch_cell_count,
        "AMR-induced particle migration upper-bound count");
    send_wire[target] = core::checkedMemoryBytesAdd(
        send_wire[target],
        core::checkedIntegralNarrow<std::uint64_t>(
            core::checkedSizeMultiply(
                core::checkedIntegralNarrow<std::size_t>(
                    patch_cell_count, "AMR-induced particle count width"),
                core::checkedIntegralNarrow<std::size_t>(
                    max_particle_wire, "maximum particle wire width"),
                "AMR-induced particle wire upper bound"),
            "AMR-induced particle wire byte width"),
        "AMR plus particle outbound wire bytes");
    send_dynamic_heap[target] = core::checkedMemoryBytesAdd(
        send_dynamic_heap[target],
        core::checkedIntegralNarrow<std::uint64_t>(
            core::checkedSizeMultiply(
                core::checkedIntegralNarrow<std::size_t>(
                    patch_cell_count, "AMR-induced dynamic particle count width"),
                core::checkedIntegralNarrow<std::size_t>(
                    max_particle_dynamic_heap, "maximum particle dynamic heap width"),
                "AMR-induced particle dynamic heap upper bound"),
            "AMR-induced particle dynamic heap byte width"),
        "AMR plus particle outbound dynamic heap bytes");
    send_state_bytes[target] = core::checkedMemoryBytesAdd(
        send_state_bytes[target],
        core::checkedIntegralNarrow<std::uint64_t>(
            core::checkedSizeMultiply(
                core::checkedIntegralNarrow<std::size_t>(
                    patch_cell_count, "AMR-induced state particle count width"),
                core::checkedIntegralNarrow<std::size_t>(
                    max_particle_state_bytes, "maximum particle state byte width"),
                "AMR-induced particle state upper bound"),
            "AMR-induced particle state byte width"),
        "AMR plus particle outbound state bytes");
    send_max_record_wire[target] = std::max(send_max_record_wire[target], max_particle_wire);
    }
  } catch (...) {
    local_sizing_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_sizing_failure, "migration admission local physical sizing");

  const std::vector<std::uint64_t> recv_particle_counts = exchangeMigrationControlU64(
      mpi_context, plan.outbound_particle_count_by_rank, "particle migration count admission exchange");
  plan.inbound_particle_count_by_rank = recv_particle_counts;
  const std::vector<std::uint64_t> recv_patch_counts = exchangeMigrationControlU64(
      mpi_context, plan.outbound_patch_count_by_rank, "AMR migration count admission exchange");
  const std::vector<std::uint64_t> recv_wire = exchangeMigrationControlU64(
      mpi_context, send_wire, "migration wire-byte admission exchange");
  const std::vector<std::uint64_t> recv_dynamic_heap = exchangeMigrationControlU64(
      mpi_context, send_dynamic_heap, "migration dynamic-heap admission exchange");
  const std::vector<std::uint64_t> recv_max_record_wire = exchangeMigrationControlU64(
      mpi_context, send_max_record_wire, "migration max-record admission exchange");
  const std::vector<std::uint64_t> recv_state_bytes = exchangeMigrationControlU64(
      mpi_context, send_state_bytes, "migration state-byte admission exchange");
  const std::vector<std::uint64_t> recv_patch_particle_upper_bound = exchangeMigrationControlU64(
      mpi_context, send_patch_particle_upper_bound,
      "AMR-induced particle-count admission exchange");

  plan.outbound_particle_count = sumU64(plan.outbound_particle_count_by_rank, "outbound particle migration count");
  plan.inbound_particle_count = sumU64(recv_particle_counts, "inbound particle migration count");
  plan.outbound_patch_count = sumU64(plan.outbound_patch_count_by_rank, "outbound AMR patch migration count");
  plan.inbound_patch_count = sumU64(recv_patch_counts, "inbound AMR patch migration count");
  plan.outbound_particle_record_capacity_count = core::checkedMemoryBytesAdd(
      plan.outbound_particle_count,
      sumU64(send_patch_particle_upper_bound, "outbound AMR-induced particle capacity count"),
      "outbound particle record capacity count");
  for (std::size_t rank = 0; rank < rank_count; ++rank) {
    plan.outbound_particle_record_capacity_count_by_rank[rank] = core::checkedMemoryBytesAdd(
        plan.outbound_particle_count_by_rank[rank],
        send_patch_particle_upper_bound[rank],
        "outbound per-rank particle record capacity count");
  }
  plan.inbound_particle_record_capacity_count = core::checkedMemoryBytesAdd(
      plan.inbound_particle_count,
      sumU64(recv_patch_particle_upper_bound, "inbound AMR-induced particle capacity count"),
      "inbound particle record capacity count");
  plan.outbound_wire_bytes = sumU64(send_wire, "outbound migration wire bytes");
  plan.inbound_wire_bytes = sumU64(recv_wire, "inbound migration wire bytes");
  plan.outbound_dynamic_heap_bytes = sumU64(send_dynamic_heap, "outbound migration dynamic heap bytes");
  plan.inbound_dynamic_heap_bytes = sumU64(recv_dynamic_heap, "inbound migration dynamic heap bytes");

  const std::uint64_t particle_records = core::checkedMemoryBytesAdd(
      plan.outbound_particle_record_capacity_count,
      plan.inbound_particle_record_capacity_count,
      "migration particle record count");
  const std::uint64_t patch_records = core::checkedMemoryBytesAdd(
      plan.outbound_patch_count, plan.inbound_patch_count,
      "migration patch record count");
  plan.record_capacity_bytes = core::checkedMemoryBytesAdd(
      core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(particle_records, "migration particle record count width"),
          sizeof(core::ParticleMigrationRecord),
          "migration particle record capacity"),
      core::checkedSizeMultiply(
          core::checkedIntegralNarrow<std::size_t>(patch_records, "migration patch record count width"),
          sizeof(core::AmrPatchMigrationRecord),
          "migration patch record capacity"),
      "migration native record capacities");
  plan.record_capacity_bytes = core::checkedMemoryBytesAdd(
      plan.record_capacity_bytes,
      core::checkedMemoryBytesAdd(
          plan.outbound_dynamic_heap_bytes,
          plan.inbound_dynamic_heap_bytes,
          "migration dynamic record heap total"),
      "migration record plus dynamic capacity");

  const std::uint64_t expected_particle_rows = core::checkedMemoryBytesAdd(
      core::checkedIntegralNarrow<std::uint64_t>(state.particles.size(), "migration current particle rows"),
      plan.inbound_particle_count,
      "migration scheduler conservative destination rows");
  plan.scheduler_remap_bytes = core::checkedMemoryBytesAdd(
      core::checkedIntegralNarrow<std::uint64_t>(
          core::checkedSizeMultiply(
              state.particles.size(), sizeof(std::uint32_t),
              "migration preserved-index capacity"),
          "migration preserved-index bytes"),
      core::checkedIntegralNarrow<std::uint64_t>(
          core::checkedSizeMultiply(
              core::checkedIntegralNarrow<std::size_t>(expected_particle_rows, "migration scheduler rows"),
              sizeof(core::TimeBinSchedulerIdentityRecord),
              "migration scheduler identity capacity"),
          "migration scheduler identity bytes"),
      "migration scheduler remap bytes");
  plan.scheduler_remap_bytes = core::checkedMemoryBytesAdd(
      plan.scheduler_remap_bytes,
      core::checkedIntegralNarrow<std::uint64_t>(
          core::checkedSizeMultiply(
              core::checkedIntegralNarrow<std::size_t>(expected_particle_rows, "migration destination ID rows"),
              sizeof(std::uint64_t),
              "migration destination particle ID capacity"),
          "migration destination particle ID bytes"),
      "migration scheduler/destination remap bytes");

  constexpr std::size_t k_unordered_node_overhead = sizeof(void*) * 3U;
  const std::size_t local_index_node_bytes =
      sizeof(std::pair<const std::uint64_t, std::uint32_t>) + k_unordered_node_overhead;
  const std::size_t target_map_node_bytes =
      sizeof(std::pair<const std::uint32_t, int>) + k_unordered_node_overhead;
  plan.index_map_bytes = core::checkedMemoryBytesAdd(
      core::checkedIntegralNarrow<std::uint64_t>(
          core::checkedSizeMultiply(state.particles.size(), local_index_node_bytes,
                                    "migration local identity map capacity"),
          "migration local identity map bytes"),
      core::checkedIntegralNarrow<std::uint64_t>(
          core::checkedSizeMultiply(
              core::checkedIntegralNarrow<std::size_t>(plan.outbound_particle_count, "migration outbound map rows"),
              target_map_node_bytes,
              "migration outbound target map capacity"),
          "migration outbound target map bytes"),
      "migration index-map bytes");

  const std::uint64_t inbound_state_bytes = sumU64(recv_state_bytes, "migration inbound state bytes");
  plan.commit_coexistence_bytes = core::checkedMemoryBytesAdd(
      transaction_baseline_before,
      inbound_state_bytes,
      "migration old/new commit coexistence");

  const std::uint64_t max_record_assembly = core::checkedMemoryBytesAdd(
      sumU64(send_max_record_wire, "migration outbound max-record staging"),
      sumU64(recv_max_record_wire, "migration inbound max-record staging"),
      "migration record fragment assembly staging");
  const std::uint64_t round_limit = core::checkedIntegralNarrow<std::uint64_t>(
      parallel::mpiTransportRoundLimitBytes(), "migration MPI transport round limit");
  plan.packet_staging_bytes = core::checkedMemoryBytesAdd(
      core::checkedMemoryBytesAdd(round_limit, round_limit, "migration MPI internal round buffers"),
      core::checkedMemoryBytesAdd(round_limit, round_limit, "migration packet send/receive vectors"),
      "migration bounded packet staging");
  plan.packet_staging_bytes = core::checkedMemoryBytesAdd(
      plan.packet_staging_bytes, max_record_assembly,
      "migration packet plus single-record assembly staging");

  if (state.hasHomogeneousDmoMetadata()) {
    const std::uint64_t materialization_bytes_per_particle =
        sizeof(std::uint8_t) + sizeof(std::uint32_t) * 3U + sizeof(double) * 2U +
        sizeof(std::uint32_t) * 2U +
        (state.particle_sidecar.sfc_key.empty() ? sizeof(std::uint64_t) : 0U);
    plan.metadata_materialization_bytes = core::checkedIntegralNarrow<std::uint64_t>(
        core::checkedSizeMultiply(
            state.particles.size(),
            static_cast<std::size_t>(materialization_bytes_per_particle),
            "migration compact-DMO metadata materialization"),
        "migration compact-DMO metadata materialization bytes");
  }

  plan.requested_extra_bytes = core::checkedMemoryBytesAdd(
      plan.commit_coexistence_bytes, plan.metadata_materialization_bytes,
      "migration commit/materialization coexistence");
  for (const auto [bytes, label] : std::array{
           std::pair{plan.record_capacity_bytes, std::string_view{"migration record capacities"}},
           std::pair{plan.scheduler_remap_bytes, std::string_view{"migration scheduler remap"}},
           std::pair{plan.index_map_bytes, std::string_view{"migration index maps"}},
           std::pair{plan.packet_staging_bytes, std::string_view{"migration bounded packets"}},
       }) {
    plan.requested_extra_bytes = core::checkedMemoryBytesAdd(
        plan.requested_extra_bytes, bytes, label);
  }
  return plan;
}

[[nodiscard]] parallel::LocalOwnershipIdentitySummary reduceLocalParticleIdentitySummary(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context) {
  const parallel::LocalOwnershipIdentitySummary local =
      parallel::summarizeLocalOwnedParticleIds(
          state.particle_sidecar.particle_id, state.validateUniqueParticleIds());
  const std::uint64_t unique_rank_count =
      mpi_context.allreduceSumUint64(local.local_particle_ids_unique ? 1ULL : 0ULL);
  return parallel::LocalOwnershipIdentitySummary{
      .local_owned_count = mpi_context.allreduceSumUint64(local.local_owned_count),
      .local_particle_id_sum = mpi_context.allreduceSumUint64(local.local_particle_id_sum),
      .local_particle_id_square_sum = mpi_context.allreduceSumUint64(local.local_particle_id_square_sum),
      .local_particle_id_xor = mpi_context.allreduceXorUint64(local.local_particle_id_xor),
      .local_particle_ids_unique =
          unique_rank_count == static_cast<std::uint64_t>(mpi_context.worldSize()),
  };
}

void requireGlobalOwnedParticlePartitionIdentity(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> expected_global_particle_ids,
    std::string_view caller) {
  std::vector<std::uint64_t> generic_local_owned_ids;
  std::span<const std::uint64_t> local_owned_ids = state.particle_sidecar.particle_id;
  if (!state.hasHomogeneousDmoMetadata()) {
    generic_local_owned_ids.reserve(state.particles.size());
    const std::uint32_t world_rank = static_cast<std::uint32_t>(mpi_context.worldRank());
    for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
      if (state.particleOwningRank(particle_index) == world_rank) {
        generic_local_owned_ids.push_back(state.particle_sidecar.particle_id[particle_index]);
      }
    }
    local_owned_ids = generic_local_owned_ids;
  }
  const parallel::ExactOwnershipPartitionReport report =
      parallel::validateExactGlobalOwnershipPartition(mpi_context, local_owned_ids, expected_global_particle_ids);
  if (!report.valid()) {
    throw std::runtime_error(std::string(caller) +
        ": exact distributed authoritative ownership table has duplicate=" +
        std::to_string(report.duplicate_count) +
        ", missing=" + std::to_string(report.missing_count) +
        ", extra=" + std::to_string(report.extra_count) + " particle IDs");
  }
}

void recordRuntimeRebalanceDecision(
    core::ProfilerSession* profiler,
    const parallel::RuntimeRebalancePlan& rebalance,
    std::uint64_t step_index) {
  if (profiler == nullptr) {
    return;
  }
  const auto join_u64 = [](std::span<const std::uint64_t> values) {
    std::ostringstream stream;
    for (std::size_t i = 0; i < values.size(); ++i) {
      if (i != 0U) {
        stream << ',';
      }
      stream << values[i];
    }
    return stream.str();
  };
  profiler->recordEvent(core::RuntimeEvent{
      .event_kind = "parallel.decomposition.runtime_rebalance",
      .severity = rebalance.should_rebalance ? core::RuntimeEventSeverity::kWarning : core::RuntimeEventSeverity::kInfo,
      .subsystem = "parallel.domain_decomposition",
      .step_index = step_index,
      .message = rebalance.should_rebalance
          ? "runtime work-weighted decomposition produced migration intent"
          : "runtime work-weighted decomposition remained below rebalance threshold",
      .payload = {{"reason", rebalance.reason},
                  {"should_rebalance", rebalance.should_rebalance ? "true" : "false"},
                  {"particle_migration_count", std::to_string(rebalance.particle_migrations.size())},
                  {"amr_patch_ownership_update_count", std::to_string(rebalance.amr_patch_ownership_updates.size())},
                  {"migrated_load_fraction", std::to_string(rebalance.migrated_load_fraction)},
                  {"used_distributed_sfc_cuts", rebalance.used_distributed_sfc_cuts ? "true" : "false"},
                  {"used_compact_planner", rebalance.used_compact_planner ? "true" : "false"},
                  {"planner_record_bytes", std::to_string(rebalance.planner_record_bytes)},
                  {"planner_sample_bytes", std::to_string(rebalance.planner_sample_bytes)},
                  {"planner_prefix_bytes", std::to_string(rebalance.planner_prefix_bytes)},
                   {"planner_migration_intent_bytes", std::to_string(rebalance.planner_migration_intent_bytes)},
                   {"planner_other_known_scratch_bytes", std::to_string(rebalance.planner_other_known_scratch_bytes)},
                   {"planner_peak_temporary_bytes", std::to_string(rebalance.planner_peak_temporary_bytes)},
                   {"planner_local_peak_temporary_bytes", std::to_string(rebalance.planner_local_peak_temporary_bytes)},
                   {"planner_bytes_per_entity", std::to_string(rebalance.planner_bytes_per_entity)},
                   {"planner_local_entity_count", std::to_string(rebalance.planner_local_entity_count)},
                   {"planner_bytes_per_local_entity", std::to_string(rebalance.planner_bytes_per_entity)},

                  {"exact_debug_audit_enabled", rebalance.exact_debug_audit_enabled ? "true" : "false"},
                  {"local_entities_considered", std::to_string(rebalance.local_entities_considered)},
                  {"global_entities_considered", std::to_string(rebalance.global_entities_considered)},
                  {"local_entities_moved", std::to_string(rebalance.local_entities_moved)},
                  {"global_entities_moved", std::to_string(rebalance.global_entities_moved)},
                  {"local_bytes_moved", std::to_string(rebalance.local_bytes_moved)},
                  {"global_bytes_moved", std::to_string(rebalance.global_bytes_moved)},
                  {"global_control_bytes", std::to_string(rebalance.global_control_bytes)},
                  {"peak_temporary_bytes", std::to_string(rebalance.peak_temporary_bytes)},
                  {"cut_displacement_fraction", std::to_string(rebalance.cut_displacement_fraction)},
                  {"sfc_cut_count", std::to_string(rebalance.sfc_cut_keys.size())},
                  {"sfc_cut_keys", join_u64(rebalance.sfc_cut_keys)},
                  {"sfc_cut_entity_ids", join_u64(rebalance.sfc_cut_entity_ids)},
                  {"current_weighted_imbalance_ratio", std::to_string(rebalance.current_metrics.weighted_imbalance_ratio)},
                  {"target_weighted_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.weighted_imbalance_ratio)},
                  {"current_memory_imbalance_ratio", std::to_string(rebalance.current_metrics.memory_imbalance_ratio)},
                  {"target_memory_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.memory_imbalance_ratio)}}});
}

void exchangeAndValidateAmrPatchPayloads(
    core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    std::uint64_t step_index,
    core::ProfilerSession* profiler) {
  const std::vector<parallel::AmrPatchPayloadRecord> local_records =
      buildMigrationAmrPatchPayloadRecords(state, world_rank);
  std::unordered_map<std::uint64_t, std::uint32_t> expected_cell_count_by_patch_id;
  expected_cell_count_by_patch_id.reserve(local_records.size());
  for (const parallel::AmrPatchPayloadRecord& record : local_records) {
    const auto [it, inserted] = expected_cell_count_by_patch_id.emplace(record.patch_id, record.cell_count);
    if (!inserted) {
      throw std::runtime_error("AMR owner-local validation detected duplicate authoritative patch payloads");
    }
    if (record.owner_rank != world_rank) {
      throw std::runtime_error("AMR owner-local validation found stale patch owner metadata");
    }
  }
  for (const parallel::AmrPatchPayloadRecord& record : local_records) {
    if (static_cast<std::uint64_t>(record.first_cell) + record.cell_count > state.cells.size()) {
      throw std::runtime_error(
          "AMR owner-local validation patch range exceeds authoritative cell storage for patch_id=" +
          std::to_string(record.patch_id));
    }
    const std::uint64_t geometry_count = static_cast<std::uint64_t>(record.cell_dim_x) *
        record.cell_dim_y * record.cell_dim_z;
    if (geometry_count != record.cell_count) {
      throw std::runtime_error(
          "AMR owner-local validation patch geometry/cell-count mismatch for patch_id=" +
          std::to_string(record.patch_id));
    }
  }

  parallel::DirectedAmrExchangeDiagnostics diagnostics;
  if (mpi_context.isEnabled() && mpi_context.worldSize() > 1) {
    std::set<std::pair<std::uint64_t, std::uint32_t>> remote_cell_keys;
    const parallel::DirectedAmrPatchPayloadExchange directed = parallel::executeBlockingDirectedAmrPatchPayloadExchange(
        mpi_context,
        local_records,
        [&state, world_rank](
            std::span<const parallel::AmrPatchBoundaryCellRequest> requests,
            std::uint64_t first_record,
            std::size_t max_records,
            std::vector<parallel::AmrPatchCellPayloadRecord>& output) {
          fillMigrationAmrPatchBoundaryCellPayloadChunk(
              state, world_rank, requests, first_record, max_records, output);
        },
        [&remote_cell_keys](
            int,
            std::span<const parallel::AmrPatchCellPayloadRecord> records) {
          for (const auto& record : records) {
            if (!remote_cell_keys.emplace(record.patch_id, record.local_cell_offset).second) {
              throw std::runtime_error(
                  "directed AMR validation received duplicate remote boundary-cell payload");
            }
          }
        },
        {},
        0U,
        step_index);
    diagnostics = directed.diagnostics;

    std::unordered_map<std::uint64_t, parallel::AmrPatchPayloadRecord> remote_patch_by_id;
    remote_patch_by_id.reserve(directed.patch_payloads_received.size());
    for (const parallel::AmrPatchPayloadRecord& record : directed.patch_payloads_received) {
      if (record.owner_rank == world_rank) {
        throw std::runtime_error("directed AMR validation received local authoritative patch as a remote ghost");
      }
      const auto [it, inserted] = remote_patch_by_id.emplace(record.patch_id, record);
      if (!inserted) {
        throw std::runtime_error("directed AMR validation received duplicate remote patch metadata");
      }
    }
  }

  if (profiler != nullptr) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "amr.directed_patch_payload_exchange",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "amr.patch_exchange",
        .step_index = step_index,
        .message = "validated owner-local AMR state and interface-scoped directed AMR payload coverage",
        .payload = {{"local_patch_payloads", std::to_string(local_records.size())},
                    {"local_patch_cell_payloads", "lazy_interface_only"},
                    {"candidate_peer_count", std::to_string(diagnostics.candidate_peer_count)},
                    {"neighbor_peer_count", std::to_string(diagnostics.neighbor_peer_count)},
                    {"directed_patch_descriptor_records_sent", std::to_string(diagnostics.directed_patch_descriptor_records_sent)},
                    {"directed_patch_descriptor_records_received", std::to_string(diagnostics.directed_patch_descriptor_records_received)},
                    {"directed_patch_cell_records_sent", std::to_string(diagnostics.directed_patch_cell_records_sent)},
                    {"directed_patch_cell_records_received", std::to_string(diagnostics.directed_patch_cell_records_received)},
                    {"control_plane_bytes", std::to_string(diagnostics.control_plane_bytes)},
                    {"patch_descriptor_bytes", std::to_string(diagnostics.patch_descriptor_bytes)},
                    {"patch_cell_payload_bytes", std::to_string(diagnostics.patch_cell_payload_bytes)},
                    {"patch_cell_send_capacity_high_water_bytes", std::to_string(diagnostics.patch_cell_send_capacity_high_water_bytes)},
                    {"patch_cell_receive_capacity_high_water_bytes", std::to_string(diagnostics.patch_cell_receive_capacity_high_water_bytes)},
                    {"communication_workspace_high_water_bytes", std::to_string(diagnostics.communication_workspace_high_water_bytes)},
                    {"patch_cell_transport_round_count", std::to_string(diagnostics.patch_cell_transport_round_count)},
                    {"remote_patch_ghost_count", std::to_string(diagnostics.remote_patch_ghost_count)},
                    {"remote_interface_count", std::to_string(diagnostics.remote_interface_count)}}});
  }
}


struct CompactDmoMigrationCommitResult {
  MigrationExchangeStats exchange_stats{};
  std::size_t outbound_count = 0U;
  std::size_t inbound_count = 0U;
};

[[nodiscard]] CompactDmoMigrationCommitResult executeCompactDmoMigrationTransaction(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& scheduler,
    const parallel::RuntimeRebalancePlan& rebalance,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    const MigrationAdmissionPlan& plan,
    std::span<const std::uint64_t> expected_global_particle_ids,
    bool exact_ownership_audit_enabled) {
  if (plan.representation != MigrationTransactionRepresentation::kCompactHomogeneousDmo) {
    throw std::logic_error("compact DMO transaction invoked with generic admission plan");
  }
  const bool transfer_persistent_id_certificate =
      mpi_context.allreduceSumUint64(
          state.hasCertifiedPersistentParticleIds() ? 1ULL : 0ULL) ==
      static_cast<std::uint64_t>(mpi_context.worldSize());
  const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
  std::vector<std::vector<std::uint32_t>> selections_by_rank;
  std::vector<std::uint32_t> outbound_rows;
  core::ParticleSoa candidate_particles;
  core::ParticleSidecar candidate_sidecar;
  std::optional<core::HierarchicalTimeBinScheduler> candidate_scheduler;
  std::size_t kept_count = 0U;
  std::size_t final_count = 0U;

  std::exception_ptr local_preparation_failure;
  try {
    selections_by_rank.resize(rank_count);
    for (std::size_t peer = 0; peer < rank_count; ++peer) {
      selections_by_rank[peer].reserve(core::checkedIntegralNarrow<std::size_t>(
          plan.outbound_particle_count_by_rank[peer], "compact DMO per-peer selection reserve"));
    }
    for (const auto& intent : rebalance.particle_migrations) {
      if (intent.old_owner_rank != world_rank || intent.new_owner_rank == world_rank) continue;
      if (intent.new_owner_rank < 0 || intent.new_owner_rank >= mpi_context.worldSize() ||
          intent.item_index >= state.particles.size()) {
        throw std::runtime_error("compact DMO migration intent is outside authoritative local row/rank bounds");
      }
      const auto row = core::checkedIntegralNarrow<std::uint32_t>(
          intent.item_index, "compact DMO migration intent local row");
      if (state.particle_sidecar.particle_id[row] != intent.particle_id ||
          state.particleOwningRank(row) != static_cast<std::uint32_t>(world_rank)) {
        throw std::runtime_error("compact DMO migration intent does not match authoritative local identity");
      }
      selections_by_rank[static_cast<std::size_t>(intent.new_owner_rank)].push_back(row);
    }
    outbound_rows.reserve(core::checkedIntegralNarrow<std::size_t>(
        plan.outbound_particle_count, "compact DMO outbound row reserve"));
    for (auto& peer_rows : selections_by_rank) {
      std::sort(peer_rows.begin(), peer_rows.end());
      if (std::adjacent_find(peer_rows.begin(), peer_rows.end()) != peer_rows.end()) {
        throw std::runtime_error("compact DMO migration contains a duplicate local row for one destination");
      }
      outbound_rows.insert(outbound_rows.end(), peer_rows.begin(), peer_rows.end());
    }
    std::sort(outbound_rows.begin(), outbound_rows.end());
    if (std::adjacent_find(outbound_rows.begin(), outbound_rows.end()) != outbound_rows.end()) {
      throw std::runtime_error("compact DMO migration contains one local row with multiple destinations");
    }
    if (outbound_rows.size() != plan.outbound_particle_count) {
      throw std::runtime_error("compact DMO physical selection count differs from admission count");
    }

    kept_count = state.particles.size() - outbound_rows.size();
    final_count = core::checkedSizeAdd(
        kept_count,
        core::checkedIntegralNarrow<std::size_t>(plan.inbound_particle_count,
            "compact DMO inbound final count"),
        "compact DMO final candidate count");
    candidate_particles.position_x_comoving.resize(final_count);
    candidate_particles.position_y_comoving.resize(final_count);
    candidate_particles.position_z_comoving.resize(final_count);
    candidate_particles.velocity_x_peculiar.resize(final_count);
    candidate_particles.velocity_y_peculiar.resize(final_count);
    candidate_particles.velocity_z_peculiar.resize(final_count);
    candidate_particles.mass_code.resize(final_count);
    candidate_sidecar.particle_id.resize(final_count);
    if (plan.candidate_requires_sfc_lane) candidate_sidecar.sfc_key.resize(final_count);
    candidate_scheduler.emplace(0U);
    candidate_scheduler->reset(
        core::checkedIntegralNarrow<std::uint32_t>(final_count,
            "compact DMO scheduler candidate population"),
        0U, scheduler.currentTick());

    std::size_t outbound_cursor = 0U;
    std::size_t destination = 0U;
    for (std::size_t row = 0; row < state.particles.size(); ++row) {
      if (outbound_cursor < outbound_rows.size() && outbound_rows[outbound_cursor] == row) {
        ++outbound_cursor;
        continue;
      }
      candidate_sidecar.particle_id[destination] = state.particle_sidecar.particle_id[row];
      if (plan.candidate_requires_sfc_lane) candidate_sidecar.sfc_key[destination] = state.particleSfcKey(row);
      candidate_particles.position_x_comoving[destination] = state.particles.position_x_comoving[row];
      candidate_particles.position_y_comoving[destination] = state.particles.position_y_comoving[row];
      candidate_particles.position_z_comoving[destination] = state.particles.position_z_comoving[row];
      candidate_particles.velocity_x_peculiar[destination] = state.particles.velocity_x_peculiar[row];
      candidate_particles.velocity_y_peculiar[destination] = state.particles.velocity_y_peculiar[row];
      candidate_particles.velocity_z_peculiar[destination] = state.particles.velocity_z_peculiar[row];
      candidate_particles.mass_code[destination] = state.particles.mass_code[row];
      ++destination;
    }
    if (destination != kept_count || outbound_cursor != outbound_rows.size()) {
      throw std::logic_error("compact DMO kept-row streaming did not consume the admitted selection exactly");
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure, "compact DMO candidate preparation");

  CompactDmoExchangeResult exchange = exchangeRuntimeCompactDmoParticles(
      mpi_context, state, selections_by_rank, plan.inbound_particle_count_by_rank,
      candidate_particles, candidate_sidecar, kept_count, plan.candidate_requires_sfc_lane);

  std::exception_ptr local_validation_failure;
  try {
    std::sort(exchange.inbound_particle_ids.begin(), exchange.inbound_particle_ids.end());
    if (std::adjacent_find(exchange.inbound_particle_ids.begin(), exchange.inbound_particle_ids.end()) !=
        exchange.inbound_particle_ids.end()) {
      throw std::runtime_error("compact DMO migration received duplicate inbound particle IDs");
    }
    for (std::size_t row = 0; row < kept_count; ++row) {
      if (std::binary_search(exchange.inbound_particle_ids.begin(),
                             exchange.inbound_particle_ids.end(),
                             candidate_sidecar.particle_id[row])) {
        throw std::runtime_error(
            "compact DMO migration received an inbound particle ID that duplicates a kept local particle");
      }
    }
    if (!candidate_particles.isConsistent() || !candidate_sidecar.isConsistent()) {
      throw std::runtime_error("compact DMO migration candidate failed canonical lane validation");
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_validation_failure, "compact DMO candidate local validation");

  if (exact_ownership_audit_enabled) {
    const auto report = parallel::validateExactGlobalOwnershipPartition(
        mpi_context, candidate_sidecar.particle_id, expected_global_particle_ids);
    if (!report.valid()) {
      throw std::runtime_error(
          "compact DMO candidate exact ownership validation failed: duplicate=" +
          std::to_string(report.duplicate_count) + ", missing=" +
          std::to_string(report.missing_count) + ", extra=" +
          std::to_string(report.extra_count));
    }
  }

  // From here to publication there are no population-sized allocations. The
  // replacement rung-zero scheduler was prepared under the same reservation.
  static_assert(std::is_nothrow_move_assignable_v<core::HierarchicalTimeBinScheduler>);
  if (!candidate_scheduler.has_value()) {
    throw std::logic_error("compact DMO scheduler candidate was not prepared");
  }
  const std::uint32_t uniform_flags = state.homogeneousDmoParticleFlags();
  const double drift_time_code = state.homogeneousDmoLastDriftTimeCode();
  const double drift_scale_factor = state.homogeneousDmoLastDriftScaleFactor();
  state.commitCompactHomogeneousDmoCandidate(
      std::move(candidate_particles), std::move(candidate_sidecar),
      static_cast<std::uint32_t>(world_rank), uniform_flags,
      drift_time_code, drift_scale_factor, transfer_persistent_id_certificate);
  scheduler = std::move(*candidate_scheduler);
  return CompactDmoMigrationCommitResult{
      .exchange_stats = exchange.stats,
      .outbound_count = outbound_rows.size(),
      .inbound_count = exchange.inbound_particle_ids.size(),
  };
}

[[nodiscard]] bool applyMeasuredRuntimeRebalancePlan(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::SimulationConfig& config,
    const RuntimeServices& services,
    int world_rank,
    const parallel::DecompositionRuntimeMeasurements& measurements,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint64_t> expected_global_particle_ids,
    core::ProfilerSession* profiler,
    std::uint64_t step_index) {
  const parallel::MpiContext& mpi_context = services.mpi_context;
  if (!config.parallel.decomposition_runtime_rebalance_enabled || mpi_context.worldSize() <= 1) {
    return false;
  }
  if (world_rank < 0 || world_rank >= mpi_context.worldSize()) {
    throw std::invalid_argument("runtime rebalance world_rank is outside MPI world");
  }
  parallel::RuntimeRebalancePlan rebalance;
  core::MemoryReservation decomposition_reservation;
  {
    std::exception_ptr decomposition_admission_failure;
    try {
      if (services.memory_governor != nullptr) {
        // One authoritative source-view memory estimate shared by every
        // MemoryGovernor admission: active mask, compact patch index,
        // gas-incidence offsets/indices/construction scratch, and the bounded
        // species table. Gas-empty/DMO states charge zero gas-incidence bytes
        // through parallel::runtimeDecompositionHasGasIncidenceSource, the
        // same condition RuntimeDecompositionSourceStorage uses.
        const parallel::RuntimeDecompositionSourceMemoryEstimate source_estimate =
            parallel::estimateRuntimeDecompositionSourceStorage(
                state, active_particle_indices);
        const std::size_t entity_count = core::checkedSizeAdd(
            state.particles.size(), state.patches.size(),
            "runtime compact decomposition entity count");
        const std::size_t block_count = entity_count / 256U +
            (entity_count % 256U != 0U ? 1U : 0U);
         const std::size_t sample_count = 256U;
         const std::size_t rank_count = core::checkedIntegralNarrow<std::size_t>(
             mpi_context.worldSize(),
             "runtime compact decomposition rank count");
          const std::size_t sample_lane_count = core::checkedSizeAdd(
              core::checkedSizeMultiply(
                  rank_count, 3U, "runtime compact decomposition sample lane count"),
              1U, "runtime compact decomposition sample lane count");

         std::size_t planning_bytes = core::checkedIntegralNarrow<std::size_t>(
             source_estimate.total_bytes,
             "runtime compact source estimate byte width");
         planning_bytes = core::checkedSizeAdd(
             planning_bytes,
             core::checkedSizeMultiply(
                 entity_count, sizeof(parallel::CompactRuntimeDecompositionRecord),
                 "runtime compact decomposition record bytes"),
             "runtime compact decomposition planning bytes");
          const std::size_t sample_bytes = core::checkedSizeMultiply(
              core::checkedSizeMultiply(
                  sample_lane_count, sample_count,
                  "runtime compact decomposition sample count"),
              sizeof(std::uint64_t) * 3U,
              "runtime compact decomposition sample bytes");
          const std::size_t sample_payload_bytes = core::checkedSizeMultiply(
              core::checkedSizeMultiply(
                  rank_count, sample_count,
                  "runtime compact decomposition sample payload count"),
              sizeof(std::uint64_t) * 3U,
              "runtime compact decomposition sample payload bytes");
          const std::size_t round_limit = std::max<std::size_t>(
              1U, parallel::mpiTransportRoundLimitBytes());
          const std::size_t gather_round_count = std::max<std::size_t>(
              1U,
              (sample_payload_bytes + round_limit - 1U) / round_limit);
          const std::size_t gather_metadata_bytes = core::checkedSizeAdd(
              core::checkedSizeMultiply(
                  rank_count,
                  sizeof(std::uint64_t) * 4U +
                      sizeof(int) * 2U + sizeof(std::size_t),
                  "runtime compact bounded all-gather metadata bytes"),
              core::checkedSizeMultiply(
                  core::checkedSizeMultiply(
                      gather_round_count, rank_count,
                      "runtime compact bounded all-gather round rank count"),
                  sizeof(int) * 2U + sizeof(std::size_t),
                  "runtime compact bounded all-gather round metadata bytes"),
              "runtime compact bounded all-gather metadata bytes");
          const std::size_t repair_candidate_bytes = core::checkedSizeMultiply(
              rank_count, sizeof(std::uint64_t) * 3U,
              "runtime compact repair candidate bytes");
          planning_bytes = core::checkedSizeAdd(
              planning_bytes, gather_metadata_bytes,
              "runtime compact decomposition planning bytes");
          planning_bytes = core::checkedSizeAdd(
              planning_bytes, repair_candidate_bytes,
              "runtime compact decomposition planning bytes");

         planning_bytes = core::checkedSizeAdd(

            planning_bytes, sample_bytes,
            "runtime compact decomposition planning bytes");
        planning_bytes = core::checkedSizeAdd(
            planning_bytes,
            core::checkedSizeMultiply(
                core::checkedSizeAdd(
                    block_count, 1U, "runtime compact prefix block count"),
                sizeof(std::uint64_t),
                "runtime compact decomposition prefix bytes"),
            "runtime compact decomposition planning bytes");
         const std::size_t metric_lane_bytes = core::checkedSizeAdd(
             core::checkedSizeMultiply(
                 5U, sizeof(std::uint64_t),
                 "runtime compact decomposition metric lane width"),
             core::checkedSizeMultiply(
                 13U, sizeof(double),
                 "runtime compact decomposition metric lane width"),
             "runtime compact decomposition metric lane width");
         planning_bytes = core::checkedSizeAdd(
             planning_bytes,
             core::checkedSizeMultiply(
                 core::checkedSizeMultiply(
                     rank_count, 2U,
                     "runtime compact decomposition metric vector count"),
                 metric_lane_bytes,
                 "runtime compact decomposition metric bytes"),
             "runtime compact decomposition planning bytes");
         planning_bytes = core::checkedSizeAdd(
             planning_bytes,
             core::checkedSizeMultiply(
                 rank_count, sizeof(std::uint64_t) * 10U,
                 "runtime compact decomposition control bytes"),
             "runtime compact decomposition planning bytes");

        planning_bytes = core::checkedSizeAdd(
            planning_bytes,
            core::checkedSizeMultiply(
                state.particles.size(),
                sizeof(parallel::ParticleMigrationIntent),
                "runtime compact particle intent bound"),
            "runtime compact decomposition planning bytes");
        planning_bytes = core::checkedSizeAdd(
            planning_bytes,
            core::checkedSizeMultiply(
                state.patches.size(),
                sizeof(parallel::AmrPatchOwnershipUpdate),
                "runtime compact patch intent bound"),
            "runtime compact decomposition planning bytes");
        decomposition_reservation = services.memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            core::checkedIntegralNarrow<std::uint64_t>(
                planning_bytes, "runtime compact decomposition planning byte width"),
            "parallel.decomposition.plan");
        decomposition_reservation.commit();
      }
    } catch (...) {
      decomposition_admission_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        decomposition_admission_failure, "runtime compact decomposition memory admission");

    std::optional<RuntimeDecompositionSourceStorage> source_storage;
    parallel::RuntimeDecompositionSourceView source_view;
    std::uint64_t local_persistent_bytes = 0U;
    std::exception_ptr source_preparation_failure;
    try {
      source_storage.emplace(state, world_rank, active_particle_indices);
      source_view = source_storage->view();
      for (std::size_t particle = 0; particle < source_view.particle_count; ++particle) {
        const std::uint32_t species = source_view.particleSpeciesTag(particle);
        if (species >= source_view.particle_memory_bytes_by_species.size()) {
          throw std::invalid_argument("runtime compact source species exceeds memory table");
        }
        local_persistent_bytes = core::checkedMemoryBytesAdd(
            local_persistent_bytes,
            source_view.particle_memory_bytes_by_species[species],
            "runtime compact persistent particle memory");
      }
      for (const std::uint32_t patch : source_view.compact_patch_indices) {
        const std::uint64_t patch_bytes = core::checkedSizeMultiply(
            source_view.patch_cell_counts[patch],
            sizeof(double) * 8U + sizeof(std::uint32_t) * 2U,
            "runtime compact persistent patch memory");
        local_persistent_bytes = core::checkedMemoryBytesAdd(
            local_persistent_bytes, patch_bytes,
            "runtime compact persistent decomposition memory");
      }
    } catch (...) {
      source_preparation_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        source_preparation_failure, "runtime compact decomposition source preparation");

    const parallel::RuntimeRebalanceConfig rebalance_config{
        .world_size = mpi_context.worldSize(),
        .imbalance_trigger_ratio = config.parallel.decomposition_rebalance_imbalance_trigger,
        .memory_trigger_ratio = config.parallel.decomposition_rebalance_memory_trigger,
        .max_migrated_load_fraction = config.parallel.decomposition_rebalance_max_migrated_load_fraction,
        .allow_particle_migration = true,
        .allow_amr_patch_reassignment = true,
    };
    parallel::DecompositionConfig decomposition_config =
        makeWorkflowDecompositionConfig(config, mpi_context.worldSize());
    std::exception_ptr decomposition_envelope_failure;
    try {
      applyMemoryAwareDecompositionEnvelope(
          decomposition_config, local_persistent_bytes, config, services);
    } catch (...) {
      decomposition_envelope_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        decomposition_envelope_failure, "runtime compact decomposition memory envelope");
    rebalance = parallel::buildCompactDistributedRuntimeRebalancePlan(
        mpi_context,
        source_view,
        decomposition_config,
        rebalance_config,
        measurements,
        makeWorkflowFeedbackCoefficients(config));
  }
  rebalance.exact_debug_audit_enabled = config.parallel.decomposition_debug_exact_ownership_audit;
  if (!rebalance.should_rebalance) {
    exchangeAndValidateAmrPatchPayloads(state, mpi_context, world_rank, step_index, profiler);
    recordRuntimeRebalanceDecision(profiler, rebalance, step_index);
    return false;
  }

  const bool restore_compact_dmo_after_commit = state.hasHomogeneousDmoMetadata();
  const bool local_compact_dmo_eligible = compactDmoMigrationEligible(state, scheduler, rebalance);
  const bool all_ranks_compact_dmo_eligible =
      mpi_context.allreduceSumUint64(local_compact_dmo_eligible ? 1ULL : 0ULL) ==
      static_cast<std::uint64_t>(mpi_context.worldSize());
  const bool use_compact_dmo_migration = all_ranks_compact_dmo_eligible &&
      compactDmoMetadataConsensus(state, scheduler, mpi_context);
  const std::array transaction_reports{
      core::collectSimulationMemoryReport(state),
      core::collectSchedulerMemoryReport(scheduler, gas_cell_scheduler)};
  const std::uint64_t transaction_baseline_before =
      core::memoryReportBaselineOwnedBytes(
          core::mergeMemoryReports(transaction_reports));
  std::uint64_t process_baseline_before = 0U;
  MigrationAdmissionPlan migration_admission_plan;
  core::MemoryReservation migration_reservation;
  std::exception_ptr migration_admission_failure;
  try {
    migration_admission_plan = use_compact_dmo_migration
        ? buildCompactDmoMigrationAdmissionPlan(
              state, rebalance, mpi_context, world_rank,
              config.parallel.decomposition_debug_exact_ownership_audit)
        : buildMigrationAdmissionPlan(
              state, rebalance, mpi_context, world_rank, transaction_baseline_before);
    if (services.memory_governor != nullptr) {
      process_baseline_before =
          services.memory_governor->snapshot().baseline_owned_bytes;
      if (transaction_baseline_before > process_baseline_before) {
        throw std::logic_error(
            "migration governor baseline is stale: state/scheduler capacity exceeds process baseline");
      }
      migration_reservation = services.memory_governor->reserve(
          core::MemoryClass::kCommunication,
          migration_admission_plan.requested_extra_bytes,
          "parallel.migration.transaction");
      migration_reservation.commit();
    }
  } catch (...) {
    migration_admission_failure = std::current_exception();
  }
  FailureCoordinator(services).rethrowCollectiveFailure(
      migration_admission_failure, "runtime migration memory admission");

  const std::uint64_t particle_index_generation_before = state.particleIndexGeneration();
  const std::uint64_t cell_index_generation_before = state.cellIndexGeneration();
  const std::uint64_t gas_identity_generation_before = state.gasCellIdentityGeneration();

  if (use_compact_dmo_migration) {
    CompactDmoMigrationCommitResult compact_result;
    std::exception_ptr compact_failure;
    try {
      compact_result = executeCompactDmoMigrationTransaction(
          state, scheduler, rebalance, mpi_context, world_rank,
          migration_admission_plan, expected_global_particle_ids,
          config.parallel.decomposition_debug_exact_ownership_audit);
    } catch (...) {
      compact_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        compact_failure, "compact homogeneous-DMO migration transaction");

    recordRuntimeRebalanceDecision(profiler, rebalance, step_index);
    if (profiler != nullptr) {
      profiler->recordEvent(core::RuntimeEvent{
          .event_kind = "parallel.decomposition.runtime_migration_commit",
          .severity = core::RuntimeEventSeverity::kWarning,
          .subsystem = "parallel.domain_decomposition",
          .step_index = step_index,
          .message = "runtime rebalance committed compact homogeneous-DMO migration at a closed rung-zero boundary",
          .payload = {
              {"migration_representation", "compact_homogeneous_dmo"},
              {"outbound_particle_count", std::to_string(compact_result.outbound_count)},
              {"inbound_particle_count", std::to_string(compact_result.inbound_count)},
              {"candidate_state_capacity_bytes", std::to_string(migration_admission_plan.candidate_state_capacity_bytes)},
              {"selection_capacity_bytes", std::to_string(migration_admission_plan.selection_capacity_bytes)},
              {"packet_staging_bound_bytes", std::to_string(migration_admission_plan.packet_staging_bytes)},
              {"record_assembly_bound_bytes", std::to_string(migration_admission_plan.record_assembly_bound_bytes)},
              {"uniform_scheduler_capacity_bytes", std::to_string(migration_admission_plan.uniform_scheduler_capacity_bytes)},
              {"exact_ownership_workspace_bound_bytes", std::to_string(migration_admission_plan.exact_ownership_workspace_bound_bytes)},
              {"migration_transaction_reserved_bytes", std::to_string(migration_admission_plan.requested_extra_bytes)},
              {"wire_bytes_sent", std::to_string(compact_result.exchange_stats.wire_bytes_sent)},
              {"wire_bytes_received", std::to_string(compact_result.exchange_stats.wire_bytes_received)},
              {"physical_exchange_rounds", std::to_string(compact_result.exchange_stats.physical_exchange_rounds)},
          },
      });
    }
    if (migration_reservation.valid()) {
      const std::array reconciled_transaction_reports{
          core::collectSimulationMemoryReport(state),
          core::collectSchedulerMemoryReport(scheduler, gas_cell_scheduler)};
      const std::uint64_t transaction_baseline_after =
          core::memoryReportBaselineOwnedBytes(core::mergeMemoryReports(reconciled_transaction_reports));
      const std::uint64_t non_transaction_baseline = process_baseline_before - transaction_baseline_before;
      const std::uint64_t reconciled_process_baseline = core::checkedMemoryBytesAdd(
          non_transaction_baseline, transaction_baseline_after,
          "compact DMO migration retained-capacity baseline reconciliation");
      migration_reservation.reconcileBaselineOwnedAndRelease(reconciled_process_baseline);
    }
    return state.particleIndexGeneration() != particle_index_generation_before ||
        state.cellIndexGeneration() != cell_index_generation_before ||
        state.gasCellIdentityGeneration() != gas_identity_generation_before;
  }

  if (restore_compact_dmo_after_commit) {
    state.materializeParticleMetadata();
  }

  std::unordered_map<std::uint64_t, std::uint32_t> local_index_by_particle_id;
  local_index_by_particle_id.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    local_index_by_particle_id.emplace(
        state.particle_sidecar.particle_id[particle_index],
        static_cast<std::uint32_t>(particle_index));
  }

  std::unordered_map<std::uint32_t, int> outbound_target_by_local_index;
  outbound_target_by_local_index.reserve(rebalance.particle_migrations.size());
  std::unordered_map<std::uint32_t, int> outbound_patch_target_by_local_index;
  outbound_patch_target_by_local_index.reserve(rebalance.amr_patch_ownership_updates.size());
  const auto add_particle_migration = [&](std::uint32_t local_index, int new_owner_rank, std::string_view source_label) {
    if (local_index >= state.particles.size()) {
      throw std::out_of_range("runtime rebalance attempted to migrate a particle index outside SimulationState");
    }
    if (new_owner_rank < 0 || new_owner_rank >= mpi_context.worldSize()) {
      throw std::invalid_argument("runtime rebalance produced particle migration target outside MPI world");
    }
    if (new_owner_rank == world_rank) {
      return;
    }
    if (state.particleOwningRank(local_index) != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime rebalance attempted to migrate a non-authoritative local particle via " + std::string(source_label));
    }
     const auto [it, inserted] = outbound_target_by_local_index.emplace(local_index, new_owner_rank);
     if (!inserted && it->second != new_owner_rank) {
       if (source_label == "amr_patch_ownership_update") {
         it->second = new_owner_rank;
         return;
       }
       throw std::runtime_error("runtime rebalance produced conflicting destinations for one particle");
     }

  };
  const auto add_patch_migration = [&](std::uint32_t local_patch_index, int new_owner_rank) {
    if (local_patch_index >= state.patches.size()) {
      throw std::out_of_range("runtime rebalance attempted to migrate an AMR patch index outside SimulationState");
    }
    if (new_owner_rank < 0 || new_owner_rank >= mpi_context.worldSize()) {
      throw std::invalid_argument("runtime rebalance produced AMR patch migration target outside MPI world");
    }
    if (new_owner_rank == world_rank) {
      return;
    }
    if (state.patches.owning_rank[local_patch_index] != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime rebalance attempted to migrate a non-authoritative local AMR patch");
    }
    const auto [it, inserted] = outbound_patch_target_by_local_index.emplace(local_patch_index, new_owner_rank);
    if (!inserted && it->second != new_owner_rank) {
      throw std::runtime_error("runtime rebalance produced conflicting destinations for one AMR patch");
    }
  };

  for (const auto& intent : rebalance.particle_migrations) {
    if (intent.old_owner_rank != world_rank || intent.new_owner_rank == world_rank) {
      continue;
    }
    const auto found = local_index_by_particle_id.find(intent.particle_id);
    if (found == local_index_by_particle_id.end()) {
      throw std::runtime_error("runtime rebalance migration intent references a particle ID missing on the owning rank");
    }
    add_particle_migration(found->second, intent.new_owner_rank, "particle_sfc_intent");
  }

  for (const auto& update : rebalance.amr_patch_ownership_updates) {
    if (update.new_owner_rank < 0 || update.new_owner_rank >= mpi_context.worldSize()) {
      throw std::runtime_error("runtime rebalance produced AMR patch owner outside MPI world");
    }
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      if (state.patches.patch_id[patch_index] != update.patch_id) {
        continue;
      }
      if (state.patches.owning_rank[patch_index] == static_cast<std::uint32_t>(world_rank) &&
          update.new_owner_rank != world_rank) {
        add_patch_migration(static_cast<std::uint32_t>(patch_index), update.new_owner_rank);
        const std::uint32_t first_cell = state.patches.first_cell[patch_index];
        const std::uint32_t cell_count = state.patches.cell_count[patch_index];
        if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) > state.cells.size()) {
          throw std::runtime_error("runtime rebalance AMR patch cell range is outside CellSoa");
        }
        for (std::uint32_t cell_offset = 0; cell_offset < cell_count; ++cell_offset) {
          const std::uint32_t cell_index = first_cell + cell_offset;
          const auto gas_particle_index = parentParticleRowForGasCellRow(
              state,
              cell_index,
              local_index_by_particle_id,
              "runtime rebalance AMR patch ownership update");
          if (gas_particle_index.has_value() &&
              state.particleOwningRank(*gas_particle_index) == static_cast<std::uint32_t>(world_rank)) {
            add_particle_migration(*gas_particle_index, update.new_owner_rank, "amr_patch_ownership_update");
          }
        }
      }
    }
  }

  std::vector<std::uint32_t> outbound_local_indices;
  outbound_local_indices.reserve(outbound_target_by_local_index.size());
  for (const auto& [local_index, target_rank] : outbound_target_by_local_index) {
    (void)target_rank;
    outbound_local_indices.push_back(local_index);
  }
  std::sort(outbound_local_indices.begin(), outbound_local_indices.end());

  std::vector<std::uint32_t> outbound_patch_indices;
  outbound_patch_indices.reserve(outbound_patch_target_by_local_index.size());
  for (const auto& [patch_index, target_rank] : outbound_patch_target_by_local_index) {
    (void)target_rank;
    outbound_patch_indices.push_back(patch_index);
  }
  std::sort(outbound_patch_indices.begin(), outbound_patch_indices.end());

  std::vector<std::uint32_t> kept_cell_rows;
  kept_cell_rows.reserve(state.cells.size());
  for (std::uint32_t cell_row = 0; cell_row < state.cells.size(); ++cell_row) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_row];
    if (patch_index >= state.patches.size()) {
      throw std::runtime_error("runtime rebalance found a gas cell with invalid patch_index");
    }
    if (!std::binary_search(outbound_patch_indices.begin(), outbound_patch_indices.end(), patch_index)) {
      kept_cell_rows.push_back(cell_row);
    }
  }
  std::vector<core::TimeBinSchedulerIdentityRecord> gas_scheduler_records =
      core::exportGasCellSchedulerIdentityRecords(gas_cell_scheduler, state, kept_cell_rows);

  std::vector<std::uint32_t> preserved_indices;
  preserved_indices.reserve(state.particles.size() - outbound_local_indices.size());
  for (std::uint32_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    if (!std::binary_search(outbound_local_indices.begin(), outbound_local_indices.end(), particle_index)) {
      preserved_indices.push_back(particle_index);
    }
  }
  std::vector<core::TimeBinSchedulerIdentityRecord> scheduler_records =
      core::exportParticleSchedulerIdentityRecords(
          scheduler, state, preserved_indices);

  std::vector<std::vector<core::ParticleMigrationRecord>> outbound_records_by_rank(
      static_cast<std::size_t>(mpi_context.worldSize()));
  for (std::size_t rank = 0; rank < outbound_records_by_rank.size(); ++rank) {
    outbound_records_by_rank[rank].reserve(
        core::checkedIntegralNarrow<std::size_t>(
            migration_admission_plan.outbound_particle_record_capacity_count_by_rank[rank],
            "outbound particle migration record reserve"));
  }
  for (const std::uint32_t local_index : outbound_local_indices) {
    const int target_rank = outbound_target_by_local_index.at(local_index);
    auto records = state.packParticleMigrationRecords(std::span<const std::uint32_t>(&local_index, 1), scheduler);
    if (records.size() != 1U) {
      throw std::runtime_error("runtime rebalance particle migration packing returned an unexpected record count");
    }
    records[0].owning_rank = static_cast<std::uint32_t>(target_rank);
    outbound_records_by_rank[static_cast<std::size_t>(target_rank)].push_back(std::move(records[0]));
  }

  std::vector<std::vector<core::AmrPatchMigrationRecord>> outbound_patch_records_by_rank(
      static_cast<std::size_t>(mpi_context.worldSize()));
  for (std::size_t rank = 0; rank < outbound_patch_records_by_rank.size(); ++rank) {
    outbound_patch_records_by_rank[rank].reserve(
        core::checkedIntegralNarrow<std::size_t>(
            migration_admission_plan.outbound_patch_count_by_rank[rank],
            "outbound AMR migration record reserve"));
  }
  for (const std::uint32_t patch_index : outbound_patch_indices) {
    const int target_rank = outbound_patch_target_by_local_index.at(patch_index);
    auto records = state.packAmrPatchMigrationRecords(std::span<const std::uint32_t>(&patch_index, 1));
    if (records.size() != 1U) {
      throw std::runtime_error("runtime rebalance AMR patch migration packing returned an unexpected record count");
    }
    records[0].patch.owning_rank = static_cast<std::uint32_t>(target_rank);
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    std::vector<std::uint32_t> patch_cell_rows;
    patch_cell_rows.reserve(cell_count);
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      patch_cell_rows.push_back(first_cell + offset);
    }
    const std::vector<core::TimeBinSchedulerIdentityRecord> patch_scheduler_records =
        core::exportGasCellSchedulerIdentityRecords(gas_cell_scheduler, state, patch_cell_rows);
    records[0].gas_cell_scheduler_records.reserve(patch_scheduler_records.size());
    for (const auto& scheduler_record : patch_scheduler_records) {
      records[0].gas_cell_scheduler_records.push_back(core::GasCellSchedulerMigrationRecord{
          .gas_cell_id = scheduler_record.element_id,
          .bin_index = scheduler_record.bin_index,
          .next_activation_tick = scheduler_record.next_activation_tick,
          .pending_bin_index = scheduler_record.pending_bin_index,
      });
    }
    for (core::GasCellMigrationRecord& gas_record : records[0].gas_cell_records) {
      gas_record.owning_rank = static_cast<std::uint32_t>(target_rank);
      gas_record.fields.patch_index = 0U;
      gas_record.fields.owning_patch_id = records[0].patch.patch_id;
    }
    outbound_patch_records_by_rank[static_cast<std::size_t>(target_rank)].push_back(std::move(records[0]));
  }

  auto particle_exchange =
      exchangeRuntimeParticleMigrationRecords(mpi_context, outbound_records_by_rank);
  auto patch_exchange =
      exchangeRuntimeAmrPatchMigrationRecords(mpi_context, outbound_patch_records_by_rank);
  std::vector<core::ParticleMigrationRecord> inbound_records =
      std::move(particle_exchange.records);
  std::vector<core::AmrPatchMigrationRecord> inbound_patch_records =
      std::move(patch_exchange.records);
  const std::size_t inbound_particle_count = inbound_records.size();
  const std::size_t inbound_patch_count = inbound_patch_records.size();
  for (const core::ParticleMigrationRecord& record : inbound_records) {
    if (record.owning_rank != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime particle migration exchange delivered a record to the wrong destination rank");
    }
  }
  for (const core::AmrPatchMigrationRecord& record : inbound_patch_records) {
    if (record.patch.owning_rank != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime AMR patch migration exchange delivered a patch to the wrong destination rank");
    }
    if (record.gas_cell_scheduler_records.size() != record.gas_cell_records.size()) {
      throw std::runtime_error("runtime AMR patch migration payload is missing gas-cell scheduler identity records");
    }
    for (const core::GasCellSchedulerMigrationRecord& scheduler_record : record.gas_cell_scheduler_records) {
      gas_scheduler_records.push_back(core::TimeBinSchedulerIdentityRecord{
          .element_id = scheduler_record.gas_cell_id,
          .bin_index = scheduler_record.bin_index,
          .next_activation_tick = scheduler_record.next_activation_tick,
          .pending_bin_index = scheduler_record.pending_bin_index,
      });
    }
  }

  if (!outbound_patch_indices.empty() || !inbound_patch_records.empty()) {
    core::AmrPatchMigrationCommit patch_commit;
    patch_commit.world_rank = world_rank;
    patch_commit.outbound_local_patch_indices = outbound_patch_indices;
    patch_commit.inbound_records = std::move(inbound_patch_records);
    state.commitAmrPatchMigration(patch_commit);
    core::rebuildSchedulerFromGasCellIdentityRecords(gas_cell_scheduler, gas_scheduler_records, state);
    core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, state);
  }

  if (!outbound_local_indices.empty() || !inbound_records.empty()) {
    scheduler_records.reserve(scheduler_records.size() + inbound_records.size());
    for (const core::ParticleMigrationRecord& record : inbound_records) {
      if (!record.has_scheduler_fields) {
        throw std::runtime_error(
            "runtime particle migration payload is missing scheduler authority");
      }
      scheduler_records.push_back(core::TimeBinSchedulerIdentityRecord{
          .element_id = record.particle_id,
          .bin_index = record.scheduler_fields.bin_index,
          .next_activation_tick = record.scheduler_fields.next_activation_tick,
          .pending_bin_index = record.scheduler_fields.pending_bin_index,
      });
    }
    core::ParticleMigrationCommit commit;
    commit.world_rank = world_rank;
    commit.outbound_local_indices = outbound_local_indices;
    commit.inbound_records = std::move(inbound_records);
    commit.preserve_gas_cell_state = true;
    state.commitParticleMigration(commit);

    std::vector<std::uint64_t> destination_particle_ids(state.particle_sidecar.particle_id.begin(),
                                                        state.particle_sidecar.particle_id.end());
    core::rebuildSchedulerFromParticleIdentityRecords(
        scheduler, scheduler_records, destination_particle_ids);
    syncTimeBinsFromScheduler(scheduler, state);
    if (config.parallel.decomposition_debug_exact_ownership_audit) {
      requireGlobalOwnedParticlePartitionIdentity(
          state, mpi_context, expected_global_particle_ids, "runtime rebalance particle migration commit");
    }
  }
  exchangeAndValidateAmrPatchPayloads(state, mpi_context, world_rank, step_index, profiler);
  if (restore_compact_dmo_after_commit) {
    if (!state.compactHomogeneousDmoMetadata(static_cast<std::uint32_t>(world_rank))) {
      throw std::runtime_error(
          "runtime migration could not restore the pre-transaction homogeneous-DMO metadata representation");
    }
  }

  recordRuntimeRebalanceDecision(profiler, rebalance, step_index);
  if (profiler != nullptr && (!outbound_local_indices.empty() || inbound_particle_count != 0U)) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.runtime_migration_commit",
        .severity = core::RuntimeEventSeverity::kWarning,
        .subsystem = "parallel.domain_decomposition",
        .step_index = step_index,
        .message = "runtime rebalance committed authoritative particle migration at a safe scheduler boundary",
        .payload = {
            {"migration_representation", "generic_full_physics"},
            {"outbound_particle_count", std::to_string(outbound_local_indices.size())},
            {"inbound_particle_count", std::to_string(inbound_particle_count)},
            {"inbound_patch_count", std::to_string(inbound_patch_count)},
            {"migration_transaction_reserved_bytes", std::to_string(migration_admission_plan.requested_extra_bytes)},
            {"migration_record_capacity_bytes", std::to_string(migration_admission_plan.record_capacity_bytes)},
            {"migration_packet_staging_bound_bytes", std::to_string(migration_admission_plan.packet_staging_bytes)},
            {"migration_scheduler_remap_capacity_bytes", std::to_string(migration_admission_plan.scheduler_remap_bytes)},
            {"migration_index_map_capacity_bytes", std::to_string(migration_admission_plan.index_map_bytes)},
            {"migration_metadata_materialization_bytes", std::to_string(migration_admission_plan.metadata_materialization_bytes)},
            {"migration_outbound_wire_traffic_upper_bound_bytes", std::to_string(migration_admission_plan.outbound_wire_bytes)},
            {"migration_inbound_wire_traffic_upper_bound_bytes", std::to_string(migration_admission_plan.inbound_wire_bytes)},
            {"migration_packet_send_capacity_bytes", std::to_string(std::max(
                 particle_exchange.stats.packet_send_capacity_bytes,
                 patch_exchange.stats.packet_send_capacity_bytes))},
            {"migration_packet_receive_capacity_bytes", std::to_string(std::max(
                 particle_exchange.stats.packet_receive_capacity_bytes,
                 patch_exchange.stats.packet_receive_capacity_bytes))},
            {"migration_communication_high_water_bytes", std::to_string(std::max(
                 particle_exchange.stats.communication_high_water_bytes,
                 patch_exchange.stats.communication_high_water_bytes))},
            {"migration_wire_bytes_sent", std::to_string(core::checkedMemoryBytesAdd(
                 particle_exchange.stats.wire_bytes_sent,
                 patch_exchange.stats.wire_bytes_sent,
                 "migration combined wire bytes sent"))},
            {"migration_wire_bytes_received", std::to_string(core::checkedMemoryBytesAdd(
                 particle_exchange.stats.wire_bytes_received,
                 patch_exchange.stats.wire_bytes_received,
                 "migration combined wire bytes received"))},
            {"current_weighted_imbalance_ratio", std::to_string(rebalance.current_metrics.weighted_imbalance_ratio)},
            {"target_weighted_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.weighted_imbalance_ratio)},
            {"current_memory_imbalance_ratio", std::to_string(rebalance.current_metrics.memory_imbalance_ratio)},
            {"target_memory_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.memory_imbalance_ratio)},
        },
    });
  }
  if (migration_reservation.valid()) {
    const std::array reconciled_transaction_reports{
        core::collectSimulationMemoryReport(state),
        core::collectSchedulerMemoryReport(scheduler, gas_cell_scheduler)};
    const std::uint64_t transaction_baseline_after =
        core::memoryReportBaselineOwnedBytes(
            core::mergeMemoryReports(reconciled_transaction_reports));
    const std::uint64_t non_transaction_baseline =
        process_baseline_before - transaction_baseline_before;
    const std::uint64_t reconciled_process_baseline =
        core::checkedMemoryBytesAdd(
            non_transaction_baseline, transaction_baseline_after,
            "runtime migration retained-capacity baseline reconciliation");
    migration_reservation.reconcileBaselineOwnedAndRelease(
        reconciled_process_baseline);
  }

  const bool local_decomposition_changed =
      state.particleIndexGeneration() != particle_index_generation_before ||
      state.cellIndexGeneration() != cell_index_generation_before ||
      state.gasCellIdentityGeneration() != gas_identity_generation_before;
  const std::uint64_t changed_rank_count = mpi_context.allreduceSumUint64(
      local_decomposition_changed ? 1ULL : 0ULL);
  return changed_rank_count > 0U;
}

// Startup density/PM occupancy grids and per-entity startup weight helpers.
// These are transient per-entity computations; only the ≤64-byte compact
// record fields survive record creation. Formulas are preserved verbatim from
// the historical rich-item startup path so ownership cuts do not shift.
constexpr std::size_t k_startup_density_grid = 16;

[[nodiscard]] double wrapStartupPeriodicCoordinate(double x, double box) {
  if (box <= 0.0) {
    return x;
  }
  double wrapped = std::fmod(x, box);
  if (wrapped < 0.0) {
    wrapped += box;
  }
  return wrapped;
}

[[nodiscard]] std::size_t startupDensityCellIndex(
    const core::SimulationConfig& config, double x, double y, double z) {
  const double box_x = config.cosmology.box_size_x_mpc_comoving;
  const double box_y = config.cosmology.box_size_y_mpc_comoving;
  const double box_z = config.cosmology.box_size_z_mpc_comoving;
  const std::size_t ix = (box_x > 0.0)
      ? std::min<std::size_t>(
            k_startup_density_grid - 1U,
            static_cast<std::size_t>(
                (wrapStartupPeriodicCoordinate(x, box_x) / box_x) *
                static_cast<double>(k_startup_density_grid)))
      : 0U;
  const std::size_t iy = (box_y > 0.0)
      ? std::min<std::size_t>(
            k_startup_density_grid - 1U,
            static_cast<std::size_t>(
                (wrapStartupPeriodicCoordinate(y, box_y) / box_y) *
                static_cast<double>(k_startup_density_grid)))
      : 0U;
  const std::size_t iz = (box_z > 0.0)
      ? std::min<std::size_t>(
            k_startup_density_grid - 1U,
            static_cast<std::size_t>(
                (wrapStartupPeriodicCoordinate(z, box_z) / box_z) *
                static_cast<double>(k_startup_density_grid)))
      : 0U;
  return (ix * k_startup_density_grid + iy) * k_startup_density_grid + iz;
}

[[nodiscard]] std::size_t startupPmXIndex(
    const core::SimulationConfig& config, double x, std::size_t pm_x_bins) {
  const double box_x = config.cosmology.box_size_x_mpc_comoving;
  if (box_x <= 0.0 || pm_x_bins == 0) {
    return std::size_t{0};
  }
  const double scaled = (wrapStartupPeriodicCoordinate(x, box_x) / box_x) *
      static_cast<double>(pm_x_bins);
  return std::min<std::size_t>(pm_x_bins - 1U, static_cast<std::size_t>(scaled));
}

struct StartupDensityOccupancyGrids {
  std::vector<std::uint32_t> occupancy;
  std::vector<std::uint32_t> active_occupancy;
  std::vector<std::uint32_t> gas_occupancy;
  std::vector<std::uint32_t> pm_x_occupancy;
};

[[nodiscard]] StartupDensityOccupancyGrids buildStartupDensityOccupancyGrids(
    const core::SimulationState& state,
    const core::SimulationConfig& config) {
  constexpr std::size_t k_grid_cells =
      k_startup_density_grid * k_startup_density_grid * k_startup_density_grid;
  StartupDensityOccupancyGrids grids;
  grids.occupancy.assign(k_grid_cells, 0U);
  grids.active_occupancy.assign(k_grid_cells, 0U);
  grids.gas_occupancy.assign(k_grid_cells, 0U);
  grids.pm_x_occupancy.assign(
      static_cast<std::size_t>(std::max(config.numerics.treepm_pm_grid_nx, 1)), 0U);
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const std::size_t cell = startupDensityCellIndex(
        config,
        state.particles.position_x_comoving[particle_index],
        state.particles.position_y_comoving[particle_index],
        state.particles.position_z_comoving[particle_index]);
    ++grids.occupancy[cell];
    if (!state.particles.time_bin.empty() && state.particles.time_bin[particle_index] == 0U) {
      ++grids.active_occupancy[cell];
    }
    if (state.particleSpeciesTag(particle_index) ==
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
      ++grids.gas_occupancy[cell];
    }
    ++grids.pm_x_occupancy[startupPmXIndex(
        config, state.particles.position_x_comoving[particle_index],
        grids.pm_x_occupancy.size())];
  }
  return grids;
}

// Startup-specific conservative particle footprint. Intentionally independent
// of estimateParticleMemoryBytesForDecomposition so startup weighting
// semantics are unchanged by module-sidecar accounting.
[[nodiscard]] std::uint64_t startupParticleMemoryBytes(
    const core::SimulationState& state, std::uint32_t species_tag) {
  std::uint64_t bytes = sizeof(double) * 7U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 3U;
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    bytes += sizeof(double);
  }
  if (!state.particle_sidecar.has_gravity_softening_override.empty()) {
    bytes += sizeof(std::uint8_t);
  }
  if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
    bytes += sizeof(double) * 8U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t);
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar)) {
    bytes += sizeof(std::uint32_t) + sizeof(double) * 13U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
    bytes += sizeof(std::uint32_t) * 2U + sizeof(double) * 8U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)) {
    bytes += sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U + sizeof(double) * 3U;
  }
  return bytes;
}

// One entity's startup weight inputs, computed as a transient per-entity
// operation. Only the ≤64-byte record fields are copied into the compact
// record; work components remain transient and are recomputed for the metric
// pass once cuts are known.
struct StartupDecompositionWeight {
  parallel::DecompositionWorkComponents work_components{};
  double weighted_load = 0.0;
  std::uint64_t entity_id = 0;
  std::uint64_t sfc_key = 0;
  std::uint64_t memory_bytes = 0;
  std::uint64_t active_target_count_recent = 0;
  std::uint64_t remote_tree_interactions_recent = 0;
  int current_owner_rank = -1;
  parallel::DecompositionEntityKind kind = parallel::DecompositionEntityKind::kParticle;
};

[[nodiscard]] StartupDecompositionWeight computeInitialDecompositionWeightForParticle(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    const parallel::DecompositionConfig& decomposition_config,
    const StartupDensityOccupancyGrids& grids,
    std::span<const std::uint32_t> patch_cell_count,
    std::size_t particle_index) {
  StartupDecompositionWeight weight;
  weight.kind = parallel::DecompositionEntityKind::kParticle;
  weight.entity_id = state.particle_sidecar.particle_id[particle_index];
  weight.current_owner_rank = static_cast<int>(state.particleOwningRank(particle_index));
  const double x = state.particles.position_x_comoving[particle_index];
  const double y = state.particles.position_y_comoving[particle_index];
  const double z = state.particles.position_z_comoving[particle_index];
  weight.sfc_key = parallel::sfcKeyForPosition(x, y, z, decomposition_config);
  const std::size_t density_cell = startupDensityCellIndex(config, x, y, z);
  const std::uint32_t local_density = grids.occupancy[density_cell];
  const std::uint32_t local_active = grids.active_occupancy[density_cell];
  const std::uint32_t local_gas = grids.gas_occupancy[density_cell];
  const std::uint32_t pm_load =
      grids.pm_x_occupancy[startupPmXIndex(config, x, grids.pm_x_occupancy.size())];
  const std::uint32_t species_tag = state.particleSpeciesTag(particle_index);
  double amr_patch_cost = 0.0;
  if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) &&
      !state.cells.patch_index.empty()) {
    std::vector<std::uint32_t> seen_patch_indices;
    for (const std::uint32_t cell_row :
         state.gas_cell_identity.rowsForParentParticleId(weight.entity_id)) {
      if (cell_row >= state.cells.patch_index.size()) {
        continue;
      }
      const std::uint32_t patch_index = state.cells.patch_index[cell_row];
      if (patch_index >= patch_cell_count.size() ||
          std::find(seen_patch_indices.begin(), seen_patch_indices.end(), patch_index) !=
              seen_patch_indices.end()) {
        continue;
      }
      seen_patch_indices.push_back(patch_index);
      amr_patch_cost += static_cast<double>(patch_cell_count[patch_index]);
    }
  }
  const double local_density_d = static_cast<double>(std::max<std::uint32_t>(local_density, 1U));
  weight.active_target_count_recent = local_active;
  weight.remote_tree_interactions_recent = static_cast<std::uint64_t>(
      std::llround(local_density_d * std::log2(local_density_d + 1.0)));
  weight.memory_bytes = startupParticleMemoryBytes(state, species_tag);
  weight.work_components = parallel::DecompositionWorkComponents{
      .particle_count_cost = 1.0,
      .gas_cell_cost = (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas))
          ? (1.0 + static_cast<double>(local_gas))
          : 0.0,
      .tree_interaction_cost = static_cast<double>(weight.remote_tree_interactions_recent),
      .pm_mesh_cost = static_cast<double>(pm_load),
      .amr_patch_cost = amr_patch_cost,
      .active_fraction_cost = static_cast<double>(local_active),
      .memory_pressure_cost = static_cast<double>(weight.memory_bytes),
      .transient_memory_cost = (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas))
          ? static_cast<double>(local_active) *
                static_cast<double>(hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell)
          : 0.0,
      .source_event_cost = (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar) ||
                            species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole))
          ? 1.0
          : 0.0,
      .communication_cost = static_cast<double>(weight.remote_tree_interactions_recent) +
          static_cast<double>(pm_load),
      .gpu_occupancy_cost = 0.0,
      .generic_work_cost = 1.0 + std::sqrt(local_density_d),
      .has_explicit_components = true,
  };
  weight.weighted_load = parallel::weightedLoadFromExplicitComponents(
      weight.work_components, weight.kind, weight.active_target_count_recent,
      weight.remote_tree_interactions_recent, /*work_units=*/1.0, weight.memory_bytes,
      decomposition_config);
  return weight;
}

// Returns nullopt for patches the historical startup path skipped: zero-cell
// patches and patches without an authoritative first cell.
[[nodiscard]] std::optional<StartupDecompositionWeight> computeInitialDecompositionWeightForPatch(
    const core::SimulationState& state,
    const parallel::DecompositionConfig& decomposition_config,
    std::size_t patch_index) {
  if (state.patches.cell_count[patch_index] == 0U) {
    return std::nullopt;
  }
  const std::uint32_t first_cell = state.patches.first_cell[patch_index];
  if (first_cell >= state.cells.size()) {
    return std::nullopt;
  }
  StartupDecompositionWeight weight;
  weight.kind = parallel::DecompositionEntityKind::kAmrPatch;
  weight.entity_id = state.patches.patch_id[patch_index];
  weight.current_owner_rank = static_cast<int>(state.patches.owning_rank[patch_index]);
  const double x = state.cells.center_x_comoving[first_cell];
  const double y = state.cells.center_y_comoving[first_cell];
  const double z = state.cells.center_z_comoving[first_cell];
  weight.sfc_key = parallel::sfcKeyForPosition(x, y, z, decomposition_config);
  weight.memory_bytes = static_cast<std::uint64_t>(state.patches.cell_count[patch_index]) *
      static_cast<std::uint64_t>(sizeof(double) * 8U + sizeof(std::uint32_t));
  weight.work_components = parallel::DecompositionWorkComponents{
      .amr_patch_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
          (1.0 + static_cast<double>(std::max(state.patches.level[patch_index], 0))),
      .memory_pressure_cost = static_cast<double>(weight.memory_bytes),
      .transient_memory_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
          static_cast<double>(hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell),
      .communication_cost = 2.0 * (
          static_cast<double>(state.patches.cell_dim_x[patch_index]) * state.patches.cell_dim_y[patch_index] +
          static_cast<double>(state.patches.cell_dim_x[patch_index]) * state.patches.cell_dim_z[patch_index] +
          static_cast<double>(state.patches.cell_dim_y[patch_index]) * state.patches.cell_dim_z[patch_index]),
      .generic_work_cost = static_cast<double>(state.patches.cell_count[patch_index]),
      .has_explicit_components = true,
  };
  weight.weighted_load = parallel::weightedLoadFromExplicitComponents(
      weight.work_components, weight.kind, /*active_target_count_recent=*/0U,
      /*remote_tree_interactions_recent=*/0U, /*work_units=*/1.0, weight.memory_bytes,
      decomposition_config);
  return weight;
}

void applyInitialGravityAwareDecomposition(
    core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_size,
    int world_rank,
    core::ProfilerSession* profiler,
    const RuntimeServices& services) {
  if (world_size <= 1) {
    return;
  }

  // MemoryGovernor admission for the compact startup planner's transient
  // live set, obtained before any planner population (occupancy grids, patch
  // mapping, compact records, plan arrays) is allocated. This answers the
  // planner-safety question only; it does not replace resulting-cut
  // feasibility checks, which startup does not apply (no
  // max_rank_memory_bytes ceiling on the initial placement). The reservation
  // is declared before every planner temporary so its destruction releases
  // only after records, occupancy grids, patch mapping, and plan arrays are
  // destroyed. The reservation is additional transient planner memory, not a
  // re-charge of canonical simulation state, so it is released (never
  // baseline-reconciled) at that lifetime boundary.
  core::MemoryReservation planner_reservation;
  std::exception_ptr planner_admission_failure;
  try {
    if (services.memory_governor != nullptr) {
      const std::size_t entity_upper_bound = core::checkedSizeAdd(
          state.particles.size(), state.patches.size(),
          "startup compact planner entity upper bound");
      const parallel::CompactStartupPlannerMemoryEstimate planner_estimate =
          parallel::estimateCompactStartupPlannerTransientBytes(
              entity_upper_bound,
              state.patches.size(),
              state.cells.size(),
              static_cast<std::size_t>(world_size),
              static_cast<std::size_t>(
                  std::max(config.numerics.treepm_pm_grid_nx, 1)),
              static_cast<std::size_t>(k_startup_density_grid) *
                  k_startup_density_grid * k_startup_density_grid);
      planner_reservation = services.memory_governor->reserve(
          core::MemoryClass::kPhaseResident,
          planner_estimate.total_bytes,
          "parallel.decomposition.startup_planner");
      planner_reservation.commit();
    }
  } catch (...) {
    planner_admission_failure = std::current_exception();
  }
  FailureCoordinator(services).rethrowCollectiveFailure(
      planner_admission_failure,
      "startup compact decomposition memory admission");

  const StartupDensityOccupancyGrids grids =
      buildStartupDensityOccupancyGrids(state, config);

  std::vector<std::uint32_t> patch_cell_count(state.patches.size(), 0U);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_index];
    if (patch_index < patch_cell_count.size()) {
      ++patch_cell_count[patch_index];
    }
  }

  const parallel::DecompositionConfig decomposition_config =
      makeWorkflowDecompositionConfig(config, world_size);

  // Compact startup placement: one ≤64-byte record per decomposition unit,
  // sorted in place by the shared compact SFC planner. No
  // DecompositionItem[N] population, no geometry retention, no per-entity
  // work_components storage.
  std::vector<parallel::CompactRuntimeDecompositionRecord> records;
  records.reserve(state.particles.size() + state.patches.size());
  std::vector<std::uint32_t> included_patch_rows;
  included_patch_rows.reserve(state.patches.size());

  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const StartupDecompositionWeight weight = computeInitialDecompositionWeightForParticle(
        state, config, decomposition_config, grids, patch_cell_count, particle_index);
    records.push_back(parallel::CompactRuntimeDecompositionRecord{
        .entity_id = weight.entity_id,
        .sfc_key = weight.sfc_key,
        .memory_bytes = weight.memory_bytes,
        .weighted_load = weight.weighted_load,
        .local_index = particle_index,
        .active_target_count_recent = weight.active_target_count_recent,
        .remote_tree_interactions_recent = weight.remote_tree_interactions_recent,
        .current_owner_rank = weight.current_owner_rank,
        .kind = weight.kind,
    });
  }

  if (state.patches.owning_rank.size() != state.patches.size()) {
    state.patches.owning_rank.assign(
        state.patches.size(), static_cast<std::uint32_t>(std::max(world_rank, 0)));
  }
  for (std::size_t patch_row = 0; patch_row < state.patches.size(); ++patch_row) {
    const std::optional<StartupDecompositionWeight> weight =
        computeInitialDecompositionWeightForPatch(state, decomposition_config, patch_row);
    if (!weight.has_value()) {
      continue;
    }
    records.push_back(parallel::CompactRuntimeDecompositionRecord{
        .entity_id = weight->entity_id,
        .sfc_key = weight->sfc_key,
        .memory_bytes = weight->memory_bytes,
        .weighted_load = weight->weighted_load,
        .local_index = state.particles.size() + included_patch_rows.size(),
        .active_target_count_recent = weight->active_target_count_recent,
        .remote_tree_interactions_recent = weight->remote_tree_interactions_recent,
        .current_owner_rank = weight->current_owner_rank,
        .kind = weight->kind,
    });
    included_patch_rows.push_back(static_cast<std::uint32_t>(patch_row));
  }

  auto plan = parallel::buildMortonSfcDecompositionFromCompact(records, decomposition_config);
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    state.particle_sidecar.owning_rank[particle_index] =
        static_cast<std::uint32_t>(plan.owning_rank_by_item[particle_index]);
  }
  for (std::size_t ordinal = 0; ordinal < included_patch_rows.size(); ++ordinal) {
    const std::uint32_t patch_row = included_patch_rows[ordinal];
    if (patch_row < state.patches.owning_rank.size()) {
      state.patches.owning_rank[patch_row] = static_cast<std::uint32_t>(
          plan.owning_rank_by_item[state.particles.size() + ordinal]);
    }
  }

  // Startup records fold components into weighted_load, so the planner's
  // component span stays empty. Recompute each record's components after cuts
  // are known and accumulate them into the same per-rank lanes the rich path
  // populated, preserving recordDistributedProfiling component counters.
  for (const parallel::CompactRuntimeDecompositionRecord& record : records) {
    const std::size_t rank = static_cast<std::size_t>(plan.owning_rank_by_item[record.local_index]);
    if (record.kind == parallel::DecompositionEntityKind::kParticle) {
      const StartupDecompositionWeight weight = computeInitialDecompositionWeightForParticle(
          state, config, decomposition_config, grids, patch_cell_count, record.local_index);
      parallel::addWorkComponentsToMetrics(plan.metrics, rank, weight.work_components, 1.0);
    } else {
      const std::size_t ordinal = record.local_index - state.particles.size();
      const std::optional<StartupDecompositionWeight> weight =
          computeInitialDecompositionWeightForPatch(
              state, decomposition_config, included_patch_rows[ordinal]);
      if (weight.has_value()) {
        parallel::addWorkComponentsToMetrics(plan.metrics, rank, weight->work_components, 1.0);
      }
    }
  }

  parallel::recordDistributedProfiling(profiler, plan.metrics, 0, 0);
  if (profiler != nullptr) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.work_weighted",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.domain_decomposition",
        .message = "initial domain decomposition used explicit work-weight components",
        .payload = {{"world_size", std::to_string(world_size)},
                    {"item_count", std::to_string(records.size())},
                    {"weighted_imbalance_ratio", std::to_string(plan.metrics.weighted_imbalance_ratio)},
                    {"memory_imbalance_ratio", std::to_string(plan.metrics.memory_imbalance_ratio)}}});
  }
  compactStateToCurrentOwner(state, world_rank);
}




}  // namespace

MigrationBalanceRuntime::MigrationBalanceRuntime(
    const core::SimulationConfig& config,
    const RuntimeServices& services) noexcept
    : m_config(config), m_services(services) {}

void MigrationBalanceRuntime::initializeOwnership(
    core::SimulationState& state) const {
  const std::size_t pm_grid_nx =
      static_cast<std::size_t>(m_config.numerics.treepm_pm_grid_nx);
  const int world_size = m_services.mpi_context.worldSize();
  const int world_rank = m_services.mpi_context.worldRank();
  if (pm_grid_nx == 0U || world_size <= 0) {
    throw std::invalid_argument(
        "initial ownership requires a nonzero PM x extent and world size");
  }
  const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
  for (std::size_t index = 0; index < state.particles.size(); ++index) {
    double wrapped_x = state.particles.position_x_comoving[index];
    if (box_size_x > 0.0) {
      wrapped_x = std::fmod(wrapped_x, box_size_x);
      if (wrapped_x < 0.0) {
        wrapped_x += box_size_x;
      }
      if (wrapped_x >= box_size_x) {
        wrapped_x -= box_size_x;
      }
    }
    std::size_t global_x = 0U;
    if (box_size_x > 0.0) {
      const double scaled =
          (wrapped_x / box_size_x) * static_cast<double>(pm_grid_nx);
      global_x = std::min(
          pm_grid_nx - 1U, static_cast<std::size_t>(scaled));
    }
    state.particle_sidecar.owning_rank[index] =
        static_cast<std::uint32_t>(parallel::pmOwnerRankForGlobalX(
            pm_grid_nx, world_size, global_x));
  }
  applyInitialGravityAwareDecomposition(
      state,
      m_config,
      world_size,
      world_rank,
      &m_services.profiler,
      m_services);
}

parallel::LocalOwnershipIdentitySummary
MigrationBalanceRuntime::reduceIdentity(
    const core::SimulationState& state) const {
  return reduceLocalParticleIdentitySummary(state, m_services.mpi_context);
}

std::vector<parallel::TopDomainLeaf>
MigrationBalanceRuntime::authoritativeTopDomainLeaves(
    const core::SimulationState& state,
    std::uint64_t decomposition_epoch,
    core::MemoryReservation* retained_reservation) const {
  if (retained_reservation == nullptr) {
    throw std::invalid_argument("top-domain seed requires a retained-memory reservation sink");
  }
  const int world_rank = m_services.mpi_context.worldRank();
  core::MemoryReservation seed_reservation;
  std::exception_ptr admission_failure;
  try {
    if (m_services.memory_governor != nullptr) {
      // Shared source-view estimate plus caller-specific seed memory. The
      // total is a documented conservative upper envelope: source storage is
      // built first (its gas construction scratch is released before seed
      // records and leaves are installed), then seed records and leaves
      // coexist until the function returns.
      const parallel::RuntimeDecompositionSourceMemoryEstimate source_estimate =
          parallel::estimateRuntimeDecompositionSourceStorage(
              state, std::span<const std::uint32_t>{});
      const std::size_t entity_count = core::checkedSizeAdd(
          state.particles.size(), state.patches.size(),
          "runtime top-domain seed entity count");
      const std::size_t seed_bytes = core::checkedSizeMultiply(
          entity_count, sizeof(parallel::CompactTopDomainSeedRecord),
          "runtime top-domain seed record bytes");
      const std::size_t source_bytes = core::checkedIntegralNarrow<std::size_t>(
          source_estimate.total_bytes,
          "runtime top-domain source estimate byte width");
      const std::size_t leaf_bytes = core::checkedSizeMultiply(
          3U * 8U, sizeof(parallel::TopDomainLeaf),
          "runtime top-domain leaf installation bytes");
      const std::size_t total_bytes = core::checkedSizeAdd(
          core::checkedSizeAdd(source_bytes, seed_bytes, "runtime top-domain seed bytes"),
          leaf_bytes, "runtime top-domain seed bytes");

      seed_reservation = m_services.memory_governor->reserve(
          core::MemoryClass::kPhaseResident,
          core::checkedIntegralNarrow<std::uint64_t>(
              total_bytes, "runtime top-domain seed byte width"),
          "parallel.decomposition.top_domain_seed");
      seed_reservation.commit();
    }
  } catch (...) {
    admission_failure = std::current_exception();
  }
  FailureCoordinator(m_services).rethrowCollectiveFailure(
      admission_failure, "runtime top-domain seed memory admission");

  std::vector<parallel::TopDomainLeaf> leaves;
  std::exception_ptr preparation_failure;
  try {
    RuntimeDecompositionSourceStorage source_storage(state, world_rank, {});
    const parallel::RuntimeDecompositionSourceView source = source_storage.view();
    leaves = parallel::buildAuthoritativeTopDomainLeavesFromSource(
        source,
        makeWorkflowDecompositionConfig(m_config, m_services.mpi_context.worldSize()),
        world_rank,
        decomposition_epoch);
  } catch (...) {
    preparation_failure = std::current_exception();
  }
  FailureCoordinator(m_services).rethrowCollectiveFailure(
      preparation_failure, "runtime top-domain seed preparation");
  *retained_reservation = std::move(seed_reservation);
  return leaves;
}

bool MigrationBalanceRuntime::rebalance(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const parallel::DecompositionRuntimeMeasurements& measurements,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint64_t> expected_global_particle_ids,
    std::uint64_t step_index) const {
  return applyMeasuredRuntimeRebalancePlan(
      state,
      particle_scheduler,
      gas_cell_scheduler,
      m_config,
      m_services,
      m_services.mpi_context.worldRank(),
      measurements,
      active_particle_indices,
      expected_global_particle_ids,
      &m_services.profiler,
      step_index);
}

}  // namespace cosmosim::workflows::internal
