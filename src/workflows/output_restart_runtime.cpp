#include "cosmosim/workflows/output_restart_runtime.hpp"
#include "cosmosim/physics/baryon_composition.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/runtime_console_reporter.hpp"
#include "cosmosim/workflows/analysis_runtime.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/openmp_runtime.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "workflows/internal/output_verification.hpp"

namespace cosmosim::workflows::internal {
namespace {

[[nodiscard]] std::uint64_t simulationOwnedCapacityBytes(
    const core::SimulationState& state) {
  return core::memoryReportBaselineOwnedBytes(
      core::collectSimulationMemoryReport(state));
}

[[nodiscard]] std::uint64_t snapshotPartitionVerificationPeakBytes() {
  // verifySingleFileScienceSnapshotPartitionHdf5 reads at most 4096 rows at
  // once: Coordinates[3], Velocities[3], Masses, ParticleIDs, plus a compact
  // local-row index lane. Keep one MiB of governed headroom for vector capacity
  // and HDF5-owned call staging without charging a second simulation state.
  return 1U << 20U;
}

[[nodiscard]] std::uint64_t incrementalOutputStagingBytes(
    const workflows::RuntimeServices& services,
    std::uint64_t physical_staging_bytes) {
  if (services.memory_governor == nullptr) {
    return physical_staging_bytes;
  }
  return incrementalMemoryBeyondRetained(
      physical_staging_bytes,
      services.memory_governor->policy().planned_overlap_reserve_bytes);
}

[[nodiscard]] core::MemoryReservation reserveOutputStaging(
    const workflows::RuntimeServices& services,
    std::uint64_t physical_staging_bytes,
    std::string_view owner) {
  if (services.memory_governor == nullptr) {
    return {};
  }
  const std::uint64_t incremental_bytes =
      incrementalOutputStagingBytes(services, physical_staging_bytes);
  core::MemoryReservation reservation = services.memory_governor->reserve(
      core::MemoryClass::kDiagnostic, incremental_bytes, owner);
  reservation.commit();
  return reservation;
}

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}

[[nodiscard]] std::string treePmAssignmentSchemeName(
    core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic: return "cic";
    case core::TreePmAssignmentScheme::kTsc: return "tsc";
  }
  throw std::runtime_error("unhandled TreePM assignment scheme enum value");
}

[[nodiscard]] std::string treePmOpeningCriterionName(
    core::TreePmOpeningCriterion opening_criterion) {
  switch (opening_criterion) {
    case core::TreePmOpeningCriterion::kGeometric: return "geometric";
    case core::TreePmOpeningCriterion::kComDistance: return "com_distance";
    case core::TreePmOpeningCriterion::kRelativeForceError:
      return "relative_force_error";
  }
  throw std::runtime_error("unhandled TreePM opening criterion enum value");
}

[[nodiscard]] std::string pmDecompositionModeName(
    core::PmDecompositionMode mode) {
  switch (mode) {
    case core::PmDecompositionMode::kSlab: return "slab";
    case core::PmDecompositionMode::kPencil: return "pencil";
  }
  throw std::runtime_error("unhandled PM decomposition mode enum value");
}

[[nodiscard]] bool hasSpeciesSpecificSoftening(
    const core::SimulationConfig& config) {
  return config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
}

[[nodiscard]] std::string describeSofteningPolicy(
    const core::SimulationConfig& config) {
  return hasSpeciesSpecificSoftening(config)
      ? "comoving_species"
      : "comoving_fixed";
}

[[nodiscard]] io::StochasticPersistentState buildStochasticPersistentState(
    const core::SimulationConfig& config,
    const core::IntegratorState& integrator_state,
    std::uint32_t rank_local_seed_offset) {
  io::StochasticPersistentState stochastic_state;
  if (config.physics.enable_star_formation && config.physics.sf_stochastic_spawning) {
    stochastic_state.modules.push_back(io::StochasticModulePersistentState{
        .module_name = "star_formation",
        .schema_version = 2,
        .rng_policy =
            "stateless_splitmix64(global_seed,gas_cell_id,global_integration_tick,birth_ordinal,model_schema_version)",
        .random_seed = config.physics.sf_random_seed,
        .rank_local_seed_offset = 0U,
        .last_committed_step_index = integrator_state.step_index,
        .deterministic_from_serialized_inputs = true,
    });
  }
  if (config.physics.enable_feedback && config.physics.fb_variant == core::FeedbackVariant::kStochastic) {
    stochastic_state.modules.push_back(io::StochasticModulePersistentState{
        .module_name = "stellar_feedback",
        .schema_version = 1,
        .rng_policy = "stateless_splitmix64(seed,step_index,star_index)",
        .random_seed = config.physics.fb_random_seed,
        .rank_local_seed_offset = rank_local_seed_offset,
        .last_committed_step_index = integrator_state.step_index,
        .deterministic_from_serialized_inputs = true,
    });
  }
  return stochastic_state;
}

[[nodiscard]] std::array<std::uint64_t, 6> snapshotLocalPartCounts(
    const core::SimulationState& state) {
  std::array<std::uint64_t, 6> counts{};
  counts[0] = static_cast<std::uint64_t>(state.cells.size());
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const auto species = static_cast<core::ParticleSpecies>(state.particleSpeciesTag(i));
    switch (species) {
      case core::ParticleSpecies::kGas:
        break;  // PartType0 authority is CellSoa, not a duplicate gas-particle row.
      case core::ParticleSpecies::kDarkMatter: ++counts[1]; break;
      case core::ParticleSpecies::kTracer: ++counts[3]; break;
      case core::ParticleSpecies::kStar: ++counts[4]; break;
      case core::ParticleSpecies::kBlackHole: ++counts[5]; break;
      default: throw std::runtime_error("snapshot output encountered unsupported particle species");
    }
  }
  return counts;
}

[[nodiscard]] std::string formatThreeDigitIndex(std::uint64_t index) {
  std::ostringstream out;
  out << std::setw(3) << std::setfill('0') << index;
  return out.str();
}

[[nodiscard]] std::string snapshotGenerationId(
    const core::FrozenConfig& frozen_config,
    std::uint64_t step_index) {
  std::string hash = frozen_config.provenance.config_hash_hex;
  if (hash.size() > 12U) hash.resize(12U);
  return "snap_" + formatThreeDigitIndex(step_index) + "_" + hash;
}

[[nodiscard]] std::string snapshotMemberFilename(
    std::string_view stem,
    std::uint64_t step_index,
    int world_size,
    int world_rank) {
  std::ostringstream out;
  out << stem << '_' << formatThreeDigitIndex(step_index);
  if (world_size > 1) out << '.' << world_rank;
  out << ".hdf5";
  return out.str();
}

[[nodiscard]] bool useSingleFileScienceSnapshot(
    const core::SimulationConfig& config,
    const parallel::DistributedExecutionTopology& topology) {
  switch (config.output.snapshot_layout) {
    case core::ScienceSnapshotLayout::kSharded:
      return topology.world_size <= 1;
    case core::ScienceSnapshotLayout::kAggregated:
      throw std::runtime_error(
          "output.snapshot_layout=aggregated reserves the fixed-file-count I/O topology contract but its production backend is not yet qualified");
    case core::ScienceSnapshotLayout::kSingle:
    case core::ScienceSnapshotLayout::kAuto:
      if (topology.world_size <= 1) return true;
#if COSMOSIM_ENABLE_MPI && COSMOSIM_HDF5_PARALLEL
      return true;
#else
      throw std::runtime_error(
          "analysis-ready single-file MPI science snapshots require an MPI-enabled Parallel HDF5 build; select a Parallel-HDF5 build or explicitly set output.snapshot_layout=sharded");
#endif
  }
  throw std::logic_error("unhandled science snapshot layout");
}

[[nodiscard]] std::string singleSnapshotFilename(
    std::string_view stem, std::uint64_t step_index) {
  return std::string(stem) + "_" + formatThreeDigitIndex(step_index) + ".hdf5";
}

[[nodiscard]] std::string formatIndexedRankedFileStem(
    std::string_view stem,
    std::uint64_t index,
    int world_size,
    int world_rank) {
  std::ostringstream out;
  out << stem << '_' << std::setw(3) << std::setfill('0') << index;
  if (world_size > 1) {
    out << "_rank" << std::setw(3) << std::setfill('0') << world_rank;
  }
  out << ".hdf5";
  return out.str();
}


[[nodiscard]] core::ProvenanceRecord makeGravityAwareProvenanceRecord(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config) {
  core::ProvenanceRecord record = core::makeProvenanceRecord(
      frozen_config.provenance.config_hash_hex,
      "unknown",
      0,
      frozen_config.normalized_text);
  record.config_schema_name = "cosmosim_config";
  record.config_schema_version = std::to_string(frozen_config.config.schema_version);
  record.raw_input_config = frozen_config.raw_text;
  record.normalized_config = frozen_config.normalized_text;
  record.derived_runtime_state =
      core::serializeDerivedRuntimeConfig(core::deriveRuntimeConfig(frozen_config));
  record.derived_runtime_state += "source_identity=" + std::string(COSMOSIM_SOURCE_IDENTITY) + "\n";
  record.derived_runtime_state += "source_manifest_sha256=" + std::string(COSMOSIM_SOURCE_MANIFEST_SHA256) + "\n";
  record.derived_runtime_state += "dependency_providers=" + std::string(COSMOSIM_DEPENDENCY_PROVIDERS) + "\n";
  const core::OpenMpRuntimeInfo openmp = core::openMpRuntimeInfo();
  record.derived_runtime_state += "openmp_compiled=" + std::string(openmp.compiled ? "true" : "false") + "\n";
  record.derived_runtime_state += "openmp_requested_threads=" + std::to_string(openmp.requested_threads) + "\n";
  record.derived_runtime_state += "openmp_active_threads=" + std::to_string(openmp.configured_threads) + "\n";
  const double dx = config.cosmology.box_size_x_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nx);
  const double dy = config.cosmology.box_size_y_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_ny);
  const double dz = config.cosmology.box_size_z_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nz);
  const double mesh_spacing_mpc_comoving = std::cbrt(dx * dy * dz);
  const gravity::TreePmSplitPolicy split_policy = gravity::makeTreePmSplitPolicyFromMeshSpacing(
      config.numerics.treepm_asmth_cells,
      config.numerics.treepm_rcut_cells,
      mesh_spacing_mpc_comoving);
  record.gravity_treepm_pm_grid = config.numerics.treepm_pm_grid_nx;
  record.gravity_treepm_pm_grid_nx = config.numerics.treepm_pm_grid_nx;
  record.gravity_treepm_pm_grid_ny = config.numerics.treepm_pm_grid_ny;
  record.gravity_treepm_pm_grid_nz = config.numerics.treepm_pm_grid_nz;
  record.gravity_treepm_assignment_scheme =
      treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
  record.gravity_treepm_window_deconvolution =
      config.numerics.treepm_enable_window_deconvolution;
  record.gravity_treepm_asmth_cells = config.numerics.treepm_asmth_cells;
  record.gravity_treepm_rcut_cells = config.numerics.treepm_rcut_cells;
  record.gravity_treepm_tree_opening_criterion =
      treePmOpeningCriterionName(config.numerics.treepm_tree_opening_criterion);
  record.gravity_treepm_tree_opening_theta = config.numerics.treepm_tree_opening_theta;
  record.gravity_treepm_tree_relative_force_tolerance =
      config.numerics.treepm_tree_relative_force_tolerance;
  record.gravity_treepm_tree_relative_force_acceleration_floor =
      config.numerics.treepm_tree_relative_force_acceleration_floor;
  record.gravity_treepm_mesh_spacing_mpc_comoving = mesh_spacing_mpc_comoving;
  record.gravity_treepm_mesh_spacing_x_mpc_comoving = dx;
  record.gravity_treepm_mesh_spacing_y_mpc_comoving = dy;
  record.gravity_treepm_mesh_spacing_z_mpc_comoving = dz;
  record.gravity_treepm_split_scale_mpc_comoving = split_policy.split_scale_comoving;
  record.gravity_treepm_cutoff_radius_mpc_comoving = split_policy.cutoff_radius_comoving;
  record.gravity_treepm_update_cadence_steps = config.numerics.treepm_update_cadence_steps;
  record.gravity_treepm_pm_decomposition_mode =
      pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
  record.gravity_treepm_tree_exchange_batch_bytes =
      config.numerics.treepm_tree_exchange_batch_bytes;
  record.gravity_treepm_pm_exchange_batch_bytes =
      config.numerics.treepm_pm_exchange_batch_bytes;
  record.gravity_softening_policy = describeSofteningPolicy(config);
  record.gravity_softening_kernel = "plummer";
  record.gravity_softening_epsilon_kpc_comoving = config.numerics.gravity_softening_kpc_comoving;
  record.gravity_pm_fft_backend = gravity::PmSolver::fftBackendName();
  record.gravity_pm_backend_capability = std::string(
      gravity::pmBackendCapabilityName(gravity::pmBackendCapability()));
  record.gravity_acceptance_profile_id = "unverified_current_source";
  switch (config.mode.zoom_long_range_strategy) {
    case core::ZoomLongRangeStrategy::kDisabled:
      record.zoom_long_range_strategy = "disabled";
      break;
    case core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection:
      record.zoom_long_range_strategy = "global_coarse_plus_focused_highres_correction";
      break;
  }
  record.zoom_region_center_x_mpc_comoving = config.mode.zoom_region_center_x_mpc_comoving;
  record.zoom_region_center_y_mpc_comoving = config.mode.zoom_region_center_y_mpc_comoving;
  record.zoom_region_center_z_mpc_comoving = config.mode.zoom_region_center_z_mpc_comoving;
  record.zoom_region_radius_mpc_comoving = config.mode.zoom_region_radius_mpc_comoving;
  record.zoom_focused_pm_grid = std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_nz);
  record.zoom_contamination_radius_mpc_comoving = config.mode.zoom_contamination_radius_mpc_comoving;
  record.deterministic_mode = config.parallel.deterministic_reduction ? "deterministic_reduction" : "non_deterministic_reduction";
  return record;
}

void applyExecutionTopologyToProvenance(
    core::ProvenanceRecord* record,
    const parallel::DistributedExecutionTopology& topology) {
  if (record == nullptr) {
    throw std::invalid_argument("applyExecutionTopologyToProvenance requires a non-null record");
  }
  record->author_rank = topology.world_rank;
  record->mpi_world_size = topology.world_size;
  record->mpi_node_local_rank = topology.local_rank;
  record->mpi_summary = topology.mpi_enabled ? "enabled" : "disabled";
  if (topology.device_assignment.uses_cuda) {
    record->gpu_summary = "cuda_device=" + std::to_string(topology.device_assignment.assigned_device_index) +
        ",visible_devices=" + std::to_string(topology.device_assignment.visible_device_count);
  } else {
    record->gpu_summary = "disabled_or_unavailable";
  }
}

[[nodiscard]] std::string pmSlabSignature(const parallel::DistributedRestartState& distributed_state) {
  std::ostringstream stream;
  for (std::size_t rank = 0; rank < distributed_state.pm_slab_begin_x_by_rank.size(); ++rank) {
    if (rank > 0) {
      stream << ';';
    }
    stream << rank << ':' << distributed_state.pm_slab_begin_x_by_rank[rank] << '-' <<
        distributed_state.pm_slab_end_x_by_rank[rank];
  }
  return stream.str();
}

bool maybeWriteOutputs(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config,
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const workflows::GravityRestartStateProvider& gravity_state,
    const workflows::AnalysisRuntime& analysis_state,
    const workflows::RuntimeServices& services,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    bool write_outputs_enabled,
    bool snapshot_due,
    bool checkpoint_due,
    std::uint64_t snapshot_interval_steps,
    std::uint64_t next_snapshot_step_index,
    double snapshot_interval_time_code,
    double next_snapshot_time_code,
    double restart_resume_dt_time_code) {
  if (!write_outputs_enabled || (!snapshot_due && !checkpoint_due)) {
    return false;
  }

#if !COSMOSIM_ENABLE_HDF5
  throw std::runtime_error(
      "runtime outputs requested, but this build lacks HDF5 support. Reconfigure with COSMOSIM_ENABLE_HDF5=ON.");
#else
  if (integrator_state.step_index == 0) {
    return false;
  }
  if (checkpoint_due) {
    core::assertCanWriteCheckpointAtBoundary(integrator_state, scheduler.currentTick());
    if (gas_cell_scheduler.currentTick() != scheduler.currentTick()) {
      throw std::runtime_error("output checkpoint requires particle and gas-cell schedulers on one integer tick");
    }
  }
  if (snapshot_due) {
    core::assertCanWriteSnapshotAtBoundary(integrator_state);
  }
  if ((snapshot_due || checkpoint_due) &&
      !core::isOutputSafeBoundary(integrator_state.last_completed_boundary_kind)) {
    return false;
  }

  bool output_flushed = false;
  if (snapshot_due) {
    const auto topology = gravity_state.runtimeTopology();
    const std::array<std::uint64_t, 6> local_counts = snapshotLocalPartCounts(state);
    std::array<std::uint64_t, 6> global_counts{};
    for (std::size_t i = 0; i < global_counts.size(); ++i) {
      global_counts[i] = services.mpi_context.allreduceSumUint64(local_counts[i]);
    }
    const bool single_file_snapshot = useSingleFileScienceSnapshot(config, topology);
    const bool collective_single_file = single_file_snapshot && topology.world_size > 1;
    const std::uint32_t physical_file_count = single_file_snapshot
        ? 1U
        : core::checkedIntegralNarrow<std::uint32_t>(topology.world_size, "snapshot member count");
    const std::filesystem::path shared_run_directory =
        report.shared_run_directory.empty() ? report.run_directory : report.shared_run_directory;
    if (collective_single_file) {
      std::string local_shared_run_directory;
      std::exception_ptr path_preparation_failure;
      try {
        std::error_code absolute_error;
        std::filesystem::path absolute_shared_run_directory =
            std::filesystem::absolute(shared_run_directory, absolute_error);
        if (absolute_error) {
          throw std::runtime_error(
              "failed to resolve collective science snapshot directory '" +
              shared_run_directory.string() + "': " + absolute_error.message());
        }
        local_shared_run_directory =
            absolute_shared_run_directory.lexically_normal().generic_string();
      } catch (...) {
        path_preparation_failure = std::current_exception();
      }
      FailureCoordinator(services).rethrowCollectiveFailure(
          path_preparation_failure,
          "collective science snapshot shared-path preparation");

      std::vector<std::uint8_t> root_path_bytes;
      if (services.mpi_context.isRoot()) {
        root_path_bytes.assign(
            local_shared_run_directory.begin(), local_shared_run_directory.end());
      }
      const std::vector<std::uint8_t> agreed_path_bytes =
          services.mpi_context.broadcastBytesFromRoot(root_path_bytes, 0);
      const std::string agreed_shared_run_directory(
          agreed_path_bytes.begin(), agreed_path_bytes.end());
      std::exception_ptr path_consensus_failure;
      if (local_shared_run_directory != agreed_shared_run_directory) {
        path_consensus_failure = std::make_exception_ptr(std::runtime_error(
            "collective single-file science snapshot requires one identical shared run "
            "directory on every MPI rank; root='" + agreed_shared_run_directory +
            "' local='" + local_shared_run_directory + "'"));
      }
      FailureCoordinator(services).rethrowCollectiveFailure(
          path_consensus_failure,
          "collective science snapshot shared-path consensus");
    }
    const std::filesystem::path snapshot_directory = shared_run_directory / "snapshots";
    const std::filesystem::path final_single_path = snapshot_directory /
        singleSnapshotFilename(config.output.output_stem, integrator_state.step_index);
    const std::filesystem::path partial_single_path =
        std::filesystem::path(final_single_path.string() + ".partial");
    report.snapshot_path = single_file_snapshot
        ? final_single_path
        : snapshot_directory / snapshotMemberFilename(
              config.output.output_stem, integrator_state.step_index,
              topology.world_size, topology.world_rank);
    report.snapshot_set_path = snapshot_directory /
        (config.output.output_stem + "_" +
         formatThreeDigitIndex(integrator_state.step_index) + ".complete");
    const std::string generation_id =
        snapshotGenerationId(frozen_config, integrator_state.step_index);

    std::exception_ptr snapshot_replacement_failure;
    if (services.mpi_context.isRoot()) {
      try {
        std::error_code directory_error;
        std::filesystem::create_directories(
            snapshot_directory, directory_error);
        if (directory_error) {
          throw std::runtime_error(
              "failed to establish science snapshot directory '" +
              snapshot_directory.string() + "': " +
              directory_error.message());
        }
        const bool snapshot_path_is_directory =
            std::filesystem::is_directory(
                snapshot_directory, directory_error);
        if (directory_error) {
          throw std::runtime_error(
              "failed to validate science snapshot directory '" +
              snapshot_directory.string() + "': " +
              directory_error.message());
        }
        if (!snapshot_path_is_directory) {
          throw std::runtime_error(
              "science snapshot path exists but is not a directory: " +
              snapshot_directory.string());
        }
        if (std::filesystem::exists(report.snapshot_set_path) ||
            (single_file_snapshot && std::filesystem::exists(final_single_path))) {
          throw std::runtime_error(
              "refusing to overwrite a committed science snapshot: " +
              report.snapshot_path.string());
        }
        if (collective_single_file) {
          std::error_code ignored;
          std::filesystem::remove(partial_single_path, ignored);
        }
      } catch (...) {
        snapshot_replacement_failure = std::current_exception();
      }
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        snapshot_replacement_failure, "science snapshot replacement preflight");

    core::ProvenanceRecord snapshot_provenance =
        makeGravityAwareProvenanceRecord(frozen_config, config);
    snapshot_provenance.derived_runtime_state += report.restart_config_compatibility_record;
    applyExecutionTopologyToProvenance(&snapshot_provenance, topology);
    snapshot_provenance.gravity_treepm_decomposition_epoch =
        gravity_state.decompositionEpoch();
    if (collective_single_file) {
      std::vector<std::uint8_t> serialized;
      if (services.mpi_context.isRoot()) {
        const std::string text = core::serializeProvenanceRecord(snapshot_provenance);
        serialized.assign(text.begin(), text.end());
      }
      const std::vector<std::uint8_t> root_bytes =
          services.mpi_context.broadcastBytesFromRoot(serialized, 0);
      snapshot_provenance = core::deserializeProvenanceRecord(
          std::string_view(reinterpret_cast<const char*>(root_bytes.data()), root_bytes.size()));
    }

    io::SnapshotWritePayload snapshot_payload;
    snapshot_payload.state = &state;
    snapshot_payload.memory_governor = services.memory_governor;
    snapshot_payload.config = &config;
    snapshot_payload.normalized_config_text = frozen_config.normalized_text;
    snapshot_payload.provenance = snapshot_provenance;
    snapshot_payload.set_member.member_index = single_file_snapshot
        ? 0U
        : core::checkedIntegralNarrow<std::uint32_t>(topology.world_rank, "snapshot member index");
    snapshot_payload.set_member.num_files_per_snapshot = physical_file_count;
    snapshot_payload.set_member.global_part_count = global_counts;
    snapshot_payload.set_member.has_global_part_count = true;
    snapshot_payload.set_member.generation_id = generation_id;
    snapshot_payload.collective_single_file = collective_single_file;
    if (collective_single_file) {
      for (std::size_t type = 0; type < 6U; ++type) {
        snapshot_payload.collective_file_row_offset[type] =
            services.mpi_context.exclusiveScanSumUint64(local_counts[type]);
      }
      snapshot_payload.collective_write_particle_softening =
          services.mpi_context.allreduceSumUint64(
              state.particle_sidecar.gravity_softening_comoving.empty() ? 0U : 1U) > 0U;
      snapshot_payload.collective_write_particle_softening_override =
          services.mpi_context.allreduceSumUint64(
              state.particle_sidecar.has_gravity_softening_override.empty() ? 0U : 1U) > 0U;
    }

    io::SnapshotIoPolicy snapshot_policy;
    snapshot_policy.dialect = config.units.coordinate_frame == core::CoordinateFrame::kComoving
        ? io::SnapshotDialect::kArepoFormat3
        : io::SnapshotDialect::kChuiNative;
    snapshot_policy.durable_publication = false;
    snapshot_policy.caller_managed_publication = collective_single_file;

    std::exception_ptr local_snapshot_write_failure;
    try {
      io::writeScienceSnapshotHdf5(
          collective_single_file ? partial_single_path : report.snapshot_path,
          snapshot_payload, snapshot_policy);
    } catch (...) {
      local_snapshot_write_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        local_snapshot_write_failure, "science snapshot write");

    // Structural verification (non-collective members only) reads header
    // metadata, so it is charged like the bounded partition readback rather
    // than a second materialized simulation state.
    const bool structural_member_verification = !collective_single_file &&
        config.output.snapshot_verification == core::SnapshotVerificationPolicy::kStructural;
    core::MemoryReservation snapshot_readback_reservation;
    std::exception_ptr snapshot_readback_admission_failure;
    try {
      snapshot_readback_reservation = reserveOutputStaging(
          services,
          (collective_single_file || structural_member_verification)
              ? snapshotPartitionVerificationPeakBytes()
              : simulationOwnedCapacityBytes(state),
          collective_single_file ? "io.snapshot.partition_readback"
                                 : "io.snapshot.readback");
    } catch (...) {
      snapshot_readback_admission_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        snapshot_readback_admission_failure, "science snapshot readback memory admission");

    const auto snapshot_verification_begin = std::chrono::steady_clock::now();
    std::exception_ptr local_snapshot_readback_failure;
    std::string local_verification_detail;
    try {
      report.snapshot_roundtrip_executed = true;
      if (collective_single_file) {
        io::verifySingleFileScienceSnapshotPartitionHdf5(
            partial_single_path, snapshot_payload, snapshot_policy);
        report.snapshot_roundtrip_ok = true;
        local_verification_detail = "distributed_rank_partition_exact_readback";
      } else if (structural_member_verification) {
        // Schema, file kind, header cosmology/time, provenance hash and
        // per-type counts. The SHA-256 member sidecar and root completion
        // marker/independent set validation below remain mandatory.
        const io::SnapshotMemberInspection member =
            io::inspectSnapshotMemberHdf5(report.snapshot_path);
        const auto expected_scale_factor = state.metadata.scale_factor;
        const auto equal = [](double lhs, double rhs) {
          const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
          return std::isfinite(lhs) && std::isfinite(rhs) &&
              std::abs(lhs - rhs) <= 32.0 * std::numeric_limits<double>::epsilon() * scale;
        };
        std::string mismatch;
        if (member.schema_name != io::gadgetArepoSchemaMap().schema_name ||
            member.schema_version != io::gadgetArepoSchemaMap().schema_version) {
          mismatch = "schema identity";
        } else if (member.file_kind != io::sharedIoContractNames().science_snapshot_file_kind) {
          mismatch = "file kind";
        } else if (member.config_hash_hex != snapshot_provenance.config_hash_hex) {
          mismatch = "provenance config hash";
        } else if (member.generation_id != generation_id) {
          mismatch = "generation id";
        } else if (member.local_part_count != local_counts ||
                   member.global_part_count != global_counts ||
                   member.num_files_per_snapshot != physical_file_count ||
                   member.member_index != snapshot_payload.set_member.member_index) {
          mismatch = "member/global particle counts or member index";
        } else if (!equal(member.scale_factor, expected_scale_factor) ||
                   !equal(member.box_size_x, config.cosmology.box_size_x_mpc_comoving) ||
                   !equal(member.box_size_y, config.cosmology.box_size_y_mpc_comoving) ||
                   !equal(member.box_size_z, config.cosmology.box_size_z_mpc_comoving) ||
                   !equal(member.omega_matter, config.cosmology.omega_matter) ||
                   !equal(member.omega_lambda, config.cosmology.omega_lambda) ||
                   !equal(member.hubble_param, config.cosmology.hubble_param)) {
          mismatch = "Header cosmology/time metadata";
        }
        report.snapshot_roundtrip_ok = mismatch.empty();
        local_verification_detail = "member_structural_header_verification";
        if (!report.snapshot_roundtrip_ok) {
          throw std::runtime_error(
              "snapshot structural verification failed: " + mismatch);
        }
      } else {
        io::SnapshotReadOptions local_read_options;
        local_read_options.require_complete_chui_set = false;
        const io::SnapshotReadResult snapshot_read =
            io::readGadgetArepoSnapshotHdf5(report.snapshot_path, config, local_read_options);
        const SnapshotRoundtripVerification snapshot_verification =
            verifySnapshotRoundtrip(
                snapshot_read, state, config, frozen_config.normalized_text,
                snapshot_provenance);
        report.snapshot_roundtrip_ok = snapshot_verification.ok;
        local_verification_detail = snapshot_verification.detail;
        if (!report.snapshot_roundtrip_ok) {
          throw std::runtime_error(
              "snapshot scientific roundtrip verification failed: " +
              snapshot_verification.detail);
        }
      }
    } catch (...) {
      local_snapshot_readback_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        local_snapshot_readback_failure, "science snapshot readback");
    profiler.counters().addCount("io.snapshot.verification_us", static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - snapshot_verification_begin).count()));
    profiler.counters().addCount(structural_member_verification
        ? "io.snapshot.structural_verifications" : "io.snapshot.full_verifications");

    std::exception_ptr completion_failure;
    if (services.mpi_context.isRoot()) {
      bool single_file_published = false;
      try {
        if (collective_single_file) {
          io::publishSingleFileScienceSnapshot(
              partial_single_path, final_single_path,
              snapshot_payload.set_member, false);
          single_file_published = true;
          io::writeSnapshotSetCompletionMarker(
              final_single_path, generation_id, 1U, global_counts, false,
              io::SnapshotCompletionIntegrityMode::kDistributedScienceReadback);
        } else {
          io::writeSnapshotSetCompletionMarker(
              report.snapshot_path, generation_id, physical_file_count,
              global_counts, false);
        }
        io::SnapshotReadOptions validation_options;
        if (collective_single_file) {
          validation_options.validate_global_id_uniqueness = false;
        }
        const io::SnapshotValidationReport validation =
            io::validateSnapshotSetHdf5(report.snapshot_set_path, validation_options);
        validation.requireValid();
        if (validation.inspection.global_part_count != global_counts ||
            validation.inspection.num_files_per_snapshot != physical_file_count) {
          throw std::runtime_error(
              "independent snapshot-set validation disagrees with runtime global metadata");
        }
      } catch (...) {
        completion_failure = std::current_exception();
        if (collective_single_file && single_file_published) {
          // A failed completion/independent-validation phase must not leave a
          // final-named science product that the workflow itself rejected.
          std::error_code ignored;
          std::filesystem::remove(report.snapshot_set_path, ignored);
          std::filesystem::remove(final_single_path, ignored);
        }
      }
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        completion_failure, "science snapshot completion/publication");

    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "snapshot.write.complete",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "io.snapshot",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = structural_member_verification
            ? "science snapshot written, structurally verified, checksummed, and transactionally published"
            : "science snapshot written, scientifically read back, and transactionally published",
        .payload = {{"science_path", report.snapshot_path.string()},
                    {"set_path", report.snapshot_set_path.string()},
                    {"generation_id", generation_id},
                    {"layout", single_file_snapshot ? "single" : "sharded"},
                    {"verification", local_verification_detail}},
    });
    if (services.console_reporter != nullptr) {
      services.console_reporter->emitSnapshotCommitted(
          integrator_state.step_index, report.snapshot_path,
          report.snapshot_set_path, physical_file_count);
    }
    profiler.counters().addCount("io.snapshot.writes");
    output_flushed = true;
  }

  if (checkpoint_due) {
    // M48-08: write-time verification owns only its bounded diagnostic
    // workspace and reserves it inside verifyRestartCheckpointHdf5(). Do not
    // reserve a second canonical state/scheduler/force-cache population here.

    core::IntegratorState restart_integrator_state = integrator_state;
    if (restart_resume_dt_time_code > 0.0) {
      // The completed step may have been shortened to land on an output
      // event.  Persist the pre-clip proposal as the next-step interval so a
      // resumed workflow follows the same timeline as an uninterrupted run.
      restart_integrator_state.dt_time_code = restart_resume_dt_time_code;
    }
    io::RestartWritePayload restart_payload;
    restart_payload.persistent_state.simulation_state = &state;
    restart_payload.integrator_state = &restart_integrator_state;
    restart_payload.scheduler = &scheduler;
    restart_payload.gas_cell_scheduler = &gas_cell_scheduler;
    const io::GravityForceCachePersistentView gravity_force_cache =
        gravity_state.restartForceCacheView(state);
    restart_payload.gravity_force_cache_view = &gravity_force_cache;
    restart_payload.provenance =
        makeGravityAwareProvenanceRecord(frozen_config, config);
    restart_payload.provenance.derived_runtime_state += report.restart_config_compatibility_record;
    // Optional science cadence is explicitly best-effort and nonpersistent.
    // Persist the pending/missed summary as provenance, not restart truth.
    // A resumed run starts fresh; historical physical states are never
    // fabricated by replaying a later state under an earlier epoch.
    restart_payload.provenance.derived_runtime_state +=
        analysis_state.optionalCadenceProvenance();
    applyExecutionTopologyToProvenance(
        &restart_payload.provenance, gravity_state.runtimeTopology());
    restart_payload.normalized_config_text = frozen_config.normalized_text;
    restart_payload.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
    restart_payload.distributed_gravity_state.schema_version = 2;
    restart_payload.distributed_gravity_state.decomposition_epoch = gravity_state.decompositionEpoch();
    restart_payload.distributed_gravity_state.world_size = gravity_state.runtimeTopology().world_size;
    restart_payload.distributed_gravity_state.pm_grid_nx = gravity_state.pmGridShape().nx;
    restart_payload.distributed_gravity_state.pm_grid_ny = gravity_state.pmGridShape().ny;
    restart_payload.distributed_gravity_state.pm_grid_nz = gravity_state.pmGridShape().nz;
    restart_payload.distributed_gravity_state.pm_decomposition_mode =
        pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
    restart_payload.distributed_gravity_state.gravity_kick_opportunity =
        integrator_state.pm_sync_state.gravityKickOpportunity();
    restart_payload.distributed_gravity_state.pm_update_cadence_steps =
        integrator_state.pm_sync_state.cadenceSteps();
    restart_payload.distributed_gravity_state.long_range_field_version =
        integrator_state.pm_sync_state.fieldVersion();
    restart_payload.distributed_gravity_state.last_long_range_refresh_opportunity =
        integrator_state.pm_sync_state.lastRefreshOpportunity();
    restart_payload.distributed_gravity_state.long_range_field_built_step_index =
        integrator_state.pm_sync_state.lastRefreshStepIndex();
    restart_payload.distributed_gravity_state.long_range_field_built_scale_factor =
        integrator_state.pm_sync_state.lastRefreshScaleFactor();
    restart_payload.distributed_gravity_state.long_range_restart_policy = "deterministic_rebuild";
    restart_payload.output_cadence_state.output_enabled = write_outputs_enabled;
    restart_payload.output_cadence_state.write_restarts = config.output.write_restarts;
    restart_payload.output_cadence_state.snapshot_due = false;
    restart_payload.output_cadence_state.checkpoint_due = false;
    restart_payload.output_cadence_state.last_completed_step_index = integrator_state.step_index;
    restart_payload.output_cadence_state.snapshot_interval_steps = snapshot_interval_steps;
    restart_payload.output_cadence_state.next_snapshot_step_index = next_snapshot_step_index;
    restart_payload.output_cadence_state.snapshot_interval_time_code =
        snapshot_interval_time_code;
    restart_payload.output_cadence_state.next_snapshot_time_code =
        next_snapshot_time_code;
    restart_payload.output_cadence_state.snapshot_stem = config.output.output_stem;
    restart_payload.output_cadence_state.restart_stem = config.output.restart_stem;
    restart_payload.stochastic_state = buildStochasticPersistentState(
        config,
        integrator_state,
        static_cast<std::uint32_t>(std::max(gravity_state.runtimeTopology().world_rank, 0)));
    restart_payload.distributed_gravity_state.borrowOwningRanksForWrite(state);
    const std::size_t world_size = static_cast<std::size_t>(gravity_state.runtimeTopology().world_size);
    restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank.resize(world_size, 0);
    restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank.resize(world_size, 0);
    for (std::size_t rank = 0; rank < world_size; ++rank) {
      const parallel::PmSlabRange range = parallel::pmOwnedXRangeForRank(
          gravity_state.pmGridShape().nx,
          static_cast<int>(world_size),
          static_cast<int>(rank));
      restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank[rank] = range.begin_x;
      restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank[rank] = range.end_x;
    }
    restart_payload.provenance.gravity_treepm_decomposition_epoch =
        restart_payload.distributed_gravity_state.decomposition_epoch;
    restart_payload.provenance.gravity_treepm_restart_world_size =
        restart_payload.distributed_gravity_state.world_size;
    restart_payload.provenance.gravity_treepm_restart_pm_grid =
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nx) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_ny) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nz);
    restart_payload.provenance.gravity_treepm_restart_slab_signature =
        pmSlabSignature(restart_payload.distributed_gravity_state);
    restart_payload.provenance.gravity_treepm_restart_kick_opportunity =
        restart_payload.distributed_gravity_state.gravity_kick_opportunity;
    restart_payload.provenance.gravity_treepm_restart_field_version =
        restart_payload.distributed_gravity_state.long_range_field_version;
    restart_payload.provenance.gravity_treepm_long_range_restart_policy =
        restart_payload.distributed_gravity_state.long_range_restart_policy;

    report.restart_path = report.run_directory / formatIndexedRankedFileStem(config.output.restart_stem, integrator_state.step_index, gravity_state.runtimeTopology().world_size, gravity_state.runtimeTopology().world_rank);
    std::optional<io::RestartVerificationResult> restart_verification_result;
    std::exception_ptr restart_io_failure;
    try {
      io::writeRestartCheckpointHdf5(
          report.restart_path, restart_payload,
          io::RestartWritePolicy{.memory_governor = services.memory_governor});
      report.restart_roundtrip_executed = true;
      restart_verification_result.emplace(io::verifyRestartCheckpointHdf5(
          report.restart_path, restart_payload,
          io::RestartVerificationPolicy{.memory_governor = services.memory_governor}));
    } catch (...) {
      restart_io_failure = std::current_exception();
    }
    FailureCoordinator(services).rethrowCollectiveFailure(
        restart_io_failure, "restart checkpoint write/streaming verification");
    const io::RestartVerificationResult& restart_verification =
        *restart_verification_result;
    const auto compatibility = parallel::evaluateDistributedRestartCompatibility(
        restart_payload.distributed_gravity_state,
        gravity_state.runtimeTopology());
    const bool restart_rank_qualified_name =
        gravity_state.runtimeTopology().world_size == 1 ||
        report.restart_path.filename().string().find("_rank") != std::string::npos;
    report.restart_roundtrip_ok =
        restart_verification.ok && restart_rank_qualified_name && compatibility.compatible();
    // A distributed checkpoint is only verified when every rank's local
    // streaming field verification succeeds. This collective belongs to the
    // checkpoint verification boundary, not to console presentation.
    const std::uint64_t restart_verification_failed_ranks =
        FailureCoordinator(services).failedRankCount(!report.restart_roundtrip_ok);
    report.restart_roundtrip_ok = restart_verification_failed_ranks == 0U;
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.write.complete",
        .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.restart",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "restart checkpoint written and streaming-verified",
        .payload = {{"path", report.restart_path.string()},
                    {"restart_schema", io::restartSchema().name},
                    {"restart_schema_version", std::to_string(io::restartSchema().version)},
                    {"boundary_kind", std::string(core::stepBoundaryKindName(integrator_state.last_completed_boundary_kind))},
                    {"scheduler_current_tick", std::to_string(scheduler.currentTick())},
                    {"scheduler_max_bin", std::to_string(scheduler.maxBin())},
                    {"pm_cadence_steps", std::to_string(integrator_state.pm_sync_state.cadenceSteps())},
                    {"pm_gravity_kick_opportunity", std::to_string(integrator_state.pm_sync_state.gravityKickOpportunity())},
                    {"pm_field_version", std::to_string(integrator_state.pm_sync_state.fieldVersion())},
                    {"pm_long_range_field_valid", integrator_state.pm_long_range_field_valid ? "true" : "false"},
                    {"restart_verification_mode", "streaming_field_exact"},
                    {"verification_workspace_high_water_bytes", std::to_string(restart_verification.workspace_high_water_bytes)},
                    {"verification_bytes_read", std::to_string(restart_verification.verified_bytes)},
                    {"verification_dataset_count", std::to_string(restart_verification.verified_dataset_count)},
                    {"verification_failed_field", restart_verification.failed_field},
                    {"output_next_snapshot_step_index", std::to_string(restart_payload.output_cadence_state.next_snapshot_step_index)},
                    {"output_next_snapshot_time_code", formatRuntimeDouble(
                        restart_payload.output_cadence_state.next_snapshot_time_code)},
                    {"stochastic_module_count", std::to_string(restart_payload.stochastic_state.modules.size())},
                    {"verification_failed_ranks", std::to_string(restart_verification_failed_ranks)}},
    });
    if (services.console_reporter != nullptr) {
      if (report.restart_roundtrip_ok) {
        services.console_reporter->emitRestartCommitted(
            integrator_state.step_index, report.restart_path);
      } else {
        services.console_reporter->emitWarning(
            "io.restart", "restart checkpoint streaming verification failed");
      }
    }
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.verify.complete",
        .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.restart",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "restart checkpoint streaming verification complete",
        .payload = {{"path", report.restart_path.string()},
                    {"restart_schema", restart_verification.schema_name},
                    {"restart_schema_version", std::to_string(restart_verification.schema_version)},
                    {"boundary_kind", restart_verification.diagnostics.last_completed_boundary_kind},
                    {"restart_safe", restart_verification.diagnostics.restart_safe ? "true" : "false"},
                    {"scheduler_current_tick", std::to_string(restart_verification.diagnostics.scheduler_current_tick)},
                    {"scheduler_max_bin", std::to_string(restart_verification.diagnostics.scheduler_max_bin)},
                    {"scheduler_element_count", std::to_string(restart_verification.diagnostics.scheduler_element_count)},
                    {"scheduler_pending_transition_count", std::to_string(restart_verification.diagnostics.scheduler_pending_transition_count)},
                    {"pm_field_version", std::to_string(restart_verification.diagnostics.pm_field_version)},
                    {"payload_hash_hex", restart_verification.payload_sha256_hex},
                    {"verification_workspace_high_water_bytes", std::to_string(restart_verification.workspace_high_water_bytes)},
                    {"verification_bytes_read", std::to_string(restart_verification.verified_bytes)}},
    });
    profiler.counters().addCount("io.restart.writes");
    output_flushed = true;
  }
  return output_flushed;
#endif
}

[[nodiscard]] std::string joinRestartCompatibilityMessages(
    std::span<const std::string> mismatch_messages) {
  std::ostringstream stream;
  for (std::size_t i = 0; i < mismatch_messages.size(); ++i) {
    if (i != 0U) {
      stream << "; ";
    }
    stream << mismatch_messages[i];
  }
  return stream.str();
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

void validateRestartResumeTopologyOrThrowImpl(
    const io::RestartReadResult& restart,
    const core::SimulationConfig& config,
    const core::FrozenConfig& frozen_config,
    const parallel::MpiContext& mpi_context) {
  const std::uint8_t expected_max_bin = static_cast<std::uint8_t>(config.numerics.hierarchical_max_rung);
  const auto require_scheduler = [&](const core::TimeBinPersistentState& scheduler, std::string_view label) {
    if (scheduler.max_bin != expected_max_bin) {
      throw std::runtime_error("restart scheduler hierarchy differs from typed configuration: " + std::string(label));
    }
    if (expected_max_bin == 0U) {
      for (const auto bin : scheduler.bin_index) if (bin != 0U) throw std::runtime_error("reference restart has nonzero bin");
      for (const auto pending : scheduler.pending_bin_index) {
        if (pending != 0U && pending != core::HierarchicalTimeBinScheduler::k_unset_pending_bin) {
          throw std::runtime_error("reference restart has nonzero pending bin");
        }
      }
      return;
    }
    if (scheduler.next_activation_tick.size() != scheduler.bin_index.size() ||
        scheduler.active_flag.size() != scheduler.bin_index.size() ||
        scheduler.pending_bin_index.size() != scheduler.bin_index.size()) {
      throw std::runtime_error("hierarchical restart scheduler lanes have inconsistent extents");
    }
    const auto period = 1ULL << expected_max_bin;
    if (scheduler.current_tick % period != 0U) throw std::runtime_error("hierarchical restart is not at a coarse synchronization tick");
    for (std::size_t row = 0U; row < scheduler.bin_index.size(); ++row) {
      if (scheduler.bin_index[row] > expected_max_bin || scheduler.active_flag[row] != 0U ||
          scheduler.pending_bin_index[row] != core::HierarchicalTimeBinScheduler::k_unset_pending_bin ||
          scheduler.next_activation_tick[row] != scheduler.current_tick) {
        throw std::runtime_error("hierarchical restart contains an open particle interval or pending transition");
      }
    }
  };
  require_scheduler(restart.scheduler_state, "particles");
  require_scheduler(restart.gas_cell_scheduler_state, "gas cells");
  if (expected_max_bin > 0U) {
    if (restart.state.cells.size() != 0U || restart.integrator_state.inside_kdk_step ||
        !restart.integrator_state.last_completed_restart_safe ||
        restart.scheduler_state.current_tick != restart.gas_cell_scheduler_state.current_tick ||
        restart.integrator_state.time_bins.max_bin != expected_max_bin ||
        restart.integrator_state.last_completed_boundary_kind != core::StepBoundaryKind::kGlobalSynchronizationPoint ||
        restart.scheduler_state.bin_index.size() != restart.state.particles.size() ||
        !restart.gas_cell_scheduler_state.bin_index.empty()) {
      throw std::runtime_error("hierarchical restart requires closed collisionless synchronized KDK authority");
    }
    for (std::size_t row = 0U; row < restart.state.particles.size(); ++row) {
      if (restart.state.particleSpeciesTag(row) != static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter) ||
          restart.state.particleLastDriftTimeCode(row) != restart.integrator_state.current_time_code ||
          restart.state.particleLastDriftScaleFactor(row) != restart.integrator_state.current_scale_factor) {
        throw std::runtime_error("hierarchical restart mixes species or particle drift epochs");
      }
    }
  }
  (void)io::restartConfigurationCompatibilityRecord(restart, frozen_config);
  const core::CudaRuntimeInfo cuda_runtime = core::queryCudaRuntime();
  const parallel::DistributedExecutionTopology runtime_topology =
      parallel::buildDistributedExecutionTopology(
          static_cast<std::size_t>(config.numerics.treepm_pm_grid_nx),
          static_cast<std::size_t>(config.numerics.treepm_pm_grid_ny),
          static_cast<std::size_t>(config.numerics.treepm_pm_grid_nz),
          mpi_context,
          config.parallel.mpi_ranks_expected,
          config.parallel.gpu_devices,
          cuda_runtime.runtime_available,
          cuda_runtime.visible_device_count,
          pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode));
  const parallel::DistributedRestartCompatibilityReport compatibility =
      parallel::evaluateDistributedRestartCompatibility(
          restart.distributed_gravity_state,
          runtime_topology);
  if (!compatibility.compatible()) {
    throw std::runtime_error(
        "ReferenceWorkflow restart topology validation failed for rank " +
        std::to_string(mpi_context.worldRank()) + "/" + std::to_string(mpi_context.worldSize()) +
        ": " + joinRestartCompatibilityMessages(compatibility.mismatch_messages));
  }
  if (restart.distributed_gravity_state.owning_rank_by_item.size() !=
      restart.state.particles.size()) {
    throw std::runtime_error(
        "ReferenceWorkflow restart topology validation failed: /distributed_gravity owning_rank_by_item count " +
        std::to_string(restart.distributed_gravity_state.owning_rank_by_item.size()) +
        " does not match local particle rows " +
        std::to_string(restart.state.particles.size()));
  }
  for (std::size_t row = 0; row < restart.state.particles.size(); ++row) {
    const int observed_owner = static_cast<int>(restart.state.particleOwningRank(row));
    const int restart_owner = restart.distributed_gravity_state.owning_rank_by_item[row];
    if (observed_owner != restart_owner) {
      throw std::runtime_error(
          "ReferenceWorkflow restart topology validation failed: local particle row " +
          std::to_string(row) + " owner mismatch: state=" + std::to_string(observed_owner) +
          ", /distributed_gravity=" + std::to_string(restart_owner));
    }
    if (restart_owner < 0 || restart_owner >= mpi_context.worldSize()) {
      throw std::runtime_error(
          "ReferenceWorkflow restart topology validation failed: local particle row " +
          std::to_string(row) + " owner rank is outside runtime world: owner=" +
          std::to_string(restart_owner) + ", world_size=" + std::to_string(mpi_context.worldSize()));
    }
  }
}

}  // namespace

void validateRestartResumeTopologyOrThrow(
    const io::RestartReadResult& restart,
    const core::SimulationConfig& config,
    const core::FrozenConfig& frozen_config,
    const parallel::MpiContext& mpi_context) {
  validateRestartResumeTopologyOrThrowImpl(
      restart, config, frozen_config, mpi_context);
}

PendingOutputBoundary initializeOutputCadence(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    const core::IntegratorState& integrator_state,
    const io::OutputCadencePersistentState* restored_output) {
  PendingOutputBoundary pending;
  if (options.write_outputs && config.output.snapshot_interval_steps > 0) {
    pending.snapshot_interval_steps =
        static_cast<std::uint64_t>(config.output.snapshot_interval_steps);
    if (restored_output != nullptr && restored_output->output_enabled) {
      if (restored_output->snapshot_interval_steps !=
          pending.snapshot_interval_steps) {
        throw std::runtime_error(
            "restart output step interval does not match the active normalized configuration");
      }
      pending.next_snapshot_step_index =
          restored_output->next_snapshot_step_index;
    } else {
      const std::uint64_t completed_intervals =
          integrator_state.step_index / pending.snapshot_interval_steps;
      if (completed_intervals == std::numeric_limits<std::uint64_t>::max() ||
          completed_intervals + 1U >
              std::numeric_limits<std::uint64_t>::max() /
                  pending.snapshot_interval_steps) {
        throw std::overflow_error(
            "output step cadence could not represent its next event");
      }
      pending.next_snapshot_step_index =
          (completed_intervals + 1U) * pending.snapshot_interval_steps;
    }
    if (pending.next_snapshot_step_index <= integrator_state.step_index) {
      throw std::runtime_error(
          "restart output step cadence does not provide a future ordered event");
    }
  } else if (restored_output != nullptr && options.write_outputs &&
             restored_output->output_enabled &&
             restored_output->snapshot_interval_steps > 0U) {
    throw std::runtime_error(
        "restart contains an active step output cadence that the current configuration disables");
  }

  if (options.write_outputs &&
      config.output.snapshot_interval_time_code > 0.0) {
    pending.snapshot_interval_time_code =
        config.output.snapshot_interval_time_code;
    if (restored_output != nullptr &&
        restored_output->snapshot_interval_time_code > 0.0) {
      if (!timelineTimesEqual(
              restored_output->snapshot_interval_time_code,
              config.output.snapshot_interval_time_code)) {
        throw std::runtime_error(
            "restart output code-time interval does not match the active normalized configuration");
      }
      pending.next_snapshot_time_code =
          restored_output->next_snapshot_time_code;
    } else {
      pending.next_snapshot_time_code = nextOrderedTimeEvent(
          config.numerics.t_code_begin,
          integrator_state.current_time_code,
          config.output.snapshot_interval_time_code);
    }
    if (!std::isfinite(pending.next_snapshot_time_code) ||
        pending.next_snapshot_time_code <= integrator_state.current_time_code) {
      throw std::runtime_error(
          "restart output code-time cadence does not provide a future ordered event");
    }
  } else if (restored_output != nullptr && options.write_outputs &&
             restored_output->snapshot_interval_time_code > 0.0) {
    throw std::runtime_error(
        "restart contains an active code-time output cadence that the current configuration disables");
  }

  if (options.write_outputs && config.output.write_restarts &&
      config.output.restart_interval_steps > 0) {
    pending.restart_interval_steps =
        static_cast<std::uint64_t>(config.output.restart_interval_steps);
    pending.next_restart_step_index = nextRestartStepIndexAfter(
        integrator_state.step_index, pending.restart_interval_steps);
  }
  if (options.write_outputs && config.output.write_restarts &&
      config.output.restart_interval_time_code > 0.0) {
    pending.restart_interval_time_code = config.output.restart_interval_time_code;
    pending.next_restart_time_code = nextOrderedTimeEvent(
        config.numerics.t_code_begin,
        integrator_state.current_time_code,
        config.output.restart_interval_time_code);
  }
  return pending;
}

std::uint64_t nextRestartStepIndexAfter(
    std::uint64_t step_index, std::uint64_t interval_steps) {
  if (interval_steps == 0U) {
    throw std::invalid_argument("checkpoint step cadence requires a positive interval");
  }
  const std::uint64_t completed_intervals = step_index / interval_steps;
  if (completed_intervals == std::numeric_limits<std::uint64_t>::max() ||
      completed_intervals + 1U >
          std::numeric_limits<std::uint64_t>::max() / interval_steps) {
    throw std::overflow_error("checkpoint step cadence could not represent its next event");
  }
  return (completed_intervals + 1U) * interval_steps;
}

[[nodiscard]] bool timelineTimesEqual(double lhs, double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= 8.0 * std::numeric_limits<double>::epsilon() * scale;
}

[[nodiscard]] double nextOrderedTimeEvent(
    double timeline_begin_time_code,
    double current_time_code,
    double interval_time_code) {
  if (!std::isfinite(timeline_begin_time_code) || !std::isfinite(current_time_code) ||
      !std::isfinite(interval_time_code) || interval_time_code <= 0.0) {
    throw std::invalid_argument("output code-time cadence requires finite times and a positive interval");
  }
  const double completed_intervals = std::max(
      0.0,
      std::floor((current_time_code - timeline_begin_time_code) / interval_time_code));
  double next_time_code =
      std::fma(completed_intervals + 1.0, interval_time_code, timeline_begin_time_code);
  if (next_time_code <= current_time_code || timelineTimesEqual(next_time_code, current_time_code)) {
    next_time_code = std::fma(completed_intervals + 2.0, interval_time_code, timeline_begin_time_code);
  }
  if (!std::isfinite(next_time_code) || next_time_code <= current_time_code) {
    throw std::overflow_error("output code-time cadence could not represent its next ordered event");
  }
  return next_time_code;
}

void latchOutputRequestForCompletedStep(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    std::uint64_t next_step_index,
    double next_time_code,
    PendingOutputBoundary& pending) {
  if (!options.write_outputs) {
    return;
  }
  const bool step_event_due = pending.snapshot_interval_steps > 0 &&
      next_step_index >= pending.next_snapshot_step_index;
  const bool time_event_due = pending.snapshot_interval_time_code > 0.0 &&
      timelineTimesEqual(next_time_code, pending.next_snapshot_time_code);
  const bool independent_restart_cadence =
      pending.restart_interval_steps > 0U || pending.restart_interval_time_code > 0.0;
  const bool restart_step_due = pending.restart_interval_steps > 0U &&
      next_step_index >= pending.next_restart_step_index;
  // No clipping: the first completed step ending at or after the event.
  const bool restart_time_due = pending.restart_interval_time_code > 0.0 &&
      (next_time_code >= pending.next_restart_time_code ||
       timelineTimesEqual(next_time_code, pending.next_restart_time_code));
  const bool restart_event_due = restart_step_due || restart_time_due;
  if (step_event_due || time_event_due) {
    pending.snapshot_due = true;
    if (!independent_restart_cadence) {
      // Legacy policy: each science snapshot event also writes a checkpoint.
      pending.checkpoint_due = pending.checkpoint_due || config.output.write_restarts;
    }
    pending.step_event_due = pending.step_event_due || step_event_due;
    pending.time_event_due = pending.time_event_due || time_event_due;
  }
  if (restart_event_due) {
    pending.checkpoint_due = true;
    pending.restart_step_event_due = pending.restart_step_event_due || restart_step_due;
    pending.restart_time_event_due = pending.restart_time_event_due || restart_time_due;
  }
}

[[nodiscard]] core::StepBoundaryKind requestedBoundaryForPendingOutput(const PendingOutputBoundary& pending) {
  if (pending.checkpoint_due) {
    return core::StepBoundaryKind::kCheckpointPoint;
  }
  if (pending.snapshot_due) {
    return core::StepBoundaryKind::kSnapshotPoint;
  }
  return core::StepBoundaryKind::kGlobalSynchronizationPoint;
}

OutputRestartRuntime::OutputRestartRuntime(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const GravityRestartStateProvider& gravity_state,
    const AnalysisRuntime& analysis_state,
    const RuntimeServices& services,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    PendingOutputBoundary& pending_output,
    bool write_outputs_enabled)
    : m_frozen_config(frozen_config),
      m_config(config),
      m_scheduler(scheduler),
      m_gas_cell_scheduler(gas_cell_scheduler),
      m_gravity_state(gravity_state),
      m_analysis_state(analysis_state),
      m_services(services),
      m_report(report),
      m_profiler(profiler),
      m_pending_output(pending_output),
      m_write_outputs_enabled(write_outputs_enabled) {}

RuntimeTaskMemoryEstimate OutputRestartRuntime::estimateMemory(
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state) const {
  if (!m_write_outputs_enabled ||
      (!m_pending_output.snapshot_due && !m_pending_output.checkpoint_due) ||
      integrator_state.step_index == 0U) {
    return {0U, true, {}};
  }
  if (!core::isOutputSafeBoundary(integrator_state.last_completed_boundary_kind)) {
    return {0U, true, {}};
  }
  const std::uint64_t state_bytes = simulationOwnedCapacityBytes(state);
  std::uint64_t peak = m_pending_output.snapshot_due
      ? incrementalOutputStagingBytes(m_services, state_bytes) : 0U;
  if (m_pending_output.checkpoint_due) {
    constexpr std::uint64_t k_restart_verification_workspace_bytes = 16ULL * 1024ULL * 1024ULL;
    // Write-time verification no longer constructs a second SimulationState,
    // scheduler population, or owning gravity force-cache replica. The live
    // state/cache are borrowed; application-owned verification staging is a
    // fixed chunk workspace. Genuine restore-candidate memory remains modeled
    // by restartReadCandidateStagingBytes() on the restore path.
    peak = std::max(peak, incrementalOutputStagingBytes(
        m_services, k_restart_verification_workspace_bytes));
  }
  return {peak, false,
          "snapshot staging plus bounded 16-MiB/rank streaming restart verification workspace"};
}

void OutputRestartRuntime::execute(OutputRestartStageView& view) {
  view.requireFresh();
  core::StepContext& context = internal::RuntimeStageAccess::outputRestartContext(
      view,
      {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kSchedulerTruth, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kOutputRestartState, RuntimeResourceAccessMode::kReadWrite},
       {RuntimeResourceKey::kDiagnostics, RuntimeResourceAccessMode::kWrite}});
  if (context.stage != core::IntegrationStage::kOutputCheck) {
    throw std::logic_error(
        "output/restart runtime received an unregistered stage");
  }
  if (!m_pending_output.snapshot_due && !m_pending_output.checkpoint_due) {
    return;
  }
  if (!context.boundary.output_safe || !context.boundary.restart_safe) {
    if (context.boundary.multirate_pm_cycle_interior && !m_pending_output.time_event_due) {
      // Synchronized boundary inside a multirate PM cycle: momenta carry an
      // open outer PM half-kick. Step-cadence snapshots and checkpoints stay
      // latched and are written at the cycle end (at most M - 1 blocks later).
      // Code-time snapshot events are planned onto cycle ends and fall
      // through to the fail-closed check below if that contract is violated.
      m_profiler.counters().addCount("io.output.deferred_to_multirate_pm_cycle_end");
      return;
    }
    if (m_pending_output.time_event_due) {
      throw std::runtime_error(
          "required code-time output event reached a boundary that is not output/restart safe");
    }
    if (m_pending_output.checkpoint_due) {
      core::assertCanWriteCheckpointAtBoundary(
          context.integrator_state, m_scheduler.currentTick());
    }
    return;
  }

  const double persisted_next_snapshot_time_code =
      m_pending_output.time_event_due
      ? nextOrderedTimeEvent(
            m_config.numerics.t_code_begin,
            context.integrator_state.current_time_code,
            m_pending_output.snapshot_interval_time_code)
      : m_pending_output.next_snapshot_time_code;
  std::uint64_t persisted_next_snapshot_step_index =
      m_pending_output.next_snapshot_step_index;
  if (m_pending_output.step_event_due) {
    do {
      if (persisted_next_snapshot_step_index >
          std::numeric_limits<std::uint64_t>::max() -
              m_pending_output.snapshot_interval_steps) {
        throw std::overflow_error(
            "output step cadence overflowed its ordered timeline");
      }
      persisted_next_snapshot_step_index +=
          m_pending_output.snapshot_interval_steps;
    } while (persisted_next_snapshot_step_index <=
             context.integrator_state.step_index);
  }

  const bool snapshot_requested = m_pending_output.snapshot_due;
  // Experimental stale-window PM reuse (cadence > 1): a resumed run rebuilds
  // the long-range field at its first force opportunity, so a checkpoint is
  // written only where an uninterrupted run would refresh at that same
  // opportunity; otherwise it stays pending to the next aligned boundary.
  const auto& pm_sync = context.integrator_state.pm_sync_state;
  const bool checkpoint_pm_aligned = pm_sync.cadenceSteps() <= 1U ||
      !context.integrator_state.pm_long_range_field_valid ||
      (pm_sync.gravityKickOpportunity() + 1U - pm_sync.lastRefreshOpportunity()) >= pm_sync.cadenceSteps();
  const bool checkpoint_now = m_pending_output.checkpoint_due && checkpoint_pm_aligned;
  if (m_pending_output.checkpoint_due && !checkpoint_now) {
    m_profiler.counters().addCount("io.restart.deferred_for_pm_cadence");
  }
  if (checkpoint_now) {
    // CHUI Part 2: serialize the gas composition sidecar into its module
    // block right before the checkpoint, so the restart payload (and its
    // integrity digest) carries it. Only a serialized image; the sidecar
    // stays the single live authority.
    physics::packGasCompositionIntoModuleSidecar(context.state);
  }
  const bool output_flushed = maybeWriteOutputs(
      m_frozen_config,
      m_config,
      context.state,
      context.integrator_state,
      m_scheduler,
      m_gas_cell_scheduler,
      m_gravity_state,
      m_analysis_state,
      m_services,
      m_report,
      m_profiler,
      m_write_outputs_enabled,
      m_pending_output.snapshot_due,
      checkpoint_now,
      m_pending_output.snapshot_interval_steps,
      persisted_next_snapshot_step_index,
      m_pending_output.snapshot_interval_time_code,
      persisted_next_snapshot_time_code,
      m_pending_output.restart_resume_dt_time_code);
  if (output_flushed) {
    if (snapshot_requested) {
      m_last_committed_snapshot_step_index = context.integrator_state.step_index;
    }
    m_pending_output.snapshot_due = false;
    m_pending_output.step_event_due = false;
    m_pending_output.time_event_due = false;
    if (checkpoint_now) {
      m_pending_output.checkpoint_due = false;
      if (m_pending_output.restart_step_event_due) {
        m_pending_output.next_restart_step_index = nextRestartStepIndexAfter(
            context.integrator_state.step_index, m_pending_output.restart_interval_steps);
        m_pending_output.restart_step_event_due = false;
      }
      if (m_pending_output.restart_time_event_due) {
        m_pending_output.next_restart_time_code = nextOrderedTimeEvent(
            m_config.numerics.t_code_begin,
            context.integrator_state.current_time_code,
            m_pending_output.restart_interval_time_code);
        m_pending_output.restart_time_event_due = false;
      }
    }
    m_pending_output.next_snapshot_step_index =
        persisted_next_snapshot_step_index;
    m_pending_output.next_snapshot_time_code =
        persisted_next_snapshot_time_code;
    m_pending_output.restart_resume_dt_time_code = 0.0;
  }
}

void OutputRestartRuntime::ensureFinalEndpointSnapshot(
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state) {
  if (!m_write_outputs_enabled || integrator_state.step_index == 0U) return;
  const double endpoint_tolerance = 1.0e-12 * std::max(
      {1.0, std::abs(integrator_state.current_time_code),
       std::abs(m_config.numerics.t_code_end)});
  if (integrator_state.current_time_code + endpoint_tolerance <
      m_config.numerics.t_code_end) {
    return;  // bounded segment/restart probe, not authoritative normal endpoint
  }
  core::assertCanWriteSnapshotAtBoundary(integrator_state);
  if (m_last_committed_snapshot_step_index.has_value() &&
      *m_last_committed_snapshot_step_index == integrator_state.step_index) {
    return;
  }

  const bool flushed = maybeWriteOutputs(
      m_frozen_config, m_config, state, integrator_state, m_scheduler,
      m_gas_cell_scheduler, m_gravity_state, m_analysis_state, m_services,
      m_report, m_profiler, true,
      true, false,
      m_pending_output.snapshot_interval_steps,
      m_pending_output.next_snapshot_step_index,
      m_pending_output.snapshot_interval_time_code,
      m_pending_output.next_snapshot_time_code,
      0.0);
  if (!flushed) {
    throw std::runtime_error(
        "failed to commit required final endpoint science snapshot at an output-safe boundary");
  }
  m_last_committed_snapshot_step_index = integrator_state.step_index;
}

}  // namespace cosmosim::workflows::internal
