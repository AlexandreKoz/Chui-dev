#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/gravity_source_ownership.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/workflows/source_runtime.hpp"
#include "cosmosim/workflows/runtime_capabilities.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/time_coordinator.hpp"
#include "cosmosim/workflows/analysis_runtime.hpp"
#include "workflows/internal/initial_condition_runtime.hpp"
#include "cosmosim/workflows/migration_balance_runtime.hpp"
#include "workflows/internal/output_restart_runtime.hpp"
#include "workflows/internal/reference_runtime_composition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <cstring>
#include <exception>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/openmp_runtime.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/node_memory_admission.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace cosmosim::workflows {
namespace {

[[nodiscard]] bool runtimePhaseDiagnosticsEnabled() {
  const char* value = std::getenv("COSMOSIM_RUNTIME_PHASE_DIAGNOSTICS");
  return value != nullptr && std::string_view(value) != "0";
}

void traceRuntimePhase(std::string_view phase) {
  if (runtimePhaseDiagnosticsEnabled()) {
    std::cerr << "runtime_phase=" << phase << '\n' << std::flush;
  }
}

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}

[[nodiscard]] std::optional<double> consoleRedshift(
    bool cosmological,
    double a_scale) {
  if (!cosmological || !std::isfinite(a_scale) || a_scale <= 0.0) {
    return std::nullopt;
  }
  return 1.0 / a_scale - 1.0;
}


[[nodiscard]] std::string treePmAssignmentSchemeName(core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return "cic";
    case core::TreePmAssignmentScheme::kTsc:
      return "tsc";
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}


[[nodiscard]] std::filesystem::path computeRunDirectory(
    const core::SimulationConfig& config,
    const std::filesystem::path* output_root_override) {
  const std::filesystem::path output_root =
      (output_root_override != nullptr && !output_root_override->empty())
      ? *output_root_override
      : std::filesystem::path(config.output.output_directory);
  return output_root / config.output.run_name;
}

[[nodiscard]] std::filesystem::path rankQualifiedRunDirectory(
    std::filesystem::path run_directory,
    int world_size,
    int world_rank) {
  if (world_size <= 1) {
    return run_directory;
  }
  std::ostringstream rank_suffix;
  rank_suffix << "_rank" << std::setw(3) << std::setfill('0') << world_rank;
  return run_directory.parent_path() / (run_directory.filename().string() + rank_suffix.str());
}

[[nodiscard]] std::filesystem::path resolveConfigRelativePath(
    const core::FrozenConfig& frozen_config,
    const std::filesystem::path& candidate) {
  if (candidate.is_absolute()) {
    return candidate;
  }

  const std::filesystem::path source_path(frozen_config.provenance.source_name);
  if (!source_path.empty() && source_path.has_parent_path()) {
    return source_path.parent_path() / candidate;
  }
  return candidate;
}

[[nodiscard]] bool repeatedCanonicalOrder(const std::vector<std::string>& observed) {
  const auto canonical = core::StageScheduler::kickDriftKickOrder();
  if (observed.empty() || observed.size() % canonical.size() != 0) {
    return false;
  }

  for (std::size_t chunk = 0; chunk < observed.size(); chunk += canonical.size()) {
    for (std::size_t i = 0; i < canonical.size(); ++i) {
      if (observed[chunk + i] != core::integrationStageName(canonical[i])) {
        return false;
      }
    }
  }
  return true;
}

void fnv1aMix(std::uint64_t& hash, const void* data, std::size_t size) {
  constexpr std::uint64_t k_fnv_offset_basis = 1469598103934665603ULL;
  constexpr std::uint64_t k_fnv_prime = 1099511628211ULL;
  if (hash == 0) {
    hash = k_fnv_offset_basis;
  }
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= k_fnv_prime;
  }
}

void fnv1aMix(std::uint64_t& hash, double value) { fnv1aMix(hash, &value, sizeof(value)); }

void fnv1aMix(std::uint64_t& hash, std::uint64_t value) { fnv1aMix(hash, &value, sizeof(value)); }

[[nodiscard]] std::uint64_t entityStateHashParticle(
    const core::SimulationState& state,
    std::uint32_t particle_row) {
  if (particle_row >= state.particles.size()) {
    throw std::out_of_range("state digest particle row is outside ParticleSoa");
  }
  std::uint64_t hash = 0;
  fnv1aMix(hash, static_cast<std::uint64_t>(0x7061727469636c65ULL));  // "particle" domain tag.
  fnv1aMix(hash, state.particle_sidecar.particle_id[particle_row]);
  fnv1aMix(hash, static_cast<std::uint64_t>(
      state.particleSpeciesTag(particle_row)));
  fnv1aMix(hash, state.particles.position_x_comoving[particle_row]);
  fnv1aMix(hash, state.particles.position_y_comoving[particle_row]);
  fnv1aMix(hash, state.particles.position_z_comoving[particle_row]);
  fnv1aMix(hash, state.particles.velocity_x_peculiar[particle_row]);
  fnv1aMix(hash, state.particles.velocity_y_peculiar[particle_row]);
  fnv1aMix(hash, state.particles.velocity_z_peculiar[particle_row]);
  fnv1aMix(hash, state.particles.mass_code[particle_row]);
  return hash;
}

[[nodiscard]] std::uint64_t entityStateHashGasCell(
    const core::SimulationState& state,
    std::uint32_t cell_row) {
  if (cell_row >= state.cells.size() || cell_row >= state.gas_cells.size()) {
    throw std::out_of_range("state digest gas-cell row is outside authoritative gas storage");
  }
  std::uint64_t hash = 0;
  fnv1aMix(hash, static_cast<std::uint64_t>(0x6761735f63656c6cULL));  // "gas_cell" domain tag.
  fnv1aMix(hash, state.gas_cells.gas_cell_id[cell_row]);
  fnv1aMix(hash, state.cells.center_x_comoving[cell_row]);
  fnv1aMix(hash, state.cells.center_y_comoving[cell_row]);
  fnv1aMix(hash, state.cells.center_z_comoving[cell_row]);
  fnv1aMix(hash, state.cells.mass_code[cell_row]);
  fnv1aMix(hash, state.gas_cells.velocity_x_peculiar[cell_row]);
  fnv1aMix(hash, state.gas_cells.velocity_y_peculiar[cell_row]);
  fnv1aMix(hash, state.gas_cells.velocity_z_peculiar[cell_row]);
  fnv1aMix(hash, state.gas_cells.density_code[cell_row]);
  fnv1aMix(hash, state.gas_cells.pressure_code[cell_row]);
  fnv1aMix(hash, state.gas_cells.internal_energy_code[cell_row]);
  fnv1aMix(hash, state.gas_cells.metal_mass_code[cell_row]);
  return hash;
}

[[nodiscard]] std::uint64_t computeStateDigest(
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state,
    const parallel::MpiContext& mpi_context) {
  // Build an order- and decomposition-independent digest over authoritative
  // physical entities.  Per-entity hashes are keyed by stable identity and
  // combined commutatively; ghost copies, gas compatibility particles, and
  // covered coarse gas cells are intentionally excluded.  This makes the
  // report usable as a reproducibility artifact rather than a dense-row-local
  // checksum while avoiding a global O(N) gather/sort.
  const std::uint32_t local_rank = static_cast<std::uint32_t>(
      std::max(mpi_context.worldRank(), 0));
  const auto authoritative_rows = internal::selectAuthoritativeGravitySourceRows(
      state,
      local_rank,
      static_cast<std::uint32_t>(std::max(mpi_context.worldSize(), 1)),
      "reference workflow final-state digest");

  std::uint64_t local_sum = 0U;
  std::uint64_t local_xor = 0U;
  for (const std::uint32_t particle_row : authoritative_rows.particle_rows) {
    const std::uint64_t entity_hash = entityStateHashParticle(state, particle_row);
    local_sum += entity_hash;
    local_xor ^= entity_hash;
  }
  for (const std::uint32_t cell_row : authoritative_rows.gas_cell_rows) {
    const std::uint64_t entity_hash = entityStateHashGasCell(state, cell_row);
    local_sum += entity_hash;
    local_xor ^= entity_hash;
  }

  const std::uint64_t global_sum = mpi_context.allreduceSumUint64(local_sum);
  const std::uint64_t global_xor = mpi_context.allreduceXorUint64(local_xor);
  const std::uint64_t global_particle_count = mpi_context.allreduceSumUint64(
      static_cast<std::uint64_t>(authoritative_rows.particle_rows.size()));
  const std::uint64_t global_gas_cell_count = mpi_context.allreduceSumUint64(
      static_cast<std::uint64_t>(authoritative_rows.gas_cell_rows.size()));

  std::uint64_t hash = 0;
  fnv1aMix(hash, static_cast<std::uint64_t>(0x676c6f62616c7631ULL));  // global digest schema v1.
  fnv1aMix(hash, global_particle_count);
  fnv1aMix(hash, global_gas_cell_count);
  fnv1aMix(hash, global_sum);
  fnv1aMix(hash, global_xor);
  fnv1aMix(hash, integrator_state.current_time_code);
  fnv1aMix(hash, integrator_state.current_scale_factor);
  fnv1aMix(hash, integrator_state.step_index);
  return hash;
}

[[nodiscard]] std::uint64_t computeParticleIdSum(const core::SimulationState& state) {
  std::uint64_t sum = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    sum += particle_id;
  }
  return sum;
}

[[nodiscard]] std::uint64_t computeParticleIdSquareSum(const core::SimulationState& state) {
  std::uint64_t sum = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    sum += particle_id * particle_id;
  }
  return sum;
}

[[nodiscard]] std::uint64_t computeParticleIdXor(const core::SimulationState& state) {
  std::uint64_t value = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    value ^= particle_id;
  }
  return value;
}


void ensureRunDirectory(const std::filesystem::path& run_directory) {
  std::filesystem::create_directories(run_directory);
}

void flushCommonArtifacts(
    const core::FrozenConfig& frozen_config,
    core::ProfilerSession& profiler,
    ReferenceWorkflowReport& report) {
  ensureRunDirectory(report.run_directory);

  report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
  if (!std::filesystem::exists(report.normalized_config_snapshot_path)) {
    core::writeNormalizedConfigSnapshot(frozen_config, report.run_directory);
  }
  report.normalized_config_snapshot_written = true;

  report.profiler_json_path = report.run_directory / "profile.json";
  report.profiler_csv_path = report.run_directory / "profile.csv";
  report.operational_report_json_path = report.run_directory / "operational_events.json";
  core::writeProfilerReportJson(profiler, report.profiler_json_path);
  core::writeProfilerReportCsv(profiler, report.profiler_csv_path);
  core::writeOperationalReportJson(
      profiler,
      report.operational_report_json_path,
      frozen_config.config.output.run_name,
      frozen_config.provenance.config_hash_hex);
}


}  // namespace

ReferenceWorkflowRunner::ReferenceWorkflowRunner(core::FrozenConfig frozen_config)
    : m_frozen_config(std::move(frozen_config)) {
  core::configureOpenMpThreads(m_frozen_config.config.parallel.omp_threads);
}

const core::FrozenConfig& ReferenceWorkflowRunner::frozenConfig() const noexcept {
  return m_frozen_config;
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(const ReferenceWorkflowOptions& options) const {
  return runImpl(nullptr, options);
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(
    const std::filesystem::path& output_root_override,
    const ReferenceWorkflowOptions& options) const {
  return runImpl(&output_root_override, options);
}


ReferenceWorkflowReport ReferenceWorkflowRunner::runImpl(
    const std::filesystem::path* output_root_override,
    const ReferenceWorkflowOptions& options) const {
  // Runtime output-root overrides must cover diagnostics/analysis as well as
  // snapshots and restarts. Preserve the frozen configuration for provenance,
  // but pass a runtime copy with the effective root to output-producing modules.
  core::SimulationConfig runtime_config = m_frozen_config.config;
  if (output_root_override != nullptr && !output_root_override->empty()) {
    runtime_config.output.output_directory = output_root_override->string();
  }
  const core::SimulationConfig& config = runtime_config;
  parallel::MpiContext mpi_context;
  mpi_context.validateExpectedWorldSizeOrThrow(config.parallel.mpi_ranks_expected);

  ReferenceWorkflowReport report;
  report.world_size = mpi_context.worldSize();
  report.world_rank = mpi_context.worldRank();
  report.shared_run_directory = computeRunDirectory(config, output_root_override);
  report.run_directory = rankQualifiedRunDirectory(
      report.shared_run_directory, report.world_size, report.world_rank);
  report.config_compatible = true;
  report.schema_compatible =
      config.schema_version == 1 && io::gadgetArepoSchemaMap().schema_version >= 2 &&
      io::isRestartSchemaCompatible(io::restartSchema().version);

  core::ProfilerSession profiler(true);
  RuntimeConsoleOptions console_options = options.console;
  if (console_options.config_path.empty()) {
    console_options.config_path = m_frozen_config.provenance.source_name;
  }
  RuntimeConsoleReporter console_reporter(
      std::move(console_options), mpi_context.worldRank(), std::cout, std::cerr);
  std::uint64_t external_runtime_reserve_bytes = 0U;
  external_runtime_reserve_bytes = core::checkedMemoryBytesAdd(
      external_runtime_reserve_bytes,
      config.parallel.gravity_backend_unknown_reserve_bytes,
      "reference workflow gravity backend reserve");
  external_runtime_reserve_bytes = core::checkedMemoryBytesAdd(
      external_runtime_reserve_bytes,
      config.parallel.process_mpi_unknown_reserve_bytes,
      "reference workflow MPI reserve");
  external_runtime_reserve_bytes = core::checkedMemoryBytesAdd(
      external_runtime_reserve_bytes,
      config.parallel.process_fftw_unknown_reserve_bytes,
      "reference workflow FFTW reserve");
  external_runtime_reserve_bytes = core::checkedMemoryBytesAdd(
      external_runtime_reserve_bytes,
      config.parallel.process_hdf5_unknown_reserve_bytes,
      "reference workflow HDF5 reserve");
  external_runtime_reserve_bytes = core::checkedMemoryBytesAdd(
      external_runtime_reserve_bytes,
      config.parallel.process_allocator_unknown_reserve_bytes,
      "reference workflow allocator reserve");
  // Node-aware admission (parallel.node_memory_admission). Off preserves the
  // historical per-process semantics (0 = unlimited). shared_equal derives a
  // conservative per-rank ceiling from the node limit; every rank computes
  // the same rule locally, so no collective is needed and none can block.
  std::uint64_t rank_hard_limit_bytes = config.parallel.process_memory_budget_bytes;
  std::optional<core::NodeMemoryAdmissionResult> node_admission;
  if (config.parallel.node_memory_admission == core::NodeMemoryAdmissionMode::kSharedEqual) {
    std::exception_ptr node_admission_failure;
    try {
      node_admission = core::computeNodeMemoryAdmission(core::NodeMemoryAdmissionInput{
          .explicit_node_budget_bytes = config.parallel.node_memory_budget_bytes,
          .node_reserve_bytes = config.parallel.node_memory_reserve_bytes,
          .ranks_on_node = static_cast<std::uint32_t>(std::max(1, mpi_context.localSize())),
          .process_budget_bytes = config.parallel.process_memory_budget_bytes,
          .discovery = core::discoverNodeMemory()});
      rank_hard_limit_bytes = node_admission->rank_hard_limit_bytes;
    } catch (...) {
      node_admission_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        node_admission_failure, "node memory admission");
  }
  core::MemoryGovernor memory_governor(core::MemoryGovernorPolicy{
      .hard_limit_bytes = rank_hard_limit_bytes,
      .external_runtime_reserve_bytes = external_runtime_reserve_bytes,
      .planned_overlap_reserve_bytes =
          config.parallel.process_output_restart_overlap_bytes,
      .safety_margin_fraction =
          config.parallel.process_memory_safety_margin_fraction,
  });
  const RuntimeServices runtime_services{
      .mpi_context = mpi_context,
      .profiler = profiler,
      .memory_governor = &memory_governor,
      .console_reporter = &console_reporter,
      .deterministic_execution = true};
  const FailureCoordinator failure_coordinator(runtime_services);
  const internal::MigrationBalanceRuntime migration_balance(
      config, runtime_services);
  if (node_admission.has_value()) {
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "memory.node_admission",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "core.memory_governor",
        .step_index = options.step_index,
        .simulation_time_code = config.numerics.t_code_begin,
        .scale_factor = 1.0,
        .message = "node-aware memory admission selected the per-rank hard limit",
        .payload = {{"node_limit_bytes", std::to_string(node_admission->node_limit_bytes)},
                    {"node_limit_source",
                     std::string(core::nodeMemoryLimitSourceLabel(node_admission->node_limit_source))},
                    {"node_usable_bytes", std::to_string(node_admission->node_usable_bytes)},
                    {"ranks_on_node", std::to_string(node_admission->ranks_on_node)},
                    {"rank_share_bytes", std::to_string(node_admission->rank_share_bytes)},
                    {"rank_hard_limit_bytes", std::to_string(node_admission->rank_hard_limit_bytes)},
                    {"process_budget_binding", node_admission->process_budget_was_binding ? "true" : "false"}},
    });
  }
  const core::OpenMpRuntimeInfo omp_runtime = core::openMpRuntimeInfo();
  profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "runtime.topology",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = "parallel.runtime",
      .step_index = options.step_index,
      .simulation_time_code = config.numerics.t_code_begin,
      .scale_factor = 1.0,
      .message = "MPI/OpenMP execution topology selected",
      .payload = {{"world_size", std::to_string(mpi_context.worldSize())},
                  {"world_rank", std::to_string(mpi_context.worldRank())},
                  {"node_local_size", std::to_string(mpi_context.localSize())},
                  {"node_local_rank", std::to_string(mpi_context.localRank())},
                  {"openmp_compiled", omp_runtime.compiled ? "true" : "false"},
                  {"openmp_requested_threads", std::to_string(omp_runtime.requested_threads)},
                  {"openmp_configured_threads", std::to_string(omp_runtime.configured_threads)}},
  });
  profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "config.freeze",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = "core.config",
      .step_index = options.step_index,
      .simulation_time_code = config.numerics.t_code_begin,
      .scale_factor = 1.0,
      .message = "frozen configuration accepted for runtime workflow",
      .payload = {{"schema_version", std::to_string(config.schema_version)},
                  {"config_hash_hex", m_frozen_config.provenance.config_hash_hex},
                  {"source_name", m_frozen_config.provenance.source_name}},
  });

  try {
    ensureRunDirectory(report.run_directory);
    const RuntimeCapabilityReport runtime_capabilities =
        buildRuntimeCapabilityReport(config);
    validateRequestedRuntimeCapabilities(config, runtime_capabilities);
    report.runtime_capability_report_path =
        report.run_directory / "runtime_capabilities.json";
    writeRuntimeCapabilityReportJson(
        runtime_capabilities, report.runtime_capability_report_path);
    core::writeNormalizedConfigSnapshot(m_frozen_config, report.run_directory);
    report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
    report.normalized_config_snapshot_written = true;

    if (!report.schema_compatible) {
      throw std::runtime_error("runtime workflow schema compatibility validation failed");
    }

#if !COSMOSIM_ENABLE_HDF5
    if (options.write_outputs) {
      throw std::runtime_error(
          "this build lacks HDF5 support, so the config-driven runtime cannot emit snapshots/restarts; reconfigure with COSMOSIM_ENABLE_HDF5=ON or use a no-output internal test path");
    }
#endif

    const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
    core::validateModePolicy(config, mode_policy);

    console_reporter.emitRuntimePhase(
        "initial_conditions_begin", "materializing authoritative initial state");
    traceRuntimePhase("initial_conditions_begin");
    const internal::InitialConditionRuntime initial_conditions(
        m_frozen_config, runtime_services);
    internal::InitialConditionStartupResult startup =
        initial_conditions.materialize(options, report.run_directory);
    traceRuntimePhase("initial_conditions_complete");
    report.ic_manifest_path = startup.manifest_path;
    core::SimulationState state = std::move(startup.state);
    const std::array initial_memory_reports{
        core::collectSimulationMemoryReport(state),
        profiler.retainedEventMemoryReport()};
    core::MemoryReport startup_memory_report_value =
        core::mergeMemoryReports(initial_memory_reports);
    memory_governor.setBaselineOwnedBytes(
        core::memoryReportBaselineOwnedBytes(startup_memory_report_value));
    core::attachMemoryGovernorSnapshot(
        startup_memory_report_value, memory_governor);
    core::attachProcessMemoryObservation(
        startup_memory_report_value, core::observeProcessMemory());
    attachDistributedMemoryTelemetry(
        startup_memory_report_value, runtime_services);
    profiler.setMemoryReport(std::move(startup_memory_report_value));
    const core::MemoryReport* startup_memory_report = profiler.memoryReport();
    if (startup_memory_report != nullptr) {
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "memory.startup_snapshot",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.memory",
          .step_index = options.step_index,
          .simulation_time_code = config.numerics.t_code_begin,
          .scale_factor = 1.0,
          .message = "startup memory accounting snapshot from owned host buffers",
          .payload = {{"persistent_total_bytes", std::to_string(startup_memory_report->totals.persistent_total_bytes)},
                      {"transient_total_bytes", std::to_string(startup_memory_report->totals.transient_total_bytes)},
                      {"unknown_external_allocations", "true"}},
      });
    }
    initializeGravityState(state, config);
    const bool restoring_from_restart = startup.restoring_from_restart;
    if (restoring_from_restart) {
      std::exception_ptr local_restart_validation_failure;
      try {
        internal::validateRestartResumeTopologyOrThrow(
            *options.restart_state_override, config, m_frozen_config,
            mpi_context);
        report.restart_config_compatibility_record = io::restartConfigurationCompatibilityRecord(
            *options.restart_state_override, m_frozen_config);
      } catch (...) {
        local_restart_validation_failure = std::current_exception();
      }
      failure_coordinator.rethrowCollectiveFailure(
          local_restart_validation_failure, "restart topology validation");
    } else if (!startup.already_partitioned) {
      migration_balance.initializeOwnership(state);
    }

    std::exception_ptr local_identity_validation_failure;
    try {
      if (!state.hasCertifiedPersistentParticleIds()) {
        core::MemoryReservation identity_validation_reservation =
            memory_governor.reserve(
                core::MemoryClass::kDiagnostic,
                core::OwnershipValidationWorkspace::k_id_validation_workspace_limit_bytes,
                "core.startup.particle_id_validation");
        identity_validation_reservation.commit();
        core::OwnershipValidationWorkspace identity_validation_scratch;
        if (!state.validatePersistentParticleIds(identity_validation_scratch)) {
          throw std::runtime_error(
              "startup authoritative state contains duplicate or zero persistent particle IDs");
        }
      }
    } catch (...) {
      local_identity_validation_failure = std::current_exception();
    }
    failure_coordinator.rethrowCollectiveFailure(
        local_identity_validation_failure, "startup particle-ID validation");

    const parallel::LocalOwnershipIdentitySummary expected_global_identity =
        migration_balance.reduceIdentity(state);
    std::vector<std::uint64_t> expected_global_particle_ids(
        state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());

    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    const std::uint64_t local_particle_id_square_sum = computeParticleIdSquareSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    const auto local_identity = parallel::summarizeLocalOwnedParticleIds(
        state.particle_sidecar.particle_id, state.validateUniqueParticleIds());
    report.local_particle_ids_unique = local_identity.local_particle_ids_unique;
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    const std::uint64_t global_particle_id_square_sum =
        mpi_context.allreduceSumUint64(local_particle_id_square_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    const bool all_ranks_have_unique_local_ids =
        mpi_context.allreduceSumUint64(report.local_particle_ids_unique ? 1ULL : 0ULL) ==
        static_cast<std::uint64_t>(mpi_context.worldSize());
    const parallel::LocalOwnershipIdentitySummary reduced_identity{
        .local_owned_count = report.global_particle_count,
        .local_particle_id_sum = report.global_particle_id_sum,
        .local_particle_id_square_sum = global_particle_id_square_sum,
        .local_particle_id_xor = report.global_particle_id_xor,
        .local_particle_ids_unique = all_ranks_have_unique_local_ids,
    };
    report.global_particle_partition_identity_match = parallel::partitionIdentityMatchesGeneratedSet(
        reduced_identity,
        expected_global_identity.local_owned_count,
        expected_global_identity.local_particle_id_sum,
        expected_global_identity.local_particle_id_square_sum,
        expected_global_identity.local_particle_id_xor);
    if (!report.global_particle_partition_identity_match) {
      throw std::runtime_error(
          "distributed ownership invariant failed after initial decomposition: duplicate, missing, or extra authoritative particle IDs");
    }
    // Homogeneous-DMO metadata is a global representation contract and must
    // be established before scheduler/runtime owner construction.  GravityRuntime
    // selects source representations during composition; compacting afterward is
    // too late and can leave ranks on different MPI protocols when some ranks own
    // zero particles.  report.global_cell_count is already a collective result.
    if (report.global_cell_count == 0U) {
      std::exception_ptr local_dmo_compaction_failure;
      try {
        if (!state.hasHomogeneousDmoMetadata() &&
            !state.compactHomogeneousDmoMetadata(
                static_cast<std::uint32_t>(std::max(mpi_context.worldRank(), 0)))) {
          throw std::runtime_error(
              "globally gas-free runtime could not establish homogeneous DMO metadata on this rank");
        }
      } catch (...) {
        local_dmo_compaction_failure = std::current_exception();
      }
      failure_coordinator.rethrowCollectiveFailure(
          local_dmo_compaction_failure, "global DMO metadata compaction");
    }

    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.ownership.partition_summary",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.ownership",
        .step_index = options.step_index,
        .simulation_time_code = config.numerics.t_code_begin,
        .scale_factor = 1.0,
        .message = "initial ownership-partition summary recorded",
        .payload = {{"world_size", std::to_string(report.world_size)},
                    {"world_rank", std::to_string(report.world_rank)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_cell_count", std::to_string(report.local_cell_count)},
                    {"global_cell_count", std::to_string(report.global_cell_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)},
                    {"local_particle_ids_unique", report.local_particle_ids_unique ? "true" : "false"},
                    {"global_particle_partition_identity_match",
                     report.global_particle_partition_identity_match ? "true" : "false"}},
    });
    console_reporter.emitInitialConditions(
        report.global_particle_count,
        report.global_cell_count,
        restoring_from_restart,
        report.ic_manifest_path);

    core::CosmologyBackgroundConfig background_config;
    background_config.hubble_param = config.cosmology.hubble_param;
    background_config.omega_matter = config.cosmology.omega_matter;
    background_config.omega_lambda = config.cosmology.omega_lambda;
    const std::optional<core::LambdaCdmBackground> background = mode_policy.cosmological_comoving_frame
        ? std::optional<core::LambdaCdmBackground>(core::LambdaCdmBackground(background_config))
        : std::nullopt;

    const core::UnitSystem runtime_units = core::makeUnitSystem(
        config.units.length_unit,
        config.units.mass_unit,
        config.units.velocity_unit);
    // Hierarchical scheduler lanes are optional population-scale owners. Admit
    // their first materialization before reset/import, then reconcile actual
    // capacities. Later bin moves use RetainedCapacityTransaction.
    core::MemoryReservation scheduler_initialization_reservation;
    const std::uint64_t state_bytes_before_scheduler = core::memoryReportBaselineOwnedBytes(
        core::collectSimulationMemoryReport(state));
    std::exception_ptr scheduler_initialization_failure;
    if (config.numerics.hierarchical_max_rung > 0) {
      try {
        const std::uint64_t per_particle = sizeof(std::uint64_t) + sizeof(std::size_t) +
            5U * sizeof(std::uint8_t) + sizeof(core::TimeStepCandidateSource) + 3U * sizeof(std::uint32_t);
        const auto population_bytes = static_cast<std::uint64_t>(core::checkedSizeMultiply(
            state.particles.size(), static_cast<std::size_t>(per_particle), "hierarchical scheduler initial bytes"));
        const std::uint64_t fixed_bytes = static_cast<std::uint64_t>(core::checkedSizeMultiply(
            2U * (static_cast<std::size_t>(config.numerics.hierarchical_max_rung) + 1U),
            sizeof(std::vector<std::uint32_t>) + 2U * sizeof(std::uint32_t), "hierarchical bin table bytes"));
        scheduler_initialization_reservation = memory_governor.reserve(core::MemoryClass::kPhaseResident,
            core::checkedMemoryBytesAdd(population_bytes, fixed_bytes, "hierarchical scheduler initial peak"),
            "time.hierarchical_scheduler_initialization");
        scheduler_initialization_reservation.commit();
      } catch (...) { scheduler_initialization_failure = std::current_exception(); }
    }
    failure_coordinator.rethrowCollectiveFailure(scheduler_initialization_failure, "scheduler initial memory admission");
    std::optional<RungZeroTimeState> time_state_owner;
    scheduler_initialization_failure = {};
    try {
      time_state_owner.emplace(initializeRungZeroTimeState(
          config,
          options,
          state,
          runtime_units,
          background.has_value() ? &background.value() : nullptr,
          restoring_from_restart ? options.restart_state_override : nullptr));
      if (scheduler_initialization_reservation.valid()) {
        const auto baseline = memory_governor.snapshot().baseline_owned_bytes;
        if (baseline < state_bytes_before_scheduler) throw std::logic_error("scheduler initialization baseline is stale");
        const auto reports = std::array{core::collectSimulationMemoryReport(state),
            core::collectSchedulerMemoryReport(time_state_owner->particleScheduler(), time_state_owner->gasCellScheduler())};
        scheduler_initialization_reservation.reconcileBaselineOwnedAndRelease(core::checkedMemoryBytesAdd(
            baseline - state_bytes_before_scheduler, core::memoryReportBaselineOwnedBytes(core::mergeMemoryReports(reports)),
            "hierarchical scheduler initial baseline"));
      }
    } catch (...) { scheduler_initialization_failure = std::current_exception(); }
    failure_coordinator.rethrowCollectiveFailure(scheduler_initialization_failure, "scheduler initialization");
    RungZeroTimeState& time_state = time_state_owner.value();
    core::IntegratorState& integrator_state = time_state.integratorState();

    const std::filesystem::path zoom_region_path =
        config.mode.zoom_region_file.empty()
        ? std::filesystem::path{}
        : resolveConfigRelativePath(m_frozen_config, std::filesystem::path(config.mode.zoom_region_file));
    traceRuntimePhase("runtime_composition_begin");
    internal::ReferenceRuntimeComposition runtime_composition =
        internal::buildReferenceRuntimeComposition(
            internal::ReferenceRuntimeCompositionInputs{
                .frozen_config = m_frozen_config,
                .config = config,
                .state = state,
                .integrator_state = integrator_state,
                .mode_policy = mode_policy,
                .units = runtime_units,
                .services = runtime_services,
                .particle_scheduler = time_state.particleScheduler(),
                .gas_cell_scheduler = time_state.gasCellScheduler(),
                .report = report,
                .profiler = profiler,
                .pending_output = time_state.pendingOutput(),
                .options = options,
                .zoom_region_path = zoom_region_path,
                .world_rank = static_cast<std::uint32_t>(
                    std::max(mpi_context.worldRank(), 0)),
            });
    traceRuntimePhase("runtime_composition_complete");
    console_reporter.emitRuntimePhase(
        "runtime_composition_ready", "solver/runtime modules frozen and instantiated");
    if (restoring_from_restart) {
      const std::string& prior_runtime_state =
          options.restart_state_override->provenance.derived_runtime_state;
      const std::size_t cadence_pos = prior_runtime_state.find(
          "optional_diagnostic_cadence_policy=");
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "analysis.optional_cadence_restart_policy",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "analysis.diagnostics",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "optional science cadence starts fresh after restart; pending historical products are not replayed",
          .payload = {
              {"policy", "coalesced_nonpersistent"},
              {"prior_checkpoint_cadence", cadence_pos == std::string::npos
                  ? "legacy_unrecorded" : prior_runtime_state.substr(cadence_pos)},
          },
      });
    }
    GravityRuntime& gravity_callback = *runtime_composition.gravity;
    if (restoring_from_restart) {
      const io::RestartReadResult& restart = *options.restart_state_override;
      gravity_callback.restoreDecompositionEpoch(
          restart.distributed_gravity_state.decomposition_epoch);
      if (restart.gravity_force_cache.valid) {
        gravity_callback.importRestartForceCache(restart.gravity_force_cache, state);
      } else {
        // A pre-v20 checkpoint has no cache that can safely service the next
        // KDK pre-kick.  Do not pretend its persisted PM validity bit implies
        // a live force cache in this process; force the legal bootstrap path.
        integrator_state.pm_long_range_field_valid = false;
      }
    }
    report.treepm_pm_grid = gravity_callback.pmGridSize();
    report.treepm_pm_grid_nx = gravity_callback.pmGridShape().nx;
    report.treepm_pm_grid_ny = gravity_callback.pmGridShape().ny;
    report.treepm_pm_grid_nz = gravity_callback.pmGridShape().nz;
    report.treepm_pm_grid_shape = std::to_string(report.treepm_pm_grid_nx) + "x" +
        std::to_string(report.treepm_pm_grid_ny) + "x" + std::to_string(report.treepm_pm_grid_nz);
    report.treepm_update_cadence_steps = static_cast<int>(integrator_state.pm_sync_state.cadenceSteps());
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "gravity.treepm_setup",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "gravity.treepm",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "TreePM runtime configuration initialized",
        .payload = {
            {"pm_grid", std::to_string(config.numerics.treepm_pm_grid_nx) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_ny) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_nz)},
            {"pm_assignment_scheme", treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme)},
            {"pm_window_deconvolution", config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
            {"asmth_cells", formatRuntimeDouble(config.numerics.treepm_asmth_cells)},
            {"rcut_cells", formatRuntimeDouble(config.numerics.treepm_rcut_cells)},
            {"mesh_spacing_x_mpc_comoving", formatRuntimeDouble(config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx))},
            {"mesh_spacing_y_mpc_comoving", formatRuntimeDouble(config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny))},
            {"mesh_spacing_z_mpc_comoving", formatRuntimeDouble(config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))},
            {"split_scale_mpc_comoving", formatRuntimeDouble(config.numerics.treepm_asmth_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"cutoff_radius_mpc_comoving", formatRuntimeDouble(config.numerics.treepm_rcut_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"pm_update_cadence_steps", std::to_string(config.numerics.treepm_update_cadence_steps)},
            {"gravitational_constant_code", formatRuntimeDouble(gravity_callback.gravitationalConstantCode())},
            {"zoom_long_range_strategy",
                config.mode.zoom_long_range_strategy == core::ZoomLongRangeStrategy::kDisabled
                    ? "disabled"
                    : "global_coarse_plus_focused_highres_correction"},
            {"zoom_region_center_x_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_center_x_mpc_comoving)},
            {"zoom_region_center_y_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_center_y_mpc_comoving)},
            {"zoom_region_center_z_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_center_z_mpc_comoving)},
            {"zoom_region_radius_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_radius_mpc_comoving)},
            {"zoom_focused_pm_grid", std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_nz)},
            {"zoom_contamination_radius_mpc_comoving", formatRuntimeDouble(config.mode.zoom_contamination_radius_mpc_comoving)},
            {"softening_policy", describeGravitySofteningPolicy(config, state)},
            {"softening_kernel", "plummer"},
            {"softening_epsilon_kpc_comoving", formatRuntimeDouble(config.numerics.gravity_softening_kpc_comoving)},
            {"pm_fft_backend", gravity::PmSolver::fftBackendName()},
            {"pm_backend_capability", std::string(gravity::pmBackendCapabilityName(gravity::pmBackendCapability()))},
            {"gravity_acceptance_profile_id", "unverified_current_source"},
        },
    });
    console_reporter.emitPmSetup(
        report.treepm_pm_grid_nx,
        report.treepm_pm_grid_ny,
        report.treepm_pm_grid_nz,
        report.treepm_update_cadence_steps,
        gravity::PmSolver::fftBackendName());
    {
      const std::array startup_reports{
          core::collectSimulationMemoryReport(state),
          core::collectSchedulerMemoryReport(
              time_state.particleScheduler(), time_state.gasCellScheduler()),
          gravity_callback.memoryReport(),
          runtime_composition.hydro_amr->memoryReport(),
          runtime_composition.source->memoryReport(),
          profiler.retainedEventMemoryReport()};
      core::MemoryReport merged_startup_memory_report = core::mergeMemoryReports(startup_reports);
      const std::uint64_t governor_baseline_bytes =
          core::memoryReportBaselineOwnedBytes(merged_startup_memory_report);
      memory_governor.setBaselineOwnedBytes(governor_baseline_bytes);
      core::attachMemoryGovernorSnapshot(
          merged_startup_memory_report, memory_governor);
      core::attachProcessMemoryObservation(
          merged_startup_memory_report, core::observeProcessMemory());
      attachDistributedMemoryTelemetry(
          merged_startup_memory_report, runtime_services);
      profiler.setMemoryReport(std::move(merged_startup_memory_report));
      const core::MemoryReport* runtime_memory_report = profiler.memoryReport();
      if (runtime_memory_report != nullptr) {
        profiler.recordEvent(core::RuntimeEvent{
            .event_kind = "memory.runtime_startup_snapshot",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "core.memory",
            .step_index = integrator_state.step_index,
            .simulation_time_code = integrator_state.current_time_code,
            .scale_factor = integrator_state.current_scale_factor,
            .message = "startup memory accounting snapshot including gravity solver workspaces",
            .payload = {{"persistent_total_bytes", std::to_string(runtime_memory_report->totals.persistent_total_bytes)},
                        {"transient_total_bytes", std::to_string(runtime_memory_report->totals.transient_total_bytes)},
                        {"unknown_total_bytes", std::to_string(runtime_memory_report->totals.unknown_total_bytes)},
                        {"rank_max_owned_bytes", std::to_string(runtime_memory_report->distributed.rank_max_owned_bytes)},
                        {"global_sum_owned_bytes", std::to_string(runtime_memory_report->distributed.global_sum_owned_bytes)},
                        {"memory_imbalance_ratio", formatRuntimeDouble(runtime_memory_report->distributed.max_to_mean_imbalance_ratio)},
                        {"memory_governor_pressure",
                         runtime_memory_report->governor_snapshot.has_value()
                             ? std::string(core::memoryPressureLabel(
                                   runtime_memory_report->governor_snapshot->pressure))
                             : "unavailable"},
                        {"memory_governor_headroom_bytes",
                         runtime_memory_report->governor_snapshot.has_value()
                             ? std::to_string(
                                   runtime_memory_report->governor_snapshot->headroom_bytes)
                             : "0"}},
        });
        console_reporter.emitMemory(*runtime_memory_report, "startup");
      }
    }

    RuntimeConsoleStartupStatus startup_status;
    startup_status.run_name = config.output.run_name;
    startup_status.run_directory = report.shared_run_directory;
    startup_status.rank_directory = report.run_directory;
    startup_status.simulation_mode = core::modeToString(config.mode.mode);
    startup_status.mpi_world_size = mpi_context.worldSize();
    startup_status.openmp_compiled = omp_runtime.compiled;
    startup_status.openmp_threads = omp_runtime.configured_threads;
    startup_status.global_particle_count = report.global_particle_count;
    startup_status.global_cell_count = report.global_cell_count;
    startup_status.t_code = integrator_state.current_time_code;
    startup_status.dt_time_code = integrator_state.dt_time_code;
    if (mode_policy.cosmological_comoving_frame) {
      startup_status.a_scale = integrator_state.current_scale_factor;
      startup_status.redshift = consoleRedshift(
          true, integrator_state.current_scale_factor);
      startup_status.endpoint_a_scale = config.numerics.a_end;
    }
    startup_status.endpoint_t_code = config.numerics.t_code_end;
    startup_status.pm_grid_nx = report.treepm_pm_grid_nx;
    startup_status.pm_grid_ny = report.treepm_pm_grid_ny;
    startup_status.pm_grid_nz = report.treepm_pm_grid_nz;
    startup_status.pm_update_cadence_steps = report.treepm_update_cadence_steps;
    startup_status.snapshot_interval_steps = config.output.snapshot_interval_steps;
    startup_status.snapshot_interval_time_code =
        config.output.snapshot_interval_time_code;
    startup_status.write_restarts = options.write_outputs && config.output.write_restarts;
    startup_status.restoring_from_restart = restoring_from_restart;
    if (const core::MemoryReport* memory_report = profiler.memoryReport();
        memory_report != nullptr && memory_report->governor_snapshot.has_value()) {
      startup_status.memory_headroom_bytes =
          memory_report->governor_snapshot->headroom_bytes;
      startup_status.memory_pressure = std::string(core::memoryPressureLabel(
          memory_report->governor_snapshot->pressure));
    }
    console_reporter.emitStartup(startup_status);

    HydroAmrRuntime& hydro_callback = *runtime_composition.hydro_amr;
    auto effective_eos_table = runtime_composition.effective_eos_table;
    TimeCoordinator time_coordinator(
        runtime_services,
        time_state,
        gravity_callback,
        hydro_callback,
        *runtime_composition.source,
        std::move(runtime_composition.execution_plan),
        migration_balance,
        std::move(effective_eos_table));
    traceRuntimePhase("time_coordinator_begin");
    time_coordinator.runRungZeroSegment(
        config,
        options,
        state,
        background.has_value() ? &background.value() : nullptr,
        expected_global_particle_ids,
        report,
        profiler,
        mode_policy,
        restoring_from_restart);
    traceRuntimePhase("time_coordinator_complete");
    runtime_composition.output_restart->ensureFinalEndpointSnapshot(
        state, integrator_state);
    runtime_composition.analysis->finalizePending(
        integrator_state.step_index, integrator_state.current_time_code,
        integrator_state.current_scale_factor);

    report.final_state_digest = computeStateDigest(state, integrator_state, mpi_context);
    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    const std::uint64_t final_local_particle_id_square_sum = computeParticleIdSquareSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    const auto final_local_identity = parallel::summarizeLocalOwnedParticleIds(
        state.particle_sidecar.particle_id, state.validateUniqueParticleIds());
    report.local_particle_ids_unique = final_local_identity.local_particle_ids_unique;
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    const std::uint64_t final_global_particle_id_square_sum =
        mpi_context.allreduceSumUint64(final_local_particle_id_square_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    const bool final_all_ranks_have_unique_local_ids =
        mpi_context.allreduceSumUint64(report.local_particle_ids_unique ? 1ULL : 0ULL) ==
        static_cast<std::uint64_t>(mpi_context.worldSize());
    const parallel::LocalOwnershipIdentitySummary final_reduced_identity{
        .local_owned_count = report.global_particle_count,
        .local_particle_id_sum = report.global_particle_id_sum,
        .local_particle_id_square_sum = final_global_particle_id_square_sum,
        .local_particle_id_xor = report.global_particle_id_xor,
        .local_particle_ids_unique = final_all_ranks_have_unique_local_ids,
    };
    // Source terms can legitimately grow the authoritative particle population
    // after the initial decomposition (for example by spawning star particles).
    // The exact expected-ID lane is extended by TimeCoordinator at each birth
    // event, so derive the final scalar identity from that current authority
    // instead of comparing the evolved state against the immutable IC digest.
    const auto final_local_expected_identity =
        parallel::summarizeLocalOwnedParticleIds(expected_global_particle_ids, true);
    const parallel::LocalOwnershipIdentitySummary final_expected_global_identity{
        .local_owned_count = mpi_context.allreduceSumUint64(
            final_local_expected_identity.local_owned_count),
        .local_particle_id_sum = mpi_context.allreduceSumUint64(
            final_local_expected_identity.local_particle_id_sum),
        .local_particle_id_square_sum = mpi_context.allreduceSumUint64(
            final_local_expected_identity.local_particle_id_square_sum),
        .local_particle_id_xor = mpi_context.allreduceXorUint64(
            final_local_expected_identity.local_particle_id_xor),
        .local_particle_ids_unique = true,
    };
    std::vector<std::uint64_t> final_generic_local_owned_particle_ids;
    std::span<const std::uint64_t> final_local_owned_particle_ids = state.particle_sidecar.particle_id;
    if (!state.hasHomogeneousDmoMetadata()) {
      final_generic_local_owned_particle_ids.reserve(state.particles.size());
      for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
        if (state.particleOwningRank(particle_index) ==
            static_cast<std::uint32_t>(mpi_context.worldRank())) {
          final_generic_local_owned_particle_ids.push_back(
              state.particle_sidecar.particle_id[particle_index]);
        }
      }
      final_local_owned_particle_ids = final_generic_local_owned_particle_ids;
    }
    const parallel::ExactOwnershipPartitionReport exact_partition =
        parallel::validateExactGlobalOwnershipPartition(
            mpi_context, final_local_owned_particle_ids, expected_global_particle_ids,
            &memory_governor);
    report.global_particle_partition_identity_match = exact_partition.valid() &&
        parallel::partitionIdentityMatchesGeneratedSet(
            final_reduced_identity,
            final_expected_global_identity.local_owned_count,
            final_expected_global_identity.local_particle_id_sum,
            final_expected_global_identity.local_particle_id_square_sum,
            final_expected_global_identity.local_particle_id_xor);
    if (!report.global_particle_partition_identity_match) {
      throw std::runtime_error(
          "distributed ownership invariant failed at run completion: duplicate, missing, or extra authoritative particle IDs; "
          "actual_count=" + std::to_string(report.global_particle_count) +
          ", expected_count=" + std::to_string(final_expected_global_identity.local_owned_count) +
          ", actual_sum=" + std::to_string(report.global_particle_id_sum) +
          ", expected_sum=" + std::to_string(final_expected_global_identity.local_particle_id_sum) +
          ", actual_square_sum=" + std::to_string(final_global_particle_id_square_sum) +
          ", expected_square_sum=" + std::to_string(final_expected_global_identity.local_particle_id_square_sum) +
          ", actual_xor=" + std::to_string(report.global_particle_id_xor) +
          ", expected_xor=" + std::to_string(final_expected_global_identity.local_particle_id_xor) +
          ", duplicate_count=" + std::to_string(exact_partition.duplicate_count) +
          ", missing_count=" + std::to_string(exact_partition.missing_count) +
          ", extra_count=" + std::to_string(exact_partition.extra_count));
    }
    report.treepm_long_range_refresh_count = gravity_callback.longRangeRefreshCount();
    report.treepm_long_range_reuse_count = gravity_callback.longRangeReuseCount();
    report.multirate_rejected_pm_refreshes = gravity_callback.rejectedLongRangeRefreshCount();
    report.multirate_rejected_pm_reuses = gravity_callback.rejectedLongRangeReuseCount();
    report.multirate_rejected_cadence_records = gravity_callback.rejectedCadenceRecordCount();
    report.treepm_cadence_records.assign(
        gravity_callback.cadenceRecords().begin(),
        gravity_callback.cadenceRecords().end());
    report.canonical_stage_order = repeatedCanonicalOrder(report.stage_sequence);
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.complete",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "workflows.reference",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "runtime workflow completed",
        .payload = {{"completed_steps", std::to_string(report.completed_steps)},
                    {"final_time_code", formatRuntimeDouble(report.final_time_code)},
                    {"final_scale_factor", formatRuntimeDouble(report.final_scale_factor)},
                    {"run_directory", report.run_directory.string()},
                    {"treepm_update_cadence_steps", std::to_string(report.treepm_update_cadence_steps)},
                    {"treepm_long_range_refresh_count", std::to_string(report.treepm_long_range_refresh_count)},
                    {"treepm_long_range_reuse_count", std::to_string(report.treepm_long_range_reuse_count)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)},
                    {"local_particle_ids_unique", report.local_particle_ids_unique ? "true" : "false"},
                    {"global_particle_partition_identity_match",
                     report.global_particle_partition_identity_match ? "true" : "false"},
                    {"final_state_digest", std::to_string(report.final_state_digest)}},
    });
    traceRuntimePhase("final_artifact_flush_begin");
    flushCommonArtifacts(m_frozen_config, profiler, report);
    // Release cached graph communicators before the caller's MPI-session owner
    // can finalize MPI. This is the normal lifecycle path; destructors retain
    // only emergency containment for exception unwinding.
    gravity_callback.shutdownMpiResources();
    console_reporter.emitDone(
        report.completed_steps,
        report.final_time_code,
        mode_policy.cosmological_comoving_frame
            ? std::optional<double>(report.final_scale_factor)
            : std::nullopt,
        consoleRedshift(mode_policy.cosmological_comoving_frame, report.final_scale_factor),
        report.shared_run_directory,
        report.run_directory,
        report.normalized_config_snapshot_path,
        report.operational_report_json_path);
    traceRuntimePhase("run_complete");
    return report;
  } catch (const std::exception& ex) {
    console_reporter.emitWarning("workflows.reference", ex.what());
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.failure",
        .severity = core::RuntimeEventSeverity::kFatal,
        .subsystem = "workflows.reference",
        .step_index = std::nullopt,
        .simulation_time_code = std::nullopt,
        .scale_factor = std::nullopt,
        .message = "runtime workflow failed",
        .payload = {{"error", ex.what()}, {"run_directory", report.run_directory.string()}},
    });
    try {
      flushCommonArtifacts(m_frozen_config, profiler, report);
    } catch (...) {
    }
    throw;
  }
}

}  // namespace cosmosim::workflows
