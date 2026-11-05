#include "cosmosim/workflows/runtime_capabilities.hpp"

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/openmp_runtime.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/io/ic_reader.hpp"

#include <fstream>
#include <stdexcept>

namespace cosmosim::workflows {
namespace {

[[nodiscard]] std::string escapeJson(std::string_view value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (const char ch : value) {
    switch (ch) {
      case '\\':
        escaped += "\\\\";
        break;
      case '"':
        escaped += "\\\"";
        break;
      case '\n':
        escaped += "\\n";
        break;
      case '\r':
        escaped += "\\r";
        break;
      case '\t':
        escaped += "\\t";
        break;
      default:
        escaped.push_back(ch);
        break;
    }
  }
  return escaped;
}

}  // namespace

const RuntimeCapability& RuntimeCapabilityReport::require(
    std::string_view name) const {
  for (const RuntimeCapability& capability : capabilities) {
    if (capability.name == name) {
      return capability;
    }
  }
  throw std::out_of_range("runtime capability is absent from report: " +
                          std::string(name));
}

std::string_view runtimeCapabilityStatusName(
    RuntimeCapabilityStatus status) noexcept {
  switch (status) {
    case RuntimeCapabilityStatus::kSupported:
      return "supported";
    case RuntimeCapabilityStatus::kProvisional:
      return "provisional";
    case RuntimeCapabilityStatus::kUnsupported:
      return "unsupported";
  }
  return "unsupported";
}

std::string_view runtimeImplementationMaturityName(
    RuntimeImplementationMaturity maturity) noexcept {
  switch (maturity) {
    case RuntimeImplementationMaturity::kNotApplicable:
      return "not_applicable";
    case RuntimeImplementationMaturity::kProductionScalable:
      return "production_scalable";
    case RuntimeImplementationMaturity::kReference:
      return "reference";
  }
  return "not_applicable";
}

std::string_view runtimeScientificMaturityName(
    RuntimeScientificMaturity maturity) noexcept {
  switch (maturity) {
    case RuntimeScientificMaturity::kNotApplicable:
      return "not_applicable";
    case RuntimeScientificMaturity::kProvisional:
      return "provisional";
    case RuntimeScientificMaturity::kValidated:
      return "validated";
  }
  return "not_applicable";
}

std::string_view runtimeScalabilityClassName(
    RuntimeScalabilityClass scalability) noexcept {
  switch (scalability) {
    case RuntimeScalabilityClass::kNotApplicable:
      return "not_applicable";
    case RuntimeScalabilityClass::kLight:
      return "light";
    case RuntimeScalabilityClass::kModerate:
      return "moderate";
    case RuntimeScalabilityClass::kScalableFft:
      return "scalable_fft";
    case RuntimeScalabilityClass::kHeavyReference:
      return "heavy_reference";
  }
  return "not_applicable";
}

RuntimeCapabilityReport buildRuntimeCapabilityReport(
    const core::SimulationConfig& config) {
  const core::OpenMpRuntimeInfo omp = core::openMpRuntimeInfo();
  RuntimeCapabilityReport report;
  report.capabilities = {
      {.name = "production_treepm_pm_backend",
       .status = gravity::pmBackendProductionReady()
           ? RuntimeCapabilityStatus::kSupported
           : RuntimeCapabilityStatus::kUnsupported,
       .requested = config.numerics.gravity_solver == core::GravitySolver::kTreePm,
       .compiled = gravity::pmBackendCapability() != gravity::PmBackendCapability::kUnavailable,
       .dependency_available = gravity::pmBackendProductionReady(),
       .runtime_available = gravity::pmBackendProductionReady(),
       .active = config.numerics.gravity_solver == core::GravitySolver::kTreePm &&
           (gravity::pmBackendProductionReady() || config.numerics.treepm_allow_diagnostic_naive_dft),
       .implementation_maturity = gravity::pmBackendProductionReady()
           ? RuntimeImplementationMaturity::kProductionScalable
           : RuntimeImplementationMaturity::kReference,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .scalability = gravity::pmBackendProductionReady()
           ? RuntimeScalabilityClass::kScalableFft
           : RuntimeScalabilityClass::kHeavyReference,
       .detail = "TreePM PM backend capability=" +
           std::string(gravity::pmBackendCapabilityName(gravity::pmBackendCapability())) +
           (config.numerics.treepm_allow_diagnostic_naive_dft
                ? "; explicit diagnostic naive-DFT override is enabled."
                : "; production mode forbids diagnostic naive-DFT fallback.")},
      {.name = "fixed_global_timestep", .status = RuntimeCapabilityStatus::kSupported,
       .requested = true, .compiled = true, .dependency_available = true,
       .runtime_available = true, .active = config.numerics.hierarchical_max_rung == 0,
       .detail = "Production ReferenceWorkflow KDK with hierarchical_max_rung=0."},
      {.name = "adaptive_global_timestep", .status = RuntimeCapabilityStatus::kSupported,
       .requested = config.numerics.hierarchical_max_rung == 0, .compiled = true,
       .dependency_available = true, .runtime_available = true,
       .active = config.numerics.hierarchical_max_rung == 0,
       .detail = "TimeCoordinator resizes global KDK intervals using collective physical criteria and endpoint/output clipping; local hierarchical integration is a separate optional DMO capability."},
      {.name = "production_hierarchical_local_timestep", .status = RuntimeCapabilityStatus::kProvisional,
       .requested = config.numerics.hierarchical_max_rung != 0, .compiled = true,
       .dependency_available = true, .runtime_available = true,
       .active = config.numerics.hierarchical_max_rung != 0,
       .implementation_maturity = RuntimeImplementationMaturity::kReference,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .detail = "Source-implemented synchronized power-of-two DMO KDK with coarse PM endpoint kicks; validation pending, gas/source multirate integration unavailable."},
      {.name = "canonical_external_ic_import",
#if COSMOSIM_ENABLE_HDF5
       .status = RuntimeCapabilityStatus::kSupported, .requested = !config.mode.ic_file.empty(),
       .compiled = true, .dependency_available = true, .runtime_available = true,
       .active = !config.mode.ic_file.empty(),
       .detail = "Typed HDF5 IC import/conversion with strict audit manifest v" +
           std::to_string(io::kIcAuditManifestSchemaVersion) +
           ", source-chunk-to-canonical streaming, schema, and provenance validation is compiled in."},
#else
       .status = RuntimeCapabilityStatus::kUnsupported, .requested = !config.mode.ic_file.empty(),
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false, .detail = "This build has COSMOSIM_ENABLE_HDF5=OFF."},
#endif
      {.name = "multifile_external_ic_import",
#if COSMOSIM_ENABLE_HDF5
       .status = RuntimeCapabilityStatus::kSupported, .requested = false,
       .compiled = true, .dependency_available = true, .runtime_available = true,
       .active = false, .detail = "Validated multi-file HDF5 ingestion is compiled in."},
#else
       .status = RuntimeCapabilityStatus::kUnsupported, .requested = false,
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false, .detail = "This build has COSMOSIM_ENABLE_HDF5=OFF."},
#endif
      {.name = "distributed_ic_import",
#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI
       .status = RuntimeCapabilityStatus::kProvisional, .requested = config.parallel.mpi_ranks_expected > 1,
       .compiled = true, .dependency_available = true, .runtime_available = true,
       .active = config.parallel.mpi_ranks_expected > 1,
       .detail = "MPI+HDF5 distributed IC ingestion is implemented; production acceptance remains provisional pending the full rank matrix."},
#else
       .status = RuntimeCapabilityStatus::kUnsupported, .requested = config.parallel.mpi_ranks_expected > 1,
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false, .detail = "Distributed IC ingestion requires MPI and HDF5 in this build."},
#endif
      {.name = "production_hydro_cooling", .status = RuntimeCapabilityStatus::kSupported,
       .requested = config.physics.enable_cooling, .compiled = true,
       .dependency_available = true, .runtime_available = true,
       .active = config.physics.enable_cooling,
       .implementation_maturity = RuntimeImplementationMaturity::kProductionScalable,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .scalability = RuntimeScalabilityClass::kLight,
       .detail = "ReferenceWorkflow attaches the configured cooling/heating integrator to fixed-grid and production AMR hydro; physical density and hydrogen number density are derived explicitly from the stored density frame and unit system."},
      {.name = "distributed_gas_cell_hydro_ghosting",
#if COSMOSIM_ENABLE_MPI
       .status = RuntimeCapabilityStatus::kProvisional,
       .requested = config.parallel.mpi_ranks_expected > 1, .compiled = true,
       .dependency_available = true, .runtime_available = true,
       .active = config.parallel.mpi_ranks_expected > 1,
       .implementation_maturity = RuntimeImplementationMaturity::kProductionScalable,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .scalability = RuntimeScalabilityClass::kModerate,
       .detail = "Production distributed hydro uses the gas-cell protocol keyed by stable gas_cell_id, ownership/decomposition generation, and hydro synchronization epoch; the legacy fixed-Cartesian multi-rank fallback is rejected."},
#else
       .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = config.parallel.mpi_ranks_expected > 1, .compiled = false,
       .dependency_available = false, .runtime_available = false, .active = false,
       .detail = "Distributed gas-cell hydro ghosting requires an MPI-enabled build."},
#endif
      {.name = "black_hole_seeding_candidate_provider",
       .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = config.physics.bh_enable_seeding, .compiled = false,
       .dependency_available = false, .runtime_available = false, .active = false,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .detail = "ReferenceWorkflow has no authoritative halo/candidate provider for black-hole seeding and fails closed when bh_enable_seeding=true; pre-seeded black-hole accretion/feedback remains available independently."},
      {.name = "distributed_metal_diffusion",
#if COSMOSIM_ENABLE_MPI
       .status = config.physics.metal_species_mode == core::MetalSpeciesMode::kCoreElements
           ? RuntimeCapabilityStatus::kSupported
           : RuntimeCapabilityStatus::kUnsupported,
       .requested = config.physics.enable_metal_diffusion && config.parallel.mpi_ranks_expected > 1,
       .compiled = true, .dependency_available = true,
       .runtime_available = config.physics.metal_species_mode == core::MetalSpeciesMode::kCoreElements,
       .active = config.physics.enable_metal_diffusion && config.parallel.mpi_ranks_expected > 1 &&
           config.physics.metal_species_mode == core::MetalSpeciesMode::kCoreElements,
       .implementation_maturity = RuntimeImplementationMaturity::kProductionScalable,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .scalability = RuntimeScalabilityClass::kModerate,
       .detail = config.physics.metal_species_mode == core::MetalSpeciesMode::kCoreElements
           ? "Element-resolved MPI diffusion uses unique-owner remote faces, directed lane exchange, equal/opposite owner corrections, physical boundary topology and the configured explicit/RKL integration path."
           : "Distributed diffusion requires physics.metal_species_mode=core_elements; the total_only scalar layout has no element-resolved conservative remote-face contract."},
#else
       .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = config.physics.enable_metal_diffusion && config.parallel.mpi_ranks_expected > 1,
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false, .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .detail = "Distributed element-resolved metal diffusion requires an MPI-enabled build."},
#endif
      {.name = "rank_remappable_restart", .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = false, .compiled = false, .dependency_available = true,
       .runtime_available = false, .active = false,
       .detail = "Restart currently supports same-world-size rank-local continuation only."},
      {.name = "asynchronous_output", .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = false, .compiled = false, .dependency_available = true,
       .runtime_available = false, .active = false,
       .detail = "Snapshot/checkpoint writes are synchronous at restart-safe boundaries."},
      {.name = "openmp_execution",
       .status = omp.compiled ? RuntimeCapabilityStatus::kSupported : RuntimeCapabilityStatus::kUnsupported,
       .requested = config.parallel.omp_threads == 0 || config.parallel.omp_threads > 1,
       .compiled = omp.compiled, .dependency_available = omp.compiled,
       .runtime_available = omp.compiled, .active = omp.compiled && omp.configured_threads > 1,
       .detail = omp.compiled
           ? "OpenMP is compiled in; Tree gravity and production FFT analysis execute real shared-memory parallel loops; active_threads=" +
                 std::to_string(omp.configured_threads) + ", requested_threads=" +
                 std::to_string(omp.requested_threads) + "."
           : "This binary was built without OpenMP; thread requests above one are rejected."},
      {.name = "production_fof_halo_finder", .status = RuntimeCapabilityStatus::kSupported,
       .requested = config.analysis.halo_on_the_fly, .compiled = true,
       .dependency_available = true, .runtime_available = true,
       .active = config.analysis.halo_on_the_fly,
       .detail = "Production FOF candidate search uses a periodic cell-linked spatial hash; all-pairs remains only as an explicit small-N oracle."},
      {.name = "distributed_fof_halo_merge",
#if COSMOSIM_ENABLE_MPI
       .status = RuntimeCapabilityStatus::kProvisional,
       .requested = config.parallel.mpi_ranks_expected > 1 && config.analysis.halo_on_the_fly,
       .compiled = true, .dependency_available = true, .runtime_available = true,
       .active = config.parallel.mpi_ranks_expected > 1 && config.analysis.halo_on_the_fly,
       .detail = "Cross-rank FOF correctness uses a root-gather global spatial-hash union and broadcasts stable particle-ID membership labels; candidate search is scalable but root catalog assembly remains a distributed-scaling limit."},
#else
       .status = RuntimeCapabilityStatus::kUnsupported, .requested = false,
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false, .detail = "Distributed FOF merge requires an MPI-enabled build."},
#endif
      {.name = "bound_subhalo_finder", .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = false, .compiled = false, .dependency_available = true,
       .runtime_available = false, .active = false,
       .detail = "No physical bound-substructure finder is claimed; halo catalogs emit no fabricated subhalo entries."},
      {.name = "production_fft_power_spectrum", .status = RuntimeCapabilityStatus::kSupported,
       .requested = config.analysis.enable_diagnostics, .compiled = true,
       .dependency_available = true, .runtime_available = true,
       .eligible = config.analysis.diagnostics_execution_policy ==
               core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional,
       .scheduled = config.analysis.enable_diagnostics &&
           config.analysis.diagnostics_execution_policy ==
               core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional,
       .active = config.analysis.enable_diagnostics &&
           config.analysis.diagnostics_execution_policy ==
               core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional,
       .implementation_maturity = RuntimeImplementationMaturity::kProductionScalable,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .scalability = RuntimeScalabilityClass::kScalableFft,
       .detail =
#if COSMOSIM_ENABLE_FFTW
           "FFT-backed power spectrum is computationally production-scalable through FFTW; science validation remains provisional and scheduling requires the provisional-science policy; direct DFT remains an explicit reference oracle."},
#else
           "FFT-backed power spectrum is computationally production-scalable through the built-in radix-2 backend on power-of-two meshes; science validation remains provisional and scheduling requires the provisional-science policy; direct DFT remains an explicit reference oracle."},
#endif
      {.name = "production_diagnostics", .status = RuntimeCapabilityStatus::kSupported,
       .requested = config.analysis.enable_diagnostics, .compiled = true,
       .dependency_available = true, .runtime_available = true,
       .active = config.analysis.enable_diagnostics,
       .detail = "Run-health and scalable science diagnostics follow the typed execution policy."},
      {.name = "provisional_diagnostics",
       .status = config.analysis.diagnostics_execution_policy ==
                     core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional
           ? RuntimeCapabilityStatus::kProvisional
           : RuntimeCapabilityStatus::kUnsupported,
       .requested = config.analysis.diagnostics_execution_policy ==
                    core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional,
       .compiled = true, .dependency_available = true, .runtime_available = true,
       .active = config.analysis.diagnostics_execution_policy ==
                 core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional,
       .detail = "Only diagnostics explicitly marked provisional remain gated behind the opt-in policy."},
  };
  for (RuntimeCapability& capability : report.capabilities) {
    if (!capability.eligible && capability.runtime_available &&
        capability.status != RuntimeCapabilityStatus::kUnsupported &&
        capability.name != "production_fft_power_spectrum") {
      capability.eligible = true;
    }
    if (!capability.scheduled && capability.active &&
        capability.name != "production_fft_power_spectrum") {
      capability.scheduled = true;
    }
  }
  return report;
}

void validateRequestedRuntimeCapabilities(
    const core::SimulationConfig& config,
    const RuntimeCapabilityReport& report) {
  if (config.numerics.hierarchical_max_rung != 0 &&
      !report.require("production_hierarchical_local_timestep").runtime_available) {
    throw std::invalid_argument(
        "runtime capability production_hierarchical_local_timestep is unsupported: " +
        report.require("production_hierarchical_local_timestep").detail);
  }
  if (config.physics.bh_enable_seeding &&
      report.require("black_hole_seeding_candidate_provider").status !=
          RuntimeCapabilityStatus::kSupported) {
    throw std::invalid_argument(
        "runtime capability black_hole_seeding_candidate_provider is unsupported: " +
        report.require("black_hole_seeding_candidate_provider").detail);
  }
  if (config.physics.enable_metal_diffusion && config.parallel.mpi_ranks_expected > 1 &&
      report.require("distributed_metal_diffusion").status !=
          RuntimeCapabilityStatus::kSupported) {
    throw std::invalid_argument(
        "runtime capability distributed_metal_diffusion is unsupported: " +
        report.require("distributed_metal_diffusion").detail);
  }
  if (config.numerics.gravity_solver == core::GravitySolver::kTreePm) {
    gravity::requireTreePmSupportOrThrow(
        config.numerics.gravity_solver,
        config.numerics.treepm_allow_diagnostic_naive_dft);
  }
}

std::string serializeRuntimeCapabilityReportJson(
    const RuntimeCapabilityReport& report) {
  std::string json = "{\n  \"schema_version\": " +
      std::to_string(report.schema_version) + ",\n  \"capabilities\": [\n";
  for (std::size_t index = 0; index < report.capabilities.size(); ++index) {
    const RuntimeCapability& capability = report.capabilities[index];
    json += "    {\"name\": \"" + escapeJson(capability.name) +
        "\", \"status\": \"" +
        std::string(runtimeCapabilityStatusName(capability.status)) +
        "\", \"requested\": " + std::string(capability.requested ? "true" : "false") +
        ", \"compiled\": " + std::string(capability.compiled ? "true" : "false") +
        ", \"dependency_available\": " + std::string(capability.dependency_available ? "true" : "false") +
        ", \"runtime_available\": " + std::string(capability.runtime_available ? "true" : "false") +
        ", \"eligible\": " + std::string(capability.eligible ? "true" : "false") +
        ", \"scheduled\": " + std::string(capability.scheduled ? "true" : "false") +
        ", \"active\": " + std::string(capability.active ? "true" : "false") +
        ", \"implementation_maturity\": \"" +
        std::string(runtimeImplementationMaturityName(capability.implementation_maturity)) +
        "\", \"scientific_maturity\": \"" +
        std::string(runtimeScientificMaturityName(capability.scientific_maturity)) +
        "\", \"scalability\": \"" +
        std::string(runtimeScalabilityClassName(capability.scalability)) +
        "\", \"detail\": \"" + escapeJson(capability.detail) + "\"}";
    json += index + 1U == report.capabilities.size() ? "\n" : ",\n";
  }
  json += "  ]\n}\n";
  return json;
}

void writeRuntimeCapabilityReportJson(
    const RuntimeCapabilityReport& report,
    const std::filesystem::path& output_path) {
  std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error(
        "failed to open runtime capability report: " + output_path.string());
  }
  output << serializeRuntimeCapabilityReportJson(report);
  if (!output) {
    throw std::runtime_error(
        "failed to write runtime capability report: " + output_path.string());
  }
}

}  // namespace cosmosim::workflows
