#include "cosmosim/workflows/runtime_capabilities.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/openmp_runtime.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/io/ic_reader.hpp"

#include <fstream>
#include <string>
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
  // CHUI-MHD-REF-v1 capabilities (pass 02). The synchronized CPU reference
  // path (leaf AMR, MPI, D/H/I/H/D, v26/v7 IO) is source-implemented and
  // UNVALIDATED; it is admitted only when every predicate of
  // mhdProductionDispatchRefusal() holds, otherwise the precise refusal is
  // reported. Status stays provisional, never certified.
  // See docs/architecture/chui_mhd_implementation_pass02_status.md.
  const bool mhd_requested = config.mhd.enabled;
  const bool mhd_nonideal_requested = config.mhd.enabled && config.mhd.nonideal.enabled;
  const std::string mhd_refusal =
      mhd_requested ? mhdProductionDispatchRefusal(config, core::ModePolicy{}) : std::string();
  const bool mhd_dispatch_available = mhd_requested && mhd_refusal.empty();
  report.capabilities.push_back(
      {.name = "mhd_ideal_reference_kernels",
       .status = RuntimeCapabilityStatus::kProvisional,
       .requested = mhd_requested, .compiled = true, .dependency_available = true,
       .runtime_available = true, .active = mhd_dispatch_available,
       .implementation_maturity = RuntimeImplementationMaturity::kReference,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .scalability = RuntimeScalabilityClass::kNotApplicable,
       .detail = "SSPRK2/PLM-MC/PPM/HLLD/HLLE/UCT-HLL/CT kernels in cosmosim_mhd "
                 "(source-implemented, UNVALIDATED, NOT_RUN)."});
  report.capabilities.push_back(
      {.name = "mhd_production_dispatch",
       .status = mhd_requested && !mhd_refusal.empty() ? RuntimeCapabilityStatus::kUnsupported
                                                       : RuntimeCapabilityStatus::kProvisional,
       .requested = mhd_requested, .compiled = true, .dependency_available = true,
       .runtime_available = !mhd_requested || mhd_dispatch_available,
       .active = mhd_dispatch_available,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .detail = mhd_requested && !mhd_refusal.empty()
           ? mhd_refusal
           : std::string("Synchronized CPU reference MHD in kHydroUpdate: leaf-AMR "
                         "same-level/coarse-fine ownership, MPI halo/EMF exchange, "
                         "global-step retry ownership, v26 restart and v7 snapshot "
                         "(source-implemented, UNVALIDATED, NOT_RUN).")});
  report.capabilities.push_back(
      {.name = "mhd_nonideal_transport",
       .status = RuntimeCapabilityStatus::kProvisional,
       .requested = mhd_nonideal_requested, .compiled = true, .dependency_available = true,
       .runtime_available = !mhd_nonideal_requested ||
           (mhd_dispatch_available &&
            config.mhd.nonideal.coefficient_provider !=
                core::MhdCoefficientProviderConfigMode::kChemistry),
       .active = mhd_nonideal_requested && mhd_dispatch_available,
       .implementation_maturity = RuntimeImplementationMaturity::kReference,
       .scientific_maturity = RuntimeScientificMaturity::kProvisional,
       .detail = "Ohmic/ambipolar (SSPRK2/RKL2) and Hall-HLL (SSPRK2) with D/H/I/H/D Strang "
                 "composition; constant and verified-table providers only (UNVALIDATED)."});
  report.capabilities.push_back(
      {.name = "mhd_chemistry_coefficient_provider",
       .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = mhd_nonideal_requested &&
           config.mhd.nonideal.coefficient_provider ==
               core::MhdCoefficientProviderConfigMode::kChemistry,
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false,
       .detail = "No validated charged-carrier/collision/grain contract exists in this "
                 "repository; chemistry-derived transport fails closed (amendment R5)."});
  report.capabilities.push_back(
      {.name = "mhd_distributed_recursive_amr",
       .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = mhd_requested && config.mhd.amr.distributed_subcycling,
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false,
       .detail = "Magnetic temporal face histories and coupled EMF subcycling are not "
                 "implemented; nonmagnetic recursive AMR is unaffected."});
  report.capabilities.push_back(
      {.name = "mhd_gpu_backend",
       .status = RuntimeCapabilityStatus::kUnsupported,
       .requested = mhd_requested && config.mhd.performance.optional_gpu_backend,
       .compiled = false, .dependency_available = false, .runtime_available = false,
       .active = false,
       .detail = "No GPU MHD kernels exist; there is no silent CPU fallback."});

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

std::string mhdProductionDispatchRefusal(const core::SimulationConfig& config,
                                         const core::ModePolicy& mode_policy) {
  (void)mode_policy;  // periodic, open and reflecting (conducting) walls are all supported
  if (!config.mhd.enabled) {
    return {};
  }
  const core::MhdConfig& mhd = config.mhd;
  if (mhd.performance.optional_gpu_backend) {
    return "mhd_performance.optional_gpu_backend=true: no GPU MHD kernels exist and there is "
           "no silent CPU fallback";
  }
  if (mhd.amr.distributed_subcycling ||
      config.numerics.amr_time_stepping == core::AmrTimeSteppingConfigMode::kRecursiveSubcycling) {
    return "magnetized recursive/distributed AMR subcycling (tier mhd_ref_amr_subcycled) is not "
           "implemented: magnetic temporal face histories and coupled EMF subcycling do not "
           "exist; use numerics.amr_time_stepping=synchronized and "
           "mhd_amr.distributed_subcycling=false";
  }
  if (mhd.nonideal.enabled &&
      mhd.nonideal.coefficient_provider == core::MhdCoefficientProviderConfigMode::kChemistry) {
    return "mhd_nonideal.coefficient_provider=chemistry: no validated charged-carrier/"
           "collision/grain contract exists (amendment R5); use constant or a verified table";
  }
  if (config.numerics.hierarchical_max_rung != 0) {
    return "numerics.hierarchical_max_rung must be 0: the MHD reference advances every leaf "
           "with one synchronized global interval";
  }
  if (config.numerics.global_step_max_retries < 1) {
    return "numerics.global_step_max_retries must be >= 1: the global-step transaction is the "
           "owner of MHD interval retries (the MHD stepper never commits a shortened interval "
           "on its own)";
  }
  if (config.numerics.amr_regrid_mode != core::AmrRegridConfigMode::kOff) {
    return "numerics.amr_regrid_mode must be off: magnetic face transfer is not wired into "
           "amr_production_regrid in this source pass (static leaf hierarchies only)";
  }
  if (config.parallel.mpi_ranks_expected > 1 &&
      config.parallel.decomposition_runtime_rebalance_enabled) {
    return "parallel.decomposition_runtime_rebalance_enabled must be false for multi-rank MHD: "
           "patch migration with magnetic face payloads is not wired in this source pass";
  }
  if (config.numerics.hydro_shadow_enabled) {
    return "numerics.hydro_shadow_enabled must be false: the hydro shadow reference compares "
           "nonmagnetic Godunov updates";
  }
#if !COSMOSIM_ENABLE_HDF5
  return "MHD requires an HDF5-enabled build: cosmosim_restart_v26 and snapshot v7 are the only "
         "magnetic restart/output authorities";
#else
  return {};
#endif
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
  if (config.mhd.enabled) {
    // Capability-specific MHD admission (pass 02): every unsupported
    // topology/IO/provider/backend/source/EOS combination is refused before
    // any initial-condition materialization with its precise reason.
    const auto reject_mhd = [&](std::string_view name) {
      const RuntimeCapability& capability = report.require(name);
      if (capability.requested && !capability.runtime_available) {
        throw std::invalid_argument(
            "MHD_UNSUPPORTED_PHYSICS: runtime capability " + std::string(name) +
            " is unavailable: " + capability.detail);
      }
    };
    reject_mhd("mhd_gpu_backend");
    reject_mhd("mhd_distributed_recursive_amr");
    reject_mhd("mhd_chemistry_coefficient_provider");
    reject_mhd("mhd_nonideal_transport");
    reject_mhd("mhd_production_dispatch");
    const std::string refusal = mhdProductionDispatchRefusal(config, core::ModePolicy{});
    if (!refusal.empty()) {
      throw std::invalid_argument("MHD_UNSUPPORTED_PHYSICS: " + refusal);
    }
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
