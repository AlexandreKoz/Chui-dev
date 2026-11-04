#include "cosmosim/workflows/analysis_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/analysis/zoom_contamination.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/optional_diagnostic_cadence.hpp"

namespace cosmosim::workflows {
namespace {

class AnalysisRuntimeImpl final : public AnalysisRuntime {
 public:
  AnalysisRuntimeImpl(
      const core::SimulationConfig& config,
      std::vector<std::string>& stage_sequence,
      const RuntimeServices& services)
      : m_config(config),
        m_stage_sequence(&stage_sequence),
        m_services(services),
        m_diagnostics(config, services.memory_governor) {}

  [[nodiscard]] std::uint64_t estimateRequiredIncrementalBytes(
      const core::SimulationState& state, std::uint64_t completed_step) const override {
    if (!m_config.analysis.enable_diagnostics ||
        completed_step % static_cast<std::uint64_t>(
            m_config.analysis.run_health_interval_steps) != 0U) {
      return 0U;
    }
    return m_diagnostics.estimateBundleIncrementalBytes(
        analysis::DiagnosticClass::kRunHealth, state);
  }

  void audit(AnalysisStageView& view) override {
    view.requireFresh();
    const core::StepContext& context = internal::RuntimeStageAccess::analysisContext(
        view,
        {{RuntimeResourceKey::kDiagnostics, RuntimeResourceAccessMode::kWrite}});
    const auto canonical = core::StageScheduler::kickDriftKickOrder();
    if (context.stage != canonical[m_audited_stage_count % canonical.size()]) {
      throw std::logic_error("runtime audit observed a noncanonical KDK stage order");
    }
    // Validate every stage, retain only sixteen recent complete step traces.
    // Lifetime stage totals live in the existing counter registry.
    constexpr std::size_t k_recent_stage_limit = 128U;
    if (m_stage_sequence->size() == k_recent_stage_limit) {
      m_stage_sequence->erase(m_stage_sequence->begin(),
          m_stage_sequence->begin() + static_cast<std::ptrdiff_t>(canonical.size()));
    }
    ++m_audited_stage_count;
    m_services.profiler.counters().addCount("runtime.audited_stage_count");
    m_stage_sequence->push_back(
        std::string(core::integrationStageName(context.stage)));
  }

  // G3.4: operational zoom-contamination report at the run-health cadence
  // (zoom_in mode). Collisionless dark-matter particles are classified by
  // mass against the GLOBAL minimum particle mass (relative tolerance
  // k_high_res_mass_tolerance, a mass-equality tolerance for identical IC
  // particle masses, not a physics parameter); rank partials combine by
  // sums and a minimum. Emitted as profiler event analysis.zoom_contamination
  // and counters; positions/radius use the comoving code length convention of
  // the TreePM zoom options (config values in Mpc comoving).
  void reportZoomContamination(const core::StepContext& context, std::uint64_t step) {
    if (m_config.mode.mode != core::SimulationMode::kZoomIn) return;
    const double radius = m_config.mode.zoom_contamination_radius_mpc_comoving > 0.0
        ? m_config.mode.zoom_contamination_radius_mpc_comoving
        : m_config.mode.zoom_region_radius_mpc_comoving;
    if (!(radius > 0.0)) return;  // rank-uniform (configuration)
    constexpr double k_high_res_mass_tolerance = 1.0e-6;
    const core::SimulationState& state = context.state;
    std::vector<double> x, y, z, mass;
    double local_min_mass = std::numeric_limits<double>::infinity();
    std::exception_ptr failure;
    try {
      const bool homogeneous = state.hasHomogeneousDmoMetadata();
      for (std::size_t p = 0; p < state.particles.size(); ++p) {
        const bool dark_matter = homogeneous ||
            (p < state.particle_sidecar.species_tag.size() &&
             state.particle_sidecar.species_tag[p] == static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter));
        if (!dark_matter) continue;
        x.push_back(state.particles.position_x_comoving[p]);
        y.push_back(state.particles.position_y_comoving[p]);
        z.push_back(state.particles.position_z_comoving[p]);
        mass.push_back(state.particles.mass_code[p]);
        local_min_mass = std::min(local_min_mass, state.particles.mass_code[p]);
      }
    } catch (...) {
      failure = std::current_exception();
    }
    FailureCoordinator(m_services).rethrowCollectiveFailure(failure, "zoom contamination census");
    const double global_min_mass = m_services.mpi_context.allreduceMinDouble(local_min_mass);
    analysis::ZoomContaminationReport report;
    try {
      if (std::isfinite(global_min_mass)) {
        report = analysis::zoomContaminationDiagnostics(
            x, y, z, mass,
            {m_config.mode.zoom_region_center_x_mpc_comoving, m_config.mode.zoom_region_center_y_mpc_comoving,
             m_config.mode.zoom_region_center_z_mpc_comoving},
            radius, global_min_mass, k_high_res_mass_tolerance, true,
            {m_config.cosmology.box_size_x_mpc_comoving, m_config.cosmology.box_size_y_mpc_comoving,
             m_config.cosmology.box_size_z_mpc_comoving});
      }
    } catch (...) {
      failure = std::current_exception();
    }
    FailureCoordinator(m_services).rethrowCollectiveFailure(failure, "zoom contamination diagnostic");
    std::array<std::uint64_t, 2> counts{report.high_res_count, report.low_res_count};
    m_services.mpi_context.allreduceSumUint64sInPlace(counts);
    std::array<double, 2> masses{report.high_res_mass, report.low_res_mass};
    m_services.mpi_context.allreduceSumDoublesInPlace(masses);
    const double nearest = m_services.mpi_context.allreduceMinDouble(report.nearest_low_res_distance);
    analysis::ZoomContaminationReport global;
    global.high_res_count = counts[0];
    global.low_res_count = counts[1];
    global.high_res_mass = masses[0];
    global.low_res_mass = masses[1];
    global.nearest_low_res_distance = nearest;
    m_services.profiler.counters().setCount("analysis.zoom.low_res_count", global.low_res_count);
    m_services.profiler.counters().setCount("analysis.zoom.high_res_count", global.high_res_count);
    m_services.profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "analysis.zoom_contamination",
        .severity = global.low_res_count > 0U ? core::RuntimeEventSeverity::kWarning : core::RuntimeEventSeverity::kInfo,
        .subsystem = "analysis.zoom",
        .step_index = step,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.timeline_step.scale_factor_end,
        .message = "measured low-resolution contamination inside the zoom contamination sphere",
        .payload = {{"radius_comoving", std::to_string(radius)},
                    {"high_res_count", std::to_string(global.high_res_count)},
                    {"low_res_count", std::to_string(global.low_res_count)},
                    {"high_res_mass_code", std::to_string(global.high_res_mass)},
                    {"low_res_mass_code", std::to_string(global.low_res_mass)},
                    {"low_res_mass_fraction", std::to_string(global.lowResMassFraction())},
                    {"nearest_low_res_distance_comoving", std::isfinite(nearest) ? std::to_string(nearest) : "inf"}}});
  }

  void executeDiagnostics(AnalysisStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::analysisContext(
        view,
        {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kParticleGravitySource, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kDiagnostics, RuntimeResourceAccessMode::kWrite}});
    if (context.stage != core::IntegrationStage::kAnalysisHooks) {
      throw std::logic_error("analysis diagnostics handler received an unregistered stage");
    }
    if (!m_config.analysis.enable_diagnostics) {
      return;
    }
    if (context.integrator_state.step_index ==
        std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("analysis completed-step index overflows uint64");
    }
    // AnalysisHooks observes the fully kicked end-of-step state before the
    // integrator commits that epoch. Label cadence and science products by the
    // state being measured, not by the still-committed step-begin metadata.
    const std::uint64_t step = context.integrator_state.step_index + 1U;
    const double scale_factor = context.timeline_step.scale_factor_end;
    if (!std::isfinite(scale_factor) || scale_factor <= 0.0) {
      throw std::runtime_error(
          "analysis completed-state scale factor must be finite and positive");
    }
    const auto run = [&](analysis::DiagnosticClass diagnostic_class,
                         core::MemoryReservation* reservation = nullptr) {
      analysis::DiagnosticsBundle bundle;
      std::exception_ptr generation_failure;
      try {
        bundle = m_diagnostics.generateBundle(
            context.state, step, scale_factor, diagnostic_class,
            context.workspace, reservation,
            m_services.mpi_context.worldSize() == 1);
      } catch (...) {
        generation_failure = std::current_exception();
      }
      // No rank may enter diagnostics collectives after a peer failed local
      // bundle construction. This is the same prepare-then-collective rule used
      // by distributed I/O.
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          generation_failure, "distributed diagnostics local generation");

      analysis::reduceDiagnosticsBundleAcrossRanks(
          bundle, m_services.mpi_context);
      attachDistributedMemoryTelemetry(bundle.memory_report, m_services);

      // Shared ordinary diagnostics have one filesystem owner. Every rank
      // contributes to the global bundle; only rank zero publishes it.
      std::exception_ptr publication_failure;
      if (m_services.mpi_context.isRoot()) {
        try {
          m_diagnostics.writeBundle(bundle);
        } catch (...) {
          publication_failure = std::current_exception();
        }
      }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          publication_failure, "distributed diagnostics publication");
    };
    if (step % static_cast<std::uint64_t>(
                   m_config.analysis.run_health_interval_steps) == 0) {
      reportZoomContamination(context, step);
      // Required health is never silently skipped. Its owner admission failure
      // is a controlled error, coordinated before any later collective phase.
      std::exception_ptr health_failure;
      core::MemoryReservation health_reservation;
      try {
        if (m_services.memory_governor != nullptr) {
          health_reservation = m_services.memory_governor->reserve(
              core::MemoryClass::kDiagnostic,
              m_diagnostics.estimateBundleIncrementalBytes(
                  analysis::DiagnosticClass::kRunHealth, context.state),
              "analysis.run_health");
          health_reservation.commit();
        }
      } catch (...) { health_failure = std::current_exception(); }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          health_failure, "required analysis memory admission");
      std::exception_ptr health_execution_failure;
      try {
        run(analysis::DiagnosticClass::kRunHealth, &health_reservation);
      } catch (...) {
        health_execution_failure = std::current_exception();
      }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          health_execution_failure, "required analysis execution");
    }

    const core::MemoryPressure memory_pressure =
        m_services.memory_governor != nullptr
        ? m_services.memory_governor->snapshot().pressure
        : core::MemoryPressure::kGreen;
    const bool defer_optional_analysis =
        memory_pressure == core::MemoryPressure::kRed ||
        memory_pressure == core::MemoryPressure::kTrip;
    const bool science_light_due_now =
        step % static_cast<std::uint64_t>(
                   m_config.analysis.science_light_interval_steps) == 0 &&
        m_config.analysis.diagnostics_execution_policy !=
            core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly;
    const bool science_heavy_due_now =
        step % static_cast<std::uint64_t>(
                   m_config.analysis.science_heavy_interval_steps) == 0 &&
        m_config.analysis.diagnostics_execution_policy ==
            core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
    const auto executeOptional = [&](analysis::DiagnosticClass diagnostic_class,
                                     internal::OptionalDiagnosticCadence& pending,
                                     bool due_now, const char* class_name) {
      if (!due_now && !pending.pending) { return; }
      const auto previous = pending;
      core::MemoryReservation reservation;
      std::exception_ptr admission_failure;
      bool memory_rejected = defer_optional_analysis;
      std::string rejection_reason = defer_optional_analysis
          ? "process_memory_pressure" : "";
      std::uint64_t requested_bytes = 0U;
      // Each rank must make the same decision. A local hard-limit rejection
      // is not an error in optional science; other failures remain fatal.
      if (!memory_rejected) {
        try {
          requested_bytes = m_diagnostics.estimateBundleIncrementalBytes(
              diagnostic_class, context.state);
          if (m_services.memory_governor != nullptr) {
            reservation = m_services.memory_governor->reserve(
                core::MemoryClass::kDiagnostic, requested_bytes,
                "analysis.optional_science");
            reservation.commit();
          }
        } catch (const core::MemoryAdmissionError& error) {
          memory_rejected = true;
          rejection_reason = error.what();
        } catch (...) {
          admission_failure = std::current_exception();
        }
      }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          admission_failure, "optional analysis memory preparation");
      const bool any_rejected =
          FailureCoordinator(m_services).failedRankCount(memory_rejected) != 0U;
      if (any_rejected) {
        reservation.release();
        if (due_now) { pending.recordDue(step); }
        m_services.profiler.recordEvent(core::RuntimeEvent{
            .event_kind = "analysis.memory_pressure_deferral",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "analysis.diagnostics",
            .step_index = step,
            .simulation_time_code = context.timeline_step.time_end_code,
            .scale_factor = scale_factor,
            .message = "optional science diagnostics deferred by collective memory admission",
            .payload = {
                {"pressure", std::string(core::memoryPressureLabel(memory_pressure))},
                {"diagnostic_class", class_name},
                {"requested_bytes", std::to_string(requested_bytes)},
                {"headroom_bytes", m_services.memory_governor != nullptr
                    ? std::to_string(m_services.memory_governor->snapshot().headroom_bytes) : "unlimited"},
                {"reason", memory_rejected ? rejection_reason : "peer_memory_rejection"},
                {"science_light_due", science_light_due_now ? "true" : "false"},
                {"science_heavy_due", science_heavy_due_now ? "true" : "false"},
                {"science_light_pending", m_science_light_pending.pending ? "true" : "false"},
                {"science_heavy_pending", m_science_heavy_pending.pending ? "true" : "false"},
                {"science_light_missed_count", std::to_string(m_science_light_pending.missed_count)},
                {"science_heavy_missed_count", std::to_string(m_science_heavy_pending.missed_count)},
            },
        });
        return;
      }
      // The physical reservation survives computation, bundle materialization,
      // serialization, and the final local write. Do not clear pending on error.
      std::exception_ptr execution_failure;
      try {
        run(diagnostic_class, &reservation);
      } catch (...) {
        execution_failure = std::current_exception();
      }
      FailureCoordinator(m_services).rethrowCollectiveFailure(
          execution_failure, "optional analysis execution");
      pending.clear();
      if (previous.pending) {
        m_services.profiler.recordEvent(core::RuntimeEvent{
            .event_kind = "analysis.memory_pressure_catchup",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "analysis.diagnostics",
            .step_index = step,
            .simulation_time_code = context.timeline_step.time_end_code,
            .scale_factor = scale_factor,
            .message = "coalesced optional science diagnostic completed at the current physical epoch",
            .payload = {
                {"diagnostic_class", class_name},
                {"first_due_step", std::to_string(previous.first_due_step)},
                {"latest_due_step", std::to_string(previous.latest_due_step)},
                {"missed_count", std::to_string(previous.missed_count)},
                {"coalesced_count", std::to_string(previous.missed_count - 1U)},
                {"actual_execution_step", std::to_string(step)},
                {"science_light_catchup", diagnostic_class == analysis::DiagnosticClass::kScienceLight ? "true" : "false"},
                {"science_heavy_catchup", diagnostic_class == analysis::DiagnosticClass::kScienceHeavy ? "true" : "false"},
                {"historical_state_replayed", "false"},
            },
        });
      }
    };
    executeOptional(analysis::DiagnosticClass::kScienceLight,
                    m_science_light_pending, science_light_due_now, "science_light");
    executeOptional(analysis::DiagnosticClass::kScienceHeavy,
                    m_science_heavy_pending, science_heavy_due_now, "science_heavy");
    if (m_services.mpi_context.isRoot()) {
      m_diagnostics.enforceRetentionPolicy();
    }
  }

  [[nodiscard]] std::string optionalCadenceProvenance() const override {
    const auto describe = [](const internal::OptionalDiagnosticCadence& state) {
      return std::string("pending=") + (state.pending ? "true" : "false") +
          ",first_due_step=" + std::to_string(state.first_due_step) +
          ",latest_due_step=" + std::to_string(state.latest_due_step) +
          ",missed_count=" + std::to_string(state.missed_count);
    };
    return "optional_diagnostic_cadence_policy=coalesced_nonpersistent\n"
           "optional_diagnostic_science_light=" + describe(m_science_light_pending) + "\n" +
           "optional_diagnostic_science_heavy=" + describe(m_science_heavy_pending) + "\n";
  }

  void finalizePending(std::uint64_t completed_step, double time_code,
                       double scale_factor) override {
    const auto drop = [&](internal::OptionalDiagnosticCadence& pending,
                          const char* class_name) {
      if (!pending.pending) { return; }
      m_services.profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "analysis.optional_cadence_dropped",
          .severity = core::RuntimeEventSeverity::kWarning,
          .subsystem = "analysis.diagnostics",
          .step_index = completed_step,
          .simulation_time_code = time_code,
          .scale_factor = scale_factor,
          .message = "pending optional science cadence discarded at run termination; no historical state replay",
          .payload = {
              {"diagnostic_class", class_name},
              {"first_due_step", std::to_string(pending.first_due_step)},
              {"latest_due_step", std::to_string(pending.latest_due_step)},
              {"dropped_count", std::to_string(pending.missed_count)},
              {"reason", "run_termination"},
          },
      });
      pending.clear();
    };
    drop(m_science_light_pending, "science_light");
    drop(m_science_heavy_pending, "science_heavy");
  }

 private:
  core::SimulationConfig m_config;
  std::vector<std::string>* m_stage_sequence = nullptr;
  std::uint64_t m_audited_stage_count = 0;
  const RuntimeServices& m_services;
  analysis::DiagnosticsEngine m_diagnostics;
  internal::OptionalDiagnosticCadence m_science_light_pending;
  internal::OptionalDiagnosticCadence m_science_heavy_pending;
};

}  // namespace

std::unique_ptr<AnalysisRuntime> makeAnalysisRuntime(
    const core::SimulationConfig& config,
    std::vector<std::string>& stage_sequence,
    const RuntimeServices& services) {
  return std::make_unique<AnalysisRuntimeImpl>(config, stage_sequence, services);
}

}  // namespace cosmosim::workflows
