#include "workflows/internal/reference_runtime_composition.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/periodic_domain.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/workflows/analysis_runtime.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/source_runtime.hpp"
#include "cosmosim/workflows/output_restart_runtime.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"

namespace cosmosim::workflows::internal {

struct CompositionAssembly {
  std::shared_ptr<AnalysisRuntime> analysis;
  std::shared_ptr<GravityRuntime> gravity;
  std::shared_ptr<HydroAmrRuntime> hydro_amr;
  std::shared_ptr<SourceRuntime> source;
  std::shared_ptr<OutputRestartRuntime> output_restart;
  std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table;
};

class DriftRuntime final {
 public:
  DriftRuntime(
      const RuntimeServices& services,
      bool periodic,
      double box_size_x,
      double box_size_y,
      double box_size_z) noexcept
      : m_services(services),
        m_periodic(periodic),
        m_box_size_x(box_size_x),
        m_box_size_y(box_size_y),
        m_box_size_z(box_size_z) {}

  void execute(DriftParticleStageView& view) const {
    view.requireFresh();
    core::StepContext& context = RuntimeStageAccess::driftContext(
        view,
        {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead}});
    if (context.stage != core::IntegrationStage::kDrift) {
      throw std::logic_error("drift task received a non-drift stage");
    }
    const double drift_factor = context.timeline_step.drift_factor_code;
    const std::uint32_t world_rank =
        static_cast<std::uint32_t>(m_services.mpi_context.worldRank());
    for (const std::uint32_t particle_index : context.active_set.particle_indices) {
      if (particle_index >= context.state.particles.size()) {
        throw std::out_of_range("drift task active particle index out of range");
      }
      if (context.state.particleOwningRank(particle_index) != world_rank) {
        continue;
      }
      // Gas-tagged generic particles are compatibility/lineage mirrors. The
      // authoritative Eulerian/AMR gas geometry lives in CellSoa and must not
      // acquire a second, independently drifted physical trajectory.
      if (context.state.particleSpeciesTag(particle_index) ==
          static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
        continue;
      }
      context.state.particles.position_x_comoving[particle_index] +=
          context.state.particles.velocity_x_peculiar[particle_index] * drift_factor;
      context.state.particles.position_y_comoving[particle_index] +=
          context.state.particles.velocity_y_peculiar[particle_index] * drift_factor;
      context.state.particles.position_z_comoving[particle_index] +=
          context.state.particles.velocity_z_peculiar[particle_index] * drift_factor;

      if (context.hierarchical_kdk.enabled) {
        // Canonical p=a*u remains fixed under the drift operator.
        const double drag = context.timeline_step.hubble_drag_factor;
        context.state.particles.velocity_x_peculiar[particle_index] *= drag;
        context.state.particles.velocity_y_peculiar[particle_index] *= drag;
        context.state.particles.velocity_z_peculiar[particle_index] *= drag;
      }

      // Canonical particle coordinates in a periodic cosmological domain must
      // remain in [0,L). Gravity/decomposition already use periodic-equivalent
      // coordinates internally; keep the authoritative state canonical too.
      if (m_periodic) {
        context.state.particles.position_x_comoving[particle_index] =
            core::canonicalPeriodicPosition(
                context.state.particles.position_x_comoving[particle_index],
                m_box_size_x);
        context.state.particles.position_y_comoving[particle_index] =
            core::canonicalPeriodicPosition(
                context.state.particles.position_y_comoving[particle_index],
                m_box_size_y);
        context.state.particles.position_z_comoving[particle_index] =
            core::canonicalPeriodicPosition(
                context.state.particles.position_z_comoving[particle_index],
                m_box_size_z);
      }
    }
    context.state.requireGasCellIdentityMapCoversDenseRows("drift task");
  }

 private:
  const RuntimeServices& m_services;
  bool m_periodic = false;
  double m_box_size_x = 0.0;
  double m_box_size_y = 0.0;
  double m_box_size_z = 0.0;
};

namespace {

[[nodiscard]] RuntimeResourceAccess read(RuntimeResourceKey resource) {
  return {.resource = resource, .mode = RuntimeResourceAccessMode::kRead};
}

[[nodiscard]] RuntimeResourceAccess write(RuntimeResourceKey resource) {
  return {.resource = resource, .mode = RuntimeResourceAccessMode::kWrite};
}

[[nodiscard]] RuntimeResourceAccess readWrite(RuntimeResourceKey resource) {
  return {.resource = resource, .mode = RuntimeResourceAccessMode::kReadWrite};
}

[[nodiscard]] RuntimeTaskSchedulingProfile schedulingProfile(
    std::uint64_t estimated_peak_bytes,
    RuntimeTaskPressureClass compute_pressure,
    RuntimeTaskPressureClass memory_bandwidth_pressure,
    RuntimeTaskPressureClass communication_pressure,
    bool optional = false,
    core::MemoryClass memory_class = core::MemoryClass::kPhaseResident,
    RuntimeTaskLifetimeBoundary lifetime_boundary = RuntimeTaskLifetimeBoundary::kTaskEnd) {
  return RuntimeTaskSchedulingProfile{
      .estimated_peak_bytes = estimated_peak_bytes,
      .peak_is_known = estimated_peak_bytes != 0U,
      .memory_class = memory_class,
      .lifetime_boundary = lifetime_boundary,
      .compute_pressure = compute_pressure,
      .memory_bandwidth_pressure = memory_bandwidth_pressure,
      .communication_pressure = communication_pressure,
      .optional = optional,
  };
}

[[nodiscard]] std::uint64_t declaredHydroTaskPeakBytes(const core::SimulationConfig& config) {
  const std::uint64_t cells = config.numerics.hydro_active_batch_max_cells == 0U
      ? static_cast<std::uint64_t>(hydro::k_hydro_automatic_active_batch_max_cells)
      : config.numerics.hydro_active_batch_max_cells;
  if (cells > std::numeric_limits<std::uint64_t>::max() /
                  hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell) {
    throw std::overflow_error("declared hydro task peak overflows uint64");
  }
  return cells * hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell;
}

[[nodiscard]] RuntimeModuleDescriptor makeAnalysisDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  static constexpr std::array stages{
      core::IntegrationStage::kGravityKickPre,
      core::IntegrationStage::kDrift,
      core::IntegrationStage::kForceRefresh,
      core::IntegrationStage::kHydroUpdate,
      core::IntegrationStage::kSourceTerms,
      core::IntegrationStage::kGravityKickPost,
      core::IntegrationStage::kAnalysisHooks,
      core::IntegrationStage::kOutputCheck,
  };
  std::vector<RuntimeTaskDeclaration> declarations;
  declarations.reserve(stages.size() + 1U);
  for (const core::IntegrationStage stage : stages) {
    declarations.push_back(RuntimeTaskDeclaration{
        .task_id = "analysis.stage_audit." + std::string(core::integrationStageName(stage)),
        .stage = stage,
        .ordinal = -100,
        .view_kind = RuntimeStageViewKind::kStageAudit,
        .resources = {write(RuntimeResourceKey::kDiagnostics)},
        .scheduling = [] {
          auto profile = schedulingProfile(
              0U, RuntimeTaskPressureClass::kLow, RuntimeTaskPressureClass::kLow,
              RuntimeTaskPressureClass::kLow);
          profile.peak_is_known = true;
          return profile;
        }(),
    });
  }
  declarations.push_back(RuntimeTaskDeclaration{
      .task_id = "analysis.diagnostics",
      .stage = core::IntegrationStage::kAnalysisHooks,
      .ordinal = 100,
      .view_kind = RuntimeStageViewKind::kAnalysis,
      .resources = {read(RuntimeResourceKey::kParticlePosition),
                    read(RuntimeResourceKey::kParticleVelocity),
                    read(RuntimeResourceKey::kParticleGravitySource),
                    read(RuntimeResourceKey::kHydroPrimitiveState),
                    read(RuntimeResourceKey::kSourceMutationState),
                    read(RuntimeResourceKey::kMigrationOwnership),
                    read(RuntimeResourceKey::kIntegratorTruth),
                    write(RuntimeResourceKey::kDiagnostics)},
      .dependencies = {"gravity::gravity.gravity_kick_post"},
      .scheduling = schedulingProfile(
          0U, RuntimeTaskPressureClass::kModerate, RuntimeTaskPressureClass::kHigh,
          RuntimeTaskPressureClass::kLow, true, core::MemoryClass::kDiagnostic),
  });

  return RuntimeModuleDescriptor{
      .module_id = "analysis",
      .schema_version = 1,
      .construction_ordinal = 10,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = std::move(declarations),
      .factory = [&config = inputs.config, &report = inputs.report, &inputs_state = inputs.state,
                  &inputs_integrator_state = inputs.integrator_state, assembly](
                     const RuntimeModuleFactoryContext& context) {
        std::shared_ptr<AnalysisRuntime> owner(
            makeAnalysisRuntime(config, report.stage_sequence, context.services));
        assembly->analysis = owner;
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        static constexpr std::array factory_stages{
            core::IntegrationStage::kGravityKickPre,
            core::IntegrationStage::kDrift,
            core::IntegrationStage::kForceRefresh,
            core::IntegrationStage::kHydroUpdate,
            core::IntegrationStage::kSourceTerms,
            core::IntegrationStage::kGravityKickPost,
            core::IntegrationStage::kAnalysisHooks,
            core::IntegrationStage::kOutputCheck,
        };
        for (const core::IntegrationStage stage : factory_stages) {
          instance.stage_tasks.push_back(RuntimeStageTaskContribution{
              .task_id = "analysis.stage_audit." +
                  std::string(core::integrationStageName(stage)),
              .task = StageAuditTask{[owner](AnalysisStageView& view) {
                owner->audit(view);
              }},
          });
        }
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "analysis.diagnostics",
            .estimate_incremental_bytes = [owner, &state = inputs_state,
                                           &integrator_state = inputs_integrator_state]() {
              if (integrator_state.step_index ==
                  std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("analysis completed-step index overflows uint64");
              }
              return owner->estimateRequiredIncrementalBytes(
                  state, integrator_state.step_index + 1U);
            },
            .task = AnalysisStageTask([owner](AnalysisStageView& view) {
              owner->executeDiagnostics(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeDriftDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs) {
  const bool periodic =
      inputs.mode_policy.gravity_boundary ==
      core::GravityBoundaryModel::kPeriodicPoisson;
  const double box_size_x =
      inputs.config.cosmology.box_size_x_mpc_comoving;
  const double box_size_y =
      inputs.config.cosmology.box_size_y_mpc_comoving;
  const double box_size_z =
      inputs.config.cosmology.box_size_z_mpc_comoving;

  return RuntimeModuleDescriptor{
      .module_id = "drift",
      .schema_version = 1,
      .construction_ordinal = 20,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "time.drift_particles",
          .stage = core::IntegrationStage::kDrift,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kDriftParticles,
          .resources = {readWrite(RuntimeResourceKey::kParticlePosition),
                        readWrite(RuntimeResourceKey::kParticleVelocity),
                        read(RuntimeResourceKey::kMigrationOwnership)},
          .dependencies = {"gravity::gravity.gravity_kick_pre"},
          .scheduling = schedulingProfile(
              0U, RuntimeTaskPressureClass::kModerate, RuntimeTaskPressureClass::kHigh,
              RuntimeTaskPressureClass::kLow),
      }},
      .factory = [periodic, box_size_x, box_size_y, box_size_z](
                     const RuntimeModuleFactoryContext& context) {
        auto owner = std::make_shared<DriftRuntime>(
            context.services,
            periodic,
            box_size_x,
            box_size_y,
            box_size_z);
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "time.drift_particles",
            .task = DriftStageTask([owner](DriftParticleStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeGravityDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  std::vector<RuntimeTaskDeclaration> declarations;
  for (const core::IntegrationStage stage : {
           core::IntegrationStage::kGravityKickPre,
           core::IntegrationStage::kForceRefresh,
           core::IntegrationStage::kGravityKickPost}) {
    declarations.push_back(RuntimeTaskDeclaration{
        .task_id = "gravity." + std::string(core::integrationStageName(stage)),
        .stage = stage,
        .ordinal = 0,
        .view_kind = RuntimeStageViewKind::kGravity,
        .resources = {read(RuntimeResourceKey::kParticlePosition),
                      readWrite(RuntimeResourceKey::kParticleVelocity),
                      read(RuntimeResourceKey::kParticleGravitySource),
                      readWrite(RuntimeResourceKey::kHydroConservedState),
                      readWrite(RuntimeResourceKey::kHydroPrimitiveState),
                      read(RuntimeResourceKey::kMigrationOwnership),
                      read(RuntimeResourceKey::kSchedulerTruth),
                      write(RuntimeResourceKey::kGravityAcceleration),
                      readWrite(RuntimeResourceKey::kIntegratorTruth)},
        .dependencies = stage == core::IntegrationStage::kForceRefresh
            ? std::vector<std::string>{"drift::time.drift_particles"}
            : (stage == core::IntegrationStage::kGravityKickPost
                   ? std::vector<std::string>{"sources::sources.update"}
                   : std::vector<std::string>{}),
        .scheduling = schedulingProfile(
            0U, RuntimeTaskPressureClass::kHigh, RuntimeTaskPressureClass::kHigh,
            RuntimeTaskPressureClass::kHigh),
    });
  }
  return RuntimeModuleDescriptor{
      .module_id = "gravity",
      .schema_version = 1,
      .construction_ordinal = 30,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = std::move(declarations),
      .factory = [inputs, assembly](const RuntimeModuleFactoryContext&) {
        assembly->gravity = std::shared_ptr<GravityRuntime>(makeGravityRuntime(
            inputs.config,
            inputs.mode_policy,
            inputs.services,
            inputs.zoom_region_path));
        RuntimeModuleInstance instance;
        instance.owner_lifetime = assembly->gravity;
        for (const core::IntegrationStage stage : {
                 core::IntegrationStage::kGravityKickPre,
                 core::IntegrationStage::kForceRefresh,
                 core::IntegrationStage::kGravityKickPost}) {
          instance.stage_tasks.push_back(RuntimeStageTaskContribution{
              .task_id = "gravity." + std::string(core::integrationStageName(stage)),
              .estimate_memory = [owner = assembly->gravity, &state = inputs.state]() {
                return owner->estimateMemory(state);
              },
              .task = GravityStageTask([owner = assembly->gravity](GravityStageView& view) {
                owner->execute(view);
              }),
          });
        }
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeHydroAmrDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  return RuntimeModuleDescriptor{
      .module_id = "hydro_amr",
      .schema_version = 1,
      .construction_ordinal = 40,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "hydro_amr.update",
          .stage = core::IntegrationStage::kHydroUpdate,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kHydroAmr,
          .resources = {readWrite(RuntimeResourceKey::kSourceMutationState),
                        readWrite(RuntimeResourceKey::kParticleVelocity),
                        readWrite(RuntimeResourceKey::kHydroConservedState),
                        readWrite(RuntimeResourceKey::kHydroPrimitiveState),
                        readWrite(RuntimeResourceKey::kAmrPatchState),
                        readWrite(RuntimeResourceKey::kEffectiveIsmThermodynamics),
                        read(RuntimeResourceKey::kMigrationOwnership),
                        read(RuntimeResourceKey::kIntegratorTruth)},
          .dependencies = {"gravity::gravity.force_refresh"},
          .scheduling = schedulingProfile(
              declaredHydroTaskPeakBytes(inputs.config), RuntimeTaskPressureClass::kHigh,
              RuntimeTaskPressureClass::kHigh, RuntimeTaskPressureClass::kHigh, false,
              core::MemoryClass::kScratchArena),
      }},
      .factory = [inputs, assembly](const RuntimeModuleFactoryContext&) {
        assembly->hydro_amr = std::shared_ptr<HydroAmrRuntime>(makeHydroAmrRuntime(
            inputs.config,
            inputs.mode_policy,
            inputs.services,
            assembly->effective_eos_table));
        RuntimeModuleInstance instance;
        instance.owner_lifetime = assembly->hydro_amr;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "hydro_amr.update",
            .estimate_memory = [owner = assembly->hydro_amr, &state = inputs.state]() {
              return owner->estimateMemory(state);
            },
            .task = HydroAmrStageTask([owner = assembly->hydro_amr](HydroAmrStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeSourceDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  return RuntimeModuleDescriptor{
      .module_id = "sources",
      .schema_version = 1,
      .construction_ordinal = 50,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "sources.update",
          .stage = core::IntegrationStage::kSourceTerms,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kSourceMutation,
          .resources = {readWrite(RuntimeResourceKey::kSourceMutationState),
                        readWrite(RuntimeResourceKey::kParticlePosition),
                        readWrite(RuntimeResourceKey::kParticleVelocity),
                        readWrite(RuntimeResourceKey::kParticleIdentity),
                        readWrite(RuntimeResourceKey::kParticleSpeciesIndex),
                        readWrite(RuntimeResourceKey::kHydroConservedState),
                        read(RuntimeResourceKey::kHydroPrimitiveState),
                        read(RuntimeResourceKey::kAmrPatchState),
                        read(RuntimeResourceKey::kEffectiveIsmThermodynamics),
                        readWrite(RuntimeResourceKey::kMigrationOwnership),
                        read(RuntimeResourceKey::kIntegratorTruth)},
          .dependencies = {"hydro_amr::hydro_amr.update"},
          .scheduling = schedulingProfile(
              0U, RuntimeTaskPressureClass::kHigh, RuntimeTaskPressureClass::kHigh,
              RuntimeTaskPressureClass::kModerate),
      }},
      .factory = [&config = inputs.config, &mode_policy = inputs.mode_policy,
                  &units = inputs.units, world_rank = inputs.world_rank,
                  &mpi_context = inputs.services.mpi_context,
                  &services = inputs.services, &state = inputs.state,
                  assembly](const RuntimeModuleFactoryContext&) {
        std::shared_ptr<SourceRuntime> owner(
            makeSourceRuntime(
                config,
                mode_policy,
                units,
                world_rank,
                mpi_context,
                assembly->effective_eos_table,
                &services));
        assembly->source = owner;
        if (assembly->hydro_amr != nullptr) {
          assembly->hydro_amr->bindSourceRuntime(owner.get());
        }
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "sources.update",
            .estimate_memory = [owner, &state]() {
              return owner->estimateMemory(state);
            },
            .task = SourceMutationStageTask([owner](SourceMutationStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeOutputRestartDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  return RuntimeModuleDescriptor{
      .module_id = "output_restart",
      .schema_version = 1,
      .construction_ordinal = 60,
      .prerequisites = {"gravity"},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "output_restart.boundary",
          .stage = core::IntegrationStage::kOutputCheck,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kOutputRestart,
          .resources = {read(RuntimeResourceKey::kParticlePosition),
                        read(RuntimeResourceKey::kParticleVelocity),
                        read(RuntimeResourceKey::kHydroConservedState),
                        read(RuntimeResourceKey::kSourceMutationState),
                        read(RuntimeResourceKey::kMigrationOwnership),
                        read(RuntimeResourceKey::kSchedulerTruth),
                        read(RuntimeResourceKey::kIntegratorTruth),
                        readWrite(RuntimeResourceKey::kOutputRestartState),
                        write(RuntimeResourceKey::kDiagnostics)},
          .dependencies = {"analysis::analysis.diagnostics"},
          .scheduling = schedulingProfile(
              0U, RuntimeTaskPressureClass::kLow, RuntimeTaskPressureClass::kHigh,
              RuntimeTaskPressureClass::kHigh),
      }},
      .factory = [inputs, assembly](const RuntimeModuleFactoryContext&) {
        if (!assembly->gravity) {
          throw std::logic_error("output/restart factory requires the gravity owner");
        }
        auto owner = std::make_shared<OutputRestartRuntime>(
            inputs.frozen_config,
            inputs.config,
            inputs.particle_scheduler,
            inputs.gas_cell_scheduler,
            *assembly->gravity,
            *assembly->analysis,
            inputs.services,
            inputs.report,
            inputs.profiler,
            inputs.pending_output,
            inputs.options.write_outputs);
        assembly->output_restart = owner;
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "output_restart.boundary",
            .estimate_memory = [owner, &state = inputs.state,
                                &integrator_state = inputs.integrator_state]() {
              return owner->estimateMemory(state, integrator_state);
            },
            .task = OutputRestartStageTask([owner](OutputRestartStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

}  // namespace

ReferenceRuntimeComposition buildReferenceRuntimeComposition(
    ReferenceRuntimeCompositionInputs inputs) {
  auto assembly = std::make_shared<CompositionAssembly>();
  if (inputs.config.physics.enable_star_formation &&
      inputs.config.physics.star_formation_model ==
          core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    assembly->effective_eos_table =
        std::make_shared<const physics::EffectiveMultiphaseEosTable>(
            physics::makeEffectiveMultiphaseEosConfig(inputs.config.physics),
            inputs.units,
            physics::makeEffectiveIsmReferenceCoolingProvider(inputs.config.physics));
  }
  RuntimeModuleRegistry registry;
  registry.registerModule(makeAnalysisDescriptor(inputs, assembly));
  registry.registerModule(makeDriftDescriptor(inputs));
  registry.registerModule(makeGravityDescriptor(inputs, assembly));
  registry.registerModule(makeHydroAmrDescriptor(inputs, assembly));
  registry.registerModule(makeSourceDescriptor(inputs, assembly));
  registry.registerModule(makeOutputRestartDescriptor(inputs, assembly));
  if (inputs.options.register_runtime_modules) {
    inputs.options.register_runtime_modules(registry);
  }
  RuntimeExecutionPlan execution_plan = registry.freezeAndInstantiate(
      RuntimeModuleFactoryContext{inputs.services});
  if (!assembly->gravity || !assembly->hydro_amr || !assembly->source ||
      !assembly->output_restart) {
    throw std::logic_error("reference runtime composition omitted a required owner");
  }
  return ReferenceRuntimeComposition{
      .analysis = std::move(assembly->analysis),
      .gravity = std::move(assembly->gravity),
      .hydro_amr = std::move(assembly->hydro_amr),
      .source = std::move(assembly->source),
      .output_restart = std::move(assembly->output_restart),
      .effective_eos_table = std::move(assembly->effective_eos_table),
      .execution_plan = std::move(execution_plan),
  };
}

}  // namespace cosmosim::workflows::internal
