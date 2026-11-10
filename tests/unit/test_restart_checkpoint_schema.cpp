#include <cassert>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

template <typename T>
concept ExposesParticleSoaStorageMember = requires(T payload) {
  payload.particle_soa_storage;
};

template <typename T>
concept ExposesTransientWorkspaceMember = requires(T payload) {
  payload.transient_workspace;
};

int main() {
  static_assert(
      std::is_same_v<decltype(cosmosim::io::RestartWritePayload{}.persistent_state.simulation_state), const cosmosim::core::SimulationState*>,
      "Restart payload must serialize canonical SimulationState ownership root");
  static_assert(
      !ExposesParticleSoaStorageMember<cosmosim::io::RestartWritePayload>,
      "Restart payload must not expose ParticleSoaStorage as persistent restart truth");
  static_assert(
      !ExposesTransientWorkspaceMember<cosmosim::io::RestartWritePayload>,
      "Restart payload must not expose transient workspace state");

  const auto& schema = cosmosim::io::restartSchema();
  assert(schema.name == "cosmosim_restart_v25");
  assert(schema.version == 25);
  assert(cosmosim::io::isRestartSchemaCompatible(24));
  assert(cosmosim::io::isRestartSchemaCompatible(23));
  assert(cosmosim::io::isRestartSchemaCompatible(22));
  assert(cosmosim::io::isRestartSchemaCompatible(21));
  assert(cosmosim::io::isRestartSchemaCompatible(20));
  assert(cosmosim::io::isRestartSchemaCompatible(15));
  assert(cosmosim::io::isRestartSchemaCompatible(14));
  assert(!cosmosim::io::isRestartSchemaCompatible(11));
  const auto& checklist = cosmosim::io::exactRestartCompletenessChecklist();
  assert(!checklist.empty());
  assert(checklist.front() == "simulation_state_lanes_and_metadata");
  bool saw_softening = false;
  bool saw_gas_identity = false;
  bool saw_star_formation_state = false;
  bool saw_hydro_geometry = false;
  bool saw_species_sidecars = false;
  bool saw_output_cadence = false;
  bool saw_stochastic_state = false;
  bool saw_restart_diagnostics = false;
  bool saw_amr_pending_flux = false;
  bool saw_amr_temporal_history = false;
  bool saw_gravity_force_cache = false;
  for (const std::string_view item : checklist) {
    saw_softening = saw_softening || item == "particle_identity_softening_and_drift_epoch_lanes";
    saw_gas_identity = saw_gas_identity || item == "gas_cell_identity_lanes";
    saw_star_formation_state =
        saw_star_formation_state || item == "gas_metal_mass_and_star_birth_identity_lanes";
    saw_hydro_geometry = saw_hydro_geometry || item == "hydro_geometry_patch_state";
    saw_species_sidecars = saw_species_sidecars || item == "species_specific_sidecars";
    saw_output_cadence = saw_output_cadence || item == "output_cadence_persistent_state";
    saw_stochastic_state = saw_stochastic_state || item == "stochastic_module_persistent_state";
    saw_restart_diagnostics = saw_restart_diagnostics || item == "restart_diagnostics_summary";
    saw_amr_pending_flux = saw_amr_pending_flux || item == "amr_pending_flux_register_state";
    saw_amr_temporal_history = saw_amr_temporal_history || item == "amr_temporal_boundary_history_state";
    saw_gravity_force_cache = saw_gravity_force_cache || item == "gravity_force_cache_at_kdk_boundary";
  }
  assert(saw_softening);
  assert(saw_gas_identity);
  assert(saw_star_formation_state);
  assert(saw_hydro_geometry);
  assert(saw_species_sidecars);
  assert(saw_output_cadence);
  assert(saw_stochastic_state);
  assert(saw_restart_diagnostics);
  assert(saw_amr_pending_flux);
  assert(saw_amr_temporal_history);
  assert(saw_gravity_force_cache);
  assert(checklist.back() == "payload_integrity_hash_and_hex");

  bool missing_fields_threw = false;
  try {
    cosmosim::io::RestartWritePayload incomplete;
    (void)cosmosim::io::restartPayloadIntegrityHash(incomplete);
  } catch (const std::invalid_argument&) {
    missing_fields_threw = true;
  }
  assert(missing_fields_threw);

  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(0);
  state.resizePatches(0);
  state.particle_sidecar.particle_id[0] = 1;
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  state.species.count_by_species[0] = 1;
  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());

  cosmosim::core::IntegratorState integrator_state;
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(1);
  scheduler.reset(1, 0, 0);

  cosmosim::io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "mode = toy\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "deadbeef", 0, payload.normalized_config_text);
  payload.distributed_gravity_state.schema_version = 2;
  payload.distributed_gravity_state.world_size = 1;
  payload.distributed_gravity_state.pm_grid_nx = 4;
  payload.distributed_gravity_state.pm_grid_ny = 4;
  payload.distributed_gravity_state.pm_grid_nz = 4;
  payload.distributed_gravity_state.owning_rank_by_item = {0};
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank = {0};
  payload.distributed_gravity_state.pm_slab_end_x_by_rank = {4};

  const std::uint64_t hash_before = cosmosim::io::restartPayloadIntegrityHash(payload);
  // Integrity traversal is allocation-free in the production payload path, so
  // even a tiny budget admits it: the result is byte-identical and nothing is
  // reserved or committed through the governor.
  cosmosim::core::MemoryGovernor tight_integrity({.hard_limit_bytes = 255U});
  assert(cosmosim::io::restartPayloadIntegrityHash(payload, &tight_integrity) == hash_before);
  assert(tight_integrity.snapshot().committed_bytes == 0U);
  assert(tight_integrity.snapshot().reserved_bytes == 0U);
  cosmosim::core::MemoryGovernor integrity_governor({.hard_limit_bytes = 1024U * 1024U});
  assert(cosmosim::io::restartPayloadIntegrityHash(payload, &integrity_governor) == hash_before);
  assert(integrity_governor.snapshot().committed_bytes == 0U);
  assert(integrity_governor.snapshot().reserved_bytes == 0U);

  integrator_state.time_bins.active_bin = 1;
  const std::uint64_t hash_after = cosmosim::io::restartPayloadIntegrityHash(payload);
  assert(hash_before != hash_after);

  integrator_state.time_bins.active_bin = 0;
  const std::uint64_t hash_restored = cosmosim::io::restartPayloadIntegrityHash(payload);
  payload.output_cadence_state.output_enabled = true;
  payload.output_cadence_state.write_restarts = true;
  payload.output_cadence_state.last_completed_step_index = integrator_state.step_index;
  payload.output_cadence_state.snapshot_interval_steps = 4;
  payload.output_cadence_state.next_snapshot_step_index = 4;
  payload.output_cadence_state.snapshot_interval_time_code = 0.25;
  payload.output_cadence_state.next_snapshot_time_code = 0.25;
  payload.output_cadence_state.snapshot_stem = "snapshot";
  payload.output_cadence_state.restart_stem = "restart";
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_restored);
  payload.output_cadence_state = {};
  payload.stochastic_state.modules.push_back(cosmosim::io::StochasticModulePersistentState{
      .module_name = "star_formation",
      .schema_version = 1,
      .rng_policy = "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)",
      .random_seed = 123456789ull,
      .rank_local_seed_offset = 2,
      .last_committed_step_index = integrator_state.step_index,
      .deterministic_from_serialized_inputs = true,
  });
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_restored);
  payload.stochastic_state = {};
  cosmosim::io::GravityForceCachePersistentState force_cache;
  force_cache.valid = true;
  force_cache.particle_id = {state.particle_sidecar.particle_id[0]};
  force_cache.particle_accel_x_comoving = {0.125};
  force_cache.particle_accel_y_comoving = {-0.25};
  force_cache.particle_accel_z_comoving = {0.5};
  payload.gravity_force_cache = &force_cache;
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_restored);
  payload.gravity_force_cache = nullptr;
  state.particle_sidecar.setGravitySofteningOverride(0, 0.125);
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_restored);
  const std::uint64_t hash_with_softening = cosmosim::io::restartPayloadIntegrityHash(payload);
  state.sidecars.upsert(cosmosim::core::ModuleSidecarBlock{
      .module_name = "hash_probe",
      .schema_version = 1,
      .payload = {std::byte{0xAA}}});
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_with_softening);
  const std::uint64_t hash_with_sidecar = cosmosim::io::restartPayloadIntegrityHash(payload);
  payload.distributed_gravity_state.decomposition_epoch = 1;
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_with_sidecar);
  payload.distributed_gravity_state.decomposition_epoch = 0;

  bool unsafe_half_step_payload_threw = false;
  try {
    cosmosim::core::IntegratorState half_step_integrator_state = integrator_state;
    half_step_integrator_state.inside_kdk_step = true;
    half_step_integrator_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
    half_step_integrator_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
    half_step_integrator_state.last_completed_restart_safe = true;
    cosmosim::io::RestartWritePayload unsafe_payload = payload;
    unsafe_payload.integrator_state = &half_step_integrator_state;
    (void)cosmosim::io::restartPayloadIntegrityHash(unsafe_payload);
  } catch (const std::runtime_error& ex) {
    const std::string message = ex.what();
    unsafe_half_step_payload_threw =
        message.find("inside_kdk_step=true") != std::string::npos &&
        message.find("half-step restart is not represented") != std::string::npos;
  }
  assert(unsafe_half_step_payload_threw);

  bool missing_metadata_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_hash_hex.clear();
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    missing_metadata_threw = true;
  }
  assert(missing_metadata_threw);

  bool provenance_hash_message_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.provenance.normalized_config_hash_hex.clear();
    invalid_metadata.provenance.config_hash_hex.clear();
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument& ex) {
    const std::string message = ex.what();
    provenance_hash_message_threw =
        message.find("provenance.normalized_config_hash_hex or provenance.config_hash_hex") != std::string::npos;
  }
  assert(provenance_hash_message_threw);

  bool mismatched_metadata_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_hash_hex = "0000000000000000";
    invalid_metadata.provenance =
        cosmosim::core::makeProvenanceRecord(invalid_metadata.normalized_config_hash_hex, "deadbeef", 0, invalid_metadata.normalized_config_text);
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    mismatched_metadata_threw = true;
  }
  assert(mismatched_metadata_threw);

  bool text_hash_mismatch_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_text = "mode = altered\n";
    invalid_metadata.normalized_config_hash_hex = payload.normalized_config_hash_hex;
    invalid_metadata.provenance = payload.provenance;
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    text_hash_mismatch_threw = true;
  }
  assert(text_hash_mismatch_threw);

  {
    cosmosim::core::SimulationState gas_state;
    gas_state.resizeParticles(4);
    gas_state.resizeCells(2);
    gas_state.resizePatches(1);
    gas_state.species.count_by_species = {2, 2, 0, 0, 0};
    for (std::size_t i = 0; i < gas_state.particles.size(); ++i) {
      gas_state.particle_sidecar.particle_id[i] = 900 + i;
      gas_state.particle_sidecar.owning_rank[i] = 0;
      gas_state.particle_sidecar.species_tag[i] =
          i < 2 ? static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)
                : static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    }
    gas_state.rebuildSpeciesIndex();
    gas_state.patches.patch_id[0] = 777;
    gas_state.patches.first_cell[0] = 0;
    gas_state.patches.cell_count[0] = 2;
    gas_state.patches.extent_x_comoving[0] = 2.0;
    gas_state.patches.extent_y_comoving[0] = 1.0;
    gas_state.patches.extent_z_comoving[0] = 1.0;
    gas_state.patches.cell_dim_x[0] = 2;
    gas_state.patches.cell_dim_y[0] = 1;
    gas_state.patches.cell_dim_z[0] = 1;
    gas_state.patches.owning_rank[0] = 0;
    gas_state.cells.patch_index[0] = 0;
    gas_state.cells.patch_index[1] = 0;
    gas_state.refreshGasCellIdentityFromParticleOrder();

    cosmosim::core::HierarchicalTimeBinScheduler gas_scheduler(2);
    gas_scheduler.reset(4, 0, 0);
    gas_scheduler.setElementBin(2, 1, 0);
    gas_scheduler.setElementBin(3, 2, 0);
    cosmosim::core::syncTimeBinMirrorsFromScheduler(
        gas_scheduler,
        gas_state,
        cosmosim::core::TimeBinMirrorDomain::kParticles);
    cosmosim::core::HierarchicalTimeBinScheduler gas_cell_scheduler(2);
    gas_cell_scheduler.reset(2, 0, 0);
    gas_cell_scheduler.setElementBin(0, 1, 0);
    gas_cell_scheduler.setElementBin(1, 2, 0);
    cosmosim::core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, gas_state);

    cosmosim::io::RestartWritePayload gas_payload = payload;
    gas_payload.persistent_state.simulation_state = &gas_state;
    gas_payload.scheduler = &gas_scheduler;
    gas_payload.gas_cell_scheduler = &gas_cell_scheduler;
    gas_payload.distributed_gravity_state.owning_rank_by_item = {0, 0, 0, 0};
    assert(cosmosim::io::restartPayloadIntegrityHash(gas_payload) != 0);

    gas_state.cells.time_bin[0] = 0;
    bool stale_cell_mirror_threw = false;
    try {
      (void)cosmosim::io::restartPayloadIntegrityHash(gas_payload);
    } catch (const std::invalid_argument&) {
      stale_cell_mirror_threw = true;
    }
    assert(stale_cell_mirror_threw);

    cosmosim::core::syncTimeBinMirrorsFromScheduler(
        gas_scheduler,
        gas_state,
        cosmosim::core::TimeBinMirrorDomain::kParticles);
    cosmosim::core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, gas_state);
    assert(cosmosim::io::restartPayloadIntegrityHash(gas_payload) != 0);

    gas_state.cells.patch_index[1] = 1;
    bool mismatched_patch_index_threw = false;
    try {
      (void)cosmosim::io::restartPayloadIntegrityHash(gas_payload);
    } catch (const std::invalid_argument& ex) {
      const std::string message = ex.what();
      mismatched_patch_index_threw =
          message.find("/state/cells/patch_index") != std::string::npos ||
          message.find("/state/patches") != std::string::npos;
    }
    assert(mismatched_patch_index_threw);

    cosmosim::core::SimulationState parentless_state;
    parentless_state.resizeParticles(0);
    parentless_state.resizeCells(1);
    parentless_state.resizePatches(1);
    parentless_state.patches.patch_id[0] = 901;
    parentless_state.patches.level[0] = 0;
    parentless_state.patches.first_cell[0] = 0;
    parentless_state.patches.cell_count[0] = 1;
    parentless_state.patches.extent_x_comoving[0] = 1.0;
    parentless_state.patches.extent_y_comoving[0] = 1.0;
    parentless_state.patches.extent_z_comoving[0] = 1.0;
    parentless_state.patches.cell_dim_x[0] = 1;
    parentless_state.patches.cell_dim_y[0] = 1;
    parentless_state.patches.cell_dim_z[0] = 1;
    parentless_state.patches.owning_rank[0] = 0;
    parentless_state.cells.patch_index[0] = 0;
    parentless_state.gas_cell_identity.assign({
        {.gas_cell_id = 7001, .parent_particle_id = std::nullopt, .owning_patch_id = 901, .local_cell_row = 0},
    });
    parentless_state.gas_cells.gas_cell_id[0] = 7001;
    parentless_state.gas_cells.parent_particle_id[0] = 0;
    parentless_state.gas_cells.velocity_x_peculiar[0] = 1.0;
    parentless_state.gas_cells.velocity_y_peculiar[0] = 2.0;
    parentless_state.gas_cells.velocity_z_peculiar[0] = 3.0;
    parentless_state.gas_cells.density_code[0] = 4.0;
    parentless_state.gas_cells.pressure_code[0] = 5.0;
    parentless_state.gas_cells.internal_energy_code[0] = 6.0;
    parentless_state.gas_cells.temperature_code[0] = 7.0;
    parentless_state.gas_cells.sound_speed_code[0] = 8.0;
    parentless_state.rebuildSpeciesIndex();
    assert(parentless_state.validateOwnershipInvariants());

    cosmosim::core::HierarchicalTimeBinScheduler parentless_scheduler(0);
    parentless_scheduler.reset(0, 0, 0);
    // Payloads with gas cells must carry the authoritative gas-cell scheduler.
    cosmosim::core::HierarchicalTimeBinScheduler parentless_gas_scheduler(0);
    parentless_gas_scheduler.reset(1, 0, 0);
    cosmosim::core::syncGasCellTimeBinMirrorsFromGasCellScheduler(parentless_gas_scheduler, parentless_state);
    cosmosim::io::RestartWritePayload parentless_payload = payload;
    parentless_payload.persistent_state.simulation_state = &parentless_state;
    parentless_payload.scheduler = &parentless_scheduler;
    parentless_payload.gas_cell_scheduler = &parentless_gas_scheduler;
    parentless_payload.distributed_gravity_state.owning_rank_by_item.clear();
    assert(cosmosim::io::restartPayloadIntegrityHash(parentless_payload) != 0);

    bool zero_parent_with_has_parent_threw = false;
    try {
      parentless_state.gas_cell_identity.assign({
          {.gas_cell_id = 7001, .parent_particle_id = 0, .owning_patch_id = 901, .local_cell_row = 0},
      });
    } catch (const std::invalid_argument& ex) {
      zero_parent_with_has_parent_threw =
          std::string(ex.what()).find("present parent_particle_id must be nonzero") != std::string::npos;
    }
    assert(zero_parent_with_has_parent_threw);
  }

  return 0;
}
