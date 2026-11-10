#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <span>
#include <vector>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::core { class MemoryGovernor; struct FrozenConfig; }

namespace cosmosim::io {

struct RestartSchema {
  // v25 adds /state/amr_temporal_boundary_history/{start,end}_internal_energy_density_comoving
  // (group schema_version 3): the dual-energy auxiliary rho*e endpoints of
  // each coarse temporal-history cell, so a mid-subcycle restart forms the
  // same temporally interpolated auxiliary ghost lane as the uninterrupted
  // run (P1-55-01). Older checkpoints restore them as NaN (unknown).
  // v24 adds /gravity_force_cache/hierarchical_endpoint: the synchronized
  // hierarchical-KDK split force (Tree and PM components) and nearest-pair
  // encounter lane at the checkpoint, so a resumed run reuses the same
  // endpoint force as the uninterrupted run instead of re-evaluating it.
  // v23 adds canonical SHA-256 payload integrity while retaining the v22 legacy hash
  // only for backward compatibility and migration diagnostics.
  // v22 adds conserved gas metal mass and immutable star-birth identity.
  std::string name = "cosmosim_restart_v25";
  std::uint32_t version = 25;
};

[[nodiscard]] const RestartSchema& restartSchema();
[[nodiscard]] bool isRestartSchemaCompatible(std::uint32_t file_schema_version);
[[nodiscard]] const std::vector<std::string_view>& exactRestartCompletenessChecklist();

struct RestartWritePolicy {
  bool enable_fsync_finalize = true;
  std::string temporary_suffix = ".part";
  // Optional existing process authority for owner-local serialization scratch.
  core::MemoryGovernor* memory_governor = nullptr;
};

// Production restart reader controls: memory governance and dimension validation
// for the read side.  The writer already reserves owning-state capacity via
// RestartWritePolicy; the reader needs the same optional authority plus a flag
// that forces the strictest cross-rank dimension consistency check before any
// allocation is committed.
struct RestartReadPolicy {
  // Optional governor for read-side scratch (candidate buffers, dimension
  // preflight reservations).  When nullptr the reader falls back to
  // unmeasured heap allocation.
  core::MemoryGovernor* memory_governor = nullptr;
  // When true the reader verifies that per-rank particle/cell counts are
  // consistent with the restart's own global metadata and throws on
  // mismatch rather than silently truncating.
  bool enforce_dimension_consistency = true;
};

// Gate-5 physical readback ownership: preflight dimensions checked BEFORE any
// population allocation, admitted via policy.memory_governor as incremental
// candidate memory coexisting with the already-live canonical state. The live
// state is NOT double-charged; only the candidate is reserved.
struct RestartReadCandidateDimensions {
  std::uint64_t particle_count = 0U;
  std::uint64_t cell_count = 0U;
  std::uint64_t patch_count = 0U;
  std::uint64_t pending_flux_count = 0U;
  std::uint64_t temporal_cell_count = 0U;
};

[[nodiscard]] std::uint64_t restartReadCandidateStagingBytes(
    const RestartReadCandidateDimensions& dimensions);

struct RestartPersistentStateView {
  // Narrow restart-serialization root: persistent-only state ownership path.
  // TransientStepWorkspace/HydroScratch/PM/MPI/output scratch are intentionally
  // outside this view so restart traversal cannot accidentally serialize them.
  const core::SimulationState* simulation_state = nullptr;
};

struct OutputCadencePersistentState {
  bool output_enabled = false;
  bool write_restarts = false;
  bool snapshot_due = false;
  bool checkpoint_due = false;
  std::uint64_t last_completed_step_index = 0;
  std::uint64_t snapshot_interval_steps = 0;
  std::uint64_t next_snapshot_step_index = 0;
  double snapshot_interval_time_code = 0.0;
  double next_snapshot_time_code = 0.0;
  std::string snapshot_stem;
  std::string restart_stem;
};

struct StochasticModulePersistentState {
  std::string module_name;
  std::uint32_t schema_version = 1;
  std::string rng_policy;
  std::uint64_t random_seed = 0;
  std::uint32_t rank_local_seed_offset = 0;
  std::uint64_t last_committed_step_index = 0;
  bool deterministic_from_serialized_inputs = true;
};

struct StochasticPersistentState {
  std::vector<StochasticModulePersistentState> modules;
};

struct GravityForceCachePersistentView {
  // Non-owning write/verification view. The live SimulationState and
  // GravityRuntime remain the sole owners of these lanes.
  bool valid = false;
  std::span<const std::uint64_t> particle_id;
  std::span<const std::uint64_t> gas_cell_id;
  std::span<const double> particle_accel_x_comoving;
  std::span<const double> particle_accel_y_comoving;
  std::span<const double> particle_accel_z_comoving;
  std::span<const double> cell_accel_x_comoving;
  std::span<const double> cell_accel_y_comoving;
  std::span<const double> cell_accel_z_comoving;
  // v24 hierarchical synchronized endpoint (rows follow particle_id). Valid
  // only together with `valid`, at a hierarchical synchronization boundary
  // whose Tree/PM split force was evaluated at the checkpoint positions.
  bool hierarchical_endpoint_valid = false;
  std::uint64_t hierarchical_pm_field_version = 0U;
  std::span<const double> hierarchical_tree_accel_x_comoving;
  std::span<const double> hierarchical_tree_accel_y_comoving;
  std::span<const double> hierarchical_tree_accel_z_comoving;
  std::span<const double> hierarchical_pm_accel_x_comoving;
  std::span<const double> hierarchical_pm_accel_y_comoving;
  std::span<const double> hierarchical_pm_accel_z_comoving;
  // Squared comoving nearest direct-pair separation; +infinity when none.
  std::span<const double> hierarchical_nearest_pair_r2_comoving;
};

struct GravityForceCachePersistentState {
  // Cached accelerations are restart-authoritative at a safe KDK boundary:
  // the next pre-kick consumes them before ForceRefresh rebuilds a new field.
  // The identity lanes make the cache safe to import after a dense-row reorder,
  // compaction, or local migration: cache order is never physical identity.
  bool valid = false;
  std::vector<std::uint64_t> particle_id;
  std::vector<std::uint64_t> gas_cell_id;
  std::vector<double> particle_accel_x_comoving;
  std::vector<double> particle_accel_y_comoving;
  std::vector<double> particle_accel_z_comoving;
  std::vector<double> cell_accel_x_comoving;
  std::vector<double> cell_accel_y_comoving;
  std::vector<double> cell_accel_z_comoving;
  // v24 hierarchical synchronized endpoint; see GravityForceCachePersistentView.
  bool hierarchical_endpoint_valid = false;
  std::uint64_t hierarchical_pm_field_version = 0U;
  std::vector<double> hierarchical_tree_accel_x_comoving;
  std::vector<double> hierarchical_tree_accel_y_comoving;
  std::vector<double> hierarchical_tree_accel_z_comoving;
  std::vector<double> hierarchical_pm_accel_x_comoving;
  std::vector<double> hierarchical_pm_accel_y_comoving;
  std::vector<double> hierarchical_pm_accel_z_comoving;
  std::vector<double> hierarchical_nearest_pair_r2_comoving;
};

// Non-owning view over an owning cache (all lanes, including v24 endpoint).
[[nodiscard]] GravityForceCachePersistentView makeGravityForceCachePersistentView(
    const GravityForceCachePersistentState& cache) noexcept;

struct RestartWritePayload {
  RestartPersistentStateView persistent_state;
  const core::IntegratorState* integrator_state = nullptr;
  // Particle scheduler; particles.time_bin is a derived compatibility mirror.
  const core::HierarchicalTimeBinScheduler* scheduler = nullptr;
  // Authoritative scheduler for gas-cell rows keyed by stable gas_cell_id in the
  // workflow-level persistence contract.  Required for v19 writes whenever gas
  // cells exist.
  const core::HierarchicalTimeBinScheduler* gas_cell_scheduler = nullptr;
  // Owning compatibility input retained for direct library callers/tests.
  // Production workflow writing should prefer gravity_force_cache_view.
  const GravityForceCachePersistentState* gravity_force_cache = nullptr;
  const GravityForceCachePersistentView* gravity_force_cache_view = nullptr;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
  parallel::DistributedRestartState distributed_gravity_state;
  OutputCadencePersistentState output_cadence_state;
  StochasticPersistentState stochastic_state;
};


struct RestartDiagnosticsSummary {
  std::string restart_schema_name;
  std::uint32_t restart_schema_version = 0;
  std::string current_boundary_kind;
  std::string last_completed_boundary_kind;
  bool restart_safe = false;
  std::uint64_t step_index = 0;
  std::uint64_t scheduler_current_tick = 0;
  std::uint32_t scheduler_max_bin = 0;
  std::uint64_t scheduler_element_count = 0;
  std::uint64_t scheduler_active_count = 0;
  std::uint64_t scheduler_pending_transition_count = 0;
  std::uint64_t gas_cell_scheduler_current_tick = 0;
  std::uint32_t gas_cell_scheduler_max_bin = 0;
  std::uint64_t gas_cell_scheduler_element_count = 0;
  std::uint64_t gas_cell_scheduler_active_count = 0;
  std::uint64_t gas_cell_scheduler_pending_transition_count = 0;
  std::uint64_t pm_cadence_steps = 0;
  std::uint64_t pm_gravity_kick_opportunity = 0;
  std::uint64_t pm_field_version = 0;
  std::uint64_t pm_last_refresh_opportunity = 0;
  std::uint64_t pm_last_refresh_step_index = 0;
  bool pm_refresh_commit_pending = false;
  bool pm_long_range_field_valid = false;
  bool output_enabled = false;
  bool output_snapshot_due = false;
  bool output_checkpoint_due = false;
  std::uint64_t output_last_completed_step_index = 0;
  std::uint64_t output_next_snapshot_step_index = 0;
  double output_next_snapshot_time_code = 0.0;
  std::uint64_t stochastic_module_count = 0;
};

struct RestartReadResult {
  core::SimulationState state;
  core::IntegratorState integrator_state;
  core::TimeBinPersistentState scheduler_state;
  // v19+: state keyed by stable gas_cell_id.  For v14-v18 reads this is
  // reconstructed from the historic cell time_bin mirror for compatibility.
  core::TimeBinPersistentState gas_cell_scheduler_state;
  std::vector<std::uint64_t> gas_cell_scheduler_ids;
  GravityForceCachePersistentState gravity_force_cache;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
  parallel::DistributedRestartState distributed_gravity_state;
  OutputCadencePersistentState output_cadence_state;
  StochasticPersistentState stochastic_state;
  RestartDiagnosticsSummary diagnostics;
  // Legacy v14-v22 FNV-1a identity retained for compatibility diagnostics.
  std::uint64_t payload_hash = 0;
  std::string payload_hash_hex;
  // v23+: canonical, endian-independent integrity digest.
  std::string payload_integrity_algorithm;
  std::string payload_integrity_sha256_hex;
  // Transient read receipt, set only after completeness and payload verification.
  // Not an on-disk schema field; synthetic embedding payloads default to false.
  bool payload_integrity_verified = false;
};

// Exact admission or constrained treepm_normalized_defaults_v2 migration.
// Existing v1 compatibility receipts survive subsequent exact resumes.
// Returns a separately identified provenance record without mutating raw input.
// Migration requires a verified v23/provenance_v7 checkpoint, a closed global
// KDK boundary and the exact current normalized text minus the two new default
// TreePM lines. Insufficient/opaque historical provenance fails closed.
[[nodiscard]] std::string restartConfigurationCompatibilityRecord(
    const RestartReadResult& restart, const core::FrozenConfig& current);

[[nodiscard]] std::uint64_t restartPayloadIntegrityHash(
    const RestartWritePayload& payload, core::MemoryGovernor* governor = nullptr);
[[nodiscard]] std::string restartPayloadIntegrityHashHex(
    const RestartWritePayload& payload, core::MemoryGovernor* governor = nullptr);

void writeRestartCheckpointHdf5(
    const std::filesystem::path& output_path,
    const RestartWritePayload& payload,
    const RestartWritePolicy& policy = {});

struct RestartVerificationPolicy {
  core::MemoryGovernor* memory_governor = nullptr;
  // Application-owned verification scratch. Independent of population size.
  std::uint64_t workspace_limit_bytes = 16ULL * 1024ULL * 1024ULL;
};

struct RestartVerificationResult {
  bool ok = false;
  RestartDiagnosticsSummary diagnostics;
  std::string schema_name;
  std::uint32_t schema_version = 0;
  std::uint64_t legacy_payload_hash = 0;
  std::string payload_sha256_hex;
  std::uint64_t verified_bytes = 0;
  std::uint64_t workspace_high_water_bytes = 0;
  std::uint64_t verified_dataset_count = 0;
  std::string failed_field;
  std::string detail;
};

[[nodiscard]] RestartVerificationResult verifyRestartCheckpointHdf5(
    const std::filesystem::path& input_path,
    const RestartWritePayload& expected,
    const RestartVerificationPolicy& policy = {});

[[nodiscard]] RestartReadResult readRestartCheckpointHdf5(
    const std::filesystem::path& input_path,
    const RestartReadPolicy& policy = {});

}  // namespace cosmosim::io
