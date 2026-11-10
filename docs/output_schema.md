# Output schema guide

20261012 adds source identity, source-manifest SHA-256 and dependency-provider
text to existing build/derived-runtime provenance, without changing science
snapshot v6, checkpoint v25 or provenance v7 field schemas. Source archives
retain unavailable Git identity; a manifest checksum is explicitly a source
identity, never a fabricated commit. Normalized configuration includes five
new conservative TreePM controls; verified predecessor migration uses the
documented `treepm_normalized_defaults_v2` receipt. See
[restart admission](restart_checkpointing.md) and the
[source handoff](repair/treepm_deep_hotpath_attempt1_closure_20261012.md).

This guide documents snapshot, restart, and provenance payload contracts currently implemented in CosmoSim.

Authoritative interfaces:

- `include/cosmosim/io/snapshot_hdf5.hpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `include/cosmosim/core/provenance.hpp`
- `include/cosmosim/core/profiling.hpp` (operational run-event report)

## 0) Canonical CHUÍ initial-condition schema

The standalone converter writes header schema `chui_canonical_v1`, version 2.
It is an interchange input bundle, not a restart checkpoint and not a normal
time-series snapshot. The HDF5 member uses canonical GADGET/AREPO group and
dataset names:

- `/Header`
- `/Provenance/ConversionManifestJson`
- `/PartType0`, `/PartType1`, `/PartType4`, `/PartType5` when populated
- `Coordinates`, `Velocities`, `Masses`, `ParticleIDs`
- species fields such as `InternalEnergy`, `Density`, `StellarFormationTime`,
  `InitialMass`, `Metallicity`, `BH_Mass`, and `BH_Mdot`

The header records cosmology, box, epoch, reconstructed 64-bit counts, canonical
unit scales, comoving-coordinate and physical-peculiar-velocity conventions,
and the SHA-256 binding of the normalized audit manifest. It also records
converter evidence:

- `ChuiConverterFullStateMaterialized = 0`;
- `ChuiConverterFlow = source_chunk_to_validate_convert_append`;
- chunk capacity, peak converted batch records, and peak reader staging bytes.

The paired JSON uses `chui_ic_audit_manifest` version 4. The exact normalized
JSON is embedded in HDF5 and emitted as a sidecar. A completion marker binds the
canonical filename, sidecar filename, and digest. Runtime import recomputes and
compares both hashes and requires the marker before materialization; a merely
well-formed digest string is not accepted as verification.

HDF5, manifest, and marker are written as `.part` members and finalized as one
recoverable bundle. Any write, digest, or rename failure removes partial and
orphaned final members.

The canonical converter currently emits one HDF5 member. It fails before file
creation when any `NumPart_ThisFile` family would exceed `UINT32_MAX`, because
that header attribute has no high-word companion. Accepted total counts retain
`NumPart_Total` plus `NumPart_Total_HighWord`; the converter never silently
truncates a one-file local count.

Canonical IC, science snapshot, and restart remain deliberately separate
schemas: canonical IC is audited import/interchange state; snapshot is analysis
output; restart is exact execution continuation.

## 1) Science snapshot schema and external dialects

Current native identity:

- `schema_name = chui_science_snapshot_v6`
- `schema_version = 6`

Science snapshots are not restart checkpoints. `/PartTypeN` uses canonical phase-space
names (`Coordinates`, `Velocities`, `Masses`, `ParticleIDs`) while CHUÍ-specific science
lanes are explicitly named extensions. `/Header`, `/Config`, `/Parameters`, `/Provenance`,
and `/Units` make the file self-describing. `NamingRulesVersion` and
`FileNamingRulesVersion` are persisted.

`SnapshotDialect` separates CHUÍ-native storage from AREPO Format-3 and GADGET-4 HDF5
semantics. For comoving AREPO/GADGET-family export, the conversion boundary is explicit:
coordinates and masses carry the configured `h` convention, stored velocities are physical
peculiar velocity divided by `sqrt(a)`, and density/pressure apply the corresponding `h^2`
conversion. The inverse conversion is applied on import. Arbitrary external populated
PartType2/PartType3 groups are rejected unless a species map resolves their meaning.

`/PartType0` rows come from authoritative dense gas cells and retain stable `gas_cell_id`,
optional parent identity, and owning patch. Gas science output includes internal energy,
density, pressure, metallicity, star-formation/effective-ISM diagnostics where applicable,
and CHUÍ temperature/sound-speed diagnostics. Stars preserve formation/birth identity and
current cumulative stellar-evolution state. PartType5 preserves the current black-hole
science sidecar using documented `CHUI_` extension names for model-specific quantities.

Science snapshot topology is an I/O policy rather than an MPI-rank invariant. With `output.snapshot_layout = auto`, serial output is one `snap_###.hdf5`; MPI output also uses one collective `snap_###.hdf5` when the build has MPI-enabled Parallel HDF5. All ranks write owned hyperslabs directly into global datasets, `NumFilesPerSnapshot = 1`, and `NumPart_ThisFile` is the complete file count. Before a collective single-file publication, every rank must agree on the same normalized absolute shared run directory; disagreement fails collectively before any rank enters Parallel-HDF5 file creation. Publication is transactional: ranks write `snap_###.hdf5.partial`, perform exact bounded partition readback, agree success collectively, and only then publish the final filename. The accompanying `chui_snapshot_set_v3` `.complete` marker records `distributed_science_readback_v1` integrity without a root-only whole-file hash pass. A cadence snapshot committed by the current runtime instance may satisfy final-endpoint closure for that same step; a committed namespace left by another invocation is not treated as idempotent and remains protected by overwrite refusal.

The explicit `sharded` compatibility layout continues to use rank members plus a `chui_snapshot_set_v2` per-member SHA-256 completion marker. `aggregated`/`snapshot_num_files` reserve a future fixed-file-count backend and currently fail closed. Legacy `snapdir_###/<generation>.complete` sets remain readable. Science snapshot topology is independent from the same-world-size, rank-local restart topology.

See `docs/snapshot_hdf5_io.md` for the detailed contract.

## 2) Restart schema

Current restart identity:

- `name = cosmosim_restart_v25`
- `version = 25`

Schema v25 adds `/state/amr_temporal_boundary_history/start_internal_energy_density_comoving`
and `end_internal_energy_density_comoving` (float64, one row per history cell, group
`schema_version = 3`): the dual-energy auxiliary `rho*e` endpoints used by temporal
coarse-to-fine AMR ghosts (NaN = unknown; v24 and older restore NaN).

Per-group layout identities (`@schema_version` on the group), bound to the outer schema:

| Group | v18–v24 | v25 |
|---|---|---|
| `/state/amr_pending_flux_registers` | 2 (metal lanes optional by presence) | 2 (layout unchanged) |
| `/state/amr_temporal_boundary_history` | 1 or 2 | exactly 3 |

The schema validator, the reader and the write-time readback verifier all reject a v25
file whose temporal-history group is not schema 3 or whose pending flux-register group is
not schema 2, and the reader rejects group schema 3 inside a pre-v25 file. Files from the
unreleased patch-56 writer (pending group 3, temporal group 2) are therefore refused
rather than coerced.

Schema v24 adds `/gravity_force_cache/hierarchical_endpoint` (attributes `valid`,
`pm_field_version`; float64 datasets `tree_accel_{x,y,z}_comoving`,
`pm_accel_{x,y,z}_comoving`, `nearest_pair_r2_comoving`, rows ordered like
`/gravity_force_cache/particle_id`). It is valid only at a hierarchical KDK synchronization
endpoint whose Tree and PM split force was evaluated for every particle at the checkpoint
positions; otherwise `valid = 0` and the datasets are empty. `nearest_pair_r2_comoving` may be
`+inf` (no direct pair inside `r_cut`). The lanes are covered by both integrity digests. v23
files remain readable and import with no endpoint (the resumed hierarchical segment then
performs its synchronized bootstrap evaluation).

Schema v23 retains the v22 metal/star-birth state and all v21/v20 continuation state, and
adds a canonical little-endian SHA-256 payload integrity digest. The legacy 64-bit FNV-1a
digest remains readable/written only as a compatibility diagnostic; new v23 integrity
acceptance is bound to `sha256-canonical-le-v1`. Both digests are produced in one logical
payload traversal rather than hashing the checkpoint state twice.

Checkpoint publication uses the shared transactional-file layer. Bare filenames are valid,
temporary paths are unique same-directory siblings and can never equal the final path, HDF5
objects are closed before publication, and the implementation uses platform-specific atomic
replace plus durable file/directory synchronization when the durable policy is selected. The
writer never deletes the previous valid checkpoint before publishing its replacement.

`isRestartSchemaCompatible` retains explicit backward read support for v14 through v22.
Same-world-size rank-local restart remains the documented distributed continuation topology;
science snapshot topology is independent and may be multifile.

## H1 workflow force-cache restart note

Restart schema v20 added `/gravity_force_cache`, containing the particle and gas-cell acceleration triplets consumed by the next KDK pre-kick, plus the persisted `IntegratorState::pm_refresh_enabled` policy bit. Schema v21 retains that contract. The ReferenceWorkflow writes a valid cache only at a restart-safe completed boundary and verifies it on read. Direct low-level checkpoint callers may write an explicitly invalid/empty cache, but those artifacts do not constitute exact `ReferenceWorkflow` continuation proof.

Restart payload includes:

- `SimulationState`
- hydro geometry state under `/state/cells`, `/state/gas_cells`, and `/state/patches`
  (`/state/gas_cells/velocity_[xyz]_peculiar` is persistent cell-local hydro velocity truth)
  (`patch_id`, `level`, `first_cell`, `cell_count`, `owning_rank`)
- authoritative gas-cell identity records under `/state/gas_cell_identity`:
  `gas_cell_id`, `has_parent_particle`, `parent_particle_id`, `owning_patch_id`,
  `local_cell_row`, `@local_row_reconstruction_policy`, and
  `@identity_generation_at_write`
- `IntegratorState`
- `HierarchicalTimeBinScheduler` persistent state (`TimeBinPersistentState`)
- distributed TreePM continuation metadata (`DistributedRestartState`)
- normalized config text and normalized config hash
- provenance record
- payload integrity hash
- `/restart_diagnostics` audit metadata (schema/boundary/scheduler/PM/output/stochastic summaries; non-authoritative)

Compatibility rule is explicit through `isRestartSchemaCompatible(version)`.

Distributed restart topology is part of the executable schema contract. Current v21 behavior supports only
same-world-size, rank-local continuation: the checkpoint's normalized config hash, `/distributed_gravity`
world size, PM grid, decomposition mode, per-rank slab table, owner table, and TreePM cadence/field metadata
must match the runtime before `ReferenceWorkflow` resumes. Rank-count-changing restart and arbitrary topology
remap are not represented by this schema and must fail clearly. Per-rank HDF5 checkpoint files are serial HDF5
files with rank-qualified names; they are not a parallel-HDF5/MPIO claim.

## 2.1) Field ownership table (snapshot vs restart)

| Ownership | Fields |
|---|---|
| Shared metadata contract | normalized config text/hash, provenance payload, schema identity |
| Snapshot-only (interoperable science output) | `/Header` cosmology attrs, `/PartTypeN` particle datasets, read aliases (`Position`, `VEL`, `ID`, etc.) |
| Restart-only (exact continuation state) | compact `/restart_diagnostics` audit metadata plus full `SimulationState` hot/cold lanes, `StateMetadata`, hydro patch geometry lanes (`CellSoa::patch_index`, gas-cell mirror lanes, authoritative `/state/gas_cell_identity` records, gas-cell velocity/thermodynamic lanes, `PatchSoa` descriptors and ownership), module sidecars + schema versions, `IntegratorState`, scheduler persistent state (`bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`), distributed TreePM restart state (`decomposition_epoch`, owning-rank table, PM slab layout, cadence/long-range metadata, restart policy), payload integrity hashes. `ParticleSoa::time_bin` and `CellSoa::time_bin` are retained only as derived mirrors/diagnostics; exact continuation imports scheduler state, rejects stale mirror conflicts, and rebuilds mirrors from scheduler authority. For gas cells with a local parent, the cell mirror is validated/rebuilt against the parent gas particle's scheduler entry; parentless cells retain their cell-local mirror until a future cell scheduler authority is introduced. H1 Cartesian CFL widths are derived from persisted cell centers/config and are not serialized as scratch caches. |

Additive softening sidecar persistence:
- Snapshot `/PartTypeN/GravitySofteningComoving` (`float64`, optional; per-particle, comoving units) and `/PartTypeN/GravitySofteningOverrideMask` (`uint8`, optional) are diagnostics/interchange mirrors.
- Restart `/state/particle_sidecar/gravity_softening_comoving` (`float64`, optional; per-particle, comoving units) and `/state/particle_sidecar/has_gravity_softening_override` (`uint8`, optional) are exact-continuation lanes. The mask is the authority for per-particle override truth; a value lane without a mask is a materialized default/diagnostic mirror and does not create overrides.

Snapshot and restart intentionally remain separate contracts: snapshot is analysis/interchange oriented;
restart is execution-resume oriented.


## Stage 2 timestep-authority schema note (2026-05-11)

Historical Stage 2 scheduler-authority documentation did not change snapshot/restart/provenance schemas. H2.4 historical material referenced `cosmosim_restart_v17`; v19 introduced authoritative gas-cell scheduling, v20 added checkpoint-authoritative gravity force caches, and the active schema is now `cosmosim_restart_v25` (v24 plus the temporal-history dual-energy auxiliary endpoints; v24 was v23 plus the hierarchical KDK synchronized split-force endpoint) with ordered code-time output events. The compatibility behavior is explicit: restart payloads retain `ParticleSoa::time_bin` and `CellSoa::time_bin` as mirrors for corruption detection, reject stale mirror conflicts against scheduler truth, and rebuild valid parent-backed mirrors from scheduler state on import. Gas-cell parent lineage is optional metadata; parentless cells keep cell-local hydro velocity and timestep mirror lanes without particle velocity access.

## H1 Hydro Restart Geometry Note

Schema v15 persists the restart-authoritative hydro geometry inputs needed to rebuild
Cartesian gas-cell geometry deterministically: cell centers, gas-cell mirror lanes,
authoritative `/state/gas_cell_identity` records, cell `patch_index`, and `PatchSoa`
descriptors including patch ownership. Restart validation rejects patch ranges that
overlap, omit cells, exceed the cell arrays, disagree with per-cell `patch_index`, or
disagree with identity-record `owning_patch_id`. Transient hydro reconstruction scratch,
ghost-fill buffers, face flux arrays, and derived CFL width caches remain outside restart
payloads.

## 3) Provenance payload

`ProvenanceRecord` persists:

- schema tag (`provenance_v7`)
- config schema identity (`config_schema_name`, `config_schema_version`)
- audit payloads (`raw_input_config`, `normalized_config`, `derived_runtime_state`)
- deterministic compatibility config hash (`normalized_config_hash_hex`)
- strong normalized-config integrity (`integrity_digest_algorithm=sha256`, `normalized_config_sha256_hex`)
- build identity (`git_sha`, compiler id/version, compiler flags, build preset, feature flags)
- hardware/runtime descriptors (CPU model, logical/physical core counts when known, system RAM, host, GPU/CUDA summary, MPI summary/world size/node-local rank, deterministic mode)
- UTC timestamp and compatibility hardware summary
- rank attribution (`author_rank`)
- compatibility: valid `provenance_v6` text records remain readable; current writers emit v7 and never reinterpret v6 as v7
- gravity/TreePM reproducibility contract:
  - controls: `gravity_treepm_pm_grid`, `gravity_treepm_assignment_scheme`,
    `gravity_treepm_window_deconvolution`, `gravity_treepm_asmth_cells`,
    `gravity_treepm_rcut_cells`, `gravity_treepm_update_cadence_steps`,
    `gravity_treepm_tree_opening_criterion`, `gravity_treepm_tree_opening_theta`,
    `gravity_treepm_tree_relative_force_tolerance`,
    `gravity_treepm_tree_relative_force_acceleration_floor`,
    `gravity_treepm_pm_decomposition_mode`, `gravity_treepm_tree_exchange_batch_bytes`, `gravity_treepm_pm_exchange_batch_bytes`
  - derived scales: `gravity_treepm_mesh_spacing_mpc_comoving` (`Δmesh`),
    `gravity_treepm_split_scale_mpc_comoving` (`r_s`),
    `gravity_treepm_cutoff_radius_mpc_comoving` (`r_cut`)
  - softening/backend: `gravity_softening_policy`, `gravity_softening_kernel`,
    `gravity_softening_epsilon_kpc_comoving`, `gravity_pm_fft_backend`,
    `gravity_pm_backend_capability`, `gravity_acceptance_profile_id`. The
    acceptance profile remains `unverified_current_source` until external
    source-bound release evidence is generated; the executable does not
    self-certify.
  - restart/debug continuation metadata:
    - `gravity_treepm_decomposition_epoch`
    - `gravity_treepm_restart_world_size`
    - `gravity_treepm_restart_pm_grid`
    - `gravity_treepm_restart_slab_signature`
    - `gravity_treepm_restart_kick_opportunity`
    - `gravity_treepm_restart_field_version`
    - `gravity_treepm_long_range_restart_policy`
  - zoom-gravity metadata:
    - `zoom_long_range_strategy`
    - `zoom_region_center_{x,y,z}_mpc_comoving`
    - `zoom_region_radius_mpc_comoving`
    - `zoom_focused_pm_grid`
    - `zoom_contamination_radius_mpc_comoving`

## 4) Naming and stability conventions

- Output stems are restricted to `[A-Za-z0-9_-]` for deterministic naming.
- Normalized config snapshots are written with run outputs.
- Schema names/versions are part of compatibility behavior and must be updated intentionally.
- Reference workflow operational report (`reference_operational_events.json`) is schema-tagged (`schema_version = 1`) and includes:
  - run label,
  - provenance config hash linkage (`provenance_config_hash_hex`),
  - severity-count summary and status,
  - structured events (`event_kind`, `severity`, `subsystem`, optional step/time/scale context, message, payload map).
  - gravity runtime events include:
    - `gravity.treepm_setup` (one setup event with PMGRID, assignment, deconvolution,
      ASMTH/RCUT controls, derived `Δmesh`/`r_s`/`r_cut`, cadence, softening policy/kernel,
      FFT backend),
    - `gravity.pm_long_range_field` (refresh/reuse event per gravity kick
      opportunity carrying the same PM contract plus
      `gravity_kick_opportunity`, `field_version`, `field_built_step_index`, and
      `field_built_scale_factor`). Because the event is emitted before cadence
      commit, these four values come from the pending integrator-issued
      `PmRefreshDirective`/solver decision, not the previous committed
      `PmSynchronizationState`. This corrects stale initial-refresh payloads
      without changing the event schema or cadence ownership.
    - `gravity.zoom_force_diagnostics` (per-kick zoom decomposition norms and low-resolution contamination counters).
    - `gravity.health_check` (targeted warning/fatal event for explicit gravity-state violations;
      fatal events are not downgraded into generic diagnostics).
    - `gravity.health_summary` (per-kick gravity health counters summarizing cheap always-on checks
      and heavy reference-only checks when policy enabled).

Gravity operational-event payload values are strings. Diagnostic doubles that
carry solver evidence use scientific notation with
`std::numeric_limits<double>::max_digits10` precision, including `G_code`, the
relative-MAC acceleration floor, split/cutoff/softening scales, force L2 norms,
and field/force-evaluation scale factors. This prevents small nonzero values
from being serialized as `0.000000`. It is a value-format correction within
operational report schema version 1; it does not alter snapshot/restart schema
or numerical state.

## 4.1) Diagnostics bundle maturity metadata

Diagnostics bundles (`<diagnostics_stem>_<diagnostic_class>_step_*.json`) now include explicit maturity metadata:

- `diagnostic_class`: cadence bucket (`run_health`, `science_light`, `science_heavy`)
- `diagnostics_execution_policy`: active policy (`run_health_only`, `run_health_and_light_science`, `all_including_provisional`)
- `diagnostic_records[]`: per-output metadata
  - `name`
  - `tier` (`infrastructure_health`, `validated_science`, `reference_science`)
  - `maturity` (`production`, `validated`, `provisional`)
  - `scalability` (`cheap`, `moderate`, `heavy_reference`)
  - `executed`
  - `policy_note`

Current intended classification:

- Production infrastructure health: `run_health_counters`
- Production infrastructure health: `gravity_health_summary`
- Validated lightweight science: `star_formation_history`, `angular_momentum_budget`, `gas_xy_slice_density`, `gas_xy_projection_density`
- Provisional heavy reference-only: `power_spectrum` (disabled unless `diagnostics_execution_policy = all_including_provisional`)

Run-health payload now includes additive gravity integrity counters:

- `gravity_softening_sidecar_size_ok`
- `non_finite_gravity_softening`
- `non_positive_particle_mass`

## 5) Change procedure for schema-affecting work

When changing snapshot/restart/provenance fields:

1. Update the corresponding interface headers and implementation.
2. Update `docs/configuration.md` if config keys or normalized text semantics changed.
3. Update validation expectations in `docs/validation_plan.md`.
4. Add/update tests in `tests/unit` + `tests/integration` + `tests/validation` as applicable.
5. Record rationale in `docs/architecture/decision_log.md`.

## Compatibility notes (2026-04-20)

- Historical snapshot schema v4 was introduced as `gadget_arepo_v4`; v5 added optional per-particle softening sidecars.
  The active science-output schema is now `chui_science_snapshot_v6` (`schema_version = 6`), which separates CHUÍ-native identity from explicit AREPO/GADGET export semantics and supports logical multifile snapshot sets.
- No external `/PartType*` dataset names were changed.
- Restart schema version/name are now `cosmosim_restart_v25`, version `25` (v25 adds the temporal-history auxiliary `rho*e` endpoints; v24 adds `/gravity_force_cache/hierarchical_endpoint`; see below). It retains the v20 force cache and adds restart-authoritative code-time output cadence with an explicit v20 compatibility default of disabled.
- Restart contract enforcement was tightened: missing continuation-critical metadata, a missing or wrong root file kind, or missing output-cadence state now fails fast with explicit path-aware errors instead of producing weak checkpoints.
- Restart schema is `cosmosim_restart_v25`; distributed TreePM state, the restart-authoritative gravity force cache, and ordered output-event state are persisted under restart-only data and covered by restart integrity hashing.
- The reader accepts the documented legacy `cosmosim_restart_v14` particle-bound import path
  by materializing `/state/gas_cell_identity` from
  `/state/gas_cells/{gas_cell_id,parent_particle_id}` with
  `has_parent_particle=true`; it requires `gas_cell_id == parent_particle_id != 0` and
  does not reinterpret `parent_particle_id=0` as absent in old files.
- New v15 files require `/state/gas_cell_identity/has_parent_particle`; malformed identity
  records with duplicate or zero `gas_cell_id`, sparse `local_cell_row`, invalid parent
  flags, or invalid patch ownership are rejected.
- Restart v6 compatibility behavior is tightened without changing payload fields: non-empty
  `CellSoa::time_bin` mirrors must match the scheduler `bin_index` of each gas cell's
  parent particle, and restart import rebuilds cells from that parent-particle scheduler
  mapping. Mixed states where `cells.size() != particles.size()` are no longer an
  exact-size validation bypass.
- Diagnostics maturity metadata is additive to diagnostics JSON bundles and does not alter snapshot/restart/provenance schema compatibility.
- Provenance v6 adds the TreePM tree opening criterion, Barnes-Hut threshold, relative-force
  tolerance, and relative-force acceleration floor. Readers accept v5 records without these
  fields and apply the historical production defaults (`com_distance`, `0.7`, `0.005`, and
  `1e-30` code-acceleration units, respectively). The zoom-gravity metadata and
  contamination-radius contract keys remain additive fields.


### AMR patch geometry restart lanes (v16)

PatchSoa now persists restart-authoritative AMR patch geometry lanes: `parent_patch_id`, `morton_key`, `origin_x/y/z_comoving`, `extent_x/y/z_comoving`, and `cell_dim_x/y/z`. Production AMR hydro requires these lanes to be explicit; legacy restart inputs without them are accepted only as non-AMR/legacy patch states and do not silently enter the production AMR hydro path.


## AMR temporal restart state (v18)

`cosmosim_restart_v19` introduced `/state/amr_temporal_boundary_history` for active local AMR coarse temporal intervals. The current `cosmosim_restart_v25` retains it (group schema 3 adds the dual-energy auxiliary endpoints) and the v20 `/gravity_force_cache`, then adds ordered code-time output events. None are part of analysis snapshots. v17 pending flux-register restart state remains supported as a legacy read path.


## Initial-condition import report counters

`IcImportReport` is runtime/profiling metadata rather than a new HDF5 snapshot
schema. Its Campaign B distributed counters distinguish:

- logical consensus phases: `logical_consensus_phase_count` and
  `routing_logical_consensus_phase_count`;
- compatibility aliases: `collective_phase_count` and
  `routing_collective_phase_count` (same logical meaning);
- actual MPI calls: `mpi_collective_call_count`, routing/non-routing totals,
  and per-operation counters for `Allreduce`, `Bcast`, `Gather`, `Gatherv`,
  `Alltoall`, and `Alltoallv`;
- normalized cost: `collectives_per_million_records`;
- stable-reader, file/dataset-open, full-hash, source-identity, batch,
  exchange, byte, staging, imbalance, `distributed_id_audit_round_count`, and
  final owner-local state counters.

For a successful routing-protocol-v3 batch,
`routing_mpi_collective_call_count == 19 * routing_batch_count` on every MPI
rank. The successful non-routing protocol is also exact:

```text
nonrouting_mpi_collective_call_count
  = 37
  + (validate_runtime_cosmology ? 1 : 0)
  + source_file_count
  + 10 * distributed_id_audit_round_count
  + mpi_bcast_call_count
```

The per-operation counters sum exactly to the total, and all ranks must report
the same actual-call count. These report fields do not alter canonical
GADGET/AREPO group, dataset, or header names.

## Star formation and gas metals

Gas snapshots write canonical `PartType0/Metallicity`, derived from conserved `metal_mass_code / gas_mass_code`, plus `GasCellIDs` for stable mesh identity. On read, gas metallicity is converted back into conserved metal mass.

Star snapshots write canonical `PartType4/Coordinates`, `Velocities`, `Masses`, `ParticleIDs`, and `Metallicity`, plus:

- `StellarFormationTime`;
- `BirthMass`;
- `StarFormationBirthKey`;
- `ParentGasCellID`;
- `BirthIntegrationTick`;
- `BirthOrdinal`.

The analysis root writes `sfr_history.csv` atomically through a `.part` file. Columns are `scale_factor`, `redshift`, `formed_mass_code`, and `cumulative_stellar_birth_mass_code`. Values are reconstructed from actual star sidecar birth masses and formation scale factors, not from an analytic rate estimate.


## Effective-ISM gas diagnostics

Gas snapshots retain canonical GADGET-style fields and add derived diagnostics:

```text
/PartType0/StarFormationRate
/PartType0/ColdCloudMassFraction
/PartType0/EffectivePressure
/PartType0/EffectiveInternalEnergy
/PartType0/IsOnEffectiveEos
```

These are analysis fields, not additional conserved authorities. Effective temperature inferred from `EffectiveInternalEnergy` is a pseudo-temperature representing unresolved pressure support, not the temperature of either subgrid phase. Gas `Metallicity` remains derived from conserved metal mass.

## Metals metadata and enrichment restart lanes

The canonical scalar dataset name remains `Metallicity` for gas and stars.
`GasCellSidecar::metal_mass_code` is the gas authority; the dataset is derived at
write time. The HDF5 `/Header` records:

- `MetalSpeciesMode`;
- `MetalDiffusionModel`;
- `MetalDiffusionTimeIntegrator`;
- `MetalDiffusionCoefficient`;
- `StellarEvolutionTablePath`.

The normalized configuration and its hash remain the complete model-identity
record. Restart state additionally stores cumulative newly synthesized metals,
unresolved enrichment mass/metals/energy/momentum, and cumulative deposited
mass/metals/energy. Missing new lanes in a legacy restart default to zero under
the existing backward-compatibility policy.

## M2D-1 optional diagnostic provenance

The existing provenance derived-runtime-state text may include `optional_diagnostic_cadence_policy=coalesced_nonpersistent` and light/heavy pending, first/latest due, and missed-count summaries. The operational event stream records actual execution epochs, coalescing, and terminal drops. No snapshot or restart schema field was added. A science product generated during catch-up represents its actual current physical epoch, never a reconstructed historical state.


## M2D I/O residency follow-up (2026-09-08)

The snapshot writer's sidecar-row lookup is now a compact, phase-local sorted
index with an optional owner-held memory-governor admission. The production
workflow supplies the existing governor; standalone callers retain the
compatible null-authority behavior. Missing/duplicate sidecar indices remain
errors. Restart module-sidecar payloads are written directly from their
authoritative byte storage without changing the dataset type, shape, contents
or schema. No snapshot/restart schema version or normalized configuration key
changes are introduced. Other writer metadata and readback capacity growth
remain subject to their existing owner contracts and external-runtime reserve.

## P6 operational report version 2

Science snapshot schema and Parallel-HDF5 single-file publication are unchanged.
Optional hierarchical snapshots occur only after all coarse intervals/kicks
close and ownership migration commits; no mixed-epoch snapshot is published.
Existing step/code-time scheduling is preserved; hierarchical step cadence counts
coarse blocks. Explicit redshift/scale-factor target lists are not added.

`operational_report.json` advances `schema_version` from 1 to 2. Its summary
`event_count`, warning/error/fatal counts are lifetime totals; recent event detail
is limited to `recent_event_limit` (4096 since the October 2026 TreePM campaign;
previously 256, which evicted early hierarchical blocks from ten-block runs).
The array layout and schema version are unchanged. Added summary fields are `retained_event_count`,
`evicted_event_count`, and `recent_event_limit`. Consumers must not infer lifetime
counts from the event array. Existing event records remain readable; ordinary
profiler JSON stays version 1. These source changes require consumer/fixture
qualification; no output or restart tests were executed in P6.

## Effective-ISM sidecar ledger units (CHUI Part 1, audit F7)

The `effective_multiphase_ism` module sidecar now records integrated code
energies (`energy_ledger_version=2`, `energy_ledger_units=integrated_code_energy`).
No HDF5 group or dataset name changes; the snapshot/restart schema versions are
unchanged. Legacy (version-less) ledgers are discarded on restart, see
`docs/restart_checkpointing.md`.

## CHUI Part 3 science-state export (composition and stellar populations)

Source only; NOT RUN — USER-RESERVED. Additive optional datasets; the
snapshot `schema_version` and every canonical GADGET/AREPO dataset are
unchanged, so existing readers are unaffected. Presence is decided from the
configuration (identical on every rank, required for collective single-file
HDF5); an absent dataset means the physics is disabled, never zero.

`/PartType0` (when the run carries a composition sidecar; attribute
`CHUI_CompositionSchemaVersion = 3`, `CHUI_CompositionLayoutName`,
`CHUI_ElementNames`, `CHUI_SpeciesNames`, `CHUI_DustNames`,
`CHUI_CompositionFractionBasis`, `CHUI_DetailedChemistryStorage`):

| Dataset | Shape | Meaning |
|---|---|---|
| `CHUI_ElementMassFractions` | [N, 11] | gas-phase element lane mass / (elements + dust) (H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe) |
| `CHUI_SpeciesMassFractions` | [N, S] | species lane mass / same basis (hybrid/detailed networks only) |
| `CHUI_DustMassFractions` | [N, D] | dust reservoir mass / same basis |
| `CHUI_IonStageMassFractions` | [N, I] | C02 metal-ion ladders only (`physics.chemistry_metal_ion_elements != none`): gas-phase nuclei mass of each ionization stage / same basis; stage names in `CHUI_IonStageNames` (`ion_C_I` .. `ion_Fe_XXVII`, selected elements only). Distribution carrier: the authoritative element total is `CHUI_ElementMassFractions` (minus CO-bound C/O) |
| `CHUI_SgsSpecificEnergy` | [N] | SGS energy / gas mass (code units), when the SGS lane exists |
| `CHUI_MeanMolecularWeight` | [N] | cached mu (amu); NaN = unknown (rebuilt by chemistry) |
| `CHUI_ChemistryModelTag`, `CHUI_IsmRegimeTag` | [N] u8 | adaptive chemistry fidelity / ISM regime tags |
| `CHUI_DetailedChemistryStateMaterialized` | [N] u8 | 1 when the adaptive C+/CO history has an explicit stable-ID sparse row; 0 means the documented implicit C+=element-C, CO=0 projection (or dense-reference storage selected by the group attribute) |
| `CHUI_ChemistryModelDwellYears`, `CHUI_IsmRegimeDwellYears` | [N] | hysteresis dwell |
| `CHUI_ChemistryLastTemperatureK`, `CHUI_ChemistryLastCoolingTimeSeconds` | [N] | last accepted local thermal/timescale history; 0 means unavailable/no finite cooling time |
| `CHUI_ChemistryLastErrorEstimate`, `CHUI_ChemistryLastElementFractionResidual`, `CHUI_ChemistryLastChargeResidual` | [N] | last accepted normalized integration error and dimensionless conservation residuals used by C03 selection |
| `CHUI_ChemistryCumulativeRadiativeLossErg`, `CHUI_ChemistryCumulativeHeatingGainErg` | [N] | cumulative accepted per-cell radiative/heating ledger; extensive under AMR regrid |

Per cell, sum(`CHUI_ElementMassFractions`) + sum(`CHUI_DustMassFractions`) = 1;
`Metallicity` remains the derived total-metal mirror (gas + dust metals).

`/PartType4` (Part 2 stellar population models; attribute
`CHUI_StellarPopulationSchemaVersion = 2`, `CHUI_StellarChannelOrder = AGB,CC,Ia`):
`CHUI_StellarReturnedMassByChannel` [N, 3], `CHUI_StellarFeedbackEnergyByChannelErg`
[N, 3], `CHUI_BirthElementMassFractions` [N, 11], `CHUI_BirthTimeCode`,
`CHUI_SnCcRealized`, `CHUI_SnIaRealized` (u32), `CHUI_EventStreamKey` (u64),
`CHUI_PopulationRecordInitialized` (u8; 0 = legacy star without a record),
`CHUI_ChronologyTimeCode`, `CHUI_ChronologyStepIndex` (u64), and
`CHUI_ChronologyInitialized` (u8). The chronology triplet identifies the last
accepted SSP/event interval represented by the cumulative ordinals and is
analysis provenance; checkpoint restart still uses the complete module
sidecar rather than reconstructing it from a science snapshot.

Restart-only transient state (flux registers, temporal histories, Jacobians,
staged deltas) is never exported to science snapshots; it lives in the
checkpoint (`docs/restart_checkpointing.md`).
