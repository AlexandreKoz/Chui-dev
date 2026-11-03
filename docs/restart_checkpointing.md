# Restart and Checkpoint System

## Scope and schema

CosmoSim restart checkpoints are **exact-continuation artifacts** and intentionally richer than analysis snapshots.
The restart schema (`cosmosim_restart_v24`) persists:

- full `SimulationState` hot/cold SoA lanes (through a narrow `RestartPersistentStateView`),
- `StateMetadata` blob,
- hydro geometry ownership lanes: cell centers, gas-cell identity, cell-local velocity, and thermodynamic fields,
  cell `patch_index`, and `PatchSoa` descriptors (`patch_id`, `level`, `first_cell`,
  `cell_count`, `owning_rank`),
- authoritative gas-cell identity records under `/state/gas_cell_identity`, including
  `gas_cell_id`, `has_parent_particle`, `parent_particle_id`, `owning_patch_id`,
  `local_cell_row`, `local_row_reconstruction_policy=explicit_dense_local_cell_row`,
  and `identity_generation_at_write` as an audit stamp,
- module sidecars (`ModuleSidecarRegistry`) with per-module schema versions,
- `IntegratorState`, including the `pm_refresh_enabled` policy bit,
- restart-authoritative `/gravity_force_cache` acceleration triplets used by the next KDK pre-kick,
- v24 `/gravity_force_cache/hierarchical_endpoint`: at a closed hierarchical KDK block the
  Tree and PM split components and the nearest direct-pair separation squared for every
  particle, so the next block opens with the same force in uninterrupted and resumed runs,
- hierarchical scheduler persistent state (`TimeBinPersistentState`),
- distributed TreePM continuation state (`distributed_gravity_state`, schema-versioned),
- normalized config text/hash and provenance payload,
- output cadence persistence state for deterministic future snapshot/restart naming and cadence decisions,
- root file-kind metadata (`cosmosim_file_kind=restart_checkpoint`),
- compact `/restart_diagnostics` audit metadata,
- canonical little-endian SHA-256 payload integrity (`sha256-canonical-le-v1`) for v23,
  plus the legacy FNV-1a 64-bit digest only for backward-compatibility diagnostics.

By design, this differs from GADGET/AREPO-style analysis snapshots where scheduler internals and opaque sidecars are not mandatory.

## Transactional publication (v23)

Checkpoint publication uses the shared transactional-file layer. Bare filenames are valid.
Temporary siblings are unique and same-directory, an empty/path-like suffix is rejected before
file creation, HDF5 handles are closed before publication, and the previous valid checkpoint is
never deliberately removed before the replacement becomes visible. Durable publication fsyncs
the temporary file and containing directory on POSIX; the guarded Windows path uses
`FlushFileBuffers` plus `ReplaceFileW`/`MoveFileExW(...WRITE_THROUGH)`. Windows runtime
behavior remains a platform acceptance item rather than an unexecuted claim.

## Ordered output timeline events (v21)

Schema v21 adds `snapshot_interval_time_code` and `next_snapshot_time_code` to
`/output_cadence`, with the latter also mirrored in `/restart_diagnostics` for audit.
The interval is anchored at the configured begin time. `ReferenceWorkflow` clips a
larger KDK interval before it crosses the next event, writes at the resulting safe
boundary, advances the persisted event, and stores the pre-clip interval as the next
step proposal. This keeps uninterrupted and restart continuation schedules aligned.
Step-modulo cadence remains supported and may be enabled alongside code-time cadence.
The integrity hash covers both new fields. Legacy v20 reads set both fields to zero,
which explicitly disables code-time events rather than inventing a future event.

## H1 fixed/patch workflow restart note (v20)

The v20 addition is deliberately narrow: it persists the committed particle and gas-cell gravity acceleration cache at a restart-safe KDK boundary. `ReferenceWorkflow` restores that cache before its next pre-kick, while continuing to rebuild PM mesh work arrays at their next legal refresh surface. This closes the workflow-level restart gap without serializing transient hydro scratch, ghost descriptors, PM meshes, tree traversal workspaces, or MPI buffers.

Restart write payloads now carry `RestartPersistentStateView` (`persistent_state.simulation_state`) rather than a broad ad-hoc object pointer. This creates an explicit outer-boundary guard: transient runtime owners such as `TransientStepWorkspace`, `HydroScratchBuffers`, PM work arrays, TreePM traversal scratch, MPI send/recv buffers, and output staging buffers are out-of-scope for restart traversal by type.

## Gravity force history, PM cadence, and relative-MAC continuation

The `/gravity_force_cache` is the committed dense particle/gas-cell acceleration
state consumed by the next KDK pre-kick. Stable `particle_id` and `gas_cell_id`
lanes make it possible to remap this cache after a row reorder without treating
dense row as physical identity. On import, all components must be finite and
the stable-ID lanes must be unique and exactly cover the current state.

The same committed particle acceleration is the only production source for the
relative-force tree MAC's previous-acceleration magnitude. The workflow uses it
only when the cache is valid, its dense extents match, and its particle-index
generation matches the current `SimulationState`. Otherwise the traversal uses
the deterministic COM-distance fallback. A restart therefore does not invent a
force scale, and a migration/compaction generation change cannot silently feed
stale row-indexed history to the MAC.

The safe production default remains `numerics.hierarchical_max_rung = 0`,
with `numerics.treepm_update_cadence_steps = 1`. The provisional collisionless
hierarchical path is explicitly opt-in and is not scientifically certified.
Rung-zero resumes require both schedulers to have `max_bin=0`, all committed bins
zero, and pending bins zero or unset. Opt-in hierarchical resumes require the
configured rung, a closed synchronization boundary, matching source drift epochs,
and no open active interval. Rank-coordinated coarse endpoints must report an
actually refreshed PM field at the exact endpoint scale factor and step index.
Ranks validate these contracts before production stepping. The
integrator/distributed restart state remains authoritative for:

- gravity-kick opportunity and refresh decision;
- long-range field version and last refresh opportunity;
- PM field build step and scale factor;
- same-topology PM slab ownership.

PM real/spectral meshes, periodic unwrapped tree coordinates, tree topology,
multipoles, hierarchy packets, and MPI buffers are transient. They are never
restart truth. The PM mesh follows `deterministic_rebuild` at the next legal
refresh; the force cache preserves the already committed kick input. Tree and
hierarchy state are rebuilt from owner particles. Their actual particle-
decomposition epoch is restored from distributed restart truth, while force and
exchange identity follows the resumed field/exchange sequence.

The in-memory TreePM long-range cache validity signature is transient as well.
It covers force epoch, force-evaluation scale factor, `G_code`, split scale,
rectangular box axes, assignment, boundary, PM decomposition mode, and
deconvolution. It is not serialized as a second owner. A newly constructed or
signature-incompatible coordinator therefore rejects an explicit reuse request
coherently; it does not silently reinterpret reuse as refresh. Divergent
explicit refresh/reuse votes fail before PM collectives. Particle-decomposition
epoch is deliberately excluded because PM ownership remains the fixed FFT slab
map; tree protocol records still require that epoch.

The distributed decomposition epoch is a workflow-owned generation, not a
rank-local dense-row counter. It advances once on every rank only when a
measured rebalance commits an actual particle-ownership change, is persisted in
`/distributed_gravity/state` and provenance, and is restored before TreePM use.
The dense-row acceleration cache is invalidated immediately on such a commit,
so a checkpoint written after migration honestly stores `cache.valid=false`
rather than stale row-indexed force history.

`provenance_v6` introduced the normalized tree-opening criterion, theta,
relative-force tolerance, and relative-acceleration floor. Current writers emit
`provenance_v7`, which retains those fields and adds a SHA-256 normalized-config
integrity digest plus explicit compiler/hardware/runtime descriptors. Readers
remain backward-compatible with valid v6 records; v6 bytes are never reinterpreted
as v7. The provenance schema does not create a second live runtime authority or
change restart ownership semantics.


## Restart-safe boundary contract

Restart checkpoints may only be written from a completed, globally coherent restart boundary. The runtime predicate `core::evaluateRestartBoundary(...)` is the narrow contract used by workflow output dispatch and HDF5 restart payload validation. It rejects half-step KDK states, local active-bin substeps, non-restart-safe boundary kinds, and PM refresh transitions with an uncommitted long-range refresh event.

Failed restart requests are hard errors, not warnings or silent skips. The diagnostic records the current and last completed boundary kind, `inside_kdk_step`, `last_completed_restart_safe`, local-substep activity, PM refresh legality/commit-pending state, `step_index`, and scheduler tick when available. Intentionally represented half-step or local-substep restart is not implemented in schema v23, so those states must not be serialized as persistent truth.


## Restart diagnostics metadata

Schema v15 retains the compact `/restart_diagnostics` group. It records the schema identity, current and last completed boundary kind, restart-safe decision, scheduler tick/bin/active/pending counts, PM cadence/field-version summary, output cadence summary, and stochastic module count. These fields are audit metadata only: authoritative continuation remains the serialized `SimulationState`, `IntegratorState`, scheduler arrays, output cadence state, stochastic state, distributed TreePM state, normalized config, provenance, and integrity hash.

## File format and compatibility

- Format: HDF5 (`writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5`).
- Root file-kind gate: restart readers require `cosmosim_file_kind=restart_checkpoint` and reject ordinary `science_snapshot` files before reading runtime truth.
- Schema version gate: `isRestartSchemaCompatible(file_schema_version)`.
- Current compatibility policy: write current v24; read v23 (no hierarchical endpoint: the resumed hierarchical segment bootstraps its first block), v22/v21 plus documented v20/v19/v18/v17/v16/v15/v14 paths. v20 has no code-time output-event fields and materializes that lane as disabled. v19 also has no serialized gravity force cache and therefore follows the explicit safe-bootstrap compatibility path rather than claiming bitwise workflow continuation.
- v14 compatibility materializes `/state/gas_cell_identity` from
  `/state/gas_cells/{gas_cell_id,parent_particle_id}` with `has_parent_particle=true`
  and requires `gas_cell_id == parent_particle_id != 0` for every cell. It does not
  reinterpret `parent_particle_id=0` as absent in old files.
- `v5` and older restart files are intentionally rejected because they can omit exact sidecar truth lanes required by v6, including materialized softening override mask/value datasets.

## Atomic write semantics

Writers always emit to `<final>.part` first and only then rename into final path.
Finalization uses a direct rename of the temp artifact (no pre-remove step), so the old
restart path is never explicitly deleted before replacement. This preserves atomic replace
behavior on filesystems where `rename` is atomic.

## Integrity and provenance

- `restartPayloadIntegrityHash` hashes state+integrator+scheduler+config text/hash+provenance and first validates that derived particle `time_bin` mirrors match scheduler `bin_index` authority. For gas cells with a local parent particle, derived cell `time_bin` mirrors are validated through the parent scheduler entry; parentless gas-cell time-bin mirrors remain cell-local diagnostics until a future cell scheduler authority is introduced. Restart validation also requires hydro patch descriptors to cover every gas-cell row exactly once and requires each cell's `patch_index` to agree with the serialized patch range. H1 Cartesian CFL metadata remains derived from persistent cell centers/config and is not serialized as a scratch cache.
- The hash covers particle lanes, sidecars, gas-cell mirror lanes, authoritative gas-cell identity records (`gas_cell_id`, `has_parent_particle`, `parent_particle_id`, `owning_patch_id`, `local_cell_row`, and reconstruction policy), gas-cell velocity and thermodynamic lanes, species counts, full star/BH/tracer sidecars (including stellar-evolution cumulative lanes), integrator
  time-bin context, scheduler persistent arrays (`bin_index`, `next_activation_tick`,
  `active_flag`, `pending_bin_index`), output cadence persistence fields, and the full serialized provenance payload. State `time_bin` arrays remain hash inputs for corruption detection, but restart continuation imports scheduler state and rebuilds mirrors rather than treating mirrors as fallback authority.
- Because gravity TreePM metadata now lives in `ProvenanceRecord`, restart artifacts carry
  and integrity-protect the exact PM controls (`pm_grid`, assignment, deconvolution,
  `asmth_cells`, `rcut_cells`, cadence), derived scales (`Δmesh`, `r_s`, `r_cut`),
  softening policy/kernel/epsilon, and FFT backend name.
- `G_code` is deterministically re-derived from the integrity-protected frozen
  `units.{length_unit,mass_unit,velocity_unit}` through
  `core::newtonGravitationalConstantCode(UnitSystem)`. It is not a second
  persisted authority and requires no restart schema field. PM and short-range
  TreePM both consume that same scale-free value. Scale factor and Hubble rate
  remain integrator truth; collisionless KDK and the gas conservative source
  apply the cosmological response to their respective state.
- Distributed continuation metadata is persisted under `/distributed_gravity/state` (text-encoded
  `parallel::DistributedRestartState`) and includes:
  - decomposition epoch,
  - per-particle owning rank mapping (`owning_rank_by_item`),
  - PM grid shape and slab ownership table (`pm_slab_begin_x_by_rank`, `pm_slab_end_x_by_rank`),
  - kick-opportunity cadence state (`gravity_kick_opportunity`, `pm_update_cadence_steps`),
  - long-range field refresh metadata (`long_range_field_version`,
    `last_long_range_refresh_opportunity`, `long_range_field_built_step_index`,
    `long_range_field_built_scale_factor`),
  - restart long-range policy (`long_range_restart_policy`).
  Distributed floating-point restart metadata is text-encoded with enough precision for exact read-back comparisons.
- **Supported distributed restart topology:** same-topology, rank-local restart only. A continuation must use the
  same MPI world size, same runtime rank for each rank-qualified checkpoint file, same normalized config hash,
  same PM grid shape, same PM decomposition mode, same per-rank PM slab table, and coherent TreePM cadence/field
  metadata before mutable runtime state is used. `ReferenceWorkflow` validates these fields on resume and rejects
  mismatches with expected/observed rank, world-size, slab, config, or cadence diagnostics. Rank-qualified filenames
  are a file-mapping guard, not a remapping algorithm.
- **Unsupported topology changes:** rank-count-changing restart, arbitrary particle/cell remap, and remote-row
  assignment from another rank's checkpoint are intentionally rejected. Each rank is authoritative only for its
  local restart payload under the validated topology. Imported PM/tree/hydro ghosts are transient and are rebuilt
  from owner state after restart; they are not persisted as truth.
- **Policy:** restart continuation uses `long_range_restart_policy=deterministic_rebuild` for PM mesh-field storage. The v20 force-cache addition, retained by v23, serializes the acceleration cache that the next KDK pre-kick consumes, so `ReferenceWorkflow` does not discard already-committed force state across a restart. PM mesh arrays themselves remain non-persistent and are rebuilt on the next legal refresh opportunity.
- Tracer restart payload includes host-coupling lanes (`host_cell_index`, `mass_fraction_of_host`,
  `last_host_mass_code`, `cumulative_exchanged_mass_code`) for deterministic continuation.
- Writer stores the canonical SHA-256 digest and also the legacy integer/hex FNV digest for migration diagnostics; both are computed in one logical traversal.
- Reader recomputes hash and rejects mismatches.
- Provenance is serialized with the checkpoint for continuation auditing.

## Exact-restart completeness checklist

The checklist exposed by `exactRestartCompletenessChecklist()` and enforced by
restart write/read checks is:

1. `simulation_state_lanes_and_metadata`
2. `particle_identity_softening_and_drift_epoch_lanes`
3. `gas_cell_identity_lanes`
4. `hydro_geometry_patch_state`
5. `species_specific_sidecars`
6. `module_sidecars_with_schema_versions`
7. `integrator_state`
8. `integrator_owned_pm_sync_state`
9. `gravity_force_cache_at_kdk_boundary`
10. `scheduler_persistent_state`
11. `gas_cell_scheduler_persistent_state_keyed_by_gas_cell_id`
12. `output_cadence_persistent_state`
13. `stochastic_module_persistent_state`
14. `restart_diagnostics_summary`
15. `distributed_gravity_state`
16. `normalized_config_text_and_hash`
17. `provenance_record`
18. `payload_integrity_hash_and_hex`
19. `amr_pending_flux_register_state`
20. `amr_temporal_boundary_history_state`

`validateContinuationMetadata(...)` now requires non-empty normalized config text, a normalized config hash that matches the text, and a provenance config hash that matches the normalized config hash before hashing/writing restart payloads.

## Exactness policy

Restart checkpoints target exact continuation on the same build/feature path. Current tests
expect bitwise equality for persisted arrays and strict scalar equality, with no tolerance
window except tiny floating-point assertion slack in selected test comparisons (`< 1e-15`)
to avoid false negatives from host-side literal conversions.

## Runtime-truth round-trip and compatibility behavior

Restart round-trip invariant coverage now explicitly checks that:

- particle IDs and canonical particle order survive reload,
- species counts/index reconstruction and gas identity-by-particle-ID survive reload,
- scheduler persistent state survives and reconstructs the same active set at resume tick,
- particle `time_bin` mirrors are rejected when stale and rebuilt from scheduler state on successful read,
- particle-bound gas-cell `time_bin` mirrors are rejected when stale against their parent gas particle scheduler entry and are rebuilt from that mapping on successful read,
- hydro patch geometry state is rejected when patch ranges overlap, omit gas cells, extend beyond
  the cell arrays, or disagree with per-cell `patch_index`,
- per-particle softening overrides survive and continue to take precedence over species/global defaults,
- normalized config text/hash and provenance config hash remain aligned after reload.

Compatibility policy for legacy restart payloads is explicit:

- Missing required continuation fields (for example scheduler persistent lanes such as
  `pending_bin_index`) are rejected with clear read errors.
- Stale particle `time_bin` mirrors that disagree with `/scheduler/bin_index` are rejected by hash/write/read validation before exact continuation is accepted. Stale particle-bound gas-cell `time_bin` mirrors are checked against the parent gas particle scheduler entry. This originated as a v6 compatibility check; the current schema is v15 after H2 gas-cell identity promotion added explicit identity records without making diagnostics authoritative truth.
- v15 malformed identity records are rejected on read: duplicate or zero `gas_cell_id`,
  missing `/state/gas_cell_identity/has_parent_particle`, non-binary parent flags,
  `parent_particle_id=0` when `has_parent_particle=true`, nonzero parent values when
  `has_parent_particle=false`, sparse or duplicate `local_cell_row`, and identity
  `owning_patch_id` values that do not match each cell's owning `PatchSoa::patch_id`.
- The v6 reader requires `/state/particle_sidecar/gravity_softening_comoving` and `/state/particle_sidecar/has_gravity_softening_override` datasets to be present. The mask remains authoritative; empty datasets encode no materialized per-particle softening lane, while populated datasets must satisfy `ParticleSidecar` ownership invariants. Missing lanes are rejected instead of being interpreted as implicit defaults.

## Parallel and scale-up note

MPI workflow restart uses rank-qualified serial-HDF5 files and supports exact
same-world-size, same-rank, same-PM-layout continuation. It is not parallel
HDF5/MPIO, and it cannot change rank count or remap one rank's checkpoint rows
onto another rank. The np1--np4 DMO gate exercises direct-vs-resumed state under
the supported topology.

## Runtime-state exactness check for reference workflow restart verification

The reference workflow now treats restart verification as exact runtime-state continuation, not a count-only smoke test. The verification path compares all restart-authoritative state families: particle hot SoA lanes, particle metadata/softening override lanes, gas-cell identity and hydro fields, patch state, star/black-hole/tracer sidecars, species count ledgers, opaque module sidecar payloads, scheduler state, provenance/config metadata, and distributed gravity restart metadata.

Module sidecar payload length is explicitly included in the restart payload integrity hash before payload bytes. This preserves hash boundary clarity across adjacent module payloads and makes omitted or reshaped module truth detectable by the hash.


### AMR patch geometry restart lanes (v16)

PatchSoa now persists restart-authoritative AMR patch geometry lanes: `parent_patch_id`, `morton_key`, `origin_x/y/z_comoving`, `extent_x/y/z_comoving`, and `cell_dim_x/y/z`. Production AMR hydro requires these lanes to be explicit; legacy restart inputs without them are accepted only as non-AMR/legacy patch states and do not silently enter the production AMR hydro path.

## AMR hydro restart-equivalence coverage

The H3 hardening pass added `integration_restart_equivalence_amr_hydro`, an HDF5-gated integration test that compares a direct local AMR hydro continuation against a checkpoint/reload/continue path. The test exercises production AMR coverage, explicit PatchSoa geometry lanes, stable `GasCellIdentityMap` records, ghost-fill/reflux safety, scheduler state, integrator state, output cadence state, and stochastic persistent state.

The restart-equivalence harness now compares explicit AMR patch lanes (`parent_patch_id`, `morton_key`, origins, extents, cell dimensions), `GasCellIdentityMap` records, and the gas identity generation value. Restart read restores the serialized gas identity generation with `GasCellIdentityMap::assignWithGeneration(...)` so a checkpoint/reload path does not masquerade as a runtime identity rebuild.

This proves local HDF5 AMR hydro restart equivalence for the exercised synchronized-sweep scenario. It does not prove MPI restart, restart after AMR migration, Berger-Colella subcycling, or replay of persistent pending flux registers.

---

The narrow TreePM normalized-defaults compatibility projection accepts only an integrity-verified v23 or v24 checkpoint with matching historical provenance, an exact legacy-default projection, and a closed global KDK state. V24 is covered by the full canonical SHA-256 digest (including output-cadence state); this does not introduce open-ended configuration migration or change the restart schema.

## AMR pending flux-register restart lanes (v17)

Restart schema v17 introduced pending AMR flux-register state; current v23 checkpoints retain it under:

```text
/state/amr_pending_flux_registers
```

The group persists stable identity and coverage metadata for deferred reflux records, including:

- `register_key`
- `coarse_patch_id`
- `coarse_gas_cell_id`
- `coarse_cell_index`
- `level`, `axis`, `orientation`
- expected/coarse/fine area coverage
- interval start/end and coarse timestep
- expected/completed fine substeps and substep coverage mask
- coarse/fine face counts
- gas-cell identity generation and patch geometry generation
- coarse/fine flux-integrated mass, momentum, and total-energy contributions

Old restart files without this group are treated as legacy inputs with an empty pending store only at a synchronized point. A pre-v18 payload that does contain pending deferred-reflux records is rejected by the temporal-history reader, because it cannot prove the matching temporal boundary history needed for safe mid-subcycle continuation. Current v17 restart schema validation required the group and datasets so a new v17 restart could not accidentally omit pending deferred reflux truth.

The persisted `patch_geometry_generation` participates in the checkpoint integrity hash. On read, CHUI verifies the original complete checkpoint digest **before** rebasing this transient dense-row epoch to the reconstructed state's `cellIndexGeneration()`. The rebase requires the stored `coarse_cell_index`, stable gas-cell identity, owning patch ID, and identity generation to agree with the restored state; any mismatch is rejected. The on-disk schema and gas-cell identity generation are not modified. This also covers restart after distributed AMR migration.

The restart payload integrity hash includes the pending-register store for v17 payloads. The HDF5-gated `integration_restart_equivalence_amr_flux_registers` test writes an incomplete pending register before restart, reads it back, merges the missing fine contribution, applies completed reflux, and compares direct vs restart continuation.

This proves the exercised single-rank pending-register restart path. It does not prove MPI AMR restart, restart after AMR patch migration, or scheduler-owned Berger-Colella subcycling.

---

## AMR temporal boundary-history restart lanes (v18)

Restart schema v19 introduced `/state/amr_temporal_boundary_history`; current `cosmosim_restart_v23` retains it for open local AMR
coarse intervals. Each history record persists patch ID/level, geometry fingerprint, gas identity generation,
interval start/end, completion state, and stable-ID patch-local conserved start/end records. The restart
payload integrity hash includes this store for v18 payloads.

`integration_restart_equivalence_amr_temporal_ghosts` checkpoints after coarse end-state capture and a
first fine temporal ghost use, with an incomplete pending reflux record present. The restarted branch reads
the history, uses it for the midpoint coarse-to-fine ghost, completes the deferred flux contribution, and
matches the direct continuation. This is single-rank local coverage only.

Legacy schemas before v18 do not contain temporal history. They may be read only for states that do not
require resumption inside an active temporal AMR interval; legacy files are not a supported way to resume
one. The v18 reader conservatively rejects a pre-v18 payload that still contains pending AMR flux-register
records, because that observable deferred-synchronization state cannot prove the required time-aligned
coarse boundary history exists.

---

## Gas-cell scheduler persistence (v19)

Schema v19 introduced a separate gas-cell time-bin state; current v23 (`cosmosim_restart_v23`) persists it in
`/gas_cell_scheduler`. Its identity key is explicitly `gas_cell_id`, not
`parent_particle_id` and not a dense local row. The group contains `gas_cell_id`,
`bin_index`, `next_activation_tick`, `active_flag`, and `pending_bin_index`, together with
`current_tick` and `max_bin` attributes.

Restart validation checks that each scheduler identity agrees with the authoritative
`GasCellIdentityMap` and that `CellSoa::time_bin` is a derived mirror of `bin_index`. The
payload integrity hash and diagnostics include this state. v14--v18 files retain a narrow
compatibility reconstruction route based on persisted cell time-bin mirrors; that route never
looks up a cell through a parent particle.

A state with gas cells but no `PatchSoa` is legal for legacy/non-AMR Cartesian continuation only
when every `cells.patch_index` is the zero sentinel and identity records use
`owning_patch_id = 0`. AMR states continue to require complete patch coverage.

## H1 fixed-patch workflow restart evidence (v20)

For the H1 fixed/Cartesian-patch workflow proof, `/gravity_force_cache` stores stable
`particle_id` and `gas_cell_id` lanes beside cached acceleration components. On import,
`ReferenceWorkflow` remaps those components to the current dense rows by stable ID rather
than assuming write-time row order. A missing, duplicate, zero, or non-covering cache
identity lane is rejected. This is required because dense row is transient storage position,
not physical identity.

**Ordinary import versus transactional restore (2026-10-27).** These are two
deliberately separate paths:

* **Restart import.** It keeps the stable-ID remap described above
  (`internal::remapForceEndpointById`, shared with the unit test that
  compares both paths). Its behavior and messages are unchanged, and it
  validates fully before writing particle rows.
* **Multirate PM-cycle rollback.** It never goes through the restart
  importer. Inside a cycle the particle layout is frozen, so the cycle-opening
  endpoint is restored row by row from a `TransactionalForceEndpoint`, with no
  ID lane, no hash map and no allocation. The restore is accepted only when
  the particle count, the particle-index generation and the identity/layout
  certificate all match, on every rank (collective validation before any
  write).

Neither path is persisted, so there is no schema change.

The HDF5 workflow-equivalence test covers a checkpoint boundary at which gas-cell rows are
intentionally reordered before resume. It proves the exercised local workflow path, not
multi-rank migration/repartition restart. Validation posture remains: CI guards;
deterministic regression and restart equivalence; multi-resolution reference convergence;
cross-code comparison; then science-readiness benchmarks.

## Restart schema 22: star formation and conserved metals

Schema 22 adds persistent gas `metal_mass_code`, AMR metal flux registers, temporal coarse-boundary metal history, star birth key/parent/tick/ordinal fields, and the associated scheduler/source metadata. The stochastic continuation key uses the restored global integration tick and stable gas-cell identity.

Schema 21 remains explicitly readable. Fields absent from schema 21 are initialized through the documented backward-load defaults rather than silently reinterpreting an old payload under the new adaptive model. A historical restart retains its normalized model selector and configuration hash; changing model is an explicit configuration/restart compatibility decision.


## Star formation and effective-ISM continuation

Restart schema v22 persists gas metal mass, stellar immutable birth identity, scheduler sidecars, integration tick, module metadata, and latest/cumulative star-formation ledgers. The effective EOS itself is reconstructed deterministically from normalized configuration and cooling policy; the restart-persistent `effective_multiphase_ism` sidecar records parameter set, EOS schema, table hash, cooling reference, and latest/cumulative energy adjustments for verification. Since CHUI Part 1 (audit F7) the adjustments are integrated code energies `sum(rho_comoving * delta_u * V_comoving)` and the sidecar carries `energy_ledger_version=2`, `energy_ledger_units=integrated_code_energy` and `legacy_density_ledger_discarded=0|1`. A restart from a sidecar without the version key discards (does not reinterpret) the legacy energy-density sums and sets `legacy_density_ledger_discarded=1`; the HDF5 restart schema version is unchanged because the sidecar payload is versioned inside the block. Derived cold fractions and pressures are recomputed rather than stored as writable authorities.

The particle-ID registry is reconstructed exactly from loaded particle IDs and then updated incrementally. A restart cannot silently select a different star-formation model or EOS table.

## M2D-1 best-effort optional diagnostics at restart

Optional light/heavy science cadence is explicitly coalesced and nonpersistent. Checkpoint `provenance.derived_runtime_state` records `optional_diagnostic_cadence_policy=coalesced_nonpersistent` and per-class pending/first/latest/missed metadata. This is existing provenance text, not a new restart field or schema version. Resumption starts a new optional cadence and emits a policy event containing the prior summary. Pending historical products are not replayed from later physical states. Required integrator/scheduler/output cadence truth and exact restart continuation remain unchanged. Persistent pressure or a segment boundary may cause best-effort optional products to be dropped with an explicit event.

## M2D distributed metadata streaming continuation

`DistributedRestartState::serializeTo(std::ostream&)` emits the same canonical text as the retained `serialize()` API. `serializedSizeBytes()` counts the exact emitted bytes without retaining the text. The HDF5 writer uses a fixed 64 KiB buffer and hyperslab writes to the existing contiguous uint8 `/distributed_gravity/state` dataset; the full integrity hash consumes the same stream. No schema migration is required. The canonical FNV-1a and SHA-256 digest streams retain the exact previous byte contract. The caller must keep the source state immutable for both passes. A changed length or failed write is a hard error; the existing temporary-file publication and full readback requirements remain in force. This is a bounded serialization path, not a complete owner-local restart-memory certificate.

## M2D owner-local integrity scratch (2026-09-09)

`RestartWritePolicy::memory_governor` optionally supplies the existing process
MemoryGovernor. The production output workflow passes its authority. The
canonical integrity traversal sorts non-owning pointers to temporal-history
records and cells, rather than copying the conserved history. Its arena bound is
`P * sizeof(record_pointer) + Cmax * sizeof(cell_pointer) + 256` bytes, where
P is the history patch count and Cmax the largest history cell count per patch.
The record and cell indexes coexist; other reader and writer allocations are
not included in this narrow bound. The legacy hash and canonical SHA-256
serialization order, restart schema and full verification policy are unchanged.
The public hash helpers retain their one-argument compatibility surface and
accept an optional governor for this scratch. This is not a complete restart
readback admission certificate; full result ownership and verification-time
scheduler/metadata coexistence remain separately tracked M2D gaps.

### DMO exact-equivalence coverage

The internal `restartRuntimeStateExactlyEquivalent(...)` verification contract
includes the authoritative particle drift-history sidecar lanes
`last_drift_time_code` and `last_drift_scale_factor`. These lanes were already
serialized in the restart schema; this is comparator coverage only and does not
change the restart schema version. Mutation-sensitivity coverage in the DMO
restart-equivalence test requires either lane to make exact comparison fail.


## Write-time streaming verification (M48-08)

Checkpoint publication no longer performs a complete restore merely to verify a file that was just written. The production write path uses `verifyRestartCheckpointHdf5()` to compare the current-schema HDF5 checkpoint against live restart truth with bounded hyperslab reads and a fixed 16 MiB/rank CHUÍ-owned workspace. Compact DMO metadata and uniform rung-zero scheduler values are obtained through logical accessors; the gravity force cache is serialized and verified through non-owning spans. `/distributed_gravity/state` is compared through streaming serialization rather than a complete readback string.

This does **not** change actual resume semantics: `readRestartCheckpointHdf5()` remains the complete owning restoration path, including legacy-schema compatibility and restore-candidate MemoryGovernor admission. The restart file schema and integrity algorithms are unchanged.

## P6 optional hierarchical boundary contract

SOURCE-IMPLEMENTED / VALIDATION PENDING. Restart schema v24 adds only the
synchronized endpoint split force (below); no fine/in-flight state. Existing scheduler
bins, next activation, integer tick, integrator/common drift epoch, source
identity, PM cadence and synchronized total-force history suffice **only at a
closed coarse block boundary**. Hierarchy validates aligned ticks, inactive
flags, no pending transitions, next activation equal to current tick, matching
particle/gas scheduler authority, DMO species, and common particle drift epoch.
Fine/in-flight checkpoints are rejected. Existing checksums, streaming verifier,
transactional publication, topology/config hash and collective checks remain.

The existing `/state/particles/time_bin` dataset now writes an optional compact
DMO mirror when present; implicit uniform zero is written only when the lane is
absent. Dataset shape/type are unchanged. Migration carries existing scheduler
and drift fields through its bounded generic transaction, then recompacts
uniform species/owner/drift metadata. Rung-zero compact migration is preserved.

Endpoint force reuse (October 2026 TreePM campaign, R1). A coarse block closes
with an all-active Tree+PM solve at the synchronized positions; the next block
opens at exactly those positions. When every rank still holds that endpoint
complete (every row's Tree generation equals the source generation, PM was
captured at it from the integrator's committed field version, no migration in
between), the block opens with it instead of a second synchronized solve.
Otherwise (first block, after a migration, a pre-v24 restart) the
all-target synchronized bootstrap solve runs, as before. A re-evaluation at the
same positions is not bitwise identical to the endpoint because the adaptive
relative MAC reads the force history, which the endpoint itself updated; v24
therefore persists the endpoint Tree/PM components and nearest-pair lane, and
`GravityRuntime::importRestartForceCache` restores them with the restored source
generation and the persisted PM field version, so uninterrupted and resumed runs
take the same branch with the same forces. Per-row Tree-generation tags remain
transient (re-stamped on import). PM opportunities/versions remain persisted authority. Exact
restart equivalence still needs qualification. Optional refit retains different
transient topology across an uninterrupted run, and spatial EMA starts cold on
restart; compare these modes separately without an unsupported bitwise claim.
New normalized config fields change hashes; readers remain compatible, while
strict workflow hash matching is not relaxed.

## TreePM recovery continuation note — 20261010

This campaign changes no snapshot/restart dataset or schema version. The
hierarchical dispatcher checks the exact admitted all-source row set and the
common synchronized position epoch before opening a block; endpoint publication
and scheduler-owned alignment remain mandatory. The workflow adds a globally
reduced mesh/split displacement restriction to the coarse interval, so nonzero
rung trajectories/step counts can differ from the earlier provisional source.
No force bootstrap is removed, and no new force history is serialized.

The two normalized booleans `numerics.treepm_gaussian_pair_lookup_enabled` and
`numerics.treepm_full_mac_diagnostics` enter the existing config/provenance hash.
Old `.param.txt` files obtain defaults, but normalizing them with this version
changes their hash. Existing strict restart/config compatibility is preserved;
the campaign does not authorize bypassing it for old checkpoints. The immutable
coefficient table is derived runtime cache data and is rebuilt on demand after
restart. Fixed octant traversal changes accumulation order, so old force-vector
checksums are numerical baselines rather than new bitwise replay expectations.

## Constrained normalized-configuration migration (20261010)

20261012 extends this admission contract to `treepm_normalized_defaults_v2`.
Active snapshot/restart/provenance schemas remain v6/v23/v7. A verified,
closed global-KDK predecessor can omit the exact five new default lines
(`treepm_full_residual_lookup_enabled=false`, `treepm_leaf_simd_enabled=false`,
`treepm_radial_coefficients_enabled=false`, block size 64, leaf size 16), with or
without the previous exact two-line default group. The prior two-line-only
projection is retained. Arbitrary missing keys, changed scientific settings,
nondefault accelerators, opaque hashes and unverified payloads do not migrate.
Successful migrations write a v2 receipt; existing v1/v2 receipts survive
subsequent exact resumes. The current normalized text must be reproduced by
the authoritative loader before admission. Exact matching configurations keep
their existing admission path. Derived acceleration caches are never restart
truth. The v1 description below records the original predecessor group and its
unchanged integrity/synchronization requirements. See the
[20261012 handoff](repair/treepm_deep_hotpath_attempt1_closure_20261012.md).

`io::restartConfigurationCompatibilityRecord` is an additive admission helper.
Exact admission now checks normalized text as well as its hash. A mismatched
configuration can use **treepm_normalized_defaults_v1** only when:

- The checkpoint was completely read and its v23 canonical SHA-256 verified by
  `readRestartCheckpointHdf5`; its transient `payload_integrity_verified` receipt
  is not an on-disk field or a substitute for checkpoint validation.
- Original `provenance_v7`, config schema 1, full normalized text, stable hash,
  strong hash, and original state configuration identity all agree.
- It is a closed, restart-safe global KDK checkpoint with rung-zero schedulers.
- Current lookup is false and full MAC diagnostics true. Removing exactly the
  canonical lines for those two default keys from the current normalized text
  produces the original text byte for byte; authoritative parsing restores the
  exact current normalized text. All other settings, including output, topology,
  source and scientific settings, must remain unchanged.

Older versions, opaque hashes, partial/changed historical configurations,
nondefault new options, hierarchical migration and insufficient provenance fail
closed. Ordinary science snapshots remain ineligible for continuation. Existing
ownership, epoch, topology, force-cache and collective admission checks still run.

Raw checkpoint config/provenance remain unchanged. A separately marked record in
new `derived_runtime_state` contains original stable/strong/payload digests,
hex-encoded original normalized text and serialized provenance, and distinct
derived config fingerprints. The workflow report retains the record and later
exact resumes carry it forward. This uses existing provenance serialization and
changes no persistent schema. Migration metadata is bounded to 1 MiB per
normalized representation and 4 MiB original serialized provenance; it allocates
no per-particle state. See the [repair handoff](repair/treepm_p1p5_qualification_and_test_reliability_repair_20261010.md) for positive
production resume tests and refusal fixtures. Qualification is NOT RUN — deferred
by user.

## Independent checkpoint cadence (2026-10-25, untested)

`output.restart_interval_steps = N > 0` writes checkpoints every `N`
completed global steps independently of science snapshots; `0` keeps the
legacy snapshot-coupled policy. The next checkpoint index is derived as the
first multiple of `N` after the restart-persisted step index, so resumed and
uninterrupted runs share one checkpoint timeline without a schema change.
Checkpoints are still written only at restart-safe synchronization
boundaries, and streaming write-time verification is unchanged.

## Code-time checkpoint cadence and PM-aligned deferral (2026-10-25 round 2, untested)

`output.restart_interval_time_code = T > 0` adds a code-time checkpoint
cadence: a checkpoint is due at the first completed global step ending at or
after the next event, and after the write the next event is the first
multiple of `T` (from `numerics.time_begin_code`) strictly after the current
time. The event is derived from restart-persisted time, so no restart schema
field was added; step and time cadences may be combined. With the
experimental `numerics.treepm_pm_reuse_policy = stale_within_cadence`, a due
checkpoint is deferred (counter `io.restart.deferred_for_pm_cadence`) until
the step immediately preceding a PM refresh opportunity, so a resumed run
refreshes the long-range field at the same step as the uninterrupted run;
the pending flags clear only when a checkpoint is actually written.

## Checkpoint request semantics (2026-10-26 audit, untested)

* **Step cadence** (`restart_interval_steps = N`): due at the first completed
  global step whose index is `>= next`, `next` = first multiple of `N` after
  the last checkpoint step (derived on resume from the persisted step index).
* **Code-time cadence** (`restart_interval_time_code = T`): due at the first
  completed global step whose end time is at or after the next event.
  **Timesteps are not shortened** for checkpoint events (science-snapshot time
  events are still clipped as before); a checkpoint is therefore *not*
  exactly aligned with `t_begin + kT`. A step crossing several events raises a
  single request; the next event is re-derived strictly after the write time
  (missed events are not back-filled).
* **Pending vs completed.** A request stays latched (`checkpoint_due`) until a
  checkpoint is actually written at an output- and restart-safe boundary.
  Concurrent step and time events, and events arriving while a request is
  pending, collapse into one checkpoint; both next events are recomputed only
  after the write. Requests are never dropped or duplicated.
* **PM-aligned deferral** (`treepm_pm_reuse_policy = stale_within_cadence`
  only): the pending request is held until the step before a PM refresh
  opportunity. The delay is bounded by `treepm_update_cadence_steps - 1`
  global steps (opportunities advance once per global step in that mode,
  which requires `hierarchical_max_rung = 0`).

## Multirate outer PM cycles (2026-10-26 round 3; IMPLEMENTED — NOT COMPILED — NOT TESTED)

With `numerics.treepm_pm_reuse_policy = multirate_outer_kick`, the coarse-block
boundaries strictly inside a PM cycle are synchronized, but they are **not**
restart-safe or output-safe (`StepBoundaryState::multirate_pm_cycle_interior`),
because momenta carry an open outer PM half-kick. At those boundaries:

* Pending checkpoints and step-cadence snapshots are held until the cycle end
  (counter `io.output.deferred_to_multirate_pm_cycle_end`).
* Analysis hooks are skipped.
* Code-time snapshot events are not reached mid-cycle, because the cycle
  length is clipped exactly to the next event or `t_end`.

**Soft stop (round 3).** The segment step budget no longer shortens the PM
cycle. If the budget is reached while a cycle is open, the segment keeps
running to the cycle end. The extra blocks are reported in
`ReferenceWorkflowReport::multirate_soft_stop_extension_blocks` and counter
`time.hierarchical_multirate_soft_stop_extension_blocks`, alongside
`requested_steps` and `completed_steps`. A segment that ends with a cycle open
is a `logic_error`.

Every checkpoint therefore closes a cycle. A resumed segment opens its next
cycle exactly as the uninterrupted run does, from state that is already
persisted:

* the scheduler;
* the integrator, including the previous quantum;
* the v24 force endpoint;
* pending output.

There is no persisted interior-cycle state and **no restart schema change**.
Restarting inside a cycle is unsupported by design. A restart whose integrator
state is not restart-safe is refused when a multirate segment starts.

The transactional cycle rollback (see `docs/time_integration.md`) is
in-memory only. It never writes or reads a checkpoint.

## CHUI Part 1 normalized-config migration (20261103)

Part 1 adds nine normalized keys (`numerics.hydro_reconstruction`,
`hydro_slope_limiter`, `hydro_integration_scheme`, `hydro_riemann_solver`,
`hydro_dual_energy`, `hydro_dual_energy_eta`, `parallel.node_memory_admission`,
`node_memory_budget_bytes`, `node_memory_reserve_bytes`). Resuming an
integrity-verified v23/v24 checkpoint written before Part 1 goes through the
existing exact-default projection in `restartConfigurationCompatibilityRecord`:
the nine lines are removed only if every one of them holds its default value
in the current configuration; any other value refuses the resume. No HDF5
schema version changes. Reproducibility impact: a resumed run continues with
the corrected numerics (F1 predictor sign and `1/a`, F3 passive lanes,
periodic box faces evaluated once), so it is not bit-identical to a
continuation under the pre-Part-1 binary.

## CHUI Part 1 stabilization: hydro timestep rule (20261103)

The gas CFL candidate now follows the stability rule of the integration
family: dimension-summed for `muscl_hancock` (default), per-axis for `ctu`
(`docs/hydro_numerics.md` section 4). Resuming any checkpoint with the default
method yields smaller multi-dimensional timesteps than the binary that wrote
it; no schema or config migration is involved, and the old per-axis rule
cannot be selected for MUSCL-Hancock because it is linearly unstable for
diagonal 3D flow. Reproducibility against pre-stabilization runs is not
bitwise.
