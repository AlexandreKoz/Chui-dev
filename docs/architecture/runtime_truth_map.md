# Runtime truth map (Stage 0 freeze audit)

20261012 additive interface clarification: drift owns position/velocity
read/write through its typed view and reads migration ownership. Registry
`descriptors()` exposes an immutable inspection span, invalidated by registry
mutation/destruction, without granting resource access. TreePM bounds,
coefficient tables and worker staging are derived coordinator-owned state;
none become particle, scheduler or restart authority. Production serial memory
estimation may evict optional acceleration caches before mandatory admission;
the additive `trimOptionalAccelerationCaches` API must never overlap a solve.
Tree generation/layout/frame checks and PM/source epochs remain authoritative.
See the [source/lifetime contract](../repair/treepm_deep_hotpath_attempt1_closure_20261012.md).

_Date: 2026-04-25 (UTC)_

## Scope and method

This document is a code-first runtime-state ownership audit. It was compiled from direct repository inspection and command-backed build/test checks, before runtime-state implementation edits.

### Commands used for discovery

- `rg --files -g 'AGENTS.md'`
- `rg --files`
- `rg -n "time_bin|timestep|active_set|scheduler|reorder|species|sidecar|softening|restart|checkpoint|snapshot|normalize|derived|provenance|migration|compact|view|cache|mirror|cell_id|gas" src include tests docs`
- targeted `sed -n` reads of core/workflow/io/config/provenance files listed below

### Build/test commands executed

- `cmake --preset cpu-only-debug`
- `cmake --build --preset build-cpu-debug`
- `ctest --preset test-cpu-debug --output-on-failure`

## Exact files inspected

### Core runtime state and scheduling

- `include/cosmosim/core/simulation_state.hpp`
- `src/core/simulation_state.cpp`
- `src/core/simulation_state_structures.cpp`
- `src/core/simulation_state_ownership.cpp`
- `src/core/simulation_state_species.cpp`
- `src/core/simulation_state_active_views.cpp`
- `include/cosmosim/core/time_integration.hpp`
- `src/core/time_integration.cpp`
- `src/core/simulation_state_metadata.cpp`

### Runtime orchestration / mutation authority

- `src/workflows/reference_workflow.cpp`
- `include/cosmosim/workflows/reference_workflow.hpp`

### Config / normalization / derived values

- `include/cosmosim/core/config.hpp`
- `src/core/config.cpp`
- `include/cosmosim/core/simulation_mode.hpp`
- `src/core/simulation_mode.cpp`

### Restart/checkpoint/snapshot/provenance

- `include/cosmosim/io/restart_checkpoint.hpp`
- `src/io/restart_checkpoint.cpp`
- `include/cosmosim/io/snapshot_hdf5.hpp`
- `src/io/snapshot_hdf5.cpp`
- `include/cosmosim/core/provenance.hpp`
- `src/core/provenance.cpp`
- `src/io/internal/io_contract.cpp`

### Tests inspected for coverage mapping

- `tests/unit/test_simulation_state.cpp`
- `tests/unit/test_time_integration.cpp`
- `tests/unit/test_species_state_organization.cpp`
- `tests/unit/test_hot_cold_sidecar_layout.cpp`
- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/unit/test_snapshot_hdf5_schema.cpp`
- `tests/unit/test_units_cosmology_provenance.cpp`
- `tests/integration/test_hierarchical_time_bins.cpp`
- `tests/integration/test_hierarchical_timestep_regression.cpp`
- `tests/integration/test_time_integration_loop.cpp`
- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_species_migration_invariants.cpp`
- `tests/integration/test_soa_species_pipeline.cpp`
- `tests/integration/test_species_mixed_state.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`
- `tests/integration/test_provenance_roundtrip.cpp`
- `tests/integration/test_parallel_two_rank_restart.cpp`

---

## Domain map

## 1) Timestep and timestep-bin state

### Storage locations (current truth)

- Scheduler authoritative bin timeline: `HierarchicalTimeBinScheduler` internals (`m_hot`, `m_elements_by_bin`, `m_active_elements`, `m_current_tick`).
- Scheduler-owned PM cadence truth: `PmSynchronizationState` (`gravity_kick_opportunity`, cadence steps, field version, last refresh opportunity, field-built step/scale factor).
- During an authorized refresh but before cadence commit, the core integrator
  emits rank-local refresh evidence in `PmRefreshDirective`. Distributed
  execution globalizes that evidence in `TimeCoordinator` before materializing
  one identical transient decision on every rank. Operational events emitted in
  that interval copy its opportunity/version/build epoch; the previous
  committed `PmSynchronizationState` remains authority for completed refreshes.
- Per-particle timestep bin mirror: `SimulationState::particles.time_bin`.
- Per-cell timestep bin mirror: `SimulationState::cells.time_bin`.
- Integrator coarse state and metadata mirror: `IntegratorState::{current_time_code,current_scale_factor,dt_time_code,step_index,scheme,time_bins}`.

### Mutation paths

- Scheduler mutates bin/tick truth via `reset`, `setElementBin`, `requestBinTransition`, `beginSubstep` + `endSubstep`.
- Runtime bin adaptation submits labeled criteria candidates in `updateAdaptiveTimeBins(...)`; scheduler reconciliation owns pending-transition mutation and commit.
- State mirrors refreshed from scheduler in `syncTimeBinMirrorsFromScheduler(...)`.
- Mirror drift detection is now explicit via `timeBinMirrorsMatchScheduler(...)` / `debugAssertTimeBinMirrorAuthorityInvariant(...)`.
- Restart load mutates integrator/scheduler via `readRestartCheckpointHdf5(...)` and `importPersistentState(...)`.

### Read paths

- Active element extraction for stepping: `scheduler.beginSubstep()`.
- Stage execution uses `ActiveSetDescriptor` built from scheduler active elements, stamped with scheduler tick/state-generation provenance.
- Diagnostics/provenance read bin distribution via scheduler diagnostics and metadata payload fields.

### Serialization paths

- Restart serializes:
  - `IntegratorState` scalar fields and `time_bins` fields.
  - Full `TimeBinPersistentState` arrays in `/scheduler/*`.
  - Particle/cell `time_bin` arrays under `/state`.
- Snapshot currently writes/reads particle data but imported particles are assigned `time_bin = 0` during read; no persistent timestep-bin continuity in snapshot I/O.

### Cache/mirror/view classification

- `SimulationState::{particles.time_bin,cells.time_bin}` is a mirror of scheduler truth during workflow execution (rebuilt from scheduler after scheduler commit).
- `CellSoa::time_bin` is a gas-cell scheduler mirror keyed through `GasCellIdentityMap::gas_cell_id`; parent-particle mapping is legacy compatibility only.
- `IntegratorState::time_bins` is metadata-level mirror and not sufficient alone to reconstruct scheduler occupancy/membership.

### Reorder/resize/migration behavior

- Particle reorder/migration physically carries particle `time_bin` as a mirror with moved rows, but exact scheduler continuation requires scheduler remap/rebuild helpers keyed by stable particle identity.
- Gas-cell and AMR patch migration carries `time_bin` only as a diagnostic/compatibility mirror; gas-cell scheduler identity records keyed by `gas_cell_id` are required for authoritative continuation.

### Restart/reload behavior

- Restart restores exact scheduler persistent vectors and verifies payload integrity.
- Restart validation rejects stale particle and parent-particle-mapped cell `time_bin` mirrors before import; successful reads rebuild mirrors from scheduler authority.
- Post-read scheduler can be reconstructed exactly (subject to schema compatibility and integrity checks).

### Ownership decision

- **Single live authority**: Stage 2 assigns timestep-bin membership, activation, active-set construction, and committed PM cadence ownership to scheduler objects (`HierarchicalTimeBinScheduler` and `PmSynchronizationState`). The core integrator may emit rank-local PM-refresh evidence, but only the MPI-aware `TimeCoordinator` may convert that evidence into a distributed directive before branch-dependent TreePM work. State `time_bin` lanes, migration records, restart mirrors, and `IntegratorState::time_bins` are mirrors/metadata only. Any future path that consumes mirrors as scheduling authority must be rejected or documented as a new ADR-backed schema/interface change.
- `gravity.pm_long_range_field` is an observer, not a second cadence owner. Its
  in-flight payload is copied from `PmRefreshDirective` because the matching
  `PmSynchronizationState` commit has not happened yet.

### Tests covering this domain

- `tests/unit/test_time_integration.cpp`
- `tests/integration/test_hierarchical_time_bins.cpp`
- `tests/integration/test_hierarchical_timestep_regression.cpp`
- `tests/integration/test_time_integration_loop.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`

### Remaining tests (audit-identified)

- Snapshot read/write test asserting intentional loss (or future preservation) of timestep bins, to avoid accidental contract drift.

---

## 2) Active-set construction

### Storage locations

- Authoritative active indices: `HierarchicalTimeBinScheduler::m_active_elements` (transient per substep).
- Step-local descriptor: `ActiveSetDescriptor` in orchestrator calls.
- Workspace compact buffers: `TransientStepWorkspace` and active views (`ParticleActiveView`, `CellActiveView`, `GravityParticleKernelView`, `HydroCellKernelView`).

### Mutation/read paths

- Mutated in scheduler `rebuildActiveSet()` and workflow active-set split loop.
- Read by stage callbacks through `StepContext.active_set`.
- Gravity/hydro compact views gather/scatter from state arrays per active indices.

### Serialization

- Active-set vectors are not serialized directly.
- Restart serializes scheduler persistent state, from which active sets can be recomputed.

### Cache/mirror/view paths

- Active views are compact transient materializations (gather buffers) with explicit scatter back to state.

### Reorder/resize/migration behavior

- Active sets are index-based and therefore sensitive to reorder/migration operations; downstream paths rely on updated indices from reorder/migration routines.
- Active eligibility is species-agnostic in the scheduler contract; migration changes tags/sidecars but does not itself change active membership unless scheduler bins are reassigned through scheduler APIs.

### Restart/reload behavior

- Active sets are recomputed after restart via scheduler state rather than loaded as persisted lists.

### Ambiguity/duplication

- No second persistent active-set truth detected; active views are clearly transient.

### Tests covering this domain

- `tests/unit/test_time_integration.cpp`
- `tests/integration/test_time_integration_loop.cpp`
- `tests/integration/test_hierarchical_time_bins.cpp`

### Missing tests

- Direct tests for `build*/scatter*ActiveView` stale-index fault handling under reordered/migrated state.

---

## 3) Particle and species ordering

### Storage locations

- Particle SoA arrays in `SimulationState::particles`.
- Ordering keys/metadata in `SimulationState::particle_sidecar` (`particle_id`, `sfc_key`, `species_tag`, `owning_rank`, flags).
- Species index maps in `SimulationState::particle_species_index`.
- Species counts in `SimulationState::species.count_by_species`.

### Mutation paths

- Full reorder through `reorderParticles(...)` with `ParticleReorderMap`.
- Species map rebuild via `SimulationState::rebuildSpeciesIndex()`.
- Migration commit (`commitParticleMigration`) compacts/remaps and then rebuilds species counts/index.

### Read paths

- Species-local loops use `particle_species_index.globalIndices(...)`.
- Workflows and I/O read sidecar tags/counts for classification and packing.

### Serialization

- Restart serializes all particle arrays, particle sidecar lanes, and species counts.
- Snapshot serializes particle content by type groups; species re-derived on read from type mapping.

### Cache/mirror/view

- `ParticleSpeciesIndex` is derived/cache-like index map rebuilt from sidecar tags.

### Reorder/resize/migration behavior

- Reorder supports sidecar policies (`kMoveWithParent` or remapped parent indirection).
- Migration commit rewrites particle arrays and sidecar rows, then recomputes species ledger/index.

### Restart/reload behavior

- Restart reads species counts and sidecars, then calls `rebuildSpeciesIndex()` + invariant validation.

### Ambiguity/duplication

- Duplicated species truth exists intentionally:
  - `species.count_by_species` (ledger),
  - `particle_sidecar.species_tag` (row truth),
  - `particle_species_index` (derived index).
- Ownership invariant checks enforce consistency but mutation authority is spread across multiple functions.

### Tests covering

- `tests/unit/test_simulation_state.cpp`
- `tests/unit/test_species_state_organization.cpp`
- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_soa_species_pipeline.cpp`
- `tests/integration/test_species_mixed_state.cpp`

### Missing tests

- Fault-injection coverage for partial updates where species ledger changed without index rebuild.

---

## 4) Species sidecars

### Storage locations

- `SimulationState::{star_particles, black_holes, tracers}` as index-addressed sidecars keyed by `particle_index`.

### Mutation paths

- Reorder remaps sidecar parent indices according to policy.
- Migration commit rebuilds each sidecar using kept rows plus inbound records.
- Physics callbacks mutate sidecar payload values (stellar evolution, BH, tracer updates).

### Read paths

- Physics kernels and diagnostics read sidecars by species.
- Migration pack path reads required sidecar rows by parent particle index.

### Serialization

- Restart serializes/deserializes all sidecar fields (including stellar channel arrays).
- Snapshot serializes selected sidecar-derived datasets where schema supports them (e.g., tracer host metadata, softening sidecar).

### Cache/mirror/view

- No separate persistent mirror found; sidecars are canonical for species-specific cold data.

### Reorder/resize/migration behavior

- Explicit sidecar sync policies and stale-index guards exist.
- Migration can drop/rebuild sidecar rows as ownership changes.

### Restart/reload behavior

- Restart reads sidecars then validates one-to-one invariants vs species-tagged particles.

### Ambiguity/duplication

- Parent link duplication: species tag lane and sidecar membership both encode species identity.

### Tests covering

- `tests/unit/test_simulation_state.cpp`
- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`

### Missing tests

- Snapshot roundtrip completeness assertions for all sidecar families (currently strongest guarantees are restart-focused).

---

## 5) Gas cell identity

### Storage locations

- Gas cell structural lanes: `SimulationState::cells`.
- Hydro/cold lanes: `SimulationState::gas_cells`.
- Identity authority: `SimulationState::gas_cell_identity`, keyed by stable `gas_cell_id` with optional `parent_particle_id` compatibility/provenance.

### Mutation paths

- Eulerian/AMR gas-cell centres are geometry authority and are not drifted through generic gas compatibility particles.
- Production workflow compaction/migration rebuilds cell arrays from gas-cell/AMR patch payloads keyed by stable `gas_cell_id` and `patch_id`; legacy gas-particle keyed helpers are compatibility-only.
- Source transactions refresh compatibility particle mass/velocity mirrors only for legacy I/O/adapters; those mirrors do not control hydro or gas gravity.
- Host cell indices for BH/tracers are remapped after rebuild.

### Read paths

- Hydro kernels read from cell/gas arrays (often through active compact views).
- Gravity consumes each owned authoritative leaf gas cell directly at its cell centre and scatters acceleration back to that cell row. Generic gas compatibility particles are excluded from gravity source mass.
- Migration and diagnostics may read gas-particle↔cell lineage correspondence where compatibility/provenance is required.

### Serialization

- Restart serializes `cells` and `gas_cells` fully.
- Snapshot focuses on particle-centric schema; gas cell arrays are not a first-class persisted group there.

### Cache/mirror/view

- `parent_particle_id` and generic gas-particle rows are compatibility/provenance mirrors, not cell-position, hydro-mass, or gravity authority.
- Gravity force-cache validity tracks particle-row, cell-row, and gas-identity generations; refinement/migration invalidates stale cell accelerations.

### Reorder/resize/migration behavior

- Gas identity continuity across migration is preserved via particle IDs, not raw cell indices.

### Restart/reload behavior

- Restart restores explicit cell arrays and patch indices.

### Ambiguity/duplication

- **Ambiguous identity authority**: gas identity is partly index-based (cell index), partly particle-ID-based (migration rebuild), and partly species-order-based (gas global index list).

### Tests covering

- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_hydro_sod_like.cpp`
- `tests/integration/test_amr_static_refinement_sync.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`

### Missing tests

- Explicit invariant tests for gas cell identity continuity across reorder + migration + restart combined path.

---

## 6) Softening state

### Storage locations

- Config-level defaults/species lanes in `SimulationConfig::numerics`.
- Optional per-particle softening value sidecar: `particle_sidecar.gravity_softening_comoving`; authoritative override membership is `particle_sidecar.has_gravity_softening_override`. A value without a mask bit is a materialized default/diagnostic mirror, not override authority.
- Provenance lanes (`ProvenanceRecord` softening policy/kernel/epsilon).

### Mutation paths

- Initial population from species policy in `maybeInitializeParticleSofteningFromSpeciesPolicy(...)`.
- Reordered/migrated alongside particle rows.
- Snapshot import may populate sidecar if `GravitySofteningComoving` dataset exists.

### Read paths

- Time-step criteria gravity hook uses sidecar if present, else species/global fallback.
- Gravity modules consume policy and sidecar values for force/timestep behavior.

### Serialization

- Restart serializes full sidecar when present.
- Snapshot writes per-particle softening additively and reads it if present.

### Cache/mirror/view

- Sidecar is canonical runtime lane for per-particle softening overrides; config values are policy defaults.

### Reorder/resize/migration behavior

- Reorder/migration preserve softening lane with particle rows.

### Restart/reload behavior

- Restart restores sidecar or empty state (policy-only mode).

### Ambiguity/duplication

- Duplicate representational lanes exist: global/species config policy plus optional per-particle sidecar overrides.

### Tests covering

- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`

### Missing tests

- Coverage for conflict resolution when config species policy and snapshot softening sidecar disagree.

---

## 7) Config-derived runtime values

### Storage locations

- Canonical typed config: `FrozenConfig::config`.
- Canonical normalized text/hash: `FrozenConfig::{normalized_text, provenance.config_hash_hex}`.
- Runtime derived policy objects built downstream (`ModePolicy`, per-module derived values).

### Mutation paths

- Parsing + normalization in `normalizeValidateFreeze(...)`.
- Compatibility alias/scalar lanes collapsed into axis-aware typed fields.
- Validation enforces bounds/required relationships.

### Read paths

- Reference workflow and modules consume typed config.
- Provenance generation and output schemas consume normalized hash + selected derived fields.

### Serialization

- Normalized config text written to run directory and embedded in snapshot/restart attributes.
- Hash continuity checked by shared IO contract validation.

### Cache/mirror/view

- Normalized text acts as canonical serialized mirror of typed config.

### Reorder/resize/migration behavior

- Not applicable directly; config drives policies that influence these paths.

### Restart/reload behavior

- Restart requires normalized text/hash/provenance agreement.

### Ambiguity/duplication

- Intentional backward-compatibility duplication exists (e.g., scalar and axis-aware PM/box keys); normalized output collapses to canonical axis-aware form.

### Tests covering

- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_simulation_mode.cpp`
- `tests/integration/test_config_examples.cpp`
- `tests/integration/test_simulation_mode_toy_runs.cpp`

### Missing tests

- Cross-check tests asserting derived runtime constants used by workflow callbacks remain in lockstep with normalized config text for every mode.

---

## 8) Restart and checkpoint truth

### Storage locations

- Authoritative restart payload: `RestartWritePayload` (state + integrator + scheduler + provenance + normalized config + distributed gravity state).
- On-disk groups: `/state`, `/integrator`, `/scheduler`, `/distributed_gravity`, plus root hash/schema attributes.

### Mutation paths

- Restart writer assembles payload, computes integrity hash, writes temp file, finalizes atomically.
- Reader reconstructs structures, re-computes integrity hash, validates schema + invariants.

### Read paths

- Workflow restart path loads restart payload then resumes scheduler/integrator/state.

### Serialization overlap with snapshot

- Overlap domains: particles, IDs, softening sidecar (optional), normalized config, provenance.
- Non-overlap: restart carries scheduler persistent state, module sidecars, strict continuation integrity hash.

### Cache/mirror/view

- No separate restart cache detected; restart payload integrity hash is the main consistency guard.

### Reorder/resize/migration behavior

- Restart stores post-mutation state as-is; no special reorder metadata beyond arrays/index fields and sidecars.

### Restart/reload behavior

- Exact continuation target with strict schema compatibility and integrity checks.

### Ambiguity/duplication

- Snapshot and restart both carry normalized config and provenance but with different completeness guarantees; this is an intentional dual-contract surface.

### Tests covering

- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_parallel_two_rank_restart.cpp`
- `tests/integration/test_provenance_roundtrip.cpp`

### Missing tests

- Contract test asserting snapshot->restart->snapshot roundtrip equivalence boundaries (what is intentionally lossy vs exact).

---

## Cross-domain duplication/mirror inventory (classified)

- `scheduler hot metadata` vs `state.*.time_bin`: **mirror** (scheduler-authoritative during loop).
- `species_tag` vs `species.count_by_species` vs `ParticleSpeciesIndex`: **canonical + derived + ledger duplication** with explicit invariant checks.
- Gas identity by `cell index` vs `gas species global index` vs `gas particle_id`: **ambiguous multi-lane identity** (currently reconciled during compaction via particle IDs).
- Softening policy config vs optional per-particle softening sidecar: **policy + override dual truth**.
- Normalized config hash in metadata/provenance/restart attributes: **redundant consistency mirrors** validated by contract checks.

## Repair tickets / follow-up prompts (no broad code change in this stage)

1. **P31-STATE-TEST-DRIFT-026 (Open)**
   - Symptom: `unit_snapshot_hdf5_schema` currently asserts `gadget_arepo_v2` while runtime schema reports `gadget_arepo_v4` during `ctest --preset test-cpu-debug`.
   - Scope: test/schema contract alignment only.

2. **P32-STATE-OWNER-AMBIGUITY-027 (Open)**
   - Scope: make mutation authority explicit for scheduler bins vs state bin mirrors and gas cell identity lanes; add invariant checks/tests before behavior changes.

3. **P33-STATE-SNAPSHOT-RESTART-BOUNDARY-028 (Open)**
   - Scope: codify and test intentional overlap/loss boundaries between snapshot and restart (especially timestep bins, scheduler state, module sidecars).

4. **P34-STATE-SOFTENING-OVERRIDE-029 (Open)**
   - Scope: explicit policy for precedence/compatibility between config species softening and imported per-particle softening sidecar values.

## Build/test status for this audit

- Configure/build with `cpu-only-debug` preset succeeded.
- Full `ctest --preset test-cpu-debug --output-on-failure` did **not** complete cleanly:
  - `unit_snapshot_hdf5_schema` failed on schema name assertion (`gadget_arepo_v2` expected).
  - Test run subsequently stalled in this environment after progressing into later integration tests, so this audit records a partial run with explicit blocker instead of claiming green closure.

## Reproducibility impact statement

This stage is documentation/audit only. No runtime behavior, schema payload, deterministic scheduling mode, config normalization flow, or provenance completeness path was modified by this change.

### Stage 0 follow-up note: migration softening value vs override authority

Migration payloads now distinguish `has_gravity_softening_value` from `has_gravity_softening_override`. The former preserves a materialized numeric value when present; the latter is the only flag that carries per-particle override authority across migration.

## Stage 0 P0-05..P0-08 repair update

The P0-05..P0-08 repair pass adds the following concrete runtime-truth structures:

- `GasCellSidecar::gas_cell_id` and `GasCellSidecar::parent_particle_id` are compatibility mirror lanes for the authoritative `GasCellIdentityMap`. They are serialized in restart/checkpoint state and included in restart integrity hashing; parent ID may be absent.
- `ParticleSidecar::has_gravity_softening_override` is the explicit owner of per-particle softening override authority. Numeric softening values without a corresponding override bit are materialized values, not independent override truth.
- `DerivedRuntimeConfig` is the owner for normalized-config-derived runtime constants.
- `ActiveSetDescriptor` now carries scheduler provenance and state generation metadata; stale descriptors fail through `debugAssertActiveSetDescriptorFresh()` before solver consumption.

Remaining caveat: some legacy physics/unit tests still operate standalone hydro cells without full workflow AMR ownership. Those states are deliberately not promoted to distributed production coverage unless the state ownership layer installs stable gas-cell identity records and validates dense-row coverage.

## 2026-07-16 Campaign A runtime-composition addendum

| Runtime fact | Live authority | Transient view/mirror | Invalidation or rebuild rule |
| --- | --- | --- | --- |
| Gas derived thermodynamics (`pressure_code`, `sound_speed_code`, `temperature_code`) | `gas_cells.density_code`, `internal_energy_code`, peculiar velocity, metal mass, plus the configured closure | the three derived cache lanes | Hydro primitive store rewrites P/c_s; any other mutation (feedback today) calls `notifyGasThermodynamicMutation()` then `hydro::refreshGasThermodynamics()`; `TimeCoordinator` refuses stale caches (epochs are rank-local and non-persistent). |
| Hydro face states, CTU transverse increments, dual-energy auxiliary energy | staged `HydroConservedStateSoa` of the advance | `HydroScratchBuffers` lanes (batch-scoped) | Rebuilt every advance from step-start state; never restart truth; only the final-pass fluxes update conserved state and flux registers. |
| Particle/gas scheduler and integrator state | `workflows::RungZeroTimeState` | particle/cell `time_bin` lanes are mirrors | Fresh construction or restart import initializes the owner; mirrors refresh only from its schedulers. |
| Production stage ordering and contributions | Frozen `workflows::RuntimeExecutionPlan` built by `RuntimeModuleRegistry` | descriptor declarations and typed task functions | Registry freezes once per run. Adding/removing a descriptor adds/removes its task without a runner callback-list edit. |
| Stage resource freshness and task grant | `SimulationState` particle/cell/gas-identity generations, particle scheduler tick, `IntegratorState::step_index`, and the frozen task declaration | `RuntimeResourceLease` plus one task-scoped `RuntimeResourceAccess` span captured in a typed stage view | Rebuild for every stage and rebind for every task. Reorder/migration/compaction or tick/step advance makes guarded leases stale; the task grant is cleared after invocation. Views are never serialized and expose no `StepContext`. |
| Active-set construction | `RungZeroTimeState` schedulers during an open substep | `ActiveSetDescriptor` and typed task views | Rebuild after scheduler commit or ownership-generation change; never repair an old view locally. |
| Output cadence and legal boundary request | `RungZeroTimeState::pendingOutput()` plus committed integrator boundary metadata | `OutputRestartStageView` | Cadence is restored from checkpoint state or initialized once, latched before a step, and consumed only at a safe boundary. |

Production workflow tasks no longer register `IntegrationCallback` objects or
receive `StepContext` through their public interface. External/custom tasks can
observe only freshness and their descriptor grant. Source-private built-in owner
access validates that grant before entering legacy rung-zero numerical bodies.
Numerical modules remain below workflow composition in the dependency graph,
and owner services include other owners only through public/narrow contracts.

## 2026-07-13 gravity runtime-truth addendum

The production gravity hardening resolves these additional authority
boundaries:

| Runtime fact | Live authority | Mirrors/transient state | Refresh or invalidation rule |
| --- | --- | --- | --- |
| Physical Newton constant in code units | Frozen `UnitSystem`, converted once through `core::newtonGravitationalConstantCode` | PM/tree option fields and the `gravity.treepm_setup` operational-event value | Re-derive from frozen length/mass/velocity units at workflow construction or restart; never hard-code unity in a physical-unit workflow. |
| Cosmological force time dependence | `CosmologicalTimeline`; during a step, `StepContext.timeline_step`; after commit, `IntegratorState::{current_scale_factor,current_hubble_rate_code}` | `PmSolveOptions::scale_factor` and PM field-build scale are validity/source-time metadata; `HydroUpdateContext` is a stage-local view | PM and short-range TreePM return scale-free `A`. Collisionless KDK applies its Hubble-drag/kick factors. Post-drift gas hydro uses `timeline_step.scale_factor_end/hubble_end_code` for `A/a^2`, `-H m`, and `-H(2K+3P)` because committed `IntegratorState` is still at step begin. Hydro callbacks must not recompute SI `H(a)`. |
| Particle decomposition epoch | `GravityStageCallback::m_decomposition_epoch` | Provenance/restart fields and TreePM hierarchy/request/response record epochs | Advance once, on every rank, only after an actual committed ownership change; restore before resumed TreePM use. Rank-local dense-row generations are not this authority. Drift does not advance this epoch. |
| Gravity source generation | Authoritative `SimulationState::gravitySourceGeneration` (bumped on source mutation including drift) | `TreePmOptions.source_generation`, PM long-range validity signature, force-cache/source snapshot metadata | Advance on every source-state mutation. PM field reuse and top-domain geometry freshness both require equality with the live generation; decomposition epoch alone is never a substitute. |
| Authoritative top-domain leaf geometry | Stable seed leaves (`GravityRuntime::m_authoritative_top_domain_seed_leaves`) from `MigrationBalanceRuntime::authoritativeTopDomainLeaves` at segment/restart/ownership-change install; per-solve O(N) bound refit in `GravityRuntime` producing the published set stamped with the current source generation | `TreePmOptions.authoritative_domain_leaves` + `authoritative_geometry_source_generation` (derived routing view; not restart truth) | Install replaces seed + published together. Ordinary refit always consumes the seed set (never the previous published result), so empty-group omission cannot erode SFC partition identity; published leaves may omit currently-empty groups. Refit after compact source rebuild before each solve when the stamp is stale; TreePM uses authoritative leaves only when generation matches **and** local sources are covered, else falls back to the local tree-root packet with an explicit fallback reason. `commitParticleDecompositionChange()` invalidates freshness immediately on ownership epoch change; reinstall follows. Seed rebuild only after ownership change or restart; drift changes source generation, never the decomposition epoch or the seed set. Freshness token is source generation — revisit for future mixed-rung KDK prediction/evaluation epochs. |
| PM long-range field validity | `TreePmCoordinator::LongRangeFieldValidity` | FFT plans, mesh arrays, halos, and diagnostic build metadata | Reuse requires an explicit physical source generation plus force epoch, build scale, `G_code`, split/box/assignment/boundary/decomposition-mode/deconvolution compatibility on every rank. Missing/incompatible reuse fails. Particle decomposition alone does not invalidate a fixed-slab PM field. |
| Dense-row gravity force history | `GravityStageCallback` committed force cache plus per-row validity | Relative-MAC previous-acceleration magnitudes and restart `/gravity_force_cache` | Invalidate immediately on a particle ownership/index-generation change; rebuild all locally owned rows before a kick when required. A post-migration checkpoint records invalid rather than stale cache truth. |
| Gravity timestep acceleration | Post-step adaptive scheduler criteria using committed `IntegratorState.current_scale_factor` | TreePM force lanes store scale-free `A`; `ComovingGravityTimeStepInput` is the transient public view for `eps_com`, `|A|`, and `a` | `computeComovingGravityTimeStep` validates finite positive `a`, converts to comoving-coordinate acceleration `|A|/a^3`, then evaluates `dt_grav=eta sqrt(eps_com/(|A|/a^3))`. Do not pass `A` or `A/a^2` with `eps_com` to the generic coordinate-neutral helper. |
| Production timestep rung | Frozen config plus particle/gas schedulers constrained to bin zero | State `time_bin` lanes remain scheduler mirrors | `hierarchical_max_rung` must be zero; restart validation requires scheduler `max_bin=0`, committed bins zero, and pending bins zero/unset. Mixed-rung production is blocked until per-element kick/drift epochs exist. |

Production cadence one means each integrator-issued rank-coordinated
force-refresh surface rebuilds PM. The scheduler core can still represent and
test coordinated local boundaries, and the coordinator retains an explicit
lower-level reuse seam, but neither is a competing production authority while
`hierarchical_max_rung=0`.

All PM density/request/response and TreePM hierarchy/request/response payloads
are explicit versioned little-endian wire records. Native C++ object layout is
not wire truth. Decode validation and rank-coordinated failure voting happen
before the next collective boundary; zero-record ranks still participate.

## 2026-08-29 M1B memory-governor addendum

| Runtime fact | Live authority | Transient view/mirror | Invalidation or rebuild rule |
| --- | --- | --- | --- |
| Process runtime memory admission | One `core::MemoryGovernor` owned by the reference-workflow composition root | Borrowed pointer in `RuntimeServices`; immutable `MemoryGovernorSnapshot` in the existing memory/profiler report | Construct once from `parallel.process_memory_budget_bytes` and existing process reserve/margin keys. Reconcile baseline owned capacity after authoritative memory reports change; never introduce subsystem-local governors. |
| Governed retained scratch capacity | `MonotonicScratchAllocator::Block` backing storage plus its committed `MemoryReservation` | Logical scratch offsets/current-use counters | Reserve before block allocation, commit after allocation, retain reservation for the block lifetime. `reset()` changes logical use only and must not release physical commitment. |
| Stage capability/freshness | Existing `RuntimeResourceLease` and runtime-resource epoch machinery | Task-scoped typed stage views | Unchanged by M1B. Memory reservations and stage-view leases are separate authorities and may coexist. |

`MemoryReport` remains the ownership/capacity tree. `MemoryGovernor` is the
admission authority. Governed entries are excluded from reconciled baseline
owned bytes because their physical capacity is already represented by the
governor committed counter; the two surfaces must not be summed as independent
physical allocations. See `docs/memory_governance.md`.

## 2026-09-06 M2D runtime scheduling authority addendum

| Runtime fact | Live authority | Transient policy/view | Invalidation or release rule |
| --- | --- | --- | --- |
| Process byte admission and deterministic batch headroom | One composition-root `core::MemoryGovernor` | `MemoryGovernorSnapshot` plus `DeterministicBatchSizingPolicy/Result` | Re-evaluate from a fresh snapshot at each task/batch decision; a sizing result never owns bytes and does not replace the required reservation for actual residency. |
| Major production task order and dependency graph | Frozen `workflows::RuntimeExecutionPlan` from `RuntimeModuleRegistry` | `RuntimeTaskSchedulingProfile`, dependency keys, and `runtimeTasksMayOverlap()` | Registry freeze validates dependencies against the already-authoritative order. The current dispatcher remains serial; scheduling metadata is not restart/scientific state. |
| Per-rank decomposition safety | Existing Morton/SFC decomposition and its `DecompositionConfig` | Hard `max_rank_memory_bytes`, transient reserve, expanded work-component telemetry | Recompute from current ownership and fresh process-governor headroom before rebalance. Any memory-unsafe target is rejected before migration; rank maximum, not mean, is the safety quantity. |
| Node/rank/thread topology evidence | `parallel::MpiContext` plus OpenMP runtime configuration | `runtime.topology` profiler event | Rebuilt at process/runtime initialization. It is diagnostic/tuning evidence, not a second scheduling or memory authority. |
| Deferred optional science diagnostics | `AnalysisRuntime` owner lifetime | Two constant-size coalesced cadence records with first/latest due epoch and missed count | Record due epochs under Red/Trip; clear only after successful execution at a later safe analysis hook or explicit terminal drop. Checkpoint provenance records pending metadata but no restart authority is created. Never override the hard memory ceiling. |

The M2D task graph does not supersede `RungZeroTimeState` scheduler authority or
change accepted kick/drift/hydro/source ordering. It describes resource lifetime
and overlap legality around the existing scientific schedule. The decomposition
hard cap likewise constrains ownership placement; it does not become a second
whole-process memory ledger and must not be summed independently with governor
baseline/commitment accounting.

## M2D-1 task and optional cadence ownership

The existing `MemoryGovernor` is the only budget authority. The runtime registry owns task declarations and deterministic dispatch; dispatcher-owned task leases have task-call lifetimes, while owner-managed leases remain with the actual subsystem allocation. The graph never becomes a second physical accounting authority. An unknown peak is not zero memory and cannot authorize overlap.

`AnalysisRuntime` owns constant-size coalesced light/heavy cadence metadata (first/latest due epoch, missed count). It is not scientific/restart truth. Successful execution clears pending work; run-segment termination records and discards it. Checkpoint provenance records its pending summary, but restart starts a new cadence. No historical physical state is replayed or relabeled.

## 2026-09-30 M48-07 migration truth update

For an all-rank-proven homogeneous-DMO state with a uniform-rung-zero particle scheduler and no heterogeneous sidecar/softening/AMR exception, migration now has a distinct `compact_homogeneous_dmo` transaction representation. Decomposition `item_index` is the authoritative local-row selector; compact scalar metadata remains transaction-level truth and is communicator-consensus checked. The transaction streams fixed DMO wire records directly from canonical rows into a compact candidate, validates before publication, then commits by swap while preserving compact metadata and existing generation invalidation semantics. Unsupported states continue through `generic_full_physics` and the existing `ParticleMigrationRecord`/AMR compatibility transaction.

Exact global ownership validation remains authoritative but is now a hard-bounded derived diagnostic: deterministic hash-owner/radix buckets are refined until each exact sort/compare fits the 64 MiB/rank workspace contract. `ExactOwnershipPartitionReport` owns exact duplicate/missing/extra counts plus bounded diagnostic samples; compact-DMO callers borrow the canonical particle-ID lane rather than materializing an N-sized ID copy.

## P6 optional authority additions (source only, qualification pending)

Scheduler remains the sole bin/activation/active-set authority. Hierarchical
quantum and coarse endpoint belong to IntegratorState/CosmologicalTimeline;
transient StepContext factors describe one operator surface. Sources are all
drifted to each fine epoch, and SimulationState advances source generation once.
Uniform DMO drift epoch remains scalar. The optional particle bin byte lane is a
mirror, not an independent clock.

Gravity runtime owns derived PM/Tree component caches. Tree kicks require the
row's exact current source generation and row identity; PM endpoint kicks also
require its source generation and committed PM version. Previous total force is
only a MAC scale and cannot authorize a kick. Fine cached MAC scales may combine
fine Tree and coarse PM components. Migration invalidates row identity; mandatory
coarse bootstrap reconstructs split caches in normal and restart paths.

Tree solver retains original octant geometry/membership as a refit certificate;
coordinator owns fixed periodic unwrap/frame and source-layout checks. Refitted
bounds/moments publish a new build generation only when complete. P2 routing
geometry continues to refresh from canonical current sources independently of
ownership. Spatial work EMA and profiler rings are bounded derived diagnostics;
spatial history is cold after restart and is not force truth.
