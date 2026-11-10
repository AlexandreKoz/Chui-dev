# Repair state recap (post-repair audit snapshot)

## 2026-11-10 — post-patch-56 v10 final source closure (static, NOT RUN)

Restart v25 group identities repaired (`src/io/restart_checkpoint.cpp`);
explicit complete-coverage AMR domain contract
(`src/workflows/hydro_amr_runtime.cpp` `globalAmrHierarchyDomain`); S02 yield
coverage report (`SspYieldTable::channelCoverage`, `sources.data_assumptions`).
No new config key; restart remains v25 (writer corrected, no layout change);
AMR wire remains v5. Ledger:
`docs/repair/chui_post_patch56_v10_final_source_closure.md`.

## 2026-11-04 — design conformance A1/A2 source checkpoint

Recursive-driver working-set admission precedes allocation, and temporal-history
start/end operations publish complete candidates. Retirement is allocation-free
for record storage and reconciles freed nested capacities. This is local AMR
correctness work, not distributed multirate completion. A3–A5 remain open.
All authored/retained tests are **NOT RUN — USER-RESERVED**; no executable
validation was performed. See
[`repair/chui_conformance_closure_20261104.md`](repair/chui_conformance_closure_20261104.md)
for the H01–X01 map, ownership/API changes and continuation.

20261103 CHUI Part 1 (foundation + hydro suite): audit F1, F3, F4, F7 repaired
in source; PPM, CTU, HLLE selection, dual energy, thermodynamic refresh
contract, node-aware admission, owner telemetry and stage transactions added.
**IMPLEMENTED IN SOURCE — NOT COMPILED, NOT TESTED** (user-reserved).
F2, F5, F6 remain open and owned by Parts 2/3. See the
[Part 1 record](repair/chui_part1_foundation_hydro_20261103.md).

20261012 regression closure: short-range fixture dispatch, valid typed config
controls, field-specific negative checks, experiment deck invariants, PM owned
allocation expectations and the stale hierarchical admission diagnostic are
repaired in source. Always-active gravity assertions and an NDEBUG probe are
included. **FIXED IN SOURCE — NOT YET TESTED.** Production code is unchanged;
compilation, CPU/MPI tests, Ewald and benchmarks are deferred to the maintainer.
See the [defect table and manual sequence](repair/treepm_deep_hotpath_test_regression_closure_20261012.md).

20261012 Attempt 1 recovery: confirmed drift velocity access repaired, composed
runtime coverage added, bounded geometry and optional governed arithmetic/SIMD
paths implemented, production block/leaf controls and staged experiments wired,
launcher mode and source/toolchain provenance repaired. **IMPLEMENTED IN SOURCE
— NOT VALIDATED.** No compilation, scientific/MPI tests, benchmarks or long-run
convergence executed, by request. Exact dispositions and manual commands:
[handoff](repair/treepm_deep_hotpath_attempt1_closure_20261012.md).


20261010 P1–P5 corrective follow-up: Python 3.8-compatible complete-snapshot
reading, coincident residual/zero-vector repair, constrained verified-v23
configuration migration, endpoint PM guards, mandatory Ewald recovery profiles,
and test comparison/assertion/CTest repairs are implemented in source. All
compilation, tests, MPI, convergence and performance remain **NOT RUN — deferred
by user**. Legacy opaque-hash migration fails closed; historical external P0
tools remain absent; lookup and hierarchical modes remain opt-in/provisional.
See the [defect dispositions, exact manifest and commands](repair/treepm_p1p5_qualification_and_test_reliability_repair_20261010.md).

20261010 unified TreePM recovery: production traversal, stable analytic residual,
opt-in governed Gaussian lookup, staged MAC accounting, external split sweeps,
and hierarchical coarse-displacement/source-epoch hardening are implemented in
source. New scientific qualification tests are written. Compilation, tests and
performance qualification were explicitly deferred by the user; no new passing
evidence or capability certification is asserted. See the
[implementation and manual commands](repair/treepm_unified_performance_recovery_20261010.md).

P6 campaign label 20261006: source architecture implemented, qualification
pending. No build, tests, or numerical campaigns executed. Current details:
[implementation](repair/p6_preproduction_performance_architecture_implementation_20261006.md)
and [follow-up matrix](repair/p6_validation_campaign_plan_20261006.md).

> **Historical / non-authoritative.** This document records an earlier campaign state. Use [`CURRENT_STATUS.md`](../CURRENT_STATUS.md) for current repository truth.


## 2026-09-12 M2D memory owner-lifetime closure campaign

Targeted the remaining source-level memory ownership defects across distributed
reflux, active-level selection, regrid coexistence, source report/event paths,
and restart/readback governance.

### Changes applied

**Distributed reflux / MPI ownership (Phase B)**

- `executeBlockingAmrFluxRegisterPayloadExchange` (`distributed_mesh.hpp` +
  `distributed_memory.cpp`): added `core::MemoryGovernor*` parameter; inbound
  records buffer is now governed via `MemoryClass::kCommunication` reservation
  with checked-arithmetic element count.
- Replaced `std::unordered_set<uint64_t> inbound_keys` with sorted-vector
  adjacent-duplicate detection — eliminates a population-scale hash-table
  allocation in the MPI hot path.
- `applyCompletePendingFluxRegistersToSimulationState` (`amr_hydro_orchestrator.cpp`):
  removed `std::vector<uint64_t> applied_keys`; applies zeroing in-place and
  uses new `PendingFluxRegisterStore::eraseCompleted()` to clear consumed
  records without a second ownership domain.
- Caller in `hydro_amr_runtime.cpp` now passes `m_memory_governor` to the flux
  exchange function and wires `regrid_memory_governor` into `ProductionAmrHydroOptions`.

**Active-level and reflux materialization (Phase C)**

- `activeRowsForLevel` (`amr_hydro_orchestrator.cpp`): eliminated
  `std::unordered_set<uint32_t>` population-scale set; now iterates
  `requested_rows` directly with a level filter, pre-reserves returned vector
  capacity.

**Source report/event ownership (Phase E)**

- Added `StellarFeedbackModel::applyEventsCountersOnly` (`stellar_feedback.hpp` +
  `stellar_feedback.cpp`): production deposition path that accumulates counters
  and performs all side effects without materializing a per-star
  `star_reports` vector.
- `source_runtime.cpp`: production loop now calls `applyEventsCountersOnly`
  instead of discarding the return value of `applyEventsWithViews`.

**Restart/readback governance (Phase F)**

- Added `RestartReadPolicy` struct (`restart_checkpoint.hpp`) with
  `MemoryGovernor*` and `enforce_dimension_consistency` flag.
- `readRestartCheckpointHdf5` now accepts the policy (backward-compatible
  default); production caller in `output_restart_runtime.cpp` passes the
  governor.

**Forward-declaration fix**

- `distributed_mesh.hpp`: added `core::MemoryGovernor` forward declaration
  alongside the existing `ProfilerSession` forward declaration.

### Evidence

- `cmake --preset local-gcc-cpu-debug` configured successfully
- `cmake --build --preset build-local-gcc-cpu-debug -j 3` compiled 346/346
  targets
- `ctest --preset test-local-gcc-cpu-debug -j 3` passed 137/140 (3 pre-existing
  environment failures: missing `<span>` header in temp build, source package
  completeness, CI label test)
- `git diff --stat` confirms 11 files changed, 349 insertions, 36 deletions

### Files touched

| File | Change |
|------|--------|
| `include/cosmosim/core/simulation_state.hpp` | `eraseCompleted()` on `PendingFluxRegisterStore` |
| `include/cosmosim/parallel/distributed_mesh.hpp` | `MemoryGovernor*` param + forward decl |
| `include/cosmosim/physics/stellar_feedback.hpp` | `applyEventsCountersOnly` declaration |
| `include/cosmosim/io/restart_checkpoint.hpp` | `RestartReadPolicy` struct |
| `src/amr/amr_hydro_orchestrator.cpp` | `activeRowsForLevel` rewrite, flux register zeroing |
| `src/parallel/distributed_memory.cpp` | Governed MPI inbound buffers |
| `src/physics/stellar_feedback.cpp` | `applyEventsCountersOnly` implementation |
| `src/workflows/hydro_amr_runtime.cpp` | Pass governor to flux exchange |
| `src/workflows/source_runtime.cpp` | Use `applyEventsCountersOnly` |
| `src/workflows/output_restart_runtime.cpp` | Pass `RestartReadPolicy` |
| `src/io/restart_checkpoint.cpp` | Accept `RestartReadPolicy` |


## 2026-07-13 gravity/MPI final regression closure

- Determined that the former
  `saw_zero_local_work_with_remote_activity` failure came from a stale
  intermediate MPI test binary. The current fixture intentionally activates
  rung zero on every rank under the production `hierarchical_max_rung=0`
  fail-closed contract; rebuilding removed the obsolete assertion without a
  solver-behavior change.
- Closed three production collective defects found by the focused matrix and
  adversarial review: rank-local serial TreePM no longer enters distributed
  residual exchange; the distributed SFC actionable-migration vote is no
  longer hidden behind a rank-local short circuit; and MPI-enabled serial PM,
  TreePM, `MpiContext`, and PM-halo paths check initialization/finalization
  before communicator queries.
- PM and TreePM collective entry consensus now binds mesh shape and operator
  controls. PM halo preparation/allocation failures and protocol identity are
  coordinated before point-to-point traffic. Focused regressions cover
  divergent shapes/controls, rank-local failure, zero-width slabs, and
  np2/np3/np4 mixed-actionability SFC execution.
- Repaired the CUDA debug build for its configured `sm_52` target by using the
  standard integer-CAS double atomic below `sm_60` and resolving device symbols
  in the static gravity library. The CUDA 12.0.140 inventory passes 123/123,
  including the PM GPU smoke; CUDA remains CIC-only and outside the certified
  TSC force profile.
- Final local evidence: repository hygiene and `git diff --check` pass;
  CPU 122/122; HDF5 124/124; HDF5+FFTW 124/124; MPI+HDF5+FFTW 150/150 with
  real np1--np4 launchers; ASan+UBSan+LSan 122/122; CUDA 123/123; focused
  gravity/MPI 36/36; DMO workflow/rank comparator 5/5.
- Reproducibility impact: no snapshot/restart dataset, schema version, config
  key, normalized dump, provenance format, output naming, force split, or
  acceptance tolerance changed in these final repairs. Failure behavior is
  stricter and collective order is now uniform. The already documented
  TSC/deconvolution/cutoff default changes remain the numerical reproducibility
  change for the overall gravity hardening task.

## 2026-07-08 MPI AMR directed-exchange closure pass

- Replaced the live production AMR patch/cell import path with bounded
  rank-envelope peer discovery and directed candidate-peer patch/cell payload
  exchange. Remote AMR data are transient read-only ghosts rather than globally
  replicated authoritative state.
- Replaced AMR reflux completion with owner-routed peer exchange keyed by
  `owner_rank`, so flux-register corrections are sent only to authoritative
  coarse-cell owners.
- Relabeled the existing serial AMR boundary regression as serial-only evidence
  and added a real MPI-launched AMR directed-exchange regression target for
  dependency-enabled environments.
- Reproducibility impact: no HDF5 snapshot/restart dataset names, normalized
  config keys, output naming, solver equations, or rank-count-changing restart
  semantics changed. Distributed AMR remote patch state remains transient
  stage scratch rebuilt from owner payloads.

## 2026-06-30 MPI.7 PM/TreePM transport scale guards

- Preserved existing slab-owned periodic PM routed deposition/interpolation and
  TreePM short-range request identity checks; added profile counters for routed
  PM records, participating peers, and slab-halo force cache hits so scale runs
  can distinguish sparse payload traffic from collective completion plumbing.
- Added explicit bounded policies for the live unbounded paths: isolated/open PM
  root gather/scatter now estimates root workspace and fails before MPI/data
  allocation when `parallel.isolated_pm_root_workspace_limit_bytes` is exceeded;
  focused zoom high-resolution source all-gather is bounded by
  `parallel.zoom_high_res_allgather_limit_bytes`.
- Reproducibility impact: no solver equation, HDF5 schema, restart payload, or
  output naming change. Two normalized config keys now document the supported
  envelope for gathered PM/zoom compatibility paths.

## 2026-06-30 MPI.6 distributed SFC runtime rebalance

- Replaced routine runtime rebalance planning with distributed compact SFC cuts: ranks sort local
  authoritative entities, exchange bounded cut samples, all-reduce exact balance/movement counters,
  and migrate only particles/AMR patches crossing accepted cuts.
- Kept the exact global `DecompositionItem` all-gather planner as debug/small-test evidence, gated
  production workflow exact ownership validation behind
  `parallel.decomposition_debug_exact_ownership_audit = false` by default.
- Reproducibility impact: no HDF5/restart schema, output naming, rank-count-changing restart, or
  solver-numerics change. Runtime decomposition ownership is stricter: accepted SFC cuts are the
  owner authority, while compact summaries and exact audit payloads are read-only diagnostics.

## 2026-06-30 MPI.5 distributed AMR patch boundaries and reflux

- Added a bounded synchronized distributed AMR hydro path: owner ranks export explicit patch geometry and
  gas-cell conserved state, consumers build transient read-only remote patch ghosts only for adjacent AMR
  faces, and coarse-fine flux-register entries whose coarse cell is remote are sent back to the authoritative
  owner rank for reflux.
- Kept local AMR subcycling and temporal ghost histories local-only. Distributed AMR calls that request
  subcycling or temporal coarse-to-fine interpolation fail explicitly rather than implying a scheduler-owned
  MPI Berger-Colella timeline.
- Reproducibility impact: no config keys, normalized config dumps, HDF5 dataset names, schema versions, or
  rank-count-changing restart semantics changed. Existing restart-authoritative pending flux-register state
  remains owner-local and is preserved by the v19 restart groups.

## 2026-06-30 MPI.4 distributed hydro ghost boundary path

- Added a live fixed-grid ReferenceWorkflow hydro MPI boundary path: ranks advertise compact boundary-cell geometry, request adjacent remote gas-cell conserved state by stable `gas_cell_id`, bind it to explicit `kImportedMpi` ghost rows before reconstruction, and route deterministic face-owner conservative ghost deltas back to authoritative owners.
- Added packet validation and a two-rank MPI+HDF5 workflow regression source covering remote ghost diagnostics, global mass conservation, and restart-vs-uninterrupted continuation for the exercised live interface. Local execution remains dependency-bound if MPI/HDF5 tooling is unavailable.
- Reproducibility impact: no config keys, normalized config dumps, HDF5 schema datasets, output naming, or rank-count-changing restart semantics changed. Remote hydro ghost state is transient stage scratch keyed by epoch/gas-cell identity generation and is rediscovered after migration or restart.

## 2026-06-30 MPI.3 H2 gas-cell identity migration

- Replaced reference-workflow production compaction/rebalance gas handling with explicit gas-cell/AMR patch payloads keyed by stable `gas_cell_id` and `patch_id`; particle migration now has an explicit `preserve_gas_cell_state` mode so optional parent-particle removal cannot delete hydro truth.
- Runtime AMR patch migration carries gas-cell scheduler identity records and rebuilds the gas scheduler after a validated commit. Initial decomposition removes stale gas patches by patch ownership, not by one-cell-per-gas-particle reconstruction.
- Added a two-rank MPI+HDF5 workflow regression source for parentless cells, two cells sharing one parent, restart readback, and restart resume. Local execution remains environment-blocked here because `cmake --preset mpi-hdf5-fftw-debug` cannot find `MPI_CXX`.
- Reproducibility impact: no config keys, restart HDF5 datasets, schema version, normalized config dumps, or output naming changed. The runtime migration transaction is stricter and fails before commit for wrong-owner or incomplete scheduler/gas payloads.

## 2026-06-07 Agent interface update

- Promoted `AGENTS.md` from a repair-only guardrail into the repository-wide CHUI/CosmoSim agent contract covering audit, repair, feature implementation, and clean packaging modes.
- Added `docs/architecture/agent_task_interface.md` and `.github/copilot-instructions.md` so Codex/Copilot-style agents share one task router, evidence floor, and clean-artifact policy.
- Reproducibility impact: documentation/interface only; solver numerics, config normalization, output naming, snapshot/restart schemas, and runtime state ownership are unchanged.


## 2026-05-25 Stage 6 memory accounting infrastructure by subsystem (S6-07)

- Added central `core::MemoryReport` accounting API with subsystem categories, persistent/transient split, capacity-based owned-byte semantics, and explicit uncertainty notes for external allocations.
- Wired `ProfilerSession` optional memory report attachment so diagnostics can carry machine-readable memory budgets without making restart payloads depend on allocator internals.
- Reproducibility impact: observational diagnostics only; solver numerics, config normalization, output naming, and restart truth are unchanged.

## 2026-05-21 Stage 5 provenance schema/version audit hardening (S5-04)

- Bumped provenance payload contract to `provenance_v5` with explicit config schema identity and auditable raw/normalized/derived runtime payload fields.
- Wired snapshot provenance metadata to carry `config_schema_name`, `config_schema_version`, `normalized_config_hash_hex`, and raw/normalized/derived config/state payloads.
- Reproducibility impact: strengthens replay auditability without changing solver numerics; normalized hash remains derived from normalized config only.

## 2026-05-21 Stage 5 config-contract audit map (S5-00)

- Added `docs/repair/stage5_config_contract_audit_map.md` as a compact file-by-file implementation map for Stage 5 config contract hardening.
- Enumerated raw input boundary, normalized typed contract, derived runtime projection, runtime entrypoints, alias-boundary constraints, and riskiest mixed-lane call sites to sequence follow-up prompts without implementation churn.
- Reproducibility impact: documentation-only; no solver numerics, config keys, normalized-config emission, provenance hash logic, snapshot/restart schema, output naming, or rank coordination behavior changed.



## 2026-05-19 Stage-bound orchestrator dispatch repair

- Replaced authoritative callback broadcast with per-stage handler buckets: handlers now declare typed `IntegrationStage` sets at registration, and `StepOrchestrator` dispatches only the active stage bucket in stable registration order.
- Converted production gravity, drift, hydro, source, tracer, and diagnostics handlers away from silent off-stage self-filtering; impossible direct off-stage calls now fail as contract violations while legitimate same-stage no-op conditions remain unchanged.
- Reproducibility impact: deterministic handler order is preserved by registration order within each stage, solver numerics and KDK stage order are unchanged, and no config keys, normalized config dumps, output naming, rank coordination, or snapshot/restart schemas changed. The dispatch path is more auditable and prevents irrelevant-stage delivery by construction.


## 2026-05-11 Stage 2 authority documentation synchronization

- Synchronized implementation-facing docs, runtime truth map, restart/schema language, ADR entries, repair-state notes, release scope, and release notes around the Stage 2 owner model: scheduler objects own timestep bins, active-set construction, and PM cadence legality; public `time_bin` lanes remain mirrors.
- Added a docs-scaffold regression so future edits must keep the scheduler authority contract, mirror policy, candidate criteria flow, active-set construction flow, invariant framework, and forbidden Phase 3 production-claim language present in docs/release artifacts.
- Reproducibility impact: documentation/test-scaffold only. No solver numerics, restart schema fields/version, HDF5 dataset layout, normalized config dumps, provenance format, output naming, or deterministic scheduler behavior changed. Restart behavior remains the documented v6 compatibility tightening that rejects stale mirrors and rebuilds valid mirrors from scheduler authority.
- Command evidence in this cycle:
  - `cmake --preset cpu-only-debug` → pass.
  - `cmake --build --preset build-cpu-debug` → pass.
  - `ctest --preset test-stage1-runtime-truth-cpu-debug --output-on-failure` → pass (17/17).
  - `ctest --preset test-cpu-debug --output-on-failure` → pass (77/77).
  - `cmake --preset hdf5-debug`, `cmake --build --preset build-hdf5-debug`, `ctest --preset test-hdf5-debug --output-on-failure` → pass (78/78 HDF5 tests).
  - `cmake --preset pm-hdf5-fftw-debug`, `cmake --build --preset build-pm-hdf5-fftw-debug`, `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` → pass (78/78 PM+HDF5+FFTW tests).
  - `cmake --preset mpi-hdf5-fftw-debug` → environment-blocked by missing `fftw3_mpi`; no MPI+HDF5+FFTW closure is claimed.

## 0) HDF5 preset smoke/restart gate repair (2026-05-10 UTC)

- Fixed the HDF5 runtime smoke fixture to provide the complete axis-aware TreePM grid triplet after the config contract started rejecting partial `treepm_pm_grid_n{xyz}` input.
- Fixed the Phase 2 restart-continuation validation fixture to set `snapshot_interval_steps = 1`, ensuring the three-step validation run actually reaches the restart write/read gate instead of inheriting the default 64-step cadence.
- Reproducibility impact: this is test/config wiring only; no solver numerics, restart schema, normalized config generation, or scheduler authority semantics changed.

## 0) Stage 2 timestep duplicate-lane reclassification and restart guards (2026-05-10 UTC)

- Enforced scheduler-owned restart truth for timestep bins: restart hash/write/read paths now validate particle `time_bin` mirrors against scheduler `bin_index`, reject stale mirrors, validate scheduler schema arrays, and rebuild mirrors from scheduler state on read.
- Reclassified migration/transfer `time_bin` fields as derived mirrors with scheduler remap/rebuild contracts, and expanded Stage 2 audit/output-schema docs with the full bin/cadence/kick-opportunity/field-version lane inventory.
- Added tests for stale restart mirror rejection, restart mirror rebuild, migration/reorder mirror consistency by stable IDs, scheduler import schema validation, and diagnostics not mutating scheduler truth.
- Reproducibility impact: exact restart continuation is strengthened without a schema bump; solver numerics, output names, normalized config dumps, and PM cadence semantics are unchanged.

## 0) Active-view lifetime, generation, and workspace reuse hardening (2026-05-07 UTC)

Validation commands for this repair:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state test_unit_hot_cold_sidecar_layout test_integration_reorder_compaction_sidecars test_unit_tree_pm_split_kernel
ctest --preset test-cpu-debug --output-on-failure -R "unit_simulation_state|integration_reorder_compaction_sidecars|unit_hot_cold_sidecar_layout"
ctest --preset test-cpu-debug --output-on-failure -R "active|stale|kernel|hot_cold|workspace"
```

Observed:

- Added generation stamps to read-only particle and cell active views, matching the existing mutable kernel-view scatter generation guards.
- Kept active/kernel view storage as `TransientStepWorkspace` scratch mirrors only, with comments documenting non-ownership, allowed mutable hot lanes, and forbidden metadata lanes.
- Extended regression coverage for stale particle scatter after migration, stale hydro scatter after gas-cell identity rebuild, and workspace clear/reuse preserving capacity without exposing stale rows.

Migration note:

- Callers should continue to construct active views through `buildParticleActiveView` and `buildCellActiveView`; aggregate/manual construction must initialize the new source-generation field to preserve stale-view auditability. No config keys, snapshot/restart payloads, or solver inputs migrate.

Reproducibility impact:

- Deterministic behavior is strengthened: stale compact views now carry auditable source generations across read-only and mutable active paths, and mutable scatter continues to fail before row-index writes after structural transforms. No solver numerics, restart/HDF5 schema, output naming, normalized-config generation, rank coordination, or persistent ownership semantics changed.


## 0) Transform fuzz/invariant harness hardening (2026-05-07 UTC)

Validation commands for this repair:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_hot_cold_sidecar_layout test_integration_reorder_compaction_sidecars test_integration_species_migration_invariants test_integration_transform_fuzz_invariants
ctest --preset test-cpu-debug --output-on-failure -R "fuzz|reorder|migration|sidecar|stale"
ctest --preset test-cpu-debug --output-on-failure -R "fuzz|reorder|migration|sidecar|stale"
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug -j4 --target test_integration_transform_fuzz_invariants
ctest --preset test-hdf5-debug --output-on-failure -R "transform_fuzz"
```

Observed:

- Added a fixed-seed transform fuzz integration harness that composes random reorder, active-view materialization, deterministic compaction, sparse sidecar migration commit, gas-state rebuild by particle ID, and HDF5 restart roundtrip when HDF5 is enabled.
- The harness uses particle-ID keyed oracles for common particle lanes, nonmonotonic IDs, SFC keys, time bins, sparse gravity-softening overrides, and all current species sidecar families (stars, black holes, tracers).
- Added negative stale active-view scatter checks after reorder and compaction so stale row-index views fail loudly instead of silently corrupting transformed rows.
- Registered a short default integration run plus a longer `validation`-labeled run of the same deterministic harness.

Reproducibility impact:

- Deterministic behavior is strengthened by fixed seeds, deterministic failure context (`seed`/`trial`/stage), particle-ID keyed comparisons after every transform, and two identical CPU CTest runs. No solver numerics, restart schema version, HDF5 dataset names, output naming, normalized-config generation, scheduling semantics, or rank coordination behavior changed.


## 0) Particle-bound gas-cell contract hardening (2026-05-07 UTC)

Validation commands for this repair:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_gas_cell_identity_invariants test_integration_species_migration_invariants
ctest --preset test-cpu-debug --output-on-failure -R "gas|hydro|migration|reorder|restart"
cmake --build --preset build-hdf5-debug -j4
ctest --preset test-hdf5-debug --output-on-failure -R "restart|roundtrip|gas"
```

- Replaced hidden gas count/cell row assumptions in core active hydro views and reference workflow gas consumers with `requireParticleBoundGasCellContract(...)`, plus named row/ID helpers for parent-particle ID lookup and gas-particle row lookup.
- Strengthened gas identity tests to cover helper seams, hydro active-view contract rejection after invalid cell resize, forbidden gas-relative reorder, and migration/compaction rebuild of hydro fields by gas particle ID.
- Reproducibility impact: no solver numerics or HDF5/restart schema changed. The repair fails earlier when gas-cell identity lanes are stale, preserving deterministic gas field attachment by particle ID across reorder/migration/restart boundaries.


## 0) Per-particle gravity softening transform/restart guardrail hardening (2026-05-07 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_integration_softening_ownership_invariants test_integration_reorder_compaction_sidecars test_integration_species_migration_invariants
ctest --preset test-cpu-debug --output-on-failure -R "softening|reorder|migration|sidecar"
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug -j4
ctest --preset test-hdf5-debug --output-on-failure -R "restart|snapshot|softening"
```

Observed:

- Added deterministic randomized reorder stress coverage for sparse per-particle gravity softening overrides, using particle ID as the oracle after every permutation.
- Extended resize/compaction coverage to prove newly grown rows do not become overrides and retained particle IDs preserve `{has_override, epsilon}` through shrink/regrow cycles.
- Extended Tree softening view coverage so mismatched source/target override masks fail loudly instead of silently truncating override authority.
- Fixed the HDF5 softening ownership restart test payload metadata so the restart reader validates the normalized-config hash lane before comparing softening value/mask roundtrip state.
- Documented that `has_gravity_softening_override` is the authoritative restart lane; value-only softening lanes are diagnostics/default mirrors and do not create overrides.

Reproducibility impact:

- Deterministic behavior is preserved. This pass adds regression coverage and documentation for existing softening ownership semantics; it does not alter solver numerics, restart dataset names, output naming, normalized config generation, scheduling, or rank coordination.

## 0) Stage 1 runtime transform audit freeze (2026-05-06 UTC)

Commands:

```bash
git status --short
rg -n "SimulationState|sidecar|ownership|runtime truth|reorder|migration|gas|restart|snapshot|softening" tests include src docs/architecture docs/repair -g '*.{cpp,hpp,md,txt,cmake}'
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure -R "simulation_state|ownership|runtime|truth"
```

Observed:

- Added `docs/repair/stage1_runtime_transform_audit.md` as a documentation-only mutation-map proof for particle, gas-cell, species-sidecar, restart, snapshot, migration, and reorder paths.
- Classified Stage 1 P0 targets 1.1, 1.2, and 1.3 from the prompt scope because no separate in-repository Stage 1 target list was found by `rg`.
- No code, build registration, solver behavior, SoA layout, restart schema, snapshot schema, or HDF5 dataset names were changed.
- The required targeted `ctest` regex initially selected integration executables that had not been built by the narrower `test_unit_simulation_state` target; after a full `cmake --build --preset build-cpu-debug -j4`, the same regex passed 3/3.

Reproducibility impact:

- Deterministic behavior unchanged. This pass records audit evidence and current mutation authority only; no runtime behavior, output naming, normalized config, provenance, scheduler, restart, snapshot, or physics numerics changed.

## 0) Stage 0 consolidation gate and runtime-truth freeze (P0-10, 2026-04-26 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "runtime|truth|timestep|scheduler|active|particle|sidecar|species|gas|hydro|softening|config|restart|checkpoint|provenance"
ctest --preset test-cpu-debug --output-on-failure

cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure

cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```

Observed:

- Stage 0 targeted runtime-truth suite on `cpu-only-debug` passed (21/21) across timestep/bin, active-set, particle/sidecar, species, gas, softening, config/provenance, and restart/checkpoint families.
- Full CPU/HDF5 suite runs are currently not fully green in this environment:
  - `unit_snapshot_hdf5_schema` mismatch (`gadget_arepo_v2` expected by test while runtime schema references newer naming),
  - long-running `integration_tree_pm_coupling_periodic` path limits timely full-cycle completion in this container.
- PM+HDF5+FFTW preset run completed with explicit failures (`8/75` failed): `integration_softening_ownership_invariants`, `unit_pm_solver`, `unit_snapshot_hdf5_schema`, `integration_tree_pm_coupling_periodic`, `integration_docs_scaffold`, `integration_release_readiness_artifacts`, `integration_runtime_app_smoke`, `validation_phase2_mpi_gravity_single_rank`.
- Added Stage 0 grouped ctest preset: `test-stage0-runtime-truth-cpu-debug`.
- Added closure artifact `docs/repair/stage0_runtime_truth_freeze.md` and recorded gate ADR entry (`ADR-INFRA-STAGE0-GATE-014`).
- Updated `docs/repair_open_issues.md` with precise Stage 0 closure blockers and retest commands.

Reproducibility impact:

- Deterministic behavior unchanged. This pass records gate evidence, test-group metadata, and blocker docs only; no new physics or solver-model features were introduced.

## 0) Restart/reload runtime-truth round-trip invariant suite (2026-04-26 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_restart_checkpoint_roundtrip
ctest --preset test-cpu-debug --output-on-failure -R "restart|checkpoint|roundtrip|provenance"
```

Observed:

- Extended `tests/integration/test_restart_checkpoint_roundtrip.cpp` to harden runtime-truth continuation checks for:
  - canonical particle ID/order survival, species ledger/index reconstruction, gas identity mapping, and sidecar roundtrip,
  - scheduler persistent-state restoration with active-set reconstruction equivalence at resume tick,
  - per-particle softening override continuation and precedence behavior after reload,
  - normalized config/provenance hash continuity across restart read.
- Added explicit missing/legacy field behavior checks:
  - missing required scheduler lane (`pending_bin_index`) fails with a clear field-specific error,
  - missing legacy optional per-particle softening dataset follows documented compatibility behavior (`no overrides`).
- Updated restart/ownership documentation (`docs/restart_checkpointing.md`, `docs/architecture/adr_runtime_truth_ownership.md`) with authoritative serialized fields, derived reconstruction rules, active-set rebuild expectations, and compatibility-path constraints.
- Updated issue ledger in `docs/repair_open_issues.md` to mark P0-09 closed with command-backed evidence.

Reproducibility impact:

- Deterministic behavior is preserved; this pass adds restart-continuation invariant tests and documentation, without changing solver numerics or adding new physics.

## 0) Active-set ownership + cache-invalidation invariant floor (2026-04-26 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug --target test_unit_time_integration test_unit_simulation_state
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_time_integration|unit_simulation_state"
```

Observed:

- Added active-set authority and divergence-risk coverage in `tests/unit/test_time_integration.cpp`:
  - `testActiveSetAuthority` validates known scheduler-driven extraction from controlled bin assignments plus authorized bin mutation refresh behavior.
  - `testActiveSetNoCompetingBuilders` demonstrates that a mirror-based competing builder can diverge after non-owner mirror mutation and is caught by `debugAssertTimeBinMirrorAuthorityInvariant(...)`.
- Added mutable compact-view invalidation + solver-derivation checks in `tests/unit/test_simulation_state.cpp`:
  - stale gravity compact views fail loudly after particle reorder generation change,
  - stale hydro compact views fail loudly after cell resize generation change,
  - gravity/hydro solver-local views are verified as derived from the same authoritative active indices used by read-only active views.
- Updated `docs/architecture/adr_runtime_truth_ownership.md` with explicit active-set owner/layering rules, forbidden module-local fallback builders, generation-counter scatter rules, and lifetime/invalidation events.
- Updated issue ledger entry in `docs/repair_open_issues.md` to mark P0-08 closed with command-backed test/document evidence.

Reproducibility impact:

- Deterministic behavior is preserved; this pass adds/strengthens ownership tests and explicit invalidation policy without changing solver numerics or physics models.

## 0) Config/runtime ownership and provenance-restart consistency hardening (2026-04-26 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_config_parser test_unit_units_cosmology_provenance test_unit_restart_checkpoint_schema
ctest --preset test-cpu-debug --output-on-failure -R "config|provenance|runtime|restart"
ctest --preset test-cpu-debug --output-on-failure
```

Observed:

- Added targeted config/runtime ownership coverage in `tests/unit/test_config_parser.cpp`:
  - deterministic normalized-text + hash stability from equivalent raw inputs,
  - legacy alias-to-canonical normalization checks for runtime-owned numerics lanes,
  - normalized text reparse checks to ensure one frozen/typed authority lane.
- Added derived-constant consistency coverage in `tests/unit/test_units_cosmology_provenance.cpp`:
  - normalized typed config -> `UnitSystem` mapping consistency,
  - normalized typed cosmology -> `LambdaCdmBackground` H0 derivation consistency.
- Added restart metadata integrity rejection coverage in `tests/unit/test_restart_checkpoint_schema.cpp`:
  - explicit rejection of normalized-text/hash mismatches even when caller-provided hash/provenance strings are non-empty.
- Updated ownership/config docs (`docs/architecture/adr_runtime_truth_ownership.md`, `docs/configuration.md`) with an explicit raw/normalized/derived/mutable/provenance classification map and anti-ambiguity notes for `time_begin_code`, `time_end_code`, `current_scale_factor`, and diagnostic redshift.

Reproducibility impact:

- Deterministic behavior is preserved and tightened: this pass adds validation/tests/docs only and does not change solver numerics or physical-model equations.

## 0) Softening ownership/priority + override-preservation invariant floor (2026-04-25 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_softening_ownership_invariants
ctest --preset test-cpu-debug --output-on-failure -R "softening|gravity|sidecar"
ctest --preset test-cpu-debug --output-on-failure
```

Observed:

- Added targeted invariant coverage in `tests/integration/test_softening_ownership_invariants.cpp` for:
  - explicit source-priority checks (`per-particle override > species default > global default`) in tree softening resolution helpers,
  - preservation of per-particle override values across resize grow/shrink, reorder/permutation, and species migration commit,
  - restart checkpoint roundtrip retention of per-particle softening override values (HDF5-enabled path),
  - active-set softening extraction mirror checks that ensure gathered lanes reflect sidecar truth,
  - diagnostics observer checks proving run-health softening diagnostics do not mutate runtime softening state.
- Registered `integration_softening_ownership_invariants` in CMake and updated ownership ADR softening rules to call out resize preservation, active-set mirror status, restart serialization, and diagnostic/provenance non-authority.

Reproducibility impact:

- Deterministic behavior is preserved; this pass adds invariant tests and ownership documentation without changing gravity force physics or introducing new softening models.


## 0) Gas-cell identity and hydro-state ownership invariant tests (2026-04-25 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_gas_cell_identity_invariants test_integration_restart_checkpoint_roundtrip
ctest --preset test-cpu-debug --output-on-failure -R "gas|hydro|cell"
ctest --preset test-cpu-debug --output-on-failure
```

Observed:

- Added targeted invariants in `tests/unit/test_gas_cell_identity_invariants.cpp` for:
  - explicit local gas identity construction checks (stable gas-particle ID anchor, hydro sidecar ownership, and active hydro extraction safety),
  - reorder and resize drift detection with loud failures when the temporary local 1:1 gas-particle/gas-cell contract is violated,
  - stale hydro active-view scatter invalidation after cell index-space resize,
  - explicit separation between persistent hydro sidecar fields and active-kernel scratch views.
- Hardened runtime guardrails:
  - `core::debugAssertGasCellIdentityContract(...)` added for cheap invariant checks,
  - gas migration/compaction rebuild helpers now assert the temporary gas contract before reconstructing local gas ownership maps,
  - `reorderParticles(...)` now fails loudly when gas-relative reorder would silently drift gas-cell ownership without an ID-based rebuild.
- Extended restart roundtrip coverage in `tests/integration/test_restart_checkpoint_roundtrip.cpp` to assert gas density-by-particle-ID mapping survives restart I/O.
- Updated runtime ownership ADR to explicitly document temporary contract status, stable fields, scratch/derived fields, and forbidden gas identity assumptions.

Reproducibility impact:

- Deterministic behavior is preserved and made louder on contract violations; no hydro solver numerics or physics models were changed.


## 0) Sidecar reorder full-payload semantics hardening (2026-05-07 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_integration_reorder_compaction_sidecars test_unit_hot_cold_sidecar_layout
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure -R "reorder|sidecar|hot_cold|species"
rg "particle_index" src include tests
```

Observed:

- Refactored `reorderParticles(...)` species sidecar synchronization through typed row-lane visitors so star, black-hole, and tracer scalar lanes plus stellar channel arrays move as complete rows under `SidecarSyncMode::kMoveWithParent`.
- Preserved default `SidecarSyncMode::kUseParentIndirection` semantics: row order is stable and only `particle_index` is remapped to the post-reorder parent index.
- Added a post-reorder defensive guard (`debugAssertSpeciesSidecarOwnershipInvariants`) enforcing exactly one sidecar row for each eligible star/BH/tracer parent and no sidecar row attached to ineligible species.
- Expanded integration coverage with adversarial reverse permutation, explicit move-vs-indirection comparison, and randomized reorder-mode fuzzing over `kByTimeBin`, `kBySfcKey`, and `kBySpecies` while checking full payload identity by particle ID.
- Audited `particle_index` references with ripgrep; reorder entry points remain centralized through `reorderParticles(...)`, while physics modules only read/write sidecar parent indices during creation/evolution paths.

Reproducibility impact:

- Deterministic reorder behavior is preserved and made more auditable. No solver numerics, restart schema, HDF5 payload schema, or MPI/GPU ownership contracts were changed.


## 0) Particle ordering/resize/reorder sidecar invariants (2026-04-25 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_reorder_compaction_sidecars
ctest --preset test-cpu-debug --output-on-failure -R "reorder|sidecar|resize|particle"
ctest --preset test-cpu-debug --output-on-failure
```

Observed:

- Added deterministic integration invariants in `tests/integration/test_reorder_compaction_sidecars.cpp` for:
  - particle identity pairing across grow/shrink operations,
  - canonical reorder + adversarial permutation consistency,
  - species labels, timestep bins, and per-particle softening override continuity by particle ID,
  - full payload movement for star/black-hole/tracer sidecar rows under `SidecarSyncMode::kMoveWithParent`,
  - stale mutable compact-view invalidation via index-space generation checks.
- Repaired `ParticleSidecar::resize` so optional `gravity_softening_comoving` values are preserved across structural resize when the override sidecar is populated.
- Repaired `reorderParticles` `kMoveWithParent` handling to move complete sidecar rows (not index-only remap), then remap parent indices to post-reorder positions.
- Added generation counters on `SimulationState` index spaces and enforced generation mismatch failures in particle/cell kernel scatter paths.

Reproducibility impact:

- Deterministic behavior is preserved; this pass hardens ownership/index-safety checks and sidecar movement invariants without changing solver numerics.


## 0) Species migration + canonical grouping invariants (2026-04-25 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_species_migration_invariants test_unit_simulation_state
ctest --preset test-cpu-debug --output-on-failure -R "species|migration|sidecar"
ctest --preset test-cpu-debug --output-on-failure
```

Observed:

- Added deterministic migration-focused integration coverage in `tests/integration/test_species_migration_invariants.cpp` for:
  - single-particle species migration with stable particle identity and canonical species regrouping checks,
  - batch migration with sidecar-family transitions (star/BH/tracer attach-detach invariants),
  - migration followed by reorder/resize and scheduler active-set extraction checks,
  - softening priority continuity (`per-particle override > species default > global default`) after migration/reorder.
- Repaired `SimulationState::{packParticleMigrationRecords,commitParticleMigration}` to preserve optional per-particle softening override sidecar values across migration commits, including inbound records.
- Updated runtime ownership docs to state allowed migration path ownership, sidecar attach/detach expectations, and species-agnostic scheduler active-eligibility behavior.

Reproducibility impact:

- Deterministic behavior is preserved; this pass hardens migration-state invariants and sidecar continuity without adding new physics or changing solver numerics.


## 0) Timestep/bin authority invariant test floor (2026-04-25 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_time_integration
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_time_integration|integration_hierarchical_time_bins|integration_restart_checkpoint_roundtrip"
```

Observed:

- Added deterministic timestep/bin authority invariants in `tests/unit/test_time_integration.cpp` for:
  - scheduler-owned assignment/reassignment,
  - unauthorized mirror mutation detection,
  - reorder identity survival (bin follows particle ID),
  - scheduler persistent-state roundtrip equivalence (restart lane contract).
- Added explicit helper seam in `core/time_integration`:
  - `timeBinMirrorsMatchScheduler(...)`
  - `debugAssertTimeBinMirrorAuthorityInvariant(...)`
  so non-owner mutation of mirror lanes is detectable in debug/test paths.
- Updated ownership documentation to record these enforcement helpers and to reduce the timestep missing-test gap in the runtime truth map.

Reproducibility impact:

- Deterministic behavior is preserved. No timestep algorithm or solver numerics changed; this pass adds invariant checks/tests and documentation only.


## 0) Runtime ownership ADR freeze (2026-04-25 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

Observed:

- Added formal ownership ADR `docs/architecture/adr_runtime_truth_ownership.md` that assigns one authoritative owner per runtime domain and codifies mirror/cache/view refresh+invalidation, mutation authority, restart/reload rules, reorder/resize/migration rules, and forbidden duplicate-authority patterns.
- Updated architecture decision log with `ADR-INFRA-OWNERSHIP-013` linking the new runtime ownership policy to the Stage 0 truth map baseline.

Reproducibility impact:

- Documentation/policy-only change; no runtime/state schema or solver behavior change.


## 0) Stage 0 runtime-truth ownership map freeze (2026-04-25 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

Observed:

- Added a code-verified runtime ownership map in `docs/architecture/runtime_truth_map.md` covering timestep bins, active-set construction, particle/species ordering, sidecars, gas identity, softening, config-derived runtime values, and restart/checkpoint truth.
- Recorded explicit duplication/mirror lanes and mutation-authority ambiguity for scheduler/state time bins and gas identity mapping.
- Documented Stage 0 follow-up repair tickets without broad runtime behavior changes.
- `cpu-only-debug` configure/build succeeded; test preset run exposed an existing blocker: `test_unit_snapshot_hdf5_schema` expects schema `gadget_arepo_v2` while runtime schema reports `gadget_arepo_v4`.

Reproducibility impact:

- No runtime behavior or schema payload changed in this pass; this is an audit/documentation freeze only.

## 0) Gravity invalid-state diagnostics/escalation hardening (2026-04-20 UTC)

Observed:

- Added explicit gravity health checks in the reference workflow for PM/force finiteness, sync-state legality, and zoom/decomposition sanity.
- Added per-kick gravity runtime summaries (`gravity.health_summary`) plus targeted warning/fatal events (`gravity.health_check`).
- Fatal gravity-state violations now throw immediately and are not hidden behind generic diagnostics bundles.
- Added additive run-health gravity counters in diagnostics bundles (`gravity_softening_sidecar_size_ok`, `non_finite_gravity_softening`, `non_positive_particle_mass`).
- Extended restart roundtrip verification and distributed restart compatibility reporting to include gravity cadence/field-state consistency checks.

Reproducibility impact:

- Cheap gravity integrity checks are always-on and deterministic.
- Heavy reference checks remain opt-in behind `analysis.diagnostics_execution_policy = all_including_provisional`; default runtime behavior is unchanged.

## 0) Zoom long-range gravity contract hardening (2026-04-20 UTC)

Observed:

- Added typed zoom metadata/config validation plus authoritative file-driven high-resolution membership and focused PM correction controls.
- Implemented explicit zoom force decomposition in TreePM coupling:
  `PM_global_coarse(all) + [PM_focused(highres)-PM_coarse(highres)] + tree_short_residual(all)`.
- Added low-resolution contamination counters and force-component norms to runtime diagnostics/provenance pathways.

Reproducibility impact:

- Zoom membership and correction controls are now explicit typed config/provenance lanes; restart/snapshot continuation remains deterministic under the existing cadence policy.

## 0) Tree gravity multipole/MAC maturation repair (2026-04-20 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target \
  test_unit_tree_gravity \
  test_integration_tree_gravity_vs_direct \
  test_integration_tree_pm_coupling_periodic \
  bench_tree_gravity
ctest --test-dir build/cpu-only-debug --output-on-failure -R \
  "unit_tree_gravity|integration_tree_gravity_vs_direct"
./build/cpu-only-debug/bench_tree_gravity
```

Observed:

- Tree gravity now supports explicit multipole order selection with production path set to quadrupole (`l=2`) plus monopole fallback.
- MAC selection is now explicit (`geometric` and COM-distance-aware), with the COM-distance-aware criterion set as default and validated against direct-force references on quasi-uniform and clustered distributions.
- Traversal order was updated to push near children first from compact SoA sidecars to improve locality in hot loops without adding node-struct bloat.
- TreePM residual local and remote-source evaluation paths now use the same upgraded MAC + multipole acceptance/evaluation semantics.

Reproducibility impact:

- Self-exclusion, softening policy application, TreePM split factor multiplication, and cutoff pruning semantics are preserved; changes are in truncation order and acceptance policy only.

## 0) Axis-aware TreePM geometry contract repair (2026-04-20 UTC)

## 0) Hierarchical KDK + TreePM explicit PM synchronization contract repair (2026-04-20 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target \
  test_integration_reference_workflow_distributed_treepm_mpi \
  test_integration_time_integration_loop \
  test_integration_hierarchical_time_bins \
  test_integration_hierarchical_timestep_regression \
  test_unit_parallel_distributed_memory
ctest --test-dir build/cpu-only-debug --output-on-failure -R \
  "integration_hierarchical_time_bins|integration_hierarchical_timestep_regression|integration_time_integration_loop|unit_parallel_distributed_memory"
```

Observed:

- The PM long-range kick is now documented and recorded as an explicit kick-surface operator at `gravity_kick_pre` and `gravity_kick_post`.
- Rank-consensus checks were hardened to include PM field version and last-refresh opportunity metadata (not only refresh votes).
- Cadence records now explicitly report PM sync surface and active/inactive kick ownership per event.
- Restart continuation rule remains explicit deterministic rebuild with persisted cadence metadata.

Reproducibility impact:

- Deterministic scheduling semantics are unchanged; this pass improves auditability and explicit operator metadata for PM synchronization and restart continuation.

Observed:

- Added canonical axis-aware config lanes for cosmological box lengths and PM grid shape while preserving scalar compatibility inputs.
- Updated PM periodic spectral operator, per-axis window deconvolution factors, and cell-volume normalization to use `(Lx,Ly,Lz)` and `(Nx,Ny,Nz)` directly.
- Updated reference workflow/runtime topology/restart metadata wiring to carry full PM shape and axis-aware geometry instead of forcing cubic `(N,N,N)` and `(L,L,L)`.
- Added parser compatibility and integration coverage for non-cubic periodic PM geometry.

Reproducibility impact:

- Normalized config output is now canonicalized to axis-aware keys (`box_size_{x,y,z}`, `treepm_pm_grid_n{xyz}`), which removes ambiguity while retaining scalar input compatibility.


## 0) TreePM Phase 3 contract freeze and gravity-stack audit baseline (2026-04-20 UTC)

Commands:

```bash
sed -n '1,220p' docs/treepm_phase2_closeout.md
sed -n '1,260p' docs/treepm_phase2_distributed_contract.md
sed -n '1,320p' docs/pm_gravity_solver.md
sed -n '1,260p' docs/tree_gravity_solver.md
sed -n '1,260p' docs/tree_pm_coupling.md
sed -n '1,220p' docs/time_integration.md
sed -n '1,220p' docs/hierarchical_timestepping_audit_2026-04-06.md
sed -n '1,220p' docs/parallel_distributed_memory_contracts.md
sed -n '1,260p' docs/validation_plan.md
sed -n '1,220p' docs/validation_ladder.md
sed -n '1,220p' docs/releases/known_issues.md
sed -n '1,260p' docs/architecture/decision_log.md
```

Observed:

- Added `docs/treepm_phase3_contract.md` as the Phase 3 contract/maturity-map anchor.
- Recorded explicit audited baseline facts (periodic PM production path, slab-only PM decomposition, cubic PM grid wiring from scalar PMGRID, monopole tree MAC, single-policy softening lane, hierarchy present but Phase 3 integrator-grade PM synchronization unproven).
- Added ADR record in `docs/architecture/decision_log.md` for the Phase 3 contract freeze.
- Updated validation/release docs to keep honest language: Phase 1/2 evidence is not Phase 3 closure; scaling artifacts are not certification.

Interpretation:

- This pass is documentation/infrastructure alignment only and does not change solver behavior, config semantics, or schema payload contracts.

## 0) Phase 2 distributed TreePM closeout failure-contract hardening (2026-04-19 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory"
```

Observed:

- Added explicit restart cadence inconsistency validation (`last_long_range_refresh_opportunity` and `long_range_field_version` coherence checks) during distributed restart-state decode.
- Added explicit restart compatibility-report lane `pm_cadence_steps_match` for invalid cadence metadata.
- Added/extended Phase 2 negative-contract coverage for rank-count mismatch, unsupported decomposition mode, communicator/layout mismatch, missing distributed restart metadata, and inconsistent cadence metadata.
- Added final hard-gate closeout document `docs/treepm_phase2_closeout.md`.

Interpretation:

- This pass does not change force-split numerics; it closes coherence/failure-reporting gaps for distributed topology + restart continuity contracts.

## 0) Phase 2 distributed TreePM validation/gate suite wiring (2026-04-19 UTC)

Commands:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug -j4 --target test_validation_phase2_mpi_gravity generate_mpi_gravity_scaling_artifacts
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -R "validation_phase2_mpi_gravity_two_rank|validation_phase2_mpi_gravity_three_rank"
```

Observed:

- Added an explicit MPI gravity validation gate (`validation_phase2_mpi_gravity_*`) covering:
  - distributed PM vs one-rank force reference (`rel_L2 <= 1e-10`),
  - distributed full TreePM vs one-rank force reference (`rel_L2 <= 5e-6`, `max_rel <= 5e-5`),
  - rank-count reproducibility sweep (`np=1,2,3`),
  - communication-stress path with tiny tree exchange batches,
  - MPI restart write/read continuation verification through reference workflow roundtrip.
- Added separate PM-only and tree-only MPI scaling artifact generators:
  - `validation/artifacts/pm_only_scaling_np{1,2}.csv`
  - `validation/artifacts/tree_only_scaling_np{1,2}.csv`
- CI now surfaces the MPI gravity gate in the MPI+HDF5+FFTW matrix row and uploads scaling artifacts.

Interpretation:

- Phase 2 evidence now has an auditable distributed gravity gate suite rather than relying on pseudo two-rank scaffolding alone.

## 0) TreePM Phase 2 distributed gravity performance hardening pass (2026-04-19 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_pm_solver bench_pm_solver bench_tree_pm_coupling bench_parallel_decomposition_exchange
ctest --test-dir build/cpu-only-debug --output-on-failure -R unit_pm_solver
./build/cpu-only-debug/bench_pm_solver
./build/cpu-only-debug/bench_tree_pm_coupling
./build/cpu-only-debug/bench_parallel_decomposition_exchange
```

Observed:

- PM FFT plans and mesh-side scratch buffers remain cached by slab layout and are now paired with reused distributed interpolation send/receive buffers for force and potential reverse communication paths.
- TreePM short-range distributed exchange now reuses persistent payload/count buffers and overlaps communication with owner-local target evaluation by running local-local residual work while batched request exchange is in-flight.
- Active-set compact sidecar arrays in TreePM coupling now resize without per-step zero-fill churn.
- Bench output now includes repeated warmup/measured iteration fields and PM plan cache counters so reuse behavior is visible in profiler output instead of only one-shot timings.

Interpretation:

- Phase 2 distributed TreePM runtime overhead is reduced via explicit buffer/plan reuse and auditable overlap without changing force-split math or ownership semantics.

## 0) TreePM Phase 2 gravity-aware ownership decomposition + migration commit repair (2026-04-19 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_unit_simulation_state
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|unit_simulation_state"
```

Observed:

- Morton SFC decomposition now supports explicit gravity-aware cost terms for:
  - owned particle count,
  - recent active target count,
  - recent remote tree interaction volume,
  - memory footprint (plus optional retained generic work term).
- Per-rank decomposition metrics now expose each component to make imbalance diagnostics auditable.
- `SimulationState` now provides explicit migration pack/commit boundaries:
  - packs hot lanes + metadata + species sidecars for migrating rows,
  - commits ownership only at one synchronization call,
  - rebuilds species counts and species index,
  - removes stale local ghost/import rows so post-commit ownership/sidecar state is unambiguous.

Interpretation:

- Distributed TreePM ownership has a real migration commit contract and gravity-aware decomposition input signal.
- One-rank baseline behavior is preserved when no migration is staged.

## 0) TreePM Phase 2 distributed workflow cadence-consensus + active-set wiring repair (2026-04-19 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_integration_reference_workflow test_integration_reference_workflow_distributed_treepm_mpi test_integration_tree_pm_coupling_periodic
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_reference_workflow"
./build/cpu-only-debug/test_integration_reference_workflow_distributed_treepm_mpi
```

Observed:

- The live gravity callback path now enforces rank-consensus cadence metadata for each gravity kick:
  - all ranks advance `gravity_kick_opportunity` every gravity kick stage,
  - refresh/reuse decisions are reduced and must agree on all ranks,
  - divergence fails loudly instead of silently drifting counters or field-version metadata.
- Early return on empty local active targets was removed so distributed PM cadence bookkeeping remains coherent even for ranks with no active targets on a given kick.
- Integration coverage now includes:
  - single-rank active-subset agreement against full solve in TreePM coupling tests,
  - distributed active-subset two-rank vs one-rank reference agreement in TreePM coupling tests,
  - distributed workflow MPI smoke test for cadence-record coherence and final digest agreement.

Interpretation:

- Distributed PM cadence and metadata contracts are now wired in the real reference workflow callback path.
- TreePM source/target split is preserved: PM refresh consumes all mass sources, while active-target export/import remains short-range-only work.

## 0) TreePM Phase 2 distributed PM interpolation return path repair (2026-04-19 UTC)

## 0) TreePM Phase 2 distributed short-range tree export/import repair (2026-04-19 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_tree_pm_coupling_periodic
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_tree_pm_coupling_periodic"
```

Observed:

- `TreePmCoordinator::evaluateShortRangeResidual` now runs an explicit active-target export/import protocol when `world_size>1`:
  - owner rank computes local-local residual,
  - owner exports target batches (bounded by `tree_exchange_batch_bytes`) to each peer,
  - peer evaluates requests against its local tree/source data and returns partial accelerations,
  - owner validates response coverage (`batch_token`, `request_id`) and accumulates remote partials.
- The workflow now wires `numerics.treepm_tree_exchange_batch_bytes` into `TreePmOptions`.
- MPI integration coverage was added (`integration_tree_pm_coupling_periodic_mpi_two_rank`) for two-rank distributed-vs-single-rank agreement including cutoff-boundary cross-rank peers.

Interpretation:

- Phase 2 short-range distributed TreePM now uses real peer participation instead of rank-local-only residual evaluation.
- The one-rank numerical contract remains the reference and is used as distributed comparison baseline.

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_pm_periodic_mode
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_pm_periodic_mode"
```

Observed:

- `PmSolver::interpolateForces` now supports slab-distributed reverse communication:
  particle-owner ranks send weighted stencil requests to slab owners, slab owners return weighted force contributions, and owners accumulate in local particle order.
- `PmSolver::interpolatePotential` uses the same reverse message contract for optional potential gather.
- Distributed PM integration coverage now includes CIC and TSC one-rank vs two-rank interpolation agreement with slab-boundary particle cases.

Interpretation:

- Phase 2 PM long-range path is now distributed through deposition, solve, and particle-force/potential return without replicated mesh gather assumptions.
- One-rank PM numerical contract remains the baseline reference and is used as the distributed comparison target.

## 0) TreePM Phase 2 distributed density assignment repair (2026-04-19 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_pm_solver test_integration_pm_periodic_mode
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_pm_solver|integration_pm_periodic_mode"
```

Observed:

- `PmSolver::assignDensity` now supports slab-distributed owner-to-owner contribution routing with
  explicit slab ownership validation and no replicated full mesh requirement on each rank.
- One-rank path remains numerically consistent with prior contract and now also includes TSC wrap/mass checks.
- Integration tests include distributed density assignment agreement against one-rank reference for both CIC and TSC.

Interpretation:

- Phase 2 PM distributed deposition infrastructure is now implemented and tested at contract level.
- Distributed FFT solve support status is unchanged from prior repair state.

## 1) Preset and dependency gates

Commands:

```bash
cat CMakePresets.json
cmake --preset hdf5-debug
cmake --preset pm-hdf5-fftw-debug
```

Observed:

- Presets now include `hdf5-debug` and `pm-hdf5-fftw-debug`.
- Dependency failure behavior is explicit and actionable.
- In this environment, both dependency-enabled presets fail at configure with:
  - `COSMOSIM_ENABLE_HDF5=ON but HDF5 was not found.`

Interpretation:

- Preset-level closure for P05/P06 exists in-tree.
- Runtime closeout for HDF5/PM feature paths is blocked by missing dependency in this environment.

## 2) CPU-only baseline evidence

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure
```

Outcome:

- Configure: **PASS**
- Build: **PASS**
- Tests: **PASS** (`36/36`)

Interpretation:

- Baseline path remains healthy and no unrelated regressions were observed.

## 3) Guardrail scripts

Commands:

```bash
./scripts/ci/check_repo_hygiene.sh
./scripts/ci/guard_feature_paths.sh
```

Outcome:

- `check_repo_hygiene.sh`: **PASS**.
- `guard_feature_paths.sh`: CPU segment passes; feature segment fails at HDF5 configure gate, as expected when dependency is absent.

Interpretation:

- Hygiene and process checks are active and correctly prevent false success claims from CPU-only evidence.

## 4) IC reader gas import state evidence

Command:

```bash
nl -ba src/io/ic_reader.cpp | sed -n '562,592p'
```

Observed:

- For `PartType0`, IC import now writes:
  - `gas_cells.density_code[cell_i]`
  - `gas_cells.internal_energy_code[cell_i]`
- Unsupported reports are for unmapped metallicity/smoothing length, not thermodynamic bypass.

Related test artifact:

- `tests/unit/test_ic_reader.cpp` includes HDF5-gated checks for gas thermodynamic mapping and optional-density defaulting behavior.

Interpretation:

- P04 implementation appears closed in code/tests.
- Full runtime verification remains dependent on HDF5 availability.

## 5) Docs and naming/process hygiene

Commands:

```bash
test -f docs/build_instructions.md
test -f CONTRIBUTING.md
find . -maxdepth 1 -type f | sed 's#^./##' | sort
```

Outcome:

- `docs/build_instructions.md`: present.
- `CONTRIBUTING.md`: present.
- Root file naming is now policy-safe (no legacy violating spreadsheet at repo root).

Interpretation:

- P08/P09/P10 are closed at repository state level.

## 6) Audit conclusion pointer

See `docs/repair_closeout_report.md` for the stop/go decision.
Current state is **STOP** until HDF5 (and then PM/HDF5/FFTW) preset build+test commands pass on the intended feature paths.

## 7) Core boundary repair (reference workflow assembly)

Commands:

```bash
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_core_dependency_direction|integration_reference_workflow"
```

Observed:

- Core dependency-direction guard passes and fails fast on forbidden upward includes in `include/cosmosim/core/**` and `src/core/**`.
- Reference workflow integration smoke test still passes after moving assembly to `workflows/`.

Interpretation:

- The architecture boundary leak from `core` into analysis/I/O/physics workflow assembly is repaired while keeping smoke behavior intact.

## 8) Configuration-contract hardening (typed freeze path)

Commands:

```bash
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_simulation_mode|integration_simulation_mode_toy_runs"
```

Observed:

- Policy-like config values are parsed from param strings, then frozen into typed enums for solver selection, coordinate frame, mode boundaries, and feedback mode/variant.
- Unknown key handling and deprecated alias mapping use a centralized key/alias registry in `src/core/config.cpp`.
- Unit/integration coverage validates unknown key rejection, alias behavior, invalid enum rejection, and deterministic normalized config/hash behavior.

Interpretation:

- String-driven runtime policy drift is reduced after freeze while preserving param-style UX and deterministic normalized-config provenance semantics.

## 9) SimulationState ownership/invariants decomposition repair

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state
./build/cpu-only-debug/test_unit_simulation_state
```

Observed:

- `src/core/simulation_state.cpp` responsibility concentration was split by invariant boundary:
  - storage-lane consistency (`simulation_state_structures.cpp`)
  - ownership checks (`simulation_state_ownership.cpp`)
  - species indexing and transfer packing (`simulation_state_species.cpp`)
  - metadata/module-sidecar serialization (`simulation_state_metadata.cpp`)
  - active-view assembly/scatter (`simulation_state_active_views.cpp`)
  - reorder/scratch logic remains in `simulation_state.cpp`
- Header-level hot-field contract notes now explicitly name allowed gravity/hydro active-view lanes.
- `tests/unit/test_simulation_state.cpp` now covers:
  - ownership invariant failure/recovery
  - unique-ID invariant failure/recovery
  - species-index rebuild correctness
  - transfer packet pack/unpack-equivalence checks
  - metadata serialize/deserialize round-trip
  - active-view hot-field writeback behavior and cold-lane protection
  - static assertion guardrails for gravity/hydro active-view compactness.

Interpretation:

- State ownership and active-view invariants are now reviewable in focused files with explicit hot/cold contracts and command-backed checks.

## 10) Snapshot/restart I/O contract boundary hardening

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_snapshot_hdf5_schema test_unit_restart_checkpoint_schema test_integration_restart_checkpoint_roundtrip
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_snapshot_hdf5_schema|unit_restart_checkpoint_schema|integration_restart_checkpoint_roundtrip"
```

Observed:

- Shared continuation-metadata contract names/validation are centralized via `include/cosmosim/io/io_contract.hpp` + `src/io/internal/io_contract.cpp`.
- Restart docs and API now publish an exact-restart checklist (`exactRestartCompletenessChecklist()`), separating restart-completeness obligations from snapshot interoperability.
- New negative checks cover restart schema mismatch rejection, missing required scheduler dataset rejection, and finalize-failure behavior that leaves target path untouched while preserving the temporary artifact for diagnosis.

Interpretation:

- Snapshot and restart responsibilities are now explicitly separated and reviewable, with stronger error behavior on continuation-critical contract violations.

## 11) Operational observability repair (structured run events)

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_profiling test_integration_reference_workflow
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_profiling|integration_reference_workflow"
```

Observed:

## 12) Phase-2 PM distributed FFT infrastructure repair state

_Date captured: 2026-04-19 (UTC)_

Commands:

```bash
cmake -S . -B build/mpi-fftw-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCOSMOSIM_ENABLE_TESTS=ON -DCOSMOSIM_ENABLE_MPI=ON -DCOSMOSIM_ENABLE_FFTW=ON -DCOSMOSIM_ENABLE_HDF5=OFF
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_pm_solver test_integration_pm_periodic_mode
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_pm_solver|integration_pm_periodic_mode"
```

Observed:

- `PmSolver::solvePoissonPeriodic` now accepts slab-owned local PM grids and uses communicator-aware FFTW MPI plans when MPI+FFTW are enabled.
- PM FFT plan/scratch buffers are cached by slab layout key and reused across calls; unit/integration checks now include plan cache reuse behavior.
- Single-rank PM periodic checks remain green on CPU-only debug path.
- This environment cannot complete MPI validation because CMake cannot find MPI (`Could NOT find MPI_CXX`).

Interpretation:

- Distributed PM FFT ownership/message contract is implemented in-tree with command-backed single-rank regression evidence.
- Multi-rank runtime closure remains blocked in this container by missing MPI runtime/development toolchain.

## 12) Infrastructure gate-bundle enforcement hardening

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
bash -n scripts/ci/enforce_infra_gates.sh
```

Observed:

- Added `scripts/ci/enforce_infra_gates.sh` to run three explicit infrastructure gates with command failure propagated per path:
  - `cpu_core_boundary_and_config_contract`
  - `hdf5_schema_and_exact_restart_contract`
  - `pm_hdf5_fftw_feature_path_validation`
- Each gate writes artifacts under `ci_artifacts/infrastructure_gates/<gate_id>/`.
- A machine-readable status file, `infrastructure_gate_report.json`, now records gate ids, presets, test regexes, and pass/fail status for CI interpretation.
- `.github/workflows/ci.yml` now includes a dedicated `infrastructure_gates` job and makes `reproducibility_gate` depend on both matrix coverage and the explicit gate bundle.
- `tests/integration/test_core_dependency_direction.cmake.in` now checks both source/header include direction and `CMakeLists.txt` target-link direction for `cosmosim_core`.

## 13) Phase 2 distributed TreePM contract freeze surface

_Date captured: 2026-04-19 (UTC)_

Commands:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug -j4
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -R "unit_config_parser|integration_docs_scaffold|integration_provenance_roundtrip|unit_parallel_distributed_memory"
```

Observed:

- Typed config now freezes Phase 2 distributed gravity controls:
  - `numerics.treepm_pm_decomposition_mode` (`slab` only),
  - `numerics.treepm_tree_exchange_batch_bytes` (`>0`).
- Provenance now records those controls for audit continuity.
- New architecture contract doc exists at `docs/treepm_phase2_distributed_contract.md`.
- CI matrix now includes `mpi-hdf5-fftw-debug`, so distributed gravity development is not only gated through optional `mpi-release` smoke.

Interpretation:

- Contract/documentation/build surfaces for Phase 2 are now explicit without claiming distributed TreePM algorithm completion.

Interpretation:

- CPU-only success is no longer treated as implicit feature-path closure because HDF5 and PM/HDF5/FFTW are evaluated in one explicit gate bundle with a consolidated status report.
- Architecture-boundary enforcement now covers both include-level and build-graph-level dependency direction.

## 12) Parallel contract hardening follow-up (R07)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|integration_parallel_two_rank_restart"
```

Observed:

- Distributed-memory transfer roles now expose typed `outbound_transfers` / `inbound_transfers` descriptors in addition to legacy send/recv index vectors.
- Plan/buffer invariants were tightened for invalid ownership combinations, zero-byte ghost payload shape, and decode shape checks.
- Reduction agreement fields now distinguish deterministic baseline vs measured sum, plus a policy helper for absolute/relative tolerance gating.
- Config consensus mismatch artifacts include property-level mismatch rows with baseline/offending rank-value pairs.

Interpretation:

- Pseudo-multi-rank contract clarity and diagnosability improved for review/CI artifacts without adding new MPI features or claiming production multi-rank closure.

## 14) Pencil PM end-to-end transposed path hardening

_Date captured: 2026-04-20 (UTC)_

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_config_parser
./build/cpu-only-debug/test_unit_config_parser
```

Observed:

- `numerics.treepm_pm_decomposition_mode` has typed enum coverage for `slab` and `pencil`.
- `pencil` now executes through explicit FFTW-MPI transposed plans (`TRANSPOSED_OUT/IN`) rather than a slab alias.
- Real-space ownership remains slab for deposition/interpolation while spectral ownership is explicitly transposed by y-range in pencil mode.
- Restart/runtime topology now carry decomposition mode metadata instead of forcing `"slab"`.

Interpretation:

- This patch closes the previous prompt-violation gap by providing an end-to-end non-aliased pencil PM mode with explicit transpose ownership semantics.

## 13) Parallel contract hardening follow-up (R07c)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_integration_parallel_two_rank_restart
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|integration_parallel_two_rank_restart"
```

Observed:

- Ghost transfer descriptors now encode explicit lifecycle intent (`ghost refresh request`, `ghost refresh receive staging`) and post-transfer residency expectation, with reserved typed migration intents documented as out-of-scope scaffolding placeholders.
- `validateGhostExchangePlan()` now enforces neighbor-slot/peer-rank role correctness and exact descriptor-index equality against canonical send/recv vectors.
- Reduction agreement checks now use explicit mode selection (`absolute only`, `relative only`, `absolute and relative`, `absolute or relative`) rather than implicit boolean logic.
- Unit coverage now includes plan-drift rejection checks and per-policy-mode reduction agreement assertions.

Interpretation:

- Pseudo-multi-rank distributed-memory contracts are now more reviewable and diagnosable without claiming production MPI transfer correctness or migration commit semantics.

## 14) Diagnostics maturity-tier repair (analysis honesty/scalability guard)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_config_parser test_unit_analysis_diagnostics test_integration_analysis_bundle
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_analysis_diagnostics|integration_analysis_bundle"
```

Observed:

- Analysis config now carries typed `diagnostics_execution_policy` with normalized-config persistence.
- Diagnostics bundles now emit per-diagnostic maturity metadata (`tier`, `maturity`, `scalability`, execution policy).
- Provisional heavy reference diagnostics (power spectrum direct summation) no longer run under default policy.
- Unit/integration tests now assert that heavy provisional diagnostics only run when explicitly opted in.

Interpretation:

## 13) Distributed-memory ownership/determinism contract hardening

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_integration_parallel_two_rank_restart
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|integration_parallel_two_rank_restart"
```

Observed:

- Distributed-memory ghost planning now has an explicit typed ownership contract (`LocalGhostDescriptor` with `kOwned`/`kGhost`) and validation for invalid ownership combinations.
- Deterministic reduction agreement helpers now provide explicit deterministic reference sums and absolute/relative agreement reporting for reproducibility checks.
- Multi-rank config-freeze consensus checks now compare normalized config hash + rank-count expectation + deterministic-reduction mode across per-rank digests.
- Unit/integration tests cover typed ownership planning, deterministic reduction agreement checks, and config-consensus checks using pseudo-multi-rank vectors.

Interpretation:

- Reviewability and contract clarity improved for rank-owned vs ghost state and reduction determinism semantics without expanding scope into new decomposition/scaling features.
- Evidence remains CPU-only/pseudo-multi-rank in this environment; no claim is made that this alone closes full production MPI execution.

- Infrastructure run-health counters remain first-class and cheap.
- Validated light science diagnostics remain available in default runs.
- Reference/provisional heavy diagnostics are quarantined behind explicit non-default policy.

- `core::ProfilerSession` now carries a minimal structured runtime event model (kind, severity, subsystem, optional step/time/scale context, message, key/value payload).
- Reference workflow emits a machine-readable operational report (`reference_operational_events.json`) linked to deterministic config provenance via `provenance_config_hash_hex`.
- Key infrastructure lifecycle/failure surfaces are explicit in event records (config freeze validation, restart/snapshot write/read begin/complete/failure).

Interpretation:

- Operational troubleshooting and CI artifact review no longer depend on ad hoc text/exception surfaces alone.
- Reproducibility posture remains unchanged: operational events are additive and provenance-linked; no solver behavior or restart/snapshot schema semantics were changed.

## 13) Infrastructure gate hardening follow-up (R08 durability/diagnostics)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
bash -n scripts/ci/run_preset_pipeline.sh
bash -n scripts/ci/enforce_infra_gates.sh
cmake --preset cpu-only-debug
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_core_dependency_direction"
bash ./scripts/ci/enforce_infra_gates.sh ci_artifacts/local_infra_gates
```

Observed:

- Infrastructure gates are now defined in `scripts/ci/infrastructure_gates_manifest.tsv` with explicit fields for gate id, preset trio, intended test scope, metadata expectation, artifact subdir, and benchmark toggle.
- `run_preset_pipeline.sh` now emits `preset_pipeline_report-<preset>.json` with per-phase command metadata and failure phase (`configure`, `build`, `test`, `artifact_collection`, `benchmark`) to improve triage fidelity.
- `enforce_infra_gates.sh` now records richer per-gate entries in `infrastructure_gate_report.json`, including `build_preset`, `test_preset`, `test_scope`, `artifact_dir`, `failed_phase`, and explicit phase commands.
- `integration_core_dependency_direction` now validates target-link direction by inspecting CMake-generated Graphviz dependency metadata from the configured build tree instead of scanning only top-level `CMakeLists.txt` text.

Interpretation:

- Gate semantics are less rename-fragile and more reviewable due to a narrow declarative manifest.
- CI failure triage can identify failing phase and command per gate without log archaeology.
- Core target-link direction guard is stronger than single-file regex scanning while remaining narrow and infrastructure-focused.


## 11) 2026-04-14 targeted infrastructure repair pass

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j1 --target   test_unit_config_parser   test_unit_simulation_state   test_unit_parallel_distributed_memory   test_unit_restart_checkpoint_schema
ctest --test-dir build/cpu-only-debug --output-on-failure -R   "unit_config_parser|unit_simulation_state|unit_parallel_distributed_memory|unit_restart_checkpoint_schema"

cmake --preset hdf5-debug
ninja -C build/hdf5-debug -j8   test_unit_restart_checkpoint_schema   test_integration_snapshot_hdf5_roundtrip   test_integration_restart_checkpoint_roundtrip
ctest --test-dir build/hdf5-debug --output-on-failure -R   "unit_restart_checkpoint_schema|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip"

cmake --preset pm-hdf5-fftw-debug
```

Observed:

- CPU targeted repair tests pass after hardening config round-trip/hash semantics, species-sidecar coverage invariants, descriptor-only ghost-plan honesty, and distributed restart decode completeness checks.
- HDF5 targeted repair tests pass after fixing restart stellar-sidecar completeness, restart payload hash coverage, continuation metadata cross-checking, and `Header/MassTable` snapshot-mass fallback.
- `pm-hdf5-fftw-debug` configure is blocked in this environment because FFTW3 development files are unavailable.

Interpretation:

- The repaired invariants are command-backed on CPU and HDF5 feature paths.
- PM/FFTW path remains an environment blocker here, not a demonstrated code regression in this repair pass.

## 14) 2026-04-19 Phase 2 PM slab ownership/storage milestone

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_unit_pm_solver
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|unit_pm_solver"
```

Observed:

- Added explicit PM slab ownership typing (`parallel::PmSlabLayout`, `PmSlabRange`) with deterministic uneven partitioning, owner-rank lookup, and validated global/local index conversions.
- `gravity::PmGridStorage` now accepts explicit slab layout and allocates only local slab storage (`local_nx * Ny * Nz`), while default construction remains the one-rank full-domain slab.
- PM solver entry points now fail fast on partial slabs, preventing pseudo-distributed use until distributed FFT and remote deposition/gather are implemented.
- Unit tests now cover uneven slab partitioning, ownership queries, index round-trips, one-rank slab equivalence, and explicit rejection of partial-slab use on the single-rank solver path.

Interpretation:

- This is an infrastructure-only Phase 2 milestone: ownership/storage contracts are now auditable without claiming distributed PM algorithm completion.
- Reproducibility posture is preserved for one-rank runs because default PM storage maps to the same full-domain indexing contract.


## 15) 2026-04-19 MPI/CUDA runtime topology milestone

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_unit_config_parser
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|unit_config_parser"
```

Observed:

- Added explicit runtime topology assembly (`parallel::DistributedExecutionTopology`) that binds together MPI world size/rank, PM slab ownership, and rank-local CUDA device assignment.
- Added `core::CudaRuntimeInfo` and explicit CUDA device-selection helpers so `parallel.gpu_devices` becomes a real runtime contract rather than a dead config field.
- Reference workflow TreePM initialization now validates `parallel.mpi_ranks_expected` against the runtime world and enables the CUDA PM assignment/interpolation path only when the runtime request is valid.
- GPU requests now fail loudly when CUDA was requested but no visible devices are available; this avoids silent CPU fallback drift in distributed runs.

Interpretation:

- This remains infrastructure work, not a claim of distributed PM FFT or MPI+GPU overlap completion.
- The code path is now materially closer to an operational multi-rank/multi-device TreePM bring-up because rank-local execution intent is explicit and auditable.

## 16) 2026-04-19 distributed restart/provenance continuation-safety milestone

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_restart_checkpoint_schema test_unit_parallel_distributed_memory test_integration_restart_checkpoint_roundtrip test_integration_parallel_two_rank_restart
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_restart_checkpoint_schema|unit_parallel_distributed_memory|integration_restart_checkpoint_roundtrip|integration_parallel_two_rank_restart"
```

Observed:

- Restart schema bumped to `cosmosim_restart_v4` with explicit `/distributed_gravity/state` payload (`parallel::DistributedRestartState`, schema_version=2).
- Distributed continuation payload now persists decomposition epoch, owning-rank table, PM slab layout metadata, gravity-kick cadence state, long-range field refresh/version metadata, and explicit `long_range_restart_policy`.
- Policy is explicit and enforced: `deterministic_rebuild` (cached long-range PM field arrays are not serialized; continuation rebuilds deterministically on cadence refresh).
- Provenance schema bumped to `provenance_v5` to add axis-aware PM grid and mesh-spacing metadata while preserving legacy scalar aliases for compatibility.
- Added typed compatibility diagnostics (`evaluateDistributedRestartCompatibility`) plus negative mismatch coverage.

Interpretation:

- This is an infrastructure-repair continuation-safety improvement for Phase 2 distributed TreePM, not a new physics model.
- Reproducibility posture remains explicit: restart carries deterministic policy + auditable rank/layout metadata, and integrity hashing now includes distributed continuation state.

- Added isolated/open PM long-range path selection in `reference_workflow` and `TreePmCoordinator`; isolated now uses doubled-domain free-space convolution (padded FFT-style linear convolution) and disables minimum-image short-range wrapping. Current limitation remains single-rank only, with explicit fail-fast on multi-rank.

## 17) 2026-04-20 Tree softening maturity: species + optional per-particle sidecars

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_tree_gravity test_integration_tree_gravity_vs_direct test_integration_snapshot_hdf5_roundtrip test_integration_restart_checkpoint_roundtrip
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_tree_gravity|integration_tree_gravity_vs_direct|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip"
```

Observed:

- Tree gravity now resolves softening at source/target with optional sidecars (species table and per-particle overrides) and uses an explicit pair rule `epsilon_pair = max(epsilon_i, epsilon_j)` across leaf P2P and accepted-node multipole paths.
- TreePM short-range residual uses the same pair law for local and distributed remote request/response paths (request packet now carries target epsilon).
- Optional per-particle softening sidecar persistence is added to snapshot (`/PartTypeN/GravitySofteningComoving`) and restart (`/state/particle_sidecar/gravity_softening_comoving`) paths.
- Schema versions are intentionally bumped for additive axis-aware PM metadata: snapshot `gadget_arepo_v4`, provenance `provenance_v5`, while restart remains `cosmosim_restart_v12`.

Interpretation:

- This is an infrastructure-contract maturation pass: explicit pairwise softening semantics are now documented, tested, and consistent between tree-only and TreePM residual evaluation paths.

## 18) 2026-04-20 distributed gravity clustered-load balancing and exchange observability hardening

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_validation_phase2_mpi_gravity
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory"
```

Observed:

- Morton SFC decomposition now chooses rank cut points by nearest target-prefix boundary once crossed, reducing clustered-load overshoot while preserving contiguous deterministic ownership ranges.
- Distributed TreePM short-range peer gating keeps exact force behavior but tightens periodic source-support bounds using wrap-aware compact intervals, reducing unnecessary remote request export volume for clustered data near periodic boundaries.
- Distributed short-range diagnostics now expose request/response bytes, target request coverage, geometry-pruned pair counts, and request pressure imbalance metrics (max-peer + ratio), enabling clustered-run communication auditing.
- Existing send/recv payload/displacement buffers continue to be reused from persistent coordinator workspace across batches/steps.

Interpretation:

- This is infrastructure repair for distributed maturity and observability under clustered loads; no force-kernel heuristics, dropped contributions, or mid-kick ownership repartitioning were introduced.
- Remaining maturity limits are tracked in `docs/repair_open_issues.md` (see open entry for clustered decomposition/short-range scheduling limits).



## 2026-05-08 — Preset alias and schema-drift gate repair

- Root cause for the failed `cmake --preset cpu-debug` command was a naming drift: build/test presets used `build-cpu-debug` and `test-cpu-debug`, but the configure preset existed only as `cpu-only-debug`. Added `cpu-debug` as a compatibility configure alias that maps to the canonical CPU-only debug build directory and metadata.
- Root cause for the docs/schema gate failure was stale restart-schema literals in release artifacts and gate checks after the runtime schema moved to `cosmosim_restart_v12`. Updated the release manifest, release docs, and release readiness check to the current restart schema.
- Reproducibility impact: no solver, snapshot, restart payload, HDF5 layout, or normalized-config behavior changed; this repair only aligns preset names and schema documentation/tests with existing runtime truth.

## 2026-05-08 — Gas-cell identity seam RFC and isolated API guard

- Added `docs/architecture/gas_cell_identity_map_rfc.md` to define the proposed `GasCellIdentityMap` seam for future AMR/moving-mesh gas-cell ownership without changing production hydro or restart behavior.
- Added an isolated core API shape that validates unique stable `gas_cell_id` and unique transient `local_cell_row` while allowing multiple gas cells to share an optional parent particle or have no parent.
- Preserved the current P0 particle-bound gas contract: production hydro/restart paths still rely on `requireParticleBoundGasCellContract(...)` until an explicit schema migration and solver remap test plan replaces it.

## 2026-04-20: Phase 3 gravity evidence campaign scaffolding

- Added campaign-specific validation configs for cosmological box and zoom runs under `configs/validation_phase3/`.
- Added research-grade observable contract in `validation/campaign/phase3_gravity_observables.json`.
- Added provenance collector `scripts/validation/collect_phase3_evidence.py` and campaign guide `docs/validation_phase3_campaign.md`.
- Validation convergence now emits deterministic force-profile and time-integration artifacts under `validation/artifacts/research_grade/phase3/`.

## 0) Phase 3 final integration closeout audit (2026-04-21 UTC)

Commands:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --test-dir build/pm-hdf5-fftw-debug -R "integration_tree_gravity_vs_direct|integration_tree_pm_coupling_periodic|integration_reference_workflow_distributed_treepm_mpi|validation_phase2_mpi_gravity_single_rank|validation_convergence|integration_reference_workflow$|integration_time_integration_loop" --output-on-failure
python3 scripts/validation/collect_phase3_evidence.py
```

Observed:

- `mpi-hdf5-fftw-debug` configure is blocked in this environment by missing FFTW MPI library (`fftw3_mpi`).
- `pm-hdf5-fftw-debug` configure/build succeeds.
- Targeted gravity/Phase 3 coherence suite is not green in this cycle:
  - `integration_reference_workflow` fails with runtime workflow schema compatibility validation error.
  - `integration_tree_pm_coupling_periodic` fails the zoom low-res contamination invariant.
  - `validation_phase2_mpi_gravity_single_rank` fails communication-stress residual-cutoff expectation despite zero force-error norm.
- Manifest/provenance collector ran successfully and refreshed Phase 3 campaign metadata to current HEAD.

Interpretation:

- Phase 3 cannot be marked closed in this audit cycle.
- Closeout verdict is recorded in `docs/treepm_phase3_closeout.md` as **incomplete** with explicit blockers and forbidden claims.

Reproducibility impact:

- Documentation/evidence-manifest refresh only; no solver, schema, or config semantic changes.

## 2026-04-26: Runtime identity, sidecar, and schema-drift hardening

Commands:

```bash
cmake --preset cpu-only-debug -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS_DEBUG=-g0
cmake --build --preset build-cpu-debug --target test_integration_species_migration_invariants test_unit_simulation_state test_integration_reorder_compaction_sidecars test_integration_softening_ownership_invariants test_unit_snapshot_hdf5_schema -j1
./build/cpu-only-debug/test_unit_simulation_state
./build/cpu-only-debug/test_integration_reorder_compaction_sidecars
./build/cpu-only-debug/test_integration_species_migration_invariants
./build/cpu-only-debug/test_integration_softening_ownership_invariants
./build/cpu-only-debug/test_unit_snapshot_hdf5_schema
```

Observed:

- Particle migration commits now validate every inbound record before mutating state: species tags must be valid, sidecar payload flags must exactly match the species contract, inbound owner rank must match the receiving commit rank, and the final kept-plus-inbound particle-ID set must remain unique.
- Ownership invariants now include particle-ID uniqueness, so duplicate identity state fails loudly instead of surviving as a later restart/I/O/migration corruption.
- Active gravity/hydro kernel-view scatter now protects the global generation registry with a mutex and rejects unregistered or stale compact views rather than silently accepting a pointer not produced by the gather path.
- Snapshot HDF5 species mapping now rejects invalid species tags and invalid `PartType` indices instead of implicitly falling back to dark matter.
- Snapshot schema drift test/documentation was aligned to the runtime additive schema identity `gadget_arepo_v4` / schema version 4, closing the stale `gadget_arepo_v2` unit-test expectation.

Interpretation:

- This is a blocker-oriented corruption hardening pass. No force-kernel heuristics, cosmological parameters, gravity split scales, hydrodynamic flux formulas, or subgrid calibration constants were changed.
- HDF5 roundtrip with real HDF5 was not rerun in this CPU-only environment; the CPU-only schema contract and identity/sidecar regression binaries passed directly.


## 2026-05-07: Migration/compaction identity contract hardening

Commands run for validation:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_integration_species_migration_invariants test_integration_parallel_two_rank_restart
cmake --build --preset build-cpu-debug -j4
./build/cpu-only-debug/test_integration_species_migration_invariants && ./build/cpu-only-debug/test_integration_parallel_two_rank_restart
ctest --preset test-cpu-debug --output-on-failure -R "migration|parallel|two_rank|sidecar|species"
cmake --preset mpi-hdf5-fftw-debug
```

Observed repair results:

- Particle migration stale-ghost inputs now reject duplicate stale indices and local-owned rows, preserving the separation between outbound owned migration and remote-owned stale ghost removal.
- Species migration regression coverage now includes wrong inbound owner rank, stale ghost sidecar removal with exact generation invalidation, duplicate/owned stale ghost errors, and stale sidecar index errors.
- The distributed-memory contract now states the migration record completeness contract for hot lanes, common sidecars, softening override truth, and species-specific sidecars, and documents gas-cell rebuild by stable gas particle ID rather than old local positions.
- The targeted migration and two-rank restart binaries pass directly after rebuild. Follow-up fixture repair aligned the broader regex-selected hot/cold and mixed-species tests with the particle-bound gas-cell identity contract, so the requested CPU ctest selection is expected to exercise only valid gas-cell fixtures.
- The MPI+HDF5+FFTW preset remains environment-blocked by missing `fftw3_mpi` development support.

Reproducibility impact:

- No solver numerics, HDF5 snapshot/restart schema, config keys, or science workflow features changed. The repair tightens rank-local ownership validation and view invalidation evidence for particle index-space mutations.

## 2026-05-07 Restart sidecar exactness repair

- Restart schema is now `cosmosim_restart_v12`: readers reject v5 and older files unless a future explicit compatibility importer is added and tested.
- Restart HDF5 always writes and requires the particle softening value and authoritative override-mask datasets, preserving exact sidecar presence instead of silently defaulting missing lanes.
- Restart payload integrity hashing now length-delimits strings and vector byte lanes and continues to cover particle identity, gas identity, species sidecars, module sidecars with schema versions, scheduler/integrator state, normalized config/provenance, and distributed gravity continuation metadata.
- Reproducibility impact: restart continuation is stricter and more auditable; old files missing v6 runtime-truth lanes are rejected rather than resumed with inferred defaults. Snapshot canonical GADGET/AREPO dataset names remain unchanged.

## 2026-05-08: Distributed ownership identity reduction floor

Commands planned for validation:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_integration_reference_workflow_distributed_treepm_mpi test_integration_parallel_two_rank_restart
./build/cpu-only-debug/test_unit_parallel_distributed_memory
ctest --preset test-cpu-debug --output-on-failure -R "parallel|distributed|two_rank|restart"
cmake --preset mpi-hdf5-fftw-debug && cmake --build --preset build-mpi-hdf5-fftw-debug -j4 && ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -R "mpi|distributed|two_rank|restart"
```

Observed repair scope:

- Added a distributed ownership identity summary helper for rank-local authoritative particle IDs.
- Reference workflow reports now expose local duplicate-ID status and a reduced partition identity check for count/sum/xor auditing.
- Distributed TreePM MPI regression assertions now check the report-level identity status and independently compare reduced identity tuples to the generated IC set.
- Documentation now distinguishes correctness-first identity reductions from mature load balancing or pencil-FFT work.

Reproducibility impact:

- No solver numerics, config keys, restart schema version, or HDF5 dataset layout changed.
- Restart long-range PM policy remains deterministic rebuild; the repair adds deterministic integer identity reductions around existing rank-local state.

## 2026-05-08: Stage 1 runtime-truth CI gate consolidation

Commands planned for validation:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-stage0-runtime-truth-cpu-debug --output-on-failure
ctest --preset test-stage1-runtime-truth-cpu-debug --output-on-failure
./scripts/ci/run_stage1_runtime_truth_gate.sh ci_artifacts/stage1_runtime_truth
```

Observed repair scope:

- Added explicit `runtime_truth`/`p0` CTest labels and subsystem labels for the CPU-runnable Stage 1 suite covering softening, sidecar, gas identity, migration, restart, provenance, scheduler mirrors, and active views.
- Added `integration_runtime_truth_ctest_labels` to audit P0 registration/labels and to fail CPU-only builds if feature-gated HDF5 app-smoke, MPI multi-rank, CUDA, or Python tests are registered while their feature is disabled.
- Added `test-stage1-runtime-truth-cpu-debug` and kept `test-stage0-runtime-truth-cpu-debug` as a compatibility alias over the same exact P0 test-name set.
- Added `scripts/ci/run_stage1_runtime_truth_gate.sh` and a manifest row so local and CI runs have one documented dependency-free Stage 1 gate command.

Reproducibility impact:

- Test registration, labels, presets, scripts, and documentation only. No solver numerics, HDF5 schema names, restart schema versions, normalized config dumps, provenance payload contents, output naming, deterministic scheduling behavior, or SoA/hot-cold state layout changed.


## 2026-05-10 — Stage 2 timestep ownership audit

- Added `docs/architecture/stage2_timestep_ownership_audit.md` as the pre-repair inventory for timestep ownership lanes, serialization classifications, mutation paths, unsafe hidden state, and follow-up repair prompts.
- Recorded `ADR-INFRA-STAGE2-TIMESTEP-OWNERSHIP-014` in `docs/architecture/decision_log.md`: the hierarchical scheduler is the intended per-element timestep authority; state `time_bin` arrays are mirrors; integrator time-bin context is metadata; the pre-repair audit classified TreePM cadence as a gravity workflow/distributed restart lane; the 2026-05-11 ADR refines live cadence authority to `PmSynchronizationState` while keeping restart/provenance cadence fields serialized metadata.
- No solver math, restart schema, output schema, or runtime behavior changed in this audit patch.

## 2026-05-10 Scheduler timestep authority repair

- Hierarchical scheduler APIs now accept labelled candidate timestep/bin submissions and reconcile them before legal transition commit; production adaptive hydro/gravity criteria submit constraints instead of mutating pending bin transitions directly.
- `StepOrchestrator` now rejects expected-scheduler-tick execution when active-set descriptors lack scheduler provenance, so production solver callbacks cannot consume mirror-derived active sets silently.
- PM long-range cadence decisions are represented by scheduler-owned synchronization state (`PmSynchronizationState`), preserving cadence > 1 refresh boundaries while leaving gravity force execution in the gravity callback.
- Reproducibility impact: integer timeline authority remains in scheduler metadata; particle/cell `time_bin` lanes remain mirrors refreshed after scheduler commit. No snapshot/restart schema or config keys changed.
- Validation note: CPU-only Stage 1 runtime-truth tests and the full CPU debug preset pass after aligning tracer/star-formation/black-hole/simulation-state-pipeline fixtures with the particle-bound gas-cell identity contract.

## 2026-05-10 Restart/HDF5 regression follow-up

- Root cause for the HDF5 restart roundtrip abort was a malformed reorder fixture: `state.particle_sidecar.sfc_key` was reassigned with six entries for a seven-particle state, so `reorderAlignedVector` correctly rejected the lane/permutation length mismatch before restart assertions could complete.
- HDF5 missing-field checks now preflight dataset existence before `H5Dopen2`, preserving explicit missing-dataset exceptions for required v6 fields (`pending_bin_index`, `has_gravity_softening_override`) without emitting misleading HDF5 diagnostics in negative tests.
- Reproducibility impact: no restart schema or compatibility semantics changed; this only fixes fixture lane cardinality and makes missing required-field diagnostics deterministic.

## 2026-05-10 Hierarchical timestep invariant traps

- Added scheduler-owned runtime invariant checks around substep begin/end, candidate reconciliation, bin-transition commit, active-set descriptor consumption, restart timestep-state import/export validation, time-bin mirror consumption, and PM cadence refresh commits.
- Bad states now fail with deterministic context strings carrying the source label, element index, current integer tick, bin, next activation tick, placeholder integrator step/time context, and PM sync context where relevant.
- Added CPU debug tests for early activation, illegal bin jumps, skipped PM refresh commits, stale mirror use, invalid restart timestep payloads, active-set mismatch, and local/global synchronization-boundary violations.
- Reproducibility impact: no solver numerics, config keys, restart schema version, HDF5 dataset layout, normalized config dumps, provenance payloads, or output naming changed. The patch hardens deterministic scheduler-state validation and turns previously silent invalid multirate states into immediate test-visible failures.

## 2026-05-10 PM validation residual cutoff contract repair

- Fixed the PM+HDF5+FFTW validation residual-cutoff gate so `validation_phase2_mpi_gravity_single_rank` accepts either internal-node cutoff pruning or leaf pair cutoff skips as deterministic evidence that the TreePM residual cutoff path was exercised.
- Reproducibility impact: validation-only contract alignment; no solver numerics, restart/schema payloads, config keys, provenance, output naming, or runtime scheduling behavior changed.

## 2026-05-10 structural-transform scheduler remap repair

- Added stable-identity scheduler remap/rebuild helpers for particle reorder, compaction/migration-style row rebuilds, and scheduler identity-record imports. Gas-cell time-bin mirror refresh now resolves parent particle identity instead of assuming cell rows and scheduler rows are positionally coupled.
- Reproducibility impact: deterministic activation by `particle_id` is preserved across local structural transforms when callers use the remap/rebuild helpers. No solver numerics, config keys, restart schema version, HDF5 dataset layout, normalized config dumps, provenance payloads, or output naming changed. The repair tightens failure behavior: missing scheduler identity records now throw instead of silently trusting `time_bin` mirrors.
- Future MPI note: local helpers do not claim fully distributed scheduler migration. Multi-rank exact continuation still requires exchanging full scheduler identity records with migrating particles and coordinating scheduler tick/max-bin metadata across ranks.

## 2026-05-11 Stage 2 verification gate hardening

- Re-audited Stage 2 scheduler authority coverage for direct timestep mirror writes, scheduler hot metadata, next activation, active flags, current ticks, integrator time-bin metadata, PM cadence fields, restart hash/read validation, and docs/ADR consistency.
- Promoted `unit_time_integration` into the runtime-truth/P0 label set and added `test-stage2-runtime-truth-cpu-debug` so the scheduler authority, split-brain removal, PM sync legality, restart equivalence, and CTest label audit tests have an explicit CPU-runnable preset gate.
- Extended the HDF5 restart roundtrip negative coverage to corrupt `/scheduler/active_flag` in-file and require restart reader rejection, complementing stale mirror and payload-hash tamper checks.
- Reproducibility impact: no solver numerics, config keys, restart schema version, dataset layout, normalized config output, provenance payload, or deterministic scheduling semantics changed. The repair only tightens test registration and negative restart validation coverage.

## 32) 2026-05-11 Stage 2 adversarial timestep restart mirror audit

- Reviewed Stage 2 timestep authority surfaces with hostile searches for writable timestep references, active-set construction from mirrors, stale cache bypasses, restart/migration transforms, and AGENTS/code-review guardrails.
- Patched the blocker restart loophole where particle-bound gas-cell `CellSoa::time_bin` mirrors were validated only when scheduler `bin_index` happened to be cell-sized. Restart validation/import now maps each cell through its parent gas particle and scheduler `bin_index`, preserving scheduler authority for mixed particle/gas states.
- Reproducibility impact: no solver behavior, schema fields, normalized config, or provenance format changed. Existing restart v6 payloads with stale gas-cell timestep mirrors now fail fast instead of importing fake mirror truth; valid payloads rebuild cell mirrors deterministically from scheduler state.
- Added the adversarial review artifact at `docs/architecture/stage2_adversarial_timestep_review.md` and documented the compatibility tightening in `docs/output_schema.md`.

## 2026-05-19 PM refresh first-class stage diagnostics hardening

- Extended integrator-owned `PmRefreshDirective` metadata so PM refresh diagnostics now include explicit refresh reason (`initial_force_bootstrap` vs `scheduled_force_refresh_stage`) and explicit PM force-evaluation scale factor from the KDK stage boundary.
- Wired the new directive fields at orchestrator stage dispatch points (`kGravityKickPre`, `kForceRefresh`) and propagated them into TreePM runtime events for stage-auditable PM cadence evidence.
- Added a targeted unit test to prove both explicit PM refresh reasons are emitted through stage dispatch and carry force-evaluation context.
- Reproducibility impact: no solver numerics, no config keys, no restart schema changes, and no output naming changes. This is observability/legality hardening of existing integrator-owned PM cadence behavior.

- S5-02: normalized numerics time/cosmology keys now use canonical names (`a_*`, `z_*`, `t_code_*`, `t_phys_*`, `integrator_time_variable`); legacy ambiguous keys remain UserConfig aliases only.

## 2026-05-29 Stage 8 restart-safe boundary contract

- Restart checkpoint requests now use an explicit `RestartBoundaryDecision` contract before hash/write traversal. Unsafe KDK half-steps, local active-bin substeps, non-restart-safe boundary kinds, and uncommitted PM refresh transitions fail with deterministic diagnostics instead of being silently skipped or serialized.
- Reproducibility impact: no solver numerics, config keys, restart schema version, HDF5 dataset names, normalized config dumps, provenance payloads, or output naming changed. Exact continuation is hardened by preventing invalid phase states from becoming restart truth.

## 2026-07-13 periodic TreePM and DMO gravity hardening

Current gravity status is consolidated in
`docs/gravity_production_readiness.md`; dated Phase 2/Phase 3 closeout records
remain historical and are not rewritten.

Implemented repair outcomes:

- Periodic TreePM builds topology, COMs, corrected quadrupoles, raw second
  moments, node bounds, MAC geometry, and cutoff pruning in a coherent per-axis
  largest-gap unwrapped frame. Rectangular seam and integer-image translation
  tests cover leaf/internal and monopole/quadrupole behavior.
- Standalone and screened multipoles use the corrected acceleration sign for
  the repository displacement convention. Screened/softened nodes use the raw
  second-moment trace and derivatives of the full radial kernel.
- The third `relative_force_error` MAC is wired through typed config,
  normalization, validation, runtime, committed force-history fallback,
  `provenance_v6`, restart audit metadata, unit tests, and integration tests.
- Empty TreePM ranks emit explicit zero-source roots; independent target
  coordinates support source-empty target owners; all-empty, source-only,
  target-only, zero-record, and migration-to-empty epochs are exercised.
- Tree hierarchy and short-range request/response records use explicit
  versioned little-endian wire formats with peer/record identity,
  decomposition/force/exchange epochs, finite/bounds checks, and exact response
  coverage.
- Periodic unwrapping and local-tree construction exceptions are now caught and
  voted on the actual MPI world before any residual collective. PM-layout and
  context metadata mismatch is rejected by solve-entry consensus rather than a
  rank-local constructor throw. Reusable TreePM exchange workspace now owns
  count/displacement, payload, remote-acceleration, and response-coverage lanes;
  throwing capacity growth occurs in coordinated protocol phases and all of
  those capacities appear in transient MPI memory accounting.
- PM halo-cache commit, compact active-position/PM-force/zoom-correction target
  preparation, and worst-case residual traversal-stack reservation are likewise
  caught and voted before ranks can enter subsequent collectives. Reserving the
  full local-node stack bound removes traversal-time growth allocations.
- PM density and force/potential request/response records also use explicit
  versioned little-endian `PMW1` formats. Actual-MPI-world preflight votes API
  operation and serial/distributed layout mode before any layout branch.
- PM slab authority now matches FFTW-MPI fixed ceil-block ownership. Periodic
  halos use nonblocking receive-first/send-second completion and handle
  same-peer neighbors and zero-width slabs. Zero-width ranks use backend-safe
  dummy allocation and participate in complete FFTW-MPI solves.
- PM and short-range TreePM now use the physical Newton constant converted from
  frozen units and return a scale-free comoving kernel. The former extra `a^2`
  factor is removed; collisionless KDK and the gas conservative source apply
  the cosmological response to their respective state.
- The gas source now treats gravity lanes as scale-free `A`, applies `A/a^2`
  to momentum and gravity work, and uses
  `S_E=rho(u dot A)/a^2-H(2K+3P)`. Fixed-grid and AMR callbacks consume the
  post-drift `timeline_step.scale_factor_end/hubble_end_code`, not the
  still-step-begin committed integrator epoch or an independently recomputed SI
  Hubble rate.
- Gravity timestep proposals now keep coordinates consistent: after step
  commit, both particles and gas cells call the new public
  `computeComovingGravityTimeStep` helper with comoving softening, scale-free
  `|A|`, and committed `a`. It validates finite positive `a`, converts to
  comoving-coordinate acceleration `|A|/a^3`, and uses
  `dt_grav=eta sqrt(a^3 epsilon_com/|A|)`. The older generic helper remains
  coordinate-neutral.
- Production cadence is restricted to one and `hierarchical_max_rung` to zero.
  Every rank-coordinated production force-refresh surface rebuilds PM;
  mixed-rung execution is blocked until per-element kick/drift epochs exist.
- The coordinator now records a transient long-range cache signature containing
  force epoch, field-build scale, `G_code`, split and box geometry, assignment,
  boundary, PM decomposition mode, and deconvolution. Explicit reuse with any
  incompatible or missing rank-local cache fails coherently; divergent votes
  fail before PM collectives. Particle decomposition epoch is intentionally not
  a fixed-slab PM invalidator.
- The workflow decomposition epoch now advances only after an actual globally
  committed ownership change, is persisted/restored, and guards tree wire
  records. Dense-row force history is invalidated immediately after such a
  change so checkpoint truth cannot contain stale acceleration rows.
- Gravity observability now includes local source/target/tree counts, global
  empty-source/target rank counts, remote hierarchy packets, unique peers, PM
  solve/reuse and halo/slab dimensions, plus tree build, multipole-refresh, and
  opened-node counters.
- `gravity.pm_long_range_field` now serializes the in-flight
  `PmRefreshDirective` decision's opportunity, field version, build step, and
  build scale factor before cadence commit. It no longer reads the previous
  committed `pm_sync_state` for a current refresh event; this is an
  observability correction, not a cadence-owner or restart-schema change.
- Runtime gravity diagnostic doubles now use scientific `max_digits10`
  formatting for `G_code`, relative-MAC floors, split/cutoff/softening scales,
  force L2 norms, and scale factors. This prevents small nonzero evidence from
  appearing as `0.000000`; solver state and output schemas are unchanged.
- The independent rectangular Ewald reference certifies the FFTW-backed
  TSC+deconvolution profile. CIC remains available but is explicitly diagnostic
  for this accuracy envelope. The normalized production defaults are now
  TSC+deconvolution and `rcut_cells=6.25` (`r_cut/r_s=5` at the default
  `asmth=1.25`). The former 4.5-cell cutoff remains accepted but diagnostic
  after showing roughly nine-percent cutoff-transition error. Unchanged configs
  that relied on implicit defaults therefore have a reproducibility-visible
  force-profile and short-range-work change.
- The certified profile obtains the public/runtime `opening_theta=0.7` default
  from typed configuration. An independent-target matched-MAC fixture asserts
  geometric, COM-distance, and relative-force error/work behavior against Ewald
  rather than only printing diagnostics.
- Periodic TreePM now fails closed unless
  `r_cut < min(Lx,Ly,Lz)/2`, at both typed config loading and direct coordinator
  entry. The residual evaluates one minimum image per source, so coarse decks
  must increase PM resolution or reduce `rcut_cells`; the DMO gate uses a
  `16^3` mesh so its default `6.25`-cell cutoff remains strictly inside the
  half-box bound.
- Isolated/open assignment and gather stencils clip at the physical domain
  rather than wrapping. Typed config now rejects
  `treepm_enable_window_deconvolution=true` after mode policy resolves to
  isolated/open gravity; tracked isolated decks explicitly select `false`, and
  runtime does not silently override the frozen policy. Focused zoom correction
  also explicitly disables unsupported isolated window deconvolution.
- A production-workflow 64-particle Zel'dovich validation exercises config
  loading, KDK, periodic TreePM, HDF5 restart, stable ownership, np1--np4, mode
  growth, mass/momentum/COM invariants, and exact direct-vs-resumed state.
- The follow-up DMO spectrum gate adds a deterministic binned 3D `P(k)`
  artifact: `12^3` CIC with matched deconvolution, 32 requested bins, all 1727
  non-DC modes, four explicit empty bins, and reported-not-subtracted Poisson
  level `0.015625`. Initial/evolved fundamental powers are
  `2.1315894140899176e-6` / `2.1827523417862502e-6`, giving growth
  `1.011929959672577` versus `1.0119399775244009`; normalized np1--np4 JSON
  SHA-256 is
  `53ee518b63d21cd28790415b8f076614497536bb6a62c95276b8ed6c14984bf3`.
- The final DMO guards compare direct displacement-mode and square-root-power growth
  increments against `a_final/a_initial-1` within 7.5 percent and require the
  two increments to agree within relative `1e-3`. The zero-force Hubble-drag
  comparison additionally gates response/ballistic RMS and growing-mode
  projection in `[0.01,0.03]`, with alignment at least `0.99`. The exact
  initial `4^3` cancellation-dominated lattice is separately force-certified
  against periodic Ewald and predicts response near `0.019`; this is an
  independently anchored integration regression, not an analytic LCDM proof.
- A cross-preset non-MPI serial/MPI1 comparison covers 64 stable IDs with zero
  position and mass error, maximum velocity error
  `0`, and byte-identical power JSON SHA-256
  `8acbf3d6826250ca4fbabb0761511ff1178cbb7c20472cb2d8c8073dd16d355c`.

Command-backed focused evidence:

- PM ownership/halo, periodic TreePM, and strict Phase-2 np1--np4 matrix:
  12/12 passed in `build/mpi-hdf5-fftw-debug`.
- Ewald plus periodic TreePM serial/MPI np2--np4 matrix: 6/6 passed.
- The follow-up focused periodic/workflow regression matrix after cache and
  diagnostics hardening also passed 6/6.
- `TMPDIR=/tmp ctest --test-dir build/mpi-hdf5-fftw-debug -R
  '^validation_dmo_zeldovich_workflow' --output-on-failure`: 5/5 passed for
  MPI1/np2/np3/np4 producers plus rank comparator.
- Certified worst TSC metrics are now set by the exact cancellation-dominated
  `4^3` DMO initial lattice: `relative_L2=8.059667017e-3`,
  `p99=8.157551429e-3`, net-force fraction `5.985844233e-17`, and translation
  drift `1.523007584e-13`. Ewald/PM/tree/total force L2 is
  `10.16780756` / `8.296819318` / `1.789203432` / `10.08586460`. The separate
  rectangular seam fixture retains translation drift `3.031142721e-14`.
- DMO direct-mode amplitude grew from `0.0049999895833398901` to
  `0.0050596390031745634`, or `1.0119299088208957` versus expected
  `1.0119399775244009`. Maximum phase drift was
  `7.2112604976761041e-14` radians with coherence one; restart phase/state drift
  was zero. Response/ballistic RMS was `0.018731640935373602`, projection
  `0.01873163714589127`, and alignment `0.99999979769618985`.

Limits deliberately remain open: the distributed tree is a bounded hierarchy
all-gather plus target/response all-to-all, not LET; PM uses communicator-wide
record collectives; evidence stops at four ranks; zoom and distributed open PM
retain bounded all-gather/root-workspace paths; and GPU TSC parity is absent.
Production hierarchical KDK is intentionally rung-zero only.
The controlled DMO gate now emits the required versioned binned 3D
power-spectrum artifact, but its direct-DFT `12^3`, 64-particle fixture is not a
large-volume resolution/convergence, performance, halo-statistics, or
cross-code campaign. Publishable DMO production and AREPO/GADGET-class parity
are not claimed.


## 2026-07-17 runtime capability boundary closure

- Removed public/protected recovery of `core::StepContext` from all typed workflow
  stage views and owner interfaces.
- Added task-scoped descriptor grants with typed-view access allowlists and
  source-private built-in access validation.
- Added expected compile-failure coverage for the former analysis-owner and
  direct-view context escapes.
- Consolidated gas-cell identity, parent-row, owner-rank, and compatibility-mirror
  helpers under `workflows/internal/gas_cell_ownership.*`.
- Removed Windows ADS metadata and generated validation outputs; repository
  hygiene now detects those recursively.
- CPU focused workflow/resource/registry tests pass. MPI/FFTW/HDF5 runtime
  validation remains delegated to a dependency-complete environment.

## Campaign B acceptance repair (2026-07-24)

Source-level blockers from the adversarial IC-ingestion audit were repaired:
complete typed no-guess bridge configuration, manifest schema v3, one shared C++
conversion engine and compiled converter, pre-read HDF5 shape/type validation,
full field/family dispositions, distributed per-file hash/schema assignment,
collective phase guards, exact file/chunk and source/final reconciliation,
all-axis/sidecar validation, disk-backed bounded-memory global duplicate-ID
audit, stable tracer-parent remapping after routing, truthful I/O and staging
counters, build-aware provisional capability reporting, expanded direct,
canonical, manifest, conversion, type-policy, scaling, route-mutation, and fault
MPI tests, and clean packaging rules. CPU/HDF5/sanitizer evidence is recorded by
the repair report. MPI runtime acceptance remains external because `MPI_CXX` is unavailable
in the repair environment; this is an evidence blocker, not a source-level
waiver.

## 2026-07-25 Campaign B scientific-correctness and streaming completion pass

- Accepted real AREPO signed `NumPart_ThisFile` and `NumFilesPerSnapshot`
  attributes through checked wide intermediates while preserving datatype
  signedness/width/rank/extent in manifest schema v4; negative, overflow,
  floating, and malformed attributes fail closed.
- Replaced path-forced velocity powers with field-specific contracts:
  `Velocities` alone uses the selected snapshot velocity convention;
  `InternalEnergy` and `BH_Mdot` use independent `L^2/T^2` and `M/T`
  conversions. Independent golden tests cover nontrivial `a`, `h`, and units.
- Added typed missing-field policy/config/provenance contracts. Gas
  thermodynamics, stellar sidecars, and BH accretion reject by default;
  configured or dialect defaults are explicit and audited. Unsupported AREPO
  wind birth times and canonical tracer mapping fail early.
- Canonical schema version 2 embeds normalized manifest JSON, recomputes and
  verifies SHA-256 against the sidecar, requires a completion marker, and records
  verified status. Converter finalization is a rollback-capable HDF5 + manifest
  + marker bundle transaction with fault tests.
- Replaced converter full-state materialization with direct
  source-chunk-to-canonical hyperslab appends and bounded external duplicate-ID
  merge. Converter evidence records peak batch/chunk capacity and explicitly
  states that no complete `SimulationState` was materialized.
- Added persistent payload reader sessions with source size/hash revalidation and
  cached dataset handles. Distributed routing now groups multiple chunks into a
  configured bounded batch and exchanges once per batch.
- Split the field-conversion catalog, canonical-bundle verifier, persistent
  reader session, byte codec, endian-stable record codec, and SHA-256 into
  narrow source units; the record codec has an independent test. The remaining
  `ic_reader_file_set.cpp` is still large and is tracked as a maintainability
  follow-up rather than misreported as fully split.
- HDF5 focused unit/config/capability/converter tests pass. MPI configuration is
  blocked by missing `MPI_CXX`; distributed capability remains provisional.

## 2026-07-26 Campaign B final acceptance source closure

- Made supplied schema-v4 manifests the sole runtime scientific authority for
  source order, hashes, aliases/conversions, species policy, and normalized
  missing-field policy/value contracts. Manifest startup no longer requires
  `mode.ic_file`; an optional path is only a conflict check.
- Closed the remaining pre-collective exception seams for manifest digest,
  serialization and exchange accounting, exact source/final reconciliation,
  duplicate-audit staging, reader-session creation, and completion identity
  validation. Added finite-time rank-selective fault fixtures for each seam.
- Replaced batch-index reader rotation with deterministic per-file reader
  ownership and one persistent payload session. Full-file hashing is bounded at
  three passes per source file independent of batch count.
- Strengthened source-session identity with size/time and native device/inode
  checks where available, plus SHA-256 validation at session start and end.
  Same-descriptor hashing is not claimed.
- Added observable reader imbalance, routing/main/exact-audit exchange counts,
  logical collective phases, bytes, bounded staging, and ingestion wall time.
- Made canonical single-file output fail before HDF5 creation when any
  `NumPart_ThisFile` family count exceeds `UINT32_MAX`.
- Removed the 4,136-line IC monolith. The largest replacement source is below
  2,000 lines, and structural CI prevents a replacement monolith or private
  header leakage.
- Removed the two Windows ADS sidecars and added recursive source-only packaging.
- Focused HDF5 production-path, reader, capability, and structure tests pass.
  MPI-enabled source/test syntax passes, but runtime np1/np2/np4 acceptance is
  environment-blocked because `MPI_CXX` is unavailable; distributed capability
  therefore remains provisional.


## 2026-07-26 Campaign B distributed final-gate micro-patch

- Every supplied manifest pointer is now validated or rejected at serial,
  streaming, and distributed entry; the old non-empty-hash authority heuristic
  is removed.
- Distributed source inspection derives its dialect from the validated supplied
  manifest. Supplied canonical manifests execute embedded/sidecar/digest/marker
  verification, and bridge/canonical mismatches fail explicitly.
- Exact source/final and duplicate-ID audit capacity, path-storage, and peak
  accounting are rank-consistent guarded phases with finite-time fault seams.
- A central IC MPI wrapper counts actual `Allreduce`, `Bcast`, `Gather`,
  `Gatherv`, `Alltoall`, and `Alltoallv` calls with saturating `noexcept`
  increments. Logical phases remain separately named; routing protocol v1 costs
  exactly 30 actual collective calls per successful batch. The non-routing
  identity is versioned as `40 + runtime_cosmology_vote + source_file_count +
  10 * distributed_id_audit_round_count + mpi_bcast_call_count`.
- The distributed regression inventory now includes supplied canonical manifests,
  exact direct-versus-supplied canonical state equivalence, canonical
  tampering/missing-marker/bad-digest cases, source-path/hash mismatches,
  dialect mismatches, invalid manifests, per-counter rank agreement and exact
  cost formulas, and the remaining audit-accounting faults.
- Distributed runtime acceptance remains provisional where MPI execution is
  unavailable; serial/HDF5 and clean-package evidence is recorded in the final
  execution plan and handoff report.

## 2026-07-30 Campaign B source-package completeness closure

- Reproduced the final-gates audit against the submitted archive: two encoded
  ADS sidecars made actual-tree hygiene fail, while the packaging script's
  unanchored `core`/`core.*` exclusions removed `src/core/` and
  `include/cosmosim/core/`, yielding a clean-looking but unconfigurable archive.
- Deleted the actual ADS sidecars and replaced ambiguous crash-dump exclusion
  with fail-closed checks restricted to regular files. Directories named
  `core` are legitimate and must survive staging and extraction.
- Added shared package sentinels, regular-file SHA-256 manifests, symlink-target
  manifests, staged/extracted inventory equivalence, extracted-tree hygiene,
  fresh CPU CMake configuration, and required `cosmosim_core`/
  `cosmosim_harness` builds to the real packaging script.
- Added the release-labelled `integration_source_package_completeness`
  regression. It also proves a regular `core.1234` candidate is rejected and no
  failed archive is retained.
- Validation: CPU 133/133; HDF5 137/137; focused Campaign B HDF5 7/7. The totals
  are one higher than the prior baselines because of the new package test.
- MPI configuration remains environment-blocked because `mpi-cxx`/`MPI_CXX` is
  unavailable. Five distributed-IC translation units pass strict
  declaration-only syntax compilation, but no rank test is claimed.
- Campaign B source and archive-completeness acceptance are closed. Distributed
  runtime capability remains `provisional` until the documented np1/np2/np4
  and fault matrices execute successfully.

## 2026-08-08 `src/core` adversarial hardening campaign

- Replaced the relocatable vector-backed scratch arena with a reusable segmented monotonic arena. Allocations are pointer-stable for the scratch epoch, over-aligned requests use absolute-address alignment, and allocation/growth arithmetic is checked before wraparound.
- Strengthened cell/AMR ownership to an exact non-overlapping partition when patch ownership exists, retained an explicit non-AMR zero-patch mode, and standardized `kInvalidGasCellRow` as the hostless black-hole/tracer sentinel instead of relying on empty-container loopholes.
- Corrected periodic ghost indexing for widths larger than the interior and made wide reflective extension a deliberate repeated reflection that preserves the established first-domain ordering. `.param.txt` magnitude/unit values now reject trailing semantic tokens.
- Hardened KDK and hierarchical scheduler failure bookkeeping: tick overflow is preflighted before transition mutation, failed stages cannot remain marked inside KDK or restart-safe, and physical callback failure is explicitly treated as an abort/recover-from-checkpoint boundary rather than fictitious full-state rollback. Gas/AMR migration host remaps are prepared and validated before commit.
- Timestep hooks now fail closed on NaN, negative, zero, or `-inf`; only explicit `+inf` means no constraint. Active-view scatter prevalidates lane consistency, bounds, and duplicate destinations, while the all-particle gravity path aliases canonical hot SoA lanes instead of repacking the full state.
- Profiling output now uses complete JSON control escaping, JSON-safe non-finite handling, deterministic ordering, CSV quoting, checked `.part` publication, non-throwing RAII scope destruction, and the fixed multi-use scope macro. Memory accounting covers the audited omitted core lanes and uses checked accumulation. Metadata/provenance decoding is strict about malformed/duplicate/unknown values, and root provenance publication is staged through `.part`.
- Unit/cosmology helpers reject non-finite or non-positive conversion/scale inputs and adaptive cosmology integration reports tolerance non-convergence rather than silently returning a depth-limited estimate. Per-rank particle/cell/patch capacities are explicitly bounded by the retained 32-bit local index contract.
- MPI topology now records node-local rank for CUDA selection and the harness uses distributed fatal-abort semantics instead of allowing one rank to return into later collectives. CUDA discovery preserves the runtime diagnostic. MPI/CUDA execution was not available in the repair environment, so those runtime paths remain evidence-blocked rather than claimed green.
- Added a CI `cosmosim_core` strict-warning gate using `-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Wundef -Wformat=2 -Wnull-dereference -Werror`. The same CPU-only strict core build passes in the repair environment.

Remaining intentional boundaries: active-bin emission still sorts for deterministic order pending a deeper ordered-bin data-structure redesign; migration retains old/new authoritative buffers during prepare/commit to preserve the strong failure guarantee, so peak-RSS optimization needs measured evidence; provenance v6 retains the established 64-bit compatibility hash and limited hardware schema until a separately versioned strong-integrity/hardware-metadata migration is designed.


## 2026-08-09 `src/core` final audit-lineage closure campaign

- Added a strict reusable core-internal line codec. `StateMetadata` now writes schema v3 with reversible escaping while reading valid v2, and provenance v6/v7 decoding rejects unknown or dangling escapes.
- Introduced `provenance_v7` without mutating v6: current records carry a real SHA-256 digest of normalized configuration text plus compiler flags, CPU/RAM/host and MPI/CUDA topology descriptors where available. Valid v6 records remain readable.
- Centralized the retained 32-bit local-index contract across authoritative resize, migration, AMR and scheduler narrowing paths; `UINT32_MAX` remains the invalid gas-cell-row sentinel.
- Replaced per-active-set comparison sorting with deterministic reusable radix ordering. Release evidence shows material `beginSubstep` improvement on 100k/500k and fine-bin-heavy scheduler cases while preserving numerical order.
- Removed accidental O(N^2) gas-cell migration validation by validating dense identity once at the operation boundary, and removed a redundant full migration-field staging copy. 50k-cell peak/steady RSS measurements are ~1.91x mostly-local, ~2.04x moderate, and ~2.31x heavy while retaining prepare/validate/commit.
- Added process-memory calibration to memory-accounting benchmarks. On the recorded Linux 200k-particle/100k-cell case, known-owned capacity is within ~0.53% of the measured allocation-induced RSS delta; RSS remains explicitly diagnostic rather than treated as exact ownership truth.
- MPI/CUDA runtime execution remains environment-blocked where the relevant tools/devices are absent; no runtime promotion is claimed from source evidence alone.

## 2026-08-29 Campaign M1B memory-governor foundation

- Added one process-level `core::MemoryGovernor` owned by the reference workflow
  composition root and borrowed through `RuntimeServices`; no global memory
  singleton or `core -> workflows` dependency was introduced.
- Reused `parallel.process_memory_budget_bytes` plus the existing process
  external-reserve, output-overlap, and safety-margin policy. Budget zero keeps
  unlimited semantics while accounting and high-water telemetry remain active.
- Added move-only RAII `MemoryReservation` with distinct pending/committed
  accounting, checked 64-bit arithmetic, deterministic pressure/headroom,
  rejection diagnostics, class attribution, and control-plane mutex safety.
- Extended `MemoryEntry` with orthogonal `MemoryClass` plus a governed-commitment
  reconciliation marker. `MemoryReport` remains the ownership tree; governed
  capacity is excluded only from the governor baseline so the same physical byte
  range is not counted twice in admission arithmetic.
- Governed real `MonotonicScratchAllocator` block growth reserve-before-allocate.
  Reservations remain committed across `reset()` because physical blocks remain
  resident, and are released with block destruction.
- Closed the post-campaign production-use gap by admitting the rung-zero
  all-particle gravity `uint32` index lane through governed scratch after the
  existing collective DMO process preflight. `FailureCoordinator` now makes a
  rank-local scratch admission/allocation failure collective-safe before later
  TreePM communication.
- Represented particle and gas-cell scheduler retained state as explicit
  persistent `MemoryEntry` owners, so the governor baseline is fully derived
  from the merged ownership report instead of manually adding scheduler bytes.
- Made finite-budget pressure thresholds total for 1-3 byte diagnostic/test
  ceilings while preserving the 85%/95% policy at realistic process budgets.
- Added governor fields to the existing profiler memory-report JSON and human
  memory formatting. `RuntimeResourceLease` remains unchanged as the independent
  stage freshness/capability contract.
- Focused governor/scratch/accounting/profiler tests and the production reference
  workflow smoke cover finite-budget propagation, telemetry, lifecycle,
  exception/move semantics, concurrency, retained scratch reset/reuse, and the
  existing impossible-budget preflight rejection.
- Repaired the two AMR hydro integration regressions without loosening
  conservation tolerances: the fixtures evolve with `gamma=1.4`, and refine /
  derefine now receive the same explicit `ProductionAmrHydroOptions` instead of
  silently using the API default `gamma=5/3` during topology changes.
- Broad PM/tree/FFT/migration/I/O/hydro/analysis reservation coverage and RSS/PSS
  reconciliation remain explicit M1C handoff items.

## 2026-08-30 MPI P1 count and collective-safety closure

- Added one lower-layer checked bounded-round planner for classic MPI
  counts/displacements. Logical totals/prefixes remain wide while every payload
  round is representable; a conservative internal transport ceiling keeps
  communication scratch bounded.
- Repaired particle and AMR-patch migration, global decomposition gathering,
  hydro boundary advertisement gathering, TreePM zoom field gathering, and the
  sharded source-runtime exchange so logical traffic is not limited by one
  `MPI_*v` `int` displacement domain.
- Repaired `MpiContext::gatherBytesToRoot` and
  `broadcastBytesFromRoot` so receive allocation/preparation is globally agreed
  before payload communication.
- Added lower-layer collective preparation coordination and pairwise ghost
  readiness handshakes. Rank-local allocation/validation failure after metadata
  now suppresses the matching payload phase coherently; sequential ghost
  schedules continue servicing later peer handshakes before raising the
  distributed failure.
- The bounded adjacent review moved compact rebalance cut samples, exact
  ownership hash redistribution, tree-pseudo wire gathering, hydro conservative
  flux corrections, and hydro ghost request/payload gathers onto the same
  bounded/coordination rules.
- Added synthetic-limit planner tests and MPI integration regressions for
  bounded rounds, one-rank post-control rejection, and later-peer liveness after
  a pairwise preparation failure.
- No physics, force law, hydro method, ownership authority, configuration
  schema, HDF5 schema, or restart schema was changed by this repair. Real
  multi-rank execution remains a separate validation requirement when an MPI
  development/runtime environment is available.

## 2026-08-31 MPI P1 post-campaign sparse-participation closure

- Closed the post-campaign sparse ghost-refresh divergence: an MPI-enabled rank
  with zero payload neighbors no longer returns before the final world-wide
  protocol failure agreement. Empty ranks perform no payload `Sendrecv`, but
  remain MPI-world participants through protocol completion.
- Added three-rank sparse-topology regressions with ranks 0/1 exchanging ghost
  payloads and rank 2 carrying no payload neighbor. The suite covers both normal
  completion and deterministic rank-0 post-metadata preparation rejection, with
  all ranks observing coordinated failure and reaching a post-test barrier.
- Extended the production distributed gas-cell migration integration test to
  force an 8-byte test transport ceiling, ensuring the real bounded Alltoallv
  path executes multiple rounds while existing restart, identity, count, and
  metal-mass invariants remain unchanged.
- Added the missing `PROCESSORS 2` CTest scheduling metadata for the two-rank
  collective-safety regression.
- Migration full-wire residency remains part of the existing M1C memory-governor
  integration handoff; this closure does not claim packet-streaming migration.
  The pre-existing sparse TreePM classic-count boundedness handoff is likewise
  unchanged.

## 2026-09-03 Campaign M1C-2 M1 acceptance closure

- Reconstructed the M1A -> M1B -> M1C-1 architecture from current source after
  all later authoritative MPI and M1C-1 post-campaign repairs.
- Closed the 12-gate M1 source/architecture acceptance matrix in
  `docs/repair/m1c2_m1_acceptance_closure_20260903.md` without changing the
  accepted physics, numerical tolerances, restart schema, or configuration
  authority.
- Re-audited the PM 128 MiB gate using retained-capacity arithmetic rather than
  logical payload. The conservative 8-rank force/potential interpolation model
  is 134,151,808 B = 448 B rank metadata + two 67,075,680 B directional wire
  capacities, leaving 65,920 B beneath 128 MiB. The model is independent of
  `N_local`, and runtime checks actual retained wire capacities plus metadata.
- Strengthened `unit_pm_solver` so the certified 512^3 / 8-rank local population
  and exact 128 MiB expression are explicit regression evidence rather than an
  undocumented arithmetic assumption.
- Focused CPU M1 acceptance tests pass. Dependency-enabled distributed runtime
  status is recorded in the M1C-2 acceptance document and must remain
  environment-blocked where MPI/FFTW development support is unavailable.

## 2026-09-03 P0/P1 post-adversarial-audit runtime closure

- Repaired the rank-local versus communicator-global TreePM/PM state boundary:
  dense particle/cell index generations no longer advance global gravity-source
  generation, while PM collective fingerprints retain global solve semantics
  without incorrectly requiring rank-local target-layout enums to match.
- Preserved optional-lane absence across generic ghost wire/unpack/commit so DMO
  ghost state cannot silently acquire hydro lanes; hydro ghost truth remains on
  the `gas_cell_id`-keyed path.
- Enforced the gas optional-parent invariant at `GasCellIdentityMap`: absent
  means no parent, present means a valid nonzero parent ID. HDF5/restart and
  star-formation/AMR roundtrip coverage passes where executable.
- Corrected the periodic TreePM integration contract to rely on the independent
  periodic Ewald validation rather than an invalid minimum-image PM-only
  monotonicity assertion; no force tolerance or certified split profile changed.
- Closed narrow distributed-IC acceptance issues in collective telemetry,
  malformed-manifest fixture setup, and single-rank `already_partitioned`
  semantics without changing bounded routing/reconciliation/failure
  coordination.
- CPU and HDF5 executable matrices are green outside the long source-package
  completeness test, and sanitizer coverage shows no production corruption;
  the deliberate `SIZE_MAX` allocation test remains an ASan-runtime exception.
- Full P0/P1 runtime closure remains environment-blocked because FFTW3 and MPI
  C++ development support are unavailable. Do not promote distributed IC or
  multi-rank TreePM/DMO acceptance until the registered np2/np3/np4/np8 matrix
  passes on a dependency-complete system. Detailed evidence is in
  `docs/repair/p0_p1_post_audit_runtime_closure_20260903.md`.

## 2026-09-05 Campaign M2C full-physics sidecar/event-memory closure

- Preserved the existing species-local gas/star/BH/tracer sidecar authority and
  added acceptance coverage proving disabled species retain zero cold-sidecar
  capacity in DMO/gas-only states.
- Removed the 48 B/star duplicate stellar-feedback carry mirror; persistent
  unresolved budget remains solely in restart/migration-authoritative star
  sidecar lanes.
- Replaced 24 B/star source-step feedback delta staging with bounded 4096-event
  batches and replaced per-event O(N_gas) distance materialization with an exact
  nearest-k lookup over a 4 B/owned-leaf x-sorted spatial index.
- Removed feedback-only full gas volume/ownership staging in favor of on-demand
  AMR volume derivation and removed tracer all-cell retained index/mask scaling.
- Reference hydro, sources, and adaptive timestep logic now share one immutable
  effective-multiphase EOS table rather than constructing three copies.
- Existing typed cooling/metal configuration remains the scientific fidelity
  authority. `core_elements` still fails closed until its full abundance,
  transport, cooling, restart, and MPI contract exists.
- Detailed source/validation evidence is recorded in
  `docs/repair/m2c_full_physics_sidecars_20260905.md`.

## 2026-09-05 M2C post-merge transient-memory closure

- Closed the post-merge audit finding that star formation could retain/materialize
  population-scale all-active indices and `StarFormationCellInput` rows. The
  reference `SourceRuntime` now scans in globally coordinated <=4096-cell
  batches and governor-admits input-capacity growth.
- Removed pre-eligibility `StarBirthPlan` population reservation while preserving
  birth-key, RNG, conservation, and particle-ID precommit semantics; empty local
  MPI batch rounds still participate in collective ID preparation.
- Replaced the production 128 B/gas-cell metal-diffusion mirror plus full-cell
  volume staging with canonical spans, 49 B/cell compact phase state, reusable
  five-scalar workspace, and cell-local topology-gradient reduction. Conservative
  equal-and-opposite flux ordering is regression-tested against the legacy
  operator.
- Removed the O(N_BH) all-active accretion index vector; empty active-index span
  denotes all local BH rows. Production seeding remains fail-closed without an
  authoritative seed-candidate provider.
- `SourceRuntime` capacities/high-waters are now part of central runtime memory
  truth and persisted profiler JSON. A 4097-cell regression proves SF input
  residency plateaus at one 4096-cell batch; compact diffusion repeated-step
  capacity also plateaus.

## 2026-09-06 M2C final acceptance closure

- Repaired the CHUI-native HDF5 bounds check so Mpc-valued
  `CHUIBoxSize{X,Y,Z}_MpcComoving` metadata is converted to the configured
  code/stored coordinate unit before validating decoded/stored coordinates.
  Snapshot schema names and report units remain unchanged.
- Added a dedicated governor reservation for retained metal-diffusion face
  topology. Replacement admission accounts for old/new graph coexistence before
  releasing the prior retained reservation.
- Added historical retained-capacity high-water tracking for move-replaced
  diffusion `rho*kappa` and face containers; smaller later topologies cannot
  erase prior peak evidence.
- Added focused regressions for kpc-vs-Mpc snapshot bounds, continued invalid
  out-of-box rejection, governed diffusion face residency, and high-water >=
  current retained capacity.
- Reproducibility/science impact: no source equations, diffusion operator,
  star-formation/feedback model, precision, HDF5 schema version, dataset names,
  normalized config, or restart payload changed.

## 2026-09-06 Campaign M2D memory-aware runtime scheduling

- Extended the existing process `MemoryGovernor` with checked deterministic
  headroom-aware batch sizing; hydro, remote hydro communication, and
  stellar-evolution/feedback paths now derive bounded residency from the same
  authority instead of using independent fixed-memory assumptions.
- Added major-task dependency, declared-peak/lifetime, optionality, and coarse
  compute/DRAM/MPI pressure metadata to the existing `RuntimeModuleRegistry`.
  The reference serial stage order is unchanged; the conservative overlap
  predicate establishes legality for future concurrency without changing the
  numerical schedule.
- Expanded SFC load metrics with transient-memory, source-event, and
  communication costs while preserving the existing geometry and weight
  authority. A hard per-rank peak-memory envelope now cuts/rejects partitions
  independently of the soft memory weight, and target candidates that exceed
  the cap are vetoed before migration.
- A current hard rank-memory violation forces rebalance despite non-triggering
  soft thresholds or migration throttling. The workflow derives one
  collective-safe cap from the minimum process-governor headroom and reserves a
  conservative maximum major-task transient allowance.
- MPI topology truth now includes node-local size as well as rank; workflow
  profiling records world/node rank topology and OpenMP thread configuration.
  No cross-rank shared-memory rewrite or new memory authority was introduced.
- Optional light/heavy science diagnostics now retain bounded pending state when
  Red/Trip pressure defers a due cadence and catch up at the first later
  non-Red analysis hook. An 8-step tiny workflow regression forces a step-4
  deferral and proves step-5 catch-up without bypassing the hard ceiling.
- Final available-environment validation is green for all 136 CPU tests other
  than the separately long-running source-package-completeness regression,
  which exceeded the external command window and produced no failure. The
  HDF5 configuration detects HDF5 1.14.5 and the focused M2D/reference-workflow
  set passes 4/4. MPI C++ and FFTW development dependencies remain unavailable.
  Detailed evidence is in `docs/repair/m2d_memory_aware_runtime_20260906.md`.

## 2026-09-06 M2D-1 source repair state

The post-merge feasible-SFC rejection and soft-throttle hard-memory relief defects are repaired in source. Exact cut repair has 86,016 exhaustive synthetic feasibility cases and focused CPU tests passing. Retained-capacity-aware batching, dispatcher-owned task leases, conservative unknown-peak overlap refusal, and coalesced optional-cadence lifecycle are implemented. Full acceptance is provisional pending real MPI+FFTW, owner-specific major-task peak contracts, and broad scientific/topology validation. See `docs/repair/m2d1_memory_scheduling_closure_20260906.md`. Do not treat historical M2D reports as fresh distributed runtime evidence.

## 2026-09-07 M2D final admission repair

The post-M2D-1 optional FFT admission finding is repaired in source with an
owner-held lease, checked workspace model and collective optional deferral.
Required-health preflight is connected to the existing task registry without
double charging. Complete major-owner peak contracts and dependency-complete
MPI/FFTW acceptance remain open. See
`docs/repair/m2d_final_admission_closure_20260907.md`. This entry is a repair
record, not an unconditional M2D acceptance declaration.

## 2026-09-07 M2D acceptance-gates source repair

The owner-contract patch exposes dynamic estimates with explicit completeness,
reuses the existing gravity phase and output staging models, publishes partial
source/AMR bounds, and preserves all physical owner reservations and serial
scientific ordering. The focused registry and reference-workflow gates pass.
The full-physics aggregate peak and distributed/runtime acceptance remain open;
see `docs/repair/m2d_acceptance_gates_20260907.md`. Historical M2D results are
not promoted to current executable evidence.

## 2026-09-07 — M2D coding closure continuation (provisional)

The current source contains the verified M2D acceptance-gates patch plus
precommitted retained-capacity growth for source births and additional
collective-safe SFC preparation boundaries. New tests cover exact headroom,
partial allocation failure/retry, scheduler/species-sidecar growth, and registered
MPI hard-memory/fault scenarios. The previous SFC/analysis repairs and serial
numerical execution are preserved. The source report
`docs/repair/m2d_coding_closure_20260907.md` records fresh commands and
limitations. Do not mark M2D fully accepted until remaining owner-local
allocation contracts and dependency-complete distributed/whole-process gates
are established; historical full-suite results are not fresh evidence.

## 2026-09-16 — DMO first-light correctness closure

The current source closes the DMO first-light audit leads for cosmological
normalization, standard GADGET/AREPO stored-velocity semantics, authoritative
particle drift-history restart equivalence, and invariant PM spectral-operator
cache reuse. The legacy custom `sqrt_a_scaled_peculiar` IC convention remains
unchanged, restart schema v23 is unchanged because the drift lanes were already
persisted, and PM scale factor remains part of collective epoch consensus even
though it no longer invalidates scale-free spectral arrays. Focused CPU and
HDF5 regressions pass in the available environment; FFTW/MPI PM validation
remains dependency-blocked. See
`docs/repair/dmo_first_light_correctness_20260916.md`.

## 2026-09-17 — DMO first-light runtime, IC, periodic-state, and snapshot-layout recovery

The first real DMO run exposed a config-driven timestep/endpoint authority bug,
not an FLRW bracketing defect. Rung-zero now selects the collective physical
minimum before KDK construction, `max_global_steps` is only a safety cap, and
scale-factor authority derives and clips to the FLRW code-time endpoint. The
repair also centralizes first-light periodic modulo semantics, accepts the
monofonIC HDF5 structural forms exercised by the failure bundle, normalizes
proven contiguous zero-based source IDs, and migrates new science output from
per-snapshot directories to stem-scoped sets inside one `snapshots/` directory
without changing science snapshot schema v6. The config-driven DMO smoke and
focused HDF5 IC/snapshot tests pass locally; MPI and FFTW runtime matrices remain
dependency-blocked.

## 2026-09-21 — 64^3 three-rank external-IC and science-output production closure

The successful 64^3/3-rank first-light campaign exposed remaining interoperability rather than solver defects. The source now treats zero-containing external GADGET/AREPO/MONOFONIC ParticleIDs as an order-independent external identity domain and applies one overflow-checked `+1` mapping before internal use. Distributed source/final species masses use compensated local accumulation and one shared MPI reduction lane order, retaining the existing tight tolerance and explicit loss/duplication failures.

Science snapshots now have a typed topology policy. `auto`/`single` selects one analysis-ready `snap_###.hdf5` for MPI only when MPI-enabled Parallel HDF5 is present; every rank writes owned global hyperslabs without gathering global state to rank zero. Publication is `.partial` -> distributed exact partition readback -> collective success agreement -> atomic final name -> completion marker. Root post-publication validation keeps bounded structural/scientific checks without materializing the redundant global ID set. Explicit `sharded` remains a compatibility backend; fixed-count `aggregated` is reserved and fails closed until implemented. Restart remains rank-oriented and unchanged. Normal completion guarantees one endpoint snapshot if cadence did not already commit that exact step/epoch.

The launcher now obtains `parallel.mpi_ranks_expected` from a narrow authoritative C++ preflight for known builds, so ordinary `./chui run CONFIG` can compose MPI without a second Python config parser. Configure probes actual C++20 `std::span` support, reports compiler/dependency providers, distinguishes serial from Parallel HDF5, and warns on suspicious mixed Conda/system stacks. Serial/HDF5 focused validation passes in the repair environment; real MPI/Parallel-HDF5/FFTW qualification remains environment-blocked where those development dependencies are absent.

## 2026-09-22 — Parallel-HDF5 / MPI diagnostics qualification closure

- HDF5 source/package dependency truth is the C API only; C is explicitly enabled and production Parallel-HDF5 capability is established by MPIO compile/link probes rather than `HDF5_IS_PARALLEL` alone. Production MPI+HDF5+FFTW presets require that capability; serial-HDF5 MPI remains explicitly compatibility-only.
- The mixed-prefix/Conda classifier normalizes separators and cannot abort configure through an invalid backslash regex. `Threads::Threads` is linked to the direct `std::thread` memory-governor test.
- MPI diagnostics now reduce communicator-global counts, boolean invariant failures, angular-momentum vectors, SFR bins, quicklook fields, and distributed memory summaries before one publisher writes shared files. MPI power spectrum is explicitly deferred instead of averaging rank-local spectra.
- Qualification is layered into narrow snapshot, restart, diagnostics, endpoint, real-launcher, feature-summary/runtime-truth, and installed-package-consumer gates. The real 64^3/3-rank first-light rerun remains a dependency/host acceptance gate and is not claimed by source inspection alone.

## 2026-09-22 — 64^3 Parallel-HDF5 smoke defect closure patch

Q-01 is repaired by removing premature MPIO-probe initialization and invalidating the cached check result immediately before the real `H5Pset_fapl_mpio`/`H5Pset_dxpl_mpio` compile-link probe. Q-02 now makes external-ID normalization typed source-set semantic state: source fragments report bounded zero/max observations, rank zero chooses one `IcExternalIdMapping`, the mapping is broadcast beside (not encoded in) the public manifest, chunk conversion and `ParentParticleIDs` consume the same mapping, and `IcImportReport` exposes the applied policy. Q-03 gives the single-file publication layer explicit ownership of creating/validating the shared snapshot directory before any rank enters Parallel-HDF5 creation; committed output protection and stale-partial policy remain unchanged. The public `IcManifest` schema and restart schema are unchanged.

Source-level and serial-HDF5 compilation can be exercised without MPI, but the real 3-rank original-MONOFONIC qualification remains a dependency/artifact gate whenever the host lacks MPI/Parallel-HDF5 or the unmodified zero-based IC is not supplied.

## 2026-09-23 — TreePM observability, exchange-memory model, and domain-geometry freshness

Code-first campaign (no build/test run in this session; validation intentionally deferred by user request). Scope is observability truth, a realistic LET exchange memory estimate under the existing MemoryGovernor, and routing-geometry freshness decoupled from ownership epoch.

### Changes applied

**P1 — residual counter truth**

- Split residual traversal counters into non-overlapping `local_owned_targets` and `incoming_remote_targets` bundles in `tree_pm_coupling.cpp`.
- Exact identity: `residual_pair_evaluations == local_pair_evaluations + incoming_remote_pair_evaluations`. Combined `tree_profile` visited/accepted/opened/PPI and cutoff counters remain the sum (compat).
- Workflow feedback uses residual only (fixes double-count of local+incoming as if both were remote).
- Event fields: `local_pair_evaluations`, `incoming_remote_pair_evaluations`, `remote_pair_evaluations` (compat alias of incoming), `total_pair_evaluations`.

**P1 — timer truth**

- `PmProfileEvent.total_ms` accumulated in the coordinator from PM-phase entry through tree short-range start; `pm_solver` totals add it once.
- Workflow `tree_wall_ms_recent` = `tree_short_range_ms` alone; `pm_wall_ms_recent` = `pm_profile.total_ms` alone.
- Remote-phase split: `incoming_remote_target_compute_ms`, `protocol_validation_ms`, `protocol_consensus_ms`, `response_exchange_ms`. `let_remote_traversal_ms` is a compat alias of incoming compute; `let_communication_ms` unchanged (still includes consensus).

**P1 — TreePM exchange memory model**

- Shared wire constants `kTreePmShortRangeRequestWireBytes=96`, `kTreePmShortRangeResponseWireBytes=80`.
- `estimateTreePmExchangeMemory` in `gravity_memory.cpp` replaces the prior `2*B` term with wire + structured/mask/count/accumulator/metadata + transient codec terms under checked arithmetic and `planSparseTreePmRound` clamp. Unordered-set allocator overhead is excluded and documented as uncertainty.
- Runtime high-water split: `let_wire_buffer_high_water_bytes` (four payloads) and `let_known_workspace_high_water_bytes` (full known peak); `let_high_water_bytes` remains a wire-only compat alias.

**P2 — domain-geometry freshness lifecycle**

- `TreePmOptions.authoritative_geometry_source_generation` + diagnostics freshness/fallback fields + `TreePmDomainGeometryFallbackReason` enum/name helper.
- `GravityRuntime::installAuthoritativeTopDomainLeaves(leaves, source_generation)` and `authoritativeDomainGeometryMatches`.
- Pre-solve O(N) `refitAuthoritativeTopDomainLeaves` (in `distributed_memory.cpp`, shared SFC helpers) after compact source rebuild; seed owner/epoch/SFC retained; empty leaves omitted; no non-finite published bounds.
- TreePM selects authoritative leaves only when generation matches **and** coverage passes; otherwise conservative tree-root fallback with a specific reason (`kNoGeometryInstalled`, `kStaleSourceGeneration`, `kDecompositionEpochMismatch`, `kSourceCoverageFailure`, `kGeometryPreparationFailure`).
- Seed path remains `MigrationBalanceRuntime::authoritativeTopDomainLeaves` at segment/restart/ownership-change only. Geometry is derived state: no restart/snapshot schema change. `decomposition_epoch` does not advance on drift.
- Event fields for geometry lifecycle plus OpenMP provenance from `core::openMpRuntimeInfo()`.

### Validation status

- Validation/testing intentionally deferred by user request.
- Build/test commands were not executed; `cmake`/`ctest` permission was denied in this environment.
- No production qualification claim is made by this campaign.
- Reproducibility: no config keys, no snapshot/restart schema, no solver numerics (theta, 0.08 gate, multipole order, cutoff, split, PM mesh, assignment, softening, cadence) were changed. Diagnostic/event payloads and a memory estimate model changed; residual pair identity and timer non-overlap are the intended observability fixes.

## 2026-09-23 — TreePM P1+P2 closure repair (seed geometry, timer purity, ABI asserts, ownership invalidation)

Small closure patch over the existing P1+P2 work in response to an adversarial static review. Implementation only; validation intentionally deferred by the user. No P1/P2 concept was reverted and no solver numerics changed.

### Changes applied

**1. Stable seed / current geometry separation**

- `GravityRuntime` now holds two sets: `m_authoritative_top_domain_seed_leaves` (stable decomposition-local leaf identities/SFC intervals, may include currently-empty groups, replaced only at segment/restart/ownership-change install) and `m_authoritative_top_domain_leaves` (current published bounds for the stamped source generation, may omit empty groups).
- `refreshAuthoritativeTopDomainGeometryForSolve` refits from the seed set, never from the previously published result, so repeated empty-group disappearance can no longer erode SFC partition identity across source generations. Legacy states without a recorded seed fall back to the current set.
- `installAuthoritativeTopDomainLeaves` replaces both sets and stamps the published generation. Out-of-seed-range sources remain covered by the nearest seed leaf with existing telemetry; fallback, coverage validation, and empty-leaf omission semantics are unchanged. No snapshot/restart schema change.

**2. Pure incoming remote compute timing**

- Split the incoming remote phase into non-overlapping scopes: new `incoming_request_decode_validation_ms` (decode, record-count/peer/epoch/identity validation, duplicate hashing, finite checks, response structure preparation), narrowed `incoming_remote_target_compute_ms` (validated target force evaluation against this rank's tree plus only the inseparable direct acceleration write), and new `incoming_response_encode_pack_ms` (response encoding, size check, pack into response send payload). `protocol_validation_ms` remains response count/displacement layout and payload buffer sizing; `protocol_consensus_ms` and `response_exchange_ms` unchanged. `let_remote_traversal_ms` remains a compatibility alias of incoming compute only. Balancing continues to consume residual pair counts and `tree_short_range_ms`, unchanged.

**3. ABI-safe packet/memory accounting**

- `static_assert(sizeof(ShortRangeTargetRequestPacket) == kTreePmShortRangeRequestWireBytes)` and the response equivalent enforce host/wire identity at compile time; a padding change fails the build instead of silently using wire bytes as host object bytes. Preflight remains rank-degree aware with the four-wire-buffer model `M_wire = 2*d*b*(96+80)`, `b = floor(B/96)`, `d = R-1`, plus structured/mask/count/accumulator/metadata/codec terms under `MemoryGovernor`. The transient codec term is documented as a conservative upper envelope, not an exact simultaneous peak.

**4. Ownership-change geometry invalidation**

- `commitParticleDecompositionChange()` now sets `m_authoritative_top_domain_source_generation_valid = false` when the decomposition epoch advances, and `publishAuthoritativeTopDomainLeavesSpan` publishes an empty routing view while invalid, so between commit and reinstall `authoritativeDomainGeometryMatches(...)` is false and stale leaves are not consumable. Ownership data, migration state, and the epoch are not cleared.

**5. Type attribute cleanup**

- `[[nodiscard, maybe_unused]] struct DomainCoverageResult` → plain `struct DomainCoverageResult` (the accessor keeps `[[nodiscard]]`; no warning suppression introduced).

**Optional retained-capacity telemetry**

- `let_known_workspace_high_water_bytes` now reads `vector.capacity()` for CHUÍ-owned structured request, mask, count, accumulator, and rank-metadata storage (including nested outer capacities) instead of reconstructing from `world_size`/`batch_size`; the transient codec term stays a modeled conservative upper envelope. Preflight remains a separate conservative estimate.

### Validation status

- Validation/testing intentionally deferred by user request.
- No production qualification claim is made by this closure campaign.
- Reproducibility: no config keys, no snapshot/restart schema, no solver numerics (theta, 0.08 gate, multipole order, cutoff, split scale, PM mesh, assignment/deconvolution, softening, cadence, timestep policy) were changed.

## 2026-09-23 — TreePM P3: OpenMP residual execution

Implementation-only campaign over the existing P1/P2 worktree. Validation intentionally deferred by the user. No P1/P2 concept was reverted and no solver numerics changed (same tree, same MAC/0.08 gate, cutoff, softening, per-target traversal order, MPI protocol).

### Changes applied

- Shared immutable residual tree across OpenMP workers; one contiguous `m_worker_stack_storage` of `T * S` slots with `S = 1 + 7*D`, with `D` bounded by the builder-enforced `kMaximumTreeDepth`; MemoryGovernor terms also cover O(worker) integer counters and O(ceil(A/64)) deterministic block diagnostics.
- Residual evaluation parallelized **between targets only** (`#pragma omp for schedule(dynamic, 1)` over 64-target logical blocks); no atomics in hot node/pair loops; integer counters are worker-scale and merged exactly after join; `local_short_range_sum_sq` remains block-deterministic.
- The full-active target-softening lane was removed. Main-thread validation checks spans, species ranges, finite/non-negative values, and independent-target sidecars; workers use the prevalidated allocation-free unchecked resolver and reuse the immutable resolved source lane. `MPI_THREAD_FUNNELED` remains the only MPI threading contract.
- Worker exceptions are captured inside the work loops under named critical sections and rethrown after join. Stack, counter, and block scratch is prepared before distributed request posting; failures are collectively coordinated.
- `openmp_observed_workers` records the maximum actual team size across local, distributed overlap, and incoming-target regions; serial execution reports one. The workflow emits configured/planned worker provenance and residual scratch high-water.
- Serial fallback under `COSMOSIM_HAVE_OPENMP=0` preserved; one shared kernel; no `omp_set_num_threads` in TreePM.
- Docs: `docs/tree_pm_coupling.md` ("## OpenMP residual execution (P3)") and `docs/profiling.md` OpenMP provenance/counters.

### Validation status

- Validation/testing intentionally deferred by user request.
- No production qualification or measured speedup claim is made by this campaign.
- Reproducibility: no config keys, no snapshot/restart schema, no solver numerics changed.

## 2026-09-23 — Parallel P4: compact large-N decomposition planner

Implementation-only campaign over the existing worktree (P1/P2/P3 retained). Validation intentionally deferred by the user.

### Changes applied

- Runtime rebalance now streams canonical particle/AMR-patch spans into the source-view overload of `buildCompactDistributedRuntimeRebalancePlan`; it does not call `buildRuntimeDecompositionItems` or materialize `vector<DecompositionItem>[N]`. The separate initial ownership placement adapter retains its specialized rich-item policy and is outside this runtime-rebalance closure.
- Feedback normalization is a direct canonical-state pass, followed by transient per-entity component recomputation while appending ≤64-byte `CompactRuntimeDecompositionRecord` values. The same weighted-load policy and coefficients remain in force.
- The compact planner retains one in-place `(sfc_key, entity_id, local_index)` sort, ≤256 cut samples/rank, bounded prefix memory queries, streaming metrics, and exact owner-change migration intents. Intent storage is reserved/admitted with a conservative particle-plus-patch bound.
- Top-domain seeds use a direct minimal sortable record path over canonical coordinates and streaming AMR bounds; the committed phase reservation spans construction and installation into the gravity owner’s seed/published vectors; `refitAuthoritativeTopDomainLeaves` and the P2 seed/published lifecycle are unchanged.
- Runtime telemetry reports local planner peak bytes, local entity count, compact/sample/prefix/intent/other-scratch terms, and local bytes-per-entity. The rich planner remains reference/debug/test-only.
- Docs: `docs/parallel_distributed_memory_contracts.md`.

### Validation status

- Validation/testing intentionally deferred by user request.
- No production qualification or measured speedup claim is made by this campaign.
- Reproducibility: no config keys, no snapshot/restart schema, no solver numerics (opening_theta, MAC, cutoff, split scale, softening, PM mesh/TSC/cadence) changed; production rebalance call path switched to the compact planner with identical cut/MPI protocol semantics.

## 2026-10-15 — TreePM performance campaign

### Changes applied

Tree-ordered TreePM residual traversal, optional threaded FFTW, adaptive
`relative_monopole` estimator and `com_distance` cutoff-node policy (opt-in),
hierarchical `rms_split_scale` coarse criterion (opt-in) with coarse-limiter
diagnostics, post-bootstrap governor reconcile (B-10 site), Open MPI hybrid
placement in `chui`, measurement tooling, fast continuation deck, and repairs
of stale tests that no longer compiled or ran.

### Validation status

Compiled and tested in this session (`cpu-only-debug` full CTest,
`pm-hdf5-fftw-debug` gravity/hierarchical gates, Release benchmark accuracy
matrix). Exact results: [campaign record](repair/treepm_performance_campaign_20261015.md).
