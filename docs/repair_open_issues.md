## Round 8 corrective source changes — 2026-10-29

20261103 CHUI Part 1 — source complete, NOT compiled/run. Closed in source:
F1 (Hancock predictor sign + 1/a), F3 (HLL passive lanes), F4 (post-feedback
thermodynamic refresh before CFL), F7 (ISM ledger volume), periodic fixed-grid
double boundary flux. Open: F2 cooling interval acceptance (Part 2), F5
newborn stellar age (Part 2), F6 automatic distributed AMR hierarchy (Part 3),
multi-layer AMR ghosts for full-order PPM/CTU at patch boundaries (Part 3),
**MANDATORY PART 3: cross-patch periodic AMR topology** (now refused
instead of zero-flux box faces; cosmo_cube on decomposed root grids depends
on it), transactional AMR sweep (Part 3). Corrective pass also fixed
AMR periodic zero-flux faces, AMR pressure-cache thermal authority and the AMR
temperature cache, and redesigned dual energy (no E overwrite). See the
[Part 1 record](repair/chui_part1_foundation_hydro_20261103.md).


- Collective science snapshots require one MPI-shared run root: the source-runtime star-formation test now supplies the same initial/resumed directory on every rank. The production path-consensus guard is unchanged.
- IC reader fixed non-routing collective count updated to 37 (Round 7 instrumentation: 66 actual = 37 fixed + 1 cosmology + 2 files + 20 audit + 6 Bcast); strict per-rank protocol assertions remain enabled.
- Owned grouped TreePM pair accounting now excludes exact same-source identity, while counting physically distinct coincident particles, matching the per-target walk. MPI and force-convergence qualification still required.
- Restart read now rebases pending AMR flux-register transient patch geometry epoch only **after** verifying the persisted payload integrity digest and stable gas/patch ownership. No schema change; MPI restart validation remains required.
- The strict TreePM historical normalized-defaults migration gate now accepts integrity-verified v24 as well as v23 checkpoints, without admitting older, unverified, or non-global-KDK checkpoints.
- Outstanding unverified/high-risk failures: MPI_ERR_TRUNCATE rank-divergent collective sequences; hierarchical TreePM trajectory accuracy and verified checkpoint migration; PM mesh-refinement accuracy; SFC hard-memory collective rejection; gas migration metal conservation. Do not waive their acceptance gates.

# Repair open issues

20261027 round 5 — source freeze. Engineering code-complete; NOT compiled,
NOT tested; owner validation required.

* Row-aligned transactional force restore (no O(N) ID hash map, no
  allocation, collective two-phase validation).
* Capacity-derived rollback memory governance (147 B/particle model,
  release-before-grow, reconciled reservation).
* Identity-row and max-mass scheduler caches. This fixed an unbounded
  per-block scratch-arena growth, and adds
  `SimulationState::particleMassGeneration()`.
* New tests: `unit_gravity_force_endpoint`, `unit_particle_mass_generation`,
  workflow sections (3), (7), (8), (10), (11).

See the [round-5 record](repair/source_freeze_round5_20261027.md).

20261026 round 4: final pre-validation completion. Engineering code-complete
at source level; NOT compiled, NOT tested; owner validation required.

* Transactional multirate rollback with an explicit restoration/invalidation
  contract and a state inventory.
* Governed reusable rollback store.
* Collective memory-pressure policy
  (`numerics.treepm_multirate_rollback_memory_policy`).
* Transactional accepted-trajectory diagnostics.
* Cosmological driver suites, plus rollback, fallback, failure-injection and
  skewed-rank workflow tests.

No known unresolved source-level defect. Limitations are classified in the
[round-4 record](repair/multirate_round4_completion_20261026.md), §E, which
also contains the owner validation manifest.

20261026 round 3 (source only; NOT compiled, NOT tested; owner validation
required):

* Explicit checked AREPO ParticleID translation (`--arepo-id-offset`).
* Multirate timestep invariants enforced: per-block criteria, transactional
  cycle rollback/retry.
* Soft stop at PM cycle ends for segmentation invariance.
* Exact cycle-end rounding.
* Strengthened multirate contract.

Open:

* Interior-cycle restart is unsupported (by design; no stop source can end inside a cycle).
* Convergence order is unmeasured: driver-level analytic suites (PM-only, tree-only, combined, variable, adaptive, two-body) and a genuine-rejection workflow test are written, not run.
* Gas/source multirate is out of scope (hierarchical KDK is DMO-only).

See the [round-3 record](repair/multirate_stabilization_round3_20261026.md).

20261026 AREPO-grade campaign (source only; NOT compiled, NOT tested; owner
validation required). Repairs:

* Rank-divergent PM routing admission and stale routing plans.
* Consensus-batch lane-count and failure divergence.
* Remote group-walk coincident-pair accounting.
* O(world_size) LET bookkeeping.
* The per-tick telemetry collective.
* The FD halo allocation outside the governor.

New features:

* Multirate outer PM level (`multirate_outer_kick`).
* Exact drift-factor memo.
* Certified uniform-mass reuse.

New tests: routing admission (MPI), in-place FFT layout, remote group walk
(MPI), fine-tick caches, gravity invariants, multirate workflow gates. New
tool: `tools/compare_dmo_arepo.py`. See the
[campaign record](repair/arepo_grade_campaign_20261026.md).


20261025 TreePM optimization campaign (source only; NOT compiled, NOT
tested, owner validation required): opt-in finite-difference PM gradient,
lazy `potential_k`, governed communication-arena growth, group-walk leaf
culling/target sizing/telemetry (and repair: the group walk now also serves
owned targets in distributed runs), LET graph subset reuse, fused diagnostic
reductions, independent checkpoint cadence, structural snapshot verification
and `numerics.gravity_timestep_eta`. See the
[campaign record](repair/treepm_optimization_campaign_20261025.md) for the
implementation matrix, risks and validation checklist.
Round 2 of the same campaign (also NOT compiled, NOT tested): batched rank
consensus (one MIN allreduce of `[v, ~v]`) for gravity entry/PM decisions,
tree topology-generation caches, stable radix target ordering, in-place FFTW
PM storage, workload-aware PM routing reserves, remote LET group walk,
experimental `numerics.treepm_pm_reuse_policy = stale_within_cadence`,
`output.restart_interval_time_code`, and uniform-mass lane release.

20261015 TreePM performance campaign: tree-ordered residual traversal
(bitwise-identical forces), threaded FFTW, opt-in adaptive estimator and
cutoff-node policy, opt-in hierarchical RMS coarse criterion and the B-10
post-bootstrap governor reconcile are implemented, compiled and tested; see
the [campaign record](repair/treepm_performance_campaign_20261015.md) for
measured accuracy/timing, exact commands and the remaining open items.
Round 4 (source-level, compiled in all three local build trees, NOT yet
executed): hierarchical softened-encounter rung criterion (D1,
`numerics.hierarchical_encounter_timestep_enabled`, default on) and
hierarchical endpoint force reuse (R1) with restart schema v24 and TreePM
short-range response wire v2. Open until the owner's test run: the D1 gate
(`integration_hierarchical_treepm_workflow`, refined error < 1e-3), hierarchical
restart equivalence with reuse, and the full HDF5/FFTW and MPI suites;
trajectory/growth qualification of the fast deck and of hierarchical
`rms_split_scale` remains open.

20261012 test regression closure: the source repairs in the
[closure handoff](repair/treepm_deep_hotpath_test_regression_closure_20261012.md)
await compilation and execution of every stated gate. In particular, active
FFTW-MPI slab/transposed allocations, Release assertion effectiveness, Ewald,
and hierarchical dispatcher/restart qualification have no new runtime evidence.
The authoritative scheduler, gravity mathematics, config and memory estimators
were not changed. Any additional failure needs its exact command and log before
a numerical or backend conclusion can be drawn.

20261012 Attempt 1 source repairs and optimizations await all execution gates:
composed runtime first light, analytic/table/SIMD/shared-derivative accuracy,
MPI collective failure/rank equivalence, restart migration, actual snapshot
units/softening, adaptive 0.25 and separately 0.5, evolved trajectories,
split-dependent independent total references and measured live-set memory.
Cross-instance table pooling, detailed table-hit telemetry and PM workspace
redesign are not implemented; no performance evidence justifies extra O(N)
storage or new defaults. Compilation/tests/benchmarks were explicitly deferred
by user. See the [exact status table](repair/treepm_deep_hotpath_attempt1_closure_20261012.md).


20261010 P1–P5 corrective follow-up: Python 3.8-compatible complete-snapshot
reading, coincident residual/zero-vector repair, constrained verified-v23
configuration migration, endpoint PM guards, mandatory Ewald recovery profiles,
and test comparison/assertion/CTest repairs are implemented in source. All
compilation, tests, MPI, convergence and performance remain **NOT RUN — deferred
by user**. Legacy opaque-hash migration fails closed; historical external P0
tools remain absent; lookup and hierarchical modes remain opt-in/provisional.
See the [defect dispositions, exact manifest and commands](repair/treepm_p1p5_qualification_and_test_reliability_repair_20261010.md).

## Unified TreePM recovery qualification pending — 20261010

The [source handoff](repair/treepm_unified_performance_recovery_20261010.md)
implements P1–P5 changes and tests without executing compilation or validation,
as requested. Gaussian lookup is opt-in; full overlapping MAC counters remain
default; hierarchical KDK stays provisional with rung-zero production default.
Pending gates include full-snapshot force/timing comparisons, FFTW Ewald,
MPI/thread equivalence, governed allocation failure, hierarchical restart and
growth/power refinement. The historical P0 snapshot/oracle tools are absent from
this checkout; the new frozen-source driver does not impersonate them. Nonlinear
long-time and isolated-orbit hierarchical qualification remain open extensions.

## P6 qualification pending — campaign label 20261006

P6 is source implementation only. No P6 subsystem has new validation evidence.
Both conservative motion refit and hierarchical DMO KDK now have optional
source paths; neither is certified. Required gates are in
[p6_validation_campaign_plan_20261006.md](repair/p6_validation_campaign_plan_20261006.md).

Remaining implementation scope: gas/source multirate KDK needs a coupled
hydro/source synchronization law; motion refit deliberately rejects unknown
row identity and heterogeneous softening; spatial feedback is cold after
restart; duplicate owner discovery is retained to avoid population-scale
peer storage or a protocol rewrite. Deposition/FFTW threading and resulting-rank
process-fit prediction remain deferred. Bounded event payload memory is
reported but does not have prospective per-event admission. New defaults stay
conservative. These statements supersede earlier roadmap intent, not earlier
independent qualification evidence.

> **Historical / non-authoritative.** This document records an earlier campaign state. Use [`CURRENT_STATUS.md`](../CURRENT_STATUS.md) for current repository truth.



## 2026-10-02 P5-Q3 remaining qualification gates

The known P5-Q2 source backlog addressed by P5-Q3 is closed for the exercised
CPU and serial-HDF5 paths: compact periodic-seam top-domain geometry, scheduler
promotion for future-activated source births, newborn drift/post-kick epoch
handling, evolved-population final ownership validation, and stale CPU fixtures.
DMO and TreePM single-rank restart equivalence pass exactly. Remaining gates are
qualification rather than known source defects: run the MPI+OpenMP force-vector
matrix, event-by-event collective fallback/recovery sequence, MPI restart
publication smoke, and Parallel-HDF5 provider-family matrix on a host providing
MPI C++, Parallel HDF5, FFTW and FFTW-MPI. Do not start P5-W as a substitute for
those gates.

## 2026-09-21 first-light science-I/O closure remaining qualification gates

Source-level closure now covers order-independent zero-based external IDs, numerically consistent distributed mass auditing, one-file Parallel-HDF5 science-output topology, transactional publication, and final-endpoint output. Remaining items are qualification/scale extensions rather than permission to fall back silently:

- execute the registered 2/3/4/8-rank snapshot/IC matrix and real 64^3 three-rank TreePM smoke on a host with MPI C++, Parallel HDF5, FFTW and FFTW-MPI;
- directly open a produced science snapshot in GLnemo2 when that executable is available; ordinary HDF5/h5py compatibility is the automated floor;
- implement and qualify the reserved fixed-count `aggregated` large-rank backend before claiming scalable file-count aggregation beyond the single-file small/moderate-rank policy.

## 2026-09-17 DMO first-light recovery remaining gates

The source-level first-light blockers for config-driven rung-zero stepping,
scale-factor endpoint authority, monofonIC structural ingestion, periodic state,
and flat snapshot-set discovery are repaired. The remaining acceptance blockers
are environment gates: this repair host lacks a discoverable MPI C++ toolchain
and FFTW development library, so real multi-rank timestep/snapshot execution and
the FFTW-backed `cosmosim_harness` TreePM smoke were not executed. These are
validation gates, not permission to weaken the repaired contracts.


## 2026-07-30 Campaign B packaging closure and remaining MPI gate

The adversarial handoff defects are closed in source: the two actual ADS
sidecars are removed; the packaging path no longer excludes legitimate
`core` directories; one shared sentinel inventory is checked before and after
compression; staged/extracted paths, symlink targets, and SHA-256 values are
compared; and a fresh extracted copy passes hygiene, CPU configuration, and the
required core builds. `integration_source_package_completeness` passes directly
and under both CPU and HDF5 CTest inventories.

Current supported evidence is CPU 133/133, HDF5 137/137, and the focused
Campaign B HDF5 set 7/7. The only Campaign B acceptance item still
**environment-blocked** is real MPI execution: this environment cannot find
`mpi-cxx` or `MPI_CXX`, so np1/np2/np4 and fault-injection tests were not
launched. `distributed_ic_import` must remain `provisional` until those exact
commands pass on a dependency-complete host.

The immutable final archive SHA-256 is reported by the packaging script and
handoff response rather than embedded in source documentation, because embedding
an archive's own hash would mutate the archive.

## 2026-07-26 Campaign B distributed final-gate evidence

The final source micro-patch now validates every non-null supplied manifest,
uses `manifest.dialect` for distributed inspection, executes canonical bundle
verification on supplied canonical manifests, guards the remaining exact-audit
capacity/path accounting seams, and instruments actual production MPI
collective calls separately from logical consensus phases. The registered MPI
matrix asserts 30 calls per routing batch and the exact non-routing formula,
including duplicate-audit rounds and chunked Bcasts. Stable per-file readers,
the three-hash-pass ceiling, exact count/mass/ownership/route/ID audits,
source identity, large-count protection, and the decomposed IC structure remain
intact.

The remaining Campaign B item is external runtime evidence, not a source waiver:
this repair environment has no `MPI_CXX`, so the registered np1/np2/np4
canonical-manifest, malformed-input, fault-injection, counter-agreement, and
scaling matrix cannot execute here. Distributed capability must remain
`provisional` until those exact tests pass on a dependency-complete host. The
final handoff must also be the exact archive produced by
`scripts/ci/package_source_zip.sh`, not an independently uploaded working-tree
ZIP.

## 2026-07-08 remaining MPI acceptance limitations

The complete local `mpi-hdf5-fftw-debug` inventory now passes 150/150, and the
focused gravity/MPI matrix passes 36/36 through np4. The items below remain open
because their scientific/restart oracles are narrower than the claimed
production behavior, not because the registered local executables are red.

- The new `integration_reference_workflow_distributed_amr_mpi_two_rank` test is
  registered as a real MPI-launched test and exercises directed AMR payload and
  owner-routed flux exchange, but it does not yet prove full production
  `ReferenceWorkflowRunner` AMR restart equivalence with pending flux-register
  persistence.
- `tests/integration/test_distributed_gas_cell_migration_mpi.cpp` still needs a
  stronger named patch/gas-cell cross-rank movement oracle with stable-ID keyed
  physical, sidecar, scheduler, and restart comparisons.
- `tests/integration/test_reference_workflow_distributed_hydro_mpi.cpp` still
  needs full conserved-state and correction-sign assertions for mass, momentum
  x/y/z, total energy, and internal-energy behavior.
- Periodic PM still uses sparse payload records over global all-to-all
  collectives. This is documented as a small-cluster transport policy only.


## 2026-05-25 Stage 6 memory accounting status

| ID | Status | Area | Finding | Resolution / Remaining action |
| --- | --- | --- | --- | --- |
| S6-07-MEMORY-ACCOUNTING | Open | Runtime memory diagnostics | Core now reports capacity-based owned bytes for SimulationState/workspace lanes and keeps unknown external allocations explicit, but tree/PM/MPI/IO subsystems still need direct report hooks for concrete owned buffers. | Keep central API stable; add per-subsystem report hooks to replace unknown estimates when each module exposes ownership details. |

## 2026-06-16 H3 AMR hydro closeout status

| ID | Status | Area | Finding | Resolution / Remaining action |
| --- | --- | --- | --- | --- |
| H3-09-AMR-HYDRO-RESTART-EQUIVALENCE | Open | AMR hydro restart evidence | The production AMR hydro path is locally conservative at the implemented single-rank level, with tests for conservative refine/derefine, same-level and coarse-fine ghost fill, live flux-register generation, automatic reflux, shock/Sedov guards, synchronization stress, and AMR patch migration payloads. However, no dedicated production AMR hydro run/restart/run equivalence test currently exercises patch geometry, ghost fill, flux registers, reflux, and stable `gas_cell_id` identity across restart. | Add and run a dedicated HDF5 restart-equivalence test for `advanceProductionAmrHydro` on explicit AMR patch geometry before claiming H3 final restart closure. Current audit evidence is in `docs/architecture/h3_amr_hydro_closeout_audit.md`; local validation was environment-blocked because `cmake` was not found on PATH. |

# Repair open issues (P01–P19 freeze ledger)

_Date captured: 2026-04-14 (UTC)_

## 2026-05-21 Stage 5 config-contract audit status

| ID | Status | Area | Finding | Resolution / Remaining action |
| --- | --- | --- | --- | --- |
| S5-00-CONFIG-CONTRACT-AUDIT | Open | Config contract hardening plan | Stage 5 needed a concise implementation map naming mixed raw/normalized/derived/runtime lanes, risky call sites, and entrypoints that must remain frozen-config-only. | Added `docs/repair/stage5_config_contract_audit_map.md` with file-by-file patch map and ordered follow-up prompts (S5-01..S5-04). Next prompts should implement targeted guards/tests without broad refactor churn. |


## 2026-05-19 Stage 4 orchestrator dispatch status

| ID | Status | Area | Finding | Resolution / Remaining action |
| --- | --- | --- | --- | --- |
| P41-STAGE4-ORCHESTRATOR-DISPATCH-041 | Closed | Runtime orchestrator dispatch | The production orchestrator broadcast every stage to every callback, relying on callback-side self-filtering. | Repaired with typed per-stage handler registration and dispatch buckets, debug/contract failures for impossible off-stage direct calls, targeted unit coverage for registered-only delivery, unregistered stage absence, and deterministic per-stage order. Continue to reject future production broadcast/self-filter orchestration. |


## Stage 0 gate status update (2026-04-26 UTC)

- **Stage 0 status: NOT CLOSED.**
- Gate evidence: `docs/repair/stage0_runtime_truth_freeze.md`.
- New physics/performance work remains blocked pending closure of Stage 0 blockers below.

### Stage 0 closure blockers (this gate run)

| ID | Status | Area | Current limitation | Next-step evidence target |
|---|---|---|---|---|
| P0-10-STAGE0-GATE-028 | Closure candidate | Stage 0 closure gate | P0-09/P0-10 repair applied: missing test sources restored, schema drift repaired, runtime smoke normalized-config check aligned, and softening restart payload metadata fixed. Full closure still requires local completion of CPU/HDF5/PM preset suites. | Run the full command bundle from `docs/repair/stage0_runtime_truth_freeze.md`; mark Stage 0 closed after those suites pass locally. |
| P32-SOFTENING-PROVENANCE-029 | Closed | Softening + restart/provenance contract | Targeted CPU and HDF5 validation now pass for softening ownership, reorder/sidecar preservation, migration preservation, restart, and snapshot schema coverage; `integration_softening_ownership_invariants` builds a valid restart payload with matching normalized config hash/provenance hash and single-rank distributed restart state. | Continue to include `ctest --preset test-hdf5-debug --output-on-failure -R "restart|snapshot|softening"` in release-gate evidence. |
| P33-APP-CONFIG-ROUNDTRIP-030 | Closed pending local full-suite certification | Runtime app normalized config contract | `integration_runtime_app_smoke` now checks canonical axis-aware normalized config key `treepm_pm_grid_nx = 24`, matching the documented normalized-config contract. | Re-run runtime app smoke, config examples, and config parser tests locally. |
| P34-RELEASE-SCHEMA-DRIFT-031 | Closed pending local full-suite certification | Release/docs schema metadata drift | Release manifest, docs scaffold, and release readiness checks now agree on snapshot `gadget_arepo_v4`, restart `cosmosim_restart_v12`, and provenance `provenance_v5`. | Re-run docs scaffold, release readiness, and snapshot/restart schema tests locally. |
| P35-PM-VALIDATION-CONTRACT-032 | Closed | PM isolated/coupling/validation contract | PM+HDF5+FFTW validation now accepts either internal-node residual pruning or leaf-pair cutoff skips as cutoff traversal evidence, matching the TreePM diagnostic contract and closing the single-rank residual-cutoff expectation mismatch. | Full PM+HDF5+FFTW preset is green with `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure`; MPI+FFTW remains environment-blocked when `fftw3_mpi` is unavailable. |
| P36-STAGE2-TIMESTEP-RESTART-MIRROR-033 | Closed | Stage 2 timestep authority / restart mirror bypass | Adversarial Stage 2 review found that restart validation checked `CellSoa::time_bin` only when scheduler `bin_index` was cell-sized, allowing mixed particle/gas restarts with stale cell mirrors to pass if particle mirrors matched scheduler truth. | Patched restart validation/import to map cells through parent gas particle scheduler entries, added stale-cell-mirror rejection in `test_unit_restart_checkpoint_schema`, refreshed restart roundtrip fixtures to sync particle and cell mirrors from scheduler authority, documented compatibility behavior in `docs/output_schema.md`, and recorded the adversarial review in `docs/architecture/stage2_adversarial_timestep_review.md`. CPU and HDF5 debug gates passed locally. |


## Current blocker ledger after dependency-enabled validation

- TreePM Phase 3 implementation/closure work remains open by contract: see `docs/treepm_phase3_contract.md` hard-gate requirements and staged maturity map.
- Production rejects PM cadence other than one and
  `hierarchical_max_rung != 0`; asynchronous/mixed-rung TreePM remains
  unimplemented rather than an accepted approximation.
- Zoom focused-PM no longer relies on configured spherical membership alone; file-driven high-resolution particle masks are now consumed directly. Remaining follow-up scope is richer irregular region metadata beyond particle-ID masks.
- Heavy gravity abnormality reference checks are intentionally policy-gated (`analysis.diagnostics_execution_policy = all_including_provisional`) and remain non-default for production cadence costs.


Environment blockers currently observed in this validation environment:


- MPI gravity gate now expects `validation_phase2_mpi_gravity_two_rank` and `validation_phase2_mpi_gravity_three_rank` in MPI-enabled presets; if MPI runners have fewer than 3 ranks available to `mpiexec`, this gate will fail by construction and should be treated as an environment-capacity blocker, not a silent skip.
- `cmake --preset pm-hdf5-fftw-debug` can fail when FFTW3 development files are unavailable.
- `cmake -S . -B build/mpi-fftw-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCOSMOSIM_ENABLE_TESTS=ON -DCOSMOSIM_ENABLE_MPI=ON -DCOSMOSIM_ENABLE_FFTW=ON -DCOSMOSIM_ENABLE_HDF5=OFF` fails in this container because MPI C++ tooling is missing (`Could NOT find MPI_CXX`), preventing two-rank PM distributed runtime closure in CI-like local validation.
- The `cpu-only-debug` preset used for quick repair validation sets `COSMOSIM_ENABLE_MPI=OFF`, so MPI two-rank tests such as `integration_reference_workflow_distributed_treepm_mpi_two_rank` are not registered in that preset.

Deferred to Phase 3 follow-up after axis-aware geometry repair:

- Full 2D pencil decomposition beyond FFTW transposed ownership remains open. Current `pencil` mode is now an end-to-end transposed PM path with explicit ownership contracts, but wider process-grid scaling and non-FFTW backend parity remain follow-up scope.
- Isolated/non-periodic PM is implemented locally and as a bounded distributed
  root gather/solve/scatter compatibility path. A scalable distributed
  open-boundary PM strategy remains open.
- CUDA PM CIC deposition/interpolation now accepts axis-aware periodic box lengths; isolated/open CUDA PM remains unsupported and fails fast.

No new in-tree CPU-path code blocker was observed on the repaired PM unit/integration single-rank checks.

Additional validation limitation for this pass:

- `test_integration_tree_pm_coupling_periodic` in `cpu-only-debug` remains expensive in this environment (naive DFT PM backend with MPI/FFTW disabled), so the performance-hardening patch evidence is currently anchored on `unit_pm_solver` plus benchmark runs (`bench_pm_solver`, `bench_tree_pm_coupling`, `bench_parallel_decomposition_exchange`) rather than full integration runtime closure in this container.
- `test_integration_tree_pm_coupling_periodic` also remains runtime-heavy after quadrupole + mature MAC tree upgrades; in this container it exceeded a 300s direct invocation timeout, so distributed short-range agreement evidence continues to rely on existing MPI-enabled validation gates (`validation_phase2_mpi_gravity_*`) outside the CPU-only preset.

## Recently closed items

| ID | Status | Area | Prior symptom | Closure evidence |
|---|---|---|---|---|
| P0-09 | Closed | Restart/reload runtime-truth round-trip invariants | Restart continuation safety needed explicit invariant checks for runtime-truth domains (particles/species/order, scheduler/bin/active set, softening overrides, gas identity, normalized config/provenance hashes) plus clear missing/legacy-field behavior. | Extended `test_integration_restart_checkpoint_roundtrip` with deterministic restart-truth checks, active-set reconstruction equivalence checks, required-field clear failure checks (`pending_bin_index`), and documented legacy compatibility coverage for missing per-particle softening dataset; updated `docs/restart_checkpointing.md` and `docs/architecture/adr_runtime_truth_ownership.md`; validated with `cmake --build --preset build-cpu-debug -j4 --target test_integration_restart_checkpoint_roundtrip` and `ctest --preset test-cpu-debug --output-on-failure -R "restart|checkpoint|roundtrip|provenance"`. |
| P0-08 | Closed | Active-set construction ownership, cache invalidation, and no-competing-authority tests | Active-set authority and solver-local derivation/invalidation rules were documented at high level but lacked explicit test-floor coverage for stale mirror divergence and mutable compact-view generation invalidation after reorder/resize. | Added deterministic active-set ownership tests in `test_unit_time_integration` (`testActiveSetAuthority`, `testActiveSetNoCompetingBuilders`) and mutable compact-view invalidation/derivation checks in `test_unit_simulation_state` (stale particle/cell view rejection after reorder/resize and solver-local view derivation from authoritative indices); updated `docs/architecture/adr_runtime_truth_ownership.md` with explicit owner/invalidation/generation/lifetime rules for active-set layers and forbidden module-local builders. |
| P0-07 | Closed | Config ownership: raw/normalized/derived/runtime/provenance lanes | Config-derived runtime values and continuation metadata ownership boundaries needed explicit classification + deterministic test floor for normalization/hash/restart metadata consistency and ambiguous time/scale semantics. | Added targeted ownership tests in `test_unit_config_parser` and `test_unit_units_cosmology_provenance`, plus restart normalized-text/hash mismatch rejection in `test_unit_restart_checkpoint_schema`; updated `docs/architecture/adr_runtime_truth_ownership.md` and `docs/configuration.md` with explicit ownership map and anti-ambiguity semantics for code-time vs scale-factor vs diagnostic redshift lanes. |
| P0-06 | Closed | Softening ownership/priority and override-preservation invariants | Softening ownership order (global/species/per-particle) and transform persistence rules were under-specified for resize/reorder/migration/restart/active extraction and diagnostics observer behavior. | Added `test_integration_softening_ownership_invariants` with `test_softening_priority_invariants`, `test_softening_override_resize_reorder_preservation`, and `test_softening_override_restart_roundtrip`; verified per-particle > species > global resolution, resize/reorder/migration persistence, restart roundtrip persistence, active-lane mirror behavior, and diagnostics non-authority; updated ADR softening preservation rules and mirror-authority notes. |
| P0-05 | Closed | Gas-cell identity and hydro-state ownership invariants | Gas identity semantics and temporary 1:1 gas-particle/gas-cell coupling were insufficiently explicit under reorder/resize/restart/active-hydro extraction paths. | Added `test_unit_gas_cell_identity_invariants` plus restart mapping checks in `test_integration_restart_checkpoint_roundtrip`; hardened the contract with `requireParticleBoundGasCellContract`, explicit row/ID helper seams, hydro active-view validation, migration-by-ID coverage in `test_integration_species_migration_invariants`, and reorder guardrails preventing silent gas-relative reorders; documented temporary contract + forbidden assumptions in `docs/architecture/adr_runtime_truth_ownership.md`. |
| P0-04 | Closed | Species migration + canonical grouping invariants | Species transitions were under-tested for canonical grouping continuity across sidecars/count ledgers/offset index maps, per-particle softening override survival, and scheduler-bin/active-set continuity after migration + reorder/resize. | Added deterministic migration invariant suite in `test_integration_species_migration_invariants` (single/batch migration, reorder/resize/active-extraction adversarial sequences, sidecar attach-detach checks, timestep-bin continuity, and softening-priority checks); repaired migration pack/commit lanes to carry optional per-particle softening overrides; updated runtime ownership docs for allowed migration path + sidecar rules + species-agnostic active-eligibility policy. |
| P0-03 | Closed | Particle ordering/resize/reorder/sidecar stale-index invariants | Index-based assumptions could silently desynchronize particle identity from species tags, sidecar payload rows, softening overrides, and cached active/kernel views across resize and reorder operations. | Hardened `reorderParticles` to move complete star/BH/tracer payload rows under `kMoveWithParent`, keep default parent-indirection rows stable while remapping `particle_index`, preserve optional softening sidecar values across particle resize, add index-generation invalidation and defensive species-sidecar ownership checks, and expand `test_integration_reorder_compaction_sidecars` coverage for reverse permutations, move-vs-indirection semantics, randomized reorder modes, compaction, and stale-view invalidation; validated with targeted build/ctest runs plus `rg "particle_index" src include tests` audit. |
| P32-TIMEBIN-AUTHORITY-027 | Closed | Timestep/bin ownership invariant coverage | Scheduler bin authority vs state mirror lanes lacked deterministic adversarial tests for unauthorized mutation, reassignment consistency, reorder pairing, and persistent-state restart equivalence. | Added targeted invariant coverage in `test_unit_time_integration` (`testTimestepBinAuthorityInvariant`, `testTimestepBinReassignmentAndRestartRoundTrip`, `testTimestepBinReorderIdentitySurvival`), plus new mirror-authority helper checks in `core/time_integration`; validated with `cmake --build --preset build-cpu-debug --target test_unit_time_integration` and `ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_time_integration|integration_hierarchical_time_bins|integration_restart_checkpoint_roundtrip"`. |
| P19-STATE-OWNERSHIP-005 | Closed | Core state ownership/invariants decomposition | `src/core/simulation_state.cpp` concentrated ownership rules, species indexing, transfer packing, metadata parsing, and active-view hot/cold logic in one implementation unit. | Split into focused implementation units; explicit gravity/hydro hot-field contract notes added; targeted invariant/contract checks added in `tests/unit/test_simulation_state.cpp`; verified by `cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state` and `./build/cpu-only-debug/test_unit_simulation_state`. |
| P19-IO-CONTRACT-006 | Closed | Snapshot vs restart contract separation | Restart continuation-critical metadata/scheduler requirements and failure behavior were under-specified relative to snapshot interoperability boundaries. | Added shared continuation metadata contract helper, explicit restart completeness checklist, and negative tests for schema mismatch/missing required scheduler lane/finalize failure behavior; verified via `test_unit_snapshot_hdf5_schema`, `test_unit_restart_checkpoint_schema`, and `test_integration_restart_checkpoint_roundtrip`. |
| P19-GATE-FFTW-TEST-001 | Closed | PM FFTW analytic mode response | `unit_pm_solver` failed at `tests/unit/test_pm_solver.cpp:82` (`cosine_similarity > 0.98`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R "unit_pm_solver|integration_tree_pm_coupling_periodic"` now passes. |
| P19-GATE-FFTW-TEST-002 | Closed | TreePM periodic PM long-range coupling (FFTW path) | `integration_tree_pm_coupling_periodic` failed with `rel_l2=18129.9` (required `<= 0.75`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` now passes (`36/36`). |
| P19-ARCH-CORE-BOUNDARY-003 | Closed | Core dependency direction | `include/cosmosim/core/reference_workflow.hpp` assembled analysis/I/O/physics workflow concerns inside `core/`. | Workflow assembly moved to `workflows/` and `integration_core_dependency_direction` guard added to fail on forbidden upward includes. |
| P19-CONFIG-CONTRACT-004 | Closed | Frozen configuration contract typing | Policy keys (`gravity_solver`, `hydro_solver`, `coordinate_frame`, boundaries, feedback mode/variant) remained string contracts after parse. | `unit_config_parser`, `unit_simulation_mode`, and `integration_simulation_mode_toy_runs` now validate enum-backed freeze behavior with deterministic normalized text/hash and centralized key+alias registry handling. |
| P19-OBS-EVENTS-007 | Closed | Operational observability | Runtime diagnostics lacked a structured run-health/event surface linked to provenance; failures were harder to audit from CI artifacts. | Added compact runtime event model in `core::ProfilerSession` and `reference_operational_events.json`. The gravity PM-field event now takes pre-commit opportunity/version/build epoch from the pending refresh directive instead of stale committed cadence state, and diagnostic doubles use scientific `max_digits10` formatting so small nonzero evidence remains visible. Cadence ownership is unchanged. |
| P19-DIAGNOSTICS-TIER-008 | Closed | Diagnostics maturity/scalability labeling | Diagnostics bundles mixed run-health outputs with heavy reference algorithms and did not encode maturity tier in code/output policy. | Added typed analysis execution policy, per-diagnostic tier/maturity/scalability metadata in bundle JSON, and tests proving provisional heavy diagnostics are non-default unless explicitly opted in (`all_including_provisional`). |
| P19-PARALLEL-CONTRACT-009 | Closed (R07/R07b/R07c hardened) | Distributed-memory ownership/determinism contracts | Owned-vs-ghost residency and reduction/config consensus expectations were partially implicit in ghost-owner vectors and ad-hoc checks. | Added explicit send/receive transfer-role descriptors, lifecycle intent + post-transfer residency typing, strict plan/descriptor consistency validation, unambiguous reduction baseline/measured naming, explicit reduction policy modes, and property-level config mismatch records with rank/value payloads; validated in `test_parallel_distributed_memory` and `test_parallel_two_rank_restart`. |
| P19-INFRA-GATE-010 | Closed | Infrastructure gate clarity + false-closure prevention | CPU-only pass signals could be misread as full closure when HDF5 or PM/HDF5/FFTW path validation was not surfaced as one explicit gate report. | Added `scripts/ci/enforce_infra_gates.sh` + CI `infrastructure_gates` job to run explicit CPU/HDF5/PM+HDF5+FFTW gate sets, write `infrastructure_gate_report.json`, and fail closure when any feature path gate fails; strengthened `integration_core_dependency_direction` with CMake target-link direction checks. |
| P20-CONFIG-ROUNDTRIP-011 | Closed | Frozen config canonicalization/provenance | Normalized config emitted non-parseable empty-value/self-hash combinations and could not round-trip through the same parser. | Removed self-hash embedding from normalized text, preserved empty-value parsing for `key =` canonical lines, hardened inline-comment stripping for URLs/path-like values, and added deterministic round-trip/hash checks in `test_unit_config_parser`. |
| P20-RESTART-STELLAR-012 | Closed | Restart completeness/integrity | Restart payloads omitted stellar-evolution bookkeeping lanes and the payload hash did not cover those lanes. | Bumped restart schema to `cosmosim_restart_v3`, serialized/read all stellar-evolution sidecar lanes, extended integrity hashing coverage, and verified round-trip in `test_integration_restart_checkpoint_roundtrip` + `test_unit_restart_checkpoint_schema`. |
| P20-SNAPSHOT-MASSTABLE-013 | Closed | Snapshot import safety | HDF5 snapshot import defaulted missing `Masses` datasets to zero instead of using `Header/MassTable`. | Snapshot reader now requires a positive `Header/MassTable` entry for fallback and reports `Masses=MassTable`; verified in `test_integration_snapshot_hdf5_roundtrip`. |
| P20-PARALLEL-HONESTY-014 | Closed | Distributed-memory scaffolding | Ghost descriptor planning treated local ghost rows as outbound payload rows, and distributed restart decode accepted missing ownership entries as rank 0. | Descriptor-only ghost planning now leaves outbound payload indices empty unless owned export rows are known, duplicate/missing restart ownership entries are rejected, and tests were updated in `test_unit_parallel_distributed_memory`. |
| P20-STATE-COVERAGE-015 | Closed | Species-sidecar invariants | Species-tagged star/BH/tracer particles could exist without exactly one corresponding sidecar row. | Ownership validation now enforces one-to-one species-tag/sidecar coverage and rejects duplicates; verified in `test_unit_simulation_state`. |
| P21-TREEPM-PHASE2-CONTRACT-016 | Closed | Distributed gravity contract freeze surface | Phase 2 distributed TreePM ownership and config/build/documentation surfaces were implicit and leaned on optional MPI release smoke. | Added typed numerics keys (`treepm_pm_decomposition_mode=slab`, `treepm_tree_exchange_batch_bytes`), provenance capture, `mpi-hdf5-fftw-debug` preset, CI matrix coverage for that preset, and explicit architecture freeze doc `docs/treepm_phase2_distributed_contract.md`; validated with `unit_config_parser`, `integration_docs_scaffold`, `integration_provenance_roundtrip`, and `unit_parallel_distributed_memory`. |
| P22-TREEPM-SLAB-LAYOUT-017 | Closed | PM slab ownership/storage contract | PM mesh ownership/indexing semantics for distributed slab decomposition were implicit and PM storage assumed globally replicated cube allocation. | Added typed `parallel::PmSlabLayout` ownership helpers (partitioning, owner lookup, global/local conversions), adapted `PmGridStorage` to allocate local slabs, kept one-rank full-domain behavior as degenerate slab, and added unit coverage in `test_unit_parallel_distributed_memory` + `test_unit_pm_solver` with explicit docs updates. |
| P23-RUNTIME-TOPOLOGY-018 | Closed | MPI/CUDA runtime infrastructure wiring | MPI slab ownership, `parallel.gpu_devices`, and CUDA PM execution policy existed as disconnected surfaces with no shared runtime topology/selection contract. | Added `core::CudaRuntimeInfo`, explicit `parallel::DistributedExecutionTopology` + rank/device assignment helpers, runtime world-size validation, CUDA device selection wiring in the reference workflow, new debug presets (`cuda-debug`, `mpi-cuda-hdf5-fftw-debug`), and targeted config/topology tests. |
| P24-TREEPM-DENSITY-DEPOSIT-019 | Closed | Phase 2 distributed PM density assignment | `PmSolver::assignDensity` previously rejected partial slabs and could not route owner contributions to slab owners. | Added slab-owner routed PM density contribution exchange with batched per-destination payloads and owner-side validation in `PmSolver::assignDensity`; added CIC/TSC wrap/boundary/mass-conservation and one-rank vs distributed agreement checks in PM tests; updated PM/distributed-memory contracts docs. |
| P25-TREEPM-INTERPOLATION-020 | Closed | Phase 2 distributed PM interpolation return path | PM interpolation previously required full-domain mesh residency on every rank and rejected partial slabs. | Added explicit reverse communication gather for force/potential where particle owners send stencil requests to slab owners and receive weighted contributions; validated CIC/TSC two-rank vs one-rank agreement and slab-boundary cases in `test_integration_pm_periodic_mode`. |
| P26-TREEPM-SHORT-RANGE-021 | Closed | Phase 2 distributed tree short-range residual | TreePM short-range residual previously evaluated only local sources and lacked peer-evaluated export/import tree participation. | Added explicit active-target export/import protocol with per-peer batching (`tree_exchange_batch_bytes`), remote peer tree evaluation, response coverage checks for duplicate/missing packets, workflow config wiring, and two-rank distributed-vs-single-rank coverage in `test_integration_tree_pm_coupling_periodic`. |
| P27-TREEPM-OWNERSHIP-022 | Closed | Gravity-aware decomposition + particle migration ownership commit | Decomposition cost used only generic work/memory terms, and migration intents existed without a committed owner-tag + sidecar rebuild boundary. | Added gravity-aware weighted decomposition terms (`owned_particle`, `active_target_count_recent`, `remote_tree_interactions_recent`, `memory_bytes`), per-rank decomposition component metrics, plus explicit `SimulationState` migration pack/commit APIs that rebuild ownership/species/sidecar mappings at a synchronization point; validated by `test_unit_parallel_distributed_memory` and `test_unit_simulation_state`. |
| P28-TREEPM-RESTART-CONTINUATION-023 | Closed | Distributed TreePM restart/provenance continuation safety | Restart payloads lacked explicit versioned distributed gravity continuation state, and rank/layout mismatches were hard to debug after resume. | Bumped restart schema to `cosmosim_restart_v4`; persisted `DistributedRestartState` (epoch, owning ranks, PM slab layout, cadence + field metadata, explicit restart policy); added restart compatibility diagnostics and mismatch coverage in `test_unit_parallel_distributed_memory`, plus restart round-trip validation in `test_integration_restart_checkpoint_roundtrip` and schema checks in `test_unit_restart_checkpoint_schema`. |
| P29-TREEPM-PHASE2-CLOSEOUT-024 | Closed | Phase 2 closeout coherence + explicit failure contracts | Final Phase 2 hard-gate checklist and some distributed failure conditions were not captured as one explicit closeout artifact with command-backed checks. | Added `docs/treepm_phase2_closeout.md`; hardened restart cadence-state validation and compatibility reporting (`pm_cadence_steps_match`); added explicit negative-contract coverage in `test_validation_phase2_mpi_gravity` and `test_unit_parallel_distributed_memory` for rank-count mismatch, unsupported decomposition mode, communicator/layout mismatch, missing distributed restart metadata, and cadence inconsistency. |


## Verified evidence (this run)

- CPU targeted repair floor passes:
  - `ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_simulation_state|unit_parallel_distributed_memory|unit_restart_checkpoint_schema"`
- HDF5 targeted repair floor passes:
  - `ctest --test-dir build/hdf5-debug --output-on-failure -R "unit_restart_checkpoint_schema|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip"`
- PM HDF5+FFTW configure remains blocked by missing FFTW3 development files in this environment:
  - `cmake --preset pm-hdf5-fftw-debug`

## Outcome tags

- **CPU-only preserved**
- **HDF5 path preserved**
- **PM HDF5+FFTW path environment-blocked**
- **P19 blocker ledger cleared**

## R08 hardening residual limitations (tracked, non-blocking)

- `integration_core_dependency_direction` target-link guard now consumes CMake Graphviz output from the configured build tree, which is materially stronger than top-level text scanning; however, it still relies on target-name pattern matching (`cosmosim_(analysis|io|physics|workflows|app)`) rather than full semantic layer tags. If naming conventions change, this check must be updated in lockstep.
- Infrastructure gate test scopes are now manifest-driven and converted into a ctest regex at runtime. This removes inline script fragility but still depends on stable test identifiers.

- Isolated/open PM remains single-rank in this stage; distributed isolated PM (MPI slabs) is intentionally blocked with runtime errors until a distributed open-boundary PM strategy is implemented.

- Tree softening maturity limitations after stage 2026-04-20:
  - Node-level accepted multipoles currently use `epsilon_node_max` as a conservative aggregate for mixed-source nodes; this favors safety/coherence over minimum-bias accuracy for strongly heterogeneous softening distributions.
  - Snapshot softening sidecar field is additive and optional; external consumers that assume only canonical GADGET datasets may ignore `GravitySofteningComoving` unless upgraded.

## New open follow-up items (post P29 closeout)

| ID | Status | Area | Current limitation | Next-step evidence target |
|---|---|---|---|---|
| P30-CLUSTERED-LOAD-MATURITY-025 | Open | Distributed gravity clustered-load maturity | Current deterministic contiguous SFC cuts reduce clustered overshoot but remain one-dimensional; extreme multi-cluster anisotropy can still create per-rank remote short-range pressure skew even when weighted load appears balanced. | Add clustered MPI scaling sweeps (2/3/4+ ranks) that jointly track weighted-load imbalance and remote request packet imbalance against acceptance thresholds in `bench_tree_only_scaling_mpi` + `bench_tree_pm_coupling`. |
| P31-STATE-TEST-DRIFT-026 | Closed | Snapshot schema contract drift | `test_unit_snapshot_hdf5_schema` had a stale `gadget_arepo_v2` expectation while runtime schema advertised `gadget_arepo_v4`; the unit contract and docs are now aligned to schema v4. | `test_unit_snapshot_hdf5_schema` rebuilt and passed in the 2026-04-26 CPU-only hardening pass. HDF5 roundtrip remains environment-dependent on HDF5-enabled preset availability. |
| P32-MIGRATION-IDENTITY-027 | Closed | Migration/compaction identity hardening | Migration, stale-ghost removal, ownership compaction, and rank-local gas rebuild needed one explicit identity-preserving contract with negative coverage for wrong-owner inbound records, stale ghosts, and stale sidecar indices. | Hardened `commitParticleMigration` stale-ghost validation; added targeted species migration negative/sidecar tests; documented `ParticleMigrationRecord` completeness and gas-cell rebuild-by-particle-ID implications. |
| P33-TRANSFORM-FUZZ-HARNESS-028 | Closed | Transform invariant hardening | Resize, reorder, compaction, migration, active-view invalidation, sparse optional sidecars, and HDF5 restart paths lacked one fixed-seed composed fuzz harness keyed by particle ID rather than row counts. | Added `integration_transform_fuzz_invariants` and `validation_transform_fuzz_invariants_long` with deterministic seed/trial failure context, all sidecar-family oracles, stale scatter negatives, gas rebuild checks, and HDF5 restart roundtrip when enabled. |
| P34-ACTIVE-VIEW-LIFETIME-029 | Closed | Active-view lifetime/generation hardening | Compact active/kernel views needed explicit transient lifetime documentation, source-generation stamps on read-only views, stale scatter coverage after migration/gas rebuild, and workspace reuse evidence. | Added generation fields to read-only `ParticleActiveView`/`CellActiveView`, kept mutable kernel scatter generation guards, documented non-ownership/hot-lane contracts in `TransientStepWorkspace`, and extended `test_unit_simulation_state` plus `test_integration_reorder_compaction_sidecars` for workspace reuse, migration invalidation, and gas rebuild invalidation. |


## 2026-04-20: Remaining Phase 3 campaign blockers

- Large-rank (beyond np2) strong/weak scaling certification sweep still required for Phase 3 closeout.
- External literature-target calibration for cosmological P(k) and halo profile systematics remains outside this stage; current campaign uses explicit in-repo reference contracts and envelopes.

## 2026-04-21: Phase 3 final integration closeout blockers

- `cmake --preset mpi-hdf5-fftw-debug` fails in this environment because `fftw3_mpi` is unavailable; distributed MPI+FFTW command bundle cannot be executed end-to-end in this cycle.
- `ctest --test-dir build/pm-hdf5-fftw-debug -R "integration_reference_workflow" --output-on-failure` fails with `runtime workflow schema compatibility validation failed`.
- `ctest --test-dir build/pm-hdf5-fftw-debug -R "integration_tree_pm_coupling_periodic" --output-on-failure` fails with `expected at least one low-res contaminant`.
- `validation_phase2_mpi_gravity_single_rank` residual-cutoff expectation mismatch is repaired; cutoff traversal evidence now accepts internal-node pruning or leaf-pair cutoff skips.
- `validation/artifacts/research_grade/phase3/correctness/power_spectrum_consistency.json` remains blocked (`missing_diagnostics_inputs`), so Phase 3 correctness/force-accuracy envelope is not fully evidenced.
- Multi-rank scaling evidence in the audited artifact set remains baseline-only and non-certifying (`phase2_baseline_scaling_summary.json`); checked-in CSV artifacts currently cover `np1` only.

Phase 3 status in this cycle: **incomplete** (see `docs/treepm_phase3_closeout.md`).

## 2026-05-07 restart sidecar exactness status

- Closed P0 blocker for v6 restart sidecar exactness in the HDF5 path: missing softening mask/value datasets and older restart schema versions are rejected, and targeted reorder/migration roundtrip tests cover identity-aligned sidecars.
- Remaining compatibility work, if needed, must be an explicit legacy importer with documented migration semantics and tests; no silent v5 compatibility is allowed.

## 2026-05-08 distributed ownership floor status

- P1 distributed ownership correctness floor is partially hardened in code by rank-local ID uniqueness summaries, reduced count/sum/xor partition checks, and H2-safe gas-cell/AMR patch migration payloads keyed by stable `gas_cell_id`/`patch_id` in the reference workflow MPI path.
- Mature load balancing and pencil FFT remain out of scope for this repair cycle.
- MPI+HDF5+FFTW closure remains dependent on the validation environment providing MPI FFTW (`fftw3_mpi`) and HDF5-enabled MPI presets; if unavailable, the MPI command bundle must be reported as environment-blocked rather than silently passed.


## 2026-05-08 notes

- P2 AMR/moving-mesh gas ownership readiness now has an isolated `GasCellIdentityMap` API guard and production distributed workflow compaction/rebalance uses stable gas-cell/patch migration payloads. Remaining blocker: broader AMR neighbor exchange/reflux and rank-count-changing restart remain separate production-claim gates.

## 2026-05-08 Stage 1 runtime-truth gate status

- Closed CI-discipline gap for CPU-runnable P0 runtime-truth repairs: the P0 suite now has an exact CTest preset, labels, a label/feature-gating audit test, and one local/CI script under `scripts/ci/`.
- Feature-specific closure remains dependency-bound: HDF5, FFTW, MPI, CUDA, and Python presets are explicit and should fail at configure time when requested dependencies are absent; CPU-only Stage 1 success must not be reported as feature-path closure.


## Stage 2 timestep ownership follow-ups (opened 2026-05-10)

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P35-STAGE2-TIMESTEP-OWNERSHIP-032 | Closed | Hierarchical timestep authority | Restart writer/hash/read now reject stale particle time-bin mirrors against scheduler `bin_index`, reader rebuilds mirrors from scheduler authority after validation, migration/transfer fields are documented as mirrors, and reorder/migration tests exercise stable-ID mirror preservation. | Remaining scheduler ownership work moved to P36 for new particle registration; no additional restart mirror guard is open for existing v6 lanes. |

## 2026-05-10 Scheduler ownership residuals

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P36-SCHEDULER-SPAWN-OWNERSHIP-033 | Open | Particle creation + scheduler authority | Star-formation and black-hole creation no longer copy cell mirror bins into new particle mirror lanes, but production creation paths still need an explicit scheduler-owned element registration policy before spawned particles can participate in hierarchical activation without a rebuild. | Add scheduler element-registration/import APIs with tests covering spawned particles, mirror refresh, and restart continuity. |

## 2026-05-10 hierarchical timestep invariant status

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P37-HIERARCHICAL-TIMESTEP-INVARIANTS-034 | Closed | Scheduler/debug invariants | Runtime checks now trap early activation, illegal sync-boundary bin transitions, skipped PM refresh commits, stale mirrors/descriptors, invalid restart timestep state, active-set mismatches, and non-monotonic tick overflow with deterministic context. | Feature-path closure still depends on HDF5/FFTW/MPI presets in environments with those dependencies available; CPU-only surrogate scheduler tests are in `test_unit_time_integration`. |


## 2026-05-10 scheduler identity exchange follow-up

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P38-SCHEDULER-MPI-IDENTITY-035 | Partially closed for gas-cell H2 migration | MPI migration + scheduler authority | Particle scheduler identity exchange exists for particle migration, and H2 gas-cell/AMR patch migration now exchanges gas-cell scheduler identity records by `gas_cell_id`. Spawned-particle registration and rank-count-changing restart are still outside this closure. | Run the new two-rank MPI+HDF5 test in an environment with `MPI_CXX`; keep P36 spawned-particle scheduler registration open before broader exact-continuation claims. |

## 2026-05-11 Stage 2 verification status

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P39-STAGE2-GATE-INCLUSION-036 | Closed | CI/runtime-truth gates | Stage 2 scheduler-authority tests existed but the strongest CPU runtime-truth preset did not explicitly include `unit_time_integration`, and there was no Stage 2-named preset for scheduler authority, split-brain removal, PM sync legality, and restart equivalence. | `unit_time_integration` is now marked `runtime_truth;p0` with scheduler/restart/treepm labels, the label audit requires it, Stage 0/1 runtime-truth presets include it, and `test-stage2-runtime-truth-cpu-debug` names the Stage 2 critical gate directly. |
| P40-STAGE2-HDF5-RESTART-NEGATIVE-037 | Closed | Restart schema validation | Restart HDF5 roundtrip covered stale particle `time_bin` mirrors and payload-hash tamper, but did not corrupt a persisted scheduler hot lane directly. | `integration_restart_checkpoint_roundtrip` now mutates `/scheduler/active_flag` to an invalid value in-file and requires `readRestartCheckpointHdf5` to reject it before exact continuation is accepted. |


## 2026-05-11 Stage 2 documentation/release anti-drift status

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P41-STAGE2-DOC-RELEASE-ANTI-DRIFT-038 | Closed | Docs/release claim discipline | Stage 2 timestep authority docs, restart/schema language, runtime truth map, ADR, repair-state docs, release notes, and docs scaffold now agree that scheduler objects own timestep truth and PM cadence legality while `time_bin` lanes are mirrors. | Keep `integration_docs_scaffold` green whenever changing Stage 2/Phase 3 language; do not claim Phase 3 multirate TreePM production closure without the hard-gate evidence in `docs/treepm_phase3_contract.md`. |

Command evidence for the documentation synchronization patch:

- `cmake --preset cpu-only-debug` → pass.
- `cmake --build --preset build-cpu-debug` → pass.
- `ctest --preset test-stage1-runtime-truth-cpu-debug --output-on-failure` → pass (17/17).
- `ctest --preset test-cpu-debug --output-on-failure` → pass (77/77).
- `cmake --preset hdf5-debug`, `cmake --build --preset build-hdf5-debug`, `ctest --preset test-hdf5-debug --output-on-failure` → pass (78/78).
- `cmake --preset pm-hdf5-fftw-debug`, `cmake --build --preset build-pm-hdf5-fftw-debug`, `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` → pass (78/78).
- `cmake --preset mpi-hdf5-fftw-debug` → blocked by missing `fftw3_mpi`; no MPI+HDF5+FFTW closure is claimed.

Remaining blockers are unchanged: P36 spawned-particle scheduler registration, P38 distributed scheduler identity exchange, and dependency-enabled MPI+HDF5+FFTW feature-path evidence are still required before broader closure claims.

## 2026-05-19 PM refresh stage-diagnostic update

No issue status changes in this patch. Open blockers remain P38 distributed scheduler identity exchange and dependency-enabled MPI+HDF5+FFTW feature-path evidence.

- S5-02 follow-up: monitor downstream configs still using legacy aliases and migrate examples to canonical keys.

| S5-04-PROVENANCE-VERSIONING | Closed | Provenance/schema auditability | Provenance lacked explicit config schema identity and separate raw/normalized/derived runtime payload capture in snapshot/restart metadata. | Added `provenance_v5` fields for config schema identity, normalized hash, and audit payloads with deterministic serialization/roundtrip coverage updates. |

## 2026-05-29 Stage 8 restart-boundary status

No existing repair issue is closed by this patch. The new boundary contract blocks unsafe restart writes but does not implement intentionally represented half-step/local-substep restart; those states remain unsupported and must continue to fail loudly unless a future schema-versioned design serializes all required scheduler/phase state.

## 2026-06-29 MPI TreePM collective-order and patch-topology repair

| ID | Status | Area | Finding | Resolution / Remaining action |
| --- | --- | --- | --- | --- |
| MPI-TREEPM-COLLECTIVE-ORDER-20260629 | Pending MPI validation | Distributed TreePM short-range transport | A rank-local zero-byte early continue could skip request/response `MPI_Alltoallv` calls after the count exchange while peers still entered the collective. | Repaired so every rank enters request and response `MPI_Alltoallv` for each globally coordinated batch, with zero counts and safe non-null empty buffers. Added byte-count/displacement overflow/alignment guards and repaired the response max-peer diagnostic. Requires the registered 2-rank and 3-rank MPI gates. |
| MPI-GAS-PATCH-REBUILD-20260629 | Pending MPI validation | Workflow ownership compaction | Rebuilding local gas state fabricated one patch (`patch_id=1`, level 0) and discarded retained AMR patch topology. | Rebuilt state now carries stable patch ID and patch-local order in migration records, preserves patch metadata/ranges/ownership, remaps dense `patch_index`, and fails loudly for missing metadata or an unsupported split geometry-bearing patch. |
| MPI-PM-SAME-PEER-HALO-20260629 | Pending MPI validation | Two-rank periodic PM halo contract | The same-peer (`left_peer == right_peer`) Sendrecv orientation had no direct regression. | Added `integration_pm_slab_halo_exchange_mpi_two_rank`, using distinct boundary-plane encodings to verify that left halo receives peer right edge and right halo receives peer left edge. |

## 2026-07-13 gravity hardening status (supersedes pending gravity rows)

| ID | Status | Area | Current evidence or limitation | Required follow-up |
| --- | --- | --- | --- | --- |
| MPI-TREEPM-COLLECTIVE-ORDER-20260629 | Closed for np1--np4 correctness | Distributed TreePM collective order | Periodic/Phase-2 focused MPI matrix passes with empty, source-only, target-only, all-empty, repeated-epoch, and migration-to-empty cases. Every global request/response batch retains zero-count participation. | Keep np1--np4 regression coverage; scaling transport is tracked separately below. |
| MPI-PM-SAME-PEER-HALO-20260629 | Closed | Periodic PM halo | Halo exchange now uses global-plane owner lookup plus receive-first nonblocking completion. np2/np3/np4 pass same-peer, zero-width-slab, and single-plane cases. | Preserve tests when changing PM layout or tags. |
| GRAVITY-PERIODIC-GEOMETRY-20260713 | Closed for covered correctness matrix | Periodic tree topology and moments | Largest-gap per-axis frame, rectangular minimum images, compact roots, corrected multipoles, remote bounds, and integer-image translation tests pass. | Extend only when adding a new tree representation or higher multipole order. |
| GRAVITY-PACKET-INTEGRITY-20260713 | Closed for fail-fast MPI contract | PM and TreePM wire protocols | Versioned little-endian PM density/force/potential records and Tree hierarchy/request/response records carry runtime identity, finite/bounds/overflow validation, and exact response coverage. Actual-MPI-world preflights coordinate solve metadata, PM halo-cache commit, compact active/PM/zoom buffers, periodic unwrap/local-tree construction, worst-case traversal-stack reservation, reusable exchange-workspace growth, request preparation, and protocol decoding so zero/source-only/target-only ranks retain uniform collective order through four ranks. Workspace capacities are memory-accounted. | No process-failure recovery is claimed; retain MPI-abort as the fatal boundary unless fault tolerance becomes a separate project. |
| GRAVITY-EWAL-TSC-PROFILE-20260713 | Closed for named small-N profile | Periodic force accuracy | Independent Ewald reference converges. TSC+deconvolution with `asmth=1.25`, `rcut=6.25`, runtime-default `theta=0.7` passes `relative_L2<=1e-2`, `p99<=5e-2`; worst values are `8.059667017e-3` and `8.157551429e-3` on the exact cancellation-dominated `4^3` DMO initial lattice. Its Ewald/PM/tree/total L2 values are `10.16780756` / `8.296819318` / `1.789203432` / `10.08586460`. The independent-target calibrated MAC matrix asserts all three MACs against Ewald. | Any mesh/split/cutoff/softening/MAC profile change needs a new envelope. |
| GRAVITY-CIC-ACCURACY-20260713 | Open scientific profile | CIC periodic force accuracy | CIC remains functional, but covered uniform and split-transition cases exceed the certified target. The historical 4.5-cell cutoff also shows roughly nine-percent cutoff-transition error. | Calibrate a distinct CIC/lower-cutoff profile or keep it explicitly diagnostic; do not weaken the TSC gate. |
| GRAVITY-CADENCE-20260713 | Closed by restriction | PM time consistency | Production config rejects cadence other than one and refreshes every rank-coordinated force surface. Lower-level reuse requires an exact coordinator cache signature; missing/incompatible reuse fails coherently rather than silently changing into refresh, and mixed explicit votes fail before collectives. | A future cadence `>1` still requires a validated predictor/interpolator and restart-tested temporal semantics; the invalidation signature is not that predictor. |
| GRAVITY-HIERARCHICAL-KDK-20260713 | Closed by fail-closed restriction; feature remains open | Mixed-rung time integration | Production config and restart require `hierarchical_max_rung=0`/bin-zero schedulers because current workflow state lacks per-element kick/drift epochs. This prevents elapsed-time-skipping mixed-rung execution. | Implement per-element KDK epochs, source prediction, synchronized PM semantics, restart state, and direct-vs-resumed mixed-rung tests before enabling nonzero rungs. |
| GRAVITY-DECOMPOSITION-EPOCH-20260713 | Closed for current workflow | Ownership/cache generation | Workflow epoch advances only after actual committed particle ownership change, is persisted/restored, and guards tree wire records. Fixed-slab PM compatibility deliberately excludes it; dense-row force history invalidates immediately. | Rank-count-changing restart/topology remap remains unsupported. |
| GRAVITY-GAS-COSMOLOGY-NORMALIZATION-20260713 | Closed for implemented source contract | Cosmological gas gravity/expansion | Gas momentum/work now consumes scale-free TreePM `A` as `A/a^2`; energy expansion is `-H(2K+3P)`. Fixed-grid and AMR callbacks use the post-drift timeline's step-end `a,H_code`, not still-step-begin committed state or SI `H`. | Add a dedicated expanding hydrodynamic cosmology validation campaign before a publication-grade baryonic claim. |
| GRAVITY-TIMESTEP-COORDINATES-20260713 | Closed for implemented criterion | Gravity timestep normalization | Force lanes store scale-free `A`; post-step criteria use committed `a` to form comoving-coordinate acceleration `|A|/a^3`, giving `dt_grav=eta sqrt(a^3 epsilon_com/|A|)` for comoving softening. | Preserve the scale-factor regression when enabling nonzero rungs or changing force-variable conventions. |
| GRAVITY-LET-SCALING-20260713 | Open scalability blocker | Distributed short-range tree | Current hierarchy `MPI_Allgatherv` plus target/response `MPI_Alltoallv` is correct through np4 but is not a locally essential tree. | Implement selective remote-branch/node import with bounded memory, epochs, equivalence tests, and communication-volume scaling evidence. |
| GRAVITY-PM-SPARSE-PEERS-20260713 | Open scalability blocker | Distributed PM records | Density/interpolation records are sparse, but count/data completion is communicator-wide. | Add a cached epoch-invalidated sparse-peer plan while retaining the collective path as a reference until equivalence is proven. |
| GRAVITY-DMO-POWER-SPECTRUM-20260713 | Closed for controlled DMO gate | DMO observables | Superseding the earlier missing-artifact assessment, the np1--np4 workflow emits a versioned binned 3D `P(k)` JSON. Direct and power growth increments meet the expected-increment and cross-estimator gates. The exact initial `4^3` cancellation-dominated lattice is separately certified against periodic Ewald; its short-step trajectory requires response/ballistic RMS and growing-mode projection in `[0.01,0.03]` with alignment at least `0.99`, preventing seeded velocity from masking a wrong force. | Keep this deterministic direct-DFT fixture as the correctness reference when adding a performant large-mesh estimator; closure does not imply publication-grade DMO validation. |
| GRAVITY-ISOLATED-OPEN-20260713 | Open scalability/scientific profile | Isolated and zoom PM | Assignment/gather clips rather than wraps. Typed config rejects window deconvolution for isolated/open gravity; tracked isolated decks explicitly set it false, with no silent runtime override. Distributed open PM is bounded root gather/scatter; focused zoom all-gathers high-resolution sources and forces unsupported isolated window deconvolution off. | Add scalable open-boundary/zoom transport and an independent isolated/zoom accuracy envelope before broader production claims. |
| GRAVITY-LARGE-DMO-20260713 | Open scientific/scaling blocker | Production cosmology | No large-volume resolution/convergence, halo, long-horizon growth, cross-code, or beyond-np4 campaign is certified. The small deterministic `P(k)` artifact is now present. | Run a documented large-volume science and scaling campaign; compare any optimized spectrum estimator against the closed direct-DFT reference. |
| GRAVITY-GPU-PARITY-20260713 | Open performance/accuracy blocker | GPU PM | CUDA remains CIC-only and is outside the certified TSC profile. | Add TSC assignment/gather parity and independent force-accuracy evidence before a production GPU claim. |

Supersession note (2026-07-13): the 2026-04-21
`power_spectrum_consistency.json` `missing_diagnostics_inputs` entry above is a
preserved historical Phase 3 campaign observation. It does not describe the
current production-workflow DMO gate, whose versioned binned 3D `P(k)` JSON is
now implemented and validated as recorded in the closed row above.

The current readiness matrix and the distinction between controlled DMO
validation and publishable production are maintained in
`docs/gravity_production_readiness.md`.


## 2026-07-17 Campaign A capability follow-up

| ID | Status | Area | Resolution | Remaining boundary |
|---|---|---|---|---|
| P40-RUNTIME-RESOURCE-CAPABILITY-037 | Closed | Workflow resource ownership | Removed public `StepContext`/`ownerContext`/`stageContext` escape, bound descriptor resources per task, enforced typed-view access allowlists, added negative compile tests, and centralized gas-cell ownership helpers. | Built-in owners remain trusted source implementations; dynamic lane tracing and descriptor-owned timestep/restart/migration schemas are not claimed. |

## Campaign B external acceptance handoff (2026-07-24)

No known source-level Campaign B acceptance item is intentionally left open.
The remaining Campaign B gate is execution evidence on a dependency-complete MPI
system: configure/build the `mpi-hdf5-fftw-debug` preset and run the complete
distributed IC matrix. Keep `distributed_ic_import` provisional until that
matrix passes. Any failure must be repaired in code or protocol; do not promote
capability by changing the status test, disabling a fault scenario, or increasing
its timeout.

Supersession note (2026-09-04): the repository-wide runtime smoke issue is now
diagnosed. Production output writes a logical snapshot set under
`snapdir_###`, while the old smoke asserted an obsolete flat snapshot path.
The P2/P3 post-audit closure aligns both serial and MPI smoke checks with the
snapshot-set completion-manifest contract and keeps the production smoke gated
on FFTW. Runtime re-execution of that gate still requires a dependency-complete
FFTW environment; the historical Campaign B ingestion result is otherwise
unchanged.

## Campaign B completion follow-ups (2026-07-25)

| ID | Status | Area | Resolution in this pass | Remaining boundary |
|---|---|---|---|---|
| CAMPAIGN-B-MPI-ACCEPTANCE-20260725 | Environment-blocked | Distributed IC runtime | Persistent HDF5 sessions, bounded multi-chunk batches, route/accounting counters, and the np1/np2/np4 source test matrix are present. | Configure/build/run `mpi-hdf5-fftw-debug` on a dependency-complete machine; do not promote `distributed_ic_import` until it passes. |
| CAMPAIGN-B-COLLECTIVE-GUARDS-20260725 | Environment-blocked runtime evidence | MPI failure consensus | Final acceptance guards manifest digest, reader-session creation, serialization, post-exchange, source/final, duplicate-audit, identity-revalidation, append, and final-state phases; rank-specific finite-time fault tests are registered and MPI-enabled syntax checks pass. | Execute the registered fault matrix under real MPI; no source-level seam is intentionally waived. |
| CAMPAIGN-B-INITIAL-REBALANCE-20260725 | Open | Distributed IC ownership | Deterministic x-slab ingestion ownership is documented and existing runtime rebalance ownership was preserved. | Measure imported work imbalance and invoke the existing legal rebalance path before the first expensive production phase when a configured threshold is exceeded; validate under clustered and zoom ICs. |
| CAMPAIGN-B-IC-SOURCE-SPLIT-20260725 | Closed | IC maintainability | The 4,136-line monolith was removed and replaced by responsibility-based manifest/schema, serial, streaming, and distributed units; every IC source is below 2,000 lines and structural CI rejects reintroduction or private-header leakage. | Preserve the structural gate as the subsystem evolves. |

## 2026-08-08 `src/core` hardening follow-ups

| ID | Status | Area | Resolution in this pass | Remaining boundary |
| --- | --- | --- | --- | --- |
| CORE-SCHEDULER-ORDERED-ACTIVE-20260808 | Closed with implementation + release benchmark evidence | Hierarchical scheduler | Deterministic active ordering now uses a reusable four-pass 8-bit radix sort over collected `uint32_t` local indices, preserving O(1) swap-with-last bin updates and avoiding a full-N scan. Release benchmarks on 100k/500k and fine-bin-heavy cases reduced `beginSubstep` time by roughly 42-54% versus the prior comparison-sort baseline; the isolated comparison-sort reference costs 1.64-1.92x the complete radix-based begin-substep path. | Rebenchmark on future representative production hardware when scheduler structure changes; no open correctness/performance blocker remains from CORE-012. |
| CORE-MIGRATION-PEAK-RSS-20260808 | Closed with implementation + measured envelope | Migration memory | Removed repeated O(N) dense-identity validation from each row, retained one operation-boundary validation, and removed the redundant full `GasCellMigrationFields` staging copy from gas-cell commit while preserving prepare/validate/commit. At 50k cells measured peak/steady RSS is ~1.90x mostly-local, ~2.04x moderate, and ~2.31x heavy migration, including exchange payloads. Moderate migration commit time scales from ~7.5 ms at 5k cells to ~100 ms at 50k and ~160 ms at 100k on this audit host after removing the accidental per-row O(N) validation. | Keep the benchmark as a regression/evidence tool; chunked destructive migration is not justified unless larger production cases exceed the documented workstation envelope. |
| CORE-PROVENANCE-V7-INTEGRITY-20260808 | Closed | Provenance | Current writers emit `provenance_v7` with a real SHA-256 normalized-config digest, strict reversible text escaping, compiler flags, CPU/RAM/host and MPI/CUDA descriptor fields. Valid v6 records remain readable and malformed v6/v7 escapes fail closed. SHA-256 is covered by standard known vectors and v7 round-trip/tamper tests. | MPI/CUDA runtime values remain environment-dependent and continue under the separate runtime-evidence issue; schema integrity itself is closed. |
| CORE-MPI-CUDA-RUNTIME-EVIDENCE-20260808 | Environment-blocked | Distributed failure/CUDA mapping | Source now uses node-local MPI rank for GPU selection, preserves CUDA discovery diagnostics, and aborts distributed harness jobs on rank-local fatal failure. CPU builds of the affected core/parallel/harness paths pass. | Execute rank-selective failure and CUDA node-local mapping tests on a dependency-complete MPI/CUDA system before claiming runtime closure. |
| CORE-METADATA-TEXT-CODEC-20260809 | Closed | Metadata persistence | `StateMetadata` v3 uses a strict reversible line codec for newline, CR, tab, backslash, equals, control bytes and UTF-8 bytes; valid v2 remains readable. Unknown/dangling/truncated escapes fail closed and hostile-string round trips are tested. | No residual CORE-018 serialization gap. |
| CORE-LOCAL-INDEX-CAPACITY-20260809 | Closed | Local index capacity | Central checked local-count/row helpers now cover authoritative state resize, migration/rebuild destinations, patch cell counts, scheduler reorder rows, and checked count addition before narrowing. Raw core SoA resize paths enforce the same capacity. | `UINT32_MAX` remains the invalid gas-cell-row sentinel; highest valid row is `UINT32_MAX-1`. |
| CORE-MEMORY-ACCOUNT-CALIBRATION-20260809 | Closed with measured evidence | Memory accounting | Controlled 200k-particle/100k-cell calibration reports 34,198,583 known-owned bytes versus a representative ~33.97 MB RSS delta on this Linux audit host (known-owned/RSS-delta ratio ~1.007); unknown external allocations stay explicitly labelled rather than folded into false precision. | RSS remains a process diagnostic, not an exact ownership total; rerun on target platforms as part of performance qualification. |


## 2026-08-29 M1B/M1C memory-governance handoff

| ID | Status | Area | Resolution in M1B | Remaining boundary |
| --- | --- | --- | --- | --- |
| MEMORY-GOVERNOR-FOUNDATION-20260829 | Closed for M1B foundation | Process memory admission | One process governor, move-only reservations, pressure/headroom/high-water telemetry, report-derived scheduler baseline reconciliation, production budget ownership, and a real production all-particle gravity index allocation through governed monotonic scratch are implemented and covered by focused/production smoke tests. Scratch admission occurs only after collective process preflight and rank-local failures are coordinated before later TreePM collectives. | M1B does not claim broad large-allocation coverage or whole-process measured-memory closure. |
| MEMORY-GOVERNOR-BROAD-INTEGRATION-20260829 | Open for M1C | Runtime large-allocation coverage | The M1B authority and reservation contract now exist. | Extend coordinated reservations to PM/FFT, tree/LET and communication staging, migration, output/restart/readback, analysis/halo work, hydro/AMR phase batches, and other target-scale allocations as justified. Preserve MPI collective safety for rank-dependent rejection. |
| MEMORY-RSS-RECONCILIATION-20260829 | Open for M1C | Acceptance telemetry | Owned-capacity and governor accounting are deterministic and reconciled logically. | Add RSS/PSS/peak-RSS observation and rank max/mean reconciliation on dependency-complete production systems; diagnose measured-versus-declared discrepancies without treating RSS as exact ownership truth. |

## 2026-09-03 M1C-2 memory-acceptance boundary

| ID | Status | Area | Closure / remaining boundary |
| --- | --- | --- | --- |
| MEMORY-M1-SOURCE-ACCEPTANCE-20260903 | Closed | M1 memory architecture | All 12 mandatory M1 source/architecture gates are closed in `docs/repair/m1c2_m1_acceptance_closure_20260903.md`; the 512^3/8-rank PM routing expression is explicitly regression-tested. |
| MEMORY-M1-DISTRIBUTED-RUNTIME-20260903 | Environment-blocked | M1 distributed runtime validation | Multi-rank PM/TreePM/rank-imbalance/restart equivalence must be reported only from an MPI+FFTW-capable environment. Missing dependencies are an environment blocker, not permission to weaken source gates. |
| MEMORY-M2-FULL-PHYSICS-20260903 | Open for M2 | Hydro/AMR/subgrid/device memory architecture | Hydro active batching, sparse AMR/regrid coexistence, chemistry fidelity tiers, species-specific subgrid sidecars, large analysis products, and device budgets remain the explicit next campaign. |

## 2026-09-05 M2C full-physics memory handoff

| ID | Status | Area | Closure / remaining boundary |
| --- | --- | --- | --- |
| MEMORY-M2C-SPECIES-SIDECARS-20260905 | Closed for released species | Full-physics persistent state | Gas/star/BH/tracer state remains species-local; DMO/gas-only zero-cold-sidecar allocation is regression-tested. Dust is not implemented and allocates nothing. |
| MEMORY-M2C-FEEDBACK-EVENTS-20260905 | Closed | Stellar feedback memory | Duplicate star carry mirrors and full-star delta staging are removed; feedback uses capped event batches, exact indexed nearest-k target selection, on-demand AMR cell volume, and O(k) per-event target scratch. |
| MEMORY-M2C-MICROPHYSICS-TABLES-20260905 | Closed for effective ISM table | Shared microphysics | Reference hydro/source/timestep consumers share one immutable EOS table. Preserve shared/read-only ownership for future tables unless measured thread-local replication is justified. |
| MEMORY-M2C-CHEMISTRY-TIERS-20260905 | Released scalar tier closed; higher tier fail-closed | Chemistry/enrichment fidelity | `primordial`/`primordial_metal_line` cooling and `metal_species_mode=total_only` remain explicit typed choices. `core_elements` remains rejected until the complete element-sidecar/yield/cooling/AMR/MPI/restart contract exists. |
| MEMORY-M3-GLOBAL-SCHEDULING-20260905 | Closed for M2D deterministic source/CPU scope; distributed runtime evidence environment-blocked | Whole-runtime resource scheduling | M2D records the major-task dependency/lifetime/pressure graph in the existing registry, derives deterministic governor-aware hydro/communication/feedback batch sizes, enforces hard rank-memory decomposition ceilings, publishes node/rank/thread topology, and catches up optional science diagnostics after transient Red/Trip deferral without forcing them through unsafe pressure. The current serial dispatcher preserves accepted scientific order. All 136 executable CPU tests outside the separately long source-package-completeness regression are green on the final tree; focused HDF5 M2D/reference-workflow tests are 4/4. MPI C++ and FFTW development dependencies are unavailable here, so multi-rank topology/scientific evidence remains an environment handoff. Predictive controllers, measured selective overlap, GPU/device-residency scheduling, and cross-rank shared-memory resources remain later campaigns. |

## 2026-09-05 M2C post-merge transient-memory closure

| ID | Status | Area | Closure / remaining boundary |
| --- | --- | --- | --- |
| MEMORY-M2C-SF-TRANSIENT-BATCHING-20260905 | Closed | Star-formation source memory | Removed all-active full-cell index materialization and population-scale input staging from production `SourceRuntime`; input work is <=4096 cells per globally coordinated batch and capacity growth is governor-admitted. |
| MEMORY-M2C-DIFFUSION-TRANSIENT-20260905 | Closed for single-rank released diffusion | Metal/enrichment source memory | Production diffusion now uses canonical gas spans plus 49 B/cell compact phase state and reusable five-scalar workspace; the 128 B/cell mirror, separate full-cell volume lane, duplicate face copy, and population-scale neighbor-accumulator graph are removed. Distributed conservative diffusion remains intentionally fail-closed until its remote flux contract exists. |
| MEMORY-M2C-SOURCE-ACCOUNTING-20260905 | Closed | Runtime memory truth | `SourceRuntime` owner capacities/high-waters are merged into startup/per-step `MemoryReport` and persisted profiler JSON, including SF/feedback/diffusion source scratch. |
| MEMORY-M2C-BH-ACTIVE-STAGING-20260905 | Closed for production accretion | BH/AGN event memory | Empty active-index span denotes all local BH rows, eliminating the O(N_BH) all-active index vector. Production seeding remains fail-closed without an authoritative candidate provider rather than allocating speculative population staging. |

## 2026-09-06 M2C final acceptance closure

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| IO-M2C-CHUI-BOX-UNITS-20260906 | Closed | Snapshot HDF5 validation/readback | CHUI Mpc-valued axis box metadata was compared directly against CHUI-native code-unit coordinates. | Reader and validator now convert box bounds into the matching code/stored unit; regression covers valid 4.0965 kpc and invalid 50.1 kpc coordinates in a 0.05 Mpc box. |
| MEMORY-M2C-DIFFUSION-FACE-GOVERNOR-20260906 | Closed | Metal-diffusion topology memory | The retained O(N_face) diffusion graph was reported but not governor-admitted during replacement. | Face capacity now has a dedicated phase-resident reservation; old/new coexistence is accounted during handoff. |
| MEMORY-M2C-DIFFUSION-HIGH-WATER-20260906 | Closed | SourceRuntime memory truth | Move-replaced diffusion containers could report only current capacity after a later smaller topology. | Explicit historical high-water is retained for `rho*kappa` and face capacities and surfaced in the existing `MemoryReport`. |

## 2026-09-06 M2D-1 remaining acceptance gates

The P1 exact-feasibility and distributed soft-throttle defects have source repairs and focused CPU evidence. Open: dependency-complete MPI+FFTW rank/fault/transaction validation; complete owner-specific major-task incremental peak models before enabling concurrency; full CPU/HDF5/source-package and scientific/topology acceptance. See `docs/repair/m2d1_memory_scheduling_closure_20260906.md`. This historical ledger is not a substitute for CURRENT_STATUS.md.

## 2026-09-07 M2D final admission residuals

The optional heavy-analysis admission defect is repaired in the current
source patch. Remaining acceptance: complete gravity/source/output state-
dependent peak and lifetime contracts; verify external-runtime reserve and
whole-process memory under representative full physics; execute the real
MPI+FFTW np2/np3/np4 safety/scaling/restart matrix and full test inventories.
Nonserial full-physics scheduling and unconditional M2D closure remain blocked
until the evidence is recorded. See the dated M2D final admission repair report.

## 2026-09-07 M2D acceptance-gates residuals

- M2D-F1: source population growth/ID coordination, AMR transaction/ghost
  coexistence, and output writer/readback metadata still require complete
  owner-specific peak models. Current partial models cannot authorize overlap.
- M2D-F2: dependency-complete MPI+HDF5+FFTW np2/np3/np4 hard-memory relief,
  collective preparation faults, transaction coexistence, conservation, and
  restart equivalence remain unexecuted in this repair environment.
- M2D-F3: full CPU/HDF5/source-package, rank max/mean RSS/PSS, work imbalance,
  legal batch equivalence, and Release topology evidence remain acceptance
  requirements. Do not claim a certified 48 GiB full-physics envelope.

The former SFC and heavy-analysis source defects remain repaired. See the new
acceptance-gates report for the exact current commands and outcomes.

## 2026-09-22 — Parallel-HDF5 / MPI diagnostics qualification closure

- HDF5 source/package dependency truth is the C API only; C is explicitly enabled and production Parallel-HDF5 capability is established by MPIO compile/link probes rather than `HDF5_IS_PARALLEL` alone. Production MPI+HDF5+FFTW presets require that capability; serial-HDF5 MPI remains explicitly compatibility-only.
- The mixed-prefix/Conda classifier normalizes separators and cannot abort configure through an invalid backslash regex. `Threads::Threads` is linked to the direct `std::thread` memory-governor test.
- MPI diagnostics now reduce communicator-global counts, boolean invariant failures, angular-momentum vectors, SFR bins, quicklook fields, and distributed memory summaries before one publisher writes shared files. MPI power spectrum is explicitly deferred instead of averaging rank-local spectra.
- Qualification is layered into narrow snapshot, restart, diagnostics, endpoint, real-launcher, feature-summary/runtime-truth, and installed-package-consumer gates. The real 64^3/3-rank first-light rerun remains a dependency/host acceptance gate and is not claimed by source inspection alone.

## 2026-09-22 — 64^3 Parallel-HDF5 smoke closure residual acceptance

| ID | Status | Area | Current state | Remaining acceptance |
| --- | --- | --- | --- | --- |
| FIRSTLIGHT-MPIO-PROBE-20260922 | Source closed; dependency-host evidence required | CMake / Parallel HDF5 | Probe state is no longer preinitialized and stale cache is invalidated before the actual MPIO compile/link check; regression poisons cached state and includes serial-provider fail-closed coverage. | Execute the registered configure regression on a host with MPI + Parallel HDF5 (and serial HDF5 for the negative subcase). |
| FIRSTLIGHT-ZERO-ID-MPI-20260922 | Source closed; real-IC evidence required | Distributed IC ingestion | One typed source-set ID mapping survives fragment assembly, serialization/broadcast, and chunk conversion; warnings are diagnostic only; overflow and shuffled/equivalence fixtures are registered. | Run the registered 1/2/3/4/8 MPI matrix and the exact original zero-based 64^3 MONOFONIC IC. |
| FIRSTLIGHT-SNAPSHOT-NAMESPACE-20260922 | Source patched; Parallel-HDF5 runtime evidence required | Single-file snapshot publication | Rank zero creates/validates `snapshots/`, collective preflight gates HDF5 entry, committed snapshots remain protected, and deterministic path-is-file failure is registered at 2/3 ranks. The 2026-09-26 follow-up additionally requires exact shared-run-directory consensus before collective file creation and makes final-endpoint de-duplication runtime-local so a prior invocation cannot silently satisfy overwrite protection. | Execute the 2/3-rank Parallel-HDF5 snapshot smoke, including divergent-path coordinated failure and committed-namespace rerun refusal; retain host-level unwritable-directory qualification where permissions can be enforced. |

## 2026-09-23 — TreePM observability / geometry-refresh residual gates

| ID | Status | Area | Current state | Remaining acceptance |
| --- | --- | --- | --- | --- |
| TREEPM-P1-COUNTERS-TIMERS-MEM-20260923 | Source closed; build/test deferred by user request | TreePM residual observability + LET memory estimate | Local/incoming residual pair identity, non-overlapping remote-phase timers, workflow tree/PM wall split, shared short-range wire constants, and `estimateTreePmExchangeMemory` feeding the existing MemoryGovernor are implemented with docs. | Run `cmake --preset cpu-only-debug`, `cmake --build --preset build-cpu-debug`, `ctest --preset test-cpu-debug --output-on-failure`, plus MPI/HDF5/FFTW presets when dependencies exist. No production qualification claim from this campaign. |
| TREEPM-P2-GEOMETRY-FRESHNESS-20260923 | Source closed; build/test deferred by user request | Top-domain routing geometry lifecycle | Generation-stamped install, pre-solve O(N) leaf refit, freshness/coverage selection with explicit fallback reasons, and runtime event fields are implemented; geometry remains derived (no restart schema change). | Same build/ctest floor plus at least one two-rank TreePM periodic run exercising drift between ownership updates; confirm fallback reasons never fire on covered sources. |

## 2026-09-23 — TreePM P1+P2 closure residual gates

| ID | Status | Area | Current state | Remaining acceptance |
| --- | --- | --- | --- | --- |
| TREEPM-P1P2-CLOSURE-20260923 | Source closed; build/test deferred by user request | Seed/published geometry split, pure incoming compute timer, ABI asserts, ownership geometry invalidation, attribute cleanup | Stable seed leaf set survives ordinary refits; incoming remote timers split decode/validation vs compute vs encode/pack; `static_assert` enforces host/wire packet identity; ownership commit invalidates geometry freshness; `DomainCoverageResult` attribute placement fixed; retained-capacity high-water uses `vector.capacity()`. | Run `cmake --preset cpu-only-debug`, `cmake --build --preset build-cpu-debug`, `ctest --preset test-cpu-debug --output-on-failure`, plus MPI/HDF5/FFTW presets when dependencies exist; confirm seed SFC identity across two source generations with an empty group and pure-compute timer semantics under the upcoming OpenMP campaign. No production qualification claim from this campaign. |

## 2026-09-23 — TreePM P3 OpenMP residual + P4 compact decomposition residual gates

| ID | Status | Area | Current state | Remaining acceptance |
| --- | --- | --- | --- | --- |
| TREEPM-P3-OPENMP-RESIDUAL-20260923 | Source closed; build/test deferred by user request | TreePM short-range residual OpenMP execution | Immutable shared tree; enforced depth bound; worker-scale stack/counter scratch; deterministic 64-target floating block diagnostics; main-thread softening validation with allocation-free worker lookup; scratch prepared before request posting; actual team-size telemetry across local/distributed/incoming regions; named worker failure capture; serial fallback; docs updated. No MPI from workers; `MPI_THREAD_FUNNELED` only. | Run `cmake --preset cpu-only-debug`, `cmake --build --preset build-cpu-debug`, `ctest --preset test-cpu-debug --output-on-failure`, plus MPI/HDF5/FFTW presets when dependencies exist; confirm deterministic residual counters and observed-worker telemetry across thread counts. Validation/testing intentionally deferred by user request. No production qualification or measured speedup claim is made by this campaign. |
| PARALLEL-P4-COMPACT-DECOMPOSITION-20260923 | Source closed for runtime rebalance, seed, and startup initial-placement paths incl. startup planner admission and dense-index contract; build/test floor deferred to P5 | Large-N decomposition planning + shared source-memory estimate + startup planner admission | Runtime rebalance uses a canonical source view and the source-view compact planner overload; no runtime-rebalance `DecompositionItem[N]` or rich top-domain seed adapter; streaming feedback normalization and transient component recomputation; ≤64 B records; deterministic sort/cuts; exact owner-change intents with conservative intent admission; minimal direct top-domain seed records with construction-to-install reservation transfer; local planner telemetry. Startup initial placement streams compact records through shared `buildMortonSfcDecompositionFromCompact` (dense `local_index`, no stored coordinates/components, pass-2 component accumulation) and now reserves its modeled transient live set (`estimateCompactStartupPlannerTransientBytes`: records, worst-case memory groups, owner/sorted-index vectors, occupancy grids, patch mapping, rank-scaled plan lanes, per-entity lookup scratch) through the single `MemoryGovernor` before any planner population is allocated, releasing at planner-temporary destruction; the rich planner is a reference wrapper over the same core. The compact core enforces `{local_index}` as a permutation of `[0, N)` (`-1` sentinel owner, duplicate rejection, full-coverage verification). Runtime rebalance and top-domain seed admissions share `estimateRuntimeDecompositionSourceStorage` with one gas-incidence condition (DMO admits zero gas bytes); seed admission adds caller-specific seed/leaf memory. `startupParticleMemoryBytes` is the startup per-entity weighting authority; `estimateParticleMemoryBytesForDecomposition` is the source-view species table for storage admissions. | Validation deferred to P5 by user request: duplicate/missing/out-of-range `local_index` rejection, startup ownership equivalence, startup memory-rejection behavior, plus `cmake --preset cpu-only-debug`, `cmake --build --preset build-cpu-debug`, `ctest --preset test-cpu-debug --output-on-failure`, distributed SFC rebalance MPI tests when MPI is available, and confirmation of canonical particle/AMR migration paths, planner byte accounting, and empty rich startup/rebalance call paths. No production qualification claim is made by this campaign. |

## 2026-09-25 — P5-Q integrated P1–P4 qualification results

Full report: `docs/repair/p5q_integrated_p1_p4_qualification_20260925.md`.

| ID | Status | Area | Current state | Remaining acceptance |
| --- | --- | --- | --- | --- |
| P5Q-P1P2P3-TELEMETRY-20260925 | Gates executed green | TreePM observability, geometry freshness, OpenMP worker telemetry | CPU floor 150/150; non-OpenMP floor 150/150 (`-DCOSMOSIM_ENABLE_OPENMP=OFF`, thread matrix correctly narrows to `{1}`); pair identity `total == local + remote` event-by-event, configured/observed worker telemetry equal to the requested team, identical pair totals across thread sets; two-rank Parallel-HDF5 TreePM test green with geometry freshness/fallback fields asserted. | None for the deferred floor. No production speedup claim. |
| P5Q-P4-COMPACT-20260925 | Partially executed; suite-gated items remain open | Compact decomposition + startup admission | Dense-index contract (duplicate/missing/out-of-range, byte-exact messages over 256×24 permutations) and startup transient admission (byte model, overflow, governor reject/accept) registered in `test_parallel_distributed_memory` and green; startup-planner budget rejection asserted at 2 ranks; distributed SFC rebalance 2-rank green. | Startup ownership equivalence and canonical particle/AMR migration paths are still blocked: `integration_distributed_gas_cell_migration_mpi_two_rank` and `integration_distributed_amr_patch_migration_mpi_two_rank` fail identically on baseline `108abfe` (round-limit / AMR wire-count errors); rebalance 3-/4-rank fail identically at baseline (hard-memory collective preparation assert). |
| P5Q-MPI-COMPILE-20260925 | Closed | MPI-only compile of P1–P4 | `persistent_limit` lambda-capture break hoisted out of the `try` in `src/parallel/distributed_memory.cpp`; `mpi-hdf5-fftw-debug` builds RC 0. | None. |
| P5Q-TOP-DOMAIN-ID-COLLISION-20260926 | Source patched; validation pending | Authoritative top-domain pseudo-node identity | Rich, compact-item, and canonical-source top-domain builders previously hashed only owner rank plus SFC begin/end. Tied Morton keys can legitimately split into multiple leaves with identical intervals, deterministically producing duplicate pseudo-node IDs and tripping hierarchy-exchange validation. All builders now also mix the stable within-rank leaf ordinal; refits retain the seed ID. | Re-run phase-2 gravity at 4/8 ranks and tied-key top-domain fixtures; confirm no duplicate pseudo-node identity while rich/compact/source builders remain identity-equivalent. |
| P5Q-COLLECTIVE-SNAPSHOT-PATH-20260926 | Source patched; validation pending | Parallel-HDF5 collective publication / endpoint namespace | Collective single-file publication now broadcasts and compares the normalized absolute shared run directory before HDF5 entry. Final-endpoint de-duplication records snapshots committed by the current runtime instance rather than accepting a `.complete` marker discovered from a previous invocation. | Re-run 2/3-rank snapshot smoke, committed-namespace rerun refusal, divergent per-rank output-root fault case, and endpoint on-/off-cadence coverage. |
| P5Q-MPI-SUITE-20260925 | **Open — pre-existing failure inventory** | Integrated MPI/HDF5/FFTW suite | 232/308 green after repairs; 76 residual failures all reproduce on baseline `108abfe` in this environment (no P1–P4 regressions found). Clusters: distributed IC reader 40, TreePM split/pm-only consistency 5 (bit-identical drift numbers), star-formation MPI split-scale configs 6, gravity-gas identity mirrors 4, pm-periodic routed-density 1, rebalance 3/4-rank 2, migration 2, phase2 gravity 4/8-rank 2, runtime smokes 3, environment probes 2, DMO validation 6 (+1 cascade), launcher 0 (repaired: `chui` mode 100755). `validation_convergence` passes standalone (contention artifact). | Dedicated campaigns per cluster (see report §5). CI equivalence unverifiable without `gh`/network. #11 symptom delta (completed_steps vs baseline mixed-epochs) remains a risk hypothesis. The #276 committed-namespace rerun behavior received the 2026-09-26 source patch above and now requires the dedicated Parallel-HDF5 validation matrix. |

## 2026-10-02 P5-Q2 / M48-Q remaining gates after bounded runtime-validation repair

The current source includes collective TreePM fallback coherence, tied-SFC refit preservation, robust split-scale boundary comparison, and stale distributed-fixture corrections. Remaining campaign blockers are explicit:

- implement/qualify compact conservative periodic-seam top-domain leaf geometry rather than raw box-wide min/max bounds;
- reproduce and diagnose the uninterrupted-vs-resumed DMO digest mismatch before changing any digest definition;
- trace and qualify the MPI restart artifact/cadence contract;
- execute P3 OpenMP OFF/1/2/4 force-vector equivalence and a real MPI+OpenMP incoming-target traversal;
- add/execute the event-by-event P2 healthy -> stale -> collective fallback -> refresh -> selective-routing telemetry sequence;
- run the distributed-IC, PM routing, migration, collective output-root, and fallback regressions under a dependency-complete MPI + Parallel-HDF5 + FFTW toolchain;
- finish the HDF5 provider-family consumer matrix and distinguish serial-HDF5 rejection from toolchain failure;
- triage the CPU-floor failures exposed by the 2026-10-02 full inventory attempt (uniform rung-zero append fixtures, stage-6 memory-report fixture, unique-ID fixture semantics, and diagnostic-governor fixture) before calling the repository qualification-clean.

Do not weaken M48 compact/bounded representations or the restart/identity/memory-governor checks to close these items.

## 2026-11-03 — CHUI Part 2 (source only; NOT RUN — USER-RESERVED)

Record: `docs/repair/chui_part2_thermochemistry_stellar_physics_20261103.md`.

| ID | Status | Remaining acceptance |
|---|---|---|
| CHUI-P2-BUILD | Source written, never compiled | Human build + the Part 2 unit tests listed in the record §11. |
| CHUI-P2-RATES | Part 2.5: provenance checked against the literature. Primary: R01–R12, R14, R16 (restored from GA08), R19, plus KWH96/Wolfire/GP98 cooling. Secondary: R13, R15, R17, R20, R21, H₂ LTE and fine structure. Unverified: R18 and the 0.4 eV photodissociation heat. Corrections: b5, φ_PAH, the Cen term, O I ΔE. Backends are still gated behind `thermochemistry_allow_unqualified_rates`. | User-reserved: run `unit_chemistry_rate_references` and the one-zone validations in `docs/chemistry_rate_provenance.md` §6. Verify R18 against Jura 1975. Owner decision on the equilibrium gate. |
| CHUI-P2-DATA | BLOCKED BY DATA (converters written, NOT RUN) | Sources, checksums and licences are in `tools/chui_data/data_sources.md`. The FG20 TREECOOL and WSS09 have no licence statement, so they are not committed. NuGrid yields are BSD-3 and need tolerance ≈ 0.02. WSS09 (HM01) with FG20 is a declared UVB mismatch, refused unless accepted explicitly. A Ploeckinger & Schaye 2020 converter is blocked. |
| CHUI-P2-TOPOLOGY | Fail closed | Part 3: distributed composition ghosts/reflux, regrid composition transfer, remote feedback deposition. |
| CHUI-P2-GAPS | Open | Element-resolved diffusion (refused by config), science-snapshot export of composition/dust/species/stellar records (not implemented), benchmarks; associations refused. |

## CHUI Part 3 open items (2026-11-04; source only, NOT RUN)

* Recursive >=3-level gas/particle/source multirate (T02) not implemented.
* Multi-rank element diffusion refused (needs per-subcycle remote exchange).
* Derefinement of octets split across ranks skipped; grandparent link may be
  unrecoverable after full sibling refinement.
* Same-level AMR patch faces conserve to roundoff (two-sided evaluation).
* Global step rollback: rung zero only; full-state copy; adaptive-chemistry
  bias of a rejected attempt not reverted.
* Wind/pre-SN channel BLOCKED BY DATA (no chui_ssp_wind_table_v1 bundled).
* Associations grouped per owner rank (decomposition-dependent research mode).
* MPI tests (2/4 ranks) for periodic AMR, distributed composition, remote
  feedback/SF gradients, regrid restart and rollback bit-identity not written.
* See docs/repair/chui_part3_distributed_amr_multirate_20261103.md.

