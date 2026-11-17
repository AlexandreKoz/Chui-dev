# C04 sampled shadow hydrodynamics — final source campaign and pre-compilation handoff (2026-11-17)

**Status: C04 SOURCE IMPLEMENTED — UNVERIFIED (UNCOMPILED / UNTESTED).**
Nothing in this campaign was configured, compiled, executed, benchmarked or
scientifically validated; those stages are owner-reserved. Reported values
are sampled numerical disagreements, never truncation-error bounds or
physical uncertainty.

* Branch `claude/gallant-curie-dm9r9v`, HEAD
  `8317141468b200c5aa3c9192dcfe3aaabd48962e` (merge of PR #140, which already
  contains the revision-6 C04 campaign `2b0ccfe`). The working tree was clean
  at start; this campaign is an uncommitted working-tree delta on top of
  HEAD. Nothing was reset, reverted, committed or pushed.
* Mode (AGENTS.md §2): repair/feature completion of the mandatory C04
  contract. No new config key, no snapshot/restart/provenance schema change,
  no physics model, no `docs/research/shai_hulud/**` change.
* Predecessor ledger: `chui_c04_rev5_cloud_source_closure.md` (revision 6).
  This document supersedes its §5 items 5 and the "rejected attempts are not
  reported" behaviour; everything else in it remains valid.

## 1. Current-HEAD gap table (written before editing; guided the patch)

Revision 6 (HEAD) had already integrated most of what the revision-5 audit
listed as open. The independent re-audit of HEAD against `00_PROMPT.md` found
the following residual **definite source defects**; each is now fixed.

| # | Requirement | Current function / caller at HEAD | Missing integration | Owner | Memory / MPI consequence |
|---|---|---|---|---|---|
| G1 | C04-05 "optional shadow failure must not veto production" | `HydroShadowReference::capture()` allocated provenance/cone scratch outside its try; orchestrator built provenance/ghost-row maps unguarded (`amr_hydro_orchestrator.cpp` capture site); `observeProduction()` grew `m_pending` unguarded; `setGlobalCandidates()` could throw on one rank inside the rank-uniform selection | a local `bad_alloc` escaped into the production sweep (step failure) or, in selection, into a rank-asymmetric path | `HydroShadowReference` | rank-local exception → possible one-rank abort while peers enter the next collective |
| G2 | C04-02 private composition must follow the production rule | `ShadowCompositionSink::apply()` added deltas without the production clamp and without the per-row mass reconciliation of `CompositionFluxAccumulator::commit()` | private composition drifted from the production rule (density-floor mass, negative lanes) | `ShadowCompositionSink` | none (private rows only) |
| G3 | C04-02/03 composition across MPI interfaces (synchronized, the distributed flagship path) | remote ghost donors always fell back to receiver fractions (`composition_remote_donor_faces`) and were cone-excluded, although production assembles the exact step-start remote donor fractions right after the sweep (`completeDistributedCompositionTransport`) | no private use of the production donor table; no hold of composition samples until it exists | `HydroShadowReference`, `HydroAmrRuntimeImpl` | uses the already-exchanged production table: **no new collective, no new message** |
| G4 | C04-02/06 rejected attempts must be counted, never published | `beginEpoch()` reset the report of a rejected attempt before publication, so `discarded_rejected_attempt_samples` could never be non-zero | rejection accounting lost | `HydroAmrRuntimeImpl` | two rank-local counters; values published by the existing sum reductions |

Re-verified as **already integrated at HEAD** (no change needed): production
numerical policy forwarding (C04-01); endpoint ghosts from the production
filler for same-level/diagonal/coarse-fine/closure-mapped/MPI 8³ AMR in
synchronized and recursive stepping; production active sets (partial
activity, suppressed duplicate periodic faces); private refined Strang
chronology with the production `BaryonicSourceStage` on a private
`SimulationState` (fields read by the stage, ISM regime and thermochemistry
operator are all populated — `cells.mass_code`, `gas_cells.{density,velocity,
internal_energy,temperature,metal_mass,gas_cell_id}`, composition sidecar);
physical region keys and one global K; node-aggregate admission (every
`MemoryGovernor::reserve()` consults the attached `NodeMemoryPool`,
`memory_governor.cpp:457`); acceptance-barrier publication with continuity
check against global-step rollback (`time_coordinator.cpp` rollback re-runs
the same step index → discarded). Production installs only
`CoolingHeatingSource` and `EffectiveIsmEnergyRelaxationSource` as hydro
sources (`hydro_amr_runtime.cpp:1273–1279, 4025–4031`), both clonable; the
non-clonable `ComovingGravityExpansionSource` is test-only, so flagship
samples are not systematically refused.

Classification of non-defects (not invented as requirements): multi-rank
fixed-Cartesian hydro is refused by production itself
(`hydro_amr_runtime.cpp:1120` throws); distributed uniform grids run as
level-0 AMR patches, whose region keys are lattice-based and decomposition
independent (MPI fixture: IDs and owners depend on world size, keys do not).

## 2. Changes in this campaign

| File | Change |
|---|---|
| `include/cosmosim/hydro/hydro_shadow_reference.hpp` | `capture()` (both), `observeProduction()`, `setGlobalCandidates()` now `noexcept`; new `provideRemoteDonorFractions()`, `unpublishedResultCount()`, `noteLocalResourceRefusal()`, `finalizePendingEndpoints(bool include_awaiting_remote_donors = true)`; capture spec gains `composition_ghost_gas_cell_ids`, `remote_donor_fractions_will_be_provided`; report gains `composition_remote_donor_slots_resolved`, `composition_reconciled_mass_code`; private `captureImpl()` |
| `src/hydro/hydro_shadow_reference.cpp` | G1 containment wrapper; G2 clamp + reconciliation in `ShadowCompositionSink::apply()` (hydro-only and chronology halfsteps); G3 ghost origin ids, remote-slot count, hold in `observeProduction()`/`finalizePendingEndpoints()`, table resolution; G4 `unpublishedResultCount()` |
| `include/cosmosim/amr/amr_hydro_orchestrator.hpp` | appended `ProductionAmrHydroOptions::hydro_shadow_remote_donor_fractions` (default false) |
| `src/amr/amr_hydro_orchestrator.cpp` | synchronized sweep capture site: preparation guarded (counted refusal), ghost origin ids + promise flag forwarded |
| `src/workflows/hydro_amr_runtime.cpp` | flag set only on the synchronized distributed composition path; endpoint finalize holds remote-donor samples; production donor table handed to the shadow right after `setRemoteDonorFractions`; forced finalize after the composition exchange and before the `kPostHalfBegin` record; rejected-attempt carry; publication keys `discarded_rejected_attempts`, `composition_remote_donor_slots_resolved`, `composition_reconciled_mass_code`, extended `composition_transport` label |
| `tests/unit/test_hydro_shadow_c04_final.cpp` (new) | fixtures for G1–G4 (NOT RUN) |
| `CMakeLists.txt` | registers `unit_hydro_shadow_c04_final` |
| docs | this file; superseding notes in `chui_c04_rev5_cloud_source_closure.md`, `chui_post_patch58_final_source_closure.md`, `docs/repair_open_issues.md`, `docs/hydro_numerics.md` |

Public-interface migration (AGENTS.md §9): all additions are appended or
defaulted; existing callers compile unchanged. `noexcept` on capture/observe
is a strengthening; no caller relied on an exception (verified by search of
`src/` and `tests/`).

## 3. C04-01 … C04-06 source-evidence matrix

| ID | Verdict | Production caller → private owner (file:line at this working tree) |
|---|---|---|
| C04-01 | **IMPLEMENTED IN SOURCE — UNVERIFIED** | AMR: `amr_hydro_orchestrator.cpp:3186–3252` forwards `solver.numericsPolicy()`, `options.reconstruction_scheme`, `reconstruction_policy`, `options.riemann_solver_kind`; fixed grid: `hydro_amr_runtime.cpp:1319–1346`. Private engine `PrivateHydroEngine` (`hydro_shadow_reference.cpp:472`) forces CTU, keeps dual energy/floors/limiter/gamma, `makeHydroRiemannSolver(kind)`, per-halfstep `hydroTransportCourantFactors(dt_sub, a_sub)`; exact FLRW epochs `hydroShadowScaleFactorAtPhysicalFraction` (`:686`). Class `kTimeOnlyCtuPpm` only for production CTU/PPM, otherwise `kMixedSpatialTemporal` (`hydroClassFor`). |
| C04-02 | **IMPLEMENTED IN SOURCE — UNVERIFIED** (declared scopes §5) | Source owner `source_runtime.cpp:2247–2276` (global halves, replay lent at post-half end) and `:2009–2080` (level windows) → observer `hydro_amr_runtime.cpp:1625/1667` → `recordShadowChronologyPoint` (`:2146`) → `HydroShadowReference::finalizeChronologies` (`hydro_shadow_reference.cpp:1804`) → `runChronology` (`:1841`) using `SourceRuntimeImpl::ShadowSmoothReplay` (`source_runtime.cpp:1758`, production `BaryonicSourceStage` on a private `SimulationState`). Private composition: `ShadowCompositionSink` (`hydro_shadow_reference.cpp:184`, apply `:212`, now production clamp + reconciliation); MPI donors from the production table (`hydro_amr_runtime.cpp:4641` → `provideRemoteDonorFractions` `hydro_shadow_reference.cpp:1391`). Rejection: continuity-checked publication `publishOrDiscardCompletedShadowEpoch` (`hydro_amr_runtime.cpp:1752`) with carried loss count. Refusals: `skipped_unsupported_source`, `chronology_skipped_{no_replay_operator,replay_unavailable,replay_incomplete,missing_point,region_changed,recursive_level_sources,unsplit_compatibility}`. |
| C04-03 | **IMPLEMENTED IN SOURCE — UNVERIFIED** | Capture in the production synchronized sweep (`amr_hydro_orchestrator.cpp:3252`, provenance `hydroShadowGhostProvenance` `:4719`, step-start ghost shell loaded just above) → `observeProduction` (`:3272`, pre-reflux) → synchronized endpoint `completeSynchronizedShadowEndpoints` (`hydro_amr_runtime.cpp:2195`, called `:4371` after local+remote reflux; MPI via `exchangeDirectedHydroGhosts`) / recursive endpoint at each level clock (`amr_hydro_orchestrator.cpp:3844–3850`, MPI via `shadowEndpointSources` `hydro_amr_runtime.cpp:3574`) → `completeHydroShadowEndpointGhosts` (`amr_hydro_orchestrator.cpp:4760`, production `fillAmrHydroGhostCells` `:4915` on private copies, edge/diagonal sources via `edgeStencilSourceTouches`) → `provideEndpointGhosts` (`hydro_shadow_reference.cpp:1349`) → midpoint = linear interpolation in `PrivateHydroEngine::setGhosts` (`:501`), self slots refreshed from private state → two-halfstep CTU/PPM in `runHydroComparison` (`:1479`) over the production active set; stale slots cone-excluded (`hydroShadowDependencyConeMask` `:755`). |
| C04-04 | **IMPLEMENTED IN SOURCE — UNVERIFIED** | `hydroShadowRegionKeyForPatch` (`amr_hydro_orchestrator.cpp:4702`) and fixed-grid tiles `buildFixedGridShadowRegions` (`hydro_amr_runtime.cpp:1819`) → `hydroShadowPhysicalRegionKey` (`hydro_shadow_reference.cpp:632`); `prepareGlobalHydroShadowSelection` (`hydro_amr_runtime.cpp:1867`): one owner per AMR patch, local top-K, bounded `allgatherBytesBounded`, global `rankHydroShadowCandidates` (`:673`), two refusal votes; empty ranks enter every collective. Multi-rank uniform grids = level-0 AMR (production refuses multi-rank fixed-Cartesian, `hydro_amr_runtime.cpp:1120`). Candidate identity is decomposition independent; bitwise-equal floating aggregates are not claimed (per-region records sorted by key before summation). |
| C04-05 | **IMPLEMENTED IN SOURCE — UNVERIFIED** | Checked arithmetic + `MemoryClass::kDiagnostic` leases before allocation in `captureImpl` (`hydro_shadow_reference.cpp:959`), `openChronology` (`:1641`), endpoint fill, endpoint remote interpolation, selection, publication; node aggregate through `NodeMemoryPool` in every `reserve()` (`memory_governor.cpp:457`); new `noexcept` containment (`capture` `:943`, `observeProduction` `:1289`, `setGlobalCandidates` `:898`, orchestrator guarded preparation). Collectives only at rank-uniform points (selection, endpoint votes, publication); the donor-table hand-off adds none. Disabled shadow: no capture spec built (pointer null), no allocation, no new collective (`hydro_shadow_enabled` gates every call). |
| C04-06 | **IMPLEMENTED IN SOURCE — UNVERIFIED (fixtures written, NOT RUN)** | §6 below; ledger rows updated; publication payload labels class, provenance, time authority, chronology outcomes, composition rule and rejection loss. |

## 4. Physical operators: supported vs refused

| Operator / input | Hydro-interval classes | Coupled Strang class | Refusal (specific counter) |
|---|---|---|---|
| Hydro CTU/PPM + HLLC/HLLE, floors, dual energy | exact production settings, two halfsteps | same | — |
| Other production reconstruction/integrator | one private CTU/PPM step, `mixed_spatial_temporal_scheme` | same engine | — |
| Cosmological expansion a(t), H, z, dt/(a dx) | exact FLRW epochs | + source-owner fractions | — |
| Hydro-owned cooling, effective-ISM relaxation | private clones | private clones | non-clonable source → `skipped_unsupported_source` |
| Thermochemistry, ISM regime, chemistry fidelity, dust, SGS dissipation | outside window | private production `BaryonicSourceStage` | `chronology_skipped_replay_{unavailable,incomplete}`, `…no_replay_operator` |
| Composition advection | private donor-cell sink + production clamp/reconciliation; MPI donors from production table (synchronized) | private (global) / deferred like production (level windows) | missing row map → `skipped_external_passives`; unknown donor → `composition_remote_donor_faces` (cells cone-excluded) |
| Element diffusion, SGS shear (neighbour-coupled) | outside window | captured production increments, linear in time | — |
| Discrete events (SF, SSP returns, winds, SN, BH) | outside window | captured production increment at interval end; no RNG/queue/deposition | — |
| Gravity KDK half-kicks | outside window | outside window | — |
| Unsplit compatibility dispatch / recursive global halves | — | — | `chronology_skipped_unsplit_compatibility`, `…recursive_level_sources` (level windows used instead) |

## 5. AMR boundary / ownership cases

| Case | Time authority in halfstep 2 | Status |
|---|---|---|
| In-patch periodic wrap, open, reflective | refreshed from private state | supported |
| Same-level local / MPI | production filler endpoint, linear in time | supported (synchronized + recursive) |
| Edge/corner (diagonal) local / MPI | production filler endpoint | supported |
| Coarse→fine (spatial or temporal history) | production filler endpoint at t_end (open coarse history in recursive) | supported |
| Fine→coarse, synchronized | production filler endpoint (post-reflux) | supported |
| Fine→coarse, recursive level clock | no legal endpoint (finer level not advanced) | refused per slot: step-start only, radius-4 cone excluded, `ghost_slots_step_start_only` / `bounded_projection` |
| Closure-mapped deep layers | production filler endpoint | supported |
| Endpoint fill failed / geometry mismatch | step-start only, cone-excluded | counted `skipped_endpoint_*` |
| Partial activity, suppressed duplicate periodic faces | production active set unchanged | supported |
| Composition donors across MPI, synchronized | production step-start donor table | supported (this campaign) |
| Composition donors across MPI, recursive | receiver fractions, composition cone excluded | counted fallback (declared scope) |

Declared approximation scopes that remain (not mandatory gaps; all labelled
in the payload): captured (not re-integrated) diffusion/SGS-shear/event
forcing; imported ghosts in the coupled chronology are production endpoint
states; recursive fine→coarse cone exclusion; recursive-MPI composition
donor fallback; the final step of a run is never published (no later
acceptance barrier exists — `HydroAmrRuntime` has no run-end hook).

## 6. Conservation, sidecars, restart, MPI, memory — why production is unaffected

* The shadow only reads `SimulationState`, the production composition
  sidecar and the production donor table (spans, `const`); every private
  state, ghost shell, composition copy, snapshot and scratch is owned by
  `HydroShadowReference` under a diagnostic lease.
* No flux register, reflux, composition accumulator, ledger, RNG cursor,
  event queue, controller, CFL/acceptance decision or restart field is
  written. `provideRemoteDonorFractions` is called after the production
  `setRemoteDonorFractions` and reads the same immutable assembled table.
* Restart/snapshot schema unchanged; shadow state is never persisted.
* Disabled: every new call is behind `hydro_shadow_enabled`; the orchestrator
  flag defaults to false; no allocation or collective is added.
* MPI: the hold changes only *when* a local comparison runs (after the
  already-collective composition exchange, before the rank-uniform
  `kPostHalfBegin` record). No new collective, no per-patch collective.

## 7. Test-source fixtures (ALL NOT RUN)

| Test (CTest name) | Covers |
|---|---|
| `unit_hydro_shadow_c04_final` (**new**) | noexcept containment; remote-donor hold/resolution vs local-row resolution (identical private composition); counted fallback without/with malformed/unfulfilled table; reconciliation bounded; held work counted then dropped on rejection |
| `unit_hydro_shadow_reference`, `unit_hydro_shadow_c04_coupled`, `unit_hydro_shadow_c04_identity`, `unit_hydro_shadow_c04_chronology` | revision-6 fixtures (policy, 8³ periodic/reflective, cones, physical keys + global K under 1/2/4/uneven/empty ranks, production `BaryonicSourceStage` replay, fail-closed reasons) |
| `integration_amr_hydro_shadow_c04_production` | production 8³ coarse/fine + diagonal + same-level endpoints, partial activity + suppressed periodic faces, bit identity on/off/refused, governor refusal, rejected attempt, recursive level clock |
| `integration_hydro_shadow_c04_workflow_strang` | real `ReferenceWorkflowRunner` T01 Strang cosmological baryonic deck, digest identity |
| `integration_hydro_shadow_c04_mpi_two_rank`, `…_four_rank` | MPI endpoint exchange, decomposition-independent keys (level-0 uniform grid with world-size-dependent IDs/owners), uniform and asymmetric refusal, digest identity |

## 8. Owner compilation handoff (instructions only — nothing was executed)

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R "hydro_shadow|amr_hydro_shadow"
ctest --preset test-cpu-debug --output-on-failure
```

Then the MPI preset from `docs/build_instructions.md` for
`integration_hydro_shadow_c04_mpi_two_rank` / `_four_rank`. A distributed
composition deck (≥2 ranks, `hydro_shadow_enabled = true`, synchronized AMR,
composition on) exercises the new donor hand-off end to end; expect
`composition_remote_donor_slots_resolved > 0` and
`composition_remote_donor_faces == 0` near rank interfaces.

Likely build-integration risk points found by source reading:

1. `hydro_shadow_reference.cpp` `ShadowCompositionSink::apply` uses
   `std::min({rows, real, private_state.size()})` (all `std::size_t`) and
   `GasCompositionSidecar::scaleRow/setValue` — check they are callable on
   a dense (`emptyLike(n, false)`) sidecar.
2. `capture()` returns `captureImpl(spec)` (move-only `Sample`, guaranteed
   elision) inside a `try` in a `noexcept` function.
3. The new test's `static_assert(noexcept(shadow.setGlobalCandidates({})))`
   relies on `{}` initializing the by-value vector parameter.
4. Indentation of the guarded block in `amr_hydro_orchestrator.cpp`
   (`try {` around the existing selection block) is intentionally unchanged
   to keep the diff minimal.
5. Revision-6 risk points listed in `chui_c04_rev5_cloud_source_closure.md`
   still apply (never compiled).

## 9. Remaining blockers

* **Mandatory source blockers:** none found after this campaign.
* **External physics data:** production rate/yield/cooling tables remain
  as recorded in `chui_post_patch58_final_source_closure.md` §4 (not a C04
  item).
* **Optional research:** private re-integration of neighbour-coupled
  diffusion/SGS shear; recursive-MPI private composition donors (would need
  level-clock donor exchanges production does not perform); a run-end
  acceptance hook to publish the final epoch.
* **Post-build qualification (owner):** compilation, all fixtures above,
  MPI 1/2/4 ranks incl. empty/uneven, recursive AMR, restart equivalence,
  memory-pressure/node-pool refusal, convergence of the disagreement metric,
  performance overhead with shadow disabled/enabled.

## 10. Addendum — design-audit repairs (2026-11-17, same day)

A subsequent audit against the original design baseline (v01/v10/v02)
found and this addendum repairs (source only, UNCOMPILED / UNTESTED):

| Item | Repair | Fixture (NOT RUN) |
|---|---|---|
| C04-05 P1: shadow memory could starve mandatory physics (v10 §5.1/§5.2) | every shadow lease uses `MemoryGovernor::reserveBestEffort` (green band only); `HydroShadowReference::reclaimBestEffort()` is registered as the governor's best-effort reclaimer, so a refused mandatory `reserve()` drops pending samples/chronologies (`evicted_for_mandatory_memory`) and retries once; busy guard prevents reclaim while the shadow iterates its lists | `unit_memory_governor_best_effort` |
| C04 cost telemetry (v02 §10.3) | `private_work_seconds_{sum,max_rank}` published | — |
| Specific skip reason | cloned source needing a retry → `skipped_source_interval_incomplete` | — |
| Shadow regime gate mirror | the private solver receives the full solver-row I01 gate; chronology halfsteps follow the private (replayed) regime tags | `unit_ism_regime_gate` |
| C03/C04 chemistry shadow | key-only selection, cap truncation and reference failures counted, controller bias held on incomplete evidence, no per-row heap allocation | `unit_chemistry_shadow_sampling` (skip 77 without detailed backend) |

The C04 verdict is unchanged: **SOURCE IMPLEMENTED — UNVERIFIED.**
