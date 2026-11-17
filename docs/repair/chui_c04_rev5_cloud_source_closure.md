# C04 sampled shadow hydrodynamics — revision-6 source closure (cloud session, 2026-11-17)

> **Superseding note (2026-11-17, final C04 source campaign):** see
> `chui_c04_final_precompile_handoff.md`. It fixes four residual source
> defects of this revision (exception containment, production-rule private
> composition clamp/reconciliation, synchronized-MPI remote composition
> donors from the production donor table — superseding §5 item 5 for
> synchronized stepping — and rejected-attempt accounting). The rest of this
> note remains valid.

**Status: C04 SOURCE IMPLEMENTED — UNVERIFIED.** Every mandatory C04 source
contract listed by the revision-5 audit now has an integrated source path, with
the approximation scopes declared below. **Nothing in this campaign was
compiled, run, tested, benchmarked or scientifically validated** (owner-reserved
stages). Scientific readiness is a separate, later judgment; disagreement
metrics are sampled numerical disagreement, never physical error bars.

* Baseline: branch `claude/nifty-bell-6y17me`, HEAD
  `2ce73fbfe970147cd229b1f6a43d7fe953078166` ("shadow baryons patching"),
  identical in C04 scope to the revision-5 archive audited on 2026-11-17.
  No newer C04 commits existed; nothing was reset or reverted.
* Scope respected: no MHD, radiative transfer, AGN or SHAI-HULUD work; no
  `docs/research/shai_hulud/**` change; no new config key; no snapshot,
  restart, provenance or PM/hydro authoritative format change.
* Mode: feature implementation + localized integration repair (AGENTS.md §2).

## 1. Implementation dependency map

### 1.1 Sampled physical states (initial / terminal)

| Comparison class | Initial private state | Terminal production state compared |
| --- | --- | --- |
| `time_only_ctu_ppm_subdivision` / `mixed_spatial_temporal_scheme` (hydro interval) | production hydro input of the sampled patch (post-S1), real cells **and** step-start ghost shell, auxiliary e | production solver output of the same patch, **pre-reflux** (both trajectories are reflux-free) |
| `bounded_projection` | as above | as above, restricted to cells outside the dependency cones of stale ghost slots |
| `coupled_full_physics_strang_interval` (global Strang) | `kPreHalfBegin` snapshot (U0, before S1) | `kPostHalfEnd` snapshot (after S2 smooth operators and discrete events) |
| `coupled_full_physics_strang_interval` (recursive level window) | level-L `kPreHalfBegin` of the first level-L window | level-L post-half end of the same window (no events inside a level window) |
| hydro-only composition check | private composition rows at hydro input | committed production composition after the hydro stage (`kPostHalfBegin`, before regrid), only when the sample interval ends at that point |

### 1.2 Source-order timeline of one accepted interval (synchronized stepping)

```
KickPre | Drift | ForceRefresh |
  S1 pre half  : [kPreHalfBegin] transport [kPreHalfAfterTransport] thermochemistry [kPreHalfEnd]
                 -> attempt opened + global selection + chronology records (observer, rank-uniform)
  Hydro        : capture(t0 ghosts) -> production advance -> observeProduction (pre-reflux)
                 -> reflux (local + remote) -> post-advance directed exchange (MPI, rank-uniform)
                 -> completeHydroShadowEndpointGhosts (production filler, read-only) -> hydro comparison
                 -> composition commit -> [kPostHalfBegin] (+ hydro-only composition check) -> regrid
  S2 post half : [..] thermochemistry [kPostHalfAfterThermochemistry] transport [kPostHalfAfterTransport]
                 discrete events [kPostHalfEnd] -> private chronology replay (lent replay) -> epoch closed
| KickPost | Analysis
next attempt start: continuity check (t_begin == previous t_end, new step index) -> publish, else discard
```

Private refined chronology replayed on the sampled region:
`[S(q0,q1) H(q0,q2) S(q1,q2)] [S(q2,q3) H(q2,q4) S(q3,q4)] + events`, with
`q_k = t0 + k dt/4`, S the production smooth operator recomputed privately,
H the private CTU/PPM halfstep with the production ghost authority.

Recursive Berger–Oliger stepping: each level clock dispatches
`S_L(pre) -> [kHalfBegin, kAfterTransport, kHalfEnd]`, `H_L` sweep, endpoint
ghosts completed by the driver **before** `S_L(post)`, then
`S_L(post) -> [kHalfBegin, kAfterThermochemistry, kAfterTransport, kHalfEnd]`
with a replay over the **level** interval. Composition face transfers are
committed by production only at the coarse end; the private level window
defers them identically.

### 1.3 Ghost-history owner

| Ghost provenance | Time authority in the private halfstep 2 | Owner of the endpoint value |
| --- | --- | --- |
| self physical (open/reflective), self periodic wrap | refreshed from the PRIVATE state (`fillHydroBoundaryGhostCells`) | private reference |
| same-level local, same-level MPI, edge/diagonal, coarse→fine (spatial or temporal history), fine→coarse, closure-mapped deep layers | linear time interpolation between the production step-start value and the production endpoint value | production AMR ghost filler (`fillAmrHydroGhostCells`) re-run read-only at the interval end by `amr::completeHydroShadowEndpointGhosts` |
| recursive fine→coarse (finer level not advanced yet), failed/rejected endpoint fills | step-start only; dependency cone (Chebyshev radius 4) excluded and counted | none (no legal endpoint at that point) |

Halfstep 1 always uses the exact production step-start ghost shell.

### 1.4 MPI collective locations (all rank-uniform)

| Location | Collective | Decision uniformity |
| --- | --- | --- |
| attempt start (source pre half or hydro stage) | publish previous epoch: allreduce sums/max, allgather of region records | pending-publication flag and continuity times are identical on all ranks |
| selection | allreduce(production AMR coverage), refusal vote, bounded allgather of ≤K (key, level) per rank, second vote | config + step cadence |
| synchronized endpoint | allreduceMax(pending), `exchangeDirectedHydroGhosts` (own collective votes) | config + cadence |
| recursive endpoint | `hooks->globalMax(pending)`, `hooks->shadowEndpointSources` (directed exchange + failure vote) | every rank traverses every level clock |
| capture, compare, endpoint fill, chronology replay, level hooks | **none** (local only) | — |

## 2. Batch ledger

| Batch | Contract | Status | Main source |
| --- | --- | --- | --- |
| A | private ghost authority + true AMR active sets (C04-03) | implemented, unverified | `src/hydro/hydro_shadow_reference.cpp` (`capture`, `PrivateHydroEngine::setGhosts`, `hydroShadowDependencyConeMask`, `provideEndpointGhosts`, `finalizePendingEndpoints`), `src/amr/amr_hydro_orchestrator.cpp` (`hydroShadowGhostProvenance`, `completeHydroShadowEndpointGhosts`, synchronized sweep capture/observe, recursive driver hook) |
| B | private source chronology, composition sink, Strang/levels (C04-02) | implemented, unverified | `hydro_shadow_reference.cpp` (`ShadowCompositionSink`, `ShadowSourceProperties`, `runChronology`, `openChronology`, `recordChronologyPoint`, `finalizeChronologies`, `compareCompositionAfterCommit`), `src/workflows/source_runtime.cpp` (`ShadowSmoothReplay`, `notifyShadow`, `notifyLevelShadow`), `include/cosmosim/workflows/hydro_shadow_chronology.hpp` |
| C | physical identities + global K (C04-04) | implemented, unverified | `hydroShadowPhysicalRegionKey`, `amr::hydroShadowRegionKeyForPatch`, `HydroAmrRuntimeImpl::prepareGlobalHydroShadowSelection` / `buildFixedGridShadowRegions`, region-result wire + sorted aggregation |
| D | resource governance, fixtures, diagnostics, ledger (C04-05/06) | implemented (fixtures written, NOT run) | leases in `capture`, `openChronology`, `completeHydroShadowEndpointGhosts`, `shadowEndpointSources`, `publishHydroShadow`; tests listed in §8 |

## 3. C04-01 … C04-06 source-evidence matrix

| ID | Verdict | Source evidence |
| --- | --- | --- |
| C04-01 numerical policy equivalence | **SOURCE IMPLEMENTED — UNVERIFIED** (preserved, hardened) | `PrivateHydroEngine` builds `HydroCoreSolver` from the production `HydroNumericsPolicy` (CTU forced, dual energy preserved), production `HydroReconstructionPolicy` (limiter, floors, gamma), `makeHydroRiemannSolver(production kind)`; per-substep `hydroTransportCourantFactors(dt_sub, a_sub)`; halfstep epochs now from exact FLRW inversion (`hydroShadowScaleFactorAtPhysicalFraction`, same 64-step bisection as the source owner) between production endpoints. |
| C04-02 coupled cosmology + physical source chronology | **SOURCE IMPLEMENTED — UNVERIFIED** with declared approximations (§5) | `SourceRuntimeImpl::ShadowSmoothReplay::advanceSmooth` runs the production `BaryonicSourceStage` (ISM regime, chemistry fidelity, thermochemistry, dust, SGS) on a private single-region `SimulationState`; `HydroShadowReference::runChronology` (refined Strang, replay-identity probe, captured diffusion/SGS-shear/event forcing, private composition sink + metal mirror); level windows via `onLevelSourceShadowBoundary` + `finalizeChronologies(replay, level)`; composition-aware private temperature (`ShadowSourceProperties`, `temperatureFactorForMu`); accepted-epoch publication (`publishOrDiscardCompletedShadowEpoch`). Composition-enabled runs are no longer refused (`ShadowCompositionSink`). |
| C04-03 production 8³ AMR + MPI/coarse-fine temporal ghosts | **SOURCE IMPLEMENTED — UNVERIFIED** | No interior guard; `hydroShadowGhostProvenance` (all classes); endpoint ghosts from the production filler (`completeHydroShadowEndpointGhosts`) for synchronized (post-reflux, MPI via `exchangeDirectedHydroGhosts`) and recursive (level clock, MPI via `shadowEndpointSources`) stepping; midpoint = linear interpolation; self slots refreshed privately; stale slots cone-excluded; production active set (partial activity, suppressed duplicate lower periodic faces) passed unchanged. |
| C04-04 decomposition-independent selection | **SOURCE IMPLEMENTED — UNVERIFIED** | `HydroShadowCandidate::region_key` from (level, nearbyint(origin/width), dims); one owner per AMR patch; fixed-grid physical tiles of `k_hydro_shadow_fixed_grid_region_extent` cells on the global lattice (multi-rank fixed grid is validation-only: no candidate, counted `skipped_region_split_across_ranks`); global K by score(key, step, seed); identity decided before admission; per-region records gathered and reduced in physical-key order (`selected_region_keys`, `completed_unique_physical_regions`). |
| C04-05 bounded memory, no authoritative side effects, MPI-safe failure | **SOURCE IMPLEMENTED — UNVERIFIED** | Checked arithmetic + `MemoryClass::kDiagnostic` leases before allocation for samples (state ×3, endpoint ghosts, production result, geometry/active copies, composition rows/fractions/deltas, two scratch envelopes), chronology snapshots, endpoint fill working set, endpoint remote interpolation, selection and publication buffers; every failure local and counted; endpoint/exchange refusals caught uniformly; the canonical state is only read (fixture asserts bit-identity). |
| C04-06 fixtures, diagnostics, ledger | **SOURCE IMPLEMENTED — UNVERIFIED (fixtures written, NOT run)** | §8; publication payload documents class, provenance, time authority, chronology outcomes and approximation scope. |

## 4. Physical-operator coverage matrix

| Operator | Hydro-interval classes | Coupled Strang class |
| --- | --- | --- |
| Cosmological expansion (a, H, z, dt/(a dx)) | exact FLRW epochs per private halfstep | same + replay intervals from the source owner's inversion |
| Hydro-owned smooth sources (Part-1 cooling, effective-ISM relaxation) | private clones (`cloneForShadow`); refused if not clonable (counted) | same |
| Thermochemistry / ISM regime / chemistry fidelity / dust / SGS dissipation | not in window (external) | **recomputed privately** with the production operator |
| Composition advection | private donor-cell sink (production rule) | private sink (global) / deferred like production (level windows) |
| Element diffusion, SGS shear (neighbour-coupled) | not in window | captured production increments, apportioned linearly in time |
| Discrete events (star birth, SSP returns, winds, SN, BH) | not in window | captured production increment at the interval end (immutable forcing; no RNG, no queue, no deposition) |
| Gravity (KDK half-kicks) | outside window | outside window (kicks bracket S1·H·S2) |
| Interface synchronization halves (recursive post-children diffusion) | — | outside the level window |

## 5. Declared approximation scopes and remaining limits

1. Captured forcing (diffusion/SGS shear, discrete events) is replayed, not
   re-integrated: disagreement in those operators is not measured.
2. Imported ghosts in the coupled chronology are production hydro endpoint
   states (post-S1 at t0, post-H at t1), not private-chronology neighbours.
3. Endpoint ghosts of synchronized runs are post-reflux accepted values; the
   compared production patch is pre-reflux (both trajectories reflux-free).
4. Recursive fine→coarse slots have no endpoint at the level clock: cone
   exclusion (`bounded_projection`), counted.
5. Remote ghost donor composition is not exchanged: inflow from such donors
   uses the receiver's fractions and those cells are excluded from the
   composition metric (`composition_remote_donor_faces`).
6. The fixed grid is single-rank in production; multi-rank fixed-grid
   validation runs submit no candidate (counted).
7. Published floating aggregates are reduced in physical-key order; each
   per-region value is decomposition independent only to the extent the
   production patch state is (bitwise identity across decompositions is not
   claimed).
8. The epoch of the final step of a run is never published (no later
   acceptance barrier). Shadow state is never persisted (no restart schema).

## 6. Memory / RSS and MPI ownership

* Owner: `HydroShadowReference` (rank-local, diagnostic class). Lifetimes:
  sample leases until compared (or retained to the chronology end, at most
  one epoch); chronology leases until finalized; endpoint fill and
  publication leases scoped to the call. `beginEpoch` drops everything.
* Bounds: global K (≤1024 by config validation) × (patch storage + ghost
  shell + chronology snapshots); publication gather ≤ 3 records × K × 512 B.
* `observed_capacity_high_water_bytes` reports measured vector capacities;
  opaque allocator overhead remains outside the measurement.

## 7. Public interface changes and migration notes

* `include/cosmosim/hydro/hydro_shadow_reference.hpp`: new `HydroShadowPatchCapture`
  capture API, `observeProduction`, endpoint protocol, chronology API,
  provenance/authority/class enums, region keys, region-result wire. The
  legacy `capture(...)`/`compare(...)` overloads remain (one region, default
  provenance, immediate comparison). `HydroShadowCandidate` gained
  `region_key` (appended; `{id, level}` initializers keep the legacy key).
* `include/cosmosim/amr/amr_hydro_orchestrator.hpp`: `ProductionAmrHydroOptions`
  gained four appended C04 fields (defaults preserve behaviour);
  `RecursiveAmrDistributedHooks::shadowEndpointSources` (non-pure, default
  empty); `hydroShadowRegionKeyForPatch`, `hydroShadowGhostProvenance`,
  `completeHydroShadowEndpointGhosts`.
* `include/cosmosim/workflows/source_runtime.hpp`:
  `setHydroShadowChronologyObserver` (non-pure, default no-op).
* New `include/cosmosim/workflows/hydro_shadow_chronology.hpp`.
* Config: no new key. `hydro_shadow_max_patches_per_rank` is the GLOBAL K of
  physical regions for AMR **and** fixed grid.
* Reproducibility impact: none on evolved state (shadow is read-only;
  fixtures assert digest/bit identity); publication moved to the acceptance
  barrier of the next attempt.

## 8. Fixtures written (NOT compiled, NOT run)

| Test | Covers |
| --- | --- |
| `tests/unit/test_hydro_shadow_reference.cpp` (updated) | policy equivalence, budgets, no 15-cell guard, malformed active set |
| `tests/unit/test_hydro_shadow_c04_coupled.cpp` (updated) | 8³ periodic/reflective, sources, cosmology, dual energy, imported-ghost cone exclusion vs refusal |
| `tests/unit/test_hydro_shadow_c04_identity.cpp` (new) | physical keys, 1/2/4-rank + uneven + empty-rank selection, no ID-zero duplicates, global K, deterministic aggregation, 8³ and 2D cones |
| `tests/unit/test_hydro_shadow_c04_chronology.cpp` (new) | production `BaryonicSourceStage` replay, replay identity = 0, coupled chronology, partial activity, every fail-closed reason, composition-map requirement, no sidecar aliasing |
| `tests/integration/test_amr_hydro_shadow_c04_production.cpp` (new) | production 8³ hierarchy with coarse/fine + diagonal + same-level endpoints, partial activity + suppressed periodic faces, shadow on/off/refused bit identity, governor refusal, rejected attempt, recursive level clock with fine→coarse cones |
| `tests/integration/test_hydro_shadow_c04_workflow_strang.cpp` (new) | real `ReferenceWorkflowRunner` T01 Strang cosmological baryonic deck (thermochemistry, composition, dust, ISM, SF, SN), digest identity, published chronology outcomes |
| `tests/integration/test_hydro_shadow_c04_mpi.cpp` (new; 2 and 4 ranks) | MPI endpoint exchange, decomposition-independent selected keys, uniform and rank-asymmetric refusal, digest identity |

Owner verification order: `cmake --preset cpu-only-debug`,
`cmake --build --preset build-cpu-debug`,
`ctest --preset test-cpu-debug --output-on-failure` (then the MPI preset for
the two `integration_hydro_shadow_c04_mpi_*` tests).
