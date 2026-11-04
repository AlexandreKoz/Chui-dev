# CHUI Part 3 — distributed periodic AMR, remote sources, rollback, science output (handoff)

Date: 2026-11-03/04. **Evidence status: SOURCE ONLY. Nothing was configured,
compiled, linked, tested, benchmarked, run under MPI or simulated (explicit
owner instruction). Every test and command below is NOT RUN — USER-RESERVED.**
The Part 1/2/2.5 parent has never been compiled either; compilation errors in
this large change set are expected and must be fixed at the first human build.
Nothing here is numerically verified or scientifically calibrated.

Companion ledger: `docs/repair/chui_part3_scope_and_closure_20261103.md`.

## 1. Git record

| Item | Value |
|---|---|
| Parent / base | `f5a92555de2c199d860755a727a7f661612912aa` (Part 2 + 2.5) |
| Branch | `claude/bold-dirac-ig5u6l` (session-designated; created at the parent SHA, did not exist remotely before; `main`, Part 1 and Part 2 branches untouched) |
| Checkpoint A | `0e5f9f8` periodic cross-patch/cross-rank topology, deep ghost layers |
| Checkpoint B | `9820d93` automatic regrid, composition regrid, distributed composition |
| Checkpoint C | `aef9b9e` element diffusion, cross-rank feedback + SF gradients, wind channel, associations |
| Checkpoint G | `d34f06d` composition / stellar-population science snapshots |
| Checkpoint D | `20a5d88` collective global-step rollback (rung zero) |
| Checkpoint E | `3865b08` TreePM–AMR force fixture, zoom contamination |
| Checkpoints F/H/I | `4280cdf` global adaptive chemistry control, predictive refinement, four-mode decks, benchmark |
| Final | the commit adding this file (see `git log`) |

Ordinary `git push` only; no force push, rebase, merge or PR. Pushes reported
success by the remote.

## 2. Closed / open matrix (feature register and carried debt)

Status vocabulary per the directive. "IMPLEMENTED IN SOURCE" = connected
production path by inspection only; always NOT RUN.

| ID / debt | Source (entry point) | Modes | Unsupported / fail-closed | Tests (NOT RUN) | Status |
|---|---|---|---|---|---|
| A01/A02 periodic cross-patch + cross-rank hydro | `amr/internal/amr_periodic_topology.hpp`; `amr_hydro_orchestrator.cpp` (`boundaryClassForSide`, `validatePeriodicAmrPatchTopology`, wrapped probes/registers); `amr_ghost_fill.cpp`; `parallel/distributed_memory_amr.cpp` (`AmrInterfaceTopologyOptions`); `hydro_amr_runtime.cpp` | cosmo_cube, zoom_in (periodic) | an untiled periodic box is refused | `unit_amr_periodic_cross_patch` | IMPLEMENTED IN SOURCE |
| H01 multi-layer ghosts (PLM/PPM full order) | `amr_hydro_geometry.cpp` deep stencil-only layers; `ProductionAmrHydroOptions::ghost_layers` | all | reflective patches narrower than the mirror depth, deep layers beyond the domain edge: counted zero-gradient reduction; CTU ghost-side transverse terms still normal-only (counted) | same | IMPLEMENTED IN SOURCE |
| F6/A01 automatic refine/derefine | `amr_production_regrid.{hpp,cpp}`; `HydroAmrRuntime::runAutomaticAmrRegrid`; `numerics.amr_regrid_*` | all | octets split across ranks skipped; grandparent link may be unrecoverable (blocks further coarsening); first-order (piecewise-constant) prolongation; ops individually transactional | `unit_amr_production_regrid`, `unit_config_chui_part3` | IMPLEMENTED IN SOURCE |
| D01 composition regrid + distributed composition | `amr_hydro_orchestrator.cpp` (`CompositionRegridPlan`); `composition_transport.cpp`; `parallel/composition_exchange.{hpp,cpp}`; `HydroAmrRuntime::completeDistributedCompositionTransport` | all | lane count > 32 refused | `unit_amr_production_regrid` | IMPLEMENTED IN SOURCE |
| D01 element-resolved diffusion | `physics/element_diffusion.{hpp,cpp}`; `SourceRuntime::executeDiffusionAndSgsShear`; `DistributedElementDiffusionCoupling` | all gas modes | lane count > 32 refused by the fixed distributed record contract | `unit_chui_part3_remote_sources` | IMPLEMENTED IN SOURCE, NOT RUN |
| S01 remote SF gradients | `star_formation_gradients.hpp` remote cells; `SourceRuntime::importRemoteStarFormationNeighbors` | all distributed | — | (MPI test not written: open) | IMPLEMENTED IN SOURCE |
| F01 cross-rank feedback | `supernova_feedback.{hpp,cpp}` (remote candidates, step-start reference, per-cell kinetic closure); `SourceRuntime::importRemoteFeedbackCandidates` | all | query volume O(active sources × ranks) | `unit_chui_part3_remote_sources` (in-process equivalence) | IMPLEMENTED IN SOURCE |
| F02 winds / pre-SN | `physics/stellar_wind.{hpp,cpp}`; `physics.ssp_wind_*` | galaxy/cluster/zoom | **BLOCKED BY DATA**: no wind table bundled; refused with `fb_presn_wind_momentum_lsun_per_c` | same | IMPLEMENTED IN SOURCE, data-blocked |
| S02 associations (RESEARCH) | `mergeAssociationSources`; `physics.ssp_association_*` | all | grouping per owner rank (decomposition-dependent); internal velocity-dispersion energy of members not ledgered | same | EXPERIMENTAL/UNQUALIFIED |
| T01/F02 global step rollback | `TimeCoordinator::runRungZeroSegment`; `numerics.global_step_max_retries` | rung zero | hierarchical runs refused (they use the multirate PM transaction); full-state copy per step (governed); adaptive-chemistry bias of a rejected attempt not reverted | `unit_config_chui_part3` (keys only; rollback integration test open) | IMPLEMENTED IN SOURCE |
| T02 recursive multirate gas/particle/source | — | — | **open**: existing 2-level local subcycling unchanged; no ≥3-level recursive gas/source multirate | — | open |
| G01 fine-level gravity | `test_treepm_amr_force_resolution.cpp`; `ProductionAmrRegridPolicy::gravity_softening_code` (counts refinement below softening) | all | no correction implemented (needed only if the fixture fails) | `unit_treepm_amr_force_resolution` | REFERENCE fixture; correction open |
| zoom contamination | `analysis/zoom_contamination.{hpp,cpp}` | zoom_in | not yet emitted by the runtime report | `unit_zoom_contamination` | IMPLEMENTED IN SOURCE (API) |
| C03/C04 global adaptive chemistry + shadow | `ThermochemistryOperator::setGlobalPromotionBias`, `nextChemistryPromotionBias`; `SourceRuntime::updateGlobalChemistryControl` (sidecar `chui_adaptive_chemistry_control_v1`) | chemistry runs | no timescale/stiffness criterion yet; fixed control constants | `unit_chui_part3_adaptive_and_decks` | EXPERIMENTAL/UNQUALIFIED |
| A03 physics indicators | jump (contact/shock), comoving Jeans, cell mass | all | cooling-length, chemistry-interface, vorticity indicators open | regrid tests | partial |
| A04 predictive refinement | `numerics.amr_predictive_feedback_age_years` | all | host lookup O(stars × patches) | decks test | EXPERIMENTAL/UNQUALIFIED |
| V02 error/resource coordination | — | — | open beyond the C03 bias loop | — | open |
| V01 science snapshots | `snapshot_hdf5.cpp` composition + stellar population datasets | all | readers do not import the new datasets (science export only) | `integration_chui_part3_science_snapshot` (HDF5) | IMPLEMENTED IN SOURCE |
| M01/M02 memory, CPU | regrid/composition/diffusion/rollback admissions | all | no new RSS measurement; benchmarks only partially written | `bench_chui_part3_kernels` | partial |
| Four-mode decks | `configs/chui_part3/*_part3.param.txt` | four modes | decks use `uv_background_model = none` and no data tables | `unit_chui_part3_adaptive_and_decks` | IMPLEMENTED IN SOURCE |
| Remote auxiliary dual energy (Part 1 debt) | — | — | open (remote ghosts still fall back to E−K, counted) | — | open |
| Transactional AMR sweep (Part 1 debt) | covered by the global step rollback when enabled | rung zero | no AMR-local undo | — | partial |
| `massiveReturn` monotonicity (Part 2.5) | unchanged; 401-point check retained | — | static review found no structural fix in scope | existing | open (owner run decides) |
| X01 MHD/RT/new AGN | not implemented by design | — | — | — | BOUNDARY |

## 3. Production call graph (Part 3 additions marked +)

```
mode dispatch -> TimeCoordinator::runRungZeroSegment
  + capture GlobalStepRetryPoint (if global_step_max_retries > 0)
  -> dt selection (+ retry limit)  -> executeSingleStep
     -> gravity kicks (TreePM; unchanged)
     -> HydroAmrRuntime: + periodic topology / deep ghosts / periodic MPI planner
        -> advance[Distributed]ProductionAmrHydro -> reflux (local + remote)
        -> composition: single rank prepare/commit | + distributed donor/correction exchanges
        -> + runAutomaticAmrRegrid (planner -> id windows -> applier)
     -> SourceRuntime: baryonic stage (+ global chemistry bias) -> SF (+ remote neighbours)
        -> stellar population (+ wind channel, + associations, + remote candidate halo,
           + step-start reference, + remote delta delivery) -> refresh
        -> diffusion (+ element-resolved when composition is enabled)
  + collective rejection verdict -> restore + retry | rethrow
  -> output (science snapshot + composition / population datasets) / checkpoint
```

## 4. AMR topology and boundaries

* Periodic box `[0, L)` per axis. A patch spanning an axis wraps in-patch; any
  other face lying on a periodic box face is a same-level or coarse/fine
  interface with its wrapped neighbour (minimum-image face contact, wrapped
  probe points and fine-to-coarse ghost volumes, wrapped register keys).
  Every non-spanning periodic-axis face must have a neighbour.
* Ghost layers = reconstruction half-width (PCM 1, PLM 2, PPM 3). Layer 1
  carries the face; deeper layers are stencil-only (no face, no flux, never
  advanced), filled from the cell containing the probe at `d` widths, at any
  legal level relation; physical deep layers use zero-gradient (open) or
  mirror (reflective) sources.
* MPI: `AmrInterfaceTopologyOptions` (periodic extent, ghost layers) drives
  peer discovery (periodic envelope images), face masks (minimum image) and
  per-face source depth (`ceil(layers × width ratio)`, clamped to the patch).
* Regrid only at synchronized boundaries (planner refuses with temporal
  histories or pending registers). Same-level patch faces are evaluated from
  both sides with identical inputs: conservation to roundoff, not by a unique
  owner (limitation).

## 5. Time integration

* Rung zero (global KDK) is the only driver with global step rollback.
  Legal retry = same model, half interval; illegal (model downgrade,
  resource refusal, invalid state) fails on every rank. Collective verdict by
  two `allreduceMax`. Captured: SimulationState (all sidecars, hierarchy,
  composition, star records), integrator, both schedulers, pending output,
  gravity endpoint force (existing cycle-transaction API). Restored state gets
  bumped gravity-source and cell-index generations.
* Recursive multirate gas/particles/sources beyond the existing 2-level local
  AMR subcycling: **not implemented** (open).
* Event chronology: unchanged Part 2 cumulative-Poisson records restored by
  the rollback; wind budgets are cumulative-table differences (subdivision
  additive).

## 6. Gravity

`unit_treepm_amr_force_resolution` compares TreePM with a direct
minimum-image reference at levels 0–3. Acceptance: error never exceeds
max(2 × level-0 error, 5e-2) and the split composition error stays below
1e-12. If the owner run fails it, that is a genuine finding: a fine-force
correction would then be required (not implemented). Refinements below the
gas softening are counted as an unqualified regime and never blocked.

## 7. Adaptive research models

C03/C04: bias `b` widens promotion thresholds by `1 + b`. `b` rises by 0.5
(capped at 4) when the global shadow max relative species difference exceeds
0.1, and falls by 0.25 (floor 0) below 0.01. It is persisted for restart.
These are fixed research constants, not calibrated. A04 predictive refinement
force-refines host cells of stars younger than the configured age. Physical
budgets are unchanged; only the resolution changes.

## 8. Remote sources and conservation ledgers

* Element diffusion: one linear operator for every mass lane. Lane totals and
  linear constraints are exact up to roundoff. Positivity uses `dt_sub ≤ cfl M/ΣC`.
  Exceeding `max_subcycles` fails closed without mutation.
* Feedback: global k-NN stencil over local and remote cells (ties broken by stable id).
  Weights are normalized globally. Per-cell kinetic cross-term closure is
  ledgered (`kinetic_cross_term_to_thermal_code`, `_unclosed_code`). Remote
  deltas are delivered inside the stage transaction.
* Behaviour change: deposition now always uses the step-start reference, so
  results differ from Part 2 when stencils of several sources overlap.

## 9. Science snapshot schema

See `docs/output_schema.md` (CHUI Part 3 section). `/PartType0`:
`CHUI_ElementMassFractions [N,11]`, `CHUI_SpeciesMassFractions`,
`CHUI_DustMassFractions`, SGS, mu, tags and dwell. `/PartType4`: per-channel
return and energy `[N,3]`, birth element mixture `[N,11]`, birth time, CC/Ia
counts and stream key. Canonical datasets are unchanged; restart-only state is
not exported. No checkpoint schema change: the regrid adds no persistent
state, the chemistry bias lives in a module sidecar, and the wind table is
bound by its checksum.

## 10. Memory owners (governed admissions added)

`amrHydroGeometryCapacity(patch, layers)` (deep ghosts) ·
`compositionCandidateBytes` (regrid candidate sidecar + CSR plan) ·
`sources.element_diffusion.workspace` · `time.global_step_rollback.state_copy`
(full state copy) · distributed exchanges size their buffers by record
counts. No RSS, PSS or speedup numbers are claimed.

## 11. Tests and owner commands (all NOT RUN — USER-RESERVED)

Unit (label `chui_part3`): `unit_amr_periodic_cross_patch`,
`unit_amr_production_regrid`, `unit_chui_part3_remote_sources`,
`unit_config_chui_part3`, `unit_treepm_amr_force_resolution` (skips 77
without FFTW), `unit_zoom_contamination`, `unit_chui_part3_adaptive_and_decks`.
Integration: `integration_chui_part3_science_snapshot` (skips 77 without
HDF5). Benchmark: `bench_chui_part3_kernels [cells]`.

```bash
cmake --preset cpu-only-debug
```

```bash
cmake --build --preset build-cpu-debug
```

```bash
ctest --preset test-cpu-debug --output-on-failure -L chui_part3
```

```bash
ctest --preset test-cpu-debug --output-on-failure
```

Then the HDF5/FFTW/MPI presets of `docs/build_instructions.md` (2 and 4 ranks)
for the snapshot, force and distributed paths.

Not yet written (open): MPI tests for 2/4-rank periodic AMR, distributed
composition, remote feedback and SF gradients. Also open: regrid restart and
migration round-trip, rollback bit-identity integration, end-to-end four-mode
runs and composition-regrid benchmarks.

## 12. Final scientific truth

* Rate coefficients are UNQUALIFIED (Part 2.5 gate retained). R18 and the 0.4 eV
  heat are unverified. Case-A vs case-B and the UVB/WSS09 mismatch refusals are unchanged.
* No wind table, TREECOOL, WSS09 or production SSP yields are bundled. Data
  paths fail closed without a checksum.
* Topology blockers: multi-rank metal diffusion, octets split across ranks at
  derefine, and recursive ≥3-level multirate gas/source time integration.
* Research features (associations, adaptive chemistry bias, predictive
  refinement) are EXPERIMENTAL/UNQUALIFIED.
* MHD, RT and a new AGN model remain separate future campaigns.
