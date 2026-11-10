# CHUÍ — post-patch-55 v10 source-closure ledger

Date: 2026-11-10. Baseline: `2d6caf9` (merge of PR #134, patch 55), clean
working tree, branch `claude/busy-turing-sgpocl`. Controlling inputs:
`chui_claude_code_patch55_v10_source_closure_prompt.md`,
`chui_patch55_source_audit.md`, `chui_patch55_v10_full_spec_gap_matrix.md`
(supplied by the owner). `chui_final_implementation_route_consolidation_v10.md`
and `chui_full_physics_memory_architecture_specification_v02.md` were **not**
supplied in this session and were not read; v10 requirements are taken from
the gap matrix's restatement.

**Executed:** source reading, `git diff`, `git diff --check`, text search,
static call-site review. **Not executed (owner request):** CMake, compiler,
linker, CTest, Python harnesses, static analyzers, MPI jobs, simulations,
benchmarks, scientific validation. Every status below is static source
evidence only; none is a build, runtime, convergence, conservation or RSS
claim. All edits are uncommitted in the working tree.

Status vocabulary: `implemented_static_unverified`, `already_resolved_in_head`,
`open_source_gap`, `blocked_external_data`, `conditional_requires_contract`,
`optional_research_deferred`.

## 1. Mandatory repairs

### P1-55-01 / H02 + A02 — dual-energy AMR ghost transfer — implemented_static_unverified

Defect reproduced at HEAD: `src/amr/amr_hydro_orchestrator.cpp` (sweep, old
lines ~2930–2966) set every non-physical ghost's auxiliary lane to
`rho_ghost * e(source_gas_cell_id)`, where the id came from one geometric
probe (`refreshTraceableGhostSourceIds`, edge filler) even when the ghost's
conserved state was a fine→coarse average or a temporal blend.

Contract now implemented (single operator for both lanes):

| Relation | Conserved state | Auxiliary `rho e` (same code path) |
|---|---|---|
| same level, local | source cell | source SoA aux lane (`sourceAuxiliary`) |
| coarse→fine, synchronized | injected coarse cell | same cell's aux |
| coarse→fine, temporal, local | `(1-a) U_s + a U_e` from history | `(1-a)(rho e)_s + a (rho e)_e` from history, same `a` (`temporalStateForCell`) |
| coarse→fine, temporal, remote | recursive driver interpolates `U` | driver interpolates `rho e` then divides by interpolated `rho` (`hydro_amr_runtime.cpp` `interpolate`) → filler product is the interpolated density |
| fine→coarse, local/remote/mixed | equal-volume mean over validated complete stencil | equal-volume mean of `rho e` over the **same** sorted selection (`averageSourceAuxiliary`) |
| sparse MPI | owner canonical-`e` conserved state | `rho_remote * e_remote` |
| physical (open/reflective) | `fillHydroBoundaryGhostCells` | `fillHydroBoundaryGhostAuxiliaryInternalEnergy` inside the filler |
| deep/edge closures, in-patch periodic wrap | copied/mirrored cell (+ momentum parity) | copied cell's aux (parity does not act on `rho e`) |
| any contributor with unknown `e` | — | NaN, counted (`auxiliary_internal_energy_ghosts_invalid`) |

Files / functions:

* `src/amr/amr_ghost_fill.cpp`: `validAuxiliary`, `sourceAuxiliary`,
  `averageSourceAuxiliary`, `cellAuxiliary`, `storeGhostAuxiliary`; every
  `storeCell` of a ghost in `fillAmrGhost` now pairs with `storeGhostAuxiliary`;
  `temporalStateForCell(..., double& auxiliary, ...)`; physical aux filled in
  `fillAmrHydroGhostCells`; history capture `historyEndpointForStateRow`
  (replaces `conservedForStateRow`), `assignStart/assignEnd` store aux.
* `include/cosmosim/amr/amr_ghost_fill.hpp`: contract comment, two counters,
  `AmrHydroSparseRemoteCell::internal_energy_code` semantics.
* `include/cosmosim/core/simulation_state.hpp`:
  `AmrTemporalBoundaryHistoryCellRecord::{start,end}_internal_energy_density_comoving`
  (NaN default, appended).
* `src/amr/amr_hydro_orchestrator.cpp`: `loadCanonicalAuxiliaryInternalEnergy`;
  target and local source SoAs carry the lane during ghost preparation (lease
  `amr.ghosts.local_source_preparation` at 7 instead of 6 doubles/cell);
  `prepared_ghost_auxiliary` in the step-start snapshot arena
  (`amrPreparedGhostWorkspaceBytes(patches, ghosts, dual_energy)` overload,
  +8 B/ghost + 64 B); sweep copies every ghost's aux from the snapshot; the
  ungoverned `remote_internal_energy` sorted vector and the one-id lookup are
  deleted.
* `include/cosmosim/amr/amr_hydro_orchestrator.hpp`: overload declaration,
  `dual_energy_unresolved_ghosts` semantics.
* `src/workflows/hydro_amr_runtime.cpp`: remote conserved state uses the
  owner's canonical `e` (same positive density / `1e-14` floor as
  `loadAmrHydroConservedStateInto`); remote interval interpolation in `rho e`.
* Persistence: restart **v25** (`src/io/restart_checkpoint.cpp`,
  `include/cosmosim/io/restart_checkpoint.hpp`): history group schema 3 with
  two new datasets; reader accepts v24 (restores NaN), group 1/2/3; schema
  validation requires the datasets at v25; admissibility check (NaN or finite
  > 0); legacy FNV and SHA-256 digests include the lanes at v25 (NaN digested
  canonically); post-write verification covers them; migration projection
  accepts v25. AMR migration wire **v5** (`migration_wire.{hpp,cpp}`, +16
  B/history cell, estimate constant 112→128). `time_coordinator.cpp`: the
  "pre-endpoint restart ⇒ PM refresh" rule is pinned to v24 so the bump does
  not alter v24 resumes. Docs: `docs/restart_checkpointing.md`,
  `docs/output_schema.md`, `docs/hydro_numerics.md`; release manifest and
  three test sources updated to the v25 name (not run).

Static source trace (canonical → ghost aux):
`GasCellSidecar::internal_energy_code[row]` → (a) local:
`loadCanonicalAuxiliaryInternalEnergy` → SoA aux of target/source real cells
→ `fillAmrHydroGhostCells` → `selectSourceCells` (same selection object as
`averageSourceState`) → `storeGhostAuxiliary(ghost_cell)` →
`prepared_ghost_auxiliary[slot]` → sweep `auxiliary[ghost_cell]` →
`HydroCoreSolver` dual-energy lane; (b) temporal: capture
`historyEndpointForStateRow` → record → restart v25 / wire v5 →
`temporalStateForCell` (same `alpha`, same record cell) → (a); (c) MPI:
`makePatchCellPayloadRecord.internal_energy_code` → `sparse_consumer`
(conserved from the same `e`) → recursive `interpolate` (`rho e` weights) →
`sourceAuxiliary` → (a). No ghost aux is read from `source_gas_cell_id`.

Reproducibility impact: temporal-history conserved endpoints now use canonical
`e` instead of the `pressure_code` cache (bug fix: alpha = 0/1 equals the state
the synchronized fill reads); remote ghosts likewise. Results change wherever
the cache disagreed with `(gamma-1) rho e`. Mid-subcycle v24 resumes are no
longer bitwise in dual-energy runs (aux NaN, counted).

### P1-55-02 / H01 + A02 — in-domain deep ghosts fail closed — implemented_static_unverified

Defect reproduced at HEAD: `fillAmrGhost` continued any empty-selection deep
layer from the previous layer (or mirrored) without a domain proof.

Implementation (`src/amr/amr_ghost_fill.cpp`): `closureBounds`,
`targetFaceIsPhysicalWall`, `domainDeclared`, `alignedDomainIndex`,
`validateTargetPhysicalFacesOnDomainWalls`, `fillClosureMappedGhost`,
`deepLayerStencilIndex`, `ghostLocationText`; `fillAmrGhost(..., mapped_probe)`.
A deep layer with no source is mapped through the physical closure only on a
non-periodic axis whose wall is proven by the declared domain
(`AmrHydroGhostFillPatch::domain_{origin,extent}_comoving`; production: the
collectively reduced patch bounding box, `globalAmrHierarchyDomain` in
`hydro_amr_runtime.cpp`, forwarded by `ProductionAmrHydroOptions`) or, without
a declared domain, by a target face that is itself `kPhysical`. Otherwise:
`kMissingSource`, `in_domain_missing_source_rejections`, exception with
patch/level/axis/side/layer/epoch. Periodic axes never close. A face with no
neighbour off a declared wall is rejected
(`physical_face_inside_domain_rejections`). Distinct counters: valid closures
(`deep_layer_extrapolated_ghosts` = proven open wall, `deep_layer_reflected_ghosts`,
`edge_stencil_physical_closures`), unsupported geometry
(`unsupported_closure_geometry_rejections`, status
`kRejectedUnsupportedClosureGeometry`), missing source, stale remote
(`stale_epoch_rejections`, unchanged). Collective failure: the synchronized
distributed advance is now wrapped and voted
(`rethrowCollectivePreparationFailure(..., "distributed AMR hydro advance")`)
before the reflux exchange; the recursive driver already reduces per stage.
Ghost-only semantics: closures create no faces/registers.

Static gate: the only remaining path that writes a deep layer without a source
is `fillClosureMappedGhost` with `outcome.closed == true`; it is reachable only
when an axis index lies beyond `ClosureAxisBounds` (declared domain wall or
own `kPhysical` face). A legal in-domain stencil (incl. multi-rank, periodic)
either resolves through `selectSourceCells` or throws.

### P2-55-03 / H01 + A02 — reflective thin patches — implemented_static_unverified

* Edge ghosts: the counted "mirror deeper than the patch → copy" downgrade is
  removed; the mirror index `2w-1-i` is resolved wherever it lands (target
  cell, or beyond the opposite face through the neighbour patch with any legal
  level relation, or a second reflection when that face is also a wall), with
  momentum reversed per reflected axis. `edge_stencil_degenerate_closures` is
  retired (never incremented, kept for ABI).
* Deep reflective physical layers of patches thinner than the mirror depth are
  no longer omitted (`appendDeepGhostLayers`): created with a provisional owner
  source and `AmrHydroGhostDescriptor::physical_closure_mapped = true`; the
  AMR filler overwrites them with the exact mirror.
* Unsupported: > 8 reflections per axis (domain narrower than the stencil),
  boundary kinds other than open/reflective/periodic, probe range overflow →
  explicit rejection at fill time (before the patch advance of that step).
  There is no separate configuration-time rejection: the geometry is decided
  per step by regrid, so the fill-time collective failure is the earliest
  point with the needed hierarchy.

### P2-55-04 / M01 — compact CTU ghost scratch — implemented_static_unverified

`HydroPatchGeometry::face_ghost_slot_count` / `faceGhostSlotExtent()` (AMR
builders set it to the face-adjacent shell; default = all slots). The CTU
ghost delta/mask vectors are sized by this prefix and indexed by the
unchanged ghost slot (no lookup table, no per-face allocation). Producer
throws `logic_error` if a face references a slot beyond the prefix; consumer
treats beyond-prefix slots as "no slab entry" (identical outcome to an empty
mask, including 1-cell-axis suppression). New `HydroScratchBuffers::
workspaceBytes/prepareForActiveSet` overloads take the face-ghost bound; the
AMR sweep workspace passes `max faceGhostSlotExtent()`; the old overloads keep
the all-ghost bound. Scale: from `144 B × (all ghosts)` to `144 B × (shell
slots)`; for PPM with edge block roughly ≥ 2/3 of that scratch removed
(estimate from slot counts, not measured).

Also fixed while reviewing ownership (M01): `selectSourceCells` returned its
pmr selection **by value** (copied through the default memory resource on
every ghost and discarded); it now returns `void`.
`CompositionFluxAccumulator::resolveRow` gains an O(1) slot fast path
(was O(ghosts) per face).

## 2. H01–V02 coverage

| ID | Status | Evidence / remaining |
|---|---|---|
| H01 | implemented_static_unverified | P1-55-02, P2-55-03, P2-55-04 above. Numerical order unverified. |
| H02 | implemented_static_unverified | P1-55-01 above. |
| C01 | blocked_external_data | Fail-closed WSS09/UVB pairing gate present (`physics.metal_cooling_accept_uv_background_mismatch` default false; `metal_cooling_table.cpp:304`). HM01-consistent paired UVB/cooling table outstanding. No code change. |
| C02 | blocked_external_data | Loaders/provenance unchanged; CI/RR/DR and metal photo-rates need an approved data contract. |
| C03 | already_resolved_in_head | Adaptive fidelity in `thermochemistry_operator.cpp`; no source gap found. |
| C04 | optional_research_deferred | Chemistry shadow present; standalone higher-order shadow hydro not present. Proposed contract: shadow hydro advances a sampled patch subset with an alternative reconstruction from the same step-start ghost snapshot, never writes canonical state, reports L1 deltas. |
| I01 | already_resolved_in_head | `effective_multiphase_ism.cpp`; calibration/double-cooling qualification open (scientific). |
| I02 | already_resolved_in_head | `ism_dust_sgs.cpp`; coefficients unqualified. |
| S01 | already_resolved_in_head | Boundary-combination certification is a test item. |
| S02 | blocked_external_data | Yield coverage/IMF provenance. |
| F01 | already_resolved_in_head | Cross-level deposition conservation needs numerical qualification. |
| F02 | blocked_external_data | Single-Z / Kroupa SB99 winds vs selected Chabrier SSP. |
| D01 | conditional_requires_contract | Shattering/coagulation need a frozen grain-collision model and rates; none invented. Composition ghost origins: single id is consumed only on fine faces (single coarse source) and on coarse register faces whose update is suppressed, so the multi-source issue of P1-55-01 does not reach composition. |
| A01 | already_resolved_in_head | Duplicate planner / admission-cap fixes preserved (not touched). |
| A02 | implemented_static_unverified | Ghost consistency + fail-closed + collective vote above. MPI diagonal/periodic completeness remains a qualification risk (owner MPI tests). |
| A03 | already_resolved_in_head | Indicators unchanged. |
| A04 | already_resolved_in_head | Predictive/cost-aware refinement unchanged. |
| G01 | conditional_requires_contract | TreePM + force-resolution policy authoritative; no mandatory finer-force contract supplied. |
| T01 | already_resolved_in_head | Single time authority; default `coarse_endpoint` documented lag; `physical_event_barrier` exact via global clip. Traced: forecast (`discreteEventStepEndLimitTimeCode`) covers legacy/cursorless records with catch-up advanced on a copy, so publication stays exactly-once on the ordinal cursor; rollback re-derives from restored state; cursors persist in restart. |
| T02 | optional_research_deferred | Level-local SF/SSP/SN/wind publication not implemented; global barrier satisfies the event-time contract (performance limitation under high event rates, unmeasured). Optional IMEX/SDC absent. |
| M01 | implemented_static_unverified | CTU scratch compaction, aux-lane admission (prep lease, snapshot arena), ungoverned remote-e vector removed, per-ghost selection copy removed. FFTW/HDF5/MPI library buffers remain estimate-only by design. |
| M02 | already_resolved_in_head | CPU serial/OpenMP/MPI paths; no GPU dependency added. |
| V01 | implemented_static_unverified | Restart v25 / wire v5 versioned with compatibility (above). End-to-end qualification deferred. |
| V02 | already_resolved_in_head | No change. |
| X01 | not a v10 gap | No MHD/RT/AGN added. |

## 3. Interface / config / restart / memory impact

* Public headers (append-only fields; designated initializers keep meaning):
  `amr_ghost_fill.hpp` (2 patch fields, 5 counters), `amr_hydro_geometry.hpp`
  (1 enum value, 1 descriptor field), `amr_hydro_orchestrator.hpp` (2 option
  fields, 1 overload), `hydro_core_solver.hpp` (1 geometry field + accessor,
  2 overloads), `simulation_state.hpp` (2 history fields),
  `restart_checkpoint.hpp` (v25).
* No new configuration key.
* Restart v25 and AMR migration wire v5 (mixed-version ranks are refused by the
  existing wire-version check).
* Memory: +1 SoA lane for ghost-preparation target/sources and +8 B per
  prepared ghost when dual energy is on; −144 B per stencil-only ghost slot of
  CTU scratch; +16 B per temporal-history cell (retained, admitted through
  `historyCandidateBound` via `sizeof`).

## 4. Behaviour changes owners should expect in tests (not run)

* Non-periodic hierarchies without a declared domain whose deep layer at an
  interface face leaves every patch now throw instead of continuing.
* A face with no neighbour inside a declared domain throws (missing interface
  metadata). Non-rectangular root tilings fail closed.
* Reflective thin patches gain deep physical layers (ghost counts change).
* Temporal history / remote states use canonical `e` (values change where the
  pressure cache was stale).

## 5. Owner validation floor (not run)

`cmake --preset cpu-only-debug`, build, `ctest` at least:
`test_amr_temporal_ghost_fill`, `test_amr_patch_boundary_consistency`,
`test_restart_equivalence_amr_temporal_ghosts`, `test_amr_periodic_cross_patch`,
`test_amr_production_hydro_integration`, `test_hydro_ctu`,
`test_restart_checkpoint_schema`, docs/release cmake checks, plus the MPI
preset (`test_directed_amr_exchange_fault_injection_mpi`, distributed
recursive AMR). New tests worth adding: multi-source fine→coarse aux mean;
temporal aux interpolation; v24→v25 restart NaN restore; thin reflective
patch mirror via neighbour; in-domain missing deep source rejection.

## 6. Correction of earlier labels

`docs/repair/chui_post_patch53_source_closure.md` Campaign 2 listed the
reflective-mirror degeneracy as the "only counted downgrade" and the
fine→coarse edge aux as NaN. Both statements are superseded by §1 above
(the degeneracy is removed; fine→coarse aux is the conserved-consistent mean).
The deep-layer continuation it described is now reachable only beyond a
proven wall.
