# CHUÍ — post-patch-53 baryonic source-closure ledger

Date: 2026-11-10. Baseline: `6070b88` (merge of PR #133, patch 53), clean
working tree. Branch `claude/modest-wozniak-3zeffi`. Controlling inputs:
`chui_claude_code_source_closure_prompt.md` and
`chui_baryonic_source_audit_patch53.md` (supplied by the owner).
`chui_final_implementation_route_consolidation_v10.md` and
`chui_full_physics_memory_architecture_specification_v02.md` were **not**
supplied in this session and were not read.

**Executed:** source inspection, `git log/blame/show/diff`, text search.
**Not executed (owner request):** CMake, compiler, linker, CTest, Python
harnesses, lint, MPI jobs, simulations, benchmarks, scientific validation.
Every status below is static source evidence only; none is a build, runtime,
convergence, conservation or RSS claim.

Root cause of P0-HYDRO-01, P0-AMR-01 and P1-AMR-02: merge `efe31de`
("Merge branch 'main' into codex/chui-integral-source-closure-20261104")
combined `afc5c9d` (CTU ghost-slab completion, candidate-ranked V02 planner,
memory-aware regrid budget) with `1727c24` ("hydro patch 1": CTU boundary
proxy, two-pass mandatory-first scheduler, raw policy initialization) by
textual concatenation.

## Status table

| ID | Status | Changed files / functions | Source evidence | Unresolved dependency |
|---|---|---|---|---|
| P0-HYDRO-01 | implemented_static_unverified | `src/hydro/hydro_core_solver.cpp`: `transverseSourceCell` (3-arg, real/periodic only), new `boundaryProxySourceCell`, CTU lambdas `sumRealCellTransverse`, `halfTransverse`, `correctSideState` in `advancePatchActiveSetBatchedWithScratch`; final-pass CTU loop; header comment of `ctu_transverse_boundary_proxy_states` | Orphaned `bool& used_boundary_proxy) -> bool {`, the 3- vs 4-argument `transverseSourceCell` calls and the duplicated `if (halfTransverse(...))` fragments are removed. One provenance per side state: real/periodic → complete tangential ghost slab (`ctu_transverse_ghost_corrected_states`) → imported-ghost owner proxy (`ctu_transverse_boundary_proxy_states`) → normal-only (`ctu_transverse_degraded_states`). Inadmissible corrections keep the normal state (`ctu_transverse_inadmissible_states`). Hancock normal prediction + one transverse correction (standard CTU), no double predictor. Fluxes/registers/mass-flux sink (composition, chemistry, dust advection) are untouched. | Compile; `test_hydro_ctu`; boundary convergence tests (owner). |
| P0-AMR-01 | implemented_static_unverified | `src/amr/amr_production_regrid.cpp`: `planProductionAmrRegrid` (deleted `schedule_pass(true/false)`; duplicate-local-patch-id guard; mandatory max-level/balance blocks are counted, not thrown); `include/cosmosim/amr/amr_production_regrid.hpp` contract comment | One candidate planner visits each local leaf once ⇒ ≤ 1 refine op per patch; mandatory first by patch id and exempt from `max_operations` and V02 ceilings; optional ranked by V02 benefit/cost; `required_patch_ids`/`required_gas_cell_ids` are summed from the final `plan.operations`, so ID windows match emitted ops. Mandatory blocks go to `mandatory_blocked_by_{max_level,balance}`, which `runAutomaticAmrRegrid` reduces and fails collectively (no rank-local throw before collectives). Mandatory below gas softening with `limit_optional_reject_mandatory` still throws inside the runtime's guarded planning phase (collective rethrow). Hysteresis lives in `patchDecision` (unchanged). | Compile; `test_amr_production_regrid` (its `capped.mandatory_blocked_by_max_level == 1` expectation now matches the counted policy). |
| P1-AMR-02 | implemented_static_unverified | `src/workflows/hydro_amr_runtime.cpp`: `runAutomaticAmrRegrid` | The stale second policy initialization (raw `amr_regrid_max_operations`, `current_scale_factor`, duplicated softening) is deleted; `policy.max_operations` keeps `memoryAwareOptionalRegridBudget(...)`, equal to `m_last_regrid_optional_budget` and to the published `optional_operation_budget`. New `regrid_scale_factor` = `timeline_step.scale_factor_end` (cosmological comoving frame, else 1): the regrid runs after the hydro stage advanced the rows to the step end and before `CosmologicalTimeline::commitStep` (which is when `current_scale_factor` becomes `a_end`); it now feeds the Jeans criterion, cooling-time physical density/redshift and the predictive-feedback drift scale. V02 byte budget (`optional_control.memory_budget_bytes`) is still enforced in the planner alongside `max_operations`. The memory-aware count is rank-local by design (per-rank op budget; 2:1 balance does not need cross-rank agreement on optional subsets). | Compile; regrid telemetry check (owner). |
| P1-HYDRO-02 | implemented_static_unverified (partial) | `src/amr/amr_ghost_fill.cpp` `fillAmrGhost` deep-layer branch; `include/cosmosim/amr/amr_ghost_fill.hpp` (`AmrHydroGhostFillPatch::physical_boundary_kind`, `AmrHydroGhostFillDiagnostics::deep_layer_reflected_ghosts`); `include/cosmosim/amr/amr_hydro_geometry.hpp` (`kFilledDeepLayerReflected`); `src/amr/amr_hydro_orchestrator.cpp` (target view sets the boundary kind; diagnostics summed) | Deep stencil layers beyond a domain wall behind a thin neighbour now use the physical closure: reflective ⇒ exact mirror of column position `2m+1-L` (column ghost / owner / interior cell) with reversed normal momentum; open ⇒ zero-gradient (the outflow closure itself). CTU edge/corner: see P0-HYDRO-01 provenance chain. | **Open code work:** (1) true corner/edge ghost provenance — ghost cells displaced along two axes are not part of `HydroPatchGeometry`, so a complete slab at face edges needs new diagonal ghost storage + ghost-fill + MPI interface topology (`exchangeDirectedHydroGhosts` ghost_layers/corner requests); (2) uncovered-but-inside-domain deep probes (remote payload depth) still take the counted continuation. Both remain counted downgrades, never labelled full-order. |
| P1-TIME-01 | implemented_static_unverified (partial) | `include/cosmosim/core/config.hpp` (`SspEventPublicationModeConfig`, `k_ssp_event_barrier_tolerance_years`, `PhysicsConfig::ssp_event_publication_mode`, `sspEventPublicationModeConfigToString`); `src/core/config.cpp` (parse, default, validation, normalized dump); `src/io/restart_checkpoint.cpp` (pre-key checkpoint projection: only default accepted); `src/workflows/source_runtime.cpp` (`physicalEventBarrier`, `enforcedEventLagBoundYears`, `discreteEventStepEndLimitTimeCode`, `DiscreteEventChronologyLedger::physical_event_barrier/lagPolicy`, stellar-population stage, telemetry); `include/cosmosim/workflows/source_runtime.hpp` comment | New key `physics.ssp_event_publication_mode = coarse_endpoint | physical_event_barrier`. `coarse_endpoint` is byte-identical in behaviour to patch 53 (lag key semantics unchanged). `physical_event_barrier` reuses the guarded forecast + `allreduceMinDouble` time-authority clip with a fixed 1-yr resolution, so global coarse steps (legal recursive-AMR/multirate synchronization points) end at the next physical event; forecast-covered violations fail closed. Events remain exactly-once (threshold/ordinal cursor unchanged); rollback re-derives the forecast from restored state; restart persists cursors. | **Open code work:** fine-level (intra-coarse-interval) event barriers on recursive level clocks are not implemented; exact mode shortens the global step instead. Legacy-record (`unforecast`) events remain counted/exempt. Owner tests: config round-trip, restart projection, event-time lag ≤ 1 yr. |
| P2-MEMORY-01 | implemented_static_unverified | `src/physics/supernova_feedback.cpp` + header: `nearestInto`, phase-reused `localCandidatesWithin` overload; `src/workflows/source_runtime.cpp` remote halo responses; `src/workflows/migration_balance_runtime.cpp` `exchangeAndValidateAmrPatchPayloads` (new `memory_governor` parameter, both call sites) | Halo responses reuse one k-capacity heap and row buffer admitted in `sources.feedback.remote_halo_responses` (no per-query allocation). Migration validation `std::set` → flat vector reserved per peer in the exchange's guarded admission callback, governor lease `migration.amr_validation.remote_cell_keys` (`kDiagnostic`, old+new coexistence); consumer appends within capacity; duplicates detected after the exchange in the voted decode phase. | FFTW/HDF5/MPI library buffers stay estimate-only by M01 design. |
| G01 | not_required_conditional | — | `planProductionAmrRegrid` already applies `amr_force_resolution_policy` (limit optional / reject mandatory, or counted reference). No mandatory source contract in the supplied inputs requires a new refinement-aware gravity correction. | Becomes required only if the flagship resolution makes gas cells finer than gas softening in mandatory regions. |
| D01 | not_required_conditional | — | Sputtering/accretion/destruction/condensation lanes exist (`src/physics/ism_dust_sgs.cpp`). No accepted shattering/coagulation data/specification was supplied; none invented. | Authoritative rate data + specification. |
| C02 / UVB / yields | blocked_missing_contract (data) | — | Loaders and provenance gates unchanged. | External table qualification and rights. |

## Correction of earlier claims

`docs/repair/chui_source_implementation_closure.md` labelled patch 53
`SOURCE-COMPLETE, EXTERNAL-DATA-INCOMPLETE`. Inspected HEAD source
contradicted that label (malformed CTU translation unit, duplicate regrid
scheduling, overwritten memory-aware budget). The label is superseded by this
ledger; after these edits the accurate label is **source repaired
(static, unverified), with open P1-HYDRO-02 corner-ghost and P1-TIME-01
fine-level event-barrier work**.

## Interface / config / restart impact

* Public headers changed: `amr_ghost_fill.hpp` (appended field + counter),
  `amr_hydro_geometry.hpp` (appended enum value), `config.hpp` (new enum,
  constant, field, serializer), `supernova_feedback.hpp` (new overload and
  private helper; old overload kept), `hydro_core_solver.hpp`,
  `amr_production_regrid.hpp`, `source_runtime.hpp` (comments only). All
  additions are appended; existing designated initializers and callers keep
  their meaning.
* New config key `physics.ssp_event_publication_mode` (default
  `coarse_endpoint`): typed parse, validation, normalized dump, restart
  projection, `docs/configuration.md`, `docs/time_integration_event_chronology.md`,
  `docs/restart_checkpointing.md`. Normalized config text gains one line, so
  the config hash of every run changes (reproducibility impact: provenance
  only; behaviour unchanged at the default).
* Behaviour changes at defaults: CTU (fixes a non-compiling merge; proxy now
  only after an incomplete slab), regrid (single planner; memory-aware budget
  actually applied; indicators at `a_end`), reflective deep ghost layers
  (mirror instead of copy). These are bug fixes, not new physics.
* No checkpoint/snapshot dataset changed.

## Owner validation floor (not run)

`cmake --preset cpu-only-debug`, build, and at least `test_hydro_ctu`,
`test_amr_production_regrid`, `test_amr_recursive_subcycling`,
`test_config_chui_part3`, restart-config projection tests, the MPI preset for
the migration validation exchange and the remote feedback halo, plus a
`physical_event_barrier` config/restart round trip.

## Continuation state

All P0/P1/P2 edits above are in the working tree (uncommitted at the time of
writing). Next exact edits if continued: (1) P1-HYDRO-02 diagonal (edge)
ghost storage in `amr_hydro_geometry.cpp` + fill in `amr_ghost_fill.cpp` +
interface requests in `src/parallel/distributed_memory_amr.cpp`; (2)
P1-TIME-01 level-clock event barrier in `source_runtime.cpp`
`executeAmrLevelSmoothSources` path; (3) owner-side test updates for the new
config key and diagnostics.

## Campaign 2 (2026-11-10): P1-HYDRO-02 edge/corner infrastructure, P1-TIME-01 follow-up

Static source edits only (uncommitted). No build, test or run.

| ID | Status | Changed files / functions | Source evidence | Remaining |
|---|---|---|---|---|
| P1-HYDRO-02 | implemented_static_unverified | `hydro_core_solver.hpp` (`HydroGhostCell::transverse_axis/transverse_offset`, `HydroPatchGeometry::edge_ghost_slot_begin/edge_ghost_depth/edgeGhostSlot`, `hydroEdgeGhostCount`); `hydro_cartesian_patch.cpp`; `hydro_core_solver.cpp` (`tangentialGhostCell` edge lookup, incomplete-slab `logic_error` with the block present, proxy only without the block); `amr_hydro_geometry.{hpp,cpp}` (`AmrHydroGeometryOptions::edge_stencil_depth`, `AmrHydroGhostDescriptor::edge_stencil`, capacity, `appendEdgeStencilGhosts`); `amr_ghost_fill.{hpp,cpp}` (transverse-aware probe/bounds, `selectSourceCells(force_any_relation)`, `fillAmrGhost` probe override, `fillEdgeStencilGhost`, column search excludes edge ghosts, three counters); `amr_hydro_orchestrator.{hpp,cpp}` (`ProductionAmrHydroOptions::ctu_edge_stencils`, `productionAmrHydroEdgeStencilDepth`, `edgeStencilSourceTouches`, edge-aware `ghostProbePoint` and source tracing, diagnostics); `distributed_mesh.hpp` + `distributed_memory_amr.cpp` (`edge_stencil_layers`, widened transverse mask, reach scale); `hydro_amr_runtime.cpp` (CTU flag + topology) | Stencil-only edge block (no faces/registers); tangential slab complete at patch edges/corners for same-level, coarse/fine (temporal), fine/coarse, periodic (incl. in-patch wrap) and MPI sources; physical domain closure per crossed axis; missing data fails closed; proxy/normal-only no longer reachable for AMR CTU geometries. Admission: geometry capacity includes the block; ghost preparation/snapshot sizes follow `ghost_cells`. | Counted degeneracy only: reflective mirror depth > patch extent (`edge_stencil_degenerate_closures`). Dual-energy auxiliary of an edge ghost whose source is a fine->coarse average stays NaN (E - K, counted as for face ghosts). Owner tests: CTU smooth-wave convergence across edges/corners, MPI diagonal-rank delivery, periodic and reflective corners. |
| P1-TIME-01 legacy records | implemented_static_unverified | `source_runtime.cpp` (`legacyStellarPopulationRecord`, forecast covers missing/cursorless records with catch-up advanced, stage enforces every timed event) | Forecast and stage share one deterministic legacy record. | — |
| P1-TIME-01 level-local publication | open | — | Exact publication is achieved by global coarse-step clipping (barrier resolution 1 yr). | Level-local SSP publication inside a recursive coarse interval: per-star contiguous cursors across level substeps, level-restricted feedback stencils (cells of level >= L), deposition after the level temporal-history end capture with open flux registers, and level-local rollback. Not implemented. |
| G01 | not_required_conditional | — | unchanged | needs a mandatory resolution contract |
| D01 | not_required_conditional | — | unchanged | needs authoritative rate data/specification |
