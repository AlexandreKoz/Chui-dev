# CHUÍ — post-patch-56 v10 final source-implementation closure ledger

> **Superseded in part (2026-11-11):** §2 (volume-based domain check), the
> S02/F02 rows of §3–§4 and the "unsettled model decisions" of §6 are
> superseded by
> [`chui_post_patch58_final_source_closure.md`](chui_post_patch58_final_source_closure.md)
> (exact leaf-tiling proof; explicit yield-coverage and wind-compatibility
> policies). This record is kept unchanged below for history.

Date: 2026-11-10. Baseline: `4775f35` (merge of PR #135, patch 56), branch
`claude/gifted-cray-zjay8f`, clean working tree at start. This ledger extends,
and does not delete, `chui_post_patch55_source_closure.md`.

Normative inputs read in full this session (unlike patch 55, the originals
were supplied): `chui_final_implementation_route_consolidation_v10.md`,
`chui_full_physics_memory_architecture_specification_v02.md`,
`chui_consolidated_physics_performance_plan_v01.md`,
`CHUI_Baryonic_Physics_AMR_Audit_2026-10-30.md`; leads:
`chui_patch56_source_review.md`, `chui_patch55_v10_full_spec_gap_matrix.md`,
`chui_data_staging.md` (repository copy is current).

**Executed:** source reading, text search, `git status/diff/diff --check`,
manual call-path and reader/writer review. **Not executed (owner request):**
CMake, compile, link, CTest, Python harnesses, static analyzers, MPI, HDF5
round trips, simulations, benchmarks, scientific validation. Every status is
static source evidence only. Nothing is committed.

Status vocabulary: `implemented_static_unverified`, `repaired_in_this_patch`,
`open_source_gap`, `blocked_external_data`, `conditional_requires_contract`,
`optional_research_deferred`, `not_in_scope`.

## 1. P0 — restart v25 group-schema regression — repaired_in_this_patch

Reproduced at HEAD (`src/io/restart_checkpoint.cpp`):
`writePendingFluxRegisterGroup` wrote `@schema_version = 3` (patch 56 changed
it from 2 with no layout change); `writeAmrTemporalBoundaryHistoryGroup` still
wrote `2` while emitting the v25 auxiliary datasets; the reader demanded `>= 3`
for v25. Every v25 file written by patch 56 was unreadable by its own reader.

Repair (same file):

* Named constants `k_pending_flux_register_group_schema = 2`,
  `k_temporal_history_group_schema_legacy_max = 2`,
  `k_temporal_history_group_schema_v25 = 3`.
* `writePendingFluxRegisterGroup` writes 2; `writeAmrTemporalBoundaryHistoryGroup`
  writes 3 unconditionally (also for an empty history).
* `readAmrTemporalBoundaryHistoryGroup`: outer v25 ⇒ group must equal 3;
  outer v18–v24 ⇒ group 1 or 2; any other pairing is refused (was: any of
  1/2/3 for any outer version). The auxiliary datasets are read iff group = 3.
* `readPendingFluxRegisterGroup`: now checks its attribute (when present):
  1..2 accepted, 3 refused. The validator already requires the attribute at
  v18+.
* `validateRestartCheckpointSchema`: at v25 requires pending = 2 and
  temporal = 3 (preflight, before any payload read).
* `verifyRestartCheckpointHdf5` (write-time readback): asserts both group
  attributes with `require_u32_attr`.

Not changed (deliberately): global restart version (v25), dataset layout,
digests (`restartPayloadIntegrityDigestsImpl` with `include_dual_energy_history_v25`
— group attributes were never digested), `isRestartSchemaCompatible`, the v24
read path, migration projection (v23/v24/v25), AMR migration wire v5
(`src/workflows/internal/migration_wire.cpp` lines ~514/535 encode/decode the
two lanes symmetrically; `k_amr_patch_record_wire_version = 5` with exact
refusal at ~790), release manifest (`cosmosim_restart_v25`).

Writer/reader field matrix (static):

| Outer | pending group | temporal group | auxiliary datasets | reader outcome |
|---|---|---|---|---|
| v25 (this writer) | 2 | 3 | present | read; validator, reader and verifier agree |
| v25 (patch-56 writer) | 3 | 2 | present | refused by validator (pending≠2) — never coerced |
| v24 / v22–v23 | 2 | 1–2 | absent | read; aux restored as NaN |
| v18–v21 | 2 (metal lanes absent) | 1–2 | absent | read; aux NaN; metal lanes 0 |
| v24 with temporal 3 | — | 3 | — | refused (inconsistent writer) |
| < v18/v19 | optional | absent | — | unchanged legacy path |

Legacy v24 NaN semantics re-traced: `readAmrTemporalBoundaryHistoryGroup` fills
NaN; `validateAmrTemporalBoundaryHistoryForRestart` admits NaN or finite > 0;
`amr_ghost_fill.cpp` `temporalStateForCell` propagates NaN; the AMR sweep
(`amr_hydro_orchestrator.cpp` ~2992–3011) copies it to the solver lane and
counts it in `dual_energy_unresolved_ghosts`; the solver's
`selectDualEnergyThermalState` then uses `E − K`. The ghost fallback is
therefore counted (AMR diagnostic, reported by `hydro_amr_runtime.cpp`); the
solver's `dual_energy_unknown_auxiliary_cells` counts real cells only.
`docs/hydro_numerics.md` had a stale pre-P1-55-01 sentence ("not counted per
ghost") — corrected. Mid-subcycle v24 resume is not bitwise in dual-energy runs
(documented).

Docs: `docs/output_schema.md` (stray `version = 24` bullet → `version = 25`;
per-group identity table), `docs/restart_checkpointing.md` (group-identity rule).

## 2. New mandatory source gap fixed — explicit AMR domain contract — repaired_in_this_patch

Finding: `globalAmrHierarchyDomain` (`src/workflows/hydro_amr_runtime.cpp`)
inferred non-periodic walls from the populated-patch bounding box; holes /
non-rectangular tilings were rejected only implicitly, at the first ghost fill
that met a neighbourless in-domain face. The authoritative hydro domain owner
was traced: the root tiling is built from the gas-cell layout
(`initial_condition_runtime.cpp` `materializeRootHydroPatchIfMissing` →
`buildCartesianGasCellRowLayout`, origin = first cell − dx/2), *not*
`cosmology.box_size_*` (gravity/PM extent; isolated modes may place gas in a
sub-box). Refinement replaces a parent by its leaf children
(`refineProductionPatchInSimulationState`), so populated patches are disjoint
leaves. The bounding box of the leaves therefore *is* the root tiling iff
coverage is complete.

Implementation: `globalAmrHierarchyDomain` now also sums leaf volumes,
`allreduceSumDouble`, and throws (identically on every rank, before any ghost
fill) when covered volume ≠ bounding-box volume (relative 1e-9). Supported
restriction stated explicitly: production AMR requires rectangular, hole-free
coverage. No new key; no interface change; one extra scalar collective per
production AMR hydro step. Docs: `docs/hydro_numerics.md`.

## 3. P2 — external data contract table

| Product / key | Model assumptions | Allowed coverage | Provenance / hash | Rights | Gate in source | Execution status |
|---|---|---|---|---|---|---|
| C01 WSS09 metal cooling — `physics.metal_cooling_mode=element_table*`, `metal_cooling_table_path/_sha256` | embeds HM01 UVB; CIE-free photoionized tables | tabulated z 0–8.989, T, n; `metal_cooling_table_range` | `4bbca8a9…3410`, bit-exact vs source HDF5 | no licence stated | `MetalCoolingTable::checkUvBackgroundCompatibility` refuses declared-UVB ≠ configured TREECOOL hash unless `metal_cooling_accept_uv_background_mismatch=true`; element tables require a Part-2 network, which refuses analytic hm12/fg20 (`config.cpp` ~1853) — gate complete | research-only (no HM01 TREECOOL) |
| FG20 TREECOOL ×0.68 — `uv_background_model=treecool_table` | FG20 H/He photo-rates only | z grid of file | `4bfc5d9b…b313`, mirror only (original host unreachable) | no | every consumer behind `thermochemistry_allow_unqualified_rates` (`config.cpp` ~1836, `thermochemistry.cpp` ~120) | research-only |
| C02 C/O/Si/Fe ions — `chemistry_metal_ion_elements`, `metal_ion_rate_table_path/_sha256` | 1e4–1e8 K, optically thin, collisional; no metal photo-rates, no low-T DR | per-file | `633a66ee…2d00`, CI/RR/DR not independently verified | no | SHA-256 check (`baryonic_runtime_config.cpp`), `declaredVerified()` + opt-in, coverage of selected elements, ionization-energy contract (`thermochemistry.cpp` ~133–160) | research-only |
| S02 NuGrid+W7 yields — `stellar_population_model=chabrier_yield_table_v1`, `ssp_yield_table_path/_sha256` | per-star net yields, IMF integrated by CHUI Chabrier (no IMF conversion needed) | NuGrid 1–25 M☉, Z 1e-4–0.02; W7 single Z | `e2161646…f8e14` | no | citation/licence mandatory; edge-row clamp counted (`table_clamps`); **new**: configured-vs-tabulated coverage in `sources.data_assumptions` | diagnostic (safe-default fragment), not qualified |
| F02 SB99 winds — `ssp_wind_model=table_v1`, `ssp_wind_table_path/_sha256` | Kroupa IMF (CHUI SSP: Chabrier); single Z=0.014 | ages in file; Z clamped and counted | `ba232a1c…5da5` | no | provenance header, monotone budgets, Z clamp counted, limitation lines surfaced; no machine-readable IMF field | diagnostic, not qualified |
| HM01 TREECOOL | — | — | absent | — | — | unavailable |
| D01 grain shattering/coagulation | no approved law/rates | — | none | — | no operator exists (none invented) | unavailable |

Blocked external data: HM01-consistent TREECOOL (or a WSS09-equivalent table
computed with FG20); independently verified C02 CI/RR/DR + metal photo-rates;
licensed/qualified yields covering the configured CC range; multi-Z SN Ia;
Chabrier-consistent multi-Z wind tables; any grain-collision rate set.

Unsettled model decisions (not source gaps): whether yield edge clamping
(e.g. CC 25–40 M☉ using the 25 M☉ row) is acceptable or must be refused; a
machine-readable IMF field in `chui_ssp_wind_table_v1` and its gate.

## 4. H01–V02 matrix

| requirement | normative clause/file | implementation entry point | source-state owner | MPI/AMR/time/restart path | status | remaining code | external data / conditional decisions | deferred verification |
|---|---|---|---|---|---|---|---|---|
| H01 MUSCL–Hancock/CTU/PPM, HLLC+HLLE, passive conservation | v10 §2, §4.1; audit F1/F3 | `HydroCoreSolver::advance` (`hydro_core_solver.cpp`), `hancockHalfStepIncrement` (`hydro_reconstruction.cpp`, sign and `dt/(a dx)` re-verified), lane-generic HLL (`hydro_riemann.cpp`), `fillAmrGhost`/`fillClosureMappedGhost` (`amr_ghost_fill.cpp`) | canonical gas conserved SoA | edge/corner ghosts, diagonal MPI requests (`distributed_memory_amr.cpp`), fail-closed deep ghosts, collective vote; **this patch**: explicit domain-coverage contract | repaired_in_this_patch (domain contract); otherwise implemented_static_unverified | none found | — | convergence, positivity, diagonal/periodic MPI delivery |
| H02 dual energy, high-Mach, EOS/comoving | v10 §4.1, §3.1 | `selectDualEnergyThermalState`, auxiliary update (`hydro_core_solver.cpp`); ghost aux (`amr_ghost_fill.cpp` P1-55-01) | `GasCellSidecar::internal_energy_code` (aux lane derived per advance) | temporal history aux endpoints, restart v25 (**group schema repaired**), wire v5 | repaired_in_this_patch (restart identity) | none | — | Mach sweep, restart round trip |
| C01 hybrid H/He NEQ + metal cooling, UVB/shielding | v10 §4.2; audit F2 | `ThermochemistrySolver` (`thermochemistry.cpp`), `MetalCoolingTable` (`metal_cooling_table.cpp`), `radiation_background.cpp`, full-interval integrator (`cooling_heating.cpp`, `SourceIntervalIncompleteError`) | composition sidecar + gas sidecar | source transaction retry; restart via module sidecars | blocked_external_data | none (gate complete) | HM01-consistent UVB/table pairing | equilibrium/transient histories |
| C02 equilibrium/analytic + detailed reference with ions | v10 §4.2 | `chemical_network.cpp`, `metal_ion_network.cpp`, `thermochemistry_operator.cpp` | chemistry sidecar | same | blocked_external_data | none | verified CI/RR/DR, photo-rates, rights | solver accuracy |
| C03 adaptive fidelity, conservative transitions | v10 §4.2 (RESEARCH) | `thermochemistry_operator.cpp` adaptive selection, `gas_composition_sidecar` | regime tag + species history | restart sidecar; AMR migration payload | implemented_static_unverified | none | qualification of switching accuracy | matched-error comparisons |
| C04 sampled shadow chemistry/hydro | v10 §2 (RESEARCH/diagnostics), §3.2 | chemistry shadow in `thermochemistry_operator.cpp` (~417–560; capped by `shadow.max_cells`, writes only `ChemistryShadowReport`) | diagnostic report only | V02 controller publication (`source_runtime.cpp`) | implemented_static_unverified (chemistry); optional_research_deferred (shadow hydro) | shadow hydro absent (optional) | shadow-hydro reference scheme/metric/budget contract not frozen | shadow disagreement qualification |
| I01 effective vs resolved ISM | v10 §4.3; audit F7 | `effective_multiphase_ism.cpp` (volume-weighted ledger re-verified), `baryonic_source_stage.cpp` | ISM regime config + gas sidecar | hydro closure, CFL refresh | implemented_static_unverified | none | calibration | no-double-cooling tests |
| I02 evolved SGS turbulence | v10 §4.3 (RESEARCH) | `ism_dust_sgs.cpp`; recursive smooth sources | SGS lane | regrid transfer, restart sidecar | implemented_static_unverified | none | closure coefficients | energy transfer qualification |
| S01 stochastic SF, cross-boundary gradients | v10 §4.3; audit F5 | `star_formation.cpp`, `workflows/internal/star_formation_gradients.hpp` | star sidecar + ID registry | per-star elapsed time (`stellarElapsedYearsInInterval`) re-verified | implemented_static_unverified | none | — | patch/rank independence |
| S02 IMF/yields/events/associations | v10 §4.4 | `SspModel`/`SspYieldTable` (`stellar_population.cpp`); **this patch**: `channelCoverage` + coverage report | star sidecar + event cursors | restart checkpoint of SSP counters | blocked_external_data | none | licensed full-coverage yields; edge-clamp acceptance decision | chronology/restart equivalence |
| F01 SN thermal/mechanical/stochastic, refresh | v10 §4.4; audit F4 | `supernova_feedback.cpp`; `refreshGasThermodynamics` after source mutation (`source_runtime.cpp` ~2420, ~2895) | gas conserved + ledgers | cross-patch/rank deposition | implemented_static_unverified | none | — | conservation/resolution/MPI tests |
| F02 winds/pre-SN, event-aware refinement | v10 §4.4 | `StellarWindTable` (`stellar_wind.cpp`), wind/pre-SN double-count refusal (`config.cpp` ~2064) | event state | same as S02 | blocked_external_data | none | Chabrier-consistent multi-Z winds; machine-readable IMF field | budget invariance |
| D01 elements, advection/diffusion, dust/molecules | v10 §4.2, §3.1 | `baryon_composition.cpp`, `composition_transport.cpp`, `element_diffusion.cpp`, `metal_diffusion.cpp`, `ism_dust_sgs.cpp` | element/dust inventories | composition reflux/migration/restart | implemented_static_unverified; conditional_requires_contract (shattering/coagulation) | none mandatory | grain-collision law + rates | element/solid partition tests |
| A01 automatic regrid, conservative migration/restart | v10 §4.5; audit F6 | `amr_production_regrid.cpp`, `refine/derefineProductionPatchInSimulationState` | patch SoA + identity map | memory-admitted transaction; restart | implemented_static_unverified | none | — | real evolving-refinement runs |
| A02 recursive distributed subcycling, temporal ghosts, reflux | v10 §4.5 | recursive driver (`hydro_amr_runtime.cpp`), `amr_hydro_orchestrator.cpp` | temporal history + pending registers | restart v25 groups repaired; domain contract | repaired_in_this_patch (restart, domain); implemented_static_unverified otherwise | none | — | MPI diagonal/periodic, reflux |
| A03 physics-aware indicators | v10 §4.5 | `amr_production_regrid.cpp` indicators | regrid planner | — | implemented_static_unverified | none | force-resolution criterion tied to G01 | indicator qualification |
| A04 predictive/cost-aware refinement | v10 §4.5 (RESEARCH) | regrid planner + resource controller | planner | memory admission | implemented_static_unverified | none | — | fixed-budget comparisons |
| G01 TreePM/KDK + force/softening agreement | v10 §4.5 | `tree_pm_coupling.cpp`, `gravity_runtime.cpp` | particle state + force cache | v20/v24 force cache in restart | conditional_requires_contract | fine-level correction only if a force-error contract demands it | force-resolution contract | force-at-fine-level tests |
| T01 one time authority, Strang, retries | v10 §4.6 | `time_coordinator.cpp`, `stage_transaction.cpp`, `source_runtime.cpp` | integrator state + scheduler | collective retry/rollback; restart cursors | implemented_static_unverified | none | — | temporal convergence |
| T02 recursive multirate sources; optional IMEX/SDC | v10 §4.6 | `executeAmrLevelSmoothSources` (`hydro_amr_runtime.cpp`/`source_runtime.cpp`), `physical_event_barrier` publication | level clocks + event cursors | global event barriers (accepted semantics) | implemented_static_unverified; optional_research_deferred (level-local events, IMEX/SDC) | none mandatory | — | accuracy vs cost |
| M01 canonical state, phase scratch, node reservations | v02 §3, §6–7; v10 §5 | `node_memory_admission.cpp` (wired in `reference_workflow.cpp` ~394), `memory_governor.cpp`, CTU face-ghost prefix scratch | governor | regrid old/new coexistence admission | implemented_static_unverified | none found | node reserve policy is site-dependent | RSS/PSS/high-water measurements |
| M02 scalar/OpenMP/MPI CPU, optional GPU | v10 §1.2, §5.3 | core + `parallel/` | — | — | implemented_static_unverified | none | — | serial/OpenMP/MPI builds |
| V01 ledgers, restart, provenance, four modes | v10 §9 | `restart_checkpoint.cpp`, source ledgers, `sources.data_assumptions` | restart/provenance | restart v25 repaired | repaired_in_this_patch (restart identity) | none | — | end-to-end four-mode runs |
| V02 error/resource coordination | v10 §2 (RESEARCH) | `include/cosmosim/core/adaptive_resource_control.hpp`, V02 controller publication | versioned controller sidecar | restart sidecar | implemented_static_unverified | none | — | adaptive-policy comparisons |
| X01 new MHD / general RT / new AGN | v10 §2, §11 BOUNDARY | — | — | — | not_in_scope | — | — | — |

October-30 F1–F7 rechecked in current code: F1 (`hancockHalfStepIncrement`
negative-sign derivative added; `hydroTransportCourantFactors` = `dt/(a dx)`),
F2 (`SourceIntervalIncompleteError` → retry), F3 (lane-generic HLL), F4
(`refreshGasThermodynamics` after SF/feedback), F5
(`stellarElapsedYearsInInterval` with newborn birth tick), F6 (production
regrid + recursive driver), F7 (`delta_energy_density * cell_volume_comoving`)
— none reopened.

## 5. Interface / config / schema / reproducibility impact

* Public header (additive): `include/cosmosim/physics/stellar_population.hpp`
  `SspYieldTable::Coverage`, `SspYieldTable::channelCoverage(std::size_t)`. No
  migration needed (no existing signature changed).
* No new configuration key. No restart layout change; restart v25 writer now
  matches its reader. AMR wire v5 unchanged.
* Reproducibility: restart output differs only in group attributes; payload
  digests unaffected. `sources.data_assumptions` gains keys (diagnostic). The
  domain-coverage check can newly refuse a perforated hierarchy at step start
  that previously failed at its first offending ghost fill (earlier, clearer,
  collective failure; no physics change for valid hierarchies).

## 6. Remaining obligations (not source gaps)

* **Source blockers:** none identified for mandatory v10/v02 items.
* **External datasets:** §3.
* **Physical/model contracts:** D01 grain collisions; G01 fine-level force
  correction; yield edge-clamp policy; wind-table IMF field.
* **Optional research:** C04 shadow hydro; T02 level-local event publication;
  IMEX/SDC.
* **Deferred (owner):** build, CTest (`test_restart_checkpoint_schema`,
  `test_restart_equivalence_amr_temporal_ghosts`, AMR temporal/boundary
  tests, HDF5 preset), MPI preset, four-mode runs, scientific validation.
  Suggested new test sources: v25 write→read with empty and non-empty temporal
  history asserting group attributes 2/3; synthetic v25 file with pending=3
  refused; perforated two-patch hierarchy refused by `globalAmrHierarchyDomain`.

**Not `SOURCE-COMPLETE` in any validated sense:** no compile or test was run.
