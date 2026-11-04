# CHUI design-conformance repair ledger — 2026-11-04

Execution contract: **SOURCE ONLY. ALL EXECUTABLE VALIDATION USER-RESERVED.**
No configure, build, compiler, test, sanitizer, static-analysis executable,
benchmark, converter, simulation or CI workflow may be executed in this campaign.
Authored tests are **NOT RUN — USER-RESERVED**. Source closure never means
numerical qualification or physical calibration.

Starting live `main`: `80181a901b7eaece921c64eb48a09040ad818b92`, verified
through the GitHub commit API on 2026-11-04 and identical to the clean local HEAD.
Branch: `codex/chui-conformance-milestone-a-20261104`, rooted at that exact SHA.
The uploaded audit is evidence, not execution instructions. Design authorities
read: implementation route v1.0, physics/performance plan v0.1 and full-physics
memory architecture v0.2, supplied with the design bundles. Their historical
execution suggestions do not override this session's prohibition.

## Requirement-to-source map established before source edits

Paths below are relative to the repository. Every unchanged requirement retains
the starting SHA above as its implementation reference. Supported configurations
are source-inspected boundaries, not tested configurations. Four modes means
cosmo_cube, zoom_in, isolated_galaxy and isolated_cluster with the existing
coordinate/boundary policies; applicability depends on enabled physics.

| ID | Production owner / entry point | Actual source behavior and remaining obligation | Modes / MPI; assumptions and memory owner | Tests / execution / change commit |
|---|---|---|---|---|
| H01 | `src/hydro/hydro_core_solver.cpp`, reconstruction/riemann modules; `HydroAmrRuntime` | MUSCL, PPM, CTU and HLLC/HLLE present. Ghost-side CTU reductions remain counted; boundary convergence open. | Four gas modes; local/distributed paths. Solver batch scratch, no second hydro authority. | Existing Part 1/3 tests; NOT RUN — USER-RESERVED; unchanged |
| H02 | same; `src/amr/amr_ghost_fill.cpp` | Auxiliary dual energy and comoving CFL present; remote auxiliary lane may fall back to E−K. | Four gas modes; remote high-Mach qualification open. Canonical internal energy plus declared solver scratch. | Existing dual-energy/reflux tests; NOT RUN — USER-RESERVED; unchanged |
| C01 | `src/physics/thermochemistry.cpp`, `thermochemistry_operator.cpp`; `SourceRuntime::executeThermochemistry` | Hybrid light chemistry, cooling/data gates present; missing verified R18/photodissociation heat and production tables. | Four gas modes; explicit unqualified-rate and UVB compatibility opt-ins retained. Composition sidecar owns history. | Backend/rate tests; NOT RUN — USER-RESERVED; unchanged |
| C02 | same | Analytic/equilibrium/detailed references present; detailed inventory is limited, no case B. No scientific qualification. | Selected backend; bounded per-row chemistry solve. | Backend/rate tests; NOT RUN — USER-RESERVED; unchanged |
| C03 | `SourceRuntime::updateGlobalChemistryControl`, `ThermochemistryOperator::decideChemistryFidelity` | Heuristic global promotion bias and local hysteresis; private bias survives rollback. Timescale/stiffness and resource scheduling open. | Four chemistry modes; collective shadow reduction; sidecar plus private cached bias. | Adaptive/rollback fixtures; NOT RUN — USER-RESERVED; A5 then D |
| C04 | `ThermochemistryOperator::advanceRow` | Sampled chemical shadow; higher-order shadow hydro and matched-physics accuracy experiments open. | Bounded samples, no production mutation by shadow. | Adaptive fixtures; NOT RUN — USER-RESERVED; unchanged |
| I01 | `src/physics/effective_multiphase_ism.cpp`, `baryonic_source_stage.cpp` | Effective/resolved regimes and source transactions present; EOS edge/relaxation qualification open. | Four gas modes; composition tags and staged row scratch. | ISM/rollback tests; NOT RUN — USER-RESERVED; unchanged |
| I02 | `src/physics/ism_dust_sgs.cpp`, composition regrid | SGS reservoir/dissipation present; shear production and calibrated closure/refinement response open. | Opt-in research model; sidecar lane, no DMO lane. | ISM/dust/SGS tests; NOT RUN — USER-RESERVED; unchanged |
| S01 | `SourceRuntime` SF path; `src/physics/star_formation.cpp` | Births and remote gradients present; MPI identity/conservation/gradient comparisons open. | Four baryonic modes; conservative transfer and admitted population growth. | Existing gradient/birth tests; NOT RUN — USER-RESERVED; unchanged |
| S02 | `src/physics/stellar_population.cpp`, `supernova_feedback.cpp`; `SourceRuntime::executePart2StellarPopulation` | SSP/event records present; cumulative quadrature not structurally monotone. Associations rank-local, bucketed, omit velocity-dispersion energy. | Stellar modes; population sidecar and staged budgets. Production yields data-blocked. | Stellar population/remote sources tests; NOT RUN — USER-RESERVED; A4 |
| F01 | `MechanicalFeedbackDeposition`, `SourceRuntime::importRemoteFeedbackCandidates` | Remote deposition present; remote search precedes association centroid change. Association repair or explicit refusal required. | Stellar modes; local/remote candidates, staged deltas; MPI end-to-end open. | Remote sources tests; NOT RUN — USER-RESERVED; A4 |
| F02 | `src/physics/stellar_wind.cpp`, source runtime | Wind parser/interval budgets present; tables absent; level-local event/source timing open. | Table-selected source; immutable data, population event ordinals. | Wind/event tests; NOT RUN — USER-RESERVED; B3 |
| D01 | `src/physics/element_diffusion.cpp`, `SourceRuntime::executeMetalDiffusion`, composition exchange | Local diffusion only; caller already reserves delta. Full topology/temporary coexistence audit and distributed subcycle transport open. | Four gas modes; multiple ranks explicitly refused. Composition truth, phase scratch and face graph. | Remote sources/composition tests; NOT RUN — USER-RESERVED; A3 then C1 |
| A01 | `src/amr/amr_production_regrid.cpp`; `HydroAmrRuntime::runAutomaticAmrRegrid` | Automatic local/owned refine/derefine present; split-rank octets skipped, operations individually committed. Whole-hierarchy transaction/genealogy/migration open. | Four gas modes; conservative composition/SGS transfer; old/new admission. | Regrid tests; NOT RUN — USER-RESERVED; C2 |
| A02 | `advanceProductionAmrHydroRecursive`, `HydroAmrRuntime` | Local nested hydro recursion; MPI explicitly refused. Per-substep remote history/ghosts/flux collectives absent. | Single rank only for recursion; same-level/periodic topology in existing synchronized MPI path. | Recursive tests; NOT RUN — USER-RESERVED; B1 |
| A03 | `RefinementEvaluator`, production planner | Mass/jump/Jeans indicators; cooling length, chemistry interface, vorticity, force indicators incomplete. | Four gas modes; mandatory resolution must not silently downgrade. | Regrid tests; NOT RUN — USER-RESERVED; D |
| A04 | production planner, `amr_predictive_feedback_age_years` | Young-star age predictor; no general event forecast/error-versus-cost selection. | Research opt-in; governor admission, no changed source budgets. | Adaptive/deck tests; NOT RUN — USER-RESERVED; D |
| G01 | TreePM / `GravityRuntime` | Existing gravitational authority retained; fine-force correction conditional on evidence, not authorized as an automatic solver replacement. Zoom diagnostic not wired to report. | Four mode policies; valid force epochs/softening required. Gravity owns endpoint/cache memory. | TreePM–AMR reference fixture; NOT RUN — USER-RESERVED; B4/qualification |
| T01 | `TimeCoordinator::runRungZeroSegment`, source intervals/transactions | Optional full-state retry only at rung zero; private chemistry cache not restored; bounded comprehensive journal absent. | Rung-zero retry; full-state copy plus schedulers/forces coexist (complete budget audit open). | Config/source rollback fixtures; NOT RUN — USER-RESERVED; A5/B5 |
| T02 | coordinator, recursive hydro driver | Finest-cell global CFL; coarse source intervals; frozen gas gravity during recursion; particle/gas/source multirate coupling incomplete. | Local recursive hydro only; collisionless hierarchical KDK is separate existing authority to preserve. | Recursive tests; NOT RUN — USER-RESERVED; B2–B4 |
| M01 | `MemoryGovernor`, node pool; AMR/history/source owners | Node foundation present. Recursive rows allocate before admission; per-level history replaces canonical state during construction. Diffusion topology/coexistence and full rollback budgets need review. | Existing shared_equal/shared_pool policy; capacity accounting is not RSS/PSS. No new baryonic work in DMO. | Memory/refusal fixtures; NOT RUN — USER-RESERVED; A1–A3/B5 |
| M02 | scalar/OpenMP/MPI kernels, runtime balance planner | CPU paths/benchmark sources present; no new performance or node-fit evidence. | Scalar CPU authority, optional accelerators; per-worker and rank costs require owner measurements. | Benchmark sources; NOT RUN — USER-RESERVED; E |
| V01 | snapshot/restart/diagnostics; `src/io/snapshot_hdf5.cpp` | Science composition/star export present; import absent. Four-mode decks exist; MPI/new-physics/restart/migration qualification incomplete. | Existing HDF5 schemas; output only at legal barriers. | Four-mode/snapshot sources; NOT RUN — USER-RESERVED; E |
| V02 | adaptive chemistry/AMR runtime controllers | Research heuristics, no qualified coordinated error/resource optimizer. | Fixed reference/research identities must remain explicit; bounded diagnostic samples. | Adaptive tests; NOT RUN — USER-RESERVED; D/E |
| X01 | scope boundary | No new MHD, RT, CR, AGN, relativistic/nuclear/dark-sector engine or SUBFIND unbinding. | Deliberately absent; preserve existing guarded capabilities. | No new scope; unchanged |

Carried F1/F3/F4/F7: source repairs in Part 1, qualification still reserved.
F2/F5: Part 2 source repair; full-interval/birth/event regressions remain reserved.
F6: production regrid present, distributed recursion still absent. Additional
carried risks: cooling row mapping, EOS edges, finite relaxation plus cooling,
high-Mach cancellation, SSP monotonicity, remote dual energy, CTU boundary order,
composition donor fractions frozen over recursive substeps, science import,
zoom reporting, missing data, MPI tests and memory/RSS reconciliation.

## Checkpoints

Initial implementation checkpoint: A1/A2 AMR working-set admission and temporal
history publication. Milestone A remains PARTIAL until A3–A5 and their
production-path regressions are addressed. Milestones B–E remain OPEN; no
distributed time integration is claimed by local recursive-hydro source.

## A1/A2 source checkpoint — completed scope

**SOURCE REPAIRED; NOT RUN — USER-RESERVED.** This checkpoint closes the
specific recursive-driver admission and temporal-publication defects in source.
It does **not** close M01, A02, T01 or Milestone A as a whole.

Production integrations:

- `src/amr/amr_hydro_orchestrator.cpp::advanceProductionAmrHydroRecursive`:
  fixed-size level/count census before heap allocation; exact row capacities
  total O(active rows), outer vector owners, descriptors, level/ratio/interval
  tables and diagnostic counters admitted together as `amr.recursive.driver`.
  Capacity over-grants fail before hierarchy evolution. Descriptor capacity
  is reused at reflux rather than rebuilt.
- `hasProductionAmrHydroCoverage`: removes the O(Ncells) bitmap. Each validated
  contiguous range requires every row's patch index to equal its owner. Two
  different ranges cannot overlap without failing that equality; their summed
  lengths equal Ncells iff coverage is complete. Geometry and stable identity
  checks remain. This is an equivalent coverage proof, not a physics change.
- The reused synchronized sweep admits its separate local/remote descriptor
  copies before constructing either, including coexistence. Pending reflux
  candidates now accept and use the existing governor before allocation.
- `src/amr/amr_ghost_fill.cpp::captureHistoryStart`: one candidate contains new
  level records plus exact-capacity copies of all retained levels. No canonical
  clear, intermediate publication or allocating restoration path remains.
- `captureAmrTemporalBoundaryHistoryEndForLevel`: validates and fills a complete
  candidate before moving it into the canonical store. Failure after processing
  an earlier record cannot complete that record in canonical history.
- History candidates admit outer records, every nested cell array, geometry
  construction coexistence (start only) and an explicit 64-byte per-allocation
  allowance. The allowance is a model, not measured allocator occupancy;
  retained allocator uncertainty remains in the existing opaque process reserve.
  Actual nested/outer capacities determine the retained baseline.
- `AmrTemporalBoundaryHistoryStore::eraseLevel` compacts existing records with
  nonthrowing moves/destruction. Retirement constructs no record/cell copy.
  `AmrRetainedStateAdmission::commitShrink` computes the reduced baseline before
  canonical erase, then reconciles freed nested capacity with the same governor
  and attached node pool. Outer vector capacity remains charged while retained.
- `mergeFluxRegistersIntoPendingStore` already builds a candidate and publishes
  once. This behavior is retained, given an outer allocation allowance and
  covered by late-incompatibility and allocation-failure regressions.

Memory/transaction boundaries: driver preparation refusal/failure precedes gas,
patch, identity, history and register mutation. History operations individually
have strong exception safety. Existing state and all candidate levels coexist
under admission; no second SimulationState is introduced. All source-changing
callbacks publish with nonthrowing moves. This assumes the existing serial
workflow mutation boundary; concurrent mutation is unsupported. The small
returned diagnostic summary transfers to the caller as before. DMO dispatch
acquires no new baryonic owner or working set. The existing process/node policy
is used; no independent budget or communicator is created.

**Limit:** a refusal in a later operation after an earlier hydro substep committed
is still a whole-step rollback problem (B5/T01). This checkpoint does not make
an entire recursive advance transactional. Distributed recursive hydro remains
explicitly unavailable, and no fresh MPI numerical evidence is asserted.

### Tests written and registered

All tests below are **NOT RUN — USER-RESERVED**, including existing numerical
assertions. No assertion or tolerance was weakened.

| CTest name | Source additions | Target |
|---|---|---|
| `unit_amr_history_atomicity` | New fixture; admission refusal with an existing parent interval, invalid last cell/geometry after earlier records, duplicate descriptors, enumeration of C++ allocation failures for nested start and end capture, allocation-free retirement and baseline reconciliation | Bit-identical history fields, unchanged canonical storage on failure; no leaked reservation |
| `unit_amr_recursive_subcycling` | Single-/three-level driver memory refusal and selected row-allocation failure; overlap/missing coverage rejection | Exact gas-lane bytes/identity generations retained, empty histories/registers, all leases released before return; existing numerical/conservation/clock tests retained |
| `unit_amr_pending_flux_registers` | Candidate partially built before incompatible existing metadata; enumeration of allocation failures before a successful merge | Canonical register storage/counters/flux fields unchanged on failure and governor baseline/leases preserved |
| `unit_amr_temporal_ghost_fill` | Existing endpoint/interpolation/identity/restart-admission contracts retained | Numerical tolerances and historical API calls unchanged; NOT RUN |

`tests/support/scoped_allocation_failure.hpp` replaces C++ allocation functions
only in the three standalone fault-test executables that include it. It is
thread-local, one-shot, covers ordinary/aligned C++ new and is disabled outside
an explicit scope. It does not inject failures into opaque library malloc calls.
`unit_amr_history_atomicity` is registered in the existing `chui_part3` unit list;
the other three targets retain their existing registration. No MPI test was
written in this checkpoint, and E's 2-/4-rank matrix remains OPEN.

Owner-run commands, provided only as a standalone validation selection:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -R '^unit_amr_(history_atomicity|recursive_subcycling|pending_flux_registers|temporal_ghost_fill)$'
```

Then retain the normal complete CPU suite and relevant MPI/HDF5/restart presets
as owner gates. New convergence, performance and science acceptance are not
implied by the fault fixtures. No benchmark kernel was added: these are lifetime
and transaction changes, and existing recursive-hydro/kernel benchmarks remain
the future measurement surface. No speedup or memory-fit claim is made.

### Compatibility

No numerical equation, flux law, integration order, coordinate conversion,
configuration key/default, HDF5 dataset, restart payload or provenance schema
changed. Successful history endpoints retain the same mathematical values.
Insufficient-memory and invalid-input behavior becomes earlier/atomic.
End capture now requires candidate headroom; tightly budgeted configurations may
refuse an operation that previously mutated history in place unsafely.

Source API additions are trailing defaulted governor arguments on both history
retirement functions and `applyCompletePendingFluxRegistersToSimulationState`.
Existing calls compile without changes; callers managing a governor should pass
it for retirement/reflux. Rebuild dependents (binary ABI is not promised).
The store's `assign` is explicitly noexcept, and `eraseLevel` is a new primitive.

### Continuation — exact next work

1. **A3 first**, before B. In `SourceRuntime::executeMetalDiffusion` around its
   `sources.element_diffusion.workspace` reservation, preserve the fact that
   delta is already admitted. The primary estimator excludes delta and small
   lanes/before/scale vectors; consolidate with checked size arithmetic and
   include caller-owned retained capacity/growth. The topology builder
   `src/workflows/internal/metal_diffusion_topology.cpp::buildMetalDiffusionTopology`
   also constructs faces, strain, boundary descriptors, directed refs, plane
   map and unique-face set before caller admission. Its reported scratch is
   incomplete (map/set not included). Repair that full lifecycle, including
   old/new face/rho-kappa buffers and retained scalar-diffusion workspace, before
   claiming no known major ungoverned diffusion scratch. MPI remains refused;
   no communication buffers exist in this local path to double-count.
2. **A4:** `SourceRuntime::executePart2StellarPopulation` imports remote feedback
   candidates before `mergeAssociationSources`. Current grouping is local to
   source owner, bucket-based and loses velocity-dispersion kinetic energy.
   The authorized safe option is to reject `ssp_clustered_associations=true`
   in typed config and direct runtime construction, and keep the unsupported
   merger unavailable, with explicit refusal tests/documentation. Do not merely
   reorder the existing approximate merge and call associations conformant.
3. **A5:** `SourceRuntime::updateGlobalChemistryControl`, its
   `m_chemistry_control_restored` flag/`m_last_shadow_report`, and
   `ThermochemistryOperator::m_promotion_bias` must join the global rollback
   contract or become derivable from the restored sidecar before decisions.
   Restore restart bias before the first fidelity decision. Empty ranks must
   contribute zero current shadow samples, not stale prior-stage diagnostics.
   Test an injected later rejection after bias publication and compare private
   controller state plus authoritative sidecar/tags/counters before/after retry.
4. Complete A's remaining source/test/documentation checkpoint, then undertake
   B1 remote per-substep exchange/history/reflux under the existing coordinator.
   B2 CFL relaxation is unsafe until each affected operator has its required
   cadence. Preserve existing collisionless KDK and TreePM authorities. B3/B4/B5
   remain integrated requirements, not independent optional toggles.
5. C/D/E retain every obligation in the map above. In particular, do not mark
   any MPI requirement closed from local fault tests or historical test reports.

Commands executed in this session were source/document reads, archive extraction,
source-editing utilities, Git inspection/branch/commit/push operations. No project
code or executable validation was run. The design-bundle estimator and data
converters were never invoked. CI is not dispatched; this feature branch does not
match the repository workflow's push-to-main trigger, and commits carry skip-ci.

## Committed handoff

Implementation checkpoint SHA: **`ecb00914100b68835dc46fe00f6e685a023edfd3`**.
This is the change reference for the A1/A2 subdefects of M01 and the history
transaction portion of A02/T01. All other source obligations retain the starting
SHA as reference; they have not been relabeled complete. The subsequent
handoff-only commit records this SHA and this file's exact touched-file inventory.
Both commits are on the dedicated feature branch. No PR, main merge, rebase or
force-push is authorized or performed.

| File | Purpose |
|---|---|
| `src/amr/amr_hydro_orchestrator.cpp` | Recursive/sweep descriptor and row admission; bitmap-free coverage; governed reflux candidate preparation and governor propagation |
| `src/amr/amr_ghost_fill.cpp` | Atomic start/end candidates and no-copy retirement with baseline reconciliation |
| `src/amr/internal/amr_retained_state_admission.hpp` | Checked in-place shrink transaction |
| `include/cosmosim/core/simulation_state.hpp` | Nonthrowing temporal-store publication and level compaction primitive |
| `include/cosmosim/amr/amr_ghost_fill.hpp` | Documented strong exception contract and defaulted retirement governor |
| `include/cosmosim/amr/amr_hydro_orchestrator.hpp` | Defaulted pending-reflux governor |
| `tests/unit/test_amr_history_atomicity.cpp` | New multi-level history failure/admission/retirement fixtures |
| `tests/unit/test_amr_recursive_subcycling.cpp` | Driver refusal/allocation failure and coverage regressions; existing numerical tests unchanged |
| `tests/unit/test_amr_pending_flux_registers.cpp` | Late merge refusal and allocation-failure publication fixtures |
| `tests/support/scoped_allocation_failure.hpp` | Test-only thread-local one-shot C++ allocation interception |
| `CMakeLists.txt` | Registers `unit_amr_history_atomicity`; other targets already registered |
| `CURRENT_STATUS.md` | Current source-only checkpoint and incomplete Milestone A status |
| `docs/repair_open_issues.md` | Prioritized remaining A3–A5 defects |
| `docs/repair_state_recap.md` | Current repair state and ledger link |
| `docs/validation_plan.md` | Owner-reserved test selection and reproducibility distinctions |
| `docs/architecture/amr_temporal_ghost_interpolation.md` | Ownership/lifetime/API migration and transaction limits |
| `docs/repair/chui_conformance_closure_20261104.md` | Complete requirement map, source closure boundary, commit identity and continuation |

Git inspection confirmed the workflow's automatic push trigger is restricted to
`main`. Only the feature branch is pushed; commit messages also include
`[skip ci]`. Git hooks are disabled per commit/push command to honor the absolute
prohibition on executing repository code. This is session-local command
configuration, not a repository policy change. No workflow was dispatched.

The required CPU and MPI/HDF5 executable evidence remains blocked by the
**user-reserved execution contract**, not by a claim of missing dependencies.
There is no compile, link, runtime, numerical-accuracy, performance or science
qualification result for this checkpoint. The next source checkpoint is A3's
complete diffusion memory lifecycle; remaining budget is not being spread over
partial B/C/D/E implementations.
