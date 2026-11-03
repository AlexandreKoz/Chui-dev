# CHUI Part 1 — foundation, hydro suite, memory admission, audit F1/F3/F4/F7

Date: 2026-11-03. Base revision: `4c4a807` (the audit's pinned revision).
Initial source pass: commit `7157116`. Corrective pass: `6ca40e1`.
Stabilization gate: the commit that adds this revision of the record.

**Evidence status: SOURCE ONLY. No configure, build, compiler, test,
benchmark, sanitizer or simulation command was executed (explicit owner
instruction). Every check below is NOT RUN — USER-RESERVED.** Compilation
errors are possible and must be fixed by the first human build.

Design inputs read: `chui_final_implementation_route_consolidation_v10.md`,
`chui_consolidated_physics_performance_plan_v01.md`,
`chui_full_physics_memory_architecture_specification_v02.md`,
`CHUI_Baryonic_Physics_AMR_Audit_2026-10-30.md` (handoff zip, not committed).

## 1. Implemented capabilities (source)

| Area | Where | Notes |
|---|---|---|
| F1 Hancock predictor | `src/hydro/hydro_reconstruction.cpp` | forward signs, cell-centre evaluation along the face normal, EOS-closed face states, ordered positivity fallback; `dt/(a dx)` via `hydroTransportCourantFactors` in fixed-grid and AMR runtimes |
| PPM | same | CW84 interfaces + monotonicity + per-cell flattening, single Hancock step; stencil completion through the ghost shell |
| CTU | `src/hydro/hydro_core_solver.cpp` | simplified 3D Colella, 2 solves/face, exact donor-cell corner transport for first-order states |
| F3 Riemann | `src/hydro/hydro_riemann.cpp` | lane-generic HLL, Einfeldt HLLE bounds, HLLC/HLLE selectable |
| Dual energy | solver + `HydroConservedStateSoa` aux lane + both runtimes | see section 3 |
| F4 refresh | `gas_thermodynamic_refresh.{hpp,cpp}`, `source_runtime.cpp`, `stellar_feedback.cpp`, `time_coordinator.cpp` | feedback touch log, SF/BH whole-rank refresh, stale-cache guard, runtime-wide temperature convention |
| F7 ledger | `effective_multiphase_ism.cpp` | `rho_c du V_c`; sidecar `energy_ledger_version=2` |
| Thermal authority | `amr_hydro_geometry.cpp`, `amr_hydro_orchestrator.cpp` | AMR loads canonical `e`, not the pressure cache |
| Periodic boundaries | `hydro_amr_runtime.cpp`, `amr_hydro_geometry.cpp`, `amr_hydro_orchestrator.cpp` | fixed-grid double flux fixed; AMR in-patch wrap; cross-patch periodic topology refused |
| Stencils | `hydro_boundary_conditions.cpp` | `completeReconstructionStencils` for every boundary kind |
| Memory | `node_memory_admission.{hpp,cpp}`, `memory_governor.{hpp,cpp}`, `parallel/node_memory_pool.{hpp,cpp}`, `reference_workflow.cpp` | `shared_equal` (deterministic equal ceilings) and `shared_pool` (node-shared counter); owner telemetry |
| Transactions | `stage_transaction.{hpp,cpp}`, `workflows/internal/mpi_stage_collective_decision.hpp` | used around the fixed-grid hydro commit |
| Config / restart | `config.{hpp,cpp}`, `restart_checkpoint.cpp` | 9 keys (+ `shared_pool` value); exact-default migration projection |

## 2. Corrective-pass changes (this revision)

1. **Dual energy redesigned** (the initial pass reset `E` from a step-local
   auxiliary computed from `E - K`, which neither avoided the cancellation it
   was meant to fix nor conserved energy): auxiliary lane initialized from the
   canonical `e`, Bryan-type selection, `E` never overwritten by the solver,
   full ledger of the thermal correction, unknown-lane and floor counters, the
   reference `off` mode unchanged.
2. **AMR thermal authority** repaired (pressure cache was used as authority).
3. **AMR periodic boundaries**: confirmed zero-flux box faces; in-patch wrap +
   once-per-face evaluation; refusal for cross-patch periodic topology.
4. **Stencils** completed through the ghost shell for open/reflective/AMR
   ghosts (PLM full order at patch boundaries except the ghost side; PPM full
   order except the last two layers).
5. **Temperature cache**: one convention (`T = (gamma-1) mu m_p e / k_B`,
   mu = 0.59) in the fixed-grid store, AMR scatter and refresh (AMR wrote `T = e`;
   fixed grid never updated `T`).
6. **F4 completeness**: star formation and BH accretion also trigger a refresh
   (effective-EOS signal speed depends on density).
7. **Node-shared pool** (`shared_pool`) attached to the existing governor.
8. Documentation now states the equations actually implemented
   (`docs/hydro_numerics.md`).

## 2b. Stabilization-gate changes

1. **Hydro CFL rule**: MUSCL-Hancock (default) now uses the dimension-summed
   condition `dt = 0.4 / sum_d (|v_d|+c_s)/(L dx_d)` (stable iff
   `sum nu <= 1` for the unsplit update; the old per-axis rule allowed 1.2);
   CTU keeps the per-axis rule. Same function for the TimeCoordinator and the
   accepted-step guard; inactive (domain-spanning single-cell) axes excluded.
2. **AMR thermodynamic caches**: scatter, reflux (local + remote), refine and
   derefine writes produce closure-consistent `P`, `c_s` and `T = f e`
   (`amr::AmrThermodynamicCachePolicy`); whole-rank post-AMR closure pass
   removed; both loaders use canonical `e` (pressure cache only as legacy
   fallback); star formation refresh is touched-row based; `mu` parameterized.
3. **Shared pool**: guarantee only after world-agreed compile-time + runtime
   checks (address-free allowlist, `MPI_WIN_UNIFIED`, alignment); otherwise
   deterministic `shared_equal` fallback, reported.
4. **Dual energy at AMR**: reflux keeps pre-reflux `e` where `E - K < eta E`
   (ledgered); ghost auxiliary energy from local source cells; remote sources
   counted as degraded.

## 3. Policies as implemented

* **Dual energy:** `docs/hydro_numerics.md` section 7. Conservative state owns
  rho, rho u, E, rho Z; thermal state = `E - K` if `>= eta E`, else auxiliary,
  else floor; correction `sum (rho e_sel - (E - K)) V` reported per step.
* **Memory:** `docs/memory_governance.md`. `shared_equal`: rank ceiling
  `floor((node_limit - reserve)/ranks)`, deterministic. `shared_pool`: node sum
  bounded by a lock-free shared counter, non-deterministic under contention.
  Accounting is not RSS.
* **StageTransaction:** the fixed-grid hydro writes SimulationState only in
  `commit()`. Staged data are the solver-owned conserved copy; the effective-ISM
  sidecar ledger is persisted only after commit; profiler events are
  diagnostics. A rejection leaves gas truth, identities, thermodynamic caches
  and ledgers untouched. A failure inside the commit callable is reported as
  `kCommitFailed` and rethrown (invariant violation). The AMR path is NOT
  transactional yet (patches scatter as they advance); that is Part 3 scope.

## 4. Audit closure map

| ID | Status | Tests (NOT RUN — USER-RESERVED) |
|---|---|---|
| F1 | repaired in source | `unit_hydro_muscl_hancock_predictor` |
| F3 | repaired in source | `unit_hydro_riemann_passive_scalars` |
| F4 | repaired in source (feedback, SF, BH; stale guard) | `unit_gas_thermodynamic_refresh` |
| F7 | repaired in source (isolated + comoving conventions) | `unit_effective_ism_energy_ledger` |
| F2 | open — Part 2 | — |
| F5 | open — Part 2 | — |
| F6 | open — Part 3 | — |

Additional confirmed defects repaired: fixed-grid periodic double flux;
AMR periodic zero-flux box faces; AMR pressure-cache thermal authority; AMR
temperature cache in code units. Supplementary audit risks still open:
fixed-grid cooling row permutation, SF gradients at patch boundaries, EOS
table edges, relaxation + cooling partition.

## 5. Remaining limitations

* PPM/CTU at AMR (and non-periodic fixed-grid) patch boundaries: explicit,
  counted reductions (table in `docs/hydro_numerics.md` section 3; CTU ghost
  side normal-only). Multi-layer ghost shells and corner ghosts are Part 3.
* **MANDATORY PART 3 DELIVERABLE — periodic cross-patch AMR topology.**
  Periodic AMR across patches is refused; decomposed periodic `cosmo_cube` gas
  whose root patches do not span the box stops with an explicit error.
  `cosmo_cube` is a first-class mode: Part 3 must implement wrapped neighbour
  detection, wrapped ghost sampling, wrapped coarse/fine flux registers and
  wrapped distributed exchange across MPI ranks, then lift the refusal.
* Dual energy: weak shocks below `eta` not heated by the auxiliary branch;
  remote AMR ghost sources fall back to `E - K` (counted); reflux keeps the
  pre-reflux `e` in cold cells (ledgered non-conservative correction).
* Gravity half-kick scaling of effective-EOS caches deviates by O(H dt) for
  one step (re-derived at the next hydro store).
* BH accretion triggers a whole-rank refresh (no touched-cell list).
* Not bit-reproducible with pre-Part-1 binaries (intended numerical fixes).

## 6. Migrations

* Config: `numerics.hydro_reconstruction`, `hydro_slope_limiter`,
  `hydro_integration_scheme`, `hydro_riemann_solver`, `hydro_dual_energy`,
  `hydro_dual_energy_eta`, `parallel.node_memory_admission`
  (`off|shared_equal|shared_pool`), `node_memory_budget_bytes`,
  `node_memory_reserve_bytes`. Defaults keep the historical method and memory
  semantics; normalized text/hash change for every configuration.
* Restart: no HDF5 schema change. Pre-Part-1 checkpoints resume only with all
  nine keys at their defaults (methods never silently activate on resume).
* Canonical state: unchanged layout. New rank-local, non-persistent thermal
  epochs on `SimulationState`. The auxiliary dual-energy lane is solver
  scratch rebuilt from `e`.
* ISM sidecar: `energy_ledger_version=2`; legacy density sums discarded.

## 7. Test inventory (all NOT RUN — USER-RESERVED)

`unit_hydro_muscl_hancock_predictor`, `unit_hydro_riemann_passive_scalars`,
`unit_hydro_ppm`, `unit_hydro_ctu`, `unit_hydro_dual_energy` (rewritten:
cold contact beyond E precision, warm-flow no-op, shock heating, contracts),
`unit_gas_thermodynamic_refresh` (+ temperature convention),
`unit_effective_ism_energy_ledger` (+ comoving vs isolated),
`unit_node_memory_admission` (+ shared pool), `unit_stage_transaction`,
`unit_config_hydro_numerics` (+ `shared_pool`), `unit_hydro_cfl_stability_rule`
(new: amplification-factor proof, anisotropic, comoving, guard),
`unit_amr_reflux_thermodynamics` (new: closure/temperature at reflux, legacy
path, dual-energy reflux fallback), `integration_amr_production_hydro_integration`
(+ refinement caches), `unit_node_memory_admission` (+ serial pool),
`unit_amr_periodic_boundary`
(new: wrapped ghosts, boundary-crossing donor-cell step with metals and
velocity, conservation, ghost-shell stencils), modified
`unit_effective_multiphase_ism`; benchmark source `bench_hydro_method_matrix`.

## 8. Human validation checklist (documentation only; not executed)

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

Then the MPI preset of `docs/build_instructions.md` for `shared_pool`
(several ranks per node) and restart migration, and the Part 1 campaign items
in `docs/validation_plan.md`.

## 9. Interfaces for Parts 2 and 3

* Part 2: conserved lanes are extended in `HydroConservedState`; source stages
  mutating `rho`, `e` or composition call `notifyGasThermodynamicMutation()` +
  `refreshGasThermodynamics()`; unfinished chemistry throws
  `StageModelDowngradeRequest` inside `StageTransaction::prepare`; feedback
  reports mutated cells through `StellarFeedbackThermodynamicTouchSink`;
  source ledgers integrate with `HydroSourceContext::cell_volume_comoving`;
  thermal sources enter the dual-energy auxiliary update through `S_E`.
* Part 3 (mandatory first item): periodic cross-patch AMR topology.
* Part 3: remote auxiliary-energy ghost exchange, multi-layer AMR ghost shells (`requiredStencilHalfWidth()` layers,
  corner ghosts for CTU), wrapped periodic patch topology (ghost sampling,
  flux registers, distributed exchange), transactional AMR sweep, dimension-
  scheduler-coordinated memory borrowing. (Dimension-summed CFL is now done.)
