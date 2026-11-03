# CHUI Part 1 — foundation, hydro suite, memory admission, audit F1/F3/F4/F7

Date: 2026-11-03. Base revision: `4c4a807` (the audit's pinned revision).
Mode: feature implementation (Part 1 of 3).

**Evidence status: SOURCE ONLY. No configure, build, test, benchmark or
simulation command was executed (explicit owner instruction). Every check
below is NOT RUN — USER-RESERVED.** Compilation errors are possible and must be
fixed by the first human build.

Design documents read: `chui_final_implementation_route_consolidation_v10.md`,
`chui_consolidated_physics_performance_plan_v01.md`,
`chui_full_physics_memory_architecture_specification_v02.md`,
`CHUI_Baryonic_Physics_AMR_Audit_2026-10-30.md` (supplied in the handoff zip,
not committed).

## 1. What changed (real code paths)

| Area | Files | Interfaces |
|---|---|---|
| MUSCL-Hancock predictor (F1) | `src/hydro/hydro_reconstruction.cpp`, `include/cosmosim/hydro/hydro_reconstruction.hpp` | `MusclHancockReconstruction` (corrected), `hydroTransportCourantFactors()`, `HydroFaceStencilStates`, `reconstructFaceFromStencil()`, `stencilDegradationCount()` |
| PPM | same | `PpmReconstruction`, `ppmCellEdges()`, `HydroReconstructionScheme` |
| Riemann (F3) | `src/hydro/hydro_riemann.cpp`, `hydro_riemann.hpp` | lane-generic HLL, Einfeldt HLLE bounds, `HydroRiemannSolverKind`, `makeHydroRiemannSolver()` |
| CTU + dual energy | `src/hydro/hydro_core_solver.cpp`, `hydro_core_solver.hpp` | `HydroIntegrationScheme`, `HydroDualEnergyPolicy`, `HydroNumericsPolicy`, `HydroCoreSolver(gamma, numerics)`, `requiredWorkspaceFeatures()`, new profile/conservation fields, `HydroSourceContext::cell_volume_comoving`, `HydroFace::owner_minus2_cell/neighbor_plus2_cell` |
| Geometry stencils | `src/hydro/hydro_cartesian_patch.cpp`, `src/amr/amr_hydro_geometry.cpp`, `src/hydro/hydro_boundary_conditions.cpp` | `completePeriodicReconstructionStencils()` |
| Fixed-grid runtime | `src/workflows/hydro_amr_runtime.cpp` | per-step reconstruction with `dt/(a dx)`, configurable reconstruction/Riemann/numerics, periodic face de-duplication, stage-transaction commit, telemetry |
| AMR runtime | `src/amr/amr_hydro_orchestrator.cpp`, `amr_hydro_orchestrator.hpp` | `ProductionAmrHydroOptions::reconstruction_scheme/slope_limiter`, `dt/(a dx)`, diagnostics counters |
| Thermodynamic refresh (F4) | `include/cosmosim/hydro/gas_thermodynamic_refresh.hpp`, `src/hydro/gas_thermodynamic_refresh.cpp`, `src/workflows/source_runtime.cpp`, `src/physics/stellar_feedback.cpp`, `include/cosmosim/physics/stellar_feedback.hpp`, `src/workflows/time_coordinator.cpp`, `SimulationState` | `refreshGasThermodynamics()`, `StellarFeedbackThermodynamicTouchSink`, `notifyGasThermodynamicMutation()/markGasThermodynamicsRefreshed()/gasThermodynamicCachesFresh()` |
| ISM ledger (F7) | `src/physics/effective_multiphase_ism.cpp`, header, `hydro_amr_runtime.cpp` sidecar | `kEffectiveIsmEnergyLedgerVersion = 2` |
| Memory | `include/cosmosim/core/node_memory_admission.hpp`, `src/core/node_memory_admission.cpp`, `memory_governor.{hpp,cpp}`, `reference_workflow.cpp` | `computeNodeMemoryAdmission()`, `discoverNodeMemory()`, parsers, `MemoryGovernor::ownerTelemetry()` |
| Transactions | `include/cosmosim/core/stage_transaction.hpp`, `src/core/stage_transaction.cpp`, `src/workflows/internal/mpi_stage_collective_decision.hpp` | `StageTransaction`, `StageCollectiveDecision`, `StageNumericalRetryRequest`, `StageModelDowngradeRequest`, `StageTransactionRejected` |
| Config | `include/cosmosim/core/config.hpp`, `src/core/config.cpp` | 9 typed keys (docs/configuration.md) |
| Restart migration | `src/io/restart_checkpoint.cpp` | exact-default projection extended to the 9 keys |

## 2. Numerical contracts

See `docs/hydro_numerics.md` (equations). Summary: forward Hancock
predictor with physical signs evaluated at the cell centre along the face
normal; `lambda = dt/(a dx)` at the flux epoch shared with the update;
PPM = CW84 interface/monotonicity/flattening + one Hancock step on the mean
parabola gradient; CTU = simplified 3D Colella (2 solves/face, Courant <= 1/2
in 3D); all conserved lanes through every Riemann branch; step-local dual
energy with `E` synchronization and measured `dual_energy_sync_delta`; ISM
ledger in `rho delta_u V` code energy.

## 3. Memory and concurrency contract

See `docs/memory_governance.md`. Node rule `share = floor((node_limit -
reserve)/ranks_on_node)`; no collectives; serial/no-MPI uses one rank.
Governed: hydro active batch (incl. CTU/dual lanes), feedback touch log,
existing governed owners. Ungoverned: FFTW/MPI/HDF5 internals (explicit
reserves only), small per-step method objects, diagnostic copies. Solver
counters are per-instance (serial per patch); OpenMP use requires
thread-owned instances.

## 4. Audit closure map

| ID | Status | Source | Tests (NOT RUN) |
|---|---|---|---|
| F1 | repaired in source | `hydro_reconstruction.cpp` (`hancockHalfStepIncrement`, `finalizeSide`), `hydroTransportCourantFactors` use in `hydro_amr_runtime.cpp` and `amr_hydro_orchestrator.cpp` | `unit_hydro_muscl_hancock_predictor` (2.4 not 2.6, axes, reversed velocity/normal, a=2, pressure sign, Z) |
| F3 | repaired in source | `hydro_riemann.cpp` `hllFlux` | `unit_hydro_riemann_passive_scalars` (metal 0.2, forced fallback, orientations, contact) |
| F4 | repaired in source | `source_runtime.cpp` touch log + refresh, `stellar_feedback.cpp` sink, `time_coordinator.cpp` stale guard | `unit_gas_thermodynamic_refresh` (u x100 -> c x10, CFL /10, branch transition) |
| F7 | repaired in source | `effective_multiphase_ism.cpp`, solver volume injection, sidecar v2 | `unit_effective_ism_energy_ledger` (parent vs 8 children, unequal volumes, solver consistency) |
| F2 | open — Part 2 | — | — |
| F5 | open — Part 2 | — | — |
| F6 | open — Part 3 | — | — |

Supplementary audit risks:
1. Fixed-grid cooling row permutation (provider indexes canonical rows by
   solver row) — inspected only partially; NOT changed. Remains open.
2. Adaptive SF gradients at patch boundaries — not inspected; open (Part 2).
3. Effective EOS table boundary/derivative behaviour — not changed; open.
4. ISM relaxation plus ordinary cooling partition — not changed; open (Part 2).
5. High-Mach thermal accuracy — dual-energy option added; quantitative decision
   pending human validation item 3.
New finding repaired: periodic fixed-grid boxes applied every boundary flux
twice (both ghost copies active) — fixed, regression in `unit_hydro_ctu`.
Risk hypothesis (not reproduced, NOT repaired): AMR physical periodic ghosts use
`source_real_cell = owner` (net-zero flux at the box boundary) — needs a
dedicated reproduction; Part 3 AMR scope.

## 5. Migrations

* Config: 9 new keys; defaults preserve historical method/memory semantics.
  Normalized text and its hash change for every configuration.
* Restart: no HDF5 schema bump; pre-Part-1 checkpoints resume through the
  exact-default projection (`docs/restart_checkpointing.md`).
* ISM sidecar payload: `energy_ledger_version=2`; legacy cumulative density sums
  discarded with `legacy_density_ledger_discarded=1`.
* Public API: additive (new fields at struct ends, defaulted constructor
  argument). `MemoryReservation` gained a private owner slot (ABI only).

## 6. Test inventory (all NOT RUN — USER-RESERVED)

| ctest name | invariant |
|---|---|
| `unit_hydro_muscl_hancock_predictor` | F1 signed advection/pressure/Z, 1/a, orientation, closure |
| `unit_hydro_riemann_passive_scalars` | F3 passive lanes in all branches, contact |
| `unit_hydro_ppm` | CW84 exactness on quadratics, monotonicity, flattening, degradation, single predictor |
| `unit_hydro_ctu` | exact donor-cell CTU (a=1, a=2), 2 solves/face, free stream 3D, periodic face once, stencil wrap |
| `unit_hydro_dual_energy` | switch, sync accounting, high-Mach contact and boosted wave |
| `unit_gas_thermodynamic_refresh` | F4 refresh, epochs, CFL tightening, branch transition |
| `unit_effective_ism_energy_ledger` | F7 volume weighting |
| `unit_node_memory_admission` | parsers, 8/16/48/128 GiB synthetic shares, overflow, owner telemetry |
| `unit_stage_transaction` | commit/reject lifecycle, single vote, peer rejection |
| `unit_config_hydro_numerics` | defaults, round trip, fail-closed validation |
| existing `unit_effective_multiphase_ism` | updated: supplies `cell_volume_comoving` |
| `bench_hydro_method_matrix` | benchmark source (not a test) |

## 7. Human validation checklist (documentation only; not executed)

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
ctest --preset test-cpu-debug -R "unit_(hydro|gas_thermodynamic|effective_ism|node_memory|stage_transaction|config_hydro)" --output-on-failure
```

Then the MPI/HDF5 presets listed in `docs/build_instructions.md` for restart
migration and multi-rank admission, and the campaign items in
`docs/validation_plan.md` (Part 1 section).

## 8. Part 2 handoff (facts, not promises)

* Conserved lanes: `HydroConservedState` {mass, momentum x/y/z, total energy,
  metal mass}; primitives carry `metallicity_mass_fraction`. New species/dust
  lanes must be added to the struct, its arithmetic, SoA, flux-register and MPI
  correction records; HLL/HLLC then transport them automatically.
* Thermodynamic refresh: chemistry/cooling/feedback stages that change
  `density_code`, `internal_energy_code` or composition outside hydro must call
  `notifyGasThermodynamicMutation()` + `refreshGasThermodynamics()` (closure via
  `HydroThermodynamicClosure`).
* Transactions: chemistry intervals that cannot finish throw
  `StageModelDowngradeRequest` (or `StageNumericalRetryRequest`) inside
  `StageTransaction::prepare`; commit only after `validate()`.
* Feedback deposition reports every mutated cell through
  `StellarFeedbackThermodynamicTouchSink` before mutation.
* Source ledgers receive `HydroSourceContext::cell_volume_comoving` and must
  integrate energies over it.
* Memory: request batch lanes through `HydroActiveBatchPolicy`/governor
  reservations; per-owner labels appear in `ownerTelemetry()`.

## 9. Blockers / uncertainty

* No compilation: syntax/type errors may exist (largest new units:
  `hydro_core_solver.cpp` batched advance, `hydro_reconstruction.cpp`).
* PPM/CTU at AMR patch boundaries are reduced-order (counted), because AMR
  ghost shells are one cell wide. Full-order shells are Part 3.
* Dual energy is step-local; thermal energy below `ulp(E)` is unrecoverable
  without a persistent entropy lane.
* The CFL coefficient remains 0.4; CTU in 3D requires <= 0.5 (satisfied), but
  user-adjustable CFL does not exist and was not added.
* Node admission assumes job-shared cgroups and does not observe other jobs on
  shared nodes; use an explicit node budget there.
* Not bit-reproducible with pre-Part-1 binaries (intended bug fixes).
