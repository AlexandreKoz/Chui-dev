# CHUÍ MHD — implementation ledger, pass 04 (corrective campaign)

**Status vocabulary:** `OPEN`, `IN_PROGRESS`, `IMPLEMENTED_UNVALIDATED`, `PARTIAL`,
`FAIL_CLOSED`, `NOT_APPLICABLE`. No row is `PASSED` or `CERTIFIED`.
**Source-only campaign:** no configure, build, CTest, simulation or benchmark was
run, nothing was committed, pushed, reset or stashed. Every evidence cell is
`NOT_RUN`; compile viability of every edited file is unverified. A guard is
always reported as `FAIL_CLOSED`, never as the feature it protects.

Inputs: REF `chui_mhd_architecture_specification_v01.md`, PERF
`chui_mhd_memory_performance_specification_v01.md`, AMEND
`chui_mhd_dev_branch_compatibility_amendment_v01.md`, pass-02 ledger, audit
`chui_mhd_specification_compliance_audit_20261127.md` (added to the tree from the
pass-04 bundle, together with `agent_prompts/chui_mhd_corrective_pass04_*`).

## 0. Worktree identity

| Field | Value |
|---|---|
| Repository | `/home/user/Chui-dev` (existing local clone) |
| Branch | `claude/inspiring-einstein-33hhmz` |
| HEAD at start and end | `2795d11f12040b824fcdecd3fda749a70c9eeb0c` (no commit made) |
| Status at start | clean |
| Status at end | 46 modified tracked files, 5 new files (this ledger, the audit, two bundle prompts, one new test) — see §6 |
| Static checks run | `git diff --check` (clean); textual inspection only |

## 1. Audit findings (M03-01..15, M04-01..05)

| ID | Status | Implementation paths / symbols | Production caller | Memory / state authority | Regression test source | Remaining gap | Evidence |
|---|---|---|---|---|---|---|---|
| M03-01 global epoch | IMPLEMENTED_UNVALIDATED | `MhdRuntimeBridge::buildDirectory` (canonical sorted integer tuples, duplicate-id check, digest), `agreeHierarchyValue` (allreduce max(v), max(~v), max(disagree)), `m_hierarchy_epoch = digest (+) m_hierarchy_generation`; `MhdLevelStepper::directoryDigest` computed in `buildTopology`; `MhdLevelStepperConfig::hierarchy_epoch` doc contract | `HydroAmrRuntimeImpl::execute` / `initializeMagneticState` → bridge `bind` → stepper wire header | Bridge owns the epoch; `topologyFingerprint` is a rank-local cache key only; stage/attempt sequence stays in the phase word | `test_mhd_hierarchy_runtime.cpp::testStaleEpochRejectedCollectively`, `testTwoRankEquivalence` (now with a shared epoch) | Real `mpiexec` 1/2/4-rank production-bridge test with unequal ownership not written/registered | NOT_RUN |
| M03-02 composition transport | FAIL_CLOSED | `mhdProductionDispatchRefusal` refuses `physics::part2ThermochemistrySelected` (the exact activation condition of `ensureGasCompositionInitialized`); bridge constructor reuses the predicate; collective runtime refusal of any active `gas_composition` sidecar in the magnetized branch | capability admission, `makeHydroAmrRuntime`, `MhdRuntimeBridge` ctor, `HydroAmrRuntimeImpl::execute` | composition owner unchanged; MHD lanes: rho, m, E, legacy metal mass | `test_mhd_amr_mpi_restarts.cpp::testConfigAndCapabilityAdmission` (composition case) | Accepted-SSPRK2-mass-flux composition transport (shared faces, reflux, MPI, retry) NOT implemented | NOT_RUN |
| M03-03 run-wide checkpoint | IMPLEMENTED_UNVALIDATED | `RestartMagneticSchemaSelection`, `restartPayloadIsMagnetized`, writer/verifier selection, `MhdRestartGroupMetadata::allow_empty_rank`, `empty_rank` marker in `writeMhdRestartGroup`/`readMhdRestartGroup`, output runtime sets the selection from `config.mhd.enabled`, restart event reports v26/v27/v25 run-wide, collective resume admission in `validateRestartResumeTopologyOrThrowImpl`; `radiationRestartSchema` declared publicly | `OutputRestartRuntime` restart write; restart resume validation | Per-rank v26 files; empty rank = valid empty group | `test_restart_checkpoint_roundtrip.cpp::testMagnetizedEmptyRankRestartV26` (+ existing v26/v25 cases) | Real distributed empty-rank restart; v27 zero-rank selection still rank-local (RT campaign, recorded) | NOT_RUN |
| M03-04 pre-admission | PARTIAL | Directory/decode/row-map reservation before `allgatherBytesBounded` (counts agreed collectively); seed layout (`kCanonicalPersistent`) + edge scratch reserved before allocation, commit only after collective success, `reconcileBaselineOwnedAndRelease`; topology item metadata pre-admitted in `buildTopology` from exact neighbour hits (`MhdLevelStepperConfig::topology_governor`), post-check, handover in `admitMemory` (no double count) | bridge `buildDirectory`/`seedFreshField`/`bind` | MemoryGovernor reservations owned by bridge/stepper | none new (existing `testMemoryRefusalBeforeAllocation`) | Per-seam fault-injection tests; table provider admission before construction (A10); HDF5 provider read buffers; local descriptor vector before admission (L-scale) | NOT_RUN |
| M03-05 bounded MPI staging | IMPLEMENTED_UNVALIDATED | `MhdLevelStepper::planExchangeRounds` (rank-invariant summed per-destination demand, header allowance), `runItems` round wrapper + `runItemsRound`, `roundsFor`, unique edge-fine phase ids, `planMemory` uses max round, `admitMemory` refuses irreducible destination > budget, `exchange_rounds` telemetry | every stepper exchange | Send/recv buffers retain max-round capacity only | `testBoundedExchangeRoundsMatchUnsplit` (bitwise, 4 KiB budget) | Allgather directory remains (small, admitted); no per-peer credit pool distinct from destination ranges | NOT_RUN |
| M03-06 peak coexistence | PARTIAL | Row maps, directory, seed, topology now charged before allocation; boundary integrals admitted via `advanceOptions` | — | Global retry copy + stepper rollback + bridge staging still coexist | none | Formal lifetime reconciliation / removal of retained bridge staging; MP02/04/05/07 measurements | NOT_RUN |
| M03-07 regrid / migration | FAIL_CLOSED | Guards unchanged in `mhdProductionDispatchRefusal`; `bindFaces` refuses topology change | — | — | existing admission tests | Governed gas+face refine/coarsen/migration transaction not implemented | NOT_RUN |
| M03-08 UCT vs face flux | IMPLEMENTED_UNVALIDATED | R-UCT: `MhdFaceSignalPatch::induction_correction`, `MhdStepOperator::evaluateIdealStage` (`plmInductionHll`, ghost-row faces solved with the same solver), `MhdEmfBuilder::buildUctHll` correction; memory estimate 4 lanes; v7 `EmfScheme = uct_hll_r_uct_v1`; derivation in `mhd_induction.hpp` | ideal stage of every interval | stage scratch only | `tests/unit/mhd/test_mhd_uct_face_consistency.cpp` | Coarse/fine and shared-edge behaviour inherited from reconciliation (tested only through existing hierarchy tests); low-beta/discontinuity campaigns (I06) | NOT_RUN |
| M03-09 AMR transfer math | PARTIAL | `MhdAmrSynchronizer::prolongFaceFluxes` quadratic coefficient now aspect-ratio-correct (bug fix) | none (regrid refused) | — | `testAnisotropicLinearFieldExactness`, existing `testNonuniformTransferIdentity` (anisotropic parent) | Tóth–Roe equivalence not derived; random-curl cycles + energy/positivity fallback campaign; production integration (M03-07) | NOT_RUN |
| M03-10 composite seed | PARTIAL | `mhd::compositeEdgePotential` (pairwise sums over integer finest-lattice sub-edges); bridge wraps every potential | `seedFreshField` | — | none new (two-level stochastic seed test NOT written) | Magnetized-shock preset and certified imported face-flux HDF5 IC not implemented (refused: seed `none` + empty state) | NOT_RUN |
| M03-11 stage domains / init | IMPLEMENTED_UNVALIDATED | `RuntimeResourceKey::kMagneticFaceState` (hydro RW, output/analysis R; declarations + view table); `HydroAmrRuntime::initializeMagneticState`; TimeCoordinator calls it before the first retry capture; bridge `ensureMagneticState(..., allow_seed)` refuses seeding from candidate/advance | `TimeCoordinator::runRungZeroSegment` | magnetic faces written only by kHydroUpdate owner | none new | Hierarchical segment path not wired (MHD refuses rung > 0) | NOT_RUN |
| M03-12 boundary budgets / output options | PARTIAL | `advanceOptions` produces boundary flux integrals when any domain side is non-periodic; accepted open/wall mass, ideal total energy and Poynting work reduced into `MhdRuntimeTelemetry`; `output_divergence_diagnostics=false` refused | bridge `advance` | stepper boundary planes (admitted) | `testConfigAndCapabilityAdmission` (divergence option) | `/MHD/FaceFlux` not implemented (no requesting contract); v7 readback test; budgets not yet in a conservation report | NOT_RUN |
| M03-13 geometry / sampler | IMPLEMENTED_UNVALIDATED | `bindFaces` checks x/y/z widths and origins; `MhdSnapshotFieldSampler` sorted index + level/dims check | bind; snapshot v7 writer | metadata-scale index (not governed) | none new | Governed sampler allocation | NOT_RUN |
| M03-14 telemetry | PARTIAL | attempted/accepted/rejected totals, exchange phases/rounds, attempted face counts, hierarchy generation counters | magnetized branch counters | — | none | Per-class high water by phase, per-peer waits, RSS/PSS deltas, OpenMP determinism audit | NOT_RUN |
| M03-15 CMake labels | IMPLEMENTED_UNVALIDATED | `cosmosim_append_test_labels(... mhd)` for all 8 MHD tests + restart roundtrip | CTest | — | static inspection | real MPI test registration | NOT_RUN |
| M04-01 neutrino retry | IMPLEMENTED_UNVALIDATED (combined MHD path FAIL_CLOSED) | `NeutrinoSpectralSource::{begin,commit,rollback}Transaction` bounded journal (opening rows, popped accepted nodes, first-use status), charged in its reservation; `GravityRuntimeImpl` cycle transaction hooks + first-use flag restore; `synchronizeModuleSidecarsAfterRollback` in multirate rollback; MHD + linear_response refused | gravity cycle transaction (global step and multirate) | live history + sidecar + journal | `test_neutrino_methods.cpp::testLinearResponseTransactionRollback` | Workflow-level injected rejection after PM refresh, restart replay, multi-rank | NOT_RUN |
| M04-02 species observables | FAIL_CLOSED | refusal of observables/instruments + composition with MHD | admission | — | admission test | depends on M03-02 | NOT_RUN |
| M04-03 observable lifecycle | FAIL_CLOSED | same refusal; observable kernels untouched | admission | — | admission test | accepted-boundary publication / staging | NOT_RUN |
| M04-04 magnetic sampler | NOT_APPLICABLE (refusals preserved) | no change; CR/RT-MHD/AGN-magnetic/RM refusals intact | — | — | — | Read-only CT sampler for future radio/CR work not written | NOT_RUN |
| M04-05 cosmology seam | IMPLEMENTED_UNVALIDATED (verification source) | Source check: workflow builds one `LambdaCdmBackground` via `cosmologyBackgroundConfigFromSimulation`; timeline (`CosmologicalTimeline`) and drag-free kick quadrature use it | — | single background | `test_mhd_ideal_reference.cpp::testExpansionWithRelativisticBackground` | Coupled wave/gravity case, MHD-off hash check | NOT_RUN |

## 2. REF gates (all NOT_RUN)

| ID | Status | Test source | Remaining |
|---|---|---|---|
| G01 | PARTIAL | `testZeroFieldLimit` | hydro-equivalence order study |
| G02 | PARTIAL | discrete complex tests, R-UCT oblique divergence | MPI/AMR topologies |
| G03 | PARTIAL | per-step hierarchy checks | 1000 steps with dynamic AMR+MPI (regrid refused) |
| G04 | PARTIAL | `testTwoLevelConservationAndCt` | open-boundary budget closure using new telemetry |
| G05 | PARTIAL | fallback assertions | certified smooth problems |
| G06 | PARTIAL | in-process two-rank + stale epoch + bounded rounds | real `mpiexec` bridge test |
| G07 | PARTIAL | transfer helper tests (+ anisotropic exactness) | production regrid (FAIL_CLOSED) |
| G08 | PARTIAL | v26 roundtrip, empty-rank v26, corruption | v7 readback |
| G09 | PARTIAL | memory governance tests | long-run plateau |
| I01 | PARTIAL | `testAlfvenConvergence` | oblique |
| I02 | OPEN | — | eigenmode suite not written |
| I03 | PARTIAL | `testBrioWu` | high-res reference |
| I04 | PARTIAL | patch equivalence | spectral convergence |
| I05 | PARTIAL | field loop | refinement campaign |
| I06 | OPEN | — | blast/rotor/current sheet |
| I07 | OPEN | — | dynamo |
| I08 | OPEN | — | dynamic crossing (regrid refused) |
| I09 | IMPLEMENTED_UNVALIDATED | `testExpansionOnly`, `testExpansionWithRelativisticBackground` | production expansion run |
| I10 | OPEN | — | coupled cosmological waves/gravity |
| N01 | PARTIAL | Ohmic decay | order campaign |
| N02 | OPEN | — | ambipolar C-shock |
| N03 | PARTIAL | Hall whistler | — |
| N04 | PARTIAL | Hall energy | — |
| N05 | PARTIAL | table tests | variable-grain/restart |
| N06 | OPEN | smoke only | variable-coefficient expanding MMS |
| N07 | PARTIAL | RKL recurrence | reference comparison |
| N08 | OPEN | — | magnetized collapse |

## 3. PERF gates (all NOT_RUN)

| ID | Status | Note |
|---|---|---|
| MP01 | PARTIAL | estimator updated for R-UCT lanes |
| MP02 | PARTIAL | scratch cap check unchanged; +2 lanes per face family |
| MP03 | IMPLEMENTED_UNVALIDATED | bounded rounds; test source |
| MP04 | OPEN | no host measurement harness |
| MP05 | OPEN | 100-interval plateau not written |
| MP06 | PARTIAL | pre-admission seams; no production regrid |
| MP07 | PARTIAL | MHD-off paths untouched except new no-op virtuals and an absent resource key |
| MP08 | PARTIAL | streamed codec; empty-rank group |
| MP09 | PARTIAL | formula sources |
| MP10 | NOT_APPLICABLE | GPU refused (CPU qualification) |
| MP11 | PARTIAL | policy unchanged |
| MP12 | PARTIAL | new counters; no RSS/PSS delta or per-peer wait |
| BP01–BP10 | PARTIAL | benchmark source only, external references `not_evaluated`; nothing measured |

## 4. Derivations and decisions (source-based; no exactness claim beyond stated)

* **R-UCT (M03-08).** `E_c = E_UCT-HLL - [a_b+ D_a|_(b-) + a_b- D_a|_(b+)]/(a_b+ + a_b-) + [a_a+ D_b|_(a-) + a_a- D_b|_(a+)]/(a_a+ + a_a-)`, `D_n(B_t) = F_face(B_t) - F_HLL^PLM(B_t)`. 1-D along `a`: equal `D_a` and fans on both a-faces, `D_b = 0` (identical L/R states), `E_UCT-HLL = -F_HLL^PLM` ⇒ `E_c = -F_face`. Divergence preserved (any edge field). Mignone & Del Zanna coefficients were NOT used (not reproduced from memory).
* **Prolongation (M03-09).** X_d-linear part of `div B` in physical coordinates fixes `a_dd = -(Δ_d/2)(b_{t1,d}/Δ_t1 + c_{t2,d}/Δ_t2)`; constant part vanishes iff the parent cell is solenoidal. Only the divergence/sum/sibling properties are claimed.
* **Composite seed (M03-10).** Pairwise integer-lattice sums make coarse = left half + right half with identical association.
* **Rounds (M03-05).** Contiguous destination ranges over dst-sorted lists keep the exact application order of the unsplit exchange.
* **Epoch (M03-01).** Integer-only digest (no struct padding or floating bits) + generation incremented on every collectively decided rebuild.

## 5. Runtime modes

Entering source paths (UNVALIDATED): synchronized static leaf AMR, serial and
MPI, ideal and constant/table non-ideal, legacy cooling sources, gravity +
cosmology (incl. relativistic background, background_only/transfer_ratio
neutrinos), v26/v7 IO. Refused: Part-2 composition, observables/instruments,
linear-response neutrinos, regrid, multi-rank rebalance, subcycling, chemistry
provider, GPU, hydro shadow, effective EOS/dual energy, RT+MHD, CR, AGN magnetic
injection, `output_divergence_diagnostics=false`.

## 6. Files touched

Core/config: `src/core/config.cpp`. MHD: `include/cosmosim/mhd/{mhd_induction,mhd_initial_conditions,mhd_level_stepper,mhd_memory}.hpp`, `src/mhd/{mhd_amr,mhd_induction,mhd_initial_conditions,mhd_level_stepper,mhd_level_topology,mhd_memory,mhd_step_operator}.cpp`. IO: `include/cosmosim/io/{mhd_restart_state,restart_checkpoint}.hpp`, `src/io/{mhd_restart_state,restart_checkpoint,snapshot_hdf5}.cpp`. Neutrinos: `include/cosmosim/neutrinos/neutrino_linear_response.hpp`, `src/neutrinos/neutrino_linear_response.cpp`. Workflows: `include/cosmosim/workflows/{gravity_runtime,hydro_amr_runtime,runtime_resources}.hpp`, `src/workflows/{gravity_runtime,hydro_amr_runtime,output_restart_runtime,reference_runtime_composition,runtime_capabilities,runtime_module_registry,time_coordinator}.cpp`, `src/workflows/internal/mhd_runtime_bridge.{hpp,cpp}`. Tests: `tests/unit/mhd/test_mhd_uct_face_consistency.cpp` (new), `tests/integration/{test_mhd_amr_mpi_restarts,test_mhd_hierarchy_runtime,test_mhd_ideal_reference,test_restart_checkpoint_roundtrip}.cpp`, `tests/unit/test_neutrino_methods.cpp`. Build: `CMakeLists.txt`. Docs: configuration, output_schema, restart_checkpointing, profiling, validation_plan, repair_open_issues, repair_state_recap, architecture/{decision_log,runtime_truth_map}, this ledger, the audit, two bundle prompts.

## 7. Continuation (precise next tasks, in order)

1. Owner compile cycle (not executed here): `cmake --preset cpu-only-debug`, `cmake --build --preset build-cpu-debug`, fix compile errors; then `ctest --preset test-cpu-debug -L mhd --output-on-failure`, `ctest --preset test-cpu-debug -R "neutrino_methods|restart_checkpoint_roundtrip|runtime_module_registry"`; then HDF5 and MPI presets from `docs/build_instructions.md`.
2. M03-02: per-face time-integrated accepted mass flux lanes per stage in `MhdStepOperator`/`MhdLevelStepper` (reconciled like gas records), donor composition at each SSPRK2 stage, conservative species update, then lift the composition guard (observables stay gated by M04-03).
3. Real `mpiexec` production-bridge test (1/2/4 ranks, zero-gas rank, stale header) with MPI CTest registration.
4. M03-07 magnetic regrid/migration transaction on `amr_production_regrid`; then I08/G07 production tests.
5. Magnetized-shock preset + certified face-flux HDF5 IC; two-level stochastic seed test.
6. MP05 100-interval plateau, MP12 RSS/PSS and per-peer waits, M03-06 lifetime reconciliation; v7 readback; `/MHD/FaceFlux` with a typed request.
7. Remove the MHD + linear_response guard after a combined injected-retry run.
