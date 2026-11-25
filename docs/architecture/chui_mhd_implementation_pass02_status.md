# CHUÍ MHD — implementation ledger, pass 02 (production-closure campaign)

**Status vocabulary:** `IMPLEMENTED_UNVALIDATED`, `PARTIAL`, `FAIL_CLOSED`,
`NOT_IMPLEMENTED`, `NOT_APPLICABLE`. Nothing in this ledger is `PASSED` or
`CERTIFIED`. **No configure, build, test, benchmark or simulation was run in
this campaign** (explicit user-reserved source-only exception to the AGENTS.md
§10 evidence floor). Compile viability of every new or changed file is
unverified; every test written here is `NOT_RUN`.

Origin tags: `PRE_EXISTING` (before pass 01), `PASS01` (introduced by commit
4382183 / 7e69a3d), `PASS02` (introduced here), `UNKNOWN`.

## 0. Baseline

| Field | Value |
|---|---|
| Repository | `/home/user/Chui-dev` (local Git clone) |
| Branch | `claude/blissful-cray-gjtfju` |
| HEAD at start | `2649fa8c07bc4191bd3bba46acc27ef79d8708e9` (pass-01 work is committed in `4382183`) |
| Dirty state at start | clean |
| Inputs added (no overwrite) | `agent_prompts/chui_mhd_production_closure_prompt_v01.md`, `docs/architecture/chui_mhd_production_closure_audit_20261124.md` |
| Audit anchors | A01–A17 re-checked against the live tree (§2) |
| Commit / push / reset | none (user instruction) |

## 1. Slice log

(updated after every slice; most recent last)

- S0 — orientation: AGENTS.md, REF, PERF, amendment, pass-01 ledger,
  interruption audit, production-closure audit read in full. Live anchors
  confirmed: `mhd_level_stepper.cpp` 92–218 uniform tiling and O(P²) id scan;
  1157–1210 `canonicalizeMirrors()` writes before validation and before the
  rollback image; 1248 catches only `MhdError`; 1259–1272 internal shortened
  accepted interval; 1294–1298 estimator calls the mutating canonicalizer;
  `hydro_amr_runtime.cpp` 5883 blanket refusal; `runtime_capabilities.cpp`
  420–435 blanket refusal.

- S1 — `MhdLevelStepper` rewritten as a transactional, topology-aware
  synchronized leaf-hierarchy driver (`include/cosmosim/mhd/mhd_level_stepper.hpp`,
  `src/mhd/mhd_level_stepper.cpp`, new `src/mhd/mhd_level_topology.cpp`).
- S2 — expansion ownership (A3): drag-free source kick quadrature
  (`core::computeComovingSourceKickFactor`, `CosmologicalStepFactors::
  {first,second}_source_kick_factor_code`), magnetized gas kick
  `applyMagnetizedGasCellSourceKick` in `src/workflows/gravity_runtime.cpp`.
- S3 — runtime composition (B1–B3): `src/workflows/internal/mhd_runtime_bridge.{hpp,cpp}`,
  magnetized `kHydroUpdate` branch and `magnetizedStableTimeStep` in
  `src/workflows/hydro_amr_runtime.cpp`, scheduler candidate in
  `src/workflows/time_coordinator.cpp`, precise admission
  `mhdProductionDispatchRefusal` in `src/workflows/runtime_capabilities.cpp`,
  `[mhd_io] initial_field_*` seed keys in `include/cosmosim/core/config.hpp`,
  `src/core/config.cpp`.
- S4 — restart v26 (E1) in `src/io/restart_checkpoint.cpp`
  (`mhdRestartSchema`, `/state/mhd` write/read, SHA-256 coverage, verifier),
  codec streaming/bounds/spelling in `src/io/mhd_restart_state.cpp`.
- S5 — snapshot v7 (E2) in `src/io/snapshot_hdf5.cpp`
  (`MhdSnapshotFieldSampler`, `PartType0/MagneticField`, `MagneticPressure`,
  `/MHD/Diagnostics`, `mhdLocalMagneticEnergyComoving`), rank-agreed
  diagnostics in `src/workflows/output_restart_runtime.cpp`.

## 2. Audit findings — disposition (A01–A17)

| ID | Disposition | Source evidence (function) | Status |
|---|---|---|---|
| A01 | Blanket refusal replaced by `mhdProductionDispatchRefusal` (precise predicates) + real dispatch `HydroAmrRuntimeImpl::execute` → `MhdRuntimeBridge::advance` | runtime_capabilities.cpp, hydro_amr_runtime.cpp | IMPLEMENTED_UNVALIDATED |
| A02 | Sparse leaf hierarchy: box-intersection plans over a sorted bin index; arbitrary rectangular extents, nonzero domain origin (`MhdLevelDomain::origin_comoving`), partial fine coverage, zero-patch ranks; O(P log P) duplicate check | `buildTopology`, `buildSpatialIndex`, `queryLevel`, `planHalo`, `planReconciliation` | IMPLEMENTED_UNVALIDATED |
| A03 | Read-only `validateCanonicalConsistency` BEFORE the rollback image; normalization inside the transaction; RAII `RollbackGuard` + no-throw `restoreRollback`; all exception types caught and collectively agreed (`agree`); estimator never normalizes | `advance`, `agree`, `reconcileFaces(compare_only)` | IMPLEMENTED_UNVALIDATED |
| A04 | `MhdRetryMode::kExternal` (default, production): retryable failure → rollback, `accepted=false`; the bridge throws `StageTransactionRejected(kNumericalRetry)` into the existing global-step retry transaction; admission requires `numerics.global_step_max_retries >= 1` | `advance`, `MhdRuntimeBridge::advance`, `mhdProductionDispatchRefusal` | IMPLEMENTED_UNVALIDATED |
| A05 | MHD stage is the single expansion owner for magnetized gas; KDK gas kick uses drag-free `integral dt/a^2`; particles and MHD-off gas unchanged | `applyGasKick`, `applyMagnetizedGasCellSourceKick`, `computeComovingSourceKickFactor` | IMPLEMENTED_UNVALIDATED (see R-EXP) |
| A06 | v26 writer/reader and v7 snapshot wired | `writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5`, `writeScienceSnapshotHdf5` | IMPLEMENTED_UNVALIDATED |
| A07 | Synchronized leaf AMR reconciles coarse/fine per stage (records, edges, faces) so no pending register exists; accepted-interval boundary gas (Σ w F/a) and Poynting (Σ w F_E) integrals follow SSPRK2/explicit/Hall/RKL2 recurrences when `produce_boundary_flux_integrals` | `idealStage`, `explicitNonIdealSubstep`, `rkl2Cycle`, `combineRecords` | IMPLEMENTED_UNVALIDATED |
| A08 | Batch slots (≤ `batch_patches`, default 8) for tile-local scratch; level-resident state separated; honest per-lane plan incl. halo imports, topology, comm; admission before allocation; commit of measured capacity; growth check | `planMemory`, `admitMemory`, `ensureAllocated`, `forEachLocalPatchBatched` | IMPLEMENTED_UNVALIDATED |
| A09 | Flux lanes streamed in ≤ 8 MiB hyperslabs straight into the candidate lanes; patch-table bound; governed candidate reservation | `readFluxLaneInto`, restart reader v26 block | IMPLEMENTED_UNVALIDATED |
| A10 | Not changed in this pass (table provider materialization); admission of table bytes not yet charged | — | PARTIAL (open) |
| A11 | `cosmosim_workflows` links `cosmosim_mhd` PUBLIC | CMakeLists.txt | IMPLEMENTED_UNVALIDATED |
| A12 | Reviewed: UCT-HLL transverse/corner states are PLM-MC in both PLM and PPM modes (pass-01 R1); the exact 1-D reduction holds only for the HLL face flux with PLM-MC states, NOT for HLLD or PPM states. Declared, not claimed general | `MhdEmfBuilder::buildUctHll` (unchanged) | PARTIAL (declared reconciliation; tests NOT_RUN) |
| A13 | Not re-derived in this pass; the Balsara = Tóth–Roe ratio-2 claim in `mhd_amr.hpp` remains an unverified assertion (dynamic regrid is fail-closed, so the transfer is not on a production path) | — | PARTIAL (open) |
| A14 | Level-wide: all patches loaded, advanced, validated, then stored; temperature/sound speed refreshed by `hydro::refreshGasThermodynamicCell`; gravity source generation bumped | `MhdRuntimeBridge::advance` | IMPLEMENTED_UNVALIDATED |
| A15 | Canonical spelling frozen as `face_flux_{x,y,z}_comoving` (repository naming rule, amendment R4); `_comov` rejected, not aliased | `readMhdRestartGroup` | IMPLEMENTED_UNVALIDATED |
| A16 | MPI halo/EMF/record exchange via `MpiMhdHaloTransport` over `exchangeBoundedAlltoallBytes`; OpenMP over batch slots; telemetry struct | `runItems`, `exchangeHalo`, `MhdRuntimeTelemetry` | IMPLEMENTED_UNVALIDATED (counters.json export open) |
| A17 | Test sources added/updated (see §4) | tests/ | NOT_RUN |

## 3. Requirement coverage (source evidence; nothing validated)

### 3.1 Prompt work packages

| ID | Requirement | Status | Evidence |
|---|---|---|---|
| A1 | Kernel contracts re-audited | PARTIAL | No kernel math changed; UCT-HLL reduction scope declared (A12); A13 open; new hierarchy test exercises PPM + uniform invariance across interfaces |
| A2 | Atomic interval (read-only validation, RAII rollback, all exceptions, attempted vs accepted) | IMPLEMENTED_UNVALIDATED | `MhdLevelStepper::advance`, `validateCanonicalConsistency`, `restoreRollback` (noexcept), `agree`; tests `testEntryMirrorMismatchIsReadOnly`, `testExternalRetryCommitsNothing`, `testForeignLateFailureRollsBack` |
| A3 | Single expansion owner | IMPLEMENTED_UNVALIDATED | R-EXP; `applyMagnetizedGasCellSourceKick`, `computeComovingSourceKickFactor`; I09/I10 cosmological runtime tests NOT_WRITTEN (kernel-level I09 exists) |
| A4 | Non-ideal at stage properties, eta/a² once, caps | IMPLEMENTED_UNVALIDATED (pass-01 kernels) | coefficient evaluation per stage in `evaluateNonIdealRate`; caps in `diffusionHalf`/`hallHalf`/`forecastNonIdeal`; no kernel change in pass 02 |
| B1 | Real runtime composition | IMPLEMENTED_UNVALIDATED | `MhdRuntimeBridge` (directory, lattice, row maps from cell centres, seed vs restart, load/advance/store, thermo refresh); `HydroAmrRuntimeImpl::execute` MHD branch |
| B2 | Scheduler candidate, external retry, Strang sources | IMPLEMENTED_UNVALIDATED | `TimeCoordinator::runRungZeroSegment` MHD candidate; kExternal retry → global-step transaction; existing S(dt/2) H(dt) S(dt/2) KDK order (source owners act on gas lanes only; magnetic faces are never touched by sources — explicit policy: star formation / feedback / sinks change gas mass/energy lanes, B_c and the face flux stay, magnetic energy unchanged) |
| B3 | Capability-specific guards | IMPLEMENTED_UNVALIDATED | `mhdProductionDispatchRefusal`; positive/negative test in `test_mhd_amr_mpi_restarts.cpp` |
| C1 | Topology-aware driver | IMPLEMENTED_UNVALIDATED | `buildTopology`/`planHalo`/`planReconciliation`; sparse coverage, arbitrary extents, nonzero origin, zero-patch ranks; deterministic keys; collective error agreement before every exchange |
| C2 | Synchronized coarse/fine | IMPLEMENTED_UNVALIDATED | R-SYNC; `planCrossLevel`, `reconcileGasRecords/EnergyRecords/Edges/Faces` fine overrides; `testTwoLevelConservationAndCt`, `testTwoLevelUniformInvariance` |
| C3 | Regrid / migration | FAIL_CLOSED | refused by `mhdProductionDispatchRefusal` (`amr_regrid_mode`, multi-rank rebalance); bridge detects topology change without magnetic transfer (`bindFaces`) |
| D1 | Bounded batching, honest accounting, admission | IMPLEMENTED_UNVALIDATED | `planMemory` (per-lane), `admitMemory` (phase/scratch/topology/comm), `ensureAllocated` commit + growth check, bridge gas staging reservation, restart candidate reservation |
| D2 | OpenMP + telemetry | IMPLEMENTED_UNVALIDATED | `forEachLocalPatchBatched` (`#pragma omp parallel for`, slot-ordered reductions); counters `mhd_*` in `HydroAmrRuntimeImpl::execute` |
| E1 | v26 restart | IMPLEMENTED_UNVALIDATED | writer/reader/verifier; `testMagnetizedRestartV26` |
| E2 | v7 snapshot | IMPLEMENTED_UNVALIDATED | `writeScienceSnapshotHdf5` v7 branch; no dedicated v7 readback test written (OPEN) |
| F | Verification sources, docs | PARTIAL | §6 |

### 3.2 Supported capability matrix (source; UNVALIDATED, NOT_RUN)

| Capability | Status |
|---|---|
| Serial fixed (single-level leaf tiling) | IMPLEMENTED_UNVALIDATED |
| Serial AMR (static 2:1 leaf tiling, synchronized) | IMPLEMENTED_UNVALIDATED |
| MPI fixed | IMPLEMENTED_UNVALIDATED (protocol tested in-process only, NOT_RUN) |
| Synchronized MPI AMR (static) | IMPLEMENTED_UNVALIDATED |
| Dynamic regrid / migration with MHD | FAIL_CLOSED |
| Distributed recursive subcycling (`mhd_ref_amr_subcycled`) | FAIL_CLOSED |
| Non-ideal Ohmic / ambipolar (SSPRK2/RKL2), Hall-HLL | IMPLEMENTED_UNVALIDATED (constant and verified-table providers) |
| Chemistry provider | FAIL_CLOSED |
| Gravity (TreePM kick) + cosmological expansion | IMPLEMENTED_UNVALIDATED (R-EXP) |
| Cooling / star formation / feedback via S(dt/2) hooks | IMPLEMENTED_UNVALIDATED (sources act on gas lanes; faces untouched) |
| Effective multiphase EOS, dual energy | FAIL_CLOSED (config) |
| Restart v26 (same rank topology) | IMPLEMENTED_UNVALIDATED |
| Science snapshot v7 | IMPLEMENTED_UNVALIDATED |
| Imported face-flux HDF5 IC | NOT_IMPLEMENTED (seeds only) |
| GPU | FAIL_CLOSED |
| Hydro shadow reference with MHD | FAIL_CLOSED |

## 4. Reconciliation notes (decisions recorded, not invented)

- **R-EXP (expansion ownership).** For magnetized runs the MHD ideal stage owns
  `-H m_c` and `-H(2K_c+3p_c+M_c)` at its SSPRK2 stage endpoints using the
  timeline's `a, H` at the interval endpoints (`MhdStageContext::
  apply_expansion_source=true`). The gas KDK kick of magnetized runs adds only
  `du = A · ∫ dt/a²` (new drag-free timeline quadrature, 64-point midpoint in
  `a`, same sampling as the existing kick integral) and applies no thermal
  scaling. The existing contract `u1 = (a0/a1) u0 + A/a1 ∫dt/a` is unchanged
  for particles and MHD-off gas. Gravitational work enters through the kinetic
  energy of the kicked velocity when `E_c` is rebuilt at the next
  `kHydroUpdate`. Order: K(dt/2) D S(dt/2) [MHD incl. expansion](dt) S(dt/2) K(dt/2).
- **R-SYNC (coarse/fine).** CHUI production AMR is a leaf tiling (no covered
  coarse cells). With synchronized sweeps every coarse/fine interface lies on
  patch boundaries; per stage, coarse interface gas/Poynting records take the
  mean of the four fine records, coarse interface edges the sum of the two
  fine edges, coarse interface faces the sum of the four fine faces (fine side
  authoritative). This is the per-stage form of reflux + reflux-curl; no
  register outlives a stage, so every step boundary is MHD-checkpoint safe.
  Ghosts: same-level copy; finer → conserved restriction; coarser → minmod
  linear prolongation of primitives (positivity-preserving), ghost staggered
  faces by face sampling. Requires a 2:1 balanced leaf tiling (fail closed).
- **R-SPELL.** v26 `/state/mhd` flux datasets are `face_flux_*_comoving`.
- **R-MU.** Non-ideal property model uses the Part-1 reference μ = 0.59.
- **R-IC.** Fresh starts seed `Φ = Φ_uniform + C(A)` on the real leaf topology
  from `[mhd_io] initial_field_seed`; imported face-flux HDF5 ICs are not
  implemented (seed `none` + empty magnetic state fails closed).

## 5. Remaining guards and why

| Guard | Location | Reason it stays |
|---|---|---|
| `amr_regrid_mode != off` refused | `mhdProductionDispatchRefusal` | `amr_production_regrid` has no magnetic transfer/coexistence transaction yet (C3) |
| multi-rank runtime rebalance refused | same | migration payloads carry no `/state/mhd` face records yet |
| recursive / distributed subcycling refused | same + config | no magnetic temporal histories (REF 7.4 tier) |
| `global_step_max_retries < 1` refused | same | retry ownership must be external (A04) |
| chemistry provider, GPU, hydro shadow refused | same | no validated contract / no kernels / nonmagnetic comparator |
| fixed-Cartesian (no patch geometry) refused | `HydroAmrRuntimeImpl::execute`, `magnetizedStableTimeStep` | face-flux authority needs explicit leaf geometry |
| topology change without magnetic transfer | `MhdRuntimeBridge::bindFaces` | never reseed or silently remap a field |
| 2:1 balance / tiling gaps | `planHalo`, `buildTopology` | prolongation/restriction contract |
| v25 writer for magnetized state | routed to v26; pre-v26 `/state/mhd` rejected | amendment R1 |
| v6 snapshot for magnetized config | routed to v7; inconsistent state/config refused | amendment R6 |

## 6. Gates and test sources (all NOT_RUN)

| Gate | Pass-02 test source |
|---|---|
| G02/G03 (per step, composite leaf complex) | `test_mhd_hierarchy_runtime.cpp` `runSteps` divergence check; two-level |
| G04 closed conservation incl. coarse/fine | `testTwoLevelConservationAndCt` |
| G06 rank equivalence (protocol) | `testTwoRankEquivalence` (in-process transport, incl. zero-patch rank); real `mpiexec` run NOT_WRITTEN |
| G07 regrid | FAIL_CLOSED (no production path) |
| G08 I/O | `testMagnetizedRestartV26` (v26 roundtrip, SHA-256 corruption, spelling, v25 unchanged, unsafe boundary); v7 readback NOT_WRITTEN |
| MP02/MP06 | `testMemoryRefusalBeforeAllocation`, existing `test_mhd_memory_governance.cpp` |
| Transactions (A2) | `testEntryMirrorMismatchIsReadOnly`, `testExternalRetryCommitsNothing`, `testForeignLateFailureRollsBack` |
| Topology (A02) | `testNonUniformExtentsAndOrigin` (arbitrary extents, nonzero origin, batch-size invariance) |
| Admission (B3) | `test_mhd_amr_mpi_restarts.cpp::testConfigAndCapabilityAdmission` |
| I02, I06, I08, N02, N08, MP03, MP05 (100 intervals), MP12, BP campaigns, cloud collapse / dynamo decks | NOT_WRITTEN in pass 02 (carried open) |

Modified existing tests: `test_mhd_ideal_reference.cpp::testRetryAndRollback`
(external default + explicit internal mode), `mhd_test_fixtures.hpp`
(config overload), `test_mhd_amr_mpi_restarts.cpp` (precise admission).

## 7. Recommended owner build order (not executed)

1. `cmake --preset cpu-only-debug && cmake --build --preset build-cpu-debug`
   (fix compile errors first: pass-02 files never compiled).
2. `ctest --preset test-cpu-debug -L mhd --output-on-failure`, then
   `integration_mhd_hierarchy_runtime`, `integration_restart_checkpoint_roundtrip`.
3. MHD-off regression: config parser/examples, restart v25 roundtrip,
   snapshot v6, gravity KDK tests (the new step fields must not change MHD-off
   kick factors).
4. HDF5 + MPI presets; add an `mpiexec -n 2` run of the canonical deck.

## 8. Interruption checkpoint / next exact work item

All pass-02 source edits are coherent (no half-edited function); every
unsupported combination remains fail-closed with a precise refusal. Next
exact items, in order: (1) owner compile/test cycle (§7); (2) C3 — wire
`MhdAmrSynchronizer` prolong/restrict into `amr_production_regrid` with
old+new coexistence admission and magnetic migration payloads, then lift the
regrid/rebalance predicates; (3) imported face-flux HDF5 IC reader;
(4) v7 readback test and real `mpiexec` MHD test; (5) A13 derivation/test and
A10 table admission; (6) remaining G/I/N/MP/BP sources of §6.
