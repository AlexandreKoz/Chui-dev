# CHUÍ MHD — implementation ledger, pass 01 (continuation campaign)

> Superseded for current status by
> `docs/architecture/chui_mhd_implementation_pass02_status.md` (pass 02:
> production dispatch, leaf-AMR/MPI stepper, v26/v7). This file is kept as the
> pass-01 record.

**Status vocabulary:** `IMPLEMENTED_UNVALIDATED`, `PARTIAL`, `FAIL_CLOSED`,
`NOT_IMPLEMENTED`, `NOT_APPLICABLE`. Nothing here is `PASSED` or `CERTIFIED`.
**No configure, build, test, benchmark or simulation was run in this
campaign** (explicit user-reserved source-only exception to the AGENTS.md §10
evidence floor). Compile viability of every new or changed file is unverified.
Written tests are `NOT_RUN`.

Origin tags for defects: `PRE_EXISTING` (before the OpenCode MHD commit),
`INTRODUCED` (introduced by OpenCode commit 7e69a3d), `UNKNOWN`.

## 0. Baseline

| Field | Value |
|---|---|
| Repository | `/home/user/Chui-dev` (local Git clone) |
| Branch | `claude/upbeat-lamport-fyk0p5` |
| HEAD at start | `7e69a3d4d6f83b794a8d1207a5178c2c95f65121` |
| Dirty state at start | clean |
| Untracked inputs added (no overwrite) | `agent_prompts/chui_mhd_claude_continuation_v01.md`, `agent_prompts/chui_mhd_claude_continuation_readme_v01.md`, `docs/architecture/chui_mhd_open_code_interruption_audit_20261124.md` |
| Live tree vs audit | matched the audited export (20 MHD C++ files + 5 core files); all audit line anchors confirmed |
| Commit / push / reset | none (user instruction) |

## 1. What changed (by slice)

**A — perimeter, build, config bridge**
* `src/workflows/runtime_capabilities.cpp`: capabilities `mhd_ideal_reference_kernels`,
  `mhd_production_dispatch`, `mhd_nonideal_transport`,
  `mhd_chemistry_coefficient_provider`, `mhd_distributed_recursive_amr`,
  `mhd_gpu_backend`; `validateRequestedRuntimeCapabilities` rejects every
  magnetized request with `MHD_UNSUPPORTED_PHYSICS` before IC materialization.
* `src/workflows/hydro_amr_runtime.cpp`: `makeHydroAmrRuntime` refuses `mhd.enabled`.
* `CMakeLists.txt`: `cosmosim_mhd` library (13 TUs, core-only dependency,
  private HDF5 for the table provider), install/export/alias; IO codec TU;
  6 test targets; 1 benchmark target.
* `src/core/config.cpp`, `include/cosmosim/core/config.hpp`: `[mhd*]` emitted
  into normalized text only when enabled (MHD-off hash byte-identical to the
  predecessor); `MhdConfig` defaulted equality; inert non-default MHD keys on a
  disabled run rejected; cross-checks (dual energy, effective-ISM EOS,
  distributed subcycling, GPU, table composition/grain ids, mixed authority).

**B — ideal core** (`mhd_riemann`, `mhd_induction`, `mhd_step_operator`)
* HLLD rewritten from Miyoshi & Kusano (2005): symmetric p_T*, degenerate
  star states (vt*=vt, Bt*=Bt), vanishing-B_n collapse, **fixed inverted
  S_M<0 branch selection**, star thermal admissibility, signed fan bounds;
  HLLE with signed bounds; `signalBounds`; first-order Rusanov on
  **cell-average** states returning `valid=false` instead of a zero flux.
* Genuine 2-D UCT-HLL edge EMF (Londrillo & Del Zanna form) on ghost-extended
  face lattices: signed a±, transversely reconstructed staggered fields, MC
  corner velocities, HLL jump dissipation; exact 1-D reduction to the face HLL
  flux; E_c = E_eff/√a; conducting-wall E_t=0 and periodic edge copy rules.
* FV update with correct sign (lower cell loses, upper gains) and 1/(a dx);
  one cyclic axis convention for every sweep; PLM-MC stencil ±1 (ghost 2) and
  bounded **characteristic** PPM (7-wave MHD primitive eigenbasis, CW84
  limiting, counted PLM positivity fallback, ghost 3); shared normal B_eff on
  both face states; metal flux from the accepted mass flux.
* Cosmology: B_eff = B_c/√a in Riemann/primitives, M_c = B_c²/(2a) in recovery
  and budget, expansion source −H m, −H(2K+3p+M) at SSPRK2 stage endpoints
  (owner flag `apply_expansion_source`), isolated a=1/H=0 enforced, multi-D
  fast CFL with a·dx, no silent skips or floors.
* `MhdLevelStepper` (new): single committer of a synchronized same-level,
  same-rank interval — rollback image, D(dt/2) H(dt/2) I(dt) H(dt/2) D(dt/2),
  same-stage neighbour halos, owner reconciliation of shared faces, boundary
  gas/Poynting fluxes and edge EMFs, retryable vs fatal classification,
  in-stepper retries only with a timeline sampler (or constant timeline),
  optional time-integrated edge EMF (Φ_end = Φ_begin − C(I)), governed memory.

**C — non-ideal** (`mhd_transport`)
* Compatible dual-mesh edge current J = curl_x B_c from staggered faces;
  edge-local constitutive law with edge-averaged B and η; Hall-HLL with grid
  Nyquist whistler speed πη_H/Δ and reconstructed jumps; conservative face
  Poynting flux (E_c×B_c)_n/a (Joule heat appears only through energy
  conservation; diagnostic J·E separate); dt_D = 0.15/[max(η_O+η_A)Σ1/Δ²],
  dt_H = 0.1/[π max|η_H|Σ1/Δ²]; correct RKL2 recurrence (MBA 2014) with
  stability factor (s²+s−2)/4, deterministic cycle plan ≤128 stages/cycle,
  per-interval caps; explicit SSPRK2 substeps for R ≤ 4; per-stage coefficient
  re-evaluation and growth detection; property model as the single unit owner.
* Table provider: build-config include fixed (HDF5 path was always compiled
  out), contract/units metadata, shape checks, BZeroLimit group, log
  interpolation for positive Ohmic/ambipolar, signed Hall, strict coverage,
  SHA-256 of the parsed bytes, admission required before evaluate.

**D — AMR, state, memory**
* `mhd_amr`: Balsara 2001 divergence-free ratio-2 prolongation (fine fluxes are
  exact integrals; Φ/4 for uniform fields), zero-slope fallback, unique
  restriction with sibling verification, transactional gas+field prolong/restrict
  with conservative thermal re-allocation, register on full geometric keys with
  pre-declared edges, coverage masks, stale/duplicate/undeclared rejection,
  reflux-curl as curl correction.
* `mhd_memory`: exact per-lane retained-capacity model, 512 B/cell check,
  batch sizing, whole-host envelope, `MhdMemoryAdmission` on the existing
  `MemoryGovernor` (no second governor).
* Core: checked counts in `MagneticFaceSoa::isConsistent`; ownership
  invariants include magnetic faces; memory report entries (only when
  non-empty); `kMagneticFaces` write permission enforced in stage contracts.
* `mhd_runtime_adapter`: E_c includes B_c²/(2a) on load, no floors on store,
  configured γ, stale identity/row-map rejection, transactional face layout,
  real provider construction (constant/table; chemistry fails closed).

**E — IO, IC, config decks, docs**
* `io/mhd_restart_state`: `/state/mhd` v26 group codec (chunked ≤8 MiB writes,
  full pre-mutation validation). v25 writer and v6 snapshot writer refuse
  magnetized states.
* `mhd_initial_conditions`: Φ = Φ_uniform + C(A) seeds (field loop, circular
  Alfvén, Hermitian stochastic seed with SplitMix64), deterministic.
* Decks: `configs/cosmo_cube/cube_mhd_ref_v01.param.txt`,
  `configs/isolated_galaxy/disk_mhd_ref_v01.param.txt` (rejected at admission).
* Docs: configuration, output_schema, restart_checkpointing, validation_plan,
  runtime_truth_map, decision_log, repair_open_issues.

**F — test and benchmark source** (all `NOT_RUN`): see §4.

## 2. Audit findings — disposition

| ID | Finding | Origin | Disposition |
|---|---|---|---|
| P0-01 | No build/production integration | INTRODUCED | CMake registered; production dispatch `FAIL_CLOSED` (capability guard, factory guard) |
| P0-02 | Reversed FV flux signs, no 1/a | INTRODUCED | Fixed (`evaluateIdealStage`, `applyBoundaryFaceFluxes`); test `testAdvectionDirection` |
| P0-03 | Axis/stencil addressing, PLM ±2 reads, edge index overrun | INTRODUCED | Fixed (cyclic convention, ±1 PLM, ghost-extended signal/staggered lattices); tests axis rotation, ghost-3 PPM patch |
| P0-04 | Fake UCT-HLL, t1/t2 mismatch | INTRODUCED | Replaced by genuine UCT-HLL on signed fans; one generic accessor for all axes |
| P0-05 | Incomplete comoving energy | INTRODUCED | Fixed in kernels/adapter/budget/AMR; expansion test I09; source-owner reconciliation open (§6) |
| P0-06 | AMR flux mismatch, double-counted restriction, hash register | INTRODUCED | Fixed; tests uniform Φ/4, identity, sibling equality, register coverage |
| P0-07 | Non-ideal never executed | INTRODUCED | D/H/I/H/D executed transactionally in `MhdLevelStepper::advance` |
| P0-08 | Cell-centred current, −J·E energy | INTRODUCED | Replaced by compatible edge curl and conservative face Poynting flux |
| P1-09 | Wrong dt bounds, Euler "explicit", broken RKL2, iterator-to-vector, nonexistent member | INTRODUCED | Fixed (bounds, SSPRK2 substeps, MBA recurrence, caps); RKL2 scalar test |
| P1-10 | Provider/guard | INTRODUCED | Table provider wired via `makeConfiguredCoefficientProvider`; chemistry fail-closed |
| P1-11 | Ungoverned allocations | INTRODUCED | Workspaces reused (no per-stage allocation), exact model, governed admission; OpenMP loops `NOT_IMPLEMENTED` |
| P1-12 | Missing lifecycle integrations | INTRODUCED | Partially: state accounting/invariants, codec, guards; dispatch/regrid/migration/MPI/v7 open |
| NEW-01 | HLLD S_M<0 region inverted (F_R* vs F_R**) | INTRODUCED | Fixed; Galilean-continuity sweep test |
| NEW-02 | HLLD degenerate denom returned zero star state | INTRODUCED | Fixed (MK05 limit) |
| NEW-03 | First-order fallback returned a silent zero flux on invalid input | INTRODUCED | Fixed (`valid=false`, stage rejection) |
| NEW-04 | `COSMOSIM_ENABLE_HDF5` never defined in `mhd_transport.cpp` | INTRODUCED | Fixed (`build_config.hpp`) |
| NEW-05 | Divergence normalizer floor 1.0 hid small-field violations | INTRODUCED | Fixed (G03 ε-scale) |
| NEW-06 | `MhdFaceFluxMutableView::geometry` non-const vs const sources (compile error) | INTRODUCED | Fixed |
| NEW-07 | Unconditional `[mhd*]` normalized text changed every legacy config hash | INTRODUCED | Fixed (conditional emission + inert-key rejection) |

## 3. Requirement coverage

### Milestones
| ID | Status | Evidence / gap |
|---|---|---|
| M0 contracts | PARTIAL | config, state, keys, errors, schema numbers; production preflight is the guard only |
| M1 ideal single patch | IMPLEMENTED_UNVALIDATED | `mhd_riemann/induction/step_operator/level_stepper` |
| M2 accuracy + AMR | PARTIAL | char. PPM, transfer, register kernels; no multi-level composition, no MPI |
| M3 cosmology & coupling | PARTIAL | kernel sources; gravity/cooling/source composition `FAIL_CLOSED` (no dispatch) |
| M4 non-ideal | IMPLEMENTED_UNVALIDATED (kernels+stepper) | constant/table providers, D/H, Strang, RKL2, Hall-HLL; no production dispatch |
| P0 resource contract | IMPLEMENTED_UNVALIDATED | `mhd_memory` estimator/admission/envelope |
| P1 patch CPU kernels | PARTIAL | governed reuse; no OpenMP loops, no timers in profile_timeline |
| P2 AMR/MPI execution | NOT_IMPLEMENTED (MPI), PARTIAL (AMR kernels) | |
| P3 non-ideal stages | IMPLEMENTED_UNVALIDATED | reusable lanes, caps, forecasts |
| P4 production telemetry | NOT_IMPLEMENTED | counters exist in results; not exported to `counters.json` |

### General gates G01–G09
| Gate | Status | Source / test |
|---|---|---|
| G01 zero-field | IMPLEMENTED_UNVALIDATED | `testZeroFieldLimit`; hydro-order comparison vs legacy `NOT_IMPLEMENTED` |
| G02 discrete complex | IMPLEMENTED_UNVALIDATED (patch, periodic, wall, AMR children) | `test_mhd_discrete_complex`, AMR test; MPI topology `NOT_IMPLEMENTED` |
| G03 long-run divergence | PARTIAL | per-interval contract check; 1,000-step AMR/MPI run not written |
| G04 closed conservation | IMPLEMENTED_UNVALIDATED | field loop, conducting box, Ohmic/Hall closed box |
| G05 floors/fallback | IMPLEMENTED_UNVALIDATED | counters + fraction gate flag; no floors exist |
| G06 MPI equivalence | NOT_IMPLEMENTED | no MPI halo protocol (fail-closed) |
| G07 regrid | IMPLEMENTED_UNVALIDATED (kernels) | AMR test; production regrid not wired |
| G08 I/O | PARTIAL | v26 group codec + corruption tests; v26/v7 writers `FAIL_CLOSED` |
| G09 memory boundedness | PARTIAL | plateau mini-test, disabled-MHD test; repeated-run plateau not written |

### Ideal I01–I10
| Gate | Status | Test |
|---|---|---|
| I01 Alfvén | IMPLEMENTED_UNVALIDATED | `testAlfvenConvergence` (32/64; oblique not written) |
| I02 magnetosonic/entropy modes | NOT_IMPLEMENTED (test) | |
| I03 Brio–Wu | PARTIAL | positivity/CT smoke; reference-solution convergence not written |
| I04 Orszag–Tang | PARTIAL | used in multi-patch equivalence; spectra not measured |
| I05 field loop | IMPLEMENTED_UNVALIDATED | `testFieldLoopConservation` |
| I06 blast/rotor/current sheet | NOT_IMPLEMENTED (test) | |
| I07 dynamo | NOT_IMPLEMENTED | needs campaign |
| I08 AMR crossing | NOT_IMPLEMENTED | no multi-level stepping |
| I09 expansion-only | IMPLEMENTED_UNVALIDATED | `testExpansionOnly` |
| I10 cosmological waves/gravity | NOT_IMPLEMENTED | needs dispatch |

### Non-ideal N01–N08
| Gate | Status | Test |
|---|---|---|
| N01 Ohmic decay | IMPLEMENTED_UNVALIDATED | `testOhmicDecay` (128 cells, 1%) |
| N02 ambipolar wave/C-shock | NOT_IMPLEMENTED (test) | operator exists |
| N03 whistler | IMPLEMENTED_UNVALIDATED | `testHallWhistler` (phase 2%, polarity) |
| N04 Hall-only energy | IMPLEMENTED_UNVALIDATED | energy ≤1e-12, bounded J·E_H diagnostic |
| N05 table | IMPLEMENTED_UNVALIDATED | synthetic fixture test (HDF5 builds); restart SHA admission in codec |
| N06 combined D+H+I | PARTIAL | smoke + edge-integral identity; order study not written |
| N07 STS vs explicit | PARTIAL | RKL2 scalar stability/order; field-level comparison not written |
| N08 cloud collapse | NOT_IMPLEMENTED | needs dispatch + gravity |

### Performance MP01–MP12
| Gate | Status | Evidence |
|---|---|---|
| MP01 | IMPLEMENTED_UNVALIDATED | estimate == `ownedCapacityBytes` test |
| MP02 | IMPLEMENTED_UNVALIDATED | 32³ PPM ideal ≈388 B/cell, non-ideal+RKL ≈492 B/cell — arithmetic on the estimator formulas during source review, not a run of the code; asserted by `testExactCapacityAccounting` (NOT_RUN) |
| MP03 | NOT_IMPLEMENTED | no MPI buffers |
| MP04 | PARTIAL | `requireHostEnvelope`; not wired into workflow preflight |
| MP05 | PARTIAL | two-interval plateau check |
| MP06 | IMPLEMENTED_UNVALIDATED | admission rejection before mutation; regrid/migration coexistence not wired |
| MP07 | IMPLEMENTED_UNVALIDATED | empty magnetic state, no report entries, conditional normalized text |
| MP08 | PARTIAL | codec chunking ≤8 MiB; writer not wired |
| MP09 | IMPLEMENTED_UNVALIDATED | N_f formula test at 16/32/64 |
| MP10 | NOT_APPLICABLE (GPU rejected, no kernels) | |
| MP11 | IMPLEMENTED_UNVALIDATED | host envelope rule |
| MP12 | NOT_IMPLEMENTED | telemetry export |

### Benchmarks BP01–BP10
`bench/mhd/bench_mhd_reference_suite.cpp`: BP01, BP04, BP05 local timing
(`measured_local_only` once run), BP10 model only; BP02, BP03, BP06–BP09 report
`not_evaluated` with reasons. No result exists.

## 4. Test sources written (NOT_RUN)
`tests/unit/mhd/{mhd_test_fixtures.hpp, test_mhd_discrete_complex.cpp,
test_mhd_riemann_reconstruction.cpp, test_mhd_memory_governance.cpp}`,
`tests/integration/{test_mhd_ideal_reference.cpp, test_mhd_nonideal_reference.cpp,
test_mhd_amr_mpi_restarts.cpp}`. Assertions use `CHUI_REQUIRE` (active under
NDEBUG). Tolerances follow the specification; analytic references are derived
in the tests (no fabricated reference outputs).

## 5. Production call chain (source)
`.param.txt` → `loadFrozenConfig*` (typed MHD keys, normalized only when
enabled) → `ReferenceWorkflow::run` → `buildRuntimeCapabilityReport` →
`validateRequestedRuntimeCapabilities` → **throws `MHD_UNSUPPORTED_PHYSICS`**
for `mhd.enabled=true`. Nothing past this point runs for MHD. The intended
chain (not wired): `HydroAmrRuntime` kHydroUpdate → build
`MhdLevelPatchBinding`s from `SimulationState::magnetic_faces` + gas rows via
`loadMhdGasFromState` → `MhdLevelStepper::admitMemory` →
`estimateStableTimeStep` to the TimeCoordinator → `advance` →
`storeMhdGasToState` → thermo cache refresh → diagnostics → v26/v7 IO.

## 6. Specification reconciliation and ambiguities (recorded, not invented)
1. **UCT-HLL form:** REF fixes "UCT-HLL with consistent 1-D reduction" but not the
   transverse state construction. Chosen: Londrillo–Del Zanna 2-D HLL with
   PLM-MC transverse staggered fields and MC corner velocities in both PLM and
   PPM modes (exact 1-D reduction to the HLL face flux with PLM-MC states).
2. **Characteristic PPM:** limiting in the 7-wave primitive eigenbasis of the
   cell (left vectors by numerical inversion); passive metal by scalar PPM;
   singular basis or non-positive face state → counted PLM fallback.
3. **Conducting wall ghosts:** v_n odd, B_n even, B_t odd, plus exact E_t = 0
   on wall edges (zero mass/energy/tangential-stress flux, wall flux preserved).
4. **Hall-HLL whistler speed:** c_w = π|η_H|/Δ (grid Nyquist); the fan bound is
   therefore stricter than the diffusive-form bound by π and always wins.
5. **Expansion ownership:** existing hydro owns −H m and −H(2K+3p) through its
   comoving gravity source; MHD in-stage source includes −H M_c. Integration
   must select one owner per run (`apply_expansion_source`). **Open.**
6. **Comoving coefficient conversion for RKL2 stages:** cycle-midpoint a for all
   stages of one cycle (properties re-evaluated per stage).
7. **BZeroLimit layout:** defined as group `/BZeroLimit/{EtaOhmicCgs,EtaHallCgs,
   EtaAmbipolarCgs}` `[nrho][nT][ncomp][ngrain]` (REF names the slice only).
8. **Restart patch geometry:** REF §13.2 lists no origin/spacing datasets;
   `patch_origin_comoving` and `cell_width_comoving` added to the v26 group
   because geometry validation requires them.
9. **Edge-integral diagnostics for Hall work:** discrete J·E_H summed over edge
   families vanishes to O((k dx)²), not exactly; test bound is relative.

## 7. Remaining findings (priority)
* **P0** production dispatch absent (fail-closed); **P0** v26 writer/reader and
  v7 snapshot absent (fail-closed); **P0** compile viability unverified.
* **P1** no MPI halo/migration; no multi-level stepping/regrid wiring; expansion
  owner; gravity/cooling/source composition and Strang half-source hooks;
  `amr_temporal`/recursive magnetized mode rejected only by guard/config.
* **P1** no OpenMP parallel loops; no profile/counter export (MP12).
* **P2** tests for I02, I06, I08, N02, N08, oblique I01, 1,000-step G03, repeated
  plateau; GPU kernels (out of compulsory scope).

## 8. Next audit target
Compile and run, in order: `test_unit_mhd_riemann_reconstruction`,
`test_unit_mhd_discrete_complex`, `test_integration_mhd_ideal_reference`
(advection direction, multi-patch equivalence, axis rotation first), then
`test_integration_mhd_nonideal_reference`, `test_integration_mhd_amr_mpi_restarts`,
`test_unit_mhd_memory_governance`, and the MHD-off regression suite
(`test_config_parser`, `test_config_examples`, memory accounting, restart v25)
to confirm the guards and the conditional normalized text change nothing.
Source-audit focus: `MhdLevelStepper` reconciliation and ghost fill
(`src/mhd/mhd_level_stepper.cpp`), characteristic PPM eigenvectors
(`src/mhd/mhd_step_operator.cpp`), Balsara prolongation coefficients
(`src/mhd/mhd_amr.cpp`).

## 9. Interruption checkpoint
- All slices A–F have source deliverables as listed above. Campaign complete
  for this pass; nothing half-edited.
- Next agent: start at §8; do not remove the fail-closed guards until the
  dispatch, v26 writer/reader and v7 writer exist and their tests run.
