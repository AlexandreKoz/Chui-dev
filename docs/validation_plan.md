# Validation plan

## 2026-11-04 final-design closure — owner-run only

New or extended tests (all NOT RUN):

* `integration_distributed_element_diffusion_mpi_np{2,4}`: the decomposition
  reference is matched to 1e-12, and global lane totals are conserved to 1e-13;
* `integration_amr_family_colocation_mpi_np{2,4}`;
* `unit_amr_recursive_subcycling`: loopback distributed driver is bitwise equal
  to the local driver;
* `unit_chui_part3_remote_sources`: association partition invariance and
  energy and momentum conservation;
* `unit_workflow_metal_diffusion_topology`;
* `unit_stellar_population_feedback`: 200-way SSP subdivision.

Thresholds are round-off bounds for explicit convex updates and id-ordered
sums. None was widened. See
[`repair/chui_final_design_closure_20261104.md`](repair/chui_final_design_closure_20261104.md) §8.

## 2026-11-04 AMR admission/history checkpoint — owner-run only

New `unit_amr_history_atomicity` and extended `unit_amr_recursive_subcycling` /
`unit_amr_pending_flux_registers` cover preallocation refusal, injected C++
allocation failures, atomic multi-level history start/end publication and pending
register publication. History fields and rejection-state gas lanes are exact
bit targets; existing hydro/interpolation/conservation tolerances remain
numerical targets. All are **NOT RUN — USER-RESERVED**. No build, test,
benchmark, simulator, converter, sanitizer, static-analysis executable or CI
workflow was run. The standalone command selection and unimplemented MPI/science
families are in [the closure ledger](repair/chui_conformance_closure_20261104.md).

The [20261012 test regression closure](repair/treepm_deep_hotpath_test_regression_closure_20261012.md)
repairs invalid split-only/config fixtures and stale PM/rung expectations.
New controls require field-specific rejection, ten stable experiment decks,
actual PM retained-capacity reports, and always-active assertions under NDEBUG.
**NOT RUN — deferred to maintainer.** Follow its focused-to-full command order;
existing scientific references and tolerances retain their acceptance authority.

## Attempt 1 recovery gates — 20261012

Source regressions now cover actual composed runtime admission/drift grants,
multi-step CLI gravity progress, MPI grant agreement, periodic half ties and
multi-wrap frame certificates, independent long-double radial derivatives,
error-gated residual tables, scalar/AVX2 lanes, coincidence/self/tails,
fast/full accepted-work equality, admitted/denied/evicted governed caches,
reuse/refit generation/layout invalidation, evolved bridge unit factors and
production hierarchical optional-kernel/restart comparisons. Existing force,
Ewald, timestep-refinement and distributed assertions retain their thresholds.
All are **NOT RUN — deferred by user**. Execute the
[A–F manual sequence](repair/treepm_deep_hotpath_attempt1_closure_20261012.md)
before accepting performance/numerical claims. Adaptive 0.5, snapshot replay,
independent total periodic forces and 48 GiB fit each need separate evidence.

CosmoSim uses a four-level validation ladder: unit, integration, regression, convergence.

## Validation structure in-repo

- Unit: `tests/unit/`
- Integration: `tests/integration/`
- Validation ladder: `tests/validation/`
- Frozen tolerances: `validation/reference/validation_tolerances_v1.txt`

`validation/input_decks/` is currently documentation/reference material only. It is **not** the authoritative executable runtime path in this repository build. The real runnable path is the config-driven application entry point (`cosmosim_harness <config.param.txt>`) plus the tested release/runtime smoke configs.

## Required command sequence

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

For HDF5 snapshot/restart flows and the real runtime smoke path:

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure
```

For PM/TreePM FFTW-specific validation when FFTW is available:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```

## Gravity Phase 1 evidence classes (explicit references)

- **Open-boundary small-`N` tree check**
  - Method: direct summation, softened with the same `epsilon_comoving` as tree solver.
  - Claim: local tree operator correctness for non-periodic pair accumulation.
- **Periodic PM plane-wave mode check**
  - Method: analytic periodic spectral solution for selected Fourier mode(s).
  - Claim: PM operator/sign/phase correctness in periodic domain.
- **Periodic PM uniform-density cancellation**
  - Method: exact lattice plus deterministic parity-balanced micro-jitter (glass-like proxy) particle placement, periodic Poisson solve with mean subtraction.
  - Claim: near-zero force under uniform density.
- **Periodic TreePM consistency check (integration test)**
  - Method: minimum-image periodic direct reference.
  - Claim: seam, cutoff, split, and distributed protocol regression only.
  - Explicit limitation: the minimum-image reference is not Ewald exact and is
    not the periodic force-accuracy certification.
- **Independent periodic Ewald reference and TreePM accuracy gate**
  - Method: test-only double-precision real-plus-reciprocal Ewald sum with
    independent alpha/truncation convergence, homogeneous-background `k=0`
    convention, central-self removal, and rectangular boxes.
  - Claim: small-N periodic force accuracy for the explicitly named production
    profile. Current certification is FFTW-backed TSC with matched window
    deconvolution, `asmth_cells=1.25`, `rcut_cells=6.25`, quadrupoles, and
    COM-distance opening at the typed/public runtime default `theta=0.7`.
    Lower positive cutoffs remain compatibility inputs, not certified defaults.
  - Gate: `relative_L2 <= 1e-2`, `p99_normalized <= 5e-2`, net-force fraction
    `<= 1e-12`, and seam translation drift `<= 1e-11` on the covered fixtures.
  - Limitation: CIC is measured but diagnostic for this envelope; it is not
    silently accepted when it exceeds the TSC gate.
- **Zoom long-range correction + contamination diagnostics check (integration test)**
  - Method: synthetic high-resolution membership mask with focused PM correction enabled; compare high-res vs low-res target force deltas against no-zoom baseline and assert contamination counters.
  - Claim: explicit zoom force decomposition is wired, high-res membership is not implicit, and contamination failure modes remain visible.
- **Periodic TreePM force-error mapping (benchmark artifact)**
  - Method: fine spectral PM plus exact pairwise short-range residual periodic proxy reference.
  - Explicit limitation: not Ewald exact.
- **Split-kernel matching checks**
  - Method: Gaussian split SR/LR factor complementarity and TreePM PM-only/tree-only/split comparisons.
  - Claim: split coupling semantics are wired and monotone in expected direction.
- **Two-body orbit controlled run + energy drift monitor**
  - Method: leapfrog integration with tree accelerations and softened two-body energy accounting.
  - Claim: bounded long-run drift in a controlled setup.
- **Static halo radial-force profile**
  - Method: tree force on shell probes versus softened direct-summation reference.
  - Claim: radial-force profile agreement for a static mass distribution.
- **Regression-only checks (AMR reflux / SF mass budget)**
  - Claim: implementation stability against historical behavior, not standalone validation proof.

## Force-error mapping artifact requirement

Run:

```bash
./build/pm-hdf5-fftw-debug/bench_tree_pm_force_error_map
```

Expected deterministic artifact:

- `validation/artifacts/tree_pm_force_error_map.csv` (generated artifact; do not rely on a source-controlled copy)

The CSV is the required Phase 1 force-error map over PMGRID/ASMTH/RCUT and must be attached/referenced in review evidence.

## Validation invariants

- No silent comoving/proper or code/SI unit mixing.
- Mode and boundary policy must remain explicit and validated.
- Config unit coverage must reject an otherwise valid isolated/open deck when
  `numerics.treepm_enable_window_deconvolution` is `true` (including inherited
  periodic default `true`) and must accept the same deck when the key is
  explicitly `false`; no silent runtime override is allowed.
- The runtime application must load a real config and emit a real run directory outcome.
- Restart roundtrip must preserve state consistency and provenance.
- Snapshot schema compatibility and alias-read behavior remain tested.
- Deterministic/reproducible config normalization and hashes remain stable for unchanged inputs.

## Expected tests for feature work

- **Unit:** parser/math/layout or local invariants.
- **Integration:** real pipeline path through owning modules.
- **Regression:** preserve prior output metadata/behavior with tolerances.
- **Convergence:** when numerical order or resolution behavior is affected.

## Hydro classical CI guards

The default CPU validation ladder now includes:

- Executable: `test_validation_hydro_classics`
- CTest entry: `validation_hydro_classics`
- Labels: `validation;integration;hydro`

This executable runs small fixed-grid Cartesian versions of Sedov, Noh, Gresho
vortex, Kelvin-Helmholtz, and an Evrard-style collapse toy. The assertions are
deliberately CI-scale: every case checks finite state, positive density,
positive pressure/internal energy, and at least one physically meaningful
invariant or trend. These tests are not paper-grade reference-profile or
convergence campaigns.

The Evrard-style case uses a fixed analytic inward gravity source through the
hydro source-term interface. It is a bounded source-coupling and collapse-trend
guard until a self-gravitating hydro workflow can own a restartable Evrard
reference path with documented gravity/hydro ownership.

## AMR hydro CI guards

The default CPU integration ladder includes three AMR hydro validation guards:

- `integration_amr_hydro_shock_tube`: refines the patch adjacent to a Sod-like
  discontinuity, advances through production AMR hydro, and checks finite
  positive state, coarse-fine ghost fill, flux-register/reflux completion,
  conservation bounds, and shock/contact trend.
- `integration_amr_hydro_sedov`: refines a compact Sedov-style blast fixture,
  verifies child states are nonzero after refinement, advances the refined
  hierarchy, and checks finite positive state, approximate symmetry, outward
  shell motion, and total-energy bound.
- `integration_amr_synchronization_stress`: evolves, refines, evolves across a
  coarse-fine boundary with flux registers and reflux, derefines through
  conservative restriction, evolves again, and checks global conserved totals.

These tests are CI-scale synchronization and conservation sentinels for H3 AMR
hydro. They intentionally do not replace future high-resolution AMR convergence
decks, restart equivalence decks, or cross-code shock/blast comparisons.

Restart comparison blocker: these AMR hydro fixtures currently validate the
in-memory production patch/regrid/ghost/reflux path only. A short AMR
run/restart/run comparison remains a follow-up requirement once the HDF5-enabled
restart gate owns an AMR patch-hierarchy fixture with explicit patch, gas-cell
identity, and reflux-synchronization expectations.

## Validation reporting in PRs

Every PR should include:

- build preset(s),
- exact test commands,
- pass/fail status,
- any intentionally deferred validation and rationale.

## Runnable smoke/config examples

Use one of these config-driven paths for honest runtime smoke checks:

- `configs/release/release_smoke_zoom_in.param.txt`
- `configs/release/release_smoke_isolated_galaxy.param.txt`
- `configs/release/release_smoke_cosmo_cube.param.txt`

Or run the built-in integration smoke gate:

- `integration_runtime_app_smoke` (HDF5+FFTW-enabled builds; the production TreePM runtime fails closed when only the diagnostic naive-DFT backend is available)
- `integration_runtime_app_mpi_treepm_smoke_two_rank` (MPI+HDF5+FFTW-enabled builds; launches the real `cosmosim_harness` through `mpiexec` with `parallel.mpi_ranks_expected = 2`)



## Phase 3 validation status

- Phase 3 validation closure criteria are defined in `docs/treepm_phase3_contract.md`.
- This repository stage does not yet provide a dedicated `test_validation_phase3_*` executable/gate.
- Current Phase 2 MPI gravity validation is necessary baseline evidence, but not sufficient evidence for Phase 3 closure.

## Phase 2 distributed TreePM validation gate (MPI gravity gate)

- Distributed workflow snapshot/restart paths are rank-qualified (`..._rank###.hdf5`) so MPI restart validation does not rely on unsafe multi-rank writes to a single file.

The Phase 2 gate is now a dedicated MPI validation suite:

- Executable: `test_validation_phase2_mpi_gravity`
- Registration requirement: `COSMOSIM_ENABLE_FFTW=ON`. The diagnostic naive-DFT backend is not production-gravity acceptance evidence.
- CTest entries:
  - `validation_phase2_mpi_gravity_single_rank`
  - `validation_phase2_mpi_gravity_two_rank`
  - `validation_phase2_mpi_gravity_three_rank`
  - `validation_phase2_mpi_gravity_four_rank`
  - `validation_phase2_mpi_gravity_eight_rank`

### Numerical contracts enforced

- Distributed PM equivalence vs one-rank reference: `rel_L2 <= 1e-10`.
- Distributed full TreePM equivalence vs one-rank reference: `rel_L2 <= 5e-6` and `max_rel <= 5e-5`.
- Force-convention contract: the PM and complementary tree return one
  scale-free comoving kernel using `G_code` derived from configured units. The
  non-unit-scale regression requires force ratio `1` between otherwise
  identical `a=0.5` and `a=1` solves; cosmological response belongs to the
  particle/gas update layer rather than an extra kernel `a^2` factor.
- Gravity-timestep coordinate contract: for comoving softening, production
  converts `|A|` to comoving-coordinate acceleration `|A|/a^3` after step
  commit and evaluates `dt_grav=eta sqrt(a^3 epsilon_com/|A|)`.
- Cosmological gas-source contract: at non-unit `a`, unit coverage requires
  `rho A/a^2-Hm` momentum and
  `rho(u dot A)/a^2-H(2K+3P)` energy. A zero-gravity homogeneous-expansion
  case isolates the energy term. Fixed-grid and AMR workflow paths must both
  use the post-drift timeline's step-end `a,H_code`.
- Comoving gravity-timestep public API contract: unit coverage at `a=0.5`
  checks `computeComovingGravityTimeStep` against
  `eta sqrt(a^3 epsilon_com/|A|)` and rejects `a=0`. The production particle
  and gas criteria must use this helper; `computeGravityTimeStep` remains the
  coordinate-neutral primitive.
- Communication stress path: tiny tree exchange batches
  (`tree_exchange_batch_bytes=64`) plus explicit coordinator refresh/local-reuse
  calls, checked against the same TreePM thresholds. This does not enable
  production config cadence greater than one.
- Coordinator cache-safety contract: reuse requires matching
  force epoch, force-evaluation scale, `G_code`, split and rectangular-box
  geometry, assignment, boundary, PM decomposition mode, and deconvolution.
  Particle decomposition is deliberately not a fixed-slab PM-field invalidator.
  Missing/incompatible cache state on one rank must make explicit reuse fail
  coherently; divergent explicit refresh votes must fail before PM collectives.
- Observability contract: focused periodic/workflow tests inspect local
  source/target/tree occupancy, global empty-rank counts, hierarchy/peer
  traffic, PM solve/reuse, halo/local-slab dimensions, and tree
  build/multipole/opened-node counters. Counters are deterministic structural
  evidence, not a substitute for timing-based strong/weak scaling artifacts.
- PM-field event coherence: a refresh event emitted before cadence commit must
  carry the current `PmRefreshDirective` opportunity/version/build step/build
  scale factor, including the initial bootstrap, rather than the previous
  committed `PmSynchronizationState` values.
- Operational precision: end-to-end report coverage requires the first PM
  refresh to carry opportunity/version `1` and rejects `0.000000` formatting
  for nonzero `gravitational_constant_code` and
  `tree_relative_force_acceleration_floor`. Gravity diagnostic doubles use
  scientific `max_digits10` strings.
- Restart continuation contract in MPI mode: reference workflow restart write/read roundtrip must report `restart_roundtrip_ok=true`, preserve rank-qualified restart naming (`..._rank###.hdf5`) for multi-rank runs, preserve distributed ownership/slab metadata exactly across write/read, and pass a production direct-vs-resumed continuation test under the same topology. Rank-count-changing restart is an explicit negative test and must fail as unsupported.
- Production time-integration honesty contract: normalized config requires
  `hierarchical_max_rung=0`; resume rejects nonzero scheduler max/committed bins
  or nonzero pending bins. Mixed-rung KDK is not accepted until per-element
  kick/drift epochs exist.
- Gravity invalid-state trap contract:
  - NaN/Inf force or PM diagnostic norms must raise fatal gravity-state runtime errors.
  - Illegal gravity sync-state transitions (cadence/version/opportunity regressions) must raise fatal gravity-state runtime errors.
  - Distributed restart compatibility must reject gravity-specific metadata mismatches (cadence state and field-version/refresh coherence).
- Periodic TreePM geometry contract: typed config and direct coordinator tests
  must reject `r_cut >= min(Lx,Ly,Lz)/2`, while the certified DMO/Ewald fixture
  keeps the cutoff strictly inside that one-minimum-image bound.
- Distributed workflow honesty floor: the reference workflow MPI test must prove that the final runtime state is partitioned rather than silently replicated by checking (a) reduced local particle/cell counts equal the reported global counts and (b) the reduced particle-ID sum/xor matches the deterministic generated-IC set.
- Required MPI CI lane: `mpi-hdf5-fftw-debug` must select the registered two-rank periodic PM/TreePM workflow tests, uneven PM slab/potential routing coverage, two/three/four-rank gravity validation, distributed gas migration, hydro interface conservation, AMR boundary/reflux, and distributed restart continuation. Missing MPI/HDF5/FFTW-MPI dependencies or insufficient runner rank capacity are blockers, not silent skips.

### Phase 2 scaling artifacts

Run in MPI-enabled builds:

```bash
cmake --build --preset build-mpi-hdf5-fftw-debug --target generate_mpi_gravity_scaling_artifacts
```

Outputs (artifact files):

- `validation/artifacts/pm_only_scaling_np1.csv`
- `validation/artifacts/pm_only_scaling_np2.csv`
- `validation/artifacts/tree_only_scaling_np1.csv`
- `validation/artifacts/tree_only_scaling_np2.csv`

For final integration hard-gates and exact command bundle, see:

- `docs/treepm_phase2_closeout.md`

## 2026-07-13 gravity hardening gates

The current gravity acceptance bundle adds these deterministic gates:

- `validation_periodic_ewald_reference`: Ewald convergence, alpha independence,
  integer-box translation invariance, zero self force, symmetric cancellation,
  sign, `G` scaling, and invalid-input rejection.
- `validation_tree_pm_ewald_accuracy`: quasi-uniform, clustered, rectangular
  seam, Gaussian-transition, and cutoff-separation classifications with full
  error distributions and tree-work diagnostics. The release gate applies only
  to the named default TSC/deconvolution/quadrupole/COM profile; diagnostic
  records cross CIC/TSC and deconvolution off/on, monopole/quadrupole and all
  three MACs, and axis/diagonal separation classes without relaxing the gate.
  An asserted independent-target MAC matrix calibrates
  `theta=sqrt(alpha)=0.07071067811865475` for `alpha=0.005`, requires internal
  accept/open work with zero particle pairs, and keeps all three MAC errors
  within 25 percent while each meets fixed Ewald gates.
- `integration_tree_pm_coupling_periodic_mpi_{two,three,four}_rank`: periodic
  hierarchy geometry, empty/source-only/target-only ranks, independent targets,
  all-empty rounds, migration-to-empty epochs, repeated exchange epochs, and
  distributed agreement. Divergent layout/context metadata must be rejected on
  every rank at coordinated solve entry; periodic unwrap/local-tree and
  PM halo-cache commit, compact active/PM/zoom target preparation, worst-case
  traversal-stack reservation, and exchange-workspace/request preparation must
  remain failure-voted before the next collective.
- `integration_pm_slab_halo_exchange_mpi_{two,three,four}_rank`: FFTW-compatible
  fixed-block slab ownership, same-peer left/right orientation, zero-width
  slabs, and nonblocking halo completion.
- the periodic PM `Nx=2,np4` case: zero-width high ranks enter the complete
  FFTW-MPI plan/transform/collective solve with backend-safe dummy allocation.
- `validation_phase2_mpi_gravity_{single,two,three,four,eight}_rank`: identical
  physical fixtures across rank counts and the unchanged strict TreePM
  thresholds `relative_L2 <= 5e-6`, `max_relative <= 5e-5`.
- `validation_dmo_zeldovich_workflow_single_rank` plus MPI np2/np3/np4/np8 and
  `validation_dmo_zeldovich_workflow_rank_equivalence`: real config loader,
  production KDK/TreePM, HDF5 restart, stable ownership, periodic linear-mode
  amplitude/phase/coherence, mass/momentum/COM checks, direct-vs-resumed
  equality, and stable-ID np1/np2/np3/np4/np8 comparison. A manual cross-preset comparator
  separately checks a true non-MPI serial binary against MPI world size one.

The DMO displacement-growth gate is increment-based. With
`Delta D_expected=a_final/a_initial-1`, both the direct-mode growth and
`sqrt(P_final/P_initial)` must grow and satisfy
`abs((D_measured-1)-Delta D_expected)/Delta D_expected <= 0.075`. Their two
measured increments must agree within relative `1e-3`. The final velocity is
also compared with the exact zero-force Hubble-drag trajectory
`u_ballistic=u_initial a_initial/a_final`: response/ballistic RMS must lie in
`[0.01,0.03]`, its growing-mode projection must lie in `[0.01,0.03]`, and
alignment must be at least `0.99`. The exact initial `4^3` lattice is a named
periodic-Ewald force fixture, whose short-step kick predicts response near
`0.019` of ballistic RMS. This is an independently anchored integration
regression, not an analytic LCDM growth proof.

Command-backed results from the MPI+HDF5+FFTW build on 2026-07-13:

- the focused ownership/halo/periodic/Phase-2 np1--np4 matrix passed 12/12;
- the DMO producer/comparator matrix passed 5/5;
- the certified Ewald profile's worst observed `relative_L2` was
  `8.059667017e-3` and worst `p99` was `8.157551429e-3`, both on the exact
  cancellation-dominated `4^3` DMO initial lattice. Its net-force fraction was
  `5.985844233e-17`, translation drift `1.523007584e-13`, and
  Ewald/PM/tree/total force L2 was
  `10.16780756` / `8.296819318` / `1.789203432` / `10.08586460`. The separate
  rectangular seam fixture retained translation drift `3.031142721e-14`;
  default-setting diagnostic cutoff-classified relative L2 values
  were `6.209e-3` (axis),
  `8.072e-3` (diagonal), `6.152e-3` (immediately beyond cutoff on-axis), and
  `2.397e-3` (1.1 cutoff on the diagonal);
- DMO direct-mode amplitude grew from `0.0049999895833398901` to
  `0.0050596390031745634`, or `1.0119299088208957` versus expected
  `1.0119399775244009`. Maximum Fourier phase drift was
  `7.2112604976761041e-14` radians with coherence `1`; uninterrupted/resumed
  phase/state drift was zero. Response/ballistic RMS was
  `0.018731640935373602`, growing-mode projection `0.01873163714589127`, and
  alignment `0.99999979769618985`.

The DMO gate retains the independent direct particle-mode diagnostic and now
also runs a deterministic binned three-dimensional `P(k)` estimator. The
estimator deposits on a `12^3` CIC (historical 2026-07-13 baseline; current built-in-FFT fixture is `16^3` CIC) mesh, deconvolves the matched assignment
window, uses 32 linear-k bins to the three-dimensional mesh corner, excludes
DC, accounts for all 1727 remaining discrete modes, retains four empty bins,
and reports without subtracting the `V/N = 0.015625` Poisson level. Its
versioned JSON records normalization, code units, assignment/window/shot-noise
policies, bin edges, centers, powers, mode counts, and explicit empty-bin state.

The np1--np4 gate reported fundamental powers
`2.1315894140899176e-6` initially and `2.1827523417862502e-6` after evolution,
giving `sqrt(P_evolved/P_initial) = 1.011929959672577` versus expected
scale-factor growth `1.0119399775244009`. After normalizing the intentionally
rank-specific `world_size`, all four JSON artifacts had SHA-256
`53ee518b63d21cd28790415b8f076614497536bb6a62c95276b8ed6c14984bf3`.
The cross-preset serial/MPI1 comparison covered 64 stable IDs with periodic
position error `0`, maximum velocity error
`0`, and relative mass error `0`; their JSON artifacts
were byte-identical with SHA-256
`8acbf3d6826250ca4fbabb0761511ff1178cbb7c20472cb2d8c8073dd16d355c`.
This closes the missing-artifact blocker for the controlled gate, but the
retained direct-DFT oracle and tiny linear fixture are not by themselves a
publishable/high-dynamic-range DMO campaign; the production estimator is now
FFT-backed and must still earn large-mesh scaling and aliasing evidence. The
estimator's normalization, window, bin, unit, and empty-bin contract is
maintained in `docs/power_spectrum_diagnostics.md`.

## Phase 3 research-grade evidence campaign

For Phase 3 gravity maturity evidence generation and separation of correctness/force-accuracy/scaling classes, use:

- Observable contract: `validation/campaign/phase3_gravity_observables.json`
- Campaign command bundle: `docs/validation_phase3_campaign.md`
- Provenance collector: `scripts/validation/collect_phase3_evidence.py`

The campaign requires explicit cosmological and zoom runtime runs plus force/time-integration/scaling artifacts. Existing `np1/np2` scaling outputs remain baseline performance evidence only and are not standalone certification.


## MPI AMR directed-exchange validation update

The MPI feature lane must include the real MPI-launched AMR directed-exchange
registration `integration_reference_workflow_distributed_amr_mpi_two_rank` and
must not count the serial `integration_amr_distributed_remote_patch_boundary`
fixture as MPI coverage. The serial test remains useful as an AMR data-structure
and orchestrator guard only.

The directed-exchange test exercises patch/cell payload exchange and owner-routed
flux-register exchange with asymmetric ranks and nonzero mass/momentum/energy
payload fields. The MPI+HDF5 lane also registers
`integration_distributed_amr_patch_migration_mpi_two_rank`, which performs a
genuine cross-rank AMR patch migration using the production bounded byte
transport and migration wire. The migrated payload includes pending reflux and
temporal-boundary history; the destination commit remaps them by stable gas/patch
identity, then the test checkpoints/reloads the migrated state and performs a
post-restart migration-commit continuation check.

These are CI-scale ownership/restart contracts, not a publication-grade AMR
scaling or convergence result. Full `ReferenceWorkflowRunner` direct-vs-resumed
comparison of all five conserved totals across a runtime-triggered AMR regrid is
still a stronger future validation target.

## Campaign B IC-ingestion acceptance repair gate (2026-07-24)

Serial/HDF5 acceptance requires:

- incomplete direct-bridge conventions rejected before HDF5 access;
- manifest schema v3 dimensional contract and version-2 rejection;
- direct, manifest-driven, and canonical converter numerical equivalence for
  coordinates, velocity, mass, density, internal energy, stellar initial mass,
  BH mass, and BH accretion rate;
- malformed attribute-shape, floating/signed-ID, wrong-rank, ambiguous-alias,
  and unknown-family fixtures under normal plus ASan/UBSan execution;
- single-file/multifile and mixed-species sidecar equivalence;
- converter source-manifest hash/schema rejection and atomic output.

MPI acceptance is registered but remains provisional until executed with the
`mpi-hdf5-fftw-debug` preset. The matrix includes one/two/four-rank direct,
canonical, and supplied-manifest ingestion; mixed Gas/DM/Star/PartType5-BH
input; alternate h/a/frame/velocity conversion; PartType2-to-tracer and
PartType3-to-star mappings; more ranks than files; exact duplicate rejection;
deliberate lost/duplicated routes; increasing-global-count bounded-staging and
local-authoritative-capacity checks; and two-rank one-sided fault injection for
owner serialization, send layout, payload validation, deserialization, sidecar
append, final state, and source/final reconciliation. Fault tests have 30-second
timeouts; increasing those timeouts is not an acceptance repair.

Required commands:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure
ctest --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
  -R 'integration_distributed_ic_reader'
```

## Post-optimization source repairs (2026-10-28)

The distributed IC protocol-v3 successful routing budget is 19 measured MPI
collective calls per batch (formerly documented as 20); the exact per-rank
accounting assertion remains enabled. The installed-package consumer must
execute its own MPI/Parallel-HDF5 capability compile/link probe rather than
preinitializing its result to false. Deterministic star formation returns no
birth for a positive expected mass below the configured minimum particle mass
without invoking the stochastic target-mass feasibility search; the
`unit_star_formation` regression checks this boundary. The multi-rank gas
force fixture must execute its rank-local-serial reference solves on *every*
rank, even those with no active target. The distributed PM potential fixture
checks globally reduced routed requests, retaining per-rank upper bounds and
numerical field comparisons.

The numerical gravity invariant, hierarchical trajectory, group-walk
accounting, SFC failure-injection, and general distributed collective failures
remain separate acceptance gates and must not be waived by these fixes.

## Adaptive bound–Jeans star-formation validation

Acceptance coverage includes physical/comoving conversion, gradient norm, virial and Jeans transitions, converging/static/expanding flow, free-fall/compression times, exact exponential depletion, source timestep limits, stable RNG golden vectors, row/rank-offset invariance, unbiased stochastic ensembles, multi-particle count consistency, conservative mass/momentum/metals, internal-energy sink accounting, invalid-state rejection, disabled no-op behavior, legacy compatibility, actual `SourceRuntime` gradient/metal use, snapshot round-trip, restart schema coverage, AMR metal transport, feedback metal return, and distributed ownership ledgers.

The affordable repository tests establish implementation correctness and deterministic contracts. They do not constitute observational calibration of an isolated disk or cosmological galaxy population. Resolution studies, long-duration feedback-regulated disks, and external cross-code comparisons remain release-validation work.

### Adaptive bound--Jeans star-formation acceptance

The campaign adds `unit_star_formation`, `integration_star_formation_box`, and
`integration_star_formation_source_runtime` for the analytic model, conservative
batch mutation, real production hydro gradients/metallicity, parent gas-particle
gravity-mirror refresh, snapshot/restart fields, and disabled/legacy behavior.
The MPI-enabled HDF5 configuration additionally registers:

- `integration_star_formation_mpi_two_rank` (owner-only spawning, including a rank
  with zero eligible local cells, against a one-rank physical reference);
- `integration_star_formation_mpi_repartition_two_rank` (birth-key, mass, and ID
  invariance under ownership exchange);
- `integration_star_formation_restart_mpi_two_rank` (direct versus checkpointed
  continuation with exact birth records and gas mass).

These MPI cases are source- and syntax-validated in environments without MPI, but
acceptance may be claimed only after they execute under `mpiexec -n 2` in an
MPI+HDF5 build.


## Star-formation completion validation

`validation_star_formation_models` validates the deterministic EOS table, threshold continuity, monotonic pressure, positive closure signal speed, interpolation, adaptive resolution response, exact timestep subdivision, and stochastic ensemble expectation. Running the executable with an output directory generates machine-readable files under `validation/star_formation/`. Production-runtime integration covers explicit profile selection, authoritative hydro inputs, covered-coarse rejection, leaf-only spawning, dense-row/patch reorder invariance, conservative birth, restart, and HDF5 fields.

The HDF5 snapshot acceptance path also includes a CHUI-native unit-consistency regression with
`length_unit = kpc` and a 0.05 Mpc box. A valid 4.0965 kpc coordinate must pass even though its
raw numeric value exceeds `0.05`, while a 50.1 kpc coordinate must still fail the same 50 kpc
box bound. This protects both the direct HDF5 validator and normal snapshot readback from
mixing Mpc-valued CHUI box metadata with code-unit coordinates.

The retained data are analytic and controlled-model validation. `validation_metadata.txt` explicitly records that no long-duration isolated-disk or observational calibration was performed. MPI acceptance requires the registered two-rank tests to execute in an MPI+HDF5 environment.

## MPI P1 count/collective-safety regression gate (2026-08-30)

The distributed-memory correctness floor now includes two complementary layers.

The CPU unit test `test_unit_parallel_distributed_memory` exercises the bounded
MPI round planner with synthetic limits rather than giant allocations. It
covers zero/all-empty traffic, one participant, asymmetric peers, exact-limit
and limit-plus-one cases, an aggregate that exceeds the simulated classic-MPI
domain while individual peers fit, multi-round exact coverage, deterministic
logical offsets, and checked multiplication/prefix overflow.

When an MPI runtime is available, `test_integration_mpi_collective_safety` must
run at two ranks and three ranks. The test forces tiny bounded rounds, injects a
single-rank failure after collective control metadata, verifies coherent
failure without entering an unmatched payload phase, and exercises a sequential
ghost-peer schedule where a post-metadata rejection must not strand a later
peer. The three-rank case also uses a deliberately sparse payload topology in
which ranks 0 and 1 exchange ghosts while rank 2 has no payload neighbor. Both
normal success and an injected rank-0 payload-preparation failure must leave
rank 2 participating in the final world agreement; a post-operation barrier is
the liveness sentinel.

`integration_distributed_gas_cell_migration_mpi_two_rank` sets the test-only
`COSMOSIM_MPI_TEST_TRANSPORT_LIMIT_BYTES` ceiling to 8 bytes. The serialized
migration packets necessarily exceed that ceiling, forcing the real production
bounded `MPI_Alltoallv` transport through multiple rounds while the existing
restart, stable gas-cell identity, global cell-count, and metal-mass assertions
verify that transport chunking does not alter authoritative state.

These tests are specifically designed to prove protocol correctness without
`INT_MAX`-sized buffers or deliberate OOM. Runtime MPI acceptance is evidence
separate from source closure; an environment without `MPI_CXX` keeps this gate
environment-blocked rather than converting the source repair into a different
algorithm or weakening the tests.

## M2D final admission closure gates

Run `unit_analysis_diagnostics`, `unit_memory_governor`,
`unit_runtime_module_registry`, and `integration_reference_workflow` on the
CPU and HDF5 Debug presets. The diagnostic unit test covers checked 256^3/512^3
mesh models, below-Red insufficient headroom, exact-boundary recovery, lease
release and scientific equivalence. The workflow test covers coalesced optional
cadence. Follow with the full inventories, source-package completeness and
MPI+HDF5+FFTW rank/fault-injection/science matrix on a provisioned host. Do not
promote CPU algorithm checks into distributed runtime certification. The
remaining owner-specific peak contracts and production full-process envelope
must be measured before closing M2D or enabling nonserial full-physics work.

## M2D major-task contract acceptance (2026-09-07)

The focused registry test covers checked retained-capacity arithmetic, dynamic
estimate evaluation, complete/partial transitions, dispatcher-owned RAII
release, owner-managed preflight without double charging, incomplete-owner
batch freedom, and qualified task lookup. The end-to-end reference workflow
continues to require a deterministic hard-memory rejection for an impossible
one-byte process ceiling and preserves the scientific stage/diagnostic checks.

Required remaining acceptance: complete source/AMR/output owner models,
CPU/HDF5/source-package inventories, MPI+HDF5+FFTW np2/np3/np4 skew/empty/
infeasible and fault-injection tests, migration coexistence and restart/science
equivalence, measured rank max/mean RSS/PSS and work imbalance, and legal
batch-size and Release rank/thread comparisons. Use the authoritative presets
in `docs/build_instructions.md`; do not count unavailable dependencies as
passing or change numerical tolerances for memory work.


## M2D final coding follow-up gates (2026-09-08)

`unit_particle_id_registry` exercises deterministic collision rehashing,
invalid initial IDs, retained/released leases, empty batches, and tight-budget
rejection. `unit_sidecar_row_lookup` exercises unsorted rows, duplicates,
missing indices, exact headroom and overflow. `unit_memory_governor` now also
checks bounded arena exhaustion and lease release. The source and reference
workflow integrations preserve the accepted scientific order.

MPI-only source-ID branches require a real MPI build and np2/np3/np4 execution;
a declaration-stub syntax check is not a runtime pass. The remaining AMR
geometry/ghost and I/O metadata allocation ownership must be completed before
claiming a whole-task certificate. Full CPU/HDF5 inventories, extracted-source
completeness, distributed fault/restart/science equivalence, and measured
whole-process/topology evidence remain separate acceptance gates.

### DMO first-light closure regressions (2026-09-16)

Focused CPU/HDF5 regressions cover: flat cosmology normalization and
`H(a=1)=H0`; explicit standard GADGET/AREPO `sqrt(a)` stored-velocity semantics
at `a=0.25` and `a=0.5` with cross-check against snapshot conversion; DMO
restart mutation sensitivity for particle drift-history state; and PM
spectral-operator cache reuse when only cosmological scale factor changes.
The PM scale factor remains in distributed collective consensus.

### DMO first-light runtime/I/O recovery regressions (2026-09-17)

First-light readiness now distinguishes kernel validation from the config-driven
production path. `integration_dmo_config_driven_firstlight_smoke_single_rank`
uses a tiny DMO state, `hierarchical_max_rung=0`, scale-factor authority, and
`runner.run()` with no `dt_time_code` override; it must require multiple physical
steps, remain below the configured safety step cap, and reach the FLRW-derived
endpoint. HDF5 acceptance includes a monofonIC-like fixture with uint64 totals,
rank-1 extent-one logical scalars, zero-based IDs, periodic-equivalent negative/
upper-edge coordinates, and fixed DM mass. Snapshot HDF5 roundtrip acceptance
uses two logical sets in one `snapshots/` directory and rejects missing members
and stale/misnamed completion markers. MPI variants remain runtime evidence, not
a pass, until executed with a real MPI C++ toolchain.

## 2026-09-21 — 64^3 three-rank first-light production-closure gate

The first-light closure adds acceptance coverage for the external-IC to analysis-snapshot path without changing TreePM or restart numerics. Required focused evidence is:

- `unit_ic_reader`: zero-containing external IDs map with one source-wide `+1` policy independent of row order; shuffled zero-based IDs are accepted; duplicate IDs and `0 + UINT64_MAX` overflow cases are rejected; the compensated mass accumulator is stressed at 64^3 = 262,144 equal-mass elements and still rejects one-particle loss/duplication without relaxing tolerance.
- `integration_snapshot_hdf5_roundtrip`: ordinary external HDF5 readback accepts shuffled zero-based IDs and the single-file completion contract remains backward compatible. The production-validator mode is also exercised with a zero-byte global-ID budget to prove the post-publication structural pass does not require an O(N_global) root ID buffer.
- `integration_snapshot_set_mpi_{2,3,4,8}_rank`: on a Parallel-HDF5 build, all ranks write one `snap_007.hdf5` through non-overlapping global hyperslabs. The fixture uses uneven ownership and an empty final rank, then verifies global counts, direct readback, persistent IDs, and the one-member completion marker. A serial-HDF5 MPI build retains the explicit sharded compatibility branch only.
- `integration_runtime_app_smoke`: an on-cadence endpoint produces exactly one committed science file. The same real harness is rerun with a cadence that cannot fire on step one and must still publish exactly one final endpoint snapshot.
- `integration_runtime_app_mpi_snapshot_smoke_{2,3}_rank`: narrow Parallel-HDF5 runtime gate with diagnostics/restarts explicitly disabled; uses off-cadence endpoint forcing and requires exactly one `snap_001.hdf5`, one completion marker, no shards, no `.partial`, and one rank-zero console stream.
- `integration_runtime_app_mpi_restart_smoke_two_rank`: restart-specific runtime gate with diagnostics disabled and rank-qualified restart artifacts.
- `integration_diagnostics_mpi_global_{2,3,4}_rank` plus `integration_runtime_app_mpi_diagnostics_smoke_two_rank`: analytically verify global counts/booleans/angular momentum/SFR/quicklooks/memory semantics and single-publisher shared files; MPI power spectrum must be explicitly deferred rather than averaged from local spectra.
- `integration_runtime_app_mpi_treepm_smoke_two_rank`: combined gate retained after the narrow snapshot/restart/diagnostics surfaces pass.
- `integration_chui_launcher_mpi_runtime_two_rank`: dependency-gated real `./chui run ... --preset` execution verifies authoritative rank inference, real MPI launch, child status, rank-zero console ownership, and science output.
- existing restart roundtrip/equivalence tests remain required because science snapshot topology is not allowed to alter checkpoint schema or same-topology continuation semantics.

Dependency-complete MPI/Parallel-HDF5/FFTW execution is a mandatory production qualification gate. A host that lacks those development libraries may establish source/configure/serial-HDF5 evidence but must report the MPI matrix as blocked rather than silently substituting the sharded backend.

## 2026-09-22 — Parallel-HDF5 smoke defect regressions

The Q-01/Q-02/Q-03 closure adds narrow regressions rather than relying on the real 64^3 smoke as the first detector. Configure qualification poisons the cached MPIO capability result and requires a real compile/link probe to execute; serial HDF5 is rejected when the production parallel requirement is enabled where a serial wrapper is available. Distributed IC acceptance runs `zero_ids`, `shuffled_zero_ids`, and zero-based-versus-pre-normalized equivalence fixtures at 1/2/3/4/8 ranks, retains duplicate rejection, and adds the `0 + UINT64_MAX` overflow conflict. The zero-ID cases assert the typed `IcExternalIdMapping`, the exact normalized global ID ranges, nonzero IDs, distributed ownership, and physical-state equivalence.

The narrow Parallel-HDF5 runtime snapshot smoke at 2 and 3 ranks begins with a nonexistent `snapshots/` directory, requires first publication to create it, reruns the committed namespace to prove overwrite refusal, and injects a regular file at the `snapshots` path to prove collective namespace preflight fails before HDF5 creation. Environment-specific unwritable-directory behavior remains suitable for dependency-host qualification because root CI cannot portably express POSIX permission denial.

## P6 follow-up qualification

The P6 source campaign deliberately executed no configure/build/compiler,
tests, MPI runs, benchmarks or scientific validation. Existing build/test floors
remain required for the second campaign. The exact staged matrix and evidence
requirements are in [p6_validation_campaign_plan_20261006.md](repair/p6_validation_campaign_plan_20261006.md).
Optional modes must not become defaults until those gates are accepted.

## TreePM unified recovery qualification (2026-10-10)

New source coverage, **not executed**:

- `unit_tree_pm_split_kernel`: independent long-double Gaussian coefficient,
  zero/tiny/extreme scales, full Hermite interval/edges, multiple softenings,
  invalid inputs, squared cutoff decisions and full/fast MAC decisions.
- `integration_tree_pm_coupling_periodic`: production direct-pair interaction
  counts at inclusive cutoff, rectangular/multi-wrap/half-box fixtures,
  analytic/lookup and worker reproducibility, retained table accounting.
- `validation_tree_pm_recovery_ewald_accuracy`: strict/adaptive and analytic/
  lookup with the existing Ewald accuracy budgets and finite softening.
  Adaptive reference-history cases use explicitly converged Ewald magnitudes,
  not simulated historical KDK history.
- `integration_hierarchical_timestep_regression`: displacement law, production
  core block dispatcher free expansion/momentum, sparse activations, source
  epochs, endpoint safety and malformed drift identity rejection.
- `integration_hierarchical_treepm_workflow` (plus `_np2`, `_np4`): full production
  workflow global-reference refinement, replay, PM endpoints, restart and mass.
- `validation_hierarchical_dmo_growth_single_rank`, `_mpi_np2`, `_mpi_np4`, and
  `validation_hierarchical_dmo_rank_equivalence`: existing Zel'dovich fixture,
  growth/power refinement at its six fundamental modes and stable-ID physical
  rank comparisons using the existing equivalence tolerance.

The frozen-source sweep produces a new direct short reference for each physical
split. Its optional small-N Ewald comparison requires increasing image/mode
limits to demonstrate convergence. An external total reference must match
particle IDs, masses, cosmological state and units. Timing alone never qualifies
mesh/leaf choices. No speedup, force budget, or hierarchical scale range is newly
certified by writing these tests. Exact commands and remaining gates are in
[the campaign report](repair/treepm_unified_performance_recovery_20261010.md).

## P1–P5 follow-up reliability gates — 20261010

`unit_tree_pm_split_kernel` now asserts the exact finite-softening coincident
residual limit, approach to zero, immutable lookup consistency and invalid-input
contracts. Periodic coupling checks reject NaN/Inf comparisons and explicitly
verify zero coincident force even when the radial coefficient overflows. Existing
seam, rectangular, half-box, cutoff, mixed-softening and distributed fixtures keep
their accuracy bounds.

`validation_tree_pm_recovery_ewald_accuracy` passes the mandatory executable flag
`--recovery-qualification`; clearing an environment variable cannot reduce it to
a baseline-only run. Strict/adaptive x analytic/lookup profiles must satisfy the
existing independent Ewald accuracy targets, exercise direct pairs and separately
accept internal multipoles. Its printed experimental qualification gates do not
certify adaptive momentum conservation or promote opt-in runtime capabilities.
The certified baseline conservation gate and finite-softening cases remain.

Production hierarchical workflow tests add verified default-key restart migration
and incompatible/insufficient/open-state refusals. Core hierarchical dispatcher
fault injection checks missing/stale endpoint PM events before cadence commit.
Source generation, all-source drift, active sets, endpoint kicks, replay and
restart/refinement assertions remain active. Full nonlinear convergence and
large-snapshot performance remain manual gates.

Official `test_` CMake targets already compile with `-UNDEBUG` or `/UNDEBUG`,
guarded by `unit_assertions_enabled`. Selected gravity/config/hierarchical/HDF5/
restart sources additionally use `tests/support/test_require.hpp`, independent
of build flags. No production NDEBUG policy or scientific tolerance was weakened.
Explicit CTest timeout/working-directory properties are preserved; numeric and
`_npN` MPI cases now receive rank-count processor accounting.

The Python reader matrix uses production-shaped v2/v3 closure fixtures and named
negative categories (completion, member identity/counts, units, masses, lanes,
IDs and checkpoint refusal). Missing Python packages are an explicit skip.
All new/changed gates are **NOT RUN — deferred by user**. See the
[repair report and manual matrix](repair/treepm_p1p5_qualification_and_test_reliability_repair_20261010.md).

## TreePM fast-settings qualification (October 2026)

Force level (done, synthetic clustered 64^3 box): total-force error quantiles
of every adaptive estimator / cutoff-node / alpha / PM-grid combination
against the strict reference on the same mesh (`docs/tree_pm_coupling.md`).
Acceptance used for the fast deck: rms <= 3e-4, p99 <= 1e-3.

Still required before calling the fast deck or hierarchical `rms_split_scale`
production-qualified:

1. The same force comparison on `snap_1698` (real clustering, real softening).
2. Matched-scale-factor trajectory comparison (positions/velocities by ID) of
   the fast deck against the strict reference deck over the snap_1698 -> z=0
   continuation, plus P(k) and halo mass function at z=0.
3. Hierarchical: `rms_split_scale` with `hierarchical_max_rung` 4 and 6 versus
   global KDK at matched scale factors (same diagnostics), restart equivalence
   and the existing `validation_hierarchical_dmo_growth_*` gates.

## AREPO-grade campaign validation additions (2026-10-26; authored, NOT RUN)

Record: [`repair/arepo_grade_campaign_20261026.md`](repair/arepo_grade_campaign_20261026.md).
All thresholds below are **provisional** (scientifically motivated, not
qualified); agreement with AREPO is a consistency check, not a proof of
correctness — both codes carry force, time-integration and discretization
error.

| Area | Owner-executable evidence | Provisional gate |
| --- | --- | --- |
| MPI routing admission | `integration_pm_routing_admission_mpi_{two,four}_rank` | every rank fails together; no hang (watchdog 120 s / case) |
| In-place FFT layout | `integration_pm_in_place_fft_layout[_mpi_{two,four}_rank]` | Fourier-mode response ≤ 1e-9 rel; distributed = serial ≤ 1e-9 rel |
| Remote group walk | `integration_tree_pm_remote_group_walk_mpi_{two,four}_rank` | grouped = per-target ≤ 1e-9 rel; identical in-cutoff pair counts |
| Fine-tick caches | `integration_tree_pm_fine_tick_caches` | retained = eager ≤ 1e-12 rel; uniform scalar = explicit lane (bitwise) |
| Gravity invariants | `validation_gravity_invariant_suite` | self-force ≤ 1e-2 one-cell force; Σ m a ≤ 1e-6 Σ m\|a\|; axis covariance ≤ 1e-9; two-body vs Ewald non-increasing with mesh and ≤ 2 % at 64³; FD ≤ 5 % |
| Rollback force endpoint | `unit_gravity_force_endpoint` | row-aligned restore bitwise equal to the stable-ID restart remap of a permuted checkpoint; repeated restores; zero heap allocations per restore; refusals (count, layout generation, identity certificate, lane extents, no endpoint) write nothing; restart-remap refusals (missing/duplicate/zero ID, non-finite, extents) write nothing |
| Particle mass generation | `unit_particle_mass_generation` | advanced by resize, reorder, index bumps and notified in-place writes; unchanged by drifts; cache-key equivalence |
| Multirate PM (driver) | `test_hierarchical_multirate_convergence` (CPU-only) | observed order ∈ [1.8, 2.2] over 3 halvings for PM-only (= leapfrog at T), tree-only (= reference), combined uniform/variable blocks (+ splitting error), adaptive rungs, split-force circular two-body orbit (momentum and angular momentum to roundoff); EdS cosmology: exact canonical momentum under constant split accelerations for any tiling/rungs/cycle (detects duplicated/missing PM half-kicks), second-order positions, cosmological oscillator self-convergence |
| Multirate PM (workflow) | `test_hierarchical_treepm_workflow` (multirate block, np 1/2/4) | block invariants (quantum ≤ criterion, exact tiling, T ≤ PM displacement bound, exact t_end, soft stop); segmentation invariance (4+4, 2×4, 1×8, 3+rest, 5+rest) bitwise with equal schedules and PM refresh totals; forced rollback bitwise (roundoff with tree reuse/refit); genuine displacement rejection (calibrated, deterministic replay, rejected work excluded from trajectory diagnostics); monotone convergence and splitting error vs exact; memory fallback (refusal on all/one/several ranks, allocation failure, first-only, `fail` policy) identical on every rank with no duplicated PM kick and no retained snapshot storage; consecutive rejections bitwise with accepted cadence history equal to the uninterrupted run; rank-local rollback failure fails every rank; rollback with tree reuse/refit/group walk (roundoff, deterministic replay); skewed/empty-rank ownership |
| Checkpoint cadence | `unit_output_cadence` | no dropped/duplicated requests; one request per crossing step |
| Matched-IC DMO vs AREPO | `tools/compare_dmo_arepo.py --enforce` | mass rel. diff ≤ 1e-12; P(k) ratio within 2 % and r(k) ≥ 0.99 for k < k_Ny/4; FoF mass-function bins (≥ 20 halos) within 10 %; median displacement ≤ 0.1 mean spacing |

Recommended owner sequence for the AREPO comparison: identical ICs
(`tools/convert_ic.*`), matched softening, `treepm_asmth_cells` /
`treepm_rcut_cells` comparable to AREPO's `ASMTH` / `RCUT`, matched PM mesh,
`gravity_timestep_eta` translated as documented (`eta_CHUI = sqrt(2
eta_AREPO)`), snapshots at the same scale factors (z = 9, 3, 1, 0), first with
the exact reference policies, then with each experimental policy toggled
separately (`multirate_outer_kick`, `finite_difference_4`, group walk,
topology refit) so differences are attributable.

## CHUI Part 1 human-run acceptance (hydro suite, F1/F3/F4/F7, memory) — 20261103

All items below are **NOT RUN — USER-RESERVED**. Unit sources:
`test_hydro_muscl_hancock_predictor`, `test_hydro_riemann_passive_scalars`,
`test_hydro_ppm`, `test_hydro_ctu`, `test_hydro_dual_energy`,
`test_gas_thermodynamic_refresh`, `test_effective_ism_energy_ledger`,
`test_node_memory_admission`, `test_stage_transaction`,
`test_config_hydro_numerics` (ctest names `unit_<name>`).

Additional campaign checks to run on top of the unit tests:
1. Smooth-wave and oblique advection convergence for PLM/PPM x MH/CTU at
   N = 16, 32, 64, 128 (`bench_hydro_method_matrix` prints the L1 error);
   expect ~2nd order for PLM, >= 2nd for PPM on smooth data, and CTU error
   below MH at equal Courant for oblique flow.
2. Sod, strong-shock (pressure ratio 1e5), Sedov and contact problems with
   each reconstruction; check positivity fallback counts and absence of
   negative pressure.
3. High-Mach cold advection (Mach 1e2-1e4) over >= 100 steps with
   `hydro_dual_energy = off` vs `thermal_fraction`; decide from the measured
   pressure error and `dual_energy_sync_delta` whether a persistent entropy
   lane is needed (Part 2/3 decision).
4. a != 1 consistency: identical results for (dt, a) and (dt/2, 1) in an
   isolated (non-expanding) box; comoving Zel'dovich pancake for the combined
   1/a and expansion-drag path.
5. Feedback heating x100: verify the next accepted interval drops by ~10 in
   a production `isolated_galaxy` run (`sources.thermodynamic_refresh` event,
   `max_sound_speed_increase_ratio`).
6. AMR: synchronized two-level patch run with PPM/CTU; inspect
   `stencil_degradation_count`/`ctu_transverse_degraded_states` and reflux
   conservation (mass, momentum, energy and metals) to round-off.
7. Mode smoke for `cosmo_cube`, `zoom_in`, `isolated_galaxy`,
   `isolated_cluster` with defaults (behavior must match pre-Part-1 apart from
   the documented F1/F3/periodic-face fixes).
8. Memory: run with `node_memory_admission = shared_equal` and `shared_pool`
   and 1, 2, 4 ranks per node; confirm `memory.node_admission` shares, that the
   node pool bounds the sum, and that refusals are reported, not crashes.
9. Periodic AMR: `cosmo_cube` with a box-spanning root patch, uniform-velocity
   metal contact crossing the box boundary for one box crossing (exact
   translation check, global mass/momentum/energy/metal conservation); confirm
   that a decomposed non-spanning root grid stops with the documented refusal.
10. `unit_amr_periodic_boundary` (new).
11. Stabilization: `unit_hydro_cfl_stability_rule`, `unit_amr_reflux_thermodynamics`,
    `integration_amr_production_hydro_integration` (refinement caches). Run a 3D
    diagonal advection with `muscl_hancock` at the new summed CFL for >= 1000
    steps (must stay bounded) and compare step counts with `ctu`.
12. MPI shared pool: on an MPI build with >= 2 ranks per node, check the
    `memory.node_admission` event (`node_pool_established`,
    `node_pool_mpi_shared_window`, `node_pool_fallback_reason`) and a forced
    contention test where one rank reserves most of the node budget.
