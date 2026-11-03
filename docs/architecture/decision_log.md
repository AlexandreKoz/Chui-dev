# Architecture decision log

## 2026-07-13 — ADR-REPAIR-GRAVITY-018: Cosmological gas-source normalization and non-degenerate DMO growth gates

### Status
Accepted (scientific-correctness bug fix)

### Context

The gravity repair makes PM and the complementary tree return a scale-free
comoving kernel `A`. The collisionless KDK path consumed that convention, but
the gas source still treated `A` as an already scaled peculiar acceleration,
recomputed a background Hubble rate outside the integrator's code-unit
timeline, and used an incomplete total-energy expansion source. The controlled
DMO gate also used absolute response floors that were weak relative to the
short expected growth increment and its seeded Zel'dovich velocity.

### Decision

- Treat `HydroSourceContext::gravity_accel_*_peculiar` as
  compatibility-named views of scale-free `A`. Apply `A/a^2` exactly once in
  gas momentum and gravity work.
- Hydro runs post-drift before `IntegratorState` commits the step. Both
  fixed-grid and AMR callbacks therefore consume
  `StepContext.timeline_step.scale_factor_end` and `hubble_end_code`, already
  in inverse code time.
- For `m=rho u`, `K=rho |u|^2/2`, and total comoving gas energy, apply
  `S_m=rho A/a^2-Hm` and
  `S_E=rho(u dot A)/a^2-H(2K+3P)`. This pressure form is valid for every
  adiabatic index supported by `HydroCoreSolver`.
- Force lanes continue to store scale-free `A`. When the gravity timestep
  criterion combines that force with comoving softening, use the public
  `ComovingGravityTimeStepInput`/`computeComovingGravityTimeStep` boundary. It
  validates finite positive `a`, converts to comoving-coordinate acceleration
  `|A|/a^3`, and after step commit evaluates
  `dt_grav=eta sqrt(a^3 epsilon_com/|A|)` with the committed scale factor. The
  existing `computeGravityTimeStep` remains coordinate-neutral.
- Let `Delta D=a_final/a_initial-1`. Require the direct displacement-mode and
  `sqrt(P_final/P_initial)` growth increments to be positive and to satisfy
  `abs((D_measured-1)-Delta D)/Delta D <= 0.075`.
- Compare final velocity with the exact homogeneous Hubble-drag trajectory
  `u_ballistic=u_initial a_initial/a_final`. Require response/ballistic RMS in
  `[0.01,0.03]`, growing-mode projection in `[0.01,0.03]`, and alignment at
  least `0.99`. Add the exact initial `4^3` cancellation-dominated lattice to
  the periodic Ewald certification; its short-step kick predicts response near
  `0.019` of ballistic RMS. The trajectory window is therefore an
  independently anchored integration regression, not an analytic LCDM proof.
- Require the power-derived and direct displacement-mode growth increments to agree within
  relative `1e-3`.

### Consequences

- Positive: collisionless and gas state now consume one TreePM normalization,
  with their cosmological response applied once in the owning update layer.
- Positive: hydro source time has one code-unit timeline authority and the
  energy source remains correct away from `gamma=5/3`.
- Positive: gravity timestep length and acceleration now use the same comoving
  coordinate system instead of mixing `epsilon_com` with a scale-free or
  peculiar-velocity acceleration.
- Positive: the DMO gate rejects zero-growth/ballistic success and couples two
  independent growth estimators without overstating the response windows.
- Reproducibility change: physical-unit cosmological gas evolution is
  numerically corrected. Noncosmological `a=1`, `H=0` source behavior is
  unchanged.
- No config key, normalized dump, snapshot field, restart field, or provenance
  schema is added.

### Evidence references

- `include/cosmosim/hydro/hydro_core_solver.hpp`
- `src/hydro/hydro_core_solver.cpp`
- `src/workflows/reference_workflow.cpp`
- `tests/unit/test_hydro_core_solver.cpp`
- `tests/validation/test_dmo_zeldovich_workflow.cpp`
- `docs/hydro_core_solver.md`
- `docs/validation_plan.md`

## 2026-07-13 — ADR-REPAIR-GRAVITY-017: Periodic TreePM correctness profile and small-cluster transport boundary

### Status
Accepted (gravity correctness and reproducibility repair)

### Context
Periodic short-range trees were formerly built from ordinary wrapped Euclidean
coordinates even though traversal later used minimum-image deltas. A compact
cluster crossing a box seam could therefore acquire a box-scale root,
misleading COM and quadrupole moments, inconsistent MAC geometry, and unsafe
remote pruning bounds. Empty source ranks also lacked an explicit hierarchy
record, gravity packets relied on incomplete identity contracts, the PM slab
authority map did not exactly match FFTW-MPI's fixed input blocks at all rank
counts, and production cadence greater than one lacked a time-consistent
long-range predictor. The previous default CIC/deconvolution-off profile also
did not meet the new independent Ewald force-accuracy target for the covered
mesh/split fixtures. The workflow additionally hard-coded `G_code=1`, inserted
an `a^2` factor into both TreePM branches even though the dynamical update
layers already own the scale-factor response, and allowed production mixed-rung
scheduler configurations without per-element kick/drift epochs.

### Decision

- For periodic TreePM only, derive a per-axis largest-circular-gap unwrapped
  source frame and build topology, COMs, quadrupoles, raw second moments, node
  bounds, and MAC geometry in that coherent frame. Keep wrapped caller state and
  nonperiodic standalone tree behavior unchanged.
- Correct the quadrupole acceleration sign for the repository's
  `x_com - x_target` displacement convention. Store the raw second-moment trace
  and evaluate softened/screened second-order expansions through derivatives of
  the complete radial kernel.
- Add the relative force-error MAC
  `G M l^2 <= alpha max(|a_previous|, a_floor) r^4`, with a target-containing
  open guard and deterministic COM-distance fallback when compatible force
  history is unavailable.
- Align the public `TreeGravityOptions::opening_theta` default with the typed
  production/Ewald-certified value `0.7`. Callers that relied on the former
  aggregate-construction default `0.6` must request it explicitly.
- Encode tree hierarchy and short-range request/response packets explicitly as
  versioned little-endian records. Carry source/owner identity, exchange
  sequence, decomposition epoch, force epoch, geometry frame, target identity,
  and exact response coverage. Represent an empty local tree with a zero-source
  root sentinel rather than omitting rank participation.
- Coordinate every throwable TreePM entry/residual phase on the actual MPI
  world. Defer constructor-local layout/context rejection to solve-entry
  consensus; vote periodic unwrap/local-tree failures before residual exchange;
  vote PM halo-cache commit, compact active/PM/zoom target preparation, and
  worst-case traversal-stack reservation before later collectives; and place
  response counts/displacements plus payload/acceleration scratch in
  a reusable, memory-accounted exchange workspace whose throwing growth occurs
  inside coordinated preflight/preparation phases.
- Permit explicit active-target coordinate lanes so a zero-source rank can own
  targets without fabricating source mass.
- Match the PM x-slab owner map to FFTW-MPI fixed blocks:
  `B=ceil(Nx/world_size)`, rank `r` owns
  `[min(rB,Nx), min((r+1)B,Nx))`. Exchange force halos by posting nonblocking
  receives before sends and waiting collectively, with side/sequence tags and
  zero-width-slab support.
- Treat legal zero-width FFTW-MPI ranks as full solve participants: retain
  logical zero extent, provide backend-safe dummy allocation, and enter plan
  creation, transforms, and all PM collectives.
- Convert the physical SI Newton constant once from frozen length/mass/velocity
  units with `core::newtonGravitationalConstantCode`. PM and the complementary
  tree consume the same `G_code` and return a scale-free comoving kernel `A`.
  Collisionless KDK and the gas conservative source apply Hubble
  drag/expansion and the `A/a^2` response to their respective state.
- Make TSC with matched transfer-window deconvolution the normalized production
  default because it passes the independent small-N Ewald gate. Keep CIC
  available as a named diagnostic/compatibility profile; do not lower the gate
  to make it pass.
- Raise the normalized `treepm_rcut_cells` default from 4.5 to 6.25 while
  retaining positive explicit lower values for compatibility/diagnostics. With
  `asmth_cells=1.25`, the certified default has `r_cut/r_s=5`; the historical
  cutoff produced roughly nine-percent error at the cutoff transition.
- Require periodic `r_cut < min(Lx,Ly,Lz)/2` during typed loading and again at
  direct TreePM solve entry. The residual evaluates one minimum image per
  source, so coarse invalid decks must increase PM resolution or reduce the
  cutoff rather than relying on ambiguous multi-image reach.
- Restrict production `treepm_update_cadence_steps` to exactly one and
  `hierarchical_max_rung` to zero. Every rank-coordinated production
  force-refresh surface rebuilds PM. Retain lower-level reuse only as a
  test/future-integration seam; missing or incompatible cache state makes reuse
  fail rather than silently changing the requested operation.
- When `gravity.pm_long_range_field` is emitted before cadence commit, populate
  its opportunity, version, build step, and build scale factor from the pending
  integrator-issued `PmRefreshDirective`/solver decision. Do not pair a current
  refresh message with the previous committed `PmSynchronizationState`; the
  event remains a diagnostic observer and does not change cadence ownership.
- Serialize gravity operational-event doubles in scientific `max_digits10`
  form. Physical `G_code`, relative-MAC floors, derived split/cutoff/softening
  scales, force norms, and scale factors must not be rounded into apparent
  zeros by fixed six-decimal formatting.
- Validate the transient PM cache against force epoch, field-build scale,
  `G_code`, split scale, box axes, assignment, boundary, decomposition mode,
  and deconvolution. Deliberately exclude particle decomposition epoch because
  the mesh remains owned by fixed FFT slabs. Tree packets still carry the
  workflow epoch, which advances only after an actual ownership transition;
  dense-row force history is invalidated immediately on that transition.
- Clip isolated/open assignment and gather stencils at physical-domain
  boundaries rather than wrapping them to the opposite face.
- Reject `numerics.treepm_enable_window_deconvolution=true` during typed config
  validation whenever mode policy resolves to isolated/open gravity. The
  periodic default remains `true`; isolated decks must explicitly select
  `false`, and runtime code must not hide the incompatibility by mutating the
  frozen policy.
- Retain the hardened hierarchy-allgather plus target request/response
  all-to-all algorithm for correctness-first workstation and small-cluster use.
  Do not label it a locally essential tree or a neighbor-scalable transport.

### Consequences

- Positive: periodic translation, seam topology/moments, empty-rank
  participation, target-only ranks, and np2/np3/np4 slab/halo ownership now have
  executable regression coverage.
- Positive: the certified TSC profile meets `relative_L2 <= 1e-2` and
  `p99 <= 5e-2` against an implementation-independent Ewald reference on the
  covered fixtures; CIC failure remains visible.
- Positive: opening controls and relative-force parameters are normalized,
  validated, and recorded in `provenance_v7`; restart force history is consumed
  only when stable-ID/generation compatible.
- Positive: physical-unit workflows no longer depend on dimensionless
  `G_code=1`, and cosmological scale-factor powers are owned by the particle or
  gas update layer rather than duplicated in the force kernel.
- Positive: production mixed-rung execution now fails at config/restart
  validation instead of skipping elapsed time under incomplete per-element KDK
  state.
- Positive: occupancy, empty-rank, hierarchy/peer, PM solve/reuse/halo/slab,
  and tree build/multipole/opened-node counters make correctness and cost
  behavior directly observable without promoting diagnostics to runtime truth.
- Reproducibility change: configurations that relied on implicit defaults now
  select TSC+deconvolution and a 6.25-cell cutoff rather than CIC without
  deconvolution and a 4.5-cell cutoff. Explicit old settings retain their
  behavior, and provenance distinguishes the profiles.
- Tradeoff: TSC costs 27 mesh points per particle rather than CIC's 8; the
  larger cutoff also increases tree traversal, pair work, and distributed
  target/response volume.
- Tradeoff: communicator-wide hierarchy and record collectives remain a scaling
  ceiling. LET, sparse-peer PM transport, large-rank scaling certification,
  large DMO science validation, and GPU TSC parity remain future work.
- Positive: the controlled DMO gate now provides a versioned binned
  three-dimensional power-spectrum artifact, Fourier phase/coherence checks,
  serial-vs-MPI1 comparison, and np1--np4 stable-ID state comparison.
- Limitation: the 64-particle, two-step direct-DFT fixture is not a large-volume
  convergence, halo-statistics, cross-code, or publication-grade campaign;
  publishable DMO-production readiness is therefore not implied.

### Evidence references

- `docs/gravity_production_readiness.md`
- `docs/tree_gravity_solver.md`
- `docs/tree_pm_coupling.md`
- `docs/power_spectrum_diagnostics.md`
- `include/cosmosim/core/units.hpp`
- `tests/integration/test_tree_pm_coupling_periodic.cpp`
- `tests/integration/test_pm_slab_halo_exchange_mpi.cpp`
- `tests/validation/test_periodic_ewald_reference.cpp`
- `tests/validation/test_tree_pm_ewald_accuracy.cpp`
- `tests/validation/test_validation_phase2_mpi_gravity.cpp`
- `tests/validation/test_dmo_zeldovich_workflow.cpp`

## 2026-06-07 — ADR-INFRA-AGENT-INTERFACE-016: Repository-wide agent task contract

### Status
Accepted (documentation/interface decision)

### Context
The root `AGENTS.md` still described a narrow infrastructure-repair session even though current CHUI/CosmoSim work regularly asks AI agents to audit, patch, implement feature stages, and return clean source archives. That mismatch made it easy for agents to over-restrict valid implementation work or under-specify packaging and evidence duties.

### Decision
- Treat `AGENTS.md` as the canonical repository-wide agent contract.
- Add `docs/architecture/agent_task_interface.md` as the detailed task router for audit, repair, feature implementation, and packaging requests.
- Add `.github/copilot-instructions.md` as a short integration surface for GitHub Copilot-style agents that points back to the canonical contract.
- Preserve the CHUI user-facing name while explicitly forbidding opportunistic renames of existing `CosmoSim`/`cosmosim` compatibility identifiers.

### Consequences
- Positive: agent sessions now have one explicit mode-selection, evidence, ownership, and clean-zip policy.
- Positive: repair-only constraints no longer accidentally block requested feature implementation work.
- Positive: clean handoff artifacts now have an explicit repo-local exclusion contract.
- Tradeoff: documentation-only changes do not prove solver readiness; agents must still provide command-backed evidence for code changes.

### Evidence references
- `AGENTS.md`
- `docs/architecture/agent_task_interface.md`
- `.github/copilot-instructions.md`
- `docs/architecture/developer_workflow_contract.md`

## 2026-05-11 — ADR-INFRA-STAGE2-TIMESTEP-AUTHORITY-015: Scheduler-owned timestep truth and claim discipline

### Status
Accepted (Stage 2 infrastructure ownership decision)

### Context
Stage 2 repair hardened timestep-bin authority after audits found that scheduler hot metadata, public `ParticleSoa::time_bin` / `CellSoa::time_bin` lanes, restart payload mirrors, migration records, and TreePM cadence metadata could drift in documentation and release language. Project guardrails require public-interface documentation, explicit restart/schema behavior, and command-backed closure evidence in the same patch.

### Decision
- `HierarchicalTimeBinScheduler` is the only live authority for timestep-bin membership, next activation, active flags, pending transitions, and scheduler-built active sets.
- `PmSynchronizationState` is the scheduler-owned authority for PM kick-opportunity cadence, field-version, and refresh-commit legality.
- `ParticleSoa::time_bin`, `CellSoa::time_bin`, migration/transfer `time_bin` records, restart state copies, and `IntegratorState::time_bins` are mirrors or metadata. They may be serialized and validated for corruption, but may not be used as fallback scheduling authority.
- Particle-bound gas-cell mirrors must be validated and rebuilt through the parent gas particle scheduler entry, not by assuming cell rows and scheduler rows are equal-sized.
- Documentation and release notes must distinguish CPU-runnable Stage 2 scheduler-authority evidence from unproven Phase 3 multirate TreePM production claims.

### Consequences
- Positive: implementation, restart/schema docs, runtime truth map, and release language now present one Stage 2 owner model.
- Positive: stale mirror acceptance is treated as a restart compatibility failure, while valid v6 payload fields and schema version remain unchanged.
- Positive: Phase 3 language remains evidence-gated; Stage 2 scheduler authority does not imply production-proven hierarchical TreePM multirate synchronization.
- Tradeoff: new particle-registration, distributed scheduler identity exchange, and full PM/MPI/HDF5 evidence remain explicit follow-up blockers instead of closure claims.

### Evidence references
- `docs/time_integration.md`
- `docs/restart_checkpointing.md`
- `docs/output_schema.md`
- `docs/architecture/runtime_truth_map.md`
- `docs/repair_state_recap.md`
- `docs/repair_open_issues.md`
- `tests/integration/test_docs_scaffold.cmake.in`

## 2026-04-26 — ADR-INFRA-STAGE0-GATE-014: Stage 0 consolidation gate remains open (runtime-truth freeze not yet closed)

### Status
Accepted (gate record; Stage 0 not closed)

### Context
Stage 0 required a full consolidation gate across runtime-truth invariant suites, CPU debug/full suites, and dependency-enabled HDF5/PM+HDF5+FFTW presets before allowing new physics or performance work.

### Decision
- Run the Stage 0 targeted runtime-truth suite and full debug preset suites using repository presets.
- Record the gate closure artifact at `docs/repair/stage0_runtime_truth_freeze.md`.
- Keep Stage 0 marked **not closed** until currently failing contract/invariant tests are repaired and revalidated.
- Add named Stage 0 grouped test execution preset `test-stage0-runtime-truth-cpu-debug` in `CMakePresets.json`.

### Consequences
- Positive: Stage 0 closure status is explicit and command-backed, with a discoverable one-command Stage 0 runtime-truth suite.
- Positive: Prevents premature resumption of feature work while known runtime-truth contract failures remain.
- Tradeoff: Progression remains blocked behind existing red tests in CPU/HDF5/PM+HDF5+FFTW paths.

### Evidence references
- `docs/repair/stage0_runtime_truth_freeze.md`
- `docs/repair_open_issues.md`
- `docs/architecture/runtime_truth_map.md`
- `docs/architecture/adr_runtime_truth_ownership.md`
- `CMakePresets.json`

## 2026-04-25 — ADR-INFRA-OWNERSHIP-013: Single-source-of-truth runtime ownership policy

### Status
Accepted (infrastructure ownership policy baseline)

### Context
Stage 0 runtime-truth mapping (`docs/architecture/runtime_truth_map.md`) captured authoritative lanes and ambiguity, but ownership/mirror/invalidation/mutation authority contracts were not yet frozen in one ADR.

### Decision
- Add `docs/architecture/adr_runtime_truth_ownership.md` as the formal single-source-of-truth policy for runtime state ownership.
- Assign one authoritative owner per runtime domain (integrator time, scheduler bins/active sets, particle order/species grouping, gas identity mapping, softening policy lanes, config/provenance lanes, restart continuation payloads).
- Define explicit mirror/cache/view refresh + invalidation rules, mutation authority boundaries, restart vs snapshot semantics, and forbidden duplicate-authority patterns.

### Consequences
- Positive: reviewers and follow-up repair prompts now have one normative ownership contract tied to concrete modules/classes/fields.
- Positive: duplication lanes remain allowed only as documented mirrors/derived indices with explicit synchronization obligations.
- Tradeoff: stricter ownership language may require targeted follow-up tests when touching ambiguous domains.

### Evidence references
- `docs/architecture/runtime_truth_map.md`
- `docs/architecture/adr_runtime_truth_ownership.md`
- `include/cosmosim/core/time_integration.hpp`
- `include/cosmosim/core/simulation_state.hpp`
- `src/workflows/reference_workflow.cpp`
- `src/io/restart_checkpoint.cpp`
- `src/io/snapshot_hdf5.cpp`

## 2026-04-25 — ADR-INFRA-AUDIT-012: Stage 0 runtime-truth freeze initiated

### Status
Accepted (infrastructure audit baseline)

### Context
Runtime-state ownership across scheduler bins, species sidecars, gas identity mapping, and restart/snapshot overlap needed a code-verified truth map before additional repair edits.

### Decision
- Freeze a repository-local runtime truth map at `docs/architecture/runtime_truth_map.md` based on direct code inspection.
- Record explicit mutation/read/serialization owners for timestep bins, active sets, ordering/sidecars, gas identity, softening, config-derived runtime values, and restart/checkpoint truth.
- Record ambiguity/duplication and follow-up repair tickets without broad behavior edits in this stage.

### Consequences
- Positive: follow-up repair prompts can target one explicit audited baseline instead of inferred ownership.
- Positive: ambiguous mutation authority is called out explicitly rather than hidden in helper paths.
- Tradeoff: documentation overhead increases, but runtime behavior remains unchanged.

### Evidence references
- `docs/architecture/runtime_truth_map.md`
- `src/core/time_integration.cpp`
- `src/core/simulation_state_species.cpp`
- `src/workflows/reference_workflow.cpp`
- `src/io/restart_checkpoint.cpp`
- `src/io/snapshot_hdf5.cpp`

## 2026-04-19 — ADR-FEATURE-GRAVITY-011: Phase 2 distributed TreePM closeout hard-gate and failure-contract freeze

### Status
Accepted (feature-branch distributed TreePM closeout)

### Context
Phase 2 distributed TreePM implementation/testing was present across PM, tree exchange, workflow cadence, migration, and restart paths, but closeout lacked one command-backed hard-gate document and explicit negative-contract checks for common distributed failure modes.

### Decision
- Freeze a final closeout artifact at `docs/treepm_phase2_closeout.md` with explicit hard-gate commands.
- Extend distributed restart compatibility reporting with an explicit cadence validity lane (`pm_cadence_steps_match`).
- Add/expand test evidence for explicit failure handling:
  - rank-count/config mismatch,
  - unsupported PM decomposition mode,
  - communicator/layout mismatch,
  - missing distributed restart metadata,
  - inconsistent cadence metadata.

### Consequences
- Positive: documentation, validation, and runtime failure reporting now tell one coherent Phase 2 closeout story.
- Positive: distributed failure modes now fail loudly with typed/reportable reasons instead of implicit drift.
- Tradeoff: Phase 2 remains intentionally constrained to slab decomposition and deterministic long-range restart rebuild policy.

### Evidence references
- `docs/treepm_phase2_closeout.md`
- `include/cosmosim/parallel/distributed_memory.hpp`
- `src/parallel/distributed_memory.cpp`
- `tests/unit/test_parallel_distributed_memory.cpp`
- `tests/validation/test_validation_phase2_mpi_gravity.cpp`

## 2026-04-19 — ADR-FEATURE-GRAVITY-010: Freeze Phase 2 distributed TreePM contract and evidence floor

### Status
Accepted (feature-branch distributed gravity contract freeze)

### Context
Phase 1 TreePM contract is implemented and documented, but Phase 2 distributed ownership/communication semantics were not frozen in typed config and architecture docs. CI surface also centered on CPU/HDF5/PM paths plus an optional MPI release smoke job, which is insufficient for distributed gravity development cadence.

### Decision
- Freeze Phase 2 config-facing contract with two typed keys:
  - `numerics.treepm_pm_decomposition_mode` (currently only `slab`)
  - `numerics.treepm_tree_exchange_batch_bytes` (batch cap for tree export/import payloads)
- Add a dedicated distributed-gravity developer preset:
  - `mpi-hdf5-fftw-debug`
- Publish repository-specific Phase 2 architecture contract in `docs/treepm_phase2_distributed_contract.md`, including:
  - PM slab ownership semantics,
  - particle-owner vs slab-owner roles,
  - distributed FFT stage ordering,
  - tree export/import batch contract,
  - active-set and migration timing constraints,
  - restart continuity and deterministic limits.
- Explicitly record that pseudo-multi-rank tests are not completion evidence for distributed TreePM implementation.

### Consequences
- Positive: distributed gravity work now has an auditable typed/config/docs baseline without claiming algorithm completion.
- Positive: CI/build surface can exercise MPI+HDF5+FFTW debug configuration directly rather than only optional MPI release smoke.
- Tradeoff: only slab decomposition is valid in this phase; future modes require explicit ADR + typed/config/docs/test expansion.

### Evidence references
- `include/cosmosim/core/config.hpp`
- `src/core/config.cpp`
- `CMakePresets.json`
- `docs/treepm_phase2_distributed_contract.md`
- `.github/workflows/ci.yml`

## 2026-04-12 — ADR-REPAIR-BOUNDARY-005: Keep reference workflow assembly outside `core/`

### Status
Accepted (infrastructure boundary repair)

### Context
`include/cosmosim/core/reference_workflow.hpp` assembled diagnostics, physics callbacks, and I/O roundtrips directly from the `core/` include tree. This inverted the dependency direction guardrail (`core` foundational; higher layers depend on it).

### Decision
- Move concrete reference workflow assembly into a dedicated integration layer:
  - `include/cosmosim/workflows/reference_workflow.hpp`
  - `src/workflows/reference_workflow.cpp`
- Keep `core/` focused on abstract stage contracts and typed simulation state.
- Add a dependency-direction guard test that fails when `include/cosmosim/core/` or `src/core/` includes `analysis/`, `io/`, `physics/`, `workflows/`, or `app/` headers.

### Consequences
- Positive: Restores one-way dependency flow from higher orchestration layers into `core`.
- Positive: Makes boundary regressions auditable in CI via an explicit guard test.
- Tradeoff: Callers should migrate toward `cosmosim::workflows` naming; transitional aliases remain for now.

## 2026-04-07 — ADR-REPAIR-FREEZE-001: Freeze repair evidence before touching production code

### Status
Accepted (stabilization phase)

### Context
The repository already contains significant architecture work (SoA storage, hot/cold split, sidecars, typed config, provenance, staged integration contracts). Entering emergency stabilization without a frozen baseline risks accidental churn and invalid claims of closure.

### Decision
Before any further repair implementation:

1. Capture a command-backed repair-state recap.
2. Capture a command-backed open-issues ledger.
3. Freeze the current public-interface baseline from `include/cosmosim/**`.
4. Record both default-path green status and feature-path red status explicitly.

### Consequences
- Positive: Later prompts can cite a single factual baseline instead of re-discovering repo state.
- Positive: Prevents CPU-only green path from being mistaken for HDF5 feature-path closure.
- Positive: Supports minimal, auditable fixes by tying each open issue to file+command+symptom.
- Tradeoff: Adds documentation overhead but no production-code churn.

### Evidence references
- `docs/repair_state_recap.md`
- `docs/repair_open_issues.md`

## 2026-04-07 — ADR-REPAIR-AUDIT-002: Hold progression to P20 until dependency-enabled presets are proven

### Status
Accepted (emergency stabilization gate)

### Context
Final repair audit was run against required scope. CPU-only path is green, but this environment lacks HDF5 development libraries, preventing configure/build/test execution on required HDF5 and PM/HDF5/FFTW paths.

### Decision
Do **not** authorize progression to P20+ until all required dependency-enabled preset commands pass in an environment with HDF5 (and FFTW for PM path).

### Consequences
- Positive: Prevents invalid closure claims based on CPU-only evidence.
- Positive: Keeps dependency failure behavior explicit and actionable.
- Tradeoff: Progression is delayed until proper feature-path environment is available.

### Evidence references
- `docs/repair_state_recap.md`
- `docs/repair_closeout_report.md`

## 2026-04-07 — ADR-REPAIR-BLOCK-003: Keep P20 progression blocked until PM HDF5+FFTW validation passes

### Status
Accepted (stabilization gate enforcement)

### Context
Dependency-enabled stabilization validation was re-run with presets in this environment.

- CPU-only configure/build/test passed.
- HDF5 configure/build/test passed.
- PM HDF5+FFTW configure/build passed, but test preset failed (`2/36` failures):
  - `unit_pm_solver` assertion failure in `tests/unit/test_pm_solver.cpp:82`.
  - `integration_tree_pm_coupling_periodic` runtime failure with `rel_l2=18129.9` against required `<= 0.75`.

### Decision
Progression beyond P19 remains blocked. Do not authorize P20 until `test-pm-hdf5-fftw-debug` passes completely.

### Consequences
- Positive: Enforces dependency-enabled acceptance criteria and avoids CPU-only false closure.
- Positive: Keeps blocker scope narrow (PM/FFTW test path) without speculative refactoring.
- Tradeoff: Delivery remains gated pending PM HDF5+FFTW test stabilization.

### Evidence references
- `docs/repair_closeout_report.md`
- `docs/repair_open_issues.md`

## 2026-04-07 — ADR-REPAIR-GATE-004: Clear P20 gate only when all three preset validation paths pass

### Status
Accepted (post-repair gate policy)

### Context
P19 stabilization was previously blocked by PM FFTW-path numerical convention defects affecting:
- `unit_pm_solver` analytic mode directionality check.
- `integration_tree_pm_coupling_periodic` long-range amplitude consistency (`rel_l2` blow-up).

After targeted PM FFTW-path correction, required validation commands were re-run and all three required presets passed:
- `test-pm-hdf5-fftw-debug`
- `test-hdf5-debug`
- `test-cpu-debug`

### Decision
P20 gate is considered **CLEARED only if all three preset paths pass in the same validation cycle**:
1. CPU-only preset path
2. HDF5 preset path
3. PM HDF5+FFTW preset path

Any future regression in one path re-blocks the gate.

### Consequences
- Positive: Prevents partial-path closure and preserves dependency-enabled confidence.
- Positive: Keeps PM FFTW path first-class in stabilization acceptance.
- Tradeoff: Slightly longer validation cycle per gate decision.

### Evidence references
- `docs/repair_closeout_report.md`
- `docs/repair_open_issues.md`


## 2026-04-14 — ADR-REPAIR-CONTINUATION-005: Remove normalized-config self-hash and reject incomplete restart continuations

### Status
Accepted (targeted infrastructure repair)

### Context
The normalized config artifact embedded its own hash, which made the canonical text internally inconsistent and non-round-trippable. Separately, restart payloads/schema omitted stellar-evolution continuation lanes while the integrity hash and compatibility contract still claimed exact continuation.

### Decision
- Define the normalized config hash over the normalized config text alone; do not embed a self-hash line inside the canonical text.
- Require continuation metadata agreement: normalized text hash == stored normalized-config hash == provenance config hash.
- Bump restart schema to `cosmosim_restart_v3` and reject older incomplete restart artifacts.
- Persist and hash the full stellar-evolution star sidecar payload.

### Consequences
- Positive: normalized config dumps are deterministic, parseable, and hash-consistent.
- Positive: restart artifacts once again match their exact-continuation claim for currently serialized state.
- Tradeoff: `cosmosim_restart_v2` artifacts are intentionally incompatible because they omit required continuation state.

### Evidence references
- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`

## 2026-04-18 — ADR-FEATURE-GRAVITY-006: Feature-branch exception for PM periodic solver contract hardening

### Status
Accepted (feature-branch gravity upgrade)

### Context
Repository repair guardrails include a repair-only prohibition on solver-behavior edits. This branch is explicitly scoped as a gravity-upgrade feature branch requiring PM solver contract hardening, including potential output formalization and stricter periodic Poisson convention tests.

### Decision
- Authorize PM solver behavior/documentation/test changes in this branch only, while preserving existing architecture direction:
  - no dependency-direction inversion,
  - explicit typed/configured runtime inputs,
  - auditable solver sign/normalization conventions,
  - evidence-backed closure via preset configure/build/test commands.

### Consequences
- Positive: Enables scientifically auditable PM potential+force contract without violating boundary discipline.
- Positive: Removes ambiguity around `δρ`, zero mode, Fourier `k` mapping, and force sign/normalization.
- Tradeoff: Branch intentionally changes numerical-solver contract surface and therefore requires same-patch docs/tests updates.

### Evidence references
- `include/cosmosim/gravity/pm_solver.hpp`
- `src/gravity/pm_solver.cpp`
- `docs/pm_gravity_solver.md`
- `tests/unit/test_pm_solver.cpp`
- `tests/integration/test_pm_periodic_mode.cpp`
- `tests/validation/test_validation_integration.cpp`


## 2026-04-18 — ADR-FEATURE-GRAVITY-007: Enable PM assignment/gather upgrade (CIC+TSC) on gravity feature branch

### Status
Accepted (feature-branch gravity upgrade)

### Context
The gravity-upgrade branch requires end-to-end PM assignment/gather configurability (`cic`, `tsc`) with auditable transfer-kernel behavior. This intentionally changes solver behavior versus prior CIC-only runtime gating.

### Decision
- Keep PM assignment and gather strictly matched by selected runtime scheme (no mixed silent modes).
- Implement TSC alongside CIC in deposition and interpolation paths.
- Make k-space deconvolution scheme-aware and apply it to the combined transfer operator with a protective floor.
- Keep conservative defaults (`cic`, deconvolution off) and record behavior in same-patch docs/tests.

### Consequences
- Positive: `kTsc` is no longer a placeholder; config/runtime/docs/tests align.
- Positive: PM transfer behavior is explicit and auditable across assignment choices.
- Tradeoff: TSC increases stencil work (27-point 3D) compared to CIC (8-point 3D).
- Tradeoff: Explicit CUDA PM path remains CIC-only in this build and rejects `tsc` rather than silently diverging.

### Evidence references
- `include/cosmosim/gravity/pm_solver.hpp`
- `src/gravity/pm_solver.cpp`
- `src/workflows/reference_workflow.cpp`
- `docs/pm_gravity_solver.md`
- `docs/configuration.md`
- `tests/unit/test_pm_solver.cpp`
- `tests/integration/test_pm_periodic_mode.cpp`
- `tests/validation/test_validation_integration.cpp`

## 2026-04-18 — ADR-FEATURE-GRAVITY-008: TreePM split/cutoff contract upgrade on gravity feature branch

### Status
Accepted (feature-branch gravity upgrade)

### Context
This branch stage requires solver-behavior edits in the TreePM coupling layer to make the AREPO/GADGET-style split contract explicit and runtime-configurable from normalized config controls (`asmth_cells`, `rcut_cells`). Existing residual traversal lacked explicit RCUT truncation behavior.

### Decision
- Treat mesh-cell controls as authoritative and derive:
  - `r_s = asmth_cells * Δmesh`
  - `r_cut = rcut_cells * Δmesh`
  - `Δmesh = box_size / PMGRID`
- Keep Gaussian split family:
  - PM: `exp(-k^2 r_s^2)`
  - residual tree: Gaussian complementary short-range factor.
- Enforce explicit residual cutoff in traversal:
  - node-level AABB pruning,
  - node-acceptance guard requiring full in-cutoff enclosure,
  - leaf pair skipping beyond `r_cut`.
- Expand diagnostics/docs/tests so split parameters, composition checks, and cutoff pruning are auditable.

### Consequences
- Positive: Removes hidden split heuristics and dual-source ambiguity.
- Positive: `rcut_cells` now changes actual residual traversal work and force contribution domain.
- Positive: Improves reproducibility/audibility via runtime diagnostics and same-patch docs/tests.
- Tradeoff: Tree residual traversal now performs additional cheap AABB distance bounds to guarantee safe pruning.

### Evidence references
- `include/cosmosim/gravity/tree_pm_split_kernel.hpp`
- `include/cosmosim/gravity/tree_pm_coupling.hpp`
- `src/gravity/tree_pm_coupling.cpp`
- `src/workflows/reference_workflow.cpp`
- `docs/tree_pm_coupling.md`
- `tests/unit/test_tree_pm_split_kernel.cpp`
- `tests/integration/test_tree_pm_coupling_periodic.cpp`

## 2026-04-18 — ADR-FEATURE-GRAVITY-009: Enable explicit PM long-range cadence/reuse in reference runtime

### Status
Accepted (feature-branch gravity upgrade)

### Context
`numerics.treepm_update_cadence_steps` existed in typed/normalized config but runtime still hard-rejected values above `1`, leaving no auditable PM refresh/reuse behavior despite hierarchical stepping and active-set gravity callbacks.

### Decision
- Treat `treepm_update_cadence_steps` as an authoritative runtime control in the reference workflow.
- Apply cadence per gravity kick opportunity (`gravity_kick_pre`, `gravity_kick_post`):
  - refresh PM long-range mesh solve every `N` opportunities,
  - reuse cached long-range field otherwise,
  - always recompute short-range tree residual each kick.
- Keep default conservative (`N=1`) to preserve immediate-refresh baseline behavior.
- Record refresh/reuse metadata (field version, build step, build scale factor) in workflow report and operational events.

### Consequences
- Positive: cadence key now changes real runtime behavior and is diagnosable.
- Positive: PM reuse is explicit, deterministic, and traceable instead of hidden in callback internals.
- Tradeoff: Phase 1 cadence is single-rank and intentionally simple; it is not yet a full multirate integrator.

### Evidence references
- `src/workflows/reference_workflow.cpp`
- `include/cosmosim/workflows/reference_workflow.hpp`
- `include/cosmosim/gravity/tree_pm_coupling.hpp`
- `src/gravity/tree_pm_coupling.cpp`
- `tests/integration/test_reference_workflow_end_to_end.cpp`
- `docs/time_integration.md`
- `docs/reference_workflow.md`
- `docs/pm_gravity_solver.md`

## 2026-04-18 — ADR-FEATURE-GRAVITY-010: Persist TreePM provenance contract in snapshot/restart/log outputs

### Status
Accepted (feature-branch gravity upgrade)

### Context
TreePM controls and derived scales materially affect force accuracy/reproducibility, but output artifacts did not carry a complete gravity contract across snapshots, restarts, and runtime operational events.

### Decision
- Extend `core::ProvenanceRecord` with gravity provenance fields covering:
  - controls (`pm_grid`, assignment, deconvolution, `asmth_cells`, `rcut_cells`, cadence),
  - derived scales (`Δmesh`, `r_s`, `r_cut`),
  - softening policy/kernel/epsilon and PM FFT backend.
- Bump snapshot schema identity to `gadget_arepo_v2` / `schema_version=2` and persist the new provenance attributes in `/Provenance`.
- Keep restart schema `cosmosim_restart_v3`; serialize gravity provenance through existing provenance dataset and integrity hash path.
- Emit structured runtime events (`gravity.treepm_setup`, enriched `gravity.pm_long_range_field`) that expose the same gravity contract at setup/refresh points.

### Consequences
- Positive: Snapshot/restart/runtime logs can reconstruct the exact TreePM runtime contract used for a run.
- Positive: Derived physical/comoving scales are auditable without re-deriving from partial controls.
- Tradeoff: Snapshot schema version increased; downstream readers that hardcode schema version `1` must accept `2`.

### Evidence references
- `include/cosmosim/core/provenance.hpp`
- `src/core/provenance.cpp`
- `include/cosmosim/io/snapshot_hdf5.hpp`
- `src/io/snapshot_hdf5.cpp`
- `src/workflows/reference_workflow.cpp`
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_provenance_roundtrip.cpp`
- `tests/integration/test_reference_workflow_end_to_end.cpp`

## 2026-04-19 — ADR-FEATURE-GRAVITY-011: TreePM Phase 1 integration hard-gate closure contract

### Status
Accepted (feature-branch gravity upgrade)

### Context
Final TreePM Phase 1 integration requires coherence across runtime wiring, config normalization, provenance/output audibility, validation docs, and command-backed evidence. Earlier stage changes introduced PMGRID/ASMTH/RCUT/cadence controls and TreePM runtime behavior, but Phase 1 closure requires a single explicit contract and hard-gate evidence checklist.

### Decision
- Adopt an explicit Phase 1 integration hard gate that must be evidenced together:
  1. no hard-coded PM mesh in workflow (PMGRID is config-driven),
  2. no hidden split assumptions (ASMTH/RCUT drive runtime `r_s`/`r_cut`),
  3. deterministic repeated-run reproducibility evidence on identical config,
  4. force-error map artifact over PMGRID/ASMTH/RCUT sweep.
- Record the final contract and evidence procedure in `docs/treepm_phase1_closeout.md` and align cross-doc language (`configuration`, `tree_pm_coupling`, `validation_ladder`, `validation_plan`).

### Consequences
- Positive: Phase 1 completion claims are tied to auditable artifacts and commands rather than narrative-only status.
- Positive: Runtime/config/docs/output semantics are synchronized for TreePM controls.
- Tradeoff: Phase 1 closeout remains scoped to single-rank periodic TreePM and non-Ewald references.

### Evidence references
- `docs/treepm_phase1_closeout.md`
- `tests/integration/test_reference_workflow_end_to_end.cpp`
- `bench/bench_tree_pm_force_error_map.cpp`
- `validation/artifacts/tree_pm_force_error_map.csv`

## 2026-04-19 — ADR-INFRA-RESTART-012: Version distributed TreePM continuation state in restart/provenance

### Status
Accepted (infrastructure repair)

### Context
Phase 2 distributed TreePM restart continuation lacked an explicit, versioned contract for rank ownership, PM slab layout, cadence state, and long-range field refresh metadata. This made rank/layout mismatch debugging ambiguous and allowed implicit continuation assumptions.

### Decision
- Bump restart schema to `cosmosim_restart_v4` and require a new restart group `/distributed_gravity/state` carrying serialized `parallel::DistributedRestartState` (schema_version `2`).
- Persist explicit distributed continuation fields:
  decomposition epoch, owning-rank table, PM grid/layout metadata, kick-opportunity cadence state, long-range refresh/version metadata, and explicit long-range restart policy.
- Adopt and enforce restart policy `deterministic_rebuild` for PM long-range field continuation.
  Cached long-range arrays are not serialized; restart resumes with deterministic rebuild at the next cadence-triggered refresh.
- Extend provenance (`provenance_v4`) with axis-aware PM grid/mesh-spacing metadata in addition to distributed restart diagnostics (epoch/world/grid/slab signature/kick+field version/policy).
- Add compatibility diagnostics API `evaluateDistributedRestartCompatibility(...)` to emit explicit mismatch reasons.

### Consequences
- Positive: distributed restart semantics are explicit, versioned, and integrity-hashed.
- Positive: rank/layout/config mismatch debugging has typed reporting instead of opaque failures.
- Tradeoff: `cosmosim_restart_v3` artifacts are intentionally incompatible with `v4`.

### Evidence references
- `include/cosmosim/io/restart_checkpoint.hpp`
- `src/io/restart_checkpoint.cpp`
- `include/cosmosim/parallel/distributed_memory.hpp`
- `src/parallel/distributed_memory.cpp`
- `include/cosmosim/core/provenance.hpp`
- `src/core/provenance.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/unit/test_parallel_distributed_memory.cpp`

## 2026-04-20 — ADR-INFRA-GRAVITY-013: Freeze TreePM Phase 3 contract and maturity map before numerics work

### Status
Accepted (infrastructure planning/audit freeze)

### Context
The repository now has a documented and tested Phase 1/2 TreePM baseline (including distributed slab ownership/message contracts and Phase 2 MPI validation gates), but does not yet have an explicit Phase 3 closure contract tied to evidence boundaries. Without a frozen Phase 3 contract, later implementation prompts risk over-claiming completion or conflating scaling artifacts with production certification.

### Decision
- Add a Phase 3 contract document (`docs/treepm_phase3_contract.md`) that records:
  - audited current baseline facts,
  - required Phase 3 additions,
  - forbidden completion claims before evidence,
  - explicit hard closure gate,
  - staged requirement-to-files maturity map.
- Update validation/release planning docs so Phase 3 is tracked as pending contract work, not completed functionality.
- Do not implement Phase 3 numerics in this patch.

### Consequences
- Positive: later implementation work has a single auditable closure contract and anti-drift guardrails.
- Positive: prevents false closure claims by separating “Phase 2 correctness evidence” from “Phase 3 maturity evidence.”
- Tradeoff: this ADR adds planning/documentation overhead before any algorithmic implementation.

### Evidence references
- `docs/treepm_phase3_contract.md`
- `docs/validation_ladder.md`
- `docs/validation_plan.md`
- `docs/releases/known_issues.md`

## 2026-04-21 — ADR-FEATURE-GRAVITY-012: Record Phase 3 final integration closeout as incomplete

### Status
Accepted (closeout audit record)

### Context
A final integration/coherence audit was run against the Phase 3 contract surfaces, campaign artifacts, and targeted gravity validation lanes. The cycle included PM/HDF5+FFTW configure/build and targeted integration/validation commands, plus manifest repro collection.

### Decision
- Publish explicit closeout artifact `docs/treepm_phase3_closeout.md`.
- Record this cycle as **not Phase 3 complete** due to command-backed blockers.
- Keep claim discipline explicit: no Phase 3 completion language until blockers are cleared in one evidence cycle.

### Consequences
- Positive: avoids false completion claims and keeps evidence/limitations auditable.
- Positive: provides an exact blocker list tied to failing commands.
- Tradeoff: release/readme language must continue to use non-closure wording for Phase 3.

### Evidence references
- `docs/treepm_phase3_closeout.md`
- `docs/repair_state_recap.md`
- `docs/repair_open_issues.md`
- `docs/releases/known_issues.md`


## 2026-05-10 — ADR-INFRA-STAGE2-TIMESTEP-OWNERSHIP-014: Freeze timestep owner before Stage 2 repairs

### Status
Accepted (infrastructure audit; PM cadence ownership refined by ADR-INFRA-STAGE2-TIMESTEP-AUTHORITY-015)

### Context
Stage 2 timestep repair has multiple existing timestep-related lanes: the hierarchical scheduler hot state, state `time_bin` mirrors, integrator metadata, restart scheduler payloads, migration/reorder preservation lanes, and TreePM cadence counters. Without a pre-change ownership map, any repair could create split-brain timestep truth.

### Decision
- Adopt `core::HierarchicalTimeBinScheduler` as the intended owner for per-element discrete timestep truth: committed bin, pending transition, activation tick, active flag/cache for the current scheduler substep, and scheduler integer tick.
- Treat `ParticleSoa::time_bin` and `CellSoa::time_bin` as mirrors derived from scheduler persistent state, including in restart payloads.
- Treat `core::IntegratorState` as the owner for global continuous time/step scalars only; its time-bin context is metadata and must not become per-element authority.
- Freeze the pre-repair audit finding that TreePM PM cadence was serialized through the gravity workflow/distributed restart cadence state; ADR-INFRA-STAGE2-TIMESTEP-AUTHORITY-015 refines the live ownership decision by assigning PM kick-opportunity cadence and refresh-commit legality to `PmSynchronizationState`, while `gravity_kick_opportunity` remains a kick-stage counter rather than a scheduler tick.
- Publish the detailed inventory and unsafe hidden-state list in `docs/architecture/stage2_timestep_ownership_audit.md` before behavior-changing Stage 2 work.

### Consequences
- Positive: later Stage 2 patches have a single authority target for scheduler legality, bin commit, activation, mirror sync, restart, and transform repair.
- Positive: restart scheduler state, state mirrors, active-set descriptors, structural transforms, and PM cadence are explicitly separated.
- Required follow-up: behavior-changing repairs must add assertions/tests for mirror freshness, structural-transform scheduler rebuild/reimport, and restart read-side scheduler/mirror compatibility before claiming closure.
- Tradeoff: this ADR classifies existing ambiguous lanes as mirrors or unsafe hidden state without fixing them in this documentation-only patch.

### Evidence references
- `docs/architecture/stage2_timestep_ownership_audit.md`
- `include/cosmosim/core/time_integration.hpp`
- `src/core/time_integration.cpp`
- `src/workflows/reference_workflow.cpp`
- `src/io/restart_checkpoint.cpp`
- `include/cosmosim/parallel/distributed_memory.hpp`

## 2026-07-16 — ADR-RUNTIME-OUTPUT-TIMELINE-016: Persist ordered code-time output events

### Status
Accepted

### Context
The production workflow could schedule snapshots only by step modulo. It could
cross a requested physical/code-time event with a larger KDK interval, and its
checkpoint cadence payload was not consumed when resuming the workflow.

### Decision
- Keep step-modulo cadence and add optional `output.snapshot_interval_time_code`,
  anchored at `numerics.time_begin_code`.
- Clip the proposed KDK interval at the earliest configured endpoint or output
  event; never dispatch the event from an unsafe boundary.
- Make the next code-time event restart-authoritative in
  `cosmosim_restart_v21`, together with the interval and payload integrity hash.
- Persist the pre-clip interval as the next-step proposal so restart and direct
  continuation use the same post-event timestep.
- Read v20 with the new cadence lanes explicitly disabled; do not synthesize
  event authority absent from the old payload.

### Consequences
- Positive: required code-time snapshots/checkpoints are reached without
  overshoot and cadence survives restart.
- Positive: normalized config, restart payload, diagnostics, tests, and user
  documentation describe the same authority.
- Tradeoff: arbitrary explicit redshift/scale-factor event lists remain outside
  this narrow schema; the supported beyond-modulo event is periodic code time.

### Evidence references
- `include/cosmosim/core/config.hpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `src/core/config.cpp`
- `src/io/restart_checkpoint.cpp`
- `src/workflows/reference_workflow.cpp`
- `tests/integration/test_reference_workflow_end_to_end.cpp`


## 2026-07-17 — Task-scoped runtime resource grants replace public owner-context escape

**Decision.** Public workflow stage views contain no `StepContext` pointer and
owner interfaces expose no `stageContext()` helper. `RuntimeExecutionPlan` binds
the frozen task declaration's resource list for one invocation, clears it with
RAII, and the source-private built-in access bridge validates required keys and
modes before entering existing rung-zero code. Typed-view allowlists reject
out-of-scope declarations at registry registration.

**Rationale.** Campaign A decomposition created organizational owners but left a
protected path by which an external owner subclass could mutate unrelated state.
The new boundary prevents public/custom modules from recovering global mutable
state while preserving current numerical order inside trusted built-in owners.

**Consequences.** Descriptor resource lists are executable stage grants rather
than documentation-only metadata. They do not dynamically instrument trusted
built-in memory accesses, and descriptor scope remains stage composition—not
restart-schema, migration-schema, or timestep-criteria registration.

## 2026-08-02 — Authoritative `.param.txt` configuration format

Accepted `docs/architecture/adr_param_txt_configuration_format.md`. Simulation inputs, examples, release profiles, normalized copies, and fixtures use lower-snake-case `.param.txt`. YAML/JSON input claims in older design or naming material are superseded; generated JSON metadata remains valid.

## 2026-08-29 — ADR-RUNTIME-MEMORY-GOVERNOR-017: Establish one process memory-admission authority

### Status
Accepted for Campaign M1B foundation

### Context
Post-M1A PM/scheduler structure has bounded the dominant known memory pathologies,
but pre-run estimation alone cannot prevent a later retained workspace from
committing memory after headroom has changed. The repository already has
`MemoryReport` ownership accounting and a separate `RuntimeResourceLease` for
stage freshness; neither should be repurposed into runtime memory admission.

### Decision
- Own exactly one `core::MemoryGovernor` at the reference-workflow composition
  root and borrow it through `RuntimeServices`; do not add a global singleton.
- Use `parallel.process_memory_budget_bytes` and the existing process external
  reserve/output-overlap/safety-margin policy as the only production budget
  source.
- Represent admitted demand with move-only RAII `MemoryReservation`, with
  pending and committed bytes tracked separately.
- Reconcile ordinary owned-capacity reporting as baseline memory while excluding
  entries marked `governed_commitment`, because those bytes are already present
  in the governor committed counter.
- Govern `MonotonicScratchAllocator` block growth first. A block reservation is
  committed only after backing allocation succeeds and remains committed across
  scratch `reset()` until the block is destroyed.
- Keep `RuntimeResourceLease` unchanged and independent.
- Keep M1B pressure thresholds internal (85% Amber, 95% Red); broad adaptive
  scheduling and broad allocation coverage are deferred to M1C/later work.

### Consequences
- Positive: large target-relevant allocations can now have deterministic
  reserve-before-allocate admission without routing small/hot-loop allocation
  traffic through a central lock.
- Positive: retained physical capacity and logical scratch use no longer share
  a false lifetime in governance accounting.
- Positive: pre-run fit estimation and runtime admission share the same process
  ceiling and reserve/margin semantics while serving different questions.
- Required follow-up: M1C must extend governed reservation coverage, add
  whole-process RSS/PSS reconciliation and distributed acceptance evidence, and
  coordinate any new rank-dependent runtime rejection before later collectives.

### Evidence references
- `include/cosmosim/core/memory_governor.hpp`
- `src/core/memory_governor.cpp`
- `docs/memory_governance.md`
- `docs/state_model_memory_layout.md`
- `docs/profiling.md`

## 2026-09-07 — ADR-M2D-ACCEPTANCE-001: Explicit confidence in major-task memory estimates

**Status:** Accepted as a fail-closed runtime-contract extension; full M2D
performance/distributed acceptance remains provisional.

The runtime registry must distinguish an actual complete incremental peak from
a useful partial owner workspace model. A zero estimate or missing callback is
not proof of zero allocation. Add state-dependent estimates to the existing
registry, preserve legacy callbacks and custom owner interfaces, and expose
qualified task lookup for admission tests and future scheduling research.
Owner-managed physical leases remain authoritative; complete preflight leases
are temporary, while dispatcher-owned leases span execution. Partial models
cannot authorize overlap or force an adaptive kernel to fail at an unnecessarily
large nominal batch. No second executor, governor, configuration surface, or
scientific ordering authority is introduced. Existing gravity and output
allocation models are reused rather than creating competing formulae.

The source/hydro/output models remain explicitly incomplete where physical
population growth, AMR transactions, or I/O metadata lack a complete bound.
This is an intentional safety boundary, not acceptance closure. No solver
precision, equations, tolerances, force criteria, CFL rules, refinement policy,
restart schema, or provenance schema changes. Memory failures can occur earlier
for complete owner models; the reference-workflow budget regression accepts
both the original gravity phase and the new typed task-preflight boundary.

## 2026-09-07 — ADR-M2D-CODING-001: Source population growth and collective SFC preparation

**Status:** Source implementation; distributed runtime acceptance pending.

The canonical source owner retains the existing star/BH ID-precommit and
scheduler authorities. A phase-local checked replacement transaction reserves
canonical/species/scheduler backing growth before physical birth mutation and
reconciles retained capacity on success or partial failure. This is a narrow
core utility, not a replacement allocator or memory governor. The existing
source equations, birth identity, source iteration, and serial stage order are
unchanged. A failed admission may reject a step earlier; successful scientific
results and physical restart schema remain unchanged.

The distributed SFC planner now agrees on rank-local entry and preparation
failures before collectives and exposes test-only fault injection at the new
memory-refinement preparation boundaries. The SFC ownership geometry and
feasibility-first cut algorithm are unchanged. Full MPI fault/restart and
whole-process memory certification remain required before production M2D
acceptance. Incomplete owner contracts cannot authorize overlap.

## 2026-09-21 — ADR-FIRSTLIGHT-SCIENCE-IO-001: Decouple science-snapshot topology from compute/restart topology

**Status:** Accepted for the small/moderate-rank production path; fixed-count aggregation remains future work.

Science snapshots are analysis products, not restart partitions. The authoritative typed `output.snapshot_layout` policy therefore selects one shared HDF5 science file for serial execution and for MPI when the linked HDF5 is Parallel-HDF5 capable. MPI ranks write only owned non-overlapping global hyperslabs; no complete particle state is gathered onto rank zero. The writer uses collective Parallel-HDF5 file/metadata construction with direct MPI-IO hyperslab payload writes, bounded per-rank staging, distributed exact partition readback, collective failure agreement, and atomic publication from `.partial` to `snap_###.hdf5`. Root post-publication validation may omit its otherwise O(N_global) duplicate ID buffer only after the distributed exact readback has established file IDs equal globally unique authoritative runtime IDs.

`sharded` remains an explicit compatibility topology. `aggregated` plus `snapshot_num_files` reserves a future fixed-count I/O topology but fails closed until a real backend is qualified. Restart/checkpoint ownership, schema, rank-local files, and same-topology continuation remain unchanged. This decision changes output topology and operator interoperability only; TreePM, timesteps, particle physical state, deterministic reductions, and restart semantics are unaffected.

## 2026-10-03 — ADR-A2-TREEPM-CADENCE-001: Rank-global PM refresh materialization before TreePM branching

**Status:** Source-level repair implemented; build and runtime qualification pending.

The core integrator remains the owner of PM-refresh opportunity semantics, but
rank-local predicates are now evidence rather than committed distributed
cadence truth. `TimeCoordinator`, which already owns workflow-level MPI stage
dispatch, reduces PM enablement, refresh-boundary legality, local-substep state,
long-range-field validity, and source-generation mismatch before any
branch-dependent TreePM collective. It deterministically materializes one
`PmRefreshDirective` on every rank. Mandatory invalidation is rank-global: any
source-generation mismatch or invalid PM field requires a shared repair at a
legal boundary; refresh permission is conservative and must be legal on every
rank except the existing initial-bootstrap synchronization exception.

`PmSynchronizationState` now supports side-effect-free opportunity preview and
a separate commit. The committed kick opportunity/version therefore advances
only after the TreePM callback successfully consumes the globally authorized
directive. Gravity retains its cadence consensus guard and strengthens exact
integer/epoch checks to rank min/max agreement; it does not synthesize equality
for scientifically meaningful source generations. The repair does not force PM
refresh every step, serialize TreePM, alter force tolerances, or change restart
schema semantics.

## P6 optional accuracy, refit and DMO block integration — 20261006

Source implementation; qualification pending. Conservative defaults retain
strict residual acceptance, full tree rebuild, global KDK and existing
rank-total decomposition feedback. Numerical opt-ins belong to the existing
frozen config, not an alternate configuration surface.

Accuracy proxy and geometric convergence guards have separate authority.
Source generation, row-layout generation, ownership epoch, tree-build generation
and PM version remain distinct facts. Refit retains immutable leaf membership
and certifies it geometrically; no equality-of-N or empirical displacement
threshold permits reuse. Spatial SFC work is bounded derived policy data, not
particle truth.

Hierarchy extends the existing scheduler/orchestrator using synchronized DMO
blocks, canonical momentum drift and coarse PM endpoint kicks. Sources share
each fine epoch; only target forces/kicks are selective. Split caches have explicit
row/PM epoch checks and are reconstructed at each block start. Publication and
migration require closed all-active coarse endpoints. Gas/source multirate laws
and lazy prediction need a separate scientifically specified extension. Core's
preparation callback accepts generic exception/string_view data; workflow supplies
its existing collective failure coordinator. No dependency-boundary exception.

## CHUI Part 1 foundation, hydro method suite and audit F1/F3/F4/F7 — 20261103

Source implementation; NOT compiled, NOT run (user-reserved validation).

* Hydro numerics are three orthogonal selections resolved once per advance:
  reconstruction (PCM/PLM/PPM, `hydro_reconstruction.hpp`), integration family
  (MUSCL-Hancock/CTU, `HydroNumericsPolicy` owned by `HydroCoreSolver`) and
  Riemann solver (HLLC/HLLE, `makeHydroRiemannSolver`). Exactly one Hancock
  half step per face state; CTU adds transverse conservative increments. The
  core-owned config enums are mapped in workflows so `core` does not include
  `hydro` (no dependency exception).
* `HydroFace` gains `owner_minus2_cell`/`neighbor_plus2_cell` (appended, default
  invalid) and a documented along-normal stencil order. Additive public change;
  existing designated initializers remain valid.
* Dual energy is step-local with total-energy synchronization; no persistent
  auxiliary lane and no restart schema change. The measured non-conservation is
  reported separately from the residual.
* Single thermodynamic refresh contract `hydro::refreshGasThermodynamics()` with
  rank-local, non-persistent mutation/refresh epochs on `SimulationState`;
  timestep selection fails closed on stale caches. `hydro` now depends on
  `core::SimulationState` for this helper (core remains foundational).
* Node-aware admission is a pure per-rank rule (equal share of the node limit
  minus a reserve subtracted once); no collective, so it cannot deadlock.
  Per-owner governor telemetry is a fixed 64-slot table.
* `core::StageTransaction` (reserve/prepare/validate/commit) with one
  collective vote per transaction through `StageCollectiveDecision`; the MPI
  adapter lives in `workflows/internal` so `core` does not depend on `parallel`.
* Fixed-grid periodic boxes evaluate each box face once (upper ghost copy);
  previously the boundary flux was applied twice.
* The effective-ISM sidecar ledger becomes an integrated energy
  (`energy_ledger_version=2`); legacy density sums are discarded on restart.

## CHUI Part 1 corrective pass — 20261103

Source only; NOT compiled or run.

* Thermal authority is the canonical specific internal energy in both hydro
  paths; AMR no longer treats the pressure cache as authority.
* Dual energy redesigned: auxiliary internal-energy lane initialized from the
  canonical `e`, Bryan-type selection, no total-energy overwrite inside the
  solver, complete ledger of the thermal correction; the earlier step-local
  "reset E" variant is removed.
* Reconstruction stencils are completed through the existing ghost shell for
  all boundary kinds (no second ghost-exchange system).
* AMR periodic physical boundaries wrap in-patch; cross-patch periodic
  topology is refused (Part 3) rather than producing zero-flux box faces.
* Node-shared accounting pool (MPI-3 shared window + lock-free counter)
  attached to each rank's existing `MemoryGovernor`; no new allocator, no
  collective in the reservation path. `core` defines the pool interface; the
  MPI backend lives in `parallel`.

## CHUI Part 1 stabilization gate — 20261103

Source only; NOT compiled or run.

* Hydro CFL rule is selected by the integration family (dimension-summed for
  MUSCL-Hancock, per-axis for CTU) through one core mapping used by both the
  TimeCoordinator candidate and the accepted-step guard.
* AMR mutations (scatter, reflux, regrid) write derived thermodynamic caches
  through `amr::AmrThermodynamicCachePolicy` (closure + temperature factor);
  the whole-rank post-AMR closure pass is removed. `mu` is an interface
  parameter (Part 1 reference 0.59).
* `shared_pool` claims its node-wide guarantee only after compile-time and
  runtime conditions are agreed world-wide; otherwise all ranks fall back to
  `shared_equal`. The core no longer static-asserts lock-freedom.
* Dual energy: explicit reflux fallback (keep pre-reflux e when E-K is
  unreliable, ledgered) and ghost auxiliary fill from local sources; remote
  sources counted.
* Periodic cross-patch AMR topology remains refused (mandatory Part 3
  deliverable).
