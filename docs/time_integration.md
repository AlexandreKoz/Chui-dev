# Time Integration Framework

This document defines the baseline time-integration contract for CosmoSim and records conservative assumptions made in the first implementation.

## Stage contract

The authoritative per-step stage order is:

1. `gravity_kick_pre`
2. `drift`
3. `force_refresh`
4. `source_terms` (`kStrangPreHalf`)
5. `hydro_update`
6. `source_terms` (`kStrangPostHalf`)
7. `gravity_kick_post`
8. `analysis_hooks`
9. `output_check`

The two source dispatches define the deterministic reference
`S(dt/2) H(dt) S(dt/2)`. `SourceRuntimeImpl` derives one `SourceStagePlan` per
dispatch from the integrator-owned `CosmologicalStepFactors` (the coarse or
hierarchical-block interval) and the dispatch's `SourceSplitPhase`:

| Plan field | `kStrangPreHalf` | `kStrangPostHalf` | `kUnsplitCompatibility` |
|---|---|---|---|
| smooth sub-interval | `[t_b, t_mid]` | `[t_mid, t_e]` | `[t_b, t_e]` |
| smooth evaluation epoch | a at physical fraction 1/4 | a at 3/4 | `a_end` |
| smooth duration | `dt/2` (code = physical, code time is linear in cosmic time) | `dt/2` | `dt` |
| discrete events | none | once, at `t_e`, `a_end`, full `dt` | once, at `t_e`, `a_end`, full `dt` |
| opens / closes interval | opens | closes | both |

Smooth order is palindromic: pre = element diffusion + SGS shear, then
thermochemistry/ISM/dust/SGS dissipation; post = thermochemistry, then
transport. Discrete events (star birth → SSP/stellar returns, winds/pre-SN, SN
feedback → BH seeding/accretion/feedback) follow the post smooth half. In
proper-coordinate runs every scale factor is exactly 1. The stellar population
sidecar carries both cumulative-Poisson thresholds and an absolute
accepted-interval cursor. Its schema-2 cursor must match the next interval
begin (within round-off); gaps, overlap and replay fail before event ordinals
or physical state commit. Schema-1 records adopt their first post-restart
boundary and are rewritten as schema 2.

**Controller publication.** Thermochemistry never publishes controller state
itself: every half (and every AMR level half) accumulates its shadow/adaptive
observations into `SourceIntervalAccumulator`; `closeSourceInterval()` performs
the single collective `updateGlobalChemistryControl()` per Strang interval.
The accumulator is part of the source runtime's rollback checkpoint, so a
rejected trial (including one rejected between the halves) is replayed from a
closed, empty interval. Opening a new interval while one is open, or a post
half without its pre half, is a `std::logic_error`.

With recursive AMR subcycling, `HydroAmrRuntime::bindSourceRuntime()` enables
the level coupling and the AMR clock owner applies
`S_L(dt_L/2) H_L(dt_L) S_L(dt_L/2)` to thermochemistry/ISM/dust/SGS for every
level subinterval through `RecursiveAmrLevelSourceOperator`
(`RecursiveAmrLevelSourceAdapter` → `SourceRuntime::executeAmrLevelSmoothSources`).
The pre half precedes the level start-history capture; the post half follows
the level sweep and its remote register exchange and precedes the end-history
capture, so children interpolate complete coarse Strang states; reflux and
restriction from finer levels commit after the children. Every rank enters
each level-source transaction, including empty local spans. While the levels
own thermochemistry, the global halves skip it (same predicate on both sides:
coupling enabled, recursive mode, production AMR coverage on any rank) and the
interval close verifies matched level pre/post counts. Element diffusion, SGS
shear transfer and all discrete events stay globally synchronized at the
coarse halves/endpoint (distinct smooth cadence vs. synchronized ownership
transactions); the legacy post-sweep recursive chemistry callback was removed.

The scheduler (`StageScheduler`) owns this ordering and solver modules interact via stage-bound `IntegrationCallback` implementations. Each callback must declare its exact `integrationStages()` set during registration; `StepOrchestrator` stores handlers in per-stage buckets and dispatches only the bucket for the active stage in stable registration order. Production callbacks must treat any off-stage direct invocation as a contract violation rather than silently self-filtering. This keeps gravity/hydro/source/output sequencing explicit and auditable, prevents solvers from receiving irrelevant stages, and avoids per-stage broadcast scans in the hot path.

### Callback migration note

Any new `IntegrationCallback` implementation must provide:

- `callbackName()` for profiling and diagnostics;
- `integrationStages()` returning the small typed stage set the handler is allowed to receive;
- `stageContracts()` returning a typed `StageContract` entry for each declared stage (required inputs, mutable domains, outputs, side effects, sync requirement, active-set family, restart/output safety, and owning subsystem);
- `onStage(...)` logic that assumes the declared stage contract and throws/asserts on impossible off-stage direct calls instead of returning silently.

Existing code that previously registered a broad callback and checked `context.stage` internally must split the behavior into explicit stage declarations or return the exact small stage set needed by the handler. No config keys, snapshot/restart payloads, or solver numerics migrate for this interface repair.

### TreePM long-range force-refresh contract (distributed TreePM runtime)

The reference workflow applies an integrator-owned, config-driven TreePM
operator inside KDK. An invalid initial force cache is bootstrapped at the
globally authorized `gravity_kick_pre` boundary. After drift, the
`force_refresh` stage evaluates the force at the end-of-step source epoch; the
`gravity_kick_post` stage consumes that cache for the closing half-kick.

The synchronization decision is explicit in workflow runtime events and
cadence records (`pm_sync_surface`, `gravity_kick_opportunity`,
`field_version`, `last_refresh_opportunity`, `active_particles_kicked`,
`inactive_particles_skipped`).

For the in-flight refresh, the core integrator first records local evidence in
`PmRefreshDirective` without advancing `PmSynchronizationState`. In distributed
execution `TimeCoordinator` collectively converts that evidence into one
rank-global directive before any branch-dependent TreePM work; in serial
execution the same directive is materialized locally. That materialized
directive is the transient decision authority until the solver succeeds and
`PmSynchronizationState` commits it. The `gravity.pm_long_range_field` event is
emitted in that interval, so its opportunity, version, build step, and build
scale-factor payload must come from the directive/solver decision. Reading the
still-previous committed cadence state there would produce a
self-contradictory refresh event. This does not move committed cadence ownership
away from `PmSynchronizationState`.

#### Production refresh and lower-level reuse rule

- `numerics.treepm_update_cadence_steps` is normalized and validated as exactly
  `1` for the production `ReferenceWorkflow`.
- Every authorized, rank-coordinated production force-refresh evaluates a new
  PM field. All ranks enter the decision and collectives, including ranks with
  zero local targets.
- Lower-level scheduler and `TreePmCoordinator` interfaces retain explicit
  refresh/reuse support for focused tests and future multirate work. Reuse is
  legal only when a compatible cached field is present on every rank. A
  missing or incompatible cache throws coherently; it is never converted into
  an implicit refresh.
- Cache compatibility covers force epoch, force-evaluation scale-factor
  metadata, `G_code`, split and box geometry, assignment, boundary,
  decomposition-layout mode, and window-deconvolution policy. Particle
  ownership/decomposition epoch is intentionally excluded because the cached
  PM field is owned by the fixed FFT slab layout, not particle rows.
- Refresh/reuse votes and cadence metadata are rank-consensus facts. Rank-local
  field-validity/source-generation evidence is reduced at the workflow
  coordination boundary before the gravity runtime may branch on PM cadence.

The scale factor stored with a field is source-epoch/validity metadata. The PM
and tree kernels are scale-free; cosmological time dependence is applied only
by the KDK factors described below.

#### KDK operator order with PM sync

For each global step:

1. `K_PM+Tree^pre` on the active set, bootstrapping force truth if necessary
2. `D` drift on active set
3. rank-coordinated `force_refresh` at the post-drift source epoch
4. hydro/source stages
5. `K_PM+Tree^post` on the active set using the refreshed force cache
6. analysis/output stages

The integrator owns refresh authorization and commits the PM synchronization
event only after the TreePM callback consumes it successfully.

#### Inactive-particle treatment

- The default workflow uses `hierarchical_max_rung=0` and all-active global KDK.
- P6 adds opt-in collisionless DMO block KDK (`1..12`), described below. Only
  active targets receive Tree forces/kicks; all source positions share the
  current fine drift epoch. This source path is not numerically qualified.

#### Restart continuation rule

- Restart payload persists cadence metadata (`gravity_kick_opportunity`,
  cadence steps, field version, last refresh opportunity, and field-build
  step/scale factor).
- Restart policy is `deterministic_rebuild`: PM mesh/tree scratch is not
  checkpointed and cannot be treated as reusable field truth. Restart schema
  v20 does persist the committed gravity-force history keyed by stable particle
  ID for exact KDK continuation; import remaps compatible rows and invalidates
  stale or migrated rows. The next integrator-authorized force surface
  reconstructs transient PM/tree state.
- Cadence metadata remains auditable across restart boundaries via distributed restart state and provenance fields.


## Stage 2 timestep authority contract

Stage 2 uses a single-owner timestep model: `HierarchicalTimeBinScheduler` is the only live authority for per-element bin assignment, next activation, active flags, pending transitions, active-set construction, and PM kick cadence metadata. Solver callbacks may propose timestep candidates and consume scheduler-built active sets, but they must not treat `ParticleSoa::time_bin`, `CellSoa::time_bin`, migration records, or restart mirrors as authority.

The reference path remains `numerics.hierarchical_max_rung=0`. P6 optionally
uses the same scheduler authority for synchronized DMO blocks, with explicit
per-row Tree force epochs and common source drift epochs. It remains
SOURCE-IMPLEMENTED / VALIDATION PENDING; gas/source multirate integration is
rejected rather than assigned an unspecified operator.

### Rung-zero global physical timestep

Rung zero does **not** mean fixed timestep. Before every KDK construction the
production `TimeCoordinator` recomputes cosmology, gravity, hydro-CFL, and
active source limits, coordinates local preparation failure, and reduces the
rank-local minimum with MPI `MIN`. `+infinity` is the explicit no-limit sentinel
for an inactive criterion; NaN, negative infinity, zero, and negative proposals
fail closed. An explicit `ReferenceWorkflowOptions::dt_time_code` is a maximum
step for testing/embedding and cannot bypass a tighter physical or endpoint
limit. The accepted value is finally clipped to the next ordered output boundary
and run endpoint before stage times are constructed. `max_global_steps` only
limits how many committed steps a segment may execute.

For cosmological scale-factor authority, every committed step additionally
requires finite positive monotonic `a`, finite redshift, and
`delta_ln(a) <= cosmology_max_delta_ln_a` within the documented floating-point
tolerance. The terminal code-time interval is FLRW-derived from `a_end`, so the
last timestep is endpoint-clipped instead of allowed to overshoot into the
future.

### Scheduler authority and mirror policy

- Authoritative live lanes: `HierarchicalTimeBinScheduler` hot metadata (`bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`) and `PmSynchronizationState` cadence fields (`gravity_kick_opportunity`, cadence steps, field version, last refresh opportunity, field-built step/scale factor).
- Derived mirrors: `ParticleSoa::time_bin`, `CellSoa::time_bin`, migration/transfer `time_bin` payloads, and serialized restart state mirrors are diagnostic or compatibility mirrors only. They are refreshed from scheduler state with `syncTimeBinMirrorsFromScheduler(...)` and may be validated for corruption, but they are never fallback scheduling truth.
- Cell mirrors in production gas-cell states map through dense local gas-cell rows validated by `GasCellIdentityMap`.
  Legacy particle-bound import paths may still map through a parent gas particle, but production scheduler invariants
  do not require local gas-particle count to equal local gas-cell count.
- Active-set descriptors for hierarchical execution must be constructed from scheduler output (`makeSchedulerActiveSetDescriptor(...)`) and must carry scheduler tick/generation provenance before solver callbacks consume them.

### Candidate criteria flow

1. Physics modules compute local candidate timesteps through typed criteria hooks such as CFL, gravity acceleration, source-term, or user-clamp hooks.
2. Candidates are submitted to the scheduler with an auditable `TimeStepCandidateSource`. Optional human-readable labels are transient `std::string_view` inputs used only at the diagnostic boundary; they are not retained per element. Production scheduler state stores compact typed provenance, so diagnostic label length cannot change population-scale retained capacity.
3. `mapDtToTimeBin(...)` normalizes physical `dt` proposals into the integer power-of-two bin hierarchy using the configured `TimeStepLimits`.
4. `reconcileCandidateTransitions()` conservatively keeps the finest submitted bin per element, validates synchronization legality, and records clipping/candidate counters.
5. Pending transitions are committed only at scheduler-controlled substep boundaries; invalid destination synchronization fails fast rather than being silently clipped or delegated to mirrors.

### Active-set construction flow

1. `beginSubstep()` validates scheduler internals and rebuilds the compact active list from scheduler `bin_index` and `next_activation_tick`.
2. `RungZeroTimeState` and `TimeCoordinator` retain the scheduler-owned
   particle and gas-cell active lists as read-only spans for the open substep.
3. `makeSchedulerActiveSetDescriptor(...)` stamps scheduler tick and state-generation provenance and validates that descriptor indices still match scheduler activity.
4. `StepOrchestrator::executeSingleStepWithDispatcher(...)` runs the canonical
   KDK stage order while `TimeCoordinator` rebuilds the appropriate typed
   resource view for each stage. Existing callback-based core entry points
   remain available for lower-level compatibility tests, not production
   workflow composition.
5. `endSubstep()` reconciles candidates, applies legal pending transitions, clears active flags, advances the integer tick, and only then are public `time_bin` mirrors refreshed from scheduler truth.

### Invariant framework

The Stage 2 invariant framework intentionally traps split-brain timestep ownership before solver callbacks can consume it:

- scheduler metadata arrays must remain same-sized, in-range, and internally consistent;
- active flags are scheduler caches for an open substep, not serialized live authority;
- active-set descriptors must match scheduler active elements and source generations;
- restart import must validate scheduler lanes and reject stale particle or cell `time_bin` mirrors before rebuilding mirrors from scheduler state;
- PM refresh opportunities may be previewed without mutation for distributed
  coordination, but an accepted opportunity must be committed after successful
  solver execution, and a refresh must be committed before the next kick
  opportunity can be accepted;
- distributed PM decisions remain rank-consensus metadata;
- production restart import requires particle and gas schedulers to use
  `max_bin = 0`, with all current/pending bins at zero or unset as appropriate;
- current evidence does not make mixed-rung TreePM synchronization
  production-proven.

The fail-closed rung-zero policy prevents the scheduler scaffold from being
mistaken for production multirate integration.

## Hierarchical integer timeline bins (infrastructure)

`HierarchicalTimeBinScheduler` implements power-of-two integer bins (`dt_bin = dt_min * 2^bin`) and maintains:

- hot sidecars (`TimeBinHotMetadata`) for `bin_index`, `next_activation_tick`, and `active_flag`
- bin-local compact ownership lists (`m_elements_by_bin`) with O(1) move/erase sidecars
- per-substep compact active lists sorted for deterministic memory access

This avoids repeated full-state scans: each substep only touches bins synchronized to the current integer tick.

## State and active sets

`IntegratorState` tracks:

- `current_time_code`
- `current_scale_factor`
- `dt_time_code`
- `step_index`
- `TimeBinContext` (hierarchical-bin scaffold)

`ActiveSetDescriptor` carries compact particle/cell subsets and subset flags, avoiding global full-state sweeps in scheduler logic.

## Timestep criteria hooks

Typed helper criteria are provided and can be combined via `TimeStepCriteriaRegistry`:

- CFL limiter: `dt_CFL <= C_cfl * (Delta x / (|u| + c_s))`
- directional hydro CFL limiter:
  `dt_hydro_CFL = min_axes(C_cfl * cell_width_axis_code / (abs(v_axis_code) + c_s))`
- Gravity limiter for comoving softening and scale-free TreePM `A`:
  `dt_grav = eta * sqrt(a^3 eps_com / |A|)`
- source-term limiter and user clamp hooks

`computeGravityTimeStep(...)` itself remains coordinate-system neutral and
evaluates `eta sqrt(length/acceleration)`. Cosmological callers should instead
use the public `ComovingGravityTimeStepInput` and
`computeComovingGravityTimeStep(...)` interface. It accepts `eps_com`, the
scale-free TreePM magnitude `|A|`, and `a`, validates a finite positive scale
factor, converts to comoving-coordinate acceleration `|A|/a^3`, and delegates
the final length/acceleration calculation to the generic helper. The production
workflow uses the committed `IntegratorState.current_scale_factor` because
adaptive-bin proposals run after step commit. Using the peculiar-velocity
response `A/a^2` directly with a comoving length would mix coordinate systems.

### Public API migration: cosmological gravity timestep

Code that previously constructed `GravityTimeStepInput` from a comoving
softening and a TreePM force lane must migrate to
`ComovingGravityTimeStepInput` from
`include/cosmosim/core/time_integration.hpp`. Pass the unscaled TreePM
`|A|` and the force-epoch scale factor; do not pre-divide by `a^2` or `a^3`.
The existing `GravityTimeStepInput`/`computeGravityTimeStep(...)` API remains
available for callers whose length and acceleration are already expressed in
one coordinate system.

`combineTimeStepCriteria` takes the conservative minimum across registered hooks and a fallback dt.

Hydro CFL candidates submitted by the reference workflow use local Cartesian patch widths where gas-cell
center metadata identifies a row-ordered patch, otherwise the same near-cubic fixed-grid fallback used by the
hydro patch builder. The scheduler remains the only authority for the accepted bin; hydro submits
`TimeStepCandidateSource::kHydroCfl` candidates and the callback verifies the accepted `dt_time_code` against the
current directional CFL bound before invoking the finite-volume update. A violation is a hard runtime error, not a
silent clamp inside the hydro solver.

`HydroCflDiagnostics` records the worst checked gas cell with local row, gas-cell ID when available, patch ID/row
when available, proposed and accepted `dt`, CFL number, safety factor, velocity components, sound speed, cell widths,
and limiting axis. The reference workflow emits this through the profiler event `hydro.cfl_guard` for reporting.

## Synchronization and legal transitions

Elements may request bin promotion/demotion through the scheduler candidate
path. Refinement is committed at the active safe boundary. Coarsening is
committed only when the current tick is synchronized to the larger destination
period; an unaligned coarsening request remains persisted in
`pending_bin_index` while the element continues on its old bin. The scheduler
records each such wait in `deferred_coarsening_events`. An unexpectedly
unaligned refinement is an invariant failure and is recorded in
`illegal_transition_attempts` before the scheduler throws.

## Diagnostics and pathological collapse visibility

`TimeBinDiagnostics` reports:

- occupancy and active counts per bin
- active fraction
- promoted/demoted counts
- clipped-to-min/max counters from dt mapping
- illegal transition attempts
- deferred coarsening events
- collapse candidate count when finest-bin occupancy dominates

## Cosmology helper equations

The helper API uses physical peculiar velocity `u` and comoving position `x`:

- `dx/dt = u / a`
- `du/dt + H(a) u = A / a^2 + source_terms`
- `da/dt = a H(a)`

Here `A = -grad_x Psi` is the scale-free comoving Newtonian kernel returned by
both PM and the complementary tree residual. Neither gravity branch inserts an
additional factor of `a^2`. Collisionless particles receive the response
through KDK Hubble-drag/kick factors. Gas receives the same `A/a^2` response in
`ComovingGravityExpansionSource`, which also applies the conservative
expansion sources using the integrator-owned `a` and `H` in code units.

For hydro conserved variables `m=rho_com u`,
`E=rho_com(e+|u|^2/2)`, and kinetic energy density `K`, the explicit source is

```text
S_m = rho_com A/a^2 - H m
S_E = rho_com (u dot A)/a^2 - H (2 K + 3 P_com).
```

Because hydro runs after drift/force refresh but before the step commit,
fixed-grid and AMR callbacks consume
`StepContext.timeline_step.{scale_factor_end,hubble_end_code}`. The committed
`IntegratorState` still represents the step-begin epoch at that point; solver
callbacks must not use it or independently recompute an SI-valued Hubble rate.

The implementation provides midpoint-integrated drift/kick prefactors over scale factor and a closed-form Hubble drag factor `a_begin / a_end`.

## Assumptions (explicit)

- Baseline stage scheme is kick-drift-kick only.
- Production PM refresh cadence is exactly one at every rank-coordinated force
  evaluation. Lower-level reuse is explicit and fail-closed, but is not a
  production asynchronous multirate PM integrator.
- Time bins use strict power-of-two integer periods in tick space.
- Production config currently permits only rung zero despite the lower-level
  multi-bin scheduler implementation.
- Scale factor update in the orchestrator uses forward Euler (`advanceScaleFactorEuler`) for conservative simplicity.
- If no cosmology background is passed to `StepOrchestrator`, the orchestrator advances `time` only and leaves `current_scale_factor` unchanged.
- Output cadence policy remains external to the core integrator; output triggering is represented by the explicit `output_check` stage callback.

## Provenance and schema implications

No snapshot schema field is changed by this gravity repair. Restart schema
version 20 already persists scheduler/cadence truth needed for deterministic
continuation plus committed stable-ID-keyed force history. PM mesh/tree scratch
and coordinator-local cache structures remain transient and are
deterministically rebuilt.

## P6 synchronized hierarchical DMO KDK

Status: **SOURCE-IMPLEMENTED / VALIDATION PENDING**. Global KDK remains default.
`StepOrchestrator::executeHierarchicalBlockWithDispatcher` adds an operator to
existing core orchestration; `TimeCoordinator` retains typed stage dispatch,
collective failure, output and migration ownership. Distributed callers supply
the existing failure coordinator through the core preparation callback seam.
Core has no workflow/MPI dependency.

For the existing force convention,

\[
\dot x=u/a,\quad \dot u+Hu=A/a^2,\quad p=au,
\quad \dot x=p/a^2,\quad \dot p=A/a.
\]

`A` is scale-free Tree+PM acceleration. Hierarchy stores `u=p/a` at the common
position epoch. Its drift operator holds `p` fixed:

\[
x_1=x_0+a_0u_0\int_{t_0}^{t_1}dt/a^2,\qquad u_1=u_0a_0/a_1.
\]

A force endpoint kick adds
`delta_u = A_endpoint * integral(dt/a) / a_endpoint` over the corresponding
half interval. The existing `computeComovingDriftFactor` supplies the kick
integral; the canonical position-drift integral uses the same 64-sample midpoint
convention in scale factor, `integral da/(a^3 H)`. SI/code time conversion remains
owned by `CosmologicalTimeline`. The reference global operator is retained with
its existing factors; this optional canonical formulation needs convergence
comparison, including free expansion and small-force cases.

At a closed coarse boundary, every particle receives a freshly synchronized
Tree+PM solve. The coarse duration is at most `2^M * global_min_gravity_dt`,
further clipped by cosmological limits, explicit workflow dt, next code-time
output event, and endpoint. The existing gravity criterion uses
`eta*sqrt(eps/(|A|/a^3))` with `eta = numerics.gravity_timestep_eta` (default 0.2). The quantum is coarse dt /
`2^M`. Each bin `b` is the largest permitted power of two no greater than its
particle criterion, with interval `quantum*2^b`, `b=0..M`, `M<=12`.
Assignments and quantum are frozen inside the block. Changes occur only after
all intervals close, so no outstanding old interval is reinterpreted.

The existing scheduler candidate lane holds proposals. Exact bin-list, active
and radix-scratch capacities are admitted before membership changes. Closure
clears active flags and sets next activation to the current aligned tick without
advancing it. Gas's empty scheduler shares the integer tick. No second scheduler
or dense particle map is added. Rung-zero uses its existing uniform physical
representation when the hierarchy is disabled. Compact DMO species/ownership/
drift metadata remain scalar; the optional one-byte bin lane is a mirror.

Each fine tick audits all nine canonical dispatches. Pre-kick uses the current
scheduler subset and exact cached Tree force from that endpoint. Drift advances
**all owned DMO sources**, updates the common drift epoch, and advances source
generation once. Scheduler subsets then select closing force targets. Their
Tree forces are refreshed before closing half-kicks. The same subset receives
its next opening half-kick on the next iteration. Bins never change mid-interval.
This first architecture deliberately pays O(N) source drift per fine tick; lazy
inactive-source prediction is not implemented.

PM follows a separate coarse endpoint kick operator. All particles receive
opening/closing PM half-kicks only at coarse surfaces. Intermediate Tree solves
use `short_range_only` and do not solve or interpolate PM. This is a defined
split operator, not extrapolation of a stale mesh. The existing cadence key
still equals 1, now counting the explicit coarse PM opportunities. Mandatory
synchronized bootstrap occurs at every block start in uninterrupted and restart
runs; together with the closing endpoint it costs two PM solves per coarse
block. Eliminating repeated endpoint bootstrap needs a restart-equivalent
component-cache contract and is deferred.

Tree-only and PM-only forces are accumulated separately, without total-minus-PM
subtraction. Each row's Tree generation must exactly equal current source
generation for a kick. Nonzero PM kick factors additionally require synchronized
PM source generation and field version. Previous-acceleration history is only
an accuracy scale: within a block it is current Tree plus the coarse PM
component, not a claim that the PM field is valid at a fine source epoch. History
must have stable row identity and lie in the current coarse generation interval;
incompatible/unavailable history falls back geometrically. The default DMO
single-drift history path is preserved independently and never legitimizes a
stale kick.

At final tick all rows are active. Closing kicks complete, scheduler intervals
close, and the integrator commits the coarse endpoint once. `step_index` counts
coarse blocks; scheduler ticks count fine drifts. Cosmological endpoint checks
precede migration/publication. Only then may the existing bounded migration
run, followed by analysis, snapshot and restart. Migration invalidates row/cache
identity and reinstalls geometry; the next mandatory bootstrap reconstructs
split forces. Fine output/checkpoint is forbidden. No mixed-epoch snapshot or
in-flight restart schema is introduced.

The bounded scale-factor table contains at most 8193 doubles (65544 bytes),
admitted through the existing process governor. New component caches cost
56 bytes/local particle plus the optional bin mirror and existing generic
scheduler lanes. Restart v23 already carries bins, next activation, common
integration/drift epochs, PM cadence, source generation and synchronized total
force history. Transient split forces, meshes, trees and the table are rebuilt.
Exact restart equivalence and multirate convergence remain required tests.

## Hierarchical recovery safeguards (2026-10-10; provisional)

The existing collisionless block dispatcher remains the only hierarchical KDK
implementation. It now rejects incomplete/duplicated all-source drift identity
rows and a mismatched common position epoch before opening a block, and supplies
its scheduler pointers through the direct dispatcher seam. The in-flight marker
is set only after block-start active-set checks. The workflow retains the fresh
synchronized bootstrap in both continued and resumed blocks: removing it would
require a separately established split-cache/history restart equivalence contract.

An independent coarse PM resolution-crossing estimate augments existing
softening/gravity, expansion, output, endpoint and explicit-dt restrictions.
At synchronized scale `a`, let `v=max_rank,particle |u|/a`,
`g=max_rank,particle |A|/a^3`, and
`ell=min(dx_mesh,dy_mesh,dz_mesh,Gaussian_split_length)`.
Require `v*dt + g*dt^2/2 <= ell`, evaluated stably as
`dt=2*ell/(v+hypot(v,sqrt(2*g*ell)))`. Both zero returns infinity.
The helper uses existing comoving/peculiar conventions; existing softening
criteria still select rungs. Rank maxima use two unconditional reductions before
any solver dispatch. The frozen-epoch estimate limits PM scale crossings; it is
not a rigorous future-force-growth bound. Expansion limits remain independent.
This changes optional hierarchical coarse intervals and therefore trajectories,
scheduler ticks and output cadence; global rung-zero behavior is untouched.

### Coarse-interval criterion selection and diagnostics (October 2026)

The `max_speed_mesh` estimate above remains the default
(`numerics.hierarchical_coarse_displacement_criterion`). Because it combines
the box-wide maximum speed with the box-wide maximum total acceleration over
one PM cell, a single fast halo particle bounds every coarse block; at fine
PM meshes this made a rung-2 coarse block cover only about 1.8 global KDK
steps, so the four fine ticks were extra work rather than savings.

`rms_split_scale` instead bounds the long-range interval by the typical
motion relative to the scale on which the PM field varies. The PM force is
the Gaussian-filtered field `F_LR(k)=exp(-k^2 r_s^2)`, smooth on `r_s`, so the
coarse interval requires

\[
v_{\rm rms}\,\Delta t_{\rm coarse}\le f\,\min(r_s,\bar d),\qquad
v_{\rm rms}^2=\frac{\sum_i m_i |u_i|^2}{a^2\sum_i m_i},\qquad
\bar d=\left(\frac{L_xL_yL_z}{N}\right)^{1/3},
\]

with `f = numerics.hierarchical_coarse_rms_displacement_factor` (default
0.25). Fast particles in dense regions are still restricted by their own
softened-acceleration rung criterion for the short-range force. Mass sums,
mass-weighted squared speeds and particle counts use three unconditional
`allreduceSumDouble` reductions in both modes. This is the conventional
GADGET-2 style long-range step rule, implemented from the published formula.

Every `time.hierarchical_block` event now records `coarse_limiter` (one of
`particle_gravity`, `pm_displacement`, `cosmology_expansion`, `explicit_dt`,
`output_or_end_boundary`), the candidate intervals
(`particle_coarse_limit_time_code`, `pm_displacement_limit_time_code`,
`expansion_limit_time_code`), `coarse_displacement_criterion`, and
`rank_local_bin_occupancy` (`n0,n1,...` committed rung histogram). The rung-2
ceiling remains: a coarse block can never exceed `2^max_rung` times the
smallest particle interval, so useful savings require a larger
`hierarchical_max_rung` once the coarse limiter no longer binds.

Timestep growth limiter: a particle's hierarchical interval may at most
double from one synchronization point to the next (refinement is
unrestricted). The acceleration-only criterion is blind when a softened pair
passes through coincidence (|a| -> 0); without the limit such a pair jumped to
the coarsest rung for a whole block. The previous interval is reconstructed
from restart-persisted state only (previous quantum `integrator.dt_time_code`
times the committed rung period), so resumed and uninterrupted runs apply the
same limit; there is no limit before the first completed step.

Softened-encounter criterion (`numerics.hierarchical_encounter_timestep_enabled`,
default `true`). The growth limiter alone left the D1 pair-crossing error at
1.62e-3 against the 1e-3 gate, because the acceleration criterion still
under-resolves the crossing itself. Each synchronized Tree solve now records,
per target, the squared separation `r_nn^2` of its nearest distinct
direct-summed source inside `r_cut` (min-combined over local, group-walk and
remote MPI contributions; response wire v2). The rung candidate is

```
dt = min( eta sqrt(eps a^3 / |A|),                       acceleration
          eta sqrt(a^3 (r_nn^2 + eps^2)^{3/2} / (G m_max)) )  encounter
```

with `eta = numerics.gravity_timestep_eta` (default 0.2), comoving scale-free `A` and `G m/(r^2+eps^2)^{3/2}` (physical
values divide by `a^3`), and `m_max` the global maximum particle mass
(conservative for unequal masses). At coincidence the encounter term is the
softened two-body dynamical time `eta sqrt(eps^3 a^3/(G m))`; beyond
`r ~ eps` it is weaker than the acceleration criterion by `sqrt(r/eps)`, so
it binds only during close encounters. Assumptions: Plummer-equivalent
softening scale `eps`; the nearest neighbour is a direct pair (nodes accepted
as multipoles are not treated as neighbours, which the relative MAC makes
safe at `r <~ eps`); the group-walk path cannot distinguish an exactly
coincident distinct source from the target and omits `r = 0`. The lane is
restart-persisted with the v24 endpoint. Global KDK is unchanged. Validation
target: `integration_hierarchical_treepm_workflow` refined trajectory error
`< 1e-3` against global KDK with no change to that gate.

Endpoint force reuse: a block opens with the previous block's all-active
endpoint force when every rank still holds it complete at the same positions
(counter `time.hierarchical_endpoint_force_reuse`); otherwise the synchronized
bootstrap solve runs (`time.hierarchical_endpoint_force_bootstrap`). This saves
one all-target Tree+PM solve per coarse block. Restart equivalence is kept by
schema v24 (see `docs/restart_checkpointing.md`).

The synchronized bootstrap can grow gravity-owned retained capacity; the
workflow now reconciles the MemoryGovernor baseline collectively immediately
after it and before the rung-membership `RetainedCapacityTransaction` (repair
register B-10). The transaction invariant itself is unchanged.

The production-workflow restart/refinement test and Zel'dovich fundamental-mode
qualification are registered for dependency-complete builds. The growth gate
compares global KDK with two hierarchical temporal resolutions, retains the
existing linear-response envelope, and targets only `k=2*pi/L` (six Cartesian
modes of the 4^3 fixture). Broader growth, power, nonlinear orbit stability and
large-rung qualification remain pending. `hierarchical_max_rung=0` and all
provisional capability classifications are retained.

## Multirate outer PM level — numerical contract (2026-10-26 stabilization; IMPLEMENTED — NOT COMPILED — NOT TESTED — OWNER VALIDATION REQUIRED)

`numerics.treepm_pm_reuse_policy = multirate_outer_kick` is **experimental and
scientifically unqualified**. Hierarchical KDK only. `treepm_update_cadence_steps = M`
in `[2, 64]` is the *target* number of coarse blocks per PM cycle. The scientific
reference remains `exact` (PM at both endpoints of every coarse block).
`stale_within_cadence` is an unrelated global-KDK field-reuse approximation.

The stabilization round replaced the original fixed-interval design. That
design froze the coarse interval for a whole cycle, so it could run a block
whose quantum was coarser than the current particle criterion, and it let
the external step budget truncate the cycle.

### 1. Phase-space variables and operators

The hierarchical integrator uses the convention above: `ẋ = p/a²`,
`ṗ = A/a`. It stores `u = p/a` at the common position epoch, and `A` is the
scale-free comoving acceleration.

| Operator | Modifies | Law (code) |
| --- | --- | --- |
| Drift `D[t0,t1]` | `x` (all sources, every fine tick) and `u` | `x += a0 u ∫dt/a²`, `u ← u a0/a1` (p fixed) |
| Tree kick `K_T` | `u` of active rungs | `u += A_tree ∫dt/a / a_now` over the rung half interval |
| PM kick `K_PM[t_a,t_b]` | `u` of **all** particles | `u += A_PM(x) ∫_{t_a}^{t_b} dt/a / a_now` |

`∫dt/a` is `computeComovingDriftFactor(a_a, a_b, 64)/time_si_per_code`, the
midpoint rule in `a` of `1/(a² H)`. `a_now` is the scale factor at which `u`
is stored when the kick is applied. The PM kick modifies momenta only. Tree
kicks and drifts never use `A_PM`.

### 2. Composition actually executed

One PM cycle opening at `t0`, of length `T`, with end `tE = t0 + T`:

```
K_PM[t0, t0+T/2] ∘ B_1 ∘ B_2 ∘ … ∘ B_m ∘ K_PM[t0+T/2, tE]
```

* `A_PM` in the opening kick is evaluated at the cycle-opening positions. This
  is the endpoint of the previous cycle (endpoint reuse, restart v24) or a
  synchronized bootstrap.
* `A_PM` in the closing kick comes from a PM refresh at the closing
  positions.
* Each `B_j` is the reference power-of-two block without its PM surfaces:
  tree half-kicks per rung around all-source drifts.

In code, `StepOrchestrator::executeHierarchicalBlockWithDispatcher` receives a
`HierarchicalPmCyclePlan`:

* `opens_cycle` attaches the opening PM factor at offset 0 of the first block.
* `closes_cycle` turns on the PM refresh, the closing PM factor, the sync
  event and restart/output safety at the final tick of the last block.
* Interior force evaluations use `short_range_only`.

The opening factor uses `a(t0+T/2)` from
`advanceScaleFactorByCosmicTime(a0, T/2)`. The closing factor integrates from
that same `a_mid` to the *actual* closing scale factor, so the two PM
half-intervals cover `[a0, a_end]` with no gap or overlap.

### 3. Time coordinate and scheduling

Scheduling uses integrator code time `t`. `T` is fixed when the cycle opens,
from physical criteria only:

```
T = min( M · Δc_open,  Δt_PMdisp(open),  t_boundary − t0,  [retry bound] )
```

where `Δc_open` is the reference per-block interval:

```
Δc_open = min( 2^maxBin · min_i dt_i(particle criteria incl. encounter and growth limiter),
               expansion limit,
               explicit dt )
```

`t_boundary` is the next code-time snapshot event or `t_end`, and a cycle
clipped by it ends exactly on it. The execution step budget does **not**
enter `T` (section 8).

Inner blocks are re-planned at every block by the reference rules:
`Δc_j = min(particle·2^maxBin, expansion, explicit)`, clipped to the
remaining cycle time. A remainder shorter than a quarter block is split in
two to avoid slivers. The blocks therefore tile `[t0, tE]` exactly, and the
number of blocks `m` is determined by the dynamics (≈ M when the criteria
are steady).

Endpoint rounding: the driver ends a block at `t + quantum·2^maxBin`, and
splitting by a power of two is exact. The closing block of a cycle still
nudges its interval by at most four ulps so that `fl(t + Δc)` equals the
planned cycle end exactly. Without this, `fl(t + fl(tE − t))` could miss
`tE` by an ulp and leave an unrepresentable sliver before `t_end`. If the
nudge fails, counter `time.hierarchical_multirate_cycle_end_inexact` is
incremented.

The reference (non-multirate) path keeps its unchanged end-boundary clip.

**Scheduler-preparation caches (2026-10-27; all hierarchical segments,
reference and multirate).** Two values are cached. They hold only derived
data, a hit costs O(1) key comparisons, and the results are bitwise identical
to the uncached path. The test seam
`ReferenceWorkflowOptions::hierarchical_test_disable_preparation_caches`
provides the uncached reference.

* **Identity row list `0..N-1`.** This is the all-active index lane, held in
  the segment workspace's scratch arena.
  * Key: the row count, `particleIndexGeneration()`, and the scratch span
    (pointer and extent).
  * Invalidated by reorder, resize, compaction, migration commit, population
    replacement, or any foreign re-preparation of the scratch.
  * Not invalidated by drifts, kicks or rollbacks, because they do not change
    row identity.
  * On a miss the arena is reset before rebuilding. This also fixes a
    defect: every block used to allocate a new `4N`-byte list in a
    never-reset arena.
* **Local maximum particle mass** (encounter criterion).
  * Key: the row count, `SimulationState::particleMassGeneration()`, and the
    mass-lane address.
  * `particleMassGeneration()` advances with every dense-row change and with
    every in-place mass write. Writers must call
    `notifyParticleMassMutation()`; the star formation, black hole, tracer,
    gas-parent, ghost, view-scatter, IC-stream and Python-binding writers do.
  * The global maximum is still reduced every block, so the collective count
    is unchanged.

Counters (accepted and rejected blocks alike, i.e. execution cost):
`hierarchical_identity_row_{rebuilds,reuses}` and
`hierarchical_mass_max_{scans,reuses}`.

### 4. Timestep invariants (enforced)

1. **Particle criteria are never violated.** Every block's quantum is
   `Δc_j / 2^maxBin ≤ min_i dt_i`, the same guarantee as the reference.
   `HierarchicalBlockRecord` reports both the quantum and the global minimum
   criterion so tests can assert it.
2. **No unrepresentable intervals.** A non-finite, non-positive or
   unrepresentable quantum or cycle length throws on every rank before any
   state changes (`time.hierarchical_quantum_unrepresentable`,
   `…_recovery_failures`).
3. **Rung transitions** only happen at block synchronization points, as in
   the reference. Rungs are fixed within a block, so no KDK half-kick history
   is reinterpreted. Changing the block interval between blocks changes only
   the next block's quantum.
4. **The PM interval honours the PM displacement criterion.** It is checked
   at the opening and **re-checked at every interior block** with the
   current velocities. If `Δt_PMdisp(now) < T`, the cycle is rejected
   (section 5).
5. **Cosmological factors** are computed from the actual block (tree kicks,
   drifts) and the actual cycle endpoints (PM kicks).

### 5. Recovery: transactional rollback and retry

The full state inventory (each element's class and its restore or
invalidation rule) is in
[the round-4 completion record](repair/multirate_round4_completion_20261026.md).
This section is the contract.

**Transaction boundaries.** The coordinator agrees each transition
collectively on every rank:

* **Start.** At a multi-block cycle opening, after the endpoint force is
  available and before rung assignment and the opening PM kick:
  1. Snapshot admission is reduced with `allreduce(min)`. A governor refusal
     or a `std::bad_alloc` while copying counts as a refusal.
  2. If every rank admitted, every rank captures the rollback point and calls
     `GravityRuntime::beginCycleTransaction`.
  3. Otherwise the memory policy applies (below).
* **Accept.** At the closing block, `commitCycleTransaction` publishes the
  diagnostics journal and the rollback point is disarmed.
* **Reject and retry.** At an interior boundary, the globally reduced PM
  displacement bound (or the test seam) rejects the cycle. Every rank runs
  the rollback and reopens with `T ← min(T/2, Δt_PMdisp(now))`; the seam
  reopens with the same `T`.
* **Unrecoverable.** The cycle fails collectively on any of these:
  * a rank-local exception in snapshot, begin, rollback or reconcile (through
    `FailureCoordinator::rethrowCollectiveFailure`);
  * a fingerprint mismatch;
  * more than 64 retries;
  * an unrepresentable `T` or quantum;
  * a refusal under policy `fail`.

  No rank is left inside a collective.

**Rollback point.** The rollback point is a segment-lifetime governed store
(`time.multirate_cycle_rollback_snapshot`). Its memory model (2026-10-27) is
derived from the buffers' element types:

| Buffer | Element | Rows | Bytes/row |
| --- | --- | --- | --- |
| positions x/y/z, velocities x/y/z | `double` | particles | 48 |
| SFC key | `uint64` | particles (when materialized) | 8 |
| force endpoint: total, Tree and PM x/y/z, nearest pair r² | `double` | particles | 80 |
| particle scheduler: bin, next activation, active, pending | `u8`, `u64`, `u8`, `u8` | particles | 11 |
| cell scheduler (same lanes) | — | gas cells (0 for DMO) | 11 |

That is **147 B per particle** (139 B without a materialized SFC key), plus a
fixed part: the `StateMetadata` string capacities and a conservative
allocator allowance of 64 B for each of the 25 owned vectors.
`IntegratorState`, the drift epoch and the pending output are fixed-size
members of the store.

Allocation, admission and release:

* **Steady state.** When the retained capacities cover the current extents
  and a reservation exists, a capture allocates nothing. It writes into the
  existing capacity, and the schedulers export into it through
  `exportPersistentStateInto`.
* **Growth** (the first capture, or larger extents):
  1. Release the old storage and its reservation first, so old and new
     allocations never overlap.
  2. Reserve the exact fresh-storage bytes.
  3. Only then allocate (`reserve` per lane).
  4. If the allocator retains more than requested, grow the reservation; the
     new grant is held before the old one is dropped (transient over-count,
     never under-count).
* **Refusal.** A governor refusal or `std::bad_alloc` on a rank releases that
  rank's storage. The collective vote then picks the same branch on every
  rank. A rank that admitted while a peer refused keeps its governed
  allocation for the next opening.
* **Release.** Everything is released at segment end, and also on failure
  through an RAII guard.
* **Reported.**
  * `multirate_rollback_snapshot_peak_bytes`: the reservation peak.
  * `…_peak_capacity_bytes`: the retained-capacity peak. It is always ≤ the
    reservation peak.
  * `…_allocations`: the number of fresh allocations, which is 1 when
    capacity is reused across cycles and retries.
  * `…_retained_bytes`: 0 after the segment.
* **Limit of the model.** It accounts for buffer capacity, not process RSS.
  Allocator headers beyond the 64 B allowance and heap fragmentation are not
  bounded exactly.

The force endpoint is a row-aligned `TransactionalForceEndpoint`. It holds no
ID lane: rows follow the state row order at capture.

Everything else a cycle must not change is covered by an FNV-1a fingerprint
(word-wise) taken at the opening, which also serves as the endpoint's
row-identity certificate:

* identity, mass, species, flags and rank;
* softening and its mask;
* particle index generation and metadata representation;
* cell and gas container sizes.

If the fingerprint differs at rollback, the rollback fails closed.

**Rollback sequence.** It runs in two collective phases, so nothing is
restored on any rank unless every rank can restore.

1. **Validation (read-only).** Each rank:
   * recomputes the fingerprint and checks the extents;
   * calls `GravityRuntime::validateTransactionalForceEndpoint`. This checks
     the particle count, the particle-index generation, the certificate, every
     endpoint and runtime lane extent and the per-row metadata extents, with
     no per-row scan.

   The verdict is agreed collectively; a refusal on any rank fails every rank
   with nothing restored.
2. **Restore.**
   1. Copy the authoritative lanes back. The extents are equal, so this
      allocates nothing.
   2. Advance the gravity source generation strictly beyond every generation
      the attempt used.
   3. Call `GravityRuntime::rollbackCycleTransaction`. It:
      * restores the endpoint **row-aligned**: no ID hash map, no allocation,
        all preconditions checked before the first write. The ordinary
        restart importer still remaps by stable ID, because a checkpoint may
        come from a different row order, and both use
        `internal::remapForceEndpointById` / `restoreForceEndpointRowAligned`;
      * re-stamps the force cache and per-row hierarchical force generations
        at the new source generation;
      * moves the journal to the rejected-attempt totals;
      * forces a full tree build on the next solve;
      * invalidates the particle ghost layer and the borrowed-source
        certificates.
   4. Set `integrator.pm_source_generation`, resynchronize the time-bin
      mirrors, reinstall the domain geometry and reconcile memory.

**Retained derived caches, and why they stay valid.**

* LET domain leaves and hierarchy: keyed by a content fingerprint of the
  geometry.
* The sparse peer graph: communication topology only, and payloads are
  exchanged every solve. Ownership is frozen inside a cycle.
* The uniform-mass certificate: the mass lanes are frozen and fingerprinted.
* The PM mesh: no PM refresh happens inside a cycle, so the mesh still
  describes the restored opening sources.

**Determinism policy.** A rollback leaves derived state that depends only on
the restored authoritative state, never on how far the rejected attempt got.
A replay of the same rejection pattern is therefore deterministic. With the
reference tree settings (no reuse or refit) the retried trajectory equals the
uninterrupted one bitwise. With reuse or refit enabled it agrees to roundoff,
because the forced full build sums in a different order than a refit would.

**Memory policy** (`numerics.treepm_multirate_rollback_memory_policy`):

* `single_block_fallback` (default): every rank opens a one-block cycle. Its
  single block honours the particle criteria and, because
  `coarse_dt < T ≤ Δt_PMdisp`, the PM displacement bound. A one-block cycle
  has no interior boundary, so it needs no rollback.
* `fail`: the run fails collectively, so the schedule never depends on the
  memory budget (reproducibility mode).

Fallbacks are recorded as cycle limiter `rollback_memory` in
`HierarchicalBlockRecord::pm_cycle_length_limiter`, distinct from the physical
limiters (`target_blocks`, `pm_displacement`, `output_or_end_boundary`,
`rollback_retry`). They are also counted in `multirate_memory_limited_cycles`
and `time.hierarchical_multirate_cycle_memory_limited`. The reference `exact`
policy never takes a snapshot.

**Diagnostics.**

* Counters:
  * `time.hierarchical_coarse_limited_by.<criterion>`
  * `time.hierarchical_multirate_cycle_limited_by.<limiter>`
  * `…_cycle_shortened_by_boundary`
  * `…_cycle_rejected_retried`
  * `…_cycle_memory_limited`
  * `…_memory_policy_failures`
  * `…_cycles_completed`
  * `…_recovery_failures`
  * `…_soft_stop_extension_blocks`
  * `…_cycle_end_inexact`
  * `time.hierarchical_quantum_unrepresentable`
* **Accepted-trajectory statistics** describe only the committed trajectory
  and equal an uninterrupted run's:
  * `treepm_long_range_refresh_count`, `…_reuse_count`;
  * the TreePM cadence records (bounded window of 256, journaled per cycle,
    so a rejected attempt can never evict accepted entries);
  * `hierarchical_block_records`;
  * `completed_steps`, `multirate_soft_stop_extension_blocks`.
* **Execution-cost statistics** include rejected work:
  * `multirate_rejected_blocks`, `multirate_rejected_pm_refreshes`,
    `multirate_rejected_pm_reuses`, `multirate_rejected_cadence_records`;
  * all profiler timings, stage counters and the soft-stop counter.
* Each `HierarchicalBlockRecord` carries `pm_displacement_limit_time_code`.
  For every accepted block, `T ≤ Δt_PMdisp` within `1e-12` relative.

### 6. Interior PM values, source prediction and authoritative epochs

* Interior PM values are needed only as the accuracy scale (total history =
  current Tree + cycle-opening PM) for rung criteria and the relative MAC.
  They are never kicked.
* Source prediction is not used. Every source is canonically drifted at
  every fine tick (section above). The global-KDK predictor fails closed for
  hierarchical lanes.
* The authoritative epochs are:
  * the PM field version and source generation of the cycle opening, which
    the opening kick requires;
  * the closing refresh at the cycle end;
  * the per-row Tree generation, equal to the current source generation for
    every tree kick.

  These checks reject a double or missing PM kick.

### 7. Approximation order (what can and cannot be claimed)

* With a fixed outer `T`, fixed inner block intervals and fixed rungs, the
  composition is a symmetric r-RESPA splitting of time-symmetric inner KDK
  blocks. For autonomous Hamiltonians this is time-reversible and second
  order in `T` and in the inner steps.
* The cosmological system is non-autonomous through `a(t)`. Kick coefficients
  are exact integrals over each half interval with `A` frozen, so the
  splitting remains formally second order. Reversibility holds only in the
  extended phase space with `t` as a coordinate.
* **Adaptive inner blocks, adaptive rungs, a retried `T` and the sliver split
  all break strict time-reversibility.** The adaptive reference scheme has
  the same limitation. No reversibility or symplecticity claim is made under
  adaptive stepping.
* Accuracy is expected to remain second order in the maximum step size for
  smooth criteria. **This has not been measured.** The owner gates are:
  * `test_hierarchical_multirate_convergence` (CPU-only). It drives the
    production block driver with analytic split forces and gates the observed
    order in `[1.8, 2.2]` for:
    * PM-only, which must also equal leapfrog at step `T`;
    * tree-only, which must also equal the reference scheme;
    * combined fast/slow oscillators with uniform and variable block
      intervals, including the splitting error relative to the reference;
    * adaptive and changing rungs;
    * cosmological (Einstein-de Sitter) kick coefficients: the canonical
      momentum must be exact under constant split accelerations for any block
      tiling, rung history and PM cycle, which detects a duplicated or missing
      PM half-kick; positions must converge at second order to the analytic
      integral;
    * a cosmological split oscillator with changing rungs, against a fine
      reference;
    * a circular two-body orbit under a smooth Tree/PM split. Linear and
      angular momentum must also be conserved to roundoff.
  * `test_hierarchical_treepm_workflow` (HDF5/FFTW, np 1/2/4). It gates
    monotone convergence of the full TreePM workflow.
* The known RESPA resonance limitation applies for large `T` relative to
  short-range dynamical times. The PM displacement criterion is the
  safeguard. It bounds displacement relative to `r_s`, not resonance.

### 8. Stop and restart policy (cycle-boundary-only, enforced)

* Interior block boundaries are synchronized but carry
  `multirate_pm_cycle_interior`. They are neither output- nor restart-safe.
  Pending snapshots and checkpoints stay latched until the cycle end, and
  analysis is skipped at interior boundaries.
* **Soft stop.** A segment whose requested step budget is reached inside an
  open cycle keeps running to the cycle end
  (`multirate_soft_stop_extension_blocks`, `report.completed_steps` and
  `report.requested_steps`). A segment ending with an open cycle throws a
  `logic_error`. There is no wall-clock or operator hard stop in the
  hierarchical segment. A future one must defer to the cycle end in the same
  way.
* Every checkpoint is therefore at a closed cycle. A resumed segment opens
  its next cycle from restart-persisted state only:
  * the scheduler;
  * the integrator, including the previous quantum for the growth limiter;
  * the v24 force endpoint;
  * the pending output.

  This is exactly what the uninterrupted run does at the same boundary.
* **Interior-cycle restart is not supported**, by design. No schema change
  was made, and none is needed. A restart that is not restart-safe is
  refused at segment start.

**Segmentation invariance claim.** For a fixed configuration, the
cycle-opening times, `T`, the block intervals and the PM refresh sequence
depend only on the state at cycle boundaries and on physical criteria. They
are therefore identical for:

* one uninterrupted run;
* runs split at any cycle boundary;
* runs with different step budgets, compared at a common cycle-end
  synchronization point.

States are expected to agree bitwise when the tree is rebuilt every solve
(the default). With topology refit or identical-source reuse enabled, the
force summation order can differ after a restart or rollback. That gives a
roundoff-level difference, not a different schedule.
