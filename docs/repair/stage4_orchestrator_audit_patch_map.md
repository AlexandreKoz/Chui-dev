# Stage 4 orchestrator audit and patch map

Date: 2026-05-18
Scope: audit-only infrastructure repair note for the Stage 4 orchestrator refactor. No solver behavior is changed by this note.

## Current runtime execution model

- The production workflow loop is `ReferenceWorkflowRunner::runImpl`: it builds scheduler active elements, splits them into particle and cell vectors, latches pending snapshot/restart intent for the completed step, calls `StepOrchestrator::executeSchedulerSubstep`, updates adaptive bins, closes the scheduler substep, syncs time-bin mirrors, then writes snapshot/restart artifacts only after the step has completed.
- `StepOrchestrator::executeSingleStep` is the only production broadcast dispatch site. It validates active-set freshness, classifies the step boundary, prepares KDK/cosmological factors, schedules the canonical stages, and broadcasts every stage to every registered callback. Each callback self-filters by `context.stage`.
- The canonical stage order is currently `kGravityKickPre`, `kDrift`, `kForceRefresh`, `kSourceTerms` (Strang pre), `kHydroUpdate`, `kSourceTerms` (Strang post), `kGravityKickPost`, `kAnalysisHooks`, `kOutputCheck`.
- Active-set authority is scheduler/state-owned. `HierarchicalTimeBinScheduler::beginSubstep` creates active scheduler elements, the workflow maps those elements to particle/cell vectors, and `makeSchedulerActiveSetDescriptor` stamps state generations and scheduler tick provenance before callbacks run.

## Broadcast/self-filter dispatch sites and filters

Production broadcast site:

- `StepOrchestrator::executeSingleStep`: loops `ordered_stages`, sets `context.stage`, injects PM directives on `kGravityKickPre` and `kForceRefresh`, and calls every registered `IntegrationCallback::onStage` for every stage.

Production callback self-filters:

- `StageAuditCallback`: records every stage and does not filter.
- `DriftCallback`: returns unless `context.stage == kDrift`.
- `GravityStageCallback`: handles `kGravityKickPre`, `kForceRefresh` through the PM directive, and `kGravityKickPost`; returns on all other stages.
- `HydroStageCallback`: returns unless `context.stage == kHydroUpdate` and cells exist.
- `StarFormationCallback`: returns unless `context.stage == kSourceTerms`.
- `BlackHoleAgnCallback`: returns unless `context.stage == kSourceTerms`.
- `TracerCallback`: returns unless tracers are enabled and `context.stage == kSourceTerms`.
- `DiagnosticsCallback`: returns unless `context.stage == kAnalysisHooks` and diagnostics are enabled.

Non-callback boundary work:

- Output/restart is not a callback today. `maybeWriteOutputs` runs after the orchestrator returns and after scheduler state is closed/synced, gated by pending output flags and boundary safety checks.

## Implicit mutation sites reachable from orchestration

- Orchestrator-owned state: `inside_kdk_step`, current/last boundary fields, integrator timeline commit, PM sync commit, `pm_long_range_field_valid`, particle drift epoch sidecars, and new-particle drift epoch repair after source-term mutation.
- Drift callback: active owned particle positions and all gas-cell centers are updated from gas particle positions.
- Gravity callback: TreePM caches and diagnostics mutate; active particle velocities are kicked; particle/cell acceleration caches are updated; inactive source prediction can read drift-epoch sidecars during local PM refreshes.
- Hydro callback: gas conserved/primitive scratch mutates, then gas density/pressure/internal energy, cell mass, gas-particle mass, and gas-particle velocities are written back.
- Source-term callbacks: star formation and black-hole AGN can mutate particle/cell state and particle counts; tracer support mutates tracer sidecars.
- Analysis callback: diagnostics bundles and retention policy mutate analysis output state, not solver arrays.
- Workflow loop: adaptive bin assignments, scheduler coverage, scheduler substep closure, time-bin mirrors, state metadata, profiler events, snapshots, and restarts mutate after orchestration returns.

## Stages needing explicit Stage 4 contracts

- `kGravityKickPre`: first KDK kick; may perform the global initial TreePM force bootstrap if the long-range field cache is invalid and the boundary is not a local substep.
- `kDrift`: active-particle position drift; gas-cell center synchronization follows particle drift.
- `kForceRefresh`: legal PM refresh surface; must carry integrator-owned PM directive and, for local substeps, predicted inactive-source requirements.
- `kHydroUpdate`: hydro-only state update and gravity-source coupling via gravity callback acceleration caches.
- `kSourceTerms`: source-term mutations, including particle-count growth; orchestrator repairs drift epoch sidecars for newly created particles.
- `kGravityKickPost`: second KDK kick from cached/fresh force data; no current PM cadence directive is issued here despite `GravityStageCallback::PmSyncSurface::kKickPost` existing.
- `kAnalysisHooks`: diagnostics-only surface.
- `kOutputCheck`: marker stage only; durable output/restart writes occur after the orchestrated step.

## Hidden PM sync/refresh paths

- Initial PM bootstrap is hidden inside `kGravityKickPre` when PM refresh is enabled, the field is invalid, and the boundary is global enough for bootstrap.
- Regular PM cadence is issued at `kForceRefresh`; the orchestrator registers the kick opportunity, the gravity callback consumes the directive, rank-consensus checks the decision, and the orchestrator commits refreshed field state after callback dispatch.
- Local-bin PM refresh requests require predicted inactive source positions and are marked unsafe for output/restart boundaries.
- TreePM solver execution, cadence records, long-range refresh/reuse counters, and `m_has_long_range_field` remain callback-owned. Cadence legality and persistent PM synchronization state are integrator-owned.

## Output/restart safety boundaries

- The workflow latches output intent before the step, maps it to `kSnapshotPoint` or `kCheckpointPoint`, passes that requested boundary into the orchestrator, and attempts writes after scheduler close/sync.
- `classifyStepBoundary` marks local active-bin steps as not restart-safe and not output-safe; `assertCanWriteSnapshotAtBoundary` and `assertCanWriteCheckpointAtBoundary` reject writes while inside KDK or on unsafe boundaries.
- HDF5 output/restart remains build-gated. When enabled, snapshot and restart writes are round-tripped immediately and restart provenance includes distributed gravity/PM cadence fields and slab ownership.

## Current test coverage and gaps

Covered today:

- Stage order and basic broadcast dispatch are covered by time-integration unit/integration tests.
- Active-set freshness, scheduler provenance, local/global boundary classification, and PM synchronization persistence/cadence basics are covered in `tests/unit/test_time_integration.cpp`.
- Reference workflow tests cover end-to-end canonical stage order, TreePM cadence records, distributed TreePM rank consensus, restart compatibility, and output/restart round trips.

Gaps for Stage 4:

- No test asserts that each production callback advertises/accepts exactly the stage set it handles; filters are duplicated as private `if (context.stage ...) return` checks.
- No test can fail when a newly registered callback silently ignores or mutates an unintended stage, because the orchestrator still broadcasts every stage to every callback.
- PM bootstrap on `kGravityKickPre`, regular `kForceRefresh`, and no-op `kGravityKickPost` PM behavior are not expressed as one explicit contract table.
- Output/restart remains outside the callback stage model, so `kOutputCheck` has no direct writer contract.
- Hydro active-face selection is callback-local and not tied to an explicit stage active-set contract.

## Minimal patch map for prompts 4.1-4.4

1. **Prompt 4.1: stage contract table and callback declarations.** Patch `include/cosmosim/core/time_integration.hpp` and `src/core/time_integration.cpp` to add an internal `StageContract`/`CallbackStageMask` helper and to let callbacks declare handled stages without changing solver math. Patch tests in `tests/unit/test_time_integration.cpp` to assert the canonical table.
2. **Prompt 4.2: replace broadcast/self-filter with orchestrator-side selective dispatch.** Patch `StepOrchestrator::registerCallback`/dispatch internals and the production callbacks in `src/workflows/reference_workflow.cpp`, `src/analysis/diagnostics.cpp`, `src/physics/star_formation.cpp`, `src/physics/black_hole_agn.cpp`, and `src/physics/tracer_support.cpp` to move stage filtering to the central dispatcher while preserving callback behavior.
3. **Prompt 4.3: make PM refresh surfaces explicit.** Patch `PmRefreshDirective`, `StepOrchestrator::executeSingleStep`, and `GravityStageCallback` to encode the initial-bootstrap, force-refresh, and post-kick no-refresh contracts in one auditable table. Extend PM cadence tests in `tests/unit/test_time_integration.cpp` plus reference workflow cadence tests.
4. **Prompt 4.4: isolate output/restart boundary contract.** Patch `ReferenceWorkflowRunner::runImpl`, `maybeWriteOutputs`, and boundary helpers in `core/time_integration` so `kOutputCheck` carries an explicit writer boundary decision while writes remain after scheduler close/sync. Extend restart/snapshot boundary tests and reference workflow output tests.

Implementation guardrails for all four prompts: do not change solver equations, KDK factors, TreePM numerics, hydro numerics, snapshot/restart schemas, or config keys; each prompt should add targeted tests for the newly centralized contract it touches.
