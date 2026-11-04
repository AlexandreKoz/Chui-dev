# Campaign A runtime-decomposition execution log

Date started: 2026-07-16  
Mode: repair  
Status: M0-M8 complete

This is the living, command-backed record for Campaign A in
`runtime_decomposition_goal_spec.md`. A milestone is complete only when its
production path and focused validation pass. The broader runtime-brain log is
historical context and is not used to claim completion of this campaign.

## Scope boundary

Campaign A decomposes the existing rung-zero production runtime. It does not
implement hierarchical KDK, distributed IC ingestion, hydro wake-up, elastic
restart, asynchronous output, new physics, or unrelated cleanup. Public config
semantics and restart schema v21 remain unchanged.

## Baseline worktree and environment

- Branch: `main...origin/main`.
- The worktree began dirty and all existing changes are user-owned. At the M0
  checkpoint, 121 tracked files differed from `HEAD` with 17,338 insertions and
  2,139 deletions, in addition to untracked hardening sources, tests, and docs.
- The current hardening state is the required implementation base. No reset,
  checkout, stash, or destructive cleanup is permitted.
- `src/workflows/reference_workflow.cpp`: 8,008 physical lines.
- Toolchain and configured build trees already provide CPU-debug, HDF5-debug,
  and MPI/HDF5/FFTW-debug coverage.
- File-producing tests in this desktop environment use
  `TMPDIR=/tmp TEMP=/tmp TMP=/tmp`.
- MPI launch requires the approved unrestricted CTest invocation because PMIx
  cannot create its listener socket inside the filesystem sandbox.

## Baseline ownership and callback map

| Current production responsibility | Current source seam | Target owner |
| --- | --- | --- |
| adaptive criteria and scheduler synchronization | lines 348-897, 3580-3830 | `TimeCoordinator` |
| gas migration/compaction and decomposition records | lines 898-1449 | `MigrationBalanceRuntime` |
| particle/AMR wire exchange | lines 1450-1849 | `MigrationBalanceRuntime` |
| restart topology validation | lines 1851-1963 | `OutputRestartRuntime` |
| runtime rebalance and AMR ownership commit | lines 1964-2568 | `MigrationBalanceRuntime` |
| initial gravity-aware decomposition | lines 2569-2808 | IC-to-migration handoff |
| gravity and hydro ghost exchange | lines 2809-3269 | gravity/hydro owners |
| run paths, identity, digest and verification | lines 3272-3580 | runtime paths/output verification |
| stage audit and drift callbacks | lines 3831-3922 | analysis/time owners |
| gravity callback | lines 3923-5279 | `GravityRuntime` |
| hydro/AMR callback | lines 5280-6558 | `HydroAmrRuntime` |
| snapshot/restart verification and writes | lines 6559-7008 | `OutputRestartRuntime` |
| output cadence callback | lines 7029-7217 | `OutputRestartRuntime` |
| startup, construction, loop and finalization | lines 7218-8008 | composition root plus `TimeCoordinator` |

The canonical stage order is
`gravity_kick_pre`, `drift`, `force_refresh`, `source_terms` (Strang pre),
`hydro_update`, `source_terms` (Strang post),
`gravity_kick_post`, `analysis_hooks`, `output_check`. The current registration
order is stage audit, drift, gravity, hydro, star formation, black hole,
diagnostics, output. Both orders are behavior-preservation invariants.

## M0 command evidence

Commands run against the unedited Campaign A baseline:

```text
git status --short --branch
git diff --stat
git diff --check
wc -l src/workflows/reference_workflow.cpp
```

Outcome: branch and dirty state recorded; `git diff --check` passed; workflow
line count was 8,008.

Focused CPU command:

```text
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/cpu-only-debug \
  --output-on-failure \
  -R '^(unit_time_integration|unit_runtime_capabilities|unit_runtime_services|unit_module_registry|integration_reference_workflow|integration_reference_workflow_hydro_row_order|integration_restart_equivalence_reference_workflow_hydro|integration_core_dependency_direction)$'
```

Outcome: PASS, 8/8.

Focused HDF5 command:

```text
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/hdf5-debug \
  --output-on-failure \
  -R '^(unit_ic_reader|unit_snapshot_hdf5_schema|unit_restart_checkpoint_schema|integration_reference_workflow|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip|integration_restart_equivalence_reference_workflow_hydro)$'
```

Outcome: PASS, 7/7.

Focused MPI command:

```text
ctest --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
  -R '^(integration_reference_workflow_distributed_treepm_mpi_two_rank|integration_reference_workflow_distributed_isolated_pm_mpi_two_rank|integration_reference_workflow_distributed_hydro_mpi_two_rank|integration_reference_workflow_distributed_amr_mpi_two_rank|integration_distributed_gas_cell_migration_mpi_two_rank|integration_distributed_sfc_rebalance_mpi_two_rank)$'
```

Outcome from the pre-implementation baseline: PASS, 6/6 outside the PMIx
sandbox; inside-sandbox failures were environment-only socket errors.

## Milestone status

- [x] M0: canonical audit/spec paths, baseline evidence, callback map, and living log.
- [x] M1: shared contracts, resource epochs/views, typed descriptor registry, and owner skeletons.
- [x] M2: output/restart extraction with unchanged schemas and cadence.
- [x] M3: IC startup and migration/balance extraction.
- [x] M4: gravity extraction and force-cache lifecycle ownership.
- [x] M5: hydro/AMR, source, and analysis extraction.
- [x] M6: time-coordinator cutover and removal of unrestricted production callback access.
- [x] M7: authoritative descriptor-driven composition and probe-module proof.
- [x] M8: full matrix, structural gates, documentation, and adversarial diff review.

## Design decisions

1. Generic epoch and view primitives remain dependency-safe; workflow-specific
   factories and descriptors live in `workflows`, never `core`.
2. `core::moduleNames()` may remain as a layer catalog, but its descriptive
   string table cannot be the authoritative runtime registry.
3. Owner services communicate through public/narrow contracts, never another
   owner's private implementation header.
4. The production cutover may use temporary adapters per milestone, but all
   unrestricted `StepContext` adapters must leave production by M6.
5. Registry ordering is explicit and deterministic. Static-initialization
   self-registration is rejected because it weakens ordering and lifetime proof.
6. Extraction and semantic change are separate checkpoints. Existing numerical
   bodies move intact before their access surfaces are narrowed.
7. No config or restart schema change is expected or authorized for Campaign A.

## Files changed by Campaign A milestone

### M0

- `docs/repair/chui_runtime_brain_prompt_result_adversarial_audit.md` (canonical path)
- `docs/repair/runtime_decomposition_goal_spec.md` (canonical path)
- `docs/repair/runtime_decomposition_exec_plan.md`

### M1

- `include/cosmosim/workflows/runtime_resources.hpp`
- `src/workflows/runtime_resources.cpp`
- `include/cosmosim/workflows/runtime_module_registry.hpp`
- `src/workflows/runtime_module_registry.cpp`
- `tests/unit/test_runtime_resources.cpp`
- `tests/unit/test_runtime_module_registry.cpp`
- `CMakeLists.txt`

Focused validation:

```text
cmake --build --preset build-cpu-debug -j4 \
  --target test_unit_runtime_resources test_unit_runtime_module_registry
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/cpu-only-debug \
  --output-on-failure -R '^(unit_runtime_resources|unit_runtime_module_registry)$'
```

Outcome: PASS, 2/2. The production epoch adapter observes existing
`SimulationState`, scheduler, and integrator authorities. Typed registry tests
cover factory/task matching, deterministic construction/execution ordering,
prerequisites, incompatibilities, freeze behavior, and descriptor removal.

### M2

- `src/workflows/internal/output_restart_runtime.hpp`
- `src/workflows/output_restart_runtime.cpp`
- `src/workflows/internal/output_verification.hpp`
- `src/workflows/snapshot_verification.cpp`
- `src/workflows/restart_verification.cpp`
- `src/workflows/reference_workflow.cpp`
- `CMakeLists.txt`

Output cadence initialization, ordered-event state, snapshot/restart assembly,
field-aware verification, and the temporary output stage adapter now belong to
`OutputRestartRuntime`. It depends on `GravityRestartStateProvider`, not the
concrete gravity callback. Snapshot/restart schemas and callback ordering are
unchanged. `reference_workflow.cpp` is 7,065 lines after M2; the largest new M2
unit is `output_restart_runtime.cpp` at 775 lines.

Focused validation:

```text
cmake --build --preset build-hdf5-debug -j4 \
  --target test_integration_reference_workflow \
           test_integration_restart_equivalence_reference_workflow_hydro \
           test_integration_restart_checkpoint_roundtrip \
           test_integration_snapshot_hdf5_roundtrip \
           test_unit_restart_checkpoint_schema
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/hdf5-debug \
  --output-on-failure \
  -R '^(integration_reference_workflow|integration_restart_equivalence_reference_workflow_hydro|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip|unit_restart_checkpoint_schema)$'
```

Outcome: PASS, 5/5. A follow-up run after moving cadence initialization passed
the two reference/restart-equivalence workflows, 2/2.

### M3

- `src/workflows/internal/initial_condition_runtime.hpp`
- `src/workflows/initial_condition_runtime.cpp`
- `src/workflows/internal/migration_balance_runtime.hpp`
- `src/workflows/migration_balance_runtime.cpp`
- `src/workflows/amr_migration_payload.cpp`
- `src/workflows/reference_workflow.cpp`
- `CMakeLists.txt`

`InitialConditionRuntime` now owns the existing single-file/generated startup
path, manifest creation, state materialization, metadata finalization, and the
initial decomposition handoff. `MigrationBalanceRuntime` now owns initial
placement, measured rebalance, migration transport and commit, compaction,
scheduler remap, identity reduction, ownership-generation invalidation, and AMR
payload exchange/validation. The reusable AMR payload projection is a separate
173-line implementation, avoiding a replacement migration monolith. No
distributed IC ingestion or restart topology relaxation was introduced.

At this checkpoint, `reference_workflow.cpp` is 5,036 lines,
`migration_balance_runtime.cpp` is 1,833 lines, and
`amr_migration_payload.cpp` is 173 lines.

Focused CPU/HDF5 validation:

```text
cmake --build --preset build-cpu-debug -j2
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/cpu-only-debug \
  --output-on-failure \
  -R '^(unit_runtime_services|unit_runtime_resources|integration_reference_workflow|integration_reference_workflow_hydro_row_order|integration_restart_equivalence_reference_workflow_hydro)$'

cmake --build --preset build-hdf5-debug -j2
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/hdf5-debug \
  --output-on-failure \
  -R '^(unit_ic_reader|unit_snapshot_hdf5_schema|unit_restart_checkpoint_schema|integration_reference_workflow|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip|integration_restart_equivalence_reference_workflow_hydro)$'
```

Outcome: PASS, CPU 5/5 and HDF5 7/7.

Focused MPI validation (outside the filesystem sandbox for PMIx sockets):

```text
cmake --build --preset build-mpi-hdf5-fftw-debug -j2
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest \
  --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
  -R '^(integration_species_migration_invariants|integration_gas_cell_identity_migration|integration_amr_patch_migration|integration_distributed_sfc_rebalance_mpi_(two|three|four)_rank|integration_distributed_gas_cell_migration_mpi_two_rank|integration_reference_workflow_distributed_(hydro|amr|treepm)_mpi_two_rank)$'
```

Outcome: PASS, 10/10. This covers exact sidecar/identity migration,
multi-rank SFC rebalance, AMR ownership transfer, gas-cell migration, and the
distributed workflow paths that consume the extracted owner.

### M4

- `include/cosmosim/workflows/gravity_runtime.hpp`
- `src/workflows/gravity_runtime.cpp`
- `src/workflows/internal/particle_ghost_runtime.hpp`
- `src/workflows/particle_ghost_runtime.cpp`
- `src/workflows/internal/output_restart_runtime.hpp`
- `src/workflows/output_restart_runtime.cpp`
- `src/workflows/reference_workflow.cpp`
- `docs/architecture/overview.md`
- `docs/reference_workflow.md`
- `CMakeLists.txt`

The former gravity callback is now a factory-created `GravityRuntime` owner.
It owns TreePM/PM cadence, force-cache compatibility/import/export,
decomposition-epoch invalidation, gravity health diagnostics, zoom source
selection, and the gravity stage contribution. Hydro receives only the
read-only `GravityAccelerationProvider`; output/restart receives only
`GravityRestartStateProvider`. The shared blocking ghost refresh projection is
a 195-line workflow helper rather than duplicated owner logic. The temporary
`IntegrationCallback` base remains explicitly limited to the M6 time-dispatch
cutover.

At this checkpoint, `reference_workflow.cpp` is 3,507 lines and
`gravity_runtime.cpp` is 1,699 lines. No replacement source exceeds the
2,000-line guardrail.

Focused validation:

```text
cmake --build --preset build-cpu-debug -j2
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/cpu-only-debug \
  --output-on-failure \
  -R '^(unit_pm_solver|unit_tree_gravity|unit_tree_pm_split_kernel|integration_pm_periodic_mode|integration_pm_scale_factor_regression|integration_tree_pm_coupling_periodic|integration_tree_gravity_vs_direct|integration_softening_ownership_invariants|integration_reference_workflow|integration_restart_equivalence_reference_workflow_hydro|integration_restart_equivalence_treepm)$'

cmake --build --preset build-hdf5-debug -j2
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/hdf5-debug \
  --output-on-failure \
  -R '^(integration_reference_workflow|integration_restart_equivalence_reference_workflow_hydro|integration_restart_equivalence_treepm|integration_restart_checkpoint_roundtrip|unit_restart_checkpoint_schema)$'
```

Outcome: PASS, CPU 11/11 and HDF5 5/5.

Focused MPI validation (outside the PMIx sandbox):

```text
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest \
  --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
  -R '^(integration_reference_workflow_distributed_(treepm|isolated_pm|hydro)_mpi_two_rank|validation_phase2_mpi_gravity_two_rank)$'
```

Outcome: PASS, 4/4. The source-transfer audit found no truncation markers,
`git diff --check` passed, and restart-equivalence preserved the existing
stable-ID/tolerance policy.

### M5

- `include/cosmosim/workflows/hydro_amr_runtime.hpp`
- `src/workflows/hydro_amr_runtime.cpp`
- `include/cosmosim/workflows/source_runtime.hpp`
- `src/workflows/source_runtime.cpp`
- `include/cosmosim/workflows/analysis_runtime.hpp`
- `src/workflows/analysis_runtime.cpp`
- `src/workflows/reference_workflow.cpp`
- `docs/architecture/overview.md`
- `docs/reference_workflow.md`
- `CMakeLists.txt`

`HydroAmrRuntime` now owns fixed-Cartesian and supported AMR orchestration,
conserved/primitive/geometry workspaces, blocking ghost freshness, remote
boundary exchange, AMR patch payload validation, reflux handoff, and accepted
CFL diagnostics. Its caller sees only aggregate communication counts and
diagnostics and supplies gravity through `GravityAccelerationProvider`.
`SourceRuntime` owns construction and the preserved star-formation then
black-hole registration order. `AnalysisRuntime` owns stage-audit and
diagnostics callbacks while preserving their positions around the dynamics and
source contributions. No wake-up, local hydro stepping, or new source model was
added.

At this checkpoint, `reference_workflow.cpp` is 1,928 lines,
`hydro_amr_runtime.cpp` is 1,727 lines, `source_runtime.cpp` is 74 lines, and
`analysis_runtime.cpp` is 86 lines. No production runtime source exceeds 2,000
lines.

Focused validation:

```text
cmake --build --preset build-cpu-debug -j2
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/cpu-only-debug \
  --output-on-failure \
  -R '^(unit_hydro_reconstruction|unit_hydro_riemann|unit_hydro_core_solver|unit_star_formation|unit_black_hole_agn|unit_analysis_diagnostics|integration_reference_workflow|integration_reference_workflow_hydro_row_order|integration_restart_equivalence_reference_workflow_hydro|integration_restart_equivalence_hydro_toy|integration_restart_equivalence_amr_hydro|integration_restart_equivalence_amr_flux_registers|integration_restart_equivalence_amr_temporal_ghosts|integration_star_formation_box|integration_black_hole_agn_toy|integration_analysis_bundle)$'
```

Outcome: PASS, 16/16.

```text
cmake --build --preset build-hdf5-debug -j2
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/hdf5-debug \
  --output-on-failure \
  -R '^(integration_reference_workflow|integration_reference_workflow_hydro_row_order|integration_restart_equivalence_reference_workflow_hydro|integration_restart_equivalence_amr_hydro|integration_restart_checkpoint_roundtrip)$'
```

Outcome: PASS, 5/5.

Focused MPI validation outside the PMIx sandbox:

```text
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest \
  --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
  -R '^(integration_reference_workflow_distributed_(hydro|amr)_mpi_two_rank|integration_distributed_gas_cell_migration_mpi_two_rank)$'
```

Outcome: PASS, 3/3. The source-transfer audit found no truncation markers and
`git diff --check` passed.

### M6

- `include/cosmosim/core/time_integration.hpp`
- `src/core/time_integration.cpp`
- `include/cosmosim/workflows/runtime_resources.hpp`
- `src/workflows/runtime_resources.cpp`
- `include/cosmosim/workflows/time_coordinator.hpp`
- `src/workflows/time_coordinator.cpp`
- `include/cosmosim/workflows/migration_balance_runtime.hpp`
- `include/cosmosim/workflows/output_restart_runtime.hpp`
- `src/workflows/internal/amr_migration_payload.hpp`
- `src/workflows/analysis_runtime.cpp`
- `src/workflows/source_runtime.cpp`
- `src/workflows/reference_workflow.cpp`
- `tests/unit/test_runtime_resources.cpp`
- `tests/integration/test_core_dependency_direction.cmake.in`
- `tests/integration/test_runtime_decomposition_structure.cmake.in`
- `CMakeLists.txt`

`RungZeroTimeState` now owns both scheduler authorities, integrator truth,
fresh/restart initialization, time-bin mirror synchronization, and pending
output cadence. `TimeCoordinator` owns the bounded segment loop, endpoint and
ordered output clipping, active-set formation, adaptive criteria, migration
and rebalance invalidation, memory snapshots, and safe output dispatch. The
core orchestrator gained dependency-safe dispatcher entry points while its
legacy callback entry points remain for lower-level compatibility tests.

Production workflow tasks no longer register `IntegrationCallback` objects or
expose `onStage(StepContext&)`. Distinct typed views carry leases over the
particle/cell/gas-identity generations, scheduler tick, and step index. The
public view surface cannot recover `SimulationState`, `StepContext`, unrelated
sidecars, or epoch mutation. The resource unit test proves stale tick/step/
generation rejection and uses a real `reorderParticles(...)` commit to prove
compaction invalidates a captured workflow lease.

The owner dependency seams were also corrected: migration and output/restart
interfaces live under `include/cosmosim/workflows/`, while AMR migration
projection is a shared helper rather than a hydro include of migration's
private header. The dependency test now scans numerical modules for upward
workflow includes and owner sources for cross-owner private headers.

Focused validation:

```text
cmake --build --preset build-cpu-debug -j4
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/cpu-only-debug \
  --output-on-failure \
  -R '^(unit_time_integration|unit_runtime_resources|unit_runtime_module_registry|integration_time_integration_loop|integration_reference_workflow|integration_reference_workflow_hydro_row_order|integration_restart_equivalence_reference_workflow_hydro|integration_core_dependency_direction|integration_runtime_decomposition_structure)$'
```

Outcome: PASS. Follow-up source/analysis equivalence coverage also passed 11/11
(`unit_star_formation`, `unit_black_hole_agn`, `unit_analysis_diagnostics`,
their focused integrations, the real workflow, and restart equivalence).

### M7

- `include/cosmosim/workflows/runtime_module_registry.hpp`
- `src/workflows/runtime_module_registry.cpp`
- `src/workflows/internal/reference_runtime_composition.hpp`
- `src/workflows/reference_runtime_composition.cpp`
- `include/cosmosim/workflows/reference_workflow.hpp`
- `src/workflows/reference_workflow.cpp`
- `include/cosmosim/workflows/time_coordinator.hpp`
- `src/workflows/time_coordinator.cpp`
- `tests/unit/test_runtime_module_registry.cpp`
- `tests/integration/test_runtime_descriptor_probe.cpp`
- `CMakeLists.txt`

The registry now constructs the real analysis, drift, gravity, hydro/AMR,
source, and output/restart owners through ordered factories. Descriptors
declare prerequisites, typed stage view, deterministic ordinal, and resource
access. `TimeCoordinator` owns the frozen `RuntimeExecutionPlan` and performs
no direct stage-owner switch calls; its stage switch only selects the typed
view used to execute matching plan tasks.

`ReferenceWorkflowOptions::register_runtime_modules` is a non-persisted
test/embedding seam invoked before registry freeze. The focused probe test adds
an analysis diagnostic task solely through this seam. It executes once in the
real workflow. The registry unit test separately proves that omitting the
probe descriptor removes its action without modifying the runner or built-in
task list.

Focused validation:

```text
cmake --build --preset build-cpu-debug -j4 \
  --target test_unit_runtime_module_registry \
           test_integration_runtime_descriptor_probe \
           test_integration_reference_workflow
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/cpu-only-debug \
  --output-on-failure \
  -R '^(unit_runtime_module_registry|integration_runtime_descriptor_probe|integration_reference_workflow)$'
```

Outcome: PASS, including the real descriptor-probe execution. At the M7
checkpoint `reference_workflow.cpp` is 657 lines; the largest runtime source is
`migration_balance_runtime.cpp` at 1,834 lines. No runtime source exceeds the
2,000-line hard guard.

### M8

Campaign A closes with every structural gate satisfied. The production workflow
root is 657 physical lines, down from 8,008. It constructs the run paths,
configuration, root MPI/profiler services, IC/runtime state, descriptor registry,
and owner graph; calls `TimeCoordinator::runRungZeroSegment`; and emits the final
summary/digest. It contains no scheduler initialization, numerical stage loop,
gravity/hydro/AMR implementation, migration serialization, IC conversion, or
snapshot/restart field comparison.

Final runtime-source sizes:

| Production unit | Lines | Owned responsibility |
| --- | ---: | --- |
| `reference_workflow.cpp` | 657 | composition, lifecycle, top-level reporting |
| `time_coordinator.cpp` | 1,235 | rung-zero scheduling, clipping, active sets, legal output boundaries |
| `migration_balance_runtime.cpp` | 1,834 | decomposition, migration, rebalance, remap, invalidation |
| `gravity_runtime.cpp` | 1,672 | force lifecycle, TreePM/PM cadence, gravity diagnostics |
| `hydro_amr_runtime.cpp` | 1,711 | current hydro/AMR orchestration and ghost/patch freshness |
| `output_restart_runtime.cpp` | 885 | cadence, writes, verification, restart/provenance assembly |
| `reference_runtime_composition.cpp` | 382 | built-in typed descriptors and owner factories |
| `runtime_module_registry.cpp` | 369 | validation, ordering, freeze, typed execution plan |
| `initial_condition_runtime.cpp` | 169 | generated/single-file IC startup and handoff |
| `source_runtime.cpp` | 114 | existing source-model ordering and invocation |
| `analysis_runtime.cpp` | 82 | audit/diagnostic contributions |

The former workflow responsibilities now have these live authorities:

| Former root responsibility | Final owner/interface |
| --- | --- |
| initial state materialization and initial placement | `InitialConditionRuntime` -> `MigrationBalanceRuntime` |
| scheduler/integrator initialization and step lifecycle | `RungZeroTimeState` and `TimeCoordinator` |
| gravity stages and restart cache | `GravityRuntime` through gravity/provider interfaces |
| hydro, AMR, ghosts, reflux handoff | `HydroAmrRuntime` |
| star-formation and black-hole sequencing | `SourceRuntime` |
| audit and scheduled diagnostics | `AnalysisRuntime` |
| migration, compaction, scheduler remap and rebalance | `MigrationBalanceRuntime` |
| snapshot/restart cadence, payloads and verification | `OutputRestartRuntime` |
| collective phase failure consensus | `FailureCoordinator` through `RuntimeServices` |

The remaining root responsibilities are intentionally high-level: immutable
configuration/run-path construction, one root `MpiContext` and profiler bundle,
startup-versus-restart selection, descriptor registration/freeze, the single call
into the coordinator, and final digest/reporting. Those are composition and
lifecycle decisions and do not belong to a numerical owner.

Production dispatch uses `DriftParticleStageView`, `GravityStageView`,
`HydroAmrStageView`, `SourceMutationStageView`, `AnalysisStageView`, and
`OutputRestartStageView`. Their public surface exposes freshness validation only;
the owning service receives the private, friend-gated bridge required for its
authorized lanes. Leases validate the relevant particle/cell/gas-identity
generations, scheduler tick, and step index. `MigrationOwnershipView` and
`TimeCriteriaStageView` provide the corresponding non-stage capabilities. Unit
tests prove stale tick, step, and generation rejection, public absence of state or
sidecar access, and invalidation after a real particle reorder/compaction commit.
Views are transient derivations and are never restart truth.

The dependency graph check passes. `core` provides dependency-safe dispatcher
entry points but has no dependency on workflows or physics; numerical modules do
not include workflows; composition depends downward on public owner interfaces;
and no owner includes another owner's private header. Shared AMR migration payload
projection is an internal peer helper rather than an ownership inversion.

The frozen `RuntimeExecutionPlan` is assembled exclusively from typed descriptors.
Each built-in descriptor declares identity, prerequisites, construction order,
stage/view kind, and resource access, while its factory constructs the real owner.
`TimeCoordinator` dispatches the plan and has no central list of concrete owner
callbacks. `test_runtime_descriptor_probe.cpp` adds a diagnostic action solely
through `ReferenceWorkflowOptions::register_runtime_modules`, observes exactly one
execution in the real workflow, and requires no runner edit. The registry unit test
proves descriptor omission removes the contribution.

Final command evidence:

```text
./scripts/ci/check_repo_hygiene.sh
```

Outcome: PASS.

```text
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure
```

Outcome: all 129 CPU-debug tests passed. The desktop command wrapper reached its
60-second yield after tests 1-127; the two remaining tests were run directly and
passed (`validation_tree_pm_ewald_accuracy`, 23.72 seconds, and
`integration_dmo_empty_rank_smoke_single_rank`, 6.57 seconds). A final post-cutover
focused architecture/behavior run passed 14/14.

```text
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug -j4
env TMPDIR=/tmp TEMP=/tmp TMP=/tmp ctest --test-dir build/hdf5-debug \
  --output-on-failure \
  -R '^(unit_(ic_reader|snapshot_hdf5_schema|restart_checkpoint_schema|runtime_resources|runtime_module_registry)|integration_(reference_workflow|runtime_descriptor_probe|restart_equivalence_reference_workflow_hydro|snapshot_hdf5_roundtrip|restart_checkpoint_roundtrip|provenance_roundtrip|core_dependency_direction|runtime_decomposition_structure))$'
```

Outcome: PASS, 13/13, including IC, snapshot, checkpoint, provenance, and direct
versus resumed reference-workflow coverage. Restart schema v21 and public config
semantics are unchanged.

```text
cmake --preset asan-debug
cmake --build build/asan-debug -j4
env ASAN_OPTIONS=detect_leaks=0 TMPDIR=/tmp TEMP=/tmp TMP=/tmp \
  ctest --test-dir build/asan-debug --output-on-failure \
  -R '^(integration_(reference_workflow|runtime_descriptor_probe|restart_equivalence_reference_workflow_hydro|reorder_compaction_sidecars|species_migration_invariants|core_dependency_direction|runtime_decomposition_structure)|unit_(time_integration|runtime_resources|runtime_module_registry|star_formation|black_hole_agn|analysis_diagnostics))$'
```

Outcome: PASS, 13/13 under AddressSanitizer and UndefinedBehaviorSanitizer. Leak
detection cannot run in this desktop environment because LeakSanitizer aborts when
the process is traced; the same binaries passed with only leak detection disabled.
This is an environment limitation, not a suppressed test failure.

```text
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug -j4
ctest --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
  -R '^(integration_runtime_descriptor_probe|integration_reference_workflow_distributed_(treepm|isolated_pm|hydro|amr)_mpi_two_rank|integration_pm_slab_halo_exchange_mpi_(two|three|four)_rank|integration_distributed_sfc_rebalance_mpi_(two|three|four)_rank|integration_distributed_gas_cell_migration_mpi_two_rank|integration_core_dependency_direction|integration_runtime_decomposition_structure|validation_phase2_mpi_gravity_(two|three|four)_rank)$'
```

Outcome: PASS, 17/17 outside the filesystem sandbox, as required for PMIx local
sockets.

Final structural/adversarial checks:

```text
wc -l src/workflows/*.cpp
rg -n 'registerCallback|onStage\s*\(|IntegrationCallback' \
  include/cosmosim/workflows src/workflows
git diff --check
```

Outcome: the root is 657 lines, every runtime unit is below 2,000 lines, the legacy
unrestricted callback search produced no matches in production workflow code, and
`git diff --check` passed. The capability guard still rejects
`hierarchical_max_rung != 0`; no distributed IC, wake-up, elastic restart,
asynchronous output, or new physics path was added.

```text
ctest --test-dir build/cpu-only-debug --output-on-failure \
  -R '^unit_runtime_capabilities$'
```

Outcome: PASS, 1/1, confirming the rung-zero-only capability report and rejection
of nonzero hierarchical rungs remain truthful.

## Acceptance-gate result

| Gate | Result | Evidence |
| --- | --- | --- |
| A: workflow reduction | PASS | 8,008 -> 657 lines; largest runtime 1,834; structural CTest |
| B: dependency direction | PASS | extended dependency CTest and 17-test MPI matrix |
| C: resource restriction | PASS after follow-up | public owner-context escape removed; task-scoped descriptor grants; negative compile tests |
| D: descriptor registration | PASS for stage composition | real-workflow probe, descriptor omission, typed-view access allowlists, and task-scoped grants |
| E: behavior preservation | PASS | CPU 129/129, HDF5 13/13, MPI 17/17, sanitizer 13/13 |
| F: no scope escape | PASS | rung-zero guard and adversarial source/diff search |
| G: maintainability evidence | PASS | owner/view/dependency/line-count/test evidence in this log |

## Unresolved risks and limitations

- Follow-up closure removed `StepContext` and owner-context recovery from every
  public stage-view and owner interface. Each built-in source-private access now
  validates the current task's descriptor grant before entering the existing
  rung-zero body. Future lane-level kernel capabilities can replace this trusted
  internal bridge without changing composition.
- Runtime resource declarations are executable task grants and typed-view
  allowlists; they do not dynamically trace every memory access inside a trusted
  built-in owner. Timestep criteria, restart payloads, and migration payloads
  remain dedicated-service contracts rather than descriptor extension points.
- The worktree contains a large pre-existing hardening patch and generated/untracked
  outputs. Campaign A used only file-scoped edits and no reset, checkout, stash, or
  cleanup, so those user-owned changes remain intact.
- LeakSanitizer coverage remains blocked by the desktop tracer/ptrace environment;
  address and undefined-behavior sanitizer coverage passed.

## Current blockers

None for Campaign A. Later work must remain separate: hierarchical KDK,
distributed IC ingestion, hydro wake-up, elastic restart, asynchronous output, and
new physics were neither implemented nor claimed here.


## 2026-07-17 resource-capability follow-up

The adversarial Campaign A review demonstrated that a derived owner could use a
protected `stageContext()` helper to recover mutable global state. The follow-up
repair removes that API, makes the runtime resource grant task-scoped, validates
view-specific access modes, adds negative compile tests, and centralizes duplicate
gas-cell ownership helpers. See
`docs/repair/runtime_resource_capability_closeout.md`.
