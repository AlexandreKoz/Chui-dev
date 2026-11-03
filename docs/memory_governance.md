# Runtime memory governance

## TreePM bounded optional acceleration — 20261012

`TreePmCoordinator` owns a fixed `TreePmResidualCoefficientTable` and an exact
`TreePmLeafScratch[]` allocation of one 16-source block per planned worker.
Both reserve `MemoryClass::kPersistentCache` before allocation and retain the
committed reservation until physical destruction. `memoryReport()` charges
`sizeof(TreePmResidualCoefficientTable)` (16384 cells of five doubles plus
metadata) and `worker_capacity * sizeof(TreePmLeafScratch)` including alignment.
The latter holds eleven double lanes and a validity lane; no extra O(N) source
or O(node count) COM-offset array is introduced. Actual `sizeof` and capacity
reports, rather than a logical element count, are authoritative.

Table replacement destroys/releases the old table before reserving its
replacement. Its exact key is effective pair epsilon, Gaussian split length
and cutoff. Disabled/heterogeneous/unsupported inputs use the analytic law;
near cancellation and uncovered cells fall back per pair. Worker growth
likewise releases old storage first. Allocation denial selects original
indexed scalar traversal, independently of table admission. Scratch contents
are overwritten before every leaf block and retain no source-generation truth.
The existing Gaussian table remains separately owned/reported.

Optional storage is admitted only from remaining headroom, independently of
the mandatory phase estimate. Before the next production task admission,
`GravityRuntime::estimateMemory` may call the coordinator's additive
`trimOptionalAccelerationCaches(required_headroom_bytes)` at the serial task
boundary. It releases these two caches if mandatory headroom is insufficient;
it never changes installed tree/PM state. Calls during a solve are forbidden.
Standalone callers can use that same boundary API. No process-global table
pool bypasses governor ownership. Per-rank replication and worker count must
be included in node-wide projections.

Two six-double inline bounds records add 96 bytes per coordinator, folded into
fixed-feedback accounting. Source bounds refresh during existing preflight;
node bounds refresh after build/refit and remain valid on exact reuse. Existing
three periodic coordinate vectors still cost their actual retained capacities
(24 bytes per particle if capacity equals size); no replacement ownership is
introduced. Tree rebuild/refit, communication, FFTW and snapshot/restart peak
overlap retain their existing owner-local contracts. A 512-cubed fit in 48 GiB
is **unproven**. See the [live-set projection and qualification gates](repair/treepm_deep_hotpath_attempt1_closure_20261012.md).

Status: M1B foundation, 2026-08-29.

This document defines the runtime memory-governance behavior implemented by
Campaign M1B. It complements the pre-run gravity and composed DMO process
memory estimators; it does not replace them.

## Authority and ownership

The reference workflow composition root owns exactly one `core::MemoryGovernor`
for the process. `workflows::RuntimeServices` borrows that governor and selected
runtime components may in turn borrow it. The governor is not a singleton and
`core` does not depend on `workflows`.

The user-visible hard ceiling remains
`parallel.process_memory_budget_bytes`. A value of zero retains the existing
unlimited-budget meaning: runtime accounting and high-water telemetry continue,
but the governor does not reject a request because of a hard ceiling.

The production policy also reuses the existing process-level reserve semantics:

- `parallel.gravity_backend_unknown_reserve_bytes`;
- `parallel.process_mpi_unknown_reserve_bytes`;
- `parallel.process_fftw_unknown_reserve_bytes`;
- `parallel.process_hdf5_unknown_reserve_bytes`;
- `parallel.process_allocator_unknown_reserve_bytes`;
- `parallel.process_output_restart_overlap_bytes`;
- `parallel.process_memory_safety_margin_fraction`.

M1B adds no second configuration namespace and no duplicate memory budget.

## Accounting model

The governor separates four quantities that must not be added blindly to the
ordinary memory report:

1. **baseline owned bytes**: live CHUI-owned capacity observed through the
   existing `MemoryReport` but not represented by a governor commitment;
2. **governed committed bytes**: physical backing storage whose reservation was
   admitted and then committed;
3. **pending reserved bytes**: admitted demand for which backing storage has not
   yet been committed;
4. **opaque/planned reserves**: configured external-runtime reserve and the
   configured output/restart overlap allowance.

`MemoryEntry::governed_commitment` is the reconciliation marker. A physical byte
range represented by a governed entry remains visible in `MemoryReport`, but
`memoryReportBaselineOwnedBytes()` excludes it from the governor baseline because
that same range is already present in `committed_bytes`. Consequently:

```text
raw accounted demand =
    baseline owned
  + governed committed
  + pending reserved
  + external runtime reserve
  + planned overlap reserve
```

The configured process safety margin is then applied with the same
multiplicative/ceiling convention used by process preflight. `MemoryReport`
totals remain a description of owned container capacity; governor fields are a
policy/accounting overlay and are not an additional physical-memory total to be
summed with those report totals.

`MemorySubsystem` and `MemoryClass` are orthogonal. `MemorySubsystem` identifies
who owns a reported allocation; `MemoryClass` identifies its residency/purpose.
M1B memory classes are canonical persistent, persistent cache, phase resident,
scratch arena, communication, diagnostic, and external runtime.

## Reservation lifecycle

`MemoryReservation` is move-only RAII state with the following lifecycle:

```text
reserve -> pending -> backing allocation succeeds -> commit
        -> physical storage lifetime -> release/destruction
```

Pending destruction releases only the pending reservation. Committed
destruction releases committed accounting. Explicit `release()` is idempotent;
moved-from reservations are inert; a second `commit()` is a logic error.

Reservation precedes target physical allocation. If the backing allocator
throws, the pending RAII object unwinds and governor state returns to its prior
value. A reservation is retained for the physical lifetime of its allocation,
not merely for the logical phase that currently uses the bytes.

## Pressure and headroom

For finite budgets the internal M1B policy uses deterministic thresholds at 85%
and 95% of the configured hard ceiling:

- **Green**: below the Amber threshold;
- **Amber**: at or above 85% and below 95%;
- **Red**: at or above 95% while still within the hard ceiling;
- **Trip**: accounted demand is already above the hard ceiling, or a prospective
  reservation is rejected because it would cross the hard ceiling.

Pressure is observational in M1B. Adaptive batch shrinking, optional-task
postponement, cache eviction, and broad phase scheduling are M1C/later work.

## Governed physical path: monotonic scratch

`MonotonicScratchAllocator` can borrow a `MemoryGovernor`. Before each new
physical block it computes the exact intended block capacity, reserves
`MemoryClass::kScratchArena`, allocates backing storage, commits the reservation,
and stores that reservation inside the block.

`reset()` resets logical offsets only. It does not free blocks or release their
committed reservations. Reusing retained blocks therefore creates no second
reservation. Destruction of the block/allocator releases the physical
commitment. The allocator reports current logical use, retained capacity,
historical logical-use high-water, and retained-capacity high-water separately.

The production `TransientStepWorkspace` receives the process governor from the
reference workflow through `RuntimeServices`. On a force-refresh pre-kick, the
production gravity runtime admits the all-particle `uint32` index lane through
this arena after the collective DMO process preflight succeeds. On a cached
pre-kick, where that full process preflight is intentionally skipped, admission
occurs after the already-collective gravity phase reservation and before the
subsequent drift view is constructed. Rank-local reservation/allocation failures
are passed through `FailureCoordinator` before peers can diverge. The drift view
therefore reuses this one governed physical owner rather than first retaining an
N-sized compatibility vector. The admitted block survives workspace reset as
long as the backing block remains resident. Existing standalone and subset
callers may still use the compatibility path when governed scratch has not been
admitted.

Authoritative particle and gas-cell scheduler retained capacities are explicit
persistent `MemoryEntry` owners. Governor baseline reconciliation therefore
comes directly from the merged ownership report rather than adding hidden
scheduler byte totals outside that report.

## Runtime-resource leases remain separate

`workflows::RuntimeResourceLease` remains the stage-view freshness/capability
mechanism. It is not a memory reservation and has not been overloaded or
renamed. A task may eventually hold both a `RuntimeResourceLease` and a
`MemoryReservation` for independent reasons.

## Telemetry

The existing profiler `memory_report` object contains an additive `governor`
object when a governor snapshot is attached. It exposes hard limit, baseline,
external/planned reserves, committed, reserved, raw and safety-adjusted demand,
headroom, pressure, peak committed/reserved/accounted demand, and rejection
count. Human-readable memory reporting includes memory class and the governed
reconciliation marker.

M1B intentionally does not sample RSS/PSS or `/proc`; measured process-memory
reconciliation is an M1C acceptance task.

## M1B boundary / M1C handoff

M1B proves the central authority, accounting model, reservation lifecycle,
production ownership path, profiler integration, and one retained physical
allocation path. It does **not** yet reserve every major allocation. PM routing,
Tree/LET storage and exchange, FFT/PM grids, migration staging, output/restart
staging, halo/analysis buffers, hydro/AMR batches, and device memory remain
broad-integration work for M1C or later campaigns.

## M1C-1 runtime integration addendum (2026-08-31)

Campaign M1C-1 extends the M1B foundation without creating a second memory
authority. Major DMO/general-runtime high-water phases now use the existing
process governor:

- TreePM admits a conservative phase high-water collectively before target-scale
  cache/staging growth. Retained gravity capacity transfers atomically from the
  phase commitment into baseline ownership at phase completion.
- Runtime decomposition planning has a narrow governed lifetime and is destroyed
  before migration staging. Migration then admits transaction/staging coexistence
  collectively, and preserved scheduler state uses compact identity records
  rather than full particle-migration records.
- Snapshot/restart readback uses the existing
  `parallel.process_output_restart_overlap_bytes` allowance as predeclared
  headroom. Actual physical staging commits only the amount beyond that allowance,
  so policy reserve and physical staging are not counted twice.
- Optional science diagnostics are deferred under Red/Trip pressure. Required
  run-health diagnostics and scientific solver semantics are unchanged.

`MemoryReservation::reconcileBaselineOwnedAndRelease()` is the physical-lifetime
handoff primitive for phases whose admitted memory becomes retained capacity. It
updates baseline ownership and releases the committed reservation under one lock,
avoiding either a double-account window or an unaccounted retained-capacity gap.

Production process-memory observation now supplies optional Linux/WSL RSS,
peak-RSS/HWM, and `smaps_rollup` PSS values at workflow reporting boundaries.
The profiler distinguishes conservative governor policy demand from current
declared residency. RSS discrepancy uses baseline-owned retained capacity +
currently committed governed capacity + the configured opaque external-runtime
estimate; inactive future output/restart overlap and uncommitted reservations do
not masquerade as current residency. RSS/PSS remain diagnostic evidence; they do
not become scientific state or replace deterministic allocation admission.

Distributed telemetry reports rank-local/sum/max/mean/imbalance for declared
memory, RSS/peak RSS where all ranks expose it, and communication high-water.
Rank maximum remains the safety quantity.

Post-M1C-1 closure also replaces population-scaled migration wire
materialization with an explicit conditional field codec and reusable fragmented
packets bounded by `parallel::mpiTransportRoundLimitBytes()`. Migration governor
admission is derived from the actual transaction coexistence, native record and
dynamic sidecar capacities, scheduler/index staging, bounded packet buffers, and
maximum single-record assembly requirements.

The detailed coverage inventory and M1C-2 handoff are in
`docs/repair/m1c1_runtime_memory_integration_20260831.md`.

## M1C-2 acceptance closure addendum (2026-09-03)

M1C-2 closes the M1 source/architecture acceptance matrix in
`docs/repair/m1c2_m1_acceptance_closure_20260903.md`. The certified 8-rank PM
routing model now has an explicit regression for the 512^3 target population and
for the exact worst-case interpolation accounting expression:

```text
448 B rank metadata
+ 67,075,680 B retained send capacity
+ 67,075,680 B retained receive capacity
= 134,151,808 B
```

This is 65,920 B below the 128 MiB contract. `N_local` is absent from the
capacity model by design; source population changes round count/traffic, not the
simultaneous routing workspace. Runtime still checks actual retained wire
capacities plus rank metadata against the 128 MiB ceiling.

M1C-2 does not promote unavailable MPI/FFTW execution into runtime evidence and
does not redefine M2 hydro/AMR/full-physics memory scope.

## M2B-1 production AMR regrid transaction addendum (2026-09-05)

Production `SimulationState` refinement and derefinement are explicit governed
transactions. `ProductionAmrHydroOptions::regrid_memory_governor` must point to
the process `MemoryGovernor`; regrid fails closed when that authority is absent.
Before allocating candidate AMR authority or transaction scratch, the regrid
path computes a checked upper bound and reserves it as phase-resident memory.
Only AMR-owned cell, gas-sidecar, patch, and gas-identity state are staged;
unrelated particle/species/module populations are not copied. A reservation or
preparation failure leaves live AMR authority unchanged.

Canonical memory reporting now includes allocated capacity for the gas-cell
stable-identity records and lookup tables, pending flux-register records, and
AMR temporal-boundary histories including nested cell records. These entries
therefore participate in `memoryReportBaselineOwnedBytes()` and the existing
governor baseline reconciliation instead of remaining hidden heap state.

Regrid is permitted only at a synchronization-safe point. Active temporal
boundary history or any unresolved pending flux register causes deterministic
pre-mutation rejection. The pending-register representation does not encode all
fine-contributor topology identities, so this campaign deliberately uses the
strong global-drain contract rather than attempting an unsafe intersection-local
remap.

## M2D memory-aware runtime scheduling addendum (2026-09-06)

M2D extends the existing process governor from allocation admission into
explicit deterministic high-water scheduling without creating a second memory
authority. The governor now supplies checked headroom-derived batch sizing used
by hydro, bounded remote-gas communication, and stellar evolution/feedback.
Each selection is deterministic for a given snapshot and policy; unlimited
budget preserves the configured/requested maximum.

The frozen `RuntimeModuleRegistry` now also records major-task dependency,
release-boundary, optionality, declared peak, and coarse compute/DRAM/MPI
pressure metadata. The accepted serial stage dispatcher remains the numerical
execution authority. A conservative overlap predicate rejects resource-write
conflicts, simultaneous high pressure in the same limiting resource class, and
known memory peaks that exceed current governor headroom. This establishes a
safe concurrency contract without silently changing scientific order.

SFC decomposition now accepts a hard per-rank peak-memory envelope in addition
to its soft memory-pressure cost. The envelope is derived collectively from the
most constrained process-governor headroom and includes a conservative
major-task transient reserve. Memory-unsafe candidate partitions are rejected
before migration, and a current hard violation forces rebalance independently
of ordinary soft load thresholds. Rank maximum is therefore a safety constraint,
not an average optimization target.

`MpiContext` records node-local MPI size/rank and the workflow publishes
rank/thread topology in `runtime.topology`. M2D does not introduce MPI shared
windows or a new node-memory authority; any future cross-rank shared resource
must demonstrate measurable bytes/runtime payoff plus NUMA and lifetime safety.

Optional science diagnostics remain subordinate to the hard memory ceiling, but
a transient Red/Trip episode no longer drops a due diagnostic permanently. The
analysis runtime retains only two owner-local pending bits (light/heavy science)
and retries them on the first later AnalysisHooks stage whose governor pressure
is below Red. A failed/unsafe pressure state still defers rather than forcing the
work through the ceiling; successful catch-up clears the pending bit. This
pending state is runtime scheduling metadata, not scientific/restart authority.

Detailed implementation and validation evidence is in
`docs/repair/m2d_memory_aware_runtime_20260906.md`.

## M2D-1 retained workspace and task admission contracts

`selectDeterministicBatchSizeForWorkspaces` uses the existing governor snapshot and policy. Retained capacity is already accounted; growth charges the full replacement allocation, and temporary workspace is charged in full. The selector returns a decision, not ownership. The caller retains the actual physical reservation until its memory is released or reconciled. Aggregate arithmetic is checked in both finite and unlimited modes.

The runtime registry distinguishes complete/unknown peaks and owner-managed/dispatcher-owned reservations. Dispatcher-owned tasks hold one RAII lease through execution. Owner-managed preflight is released before the owner's own reservation; no physical byte range is counted twice. All existing heavy built-in owners remain owner-managed. Unknown peaks cannot authorize concurrency. The production numerical order is unchanged, and completion of the remaining owner-specific peak models is a separate acceptance gate.

## M2D final diagnostic admission (2026-09-07)

The production analysis owner holds a `kDiagnostic` lease through generation
and writing. `estimatePowerSpectrumMemory` models the real/complex coexistence,
FFT worker lines and reductions with checked 64-bit arithmetic. The bundle
model also includes explicit reusable ownership-validation scratch. Optional
hard-limit refusal is collective and defers the cadence; required health fails
in a controlled manner. Opaque FFTW planner/runtime memory remains in the
existing external-runtime reserve. The 8 MiB bundle metadata allowance is an
engineering provision, not a proof of all external allocation maxima.

The task registry receives a required-health-only incremental preflight.
Optional products have their own owner-held admission. Unknown whole-phase
peaks, including incomplete gravity/source/output contracts, remain
non-overlappable; production numerical stages remain serial. Do not use the
known required-health estimate as a complete analysis-stage peak or charge it
again while the owner holds its physical reservation. Public APIs and the
remaining limitations are documented in
`docs/repair/m2d_final_admission_closure_20260907.md`.

## M2D acceptance-gates owner contracts (2026-09-07)

The existing `RuntimeExecutionPlan` now accepts an additive, state-dependent
`RuntimeTaskMemoryEstimate` from each built-in major owner. The estimate carries
an incremental byte count, an explicit completeness bit, and a static-lifetime
uncertainty description. `taskMemoryEstimate("module::task_id")` exposes the
current model without acquiring memory or a stage view. A missing implementation
in a third-party owner defaults to incomplete, preserving source compatibility.

Gravity reuses its actual conservative TreePM phase estimator and subtracts
retained gravity baseline capacity once. Its existing physical phase lease is
unchanged. Source and hydro/AMR publish checked bounded-workspace estimates;
output/restart publishes the maximum of its separately released snapshot and
restart staging demands, using the same planned-overlap policy as the physical
writer. Disabled/no-op source and output paths can report a complete zero.

Incomplete models are **not whole-task certificates**. Their byte count is
informational and cannot authorize overlap or become a mandatory preflight:
adaptive kernels must retain the ability to select a smaller legal batch.
Complete owner-managed models may be preflighted, with the temporary lease
released before the owner acquires its actual allocation lease. Dispatcher-owned
tasks require a complete model and keep their RAII lease through execution.
The production stage dispatcher remains serial. Static unknown peaks and
resource/dependency conflicts continue to reject speculative overlap.

The remaining unmodeled terms include source population growth and ID
coordination, complete AMR transaction/ghost/geometry coexistence, and output
writer/readback metadata and capacity growth. Those owners retain their
allocation-time governors. The gravity model and the new task graph do not
establish a measured whole-process ceiling, and no concurrency is enabled by
this patch. The single process governor, baseline/commitment reconciliation,
external-runtime reserve, and existing configuration keys remain authoritative.

## M2D retained population-growth admission

The source runtime now uses `core::RetainedCapacityTransaction` at the
source-precommit boundary. It plans actual replacement capacities for canonical
particle SoA, species-specific star/BH lanes, species indices and scheduler
metadata. A single aggregate reservation covers the old/new backing coexistence
before the first physical birth mutation. Each lane is replaced without
changing its logical size; the existing source model then performs the original
resize and numerical update. Successful growth transfers the measured retained
capacity into the existing governor baseline exactly once. If a later lane
throws, the already-retained capacity remains accounted and logical state is
unchanged. A consumed transaction cannot be reused; a retry must rebuild the
plan from current capacities. The production plan is phase-local, requires
nonthrowing container swap, and uses standard numeric SoA vectors. Arbitrary
third-party allocator over-allocation is not certified: the plan rejects a
replacement capacity above its declared bound.

The ID-precommit interface has an additive no-op compatibility hook. The
production distributed registry connects it to the one source owner and its
scheduler. All ranks, including ranks with no births, enter the same collective
admission boundary. The existing governor remains the only byte authority.
Unmodeled ID coordination, source metadata, and other source scratch do not
become certified merely because canonical growth is now governed. The
hydro/AMR and output contracts remain conservative/partial, and no speculative
full-physics overlap is enabled.

The distributed SFC planner also coordinates local entry, sorting/sample
preparation, and memory-refinement metadata failures before their next
collective. Test-only fault points are `sfc_local_preparation`,
`sfc_cut_sample`, `sfc_memory_preflight`, `sfc_memory_prefix`, and
`sfc_memory_repair_metadata`. The MPI test suite must verify these on real
np2/np3/np4 ranks; CPU-only compilation is not distributed certification.

## M2D final coding follow-up — bounded ID and phase metadata (2026-09-08)

The distributed source-ID registry now retains a sorted uint64 shard rather than
a node-owning hash set. A fixed local initialization batch is exchanged through
bounded MPI rounds. Collision requests, decisions, and validation work use a
`GovernedScratchArena` with a null upstream resource; the arena has one physical
allocation and one existing-governor lease. The registry's retained shard has a
separate persistent-cache lease, and precommit results remain leased until the
source birth transaction finishes. The internal factory is test-only/owner-local;
the production SourceRuntime remains the sole coordination authority. The
`ParticleIdPrecommit::finishPrecommit()` hook defaults to a no-op for existing
implementations. No ID construction or collision ordering contract changes.

The snapshot writer now uses compact sorted sidecar-row lookup records instead
of three unordered maps. Its optional governor pointer is supplied by the
production output owner; the standalone API remains compatible. The lookup
arena is released after serialization. Restart writes module-sidecar payloads
directly from canonical bytes instead of retaining a second payload. AMR
synchronized hydro uses a governed sparse active-row index and allocation-free
patch-local membership queries, preserving the original solver active order.

These are physical lifetime improvements, not a complete whole-task certificate.
The ID scratch envelope is bounded by the admitted source batch, global record
count and configured transport round size; the retained shard scales with its
authoritative ID population. Full source metadata, AMR geometry/ghost/flux
coexistence, and all writer/readback metadata remain incomplete owner models.
The existing serial dispatcher, local physical leases, incomplete-contract
overlap refusal, and external-runtime reserve remain authoritative.

## P6 owners and coexistence (source only)

The process ceiling and existing governor remain authoritative. New optional
owners coexist with compact records/migration, PM/spectral storage, tree lanes,
worker DFS scratch and the communication arena; none supplies a separate cap.

| Owner | Scale / element | Capacity and lifetime | Admission / coexistence |
| --- | --- | --- | --- |
| Gravity runtime Tree/PM split caches | six double lanes + uint64 generation = 56 bytes/local particle; two PM identity scalars | Only hierarchy; retained across fine solves, invalidated by migration | Both shared gravity estimates and phase reservations include 56N; actual vector capacities reconciled. Coexists with existing total-force history. |
| Scheduler hierarchy | existing byte/uint64/size_t/candidate/member lanes; typically 33 bytes/particle live, ABI dependent | Enabled max rung >0; exact bin occupancy reserves; retained capacities governed | Initial process reservation models lanes/bin tables; subsequent bin changes use RetainedCapacityTransaction. Migration uses existing generic bounded transaction. |
| DMO bin mirror | one uint8/local particle | Allocated for hierarchy; uniform species/owner/drift remain scalar | Existing state owner, initial and block capacity transactions; not scheduler truth. |
| Timeline | 8193 doubles maximum, 65544 bytes | Stack during one hierarchical block | PhaseResident reservation precedes entry; coexists with solve and endpoint publication. |
| Tree refit certificate | four root doubles plus unwrap/frame/identity scalars | Solver/coordinator lifetime | O(1); existing coordinate and node lanes reused; no second tree/N scratch. |
| Hot worker counters | `sizeof(TreePmTraversalCounters)` per worker, alignas(64), scalar counters and summed/max block time | Current force / retained bounded worker workspace | Shared worker estimate uses actual hot sizeof; existing bounded DFS scratch and arena. |
| Optional spatial-worker scratch | `kTreePmResidualSpatialCounterBytes` per worker, separate alignas(64) work/target arrays with 64 bins each | Allocated/grown and zeroed once per enabled solve; retained capacity reported separately | Shared worker preflight conservatively includes this width even when disabled, preserving the previous admission envelope. Default-disabled execution has no histogram allocation or updates; no per-target allocation. |
| Spatial EMA | two 64-double arrays plus bounded current-force/planner views | Coordinator lifetime | Fixed feedback estimate/report entries; no population/run-length owner. |
| Cadence detail | at most 256 records with three control-path strings each | Recent history + bounded final report copy | Gravity phase model includes warm capacity; actual record/string capacities reported. |
| Operational events | 256 recent records; bounded merged report view | Profiler lifetime | Source capacity model includes strings/map buckets/nodes and lifetime scalar counts. Per-event prospective admission is not implemented. |
| Stage audit | at most 128 stage-name strings plus lifetime scalar | Workflow report lifetime | Bounded; small control metadata uses existing process/allocator reserve, not a new allocator or cap. |

Hierarchy uses existing uint32 identity scratch for all-source drift, explicitly
admitted by the workspace arena. Reports include its governed commitment and
exclude it from baseline owned bytes, avoiding double reservation. No full tree
or PM mesh per worker, gas-zero owners or homogeneous cold lanes are restored.

Finite resulting-rank fit is still stronger than compact planner transient fit.
P6 admits new scheduler/force owners at materialization and reconciles runtime
process baseline, but does not add a predictive full resulting-rank live-set
model to SFC cuts. Oversized future ranks can still be rejected at runtime.
Qualification must exercise simultaneous planner, migration, hierarchy, PM,
Tree, arena and output/restart peaks under finite budgets; allocator/backend
unknowns retain existing reserves. These are limitations, not a 256^3 fit claim.

## Optional TreePM coefficient owner (2026-10-10)

The Gaussian pair table has one `kPersistentCache` reservation per enabled
`TreePmCoordinator`, admitted before its sole allocation and committed only after
construction. Its report entry uses the physical `sizeof` capacity and marks
`governed_commitment`; failed allocation releases pending admission. Destruction
frees the object before releasing its commitment. Disabling lookup on a later
call retains the bounded table as an explicit cache high-water. No table exists
on the default analytic path. Split/softening changes require no new allocation.

Whole-gravity preflight includes the optional capacity. The workflow phase lease
excludes it because the coordinator holds the separate persistent commitment,
just as the communication arena has its own owner. No per-particle storage is
added by this optimization; existing optional hierarchical split-force lanes
remain governed by their existing transactions. This does not certify a large
DMO simulation's whole-process fit.

## Node-aware admission, owner telemetry and stage transactions (CHUI Part 1, 20261103)

Status: source only, NOT RUN.

**Rank vs node.** `MemoryGovernor` remains the single rank-local authority.
`core::computeNodeMemoryAdmission()` (`include/cosmosim/core/node_memory_admission.hpp`)
derives its hard limit when `parallel.node_memory_admission = shared_equal`:

    node_limit  = min(node_memory_budget_bytes [if set], cgroup v2/v1 limit,
                      SLURM_MEM_PER_NODE, physical RAM)
    node_usable = node_limit - node_memory_reserve_bytes      (once)
    rank_share  = floor(node_usable / ranks_on_node)           (MPI_COMM_TYPE_SHARED size)
    hard_limit  = min(rank_share, process_memory_budget_bytes if non-zero)

`sum(hard_limit) <= node_usable` for every co-located rank set, independent of
how unevenly ranks allocate, so no node-level shared counter or MPI window is
required. Rank-local external reserves (`process_*_unknown_reserve_bytes`) stay
inside each rank's share: nothing is subtracted twice. Discovery is
best-effort and never throws; if no limit is known the admission refuses and
asks for an explicit node budget. Assumptions: a discovered cgroup/scheduler
limit is shared by all co-located ranks (conservative for per-task cgroups);
physical RAM is optimistic on shared nodes. The 36/42 GiB numbers of the
historical 48 GiB workstation are only an illustrative profile; the rule is
size-independent (unit tests use synthetic 8/16/48/128 GiB nodes). Serial and
non-MPI builds use `ranks_on_node = 1` and never call MPI. The selection is
reported in the `memory.node_admission` profiler event.

**Owner telemetry.** `MemoryGovernor::ownerTelemetry()` returns, per owner label
(first 63 labels, then `<other>`), live bytes, peak live bytes, reservation and
rejection counts. The table is fixed-size; reservation never allocates for it.
Accounting remains an estimate of CHUI-owned allocations, not OS residency
(`observeProcessMemory()` stays diagnostic evidence).

**Batch-scoped hydro scratch.** CTU (`3 x HydroConservedState` per touched
cell) and dual-energy (`2 x double` per touched cell) lanes are planned by
`HydroScratchBuffers::workspaceBytes/prepareForActiveSet` when the batch policy
requests them (`HydroCoreSolver::requiredWorkspaceFeatures`), are sized by the
active patch/batch, and are never persistent all-mesh copies. The batched
solver refuses an unplanned lane rather than growing after admission. PPM
stages up to six primitive stencil cells per face (was four).

**Stage transactions.** `core::StageTransaction` admits all transient bytes
before work (refusal is recorded, not thrown, so every rank still votes),
classifies preparation failures (numerical retry / model downgrade / resource /
invalid state), performs exactly one collective vote, commits pending
reservations before the authoritative mutation, and releases everything exactly
once on every path. Wired today around the fixed-grid hydro commit (staged
conserved state -> `SimulationState`). Feedback touch-log capacity is admitted
before each deposition batch (old + new storage coexist during growth).

**Known ungoverned allocations (unchanged by Part 1):** FFTW/MPI/HDF5 internal
buffers (covered only by the explicit reserves), per-step reconstruction /
Riemann solver objects (a few hundred bytes), the `ownerTelemetry()` result
vector at diagnostic boundaries.

### Node-shared accounting pool (`shared_pool`, corrective pass 20261103)

`shared_equal` gives each co-located rank a fixed equal ceiling: deterministic
and provably non-overcommitting, but it cannot let a heavy rank borrow from a
light one. `shared_pool` keeps the same node limit and single reserve and adds
one node-wide accounting counter:

* `parallel::NodeSharedMemoryPool` allocates one 8-byte MPI-3 shared window on
  node-local rank 0 (`MPI_Comm_split_type(MPI_COMM_TYPE_SHARED)` +
  `MPI_Win_allocate_shared`), collectively, once, at workflow start (main
  thread). Non-MPI builds and single-rank runs use a process-local counter.
* `MemoryGovernor::attachNodePool()` acquires the rank's current accounted
  bytes; afterwards every reservation acquires its bytes from the counter by
  lock-free compare-and-swap (`std::atomic_ref<uint64_t>`, compile-time checked
  lock-free) and every release returns them. Baseline reconciliation of
  already-resident memory is recorded unconditionally (`forceAcquire`), so the
  counter reflects reality even if it exceeds capacity, in which case all new
  reservations on the node fail until memory is released.
* Invariant: `sum over co-located ranks of governed bytes <= node_usable`
  whenever reservations succeed; a single rank may use up to `node_usable`
  (capped by a non-zero `process_memory_budget_bytes`).
* No MPI call or collective per reservation (works with
  `MPI_THREAD_FUNNELED`), no waiting, hence no deadlock; cost is one atomic CAS
  at reservation boundaries, never in per-cell/per-face loops.
* Failure is explicit (`MemoryAdmissionError` with `node_pool_exhausted`, used
  and capacity bytes) but, under contention, NOT run-to-run reproducible; use
  `shared_equal` when admission must be deterministic.
* What is governed: everything that goes through `MemoryGovernor` (hydro batch
  incl. CTU/dual lanes, AMR sweep workspace incl. the auxiliary lane, AMR
  regrid transactions, feedback touch log, existing owners) plus the declared
  external reserves. Not governed: allocations that bypass the governor and
  library internals beyond the declared reserves. The counter is an
  accounting limit, not OS RSS/PSS; `observeProcessMemory()` remains the
  diagnostic for actual residency.
* Dependency: borrowing that is coordinated with work redistribution (moving
  load away from a rank near the node limit) belongs to the Part 3
  distributed scheduler.

### Shared-pool guarantee conditions and fallback (stabilization pass)

`std::atomic_ref<uint64_t>::is_always_lock_free` alone does NOT prove that
atomics on memory mapped into several processes are atomic across those
processes: the C++ standard only *recommends* that lock-free atomics be
address-free. The node-wide sum guarantee of `shared_pool` is therefore
claimed only when all of the following hold, and the decision is agreed by
`MPI_Allreduce(MIN)` over the world at construction (all ranks use the pool or
none does):

1. Toolchain/ISA allowlist: GCC/Clang on x86-64, AArch64 or POWER64, where
   64-bit lock-free atomics compile to address-free hardware instructions
   (`lock cmpxchg`, LSE/LL-SC, `ldarx/stdcx`) on cache-coherent memory, and
   `std::atomic_ref<uint64_t>::is_always_lock_free`.
2. `MPI_WIN_MODEL == MPI_WIN_UNIFIED` for the shared window (direct load/store
   is coherent with the window's public copy).
3. The counter from `MPI_Win_shared_query` is 8-byte aligned.
4. Initialization: node rank 0 stores 0 (release), then `MPI_Win_sync` and a
   node barrier before any rank reserves. Updates are acquire-release CAS loops;
   no rank ever waits on another.

If any condition fails anywhere, the window is freed, `established()` is false,
every rank keeps the deterministic `shared_equal` ceiling, and the
`memory.node_admission` event reports `node_pool_established = false` with the
reason. Single-process runs use a process-local counter (always correct).

Unequal ranks: the pool only bounds the sum; a rank may hold up to the usable
node budget (capped by a non-zero process budget). A refused reservation
acquires nothing (no rollback needed); releases and governor destruction
return exactly the bytes held. Teardown: governors are destroyed before the
pool; the pool's `MPI_Win_free` is collective and is skipped during exception
unwinding so a failing rank cannot hang (MPI_Finalize/Abort reclaims it).
Baseline reconciliation records resident memory even above capacity, which
then refuses new reservations node-wide until releases happen (truthful
accounting, not hidden overcommit).

Precise guarantee now provided:
* `off`: per-rank `process_memory_budget_bytes` (0 = unlimited), as before.
* `shared_equal`: for every node, sum over co-located ranks of governed bytes
  `<= node_limit - node_reserve`, deterministic.
* `shared_pool` with `established()`: the same node-wide bound, with
  unequal per-rank use; non-deterministic admission order under contention.
* `shared_pool` without `established()`: exactly the `shared_equal` guarantee.
All of these bound GOVERNED accounting (CHUI allocations that go through the
governor plus declared external reserves), not OS RSS/PSS.
