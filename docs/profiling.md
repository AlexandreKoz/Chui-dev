# Profiling and benchmark workflow

## Attempt 1 staged experiments — 20261012

The frozen-source driver now defaults to a single mesh/leaf/block/policy/kernel
candidate. Mesh 64/96/128, leaves 4/8/16/32 and blocks 16/32/64/128/256 should
be compared in stages; candidate count is checked before launch (default cap
40). `--asmths` and `--rcuts` add bounded split studies. `--maximum-angle 0.25`
and 0.5 are separate adaptive experiments. Kernels are `analytic`, `lookup`,
`full_residual`, `simd`, `simd_full_residual` and `radial`. Reuse/refit are
independently controlled with `--reuse on|off` and `--refit on|off`; unchanged
frozen-state reuse does not represent moving production-source reuse.

Production uses the same target block setting through
`numerics.treepm_residual_block_size` (default64), with leaf size controlled by
`treepm_tree_max_leaf_size`. The typed C++ experiment generator writes isolated
normalized decks. Preserve runtime's emitted `normalized_config.param.txt`,
hash, build/source identity and event logs for every comparison.

Benchmark JSONL/CSV adds source preprocessing, explicit refit/rebuild and node
geometry preparation timing, compiler/CPU/source/provider metadata and selected
mode target counts. Production events expose requested flags, block/leaf and
actual eligible target counters. `residual_lookup_targets` counts table-eligible
targets with a table installed, not interpolation hits. `simd_leaf_targets`
counts targets selecting AVX2 staging, not vector occupancy. `skipped_geometry`
in benchmark output (`skipped_geometry_evaluations` in runtime) identifies
unevaluated COM/envelope/softening guards. Existing exact interaction/node
counts and worker imbalance remain; no new per-pair clocks/atomics were added.

Each physical split needs its own short-force oracle. The existing independent
total-reference manifest gate remains mandatory for external total comparisons;
matching short forces or fewer interactions alone cannot certify PM+Tree.
No new timing was collected. See the [exact staged commands and memory gates](repair/treepm_deep_hotpath_attempt1_closure_20261012.md).

## OpenMP concurrency contract

`ProfilerSession` and `CounterRegistry` use per-thread shards. Worker threads
record scopes, byte counts, and counters without a global mutex on each
operation. Control-path operational events use a session-wide bounded ring
under a short mutex; no event is recorded per node or pair. The first use by a thread registers its shard under a short mutex;
report generation and reset are quiescent-boundary operations that merge shards
deterministically after worker activity has joined. `AllocatorStats` uses atomic
accounting. Callers must not request a report/reset while worker threads are
still mutating the same session.


CosmoSim benchmarks are lightweight hooks for performance trend visibility, not proofs of correctness.

## Benchmark directory structure

- Core/system hooks: `bench/bench_*.cpp`
- Module-focused hooks: `bench/gravity`, `bench/hydro`, `bench/amr`, `bench/io`, `bench/mini`
- Reporting utility: `bench/reporting/bench_report.hpp`
- Baseline size presets: `bench/baselines/benchmark_sizes_v1.txt`

## Build and run

```bash
cmake --preset cpu-only-release
cmake --build --preset build-cpu-release --target bench_config_parser bench_profiling_overhead bench_pm_solver bench_tree_pm_coupling bench_parallel_decomposition_exchange
./build/cpu-only-release/bench_config_parser
./build/cpu-only-release/bench_profiling_overhead
./build/cpu-only-release/bench_pm_solver
./build/cpu-only-release/bench_tree_pm_coupling
./build/cpu-only-release/bench_parallel_decomposition_exchange
```

## Reporting expectations

Each benchmark report should include:

- build type/preset,
- hardware summary and thread count,
- enabled feature flags,
- setup vs steady-state timing where practical,
- throughput/rate metrics or effective bandwidth proxies for memory-sensitive paths.
- cache/build reuse counters for persistent infrastructure where applicable (for example PM FFT plan cache build count, warmup vs measured solve counts).

Do not claim scientific correctness from benchmark throughput.

## Distributed TreePM phase-2 performance hardening notes

- `bench_pm_solver` now reports `plan_cache_size` and `plan_build_count`; in a stable slab layout `plan_build_count` should stay flat across measured solves.
- `bench_tree_pm_coupling` runs repeated measured iterations and exercises PM cadence reuse (`refresh_long_range_field=false` on selected iterations) so PM reuse impacts are visible in timing deltas.
- `bench_parallel_decomposition_exchange` remains the communication planning baseline for ownership/decomposition costs; compare alongside TreePM coupling output when evaluating communication overhead.

## Profiling discipline

- Keep hot-path instrumentation low-overhead and explicitly gated.
- Preserve deterministic behavior when profiling is disabled.
- Compare against prior baselines using the same preset and workload size.

## TreePM residual counters, timers, and geometry telemetry

The gravity step event derived from `TreePmDiagnostics` / `TreePmProfileEvent`
reports:

- **Pair counters.** `local_pair_evaluations` and
  `incoming_remote_pair_evaluations` are non-overlapping residual bundles.
  `remote_pair_evaluations` is a compatibility alias of the incoming-remote
  value; `total_pair_evaluations` is the residual total and equals the exact
  sum of the two local/incoming components for the same solve.
- **Wall timers.** `tree_wall_ms_recent` is short-range tree wall time alone;
  `pm_wall_ms_recent` is `pm_profile.total_ms` alone (PM phase entry through
  tree short-range start, including long-range refresh when taken).
  `incoming_request_decode_validation_ms`,
  `incoming_remote_target_compute_ms`, and
  `incoming_response_encode_pack_ms` split incoming request
  decode/validation/hash, pure incoming tree compute (force evaluation of
  validated targets against this rank's tree only), and response
  encode/pack respectively;
  `protocol_validation_ms`, `protocol_consensus_ms`, and
  `response_exchange_ms` split response layout/payload sizing, failure
  consensus, and the response `MPI_Neighbor_alltoallv` call respectively.
  `let_remote_traversal_ms` is a compatibility alias of
  `incoming_remote_target_compute_ms`. OpenMP provenance
   fields `openmp_compiled`, `openmp_configured_workers`, and
   `openmp_observed_workers` come from `TreePmDiagnostics`
   (configured is the planned configured/runtime maximum; observed is the
   maximum actual team size recorded across the local, distributed overlap,
   and incoming-target regions during the last short-range solve; serial
   execution reports one).

  `residual_local_target_count` / `residual_incoming_target_count` count
  targets evaluated in each bundle for that solve;
   `residual_worker_scratch_high_water_bytes` is the retained contiguous
   worker-stack capacity (`worker_count * (1 + 7 * kMaximumTreeDepth) *
   sizeof(TreeLocalIndex)`); worker counter bundles are `O(worker_count)` and
   deterministic block floating diagnostics are `O(ceil(active_targets/64))`.
   The matching MemoryGovernor estimates are
   `gravity.estimate.treepm_residual_worker_scratch` and
   `gravity.estimate.treepm_residual_block_diagnostics`.

- **LET high-water.** `let_wire_buffer_high_water_bytes` is the four reusable
  short-range payload buffers; `let_known_workspace_high_water_bytes` is the
   known workspace peak (wire + structured + counts/masks/accumulators +
   metadata + transient codec). CHUÍ-owned retained storage is read from

  actual `vector.capacity()`; the transient codec term remains a modeled
  conservative upper envelope. These are capacity high-waters, not bytes
  communicated and not preflight estimates.
- **Domain geometry.** `domain_geometry_source_generation`,
  `current_gravity_source_generation`, `domain_geometry_fresh`,
  `domain_geometry_refreshed`, `domain_geometry_fallback_used`,
  `domain_geometry_fallback_reason` (string name of
  `TreePmDomainGeometryFallbackReason`),
  `domain_geometry_refreshed_leaf_count`,
  `domain_geometry_out_of_seed_range_source_count`, and
  `domain_geometry_uncovered_source_count` record the derived top-domain
  freshness/refit/fallback path. The published refit is seeded from the
  stable decomposition-local seed leaf set (not the previous published
  result), and an ownership commit invalidates freshness until reinstall.
  They are observational; restart/snapshot schema is unchanged.

## PM routing runtime events

The `gravity.treepm` runtime event keeps total communication traffic separate from
resident routing memory. Its PM payload includes logical density/force/potential
route counts and peer counts, total MPI bytes sent/received, MPI wait time, and
four distinct capacity high-waters: `pm_routed_send_buffer_high_water_bytes`,
`pm_routed_receive_buffer_high_water_bytes`,
`pm_routed_combined_buffer_high_water_bytes`, and
`pm_routed_workspace_high_water_bytes`. The combined value is the simultaneous
capacity of the two reusable wire buffers; the workspace value additionally
includes retained rank-scale routing count/displacement metadata and the reused
per-peer packing/response cursor. The capacity model intentionally leaves 64 KiB
below the 128 MiB/rank M1A engineering ceiling. These are capacity/high-water
metrics, not aliases for bytes communicated over the phase.

The DMO process preflight also reports hierarchical scheduler memory as three
separate values: logical live bytes, retained owned capacity, and historical
retained-capacity high-water. Candidate labels remain transient compatibility
inputs and do not contribute population-scale scheduler storage.

## M1B memory-governor telemetry

The existing profiler `memory_report` JSON object may contain a `governor`
object. It reports the authoritative process hard limit, reconciled baseline
owned bytes, configured external/planned reserves, governed committed bytes,
pending reserved bytes, raw and safety-adjusted accounted demand, raw headroom,
pressure (`green`, `amber`, `red`, or `trip`), historical committed/reserved/
accounted high-waters, and the reservation rejection count. This is an additive
view in the existing profiler format, not a second profiling stream.

Do not add `memory_report.totals` to governor `committed_bytes`: governed
physical blocks remain visible as ordinary `MemoryEntry` capacity for ownership
reporting, but are excluded from `baseline_owned_bytes` using the
`governed_commitment` reconciliation marker. The governor snapshot is the policy
view; the ordinary memory entries remain the ownership/capacity view. M1B does
not claim RSS/PSS reconciliation.

## Hydro runtime events

The reference workflow emits `hydro.conservation` events from the Godunov hydro stage when a profiler session is
present. Payload values are volume-integrated local diagnostics over the cells updated by that stage. The event includes
before/after mass, flux/source/floor mass deltas, residuals for mass, momentum x/y/z, total energy, and derived internal
energy, total-energy source/floor deltas, `internal_energy_floor_count`, and the tolerances used by the CPU closed-box
regression. Source-term deltas are reported separately from flux deltas so gravity/expansion work is not classified as
a flux-conservation error.

The solver also records these same totals in `HydroProfileEvent::conservation`, alongside the existing face count,
fallback counters, bytes moved, and stage timings.

For production AMR hydro, the reference workflow emits `hydro.amr_production_stage` after the AMR hydro synchronization
point. The payload includes patch and active cell/face counts, `flux_register_entry_count`, reflux corrected cell count,
corrected mass, corrected momentum x/y/z, corrected total energy, corrected internal energy, complete register count,
and skipped register counts for incomplete, area-mismatched, or missing-target registers. These diagnostics are
observational only; flux-register ownership remains in the AMR hydro synchronization path and is not persisted as
restart truth.

## Initial-condition ingestion events

After initial state construction, every rank emits one
`io.ic_ingestion.summary` event in subsystem `io.initial_conditions`. Its
payload records the single `provenance_authority` selected for the import and
includes:

- `files_assigned`
- `chunks_assigned`
- `logical_metadata_bytes_read`
- `hash_bytes_read`
- `logical_payload_bytes_read`
- `converted_payload_bytes`
- `bytes_serialized`
- `bytes_sent`
- `bytes_received`
- `manifest_metadata_bytes_communicated`
- `records_read`, `records_converted`, and `records_routed`
- `source_file_open_count`, `source_dataset_open_count`
- `full_file_hash_pass_count`, `source_identity_validation_count`
- `routing_batch_count`, `reader_batches_assigned`, `reader_records_assigned`, `reader_record_imbalance`
- `main_exchange_count`, `exact_audit_exchange_count`, `distributed_id_audit_round_count`
- `logical_consensus_phase_count`, `routing_logical_consensus_phase_count`
- compatibility aliases `collective_phase_count`, `routing_collective_phase_count`
  (logical phases only)
- `mpi_collective_call_count`, `routing_mpi_collective_call_count`,
  `nonrouting_mpi_collective_call_count`
- `mpi_allreduce_call_count`, `mpi_bcast_call_count`,
  `mpi_gather_call_count`, `mpi_gatherv_call_count`,
  `mpi_alltoall_call_count`, `mpi_alltoallv_call_count`
- `collectives_per_million_records`
- `wall_time_nanoseconds`
- `peak_staging_bytes`
- final local particle, gas, star, black-hole, and tracer counts
- `already_partitioned`

`bytes_read` is retained as the checked sum of metadata, SHA-256, and particle
payload reads. Metadata is the decoded logical HDF5 header payload; hashing
counts every source byte read by the assigned hashing rank; payload counts only
datasets actually read, never merely inspected or explicitly dropped fields.

Peak staging is a capacity-based high-water mark for the actual bounded
import/routing workspace, not the final authoritative owner-local state. It
includes all simultaneously live coordinate, velocity, mass, ID, gas, star,
black-hole, tracer, `ParticleRecord`, nested per-peer, flattened exchange,
count/displacement, decode, coverage, and ID-reconciliation buffers. The exact global duplicate-ID audit
uses rank-local sorted temporary runs and a bounded-memory external merge; disk
run bytes are not RAM staging and are removed before import returns. For a distributed fixture, the required evidence is that each nonempty source
file has one stable payload reader/session. In the default verified-identity mode, complete-file
SHA-256 work is one inspection pass per source file plus stable file-identity validation around
ingestion; explicit strict-full-rehash mode performs one additional completion SHA-256 pass.
The hash-pass count is independent of batch count, each source chunk is assigned once, reader-record imbalance is reported as the maximum minus minimum assigned-record count across ranks, main exchanges
scale with routing batches rather than source chunks, each source ID balances
against one final ID, and no rank allocates authoritative arrays sized to the
global particle count merely because MPI is enabled. `main_exchange_count` and the compatibility
`routing_collective_phase_count` are global protocol counters recorded on rank
zero; byte counters remain rank-local and may be reduced by the caller. The
compatibility collective fields count logical rank-consistent protocol phases,
not raw MPI calls. Actual production communicator calls are counted centrally by
the `mpi_*_call_count` fields. Their per-kind sum equals
`mpi_collective_call_count`; `routing_mpi_collective_call_count` is exactly 19
calls per successful routing batch in routing protocol version 3 (11 consensus votes,
three coverage reductions, two `Alltoall`/`Alltoallv` pairs, and one exact
reconciliation reduction). Fixed discovery, manifest, final audit, and
finalization calls are reported by `nonrouting_mpi_collective_call_count`.
For routing protocol version 3 the successful-path identity is:

```text
routing_mpi_collective_call_count
  = 19 * routing_batch_count

nonrouting_mpi_collective_call_count
  = 37
  + (validate_runtime_cosmology ? 1 : 0)
  + source_file_count
  + 10 * distributed_id_audit_round_count
  + mpi_bcast_call_count
```

The `mpi_bcast_call_count` term is explicit because length-prefixed metadata
broadcasts may require more than one 64 MiB payload chunk. The distributed MPI
acceptance test checks both identities rather than applying an arbitrary loose
ceiling. `collectives_per_million_records` uses the globally routed record
count and is zero when no records are routed.

These counters are scalability evidence, not a substitute for scientific
validation. Exact distributed duplicate-ID, count, mass, ownership, provenance,
and sidecar checks run separately before the workflow accepts the state.

## Workflow hooks for documentation/scaffolding changes

Documentation changes still need auditable developer workflow checks:

- integration doc-scaffold smoke test (`integration_docs_scaffold`)
- docs reference scan benchmark (`bench_docs_reference_scan`)

These hooks ensure core docs remain present and internally referenced as the codebase evolves.

## Star-formation profiling

`bench_star_formation_spawn` profiles three explicit adaptive-model phases: `eligibility_no_births`, `eligibility_sparse_births`, and `dense_birth_plan_and_batch_append`. Each line reports `cells_per_second`, `births_per_second`, `bytes_allocated_per_step`, `allocation_count_per_step`, and `peak_temporary_memory_upper_bound_bytes`. The peak is an upper bound because the tracker includes storage retained by the appended authoritative particle and stellar sidecars as well as temporary planning storage. The implementation uses one particle resize, one star-sidecar resize, one species-index rebuild, no collective inside the cell loop, no full existing-particle-ID scan, and no allocation node per candidate or newborn ID. Compare allocation counts across cell-count sweeps: a fixed number of vector/sidecar growth allocations is acceptable, while allocation growth proportional to candidate count is a regression.


## Star-formation and effective-ISM benchmarks

`bench_star_formation_spawn` reports no-birth, sparse-birth, and dense plan/append throughput plus total allocations, allocated bytes, and temporary-memory upper bounds. Exact ID precommit uses sorted contiguous batches, avoiding one heap node per birth. `bench_effective_multiphase_ism` separately reports table initialization, lookup, direct-equilibrium evaluation, and hydro-closure throughput. Table lookup reconstructs neither cooling tables nor EOS state and allocates nothing per cell.

## M1C-1 process-memory reconciliation fields

When the runtime memory report is attached, profiling emits a `process_memory`
object with `current_declared_residency_bytes`, the compatibility alias
`known_accounted_bytes`, optional current RSS, optional peak RSS, optional PSS,
non-negative `unexplained_resident_bytes`, and an optional observed/known ratio.
Current declared residency excludes inactive future output/restart policy
headroom and uncommitted reservations; governor `accounted_bytes` remains the
separate conservative admission/policy quantity. Unavailable OS measurements
are JSON `null`, not zero.

A `distributed_process_memory` object records rank count and per-metric local,
global-sum, rank-max, rank-mean, and max/mean imbalance for governor-accounted
demand, current RSS, peak RSS, and communication high-water. OS-derived metrics
are valid only when all ranks provide a value. These fields are observational;
the deterministic memory governor remains the allocation authority.

## M2D runtime topology telemetry

The reference workflow emits a `runtime.topology` profiler event containing MPI
world size/rank, MPI node-local size/rank where available, OpenMP compile
status, requested thread count, and configured thread count. These fields are
intended for rank/thread topology comparisons and memory/work imbalance
interpretation. They are diagnostic evidence only and do not alter numerical
configuration or restart state.

Optional science-diagnostic pressure handling emits
`analysis.memory_pressure_deferral` with current pressure plus light/heavy due
and pending state. When a previously deferred diagnostic later executes on a
non-Red analysis hook outside its original cadence, the profiler emits
`analysis.memory_pressure_catchup` with light/heavy catch-up flags. These events
make bounded deferral and starvation recovery auditable without promoting the
pending bits into restart or scientific state.

## M2D-1 optional cadence evidence

`analysis.memory_pressure_deferral` includes pending state and missed counts. `analysis.memory_pressure_catchup` records class, first/latest due step, missed/coalesced count, actual execution step, and `historical_state_replayed=false`. `analysis.optional_cadence_dropped` records pending work discarded at segment termination; `analysis.optional_cadence_restart_policy` records the previous checkpoint's cadence summary and the nonpersistent restart policy. These are operational provenance events, not scientific outputs or an exact-cadence guarantee.

## M2D diagnostic admission evidence

The existing governor snapshot/rejection counters and `RuntimeEvent` stream
are the only memory-reporting authority. `analysis.memory_pressure_deferral`
records class, requested bytes, headroom and reason; catch-up/terminal events
retain actual versus due epochs. The model distinguishes FFT mesh coexistence
from reduction workspace, and the existing process RSS/PSS report remains the
measurement authority. Compare identical numerical configurations and thread
counts; do not equate a declared reservation with measured RSS. A required
health preflight is owner-managed and does not constitute a second physical
allocation.

### PM spectral-operator cache telemetry

`PmProfileEvent::spectral_operator_rebuilds` counts reconstructions of the
small periodic PM spectral axis/window metadata used to derive the Poisson
operator. M48-09 no longer caches a full Poisson/gradient operator array.
Cosmological `scale_factor`, TreePM split scale, and gravitational constant
remain force-state/collective inputs but are applied directly and are not axis
metadata cache dependencies. Changing only those values therefore leaves this
counter at zero after the axis metadata has already been built for the same
mesh, box lengths, assignment scheme, and deconvolution policy.

## MPI global diagnostics ownership

Production `AnalysisRuntime` distinguishes rank-local operational state from shared science/run-health diagnostics. A due shared diagnostic is first prepared locally on every rank; local preparation failures are coordinated before any diagnostic reduction. Global health counters and SFR bins use communicator sums, invariant booleans use communicator-wide AND semantics with failing-rank counts, and angular-momentum **vectors** are summed component-wise before any norm is reported. Slice quicklooks reduce density sums plus sample counts before forming the global average; projection grids are element-wise sums.

Only rank 0 publishes ordinary shared diagnostic JSON/CSV files after those reductions complete, using `.part` then rename transactional publication. Memory reporting preserves publisher-rank local ownership while also exposing distributed rank-sum, rank-max, and imbalance fields; RSS/local ownership is never relabeled as one global process value. Large particle state is not gathered to the publisher.

A correct MPI power spectrum requires a global density field and global/distributed FFT. Until that contract exists, distributed diagnostics record `unsupported_under_mpi_requires_global_density_fft` and publish no averaged rank-local spectrum. Serial power-spectrum behavior is unchanged. Floating reductions retain the repository tolerance-equivalence contract rather than promising bitwise rank-count invariance; integer counts and boolean outcomes are exact.

## P6 force observability and bounded retention (unqualified)

Each `gravity.treepm_let` event reports source preprocessing, identity/hash,
softening resolution, Morton order, topology, multipoles, node count, depth,
reuse attempt/result/reason, membership-check time and refit time. Full rebuild
and refit timings are inclusive; multipole time is a child. They must not be
added again to inclusive Tree residual phase time.

Local-owned and incoming-target counter families are separate. Accepted
internal multipoles and accepted leaves are separate; `accepted_nodes` remains
a compatibility sum. Direct pairs count actual kernel evaluations. Node visits,
opened nodes, cutoff-pruned nodes/skipped pairs, selected/relative MAC,
maximum-angle/strict-envelope, softening, near-node and cutoff-containment
rejections are exposed. Rejection families overlap; opened nodes count unique
descents. Pair identity remains `combined = local + incoming`.

Aligned worker bundles avoid hot-path atomics and intentional false sharing.
Current-force nonempty worker-region summaries report target min/max, visit,
pair and multipole maxima, work-time min/max/sum and region count. Means can be
computed from aggregate totals / region count. Work time sums logical-block
durations, including target preparation, force-slot writes and local block
diagnostics; scheduling waits and MPI are excluded. P6.1 uses two clock reads
per existing 64-target block in all three traversal regions and none in the
target evaluator. Existing `worker_work_ms_*` event fields summarize these
block-duration sums. Native `local_traversal.block_work_ms_max` and
`incoming_traversal.block_work_ms_max` expose the largest block duration for
tail inspection; they add no serialized event field or block table.

Optional fixed 64-bin work/target histograms live in separate aligned worker
scratch, prepared once per enabled solve and reduced in fixed worker/bin order
before the unchanged global spatial history update. Default-disabled P5-W
does not allocate, clear, update or reduce histogram scratch. Hot scalar
counters retain cache-line isolation. This reduces structural instrumentation
overhead; no benchmark or validation evidence is claimed.

PM `total_ms` is inclusive. `assign_ms` includes deposition, routing and
normalization; `density_routing_wait_ms` is its nested MPI wait. Forward FFT and
spectral operator are distinct. Each axis exposes factor sweep, inverse FFT and
normalization; these lie inside `gradient_ms`. `fft_inverse_ms` also includes
an optional potential inverse. Halo exchange is separately labeled inside PM
phase; interpolation is inclusive. FFT transpose counters/timers may be nested
inside FFT time. This is not an additive ledger.

LET retains candidate/export/import peers and targets, geometry fallback reason,
domain/graph cache hits, wire bytes, communication arena capacity/high-water,
and exported targets per requested target. Fine hierarchical events use
`gravity.short_range_force`, `short_range_only=true`, zero PM solve/reuse counts,
and the actual force-evaluation time/scale factor. Cadence records add
`active_force_targets` and `short_range_only`; their legacy kick-count field is
zero for hierarchical force-only surfaces. In the reference force records it
retains its historical planned-kick-target meaning; traversal target counters
and `active_force_targets` are the actual force work authority.

Recent operational events are capped at 256, force-cadence records at 256 and
stage traces at 128 (16 canonical eight-stage sequences). Every stage is still
audited. Lifetime severity totals and fixed-name work/stage counters survive
eviction. Operational-report JSON is version 2: `event_count` is lifetime total,
with retained/evicted counts and retention limit. Ordinary profiler JSON remains
version 1. `ProfilerSession::events()` now returns recent ordered detail only;
consumers needing lifetime counts must use severity totals/counters. Retained
payload/string capacity is included in memory reports; prospective admission of
each dynamic event payload remains a follow-up limitation.

## Frozen-source split/leaf/block sweeps

`bench_tree_pm_sweep` and `python3 tools/treepm_sweep.py` provide external
benchmark controls without another simulation-input language. The driver reads
`id,x,y,z,mass,epsilon` CSV in increasing unique ID order, or generates a fixed
synthetic fixture. Full-population force refreshes are timed; fixed sampled IDs
are used for accuracy. It always refreshes PM, rebuilds the tree by default,
and reports warmups separately from repeated measured evaluations.

Sweep mesh `64,96,128,192`, leaf `4,8,16,32`, optional block/thread axes, strict or
adaptive acceptance, analytic or opt-in lookup. Cell-valued `asmth`/`rcut` derive
the candidate physical split from its own mesh. CSV softenings are preserved;
`--epsilon` applies only to generated sources. `--history bootstrap` uses a
strict analytic same-snapshot total force and labels it explicitly; `fallback`
supplies no history. Neither reconstructs the simulation's historical KDK state.

JSONL records distinguish wall, PM kernel phases, tree build/multipoles and
traversal; wall is the authoritative end-to-end measurement. Counts, actual
workers, solver-retained capacities, current RSS and process peak RSS are
reported. Driver arrays contribute to RSS, not solver-owned capacities. Fixed
force CSV contains total, Tree-only, and independently accumulated matching
short residuals. Relative L2, P95/P99/max normalized errors, absolute tails and
small-force absolute maxima accompany each oracle comparison. The normalization
floor is explicitly `1e-3*RMS(reference)`.

The Python workflow saves exact commands, binary/source/target/reference hashes,
stdout/stderr, candidate status and CSV/Markdown tradeoff reports. Each candidate
runs in a fresh process. A failed backend/geometry is a saved failure, never a
silent resolution substitution. `--ewald-level` is limited to <=4096 uniform-
softening sources and uses the independent existing periodic Ewald reference;
increase the level before accepting its result. Large meshes require FFTW.
The diagnostic is a single-owner benchmark; MPI qualification uses registered
force/workflow tests. See [manual commands](repair/treepm_unified_performance_recovery_20261010.md).

With `treepm_full_mac_diagnostics=false`, events explicitly say
`mac_accounting=fast_guard_short_circuit`, include skipped MAC counts for both
local and incoming targets, and continue to report unique descent and actual
pair counts. Selected-MAC rejection/fallback counts then describe evaluated
MACs only. Full mode retains overlapping failed guards; skipped is never passing.

### Complete snapshot export and reference identity

P1–P5 follow-up export consumes the entire production `.complete` snapshot set,
including all members, rather than interpreting one rank member as a population.
The CSV sidecar binds member/manifest checksums, persistent IDs, epoch, rectangular
geometry, declared dialect and conversion factors. Keep that sidecar with the
CSV and target list. `--fixed-target-ids` preserves the historical oracle IDs;
changing their file after export deliberately fails its checksum check.
`--total-reference-manifest` is mandatory for external total-force comparisons.
Its evidence and frozen-state identities must match; matching-split direct
short-force comparisons are generated separately per candidate. Python >=3.8,
NumPy and h5py are required for these offline workflows. Exact commands and the
external reference sidecar fields are in the
[20261010 repair handoff](repair/treepm_p1p5_qualification_and_test_reliability_repair_20261010.md). No new timings were collected.

## Clustered-box TreePM benchmark (October 2026)

`tools/generate_clustered_dmo_box.py OUT.csv --order random|morton` writes a
clustered 64^3 DMO workload for `bench_tree_pm_sweep --particles OUT.csv`.
The benchmark accepts `--alpha`, `--estimator residual_derivative|relative_monopole`,
`--cutoff-node-policy contain|com_distance`, `--group-walk on|off` and `--maximum-angle` up to 1.0;
`--forces FILE` writes total/tree/direct-residual forces for the selected
targets, so two runs (candidate and `--policy strict --kernel analytic` on the
same mesh) give a total-force error distribution. Report pairs/target,
multipoles/target, traversal and PM milliseconds together with the error
quantiles; the measured table is in `docs/tree_pm_coupling.md`. Hierarchical
runs additionally record `coarse_limiter` and `rank_local_bin_occupancy` in
every `time.hierarchical_block` event of `operational_report.json`.

## Hydro method matrix benchmark (CHUI Part 1, source only)

`bench/hydro/bench_hydro_method_matrix.cpp` (target `bench_hydro_method_matrix`,
built with `COSMOSIM_ENABLE_BENCHMARKS=ON`) runs PCM/PLM/PPM x MH/CTU x
HLLC/HLLE on a periodic N^3 smooth oblique wave and prints, per
configuration: `wall_ms`, `cells_per_s`, `faces_per_s`, `riemann_solves`,
`riemann_fallbacks`, `positivity_fallbacks`, `stencil_degradations`,
`ctu_degraded_states`, `scratch_high_water_bytes`, `moved_bytes`,
`l1_density_error`. Arguments: `[cells_per_axis=32] [steps=8] [batch_cells=0]`.
Compare methods by error at fixed wall time; sweep batch_cells (e.g. 512,
2048, 8192, 16384) as candidates. Targets (not claims): <= 5% overhead of the
default path versus pre-Part-1 with advanced options disabled. NOT RUN.
