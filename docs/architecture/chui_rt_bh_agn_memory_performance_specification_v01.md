# CHUÍ — RT/BH/AGN Memory, Performance and HPC Specification

**Version:** 0.1 · **Date:** 2026-11-26 · **Status:** proposed engineering contract  
**Parent:** [Master design](chui_rt_bh_agn_master_design_v01.md)

## 1. Objective and platform doctrine

Desktop-first means optimize scientific return on constrained hardware, not impose a maximum cluster size or require simplified physics on desktops. Same selected equations and model IDs across CPU, GPU, MPI and supercomputers, subject to documented numerical tolerance. The user may select different physics profiles but hardware does not silently choose them. Minimize wall time at fixed scientific accuracy, within host/node/device/rank budgets and deterministic restart constraints. This specification extends CHUÍ Memory Architecture and Resource Governance v0.1 (2026-08-23); the 48 GiB example has 36 GiB nominal, 42 GiB hard safe simulation ceiling and 6 GiB external reserve; not universal thresholds. Actual machine allocation, cgroup/job limits, runtime libraries and per-node topology define budgets. In HPC account for sums of processes sharing a node, not merely independent rank limits.

## 2. Memory ownership and quantitative scale law

| Owner | Memory class | Scaling and lifetime |
|---|---|---|
| RT canonical | CanonicalPersistent | 4 G b per stored cell for M1 (G groups, b scalar bytes); patch geometry owns identity; persists through restart |
| RT stage arrays | PhaseResident/ScratchArena | stage state and transport/implicit intermediates, bounded by active patches and required stage history |
| RT ghosts/registers | Communication / PhaseResident | ghosted surface and coarse-fine faces, all active groups, stage/time history where necessary |
| Opacity/spectral tables | PersistentCache | table dimensions, one immutable version per memory domain where possible; not per element |
| Chemistry interface | ScratchArena | coupled Jacobian/iteration scratch per active worker or batch; existing species state not duplicated |
| BH/disk | CanonicalPersistent | N_BH * measured field stride; no cost per DM/gas cell |
| BH queries | ScratchArena/Communication | active BH batch, bounded patch intersections and sufficient-statistic reductions |
| AGN event queue | CanonicalPersistent for pending physical events; ScratchArena for plans | O(active/pending events), bounded staged recipient count, no permanent dense BH-cell map |
| AMR changes | PhaseResident | coexistence of old/new state plus prolongation, CT, ghosts and registers |
| Device mirrors | separately budgeted | only necessary resident arrays, explicit ownership and synchronization; host staging reported |

Persistent uniform-equivalent lower bounds (fp64, 3 groups, no ghosts/covered parents/scratch):

| Grid | FLD (1/group) | M1 (4/group) | S_N (48/group) | M1 with 9 groups |
|---|---:|---:|---:|---:|
| 128³ | 0.046875 GiB | 0.1875 GiB | 2.25 GiB | 0.5625 GiB |
| 256³ | 0.375 GiB | 1.5 GiB | 18 GiB | 4.5 GiB |
| 512³ | 3 GiB | 12 GiB | 144 GiB | 36 GiB |
| 1024³ | 24 GiB | 96 GiB | 1152 GiB | 288 GiB |

These are storage estimates, not benchmarked footprints or hardware feasibility guarantees. M1 at 512³ with 6 groups is 24 GiB; an independent photon-number four-moment representation doubles moment storage. For block side p and ghost width h, fully ghosted patch inflation is ((p+2h)/p)^3, e.g. p=8,h=2 gives 3.375; p=16,h=2 gives 1.953125. Do not allocate separate full ghosted copies for every stage/group if bounded tiles or shared storage suffice. Covered parents required for temporal interpolation count explicitly. `N_alloc` includes authoritative cells, retained covered parents and capacity rounding. Ideal per-rank average is total/ranks; hard safety uses rank max and node sum.

For fp64 M1 stage planning, canonical bytes=32 G N_alloc, one full stage copy adds another 32 G N_stage; face flux scratch upper estimate=32 G N_face_active per stored flux set (plus time integrals/metadata). Actual method may fuse/reuse but must prove stage data remains available. A dense local Jacobian with K independent species and G scalar radiation energies plus thermal variable costs 8 (K+G+1)^2 bytes per solve before factorization scratch; if flux/momentum also unknown include those dimensions. Bound concurrency, exploit exact block structure and elimination, avoid allocating one dense matrix per gas cell. Do not introduce silent mixed precision. Float32 and compression require separate range/error budgets including faint radiation, high dynamic range, diffusion cancellations and deterministic ledgers.

## 3. Admission and memory pressure

Before run, estimate maximum canonical state from resolution/refinement, groups, disk profile, BH count bounds, patch shape, halo geometry, scheduler histories, MPI and external libraries. At runtime reserve all large allocations before commit, including AMR old/new state and GPU mirrors. Track `committed_bytes`, `reserved_bytes`, measured RSS/PSS, device usage, capacity high-water, per-owner/per-rank/node maxima and expected lifetime. Overflow-check all products/sums and MPI counts. Reservations include failure/rollback workspace. No OOM as normal control flow. Phase-aware schedule reuses scratch across mutually exclusive operations; do not retain `N_cells*N_groups*N_faces` records when direct accumulation sufficient.

Pressure response order: defer optional diagnostics and nonessential concurrent work; shrink batch/exchange size within certified range; release unproductive caches; use explicitly certified numerically equivalent low-memory execution; then reject/defer according to model contract or fail before allocation. Never change group count, c_hat, angular backend, physical feedback mode, or required resolution silently. Mandatory regrid admission includes transient coexistence; optional refinement within declared policy records skipped regions and residual resolution indicators. An accepted plan must have sufficient memory for completion. No active physical state spilling to disk by default. Checkpoint output streams rather than building duplicate full state.

## 4. CPU execution

RT regular stencil kernels: patch-local SoA, aligned lanes, spatially coherent ordering, explicit vectorization diagnostics, thread-local scratch, first-touch/NUMA binding. Benchmark group-major versus cell-tiled layout for transport and local coupling; maintain one canonical format and avoid whole-mesh transposes per substep. Reuse derived closure tensors only if measured time savings justify bandwidth and bytes; compute cheap M1 algebra on demand. Tile face fluxes and preserve unique flux ownership. Use per-thread accumulation/reduction to avoid hot atomics; deterministic mode fixes summation order, production mode declares tolerances. Nonlinear solve batches grouped by network/group count or stiffness only if sorting cost and state mapping are bounded. Do not materialize all-cell active lists when views suffice. BH operations parallelize over sparse active BHs and patch queries; use existing spatial indexes, not O(N_BH*N_cells) global scans.

## 5. Spatial MPI and load balance

Persistent patch and particle ownership are the authority. RT halo exchanges use existing patch topology, bounded peer buffers, explicit stage and epoch IDs. Prefer packed SoA/contiguous slabs, receive directly into ghosts where safe, no multiple wire copies. MPI count limits use bounded rounds or verified large-count APIs; progress cannot depend on empty ranks participating incorrectly. Interior stencil work overlaps exchange if measured beneficial; boundary fluxes and reflux are deterministic or tolerance-equivalent. Regrid migrates canonical radiation and associated scheduler history once. Physical source events use owner-target protocol and aggregate sufficient statistics, not root-global catalogs or all-rank gas arrays. Distributed halo seeding provider must avoid unbounded global truth replication.

Work weighting includes cell transport cost G, RT subcycle count, nonlinear iterations/retries, hydro/MHD, BH query/deposition load, angular/packet work and memory pressure. Rebalance only when predicted gain outweighs migration; record imbalance and rank max. Pin active patch and source event to topology epoch; no regrid mid-commit. Cross-rank local subcycling is opt-in only after temporal ghosts and flux-register synchronization pass. Level-synchronous RT is safe fallback when scientifically equivalent at a smaller timestep. Global timestep reductions are explicit measured costs. Avoid collectives per BH when neighborhood reductions suffice.

## 6. Accelerators and large HPC

P0 reference CPU/OpenMP+MPI; P1 GPU production transport and local coupled solves with same physics equations, stage ordering and versioned datasets. Keep canonical RT patch data device-resident over subcycles where feasible; exchange only required halos and changed coupling state. No host roundtrip of complete radiation state per stage. GPU kernel boundaries preserve stage synchronization; fuse local operations if memory traffic improves and numerical equivalence survives. Batch stiffness groups to mitigate warp divergence; bounded fallback for exceptional cells must not change chemistry or transport model. Use GPU-aware MPI or explicit pinned staging, measured independently. Device budget includes library workspaces, pinned host staging and asynchronous buffers. Multi-GPU placement follows actual NUMA/interconnect topology. CPU reference always available, no separate unqualified GPU physics.

Supercomputer execution must support distributed spatial decomposition without rank-count assumptions. Use 64-bit global IDs and sizes; local index narrowing only with range proof. As scale increases, optional frequency and angular decomposition become promising for directional/VET/Monte Carlo backends, but require matter-coupling reductions and cross-axis memory accounting. M1 primary spatial decomposition is simpler and preferred until benchmark evidence supports extra axes. P1 VET needs bounded nonlocal closure/sweep workspace, explicit communication schedule and convergence counters. Ray/packet queues use backpressure, bounded rounds and migration accounting; no unbounded per-source fanout. Do not claim GPU or massive-MPI speedups until measured. Backend registration exposes numerical capabilities, precision, memory contract and required communication rather than hardcoded machine classes.

## 7. Diagnostics and performance metrics

Per RT step: leaf cell/group updates, ghost bytes, face flux bytes, IMEX iterations and failures, subcycles, c_hat, transport/interaction/reflux wall time, achieved bandwidth and throughput. Per BH: active count, neighborhood queries/visited patches, candidates, MPI traffic, supply/actual withdrawal, disk integration steps and event delays. Per AGN: recipient counts, plans/commits/retries, queue capacity, refinement requested/admitted/denied, high-water by owner and regrid time. Per run: rank max/mean wall time, node memory sum, per-device peak, checkpoint size/time, calibration and backend IDs. Performance diagnostics should be compact counters, not per-element strings; detailed samples optional and budgeted.

Evaluate strong scaling at fixed problem and weak scaling at fixed workload per rank/device, but also time-to-accuracy and memory-to-accuracy for same physical target. Report compiler flags, commit, machine, rank/thread affinity, GPU model, spectral/chemistry profile, AMR levels, c_hat, tolerances and iteration counts. No comparison of different accuracy/physics profiles as speedup. Include cold-start, steady-state and regrid/checkpoint costs. Define baseline before optimization.

## 8. Proposed benchmark ladder and acceptance

| Gate | Proposed experiment |
|---|---|
| Disabled cost | identical DMO/hydro/MHD output within existing deterministic contract; zero RT arrays/BH optional sidecars when disabled |
| Memory slope | 64³,128³,256³ and safe 512³ probes with 1/3/6/9 groups; compare measured canonical capacity to scale law; report covered parents and ghosts |
| Stage scratch | repeated RT steps, no monotonic growth; bounded batches/leases; peak matches reservation within declared allocator reserve |
| Source sparsity | fixed cells with increasing BH count, fixed BH with increasing cells; demonstrate indexed local queries, bounded communication |
| AMR pressure | jet refinement and regrid under tight budget; verify no partial mutation or silent model change |
| Rank equivalence | serial and 2/4/8 MPI plus uneven/empty ranks; same conserved ledgers and numerical tolerance, restart rank changes when supported |
| CPU/GPU | same physics and input, compare convergence, residual distributions, wall time and peak memory |
| Single workstation | reference problem fits configured safe budget, no swap dependence, record measured speed and memory |
| Small/medium HPC | strong/weak scaling series with representative radiation and clustered source loads |
| Supercomputer | staged 64→256→1024+ ranks/devices as access permits, no root-global N state, communication scaling and imbalance reported |

Numeric performance targets shall be set after instrumented baselines; do not invent speedup claims. Proposed mandatory structural targets: no persistent full-population duplicated gas representation, no unbounded BH×cell or source×ray exchange staging, disabled-module memory independent of gas population, all large allocations accounted and admitted, no untracked physical fallback. Scientific residuals and convergence gates in master and companion documents must pass before claiming any optimization. Hardware profiles express workload size and resource budget, not exclusive access to physical features.

## 9. Integration preflight and risks

Inspect actual `MemoryGovernor`, phase arena, MPI protocol, `SimulationState`, device backend, and source runtime contracts in the pinned development commit. Extend existing governor rather than add a competing allocator. Reconcile existing MHD buffers and dual-energy paths; avoid double staging chemistry. Candidate new benchmark files follow live `.cpp` naming, e.g. `bench/bench_rt_transport_throughput.cpp`, `bench/bench_rt_coupling_batch.cpp`, `bench/bench_agn_sparse_deposition.cpp`. Use existing bench/test registration. Major risks: radiation CFL dominates tiny fine patches; IMEX stage storage; stiff chemistry divergence; nonlocal closure scaling; hot BH neighborhoods; overlapping AMR regrid buffers; GPU data movement; concurrent campaigns changing schemas. Mitigate with measured subcycling, batching, admission, topology epochs, explicit interface amendments and separate qualification. No promise that uniform 512³ full-physics RT fits 48 GiB: 9-group moments alone are 36 GiB before any other state.

## References

- CHUÍ Memory Architecture and Resource Governance Specification v0.1, 2026-08-23, particularly sections 2, 4–12, 14–17.
- CHUÍ revised Master Design Document, 2026-03-11, workstation/small-cluster and scale-up objectives.
- He et al. (2024), QUOKKA multigroup GPU radiation hydrodynamics: https://arxiv.org/abs/2407.18304
- He et al. (2024), local IMEX: https://arxiv.org/abs/2404.08247

All costs not explicitly derived from field counts are prospective engineering requirements, not measured CHUÍ performance.
