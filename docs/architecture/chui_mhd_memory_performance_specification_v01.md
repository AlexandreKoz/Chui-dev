---
title: "CHUÍ | MHD Memory and Performance Architecture"
subtitle: "Normative Resource, Execution, MPI/GPU, Profiling and Pareto Benchmark Specification — CHUI-MHD-PERF v1.0"
author: "CHUÍ / CosmoSim technical design"
date: "22 November 2026"
lang: en-US
---

# Document control

| Item | Fixed contract |
|---|---|
| Identifier | `chui_mhd_memory_performance_specification_v01` |
| Version and status | v1.0; **proposed normative design**, not an implemented or benchmark-certified subsystem |
| Parent physics design | `chui_mhd_architecture_specification_v01` (CHUI-MHD-REF v1.0; 22 November 2026) |
| Parent resource policy | `chui_memory_architecture_and_resource_governance_specification_v01` (23 August 2026) |
| Runtime baseline | `AlexandreKoz/CHUI-Cosmology-Simulator`, `main` inspected 22–23 November 2026 |
| Scope | Ideal and non-ideal MHD, patch-AMR, cosmological/isolated modes, CPU/OpenMP, MPI, optional GPU |
| Primary hardware | 48 GiB workstation; reference small clusters; constrained 8 GiB-class GPU backend |
| Scope boundary | Does not change MHD equations, precision, reference solvers, cosmological conventions, or scientific acceptance criteria |

**Precedence.** The MHD physics specification controls equations, numerical methods, face/edge authority, SSPRK2, UCT-HLL, Hall-HLL, RKL2, AMR synchronisation, reference profiles and physical failure semantics. The overarching CHUÍ memory-governance specification controls global budgeting, reservation authority, hard ceilings, and allocation classes. This companion specifies how the **MHD subsystem implements and measures** those contracts. Where an existing run has a stricter externally imposed machine limit, the stricter limit wins. Changes to any frozen scientific method require an amendment to the parent physics document, never an unrecorded “performance” patch.

**Normative language.** MUST/SHALL means a release requirement; SHOULD indicates a default for which a documented and tested exception is possible; MAY denotes an optional mode. Numerical budgets in this document are *engineering acceptance targets* rather than claims of observed performance. A method can be numerically certified but **not** performance-qualified, and vice versa; neither status may be represented as the other.

# 1. Mission and measurable objective

CHUÍ's aim is to maximize **scientifically validated fidelity**, predictive value, and freedom to test alternative physics under finite hardware resources. MHD performance engineering is therefore not a contest for the most simulated cell-updates per second; it is a means of making credible simulations possible on workstations and small clusters while retaining a scalable HPC path.

For a fixed physical problem $p$, an accuracy specification $\varepsilon_\star$, and hardware $h$, define

$$
T_{\mathrm{accuracy}}(p,\varepsilon_\star,h)
=\min_{\theta\in \Theta_{\mathrm{valid}}} T_{\mathrm{wall}}(p,\theta,h),
$$

subject to $\varepsilon_{\mathrm{numeric}}\le \varepsilon_\star$, accepted physical and conservation criteria, $M_{\mathrm{peak}}\le M_{\mathrm{safe}}$, valid provenance and restart behavior, and no hidden change of physical equations. $\theta$ includes batch size, patch scheduling, placement, permitted numerically equivalent kernel variants and launch configuration; reference physics choices such as HLLD/HLLE, UCT-HLL, floating-point authority, and tolerances are fixed in their certified profiles.

The result for each candidate is a **vector**, not a scalar prestige score:

$$
\mathcal P=(\boldsymbol\varepsilon_{\mathrm{science}},\ T_{\mathrm{wall}},\ M_{\mathrm{RSS,max}},\ M_{\mathrm{device,max}},\ V_{\mathrm{MPI}},\ E_{\mathrm{joule}},\ N_{\mathrm{failures}}).
$$

Energy in joules is included only when measured using a documented instrument or counter; otherwise it must be `unavailable` rather than estimated from TDP. Pareto comparisons are meaningful only within a fixed problem, matched physics, accuracy protocol and hardware class.

# 2. What is already real in CHUÍ

This design attaches to the **existing** `cosmosim::...` source tree, not a hypothetical fresh code. At inspection, the following interfaces were verified:

- `core::MemoryGovernor`, `MemoryReservation`, `DeterministicBatchSizingPolicy`, `selectDeterministicBatchSizeForWorkspaces`, `MemoryGovernorSnapshot` and `core::GovernedScratchArena` already implement admission and explicitly bounded PMR arenas. In particular, an arena uses a null upstream allocator; underestimated scratch **fails** rather than escaping to the heap.
- The production governor has `kCanonicalPersistent`, `kPersistentCache`, `kPhaseResident`, `kScratchArena`, `kCommunication`, `kDiagnostic`, and `kExternalRuntime`, along with `kGreen/kAmber/kRed/kTrip` pressure states. `hard_limit_bytes=0` currently means unlimited accounting; MHD certified resource profiles may **not** rely on this unlimited behavior.
- `HydroConservedState` is a six-double cell-centered record (mass, three momentum components, total fluid energy and metal density); MHD changes energy semantics through the parent specification but must not silently allocate a second independent population-scale gas state.
- `amr::PatchHierarchy`, `FluxRegisterAccumulator`, `RefluxSynchronizer`, local AMR subcycling, and stable gas-cell IDs exist. Magnetic edge-reflux authority is new. MPI-distributed AMR subcycling and remote temporal interpolation are **not** treated as existing production capabilities; synchronized MPI AMR is the reference starting point.
- `parallel::DecompositionWorkComponents` already includes AMR cost, active fraction, memory pressure, transient memory, communication and accelerator occupancy cost components; MHD must augment those measurements rather than create a competing global decomposition system.
- The accepted simulation configuration is typed **`.param.txt`**, not YAML. The parent specification reserves `[mhd]`, `[mhd_nonideal]`, `[mhd_amr]`, `[mhd_io]` and versioned science/checkpoint schemas.

None of the above is evidence that CHUÍ currently contains working MHD kernels. This document describes additions necessary for that implementation.

# 3. Operating envelopes, ceilings and priority rules

## 3.1 Fixed memory profile

For the reference workstation, retain the August CHUÍ memory-governance policy:

| Resource | Proposed default | Interpretation |
|---|---:|---|
| Host physical memory | 48 GiB | Illustrative reference machine, **not** allocatable capacity |
| Host OS/library/runtime reserve | 6 GiB | Must include measured MPI/FFTW/HDF5/allocator needs |
| Safe simulation process ceiling | **42 GiB** | Hard admission threshold across all simultaneously resident ranks/processes on the host |
| Nominal sustained residency | **36 GiB** | Preferred operating region; non-fatal temporary peaks only when admitted |
| Contingency band | 36–42 GiB | Short, justified peaks; growth forbidden without admission |
| Per-rank MHD communication budget | **64 MiB** | Default combined send, receive and reconciliation **working buffers**; not per peer |
| Per-rank MHD communication max | **128 MiB** | Explicit high-memory option, subject to whole-node governor approval |
| MHD temporary active-cell budget | **512 B per active cell** | Initial PPM reference engineering cap for MHD-owned phase resident + scratch, excluding existing canonical state |
| Optional 8 GiB-class GPU | **6 GiB** governed live allocation target | Initial working ceiling; reserve at least 2 GiB for display/runtime/other use unless device discovery proves a stricter ceiling |

**Important:** 42 GiB is a **whole-host** MHD+hydro+gravity+communication+external-runtime limit, not 42 GiB per MPI rank. Across multiple ranks on one host, the sum of rank RSS and explicitly shared buffers must be reconciled without incorrectly double-counting shared mappings. Governor accounting is rank-local; a new node-level admission coordinator or conservative aggregate preflight is therefore mandatory for multi-rank workstation certification. The 64/128 MiB communication targets are the MHD subsystem's buffers, independent from the older **128 MiB/rank PM-deposition objective**, and cannot be assumed free when PM and MHD overlap.

For other machines, set $M_{\mathrm{safe}}=M_{\mathrm{physical}}-M_{\mathrm{OS/runtime\ reserve}}$ with hardware-dependent reserve and measured external-runtime costs. **Do not charge the same reserve twice:** the 6 GiB workstation reserve lies outside the 42 GiB simulation-process ceiling; `external_runtime_reserve_bytes` in `MemoryGovernorPolicy` must describe opaque library/runtime demand *inside* the admitted process accounting, not be a second copy of the entire OS reserve. `MemoryGovernorPolicy::hard_limit_bytes=0` is supported for legacy behavior but disallowed for MHD performance-qualified runs.

For GPU operation, device and host budgets are **separate**. Pinned host buffers remain counted in the host budget; device allocations must be admitted to a device governor before allocation. GPU driver memory, CUDA contexts, allocator pools and other external reservations require measured headroom. No out-of-core oversubscription of live MHD state is accepted as a reference fallback.

## 3.2 Pressure and degradation

At green pressure, admit performance-positive caches and multiple independent batches. At amber pressure (existing governor threshold: 85% of admitted budget), pause optional analysis and speculative caches. At red pressure (95%), admit only mandatory work already proved to fit, reduce future batches within certified bounds, and avoid overlapping HDF5 output with an MHD peak. At trip, reject before allocation, identify the owner and request, and preserve the last accepted state. If a physical timestep is too stiff, fail with the parent specification's stiffness error; **never** hide stiffness by lowering physical coefficients or magnetic precision.

Pressure state is not a scientific control variable: it can change *scheduling, concurrency, scratch sizes and output timing*, but never the physical equations, spatial order of the reference profile, or accepted conservation tolerances.

# 4. Exact geometric memory equations

Let an active Cartesian patch contain $n_x n_y n_z$ real cells, $N_c=n_xn_yn_z$. The face-centered *persistent* double-precision flux authority is:

$$
N_f=(n_x+1)n_yn_z+n_x(n_y+1)n_z+n_xn_y(n_z+1),\qquad M_{\Phi}=8N_f\ \mathrm{bytes}.
$$

The stage-local edge EMF arrays contain:

$$
N_e=n_x(n_y+1)(n_z+1)+(n_x+1)n_y(n_z+1)+(n_x+1)(n_y+1)n_z,\qquad M_{\mathcal E}=8N_e\ \mathrm{bytes}.
$$

For cubic $n$ patches: $N_f=3(n+1)n^2$, $N_e=3n(n+1)^2$. These counts describe *full contiguous patch arrays*, without halo or MPI wire duplication. Magnetic flux has dimensions of field times **comoving** face area; the EMF holds an oriented **line integral** and is accumulated in the same physical stage/scale-factor convention as the parent MHD specification.

| Cubic patch | Real cells | Face flux f64 | Edge EMF f64 | Full-cell ghost expansion with PPM $g=3$ |
|---|---:|---:|---:|---:|
| $16^3$ | 4,096 | 0.100 MiB | 0.106 MiB | 2.60x |
| $32^3$ | 32,768 | **0.773 MiB** | **0.798 MiB** | 1.675x |
| $64^3$ | 262,144 | 6.094 MiB | 6.189 MiB | 1.308x |

Ghost expansion factor for a *complete cell-centered ghost cube* is $[(n+2g)/n]^3$, not the mandatory cost of every algorithm. In particular, do not allocate that full cube when only boundary tiles are needed. PLM with $g=2$ gives corresponding factors 1.953x, 1.424x and 1.199x for $n=16,32,64$.

A uniform global $512^3$ magnetic grid requires approximately **3.006 GiB** of face flux in double precision if its three axes are each stored as one global array (global face-boundary overhead included). If represented as 4,096 distinct $32^3$ patches with mirrored patch faces, the same magnetic storage is **3.094 GiB** before any ghosts, parent-level state and allocator overhead. This ~3% overhead is a *calculated geometric bound*, not a benchmark.

At $512^3$, **one** f64 lane spanning all cells costs exactly 1 GiB. No persistent full-gas-population vectors of `B_cell`, current, magnetic gradients, three transport coefficients, or per-face wave states are authorized in the reference. A persistent global derived lane requires a separately versioned memory justification and demonstrated time-to-accuracy advantage.

# 5. Canonical, phase-local and scratch ownership

| Allocation | Authority and lifetime | Memory class | Policy |
|---|---|---|---|
| Oriented `face_flux_*_comov` | Persistent per active patch; checkpoint truth | `kCanonicalPersistent` | f64 only; no second canonical magnetic representation |
| Patch geometric key, owner/epoch, shape/offset | Persistent structural truth | `kCanonicalPersistent` | Stable identity; compact integers, no strings per cell/face |
| Cell mass/momentum/total energy/metals | Existing authoritative hydro/gas state | Existing owner | No duplicate MHD gas container at rest |
| Cell-centered $B$ and primitive state | Recomputed for active tile/stage | `kPhaseResident` | No full-population permanent cache |
| HLLD left/right reconstructed states and wave speeds | Current face tile | `kScratchArena` | Stream and recycle; no persistent per-face objects |
| Edge EMF and discrete current | Current stage/patch | `kScratchArena` | Unique edge ownership, no per-face duplicate EMF list |
| SSPRK2 trial states and rollback snapshot | Current active interval | `kPhaseResident` / `kScratchArena` | Atomic accept/reject; do not modify persistent truth before commit |
| Ohmic/Hall/ambipolar coefficients | Stage-local derived from provider | `kPhaseResident` | Avoid three full-grid coefficient lanes |
| RKL2 intermediate magnetic state | One diffusion cycle | `kScratchArena` | Recycle stage buffers; stage count changes work, not persistent residency |
| AMR edge-integral registers | Until coarse/fine interval accepted | `kCanonicalPersistent` while pending | Sparse boundary-scaled and epoch keyed; drained before safe checkpoint |
| MPI face/edge halo and migration packets | Bounded rounds | `kCommunication` | Count send and receive simultaneously; fail on stale epochs |
| Device buffers and pinned host staging | Explicit device lease + host lease | Device/host budgets | No silent host/device double authority |
| Diagnostics, spectra, visualization | Episodic | `kDiagnostic` | Never causes mandatory physics to exceed hard ceiling |

**Authority invariant:** one logical owner exists for each magnetic face and time-integrated edge; shared physical boundaries may have mirrored read-only copies, but only the designated owner commits. A view may alias this storage; an out-of-date copy does not become a second authority. All retained `std::vector`/aligned-vector **capacities**, not just logical sizes, count toward ownership.

**Phase lifetime:** prepare > admit leases > read canonical input > construct one active patch/tile view > reconstruct > Riemann / EMF > update trial gas + faces > reconcile MPI/AMR > validate energy/divergence/positivity > commit atomically > release/reuse scratch. Every pointer acquired from the phase arena becomes invalid after the arena's release. `GovernedScratchArena` and null upstream resource are the default allocation mechanism; an unexpected `std::bad_alloc`/overflow is a correctness failure and **cannot** be misinterpreted as optional-task deferral.

# 6. Patch sizes, tiles, and batching

The **reference performance probe** uses $32^3$ Cartesian patches with PPM's $g=3$ and 64-byte host alignment. This is a benchmark geometry, **not** permission to change scientific AMR refinement geometry to fit a kernel. Existing patch dimensions remain authoritative. Allow $16^3$ and $64^3$ as prescribed benchmarking alternatives; reject unsupported dimensions only if they violate an actual MHD stencil or topology contract, not because performance is poor.

Default MHD-owned `PhaseResident + ScratchArena` active capacity is **at most 512 bytes per active real cell**, including: two-stage RK trial state/rollback, cell primitive cache, face reconstruction buffer, edge/current work, non-ideal stage data, vector growth headroom, padding, and temporary index maps. Exclude previously existing persistent gas state, authoritative persistent face flux, separately governed MPI communication and optional diagnostics. For a $32^3$ batch this is **16 MiB**; eight concurrently active independent batches require a maximum of **128 MiB** of such scratch plus their other classes. This limit is deliberately conservative and should be benchmarked; any increase is a resource-contract amendment, not a silent allocator expansion.

Stage memory has a particularly important restriction. **SSPRK2 and rejection need explicit old/stage/new live-set accounting.** Each trial must keep an immutable recoverable image of the currently active canonical gas-and-face data, or use copy-on-write/in-place algorithmic proof that delivers the same atomicity. At the basic six-double hydro state width, one $32^3$ gas trial image is 1.5 MiB and one face trial image is 0.773 MiB. Two independently held gas-and-face trial images alone therefore occupy about **4.55 MiB**; this is already part of the 16 MiB batch allowance, not an additional omitted cost. New optional chemical scalar lanes increase the actual trial width and must reduce permitted batch size or justify a separate reviewed contract.

**Tile policy:** boundary ghost staging is tile/face-shell-based by default. The producer may construct a complete $(n+2g)^3$ temporary cube only when its admitted retained capacity fits within the same 512 B/cell cap and its time-to-accuracy is demonstrably better. Edge EMF arrays may be kept full-patch if they avoid expensive recomputation and fit within that cap. Stream per-face HLLD results directly into gas-flux accumulators plus the adjacent edge-EMF builder; there is no persistent all-face primitive-state or flux object graph.

Batch selection uses the existing deterministic workspace-sizing API. Begin at the maximum of **eight $32^3$ equivalent patches per task group** (262,144 active cells) and shrink to the largest capacity that the available headroom admits, rounded down to full patch or valid tile units. The minimum admitted unit is a whole patch unless the patch's measured scratch would exceed 16 MiB; in that case, explicit tile splitting is required. Zero-work MPI ranks participate in synchronization without allocating synthetic scratch arrays. Host threads use task-local arenas or disjoint slices of a single governed arena; multiply per-worker scratch in preflight, never count it once and allocate it once per thread.

The project memory governor may retain performance-positive caches **only** when their amortized reuse beats recomputation within available headroom. Each cache has `owner`, capacity, last-use epoch, invalidation events, high-water and measured wall-time benefit. Never use stale cached geometry or coefficients after regridding, migration or source-state changes.

# 7. CPU execution and data-movement architecture

## 7.1 Reference kernel sequence

The reference spatial solver is fixed by CHUI-MHD-REF v1: second-order unsplit finite volume, SSPRK2, PPM (or the separately certified PLM minimal mode), HLLD with HLLE fallback, UCT-HLL and face-flux CT. The execution backend may reorder **mathematically equivalent** patch-local work but must not alter numerical stage boundaries, face fluxes, edge orientation or accepted reduction semantics beyond verified floating-point tolerances.

For each active patch and SSPRK2 stage, execute the following dependency graph:

```text
admit/prepare patch + stage-time halos
    |
    +--> recover gas primitives and B_cell from authoritative face flux
    |       |
    |       +--> PPM reconstruction (small contiguous tiles)
    |                 |
    |                 +--> HLLD/HLLE flux and face electric signals
    |                            |
    |                            +--> cell-conserved updates (unique owners)
    |                            +--> UCT-HLL edge EMFs (unique edge owners)
    |                                        |
    +---------------------------------------+--> CT face updates
                                                 |
                                    stage validation + registered fluxes
                                                 |
                                    MPI/AMR agreement + atomic commit
```

The non-ideal stages reuse the magnetic mesh/face authority but compute stage-local discrete currents, coefficients, Poynting fluxes and Hall/diffusive EMFs. Each Strang component owns its own measured interval; stage work cannot be silently deferred until after a supposedly accepted physics step. The CPU implementation MUST avoid per-face heap allocations and per-face virtual calls, dynamic strings, hash map lookups or MPI requests in the hot Riemann/EMF loops. Per-patch runtime selection outside the inner loop is permitted; the general `HydroRiemannSolver` virtual interface is not a mandate to replicate virtual dispatch for MHD.

## 7.2 Parallel ownership and scheduling

- Parallelize across active patches first; when insufficient patches exist, parallelize over disjoint face/edge tiles with explicit ownership. At shared edges, compute or reduce in **one** deterministic owner task; no write races or unbounded atomic contention.
- Stage barriers occur at actual mathematical dependencies (e.g. after complete face signals before edge EMFs, after owner edge reconciliation before CT), not after each face or component.
- At green memory pressure, use task parallelism among independent patches. At amber/red pressure, shrink simultaneous patch count without changing equations or per-cell reference order.
- First-touch placement should follow the task's NUMA locality. Avoid cross-socket cache bouncing of face/edge arrays; use measured affinity profiles, not hard-coded P-core/E-core or socket assumptions.
- Align vectorizable SoA arrays to at least 64 bytes. Use continuous axis-indexed loops, precomputed stride metadata, and branch-minimized fallbacks. HLLD fallback faces are counted and handled correctly; a branch-free approximation is not acceptable if it changes the numerical flux.
- Build microbench variants for **tile-local recompute** versus **short-lived cache of primitives**. Choose only on measured time-to-accuracy; the fastest raw HLLD loop need not win the full stage.

**Floating-point policy.** All canonical conserved magnetic fluxes, gas total-energy state and non-ideal diffusion coefficients retain the parent reference precision (f64). Fused multiply-add, vector reductions or fast-math flags are allowed in research/performance trials only if the reference scientific gates and deterministic-mode requirements hold. No implicit f32 demotion of magnetic authority; a precision campaign is a different versioned scientific study. Reproducibility-mode arithmetic flags are recorded.

## 7.3 Data movement and roofline

The performance profile MUST report (a) algorithmic **minimum** input/output bytes computed from kernel storage fields, (b) observed hardware memory bytes where available, (c) stage wall time, (d) estimated and measured arithmetic intensity with caveats, and (e) patch size and cache residency. Hardware counters are platform dependent; unsupported measures become `not_available` rather than assumed numbers.

**No fabricated FLOP or bandwidth constants.** A roofline upper bound is $P_{\mathrm{attainable}}\le \min(P_{\mathrm{peak}},I\,B_{\mathrm{measured}})$, where $I$ is measured/estimated FLOPs per byte and $B_{\mathrm{measured}}$ is a benchmarked sustained memory-bandwidth bound. Streamable stages (primitive conversion, flux accumulation, CT) and arithmetic-heavy HLLD or nonlinear coefficient evaluation must be profiled separately. Avoid presenting an aggregate cells/sec number as proof of FLOP efficiency.

# 8. MPI ownership, message bounds and communication overlap

## 8.1 Message semantics

Each logical shared face and edge has exactly one deterministic owner identified by the parent specification's stable geometric key and geometry epoch. The domain decomposition system owns rank placement. Magnetic flux state, time-integrated edge EMFs, coarse-fine edge-register contributions, and regrid/migration metadata are exchanged in **separate typed bounded messages** (or compact wire sections), carrying only non-derivable fields. Do not send duplicate complete patch copies just to obtain one face shell.

Messages must contain schema/epoch, source/target rank, patch/face or edge key range, level, interval, numeric datatype and payload count; a batch-level header can amortize shared metadata. The wire format is distinct from in-memory SoA, but permitted transformation must be at most **one bounded compact packing buffer** per direction. The receiver must consume into the authoritative or transient destination without an additional full decoded population.

Default **aggregate** MHD communication staging per MPI rank is 64 MiB (including all simultaneous send, receive, temporary reductions and unpack buffers); the explicit maximum is 128 MiB. Each peer's share is allocated from this one pool, and dynamic credit-based chunking prevents total simultaneous allocation from exceeding it. The pool must reserve send and receive coexistence, so a naive implementation cannot claim two 64 MiB buffers as one. If one patch message exceeds the budget, split by planes/edge key ranges; never bypass the limit. A chunk must respect MPI count and datatype constraints even for very large meshes.

Communication growth model for a regular partition is $O(g A_{\mathrm{boundary}})$ for ghost state plus $O(E_{\mathrm{shared}})$ for shared EMFs and AMR registers, not $O(N_{\mathrm{local}}\times\mathrm{stencil})$. In a highly fragmented AMR partition, boundary measure can itself approach population scale; the cost must be measured honestly and inform partitioning. **“Surface scaled” is a design preference with geometric conditions, not a claim for arbitrary AMR partitions.**

## 8.2 Nonblocking exchange without excess residency

Communications may be nonblocking only when the owner and epoch are frozen, the destination buffers are admitted, and progress does not race with stage mutation. Overlap reconstruction of interiors while boundary messages are in flight; after exchange completion perform deterministic owner reconciliation for shared EMFs, then CT. No full additional ghost copy is allowed merely to enable overlap.

The critical path is the maximum rank time to stage completion, not the mean. Report: send/receive bytes, message count, min/median/p95 message sizes, wait time, packing time, unpacking time, edge reconciliation, AMR register time, rank-idle imbalance and memory high-water. Rank-global reductions must include zero-cell ranks without invalid collectives. MPI failure or epoch disagreement fails collectively and leaves the previous accepted global state intact.

## 8.3 Load balancing and migration

Augment existing `DecompositionWorkComponents` with **measured** MHD weights: active gas cells, active face/edge counts, reconstruction order, fallback count, Hall subcycle count, RKL stage count, boundary work, recent cache-miss or GPU occupancy proxy and memory high-water. Use a smoothed weighted estimate from previous intervals (avoiding oscillation on one transient outlier). Weights are *learned from profiling* and stored in provenance; do not invent a fixed “one Hall cell = X ideal cells” constant.

Migration uses a two-phase admission/commit transaction: reserve sender retained old state + receiver new state + both bounded wire buffers + any active AMR register/epoch overhead; validate all destination capacity before any authoritative ownership change; transfer in bounded chunks; collectively validate shape, flux totals and epoch; then publish new ownership and free old capacity. When safe coexistence cannot fit, postpone migration to the next sync or fail before modifying state. Rank-count-changing restarts remain outside the parent v1 contract.

# 9. AMR, regrid, conservation and topology memory peaks

AMR operations require more than storing the leaf field once. Before a regrid transaction, the admission estimate is

$$
M_{\mathrm{regrid,peak}}=M_{\mathrm{old\ state}}+M_{\mathrm{new\ state}}+M_{\mathrm{transfer\ scratch}}+M_{\mathrm{boundary\ registers}}+M_{\mathrm{wire}}+M_{\mathrm{external\ reserve}}.
$$

The `old` and `new` terms describe **only simultaneously coexisting allocations**; existing canonical memory already charged to the governor must not be double counted, while a reallocating vector must reserve its **replacement alongside the retained old capacity**. A candidate transaction is rejected if any rank or whole-host admission fails. Gas and magnetic topology are published together; stale flux keys or partial face interpolation cause rollback.

**Sparse registers.** Ordinary fluid flux registers and magnetic EMF reflux-curl registers may grow only with coarse/fine interfaces and integrated stage count, **not** in proportion to every active face or every Hall/RKL substep held separately. During a substep, accumulate immediately into the one keyed, time-integrated edge authority; do not materialize a `N_substeps × N_edges` history. For each boundary register record, account for key, integer interval, integrated EMF, owner and geometry epoch, coverage/completeness metadata, vector capacity and allocator overhead. A bound must be computable from the number of legal interface edges **before** stepping.

**Temporal histories and subcycling.** The initial multi-rank reference uses synchronized AMR level sweeps. The future `mhd_ref_amr_subcycled` certification tier must budget coarse temporal magnetic-face state, conservative gas history and pending EMF registers along with remote temporal ghosts; those objects cannot be hidden inside ordinary hydro history counts. Any incomplete interval is forbidden at an MHD checkpoint boundary under the parent physics specification. During regrid, reject if pending reflux-curl or external interval histories remain active.

**Cancellation/rollback.** After an invalid Riemann state, divergence check, stale MPI geometry or failed memory admission, every trial patch returns to the last accepted canonical gas+magnetic state. Rollback snapshots are scoped to active batches; never copy the entire global gas/magnetic population merely to undo one local stage. A failure following partial cross-rank acceptance is a collective consistency error, not a reason to continue on surviving ranks.

# 10. Non-ideal operator costs and stiffness control

## 10.1 Operator work and residency

The parent method is `D(dt/2) -> H(dt/2) -> I(dt) -> H(dt/2) -> D(dt/2)`, with coupled total-energy and EMF accounting. Diffusion $D$ uses explicit SSPRK2 below the defined stiffness ratio and RKL2 otherwise, with a **128-stage cap per RKL cycle**. Hall $H$ uses the reference explicit HLL/SSPRK2 and bounded subcycling; its actual stable timestep comes from the parent physical and whistler-wave criteria.

- **Diffusion:** stage scratch is reusable independent of the number of RKL stages. Recompute required field-dependent coefficient bounds/coefficients each stage; report coefficient provider time and cache hit rate, not merely diffusion update time.
- **Hall:** two SSPRK2 stages are performed for each accepted Hall subcycle; update the same bounded magnetic/energy trial buffers. Substep count increases execution time but **cannot** multiply resident buffer count.
- **Chemistry/table:** cache a versioned, read-only HDF5 transport table at rank or node scope where safe. Its measured storage must be admitted as `kPersistentCache` or read-only external runtime. The hot loop queries a small batch of cells with consistent time and physical units; no per-cell heap allocation or table object construction is allowed.
- **Mixed composition:** physical coefficient providers may request auxiliary charge-carrier state; their persistent chemistry lanes are owned by the chemistry subsystem and counted once. MHD cannot duplicate a complete chemical state in a private cache.

## 10.2 Stiffness forecasting and cost refusal

For an attempted physical interval $\Delta t$ the scheduler computes conservative candidate substeps for ideal, diffusion and Hall processes using the parent equations. Before committing to the interval, preflight an upper-bound forecast of $N_H$, $N_{\mathrm{RKL\ cycles}}$, and $N_{\mathrm{RKL\ stages}}$ with observed active-cell distribution; periodically revise the forecast when coefficients change.

**Reference operational policy:** by default, permit up to **4,096 Hall SSPRK2 subcycles per accepted *level-local* MHD interval**. When the forecast exceeds this cap, reject the candidate interval and ask the existing scheduler to shorten it, preserving the parent stability constraints. If repeated shortening cannot satisfy both physics and global scheduling rules within the parent retry policy, return `MHD_STIFFNESS_BUDGET_EXCEEDED` with the predicted work and the dominant cells/levels. The cap is an *operational cost ceiling*, not a physics coefficient clamp and not a reason to truncate a physical interval. This new cost guard is subordinate to the parent eight-retry policy and must be included in the normalized `.param.txt` configuration. Deliberately increasing the cap is an explicit performance-profile change, not a silent dynamic override.

RKL cycles obey the parent 128-stage cap; a longer diffusion interval is partitioned into multiple complete cycles rather than dropping stages. An equivalent preflight operational ceiling is **4,096 total RKL stages per level-local interval**; exceeding it triggers scheduler shortening or a reported stiffness refusal. It is not a modification to the 128-stage-per-cycle physical integrator. Both caps may be raised only with a versioned memory/performance-profile override and measured resource justification; the certified physics is unchanged.

Keep separate counters for **attempted**, **accepted** and **rejected** stages so that a repeated failed interval does not masquerade as numerical throughput. The total cost model is

$$
W_{\mathrm{MHD}}\approx 2N_c C_I + 2N_H N_c C_H+\sum_{j=1}^{N_D} N_{\mathrm{RKL},j}N_c C_D
+ W_{\mathrm{halo}}+W_{\mathrm{AMR}}+W_{\mathrm{sources}},
$$

where $C_I,C_H,C_D$ are **measured** per-cell stage costs and $N_H$ is the number of Hall subcycles per *half-interval* (the factor of two accounts for the two Hall operators; each $C_H$ includes its SSPRK2 stages). The formula is a performance estimator, not a scientific time integration rule, and can be refined per level, cell count and communication regime. Record the estimator residual over completed runs.

# 11. GPU and accelerator backend

The CPU reference is the scientific authority. The GPU implementation reproduces **the same** HLLD/HLLE, reconstruction, UCT-HLL, SSPRK2, Hall and RKL algorithms, discrete operators and total-energy convention. A separate GPU-specific physical solver is prohibited in the reference. Device kernel fusion is permitted only when it preserves edge-owner and stage synchronization semantics and passes the parent numerical equivalence gates.

**Data placement:** for active GPU-resident patch batches, use SoA device fields and one versioned ownership/dirty-state tracker. Copy persistent magnetic state and necessary gas lanes to the device once per resident interval; avoid copying entire patch populations each RK/STS/Hall substep. Reuse device-resident EMFs and stage buffers within the governed active batch. Halo exchange may use staged pinned host buffers unless GPUDirect or GPU-aware MPI is verified; both host and device staging are governed. Never assume GPU-aware MPI from a compiler flag alone.

**Initial 8 GiB-class budget:** target at most **6 GiB** MHD-related device live set, subject to current free memory, driver/runtime reserve and memory accounting. Because the CPU-only reference may already require the whole host RAM ceiling, device allocation does **not** authorize host duplication without admission. The effective device safe ceiling is the minimum of 6 GiB and physically free device memory minus a measured runtime reserve; if this is insufficient for even one patch/tile, the optional GPU backend refuses launch or, if explicitly requested by the user and scientifically certified, uses a CPU backend at a legal boundary. No automatic hidden cross-device precision change or mid-stage fallback.

**Execution geometry:** benchmark tile sizes and GPU occupancy separately; initial tuning trials use contiguous groups of **32³-equivalent cells**, with kernel work decomposed into face and edge tiles. Evaluate one kernel per reconstruction / HLLD / EMF / update stage versus fused variants. Record launch overhead, per-stage device allocations, H2D/D2H bytes, transfer-overlap time, occupancy, achieved bandwidth and device-to-host synchronization count. Avoid global repeated allocations of `B_cell`, `J` and three coefficient vectors; stage-local derived work remains chunked.

**GPU comparability gate:** with the same input physics and reference profile, GPU and CPU must satisfy all parent scientific tolerances and match independent integrals within their declared floating-point equivalence tolerance. The default reference f64 authority cannot become f32 because a particular GPU delivers low f64 throughput. Lower-precision research modes require a separate version and a full accuracy/error-budget campaign.

# 12. I/O, checkpoints, analysis and provenance cost

Use the parent **science snapshot v7** and **checkpoint v25** MHD schemas. Science data may carry cell-centered derived `MagneticField`, while checkpoints carry authoritative face fluxes, patch offsets, geometry epochs and the physics profile. They are **not** interchangeable.

Chunked/collective HDF5 writing MUST consume canonical magnetic face arrays directly through hyperslabs or bounded chunks, without creating a second contiguous full-field mirror. Choose initial MHD I/O working-chunk target **8 MiB per rank**, admitted and counted within I/O phase resources (not the MHD 64 MiB communication pool unless the same allocation is actually reused). Actual chunks can shrink for topology or library restrictions. Publication retains existing CHUÍ `.partial`/`.complete` transaction rules and must remain checkpoint-safe at global magnetic synchronization boundaries.

Before snapshot creation, the memory governor must reserve HDF5 staging and library opaque runtime overhead, and confirm the intended I/O phase does not overlap an already saturated MHD stage. If the solver is near its hard ceiling, optional science-output staging may be deferred to the next legal output boundary; **mandatory checkpoint safety and defined output cadence must never be silently discarded**. If no legal safe output can be admitted, fail with an explicit budget error rather than writing a partial authoritative checkpoint.

Optional $B$-field projections, spectra and Faraday-rotation intermediates should use streaming or distributed analysis and release their scratch before the next peak MHD stage. Global root-gather of all magnetic faces is prohibited in production workloads. For observables requiring additional microphysical populations, the analysis provenance records their assumptions separately from the magnetic solver's fidelity claims.

All performance output extends the existing `profile_timeline.json`, `counters.json`, `provenance.json`/normalized configuration conventions where those are already used. No second isolated profiler format is made authoritative. Performance evidence must bind a source identity/commit or manifest hash, config hash, physics provider hashes, compiler, CPU/GPU model, MPI topology, patch shapes and numeric mode. Output JSON is permitted as machine-readable telemetry; simulation inputs remain `.param.txt`.

# 13. Performance telemetry: names, sampling, and accounting

**Required event identity:** `(run_hash, rank, node, level, patch_id_or_batch_id, accepted_interval_id, physics_stage, stage_attempt, accepted_flag, start_time, end_time, geometry_epoch)`. Counters are additive by phase, but a duration summed across ranks is **not** wall time; report both summed work time and maximum critical-path wall time.

| Metric name | Units | Definition |
|---|---|---|
| `mhd_active_cells` | count | Unique real cells touched by accepted stage |
| `mhd_active_faces` / `mhd_active_edges` | count | Unique geometric keys processed |
| `mhd_ideal_stage_count` | count | SSPRK2 operator evaluations, separately attempted/accepted |
| `mhd_hall_subcycle_count` | count | Accepted Hall explicit subcycles; also reject count |
| `mhd_rkl_stage_count` / `mhd_rkl_cycle_count` | count | Actual stages and cycles, not estimates |
| `mhd_hlld_fallback_faces` | count | Faces where HLLE/first-order recovery required |
| `mhd_coeff_eval_seconds` | s | Transport provider time, including table interpolation |
| `mhd_flux_seconds` / `mhd_emf_seconds` | s | Riemann flux and CT edge construction exclusive time |
| `mhd_ideal_seconds` / `mhd_nonideal_seconds` | s | Stage-inclusive and explicitly attributed time |
| `mhd_amr_reflux_seconds` / `mhd_regrid_seconds` | s | Boundary sync and topology transaction |
| `mhd_mpi_pack_seconds` / `mhd_mpi_wait_seconds` | s | Separate serialization and blocking critical path |
| `mhd_mpi_bytes_sent` / `mhd_mpi_bytes_recv` | bytes | Actual payload plus headers, not estimated field volume |
| `mhd_face_flux_capacity_bytes` | bytes | Actual persistent allocated magnetic capacity |
| `mhd_phase_scratch_high_water_bytes` | bytes | Sum of co-live MHD-owned phase/scratch allocations |
| `mhd_comm_high_water_bytes` | bytes | Per-rank aggregate simultaneous MHD communication buffers |
| `mhd_peak_rss_bytes` / `mhd_rank_max_rss_bytes` | bytes | Measured process/rank max; document RSS/PSS limitations |
| `mhd_device_high_water_bytes` | bytes | Governed GPU live set plus identified allocator-pool overhead |
| `mhd_host_device_bytes` | bytes | Directional H2D/D2H totals including pinned staging |
| `mhd_rejection_count` / `mhd_memory_admission_failures` | count | Distinguish numerical retries from memory rejections |
| `mhd_time_to_accuracy_seconds` | s | Full successful run wall-time at predetermined error target |
| `mhd_source_energy_residual` / `mhd_divergence_max_norm` | unit-declared | Existing physics-gate metrics bound to exact run |

Profiles MUST attribute the same physical byte range once to a memory class even if multiple views overlap. Where a process RSS is larger than tracked CHUÍ capacities, report the **unreconciled difference** (allocator fragmentation, C/C++ runtime, MPI/FFTW, mapping, GPU driver) rather than conceal it. Node-level PSS should be used when available to avoid falsely summing genuinely shared pages. Host sampling cannot guarantee observing every instantaneous peak; the governor reservation high-water and allocator-level capacity reports are independently mandatory.

**Sampling:** stage timer start/end records for **every** accepted and rejected stage; high-water governor counters captured on every allocation/lease transition; OS RSS at least at each accepted level interval and each communication/regrid/output transition. Higher-rate OS memory sampling may be enabled in profiling mode only; log its measurement overhead. Profiling itself must not allocate unbounded per-stage event vectors: use fixed-capacity ring buffers or bounded chunked streaming of events with deterministic dropped-event reporting. Aggregate telemetry at rank and node level, and retain fine-grained traces only for explicitly selected campaigns.

**Normalized throughput:** publish wall-seconds per accepted *physical interval*, accepted real-cell **operator** updates per second (ideal, Hall, RKL separately), successful simulated physical time per real wall hour, and the actual numerical accuracy. Never count a Hall subcycle or failed retry as a fully advanced cosmological timestep. Record each distinct stage family before constructing blended throughput measures.

# 14. Benchmarks and the Pareto frontier

## 14.1 Required workload suite

| ID | Workload and fixed geometry | Principal objective |
|---|---|---|
| `BP01` | 1D/3D smooth Alfvén wave; 32/64/128 cells per wavelength | Error/order versus wall-time and f64 bandwidth |
| `BP02` | Brio-Wu and low-beta magnetized blast, $32^3$–$128^3$ | Robust HLLD/HLLE stage time; fallback and energy cost |
| `BP03` | Orszag–Tang and field-loop advection; fixed/uniform meshes | Numerical dissipation at time-to-accuracy |
| `BP04` | $32^3$ patches tiled to $128^3$, ratio-2 AMR boundary crossing | Ghost, CT exchange, reflux-curl and AMR overhead |
| `BP05` | Ohmic sinusoidal decay; constant and spatially varying coefficients | Explicit/RKL2 crossover and coefficient evaluation |
| `BP06` | Hall whistler waves; signed coefficients, rising stiffness | Stable Hall cost curve and scheduler refusal |
| `BP07` | Full `D/H/I/H/D` non-ideal problem, valid table provider | Coupled source cost, stages, memory high-water |
| `BP08` | Representative magnetized cloud collapse | Science accuracy, source, AMR and timestep interactions |
| `BP09` | Magnetized cosmological $128^3$ and zoom-in toy workflows | End-to-end gravity/hydro/MHD overlap and resource fit |
| `BP10` | Allocation/regrid/migration stress at 16³/32³/64³ patches | Retained capacity, OOM avoidance, MPI worst rank |

Use analytic references or frozen high-resolution trustworthy numerical references as appropriate; inherit **G01–G09, I01–I10 and N01–N08** scientific gates from the parent physics specification. Each workload has both (a) same-physics, matched-resolution and (b) fixed-accuracy benchmarking. At least three mesh resolutions are required for empirical convergence where an analytic order claim is made. All run descriptors, initial conditions, input tables, comparators and statistical procedures must be versioned before timing starts.

**Cross-code controls:** use Athena++/AthenaK for ideal FV/AMR and IDEFIX for non-ideal uniform-mesh studies. AREPO is a separate moving-mesh physics comparator; do not claim its irregular-mesh advection advantage or cost has been isolated unless equations, initial data, boundary conditions, cooling, transport coefficients and error measurements are matched. Quokka's GPU MHD is an optional backend comparator only when the same numerical problem/precision has an established supported implementation.

## 14.2 Fixed environment and statistical method

- Build **release binaries** with compiler/version/flags, SIMD/FMA mode, MPI library, FFT backend, Kokkos/CUDA if applicable, and hash recorded. No debug/profiling overhead in primary performance number; collect profiling evidence in a paired instrumented run and quantify observer overhead.
- Record CPU model, P/E-core or NUMA topology, thread placement, governor/clock policy, RAM speed, accelerator, free memory, thermal conditions, rank mapping and I/O device. Avoid unexplained comparison of a GPU run with a different physical solution to a CPU run.
- Run one untimed warmup for each stable configuration, then at least **seven** independent wall-time samples for kernel/medium benchmarks with median and median absolute deviation. For long astrophysical campaigns, use at least **three** independent wall times when feasible and report uncertainty; never claim a tight confidence interval from fewer independent runs.
- For a claimed >10% improvement or regression, run at least five independent paired repetitions; compute a 95% confidence interval using a documented paired/bootstrap protocol. Promotion still follows parent Section 17: no safety regression, improvement above measurement uncertainty, no unexplained >5% peak memory cost on workstation class.
- Archive the full error vector, not just a visually attractive plot. Any run that fails physical validity, convergence, restart consistency or its budget is **ineligible** for Pareto ranking.

## 14.3 Reference performance qualification, distinct from scientific certification

The following are predeclared *performance engineering targets*. Failure leaves the numerically correct solver available but marks it `PERFORMANCE_NOT_QUALIFIED` on the corresponding hardware/workload; it does **not** allow relaxed physics tests or false claims of equivalence.

- **Scientific methods unchanged:** all parent REF-I/REF-NI numerical gates still pass with every enabled optimization.
- **Competitive CPU baseline:** on the same physical CPU configuration and with matched PDEs/accuracy, median time-to-accuracy is at most **2.0 times** Athena-family reference on at least two of the three ideal benchmark families (`BP01`, `BP02`, `BP03`). The end-to-end run includes solver kernel, halo and boundary costs; external code integration setup and I/O are excluded identically. If an external reference cannot be run on matched hardware, mark this item `not_evaluated`, never pass by assertion.
- **Competitive non-ideal CPU baseline:** median time-to-accuracy is at most **2.5 times** IDEFIX on at least two of `BP05`, `BP06`, `BP07` with genuinely matching transport models and resolution. This is a stretch target, not a statement that both codes already match.
- **Four-thread CPU parallelism:** for at least eight independent $32^3$ patches per rank and a stable compute-dominated ideal workload, measured four-thread speedup is **at least 2.5x** versus one thread on the same hardware, with no significant numerical regression. Failure requires a categorized bottleneck report (e.g. bandwidth, lock contention, insufficient independent work).
- **MPI weak scaling:** on 2/4/8 ranks with at least 32 independent $32^3$-equivalent patches per rank and matched communication topology, strong or weak-scaling efficiency must be reported. An initial **65% weak-scaling efficiency at 8 ranks** is the performance target, conditional on available homogeneous cluster hardware and stable problem geometry; otherwise mark `not_evaluated`. No performance claim substitutes for MPI scientific equivalence.
- **GPU option:** if a GPU backend exists, deliver either a measured **1.5x** or better median time-to-accuracy advantage over the same CPU node on a sufficiently large GPU-fitting workload, **or** document that the hardware's f64 capability/data-transfer requirements prevent it. GPU certification still requires scientific equivalence and 6 GiB target adherence; lack of speedup is a performance result, not a reason to demote precision.

External-code license and access restrictions are respected; the reference benchmark harness must support independently generated common ICs/results without copying proprietary source code. CPU/GPU qualifications are per named hardware class and cannot be generalized from one desktop or cluster.

# 15. Quantitative memory acceptance gates

| Gate | Mandatory check | Acceptance rule |
|---|---|---|
| `MP01` | Geometric accounting | Persistent f64 face capacities and stage f64 edge capacities agree with exact patch formulas (including mirrored boundary faces) within **1%** of separately measured *owned array* capacity after alignment; any remaining bytes itemized |
| `MP02` | Stage scratch | For $32^3$ PPM reference probes, owned phase+scratch peak **≤512 B per active cell**, counted over simultaneously active batches; no untracked unbounded per-worker buffers |
| `MP03` | MHD communication | Default combined MHD send+receive+reduction high-water **≤64 MiB/rank**; explicit approved high-memory profile **≤128 MiB/rank** |
| `MP04` | Node fit | Actual observed whole-host peak, conservatively reconciled with rank RSS/PSS and external allocations, is **≤42 GiB** for the 48 GiB profile, with the required runtime reserve present |
| `MP05` | Steady allocation plateau | After warmup and 100 repeated fixed-topology accepted intervals, capacity high-water does not grow by >**2%** absent documented cache tuning or intended state growth |
| `MP06` | Transaction coexistence | Regrid and migration tests admit exact old+new+scratch+wire coexistence; a deliberately insufficient budget rejects **before** any authoritative mutation |
| `MP07` | Legacy compatibility | With MHD disabled, no MHD face/edge/coefficient buffers exist; on identical DMO/hydro tests, peak resident memory increases by **≤1%** and median wall-time by **≤3%** versus pinned pre-MHD baseline, subject to measurement variance |
| `MP08` | I/O staging | MHD HDF5 extra per-rank staging **≤8 MiB** initial chunk target, or smaller if reserve dictates; no full-population duplicate output representation |
| `MP09` | Scaling laws | Patch tests at 16³/32³/64³ show persistent face cost follows exact $N_f$ slope; communication buffers plateau at configured ceiling under increased local populations |
| `MP10` | GPU budget (if enabled) | Governed device live set **≤6 GiB** for 8 GiB-class target and host pinned buffers included in host ceiling; CPU/GPU numerical gates pass |
| `MP11` | Disabled unlimited governor | Every `PERFORMANCE_QUALIFIED` run has a positive finite governor hard limit and documented external-runtime reserve; no `hard_limit_bytes=0` bypass |
| `MP12` | State/telemetry truth | Rank maxima, per-class capacities, reservation high-water, peak measured process memory, and untracked delta are present; no fake measurements or silently ignored telemetry drops |

**A subtle but essential caveat:** `MP01` compares **owned array capacity**, not whole-process RSS, because allocator overhead and unrelated runtime state make a 1% RSS equality scientifically meaningless. `MP04` checks the actual process/host envelope, while `MP01` checks the correctness of our analytical MHD array model.

When the 512 B/cell cap or 64 MiB/rank default proves suboptimal, we may benchmark smaller tiles, different scratch reuse, or a versioned high-memory profile. We cannot simply replace the budget with an unmeasured “auto” allocation that consumes available memory without bounded predictions.

# 16. Performance optimization and change-control policy

The order of optimization work is deliberate:

1. Remove redundant full-population representations, hash-driven per-face operations and repeated allocations. Establish correct cache- and capacity-accounting first.
2. Improve contiguous patch traversal, vectorization, stage-local fusion, edge ownership, and false-sharing behavior. Keep exact scientific operators.
3. Benchmark alternative patch/tile widths and allowed scratch/caching policies; choose time-to-accuracy winners within fixed resource envelopes.
4. Improve MPI boundary packing and overlap, controlled concurrent batches, and memory-aware decomposition.
5. Port existing kernels to GPU with exact f64 reference semantics and separately measured transfer costs.
6. Only in a **research profile**, investigate compressed/mixed precision, asynchronous numerical stages, high-order compatible-DG, vector potential or implicit Hall. Such experiments have different schema IDs and may replace the reference **only** by meeting the parent scientific and Pareto-promotion rules.

Every performance PR must attach a **before/after** table specifying source hashes, exact numerical mode, hardware, physics error vector, magnetic-divergence maximum, energy residual, per-stage wall-time, rank-max RSS, stage scratch, MPI bytes and reproducibility caveats. Faster at the cost of unaccounted heating, extra divergence, or a different physical coefficient model is a regression, not an optimization.

**Reproducibility:** reference deterministic mode uses stable geometric-key ordering for rank-local shared-edge reductions and records execution flags; it does not promise bitwise identical answers across different GPU/CPU arithmetic architectures. Within the same supported execution environment, restart equivalence and cross-rank scientific tolerance tests remain required. Fast non-deterministic accumulation is a research/performance variant that must be measured and explicitly identified; it cannot be used as the canonical reference silently.

# 17. Configuration contract and proposed defaults

Extend the parent `.param.txt` format with **one** `[mhd_performance]` section and, only where the existing parser already accepts its established keys, reuse existing `[memory]` and `[parallel]` policy. All new keys are **proposed parser additions**, not current runnable settings. All quantities in this section are normalized into a frozen run-config hash and included in provenance and restart compatibility checks where they affect execution semantics.

```ini
[mhd_performance]
profile = workstation_48gib_v01
patch_probe_size = 32
patch_tile_target_cells = 32768
ppm_ghost_cells = 3
plm_ghost_cells = 2
active_batch_target_patches = 8
mhd_scratch_max_bytes_per_active_cell = 512
mhd_comm_budget_mib_per_rank = 64
mhd_comm_max_mib_per_rank = 128
mhd_io_chunk_mib_per_rank = 8
hall_subcycle_cost_cap_per_level_interval = 4096
rkl_stage_cost_cap_per_level_interval = 4096
stage_arena_null_upstream = true
track_rank_max_memory = true
track_phase_high_water = true
track_mpi_bytes = true
record_time_to_accuracy = true
optional_gpu_backend = false
gpu_device_target_mib = 6144
allow_untracked_allocations = false
```

**Interpretation:** `patch_probe_size` is only a benchmark/test shape and does **not** override existing AMR geometry. `patch_tile_target_cells` is an advisory tile target; actual patch sizes and stencil dimensions win. `active_batch_target_patches` is an upper bound for a $32^3$-equivalent group, not a required allocation when work or headroom is smaller. `hall_subcycle_cost_cap_per_level_interval` and `rkl_stage_cost_cap_per_level_interval` constrain scheduling cost, not physics. `optional_gpu_backend=false` keeps the CPU reference runnable even on systems without CUDA/Kokkos.

**Validation rules:** any missing mandatory numeric value, negative budget, contradictory cap, unknown key, unsupported profile or mismatch between memory reserve and GPU/host discovery fails at preflight. A finite node/host hard ceiling is mandatory in qualification profiles; the existing general governor remains the sole admission authority. Existing CHUÍ configuration may use different section/key spelling; implementation must map to typed, versioned fields without creating a second ungoverned policy object. No duplicate YAML/JSON config parser.

Profile variants are fixed as follows: `workstation_48gib_v01` uses 42 GiB host safe/36 GiB nominal targets and the caps above; `cluster_v01` derives the host safe limit from the machine plus explicit runtime reserve and otherwise retains the same MHD scratch/communication defaults unless amended; `gpu_8gib_v01` enables GPU with a 6 GiB live device target and the same f64 physics, while host admission remains independent. Physical MHD profile is always identified separately (`mhd_ref_ideal_v01` or `mhd_ref_nonideal_v01`).

# 18. Failure behavior and nonnegotiable limits

| Failure | Required response |
|---|---|
| MHD governor reservation denied | Reduce *next* batch or optional overlap before work begins; otherwise `MHD_MEMORY_BUDGET_EXCEEDED` with owner and byte details |
| One patch/tile cannot fit | Fail preflight with minimum required owned, scratch, ghost, wire, temporary coexistence and external-reserve breakdown |
| Unexpected allocation growth | Fail with owner and call site; prohibit ungoverned heap fallback |
| MPI buffer/record exceeds bound | Chunk deterministically; if a single indivisible semantic item is oversized, reject with structured diagnostic |
| Rank/edge geometry epoch mismatch | Collective rollback/abort; `MHD_MPI_EPOCH_MISMATCH`, no corrupted partial authority |
| Hall or RKL predicted cost exceeds operational limit | Ask global scheduler to shorten legal interval; if not possible, `MHD_STIFFNESS_BUDGET_EXCEEDED`; no physics coefficient alteration |
| GPU allocation/launch failure | At legal boundary only, use explicitly allowed certified CPU backend if policy permits; otherwise fail safely, never change precision or continue a partial stage |
| Invalid energy/divergence/positivity | Parent scientific failure and retry policy; performance profile cannot override |
| Mandatory checkpoint cannot fit | Controlled failure at safe boundary; no incomplete or fabricated restart authority |
| Optional diagnostic cannot fit | Defer it with explicit event and provenance; never silently disable mandatory scientific evidence |

Every error includes `run_hash`, owner/rank/node, geometry and time/scale-factor epoch, requested/committed/reserved bytes, safe limit, predicted stage count where relevant, and the last accepted-state identity. A resource failure is **not** grounds to silently select a lower-order reconstruction or a different physical transport model.

# 19. Implementation and qualification sequence

| Milestone | Deliverable | Gate |
|---|---|---|
| `P0` Resource contract | Complete `MhdMemoryEstimator`, typed perf config, device/host budget and topology-aware admission | `MP01`, `MP04`, `MP11` preflight and failure-injection tests |
| `P1` Patch CPU kernels | Governed stage arena, f64 SoA, contiguous HLLD/UCT/CT trial update, timers | `MP02`, `MP05`, `MP07`, parent ideal physics gates |
| `P2` AMR/MPI execution | Bounded face/edge communication, owner keys, edge registers, migration/regrid transactions | `MP03`, `MP06`, `MP09`, parent MPI/AMR gates |
| `P3` Non-ideal stages | Reusable RKL/Hall stage state, coefficient provider caching, stiffness forecast | `MP02`, `MP05`, cost cap, parent non-ideal gates |
| `P4` Production telemetry | Integrated profiler, rank maxima, I/O admission, trace retention | `MP08`, `MP12`, checkpoint integrity, full benchmark suite |
| `P5` Competitive CPU review | Athena-family and IDEFIX matched baseline campaigns | CPU qualification and published accuracy/time/memory plots |
| `P6` GPU and research | Device governor, transfer accounting, optional optimized/novel kernels | `MP10`, CPU/GPU equivalence, benchmark/Pareto promotion rules |

**Release identities:** `REF-I NUMERICALLY CERTIFIED` and `REF-NI NUMERICALLY CERTIFIED` retain their precise parent meanings. This companion adds `MHD-PERF-WORKSTATION-v1 QUALIFIED`, `MHD-PERF-CLUSTER-v1 QUALIFIED`, and `MHD-PERF-GPU8-v1 QUALIFIED`. A label may be emitted only when its relevant `MP` and comparative performance gates have archived passing evidence on the named hardware class. At issuance of this design, **none** of those labels is asserted to be achieved.

# 20. Requirements traceability and research boundaries

| Contract | Owning module/interface | Acceptance evidence |
|---|---|---|
| Magnetic/edge geometric costs | `mhd_state`, `mhd_geometry`, `MhdMemoryEstimator` | Exact array counts and MP01/MP09 |
| Single authority and stage snapshots | `SimulationState`, MHD stage runtime, governor | MP02/MP05, restart/retry equivalence |
| Memory admission and batching | `MemoryGovernor`, `GovernedScratchArena`, MHD workflow | MP04/MP06/MP11, failure injection |
| MPI ownership and bounded exchange | `parallel`, `mhd_amr`, `mhd_induction` | MP03, parent G06/G08, rank-imbalance probes |
| Stiff transport cost controls | `mhd_transport`, existing time scheduler | Predict/actual stage-cost report, N01–N07 |
| Device memory and parity | `mhd` GPU backend, device governor | MP10, parent physics gates |
| Science I/O and restart | Existing I/O workflow, new v7/v25 MHD schema | MP08, parent G08 |
| Pareto comparisons | `bench/mhd`, versioned reference runner | BP01–BP10, paired timings and archived errors |
| Experimental numerical methods | Independent research profile owner | Parent Section 17 promotion; never implicit in REF |

**Scientific boundaries not changed here:** This document does not provide a new Hall integrator, DG solver, vector-potential gauge, or transport-chemistry closure. Those remain separate *research directions* explicitly excluded from reference-v1 physics. It does completely specify the **memory and performance contracts for the selected reference algorithms**, including default budget, ownership, stage coexistence, MPI/GPU admission, output costs, telemetry, benchmark procedures, change management and failure handling. No performance advantage is claimed without measurements.

# Appendix A. Worked host sizing examples

### A1. One $32^3$ active PPM reference batch

- Persistent magnetic flux: $3(33)(32^2)\times8=811{,}008$ B = **0.773 MiB**.
- Complete stage edge-EMF array: $3(32)(33^2)\times8=836{,}352$ B = **0.798 MiB**.
- Basic six-double gas trial image: $6(32^3)\times8=1{,}572{,}864$ B = **1.500 MiB**.
- Two separate trial images for gas + face flux: $2(1.500+0.773)=4.547$ MiB, included **inside** the stage scratch budget.
- Full PPM $g=3$ ghost cell cube if explicitly materialized: $38^3/32^3=1.675$ times the real cell count; tile-shell reconstruction is preferred.
- MHD phase resident + scratch engineering cap: $32768\times512=16$ MiB, subject to actual cell auxiliary field costs.

These figures do **not** include canonical gas sidecars, geometry metadata, MPI wire buffers, persistent flux registers, table caches, FFT/PM, TreePM or the OS. They must not be added to the whole-node total without checking co-live ownership (for example, persistent face flux is already counted separately from the two **trial** face images).

### A2. A uniform $512^3$ magnetic population

One global staggered f64 magnetic representation is ~**3.006 GiB**; 4,096 independent $32^3$ patch representations total **3.094 GiB** due to patch face duplication. A single extra persistent three-component cell-centered $B$ cache would add **3 GiB** on top, with no new scientific authority. A full-population 512 B/cell stage buffer would be **64 GiB**—already impossible on the 48 GiB reference workstation—so **only active patch/tile batches** may hold that scratch. The 48 GiB program's long-term success does **not** imply that a $512^3$ full-physics uniform hydro+MHD run is feasible; zoom-in and AMR are the intended scientific strategy.

### A3. One eight-patch active CPU group

Eight $32^3$ patches have 262,144 active real cells. MHD phase+scratch at 512 B/cell is **128 MiB**. Magnetic persistent state across those eight patches is **6.19 MiB** (charged once to canonical ownership), stage edge EMFs would be **6.38 MiB** if all eight full arrays are simultaneously live (already inside the 128 MiB phase allowance), while the default MPI communication pool remains separately capped at **64 MiB/rank**. The system must still admit other co-live gravity/hydro, source, runtime and host-node memory before starting the group.

# Appendix B. Performance-code-review checklist

Every source change materially affecting MHD allocations, stage execution, communication or device dataflow MUST answer:

1. What is the owner, `MemoryClass`, lifetime, and exact shape/type/stride of each new allocation?
2. Which allocations can coexist in the worst legal RK, RKL, Hall, regrid, migration, MPI or checkpoint state?
3. What is the calculated byte formula at 16³, 32³, 64³ and a 512³ global magnetic scale?
4. Which capacity or scratch lease is reusable, and what invalidates it?
5. Is any data another copy of canonical gas or magnetic authority? If so, why and for how long?
6. Does the communication path stay bounded for zero-work ranks, fragmented AMR, count overflow and uneven ownership?
7. What are the exact before/after wall times, physics-error vectors, peak capacities, RSS, MPI bytes and failure counts?
8. Does the code preserve the same f64 reference equations and energy/AMR/restart contracts?
9. Is an optional GPU/device allocation budgeted alongside its host staging and runtime overhead?
10. What tests prove both the predicted resource ceiling and unchanged scientific fidelity?

# References and source basis

- **[P1]** `chui_mhd_architecture_specification_v01.md`, CHUI-MHD-REF v1.0 (22 Nov 2026), companion governing equations, time integration, AMR, I/O and physical acceptance gates.
- **[P2]** `chui_memory_architecture_and_resource_governance_specification_v01(1).docx`, 23 Aug 2026. Governs 48 GiB host envelope, budget classes, hard admission, low-overhead stage scratch, MPI and precision doctrine.
- **[P3]** CHUÍ public repository, `main`, inspected 22–23 Nov 2026: https://github.com/AlexandreKoz/CHUI-Cosmology-Simulator . Relevant live interfaces: `include/cosmosim/core/memory_governor.hpp`, `include/cosmosim/core/governed_scratch_arena.hpp`, `include/cosmosim/parallel/distributed_memory.hpp`, `include/cosmosim/hydro/hydro_core_solver.hpp`, `docs/architecture/amr_hydro_subcycling_and_flux_registers.md`, `docs/architecture/adr_param_txt_configuration_format.md`.
- **[P4]** `arepo lessons.pdf`; architectural lessons on data motion, hot/cold separation, irregular mesh geometry, cache locality and multiphysics coupling.
- **[P5]** `CosmoSim_Master_Design_Document_Revised-1.pdf`, 11 Mar 2026; desktop-first, modular patch AMR, verification and HPC scale-up priorities.
- **[R1]** Stone et al. 2020, *The Athena++ Adaptive Mesh Refinement Framework*, ApJS 249:4, https://doi.org/10.3847/1538-4365/ab929b . Reference class for FV/CT/AMR.
- **[R2]** Lesur et al. 2023, *IDEFIX*, A&A 677:A9, https://doi.org/10.1051/0004-6361/202346005 . Non-ideal MHD and performance portability.
- **[R3]** Gardiner & Stone 2008, *An unsplit Godunov method for ideal MHD via CT in 3D*, JCP, https://doi.org/10.1016/j.jcp.2007.12.017 . Reference unsplit MHD/CT design.
- **[R4]** Mignone & Del Zanna 2020, *Systematic construction of upwind constrained transport schemes for MHD*, JCP, https://doi.org/10.1016/j.jcp.2020.109748 . Upwind EMF design.

**End of CHUI-MHD-PERF v1.0.** This is a proposed contract for the next engineering campaign, not a report of completed measurements or an assertion that full-physics workstation scaling has been demonstrated.
