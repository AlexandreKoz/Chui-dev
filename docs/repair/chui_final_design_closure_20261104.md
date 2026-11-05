# CHUI final design v1.0: implementation-closure ledger and handoff (2026-11-04)

**Evidence status: SOURCE ONLY. Nothing was configured, compiled, linked, tested, run,
benchmarked or sent to CI. Every test below is NOT RUN — USER-RESERVED.** "Implemented in
source" means only that a static read confirms the production call path, the memory
owners, the communication protocol and the refusal paths. It says nothing about
numerical validity or scientific qualification.

| Item | Value |
|---|---|
| Base | `main` @ `de1eac50dccbb93a9480956ed351876c548ad827`, re-fetched and confirmed as the latest remote HEAD before work began |
| Branch | `claude/epic-fermi-2bcq65`, rooted at the base; only ordinary commits and pushes, no force-push, rebase, merge or PR |
| Prior ledger | [`chui_conformance_closure_20261104.md`](chui_conformance_closure_20261104.md): its A1/A2 Codex repairs are preserved, and its A3–A5 continuation is closed below |

## 1. Commit sequence

| SHA | Gate | Content |
|---|---|---|
| `3ea0889` | G0.1, G0.4 | `SimulationConfig` member is now ctor-initialized, and so is a test that default-constructed it. The chemistry promotion bias is derived from the canonical sidecar. |
| `ce8c3da` | G0.3 | Associations are found by global FoF before the remote stencil query. Dispersion kinetic energy (KE) is conserved. Partition-invariance tests added. |
| `3e5c694` | G2.1, G0.5 | Distributed element-resolved diffusion with a unique flux authority. One complete diffusion estimator. Topology construction is admitted. |
| `b15597b` | G1.1, G1.2, G2.3 | Distributed recursive AMR, level-relaxed CFL with a limiter record, and the split-interface reflux defect fix. |
| `6102e6c` | G1.1/G2.3 tests | Loopback-hook protocol, peer refusal, and register-merge regressions. |
| `afd75fd` | G1.6 | Source-controller rollback checkpoint added to the whole-step transaction. |
| `2897999` | G2.2 | Octet families are co-located in every migration plan. |
| `7495941` | G3.2 | SSP cumulative return is structurally monotone (fixed nodes). |
| `4c70b0d` | G3.4 | Zoom contamination is emitted in operational reporting. |
| this commit | docs | This ledger, configuration and status docs. |

## 2. Requirement-closure matrix (H01–X01)

| ID | Status | Basis and exact remaining work |
|---|---|---|
| H01 Godunov/PPM/CTU | **PARTIAL** | The solver suite is unchanged. G3.1 CTU transverse reconstruction at MPI/coarse-fine boundaries is not reworked. Next edit: `amr_ghost_fill.cpp` deep-layer extrapolation (`deep_layer_extrapolated_ghosts`) for remote sparse sources. Owner evidence is needed first. |
| H02 dual energy / comoving | **IMPLEMENTED IN SOURCE — NOT RUN** (remote auxiliary lane) | The remote sparse ghost record now carries the owner's specific internal energy (`AmrHydroSparseRemoteCell::internal_energy_code`, filled from `AmrPatchCellPayloadRecord::internal_energy_code`, and interpolated with the conserved weights for recursive remote coarse ghosts). The AMR sweep seeds rho_ghost * e for remote-sourced ghosts, so the E−K fallback (counted in `dual_energy_unresolved_ghosts`) remains only for genuinely unknown sources. Missing test: a dual-energy remote-boundary assertion in `test_amr_distributed_remote_patch_boundary.cpp` that `dual_energy_unresolved_ghosts == 0`. |
| C01 thermochemistry | **BLOCKED BY EXTERNAL DATA** (R18 / photodissociation heating tables) | Source contract unchanged. |
| C02 network references | **EXPERIMENTAL/UNQUALIFIED** | Unchanged. |
| C03 adaptive chemistry control | **IMPLEMENTED IN SOURCE — NOT RUN** (transaction); control law **EXPERIMENTAL/UNQUALIFIED** | G0.4: the single authority is the `chui_adaptive_chemistry_control_v1` sidecar. It is re-read before every fidelity decision. Empty ranks contribute zero current samples, and a malformed sidecar is an error. Rollback and restart restore it with `SimulationState`. Timescale/stiffness indicators are still heuristic (G4.3 PARTIAL). |
| C04 shadow verification | **EXPERIMENTAL/UNQUALIFIED** | The bounded sampled shadow never mutates state (unchanged). |
| I01 effective/resolved ISM | **EXPERIMENTAL/UNQUALIFIED** | Unchanged. |
| I02 dust/SGS | **EXPERIMENTAL/UNQUALIFIED** | Dust lanes use the distributed operator. SGS now includes staged traceless resolved-stress production/transfer on the same AMR/MPI face topology; calibration and numerical qualification remain pending. |
| S01 star formation | **IMPLEMENTED IN SOURCE — NOT RUN** | Unchanged. The particle-ID occupancy cache is invalidated on rollback (G1.6). |
| S02 SSP / events | **IMPLEMENTED IN SOURCE — NOT RUN** (yields **BLOCKED BY EXTERNAL DATA**) | G3.2: the fixed-node non-negative linear interpolant makes C(age) monotone by construction and lets interval returns telescope. The clamp is now only a counted round-off guard. Behaviour change: returns move at the quadrature-error level. |
| F01 mechanical feedback | **IMPLEMENTED IN SOURCE — NOT RUN** | G0.3 (see §4). The remote stencil query now also covers wind-only and dispersion-only sources. |
| F02 winds | **BLOCKED BY EXTERNAL DATA** (wind table) | The coupling path exists. Its energy and momentum enter the deposition energy closure. |
| D01 element diffusion | **IMPLEMENTED IN SOURCE — NOT RUN** | G2.1 (see §4). The legacy total-Z scalar operator is explicitly single-rank (refused on more than one rank). |
| A01 regrid | **IMPLEMENTED IN SOURCE — NOT RUN** | G2.2: octet families are co-located by every migration plan, so the split-sibling skip is transient. Parent-link recovery for twice-derefined grandparents keeps its documented limitation (link 0 means it is not derefined further; the id is never wrong). Whole-hierarchy regrid transaction: the rung-zero global rollback covers it. |
| A02 recursive AMR | **IMPLEMENTED IN SOURCE — NOT RUN** | G1.1 (see §4). |
| A03 physics indicators | **PARTIAL** | G4.1 is not implemented this session: the cooling-length, chemical-interface, vorticity and force-resolution indicators are missing. Next edit: `RefinementEvaluator` / `amr_production_regrid.cpp::patchDecision`, keeping a mandatory/optional split. |
| A04 predictive AMR | **EXPERIMENTAL/UNQUALIFIED / PARTIAL** | G4.2: the young-star predictor exists. The event-aware forecast against `SspModel` event times is not implemented. |
| G01 TreePM | **IMPLEMENTED IN SOURCE — NOT RUN** | TreePM remains the sole gravity authority. G3.4: zoom contamination is reported through event `analysis.zoom_contamination`. No fine-gravity solver was added because no evidence justifies one. |
| T01 rollback | **IMPLEMENTED IN SOURCE — NOT RUN** (full-copy reference; not a bounded journal) | G1.6 (see §5). A bounded delta journal remains an optimisation, not a correctness gap. |
| T02 multirate | **IMPLEMENTED IN SOURCE — NOT RUN**, under the coarse-cadence operator contract of §5 | G1.2–G1.5. Per-level source and event integration is deliberately NOT introduced (§5). |
| M01 memory | **IMPLEMENTED IN SOURCE — NOT RUN** | G0.5 (see §6). Capacity accounting is not RSS. |
| M02 performance | **PARTIAL** | No new benchmark kernels this session (§8). |
| V01 outputs/restarts | **PARTIAL** | No schema change this session. Science composition/star snapshots remain **export-only**; they are NOT interchangeable with restart. Restart covers `SimulationState`, including composition, populations and sidecars. Snapshot import is absent. |
| V02 research control | **EXPERIMENTAL/UNQUALIFIED / PARTIAL** | Unchanged apart from the C03 transaction fix. |
| X01 scope | **OUT OF SCOPE BY DESIGN** | No MHD, RT, CR, AGN framework, relativistic or nuclear physics, dark-sector model, or SUBFIND was added. |

## 3. Gate summary

* **Gate 0**: closed in source. G0.1, G0.3, G0.4 and G0.5 are implemented. G0.2: the A1/A2 Codex
  repairs were re-read and are unchanged. The distributed recursive driver reuses
  their admission (`amr.recursive.driver`) and their atomic history and register publication.
  All new driver preparation is reduced collectively before mutation.
* **Gate 1**: G1.1 and G1.2 are implemented. G1.3–G1.5 are implemented under the coarse-cadence contract
  (§5). G1.6 has the full-copy transaction plus the source-controller checkpoint.
* **Gate 2**: G2.1, G2.2 and G2.3 are implemented, including the split-interface reflux defect.
* **Gate 3**: G3.2 and G3.4 are implemented. G3.1 and G3.3 are PARTIAL (see the matrix).
* **Gate 4**: PARTIAL. No new controllers. Existing research controllers stay marked research.
* **Gate 5**: the tests in §8 are written and registered. Four-mode decks predate this session.
  No recursive-AMR four-mode deck and no benchmark sources were added (PARTIAL).

## 4. Production call-graph changes

```
TimeCoordinator::runRungZeroSegment
  updateAdaptiveTimeBins: allreduce min AMR level -> per-cell hydro-CFL relaxation r^(L-Lmin)
                          (recursive only) -> TimestepLimiterRecord -> event time.timestep_limiter
  GlobalStepRetryPoint capture: + SourceRuntime::captureRollbackCheckpoint
  rollback: + restoreRollbackCheckpoint (ID occupancy cache invalidated); accept: discard
HydroAmrRuntime::runProductionAmrHydroPath
  recursive: verifyAcceptedRecursiveAmrHydroCfl (collective)
  MPI && recursive  -> advanceDistributedProductionAmrHydroRecursive(..., DistributedRecursiveHooks)
       prepare (all-level directed exchange; cross-rank 2:1 check)
       per level clock: beginLevelSweep (levels {L, L+1} exchange; retain L start;
                        interpolate L-1 remote coarse to t) -> level sweep
                        -> exchangeLevelFluxRegisters (merge into owner's pending store)
                        -> captureLevelEnd -> children -> reflux -> retireLevel
  MPI synchronized -> exchangeDirectedHydroGhosts (extracted) -> sweep
                      -> deferred local + inbound halves -> mergeFluxRegisterEntriesByKey -> apply
SourceRuntime
  executeThermochemistry: loadCommittedChemistryControl(state) before staging
  executePart2StellarPopulation: budgets -> mergeAssociationsDistributed
                                 -> importRemoteFeedbackCandidates -> stage -> remote deltas
  executeDiffusionAndSgsShear: gatherRemoteDiffusionCells -> topology (remote faces, global envelope)
                         -> DistributedElementDiffusionCoupling::setup -> advanceElementResolvedDiffusion
MigrationBalanceRuntime::applyMeasuredRuntimeRebalancePlan
  -> enforceAmrFamilyColocation (also runs when load rebalance is disabled)
AnalysisRuntime::executeDiagnostics (run-health cadence) -> reportZoomContamination
```

**G1.1 details.** `RecursiveAmrDistributedHooks` lives in `amr_hydro_orchestrator.hpp` and is
MPI-free. `DistributedAmrRemotePatch::state_time_code` labels each remote state with its level
clock. The ghost filler accepts a sparse remote coarse source under temporal coarse-to-fine only
when the source's label equals the fill time exactly. A stale end-state remote coarse ghost
therefore throws instead of silently filling.

**G0.3 details.** `physics::discoverAssociationLeaders` (hash grid plus union-by-min-id) and
`mergeAssociationMembers` (id-ordered sums) operate on an admitted global descriptor all-gather.
Each member is routed to its leader's owner rank. The resulting merged budget is bitwise
independent of the partition.

**G2.1 details.** The flux authority is the side with the smaller `gas_cell_id`. A claim and
conductance handshake validates a one-to-one interface; any asymmetry is a collective error.
Each subcycle performs two exchanges: fresh values go to the authority, and transfers come back.
Stable dt, subcycle count, positivity and conservation verdicts are global, and all ranks commit
together. The Smagorinsky strain includes remote neighbour velocities.

## 5. Conservation, time-integration and rollback contracts

* **Hydro.** Each level-L cell integrates exactly once over its own interval. The (L, L+1) reflux
  is applied once, after all fine substeps, from register halves on any ranks. Pending-store
  substep coverage masks prevent double counting. Split synchronized interfaces are now refluxed
  (G2.3 fix).
* **Bitwise versus round-off.**
  * Association merges: bitwise partition independent.
  * Recursive hydro: bitwise equal to the local driver when no remote data exist (loopback test).
    Across decompositions it differs at round-off level, because register halves are summed in
    source-rank order.
  * Diffusion: lane totals are conserved to round-off, and the result is decomposition
    independent to round-off. Bitwise identity is NOT claimed.
* **Gravity (G1.3).** TreePM KDK half-kicks bracket the whole coarse interval for every level.
  This is a Strang gravity/hydro split at the coarse cadence. No acceleration is used inside a
  fine substep, so no stale acceleration is consumed. The gravity, expansion and acceleration dt
  candidates are not level-relaxed, so the coarse interval respects them for the finest cell.
  Collisionless hierarchical DMO KDK is unchanged.
* **Sources and events (G1.4/G1.5).**
  * Hydro source terms passed to the sweep (cooling, effective-ISM relaxation) run per level
    substep with that substep's update.
  * `SourceRuntime` operators run once per accepted coarse interval after the hierarchy
    synchronizes: Part 2 thermochemistry, dust/SGS, star formation, SSP return/events, feedback
    and diffusion.
  * Every cell's interval is covered exactly once with a positive proper-time duration
    (`makeSourceInterval`).
  * Events use the cumulative-Poisson chronology on the restored population sidecar, so
    duplicates across retries are impossible. Deposition uses the synchronized step-end gas state.
  * Their stability constraints are not relaxed. This is a deliberate first-order coarse-cadence
    contract, not a per-level split.
* **Rollback (G1.6).**
  * Rejection reasons that roll back the step: the collective `kNumericalRetry` or
    `kPeerRejected` verdict. Any other kind is rethrown on every rank.
  * What is restored: the full governed `SimulationState` copy (composition, populations and
    ordinals, chemistry sidecar, histories, registers), the integrator, both schedulers, pending
    output, the force endpoint, and the source-controller checkpoint.
  * HydroAmrRuntime holds only diagnostics and key-validated caches.
  * A capture refusal on any rank is collective, before any mutation.
  * Retry halves dt and keeps the same physical model.

## 6. Memory ownership summary (governor owners added)

| Owner | Covers |
|---|---|
| `sources.associations.global_discovery` | Gathered descriptors plus discovery scratch |
| `sources.element_diffusion.workspace` | Single estimator including delta |
| `sources.metal_diffusion.topology_construction` | Exact-count bound with node-container allowances, admitted before construction |
| `sources.element_diffusion.remote_cells` | Per-peer receive |
| `sources.element_diffusion.remote` | Claims, conductances, one subcycle's state/delta records |
| `hydro.recursive.remote_interval` | Retained remote coarse start/end states |
| `hydro.recursive.remote_interpolated` | Per-sweep interpolated copies |
| `migration.amr_family_colocation` | Family triples |

* The existing node-wide `shared_equal` / `shared_pool` policy is used unchanged. No independent
  quota is introduced.
* **Bounded unadmitted allowances** (documented, O(interface) per sweep or step): the per-sweep
  outbound register-half vectors in the recursive driver, the association member-record exchange
  vectors, and the O(lanes) vectors in the diffusion coupling.
* Allocator overhead stays in the opaque process reserve.

## 7. Compilation-risk notes (static only; nothing compiled)

* **Resolved:** the deleted `SimulationConfig()` used as a member (`source_runtime.cpp`) and in
  `test_hydro_cfl_stability_rule.cpp`.
* **Watch at the first build:**
  * `HydroAmrRuntimeImpl::DistributedRecursiveHooks` is a nested class that uses enclosing
    private members declared later. This is legal in complete-class contexts.
  * Designated-initializer order in the new `TimestepLimiterRecord`,
    `AmrInterfaceTopologyOptions` and `PatchBoundaryFace` literals.
  * `std::erase_if` on `std::vector` (C++20).
  * `allreduceSumUint64sInPlace` with a `std::array<uint64_t,256>` span.
  * Defaulted trailing parameters of the anonymous-namespace functions in `time_coordinator.cpp`.

## 8. Tests written and registered (ALL NOT RUN — USER-RESERVED)

| CTest name | Purpose |
|---|---|
| `integration_distributed_element_diffusion_mpi_np2`, `_np4` | Contiguous, interleaved and empty-rank decompositions, open and periodic, against the single-domain reference (1e-12). Global lane conservation (1e-13). One authority per interface. Collective refusal without deadlock leaves state bitwise unchanged. |
| `integration_amr_family_colocation_mpi_np2`, `_np4` | Split-family consolidation. The family follows its minimum-id child. Co-located families are untouched. |
| `unit_amr_recursive_subcycling` (extended) | Loopback distributed driver is bitwise equal to the local driver, with the hook sequence checked. Peer-only level traversal. Peer refusal before mutation. Split-register merge and refusal. |
| `unit_chui_part3_remote_sources` (extended) | FoF chaining across the old lattice boundary. Permutation and partition invariance. Mass, momentum and lab-frame KE conservation. Duplicate refusal. Periodic wrap. |
| `unit_workflow_metal_diffusion_topology` (extended) | Remote face geometry, side and owner. Strain from remote velocity. Bound ≥ measured scratch. Direct and periodic faces kept distinct. |
| `unit_stellar_population_feedback` (extended) | 200-way SSP subdivision with zero clamps and per-element positivity. |

**Not written this session** (owner priority list):

* 2-/4-rank distributed recursive AMR end-to-end against a 1-rank run. This needs an AMR deck with
  `numerics.amr_time_stepping = recursive_subcycling`.
* Injected numerical-retry workflow test comparing the source checkpoint bit-for-bit.
* Association MPI end-to-end test.
* Benchmark sources for recursive overhead and remote diffusion.

Owner commands (documentation only):

```bash
cmake --preset cpu-only-debug
```

```bash
cmake --build --preset build-cpu-debug
```

```bash
ctest --preset test-cpu-debug --output-on-failure -L chui_part3
```

Then run the MPI presets from `docs/build_instructions.md` for the `*_np2` / `*_np4` tests.

## 9. External data that genuinely remains unresolved

* Production SSP yield tables (`chui_ssp_yield_table_v1`).
* Stellar wind tables (`chui_ssp_wind_table_v1`).
* R18 / photodissociation-heating rate data.

All are checksum- and reference-bound. Nothing was fabricated.

## 10. Precise continuation (same branch)

1. **G3.1:** CTU transverse reconstruction at remote and coarse/fine boundaries. Add the H02
   remote dual-energy regression.
2. **G4.1:** cooling-length, chemical-interface, vorticity and force-resolution indicators in
   `patchDecision`. Keep mandatory accuracy and resource priority separate. Never coarsen below
   a mandatory criterion.
3. **G4.2:** event-aware predictor from `SspModel` cumulative event times, at legal regrid
   barriers only.
4. **Tests:** the recursive end-to-end MPI deck, plus injected-retry and association MPI tests.
5. **Benchmarks:** recursive-driver overhead, remote diffusion subcycles, family co-location
   traffic.

## 11. Readiness for the owner's first compilation

The source is coherent and ready for the owner's first compile-and-debug pass. Mandatory
items that remain PARTIAL are listed in §2 and §10. They are known and bounded, and none of them
is silently bypassed: unsupported conditions throw. No feature here is numerically validated or
scientifically qualified.
