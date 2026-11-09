# CHUÍ — patch 53 source-completion record (supersedes the patch-52 matrix)

Date: 2026-11-09. Branch `claude/nifty-ritchie-x69law`, restarted from
`origin/main` `f9ba682` (merge of PR #132 = patch 52). The supplied
`Cosmology-Simulator-main(52).zip` is byte-identical to that tree (`diff -r`:
0 differences). Input audit: `03_AUDITORIA_PATCH52.md`.

**Executed:** source inspection, `git diff --check`, SHA-256 checks, and
data-only Python commands (download, conversion, validation; see §WP-C).
**Not executed:** CMake, compiler, linker, any CHUI test or executable, MPI
job, simulation, benchmark or sanitizer. No statement below is a runtime,
deadlock-freedom, convergence or RSS claim.

**Final label: SOURCE-COMPLETE, EXTERNAL-DATA-INCOMPLETE** (rationale at the
end). Scientific validation remains a future phase.

## WP-A (P0-52-01) — directed AMR MPI protocol

`parallel::executeBlockingDirectedAmrPatchPayloadExchange`
(`src/parallel/distributed_memory_amr.cpp`):

| Epoch | Rank-local fallible work (all guarded) | Agreement | Next communication |
|---|---|---|---|
| E0 | local patch index (admitted `parallel.directed_amr.local_patch_index`), record validation, owner check, callback presence; envelope (validates records), control-plane admission + allocation | folded into vote "directed AMR envelope preflight" | `MPI_Allgather` envelopes (return code checked) |
| E2a | envelope rank check; symmetric peer graph (each pair evaluated on the rank-ordered envelope pair, so both members compute the identical floating-point expression); remote descriptor upper bound from gathered patch counts | vote "directed AMR candidate peer graph" | — |
| E2b | staged-result admission (`…remote_patch_descriptors`) + reserve at the bound | vote "directed AMR staged descriptor admission" | first pair epoch |
| E3 per pair, ascending peer rank (= lexicographic (min,max) edge order on both sides; no cyclic wait) | pair-health handshake first: a rank that already failed keeps visiting every remaining peer and sends 0, both sides skip that pair | pairwise | descriptor count |
| E3 (a) | descriptor receive admission + resize (`exchangePodRecordsWithPeer`) | pairwise readiness | descriptor payload rounds (int-safe, bounded) |
| E3 (b) | descriptor validation, owner check, peer-count bound, planning admission, request planning both directions, request/patch indices, coverage counters, staging append within reserved capacity, traffic arithmetic | pairwise "plan ready" | symmetric skip if both request sets are empty, else cell stream |
| E3 (c) | limit narrowing, transfer plans, caller admission callback, chunk admission (only when the caller did not size the window) + chunk capacity, diagnostics | pairwise admission readiness | per round: producer callback + count check → producer readiness → payload → consumer callback → consumer readiness |
| E3 (d) | coverage checks after the pair's last message | marks this rank failed (announced at its next pair entries) | next pair |
| E4 | — | vote "directed AMR patch exchange" | result (with lease) published only on unanimous success; otherwise every rank throws and the staging is discarded |

Call-path proof sketch: inside an epoch every throwing statement either
precedes a pairwise readiness message (both ranks then stop at the same
message) or follows the epoch's last message; between epochs only the health
handshake runs; a failed rank never skips an epoch its peer will enter.
Remaining exception: an MPI error return is communicator-fatal (default MPI
error handler) and is not claimed recoverable.

`executeGovernedAmrFluxRegisterPayloadExchange` (new; the legacy vector-return
wrapper is kept) follows the same pattern: guarded validation + count control
plane → vote → `MPI_Alltoall` → guarded active-peer/inbound admission → vote →
pair-health epochs with a staging readiness message → final vote. The result
carries its inbound lease.

Callers (pre-exchange descriptor builds and post-exchange decodes are now
guarded + voted): `HydroAmrRuntimeImpl::exchangeDirectedHydroGhosts` (also
reservations declared before storage), the non-recursive distributed branch
(remote neighbour index, reflux packing, remote reflux decode/apply),
`DistributedRecursiveHooks` register exchange (governed result),
`gatherRemoteDiffusionCells`, `importRemoteStarFormationNeighbors` (now
admitted per peer), `exchangeAndValidateAmrPatchPayloads`.
Future test source (not run): `tests/integration/test_directed_amr_exchange_fault_injection_mpi.cpp`
(3 ranks: baseline, consumer/producer/admission failure, unilateral governor
refusal, one-record rounds, periodic topology, flux stale metadata; the ctest
TIMEOUT is the deadlock oracle).

## WP-B (P1-52-02) and M01 ownership

| Buffer | Admission (owner label) | Lifetime |
|---|---|---|
| association routing outbound/destinations/kept | `sources.associations.routing` (reserved at local size) | routing → merge |
| owner_of / id_to_leader / local_leader | `sources.associations.routing_scratch` | routing scope |
| kept regrowth, members, merged | `sources.associations.merge` (old+new capacity) | merge scope; merged supersedes the caller's larger staged_budgets |
| received member records | exchange lease (not charged again) | until appended |
| feedback candidate index (sorted rows/x, row-slot) | `physics.feedback.candidate_index` (old/new coexistence, strong exception safety) | deposition object |
| remote candidate storage | exchange lease adopted by the deposition object | deposition object |
| staged per-cell deltas | `physics.feedback.staged_rows` (sources × k, + inbound deltas) | deposition object |
| SSP stage per-star staging | `sources.stellar_population.stage` | stage |

M01 sweep (source-derived; not measured): canonical state, AMR geometry/
histories/flux registers, hydro batches, remote ghosts, composition exchanges,
directed AMR exchange, rollback state copy (`time.global_step_rollback.state_copy`),
treePM/PM owned fields and communication arenas (`gravity.estimate.*`),
restart/snapshot streaming, chemistry adaptive rows and the node pool are
admitted with named owners. **Estimate-only by design (M01):** FFTW plans and
internal buffers, HDF5 library chunk caches, MPI library internal buffers
(`external_runtime_reserve_bytes`). **Known not separately admitted:** the
`std::set` of remote cell keys in the migration validation exchange
(diagnostic path, bounded by received boundary cells) and the transient
`localCandidatesWithin` result vectors (≤ k entries each).

## WP-C — data products

All four expected SHA-256 values matched the supplied bytes. Data-only
reproduction: `fetch_and_convert.sh` regenerated S02 `e2161646…`, F02
`ba232a1c…`, UVB `4bfc5d9b…` and **C01 `4bbca8a9…`** byte-identically from
pinned upstream sources (WSS09 tarball `228ec928…` downloaded in this
session). C02 `633a66ee…` regenerates byte-identically from its own package;
its 10 962 cooling values match an independent Gnat & Ferland 2012 mirror
within 0.0005 dex; its CI/RR/DR coefficients are **NOT INDEPENDENTLY
SOURCE-VERIFIED** (primary archive blocked). No data file is committed.
Fact fields per product: `tools/chui_data/dataset_registry.json`; staging:
`docs/chui_data_staging.md`.

Source integration added in this patch: per-module collective agreement of
data-load/checksum failures (`RuntimeModuleRegistry::freezeAndInstantiate`);
adaptive + ladder requires an element table (config + operator); wind-table
limitations retained; one-time `sources.data_assumptions` record. Verified
existing seams: per-rank SHA-256 enforcement for every table (config requires
a 64-hex hash whenever a path is set), unverified-C02 opt-in gate, energy
contract gate, photo-rate/UVB coherence gate, ladder exclusion from the C01
query, C/O fine-structure exclusion under ladders, `-inf` as exact zero in
interpolation, HM01/FG20 mismatch refusal by default, restart compatibility of
all data keys through the normalized configuration.

C02 labels: **C02_DATA_ASSEMBLED_UNQUALIFIED, C02_LOADER_INTEGRATED**.
Not supported: C02_SOURCE_CROSS_CHECKED (only the cooling channel is
cross-checked), C02_SCIENCE_QUALIFIED, C02_PRODUCTION_READY.
C01: bytes present and reproducible (user-local), loader integrated; **no
HM01-consistent UVB** exists, so a consistent C01 + UVB model is
EXTERNAL-DATA-BLOCKED (only the explicit mismatch opt-in runs it).

## WP-D — T02 chronology

`docs/time_integration_event_chronology.md`. Default `0` keeps its meaning
("extra bound disabled"); the ledger now records the policy label, earliest
and latest physical event times, unforecast timed events and catch-up events,
and states that deposition happens at the publication epoch. The forecast is
MIN-reduced across ranks before the clip; the publisher enforces the bound for
forecast-covered records. Face-identity ledger and phase-tagged handshake:
re-read, no new source gap found.

## H01–V02/X01 evidence matrix

Columns: SG = source gap, ED = external-data gap, CQ = conditional
qualification, FV = future validation. "Fixed" = a previous finding fixed in
source in this patch.

| ID | Obligation (modes) | Production entry → state → MPI/AMR → rollback/restart/provenance | Source evidence | SG | ED | CQ | FV | Fixed |
|---|---|---|---|---|---|---|---|---|
| H01 | MUSCL/CTU/PPM, HLLC→HLLE (gas modes) | hydro stage → `HydroCoreSolver` (`src/hydro/hydro_core_solver.cpp`) → AMR sweeps / ghost exchange → step rollback copy | unchanged | no | — | — | convergence, positivity | — |
| H02 | dual energy, EOS/CFL consistency | `hydro_core_solver.cpp` dual-energy paths, effective EOS closure | unchanged | no | — | — | Mach sweep | — |
| C01 | hybrid network + element-table metal cooling | `ThermochemistrySolver::coolingHeating` → `MetalCoolingTable::evaluate` → per-rank load + SHA → normalized config | loader, UVB-mismatch gate, ladder exclusion | no | HM01-consistent UVB absent; WSS09 not redistributable | rates unqualified (opt-in) | cooling histories | collective load failure (this patch) |
| C02 | detailed reference incl. metal ions | `MetalIonRateTable::parse/loadFromFile`, `MetalIonLadderSolver`, `ThermochemistryOperator::advanceRow` → ion sidecar lanes → composition exchange/migration/snapshot | grammar, opt-in gate, energy contract, adaptive guard (this patch), validator | no | CI/RR/DR primary-source check; rights | research only | 58×3 coverage, subsets, restart, ranks | adaptive-row cooling loss (this patch) |
| C03 | adaptive fidelity | `decideChemistryFidelity`, hybrid solver for demoted rows | ladder + table requirement added | no | — | research | switching error | yes (guard) |
| C04 | sampled shadow | `ThermochemistryOperator` shadow, `updateGlobalChemistryControl` | unchanged | no | — | research | shadow thresholds | — |
| I01 | effective / resolved ISM | `EffectiveIsmThermodynamicClosure`, regime config | unchanged | no | — | — | no double cooling | — |
| I02 | SGS reservoir | SGS shear on level phases, face ledger | unchanged (patch 52) | no | — | research | energy transfer | — |
| S01 | conservative stochastic births | `StarFormation::applyFromInputs` at the interval end, age 0 | unchanged | no | — | — | rank independence | — |
| S02 | yields, event chronology | `SspYieldTable`, `executePart2StellarPopulation`, population sidecar | staging admission (this patch) | no | rights review | yield-set qualification | restart chronology | — |
| F01 | mechanical feedback, closure refresh | `MechanicalFeedbackDeposition` stage/commit, halo + delta exchanges | index/staging governance, strong exception safety (this patch) | no | — | — | conservation, MPI | yes (P1-52-02) |
| F02 | winds / pre-SN | `StellarWindTable::interval` | limitations retained + reported | no | Chabrier/multi-Z wind data; rights | IMF mismatch decision | budgets | — |
| D01 | elements, transport, diffusion | composition exchange, `DistributedElementDiffusionCoupling`, directed AMR gather | directed exchange protocol (this patch) | no | — | — | patch/level/MPI transport | yes (P0-52-01) |
| A01 | production regrid | `runAutomaticAmrRegrid` | unchanged | no | — | — | evolving AMR | — |
| A02 | distributed recursive subcycling | `advanceProductionAmrHydroRecursive`, `DistributedRecursiveHooks` | governed flux result; ghost decode vote (this patch) | no | — | — | recursive vs synchronized | yes (P0-52-01) |
| A03 | physics-aware indicators | regrid planner incl. cooling times | unchanged | no | — | — | resolved error | — |
| A04 | predictive refinement | `forecastStellarEvent` | unchanged | no | — | research | — | — |
| G01 | TreePM/KDK, fine-gas force | `src/gravity/tree_pm_coupling.cpp`, `GravityRuntime` | unchanged | no | — | **conditional (v10)** | force accuracy | — |
| T01 | one time authority, Strang, events | `TimeCoordinator` global loop, `SourceStagePlan` | event clip ordered with outputs (patch 52) | no | — | — | temporal convergence | — |
| T02 | multirate coupling, event semantics | level clocks, face ledger, event ledger | policy record (this patch) | no | — | lag value is a scientific choice | event timing, face conservation | P1-52-03 made explicit |
| M01 | governed memory, §7.2 collective safety | governor + node pool; all exchanges above | directed AMR + flux protocol, routing/index/staging admission, module-construction vote | no known unprotected rank-local hazard in the directed AMR exchange and its callers | — | library internals estimate-only | RSS/PSS, low-memory fault injection | yes (P0-52-01, P1-52-02) |
| M02 | CPU scalar/OpenMP/MPI | unchanged | — | no | — | — | scaling | — |
| V01 | ledgers, restart, provenance | normalized config (paths + SHA), `sources.data_assumptions`, event ledger | this patch | no | — | — | restart equivalence | — |
| V02 | error/resource coordination | unchanged | — | no | — | research | — | — |
| X01 | out of scope | no MHD/GRMHD, moving mesh, generalized RT, new AGN, GPU requirement | — | — | — | — | — | — |

## Open items

* External data: HM01-consistent H/He UVB; primary-source verification of the
  C02 CI/RR/DR coefficients and NIST potentials; low-T DR resonances; metal
  photo-rates from a spectral UVB; redistribution rights for every product;
  Chabrier-IMF / multi-Z wind tables.
* Test-source defect observed (pre-existing, not run):
  `tests/unit/test_baryonic_source_stage_rollback.cpp::testCoolingTimePromotion`
  builds an adaptive operator with `fixed_hybrid`, which the operator refuses
  ("adaptive chemistry fidelity requires the fixed_detailed backend"); it will
  throw when executed and needs a detailed-layout fixture.
* Validation campaign: everything in `docs/validation_plan.md` (patch 52/53
  sections).

## Label rationale

Every mandatory production-source contract named by the patch-52 audit and
found in this sweep is implemented and wired (directed AMR/flux MPI protocol,
governed routing/index/staging, data loaders, gates and provenance). The
selected full model is nevertheless missing required external bytes (an
HM01-consistent UVB for C01; primary-source-verified C02 coefficients), and
the available data are unqualified and not rights-cleared. Hence
**SOURCE-COMPLETE, EXTERNAL-DATA-INCOMPLETE** — not DATA-READY, and not a
scientific-qualification claim.

---

# CHUÍ — patch 52 source-completion record (supersedes the patch-51 matrix)

Date: 2026-11-08. Branch `claude/nifty-ritchie-x69law` (base: `180b5bb`, the
merge of PR #131 carrying the patch-51 commits the audit describes; the
`main(51)` zip itself was not supplied, so byte identity was not verified). Input:
`CHUI_patch51_auditoria_independente_fonte_2026-11-08.md`. Static edits and
static cross-file tracing only. **Nothing was compiled, configured, linked,
tested, benchmarked, simulated or validated; only `git diff --check` ran.** No
statement below is a runtime, deadlock-freedom, convergence, positivity,
RSS/PSS or machine-fit claim.

**Status: NOT SOURCE-COMPLETE.** Work packages A–D are implemented in source
(unvalidated). One codebase-wide M01 §7.2 property is still open in source:
`parallel::executeBlockingDirectedAmrPatchPayloadExchange`
(`src/parallel/distributed_memory_amr.cpp`) can throw rank-locally between its
point-to-point peer messages (record validation, owner checks, request
planning, consumer callbacks), so a local failure there can still leave a peer
blocked. It is a different primitive from the destination exchange repaired
here; closing it needs a guarded per-peer protocol (see §Open source items).

## WP-A (P1-51-01) — collective-safe destination exchange

`parallel::exchangeDestinationRecordsInto` (`src/parallel/composition_exchange.cpp`)
replaces the anonymous `exchangeIntoDestination`; all front ends route through it.

| Step | Rank-local fallible operations (guarded) | Failure agreement | Next collective |
|---|---|---|---|
| 1 | argument/allocator/destination validation; control-plane admission (`kCommunication`, counts/offsets/displacements/cursor, `(ranks+1)·(2·8+3·sizeof(size_t)+4·4)` B) and allocation; per-destination tally | `rethrowCollectivePreparationFailure("destination record exchange preparation")` | `MPI_Alltoall` counts (no throwing statement in between) |
| 4 | none: saturating round plan (`noexcept` arithmetic) | — | `allreduceMaxUint64(local_rounds)` |
| 5 | MPI-int range of the maximal round (`overflow_error`), receive lease + scratch (bucket + both maximal stagings) admission, checked offsets, receive-buffer allocation (caller allocator), bucket fill, staging `reserve` | `rethrowCollectivePreparationFailure("destination record exchange admission")` | `MPI_Alltoallv` rounds |
| 7 | none: `resize` within reserved capacity, `memcpy`, pre-checked `int` narrowing | — | next round |

Front-end phases: `exchangeCompositionDonorRequests` (destination list admitted
and allocated inside its validation `try`, vote "composition donor request
validation"), `exchangeCompositionLaneChunks` (vote "composition lane chunk
validation"), then the table above. Serial/one-rank path: no collective, local
throws propagate. MPI error returns inside steps 3/4/7 remain communicator-fatal
(default handler), as before. Zero-work ranks traverse every collective with
empty streams.

Production callers made collective-safe (decode/merge after a collective is a
guarded phase voted before the next collective):

| Caller | Votes added |
|---|---|
| `DistributedElementDiffusionCoupling::setup` (claim/conductance) | typed governed receive replaces byte receive + unguarded typed copy that preceded "interface validation" |
| `SourceRuntimeImpl::importRemoteFeedbackCandidates` | "remote feedback halo queries", "… responses", "… import" |
| `SourceRuntimeImpl::mergeAssociationsDistributed` | decode now guarded (reported at "association discovery"), "association member routing", "association member merge" |
| `executePart2StellarPopulation` feedback deltas | "mechanical feedback remote delta export", "… import" |
| `executePart2StellarPopulation` prologue | "stellar population stage prologue" |
| `importRemoteStarFormationNeighbors` | "star-formation remote neighbour import" |

## WP-B (P1-51-02) — ownership and coexistence (source-derived, not measured)

`GovernedRecords`/`GovernedRecordBuffer` now declare `lease` before the storage
(implicit destruction frees storage first), move assignment frees the old
storage before releasing the old lease, `releaseStorage()` does both in that
order; copies are deleted. The primitive charges its own memory;
`compositionLaneExchangeAdmissionBytes` now covers caller-owned buffers only
(outbound records + destinations + reassembled lanes + sort order), removing the
patch-51 double charge of received records and bucket.

| Stream | Caller-owned (reservation, owner) | Primitive-owned transient | Returned | Coexistence / release |
|---|---|---|---|---|
| donor requests | requests (small, ungoverned control), destination list (`parallel.composition_donor_requests.destinations`) | control, bucket, staging | records + lease | destinations die at return; received requests live through response build |
| donor fractions | `response_reservation` before `responses`/destinations | same | chunk records + lease | outbound freed and lease released right after the exchange; assembly lease scoped before `fractions` |
| corrections | `correction_reservation` before `corrections`/destinations | same | chunk records + lease | inbound records + assembled corrections until commit |
| diffusion claims / conductances | `m_reservation` (now declared before the tables): tables, outbound claim/reply + destinations, export/dependent rows | same | typed records + lease (`…claim_exchange`, `…conductance_exchange`) | inbound vectors live to the end of `setup` |
| diffusion state / delta chunks | `m_reservation`: outbound + reassembled + accumulation scratch (slot/id/rank/flat delta) per subcycle | same | chunk records + lease | per subcycle |
| feedback queries | `query_lease` (bounded active sources × (ranks−1)) | same | queries + lease | released after the responses are built |
| feedback candidates | `response_lease` (≤ k per query) | same | candidates + lease | storage moved into `deposition`; lease returned to `remote_candidate_lease`, declared before `deposition` |
| association members | outbound/kept (covered by `global_reservation` scale only — see open items) | same | records + lease | appended to `kept`, lease dies at scope end |
| feedback deltas | `deposition.remoteDeltas` (owned by deposition scratch) | same | records + lease | merged before the vote, then freed |
| interface face ledger | `sources.amr_interface_face_ledger` per level (`m_interface_face_ledger_leases`) | — | — | pre half → post-children closure; dropped at interval open/close and rollback |

External MPI/FFTW/HDF5 internals remain M01 estimates. The step-rollback state
copy, HDF5 chunk buffers and PM meshes are separate phases admitted by their
own owners (unchanged); no new overlap with them was introduced.

## WP-C (T02/S01/S02/F01/F02/T01) — physical event time vs publication epoch

| Event family | Physical-time meaning | Time authority | Publication epoch | Conserved transfer | Rollback / restart owner |
|---|---|---|---|---|---|
| Star birth (S01) | model: birth at the interval end (probability integrated over the interval) | plan `event_time_code` = coarse end | same instant (age 0, F5; `stellarElapsedYearsInInterval`) | gas → star mass/momentum/elements | population record created at birth (`chronology_time_code` = birth); sidecar restart |
| Continuous SSP return (AGB/CC mass, elements, dust), table winds/pre-SN | interval-integrated over absolute age `[age_b, age_e]`; no intra-step instant exists in the model | accepted interval | interval end | ejecta from star, exact element/dust sums | `stellar_age_years_last`, record cursor (`advanceStellarPopulationChronology`, gap/replay refusal) |
| Stochastic CC/Ia SN | **exact** model age `a_k = inf{a: Λ_c(a) ≥ Λ_k}` (`physics::stellarEventAgeYears`, `realizedEventChronology`) | event forecast `SourceRuntime::discreteEventStepEndLimitTimeCode` → MIN over ranks → step clip when `physics.ssp_max_event_publication_lag_years > 0` | interval end; lag `age_e − a_k` recorded (`DiscreteEventChronologyLedger`), enforced fail-closed when bounded | SN energy/momentum/ejecta via mechanical deposition | thresholds/ordinals in the record (restart); ledger in the rollback checkpoint; forecast is a pure function of checkpointed state |
| Expectation-mode SN | fractional expected count, interval-integrated (no discrete time) | interval | interval end | same | same |
| BH seeding/accretion/feedback | seeding evaluated on the synchronized end state; accretion/feedback rate × full interval | interval | interval end | BH mass/energy budgets | BH sidecar (unchanged) |

Trace of an SN whose threshold lies inside a coarse interval: the record holds
`Λ_k` (pre-drawn); at step start every rank evaluates
`nextStellarEventAgeYears` for its stars, `time_coordinator.cpp` (global step
loop) reduces MIN and clips `ordered_step_limit_time_code` to
`t_k + L` (event `time.source_event_publication_clip`, ordered with output
events); the source post half of that interval realizes the event
(`realizeEvents`), computes `a_k`, records the lag and refuses lag > L·(1+10⁻⁶)
for records covered by the forecast; deposition is at the synchronized end, with
no pending reflux/history/split (patch-51 epoch guard). With `L = 0` (default)
behaviour is unchanged and the lag (≤ coarse dt) is reported. The criterion never
shortens a step below `L`. Legacy records catching up past chronology are counted
(`catch_up_discrete_events`), not timed. The event's identity, ordinal, age and
budget do not depend on refinement, rank, retry or restart.

Limitation (explicit): deposition remains first-order in time with respect to the
physical event time unless `L` is chosen small; no claim of exact-time impulsive
deposition is made. Not numerically qualified.

Also fixed: `executePart2StellarPopulation` referenced `level_selection`, an
identifier not in scope (introduced by `863f56a`); discrete events use the
accepted coarse step index.

## WP-D (T02/I02/D01) — exact coarse/fine face coverage

Each rank builds, in `restrictDiffusionFacesToLevelSelection`, the sorted
`CanonicalInterfaceFace` set (lo/hi stable `gas_cell_id`, axis, side relative to
lo, periodic flag, remote flag) of the coarse/fine faces it integrates. The
level pre half parks the set in `SourceIntervalAccumulator::pending_interface_faces`
(governed lease, vote "recursive AMR interface face ledger open"); the
post-children closure compares its own set element-for-element (vote "… close")
before the existing global count check. Because hierarchy, leaf identity and
rank ownership are frozen inside the coarse interval, per-rank set equality on
all ranks is equality of the global face multiset: a missing face balanced by a
duplicate is detected, and duplicates inside one phase are refused. Mirrored MPI
faces are matched key-for-key across owners by the claim handshake, which now
also carries the T02 phase tag (`ClaimRecord::phase_check`). Ledgers are dropped
at interval open/close and step rollback (`resetInterfaceFaceLedger`); none
survives a legal restart boundary. Chronology (unchanged from patch 51): pre
half at `t_b` with finer levels synchronized; closure after all children,
complete reflux and level-history retirement; each coarse/fine face integrates
`dt_L/2 + dt_L/2`, same-level faces `dt/2 + dt/2`.

## H01–V02/X01 matrix (patch 52)

| ID | Classification | Production evidence |
|---|---|---|
| H01 | source implemented (unvalidated) | `src/hydro/hydro_core_solver.cpp` (unchanged) |
| H02 | source implemented (unvalidated) | dual energy / EOS refresh (unchanged) |
| C01 | external production-data dependency | loaders + refusals; tables absent |
| C02 | source implemented (unvalidated) + external data | ion ladder; charge exchange a declared limitation |
| C03 | source implemented (unvalidated, research) | `fidelityModelForRow`, hybrid shadow |
| C04 | source implemented (unvalidated, research) | `ThermochemistryOperator::advanceRow` shadow |
| I01 | source implemented (unvalidated) | unchanged |
| I02 | source implemented (unvalidated, research) | SGS shear on level phases, exact face ledger |
| S01 | source implemented (unvalidated) | endpoint birth, age 0 |
| S02 | source implemented (unvalidated) + external data | exact event ages, publication-lag ledger/bound; yields absent |
| F01 | source implemented (unvalidated) | mechanical deposition, collective-safe halo/deltas |
| F02 | external production-data dependency | wind table absent; loader present |
| D01 | source implemented (unvalidated) | collective-safe governed exchange, exact coverage |
| A01 | source implemented (unvalidated) | unchanged |
| A02 | source implemented (unvalidated) | post-children phase + face ledger |
| A03 | source implemented (unvalidated) | unchanged |
| A04 | source implemented (unvalidated, research) | `forecastStellarEvent` (unchanged) |
| G01 | conditional pending qualification | by specification |
| T01 | source implemented (unvalidated) | Strang plan, exactly-once events, event clock |
| T02 | source implemented (unvalidated) | explicit physical event time + time-authority clip; exact face coverage |
| M01 | **source incomplete** | destination-exchange streams closed; directed AMR peer exchange still has rank-local throw windows; library internals estimate-only |
| M02 | source implemented (unvalidated) | unchanged |
| V01 | source implemented (unvalidated) | ledgers, restart compat for the new key |
| V02 | source implemented (unvalidated, research) | unchanged |
| X01 | out of scope | no MHD/GRMHD, moving mesh, CHIMES/Grackle, new AGN, GPU requirement |

Config/schema/reproducibility: new key `physics.ssp_max_event_publication_lag_years`
(default 0; typed validation, normalized dump, restart projection accepts older
checkpoints only at the default). No snapshot/restart field change. Default runs:
unchanged physics and step sequence (the new ledger is diagnostic). Runs with
`L > 0` take shorter steps (deterministic, rank-independent).

Tests written but NOT RUN: `tests/unit/test_composition_lane_chunk_wire.cpp`
(`testGovernedLeaseLifetime`), `tests/unit/test_stellar_population_feedback.cpp`
(`testEventPhysicalTimeIsSubdivisionInvariant`).

## Public interface changes and migration (`include/cosmosim/**`)

* `parallel::GovernedRecords<T>` / `GovernedRecordBuffer`: member order is now
  `lease, records|bytes`; types are move-only (copy deleted) and no longer
  aggregates. Brace-initialization by position must be replaced by member
  assignment; `.records`/`.bytes`/`.lease` access is unchanged. New
  `releaseStorage()`.
* New `parallel::exchangeDestinationRecordsInto`, `DestinationReceiveAllocator`
  and template `parallel::exchangeGovernedRecords<T>`; prefer them over
  `executeBlockingDestinationRecordExchange` (kept, ungoverned, same protocol).
* `parallel::compositionLaneExchangeAdmissionBytes` now returns caller-owned
  bytes only (smaller than before); callers that also counted the primitive's
  receive buffer/bucket must not add them back.
* `workflows::SourceRuntime::discreteEventStepEndLimitTimeCode` (virtual,
  default +infinity: other implementations need no change).
* `physics::{stellarEventCumulativeCount, stellarEventChannelHorizonYears,
  stellarEventAgeYears, RealizedEventChronology, realizedEventChronology,
  nextStellarEventAgeYears}` (additive).
* `core::PhysicsConfig::ssp_max_event_publication_lag_years` (additive, default 0).

## Open source items

1. `executeBlockingDirectedAmrPatchPayloadExchange`: guarded per-peer protocol
   (validate/plan/allocate per peer, then agree before the next peer message,
   or a two-phase collective plan) — not done in patch 52.
2. Association member routing scratch (`owner_of`, `local_leader`, outbound) is
   allocated in a guarded phase but not separately admitted (bounded by the
   admitted global descriptor set).
3. Feedback-candidate index growth inside `setRemoteCandidates` is not
   separately admitted.

Later validation (not source gaps): build matrix; multi-rank low-memory fault
injection on every vote above; event timestamp / lag-bound tests; face-ledger
mismatch injection and recursive-vs-synchronized conservation; restart
equivalence with `L > 0`; RSS/PSS reconciliation; four-mode workflows.

---

# CHUÍ — patch 51 source-completion record (supersedes the patch-50 matrix)

Date: 2026-11-08. Branch `claude/intelligent-cray-m8rxdv`. Input audit:
`CHUI_patch50_auditoria_independente_fonte_2026-11-08.md`. Static edits and
static cross-file tracing only. **Nothing was compiled, configured, linked,
tested, benchmarked, simulated or validated.** No statement below is a runtime,
convergence, MPI-deadlock-freedom, performance or memory-fit claim.

**Status: source work packages A–D implemented; NOT SOURCE-COMPLETE only in the
sense that T02 discrete events follow the specification-permitted
synchronized-epoch exception (owner confirmation requested, see below).**
External datasets (C01/C02/S02/F02) and G01 qualification remain separate.

## Audit findings → source changes

| Finding | Resolution (path::symbol) | Guaranteed by source / left to qualification |
|---|---|---|
| P0/P1-50-A coarse/fine time | `amr::RecursiveAmrSourcePhase::kPostChildrenSynchronization` (called by `advanceProductionAmrHydroRecursiveImpl` after children, reflux and history retirement, every rank); `RecursiveAmrLevelSourceAdapter::apply`; `SourceRuntime::executeAmrLevelInterfaceSynchronization`; `SourceRuntimeImpl::{restrictDiffusionFacesToLevelSelection, levelPhaseTag, executeAmrLevelSmoothSources}`; phase-tagged `gatherRemoteDiffusionCells` and `DistributedElementDiffusionCoupling(…, phase_tag)` | Source: every face owned by its COARSER side L; same-level faces dt/2 in L pre + dt/2 in L post; coarse/fine faces dt/2 in L pre (all finer levels at t_b, before L history start) + dt/2 at the post-children synchronization (all levels ≥ L at t_e, after reflux/history retirement). Both sides read at the same instant, same interval, equal/opposite transfers applied once; pre/sync pairing verified by global coarse/fine face counts; open splits refused at interval close and before events; remote records of another phase refused (phase-tagged sequences). Qualification: temporal order of the S(dt/2)[subtree]S(dt/2) split, positivity margin, recursive vs synchronized equivalence. |
| A — composition ordering risk | `composition_transport.hpp` contract | Source/argument: donor fractions are read at prepare time (post-source); diffusion is mass-neutral per cell and maximum-principle bounded, so outflow F·X ≤ M_s iff F ≤ M (hydro-only condition) — diffusion adds no negativity risk. Deferred-donor CFL remains checked exactly and refused collectively. |
| P1-50-B hybrid shadow / cooling authority | `ThermochemistryOperator::advanceRow` shadow (private detailed ladder reconstruction for demoted rows, reference never ladder-less), `fidelityModelForRow` used by `decideChemistryFidelity` and `coolingTimeSeconds` | Source: reduced rows stay reduced; reference measures demotion error (stage fractions, metal electrons / reference total n_e, T, stored potential); published once per interval through the existing C03/V02 path (rollback-restored). Cooling estimate uses the row's own model: no metal cooling vanishes, no double counting. Qualification: thresholds. |
| P1-50-C MPI lease / allocation | `parallel::{GovernedRecords, GovernedRecordBuffer, executeGovernedDestinationRecordExchange}`, internal `exchangeIntoDestination`; callers in `completeDistributedCompositionTransport`, `DistributedElementDiffusionCoupling::addRemoteFaceTransfers` | Source: receive directly into the final typed buffer (no decode copy); its lease is returned with it; staging + bucket admitted separately and pre-reserved at the maximal round; no allocation between the admission vote and the payload rounds. An MPI error return is still rank-local and fatal for the communicator (default MPI error handler) — not claimed recoverable. |
| P1-50-D discrete events | `SourceRuntimeImpl::executeDiscreteSourceEvents` epoch guard + `SourceIntervalAccumulator::discrete_event_dispatches` | Specification-permitted exception (route v10 §4.6 / T02 “stochastic-event separation”, “do not force stochastic events into smooth integrators”, F5 absolute birth time): events run once per interval at the absolute interval end, the only legal publication epoch (hierarchy synchronized, verified collectively). Budgets are age-integrated (cumulative cursors, gap/replay refusal), invariant under level cadence/refinement/rank/retry; unresolved deposits carried. **Owner decision requested** if per-substep event dispatch is required instead. |
| P2-50-E data / charge exchange | unchanged loaders and refusals | Metal–H/He charge exchange is not named by C02 in route v10 (“selected metal ions/molecules”, qualification “exact network inventory, reaction data”); it is recorded as a model limitation, not implemented, and no rate is invented. |

## H01–V02 source matrix (patch 51)

| ID | Status | Production evidence / exact remaining gap |
|---|---|---|
| H01 | implemented in source | `hydro_core_solver.cpp` (unchanged) |
| H02 | implemented in source | dual energy / EOS refresh (unchanged) |
| C01 | external-data dependency | loaders + refusals; production tables absent |
| C02 | implemented in source + external-data dependency | ladder physics (patch 50); charge exchange a declared limitation; verified dataset with energy contract absent |
| C03 | implemented in source (research) | hybrid-row ion error now measured; `fidelityModelForRow` |
| C04 | implemented in source (research) | detailed reference for every sampled row, stack-only |
| I01 | implemented in source | unchanged |
| I02 | implemented in source (research) | SGS shear on causal level phases |
| S01 | implemented in source | synchronized-epoch events (exception above) |
| S02 | external-data dependency | age-integrated returns; yields absent |
| F01 | implemented in source | synchronized-epoch deposition |
| F02 | external-data dependency | wind dataset absent |
| D01 | implemented in source | schema-2 wire, governed leases, causal multilevel diffusion |
| A01 | implemented in source | unchanged |
| A02 | implemented in source | new post-children phase in recursive driver |
| A03 | implemented in source | cooling time from the row's fidelity model |
| A04 | implemented in source (research) | unchanged |
| G01 | conditional | qualification pending by specification |
| T01 | implemented in source | Strang plan, exactly-once events |
| T02 | implemented in source with specification-permitted event exception | smooth sources/diffusion/SGS causal per level; events at synchronized epochs |
| M01 | implemented in source for traced owners | MPI receive leases; FFTW/HDF5/MPI-library internals estimate-only |
| M02 | implemented in source | unchanged |
| V01 | implemented in source | ledgers/restart (no schema change this patch) |
| V02 | implemented in source (research) | hybrid ion disagreement reaches the joint controller |

X01 out of scope (no MHD, GRMHD, moving mesh, SPH, CHIMES, Grackle, new AGN).

Restart/schema impact: none (all new state is runtime controller state inside
`SourceIntervalAccumulator`, captured by the existing rollback checkpoint and
always closed at legal restart boundaries). Reproducibility impact: recursive
AMR runs with diffusion/SGS change (coarse/fine faces now integrated on the
coarser clock); adaptive runs with shadow sampling change their promotion
bias (hybrid ion disagreement).

Tests written/updated but NOT RUN: `tests/unit/test_amr_recursive_subcycling.cpp`
(post-children phase counts), `tests/unit/test_composition_lane_chunk_wire.cpp`
(governed return type).

---

# CHUÍ — patch 50 source-completion record (supersedes the patch-49 matrix)

Date: 2026-11-08. Branch `claude/intelligent-cray-m8rxdv`. Scope: static source
edits and static cross-file tracing only.
**No compilation, configuration, tests, benchmarks, simulations, sanitizers or
scientific validation were performed.** Nothing below is COMPILED, TESTED,
SCIENTIFICALLY VALIDATED or PERFORMANCE VERIFIED.

**Status: NOT SOURCE-COMPLETE.** Remaining mandatory source gap: T02 discrete
events (SF/SSP/winds/SN/BH) are still a coarse-barrier synchronization
operator, not dispatched per AMR level substep (exact entry point below).
External datasets (C01/C02/S02/F02) and G01 qualification remain separate.

## Audit findings closed in source

| ID | Fix (path::symbol) | Contract chosen |
|---|---|---|
| P0-49-01 | `parallel/composition_exchange.{hpp,cpp}`: `CompositionLaneChunkRecord`, `appendCompositionLaneChunks`, `assembleCompositionLaneChunks`, `exchangeCompositionLaneChunks`, bounded-round `executeBlockingDestinationRecordExchange`; callers `HydroAmrRuntimeImpl::completeDistributedCompositionTransport`, `DistributedElementDiffusionCoupling::{setup,addRemoteFaceTransfers}`; bound `physics::k_max_constructible_composition_lane_count` checked in `CompositionLayout::make` and `static_assert` in `hydro_amr_runtime.cpp` | Wire schema 2: fixed 32-double chunks, ceil(L/32) per cell; header = schema, stream kind, layout signature, lane count, chunk index/count/range, exchange sequence, stable `gas_cell_id`, source/destination rank. Reassembly sorts by (id, source) and refuses truncation, duplicates, foreign layout, stale sequence, misaddressing. Receive + staging capacity reserved (`kCommunication`) before allocation; argument/admission failures are collective before the payload `Alltoallv`; per-peer round cap keeps MPI int counts below `INT_MAX/2`. The 32-lane constant was removed, not raised. Capability gate unchanged (`core_elements` allowed, scalar-Z refused on >1 rank). |
| P1-49-02 | `ThermochemistrySolver::metalIonEnvironment`; used by `ThermochemistryOperator::attachMetalIonState`, `integrateCoupled::advance_ions`, `ThermochemistrySolver::coolingHeating`; refusal in `ThermochemistrySolver` ctor | One radiation context. Grey approximation: HI Rahmati factor applied to every ion photo/photoheating rate (`MetalIonPhotoAttenuationModel::kHiSelfShieldingGrey`, documented as not frequency resolved). Photo-rate datasets without a configured UVB are refused. Equilibrium initialization iterates the coupled n_e (≤ 8 bounded iterations). |
| P1-49-03 | `transferLadderNuclei`, `ThermochemistryOperator::integrateCoupled`, staging in `advanceRow`, `ThermochemistryOperator::reconcileMetalIonStageLanes` (called after dust in `BaryonicSourceStage::prepare`) | Invariant: Σ_k stage nuclei = element − CO-bound nuclei after every chemistry commit. CO formation/destruction moves nuclei through the NL97 channel stages (C II, O I); remainder proportional only if the channel is exhausted; no hidden creation. `kCII` = atomic-carbon (non-CO) reservoir of the light network; ionization distribution owned by the ladder. Transport, diffusion, prolong/restrict are linear in all lanes, so the invariant is carried; dust changes are reconciled on the private row. |
| P1-49-04 | `ThermochemistrySolver::advanceCell` charge residual; closing check in `integrateCoupled`; `BaryonicSourceStage::prepare` history | Charge = light-species charges + metal-ion electrons beyond NL97 C+ (identical inventory to n_e). Ladder normalization goes to `metal_ion_normalization_residual` → element-fraction history, never charge history. |
| P1-49-05 | EOS refresh in `integrateCoupled` (μ, T, n_e/n_H from final ions at `u_end`); stage commit returns the rows to the caller's thermodynamic/CFL refresh | Conserved: u (thermal), element nuclei, CO/atomic partition. Derived: n_e, μ, T, P, c_s. No energy is manufactured to hold T. |
| P1-49-06 | `MetalIonRateTable::{hasIonizationEnergyContract, ionizationPotentialErg, storedIonizationEnergyPerNucleusErg}`, `energy_convention`/`ionization_potential` parsing; `SourceLedger::metal_ion_potential_delta_erg`; CII term moved out of `bindingEnergyDensity` when C ladder runs | Convention `collisional_ionization_work_in_cooling`: L_k contains collisional work, photoheating is excess energy; the stored-potential change is ledgered separately and not subtracted again from u. Missing contract ⇒ refusal, never zero-fill. No potentials compiled in. |
| P1-49-07 | `ThermochemistryOperator::integrateCoupled` (shared by production and reference), shadow block in `advanceRow`; `ChemistryShadowReport::{ion_samples, max_ion_stage_fraction_difference, …}`; `SourceRuntimeImpl::{accumulateSourceIntervalObservation, updateGlobalChemistryControl}`; `V02ControlState::{last_shadow_ion_error, last_shadow_ion_samples}` | Reference = same chronology with the reference solver's own ladder (10× tighter rtol, 4× substeps), on stack copies only. Global MPI-max ion disagreement enters the promotion-bias law and V02 restart state. Local sampled diagnostic, never a certified global error; no memory-pressure downgrade path added. |
| T02 (diffusion/SGS) | `SourceRuntimeImpl::{restrictDiffusionFacesToClockLevel, leafLevelOfRow, executeAmrLevelSmoothSources, executeDiffusionAndSgsShear(face_clock_level)}`; `makeSourceStagePlan` (`smooth_transport = !recursive_levels_active`); `RemoteDiffusionCell::level`, `SgsRemoteShearFace::remote_level`; check in `closeSourceInterval` | Face clock = max(level of both leaves). Each face integrated on one clock in palindromic level halves (pre: transport→chemistry; post: chemistry→transport), total integrated time = coarse dt; transfers equal/opposite and applied immediately ⇒ no coarse/fine mismatch, no reflux register needed. Rank-crossing faces use the symmetric clock so the claim/authority handshake stays one-to-one. |
| M01 (gap found) | `core::collectSimulationMemoryReport` (new `gas_composition.lanes_sparse_rows_and_history` entry); `TimeCoordinator` retry point reservation | The composition sidecar was absent from the report ⇒ governor baseline, `RetainedCapacityTransaction` baseline reconciliation and rollback state-copy admission undercounted it. Rollback now also admits both scheduler copies. |

## Production call flow — coarse level with a finer level and a rank-spanning feedback event

1. Clocks: `TimeCoordinator` agrees `accepted_dt` (MPI min) → optional retry
   point (state + composition + schedulers admitted, `m_source.captureRollbackCheckpoint`).
2. Pre-source: `StepOrchestrator` → `SourceRuntimeImpl::execute(kStrangPreHalf)`
   → `makeSourceStagePlan` (recursive levels: no global smooth chemistry or
   transport) → `openSourceInterval`.
3. Hydro: `runProductionAmrHydroPath` → `advanceDistributedProductionAmrHydroRecursive`
   → level L: `apply(kPreHydroHalf)` → `executeAmrLevelSmoothSources` →
   level-clock diffusion/SGS (faces with clock L, collective gather/handshake,
   chunked lane exchange) → thermochemistry (coupled ion chronology) →
   history start → sweep → register exchange → `apply(kPostHydroHalf)`
   (chemistry → level-clock transport) → history end.
4. Fine substeps: same per child level k = 0..r-1 with `subintervalUpdate`.
5. Reflux: `applyCompletePendingFluxRegistersToSimulationState`; temporal
   history retired; synchronization check.
6. Composition transport: `completeDistributedCompositionTransport` (schema-2
   chunks, collective validation, single commit).
7. Post-source: `execute(kStrangPostHalf)` → `executeDiscreteSourceEvents`
   (SF → SSP/winds/pre-SN/SN → BH, once, at the coarse end epoch; remote
   deposits through the existing owner-addressed exchanges) → `closeSourceInterval`
   (level-half and level-transport counters verified; one controller publication).
8. Regrid at the synchronized boundary; rollback restores state, schedulers,
   force endpoint and source checkpoint on collective numerical rejection;
   restart only at closed intervals.

## H01–V02 source matrix (patch 50)

| ID | Normative requirement | Production evidence (`path::symbol`) | Status | Exact remaining code gap | Deferred obligations |
|---|---|---|---|---|---|
| H01 | MUSCL–Hancock/CTU/PPM, HLLC+HLLE | `hydro_core_solver.cpp`, CTU ghost slabs (patch 48/49, unchanged) | implemented in source | none known | order/convergence, multi-axis positivity |
| H02 | Dual energy, EOS/comoving CFL | `hydro_core_solver.cpp`, `gas_thermodynamic_refresh.cpp` | implemented in source | none known | Mach sweep |
| C01 | Hybrid light network + element cooling | `thermochemistry.cpp`, `metal_cooling_table.cpp` | external-data dependency | none in source | production UVB/metal tables |
| C02 | Detailed reference incl. selected metal ions | `metal_ion_network.cpp`, `ThermochemistryOperator::integrateCoupled`, P1-49-02..06 above | implemented in source + external-data dependency | none known; metal–H/He charge exchange unrepresented (named omission) | verified dataset with energy contract; conservation/equilibrium tests |
| C03 | Adaptive fidelity, hysteresis | `decideChemistryFidelity`, `BaryonicSourceStage`, ion history fixes | implemented in source (research) | none known | matched-error qualification |
| C04 | Sampled shadow | `advanceRow` shadow + coupled ion reference | implemented in source (research) | none known | bounded-cost measurement |
| I01 | Effective vs resolved ISM | `effective_multiphase_ism.cpp`, `BaryonicSourceStage` | implemented in source | none known | consistency tests |
| I02 | SGS reservoir | `ism_dust_sgs.cpp`, level-clock shear (T02) | implemented in source (research) | none known | energy-transfer calibration |
| S01 | SF eligibility, conservative births | `star_formation.cpp`, `executeDiscreteSourceEvents` | implemented in source | per-level birth cadence (T02) | patch/rank independence |
| S02 | SSP returns, event timeline | `stellar_population.cpp`, cumulative age bookkeeping | external-data dependency | per-level deposition cadence (T02) | licensed yields |
| F01 | Mechanical SN deposition | `supernova_feedback.cpp` | implemented in source | per-level deposition cadence (T02) | MPI/resolution tests |
| F02 | Winds/pre-SN | `stellar_wind.cpp` | external-data dependency | per-level cadence (T02) | wind dataset |
| D01 | Elements, advection/diffusion, dust, molecules | `CompositionFluxAccumulator`, schema-2 wire, `DistributedElementDiffusionCoupling` | implemented in source | none known (P0-49-01 closed) | ≥ 37-lane multi-rank tests |
| A01 | Production regrid | `runAutomaticAmrRegrid`, `planProductionAmrRegrid` | implemented in source | none known | evolving-workflow tests |
| A02 | Distributed recursive subcycling | `advanceProductionAmrHydroRecursiveImpl`, pending registers | implemented in source | none known | nested-history MPI tests |
| A03 | Physics-aware indicators | `patchDecision`, cooling-time path | implemented in source | none known | resolved-error studies |
| A04 | Predictive/cost-aware AMR | predictive block, joint ordering | implemented in source (research) | none known | fixed-budget comparisons |
| G01 | TreePM/KDK, fine-gas force check | TreePM, AMR softening policy | conditional | correction only if qualification demands | force-accuracy qualification |
| T01 | One time authority, Strang, chronology | `SourceStagePlan`, once-only events, retry | implemented in source | none known | temporal convergence |
| T02 | Recursive multirate particle/AMR/source | level chemistry + level-clock diffusion/SGS (this patch); KDK unchanged | **partial** | `SourceRuntimeImpl::executeDiscreteSourceEvents` still runs at the coarse post half: missing per-level-substep SF/SSP/wind/SN/BH dispatch with level-timestamped births, cross-level/cross-rank deposit aggregation and deferred particle materialization at a legal barrier | mixed-source convergence, exactly-once under retry |
| M01 | Canonical state, governed peaks | ownership table below; composition now in baseline | implemented in source for traced owners; FFTW/HDF5/MPI internals estimate-only | none known in traced owners | RSS/PSS measurement, 48 GiB fit not claimed |
| M02 | CPU serial/OpenMP/MPI, optional GPU | unchanged runtimes | implemented in source | none known | rank/NUMA measurements |
| V01 | Ledgers, four modes, restart/provenance | ledgers extended (ion potential/CO/inventory); restart notes | implemented in source | none known | end-to-end four-mode tests |
| V02 | Error/resource coordination | joint controller + ion shadow keys | implemented in source (research) | none known | policy comparisons |

X01 remains outside scope (no MHD/RT/new AGN/CHIMES/Grackle added).

## M01 ownership / admission table (statically traced)

| Phase / buffer | Owner | Capacity formula | Allocation site | Admission site | Lifetime / release |
|---|---|---|---|---|---|
| Composition sidecar (dense lanes, sparse rows, history) | `GasCompositionSidecar` | `ownedCapacityBytes()` (Σ lane capacities, sparse ids/values, 10 history vectors) | configure/growth | baseline via `collectSimulationMemoryReport` (new); growth via `RetainedCapacityTransaction` | persistent |
| Global-step rollback copy | `TimeCoordinator` retry point | baseline(state) + baseline(schedulers) | `retry_point.state.emplace` | `time.global_step_rollback.state_copy` (kPhaseResident) | until accept/restore |
| Multirate PM rollback | `TimeCoordinator` | rows × per-row + metadata | rollback store | `time.multirate_cycle_rollback_snapshot` | cycle |
| Composition donor/correction exchange | `completeDistributedCompositionTransport` | `compositionLaneExchangeAdmissionBytes(n, L)` + primitive recv/staging | response vectors / primitive | `hydro.distributed_amr.composition_*` + `parallel.*` (kCommunication) | call frame |
| Diffusion interface exchange | `DistributedElementDiffusionCoupling` | tables + 2× lane-exchange bytes per face | setup / per subcycle | `sources.element_diffusion.remote` + primitive | coupling lifetime |
| Diffusion workspace | `advanceElementResolvedDiffusion` | `elementDiffusionWorkspaceBytes(L, rows)` | operator | `sources.element_diffusion.workspace` | phase |
| Diffusion topology / faces | `executeDiffusionAndSgsShear` | construction bound + retained faces | builder | `sources.metal_diffusion.*` | phase / retained graph |
| SGS shear staging | `SgsShearTransferWorkspace` | `requiredBytes(cells)` | prepare | `sources.metal_diffusion.cell_phase` | phase |
| Detailed chemistry Jacobian / ion ladder / shadow | `ThermochemistrySolver`, `MetalIonLadderSolver`, `integrateCoupled` | fixed stack arrays (≤ 11² doubles, ≤ 27 stages, two `CoupledRowOutcome`) | stack | none required (no heap) | row |
| Source-stage staging | `BaryonicSourceStage` | rows × lanes | `m_lanes.reserve` | source-stage transaction (existing) | stage |
| AMR regrid old+new hierarchy | `runAutomaticAmrRegrid` | candidate state copy | candidate | `amr.production_*.transaction` | until publish |
| AMR temporal history / registers | AMR | per-interface | ghost fill | `amr.temporal_history.*`, `amr.pending_flux_registers.*` | level interval |
| TreePM / FFTW | gravity runtime | worker scratch + declared reserve | solver | `treepm.worker_leaf_scratch`, `parallel.process_fftw_unknown_reserve_bytes` (estimate-only) | step |
| HDF5 staging | output runtime | declared reserve + overlap | library | `parallel.process_hdf5_unknown_reserve_bytes`, `process_output_restart_overlap_bytes` (estimate-only) | output |

Node aggregation: every reservation above goes through the rank governor,
which charges the attached `NodeMemoryPool` when shared-pool admission is
configured (`MemoryGovernor::attachNodePool`). FFTW/HDF5/MPI-library internals
remain opaque estimate-only reserves; no RSS/PSS or 48 GiB fit is claimed.

## Regression review (static)

H01/H02, A01–A04, G01, M02: untouched. I01/I02: SGS shear unchanged in uniform
runs; recursive runs integrate shear per level clock. S01/S02/F01/F02: event
code unchanged. D01: wire protocol replaced (all three streams + diffusion);
other users of `executeBlockingDestinationRecordExchange` keep their arguments
(two new defaulted parameters) and gain one collective argument check and one
round-count reduction per call (rank-uniform). T01: plan ordering unchanged in
uniform runs. V01/V02: additive ledger/state fields. C01: binding ledger
unchanged unless the C ladder is selected. New refusals (deliberate): ion
datasets without the energy contract; photo-rate ion datasets without a UVB.

## Tests written but not run

`tests/unit/test_composition_lane_chunk_wire.cpp` (registered as
`unit_composition_lane_chunk_wire`, NOT RUN — USER-RESERVED).

---

# Patch-49 record (history, superseded by the patch-50 matrix above)

# CHUÍ — patch 49 source-completion record (supersedes the patch-48 self-assessment)

Date: 2026-11-05. Scope: static source edits only.
**No compilation, tests, benchmarks, simulations, or scientific validation were performed in this campaign.**

**Status: NOT SOURCE-COMPLETE.** The patch-48 section below is retained as
history; its matrix was an unverified self-declaration and is superseded by the
matrix in this section (independent audit
`CHUI_patch48_auditoria_independente_fonte_2026-11-04.md`).

## Contradictions closed (P0 + additional findings)

| ID | Defect | Resolution (symbol) |
|---|---|---|
| P0-01 | stale `source_evaluation_scale_factor` / undefined `smooth_source_dt_code` | `SourceStagePlan` + `makeSourceStagePlan()`; every operator reads its clock from the plan |
| P0-02 | `executeThermochemistry` returned a report without value; caller arity | returns `void`; observations only accumulate (`accumulateSourceIntervalObservation`), publication in `closeSourceInterval()` |
| P0-03/04 | `const amr_options` mutated; undefined adapter/context types | mutable options; `RecursiveAmrLevelSourceAdapter` (hydro_amr_runtime.cpp); legacy callback API removed |
| P0-05 | `setMultilevelSourceRuntime ... override` with no base | removed; `bindSourceRuntime()` (already called by `reference_runtime_composition.cpp`) binds and calls `setAmrLevelSourceCouplingEnabled()` |
| P0-06 | CTU ghost telemetry fields undeclared | declared in `HydroProfileEvent` and `ProductionAmrHydroDiagnostics`; patch/level aggregation and global rank sums |
| P0-07 | private `invalidateOccupancyCache()` called from rollback | public `DistributedParticleIdRegistry::invalidateForRollback()` |
| new | `persist_v02_regrid_state` called, never defined; applied regrid never published V02 | defined; AMR half of V02 published on no-op and applied passes |
| new | `ProductionAmrRegridPlan::optional_admitted_*`, `optional_min_admitted_benefit` undeclared | declared |
| new | predictive AMR block: `catch` body swallowed the local forecast; `predictive_enabled` undeclared | restructured, admitted before allocation |
| new | `admitCompositionSparseActivation` called at three sites, never defined | defined (retained-capacity transaction, canonical class) |

## Production reachability (traced call flow)

* `.param.txt` → `SimulationConfig` → `reference_runtime_composition.cpp`
  builds hydro then source owners and calls `HydroAmrRuntime::bindSourceRuntime()`
  → coupling enabled iff `numerics.amr_time_stepping = recursive_subcycling`.
* `StepOrchestrator` dispatches `kSourceTerms` twice (`kStrangPreHalf`,
  `kStrangPostHalf`) in both the single-step and hierarchical-block paths →
  `SourceRuntimeImpl::execute()` → plan → smooth halves → (post only)
  `executeDiscreteSourceEvents()` → `closeSourceInterval()`.
* `runProductionAmrHydroPath()` → recursive driver
  (`advanceProductionAmrHydroRecursive[Impl]`, serial and distributed) → per level
  `RecursiveAmrLevelSourceOperator::apply(kPreHydroHalf)` before history start
  and `apply(kPostHydroHalf)` after sweep/register exchange →
  `executeAmrLevelSmoothSources()` → `executeThermochemistryInterval()`.
* Regrid cooling-length criterion: accepted chemistry history, else
  `SourceRuntime::fillAmrCoolingTimes()` (same rate set), else legacy provider.

## H01–V02 source-coverage matrix (patch 49)

| ID | Status | Evidence / exact remaining gap |
|---|---|---|
| H01 | Implemented in source, not tested | CTU ghost-slab increments + counted fallback; telemetry now declared/aggregated. Formal order/convergence: later validation |
| H02 | Implemented in source, not tested | dual-energy, floors, comoving factors (unchanged) |
| C01 | External data dependency | hybrid network + table loaders; production UVB/metal tables not supplied |
| C02 | Implemented in source, not tested + **External data dependency** | `metal_ion_network.*`: C/O/Si/Fe stage inventory, validated `chui_metal_ion_rates_v1` loader, ladder solver, layout/transport/diffusion/snapshot/restart-signature/config. No licensed, source-verified rate dataset exists in the repository; charge exchange with H/He not represented (schema v1) |
| C03 | Implemented in source, not tested | residual/error history, hysteresis/dwell (physical years), sparse mapping; ion ladders initialized from equilibrium on promotion (counted) |
| C04 | Implemented in source, not tested | shadow non-committing; controller published once per Strang interval |
| I01 | Implemented in source, not tested | regime split unchanged |
| I02 | Implemented in source, not tested | SGS shear transfer now receives the defined half duration |
| S01 | Implemented in source, not tested | births only in post half at `t_end`; global synchronized (not per level) |
| S02 | External data dependency | cursor/chronology once per interval; production yield tables absent |
| F01 | Implemented in source, not tested | feedback once per interval (post half) |
| F02 | External data dependency | winds/pre-SN once per interval; production wind dataset absent |
| D01 | Implemented in source, not tested | capability authorizes `core_elements`, refuses `total_only` on >1 rank; ion lanes diffuse with elements |
| A01 | Implemented in source, not tested | mandatory/optional budgets, transactional publication (patch 48, preserved) |
| A02 | Implemented in source, not tested | recursive distributed reflux (patch 48, preserved); level Strang ordering w.r.t. histories now explicit |
| A03 | Implemented in source, not tested | cooling-length fallback through `fillAmrCoolingTimes` |
| A04 | Implemented in source, not tested (research) | predictive block repaired |
| G01 | Conditional on later qualification | TreePM/KDK retained; no multilevel correction without evidence |
| T01 | Implemented in source, not tested | `SourceStagePlan`, once-only events, single publication, rollback of open interval |
| T02 | **Partial, exact missing code** | smooth thermochemistry/ISM/dust/SGS per level clock (pre+post). Diffusion, SGS shear transfer, SF, SSP, winds, SN, BH remain globally synchronized at coarse halves/endpoint. Missing: level-local cadence for element diffusion/SGS transfer with coarse/fine flux registers for composition-diffusion fluxes, and per-level event ownership for SF/feedback |
| M01 | **Partial — not exhaustively traced** | sparse species/ion storage, transport sparse activation (`admitCompositionSparseActivation`), element-diffusion workspace (`elementDiffusionWorkspaceBytes` over all lanes, including ion lanes) and predictive-AMR points are admitted before allocation. Not traced end-to-end in this campaign: that a single admission decision accounts for the coexisting rollback copy + PM + HDF5 + MPI buffer peak, and that `parallel/node_memory_pool` node-aggregate admission covers every phase-resident owner |
| M02 | Implemented in source, not tested | CPU serial/OpenMP/MPI unchanged; GPU optional |
| V01 | Implemented in source, not tested | ledgers/restart; four modes untouched |
| V02 | Implemented in source, not tested (research) | chemistry + AMR halves of the joint controller both publish; local estimators, not global error |
| X01 | Out of scope | no MHD/RT/new AGN |

## Remaining validation gates (owner-reserved)

Build all presets; unit tests (`test_amr_recursive_subcycling` updated for the
operator seam); Strang temporal convergence with mixed sources; exactly-once
event counts per interval under retry; recursive vs synchronized AMR source
equivalence; metal-ion ladder conservation/equilibrium against a qualified
dataset; restart round-trip with ion lanes; MPI rank-count independence.

---

# Patch-48 history (unverified self-assessment, superseded)

# CHUÍ — fechamento estático de implementação H01–V02

Data: 2026-11-04. Escopo: inspeção e edição exclusivamente estáticas do patch
47. Conforme determinação do proprietário, **não houve compilação, teste,
benchmark, simulação, sanitizer ou validação numérica** nesta campanha.

## Resultado executivo

Esta campanha fechou em fonte as contradições concretas do patch 47 e os
blocos subsequentes de integração:

1. **A02:** o refluxo recursivo distribuído deixou de ser proibido pelo sweep.
   As contribuições remotas entram no `PendingFluxRegisterStore` do proprietário
   grosso com intervalo/subpasso explícito, validação de duplicidade por origem,
   área espacial somada entre ranks no mesmo subpasso e aplicação somente após
   cobertura temporal completa. O store já é migrado e serializado no restart.
2. **T01:** a ordem canônica passou a
   `S(dt/2) -> H(dt) -> S(dt/2)`. Termoquímica, dust, SGS e difusão usam os
   intervalos físicos de meia etapa; nascimento, evolução/eventos, feedback,
   ventos e BH executam uma única vez no post-step no epoch absoluto final.
3. **D01:** a capability agora reflete o operador MPI elementar realmente
   existente e aceita apenas layout `core_elements` em build MPI. O passo de
   difusão recebe `dt/2` em cada metade Strang. O layout escalar `total_only`
   continua explicitamente recusado em múltiplos ranks.
4. **A01:** refinamentos mandatórios não consomem o orçamento opcional. O
   runtime prepara a hierarquia inteira numa cópia admitida pelo governador,
   vota coletivamente e publica por move no-throw somente após sucesso global.
5. **H01/I02:** CTU transversal em halos AMR/MPI e transferência SGS por
   cisalhamento foram ligados aos caminhos conservativos, scratch governado,
   reflux/regrid e telemetria existentes.
6. **M01:** as lanes químicas detalhadas adaptativas passaram a armazenamento
   esparso por `gas_cell_id`, com referência fixa densa, capacidade governada,
   transporte, migração, regrid, rollback e restart schema 3.
7. **A03/A04/V02:** a política força/resolução agora toma ação explícita; o
   forecast de eventos cruza ranks; e um controlador persistente selecionável
   coordena shadow químico, benefício AMR e custos tipados de cell-step,
   comunicação, churn e RAM sem substituir as referências fixas.

A inspeção pós-edição não encontrou chamadas desconectadas nas novas
assinaturas, gates de capability contraditórios, whitespace inválido ou a
antiga guarda A02. Ela encontrou e corrigiu ainda o teste de ordem que assumia
Lie splitting e a interseção geométrica que confundia contato de aresta/canto
com área positiva de face.

## Matriz normativa (24 requisitos)

`IMPLEMENTADO EM FONTE, NÃO TESTADO` significa apenas cobertura estática do
caminho normal de produção; não é qualificação científica. `PARCIAL` significa
que permanece código obrigatório ausente.

| ID | Status nesta árvore | Evidência e limite preciso |
|---|---|---|
| H01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | MUSCL–Hancock/PPM/CTU/HLLC-HLLE permanecem em `src/hydro/hydro_core_solver.cpp`. O CTU agora reconstrói, com o mesmo PCM/PLM/PPM e Riemann do avanço, os dois fluxos de cada eixo tangencial do ghost a partir do halo imutável já preenchido por AMR/MPI; o incremento usa o volume-alvo e não cria autoridade persistente. Apenas a região de aresta/canto cujo stencil cruza uma segunda interface usa fallback normal conservativo explícito, contado em `ctu_transverse_degraded_states`; correções/axes remotos têm telemetria própria. Scratch do halo foi incluído em `HydroScratchBuffers::workspaceBytes`, na admissão AMR e nos diagnósticos runtime. Prova de ordem e estabilidade remota pertencem à campanha numérica posterior. |
| H02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Dual energy, autoridade de energia interna, recuperação de primitivas, floors e conversões comoventes: `src/hydro/hydro_core_solver.cpp`, `src/amr/amr_ghost_fill.cpp`, `src/hydro/gas_thermodynamic_refresh.cpp`. Falta qualificação numérica, não uma rota conhecida de fonte. |
| C01 | **BLOQUEADO POR DADOS EXTERNOS** | Operador, loader SHA-256/provenance, H/He não-equilíbrio, heating/shielding e cooling elementar estão em `src/physics/thermochemistry*.cpp` e `metal_cooling_table.cpp`. As tabelas físicas de produção exigidas pelos perfis não acompanham o pacote. |
| C02 | **PARCIAL / BLOQUEADO POR DADOS EXTERNOS** | Backends analítico/equilíbrio existem; `fixed_detailed` em `src/physics/chemical_network.cpp` permanece rede reduzida primordial + C+/CO. Faltam inventário multiiônico aprovado e taxas/tabelas verificáveis; sem esses dados não é lícito inventar reações, mas a ligação do inventário definitivo também permanece código futuro. |
| C03 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Temperatura, densidade, H2/CO, cooling time, erro/resíduos, histerese, dwell e transição conservativa estão em `ThermochemistryOperator::decideChemistryFidelity`/`BaryonicSourceStage`. `GasCompositionSidecar` schema 3 preserva o histórico e a extensão esparsa por identidade através de commit, transporte MPI, migração, regrid, rollback/restart e HDF5. `SourceRuntimeImpl::updateGlobalChemistryControl` reduz o shadow global e, em modo V02 conjunto, grava bias, amostras e promoções multinível no controlador canônico. A qualificação dos limiares continua pendente de validação, não de código conhecido. |
| C04 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Shadow detalhado amostrado, determinístico, limitado e sem commit: `ThermochemistryOperator::advanceRow`; redução global em `SourceRuntimeImpl::updateGlobalChemistryControl`. O controlador agora avança uma vez por intervalo Strang, não entre as duas metades. |
| I01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Regime efetivo versus frio resolvido, ownership térmico único e ausência de double cooling: `effective_multiphase_ism.cpp`, `baryonic_source_stage.cpp`, `gas_thermodynamic_refresh.cpp`. |
| I02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `prepareSgsShearTransfer` em `ism_dust_sgs.cpp` usa `sgs_eddy_viscosity_coefficient` no stress traceless `nu_t=C_v Delta^2 |S*|`: impulsos e trabalho são antissimétricos por face, a produção sai da energia resolvida e entra na lane SGS no mesmo commit, e o limite viscoso pede retry numérico. `SourceRuntimeImpl::executeDiffusionAndSgsShear` reutiliza o grafo físico local/coarse-fine/periódico/MPI, valida interfaces remotas coletivamente, governa staging e refresca caches. Advecção/reflux/regrid/restart, feedback e dissipação permanecem na autoridade SGS existente. O modelo continua EXPERIMENTAL/UNQUALIFIED e requer validação numérica futura, não código obrigatório conhecido. |
| S01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Elegibilidade Jeans/fluxo/boundedness, vizinhos AMR/MPI, nascimento conservativo, ID determinístico e birth epoch final: `star_formation.cpp`, `star_formation_gradients.hpp`, `SourceRuntimeImpl::execute`. |
| S02 | **BLOQUEADO POR DADOS EXTERNOS** | SSP, idade/metallicidade, cronologia cumulativa, associações e loader de yields existem em `stellar_population.cpp` e `source_runtime.cpp`. O schema 2 da sidecar populacional preserva também o cursor absoluto do último intervalo aceito, recusa gaps/replays e migra schema 1; ordinais e cursor só avançam no mesmo commit dos depósitos. Faltam tabelas licenciadas de produção; nenhuma equivalência FIRE-3/COLIBRE é declarada. O depósito em barreiras AMR multinível continua limitado por T02. |
| F01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Ramos mecânico resolved/unresolved, orçamento exclusivo térmico/SGS, depósito conservativo e refresh termodinâmico: `supernova_feedback.cpp` e `SourceRuntimeImpl::executePart2StellarPopulation`. |
| F02 | **PARCIAL / BLOQUEADO POR DADOS EXTERNOS** | Vento e pre-SN, loader e orçamento não aditivo existem em `stellar_wind.cpp`/`source_runtime.cpp`; falta dataset físico de vento de produção e a cronologia por nível exigida por T02. |
| D01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Transporte conservativo de elementos/dust/moléculas/SGS por `CompositionFluxAccumulator`; difusão MPI por `DistributedElementDiffusionCoupling`; capability coerente em `buildRuntimeCapabilityReport`; integração Strang usa `source_dt_code`. |
| A01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `planProductionAmrRegrid` separa `admitted_optional_operations`; `runAutomaticAmrRegrid` usa candidato integral governado + voto coletivo + publicação no-throw. Prolong/restrict de composição e identidades continuam nos operadores existentes. |
| A02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `advanceProductionAmrHydroRecursiveImpl`, `mergeFluxRegistersIntoPendingStore`, `DistributedRecursiveHooks::exchangeLevelFluxRegisters`, histories start/end e restart de `PendingFluxRegisterRecord`. Intervalos incompatíveis e duplicatas por origem falham fechados; área fina repartida é recomposta no subpasso. |
| A03 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `patchDecision` mantém massa/Jeans/cooling length como mandatórios e jump, metal/species/vorticidade/eventos como opcionais tipados. Todo backend Part 2 publica `chemistry_last_cooling_time_seconds`; `runAutomaticAmrRegrid` o converte ao tempo-code autoritativo. `AmrForceResolutionPolicy` oferece a referência histórica contabilizada e a política produtiva que bloqueia refinamento opcional abaixo do softening e falha coletivamente se requisito mandatório violar softening, `amr_max_level` ou a hierarquia 2:1. G01 continua condicional à validação. |
| A04 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `planProductionAmrRegrid` calcula severidade normalizada por critério e custos separados de cell-steps incrementais, faces MPI, churn e RAM; o modo conjunto ordena benefício/cell-step e usa MPI/churn como desempate, aplica piso com histerese e orçamento, enquanto mandatórios bypassam os tetos. `runAutomaticAmrRegrid` troca posições atual/driftada por `allgatherBytesBounded`, resolve o host somente no owner e propaga falhas/admissão coletivamente. O modelo permanece RESEARCH/UNQUALIFIED até validação numérica. |
| G01 | **CONDICIONAL À VALIDAÇÃO** | TreePM/KDK e auditoria AMR-softening existem. A especificação condiciona uma correção gravitacional multinível adicional ao resultado de qualificação futura; esta rodada não executou essa qualificação e não inventou novo solver. |
| T01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `SourceSplitPhase`, ordem de nove dispatches em `time_integration.cpp`, inversão determinística de tempo cosmológico, `makeSourceInterval` por metade e eventos somente post. Rollback global restaura estado, schedulers, force endpoint e checkpoint privado da fonte. |
| T02 | **PARCIAL** | `advanceProductionAmrHydroRecursiveImpl` agora invoca `SourceRuntimeImpl::executeAmrLevelSmoothSources` em ambas as metades de cada subintervalo, inclusive nos ranks sem linhas locais; termoquímica/dust/SGS usam o relógio `dt_L`, entram nos históricos temporais e shadow, promoções adaptativas e controlador são agregados uma vez por passo com rollback privado. A autoridade de eventos permanece O(1) na sidecar populacional e ganhou cursor absoluto persistente com rejeição de gap/replay e migração de schema. SF/SSP/feedback/ventos/difusão ainda usam a barreira global enquanto faltam depósitos cross-level e difusão no grafo multinível. |
| M01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `GasCompositionSidecar` mantém apenas as lanes básicas densas e armazena C+/CO adaptativos em vetores contíguos ordenados por `gas_cell_id`, com projeção implícita conservativa e retenção do histórico após rebaixamento. Crescimento em química, transporte Godunov MPI, difusão e demais commits usa `RetainedCapacityTransaction`/governador; regrid e migração dimensionam exatamente o candidato já admitido. Pack schema 3, restart com migração schema 1/2, HDF5, AMR prolong/restrict, ownership, rollback e depósitos de feedback/poeira preservam a mesma autoridade. O modo detalhado fixo continua denso como referência determinística. Qualificação de pico MPI fica para a campanha de testes reservada. |
| M02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | CPU scalar/SIMD/OpenMP/MPI, decomposição/migração e balanceamento ponderado permanecem nos runtimes existentes. GPU segue opcional e não foi transformada em dependência. |
| V01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Contratos de configuração/provenance, ledgers, snapshots HDF5, restart versionado e quatro decks científicos permanecem conectados em `config.cpp`, `snapshot_hdf5.cpp`, `restart_checkpoint.cpp` e `configs/chui_part3`. Aceitação end-to-end foi deliberadamente adiada. |
| V02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `V02ControlMode` seleciona `fixed_reference` ou `joint_error_resource_v1`. O estado versionado em `adaptive_resource_control.hpp` reside na sidecar de `SimulationState`, portanto acompanha rollback/restart sem segunda autoridade. `SourceRuntimeImpl` publica shadow, bias e demanda química multinível; `runAutomaticAmrRegrid` reserva crescimento esparso antes de expor headroom opcional, persiste decisões admitidas ou no-op e usa estimadores AMR tipados sem inventar erro físico universal. O modo fixo preserva bias químico e ordem AMR por `patch_id`; física mandatória nunca é reduzida pelo controlador. O modo conjunto permanece RESEARCH/UNQUALIFIED até qualificação. |

## Contratos alterados nesta campanha

- `StepContext::source_split_phase` é transitório e não cria nova autoridade de
  tempo. `IntegratorState` e `CosmologicalTimeline` continuam autoritativos.
- `PendingFluxRegisterRecord` continua a única autoridade persistente do reflux
  incompleto; não foi criado store paralelo. A semântica de
  `fine_area_accumulated_comov` passou a ser a cobertura espacial do subpasso
  corrente/completo, e não o máximo de uma contribuição por rank.
- O regrid automático conserva uma única publicação autoritativa de
  `SimulationState`; o candidato é scratch transacional governado e destruído
  após commit/rejeição.
- Difusão distribuída somente é elegível com composição elementar. A capability
  não promove o operador escalar legado nem oculta a dependência de MPI.
- `chui_v02_joint_error_resource_control_v1` é o único estado persistente do
  controlador conjunto. Química continua dona de composição/fidelidade e AMR
  continua dono da hierarquia; o sidecar guarda apenas histórico de decisão e
  reservas, participa da cópia de rollback/restart e não duplica campos físicos.

## Pendências obrigatórias verdadeiras

Não é correto declarar fechamento integral H01–V02 nesta árvore. Permanecem
obrigatórias a implementação T02 multinível bariônica e a ligação do
inventário definitivo de C02. C02 requer primeiro a especificação/dados verificáveis do
inventário detalhado; S02/F02 requerem tabelas físicas de produção. G01 é
deliberadamente condicional à campanha de validação prevista pela norma.

X01 permanece fora de escopo: nenhuma implementação de MHD, moving mesh, SPH,
CHIMES, Grackle, RT multifrequência ou nova família AGN foi introduzida.
