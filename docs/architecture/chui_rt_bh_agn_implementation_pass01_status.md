# CHUÍ RT/BH/AGN — implementation ledger, pass 01

**NO BUILD, NO TESTS, NO BENCHMARKS, NO SCIENTIFIC VALIDATION WERE RUN BY THIS AGENT.**
Compile viability of every new or changed file is unverified. Every test is
`NOT_RUN`. Statuses: `IMPLEMENTED_UNVALIDATED`, `PARTIAL`, `FAIL_CLOSED`,
`NOT_IMPLEMENTED`, `EXCLUDED_RESEARCH`. Qualification of every row: `NOT_RUN`
(or `BLOCKED` where a dependency is missing). Nothing is VERIFIED, VALIDATED,
CALIBRATED or PRODUCTION_QUALIFIED.

## 1. Baseline and repository

| Field | Value |
|---|---|
| Repository | `AlexandreKoz/Chui-dev` |
| Base | `main` @ `6310aba6befea79a4e13e073d07ad572696bb7c8` (verified with `git ls-remote`) |
| Feature branch | `claude/funny-volta-bcqf37` (session-mandated isolated branch) |
| Pre-existing changes | none (clean worktree) |
| Specifications | five v0.1 files copied unchanged into `docs/architecture/` (SHA-256 in the preflight) |
| Amendments | `chui_rt_bh_agn_compatibility_amendment_v01.md` (A1–A15) |
| Remote CI on push | none (`ci.yml` triggers on `pull_request` and pushes to `main` only) |
| Commits | `14cf3a7` shared contracts · `c38fedf` BH engine · `3187f23` AGN · `14fc64f` chemistry interface (A5) · `3c19832` radiation · `d5c555a` P1 · `c22728c` integration · `d4e3599` restart v27 · `a449019` tests/docs · final ledger commit |

## 2. Implemented source (grouped)

* **Shared (B01–B04):** `include/cosmosim/physics/rt_bh_agn_frame.hpp` (proper-frame table, `PeriodicBox`); `include/cosmosim/physics/source_transaction.hpp`, `src/physics/source_transaction.cpp` (prepare/validate/authorize/commit batch); `include/cosmosim/core/config_rt_bh_agn.hpp`, `src/core/config_rt_bh_agn.cpp`, `src/core/internal/config_rt_bh_agn_{tables,validation}.hpp`, `src/core/config.cpp` (registry, parse, inert-key rejection, admission, hash-preserving dump), `include/cosmosim/core/config.hpp` (three nested structs).
* **BH engine:** `include/cosmosim/physics/black_hole/` and `src/physics/black_hole/`: `bh_engine_state` (module sidecar store), `bh_supply` (kernel, Bondi, Rosas-Guevara, Hopkins–Quataert torque, decomposition, resolved-flux P1 measurement), `bh_disk` (exact piecewise reservoir), `bh_budget` (Eddington, ledgers, emission budget), `bh_dynamics` (Chandrasekhar/Tremmel), `bh_merger` (bound-pair criteria, disjoint selection, record merge, Peters time), `bh_spin` (Kerr ISCO, thin-disk spin-up, BZ extraction; P1), `bh_neighborhood` (proximity filter, owner-target exchange interface, serial exchange), `bh_engine` (interval coordinator), `bh_seed_provider` (halo seeding transaction).
* **AGN:** `include/cosmosim/physics/agn/agn_injection.hpp`, `src/physics/agn/agn_injection.cpp` (thermal continuous/reservoir, kicks, bipolar jets, winds, radiation disposition, refinement requests); `agn_magnetic.*` (P1 CT-compatible edge-potential injection and energy quadratic).
* **Chemistry interface (shared file, A5):** `thermochemistry.hpp/.cpp`, `thermochemistry_operator.hpp/.cpp`, `source_interval.hpp/.cpp` — additive `TransportedPhotoRates` input, provider hook, `SourceLedger::transported_photoionizations`. Necessary to couple RT without a second network; defaults unchanged.
* **Radiation (`cosmosim_radiation`):** `include/cosmosim/radiation/` + `src/radiation/`: `radiation_closure`, `radiation_spectrum`, `radiation_transport` (topology, ghosts, HLL/GLF, MC joint limiting, reflux), `radiation_matter_coupling` (grey Newton, photoionization, photo-rates), `radiation_cosmology`, `radiation_runtime` (allocation, PD-ARS subcycling, retries, staged commit), `radiation_amr`, `radiation_sources` (AGN), `radiation_angular` (P1 rays, VET, Monte Carlo).
* **Core storage:** `RadiationMomentSoa` in `simulation_state.hpp`, consistency in `simulation_state_structures.cpp`, memory accounting in `memory_accounting.cpp`.
* **Integration:** `src/workflows/internal/rt_bh_agn_bridge.{hpp,cpp}` (MPI support exchange, serial FoF provider, epoch composition); `src/workflows/source_runtime.cpp` (constructor guard mapping, bridge construction, discrete-epoch branch, touched-row refresh).
* **Persistence:** `src/io/restart_checkpoint.cpp` (schema v27, `/state/radiation`, digest coverage, schema-group validation).
* **Build:** `CMakeLists.txt` (new sources, `cosmosim_radiation` target/alias/install, tests, benchmark).
* **Reused unchanged:** legacy `BlackHoleAgnModel` (comparator), particle-ID precommit registry, module-sidecar migration/restart, `FofHaloFinder`, `exchangeBoundedAlltoallBytes`, `MpiContext` reductions, `GasCompositionSidecar::scaleRow`, `synchronizeMetalMirrorFromComposition`, `remapSchedulerByParticleId`, `commitParticleMigration`.

## 3. Coverage matrix

| ID | Spec | Source evidence | Runtime entry | Admission | Persistence | Memory/parallel | Tests | Status | Missing / next repair |
|---|---|---|---|---|---|---|---|---|---|
| B01 | master §4 | `RadiationMomentSoa`, `BlackHoleEngineStore`, transactions on canonical lanes | bridge | — | v27 / module block | lazy, zero when disabled | radiation_core, bh_engine_interval | IMPLEMENTED_UNVALIDATED | audit duplicate-authority claims |
| B02 | master §5 | `SourceTransactionBatch` | engine, AGN, seeding, RT commit | — | n/a | vote via `BhSupportExchange::allAgree` | source_transaction | IMPLEMENTED_UNVALIDATED | fatal mid-commit -> checkpoint recovery is a workflow policy, not coded here |
| B03 | master §5 | `ProperFrameContext` | all | — | — | — | numerics | IMPLEMENTED_UNVALIDATED | RT a(t) interpolation is linear in proper time (approx.) |
| B04 | master §7 | `config_rt_bh_agn.*`, `config.cpp` | `validateConfig` | yes | normalized dump (enabled only) | — | rt_bh_agn_config | IMPLEMENTED_UNVALIDATED | docs key table complete; examples test passes only if decks parse (NOT_RUN) |
| BH01 | BH §3.1 | `bh_seed_provider`, bridge FoF provider | bridge (cadence) | `halo_fof_serial`, 1 rank | module block | O(N) FoF per cadence | config refusal only | PARTIAL | distributed provider FAIL_CLOSED; no seeding unit test; seed composition not retained in BH |
| BH02 | BH §3.2 | `bh_dynamics`, engine DF block | engine | `unresolved_chandrasekhar` | cumulative ledger | rank-local samples | numerics | PARTIAL | cross-rank background; DF timestep limit not fed to scheduler; gas drag NOT_IMPLEMENTED |
| BH03 | BH §4 | `bh_supply`, engine rounds 1–3 | engine | supply models | last_* diagnostics | owner-target rounds | numerics, interval | IMPLEMENTED_UNVALIDATED | torque calibration is a user input; decomposition classifier unverified vs SIMBA code |
| BH04 | BH §5 | `evaluateWithdrawal` | engine round 4 | — | disk lanes | joint per-row limit, vote | interval, transaction | IMPLEMENTED_UNVALIDATED | magnetic sink = retained field (policy text only); max 0.5 per-row fraction is a numerical limit |
| BH05 | BH §5 | `bh_disk`, engine substeps | engine | `scalar_reservoir`/`instantaneous` | disk lanes | — | numerics | IMPLEMENTED_UNVALIDATED | — |
| BH06 | BH §6 | `bh_budget`, engine | engine | `bh_mass_ledger` | cumulative ledgers | — | numerics | IMPLEMENTED_UNVALIDATED | L/L_Edd uses eps_r of the config (no state-dependent efficiency) |
| BH07 | BH §8 | `bh_merger`, engine merger block | engine | `bound_pair` | lineage lanes | local only | numerics | PARTIAL | cross-rank pairs; tidal criterion; ownership-ledger compatibility of removal |
| RT01 | RT §3 | `radiation_spectrum`, `RadiationMomentSoa` | runtime | profiles | v27 | 4G·8·N | radiation_core | IMPLEMENTED_UNVALIDATED | broad/agn_hard FAIL_CLOSED |
| RT02 | RT §2, §5 | `radiation_closure`, `radiation_transport` | runtime | `rt_riemann_solver`, `rt_reconstruction` | — | single rank | radiation_core | IMPLEMENTED_UNVALIDATED | topology build O(P x ghosts x P) |
| RT03 | RT §6 | `RadiationRuntime::advanceInterval` | bridge | `rt_integrator` | — | work copies (see PERF01) | radiation_core | IMPLEMENTED_UNVALIDATED | first-order fallback substeps; no asymptotic test |
| RT04 | RT §6 | `radiation_matter_coupling`, A5 hooks | runtime | grey reference only | — | — | radiation_core | PARTIAL | photoionization reconciliation not wired -> FAIL_CLOSED; mixed-frame coefficients to verify |
| RT05 | RT §4, §6 | `radiation_cosmology`, runtime a^4 handling | runtime | reduced c refused in cosmology | v27 scale factor | — | radiation_core | PARTIAL | expansion-work ledger not reported separately; cosmological RSLA refused |
| RT06 | RT §7 | topology ghosts/reflux, `radiation_amr` | runtime | MPI/regrid/subcycling refused | v27 | single rank | radiation_core | PARTIAL | MPI halo exchange, regrid transfer, migration FAIL_CLOSED |
| RT07 | RT §7 | `radiation_sources` | bridge | AGN sources (SED table) | pending budget debit | — | none | PARTIAL | stellar sources FAIL_CLOSED; no AGN-source unit test |
| AGN01 | AGN §2 | engine pending budgets, injector disposition | engine/AGN | `agn_*` | module block | — | interval | IMPLEMENTED_UNVALIDATED | spin-powered reserve absent (P1) |
| AGN02 | AGN §3 | injector statistics/commit rounds | AGN | — | — | owner-target, vote | interval | IMPLEMENTED_UNVALIDATED | events occur at the synchronized epoch (no inactive-target latency path needed) |
| AGN03 | AGN §4.1–4.2 | thermal + kick quadratic with recoil | AGN | modes | pending lanes | — | numerics, interval | IMPLEMENTED_UNVALIDATED | stochastic targeting not implemented (deterministic only) |
| AGN04 | AGN §4.3 | wind/jet blocks | AGN | velocity domain | pending wind mass | — | numerics (quadratic) | IMPLEMENTED_UNVALIDATED | no wind/jet integration test |
| AGN05 | AGN §7 | `AgnRefinementRequest` | AGN report | `agn_refinement_*` | — | — | none | PARTIAL | requests not consumed by `amr_production_regrid` tagging |
| P1-BH01 | BH §3.1, §4.3 | `resolvedControlSurfaceFlux` | none | refused | — | — | none | PARTIAL / local seeds NOT_IMPLEMENTED | no qualified seed criteria (no invented thresholds) |
| P1-BH02 | BH §7 | `accreteKerrThinDisk`, Kerr functions | none | refused | J lanes exist | — | numerics | PARTIAL | warp alignment (Fiacconi/Cenci) NOT_IMPLEMENTED; not wired |
| P1-BH03 | BH §7–8 | `extractBlandfordZnajek`, `petersCircularInspiralTime` | none | refused | — | — | numerics | PARTIAL | thick/thin/slim states, super-Eddington, remnant/recoil fits NOT_IMPLEMENTED |
| P1-RT01 | RT §7 | `radiation_angular` (rays, VET tensor, MC) | none | `rt_transport != m1` refused | — | bounded queue | radiation_core (VET, RNG) | PARTIAL | not wired; VET transport backend, hard UV/X-ray, dust IR NOT_IMPLEMENTED |
| P1-AGN01 | AGN §4.4, §6 | `agn_magnetic` | none | refused | — | — | numerics | PARTIAL | CT-authority wiring + AMR reflux-curl; nozzle NOT_IMPLEMENTED |
| P1-HPC01 | perf §6 | none | — | GPU not offered for RT | — | — | — | FAIL_CLOSED | no device kernels (no fake backend) |
| INT01 | master §6 | bridge + source_runtime | discrete epoch | — | — | collective rounds | interval (engine level) | PARTIAL | RT not Strang-split with the smooth chemistry stage; end-to-end runtime test NOT_WRITTEN |
| INT02 | master §7 | restart v27, module block | writer/reader | — | yes | — | none | PARTIAL | science-snapshot export, restart-verifier field comparison and v27 roundtrip test NOT_WRITTEN |
| PERF01 | perf §2–3 | memory accounting `rt.*`, `radiationCanonicalBytes` | — | — | — | — | bench | PARTIAL | NO MemoryGovernor reservations for RT work copies (~5 x canonical), BH/AGN scratch, ray/MC buffers |
| PERF02 | perf §4–5 | SoA lanes, bounded exchange, rank boxes | — | — | — | — | — | PARTIAL | no OpenMP; per-round O(N_owned) scans; topology O(P^2) |
| RES | master §3 | — | — | — | — | — | — | EXCLUDED_RESEARCH | S_N, adaptive angular/spectral, global implicit, level-dependent c, radial disk PDE, PN/GRMHD, relativistic jets, CRs |

Counts (35 rows): IMPLEMENTED_UNVALIDATED 17 · PARTIAL 16 · FAIL_CLOSED 1 · EXCLUDED_RESEARCH 1 (NOT_IMPLEMENTED sub-items are listed inside the PARTIAL P1 rows).

## 4. Runtime trace (source evidence; not executed)

1. `.param.txt` → `normalizeValidateFreeze` (`src/core/config.cpp`, RT/BH/AGN parse block after `bh_speed_of_light_si`) → `internal::validateRtBhAgnConfig` (`src/core/config_rt_bh_agn.cpp`).
2. `SourceRuntimeImpl` constructor (`source_runtime.cpp`): seeding guard mapping; `m_rt_bh_agn = RtBhAgnSourceBridge(...)` when `reference_v1` or `rt_enabled`.
3. `executeDiscreteSourceEvents` (coarse end, hierarchy synchronized, units: code time + proper seconds from the cosmology background, owner: SimulationState) → legacy `m_black_hole.apply` skipped for `reference_v1` → `RtBhAgnSourceBridge::executeDiscreteEpoch`.
4. `BlackHoleEngine::advanceInterval` (proper dt, `ProperFrameContext`): rounds `kRadiusBins` → `kKernelStatistics` → `kCounterRotation` → `kWithdrawal` (`SourceTransactionBatch` prepare/validate → `allAgree` → commit on gas lanes) → owner update (disk ODE, ledger, budgets → `chui_bh_engine_v1`) → DF → mergers (`commitParticleMigration`, `remapSchedulerByParticleId`) → `store.store`.
5. `applyHaloSeeding` (cadence; precommit IDs; transaction) → caller appends scheduler elements and created IDs.
6. `AgnFeedbackInjector::apply`: radiation disposition → `kAgnStatistics` → amplitudes → `kAgnCommit` (batch, vote, commit) → BH recoil/mass update → store.
7. RT: `allocateRadiationState` (first use) → `radiationRecordsMatchPatches` → `buildAgnRadiationSources` (debit pending radiation) → `RadiationRuntime::advanceInterval` (PD-ARS subcycles on comoving a^4 moments, grey implicit exchange, edge transfer) → staged gas transaction + lane commit.
8. Back in SourceRuntime: compatibility mirrors, touched-row thermodynamic refresh, scheduler append for seeds, one collective gravity-source generation bump.
9. Checkpoint: `writeRestartCheckpointHdf5` → v27 if radiation non-empty (+ module blocks); reader validates and restores.
Paths reachable only from tests/API: Kerr/BZ, rays/VET/MC, magnetic injection, resolved flux, photoionization coupling, explicit `kExplicitReference` comparator (config-reachable).

## 5. Physical and numerical invariants (implemented formulas)

* Gas transactions: total-energy convention `u' = (E + dE - |p+dp|^2/2m')/m'`; withdrawal leaves `u, v` unchanged; inelastic mixing heats (test_source_transaction).
* Accretion: gas → disk mass/momentum exact; captured thermal+kinetic energy and capture dissipation in named ledgers; particle mass `M_bh + M_d + M_wind,pending`.
* Disk: exact exponential/linear pieces with exact crossing time; `D = M + S dt - M'` closes to round-off.
* Ledgers: legacy `dM_bh = dm0`; relativistic `dM_bh = (1 - eps_acc) dm0 - E_spin/c^2`, wind mass excluded from dm0; never mixed.
* AGN: budgets debited by released amounts; kick/jet amplitude from `C A^2 + B A = Q` with recoil; jet lobes normalized to equal and opposite momenta; wind energy `dm(u_w + v_w^2/2) + v_bh·dp_rest + dm v_bh^2/2`; owner closure residual reported.
* RT: M1 (Levermore), admissible cone, joint theta, exact eigen-bounds, PD-ARS eqs. 28–29, RSLA via (E, G = F c_hat/c), matter receives c G^0 and G (c/c_hat invariant), edge transfer `H E_l phi_low E_g` exact exponential, a^4 storage.
* Deviations from spec: frame velocity frozen in the grey Newton; first-order fallback; Plummer-equivalent merger potential; DF rank-local; RT single rank.

## 6. Verification inventory (all NOT_RUN)

| Gate | Test source (registered name) | Covered | Missing |
|---|---|---|---|
| V01 | `unit_radiation_core` (closure limits, wave speeds, rotational invariance of bounds, joint limiter, periodic pulse conservation) | partial | oblique beam, smooth convergence slope >= 1.8, HLL/GLF comparison |
| V02 | `unit_radiation_core` (grey implicit LTE conservation) | partial | Marshak, dynamic diffusion, moving pulse, retry/rollback, RSLA invariant run |
| V03 | `unit_radiation_core` (spectral identities, redshift transfer conservation) | partial | Strömgren (needs photoionization), shadow/crossing beams, homogeneous expansion run, emission-once test |
| V04 | `unit_rt_bh_agn_config` (disabled equivalence of normalized text, refusals) | partial | AMR/rank fronts, restart v27 roundtrip, empty ranks |
| V05 | `unit_bh_agn_numerics`, `integration_bh_engine_interval`, `unit_source_transaction` | partial | seeding uniqueness/restart, torque decomposition on particle data, orbit decay |
| V06 | `unit_bh_agn_numerics` (Kerr limits, spin-up, misalignment refusal, BZ area theorem) | partial | vector-J alignment, efficiency transitions, remnant fits (blocked: no pinned fits) |
| V07 | `unit_bh_agn_numerics` (quadratic, kinetic increment, magnetic div/energy), `integration_bh_engine_interval` (thermal closure) | partial | jets/winds integration, concurrent events, rollback, refinement |
| V08 | `integration_bh_engine_interval` (engine-level interval) | partial | full SourceRuntime interval, rank equivalence, galaxy/cluster/cosmological decks (`configs/rt_bh_agn/*.param.txt` written, not run) |
| V09 | `bench_rt_transport_throughput` (driver only) | partial | memory slope, scaling ladder, CPU/GPU (GPU blocked) |

Planned commands (execution campaign): `cmake --preset cpu-only-debug`, `cmake --build --preset build-cpu-debug`, `ctest --preset test-cpu-debug -L rt_bh_agn --output-on-failure`, then the full suite (disabled-path regression), HDF5 preset for v27.

## 7. Severity-ranked gap ledger

| Sev | Gap | Evidence | Consequence / current behaviour | Proposed narrow repair |
|---|---|---|---|---|
| CRITICAL | No compilation or test evidence for ~11k new lines | this pass | unknown build breaks | execution campaign: build, fix, run `-L rt_bh_agn` and full suite |
| HIGH | Photoionization reconciliation not wired | A5; `config_rt_bh_agn.cpp` refusal | production ionizing RT refused (safe) | store per-row RT absorbed energy + rates, consume `transported_photoionizations` in `BaryonicSourceStage` commit, correct E_g, then lift refusal |
| HIGH | RT / BH / AGN scratch not reserved through MemoryGovernor | runtime work copies, exchange buffers, filter maps | allocation outside admission | reserve `kPhaseResident`/`kScratchArena` before each round/interval; bound filter maps |
| HIGH | Merger particle removal vs workflow ownership ledger | `bh_engine.cpp` merger block | possible ledger/count assertion at the next validation | add removed-ID sink to `StepContext`/ownership ledger; test |
| HIGH | Multi-rank RT and RT regrid | admission refusals A14 | RT single-rank only (safe) | patch halo exchange per stage; moment transfer in regrid/migration payloads |
| MEDIUM | Cross-rank mergers and DF backgrounds | engine (local only) | pairs/backgrounds split by rank boundaries are missed (biased, not unsafe) | route BH/particle queries through the support exchange |
| MEDIUM | Refinement requests not consumed | `AgnInjectionReport::refinement_requests` | unresolved injection recorded only | map requests to AMR tags in `amr_production_regrid` |
| MEDIUM | Science-snapshot export and restart-verifier comparison of new state | INT02 | no output of radiation fields/engine lanes | add v8 snapshot groups; compare in `restart_verification.cpp` |
| MEDIUM | Mixed-frame coefficients, He I Verner coefficients, Tchekhovskoy BZ coefficients are transcriptions | A6, A13, `bh_spin.hpp` | possible coefficient errors | audit 02 against primary papers |
| MEDIUM | Topology/locator O(P^2) and per-round O(N_cells) scans | `radiation_transport.cpp`, `radiation_angular.cpp`, filter rounds | scaling cost at many patches | spatial hash of records; persistent leaf index |
| LOW | Seed composition not tracked in the BH; DF limiter not fed to the scheduler | `bh_seed_provider.cpp`, `bh_dynamics.cpp` | minor ledger detail / timestep | extend record; submit candidate timestep |
| LOW | RT a(t) and H(t) linear interpolation within an interval | `radiation_runtime.cpp` | small cosmological timing error | query the cosmology service per substep |

## 8. Handoff

First audit targets: `src/physics/black_hole/bh_engine.cpp` (rounds, owner update, merger removal), `src/physics/agn/agn_injection.cpp` (amplitude and closure), `src/radiation/radiation_transport.cpp` and `radiation_runtime.cpp` (fluxes, IMEX staging, commit), `src/workflows/internal/rt_bh_agn_bridge.cpp` and the `executeDiscreteSourceEvents` branch, `src/io/restart_checkpoint.cpp` v27 edits (all `== v26` changes), `src/core/config.cpp` insertion points. Cross-campaign risks: shared thermochemistry files (A5), restart schema numbering with concurrent MHD work, SourceRuntime discrete-epoch ordering with stellar feedback, `SimulationState` layout change (new `radiation` member). No merge, deployment or next campaign was started.

## Slice log

- S0 preflight; S1 shared contracts; S2 BH engine; S3 AGN; S4 chemistry interface; S5 radiation; S6 P1 extras; S7 integration; S8 restart v27; S9 tests/docs; S10 ledger and final sweep (`git diff --check` clean outside the verbatim specification files, which carry Markdown hard-break trailing spaces and were copied unchanged).
