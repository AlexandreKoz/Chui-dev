# CHUI final implementation and scientific-data closure, pass 2 (2026-11-04)

**Evidence status.**

* CHUI was **not** configured, compiled, linked, tested, benchmarked or run. Every CHUI test
  and benchmark named below is **NOT RUN — USER-RESERVED**.
* Executed in this pass: data downloads, the Python converters and validators in
  `tools/chui_data/` (they read text and HDF5 only), SHA-256 checksums, and git.

This file supersedes the status columns of
[`chui_final_design_closure_20261104.md`](chui_final_design_closure_20261104.md) (pass 1). That
file's call graph and contracts remain valid.

| Item | Value |
|---|---|
| Branch | `claude/epic-fermi-2bcq65` (verified HEAD at start: `408b6c88e2edae99953ed8375c6f0f83ba75f3f9`) |
| Data workflow | `tools/chui_data/fetch_and_convert.sh <dir>`. A second independent run reproduced every derived file bit for bit. |
| Data record | `tools/chui_data/data_sources.md` §"2026-11-04 acquisition and conversion" (checksum manifest, residuals, limitations) |

## 1. Commits of this pass

| SHA | Item | Content |
|---|---|---|
| `ae1424a` | D1–D5 | NuGrid/W7 yield table, SB99 wind table, 49-redshift WSS09 table, validated FG20 TREECOOL, R18/C06 provenance, single-Z wind clamp reporting, `configs/chui_final` decks. |
| `7b2aca1` | R1 | Physics-aware AMR indicators: cooling length (mandatory); metallicity, species and vorticity (optional); mandatory-first budget. |
| `7e521f7` | R2 | Event-aware predictive refinement from the canonical SSP chronology. |
| `a2a0dfc` | R3 | Cooling-time / stiffness criterion for adaptive chemistry fidelity. |
| `c9fec79` | V5 | Kernel and MPI benchmark sources. |
| this commit | docs | This ledger and the status pointers. |

## 2. Scientific data (results of the executed converters and validators)

| Product | Status | Derived sha256 |
|---|---|---|
| SSP yields (NuGrid Fryer12-delay + W7) | **AVAILABLE — converted, VALID** | `e2161646eadc4119e1b4976045fd6dad7876c7b4940a0c9174e94e1fa56e8f14` |
| SSP winds (Starburst99, Z = 0.014 Geneva v00) | **AVAILABLE — converted** (single Z) | `ba232a1cfe2c1996cadf55f14df2476280f8f56e81fe25e2e06c90464b655da5` |
| WSS09 element cooling, z = 0–8.989 (49 tables) | **AVAILABLE — converted, VALID**. UVB-paired use remains refused (HM01 vs FG20). | `4bbca8a9e5123226392f077f852b89ec551ef4d233b366b3d9c8ee443b143410` |
| FG20 TREECOOL (×0.68 heating) | **AVAILABLE — VALID** | `4bfc5d9bd6f0541addd01256ac20b5e6bc8a7aa5fc8849b9cbdf34cb0882b313` |

Key quantitative facts (details in `data_sources.md`):

* **Yields.**
  * Ejecta closure: 2.1e-4.
  * Tracked net residual: ≤ 1.656e-2, against a declared tolerance of 0.02.
  * Every tracked gross ejecta is positive.
  * Untracked mass stays in the star particle: a median of 2e-4 and a maximum of 1.7e-2 of
    massive or AGB ejecta.
  * SN Ia: 11.8% of the W7 ejecta (Ni, Ar, Mn, Cr) is excluded, explicitly.
  * Ia has one metallicity only.
* **Winds.**
  * Cumulative energy: E = 1.61e48 erg/M☉. Cumulative momentum: P = 77.5 M☉ km/s per M☉.
  * The momentum is SB99's own Σ Ṁ v∞, integrated in time; it is not derived from the energy.
  * The energy quadrature check agrees to 0.78%.
  * Limitations: Kroupa IMF rather than Chabrier, and one metallicity (other Z are now reported
    as clamps).
* **Cooling.**
  * Values round-trip float32 bit-exactly.
  * H, He and Compton are not converted.
  * Per-element non-additivity relative to `Total_Metals`: median 3–7% for T > 10⁴ K.

**Remaining genuine data blockers.**

* An HM01 H/He photorate table consistent with WSS09, or FG20-consistent element tables. Without
  one, WSS09 + FG20 is usable only through the explicit mismatch override.
* Multi-metallicity SB99 runs. stsci.edu was blocked by the proxy (403), and only Z = 0.014
  outputs were found.
* Multi-metallicity SN Ia yields.
* The primary texts of Jura 1975 and Black & Dalgarno 1977.

## 3. Requirement matrix (H01–X01, updated)

| ID | Status | Basis / exact remaining work |
|---|---|---|
| H01 CTU at interfaces | **PARTIAL** | No change to CTU transverse reconstruction. At MPI and coarse/fine faces, deep stencil layers beyond the delivered ghost depth are still zero-gradient extrapolated. They are counted (`deep_layer_extrapolated_ghosts`, `ctu_transverse_degraded_states`), never silent. Next edit: request `ghost_layers = productionAmrHydroGhostLayers(...) + 1` for CTU in `exchangeDirectedHydroGhosts` and fill the corner and transverse ghosts in `amr_ghost_fill.cpp`. The N1 regression set (smooth-wave convergence, shock across levels, high-Mach remote) is not written. |
| H02 dual energy | **IMPLEMENTED IN SOURCE — NOT RUN** | The owner's specific e is transported to remote ghosts and consumed for the auxiliary lane (pass 1, `408b6c8`). Its dedicated regression is not written. |
| C01 thermochemistry | **IMPLEMENTED IN SOURCE — NOT RUN**; UVB pairing **BLOCKED BY EXTERNAL DATA** | R18 and C06 are now secondary-verified. The gate is unchanged. |
| C02 network references | **EXPERIMENTAL/UNQUALIFIED** | Unchanged. |
| C03 adaptive chemistry | **IMPLEMENTED IN SOURCE — NOT RUN**, **EXPERIMENTAL/UNQUALIFIED** | Composition thresholds plus the R3 cooling-time/stiffness criterion, canonical-sidecar bias, hysteresis and dwell, all logged. There is still no solver-error-estimate criterion, because per-row error estimates are not retained. |
| C04 shadow | **EXPERIMENTAL/UNQUALIFIED** | Unchanged; it never mutates state. |
| I01, I02 | **EXPERIMENTAL/UNQUALIFIED** | Unchanged. |
| S01 | **IMPLEMENTED IN SOURCE — NOT RUN** | Unchanged. |
| S02 SSP | **IMPLEMENTED IN SOURCE — NOT RUN**; yields **AVAILABLE** (D1) | Monotone quadrature from pass 1. |
| F01 feedback | **IMPLEMENTED IN SOURCE — NOT RUN** | Pass 1 global associations. The MPI end-to-end test is not written (V2). |
| F02 winds | **IMPLEMENTED IN SOURCE — NOT RUN**; data **AVAILABLE** for Z = 0.014 only | Other metallicities: **BLOCKED BY EXTERNAL DATA**. |
| D01 diffusion | **IMPLEMENTED IN SOURCE — NOT RUN** | Pass 1. The legacy scalar operator remains single-rank and is refused explicitly. |
| A01 regrid | **IMPLEMENTED IN SOURCE — NOT RUN** | Pass 1 octet co-location. |
| A02 recursive AMR | **IMPLEMENTED IN SOURCE — NOT RUN** | Pass 1. The end-to-end MPI test (V1) is not written. |
| A03 indicators | **IMPLEMENTED IN SOURCE — NOT RUN**, with one restriction | Mandatory: Jeans, mass and cooling length. Optional: jump, metallicity, species and vorticity. Mandatory requests are scheduled first, and mandatory truncation is reported. The cooling-length criterion is supported **only** with the legacy in-solver cooling. With Part 2 backends the key is refused explicitly (**PARTIAL** for that combination): SourceRuntime must publish per-row cooling times, for example `ThermochemistryOperator::coolingHeating` evaluated in the regrid planner through a source-owned callback. Feedback structure is covered by the R2 forcing. Gravity adequacy keeps the counted `refined_below_gravity_softening` regime; no new solver. |
| A04 predictive AMR | **IMPLEMENTED IN SOURCE — NOT RUN**, **EXPERIMENTAL/UNQUALIFIED** | Forecasts the actual chronology (stochastic mode). Expectation mode uses a ≥ 0.5-event heuristic. Limitation: a drifted host on another rank is not forced. |
| G01 TreePM / zoom | **IMPLEMENTED IN SOURCE — NOT RUN** | Unchanged. |
| T01 rollback | **IMPLEMENTED IN SOURCE — NOT RUN** (governed full-copy reference) | §5. A bounded journal is not required for correctness; it remains a memory optimisation (**PARTIAL** against memory architecture v0.2 only for very large runs). |
| T02 multirate | **PARTIAL** | §4. Source operators run at coarse cadence, which is not the fully level-local coupling the original design requires. |
| M01 memory | **IMPLEMENTED IN SOURCE — NOT RUN** | New admitted owners: `hydro.amr_regrid.cooling_times`. |
| M02 performance | **IMPLEMENTED IN SOURCE (benchmarks) — NOT RUN** | `bench_chui_final_kernels`, `bench_chui_final_mpi`. No results. |
| V01 outputs / restart | **OUT OF SCOPE BY DESIGN** (snapshot restart) / **IMPLEMENTED IN SOURCE — NOT RUN** (checkpoints) | Science snapshots are analysis-only. `docs/output_schema.md` §1 says so, and the restart reader rejects `science_snapshot` files. The authoritative restart is the checkpoint, which carries the full `SimulationState`: composition, populations, sidecars, AMR histories and registers. |
| V02 error/resource control | **PARTIAL** / **EXPERIMENTAL/UNQUALIFIED** | A fixed-fidelity reference stays the default (all new keys default to 0). The mandatory/optional priority now gives a resource-aware refinement budget. There is still no coordinated controller using measured error and node memory to pick budgets. Next: feed `MemoryGovernorSnapshot::headroom_bytes` into `amr_regrid_max_operations` for optional requests only. |
| X01 scope | **OUT OF SCOPE BY DESIGN** | No MHD, RT, CR, new AGN engine, GR or dark-sector physics added. |

## 4. N2: operator time authorities under recursive subcycling

| Operator | Time authority | Cadence | Proper-time coverage | Interval-limiting constraint | State valid at intermediate hydro times? | Rejected step | Splitting |
|---|---|---|---|---|---|---|---|
| Collisionless gravity (TreePM KDK) | `TimeCoordinator` / `IntegratorState` | Coarse interval (rung zero); the DMO hierarchical segment is unchanged | [t, t+Δt] once | Gravity η√(ε/\|a\|) per particle, not relaxed | Particles drift across the interval | Full-copy restore plus force endpoint | KDK, second order |
| Gas gravity kicks | Same KDK | Half-kicks at the coarse boundaries, all levels | Once per interval | Cell acceleration criterion, not relaxed | Gas velocity inside substeps has no gravity impulse | Same | Strang (gravity ⊕ hydro), second order in coarse Δt |
| Hydro transport | `advance*Recursive` level clocks | Level L: Δt/r^(L−Lmin) | Each cell exactly once over its own interval | Level-relaxed CFL with a collective guard | Yes, at each level clock (temporal ghosts) | Same | Conservative FV with reflux |
| Hydro source terms (legacy cooling, effective-ISM relaxation) | Sweep update context | Per level substep | Each substep | Internal subcycling | Yes | Same | Level-local |
| Part 2 thermochemistry, dust, SGS | `SourceRuntime` (`makeSourceInterval`) | Two coarse-cadence half steps around hydro | Full interval as two physical half intervals | Internal error-controlled subcycling; rejection leads to retry | **No.** Recursive level clocks are still not source-local. | Stage transaction plus full-copy | Strang reference at coarse cadence |
| Star formation | `SourceRuntime` | Once in the post half | Full accepted interval; birth epoch is the absolute endpoint | SF timestep candidate, not relaxed | No | Transaction and full copy | Discrete endpoint event, not duplicated |
| SSP return, events, feedback, winds | `SourceRuntime` plus cumulative chronology | Once in the post half | Full interval; events by cumulative Poisson, never duplicated | none (cumulative) | Deposits at the synchronized step end | Population sidecar restored | Discrete endpoint event |
| Element diffusion | `SourceRuntime` | Two coarse-cadence half steps, each with internal global subcycles | Full interval as two half intervals | Parabolic candidate plus internal stable Δt | No | Collective all-or-nothing | Strang reference |
| Regrid, output, checkpoint | Coordinator | Synchronized boundaries only | — | — | — | Full copy | — |

**Assessment against the original T02.** Hydro, its source terms and the reflux are level-local,
and every interval is covered exactly once with positive duration and correct chronology. Nothing
is double counted, skipped or stale-accelerated.

The canonical `SourceRuntime` operators now use deterministic Strang splitting relative to the
fine hydro clocks. The coarse Δt is **not** bounded by any thermochemical cooling time: no such
candidate exists in `updateAdaptiveTimeBins`. In deeply refined, rapidly cooling gas, the
splitting error therefore scales with the coarse interval, not with the level clock. This does
**not** meet a fully level-local T02. **T02 remains PARTIAL.**

Exact continuation:

1. Add a `ProductionAmrHydroOptions` per-level source callback, called by
   `advanceProductionAmrHydroRecursiveImpl` after each level sweep with that level's rows and
   `subintervalUpdate`.
2. Implement it in `SourceRuntime` as a row-subset `executeThermochemistry` (staged transaction,
   collective vote per level clock), and skip those rows in the coarse-interval call.
3. Until then, the owner can bound the error by setting `numerics.source_max_fractional_change`.
   Alternatively, add a cooling-time timestep candidate. That would be a new key; not done here.

No duplicate time authority is introduced either way.

## 5. N3: rollback and memory review

* **Coverage.**
  * `GlobalStepRetryPoint` restores the full `SimulationState`: gas, particles, composition,
    populations and event ordinals, the chemistry control sidecar, AMR hierarchy, histories and
    registers.
  * It also restores the integrator, both schedulers, pending output and the gravity force
    endpoint and epochs (generation bumped past the reached value).
  * Since pass 1 it also covers the source-controller checkpoint: counters and ledgers, plus
    invalidation of the ID occupancy cache.
  * HydroAmrRuntime keeps only diagnostics, monotone high-water values and key-validated caches.
    `m_last_predictive_forecast_events` is a diagnostic.
  * Profiler counters deliberately count attempts (`time.global_step_rollbacks`); they are not
    restored, by design.
* **Admission.** The copy is reserved at the accounted baseline before the copy is made. A
  refusal on any rank is collective, before any mutation.
* **Retries.** Collective `kNumericalRetry` or `kPeerRejected` halves Δt and keeps the same
  physics. Any other rejection is rethrown on every rank.
* **Scope.** Enabled only with `numerics.global_step_max_retries > 0`, and only at rung zero.
* **Memory.** The full copy doubles `SimulationState` residency while enabled. This is governed
  and acceptable as the reference. A bounded journal would reduce memory, not fix correctness.

## 6. Tests and benchmarks written in this pass (NOT RUN)

| Name | Content |
|---|---|
| `unit_amr_production_regrid` | Metallicity, vorticity and species indicators; mandatory cooling length and its guard; mandatory truncation reported; missing cooling times refused. |
| `unit_config_chui_part3` | R1/R2 keys round-trip; refusal of Part 2 + cooling length, of species without a network, and of lookahead without regrid or populations. |
| `unit_stellar_population_feedback` | Forecast equals the realized chronology over 60 steps; no mutation; expectation heuristic. |
| `unit_baryonic_source_stage_rollback` | Cooling-time promotion flagged; ratio 0 leaves behaviour unchanged. |
| `unit_chui_part3_remote_sources` | Single-metallicity wind table reports clamps. |
| `unit_chemistry_rate_references` | R18 is secondary-verified and the gate holds. |
| `unit_chui_part3_adaptive_and_decks` | All five `configs/chui_final` decks parse, are recursive, carry checksums, and only the experimental deck has the UVB override. |
| `bench_chui_final_kernels`, `bench_chui_final_mpi` | Executables only (not CTest). |

**Still not written** (owner priority):

* **V1:** 2- and 4-rank versus 1-rank recursive AMR through `ReferenceWorkflowRunner`. This needs
  an AMR initial-condition deck.
* **V2:** association and feedback MPI end-to-end.
* **V3:** injected-retry integration comparing a retried step bit-for-bit against a clean one.
* **N1:** the CTU and dual-energy regression set.

Pass 1 already has loopback-protocol tests for the recursive driver and the distributed
diffusion MPI tests.

## 7. Known compile-risk locations (static only)

* `hydro_amr_runtime.cpp`, `runAutomaticAmrRegrid`:
  * `std::optional<physics::SspModel>` / `StellarPopulationStore` emplacement;
  * designated `CoolingRateQuery{...}` order;
  * `policy.physics.cooling_time_code` span lifetime (the vector lives in the same scope).
* `amr_production_regrid.cpp`: structured binding `[mandatory, patch_ptr]` over pairs, and
  `std::stable_partition` on a pair vector.
* `config.cpp`: the `std::initializer_list<std::pair<const char*, double*>>` parse loop.
* `thermochemistry_operator.cpp`: `m_solver.coolingHeating(input, input.species_cm3, t)`, a const
  member call.
* Pass-1 items listed in `chui_final_design_closure_20261104.md` §7.

## 8. Owner commands

Build the data:

```bash
python3 -m venv .venv-data
```

```bash
.venv-data/bin/pip install h5py numpy
```

```bash
PATH=$PWD/.venv-data/bin:$PATH tools/chui_data/fetch_and_convert.sh "$PWD/run/chui_data_build"
```

Then link the derived files into `<run_dir>/chui_data/` for the `configs/chui_final` decks.

First build and tests:

```bash
cmake --preset cpu-only-debug
```

```bash
cmake --build --preset build-cpu-debug
```

```bash
ctest --preset test-cpu-debug --output-on-failure -L chui_part3
```

```bash
ctest --preset test-cpu-debug --output-on-failure -R 'unit_(amr_production_regrid|config_chui_part3|stellar_population_feedback|baryonic_source_stage_rollback|chemistry_rate_references|chui_part3_adaptive_and_decks|chui_part3_remote_sources)$'
```

After that, run the MPI and benchmark presets from `docs/build_instructions.md`, including
`integration_distributed_element_diffusion_mpi_np{2,4}`,
`integration_amr_family_colocation_mpi_np{2,4}` and `bench_chui_final_*`.

## 9. Readiness

The source is ready for the owner's first compile-and-debug pass. Mandatory items still
**PARTIAL**:

* T02 level-local source coupling (§4);
* H01 CTU interface stencils;
* V02 coordinated controller;
* the V1–V3 end-to-end tests;
* A03 cooling length under Part 2 backends.

Each has an exact continuation above, and unsupported combinations are refused explicitly.
Nothing is reported as compiled, passing, validated or qualified.
