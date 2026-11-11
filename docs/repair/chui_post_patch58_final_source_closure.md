# CHUÍ — post-patch-58 final source-implementation closure ledger

Date: 2026-11-11. Baseline: `f2f9fc4` (merge of PR #136, patch 56/58 closure),
branch `claude/eager-galileo-92u2k2`, clean working tree at start. This ledger
extends, and does not delete, `chui_post_patch56_v10_final_source_closure.md`
and `chui_post_patch55_source_closure.md`. Where a status below differs from
those ledgers, **this ledger supersedes it** (marked "supersedes").

Normative inputs read in full: `00_PROMPT.md` (patch 58 closure pack),
`chui_patch58_independent_source_audit.md`,
`chui_final_implementation_route_consolidation_v10.md`,
`chui_full_physics_memory_architecture_specification_v02.md`,
`chui_data_staging.md`; repository `AGENTS.md`, the patch-56 ledger and the
module docs touched below.

**Executed:** source reading, text search, manual call-path and
reader/writer review, `git status`, `git diff --stat`, `git diff --check`.
**NOT RUN (owner request):** CMake configure, compile, link, CTest, Python
tools, static analyzers, MPI, HDF5 round trips, simulations, benchmarks,
numerical or scientific validation. Every status below is static source
evidence only. Nothing is committed.

Status vocabulary (unchanged): `implemented_static_unverified`,
`repaired_in_this_patch`, `open_source_gap`, `blocked_external_data`,
`conditional_requires_contract`, `optional_research_deferred`,
`not_in_scope`, plus `open_owner_science_decision`.

---

## 1. P1-GEOMETRY — exact AMR leaf-tiling proof — repaired_in_this_patch

**Defect (confirmed, audit P1).** `globalAmrHierarchyDomain`
(`src/workflows/hydro_amr_runtime.cpp`, patch 56) compared the allreduced sum
of populated-patch volumes with the bounding-box volume. Counterexample
(x-slabs `[0,.5) [.4,.7) [.8,1)`): equal volumes, yet an overlap and a hole.
Malformed populated patches (NaN/inf/non-positive extent) contributed volume
0 and were otherwise skipped, and `level < 0` populated patches were silently
ignored. Disjointness is *not* independently enforced elsewhere:
`hasProductionAmrHydroCoverage` proves row/patch-ID ownership (each row names
one patch), not geometric disjointness; refine/derefine construct children by
halving (`refineProductionPatchInSimulationState`), but migration
(`commitAmrPatchMigration`), restart restoration and per-rank root
materialization (`materializeRootHydroPatchIfMissing`, patch ID = rank + 1 for
non-generated ICs) can each deliver geometry that no construction invariant
re-checks. Option (ii) of the prompt was therefore implemented.

**Mathematical basis (what establishes non-overlap and coverage).** Every
populated leaf is mapped to a half-open box `[a,b)` on the integer lattice
whose unit is the globally finest populated cell width per axis (faces must
lie on lattice planes within `1e-6 + 64·eps·span` lattice units; each cell
width must be an integer number of lattice units, i.e. mixed levels must be
index/dyadically aligned; span ≤ 2^40 per axis). By inclusion–exclusion

    1_[a,b)(x) = Σ_{c ∈ corners} s(c) H(x − c),  H(y) = Π_axis [y_axis ≥ 0],
    s(c) = (−1)^(number of axes where c takes b).

The coverage multiplicity `m = Σ_P 1_P` is therefore the image under `H` of
the summed signed corner measure `W`. Orthant functions of distinct lattice
points are linearly independent (evaluate at a support point minimal in the
product order: only its own term is non-zero), hence

    m ≡ 1_B  ⇔  W = W_B (the 8 signed corners of the bounding box B).

`m ≡ 1_B` is exactly "every lattice cell of B covered once, nothing outside
covered": complete coverage **and** non-overlap, proved in integer arithmetic,
for any number of levels, ranks and patch shapes. Corners shared by
neighbouring leaves cancel; a valid tiling aggregates to exactly 8 points. The
summed lattice volume is used only to *classify* a failure (overlap / hole /
overlap compensating a hole) in the diagnostic, never as proof.

**Implementation.**

| File | Content |
|---|---|
| `include/cosmosim/amr/amr_geometry_coverage.hpp`, `src/amr/amr_geometry_coverage.cpp` (new, MPI-free) | `AmrPatchGeometryRecord`, `diagnoseAmrPatchGeometry` (level ≥ 0, cell count = checked dims product, finite origin, finite positive extent, finite upper bound strictly above origin — catches absorption overflow — finite positive width; message names patch ID, rank, level, axis, value), `isLegacyUndeclaredAmrPatchGeometry`, `makeAmrCoverageLattice`, `amrCoverageSnapTolerance`, `appendAmrPatchSignedCorners` (representability, containment, integer cell width), `aggregateAmrSignedCorners` (sort/merge/drop zero, checked int64), `verifyAmrExactTiling`, `verifyAmrPatchSetTiling` (serial convenience). |
| `src/workflows/internal/amr_geometry_assurance.{hpp,cpp}` (new) | `AmrGeometryAssurance::assure` collective protocol and fingerprint cache; `AmrGeometryAssurancePolicy` (explicit-geometry requirement, periodic box). |
| `src/workflows/hydro_amr_runtime.cpp` | `globalAmrHierarchyDomain(state, amr_options, phase)` now delegates to the assurance (policy from `physical_boundary_kind` / `periodic_box_extent_comoving`); called (a) in `runProductionAmrHydroPath` before any ghost exchange/fill, (b) in `runAutomaticAmrRegrid` right after the regrid application vote. Telemetry: `geometry_proof_runs`, `geometry_proof_global_leaves`, `geometry_proof_gathered_corner_points` in `hydro.amr_production_stage`. |
| `src/workflows/source_runtime.cpp` | `SourceRuntimeImpl::execute`: assurance before the Strang pre-half / unsplit dispatch (that half runs before the hydro stage, so on the first step after IC materialization, restart restore or migration the source-stage patch-topology consumers — element diffusion, SF gradients, remote neighbour import — would otherwise read an unproven hierarchy). Legacy all-zero rows are non-AMR here (lenient policy). |
| `CMakeLists.txt` | new sources in `cosmosim_amr` / `cosmosim_workflows`; test `unit_amr_geometry_coverage`. |

**Entry points and invariant survival.** IC/root tiling: proved at the first
hydro step / pre-half source stage (cache empty). Refine/derefine: proved
after every applied automatic regrid and again (cached) before the next
ghost fill. Migration (`migration_balance_runtime.cpp`
`commitAmrPatchMigration`, ownership compaction): changes the local populated
geometry fingerprint on both ranks → full proof before the next ghost fill or
pre-half source stage. Restart restore: new process, empty cache → full proof
before the first ghost fill. Synchronized and recursive hydro: both enter
through `runProductionAmrHydroPath`; recursive level clocks run inside a
step whose hierarchy was proved and that cannot regrid mid-step (regrid is at
the synchronized hierarchy boundary). The proof is per step boundary, never
per face/ghost. Fail-closed missing interior deep/edge CTU stencils, CTU
validity checks and the flux-face set are untouched (no diagonal faces).

**Collective MPI behaviour.** C1 `allreduceMaxUint64` (cache gate) → V1 vote
`rethrowCollectivePreparationFailure` (malformed leaves) → C2 populated count
→ C3 nine `allreduceMinDouble` (box and finest width) → V2 vote
(representability, periodic box, local aggregation) → C4 residual count and
volume → V3 vote (gather admission) → C5 bounded `gatherBytesToRoot` (round
-limited transport) → root verdict (any root exception becomes the verdict
text) → C6 `broadcastBytesFromRoot` → identical `std::runtime_error` on every
rank. No rank throws before a collective its peers are in; ranks owning no
patch participate in every collective. Cache validity is only ever set after
a collectively accepted proof, so it is rank-uniform.

**Memory ownership and complexity.** Local scratch `P_local × (64 B record +
8 × 32 B corners)` reserved as `kScratchArena`
(`hydro.amr.geometry_proof.local_scratch`) before allocation; root receive
`2 × 32 B × R_gathered` as `kCommunication`
(`hydro.amr.geometry_proof.root_gather`); all released at return
(reservation declared before the buffers it covers). Work O(P_local log
P_local) per rank, O(R log R) at the root. After local cancellation a rank
keeps only the corners of its own owned region (surface-vertex proportional),
so `R_gathered ≤ 8 P_total` worst case and typically ∝ Σ rank-region boundary
vertices; never ∝ N_cells, and no rank ever holds all patch records or a
dense finest-level bitmap. Steady state (unchanged hierarchy): one O(P_local)
FNV fingerprint + one scalar collective per call.

**Supported / rejected operating modes.** Supported: rectangular,
hole-free, non-overlapping leaf tilings at any mix of levels and ranks, legal
fine/coarse adjacencies, zero-patch ranks, inactive (cell_count = 0) parents
and metadata-only patches (skipped), isolated sub-box hydro domains inside a
larger gravity/PM box, periodic runs whose leaves tile `[0, L)` on every
periodic axis (box convention of `amr/internal/amr_periodic_topology.hpp`).
Rejected collectively with a diagnostic: holes, overlaps (incl. compensating
ones), non-rectangular unions, malformed/non-finite/zero/negative geometry,
negative levels, non-lattice-aligned faces or non-integer width ratios,
populated legacy all-zero descriptors on the production AMR path, periodic
sub-box domains. No physical open boundary is inferred from a hole; no
non-rectangular root-domain boundary model exists in source, so none is
admitted.

**Config / restart / schema.** No new key; no restart, snapshot or provenance
schema change; restart v25 (pending flux group 2, temporal history group 3)
and AMR migration wire v5 unchanged. Reproducibility: physics of valid
hierarchies unchanged; invalid hierarchies previously accepted by the volume
sum (or failing later at a ghost) now fail at the step boundary.

**Unexecuted verification (NOT RUN):** `tests/unit/test_amr_geometry_coverage.cpp`
cases 1–9 (see `docs/validation_plan.md`, patch 58); MPI fault-injection runs
listed there.

## 2. P2-S02 — explicit, enforceable yield-coverage policy — repaired_in_this_patch (science approval OPEN)

**Defect.** Coverage was reported (`channelCoverage`,
`sources.data_assumptions`) but the table mode always edge-clamped silently
in model terms (default CC range 8–40 M☉ vs staged NuGrid 1–25 M☉; W7 Ia
single-Z).

**Trace.** AGB channel integrates `[imf_min, min(imf_max, cc_min)]`, CC
`[max(imf_min, cc_min), min(imf_max, cc_max)]` (`cumulativeMassiveReturn`,
64 fixed log-mass nodes from the turnoff node upward, birth metallicity of
the star's record); table queries happen only in `SspModel::ejectaElements`
(`SspYieldTable::massive`) and `iaReturnPerEvent` (`snia`, metallicity-only);
both run inside the stellar-population `StageTransaction::prepare`
(`SourceRuntimeImpl::executePart2StellarPopulation`) whose exceptions are
voted collectively and recorded as a rejection; commit (masses, ages, records,
event ordinals, cursors, deposition) runs only after a unanimous accept. The
predictive-AMR `SspModel` uses no table. Counters are runtime diagnostics
(rollback checkpoint), not restart state.

**Implementation.**

* `include/cosmosim/core/config.hpp` / `src/core/config.cpp`: typed key
  `physics.ssp_yield_table_coverage = reject | edge_clamp_research`
  (default `reject`): parse, key registry default, validation (a research
  value without `chabrier_yield_table_v1` is refused), normalized dump,
  `sspYieldTableCoverageConfigToString`.
* `include/cosmosim/physics/stellar_population.hpp`,
  `src/physics/stellar_population.cpp`: `SspYieldCoveragePolicy`,
  `SspConfig::yield_coverage`, `k_chui_ssp_imf_identity =
  "chabrier2003_individual"`, `SspYieldChannelCoverage` /
  `SspYieldCoverageAssessment` / `assessSspYieldCoverage`,
  `sspYieldCoveragePolicyName`, `SspYieldTable::massive(…, mass_clamped,
  metallicity_clamped)` overload, `metallicityCount` / `massCount`,
  `SspModel::table_mass_clamps` / `table_metallicity_clamps`.
  Construction (strict): refuse if the AGB upper end or the whole CC range is
  outside the tabulated mass axis, a mass axis is a singleton, or any channel
  (AGB/CC/Ia) is single-metallicity. At use (strict): any clamped query throws
  `SourceModelUnavailableError` naming channel, mass, Z and the grid. Research:
  edge row, counted per query and per axis. Ia has no mass-domain check.
* **Demonstrated source fix in the same trace:** clamp detection now treats a
  query within 1e-12 (relative) of an axis end as on-grid. The fixed nodes are
  `pow(10, log10 m)`, a few ulp outside exact channel limits; without this the
  strict gate would refuse an exactly covering table, and the research
  counters over-counted. Interpolated values are unchanged.
* `src/workflows/source_runtime.cpp`: policy mapped into `SspConfig`;
  split counters in the rollback checkpoint, in the stage-rejection restore
  (the wind clamp counter, previously not restored on rejection, is now
  restored too) and in `sources.stellar_population`; `sources.data_assumptions`
  records `ssp_yield_table_coverage`,
  `ssp_yield_table_coverage_research_approximation`,
  `ssp_yield_strict_construction_refusal`, per-channel tabulated/configured
  ranges, `_single_metallicity`, `_mass_range_clamped`,
  `_lower_mass_enforced_at_use` (existing key names kept).
* `src/io/restart_checkpoint.cpp`: key added to `part3PhysicsKeysAtDefaults`;
  pre-patch-58 projection accepts only the default.

**Startup vs runtime.** Construction-time checks cover what every run queries
(CC range, AGB upper end) and what can never be covered (single-Z channel).
The AGB lower end (reached after ~lifetime(m_min)) and the actual birth
metallicities are known only at use and are enforced there; the strict gate
fires up to one quadrature node (~0.03 dex) early (conservative). No claim of
coverage is derived from current stellar metallicities.

**MPI / rollback.** Construction refusal: SourceRuntime construction is
agreed collectively ("runtime module construction 'sources'"). At-use
refusal: inside `prepare`, voted (`stellar population budgets`) → transaction
rejected on every rank → nothing commits; counters restored.

**Restart / config.** New normalized key ⇒ new config hash/provenance. A
pre-patch-58 checkpoint resumes only through the legacy projection, which
requires table physics disabled and the key at default; a checkpoint whose
run used table yields is refused (never silently re-interpreted).

**Accepted / rejected combinations.** `legacy_table` / `chabrier_reference_v1`
+ default: unaffected. `chabrier_yield_table_v1` + `reject`: admitted only
for a table covering CC and AGB-upper ranges with multi-Z channels; refused at
use outside the grid. `chabrier_yield_table_v1` + `edge_clamp_research`:
admitted, named, counted. Research value without table mode: config error.
**Staged NuGrid/W7 data:** refused by `reject` (CC 25–40 M☉ uncovered; W7
single Z; NuGrid AGB lower bound 1 M☉ would also refuse at use late in a
cosmological run); requires `edge_clamp_research` — added explicitly to the
four `configs/chui_final/*` decks that use it and to
`profiles/SAFE_DEFAULT_AVAILABLE_DATA.param.fragment` (+ generator).

**OPEN owner science decisions:** (a) is CC edge-row use for 25–40 M☉
acceptable, or must the CC upper mass / data change; (b) is a
single-metallicity Ia (W7) table acceptable (a per-axis strict exemption
would be a scientific decision and is not implemented); (c) is the AGB
1 M☉ lower bound acceptable for long runs.

## 3. P2-F02 — machine-checkable wind IMF/metallicity contract — repaired_in_this_patch (science approval OPEN)

**Defect.** `chui_ssp_wind_table_v1` carries the IMF only as free text; the
staged SB99 table is Kroupa (0.1–0.5–100), single Z = 0.014, while the CHUI
SSP is Chabrier.

**Implementation.**

* `include/cosmosim/physics/stellar_wind.hpp`, `src/physics/stellar_wind.cpp`:
  `chui_ssp_wind_table_v2` = v1 data layout + required header declarations,
  each exactly once: `source`, `license`, `imf` (token `[a-z0-9_]+`,
  `unknown` allowed and never compatible), `imf_mass_min_msun`,
  `imf_mass_max_msun` (finite, positive, max > min, full-string numeric
  parse), `normalization = per_msun_initial_stellar_mass`,
  `energy_unit = erg_per_msun`, `momentum_unit = msun_km_s_per_msun`,
  `metallicity_definition = metal_mass_fraction`. Duplicates, malformed values,
  two schema headers, an unknown schema version, or a v2 key inside a v1 table
  refuse the table. v1 loading is unchanged otherwise (IMF undeclared; free
  text never parsed as metadata; `# license` tolerated). SHA-256 binding
  unchanged; the verified hash is now kept (`sha256()`).
  `StellarWindTableMetadata`, `StellarWindCompatibilityPolicy`,
  `StellarWindCompatibilityReport`, `assessStellarWindCompatibility`,
  `stellarWindCompatibilityPolicyName`, `metallicityMin/Max`.
* Config: `physics.ssp_wind_table_compatibility = strict |
  approximate_research` (default `strict`), same parse/registry/validation/
  dump/restart-projection pattern; research value requires
  `ssp_wind_model = table_v1`.
* `SourceRuntimeImpl` constructor: assessment against
  `k_chui_ssp_imf_identity` and the configured `ssp_imf_mass_{min,max}_msun`;
  strict refuses v1 (undeclared), `unknown`/mismatched IMF, mismatched IMF
  limits (per-unit-mass normalization differs) and single-Z tables. At use
  (strict) an out-of-grid birth metallicity throws inside the stellar
  population transaction; research counts the clamp as before.
  `sources.data_assumptions` gains `ssp_wind_table_sha256`, `_schema`, `_imf`
  (`undeclared` for v1), `_imf_mass_msun`, `_normalization`,
  `_metallicity_range`, `ssp_imf_identity`, `ssp_wind_table_compatibility`,
  `_imf_matches_ssp`, `_imf_limits_match_ssp`, `_accepted_approximations`.
* `tools/chui_data/convert_sb99_winds.py`: `--schema v2` with owner-stated
  `--imf-identity`, `--imf-mass-min-msun`, `--imf-mass-max-msun` (not inferred;
  default `v1` keeps the pinned bytes/SHA-256 of the staged product).

**Ownership / double counting.** Unchanged: the wind channel adds energy and
momentum only; mass stays in the SSP CC return; `fb_presn_wind_momentum`
remains mutually exclusive. No IMF reweighting, no multi-Z synthesis.

**Accepted / rejected combinations.** Strict: v2 + `chabrier2003_individual`
+ matching limits + ≥ 2 metallicities, in-grid queries. Research: anything
loadable, approximations named and recorded. **Staged SB99 v1 table:**
refused by strict (undeclared IMF; single Z); requires
`approximate_research` (decks and SAFE_DEFAULT fragment updated). A v2
re-conversion of it must declare `kroupa2001`-type identity with 0.1–100 limits
and would still be refused by strict (IMF mismatch, single Z) — correctly.

**OPEN owner science decision:** whether the Kroupa single-Z SB99 wind budget
is acceptable with the Chabrier SSP for any production claim. Data need:
Chabrier-consistent multi-Z v2 tables.

## 4. External data contract table (supersedes the S02/F02 rows of patch 56)

| Product | Gate in source (patch 58) | Strict status | Execution status |
|---|---|---|---|
| S02 NuGrid+W7 (`e2161646…f8e14`) | `ssp_yield_table_coverage` construction + at-use gate | refused (CC 25–40 uncovered, single-Z Ia) | research only (`edge_clamp_research`) |
| F02 SB99 v1 (`ba232a1c…5da5`) | `ssp_wind_table_compatibility` + v2 metadata contract | refused (IMF undeclared, single Z) | research only (`approximate_research`) |
| C01 WSS09, FG20 UVB, C02, HM01, D01 | unchanged (patch 56 table) | unchanged | unchanged |

Blocked external data (unchanged list, now machine-gated for S02/F02):
full-range multi-Z yields; multi-Z Ia; Chabrier-consistent multi-Z wind v2
tables; HM01-consistent UVB/cooling; verified C02 CI/RR/DR + photo-rates;
grain-collision rates.

## 5. H01–V02 matrix (supersedes the patch-56 matrix where marked)

| req | entry point(s) | owner | MPI/AMR/time/restart path | status | remaining mandatory code | external data / decisions | deferred verification (NOT RUN) |
|---|---|---|---|---|---|---|---|
| H01 | `HydroCoreSolver::advance`, `hancockHalfStepIncrement`, lane-generic HLL, `fillAmrGhost` | canonical gas SoA | **supersedes:** exact leaf-tiling proof before ghost fill | repaired_in_this_patch (geometry proof); else implemented_static_unverified | none found | — | convergence, positivity, MPI geometry fault injection |
| H02 | `selectDualEnergyThermalState`, ghost aux | gas sidecar | restart v25 unchanged | implemented_static_unverified | none | — | Mach sweep, restart round trip |
| C01 | `ThermochemistrySolver`, `MetalCoolingTable` | composition + gas sidecar | source transaction | blocked_external_data | none | HM01-consistent UVB/table | equilibrium/transient histories |
| C02 | `metal_ion_network.cpp` | chemistry sidecar | same | blocked_external_data | none | verified rates | solver accuracy |
| C03 | adaptive selection | regime tag + history | restart/migration | implemented_static_unverified | none | qualification | matched-error comparisons |
| C04 | chemistry shadow | diagnostic | V02 publication | implemented_static_unverified; shadow hydro optional_research_deferred | none mandatory | shadow-hydro contract | disagreement qualification |
| I01 | `effective_multiphase_ism.cpp` | ISM config | closure/CFL refresh | implemented_static_unverified | none | calibration | no double cooling |
| I02 | `ism_dust_sgs.cpp` | SGS lane | regrid/restart | implemented_static_unverified | none | coefficients | energy transfer |
| S01 | `star_formation.cpp`, gradients | star sidecar | pre-half source stage now runs on a proven hierarchy | implemented_static_unverified | none | — | patch/rank independence |
| S02 | `SspModel`, `assessSspYieldCoverage`, at-use gate | star sidecar + cursors | transaction rejection; config key + restart projection | **supersedes:** repaired_in_this_patch (policy); open_owner_science_decision (edge clamp, single-Z Ia) | none | full-range multi-Z yields | §2 fixtures; staged-data refusal run |
| F01 | `supernova_feedback.cpp`, thermodynamic refresh | gas + ledgers | cross-rank deposition | implemented_static_unverified | none | — | conservation/MPI |
| F02 | `StellarWindTable` v1/v2, `assessStellarWindCompatibility` | event state | construction gate + at-use gate | **supersedes:** repaired_in_this_patch (contract); open_owner_science_decision (Kroupa/single-Z) | none | Chabrier multi-Z v2 tables | §3 fixtures |
| D01 | composition/diffusion/dust | element/dust inventories | composition reflux/migration | implemented_static_unverified; conditional_requires_contract (shattering/coagulation) | none mandatory | grain-collision law + rates | partition tests |
| A01 | `amr_production_regrid.cpp`, refine/derefine | patch SoA + identity | **supersedes:** post-regrid exact tiling proof | implemented_static_unverified (+ proof repaired_in_this_patch) | none | — | evolving-refinement runs |
| A02 | recursive driver, orchestrator | temporal history + registers | proof before every ghost fill epoch | implemented_static_unverified (+ proof) | none | — | MPI diagonal/periodic, reflux |
| A03 | regrid indicators | planner | — | implemented_static_unverified | none | force-resolution criterion (G01) | indicator qualification |
| A04 | predictive/cost-aware regrid | planner | memory admission | implemented_static_unverified | none | — | fixed-budget comparisons |
| G01 | `tree_pm_coupling.cpp` | particles + force cache | restart force cache | conditional_requires_contract | fine-level correction only if a force-error contract demands it | force-resolution contract | force at fine levels |
| T01 | `time_coordinator.cpp`, `stage_transaction.cpp` | integrator + scheduler | collective rollback | implemented_static_unverified | none | — | temporal convergence |
| T02 | level smooth sources, event barrier | level clocks + cursors | global event barriers | implemented_static_unverified; level-local events / IMEX/SDC optional_research_deferred | none mandatory | — | accuracy vs cost |
| M01 | governor, node admission | governor | geometry proof scratch/gather leases (new) | implemented_static_unverified | none found | site reserve policy | RSS/high-water |
| M02 | core + parallel | — | — | implemented_static_unverified | none | — | serial/OpenMP/MPI builds |
| V01 | restart, ledgers, `sources.data_assumptions` | restart/provenance | two new normalized keys; v25/wire v5 unchanged | implemented_static_unverified (+ S02/F02 provenance) | none | — | four-mode end-to-end |
| V02 | adaptive resource control | controller sidecar | restart sidecar | implemented_static_unverified | none | — | adaptive-policy comparisons |
| X01 | — | — | — | not_in_scope | — | — | — |

## 6. Interface, config, schema and reproducibility impact (migration notes)

* **New public header** `include/cosmosim/amr/amr_geometry_coverage.hpp`
  (MPI-free; `core`-independent; dependency direction respected).
* **Additive public changes:** `stellar_population.hpp` (policy enum,
  `SspConfig::yield_coverage` appended last — aggregate initializers
  unaffected; `k_chui_ssp_imf_identity`; assessment types; `massive` overload;
  axis counts; split counters), `stellar_wind.hpp` (metadata/compat types,
  accessors), `config.hpp` (two enums, two `PhysicsConfig` fields appended,
  two `toString`). No existing signature changed.
* **Behaviour changes:** (1) S02/F02 strict defaults refuse data that patch 56
  edge-clamped/admitted — explicit opt-ins restore the previous numerics, now
  named and recorded; (2) wind-table header recognition is an exact
  (whitespace/CR-trimmed) match instead of a prefix match, and v2-reserved
  keys are refused in v1 files; (3) yield clamp classification within 1e-12 of
  an axis end; (4) malformed/overlapping AMR hierarchies refused at the step
  boundary.
* **Config:** two new typed keys with validation, normalized dump, docs
  (`docs/configuration.md`), tests; no second config system; `.param.txt`
  only.
* **Schema:** no snapshot/restart/provenance schema change; restart v25 with
  pending flux group 2 and temporal history group 3; AMR wire v5. Wind table
  file schema v2 added (reader, docs, converter, fixtures in one change).
* **Reproducibility:** normalized configuration (hash) gains two lines for
  every run; runs using table yields/winds must now state the research
  policies; numerics under the research policies equal patch 56 except the
  clamp classification counters.

## 7. Remaining obligations

* **Mandatory source gaps:** none identified for H01–V02 after this patch
  (static evidence only).
* **Owner science decisions (OPEN):** §2 (a)–(c), §3.
* **External data:** §4.
* **Conditional / deferred research (not failed mandatory fixes):** D01 grain
  shattering/coagulation; G01 multilevel gravity correction; C04 shadow hydro;
  T02 independent fine-level events, IMEX/SDC; X01 MHD, general RT, new AGN.
* **Known limitation (documented, not a gap):** non-rectangular root hydro
  domains have no boundary model and are refused; the strict yield gate may
  fire one quadrature node early at the AGB lower bound.
* **NOT RUN:** compilation, unit/integration tests (incl.
  `unit_amr_geometry_coverage`, `unit_ssp_wind_data_contracts`,
  `unit_config_chui_part3`, `unit_chui_part3_adaptive_and_decks`), MPI, HDF5
  restart round trips, simulations, numerical validation, benchmarks.

**Not `SOURCE-COMPLETE` in any validated sense:** no compile or test was run.
