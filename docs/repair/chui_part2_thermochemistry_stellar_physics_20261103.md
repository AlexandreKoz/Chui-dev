# CHUI Part 2 — thermochemistry, composition, dust, ISM/SGS, star formation, stellar populations, feedback

Date: 2026-11-03. Part 1 source baseline: `ef3c7f1` (tree identical to the
merge `75684d6` that this branch started from). Branch:
`claude/dreamy-ritchie-6rm8b1` (the session's designated Part 2 feature branch).
Commits: listed by `git log ef3c7f1..HEAD` (checkpoints A–F);
every commit was pushed to `origin` with ordinary `git push` (no force push,
no merge, no PR).

**Evidence status: SOURCE ONLY. No configure, build, compiler, test,
benchmark, sanitizer or simulation was executed (explicit owner instruction).
Every check below is NOT RUN — USER-RESERVED.** Compilation errors are
expected in a change of this size and must be fixed by the first human build.

Vocabulary: IMPLEMENTED IN SOURCE (connected production path, by inspection
only) · REFERENCE ONLY · EXPERIMENTAL/UNQUALIFIED · BLOCKED BY DATA/TOPOLOGY ·
NOT RUN — USER-RESERVED.


## P2.5 Scientific provenance and physics-data qualification (2026-11-03, from `ed18d36`)

**Still SOURCE ONLY.** Nothing was configured, compiled, built, tested,
benchmarked or simulated. Every test is **NOT RUN — USER-RESERVED**.
Production qualification is **not** declared.

The literature was read from arXiv in this pass, so the corrective-gate
statement that no coefficient could be verified (§0.3) is superseded.
Section §0.3 is retained for history.

### P2.5.1 Rates

Full record: `docs/chemistry_rate_provenance.md`.

- **Primary-verified:**
  - R01–R06: KWH96 Table 2.
  - R07–R09: TREECOOL contract plus Rahmati 2013 eqs. 13–14.
  - R10–R12, R14: GP98 Table 1.
  - R16: GA08 reaction 9.
  - R19: GML07 eqs. 32/34 and DB96 eq. 37.
  - Cooling: KWH96 Table 1 and eq. 23; Wolfire 2003 eqs. 23–24; GP98 A5–A7.
  - Compton: derived from constants.
- **Secondary-verified** (primary inaccessible, reproduced by a compilation):
  - R13, R15, R17 via Abel et al. 1997.
  - R20, R21 via Glover & Clark 2012.
  - HM79 H₂ LTE and HM89 fine structure via KROME.
- **Unverified:** R18 (Jura 1975) and the 0.4 eV H₂ photodissociation heat.
- **Corrections:**
  - R16 restored from GA08.
  - H₂ self-shielding b5 fixed; it was 10× too small.
  - φ_PAH = 0.5 added to photoelectric heating.
  - The Cen 1992 He 2³S term, which is not in KWH96, was removed.
  - O I ΔE set to 227.7 K.
  - R12's 3.2e4 K cap is now documented as GP98's published validity limit.
- **Case A vs case B.** The KWH96 fit is case-A-like. Case B is not
  implemented. Ferland 1992 cross-check functions were added.
- **Gate.** The gate is retained for all three backends. The hybrid and
  detailed manifests are not `qualified()`. For the equilibrium backend,
  provenance is complete but validation is user-reserved.

### P2.5.2 Public interface changes (migration notes)

- `include/cosmosim/physics/chemical_network.hpp`:
  - New `RateVerification::kVerifiedAgainstSecondarySource`, plus
    `rateVerificationQualifies` and `rateVerificationLabel`.
  - **`NetworkManifest::cooling_verification` was removed** and replaced by
    `std::vector<CoolingChannelRecord> cooling_records`. `qualified()` is
    false when the record list is empty.
  - New `equilibriumPrimordialManifest()`.
  - `RateEnvironment::helium_nuclei_cm3`.
  - New rate functions: `recombinationHIICaseA/BFerland1992`,
    `collisionalDissociationH2ByH*`, `h2DissociationCriticalDensity` and
    `h2ThermalDopplerB5`.
  - `k_wolfire2003_phi_pah`.
  - Migration: callers of `cooling_verification` iterate `cooling_records`.
- `include/cosmosim/physics/metal_cooling_table.hpp`:
  - `MetalCoolingTableProvenance::{uv_background, uv_background_treecool_sha256}`.
  - `UvBackgroundCompatibility`, `checkUvBackgroundCompatibility()` and
    `referenceHydrogenMassFraction()`.
  - The `chui_metal_cooling_table_v1` **input** format now *requires*
    `# uv_background = ...` and a `reference_hydrogen_mass_fraction` block.
    Older tables are refused loudly; regenerate them with
    `tools/chui_data/convert_wss09_metal_cooling.py`. No production table
    existed before this pass.
- `include/cosmosim/physics/baryonic_runtime_config.hpp`: new
  `BaryonicDataProducts::uv_background_compatibility`.
- `include/cosmosim/core/config.hpp`: new
  `PhysicsConfig::metal_cooling_accept_uv_background_mismatch`.

### P2.5.3 Consistency audit (§C)

Confirmed defects, each fixed with a regression test:

1. **Metal-table abundance scaling.** It used X_e/X_e,ref, while WSS09 eq. 4
   defines the scaling relative to hydrogen. Now
   (X_e/X_H)/(X_e,ref/X_H,ref). Test: `unit_thermochemistry_backends` in
   `testMetalTable`.
2. **Undeclared UVB in metal tables.** WSS09/HM01 against an FG20 TREECOOL
   could mix silently. A declared UVB is now required, and a mismatch is
   refused unless `metal_cooling_accept_uv_background_mismatch = true`.
   Tests: `testMetalTableUvBackgroundPolicy` and
   `testUvBackgroundMismatchKeyRoundTrip`.
3. **Data files not bound to the run.** A `*_table_path` without a checksum
   recorded only the path in the config hash, so a replaced file could
   change a restarted run silently. Every configured data path now requires
   its 64-hex SHA-256, which is verified at load. Test:
   `unit_config_part2_baryonic_physics`.
4. **Rate defects in §P2.5.1** (R16, b5, φ_PAH, Cen term, O I ΔE). Test:
   `unit_chemistry_rate_references`.

Checked and found consistent (no change):

- **Comoving to physical.** Density is multiplied by a⁻³ at the evaluation
  epoch, volume by a³ in the stager, and isolated modes use a = 1, z = 0.
  Isolated runs therefore evaluate the TREECOOL and metal table at z = 0.
- **Compton.** It is added in `primordial()` and removed exactly for
  non-cosmological modes, with the sign handled.
- **Mean molecular weight.** μ (amu units) is consistent between the
  thermochemistry EOS (`particleDensity`: metals neutral, CO merges C+O)
  and the ISM Jeans length (`k_amu_g`).
  - Risk: the effective multiphase ISM uses its own configured μ with m_p
    (a 0.7% unit difference). It does not consume the sidecar μ.
- **Composition and dust.** Thermochemistry inputs use gas-phase element
  lanes (dust held separately), consistent with the WSS09 scaling by
  gas-phase abundances.
- **Stellar mass return.** Table mode uses the table ejecta (M − M_final)
  and net yields, and does not rebalance. Untracked mass stays in the star
  as a diagnostic, and ejecta are re-capped at the star mass (corrective
  gate).
- **Effective vs resolved ISM.** Ownership is unchanged (regime decision
  staged in `BaryonicSourceStage`). Element-table metal cooling is not
  evaluated for effective-regime cells.

Risk hypotheses (not confirmed defects):

- The Rahmati factor uses σ̄ = 2.49e-18 (HM01, z = 3) for any UVB.
- The same factor is applied to He rates.

### P2.5.4 Configuration, restart and schema

- **New key** `physics.metal_cooling_accept_uv_background_mismatch`
  (default `false`).
  - Typed bool, with a no-op refusal when no table mode is set.
  - Included in the normalized dump and `part2PhysicsKeysAtDefaults`.
  - A restart projection removes the line only at the default value.
- **Reproducibility impact.**
  - The normalized config text, and therefore the config hash, gains one
    line.
  - Part 2 checkpoints written before Part 2.5 restart only with the key
    at its default.
  - Configs that set a data path without a checksum are now refused.
- **No snapshot or restart payload schema change.** The baryonic sidecars
  are unchanged.

### P2.5.5 Data products

Recorded in `tools/chui_data/data_sources.md`. Nothing was committed.

- **TREECOOL FG20.** Both files were obtained from a GitHub mirror; the
  original host was unreachable. SHA-256 values are recorded. The columns
  match the reader. No licence is stated.
- **WSS09.**
  - `z_0.000.hdf5` was obtained from the pyathena mirror; the original
    tarballs returned 403.
  - Its layout and sign convention were inspected.
  - No licence is stated.
- **NuPyCEE/NuGrid.**
  - Licence: BSD-3.
  - Gross isotope yields; Σ yields = M − M_final within 2.1e-4.
  - The untracked net residual reaches 1.66% of the ejecta for massive
    stars, so tolerance ≈ 0.02 must be declared explicitly.
- **Tools** (written, NOT RUN):
  - `validate_treecool.py`;
  - `convert_wss09_metal_cooling.py` (streaming, bounded memory, never
    converts `Metal_free`);
  - `convert_nugrid_yields.py` (net = yield − X0·ejecta, refuses rather
    than rebalances, explicit Ia untracked policy).
- **Blocked:**
  - a Ploeckinger & Schaye 2020 converter (element-resolved semantics and
    shielding convention unspecified);
  - an HM01 TREECOOL that would make WSS09 UVB-compatible (none verified).

### P2.5.6 Tests (NOT RUN — USER-RESERVED)

- New: `unit_chemistry_rate_references`
  (`tests/unit/test_chemistry_rate_references.cpp`, registered in
  `CMakeLists.txt`).
- Extended:
  - `unit_thermochemistry_backends`: metal H-normalization, mandatory
    fields, UVB policy.
  - `unit_config_part2_baryonic_physics`: new key round trip and no-op
    refusal, checksum requirement.
- Expected behavioural shifts in existing tests:
  - Larger b5 means weaker H₂ self-shielding, so photodissociation is
    stronger when G₀ > 0.
  - R16 adds H₂ destruction in warm gas.
  - φ_PAH lowers ψ, raising photoelectric heating efficiency at fixed n_e.
  - Removing the Cen term changes HeII cooling slightly at high density.

  Thresholds in `testDustCatalysedH2Formation` and
  `testDetailedCarbonConservation` may need review on the first run. They
  were not weakened.

### P2.5.7 Compile and validation risks

- New aggregate initializations go through `withStatus(...)`.
- The `CoolingChannelRecord` vector literals.
- `std::all_of` in the config checksum lambda.
- These must compile on the first build.

All remaining numerical validation is listed in
`docs/chemistry_rate_provenance.md` §6.

### P2.5.8 Interpolation-safety correction (after review of `618b276`)

**Status: SOURCE ONLY.** Nothing was compiled or run; every test is NOT RUN
— USER-RESERVED.

- **Confirmed defect.** In `src/physics/metal_cooling_table.cpp`, `locate()`
  accepted any query on a singleton axis. A table converted from WSS09
  `z_0.000.hdf5` alone therefore supplied z = 0 rates at every redshift,
  even under `kReject`.
- **Fix.** A singleton axis now matches only its own coordinate.
  - Tolerance is 1e-8 × max(1, |value|). This is at least 20× the %.9g
    ASCII representation error and the double round-off of log10 T,
    log10 n_H and 1/a − 1, and about 1e6× below a WSS09 grid step.
  - A query on a multi-point axis that lies just past an end, but within
    the same tolerance, is snapped onto the grid.
  - Under `kReject`, out-of-range T, n_H and z throw.
  - Under `kClampDensityAndRedshift`, n_H and z are clamped and
    `Result::clamped` is set. T is never clamped.
  - Clamps are now **reported**: `CoolingHeatingTerms::metal_table_clamped`
    feeds `SourceLedger::metal_table_clamps`, counted once per accepted
    substep and logged by the source stage.
  - Table values are unchanged.
- **Related helper.** `SspYieldTable::massive()` and `snia()` had the same
  silent singleton behaviour. Off-grid queries now set the existing
  `clamped` flag, counted as `table_clamps`. Interpolated values are
  unchanged.
- **Converter.** The WSS09 converter refuses one-redshift output unless
  `--allow-single-redshift` is given, and it writes `redshift_coverage`.
  The restriction is documented in `tools/chui_data/data_sources.md`.
- **Step-subdivision review.**
  - Findings:
    - `massiveReturn` clamped a negative quadrature difference independently
      for the mass and for each element, silently. That breaks element
      closure and invalidates the additivity check whenever it fires.
    - The cumulative return uses age-dependent quadrature nodes, so
      monotonicity is numerical, not structural.
  - Fix: a negative difference is still clamped to zero, but it is now
    counted (`SspModel::nonmonotone_return_clamps`, restored on stage
    rejection and logged). The mass is then re-derived as the sum of the
    clamped elements. Results without a negative difference are unchanged.
  - Test changes:
    - The subdivision test now asserts that no such clamp fired.
    - It checks per-interval element closure and element-wise additivity.
    - It checks that the cumulative return is non-decreasing on a 401-point
      log-age grid for both channels, for the total and for every element.
    - A failure of that monotonicity check would be a genuine quadrature
      finding. The tolerance must not be widened to hide it.
- **Tests.**
  - `unit_thermochemistry_backends`: `testMetalTableSingletonAxes`. It covers
    z = 0 evaluated at z = 3 under reject, plus T, n_H and z singletons,
    clamp reporting and edge snapping.
  - `unit_stellar_population_feedback`: the extended subdivision test and
    singleton yield-axis clamp reporting.

## 0. Corrective gate (2026-11-03, after independent review of `5625b20`)

**Still SOURCE ONLY — nothing compiled or run; every check NOT RUN — USER-RESERVED.**

### 0.1 Class-ownership and integration defects corrected
* `executePart2StellarPopulation()` had been inserted inside
  `DistributedParticleIdRegistry` (it used `SourceRuntimeImpl` state and was
  called from `SourceRuntimeImpl`). The complete implementation now lives in
  `SourceRuntimeImpl`; state ownership, memory admission, collective ordering
  and transaction semantics are unchanged. No globals/friends/duplicated state.
* Static review fixes: one-cell physical dimension in the SF neighbour index
  returned NaN (now 0, missing remote neighbours still NaN); deposited ejecta
  mass now equals the element + dust budget and is re-capped at the star mass
  (gas mass and composition change identically; no mass created by a
  commit-time clamp); SSP diagnostics restored on stage rejection.

### 0.2 Transactional guarantees (local source-stage atomicity)
* New `physics::BaryonicSourceStage` (`baryonic_source_stage.{hpp,cpp}`):
  ISM-regime decision + dwell, chemistry-fidelity decision + dwell,
  thermochemistry, live dust and SGS dissipation are computed from the
  step-start state by PURE kernels (`decideIsmRegime`,
  `ThermochemistryOperator::decideChemistryFidelity/advanceRow`,
  `advanceDustCell` on a private one-row scratch, `sgsDissipate`) and staged
  per row (all composition lanes, e, mu, tags, dwell). SimulationState,
  composition metadata and the published diagnostics change only in
  `commit()` after the StageTransaction collective vote; a rejection leaves
  them bit-identical. Staging is bounded ((lanes + 8) doubles per active row)
  and admitted through the transaction; no full-state copy.
* Stellar-population stage: event ordinals, star masses/ages, population
  records and gas deltas were already staged; the mutable SSP diagnostic
  counters are now restored on rejection.
* NOT covered (Part 3): a rejection after earlier stages of the same step
  (hydro, AMR sweep) committed — global distributed step rollback.

### 0.3 Rate provenance
* The primary sources could not be consulted: the session network policy
  blocks `arxiv.org`. Therefore **no coefficient is verified**. Every reaction
  and cooling fit carries `RateVerification::kUnverifiedTranscription`;
  `equilibrium_primordial`, `fixed_hybrid`, `fixed_detailed` refuse to run
  unless `physics.thermochemistry_allow_unqualified_rates = true` (config
  validation and solver construction). Only the Compton coefficient is derived
  from CODATA constants; `analytic_reference` uses no rate data.
* R16 (H₂ + H → 3H) removed: the implemented form decreased with T above
  1e4 K and vanished by 1e6 K (transcription defect); listed as omitted.
* `docs/chemistry_rate_provenance.md`: per-rate implemented form, claimed
  source, regime risks (e.g. R14 growth above 1e5 K, R12 clamp is an
  implementation choice), owner verification procedure, numerical reference
  values of the implemented formulas, missing channels and consequences.

### 0.4 Stellar yields
* No silent rebalancing: a net-yield residual above the table's declared
  `net_yield_sum_tolerance` (default 1e-3 of the ejecta) rejects the table;
  no renormalization of per-element ejecta; ejected mass = sum of tracked
  element ejecta; untracked residual and negative-ejecta clamps are counted.
* `chabrier_reference_v1` remains REFERENCE ONLY (birth-mixture return), not a
  calibrated enrichment model.
* SSP return is now the difference of cumulative channel returns: exactly
  additive under step subdivision; discrete events were already subdivision-
  and restart-invariant (cumulative-Poisson thresholds in the persistent
  particle-indexed record).

### 0.5 Options and outputs
* `ssp_clustered_associations = true` is now refused (was a no-op).
* **Science snapshots do NOT export composition, species, dust, SGS, ISM tags or
  stellar population records.** Only restart checkpoints carry them (module
  blocks). Checkpoint persistence is not science-snapshot support.

### 0.6 External data requirements (nothing bundled or synthesized)
| Data | Key | Needed for |
|---|---|---|
| TREECOOL UV background | `uv_background_table_path` (+ `_sha256`) | `uv_background_model = treecool_table` |
| element-resolved metal cooling, `chui_metal_cooling_table_v1` | `metal_cooling_table_path` | `metal_cooling_mode = element_table*` |
| SSP yields, `chui_ssp_yield_table_v1` | `ssp_yield_table_path` | `stellar_population_model = chabrier_yield_table_v1` |
| verified rate coefficients (publications) | — | qualifying any network/equilibrium backend |

### 0.7 Tests added (NOT RUN — USER-RESERVED)
`unit_baryonic_source_stage_rollback` (failure injection, bit-identical
preservation), `unit_star_formation_gradients` (patch boundary, periodic,
coarse/fine, missing remote neighbour), plus updates to
`unit_thermochemistry_backends` (unqualified-rate refusal),
`unit_stellar_population_feedback` (exact subdivision additivity, yield table
never rebalanced), `unit_config_part2_baryonic_physics` (opt-in and
association refusal).

### 0.8 Known compilation risks
Large unbuilt change set: new aggregate fields / designated-initializer
orders in tests, `std::optional`/`<string>` includes, mutable members in
`SspModel`, the moved 250-line method, new headers in `cosmosim_physics`,
and `src/`-relative includes in tests (`target_include_directories` added).

## 1. Production call graph (one owner per source term)

```
mode dispatch (four modes unchanged)
 -> TimeCoordinator accepted interval [t_begin, t_end], timeline a_begin/a_end
 -> SourceRuntime (Strang pre half, smooth sources only, physical midpoint dt/4)
      thermochemistry/dust/SGS + metal diffusion for dt/2
 -> HydroAmrRuntime (hydro stage)
      ensureGasCompositionInitialized (fresh init / restart restore by gas_cell_id)
      HydroCoreSolver transport (+ effective relaxation; + legacy cooling ONLY
        when thermochemistry_backend = legacy_cooling)
      HydroFaceMassFluxSink -> CompositionFluxAccumulator (composition advection)
      fixed grid: StageTransaction (hydro + composition commit together)
      AMR: prepare/commit right after the (non-transactional, Part 3) sweep
 -> SourceRuntime (Strang post half; smooth midpoint 3dt/4, events at a_end):
      1. ISM regime manager (sync point)            [ism_regime_model]
      2. local chemistry fidelity selection         [chemistry_adaptive_selection]
      3. thermochemistry StageTransaction (one MPI vote) -> commit e, species, mu
      4. dust (exact accretion/sputtering) + SGS dissipation
      5. refresh P/T/c_s on touched rows (F4 contract) before any consumer
      6. star formation (neighbour-index gradients, composition transfer,
         population record, birth time = t_end)
      7. stellar populations: Part 2 SSP + events + mechanical feedback
         (StageTransaction)  XOR  legacy table evolution + legacy deposition
      8. metal diffusion for dt/2; BH once at t_end
 -> CFL guard of the next step reads refreshed caches
 -> restart: composition packed into module block chui_gas_composition_v1;
    stellar records in particle-indexed chui_stellar_population_v1
```

Cooling is advanced exactly once per accepted interval: in-solver for
`legacy_cooling`, in the source stage for every Part 2 backend (the hydro
runtime does not register `CoolingHeatingSource` then).

## 2. Audit closure

| ID | Status | Where | Tests (NOT RUN) |
|---|---|---|---|
| F2 incomplete cooling interval | IMPLEMENTED IN SOURCE | `cooling_heating.cpp` (implicit BE + step doubling, `completed` flag, unchanged energy on refusal), `CoolingHeatingSource` throws `SourceIntervalIncompleteError`, fixed-grid hydro converts to `kNumericalRetry`; Part 2 solver same contract (`thermochemistry.cpp`), source stage rejects collectively | `unit_cooling_f2_complete_interval`, `unit_thermochemistry_backends` (network F2) |
| F5 newborn star aged a full step | IMPLEMENTED IN SOURCE | `stellarElapsedYearsInInterval`, legacy governed path and Part 2 path skip stars whose birth tick equals the step tick | `unit_stellar_population_feedback` |
| cooling geometry row permutation | CONFIRMED DEFECT, REPAIRED | fixed-grid property provider maps solver geometry rows to dense rows; composition sidecar has a gas_cell_id identity guard | `unit_gas_composition_substrate` |
| SF patch-boundary gradients | REPAIRED | `src/workflows/internal/star_formation_gradients.hpp` (geometric neighbours, periodic, coarse/fine, NaN on missing remote) | (no dedicated test yet — gap) |
| effective EOS table edges | REPAIRED | power-law asymptote above table (flagged), c² = dP/dρ of the interpolant | `unit_ism_dust_sgs_adaptive` |
| relaxation + cooling double count | REPAIRED | relaxation acts only on effective-regime cells; regime tag gates closure/refresh/relaxation/thermochemistry | `unit_ism_dust_sgs_adaptive` |
| high-Mach thermal cancellation | inherited from Part 1 dual-energy design; unchanged | — | — |
| F1/F3/F4/F7 | Part 1 (unchanged) | | |
| F6 runtime AMR / recursive subcycling | Part 3 | | |

## 3. Species, reactions, data

Element vector (11): H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe; A09 solar
pattern; Y_p = 0.2449 (Aver et al. 2015). Species slots: HI, HII, H⁻, H₂, H₂⁺,
HeI, HeII, HeIII (hybrid), + C⁺, CO (detailed, NL97 convention: non-CO C is
C⁺, atomic O derived). Electrons derived by charge neutrality. Reactions R01–R21
with provenance in `chemical_network.cpp` (KWH96, Galli & Palla 1998, Karpas
1979, Abel et al. 1997, Jura 1975, Draine & Bertoldi 1996, Glover & Mac Low
2007, Nelson & Langer 1997). Integer element/charge balance validated at solver
construction. **Coefficients were transcribed from memory of the cited papers
and must be checked against the originals by the owner.** Omitted channels are
listed in each manifest (`omitted_channels`), e.g. HD, three-body H₂, CO
cooling, H⁻ photodetachment.

Cooling/heating: KWH96 primordial terms, Compton (cosmological modes only),
H₂ (HM79 LTE + GP98 low-density), C⁺ 158 µm / O I 63 µm two-level (HM89),
photoelectric (BT94/Wolfire 2003, needs G0>0 and dust), UVB photoheating.

Backends: `analytic_reference` (REFERENCE ONLY), `equilibrium_primordial`
(KWH96 ionization equilibrium), `fixed_hybrid`, `fixed_detailed`
(IMPLEMENTED IN SOURCE). Solver: backward Euler + Newton (FD Jacobian, dense
LU ≤11 unknowns), step-doubling error, Richardson correction, long-double
clock closing the interval exactly, element projection ledgered.

External data — **BLOCKED BY DATA (nothing bundled or synthesized)**:
* UV background: TREECOOL file (`uv_background_table_path`, SHA-256 recorded).
  Without it Part 2 runs must use `uv_background_model = none`; `hm12`/`fg20`
  are legacy analytic proxies and are refused by Part 2 backends.
* Rate coefficients: UNQUALIFIED (see §0.3); network/equilibrium backends need
  the explicit opt-in.
* Hot metal cooling: `chui_metal_cooling_table_v1` converted by the owner from
  a licensed element-resolved table (e.g. Wiersma et al. 2009).
* Production SSP yields: `chui_ssp_yield_table_v1` (AGB/CC net yields vs mass
  and Z, Ia per-event ejecta). Without it only `chabrier_reference_v1`
  (birth-mixture return, REFERENCE ONLY) is available.

## 4. Ownership matrix

| Quantity | Owner | Notes |
|---|---|---|
| gas nuclei per element | composition element lanes (mass) | Z (`metal_mass_code`) is a derived mirror incl. dust metals |
| species | species lanes (mass); sum per element == element lane | projection at solver end |
| dust solids | dust lanes; stoichiometry MgFeSiO4 / C | transfers element-conserving |
| thermal energy | `internal_energy_code` (canonical) | chemical binding energy ledgered, not in e |
| SGS energy | SGS lane (extensive) | only SF/Jeans support, not hydro pressure |
| mu | derived cache (amu), NaN = unknown | rebuilt by chemistry; used by refresh |
| ISM regime / chemistry tag + dwell | composition metadata lanes | persisted, migrated |
| star birth mixture, event thresholds/ordinals, birth time | `chui_stellar_population_v1` | migrates and checkpoints with the particle |

## 5. Transport and topology

* Fixed grid (single rank): composition advected with the final hydro mass
  flux, donor-concentration, step-start donors, periodic wrap; inside the
  hydro StageTransaction.
* Local AMR (single rank): same; composition reflux by fine-face mirroring.
* **Fail closed (Part 3)**: distributed production AMR with composition;
  imported MPI ghosts; AMR regrid with composition (prolong/restrict kernels
  exist: `prolongComposition`/`restrictComposition`); legacy particle-bound
  gas rebuild; mechanical feedback across MPI ranks.
* MPI migration of gas cells carries the composition payload (wire v3/v4).
* Metal diffusion diffuses only total Z; element-resolved diffusion is not
  implemented, so config REFUSES `enable_metal_diffusion` with a Part 2
  backend (fail closed, see §10).

## 6. ISM / SGS
Regimes by Truelove Jeans criterion + effective threshold with hysteresis
and dwell; effective cells skip resolved thermochemistry; closure, refresh and
relaxation honour the tag. SGS `budgeted_v1`: feedback-fraction source, exact
dissipation to heat; no shear production (EXPERIMENTAL/UNQUALIFIED).

## 7. Star formation and chronology
Gradients via geometric neighbour index; SF removes the exact proportional
share of every lane; birth mixture recorded; birth time = t_end; birth tick =
step index; newborn stars are not evolved in their birth interval; event
thresholds are restart-, subdivision- and order-invariant.

## 8. SSP and feedback
Chabrier 2003 IMF, Raiteri 1996 lifetimes, Kalirai 2008 IFMR, NS 1.4 Msun,
direct collapse above `ssp_cc_mass_max_msun`, Ia t⁻¹ DTD with N_Ia = 1.3e-3.
Continuous AGB/CC mass return; discrete CC/Ia energy events. Mechanical
deposition: Hopkins 2018 coupling with KO15 terminal momentum, net vector
removed mass-weighted (radial deficit ledgered), energy closure by momentum
scaling, SGS share, dust condensation (AGB/CC efficiencies) and SN
destruction. Pre-SN momentum is an EXPERIMENTAL user closure (default off).
Associations (`ssp_clustered_associations`) are not implemented and the option
is refused by config validation.

## 9. Config / schema
55 new `physics.*` keys (see `docs/configuration.md` §CHUI Part 2), typed,
validated, dumped; `uv_background_model = treecool_table` added. Restart:
pre-Part-2 checkpoints resume only with all Part 2 keys at defaults (exact
projection in `restart_checkpoint.cpp`). New module blocks
`chui_gas_composition_v1`, `chui_stellar_population_v1` (no change to the
HDF5 restart schema version itself; snapshots do not yet export composition —
gap). Migration wire versions bumped (particle 3, AMR patch 4).

## 10. Known gaps / risks (honest list)
1. Nothing compiled; expect type/signature errors (e.g. new aggregate fields,
   `std::optional` includes, designated-initializer orders).
2. Rate coefficients are UNQUALIFIED (arxiv.org blocked during the gate); owner verification required.
3. Element-resolved metal diffusion not implemented; metal diffusion is
   refused by config when a Part 2 backend is selected.
4. Snapshot HDF5 export of composition/dust/species not implemented.
5. Composition transport is first-order in concentration.
6. (fixed by the gate) regime/fidelity/dwell/dust/SGS are staged and atomic.
7. (fixed by the gate) see §0.2.
8. Association mode refused; benchmark matrices not implemented.
9. (fixed by the gate) dedicated SF gradient test added.
10. Equilibrium backend is expensive (bisection × KWH iteration) — reference use.

## 11. Tests and fixtures (all NOT RUN — USER-RESERVED)
`unit_gas_composition_substrate`, `unit_composition_transport`,
`unit_cooling_f2_complete_interval`, `unit_thermochemistry_backends`,
`unit_ism_dust_sgs_adaptive`, `unit_stellar_population_feedback`,
`unit_config_part2_baryonic_physics`, `unit_baryonic_source_stage_rollback`,
`unit_star_formation_gradients`; updated `unit_migration_wire`.
Fixtures: `configs/chui_part2/{cosmo_cube,zoom_in,isolated_galaxy,isolated_cluster}_part2.param.txt`.

Human commands (documentation only):
```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

## 12. Part 3 obligations
1. **First:** cross-patch / cross-rank periodic AMR for decomposed `cosmo_cube`.
2. Distributed composition ghosts, composition reflux payloads, regrid wiring of
   `prolongComposition`/`restrictComposition`, all-field multilevel transport.
3. Remote feedback deposition and SF neighbour data across ranks.
4. Transactional AMR + global event/source rollback; transactional regime/
   fidelity decisions.
5. Global adaptive chemistry scheduling and statistical shadow qualification.
6. Element-resolved diffusion, snapshot export, benchmarks, four-mode science
   validation, CPU/MPI performance qualification.
