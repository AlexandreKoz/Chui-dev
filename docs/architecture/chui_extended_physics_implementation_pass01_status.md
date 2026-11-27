# CHUI extended physics — implementation pass 01 status ledger

Campaign: CHUI-EXT-PHYS-REF-v1 (charter `chui_extended_physics_bulk_implementation_claude_v01.md`,
master `chui_extended_physics_master_architecture_v01.md`). Date: 2026-11-27.

**SOURCE ONLY.** Nothing in this pass was configured, compiled (including syntax-only), linked, run,
tested, benchmarked or scientifically qualified. Every row below has execution `NOT_RUN` and
qualification `NOT_QUALIFIED`. A row marked `IMPLEMENTED_UNVALIDATED` means the source path exists
end to end as described in its evidence column; it does **not** mean it compiles or is correct.

* Repository: `AlexandreKoz/Chui-dev`; base `main` @ `ada3564d2af3b4219b92bc154b97966ee9cd209a`;
  branch `feature/extended_physics_bulk_pass01`.
* Companion documents: preflight `chui_extended_physics_preflight_pass01.md`; numerical contracts,
  derivations, recalled-lead labelling and per-test expected values
  `chui_extended_physics_numerical_contracts_pass01.md` (cited below as **NC §…**).
* Status vocabulary: `EXISTING_UNQUALIFIED`, `IMPLEMENTED_UNVALIDATED`, `PARTIAL`, `FAIL_CLOSED`,
  `NOT_STARTED`, `EXCLUDED_RESEARCH`.

---

## 1. Requirement counts

Denominator: the rows of §3 (mixed `P0/P1` charter items are split into one P0 row and one P1 row;
test-source items N06/R06/CR08/L04/C04/RA06/O05 are tracked in §4 and not counted here; research
exclusions in §5).

| Tier | Rows | EXISTING_UNQUALIFIED | IMPLEMENTED_UNVALIDATED | PARTIAL | FAIL_CLOSED | NOT_STARTED |
|---|---|---|---|---|---|---|
| P0 | 32 | 2 | 9 | 14 | 6 | 1 |
| P1 | 24 | 0 | 1 | 12 | 6 | 5 |
| Total | 56 | 2 | 10 | 26 | 12 | 6 |

Execution: 56/56 `NOT_RUN`. Qualification: 56/56 `NOT_QUALIFIED`. Research exclusions: 10 (§5).

---

## 2. End-to-end paths

Traced as `param.txt → config → validation → allocation → initialization → stage → kernel → coupling →
output/restart → test registration` (I05).

| Path | State | Trace |
|---|---|---|
| Relativistic background (S01/N01) in any cosmological run | Implemented (unvalidated) | `[cosmology] radiation_background, t_cmb_k, n_ur, ncdm_masses_ev, ncdm_temperature_ratio` → `parseExtendedCosmologyKeys` (`src/core/config_extended_physics.cpp`) → `validateExtendedPhysicsConfig` + flat closure (`src/core/config.cpp`) → `cosmologyBackgroundConfigFromSimulation` (single factory) → `LambdaCdmBackground` at all three construction sites (`config.cpp` integration authority, `time_coordinator.cpp`, `reference_workflow.cpp`) → H(a) in every drift/kick/dt. Deck `configs/extended_physics/neutrino_background_dmo.param.txt`; test `unit_cosmological_background`. |
| Neutrino-aware DMO box, transfer ratio / linear response (N03) | Implemented (unvalidated) | `[neutrinos] method, table_*, response_*` → parse/dump/admission → `GravityRuntimeImpl` builds `neutrinos::NeutrinoSpectralSource::fromConfig` (tables loaded with SHA-256, gauge and cosmology checks) → pointer in `PmSolveOptions::spectral_source_hook` → `pm_solver.cpp` `requireWavenumberDomain` precheck, serial observation, multiplier inside `applyDerivedPoissonOperator` → forces → history stored in sidecar `chui_neutrino_response_v1` after each refreshing PM solve, restored on first execute. Deck `neutrino_linear_response_dmo.param.txt` (table SHA placeholders, must be replaced by real CLASS/CAMB output); test `unit_neutrino_methods`. |
| Stellar RT sources (F02) | Implemented (unvalidated) | `[physics] rt_stellar_sources, rt_stellar_luminosity_table_path/_sha256` → `validateRadiation` (requires table + cosmological + ionizing3) → `RtBhAgnSourceBridge` loads `chui_ssp_group_luminosity_v1` → `buildStellarRadiationSources` (age integral per interval, bucket-located deposits) → `RadiationRuntime` source lane → existing M1/PD-ARS. Test `unit_rhd_reference_backends` (EXT-F02-T01). |
| Isolated galaxy HI cube + direct RIME (RA01, O01, O02) | Implemented (unvalidated) | `[observables] hi_cube…`, `[instruments] …` → admission → `AnalysisRuntime::executeDiagnostics` → `observables::ObservablePipeline` (read-only `GasScienceAccessor`, `LeafGeometry::fromPatches`, `traceStraightRay`, `hiBrightnessTemperature`) → FITS `hi_cube` → instrument coupling (`simulateVisibilitiesDirect`) → `NativeVisibilityWriter` (`chui_visibilities_v1`) or `exportOskarInputs`. Deck `isolated_galaxy_hi_cube_rime.param.txt`; tests `unit_observables_kernels`, `unit_instrument_rime`, `unit_extended_physics_config`. |
| IGM Ly-α forest (L01 snapshot sightlines) | Implemented (unvalidated), partial scope | `[observables] lyman_alpha_forest, forest_*` → pipeline Voigt τ along parallel sightlines → FITS `lya_flux` + summary (`forestStatistics`: mean flux, PDF, 1-D power). |
| tSZ/kSZ/τ_e maps (C03 non-relativistic), DM map (RA05) | Implemented (unvalidated) | `[observables] sz_maps, dispersion_map` → `integrateSz`, `integrateDispersion` → FITS `compton_y`, `tau_e`, `ksz_dt_over_t`, `dispersion_measure`. |
| RT photoionization (F02 reconciliation) | **Refused** | inherited `validateRadiation` refusal retained (no reconciliation with thermochemistry, §6 G-F02-1). |
| CR-MHD gas+CR budget | **Refused** | `cosmic_rays.model ≠ none` refused with exact missing edge (`src/core/config_extended_physics_sections.cpp`). |
| Radio full-Stokes transfer | **No runtime path** | kernels only (`radio_transfer.hpp`); B is not attached to the pipeline accessor (§6 G-RA-1). |
| Neutrino particles / δf (N04) | **Refused** | `neutrinos.method = particles | delta_f` refused (no neutrino particle species). |
| VET / FLD runtime transport (R03), causal rays (R04) | **Refused** | `rt_transport` other than M1 still refused; backends are reference kernels. |
| Accelerated RIME (O03), UVH5/MS (O04) | **Refused / absent** | `instruments.backend = accelerated` refused; no UVH5/MS key exists. |
| Coupled RHD 21 cm, spectral-CR synchrotron, multi-plane lensing driver | **Absent** | dependencies incomplete (RA02, CR runtime, C02 plane builder). |

---

## 3. Requirement rows

Columns: **ID** · **Tier** · **Contract / reference** · **Base evidence** · **Code + call-site evidence** ·
**Config admission** · **State / restart / I/O** · **Memory model** · **Tests** · **Status** ·
**Blocker / next action**. Execution `NOT_RUN` and qualification `NOT_QUALIFIED` apply to every row.

### 3.1 Foundations (§4)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| F01 | P0 | NC §F01: interval workspace formula `8·25·n_values + coupled·3·n_cells·52` | RT runtime unguarded vectors; prior ledger "~5×" estimate | `radiationIntervalWorkspaceBytes`, reservation in `RadiationRuntime::advanceInterval` (`src/radiation/radiation_runtime.cpp`), governor forwarded `SourceRuntimeImpl` → `RtBhAgnSourceBridge::setMemoryGovernor` | none (governor always used when present) | none | `kPhaseResident`, reserve→allocate→commit, RAII release | — (no new test; formula exercised by bench print) | PARTIAL | Level-solver/ghost buffers not in formula; O(P²) topology build; MPI halo/reflux refused (G-F01-1..3) |
| F02 | P0 | NC §F02: `E_g = M_* ∫ L_g(τ,Z) dτ`, external SSP table | `rt_stellar_sources` FAIL_CLOSED | `radiation_stellar_sources.{hpp,cpp}` (`parseStellarLuminosityTable`, `loadStellarLuminosityTable`, `buildStellarRadiationSources`); bridge call in `src/workflows/internal/rt_bh_agn_bridge.cpp` | admitted with table+SHA, cosmological, ionizing3 (`config_rt_bh_agn.cpp`) | none (stateless in age) | table shared read-only; deposits bounded by star count | EXT-F02-T01 | PARTIAL | Stellar sources done; photoionization/thermochemistry reconciliation and transported/external UV partition not done (G-F02-1) |
| F03 | P0 | — | dynamic magnetic/radiation regrid, MHD+RT refused | untouched | refusals retained | — | — | — | FAIL_CLOSED | Out of reach without executing MHD regrid; inherited (G-F03-1) |
| F04 | P0 | — | restart v27 `/state/radiation`; no radiation science export | none for radiation/MHD science fields (N03 state persisted, see N03) | — | no snapshot schema change | — | — | NOT_STARTED | Radiation science-field export + schema bump (G-F04-1) |

### 3.2 Shared contracts (§5)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| S01 | P0 | NC §S01 (photons from T_cmb, massless N_eff, FD eigenstates, strict closure) | `LambdaCdmBackground` with 3 independent construction sites, no photons/ν | `include/cosmosim/core/cosmology.hpp`, `src/core/cosmology.cpp`; factory `cosmologyBackgroundConfigFromSimulation` used at all sites | `[cosmology] radiation_background…`; disabled ⇒ legacy text/hash | config-only (hash) | 4001-point table, constant | EXT-S01-T01..T04 | IMPLEMENTED_UNVALIDATED | Verify closure vs CLASS background (execution campaign) |
| S02 | P0 | NC §S02 (tables, gauge, k units, no extrapolation; provider status) | none | `spectral/{boltzmann_tables,atomic_data,frequency_axis,provider_status,reference_fft}` | via `[neutrinos] table_*`, `[observables] hi_collision_table_*` | provenance strings into products | tables read-only | EXT-S02 cases in `unit_spectral_data` | IMPLEMENTED_UNVALIDATED | Recalled-lead constants to verify against primary sources (G-SCI-1) |
| S03-a | P0 | NC §S03 snapshot rays (periodic/open, proper path, z, ties, parents excluded) | none | `lightcone/ray_geometry.{hpp,cpp}` (`LeafGeometry::fromPatches/uniformGrid`, CSR bucket index, `traceStraightRay`); used by `ObservablePipeline` | implicit through observables | read-only | bucket index bytes `indexBytes()` (not governed, G-S03-2) | EXT-S03-T01..T03 | IMPLEMENTED_UNVALIDATED | Pipeline single rank; `mergeRankSegments` not exercised by an MPI exchange |
| S03-b | P1 | NC §S03 lightcone shells, temporal weights, replication, streaming | none | `lightcone/lightcone_plan.{hpp,cpp}` (`ComovingDistanceTable`, `planLightconeShells`, `temporalWeights`, `replicationMetadata`, `observedRedshift`) | no key (no driver) | no resumable product | — | EXT-S03-T04/T05 | PARTIAL | Streaming lightcone driver, MPI ray migration, HEALPix (G-S03-1) |

### 3.3 Neutrinos (§6)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| N01 | P0 | NC §S01/N01 FD integrals (GL16, asymptotics), validity bound | none | `cosmology.cpp` `massiveNeutrinoDensity/PressureParameter` | `ncdm_masses_ev`, `ncdm_temperature_ratio`; `response_max_mean_velocity_c` validity | — | — | EXT-N01-T01/T02 | IMPLEMENTED_UNVALIDATED | — |
| N02 | P0 | transfer-table consistency | none | `BoltzmannTransferSet` (cosmology/gauge/species checks) | tables admitted only for N03 | — | — | EXT-N02-T01 | PARTIAL | Neutrino species IC adapter not implemented (no species) (G-N02-1) |
| N03 | P0 | NC §N03 ratio (N03a) and derived linear response (N03b); perturbation source, k=0 excluded, no CIC deconvolution of the multiplier | single periodic Poisson operator | `neutrino_linear_response.{hpp,cpp}`; `gravity/pm_spectral_hook.hpp`; `pm_solver.cpp`; `gravity_runtime.cpp` | `neutrinos.method = transfer_ratio | linear_response` admitted (periodic, CPU PM, tables, k-domain ⊇ Nyquist) | sidecar `chui_neutrino_response_v1` (whole block); fingerprint in split-scale lane | `kCanonicalPersistent` `(2+2·k_bins)·8·capacity` (4.06 MiB default) | EXT-N03-T01..T06 | IMPLEMENTED_UNVALIDATED | CUDA PM and isolated PM refused; CLASS P(k) suppression comparison pending |
| N04 | P1 | NC §N04 FD sampling, drift/kick, δf weights | none | `neutrino_particles.{hpp,cpp}` (`sampleThermalMarker`, `driftMarker`, `kickMarker`, `deltaFWeight`, `requireNonNegativeTreeMasses`) | `particles`, `delta_f` refused | — | — | EXT-N04-T01 | FAIL_CLOSED | Neutrino particle species axis (snapshot/restart) (G-N04-1) |
| N05-a | P0 | species auto/cross, total matter | none | `species_power_spectrum.{hpp,cpp}` (`crossPowerSpectrum`, `combineSpeciesSpectra`, `crossSpectraFromRatio`) | no output key | — | per-call mesh | EXT-N05-T01 | PARTIAL | Not called from the analysis runtime (G-N05-1) |
| N05-b | P1 | δf noise, halo profiles | none | `deltaFNoisePower`, `shotNoisePower` | — | — | — | — | PARTIAL | Halo profiles NOT_STARTED |

### 3.4 Radiation hydrodynamics (§7)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| R01 | P0 | existing M1, HLL/GLF, PD-ARS | prior pass (v27) | unchanged except F01/F02 lanes | unchanged | v27 unchanged | F01 | existing | EXISTING_UNQUALIFIED | Inherited gaps of the RT/BH/AGN ledger |
| R02-a | P0 | three ionizing bands | `radiation_spectrum.hpp` ionizing3 | unchanged | unchanged | — | — | existing | EXISTING_UNQUALIFIED | — |
| R02-b | P1 | FUV/LW, X-ray bands with real cross sections | none | none | — | — | — | — | NOT_STARTED | Cross-section providers |
| R03 | P1 | NC §R03 discrete ordinates, normalized Eddington tensor, VET, FLD (LP limiter) | tensor constructor only | `radiation_reference_transport.{hpp,cpp}` (`VetTransport`, FLD CG) | `rt_transport = vet` refused | — | patch-local | EXT-R03-T01..T03 | PARTIAL | Not dispatched in `RadiationRuntime` (G-R03-1) |
| R04 | P1 | NC §R04 causal rays with travel time, ledger | steady rays | `radiation_causal_transport.{hpp,cpp}` `CausalRayTransport` (bounded queue, `serialize/restore`) | no key | serialize exists, no checkpoint lane | bounded queue capacity | EXT-R04-T01 | PARTIAL | Direct/diffuse handoff, runtime stage (G-R04-1) |
| R05 | P1 | HG Monte Carlo reference, counter RNG | isotropic elastic MC | `HgMonteCarloSettings`, HG sampler | no key | — | bounded packets | EXT-R05-T01 | PARTIAL | AP/implicit MC FAIL_CLOSED (no complete discrete contract) |

### 3.5 Cosmic rays (§8)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| CR01 | P0 | NC §CR grey p/e, P=e/3 | none | `cosmic_ray_transport` `advectGrey`, `CrGasExchange` | `cosmic_rays.model` refused | none | lazy (no allocation when off) | EXT-CR01-T01 | FAIL_CLOSED | MHD level-stepper stage (G-CR-1) |
| CR02-a | P0 | anisotropic diffusion (Sharma–Hammett limiter), streaming | none | `anisotropicDiffusionRate`, `diffuseExplicit`, `streamGrey` | refused | — | patch | EXT-CR02-T01/T02/T04 | FAIL_CLOSED | as CR01 |
| CR02-b | P1 | super-time-stepping | none | `diffuseRkl2` via `mhd::makeRkl2Coefficients` | refused | — | — | EXT-CR02-T03 | FAIL_CLOSED | implicit diffusion NOT_STARTED |
| CR03 | P1 | staggered two-moment, exact relaxation | none | `advanceTwoMoment` | refused | — | — | EXT-CR03-T01 | FAIL_CLOSED | as CR01 |
| CR04-a | P0 | fixed momentum bins | none | `cosmic_ray_spectrum` `CrSpectralState`, `advanceMomentumSpace` | refused | — | `bytesFor = lanes·bins·cells·8` | EXT-CR04-T02 | FAIL_CLOSED | as CR01 |
| CR04-b | P1 | piecewise power law; spectral two-moment | none | `reconstructPowerLaw` | refused | — | 2 lanes | EXT-CR04-T01 | PARTIAL | Spectral two-moment NOT_STARTED (no verified discrete contract) |
| CR05-a | P0 | injection normalization, electron sync+IC with S01 u_cmb, table-only proton losses | none | `cosmic_ray_sources` `injectCosmicRays`, `electronRadiativeLossCoefficient`, `CrLossProvider` | refused | — | — | EXT-CR05-T01 | FAIL_CLOSED | Feedback budget debit edge |
| CR05-b | P1 | secondaries, γ, CR ionization | none | providers return UNAVAILABLE | refused | — | — | — | FAIL_CLOSED | Versioned rate data |
| CR06 | P1 | Eulerian electron spectra; tracer spectra | none | CR04+CR05 kernels | refused | — | — | EXT-CR06-T01 | PARTIAL | Tracer spectra NOT_STARTED |
| CR07-a | P0 | memory formula, no group dropping | none | `CrSpectralState::bytesFor`, governor reserve | refused | none | see §7 | — | FAIL_CLOSED | AMR/MPI flux registers |
| CR07-b | P1 | AMR/MPI registers, regrid, restart, RNG | none | none | — | — | — | — | NOT_STARTED | after CR runtime stage |

### 3.6 Ly-α, CMB, radio (§9–§11)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| L01 | P0 | NC §L01 Voigt (Weideman N=32), Hubble velocity axis, ordered segments | none | `line_transfer` `lymanAlphaOpticalDepth`, `segmentAveragedProfile`, `forestStatistics`; pipeline | `[observables] lyman_alpha_forest` admitted | FITS + summary | product vectors (ungoverned, G-OBS-1) | EXT-L01-T01..T03 | PARTIAL | FGPA comparator, 3-D flux statistics, lightcone sightlines NOT_STARTED |
| L02 | P1 | NC §L02 resonant MC, exact two-component atom-velocity sampler, recoil, dust | none | `resonant_lyman_alpha.{hpp,cpp}`, `_sampler.hpp` | no key | — | bounded packets | EXT-L02-T01 | PARTIAL | No runtime/CLI consumer |
| L03 | P1 | scattering-rate & momentum estimators | none | `ResonantResult::scattering_rate_per_cell`, `momentum_rate_per_cell` | no key; pipeline uses x_α = 0 (collisional-only model, declared) | — | — | — | PARTIAL | Dynamic coupling FAIL_CLOSED (no source/rollback contract) |
| C01 | P0 | external C_ℓ tables, S01 T_cmb | none | `cmb_secondaries` `parseCmbSpectraTable`, `clTT`, `gaussianRealizationFlatSky` | no key | — | flat-sky map | EXT-C01-T01 | PARTIAL | TE/EE realization (G-C-1) |
| C02-a | P0 | Born convergence | none | `rayShootLensPlanes` (single plane = Born) | no key | — | — | EXT-C02-T01 | PARTIAL | Plane-construction driver from snapshots (G-C-2) |
| C02-b | P1 | multi-plane with lens–lens coupling | none | `rayShootLensPlanes` | no key | — | — | EXT-C02-T01 | PARTIAL | as C02-a |
| C03-a | P0 | NR tSZ y, kSZ, τ_e; ISW kernel | none | `integrateSz` (pipeline), `integrateIsw` (kernel) | `sz_maps` admitted | FITS | — | EXT-C03-T01 | IMPLEMENTED_UNVALIDATED | ISW/RS potential history not wired (kernel only) |
| C03-b | P1 | relativistic tSZ | none | `tszSpectralFactor` NR only | refused (no key) | — | — | — | FAIL_CLOSED | Verified relativistic kernel |
| RA01 | P0 | NC §RA01 full τ HI transfer, spin temperature, ν₁₀ = 1420.40575 MHz | none | `hiBrightnessTemperature`, `spinTemperature`, `parseCollisionalCouplingTable`; pipeline | `hi_cube`, `hi_spin_model`, `hi_transfer` admitted | FITS `hi_cube` | cube in memory (G-OBS-1) | EXT-RA01-T01/T02 | IMPLEMENTED_UNVALIDATED | Ly-α colour coupling input absent (RA02) |
| RA02 | P1 | cosmic dawn coupling | none | none | — | — | — | — | NOT_STARTED | Ly-α background and X-ray heating providers |
| RA03-a | P0 | free-free (Gaunt), synchrotron power law, Faraday | none | `radio_transfer` `freeFree*`, `synchrotronPowerLaw*`, `faradayRotationRate` | no key | — | — | EXT-RA03-T01/T02 | PARTIAL | B not attached to accessor (G-RA-1) |
| RA03-b | P1 | synchrotron from CR spectra | none | `synchrotronFromSpectrum`, `synchrotronKernelF/G` | no key | — | — | EXT-RA03-T02 | PARTIAL | CR runtime |
| RA04 | P1 | Stokes `dS/ds = j − KS`, Padé(6,6) augmented exponential (singular K safe) | none | `stokesStep`, `propagationMatrix` | no key | — | — | EXT-RA04-T01 | PARTIAL | No runtime consumer |
| RA05-a | P0 | DM/RM with redshift factors | none | `integrateDispersion` (pipeline DM) | `dispersion_map` admitted | FITS | — | EXT-RA05-T01 | PARTIAL | RM requires B (G-RA-1) |
| RA05-b | P1 | Stokes maps, spectral index, FRB | none | none | — | — | — | — | NOT_STARTED | after RA04 runtime |

### 3.7 Instruments (§12)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| O01 | P0 | FITS/WCS products; OSKAR input export | none | `fits_writer`, `exportOskarInputs` | `instruments.backend = oskar_export` admitted | FITS, OSKAR dir | — | EXT-O01-T01 | PARTIAL | HEALPix, pyuvsim/CASA/WSClean adapters |
| O02 | P1 | NC §O02 direct full-pol RIME, Jones chain G·D·Ionosphere·E, w-term, smearing, noise, flags | none | `interferometry.{hpp,cpp}` `simulateVisibilitiesDirect`, `NativeVisibilityWriter`; pipeline coupling | `[instruments]` admitted (direct_rime) | `chui_visibilities_v1` streamed per time | batch per time step | EXT-O05-T01..T06 | IMPLEMENTED_UNVALIDATED | Foreground product separation not represented (signal-only sky) |
| O03 | P1 | accelerated backends | none | none | `accelerated` refused | — | — | — | FAIL_CLOSED | Choose dependency (GPU direct / NUFFT) |
| O04 | P1 | UVH5 / MS | none | none | no key | — | — | — | NOT_STARTED | Adapter dependency |

### 3.8 Cross-cutting (§13)

| ID | Tier | Contract | Base | Code / call site | Config | State / restart | Memory | Tests | Status | Blocker / next |
|---|---|---|---|---|---|---|---|---|---|---|
| I01 | P0 | stage placement | existing runtimes | N03 in gravity runtime; F02 in source bridge; observables in analysis runtime | — | — | — | — | PARTIAL | No new registry stage entries; CR stage absent |
| I02 | P0 | single `.param.txt` path, inert-key rejection, disabled = unchanged hash | registry | `configKeyRegistry` + `extendedPhysicsKeyDefaults`, `validateExtendedPhysicsConfig`, dumps | all sections | config hash | — | `unit_extended_physics_config` | IMPLEMENTED_UNVALIDATED | — |
| I03 | P0 | provenance in products; single authority | — | product headers carry provenance, table SHA, seeds, axes, `NOT_QUALIFIED` | — | sidecar for N03 | — | — | PARTIAL | Input commit not embedded in products (config provenance only) |
| I04 | P0 | governed allocations | `MemoryGovernor` | radiation workspace, N03 history, CR spectral state reserved | — | — | §7 | bench | PARTIAL | Observable products and ray index ungoverned (G-OBS-1, G-S03-2) |
| I05 | P0 | vertical paths | — | §2 | decks | — | — | decks + unit tests | PARTIAL | RT photoionization, CR-MHD, radio Stokes paths refused/absent |
| I06 | P0 | docs + CMake registration | — | 6 new libraries, 9 tests, 1 bench registered; docs updated | — | — | — | — | IMPLEMENTED_UNVALIDATED | Not configured or built |

---

## 4. Verification sources (all NOT_RUN)

Registered in `CMakeLists.txt` (foreach `cosmosim_ext_phys_unit`, labels `ext_phys;unit`,
`COSMOSIM_SOURCE_DIR` defined): `unit_cosmological_background`, `unit_spectral_data`,
`unit_neutrino_methods`, `unit_lightcone_geometry`, `unit_rhd_reference_backends`,
`unit_cosmic_ray_transport`, `unit_observables_kernels`, `unit_instrument_rime`,
`unit_extended_physics_config`. Benchmark: `bench_extended_physics_kernels`. Gate mapping:
`docs/validation_plan.md`. Expected values: NC per-section "Tests" lists.

| Item | Coverage authored | Not authored |
|---|---|---|
| N06 | FD limits, table gauge rejection, ratio/response analytic history, sidecar round trip, δf weights | IC/restart/rank equivalence, PM/tree force split integration |
| R06 | ordinates/VET/FLD diffusion limit, causal ledger, HG moments, stellar age integral | Strömgren with stellar sources, MPI/regrid, RT-MHD |
| CR08 | advection, anisotropic Gaussian, RKL2, streaming, two-moment diffusive limit, spectral adiabatic, cooling | CR shock, coupled MHD, restart/regrid, no-CR equivalence run |
| L04 | slab τ, Voigt wings, resonant slab escape | packet convergence vs RASCAS/COLT, Nyx statistics |
| C04 | uniform slab y/τ, Born point lens | relativistic tSZ, Agora lightcones |
| RA06 | thin/thick HI, free-free, synchrotron closed forms, Stokes rotation | full-Stokes radio end-to-end |
| O05 | phase, polarization, smearing, beams, Jones, noise, flags, polar uvw | accelerated vs direct, external fixtures |
| V00 | disabled defaults dump nothing; inert keys rejected | full-suite disabled regression run |

---

## 5. Research exclusions (EXCLUDED_RESEARCH, no keys created)

Adaptive neutrino phase-space switching; spectral CR AMR; adaptive radiation group changes; adaptive
angular closure switching; hybrid Ly-α MC/diffusion; observable-driven AMR; unified adaptive fidelity
controller; native Einstein–Boltzmann/GR engine; global kinetic/PIC; full receiver hardware.

---

## 6. Gap list (severity P0 critical … P3 minor)

Kind: **D** confirmed source defect/omission, **R** inferred risk, **U** unavailable evidence, **X**
research exclusion. Origin: **inh** inherited from base, **int** introduced or left by this pass.

| Gap | Sev | Req | Kind/Origin | Evidence | Missing edge | Current guard | Fix | Test |
|---|---|---|---|---|---|---|---|---|
| G-BUILD-1 | P0 | all | U/int | whole diff | never compiled; 14k lines unverified | none (branch unmerged) | execution campaign: configure, build, fix | full suite |
| G-F02-1 | P1 | F02/R01 | D/inh | `config_rt_bh_agn.cpp::validateRadiation` photoionization refusal | thermochemistry reconciliation, UV partition | refusal retained | derive H/He photo-rates coupling with rollback | Strömgren ladder |
| G-F01-1 | P1 | F01 | D/inh | RT topology build O(P²) in `radiation_runtime.cpp` | patch spatial index | none (cost only) | reuse `LeafGeometry`-style bucket index | topology bench |
| G-F01-2 | P1 | F01 | D/inh | multi-rank RT refusal (A14) | owned-face halo, synchronized reflux | refusal | MPI halo exchange with epochs | 2-rank equivalence |
| G-F01-3 | P2 | F01 | D/int | `radiationIntervalWorkspaceBytes` | level-solver ghost buffers excluded | formula undercount | extend formula | governor rejection test |
| G-F03-1 | P1 | F03 | D/inh | MHD regrid refusals | CT prolongation/EMF reflux | refusal | dedicated MHD campaign | — |
| G-F04-1 | P1 | F04 | D/inh | snapshot schema 6/7 | radiation science fields | none | schema bump + readers | schema test |
| G-N02-1 | P1 | N02/N04 | D/inh | `core::ParticleSpecies` lacks neutrinos | species axis in snapshot/restart | `particles/delta_f` refused | species axis migration project | IC round trip |
| G-N04-1 | P1 | N04 | D/int | `requireNonNegativeTreeMasses` | signed-multipole tree | refusal | PM-only scale-limited path first | force split test |
| G-N05-1 | P2 | N05 | D/int | `species_power_spectrum` unused at runtime | analysis output key | none | add output key + dump | output schema test |
| G-N03-1 | P2 | N03 | R/int | `pm_solver.cpp` hook | CUDA PM / isolated PM | refusal | port observation to CUDA | — |
| G-N03-2 | P2 | N03 | R/int | linear-response replay after restart | replay semantics depend on deterministic PM refresh cadence | sidecar method/binning checks | verify by restart equivalence run | restart test |
| G-S03-1 | P2 | S03 | D/int | `lightcone_plan.cpp` | streaming driver, HEALPix, ray migration | no key | implement driver in analysis runtime | lightcone resume test |
| G-S03-2 | P3 | S03/I04 | D/int | `LeafGeometry` bucket index | governor reservation | none | reserve `indexBytes()` | — |
| G-OBS-1 | P2 | I04/L01/RA01 | D/int | `observable_pipeline.cpp` product buffers | governor reservation (HI cube = n_side²·channels·8 B) | none | reserve before allocation | governor rejection test |
| G-OBS-2 | P2 | I05 | D/int | pipeline single-rank sightlines | MPI segment exchange (`mergeRankSegments`) | admission refuses `parallel.mpi_ranks_expected > 1` (`config_extended_physics_sections.cpp`) | collective segment gather | 2-rank equivalence |
| G-RA-1 | P2 | RA03/RA05 | D/int | `GasScienceAccessor` B unavailable | CT B sample attachment | RM/Faraday not produced | attach MHD face-centred B | RM screen test |
| G-R03-1 | P2 | R03 | D/int | `VetTransport` unused by `RadiationRuntime` | runtime dispatch, MPI angular comm | `rt_transport` refusal | dispatch behind unqualified mode | beam-crossing test |
| G-R04-1 | P2 | R04 | D/int | `CausalRayTransport` | direct/diffuse handoff, checkpoint lane | no key | ledgered handoff | budget test |
| G-CR-1 | P1 | CR01–CR07 | D/int | `config_extended_physics_sections.cpp` | MHD stage, AMR/halo/restart lanes, budget debit | admission refusal | CR stage in MHD level stepper | coupled CR-MHD budget |
| G-C-1 | P3 | C01 | D/int | `gaussianRealizationFlatSky` TT only | TE/EE correlated realization | — | Cholesky per ℓ | spectrum recovery |
| G-C-2 | P2 | C02 | D/int | `rayShootLensPlanes` | plane construction from snapshots | no key | driver | convergence test |
| G-SCI-1 | P1 | S02/L/RA/CR/R | U/int | NC "recalled lead" labels: Lyα f and A, A₁₀ = 2.85e-15 s⁻¹, Gaunt fit, RL synchrotron forms, LP limiter, HG, TMS/Sharma–Hammett, Weideman, Field spin, AHB, Elbers response | citations recalled, not verified against primary sources | labelled in NC | verify each against the paper/data source | value tests already authored |
| G-INH-1 | P3 | — | R/inh | stellar_evolution default-background overload | second background instance in stellar evolution | none | route through factory | — |
| G-INT-1 | P3 | F02 | R/int | `RtBhAgnSourceBridge` reference member | bridge not move-assignable | none | pointer member | — |
| G-HYG-1 | P3 | docs | R/int | copied spec markdown | trailing whitespace in verbatim copies (`git diff --check`) | intentional (unchanged copies) | none | — |

---

## 7. Memory models (formulas, not measurements)

| Owner | Class | Formula | 128³ | 256³ | 512³ |
|---|---|---|---|---|---|
| Radiation canonical state (ionizing3) | canonical | 4·3·N·8 | 0.19 GiB | 1.50 GiB | 12.0 GiB |
| Radiation interval workspace (coupled) | phase-resident | 25·3N·8 + 3N·52 | 1.48 GiB | 11.81 GiB | 94.5 GiB |
| CR spectral, 2 species × 16 bins × 2 lanes | canonical | 2·16·2·N·8 | 1.00 GiB | 8.00 GiB | 64.0 GiB |
| CR spectral two-moment (8 scalars, not implemented) | — | 2·16·8·N·8 | 4.00 GiB | 32.0 GiB | 256 GiB |
| N03 response history (64 bins, capacity 4096) | canonical | (2+2·64)·8·4096 | 4.06 MiB | 4.06 MiB | 4.06 MiB |

N = uniform-equivalent leaf cells per rank. Against the 42 GiB hard / 36 GiB nominal profile, a coupled
256³ RT interval workspace alone needs ≈ 13.3 GiB together with its canonical state; 512³ does not fit on
one rank. No feasibility claim is made. Groups, precision and models are never reduced to fit: the
reservation refuses before allocation.

---

## 8. Config, schema, migration notes

* New sections `[neutrinos]`, `[cosmic_rays]`, `[observables]`, `[instruments]`, new `[cosmology]` keys and
  `rt_stellar_luminosity_table_path/_sha256` (`docs/configuration.md`). All disabled by default; disabled
  sections emit nothing in the normalized dump, so config hashes and restart admission of existing runs are
  unchanged. Non-default values on a disabled model are rejected (inert-key rule).
* `rt_stellar_sources = true` changes from unconditional refusal to admission with explicit requirements
  (migration: previously refused decks now need the table keys).
* No snapshot or restart schema version changes. N03 state uses the existing module-sidecar block.
* New product schemas: FITS products, `chui_observables_summary_v1`, `chui_visibilities_v1`
  (`docs/output_schema.md`).
* Reproducibility impact: disabled path unchanged; N03 deterministic per rank topology; RIME noise from
  counter-based seeds.

## 9. New files and changed interfaces (index)

* New public headers: `core/config_extended_physics.hpp`; `gravity/pm_spectral_hook.hpp`;
  `spectral/{provider_status,boltzmann_tables,atomic_data,frequency_axis,reference_fft}.hpp`;
  `neutrinos/{neutrino_linear_response,neutrino_particles,species_power_spectrum}.hpp`;
  `lightcone/{ray_geometry,lightcone_plan}.hpp`;
  `radiation/{radiation_stellar_sources,radiation_reference_transport,radiation_causal_transport}.hpp`;
  `cosmic_rays/{cosmic_ray_transport,cosmic_ray_spectrum,cosmic_ray_sources}.hpp`;
  `observables/{gas_science_view,line_transfer,resonant_lyman_alpha,resonant_lyman_alpha_sampler,cmb_secondaries,radio_transfer,fits_writer,observable_pipeline}.hpp`;
  `instruments/interferometry.hpp`.
* Changed public interfaces (additive): `core/cosmology.hpp` (`RadiationSpeciesConfig`, new background
  methods); `core/config.hpp` (`CosmologyConfig::radiation_species`, `SimulationConfig::extended`,
  `cosmologyBackgroundConfigFromSimulation`, `validateExtendedPhysicsConfig`); `core/constants.hpp`
  (physical constants); `core/config_rt_bh_agn.hpp` (stellar table keys); `gravity/pm_solver.hpp`
  (`PmSolveOptions::spectral_source_hook`, default null); `radiation/radiation_runtime.hpp`
  (`governor` setting, `radiationIntervalWorkspaceBytes`). Migration: none needed for existing callers
  (all additions default to the previous behaviour).
* New libraries: `cosmosim_spectral`, `cosmosim_lightcone`, `cosmosim_neutrinos`, `cosmosim_cosmic_rays`,
  `cosmosim_instruments`, `cosmosim_observables` (dependency direction in `decision_log.md`).

## 10. Commits

| SHA | Slice |
|---|---|
| 87e6058 | docs: preflight, unchanged master/charter/audit-cycle copies |
| 1a8fc54 | S01/N01 background authority |
| 2ace5bc | S02 spectral data |
| 4a2854f | N03–N05 neutrinos |
| 21786ee | S03 lightcone |
| 2eed4bf | F01/F02, R03–R05 |
| 07ace24 | CR01–CR07 |
| bfad43c | L/C/RA observables and pipeline |
| 7c508dc | O01–O05 instruments |
| 2dde9d7 | decks, config test, bench, configuration docs |
| (this commit) | schema/restart/validation/profiling/decision docs, status, this ledger |

Unresolved design decisions: neutrino particle species axis (N04); CR stage placement inside the MHD level
stepper; whether observables move to a standalone analysis CLI for multi-snapshot lightcones.
