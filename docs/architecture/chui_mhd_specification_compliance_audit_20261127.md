# CHUÍ MHD — Specification compliance audit, pass 04 (2026-11-27)

**Status:** SOURCE_ONLY / NOT_CERTIFIED. Archive `Chui-dev-main(13)(1).zip`; predecessor `Chui-dev-main(7)(2).zip`. Reference documents: the repository's CHUI-MHD-REF-v1, CHUI-MHD-PERF-v1, development compatibility amendment, previous pass-02 ledger and 2026-11-26 independent audit. No repository changes, build, test, benchmark or simulation. Archive has no verified live commit identity. Focused source audit, not a mathematical proof or complete audit of the new extended-physics project.

## Executive decision

The new archive is an extended-physics campaign, not a magnetic corrective campaign. Full byte-for-byte comparison of common files confirms 74 new, 29 changed, zero removed files. All 44 selected MHD kernels, headers, tests, benchmarks, MHD bridge, hydro admission/dispatch and magnetic I/O files are byte-identical. The shared CMake/config/gravity/source/scheduling/analysis files did change. None of the 15 findings of the preceding magnetic audit has evidence of closure. Keep the old IDs and do not repeat the bulk implementation. New cross-module concerns below must be closed or explicitly gated before combined physics can be admitted.

## Newly identified or expanded concerns

| ID | Severity | Evidence and implication | Required source closure |
|---|---|---|---|
| M04-01 | P0, conditional on combined linear-response neutrinos + retry | `src/core/config_extended_physics.cpp:363–477` admits linear response without a magnetic exclusion; `runtime_capabilities.cpp:435–438` requires global retries for MHD. `time_coordinator.cpp:1193–1220,1475–1496` copies/restores state and calls gravity transaction rollback. `gravity_runtime.cpp:2370–2373` restores live neutrino source only on first use; `:2430–2436` writes response to canonical sidecar after refresh. `neutrino_linear_response.cpp:452–479` appends live history. `gravity_runtime.cpp:642–704` rollback restores forces and cadence but not neutrino history; `tree_pm_coupling.cpp:1994–2000` invalidates tree source maps only. Restored canonical sidecar and live history can diverge after a rejected interval containing a PM refresh; subsequent history can contain rejected or nonmonotonic time nodes. | Make the live spectral hook transaction-aware, including current multiplier, ratio, history append state, initialization status and derived cache fingerprint. Restore from the opening accepted state or a bounded journal, including the empty-history case (`restoreFrom` returns false on absence and does not reset). Charge all copies. Temporarily refuse MHD+linear_response until this path is complete; preserve unrelated nonmagnetic configurations while recording the general retry defect. Add forced post-refresh rejection/retry, direct shortened trajectory, restart replay and no-history cases. |
| M04-02 | P0 inherited composition risk now exposed through observables | `config_extended_physics_sections.cpp:337–358` admits observables with species-resolved thermochemistry; `hydro_amr_runtime.cpp:1086–1150` MHD branch still returns before composition advection (`:1359ff`). New `observable_pipeline.cpp:108` samples the gas composition through `GasScienceAccessor`. Even nonmagnetic-dependent HI, DM, Ly-alpha and SZ products may reflect incorrectly transported species in a magnetized run. | Close M03-02 using accepted MHD mass fluxes. Until then refuse any MHD run with active unsupported composition (including observables and loaded restart state). Do not merely refuse RM; this affects species-dependent products. Add manufactured passive species transport plus independent budget and derived observable regression. |
| M04-03 | P1 transaction / resource risk, not an executed failure | `analysis_runtime.cpp:39–43,46–55,176–181` constructs the pipeline and runs it inside AnalysisHooks, while the admission estimator only considers existing diagnostics. `observable_pipeline.cpp:251–273,387–390,411–419` publishes science files/visibility products/summary inside that call; no magnetic transaction is passed. `time_coordinator.cpp:1482–1496` accepts/rejects the global interval later. Thus observable allocations/publication are outside the explicit MHD retry/memory contract; a later rejection or partial output failure needs an atomic policy. | Make optional observables run after accepted synchronization or stage products transactionally and publish only on commit. Account peak product/geometry/FFT buffers and enforce configured actual MPI size, not expected size alone. Keep combined MHD+observables gated until composition, accepted-state timing, resource and output contracts are established. No broad rewrite of the new observable kernels. |
| M04-04 | P1 integration boundary, deliberately absent | `gas_science_view.hpp:54–57` defines optional physical-Tesla B attachment; `gas_science_view.cpp:132–150` refuses unavailable magnetic samples; `observable_pipeline.cpp:108,325` never attaches CT fields and explicitly reports RM unavailable. New CR kernels depend on magnetic geometry/RKL2 but `config_extended_physics_sections.cpp:294–305` explicitly rejects all enabled CR models. `config_rt_bh_agn.cpp:240–241` still refuses RT+MHD. | Preserve these refusals. Define a read-only accepted-state CT sampler for later radio/CR work; derive physical B with declared unit conversion and a^-2, bounded views and stable gas row mapping. Do not enable RM, synchrotron, CR backreaction, RT-MHD or AGN magnetic injection on the strength of standalone kernels. If the new sampler is implemented, test uniform-field conversion and regrid/ownership freshness; no new physical model is required in this pass. |
| M04-05 | P1 cross-cosmology verification | New extended background is admitted by `config_extended_physics.cpp:325–360`; current MHD source uses existing timeline (unchanged numerical kernels), and gravity still branches to magnetized kick at `gravity_runtime.cpp:3092–3099`. No new MHD-specific test was added. The new background/PM hook changes the integration environment, not the MHD equations. | Check all scale-factor/H(a) providers use the single configured background. Add homogeneous expansion and coupled wave/gravity source cases with radiation background and background-only neutrinos, including legacy-LambdaCDM equality and MHD-off hash checks. Treat as an unvalidated combination, not proof of a wrong H(a) formula. Linear response is conditional on M04-01. |

## Evidence and scope qualification

The new extended-physics ledger itself lists CR-MHD as refused, missing magnetic observable sampling (G-RA-1), and MHD regrid as inherited (G-F03-1). The previous MHD pass-02 ledger still explicitly admits the HLLD/PPM UCT reduction limitation (A12) and unverified Balsara/Toth–Roe equivalence (A13). The current CMake MHD registrations moved to 1226–1262; only the hierarchy test has the explicit `mhd` label in that block. No new MHD test was added. M03 source locations in unchanged files remain valid. Shared-file line numbers must be re-resolved by symbol before patching; the CMake position is updated here. M04-03 is a lifecycle risk requiring a regression rather than a demonstrated rejected-output event. No claim that every other extended-physics feature is sound.

## Inherited open findings (M03 IDs retained)

| ID | Severity | Evidence / source finding | Required correction and regression |
|---|---|---|---|
| M03-01 | P0 | `src/workflows/internal/mhd_runtime_bridge.cpp:159–182,531–598`: local patch/row hash passed as `hierarchy_epoch`. `src/mhd/mhd_level_stepper.cpp:382,415` serializes and compares it with the receiver's local value. Ordinary distinct rank partitions therefore have incompatible wire epochs. | Use deterministic global descriptor hash/generation common to all ranks (separate local cache key); validate collectively; production 2/4-rank unequal partitions, zero-patch rank, intentional stale epoch. REF 7.5, G06. |
| M03-02 | P0 | `src/workflows/hydro_amr_runtime.cpp:1086,1110–1150`: initializes composition then returns from MHD branch; legacy `CompositionFluxAccumulator` path at 1359ff and AMR counterpart are bypassed. `MhdRuntimeBridge::advanceOptions` at 148–156 disables boundary gas integrals; MHD staging contains six conserved lanes. MHD admission does not reject active composition. Species/element/dust mass therefore lacks matching conservative advection with MHD mass flux. | Route accepted per-face mass fluxes through composition transport, including coarse/fine/MPI and SSPRK2 stages, preserve per-cell and global budgets. Until integrated, explicitly refuse any active composition-transport combination, including on restart. REF 6,9 and existing composition contract. |
| M03-03 | P0 | `src/io/restart_checkpoint.cpp:4031–4046,4611–4617` picks magnetic schema by rank-local nonempty magnetic arrays. Zero-gas rank selects v25 in otherwise magnetized v26 run. `src/workflows/output_restart_runtime.cpp:918–922` repeats local selection for metadata; collective science output already explicitly accommodates zero-gas ranks at `snapshot_hdf5.cpp:1091ff`. | Run-wide schema decision from normalized config/collectively agreed capability; empty valid v26 group with provider provenance on zero-leaf ranks, reader/verifier and restart-set manifest agree. Test real distributed restart with empty rank; also verify RT v27 metadata, avoid silent schema collision. REF 13, G08. |
| M03-04 | P0/P1 | `mhd_runtime_bridge.cpp:410–462`: magnetic layout and seed-edge scratch allocated and state committed before MHD stepper governor admission; `:185–239,531–621` directory allgather, row maps, stepper topology allocation precede corresponding leases; `mhd_level_stepper.cpp:2045ff` explicitly charges already-built topology. `mhd_runtime_bridge.cpp:124–125` table provider constructed before table admission. | Preflight checked bytes before allocations, transaction candidate until collective acceptance, charge topology, row maps, seed scratch, provider and replacement coexistence; memory rejection tests at each seam. PERF 5,8,9, MP06. |
| M03-05 | P1 | `mhd_level_stepper.cpp:349–461,1980–2034`: packs entire per-peer phase in double vectors, then byte vectors, materializes receive byte and decoded double vectors. `planMemory` refuses over-budget phases; no MHD credit-based subdivision into planes/edge ranges. Bounded underlying MPI transport alone is not bounded application staging. | Implement compact per-direction chunk/credit pool, include retained capacities, reductions and simultaneous packing/unpacking in cap; scale fragmented layouts; MP03/09. |
| M03-06 | P1 | `mhd_runtime_bridge.cpp:547–559` retains six gas lanes per local cell, stepper `planMemory:1890ff` retains multiple full-level copies; time coordinator `1152–1225` separately captures complete global state for mandatory retry. REF/PERF permits necessary trial and rollback state but prohibits duplicate gas container at rest and requires peak coexistence accounting. | Charge all overlapping images; avoid retaining staging across inactive intervals or formally reconcile lifetime; eliminate redundant copies where transaction guarantees permit. MP02/04/05/07. |
| M03-07 | P1 | REF 7.2/7.3, PERF 9: `runtime_capabilities.cpp:439ff` rejects dynamic regrid and multi-rank rebalance. Transfer helpers exist, but `amr_production_regrid` / migration payloads do not integrate magnetic authority. | Implement governed gas+field refine/coarsen and migration, ownership/epoch updates and collective failure rollback; G07/I08/MP06. Keep guards until then. |
| M03-08 | P1 numerical | Pass-02 ledger A12 explicitly admits UCT-HLL exact 1D reduction only for HLL face flux with PLM-MC states, not the primary HLLD/PPM profile; induction implementation unchanged from pass01. | Derive/implement an upwind CT construction consistent with actual face flux and reconstruction or approved equivalent; test oblique, axis and strict 1D consistency with HLLD/PPM. REF 3,6. |
| M03-09 | P1 numerical | Pass-02 A13 does not derive Balsara/Tóth–Roe equivalence. Dynamic regrid not live. | Verify divergence-and-curl constraints, anisotropic rectangular patches, limited energy-compatible fallback; do not relabel algorithm as certified. REF 7.2. |
| M03-10 | P1 | `seedFreshField` builds field independently per patch and checks local divergence; `integrateEdgePotential` uses three-point quadrature, stochastic seed integrates each edge separately. Quadrature over a coarse edge is not generally identical to sum over fine edges. Entry consistency is strict (`mhd_level_stepper.cpp:1042ff`). | Construct common composite edge potential and restrict additive line integrals before curl, or exact analytic integrals; test nonuniform stochastic seed on two levels. Also imported face-flux IC and magnetized shock preset are missing REF 10 requirements. |
| M03-11 | P1 | `src/core/time_integration.cpp:111–119` introduces magnetic-domain owner validation, but `kMagneticFaces` does not appear in any actual registered stage contracts beyond validator and enum. Initialization also mutates magnetic truth from timestep candidate. | Register actual magnetic inputs/outputs/mutation and read-only consumers; initialize before accepted timeline/output/rollback, with collective admission. Preserve MHD-off contracts. |
| M03-12 | P1 | `MhdRuntimeBridge::advanceOptions:148–156` disables production boundary flux integrals. Snapshot v7 has derived B, magnetic pressure and diagnostics but no topology-resolved optional `/MHD/FaceFlux` path; `mhd_io.output_divergence_diagnostics` is parsed/normalized but not consumed by snapshot code. | Record open-boundary gas/Poynting work for accepted conservation budgets without unnecessary global storage; implement requested topology output or reject it explicitly; resolve config option; v7 readback and provenance/marker test. REF 7.5/13. |
| M03-13 | P1 | `mhd_runtime_bridge.cpp:475–493` record geometry compares level/dims and only x origin/width directly (geometry generation may catch changes elsewhere). Snapshot `MhdSnapshotFieldSampler` at 981ff scans all records for each patch, O(P²). | Validate full 3-axis geometry and restart provider/epoch; indexed sampler and checked memory. |
| M03-14 | P1 | `hydro_amr_runtime.cpp:1123–1143` publishes a subset of counters only after success; pass-02 ledger says full rank telemetry is open. | Distinguish attempted/rejected/accepted, phase and per-class high water, MPI waits and messages, actual rank RSS/PSS and untracked delta; MP12 and PERF 13. |
| M03-15 | P2 | `CMakeLists.txt:1144–1180`: seven MHD-specific test executables registered but only hierarchy test in this block explicitly assigned mhd label. Pass-02 recommended `ctest -L mhd` can miss the other six. | Add labels to all MHD tests and v26 suite or use explicit names; inspect test discovery. Add missing science and MPI gates. |


## Specification coverage (source only)

| REF section | Status | Remaining |
|---|---|---|
| 3–6 equations, face authority, ideal solver | Implemented / partial consistency | UCT-HLL consistency and positivity/fallback/energy tests; no certification |
| 7 AMR/MPI | Partial | Static synchronized stepper exists; real MPI epoch defect; dynamic regrid/migration missing; recursive tier gated |
| 8 nonideal | Implemented in source | Constant/table, Ohmic/AD, Hall, RKL2, symmetric split; numerical verification and table memory; live chemistry correctly gated |
| 9 coupling | Partial | Gravity expansion separation implemented, gas thermal source hooks available; composition bypass, full coupled tests; effective EOS/dual energy gated |
| 10 IC | Partial | Uniform/loop/Alfven/stochastic source; composite seed consistency, imported face flux and shock preset missing |
| 11–12 interfaces/config | Partial | Production bridge and typed config exist; stage mutation declaration and ignored/missing option audit |
| 13 I/O | Partial | v26/v7 integrated; zero-rank schema defect, full topology export, v7 readback, restart integrity/metadata checks |
| 14–15 memory/failure | Partial | Batching/reservation and external retry present; early ungoverned allocations and communication cap violations, telemetry gaps |
| 16–20 validation/release | Not certified | All gates unexecuted; several test sources missing; no cross-code comparison |

## Gate inventory

G01 partial zero-field test; G02 partial discrete complex tests; G03 short hierarchy runs only (1000 steps/dynamic AMR not covered); G04 partial closed periodic test; G05 partial fallback assertions; G06 in-process two-rank stepper tests, real production MPI missing; G07 helper tests but production regrid missing; G08 v26 codec + roundtrip tests, v7 readback missing; G09 partial memory tests, long-run plateau missing. All NOT_RUN.

I01 Alfvén convergence source; I02 missing eigenmode suite; I03 Brio–Wu source; I04 Orszag–Tang used for patch equivalence, not full spectral convergence; I05 field-loop source, full Galilean/refinement campaign partial; I06 missing low-beta/blast/rotor suite; I07 no dynamo campaign; I08 static hierarchy tests only; I09 kernel expansion test, no production expansion case; I10 missing coupled cosmological wave/gravity verification. All NOT_RUN.

N01 Ohmic decay test (full spatial order campaign not established); N02 missing ambipolar C-shock; N03 Hall whistler source; N04 Hall energy checks source; N05 table tests but full variable-grain convergence/restart case incomplete; N06 combined smoke rather than required variable-coefficient expanding manufactured convergence; N07 RKL recurrence/source tests not full reference comparison; N08 cloud-collapse missing. All NOT_RUN.

MP01/02 partial source tests; MP03 application-level bounded rounds missing; MP04 no actual host measurements; MP05 no 100-interval plateau; MP06 no production regrid/migration tests; MP07 disabled allocations test but no pinned runtime/memory comparison; MP08 streamed codec present, full peak/zero-rank I/O audit missing; MP09 formula/bench source not measured scaling; MP10 GPU gated, not applicable to CPU qualification; MP11 policy requires finite positive budget for performance qualification (unmeasured); MP12 partial counters. BP01–BP10 benchmark source only, several report `not_evaluated`. No measured frontier.


## Updated closure priorities

1. Fail-closed admission for unsupported composition and the newly unsafe retry combinations; fix global epoch, run-wide restart selection, pre-admission allocations and spectral-source rollback. Add focused regression source.
2. Close the frozen numerical contracts (UCT actual face-flux reduction, AMR geometry/energy transfers, composite initialization). Implement conservative composition transport, static MPI bounded rounds and all shared budgets. Reuse current source owners.
3. Wire dynamic magnetic refine/coarsen and migration, guarded transactions and legal initialization; complete accepted-state output, boundary budgets and stage access. Preserve the separate RT v27 schema and existing neutrino sidecar semantics.
4. Complete verification and performance source, CMake discovery, profiling and status ledger. No compile/test execution. No new CR/RHD/GRMHD, GPU or live unverified chemistry coefficient campaign.

## Changed-file appendix

### Modified (29)
- `CMakeLists.txt`
- `CURRENT_STATUS.md`
- `README.md`
- `docs/architecture/decision_log.md`
- `docs/configuration.md`
- `docs/output_schema.md`
- `docs/profiling.md`
- `docs/repair_open_issues.md`
- `docs/repair_state_recap.md`
- `docs/restart_checkpointing.md`
- `docs/validation_plan.md`
- `include/cosmosim/core/config.hpp`
- `include/cosmosim/core/config_rt_bh_agn.hpp`
- `include/cosmosim/core/constants.hpp`
- `include/cosmosim/core/cosmology.hpp`
- `include/cosmosim/gravity/pm_solver.hpp`
- `include/cosmosim/radiation/radiation_runtime.hpp`
- `src/core/config.cpp`
- `src/core/config_rt_bh_agn.cpp`
- `src/core/cosmology.cpp`
- `src/gravity/pm_solver.cpp`
- `src/radiation/radiation_runtime.cpp`
- `src/workflows/analysis_runtime.cpp`
- `src/workflows/gravity_runtime.cpp`
- `src/workflows/internal/rt_bh_agn_bridge.cpp`
- `src/workflows/internal/rt_bh_agn_bridge.hpp`
- `src/workflows/reference_workflow.cpp`
- `src/workflows/source_runtime.cpp`
- `src/workflows/time_coordinator.cpp`

### New (74)
- `bench/bench_extended_physics_kernels.cpp`
- `configs/extended_physics/isolated_galaxy_hi_cube_rime.param.txt`
- `configs/extended_physics/neutrino_background_dmo.param.txt`
- `configs/extended_physics/neutrino_linear_response_dmo.param.txt`
- `configs/extended_physics/synthetic_eight_station_layout.txt`
- `docs/architecture/chui_extended_physics_audit_patch_cycle_v01.md`
- `docs/architecture/chui_extended_physics_bulk_implementation_claude_v01.md`
- `docs/architecture/chui_extended_physics_implementation_pass01_status.md`
- `docs/architecture/chui_extended_physics_master_architecture_v01.md`
- `docs/architecture/chui_extended_physics_numerical_contracts_pass01.md`
- `docs/architecture/chui_extended_physics_preflight_pass01.md`
- `include/cosmosim/core/config_extended_physics.hpp`
- `include/cosmosim/cosmic_rays/cosmic_ray_sources.hpp`
- `include/cosmosim/cosmic_rays/cosmic_ray_spectrum.hpp`
- `include/cosmosim/cosmic_rays/cosmic_ray_transport.hpp`
- `include/cosmosim/gravity/pm_spectral_hook.hpp`
- `include/cosmosim/instruments/interferometry.hpp`
- `include/cosmosim/lightcone/lightcone_plan.hpp`
- `include/cosmosim/lightcone/ray_geometry.hpp`
- `include/cosmosim/neutrinos/neutrino_linear_response.hpp`
- `include/cosmosim/neutrinos/neutrino_particles.hpp`
- `include/cosmosim/neutrinos/species_power_spectrum.hpp`
- `include/cosmosim/observables/cmb_secondaries.hpp`
- `include/cosmosim/observables/fits_writer.hpp`
- `include/cosmosim/observables/gas_science_view.hpp`
- `include/cosmosim/observables/line_transfer.hpp`
- `include/cosmosim/observables/observable_pipeline.hpp`
- `include/cosmosim/observables/radio_transfer.hpp`
- `include/cosmosim/observables/resonant_lyman_alpha.hpp`
- `include/cosmosim/observables/resonant_lyman_alpha_sampler.hpp`
- `include/cosmosim/radiation/radiation_causal_transport.hpp`
- `include/cosmosim/radiation/radiation_reference_transport.hpp`
- `include/cosmosim/radiation/radiation_stellar_sources.hpp`
- `include/cosmosim/spectral/atomic_data.hpp`
- `include/cosmosim/spectral/boltzmann_tables.hpp`
- `include/cosmosim/spectral/frequency_axis.hpp`
- `include/cosmosim/spectral/provider_status.hpp`
- `include/cosmosim/spectral/reference_fft.hpp`
- `src/core/config_extended_physics.cpp`
- `src/core/config_extended_physics_sections.cpp`
- `src/core/internal/config_extended_physics.hpp`
- `src/cosmic_rays/cosmic_ray_sources.cpp`
- `src/cosmic_rays/cosmic_ray_spectrum.cpp`
- `src/cosmic_rays/cosmic_ray_transport.cpp`
- `src/instruments/interferometry.cpp`
- `src/lightcone/lightcone_plan.cpp`
- `src/lightcone/ray_geometry.cpp`
- `src/neutrinos/internal/gauss_legendre.hpp`
- `src/neutrinos/neutrino_linear_response.cpp`
- `src/neutrinos/neutrino_particles.cpp`
- `src/neutrinos/species_power_spectrum.cpp`
- `src/observables/cmb_secondaries.cpp`
- `src/observables/fits_writer.cpp`
- `src/observables/gas_science_view.cpp`
- `src/observables/line_transfer.cpp`
- `src/observables/observable_pipeline.cpp`
- `src/observables/radio_transfer.cpp`
- `src/observables/resonant_lyman_alpha.cpp`
- `src/radiation/radiation_causal_transport.cpp`
- `src/radiation/radiation_reference_transport.cpp`
- `src/radiation/radiation_stellar_sources.cpp`
- `src/spectral/atomic_data.cpp`
- `src/spectral/boltzmann_tables.cpp`
- `src/spectral/frequency_axis.cpp`
- `src/spectral/reference_fft.cpp`
- `tests/unit/test_cosmic_ray_transport.cpp`
- `tests/unit/test_cosmological_background.cpp`
- `tests/unit/test_extended_physics_config.cpp`
- `tests/unit/test_instrument_rime.cpp`
- `tests/unit/test_lightcone_geometry.cpp`
- `tests/unit/test_neutrino_methods.cpp`
- `tests/unit/test_observables_kernels.cpp`
- `tests/unit/test_rhd_reference_backends.cpp`
- `tests/unit/test_spectral_data.cpp`
