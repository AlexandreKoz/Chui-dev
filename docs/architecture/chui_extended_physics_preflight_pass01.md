# CHUI extended physics — preflight, pass 01

**Campaign:** CHUI-EXT-PHYS-REF-v1 · **Pass:** 01 (bulk source implementation; NOT qualification)
**Date:** 2026-11-27 · **Mode:** feature implementation (AGENTS.md §2) under the owner's
source-only execution exception (charter §0).

**NO CMAKE CONFIGURATION, COMPILATION, TEST, BENCHMARK, SIMULATION OR CI WAS RUN BY THIS AGENT.**

## 1. Repository and base

| Field | Value | Evidence |
|---|---|---|
| Repository | `AlexandreKoz/Chui-dev` (origin `https://github.com/AlexandreKoz/Chui-dev`) | `git remote -v` |
| Base branch | `main` | `git fetch origin main` |
| Base SHA | `ada3564d2af3b4219b92bc154b97966ee9cd209a` (merge of PR 145, RT/BH/AGN pass 01) | `git rev-parse origin/main` |
| Session branch at start | `claude/nice-lovelace-a5o24z` == base SHA, clean worktree | `git status`, `git rev-parse` |
| Feature branch | `feature/extended_physics_bulk_pass01` created from the base SHA (name did not exist locally or remotely) | `git checkout -b ... ada3564` |
| Divergence vs. design snapshot | The master was reviewed against `Chui-dev-main(8)(1).zip`. The live base contains the RT/BH/AGN pass-01 merge (ledger `chui_rt_bh_agn_implementation_pass01_status.md`, restart v27) and MHD pass 02. No later commit exists on `main`. The snapshot findings in master §2 were re-checked against live source below. | `git log origin/main` |
| CI triggers | `.github/workflows/ci.yml`: `pull_request`, `push` to `main`, nightly `schedule`, `workflow_dispatch`. A push of a feature branch triggers no job. | file read |

## 2. Design inputs

| Input | Location | SHA-256 |
|---|---|---|
| Master architecture v0.1 (unchanged copy) | `docs/architecture/chui_extended_physics_master_architecture_v01.md` | `6d97a48469f7b0bcf89201bd64e4ff519390c62e8ec71b9a8aac5e65e30e9dc9` |
| Bulk implementation charter v0.1 (unchanged copy) | `docs/architecture/chui_extended_physics_bulk_implementation_claude_v01.md` | `ee4f57b063556fbdeb5711d5ed95c017c8152b51ed5b8e85c1c6fcae9b3c13fc` |
| Audit/repair cycle template v0.1 (unchanged copy) | `docs/architecture/chui_extended_physics_audit_patch_cycle_v01.md` | `bab79dec1188728ea39e3098042b845d3df41f30e492fe64080c255039e03bca` |

No earlier accepted copy of these files existed in the repository. The `.docx` rendering was not supplied
and is not required. The eleven detailed specifications listed in master §17 do **not** exist; this pass
therefore writes its own equation-to-code notes in `chui_extended_physics_numerical_contracts_pass01.md`.

### Governance documents inspected (live; targeted sections, not every addendum line)

`AGENTS.md`, `CURRENT_STATUS.md`, `README.md`, `CONTRIBUTING.md`, `docs/architecture/overview.md`,
`developer_workflow_contract.md`, `runtime_truth_map.md`, `adr_runtime_truth_ownership.md`,
`adr_param_txt_configuration_format.md`, `decision_log.md`, `docs/repair_state_recap.md`,
`docs/repair_open_issues.md`, `docs/naming_conventions.md`, `docs/file_naming_conventions.md`,
`docs/configuration.md`, `docs/output_schema.md`, `docs/restart_checkpointing.md`,
`docs/validation_plan.md`, `docs/memory_governance.md`, RT/BH/AGN design + ledger + amendment,
MHD pass-01/02 ledgers. Missing documents: the eleven master §17 subsystem specifications;
`chui_extended_physics_memory_performance_specification_v01.md` (the memory model of this pass is the
live `MemoryGovernor` contract plus the formulas in the ledger, estimates only).

### Precedence and conflicts recorded

1. The master's YAML is illustrative; live `AGENTS.md` and `adr_param_txt_configuration_format.md` win:
   every new option is a `.param.txt` key in the single typed `SimulationConfig` path.
2. Master directory names (`neutrinos`, `cosmic_rays`, `lightcone`, `observables`, `instruments`) are
   adopted under `include/cosmosim/` and `src/` (no existing module conflicts). The shared cosmological
   background authority is placed in `core` because `LambdaCdmBackground` (core) is already the
   single H(a) authority consumed by the integrator; extending it avoids a second authority (S01).
3. Restart/snapshot schema: live current versions are restart v27 (`/state/radiation`), science snapshot v6 with MHD v7 (`snapshot_hdf5.hpp`), not an
   assumed number. Any new evolved state is versioned from the live value.

## 3. Live implementation map (base evidence)

| Area | Live owner / source | Finding |
|---|---|---|
| H(a) authority | `include/cosmosim/core/cosmology.hpp` `LambdaCdmBackground`; constructed in `src/core/config.cpp` (integration authority), `src/workflows/time_coordinator.cpp` (cosmology dt), `src/workflows/reference_workflow.cpp` (runtime) — each copies only `hubble_param, omega_matter, omega_lambda`. `omega_radiation`/`omega_curvature` exist in the struct but no config key feeds them. | Three independent construction sites = latent duplicate-authority risk; no photons, no neutrinos, no T_CMB. S01/N01 must route all sites through one factory. |
| Cosmology config | `CosmologyConfig` (`omega_matter`, `omega_lambda`, `omega_baryon`, `hubble_param`, `sigma8`, `scalar_index_ns`, box) | No neutrino, CMB temperature, N_eff or curvature keys. |
| Particle species | `core::ParticleSpecies` {DM, gas, star, BH, tracer}; `k_particle_species_count` sizes many arrays (softening, I/O, accounting) | A neutrino particle species does not exist; adding one changes the restart/snapshot species axis (N04 integration blocker). |
| PM gravity | `src/gravity/pm_solver.cpp` `applyDerivedPoissonOperator`: single periodic Poisson operator `-4 pi G W^-2 S(k) / k^2`; options fingerprinted collectively | Single insertion point for an isotropic spectral source multiplier (N03); isolated/open PM is separate. |
| Radiation | `src/radiation/*` (M1, HLL/GLF, PD-ARS, grey coupling, edge redshift transfer, AMR operators, AGN sources, P1 rays/VET tensor/MC) | Matches master §2: admission (`src/core/config_rt_bh_agn.cpp::validateRadiation`) refuses non-M1 transport, stellar sources, MPI>1, regrid, async AMR, MHD+RT. |
| Radiation spectrum | `radiation_spectrum.hpp`: `ionizing3` uses exact thresholds 13.598434 / 24.587389 / 54.417765 eV with explicit upper edge | R02 three-band baseline already EXISTING_UNQUALIFIED. |
| Angular | `radiation_angular.hpp`: steady-state rays, `vetFromAngularSamples` (tensor only), isotropic elastic MC | Confirms master: no VET transport, no causal rays, no anisotropic scattering. |
| Stellar sources | `validateRadiation`: `rt_stellar_sources` → FAIL_CLOSED ("no age/metallicity-dependent provider") | F02 open. |
| Memory | `core::MemoryGovernor` (`reserve`, `MemoryReservation::commit/release`, classes), `GovernedScratchArena` | Reused by every new module for population/mesh/packet allocations. |
| MHD | `cosmosim_mhd`, restricted synchronized CPU leaf-AMR dispatch; dynamic magnetic regrid refused | F03: refusals retained; CR transport needs B from this authority. |
| Restart / snapshot | `src/io/restart_checkpoint.cpp` v27 (`/state/radiation`), MHD v26; science snapshot schema 6, MHD snapshot 7 | F04 open (no radiation science export). |
| Observables | `src/analysis/diagnostics.cpp` (power spectrum diagnostics, mesh n) | No ray/lightcone, Ly-a, CMB, radio or instrument modules. |
| Cosmic rays | none | Not established (master §2 confirmed). |

## 4. Slice plan (dependency order)

S0 preflight/spec → S1 shared background authority (S01/N01) → S2 tables and atomic/spectral data
(S02) → S3 neutrino methods (N02–N05) → S4 ray/lightcone service (S03) → S5 RHD extensions and the
F01/F02 foundation items that are reachable in source → S6 cosmic rays (CR01–CR07) → S7 observables
(L, C, RA) → S8 instruments (O) → S9 integration, config, CMake, docs, ledger.

## 5. Explicit pass-01 decisions

* No new restart schema version is introduced unless a new **evolved authoritative** field is added to
  canonical runtime state; observables and instruments are read-only consumers with their own versioned
  product schemas.
* Every new runtime-physics option is admitted only through `validateExtendedPhysicsConfig`; options
  whose integration edge is missing are rejected with the exact missing edge (FAIL_CLOSED), never
  silently ignored. Disabled defaults emit no normalized text, preserving config hashes and restart
  admission of existing runs.
* External data (CLASS/CAMB tables, SSP luminosity tables, collisional-coupling tables, CR loss-rate
  tables) are never bundled or invented: loaders verify schema, SHA-256 and domain and refuse otherwise.
