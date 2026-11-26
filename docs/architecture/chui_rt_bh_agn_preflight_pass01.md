# CHUÍ RT/BH/AGN — source preflight, pass 01

**Campaign:** CHUI-RT-BH-AGN-REF-v1 · **Pass:** 01 (source implementation, NOT qualification)
**Date:** 2026-11-26 · **Mode:** feature implementation, source only (user-reserved validation exception:
no configure, build, compile, test, benchmark or simulation was run).

## 1. Baseline

| Field | Value |
|---|---|
| Repository | `AlexandreKoz/Chui-dev` (origin `https://github.com/AlexandreKoz/Chui-dev`) |
| Base branch | `main` |
| Base SHA | `6310aba6befea79a4e13e073d07ad572696bb7c8` (remote `refs/heads/main` verified with `git ls-remote`) |
| Feature branch | `claude/funny-volta-bcqf37` (session-mandated isolated branch; the charter's suggested `feature/rt_bh_agn_bulk_pass01` name was not used because the execution environment restricts pushes to this branch) |
| Worktree at start | clean, HEAD == base SHA |
| CI triggers (`.github/workflows/ci.yml`) | `pull_request`, `push` to `main`, nightly schedule, manual dispatch. A push of the feature branch does **not** trigger CI. No PR is opened. |
| Specifications | five files copied unchanged from `chui_rt_bh_agn_design_v01.zip` (no pre-existing copies) |

Specification SHA-256 (authoritative copies in `docs/architecture/`):

```
bfc50a80…cbc39c  chui_rt_bh_agn_master_design_v01.md
c93a4623…a101d16 chui_radiative_transfer_architecture_specification_v01.md
ac5f595c…14b7d   chui_black_hole_accretion_spin_specification_v01.md
eae18b9e…2dd     chui_agn_feedback_amr_specification_v01.md
160216c6…1d8     chui_rt_bh_agn_memory_performance_specification_v01.md
```

The implementation charter is archived as `agent_prompts/chui_rt_bh_agn_bulk_implementation_claude_v01.md`.

## 2. Live source anchors (read, not inferred from names)

| Area | Anchor | Finding |
|---|---|---|
| Legacy BH | `include/cosmosim/physics/black_hole_agn.hpp`, `src/physics/black_hole_agn.cpp` | `BlackHoleAgnModel::applyAccretionFromView`: host-cell Bondi, `min(Bondi, Edd)` cap, Newtonian `dm_bh = dm_acc`, thermal deposit `eps_f eps_r mdot c^2` into the host specific internal energy, metals removed proportionally. Seeding converts host gas → BH using `ParticleIdPrecommit`. Preserved unchanged as `legacy_host_cell`. |
| Seeding guard | `src/workflows/source_runtime.cpp` ~1121–1129 (constructor of `SourceRuntimeImpl`) | `physics.bh_enable_seeding=true` throws: no authoritative halo provider. Production seeding stays fail-closed for distributed runs (see §4). |
| BH callsite | `src/workflows/source_runtime.cpp` `executeDiscreteSourceEvents` ~2617 | BH runs once per source interval at the coarse-interval end epoch (`plan.event_time_code == plan.coarse_time_end_code`), after star formation and stellar returns/feedback; AMR registers and temporal histories must be empty. This is the single legal mutation epoch the new coordinator uses. |
| Halo provider | `include/cosmosim/analysis/halo_workflow.hpp` | `FofHaloFinder::buildCatalog` (local) and `buildDistributedCatalogFromView` (packs all particles to root → root-global truth, inadmissible at production scale). No distributed in-situ provider exists. |
| Thermochemistry | `include/cosmosim/physics/thermochemistry_operator.hpp`, `thermochemistry.hpp`, `chemical_network.hpp` | Pure per-row `advanceRow`, staged row result, `SourceLedger`; photo-rates enter only via `RateEnvironment::gamma_*`/`heat_*` built from the UVB table (+Rahmati shielding). No transported-radiation input and no integrated photoionization counter existed. |
| State | `include/cosmosim/core/simulation_state.hpp` | `BlackHoleParticleSidecar` (10 lanes), `MagneticFaceSoa` (core-owned, patch-keyed, empty when disabled — precedent for mesh-owned optional state), `ModuleSidecarBlock` particle-indexed rows keyed by stable particle ID (migrate with particles, persisted in restart; precedent `StellarPopulationStore`). |
| Scheduler | `include/cosmosim/core/time_integration.hpp`, `time_scheduler.hpp` | Stages `kHydroUpdate`, `kSourceTerms`; `SourceSplitPhase` Strang halves; `HierarchicalTimeBinScheduler::appendElements` for births; no element removal, but `remapSchedulerByParticleId` exists. |
| Units | `core/units.hpp`, `physics/source_interval.hpp` | Positions comoving code length; velocities physical peculiar; `rho_phys = rho_code/a^3` in comoving frame; specific internal energy is frame invariant; `SourceInterval::physical_seconds` owns `dt_phys`. |
| Restart | `src/io/restart_checkpoint.cpp` | Active constants: nonmagnetic `v25`, magnetized `v26` (= v25 + `/state/mhd`). Module sidecar blocks are persisted generically with their own schema versions. |
| Memory | `core/memory_governor.hpp`, `governed_scratch_arena.hpp` | `MemoryGovernor::reserve(MemoryClass, bytes, owner)`, classes `kCanonicalPersistent … kExternalRuntime`; `GovernedScratchArena` PMR arena. |
| MPI | `parallel/distributed_memory.hpp` | `MpiContext` reductions, `rethrowCollectivePreparationFailure`, `allgatherBytesBounded`, `exchangeBoundedAlltoallBytes`. |
| Stellar feedback | `physics/supernova_feedback.hpp` `MechanicalFeedbackDeposition::commit` | Conservative cell update convention: `m' = m+dm`, `v' = (m v + dp)/m'`, specific `u' = (m u + dE_th)/m'`, density scales with mass, element lanes `addValue`. AGN commit follows this convention. |
| MHD | `docs/architecture/chui_mhd_implementation_pass02_status.md` | Lorentz–Heaviside normalized face flux; source owners never touch magnetic faces (explicit policy). |
| GPU | `core/cuda_runtime.hpp`, `device_buffer.hpp` | CUDA runtime probing and device buffers exist; no RT kernels. |

## 3. Conflicts and compatibility amendments

See `chui_rt_bh_agn_compatibility_amendment_v01.md` (written in this pass). Summary:

- **A1 config keys:** the spec's proposed `physics.rt_*`, `physics.bh_*`, `physics.agn_*` spellings are adopted under `[physics]`; legacy keys unchanged; `physics.bh_engine_model = legacy_host_cell` (default) keeps the old callback exactly.
- **A2 BH state:** advanced BH/disk/budget truth is a particle-indexed module sidecar `chui_bh_engine_v1` keyed by stable particle ID (existing migration/restart machinery), not new `BlackHoleParticleSidecar` lanes.
- **A3 radiation state:** core-owned `RadiationMomentSoa` (patch-keyed, empty when disabled) and restart schema `v27` (= v25/v26 layout + `/state/radiation`), written only when radiation is non-empty.
- **A4 seeding:** serial FoF provider admitted only for one rank; distributed seeding fail-closed.
- **A5 chemistry interface:** additive transported photo-rate input and integrated photoionization counters in the thermochemistry ledger; no second network.
- **A6 cross sections:** H I/He II exact hydrogenic (Osterbrock & Ferland eq. 2.4); He I Verner et al. (1996) fit labelled `kUnverifiedTranscription`, gated by `physics.thermochemistry_allow_unqualified_rates`.

## 4. Feature coverage at start (before this pass)

Every P0/P1 feature was **missing** except: legacy host-cell Bondi/Eddington/thermal (`legacy_host_cell` comparator, existing), fail-closed seeding guard (existing), particle-ID precommit registry (existing, reused), module sidecar persistence (existing, reused), MHD CT authority (existing, consumed read-only).

## 5. Dependency map (implementation order)

`B03 frame` → `B02 transactions` → `B04 config` → `BH state` → `BH03 supply` → `BH04 withdrawal` → `BH05 disk` → `BH06 ledger/budget` → `AGN01 budget` → `AGN02–04 injection` → `RT01 spectrum/state` → `RT02 closure/transport` → `RT03 IMEX` → `RT04 coupling` → `RT05 cosmology` → `RT06 AMR` → `RT07 sources` → P1 → integration → persistence → tests/docs.

Truth map for new state: see §4 of `chui_rt_bh_agn_implementation_pass01_status.md`.
