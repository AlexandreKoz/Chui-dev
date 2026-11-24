# CHUÍ MHD — OpenCode interruption source audit

**Date:** 2026-11-24. **Scope:** focused comparative source audit, not a complete mathematical proof, build or execution. No modifications to input archives. Findings refer to the inspected post-OpenCode ZIP; verify against local Git worktree before patching.

## Evidence and scope

- Before: `Chui-dev-main(7).zip`, SHA-256 `5918ed37b7db4479204baef80400ed4bdafc1c88451873272e802d5d7599985b`.
- After: `Chui-dev-main(1)(1).zip`, SHA-256 `76850ce45a7669fcf34c3000aedffc59675f7b7f0505e522cd631fcb91ce6924`.
- Compared archive-relative regular-file paths and CRCs. No removals; 27 new files (20 MHD C++ files and seven handoff/spec/reference documents), five changed existing files. Archive entries 1162→1194 (includes directories). 6,357 physical lines in nine main .cpp and nine public .hpp files, excluding two internal SHA files. No Git metadata; commit and dirty state cannot be inferred.
- Exact normative basis: included `docs/architecture/chui_mhd_architecture_specification_v01.md`, `chui_mhd_memory_performance_specification_v01.md`, `chui_mhd_dev_branch_compatibility_amendment_v01.md`, and `agent_prompts/chui_mhd_bulk_implementation_opencode_v02.md`. Source paths/line anchors below refer to after snapshot, not future patches.
- Static inspection of key MHD kernels, core config/state and integration absence; not exhaustive inspection of all HLLD branches, SHA/table interpolation, other existing hydro/gravity modules, all performance data or every possible defect. No builds, tests, benchmarks, external simulations, or performance certification. Prior unresolved hydro/AMR defects not reassessed.

## Changed files and status

All 20 MHD C++ files are new: headers `mhd_amr`, `mhd_geometry`, `mhd_induction`, `mhd_riemann`, `mhd_runtime_adapter`, `mhd_state`, `mhd_step_operator`, `mhd_transport`, `mhd_validation` under `include/cosmosim/mhd`; matching implementations under `src/mhd`; `src/mhd/internal/mhd_sha256.cpp/.hpp`. Five existing files changed: `include/cosmosim/core/config.hpp`, `include/cosmosim/core/simulation_state.hpp`, `include/cosmosim/core/time_integration.hpp`, `src/core/config.cpp`, `src/core/simulation_state_structures.cpp`. No changes to top-level CMakeLists, workflows, production AMR, parallel, IO, tests, bench, configs. This establishes lack of new production integration, not absence of preexisting infrastructure.

| Area | Source status | Evidence |
|---|---|---|
| Config | PARTIAL | Enum/key parsing/normalization and validation in config.cpp; no runtime consumer outside MHD adapter for actual magnetic stepping; compatibility of old normalized hashes untested |
| State | PARTIAL | MagneticFaceSoa and kMagneticFaces added; lifecycle/accounting/IO not wired |
| Ideal | PARTIAL | Local SSPRK2 and HLLD/HLLE, PLM/scalar PPM, CT source, but defects below |
| AMR | PARTIAL | Standalone transfer/reflux helpers, not production regrid/level/MPI integrated |
| Nonideal | PARTIAL | Rate operators, provider classes, RKL2/Hall source exist; no actual top-level D/H invocation; constant only configured |
| Memory/PERF | PARTIAL | Some ownedCapacityBytes helpers, no complete governed allocation/parallel path; no qualification evidence |
| IO/IC/scheduler | NOT_IMPLEMENTED for MHD production | No corresponding existing owners changed |
| Validation | NOT_IMPLEMENTED for MHD test/benchmark sources | No MHD test/bench files added; diagnostic helpers are not tests |
| Ledger | NOT_IMPLEMENTED | Requested pass01 ledger absent |

## Prioritized defects and closure targets

**P0-01 — no build/production integration.** CMakeLists unchanged, no mhd source registrations; src/workflows and existing AMR/parallel/IO unchanged. Parser accepts MHD options without a demonstrated production dispatch. Add explicit unavailable-mode admission guard first, then real call path and build source registration. No claim of compiling. Ref REF §§11–13, 18; amendment R1/R2/R6.

**P0-02 — reversed FV flux signs.** `mhd_step_operator.cpp:648,707–731`: for face between f-1 and f, lower-index cell receives +F/dx and upper-index cell -F/dx although Riemann flux is in positive normal direction. Correct all lanes; apply 1/a. Verify signed two-cell flux and advection. REF §§4.3,6.1.

**P0-03 — axis/stencil addressing.** Same file:48–68,600–607,642–675: cyclic Y extents (z,x), but index/reconstruction uses (x,z), wrong on rectangular shapes. PLM two ghosts with unconditional +/-2 loads at exterior ghost center can access outside primitive storage. Fix consistent mappings and mode-specific stencil. `mhd_induction.cpp:108–207`: x-edge sweep e.g. calls faceIndexY(i,j,k) for k=nz, beyond Y-face tangential range; analogous other orientations. Correct halos/indices and shared periodic edges, not clamping incorrectly. Add noncubic and all-axis tests.

**P0-04 — incorrect UCT-HLL.** `mhd_induction.cpp:30–80` inverse-speed blend instead of signed directional HLL fan with field jumps; `mhd_riemann.cpp:128–148` provides generic nonnegative speed rather than actual signed fan bounds; `mhd_step_operator.cpp:735–770` stores x-face t1=Ey,t2=Ez, but induction edge mapping reads x t2 for Ey and t1 for Ez. Fix true UCT-HLL and one-dimensional reduction; verify orientation. Do not rely on divergence preservation alone. REF §§3,6,7.

**P0-05 — incomplete comoving energy.** `mhd_step_operator.cpp:379,477,519–549,648,1005`: a ignored in primitive recovery/budget, gas flux missing 1/a, expansion after acceptance momentum-only; B_eff not applied to Riemann. `mhd_runtime_adapter.cpp:203–289`: load omits magnetic term pending caller; store subtracts B from six-lane MhdGasStateSoa::load (`mhd_step_operator.cpp:101–110`) that doesn't populate B; uses fixed gamma and silent floor. Repair end-to-end a, EOS, energy, source ownership, stage timing, fast-CFL and state recovery; verify physical field and expansion test. REF §4.3,6.3,9. The right approach is not merely to add more B lanes as authority.

**P0-06 — AMR integrated flux mismatch.** `mhd_amr.cpp:27–48,255–352` directly copies parent integrated flux to child face or mean; for uniform B, child area is quarter of coarse area, so each of four subfaces carries parent Phi instead of Phi/4. `:158–253` restriction sums sibling block-boundary coarse faces twice (inclusive normal loops and +=). Correct true compatible prolongation, unique restriction, conservative energy adjustment and transactional topology. `MhdEmfRegister` (`:61–128`) dynamic append/hash duplicate logic does not establish full geometric contributor and temporal coverage; reflux-curl and accumulation need end-to-end review and runtime integration. REF §7.

**P0-07 — nonideal not executed.** `mhd_step_operator.cpp:795,828–1029`: coefficient provider only determines allocation; no D/H integrator invocation. Add transactional D/H/I/H/D, actual edge/gas integral export, stage exchange and source chronology. Current independently accepted dt halving is not sufficient for level-wide common time. Post-step divergence exception and source update after budget require full accepted-state rollback. REF §8.5.

**P0-08 — nonideal energy and current.** `mhd_transport.cpp:761–885` centered cell field/current then average to edges, not compatible staggered curl; total energy rate set to -J dot E at 808/870, not -div(E_ni cross B) with correct comoving factor. Ohmic/ambipolar heat must remain in closed-box total energy; Hall physical heating zero; preserve discrete energy/Poynting flux and distinguish numerical dissipation. Hall uses no whistler HLL stabilization. REF §§4.2,8.2–8.4.

**P1-09 — incorrect D/H integration and CFL; visible source compile issue.** `mhd_transport.cpp:887–934`: min-spacing squared rather than sum-of-inverse-square bound; Hall eta multiplied by field in gauss despite eta being diffusivity; audit enabled flags/units and actual fan. `:938–1090`: explicit branch at 1000 is forward Euler; RKL initialization with full dt, recurrence lacks standard L(Y0) contribution and uses inconsistent weights; at 1046 vector initialized directly from iterator (not valid vector construction). Hall `:1095–1185` copies old physical properties into stages and repeats allocations; no fresh B/T coefficients from stage state, no correct HLL fan, no conservative integral. Fix SSPRK2/RKL2, stability/accuracy, caps, conservative integrated EMFs and source-state dependencies. These are static findings, not compiler results. REF §§8.3–8.5.

**P1-10 — provider and guard.** `mhd_runtime_adapter.cpp:361–389` config constructs constant only; table and chemistry fail closed. Table provider class exists but is not production selected. Wire verified table, hash, strict boundaries and BZeroLimit, no synthetic physical data. Chemistry remains separate gated scope with unavailable validated data. Audit all config options actual use; preserve predecessor normalized config hash exact-equivalence and fail-closed runtime. REF §§8.1,12–13, amendment R5.

**P1-11 — ungoverned allocations and incomplete ownership.** `mhd_step_operator.cpp:193–231,829–835` local workspace/snapshots allocate without existing governor; ownedCapacityBytes omits cell_properties and coefficient_cache. `mhd_transport.cpp:987–1185` many temporary full-patch vectors per stage; `mhd_runtime_adapter.cpp:98–158` append canonical face vectors without governed transactional lifecycle. `simulation_state_structures.cpp:352–390` basic offset check, unchecked multiplication and incomplete ID/geometry validation; rest of core state lifecycle/accounting unchanged. No MHD OpenMP/production MPI. Repair real governance, measured capacity model, stage rollback, fixed memory caps and communication. PERF §§4–13,15–19; no measurements performed.

**P1-12 — missing lifecycle integrations.** No production state init, scheduler magnetic domains consumed, same-level MPI stage ghosts, owner EMF/gas flux synchronization, topology remap, regrid/migration, source/cooling coupling, snapshot v7, restart v26, IC parser, config profiles, test/bench registrations or ledger. Reuse existing eight-stage integration, recursive AMR and Strang source hooks; v25 is nonmagnetic dual-energy restart, v6 nonmagnetic snapshot; avoid collisions. MHD synchronized CPU/MPI reference compulsory; optional recursive MHD/GPU/chemistry fail closed if incomplete. REF §§7–18, amendment, PERF.

## Handoff and limitations

The recommendation is to continue rather than rewrite the module, prioritizing P0 corrections and runtime admission before full integration. Preserve valid topology/type work and inspect all numerical branches. No percentage completion estimate is defensible given absent production path and unexecuted verification. The audit is a targeted continuation map, not an exhaustive proof of source correctness; later source audits and build/test/science campaigns remain mandatory. No source was modified. Prompt includes detailed implementation ordering and source-only guardrails, and requires the agent to record its actual local Git identity.
