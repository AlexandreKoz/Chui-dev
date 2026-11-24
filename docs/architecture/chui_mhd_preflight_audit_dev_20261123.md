# CHUÍ MHD — Pre-implementation source compatibility audit

**Date:** 2026-11-23  
**Source inspected:** user-supplied `Chui-dev-main(7).zip`  
**Archive SHA-256:** `5918ed37b7db4479204baef80400ed4bdafc1c88451873272e802d5d7599985b`  
**Method:** extracted source/headers/docs and statically inspected; **no CMake configuration, compile, unit tests, science tests or runtime benchmarks**.  
**Scope:** suitability of the proposed CHUÍ-MHD-REF-v1 / CHUÍ-MHD-PERF-v1 bulk-implementation prompt against this actual development snapshot.

## Verdict

**Proceed after reconciliation, NOT with the earlier prompt unmodified.** The package is a significantly expanded development snapshot compared with the `main` tree used in the earlier planning conversation. It already contains source-level PPM/CTU, optional dual-energy, distributed recursive AMR and level-owned source hooks, and restart schema v25. There is **no dedicated MHD solver**. The added capability helps MHD integration, but the existing subsystems are substantially more complex and the repository repeatedly marks recent changes as **source only / not compiled or run**.

The ZIP has no `.git` directory, so the exact branch, Git HEAD and dirty worktree cannot be established from the archive. OpenCode must take that baseline directly from the **user's real local Git clone**.

## Confirmed facts and their consequences

| Priority | Source evidence | Confirmed source state | MHD implementation consequence |
|---|---|---|---|
| **P0** | `docs/output_schema.md` restart section; `src/io/restart_checkpoint.cpp:58–63` | Native restart **v25 already exists** and encodes AMR temporal dual-energy auxiliary history. | The MHD specification's planned **v24→v25** migration is obsolete for this snapshot. Preserve v25; allocate **v26** for MHD authoritative face fluxes, add v26 name/validator/reader/writer/readback, and test v25 legacy MHD-off compatibility. |
| **P0** | `AGENTS.md:1–113` | Repository-wide AI agent contract, architecture dependencies and v2 naming rules. | The agent must read this before editing. User's explicit no-build campaign is an acknowledged exception to normal repository evidence floor, never grounds for claiming tests passed. |
| **P0** | `src/amr/amr_hydro_orchestrator.cpp:3573–3900`; `include/cosmosim/amr/amr_hydro_orchestrator.hpp` | A recursive AMR hierarchy driver and distributed hooks are already in source; earlier architecture documentation still describes older restrictions. | Reuse the existing driver and interfaces. Initially certify synchronized MHD-MPI; implement magnetic recursive support only with complete face/edge histories and global ownership. Do **not** disable or regress existing non-MHD recursive mode. Unsupported MHD-specific combinations must fail closed. |
| **P0** | `src/workflows/source_runtime.cpp:1903–2207`; `src/workflows/hydro_amr_runtime.cpp:733–760` | Source half-step and recursive level-source hooks already exist. | Integrate MHD with existing Strang chronology; do not build a second source stage or silently apply gravity/cooling twice. |
| **P0** | `include/cosmosim/hydro/hydro_core_solver.hpp`; `docs/hydro_numerics.md`; `include/cosmosim/core/config.hpp:30–65` | Hydro supports PPM, CTU option, HLLE/HLLC, dual-energy option; these are not necessarily validated. | Keep legacy hydro untouched when MHD disabled; reuse reconstructor/metadata where mathematically appropriate, but SSPRK2+HLLD+UCT is a distinct MHD update and must not be incorrectly presented as preexisting hydro functionality. |
| **P1** | `docs/architecture/runtime_truth_map.md`; `docs/architecture/overview.md`; `src/workflows/hydro_amr_runtime.cpp`; `src/workflows/runtime_module_registry.cpp` | Production module stages are orchestrated by workflow ownership rather than standalone hydro calls. | Route config → typed runtime plan → `HydroAmrRuntime` → MHD trial/state/AMR/output; core must not depend on mhd/workflow. |
| **P1** | `include/cosmosim/core/simulation_state.hpp`; `docs/hydro_numerics.md` | Gas truth belongs to `SimulationState` and hydro operates on transient conserved SoAs; dual-energy auxiliary has own semantics. | Add authoritative staggered magnetic state with explicit patch identity/epochs and controlled energy recovery. Preserve gas sidecar truth and existing auxiliary semantics. |
| **P1** | `src/amr/amr_hydro_orchestrator.cpp`; `src/workflows/hydro_amr_runtime.cpp` | Several high-complexity AMR/hydro components already contain thousands of lines of source. | Minimize unbounded edits to existing orchestrators; use narrow interfaces and source-level integration tracing. |
| **P1** | `chui_data/chui_metal_ion_rates_v1_C_O_Si_Fe_unverified.txt`; `chui_data/README_C02_provenance.md`; `docs/repair_open_issues.md` | A detailed ion-rate table is expressly **unverified** and existing data/physical closure gaps are documented. | Constant-η and properly validated tabulated providers may be reference paths; never claim microphysical Hall/ambipolar coefficient certification from this dataset. Unsupported detailed providers fail closed. |
| **P1** | `docs/architecture/adr_param_txt_configuration_format.md`; `docs/configuration.md`; `docs/naming_conventions.md` | Typed `.param.txt` grammar and repository-authoritative v2 naming exist. | Add normalized, validated keys in existing config path, not YAML/JSON; follow current v2 naming. |
| **P1** | `include/cosmosim/core/memory_governor.hpp`; `include/cosmosim/core/governed_scratch_arena.hpp` | Governed reservations and null-upstream PMR arena already exist. | Reuse; no duplicate MHD governor, no hidden fallback heap, account co-live state and AMR transaction overlap. |
| **P1** | `include/cosmosim/io/snapshot_hdf5.hpp:42–43`; `docs/output_schema.md` | Native science snapshot is still **v6**. | Proposed **v7** MHD science snapshot does not collide. Preserve v6 on MHD-off runs. |
| **P2** | `CURRENT_STATUS.md`; `docs/repair_open_issues.md` | Multiple recent substantial code changes are described as source-only and unvalidated. | MHD-first source work may be productive, but correctness and build compatibility are unproven; the agent ledger must distinguish pre-existing defects from changes introduced in the MHD campaign. |

## Source architecture implications

The MHD subsystem should add geometrically authoritative face-flux storage and edge-EMF operators without making MHD a dependency of `core`. The ideal reference must wire through the existing workflow stage while preserving all existing hydro/thermochemical branches when `mhd.enabled=false`. The reference MHD integrator can use modern SSPRK2/UCT and HLLD independently of the existing hydro's optional CTU. Existing AMR fluid flux registers do **not** replace magnetic reflux-curl.

Versioned checkpoint migration must be **v25 (legacy non-MHD) → v26 (MHD extension)**. Do not globally change ordinary MHD-disabled runs to a new checkpoint format unless a separate migration is deliberately approved; the two schema branches must be tested for compatibility. Current v25 already has nontrivial group-specific schema validation that must remain intact.

Do not treat the recursive MPI code as certified merely because it is present. For the MHD reference, synchronized multi-rank operation is the minimum source target; distributed MHD recursive stepping can be added only with correct temporal face histories, shared-edge ownership, reflux-curl, coherent rollback and tests.

## Evidence limitations

This is a **focused source compatibility audit**, not a mathematical correctness audit of every existing hydro/AMR/chemistry implementation. The archive lacks `.git` metadata. No compilation, test execution, benchmark or science validation was performed. Source presence establishes interfaces, not functional or numerical certification.

## Handoff recommendation

Use the updated `chui_mhd_bulk_implementation_opencode_v02.md` together with `chui_mhd_dev_branch_compatibility_amendment_v01.md` and the two normative MHD specifications. Read the amendment as a **branch-specific compatibility override** of the v1 specs' obsolete repository-state facts—not as a modification of their physical equations or numerical acceptance requirements. Inspect the actual local checkout afresh before editing; if it differs from this ZIP, reconcile rather than blindly applying the amendment.
