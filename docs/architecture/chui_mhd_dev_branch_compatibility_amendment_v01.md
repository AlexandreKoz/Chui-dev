# CHUÍ MHD — Development-snapshot compatibility amendment v1

**Status:** normative compatibility amendment for the `Chui-dev-main(7).zip` snapshot; applies to OpenCode first-pass implementation only when the actual local clone matches these source interfaces.  
**Purpose:** repair dated source-state assumptions in CHUÍ-MHD-REF-v1 and CHUÍ-MHD-PERF-v1 **without changing their scientific/numerical target equations**.

## Authority and precedence

1. The real **local Git clone**, repository `AGENTS.md`, current accepted ADRs, the live source and its existing schemas take precedence for existing runtime ownership, public API and compatibility.
2. This amendment supersedes only outdated repo-state statements and explicit conflicting schema-number instructions in CHUÍ-MHD-REF-v1 and CHUÍ-MHD-PERF-v1 and the old OpenCode prompt.
3. CHUÍ-MHD-REF-v1 remains authoritative for magnetic equations, face-integrated flux authority, SSPRK2/UCT, HLLD/HLLE, conservative gas/magnetic energy coupling, Ohmic/ambipolar/Hall physics, AMR mathematics and scientific gates.
4. CHUÍ-MHD-PERF-v1 remains authoritative for MHD resource budgets, allocation lifetimes, performance gates and scientific/compute measurements.
5. If the actual local clone differs materially from this snapshot, record a new versioned reconciliation before changing APIs or schemas.

## Mandatory corrections

### R1. Checkpoint schema collision — replace all MHD v25 references with v26

The inspected repository **already writes/reads `cosmosim_restart_v25`**, whose own AMR temporal-history group schema 3 persists dual-energy auxiliary endpoints. Therefore:

- Preserve the complete existing v25 contract, validator, readback and legacy-reader behavior.
- The MHD-bearing checkpoint variant SHALL be `cosmosim_restart_v26`, version 26. Its authoritative `/state/mhd` patch-flux group shall carry the layout/metadata specified in CHUÍ-MHD-REF-v1 §13, preserving all previously required v25 data and group-layout identities.
- v26 may be read only when its full magnetic authority and integrity contracts are present and valid; non-MHD checkpoints preserve existing v25 behavior. Do not interpret a v25 checkpoint as magnetized by inferring cell-centered fields from science output.
- Revise references to `v25 MHD checkpoint` in gates G08, MP08, tests, prompt and comments to `v26 MHD checkpoint`. Revise `v24 legacy MHD-disabled checkpoint` to `v25 legacy MHD-disabled checkpoint`, while preserving older formats according to existing readers.
- Add explicit regressions: old v25 non-MHD write/read and source-on MHD v26 write/read; deliberate v26 missing/corrupt magnetic dataset, topology mismatch, dual-energy auxiliary preservation, normalized config mismatch and hash mismatch.
- **Do not reuse the v25 name, version, or existing v25 group schema identities for new magnetic fields.**

### R2. Pre-existing distributed recursive AMR and source coupling

`advanceProductionAmrHydroRecursiveImpl` and its distributed hooks are already present in `src/amr/amr_hydro_orchestrator.cpp`; `RecursiveAmrLevelSourceOperator`, Strang half-step hooks and source-interval bookkeeping are present in the development tree. They are **source implementations, not evidence of validated production capability**.

- Reuse these owners where mathematically compatible. Do not write a duplicate recursive level driver, temporal scheduler, MPI failure protocol, source-splitting authority or shadow state.
- MHD enabled with synchronized distributed AMR is the compulsory reference tier. Existing unmagnetized recursive AMR remains unchanged and available.
- A new MHD-distributed recursive mode may be enabled only after full magnetic histories, edge EMF ownership, reflux-curl, source chronology, MPI failure consensus and deterministic restart/migration rules exist and have been source-integrated. In this source-only pass, mark it `IMPLEMENTED_UNVALIDATED` or fail closed for MHD-specific mode if incomplete; do not claim scientific certification.

### R3. Existing hydrodynamic numerical features

`hydro` already includes source PPM, optional CTU, HLLE/HLLC, and dual-energy. The MHD reference SSPRK2/HLLD/UCT path is new.

- Preserve unmagnetized defaults and numerical-selection semantics. Reuse physical reconstruction and ghost-stencil infrastructure only after checking the full MHD prediction and compatibility requirements.
- Do not equate hydro `ctu` with constrained transport or with SSPRK2. They are different algorithms.
- Extend gas-magnetic total-energy authority without silently repurposing the existing dual-energy auxiliary. Document magnetic-energy exclusion when recovering gas internal energy. Maintain gas-only existing thermodynamic closure semantics where MHD disabled.
- Avoid source rewrites in `hydro_core_solver.cpp` and `amr_hydro_orchestrator.cpp` beyond carefully scoped integration seams.

### R4. Repository agent policy and naming

The local repository's `AGENTS.md`, `docs/naming_conventions.md`, `docs/file_naming_conventions.md`, `docs/architecture/runtime_truth_map.md`, `docs/architecture/developer_workflow_contract.md` and `docs/repair_open_issues.md` are mandatory orientation material. Existing namespaces `cosmosim::`, current `.cpp/.hpp` and `.param.txt` win over older external naming examples. The explicit user-reserved source-only campaign means **no build/test execution in the agent's first pass**; report this as validation deferred, not as a waiver of eventual verification.

### R5. Chemistry and coefficient-provider honesty

The included metal-ion-rate dataset is labelled unverified and the source records multiple external-data limitations. Implement genuine constant-coefficient and verified table-lookup paths as specified. Do not derive physical Ohmic/Hall/ambipolar coefficients from incomplete metal-ion fractions alone. A missing validated carrier/grain/neutral collision contract must explicitly disable or fail the requested detailed microphysics provider, never use fabricated rates.

### R6. Snapshot schema remains v7 for MHD

The local science snapshot is v6, so the v1 specification's proposed MHD v7 extension remains valid. The existing v6 output on MHD-disabled runs and current marker/transaction policy must remain unchanged.

### R7. Source-only campaign and feature status

The source snapshot includes substantial as-yet-unexecuted hydro/AMR/source enhancements. Treat them as **current source interfaces**, not proven capabilities. Record whether an issue is PRE-EXISTING, INTRODUCED, or UNKNOWN in the MHD implementation ledger. Preserve rollback and MPI coherence semantics; verify later through owner-run builds and science tests.

## Minimum current-code touchpoints

- `include/cosmosim/core/config.hpp`, `src/core/config.cpp`, `docs/configuration.md` — typed `.param.txt` configuration and normalization.
- `include/cosmosim/core/simulation_state.hpp`, `src/core/simulation_state*.cpp` — authoritative state and registry, without introducing core→mhd dependencies.
- `include/cosmosim/hydro/hydro_core_solver.hpp`, `include/cosmosim/hydro/hydro_riemann.hpp`, `include/cosmosim/hydro/hydro_reconstruction.hpp`, `docs/hydro_numerics.md` — current hydro semantics.
- `include/cosmosim/amr/amr_hydro_orchestrator.hpp`, `src/amr/amr_hydro_orchestrator.cpp`, `src/amr/amr_production_regrid.cpp`, `src/amr/amr_ghost_fill.cpp` — topology, source history, transfer, recursive stepping.
- `src/workflows/hydro_amr_runtime.cpp`, `src/workflows/source_runtime.cpp`, `src/workflows/runtime_module_registry.cpp`, `src/workflows/time_coordinator.cpp` — runtime composition and Strang/KDK stage authority.
- `src/io/restart_checkpoint.cpp`, `include/cosmosim/io/restart_checkpoint.hpp`, `src/io/snapshot_hdf5.cpp`, `docs/output_schema.md` — versioned state serialization and science output.
- `include/cosmosim/core/memory_governor.hpp`, `include/cosmosim/core/governed_scratch_arena.hpp` — physical reservation/admission.

## Admission rule

Any MHD feature that cannot be fully integrated while preserving existing invariants **fails closed only for MHD-enabled configurations**, with a concrete error and ledger entry. Do not regress baseline DMO, nonmagnetized hydro, or existing experimental modes merely to simplify MHD implementation.
