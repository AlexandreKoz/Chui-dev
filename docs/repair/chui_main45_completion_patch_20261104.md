# CHUÍ main(45) completion patch — 2026-11-04

**Status: SOURCE EDITS ONLY — NOT CONFIGURED, NOT COMPILED, NOT TESTED, NOT RUN.**

This patch is a focused follow-up to `docs/repair/chui_final_data_and_closure_20261104.md`.
It addresses the concrete source-level gaps found in the post-Claude static audit of main(45),
without opening a new physics campaign. All executable validation remains owner-reserved.

## 1. Automatic AMR regrid safety and resource closure

The production regrid path now treats the scientific-resolution and resource-control contracts
separately:

- mandatory refinement requests (cell mass, Jeans and cooling-length criteria) are scheduled in a
  deterministic first pass and no longer consume the optional `max_operations` budget;
- optional structure/predictive requests consume a second, memory-aware budget derived from the
  current `MemoryGovernor` headroom;
- candidate `stable_partition` scratch has been removed; mandatory and optional requests are
  selected by two deterministic scans;
- planning, cooling-time preparation and predictive-feedback preparation are surrounded by
  collective preparation-failure gates before any distributed collective that depends on the
  prepared buffers;
- mandatory requests blocked by `amr_max_level` or the current 2:1 hierarchy are explicit
  collective step failures instead of silently advancing under-resolved gas;
- predictive-feedback points are exchanged globally, so a star whose drifted host moves across an
  MPI decomposition boundary can refine the destination rank;
- predictive host lookup uses admitted flat patch/row indices instead of an O(Nstar*Npatch) scan.

The planner still applies one refine/derefine pass at a synchronized hierarchy barrier. A hard
request that cannot be made legal under the current maximum level or current balance is refused;
this patch does not silently relax that physical requirement.

## 2. Recursive AMR thermochemistry / ISM / dust / SGS coupling

`ProductionAmrHydroOptions` now exposes a workflow-owned recursive-level source callback. The
production hydro runtime binds the existing `SourceRuntime` (non-owning; there is still one source
authority) and invokes the callback after each level hydro sweep, including empty-rank level clocks,
before the level end-state history is published.

For Part-2 thermochemistry this means:

1. hydro advances level L over its accepted subinterval;
2. thermochemistry/effective-ISM/dust/SGS advances exactly the rows of level L over the same
   `HydroUpdateContext` interval;
3. the resulting level end state is captured;
4. finer children interpolate from that source-updated coarse endpoint.

The coarse `SourceTerms` stage detects that recursive thermochemistry already covered the current
step and does not apply the same interval again. Adaptive chemistry shadow evidence is accumulated
over recursive substeps and the global promotion-bias decision is published once at the synchronized
coarse barrier, so the controller does not change fidelity policy halfway through one accepted
coarse interval. The recursive accumulator is part of the source rollback checkpoint.

Star formation, SSP evolution/feedback/winds, metal diffusion and BH/AGN remain synchronized-barrier
operators in this patch. They retain their own timestep/event/subcycling contracts and are not
silently duplicated at every AMR level. Consequently, this patch materially closes the stiff gas
thermochemistry part of T02 but does **not** claim that every baryonic event/source operator is now a
fully level-local multirate method. That stronger interpretation still requires numerical design and
owner-run qualification rather than an untested mechanical relocation of event operators.

## 3. Cooling-length AMR with Part-2 thermochemistry

The AMR cooling-length criterion is no longer restricted to legacy in-solver cooling.
`ThermochemistryOperator::coolingTimeSeconds` evaluates the same configured Part-2 rate set used by
the actual source operator. `SourceRuntime::fillAmrCoolingTimes` publishes per-row code-time cooling
times to the hydro/AMR owner at the synchronized regrid barrier.

Supported table/network restrictions and UV-background compatibility checks remain those of the
thermochemistry owner; the AMR planner does not create a second cooling model.

## 4. CTU imported-boundary closure and remote dual energy

CTU no longer treats every imported AMR/MPI ghost transverse increment as unknowable. When an
imported read-only ghost has a validated adjacent real owner, the solver uses that owner's
transverse increment as a counted one-sided boundary proxy. Periodic ghosts retain their wrapped
real-cell behavior; physical open/reflective boundaries retain the existing normal-only fallback.

New telemetry `ctu_transverse_boundary_proxy_states` separates this closure from genuinely degraded
or inadmissible states. This is a source-level boundary closure, **not** a claim of measured full-order
CTU convergence at AMR/MPI boundaries; the owner must run the smooth-wave/shock convergence tests.

The distributed remote-patch integration regression now supplies the remote owner's specific
internal energy and asserts that a resolvable remote boundary produces zero
`dual_energy_unresolved_ghosts`.

## 5. Tests edited/written (all NOT RUN)

- `tests/unit/test_amr_production_regrid.cpp`
  - mandatory refinement bypasses a zero optional-operation budget;
  - max-level obstruction is reported explicitly.
- `tests/unit/test_amr_recursive_subcycling.cpp`
  - three-level recursive source callback cadence is 1/2/4 and invokes every level.
- `tests/unit/test_config_chui_part3.cpp`
  - Part-2 thermochemistry + cooling-length AMR configuration is accepted; disabled cooling remains
    rejected.
- `tests/unit/test_hydro_ctu.cpp`
  - imported read-only ghost uses the counted one-sided CTU closure.
- `tests/integration/test_amr_distributed_remote_patch_boundary.cpp`
  - resolvable remote auxiliary internal energy does not fall back to E-K.

Existing distributed recursive driver, distributed diffusion and association tests from the previous
passes are preserved unchanged.

## 6. Remaining qualification / intentionally unclaimed scope

The following are not honestly closable by static editing alone and remain owner-validation items:

- full numerical convergence/order of the one-sided CTU AMR/MPI closure;
- 1-rank versus 2/4-rank end-to-end recursive-workflow comparison over realistic AMR initial
  conditions;
- natural numerical-retry/rejection campaigns and bitwise rollback comparison under the compiled
  runtime;
- calibration of adaptive chemistry, predictive AMR and memory/error control;
- performance and RSS/PSS measurements;
- scientific sensitivity to the documented yield/wind/UVB approximations.

The coordinated V02 controller remains **EXPERIMENTAL/UNQUALIFIED**. This patch makes the optional AMR
operation cap responsive to node-memory headroom while preserving hard refinement requirements, but
it does not claim an empirically calibrated global accuracy/cost optimizer.

## 7. First owner validation after applying this patch

Start with a CPU debug configure/build, then the focused AMR/CTU/source tests before broad CTest and
MPI. Any compile error or failed assertion should be treated as evidence for a corrective debugging
pass, not as a reason to add more features.
