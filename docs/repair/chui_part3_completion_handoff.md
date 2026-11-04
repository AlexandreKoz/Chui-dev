# CHUI Part 3 completion — handoff (Milestone I session 1)

**Evidence status: SOURCE ONLY. Nothing configured, compiled, linked, tested,
run or benchmarked. Every item is NOT RUN — USER-RESERVED.**

| Item | Value |
|---|---|
| Base | `d34508cc23d29d7bdb52ab9da57917019e926d27` (first Part 3 pass) |
| Branch | `claude/chui-part3-closure-t02` (new, rooted exactly at the base; `claude/bold-dirac-ig5u6l`, `main`, Part 1/2 branches untouched) |
| Milestone I subcommit 1 | `bd5ea29` recursive ≥3-level AMR subcycling in the production hydro path |
| This handoff | the commit adding this file |

**Milestone I status: PARTIAL.** The recursive gas-hierarchy clock is in the
production call graph. Distributed recursion, level-relaxed CFL,
per-level sources/events, gravity coupling per level and bounded rollback
are OPEN (§3).

## 1. Time-state ownership inventory (live source, this revision)

| Time state | Single owner | Mirrors / derived | Persisted | Legal sync point |
|---|---|---|---|---|
| Global accepted interval `[t, t+dt]`, `a(t)` | `core::IntegratorState` via `TimeCoordinator` (`runRungZeroSegment` / `runHierarchicalSegment`) | `HydroUpdateContext`, `SourceInterval` | restart | step boundary |
| Particle rungs / DMO block KDK | `HierarchicalTimeBinScheduler` (particles) + multirate PM cycle transaction | particle `time_bin` lanes (mirror) | restart | closed coarse cycle |
| Gas-cell scheduler bins | `HierarchicalTimeBinScheduler` (gas cells) | cell `time_bin` (mirror) | restart | step boundary |
| AMR level clocks (NEW) | `advanceProductionAmrHydroRecursive` local recursion (`t_i`, `dt_i = dt / r^(L-Lmin)`) | `HydroUpdateContext` per substep (`subintervalUpdate`) | no — opens and closes inside one hydro stage | start/end of the coarse interval (enforced: empty history + registers) |
| Coarse/fine temporal ghost history | `SimulationState::amr_temporal_boundary_history`, per-level records (NEW `*ForLevel`) | — | restart (but always empty at step boundaries in recursive mode) | retired per level after children |
| Coarse/fine flux registers | `SimulationState::pending_flux_registers`, interval per level pair (NEW `reflux_intervals_by_coarse_level`) | — | restart (empty at step boundary) | applied when complete |
| Gravity source generation / forces | `GravityRuntime` (TreePM), generation on `SimulationState` | cell/particle accel caches | partially | KDK stages of the global interval |
| Source operator interval | `SourceRuntime` over the global interval (Lie/Strang as Part 2) | `SourceInterval` | n/a | after hydro hierarchy synchronizes |
| Stellar event chronology | `chui_stellar_population_v1` records (cumulative-Poisson thresholds) | star particle cumulative lanes | restart | source stage |
| Step rollback point | `GlobalStepRetryPoint` (rung zero, full copy, governed) | — | no | step boundary |

## 2. What genuinely changed in the production call graph

```
TimeCoordinator (rung zero, unchanged dt authority)
 -> HydroAmrRuntime AMR stage
     amr_time_stepping == recursive_subcycling (single rank):
       advanceProductionAmrHydroRecursive            [NEW]
         validate synchronized start + 2:1 balance (periodic wrap)   [NEW]
         level i step: open level interval (StartForLevel)           [NEW]
           -> advanceProductionAmrHydro(level-i rows only,
                temporal coarse->fine ghosts from level i-1,
                registers tagged per level pair)                     [kernel reused]
           -> close interval (EndForLevel)                           [NEW]
           -> r substeps of level i+1 (recursion)                    [NEW]
           -> complete (i,i+1) reflux, RetireLevel                   [NEW]
       verify no open interval / pending register at the end         [NEW]
     synchronized (default) / distributed paths: unchanged
 -> composition commit, regrid, sources, output: unchanged order
```

## 3. Milestone I — OPEN items, in recommended order

1. **Distributed recursive sweep (P0).** `hydro_amr_runtime.cpp` refuses
   `recursive_subcycling` on more than one rank. Next edit: move the directed
   exchange block (`executeBlockingDirectedAmrPatchPayloadExchange` through
   the remote flux payload exchange, around the `if (m_mpi_context.isEnabled()
   && worldSize > 1)` branch) into a callable that the recursive driver
   invokes per level substep. Pass it via a new
   `ProductionAmrHydroOptions::per_sweep_exchange` callback, restricted to the
   level's patches and to remote patches whose time equals `t_i`. Remote
   temporal history needs the coarse owner to send start/end states:
   `kRemoteReadOnly` with `temporal_boundary_history` is already supported by
   the ghost filler. The pending-register merge must accept remote registers
   for the same per-level table.
2. **Level-relaxed CFL (P0 for the cost benefit).** `time_coordinator.cpp`
   gas CFL candidate (`computeDirectionalCflTimeStep`, around line 706): in
   recursive mode multiply each cell's CFL dt by `r^(L - Lmin)`. Use a
   globally reduced `Lmin` (allreduce min over patch levels) so rank
   decomposition cannot change dt. The HydroAmrRuntime accepted-step guard
   must apply the same rule.
3. **Per-level source/event coupling (P0).** Sources currently integrate the
   full coarse interval once, after the hierarchy synchronizes (exact interval
   coverage, Lie splitting at the coarse level). Required next:
   - an optional per-level source split with `SourceInterval` per level substep;
   - events within `[t_i, t_i+dt_i]` deposited at their own level's barrier.
   Test events crossing level boundaries.
4. **Gravity per level (P1).** Hydro gravity source terms use the global
   step's frozen acceleration in every substep. Required next:
   - interpolate the acceleration between KDK endpoints, or document the
     first-order coupling;
   - keep TreePM as the authority;
   - keep the DMO hierarchical mode untouched.
5. **Bounded rollback (P0, T01/M01).** Replace the full-state
   `GlobalStepRetryPoint` with a journal of the lanes the step mutates:
   gas rows (all lanes plus composition), touched star rows, population
   sidecar, histories, registers, the chemistry bias (`SourceRuntime::
   updateGlobalChemistryControl` mutates `ThermochemistryOperator` and is
   currently not reverted) and pending output.
6. **Composition under subcycling.** `CompositionFluxAccumulator` uses
   step-start donor fractions for all substeps (conservative, first order in
   time). Validation of negativity at large substep counts is owner-run.
7. **4-level / cosmological / mixed-bin tests**, restart at the legal
   boundary, low-memory refusal of the recursive level-row admission.

## 4. Source tests (NOT RUN — USER-RESERVED)

`unit_amr_recursive_subcycling`:
- clock structure {1, 2, 4};
- uniform state preserved at 1e-12;
- mass/momentum conservation at 1e-12 across two nested interfaces with a
  periodic wrap;
- illegal mid-cycle start refused without mutation;
- non-2:1 hierarchy refused without mutation;
- single-level equivalence with the synchronized reference at 1e-14.

Owner commands (documentation only):

```bash
cmake --preset cpu-only-debug
```

```bash
cmake --build --preset build-cpu-debug
```

```bash
ctest --preset test-cpu-debug --output-on-failure -L chui_part3
```

## 5. Follow-on milestones II–V

Still OPEN, unchanged from `CHUI_Part3_Completion_and_Design_Closure_v1.md`
§3:
- MPI element diffusion and the confirmed `elementDiffusionWorkspaceBytes`
  under-count (the `delta` buffer of `lanes x rows` doubles is not admitted);
- distributed octet derefinement;
- association grouping before candidate discovery, and internal
  velocity-dispersion energy;
- remote dual energy and CTU boundaries;
- adaptive research completion;
- 2-/4-rank MPI fixtures.

## 6. How the next session resumes

Check out `claude/chui-part3-closure-t02` at the latest commit. Read this file,
then the Milestone I section of the completion directive. Continue with
§3 item 1 (distributed recursive sweep), then items 2 and 5, keeping
Milestone I PARTIAL until all of §3 is in the production path.
