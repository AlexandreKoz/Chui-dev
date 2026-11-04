# CHUI Part 3 — requirement-by-requirement scope and closure ledger

Date: 2026-11-03/04. Parent: `f5a92555de2c199d860755a727a7f661612912aa`
(Part 2 + Part 2.5). Branch: `claude/bold-dirac-ig5u6l` (the session's
designated feature branch, created at exactly the parent SHA; `main`, the
Part 1 branch `claude/chui-foundation-impl-hlxp2t` and the Part 2 branch
`claude/dreamy-ritchie-6rm8b1` are untouched).

**Evidence status: SOURCE ONLY. No configure, build, compiler, test,
benchmark, CI or simulation was executed (explicit owner instruction). Every
test below is NOT RUN — USER-RESERVED.** The accumulated Part 1/2/2.5 sources
have never been compiled either; compilation errors are expected.

Labels: IMPLEMENTED IN SOURCE · REFERENCE ONLY · EXPERIMENTAL/UNQUALIFIED ·
BLOCKED BY DATA/TOPOLOGY · NOT RUN — USER-RESERVED. "Open" = not yet
implemented in this campaign (never a completion claim).

This ledger is updated at every checkpoint commit; the checkpoint column
names the commit subject prefix (SHAs are listed in the final handoff
`docs/repair/chui_part3_distributed_amr_multirate_20261103.md`).

## 1. Feature register (route consolidation v10, H01–X01)

| ID | Promised behaviour | Current source path(s) | Limitation / fail-closed condition | Owner | Acceptance source tests | Status | Data restrictions | Checkpoint |
|---|---|---|---|---|---|---|---|---|
| H01 | MUSCL–Hancock, CTU, PLM/PPM, HLLC+HLLE | `src/hydro/*` (Part 1); AMR: `amr_hydro_geometry.cpp` deep ghost layers | CTU transverse terms from non-periodic AMR ghosts still normal-only (counted) | hydro | Part 1 tests; `unit_amr_periodic_cross_patch` (PLM/PPM across patches) | IMPLEMENTED IN SOURCE; NOT RUN | — | P1, A |
| H02 | dual energy, comoving CFL | Part 1 | remote AMR ghosts: auxiliary energy from local source cells only (remote counted as unresolved) | hydro | Part 1 tests | IMPLEMENTED IN SOURCE (Part 1); remote aux open | — | P1 |
| C01 | fixed hybrid network | Part 2 `thermochemistry.cpp` | rates UNQUALIFIED gate (Part 2.5) retained | physics | Part 2 tests | EXPERIMENTAL/UNQUALIFIED | rate provenance; TREECOOL/WSS09 tables | P2 |
| C02 | equilibrium/analytic + detailed reference | Part 2 | as C01 | physics | Part 2 tests | EXPERIMENTAL/UNQUALIFIED / REFERENCE ONLY | as C01 | P2 |
| C03 | adaptive fidelity, hysteresis, global sync | Part 2 local tagging (`baryonic_source_stage.cpp`) | global legal-sync scheduling: open (Checkpoint F) | physics/workflows | — | open | — | F |
| C04 | sampled shadow reference | — | open (Checkpoint F) | physics | — | open | — | F |
| I01 | effective vs resolved ISM regimes | Part 2 | — | physics | Part 2 tests | IMPLEMENTED IN SOURCE (P2) | — | P2 |
| I02 | SGS turbulence reservoir | Part 2 `budgeted_v1` | no shear production; refinement behaviour open (B/F) | physics | Part 2 | EXPERIMENTAL/UNQUALIFIED | — | P2, B, F |
| S01 | SF eligibility, conservative births | Part 2 | remote SF gradients across ranks: open (C) | physics/workflows | `unit_star_formation_gradients` | IMPLEMENTED IN SOURCE (local); remote open | — | P2, C |
| S02 | SSP returns, event timeline, associations (research) | Part 2 `stellar_population.cpp` | associations refused by config (open, C) | physics | Part 2 | IMPLEMENTED IN SOURCE; associations open | production yields need `chui_ssp_yield_table_v1` | P2, C |
| F01 | mechanical SN deposition | Part 2 | cross-rank deposition refused (open, C) | physics/workflows | Part 2 | IMPLEMENTED IN SOURCE (local) | — | P2, C |
| F02 | winds / pre-SN channels, event handling | Part 2 (pre-SN momentum EXPERIMENTAL) | time-dependent winds open (C) | physics | — | open | wind-rate tables | C |
| D01 | elements, advection, diffusion, dust | Part 2 (advection, dust); Part 3 distributed composition (`composition_transport.cpp`, `parallel/composition_exchange.cpp`) and regrid | element diffusion refused (C) | physics/amr | Part 2, `unit_amr_production_regrid` | advection + regrid + distributed IMPLEMENTED IN SOURCE (B); diffusion open (C) | — | P2, B, C |
| A01 | production automatic refine/derefine, prolong/restrict, ghosts, reflux, migration, restart | `amr_production_regrid.{hpp,cpp}`, `amr_hydro_orchestrator.cpp` (composition regrid), `hydro_amr_runtime.cpp::runAutomaticAmrRegrid` | first-order (piecewise-constant) prolongation; octets split across ranks skipped; grandparent link may be unrecoverable (blocks further coarsening); regrid ops individually transactional | amr/workflows | `unit_amr_periodic_cross_patch`, `unit_amr_production_regrid` | IMPLEMENTED IN SOURCE (A, B); NOT RUN | — | A, B |
| A02 | distributed recursive subcycling + synchronized reference | existing 2-level local subcycling | ≥3 levels / distributed: open (B/D) | amr/workflows | — | open | — | B, D |
| A03 | physics-aware indicators | `RefinementEvaluator` (mass/gradient/particles/Jeans) | shock/cooling/chemistry/feedback/gravity indicators open (B/F) | amr | — | open | — | B, F |
| A04 | predictive / cost-aware refinement | — | open (F) | amr/workflows | — | open | — | F |
| G01 | TreePM fine-gas force resolution, correction if needed | `src/gravity/*` | audit + reference tests open (E) | gravity | — | open | — | E |
| T01 | one time authority, Strang, full-interval stiff integration, event chronology | `time_coordinator.cpp`, Part 2 source stage | global all-stage rollback open (D) | workflows | Part 2 rollback test | partial; open (D) | — | P2, D |
| T02 | recursive multirate particle/AMR/source | DMO block KDK (P6) | gas/source multirate open (D) | workflows | — | open | — | D |
| M01 | canonical state, phase scratch, node-wide MPI reservations | `memory_governor`, `node_memory_admission` | deep-ghost bytes admitted via `amrHydroGeometryCapacity(patch, layers)` (A) | core | Part 1 tests | IMPLEMENTED IN SOURCE (P1); Part 3 additions per checkpoint | — | P1, A, H |
| M02 | scalar/OpenMP/MPI CPU, load balance | existing | instrumentation additions open (H) | all | — | partial | — | H |
| V01 | ledgers, mode tests, restart, provenance, science snapshots | checkpoint persistence only | science snapshot export of composition/stars open (G) | io | — | open | — | G, I |
| V02 | error/resource coordination research | — | open (F) | workflows | — | open | — | F |
| X01 | MHD/RT/new AGN boundary | — | out of scope by directive; extension seams only | — | — | BOUNDARY (not implemented by design) | — | — |

## 2. Carried debt (Part 1, Part 2, Part 2.5 handoffs and audit F1–F7)

| Debt | Origin | Status | Checkpoint |
|---|---|---|---|
| Cross-patch periodic AMR refused | Part 1 §5, Part 2 §12.1 | **IMPLEMENTED IN SOURCE** — refusal replaced by wrapped topology; only an untiled periodic box is refused | A |
| Multi-layer AMR ghost shells / PPM-CTU full order at patch boundaries | Part 1 §5 | **IMPLEMENTED IN SOURCE** for PLM/PPM (deep stencil-only layers, local + remote + wrapped); CTU ghost-side transverse terms still normal-only (counted) | A |
| Distributed composition ghosts + composition reflux payloads | Part 2 §5/§12.2 | **IMPLEMENTED IN SOURCE** — donor-fraction request/response and owner-targeted remote corrections (fine-side mirrored transfers onto remote coarse cells); refusal removed | B |
| Regrid wiring of composition prolong/restrict | Part 2 §5 | **IMPLEMENTED IN SOURCE** — CSR transfer plan in the regrid candidate commit (same weights as hydro; lane totals verified before commit) | B |
| Transactional AMR sweep | Part 1 §3 | open | B/D |
| Remote auxiliary dual energy | Part 1 §5 | open | B |
| Remote SF gradients, remote feedback | Part 2 §12.3 | open | C |
| Element-resolved diffusion (Z-only refused) | Part 2 §10.3 | open | C |
| Associations (refused) | Part 2 §0.5 | open | C/F |
| Global distributed step rollback | Part 2 §0.2 | open | D |
| Global adaptive chemistry + shadow | Part 2 §12.5 | open | F |
| Science snapshot export of composition / stars | Part 2 §0.5, §10.4 | open | G |
| Benchmark matrices, four-mode validation sources | Part 2 §10.8 | open | H/I |
| F6 runtime AMR | audit | **IMPLEMENTED IN SOURCE** — `numerics.amr_regrid_mode = automatic` drives refine/derefine in ordinary production evolution | A, B |
| `massiveReturn` quadrature monotonicity | Part 2.5 §P2.5.8 | retained 401-point check; static review in C | C |

## 3. Static integration-risk ledger (parent never compiled)

* Part 1/2/2.5 sources are uncompiled; this campaign adds to them without
  rewriting existing tests to avoid anticipated failures.
* Checkpoint A public interface changes (migration notes in the final
  handoff): `AmrHydroGeometryOptions::ghost_layers`,
  `AmrHydroGhostDescriptor::ghost_layer`, `amrHydroGeometryCapacity(patch,
  ghost_layers = 1)`, `AmrHydroGhostFillPatch::periodic_box_extent_comoving`,
  `AmrHydroGhostFillDiagnostics::deep_layer_extrapolated_ghosts`,
  `AmrHydroGhostFillStatus::kFilledDeepLayerExtrapolated`,
  `ProductionAmrHydroOptions::ghost_layers`, `productionAmrHydroGhostLayers()`,
  periodic overload of `populateAmrHydroFluxRegisterFaces`,
  `parallel::AmrInterfaceTopologyOptions` and the trailing defaulted
  `topology` parameter of `planDirectedAmrPatchBoundaryCellRequests` and
  `executeBlockingDirectedAmrPatchPayloadExchange`. All additions are
  defaulted, so existing call sites keep their historical behaviour except
  where the production runtime opts in.
* Behavioural change by design: production AMR now builds
  `stencil half-width` ghost layers (PLM 2, PPM 3) by default, so PLM/PPM keep
  full order across patch interfaces (previously reduced and counted). Results
  near patch boundaries change; not bit-reproducible with Part 2 binaries.

## 4. Owner-run commands (documentation only — NOT RUN)

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure -L chui_part3
```
