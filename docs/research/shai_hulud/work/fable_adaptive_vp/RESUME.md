# RESUME — how to continue this project without repeating anything

State at hand-over (2026-11-17): all planned experiments are complete and committed; nothing is
running; 71/71 regression tests pass (`results/pytest.log`). Verdict and numbers:
`FINAL_REPORT.md`. Mathematics and contracts: `research/local_refinement_contract.md`.

## 1. Where everything is

* Repository `AlexandreKoz/Chui-dev`, branch `claude/blissful-feynman-udaf5r`, directory
  `docs/research/shai_hulud/work/fable_adaptive_vp/` (this tree only; production CHUÍ untouched).
* Committed: all source, tests, drivers, every result JSON, `results/SUMMARY.md`, the reports.
* **Not committed** (≈ 280 MB, regenerable): the particle controls `vp_ref_ng<NG>_n<N>.npz`. They
  live in `$SHADAPT_DATA` (default `results/data/`, git-ignored). Their SHA-256, wall time and
  conserved quantities are in `results/reference_ng*_n*.json`. Regenerate only if you need to score
  *new* arms: `python3 experiments/exp01_reference.py 256 4096` (≈ 20 min), same for 512; the
  generator checkpoints every 50 steps and resumes.
* The pristine Fable breakthrough tree is needed only by `exp00_baseline_repro.py`; it is inside
  `original_sources/shai-hulud-fable-breakthrough-results(1).zip` of the continuity pack v1.1.

## 2. Environment (two minutes)

```bash
python3 -m venv --system-site-packages venv && . venv/bin/activate
pip install numpy scipy numba pytest          # used: numpy 2.5.3, scipy 1.18.1, numba 0.68.1, pytest 9.2.0
cd docs/research/shai_hulud/work/fable_adaptive_vp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=/tmp/numba_cache
python3 -m pytest -q tests -p no:cacheprovider          # 71 passed, ~10 s after the first JIT compile
```

`PYTHONDONTWRITEBYTECODE` / `-p no:cacheprovider` / `NUMBA_CACHE_DIR` keep `__pycache__` and
`.pytest_cache` out of the tree (the CHUÍ hygiene script rejects them).

## 3. Do NOT repeat (already established, with the file that proves it)

| Question | Answer | Evidence |
|:--|:--|:--|
| Does the baseline reproduce? | bit-identically | `results/exp00_baseline_repro.json` |
| Is the composite carrier smooth / interpolatory / exact in the uniform limit / identity under splits? | yes (C⁴ quintic; leak ≤ 2λ⁸·surplus) | `tests/test_carrier.py`, contract §2 |
| Is mass exact, is there interface density, is the deposit refinement-invariant? | ledger exactly 1.0; none; yes | `tests/test_deposit.py` |
| Gate B at t = 3 on PM 256² and 512²? | **not passed** | `results/exp03_gate_ng{256,512}.json`, `SUMMARY.md` §2 |
| Equal wall time before heavy folding (t ≤ 2.5)? | sheet wins 1.1–3.6× (position) | `SUMMARY.md` "Same criteria at t = 2.0 / 2.5" |
| How many nodes after crossing? | n ≈ N_p/2 per dimension, any degree | `exp07`, `exp08` |
| Which trigger timing works? | fire by t ≈ 1.5; first-crossing trigger is 30× worse | `results/exp04_triggers_ng256.json` |
| `d6` vs `pred`? | `pred` better at equal tolerance; both blind once h ≈ Δ_pm | `exp04`, `exp06` |
| Adaptive (spacing-driven) quadrature, Gauss panels? | both lose to the uniform label lattice | `results/exp05_quadrature_ng256.json` |
| Refining beyond n = ng with truncation monitors / budgets? | no gain (budget spent in the wrong blocks) | `results/exp07_deep_local_ng512.json` |
| Multi-stream proximity trigger? | reproduces uniform n = 2ng with 84 % of the nodes, 37 % of the node-steps | `results/exp08_multistream_trigger_ng*.json` |
| Wider halo or more solver sweeps? | no effect on any result | `exp04` (halo2, sweeps2) |
| PM field convergence | 5.7e-4 between PM 256² and 512² at t = 3 | `results/exp02_field_convergence.json` |

## 4. API map (all in `src/shadapt/`)

```
forest.Forest(n0, B, lmax)            .refine(l,I,J,reason,step,t,info) .coarsen(...) .can_coarsen(...) .check()
                                      .leaves() .node_keys(l) .events .lineage .to_dict()/.from_dict()
hcarrier.HierarchicalCarrier(forest, k=5, halo_tiles=1, solver="adi"|"adi-tol"|"splu", sweeps=1)
                                      .update(zs)->[a_m]   (zs[m]: nodal data of level m, periodic fields)
                                      .eval_points(a, xi, deriv) .block_values_cf(a, L, d, W1, W2) .groups .leaf_le
sheet.AdaptiveSheet.from_ic(ic, n0, B, lmax, k)   state: .keys[l] .xi[l] .d[l] .p[l] .S[l]
                                      .coeffs("x"|"xp"|"xpS") .evaluate(xi, which, deriv) .refine(blocks) .coarsen(blocks)
                                      .rebase() .apply_topology_change(fn)  (.auto_rebase)
solver.Quadrature(mode="lattice"|"adaptive", sub0, target, nq_min, nq_max, rule="midpoint"|"gauss2")
solver.AdaptiveVP(sheet, ng, quadrature, controller)  .step(dt) .sync() .solve_field() .conserved()
                                      .sheet_integrals() .save_checkpoint(path) / AdaptiveVP.load_checkpoint(path, Controller.from_state)
solver.NBodyVP(ic, n, ng)             .step(dt) .sync() .conserved()        (lean control)
controller.Controller(indicator="pred"|"d6"|"max"|"ms", tol, tau, tau_p, interval, cooldown, max_nodes,
                      floor_cells, enable_on="always"|"precaustic"|"crossing", t_enable, buffer, coarsen, hyst,
                      auto_rebase, ms_dilate, ms_lead)   .indicators(run) .maybe_adapt(run, dt) .history
diagnostics.compare_sheet / compare_nbody (OFFLINE ONLY), loop_actions, continuity_jumps, nodal_action_defect
```

Rules that are easy to break: call `run.sync()` before reading momenta (kicks are merged);
`sub0 ≥ 2^lmax` keeps the deposit lattice globally uniform (anything else is second order at level
interfaces); the controller must never import `diagnostics.compare_*`; `floor_cells` defaults to
0.5 PM cells and **blocks useful refinement at PM 512² — pass `floor_cells=0.0` for `ms`**.

## 5. Adding an experiment

`experiments/bench.py`: `bench.sheet(ng, n0, lmax, k, quad=..., ctrl=...)` / `bench.nbody(ng, n)`
build arm specs; `bench.run_arms("name.json", arms, jobs=1)` runs each arm in its own process
(`arm.py`: timed phase without the control in memory, then offline scoring), saves after every arm
and skips finished ones. Use `jobs=1` for any number that enters an equal-wall-time statement.
Add a table to `experiments/analyze.py` and re-run it to refresh `results/SUMMARY.md`.

## 6. Known defects and open points (start here, in this order)

1. **Lead anomaly of the `ms` trigger** (PM 512²: 6.2e-9 / 7.9e-9 / 1.5e-8 for lead 0.3 / 0.6 / 1.0).
   Hypothesis: nodes created near a forming caustic inherit the coarse carrier's momentum error,
   which then grows linearly in time. Test: log `|p_new − p_ref|` at creation (offline) per split
   and compare with the final error of those blocks; one PM 512² arm (≈ 2.5 min) is enough.
2. **Pre-caustic seed for `ms`**: the mask is empty until the first fold, so bases coarser than ng
   fail (2.4e-8, 3.2e-6). In `Controller.indicators`, also mark the cells of cell centres with
   `tau_c ≤ ms_lead` (the array `tc` is already computed there). Expect first split at t ≈ 1.3.
3. **Hybrid release** (the recommended next research step): convert a block's samples to particles
   when `eta_ms` fires instead of splitting it. `AdaptiveVP.deposit` already loops over blocks; the
   released block's sample lattice (positions from the carrier, momenta from `coeffs("xp")`) is
   the particle set; its mass is the block's label mass (ledger unchanged).
4. **Adjoint force / consistent mass matrix** (Call 2): the carrier is a linear operator
   `nodal data → coefficients`; its transpose applied to `Σ_s m_s ∇V(x_s) ⊗ (evaluation weights)`
   is the energy-derived nodal force. `block_values_cf` and the sparse `P`, `Mn`, line factors give
   every piece of that transpose.
5. Prototype overheads to remove before any larger run: operators are rebuilt from scratch at each
   regrid and stored as scipy CSR (peak RSS up to 1.2 kB per node in partially refined runs);
   `Forest.n_nodes()` is recomputed per budget trial.
6. `results/exp03_gate_ng512.json`, arm `local-quin/n512L1 pred budget1.5`: wall 632 s is an
   artefact of a since-fixed controller cost bug (budget roll-back retried for every denied
   request); its accuracy (5.8e-8) is valid, its timing is not. The same configuration without the
   floor ran in 108 s in `exp07`.
7. Source hashes differ between result files (the controller gained the `ms` indicator and the
   budget fix after `exp03`); per-arm hashes are recorded from `exp07` on. Indicators `d6`/`pred`
   and all numerics of `exp03` are unaffected by those edits.

## 7. What Codex Astra needs for 3D

See `research/local_refinement_contract.md` §8: what carries over verbatim (octree, tensor
hierarchical-surplus carrier with C⁴ continuity and no hanging constraints, line solves, lattice
quadrature, `ms` trigger) and what must be redesigned (hashed tile storage, halo cost, node budget
n ≈ N_p/2 in multi-stream regions, 64-point deposit cost, adjoint force, comoving time stepping,
parallel ownership, scoring without a single converged control).
