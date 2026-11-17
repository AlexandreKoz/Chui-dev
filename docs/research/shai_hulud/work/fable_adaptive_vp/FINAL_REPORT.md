# OPERATION SHAI-HULUD — THE ADAPTIVE MAKER: final report (Call 1)

Date: 2026-11-17. Scope: a locally h-adaptive, high-order, self-gravitating 2D Vlasov–Poisson
phase-sheet solver, built on the Fable breakthrough package, and its comparison with N-body PM at
equal wall time and equal memory. Research code only; **no production CHUÍ file was touched**.
Every number below is in `results/SUMMARY.md` (generated from `results/*.json` by
`experiments/analyze.py`). Evidence tags: PROVED (argument in
`research/local_refinement_contract.md`), MEASURED (result file / regression test), CONJECTURE,
FAILED, OPEN. Nothing here is peer reviewed and no novelty claim is made.

## 1. Verdict

**Gate B (local sheet ≤ N-body at equal wall time in position *and* density, and ≤ ½ N-body at
equal memory, at t = 3 with ≥ 30 % folded labels, on PM 256² and 512²): NOT PASSED → NO-GO** for
the adaptive sheet as a replacement of N-body PM in the multi-stream regime.

**CONDITIONAL GO** for three narrower things that the evidence does support:

1. the sheet as the better representation **up to about 20 % folded labels** (t ≤ 2.5 here): it
   (with n between ng/2 and ng) beats N-body at equal wall time by 1.1–3.6× in position and
   1.4–3.7× in density on both PM grids, and by 15–300× at equal peak memory;
2. the **composite hierarchical-surplus carrier** as the local-refinement kernel (globally C⁴,
   no hanging-node constraints, split = identity, exact mass ledger): it does what it was designed
   to do and extends to 3D by tensor structure;
3. the **Eulerian multi-stream proximity trigger** as the refinement criterion after shell
   crossing (it reproduces the uniform fine sheet; truncation monitors do not).

Conditions: the sheet must hand over to particles (or accept N-body-like node counts) once most
mass is multi-stream; a 3D investment is justified only for regimes where the multi-stream mass
fraction is small. This is the dossier's option (B): *retain the sheet as a memory-efficient and
geometric-analysis specialised mode while N-body remains the gravity default.*

| | strongest result | evidence |
|:--|:--|:--|
| **Positive** | At t = 2.0 (5.8 % folded), PM 512²: a quintic sheet with 65 536 nodes equals 2048² N-body (4.19 M particles) in position (2.02e-9 vs 1.92e-9) with **64× fewer degrees of freedom, 1.7× less wall time (57 s vs 97 s) and 17× less peak memory (24 vs 416 MB)**; at equal wall time it is 3.6× more accurate in position and 2.4× in density. Still true at t = 2.5 (19.5 % folded) for n = ng. | `SUMMARY.md` §2.2 and "Same criteria at t = 2.0 / 2.5" |
| **Negative** | At t = 3.0 (37 % folded), PM 512²: every sheet arm at n ≤ ng is **5–7× worse than equal-wall-time N-body in position and 3.5× in density**; matching 2048² N-body needs a sheet with n = 2·ng = N_p/2 nodes per dimension (cubic and quintic alike) at 1.2–1.4× the N-body wall time, and local refinement does not change this. | `SUMMARY.md` §2.2, §6, §7 |

## 2. What was implemented (all under `src/shadapt/`, 71 regression tests)

**Topology** (`forest.py`). Periodic quadtree of B×B-cell blocks (B = 8) with nesting, 2:1 balance
across edges, corners and the periodic seam, exact dyadic tiling, unique deterministic node
ownership (a node belongs to the coarsest level whose grid contains it; hanging nodes are ordinary
active nodes), parent/child lineage and an event log with the reason of every split
(`indicator:<name>`, `balance`, `buffer`) or merge.

**Carrier** (`hcarrier.py`) — the mathematical core. With β the cardinal B-spline of degree
k ∈ {3, 5} and T_m the level-m grid points in the level-m region plus a halo of one block ring:

```
F = Σ_m S_m ,   S_0 = periodic interpolating spline of the level-0 nodes (the baseline carrier),
S_m ∈ span{β_{m,j}: j ∈ T_m}:  S_m = surplus (node value − coarser levels) at the active level-m nodes,
                                S_m = 0 at every other grid point of T_m.
```

* P1 PROVED: each level is a symmetric positive definite problem (principal submatrix of A⊗A),
  condition number ≤ 9 (cubic), ≤ 56.25 (quintic) for any refinement shape; F is linear in the data.
* P2 PROVED + MEASURED: F ∈ C^{k−1} on the whole torus — **C⁴ for quintics across every coarse–fine
  interface** (measured normal-derivative jumps 5e-15, 4e-12, 5e-10, 5e-8, 2e-6 relative for orders
  0–4; the 5th derivative jumps by O(1)). There is no trace mismatch, hence no hanging-node
  constraint or mortar is needed.
* P3 PROVED up to a measured leak: F interpolates every active node; nodes within two cells outside
  a finer solve region are off by ≤ 2·λ⁸·surplus (λ = 0.431 quintic, 0.268 cubic).
* P4 PROVED + MEASURED: if a level is complete the carrier **is** the baseline fine spline (2e-15);
  the coupled solver then reproduces the baseline `VPRun` to 1e-12–1e-13.
* P5 PROVED up to the leak + MEASURED: **a split leaves the represented sheet unchanged** (new
  nodes are read from the old carrier: zero surplus) — density, momentum, energy and loop actions
  are conserved by the transfer to 1e-10; mass exactly. No reference value is ever used.
* P6 MEASURED: order k+1 in the fine interior, in the interface strip and in the coarse region
  (4.0 cubic, 6.2 quintic).
* Evaluation: accumulated two-scale representation; a leaf block's sample lattice is two small
  dense matrix products. Solver: one sweep of exact 1D line inverses (exact on rectangles, bands,
  full torus; defect only at re-entrant corners ≥ B cells from the data).

**Mass and deposit** (`solver.py`). Each leaf block deposits its label mass exactly once through a
tensor midpoint lattice of equal positive weights; masses are dyadic, so the **global ledger is
exactly 1.0** (tested with `==`), the grid mass 1 ± 2e-13. In lattice mode all levels share one
uniform label lattice: a cold sheet at rest gives exactly uniform density on any mesh, and a split
does not change ρ.

**Dynamics**. KDK with one field solve per step, PCS, exact-gradient force at the nodes of all
levels, exact nodal action transport, merged kicks, bitwise checkpoint/restart. *The nodal force
is the gradient of the PM potential at the node, not the derivative of the discrete energy with
respect to the node (adjoint of the quadrature): energy is not an invariant of the scheme* —
measured (E−E0)/K = −3.6e-5 at t = 3 for N-body and resolved sheets alike, momentum drift
2e-7 / 2.8e-6 (self-force of the "grad" law), carrier loop-action deviation 9e-7 (n = ng) to
1e-4 (n = ng/4).

**Triggers** (`controller.py`), all state-only (a test runs them with `open` and `numpy.load`
disabled): `d6` truncation monitor κ_k|Δ^{k+1}c| (κ₅ = 1/15360, measured to 0.2 %); `pred`
look-ahead monitor on x + τp + ½τ²a (equals `d6` of the future state for ballistic motion —
tested); strain diagnostics det Dx and time-to-caustic with `precaustic` / `crossing` gates; `ms`
Eulerian multi-stream proximity. Controller: threshold, hysteresis, hard node budget (never
exceeded), cooldown, floor, 2:1 balance, recorded reason and indicator values for every split,
automatic exact re-base when a level fills.

**Fair control** (`NBodyVP`): the same compiled kernels and force law, in-place kick/drift, one
force gather per step, two state arrays — a first version with two gathers per step and action
tracking was 2.7× slower and would have biased every equal-wall-time conclusion in favour of the
sheet (65.8 s → 24.4 s for 1024² at PM 256² in the exploratory runs; 21.8 s clean).

## 3. Baseline, controls, field convergence

* Original code, this environment: 30/30 tests; the 512² control at PM 64² and four saved arms
  reproduce **bit-identically** (`exp00`). Frozen: time variable, KDK order, IC, PM force law,
  spline conventions, metrics. New dependencies: numba 0.68.1, scipy.sparse/linalg.
* Scoring controls: 4096² particles on the same PM grid. Their own error (2048² vs 4096²):
  4.9e-10 in x at PM 256², 5.2e-9 at PM 512² (so arms below ~3e-9 at PM 512² are not resolved by
  the control). Folded label fraction at t = 3: 0.372 / 0.368.
* **PM field convergence**: positions still move by 5.7e-4 (t = 3) between PM 256² and 512²
  (1.6e-3 between 128² and 256²). Every sheet-vs-N-body difference in this report is 10³–10⁵ times
  smaller: the comparison is about the efficiency of two discretisations of the *same discrete
  model*, not about physical accuracy.

## 4. Gate results at t = 3 (clean timings: one arm at a time, one thread; repeatability ±5 %)

ratio = sheet error / N-body error read off the N-body ladder at equal cost (log-log
interpolation). Required: ≤ 1 at equal wall time, ≤ 0.5 at equal memory. RSS+ = peak resident
memory above the post-warm-up baseline; state = bytes of evolving state.

**PM 256²** (N-body ladder x_rms / wall: 512²: 4.63e-7 / 5.9 s; 768²: 5.64e-8 / 12.7 s; 1024²: 1.77e-8 / 21.8 s; 1536²: 2.13e-9 / 49.1 s)

| arm | nodes | wall s | x_rms | x ratio (wall) | rho ratio (wall) | F ratio (wall) | RSS+ MB | x ratio (RSS) | x ratio (state) | wall gate | memory gate RSS / state |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|:-:|
| uniform quintic n256, 1024² lattice | 65 536 | 15.8 | 2.57e-8 | **0.73** | **1.14** | 1.60 | 17.1 | 0.018 | 0.010 | position only | PASS / PASS |
| uniform cubic n256, 1024² lattice | 65 536 | 15.1 | 3.29e-8 | 0.85 | 1.17 | 2.02 | 16.4 | 0.021 | 0.013 | position only | PASS / PASS |
| uniform early refinement 128→256 at t = 1.5 | 65 536 | 14.2 | 2.70e-8 | 0.62 | 1.16 | 1.49 | 25.8 | 0.058 | 0.011 | position only | PASS / PASS |
| uniform refinement at t = 2.0 (late) | 65 536 | 12.8 | 1.22e-6 | 21.9 | 1.48 | 1.87 | 26.2 | 2.7 | 0.48 | FAIL | FAIL / PASS |
| **local** quintic 128→256, `pred` | 65 536 | 16.5 | 2.85e-8 | **0.89** | **1.13** | 1.65 | 81.8 | 0.99 | 0.011 | position only | FAIL / PASS |
| **local** quintic 64→256, `pred` | 65 536 | 16.8 | 2.89e-8 | 0.94 | 1.13 | 1.68 | 82.2 | 1.02 | 0.011 | position only | FAIL / PASS |
| local quintic, `d6` (same tol) | 61 032 | 16.3 | 5.13e-8 | 1.56 | 1.14 | 1.64 | 77.2 | 1.59 | 0.020 | FAIL | FAIL / PASS |
| local quintic, gate = first crossing | 65 536 | 14.4 | 7.67e-7 | 17.9 | 1.24 | 1.55 | 74.5 | 22 | 0.30 | FAIL | FAIL / PASS |
| uniform quintic n512, 2048² lattice | 262 144 | 58.8 | 9.58e-9 | 7.2† | 2.8 | 5.7 | 67.2 | 0.23 | 0.13 | FAIL | PASS / PASS |

**PM 512²** (N-body ladder: 1024²: 1.36e-7 / 21.7 s; 1536²: 2.10e-8 / 55.4 s; 2048²: 5.16e-9 / 96.8 s)

| arm | nodes | wall s | x_rms | x ratio (wall) | rho ratio (wall) | F ratio (wall) | RSS+ MB | x ratio (RSS) | x ratio (state) | wall gate | memory gate RSS / state |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|:-:|
| uniform quintic n256, 1024² lattice | 65 536 | 19.0 | 2.09e-7 | 0.98 | 0.68 | 1.12 | 25.3 | 0.055† | 0.015† | PASS | PASS / PASS |
| uniform quintic n512, 1024² lattice | 262 144 | 31.0 | 1.42e-7 | 2.12 | 4.04 | 3.22 | 67.3 | 0.29 | 0.15 | FAIL | PASS / PASS |
| uniform quintic n512, 2048² lattice | 262 144 | 70.6 | 5.74e-8 | **5.02** | **3.50** | 7.53 | 67.2 | 0.115 | 0.060 | FAIL | PASS / PASS |
| uniform cubic n512, 2048² lattice | 262 144 | 68.1 | 5.71e-8 | 4.57 | 3.43 | 7.29 | 66.2 | 0.109 | 0.060 | FAIL | PASS / PASS |
| **local** quintic 256→512, `pred` | 254 992 | 79.8 | 5.74e-8 | **6.84** | **3.45** | 8.20 | 303.8 | 5.2 | 0.065 | FAIL | FAIL / PASS |
| **local** quintic 128→512, `pred` | 254 992 | 81.4 | 5.74e-8 | 7.20 | 3.45 | 8.32 | 311.1 | 5.5 | 0.065 | FAIL | FAIL / PASS |
| local quintic 256→512, `d6` | 222 592 | 72.9 | 5.96e-8 | 5.65 | 3.50 | 7.78 | 252.8 | 3.4 | 0.055 | FAIL | FAIL / PASS |
| uniform quintic n1024, 2048² lattice (‡) | 1 048 576 | 118.7 | 6.20e-9 | ≈2.0† | ≈1.1 | — | 266 | — | — | FAIL | — |
| local `ms` 512→1024 (‡) | 876 680 | 137.1 | 6.18e-9 | ≈2.9† | ≈1.1 | — | 1031 | — | — | FAIL | — |

† N-body value extrapolated beyond the ladder. ‡ exp07/exp08: two arms ran concurrently, timings
indicative (uniform n512 measured 68.4 s there vs 70.6 s clean). The complete tables of all gate arms
(27 at PM 256², 24 at PM 512², plus repeats), with CPU time, state bytes, spectra, conservation and interface diagnostics, are in
`results/SUMMARY.md` §2 and reproduced in `reports/decision.md`.

**Equal-cost conclusions.**
* *Equal wall time, t = 3:* PM 256² — position 0.62–0.94 of N-body (sheet better), density
  1.13–1.17 (sheet worse), force error 1.5–2.0 (sheet worse): no advantage. PM 512² — sheets with
  n ≤ ng are 5–7× worse in position and 3.5× in density. The gate fails on both grids.
* *Equal memory, t = 3:* uniform sheets pass by a wide margin on measured peak RSS (position ratio
  0.018 at PM 256², 0.06–0.12 at PM 512²). Locally refined arms **fail on measured RSS** (ratio
  1.0 and 5) because this prototype stores the refinement operators as scipy sparse matrices and
  rebuilds them at every regrid; by bytes of evolving state they pass (0.011, 0.065).
* *Per-sample cost:* on the same deposit lattice the sheet is 1.37–1.38× faster than N-body
  (15.8 vs 21.8 s, 70.6 vs 96.8 s) because forces are gathered at nodes only.

**The same criteria earlier in the collapse** (supplementary; `SUMMARY.md`):

| time, folded labels | PM | arm | x ratio (wall) | rho ratio (wall) | x ratio (RSS) |
|:--|:--|:--|--:|--:|--:|
| t = 2.0, 5.8 % | 256² | uniform n256, 1024² lattice | 0.37 | 0.54 | 0.009 |
| | 256² | local 128→256 `pred` | 0.67 | 0.59 | 0.74 |
| | 512² | uniform n256, 2048² lattice | 0.28 | 0.42 | 0.003 |
| | 512² | local 256→512 `pred` | 0.64 | 0.70 | 0.49 |
| t = 2.5, 19.5 % | 256² | uniform n256, 1024² lattice | 0.38 | 0.32 | 0.010 |
| | 512² | uniform n512, 2048² lattice | 0.94 | 0.71 | 0.021 |
| | 512² | uniform n256, 2048² lattice | 1.51 | 0.66 | 0.013 |

## 5. Why the gate fails — findings

1. **After shell crossing the sheet needs n ≈ N_p/2 nodes per dimension, independent of the spline
   degree** (MEASURED). PM 512², 2048² lattice, t = 3: n256 1.69e-7, n512 5.74e-8 (cubic
   5.71e-8), n1024 6.20e-9, N-body 2048² 5.16e-9. PM 256², 1024² lattice: n256 2.57e-8, n512
   1.58e-8, N-body 1024² 1.77e-8. Before crossing the same sheets are far over-resolved
   (n256 = n1024 = N-body 2048² at t = 2.0 to 5 %). The multi-stream region has structure at the
   force-resolution scale that interpolation between nodes cannot supply (CONJECTURE: kicks
   received while crossing PM-regularised caustics; the PM kernel kinks alone are ruled out by the
   pre-crossing agreement).
2. **N-body converges fast**: x_rms ∝ wall^(−2.0…−2.4), because equal-mass particles on a label
   lattice are a lattice rule for the density functional. The uniform label lattice is also the
   best sheet quadrature: spacing-adaptive block quadrature is 30–100× worse at equal sample count
   (Euler–Maclaurin edge terms no longer cancel) and Gauss panels 3–6× worse (FAILED, `exp05`).
3. **The fold is not local on the gate problem.** 37 % of labels are folded, half of the blocks
   contain a fold and 78–86 % of the blocks lie in or next to multi-stream regions (the `ms` mask);
   at the tolerances that matter the refined region is 85–100 % of the sheet. Locality saves at most 16 % of the nodes at t = 3 (877k vs 1049k) and
   ~60 % of the node-steps, which is not where the time goes (samples dominate the cost).
4. **Truncation monitors stop identifying the limiting blocks once h ≈ Δ_pm** (MEASURED, `exp06`):
   log-correlation between the `d6` map and the true error map 0.84–0.92 for n = ng/2 but
   −0.03…0.33 for n = ng at PM 256² after t = 2 (0.2–0.6 at PM 512²); the worst blocks are at the edge of the
   folded region where the sheet is still smooth. Budget-limited refinement with these monitors
   spends the budget elsewhere: 1.5× and 2× the nodes leave the error at 5.7e-8 (`exp07`).
5. **Late refinement cannot repair committed error** (MEASURED, `exp04`): uniform refinement at
   t = 1.0 / 1.5 / 1.75 / 2.0 gives 2.58e-8 / 2.70e-8 / 5.12e-8 / 1.22e-6; a trigger gated on the
   first crossing is 30× worse than one firing by t = 1.5. State-only predictive triggers do fire
   in time (first split t = 1.4–1.6) and reproduce the uniform result (2.85e-8) with half the
   node-steps; `pred` beats `d6` at equal tolerance (2.85e-8 vs 5.13e-8).
6. **A hard budget fails ungracefully**: capped at 0.6·ng² nodes the error is 3.5e-6 (120× the
   uncapped run): the budget is exhausted before the late-time needs appear.

## 6. Alternative investigated: multi-stream proximity trigger (`ms`)

Refine every block that lies in, or within a ballistic lead of, the PM cells containing folded
mass (dilated by two cells). MEASURED (`exp08`, PM 512², 512→1024): **6.18e-9 with 877k nodes at
t = 3 and 117 M node-steps — the accuracy of the uniform n = 1024 sheet (6.20e-9, 1049k nodes,
315 M node-steps)**, where tolerance-driven truncation refinement with 79 % of the nodes gave
2.1e-8. At PM 256² it reproduces uniform n = 512 (1.56e-8 vs 1.58e-8) with 86 % of the nodes.
Open points: the result is not monotone in the lead (6.2e-9 / 7.9e-9 / 1.5e-8 for lead 0.3 / 0.6 /
1.0 at PM 512²; halo width and solver sweeps ruled out; CONJECTURE: nodes created near a forming
caustic freeze the coarse carrier's momentum error, which then grows linearly); the mask is empty
before the first fold, so a base mesh coarser than ng fails (2.4e-8, 3.2e-6) — it needs a
pre-caustic seed.

## 7. Limitations of this evidence

One test problem (single large collapse, static 2D background, three long modes); the gate's
≥ 30 % fold requirement makes it a worst case for locality. Python/numba prototype, one thread:
timings rank methods, they do not predict C++; the N-body arm is lean but not tuned beyond that.
The density error at t = 3 is dominated by sampling ripple in voids and is not monotone in N
(768²: 9.4e-5, 1024²: 1.01e-4), so density ratios within ±20 % of 1 are ties. Peak RSS of locally
refined arms measures prototype operator storage (up to 1.2 kB per node), not the method. Regrid
rebuilds operators from scratch. No adjoint force, no energy conservation, no symplecticity claim,
no expanding background, no 3D, no literature check. Control accuracy at PM 512² is 5e-9 (2048² vs
4096²). The lead anomaly of §6 is unexplained.

## 8. Next steps (in order)

1. **Hybrid release (dossier Gate F), 2D first.** Keep the sheet where the multi-stream mask is
   empty and convert blocks to particles when they enter it (the `ms` mask and the per-block
   lineage are already there; mass is per label cell, so the ledger stays exact). Pass criterion:
   N-body accuracy at t = 3 with the node count of the single-stream region.
2. **Pre-caustic seeding of the `ms` trigger** (flag blocks with time-to-caustic ≤ lead before any
   fold exists) and resolution of the lead anomaly (measure the frozen transfer error directly).
3. **Call 2 (IX-GRAVITAS)**: derive the adjoint nodal force and consistent mass matrix for this
   carrier (the cardinal functions are available through the linear carrier operator) and test
   whether a smoother/compatible field changes finding 1 — the one result that could reopen the
   equal-cost question.
4. Only then 3D (requirements listed in `research/local_refinement_contract.md` §8).

Reproduction, file map and continuation instructions: `README.md`, `RESUME.md`.
