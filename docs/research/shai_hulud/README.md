# SHAI-HULUD — structure-preserving adaptive remapping of cold phase-space sheets

Mathematical research note for CHUI. This directory is self-contained research material
(manuscript, reference prototype, measured results). **Nothing here is wired into the CHUI
build, configuration, schema or runtime**; no production interface is changed.

## Contents

| Path | What it is |
|---|---|
| `manuscript/shai_hulud.pdf` | Compiled manuscript (31 pages): definitions, theorems, proofs, counterexamples, algorithm, error analysis, experiments, cosmological extension, comparison with prior work, ledger, verdict. |
| `manuscript/shai_hulud.tex`, `sec_*.tex` | LaTeX sources. `tables_generated.tex` is produced from the result files. |
| `proof_ledger.md` | Status of every statement: established / new derivation / numerical observation / conditional / open. |
| `prototype/sheet1d.py` | d = 1 reference implementation: exact-symplectic KDK with exact action transport, action-corrected interpolation AC(m,r), equi-affine monitor, monotone-rearrangement (1-D OT) remap, invariants. |
| `prototype/sheet2d.py` | d = 2 single-chart implementation: periodic splines, C¹/C² Hermite generating-function elements, Newton–Krylov Monge–Ampère solver on the torus, OT remap, circulation diagnostics. |
| `prototype/experiments_1d.py`, `experiments_2d.py` | Experiment drivers (write `results/*.json`). |
| `prototype/verify_propositions.py` | Symbolic/numeric checks of the closed-form statements (11 checks). |
| `prototype/make_tables.py`, `make_figures.py` | Tables and figures of the manuscript, generated from `results/`. |
| `results/` | Measured results (`*.json`) and figures. |

## The method in five lines

1. Each vertex carries its label `q`, phase-space position `(x, p)` and the **action** `Σ` (`dΣ = p·dx` along the sheet).
2. A kick–drift–kick step with a common, exact-gradient potential transports `(x, p, Σ)` exactly on a carrier manifold: mass, Lagrangian property and all Poincaré invariants are exactly conserved between remaps.
3. A monitor built from symplectic invariants (edge action defects / cubic form) defines a sampling density on **label space**.
4. An optimal-transport (Monge–Ampère) map from a fixed reference mesh to label space equidistributes it. Mass is not transported; it is recovered from the labels.
5. New vertex data are read off a reconstruction of the **generating function** (action-corrected interpolation in d = 1; C¹ Hermite elements chart-wise in d ≥ 2), whose Liouville form is exact by construction.

## Headline results (all measured; see the manuscript for tables)

- d = 1 (2-D phase space), external Hamiltonians, errors against exact samples of the discretely transported sheet:
  - Poincaré invariant of the carrier: **round-off (≤ 6e-15)** for the action-corrected remap in every test, for 1–399 remaps; interpolation remaps drift (linear up to −0.37, cubic up to 3.8e-4, quintic 6.7e-6 at 399 remaps).
  - AC(4,1) is 5th-order in the normal direction (cubic: 4th); with 39 remaps it is 70–750× more accurate than cubic remapping at equal N (N = 256–4096).
  - Collapse and chaotic tests at t = 12: fixed Lagrangian mesh has O(1) error; adaptive AC(4,1) keeps 4e-6 (N = 1024) / 3e-6 (N = 2048).
  - Smooth phase mixing ("twist"): adaptivity gives **no** benefit (stated as a negative result).
  - Total mass error: 0 in all runs.
- d = 2 (4-D phase space), single-stream regime, three Monge–Ampère remaps:
  - generating-function elements: circulation at quadrature round-off (1e-12–1e-14 for the C² element); spline carriers have genuine defects (bicubic 3e-7 → 6e-10 for n = 32 → 128).
  - exactness costs about one order of momentum accuracy relative to an unconstrained spline of equal degree.
  - Monge–Ampère Newton–Krylov: 4–5 iterations to 1e-13, convexity certificate λ_min ≥ 0.86; total mass 1 to round-off.

## What is not established

- A globally conforming, exactly Lagrangian discrete carrier across caustics in d ≥ 2 (open problem; the method is exact chart-wise there).
- Stability for unboundedly many remaps (numerical evidence only), and anything about the self-consistent Vlasov–Poisson coupling (no self-gravitating run was made).
- Priority: the literature check was limited; only the abstract of the ColDICE paper could be re-verified in this session.

## Reproducing

Requirements: Python ≥ 3.10 with `numpy`, `scipy`, `sympy`, `matplotlib`; a LaTeX distribution for the PDF.

```bash
cd docs/research/shai_hulud/prototype
python3 verify_propositions.py        # ~40 s, expects "11/11 checks passed"
python3 experiments_1d.py             # ~25 min single core (or pass job names to run in parallel)
python3 experiments_2d.py             # ~5 min
python3 make_tables.py && python3 make_figures.py
cd ../manuscript && pdflatex shai_hulud.tex && pdflatex shai_hulud.tex
```

Individual jobs: `python3 experiments_1d.py history_twist history_spiral history_driven history_comoving
convergence_twist accumulation_twist monitor_ablation` (any subset).
