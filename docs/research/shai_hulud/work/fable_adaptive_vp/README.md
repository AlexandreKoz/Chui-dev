# SHAI-HULUD — Call 1 "The Adaptive Maker": locally h-adaptive high-order phase sheet (2D Vlasov–Poisson)

Self-contained research package (Python + numba). **Not part of CHUÍ**: nothing here is wired into
the CHUÍ build, configuration, schema or runtime, and no production source file was modified.
The original Fable breakthrough sources are vendored **verbatim** under `vendor/` (checksums in
`vendor/VENDOR_SHA256.txt`); all new code is in `src/shadapt/`.

Read first: `reports/decision.md` (verdict and every arm), then
`research/local_refinement_contract.md` (definitions, proofs, invariants, known failures, 3D needs),
`reports/CHECKPOINT_FOR_CALL_2.md` (API map and exact next experiment), `RELEASE_NOTES.md`.
Every number quoted anywhere is in `results/SUMMARY.md`, which is generated from `results/*.json`.

## What is in the box

| Path | Content |
|---|---|
| `src/shadapt/forest.py` | periodic quadtree of blocks: nesting, 2:1 balance (edges, corners, seam), exact tiling, deterministic node ownership, lineage, event log, serialisation |
| `src/shadapt/hcarrier.py` | **composite hierarchical-surplus spline carrier** (cubic / quintic): globally `C^{k-1}`, interpolatory, linear, local; accumulated two-scale representation; exact line solves; tensor block evaluation |
| `src/shadapt/sheet.py` | adaptive sheet state `(d, p, Σ)` on the active nodes, carrier-based transfer (split = identity, merge, exact re-base) |
| `src/shadapt/solver.py` | leaf quadrature with exact mass ledger (lattice / adaptive, midpoint / Gauss), KDK Vlasov–Poisson driver with merged kicks, bitwise checkpoint/restart, lean N-body PM control |
| `src/shadapt/controller.py` | state-only indicators (`d6` truncation, `pred` look-ahead, `ms` multi-stream proximity, strain time-to-caustic) and the controller (threshold, hysteresis, hard node budget, cooldown, floor, gates, reasons) |
| `src/shadapt/kernels.py` | B-spline tables, compiled PCS deposit / gather / kick kernels (tested equal to the baseline PM) |
| `src/shadapt/diagnostics.py` | loop actions, interface-jump probe, offline comparison with a control (never used online) |
| `tests/` | 71 regression tests (topology, carrier, deposit, coupled solver, controller) |
| `experiments/` | `exp00` baseline reproduction … `exp08` multi-stream trigger; `arm.py` (one arm per process), `bench.py` (resumable runner), `analyze.py` (tables + gate evaluation) |
| `results/` | one JSON per experiment (command, versions, source hash, wall / CPU / RSS per arm), `SUMMARY.md`, `gate_evaluation.json`, `pytest.log` |
| `research/local_refinement_contract.md` | the technical derivation and contract |
| `reports/` | `decision.md`, `CHECKPOINT_FOR_CALL_2.md` |
| `vendor/fable_breakthrough/` | verbatim baseline `src/shsheet/*.py` (+ `sheet2d.py`), imported for the frozen conventions and for equivalence tests |

## Dependencies

Python ≥ 3.10, `numpy`, `scipy`, `numba`, `pytest` (`matplotlib` is not needed). Developed with
Python 3.13.18, numpy 2.5.3, scipy 1.18.1, numba 0.68.1, pytest 9.2.0 on a 4-core Linux VM, one
thread per process. New relative to the baseline: **numba** (compiled PM kernels for *both* the sheet
and the N-body arms) and **scipy.sparse / scipy.linalg** (carrier operators, banded Cholesky).
No GPU, no network, no compiled extension is shipped.

## One-command smoke test

```bash
cd docs/research/shai_hulud/work/fable_adaptive_vp
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q tests -p no:cacheprovider      # 71 tests, ~10 s after the first JIT compile
```

## Reproducing the experiments

Large controls are regenerated, not committed (≈45–52 MB each); they go to `$SHADAPT_DATA`
(default `results/data/`, git-ignored). Every driver is resumable: results are saved after each arm
and finished arms are skipped; the control generator checkpoints every 50 steps.

```bash
cd experiments
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 exp00_baseline_repro.py /path/to/pristine/work/fable_breakthrough   # baseline first (2 min)
for ng in 64 128 256 512; do python3 exp01_reference.py $ng 2048; done       # controls (4-5 min each)
python3 exp01_reference.py 256 4096 ; python3 exp01_reference.py 512 4096    # scoring controls (~20 min each)
python3 exp02_field_convergence.py                                           # seconds
python3 exp03_gate.py 256 --repeat ; python3 exp03_gate.py 512 --repeat      # gate, one arm at a time (~10 + ~45 min)
python3 exp04_triggers.py 256        # trigger-timing study           (3 arms at a time, ~5 min)
python3 exp05_quadrature.py 256      # quadrature study               (~3 min)
python3 exp06_indicator_validation.py 256 128 ; python3 exp06_indicator_validation.py 256 256
python3 exp06_indicator_validation.py 512 512
python3 exp07_deep_local.py 512      # refinement beyond n = ng       (~8 min)
python3 exp08_multistream_trigger.py 512 ; python3 exp08_multistream_trigger.py 256
python3 analyze.py                   # results/SUMMARY.md, results/gate_evaluation.json
```

All initial data are deterministic (no random numbers enter any run). Each JSON records the command
line, package versions, thread settings, a hash of `src/shadapt/*.py` and of the vendored baseline.

## Minimal use

```python
import sys; sys.path[:0] = ["src", "vendor/fable_breakthrough/src"]
from shadapt.sheet import AdaptiveSheet
from shadapt.solver import AdaptiveVP, Quadrature
from shadapt.controller import Controller
from shsheet.vp import GrowingModeIC

sheet = AdaptiveSheet.from_ic(GrowingModeIC(), n0=128, B=8, lmax=1, k=5)     # quintic, one local level
ctrl = Controller(indicator="pred", tol=1e-6, tau=0.6, interval=10, cooldown=10)
run = AdaptiveVP(sheet, ng=256, quadrature=Quadrature("lattice", sub0=8), controller=ctrl)
for _ in range(300):
    run.step(0.01)
print(run.conserved(), sheet.nodes_per_level(), sheet.f.events[:3])
run.save_checkpoint("state.npz")                                              # bitwise restartable
```

## Supported / not supported

Supported: 2D periodic static-background Vlasov–Poisson; cubic and quintic composite carriers on a
2:1-balanced block quadtree; split, merge, re-base; uniform-lattice quadrature with exact mass
ledger; state-only triggers with budget, cooldown, hysteresis; checkpoint/restart; N-body control.

Not supported / not done: 3D; expanding background; adjoint (energy-derived) nodal force;
anisotropic or p-refinement; block-wise time steps; sheet-to-particle release; MPI; any C++;
incremental (non-rebuilding) regrid; a literature / novelty check.
