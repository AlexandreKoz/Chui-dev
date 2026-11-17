# RELEASE NOTES — work/fable_adaptive_vp (2026-11-17)

Facts only: what exists, what was run, what failed, what was not done.

## Implemented

`src/shadapt/`: `forest.py` (block quadtree), `hcarrier.py` (composite hierarchical-surplus
carrier), `sheet.py` (state and transfer), `solver.py` (quadrature, KDK driver, checkpoints, N-body
control), `controller.py` (indicators `d6`, `pred`, `ms`, strain gates; budget, cooldown,
hysteresis, floor, reasons), `kernels.py` (B-spline tables, compiled PM kernels), `diagnostics.py`.
`tests/`: 71 tests. `experiments/`: `exp00`–`exp08`, `arm.py`, `bench.py`, `analyze.py`.
Vendored verbatim: baseline `src/shsheet/*.py` and `vendor/fable_prototype/sheet2d.py`
(`vendor/VENDOR_SHA256.txt`).

## Run and passed

* Baseline (original code): `pytest -q tests` → 30 passed; `exp00`: 512² control at PM 64² and four
  saved arms bit-identical.
* `python3 -m pytest -q tests -p no:cacheprovider` → **71 passed** (`results/pytest.log`).
* Controls: 2048² particles at PM 64², 128², 256², 512²; 4096² at PM 256², 512²
  (`results/reference_*.json` hold SHA-256, timings, conserved quantities).
* `exp02` field convergence; `exp03` gate at PM 256² (33 arms incl. repeats) and PM 512² (30 arms incl.
  repeats), one arm at a time; `exp04` trigger study (25 arms); `exp05` quadrature study (14 arms); `exp06`
  indicator validation (3 cases); `exp07` deep local refinement (6 arms); `exp08` multi-stream
  trigger (8 arms at each of PM 512², 256²); `analyze.py`.
* No arm failed to run. No file of production CHUÍ was created or modified.

## Negative outcomes recorded (not bugs)

Gate B not passed at either PM grid; spacing-adaptive and Gauss-panel quadratures lose to the
uniform label lattice; first-crossing trigger 30× worse; truncation monitors blind to the limiting
blocks at h ≈ Δ_pm; budget-limited deep refinement ineffective; hard budget below the uniform
count fails ungracefully; `ms` trigger with a base coarser than ng fails; `ms` accuracy not
monotone in the lead (open).

## Incidents

* The container restarted during the PM 512² gate run. Scratch data and the committed tree
  survived; the driver resumed and skipped the three finished arms. Nothing was recomputed.
* First exploratory N-body timings were 2.7× too slow (two force gathers per step, action
  tracking, temporaries). The N-body control was made lean before any gate run; exploratory
  results were discarded (they are not in `results/`).
* Controller cost bug: with an exhausted node budget every denied request was trial-split and
  rolled back. Arm `local-quin/n512L1 pred budget1.5` of `exp03` at PM 512² therefore shows
  wall 632 s (accuracy valid). Fixed (denial without trial when the request cannot fit); the
  arm was **not** re-run; `exp07` contains the same configuration without the floor (108 s).
* The source hash stored in `exp03` (`43540d3323f00548`) predates two controller edits (budget
  fix, `ms` indicator). They do not change `d6`/`pred` decisions or any numerics of `exp03`.
  From `exp07` on each arm records the hash of the source it ran with.
* `exp04`–`exp08` ran two or three arms concurrently: their wall times are indicative only and
  are labelled as such. Peak RSS of the locally refined arms is dominated by scipy sparse
  operators that are rebuilt at every regrid.
* Timing repeatability of the gate (three repetitions): N-body 768² 12.7 / 11.9 / 11.6 s; quintic
  n256 15.8 / 14.8 / 14.8 s; local 128→256 16.5 / 17.9 / 17.3 s (PM 256²); N-body 1536²
  55.4 / 54.6 / 55.7 s; quintic n512 70.6 / 70.3 / 68.7 s; local 256→512 79.8 / 78.1 / 77.1 s (PM 512²).

## Not implemented / not attempted

Adjoint (energy-derived) force and consistent mass matrix; sheet-to-particle release; pre-caustic
seeding of the `ms` mask; anisotropic or p-refinement; incremental regrid and matrix-free
operators; block time steps; expanding background and Layzer–Irvine; 3D; MPI; C++; any CHUÍ
integration; figures; a literature / novelty check; a second test problem (the gate problem has a
large fold fraction by design); re-run of the one arm with the timing artefact.
