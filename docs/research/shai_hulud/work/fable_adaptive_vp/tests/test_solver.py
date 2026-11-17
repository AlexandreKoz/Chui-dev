"""Coupled Vlasov-Poisson driver: uniform-limit equivalence with the baseline,
conservation ledgers, determinism, checkpoint/restart, mid-run splits,
stress of under-resolved folds."""
import os

import numpy as np
import pytest

from shadapt import diagnostics as dg
from shadapt.controller import Controller
from shadapt.sheet import AdaptiveSheet
from shadapt.solver import AdaptiveVP, NBodyVP, Quadrature
from shsheet.pm import PM
from shsheet.sheet import Sheet
from shsheet.vp import GrowingModeIC, VPRun

IC = GrowingModeIC()
DT = 0.01


def test_nbody_matches_baseline_nodes_arm():
    b = VPRun(Sheet.from_ic(IC, 32), PM(32, 4, "grad"), "nodes")
    nb = NBodyVP(IC, 32, 32, track_action=True)
    for _ in range(40):
        b.step(DT)
        nb.step(DT)
    nb.sync()
    lean = NBodyVP(IC, 32, 32)                    # lean production path (no action, in-place kernels)
    for _ in range(40):
        lean.step(DT)
    lean.sync()
    assert np.abs(lean.x - nb.x).max() < 1e-13 and np.abs(lean.p - nb.p).max() < 1e-13
    assert np.abs(nb.x.reshape(32, 32, 2) - b.sh.x).max() < 1e-13
    assert np.abs(nb.p.reshape(32, 32, 2) - b.sh.p).max() < 1e-13
    assert np.abs(nb.S.reshape(32, 32) - b.sh.S).max() < 1e-13
    cb, cn = b.conserved(), nb.conserved()
    assert abs(cb["E"] - cn["E"]) < 1e-14 and abs(cn["mass"] - 1) < 1e-13


@pytest.mark.parametrize("car,k", (("cub", 3), ("quin", 5)))
def test_unrefined_sheet_matches_baseline_sheet_arm(car, k):
    b = VPRun(Sheet.from_ic(IC, 16), PM(32, 4, "grad"), car, 4)
    run = AdaptiveVP(AdaptiveSheet.from_ic(IC, 16, B=8, lmax=0, k=k), 32, Quadrature("lattice", 4))
    for _ in range(40):
        b.step(DT)
        run.step(DT)
    run.sync()
    sh = run.sh
    assert np.abs((sh.xi[0] + sh.d[0]).reshape(16, 16, 2) - b.sh.x).max() < 1e-13
    assert np.abs(sh.p[0].reshape(16, 16, 2) - b.sh.p).max() < 1e-13
    assert np.abs(sh.S[0].reshape(16, 16) - b.sh.S).max() < 1e-13
    assert np.abs(run.rho - b.rho).max() < 1e-11


def test_fully_refined_sheet_matches_baseline_fine_sheet():
    """Uniform-limit equivalence of the DYNAMICS: a sheet whose every block is
    split at t = 0 (two levels of nodes, composite carrier) evolves like the
    baseline uniform sheet with twice the nodes and the same deposit lattice."""
    b = VPRun(Sheet.from_ic(IC, 32), PM(32, 4, "grad"), "quin", 2)
    sh = AdaptiveSheet.from_ic(IC, 16, B=8, lmax=1, k=5)
    sh.refine(sh.f.leaves(), reason="test")
    sh.set_from_ic(IC)                                    # exact IC on all nodes (as the baseline)
    run = AdaptiveVP(sh, 32, Quadrature("lattice", 4))
    for _ in range(30):
        b.step(DT)
        run.step(DT)
    run.sync()
    xb = b.sh.x
    x0 = (sh.xi[0] + sh.d[0]).reshape(16, 16, 2)
    assert np.abs(x0 - xb[::2, ::2]).max() < 1e-12
    k1 = sh.keys[1]
    assert np.abs((sh.xi[1] + sh.d[1]) - xb[k1 // 32, k1 % 32]).max() < 1e-12
    assert np.abs(run.rho - b.rho).max() < 1e-10
    # and re-basing the fully refined sheet is exact
    xi = np.random.default_rng(0).random((2000, 2))
    before = sh.evaluate(xi, "xpS")
    assert sh.rebase() and sh.f.n0 == 32 and sh.top == 0
    after = sh.evaluate(xi, "xpS")
    assert max(np.abs(a - c).max() for a, c in zip(before, after)) < 1e-13


def _adaptive_run(nsteps=60, tol=3e-4, lmax=2, max_nodes=None, n0=16, ng=64, coarsen=False, quad=None):
    sh = AdaptiveSheet.from_ic(GrowingModeIC(eps=(0.30, 0.21, 0.09)), n0, B=8, lmax=lmax, k=5)
    ctrl = Controller(indicator="pred", tol=tol, interval=5, cooldown=10, max_nodes=max_nodes, coarsen=coarsen)
    run = AdaptiveVP(sh, ng, quad or Quadrature("lattice", 8), ctrl)
    for _ in range(nsteps):
        run.step(0.02)
    return run


def test_adaptive_run_refines_locally_and_conserves():
    run = _adaptive_run(nsteps=60, n0=32)
    sh, f = run.sh, run.sh.f
    f.check()
    assert f.events, "the controller never split anything"
    # genuinely local: several leaf levels coexist, the finest level covers only
    # part of the torus, and the node count is far below the uniform fine mesh
    assert len({l for (l, _, _) in f.leaves()}) >= 2
    top = f.top_level()
    assert top == 2 and 0.0 < len(f.exists[top]) / f.nbl(top) ** 2 < 0.5
    assert sh.n_nodes() < 0.6 * f.ncell(top) ** 2
    # every split is explained
    assert all(e["reason"].startswith("indicator") or e["reason"] in ("balance", "buffer") for e in f.events if e["kind"] == "split")
    assert all("value" in e and "tol" in e and e["value"] > e["tol"] for e in f.events if e["reason"].startswith("indicator"))
    c = run.conserved()
    assert run.mass_ledger == 1.0 and abs(c["grid_mass"] - 1) < 1e-12 and abs(c["mass"] - 1.0) < 1e-12
    assert abs(c["Px"]) < 1e-4 and abs(c["Py"]) < 1e-4          # "grad" law: truncation-level momentum drift
    # interfaces stay C^4 during the evolution
    j = dg.continuity_jumps(sh, "xp")
    assert all(j["jump"][m] <= 1e-6 * j["scale"][m] + 1e-10 for m in range(5))


def test_determinism():
    a, b = _adaptive_run(60), _adaptive_run(60)
    assert a.sh.f.signature() == b.sh.f.signature() and a.sh.f.events == b.sh.f.events
    for l in range(len(a.sh.keys)):
        assert np.array_equal(a.sh.keys[l], b.sh.keys[l])
        assert np.array_equal(a.sh.d[l], b.sh.d[l]) and np.array_equal(a.sh.p[l], b.sh.p[l]) and np.array_equal(a.sh.S[l], b.sh.S[l])
    assert np.array_equal(a.rho, b.rho)
    assert [h["n_split"] for h in a.ctrl.history] == [h["n_split"] for h in b.ctrl.history]


@pytest.mark.parametrize("mode", ("lattice", "adaptive"))
def test_checkpoint_restart_is_bitwise(tmp_path, mode):
    def make():
        sh = AdaptiveSheet.from_ic(GrowingModeIC(eps=(0.30, 0.21, 0.09)), 16, B=8, lmax=2, k=5)
        return AdaptiveVP(sh, 64, Quadrature(mode, 8, target=0.5), Controller(indicator="pred", tol=3e-4, interval=5, cooldown=10))
    a = make()
    for _ in range(70):
        a.step(0.02)
    b = make()
    for _ in range(45):
        b.step(0.02)
    path = os.path.join(tmp_path, "ck.npz")
    b.save_checkpoint(path)
    c = AdaptiveVP.load_checkpoint(path, controller_factory=Controller.from_state)
    for _ in range(25):
        c.step(0.02)
    assert c.step_count == a.step_count and c.sh.f.signature() == a.sh.f.signature()
    assert len(c.sh.f.events) == len(a.sh.f.events) > 0
    for l in range(len(a.sh.keys)):
        assert np.array_equal(a.sh.d[l], c.sh.d[l]) and np.array_equal(a.sh.p[l], c.sh.p[l]) and np.array_equal(a.sh.S[l], c.sh.S[l])
    assert np.array_equal(a.rho, c.rho) and np.array_equal(a.quad.nq, c.quad.nq)


def test_mid_run_split_conserves_mass_momentum_energy_and_field():
    sh = AdaptiveSheet.from_ic(GrowingModeIC(eps=(0.30, 0.21, 0.09)), 16, B=8, lmax=2, k=5)
    run = AdaptiveVP(sh, 64, Quadrature("lattice", 8))
    for _ in range(40):
        run.step(0.02)
    c0, rho0 = run.conserved(), run.rho.copy()
    la0 = dg.loop_actions(sh, nlines=4)
    rec = sh.refine([(0, 0, 0), (0, 1, 1)], reason="test")
    run.quad.sync(sh.car)
    run.solve_field()
    c1 = run.conserved()
    la1 = dg.loop_actions(sh, nlines=4)
    assert rec["created"] > 0
    assert c1["ledger_mass"] == 1.0 == c0["ledger_mass"]
    for q in ("Px", "Py", "K", "U", "E"):
        assert abs(c1[q] - c0[q]) < 1e-10 * max(abs(c0["K"]), 1e-3), q
    assert np.abs(run.rho - rho0).max() < 1e-8 * rho0.max()
    assert abs(la1["max_deviation"] - la0["max_deviation"]) < 1e-10


def test_hard_node_budget_and_cooldown():
    run = _adaptive_run(nsteps=90, tol=1e-5, max_nodes=900)
    assert run.sh.n_nodes() <= 900
    assert all(h["nodes"] <= 900 for h in run.ctrl.history)
    assert sum(h["n_denied_budget"] for h in run.ctrl.history) > 0, "budget was never the binding constraint"
    # cooldown: no block is split within `cooldown` steps of its creation
    born = {}
    for e in run.sh.f.events:
        if e["kind"] != "split":
            continue
        blk = (e["level"], e["I"], e["J"])
        if blk in born and e["reason"].startswith("indicator"):
            assert e["step"] - born[blk] >= run.ctrl.cooldown
        for a in (0, 1):
            for b in (0, 1):
                born[(e["level"] + 1, 2 * e["I"] + a, 2 * e["J"] + b)] = e["step"]


def test_underresolved_fold_stress_stays_finite_and_conservative():
    """n0 = 16 through shell crossing with a tight level cap: badly
    under-resolved folds.  The solver must stay finite, keep the mass ledger
    exact and the topology valid; accuracy is NOT claimed here."""
    run = _adaptive_run(nsteps=150, tol=1e-3, lmax=1)
    c = run.conserved()
    assert all(np.isfinite(v) for v in c.values())
    assert c["folded_fraction"] > 0.05 and c["det_min"] < 0          # it did fold
    assert run.mass_ledger == 1.0 and abs(c["grid_mass"] - 1) < 1e-12
    run.sh.f.check()
    assert abs(c["E"] - 0.0) < 0.5 * c["K"]                            # no blow-up
