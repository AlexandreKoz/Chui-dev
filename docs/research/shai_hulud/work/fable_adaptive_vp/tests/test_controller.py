"""State-only indicators and the refinement controller."""
import numpy as np

from helpers import TWO_PI
from shadapt.controller import KAPPA, Controller
from shadapt.forest import Forest
from shadapt.hcarrier import HierarchicalCarrier
from shadapt.kernels import eval_weights, prefilter_periodic
from shadapt.sheet import AdaptiveSheet
from shadapt.solver import AdaptiveVP, Quadrature
from shsheet.vp import GrowingModeIC


def test_error_constants_of_the_periodic_splines():
    """kappa_k |Delta^{k+1} c| is the leading interpolation error."""
    for k in (3, 5):
        n = 64
        x = np.arange(n) / n
        f = (np.sin(TWO_PI * x)[:, None] * np.ones((1, n)))[:, :, None]
        c = prefilter_periodic(f, k)[:, 0, 0]
        t = np.linspace(0, 1, 41)[:-1]
        js, w = eval_weights(t, k)
        err = 0.0
        for i in range(n):
            v = sum(w[0][a] * c[(i + js[a]) % n] for a in range(k + 1))
            err = max(err, np.abs(v - np.sin(TWO_PI * (i + t) / n)).max())
        dk = c.copy()
        for _ in range(k + 1):
            dk = np.roll(dk, -1) - dk
        assert abs(err / np.abs(dk).max() / KAPPA[k] - 1.0) < 0.02


class _Plane:
    """Ballistic plane-wave sheet x = q + (a0 + t) A sin(2 pi q1) e1 (no force)."""
    periods = np.zeros(2)

    def __init__(self, t, A=0.12, mode=2):
        self.t, self.A, self.m = t, A, mode

    def __call__(self, q):
        s = np.sin(TWO_PI * self.m * q[..., 0])
        p = np.stack([self.A * s, 0.02 * np.cos(TWO_PI * q[..., 1])], axis=-1)
        return q + self.t * p, p, np.zeros(q.shape[:-1])


def _ind(ic, tau, n0=32, ng=64):
    sh = AdaptiveSheet.from_ic(ic, n0, B=8, lmax=2, k=5)
    run = AdaptiveVP(sh, ng, Quadrature("lattice", 4))
    run.phi = np.zeros((ng, ng))                 # prescribed: no gravity
    run.acc = np.zeros((sh.n_nodes(), 2))
    ctrl = Controller(indicator="pred", tau=tau, tau_p=0.0)
    return ctrl.indicators(run), ctrl, run


def test_predictive_indicator_is_the_d6_indicator_of_the_future_state():
    """For ballistic motion the look-ahead monitor at time t equals the plain
    truncation monitor at t + tau exactly: it sees the under-resolution
    `tau` earlier."""
    tau = 0.25
    now, _, _ = _ind(_Plane(0.3), tau)
    fut, _, _ = _ind(_Plane(0.3 + tau), tau)
    assert np.allclose(now["eta_pred"], fut["eta_x"], rtol=1e-9, atol=1e-13)
    assert now["eta_pred"].max() > 1.5 * now["eta_x"].max()
    # strain diagnostic: the first caustic of this map is at t_c = 1/(2 pi m A)
    tc = 1.0 / (TWO_PI * 2 * 0.12)
    assert abs(now["tau_c"].min() - (tc - 0.3)) < 0.05 * tc
    assert now["det_min"].min() > 0
    late, _, _ = _ind(_Plane(tc + 0.1), tau)
    assert late["det_min"].min() < 0 and late["tau_c"].min() == 0.0


def test_indicator_predicts_the_measured_interpolation_error():
    """Unrefined mesh: kappa |Delta^6 c| / Delta_pm agrees with the measured
    max carrier error (in PM cells) within a factor 2, block by block maxima."""
    ic = _Plane(0.35)
    ind, ctrl, run = _ind(ic, 0.0)
    sh = run.sh
    xi = np.random.default_rng(0).random((40000, 2))
    (x,) = sh.evaluate(xi, "x")
    err = np.abs(x - ic(xi)[0]).max() * run.pm.ng
    assert 0.5 < ind["eta_x"].max() / err < 2.0


def test_indicators_use_only_the_state(monkeypatch):
    """No file is opened and no reference object exists when indicators and
    decisions are computed."""
    import builtins
    sh = AdaptiveSheet.from_ic(GrowingModeIC(eps=(0.3, 0.21, 0.09)), 16, B=8, lmax=2, k=5)
    run = AdaptiveVP(sh, 64, Quadrature("lattice", 8), Controller(indicator="max", tol=1e-4, interval=1, cooldown=0))
    run.step(0.02)

    def no_open(*a, **k):
        raise AssertionError("refinement decision tried to read a file")
    monkeypatch.setattr(builtins, "open", no_open)
    monkeypatch.setattr(np, "load", no_open)
    for _ in range(30):
        run.step(0.02)
    assert run.sh.f.events


def test_gates_and_trigger_ordering():
    """Predictive trigger fires before the strain (pre-caustic) gate, which
    fires before the crossing gate; all three are functions of the state."""
    first = {}
    for gate in ("always", "precaustic", "crossing"):
        sh = AdaptiveSheet.from_ic(GrowingModeIC(eps=(0.3, 0.21, 0.09)), 16, B=8, lmax=1, k=5)
        ctrl = Controller(indicator="pred", tol=1e-4, interval=2, cooldown=4, enable_on=gate, tau_lead=0.2)
        run = AdaptiveVP(sh, 64, Quadrature("lattice", 8), ctrl)
        for _ in range(70):
            run.step(0.02)
        first[gate] = ctrl.first_split_t
        if gate == "crossing":
            lat = [h for h in ctrl.history if h["enabled"]][0]
            assert lat["det_min"] <= 0.0
    assert first["always"] < first["precaustic"] < first["crossing"]


def test_patch_count_is_stable_with_hysteresis():
    """Over-refine a smooth, slowly evolving sheet by hand, then let the
    controller run with merging enabled: it merges once and does not
    oscillate (no split follows a merge of the same block)."""
    sh = AdaptiveSheet.from_ic(GrowingModeIC(eps=(0.02, 0.014, 0.006)), 32, B=8, lmax=2, k=5)
    sh.refine([(0, 1, 1), (0, 2, 2)], reason="manual")
    sh.set_from_ic(GrowingModeIC(eps=(0.02, 0.014, 0.006)))
    ctrl = Controller(indicator="max", tol=1e-3, hyst=0.1, interval=2, cooldown=4, coarsen=True)
    run = AdaptiveVP(sh, 64, Quadrature("lattice", 4), ctrl)
    for _ in range(40):
        run.step(0.02)
    kinds = [e["kind"] for e in sh.f.events if e["reason"] != "manual" and e["kind"] in ("split", "merge")]
    assert "merge" in kinds and "split" not in kinds
    counts = [len(h["levels"]) for h in ctrl.history]
    assert counts[-1] == 1 and sh.n_nodes() == 32 * 32          # back to the uniform mesh, and it stays there
    assert sum(1 for a, b in zip(counts, counts[1:]) if b > a) == 0
    assert run.mass_ledger == 1.0


def test_resolution_floor_blocks_refinement_below_the_pm_scale():
    f = Forest(16, 8, 3)
    sh = AdaptiveSheet.from_ic(GrowingModeIC(eps=(0.3, 0.21, 0.09)), 16, B=8, lmax=3, k=5)
    ctrl = Controller(indicator="pred", tol=1e-9, interval=1, cooldown=0, floor_cells=1.6)
    run = AdaptiveVP(sh, 32, Quadrature("lattice", 8), ctrl)     # PM cell = 1/32: level-1 cells (1 PM cell x stretch <= 1.4) are below the floor
    for _ in range(10):
        run.step(0.02)
    assert sh.top <= 1
    assert sum(h["n_denied_floor"] for h in ctrl.history) > 0
    assert isinstance(HierarchicalCarrier(f, k=5).memory_floats(2)["coefficients"], float)
