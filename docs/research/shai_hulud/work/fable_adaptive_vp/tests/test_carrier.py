"""Composite hierarchical carrier: exactness, interpolation, uniform limit,
continuity and order at coarse-fine interfaces, locality, linear solver."""
import numpy as np
import pytest

from helpers import nodal, nodes_xi, patch_forest, random_forest, smooth
from shadapt.forest import Forest
from shadapt.hcarrier import HierarchicalCarrier
from shadapt.kernels import DECAY
from shadapt.sheet import AdaptiveSheet
from shadapt import diagnostics as dg
from shsheet.pspline import PSpline2, grid


@pytest.mark.parametrize("k", (3, 5))
def test_constants_are_reproduced_exactly(k):
    f = random_forest(32, 8, 3, nsplit=20, seed=1)
    car = HierarchicalCarrier(f, k=k)
    zs = [np.full((len(f.node_keys(l)), 2), [0.37, -1.25]) for l in range(f.top_level() + 1)]
    a = car.update(zs)
    xi = np.random.default_rng(0).random((4000, 2))
    assert np.abs(car.eval_points(a, xi) - [0.37, -1.25]).max() < 2e-15


def test_identity_sheet_is_affine_exact():
    """x(xi) = xi + const, p = const: the degree-one (affine) part is carried
    outside the spline, so affine sheets are exact on any topology."""
    f = random_forest(32, 8, 3, nsplit=20, seed=2)
    sh = AdaptiveSheet(f, k=5)
    for l in range(len(sh.keys)):
        sh.d[l][:] = [0.3, -0.2]
        sh.p[l][:] = [0.1, 0.05]
    xi = np.random.default_rng(3).random((3000, 2))
    (x, gx), (p, gp) = sh.evaluate(xi, "xp", deriv=1)
    assert np.abs(x - xi - [0.3, -0.2]).max() < 1e-15
    assert np.abs(gx - np.eye(2)).max() < 1e-12 and np.abs(gp).max() < 1e-12


@pytest.mark.parametrize("k", (3, 5))
def test_level0_is_the_baseline_spline(k):
    f = Forest(32, 8, 2)
    car = HierarchicalCarrier(f, k=k)
    a = car.update(nodal(f))
    xi = np.random.default_rng(0).random((3000, 2))
    base = PSpline2(smooth(grid(32).reshape(-1, 2)).reshape(32, 32, 2), k)
    assert np.abs(car.eval_points(a, xi) - base(xi)).max() < 5e-15


@pytest.mark.parametrize("k", (3, 5))
def test_uniform_limit_equals_fine_baseline_spline(k):
    """Refine everything: the composite carrier IS the periodic interpolating
    spline of all nodes on the fine grid (two levels deep)."""
    f = Forest(16, 8, 2)
    for (I, J) in sorted(f.exists[0]):
        f.refine(0, I, J)
    for (I, J) in sorted(f.exists[1]):
        f.refine(1, I, J)
    car = HierarchicalCarrier(f, k=k)
    a, sur = car.update(nodal(f), return_surplus=True)
    xi = np.random.default_rng(1).random((3000, 2))
    fine = PSpline2(smooth(grid(64).reshape(-1, 2)).reshape(64, 64, 2), k)
    assert np.abs(car.eval_points(a, xi) - fine(xi)).max() < 2e-14
    assert max(np.abs(s).max() for s in sur[1:]) > 1e-6     # the surpluses are not trivially zero


@pytest.mark.parametrize("k", (3, 5))
@pytest.mark.parametrize("seed", (0, 1))
def test_interpolates_every_active_node(k, seed):
    f = random_forest(32, 8, 3, nsplit=30, seed=seed)
    car = HierarchicalCarrier(f, k=k)
    zs = nodal(f)
    a, sur = car.update(zs, return_surplus=True)
    smax = max(np.abs(s).max() for s in sur[1:])
    err = max(np.abs(car.eval_points(a, nodes_xi(f, l)) - zs[l]).max() for l in range(f.top_level() + 1))
    # exact up to the documented skirt leak  lambda^(halo*B) * surplus
    assert err <= 2.0 * DECAY[k] ** (car.ht * car.B) * smax
    car2 = HierarchicalCarrier(f, k=k, halo_tiles=2)
    a2 = car2.update(zs)
    err2 = max(np.abs(car2.eval_points(a2, nodes_xi(f, l)) - zs[l]).max() for l in range(f.top_level() + 1))
    assert err2 <= 2.0 * DECAY[k] ** (2 * car.B) * smax + 1e-14


@pytest.mark.parametrize("k", (3, 5))
def test_fixed_sweep_solver_matches_exact_solves(k):
    f = random_forest(32, 8, 3, nsplit=30, seed=3)
    zs = nodal(f)
    ref, sur = HierarchicalCarrier(f, k=k, solver="splu").update(zs, return_surplus=True)
    smax = max(np.abs(s).max() for s in sur[1:])
    tol = HierarchicalCarrier(f, k=k, solver="adi-tol", tol=1e-14)
    at = tol.update(zs)
    assert max(np.abs(x - y).max() for x, y in zip(at, ref)) < 1e-12 * max(1.0, smax / 1e-6)
    assert tol.stats.get("cg_fallbacks", 0) == 0
    for sweeps, bound in ((1, 3e-3), (2, 1e-3)):
        a = HierarchicalCarrier(f, k=k, solver="adi", sweeps=sweeps).update(zs)
        assert max(np.abs(x - y).max() for x, y in zip(a, ref)) <= bound * smax


def test_carrier_is_linear_in_the_nodal_data():
    f = random_forest(32, 8, 2, nsplit=15, seed=5)
    car = HierarchicalCarrier(f, k=5)
    rng = np.random.default_rng(0)
    z1 = [rng.standard_normal((len(f.node_keys(l)), 2)) for l in range(f.top_level() + 1)]
    z2 = [rng.standard_normal((len(f.node_keys(l)), 2)) for l in range(f.top_level() + 1)]
    a1, a2 = car.update(z1), car.update(z2)
    a3 = car.update([2.0 * u - 3.0 * v for u, v in zip(z1, z2)])
    assert max(np.abs(c - (2.0 * a - 3.0 * b)).max() for a, b, c in zip(a1, a2, a3)) < 1e-12


@pytest.mark.parametrize("k", (3, 5))
@pytest.mark.parametrize("kind", ("coarse-fine", "all"))
def test_interface_continuity_is_k_minus_1(k, kind):
    """Measured two-sided jumps of F and of its normal derivatives up to order
    k-1 across every coarse-fine (hanging) interface -- and across every
    other block edge -- are at round-off; the k-th derivative genuinely
    jumps: the carrier is C^{k-1} (C^4 for quintic, C^2 for cubic), not C^k."""
    f = patch_forest(32, 8, 2)
    sh = AdaptiveSheet(f, k=k)
    for l in range(len(sh.keys)):
        z = smooth(sh.xi[l])
        sh.d[l], sh.p[l] = 0.05 * z, 0.3 * z[:, ::-1]
    j = dg.continuity_jumps(sh, "xp", kind=kind)
    assert j["n"] > 50
    for m in range(k):
        assert j["jump"][m] <= 1e-6 * j["scale"][m] + 1e-10, (m, j)
    assert j["jump"][k] > 1e-2 * j["scale"][k]


@pytest.mark.parametrize("k,rate_min", ((3, 3.6), (5, 5.5)))
def test_convergence_order_in_patch_and_at_interface(k, rate_min):
    """Same refinement pattern on n0 = 16 and 32 (so every h halves): errors
    in the fine interior, in a strip around the coarse-fine interface and in
    the coarse region all converge at order k+1."""
    errs = []
    rng = np.random.default_rng(0)
    for n0 in (16, 32):
        f = Forest(n0, 8, 1)
        nb = f.nb0
        for I in range(nb // 4, 3 * nb // 4):
            for J in range(nb // 4, 3 * nb // 4):
                f.refine(0, I, J)
        car = HierarchicalCarrier(f, k=k)
        a = car.update(nodal(f))
        xi_in = 0.35 + 0.3 * rng.random((4000, 2))                       # deep inside the patch
        xi_if = np.stack([0.25 + (rng.random(4000) - 0.5) * 0.5 / n0, 0.3 + 0.4 * rng.random(4000)], axis=1)  # +-1/4 coarse cell around the interface
        xi_out = np.stack([0.9 * rng.random(4000) % 0.2, rng.random(4000)], axis=1)  # coarse region
        errs.append([np.abs(car.eval_points(a, z) - smooth(z)).max() for z in (xi_in, xi_if, xi_out)])
    rates = np.log2(np.array(errs[0]) / np.array(errs[1]))
    assert np.all(rates > rate_min), (errs, rates)
    # and the fine interior is (2^(k+1))-ish more accurate than the coarse region
    assert errs[1][0] < errs[1][2] / 2 ** k


def test_split_leaves_the_carrier_unchanged():
    """Transfer = identity: nodes created with the old carrier's values have
    zero surplus, so x(xi), p(xi), S(xi) are the same functions before and
    after a split; old nodes are untouched bit for bit."""
    f = patch_forest(32, 8, 3)
    sh = AdaptiveSheet(f, k=5, periods=(0.2, -0.1))
    for l in range(len(sh.keys)):
        z = smooth(sh.xi[l])
        sh.d[l], sh.p[l] = 0.05 * z, 0.3 * z[:, ::-1]
        sh.S[l] = 0.1 * z[:, 0] * z[:, 1] + sh.xi[l] @ sh.periods
    xi = np.random.default_rng(7).random((6000, 2))
    before = sh.evaluate(xi, "xpS")
    _, sur = sh.coeffs("xpS", return_surplus=True)
    smax = max(np.abs(s).max() for s in sur[1:])
    old = [(k.copy(), d.copy(), p.copy(), S.copy()) for k, d, p, S in zip(sh.keys, sh.d, sh.p, sh.S)]
    rec = sh.refine([(0, 0, 3), (1, 2, 3), (2, 7, 7)], reason="test")
    assert rec["created"] > 0 and rec["removed"] == 0
    after = sh.evaluate(xi, "xpS")
    # identity up to the skirt leak lambda^(halo*B) * surplus (see contract); measured ~1e-11
    bound = 2.0 * DECAY[5] ** (sh.car.ht * sh.car.B) * smax
    assert bound < 1e-8
    for b, a in zip(before, after):
        assert np.abs(a - b).max() <= bound
    for l, (k0, d0, p0, S0) in enumerate(old):
        pos = np.searchsorted(sh.keys[l], k0)
        assert np.array_equal(sh.keys[l][pos], k0)
        assert np.array_equal(sh.d[l][pos], d0) and np.array_equal(sh.p[l][pos], p0) and np.array_equal(sh.S[l][pos], S0)


def test_correction_support_is_local():
    """Storage and influence of a patch are local: coefficient tiles of level
    1 cover a bounded neighbourhood of the patch, and perturbing a fine node
    changes the carrier only near it (exponential decay), which a global
    regrid / global fine spline would not do with this storage."""
    f = Forest(64, 8, 2)                          # 8 x 8 blocks
    f.refine(0, 3, 3)
    car = HierarchicalCarrier(f, k=5)
    lv = car.lev[1]
    assert lv.nT == (4 * 8) ** 2                  # patch (2 x 2 fine tiles) + one halo ring of tiles
    assert lv.nD == (8 * 8) ** 2 and lv.nD <= (2 * 64) ** 2 // 4   # bounded neighbourhood, not the fine grid
    zs = nodal(f)
    a0 = car.update(zs)
    zs[1] = zs[1].copy()
    zs[1][len(zs[1]) // 2] += 1.0                 # unit bump on one fine node
    a1 = car.update(zs)
    xi = np.random.default_rng(0).random((20000, 2))
    dv = np.abs(car.eval_points(a1, xi) - car.eval_points(a0, xi)).max(axis=1)
    N = f.ncell(1)
    kk = f.node_keys(1)[len(zs[1]) // 2]
    c = np.array([kk // N, kk % N]) / N
    dist = np.abs((xi - c + 0.5) % 1.0 - 0.5).max(axis=1) * N       # in fine cells
    assert dv[dist > 12].max() < 1e-4 and dv[dist > 24].max() == 0.0
    assert np.abs(a1[0] - a0[0]).max() == 0.0     # a fine node never changes the coarse level


def test_block_tensor_evaluation_matches_point_evaluation():
    f = random_forest(32, 8, 3, nsplit=25, seed=6)
    car = HierarchicalCarrier(f, k=5)
    a = car.update(nodal(f))
    B = car.B
    for (L, d), (ids, _) in car.groups.items():
        W = car.weights(d, 2)
        ids2, V = car.block_values(a, L, d, W, W)
        N = f.ncell(L)
        u = (np.arange(2 * B) + 0.5) / 2
        x1 = (car.leafI[ids2][:, None] * B + u[None, :]) / N
        x2 = (car.leafJ[ids2][:, None] * B + u[None, :]) / N
        xi = np.stack([np.broadcast_to(x1[:, :, None], V.shape[:3]), np.broadcast_to(x2[:, None, :], V.shape[:3])], axis=-1).reshape(-1, 2)
        assert np.abs(V.reshape(-1, 2) - car.eval_points(a, xi)).max() < 1e-13
