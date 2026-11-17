"""Shared fixtures for the regression tests (no reference data anywhere)."""
import numpy as np

from shadapt.forest import Forest

TWO_PI = 2.0 * np.pi


def smooth(xi):
    """Smooth periodic 2-component test field."""
    return np.stack([np.sin(TWO_PI * xi[:, 0]) * np.cos(TWO_PI * 2 * xi[:, 1] + 0.3),
                     np.exp(np.cos(TWO_PI * (xi[:, 0] + xi[:, 1])))], axis=1)


def nodes_xi(f, l):
    N = f.ncell(l)
    k = f.node_keys(l)
    return np.stack([k // N, k % N], axis=1) / N


def nodal(f, fun=smooth):
    return [fun(nodes_xi(f, l)) for l in range(f.top_level() + 1)]


def random_forest(n0=32, B=8, lmax=3, nsplit=25, seed=0):
    rng = np.random.default_rng(seed)
    f = Forest(n0, B, lmax)
    for _ in range(nsplit):
        lv = f.leaves()
        l, I, J = lv[rng.integers(len(lv))]
        if l < f.lmax:
            f.refine(l, I, J, reason="test")
    f.check()
    return f


def patch_forest(n0=32, B=8, lmax=2):
    """One refined patch of 2x2 level-0 blocks with a nested level-2 block."""
    f = Forest(n0, B, lmax)
    for (I, J) in ((1, 1), (2, 1), (1, 2), (2, 2)):
        f.refine(0, I, J, reason="test")
    if lmax >= 2:
        f.refine(1, 3, 3, reason="test")
        f.refine(1, 4, 4, reason="test")
    f.check()
    return f
