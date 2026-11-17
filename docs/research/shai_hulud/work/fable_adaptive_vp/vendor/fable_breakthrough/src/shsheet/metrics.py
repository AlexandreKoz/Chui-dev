"""Offline diagnostics.  The oracle appears ONLY here and in experiment
scripts, never in reconstruction/adaptation code."""
from __future__ import annotations

import numpy as np

from .conventions import omega


def _rms(v):
    return float(np.sqrt(np.mean(np.sum(v * v, axis=-1)))) if v.ndim > 1 else float(np.sqrt(np.mean(v * v)))


def physical_errors(x, p, S, xo, po, So, T=None, lp=1.0):
    """Position / momentum / action errors against oracle samples at the same
    labels, plus normal / tangential split in the scaled metric (x, p/lp) when
    the true tangent T (M, 4, 2) is given."""
    bad = ~(np.isfinite(x).all(axis=1) & np.isfinite(p).all(axis=1))
    x, p = np.where(bad[:, None], 1e30, x), np.where(bad[:, None], 1e30, p)
    out = dict(x_rms=_rms(x - xo), x_max=float(np.abs(x - xo).max()), p_rms=_rms(p - po), p_max=float(np.abs(p - po).max()),
               n_bad=int(bad.sum()))
    if S is not None:
        dS = S - So
        out["S_rms"] = _rms(dS - dS.mean())
    if T is not None:
        e = np.concatenate([x - xo, (p - po) / lp], axis=1)
        Ts = T.copy()
        Ts[:, 2:, :] /= lp
        Q, _ = np.linalg.qr(Ts)
        et = np.einsum("mca,mc->ma", Q, e)
        en = e - np.einsum("mca,ma->mc", Q, et)
        out.update(normal_rms=_rms(en), normal_max=float(np.sqrt(np.sum(en * en, axis=1)).max()), tangential_rms=_rms(et))
    return out


# -- spectral line integrals on smooth sub-segments --------------------------------
def _cheb(m):
    k = np.arange(m)
    u = np.cos(np.pi * k / (m - 1))
    c = np.ones(m)
    c[0] = c[-1] = 2.0
    c *= (-1.0) ** k
    U = np.tile(u, (m, 1)).T
    dU = U - U.T
    D = np.outer(c, 1.0 / c) / (dU + np.eye(m))
    D -= np.diag(D.sum(axis=1))
    # Clenshaw-Curtis weights
    w = np.zeros(m)
    N = m - 1
    for j in range(m):
        s = 0.0
        for r in range(1, N // 2 + 1):
            b = 1.0 if 2 * r == N else 2.0
            s += b / (4 * r * r - 1) * np.cos(2 * r * j * np.pi / N)
        w[j] = (2.0 / N) * (1.0 - s) * (0.5 if j in (0, N) else 1.0)
    return u, D, w


_U, _D, _WCC = _cheb(18)


def _segment_breaks(a, b, n):
    """Parameters in [0,1] where the segment a->b crosses grid lines k/n."""
    ts = [0.0, 1.0]
    for ax in range(2):
        lo, hi = sorted((a[ax], b[ax]))
        for k in range(int(np.ceil(lo * n - 1e-12)), int(np.floor(hi * n + 1e-12)) + 1):
            if abs(b[ax] - a[ax]) > 1e-15:
                t = (k / n - a[ax]) / (b[ax] - a[ax])
                if 1e-12 < t < 1 - 1e-12:
                    ts.append(t)
    ts = np.unique(np.round(ts, 14))
    return ts


def line_integral(car, a, b, n):
    """int p.dx along the carrier image of the straight label segment a->b,
    split at cell boundaries; spectral (Chebyshev) on each smooth piece.
    Also returns Sigma_car(b) - Sigma_car(a) accumulated piecewise."""
    ts = _segment_breaks(a, b, n)
    tot, dS = 0.0, 0.0
    for t0, t1 in zip(ts[:-1], ts[1:]):
        # shrink by a hair so that every point lies in the same cell
        tau = t0 + (t1 - t0) * (0.5 * (_U + 1.0) * (1 - 2e-13) + 1e-13)
        xi = a[None, :] + tau[:, None] * (b - a)[None, :]
        x, p, S = car(xi)
        tot += float(np.sum(_WCC[:, None] * p * (_D @ x)))
        if S is not None:
            dS += float(S[0] - S[-1])  # _U runs from +1 to -1
    return tot, dS


def loop_circulation(car, verts, n):
    """oint p.dx around the label polygon verts (closed automatically)."""
    tot, dS = 0.0, 0.0
    for a, b in zip(verts, np.roll(verts, -1, axis=0)):
        c, s = line_integral(car, a, b, n)
        tot += c
        dS += s
    return tot, dS


def random_loops(rng, count, n, rmin=0.6, rmax=6.0):
    """Random label quadrilaterals with radii between rmin and rmax cells."""
    out = []
    for _ in range(count):
        c = rng.random(2)
        r = (rmin + (rmax - rmin) * rng.random()) / n
        th = np.sort(rng.random(4) * 2 * np.pi)
        out.append(c[None, :] + r * np.stack([np.cos(th), np.sin(th)], axis=1))
    return out


def two_form_defect(car, n, rng, count=400, m=10, size=0.8):
    """Pointwise omega(d_1 F, d_2 F) of the carrier at the centres of random
    sub-squares strictly inside cells, by Chebyshev spectral differentiation
    (m x m probe).  Normalised by |d_1 F| |d_2 F|."""
    u, D, _ = _cheb(m)
    cells = rng.integers(0, n, size=(count, 2))
    off = 0.5 * (1 - size) + 0.5 * size * (u + 1.0)  # local coords in (0,1)
    s1, s2 = np.meshgrid(off, off, indexing="ij")
    xi = (cells[:, None, None, :] + np.stack([s1, s2], axis=-1)[None]) / n
    x, p, _ = car(xi.reshape(-1, 2))
    z = np.concatenate([x, p], axis=1).reshape(count, m, m, 4)
    sc = 2.0 / (size / n)  # d/dxi = (2 / width) d/du
    d1 = np.einsum("ab,cbjk->cajk", D, z) * sc
    d2 = np.einsum("ab,cjbk->cjak", D, z) * sc
    w = omega(d1, d2)
    nrm = np.linalg.norm(d1, axis=-1) * np.linalg.norm(d2, axis=-1)
    mid = m // 2
    return np.abs(w[:, 1:-1, 1:-1] / nrm[:, 1:-1, 1:-1]).max(axis=(1, 2)), np.abs(w[:, mid, mid])


def seam_jump(car, n, rng, count=4000):
    """Max jump of (x, p, Sigma) across cell edges evaluated from both sides
    (only for carriers exposing eval_local)."""
    I = rng.integers(0, n, size=(count, 2))
    t = rng.random(count)
    jx = jp = jS = 0.0
    for ax in (0, 1):
        s_hi = np.zeros((count, 2)); s_hi[:, ax] = 1.0; s_hi[:, 1 - ax] = t
        s_lo = np.zeros((count, 2)); s_lo[:, 1 - ax] = t
        e = np.zeros(2, dtype=np.int64); e[ax] = 1
        xa, pa, Sa = car.eval_local(I, s_hi)
        xb, pb, Sb = car.eval_local(I + e, s_lo)
        jx = max(jx, float(np.abs(xa - xb).max())); jp = max(jp, float(np.abs(pa - pb).max()))
        jS = max(jS, float(np.abs(Sa - Sb).max()))
    return dict(x=jx, p=jp, S=jS)
