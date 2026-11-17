"""Offline diagnostics: carrier loop actions, interface continuity probes,
and errors against a grid-matched control.

Functions that take a `ref` argument are EVALUATION ONLY and are never called
by the solver or the controller.
"""
from __future__ import annotations

from math import gcd

import numpy as np


# ---------------------------------------------------------------------------
# state-only structural diagnostics
# ---------------------------------------------------------------------------
def loop_actions(sheet, nlines=16, oversample=4):
    """Carrier loop integrals  c_a(line) = oint p_h . d x_h  along `nlines`
    equally spaced label lines in each direction (the fundamental cycles of
    the label torus).  For the exact sheet every such integral equals the
    conserved Liouville period.  Returns dict(c1, c2 arrays, max deviation
    from sheet.periods, spread)."""
    f = sheet.f
    Nf = f.ncell(f.top_level())
    m = Nf * oversample
    s = (np.arange(m) + 0.5) / m
    out = {}
    for ax in (0, 1):
        vals = []
        for j in range(nlines):
            xi = np.empty((m, 2))
            xi[:, ax] = s
            xi[:, 1 - ax] = j / nlines
            (x, gx), (p, _) = sheet.evaluate(xi, "xp", deriv=1)
            vals.append(float(np.mean(np.sum(p * gx[:, :, ax], axis=1))))
        out[f"c{ax + 1}"] = np.array(vals)
    dev = max(float(np.abs(out["c1"] - sheet.periods[0]).max()), float(np.abs(out["c2"] - sheet.periods[1]).max()))
    out["max_deviation"] = dev
    out["spread"] = max(float(np.ptp(out["c1"])), float(np.ptp(out["c2"])))
    return out


def nodal_action_defect(sheet, nsamp=4000, seed=0):
    """max |S_h(xi) - S_node-consistent| proxy: RMS over random label segments
    of  [S_h(b) - S_h(a)] - int_a^b p_h . dx_h  (carrier action vs carrier
    Liouville form along short straight label segments)."""
    rng = np.random.default_rng(seed)
    f = sheet.f
    h = 1.0 / f.ncell(0)
    a = rng.random((nsamp, 2))
    dirs = rng.integers(0, 2, nsamp)
    gl, gw = np.polynomial.legendre.leggauss(6)
    gl, gw = 0.5 * (gl + 1.0), 0.5 * gw
    b = a.copy()
    b[np.arange(nsamp), dirs] += h
    (Sa,) = sheet.evaluate(a, "S")
    (Sb,) = sheet.evaluate(b, "S")
    integ = np.zeros(nsamp)
    for t, w in zip(gl, gw):
        xi = a + t * (b - a)
        (x, gx), (p, _) = sheet.evaluate(xi, "xp", deriv=1)
        integ += w * h * np.sum(p * gx[np.arange(nsamp), :, dirs], axis=1)
    dfc = (Sb - Sa) - integ
    return dict(rms=float(np.sqrt(np.mean(dfc ** 2))), max=float(np.abs(dfc).max()))


def interface_points(forest, per_edge=3, kind="coarse-fine"):
    """Label points on block edges: returns (xi, normal axis, level).
    kind "coarse-fine": edges between a leaf block of level l and a coarser
    region (the hanging interfaces); "all": every edge of every leaf block
    (also probes the places where the evaluation level of the composite
    carrier changes)."""
    pts, axes, lev = [], [], []
    B = forest.B
    for (l, I, J) in forest.leaves():
        if l == 0 and kind == "coarse-fine":
            continue
        nb = forest.nbl(l)
        N = forest.ncell(l)
        for ax, sgn in ((0, -1), (0, 1), (1, -1), (1, 1)):
            q = [I, J]
            q[ax] = (q[ax] + sgn) % nb
            if kind == "coarse-fine" and tuple(q) in forest.exists[l]:
                continue
            for tt in (np.arange(per_edge) + 0.37) / per_edge:
                p = [0.0, 0.0]
                p[ax] = ((I if ax == 0 else J) + (1 if sgn > 0 else 0)) * B / N
                p[1 - ax] = ((J if ax == 0 else I) + tt) * B / N
                pts.append(p)
                axes.append(ax)
                lev.append(l)
    return np.array(pts).reshape(-1, 2), np.array(axes, dtype=int), np.array(lev, dtype=int)


def continuity_jumps(sheet, which="xp", kind="coarse-fine", per_edge=3):
    """Two-sided jumps of the carrier and of its normal derivatives of order
    0..k across block edges, WITHOUT difference-quotient truncation error:
    along the normal line, inside one cell of the finest level on either side
    of the edge, every component of F is exactly a polynomial of degree k, so
    it is recovered from k+1 samples per side and the two polynomials are
    compared at the edge.  Returns, for m = 0..k, max |[d^m F/dn^m]| h^m and
    max |d^m F/dn^m| h^m (h = finest cell).  A C^{k-1} carrier has jumps at
    round-off for m <= k-1 and a genuine jump for m = k."""
    from math import factorial
    f = sheet.f
    k = sheet.k
    xi, axes, _ = interface_points(f, per_edge, kind)
    if len(xi) == 0:
        return dict(n=0)
    h = 1.0 / f.ncell(f.top_level())
    s = (np.arange(k + 1) + 0.5) / (k + 1) * 0.9             # sample abscissae in units of h
    coef = {}
    for side in (-1, 1):
        vals = []
        for sj in s:
            z = xi.copy()
            z[np.arange(len(xi)), axes] += side * sj * h
            vals.append(np.concatenate([np.asarray(c).reshape(len(xi), -1) for c in sheet.evaluate(z, which)], axis=1))
        V = np.vander(side * s, k + 1, increasing=True)        # same signed coordinate on both sides
        coef[side] = np.linalg.solve(V, np.stack(vals, axis=0).reshape(k + 1, -1)).reshape(k + 1, len(xi), -1)
    out = dict(n=int(len(xi)), jump=[], scale=[])
    for m in range(k + 1):
        out["jump"].append(float(factorial(m) * np.abs(coef[1][m] - coef[-1][m]).max()))
        out["scale"].append(float(factorial(m) * np.abs(coef[1][m]).max()))
    return out


# ---------------------------------------------------------------------------
# offline comparison with a grid-matched control
# ---------------------------------------------------------------------------
def power_spectrum(rho, nbins=None):
    ng = rho.shape[0]
    dk = np.fft.rfft2(rho / rho.mean() - 1.0) / (ng * ng)
    kx = np.fft.fftfreq(ng, d=1.0 / ng)
    ky = np.fft.rfftfreq(ng, d=1.0 / ng)
    kk = np.sqrt(kx[:, None] ** 2 + ky[None, :] ** 2)
    w = np.ones_like(kk)
    w[:, 1:-1] = 2.0                       # rfft half-plane multiplicity
    edges = np.arange(0.5, ng // 2 + 1.0, 1.0)
    idx = np.digitize(kk.ravel(), edges)
    P = np.bincount(idx, weights=(np.abs(dk) ** 2 * w).ravel(), minlength=len(edges) + 1)[1:len(edges)]
    cnt = np.bincount(idx, weights=w.ravel(), minlength=len(edges) + 1)[1:len(edges)]
    return 0.5 * (edges[:-1] + edges[1:]), P / np.maximum(cnt, 1.0)


def _lattice_common(n, nref):
    g = gcd(n, nref)
    return n // g, nref // g, g


def compare_fields(pm, rho, ref_rho):
    """Density, potential-gradient (force) and spectrum errors on the PM grid."""
    out = dict(rho_rel_l2=float(np.linalg.norm(rho - ref_rho) / np.linalg.norm(ref_rho - 1.0)),
               rho_max=float(rho.max()), rho_max_ratio=float(rho.max() / ref_rho.max()))
    a = pm.grid_acceleration(pm.potential(rho))
    ar = pm.grid_acceleration(pm.potential(ref_rho))
    da = np.sqrt(np.sum((a - ar) ** 2, axis=-1))
    nr = np.sqrt(np.sum(ar ** 2, axis=-1))
    out["force_rel_l2"] = float(np.sqrt(np.mean(da ** 2)) / np.sqrt(np.mean(nr ** 2)))
    out["force_max_over_rms"] = float(da.max() / np.sqrt(np.mean(nr ** 2)))
    hi = ref_rho > 4.0                      # collapsed structure (caustics / halo)
    out["force_rel_l2_dense"] = float(np.sqrt(np.mean(da[hi] ** 2)) / np.sqrt(np.mean(nr[hi] ** 2))) if hi.any() else None
    kc, P = power_spectrum(rho)
    _, Pr = power_spectrum(ref_rho)
    sel = kc <= pm.ng / 4
    rel = np.abs(P[sel] - Pr[sel]) / Pr[sel]
    out["pk_rel_max_to_quarter_nyq"] = float(rel.max())
    out["pk_rel_rms_to_quarter_nyq"] = float(np.sqrt(np.mean(rel ** 2)))
    sel2 = kc <= pm.ng / 2
    rel2 = np.abs(P[sel2] - Pr[sel2]) / Pr[sel2]
    out["pk_rel_max_to_nyq"] = float(rel2.max())
    return out


def compare_nbody(run, ref, k):
    """N-body arm vs control at the common label lattice (particle labels)."""
    n = run.n
    nref = ref[f"x{k}"].shape[0]
    s, sr, g = _lattice_common(n, nref)
    x = run.x.reshape(n, n, 2)[::s, ::s]
    p = run.p.reshape(n, n, 2)[::s, ::s]
    ex, ep = x - ref[f"x{k}"][::sr, ::sr], p - ref[f"p{k}"][::sr, ::sr]
    out = dict(labels=int(g * g), x_rms=float(np.sqrt(np.mean(np.sum(ex ** 2, -1)))), x_max=float(np.abs(ex).max()),
               p_rms=float(np.sqrt(np.mean(np.sum(ep ** 2, -1)))), p_max=float(np.abs(ep).max()))
    out.update(compare_fields(run.pm, run.rho, ref[f"rho{k}"]))
    return out


def compare_sheet(run, ref, k, lattice=256):
    """Sheet arm vs control: (a) the CARRIER evaluated on a uniform label
    lattice (the sheet as a continuum representation: the number used for
    the gate), (b) the level-0 nodes only."""
    sh = run.sh
    nref = ref[f"x{k}"].shape[0]
    lat = min(lattice, nref)
    sr = nref // lat
    g = np.arange(lat) / lat
    xi = np.stack(np.meshgrid(g, g, indexing="ij"), axis=-1).reshape(-1, 2)
    x, p = sh.evaluate(xi, "xp")
    xr, pr = ref[f"x{k}"][::sr, ::sr].reshape(-1, 2), ref[f"p{k}"][::sr, ::sr].reshape(-1, 2)
    ex, ep = x - xr, p - pr
    out = dict(labels=int(lat * lat), x_rms=float(np.sqrt(np.mean(np.sum(ex ** 2, -1)))), x_max=float(np.abs(ex).max()),
               p_rms=float(np.sqrt(np.mean(np.sum(ep ** 2, -1)))), p_max=float(np.abs(ep).max()))
    n0 = sh.f.n0
    s, sr0, g0 = _lattice_common(n0, nref)
    x0 = (sh.xi[0] + sh.d[0]).reshape(n0, n0, 2)[::s, ::s]
    p0 = sh.p[0].reshape(n0, n0, 2)[::s, ::s]
    e0, q0 = x0 - ref[f"x{k}"][::sr0, ::sr0], p0 - ref[f"p{k}"][::sr0, ::sr0]
    out.update(node0_labels=int(g0 * g0), node0_x_rms=float(np.sqrt(np.mean(np.sum(e0 ** 2, -1)))),
               node0_p_rms=float(np.sqrt(np.mean(np.sum(q0 ** 2, -1)))))
    out.update(compare_fields(run.pm, run.rho, ref[f"rho{k}"]))
    return out
