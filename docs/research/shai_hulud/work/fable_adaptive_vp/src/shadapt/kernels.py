"""Low-level numerical kernels shared by every arm (sheet and N-body alike).

* B-spline tables (degree 3 and 5): collocation taps, two-scale (subdivision)
  masks, point-evaluation weights.  The polynomial pieces are taken from the
  vendored baseline (`shsheet.pspline._BASIS`, exact rational arithmetic), so
  the level-0 carrier is the baseline carrier.
* PCS particle-mesh deposit and exact-gradient ("grad" law) gather, compiled
  with numba.  They implement the formulas of the baseline `shsheet.pm.PM`
  (order 4, force="grad") and are tested against it to round-off; they are
  serial and deterministic (fixed summation order = input order).
"""
from __future__ import annotations

from math import comb

import numba
import numpy as np

from shsheet.pspline import _BASIS

TWO_PI = 2.0 * np.pi

# values of the centred cardinal B-spline at the integers (collocation taps)
COLLOC = {3: np.array([1.0, 4.0, 1.0]) / 6.0, 5: np.array([1.0, 26.0, 66.0, 26.0, 1.0]) / 120.0}
# dominant root of the collocation symbol: decay rate of the cardinal spline per knot
DECAY = {3: 2.0 - np.sqrt(3.0), 5: 0.43057535}


def two_scale_mask(k):
    """h_j, j = -(k+1)/2 .. (k+1)/2, with  beta(u) = sum_j h_j beta(2u - j)."""
    return np.array([comb(k + 1, j) for j in range(k + 2)], dtype=float) / 2.0 ** k


def eval_weights(t, k, nder=0):
    """Weights w[a] (and derivatives w.r.t. the grid coordinate) of the k+1
    B-splines that are non-zero at fractional offset t in [0, 1):
    s(i0 + t) = sum_a w[a](t) c[i0 + js[a]].  Returns (js, [w, dw, ...]) with
    each array of shape (k+1, len(t))."""
    js, P = _BASIS[k]
    t = np.asarray(t, dtype=float)
    pw = np.ones((k + 1,) + t.shape)
    for m in range(1, k + 1):
        pw[m] = pw[m - 1] * t
    out = []
    C = P.copy()
    for _ in range(nder + 1):
        out.append(np.tensordot(C, pw, axes=(1, 0)))
        Cn = np.zeros_like(C)
        Cn[:, :-1] = C[:, 1:] * np.arange(1, k + 1)
        C = Cn
    return js, out


def prefilter_periodic(f, k):
    """Coefficients of the periodic interpolating tensor spline of degree k
    through f[i, j, ...] (the baseline `_prefilter` operator, evaluated with a
    real 2D FFT; agrees with it to round-off)."""
    if k == 1:
        return f.copy()
    n = f.shape[0]
    taps = COLLOC[k]
    h = (len(taps) - 1) // 2

    def symbol(m):
        lam = np.full(len(m), taps[h])
        for j in range(1, h + 1):
            lam = lam + 2.0 * taps[h + j] * np.cos(TWO_PI * j * m / n)
        return lam
    lam2 = symbol(np.arange(n))[:, None] * symbol(np.arange(n // 2 + 1))[None, :]
    c = np.fft.rfft2(f, axes=(0, 1))
    c /= lam2.reshape(lam2.shape + (1,) * (f.ndim - 2))
    return np.fft.irfft2(c, s=(n, n), axes=(0, 1))


# ---------------------------------------------------------------------------
# particle-mesh kernels (PCS = order 4), 2D periodic unit box
# ---------------------------------------------------------------------------
@numba.njit(cache=True, fastmath=False)
def _pcs_w(t, w):
    w[0] = (1.0 - t) ** 3 / 6.0
    w[1] = (3.0 * t ** 3 - 6.0 * t ** 2 + 4.0) / 6.0
    w[2] = (-3.0 * t ** 3 + 3.0 * t ** 2 + 3.0 * t + 1.0) / 6.0
    w[3] = t ** 3 / 6.0


@numba.njit(cache=True, fastmath=False)
def _pcs_dw(t, dw):
    dw[0] = -(1.0 - t) ** 2 / 2.0
    dw[1] = (9.0 * t ** 2 - 12.0 * t) / 6.0
    dw[2] = (-9.0 * t ** 2 + 6.0 * t + 3.0) / 6.0
    dw[3] = t ** 2 / 2.0


@numba.njit(cache=True, fastmath=False)
def deposit_pcs(rho, x, mass):
    """rho[ng, ng] += sum of PCS-assigned masses (no ng^2 normalisation).
    x (M, 2) any real values (periodic); mass scalar applied to every point."""
    ng = rho.shape[0]
    wa = np.empty(4)
    wb = np.empty(4)
    ia = np.empty(4, dtype=np.int64)
    ib = np.empty(4, dtype=np.int64)
    for m in range(x.shape[0]):
        u = x[m, 0] * ng
        v = x[m, 1] * ng
        fu = np.floor(u)
        fv = np.floor(v)
        _pcs_w(u - fu, wa)
        _pcs_w(v - fv, wb)
        i0 = int(fu) - 1
        j0 = int(fv) - 1
        for a in range(4):
            ia[a] = (i0 + a) % ng
            ib[a] = (j0 + a) % ng
        for a in range(4):
            wm = wa[a] * mass
            for b in range(4):
                rho[ia[a], ib[b]] += wm * wb[b]


@numba.njit(cache=True, fastmath=False)
def deposit_pcs_w(rho, x, mass):
    """As deposit_pcs with one mass per point (mass (M,))."""
    ng = rho.shape[0]
    wa = np.empty(4)
    wb = np.empty(4)
    ia = np.empty(4, dtype=np.int64)
    ib = np.empty(4, dtype=np.int64)
    for m in range(x.shape[0]):
        u = x[m, 0] * ng
        v = x[m, 1] * ng
        fu = np.floor(u)
        fv = np.floor(v)
        _pcs_w(u - fu, wa)
        _pcs_w(v - fv, wb)
        i0 = int(fu) - 1
        j0 = int(fv) - 1
        for a in range(4):
            ia[a] = (i0 + a) % ng
            ib[a] = (j0 + a) % ng
        for a in range(4):
            wm = wa[a] * mass[m]
            for b in range(4):
                rho[ia[a], ib[b]] += wm * wb[b]


@numba.njit(cache=True, fastmath=False)
def gather_pcs_grad(phi, x, V, acc):
    """V[m] = sum_g phi_g W(x_m - x_g),  acc[m] = -grad V(x_m)  (exact gradient
    of the PCS-interpolated potential: the baseline "grad" force law)."""
    ng = phi.shape[0]
    wa = np.empty(4)
    wb = np.empty(4)
    da = np.empty(4)
    db = np.empty(4)
    ia = np.empty(4, dtype=np.int64)
    ib = np.empty(4, dtype=np.int64)
    for m in range(x.shape[0]):
        u = x[m, 0] * ng
        v = x[m, 1] * ng
        fu = np.floor(u)
        fv = np.floor(v)
        _pcs_w(u - fu, wa)
        _pcs_w(v - fv, wb)
        _pcs_dw(u - fu, da)
        _pcs_dw(v - fv, db)
        i0 = int(fu) - 1
        j0 = int(fv) - 1
        for a in range(4):
            ia[a] = (i0 + a) % ng
            ib[a] = (j0 + a) % ng
        val = 0.0
        g0 = 0.0
        g1 = 0.0
        for a in range(4):
            for b in range(4):
                f = phi[ia[a], ib[b]]
                val += wa[a] * wb[b] * f
                g0 += da[a] * wb[b] * f
                g1 += wa[a] * db[b] * f
        V[m] = val
        acc[m, 0] = -g0 * ng
        acc[m, 1] = -g1 * ng


@numba.njit(cache=True, fastmath=False)
def gather_pcs_hess(phi, x, H):
    """Hessian of the PCS-interpolated potential at x: H[m] = (V_xx, V_xy, V_yy).
    Piecewise linear (second derivative of a cubic kernel)."""
    ng = phi.shape[0]
    wa = np.empty(4)
    wb = np.empty(4)
    da = np.empty(4)
    db = np.empty(4)
    ha = np.empty(4)
    hb = np.empty(4)
    for m in range(x.shape[0]):
        u = x[m, 0] * ng
        v = x[m, 1] * ng
        fu = np.floor(u)
        fv = np.floor(v)
        tu = u - fu
        tv = v - fv
        _pcs_w(tu, wa)
        _pcs_w(tv, wb)
        _pcs_dw(tu, da)
        _pcs_dw(tv, db)
        ha[0] = 1.0 - tu
        ha[1] = 3.0 * tu - 2.0
        ha[2] = -3.0 * tu + 1.0
        ha[3] = tu
        hb[0] = 1.0 - tv
        hb[1] = 3.0 * tv - 2.0
        hb[2] = -3.0 * tv + 1.0
        hb[3] = tv
        i0 = int(fu) - 1
        j0 = int(fv) - 1
        hxx = 0.0
        hxy = 0.0
        hyy = 0.0
        for a in range(4):
            ii = (i0 + a) % ng
            for b in range(4):
                f = phi[ii, (j0 + b) % ng]
                hxx += ha[a] * wb[b] * f
                hxy += da[a] * db[b] * f
                hyy += wa[a] * hb[b] * f
        H[m, 0] = hxx * ng * ng
        H[m, 1] = hxy * ng * ng
        H[m, 2] = hyy * ng * ng


@numba.njit(cache=True, fastmath=False)
def deposit_pcs_blocks(rho, V, mass):
    """Deposit the tensor sample lattices of a group of blocks without copies:
    V[nblk, 2, n1, n2] holds x (component 0) and y (component 1)."""
    ng = rho.shape[0]
    wa = np.empty(4)
    wb = np.empty(4)
    ia = np.empty(4, dtype=np.int64)
    ib = np.empty(4, dtype=np.int64)
    for q in range(V.shape[0]):
        for r in range(V.shape[2]):
            for c in range(V.shape[3]):
                u = V[q, 0, r, c] * ng
                v = V[q, 1, r, c] * ng
                fu = np.floor(u)
                fv = np.floor(v)
                _pcs_w(u - fu, wa)
                _pcs_w(v - fv, wb)
                i0 = int(fu) - 1
                j0 = int(fv) - 1
                for a in range(4):
                    ia[a] = (i0 + a) % ng
                    ib[a] = (j0 + a) % ng
                for a in range(4):
                    wm = wa[a] * mass
                    for b in range(4):
                        rho[ia[a], ib[b]] += wm * wb[b]


@numba.njit(cache=True, fastmath=False)
def block_spacing(V, out):
    """out[q, 0/1] = largest max-norm distance between neighbouring samples of
    block q along the first / second lattice direction."""
    for q in range(V.shape[0]):
        m1 = 0.0
        m2 = 0.0
        for r in range(V.shape[2]):
            for c in range(V.shape[3]):
                if r + 1 < V.shape[2]:
                    d = max(abs(V[q, 0, r + 1, c] - V[q, 0, r, c]), abs(V[q, 1, r + 1, c] - V[q, 1, r, c]))
                    if d > m1:
                        m1 = d
                if c + 1 < V.shape[3]:
                    d = max(abs(V[q, 0, r, c + 1] - V[q, 0, r, c]), abs(V[q, 1, r, c + 1] - V[q, 1, r, c]))
                    if d > m2:
                        m2 = d
        out[q, 0] = m1
        out[q, 1] = m2


@numba.njit(cache=True, fastmath=False)
def kick_pcs(phi, x, p, kappa):
    """In-place kick  p += kappa * (-grad V)(x)  with the exact gradient of
    the PCS-interpolated potential (same force as gather_pcs_grad; the
    potential value itself is not accumulated)."""
    ng = phi.shape[0]
    wa = np.empty(4)
    wb = np.empty(4)
    da = np.empty(4)
    db = np.empty(4)
    ia = np.empty(4, dtype=np.int64)
    ib = np.empty(4, dtype=np.int64)
    kn = kappa * ng
    for m in range(x.shape[0]):
        u = x[m, 0] * ng
        v = x[m, 1] * ng
        fu = np.floor(u)
        fv = np.floor(v)
        _pcs_w(u - fu, wa)
        _pcs_w(v - fv, wb)
        _pcs_dw(u - fu, da)
        _pcs_dw(v - fv, db)
        i0 = int(fu) - 1
        j0 = int(fv) - 1
        for a in range(4):
            ia[a] = (i0 + a) % ng
            ib[a] = (j0 + a) % ng
        g0 = 0.0
        g1 = 0.0
        for a in range(4):
            for b in range(4):
                f = phi[ia[a], ib[b]]
                g0 += da[a] * wb[b] * f
                g1 += wa[a] * db[b] * f
        p[m, 0] -= kn * g0
        p[m, 1] -= kn * g1


@numba.njit(cache=True, fastmath=False)
def drift_inplace(x, p, dt):
    for m in range(x.shape[0]):
        x[m, 0] += dt * p[m, 0]
        x[m, 1] += dt * p[m, 1]


class PMGrid:
    """2D periodic PM field with the baseline conventions (box 1, total mass 1,
    4 pi G = 1:  lap(phi) = rho - 1), PCS assignment, exact-gradient force.
    Same Green function as `shsheet.pm.PM` (continuous -1/k^2, zero mode 0)."""

    order = 4
    force = "grad"

    def __init__(self, ng):
        self.ng = ng
        k = TWO_PI * np.fft.fftfreq(ng, d=1.0 / ng)
        kr = TWO_PI * np.fft.rfftfreq(ng, d=1.0 / ng)
        k1, k2 = np.meshgrid(k, kr, indexing="ij")
        k2s = k1 ** 2 + k2 ** 2
        k2s[0, 0] = 1.0
        self.green_r = -1.0 / k2s
        self.green_r[0, 0] = 0.0

    def new_density(self):
        return np.zeros((self.ng, self.ng))

    def finalize_density(self, acc):
        """Accumulated masses -> density with mean = total mass (ng^2 factor)."""
        return acc * float(self.ng * self.ng)

    def potential(self, rho):
        return np.fft.irfft2(np.fft.rfft2(rho - rho.mean()) * self.green_r, s=rho.shape)

    def field(self, phi, x):
        V = np.empty(x.shape[0])
        a = np.empty((x.shape[0], 2))
        gather_pcs_grad(phi, np.ascontiguousarray(x), V, a)
        return V, a

    def energy_potential(self, rho, phi):
        return 0.5 * float(np.mean((rho - rho.mean()) * phi))

    def grid_acceleration(self, phi):
        """Exact gradient of the interpolated potential AT the grid points
        (diagnostic force field for force-error norms): (ng, ng, 2)."""
        ng = self.ng
        g = (np.arange(ng) / ng)
        X = np.stack(np.meshgrid(g, g, indexing="ij"), axis=-1).reshape(-1, 2)
        return self.field(phi, X)[1].reshape(ng, ng, 2)
