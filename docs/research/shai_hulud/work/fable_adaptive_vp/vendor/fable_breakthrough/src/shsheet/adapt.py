"""Label maps, state-only monitors, carrier-based remap and h-refinement.

Nothing in this module may touch the oracle: monitors are functions of the
nodal state, maps are functions of monitors, and new nodal values are read
from a carrier built from the old nodal state.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from .carriers import make_carrier
from .pspline import PSpline2, grid
from .sheet import Sheet

_VENDOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "vendor", "fable_prototype")
_GL = np.polynomial.legendre.leggauss(6)
_GLX, _GLW = 0.5 * (_GL[0] + 1.0), 0.5 * _GL[1]


# -- label maps -----------------------------------------------------------------
class ShiftMap:
    """T(xi) = xi + s (measure preserving)."""

    def __init__(self, s):
        self.s = np.asarray(s, float)

    def __call__(self, xi):
        return xi + self.s

    def invert(self, q):
        return q - self.s

    def det(self, xi):
        return np.ones(xi.shape[0])


class SeparableMap:
    """T(xi) = (T_1(xi_1), T_2(xi_2)) with analytic inverse
        T_a^{-1}(q) = C_a(q) = q + sum_m [al_m sin(2 pi m q) - be_m (cos(2 pi m q) - 1)] / (2 pi m),
    i.e. node density per unit label rho_a = C_a' = 1 + sum_m (al_m cos + be_m sin) > 0.
    Smooth by construction (band limit M), exactly invertible."""

    def __init__(self, alpha, beta):
        self.al, self.be = np.asarray(alpha, float), np.asarray(beta, float)  # (2, M)
        self.m = np.arange(1, self.al.shape[1] + 1)

    def _C(self, q, a):
        ph = 2 * np.pi * q[:, None] * self.m
        C = q + np.sum((self.al[a] * np.sin(ph) - self.be[a] * (np.cos(ph) - 1.0)) / (2 * np.pi * self.m), axis=1)
        rho = 1.0 + np.sum(self.al[a] * np.cos(ph) + self.be[a] * np.sin(ph), axis=1)
        return C, rho

    def invert(self, q):
        return np.stack([self._C(q[:, a], a)[0] for a in (0, 1)], axis=1)

    def __call__(self, xi):
        out = []
        for a in (0, 1):
            q = xi[:, a].copy()
            for _ in range(60):
                C, rho = self._C(q, a)
                step = (C - xi[:, a]) / rho
                q -= np.clip(step, -0.2, 0.2)
                if np.max(np.abs(step)) < 1e-15:
                    break
            out.append(q)
        return np.stack(out, axis=1)

    def det(self, xi):
        q = self(xi)
        return 1.0 / (self._C(q[:, 0], 0)[1] * self._C(q[:, 1], 1)[1])

    def density(self, q, a):
        return self._C(q, a)[1]

    @classmethod
    def from_density(cls, rho1, rho2, modes=6, floor=0.35):
        """rho_a: positive samples on a uniform label grid (any length)."""
        co = []
        for rho in (rho1, rho2):
            r = np.maximum(rho / rho.mean(), floor)
            r /= r.mean()
            f = np.fft.rfft(r) / len(r)
            al, be = 2 * f[1:modes + 1].real, -2 * f[1:modes + 1].imag
            # guarantee positivity of the band-limited density
            qq = np.arange(2048) / 2048
            ph = 2 * np.pi * qq[:, None] * np.arange(1, modes + 1)
            osc = np.sum(al * np.cos(ph) + be * np.sin(ph), axis=1)
            sc = min(1.0, (1.0 - floor) / max(-osc.min(), 1e-300))
            co.append((al * sc, be * sc))
        return cls(np.stack([co[0][0], co[1][0]]), np.stack([co[0][1], co[1][1]]))


class OTMap:
    """Monge-Ampere optimal-transport label map T = id + grad(phi) computed by
    Fable's (vendored, unmodified) Newton-Krylov solver."""

    def __init__(self, M_grid, n):
        if _VENDOR not in sys.path:
            sys.path.insert(0, _VENDOR)
        import sheet2d as F
        self.info = F.solve_ot(M_grid, n)
        self.tm = F.TorusMap.from_potential(self.info["phi"])

    def __call__(self, xi):
        return self.tm(xi)

    def invert(self, q):
        xi, err = self.tm.invert(q)
        self.last_invert_err = err
        return xi

    def det(self, xi):
        return np.linalg.det(self.tm(xi, 1)[1])


def cell_masses(sh):
    """Label-space areas of the reference cells (uniform mu, total mass 1)."""
    n = sh.n
    q = sh.labels()
    d = q - grid(n)
    ge = np.stack(np.meshgrid(np.arange(n + 1) / n, np.arange(n + 1) / n, indexing="ij"), axis=-1)
    idx = np.arange(n + 1) % n
    Q = ge + d[idx][:, idx]
    a, b, c, e = Q[:-1, :-1], Q[1:, :-1], Q[:-1, 1:], Q[1:, 1:]
    cr = lambda u, v: u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]
    return 0.5 * (cr(b - a, e - a) + cr(e - a, c - a))


# -- monitors (functions of the nodal state only) ---------------------------------
def edge_action_defects(sh):
    """|eps_e| on horizontal / vertical edges, eps_e = int_e p_h.dx_h - dSigma,
    with the periodic cubic spline as (p_h, x_h).  (2, n, n).  This is the
    action-aware indicator: it measures the omega-pairing of the interpolation
    error with the edge tangent and needs no derivative estimates."""
    n = sh.n
    G = grid(n)
    sp = PSpline2(np.concatenate([sh.x - G, sh.p], axis=-1), 3)
    Sper = sh.S - G @ sh.periods
    out = np.zeros((2, n, n))
    for ax in (0, 1):
        E = np.zeros(n * n)
        for gx, gw in zip(_GLX, _GLW):
            xi = G.reshape(-1, 2).copy()
            xi[:, ax] += gx / n
            v, g = sp(xi, 1)
            dx = g[:, :2, ax] / n
            dx[:, ax] += 1.0 / n
            E += gw * np.sum(v[:, 2:] * dx, axis=1)
        out[ax] = np.abs(E.reshape(n, n) - (np.roll(Sper, -1, ax) - Sper + sh.periods[ax] / n))
    return out


def d4_indicator(sh, lp=1.0):
    """Conventional indicator: norm of the fourth undivided difference of
    (x, p/lp) along each reference direction.  (2, n, n)."""
    z = np.concatenate([sh.x - grid(sh.n), sh.p / lp], axis=-1)
    out = []
    for ax in (0, 1):
        r = lambda k: np.roll(z, -k, ax)
        d4 = r(2) - 4 * r(1) + 6 * z - 4 * r(-1) + r(-2)
        out.append(np.linalg.norm(d4, axis=-1))
    return np.stack(out)


def separable_map_from_indicator(sh, E, order_m, modes=6, floor=0.35):
    """Cheap graded map: per-direction equidistribution of the marginal
    indicator.  E (2, n, n) measured on the CURRENT mesh, local error ~ C h^m
    => new node density  rho_new ~ E^(1/m) rho_cur."""
    n = sh.n
    q = sh.labels()
    cur = sh.label_map if isinstance(sh.label_map, SeparableMap) else None
    rhos = []
    for a in (0, 1):
        e = np.sqrt(np.mean(E[a] ** 2, axis=1 - a))  # marginal along direction a
        qa = (q[:, 0, 0] if a == 0 else q[0, :, 1]) % 1.0
        rc = cur.density(qa, a) if cur is not None else np.ones(n)
        r = np.maximum(e, 1e-300) ** (1.0 / order_m) * rc
        o = np.argsort(qa)
        qq = (np.arange(512) + 0.5) / 512
        rhos.append(np.interp(qq, qa[o], r[o], period=1.0))
    return SeparableMap.from_density(rhos[0], rhos[1], modes=modes, floor=floor)


def ot_map_from_indicator(sh, E, order_m, width=0.05, floor=0.35, mgrid=64):
    """Isotropic OT map: density ~ |E|^(2/m) rho_cur (2D: h ~ rho^(-1/2)),
    mollified at a fixed label scale (cf. Fable Rem. 3.6)."""
    n = sh.n
    e = np.sqrt(E[0] ** 2 + E[1] ** 2)
    q = sh.labels().reshape(-1, 2)
    xi = grid(n).reshape(-1, 2)
    rc = 1.0 / sh.label_map.det(xi) if sh.label_map is not None else np.ones(n * n)
    r = np.maximum(e.ravel(), 1e-300) ** (2.0 / order_m) * rc
    # scatter to a uniform label grid (nearest cell average), then mollify
    idx = (np.floor(q[:, 0] % 1.0 * mgrid).astype(int) % mgrid) * mgrid + (np.floor(q[:, 1] % 1.0 * mgrid).astype(int) % mgrid)
    num = np.bincount(idx, weights=r, minlength=mgrid * mgrid)
    den = np.bincount(idx, minlength=mgrid * mgrid).astype(float)
    k = 2 * np.pi * np.fft.fftfreq(mgrid, d=1.0 / mgrid)
    k1, k2 = np.meshgrid(k, k, indexing="ij")
    filt = np.exp(-0.5 * width ** 2 * (k1 ** 2 + k2 ** 2))
    sm = lambda f: np.real(np.fft.ifft2(np.fft.fft2(f.reshape(mgrid, mgrid)) * filt))
    Mq = sm(num) / np.maximum(sm(den), 1e-12)
    Mq = np.maximum(Mq / Mq.mean(), floor)
    return OTMap(Mq / Mq.mean(), n)


# -- transfer -------------------------------------------------------------------
def remap(sh, carrier, new_map):
    """New sheet on the same reference grid with labels q = new_map(xi), nodal
    values read from `carrier` (name) built from the old nodal state."""
    n = sh.n
    car, _ = make_carrier(carrier, sh) if isinstance(carrier, str) else (carrier, 0.0)
    xi = grid(n).reshape(-1, 2)
    q_new = xi if new_map is None else new_map(xi)
    xi_old = q_new if sh.label_map is None else sh.label_map.invert(q_new)
    x, p, S = car(xi_old)
    # the carrier returns x(xi_old) = xi_old + periodic = xi_new + periodic, because
    # xi_old - xi_new is periodic: the degree-one convention of Sheet is preserved.
    return Sheet(x.reshape(n, n, 2), p.reshape(n, n, 2), S.reshape(n, n), sh.periods, new_map), xi_old


def refine2(sh, carrier):
    """Uniform h-refinement n -> 2n in the reference domain; the label map is
    unchanged.  Old nodes are retained exactly by interpolatory carriers."""
    n = sh.n
    car, _ = make_carrier(carrier, sh) if isinstance(carrier, str) else (carrier, 0.0)
    xi = grid(2 * n).reshape(-1, 2)
    x, p, S = car(xi)
    return Sheet(x.reshape(2 * n, 2 * n, 2), p.reshape(2 * n, 2 * n, 2), S.reshape(2 * n, 2 * n), sh.periods, sh.label_map)
