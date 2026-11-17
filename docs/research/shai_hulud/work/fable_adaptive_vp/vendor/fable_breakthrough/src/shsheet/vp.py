"""Self-consistent 2D periodic Vlasov-Poisson evolution of a cold phase-space
sheet (static background, no expansion) with PM gravity and KDK.

State: Sheet (x, p, Sigma) on the n x n reference grid; mass measure uniform
in labels (total mass 1).  Units as in pm.py (4 pi G = 1, mean density 1).

Density assignment ("deposit"):
  "nodes"         every node is a particle of mass 1/n^2 (conventional N-body PM)
  "lin|cub|quin|ac3|csa3"  sheet quadrature: each reference cell is sampled on
                  an sub x sub midpoint lattice of labels, positions are read
                  from the named carrier of the CURRENT nodal state, and each
                  sample carries the pulled-back mass 1/(n sub)^2.  All streams
                  are summed automatically (the map need not be injective).
The potential is then a single function V(x, t) on the box and nodes (and any
massless tracers) are kicked with the PM force law selected in pm.PM.
Action transport: Sigma -= kick * V(x), Sigma += p.dx - drift |p|^2/2.
"""
from __future__ import annotations

import time

import numpy as np

from .carriers import make_carrier
from .pspline import grid
from .sheet import Sheet

TWO_PI = 2.0 * np.pi


class GrowingModeIC:
    """x0 = q + grad chi,  p0 = grad chi  (linear growing mode, rate exp(t)),
    chi = [e1 cos(k q1) + e2 cos(k q2) + e3 cos(k (q1 + q2) + ph)] / k^2.
    Exactly Lagrangian, Sigma0 = chi + |grad chi|^2 / 2, zero periods.
    Linear density contrast at t=0: delta = e1 cos + e2 cos + 2 e3 cos."""

    periods = np.zeros(2)

    def __init__(self, eps=(0.10, 0.07, 0.03), phase=0.9):
        self.eps, self.phase = eps, phase

    def __call__(self, q):
        k = TWO_PI
        e1, e2, e3 = self.eps
        a3 = k * (q[..., 0] + q[..., 1]) + self.phase
        chi = (e1 * np.cos(k * q[..., 0]) + e2 * np.cos(k * q[..., 1]) + e3 * np.cos(a3)) / k ** 2
        g = -np.stack([e1 * np.sin(k * q[..., 0]) + e3 * np.sin(a3), e2 * np.sin(k * q[..., 1]) + e3 * np.sin(a3)], axis=-1) / k
        return q + g, g.copy(), chi + 0.5 * np.sum(g * g, axis=-1)


def _bilinear_points(x, sub):
    """Positions at the sub x sub midpoint lattice of every cell (bilinear)."""
    n = x.shape[0]
    d = x - grid(n)
    t = (np.arange(sub) + 0.5) / sub
    out = np.empty((n, sub, n, sub, 2))
    d10, d01, d11 = np.roll(d, -1, 0), np.roll(d, -1, 1), np.roll(np.roll(d, -1, 0), -1, 1)
    for a, s in enumerate(t):
        for b, u in enumerate(t):
            out[:, a, :, b] = (1 - s) * (1 - u) * d + s * (1 - u) * d10 + (1 - s) * u * d01 + s * u * d11
    xi = grid(n * sub) + 0.5 / (n * sub)
    return out.reshape(n * sub, n * sub, 2) + xi


class VPRun:
    def __init__(self, sh, pm, deposit="nodes", sub=1, tracers=None):
        self.sh, self.pm, self.deposit, self.sub = sh, pm, deposit, sub
        self.t = 0.0
        self.tr = tracers  # dict(x, p, S) of massless tracers or None
        self.phi = None
        self.rho = None
        self.cpu = dict(deposit=0.0, poisson=0.0, gather=0.0, transfer=0.0)
        self.n_field_solves = 0

    # -- density ------------------------------------------------------------------
    def deposit_points(self):
        sh, n, s = self.sh, self.sh.n, self.sub
        if self.deposit == "nodes":
            return sh.x.reshape(-1, 2), 1.0 / (n * n)
        if self.deposit == "lin":
            return _bilinear_points(sh.x, s).reshape(-1, 2), 1.0 / (n * s) ** 2
        car, _ = make_carrier(self.deposit, sh)
        xi = (grid(n * s) + 0.5 / (n * s)).reshape(-1, 2)
        return car(xi)[0], 1.0 / (n * s) ** 2

    def solve_field(self):
        t0 = time.perf_counter()
        xs, m = self.deposit_points()
        self.rho = self.pm.deposit(xs, m)
        t1 = time.perf_counter()
        self.phi = self.pm.potential(self.rho)
        self.cpu["deposit"] += t1 - t0
        self.cpu["poisson"] += time.perf_counter() - t1
        self.n_field_solves += 1

    def _kick(self, k):
        t0 = time.perf_counter()
        sh = self.sh
        V, a = self.pm.field(self.phi, sh.x.reshape(-1, 2))
        sh.p = sh.p + k * a.reshape(sh.p.shape)
        sh.S = sh.S - k * V.reshape(sh.S.shape)
        if self.tr is not None:
            V, a = self.pm.field(self.phi, self.tr["x"])
            self.tr["p"] = self.tr["p"] + k * a
            self.tr["S"] = self.tr["S"] - k * V
        self.cpu["gather"] += time.perf_counter() - t0

    def _drift(self, d):
        for o in ((self.sh,) if self.tr is None else (self.sh, self.tr)):
            x, p, S = (o.x, o.p, o.S) if isinstance(o, Sheet) else (o["x"], o["p"], o["S"])
            xn = x + d * p
            Sn = S + np.sum(p * (xn - x), axis=-1) - 0.5 * d * np.sum(p * p, axis=-1)
            if isinstance(o, Sheet):
                o.x, o.S = xn, Sn
            else:
                o["x"], o["S"] = xn, Sn

    def step(self, dt):
        if self.phi is None:
            self.solve_field()
        self._kick(0.5 * dt)
        self._drift(dt)
        self.solve_field()
        self._kick(0.5 * dt)
        self.t += dt

    # -- diagnostics --------------------------------------------------------------
    def conserved(self):
        if self.phi is None:
            self.solve_field()
        n = self.sh.n
        K = 0.5 * float(np.sum(self.sh.p ** 2)) / (n * n)
        U = self.pm.energy_potential(self.rho, self.phi)
        P = self.sh.p.reshape(-1, 2).sum(axis=0) / (n * n)
        return dict(t=self.t, K=K, U=U, E=K + U, Px=float(P[0]), Py=float(P[1]), mass=float(self.rho.mean()))


def ring_tracers(ic, centre, radius, count):
    th = TWO_PI * np.arange(count) / count
    q = np.asarray(centre)[None, :] + radius * np.stack([np.cos(th), np.sin(th)], axis=1)
    x, p, S = ic(q)
    return dict(x=x, p=p, S=S)


def ring_circulation(tr):
    """oint p.dx along the closed tracer ring (spectral in the ring angle) and
    the max residual of dS/dtheta - p.dx/dtheta (action consistency)."""
    x, p, S = tr["x"], tr["p"], tr["S"]
    N = x.shape[0]
    ik = 1j * np.fft.fftfreq(N, d=1.0 / N)
    dx = np.real(np.fft.ifft(ik[:, None] * np.fft.fft(x, axis=0), axis=0))
    dS = np.real(np.fft.ifft(ik * np.fft.fft(S)))
    integrand = np.sum(p * dx, axis=1)
    spec = np.abs(np.fft.fft(integrand))
    tail = float(spec[N // 2 - N // 16: N // 2].max() / max(spec.max(), 1e-300))  # resolution indicator
    return float(TWO_PI * integrand.mean()), float(np.abs(dS - integrand).max()), tail
