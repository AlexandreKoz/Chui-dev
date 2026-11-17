"""2D periodic particle-mesh gravity on the unit torus.

Units: box L = 1, total mass 1 (mean density 1), 4 pi G = 1, so that
    lap(phi) = rho - 1 ,
and a cold uniform medium is Jeans-unstable with growth rate exp(t).

Two force laws are provided, because the choice matters for a phase sheet:

 * "grad"  (common-potential / energy-conserving):  V(x) = sum_g phi_g W(x - x_g)
           is a single C^(r-2) function of position (B-spline kernel of order
           r) and the acceleration is its EXACT gradient, a(x) = -grad V(x).
           Every sheet point feels the same Hamiltonian H = |p|^2/2 + V(x, t):
           KDK is an exact symplectic map of R^4, action transport is exact,
           and all Poincare invariants of the sheet are conserved to round-off.
           Price: a small self-force, total momentum conserved only to
           truncation.
 * "mc"    (momentum-conserving, the standard PM choice): a grid acceleration
           a_g = -D_h phi (centred differences) interpolated with the same
           kernel W used for deposition (adjoint pair): zero self-force and
           exact momentum conservation, but a(x) is NOT a gradient, so the
           sheet's Lagrangian property / loop actions drift.
"""
from __future__ import annotations

import numpy as np


def _w1d(u, order):
    """B-spline assignment weights of the given order (2 CIC, 3 TSC, 4 PCS)
    at grid-unit positions u.  Returns base index i0 (M,), weights w (order, M)
    and derivative weights dw (order, M) w.r.t. u."""
    if order == 2:
        i0 = np.floor(u).astype(np.int64)
        t = u - i0
        w = np.stack([1 - t, t])
        dw = np.stack([-np.ones_like(t), np.ones_like(t)])
    elif order == 3:
        i0 = np.floor(u + 0.5).astype(np.int64) - 1
        t = u - (i0 + 1)  # in [-1/2, 1/2)
        w = np.stack([0.5 * (0.5 - t) ** 2, 0.75 - t * t, 0.5 * (0.5 + t) ** 2])
        dw = np.stack([-(0.5 - t), -2 * t, 0.5 + t])
    elif order == 4:
        i0 = np.floor(u).astype(np.int64) - 1
        t = u - (i0 + 1)  # in [0, 1)
        w = np.stack([(1 - t) ** 3 / 6, (3 * t ** 3 - 6 * t ** 2 + 4) / 6, (-3 * t ** 3 + 3 * t ** 2 + 3 * t + 1) / 6, t ** 3 / 6])
        dw = np.stack([-(1 - t) ** 2 / 2, (9 * t ** 2 - 12 * t) / 6, (-9 * t ** 2 + 6 * t + 3) / 6, t ** 2 / 2])
    else:
        raise ValueError(order)
    return i0, w, dw


class PM:
    def __init__(self, ng, order=4, force="grad", fd_order=4):
        self.ng, self.order, self.force, self.fd_order = ng, order, force, fd_order
        k = 2 * np.pi * np.fft.fftfreq(ng, d=1.0 / ng)
        k1, k2 = np.meshgrid(k, k, indexing="ij")
        k2s = k1 ** 2 + k2 ** 2
        k2s[0, 0] = 1.0
        self.green = -1.0 / k2s
        self.green[0, 0] = 0.0

    # -- deposition -------------------------------------------------------------
    def deposit(self, x, m):
        """Mass density (mean 1 for total mass 1) from points x (M,2), masses m."""
        ng = self.ng
        i0, wa, _ = _w1d(x[:, 0] * ng, self.order)
        j0, wb, _ = _w1d(x[:, 1] * ng, self.order)
        rho = np.zeros(ng * ng)
        m = np.broadcast_to(m, (x.shape[0],))
        for a in range(self.order):
            ia = (i0 + a) % ng
            for b in range(self.order):
                rho += np.bincount(ia * ng + (j0 + b) % ng, weights=m * wa[a] * wb[b], minlength=ng * ng)
        return rho.reshape(ng, ng) * ng * ng

    def potential(self, rho):
        return np.real(np.fft.ifft2(np.fft.fft2(rho - rho.mean()) * self.green))

    # -- gather -----------------------------------------------------------------
    def _gather(self, f, x, grad=False):
        ng = self.ng
        i0, wa, dwa = _w1d(x[:, 0] * ng, self.order)
        j0, wb, dwb = _w1d(x[:, 1] * ng, self.order)
        val = np.zeros(x.shape[0])
        g = np.zeros((x.shape[0], 2)) if grad else None
        for a in range(self.order):
            ia = (i0 + a) % ng
            for b in range(self.order):
                fg = f[ia, (j0 + b) % ng]
                val += wa[a] * wb[b] * fg
                if grad:
                    g[:, 0] += dwa[a] * wb[b] * fg * ng
                    g[:, 1] += wa[a] * dwb[b] * fg * ng
        return (val, g) if grad else val

    def field(self, phi, x):
        """Returns V(x) and acceleration a(x) at points x for the chosen law."""
        if self.force == "grad":
            V, g = self._gather(phi, x, grad=True)
            return V, -g
        ng = self.ng
        acc = []
        for ax in (0, 1):
            r = lambda s: np.roll(phi, -s, ax) - np.roll(phi, s, ax)
            d = r(1) * (0.5 * ng) if self.fd_order == 2 else (8 * r(1) - r(2)) * (ng / 12.0)
            acc.append(-self._gather(d, x))
        return self._gather(phi, x), np.stack(acc, axis=1)

    def energy_potential(self, rho, phi):
        """U = 1/2 int (rho - 1) phi  (grid quadrature)."""
        return 0.5 * float(np.mean((rho - rho.mean()) * phi))
