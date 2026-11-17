"""Sheet state, external Hamiltonians, exact-symplectic KDK with exact action
transport, and the offline oracle (never used by any online algorithm).

Sheet F(xi) = (x(xi), p(xi)) on the reference torus [0,1)^2, with
  x(xi + e_a) = x(xi) + e_a            (degree-one position map, box L = 1)
  p periodic
  S(xi + e_a) = S(xi) + c_a            (c = Liouville periods, conserved)
Labels q = T(xi) (label map; identity unless r-adapted).  Mass measure is
mu = Lebesgue(dq) with total mass 1; the mass of a reference cell is the
label-space area of its image under T.
"""
from __future__ import annotations

import numpy as np

from .pspline import grid

TWO_PI = 2.0 * np.pi


class TwoModePotential:
    """V(x) = -A [cos k x1 + cos k x2 + g cos k (x1 - x2)], k = 2 pi.
    (Same family as Fable's Potential2; A = 1/k^2, g = 0.5 by default.)"""

    def __init__(self, A=1.0 / TWO_PI ** 2, gamma=0.5):
        self.A, self.gamma = A, gamma

    def V(self, x, t=0.0):
        k = TWO_PI
        return -self.A * (np.cos(k * x[..., 0]) + np.cos(k * x[..., 1])
                          + self.gamma * np.cos(k * (x[..., 0] - x[..., 1])))

    def dV(self, x, t=0.0):
        k = TWO_PI
        s12 = self.gamma * np.sin(k * (x[..., 0] - x[..., 1]))
        return self.A * k * np.stack([np.sin(k * x[..., 0]) + s12, np.sin(k * x[..., 1]) - s12], axis=-1)


class DrivenPotential(TwoModePotential):
    """Time-dependent, non-separable:  adds eps cos(k (x1 + 2 x2) - Om t)."""

    def __init__(self, A=1.0 / TWO_PI ** 2, gamma=0.5, eps=0.3, Om=2.0):
        super().__init__(A, gamma)
        self.eps, self.Om = eps, Om

    def V(self, x, t=0.0):
        k = TWO_PI
        return super().V(x) - self.A * self.eps * np.cos(k * (x[..., 0] + 2 * x[..., 1]) - self.Om * t)

    def dV(self, x, t=0.0):
        k = TWO_PI
        s = self.A * self.eps * k * np.sin(k * (x[..., 0] + 2 * x[..., 1]) - self.Om * t)
        return super().dV(x) + np.stack([s, 2 * s], axis=-1)


def kdk(pot, x, p, S, dt, nsteps, t0=0.0):
    """Kick-drift-kick.  Each sub-map is an exact symplectic map of R^4 with
    an explicit generating function; S is advanced by exactly that function:
        kick  : p -= k dV(x),          S -= k V(x)
        drift : x += d p,              S += p.(x_new - x) - d |p|^2 / 2
    (= +d|p|^2/2; the first form uses the displacement actually applied in
    floating point).  Hence dS = p.dx on the transported sheet for any dt."""
    t = t0
    for _ in range(nsteps):
        S = S - 0.5 * dt * pot.V(x, t)
        p = p - 0.5 * dt * pot.dV(x, t)
        xn = x + dt * p
        S = S + np.sum(p * (xn - x), axis=-1) - 0.5 * dt * np.sum(p * p, axis=-1)
        x = xn
        t += dt
        S = S - 0.5 * dt * pot.V(x, t)
        p = p - 0.5 * dt * pot.dV(x, t)
    return x, p, S


class ColdIC:
    """Initial sheet  x0 = q + a grad(chi)(q),  p0 = u0 + b grad(chi)(q),
    chi(q) = (cos 2 pi q1 + cos 2 pi (q2 + 0.3) ) / (2 pi)^2 scaled by amp.
    Exactly Lagrangian (dp0^dx0 = 0 because both are gradients of the same
    chi up to constants) with action
        S0 = u0 . x0 + b chi + a b |grad chi|^2 / 2,
    Liouville periods c = u0.  amp = 0, u0 = 0 is the cold sheet at rest."""

    def __init__(self, amp=0.0, a=0.0, b=0.0, u0=(0.0, 0.0)):
        self.amp, self.a, self.b, self.u0 = amp, a, b, np.asarray(u0, float)

    def _chi(self, q):
        k = TWO_PI
        c = self.amp / k ** 2
        chi = c * (np.cos(k * q[..., 0]) + np.cos(k * (q[..., 1] + 0.3)))
        g = -c * k * np.stack([np.sin(k * q[..., 0]), np.sin(k * (q[..., 1] + 0.3))], axis=-1)
        return chi, g

    def __call__(self, q):
        chi, g = self._chi(q)
        x = q + self.a * g
        p = self.u0 + self.b * g
        S = x @ self.u0 + self.b * chi + 0.5 * self.a * self.b * np.sum(g * g, axis=-1)
        return x, p, S

    @property
    def periods(self):
        return self.u0.copy()


class Oracle:
    """Exact samples of the discretely transported sheet at arbitrary labels.
    OFFLINE EVALUATION ONLY: never called by reconstruction, monitors, remaps
    or forces."""

    def __init__(self, pot, ic, dt):
        self.pot, self.ic, self.dt = pot, ic, dt

    def __call__(self, q, nsteps):
        x, p, S = self.ic(q)
        return kdk(self.pot, x, p, S, self.dt, nsteps)

    def tangent(self, q, nsteps, eps=1e-6):
        """Columns d_a (x, p) / d q_a by central differences of the oracle."""
        out = []
        for a in range(2):
            e = np.zeros(2)
            e[a] = eps
            xp, pp, _ = self(q + e, nsteps)
            xm, pm, _ = self(q - e, nsteps)
            out.append(np.concatenate([xp - xm, pp - pm], axis=-1) / (2 * eps))
        return np.stack(out, axis=-1)  # (..., 4, 2)


class Sheet:
    """Nodal state on the uniform n x n reference grid."""

    def __init__(self, x, p, S, periods=(0.0, 0.0), label_map=None):
        self.x, self.p, self.S = x, p, S
        self.n = x.shape[0]
        self.periods = np.asarray(periods, float)
        self.label_map = label_map  # None = identity; else object with __call__(xi)->q and invert(q)->xi

    @classmethod
    def from_ic(cls, ic, n, label_map=None):
        xi = grid(n)
        q = xi if label_map is None else label_map(xi.reshape(-1, 2)).reshape(n, n, 2)
        x, p, S = ic(q)
        return cls(x, p, S, ic.periods, label_map)

    def labels(self):
        xi = grid(self.n)
        return xi if self.label_map is None else self.label_map(xi.reshape(-1, 2)).reshape(self.n, self.n, 2)

    def advance(self, pot, dt, nsteps, t0=0.0):
        self.x, self.p, self.S = kdk(pot, self.x, self.p, self.S, dt, nsteps, t0)
        return self
