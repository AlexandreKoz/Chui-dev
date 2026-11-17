"""d = 3 local symplectic-frame construction (3D labels, 6D phase space).

CHART-WISE (non-conforming) version of the atlas carrier:
  z_h  : tricubic Hermite interpolant of nodal (x, p), nodal derivatives by
         4th-order centred differences;
  W_h  : tricubic Hermite interpolant of the nodal action Sigma with nodal
         gradient beta_a = p . d_a x (frame independent) and symmetrised
         difference estimates of the mixed derivatives;
  kappa = grad W_h - p_h . D x_h   (frame-independent action-mismatch one-form);
  cell frame (Ur, Ui) in U(3):  x = x_h - Ui dP,  p = p_h + Ur dP,
         dP = (D Y)^{-T} kappa,   D Y = Ur^T D x_h + Ui^T D p_h.
Inside each cell P.dY = d(W_h + Q_M o z_h) exactly (cell-wise exactly
Lagrangian).  With Ur = I, Ui = 0 this is Fable's configuration-chart front,
p = (Dx)^{-T} grad W_h, which is singular at caustics.  kappa does NOT vanish
on cell faces here (no edge/face action closure in 3D yet), so neighbouring
cells with different frames differ by a measurable seam jump.
"""
from __future__ import annotations

import itertools

import numpy as np

from .conventions import unitary_frame
from .pspline import hermite

TWO_PI = 2.0 * np.pi


class Potential3:
    """V = -A [sum_a w_a cos(k x_a) + g cos(k (x1 - x2 + x3))]."""

    def __init__(self, w=(1.0, 0.6, 0.35), gamma=0.4, A=1.0 / TWO_PI ** 2):
        self.w, self.g, self.A = np.asarray(w), gamma, A

    def V(self, x, t=0.0):
        k = TWO_PI
        return -self.A * (np.sum(self.w * np.cos(k * x), axis=-1) + self.g * np.cos(k * (x[..., 0] - x[..., 1] + x[..., 2])))

    def dV(self, x, t=0.0):
        k = TWO_PI
        s = self.g * np.sin(k * (x[..., 0] - x[..., 1] + x[..., 2]))
        return self.A * k * (self.w * np.sin(k * x) + s[..., None] * np.array([1.0, -1.0, 1.0]))


def grid3(shape):
    return np.stack(np.meshgrid(*[np.arange(m) / m for m in shape], indexing="ij"), axis=-1)


def _fd(f, ax, m):
    r = lambda k: np.roll(f, -k, ax) - np.roll(f, k, ax)
    return (8 * r(1) - r(2)) * (m / 12.0)


def jets(f, shape, first=None):
    """Nodal derivative table {alpha in {0,1}^3: array} w.r.t. xi; `first`
    optionally supplies exact first derivatives (list of 3 arrays)."""
    J = {(0, 0, 0): f}
    g = [(_fd(f, a, shape[a]) if first is None else first[a]) for a in range(3)]
    for a in range(3):
        J[tuple(int(a == b) for b in range(3))] = g[a]
    for a, b in ((0, 1), (0, 2), (1, 2)):
        al = tuple(int(c in (a, b)) for c in range(3))
        J[al] = 0.5 * (_fd(g[a], b, shape[b]) + _fd(g[b], a, shape[a]))
    J[(1, 1, 1)] = (_fd(_fd(g[0], 1, shape[1]), 2, shape[2]) + _fd(_fd(g[1], 0, shape[0]), 2, shape[2])
                    + _fd(_fd(g[2], 0, shape[0]), 1, shape[1])) / 3.0
    return J


def herm3(J, shape, I, s):
    """Tricubic Hermite value and LOCAL gradient (d/ds) at cells I, coords s."""
    B = [hermite(s[:, a]) for a in range(3)]
    val, grad = 0.0, [0.0, 0.0, 0.0]
    for c in itertools.product((0, 1), repeat=3):
        idx = tuple((I[:, a] + c[a]) % shape[a] for a in range(3))
        for al in itertools.product((0, 1), repeat=3):
            f = J[al][idx]
            sc = np.prod([(1.0 / shape[a]) ** al[a] for a in range(3)])
            b = [B[a][0][2 * al[a] + c[a]] for a in range(3)]
            db = [B[a][1][2 * al[a] + c[a]] for a in range(3)]
            ex = (lambda w: w[:, None]) if f.ndim > 1 else (lambda w: w)
            val = val + sc * ex(b[0] * b[1] * b[2]) * f
            grad[0] = grad[0] + sc * ex(db[0] * b[1] * b[2]) * f
            grad[1] = grad[1] + sc * ex(b[0] * db[1] * b[2]) * f
            grad[2] = grad[2] + sc * ex(b[0] * b[1] * db[2]) * f
    return val, np.stack(grad, axis=-1)


class Carrier3:
    def __init__(self, x, p, S, frame="osc", lp=1.0):
        self.shape = x.shape[:3]
        self.lp, self.frame = lp, frame
        z = np.concatenate([x - grid3(self.shape), p / lp], axis=-1)
        self.Jz = jets(z, self.shape)
        dx = [self.Jz[tuple(int(a == b) for b in range(3))][..., :3].copy() for a in range(3)]
        for a in range(3):
            dx[a][..., a] += 1.0
        beta = [np.sum(z[..., 3:] * dx[a], axis=-1) for a in range(3)]
        self.JW = jets(S / lp, self.shape, first=beta)

    def zh(self, I, s):
        v, g = herm3(self.Jz, self.shape, I, s)
        z = v.copy()
        dz = g.copy()  # (M, 6, 3) local derivatives
        for a in range(3):
            z[:, a] += (I[:, a] + s[:, a]) / self.shape[a]
            dz[:, a, a] += 1.0 / self.shape[a]
        return z, dz

    def frames(self, I):
        if self.frame == "x":
            return np.broadcast_to(np.eye(3), (len(I), 3, 3)), np.zeros((len(I), 3, 3))
        _, dz = self.zh(I, np.full((len(I), 3), 0.5))
        sc = np.array(self.shape, float)
        Tx, Tp = dz[:, :3, :] * sc, dz[:, 3:, :] * sc
        if self.frame == "osc":
            return unitary_frame(Tx, Tp)
        best = None  # Arnold: best of the 8 coordinate charts
        for mask in range(8):
            u = np.array([1j if (mask >> k) & 1 else 1.0 for k in range(3)])
            r, i_ = np.diag(u.real), np.diag(u.imag)
            smin = np.linalg.svd(r.T @ Tx + i_.T @ Tp, compute_uv=False)[:, -1]
            if best is None:
                best, Ur, Ui = smin.copy(), np.tile(r, (len(I), 1, 1)), np.tile(i_, (len(I), 1, 1))
            else:
                sel = smin > best
                best[sel] = smin[sel]
                Ur[sel], Ui[sel] = r, i_
        return Ur, Ui

    def eval_local(self, I, s, lagrangian=True, diag=False):
        z, dz = self.zh(I, s)
        if not lagrangian:
            return z[:, :3], z[:, 3:] * self.lp
        _, dW = herm3(self.JW, self.shape, I, s)
        kappa = dW - np.einsum("mc,mca->ma", z[:, 3:], dz[:, :3, :])
        Ur, Ui = self.frames(I)
        DY = np.einsum("mck,mca->mka", Ur, dz[:, :3, :]) + np.einsum("mck,mca->mka", Ui, dz[:, 3:, :])
        dP = np.linalg.solve(np.swapaxes(DY, 1, 2), kappa[..., None])[..., 0]
        x = z[:, :3] - np.einsum("mck,mk->mc", Ui, dP)
        p = (z[:, 3:] + np.einsum("mck,mk->mc", Ur, dP)) * self.lp
        if diag:
            sc = np.array(self.shape, float)
            sv = np.linalg.svd(DY * sc, compute_uv=False)
            svx = np.linalg.svd(dz[:, :3, :] * sc, compute_uv=False)
            return x, p, dict(cond=sv[:, 0] / sv[:, -1], sig_min=sv[:, -1], sig_min_x=svx[:, -1], detDx=np.linalg.det(dz[:, :3, :] * sc))
        return x, p

    def __call__(self, xi, **kw):
        u = np.asarray(xi) * np.array(self.shape)
        I = np.floor(u).astype(np.int64)
        return self.eval_local(I, u - I, **kw)
