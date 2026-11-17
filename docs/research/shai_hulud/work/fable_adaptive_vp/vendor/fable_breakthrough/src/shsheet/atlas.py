"""Conforming symplectic-frame atlas carrier (CSA carrier).

A continuous, node-interpolating, cell-wise EXACTLY Lagrangian reconstruction
of a 2D sheet in 4D phase space from nodal (x, p, Sigma) on a uniform
reference grid, valid through configuration-space caustics.

Construction (see research/atlas_derivation.md for the theorem and proof):

 1. z_h   : periodic spline (degree 3 or 5) of the nodal (x, p).
 2. nodal tangent fix : nodal first derivatives are moved (closed form) onto
    Lagrangian planes, omega(d_1 z, d_2 z) = 0 at every node, by a C^1
    bicubic-Hermite correction that changes no nodal value.
 3. edge-action bubbles : on every grid edge a quartic bubble along the
    omega-conjugate of the chord makes  int_edge p.dx = Sigma_j - Sigma_i
    exactly (the equation is exactly linear in the amplitude).  Bubbles are
    extended into the two adjacent cells with cubic Hermite blends, so that
    z_h' (steps 1-3) is C^1 and edge traces are single valued.
 4. generating potential : beta = z_h'^* theta = p_h'.dx_h' is a one-form on
    label space.  W_h is the bicubically blended Coons patch with boundary
    traces Sigma_i + int beta (single valued by 3) and boundary normal
    derivatives beta_n (twist-compatible by 2).  Then kappa := dW_h - beta
    vanishes identically on every cell boundary.
 5. frame correction : in each cell choose ANY linear symplectic frame
    (Y, P) = M (x, p) in which D_xi Y is invertible on the cell (default: the
    osculating unitary frame of the tangent plane at the cell centre), and
    set   Y = Y_h',  P = P_h' + (D_xi Y)^{-T} kappa.
    Then P.dY = d(W_h + Q_M o z_h') exactly in the cell, the carrier equals
    z_h' on all cell edges (so it is continuous for arbitrary, mutually
    unrelated frames in neighbouring cells), and
        p.dx = d Sigma_car,  Sigma_car = W_h + Q_M(z_h') - Q_M(z_car),
    with Sigma_car continuous across cells: a global, conforming, exactly
    Lagrangian C^0 immersion with the nodal values and periods of the data.

No reference/oracle information is used anywhere.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial import chebyshev as Ch

from .conventions import unitary_frame
from .pspline import PSpline2, bubble, grid, hermite

_GL = np.polynomial.legendre.leggauss(8)
_GLX, _GLW = 0.5 * (_GL[0] + 1.0), 0.5 * _GL[1]


def _K(v):
    """omega-conjugation: omega(a, b) = a . K(b);  K(b) = (-b_p, b_x)."""
    return np.concatenate([-v[..., 2:], v[..., :2]], axis=-1)


def _om(a, b):
    return np.sum(a * _K(b), axis=-1)


def momentum_scale(x, p):
    """Momentum unit lp balancing the scaled tangent blocks, |Dp| / lp ~ |Dx|
    (state-only; defines the metric (x, p / lp) used for 'normal' directions)."""
    n = x.shape[0]
    xp = x - grid(n)
    d = lambda f, ax: (np.roll(f, -1, ax) - np.roll(f, 1, ax)) * (0.5 * n)
    gx2 = 2.0 + np.mean(d(xp, 0) ** 2 + d(xp, 1) ** 2) * 2
    gp2 = np.mean(d(p, 0) ** 2 + d(p, 1) ** 2) * 2
    return float(np.sqrt(gp2 / gx2)) if gp2 > 1e-24 else 1.0


class AtlasCarrier:
    def __init__(self, x, p, S, order=3, frame="osc", periods=(0.0, 0.0), tangent_fix=True,
                 edge_action=True, lagrangian=True, lp=None, block=1, interior=4, rho_min=0.5, gauge="normal", x_thr=0.5):
        n = x.shape[0]
        self.n, self.h, self.order = n, 1.0 / n, order
        self.frame_kind, self.lagrangian = frame, lagrangian
        self.interior, self.amode = 0, None
        G = grid(n)
        zper = np.concatenate([x - G, p], axis=-1)
        if lp is None:
            lp = momentum_scale(x, p)
        self.lp = lp
        zper[..., 2:] /= lp
        self.sp = PSpline2(zper, order)
        self.c = np.asarray(periods, float) / lp
        self.Sp = S / lp - G @ (np.asarray(periods, float) / lp)  # periodic part of the scaled action
        self.dT = np.zeros((n, n, 4, 2))
        self.ah = np.zeros((n, n, 4))
        self.av = np.zeros((n, n, 4))
        ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        self._I = np.stack([ii.ravel(), jj.ravel()], axis=1)
        N = n * n
        zero = np.zeros((N, 2))

        # -- 2. nodal tangent fix ------------------------------------------------
        z0, dz0 = self._zh(self._I, zero)
        self.node_defect_before = float(np.abs(_om(dz0[:, :, 0], dz0[:, :, 1])).max()) * n * n
        if tangent_fix:
            t1, t2 = dz0[:, :, 0], dz0[:, :, 1]
            c = _om(t1, t2)
            g = np.sum(t1 * t1, axis=-1) + np.sum(t2 * t2, axis=-1)
            lam = 2.0 * c / (g + np.sqrt(np.maximum(g * g - 4.0 * c * c, 0.0)))
            self.dT[..., 0] = (-lam[:, None] * _K(t2)).reshape(n, n, 4)
            self.dT[..., 1] = (lam[:, None] * _K(t1)).reshape(n, n, 4)
            z0, dz0 = self._zh(self._I, zero)
        self.node_defect = float(np.abs(_om(dz0[:, :, 0], dz0[:, :, 1])).max()) * n * n

        # -- 3. edge-action bubbles ----------------------------------------------
        self.edge_defect = np.zeros((2, n, n))
        for ax in (0, 1):
            E = np.zeros(N)
            A1 = np.zeros(N)
            e = np.zeros(2, dtype=np.int64)
            e[ax] = 1
            z1, _ = self._zh(self._I + e, zero)
            ch = z1 - z0
            nvec = _K(ch) / np.linalg.norm(ch, axis=-1, keepdims=True)
            if gauge == "x":
                # position-preserving gauge: wherever the edge is not close to a
                # caustic direction, correct the action with a pure-momentum
                # bubble (x keeps its spline value); otherwise use the
                # omega-conjugate of the chord.
                cx = np.linalg.norm(ch[:, :2], axis=-1)
                use_p = cx >= x_thr * np.linalg.norm(ch, axis=-1)
                nvec[use_p, :2] = 0.0
                nvec[use_p, 2:] = ch[use_p, :2] / cx[use_p, None]
            for gx, gw in zip(_GLX, _GLW):
                s = zero.copy()
                s[:, ax] = gx
                z, dz = self._zh(self._I, s)
                E += gw * np.sum(z[:, 2:] * dz[:, :2, ax], axis=-1)
                A1 += gw * bubble(gx)[0] * _om(nvec, dz[:, :, ax])
            Spn = self.Sp.ravel()
            Sp1 = np.roll(self.Sp, -1, ax).ravel()
            eps = E - (Sp1 - Spn + self.c[ax] * self.h)
            self.edge_defect[ax] = eps.reshape(n, n)
            if edge_action:
                amp = (-eps / A1)[:, None] * nvec
                if ax == 0:
                    self.ah = amp.reshape(n, n, 4)
                else:
                    self.av = amp.reshape(n, n, 4)

        # -- 4. edge polynomial tables (Chebyshev) for beta = p.dx ---------------
        m = 9 if order == 3 else 12
        uu = np.cos(np.pi * np.arange(m) / (m - 1))
        Vinv = np.linalg.inv(Ch.chebvander(uu, m - 1))
        self.tab = {}
        for ax in (0, 1):
            bt = np.zeros((m, N))  # tangential component of beta along the edge
            bn = np.zeros((m, N))  # normal component
            for k_, u in enumerate(uu):
                s = zero.copy()
                s[:, ax] = 0.5 * (u + 1.0)
                z, dz = self._zh(self._I, s)
                bt[k_] = np.sum(z[:, 2:] * dz[:, :2, ax], axis=-1)
                bn[k_] = np.sum(z[:, 2:] * dz[:, :2, 1 - ax], axis=-1)
            ct = (Vinv @ bt).reshape(m, n, n)
            cn = (Vinv @ bn).reshape(m, n, n)
            self.tab[ax] = dict(t=ct, n=cn, I=Ch.chebint(ct, lbnd=-1, scl=0.5, axis=0), dn=Ch.chebder(cn, scl=2.0, axis=0))
        sgn = (-1.0) ** np.arange(m)
        self.Ws = np.tensordot(sgn, self.tab[0]["t"], axes=(0, 0))  # beta_s at nodes
        self.Wt = np.tensordot(sgn, self.tab[1]["t"], axes=(0, 0))  # beta_t at nodes
        sgd = (-1.0) ** np.arange(m - 1)
        tw_a = np.tensordot(sgd, self.tab[0]["dn"], axes=(0, 0))  # d_s beta_t at nodes
        tw_b = np.tensordot(sgd, self.tab[1]["dn"], axes=(0, 0))  # d_t beta_s at nodes
        self.twist_mismatch = float(np.abs(tw_a - tw_b).max())
        self.Wst = 0.5 * (tw_a + tw_b)
        # edge-action closure after the bubbles (should be round-off)
        cl = []
        for ax in (0, 1):
            Iend = self.tab[ax]["I"].sum(axis=0)  # value at u = +1
            cl.append(np.abs(Iend - (np.roll(self.Sp, -1, ax) - self.Sp + self.c[ax] * self.h)).max())
        self.edge_closure = float(max(cl))

        # -- 5. frames ------------------------------------------------------------
        nb = n // block
        cI = (self._I.reshape(n, n, 2)[::block, ::block] + (block // 2)).reshape(-1, 2)
        _, dzc = self._zh(cI, np.full((cI.shape[0], 2), 0.5 * (block % 2)))
        Tx, Tp = dzc[:, :2, :], dzc[:, 2:, :]
        if frame == "x":
            Ur = np.broadcast_to(np.eye(2), (cI.shape[0], 2, 2)).copy()
            Ui = np.zeros_like(Ur)
        elif frame in ("osc", "auto"):
            Ur, Ui = unitary_frame(Tx, Tp)
        elif frame == "arnold":
            best = None
            for mask in range(4):
                u = np.array([1j if (mask >> k_) & 1 else 1.0 for k_ in range(2)])
                r, i_ = np.diag(u.real), np.diag(u.imag)
                smin = np.linalg.svd(r.T @ Tx + i_.T @ Tp, compute_uv=False)[:, -1]
                if best is None:
                    best, Ur, Ui = smin.copy(), np.tile(r, (len(smin), 1, 1)), np.tile(i_, (len(smin), 1, 1))
                else:
                    sel = smin > best
                    best[sel] = smin[sel]
                    Ur[sel], Ui[sel] = r, i_
        else:
            raise ValueError(frame)
        rep = lambda a: np.repeat(np.repeat(a.reshape(nb, nb, 2, 2), block, 0), block, 1)
        self.Ur, self.Ui = rep(Ur), rep(Ui)
        self.x_chart = np.zeros((n, n), dtype=bool)
        if frame == "auto":
            # configuration chart (p-only correction, x untouched) wherever
            # sigma_min(Dx) / sigma_min(Dz) >= x_thr on the probe lattice,
            # osculating frame elsewhere.  Mixed neighbours are conforming.
            rx = np.full(N, np.inf)
            for a_ in (0.02, 0.5, 0.98):
                for b_ in (0.02, 0.5, 0.98):
                    _, dz = self._zh(self._I, np.tile([a_, b_], (N, 1)))
                    sx = np.linalg.svd(dz[:, :2, :], compute_uv=False)[:, -1] * np.sign(np.linalg.det(dz[:, :2, :]))
                    rx = np.minimum(rx, sx / np.linalg.svd(dz, compute_uv=False)[:, -1])
            self.x_chart = (rx >= x_thr).reshape(n, n)
            self.Ur[self.x_chart] = np.eye(2)
            self.Ui[self.x_chart] = 0.0

        # -- 5b. chart validity: rho = sigma_min(D Y) / sigma_min(D z) is the cosine
        # of the largest principal angle between the local tangent plane and the
        # frame plane.  Cells with rho < rho_min (or det DY <= 0) anywhere on a
        # 3x3 probe lattice FALL BACK to z_h' (not exactly Lagrangian inside the
        # cell, still continuous, still zero net omega-flux through the cell).
        rho = np.full(N, np.inf)
        for a_ in (0.02, 0.5, 0.98):
            for b_ in (0.02, 0.5, 0.98):
                _, dz = self._zh(self._I, np.tile([a_, b_], (N, 1)))
                DY = (np.einsum("mck,mca->mka", self.Ur.reshape(N, 2, 2), dz[:, :2, :])
                      + np.einsum("mck,mca->mka", self.Ui.reshape(N, 2, 2), dz[:, 2:, :]))
                sv = np.linalg.svd(DY, compute_uv=False)[:, -1] * np.sign(np.linalg.det(DY))
                rho = np.minimum(rho, sv / np.linalg.svd(dz, compute_uv=False)[:, -1])
        self.rho = rho.reshape(n, n)
        self.valid = self.rho >= rho_min
        self.fallback_fraction = float(1.0 - self.valid.mean())

        # -- 4b. optional interior modes: W_h += sum_k a_k phi_k with phi_k and
        # grad phi_k vanishing on the cell boundary, a = argmin of the
        # G^{-1}-weighted L2 norm of kappa (size of the normal correction).
        if not lagrangian:
            self.valid[:] = False
            self.fallback_fraction = 1.0
        if interior and lagrangian:
            K = interior
            gl = np.polynomial.legendre.leggauss(6)
            gx, gw = 0.5 * (gl[0] + 1.0), 0.5 * gl[1]
            A = np.zeros((N, K, K))
            r = np.zeros((N, K))
            for a_, wa in zip(gx, gw):
                for b_, wb in zip(gx, gw):
                    s = np.tile([a_, b_], (N, 1))
                    z, dz = self._zh(self._I, s)
                    _, dW = self._potential(self._I, s)
                    kap = dW - np.einsum("mc,mca->ma", z[:, 2:], dz[:, :2, :])
                    Gi = np.linalg.inv(np.einsum("mca,mcb->mab", dz, dz))
                    gphi = self._modes(s, K)[1]  # (M, K, 2)
                    Gg = np.einsum("mab,mkb->mka", Gi, gphi)
                    A += wa * wb * np.einsum("mka,mla->mkl", gphi, Gg)
                    r += wa * wb * np.einsum("mka,ma->mk", Gg, kap)
            self.amode = (-np.linalg.solve(A, r[..., None])[..., 0]).reshape(n, n, K)
            self.interior = K

    # ---------------------------------------------------------------------------
    def _zh(self, I, s):
        """z_h' = spline + tangent fix + edge bubbles, in scaled variables.
        I (M,2) integer cell index (unwrapped), s (M,2) local coords in [0,1].
        Returns z (M,4) and local derivatives dz (M,4,2) = h * d z / d xi."""
        n, h = self.n, self.h
        xi = (I + s) * h
        val, gr = self.sp(xi, 1)
        z = val.copy()
        z[:, :2] += xi
        dz = gr * h
        dz[:, 0, 0] += h
        dz[:, 1, 1] += h
        i, j = I[:, 0] % n, I[:, 1] % n
        ip, jp = (i + 1) % n, (j + 1) % n
        Hs, dHs = hermite(s[:, 0])
        Ht, dHt = hermite(s[:, 1])
        if self.dT.any():
            for a, ia in ((0, i), (1, ip)):
                for b, jb in ((0, j), (1, jp)):
                    d0, d1 = self.dT[ia, jb, :, 0], self.dT[ia, jb, :, 1]
                    z += d0 * (Hs[2 + a] * Ht[b])[:, None] + d1 * (Hs[a] * Ht[2 + b])[:, None]
                    dz[:, :, 0] += d0 * (dHs[2 + a] * Ht[b])[:, None] + d1 * (dHs[a] * Ht[2 + b])[:, None]
                    dz[:, :, 1] += d0 * (Hs[2 + a] * dHt[b])[:, None] + d1 * (Hs[a] * dHt[2 + b])[:, None]
        if self.ah.any() or self.av.any():
            bs, dbs = bubble(s[:, 0])
            bt, dbt = bubble(s[:, 1])
            a0, a1 = self.ah[i, j], self.ah[i, jp]
            z += bs[:, None] * (Ht[0][:, None] * a0 + Ht[1][:, None] * a1)
            dz[:, :, 0] += dbs[:, None] * (Ht[0][:, None] * a0 + Ht[1][:, None] * a1)
            dz[:, :, 1] += bs[:, None] * (dHt[0][:, None] * a0 + dHt[1][:, None] * a1)
            v0, v1 = self.av[i, j], self.av[ip, j]
            z += bt[:, None] * (Hs[0][:, None] * v0 + Hs[1][:, None] * v1)
            dz[:, :, 0] += bt[:, None] * (dHs[0][:, None] * v0 + dHs[1][:, None] * v1)
            dz[:, :, 1] += dbt[:, None] * (Hs[0][:, None] * v0 + Hs[1][:, None] * v1)
        return z, dz

    @staticmethod
    def _modes(s, K):
        """Interior modes phi_k = b(s) b(t) m_k(s, t) and local gradients."""
        bs, dbs = bubble(s[:, 0])
        bt, dbt = bubble(s[:, 1])
        a, b = s[:, 0] - 0.5, s[:, 1] - 0.5
        one, zero = np.ones_like(a), np.zeros_like(a)
        mons = [(one, zero, zero), (a, one, zero), (b, zero, one), (a * b, b, a),
                (a * a, 2 * a, zero), (b * b, zero, 2 * b), (a * a * b, 2 * a * b, a * a), (a * b * b, b * b, 2 * a * b),
                (a * a * b * b, 2 * a * b * b, 2 * a * a * b)][:K]
        phi = np.stack([bs * bt * m for m, _, _ in mons], axis=1)
        g = np.stack([np.stack([dbs * bt * m + bs * bt * ms, bs * dbt * m + bs * bt * mt], axis=1) for m, ms, mt in mons], axis=1)
        return phi, g

    # ---------------------------------------------------------------------------
    def _potential(self, I, s):
        """Coons potential W_h and its local gradient at (cell I, local s)."""
        n, h = self.n, self.h
        i, j = I[:, 0] % n, I[:, 1] % n
        ip, jp = (i + 1) % n, (j + 1) % n
        us, ut = 2.0 * s[:, 0] - 1.0, 2.0 * s[:, 1] - 1.0
        ev = lambda c, u, a, b: Ch.chebval(u, c[:, a, b], tensor=False)
        T0, T1 = self.tab[0], self.tab[1]
        c0 = self.c[0] * h
        c1 = self.c[1] * h
        base = self.c[0] * I[:, 0] * h + self.c[1] * I[:, 1] * h
        S00 = self.Sp[i, j] + base
        S10 = self.Sp[ip, j] + base + c0
        S01 = self.Sp[i, jp] + base + c1
        S11 = self.Sp[ip, jp] + base + c0 + c1
        # horizontal edges (functions of s): bottom (i, j), top (i, jp)
        Wb, Wt_ = S00 + ev(T0["I"], us, i, j), S01 + ev(T0["I"], us, i, jp)
        Bb, Bt_ = ev(T0["t"], us, i, j), ev(T0["t"], us, i, jp)
        Nb, Nt_ = ev(T0["n"], us, i, j), ev(T0["n"], us, i, jp)
        dNb, dNt_ = ev(T0["dn"], us, i, j), ev(T0["dn"], us, i, jp)
        # vertical edges (functions of t): left (i, j), right (ip, j)
        Wl, Wr = S00 + ev(T1["I"], ut, i, j), S10 + ev(T1["I"], ut, ip, j)
        Bl, Br = ev(T1["t"], ut, i, j), ev(T1["t"], ut, ip, j)
        Nl, Nr = ev(T1["n"], ut, i, j), ev(T1["n"], ut, ip, j)
        dNl, dNr = ev(T1["dn"], ut, i, j), ev(T1["dn"], ut, ip, j)
        Hs, dHs = hermite(s[:, 0])
        Ht, dHt = hermite(s[:, 1])
        W = (Ht[0] * Wb + Ht[1] * Wt_ + Ht[2] * Nb + Ht[3] * Nt_
             + Hs[0] * Wl + Hs[1] * Wr + Hs[2] * Nl + Hs[3] * Nr)
        dWs = (Ht[0] * Bb + Ht[1] * Bt_ + Ht[2] * dNb + Ht[3] * dNt_
               + dHs[0] * Wl + dHs[1] * Wr + dHs[2] * Nl + dHs[3] * Nr)
        dWt = (dHt[0] * Wb + dHt[1] * Wt_ + dHt[2] * Nb + dHt[3] * Nt_
               + Hs[0] * Bl + Hs[1] * Br + Hs[2] * dNl + Hs[3] * dNr)
        for a, ia, Sa in ((0, i, (S00, S01)), (1, ip, (S10, S11))):
            for b, jb in ((0, j), (1, jp)):
                w, ws, wt, wst = Sa[b], self.Ws[ia, jb], self.Wt[ia, jb], self.Wst[ia, jb]
                W -= w * Hs[a] * Ht[b] + ws * Hs[2 + a] * Ht[b] + wt * Hs[a] * Ht[2 + b] + wst * Hs[2 + a] * Ht[2 + b]
                dWs -= w * dHs[a] * Ht[b] + ws * dHs[2 + a] * Ht[b] + wt * dHs[a] * Ht[2 + b] + wst * dHs[2 + a] * Ht[2 + b]
                dWt -= w * Hs[a] * dHt[b] + ws * Hs[2 + a] * dHt[b] + wt * Hs[a] * dHt[2 + b] + wst * Hs[2 + a] * dHt[2 + b]
        dW = np.stack([dWs, dWt], axis=1)
        if self.interior:
            phi, g = self._modes(s, self.interior)
            am = self.amode[i, j]
            W = W + np.sum(am * phi, axis=1)
            dW = dW + np.einsum("mk,mka->ma", am, g)
        return W, dW

    # ---------------------------------------------------------------------------
    def eval_local(self, I, s, diag=False):
        """Carrier at (cell I, local s).  Returns x, p, Sigma (physical units)
        and optionally a diagnostics dict."""
        n = self.n
        z, dz = self._zh(I, s)
        i, j = I[:, 0] % n, I[:, 1] % n
        W, dW = self._potential(I, s)
        if not self.lagrangian:
            # skeleton-exact carrier: z_h' with the Coons action (exact edge
            # actions, zero net omega-flux per cell, NOT pointwise Lagrangian)
            out = (z[:, :2], z[:, 2:] * self.lp, W * self.lp)
            return out + ({},) if diag else out
        beta = np.einsum("mc,mca->ma", z[:, 2:], dz[:, :2, :])
        kappa = dW - beta
        Ur, Ui = self.Ur[i, j], self.Ui[i, j]
        DY = np.einsum("mck,mca->mka", Ur, dz[:, :2, :]) + np.einsum("mck,mca->mka", Ui, dz[:, 2:, :])
        ok = self.valid[i, j]
        DY[~ok] = np.eye(2)
        dP = np.linalg.solve(np.swapaxes(DY, 1, 2), kappa[..., None])[..., 0]
        dP[~ok] = 0.0  # fallback cells: carrier = z_h' (Sigma = W_h there)
        x = z[:, :2] - np.einsum("mck,mk->mc", Ui, dP)
        p = z[:, 2:] + np.einsum("mck,mk->mc", Ur, dP)
        # Q_M = (P.Y - p.x)/2 ; Y is unchanged by the correction, P -> P + dP
        Y = np.einsum("mck,mc->mk", Ur, z[:, :2]) + np.einsum("mck,mc->mk", Ui, z[:, 2:])
        Ph = -np.einsum("mck,mc->mk", Ui, z[:, :2]) + np.einsum("mck,mc->mk", Ur, z[:, 2:])
        Qh = 0.5 * (np.sum(Ph * Y, axis=1) - np.sum(z[:, 2:] * z[:, :2], axis=1))
        Qc = 0.5 * (np.sum((Ph + dP) * Y, axis=1) - np.sum(p * x, axis=1))
        Sig = W + Qh - Qc
        out = (x, p * self.lp, Sig * self.lp)
        if diag:
            sv = np.linalg.svd(DY, compute_uv=False)
            svx = np.linalg.svd(dz[:, :2, :], compute_uv=False)
            d = dict(sig_min=sv[:, -1] * n, cond=sv[:, 0] / np.maximum(sv[:, -1], 1e-300),
                     detDY=np.linalg.det(DY) * n * n, detDx=np.linalg.det(dz[:, :2, :]) * n * n,
                     sig_min_x=svx[:, -1] * n, corr=np.linalg.norm(dP, axis=1), zh=z)
            return out + (d,)
        return out

    def __call__(self, xi, diag=False):
        u = np.asarray(xi, float) * self.n
        I = np.floor(u).astype(np.int64)
        return self.eval_local(I, u - I, diag)

    def chart_report(self, sub=3):
        """Per-cell chart conditioning sampled on a sub x sub lattice per cell."""
        n = self.n
        g = (np.arange(sub) + 0.5) / sub
        sig, cond, sx, det = [], [], [], []
        for a in g:
            for b in g:
                s = np.tile([a, b], (n * n, 1))
                *_, d = self.eval_local(self._I, s, diag=True)
                sig.append(d["sig_min"]); cond.append(d["cond"]); sx.append(d["sig_min_x"]); det.append(d["detDx"])
        sig, cond, sx, det = map(np.array, (sig, cond, sx, det))
        return dict(sig_min=sig.min(axis=0).reshape(n, n), cond_max=cond.max(axis=0).reshape(n, n),
                    sig_min_x=sx.min(axis=0).reshape(n, n),
                    folded=((det.min(axis=0) < 0) & (det.max(axis=0) > 0)).reshape(n, n),
                    detDx_min=det.min(axis=0).reshape(n, n))

    def memory_floats(self):
        """Number of stored float64 values (carrier tables + spline)."""
        t = sum(v.size for ax in (0, 1) for v in self.tab[ax].values())
        return int(self.sp.c.size + t + self.dT.size + self.ah.size + self.av.size + self.Ur.size * 2 + 4 * self.n ** 2)
