"""SHAI-HULUD reference prototype, d = 1 (two-dimensional phase space).

A cold phase-space sheet in T*S^1 is a closed curve
    q -> (x(q), p(q)),  q in [0, 2 pi),  x(q + 2 pi) = x(q) + 2 pi.
The discrete state carries, per vertex, the *Legendrian* data
    q_i  (Lagrangian label; carries the mass measure),
    x_i, p_i  (phase-space position),
    S_i  (action / generating function: dS = p dx along the sheet),
with wrap constants  q -> q + 2 pi, x -> x + 2 pi, S -> S + c,  c = oint p dx.

Everything in this file is research/reference code for the manuscript in
../manuscript/shai_hulud.tex.  It is intentionally not part of the CHUI
production tree and is not wired into the CMake build.

Conventions
-----------
omega(a, b) = a_x b_p - a_p b_x  (canonical symplectic form dx ^ dp),
theta = p dx,  oint_{boundary D} theta = - int_D omega  for counter-clockwise loops.
"""
from __future__ import annotations

import numpy as np

LD = np.longdouble
TWO_PI = 2.0 * np.pi


# --------------------------------------------------------------------------
# Hamiltonians  H(x, p, t) = alpha(t) p^2 / 2 + beta(t) V(x, t)
# --------------------------------------------------------------------------
class Pendulum:
    """Autonomous pendulum H = p^2/2 - kappa cos x."""

    name = "pendulum"

    def __init__(self, kappa=1.0):
        self.kappa = kappa

    def V(self, x, t):
        return -self.kappa * np.cos(x)

    def dV(self, x, t):
        return self.kappa * np.sin(x)

    def d2V(self, x, t):
        return self.kappa * np.cos(x)

    def drift_coef(self, t0, t1):
        return t1 - t0

    def kick_coef(self, t0, t1):
        return t1 - t0

    def energy(self, x, p, t):
        return 0.5 * p * p + self.V(x, t)


class DrivenPendulum(Pendulum):
    """Non-autonomous, chaotic: H = p^2/2 - kappa cos x - eps cos(x - Om t)."""

    name = "driven"

    def __init__(self, kappa=1.0, eps=0.25, Om=1.7):
        super().__init__(kappa)
        self.eps = eps
        self.Om = Om

    def V(self, x, t):
        return -self.kappa * np.cos(x) - self.eps * np.cos(x - self.Om * t)

    def dV(self, x, t):
        return self.kappa * np.sin(x) + self.eps * np.sin(x - self.Om * t)

    def d2V(self, x, t):
        return self.kappa * np.cos(x) + self.eps * np.cos(x - self.Om * t)


class ComovingEdS:
    """Comoving Hamiltonian with explicit time dependence (Einstein-de Sitter),

        H = p^2 / (2 a(t)^2) + phi(x),   a(t) = (t / t0)^(2/3),

    in cosmic time t with canonical momentum p = a^2 dx/dt.  The peculiar
    potential phi is prescribed (external), phi(x) = -kappa cos x; a
    time-independent phi is what linear growth (delta ~ a) gives in these
    variables.  The drift coefficient is the exact integral of a^-2.
    """

    name = "comoving_eds"

    def __init__(self, kappa=0.25, t0=1.0):
        self.kappa = kappa
        self.t0 = t0

    def a(self, t):
        return (t / self.t0) ** (2.0 / 3.0)

    def V(self, x, t):
        return -self.kappa * np.cos(x)

    def dV(self, x, t):
        return self.kappa * np.sin(x)

    def d2V(self, x, t):
        return self.kappa * np.cos(x)

    def drift_coef(self, t0, t1):
        # int a^-2 dt = int (t/t0)^(-4/3) dt = -3 t0^(4/3) t^(-1/3)
        f = lambda t: -3.0 * self.t0 ** (4.0 / 3.0) * t ** (-1.0 / 3.0)
        return f(t1) - f(t0)

    def kick_coef(self, t0, t1):
        return t1 - t0

    def energy(self, x, p, t):
        return 0.5 * p * p / self.a(t) ** 2 + self.V(x, t)


# --------------------------------------------------------------------------
# Exact-symplectic KDK step with exact action transport
# --------------------------------------------------------------------------
def kdk_step(model, x, p, S, t, dt):
    """One kick-drift-kick step.  Each sub-map is an exact symplectic map of
    T*X with an explicit generating function, and S is advanced by exactly
    that generating function, so that dS = p dx holds on the transported
    sheet to round-off (Theorem 'action transport' of the manuscript)."""
    th = t + 0.5 * dt
    t1 = t + dt
    K = model.kick_coef(t, th)
    S = S - K * model.V(x, t).astype(LD)
    p = p - K * model.dV(x, t)
    D = model.drift_coef(t, t1)
    # Round-off-consistent action transport: use the displacement actually
    # applied in floating point, S += p (x_new - x) - D p^2 / 2  (equal to
    # D p^2 / 2 in exact arithmetic).  Otherwise the rounding of x enters the
    # edge actions as p * ulp(x), uncorrelated between neighbouring vertices.
    x_new = x + D * p
    pl = p.astype(LD)
    S = S + pl * (x_new.astype(LD) - x.astype(LD)) - (0.5 * D) * (pl * pl)
    x = x_new
    K = model.kick_coef(th, t1)
    S = S - K * model.V(x, t1).astype(LD)
    p = p - K * model.dV(x, t1)
    return x, p, S


def integrate(model, x, p, S, t, dt, nsteps):
    for _ in range(nsteps):
        x, p, S = kdk_step(model, x, p, S, t, dt)
        t += dt
    return x, p, S, t


# --------------------------------------------------------------------------
# Initial condition (exactly Lagrangian, explicit generating function)
# --------------------------------------------------------------------------
class SineIC:
    """x0(q) = q,  p0(q) = u0 + A sin q   =>   S0(q) = u0 q - A cos q,
    c = oint p dx = 2 pi u0."""

    def __init__(self, u0=0.1, A=0.2):
        self.u0 = u0
        self.A = A
        self.c = TWO_PI * u0

    def __call__(self, q):
        q = np.asarray(q, dtype=float)
        x = q.copy()
        p = self.u0 + self.A * np.sin(q)
        S = (self.u0 * q.astype(LD)) - self.A * np.cos(q).astype(LD)
        return x, p, S

    def tangent(self, q):
        q = np.asarray(q, dtype=float)
        return np.ones_like(q), self.A * np.cos(q)


def oracle_tangent(model, ic, q, t0, dt, nsteps):
    """Oracle sample together with the exact tangent d(x, p)/dq of the
    transported carrier (the tangent map of the discrete flow)."""
    x, p, S = ic(q)
    tx, tp = ic.tangent(q)
    t = t0
    for _ in range(nsteps):
        th, t1 = t + 0.5 * dt, t + dt
        K = model.kick_coef(t, th)
        tp = tp - K * model.d2V(x, t) * tx
        S = S - K * model.V(x, t).astype(LD)
        p = p - K * model.dV(x, t)
        D = model.drift_coef(t, t1)
        tx = tx + D * tp
        x_new = x + D * p
        pl = p.astype(LD)
        S = S + pl * (x_new.astype(LD) - x.astype(LD)) - (0.5 * D) * (pl * pl)
        x = x_new
        K = model.kick_coef(th, t1)
        tp = tp - K * model.d2V(x, t1) * tx
        S = S - K * model.V(x, t1).astype(LD)
        p = p - K * model.dV(x, t1)
        t = t1
    return x, p, S, tx, tp


def oracle(model, ic, q, t0, dt, nsteps):
    """Exact sample of the transported carrier: integrate labels q from the
    initial condition with the *same* integrator and step sequence."""
    x, p, S = ic(q)
    x, p, S, _ = integrate(model, x, p, S, t0, dt, nsteps)
    return x, p, S


# --------------------------------------------------------------------------
# Sheet state
# --------------------------------------------------------------------------
class Sheet:
    def __init__(self, q, x, p, S, c):
        self.q = np.asarray(q, dtype=float)
        self.x = np.asarray(x, dtype=float)
        self.p = np.asarray(p, dtype=float)
        self.S = np.asarray(S, dtype=LD)
        self.c = float(c)

    @property
    def n(self):
        return self.q.size

    def copy(self):
        return Sheet(self.q.copy(), self.x.copy(), self.p.copy(), self.S.copy(), self.c)


def omega(ax, ap, bx, bp):
    return ax * bp - ap * bx


def _ext(a, wrap, g):
    """Periodic extension with g ghost entries on each side."""
    return np.concatenate([a[-g:] - wrap, a, a[:g] + wrap])


def _lagrange_weights(Q, qn):
    """Lagrange basis weights.  Q: (M, m) abscissae, qn: (M,)."""
    m = Q.shape[1]
    W = np.ones_like(Q)
    for k in range(m):
        for l in range(m):
            if l != k:
                W[:, k] *= (qn - Q[:, l]) / (Q[:, k] - Q[:, l])
    return W


# --------------------------------------------------------------------------
# Reconstructions (carriers).
#
# Method names:  "lin" / "cub" / "quin"  : Lagrange interpolation in the label
#                                           q on 2 / 4 / 6 nodes (no action).
#                "acMR" (ac21, ac41, ac43, ac61, ac65): action-corrected
#                interpolation.  Base curve = Lagrange interpolant z^L(q) on
#                M nodes; correction = beta(q) n0 with a nodal "bubble"
#                    beta(q) = pi_M(q) * sum_{k<R} c_k qhat^k ,
#                pi_M the nodal polynomial of the stencil (vanishing at all
#                stencil nodes) and n0 the chord conormal, omega(dz, n0) = 1.
#                The R coefficients are fixed by requiring that the
#                reconstructed arc reproduces the *exact* edge actions
#                    int_{q_a}^{q_b} p x_q dq = S_b - S_a
#                on R cells of the stencil (R = 1: the evaluated cell only;
#                R = M - 1: every cell of the stencil).
# --------------------------------------------------------------------------
METHODS = {
    "lin": (2, 0), "cub": (4, 0), "quin": (6, 0),
    "ac21": (2, 1), "ac41": (4, 1), "ac43": (4, 3), "ac61": (6, 1), "ac65": (6, 5),
}
G = 3  # ghost width used everywhere

_GL_X, _GL_W = np.polynomial.legendre.leggauss(10)
_GL_X = 0.5 * (_GL_X + 1.0)
_GL_W = 0.5 * _GL_W


def _lagrange_dweights(Q, qn):
    """Derivatives of the Lagrange basis polynomials at qn."""
    m = Q.shape[1]
    D = np.zeros_like(Q)
    for k in range(m):
        acc = np.zeros(Q.shape[0])
        for i in range(m):
            if i == k:
                continue
            term = 1.0 / (Q[:, k] - Q[:, i])
            for l in range(m):
                if l != k and l != i:
                    term = term * (qn - Q[:, l]) / (Q[:, k] - Q[:, l])
            acc += term
        D[:, k] = acc
    return D


def _pihat(Uh, uh):
    """Nodal polynomial prod_m (uh - Uh_m) and its derivative."""
    m = Uh.shape[1]
    d = uh[:, None] - Uh
    pi = np.prod(d, axis=1)
    dpi = np.zeros_like(uh)
    for k in range(m):
        dpi += np.prod(np.delete(d, k, axis=1), axis=1)
    return pi, dpi


class CellData:
    """Per-cell reconstruction data for every cell [q_j, q_{j+1}]."""

    def __init__(self, sheet, method, sig=(1.0, 1.0)):
        m, r = METHODS[method]
        self.m, self.r, self.h = m, r, m // 2
        h = self.h
        n = sheet.n
        qe = _ext(sheet.q, TWO_PI, G)
        xe = _ext(sheet.x, TWO_PI, G)
        pe = _ext(sheet.p, 0.0, G)
        Se = _ext(sheet.S, LD(sheet.c), G)
        j = np.arange(n) + G
        idx = j[:, None] + np.arange(1 - h, h + 1)[None, :]
        self.U = qe[idx] - qe[j][:, None]
        self.X = xe[idx] - xe[j][:, None]
        self.P = pe[idx]
        self.Srel = (Se[idx] - Se[j][:, None]).astype(float)
        self.xj, self.Sj = xe[j], Se[j]
        self.dq = self.U[:, h]
        dx, dp = self.X[:, h], self.P[:, h] - self.P[:, h - 1]
        sx2, sp2 = sig[0] ** 2, sig[1] ** 2
        nrm = dx * dx / sx2 + dp * dp / sp2
        self.n0x, self.n0p = (-dp / sp2) / nrm, (dx / sx2) / nrm
        self.Uh = self.U / self.dq[:, None]
        self.coef = np.zeros((n, r))
        if r:
            cells = [(h - 1, h)] if r == 1 else [(k, k + 1) for k in range(m - 1)]
            A = np.zeros((n, r, r))
            rhs = np.zeros((n, r))
            for ci, (a, b) in enumerate(cells):
                ua, ub = self.U[:, a], self.U[:, b]
                acc = np.zeros(n)
                for tg, wg in zip(_GL_X, _GL_W):
                    u = ua + tg * (ub - ua)
                    xL, pL, dxL, dpL = self._base(np.arange(n), u)
                    w = wg * (ub - ua)
                    acc += w * pL * dxL
                    pi, _ = _pihat(self.Uh, u / self.dq)
                    om = dxL * self.n0p - dpL * self.n0x
                    for k in range(r):
                        A[:, ci, k] += w * pi * (u / self.dq) ** k * om
                rhs[:, ci] = (self.Srel[:, b] - self.Srel[:, a]) - acc
            self.coef = np.linalg.solve(A, rhs[..., None])[..., 0]

    def _base(self, c, u):
        """Lagrange base curve (relative x, p) and u-derivatives in cells c."""
        Wt = _lagrange_weights(self.U[c], u)
        dWt = _lagrange_dweights(self.U[c], u)
        return (np.sum(Wt * self.X[c], axis=1), np.sum(Wt * self.P[c], axis=1),
                np.sum(dWt * self.X[c], axis=1), np.sum(dWt * self.P[c], axis=1))

    def curve(self, c, u):
        """Corrected carrier (relative x, p) and u-derivatives in cells c."""
        xL, pL, dxL, dpL = self._base(c, u)
        if self.r:
            dq = self.dq[c]
            uh = u / dq
            pi, dpi = _pihat(self.Uh[c], uh)
            poly = np.zeros_like(u)
            dpoly = np.zeros_like(u)
            for k in range(self.r):
                poly += self.coef[c, k] * uh**k
                if k:
                    dpoly += k * self.coef[c, k] * uh ** (k - 1)
            beta = pi * poly
            dbeta = (dpi * poly + pi * dpoly) / dq
            xL, pL = xL + beta * self.n0x[c], pL + beta * self.n0p[c]
            dxL, dpL = dxL + dbeta * self.n0x[c], dpL + dbeta * self.n0p[c]
        return xL, pL, dxL, dpL

    def action(self, c, u):
        """int_0^u p x_u du along the corrected carrier of cells c."""
        acc = np.zeros(u.shape, dtype=LD)
        for tg, wg in zip(_GL_X, _GL_W):
            _, p, dx, _ = self.curve(c, tg * u)
            acc += (wg * u * p * dx).astype(LD)
        return acc

    def cell_actions(self):
        c = np.arange(self.dq.size)
        return self.action(c, self.dq)


def evaluate(sheet, qn, method, sig=(1.0, 1.0), cd=None):
    """Sample the reconstructed extended sheet (x, p, S) at labels qn."""
    qn = np.asarray(qn, dtype=float)
    cd = CellData(sheet, method, sig) if cd is None else cd
    qe = _ext(sheet.q, TWO_PI, G)
    c = np.clip(np.searchsorted(qe, qn, side="right") - 1 - G, 0, sheet.n - 1)
    u = qn - sheet.q[c]
    x, p, _, _ = cd.curve(c, u)
    x_new = cd.xj[c] + x
    if cd.r:
        # action along the corrected carrier, made consistent with the
        # rounding of the stored position (see kdk_step)
        rho = (x_new.astype(LD) - cd.xj[c].astype(LD)) - x.astype(LD)
        S = cd.Sj[c] + cd.action(c, u) + p.astype(LD) * rho
    else:  # plain Lagrange family: the action is interpolated like any field
        Wt = _lagrange_weights(cd.U[c], u)
        S = cd.Sj[c] + np.sum(Wt * cd.Srel[c], axis=1).astype(LD)
    return x_new, p, S


# --------------------------------------------------------------------------
# Monitor and one-dimensional optimal transport (monotone rearrangement)
# --------------------------------------------------------------------------
def chord_defects(sheet):
    """Symplectically invariant edge action defects
        eps_j = (S_{j+1} - S_j) - (p_j + p_{j+1}) (x_{j+1} - x_j) / 2,
    i.e. the signed phase-space area between the carrier arc and its chord."""
    xe = np.append(sheet.x, sheet.x[0] + TWO_PI)
    pe = np.append(sheet.p, sheet.p[0])
    Se = np.append(sheet.S, sheet.S[0] + LD(sheet.c))
    trap = (0.5 * (pe[:-1] + pe[1:]) * (xe[1:] - xe[:-1])).astype(LD)
    return ((Se[1:] - Se[:-1]) - trap).astype(float)


def affine_curvature(sheet):
    """kappa = omega(z_q, z_qq) at nodes from inscribed-triangle areas
    (a symplectically invariant three-point estimate)."""
    qe, xe, pe = (_ext(sheet.q, TWO_PI, 1), _ext(sheet.x, TWO_PI, 1), _ext(sheet.p, 0.0, 1))
    dm, dp_ = qe[1:-1] - qe[:-2], qe[2:] - qe[1:-1]
    area2 = omega(xe[1:-1] - xe[:-2], pe[1:-1] - pe[:-2], xe[2:] - xe[1:-1], pe[2:] - pe[1:-1])
    kap = area2 / (0.5 * dm * dp_ * (dm + dp_))
    speed = np.hypot(xe[2:] - xe[:-2], pe[2:] - pe[:-2]) / (dm + dp_)
    return kap, speed


def monitor_density(sheet, mode="tri", floor=0.2, arc=0.3, turn=0.3, smooth=4):
    """Monitor density M(q_j) > 0 on label space, unit mean.

    M = Maff + floor + arc * Marc + turn * Mturn  (each normalised to unit
    mean), with Maff = |kappa|^(1/3) the equi-affine arclength density,
    Marc = |z_q| the (scaled) Euclidean arclength density and
    Mturn = |kappa| / |z_q|^2 the turning-angle density."""
    q = sheet.q
    dq = np.diff(np.append(q, q[0] + TWO_PI))
    wq = 0.5 * (dq + np.roll(dq, 1))  # nodal quadrature weights in q
    kap, speed = affine_curvature(sheet)
    if mode == "tri":
        Maff = np.abs(kap) ** (1.0 / 3.0)
    elif mode == "act":
        eps = np.abs(chord_defects(sheet))
        cell = (12.0 * eps) ** (1.0 / 3.0) / dq
        Maff = 0.5 * (cell + np.roll(cell, 1))
    else:
        raise ValueError(mode)
    mean = lambda f: np.sum(f * wq) / TWO_PI
    Mturn = np.abs(kap) / speed**2
    M = Maff / mean(Maff) + floor + arc * speed / mean(speed) + turn * Mturn / mean(Mturn)
    for _ in range(smooth):
        M = 0.25 * np.roll(M, 1) + 0.5 * M + 0.25 * np.roll(M, -1)
    return M / mean(M)


def ot_labels(sheet, M, n_new=None):
    """Monotone rearrangement (the d = 1 optimal transport map with the
    label 0 pinned): returns labels q'_i = T(i / n_new) with
    T_# Lebesgue = M dq / (2 pi)."""
    from scipy.interpolate import PchipInterpolator

    n_new = sheet.n if n_new is None else n_new
    q = np.append(sheet.q, sheet.q[0] + TWO_PI)
    Me = np.append(M, M[0])
    C = np.concatenate([[0.0], np.cumsum(0.5 * (Me[1:] + Me[:-1]) * np.diff(q))])
    inv = PchipInterpolator(C / C[-1], q - q[0])
    lev = np.arange(n_new) / n_new
    return q[0] + inv(lev)


def equidistribution_quality(sheet, M):
    dq = np.diff(np.append(sheet.q, sheet.q[0] + TWO_PI))
    Mc = 0.5 * (M + np.roll(M, -1))
    w = Mc * dq
    return w.max() / w.mean()


def remap(sheet, method, mode="tri", n_new=None, sig=(1.0, 1.0), **mon):
    M = monitor_density(sheet, mode=mode, **mon)
    qn = ot_labels(sheet, M, n_new)
    x, p, S = evaluate(sheet, qn, method, sig)
    return Sheet(qn, x, p, S, sheet.c)


# --------------------------------------------------------------------------
# Invariants and diagnostics
# --------------------------------------------------------------------------
def carrier_period(sheet, method, sig=(1.0, 1.0)):
    """oint p dx over the reconstructed closed carrier, by Gauss quadrature of
    the piecewise-polynomial reconstruction (geometry only: the stored action
    enters only through the reconstruction itself)."""
    cd = CellData(sheet, method, sig)
    return float(np.sum(cd.cell_actions()))


def polygon_period(sheet):
    xe = np.append(sheet.x, sheet.x[0] + TWO_PI)
    pe = np.append(sheet.p, sheet.p[0])
    return float(np.sum((0.5 * (pe[:-1] + pe[1:]) * (xe[1:] - xe[:-1])).astype(LD)))


def total_mass(sheet):
    """mu = dq / (2 pi) on label space: cell masses are label increments."""
    dq = np.diff(np.append(sheet.q, sheet.q[0] + TWO_PI))
    return float(np.sum(dq.astype(LD)) / LD(TWO_PI)), float(dq.min())
