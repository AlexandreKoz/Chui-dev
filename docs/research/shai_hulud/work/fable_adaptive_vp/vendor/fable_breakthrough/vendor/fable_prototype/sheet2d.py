"""SHAI-HULUD reference prototype, d = 2 (four-dimensional phase space),
single-chart regime (the sheet is a graph over configuration space, i.e.
before shell crossing).

Reference domain Xi = T^2 = [0,1)^2 with a uniform n x n node grid.
Per node:  q (label, 2), x (position, 2), p (momentum, 2), S (action).

Contents
--------
* PSpline2      periodic bicubic interpolating splines with derivatives
* hermite_eval  tensor-product Hermite elements (bicubic C^1, biquintic C^2)
* solve_ot      Monge-Ampere optimal transport on the torus (Newton-Krylov)
* Sheet2 + remap with five reconstructions (see Carrier2):
      "lin" / "cub" / "quin"   direct interpolation of (x, p, S) in xi
      "front3" / "front5"      x by a spline, generating function
               W = S o x^{-1} by a C^1 / C^2 Hermite element with the exact
               nodal gradients (D x)^T p, and p := (D x)^{-T} grad_xi W
               (exactly Lagrangian carrier).

Research/reference code for ../manuscript/shai_hulud.tex; not part of the
CHUI production tree.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.sparse.linalg import LinearOperator, gmres

TWO_PI = 2.0 * np.pi


# --------------------------------------------------------------------------
# Periodic bicubic interpolating spline on a uniform n x n grid of [0,1)^2
# --------------------------------------------------------------------------
def _beta(u, k):
    """Centred uniform B-spline of degree k."""
    from math import comb, factorial
    out = np.zeros_like(u)
    for m in range(k + 2):
        out += (-1) ** m * comb(k + 1, m) * np.maximum(u + 0.5 * (k + 1) - m, 0.0) ** k
    return out / factorial(k)


def _bs(t, k):
    """Weights, first and second derivatives of the k+1 B-splines of odd
    degree k that are non-zero at fractional offset t in [0, 1)."""
    js = np.arange(-(k - 1) // 2, (k + 1) // 2 + 1)
    u = t[None, :] - js[:, None]
    w = _beta(u, k)
    d = _beta(u + 0.5, k - 1) - _beta(u - 0.5, k - 1)
    d2 = _beta(u + 1.0, k - 2) - 2 * _beta(u, k - 2) + _beta(u - 1.0, k - 2)
    return js, w, d, d2


class PSpline2:
    """Interpolating periodic tensor-product spline (degree 3 or 5) of
    f[i, j, ...] = f(i/n, j/n)."""

    def __init__(self, f, order=3):
        f = np.asarray(f, dtype=float)
        self.n = f.shape[0]
        self.k = order
        self.trail = f.shape[2:]
        f = f.reshape(self.n, self.n, -1)
        self.c = np.stack([ndimage.spline_filter(f[..., k], order=order, mode="grid-wrap") for k in range(f.shape[2])], axis=-1)

    def __call__(self, xi, deriv=0):
        """xi: (M, 2).  Returns value (M, ...); with deriv>=1 also the
        gradient (M, ..., 2); with deriv>=2 also the Hessian (M, ..., 2, 2)."""
        n = self.n
        u = np.asarray(xi) * n
        i0 = np.floor(u).astype(int)
        t = u - i0
        js, w1, d1, e1 = _bs(t[:, 0], self.k)
        _, w2, d2, e2 = _bs(t[:, 1], self.k)
        M = u.shape[0]
        K = self.c.shape[2]
        val = np.zeros((M, K))
        g = np.zeros((M, K, 2)) if deriv >= 1 else None
        H = np.zeros((M, K, 2, 2)) if deriv >= 2 else None
        for a, ja in enumerate(js):
            ia = (i0[:, 0] + ja) % n
            for b, jb in enumerate(js):
                ib = (i0[:, 1] + jb) % n
                cc = self.c[ia, ib]
                val += (w1[a] * w2[b])[:, None] * cc
                if deriv >= 1:
                    g[:, :, 0] += (d1[a] * w2[b])[:, None] * cc * n
                    g[:, :, 1] += (w1[a] * d2[b])[:, None] * cc * n
                if deriv >= 2:
                    H[:, :, 0, 0] += (e1[a] * w2[b])[:, None] * cc * n * n
                    H[:, :, 1, 1] += (w1[a] * e2[b])[:, None] * cc * n * n
                    m = (d1[a] * d2[b])[:, None] * cc * n * n
                    H[:, :, 0, 1] += m
                    H[:, :, 1, 0] += m
        out = [val.reshape((M,) + self.trail)]
        if deriv >= 1:
            out.append(g.reshape((M,) + self.trail + (2,)))
        if deriv >= 2:
            out.append(H.reshape((M,) + self.trail + (2, 2)))
        return out[0] if deriv == 0 else tuple(out)


def grid(n):
    g = np.arange(n) / n
    return np.stack(np.meshgrid(g, g, indexing="ij"), axis=-1)


class TorusMap:
    """A degree-one map of the torus, m(xi) = xi + periodic(xi), by splines."""

    def __init__(self, m_nodes, order=3):
        n = m_nodes.shape[0]
        self.n = n
        self.sp = PSpline2(m_nodes - grid(n), order)

    @classmethod
    def from_potential(cls, phi, oversample=4):
        """Label map T = id + grad(phi) of an optimal-transport potential,
        evaluated spectrally on an oversampled grid (so that the mesh map is
        known to far better accuracy than the fields it parameterises)."""
        n = phi.shape[0]
        N = oversample * n
        ph = np.fft.fft2(phi)
        big = np.zeros((N, N), dtype=complex)
        kk = np.fft.fftfreq(n, d=1.0 / n).astype(int)
        big[np.ix_(kk % N, kk % N)] = ph * oversample**2
        k = 2 * np.pi * np.fft.fftfreq(N, d=1.0 / N)
        k1, k2 = np.meshgrid(k, k, indexing="ij")
        g = np.stack([np.real(np.fft.ifft2(1j * k1 * big)), np.real(np.fft.ifft2(1j * k2 * big))], axis=-1)
        self = cls.__new__(cls)
        self.n = N
        self.sp = PSpline2(g, 5)
        return self

    def __call__(self, xi, deriv=0):
        if deriv == 0:
            return xi + self.sp(xi)
        r = self.sp(xi, deriv)
        return (xi + r[0], r[1] + np.eye(2)) + tuple(r[2:])

    def invert(self, target, guess=None, iters=40, tol=1e-13):
        """Solve m(xi) = target by damped Newton iteration."""
        xi = target - self.sp(target) if guess is None else guess.copy()
        for _ in range(iters):
            val, J = self(xi, 1)
            res = val - target
            if np.max(np.abs(res)) < tol:
                break
            step = np.linalg.solve(J, res[..., None])[..., 0]
            nrm = np.linalg.norm(step, axis=1, keepdims=True)
            xi = xi - step * np.minimum(1.0, 0.1 / np.maximum(nrm, 1e-300))
        return xi, float(np.max(np.abs(self(xi) - target)))


# --------------------------------------------------------------------------
# Tensor-product Hermite elements on the uniform grid:
#   degree 3 (Bogner-Fox-Schmit, C^1) and degree 5 (C^2).
# --------------------------------------------------------------------------
def _herm(s, degree):
    """Hermite basis on [0,1]: B[alpha][end](s) reproduces the alpha-th
    derivative at end in {0, 1}; returns values and first derivatives."""
    if degree == 3:
        B = [[2 * s**3 - 3 * s**2 + 1, -2 * s**3 + 3 * s**2], [s**3 - 2 * s**2 + s, s**3 - s**2]]
        D = [[6 * s**2 - 6 * s, -6 * s**2 + 6 * s], [3 * s**2 - 4 * s + 1, 3 * s**2 - 2 * s]]
    else:
        B = [[1 - 10 * s**3 + 15 * s**4 - 6 * s**5, 10 * s**3 - 15 * s**4 + 6 * s**5],
             [s - 6 * s**3 + 8 * s**4 - 3 * s**5, -4 * s**3 + 7 * s**4 - 3 * s**5],
             [0.5 * s**2 - 1.5 * s**3 + 1.5 * s**4 - 0.5 * s**5, 0.5 * s**3 - s**4 + 0.5 * s**5]]
        D = [[-30 * s**2 + 60 * s**3 - 30 * s**4, 30 * s**2 - 60 * s**3 + 30 * s**4],
             [1 - 18 * s**2 + 32 * s**3 - 15 * s**4, -12 * s**2 + 28 * s**3 - 15 * s**4],
             [s - 4.5 * s**2 + 6 * s**3 - 2.5 * s**4, 1.5 * s**2 - 4 * s**3 + 2.5 * s**4]]
    return B, D


def hermite_eval(Dn, xi):
    """Evaluate the tensor-product Hermite interpolant with nodal derivative
    data Dn[alpha][beta] = d_1^alpha d_2^beta W at the nodes (alpha, beta in
    0..m, m = 1: bicubic C^1, m = 2: biquintic C^2) at points xi (M, 2).
    Returns value (M,) and xi-gradient (M, 2)."""
    m = len(Dn) - 1
    n = Dn[0][0].shape[0]
    h = 1.0 / n
    u = np.asarray(xi) * n
    i0 = np.floor(u).astype(int)
    Bs, dBs = _herm(u[:, 0] - i0[:, 0], 2 * m + 1)
    Bt, dBt = _herm(u[:, 1] - i0[:, 1], 2 * m + 1)
    val = np.zeros(u.shape[0])
    g1 = np.zeros(u.shape[0])
    g2 = np.zeros(u.shape[0])
    for a in range(2):
        ia = (i0[:, 0] + a) % n
        for b in range(2):
            ib = (i0[:, 1] + b) % n
            for al in range(m + 1):
                for be in range(m + 1):
                    f = Dn[al][be][ia, ib] * h ** (al + be)
                    val += f * Bs[al][a] * Bt[be][b]
                    g1 += f * dBs[al][a] * Bt[be][b] * n
                    g2 += f * Bs[al][a] * dBt[be][b] * n
    return val, np.stack([g1, g2], axis=1)


def _dfd(f, axis, n, order=6):
    """Periodic central difference with respect to xi (order 4 or 6)."""
    r = lambda k: np.roll(f, -k, axis) - np.roll(f, k, axis)
    if order == 4:
        return (8 * r(1) - r(2)) * n / 12.0
    return (45 * r(1) - 9 * r(2) + r(3)) * n / 60.0


def front_nodal_data(S, g, degree):
    """Nodal derivative table for the generating function from its exact
    nodal values S and exact nodal gradients g = (D x)^T p; the higher mixed
    derivatives are finite differences of g (they only enter at higher order)."""
    n = S.shape[0]
    g1, g2 = g[..., 0], g[..., 1]
    d = lambda f, ax: _dfd(f, ax, n)
    w11 = 0.5 * (d(g1, 1) + d(g2, 0))
    if degree == 3:
        return [[S, g2], [g1, w11]]
    w20, w02 = d(g1, 0), d(g2, 1)
    w21, w12 = d(w20, 1), d(w02, 0)
    w22 = 0.5 * (d(w21, 1) + d(w12, 0))
    return [[S, g2, w02], [g1, w11, w12], [w20, w21, w22]]


# --------------------------------------------------------------------------
# Optimal transport on the torus:  T = id + grad(phi),  T_# Leb = M dq,
#   M(T(xi)) det(I + D^2 phi(xi)) = 1,   xi |-> |xi|^2/2 + phi convex.
# --------------------------------------------------------------------------
def solve_ot(M_grid, n, tol=1e-10, max_newton=40, verbose=False):
    """M_grid: positive density on the regular m x m label grid (any m),
    normalised internally to unit mean.  Returns dict with the map at the
    n x n nodes, the potential and convergence information."""
    M_grid = M_grid / M_grid.mean()
    logM = PSpline2(np.log(M_grid))
    k = 2 * np.pi * np.fft.fftfreq(n, d=1.0 / n)
    k1, k2 = np.meshgrid(k, k, indexing="ij")
    lap = -(k1**2 + k2**2)
    lap[0, 0] = 1.0
    xi = grid(n).reshape(-1, 2)

    def derivs(phi):
        ph = np.fft.fft2(phi)
        g1 = np.real(np.fft.ifft2(1j * k1 * ph))
        g2 = np.real(np.fft.ifft2(1j * k2 * ph))
        h11 = np.real(np.fft.ifft2(-k1 * k1 * ph))
        h12 = np.real(np.fft.ifft2(-k1 * k2 * ph))
        h22 = np.real(np.fft.ifft2(-k2 * k2 * ph))
        return g1, g2, h11, h12, h22

    def residual(phi):
        g1, g2, h11, h12, h22 = derivs(phi)
        det = (1 + h11) * (1 + h22) - h12**2
        if det.min() <= 0 or (1 + h11).min() <= 0:
            return None
        Tq = (xi + np.stack([g1.ravel(), g2.ravel()], axis=1)) % 1.0
        lm, glm = logM(Tq, 1)
        F = np.log(det).ravel() + lm
        return F.reshape(n, n), (h11, h12, h22, det, glm.reshape(n, n, 2))

    def invlap(r):
        rh = np.fft.fft2(r)
        rh[0, 0] = 0.0
        return np.real(np.fft.ifft2(rh / lap))

    phi = np.zeros((n, n))
    F, aux = residual(phi)
    hist = [float(np.abs(F - F.mean()).max())]
    for it in range(max_newton):
        h11, h12, h22, det, glm = aux
        a11, a12, a22 = (1 + h22) / det, -h12 / det, (1 + h11) / det  # inverse of I + D^2 phi

        def matvec(v):
            psi = v.reshape(n, n)
            g1, g2, p11, p12, p22 = derivs(psi)
            out = a11 * p11 + 2 * a12 * p12 + a22 * p22 + glm[..., 0] * g1 + glm[..., 1] * g2
            return (out - out.mean()).ravel()

        A = LinearOperator((n * n, n * n), matvec=matvec, dtype=float)
        P = LinearOperator((n * n, n * n), matvec=lambda r: invlap(r.reshape(n, n)).ravel(), dtype=float)
        rhs = -(F - F.mean()).ravel()
        dphi, _ = gmres(A, rhs, M=P, rtol=1e-8, restart=60, maxiter=20)
        dphi = dphi.reshape(n, n)
        lam = 1.0
        while True:
            r = residual(phi + lam * dphi)
            if r is not None and np.abs(r[0] - r[0].mean()).max() < (1 - 1e-4 * lam) * hist[-1] + 1e-15:
                break
            lam *= 0.5
            if lam < 1e-6:
                break
        phi = phi + lam * dphi
        F, aux = r if r is not None else residual(phi)
        hist.append(float(np.abs(F - F.mean()).max()))
        if verbose:
            print(f"   MA newton {it}: residual {hist[-1]:.3e} (lam={lam:.3f})")
        if hist[-1] < tol:
            break
    g1, g2, h11, h12, h22 = derivs(phi)
    T = grid(n) + np.stack([g1, g2], axis=-1)
    det = (1 + h11) * (1 + h22) - h12**2
    tr = 2 + h11 + h22
    lam_min = 0.5 * (tr - np.sqrt(np.maximum(tr**2 - 4 * det, 0.0)))
    return dict(T=T, phi=phi, residual_history=hist, det_min=float(det.min()), det_max=float(det.max()),
                hess_eig_min=float(lam_min.min()), cost=float(np.mean(g1**2 + g2**2)),
                curl_displacement=0.0)


# --------------------------------------------------------------------------
# External Hamiltonian H = |p|^2/2 + V(x) on T^2 and exact-symplectic KDK
# --------------------------------------------------------------------------
class Potential2:
    """V(x) = -A [cos k x1 + cos k x2 + gamma cos k (x1 - x2)],  k = 2 pi."""

    def __init__(self, A=1.0 / TWO_PI**2, gamma=0.5):
        self.A, self.gamma = A, gamma

    def V(self, x):
        k = TWO_PI
        return -self.A * (np.cos(k * x[..., 0]) + np.cos(k * x[..., 1]) + self.gamma * np.cos(k * (x[..., 0] - x[..., 1])))

    def dV(self, x):
        k = TWO_PI
        s12 = self.gamma * np.sin(k * (x[..., 0] - x[..., 1]))
        return self.A * k * np.stack([np.sin(k * x[..., 0]) + s12, np.sin(k * x[..., 1]) - s12], axis=-1)


def kdk(pot, x, p, S, dt, nsteps):
    """Kick-drift-kick with exact action transport (cf. sheet1d.kdk_step)."""
    for _ in range(nsteps):
        S = S - 0.5 * dt * pot.V(x)
        p = p - 0.5 * dt * pot.dV(x)
        x_new = x + dt * p
        S = S + np.sum(p * (x_new - x), axis=-1) - 0.5 * dt * np.sum(p * p, axis=-1)
        x = x_new
        S = S - 0.5 * dt * pot.V(x)
        p = p - 0.5 * dt * pot.dV(x)
    return x, p, S


def oracle2(pot, q, dt, nsteps):
    """Exact samples of the transported carrier for the cold-at-rest sheet
    x0 = q, p0 = 0, S0 = 0 (exactly Lagrangian, zero Liouville class)."""
    return kdk(pot, q.copy(), np.zeros_like(q), np.zeros(q.shape[:-1]), dt, nsteps)


# --------------------------------------------------------------------------
# Sheet state and reconstructions
# --------------------------------------------------------------------------
class Sheet2:
    """phi is the optimal-transport potential of the current label map
    q = T(xi) = xi + grad(phi)(xi)  (None for the identity)."""

    def __init__(self, q, x, p, S, phi=None):
        self.q, self.x, self.p, self.S = q, x, p, S
        self.n = q.shape[0]
        self.phi = phi

    def label_map(self):
        return TorusMap(self.q) if self.phi is None else TorusMap.from_potential(self.phi)


METHODS2 = ("lin", "cub", "quin", "front3", "front5")


class Carrier2:
    """Continuous reconstruction (carrier) of a Sheet2.

    "lin"    bilinear interpolation of (x, p, S) in xi
    "cub"    bicubic spline interpolation of (x, p, S)
    "quin"   biquintic spline interpolation of (x, p, S)
    "front3" x: bicubic spline;  W: bicubic Hermite (C^1);  p = (Dx)^-T grad W
    "front5" x: biquintic spline; W: biquintic Hermite (C^2); p = (Dx)^-T grad W
    For the front carriers p . dx = dW identically: exactly Lagrangian."""

    def __init__(self, sh, method):
        self.method = method
        n = sh.n
        self.n = n
        self.qmap = sh.label_map()  # the same (spectral) label map for every method
        if method == "lin":
            self.sh = sh
            return
        order = 5 if method in ("quin", "front5") else 3
        self.xmap = TorusMap(sh.x, order)
        if method in ("cub", "quin"):
            self.psp = PSpline2(sh.p, order)
            self.Ssp = PSpline2(sh.S, order)
        elif method in ("front3", "front5"):
            xi = grid(n).reshape(-1, 2)
            _, Dx = self.xmap(xi, 1)  # Dx[m, i, a] = d x_i / d xi_a
            g = np.einsum("mia,mi->ma", Dx, sh.p.reshape(-1, 2)).reshape(n, n, 2)
            self.Dn = front_nodal_data(sh.S, g, order)
        else:
            raise ValueError(method)

    def _bilinear(self, f, xi):
        n = self.n
        u = xi * n
        i0 = np.floor(u).astype(int)
        s, t = u[:, 0] - i0[:, 0], u[:, 1] - i0[:, 1]
        if f.ndim == 3:
            s, t = s[:, None], t[:, None]
        i, j, ip, jp = i0[:, 0] % n, i0[:, 1] % n, (i0[:, 0] + 1) % n, (i0[:, 1] + 1) % n
        return (1 - s) * (1 - t) * f[i, j] + s * (1 - t) * f[ip, j] + (1 - s) * t * f[i, jp] + s * t * f[ip, jp]

    def __call__(self, xi):
        """Returns x (unwrapped), p, S at reference points xi (M, 2)."""
        if self.method == "lin":
            g = grid(self.n)
            x = xi + self._bilinear(self.sh.x - g, xi)
            return x, self._bilinear(self.sh.p, xi), self._bilinear(self.sh.S, xi)
        x, Dx = self.xmap(xi, 1)
        if self.method in ("cub", "quin"):
            return x, self.psp(xi), self.Ssp(xi)
        W, gW = hermite_eval(self.Dn, xi)
        p = np.linalg.solve(np.swapaxes(Dx, 1, 2), gW[..., None])[..., 0]  # (Dx)^T p = grad W
        return x, p, W


def cubic_form_norm2(sh):
    """Squared Frobenius norm of the symplectic cubic form
    C_abc = omega(d_a z, d_bc z) in label coordinates, at the nodes."""
    n = sh.n
    xi = grid(n).reshape(-1, 2)
    xm, pm, qm = TorusMap(sh.x), PSpline2(sh.p), sh.label_map()
    _, Dx, Hx = xm(xi, 2)
    _, Dp, Hp = pm(xi, 2)
    _, Dq = qm(xi, 1)
    C = np.einsum("mia,mibc->mabc", Dx, Hp) - np.einsum("mia,mibc->mabc", Dp, Hx)
    Ji = np.linalg.inv(Dq)  # d xi / d q
    Cq = np.einsum("mabc,mad,mbe,mcf->mdef", C, Ji, Ji, Ji)
    return np.sum(Cq**2, axis=(1, 2, 3)).reshape(n, n)


def monitor_on_labels(sh, m=None, floor=0.5, power=0.4, width=0.05):
    """Monitor density M(q) on the regular m x m label grid,
        M ~ (|C|^2 + (floor <|C|>)^2)^(power/2),  power = d/(3+d) = 2/5,
    mollified by a Gaussian of fixed label-space width (independent of the
    mesh): the regularity of the mesh map T is what the remap error sees."""
    n = sh.n
    m = n if m is None else m
    C2 = cubic_form_norm2(sh)
    Mn = (C2 + (floor * np.sqrt(C2).mean()) ** 2) ** (0.5 * power)
    xi_of_q, err = sh.label_map().invert(grid(m).reshape(-1, 2))
    Mq = np.exp(PSpline2(np.log(Mn))(xi_of_q % 1.0)).reshape(m, m)
    k = 2 * np.pi * np.fft.fftfreq(m, d=1.0 / m)
    k1, k2 = np.meshgrid(k, k, indexing="ij")
    Mq = np.real(np.fft.ifft2(np.fft.fft2(Mq) * np.exp(-0.5 * width**2 * (k1**2 + k2**2))))
    return Mq / Mq.mean(), err


def remap2(sh, method, T_new, phi_new):
    """Resample the sheet at the new labels T_new = xi + grad(phi_new) using
    the given reconstruction.  The label map used for point location is the
    same (spectral) one for every method, so only the field interpolation
    differs between methods."""
    n = sh.n
    car = Carrier2(sh, method)
    xi_old, err = car.qmap.invert(T_new.reshape(-1, 2))
    x, p, S = car(xi_old)
    return Sheet2(T_new.copy(), x.reshape(n, n, 2), p.reshape(n, n, 2), S.reshape(n, n), phi_new), err


def loop_circulation(car, centre, radius, npts=4096):
    """oint p . dx along the image of a circle in the reference domain.
    Composite trapezoid on npts and 2 npts chords, Richardson-extrapolated."""
    def trap(m):
        th = TWO_PI * np.arange(m) / m
        xi = centre[None, :] + radius * np.stack([np.cos(th), np.sin(th)], axis=1)
        x, p, _ = car(xi)
        dx = np.roll(x, -1, axis=0) - x
        return float(np.sum(0.5 * (p + np.roll(p, -1, axis=0)) * dx))
    return (4.0 * trap(2 * npts) - trap(npts)) / 3.0


def label_mass(q):
    """Total mass of the pulled-back measure: sum of the signed areas of the
    label-space images of the reference cells (uniform mu, unit total mass).
    Returns (total, min cell mass, max cell mass)."""
    n = q.shape[0]
    d = q - grid(n)  # periodic part
    ge = np.stack(np.meshgrid(np.arange(n + 1) / n, np.arange(n + 1) / n, indexing="ij"), axis=-1)
    idx = np.arange(n + 1) % n
    Q = ge + d[idx][:, idx]
    a, b, c, e = Q[:-1, :-1], Q[1:, :-1], Q[:-1, 1:], Q[1:, 1:]
    cross = lambda u, v: u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]
    area = 0.5 * (cross(b - a, e - a) + cross(e - a, c - a))
    return float(area.sum()), float(area.min()), float(area.max())
