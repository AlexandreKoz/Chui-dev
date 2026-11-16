"""Symbolic and numerical verification of the closed-form statements used in
the SHAI-HULUD manuscript (../manuscript/shai_hulud.tex).

Each check prints PASS/FAIL (or the measured numbers for the statements that
are quantitative) and everything is written to ../results/verification.json.

Usage:  python3 verify_propositions.py
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np
import sympy as sp

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
report = {}


def record(key, ok, **info):
    report[key] = {"pass": bool(ok), **info}
    print(f"[{'PASS' if ok else 'FAIL'}] {key}: " + ", ".join(f"{k}={v}" for k, v in info.items()))


# ---------------------------------------------------------------------------
# 1. P1 edge action defect:  trapezoid - int p dx = omega(z', z'') h^3 / 12 + O(h^4)
# ---------------------------------------------------------------------------
def check_edge_defect():
    h, s = sp.symbols("h s")
    a = sp.symbols("x1:5")
    b = sp.symbols("p0:5")
    # generic curve z(q) Taylor-expanded at q = 0, evaluated on [0, h]
    q = sp.symbols("q")
    x = sum(a[k - 1] * q**k / sp.factorial(k) for k in range(1, 5))
    p = b[0] + sum(b[k] * q**k / sp.factorial(k) for k in range(1, 5))
    exact = sp.integrate(p * sp.diff(x, q), (q, 0, h))
    trap = sp.Rational(1, 2) * (p.subs(q, 0) + p.subs(q, h)) * (x.subs(q, h) - x.subs(q, 0))
    d = sp.series(sp.expand(trap - exact), h, 0, 4).removeO()
    lead = sp.Rational(1, 12) * (a[0] * b[2] - b[1] * a[1]) * h**3  # omega(z', z'') h^3 / 12
    record("edge_defect_p1_leading_term", sp.simplify(d - lead) == 0, formula="omega(z',z'') h^3/12")


# ---------------------------------------------------------------------------
# 2. P1 triangle flux for a Lagrangian graph p = grad phi(x):
#    oint_{boundary tau} theta_PL = (1/4) [C(u,u,v) - C(u,v,v)] + O(h^4),  C = D^3 phi
# ---------------------------------------------------------------------------
def check_triangle_flux():
    x1, x2, h = sp.symbols("x1 x2 h")
    # generic quartic generating function
    cs = sp.symbols("c0:15")
    monos = [x1**i * x2**j for i in range(5) for j in range(5 - i)]
    phi = sum(c * m for c, m in zip(cs, monos))
    grad = [sp.diff(phi, x1), sp.diff(phi, x2)]
    u1, u2, v1, v2 = sp.symbols("u1 u2 v1 v2")
    V0 = (0, 0)
    V1 = (h * u1, h * u2)
    V2 = (h * v1, h * v2)

    def pv(V):
        return [g.subs({x1: V[0], x2: V[1]}) for g in grad]

    def trap(A, B):
        pa, pb = pv(A), pv(B)
        return sp.Rational(1, 2) * sum((pa[k] + pb[k]) * (B[k] - A[k]) for k in range(2))

    loop = sp.expand(trap(V0, V1) + trap(V1, V2) + trap(V2, V0))
    lead = sp.series(loop, h, 0, 4).removeO()

    def C(a, b, c):  # third derivative tensor at the origin
        tot = 0
        for i, j, k in itertools.product(range(2), repeat=3):
            d = sp.diff(phi, *[(x1, x2)[m] for m in (i, j, k)]).subs({x1: 0, x2: 0})
            tot += d * a[i] * b[j] * c[k]
        return tot

    U, Vv = (u1, u2), (v1, v2)
    pred = sp.Rational(1, 4) * (C(U, U, Vv) - C(U, Vv, Vv)) * h**3
    record("triangle_flux_p1_cubic_form", sp.expand(lead - pred) == 0,
           formula="(1/4)[C(u,u,v)-C(u,v,v)] h^3")
    # explicit example of the manuscript: phi = x1^2 x2, unit right triangle
    phi_ex = x1**2 * x2
    g = [sp.diff(phi_ex, x1), sp.diff(phi_ex, x2)]
    P = lambda V: [gg.subs({x1: V[0], x2: V[1]}) for gg in g]
    T = lambda A, B: sp.Rational(1, 2) * sum((P(A)[k] + P(B)[k]) * (B[k] - A[k]) for k in range(2))
    val = sp.simplify(T((0, 0), (h, 0)) + T((h, 0), (0, h)) + T((0, h), (0, 0)))
    record("triangle_flux_example_x1sq_x2", sp.simplify(val - h**3 / 2) == 0, value=str(val))


# ---------------------------------------------------------------------------
# 3. Rigidity of exactly isotropic P1 complexes on the Kuhn-triangulated lattice
# ---------------------------------------------------------------------------
def check_rigidity_symbol():
    t1, t2 = sp.symbols("tau1 tau2")  # tau_m = tan(theta_m / 2)
    t3 = (t1 + t2) / (1 - t1 * t2)  # tan((theta1 + theta2) / 2)
    cond = sp.factor(sp.simplify((t1 + t2 - t3) * (1 - t1 * t2)))
    record("rigidity_symbol_factorisation", sp.simplify(cond + t1 * t2 * (t1 + t2)) == 0, symbol=str(cond))


def lattice_kernel_dim(n, d, jitter=0.0, seed=0):
    """Kernel dimension of the linear system 'every edge action of the P1
    complex is exact':  (p_i + p_j) . (x_j - x_i) / 2 = S_j - S_i + c . dxi
    for unknowns (p_i, S_i, c) on the periodic Kuhn-triangulated n^d lattice."""
    rng = np.random.default_rng(seed)
    idx = np.arange(n**d).reshape((n,) * d)
    grid = np.stack(np.meshgrid(*[np.arange(n)] * d, indexing="ij"), axis=-1).reshape(-1, d)
    xpos = grid / n + jitter * (rng.random((n**d, d)) - 0.5) / n
    dirs = [np.array(a) for a in itertools.product((0, 1), repeat=d) if any(a)]
    nv = n**d
    nunk = nv * (d + 1) + d
    rows = []
    for a in dirs:
        nb = np.roll(idx, shift=tuple(-a), axis=tuple(range(d))).reshape(-1)
        for i in range(nv):
            j = nb[i]
            wrap = (grid[i] + a) // n  # crossing of the periodic boundary
            dx = xpos[j] + wrap - xpos[i]
            row = np.zeros(nunk)
            row[i * (d + 1): i * (d + 1) + d] += 0.5 * dx
            row[j * (d + 1): j * (d + 1) + d] += 0.5 * dx
            row[i * (d + 1) + d] += 1.0
            row[j * (d + 1) + d] -= 1.0
            row[nv * (d + 1):] -= a / n
            rows.append(row)
    A = np.array(rows)
    sv = np.linalg.svd(A, compute_uv=False)
    tol = 1e-9 * sv[0]
    rank = int(np.sum(sv > tol))
    return nunk - rank, nunk, len(rows)


def check_rigidity_numeric():
    out = {}
    ok = True
    for n in (5, 6, 8):
        k, nunk, nrow = lattice_kernel_dim(n, 2)
        out[f"d2_n{n}_regular"] = [k, nunk, nrow]
        ok &= (k == 3 * n)
        k2, _, _ = lattice_kernel_dim(n, 2, jitter=0.3)
        out[f"d2_n{n}_jittered"] = [k2, nunk, nrow]
    for n in (4, 5):
        k, nunk, nrow = lattice_kernel_dim(n, 3)
        out[f"d3_n{n}_regular"] = [k, nunk, nrow]
        k2, _, _ = lattice_kernel_dim(n, 3, jitter=0.3)
        out[f"d3_n{n}_jittered"] = [k2, nunk, nrow]
    record("rigidity_kernel_dimension_d2_equals_3n", ok, **out)


# ---------------------------------------------------------------------------
# 4. Sliding vertices along a convex polygon strictly loses action:
#    the midpoint polygon of the regular N-gon has area cos^2(pi/N) * area.
# ---------------------------------------------------------------------------
def check_polygon_area_loss():
    ok = True
    vals = {}
    for N in (8, 64, 512):
        th = 2 * np.pi * np.arange(N) / N
        z = np.stack([np.cos(th), np.sin(th)], axis=1)
        area = lambda w: 0.5 * np.sum(w[:, 0] * np.roll(w[:, 1], -1) - w[:, 1] * np.roll(w[:, 0], -1))
        mid = 0.5 * (z + np.roll(z, -1, axis=0))
        ratio = area(mid) / area(z)
        vals[f"N{N}"] = ratio
        ok &= abs(ratio - np.cos(np.pi / N) ** 2) < 1e-12
    record("polygon_midpoint_area_factor_cos2", ok, **vals)


# ---------------------------------------------------------------------------
# 5. Cubic Hermite front: Lebesgue constant in W^{1,infty} equals 2;
#    sub-arc area of the area-exact parabola equals eps * s^3 (Archimedes).
# ---------------------------------------------------------------------------
def check_hermite_constants():
    s = np.linspace(0, 1, 200001)
    a = (1 - s) * (1 - 3 * s)
    b = s * (3 * s - 2)
    c = 6 * s * (1 - s)
    lam = np.max(np.abs(a) + np.abs(b) + np.abs(c))
    record("hermite_front_lebesgue_constant", abs(lam - 2.0) < 1e-9, value=float(lam))
    ss, eps = sp.symbols("s eps")
    W = eps * (3 * ss**2 - 2 * ss**3)
    eta = sp.diff(W, ss)
    sub = sp.simplify(W - sp.Rational(1, 2) * ss * eta)
    record("parabola_subarc_area_s_cubed", sp.simplify(sub - eps * ss**3) == 0, value=str(sub))


# ---------------------------------------------------------------------------
# 6. Quadratic (P2) edge with a true mid-parameter tracer: action error O(h^5)
# ---------------------------------------------------------------------------
def check_p2_edge_order():
    f = lambda q: (q + 0.3 * np.sin(1.3 * q), np.cos(0.7 * q) + 0.2 * q**2)
    qq = np.linspace(0, 1, 200001)

    def exact(h):
        q = qq * h
        x, p = f(q)
        return np.trapezoid(p[:-1] * 0 + 0.5 * (p[1:] + p[:-1]), x[1:] - x[:-1]) if False else float(
            np.sum(0.5 * (p[1:] + p[:-1]) * (x[1:] - x[:-1])))

    def simpson(h):
        (x0, p0), (xm, pm), (x1, p1) = f(0.0), f(0.5 * h), f(h)
        return (p0 * (-3 * x0 + 4 * xm - x1) + 4 * pm * (x1 - x0) + p1 * (x0 - 4 * xm + 3 * x1)) / 6.0

    hs = np.array([0.4, 0.2, 0.1, 0.05])
    # reference by Richardson-clean composite rule on a very fine grid
    err = np.array([abs(simpson(h) - _fine_action(f, h)) for h in hs])
    orders = np.log2(err[:-1] / err[1:])
    record("p2_edge_action_superconvergence", np.all(orders > 4.7), orders=[round(float(o), 2) for o in orders])


def _fine_action(f, h, m=4096):
    # composite Gauss-Legendre for int_0^h p x' dq with a spectral-accurate derivative
    xg, wg = np.polynomial.legendre.leggauss(8)
    tot = 0.0
    for k in range(m // 64):
        a, b = k * h / (m // 64), (k + 1) * h / (m // 64)
        q = 0.5 * (a + b) + 0.5 * (b - a) * xg
        dq = 1e-6
        x1, _ = f(q + dq)
        x0, _ = f(q - dq)
        _, p = f(q)
        tot += 0.5 * (b - a) * np.sum(wg * p * (x1 - x0) / (2 * dq))
    return tot


# ---------------------------------------------------------------------------
# 7. A grid force that is not an exact gradient is not a symplectic kick:
#    curl of the CIC-interpolated finite-difference force vs gradient of a
#    C^1 (quadratic B-spline) interpolated potential.
# ---------------------------------------------------------------------------
def check_force_curl():
    out = {}
    rng = np.random.default_rng(1)
    for n in (16, 32, 64):
        h = 1.0 / n
        g = np.arange(n) * h
        X, Y = np.meshgrid(g, g, indexing="ij")
        phi = np.cos(2 * np.pi * X) * np.cos(2 * np.pi * Y) + 0.3 * np.sin(2 * np.pi * (X + 2 * Y))
        gx = -(np.roll(phi, -1, 0) - np.roll(phi, 1, 0)) / (2 * h)
        gy = -(np.roll(phi, -1, 1) - np.roll(phi, 1, 1)) / (2 * h)
        pts = rng.random((4000, 2))
        i = np.floor(pts[:, 0] / h).astype(int)
        j = np.floor(pts[:, 1] / h).astype(int)
        s = pts[:, 0] / h - i
        t = pts[:, 1] / h - j
        ip, jp = (i + 1) % n, (j + 1) % n
        # CIC (bilinear) interpolation of each force component; analytic curl
        dgx_dy = ((1 - s) * (gx[i, jp] - gx[i, j]) + s * (gx[ip, jp] - gx[ip, j])) / h
        dgy_dx = ((1 - t) * (gy[ip, j] - gy[i, j]) + t * (gy[ip, jp] - gy[i, jp])) / h
        dgx_dx = ((1 - t) * (gx[ip, j] - gx[i, j]) + t * (gx[ip, jp] - gx[i, jp])) / h
        curl = np.abs(dgx_dy - dgy_dx)
        out[f"n{n}_cic_fd_curl_over_grad"] = float(curl.max() / np.abs(dgx_dx).max())
    vals = list(out.values())
    orders = [float(np.log2(vals[k] / vals[k + 1])) for k in range(len(vals) - 1)]
    record("cic_fd_force_is_not_a_gradient", all(v > 1e-6 for v in vals), orders=orders, **out)


# ---------------------------------------------------------------------------
# 8. von Neumann analysis of the action-corrected remaps in the frozen
#    (linear graph) setting: uniform nodes X_j = j h, data (W_j, P_j), shift
#    by sigma * h.  Spectral radius of the 2x2 amplification matrix.
# ---------------------------------------------------------------------------
def amplification(m, sigma, theta):
    """Amplification matrix acting on (W / h, P) Fourier amplitudes for the
    scheme: base = Lagrange interpolation of P on m nodes, one bubble
    c * pi_m(X) fixed by the exact cell action W_{j+1} - W_j."""
    hw = m // 2
    nodes = np.arange(1 - hw, hw + 1).astype(float)
    xg, wg = np.polynomial.legendre.leggauss(12)

    def lag(x):
        L = np.ones((np.size(x), m))
        x = np.atleast_1d(x)
        for k in range(m):
            for l in range(m):
                if l != k:
                    L[:, k] *= (x - nodes[l]) / (nodes[k] - nodes[l])
        return L

    def pim(x):
        return np.prod(np.atleast_1d(x)[:, None] - nodes[None, :], axis=1)

    def integ(fun, a, b):
        x = 0.5 * (a + b) + 0.5 * (b - a) * xg
        return 0.5 * (b - a) * (wg[:, None] * fun(x)).sum(axis=0) if fun(x).ndim > 1 else 0.5 * (b - a) * np.sum(wg * fun(x))

    IL = integ(lag, 0.0, 1.0)  # int_0^1 L_k
    Ipi = integ(pim, 0.0, 1.0)
    ILs = integ(lag, 0.0, sigma)
    Ipis = integ(pim, 0.0, sigma)
    Ls = lag(sigma)[0]
    pis = pim(sigma)[0]
    ph = np.exp(1j * theta * nodes)  # P_{j+k} = Phat e^{i theta k}
    # c = [(W_{j+1} - W_j)/h - sum_k IL_k P_k] / Ipi
    cW = (np.exp(1j * theta) - 1.0) / Ipi
    cP = -np.sum(IL * ph) / Ipi
    A = np.zeros((2, 2), dtype=complex)
    # new W/h = W_j/h + sum ILs_k P_k + c Ipis
    A[0, 0] = 1.0 + cW * Ipis
    A[0, 1] = np.sum(ILs * ph) + cP * Ipis
    # new P = sum Ls_k P_k + c pis
    A[1, 0] = cW * pis
    A[1, 1] = np.sum(Ls * ph) + cP * pis
    return A


def check_von_neumann():
    out = {}
    for m in (2, 4, 6):
        rho = 0.0
        arg = None
        for sigma in np.linspace(0.0, 1.0, 101):
            for theta in np.linspace(1e-3, np.pi, 240):
                r = np.max(np.abs(np.linalg.eigvals(amplification(m, sigma, theta))))
                if r > rho:
                    rho, arg = r, (float(sigma), float(theta))
        out[f"ac{m}1_max_spectral_radius"] = float(rho)
        out[f"ac{m}1_argmax_sigma_theta"] = arg
    record("von_neumann_spectral_radius", all(v <= 1 + 1e-9 for k, v in out.items() if k.endswith("radius")), **out)


if __name__ == "__main__":
    check_edge_defect()
    check_triangle_flux()
    check_rigidity_symbol()
    check_rigidity_numeric()
    check_polygon_area_loss()
    check_hermite_constants()
    check_p2_edge_order()
    check_force_curl()
    check_von_neumann()
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "verification.json"), "w") as f:
        json.dump(report, f, indent=1, default=str)
    n_fail = sum(not v["pass"] for v in report.values())
    print(f"{len(report) - n_fail}/{len(report)} checks passed")
