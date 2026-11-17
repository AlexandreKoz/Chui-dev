"""Periodic tensor-product interpolating splines on uniform n x n grids of the
unit torus (degree 1, 3 or 5), with exact rational basis polynomials and FFT
prefiltering.  Independent implementation (not derived from Fable's
sheet2d.PSpline2, which relies on scipy.ndimage)."""
from __future__ import annotations

from fractions import Fraction
from math import comb, factorial

import numpy as np


def _basis_poly(k):
    """Coefficients of the k+1 cardinal B-splines of odd degree k that are
    non-zero at fractional offset t in [0, 1):  w_j(t) = beta^k(t - j),
    j = -(k-1)/2 .. (k+1)/2.  Returns (js, P) with P[j, m] the coefficient of
    t^m (exact rational arithmetic, converted to float)."""
    js = list(range(-(k - 1) // 2, (k + 1) // 2 + 1))
    P = np.zeros((k + 1, k + 1))
    half = (k + 1) // 2
    for a, j in enumerate(js):
        coef = [Fraction(0)] * (k + 1)
        for m in range(k + 2):
            c = half - j - m  # term (t + c)_+^k, active on (0,1) iff c >= 0
            if c < 0:
                continue
            sgn = Fraction((-1) ** m * comb(k + 1, m), factorial(k))
            for r in range(k + 1):
                coef[r] += sgn * comb(k, r) * Fraction(c) ** (k - r)
        P[a] = [float(v) for v in coef]
    return np.array(js), P


_BASIS = {k: _basis_poly(k) for k in (1, 3, 5)}


def _weights(t, k, nder):
    """Basis weights and derivatives (w.r.t. t) up to order nder at offsets t."""
    js, P = _BASIS[k]
    pw = np.ones((k + 1,) + t.shape)
    for m in range(1, k + 1):
        pw[m] = pw[m - 1] * t
    out = []
    C = P.copy()
    for _ in range(nder + 1):
        out.append(np.tensordot(C, pw, axes=(1, 0)))  # (k+1, M)
        # differentiate the polynomial coefficients
        Cn = np.zeros_like(C)
        Cn[:, :-1] = C[:, 1:] * np.arange(1, k + 1)
        C = Cn
    return js, out


def _prefilter(f, k):
    """Interpolation prefilter along the first two axes (circulant, FFT)."""
    if k == 1:
        return f.copy()
    n = f.shape[0]
    taps = {3: {0: 4 / 6, 1: 1 / 6}, 5: {0: 66 / 120, 1: 26 / 120, 2: 1 / 120}}[k]
    m = np.arange(n)
    lam = np.zeros(n)
    for j, v in taps.items():
        lam += v * np.cos(2 * np.pi * j * m / n) * (1.0 if j == 0 else 2.0)
    c = np.fft.fft(f, axis=0) / lam.reshape((n,) + (1,) * (f.ndim - 1))
    c = np.fft.fft(c, axis=1) / lam.reshape((1, n) + (1,) * (f.ndim - 2))
    return np.real(np.fft.ifft(np.fft.ifft(c, axis=1), axis=0))


class PSpline2:
    """Interpolating periodic tensor spline of f[i, j, ...] = f(i/n, j/n)."""

    def __init__(self, f, order=3):
        f = np.asarray(f, dtype=float)
        self.n = f.shape[0]
        self.k = order
        self.trail = f.shape[2:]
        self.c = _prefilter(f.reshape(self.n, self.n, -1), order)

    def __call__(self, xi, deriv=0):
        """xi (M, 2) (any real values; periodic).  Returns value (M, ...);
        deriv>=1 adds the xi-gradient (M, ..., 2); deriv>=2 the Hessian."""
        n, k = self.n, self.k
        u = np.asarray(xi, dtype=float) * n
        i0 = np.floor(u).astype(np.int64)
        t = u - i0
        js, w1 = _weights(t[:, 0], k, deriv)
        _, w2 = _weights(t[:, 1], k, deriv)
        M, K = u.shape[0], self.c.shape[2]
        val = np.zeros((M, K))
        g = np.zeros((M, K, 2)) if deriv >= 1 else None
        H = np.zeros((M, K, 2, 2)) if deriv >= 2 else None
        ia = [(i0[:, 0] + j) % n for j in js]
        ib = [(i0[:, 1] + j) % n for j in js]
        for a in range(k + 1):
            for b in range(k + 1):
                cc = self.c[ia[a], ib[b]]
                val += (w1[0][a] * w2[0][b])[:, None] * cc
                if deriv >= 1:
                    g[:, :, 0] += (w1[1][a] * w2[0][b])[:, None] * cc
                    g[:, :, 1] += (w1[0][a] * w2[1][b])[:, None] * cc
                if deriv >= 2:
                    H[:, :, 0, 0] += (w1[2][a] * w2[0][b])[:, None] * cc
                    H[:, :, 1, 1] += (w1[0][a] * w2[2][b])[:, None] * cc
                    mm = (w1[1][a] * w2[1][b])[:, None] * cc
                    H[:, :, 0, 1] += mm
                    H[:, :, 1, 0] += mm
        out = [val.reshape((M,) + self.trail)]
        if deriv >= 1:
            out.append((g * n).reshape((M,) + self.trail + (2,)))
        if deriv >= 2:
            out.append((H * n * n).reshape((M,) + self.trail + (2, 2)))
        return out[0] if deriv == 0 else tuple(out)


def grid(n):
    g = np.arange(n) / n
    return np.stack(np.meshgrid(g, g, indexing="ij"), axis=-1)


# Cubic Hermite blends on [0, 1] and the quartic edge bubble (unit integral)
def hermite(t):
    """Returns (H0, H1, K0, K1) and their derivatives (dH0, dH1, dK0, dK1)."""
    t2, t3 = t * t, t * t * t
    H = (1 - 3 * t2 + 2 * t3, 3 * t2 - 2 * t3, t - 2 * t2 + t3, -t2 + t3)
    dH = (-6 * t + 6 * t2, 6 * t - 6 * t2, 1 - 4 * t + 3 * t2, -2 * t + 3 * t2)
    return H, dH


def bubble(s):
    """b(s) = 30 s^2 (1-s)^2 and b'(s).  b, b' vanish at both ends."""
    return 30 * s * s * (1 - s) ** 2, 60 * s * (1 - s) * (1 - 2 * s)
