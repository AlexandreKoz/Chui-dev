"""Conventions (single source of truth) and linear symplectic-frame utilities.

CONVENTIONS  (used by every module of this package)
---------------------------------------------------
* Momentum is the *specific* momentum p = dx/dt (velocity; per unit mass).
  Mass lives on label space as a measure mu(dq); it never multiplies the
  single-particle symplectic form.
* Liouville one-form        theta = p . dx = sum_i p_i dx_i
* Symplectic two-form       omega = d theta = sum_i dp_i ^ dx_i
  As a bilinear form on phase-space vectors a = (a_x, a_p), b = (b_x, b_p):
        omega(a, b) = a_p . b_x - a_x . b_p .
  (Fable's sheet1d.py uses the opposite sign, dx ^ dp; all statements below
  are sign-consistent within this package and are tested in
  tests/test_conventions.py.)
* Hamilton's equations      xdot = dH/dp,  pdot = -dH/dx   <=>   i_{X_H} omega = -dH.
* Action transport          dSigma/dt = p . xdot - H  along trajectories, so
  that dSigma_t = F_t^* theta on a sheet F_t(q) = (x_t(q), p_t(q)) whenever
  dSigma_0 = F_0^* theta (plus conserved periods on non-contractible loops).
* A sheet is Lagrangian iff F^* omega = 0, i.e. (Dx)^T Dp is symmetric,
  i.e. omega(d_a F, d_b F) = 0 for all label directions a, b.
* A linear symplectic frame is M = [[A, B], [C, D]] with  M^T Om M = Om,
        Y = A x + B p,   P = C x + D p,
  and Om = [[0, I], [-I, 0]] (so that omega(a, b) = -a^T Om b ... see om()).
  Equivalent block conditions:  A^T C, B^T D symmetric,  A^T D - C^T B = I.
* Frame action transition (exact, for ANY map into phase space):
        P . dY - p . dx = d Q_M(x, p),
        Q_M(x, p) = 1/2 x^T (C^T A) x + x^T (C^T B) p + 1/2 p^T (D^T B) p .
  Hence if dSigma = p . dx on a sheet, W_M := Sigma + Q_M(x, p) satisfies
  dW_M = P . dY, and wherever DY is invertible the sheet is the graph
  P = grad_Y W_M  (a generating function in the frame M).
"""
from __future__ import annotations

import numpy as np


def omega(a, b):
    """omega(a, b) = a_p . b_x - a_x . b_p for vectors stacked as (..., 2d)."""
    d = a.shape[-1] // 2
    return np.sum(a[..., d:] * b[..., :d], axis=-1) - np.sum(a[..., :d] * b[..., d:], axis=-1)


def Jmat(d):
    """Matrix Om with  a^T Om b = a_x . b_p - a_p . b_x  (= -omega(a, b))."""
    Z, I = np.zeros((d, d)), np.eye(d)
    return np.block([[Z, I], [-I, Z]])


def is_symplectic(M, tol=1e-12):
    d = M.shape[0] // 2
    Om = Jmat(d)
    return bool(np.max(np.abs(M.T @ Om @ M - Om)) < tol)


def blocks(M):
    d = M.shape[0] // 2
    return M[:d, :d], M[:d, d:], M[d:, :d], M[d:, d:]


def frame_inverse(M):
    """Inverse of a symplectic matrix: [[D^T, -B^T], [-C^T, A^T]]."""
    A, B, C, D = blocks(M)
    return np.block([[D.T, -B.T], [-C.T, A.T]])


def Q_frame(M, x, p):
    """Quadratic generating correction Q_M(x, p) (see module docstring)."""
    A, B, C, D = blocks(M)
    return (0.5 * np.einsum("...i,ij,...j->...", x, C.T @ A, x)
            + np.einsum("...i,ij,...j->...", x, C.T @ B, p)
            + 0.5 * np.einsum("...i,ij,...j->...", p, D.T @ B, p))


def unitary_frame(Tx, Tp):
    """Osculating unitary frame of a (nearly) Lagrangian tangent plane.

    Tx, Tp : (..., d, d) with columns d_a x, d_a p (a = label direction).
    Z = Tx + i Tp = U H (polar decomposition, U unitary).  Returns (Ur, Ui) =
    (Re U, Im U).  The frame is
        A = Ur^T, B = Ui^T, C = -Ui^T, D = Ur^T        (M in Sp(2d) & O(2d))
    which is exactly symplectic for ANY input (U is unitary by construction).
    In this frame DY = Re(H) and DP = Im(H); Im(H) = 0 iff the plane is
    Lagrangian, and then DY = H is symmetric positive definite with the
    singular values of the full tangent map [Tx; Tp].
    """
    Z = Tx + 1j * Tp
    V, _, Wh = np.linalg.svd(Z)
    U = V @ Wh
    return U.real.copy(), U.imag.copy()


def frame_matrix(Ur, Ui):
    """Full 2d x 2d matrix of the unitary frame (for tests / diagnostics)."""
    return np.block([[Ur.T, Ui.T], [-Ui.T, Ur.T]])


def arnold_frames(d):
    """The 2^d coordinate ('Arnold') charts: U = diag(u), u_k in {1, i}."""
    out = []
    for mask in range(2 ** d):
        u = np.array([1j if (mask >> k) & 1 else 1.0 for k in range(d)])
        U = np.diag(u)
        out.append((U.real.copy(), U.imag.copy()))
    return out
