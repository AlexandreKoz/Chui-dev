"""Composite hierarchical-surplus spline carrier on a block quadtree.

Definition (degree k = 3 or 5; B-splines beta of level m on the grid 1/N_m):

    F(xi) = sum_{m=0}^{top} S_m(xi),      S_m = sum_{j in T_m} c_{m,j} beta_m(xi - xi_j)

  S_0  the periodic interpolating spline through the level-0 nodes (this IS
       the baseline carrier `SplineCarrier`);
  S_m  (m >= 1) the unique spline in span{beta_m(. - xi_j): j in T_m} with
           S_m(xi_j) = z_j - sum_{l<m} S_l(xi_j)   at the active level-m nodes,
           S_m(xi_j) = 0                            at every other grid point of T_m,
       where T_m ("solve tiles") is the set of level-m grid points in the
       region covered by level-m blocks dilated by `halo_tiles` block rings.

Properties (proved in research/local_refinement_contract.md, tested in tests/):
  * F is C^{k-1} on the whole torus: every S_m is a finite sum of C^{k-1}
    B-splines.  There is no trace mismatch at coarse-fine interfaces at all.
  * F interpolates every active node of every level up to the "skirt leak"
    lambda^(halo_tiles*B) * surplus (lambda = 0.43 quintic, 0.27 cubic).
  * F is a linear function of the nodal data, and the linear system for S_m
    is symmetric positive definite with condition number <= 56 (k = 5).
  * Uniform limit: if all blocks of level 1..m exist, F is exactly the
    periodic interpolating spline of all nodes on the level-m grid.
  * Creating nodes with values F_old(xi_new) leaves F unchanged (surplus 0).

Evaluation uses the accumulated representation  A_m = sum_{l<=m} S_l  written
in level-m B-spline coefficients (a_m = P_m a_{m-1} + c_m, P_m = two-scale
subdivision) on tile sets D_m; a leaf block evaluates from A_{le}, le >= its
own level being the finest level whose correction can be non-zero on it.
Evaluation at the tensor quadrature lattice of a block is two small dense
matrix products.  Everything that depends only on the topology (index
windows, P_m, the SPD matrices, banded factors) is built once per regrid.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .kernels import COLLOC, DECAY, eval_weights, prefilter_periodic, two_scale_mask


def _dilate(mask, r):
    """Periodic Chebyshev dilation of a boolean (nb, nb) map by r tiles."""
    out = mask.copy()
    for _ in range(r):
        o = out
        for ax in (0, 1):
            o = o | np.roll(o, 1, ax) | np.roll(o, -1, ax)
        out = o
    return out


def _coarsen(mask):
    nb = mask.shape[0] // 2
    return mask.reshape(nb, 2, nb, 2).any(axis=(1, 3))


class _Level:
    """Topology-dependent data of one level m >= 1."""
    __slots__ = ("m", "N", "nb", "dkeys", "nD", "tsel", "nT", "P", "M", "Mn", "node_rows", "node_keys",
                 "perm", "chol", "crow", "lam", "lu", "tmap", "full", "dense_of_T")


class HierarchicalCarrier:
    def __init__(self, forest, k=5, halo_tiles=1, solver="adi", sweeps=1, tol=1e-13, max_iter=60):
        if k not in (3, 5):
            raise ValueError("degree must be 3 or 5")
        self.f, self.k, self.ht = forest, k, int(halo_tiles)
        # solver: "adi"   fixed number of alternating line sweeps from a zero guess
        #                   (a deterministic LINEAR operator; default, 1 sweep)
        #         "adi-tol" line sweeps to relative residual `tol` (verification)
        #         "splu"  sparse LU (verification)
        self.solver, self.sweeps, self.tol, self.max_iter = solver, int(sweeps), tol, max_iter
        self.B = forest.B
        if self.B < 8:
            raise ValueError("block size B >= 8 required by the stencil/tile logic")
        self.h = (k - 1) // 2
        self.version = None
        self.stats = dict(layout_s=0.0, update_s=0.0, solve_iters=[], n_layout=0)
        self._wcache = {}
        self.build_layout()

    # ------------------------------------------------------------------ layout
    def _flat(self, m, i, j):
        """Flat coefficient index of level-m grid point (i, j) (arrays, any
        integers; periodic).  Returns (flat, valid)."""
        if m == 0:
            n0 = self.f.n0
            return (i % n0) * n0 + (j % n0), np.ones(np.shape(i), dtype=bool)
        lv = self.lev[m]
        B = self.B
        i = i % lv.N
        j = j % lv.N
        pos = lv.tmap[i // B, j // B]                      # tile slot, -1 if the tile is not stored
        valid = pos >= 0
        return np.where(valid, pos, 0) * (B * B) + (i % B) * B + (j % B), valid

    def build_layout(self):
        t0 = time.perf_counter()
        f, B, k, h = self.f, self.B, self.k, self.h
        top = f.top_level()
        self.top = top
        self.lev = [None] * (top + 1)
        # -- tile maps -----------------------------------------------------------
        exist, T, need = [None] * (top + 1), [None] * (top + 1), [None] * (top + 1)
        for m in range(1, top + 1):
            nb = f.nbl(m)
            e = np.zeros((nb, nb), dtype=bool)
            blk = np.array(sorted(f.exists[m]), dtype=np.int64)
            e[blk[:, 0], blk[:, 1]] = True
            exist[m] = e
            T[m] = _dilate(e, self.ht)
            need[m] = _dilate(T[m], 1)
        # -- leaf table and evaluation levels --------------------------------------
        leaves = f.leaves()
        nl = len(leaves)
        self.leafL = np.array([q[0] for q in leaves], dtype=np.int64)
        self.leafI = np.array([q[1] for q in leaves], dtype=np.int64)
        self.leafJ = np.array([q[2] for q in leaves], dtype=np.int64)
        self.leaf_le = self.leafL.copy()
        evalreq = [None] + [np.zeros((f.nbl(m), f.nbl(m)), dtype=bool) for m in range(1, top + 1)]
        for q in range(nl):
            L, I, J = leaves[q]
            for m in range(top, L, -1):
                s = 1 << (m - L)
                if need[m][I * s:(I + 1) * s, J * s:(J + 1) * s].any():
                    self.leaf_le[q] = m
                    break
            m = int(self.leaf_le[q])
            if m >= 1:
                s = 1 << (m - L)
                evalreq[m][I * s:(I + 1) * s, J * s:(J + 1) * s] = True
        D = [None] * (top + 2)
        for m in range(top, 0, -1):
            base = evalreq[m] | T[m]
            if m < top:
                base = base | _coarsen(D[m + 1])
            D[m] = _dilate(base, 1)
        # -- per-level structures --------------------------------------------------
        for m in range(1, top + 1):
            lv = _Level()
            lv.m, lv.N, lv.nb = m, f.ncell(m), f.nbl(m)
            dI, dJ = np.nonzero(D[m])
            lv.dkeys = dI.astype(np.int64) * lv.nb + dJ          # sorted (row-major nonzero)
            lv.nD = len(lv.dkeys) * B * B
            lv.tmap = np.full((lv.nb, lv.nb), -1, dtype=np.int64)  # block-resolution lookup (not a fine-grid array)
            lv.tmap[dI, dJ] = np.arange(len(dI))
            self.lev[m] = lv
        r = np.arange(B, dtype=np.int64)
        for m in range(1, top + 1):
            lv = self.lev[m]
            N = lv.N
            tI, tJ = lv.dkeys // lv.nb, lv.dkeys % lv.nb
            gi = (tI[:, None, None] * B + r[None, :, None]) + 0 * r[None, None, :]
            gj = (tJ[:, None, None] * B + r[None, None, :]) + 0 * r[None, :, None]
            gi, gj = gi.ravel(), gj.ravel()                       # global indices of all D points
            # subdivision P_m : a_{m-1} -> level-m coefficients on D_m (CSR assembled row by row)
            mask = two_scale_mask(k)
            hw = (k + 1) // 2
            nst = hw + 1 if k == 5 else hw + 1                    # max coarse taps per dimension (4 for k=5, 3 for k=3)
            offs = np.zeros((2, nst), dtype=np.int64)
            wts = np.zeros((2, nst))
            for par in (0, 1):
                q = 0
                for o in range(-hw, hw + 1):                      # o = j - 2 i
                    if (o - par) % 2 == 0:
                        offs[par, q] = (par - o) // 2             # coarse index = q + (par - o)/2
                        wts[par, q] = mask[o + hw]
                        q += 1
            ncoarse = f.n0 * f.n0 if m == 1 else self.lev[m - 1].nD
            CH = 1 << 15                                          # rows per chunk: bounds the transient memory
            dat, col, cnt = [], [], []
            for r0 in range(0, lv.nD, CH):
                ci_, cj_ = gi[r0:r0 + CH], gj[r0:r0 + CH]
                n = len(ci_)
                pi, pj = ci_ & 1, cj_ & 1
                ci = (ci_ >> 1)[:, None] + offs[pi]                # (n, nst)
                cj = (cj_ >> 1)[:, None] + offs[pj]
                vals = wts[pi][:, :, None] * wts[pj][:, None, :]
                cols = np.empty((n, nst, nst), dtype=np.int32)
                for qa in range(nst):
                    fl, ok = self._flat(m - 1, np.broadcast_to(ci[:, qa, None], (n, nst)), cj)
                    if not ok.all():
                        raise RuntimeError("layout error: coarse support missing for subdivision")
                    cols[:, qa, :] = fl
                keep = vals != 0.0
                dat.append(vals[keep])
                col.append(cols[keep])
                cnt.append(keep.reshape(n, -1).sum(axis=1))
            indptr = np.concatenate([[0], np.cumsum(np.concatenate(cnt))]).astype(np.int32)
            lv.P = sp.csr_matrix((np.concatenate(dat), np.concatenate(col), indptr), shape=(lv.nD, ncoarse))
            del dat, col, cnt
            # solve region T_m (subset of D_m points, in D flat order)
            inT = T[m][tI, tJ]
            tsel = np.nonzero(np.repeat(inT, B * B))[0]
            lv.tsel, lv.nT = tsel, len(tsel)
            tindex = np.full(lv.nD, -1, dtype=np.int32)
            tindex[tsel] = np.arange(lv.nT, dtype=np.int32)
            ti, tj = gi[tsel], gj[tsel]
            taps = COLLOC[k]
            w = np.arange(-h, h + 1, dtype=np.int64)
            tt = np.broadcast_to(taps[:, None] * taps[None, :], (1, k, k))

            def colloc_rows(ri, rj):
                """CSR rows of the collocation operator restricted to T_m for the
                level-m grid points (ri, rj) (columns outside T_m dropped)."""
                dat, col, cnt = [], [], []
                for r0 in range(0, len(ri), CH):
                    a_, b_ = ri[r0:r0 + CH], rj[r0:r0 + CH]
                    n = len(a_)
                    cols = np.empty((n, k, k), dtype=np.int32)
                    for a in range(k):
                        fl, ok = self._flat(m, np.broadcast_to((a_ + w[a])[:, None], (n, k)), b_[:, None] + w[None, :])
                        cols[:, a, :] = np.where(ok, tindex[fl], -1)
                    keep = cols >= 0
                    dat.append(np.broadcast_to(tt, (n, k, k))[keep])
                    col.append(cols[keep])
                    cnt.append(keep.reshape(n, -1).sum(axis=1))
                indptr = np.concatenate([[0], np.cumsum(np.concatenate(cnt))]).astype(np.int32)
                return sp.csr_matrix((np.concatenate(dat), np.concatenate(col), indptr), shape=(len(ri), lv.nT))

            # full-torus level: the solve is the periodic prefilter (2D FFT)
            lv.full = lv.nT == N * N
            lv.dense_of_T = (ti % N) * N + (tj % N) if lv.full else None
            # the full SPD matrix is only needed by the verification solvers / extra sweeps
            need_M = self.solver != "adi" or (self.sweeps > 1 and not lv.full)
            lv.M = colloc_rows(ti, tj) if need_M else None
            # active nodes of this level -> rows of the solve region
            keys = f.node_keys(m)
            lv.node_keys = keys
            fl, ok = self._flat(m, keys // N, keys % N)
            if not ok.all() or (tindex[fl] < 0).any():
                raise RuntimeError("layout error: active node outside its solve region")
            lv.node_rows = tindex[fl]
            lv.Mn = colloc_rows(keys // N, keys % N)
            lv.lu = None
            lv.perm, lv.chol, lv.crow, lv.lam = [], [], [], None
            if lv.full and self.solver != "splu":
                continue
            # line (ADI) operators: exact inverses of the 1D collocation matrix restricted
            # to each grid line of T_m (banded Cholesky; complete periodic lines by FFT)
            lv.perm, lv.chol, lv.crow = [], [], []
            ti_w, tj_w = ti % N, tj % N
            q = np.arange(N // 2 + 1)
            lam = np.full(N // 2 + 1, taps[h])
            for d in range(1, h + 1):
                lam = lam + 2.0 * taps[h + d] * np.cos(2.0 * np.pi * d * q / N)
            lv.lam = lam
            for axis in (0, 1):                                   # lines along i (fixed j), then along j
                pos, line = (ti_w, tj_w) if axis == 0 else (tj_w, ti_w)
                fl, ok = self._flat(m, ti - (1 if axis == 0 else 0), tj - (1 if axis == 1 else 0))
                haspred = ok & (tindex[np.where(ok, fl, 0)] >= 0)
                origin = np.full(N, N, dtype=np.int64)
                np.minimum.at(origin, line[~haspred], pos[~haspred])
                complete = origin == N                            # no segment start: full periodic line
                origin[complete] = 0
                posr = (pos - origin[line]) % N                   # rotate so that no segment straddles the seam
                order = np.lexsort((posr, line))
                ls, ps = line[order], posr[order]
                ab = np.zeros((h + 1, lv.nT))
                ab[0] = taps[h]
                for d in range(1, h + 1):
                    same = (ls[d:] == ls[:-d]) & (ps[d:] - ps[:-d] == d)
                    ab[d, :-d] = np.where(same, taps[h + d], 0.0)
                lv.perm.append(order)
                lv.chol.append(sla.cholesky_banded(ab, lower=True))
                lv.crow.append(np.nonzero(complete[ls])[0])
        # -- evaluation windows per (L, d) group -------------------------------------
        self.groups = {}
        for (L, d) in sorted(set(zip(self.leafL.tolist(), (self.leaf_le - self.leafL).tolist()))):
            ids = np.nonzero((self.leafL == L) & (self.leaf_le - self.leafL == d))[0]
            m = L + d
            nc = B << d
            w = np.arange(-h, nc + h + 1, dtype=np.int64)           # nc + k window points
            wi = (self.leafI[ids] * nc)[:, None] + w[None, :]
            wj = (self.leafJ[ids] * nc)[:, None] + w[None, :]
            fl, ok = self._flat(m, np.broadcast_to(wi[:, :, None], (len(ids), len(w), len(w))),
                                np.broadcast_to(wj[:, None, :], (len(ids), len(w), len(w))))
            if not np.all(ok):
                raise RuntimeError("layout error: evaluation window outside coefficient tiles")
            self.groups[(L, d)] = (ids, fl.astype(np.int32))
        # -- sorted block tables for point location ------------------------------------
        self._exist_keys = []
        self._leaf_tab = []
        for l in range(top + 1):
            nb = f.nbl(l)
            ek = np.array(sorted(I * nb + J for (I, J) in f.exists[l]), dtype=np.int64)
            self._exist_keys.append(ek)
            sel = self.leafL == l
            lk = self.leafI[sel] * nb + self.leafJ[sel]
            o = np.argsort(lk)
            self._leaf_tab.append((lk[o], self.leaf_le[sel][o]))
        self.version = f.version
        self.stats["layout_s"] += time.perf_counter() - t0
        self.stats["n_layout"] += 1

    # ------------------------------------------------------------------ memory
    def memory_floats(self, ncomp):
        """Floats held for `ncomp` coefficient components plus the topology
        operators (sparse matrices counted as 1.5 floats per non-zero)."""
        n = self.f.n0 ** 2 * ncomp
        op = 0.0
        for m in range(1, self.top + 1):
            lv = self.lev[m]
            n += lv.nD * ncomp
            op += 1.5 * (lv.P.nnz + (lv.M.nnz if lv.M is not None else 0) + lv.Mn.nnz)
            op += (2 * (self.h + 1) * lv.nT + 2 * lv.nT) if lv.chol else 0
        win = sum(v[1].size for v in self.groups.values())
        return dict(coefficients=float(n), operators=float(op), windows=float(win))

    # ------------------------------------------------------------------ solve
    def _line_solve(self, lv, q, y):
        """Exact solve with the 1D collocation operator along direction q on T_m."""
        o = lv.perm[q]
        t = y[o]
        s = sla.cho_solve_banded((lv.chol[q], True), t, check_finite=False)
        cr = lv.crow[q]
        if cr.size:
            Y = t[cr].reshape(-1, lv.N, t.shape[1])
            s[cr] = np.fft.irfft(np.fft.rfft(Y, axis=1) / lv.lam[None, :, None], n=lv.N, axis=1).reshape(-1, t.shape[1])
        out = np.empty_like(y)
        out[o] = s
        return out

    def _cg(self, lv, r, c, scale):
        """Plain conjugate gradients on the SPD system (guaranteed fallback)."""
        res = r - lv.M @ c
        p = res.copy()
        rs = np.sum(res * res, axis=0)
        it = 0
        while float(np.abs(res).max()) > self.tol * scale and it < 2000:
            Mp = lv.M @ p
            den = np.sum(p * Mp, axis=0)
            al = np.where(den > 0, rs / np.where(den > 0, den, 1.0), 0.0)
            c = c + al * p
            res = res - al * Mp
            rs2 = np.sum(res * res, axis=0)
            p = res + np.where(rs > 0, rs2 / np.where(rs > 0, rs, 1.0), 0.0) * p
            rs = rs2
            it += 1
        return c, it

    def _solve(self, lv, r, tag):
        """Solve M c = r on T_m.

        The product of the two exact line inverses is the exact inverse of M
        whenever T_m is a rectangle / band / the full torus; otherwise it
        differs from it only through couplings across re-entrant corners of
        T_m, which lie >= halo_tiles*B points away from every non-zero of r,
        where c has already decayed by lambda^(halo_tiles*B).  "adi" applies a
        fixed number of alternating sweeps (default 1: the interpolation defect
        at the nodes stays at the skirt-leak level, measured in tests);
        "adi-tol" iterates to `tol` with a CG fallback.  A level whose solve
        region is the whole torus uses the periodic FFT prefilter (exact)."""
        if self.solver == "splu":
            if lv.lu is None:
                lv.lu = spla.splu(lv.M.tocsc(), permc_spec="MMD_AT_PLUS_A")
            return lv.lu.solve(r)
        if lv.full and self.solver != "splu":
            N = lv.N
            dense = np.empty((N * N, r.shape[1]))
            dense[lv.dense_of_T] = r
            return prefilter_periodic(dense.reshape(N, N, -1), self.k).reshape(N * N, -1)[lv.dense_of_T]
        if self.solver == "adi":
            c = self._line_solve(lv, 1, self._line_solve(lv, 0, r))
            for it in range(1, self.sweeps):
                first, second = (0, 1) if it % 2 == 0 else (1, 0)
                c = c + self._line_solve(lv, second, self._line_solve(lv, first, r - lv.M @ c))
            return c
        scale = float(np.abs(r).max())
        if scale == 0.0:
            self.stats["solve_iters"].append(0)
            return np.zeros_like(r)
        c = np.zeros_like(r)
        res = r
        it = 0
        err = float(np.abs(res).max())
        while err > self.tol * scale:
            if it >= self.max_iter:
                c, n_cg = self._cg(lv, r, c, scale)
                self.stats["cg_fallbacks"] = self.stats.get("cg_fallbacks", 0) + 1
                it += n_cg
                break
            first, second = (0, 1) if it % 2 == 0 else (1, 0)
            c = c + self._line_solve(lv, second, self._line_solve(lv, first, res))
            res = r - lv.M @ c
            new = float(np.abs(res).max())
            it += 1
            if new > 0.9 * err and it >= 4:      # not contracting: hand over to CG
                c, n_cg = self._cg(lv, r, c, scale)
                self.stats["cg_fallbacks"] = self.stats.get("cg_fallbacks", 0) + 1
                it += n_cg
                break
            err = new
        self.stats["solve_iters"].append(it)
        return c

    # ------------------------------------------------------------------ update
    def update(self, zs, tag="x", return_surplus=False):
        """zs[m]: nodal data (n_m, ncomp) of level m in `forest.node_keys(m)`
        order (periodic fields only: displacement, momentum, periodic action).
        Returns the list of accumulated coefficient arrays a[m] (flat, ncomp)."""
        if self.version != self.f.version:
            self.build_layout()
        t0 = time.perf_counter()
        n0 = self.f.n0
        z0 = np.asarray(zs[0], dtype=float)
        nc = z0.shape[1]
        a = [prefilter_periodic(z0.reshape(n0, n0, nc), self.k).reshape(n0 * n0, nc)]
        surplus = [None]
        for m in range(1, self.top + 1):
            lv = self.lev[m]
            b = lv.P @ a[m - 1]
            rn = np.asarray(zs[m], dtype=float) - lv.Mn @ b[lv.tsel]
            r = np.zeros((lv.nT, nc))
            r[lv.node_rows] = rn
            c = self._solve(lv, r, tag)
            b[lv.tsel] += c
            a.append(b)
            surplus.append(rn)
        self.stats["update_s"] += time.perf_counter() - t0
        return (a, surplus) if return_surplus else a

    # ------------------------------------------------------------------ evaluation
    def weights(self, d, s, deriv=False):
        """Tensor evaluation matrix for a block evaluated d levels finer with
        s midpoint samples per (own-level) cell: (B*s, B*2^d + k)."""
        key = (d, s, deriv)
        W = self._wcache.get(key)
        if W is None:
            B, k = self.B, self.k
            nc = B << d
            u = (np.arange(B * s) + 0.5) / s * (1 << d)
            i0 = np.floor(u).astype(np.int64)
            _, ws = eval_weights(u - i0, k, 1 if deriv else 0)
            W = np.zeros((B * s, nc + k))
            rr = np.arange(B * s)
            for q in range(k + 1):
                W[rr, i0 + q] = ws[1 if deriv else 0][q]
            self._wcache[key] = W
        return W

    def node_weights(self, d):
        """Evaluation matrix at the B+1 own-level grid lines of a block (closed)."""
        key = (d, "nodes")
        W = self._wcache.get(key)
        if W is None:
            B, k = self.B, self.k
            nc = B << d
            taps = COLLOC[k]
            W = np.zeros((B + 1, nc + k))
            for c in range(B + 1):
                W[c, (c << d): (c << d) + k] = taps      # value at integer knot: taps on i-h..i+h
            self._wcache[key] = W
        return W

    def block_values(self, a, L, d, W1, W2, sel=None):
        """Tensor evaluation for the (L, d) group: returns (ids, V) with
        V[nblk, n1, n2, ncomp] = W1 . A . W2^T on each block window."""
        ids, idx = self.groups[(L, d)]
        if sel is not None:
            ids, idx = ids[sel], idx[sel]
        A = a[L + d][idx]                                          # (nblk, mw, mw, nc)
        t1 = np.tensordot(A, W1, axes=([1], [1]))                  # (nblk, mw, nc, n1)
        t2 = np.tensordot(t1, W2, axes=([1], [1]))                 # (nblk, nc, n1, n2)
        return ids, np.moveaxis(t2, 1, 3)

    def block_values_cf(self, a, L, d, W1, W2, sel=None, ncomp=None):
        """As block_values but components-first and C-contiguous:
        V[nblk, ncomp, n1, n2] (the layout consumed by the deposit kernel)."""
        ids, idx = self.groups[(L, d)]
        if sel is not None:
            ids, idx = ids[sel], idx[sel]
        A = a[L + d][idx] if ncomp is None else a[L + d][idx, :ncomp]
        t1 = np.tensordot(A, W1, axes=([1], [1]))                  # (nblk, mw, nc, n1)
        return ids, np.tensordot(t1, W2, axes=([1], [1]))          # (nblk, nc, n1, n2)

    def locate(self, xi):
        """Leaf block level and its evaluation level for points xi (M, 2)."""
        f = self.f
        xi = np.asarray(xi, dtype=float)
        w = xi - np.floor(xi)
        lev = np.full(len(xi), -1, dtype=np.int64)
        le = np.zeros(len(xi), dtype=np.int64)
        for l in range(self.top, -1, -1):
            todo = lev < 0
            if not todo.any():
                break
            nb = f.nbl(l)
            bi = np.minimum((w[todo] * nb).astype(np.int64), nb - 1)
            key = bi[:, 0] * nb + bi[:, 1]
            ek = self._exist_keys[l]
            pos = np.minimum(np.searchsorted(ek, key), len(ek) - 1)
            hit = ek[pos] == key
            lk, lle = self._leaf_tab[l]
            idx = np.nonzero(todo)[0][hit]
            lev[idx] = l
            p2 = np.searchsorted(lk, key[hit])
            le[idx] = lle[p2]
        return lev, le

    def eval_points(self, a, xi, deriv=0):
        """F-components at arbitrary reference points xi (M, 2): values (M, nc)
        and, if deriv >= 1, the xi-gradient (M, nc, 2)."""
        xi = np.asarray(xi, dtype=float)
        M, nc = len(xi), a[0].shape[1]
        val = np.zeros((M, nc))
        g = np.zeros((M, nc, 2)) if deriv >= 1 else None
        _, le = self.locate(xi)
        k = self.k
        for m in np.unique(le):
            sel = np.nonzero(le == m)[0]
            N = self.f.ncell(int(m))
            u = xi[sel] * N
            i0 = np.floor(u).astype(np.int64)
            t = u - i0
            js, w1 = eval_weights(t[:, 0], k, deriv)
            _, w2 = eval_weights(t[:, 1], k, deriv)
            am = a[int(m)]
            v = np.zeros((len(sel), nc))
            g0 = np.zeros((len(sel), nc)) if deriv >= 1 else None
            g1 = np.zeros((len(sel), nc)) if deriv >= 1 else None
            for qa in range(k + 1):
                for qb in range(k + 1):
                    fl, ok = self._flat(int(m), i0[:, 0] + js[qa], i0[:, 1] + js[qb])
                    if not ok.all():
                        raise RuntimeError("eval_points: stencil outside coefficient tiles")
                    cc = am[fl]
                    v += (w1[0][qa] * w2[0][qb])[:, None] * cc
                    if deriv >= 1:
                        g0 += (w1[1][qa] * w2[0][qb])[:, None] * cc
                        g1 += (w1[0][qa] * w2[1][qb])[:, None] * cc
            val[sel] = v
            if deriv >= 1:
                g[sel, :, 0] = g0 * N
                g[sel, :, 1] = g1 * N
        return val if deriv == 0 else (val, g)

    def decay(self):
        return DECAY[self.k]
