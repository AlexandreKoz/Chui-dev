"""Adaptive sheet state: nodal (displacement, momentum, action) on the active
nodes of a block quadtree, the composite carrier built from them, the leaf
quadrature, and the transfer operators (split / merge).

Conventions are the baseline ones (`shsheet.conventions`): specific momentum
p = dx/dt, theta = p.dx, positions x(xi) = xi + d(xi) with d periodic (degree
one map of the unit torus), action S(xi + e_a) = S(xi) + c_a.  The mass
measure is the Lebesgue measure of the LABEL torus (total 1): a leaf cell of
level l carries mass 4^-l / n0^2, whatever the state.  Mass is a property of
label cells, never of nodes, and is not touched by any transfer.

State-only rule: nothing in this module sees a reference solution.
"""
from __future__ import annotations

import numpy as np

from .forest import Forest
from .hcarrier import HierarchicalCarrier


class AdaptiveSheet:
    def __init__(self, forest, k=5, periods=(0.0, 0.0), halo_tiles=1, solver="adi", sweeps=1):
        self.f = forest
        self.k = k
        self.periods = np.asarray(periods, dtype=float)
        self.car = HierarchicalCarrier(forest, k=k, halo_tiles=halo_tiles, solver=solver, sweeps=sweeps)
        self.keys, self.xi, self.d, self.p, self.S = [], [], [], [], []
        self.transfer_log = []
        self._sync_nodes_new()

    # ------------------------------------------------------------------ nodes
    def _xi_of(self, l, keys):
        N = self.f.ncell(l)
        return np.stack([keys // N, keys % N], axis=1) / float(N)

    def _sync_nodes_new(self):
        self.keys = [self.f.node_keys(l) for l in range(self.f.top_level() + 1)]
        self.xi = [self._xi_of(l, k) for l, k in enumerate(self.keys)]
        self.d = [np.zeros((len(k), 2)) for k in self.keys]
        self.p = [np.zeros((len(k), 2)) for k in self.keys]
        self.S = [np.zeros(len(k)) for k in self.keys]

    @classmethod
    def from_ic(cls, ic, n0, B=8, lmax=3, k=5, **kw):
        sh = cls(Forest(n0, B, lmax), k=k, periods=getattr(ic, "periods", (0.0, 0.0)), **kw)
        sh.set_from_ic(ic)
        return sh

    def set_from_ic(self, ic):
        for l in range(len(self.keys)):
            x, p, S = ic(self.xi[l])
            self.d[l], self.p[l], self.S[l] = x - self.xi[l], p.copy(), S.copy()

    @property
    def top(self):
        return len(self.keys) - 1

    def n_nodes(self):
        return int(sum(len(k) for k in self.keys))

    def nodes_per_level(self):
        return [int(len(k)) for k in self.keys]

    def x_nodes(self):
        """Unwrapped positions of all active nodes, levels concatenated."""
        return np.concatenate([self.xi[l] + self.d[l] for l in range(len(self.keys))], axis=0)

    def level_slices(self):
        out, o = [], 0
        for k in self.keys:
            out.append(slice(o, o + len(k)))
            o += len(k)
        return out

    # ------------------------------------------------------------------ carrier
    def nodal_fields(self, which):
        """Per-level nodal arrays of periodic fields.  which: sequence drawn
        from 'x' (displacement, 2), 'p' (2), 'S' (periodic action, 1)."""
        out = []
        for l in range(len(self.keys)):
            cols = []
            for w in which:
                if w == "x":
                    cols.append(self.d[l])
                elif w == "p":
                    cols.append(self.p[l])
                elif w == "S":
                    cols.append((self.S[l] - self.xi[l] @ self.periods)[:, None])
                else:
                    raise ValueError(w)
            out.append(np.concatenate(cols, axis=1))
        return out

    def coeffs(self, which="x", return_surplus=False):
        return self.car.update(self.nodal_fields(which), tag=which, return_surplus=return_surplus)

    def evaluate(self, xi, which="xp", deriv=0):
        """Carrier at reference points xi (M, 2).  Returns arrays in the order
        of `which`; 'x' is returned as the unwrapped position xi + d(xi) and
        'S' as the full (non-periodic) action."""
        a = self.coeffs(which)
        r = self.car.eval_points(a, xi, deriv)
        v, g = (r, None) if deriv == 0 else r
        out, o = [], 0
        for w in which:
            n = 1 if w == "S" else 2
            val = v[:, o:o + n]
            if w == "x":
                val = val + xi
            elif w == "S":
                val = val[:, 0] + xi @ self.periods
            item = val
            if deriv >= 1:
                gg = g[:, o:o + n].copy()
                if w == "x":
                    gg[:, 0, 0] += 1.0
                    gg[:, 1, 1] += 1.0
                elif w == "S":
                    gg = gg[:, 0] + self.periods
                item = (val, gg)
            out.append(item)
            o += n
        return out

    # ------------------------------------------------------------------ transfer
    def _retopologize(self, a_old, old_car_eval):
        """After a topology change: keep old nodes bit-for-bit, read new nodes
        from the OLD composite carrier (never from a reference)."""
        f = self.f
        new_keys = [f.node_keys(l) for l in range(f.top_level() + 1)]
        nd, npp, nS, nxi = [], [], [], []
        created, kept = 0, 0
        for l, nk in enumerate(new_keys):
            xi = self._xi_of(l, nk)
            d = np.empty((len(nk), 2))
            p = np.empty((len(nk), 2))
            S = np.empty(len(nk))
            have = np.zeros(len(nk), dtype=bool)
            if l < len(self.keys) and len(self.keys[l]):
                ok = self.keys[l]
                pos = np.minimum(np.searchsorted(ok, nk), len(ok) - 1)
                have = ok[pos] == nk
                d[have], p[have], S[have] = self.d[l][pos[have]], self.p[l][pos[have]], self.S[l][pos[have]]
            new = ~have
            if new.any():
                v = old_car_eval(a_old, xi[new])
                d[new], p[new] = v[:, 0:2], v[:, 2:4]
                S[new] = v[:, 4] + xi[new] @ self.periods
            created += int(new.sum())
            kept += int(have.sum())
            nd.append(d), npp.append(p), nS.append(S), nxi.append(xi)
        removed = self.n_nodes() - kept
        self.keys, self.xi, self.d, self.p, self.S = new_keys, nxi, nd, npp, nS
        self.car.build_layout()
        return created, removed

    def apply_topology_change(self, change):
        """Run `change(forest)` (which may split and merge blocks) and transfer
        the state.  Returns dict(created, removed)."""
        a_old = self.coeffs("xpS")
        v0 = self.f.version
        change(self.f)
        if self.f.version == v0:
            return dict(created=0, removed=0)
        created, removed = self._retopologize(a_old, self.car.eval_points)
        rec = dict(created=created, removed=removed, nodes=self.n_nodes(), version=self.f.version)
        self.transfer_log.append(rec)
        return rec

    def refine(self, blocks, reason="manual", step=-1, t=float("nan"), info=None):
        """Split the listed leaf blocks [(l, I, J), ...] (2:1 balance enforced)."""
        def change(f):
            for (l, I, J) in blocks:
                if f.is_leaf(l, I, J) and l < f.lmax:
                    f.refine(l, I, J, reason=reason, step=step, t=t, info=(info or {}).get((l, I, J)))
        return self.apply_topology_change(change)

    def coarsen(self, blocks, reason="manual", step=-1, t=float("nan"), info=None):
        """Merge the children of the listed refined blocks where permitted."""
        def change(f):
            for (l, I, J) in blocks:
                if f.can_coarsen(l, I, J):
                    f.coarsen(l, I, J, reason=reason, step=step, t=t, info=(info or {}).get((l, I, J)))
        return self.apply_topology_change(change)

    def rebase(self):
        """If every level-0 block is refined, relabel level 1 as the new level
        0 (n0 -> 2 n0).  Exact: the composite carrier is unchanged (uniform
        limit), only bookkeeping changes.  Returns True if rebased."""
        f = self.f
        if f.lmax < 1 or len(f.refined[0]) != len(f.exists[0]):
            return False
        nf = Forest(2 * f.n0, f.B, f.lmax - 1)
        # blocks of old level l+1 are blocks of new level l
        for l in range(1, f.lmax):
            for (I, J) in sorted(f.refined[l]):
                nf.refined[l - 1].add((I, J))
                for a in (0, 1):
                    for b in (0, 1):
                        nf.exists[l].add((2 * I + a, 2 * J + b))
        nf.lineage = {(k[0] - 1, k[1], k[2]): v for k, v in f.lineage.items() if k[0] >= 2}
        nf.events = f.events + [dict(kind="rebase", n0=nf.n0)]
        nf.version = f.version + 1
        nf.check()
        # merge node arrays of old levels 0 and 1 into the new level 0
        N1 = f.ncell(1)
        k0 = self.keys[0]
        k0n = (k0 // f.n0) * 2 * N1 + (k0 % f.n0) * 2
        keys = np.concatenate([k0n, self.keys[1]])
        o = np.argsort(keys)
        assert len(keys) == N1 * N1 and np.array_equal(keys[o], np.arange(N1 * N1))
        cat = lambda arr: np.concatenate([arr[0], arr[1]], axis=0)[o]
        d0, p0, S0 = cat(self.d), cat(self.p), cat(self.S)
        rest = (self.d[2:], self.p[2:], self.S[2:])
        self.f = nf
        self.car = HierarchicalCarrier(nf, k=self.k, halo_tiles=self.car.ht, solver=self.car.solver, sweeps=self.car.sweeps)
        self.keys = [nf.node_keys(l) for l in range(nf.top_level() + 1)]
        self.xi = [self._xi_of(l, k) for l, k in enumerate(self.keys)]
        self.d, self.p, self.S = [d0] + list(rest[0]), [p0] + list(rest[1]), [S0] + list(rest[2])
        for l in range(len(self.keys)):
            assert len(self.d[l]) == len(self.keys[l])
        return True

    # ------------------------------------------------------------------ memory
    def state_floats(self):
        """Persistent floats of the evolving state: (d, p, S) = 5 per node,
        plus one int64 key per node."""
        return 6 * self.n_nodes()
