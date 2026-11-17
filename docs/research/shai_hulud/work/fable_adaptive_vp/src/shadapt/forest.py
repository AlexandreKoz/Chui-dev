"""Label-space topology: a periodic quadtree of blocks.

Reference (label) torus [0,1)^2.  Level l has N_l = n0 * 2^l cells per
dimension and nb_l = N_l / B blocks per dimension; block (l, I, J) covers the
level-l cells [I B, (I+1) B) x [J B, (J+1) B).  All level-0 blocks exist.  A
block may be *refined*: its four children (l+1, 2I+a, 2J+b) then exist.  A
block that exists and is not refined is a *leaf*; the leaf blocks tile the
torus exactly once and their B x B cells are the leaf cells that carry mass.

Invariants (checked by `check()`):
  nesting   a block at level l >= 1 exists iff its parent exists and is refined
  2:1       every 8-neighbour (edges and corners, periodic) of a leaf at level
            l is covered by blocks of level >= l-1  (equivalently: all 8
            neighbours of a refined block exist at the same level)
  tiling    sum of leaf areas == 1 exactly (dyadic arithmetic)

Nodes.  The grid points of level l are (i, j) / N_l.  The *active nodes of
level l* are the level-l grid points that lie in the CLOSED region covered by
existing level-l blocks and are not grid points of level l-1 (at least one
odd index); level 0 owns all its grid points.  Hence every active node has
exactly one owner level and one integer key  i * N_l + j  there (deterministic
ownership; nodes on shared block edges are stored once; a node on a
coarse-fine interface -- a "hanging" node -- is an ordinary active node).
"""
from __future__ import annotations

import numpy as np


def _is_pow2(v):
    return v >= 1 and (v & (v - 1)) == 0


class Forest:
    def __init__(self, n0, B=8, lmax=3):
        if not (_is_pow2(n0) and _is_pow2(B) and n0 % B == 0 and n0 // B >= 1):
            raise ValueError("n0 and B must be powers of two with B | n0")
        self.n0, self.B, self.lmax = int(n0), int(B), int(lmax)
        self.nb0 = self.n0 // self.B
        self.exists = [set() for _ in range(self.lmax + 1)]
        self.refined = [set() for _ in range(self.lmax + 1)]
        self.exists[0] = {(i, j) for i in range(self.nb0) for j in range(self.nb0)}
        self.lineage = {}   # (l, I, J) -> dict(parent, step, t, reason, ...)
        self.events = []    # chronological record of every split / merge
        self.version = 0    # bumped on every topology change

    # -- sizes ----------------------------------------------------------------
    def nbl(self, l):
        return self.nb0 << l

    def ncell(self, l):
        return self.n0 << l

    def is_leaf(self, l, I, J):
        return (I, J) in self.exists[l] and (I, J) not in self.refined[l]

    def leaves(self):
        """Sorted list of leaf blocks (l, I, J): deterministic order."""
        out = []
        for l in range(self.lmax + 1):
            out.extend((l, I, J) for (I, J) in sorted(self.exists[l] - self.refined[l]))
        return out

    def top_level(self):
        return max(l for l in range(self.lmax + 1) if self.exists[l])

    def neighbours(self, l, I, J):
        nb = self.nbl(l)
        out = []
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                if a or b:
                    q = ((I + a) % nb, (J + b) % nb)
                    if q != (I, J) and q not in out:
                        out.append(q)
        return out

    # -- refinement -------------------------------------------------------------
    def refine(self, l, I, J, reason="indicator", step=-1, t=float("nan"), info=None, _depth=0):
        """Split the leaf block (l, I, J).  Enforces 2:1 balance by first
        splitting coarser neighbours (recorded with reason "balance").
        Returns the list of blocks split by this call, in order."""
        if not self.is_leaf(l, I, J):
            raise ValueError(f"block {(l, I, J)} is not a leaf")
        if l >= self.lmax:
            raise ValueError("level cap reached")
        done = []
        for (Ia, Ja) in self.neighbours(l, I, J):
            if (Ia, Ja) not in self.exists[l]:
                # by the 2:1 invariant the covering block is the leaf parent
                par = (l - 1, Ia // 2, Ja // 2)
                if self.is_leaf(*par):
                    done += self.refine(*par, reason="balance", step=step, t=t,
                                        info=dict(caused_by=[l, I, J]), _depth=_depth + 1)
        self.refined[l].add((I, J))
        for a in (0, 1):
            for b in (0, 1):
                c = (2 * I + a, 2 * J + b)
                self.exists[l + 1].add(c)
                self.lineage[(l + 1,) + c] = dict(parent=[l, I, J], step=step, t=t, reason=reason)
        ev = dict(kind="split", level=l, I=I, J=J, reason=reason, step=step, t=t)
        if info:
            ev.update(info)
        self.events.append(ev)
        self.version += 1
        done.append((l, I, J))
        return done

    def can_coarsen(self, l, I, J):
        """True if the refined block (l, I, J) may be merged back into a leaf."""
        if (I, J) not in self.refined[l]:
            return False
        kids = [(2 * I + a, 2 * J + b) for a in (0, 1) for b in (0, 1)]
        if any(k in self.refined[l + 1] for k in kids):
            return False
        for k in kids:  # no neighbour of a child may be refined (would be 2 levels finer)
            for q in self.neighbours(l + 1, *k):
                if q in self.refined[l + 1]:
                    return False
        return True

    def coarsen(self, l, I, J, reason="indicator", step=-1, t=float("nan"), info=None):
        if not self.can_coarsen(l, I, J):
            raise ValueError(f"block {(l, I, J)} cannot be coarsened")
        for a in (0, 1):
            for b in (0, 1):
                c = (2 * I + a, 2 * J + b)
                self.exists[l + 1].discard(c)
                self.lineage.pop((l + 1,) + c, None)
        self.refined[l].discard((I, J))
        ev = dict(kind="merge", level=l, I=I, J=J, reason=reason, step=step, t=t)
        if info:
            ev.update(info)
        self.events.append(ev)
        self.version += 1

    # -- invariants -----------------------------------------------------------
    def check(self):
        """Raise AssertionError if any topological invariant is violated."""
        nb0 = self.nb0
        assert self.exists[0] == {(i, j) for i in range(nb0) for j in range(nb0)}
        for l in range(self.lmax + 1):
            assert self.refined[l] <= self.exists[l]
            if l == self.lmax:
                assert not self.refined[l]
            kids = {(2 * I + a, 2 * J + b) for (I, J) in self.refined[l] for a in (0, 1) for b in (0, 1)}
            if l < self.lmax:
                assert kids == self.exists[l + 1], "nesting violated"
            for (I, J) in self.refined[l]:
                for q in self.neighbours(l, I, J):
                    assert q in self.exists[l], f"2:1 balance violated at {(l, I, J)}"
        # exact tiling: leaf areas in units of the finest block area
        top = self.lmax
        tot = sum(len(self.exists[l] - self.refined[l]) * 4 ** (top - l) for l in range(self.lmax + 1))
        assert tot == (self.nb0 << top) ** 2, "leaf blocks do not tile the torus exactly once"
        return True

    # -- cells and mass -------------------------------------------------------
    def cell_mass(self, l):
        """Label-space area (= mass, total 1) of one leaf cell at level l; exact."""
        return 1.0 / float(self.ncell(l)) ** 2

    def n_leaf_cells(self):
        return sum(len(self.exists[l] - self.refined[l]) for l in range(self.lmax + 1)) * self.B ** 2

    def leaf_level_map(self):
        """Level of the leaf block covering each finest-level block slot
        (dense array of size nbl(top)^2; diagnostic / test helper only)."""
        top = self.top_level()
        nb = self.nbl(top)
        out = np.full((nb, nb), -1, dtype=np.int8)
        cnt = np.zeros((nb, nb), dtype=np.int16)
        for (l, I, J) in self.leaves():
            s = 1 << (top - l)
            out[I * s:(I + 1) * s, J * s:(J + 1) * s] = l
            cnt[I * s:(I + 1) * s, J * s:(J + 1) * s] += 1
        assert np.all(cnt == 1), "leaf ownership is not unique"
        return out

    # -- nodes ----------------------------------------------------------------
    def node_keys(self, l):
        """Sorted unique keys i * N_l + j of the active nodes owned by level l."""
        N = self.ncell(l)
        if l == 0:
            return np.arange(N * N, dtype=np.int64)
        if not self.exists[l]:
            return np.zeros(0, dtype=np.int64)
        B = self.B
        blk = np.array(sorted(self.exists[l]), dtype=np.int64)
        r = np.arange(B + 1, dtype=np.int64)
        i = (blk[:, 0, None] * B + r[None, :]) % N            # (nblk, B+1)
        j = (blk[:, 1, None] * B + r[None, :]) % N
        ii = np.broadcast_to(i[:, :, None], (len(blk), B + 1, B + 1))
        jj = np.broadcast_to(j[:, None, :], (len(blk), B + 1, B + 1))
        odd = ((ii | jj) & 1) == 1
        return np.unique(ii[odd] * N + jj[odd])

    def n_nodes(self):
        return sum(len(self.node_keys(l)) for l in range(self.lmax + 1))

    # -- serialisation --------------------------------------------------------
    def to_dict(self):
        return dict(n0=self.n0, B=self.B, lmax=self.lmax, version=self.version,
                    refined=[sorted(list(map(list, s))) for s in self.refined],
                    lineage=[[list(k), v] for k, v in sorted(self.lineage.items())], events=self.events)

    @classmethod
    def from_dict(cls, d):
        f = cls(d["n0"], d["B"], d["lmax"])
        for l, blocks in enumerate(d["refined"]):
            for (I, J) in blocks:
                f.refined[l].add((I, J))
                for a in (0, 1):
                    for b in (0, 1):
                        f.exists[l + 1].add((2 * I + a, 2 * J + b))
        f.lineage = {tuple(k): v for k, v in d["lineage"]}
        f.events = list(d["events"])
        f.version = d["version"]
        f.check()
        return f

    def signature(self):
        """Hashable, order-independent description of the topology."""
        return tuple(tuple(sorted(s)) for s in self.refined)
