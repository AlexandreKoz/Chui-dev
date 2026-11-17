"""Topology: unique leaf ownership, exact tiling, 2:1 balance, periodic seams,
deterministic node ownership, lineage, serialisation."""
import numpy as np
import pytest

from helpers import random_forest
from shadapt.forest import Forest


def test_initial_forest_is_uniform():
    f = Forest(32, 8, 3)
    f.check()
    assert len(f.leaves()) == 16 and f.n_nodes() == 32 * 32 and f.n_leaf_cells() == 32 * 32


@pytest.mark.parametrize("seed", range(6))
def test_random_refinement_keeps_all_invariants(seed):
    f = random_forest(32, 8, 3, nsplit=40, seed=seed)
    f.check()                                    # nesting, 2:1 (edges AND corners), exact tiling
    lm = f.leaf_level_map()                      # asserts every slot is owned by exactly one leaf
    assert lm.min() >= 0
    # exact mass partition: leaf-cell masses sum to exactly 1 (dyadic)
    tot = sum(f.cell_mass(l) * f.B ** 2 for (l, _, _) in f.leaves())
    assert tot == 1.0
    # 2:1 on the dense map, including corners and the periodic seam
    for ax in (0, 1):
        for sh in (1, -1):
            assert np.abs(lm - np.roll(lm, sh, ax)).max() <= 1
    for s0 in (1, -1):
        for s1 in (1, -1):
            assert np.abs(lm - np.roll(np.roll(lm, s0, 0), s1, 1)).max() <= 1


def test_balance_cascade_across_the_periodic_seam():
    f = Forest(32, 8, 3)                         # 4 x 4 level-0 blocks
    f.refine(0, 0, 0)
    done = f.refine(1, 0, 0)                     # corner block: neighbours live across both seams
    f.check()
    split0 = {(e["I"], e["J"]) for e in f.events if e["level"] == 0}
    # the level-1 block (0,0) touches level-0 blocks (3,0), (0,3), (3,3) through the seam
    assert {(3, 0), (0, 3), (3, 3)} <= split0
    assert all(e["reason"] in ("indicator", "balance") for e in f.events)
    assert sum(e["reason"] == "balance" for e in f.events) == 3 and done[-1] == (1, 0, 0)


@pytest.mark.parametrize("seed", range(3))
def test_node_ownership_is_unique_and_deterministic(seed):
    f = random_forest(32, 8, 3, nsplit=30, seed=seed)
    seen = set()
    for l in range(f.top_level() + 1):
        k = f.node_keys(l)
        assert np.all(np.diff(k) > 0)            # sorted, no duplicate inside a level
        N = f.ncell(l)
        i, j = k // N, k % N
        if l > 0:
            assert np.all(((i | j) & 1) == 1)    # not a grid point of the coarser level
        top = f.top_level()
        s = 1 << (top - l)
        for a, b in zip((i * s).tolist(), (j * s).tolist()):   # label position on the finest lattice
            assert (a, b) not in seen            # no node owned by two levels
            seen.add((a, b))
    g = random_forest(32, 8, 3, nsplit=30, seed=seed)
    assert g.signature() == f.signature() and g.events == f.events
    assert all(np.array_equal(f.node_keys(l), g.node_keys(l)) for l in range(f.top_level() + 1))


def test_one_split_is_local():
    """A test that fails if 'local refinement' were a global regrid."""
    f = Forest(64, 8, 3)                         # 8 x 8 blocks, 4096 nodes
    n0, c0 = f.n_nodes(), f.n_leaf_cells()
    f.refine(0, 3, 3)
    B = f.B
    assert f.n_leaf_cells() - c0 == 3 * B * B                     # one block -> four
    assert f.n_nodes() - n0 == (2 * B + 1) ** 2 - (B + 1) ** 2      # new nodes of the closed patch only
    assert len(f.leaves()) == 64 - 1 + 4
    assert f.n_nodes() < 1.06 * n0                                  # nowhere near a 4x global regrid


def test_hanging_nodes_are_active_and_counted_once():
    f = Forest(32, 8, 2)
    f.refine(0, 1, 1)
    f.refine(0, 2, 1)                            # two adjacent fine blocks share an interior edge
    k = f.node_keys(1)
    N = f.ncell(1)
    i, j = k // N, k % N
    # nodes on the outer boundary (hanging) exist ...
    assert np.any((i == 16) & (j % 2 == 1)) and np.any((i == 48) & (j % 2 == 1))
    # ... and the shared edge i = 32 is stored once
    shared = (i == 32) & (j >= 16) & (j <= 32)
    assert shared.sum() == 8                     # odd j only (even j are level-0 nodes)
    assert len(k) == len(np.unique(k)) == (33 * 17) - (17 * 9)


def test_merge_restores_topology_and_lineage():
    f = Forest(32, 8, 3)
    sig0 = f.signature()
    f.refine(0, 1, 1, step=5, t=0.05)
    f.refine(1, 2, 2, step=7, t=0.07)
    assert f.lineage[(2, 4, 4)]["parent"] == [1, 2, 2] and f.lineage[(2, 4, 4)]["step"] == 7
    assert not f.can_coarsen(0, 1, 1)            # has a refined child
    balance_parents = [(e["level"], e["I"], e["J"]) for e in f.events if e["reason"] == "balance"]
    assert f.can_coarsen(1, 2, 2)
    f.coarsen(1, 2, 2)
    for (l, I, J) in reversed(balance_parents + [(0, 1, 1)]):
        if f.can_coarsen(l, I, J):
            f.coarsen(l, I, J)
    f.check()
    assert f.signature() == sig0 and not f.lineage


def test_merge_is_refused_when_it_would_break_balance():
    f = Forest(32, 8, 3)
    f.refine(0, 1, 1)
    f.refine(0, 2, 1)
    f.refine(1, 3, 2)                            # level-2 block adjacent to children of (0,2,1)
    assert not f.can_coarsen(0, 2, 1)
    with pytest.raises(ValueError):
        f.coarsen(0, 2, 1)


def test_serialisation_round_trip():
    f = random_forest(32, 8, 3, nsplit=20, seed=4)
    g = Forest.from_dict(f.to_dict())
    assert g.signature() == f.signature() and g.lineage == f.lineage and g.events == f.events
    assert all(np.array_equal(f.node_keys(l), g.node_keys(l)) for l in range(f.top_level() + 1))
