"""Quadrature study (PM 256, quintic n = 256): is a spacing-adaptive leaf
quadrature better than the uniform label lattice?  (It is not: recorded
negative result.)

  lattice / midpoint   uniform label lattice (the baseline rule)
  lattice / gauss2     two Gauss points per panel, uniform panels
  adaptive / midpoint  per-block sample counts from the block's image size
                       (target = largest sample spacing in PM cells)
  adaptive / gauss2    same counts, Gauss panels

Usage: python3 exp05_quadrature.py [NG]   -> results/exp05_quadrature_ng<NG>.json
(indicative timings: arms run concurrently)."""
from __future__ import annotations

import sys

import bench

if __name__ == "__main__":
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    A = {}
    for s in (2, 4):
        for rule in ("midpoint", "gauss2"):
            A[f"lattice-{rule}/s{s}"] = bench.sheet(ng, ng, quad=dict(mode="lattice", sub0=s, rule=rule))
    for tg in (1.0, 0.7, 0.5, 0.35):
        for rule in ("midpoint", "gauss2"):
            A[f"adaptive-{rule}/target{tg}"] = bench.sheet(ng, ng, quad=dict(mode="adaptive", sub0=2, target=tg, nq_min=8, nq_max=256, rule=rule))
    for tg in (1.0, 0.7):
        A[f"adaptive-midpoint/target{tg} lumped(nq_min=2)"] = bench.sheet(ng, ng, quad=dict(mode="adaptive", sub0=2, target=tg, nq_min=2, nq_max=256))
    bench.run_arms(f"exp05_quadrature_ng{ng}.json", A, jobs=3, meta_extra=dict(ng=ng), script=__file__)
