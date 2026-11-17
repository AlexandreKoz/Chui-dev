"""Gate-B benchmark: fixed cubic/quintic sheets, uniform early refinement,
local h-refinement and N-body PM on one PM grid, scored against the 4096^2
grid-matched control.  Arms run ONE AT A TIME, each in its own process, so
wall / CPU / peak-RSS numbers are clean.

Usage: python3 exp03_gate.py NG [--repeat]      (NG = 256 or 512)
       -> results/exp03_gate_ng<NG>.json        (resumable)
`--repeat` adds two repetitions of three key arms for timing variance."""
from __future__ import annotations

import sys

import bench


def ctrl(**kw):
    d = dict(indicator="pred", tol=1e-6, tau=0.6, tau_p=0.3, interval=10, cooldown=10)
    d.update(kw)
    return d


def arms(ng, repeat=False):
    A = {}
    lat = lambda sub0: dict(mode="lattice", sub0=sub0)
    # -- N-body ladder -------------------------------------------------------------------
    for f in ((1, 1.5, 2, 3, 4, 6) if ng <= 256 else (1, 1.5, 2, 3, 4)):
        n = int(round(f * ng))
        A[f"nbody/n{n}"] = bench.nbody(ng, n)
    # -- fixed sheets: 2ng and 4ng deposit lattices --------------------------------------
    for (n0, s) in ((ng // 4, 8), (ng // 2, 4), (ng, 2), (ng // 2, 8), (ng, 4)):
        A[f"sheet-quin/n{n0}s{s}"] = bench.sheet(ng, n0, k=5, quad=lat(s))
    for (n0, s) in ((ng, 2), (ng, 4)):
        A[f"sheet-cub/n{n0}s{s}"] = bench.sheet(ng, n0, k=3, quad=lat(s))
    if ng <= 256:
        A[f"sheet-quin/n{ng}s8"] = bench.sheet(ng, ng, k=5, quad=lat(8))
        A[f"sheet-quin/n{2 * ng}s4"] = bench.sheet(ng, 2 * ng, k=5, quad=lat(4))
    # -- uniform early refinement ng/2 -> ng (the baseline refine2 control), 4ng lattice ---
    for st in (100, 150, 200):
        A[f"href-quin/{ng // 2}to{ng}@t{st / 100:.1f}"] = bench.sheet(ng, ng // 2, lmax=1, k=5, quad=lat(8), refine_all_at=[st])
    # -- local h-refinement, 4ng lattice -----------------------------------------------------
    A[f"local-quin/n{ng // 2}L1 pred"] = bench.sheet(ng, ng // 2, lmax=1, k=5, quad=lat(8), ctrl=ctrl())
    A[f"local-quin/n{ng // 4}L2 pred"] = bench.sheet(ng, ng // 4, lmax=2, k=5, quad=lat(16), ctrl=ctrl())
    A[f"local-quin/n{ng // 2}L1 d6"] = bench.sheet(ng, ng // 2, lmax=1, k=5, quad=lat(8), ctrl=ctrl(indicator="d6"))
    A[f"local-quin/n{ng // 2}L1 pred tol1e-5"] = bench.sheet(ng, ng // 2, lmax=1, k=5, quad=lat(8), ctrl=ctrl(tol=1e-5))
    A[f"local-quin/n{ng // 2}L1 pred gate=crossing"] = bench.sheet(ng, ng // 2, lmax=1, k=5, quad=lat(8), ctrl=ctrl(tol=1e-7, enable_on="crossing"))
    A[f"local-cub/n{ng // 2}L1 pred"] = bench.sheet(ng, ng // 2, lmax=1, k=3, quad=lat(8), ctrl=ctrl())
    # local refinement BEYOND the PM-cell spacing, under a hard node budget of 1.5 ng^2
    A[f"local-quin/n{ng}L1 pred budget1.5"] = bench.sheet(ng, ng, lmax=1, k=5, quad=lat(4), ctrl=ctrl(max_nodes=int(1.5 * ng * ng)))
    # -- local h-refinement, 2ng lattice (looser tolerance matched to the coarser quadrature) --
    A[f"local-quin/n{ng // 2}L1 pred tol1e-5 (2ng lattice)"] = bench.sheet(ng, ng // 2, lmax=1, k=5, quad=lat(4), ctrl=ctrl(tol=1e-5))
    A[f"local-quin/n{ng // 4}L2 pred tol3e-5 (2ng lattice)"] = bench.sheet(ng, ng // 4, lmax=2, k=5, quad=lat(8), ctrl=ctrl(tol=3e-5))
    if repeat:
        for r in (2, 3):
            A[f"nbody/n{3 * ng}#r{r}"] = bench.nbody(ng, 3 * ng)
            A[f"sheet-quin/n{ng}s4#r{r}"] = bench.sheet(ng, ng, k=5, quad=lat(4))
            A[f"local-quin/n{ng // 2}L1 pred#r{r}"] = bench.sheet(ng, ng // 2, lmax=1, k=5, quad=lat(8), ctrl=ctrl())
    return A


if __name__ == "__main__":
    ng = int(sys.argv[1])
    bench.run_arms(f"exp03_gate_ng{ng}.json", arms(ng, "--repeat" in sys.argv), jobs=1, meta_extra=dict(ng=ng), script=__file__)
