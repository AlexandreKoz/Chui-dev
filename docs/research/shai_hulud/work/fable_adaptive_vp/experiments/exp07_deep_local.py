"""Local refinement BEYOND the PM-cell node spacing (PM 512 by default).

At PM 512^2 a uniform sheet with n = ng falls behind N-body after t = 2
although it shares its deposit lattice: the position error at t = 3 is the
same for cubic and quintic carriers, i.e. it is under-resolution of a
localised structure, not truncation of a smooth one.  This experiment asks
whether LOCAL levels on top of n = ng repair it, and what the resolution
floor (no split of blocks whose cells are already smaller than `floor_cells`
PM cells) does.  Hard node budgets in units of ng^2.

Usage: python3 exp07_deep_local.py [NG] [--clean]   -> results/exp07_deep_local_ng<NG>.json
(--clean: one arm at a time).  The budget arm WITH the default floor is in exp03."""
from __future__ import annotations

import sys

import bench

if __name__ == "__main__":
    ng = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 512
    clean = "--clean" in sys.argv
    lat = lambda s: dict(mode="lattice", sub0=s)

    def C(**kw):
        d = dict(indicator="pred", tol=1e-6, tau=0.6, tau_p=0.3, interval=10, cooldown=10)
        d.update(kw)
        return d
    A = {}
    A[f"sheet-quin/n{ng}s4 (uniform, n = ng)"] = bench.sheet(ng, ng, quad=lat(4))
    A[f"sheet-quin/n{2 * ng}s2 (uniform, n = 2 ng)"] = bench.sheet(ng, 2 * ng, quad=lat(2))
    # tolerance-driven (no budget, no floor): refine wherever the monitor exceeds a LOOSE tolerance
    for tol in (1e-3, 1e-4):
        A[f"local/n{ng}L1 tol{tol:g} floor0 nobudget"] = bench.sheet(ng, ng, lmax=1, quad=lat(4), ctrl=C(tol=tol, floor_cells=0.0))
    # budget-driven with the tight default tolerance: the budget is spent in the order of the monitor
    A[f"local/n{ng}L1 tol1e-6 floor0 budget1.5"] = bench.sheet(ng, ng, lmax=1, quad=lat(4), ctrl=C(floor_cells=0.0, max_nodes=int(1.5 * ng * ng)))
    A[f"local/n{ng}L2 tol1e-6 floor0 budget2.0"] = bench.sheet(ng, ng, lmax=2, quad=lat(4), ctrl=C(floor_cells=0.0, max_nodes=int(2.0 * ng * ng)))
    bench.run_arms(f"exp07_deep_local_ng{ng}.json", A, jobs=1 if clean else 2, meta_extra=dict(ng=ng, clean=clean), script=__file__)
