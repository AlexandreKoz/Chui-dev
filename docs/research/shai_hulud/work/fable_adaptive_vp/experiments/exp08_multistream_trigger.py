"""Alternative trigger investigated after the truncation monitors failed at
PM 512 (exp03, exp07): Eulerian multi-stream proximity ("ms").

Finding that motivated it: after shell crossing a sheet needs nodes at about
half a PM cell (n = 2 ng) to match N-body on the same deposit lattice, and
the blocks that limit the accuracy are those that are about to enter the
multi-stream region -- where the sheet is still smooth and every truncation
monitor is small.  The "ms" trigger refines every block that is inside, or
within a ballistic lead time of, the PM cells that contain folded mass.

Usage: python3 exp08_multistream_trigger.py [NG]   -> results/exp08_multistream_trigger_ng<NG>.json
(two arms at a time: timings indicative)."""
from __future__ import annotations

import sys

import bench

if __name__ == "__main__":
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 512
    lat = lambda s: dict(mode="lattice", sub0=s)

    def M(**kw):
        d = dict(indicator="ms", tol=0.5, tau=0.6, tau_p=0.3, interval=10, cooldown=10, floor_cells=0.0, ms_dilate=2)
        d.update(kw)
        return d
    A = {}
    # deposit lattice 4 ng, target n = 2 ng inside the multi-stream region
    A[f"local-ms/n{ng}L1 lead0.6 dilate2"] = bench.sheet(ng, ng, lmax=1, quad=lat(4), ctrl=M())
    A[f"local-ms/n{ng}L1 lead0.3 dilate2"] = bench.sheet(ng, ng, lmax=1, quad=lat(4), ctrl=M(tau=0.3))
    A[f"local-ms/n{ng}L1 lead1.0 dilate4"] = bench.sheet(ng, ng, lmax=1, quad=lat(4), ctrl=M(tau=1.0, ms_dilate=4))
    A[f"local-ms/n{ng // 2}L2 lead0.6 dilate2"] = bench.sheet(ng, ng // 2, lmax=2, quad=lat(8), ctrl=M())
    A[f"local-ms/n{ng // 4}L3 lead0.6 dilate2"] = bench.sheet(ng, ng // 4, lmax=3, quad=lat(16), ctrl=M())
    # references of the same deposit lattice (re-run here so that the table is self-contained)
    A[f"sheet-quin/n{2 * ng}s2 (uniform, n = 2 ng)"] = bench.sheet(ng, 2 * ng, quad=lat(2))
    A[f"sheet-quin/n{ng}s4 (uniform, n = ng)"] = bench.sheet(ng, ng, quad=lat(4))
    A[f"nbody/n{4 * ng}"] = bench.nbody(ng, 4 * ng)
    bench.run_arms(f"exp08_multistream_trigger_ng{ng}.json", A, jobs=2, meta_extra=dict(ng=ng), script=__file__)
