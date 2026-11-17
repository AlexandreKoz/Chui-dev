"""Trigger-timing study (PM 256, 4ng deposit lattice, quintic, n0 = 128 -> 256).

Which state-only trigger refines early enough?  Compared: the truncation
monitor "d6" and the predictive (look-ahead) monitor at several tolerances
and horizons; the gates "precaustic" (strain time-to-caustic), "crossing"
(first det Dx <= 0: the late trigger) and a clock gate at t = 1.5; and
prescribed uniform refinement at t = 1.0, 1.5, 2.0.  Accuracy is scored
offline against the 4096^2 control; timings here are indicative only (several
arms run concurrently) -- clean timings are in exp03.

Usage: python3 exp04_triggers.py [NG]     -> results/exp04_triggers_ng<NG>.json
"""
from __future__ import annotations

import sys

import bench

if __name__ == "__main__":
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    q = dict(mode="lattice", sub0=8)

    def C(**kw):
        d = dict(indicator="pred", tol=1e-6, interval=10, cooldown=10, tau=0.3, tau_p=0.3)
        d.update(kw)
        return d
    A = {}
    A[f"sheet-quin/n{ng}s4 (uniform from t=0)"] = bench.sheet(ng, ng, quad=dict(mode="lattice", sub0=4))
    A[f"sheet-quin/n{ng // 2}s8 (never refined)"] = bench.sheet(ng, ng // 2, quad=q)
    for st in (100, 150, 175, 200):
        A[f"href/uniform@t{st / 100:.2f}"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, refine_all_at=[st])
    for tau in (0.3, 0.6):
        for tol in (1e-5, 1e-6, 1e-7):
            A[f"local/pred tau{tau} tol{tol:g}"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tau=tau, tol=tol))
    for tol in (1e-5, 1e-6, 1e-7):
        A[f"local/d6 tol{tol:g}"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(indicator="d6", tol=tol))
    A["local/pred tol1e-7 gate=precaustic(0.3)"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tol=1e-7, enable_on="precaustic", tau_lead=0.3))
    A["local/pred tol1e-7 gate=precaustic(0.6)"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tol=1e-7, enable_on="precaustic", tau_lead=0.6))
    A["local/pred tol1e-7 gate=crossing"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tol=1e-7, enable_on="crossing"))
    A["local/pred tol1e-7 gate=t1.5"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tol=1e-7, t_enable=1.5))
    A["local/pred tol1e-7 gate=t1.0"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tol=1e-7, t_enable=1.0))
    A["local/pred tau0.6 tol1e-6 buffer1"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tau=0.6, buffer=1))
    A["local/pred tau0.6 tol1e-6 halo2"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tau=0.6), halo_tiles=2)
    A["local/pred tau0.6 tol1e-6 sweeps2"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tau=0.6), sweeps=2)
    A["local/pred tau0.6 tol1e-6 budget0.6"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tau=0.6, max_nodes=int(0.6 * ng * ng)))
    A["local/pred tau0.6 tol1e-6 no-rebase"] = bench.sheet(ng, ng // 2, lmax=1, quad=q, ctrl=C(tau=0.6, auto_rebase=False))
    bench.run_arms(f"exp04_triggers_ng{ng}.json", A, jobs=3, meta_extra=dict(ng=ng), script=__file__)
