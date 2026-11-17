"""Grid-matched particle controls ("references") for the 2D self-gravitating
sheet: N-body PM with n_ref^2 equal-mass particles on the SAME PM grid,
kernel, force law, time step and initial data as every test arm.  The
continuum-in-label limit of the discrete model on that PM grid is what both
the sheet and the N-body arms approximate, so comparing against it isolates
the sheet / sampling error from the PM field error.

Usage:  python3 exp01_reference.py NG N_REF          (e.g. 256 2048)

Writes  <DATA>/vp_ref_ng<NG>_n<N_REF>.npz  (x, p on the 512^2 label lattice,
PM density, at steps 100..300) and results/reference_ng<NG>_n<N_REF>.json
(conserved quantities, folded fractions, timing, sha256 of the npz).
OFFLINE EVALUATION DATA ONLY.  Resumable: a checkpoint is written every 50
steps and picked up on restart.
"""
from __future__ import annotations

import hashlib
import os
import sys
import time

import numpy as np

import common
from shadapt.solver import NBodyVP
from shsheet.vp import GrowingModeIC

DT, NSTEPS = 0.01, 300
SNAP = (100, 150, 200, 250, 300)
LATTICE = 512

if __name__ == "__main__":
    ng, n = int(sys.argv[1]), int(sys.argv[2])
    os.makedirs(common.DATA, exist_ok=True)
    tag = f"ng{ng}_n{n}"
    final = os.path.join(common.DATA, f"vp_ref_{tag}.npz")
    ck = os.path.join(common.DATA, f"vp_ref_{tag}.ckpt.npz")
    if os.path.exists(final):
        print("exists:", final)
        sys.exit(0)
    run = NBodyVP(GrowingModeIC(), n, ng)
    out = dict(n=n, ng=ng, dt=DT, snaps=np.array(SNAP), lattice=LATTICE)
    cons = {0: run.conserved()}
    k0 = 0
    wall_prev = 0.0
    if os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        run.x, run.p, run.S = z["x"], z["p"], z["S"]
        run.t, run.step_count, k0 = float(z["t"]), int(z["step"]), int(z["step"])
        out.update({k: z[k] for k in z.files if k[0] in "xpr" and k[1:].isdigit() or k.startswith(("rho", "negdet", "detmin"))})
        cons = z["cons"].item()
        wall_prev = float(z["wall"])
        run.solve_field()
        print(f"resumed from step {k0}", flush=True)
    st = n // LATTICE
    t0 = time.time()
    for k in range(k0 + 1, NSTEPS + 1):
        run.step(DT)
        if k in SNAP:
            out[f"x{k}"] = run.x.reshape(n, n, 2)[::st, ::st].copy()
            out[f"p{k}"] = run.p.reshape(n, n, 2)[::st, ::st].copy()
            out[f"rho{k}"] = run.rho.copy()
            cons[k] = run.conserved()
            fr, dm = run.folded_fraction()
            out[f"negdet{k}"], out[f"detmin{k}"] = fr, dm
            print(f"ng={ng} n={n} t={run.t:.2f} E={cons[k]['E']:+.6e} P=({cons[k]['Px']:+.1e},{cons[k]['Py']:+.1e}) "
                  f"folded label fraction {fr:.4f} rho_max {run.rho.max():.1f} elapsed {wall_prev + time.time() - t0:.0f}s", flush=True)
        if k % 50 == 0 and k < NSTEPS:
            np.savez(ck, x=run.x, p=run.p, S=run.S, t=run.t, step=k, cons=np.array(cons, dtype=object), wall=wall_prev + time.time() - t0,
                     **{kk: v for kk, v in out.items() if kk not in ("n", "ng", "dt", "snaps", "lattice")})
    wall = wall_prev + time.time() - t0
    out["wall_s"] = wall
    np.savez(final, **out)
    if os.path.exists(ck):
        os.remove(ck)
    h = hashlib.sha256(open(final, "rb").read()).hexdigest()
    common.save(f"reference_{tag}.json", dict(meta=common.metadata(__file__, ng=ng, n=n, dt=DT, nsteps=NSTEPS), file=os.path.basename(final), sha256=h,
                                              wall_s=wall, peak_rss_mb=common.peak_rss_mb(), cpu=run.cpu,
                                              conserved={str(k): v for k, v in cons.items()},
                                              folded={str(k): float(out[f"negdet{k}"]) for k in SNAP},
                                              rho_max={str(k): float(out[f"rho{k}"].max()) for k in SNAP}))
    print("[saved]", final, h[:16])
