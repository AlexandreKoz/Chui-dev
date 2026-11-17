"""Resumable multi-arm benchmark runner.  Every arm runs in its own process
(`arm.py`), results are saved after each arm, and arms already present in
the output file are skipped (delete the entry to re-run it)."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import common

DT, NSTEPS = 0.01, 300
SNAP = [100, 150, 200, 250, 300]


def ref_path(ng, n=None):
    """Finest available grid-matched control for this PM grid."""
    for nn in ([n] if n else [4096, 2048]):
        p = os.path.join(common.DATA, f"vp_ref_ng{ng}_n{nn}.npz")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"no control for PM {ng}: run exp01_reference.py {ng} 2048")


def _one(name, spec):
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("NUMBA_NUM_THREADS", "1")
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(common.HERE, "arm.py"), json.dumps(spec)], capture_output=True, text=True, env=env)
    for line in p.stdout.splitlines():
        if line.startswith("RESULT "):
            r = json.loads(line[7:])
            r["total_elapsed_s"] = time.time() - t0
            return name, r
    return name, dict(failed=True, spec=spec, returncode=p.returncode, stderr=p.stderr[-3000:], stdout=p.stdout[-1000:])


def fmt(name, r, snaps=("150", "200", "300")):
    if r.get("failed"):
        return f"{name:44s} FAILED: {r['stderr'][-300:]}"
    s = "  ".join(f"t{int(k) / 100:.1f}: x {r['snaps'][k]['x_rms']:.1e} rho {r['snaps'][k]['rho_rel_l2']:.1e} F {r['snaps'][k]['force_rel_l2']:.1e}"
                  for k in snaps if k in r["snaps"])
    return (f"{name:44s} {s} | wall {r['wall_s']:6.1f}s rss+ {r['peak_rss_incr_mb']:6.1f}MB nodes {r['nodes_final']:7d} "
            f"samp/step {r['deposit_samples_per_step']:9.0f}")


def run_arms(out_name, arms, jobs=1, meta_extra=None, script=__file__):
    res = common.load_partial(out_name, dict(meta=common.metadata(script, dt=DT, nsteps=NSTEPS, jobs=jobs, **(meta_extra or {})), runs={}))
    todo = [(n, s) for n, s in arms.items() if n not in res["runs"] or res["runs"][n].get("failed")]
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        for name, r in ex.map(lambda a: _one(*a), todo):
            r["concurrent_jobs"] = jobs
            res["runs"][name] = r
            common.save(out_name, res)
            print(fmt(name, r), flush=True)
    return res


def nbody(ng, n, **kw):
    return dict(kind="nbody", ng=ng, n=n, dt=DT, nsteps=NSTEPS, snaps=SNAP, ref=ref_path(ng), **kw)


def sheet(ng, n0, lmax=0, k=5, quad=None, ctrl=None, **kw):
    return dict(kind="sheet", ng=ng, n0=n0, lmax=lmax, k=k, B=kw.pop("B", 8), quad=quad or dict(mode="lattice", sub0=4), ctrl=ctrl,
                dt=DT, nsteps=NSTEPS, snaps=SNAP, ref=ref_path(ng), **kw)
