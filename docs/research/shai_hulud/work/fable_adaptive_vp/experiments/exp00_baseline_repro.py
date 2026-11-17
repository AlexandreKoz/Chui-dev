"""Step 1 of Call 1 ("baseline first"): reproduce saved Vlasov-Poisson arms of
the Fable breakthrough package and their fixed-PM-grid control with the
ORIGINAL, unmodified baseline code, in the current environment.

Usage:
  python3 exp00_baseline_repro.py /path/to/pristine/work/fable_breakthrough

Re-runs  (i) the 512^2 particle control at PM 64^2 (exp_c_reference.py 512)
and compares it with the saved results/vp_reference_n512.npz, and
(ii) the arms nbody/n128, sheet-cub/n32s4, sheet-quin/n64s2, href-cub/32to64@150
of exp_c_vp.py against the saved 1024^2 control, and compares every number
with the saved results/exp_c_vp.json.  Writes results/exp00_baseline_repro.json.
Nothing is written into the pristine tree (the control is re-run in memory)."""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np

import common  # noqa: F401  (metadata helpers; does NOT put the baseline on the path)

ARMS = ("nbody/n128", "sheet-cub/n32s4", "sheet-quin/n64s2", "href-cub/32to64@150")


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


if __name__ == "__main__":
    base = os.path.abspath(sys.argv[1])
    # the pristine baseline package must win over the vendored copy for this script
    sys.path.insert(0, os.path.join(base, "src"))
    sys.path.insert(0, os.path.join(base, "experiments"))
    for m in [k for k in sys.modules if k == "common" or k.startswith("shsheet")]:
        del sys.modules[m]
    expc = _load(os.path.join(base, "experiments", "exp_c_vp.py"), "exp_c_vp_baseline")
    from shsheet.pm import PM
    from shsheet.sheet import Sheet
    from shsheet.vp import GrowingModeIC, VPRun

    import common as base_common  # the baseline's own helper (for its source hash)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    out = dict(meta=dict(command="python3 " + " ".join(sys.argv), python=sys.version.split()[0], numpy=np.__version__,
                         baseline_source_sha256_16=base_common.source_hash(), timestamp=time.strftime("%Y-%m-%dT%H:%M:%S")),
               frozen=dict(time_variable="t (static background, no expansion)", integrator="KDK, field solved after the drift",
                           dt=expc.DT, nsteps=expc.NSTEPS, pm_grid=expc.NG, kernel="PCS (order 4)", force_law=expc.FORCE,
                           ic="GrowingModeIC(eps=(0.10,0.07,0.03), phase=0.9)", units="L=1, total mass 1, 4 pi G = 1"),
               control={}, arms={})
    saved = json.load(open(os.path.join(base, "results", "exp_c_vp.json")))
    out["meta"]["saved_meta"] = saved["meta"]

    # (i) control: 512^2 particles at PM 64^2, compared with the saved snapshot file
    ref512 = np.load(os.path.join(base, "results", "vp_reference_n512.npz"))
    run = VPRun(Sheet.from_ic(GrowingModeIC(), 512), PM(expc.NG, expc.ORDER, expc.FORCE), "nodes")
    t0 = time.perf_counter()
    for k in range(1, expc.NSTEPS + 1):
        run.step(expc.DT)
        if k in expc.SNAP:
            x = run.sh.x[::4, ::4]
            p = run.sh.p[::4, ::4]
            out["control"][str(k)] = dict(x_maxdiff=float(np.abs(x - ref512[f"x{k}"]).max()), p_maxdiff=float(np.abs(p - ref512[f"p{k}"]).max()),
                                          rho_maxdiff=float(np.abs(run.rho - ref512[f"rho{k}"]).max()))
            print(f"control n=512 step {k}: {out['control'][str(k)]}", flush=True)
    out["control"]["wall_s"] = time.perf_counter() - t0
    out["control"]["saved_wall_s"] = float(ref512["wall_s"])

    # (ii) arms against the saved 1024^2 control
    ref = np.load(os.path.join(base, "results", "vp_reference_n1024.npz"))
    todo = expc.arms(False)
    worst = 0.0
    for name in ARMS:
        r = expc.run_arm(todo[name], ref)
        cmp = {}
        for s, v in r["snaps"].items():
            sv = saved["runs"][name]["snaps"][s]
            cmp[s] = {q: dict(now=v[q], saved=sv[q], rel=abs(v[q] - sv[q]) / max(abs(sv[q]), 1e-300)) for q in ("x_rms", "p_rms", "rho_rel_l2", "E")}
            worst = max(worst, max(c["rel"] for c in cmp[s].values()))
        out["arms"][name] = dict(compare=cmp, wall_s=r["wall_s"], saved_wall_s=saved["runs"][name]["wall_s"])
        print(f"{name}: max rel deviation from saved numbers {max(c['rel'] for s in cmp.values() for c in s.values()):.2e}  wall {r['wall_s']:.1f}s (saved {saved['runs'][name]['wall_s']:.1f}s)", flush=True)
    out["max_relative_deviation_arms"] = worst
    out["control_bit_identical"] = bool(all(v["x_maxdiff"] == 0 and v["p_maxdiff"] == 0 and v["rho_maxdiff"] == 0 for k, v in out["control"].items() if k.isdigit()))
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"), exist_ok=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "exp00_baseline_repro.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print("control bit-identical:", out["control_bit_identical"], " worst arm deviation:", worst)
