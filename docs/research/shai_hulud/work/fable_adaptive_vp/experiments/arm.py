"""One benchmark arm in its own process (clean wall / CPU / peak-RSS numbers).

Usage:  python3 arm.py '<json spec>'      -> prints one JSON line "RESULT {...}"

spec (all arms): kind, ng, dt, nsteps, snaps, ref (npz path), scratch (dir)
  kind = "nbody": n
  kind = "sheet": n0, lmax, k, B, halo_tiles, sweeps,
                  quad = {mode, sub0, target, nq_min, nq_max},
                  ctrl = {Controller kwargs} or null,
                  refine_all_at = [steps]   (uniform early refinement control)

Protocol.  Phase 1 (timed): the run.  At snapshot steps the state needed for
evaluation is written to the scratch directory (the time spent writing is
excluded from the wall clock; peak RSS is read at the end of phase 1).
Phase 2 (untimed): the control is loaded and every metric is computed from
the saved states.  The control is therefore never in memory while the arm
runs, and never visible to the solver or the controller.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

import common
from shadapt import diagnostics as dg
from shadapt.controller import Controller
from shadapt.sheet import AdaptiveSheet
from shadapt.solver import AdaptiveVP, NBodyVP, Quadrature
from shsheet.vp import GrowingModeIC


def warmup():
    """Compile every numba kernel and touch the code paths before timing."""
    ic = GrowingModeIC()
    nb = NBodyVP(ic, 16, 16)
    nb.step(0.01)
    sh = AdaptiveSheet.from_ic(ic, 16, B=8, lmax=1, k=5)
    run = AdaptiveVP(sh, 16, Quadrature("lattice", 2))
    run.step(0.01)
    sh.refine([(0, 0, 0)])
    run.step(0.01)
    run.conserved()


def run_nbody(spec, scratch):
    run = NBodyVP(GrowingModeIC(), spec["n"], spec["ng"])
    c0 = run.conserved()
    wall = 0.0
    cpu0 = time.process_time()
    io_cpu = 0.0
    saved = {}
    for k in range(1, spec["nsteps"] + 1):
        t0 = time.perf_counter()
        run.step(spec["dt"])
        wall += time.perf_counter() - t0
        if k in spec["snaps"]:
            c1 = time.process_time()
            f = os.path.join(scratch, f"nb_{k}.npz")
            run.sync()
            np.savez(f, x=run.x, p=run.p, rho=run.rho)
            saved[k] = f
            io_cpu += time.process_time() - c1
    res = dict(wall_s=wall, cpu_s=time.process_time() - cpu0 - io_cpu, peak_rss_mb=common.peak_rss_mb(), cpu=dict(run.cpu),
               memory=run.memory_report(), deposit_samples_per_step=spec["n"] ** 2, nodes_final=spec["n"] ** 2,
               node_steps=spec["n"] ** 2 * spec["nsteps"], E0=c0["E"], K0=c0["K"])
    # phase 2
    ref = np.load(spec["ref"])
    snaps = {}
    for k, f in saved.items():
        z = np.load(f)
        run.x, run.p, run.rho = z["x"], z["p"], z["rho"]
        run.phi = run.pm.potential(run.rho)
        r = dg.compare_nbody(run, ref, k)
        c = run.conserved()
        fr, dm = run.folded_fraction()
        r.update(E=c["E"], K=c["K"], U=c["U"], Px=c["Px"], Py=c["Py"], mass=c["mass"], folded_fraction=fr, det_min=dm,
                 ref_folded_fraction=float(ref[f"negdet{k}"]))
        snaps[str(k)] = r
        os.remove(f)
    res["snaps"] = snaps
    return res


def run_sheet(spec, scratch):
    ic = GrowingModeIC()
    sh = AdaptiveSheet.from_ic(ic, spec["n0"], B=spec.get("B", 8), lmax=spec.get("lmax", 0), k=spec.get("k", 5),
                               halo_tiles=spec.get("halo_tiles", 1), sweeps=spec.get("sweeps", 1))
    q = spec.get("quad", dict(mode="lattice", sub0=4))
    quad = Quadrature(q.get("mode", "lattice"), q.get("sub0", 4), q.get("target", 0.5), q.get("nq_min"), q.get("nq_max"),
                      rule=q.get("rule", "midpoint"))
    ctrl = Controller(**spec["ctrl"]) if spec.get("ctrl") else None
    run = AdaptiveVP(sh, spec["ng"], quad, ctrl)
    c0 = run.conserved()
    wall = 0.0
    cpu0 = time.process_time()
    io_cpu = 0.0
    saved = {}
    refine_all = set(spec.get("refine_all_at", []))
    samples_hist, nodes_hist = [], []
    manual = []
    for k in range(1, spec["nsteps"] + 1):
        t0 = time.perf_counter()
        run.step(spec["dt"])
        if k in refine_all:
            # uniform early refinement control (the baseline `refine2`): split everything, then relabel
            run.sync()
            sh.refine(sh.f.leaves(), reason="uniform", step=k, t=run.t)
            if sh.rebase():
                run.quad.sub0 = max(1, run.quad.sub0 // 2)
            run.quad.sync(sh.car)
            run.solve_field()
            manual.append(dict(step=k, t=run.t, nodes=sh.n_nodes()))
        wall += time.perf_counter() - t0
        if k % 10 == 0 or k == spec["nsteps"]:
            samples_hist.append(run.last_samples)
            nodes_hist.append(sh.n_nodes())
        if k in spec["snaps"]:
            c1 = time.process_time()
            f = os.path.join(scratch, f"sh_{k}.npz")
            run.sync()
            run.save_checkpoint(f)
            saved[k] = f
            io_cpu += time.process_time() - c1
    res = dict(wall_s=wall, cpu_s=time.process_time() - cpu0 - io_cpu, peak_rss_mb=common.peak_rss_mb(), cpu=dict(run.cpu),
               memory=run.memory_report(), deposit_samples_total=run.deposit_samples,
               deposit_samples_per_step=run.deposit_samples / max(run.n_field_solves, 1), nodes_final=sh.n_nodes(),
               nodes_per_level_final=sh.nodes_per_level(), node_steps=run.node_steps, E0=c0["E"], K0=c0["K"],
               samples_every_10=samples_hist, nodes_every_10=nodes_hist, n_field_solves=run.n_field_solves,
               carrier_layouts=sh.car.stats["n_layout"], carrier_layout_s=sh.car.stats["layout_s"], manual_refinements=manual,
               leaf_blocks_final=len(sh.f.leaves()), leaf_cells_final=sh.f.n_leaf_cells(),
               refined_area_fraction=[len(sh.f.exists[l]) / float(sh.f.nbl(l) ** 2) for l in range(1, sh.f.top_level() + 1)])
    if ctrl is not None:
        ev = sh.f.events
        res["controller"] = dict(params=ctrl.params(), first_split_t=ctrl.first_split_t, latch_step=ctrl.latch_step,
                                 n_events=len(ev), n_indicator=sum(1 for e in ev if e["reason"].startswith("indicator")),
                                 n_balance=sum(1 for e in ev if e["reason"] == "balance"),
                                 n_buffer=sum(1 for e in ev if e["reason"] == "buffer"),
                                 n_merge=sum(1 for e in ev if e["kind"] == "merge"),
                                 denied_budget=sum(h["n_denied_budget"] for h in ctrl.history),
                                 denied_floor_evals=sum(1 for h in ctrl.history if h["n_denied_floor"]),
                                 denied_level_evals=sum(1 for h in ctrl.history if h["n_denied_level"]),
                                 all_splits_have_reason=all("reason" in e for e in ev),
                                 history=[{kk: h[kk] for kk in ("step", "t", "eta_max", "eta_d6_max", "eta_pred_max", "det_min", "tau_c_min",
                                                                 "enabled", "nodes", "n_split", "n_balance", "n_merge", "n_denied_budget")}
                                          for h in ctrl.history])
        res["events_head"] = ev[:40]
    # phase 2 (untimed, control loaded only now)
    ref = np.load(spec["ref"])
    snaps = {}
    for k, f in saved.items():
        r2 = AdaptiveVP.load_checkpoint(f)       # carries the density/potential the run actually used
        r = dg.compare_sheet(r2, ref, k, lattice=spec.get("lattice", 256))
        c = r2.conserved()
        r.update(E=c["E"], K=c["K"], U=c["U"], Px=c["Px"], Py=c["Py"], mass=c["grid_mass"], ledger_mass=c["ledger_mass"],
                 folded_fraction=c["folded_fraction"], det_min=c["det_min"], ref_folded_fraction=float(ref[f"negdet{k}"]),
                 nodes=r2.sh.n_nodes(), nodes_per_level=r2.sh.nodes_per_level(), samples=r2.last_samples)
        if spec.get("action_diag", True):
            la = dg.loop_actions(r2.sh, nlines=8, oversample=2)
            r.update(loop_action_max_dev=la["max_deviation"], loop_action_spread=la["spread"])
            if k == max(saved):
                r.update(nodal_action_defect=dg.nodal_action_defect(r2.sh, nsamp=2000))
                r.update(interface_jumps=dg.continuity_jumps(r2.sh, "xp") if r2.sh.top >= 1 else None)
        snaps[str(k)] = r
        os.remove(f)
    res["snaps"] = snaps
    return res


if __name__ == "__main__":
    spec = json.loads(sys.argv[1])
    scratch = spec.get("scratch") or os.path.join(common.DATA, "scratch")
    scratch = os.path.join(scratch, f"arm_{os.getpid()}")
    os.makedirs(scratch, exist_ok=True)
    warmup()
    base = common.current_rss_mb()
    res = run_nbody(spec, scratch) if spec["kind"] == "nbody" else run_sheet(spec, scratch)
    res["rss_baseline_mb"] = base
    res["peak_rss_incr_mb"] = res["peak_rss_mb"] - base
    res["spec"] = spec
    try:
        os.rmdir(scratch)
    except OSError:
        pass
    print("RESULT " + json.dumps(res, default=float))
