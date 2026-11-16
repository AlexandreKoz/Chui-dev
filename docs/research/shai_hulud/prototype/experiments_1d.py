"""Numerical experiments for the d = 1 SHAI-HULUD reference prototype.

Usage:  python3 experiments_1d.py [quick]

Writes ../results/exp1d_*.json.  All errors are measured against the
*oracle*: labels integrated from the initial condition with the same
integrator and step sequence, i.e. exact samples of the transported carrier
Psi^n(Lambda_0).  Time-integration error therefore cancels identically and
what is reported is purely representation + remapping error.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from sheet1d import (TWO_PI, CellData, ComovingEdS, DrivenPendulum, Pendulum, Sheet, SineIC,
                     affine_curvature, carrier_period, chord_defects, evaluate, integrate,
                     monitor_density, omega, oracle, oracle_tangent, ot_labels, polygon_period,
                     total_mass)

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
SUB = np.array([0.125, 0.375, 0.625, 0.875])  # sub-cell sampling fractions


def diagnostics(model, ic, sheet, method, t0, dt, nsteps, t, n_dense=8192):
    """All diagnostics of a sheet state against the oracle (one oracle pass)."""
    n = sheet.n
    out = {"t": t, "nsteps": nsteps, "n": n}
    dq = np.diff(np.append(sheet.q, sheet.q[0] + TWO_PI))
    qs = (sheet.q[:, None] + SUB[None, :] * dq[:, None]).ravel()
    qd = (np.arange(n_dense) + 0.5) / n_dense * TWO_PI
    qa = np.concatenate([sheet.q, qs, qd])
    xa, pa, Sa, txa, tpa = oracle_tangent(model, ic, qa, t0, dt, nsteps)
    i1, i2 = n, n + qs.size
    # (i) nodal errors against the oracle at the node labels
    xo, po, So = xa[:i1], pa[:i1], Sa[:i1]
    out["nodal_err"] = float(np.max(np.hypot(sheet.x - xo, sheet.p - po)))
    dS = (sheet.S - So).astype(float)
    out["action_err"] = float(np.max(np.abs(dS - dS[0])))
    # (ii) representation error, sub-sampled in every cell of the method's mesh
    xs, ps, tx, tp = xa[i1:i2], pa[i1:i2], txa[i1:i2], tpa[i1:i2]
    cd = CellData(sheet, method)
    xr, pr, _ = evaluate(sheet, qs, method, cd=cd)
    ex, ep = xr - xs, pr - ps
    spd = np.hypot(tx, tp)
    en = omega(tx, tp, ex, ep) / spd
    et = (tx * ex + tp * ep) / spd
    w = np.repeat(dq, SUB.size) / SUB.size / TWO_PI  # mass weights
    out["en_max"] = float(np.max(np.abs(en)))
    out["et_max"] = float(np.max(np.abs(et)))
    out["en_l2m"] = float(np.sqrt(np.sum(w * en**2)))  # mass-weighted L2
    out["et_l2m"] = float(np.sqrt(np.sum(w * et**2)))
    # (iii) Poincare invariant of the method's own carrier
    out["period_err"] = carrier_period(sheet, method) - sheet.c
    out["polygon_period_err"] = polygon_period(sheet) - sheet.c
    # (iv) mass
    m, dqmin = total_mass(sheet)
    out["mass_err"] = m - 1.0
    out["dq_min"] = dqmin
    out["dq_max"] = float(dq.max())
    # (v) chord (P1) action defects: equidistribution outcome
    kap, _ = affine_curvature(sheet)
    out["chord_defect_max_geom"] = float(np.max(np.abs(kap) * dq**3 / 12.0))
    if METH_HAS_ACTION.get(method, False):
        out["chord_defect_max_action"] = float(np.abs(chord_defects(sheet)).max())
    # (vi) energy of the carrier vs oracle on dense uniform labels
    xd, pd = xa[i2:], pa[i2:]
    xe, pe, _ = evaluate(sheet, qd, method, cd=cd)
    out["energy_err"] = float(np.mean(model.energy(xe, pe, t)) - np.mean(model.energy(xd, pd, t)))
    out["energy_ref"] = float(np.mean(model.energy(xd, pd, t)))
    # mass-weighted errors on dense uniform labels (measure-level accuracy)
    out["dense_err_max"] = float(np.max(np.hypot(xe - xd, pe - pd)))
    return out


METH_HAS_ACTION = {"ac21": True, "ac41": True, "ac43": True, "ac61": True, "ac65": True}


def sheet_complexity(model, ic, t0, dt, nsteps, n_dense=1 << 16):
    """Intrinsic complexity of the exact carrier: equi-affine length
    A = int |kappa|^(1/3) dq, Euclidean length, and the adaptivity gain
    G = (2 pi)^3 max|kappa| / A^3 (ratio of the max P1 chord defect of
    uniform-in-q sampling to that of affine-equidistributed sampling)."""
    q = np.arange(n_dense) / n_dense * TWO_PI
    x, p, S = oracle(model, ic, q, t0, dt, nsteps)
    sh = Sheet(q, x, p, S, ic.c)
    kap, spd = affine_curvature(sh)
    h = TWO_PI / n_dense
    A = float(np.sum(np.abs(kap) ** (1 / 3)) * h)
    return {"affine_length": A, "euclid_length": float(np.sum(spd) * h),
            "kappa_max": float(np.abs(kap).max()),
            "gain": float(TWO_PI**3 * np.abs(kap).max() / max(A, 1e-300) ** 3),
            "speed_max": float(spd.max()), "speed_min": float(spd.min())}


def run(model, ic, n, t0, dt, nsteps, method, adaptive, remap_every, out_every,
        mode="tri", mon=None, verbose=False):
    mon = {} if mon is None else mon
    q = np.arange(n) / n * TWO_PI
    x, p, S = ic(q)
    sh = Sheet(q, x, p, S, ic.c)
    t = t0
    hist = []
    nrem = 0
    for k in range(1, nsteps + 1):
        sh.x, sh.p, sh.S, t = integrate(model, sh.x, sh.p, sh.S, t, dt, 1)
        if adaptive and k % remap_every == 0 and k < nsteps:
            M = monitor_density(sh, mode=mode, **mon)
            qn = ot_labels(sh, M)
            xn, pn, Sn = evaluate(sh, qn, method)
            sh = Sheet(qn, xn, pn, Sn, sh.c)
            nrem += 1
        if k % out_every == 0 or k == nsteps:
            d = diagnostics(model, ic, sh, method, t0, dt, k, t)
            d["remaps"] = nrem
            hist.append(d)
            if verbose:
                print(f"   t={t:6.2f} rem={nrem:3d} en_max={d['en_max']:.2e} et_max={d['et_max']:.2e} "
                      f"per={d['period_err']:+.1e} nodal={d['nodal_err']:.1e} dqmin={d['dq_min']:.1e}", flush=True)
    return hist, sh


TESTS = {
    # smooth phase mixing in the circulating region (no hyperbolic point on the sheet)
    "twist": dict(model=Pendulum(1.0), ic=SineIC(2.6, 0.5), t0=0.0, dt=0.02, nsteps=2000),
    # cold collapse into the pendulum island; the hyperbolic point sits on the sheet
    "spiral": dict(model=Pendulum(1.0), ic=SineIC(0.05, 0.0), t0=0.0, dt=0.02, nsteps=600),
    # non-autonomous, chaotic separatrix layer
    "driven": dict(model=DrivenPendulum(1.0, 0.25, 1.7), ic=SineIC(0.0, 0.3), t0=0.0, dt=0.02, nsteps=600),
    # explicitly time-dependent comoving Hamiltonian (EdS), a from 1 to ~ 11.7
    "comoving": dict(model=ComovingEdS(0.25, 1.0), ic=SineIC(0.02, 0.0), t0=1.0, dt=0.04, nsteps=975),
}


def exp_history(test, n, methods, remap_every, n_out=20):
    cfg = TESTS[test]
    res = {"test": test, "n": n, "remap_every": remap_every, "dt": cfg["dt"], "nsteps": cfg["nsteps"], "runs": {}}
    out_every = max(cfg["nsteps"] // n_out, 1)
    for label, method, adaptive in methods:
        t_ = time.time()
        hist, _ = run(cfg["model"], cfg["ic"], n, cfg["t0"], cfg["dt"], cfg["nsteps"], method, adaptive,
                      remap_every, out_every)
        res["runs"][label] = hist
        f = hist[-1]
        print(f"[{test}] {label:14s} N={n} en_max={f['en_max']:.2e} en_l2m={f['en_l2m']:.2e} et_max={f['et_max']:.2e} "
              f"period_err={f['period_err']:+.2e} mass_err={f['mass_err']:+.1e} E_err={f['energy_err']:+.1e} "
              f"({time.time() - t_:.0f}s)", flush=True)
    res["complexity"] = [dict(t=cfg["t0"] + k * cfg["dt"], **sheet_complexity(cfg["model"], cfg["ic"], cfg["t0"], cfg["dt"], k))
                         for k in range(0, cfg["nsteps"] + 1, max(cfg["nsteps"] // 10, 1))]
    return res


def exp_convergence(test, ns, methods, remap_every):
    cfg = TESTS[test]
    res = {"test": test, "ns": ns, "remap_every": remap_every, "runs": {}}
    for label, method, adaptive in methods:
        rows = []
        for n in ns:
            hist, _ = run(cfg["model"], cfg["ic"], n, cfg["t0"], cfg["dt"], cfg["nsteps"], method, adaptive,
                          remap_every, cfg["nsteps"])
            rows.append(hist[-1])
        res["runs"][label] = rows
        en = [r["en_max"] for r in rows]
        orders = [np.log2(en[i] / en[i + 1]) for i in range(len(en) - 1)]
        print(f"[conv {test}] {label:14s} en_max=" + " ".join(f"{e:.2e}" for e in en) +
              "  orders=" + " ".join(f"{o:.2f}" for o in orders) +
              "  period_err=" + " ".join(f"{r['period_err']:+.1e}" for r in rows), flush=True)
    return res


def exp_accumulation(test, n, methods, intervals):
    cfg = TESTS[test]
    res = {"test": test, "n": n, "intervals": intervals, "runs": {}}
    for label, method, adaptive in methods:
        rows = []
        for iv in intervals:
            hist, _ = run(cfg["model"], cfg["ic"], n, cfg["t0"], cfg["dt"], cfg["nsteps"], method, adaptive,
                          iv, cfg["nsteps"])
            rows.append(hist[-1])
        res["runs"][label] = rows
        print(f"[accum {test}] {label:10s} remaps=" + " ".join(str(r["remaps"]) for r in rows) +
              "  period_err=" + " ".join(f"{r['period_err']:+.1e}" for r in rows) +
              "  en_max=" + " ".join(f"{r['en_max']:.1e}" for r in rows), flush=True)
    return res


def save(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w") as f:
        json.dump(obj, f, indent=1)


ALL_METHODS = [
    ("fixed-lin", "lin", False), ("fixed-cub", "cub", False), ("fixed-ac41", "ac41", False),
    ("adapt-lin", "lin", True), ("adapt-ac21", "ac21", True),
    ("adapt-cub", "cub", True), ("adapt-ac41", "ac41", True),
    ("adapt-quin", "quin", True), ("adapt-ac61", "ac61", True),
]

def exp_monitor_ablation():
    """AC(4,1) with the monitor driven by inscribed-triangle areas ("tri")
    versus by the measured edge action defects ("act")."""
    out = {}
    for test, n, every in (("spiral", 1024, 20), ("driven", 2048, 20)):
        cfg = TESTS[test]
        for mode in ("tri", "act"):
            hist, _ = run(cfg["model"], cfg["ic"], n, cfg["t0"], cfg["dt"], cfg["nsteps"], "ac41", True, every,
                          cfg["nsteps"], mode=mode)
            f = hist[-1]
            out[f"{test}_{mode}"] = {k: f[k] for k in ("en_max", "et_max", "period_err", "en_l2m", "dq_min")}
            print(f"[monitor ablation] {test} {mode}: en_max={f['en_max']:.2e} period_err={f['period_err']:+.1e}", flush=True)
    return out


JOBS = {
    "monitor_ablation": lambda: save("exp1d_monitor_ablation.json", exp_monitor_ablation()),
    "history_twist": lambda: save("exp1d_history_twist.json", exp_history("twist", 1024, ALL_METHODS, 50)),
    "history_spiral": lambda: save("exp1d_history_spiral.json", exp_history("spiral", 1024, ALL_METHODS, 20)),
    "history_driven": lambda: save("exp1d_history_driven.json", exp_history("driven", 2048, ALL_METHODS, 20)),
    "history_comoving": lambda: save("exp1d_history_comoving.json", exp_history("comoving", 1024, ALL_METHODS, 25)),
    "convergence_twist": lambda: save("exp1d_convergence_twist.json", exp_convergence(
        "twist", [256, 512, 1024, 2048, 4096], ALL_METHODS, 50)),
    "accumulation_twist": lambda: save("exp1d_accumulation_twist.json", exp_accumulation(
        "twist", 1024, [m for m in ALL_METHODS if m[2]], [1000, 500, 250, 100, 50, 25, 10, 5])),
}

if __name__ == "__main__":
    names = sys.argv[1:] or list(JOBS)
    for name in names:
        JOBS[name]()
