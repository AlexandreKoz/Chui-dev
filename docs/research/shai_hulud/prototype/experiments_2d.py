"""d = 2 experiments (four-dimensional phase space, single-chart regime).

(a) Monge-Ampere optimal transport of the *sampling* on label space:
    convergence, convexity, exact total mass, variable cell masses.
(b) OT-guided remapping of a cold sheet under an external Hamiltonian, with
    five reconstructions (bilinear, bicubic and biquintic splines, and the
    C^1 / C^2 generating-function "front" elements).  Reported: error against the oracle and the Kelvin /
    Poincare circulation of each method's own carrier on fixed loops.

Usage:  python3 experiments_2d.py
Writes ../results/exp2d.json
"""
from __future__ import annotations

import json
import os

import numpy as np

from sheet2d import (Carrier2, Potential2, Sheet2, grid, kdk, label_mass, loop_circulation,
                     monitor_on_labels, oracle2, remap2, solve_ot)

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
DT, NSTEPS = 0.01, 90
REMAP_STEPS = (30, 50, 70)


def loops(seed=3, count=24):
    rng = np.random.default_rng(seed)
    return [(rng.random(2), 0.05 + 0.2 * rng.random()) for _ in range(count)]


def final_diagnostics(pot, sh, method, nsteps):
    xo, po, So = oracle2(pot, sh.q, DT, nsteps)
    car = Carrier2(sh, method)
    circ = np.array([loop_circulation(car, c, r) for c, r in loops()])
    tot, amin, amax = label_mass(sh.q)
    return {
        "x_err_max": float(np.abs(sh.x - xo).max()), "p_err_max": float(np.abs(sh.p - po).max()),
        "p_err_l2": float(np.sqrt(np.mean((sh.p - po) ** 2))),
        "S_err_max": float(np.abs((sh.S - So) - (sh.S - So).mean()).max()),
        "circulation_max": float(np.abs(circ).max()), "circulation_rms": float(np.sqrt(np.mean(circ**2))),
        "mass_err": tot - 1.0, "cell_mass_ratio": amax / amin,
    }


def run(n, method, adaptive, verbose=False):
    pot = Potential2()
    q = grid(n)
    sh = Sheet2(q.copy(), q.copy(), np.zeros_like(q), np.zeros((n, n)))
    ot_info = []
    done = 0
    for step in list(REMAP_STEPS) + [NSTEPS]:
        sh.x, sh.p, sh.S = kdk(pot, sh.x, sh.p, sh.S, DT, step - done)
        done = step
        if adaptive and step < NSTEPS:
            M, inv_err = monitor_on_labels(sh, floor=0.15)
            ot = solve_ot(M, n)
            sh, loc_err = remap2(sh, method, ot["T"], ot["phi"])
            ot_info.append({"step": step, "monitor_contrast": float(M.max() / M.min()),
                            "ma_residual": ot["residual_history"][-1], "ma_newton_iters": len(ot["residual_history"]) - 1,
                            "hess_eig_min": ot["hess_eig_min"], "det_min": ot["det_min"], "det_max": ot["det_max"],
                            "w2_cost": ot["cost"], "locate_err": loc_err, "label_inverse_err": inv_err})
            if verbose:
                print("    ", ot_info[-1])
    d = final_diagnostics(pot, sh, method, NSTEPS)
    d["ot"] = ot_info
    return d


if __name__ == "__main__":
    res = {"dt": DT, "nsteps": NSTEPS, "remap_steps": REMAP_STEPS, "runs": {}}
    for n in (32, 64, 128):
        for adaptive in (False, True):
            for method in ("lin", "cub", "front3", "quin", "front5"):
                key = f"{'adapt' if adaptive else 'fixed'}-{method}-n{n}"
                d = run(n, method, adaptive)
                res["runs"][key] = d
                print(f"{key:20s} p_err_max={d['p_err_max']:.2e} x_err_max={d['x_err_max']:.2e} "
                      f"circ_max={d['circulation_max']:.2e} circ_rms={d['circulation_rms']:.2e} "
                      f"mass_err={d['mass_err']:+.1e} cell_mass_ratio={d['cell_mass_ratio']:.1f}", flush=True)
                if adaptive and method == "front5":
                    for o in d["ot"]:
                        print(f"      OT step {o['step']}: contrast={o['monitor_contrast']:.1f} newton={o['ma_newton_iters']} "
                              f"res={o['ma_residual']:.1e} eigmin={o['hess_eig_min']:.3f} W2^2={o['w2_cost']:.2e}")
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "exp2d.json"), "w") as f:
        json.dump(res, f, indent=1)
