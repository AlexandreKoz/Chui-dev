"""Explicit field-convergence test (offline, controls only).

(a) Control self-convergence at fixed PM grid: 2048^2 vs 4096^2 particles
    (how good is the 'continuum-in-label' control the arms are scored against).
(b) PM field convergence: the SAME label lattice evolved on PM grids
    64, 128, 256, 512 (2048^2 particles each): how much the physical answer
    still moves with the force resolution.  Sheet-vs-N-body differences that
    are far below (b) are differences between two discretisations of a model
    whose own field error is larger -- they are reported, but (b) is the
    scale of physical significance.

Usage: python3 exp02_field_convergence.py        -> results/exp02_field_convergence.json
"""
from __future__ import annotations

import os

import numpy as np

import common

SNAP = (100, 150, 200, 250, 300)


def load(ng, n):
    p = os.path.join(common.DATA, f"vp_ref_ng{ng}_n{n}.npz")
    return np.load(p) if os.path.exists(p) else None


def rms(a):
    return float(np.sqrt(np.mean(np.sum(a ** 2, axis=-1))))


if __name__ == "__main__":
    out = dict(meta=common.metadata(__file__), self_convergence={}, field_convergence={}, folded={})
    for ng in (256, 512):
        a, b = load(ng, 2048), load(ng, 4096)
        if a is None or b is None:
            continue
        out["self_convergence"][str(ng)] = {
            str(k): dict(x_rms=rms(a[f"x{k}"] - b[f"x{k}"]), p_rms=rms(a[f"p{k}"] - b[f"p{k}"]),
                         x_max=float(np.abs(a[f"x{k}"] - b[f"x{k}"]).max()),
                         rho_rel_l2=float(np.linalg.norm(a[f"rho{k}"] - b[f"rho{k}"]) / np.linalg.norm(b[f"rho{k}"] - 1.0)))
            for k in SNAP}
        out["folded"][str(ng)] = {str(k): float(b[f"negdet{k}"]) for k in SNAP}
        out.setdefault("rho_max", {})[str(ng)] = {str(k): float(b[f"rho{k}"].max()) for k in SNAP}
    grids = [g for g in (64, 128, 256, 512) if load(g, 2048) is not None]
    for g0, g1 in zip(grids[:-1], grids[1:]):
        a, b = load(g0, 2048), load(g1, 2048)
        out["field_convergence"][f"{g0}->{g1}"] = {
            str(k): dict(x_rms=rms(a[f"x{k}"] - b[f"x{k}"]), p_rms=rms(a[f"p{k}"] - b[f"p{k}"]),
                         x_max=float(np.abs(a[f"x{k}"] - b[f"x{k}"]).max()),
                         x_rms_in_fine_cells=rms(a[f"x{k}"] - b[f"x{k}"]) * g1) for k in SNAP}
    common.save("exp02_field_convergence.json", out)
    for k, v in out["self_convergence"].items():
        print(f"PM {k}: control 2048^2 vs 4096^2 at t=3: x_rms {v['300']['x_rms']:.2e}  rho {v['300']['rho_rel_l2']:.2e}   folded {out['folded'][k]['300']:.3f}")
    for k, v in out["field_convergence"].items():
        print(f"PM {k}: x_rms at t=1.5/2/3 = {v['150']['x_rms']:.2e} / {v['200']['x_rms']:.2e} / {v['300']['x_rms']:.2e}")
