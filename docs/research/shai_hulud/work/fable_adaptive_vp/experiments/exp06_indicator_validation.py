"""Offline validation of the state-only indicators against the true error map.

A uniform quintic sheet (n = N0, default NG/2; deposit lattice 4 NG) is evolved; at t = 1.5, 2.0,
2.5, 3.0 the indicators are computed from the state exactly as the controller
would, and -- separately, offline -- the carrier error against the 4096^2
control is binned per block.  The control is never seen by the indicator
code.  Answers: where does the sheet error live (folded labels or
everywhere)?  do the indicators find it?  by how much do they under-estimate
the accumulated error?

Usage: python3 exp06_indicator_validation.py [NG [N0]]  -> results/exp06_indicator_validation_ng<NG>_n<N0>.json
"""
from __future__ import annotations

import sys

import numpy as np

import bench
import common
from shadapt.controller import Controller
from shadapt.sheet import AdaptiveSheet
from shadapt.solver import AdaptiveVP, Quadrature
from shsheet.vp import GrowingModeIC

if __name__ == "__main__":
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    n0 = int(sys.argv[2]) if len(sys.argv) > 2 else ng // 2
    sh = AdaptiveSheet.from_ic(GrowingModeIC(), n0, B=8, lmax=1, k=5)
    run = AdaptiveVP(sh, ng, Quadrature("lattice", 4 * ng // n0))
    ctrl = Controller(indicator="pred", tau=0.6, tau_p=0.3)       # used only to evaluate indicators
    out = dict(meta=common.metadata(__file__, ng=ng, n0=n0), snaps={})
    state_ind = {}
    for k in range(1, bench.NSTEPS + 1):
        run.step(bench.DT)
        if k in (150, 200, 250, 300):
            run.sync()
            ind = ctrl.indicators(run)                             # state only
            state_ind[k] = (ind, sh.car.leafI.copy(), sh.car.leafJ.copy(), sh.evaluate(
                np.stack(np.meshgrid(np.arange(256) / 256, np.arange(256) / 256, indexing="ij"), -1).reshape(-1, 2), "x")[0])
    # ---- offline part: the control is opened only now ----
    ref = np.load(bench.ref_path(ng))
    nb = n0 // 8
    cb = max(1, nb // 32)                                   # coarsen the stored maps to at most 32 x 32
    for k, (ind, LI, LJ, x) in state_ind.items():
        xr = ref[f"x{k}"]
        st = xr.shape[0] // 256
        err = np.sqrt(np.sum((x - xr[::st, ::st].reshape(-1, 2)) ** 2, axis=-1)).reshape(256, 256) * ng     # PM cells
        if nb > 32:                                         # blocks finer than the 256^2 evaluation lattice allows: use 32 x 32 super-blocks
            sb = 32
        else:
            sb = nb
        emax = err.reshape(sb, 256 // sb, sb, 256 // sb).max(axis=(1, 3))
        erms = np.sqrt((err ** 2).reshape(sb, 256 // sb, sb, 256 // sb).mean(axis=(1, 3)))
        nr = xr.shape[0]
        a = np.roll(xr, -1, 0) - xr
        a[-1, :, 0] += 1.0
        b = np.roll(xr, -1, 1) - xr
        b[:, -1, 1] += 1.0
        det = (a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]) * nr * nr
        fold = (det < 0).reshape(sb, nr // sb, sb, nr // sb).mean(axis=(1, 3))

        def tomap(v):
            m = np.zeros((nb, nb))
            m[LI, LJ] = v
            return m.reshape(sb, nb // sb, sb, nb // sb).max(axis=(1, 3))      # super-block maximum
        e = emax.ravel()
        rec = dict(t=k * bench.DT, err_max_cells=float(err.max()), err_rms_cells=float(np.sqrt(np.mean(err ** 2))),
                   folded_label_fraction=float((det < 0).mean()), blocks_with_fold=float((fold > 0).mean()),
                   rms_err_in_folded_blocks=float(np.sqrt(np.mean(erms[fold > 0] ** 2))) if (fold > 0).any() else None,
                   rms_err_in_unfolded_blocks=float(np.sqrt(np.mean(erms[fold == 0] ** 2))) if (fold == 0).any() else None,
                   share_of_squared_error_in_folded_blocks=float((erms[fold > 0] ** 2).sum() / (erms ** 2).sum()) if (fold > 0).any() else 0.0,
                   block_fraction_err_gt={f"{th:g}": float((e > th).mean()) for th in (1e-3, 1e-4, 1e-5, 1e-6)},
                   err_block_max=emax.tolist(), fold_block=fold.tolist())
        rec["extent_block"] = tomap(ind["extent"]).tolist()
        rec["floor_would_deny_fraction_of_top_error_blocks"] = float(np.mean(tomap(ind["extent"]).ravel()[np.argsort(-e)[: max(1, len(e) // 20)]] <= 0.5))
        for name in ("eta_x", "eta_p", "eta_pred", "eta_d6"):
            v = tomap(ind[name]).ravel()
            rec[name + "_block"] = tomap(ind[name]).tolist()
            rec[name + "_logcorr_with_err"] = float(np.corrcoef(np.log(e + 1e-300), np.log(v + 1e-300))[0, 1])
            rec[name + "_over_err_median"] = float(np.median(v / e))
            rec[name + "_over_err_p10_p90"] = [float(np.percentile(v / e, 10)), float(np.percentile(v / e, 90))]
        out["snaps"][str(k)] = rec
        print(f"t={rec['t']:.1f}: err rms {rec['err_rms_cells']:.2e} max {rec['err_max_cells']:.2e} PM cells; folded labels {rec['folded_label_fraction']:.3f}, "
              f"blocks with fold {rec['blocks_with_fold']:.2f}, share of err^2 in folded blocks {rec['share_of_squared_error_in_folded_blocks']:.3f}; "
              f"corr(log eta_d6, log err) {rec['eta_d6_logcorr_with_err']:.2f}, eta_d6/err median {rec['eta_d6_over_err_median']:.3f}; "
              f"corr pred {rec['eta_pred_logcorr_with_err']:.2f}, pred/err median {rec['eta_pred_over_err_median']:.3f}", flush=True)
    common.save(f"exp06_indicator_validation_ng{ng}_n{n0}.json", out)
