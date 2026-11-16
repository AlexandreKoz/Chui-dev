"""Figures for the SHAI-HULUD manuscript, drawn from ../results/*.json.

Usage:  python3 make_figures.py      (writes ../results/figures/*.pdf)

Colour is by reconstruction (fixed categorical order, validated palette),
line style by sampling strategy (solid = adaptive, dashed = fixed Lagrangian
mesh), marker shape repeats the reconstruction; every series is also directly
labelled and all plotted numbers are tabulated in the manuscript.
"""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(RES, "figures")

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
COL = {"lin": "#2a78d6", "cub": "#eb6834", "ac41": "#1baf7a", "quin": "#eda100", "ac61": "#e87ba4",
       "ac21": "#6250d6", "front3": "#1baf7a", "front5": "#e87ba4"}
MRK = {"lin": "o", "cub": "s", "ac41": "D", "quin": "^", "ac61": "v", "ac21": "P", "front3": "D", "front5": "v"}
NAME = {"lin": "linear", "cub": "cubic", "ac41": "AC(4,1)", "quin": "quintic", "ac61": "AC(6,1)", "ac21": "AC(2,1)",
        "front3": "front-3", "front5": "front-5"}

plt.rcParams.update({
    "font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.grid": True, "axes.axisbelow": True, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "lines.linewidth": 1.6, "lines.markersize": 4.5, "pdf.fonttype": 42,
})


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def style(label):
    kind, rec = label.split("-")
    return dict(color=COL[rec], marker=MRK[rec], linestyle="-" if kind == "adapt" else (0, (4, 2)),
                markerfacecolor=COL[rec] if kind == "adapt" else SURF, markeredgecolor=COL[rec], markeredgewidth=1.0)


def end_label(ax, x, y, text, dy=0.0):
    ax.annotate(text, (x, y), xytext=(5, dy), textcoords="offset points", va="center", fontsize=7.5, color=INK)


def fig_history():
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2), constrained_layout=True)
    shown = ["fixed-lin", "fixed-ac41", "adapt-lin", "adapt-cub", "adapt-ac41"]
    offs = {"spiral": {"fixed-lin": 4, "fixed-ac41": -5}, "driven": {"fixed-lin": 4, "fixed-ac41": -5}}
    for ax, test, title in zip(axes[:2], ("spiral", "driven"),
                               ("Cold collapse in a pendulum well (N = 1024)", "Driven, chaotic pendulum (N = 2048)")):
        d = load(f"exp1d_history_{test}.json")
        for lab in shown:
            h = d["runs"][lab]
            t = [r["t"] for r in h]
            e = [max(r["en_max"], 1e-16) for r in h]
            ax.semilogy(t, e, markevery=3, **style(lab))
            kind, rec = lab.split("-")
            end_label(ax, t[-1], e[-1], f"{'adaptive' if kind == 'adapt' else 'fixed'} {NAME[rec]}", offs[test].get(lab, 0))
        ax.set_xlim(0, t[-1] * 1.42)
        ax.set_xlabel("time t")
        ax.set_title(title, fontsize=8.5, loc="left")
    axes[0].set_ylabel("max normal distance to the exact sheet")
    ax = axes[2]
    for test, ls, mk in (("spiral", "-", "o"), ("driven", (0, (4, 2)), "s")):
        c = load(f"exp1d_history_{test}.json")["complexity"]
        t = [r["t"] for r in c if r["gain"] == r["gain"]]
        g = [r["gain"] for r in c if r["gain"] == r["gain"]]
        ax.semilogy(t, g, color=INK2, linestyle=ls, marker=mk, markerfacecolor=SURF)
        end_label(ax, t[-1], g[-1], "collapse" if test == "spiral" else "driven", 5 if test == "spiral" else -5)
    ax.set_xlim(0, 15.5)
    ax.set_xlabel("time t")
    ax.set_ylabel("adaptivity gain G(t) of the exact sheet")
    ax.set_title("Intrinsic gain of equi-affine sampling", fontsize=8.5, loc="left")
    fig.savefig(os.path.join(FIG, "fig_history.pdf"))
    plt.close(fig)


def fig_convergence():
    d = load("exp1d_convergence_twist.json")
    ns = np.array(d["ns"])
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), constrained_layout=True)
    ax = axes[0]
    for lab in ("adapt-lin", "adapt-cub", "adapt-quin", "adapt-ac41", "adapt-ac61", "fixed-cub", "fixed-ac41"):
        e = np.array([r["en_max"] for r in d["runs"][lab]])
        ax.loglog(ns, e, **style(lab))
        kind, rec = lab.split("-")
        dy = {"adapt-quin": 13, "adapt-ac41": 4, "adapt-ac61": -5, "fixed-ac41": -6}.get(lab, 0)
        end_label(ax, ns[-1], e[-1], f"{'adaptive' if kind == 'adapt' else 'fixed'} {NAME[rec]}", dy)
    for order, c0 in ((2, 3e-1), (4, 3e-3), (5, 4e-5)):
        ax.loglog(ns, c0 * (ns / ns[0]) ** (-float(order)), color=GRID, linewidth=1.0, zorder=0)
        ax.annotate(f"$N^{{-{order}}}$", (ns[0], c0), xytext=(-2, 5), textcoords="offset points", fontsize=7, color=INK2)
    ax.set_xlim(ns[0] * 0.8, ns[-1] * 4.2)
    ax.set_xlabel("vertices N")
    ax.set_ylabel("max normal distance to the exact sheet at t = 40")
    ax.set_title("Accuracy (twist test, 39 remaps)", fontsize=8.5, loc="left")
    ax = axes[1]
    for lab in ("adapt-lin", "adapt-cub", "adapt-quin", "adapt-ac41", "adapt-ac61"):
        e = np.array([max(abs(r["period_err"]), 1e-16) for r in d["runs"][lab]])
        ax.loglog(ns, e, **style(lab))
        rec = lab.split("-")[1]
        end_label(ax, ns[-1], e[-1], NAME[rec], 5 if rec == "ac61" else (-5 if rec == "ac41" else 0))
    ax.axhline(1e-16, color=INK2, linewidth=0.6)
    ax.annotate("plotted floor $10^{-16}$ (values at or below round-off)", (ns[0], 1e-16), xytext=(0, 4),
                textcoords="offset points", fontsize=7, color=INK2)
    ax.set_xlim(ns[0] * 0.8, ns[-1] * 2.6)
    ax.set_xlabel("vertices N")
    ax.set_ylabel(r"$|\oint p\,dx - c|$ of the method's own carrier")
    ax.set_title("Poincaré invariant (exact value c = 5.2π)", fontsize=8.5, loc="left")
    fig.savefig(os.path.join(FIG, "fig_convergence.pdf"))
    plt.close(fig)


def fig_accumulation():
    d = load("exp1d_accumulation_twist.json")
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4), constrained_layout=True)
    for ax, key, yl, title in zip(
            axes, ("period_err", "en_max"),
            (r"$|\oint p\,dx - c|$ of the method's own carrier", "max normal distance to the exact sheet at t = 40"),
            ("Poincaré invariant vs number of remaps", "Accuracy vs number of remaps")):
        for lab in ("adapt-lin", "adapt-cub", "adapt-quin", "adapt-ac41", "adapt-ac61"):
            rows = d["runs"][lab]
            nr = np.array([r["remaps"] for r in rows])
            e = np.array([max(abs(r[key]), 1e-16) for r in rows])
            ax.loglog(nr, e, **style(lab))
            rec = lab.split("-")[1]
            dy = 0
            if key == "period_err":
                dy = 5 if rec == "ac61" else (-5 if rec == "ac41" else 0)
            end_label(ax, nr[-1], e[-1], NAME[rec], dy)
        ax.set_xlim(0.8, 1800)
        ax.set_xlabel("number of remaps in [0, 40]")
        ax.set_ylabel(yl)
        ax.set_title(title + " (N = 1024)", fontsize=8.5, loc="left")
    axes[0].axhline(1e-16, color=INK2, linewidth=0.6)
    axes[0].annotate("plotted floor $10^{-16}$", (1, 1e-16), xytext=(0, 4), textcoords="offset points", fontsize=7, color=INK2)
    fig.savefig(os.path.join(FIG, "fig_accumulation.pdf"))
    plt.close(fig)


def fig_sheet():
    """Phase portrait of the collapse test: exact sheet, fixed Lagrangian
    vertices and OT-equidistributed vertices at equal count."""
    from sheet1d import TWO_PI, Sheet, monitor_density, ot_labels, oracle
    from experiments_1d import TESTS, run

    cfg = TESTS["spiral"]
    n, nst = 256, 450
    hist, sh = run(cfg["model"], cfg["ic"], n, cfg["t0"], cfg["dt"], nst, "ac41", True, 20, nst)
    # exact sheet, sampled densely in the adaptive labels so the stretched part is drawn
    M = monitor_density(Sheet(sh.q, sh.x, sh.p, sh.S, sh.c))
    qfine = ot_labels(sh, M, 40000)
    xd, pd, _ = oracle(cfg["model"], cfg["ic"], qfine, cfg["t0"], cfg["dt"], nst)
    qf = np.arange(n) / n * TWO_PI
    xf, pf, _ = oracle(cfg["model"], cfg["ic"], qf, cfg["t0"], cfg["dt"], nst)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), constrained_layout=True, sharey=True)
    wrap = lambda x: np.mod(x + np.pi, TWO_PI) - np.pi
    for ax, xs, ps, rec, title in ((axes[0], xf, pf, "cub", f"Fixed Lagrangian labels (N = {n})"),
                                   (axes[1], sh.x, sh.p, "lin", f"Equidistributed labels after remapping (N = {n})")):
        xw = wrap(xd)
        xw[np.abs(np.diff(xw, prepend=xw[0])) > np.pi] = np.nan
        ax.plot(xw, pd, "-", color="#a5a49e", linewidth=0.5, zorder=1)
        ax.plot(wrap(xs), ps, "o", color=COL[rec], markersize=2.6, markeredgecolor=SURF, markeredgewidth=0.3, zorder=2)
        ax.set_xlabel("x")
        ax.set_title(title, fontsize=8.5, loc="left")
        ax.grid(False)
    axes[0].set_ylabel("p")
    axes[0].annotate("exact sheet (grey), vertices (colour); t = 9", (0.02, 0.04), xycoords="axes fraction", fontsize=7.5, color=INK2)
    fig.savefig(os.path.join(FIG, "fig_sheet.pdf"), dpi=300)
    plt.close(fig)


def fig_2d():
    from sheet2d import Potential2, Sheet2, grid, kdk, monitor_on_labels, solve_ot

    d = load("exp2d.json")
    n = 64
    q = grid(n)
    sh = Sheet2(q.copy(), q.copy(), np.zeros_like(q), np.zeros((n, n)))
    sh.x, sh.p, sh.S = kdk(Potential2(), sh.x, sh.p, sh.S, 0.01, 70)
    M, _ = monitor_on_labels(sh, floor=0.15)
    T = solve_ot(M, n)["T"]
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.5), constrained_layout=True)
    ax = axes[0]
    im = ax.imshow(M.T, origin="lower", extent=(0, 1, 0, 1), cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
        "seq", ["#eef4fc", "#2a78d6", "#123a70"]), interpolation="bilinear")
    Te = np.concatenate([T, T[:1] + np.array([1.0, 0.0])], axis=0)
    Te = np.concatenate([Te, Te[:, :1] + np.array([0.0, 1.0])], axis=1)
    for k in range(0, n + 1, 2):
        ax.plot(Te[k, :, 0], Te[k, :, 1], color=INK, linewidth=0.35)
        ax.plot(Te[:, k, 0], Te[:, k, 1], color=INK, linewidth=0.35)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("label $q_1$")
    ax.set_ylabel("label $q_2$")
    ax.set_title("Monitor M(q) and the optimal-transport mesh", fontsize=8.5, loc="left")
    fig.colorbar(im, ax=ax, shrink=0.8, label="M (unit mean)")
    ns = np.array([32, 64, 128])
    for ax, key, yl, title in ((axes[1], "p_err_max", "max momentum error vs oracle", "Accuracy after 3 OT remaps"),
                               (axes[2], "circulation_max", r"max $|\oint p\cdot dx|$ over 24 loops", "Kelvin/Poincaré circulation of the carrier")):
        for rec in ("lin", "cub", "quin", "front3", "front5"):
            e = np.array([max(d["runs"][f"adapt-{rec}-n{m}"][key], 1e-16) for m in ns])
            st = dict(color=COL[rec], marker=MRK[rec], markeredgecolor=COL[rec])
            if rec.startswith("front"):
                st.update(linestyle=(0, (1, 1.5)))
            ax.loglog(ns, e, **st)
            end_label(ax, ns[-1], e[-1], NAME[rec])
        ax.set_xlim(27, 330)
        ax.set_xticks(ns)
        ax.set_xticklabels([str(m) for m in ns])
        ax.minorticks_off()
        ax.set_xlabel("nodes per dimension n")
        ax.set_ylabel(yl)
        ax.set_title(title, fontsize=8.5, loc="left")
    fig.savefig(os.path.join(FIG, "fig_2d.pdf"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_history()
    fig_convergence()
    fig_accumulation()
    fig_sheet()
    fig_2d()
    print("figures written to", os.path.normpath(FIG))
