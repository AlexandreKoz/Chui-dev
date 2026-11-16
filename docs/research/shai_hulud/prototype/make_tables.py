"""Generate the LaTeX tables of the manuscript from ../results/*.json, so
that every number printed in the paper is traceable to a result file.

Usage:  python3 make_tables.py     (writes ../manuscript/tables_generated.tex)
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
OUT = os.path.join(HERE, "..", "manuscript", "tables_generated.tex")

LAB = {"fixed-lin": "fixed, linear", "fixed-cub": "fixed, cubic", "fixed-ac41": "fixed, AC(4,1)",
       "adapt-lin": "adaptive, linear", "adapt-ac21": "adaptive, AC(2,1)", "adapt-cub": "adaptive, cubic",
       "adapt-ac41": "adaptive, AC(4,1)", "adapt-quin": "adaptive, quintic", "adapt-ac61": "adaptive, AC(6,1)"}


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def sci(v, digits=1, signed=False):
    if v == 0:
        return "$0$"
    s = f"{v:+.{digits}e}" if signed else f"{abs(v):.{digits}e}"
    m, e = s.split("e")
    return f"${m}\\!\\times\\!10^{{{int(e)}}}$"


def table_history():
    out = []
    names = {"twist": "Twist ($N=1024$, $t=40$, 39 remaps)", "spiral": "Collapse ($N=1024$, $t=12$, 29 remaps)",
             "driven": "Driven ($N=2048$, $t=12$, 29 remaps)", "comoving": "Comoving EdS ($N=1024$, $a\\simeq 11.7$, 38 remaps)"}
    out.append("\\begin{tabular}{llrrrr}\\toprule")
    out.append("test & sampling, reconstruction & $e_\\perp^{\\max}$ & $e_\\parallel^{\\max}$ & $\\oint p\\,dx-c$ & $E_h-E$\\\\\\midrule")
    for test in ("twist", "spiral", "driven", "comoving"):
        d = load(f"exp1d_history_{test}.json")
        first = True
        for lab in LAB:
            f = d["runs"][lab][-1]
            out.append(f"{names[test] if first else ''} & {LAB[lab]} & {sci(f['en_max'])} & {sci(f['et_max'])} & "
                       f"{sci(f['period_err'], signed=True)} & {sci(f['energy_err'], signed=True)}\\\\")
            first = False
        out.append("\\midrule" if test != "comoving" else "\\bottomrule")
    out.append("\\end{tabular}")
    return "\n".join(out)


def table_convergence():
    d = load("exp1d_convergence_twist.json")
    ns = d["ns"]
    out = ["\\begin{tabular}{l" + "r" * len(ns) + "r}\\toprule",
           "sampling, reconstruction & " + " & ".join(f"$N={n}$" for n in ns) + " & order\\\\\\midrule"]
    for lab in LAB:
        e = [r["en_max"] for r in d["runs"][lab]]
        k = 3 if lab == "adapt-ac61" else 4  # AC(6,1) reaches the round-off conditioning floor at N = 4096
        order = np.log2(e[k - 1] / e[k])
        out.append(f"{LAB[lab]} & " + " & ".join(sci(v) for v in e) + f" & {order:.2f}\\\\")
    out.append("\\midrule")
    out.append("\\multicolumn{" + str(len(ns) + 2) + "}{l}{Poincar\\'e invariant error $\\oint p\\,dx-c$ of the carrier}\\\\\\midrule")
    for lab in LAB:
        e = [r["period_err"] for r in d["runs"][lab]]
        out.append(f"{LAB[lab]} & " + " & ".join(sci(v, signed=True) for v in e) + " & \\\\")
    out.append("\\bottomrule\\end{tabular}")
    return "\n".join(out)


def table_accumulation():
    d = load("exp1d_accumulation_twist.json")
    rem = [r["remaps"] for r in d["runs"]["adapt-ac41"]]
    out = ["\\begin{tabular}{l" + "r" * len(rem) + "}\\toprule",
           "remaps in $[0,40]$ & " + " & ".join(str(r) for r in rem) + "\\\\\\midrule",
           "\\multicolumn{" + str(len(rem) + 1) + "}{l}{$\\oint p\\,dx-c$}\\\\"]
    for lab in d["runs"]:
        out.append(f"{LAB[lab]} & " + " & ".join(sci(r["period_err"], signed=True) for r in d["runs"][lab]) + "\\\\")
    out.append("\\midrule\\multicolumn{" + str(len(rem) + 1) + "}{l}{$e_\\perp^{\\max}$ at $t=40$}\\\\")
    for lab in d["runs"]:
        out.append(f"{LAB[lab]} & " + " & ".join(sci(r["en_max"]) for r in d["runs"][lab]) + "\\\\")
    out.append("\\bottomrule\\end{tabular}")
    return "\n".join(out)


def table_complexity():
    out = ["\\begin{tabular}{lrrrrr}\\toprule",
           "test & $t$ & $\\mathcal A(t)$ & Euclidean length & $\\max|\\kappa|$ & gain $G(t)$\\\\\\midrule"]
    for test in ("twist", "spiral", "driven"):
        c = load(f"exp1d_history_{test}.json")["complexity"]
        rows = [r for r in c if r["affine_length"] > 0]
        pick = [rows[i] for i in sorted(set(np.linspace(0, len(rows) - 1, 4).astype(int)))]
        for i, r in enumerate(pick):
            out.append(f"{test if i == 0 else ''} & {r['t']:.1f} & {r['affine_length']:.2f} & {r['euclid_length']:.1f} & "
                       f"{sci(r['kappa_max'])} & {sci(r['gain'])}\\\\")
        out.append("\\midrule" if test != "driven" else "\\bottomrule")
    out.append("\\end{tabular}")
    return "\n".join(out)


def table_2d():
    d = load("exp2d.json")["runs"]
    name = {"lin": "bilinear", "cub": "bicubic spline", "quin": "biquintic spline", "front3": "front-3 ($C^1$ Hermite $W$)",
            "front5": "front-5 ($C^2$ Hermite $W$)"}
    out = ["\\begin{tabular}{lrrrrrr}\\toprule",
           " & \\multicolumn{3}{c}{$\\max|p-p^{\\rm oracle}|$} & \\multicolumn{3}{c}{$\\max_\\gamma|\\oint_\\gamma p\\cdot dx|$}\\\\",
           "reconstruction & $n=32$ & $n=64$ & $n=128$ & $n=32$ & $n=64$ & $n=128$\\\\\\midrule",
           "\\multicolumn{7}{l}{adaptive (3 optimal-transport remaps)}\\\\"]
    for rec in ("lin", "cub", "quin", "front3", "front5"):
        out.append(name[rec] + " & " + " & ".join(sci(d[f"adapt-{rec}-n{n}"]["p_err_max"]) for n in (32, 64, 128)) + " & " +
                   " & ".join(sci(d[f"adapt-{rec}-n{n}"]["circulation_max"]) for n in (32, 64, 128)) + "\\\\")
    out.append("\\midrule\\multicolumn{7}{l}{fixed Lagrangian mesh (exact nodal data; reconstruction defect only)}\\\\")
    for rec in ("lin", "cub", "quin", "front3", "front5"):
        out.append(name[rec] + " & -- & -- & -- & " +
                   " & ".join(sci(d[f"fixed-{rec}-n{n}"]["circulation_max"]) for n in (32, 64, 128)) + "\\\\")
    out.append("\\bottomrule\\end{tabular}")
    ot = d["adapt-front5-n128"]["ot"]
    out2 = ["\\begin{tabular}{rrrrrrr}\\toprule",
            "remap at step & $\\max M/\\min M$ & Newton its. & MA residual & $\\lambda_{\\min}(I+D^2\\varphi)$ & $\\int|\\nabla\\varphi|^2$ & point-location residual\\\\\\midrule"]
    for o in ot:
        out2.append(f"{o['step']} & {o['monitor_contrast']:.2f} & {o['ma_newton_iters']} & {sci(o['ma_residual'])} & "
                    f"{o['hess_eig_min']:.3f} & {sci(o['w2_cost'])} & {sci(o['locate_err'])}\\\\")
    out2.append("\\bottomrule\\end{tabular}")
    m = max(abs(v["mass_err"]) for v in d.values())
    r = max(v["cell_mass_ratio"] for v in d.values())
    return "\n".join(out), "\n".join(out2), m, r


def table_verification():
    v = load("verification.json")
    rig = v["rigidity_kernel_dimension_d2_equals_3n"]
    out = ["\\begin{tabular}{llrrr}\\toprule", "$d$ & vertex positions & $n$ & unknowns $(p,\\Sigma,c)$ & kernel dimension\\\\\\midrule"]
    for key, lab in (("d2_n5", "2"), ("d2_n6", "2"), ("d2_n8", "2"), ("d3_n4", "3"), ("d3_n5", "3")):
        n = key.split("n")[1]
        for kind, kn in (("regular", "regular lattice"), ("jittered", "generic (jittered)")):
            k, nunk, _ = rig[f"{key}_{kind}"]
            out.append(f"{lab} & {kn} & {n} & {nunk} & {k}\\\\")
    out.append("\\bottomrule\\end{tabular}")
    return "\n".join(out), v


if __name__ == "__main__":
    t2d, tot, mass2, ratio2 = table_2d()
    trig, v = table_verification()
    curl = v["cic_fd_force_is_not_a_gradient"]
    vn = v["von_neumann_spectral_radius"]
    parts = [
        "% Generated by prototype/make_tables.py from results/*.json -- do not edit by hand.",
        "\\newcommand{\\TableHistory}{" + table_history() + "}",
        "\\newcommand{\\TableConvergence}{" + table_convergence() + "}",
        "\\newcommand{\\TableAccumulation}{" + table_accumulation() + "}",
        "\\newcommand{\\TableComplexity}{" + table_complexity() + "}",
        "\\newcommand{\\TableTwoD}{" + t2d + "}",
        "\\newcommand{\\TableOT}{" + tot + "}",
        "\\newcommand{\\TableRigidity}{" + trig + "}",
        "\\newcommand{\\TwoDMassErr}{" + sci(mass2) + "}",
        "\\newcommand{\\TwoDMassRatio}{" + f"{ratio2:.1f}" + "}",
        "\\newcommand{\\CurlRatios}{" + ", ".join(f"{curl[k]:.2f}" for k in ("n16_cic_fd_curl_over_grad", "n32_cic_fd_curl_over_grad", "n64_cic_fd_curl_over_grad")) + "}",
        "\\newcommand{\\VNRadius}{" + f"{max(vn['ac21_max_spectral_radius'], vn['ac41_max_spectral_radius'], vn['ac61_max_spectral_radius']) - 1:.1e}" + "}",
    ]
    with open(OUT, "w") as f:
        f.write("\n".join(parts) + "\n")
    print("wrote", os.path.normpath(OUT))
