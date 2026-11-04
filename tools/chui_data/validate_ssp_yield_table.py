#!/usr/bin/env python3
"""Independent integrity check of a chui_ssp_yield_table_v1 file against its
NuPyCEE sources (D1, 2026-11-04 data closure).

Re-derives every quantity from the source files without the converter's
code path and checks:
  * header provenance (magic, citation, licence, tolerance, source sha256);
  * element order equals the CHUI catalog (H He C N O Ne Mg Si S Ca Fe);
  * channel assignment: agb masses < cc_mass_min <= cc masses; rectangular grid;
  * ejecta == M - M_final per model and == sum of gross isotope yields
    within --ejecta-tolerance (default 5e-4 relative);
  * per element: gross tracked ejecta = X0_e * m_ej + net_e >= -1e-12 M_sun
    (a negative gross ejecta is physically inconsistent and REJECTED);
  * |sum net_e| / m_ej <= declared net_yield_sum_tolerance (CHUI loader rule)
    and equals minus the untracked net residual (mass conservation);
  * Ia: element masses >= 0, sum == ejecta within 1e-3 (CHUI loader rule),
    tracked + untracked == total W7 ejecta;
  * mass and metallicity coverage reported.
Exit status 0 only if every check passes. Pure Python; reads text only.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

CHUI = ["H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def element(label: str) -> str:
    m = re.fullmatch(r"([A-Z][a-z]?)-?(\d+)", label.strip())
    if not m:
        raise SystemExit(f"bad isotope {label}")
    return m.group(1)


def read_massive(path: Path):
    models, cur = [], None
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("H Table:"):
            m = re.search(r"M\s*=\s*([0-9.eE+-]+)\s*,\s*Z\s*=\s*([0-9.eE+-]+)", line)
            cur = {"M": float(m.group(1)), "Z": float(m.group(2)), "Mf": None, "iso": []}
            models.append(cur)
        elif line.startswith("H Mfinal:") and cur is not None:
            cur["Mf"] = float(line.split(":", 1)[1])
        elif line.startswith("&") and cur is not None and not line.startswith("&Isotopes"):
            f = [x.strip() for x in line.strip().strip("&").split("&")]
            cur["iso"].append((f[0], float(f[1]), float(f[2])))
    return models


def read_table(path: Path):
    header, channels, ch = {}, {}, None
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("#"):
            body = line[1:].strip()
            if "=" in body and not body.startswith("model"):
                k, v = body.split("=", 1)
                header[k.strip()] = v.strip()
            elif body == "chui_ssp_yield_table_v1":
                header["magic"] = True
            continue
        tok = line.split()
        if not tok:
            continue
        if tok[0] == "elements":
            header["elements"] = tok[1:]
        elif tok[0] == "channel":
            ch = tok[1]
            channels[ch] = {"metallicity": [], "mass": [], "rows": {}}
        elif tok[0] in ("metallicity", "mass"):
            channels[ch][tok[0]] = [float(x) for x in tok[1:]]
        elif tok[0] == "row":
            if ch == "ia":
                channels[ch]["rows"][(int(tok[1]), 0)] = [float(x) for x in tok[2:]]
            else:
                channels[ch]["rows"][(int(tok[1]), int(tok[2]))] = [float(x) for x in tok[3:]]
    return header, channels


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("table", type=Path)
    p.add_argument("--massive", required=True, type=Path)
    p.add_argument("--ia", required=True, type=Path)
    p.add_argument("--cc-mass-min", type=float, default=8.0)
    p.add_argument("--ejecta-tolerance", type=float, default=5e-4)
    a = p.parse_args()
    errors = []
    header, ch = read_table(a.table)
    if not header.get("magic"):
        errors.append("missing magic")
    for key in ("citation", "license", "net_yield_sum_tolerance"):
        if not header.get(key):
            errors.append(f"missing header {key}")
    if header.get("elements") != CHUI:
        errors.append(f"element order {header.get('elements')} != CHUI catalog")
    if f"sha256={sha(a.massive)}" not in header.get("source_massive", ""):
        errors.append("source_massive sha256 mismatch")
    if f"sha256={sha(a.ia)}" not in header.get("source_ia", ""):
        errors.append("source_ia sha256 mismatch")
    tol = float(header.get("net_yield_sum_tolerance", "0"))
    if any(m >= a.cc_mass_min for m in ch["agb"]["mass"]) or any(m < a.cc_mass_min for m in ch["cc"]["mass"]):
        errors.append("channel assignment violates cc_mass_min")
    worst_ej, worst_net, min_gross, worst_closure = 0.0, 0.0, float("inf"), 0.0
    for model in read_massive(a.massive):
        c = "agb" if model["M"] < a.cc_mass_min else "cc"
        iz = ch[c]["metallicity"].index(model["Z"])
        im = ch[c]["mass"].index(model["M"])
        row = ch[c]["rows"][(iz, im)]
        ej, net = row[0], row[1:]
        ej_src = model["M"] - model["Mf"]
        gross_sum = sum(y for _, y, _ in model["iso"])
        x0 = defaultdict(float)
        untracked = 0.0
        for iso, y, x in model["iso"]:
            e = element(iso)
            x0[e] += x
            if e not in CHUI:
                untracked += y - x * ej_src
        worst_ej = max(worst_ej, abs(ej - ej_src) / ej_src, abs(gross_sum - ej_src) / ej_src)
        s = sum(net)
        worst_net = max(worst_net, abs(s) / ej)
        if abs(s) / ej > tol:
            errors.append(f"M={model['M']} Z={model['Z']}: net sum {s:.3e} above tolerance")
        # Exact identity: tracked net + untracked net = sum(gross) - ej * sum(X0).
        # Its right side is the SOURCE non-closure (rounded yields, sum X0 != 1).
        closure = gross_sum - ej_src * sum(x for _, _, x in model["iso"])
        worst_closure = max(worst_closure, abs(closure) / ej_src)
        if abs(s + untracked - closure) > 1e-6 * ej + 1e-8:
            errors.append(f"M={model['M']} Z={model['Z']}: net identity violated ({s + untracked:.4e} vs {closure:.4e})")
        for k, e in enumerate(CHUI):
            gross = x0[e] * ej + net[k]
            min_gross = min(min_gross, gross)
            if gross < -1e-12:
                errors.append(f"M={model['M']} Z={model['Z']} {e}: negative gross ejecta {gross:.3e}")
    if worst_ej > a.ejecta_tolerance:
        errors.append(f"ejecta inconsistency {worst_ej:.3e} > {a.ejecta_tolerance}")
    ia = ch["ia"]["rows"][(0, 0)]
    ia_total = 0.0
    for line in a.ia.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("&") and not line.startswith("&Isotopes"):
            ia_total += float(line.strip().strip("&").split("&")[1])
    m = re.search(r"#\s*ia_untracked_mass_msun Z=[^=]*=\s*([0-9.eE+-]+)", a.table.read_text(encoding="ascii"))
    untracked_ia = float(m.group(1)) if m else float("nan")
    if not m:
        errors.append("header lacks ia_untracked_mass_msun")
    if any(v < 0 for v in ia[1:]) or abs(sum(ia[1:]) - ia[0]) > 1e-3 * ia[0]:
        errors.append("Ia row violates non-negativity or sum == ejecta")
    if not abs(ia[0] + untracked_ia - ia_total) <= 1e-6 * ia_total:
        errors.append(f"Ia tracked {ia[0]} + untracked {untracked_ia} != W7 total {ia_total}")
    print(f"table sha256 = {sha(a.table)}")
    print(f"agb masses {ch['agb']['mass']}  cc masses {ch['cc']['mass']}  Z {ch['agb']['metallicity']}")
    print(f"max |ejecta - (M - M_final)|, |sum gross - ejecta| relative = {worst_ej:.3e}")
    print(f"max |sum tracked net| / ejecta = {worst_net:.4e} (declared tolerance {tol:g})")
    print(f"max source non-closure |sum gross - ej*sum X0| / ejecta = {worst_closure:.3e}")
    print(f"min tracked gross element ejecta = {min_gross:.3e} Msun")
    print(f"Ia: total W7 ejecta {ia_total:.6f}, tracked {ia[0]:.6f}, untracked (excluded) {untracked_ia:.6f} Msun "
          f"({untracked_ia / ia_total:.2%}); Fe {ia[11]:.4f}")
    for e in errors:
        print("ERROR:", e)
    print("VALID" if not errors else f"INVALID ({len(errors)} errors)")
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
