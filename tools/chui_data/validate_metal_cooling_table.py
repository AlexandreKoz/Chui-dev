#!/usr/bin/env python3
"""Integrity check of a chui_metal_cooling_table_v1 file (D3, 2026-11-04).

Reader-grammar checks (mirrors src/physics/metal_cooling_table.cpp):
  magic; citation/license/solar_reference/uv_background present; known keys
  only; finite numbers; strictly increasing log10_T, log10_nH, redshift;
  0 < reference_hydrogen_mass_fraction < 1; reference_mass_fraction per
  element > 0; every lambda block and ne_over_nH hold N_z*N_T*N_nH values.
Physics checks (with --sources, h5py required):
  * the converted value equals the source Net_Cooling bit-for-bit after the
    9-significant-digit float32 round trip, for every element and z;
  * sum over the nine element increments vs WSS09 Total_Metals/Net_cooling
    at solar abundance (reported, not a pass criterion: Total_Metals is an
    independent all-metals CLOUDY run, and the per-element increments are
    defined as differences; non-additivity is physical);
  * Metal_free (H, He, Compton) is NOT present (no double counting);
  * redshift coverage and the signed (heating) fraction are reported.
"""
from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

KNOWN = {"elements", "log10_T", "log10_nH", "redshift", "reference_mass_fraction",
         "reference_hydrogen_mass_fraction", "ne_over_nH", "lambda"}
GROUP = {"C": "Carbon", "N": "Nitrogen", "O": "Oxygen", "Ne": "Neon", "Mg": "Magnesium",
         "Si": "Silicon", "S": "Sulphur", "Ca": "Calcium", "Fe": "Iron"}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("table", type=Path)
    p.add_argument("--sources", nargs="*", type=Path, default=[])
    a = p.parse_args()
    errors, header, blocks, elements = [], {}, {}, []
    digest = hashlib.sha256()
    magic = False
    with a.table.open("rb") as fh:
        for raw in fh:
            digest.update(raw)
            line = raw.decode("ascii").strip()
            if not line:
                continue
            if line.startswith("#"):
                body = line[1:].strip()
                magic = magic or body == "chui_metal_cooling_table_v1"
                if "=" in body:
                    k, v = body.split("=", 1)
                    header.setdefault(k.strip(), v.strip())
                continue
            tok = line.split()
            key = tok[0]
            if key not in KNOWN:
                errors.append(f"unknown key {key}")
                continue
            if key == "elements":
                elements = tok[1:]
                continue
            name = tok[1] if key == "lambda" else key
            vals = [float(x) for x in (tok[2:] if key == "lambda" else tok[1:])]
            if not all(math.isfinite(v) for v in vals):
                errors.append(f"non-finite in {name}")
            blocks[name] = vals
    if not magic:
        errors.append("missing magic")
    for k in ("citation", "license", "solar_reference", "uv_background"):
        if not header.get(k):
            errors.append(f"missing {k}")
    lt, ln, lz = blocks.get("log10_T", []), blocks.get("log10_nH", []), blocks.get("redshift", [])
    for nm, ax in (("log10_T", lt), ("log10_nH", ln), ("redshift", lz)):
        if not ax or any(not b > c for c, b in zip(ax, ax[1:])):
            errors.append(f"{nm} empty or not strictly increasing")
    xh = blocks.get("reference_hydrogen_mass_fraction", [0])[0]
    if not 0 < xh < 1:
        errors.append("reference_hydrogen_mass_fraction out of (0,1)")
    rmf = blocks.get("reference_mass_fraction", [])
    if len(rmf) != len(elements) or any(v <= 0 for v in rmf):
        errors.append("reference_mass_fraction must be positive, one per element")
    n = len(lt) * len(ln) * len(lz)
    for e in elements + ["ne_over_nH"]:
        name = e if e == "ne_over_nH" else e
        if len(blocks.get(name, [])) != n:
            errors.append(f"block {name} has {len(blocks.get(name, []))} values, expected {n}")
    if any(v <= 0 for v in blocks.get("ne_over_nH", [1])):
        errors.append("ne_over_nH must be positive")
    for forbidden in ("H", "He", "Metal_free"):
        if forbidden in blocks:
            errors.append(f"block {forbidden} present: H/He/Compton belong to the CHUI network")
    neg = sum(sum(1 for v in blocks[e] if v < 0) for e in elements) / max(1, n * len(elements))
    print(f"table sha256 = {digest.hexdigest()}")
    print(f"elements = {elements}; grid N_z={len(lz)} N_T={len(lt)} N_nH={len(ln)}")
    print(f"z = [{lz[0]:g}, {lz[-1]:g}]; T = [{10**lt[0]:.3g}, {10**lt[-1]:.3g}] K; "
          f"n_H = [{10**ln[0]:.3g}, {10**ln[-1]:.3g}] cm^-3")
    print(f"X_H,ref = {xh}; signed (net heating) fraction = {neg:.3f}")
    if a.sources:
        import h5py  # noqa: PLC0415
        import numpy as np  # noqa: PLC0415
        nt, nn = len(lt), len(ln)
        worst_rt, worst_add = 0.0, 0.0
        by_z = {}
        for path in a.sources:
            with h5py.File(path, "r") as f:
                by_z[float(f["Header/Redshift"][()][0])] = path
        for iz, z in enumerate(lz):
            path = min(by_z, key=lambda zz: abs(zz - z))
            with h5py.File(by_z[path], "r") as f:
                total = np.zeros((nt, nn))
                for e in elements:
                    src = f[GROUP[e] + "/Net_Cooling"][()].astype(np.float64)
                    conv = np.array(blocks[e][iz * nt * nn:(iz + 1) * nt * nn]).reshape(nt, nn)
                    if not np.array_equal(conv.astype(np.float32), src.astype(np.float32)):
                        errors.append(f"z={z}: {e} differs from the source after float32 round trip")
                    worst_rt = max(worst_rt, float(np.max(np.abs(conv - src) / np.maximum(np.abs(src), 1e-300))))
                    total += conv
                tm = f["Total_Metals"]
                ref = tm["Net_Cooling" if "Net_Cooling" in tm else "Net_cooling"][()].astype(np.float64)
                mask = np.abs(ref) > 1e-26
                worst_add = max(worst_add, float(np.max(np.abs(total[mask] - ref[mask]) / np.abs(ref[mask]))))
        print(f"max relative round-trip difference vs source = {worst_rt:.2e}")
        print(f"max |sum_i Lambda_i - Total_Metals| / |Total_Metals| (|ref| > 1e-26) = {worst_add:.3f} "
              "(non-additivity of per-element increments; informational)")
    for e in errors[:20]:
        print("ERROR:", e)
    print("VALID" if not errors else f"INVALID ({len(errors)} errors)")
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
