#!/usr/bin/env python3
"""Validate a TREECOOL-format UV background file against CHUI's reader
contract and print a provenance manifest (SHA-256, redshift range, rows).

STATUS: written in CHUI Part 2.5; NOT RUN by the author (executable
validation is USER-RESERVED).

Reader contract (src/physics/radiation_background.cpp,
UvBackgroundTable::loadTreecool), verified against the Faucher-Giguere (2020,
MNRAS 493, 1614) TREECOOL README:
  column 1  log10(1 + z)                 strictly increasing
  column 2  Gamma_HI    [s^-1]           photoionization rate per HI
  column 3  Gamma_HeI   [s^-1]
  column 4  Gamma_HeII  [s^-1]
  column 5  qdot_HI     [erg s^-1]       photoheating rate per HI
  column 6  qdot_HeI    [erg s^-1]
  column 7  qdot_HeII   [erg s^-1]
  '#' lines and blank lines are ignored; all values finite and >= 0.
CHUI interpolates log10 of the rates in log10(1+z) and returns ZERO rates
above the last row (before reionization), so the file's z_max is physically
significant. Tiny placeholder values (FG20 writes 1e-40 where a rate is
undefined) are legal; this tool reports them so the owner can see where the
background is effectively off.

The tool never edits the file.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

COLUMNS = ["log10(1+z)", "Gamma_HI", "Gamma_HeI", "Gamma_HeII", "qdot_HI", "qdot_HeI", "qdot_HeII"]
PLACEHOLDER = 1.0e-39


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("treecool", type=Path)
    parser.add_argument("--expect-sha256", default="")
    args = parser.parse_args()

    digest = hashlib.sha256()
    rows = []
    errors = []
    with args.treecool.open("rb") as handle:
        for number, raw in enumerate(handle, start=1):
            digest.update(raw)
            line = raw.decode("ascii", "replace").strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 7:
                errors.append(f"row {number}: {len(fields)} columns, expected 7")
                continue
            try:
                values = [float(x) for x in fields[:7]]
            except ValueError:
                errors.append(f"row {number}: non-numeric value")
                continue
            if any(not math.isfinite(v) or v < 0.0 for v in values):
                errors.append(f"row {number}: non-finite or negative value")
            if len(fields) > 7:
                errors.append(f"row {number}: {len(fields)} columns; CHUI reads only the first 7 (check the format)")
            rows.append((number, values))
    sha = digest.hexdigest()
    if len(rows) < 2:
        errors.append("fewer than 2 data rows")
    for (n0, a), (n1, b) in zip(rows, rows[1:]):
        if not b[0] > a[0]:
            errors.append(f"rows {n0}->{n1}: log10(1+z) not strictly increasing")
    if args.expect_sha256 and args.expect_sha256 != sha:
        errors.append(f"SHA-256 {sha} != expected {args.expect_sha256}")

    print(f"file            = {args.treecool}")
    print(f"sha256          = {sha}")
    print(f"rows            = {len(rows)}")
    if rows:
        z_lo = 10.0 ** rows[0][1][0] - 1.0
        z_hi = 10.0 ** rows[-1][1][0] - 1.0
        print(f"redshift_range  = [{z_lo:.4g}, {z_hi:.4g}]  (CHUI: zero UVB above z = {z_hi:.4g})")
        for c in range(1, 7):
            placeholders = [n for n, v in rows if 0.0 < v[c] <= PLACEHOLDER]
            zeros = [n for n, v in rows if v[c] == 0.0]
            first = rows[0][1][c]
            print(f"{COLUMNS[c]:<10} z_min value = {first:.4e}; placeholder(<=1e-39) rows = {len(placeholders)};"
                  f" zero rows = {len(zeros)}")
        # Physical sanity (reported, not enforced): mean photoheating energy
        # per ionization qdot/Gamma should be a few eV for HI.
        number, v = rows[0]
        if v[1] > 0.0:
            print(f"qdot_HI/Gamma_HI at z_min = {v[4] / v[1] / 1.602176634e-12:.3f} eV (expected O(1-10) eV)")
    print("config:\n  uv_background_model = treecool_table\n"
          f"  uv_background_table_path = {args.treecool}\n  uv_background_table_sha256 = {sha}")
    if errors:
        print("INVALID for CHUI:")
        for e in errors:
            print("  - " + e)
        return 1
    print("valid for the CHUI TREECOOL reader contract")
    return 0


if __name__ == "__main__":
    sys.exit(main())
