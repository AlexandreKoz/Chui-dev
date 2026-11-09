#!/usr/bin/env python3
"""Data-only integrity check of a ``chui_metal_ion_rates_v1`` file (C02).

Patch 53. This is NOT a CHUI runtime test and NOT scientific qualification.

Reader-grammar checks mirror ``MetalIonRateTable::parse`` in
``src/physics/metal_ion_network.cpp``:
  * keywords: schema, provenance, license, verification, energy_convention,
    ionization_potential, temperature_grid_log10_k, redshift_grid, rate,
    photo; anything else is an error;
  * numbers parse completely and are finite, except the literal ``-inf`` in
    log10 series (an exact zero rate);
  * grids have a count >= 2 and are strictly increasing;
  * every present element covers every stage (C 7, O 9, Si 15, Fe 27) with
    collisional_ionization, recombination and cooling, each with one value
    per grid point; duplicates are refused;
  * the bare nucleus cannot ionize (all -inf) and the neutral atom cannot
    recombine (all -inf);
  * photo rates, if any, cover every ionizable stage of that element;
  * ionization potentials are positive, one per ionizable stage (stage < Z),
    no duplicates; the energy contract (energy_convention + all potentials)
    is reported per element.
Optional ``--gf12-dir`` compares every cooling series with Gnat & Ferland
(2012) ion-by-ion tables in VizieR/CDS byte-by-byte format (files C.txt,
O.txt, Si.txt, Fe.txt): the C02 values must equal log10 of the tabulated
efficiencies to the declared quantization tolerance.

The verdict never upgrades ``verification unverified``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

STAGES = {"C": 7, "O": 9, "Si": 15, "Fe": 27}
TEMPERATURE_PROCESSES = ("collisional_ionization", "recombination", "cooling")
PHOTO_PROCESSES = ("ionization", "heating")
KEYWORDS = {"schema", "provenance", "license", "verification", "energy_convention", "ionization_potential",
            "temperature_grid_log10_k", "redshift_grid", "rate", "photo"}


def parse_log10(token: str, context: str) -> float:
    if token == "-inf":
        return -math.inf
    try:
        value = float(token)
    except ValueError:
        raise SystemExit(f"{context}: invalid number '{token}'")
    if not math.isfinite(value):
        raise SystemExit(f"{context}: invalid number '{token}'")
    return value


def parse(path: Path) -> dict:
    out = {"schema": None, "provenance": "", "license": "", "verification": None, "energy_convention": "",
           "potentials": {}, "grids": {}, "series": {}}
    for number, raw in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        line = raw.split("#", 1)[0]
        tokens = line.split()
        if not tokens:
            continue
        key, rest = tokens[0], tokens[1:]
        context = f"{path.name} line {number}"
        if key not in KEYWORDS:
            raise SystemExit(f"{context}: unknown keyword '{key}'")
        if key == "schema":
            if rest[:1] != ["chui_metal_ion_rates_v1"]:
                raise SystemExit(f"{context}: unsupported schema {rest[:1]}")
            out["schema"] = rest[0]
        elif key in ("provenance", "license"):
            out[key] = " ".join(rest)
        elif key == "verification":
            if rest[:1] not in (["verified_against_source"], ["unverified"]):
                raise SystemExit(f"{context}: verification must be verified_against_source or unverified")
            out["verification"] = rest[0]
        elif key == "energy_convention":
            if rest[:1] != ["collisional_ionization_work_in_cooling"]:
                raise SystemExit(f"{context}: unsupported energy_convention {rest[:1]}")
            if out["energy_convention"]:
                raise SystemExit(f"{context}: energy_convention declared twice")
            out["energy_convention"] = rest[0]
        elif key == "ionization_potential":
            if len(rest) != 3:
                raise SystemExit(f"{context}: malformed ionization_potential")
            symbol, stage, token = rest[0], int(rest[1]), rest[2]
            if symbol not in STAGES:
                raise SystemExit(f"{context}: element '{symbol}' is not one of C, O, Si, Fe")
            if stage + 1 >= STAGES[symbol]:
                raise SystemExit(f"{context}: ionization_potential stage must be < Z")
            ev = float(token)
            if not (math.isfinite(ev) and ev > 0.0):
                raise SystemExit(f"{context}: ionization potential must be a positive finite eV value")
            if (symbol, stage) in out["potentials"]:
                raise SystemExit(f"{context}: duplicate ionization_potential")
            out["potentials"][(symbol, stage)] = ev
        elif key in ("temperature_grid_log10_k", "redshift_grid"):
            count = int(rest[0])
            if count < 2 or len(rest) - 1 != count:
                raise SystemExit(f"{context}: grid count mismatch")
            values = [parse_log10(t, context) for t in rest[1:]]
            if any(not math.isfinite(v) for v in values) or any(not b > a for a, b in zip(values, values[1:])):
                raise SystemExit(f"{context}: grid must be finite and strictly increasing")
            if key in out["grids"]:
                raise SystemExit(f"{context}: grid declared twice")
            out["grids"][key] = values
        else:  # rate / photo
            symbol, stage, process = rest[0], int(rest[1]), rest[2]
            if symbol not in STAGES or stage >= STAGES[symbol]:
                raise SystemExit(f"{context}: bad element/stage")
            allowed = TEMPERATURE_PROCESSES if key == "rate" else PHOTO_PROCESSES
            if process not in allowed:
                raise SystemExit(f"{context}: unknown {key} process '{process}'")
            grid = out["grids"].get("temperature_grid_log10_k" if key == "rate" else "redshift_grid")
            if grid is None:
                raise SystemExit(f"{context}: rate line precedes its grid declaration")
            values = [parse_log10(t, context) for t in rest[3:]]
            if len(values) != len(grid):
                raise SystemExit(f"{context}: {len(values)} values for {len(grid)} grid points")
            series_key = (key, symbol, stage, process)
            if series_key in out["series"]:
                raise SystemExit(f"{context}: duplicate rate line")
            out["series"][series_key] = values
    if out["schema"] is None or not out["provenance"] or not out["license"] or out["verification"] is None:
        raise SystemExit("missing schema / provenance / license / verification")
    if "temperature_grid_log10_k" not in out["grids"]:
        raise SystemExit("missing temperature grid")
    return out


def structural_report(table: dict) -> dict:
    present = sorted({k[1] for k in table["series"]}, key=list(STAGES).index)
    report = {"elements": present, "stages": {}, "energy_contract": {}, "minus_inf_counts": {}}
    for symbol in present:
        stages = STAGES[symbol]
        report["stages"][symbol] = stages
        for stage in range(stages):
            for process in TEMPERATURE_PROCESSES:
                if ("rate", symbol, stage, process) not in table["series"]:
                    raise SystemExit(f"{symbol} stage {stage} is missing {process}")
        bare = table["series"][("rate", symbol, stages - 1, "collisional_ionization")]
        if any(v != -math.inf for v in bare):
            raise SystemExit(f"bare {symbol} nucleus cannot ionize")
        neutral = table["series"][("rate", symbol, 0, "recombination")]
        if any(v != -math.inf for v in neutral):
            raise SystemExit(f"neutral {symbol} cannot recombine")
        photo = any(k[0] == "photo" and k[1] == symbol for k in table["series"])
        if photo:
            for stage in range(stages - 1):
                for process in PHOTO_PROCESSES:
                    if ("photo", symbol, stage, process) not in table["series"]:
                        raise SystemExit(f"photo rates for {symbol} must cover every ionizable stage")
        report["energy_contract"][symbol] = bool(table["energy_convention"]) and all(
            (symbol, s) in table["potentials"] for s in range(stages - 1))
        for process in TEMPERATURE_PROCESSES:
            report["minus_inf_counts"][f"{symbol}.{process}"] = sum(
                sum(1 for v in table["series"][("rate", symbol, s, process)] if v == -math.inf) for s in range(stages))
        cooling_inf = report["minus_inf_counts"][f"{symbol}.cooling"]
        if cooling_inf:
            report.setdefault("warnings", []).append(f"{symbol}: {cooling_inf} cooling entries are -inf (exact zero)")
    report["rate_rows"] = sum(1 for k in table["series"] if k[0] == "rate")
    report["photo_rows"] = sum(1 for k in table["series"] if k[0] == "photo")
    report["ionization_potentials"] = len(table["potentials"])
    grid = table["grids"]["temperature_grid_log10_k"]
    report["temperature_grid"] = {"count": len(grid), "min_K": 10 ** grid[0], "max_K": 10 ** grid[-1]}
    return report


def read_gf12(path: Path, stages: int) -> list[tuple[float, list[float]]]:
    rows = []
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        cells = raw.split()
        if len(cells) != stages + 2:  # T, Z+1 ions, CIE
            continue
        try:
            values = [float(c) for c in cells]
        except ValueError:
            continue
        rows.append((values[0], values[1:1 + stages]))
    return rows


def gf12_crosscheck(table: dict, directory: Path, tolerance_dex: float) -> dict:
    grid = table["grids"]["temperature_grid_log10_k"]
    result = {}
    for symbol, stages in STAGES.items():
        if ("rate", symbol, 0, "cooling") not in table["series"]:
            continue
        source = directory / f"{symbol}.txt"
        if not source.exists():
            result[symbol] = {"status": "source file absent"}
            continue
        rows = read_gf12(source, stages)
        if len(rows) != len(grid):
            result[symbol] = {"status": f"row count {len(rows)} != grid {len(grid)}"}
            continue
        worst_t = 0.0
        worst = 0.0
        compared = 0
        zero_mismatch = 0
        for i, (t, values) in enumerate(rows):
            worst_t = max(worst_t, abs(math.log10(t) - grid[i]))
            for stage, v in enumerate(values):
                c02 = table["series"][("rate", symbol, stage, "cooling")][i]
                if v <= 0.0:
                    zero_mismatch += 0 if c02 == -math.inf else 1
                    continue
                worst = max(worst, abs(math.log10(v) - c02))
                compared += 1
        result[symbol] = {"status": "ok" if worst <= tolerance_dex and zero_mismatch == 0 and worst_t <= 3e-3
                          else "MISMATCH",
                          "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                          "values_compared": compared, "max_abs_log10_difference": worst,
                          "max_abs_log10_temperature_difference": worst_t, "zero_mismatches": zero_mismatch}
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table", type=Path)
    parser.add_argument("--expect-sha256", default="")
    parser.add_argument("--gf12-dir", type=Path, default=None,
                        help="directory with Gnat & Ferland 2012 ion-by-ion tables C.txt, O.txt, Si.txt, Fe.txt")
    parser.add_argument("--gf12-tolerance-dex", type=float, default=0.0026,
                        help="3-significant-digit source rounding (<= 0.0022 dex) + 0.001 dex quantization")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    digest = hashlib.sha256(args.table.read_bytes()).hexdigest()
    if args.expect_sha256 and digest != args.expect_sha256:
        print(f"SHA-256 MISMATCH: {digest} != {args.expect_sha256}", file=sys.stderr)
        return 2
    table = parse(args.table)
    report = {"file": args.table.name, "sha256": digest, "bytes": args.table.stat().st_size,
              "schema": table["schema"], "verification": table["verification"],
              "energy_convention": table["energy_convention"], "license": table["license"],
              "provenance": table["provenance"]}
    report.update(structural_report(table))
    if args.gf12_dir is not None:
        report["gf12_cooling_crosscheck"] = gf12_crosscheck(table, args.gf12_dir, args.gf12_tolerance_dex)
    report["verdict"] = "STRUCTURALLY_VALID (data-only; not scientific qualification)"
    if table["verification"] != "verified_against_source":
        report["verdict"] += "; dataset declares verification = unverified (research opt-in required)"
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        for key, value in report.items():
            print(f"{key} = {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
