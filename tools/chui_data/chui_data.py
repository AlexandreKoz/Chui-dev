#!/usr/bin/env python3
"""CHUI external physics data: one data-only command interface (patch 53).

Subcommands (none of them builds or runs CHUI):

  inspect  <dir>                 list every registered product found in <dir>
  verify   <dir>                 SHA-256 check of every registered file present
                                 (exit 2 on a mismatch, 3 if --require-all and one is missing)
  stage    --to <dir> <src>...   copy registered files from bundle directories
                                 or zip archives into <dir>, after verifying the
                                 recorded SHA-256 (never overwrites a differing file)
  fetch    <dir>                 run tools/chui_data/fetch_and_convert.sh <dir>
                                 (pinned upstream hashes; needs network, curl,
                                 git, and h5py+numpy for the WSS09 conversion)
  emit-config --profile {safe_default,research_c02} [--data-dir chui_data]
                                 print a .param fragment using the registry hashes
  sha256sums <dir>               write <dir>/SHA256SUMS for `sha256sum -c`
  report   <dir> [--json]        machine-readable fact table (registry + local bytes)

Policy: data files live in a user-local directory (default ``chui_data/``,
git-ignored). Nothing here downloads executable code or runs anything from a
data archive; zip members are read as bytes only. Relative paths are resolved
against the current directory, so run CHUI from the directory that holds the
data directory named in the configuration.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "dataset_registry.json"


def load_registry() -> dict:
    return json.loads(REGISTRY.read_text(encoding="utf-8"))


def files(registry: dict) -> list[dict]:
    return [d for d in registry["datasets"] if d.get("file")]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cmd_inspect(args) -> int:
    root = Path(args.dir)
    for d in files(load_registry()):
        path = root / d["file"]
        if path.exists():
            actual = sha256_path(path)
            state = "OK" if actual == d["sha256"] else f"SHA-256 MISMATCH ({actual})"
            print(f"{d['subsystem']:<10} {d['file']}: present, {path.stat().st_size} bytes, {state}")
        else:
            print(f"{d['subsystem']:<10} {d['file']}: absent")
    return 0


def cmd_verify(args) -> int:
    root = Path(args.dir)
    status = 0
    for d in files(load_registry()):
        path = root / d["file"]
        if not path.exists():
            print(f"MISSING  {d['file']}")
            if args.require_all:
                status = max(status, 3)
            continue
        actual = sha256_path(path)
        if actual != d["sha256"]:
            print(f"MISMATCH {d['file']}: {actual} != {d['sha256']}")
            status = 2
        else:
            print(f"OK       {d['file']}")
    return status


def iter_candidates(source: Path):
    if source.is_dir():
        for path in source.rglob("*"):
            if path.is_file():
                yield path.name, path.read_bytes
    elif zipfile.is_zipfile(source):
        with zipfile.ZipFile(source) as archive:
            for info in archive.infolist():
                if info.is_dir():
                    continue
                name = Path(info.filename).name
                yield name, (lambda info=info: zipfile.ZipFile(source).read(info))
    elif source.is_file():
        yield source.name, source.read_bytes


def cmd_stage(args) -> int:
    target = Path(args.to)
    target.mkdir(parents=True, exist_ok=True)
    wanted = {d["file"]: d for d in files(load_registry())}
    staged, status = set(), 0
    for source in map(Path, args.sources):
        for name, read in iter_candidates(source):
            if name not in wanted or name in staged:
                continue
            data = read()
            if sha256_bytes(data) != wanted[name]["sha256"]:
                print(f"REFUSED  {name} from {source}: SHA-256 differs from the registry")
                status = 2
                continue
            destination = target / name
            if destination.exists() and sha256_path(destination) != wanted[name]["sha256"]:
                print(f"REFUSED  {destination} exists with different bytes (not overwritten)")
                status = 2
                continue
            if not destination.exists():
                destination.write_bytes(data)
            staged.add(name)
            print(f"STAGED   {name} -> {destination}")
    for name in sorted(set(wanted) - staged):
        print(f"ABSENT   {name} (not found in the given sources)")
    return status


def cmd_fetch(args) -> int:
    script = HERE / "fetch_and_convert.sh"
    return subprocess.call(["bash", str(script), args.dir])


SAFE_DEFAULT_HEADER = """# SAFE_DEFAULT_AVAILABLE_DATA -- PARTIAL, DIAGNOSTIC fragment (patch 53).
# Generated by tools/chui_data/chui_data.py emit-config --profile safe_default.
# Merge into the [physics] section of a mode deck. It enables ONLY data whose
# consumers run WITHOUT the unqualified-rates research opt-in:
#   S02 NuGrid+W7 yields and F02 SB99 winds with the analytic_reference
#   thermochemistry backend (the only Part 2 backend exempt from the
#   unqualified-rate gate).
# Deliberately EXCLUDED (cannot be assembled consistently without opt-ins):
#   * FG20 UVB: every backend that consumes it (equilibrium_primordial,
#     fixed_hybrid, fixed_detailed) is gated by
#     thermochemistry_allow_unqualified_rates; analytic_reference ignores it.
#   * C01 WSS09: embeds HM01; no HM01 TREECOOL exists, so FG20 + WSS09 is a
#     declared mismatch (metal_cooling_accept_uv_background_mismatch).
#   * C02: declared unverified (research opt-in only).
# Scientific caveats that remain even here: SB99 winds assume a Kroupa IMF
# (CHUI SSP: Chabrier) and a single metallicity Z = 0.014; NuGrid yields cover
            f"ssp_wind_table_sha256 = {f02['sha256']}",
            "ssp_yield_table_coverage = edge_clamp_research",
            "ssp_wind_table_compatibility = approximate_research",
            "fb_supernova_model = mechanical_resolved_unresolved",# Requires a mode deck with star formation and feedback enabled.
[physics]
"""

RESEARCH_C02_HEADER = """# RESEARCH_C02_EXPLICIT_OPT_IN -- RESEARCH-ONLY fragment (patch 53).
# Generated by tools/chui_data/chui_data.py emit-config --profile research_c02.
# The C02 dataset declares `verification unverified`; this fragment therefore
# REQUIRES thermochemistry_allow_unqualified_rates = true. It is NOT
# scientifically qualified and NOT production-ready.
# Physical domain of the dataset (README_C02 / registry):
#   * Gnat & Ferland 2012 ion-by-ion efficiencies, 1e4-1e8 K, low density,
#     optically thin; collisional ionization (Voronov 1997) and RR/DR fits;
#   * NO metal photoionization/photoheating -> uv_background_model = none
#     (collisional regime). Using the FG20 UVB here would photoionize H/He but
#     not C/O/Si/Fe (incoherent radiation context);
#   * metal_cooling_table_range = reject: a cell below 1e4 K is refused (stage
#     rejection) instead of extrapolating the ladder; switch to
#     clamp_density_redshift only as an acknowledged, counted clamp;
#   * no other element cooling is configured (metal_cooling_mode = none):
#     N, Ne, Mg, S, Ca do not cool in this fragment;
#   * do NOT enable chemistry_adaptive_selection without an element table
#     (refused by the configuration validator).
# Energy convention: collisional_ionization_work_in_cooling, consistent with
# Gnat & Ferland 2012 sect. 2; the stored ion potential is a diagnostic
# partition only. CI/RR/DR coefficients are NOT independently source-verified.
[physics]
"""


def cmd_emit_config(args) -> int:
    registry = {d["id"]: d for d in load_registry()["datasets"]}
    data_dir = args.data_dir.rstrip("/")
    lines = []
    if args.profile == "safe_default":
        lines.append(SAFE_DEFAULT_HEADER)
        s02 = registry["S02_nugrid_fryer12delay_w7"]
        f02 = registry["F02_sb99_z014_genevav00"]
        lines += [
            "thermochemistry_backend = analytic_reference",
            "uv_background_model = none",
            "metal_cooling_mode = none",
            "thermochemistry_allow_unqualified_rates = false",
            "stellar_population_model = chabrier_yield_table_v1",
            f"ssp_yield_table_path = {data_dir}/{s02['file']}",
            f"ssp_yield_table_sha256 = {s02['sha256']}",
            "ssp_wind_model = table_v1",
            f"ssp_wind_table_path = {data_dir}/{f02['file']}",
            f"ssp_wind_table_sha256 = {f02['sha256']}",
            "ssp_yield_table_coverage = edge_clamp_research",
            "ssp_wind_table_compatibility = approximate_research",
            "fb_supernova_model = mechanical_resolved_unresolved",
            "ssp_max_event_publication_lag_years = 0.0",
        ]
    else:
        lines.append(RESEARCH_C02_HEADER)
        c02 = registry["C02_metal_ion_rates_C_O_Si_Fe_unverified"]
        lines += [
            "thermochemistry_backend = fixed_detailed",
            "thermochemistry_allow_unqualified_rates = true",
            "chemistry_metal_ion_elements = C,O,Si,Fe",
            f"metal_ion_rate_table_path = {data_dir}/{c02['file']}",
            f"metal_ion_rate_table_sha256 = {c02['sha256']}",
            "metal_cooling_mode = none",
            "metal_cooling_table_range = reject",
            "uv_background_model = none",
            "chemistry_adaptive_selection = off",
        ]
    print("\n".join(lines))
    return 0


def cmd_sha256sums(args) -> int:
    root = Path(args.dir)
    out = []
    for d in files(load_registry()):
        if (root / d["file"]).exists():
            out.append(f"{d['sha256']}  {d['file']}")
    (root / "SHA256SUMS").write_text("\n".join(out) + "\n", encoding="ascii")
    print(f"wrote {root / 'SHA256SUMS'} ({len(out)} entries); verify with: (cd {root} && sha256sum -c SHA256SUMS)")
    return 0


def cmd_report(args) -> int:
    registry = load_registry()
    root = Path(args.dir) if args.dir else None
    rows = []
    for d in registry["datasets"]:
        facts = dict(d.get("facts", {}))
        if root is not None and d.get("file"):
            path = root / d["file"]
            facts["bytes_present"] = path.exists() and sha256_path(path) == d["sha256"]
        rows.append({"id": d["id"], "subsystem": d["subsystem"], "file": d.get("file"), "sha256": d.get("sha256"),
                     "facts": facts, "blocked_reason": d.get("blocked_reason", ""), "labels": d.get("labels", [])})
    if args.json:
        print(json.dumps({"assessed": registry["assessed"], "datasets": rows}, indent=2))
    else:
        keys = list(rows[0]["facts"])
        print("id".ljust(42) + " ".join(k[:12].ljust(12) for k in keys))
        for row in rows:
            print(row["id"].ljust(42) + " ".join(str(row["facts"].get(k, "")).ljust(12) for k in keys))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("inspect"); p.add_argument("dir"); p.set_defaults(func=cmd_inspect)
    p = sub.add_parser("verify"); p.add_argument("dir"); p.add_argument("--require-all", action="store_true")
    p.set_defaults(func=cmd_verify)
    p = sub.add_parser("stage"); p.add_argument("--to", required=True); p.add_argument("sources", nargs="+")
    p.set_defaults(func=cmd_stage)
    p = sub.add_parser("fetch"); p.add_argument("dir"); p.set_defaults(func=cmd_fetch)
    p = sub.add_parser("emit-config"); p.add_argument("--profile", choices=["safe_default", "research_c02"],
                                                      required=True)
    p.add_argument("--data-dir", default="chui_data"); p.set_defaults(func=cmd_emit_config)
    p = sub.add_parser("sha256sums"); p.add_argument("dir"); p.set_defaults(func=cmd_sha256sums)
    p = sub.add_parser("report"); p.add_argument("dir", nargs="?"); p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_report)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
