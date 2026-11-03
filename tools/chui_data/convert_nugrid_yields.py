#!/usr/bin/env python3
"""Convert NuPyCEE/NuGrid-format stellar yield tables into CHUI's
``chui_ssp_yield_table_v1`` (src/physics/stellar_population.cpp).

STATUS: written in CHUI Part 2.5; NOT RUN by the author (executable
validation is USER-RESERVED). The input layout was established by reading
NuPyCEE ``yield_tables/agb_and_massive_stars_nugrid_MESAonly_fryer12delay.txt``
and ``yield_tables/sn1a_i99_W7.txt`` (see tools/chui_data/data_sources.md).

Input contract (NuPyCEE ASCII):
  'H Table: (M=<Msun>,Z=<Z>)' starts a model; 'H Mfinal: <Msun>' gives the
  final (remnant) mass; '&Isotopes &Yields &X0 &Z &A' header, then rows
  '&<El>-<A> &<yield Msun> &<X0> ...'. 'Yields' are TOTAL (gross) ejected
  masses per isotope; X0 is the initial mass fraction of the star.
  SN Ia file: '&Isotopes &Z=<Z> ...' header, rows '&<El>-<A> &<Msun> ...'.

Conversion (no rebalancing, ever):
  ejecta  m_ej = M_initial - M_final                      (per star, Msun)
  net     p_e  = sum_{isotopes of e} (yield_iso - X0_iso * m_ej)
  for the 11 CHUI elements H He C N O Ne Mg Si S Ca Fe. Net yields of all
  OTHER elements (Na, Al, Ar, Ni, s-process, ...) are UNTRACKED: they are
  reported per model and their sum makes the tracked net sum non-zero. The
  CHUI loader refuses a table whose |sum p_e| / m_ej exceeds the declared
  net_yield_sum_tolerance; this tool writes the tolerance the user passes and
  REFUSES to write if any model exceeds it (it never moves the residual into
  hydrogen or any other element). Reading the Fryer (2012) delayed-remnant
  set showed untracked net fractions up to ~1.7% of the ejecta for 12-20 Msun
  models (mostly Ar, Ni, Cr, Ti), so a tolerance of ~0.02 is needed and must
  be chosen explicitly.

  Channel assignment: models with M < --cc-mass-min are 'agb', others 'cc'.
  The CHUI grid must be rectangular: every metallicity must have the same
  mass list per channel, otherwise the tool refuses (no regridding).

  SN Ia: per-event ELEMENT masses. CHUI requires them to sum to the ejecta
  within 1e-3, so untracked Ia elements (W7: ~0.15 Msun of Ni, Mn, Cr, Ar,
  Ti, ...) need an explicit policy: --ia-untracked exclude defines the Ia
  ejecta as the tracked element mass (the untracked mass is NOT returned to
  the gas and is reported in the header); there is no default. If the Ia
  file has one metallicity it is used for all metallicities (the CHUI Ia
  axis then has one entry).

Licence: NuPyCEE is distributed under a BSD 3-clause licence (repository
LICENSE); the yields themselves come from the cited publications (Pignatari
et al. 2016; Ritter et al. 2018; Iwamoto et al. 1999). --license and
--citation are REQUIRED and copied into the header.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

CHUI_ELEMENTS = ["H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"]
CONVERTER_VERSION = "chui_nugrid_yield_converter_v1"


def fail(message: str) -> None:
    raise SystemExit(f"convert_nugrid_yields: {message}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def element_of(isotope: str) -> str:
    match = re.fullmatch(r"([A-Z][a-z]?)-?(\d+)", isotope.strip())
    if not match:
        fail(f"unrecognized isotope label '{isotope}'")
    return match.group(1)


def parse_massive(path: Path):
    """Yield (mass, Z, m_final, [(isotope, yield, x0)]) one model at a time."""
    model = None
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("H Table:"):
                if model is not None:
                    yield model
                match = re.search(r"M\s*=\s*([0-9.eE+-]+)\s*,\s*Z\s*=\s*([0-9.eE+-]+)", line)
                if not match:
                    fail(f"{path}: cannot parse '{line.strip()}'")
                model = [float(match.group(1)), float(match.group(2)), None, []]
            elif line.startswith("H Mfinal:") and model is not None:
                model[2] = float(line.split(":", 1)[1])
            elif line.startswith("&") and model is not None and not line.startswith("&Isotopes"):
                fields = [f.strip() for f in line.strip().strip("&").split("&")]
                if len(fields) < 3:
                    fail(f"{path}: malformed row '{line.strip()}'")
                model[3].append((fields[0], float(fields[1]), float(fields[2])))
    if model is not None:
        yield model


def parse_ia(path: Path):
    metallicities = []
    masses = defaultdict(lambda: defaultdict(float))  # z index -> element -> Msun
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("&Isotopes"):
                metallicities = [float(f.split("=", 1)[1]) for f in line.strip().strip("&").split("&")[1:]
                                 if "=" in f]
                if not metallicities:
                    fail(f"{path}: Ia header lists no metallicity")
            elif line.startswith("&"):
                fields = [f.strip() for f in line.strip().strip("&").split("&")]
                if len(fields) < 1 + len(metallicities):
                    fail(f"{path}: malformed Ia row '{line.strip()}'")
                for iz in range(len(metallicities)):
                    value = float(fields[1 + iz])
                    if value < 0.0:
                        fail(f"{path}: negative Ia mass")
                    masses[iz][element_of(fields[0])] += value
    if not metallicities:
        fail(f"{path}: no Ia table found")
    return metallicities, masses


def fmt(v: float) -> str:
    return f"{v:.8e}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--massive", required=True, type=Path, help="AGB + massive-star NuPyCEE table")
    parser.add_argument("--ia", required=True, type=Path, help="SN Ia NuPyCEE table")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--citation", required=True)
    parser.add_argument("--license", required=True)
    parser.add_argument("--net-yield-sum-tolerance", required=True, type=float,
                        help="relative |sum of tracked net yields|/ejecta allowed per model (CHUI max 0.05)")
    parser.add_argument("--cc-mass-min", type=float, default=8.0)
    parser.add_argument("--ia-untracked", required=True, choices=["exclude"],
                        help="'exclude': Ia ejecta := tracked element mass (untracked mass reported, not returned)")
    args = parser.parse_args()
    if not 0.0 <= args.net_yield_sum_tolerance <= 0.05:
        fail("--net-yield-sum-tolerance must be in [0, 0.05] (CHUI loader limit)")

    grid = {"agb": defaultdict(dict), "cc": defaultdict(dict)}  # channel -> Z -> M -> row
    worst = 0.0
    report = []
    for mass, z, m_final, isotopes in parse_massive(args.massive):
        if m_final is None:
            fail(f"model M={mass} Z={z} has no 'H Mfinal'")
        ejecta = mass - m_final
        if not ejecta >= 0.0:
            fail(f"model M={mass} Z={z}: negative ejecta {ejecta}")
        net = dict.fromkeys(CHUI_ELEMENTS, 0.0)
        untracked = defaultdict(float)
        for isotope, gross, x0 in isotopes:
            el = element_of(isotope)
            p = gross - x0 * ejecta
            if el in net:
                net[el] += p
            else:
                untracked[el] += p
        tracked_sum = sum(net.values())
        residual = abs(tracked_sum) / max(ejecta, 1e-30)
        worst = max(worst, residual)
        top = sorted(untracked.items(), key=lambda kv: -abs(kv[1]))[:3]
        report.append(f"# model M={mass:g} Z={z:g} ejecta={ejecta:.6g} tracked_net_sum/ejecta={tracked_sum / max(ejecta, 1e-30):+.3e}"
                      f" largest_untracked=" + ",".join(f"{k}:{v:+.2e}" for k, v in top))
        if residual > args.net_yield_sum_tolerance:
            fail(f"model M={mass} Z={z}: |tracked net sum|/ejecta = {residual:.3e} exceeds "
                 f"--net-yield-sum-tolerance {args.net_yield_sum_tolerance}; refusing (yields are never rebalanced)")
        channel = "agb" if mass < args.cc_mass_min else "cc"
        if mass in grid[channel][z]:
            fail(f"duplicate model M={mass} Z={z}")
        grid[channel][z][mass] = (ejecta, [net[e] for e in CHUI_ELEMENTS])

    axes = {}
    for channel, by_z in grid.items():
        if not by_z:
            fail(f"no '{channel}' models (check --cc-mass-min)")
        zs = sorted(by_z)
        masses = sorted(by_z[zs[0]])
        for z in zs:
            if sorted(by_z[z]) != masses:
                fail(f"channel {channel}: mass grid at Z={z} differs from Z={zs[0]} (no regridding)")
        axes[channel] = (zs, masses)

    ia_z, ia_masses = parse_ia(args.ia)
    order = sorted(range(len(ia_z)), key=lambda i: ia_z[i])
    ia_rows = []
    ia_untracked = []
    for iz in order:
        tracked = [ia_masses[iz].get(e, 0.0) for e in CHUI_ELEMENTS]
        untracked = sum(v for el, v in ia_masses[iz].items() if el not in CHUI_ELEMENTS)
        ia_rows.append((ia_z[iz], sum(tracked), tracked))
        ia_untracked.append(untracked)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii") as out:
        out.write("# chui_ssp_yield_table_v1\n")
        out.write(f"# citation = {args.citation}\n")
        out.write(f"# license = {args.license}\n")
        out.write(f"# net_yield_sum_tolerance = {args.net_yield_sum_tolerance:g}\n")
        out.write(f"# converter = {CONVERTER_VERSION}\n")
        out.write(f"# source_massive = {args.massive.name} sha256={sha256_file(args.massive)}\n")
        out.write(f"# source_ia = {args.ia.name} sha256={sha256_file(args.ia)}\n")
        out.write("# convention = ejecta = M - M_final; net_e = sum_iso(yield - X0 * ejecta); no rebalancing\n")
        out.write(f"# max_tracked_net_residual = {worst:.4e}\n")
        out.write("# ia_untracked_policy = exclude (Ia ejecta = tracked element mass; untracked mass not returned)\n")
        for (z, _, _), u in zip(ia_rows, ia_untracked):
            out.write(f"# ia_untracked_mass_msun Z={z:g} = {u:.6e}\n")
        for line in report:
            out.write(line + "\n")
        out.write("elements " + " ".join(CHUI_ELEMENTS) + "\n")
        for channel in ("agb", "cc"):
            zs, masses = axes[channel]
            out.write(f"channel {channel}\n")
            out.write("metallicity " + " ".join(f"{z:g}" for z in zs) + "\n")
            out.write("mass " + " ".join(f"{m:g}" for m in masses) + "\n")
            for iz, z in enumerate(zs):
                for im, m in enumerate(masses):
                    ejecta, net = grid[channel][z][m]
                    out.write(f"row {iz} {im} {fmt(ejecta)} " + " ".join(fmt(v) for v in net) + "\n")
        out.write("channel ia\n")
        out.write("metallicity " + " ".join(f"{z:g}" for z, _, _ in ia_rows) + "\n")
        for iz, (_, ejecta, tracked) in enumerate(ia_rows):
            out.write(f"row {iz} {fmt(ejecta)} " + " ".join(fmt(v) for v in tracked) + "\n")
    print(f"wrote {args.output}; max tracked net residual {worst:.3e}; sha256 {sha256_file(args.output)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
