#!/usr/bin/env python3
"""Convert Wiersma, Schaye & Smith (2009) element-resolved cooling tables
(HDF5, one file per redshift) into CHUI's ``chui_metal_cooling_table_v1``.

STATUS: written in CHUI Part 2.5; NOT RUN by the author (all executable
validation is USER-RESERVED). The source layout was established by reading
``z_0.000.hdf5`` (Header/Version 100422) with h5py, see
``tools/chui_data/data_sources.md``.

Source contract (WSS09 sect. 2, eq. 4 and the table README):
  * ``<Element>/Net_Cooling``   shape (N_T, N_nH): Lambda_i,sun / n_H^2
    [erg cm^3 s^-1] at the CLOUDY 07.02 solar abundances, positive = net
    cooling, negative = net heating, computed as the difference between the
    all-element run and the run without element i.
  * ``<Element>/Temperature_bins`` [K], ``<Element>/Hydrogen_density_bins``
    [cm^-3] (LINEAR values; converted here to log10).
  * ``Solar/Electron_density_over_n_h`` (N_T, N_nH): (n_e/n_H)_sun, the
    denominator of the eq. 4 electron correction.
  * ``Header/Abundances/Solar_mass_fractions`` and ``Abund_names``.
  * ``Metal_free/*`` holds the H+He contribution (and Compton cooling): it is
    deliberately NOT converted, because CHUI's non-equilibrium network owns
    H/He and Compton; writing it would double count.
  * ``Header/Redshift``; tables assume photoionization by the Haardt & Madau
    (2001) background + CMB, optically thin gas, ionization equilibrium.

Output contract (src/physics/metal_cooling_table.cpp):
  index order z slowest, then T, then n_H; one block per line; sign and
  values unchanged (float32 values written with 9 significant digits, which
  round-trips float32 exactly). ``reference_mass_fraction`` and
  ``reference_hydrogen_mass_fraction`` come from Solar_mass_fractions so the
  runtime applies 10^[i/H] = (X_i/X_H)/(X_i,sun/X_H,sun) exactly as eq. 4.

Memory: bounded by one (N_T, N_nH) dataset (~114 kB) at a time; the output is
streamed. The output file is large (~4.6 MB of text per element per
redshift); use --zmax to restrict the redshift range.

Licence: the WSS09 tables carry no licence statement on their distribution
page. The converter therefore REQUIRES --license; it never invents one.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

CONVERTER_VERSION = "chui_wss09_converter_v1"
ELEMENTS = [  # (WSS09 group, CHUI symbol)
    ("Carbon", "C"), ("Nitrogen", "N"), ("Oxygen", "O"), ("Neon", "Ne"), ("Magnesium", "Mg"),
    ("Silicon", "Si"), ("Sulphur", "S"), ("Calcium", "Ca"), ("Iron", "Fe"),
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(f"convert_wss09_metal_cooling: {message}")


def decode(value) -> str:
    if hasattr(value, "tobytes"):
        value = value.tobytes()
    if isinstance(value, bytes):
        return value.decode("ascii", "replace").strip("\x00 ").strip()
    return str(value).strip()


def read_scalar_string(f, path: str) -> str:
    data = f[path][()]
    return decode(data[0] if getattr(data, "shape", ()) else data)


def same_axis(a, b, rel: float = 1.0e-6) -> bool:
    return len(a) == len(b) and all(abs(x - y) <= rel * max(abs(x), abs(y)) for x, y in zip(a, b))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", type=Path, help="WSS09 z_*.hdf5 files (any order)")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--license", required=True,
                        help="licence/redistribution terms under which YOU obtained the tables (no default)")
    parser.add_argument("--uv-background", required=True,
                        help="identity of the UVB embedded in the tables, e.g. "
                             "'HM01+CMB (Haardt & Madau 2001), CLOUDY 07.02'; 'none' only for collisional tables")
    parser.add_argument("--uv-background-treecool-sha256", default="",
                        help="SHA-256 of the TREECOOL file that reproduces the embedded UVB, if one exists")
    parser.add_argument("--zmax", type=float, default=math.inf)
    parser.add_argument("--elements", default=",".join(s for _, s in ELEMENTS),
                        help="comma-separated CHUI symbols to convert (default: all nine)")
    args = parser.parse_args()

    try:
        import h5py  # noqa: PLC0415  (optional dependency, only for this tool)
    except ImportError:
        fail("h5py is required (pip install h5py)")

    wanted = [s.strip() for s in args.elements.split(",") if s.strip()]
    groups = [(g, s) for g, s in ELEMENTS if s in wanted]
    if len(groups) != len(wanted):
        fail(f"unknown element in --elements {wanted}; allowed: {[s for _, s in ELEMENTS]} (H/He are refused)")

    # Pass 1: headers, axes and redshift ordering (no table data retained).
    files = []
    t_axis = n_axis = None
    solar = None
    for path in args.inputs:
        with h5py.File(path, "r") as f:
            z = float(f["Header/Redshift"][()][0])
            note = read_scalar_string(f, "Header/Note") if "Header/Note" in f else ""
            version = read_scalar_string(f, "Header/Version") if "Header/Version" in f else "unknown"
            names = [decode(x) for x in f["Header/Abundances/Abund_names"][()]]
            fractions = [float(x) for x in f["Header/Abundances/Solar_mass_fractions"][()]]
            this_solar = dict(zip(names, fractions))
            for group, _ in groups:
                if group not in f or "Net_Cooling" not in f[group]:
                    fail(f"{path}: missing {group}/Net_Cooling")
                tg = [float(x) for x in f[f"{group}/Temperature_bins"][()]]
                ng = [float(x) for x in f[f"{group}/Hydrogen_density_bins"][()]]
                if t_axis is None:
                    t_axis, n_axis = tg, ng
                if not same_axis(tg, t_axis) or not same_axis(ng, n_axis):
                    fail(f"{path}: {group} axes differ from the first file/group (refusing to regrid)")
                if f[f"{group}/Net_Cooling"].shape != (len(t_axis), len(n_axis)):
                    fail(f"{path}: {group}/Net_Cooling shape {f[group + '/Net_Cooling'].shape} is not (N_T, N_nH)")
            if f["Solar/Electron_density_over_n_h"].shape != (len(t_axis), len(n_axis)):
                fail(f"{path}: Solar/Electron_density_over_n_h has an unexpected shape")
            if solar is None:
                solar = this_solar
            elif any(abs(solar[k] - this_solar.get(k, -1.0)) > 1e-7 for k in solar):
                fail(f"{path}: solar abundances differ between files")
            if "Haardt" in note and args.uv_background.strip().lower() == "none":
                fail(f"{path}: Header/Note says '{note}' but --uv-background none was given")
            if z <= args.zmax:
                files.append((z, path, version, note))
    if not files:
        fail("no input file at or below --zmax")
    files.sort(key=lambda item: item[0])
    for a, b in zip(files, files[1:]):
        if not b[0] > a[0]:
            fail(f"duplicate redshift {a[0]} ({a[1]}, {b[1]})")
    if any(t <= 0.0 for t in t_axis) or any(n <= 0.0 for n in n_axis):
        fail("non-positive axis value")
    log_t = [math.log10(t) for t in t_axis]
    log_n = [math.log10(n) for n in n_axis]
    for axis, name in ((log_t, "log10_T"), (log_n, "log10_nH")):
        if any(not b > a for a, b in zip(axis, axis[1:])):
            fail(f"{name} is not strictly increasing")
    if "Hydrogen" not in solar:
        fail("Solar_mass_fractions lacks Hydrogen")
    for group, _ in groups:
        if not solar.get(group, 0.0) > 0.0:
            fail(f"Solar_mass_fractions lacks a positive {group} entry")

    def fmt(values) -> str:
        out = []
        for v in values:
            v = float(v)
            if not math.isfinite(v):
                fail("non-finite value in a source table")
            out.append(f"{v:.9g}")
        return " ".join(out)

    def stream_block(handle, prefix: str, dataset: str) -> None:
        handle.write(prefix)
        for _, path, _, _ in files:
            with h5py.File(path, "r") as f:
                data = f[dataset][()]  # one (N_T, N_nH) array in memory
                for row in data:      # T slowest, n_H fastest
                    handle.write(" " + fmt(row))
        handle.write("\n")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii") as out:
        out.write("# chui_metal_cooling_table_v1\n")
        out.write(f"# citation = Wiersma, Schaye & Smith 2009, MNRAS 393, 99 (table version {files[0][2]})\n")
        out.write(f"# license = {args.license}\n")
        out.write("# solar_reference = CLOUDY 07.02 default solar abundances (WSS09 Table 1), "
                  "Header/Abundances/Solar_mass_fractions\n")
        out.write(f"# uv_background = {args.uv_background}\n")
        if args.uv_background_treecool_sha256:
            out.write(f"# uv_background_treecool_sha256 = {args.uv_background_treecool_sha256}\n")
        out.write(f"# converter = {CONVERTER_VERSION}\n")
        out.write("# units = Lambda/n_H^2 [erg cm^3 s^-1]; sign: positive = net cooling (WSS09 Net_Cooling, unchanged)\n")
        out.write("# electron_reference = Solar/Electron_density_over_n_h (WSS09 eq. 4 denominator)\n")
        out.write("# excluded = Metal_free (H, He, Compton: owned by the CHUI network; never converted)\n")
        out.write("# assumptions = optically thin, ionization equilibrium, dust free (WSS09 sect. 2)\n")
        for z, path, version, note in files:
            out.write(f"# source z={z:g} file={path.name} sha256={sha256_file(path)} version={version} note={note}\n")
        out.write("elements " + " ".join(s for _, s in groups) + "\n")
        out.write("log10_T " + fmt(log_t) + "\n")
        out.write("log10_nH " + fmt(log_n) + "\n")
        out.write("redshift " + fmt(z for z, _, _, _ in files) + "\n")
        out.write("reference_mass_fraction " + fmt(solar[g] for g, _ in groups) + "\n")
        out.write("reference_hydrogen_mass_fraction " + fmt([solar["Hydrogen"]]) + "\n")
        stream_block(out, "ne_over_nH", "Solar/Electron_density_over_n_h")
        for group, symbol in groups:
            stream_block(out, f"lambda {symbol}", f"{group}/Net_Cooling")
    print(f"wrote {args.output} ({len(files)} redshifts, {len(log_t)} x {len(log_n)} grid, "
          f"{len(groups)} elements); sha256 {sha256_file(args.output)}")
    print("Set physics.metal_cooling_table_path/_sha256 accordingly. The tables embed "
          f"'{args.uv_background}': see docs/chemistry_rate_provenance.md for the UVB compatibility policy.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
