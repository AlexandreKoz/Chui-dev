#!/usr/bin/env python3
"""Derive TREECOOL-format H/He photoionization and photoheating rates from a
Haardt & Madau UV-background SPECTRUM, by direct frequency integration.

Purpose (C01 UVB consistency). The WSS09 element cooling tables were computed
with CLOUDY 07.02 photoionized by the Haardt & Madau galaxies+quasars
background (CLOUDY ``table HM05``, data file ``haardt_madau_galaxy.dat``; the
revised HM01 continuum supplied by F. Haardt in 2005, see the file header).
CHUI takes its H/He photo-rates from a TREECOOL file, so a CONSISTENT pairing
needs the H/He rates of that same spectrum. This tool computes them; it never
rescales, smooths or extrapolates the source spectrum.

Inputs (``--format``):
  * ``cloudy_hm05``: CLOUDY data file (``#`` comments, magic number 50808, a
    ``wavelength(A) z=... z=...`` header line, then rows lambda[A] J_nu(z_k)).
  * ``cuba_uvb``: CUBA ``UVB.out`` (``#`` comments, one line of redshifts,
    then rows lambda[A] J_nu(z_k)). Used to VALIDATE the method against the
    authors' own ``photorates.out`` (``--reference-photorates``).
  J_nu is in erg s^-1 cm^-2 Hz^-1 sr^-1 at the rest-frame wavelength.

Rates (optically thin, per ion):
  Gamma_i = 4 pi \\int_{nu_i} J_nu sigma_i(nu) / (h nu) d nu           [s^-1]
  qdot_i  = 4 pi \\int_{nu_i} J_nu sigma_i(nu) (h nu - h nu_i)/(h nu) d nu [erg s^-1]
  i = H I, He I, He II; sigma_i from the analytic fits of Verner, Ferland,
  Korista & Yakovlev 1996 (ApJ 465, 487), Table 1, ground-state shells. The
  integral runs from each threshold to min(E_max of the fit, table maximum).

Numerics: a repeated wavelength is a spectral edge (step). In ln E, the
spectrum is interpolated log-log between tabulated points (linearly where an endpoint is zero, so zeros stay zero), sigma is
evaluated exactly, each native interval is integrated with composite Simpson
(``--substeps`` points), and the threshold is inserted exactly. The tool
re-integrates at 4x resolution and refuses output whose relative change
exceeds ``--convergence-tolerance`` (default 1e-4).

Output contract (src/physics/radiation_background.cpp,
UvBackgroundTable::loadTreecool): ``#`` provenance header, then rows
  log10(1+z) Gamma_HI Gamma_HeI Gamma_HeII qdot_HI qdot_HeI qdot_HeII
strictly increasing in log10(1+z), all finite and >= 0. CHUI interpolates
log10 of positive rates and returns zero above the last row.

Standard library only; run with ``python3 -I``. Nothing here builds CHUI.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

H_PLANCK_ERG_S = 6.62607015e-27
EV_ERG = 1.602176634e-12
HC_EV_ANGSTROM = 12398.419843320026
MAGIC_HM05 = 50808

# Verner et al. 1996 Table 1: E_th, E_max, E0 [eV], sigma0 [Mb], ya, P, yw, y0, y1.
VERNER = {
    "HI": (13.60, 5.0e4, 0.4298, 5.475e4, 32.88, 2.963, 0.0, 0.0, 0.0),
    "HeI": (24.59, 5.0e4, 13.61, 949.2, 1.469, 3.188, 2.039, 0.4434, 2.136),
    "HeII": (54.42, 5.0e4, 1.720, 1.369e4, 32.88, 2.963, 0.0, 0.0, 0.0),
}
SPECIES = ("HI", "HeI", "HeII")
# Threshold cross sections implied by the fits (Mb); a self-check against
# transcription errors in the table above.
THRESHOLD_SIGMA_MB = {"HI": 6.35, "HeI": 7.44, "HeII": 1.59}


def sigma_cm2(species: str, energy_ev: float) -> float:
    e_th, e_max, e0, s0, ya, p, yw, y0, y1 = VERNER[species]
    if energy_ev < e_th or energy_ev > e_max:
        return 0.0
    x = energy_ev / e0 - y0
    y = math.sqrt(x * x + y1 * y1)
    f = ((x - 1.0) ** 2 + yw * yw) * y ** (0.5 * p - 5.5) * (1.0 + math.sqrt(y / ya)) ** (-p)
    return s0 * 1.0e-18 * f


def self_check_cross_sections() -> None:
    for species, expected in THRESHOLD_SIGMA_MB.items():
        got = sigma_cm2(species, VERNER[species][0]) / 1.0e-18
        if abs(got - expected) > 0.01 * expected:
            raise SystemExit(f"cross-section self-check failed for {species}: {got:.4f} Mb != {expected} Mb")


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_spectrum(path: Path, fmt: str) -> tuple[list[float], list[float], list[list[float]]]:
    """Returns (redshifts, wavelengths_A ascending, J[z][lambda])."""
    rows: list[list[float]] = []
    redshifts: list[float] | None = None
    magic_seen = False
    with path.open("r", encoding="latin-1") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if fmt == "cloudy_hm05":
                if not magic_seen:
                    if int(line.split()[0]) != MAGIC_HM05:
                        raise SystemExit(f"{path}: magic number {line.split()[0]} != {MAGIC_HM05}")
                    magic_seen = True
                    continue
                if redshifts is None:
                    fields = line.split()
                    if not fields[0].lower().startswith("wavelength"):
                        raise SystemExit(f"{path}: expected the 'wavelength(A) z=...' header line")
                    redshifts = [float(f.split("=", 1)[1]) for f in fields[1:]]
                    continue
            elif fmt == "cuba_uvb":
                if redshifts is None:
                    redshifts = [float(f) for f in line.split()]
                    continue
            else:
                raise SystemExit(f"unknown --format {fmt}")
            values = [float(f) for f in line.split()]
            rows.append(values)
    if redshifts is None or not rows:
        raise SystemExit(f"{path}: no spectrum found")
    nz = len(redshifts)
    for k, values in enumerate(rows):
        if len(values) != nz + 1:
            raise SystemExit(f"{path}: data row {k} has {len(values) - 1} intensities, expected {nz}")
        if any(not math.isfinite(v) or v < 0.0 for v in values):
            raise SystemExit(f"{path}: data row {k} has a non-finite or negative value")
    wavelengths = [r[0] for r in rows]
    edges = 0
    for k in range(1, len(wavelengths)):
        if wavelengths[k] < wavelengths[k - 1]:
            raise SystemExit(f"{path}: wavelengths decrease at row {k} ({wavelengths[k]})")
        if wavelengths[k] == wavelengths[k - 1]:
            # A repeated wavelength encodes a discontinuity (absorption edge):
            # the zero-width interval contributes nothing and each side keeps
            # its own tabulated value.
            edges += 1
    if edges:
        print(f"{path.name}: {edges} repeated wavelength(s) treated as spectral edges")
    if any(not redshifts[k] > redshifts[k - 1] for k in range(1, nz)):
        raise SystemExit(f"{path}: redshifts not strictly increasing")
    intensities = [[rows[i][k + 1] for i in range(len(rows))] for k in range(nz)]
    return redshifts, wavelengths, intensities


def interp_loglog(e: float, e_a: float, e_b: float, j_a: float, j_b: float) -> float:
    if e_b == e_a:
        return j_a
    t = (math.log(e) - math.log(e_a)) / (math.log(e_b) - math.log(e_a))
    if j_a > 0.0 and j_b > 0.0:
        return math.exp((1.0 - t) * math.log(j_a) + t * math.log(j_b))
    return (1.0 - t) * j_a + t * j_b


def rates_for(energies: list[float], j: list[float], species: str, substeps: int) -> tuple[float, float]:
    """energies ascending [eV], j aligned. Returns (Gamma [s^-1], qdot [erg s^-1])."""
    e_th, e_max = VERNER[species][0], VERNER[species][1]
    upper = min(e_max, energies[-1])
    if substeps % 2 == 1:
        substeps += 1
    gamma = 0.0
    heat = 0.0
    for k in range(len(energies) - 1):
        e_a, e_b = energies[k], energies[k + 1]
        if e_b == e_a:
            continue  # spectral edge (zero-width interval)
        lo, hi = max(e_a, e_th), min(e_b, upper)
        if not hi > lo:
            continue
        # Composite Simpson in ln E over [lo, hi] (exact threshold insertion).
        u_lo, u_hi = math.log(lo), math.log(hi)
        h = (u_hi - u_lo) / substeps
        g_sum = 0.0
        q_sum = 0.0
        for m in range(substeps + 1):
            weight = 1.0 if m in (0, substeps) else (4.0 if m % 2 == 1 else 2.0)
            e = math.exp(u_lo + m * h)
            jj = interp_loglog(e, e_a, e_b, j[k], j[k + 1])
            s = sigma_cm2(species, e)
            g_sum += weight * jj * s
            q_sum += weight * jj * s * (e - e_th)
        gamma += g_sum * h / 3.0
        heat += q_sum * h / 3.0
    # d nu / (h nu) = d ln E / h ; 4 pi sr.
    prefactor = 4.0 * math.pi / H_PLANCK_ERG_S
    return prefactor * gamma, prefactor * heat * EV_ERG


def all_rates(wavelengths: list[float], intensities: list[list[float]], substeps: int) -> list[list[float]]:
    energies = [HC_EV_ANGSTROM / w for w in reversed(wavelengths)]
    out = []
    for j_by_lambda in intensities:
        j = list(reversed(j_by_lambda))
        gammas, heats = [], []
        for species in SPECIES:
            g, q = rates_for(energies, j, species, substeps)
            gammas.append(g)
            heats.append(q)
        out.append(gammas + heats)
    return out


def read_cuba_photorates(path: Path) -> dict[float, list[float]]:
    """CUBA photorates.out: z HI_i HI_h HeI_i HeI_h HeII_i HeII_h Compton (heating in eV/s)."""
    table: dict[float, list[float]] = {}
    with path.open("r", encoding="latin-1") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            v = [float(f) for f in line.split()]
            if len(v) < 7:
                continue
            table[v[0]] = [v[1], v[3], v[5], v[2] * EV_ERG, v[4] * EV_ERG, v[6] * EV_ERG]
    return table


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("--format", required=True, choices=("cloudy_hm05", "cuba_uvb"))
    parser.add_argument("--output", type=Path, help="TREECOOL file to write (omit in validation mode)")
    parser.add_argument("--expect-spectrum-sha256", default="")
    parser.add_argument("--source", default="", help="provenance text: archive URL and checksum")
    parser.add_argument("--identity", default="", help="UV background identity written to the header")
    parser.add_argument("--license", default="", help="licence statement (never invented)")
    parser.add_argument("--zmax", type=float, default=math.inf, help="drop spectrum redshifts above this")
    parser.add_argument("--substeps", type=int, default=16)
    parser.add_argument("--convergence-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--reference-photorates", type=Path,
                        help="CUBA photorates.out: compare instead of / in addition to writing")
    parser.add_argument("--reference-tolerance", type=float, default=0.05)
    args = parser.parse_args()

    self_check_cross_sections()
    spectrum_sha = sha256_of(args.spectrum)
    if args.expect_spectrum_sha256 and spectrum_sha != args.expect_spectrum_sha256:
        raise SystemExit(f"spectrum sha256 {spectrum_sha} != expected {args.expect_spectrum_sha256}")
    redshifts, wavelengths, intensities = read_spectrum(args.spectrum, args.format)
    keep = [k for k, z in enumerate(redshifts) if z <= args.zmax]
    redshifts = [redshifts[k] for k in keep]
    intensities = [intensities[k] for k in keep]

    rates = all_rates(wavelengths, intensities, args.substeps)
    fine = all_rates(wavelengths, intensities, 4 * args.substeps)
    worst = 0.0
    for coarse_row, fine_row in zip(rates, fine):
        for a, b in zip(coarse_row, fine_row):
            if b > 0.0:
                worst = max(worst, abs(a - b) / b)
    print(f"convergence: max relative change at 4x resolution = {worst:.3e}")
    if worst > args.convergence_tolerance:
        raise SystemExit("integration not converged; increase --substeps")

    if args.reference_photorates is not None:
        reference = read_cuba_photorates(args.reference_photorates)
        names = ["Gamma_HI", "Gamma_HeI", "Gamma_HeII", "qdot_HI", "qdot_HeI", "qdot_HeII"]
        worst_by = [0.0] * 6
        compared = 0
        for z, row in zip(redshifts, rates):
            match = [zr for zr in reference if abs(zr - z) <= 1.0e-3 * max(1.0, z)]
            if not match:
                continue
            ref = reference[match[0]]
            if ref[0] <= 0.0:
                continue
            compared += 1
            for c in range(6):
                if ref[c] > 0.0:
                    worst_by[c] = max(worst_by[c], abs(row[c] - ref[c]) / ref[c])
        print(f"reference comparison over {compared} common redshifts (max relative difference):")
        for name, w in zip(names, worst_by):
            print(f"  {name:10s} {w:.3e}")
        if compared == 0:
            raise SystemExit("no common redshift with the reference")
        if max(worst_by) > args.reference_tolerance:
            raise SystemExit(f"reference comparison exceeds tolerance {args.reference_tolerance}")
        print("reference comparison: PASS")

    if args.output is not None:
        if not args.license:
            raise SystemExit("--license is required when writing (never invented)")
        lines = [
            "# chui_treecool_from_spectrum_v1",
            "# TREECOOL: log10(1+z) Gamma_HI Gamma_HeI Gamma_HeII [s^-1] qdot_HI qdot_HeI qdot_HeII [erg s^-1]",
            f"# uv_background: {args.identity}",
            f"# spectrum_file: {args.spectrum.name} (format {args.format})",
            f"# spectrum_sha256: {spectrum_sha}",
            f"# source: {args.source}",
            f"# license: {args.license}",
            "# method: optically thin per-ion rates, 4 pi int J_nu sigma/(h nu) dnu; Verner et al. 1996 ApJ 465 487",
            "#   cross sections (H I, He I, He II); log-log J interpolation; Simpson in ln E; "
            f"{args.substeps} substeps/interval; 4x convergence {worst:.2e}",
            f"# redshifts: {len(redshifts)} rows, z = {redshifts[0]:g} .. {redshifts[-1]:g}; "
            "CHUI returns zero above the last row",
            "# generator: tools/chui_data/convert_hm_spectrum_treecool.py",
        ]
        for z, row in zip(redshifts, rates):
            lines.append(f"{math.log10(1.0 + z):.6f} " + " ".join(f"{v:.9e}" for v in row))
        args.output.write_text("\n".join(lines) + "\n", encoding="ascii")
        print(f"wrote {args.output} ({len(redshifts)} rows); sha256 {sha256_of(args.output)}")
        z0 = rates[0]
        print(f"z={redshifts[0]:g}: Gamma_HI={z0[0]:.3e} s^-1  qdot_HI/Gamma_HI={z0[3] / z0[0] / EV_ERG:.2f} eV")
        for target in (2.0, 3.0, 6.0):
            k = min(range(len(redshifts)), key=lambda i: abs(redshifts[i] - target))
            r = rates[k]
            print(f"z={redshifts[k]:g}: Gamma_HI={r[0]:.3e} Gamma_HeI={r[1]:.3e} Gamma_HeII={r[2]:.3e} s^-1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
