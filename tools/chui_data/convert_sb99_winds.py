#!/usr/bin/env python3
"""Convert a Starburst99 ``*.power1`` (stellar-wind power, energy and
momentum-flux) output into CHUI's ``chui_ssp_wind_table_v1``
(src/physics/stellar_wind.cpp). D2 data closure, 2026-11-04.

Physics:
  * SB99 ``*.power1`` columns are POPULATION-INTEGRATED quantities of an
    instantaneous burst of total initial mass M_tot (input deck TOMA, in
    1e6 Msun): wind power L_w(t) [erg/s] (log10), cumulative wind energy
    E_w(t) [erg] (log10), wind momentum flux pdot_w(t) [dyn] (log10), each
    total and per stellar class. They are SB99 sums over the track stars of
    0.5 Mdot v_inf^2 and Mdot v_inf with its own mass-loss/terminal-velocity
    prescriptions (IWIND), i.e. an INDEPENDENT momentum prescription (the
    momentum is never inferred from the energy here). Supernovae are NOT in
    ``*.power1`` (they are ``*.snr1``), so there is no double counting with
    the CHUI SN channel. Wind MASS is not used: CHUI returns wind mass through
    the SSP continuous return.
  * Output per unit initial stellar mass of the SB99 population:
      E(t)/M_tot [erg/Msun]  = SB99 cumulative energy column / M_tot;
      P(t)/M_tot [Msun km/s / Msun] = int_0^t pdot dt / M_tot, trapezoidal in
        linear time on the SB99 time grid, plus pdot(t_0) * t_0 for [0, t_0]
        (the reader's linear-ramp convention below the first age).
    The pdot integral is non-decreasing by construction (pdot >= 0); the
    energy column is checked to be non-decreasing. As a quadrature check
    the same trapezoid applied to L_w is compared with the SB99 energy column
    (reported).
  * log10 value -30 is SB99's zero and is converted to 0.
Limitations recorded in the header: the SB99 IMF (from the input deck) may
differ from the CHUI SSP IMF (Chabrier); one metallicity per run (CHUI reports
metallicity extrapolation, see stellar_wind.cpp).
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

MSUN_G = 1.98847e33
SEC_PER_YR = 3.15576e7
KMS = 1.0e5


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_input(path: Path) -> dict:
    lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    out = {}
    for i, line in enumerate(lines):
        for key in ("ISF", "TOMA", "NINTERV", "XPONENT", "XMASLIM", "SNCUT", "BHCUT", "IZ", "IWIND"):
            if f"[{key}]" in line:
                j = i + 1
                while j < len(lines) and (lines[j].strip() == "" or re.match(r"^[A-Z]", lines[j].strip())):
                    j += 1  # skip option-list lines (e.g. track table under [IZ])
                out[key] = lines[j].strip()
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--power1", required=True, type=Path)
    p.add_argument("--input1", required=True, type=Path)
    p.add_argument("--metallicity", required=True, type=float, help="mass-fraction Z of the SB99 tracks used")
    p.add_argument("--citation", required=True)
    p.add_argument("--license", required=True)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    deck = parse_input(a.input1)
    if float(deck.get("ISF", "0")) > 0:
        raise SystemExit("convert_sb99_winds: continuous star formation run; an instantaneous burst is required")
    m_tot = float(deck["TOMA"]) * 1.0e6
    rows = []
    for line in a.power1.read_text(encoding="ascii", errors="replace").splitlines():
        tok = line.split()
        if len(tok) == 12 and re.match(r"^[.0-9]", tok[0]):
            rows.append([float(x) for x in tok])
    if len(rows) < 2:
        raise SystemExit("convert_sb99_winds: no data rows")
    t = [r[0] for r in rows]
    if any(not b > c for c, b in zip(t, t[1:])):
        raise SystemExit("convert_sb99_winds: time axis not increasing")
    lin = lambda v: 0.0 if v <= -29.999 else 10.0 ** v  # noqa: E731
    lw = [lin(r[1]) for r in rows]
    ew = [lin(r[6]) for r in rows]
    pdot = [lin(r[7]) for r in rows]
    # cumulative momentum (g cm/s) and a quadrature check of the energy
    p_cum, e_trap = [pdot[0] * t[0] * SEC_PER_YR], [lw[0] * t[0] * SEC_PER_YR]
    for i in range(1, len(t)):
        dt = (t[i] - t[i - 1]) * SEC_PER_YR
        p_cum.append(p_cum[-1] + 0.5 * dt * (pdot[i] + pdot[i - 1]))
        e_trap.append(e_trap[-1] + 0.5 * dt * (lw[i] + lw[i - 1]))
    for i in range(1, len(ew)):
        if ew[i] < ew[i - 1]:
            raise SystemExit(f"convert_sb99_winds: SB99 cumulative energy decreases at t={t[i]:.4e} yr")
    e_per = [e / m_tot for e in ew]
    p_per = [pc / m_tot / (MSUN_G * KMS) for pc in p_cum]
    quad = abs(e_trap[-1] - ew[-1]) / ew[-1]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open("w", encoding="ascii") as out:
        out.write("# chui_ssp_wind_table_v1\n")
        out.write(f"# source = Starburst99 {a.power1.name} sha256={sha(a.power1)}; deck {a.input1.name} "
                  f"sha256={sha(a.input1)}; {a.citation}\n")
        out.write(f"# license = {a.license}\n")
        out.write("# converter = chui_sb99_wind_converter_v1\n")
        out.write(f"# sb99_total_mass_msun = {m_tot:g} (TOMA); IMF exponents {deck.get('XPONENT')} "
                  f"mass limits {deck.get('XMASLIM')}; SN cut {deck.get('SNCUT')}; BH cut {deck.get('BHCUT')}; "
                  f"tracks IZ={deck.get('IZ')}; wind model IWIND={deck.get('IWIND')}\n")
        out.write("# quantities = cumulative wind energy [erg/Msun] and momentum [Msun km/s per Msun] of the "
                  "SB99 population per unit initial stellar mass (population integrated, winds only, no SNe)\n")
        out.write(f"# energy_quadrature_check = |trapz(L_w) - E_SB99| / E_SB99 at t_max = {quad:.3e}\n")
        out.write(f"# totals_at_t_max = E {e_per[-1]:.6e} erg/Msun, P {p_per[-1]:.6e} Msun km/s/Msun\n")
        out.write("# limitation = SB99 IMF differs from the CHUI Chabrier SSP IMF; single metallicity\n")
        out.write(f"metallicities 1 {a.metallicity:g}\n")
        out.write(f"ages {len(t)}\n")
        for ti, e, pp in zip(t, e_per, p_per):
            from math import log10
            out.write(f"{log10(ti):.8f} {e:.8e} {pp:.8e}\n")
    print(f"wrote {a.output} ({len(t)} ages, {t[0]:.3g}-{t[-1]:.3g} yr); E={e_per[-1]:.3e} erg/Msun "
          f"P={p_per[-1]:.3e} Msun km/s/Msun; quadrature check {quad:.2e}; sha256 {sha(a.output)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
