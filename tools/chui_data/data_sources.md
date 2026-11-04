# CHUI physics-data sources (Part 2.5 record)

This file records where each external data product came from during the
Part 2.5 scientific-provenance pass on 2026-11-03: its URL, checksum,
licence status and the format facts that were checked. **No third-party
data file is committed to this repository.** None of the sources below
carries an explicit redistribution licence that covers the data files
themselves, except the NuPyCEE repository's BSD-3-clause licence. The owner
downloads each file and records the checksum in the run configuration
(`*_sha256` keys).

The converters and the validator in this directory were written but **NOT
RUN** by the author. All executable validation is user-reserved.

## 1. UV background: Faucher-Giguère (2020), TREECOOL

| Item | Value |
|---|---|
| Publication | Faucher-Giguère 2020, MNRAS 493, 1614 (arXiv:1903.08657) |
| Original host | https://galaxies.northwestern.edu/uvb (DNS failure/403 from this environment; not reachable) |
| Mirror used | `https://raw.githubusercontent.com/yipihey/ChemistryKernels.jl/main/data/` (README also mirrored by `jeonggyukim/pyathena`) |
| `fg20_treecool_eff_rescaled_heating_rates_068.dat` | 195 rows; sha256 `4bfc5d9bd6f0541addd01256ac20b5e6bc8a7aa5fc8849b9cbdf34cb0882b313`. **Recommended**: the README says the heating rates are rescaled by 0.68 to match the IGM temperature (Gaikwad et al. 2020). |
| `fg20_treecool_eff_default.dat` | 195 rows; sha256 `c6e6030be92c6965c019c53915f31661782eb2cafed058cf9672b7efefc24aaa` |
| Columns (README, verbatim) | `log(1+z), Gamma_HI, Gamma_HeI, Gamma_HeII, qdot_HI, qdot_HeI, qdot_HeII`. Gamma in s⁻¹, qdot in erg s⁻¹. This is identical to the CHUI 7-column reader (`UvBackgroundTable::loadTreecool`). |
| Redshift convention | First column is log10(1+z), 0.000–0.970 (z 0–8.33). CHUI returns zero rates above the last row. |
| Placeholders | 1e-40 (Gamma_HeII) and 6.8e-41 (qdot_HeII) near z_max. They are legal and positive, and log interpolation keeps them negligible. |
| Licence | **None stated** in the README or on the mirrors. **Do not commit.** |
| Cross-check | The ratio of the two files' z=0 qdot_HI values is 1.3784e-25 / 2.027e-25 = 0.680, consistent with the README. |

Validator: `validate_treecool.py <file> [--expect-sha256 ...]` (NOT RUN).

## 2. Metal cooling: Wiersma, Schaye & Smith (2009)

| Item | Value |
|---|---|
| Publication | Wiersma, Schaye & Smith 2009, MNRAS 393, 99 (arXiv:0807.3748) |
| Original host | https://www.strw.leidenuniv.nl/WSS08/. The index page is readable, but its tar archives (`coolingtables.tar.gz`, `z_*.tar.gz`) returned HTTP 403 here. |
| Mirror used | `https://raw.githubusercontent.com/jeonggyukim/pyathena/master/data/microphysics/Wiersma09/tables/z_0.000.hdf5` |
| `z_0.000.hdf5` | 4 741 335 bytes; sha256 `ef31e9e4c5cc11e814a89d638e9af9a37a0f1c627bc7e7575773f5ae3906c968`; Header/Version `100422` |
| Layout (read with h5py) | Per-element groups `Carbon`, `Nitrogen`, `Oxygen`, `Neon`, `Magnesium`, `Silicon`, `Sulphur`, `Calcium`, `Iron`, each holding `Net_Cooling` (352×81, float32), `Temperature_bins` (352, linear K, 100 K–1e9 K) and `Hydrogen_density_bins` (81, linear cm⁻³, 1e-8–1). `Solar/Electron_density_over_n_h` (352×81). `Metal_free/{Net_Cooling, Electron_density_over_n_h, Mean_particle_mass}` (7×352×81, He mass fraction 0.238–0.298). `Total_Metals/Net_cooling`. `Header/Abundances/{Abund_names, Solar_mass_fractions, Solar_number_ratios}` (11 elements; X_H = 0.70649785). `Header/Note` = "Tables calculated with a Haardt and Madau (2001)". `Header/Reference` = WSS09. |
| Units and sign | Λ/n_H² [erg cm³ s⁻¹]. Positive is net cooling. At 100 K and n_H = 1e-8 the values are **negative** (photoheating), so the sign must be preserved. CHUI interpolates in asinh space and never takes a log of a signed quantity. |
| Normalization | WSS09 eq. 4: Λ = Λ_H,He + Σ_i Λ_i,⊙ [(n_e/n_H)/(n_e/n_H)⊙] 10^[i/H], with 10^[i/H] = (n_i/n_H)/(n_i/n_H)⊙. The electron density is tabulated without metal electrons. |
| Assumptions | HM01 (quasars + galaxies) plus CMB photoionization, CLOUDY 07.02, optically thin, ionization equilibrium, dust free. |
| Licence | **None stated** on the WSS08 page or the mirror. **Do not commit.** |

**Redshift coverage restriction (Part 2.5 interpolation-safety correction).**
Only `z_0.000.hdf5` was obtained in this pass. A table converted from it
alone has a **single redshift (z = 0)** and is **not** a cosmological cooling
table.

- The CHUI reader treats every singleton axis (T, n_H or z) as one
  coordinate, matched within 1e-8 × max(1, |value|).
- Under `metal_cooling_table_range = reject`, any other redshift (for
  example z = 3) is **refused**.
- Under `clamp_density_redshift`, the z = 0 rates are used and every such
  evaluation is **reported** as a clamp (thermochemistry ledger
  `metal_table_clamps`, logged per source stage).
- Temperature is never clamped.
- Before this correction a singleton axis matched any query, so z = 0 rates
  were applied silently at all redshifts.
- The converter refuses a one-redshift output unless
  `--allow-single-redshift` is given, and it records `redshift_coverage` in
  the header.
- Cosmological runs need the full WSS09 redshift set (z = 0–8.989, one HDF5
  file per redshift). Those files were not obtainable here.

Converter: `convert_wss09_metal_cooling.py` (NOT RUN). It converts
per-element `Net_Cooling` and `Solar/Electron_density_over_n_h`, and never
`Metal_free` (CHUI owns H, He and Compton).

**UVB compatibility.** The WSS09 tables embed HM01. No TREECOOL file
verified to reproduce that spectrum was available. A WSS09 table combined
with an FG20 TREECOOL is therefore a *declared mismatch*, refused unless
`physics.metal_cooling_accept_uv_background_mismatch = true` (see
`docs/chemistry_rate_provenance.md` §5).

Alternatives found but **not converted**:

- Ploeckinger & Schaye 2020 (FG20 UVB, self-shielding, dust and cosmic rays).
  The tables include their own H/He and ISM shielding assumptions, so their
  semantics differ from a per-element optically-thin increment. The converter
  is blocked until the element-resolved quantities and their shielding
  convention are specified.
- Gnat & Ferland 2012 CIE/ion-by-ion data (collisional only).

## 3. Stellar yields: NuPyCEE / NuGrid

| Item | Value |
|---|---|
| Repository | https://github.com/NuGrid/NuPyCEE (`LICENSE`: BSD 3-clause, "Copyright (c) 2014 - 2016, NuGrid Team"; sha256 of LICENSE `1ef5f600…a3b5`) |
| `yield_tables/agb_and_massive_stars_nugrid_MESAonly_fryer12delay.txt` | 682 225 bytes; sha256 `f2edf5ea29ed12d602022ae3605892b7cd0d4642eaf46945d813b2bdc9f002a1`. "NuGrid yields Set1 with set1.2, set1.1, set1.3a, set1.4a, set1.5a", Fryer et al. 2012 delayed remnants. 60 models: M = 1, 1.65, 2, 3, 4, 5, 6, 7, 12, 15, 20, 25 M☉ × Z = 0.02, 0.01, 0.006, 0.001, 0.0001. |
| Format facts (read) | `H Table: (M=..,Z=..)`, `H Lifetime`, `H Mfinal`, `H CCSN energy`, then `&Isotopes &Yields &X0 &Z &A` rows. Yields are **gross** ejected isotope masses (M☉) and X0 is the initial mass fraction. Check: Σ yields = M − M_final within 2.1e-4 and Σ X0 = 1 to ≤ 6e-5 for all 60 models. |
| Net-yield residual | Net yields computed as yield − X0·(M − M_final). The 11 CHUI elements' net sum is −(untracked net). The untracked part is ≤ 2.5e-4 of the ejecta for AGB models and up to **1.66e-2** for massive models (Z=0.01, 20 M☉), dominated by Ar, Ni, Cr and Ti. The CHUI loader default tolerance (1e-3) therefore refuses a direct conversion. The owner must declare `net_yield_sum_tolerance ≈ 0.02` explicitly. Nothing is rebalanced. |
| `yield_tables/sn1a_i99_W7.txt` | sha256 `352bb0b3d885c411a217f8242c7be792f5055d7444c2c802ed6601b11a33b39f`. Iwamoto et al. 1999 W7, single Z = 0.02 column, isotope masses in M☉. |
| Licence | Repository BSD-3-clause. The yields are from Pignatari et al. 2016 (ApJS 225, 24), Ritter et al. 2018 (MNRAS 480, 538) and Iwamoto et al. 1999 (ApJS 125, 439), which must be cited. Redistribution is permitted under BSD-3 with notice, but the files are still **not committed** (size, and no need). |

Converter: `convert_nugrid_yields.py` (NOT RUN). Massive models go to the
`agb` channel when M < `--cc-mass-min` and to `cc` otherwise. The grid must
be rectangular. SN Ia untracked mass requires `--ia-untracked exclude` (no
default). The choice and the untracked masses are written to the header.

Scientific limitations of this set in CHUI:

- Masses above 25 M☉ and below 1 M☉ are clamped. The clamps are counted
  (`table_clamps`).
- CHUI's own remnant model is not used in table mode. The table's
  M_final (Fryer 2012 delayed) defines the ejecta.
- The 7–8 M☉ and 8–12 M☉ gaps take the nearest model.

## 2026-11-04 acquisition and conversion (D1–D5, final data closure)

The data in this section was downloaded and converted by the agent. Only the
permitted Python converters and validators were executed. CHUI was not built
or run. The workflow is `tools/chui_data/fetch_and_convert.sh <dir>`, which
pins every source SHA-256. A second independent run produced bit-identical
derived files. No third-party data is committed.

### Checksum manifest

| Product | Source / mirror | Source sha256 | Derived file | Derived sha256 |
|---|---|---|---|---|
| NuGrid massive + AGB (Fryer12 delay) | raw.githubusercontent.com/NuGrid/NuPyCEE/master/yield_tables/ | `f2edf5ea…f002a1` | `nugrid_fryer12delay_w7_chui_ssp_yield_table_v1.txt` | `e2161646eadc4119e1b4976045fd6dad7876c7b4940a0c9174e94e1fa56e8f14` |
| SN Ia W7 (Iwamoto 1999) | same repository | `352bb0b3…b33b39f` | (same file, channel `ia`) | — |
| Starburst99 winds, Z = 0.014 Geneva v00, 10⁶ M☉ burst | github.com/jeonggyukim/pyathena `data/sb99/Z014_M1E6_GenevaV00_logdt_10Gyr` (MIT); stsci.edu blocked (proxy 403) | power1 `24aaa7b2…a0397f`, input1 `33f66144…6bac854` | `sb99_z014_genevav00_chui_ssp_wind_table_v1.txt` | `ba232a1cfe2c1996cadf55f14df2476280f8f56e81fe25e2e06c90464b655da5` |
| FG20 TREECOOL (×0.68 heating) | ChemistryKernels.jl mirror | `4bfc5d9b…882b313` | used as is | `4bfc5d9bd6f0541addd01256ac20b5e6bc8a7aa5fc8849b9cbdf34cb0882b313` |
| WSS09 cooling, 49 redshifts z = 0–8.989 (+ collis, photodis, nocompton, not converted) | virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/EAGLE/coolingtables.tar.gz (strw.leidenuniv.nl 403) | tarball `228ec92871e8170a181ca5fb80a591b33daec4475a30b9482419f127670a1543` (per-file sha256 in the table header) | `wss09_hm01_z0-8.989_chui_metal_cooling_table_v1.txt` (52 MB) | `4bbca8a9e5123226392f077f852b89ec551ef4d233b366b3d9c8ee443b143410` |

### D1 yields: results of `validate_ssp_yield_table.py` (VALID)

**Grid.**
- AGB channel: 1, 1.65, 2, 3, 4, 5, 6, 7 M☉.
- CC channel: 12, 15, 20, 25 M☉.
- Metallicities: Z = 1e-4, 1e-3, 6e-3, 0.01, 0.02.
- Ia: Z = 0.02 only; it is used at every Z, so Ia metallicity dependence is
  absent.
- Masses of 8–12 M☉ take the nearest model, and masses outside 1–25 M☉ are
  clamped and counted.

**Conservation checks.**
- Ejecta equal M − M_final, and match the summed gross yields within
  2.1e-4 (relative).
- The source's own non-closure |Σ gross − ej·ΣX0|/ej is ≤ 2.5e-4, from
  rounded yields and ΣX0 ≠ 1.
- The exact identity tracked + untracked net = Σ gross − ej·ΣX0 holds for all
  60 models.
- |Σ tracked net|/ej ≤ 1.656e-2, against a declared tolerance of 0.02 (the
  loader limit is 0.05). The worst case is the 20 M☉, Z = 0.01 model, driven
  by Ar, Ni, Cr and Ti.
- The minimum tracked gross element ejecta is 5.6e-8 M☉; no element is
  negative.

**Untracked-mass policy (explicit).**
- CHUI deposits only the 11 tracked elements, plus dust from the tracked
  gross yields. The star loses exactly the deposited mass. The untracked
  gross ejecta (Na, Al, Ar, Ni, s-process and others) therefore stay locked
  in the star particle.
- Massive and AGB stars: that locked fraction is a median of 2e-4 and a
  maximum of 1.68e-2 of the ejecta.
- SN Ia: the W7 total ejecta are 1.368 M☉. Of that, 0.162 M☉ (**11.8%**) is
  untracked and excluded: Ni 0.126 (mainly ⁵⁸Ni), Ar 0.016, Mn 0.009,
  Cr 0.008, Co 0.001.
- Fe (0.746 M☉) and the other tracked elements are exact.

**Scientific impact.**
- The Ia contribution to total metallicity is low by about 12% by mass.
- Ni, Mn and Cr enrichment is absent.
- Nothing is moved into H or renormalized.

### D2 winds: SB99 → `chui_ssp_wind_table_v1`

The momentum used is SB99's own population momentum flux Σ Ṁ v∞ (its
mass-loss and terminal-velocity prescriptions with wind model 0). It is
integrated in time and is **not** inferred from the energy. The energy is
SB99's cumulative wind energy column. A trapezoid integration of the wind
power reproduces it to 0.78%. Supernovae are excluded, because they are in
`.snr1`. Wind mass is not used, because CHUI's SSP return already contains it.

| Quantity | Value | Basis |
|---|---|---|
| Cumulative energy at 10 Gyr | 1.61e48 erg/M☉ | SB99 cumulative wind energy |
| Cumulative momentum at 10 Gyr | 77.5 M☉ km/s/M☉ | Integrated momentum flux |
| Effective speed 2E/P | about 2,070 km/s | Ratio of the two (consistency check only) |

Both cumulative budgets are non-negative and non-decreasing, over 1001 ages
from 10⁴ to 10¹⁰ yr.

**Limitations.**
- The SB99 IMF is Kroupa (0.1–0.5–100, slopes 1.3 and 2.3), not the CHUI
  Chabrier IMF. The normalization is per M☉ of the SB99 population.
- There is a single metallicity (Z = 0.014). The wind reader now
  **reports** every use at another Z as a metallicity clamp; before this
  change, single-Z use was silent.
- The SB99 deck's model-designation string says "M1E5", but TOMA = 1.000
  (×10⁶ M☉) and the directory name both give 10⁶ M☉. The normalization uses
  TOMA.
- More metallicities require additional SB99 runs, which the owner must run
  or obtain.

### D3 cooling and UVB

**WSS09 table conversion.**
- 49 redshifts, 176 T × 41 n_H. This is the standard-resolution set; the
  earlier pyathena z = 0 file was the 352 × 81 variant, which explains its
  different checksum.
- Values round-trip bit-exactly against float32.
- The table holds the nine metals (C N O Ne Mg Si S Ca Fe).
- `Metal_free` is not converted, so the H/He and Compton terms are owned by
  the CHUI network and nothing is double counted.
- The CMB is included in the WSS09 photoionization; the per-element
  increments contain no Compton term.

**Per-element additivity (informational).**
- The sum of the nine per-element increments differs from WSS09
  `Total_Metals`. For T > 10⁴ K and |Λ| > 1e-23 the median difference is
  3–7% and the 95th percentile is 24–51%.
- This is non-additivity intrinsic to the difference-run definition, and is
  recorded for qualification.

**UVB compatibility: unresolved.**
- WSS09 embeds HM01. No HM01 H/He photorate (TREECOOL-type) table was found.
- FG20 TREECOOL plus WSS09 therefore remains a declared mismatch, refused
  unless `physics.metal_cooling_accept_uv_background_mismatch = true`, which
  marks an acknowledged experiment.
- A consistent pairing requires either:
  - HM01 photoionization and photoheating rates derived from the CUBA HM01
    spectra, or
  - element tables recomputed with FG20.
- Both are outside this session; the second is not available.

**TREECOOL.** `validate_treecool.py`: VALID. 195 rows, z = 0–8.333,
qdot_HI/Γ_HI = 2.38 eV at z = 0.

### Example configurations

The decks are `configs/chui_final/*.param.txt`. They contain these data paths
and checksums, and expect the derived files in `chui_data/` next to the run
directory.
