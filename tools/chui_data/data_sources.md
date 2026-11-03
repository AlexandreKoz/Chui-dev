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
