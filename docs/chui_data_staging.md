# Staging the CHUI external physics data (user-local)

Patch 53. No third-party table is committed to this repository. Each run binds
its tables by path **and** SHA-256 in the `.param.txt` configuration; the
loaders recompute the checksum on every rank and refuse a mismatch, and a
failure on any rank is agreed collectively during runtime-module
construction.

## Products and status (2026-11-09 assessment)

`tools/chui_data/dataset_registry.json` is the machine-readable source of the
table below; `python3 -I tools/chui_data/chui_data.py report <dir>` evaluates
`bytes_present` against a local directory.

| Product | File | SHA-256 | Bytes present | Cross-checked | Rights confirmed | Qualified |
|---|---|---|---|---|---|---|
| S02 NuGrid + W7 yields | `nugrid_fryer12delay_w7_chui_ssp_yield_table_v1.txt` | `e2161646…f8e14` | yes | regenerated from pinned NuPyCEE upstream | no | no |
| F02 SB99 winds Z=0.014 | `sb99_z014_genevav00_chui_ssp_wind_table_v1.txt` | `ba232a1c…5da5` | yes | regenerated from pinned pyathena SB99 outputs | no | no |
| FG20 TREECOOL ×0.68 | `fg20_treecool_eff_rescaled_heating_rates_068.dat` | `4bfc5d9b…b313` | yes | mirror only (original host unreachable) | no | no |
| C01 WSS09 (HM01) | `wss09_hm01_z0-8.989_chui_metal_cooling_table_v1.txt` | `4bbca8a9…3410` | yes, user-local (acquired 2026-11-09) | bit-exact vs source HDF5; pinned tarball hash | no (no licence stated) | no |
| C02 C/O/Si/Fe ions | `chui_metal_ion_rates_v1_C_O_Si_Fe_unverified.txt` | `633a66ee…2d00` | yes | cooling vs GF12 mirror only; CI/RR/DR NOT independently verified | no | no (`unverified`) |
| HM01 TREECOOL | — | — | **no** | — | — | — |

## Procedure

```bash
mkdir -p chui_data
python3 -I tools/chui_data/chui_data.py stage --to chui_data <bundle-dir-or-zip> [...]
python3 -I tools/chui_data/chui_data.py verify chui_data --require-all
python3 -I tools/chui_data/chui_data.py sha256sums chui_data
(cd chui_data && sha256sum -c SHA256SUMS)
```

To regenerate S02, F02, UVB and C01 from their pinned upstream sources
(network, `curl`, `git`, and `h5py` + `numpy` for WSS09):

```bash
python3 -I tools/chui_data/chui_data.py fetch chui_data
```

`stage` never overwrites a file with different bytes and refuses any member
whose SHA-256 differs from the registry. `chui_data/` and `chui_data_local/`
are git-ignored.

## Configuration fragments

* `profiles/SAFE_DEFAULT_AVAILABLE_DATA.param.fragment`: **partial,
  diagnostic**. S02 + F02 with `analytic_reference`, no research opt-ins. The
  FG20 UVB is excluded because every backend that consumes it is gated as
  unqualified; C01 is excluded because no HM01-consistent UVB exists; C02 is
  excluded because it is unverified.
* `profiles/RESEARCH_C02_EXPLICIT_OPT_IN.param.fragment`: **research only**.
  `fixed_detailed` with the C/O/Si/Fe ladders, the checked C02 hash, the
  `thermochemistry_allow_unqualified_rates = true` opt-in, collisional regime
  (`uv_background_model = none`), `metal_cooling_table_range = reject`.

Regenerate both with `chui_data.py emit-config --profile {safe_default,research_c02}`.

## Known data/model mismatches

* WSS09 embeds HM01; FG20 is not HM01. The pair is refused unless
  `metal_cooling_accept_uv_background_mismatch = true` (recorded in the
  normalized configuration and in the `sources.data_assumptions` event).
* SB99 winds assume a Kroupa IMF; the CHUI SSP uses Chabrier. Single Z.
* NuGrid yields cover 1–25 M☉ and Z = 10⁻⁴–0.02; W7 is single-Z; 11.8% of
  the W7 ejecta are untracked.
* C02: 10⁴–10⁸ K, optically thin, collisional; no metal photo-rates; no
  low-T dielectronic resonances; Voronov-rounded potentials. With
  `chemistry_adaptive_selection` an element table is required (demoted rows).
