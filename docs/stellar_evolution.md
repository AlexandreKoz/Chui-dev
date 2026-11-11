# Stellar evolution and enrichment bookkeeping

The production implementation and scientific rationale are documented in
[`metals_enrichment_and_mixing.md`](metals_enrichment_and_mixing.md).

## Configuration

Under `[physics]`:

- `enable_stellar_evolution` enables delayed SSP bookkeeping and source-runtime deposition;
- `stellar_evolution_table_path` selects a v2 age-by-birth-metallicity table;
- `stellar_evolution_require_production_table=true` rejects the built-in zero-yield compatibility table;
- `stellar_evolution_hubble_time_years` is retained only for configuration compatibility and legacy isolated API calls. Production cosmological ages use the FLRW integral.

## V2 table rows

Each whitespace-delimited row contains 22 values:

1. birth metallicity;
2. SSP age in years;
3. cumulative returned mass fraction;
4. cumulative total ejected metal fraction;
5. cumulative newly synthesized metal fraction;
6. cumulative event count per initial code mass;
7. cumulative energy in erg per initial code mass;
8-10. returned mass fractions for winds/AGB, CCSN, SNIa;
11-13. total ejected metal fractions by channel;
14-16. newly synthesized metal fractions by channel;
17-19. event counts by channel;
20-22. energies by channel.

Metadata comments use `# key = value`; supported keys include table ID/version,
source papers/repository, redistribution license, SHA-256, IMF, stellar mass
range, solar abundance reference, and `production_calibrated`.

The deterministic test table is `resources/stellar_evolution/test_synthetic_v2.txt`.
It is synthetic and non-production. `resources/stellar_evolution/default_v1.txt`
is retained only as historical source data and is not accepted by the v2 production loader.

## Part 2 SSP yield-table coverage and wind-table compatibility (patch 58)

The Part 2 SSP (`stellar_population_model = chabrier_yield_table_v1`,
`SspModel` in `src/physics/stellar_population.cpp`) integrates AGB returns over
`[ssp_imf_mass_min_msun, min(ssp_imf_mass_max_msun, ssp_cc_mass_min_msun)]`
and CC returns over `[ssp_cc_mass_min_msun, ssp_cc_mass_max_msun]` (capped by
the IMF limits) on fixed quadrature nodes from the turnoff mass upward, at the
star's birth metallicity; Ia uses a per-event, metallicity-only row.
`physics.ssp_yield_table_coverage` decides what happens outside the tabulated
`chui_ssp_yield_table_v1` grid:

* `reject` (default): construction refuses a table whose mass axes do not
  contain the AGB upper end and the whole CC range, or whose channel has a
  single metallicity; at use, any remaining out-of-grid query (AGB turnoff
  below the tabulated minimum mass — conservatively up to one quadrature node,
  ~0.03 dex, early — or a birth metallicity outside the grid) throws inside
  the stellar-population `StageTransaction`, so the stage is rejected
  collectively and no return, event ordinal, cursor or ledger change commits.
  Queries within 1e-12 (relative) of an axis end are on the grid.
* `edge_clamp_research`: the historical edge row (no extrapolation), counted
  per query (`yield_table_clamps`) and per axis
  (`yield_table_mass_clamps`, `yield_table_metallicity_clamps`) and recorded.

Wind tables (`ssp_wind_model = table_v1`) are checked against the CHUI SSP
IMF identity `chabrier2003_individual` and the configured IMF limits under
`physics.ssp_wind_table_compatibility` (`strict` default /
`approximate_research`); see `include/cosmosim/physics/stellar_wind.hpp` for
the `chui_ssp_wind_table_v2` header contract. The wind channel still returns
energy and momentum only; its mass remains in the SSP CC return.
