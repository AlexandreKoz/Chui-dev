# Tools

Place standalone command-line tools here. Keep ownership explicit; avoid mixing solver logic into tools.

- `convert_ic_manifest.py`: emit JSON manifest describing supported GADGET/AREPO IC aliases and required header fields.

- `treepm_sweep.py`: frozen-source TreePM export and split-aware diagnostic sweeps.
- `generate_clustered_dmo_box.py`: synthetic clustered DMO box (NFW-like halos
  plus uniform background) written as the `bench_tree_pm_sweep --particles`
  CSV. `--order random|morton` selects row memory order. Workload fixture for
  interaction-count/timing comparisons only; not a cosmological IC.
- `hierarchical_dmo_qualification.py`: matched global/hierarchical/refined DMO
  trajectory and low-k growth/power comparison.
- `chui_snapshot.py`: shared, offline, Python 3.8-compatible complete logical
  science snapshot reader used by both tools; not a simulation-input parser.

The two qualification workflows require Python >=3.8, NumPy and h5py compatible
with that interpreter. This does not redefine the minimum version of unrelated
bindings or repository tools. Pass a `.complete` marker, one member of its set,
or a directory containing exactly one marker. Modern CHUI v6 science requires
its production v2/v3 completion manifest, even for one HDF5 file. All listed
members, counts, cosmology, rectangular geometry, units, generation and IDs are
validated before export. v2 member checksums are enforced; v3 readback evidence
is honored and current file checksums are recorded separately. Checkpoints are
refused. External/CHUI <=v5 legacy input requires `--legacy-dialect`; absent
units also require the complete `--legacy-length-unit`, `--legacy-mass-unit`,
`--legacy-velocity-unit` triple. Unsupported revisions and ambiguous units fail.

Export defaults to declared stored-to-internal dialect conversions: interoperable
positions/masses divide by h, native values retain their declared code units.
`--position-scale`/`--mass-scale` explicitly override these factors. Box axes,
source/member/marker checksums, ID fingerprints, epoch, units, config identifiers
and conversion factors are saved in the CSV `.metadata.json` sidecar. Supply the
original comoving `--epsilon` in those solver length units; no mesh-dependent
softening change is implied. `--fixed-target-ids FILE` preserves an existing
oracle's persistent-ID order; `--target-ids FILE` names the exported list.

External `--total-reference` comparisons require `--total-reference-manifest`
with matching CSV/source/target hashes, epoch, box, epsilon and G, a reference
kind of `independent_periodic_ewald` or `converged_total`, and explicit
`qualification_evidence`. A short-range oracle must match each candidate split;
it cannot be reused as a total-force reference. Hierarchical comparisons require
matching normalized scientific provenance, ignoring only documented execution,
output and time-resolution keys. They report the declared inverse stored-length
k units and retain the quarter-Nyquist estimator restriction.

Reader staging is chunked (65,536 rows by default). Only selected-species
position/mass/(optional velocity) lanes persist; one all-species uint64 ID lane
is sorted to reject cross-member/type collisions and freed before selected
sorting. These are external-tool allocations, not production MemoryGovernor
buffers. Large snapshots still require memory for the selected state and FFT
estimator; no large-problem fit is claimed.

Manual Python regression (not executed during this repair):

```bash
python3.8 tests/python/test_treepm_snapshot_tools.py
```

Missing NumPy/h5py gives explicit exit 77 (CTest SKIP); all data/semantic
rejections use assertions on named error categories. See the
[repair report](../docs/repair/treepm_p1p5_qualification_and_test_reliability_repair_20261010.md) for the full matrix and commands.

## CHUI physics-data tools (`tools/chui_data/`, Part 2.5)

These tools were written but **not executed** by the author. Executable
validation is user-reserved.

- `validate_treecool.py`: checks a TREECOOL file against the CHUI 7-column
  reader contract and prints its SHA-256 and config lines.
- `convert_wss09_metal_cooling.py`: streams the Wiersma, Schaye & Smith (2009)
  per-redshift HDF5 tables into `chui_metal_cooling_table_v1`. It needs
  h5py, preserves the sign and axes, never converts `Metal_free`, and
  requires `--license` and `--uv-background`.
- `convert_nugrid_yields.py`: converts NuPyCEE/NuGrid yields into
  `chui_ssp_yield_table_v1` (net = yield − X0·ejecta). It refuses rather
  than rebalances and requires an explicit Ia untracked-mass policy.
- `data_sources.md`: URLs, checksums, licence status and the format facts
  checked for each source. No third-party data is committed.
