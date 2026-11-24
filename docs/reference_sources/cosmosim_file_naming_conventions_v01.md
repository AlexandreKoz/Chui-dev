# CosmoSim — **Master File Naming Conventions**
*Version 1.0 (2025-11-06)*

This document defines **authoritative file and directory naming rules** for the CosmoSim codebase and adjacent tooling (C++ engine, YAML configs, Python bindings, HDF5 snapshots, analysis outputs, docs). It complements but does not duplicate variable naming rules.

---

## 0) General Principles
- **Lower snake case for all files and directories:** `lower_snake` (e.g., `tree_builder.cc`, `power_spectrum.py`).  
- **ASCII only, no spaces, no accents, no mixed case.**
- **Use `_` not `-`** except for human-facing docs/images in `docs/` where `kebab-case` is allowed.
- **Stable stems over time.** Prefer changing parameters inside files vs. renaming stems.
- **Semantic and search-friendly.** Include module and role in the name (`pm_solver`, `halo_finder`).
- **Avoid collisions across OSes.** Don’t rely on case sensitivity or reserved names.
- **Line endings:** Use LF (`\n`) across all tracked files; `.gitattributes` enforces this.

---

## 1) Directory Layout (top-level skeleton)
```
cosmosim/
  CMakeLists.txt
  cmake/                 # CMake modules & toolchain files
  external/              # vendored/third-party (read-only)
  include/               # public headers (PascalCase types ok; files still lower_snake.hh)
  src/                   # C++ sources
    gravity/             # tree, pm, treepm
    hydro/               # solver, riemann, limiters
    amr/                 # patch mgmt, prolongation, restriction
    io/                  # hdf5, readers/writers
    analysis/            # in-situ & post-processing kernels
    utils/               # math, allocators, timing, logging
  python/                # pycosmosim bindings & tools
  configs/               # YAML profiles and examples
  tests/                 # unit & integration tests
    unit/
    integration/
  bench/                 # micro & macro benchmarks
  examples/              # runnable minimal examples
  output/                # run artifacts (snapshots, logs, figs)
  docs/                  # user/dev documentation, figures
  scripts/               # CLI helpers, CI scripts
  tools/                 # small utilities, formatters, converters
```

---

## 2) Source Code Files
### 2.1 C++ Engine
- **Source:** `.cc` (or `.cpp`, choose one globally; default: `.cc`)
- **Headers:** `.hh` (internal), `.hpp` (public headers in `include/` only if you prefer; default: `.hh`)
- **Unit under module folder:**  
  `src/<module>/<unit>_<role>.<ext>`

**Examples**
```
src/gravity/pm_solver.cc
src/gravity/tree_builder.cc
src/gravity/tree_node.hh
src/hydro/hydro_solver.cc
src/amr/patch_manager.cc
src/io/hdf5_writer.cc
src/utils/thread_pool.hh
```

### 2.2 Public Headers (installed API)
- Location: `include/cosmosim/<area>/<name>.hpp`
- File name: `lower_snake.hpp`, class names inside may be `PascalCase`

**Example**
```
include/cosmosim/gravity/tree_pm.hpp
include/cosmosim/io/snapshot_reader.hpp
```

### 2.3 Python Bindings & Tools
- Modules: `lower_snake.py`
- Package layout: `python/pycosmosim/<area>/<name>.py`

**Examples**
```
python/pycosmosim/io/snapshot.py
python/pycosmosim/analysis/power_spectrum.py
```

---

## 3) Configurations & Presets
- **YAML profiles:** `configs/<profile>/<stem>.yaml`  
  Profiles: `cosmo_cube`, `zoom_in`, `isolated_cluster`, `isolated_galaxy`

**Grammar**
```
<stem> := <scenario>_<size>_v<major><minor>
<scenario> := mw_zoom | cluster_beta | cube_dm_only | disk_isolated
<size> := n256 | n512 | n1024         # indicative resolution/grid
```

**Examples**
```
configs/zoom_in/mw_zoom_n512_v01.yaml
configs/cosmo_cube/cube_dm_only_n256_v02.yaml
configs/isolated_galaxy/disk_isolated_n256_v01.yaml
```

- **CMake presets/toolchains:** `cmake/<platform>_<compiler>.cmake`  
  e.g., `cmake/linux_gcc.cmake`, `cmake/windows_msvc.cmake`

- **CI scripts:** `scripts/ci/<task>.sh` → `scripts/ci/run_tests.sh`

---

## 4) Runs, Outputs, and Artifacts
All run artifacts live under `output/<run_id>/` with **symlink** `output/latest -> output/<run_id>`.

### 4.1 Run Identifier
**Grammar**
```
<run_id> := <date8>_<tag>_<git7>
<date8>  := YYYYMMDD
<tag>    := short_kebab descriptor (e.g., "mw-zoom-n512")
<git7>   := short SHA (7 chars)
```
**Example:** `20251106_mw-zoom-n512_a1b2c3d`

### 4.2 Snapshots & Restarts
- **Snapshots (HDF5):** `snap_<index3>.hdf5` (zero-padded, 000…999)  
- **Restarts (HDF5):** `restart_<index3>.hdf5`  
- **Lightcones:** `lightcone_<zfixed>.hdf5` → `lightcone_z0.50.hdf5`

**Examples**
```
output/20251106_mw-zoom-n512_a1b2c3d/snap_000.hdf5
output/20251106_mw-zoom-n512_a1b2c3d/snap_120.hdf5
output/20251106_mw-zoom-n512_a1b2c3d/restart_060.hdf5
```

### 4.3 Logs, Metrics, and Profiles
- **Main log:** `run.log`  
- **JSON profile:** `profile_timeline.json` (timings per phase)  
- **Counters:** `counters.json` (kernel-level metrics)  
- **Config copy:** `config.yaml` (exact config used)  
- **Provenance:** `provenance.json` (compiler flags, CPU/GPU, git SHA, seeds)

### 4.4 In‑situ & Analysis Outputs
- **Power spectrum:** `pk_z<zval>.csv` → `pk_z0.000.csv`
- **2PCF:** `xi_r_z<zval>.csv`
- **Halo catalogs:** `halos_fof_z<zval>.hdf5`, `halos_so_z<zval>.hdf5`
- **Quicklooks (PNG):** `slice_<field>_z<zval>_<axis>.png` → `slice_rho_z0.000_xy.png`
- **Projections (PNG):** `proj_<field>_z<zval>_<axis>.png`
- **SFR history:** `sfr_history.csv`
- **Rotation curve:** `vrot_<radius_spec>_z<zval>.csv` → `vrot_2kpc_z0.000.csv`

### 4.5 Temporary/Checkpoint Conventions
- **Checkpoints:** `chkpt_<phase>_<index3>.bin` (internal binary format)
- **Scratch:** keep in `output/<run_id>/_scratch/` (leading underscore marks nonversioned)
- **Partial files:** suffix `.part` during write; finalize atomically to target name

---

## 5) Tests & Benchmarks
### 5.1 Unit Tests
- **File:** `tests/unit/<module>/test_<feature>_<case>.cc`
- **Golden data:** `tests/golden/<suite>/<case>/*`

**Examples**
```
tests/unit/gravity/test_pm_poisson_solution.cc
tests/unit/hydro/test_sod_1d.cc
tests/golden/hydro/sod_1d/expected_solution.csv
```

### 5.2 Integration Tests
- **File:** `tests/integration/<scenario>/it_<scenario>_<size>.cc`  
  e.g., `tests/integration/cosmo/it_cube_dm_only_n128.cc`

### 5.3 Benchmarks
- **File:** `bench/<module>/bench_<metric>.cc` → `bench/gravity/bench_tree_build_mpps.cc`

---

## 6) Documentation & Figures
- **Markdown pages:** `docs/<topic>.md` → `docs/build_instructions.md`
- **Design proposals (RFCs):** `docs/rfc/<rfc_id>_<slug>.md` → `docs/rfc/0001_time_integrator.md`
- **Figures/images:** `docs/img/<topic>-<short>.png` (kebab-case)  
  e.g., `docs/img/tree-layout-overview.png`

- **Release notes:** `docs/releases/v<semver>.md` → `docs/releases/v0.1.0.md`

---

## 7) Script & Tooling Conventions
- **Shell:** `scripts/<verb>_<noun>.sh` → `scripts/build_release.sh`
- **Python tools:** `tools/<verb>_<noun>.py` → `tools/convert_music_ic.py`
- **CMake Presets:** `CMakePresets.json` (standard name)

- **Environment example files:** `.env.example`, never `.env` committed.

---

## 8) HDF5 Group & Dataset Naming (recap)
- **Groups:** `/Header`, `/PartType0`, `/PartType1`, `/PartType4`, `/PartType5`
- **Datasets:** `Coordinates`, `Velocities`, `Masses`, `InternalEnergy`, `Density`, `Metallicity`, `SmoothingLength`, `ParticleIDs`
- **Header attrs:** `NumPart_ThisFile`, `NumPart_Total`, `BoxSize`, `Time`, `Redshift`, `Omega0`, `OmegaLambda`, `HubbleParam`, plus provenance keys.

> These dataset names are **canonical** for compatibility (yt/GADGET-style tooling). Do not rename.

---

## 9) File Stem Grammars (Regex)
Useful for validation scripts / CI.

- **Run ID:**  
  `^[0-9]{8}_[a-z0-9]+(?:-[a-z0-9]+)*_[0-9a-f]{7}$`

- **Snapshot/Restart:**  
  `^snap_[0-9]{3}\.hdf5$`  
  `^restart_[0-9]{3}\.hdf5$`

- **Config:**  
  `^[a-z0-9]+(?:_[a-z0-9]+)*_n[0-9]{3,4}_v[0-9]{2}\.ya?ml$`

- **PNG quicklooks:**  
  `^(slice|proj)_[a-z0-9]+_z[0-9]+\.[0-9]{3}_(xy|xz|yz)\.png$`

- **CSV metrics:**  
  `^(pk|xi_r|vrot|sfr_history)(_[a-z0-9]+)*(_z[0-9]+\.[0-9]{3})?\.csv$`

---

## 10) Reserved Names & OS Caveats
- Avoid Windows-reserved names: `CON, PRN, AUX, NUL, COM1..COM9, LPT1..LPT9`.
- Avoid trailing dots/spaces and case-only disambiguation.
- Keep path length under **240 chars** for Windows compatibility.
- Use forward slashes in code; C++ joins with portable utilities.

---

## 11) Symlinks & “Latest”
- Maintain `output/latest` → last completed `<run_id>`.
- Never point `latest` at a **running** run-id; update only on clean shutdown.

---

## 12) Governance
- Changes to these conventions require a PR labeled `file-naming-rules`.
- Add CI checks that validate names against the regex in §9.
- On approved changes, bump `FileNamingRulesVersion` in `provenance.json` within run folders.

---

## 13) Quick Reference (Cheat Sheet)
- **Everything lower_snake** (code, configs, scripts).  
- **Run folders:** `YYYYMMDD_tag_git7` (e.g., `20251106_mw-zoom-n512_a1b2c3d`).  
- **Snapshots:** `snap_###.hdf5` (zero-padded).  
- **Configs:** `<scenario>_<size>_vNN.yaml`.  
- **Logs/metrics:** `run.log`, `profile_timeline.json`, `pk_z0.000.csv`.  
- **Docs images:** `docs/img/kebab-case.png` only in docs.

---

**End of document — CosmoSim File Naming v1.0**
