# CosmoSim — **Master Variable Naming Conventions**
*Version 1.0 (2025‑11‑06)*

This document defines **authoritative naming rules** for the CosmoSim codebase and adjacent tools (C++ kernel, YAML configs, Python API, HDF5 I/O, analysis scripts). It is designed for **clarity, reproducibility, and performance‑aware code** that scales from desktop to HPC.

---

## 0) General Principles
- **Consistency beats preference.** If a rule exists here, follow it everywhere.
- **Semantic, readable names.** Avoid cryptic abbreviations (except the approved list below).
- **Unit clarity is mandatory.** Use suffixes (`_si`, `_cgs`, `_code`, `_comov`, `_phys`) or structured types to avoid unit ambiguity.
- **Coordinate frame clarity.** Encode frame/time in names where ambiguity is possible (e.g., `pos_comov`, `vel_phys`, `a_scale`).
- **Scope & lifetime signals.** Prefixes for globals, statics, members, device buffers, etc. (see §2.6).
- **Match external schemas.** When interoperating (e.g., GADGET HDF5), keep **their** canonical names, even if style differs (see §4).

**Approved Abbreviations** (use lower case in variables and namespaces; PascalCase in types):  
`dm, gas, star, bh, pm, fft, fftw, ic, io, id, rng, simd, aabb, bbox, mhd, rt, sfr, imf, fof, so, amr, mpc, kpc, pc, Myr, Gyr, H0`

---

## 1) C++ Core (Engine & Modules)
### 1.1 Naming Styles (by category)
- **Namespaces:** `lower_snake` → `cosmo`, `cosmo::io`, `cosmo::hydro`
- **Files:** `lower_snake.cc/.hh` → `tree_builder.cc`, `pm_solver.hh`
- **Classes / Structs / EnumTypes:** `PascalCase` → `TreeNode`, `ParticleBlock`, `HaloCatalog`
- **Functions / Methods:** `camelCase` → `buildTree()`, `computeForces()`, `advanceKDK()`
- **Variables (local / parameters):** `lower_snake` → `cell_count`, `theta_open`, `time_step`
- **Data members:** `m_lower_snake` (prefix `m_`) → `m_mass`, `m_child_index`
- **Constants / constexpr / compile‑time params:** `UPPER_SNAKE` → `MAX_LEVELS`, `DEFAULT_THETA`
- **Macros (rare):** `UPPER_SNAKE` and `COSMOSIM_` prefix → `COSMOSIM_LIKELY`
- **Template parameters:** `T`, `TScalar`, `TIndex`, `TPolicy` (PascalCase with `T` prefix)
- **Enum values (scoped `enum class`):** `PascalCase` → `enum class BcType { Periodic, Open, Reflective };`

### 1.2 Units & Frames
- **Units suffix (if using raw scalars):**  
  - `_si` (SI), `_cgs` (CGS), `_code` (normalized code units)  
  - `_phys` (proper/physical), `_comov` (comoving)  
  - Examples: `density_si`, `length_code`, `pos_comov`, `vel_phys`
- **Time/redshift/scale:** `t_phys`, `t_code`, `z`, `a_scale` (scale factor), `dt_code`
- Prefer **strong types** for units where possible (e.g., wrappers), else enforce suffix rules.

### 1.3 Spatial & Index Conventions
- **Indices:** `i, j, k` for loops; `pidx` (particle index), `cidx` (cell index), `lidx` (level index)
- **Components:** suffix `_x/_y/_z` for SoA arrays → `pos_x`, `pos_y`, `pos_z`
- **Ranges:** `i_begin`, `i_end` (half‑open). Iteration blocks: `start_idx`, `end_idx`.

### 1.4 Physics Objects & Data Layout
- **SoA arrays:** `pos_x`, `pos_y`, `pos_z`, `vel_x`, `vel_y`, `vel_z`, `mass`, `id`
- **Gas attributes:** `rho`, `u_int` (specific internal energy), `press`, `h_smooth`, `Z_metal`
- **Star/BH:** `birth_time`, `sfr_inst`, `mdot_bh`, `e_feedback`
- **Flags/Tags:** predicate style → `is_active`, `is_leaf`, `needs_refine`

### 1.5 Concurrency & GPU
- **Thread‑locals:** `tls_` prefix → `tls_rng`, `tls_scratch`
- **Atomics:** `atomic_` prefix → `atomic_task_counter`
- **Mutexes / Locks:** `mtx_...`, `lock_...` → `mtx_tree`, `lock_queue`
- **OpenMP:** `omp_` prefix where needed → `omp_chunk_size`
- **Device / Host buffers:** `d_` for device, `h_` for host, `p_` for pinned → `d_pos_x`, `h_pos_x`, `p_sendbuf`
- **CUDA kernels:** `camelCase` with `Kernel` suffix → `computeForcesKernel<<<>>>`
- **SIMD lanes:** avoid embedding ISA names; use traits types (`SimdVec<T>`) and keep variables generic.

### 1.6 Error Handling / Logging
- **Error codes:** `ErrCode::IoFailure`, `ErrCode::InvalidParam`
- **Status booleans:** `ok`, `has_error`
- **Logger categories:** `log_io`, `log_tree`, `log_hydro` (lower_snake); messages begin with verb.

### 1.7 Configuration Objects
- **Config holders:** `SimulationConfig`, `GravityParams`, `HydroParams`
- Members use `m_` and `lower_snake` → `m_theta_open`, `m_pm_grid`, `m_softening_dm`

---

## 2) YAML / JSON Configuration
- **Top‑level keys:** `lower_snake` → `run`, `cosmology`, `physics`, `numerics`, `parallel`
- **Profiles / modes:** `cosmo_cube`, `zoom_in`, `isolated_cluster`, `isolated_galaxy`
- **Booleans:** `true/false` (never `yes/no`)
- **Enums (string values):** `lower_snake` → `gravity: treepm`, `time_integrator: kdk`
- **Durations & units:** attach space with unit → `snapshot_interval: 50 Myr`
- **Vectors:** YAML lists → `softening: { dm: 0.8 kpc, gas: adaptive }`
- **Validation:** missing/invalid keys raise `ConfigError` with path, e.g., `physics.cooling`

**Reserved keys (do not rename):** `H0`, `Omega_m`, `Omega_b`, `Omega_L`, `sigma8`, `ns` (match literature).

---

## 3) Python (pycosmosim Bindings & Notebooks)
- **Modules / packages:** `lower_snake` → `pycosmosim.io`
- **Classes:** `PascalCase` → `Snapshot`, `HaloFinder`
- **Functions:** `lower_snake` → `load_snapshot()`, `power_spectrum()`
- **Variables:** `lower_snake` → `box_size`, `redshift`
- **Constants:** `UPPER_SNAKE` → `G_NEWTON`, `SOLAR_MASS`
- **Jupyter artifacts:** temp vars start with `_` → `_fig`, `_ax`
- **Avoid shadowing C++ names.** Provide thin, pythonic façades over core C++ types.

---

## 4) HDF5 Snapshot Schema
To ensure compatibility with GADGET‑style tools and yt, **use canonical dataset names exactly** (PascalCase & pluralization as defined).

- **Groups:** `/Header`, `/PartType0`, `/PartType1`, `/PartType4`, `/PartType5`
- **Datasets (examples):** `Coordinates`, `Velocities`, `Masses`, `InternalEnergy`, `Density`, `Metallicity`, `SmoothingLength`, `ParticleIDs`
- **Header attributes:** `NumPart_ThisFile`, `NumPart_Total`, `BoxSize`, `Time`, `Redshift`, `Omega0`, `OmegaLambda`, `HubbleParam`, plus provenance (`GitCommit`, `ConfigHash`).
- **Compression & chunking keys** are implementation details; keep dataset names stable across versions.

> **Note:** Even though this differs from our C++ naming style, **do not change these HDF5 names**.

---

## 5) Analysis & In‑situ Outputs
- **File stems:** `lower_snake` → `pk_z1.0.csv`, `halo_catalog_z0.5.hdf5`
- **Columns / fields:** `lower_snake` with unit suffixes where needed → `k_hmpc`, `pk_mpk3`
- **Figures:** file names include metric and redshift → `power_spectrum_z0.0.png`

---

## 6) Test Fixtures & Benchmarks
- **Tests:** `test_<feature>_<case>` → `test_pm_poisson_solution`, `test_tree_force_error`
- **Golden data files:** `golden/<suite>/<case>/*`
- **Benchmarks:** `bench_<module>_<metric>` → `bench_tree_build_mpps`

---

## 7) Prefix & Suffix Quick Reference
- `m_`   → class/struct data members (`m_mass`, `m_child_idx`)
- `d_`   → device (GPU) buffer (`d_pos_x`), `h_` host, `p_` pinned
- `tls_` → thread‑local storage (`tls_rng`)
- `mtx_` / `lock_` → synchronization primitives
- `_si` / `_cgs` / `_code` → units
- `_phys` / `_comov` → frame
- `_x/_y/_z` → Cartesian components
- `_begin/_end` → half‑open ranges
- `is_` / `has_` / `can_` → boolean predicates

---

## 8) Examples (Canonical)
```cpp
namespace cosmo::gravity {
  struct TreeNode {
    float m_mass;
    float m_com_x, m_com_y, m_com_z;   // center of mass (code units, comoving)
    int   m_child_index[8];
    bool  m_is_leaf;
  };

  // Build & traverse
  void buildTree(const ParticleBlock& block, float theta_open);
  void computeForces(const ParticleBlock& block,
                     float a_scale, float softening_code,
                     float* d_acc_x, float* d_acc_y, float* d_acc_z);
}
```

```yaml
run:
  mode: zoom_in
  end_time: z=0
cosmology:
  H0: 67.7
  Omega_m: 0.31
  Omega_b: 0.048
  Omega_L: 0.69
physics:
  hydro: true
  star_formation: sk_law
numerics:
  gravity: treepm
  theta: 0.6
  pm_grid: 512
parallel:
  openmp_threads: auto
```
```python
# pycosmosim example
from pycosmosim import io, analysis

snap = io.load_snapshot("snap_050.hdf5")
k, pk = analysis.power_spectrum(snap.positions_comov, box_size_mpc)
```
---

## 9) Anti‑Patterns (Do Not Do)
- Mixing frames without suffixes: `pos` vs `pos_phys` — **always tag**.
- Ambiguous one‑letter vars in public APIs (`x`, `y`) — use in local scopes only.
- Reinventing HDF5 dataset names — keep canonical GADGET names.
- Encoding SIMD/ISA in variable names (`acc_x_avx2`) — use type traits instead.

---

## 10) Changelog & Governance
- Changes to this document require a PR labeled `naming-rules` and approval by maintainers.
- Upon merge, bump the `NamingRulesVersion` field in `/Header` attributes of snapshots.

---

**End of document — CosmoSim Naming v1.0**