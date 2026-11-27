# CHUÍ

2026-10-12 TreePM recovery repairs the actual drift resource contract and adds
bounded periodic arithmetic, early descent, optional homogeneous residual/SIMD
kernels, production block/leaf tuning and explicit source/toolchain provenance.
**Compilation, tests, MPI, benchmarks and convergence were not run, by request.**
New arithmetic modes remain opt-in; global KDK and conservative numerical
defaults remain. The [source handoff](docs/repair/treepm_deep_hotpath_attempt1_closure_20261012.md)
contains the exact changes, memory contracts and manual qualification sequence.
Generate separate normalized experiment decks with
`cosmosim_treepm_experiment_configs BASE.param.txt NEW_DIRECTORY` after building.
The evolved snapshot example is conditional external-IC replay, not checkpoint
continuation. If extraction loses launcher permissions, use
`python3 ./chui run CONFIG --preset PRESET`.

P6 optional performance/integration architecture is **SOURCE-IMPLEMENTED /
VALIDATION PENDING**. Conservative defaults remain in effect. Read the
[implementation report](docs/repair/p6_preproduction_performance_architecture_implementation_20261006.md)
before enabling new modes; the [qualification campaign](docs/repair/p6_validation_campaign_plan_20261006.md)
has not been executed.

CHUÍ is a desktop-first and small-cluster-first cosmological simulation framework written in C++20. 

The project is designed to maximize scientifically useful cosmological simulation capability under explicit memory and reproducibility constraints while preserving a path to MPI/HPC execution. The long-term science target is galaxy formation from cosmological initial conditions; the current operational milestone is the dark-matter-only (DMO) TreePM path.

## Current milestone: real DMO first light

CHUÍ has completed a real short **32³ periodic cosmological dark-matter-only first-light smoke run** through the normal config-driven `cosmosim_harness` path using an unmodified monofonIC/GADGET-style PartType1 HDF5 initial-condition file. The exercised path covered external IC ingestion, runtime composition, TreePM gravity, adaptive rung-zero cosmological stepping, scale-factor endpoint integration, HDF5 science snapshots, final artifact flush, and clean process exit.

That result establishes an **operational DMO first-light path**, not scientific production certification. Rank-equivalence, long-interval growth, power-spectrum convergence, larger-target runtime/memory qualification, and full-physics validation remain separate gates.

## What works today

### Supported and exercised

- C++20 core with typed `.param.txt` configuration, normalization, provenance, and deterministic config hashing.
- Periodic cosmological DMO execution with TreePM gravity.
- Adaptive **global rung-zero** timestep selection and scale-factor endpoint authority.
- GADGET/AREPO-style HDF5 IC bridging, including the monofonIC structural forms exercised by first light.
- Periodic coordinate canonicalization and checked external-ID normalization.
- HDF5 science snapshots with transactional completion markers and readback validation.
- Same-world-size restart continuation where the documented restart path is enabled.
- Memory-governed runtime architecture with explicit ownership/reservation contracts for major live sets.
- MPI/FFTW code paths when the required development dependencies are available.

### Provisional or intentionally fail-closed

- Collisionless hierarchical KDK is available for opt-in qualification with `hierarchical_max_rung > 0`; it remains provisional. The safe production default is rung zero.
- Distributed IC import and multi-rank DMO execution remain provisional until the dependency-complete rank matrix and rank-equivalence gates are closed.
- Rank-count-changing restart is unsupported.
- Power-spectrum output is computationally available on supported FFT backends but is not yet a scientific production-certification claim.
- Full-physics hydro/AMR/subgrid modules have their own validation status and are not implied to be production-certified by DMO first light.

### Future acceptance gates

The next DMO qualification stages are an external-IC clean-tree replay with the repaired normalized artifact, multi-rank equivalence, linear-growth validation, power-spectrum validation/convergence, and larger-target memory/runtime qualification before any 512³ production claim.

The [20261010 TreePM recovery handoff](docs/repair/treepm_unified_performance_recovery_20261010.md) documents the deterministic traversal update, stable analytic Gaussian residual, experimental direct-pair lookup, frozen-source mesh/leaf/block sweep, and hierarchical qualification tests. These source changes have not been compiled, tested, or benchmarked in that campaign; lookup stays disabled by default.

## Quickstart

### CPU-only development build

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

For HDF5 I/O work:

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure
```

For the first-light-style TreePM path, use the repository's HDF5+FFTW preset:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
```

The shipped smoke deck is:

```text
configs/chui_firstlight_32_smoke.param.txt
```

It references:

```text
../ics/chui_ic_32_10mpch_z24.hdf5
```

IC paths are resolved relative to the config file, so with the repository layout above the external file should be placed at:

```text
ics/chui_ic_32_10mpch_z24.hdf5
```

The external monofonIC file is **not** bundled with the repository. After supplying a compatible IC, the normal human-facing entry point is:

```bash
./chui run configs/chui_firstlight_32_smoke.param.txt
```

If more than one runnable preset build exists, the launcher deliberately refuses to guess. Select the intended build explicitly. Unless `--quiet` is used, `[CHUI][LAUNCH]` records the selected executable, build directory, launch mode, and requested MPI rank count before process replacement:

```bash
./chui run configs/chui_firstlight_32_smoke.param.txt \
  --preset pm-hdf5-fftw-debug
```

The runtime now emits compact rank-0 native status records such as `[CHUI][START]`, `[CHUI][IC]`, `[CHUI][STEP]`, `[CHUI][SNAPSHOT]`, `[CHUI][RESTART]`, and `[CHUI][DONE]`. In MPI runs, startup/completion identify the logical shared `run_directory` separately from rank 0's `rank_directory`, and snapshot success reports the committed logical set (`members` + `.complete` manifest) rather than presenting `.0.hdf5` as the whole snapshot. Console cadence is presentation-only and can be controlled without changing the scientific configuration:

```bash
./chui run configs/chui_firstlight_32_smoke.param.txt --status-every 5 --status-seconds 20
./chui run configs/chui_firstlight_32_smoke.param.txt --quiet
```

For distributed qualification, use the documented `mpi-hdf5-fftw-*` presets only on a machine with MPI, **Parallel HDF5**, FFTW, and FFTW-MPI development support. For known CMake builds, `./chui run CONFIG` asks the C++ harness preflight for the authoritative typed `parallel.mpi_ranks_expected` value and composes MPI automatically when it is greater than one; Python does not parse or rewrite the scientific parameter file. An explicit `--mpi N` remains available and must agree with the typed config. The launcher rejects a selected build when its CMake cache explicitly records `COSMOSIM_ENABLE_MPI=OFF`:

```bash
./chui run configs/production.param.txt --preset mpi-hdf5-fftw-release --mpi 8
```

Direct harness execution remains supported for CI, debugging, and low-level use:

```bash
./build/pm-hdf5-fftw-debug/cosmosim_harness \
  configs/chui_firstlight_32_smoke.param.txt
```

`chui_telemetry.py` remains an optional process/resource diagnostic tool for suspicious runs; it is not part of the normal launch path. See [`docs/build_instructions.md`](docs/build_instructions.md).

## Output layout

New science snapshots are published into one flat snapshot directory. A committed serial run resembles:

```text
outputs/<run_name>/
├── normalized_config.param.txt
├── provenance...
├── runtime...
└── snapshots/
    ├── snap_001.hdf5
    ├── snap_001.complete
    ├── snap_002.hdf5
    ├── snap_002.complete
    └── ...
```

For current small/moderate-rank production MPI runs, `output.snapshot_layout = auto` selects the same analysis-ready **single-file** science product when the build has MPI-enabled Parallel HDF5:

```text
snap_042.hdf5
snap_042.complete
```

All ranks write their owned rows directly into global HDF5 datasets; the implementation does not gather the universe onto rank 0. The final `.hdf5` name is published atomically only after distributed exact partition readback succeeds. The `.complete` sidecar remains transactional CHUÍ evidence, but ordinary analysis opens `snap_042.hdf5` directly. `snapshot_layout = sharded` retains the legacy one-member-per-rank compatibility backend explicitly; `aggregated` is a reserved future fixed-file-count topology and currently fails closed rather than pretending to be implemented. Restart/checkpoint topology remains independent and rank-oriented. See [`docs/output_schema.md`](docs/output_schema.md) and [`docs/snapshot_hdf5_io.md`](docs/snapshot_hdf5_io.md).

## Architecture at a glance

```text
Typed config + provenance
          │
          v
Cosmology / TimeCoordinator
          │
          v
      TreePM gravity
          │
          v
 canonical SimulationState
          │
          v
 HDF5 snapshot / restart
          │
          v
 analysis + validation
```

The implementation is SoA-oriented and separates canonical scientific state from caches, phase scratch, communication buffers, and diagnostics. Runtime ownership follows a single-authority model; periodic/open geometry is an explicit mode-policy contract; large allocations are governed rather than treated as unbounded workspace. See [`docs/architecture/overview.md`](docs/architecture/overview.md), [`docs/architecture/runtime_truth_map.md`](docs/architecture/runtime_truth_map.md), and [`docs/architecture/adr_runtime_truth_ownership.md`](docs/architecture/adr_runtime_truth_ownership.md).

## Validation status

| Gate | Current status |
|---|---|
| Operational 32³ single-rank DMO first light | **Achieved** |
| Focused configuration and HDF5 IC regressions | **Passing on exercised CPU/HDF5 paths** |
| Normalized-config floating-point roundtrip | **Passing focused regression** |
| Normalized config reusable for HDF5 IC compatibility checks | **Passing focused regression** |
| Optional HDF5 `NumPart_Total_HighWord` absent/present-malformed handling | **Passing focused regression** |
| MPI distributed IC / rank equivalence | **Pending dependency-complete execution** |
| Linear growth validation | **Pending** |
| Power-spectrum scientific validation/convergence | **Pending** |
| 512³ production qualification | **Not certified** |
| Full-physics production certification | **Not certified** |
| Extended physics (neutrinos, observables, instruments, CR kernels; pass 01) | **Source only, not built or run** — see [ledger](docs/architecture/chui_extended_physics_implementation_pass01_status.md) |

The authoritative current status, including environment-specific validation evidence and blockers, is [`CURRENT_STATUS.md`](CURRENT_STATUS.md).

## Reproducibility discipline

CHUÍ treats reproducibility as a first-order runtime contract:

- authoritative typed configuration is normalized before execution;
- finite floating-point values in normalized config use a round-trip-safe textual representation;
- normalized configuration is hashed and preserved in run artifacts;
- IC convention, unit, frame, velocity, and species policies are explicit rather than guessed from filenames;
- snapshot/restart schemas are versioned and transactional publication prevents partial artifacts from being mistaken for completed state.

See [`docs/configuration.md`](docs/configuration.md), [`docs/output_schema.md`](docs/output_schema.md), and [`docs/restart_checkpointing.md`](docs/restart_checkpointing.md).

## Documentation map

- **Current snapshot truth:** [`CURRENT_STATUS.md`](CURRENT_STATUS.md)
- **Build and dependencies:** [`docs/build_instructions.md`](docs/build_instructions.md)
- **Runtime configuration:** [`docs/configuration.md`](docs/configuration.md)
- **Architecture and ownership:** [`docs/architecture/overview.md`](docs/architecture/overview.md)
- **Distributed IC ingestion:** [`docs/architecture/distributed_ic_ingestion.md`](docs/architecture/distributed_ic_ingestion.md)
- **Snapshot/output schema:** [`docs/output_schema.md`](docs/output_schema.md)
- **Restart contract:** [`docs/restart_checkpointing.md`](docs/restart_checkpointing.md)
- **Validation ladder:** [`docs/validation_plan.md`](docs/validation_plan.md)
- **Profiling:** [`docs/profiling.md`](docs/profiling.md)
- **Contribution workflow:** [`CONTRIBUTING.md`](CONTRIBUTING.md)
- **Agent/repository contract:** [`AGENTS.md`](AGENTS.md)

Full-physics modules such as cooling, star formation, feedback, metals, AMR, and black-hole physics retain their dedicated documentation and validation boundaries; their presence in the source tree does not broaden the DMO first-light certification claim.
