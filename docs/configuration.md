# Configuration reference

CosmoSim uses a GADGET/AREPO-style `.param.txt` workflow and normalizes it into a typed, validated `SimulationConfig` before execution. `.param.txt` is the authoritative simulation-input format for runtime decks, examples, release profiles, and test fixtures. YAML and JSON are not accepted alternative simulation-input formats; JSON remains appropriate for generated provenance, manifests, counters, profiles, and capability reports. See `docs/architecture/adr_param_txt_configuration_format.md`.

Authoritative structures: `include/cosmosim/core/config.hpp` and `src/core/config.cpp`.

## CLI runtime entry point

The shipped executable is now a real config-driven runtime application:

```bash
cosmosim_harness <config.param.txt>
```

That path:

- loads the config via `loadFrozenConfigFromFile(...)`
- validates the typed runtime contract
- constructs the reference workflow
- runs the live scheduler-driven execution path
- writes runtime artifacts into the config-driven run directory

The concrete run directory is:

```text
<output.output_directory>/<output.run_name>
```

## Parsing model

- Input supports `key=value` and `key value` forms.
- Comments: `#`, `;`, `//`.
- Sections (for example `[physics]`) prefix unqualified keys (`cooling_model` becomes `physics.cooling_model`).
- Unknown keys fail by default; compatibility opt-out is explicit (`compatibility.allow_unknown_keys=true`).

## Reproducibility invariants

- Config is normalized to canonical text.
- Every finite double-valued normalized field uses a shortest round-trip-safe decimal representation; reparsing `normalized_config.param.txt` therefore preserves the frozen binary value rather than a display-rounded approximation.
- Canonical text is hashed (stable FNV-1a) for provenance.
- The normalized config snapshot is written to `normalized_config.param.txt` inside the run directory once config loading succeeds.
- Stable naming constraints apply to `output.output_stem`, `output.restart_stem`, and diagnostics/halo stems.
- Mode policy is validated before runtime (`mode.mode`, boundary selection, zoom requirements).

## External strings vs internal typed contract

- User-facing param files remain string-based (`gravity_solver=treepm`, `fb_mode=momentum`, etc.).
- During freeze, policy-like values are converted to enums in `core::SimulationConfig`:
  - `units.coordinate_frame` -> `core::CoordinateFrame`
  - `numerics.gravity_solver` / `numerics.hydro_solver` -> `core::GravitySolver` / `core::HydroSolver`
  - `mode.hydro_boundary` / `mode.gravity_boundary` -> `core::ModeHydroBoundary` / `core::ModeGravityBoundary`
  - `physics.fb_mode` / `physics.fb_variant` -> `core::FeedbackMode` / `core::FeedbackVariant`
  - `physics.uv_background_model` / `physics.self_shielding_model` / `physics.cooling_model` -> `core::UvBackgroundModel` / `core::SelfShieldingModel` / `core::CoolingModel`
- Free-form strings remain strings only where the value is path-like, label-like, or intentionally open-ended (`ic_file`, output stems, table paths, `reionization_model`).

## Central key registry

- `src/core/config.cpp` contains a single authoritative key registry that records canonical keys and defaults.
- The same source file also defines deprecated alias mappings.
- Unknown-key checks and alias translation both use this registry path during freeze.

## Key groups

## `schema_version`

- `schema_version` (int, required to be `1` in this build)

## `[mode]`

- `mode` (`cosmo_cube`, `zoom_in`, `isolated_galaxy`, `isolated_cluster`)
- `ic_file`
- `ic_convention` (`generated`, `chui_canonical_v1`, `gadget_arepo_bridge_v1`, `manifest_v1`)
- `ic_manifest_file` (required for `manifest_v1`)
- `ic_chunk_particle_count` (positive bounded HDF5 read chunk; default `65536`)
- `ic_staging_particle_count` (positive bounded MPI routing/validation staging limit; default `65536`)
- `ic_part_type2_policy`, `ic_part_type3_policy` (`reject`, `dark_matter`, `star`, `black_hole`, `tracer`)
- missing-field policy keys: `ic_gas_internal_energy_policy`,
  `ic_gas_density_policy`, `ic_star_formation_time_policy`,
  `ic_star_initial_mass_policy`, `ic_star_metallicity_policy`, and
  `ic_bh_mdot_policy`
- missing-field policy values: `reject`, `reconstruct`, `use_config_value`,
  `dialect_defined_default`, `preserve_unavailable`; this build rejects
  `reconstruct` and `preserve_unavailable` until their scientific storage paths
  exist
- configured-value companions: `ic_gas_internal_energy_value_code`,
  `ic_gas_density_value_code`, `ic_star_formation_time_value`,
  `ic_star_initial_mass_value_code`, `ic_star_metallicity_value`, and
  `ic_bh_mdot_value_code`
- `zoom_high_res_region` (bool)
- `zoom_region_file` (required when `zoom_high_res_region=true`)
- `zoom_long_range_strategy` (`disabled`, `global_coarse_plus_focused_highres_correction`)
- `zoom_region_center_x`, `zoom_region_center_y`, `zoom_region_center_z` (comoving Mpc)
- `zoom_region_radius` (comoving Mpc; required `>0` when `zoom_high_res_region=true`)
- `zoom_focused_pm_grid_nx`, `zoom_focused_pm_grid_ny`, `zoom_focused_pm_grid_nz` (required when focused correction is enabled)
- `zoom_contamination_radius` (comoving Mpc; defaults to `zoom_region_radius` when `<=0`)
- `hydro_boundary` (`auto`, `periodic`, `open`, `reflective`)
- `gravity_boundary` (`auto`, `periodic`, `isolated_monopole`)
  - `isolated_monopole` activates non-periodic/open PM long-range gravity.
    Single-rank execution is the normal path. Multi-rank execution is accepted
    only by the bounded root gather/solve/scatter compatibility route governed
    by `parallel.isolated_pm_root_workspace_limit_bytes`; it is not a scalable
    distributed-open-PM claim.

### Initial-condition convention contract

External initial conditions are never selected by filename alone. The typed
`mode.ic_convention` value owns the scientific interpretation:

- `generated`: requires `mode.ic_file=generated` and uses the in-repository generator.
- `chui_canonical_v1`: reads an IC already expressed in the canonical CHUÍ v1 units,
  frame, velocity, species, and provenance contract.
- `gadget_arepo_bridge_v1`: requires a complete direct-bridge scientific
  convention. The source length, mass, and velocity SI scales, coordinate frame,
  velocity convention, and all length/mass/velocity h and scale-factor exponents
  are mandatory. There are no kpc/Msun/km/s or frame defaults. Standard
  GADGET/AREPO cosmological snapshot velocities must declare
  `gadget_arepo_stored_peculiar`, which applies
  `v_peculiar = v_stored * sqrt(a)`. The legacy/custom
  `sqrt_a_scaled_peculiar` convention is intentionally distinct and retains its
  existing inverse-`sqrt(a)` import semantic.
- `manifest_v1`: loads the strict versioned audit manifest named by
  `mode.ic_manifest_file`; relative source paths are resolved from the manifest.
  `mode.ic_file` is not required. If it is supplied, it is only a compatibility
  consistency check and must resolve to the same first source member or runtime
  startup fails as a conflicting authority. Source order, hashes, field
  contracts, species policy, and missing-field policy/value pairs are owned by
  the manifest; direct-bridge policy keys are not duplicate scientific
  authorities in this mode.

An external `mode.ic_file` without an explicitly written `mode.ic_convention`
fails validation. This prevents the runtime from guessing h scaling, coordinate
frame, scale-factor powers, velocity convention, or species policy. The normalized
configuration records the resolved convention, manifest path, chunk/staging limits,
and Type2/Type3 policies.

Example direct bridge import:

```ini
[mode]
mode = cosmo_cube
ic_file = ../ics/cube.0.hdf5
ic_convention = gadget_arepo_bridge_v1
ic_chunk_particle_count = 65536
ic_staging_particle_count = 65536
ic_part_type2_policy = reject
ic_part_type3_policy = reject
ic_bridge_source_length_unit_to_si = 3.0856775814913673e22
ic_bridge_source_mass_unit_to_si = 1.98847e30
ic_bridge_source_velocity_unit_to_si = 1000
ic_bridge_coordinate_frame = comoving
ic_bridge_velocity_convention = gadget_arepo_stored_peculiar
ic_bridge_length_hubble_exponent = -1
ic_bridge_length_scale_factor_exponent = 0
ic_bridge_mass_hubble_exponent = -1
ic_bridge_mass_scale_factor_exponent = 0
ic_bridge_velocity_hubble_exponent = 0
ic_bridge_velocity_scale_factor_exponent = 0
```


Missing physical fields fail closed by default. Gas `InternalEnergy` and
`Density`, stellar formation time/initial mass/metallicity, and black-hole
`BH_Mdot` are imported only when present or when an explicit typed policy is
selected. Every non-reject resolution is recorded in normalized configuration,
the IC manifest, runtime provenance, and conversion diagnostics. Negative
stellar formation times are diagnosed as unsupported AREPO wind particles
rather than reinterpreted as ordinary stars.

Direct-bridge values are part of the normalized configuration and provenance hash.
Accepted coordinate frames are `comoving` and `physical`. Accepted velocity
conventions are `physical_peculiar`, `sqrt_a_scaled_peculiar`, and
`comoving_coordinate_rate`. Unit scales must be finite and positive; all six
h/a exponents must be finite. Omitting any member of the contract is a
configuration error before HDF5 is opened.

Example manifest-driven import:

```ini
[mode]
mode = cosmo_cube
ic_file = ../ics/cube.0.hdf5
ic_convention = manifest_v1
ic_manifest_file = ../ics/cube.ic_manifest.json
```

## `[cosmology]`

- `omega_matter`, `omega_lambda`, `omega_baryon`
- `hubble_param`, `sigma8`, `scalar_index_ns`

The production `.param.txt` cosmology is currently a flat, radiation-free
matter-plus-Lambda model. Validation therefore requires
`omega_matter + omega_lambda = 1` within `1e-10`. This makes the configured
`hubble_param` the actual present-day Hubble constant, i.e. `H(a=1) = H0`.
Einstein-de Sitter (`omega_matter=1`, `omega_lambda=0`) is accepted.
- canonical axis-aware box lengths:
  - `box_size_x`, `box_size_y`, `box_size_z` (recommended; normalized output always emits these keys)
- backward-compatible scalar input:
  - `box_size` / `box_size_mpc_comoving` (legacy alias; when provided alone it maps to `box_size_x=y=z`)

## `[numerics]`

- `t_code_begin`, `t_code_end`
- `a_begin`, `a_end`, `z_begin`, `z_end`, `t_phys_begin`, `t_phys_end`, `integrator_time_variable`
- `max_global_steps`, `hierarchical_max_rung`, `amr_max_level`
  - `hierarchical_max_rung=0` retains the default global KDK path.
    Values `1..12` explicitly select source-implemented, unqualified DMO block
    KDK in `cosmo_cube`. Gas, stars, black holes, tracers, and enabled source
    physics are rejected. Bins are fixed inside each coarse block; code-time
    output events and cosmological limits clip the coarse interval. See the
    P6 operator and restart contract in `docs/time_integration.md`.
- `gravity_softening` / `gravity_softening_kpc_comoving`
  - the production gravity timestep combines this comoving length with
    comoving-coordinate acceleration `|A|/a^3`, not directly with scale-free
    TreePM `|A|` or the peculiar-velocity response `|A|/a^2`; after step commit
    it evaluates `dt_grav=eta sqrt(a^3 epsilon_com/|A|)`.
- `gravity_solver`, `hydro_solver`
- `hydro_active_batch_max_cells` (uint64, default `0`)
  - `0` selects the deterministic automatic ceiling (`16384` cells per active batch)
    before runtime headroom is considered.
  - a positive value caps the active hydro reconstruction/Riemann batch directly.
  - the existing process `MemoryGovernor` may only reduce this cap from current
    headroom; it does not create a second hydro memory authority or adaptive AIMD
    controller. The selected batch is aligned to eight cells when headroom permits.
- Hydro numerical-method selectors (CHUI Part 1; equations and compatibility
  matrix in `docs/hydro_numerics.md`). Defaults reproduce the historical
  PLM/MC/MUSCL-Hancock/HLLC path apart from the documented F1/F3 bug fixes.
  - `hydro_reconstruction` (`piecewise_constant` | `piecewise_linear` |
    `piecewise_parabolic`, default `piecewise_linear`). PPM needs a 3-cell
    stencil per side; near non-periodic patch boundaries it reduces to PLM and
    the reduction is counted (`stencil_degradation_count`).
  - `hydro_slope_limiter` (`minmod` | `mc` | `van_leer`, default `mc`). Used by
    PLM and by the PPM interface slopes/flattening fallback.
  - `hydro_integration_scheme` (`muscl_hancock` | `ctu`, default
    `muscl_hancock`). `ctu` is the simplified 3D corner-transport-upwind variant
    (2 Riemann solves per face). It is stable for Courant <= 1/2 in 3D; the
    runtime hydro CFL coefficient is 0.4.
  - `hydro_riemann_solver` (`hllc` | `hlle`, default `hllc`). HLLC falls back to
    HLLE (Einfeldt bounds) for degenerate/inadmissible star states.
  - `hydro_dual_energy` (`off` | `thermal_fraction`, default `off`) and
    `hydro_dual_energy_eta` (double in (0, 1), default `1.0e-3`). With
    `thermal_fraction`, cells whose conservative thermal energy is below
    `eta` of the total energy take their thermal state from an auxiliary
    internal-energy lane initialized from the canonical `e`; the solver never
    overwrites `E`, and the resulting non-conservative thermal correction is
    ledgered in `dual_energy_sync_delta` (`docs/hydro_numerics.md` section 7).
  - Validation is fail-closed: unknown enum spellings and eta outside (0, 1)
    are rejected. All six keys appear in the normalized dump, so enabling them
    changes the normalized-config hash (reproducibility: intended).

Time/scale semantics (anti-ambiguity contract):

- `integrator_time_variable` is runtime authority, not descriptive metadata.
  `scale_factor` and `ln_a` make `a_begin/a_end` authoritative and derive the
  code-time interval from the configured flat FLRW background. `code_time`
  makes `t_code_begin/t_code_end` authoritative and derives the matching final
  scale factor. `physical_time` is rejected at config freeze until a complete
  production conversion path exists.
- A cosmological config that explicitly supplies contradictory scale-factor and
  code-time endpoint families fails before IC ingestion. Normalized config emits
  the resolved `a`, `z`, and code-time endpoints actually used by the runtime.
- If no time authority is explicitly named, a supplied `a_end`/`z_end` selects
  scale-factor authority; a code-time-only endpoint selects code time.
  Non-cosmological modes use code time.
- `max_global_steps` is a termination/safety cap. It does not generate the
  physical timestep. In the supported rung-zero path the accepted global step is
  the MPI minimum of active physical criteria, optional explicit test/embedding
  maximum, ordered output boundary, and authoritative run endpoint.
- Legacy `time_begin_code`/`time_end_code` and `initial_scale_factor`/`initial_redshift` are accepted only as user-input aliases.
- Committed/restart cosmological scale-factor authority is
  `IntegratorState.current_scale_factor`, while redshift remains derived
  (`z=1/a-1`). Within an in-flight step, `StepContext.timeline_step` owns stage
  epochs; post-drift hydro uses its `scale_factor_end/hubble_end_code` before
  `IntegratorState` commits.
- TreePM runtime controls (typed + normalized, no hidden workflow defaults):
  - canonical axis-aware PM grid:
    - `treepm_pm_grid_nx`, `treepm_pm_grid_ny`, `treepm_pm_grid_nz` (normalized output always emits these keys)
  - backward-compatible scalar input:
    - `treepm_pm_grid` (legacy alias; when provided alone it maps to `treepm_pm_grid_nx=ny=nz`)
  - `treepm_asmth_cells` (float, default `1.25`; split scale in mesh-cell units)
  - `treepm_rcut_cells` (finite positive float, default `6.25`; user-facing short-range cutoff control in mesh-cell units; this drives explicit residual-traversal pruning in the TreePM coupling path)
  - `treepm_tree_opening_criterion` (`geometric`, `com_distance`, or `relative_force_error`; default `com_distance`)
  - `treepm_tree_opening_theta` (finite positive float, default `0.7`; Barnes-Hut threshold and deterministic first-evaluation fallback threshold for the relative-force criterion)
  - `treepm_tree_relative_force_tolerance` (finite positive float, default `0.005`; relative-force MAC tolerance)
  - `treepm_tree_relative_force_acceleration_floor` (finite positive float in code-acceleration units, default `1e-30`; lower bound for the previous-acceleration scale used by the relative-force MAC)
  - `treepm_assignment_scheme` (`cic` or `tsc`; default `tsc`; the independent Ewald gate certifies TSC, while CIC remains available as a diagnostic/compatibility variant)
  - `treepm_enable_window_deconvolution` (bool, default `true`; applies
    scheme-aware PM transfer deconvolution for both `cic` and `tsc` on the
    periodic solver path). Isolated/open gravity requires this value to be
    explicitly `false`; typed config validation rejects `true` because that
    solver path has no periodic assignment window to deconvolve.
  - `treepm_update_cadence_steps` (int, default and currently required value
    `1`; every integrator-issued production force-refresh surface rebuilds the
    PM long-range field)
  - `treepm_allow_diagnostic_naive_dft` (bool, default `false`; explicit
    diagnostic/test override permitting the tiny-grid non-FFTW DFT backend.
    Production TreePM fails before evolution when no production FFT backend is
    available and this override is false)
  - `treepm_pm_decomposition_mode` (`slab` or `pencil`; default `slab`; `pencil` selects FFTW-MPI transposed spectral ownership while real-space particle deposition/interpolation retains x-slab ownership)
  - `treepm_tree_exchange_batch_bytes` (uint64 bytes, default `4194304`; cap for tree export/import payload chunking in distributed gravity paths)
  - `treepm_pm_exchange_batch_bytes` (uint64 bytes, default `16777216`; per-peer wire-payload high-water for distributed PM density/force routing; the implementation communicates in globally coordinated bounded rounds and does not retain TSC cell records proportional to `27 * N_local`)

TreePM split/cutoff semantics in this phase:

- `Δx = cosmology.box_size_x / numerics.treepm_pm_grid_nx`
- `Δy = cosmology.box_size_y / numerics.treepm_pm_grid_ny`
- `Δz = cosmology.box_size_z / numerics.treepm_pm_grid_nz`
- representative split spacing for TreePM scalar split controls remains
  `Δmesh = cbrt(Δx * Δy * Δz)` in this phase (pencil/isotropic split redesign deferred)
- `r_s = numerics.treepm_asmth_cells * Δmesh`
- `r_cut = numerics.treepm_rcut_cells * Δmesh`

For periodic gravity, typed validation additionally requires
`r_cut < 0.5 * min(Lx,Ly,Lz)`. The short-range residual evaluates one periodic
minimum image per source, so equality or a larger cutoff would make that
single-image traversal ambiguous. Invalid decks fail with the derived cutoff,
half-shortest-axis bound, and PM grid shape. Increase PM grid resolution or
reduce `treepm_rcut_cells`; the restriction does not apply to the isolated/open
operator.

Normalization emits the dimensionless controls exactly as provided and these are the same values consumed by runtime mapping.

Tree opening-criterion semantics are:

- `geometric`: accept an internal node when `l / r < theta`, where `l` is node width and `r` is target-to-node-COM distance.
- `com_distance` (default): tighten the Barnes-Hut test to account for the node center-to-COM offset, `(l + delta_com) / r < theta`.
- `relative_force_error`: accept when `G M l^2 <= alpha max(|a_previous|, a_floor) r^4`, using `treepm_tree_relative_force_tolerance` for `alpha`. When a finite previous-acceleration magnitude is unavailable, traversal deterministically falls back to the `com_distance` rule with `treepm_tree_opening_theta`.

All three criteria, the opening threshold, relative tolerance, and acceleration floor are typed fields in the frozen config and are emitted in normalized config text. This makes criterion selection and its numerical thresholds part of the reproducibility hash.

`treepm_update_cadence_steps` is production-restricted to `1`:

- Values other than `1` fail config validation with an explicit production-safety error.
- In default `hierarchical_max_rung=0`, every integrator-issued
  rank-coordinated force-refresh surface rebuilds PM. The optional hierarchy
  instead has explicit coarse PM endpoint surfaces; fine Tree-only solves
  perform no PM interpolation and are not PM cache-reuse calls.
- The coordinator's explicit lower-level reuse API requires a transient cache
  signature match for force epoch, force-evaluation scale factor, `G_code`,
  split scale, all three box axes, assignment scheme, boundary condition, PM
  decomposition mode, and window-deconvolution policy. A missing or
  incompatible cache makes reuse fail coherently; it is not silently converted
  into a solve. Divergent explicit refresh votes fail before PM collectives.
- Particle decomposition epoch is deliberately not part of PM compatibility:
  particles may change owner while the PM field remains owned by the same fixed
  FFT slabs. Tree packets still carry the actual workflow decomposition epoch,
  and dense-row force history is invalidated on an ownership change.
- That signature prevents accidental reuse across known solver/runtime changes,
  but it is not a long-range time predictor/interpolator for arbitrary
  multi-kick evolution.

`treepm_rcut_cells` is normalized, carried into provenance, and consumed by runtime residual traversal pruning (node-level AABB pruning + acceptance guard + pair culling beyond `r_cut`).

The certified default combines `treepm_asmth_cells=1.25` with
`treepm_rcut_cells=6.25`, so `r_cut/r_s=5`. The former `4.5`-cell cutoff and
other positive explicit values remain available for compatibility and
parameter sweeps when the periodic half-shortest-axis geometry guard is also
satisfied, but are not covered by the default Ewald certification: the
4.5-cell profile showed about nine-percent force error at the cutoff transition. Raising
the cutoff increases short-range traversal radius, pair work, and distributed
target/response volume; this is a deliberate accuracy/cost tradeoff and is part
of normalized config/provenance.
Zoom long-range strategy semantics:

- `disabled`: standard TreePM split force only.
- `global_coarse_plus_focused_highres_correction`:  
  `a_total = a_PM_global_coarse(all sources) + a_PM_zoom_correction(high-res targets) + a_tree_short_residual(all sources)`  
  where `a_PM_zoom_correction = a_PM_focused(high-res sources) - a_PM_coarse(high-res sources)` to avoid double counting while retaining global tidal content from the coarse PM solve.
  The focused correction uses the isolated/open convolution and forces window
  deconvolution off; isolated deconvolution is unsupported even when the global
  periodic PM profile enables TSC deconvolution.
`treepm_assignment_scheme` now maps directly to matched PM assignment+gather behavior:

- `cic`: 2-point/axis stencil, lower cost, stronger smoothing.
- `tsc`: 3-point/axis stencil, higher cost, smoother transfer and typically lower anisotropy.

`treepm_enable_window_deconvolution=true` deconvolves the matched deposit/gather
transfer window for the selected scheme with a safety floor in k-space and is
enabled by default for the certified periodic TSC profile. It is invalid with
an isolated/open gravity boundary. There is no silent mode-dependent override:
tracked isolated decks set the key to `false`, and user isolated decks must do
the same.
`treepm_pm_decomposition_mode`, `treepm_tree_exchange_batch_bytes`, and `treepm_pm_exchange_batch_bytes` are live
distributed runtime controls, not ignored Phase 2 placeholders. Both slab and
pencil FFTW-MPI paths are implemented; their current np1--np4 correctness
evidence remains a small-cluster validation boundary, not a large-scale
performance claim.

## Migration notes (axis-aware TreePM geometry)

- Existing scalar configs remain valid:
  - `cosmology.box_size = L` maps to `box_size_x=box_size_y=box_size_z=L`.
  - `numerics.treepm_pm_grid = N` maps to `treepm_pm_grid_nx=ny=nz=N`.
- Mixed partial axis input is rejected:
  - either provide all three axis keys or provide the scalar compatibility key.
- Normalized config dumps now always emit canonical axis-aware keys (`box_size_{x,y,z}`, `treepm_pm_grid_n{xyz}`).

## Migration note (isolated/open window deconvolution)

The global default for `numerics.treepm_enable_window_deconvolution` remains
`true` for the periodic TSC profile. Consequently, an isolated/open deck that
omits this key inherits `true` and now fails during typed config loading with a
key-specific error. Set
`numerics.treepm_enable_window_deconvolution = false` explicitly. Runtime code
does not silently rewrite the frozen config, so normalized dumps and
reproducibility hashes continue to describe the solver policy actually used.

## Migration note (periodic cutoff geometry)

Periodic decks must satisfy
`treepm_rcut_cells * cbrt((Lx/Nx)(Ly/Ny)(Lz/Nz)) < min(Lx,Ly,Lz)/2`.
Older coarse-mesh decks that were accepted at or above this boundary now fail
during typed loading, and direct API callers are rejected again at TreePM solve
entry. Increase `treepm_pm_grid_n{x,y,z}` as appropriate or explicitly lower
`treepm_rcut_cells`; either choice changes the normalized numerical profile and
therefore requires fresh force-accuracy evidence.

## `[physics]`

Core toggles:

- `enable_cooling`, `enable_star_formation`, `enable_feedback`, `enable_stellar_evolution`

Cooling/heating:

- `reionization_model`, `uv_background_model`, `self_shielding_model`
- `cooling_model`, `metal_line_table_path`, `temperature_floor_k`

Star formation:

- `sf_density_threshold_code`, `sf_temperature_threshold_k`
- `sf_min_converging_flow_rate_code`, `sf_epsilon_ff`
- `sf_min_star_particle_mass_code`, `sf_stochastic_spawning`, `sf_random_seed`

Stellar feedback:

- `fb_mode` (`thermal`, `kinetic`, `momentum`, `thermal_kinetic_momentum`)
- `fb_variant` (`none`, `delayed_cooling`, `stochastic`)
- `fb_use_returned_mass_budget`
- `fb_epsilon_thermal`, `fb_epsilon_kinetic`, `fb_epsilon_momentum`
- `fb_sn_energy_erg_per_mass_code`, `fb_momentum_code_per_mass_code`
- `fb_neighbor_count`, `fb_delayed_cooling_time_code`
- `fb_stochastic_event_probability`, `fb_random_seed`

Stellar evolution + AGN:

- `stellar_evolution_table_path`, `stellar_evolution_hubble_time_years`
- `enable_black_hole_agn`
- `bh_seed_halo_mass_threshold_code`, `bh_seed_mass_code`, `bh_seed_max_per_cell`
- `bh_alpha_bondi`, `bh_use_eddington_cap`, `bh_epsilon_r`, `bh_epsilon_f`
- `bh_feedback_coupling_efficiency`, `bh_duty_cycle_active_edd_ratio_threshold`
- `bh_proton_mass_si`, `bh_thomson_cross_section_si`, `bh_newton_g_si`, `bh_speed_of_light_si`

Tracers:

- `enable_tracers`, `tracer_track_mass`, `tracer_min_host_mass_code`

## `[output]`

- `run_name`
- `output_directory` (output root; the runtime appends `run_name` to form the concrete run directory)
- `output_stem`, `restart_stem` (stable character set only)
- `snapshot_interval_steps` (zero disables step-modulo events)
- `snapshot_interval_time_code` (zero disables code-time events; positive values are
  anchored at `numerics.time_begin_code` and steps are clipped rather than crossing them)
- `snapshot_layout` (`auto`, `single`, `aggregated`, or `sharded`; default `auto`)
  - serial output is one ordinary HDF5 science file;
  - MPI `auto`/`single` requires an MPI-enabled Parallel-HDF5 build and writes one collective `snap_NNN.hdf5` without a root gather;
  - `sharded` explicitly selects the legacy one-file-per-rank compatibility topology;
  - `aggregated` reserves a future fixed-file-count topology and currently fails closed at runtime.
- `snapshot_num_files` (default `0`, automatic; values greater than one are accepted only with `snapshot_layout = aggregated` and do not currently activate an unimplemented backend)
- `write_restarts`
- `restart_interval_steps` (default `0`: checkpoints follow snapshot events; `N > 0`: independent checkpoint cadence)
- `restart_interval_time_code` (default `0`; `T > 0`: independent code-time checkpoint cadence, combinable with `restart_interval_steps`)
- `snapshot_verification` (`full_readback` default, or `structural`)

At least one snapshot interval must be positive. Code-time cadence is part of the
normalized config and its next ordered event is restart-authoritative in restart
schema v21.

## `[analysis]`

- `enable_diagnostics`, `enable_halo_workflow`, `halo_on_the_fly`
- `diagnostics_execution_policy` (`run_health_only`, `run_health_and_light_science`, `all_including_provisional`)
- `run_health_interval_steps`, `science_light_interval_steps`, `science_heavy_interval_steps`
- `retention_bundle_count`
- `power_spectrum_mesh_n`, `power_spectrum_bin_count`, `sf_history_bin_count`, `quicklook_grid_n`
  - with the built-in FFT backend (`COSMOSIM_ENABLE_FFTW=OFF`), `power_spectrum_mesh_n` must be a power of two; the parser rejects unsupported meshes before execution;
  - FFTW-enabled builds accept FFTW-supported mesh sizes.
- `diagnostics_stem`, `halo_catalog_stem`, `merger_tree_stem`
- `halo_fof_linking_length_factor`, `halo_fof_min_group_size`
- `halo_include_gas`, `halo_include_stars`, `halo_include_black_holes`

Diagnostics maturity policy:

- `run_health_only`: only infrastructure health counters run.
- `run_health_and_light_science` (default): run-health + validated lightweight science diagnostics.
- `all_including_provisional`: also schedules scientifically provisional diagnostics. The power-spectrum estimator is computationally scalable FFT-backed; its direct DFT exists only as an explicit validation oracle.

## `[parallel]`

- `mpi_ranks_expected`, `omp_threads`, `gpu_devices`

CPU thread execution contract:

- `omp_threads = 0` selects the OpenMP runtime default when OpenMP is compiled in; positive values request an explicit team size. Builds without OpenMP reject values greater than one. Tree-gravity active-set traversal and production FFT analysis contain real OpenMP work when available.

GPU execution contract:

- `gpu_devices = 0` keeps PM execution on the host.
- `gpu_devices > 0` explicitly requests the CUDA PM path for CIC assignment/interpolation. The runtime validates that CUDA devices are visible and assigns ranks to devices round-robin across the requested device pool.
- `gpu_devices` must be `>= 0`; requesting more devices than are visible is a configuration/runtime error rather than a silent fallback.
- `deterministic_reduction`
- `decomposition_runtime_rebalance_enabled`
- `decomposition_debug_exact_ownership_audit` (default `false`): keep the routine runtime rebalance path on distributed compact SFC-cut metadata, but run the exact global ownership partition audit after a migration commit in debug/small MPI runs. This exact audit is read-only evidence and is not a production ownership authority.
- `isolated_pm_root_workspace_limit_bytes` (default `268435456`): maximum rank-0 transient workspace allowed for the current isolated/open PM root-gather path. Multi-rank isolated PM is accepted only inside this explicit small-grid budget and diagnostics must not describe it as scalable.
- `gravity_memory_budget_bytes` (default `0`, meaning unlimited): backward-compatible gravity-only preflight ceiling. CHUI-owned PM/FFT `PlanResources` arrays are now explicit known memory; only genuinely library-owned FFTW/CUDA plan/runtime allocations remain external.
- `process_memory_budget_bytes` (default `0`, meaning no total-process ceiling): authoritative composed per-rank process ceiling when set. The existing DMO preflight combines live canonical state, retained scheduler state, gravity/tree staging, PM fields, explicit FFT plan arrays, bounded PM/Tree communication, configured output/restart overlap, external reserves, and the process safety margin before large gravity staging allocations. M1B also constructs the process `MemoryGovernor` from this same ceiling for reserve-before-allocate runtime admission; zero remains unlimited while accounting/high-water stays active. This does not silently change the legacy gravity-only key and no second memory-budget configuration surface exists.
- Node-aware admission (CHUI Part 1; `docs/memory_governance.md`):
  - `node_memory_admission` (`off` | `shared_equal` | `shared_pool`, default
    `off`). `shared_pool` uses the same node limit/reserve but a node-shared
    accounting counter (ranks may use unequal amounts; the node-wide sum is
    bounded; admission under contention is not run-to-run reproducible). `off`
    preserves the historical per-process semantics exactly
    (`process_memory_budget_bytes = 0` stays unlimited). `shared_equal` derives
    the node limit as the minimum of `node_memory_budget_bytes` (if non-zero),
    the job cgroup limit, `SLURM_MEM_PER_NODE`, and physical RAM; subtracts
    `node_memory_reserve_bytes` once; and gives each co-located MPI rank
    `floor(usable / ranks_on_node)`. The per-rank governor ceiling is the
    smaller of that share and a non-zero `process_memory_budget_bytes`.
  - `node_memory_budget_bytes` (uint64, default `0` = discover) and
    `node_memory_reserve_bytes` (uint64, default `0`). Setting either while
    `node_memory_admission = off` is rejected (it would be silently ignored).
    A reserve that consumes the whole node limit is rejected.
- `process_mpi_unknown_reserve_bytes`, `process_fftw_unknown_reserve_bytes`, `process_hdf5_unknown_reserve_bytes`, `process_allocator_unknown_reserve_bytes` (defaults `0`): explicit per-rank reserves for library/allocator memory CHUI cannot derive from owned container geometry. Known CHUI arrays must not be hidden in these reserves.
- `process_output_restart_overlap_bytes` (default `0`): explicit owned I/O/restart staging allowance that may overlap the resident DMO+gravity state. IC import staging is excluded from the gravity-phase peak because importer ownership ends before the production integrator takes over the canonical state.
- `process_memory_safety_margin_fraction` (default `0`, range `[0,1]`): multiplicative margin applied after known owned bytes plus external reserves. The M1B governor applies the same multiplicative/ceiling rule to reconciled runtime demand so preflight and runtime admission use compatible policy arithmetic. Swap and host overcommit are not part of this budget.
- `zoom_high_res_allgather_limit_bytes` (default `268435456`): maximum transient four-field high-resolution source payload allowed for the current focused zoom PM correction all-gather. Exceeding it fails before allocation/communication.

Distributed-memory determinism contract (infrastructure scope):

- `mpi_ranks_expected` must match across all ranks for a valid distributed run contract.
- `deterministic_reduction=true` means rank-local contributions are compared against a rank-ordered deterministic reference sum in tests and smoke paths.
- `deterministic_reduction=false` is allowed only when reduction error is evaluated explicitly against a documented policy mode in code/tests.
- The normalized frozen config hash (`provenance.config_hash_hex`) is computed over the normalized config text itself and must match both the stored normalized-config hash field and the provenance hash field in continuation artifacts.

## `[units]`

- `length_unit`, `mass_unit`, `velocity_unit`, `coordinate_frame`

The three unit strings define one `core::UnitSystem` and its derived code-time
unit `T_code=L_code/V_code`. Production gravity derives

```text
G_code = G_SI * M_SI_per_code * T_SI_per_code^2 / L_SI_per_code^3
```

through `core::newtonGravitationalConstantCode(...)`; it does not assume
`G_code=1`. PM and the TreePM residual receive the same value. The normalized
unit strings and config hash remain reproducibility authority; the derived
value is also reported in the `gravity.treepm_setup` operational event. This is
an additive public helper, not a new config key or snapshot/restart schema
field.

## `[compatibility]`

- `allow_unknown_keys`

## Example configs

Canonical examples are in `configs/`:

- `minimal_cosmosim.param.txt`
- `cosmo_cube.param.txt`
- `zoom_in.param.txt`
- `isolated_galaxy.param.txt`
- `isolated_cluster.param.txt`
- `cooling_relaxation.param.txt`
- `configs/release/release_smoke_*.param.txt`


## Campaign B supplied-manifest authority

For `mode.ic_convention = manifest_v1`, `mode.ic_manifest_file` is the only
required source authority. `mode.ic_file` is optional and, when present, must
resolve to the same first source member. The validated manifest controls source
order, dialect, hashes, schema aliases, conversions, missing-field contracts and
replacement values, and species policy. Runtime configuration does not reapply
parallel scientific policies over the manifest.

The equivalent C++ APIs follow the same rule: any non-null supplied manifest is
validated or rejected; an incomplete manifest never falls back to source
discovery or live configuration. Distributed canonical manifests select
canonical inspection and therefore require a valid embedded/sidecar/marker
bundle.

## Campaign B capability and source-handoff status

No user configuration key or scientific ingestion policy changed in the final
source-package repair. Supplied-manifest authority and canonical bundle
validation remain as documented above.

Release packaging is now a build-validating source-distribution gate:

```bash
bash scripts/ci/package_source_zip.sh /absolute/path/to/source.zip archive_root
```

The command fails unless the actual source tree, staged tree, and extracted
archive are free of ADS/generated artifacts; required source sentinels exist;
staged and extracted file inventories and hashes match; and the extracted copy
configures with `cpu-only-debug` and builds `cosmosim_core` plus
`cosmosim_harness`.

This packaging closure does not promote distributed runtime capability. With
MPI/HDF5 enabled, `distributed_ic_import` remains `provisional` until the real
np1/np2/np4 and fault-injection acceptance matrix passes on a host with a
working MPI C++ toolchain.

## Adaptive bound–Jeans star formation

`physics.star_formation_model` is typed and accepts `legacy_schmidt_threshold` or `adaptive_bound_jeans`. Existing configurations without an explicit selector retain the legacy default; new reference profiles must select the adaptive model explicitly.

Adaptive keys:

- `sf_epsilon_ff`: local collapse efficiency; the reference profile uses `1.0` only after bound/Jeans/convergence filtering.
- `sf_bound_alpha_vir_max`: strict upper bound on the local support parameter.
- `sf_require_converging_flow`: require negative physical peculiar-flow divergence.
- `sf_collapse_timescale`: `free_fall` or `minimum_free_fall_or_compression`.
- `sf_jeans_mass_floor_code`: Jeans comparison floor.
- `sf_star_particle_mass_policy`: `fixed` or `gas_resolution_fraction`.
- `sf_target_star_particle_mass_code`, `sf_min_star_particle_mass_code`, `sf_max_star_particle_mass_code`.
- `sf_target_star_particle_mass_fraction`.
- `sf_max_spawn_particles_per_cell_step`.
- `sf_max_fractional_mass_conversion`.
- `sf_min_remaining_gas_fraction`, `sf_min_remaining_gas_mass_code`.
- `sf_random_seed`.
- `sf_temperature_safety_ceiling_k`; zero disables it.

Mass-valued keys accept raw code mass or explicit `kg`, `g`, and `msun` suffixes. Normalized configuration stores code-unit numbers and includes every key in the configuration hash. Unknown or misspelled keys fail parsing. `enable_star_formation=false` bypasses the production scan and consumes no star-formation random draw.


## Star-formation model selection

All star-forming profiles must explicitly set one of:

```text
physics.star_formation_model = legacy_schmidt_threshold
physics.star_formation_model = adaptive_bound_jeans
physics.star_formation_model = effective_multiphase_tng_like
```

The effective model provides typed keys for physical `sf_effective_n_h_threshold`, optional baryon overdensity, threshold timescale, evaporation law, cold/SN energies, massive-star fraction, `q_eos`, isothermal reference, hot-excess tolerance, relaxation policy, table resolution/range, SSP birth-mass convention, and feedback coupling. Number density accepts `cm^-3` or `m^-3`; temperatures accept `K`; times accept code, seconds, years, Myr, or Gyr; specific energy accepts code, J/kg, or erg/g. All resolved values enter normalized configuration and the config hash.

Known-inconsistent combinations, including `effective_eos_only` with explicit feedback enabled, fail validation. Shared particle-mass, stochastic, conversion, and remnant-gas keys are not duplicated for the effective model.

### Metals, delayed enrichment, and turbulent mixing

The scalar total-metal mode is the compatibility and memory-efficient default.

- `metal_species_mode`: `total_only` (default). `core_elements` is reserved and currently rejected fail-closed.
- `stellar_evolution_require_production_table`: default `false`; when true, a non-empty production-calibrated v2 table is required.
- `enable_metal_diffusion`: default `false`.
- `metal_diffusion_model`: `none` (default) or `smagorinsky`.
- `metal_diffusion_time_integrator`: `explicit_subcycling` (default) or `rkl2`.
- `metal_diffusion_coefficient`: default `0.05`; dimensionless Smagorinsky calibration coefficient.
- `metal_diffusion_cfl`: default `0.4`, accepted range `(0, 0.5]`.
- `metal_diffusion_max_subcycles`: default `128`.
- `metal_diffusion_max_rkl_stages`: default `64`, minimum `2`.
- `metal_diffusion_coefficient_floor_code`: default `0`.
- `metal_diffusion_coefficient_ceiling_code`: default `1e30`.

Diffusion requires `hierarchical_max_rung = 0`. Enabling diffusion with model
`none`, a non-positive coefficient, invalid stage limits, or invalid coefficient
bounds is rejected during typed configuration loading. See
`docs/metals_enrichment_and_mixing.md` for current AMR/MPI scope.

## P6 optional numerical and performance policies

All options use the existing `.param.txt` typed freeze, validation and normalized
dump. No new parser or budget authority is introduced. **Validation pending.**

| Key | Default | Authority / constraints |
| --- | --- | --- |
| `numerics.treepm_adaptive_acceptance_enabled` | `false` | Select adaptive relative residual accuracy plus independent guards; false retains strict reference policy. |
| `numerics.treepm_adaptive_maximum_opening_angle` | `0.25` | Finite `(0,1]` (widened from `(0,0.5]` in the October 2026 TreePM campaign; existing values remain valid); adaptive convergence safeguard, not the primary accuracy tolerance. This is an unqualified CHUI policy default. |
| `numerics.treepm_adaptive_error_estimator` | `residual_derivative` | `residual_derivative` (previous behavior: second-order screened-kernel derivative proxy) or `relative_monopole` (classic `G M l^2 <= alpha |a| r^4`, no per-node transcendental). Adaptive policy only. |
| `numerics.treepm_cutoff_node_policy` | `contain` | `contain` (previous behavior: open every node straddling `r_cut`) or `com_distance` (evaluate a straddling node at its COM distance, zero beyond `r_cut`). Adaptive policy only; the strict reference always contains. |
| `numerics.treepm_group_walk_enabled` | `false` | Opt-in group walk (one tree walk per leaf of targets, vectorized pair loop); requires `treepm_adaptive_acceptance_enabled = true`, `treepm_adaptive_error_estimator = relative_monopole`, `treepm_gaussian_pair_lookup_enabled = true`, leaf SIMD and full-residual lookup off. See `docs/tree_pm_coupling.md`. |
| `numerics.hierarchical_coarse_displacement_criterion` | `max_speed_mesh` | `max_speed_mesh` (previous behavior: fastest particle with the largest total acceleration may move at most one PM cell/split length per coarse block) or `rms_split_scale` (mass-weighted RMS speed may move `factor x min(r_s, mean interparticle spacing)`). Hierarchical KDK only. |
| `numerics.hierarchical_coarse_rms_displacement_factor` | `0.25` | Finite `(0,1]`; used only by `rms_split_scale`. |
| `numerics.hierarchical_encounter_timestep_enabled` | `true` | Hierarchical KDK only. Adds the softened-encounter rung criterion `dt <= eta sqrt(a^3 (r_nn^2+eps^2)^{3/2} / (G m_max))` (`eta = numerics.gravity_timestep_eta`, default 0.2), `r_nn` the nearest direct-pair separation from the synchronized Tree solve. `false` restores the acceleration-only rung criterion. See `docs/time_integration.md`. |
| `numerics.treepm_identical_source_tree_reuse_enabled` | `false` | Exact generation/count/build/frame/ownership/row-identity match required. |
| `numerics.treepm_topology_refit_enabled` | `false` | Certify every source inside its original leaf cell, fixed unwrap frame and row membership; rebuild on rejection. |
| `parallel.decomposition_spatial_work_enabled` | `false` | Fixed periodic `cosmo_cube`, 64 bounded SFC bins; reference feedback retained when false. |
| `numerics.hierarchical_max_rung` (existing) | `0` | `1..12`: optional collisionless DMO block KDK; 0 remains reference. |

Adaptive mode reuses `treepm_tree_relative_force_tolerance`,
`treepm_tree_relative_force_acceleration_floor`, and
`treepm_tree_opening_theta`. Missing, incompatible, tiny, zero or nonfinite
history selects geometric fallback without manufacturing an acceleration.
The kernel inequality is documented in `docs/tree_pm_coupling.md`.

The added normalized fields change the config hash even when defaults are
selected. Legacy restart files still use their existing reader/schema path,
but strict workflow config-hash checks are retained: old normalized provenance
is not silently treated as the new config. Qualification must establish an
explicit fixture/config migration if needed. No config tests were added or run
in this source-only campaign; they are required in the second campaign.

### October 2026 additions: compatibility and measured accuracy

All four keys above default to the previous behavior, so existing decks keep
their numerical results; the normalized dump gains the new lines and the
config hash changes accordingly. Measured total-force errors and timings of
each combination (clustered 64^3 box, strict-reference comparison) are in
`docs/tree_pm_coupling.md` ("Adaptive estimator and cutoff-node policy").
`configs/production_dmo_fast.param.txt` collects the qualified fast settings.

## TreePM recovery controls (2026-10-10; qualification pending)

Both keys use the typed `[numerics]` boolean parser and normalized `.param.txt`
dump; their values participate in the normalized configuration hash.

| Key | Default | Contract |
|---|---|---|
| `treepm_gaussian_pair_lookup_enabled` | `false` | Experimental cubic-Hermite Gaussian long-range radial coefficient for direct leaf pairs. Plummer softening is evaluated independently. Accepted-node monopoles/quadrupoles and MAC derivatives remain analytic. Enable only for controlled qualification. |
| `treepm_full_mac_diagnostics` | `true` | Preserve overlapping rejection accounting. `false` short-circuits the selected MAC when an independent guard already rejects, and reports `skipped_mac_evaluations`; a skipped MAC is neither passing nor failing evidence. Force acceptance policies are unchanged. |
| `treepm_full_residual_lookup_enabled` | `false` | Experimental bounded Hermite table for the complete Plummer-minus-Gaussian direct coefficient. Requires proven uniform effective pair softening, supported scales and an error-gated cell; otherwise evaluates the analytic law. Governed allocation denial also falls back. |
| `treepm_leaf_simd_enabled` | `false` | Experimental 16-source worker staging for uniform source softening; runtime-selected AVX2 distance/contribution arithmetic with portable scalar fallback. No particle-sized packed cache. |
| `treepm_radial_coefficients_enabled` | `false` | Experimental shared analytic radial value/derivatives for accepted monopoles/quadrupoles, retaining the full second moment. Adaptive MAC remains on its established analytic proxy. |
| `treepm_residual_block_size` | `64` | Production target scheduling block: exactly 16, 32, 64, 128 or 256. Memory estimates and all local/incoming traversal regions use this setting. |
| `treepm_tree_max_leaf_size` | `16` | Production tree leaf limit, integer in [1,1024]. This changes tree partition/approximation, so qualify force errors as well as timing. |

All five keys are under `[numerics]`, participate in typed validation and
normalized dumps/hashes, and preserve conservative defaults. Fast diagnostics
can skip COM geometry after independently mandatory descent;
`skipped_geometry_evaluations` identifies those unevaluated guards. Installed
lookup and AVX2 target counts describe eligible targets, not table hits or
vectorized pair counts. Requested flags alone do not prove acceleration ran.

`cosmosim_treepm_experiment_configs BASE.param.txt NEW_DIRECTORY` loads and
revalidates canonical text with the C++ loader, then writes separately named,
normalized experiments and hashes. It preserves the base physical split,
softening, IC and cosmology. `reference` and `optimized_matched` have identical
strict acceptance; the latter selects fast diagnostics. Separate `matched_*`
files attribute each kernel/reuse/refit switch. `adaptive_025` uses alpha 0.005
and angle 0.25; `adaptive_050_experimental` requires separate qualification.
Every experiment uses global KDK. Neither angle nor any new kernel is certified
by successful parsing. See the [source handoff](repair/treepm_deep_hotpath_attempt1_closure_20261012.md).

These additions change newly normalized text/hashes even when defaults are
selected. Existing strict restart/config compatibility checks remain in force;
no checkpoint hash is rewritten or grandfathered silently. Snapshot/restart
schemas are unchanged. New code remains replay-deterministic under its existing
execution contract, but reverse-octant traversal changes accumulation order
relative to earlier revisions. Compare older checksums with numerical tolerances.
No production mesh, leaf-capacity, OpenMP scheduling, acceptance, or hierarchical
rung default changes. `hierarchical_max_rung=0` remains the safe global-KDK default.

## October 25 2026 optimization campaign keys (unqualified; NOT compiled or tested)

All keys use the existing `.param.txt` typed freeze, validation and normalized
dump. Every default reproduces the previous behavior; every faster policy is
opt-in. See the [campaign record](repair/treepm_optimization_campaign_20261025.md).

| Key | Default | Authority / constraints |
| --- | --- | --- |
| `numerics.gravity_timestep_eta` | `0.2` | Finite `(0,1]`. Coefficient `eta` of the gravity criterion `dt = eta sqrt(eps a^3 / |A|)` (scale-free comoving `A`) used by the global KDK particle and gas-cell candidates and by the hierarchical rung criterion, and of the softened-encounter bound `dt <= eta sqrt(a^3 (r_nn^2+eps^2)^{3/2}/(G m_max))`. `0.2` is the value previously hardcoded in `time_coordinator.cpp`. This is **not** AREPO's `ErrTolIntAccuracy` (AREPO uses `sqrt(2 eta eps/|a|)`; equal step for `eta_CHUI = sqrt(2 eta_AREPO)`). The hydro CFL coefficient and source-term limits are unchanged. |
| `numerics.treepm_group_walk_target_size` | `0` | `0` keeps the legacy group (active targets of one source leaf). `N` in `[2, treepm_residual_block_size]` groups up to `N` consecutive tree-order targets inside one residual block, closed early when the group AABB would exceed `0.5 r_cut` on any axis. Requires `treepm_group_walk_enabled = true`. |
| `numerics.treepm_pm_force_gradient` | `spectral` | `spectral` (reference: three inverse FFTs of `-i k phi_k`) or `finite_difference_4` (one inverse FFT plus the fourth-order centered difference `-[8(phi_{i+1}-phi_{i-1})-(phi_{i+2}-phi_{i-2})]/(12h)`; no retained `phi_k` copy). Periodic host PM only, `>= 4` cells per axis, every MPI rank must own two `x` planes. Forces near the mesh scale are weaker by `T(kh) = (8 sin kh - sin 2kh)/(6kh)`; requires owner force validation. |
| `output.restart_interval_steps` | `0` | `0` keeps checkpoints coupled to snapshot events (legacy `write_restarts` behavior). `N > 0` writes a checkpoint every `N` completed global steps (coarse blocks under hierarchical KDK) and snapshot events no longer imply checkpoints. Requires `write_restarts = true`. The next checkpoint is derived from the restart-persisted step index (no restart schema change). |
| `output.restart_interval_time_code` | `0` | Finite `>= 0`. `T > 0` writes a checkpoint at the first completed global step (coarse block) ending at or after each event `t_begin + kT`; after a write the next event is the first multiple strictly after the current time. Either independent interval (step or time) disables snapshot-coupled checkpoints. Requires `write_restarts = true`. Derived from restart-persisted time on resume (no restart schema change). |
| `numerics.treepm_pm_reuse_policy` | `exact` | **Experimental, unqualified.** `exact` (reference: PM field refreshed for the exact source state; requires `treepm_update_cadence_steps = 1`), `multirate_outer_kick` (hierarchical KDK only: r-RESPA outer PM level. `treepm_update_cadence_steps` in `[2, 64]` is the *target* number of coarse blocks per PM cycle. The cycle length is fixed at opening from that target, the PM displacement criterion and the next output/end boundary. Inner blocks are re-planned by the reference criteria. An interior PM displacement violation rolls the cycle back and retries it shorter. PM is solved only at cycle ends. Output, checkpoints and the step budget are held to cycle ends (soft stop). See `docs/time_integration.md`) or `stale_within_cadence` (the long-range field is reused for up to `treepm_update_cadence_steps - 1` global steps after a refresh, while particles drift; requires cadence in `[2, 64]` and `hierarchical_max_rung = 0`). Checkpoints are deferred to the step before a refresh opportunity so resume reproduces the uninterrupted field age (`io.restart.deferred_for_pm_cadence`). Force error grows with field age; not equivalent to the reference. |
| `numerics.treepm_multirate_rollback_memory_policy` | `single_block_fallback` | Used only with `treepm_pm_reuse_policy = multirate_outer_kick`. Response when the MemoryGovernor (or an allocation) refuses the PM-cycle rollback snapshot (147 B/particle plus a small fixed part, derived from the buffer types; see `docs/time_integration.md` §5) on **any** rank; admission is agreed collectively, so every rank takes the same branch. `single_block_fallback`: that cycle opens as one coarse block on every rank. The block cannot be rejected and still honours the particle criteria and the PM displacement bound. It is recorded as cycle limiter `rollback_memory` in the block records, `multirate_memory_limited_cycles` and counter `time.hierarchical_multirate_cycle_memory_limited`. The numerical schedule then depends on the memory budget. `fail`: reproducibility mode; the run fails collectively instead, so the schedule never depends on memory. Added 2026-10-26 (round 4): changes the normalized config text/hash; no restart schema change. |
| `output.snapshot_verification` | `full_readback` | `full_readback` (materialize and compare every field) or `structural` (header-only: schema, file kind, provenance config hash, generation id, member/global per-type counts, header time/cosmology). Applies to non-collective members; collective single-file snapshots always use the bounded rank-partition readback. The SHA-256 member sidecar, completion marker and root snapshot-set validation run in both modes. |

The new keys change normalized text and therefore config hashes even when
defaults are selected. Existing restarts written before this change no longer
match the current normalized text and fail closed at resume; the narrow v23
migration projection does not cover these keys. No snapshot or restart schema
is changed.

A complete performance-oriented profile that combines these keys with the
audited TreePM settings is
`configs/performance/dmo_128_pm256_performance_unqualified.param.txt`; it is
explicitly marked as requiring scientific qualification.

## CHUI Part 2 baryonic physics keys (2026-11-03; source only, NOT RUN)

All keys live in `[physics]`. Defaults reproduce the Part 1 behaviour
(legacy cooling, total-metal tracking, legacy stellar table and deposition).
Pre-Part-2 checkpoints resume only with every key below at its default.
Handoff: `docs/repair/chui_part2_thermochemistry_stellar_physics_20261103.md`.

| Key | Values / default | Meaning |
|---|---|---|
| `thermochemistry_backend` | `legacy_cooling` \| `analytic_reference` \| `equilibrium_primordial` \| `fixed_hybrid` \| `fixed_detailed` | Part 2 backends run in the source stage and require `metal_species_mode = core_elements` and `enable_cooling = true`. |
| `metal_species_mode` | `total_only` \| `core_elements` | `core_elements` enables the 11-element composition sidecar (only with a Part 2 backend). |
| `uv_background_model` | adds `treecool_table` | Part 2 backends accept only `none` or `treecool_table` (`hm12`/`fg20` are legacy analytic proxies). |
| `uv_background_table_path`, `uv_background_table_sha256` | "" | TREECOOL file and its SHA-256. Since Part 2.5 the checksum is **required** whenever a path is set (64 lower-case hex digits), for this pair and likewise for `metal_cooling_table_*` and `ssp_yield_table_*`. Data files are thereby bound to the config hash and restart identity. `tools/chui_data/validate_treecool.py` prints the value. |
| `thermochemistry_relative_tolerance` | 1e-4, (0, 0.1) | stiff-solver relative tolerance |
| `thermochemistry_abundance_absolute_tolerance` | 1e-12, (0, 1e-4] | per-H-nucleus abundance atol |
| `thermochemistry_max_substeps`, `thermochemistry_max_newton_iterations` | 4096, 12 | step-attempt and Newton budgets; exhaustion rejects the stage (F2) |
| `thermochemistry_analytic_relaxation_time_years`, `thermochemistry_analytic_equilibrium_temperature_k` | 1e6, 1e4 | REFERENCE backend only |
| `isrf_g0_habing` | 0 | FUV field for H₂/CO photodissociation and photoelectric heating |
| `chemistry_shielding_length_max_pc` | 100 | cap on the Jeans shielding length |
| `metal_cooling_mode` | `none` \| `fine_structure_cold` \| `element_table` \| `element_table_and_fine_structure` | network backends only; table modes need `metal_cooling_table_path` |
| `metal_cooling_table_path`, `metal_cooling_table_sha256`, `metal_cooling_table_range` (`reject` \| `clamp_density_redshift`), `metal_cooling_switch_temperature_k` (1e4) | | element-resolved table (`chui_metal_cooling_table_v1`) |
| `metal_cooling_accept_uv_background_mismatch` | false | Part 2.5. A metal table must declare its embedded `uv_background`. A table is compatible only if it declares `none` and no CHUI UVB is configured, or if its `uv_background_treecool_sha256` equals the configured TREECOOL SHA-256. Any other pairing (for example WSS09/HM01 tables with an FG20 TREECOOL) is refused unless this key is `true`. Setting it `true` without an element-table `metal_cooling_mode` is refused as a no-op. The key is part of the normalized config and the config hash. Part 2 checkpoints that predate it restart only at the default. |
| `chemistry_primordial_helium_mass_fraction` | 0.2449 | initialization policy |
| `chemistry_metallicity_scaled_dust_proxy` | false | D/D_MW = Z/Z_sun when no live dust lanes |
| `dust_model` | `none` \| `two_reservoir_v1` | live silicate + carbon dust (network backends) |
| `dust_accretion_timescale_ref_years`, `dust_sputtering_timescale_ref_years`, `dust_sn_swept_mass_msun`, `dust_sn_destruction_efficiency`, `dust_agb_condensation_efficiency`, `dust_cc_condensation_efficiency` | 2e7, 1.7e8, 6800, 0.3, 0.2, 0.1 | uncalibrated research defaults |
| `ism_regime_model` | `legacy` \| `effective_or_resolved` | required with Part 2 + effective SF |
| `ism_jeans_resolution_cells`, `ism_regime_hysteresis_fraction`, `ism_regime_min_dwell_years` | 4, 0.2, 1e6 | regime criterion |
| `sgs_turbulence_model` | `none` \| `budgeted_v1` | research closure |
| `sgs_dissipation_coefficient`, `sgs_eddy_viscosity_coefficient`, `sgs_feedback_energy_fraction` | 1.58, 0.05, 0 | |
| `chemistry_adaptive_selection` | `off` \| `local_v1` | requires `fixed_detailed` |
| `chemistry_adaptive_hysteresis_fraction`, `chemistry_adaptive_min_dwell_years`, `chemistry_adaptive_cold_temperature_k`, `chemistry_adaptive_dense_nh_cm3` | 0.3, 1e6, 3e3, 1 | |
| `chemistry_shadow_sample_fraction`, `chemistry_shadow_max_cells` | 0, 0 | bounded shadow reference sampling |
| `stellar_population_model` | `legacy_table` \| `chabrier_reference_v1` \| `chabrier_yield_table_v1` | Part 2 SSP; table mode needs `ssp_yield_table_path` |
| `ssp_yield_table_path`, `ssp_yield_table_sha256` | "" | `chui_ssp_yield_table_v1` |
| `ssp_imf_mass_min_msun`, `ssp_imf_mass_max_msun`, `ssp_cc_mass_min_msun`, `ssp_cc_mass_max_msun` | 0.1, 100, 8, 40 | |
| `ssp_snia_per_msun`, `ssp_snia_delay_min_years`, `ssp_snia_delay_max_years`, `ssp_energy_per_sn_erg` | 1.3e-3, 4e7, 1.37e10, 1e51 | |
| `ssp_reference_net_metal_yield` | 0 | REFERENCE yields only |
| `ssp_stochastic_sn_events` | true | discrete events vs expectation |
| `ssp_clustered_associations` | false | RESEARCH; see the Part 3 associations table below (global friends-of-friends, decomposition independent since 2026-11-04) |
| `thermochemistry_allow_unqualified_rates` | false | Required for `equilibrium_primordial`, `fixed_hybrid` and `fixed_detailed`. After Part 2.5 the equilibrium rate set is primary-verified, but numerical validation is USER-RESERVED and pending. The network manifests still contain secondary-verified (R13/R15/R17/R20/R21) and unverified (R18) entries. See `docs/chemistry_rate_provenance.md`. |
| `fb_supernova_model` | `legacy` \| `mechanical_resolved_unresolved` | required by Part 2 populations when feedback is on |
| `fb_presn_wind_momentum_lsun_per_c` | 0 | EXPERIMENTAL pre-SN momentum closure |

Additional rule: `enable_metal_diffusion = true` is refused with a Part 2
backend (element-resolved diffusion not implemented).

## CHUI Part 3 — automatic AMR regrid (`numerics.*`)

Source only; NOT RUN — USER-RESERVED. All keys are typed, validated and
emitted in the normalized dump. Checkpoints written before Part 3 resume only
with every key at its default (exact restart projection).

| Key | Type / range | Default | Meaning |
|---|---|---|---|
| `amr_regrid_mode` | `off` \| `automatic` | `off` | `automatic` refines/derefines production AMR at synchronized hierarchy boundaries (after hydro sweep, reflux and composition commit). `off` keeps the static hierarchy. |
| `amr_regrid_interval_steps` | int >= 1 | `1` | Regrid every N accepted steps. Passes alternate refine-only / derefine-only (2:1 face balance without exchanging decisions across ranks). |
| `amr_regrid_max_operations` | int >= 1 | `64` | Deterministic per-pass budget (ascending patch id). |
| `amr_refine_cell_mass_code` | double >= 0 | `0.0` | Refine when the comoving cell mass `rho_com dx^3` exceeds it (0 disables). |
| `amr_refine_jump_threshold` | double >= 0 | `0.0` | Refine when `max(|drho|/rho_min, |dP|/P_min)` over in-patch face neighbours exceeds it (contact/shock indicator; 0 disables). |
| `amr_refine_jeans_cells` | double >= 0 | `4.0` | Truelove criterion on the comoving Jeans length `c_s sqrt(pi a/(G rho_com))` (0 disables). |
| `amr_derefine_hysteresis` | (0, 1] | `0.5` | Derefinement requires every criterion to be below `hysteresis x threshold` (Jeans: above `threshold / hysteresis`). |
| `amr_refine_cooling_length_cells` | double >= 0 | `0.0` | **Mandatory** (2026-11-04). Refine when the cooling length `c_s t_cool` spans fewer than N cells (physical lengths). `t_cool = e / (Lambda_net n_H^2 / rho)` comes from the legacy in-solver cooling-rate provider; net heating imposes no constraint. Requires `enable_cooling = true` and `thermochemistry_backend = legacy_cooling`; any other combination is refused (Part 2 backends do not publish cooling times). |
| `amr_refine_metallicity_jump` | double >= 0 | `0.0` | Optional. Refine when `max |Z_a - Z_b| / max(Z_a, Z_b, 1e-6)` over in-patch face neighbours exceeds it. |
| `amr_refine_species_jump` | double >= 0 | `0.0` | Optional. Refine when the maximum absolute jump of a species mass fraction (composition gas mass) over in-patch face neighbours exceeds it. Requires `fixed_hybrid` or `fixed_detailed`. |
| `amr_predictive_event_lookahead_years` | double >= 0 | `0.0` | R2 (RESEARCH). A star's host cell, both now and at its drifted position `x + v dt / a`, is refined when that star's own event chronology produces a core-collapse or Ia event within this window. The chronology is the cumulative-Poisson thresholds of its population record, evaluated on a copy; no ordinal is consumed and no budget changes. In expectation mode (no stochastic events) the criterion is at least 0.5 expected events, which is a heuristic. Requires `amr_regrid_mode = automatic` and a Part 2 population model. Forced rows are optional-priority refinements, and they also block derefinement. Limitation: a drifted host cell on another rank is not forced. |
| `amr_refine_vorticity` | double >= 0 | `0.0` | Optional. Refine when `|curl v| dx / c_s` (comoving gradients and widths, so the value is scale-free) exceeds it. |

**Mandatory versus optional (R1).**
- Mandatory criteria: Jeans, cell mass and cooling length. They are
  scheduled before any optional request. A cell inside their hysteresis band
  is never derefined.
- Optional criteria: the jump indicator, metallicity, species and vorticity.
  The per-pass budget truncates these first.
- A truncated mandatory request is reported in the
  `hydro.amr_automatic_regrid` event as `mandatory_truncated_by_budget`.
- Every stencil uses in-patch neighbours, so decisions are identical on any
  MPI decomposition.

`automatic` requires at least one non-zero criterion and `amr_max_level` in
[1, 30]. Regrid transfers mass, momentum, total energy and every composition
lane (elements, species, dust, SGS energy) conservatively; a refused memory
admission leaves the hierarchy untouched.

## CHUI Part 3 — stellar winds and associations (`physics.*`)

| Key | Type / range | Default | Meaning |
|---|---|---|---|
| `ssp_wind_model` | `none` \| `table_v1` | `none` | `table_v1`: age/metallicity-dependent wind/pre-SN mechanical energy and momentum from a `chui_ssp_wind_table_v1` file (no mass: wind mass is in the SSP CC continuous return). Refused together with `fb_presn_wind_momentum_lsun_per_c > 0` (double count). |
| `ssp_wind_table_path` / `ssp_wind_table_sha256` | path / 64-hex | empty | Required for `table_v1`; checksum verified at load. BLOCKED BY DATA: no table is bundled. |
| `ssp_clustered_associations` | bool | `false` | RESEARCH (EXPERIMENTAL/UNQUALIFIED): merge active feedback sources into associations BEFORE the remote stencil query. Discovery is a GLOBAL friends-of-friends over all ranks (link: minimum-image comoving separation <= radius AND birth-time difference <= window); identity = minimum member source id; members merge in id order on the leader's owner rank. Mass, elements, dust, SN/wind/pre-SN budgets and vector momentum are conserved; the members' internal velocity-dispersion kinetic energy is carried as `association_dispersion_energy_code` and deposited. Independent of the MPI decomposition. **Migration (2026-11-04):** grouping changed from rank-local lattice buckets to global FoF, so results differ from earlier revisions with this key enabled. |
| `ssp_association_radius_code` | double >= 0 | `0.0` | Comoving FoF linking length (code length); must be > 0 when associations are enabled. Semantics changed from lattice cell size to linking length (2026-11-04). |
| `ssp_association_time_window_years` | double > 0 | `3.0e6` | Birth-time linking window. |

`physics.enable_metal_diffusion` with an element-resolved composition sidecar
now selects the element-resolved mixing operator (Part 2 refused it).
Multi-rank metal diffusion remains refused (named blocker: no per-subcycle
remote interface exchange).

## CHUI Part 3 Milestone I — `numerics.amr_time_stepping`

`synchronized` (default; every AMR level advances the coordinator interval)
or `recursive_subcycling` (Berger-Oliger: level L advances dt / r^(L - L_min)
with nested temporal ghost histories and per-level reflux; requires a 2:1
balanced hierarchy, validated across ranks). **2026-11-04:** supported on any
number of ranks (directed ghost exchange per level clock, remote coarse ghosts
interpolated in time, remote register halves merged into the coarse owner's
pending store; ranks without cells traverse every level). The hydro CFL
candidate of a level-L cell is relaxed by r^(L - L_min,global); gravity,
expansion, source, diffusion and chemistry candidates are NOT relaxed (those
operators run at the coarse cadence). A level-aware accepted-step CFL guard
runs before the sweep. Event `time.timestep_limiter` records the limiting
criterion and level offset. Source only, NOT RUN — USER-RESERVED.

Metal diffusion (`physics.enable_metal_diffusion`) is supported on several
ranks for the element-resolved composition operator (unique-owner interface
fluxes); the legacy total-Z scalar operator without the composition sidecar
remains single-rank and is refused explicitly on several ranks.
