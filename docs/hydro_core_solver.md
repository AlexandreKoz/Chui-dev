# Hydrodynamics core solver (finite-volume Godunov scaffold)

> Part 1 (20261103): the reconstruction/integration/Riemann method matrix, the
> corrected Hancock predictor, CTU, PPM, dual energy and the passive-scalar
> contract are specified in [`hydro_numerics.md`](hydro_numerics.md), which
> supersedes conflicting statements below.

## Variable definitions

The hydro core uses comoving conserved variables per cell:

- `mass_density_comoving = rho`
- `momentum_density_*_comoving = rho * u_*_peculiar`
- `total_energy_density_comoving = rho * (e_internal + 0.5 |u_peculiar|^2)`

The explicit update implemented in `HydroCoreSolver::advancePatch` is:

\[
U^{n+1}_i = U^n_i - \frac{\Delta t}{a V_i} \sum_f A_f \hat{F}_f + \Delta t\, S_i.
\]

with `a = scale_factor`, `V_i = cell_volume_comoving`, and `A_f = face.area_comoving`.

## Conservation diagnostics

`HydroConservationTotals` reports volume-integrated diagnostics over the selected real cells:

- mass: `sum rho * V`
- momentum x/y/z: `sum rho * u_* * V`
- total energy: `sum rho * (e_internal + 0.5 |u|^2) * V`
- internal energy: `sum (E_total - 0.5 |m|^2 / rho) * V`

The internal-energy entry is derived from the conserved total-energy and momentum densities. It is not a second
evolved energy authority.

`HydroProfileEvent::conservation` records the selected-cell totals before and after the hydro update, plus separate
volume-integrated deltas from conservative face fluxes, explicit source terms, and internal-energy floor corrections.
The residual is:

\[
U_\mathrm{after} - U_\mathrm{before} - \Delta U_\mathrm{flux} - \Delta U_\mathrm{source} - \Delta U_\mathrm{floor}.
\]

This keeps deliberate source-term work and momentum injection out of the flux-conservation error budget. For a closed
periodic full-patch update without source terms or floor corrections, the flux deltas and residuals for mass, momentum,
total energy, and the derived internal-energy diagnostic are expected to remain within the test tolerance. Active-set
updates report the budget over the caller-selected active cells; flux across the active/inactive boundary is therefore
physical active-set exchange, not automatically a conservation error.

If the update would leave non-positive internal energy, the solver raises total energy to the kinetic energy plus the
small positive floor and increments `internal_energy_floor_count`; the resulting energy addition is reported in
`floor_delta` rather than hidden inside flux or source accounting.

## Stage structure and ownership

The update is split into stages to keep AMR reuse and future GPU kernels straightforward:

1. Reconstruction (`HydroReconstruction`): produce left/right primitive states per face.
2. Riemann (`HydroRiemannSolver`): compute face-centered numerical fluxes.
3. Flux accumulation (`HydroCoreSolver`): apply conservative owner/neighbor updates.
4. Source application (`HydroSourceTerm`): apply gravity/expansion and optional future sources.

During the Riemann stage, callers may pass a `HydroFluxRegisterSink`. `HydroCoreSolver` remains AMR-neutral: it only
checks `HydroPatchGeometry::flux_register_faces` for precomputed per-face metadata and emits the actual numerical flux
returned by the configured `HydroRiemannSolver`. Faces without flux-register metadata do not call the sink, so same-level
and ordinary physical-boundary faces do not produce reflux records.

The flux-register sign convention is oriented to the coarse cell face named by the metadata. `HydroFace::normal_*`
continues to point from the hydro face owner toward its neighbor/ghost, and the finite-volume owner update subtracts
`F * A * dt / (a V)`. For reflux diagnostics, `coarse_orientation_sign = +1` means the Riemann flux is already aligned
with the outward normal of the registered coarse face; `-1` flips the flux before it reaches the sink. Fine faces that
are swept from the fine patch side of a coarse-fine boundary typically use `-1` when their local outward normal points
back toward the coarse ghost. The AMR accumulator area-weights all fine-face fluxes sharing the same register key and
emits one deterministic `FluxRegisterEntry` ordered by that key.

In the production AMR hydro path, refluxing is applied automatically by `advanceProductionAmrHydro` after all active
coarse and fine patch sweeps for the current hydro update have contributed to the shared flux register accumulator and
after patch-local conserved states have been scattered back to authoritative `SimulationState` gas rows. The reflux
application happens before AMR restriction/derefine decisions, output checks, or restart-boundary writes, so the
coarse-fine flux mismatch is corrected in the live owner cells rather than left for test-only callers. Reflux targets
are validated by patch ID, stable `gas_cell_id`, and patch-local cell index; stale coarse patch/gas-cell mappings fail
before any correction is applied.

Implementation is now split along these boundaries:

- `src/hydro/hydro_reconstruction.cpp` + `include/cosmosim/hydro/hydro_reconstruction.hpp`:
  limiter library, piecewise-constant reconstruction, and axis-aware MUSCL-Hancock predictor/floors.
- `src/hydro/hydro_riemann.cpp` + `include/cosmosim/hydro/hydro_riemann.hpp`:
  HLLE and HLLC flux construction with HLL fallback on HLLC degeneracy.
- `src/hydro/hydro_core_solver.cpp` + `include/cosmosim/hydro/hydro_core_solver.hpp`:
  patch update scaffold (validation, cache fill, conservative flux accumulation, source-term pass, profiling).
- `src/hydro/hydro_boundary_conditions.cpp` + `include/cosmosim/hydro/hydro_boundary_conditions.hpp`:
  explicit Cartesian physical-boundary ghost metadata and transient ghost-state fill.

Hydro hot data remains SoA in `HydroConservedStateSoa`. Cold patch metadata is isolated in `HydroPatchColdData`.
For hierarchical timestepping, `HydroActiveSetView` allows explicit active cell/face subsets via
`advancePatchActiveSet` without rewriting geometry ownership.
For memory-sensitive loops, `HydroScratchBuffers` and `HydroPrimitiveCacheSoa` can be reused across calls via
`advancePatchWithScratch` / `advancePatchActiveSetWithScratch` to avoid repeated allocations and repeated
primitive reconstruction for piecewise-constant paths.

Production workflow hydro validates dense `GasCellIdentityMap` coverage before loading or storing gas-cell state.
Hydro state is cell-local: parentless cells load velocity and thermodynamics from `GasCellSidecar`, and split cells
may share one optional parent particle without making that parent the gas-cell identity authority. Particle mass and
velocity lanes are updated only as compatibility mirrors through explicit optional-parent lookup, and a shared parent
is mirrored at most once per hydro store pass. Conservative ghost flux corrections are keyed by stable `gas_cell_id`;
`parent_particle_id` is lineage metadata only.

AMR patch hydro geometry is supplied by `cosmosim/amr/amr_hydro_geometry.hpp` as an adapter layer around
`HydroPatchGeometry`. The adapter maps one `amr::PatchDescriptor` to patch-local real cells, ghost placeholders, and
face metadata keyed by stable `gas_cell_id`, while `HydroCoreSolver` continues to consume only the geometry-neutral
hydro types. Persistent hydro truth remains in `core::SimulationState`; AMR patch geometry views capture the
`GasCellIdentityMap` generation and must be rebuilt after gas-cell identity remaps.

AMR same-level and coarse-fine ghost state is filled by `cosmosim/amr/amr_ghost_fill.hpp` after patch-local conserved
state has been loaded from `SimulationState`. Physical AMR boundaries delegate to the H1 Cartesian fill rules
(`periodic`, `open`, and `reflective`). Same-level AMR ghosts are copied from the geometrically adjacent source patch
cell at the same refinement level. Fine-side coarse-fine ghosts use piecewise-constant coarse-to-fine injection from
the adjacent coarse cell. Coarse-side coarse-fine ghosts use a monotone arithmetic average of the adjacent fine boundary
cells that tile the complete coarse ghost volume, including the normal case where that stencil spans multiple fine child
patches. Contributors are validated as one aligned, non-overlapping, gap-free refinement stencil and accumulated in
geometric stencil order; this makes the fine boundary state available for reconstruction/restriction tests without
replacing flux registers or refluxing. AMR ghost writes are scratch-only writes into the patch-local
`HydroConservedStateSoa` ghost rows. Same-level, coarse-fine, and remote imported ghost sources are read-only; stale
remote ghost epochs are rejected before fill.

## Cartesian patch geometry

H1 production fixed-grid hydro uses `makeCartesianPatchGeometry(...)` from
`include/cosmosim/hydro/hydro_cartesian_patch.hpp` to build a dense local Cartesian patch. The patch contract is:

- row order is `linearCellIndex(i, j, k) = i + nx * (j + ny * k)`;
- `cellIjk(row)` and `neighborCell(row, di, dj, dk)` are the explicit indexing helpers;
- `HydroPatchGeometry` carries `nx`, `ny`, `nz`, origin, axis cell widths, and one uniform cell volume;
- internal faces are built in x, y, and z only, with axis metadata, explicit owner-minus and neighbor-plus
  stencil cells, and normals pointing from owner to neighbor;
- face areas are transverse products: `dy dz`, `dx dz`, and `dx dy`.

The reference workflow first tries to derive `nx`, `ny`, `nz`, origin, and spacing from rectilinear 3D
`CellSoa` center metadata with row order matching the helper. If that metadata is not available, fixed-grid toy
states use exact near-cubic factors of `cell_count` over the configured axis-aware box lengths. Geometry caching is
invalidated by cell count, dimensions, origin, widths, timestep used by the reconstruction predictor, hydro boundary
policy, and patch metadata signature.

`makeCartesianPatchGeometry(...)` builds only internal owner-neighbor faces. Runtime physical-boundary faces are added
explicitly by `appendCartesianBoundaryGhostFaces(...)`, which appends one transient ghost row per lower/upper Cartesian
boundary face in x, y, and z. Boundary faces keep the real cell as `owner_cell`; `neighbor_cell` points at a transient
ghost row after the real-cell range, and `ghost_cell_slot` indexes `HydroPatchGeometry::ghost_cells`.

Each `HydroGhostCell` records:

- the authoritative owner real cell;
- the source real cell used to fill the ghost state;
- the transient ghost row and slot;
- boundary kind (`periodic`, `open`, `reflective`, or `imported_mpi`);
- face axis and lower/upper side;
- mutation rights.

Physical-boundary ghost rows are `kWritablePhysicalBoundaryScratch`: they may be filled from real-cell state before
reconstruction, but they are not authoritative cells and are never written back to `SimulationState`. Imported MPI ghosts
remain distinguished as `kImportedMpi` / `kReadOnlyImported`; the physical-boundary fill helper intentionally skips them.
Imported ghosts carry the remote `origin_rank`, stable `origin_gas_cell_id`, hydro sync epoch, and gas-cell identity
generation used as the decomposition epoch. They are valid for one hydro stage only. A later stage, migration, regrid, or
gas-cell identity remap must rediscover and refill them before reconstruction.

The fill rules are:

- periodic ghosts copy primitive/conserved state from the opposite interior side;
- open/outflow ghosts copy the boundary-adjacent real cell;
- reflective ghosts copy density, pressure, and tangential velocity while reversing only the velocity component normal
  to the face axis.

`HydroCoreSolver` validates that active cells and face owners are real cells and that any explicit face stencil rows are
inside the conserved/primitive storage. Ghost rows may be consumed as reconstruction neighbors through active boundary
faces, but flux deltas accumulated into ghost rows are scratch only because source/update writeback iterates active real
cells. A single-cell patch therefore receives six physical-boundary ghost faces after boundary append and remains valid
without persistent ghost cells.

## MUSCL-Hancock reconstruction contract

`MusclHancockReconstruction` reconstructs the full primitive vector
`{rho_comoving, vel_x_peculiar, vel_y_peculiar, vel_z_peculiar, pressure_comoving}`. It does not infer Cartesian
direction from row contiguity. Instead, each `HydroFace` carries:

- `axis`, selecting x, y, or z directional reconstruction;
- `owner_minus_cell`, the owner-side upwind stencil row when available;
- `neighbor_plus_cell`, the neighbor-side downwind stencil row when available.

Missing stencil rows are treated as a zero local slope on that side, which is the deliberate boundary/degenerate-patch
fallback. Boundary faces still use their ghost-filled `neighbor_cell` state from `appendCartesianBoundaryGhostFaces(...)`
and `fillHydroBoundaryGhostCells(...)`.

For distributed fixed-grid workflow hydro, `HydroStageCallback` can convert a physical boundary ghost face into an
imported MPI face when compact boundary-cell geometry advertisements show an adjacent remote owner cell. The imported
row is filled from a requested conserved payload containing only the required hydro stage fields:
`mass_density_comoving`, all three momentum-density components, and `total_energy_density_comoving`. Face authority is
deterministic: the lower owner rank advances the inter-rank face, writes the local owner delta, and accumulates the
opposite-side delta in the imported ghost row. The workflow restores that read-only ghost row after the solve and sends
the conservative delta back to the authoritative owner keyed by stable `gas_cell_id`. The higher-rank peer may import the
same boundary state for lifecycle/diagnostic coverage, but it does not advance that interface face.

The predictor CFL is axis-aware through
`HydroReconstructionPolicy::dt_over_cell_width_code = {dt/dx, dt/dy, dt/dz}`. The legacy scalar
`dt_over_dx_code` remains a compatibility fallback when the per-axis entry is zero. Directional evolution uses the
velocity component normal to the face for density, normal momentum, transverse velocity advection, and pressure updates;
transverse velocity components are reconstructed and predicted rather than copied from cell centers.

The reference workflow rebuilds the reconstruction policy whenever the accepted timestep or local Cartesian geometry
changes, so `dt_over_cell_width_code` cannot remain cached from an older scheduler bin. Before
`HydroCoreSolver::advancePatchActiveSetWithScratch(...)` is called, `HydroStageCallback` checks the accepted
integrator timestep against the local directional CFL bound
`min_axes(C_cfl * cell_width_axis / (abs(v_axis) + sound_speed))` for the active gas cells. If the accepted timestep is
larger than the current hydro CFL bound beyond tolerance, the callback throws before mutating persistent gas state.
The hydro solver itself does not clamp `dt`; scheduler-visible `HydroCfl` candidates and profiler-visible
`HydroCflDiagnostics` carry the proposal and worst-cell identity.

## Source-term convention

`ComovingGravityExpansionSource` supplies conservative-source hooks using
`gravity_accel_*_peculiar` and `hubble_rate_code`:

- `S_rho = 0`
- `S_m = rho A/a^2 - H m`
- `S_E = rho (u · A)/a^2 - H (2 K + 3 P)`

Here `rho` is comoving mass density, `m=rho u`,
`K=0.5 rho |u|^2`, `P` is comoving pressure, and `u` is physical peculiar
velocity. Despite the compatibility field stem, `gravity_accel_*_peculiar`
contains the scale-free comoving TreePM kernel `A`, not the already-scaled
peculiar acceleration. The source applies `A/a^2` exactly once for both
momentum and gravity work.

`HydroUpdateContext::scale_factor` and `hubble_rate_code` are derived views of
integrator timeline truth. Hydro runs post-drift before `IntegratorState`
commits the step, so both fixed-grid and AMR callbacks use
`StepContext.timeline_step.scale_factor_end` and `hubble_end_code`, already in
inverse code time. They do not read the still-step-begin committed epoch or
recompute an SI-valued `H(a)` independently. The expansion energy term above
is valid for every supported adiabatic index. It must not be replaced by
`-2 H U` without first proving the corresponding variable and equation
convention.

This is intentionally conservative and modular; additional baryonic sources should be implemented through extra `HydroSourceTerm` instances.

M2A.1 also permits `HydroSourceContext` to provide thermodynamic/cooling inputs through a read-only `HydroCellSourcePropertyProvider`. The fixed-grid runtime uses this seam to derive physical density, hydrogen number density, metallicity, and temperature for the cell currently being cooled rather than allocating four full-grid derived vectors. The legacy span fields remain a compatibility API, and the cooling source gives the provider precedence when one is supplied. The provider changes data lifetime only; it does not change the cooling-rate model or source integration.

## Assumptions

- Face normals are unit vectors.
- Cartesian patch geometry uses a uniform `cell_volume_comoving` per dense local patch.
- `advancePatch` updates the full patch; `advancePatchActiveSet` updates only caller-selected cell/face subsets.
- Boundary handling is explicit for Cartesian physical boundaries: append boundary ghost faces, fill ghost rows, then
  reconstruct using those ghost-filled states. A face with no neighbor remains a legacy/internal fallback that mirrors
  owner state.
- `MusclHancockReconstruction` consumes explicit axis/stencil metadata from `HydroFace`; hand-built geometries that omit
  stencil rows get boundary-style zero slopes on the omitted side rather than row-difference inference.
- The provided integration test uses a closed 3D Cartesian Sod-like setup to exercise x/y/z face construction while
  preserving conservation.
- The classical hydro validation executable (`validation_hydro_classics`) uses the same Cartesian patch, scratch/cache,
  reconstruction, Riemann, conservation-total, and source-term paths for CI-scale Sedov, Noh, Gresho vortex,
  Kelvin-Helmholtz, and Evrard-style collapse guards. These cases are runtime-credibility sentinels for positivity,
  finite state, conservation/source-energy bounds, symmetry, and qualitative physical trends; they are not final
  publication-grade convergence or reference-profile campaigns.
