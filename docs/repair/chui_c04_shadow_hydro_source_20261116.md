# CHUÍ — C04 sampled shadow hydrodynamics: revision-4 source handoff (2026-11-17)

**Status: mandatory v10, PARTIAL source implementation, UNCOMPILED and UNVERIFIED.**
This addendum supersedes the 2026-11-16 source eligibility narrative, not the
pre-existing AMR geometry proof. Input is the user-selected
`Chui-dev-main(4)(2).zip`; the attached instruction pack still names archive
`Chui-dev-main(3).zip`, which was not used as the implementation baseline.
No compilation, tests, simulations, benchmark, static analysis or scientific
validation occurred during this code handoff.

## P1-C04-01 — Production numerical policy (source implemented)

Fixed-grid (`HydroAmrRuntimeImpl`) and AMR (`ProductionAmrHydroOptions`/
`advanceProductionAmrHydro`) pass the actual `HydroReconstructionPolicy`,
`HydroNumericsPolicy`, integration/reconstruction choice and Riemann solver kind.
The private reference uses CTU+PPM, with identical limiter, adiabatic gamma,
density/pressure floors, MUSCL predictor contract, dual-energy controls,
HLLC/HLLE choice and fallback machinery. Its CTU+PPM versus CTU+PPM path
advances two halfsteps, recomputing `dt/(a dx)` and a/H at substep midpoints.
CTU/PLM and MUSCL/PLM versus reference CTU/PPM are **mixed spatial/temporal**
comparisons and must not be interpreted as pure Richardson error estimates.
The seven reported channels include the conservative metal lane and optional
auxiliary internal-energy lane. Comparison is of local volume-integrated
*disagreement*, not a global conservation residual or external error bound.

## P1-C04-02 — Coupled sources/cosmology (partial)

`HydroSourceTerm::cloneForShadow()` is a restricted private replay protocol.
Cooling/heating uses a read-only adapter over the production's immutable,
qualified lookup source; effective-ISM energy relaxation gets a private source
object and ledger. Nonclonable sources return null, with
`skipped_unsupported_source`. The sampled source context reuses the selected
thermodynamic closure read-only and snapshots the resolved/effective gate.
A private `HydroCellSourcePropertyProvider` derives changing mass/hydrogen
physical densities, thermal energy/temperature and metallicity from the
reference conservative state and sampled initial conversion, including
cosmological a⁻³ dilution. This is an **approximation** where the production
microphysical source uses nonlocal chemical/ion histories or temperature
conversion dependent on evolving molecular weights; those contexts must NOT
be advertised as fully source-equivalent. Substep redshifts are determined by
the independently advanced scale-factor midpoint, not the production endpoint.

**OPEN:** external Strang halfsteps, independently evolving element/ion
composition, discrete/stochastic star/feedback events, global gravity coupling (including `ComovingGravityExpansionSource`,
which is refused until private gravity acceleration histories exist),
and stiff source-driven event/barrier chronology do not have a clone/replay
contract. The existing composition face-mass-flux sink explicitly refuses the
sample (`skipped_external_passives`). An operator without a private clone
refuses a sample instead of mutating its live state. External stage effects
are *outside the hydro interval* and are not separately replayed. The
currently implemented reference must not be claimed to cover full baryonic
stage-level source fidelity.

## P1-C04-03 — Small AMR patches and ghosts (partial)

For self-contained periodic/open/reflective Cartesian patches, all ghosts
have local real-cell source/owner indices and writable-physical status; the
reference refreshes conserved and auxiliary-energy ghosts from private
state at **every reference substep**. This provides an 8³ sampled path
without falsely relaxing the stencil guard or making artificial diagonal
flux-register faces. All other patches continue with the existing 4/7-cell
interior guard. An 8³ patch whose exterior is imported through same-level
MPI/diagonals or coarse/fine AMR is **not supported** without an independent
midpoint temporal history. If no interior remains, report
`skipped_temporal_ghost`, not a zero error or low-order substitute.
`AmrTemporalBoundaryHistoryRecord` from the live hierarchy is not reused
as an evolving private-neighborhood authority.

**OPEN:** bounded multi-patch neighborhood and/or independently reconstructed
midpoint coarse/fine and diagonal MPI ghost histories and their conservative
auxiliary/face bookkeeping. This is mandatory to close C04 for ordinary
multi-patch 8³ AMR hierarchies; the current source does NOT satisfy that scope.

## P2-C04-04 — Distributed sampling (partial)

At a rank-uniform hydro boundary, each rank holds at most K candidates,
ranked by a stable 64-bit score of (patch ID, level, step, seed), with a
stable ID/level tie-breaker. Only these bounded candidates are exchanged,
then globally top-K ranked; empty ranks join every collective. A selected
identity is advanced at most once per epoch; local byte/cell limits can
produce `skipped_resource` rather than substitute a different patch.
`hydro_shadow_max_patches_per_rank` retains its old key name but now means
global AMR K; per-rank memory/cell caps retain their old scope. Fixed-grid
sampling still uses local region ID 0 on all ranks, so it is **not** invariant
to an alternative fixed-grid decomposition. Patch geometry/migration identity
uniqueness still requires owner MPI tests.

## P2-C04-05 — Memory and non-invasiveness (source implemented; unverified)

Scratch/clone/candidate capacity is computed with checked 64-bit products and
sums before governor leases. The global selection's byte cap, rank-wide
resource refusal vote, and bounded candidate wire prevent a hidden full
patch-inventory allgather. Local shadow failures are recorded; all MPI
reductions/publish remain at rank-uniform stage boundaries, not in
conditional per-patch reference advances. Reference copies and scratch are
phase-local and discarded; no authoritative gas, flux register, scheduler,
AMR history, stochastic source, restart or SHAI-HULUD research file is edited.
When shadow is disabled, the path returns before candidate allocation.
`observed_capacity_high_water_bytes` and `workspace_high_water_bytes` have
different meanings (observed local reference capacity vs admitted reservation).
Budget/fallback behavior requires the owner's eventual MPI/memory execution.

## P2-C04-06 — Test sources, provenance and verification (partial)

`tests/unit/test_hydro_shadow_reference.cpp` covers MUSCL/CTU paths, cadence,
budget and disabled/no-sample cases. The new
`tests/unit/test_hydro_shadow_c04_coupled.cpp` adds 3D 8³ periodic/reflective
coverage, non-default limiter and floor, smooth source/closure, evolving
comoving background, dual-energy ghost lanes, decomposed top-K identity
simulation and imported-ghost refusal. CMake registers the new source-only
test. These have **not been compiled or run**. True MPI diagonal and
coarse/fine ghost fixtures, external Strang/discrete-event, restart and
scientific convergence studies remain unimplemented/unrun and must not be
listed as passed.

No new restart or migration fields, HDF5 datasets, user-facing configuration
keys, flux-register owner or physical model were added. The
`hydro.shadow_reference` event documents selected schemes, limiter, floors,
Riemann solver, global K/seed, comparison axis, completion and separate skip
reasons; `no_samples` remains explicit. Source-level status remains PARTIAL
until the missing C04 interfaces, fixtures and owner verification are closed.

## AMR source-contract closure

- `src/workflows/internal/amr_geometry_assurance.cpp`: the pre-source lenient policy collectively refuses any mixed population of legacy all-zero and explicit geometric patches, even when the mixture is split across ranks. It continues accepting all-legacy non-AMR states and ignores inactive parents.
- `include/cosmosim/amr/amr_geometry_coverage.hpp` and `src/amr/amr_geometry_coverage.cpp`: maximum admitted lattice span is now `2^35`. The existing `1e-6 + 64 ε span` tolerance is thereby strictly less than `0.001` finest-cell widths across all admitted spans, keeping integer signed-corner proof and fail-closed oversized-span handling intact.

**Classification:** AMR P2 source defects addressed, subject to tests. C04 contains a real independently advanced research reference and production hooks, with a delimited supported physics/geometry envelope; full coupled/complex-AMR shadow source-completion remains open if the full v10 requirement is interpreted to require every mode/physics combination.
