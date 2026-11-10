# Hydro numerics contract (CHUI Part 1)

Status: **source implementation complete; numerical qualification NOT RUN.**
Every statement about accuracy below is a design contract backed by source
tests that were written but not executed by the authoring agent. Mark any
claim "verified" only after the human-run checks in `docs/validation_plan.md`.

The method has three orthogonal selections (`docs/configuration.md`):

| dimension | options | default |
|---|---|---|
| reconstruction | `piecewise_constant`, `piecewise_linear` (limited MUSCL), `piecewise_parabolic` (CW84 PPM) | `piecewise_linear` |
| time integration | `muscl_hancock` (directionally unsplit, 1 solve/face), `ctu` (corner transport upwind, 2 solves/face) | `muscl_hancock` |
| Riemann solver | `hllc` (HLLE fallback), `hlle` | `hllc` |

All combinations are mathematically compatible: reconstruction supplies
spatial interface values and applies exactly one forward Hancock half step;
CTU adds transverse conservative increments on top of those states; the
Riemann solver turns the final states into a flux. No combination stacks two
half-time predictors.

## 1. State, units and the update

Authoritative gas state (owner `core::SimulationState`, see
`docs/architecture/runtime_truth_map.md`): comoving mass density `rho`,
peculiar velocity `u`, specific internal energy `e`, metal mass. The solver
works on its own staged copy `HydroConservedStateSoa` with lanes
`(rho, rho u, E = rho (e + |u|^2/2), rho Z)` (field schema in
`hydro_core_solver.hpp`; every lane goes through every flux branch).

Conservative update with `a` the flux epoch (interval midpoint in cosmological
runs, 1 otherwise), `V` the comoving cell volume, `A_f` the comoving face area:

    U_i^{n+1} = U_i^n - dt/(a V_i) sum_f A_f F_f(W_L, W_R) + dt S_i(U^n)

Expansion drag and gravity are gravity-owned KDK half kicks; hydro sources are
the local thermochemical terms (cooling, effective-ISM relaxation).

**Transport Courant factor.** Every predictor uses
`lambda_axis = dt / (a dx_axis)` built by `hydroTransportCourantFactors()` at
the same `a` as the flux update. (Before Part 1 the runtime passed `dt/dx`
without `1/a`: audit F1 follow-up, now fixed in the fixed-grid runtime and in
`amr_hydro_orchestrator.cpp`.) The CFL selection uses the proper crossing
length `a dx` in comoving runs (`computeDirectionalCflTimeStep`).

## 2. MUSCL-Hancock predictor (audit F1)

Stencils are ordered along the face normal `n` (owner -> neighbor = +n;
`HydroFace` documents the convention). With `Delta q` the limited slope along
`n`, `u_n = u . n`, and the cell-centre state `W_i`:

    W_face = W_i +/- Delta W / 2 + (lambda/2) D(W_i, Delta W)
    D rho = -(u_n Delta rho + rho Delta u_n)
    D u   = -u_n Delta u - n Delta P / rho
    D P   = -(u_n Delta P + Gamma1 P Delta u_n)
    D e   = -(u_n Delta e + (P/rho) Delta u_n)
    D Z   = -u_n Delta Z

`+` is the owner side, `-` the neighbor side, `Gamma1 = rho c^2 / P` from the
cell's closure. The defective version applied `-0.5 lambda` to an already
negative derivative (time-reversed advection) and evaluated the pressure
acceleration with the wrong sign; for the audit ramp `[1,2,3,4]`, `u_n = 1`,
`lambda = 0.2` the owner face density is now `2.4` (was `2.6`). Face states are
thermodynamically closed: ideal-gas states get `e = P/((gamma-1) rho)` and
`c^2 = Gamma1 P / rho`.

**Positivity.** If a predicted state is inadmissible (non-finite, `rho` or `P`
below floor) the side falls back to the unpredicted limited state, then to the
first-order cell state, and only if the cell average itself violates the
floors is it clamped. Each fallback increments `positivityFallbackCount()`.
Fluxes are never zeroed.

## 3. PPM (Colella & Woodward 1984)

Interface values `a_{i+1/2} = q_i + (q_{i+1}-q_i)/2 - (dq_{i+1}-dq_i)/6` with
MC-limited `dq` (fourth-order for smooth data), the CW84 monotonicity
constraint (extrema -> constant; overshoots steepened to `3q - 2a`), and
per-cell CW84 flattening `chi = clamp(10(|dp1/dp2| - 0.75), 0, 1)` when
`|dp1|/min(p) > 0.33` and the flow is compressive, blending toward the limited
PLM edges. Omitted by design: contact steepening (multi-D artifacts) and the
upwind `max(chi_i, chi_{i+-1})` (needs a 7-cell stencil). Time evolution is the
same single Hancock step using the parabola's mean gradient `a_R - a_L`.

Stencil: cells `i-2 .. i+3` for face `i+1/2` (`HydroFace::owner_minus2_cell`,
`neighbor_plus2_cell`). Production AMR completes these stencils through a
normal ghost shell whose depth is selected from the reconstruction; periodic
ghosts resolve to their wrapped real cell, so periodic boxes (fixed grid and
box-spanning AMR patches) are full order everywhere. With the required shell:

| face position | PLM | PPM |
|---|---|---|
| interior, >= 2 cells from a non-periodic patch boundary | full | full |
| interior face adjacent to the boundary cell | full | full |
| patch-boundary face (owner real, neighbor ghost) | full | full |

The production AMR/MPI planner requests `requiredStencilHalfWidth()` normal
layers (one for PCM, two for PLM and three for PPM), and the ghost filler
resolves each layer at the step-start/temporally interpolated level clock.
**Domain versus missing data (P1-55-02, P2-55-03).** A stencil-only ghost
(deep normal layer, edge stencil ghost, or deep reflective layer of a patch
thinner than the mirror depth) is described by the target-relative index of
the stencil cell it represents. On every non-periodic axis on which that
index leaves the patch, the global boundary is classified from the declared
hierarchy domain (`ProductionAmrHydroOptions::domain_{origin,extent}_comoving`,
production: the collectively reduced bounding box of all populated leaf
patches) or, when no domain is declared, from a target face that is itself
`kPhysical`. **Supported domain contract (patch 56 closure):** the hydro domain
is the root tiling built from the gas-cell layout, not `cosmology.box_size_*`
(the gravity/PM extent of isolated modes); production AMR therefore requires
the populated leaves to tile their bounding box completely.
`globalAmrHierarchyDomain` (`hydro_amr_runtime.cpp`) sums the leaf volumes
collectively and refuses, before any ghost fill and identically on every
rank, a hierarchy whose covered volume differs from the bounding-box volume
(relative 1e-9): internal holes or non-rectangular tilings are never treated as
physical walls. Only beyond
such a wall is the physical closure applied: open -> zero-gradient onto the
last in-domain row (`deep_layer_extrapolated_ghosts`); reflective -> mirror
`i' = 2w - 1 - i` with that axis' momentum reversed, repeated for a domain
narrower than the stencil (`deep_layer_reflected_ghosts`). The mapped
in-domain cell is then resolved as a real target cell or through the
ordinary local / coarse-fine (temporal) / fine-coarse / MPI selection, so a
thin-patch mirror that lands beyond the opposite face is read from the
neighbour patch there. A stencil cell that is (or is not provably out of)
inside the domain and has no source fails closed
(`in_domain_missing_source_rejections`); it is never continued from a
shallower ghost. A face with no neighbour that is not on a declared domain wall
is rejected as missing interface metadata
(`physical_face_inside_domain_rejections`); closures that need more than eight
reflections per axis or an unsupported hierarchy boundary kind are rejected
(`unsupported_closure_geometry_rejections`). Reflective deep physical layers of
thin patches are no longer omitted. `edge_stencil_degenerate_closures` is
retired (always 0).

## 4. CTU (simplified Colella 1990 in 3D)

1. Pass 1: reconstruct every active face with the normal-only Hancock
   predictor and solve the Riemann problem.
2. Accumulate per touched cell and per axis the conservative increment
   `dU_d(c) = -dt/(a V) sum_{faces of c normal to d} (+/-) A F` (equal and
   opposite on shared faces; periodic ghost faces land on the wrapped cell).
3. Pass 2: for a face normal to `d`, add `0.5 sum_{d' != d} dU_{d'}(c)` to each
   side's conserved state (`c` = owner or neighbor; periodic ghosts use their
   source cell), convert back with the same EOS (effective-ISM states keep the
   barotropic closure linearized: `P' = P + c^2 (rho' - rho)`), and solve again.
   Only pass-2 fluxes update the state and the AMR flux registers.

Base Riemann calls: 2 per face (6 per cell in 3D, 4 in 2D). A non-periodic
interface ghost also evaluates two transient Riemann problems per available
tangential direction to reconstruct its pass-1 divergence. These extra solves
are included in `riemann_solve_count` and `ctu_transverse_ghost_axes`. Stability: Courant
`<= 1` in 2D, `<= 1/2` in 3D (no Saltzman 12-solve corner terms); the gas CFL
candidate is `dt = 0.4 min_d a dx_d/(|v_d|+c)`, which satisfies it. With
first-order states the scheme reproduces the exact donor-cell corner-transport
update (unit test `test_hydro_ctu`).

At patch boundaries the owner side uses its ordinary pass-1 divergence. For a
non-periodic ghost side, the solver uses the immutable AMR/MPI face slab and
the selected PCM/PLM/PPM reconstruction to evaluate the lower/upper flux in
each tangential direction on the target control volume. The term is accepted
only if the complete selected-order stencil exists.

**Edge stencil ghosts (P1-HYDRO-02).** Production AMR with CTU builds, per
patch, a stencil-only *edge ghost block*: for every face (axis `n`, side) and
transverse axis `t` with extent > 1, the cells one layer beyond the face and
`1..r` cells beyond the transverse patch face (`r` = stencil half-width =
ghost layers; PCM 1, PLM 2, PPM 3), keyed by position
(`HydroPatchGeometry::edgeGhostSlot`). They own no face and no flux register
and are never advanced (read-only imported scratch). The tangential slab of a
face ghost resolves `g + o e_t` to the standard shell inside the patch and to
the edge block beyond it, so the slab is complete at every patch edge and
corner (the corner region is covered by the two one-dimensional slabs; no
three-axis corner cell is read). The AMR ghost filler resolves each edge
ghost with any legal level relation (same level, coarse→fine with the
temporal history, fine→coarse conservative average), local or sparse MPI,
periodic wrap (including in-patch wraps), and applies the physical domain
closure on every crossed non-periodic domain face (open: zero-gradient;
reflective: mirror with reversed normal momentum). Missing source data fails
closed. Local ghost sources include edge/corner-adjacent patches, and the
directed MPI exchange widens the target's transverse intervals by `r` cells
(`AmrInterfaceTopologyOptions::edge_stencil_layers`) so diagonal remote cells
are delivered. Telemetry: `edge_stencil_ghosts_filled`,
`edge_stencil_physical_closures` (closures use the domain-versus-missing-data
rule of section 3; `edge_stencil_degenerate_closures` is retired). The CTU
ghost-slab scratch is sized by the face-adjacent ghost prefix only
(`HydroPatchGeometry::face_ghost_slot_count`, P2-55-04); deep and edge slots
never receive a slab correction and hold no entry.

With the edge block present an incomplete slab is a geometry defect
(`std::logic_error`), and the one-sided owner proxy is disabled. Geometries
without the block (fixed-grid / custom, or a ghost-layer override below the
stencil width) keep the counted closure: an imported ghost takes the adjacent
real owner cell's transverse increment
(`ctu_transverse_boundary_proxy_states`), a physical open/reflective ghost
keeps the normal-only predictor (`ctu_transverse_degraded_states`). Successful slab corrections are reported by
`ctu_transverse_ghost_corrected_states`; the slab and proxy counts are
disjoint subsets of `ctu_transverse_corrected_states`, so the CTU order is
full only where neither proxy nor degraded counts are non-zero. Periodic ghosts use their wrapped
real cell, so periodic boxes are corrected on both sides. Inadmissible corrected states keep the normal-only state
(`ctu_transverse_inadmissible_states`). For partial active sets the
transverse increment of a cell uses exactly the active faces, i.e. the same
faces the final update uses. Scratch: three conserved states per touched cell
plus three per ghost-shell slot and an axis-validity byte, all transient and
admitted before the advance (`ctu_transverse_workspace`).

**Hydro CFL rule (stabilization pass).** The gas timestep is selected by the
rule of the configured integration family, identically in `TimeCoordinator`
and in the accepted-step guard (`core::hydroCflStabilityRuleFor`,
`core::computeDirectionalCflTimeStep`, `core::makeHydroCflDiagnostics`):

* `muscl_hancock` (directionally unsplit, no transverse terms):
  `dt = C / sum_{active d} (|v_d| + c_s) / (L dx_d)`, `C = 0.4`. Linear
  analysis of the implemented update with first-order states gives the unsplit
  donor-cell amplification factor `g = 1 - sum_d nu_d (1 - e^{-i theta_d})`,
  which is stable iff `sum_d nu_d <= 1`; the limited, Hancock-predicted scheme
  is used with the same bound plus margin. The previous per-axis rule admitted
  `sum_d nu_d = 1.2` for 3D diagonal flow and is removed for this family.
* `ctu`: `dt = C min_{active d} L dx_d / (|v_d| + c_s)`, `C = 0.4` (per-axis
  `nu_d <= 1/2` for the simplified 3D variant).
* `L = a` for comoving coordinates (comoving widths, peculiar speeds), `1`
  otherwise. An axis is inactive only when one cell spans the whole configured
  domain along it (`core::activeHydroAxes`, rank-independent), so 1D/2D setups
  keep their 1D/2D steps (identical to before in 1D).
* Consequence: default (`muscl_hancock`) 3D runs take steps up to 3x smaller
  (isotropic sound-speed limited) than before; results and step counts are not
  reproducible against earlier binaries (intended stability fix, no config key
  restores the unsafe rule).

## 5. Riemann solvers and passive scalars (audit F3)

HLL intermediate fluxes are evaluated with lane-generic `HydroConservedState`
arithmetic, so every lane - including `metal_mass_density_comoving` and any lane
added later - is transported in the supersonic, HLL, HLLC-star and fallback
branches. For uniform `Z`, `F(rho Z) = Z F(rho)` exactly (audit case: mass 1,
metal 0.2). HLLE wave speeds are the min/max of Davis and Einfeldt (1988,
EOS-agnostic `d^2` form) estimates. HLLC uses Davis speeds and falls back to
HLLE when the contact denominator vanishes, `p* <= 0`, or a star state is
non-positive.

Extension contract for Part 2 species/dust: add the lane to
`HydroConservedState` (and its arithmetic, SoA, flux-register record and MPI
correction record), and to `HydroPrimitiveState` as a mass fraction. No
per-face heap allocation is needed and DMO/gravity paths never see the lanes.

## 6. EOS / thermodynamic authority, floors, refresh (audit F4)

* Authority (both fixed-grid and AMR paths): comoving `rho`, peculiar `u`,
  canonical specific internal energy `e = gas_cells.internal_energy_code`,
  metal mass. `E` is rebuilt as `rho e + K` at the start of each hydro advance.
  `P`, `c_s`, `T` are derived caches. **Repair:** the AMR loader and
  `primitiveForRow` used the pressure cache (or `max(P_cache, (gamma-1) rho e)`)
  as thermal authority, which ignored internal-energy changes that did not
  refresh the cache and mis-built `E` for effective-EOS cells; they now use `e`
  (the cache is only a fallback when `e` was never set).
* Closure: ideal gas (`gamma = 5/3`) or `HydroThermodynamicClosure`
  (effective multiphase ISM) evaluated at the flux epoch.
* Temperature cache: `T = f e` with `f = (gamma-1) mu m_p v_unit^2 / k_B`,
  `mu = 0.59` (the cooling module's convention), applied by the fixed-grid
  store, the AMR scatter, AMR reflux (local and remote registers), AMR
  refine/derefine writes and the source-stage refresh. **Repair:** these paths
  wrote `T = e` (code units) or never updated `T`. The AMR paths now go through
  `amr::AmrThermodynamicCachePolicy` (temperature factor + active closure), so
  effective-EOS `P`/`c_s` are also produced at the mutation itself (the former
  whole-rank post-AMR closure pass, which read the stale pressure cache, is
  removed). Focused callers that pass no policy keep the legacy values.
  `mu` is a parameter of `hydro::temperatureKelvinPerInternalEnergyCode`
  (Part 1 reference value 0.59); Part 2 supplies composition/ionization-
  dependent `mu` through the closure interface.
* Gravity half-kicks rescale `e`, `P`, `T` by the adiabatic expansion factor
  and `c_s` by its square root: exact for the ideal closure; for effective-EOS
  cells (whose `P`, `c_s` depend on density, not `e`) the cache deviates by
  O(H dt) for at most one step and is re-derived at the next hydro store.
* Floors: reconstructed-state floors (above); the conserved update applies the
  internal-energy floor `E >= K + epsilon`, measured in `floor_delta`.
* Refresh contract: `hydro::refreshGasThermodynamics()` /
  `refreshAllGasThermodynamics()`. Any stage that mutates `rho`, `e` or
  composition outside the hydro store calls
  `SimulationState::notifyGasThermodynamicMutation()` then a refresh. Wired:
  feedback deposition (touched-cell log, before mutation), star formation
  (its reported host rows), BH accretion (whole-rank refresh, only on steps
  with accretion, because the BH model exposes no touched-cell list). Both
  hydro loaders take the canonical `e` and use the pressure cache only when
  `e` was never set. `TimeCoordinator` refuses to build
  timestep criteria from stale caches (`std::logic_error`).

## 7. Dual energy (high-Mach thermal accuracy)

Selectable: `hydro_dual_energy = off` (reference: thermal energy is always
`E - K`) or `thermal_fraction`.

* **Lanes and authority.** The conservative update owns `rho`, `rho u`, `E`,
  `rho Z` in every mode. With dual energy on, `HydroConservedStateSoa` carries
  an auxiliary lane `rho e_aux` (not a flux lane) initialized each advance
  from the canonical `e` - never from `E - K`, whose cancellation error is the
  problem being solved. Since P1-55-01 every AMR ghost row (physical,
  same-level, coarse/fine, temporal, MPI, deep/edge stencil) carries the
  auxiliary value formed with the same source coverage, weights and time as
  its conserved state (see the AMR paragraph below); a ghost row whose source
  `e` is unknown (for example a temporal ghost interpolated from a pre-v25
  restored history) is NaN, its primitive state uses `E - K`, and it is
  counted once per prepared ghost in the AMR diagnostic
  `dual_energy_unresolved_ghosts`. A real cell whose lane is unknown at step
  start is counted in the solver's `dual_energy_unknown_auxiliary_cells`.
* **Primitive selection (reconstruction, Riemann, sources).** For a cell with
  `E - K <= 0` or `E - K < eta E` and a positive auxiliary value, `e`, `P` and
  `c^2` come from `rho e_aux`; otherwise from `E - K`. Effective-ISM closures
  are applied after this selection.
* **Auxiliary update (non-conservative):**
  `(rho e)^{n+1} = (rho e)^n - dt/(aV) sum_f A F_rho e_up - P^n dt/(aV) sum_f A u_f.n + dt (S_E - u.S_m)`,
  with `e_up` the reconstructed face specific energy upwind of the mass flux,
  `u_f` the mean of the face states' normal velocities, and `S_E - u.S_m` the
  thermal share of the hydro sources (cooling, effective-ISM relaxation,
  `-3HP` for expansion).
* **Final selection per cell:** `E - K` if positive and `>= eta E` (shocks
  heat through `E` and are captured; the auxiliary lane is resynchronized to
  `E - K`), else `rho e_aux` if positive (cold hypersonic flow, strong
  rarefactions), else a thermal floor (counted, `dual_energy_floor_cells`).
* **No total-energy overwrite.** The solver never modifies `E`. The runtime
  stores `e` from the selected lane value; the resulting change of the
  canonical thermal energy, `sum (rho e_sel - (E - K)) V`, is the complete
  non-conservative ledger (`HydroConservationReport::dual_energy_sync_delta`,
  `dual_energy_sync_count`) and is emitted per step. The conservative residual
  excludes nothing new.
* **Persistence.** No new persistent field: the canonical `e` is already
  restart truth, and the auxiliary lane is rebuilt from it each advance.
* **Interactions.** Cooling/ISM sources enter through `S_E`; feedback and other
  source stages mutate `e` between advances (the lane picks them up).
* **AMR reflux (explicit fallback).** Reflux corrects mass, momentum, `E` and
  metals conservatively. If the corrected `E - K` of a coarse cell is below
  `eta E` (or non-positive) the cell keeps its pre-reflux specific internal
  energy instead of the cancellation-dominated `E - K`; the thermal
  correction is ledgered in `RefluxDiagnostics::dual_energy_reflux_*` and
  reported in the `hydro.amr_production_stage` event.
* **AMR ghosts (P1-55-01).** The auxiliary lane of every AMR ghost is the
  comoving density `rho e` formed by the ghost filler with *exactly* the
  operator that forms the ghost's conserved state: single-source injection
  (same level, coarse to fine), the equal-volume mean over the complete fine
  stencil (fine to coarse), the temporal blend
  `(1 - alpha) (rho e)_start + alpha (rho e)_end` with the same `alpha` from
  the persistent temporal history (whose v25 records carry the `rho e`
  endpoints), the owner's `rho e` of a sparse MPI cell (a remote recursive
  interval is interpolated in `rho e` and divided by the interpolated `rho`),
  and the copied/mirrored cell for physical, deep-layer and edge closures.
  Local sources, the target's real cells and the step-start ghost snapshot
  carry the lane (admitted as a seventh SoA lane and one double per prepared
  ghost). A ghost touching any source with unknown `e` is stored as NaN
  (explicitly invalid, never an unrelated finite value) and is counted in
  `dual_energy_unresolved_ghosts`; the solver then uses `E - K` there.
  `AmrHydroGhostFillDiagnostics::auxiliary_internal_energy_ghosts_{transferred,invalid}`
  report the transfer. Remote conserved states use the owner's canonical `e`
  as thermal authority (same as local patches). Per-step AMR thermal
  corrections are summed in `dual_energy_thermal_correction_code`.
* **Limits.** Weak shocks whose post-shock thermal fraction stays below `eta`
  are not heated by the auxiliary branch (same as Enzo-type schemes);
  `eta = 1e-3` is the default; the pdV face velocity is first-order accurate.

## 8. Diagnostics (per hydro step)

`HydroProfileEvent`: `riemann_solve_count`, `riemann_fallback_count`,
`positivity_fallback_count`, `stencil_degradation_count`,
`ctu_transverse_{corrected,degraded,inadmissible}_states`,
`dual_energy_auxiliary_cells`, `bytes_moved`, timings, and the conservation
report (`flux/source/floor/dual_energy_sync` deltas, residual). Fixed-grid runs
emit them in the `hydro.conservation` event, AMR runs in the AMR diagnostics
event.

## 9. Periodic boundaries (supplementary repairs)

* **Fixed grid:** a periodic box face appears twice (lower ghost face of cell
  0, upper ghost face of cell n-1) and both deposit on the same two real cells;
  the runtime activated both, applying every box-boundary flux twice. Only the
  upper copy is evaluated now; the wrapped source drives activity.
* **AMR:** physical periodic ghosts were built with `source_real_cell = owner`,
  so each ghost copied the boundary cell itself and the flux was deposited back
  on that cell: every periodic box face of an AMR patch was a zero-flux wall.
  (Zero NET global boundary flux is physical for a periodic box; zero LOCAL
  flux is not.) Now a patch that spans the box along an axis wraps in-patch
  (ghost = opposite cell, flux deposited there, face evaluated once).
  A patch that touches the box boundary WITHOUT spanning the box would need
  wrapped neighbour patches (ghost sampling, coarse/fine registers,
  distributed exchange across the wrap); that topology is refused with a
  `std::runtime_error` (Part 3 scope) instead of silently closing the box.
  **Operational consequence:** periodic `cosmo_cube` gas on a decomposed root
  grid whose patches do not span the box now stops with this error; the
  previous runs with such topology had closed boxes.
* Comoving coordinates: the periodic box extent is
  `cosmology.box_size_{x,y,z}_mpc_comoving` in the patch coordinate frame,
  the same convention the PM/tree gravity and the fixed-grid layout use.
  Periodic TreePM is untouched.

## 10. Compromises stated explicitly

* Hydro sources (cooling, ISM relaxation) are evaluated on the old state and
  added to the flux increment (first-order operator coupling); the Hancock
  predictor carries no source terms.
* "Second order" applies to smooth flow away from limiter activity and with
  PLM/PPM stencils available; measured convergence is a human-run target.
* Reconstruction counters are per-instance and not thread-safe; the solver is
  serial per patch (OpenMP parallelism, if added, must use thread-owned
  reconstruction/Riemann instances).
