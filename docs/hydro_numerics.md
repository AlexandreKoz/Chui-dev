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
`neighbor_plus2_cell`). Every geometry builder now completes stencils through
the existing one-layer ghost shell (`completeReconstructionStencils`): real
cells first, then the adjacent boundary ghost row; periodic ghosts resolve to
their wrapped real cell, so periodic boxes (fixed grid and box-spanning AMR
patches) are full order everywhere. With one ghost layer:

| face position | PLM | PPM |
|---|---|---|
| interior, >= 2 cells from a non-periodic patch boundary | full | full |
| interior face adjacent to the boundary cell | full (ghost row used) | owner side full, neighbor side PLM |
| patch-boundary face (owner real, neighbor ghost) | owner full, ghost side first order | owner PLM, ghost side first order |

Every reduction is counted in `stencil_degradation_count` (fixed-grid
`hydro.conservation` event, AMR diagnostics). Removing the remaining
reductions needs a second/third ghost layer from the AMR ghost fill (Part 3
contract: `HydroGhostCell` gains a layer index and the fill must provide
`requiredStencilHalfWidth()` layers).

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

Riemann calls: 2 per face (6 per cell in 3D, 4 in 2D). Stability: Courant
`<= 1` in 2D, `<= 1/2` in 3D (no Saltzman 12-solve corner terms); the gas CFL
candidate is `dt = 0.4 min_d a dx_d/(|v_d|+c)`, which satisfies it. With
first-order states the scheme reproduces the exact donor-cell corner-transport
update (unit test `test_hydro_ctu`).

At patch boundaries the owner side of every boundary face is corrected (the
owner is a real cell with all its faces); the ghost side has no transverse
faces and keeps the normal-only state - counted per face side as
`ctu_transverse_degraded_states` (the corrector is not silently dropped).
Periodic ghosts use their wrapped real cell, so periodic boxes are corrected
on both sides. Inadmissible corrected states keep the normal-only state
(`ctu_transverse_inadmissible_states`). For partial active sets the
transverse increment of a cell uses exactly the active faces, i.e. the same
faces the final update uses. Scratch: three conserved states per touched cell,
batch-scoped, admitted before the advance (`ctu_transverse_workspace`).

**CFL caveat (pre-existing, unchanged):** the directionally unsplit
MUSCL-Hancock family without transverse terms is linearly stable only for
`sum_d nu_d <= 1`; the per-axis CFL coefficient 0.4 allows `sum_d nu_d` up to
1.2 for 3D diagonal flow. CTU does not have this problem. A dimension-summed
CFL option for MUSCL-Hancock is recommended for the Part 3 scheduler work; it
was not changed here because it alters every run's timestep.

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
  store, the AMR scatter and the source-stage refresh. **Repair:** the AMR
  scatter wrote `T = e` (code units) and the fixed-grid store never updated `T`.
  AMR reflux and regrid still write the legacy `T = e` for the cells they touch
  until the next hydro scatter (documented limitation).
* Floors: reconstructed-state floors (above); the conserved update applies the
  internal-energy floor `E >= K + epsilon`, measured in `floor_delta`.
* Refresh contract: `hydro::refreshGasThermodynamics()` /
  `refreshAllGasThermodynamics()`. Any stage that mutates `rho`, `e` or
  composition outside the hydro store calls
  `SimulationState::notifyGasThermodynamicMutation()` then a refresh. Wired:
  feedback deposition (touched-cell log, before mutation), star formation and
  BH accretion (whole-rank refresh on steps where they removed gas, because
  those models have no touched-cell lists). `TimeCoordinator` refuses to build
  timestep criteria from stale caches (`std::logic_error`).

## 7. Dual energy (high-Mach thermal accuracy)

Selectable: `hydro_dual_energy = off` (reference: thermal energy is always
`E - K`) or `thermal_fraction`.

* **Lanes and authority.** The conservative update owns `rho`, `rho u`, `E`,
  `rho Z` in every mode. With dual energy on, `HydroConservedStateSoa` carries
  an auxiliary lane `rho e_aux` (not a flux lane) initialized each advance
  from the canonical `e` - never from `E - K`, whose cancellation error is the
  problem being solved. Physical ghosts copy their source; same-level and
  coarse/fine AMR ghosts and imported MPI ghosts are "unknown" (NaN): their
  primitive state uses `E - K` (not counted per ghost); a real cell whose lane
  is unknown at step start is counted in `dual_energy_unknown_auxiliary_cells`.
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
  source stages mutate `e` between advances (the lane picks them up); AMR
  reflux applies flux corrections to the stored state after the scatter and
  re-derives `e` through `E - K` for the corrected coarse cells (cancellation
  there is not protected: limitation).
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
