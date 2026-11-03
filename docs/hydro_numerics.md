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
`neighbor_plus2_cell`). Fixed-grid periodic boxes complete all stencils by
wrap (`completePeriodicReconstructionStencils`). Where the wide cells are not
available (non-periodic patch boundary, AMR patch boundary: the AMR ghost shell
is one cell wide) the side reduces to PLM, or to first order if the near cell is
also missing; every reduction is counted in `stencil_degradation_count`
(fixed-grid `hydro.conservation` event, AMR `hydro.amr` diagnostics). This is
an explicit, measured order reduction, not a silent selector change.
Multi-layer AMR ghost shells are Part 3 scope.

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
`<= 1` in 2D, `<= 1/2` in 3D (no Saltzman 12-solve corner terms). With
first-order states the scheme reproduces the exact donor-cell corner-transport
update (unit test `test_hydro_ctu`). Sides without transverse information
(open/reflective/imported-MPI ghosts) use the normal-only state and are counted
as `ctu_transverse_degraded_states`; inadmissible corrected states keep the
normal-only state (`ctu_transverse_inadmissible_states`). For partial active
sets the transverse increment of a cell uses exactly the active faces, i.e. the
same faces the final update uses. Scratch: three conserved states per touched
cell, batch-scoped, admitted before the advance (`ctu_transverse_workspace`).

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

* Authority: `rho`, `e` (via `E - K` inside hydro), `u`, metal mass. `P`, `c_s`,
  `T` are derived caches.
* Closure: ideal gas (`gamma = 5/3`) or `HydroThermodynamicClosure`
  (effective multiphase ISM) evaluated at the flux epoch.
* Floors: reconstructed-state floors (above); the conserved update applies the
  existing internal-energy floor `E >= K + epsilon`, measured in `floor_delta`.
* Refresh contract: `hydro::refreshGasThermodynamics()`
  (`include/cosmosim/hydro/gas_thermodynamic_refresh.hpp`). A stage that
  mutates the authoritative inputs outside the hydro primitive store calls
  `SimulationState::notifyGasThermodynamicMutation()` and then the refresh.
  Wired today: stellar feedback deposition in `SourceRuntime` (touched cells are
  logged before mutation with their pre-stage `e`; `T` is rescaled by
  `e_new/e_before`; the effective-ISM closure decides branch transitions).
  `TimeCoordinator` refuses to build the next timestep criteria from stale
  caches (fail-closed `std::logic_error`).

## 7. Dual energy (high-Mach thermal accuracy)

`hydro_dual_energy = thermal_fraction` enables a step-local auxiliary
internal-energy update (Bryan et al. 1995 style switch) with total-energy
synchronization:

    (rho e)_aux = (rho e)^n - dt/(aV) sum A F_rho e_upwind
                  - P^n dt/(aV) sum A u_face.n + dt (S_E - u.S_m)

If `E - K < eta E` (or `<= 0`) after the conservative update and
`(rho e)_aux > 0`, `E <- K + (rho e)_aux`. The change is reported in
`HydroConservationReport::dual_energy_sync_delta` (excluded from the residual)
and counted. Because `E` carries the synchronized value, the next step
re-derives `e` from `E` exactly: no extra persistent field, no restart schema
change. Limits: (i) the total-energy change is a measured non-conservation;
(ii) thermal energy below `ulp(E)` (`K/(rho e) >~ 1e13`) cannot be represented
by this synchronized form and would require a persistent entropy lane (a
Part 2/3 decision if the validation shows it is needed); (iii) the face
velocity in the pdV term is the arithmetic mean of the reconstructed normal
velocities (first-order accurate source).

## 8. Diagnostics (per hydro step)

`HydroProfileEvent`: `riemann_solve_count`, `riemann_fallback_count`,
`positivity_fallback_count`, `stencil_degradation_count`,
`ctu_transverse_{corrected,degraded,inadmissible}_states`,
`dual_energy_auxiliary_cells`, `bytes_moved`, timings, and the conservation
report (`flux/source/floor/dual_energy_sync` deltas, residual). Fixed-grid runs
emit them in the `hydro.conservation` event, AMR runs in the AMR diagnostics
event.

## 9. Periodic fixed-grid boundaries (supplementary repair)

A periodic box face appears twice in the fixed-grid geometry (lower ghost face
of cell 0 and upper ghost face of cell n-1) and both deposit on the same two
real cells. The fixed-grid runtime previously activated both, applying every
box-boundary flux twice. It now evaluates only the upper-side copy and treats
the wrapped source cell as the neighbor for activity
(`HydroAmrRuntimeImpl::buildActiveFaceView`; regression in `test_hydro_ctu`).

## 10. Compromises stated explicitly

* Hydro sources (cooling, ISM relaxation) are evaluated on the old state and
  added to the flux increment (first-order operator coupling); the Hancock
  predictor carries no source terms.
* "Second order" applies to smooth flow away from limiter activity and with
  PLM/PPM stencils available; measured convergence is a human-run target.
* Reconstruction counters are per-instance and not thread-safe; the solver is
  serial per patch (OpenMP parallelism, if added, must use thread-owned
  reconstruction/Riemann instances).
