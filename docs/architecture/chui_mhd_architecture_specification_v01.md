---
title: "CHUÍ | Magnetohydrodynamics and Non-Ideal MHD"
subtitle: "Normative Mathematical, Numerical, AMR, Software and Verification Specification — CHUI-MHD-REF v1.0"
author: "CHUÍ / CosmoSim technical design"
date: "22 November 2026"
lang: en-US
---

# Contents

- [Document control and mission](#document-control-and-decision-status)
- [1–3. Scope, source basis, and frozen reference decisions](#1-mission-scope-and-normative-language)
- [4–6. Governing equations, discrete geometry and ideal solver](#4-governing-physical-model-and-units)
- [7–9. AMR, non-ideal transport and physical coupling](#7-constrained-transport-and-patch-amr-mathematics)
- [10–12. ICs, software interfaces and configuration](#10-initial-conditions-boundaries-and-magnetic-seeding)
- [13–15. HDF5, memory governance and failure behavior](#13-hdf5-data-snapshots-restart-and-provenance)
- [16–20. Validation, Pareto studies, implementation and release](#16-validation-framework-and-measurable-acceptance-gates)
- [Appendices A–C and bibliography](#appendix-a-full-reference-paramtxt-example-and-parser-rejection-rules)

---

# Document control and decision status

| Field | Frozen value |
|---|---|
| Document identity | `chui_mhd_architecture_specification_v01` |
| Version | 1.0, design baseline; 22 November 2026 |
| Status | **Normative proposed implementation contract**; neither implemented nor validated by this document |
| Applies to | Newtonian ideal MHD and single-fluid Ohmic–ambipolar–Hall MHD; cosmological cubes, zoom-ins, isolated galaxies and clusters |
| Reference profiles | `mhd_ref_ideal_v01`, `mhd_ref_nonideal_v01` |
| Existing repository | `AlexandreKoz/CHUI-Cosmology-Simulator`, `main`, inspected 22 November 2026 |
| Design authority | This specification controls new MHD code; existing hydro/AMR behavior is preserved for runs with MHD disabled |
| Change control | Numeric contracts, state ownership, cosmological normalization and output schemas change only by versioned design amendment and migration tests |

**Meaning of “complete.”** All architectural selections, defaults, ownership rules, compatibility policies, error behavior and validation gates for the *reference MHD v1* are fixed here. Where a scientific outcome is not yet measured, this document fixes the experiment and acceptance rule; it does **not** falsely declare validation complete. Advanced DG, vector-potential and implicit Hall research are explicitly separate from the v1 reference and cannot alter its outputs silently.

## 1. Mission, scope and normative language

CHUÍ's objective is not the largest collection of solvers. Its objective is the most physically faithful and scientifically predictive simulation attainable within a declared computational budget, with an equally rigorous facility for investigating new physics and numerical algorithms. The reference MHD module must therefore offer: reproducible ideal and non-ideal magnetic evolution; robust conservation; physically interpretable transport coefficients; AMR/MPI correctness; low, accountable memory use; and quantitative comparisons against independent codes.

**MUST / SHALL** denotes a non-negotiable reference acceptance condition; **SHOULD** denotes an engineering preference for which an alternative requires written evidence; **MAY** denotes a supported optional path. A claimed *certified* profile MUST satisfy every applicable gate in Section 20. A failed gate keeps the implementation in development status.

**In scope, v1:** 3D Newtonian single-fluid MHD, ideal induction, Ohmic resistivity, ambipolar diffusion, Hall term, energy-consistent transport, gas cooling/source compatibility, static and dynamic ratio-2 patch AMR, MPI ghost exchange, local AMR subcycling after its prerequisite scheduler contract is validated, cosmological expansion, initial conditions, science snapshots, exact restarts and bounded memory. All three non-ideal coefficients may be spatially variable and Hall coefficients signed.

**Explicitly excluded from REF-v1:** GRMHD, anisotropic Braginskii viscosity/conduction, radiation-MHD, cosmic-ray MHD, electron-inertia extended MHD, two-fluid plasma, shearing-periodic boxes, embedded unstructured Voronoi dynamics, adaptive switching between discontinuous Galerkin and finite volume, and automatic region-by-region switching of the physical governing equations. The runtime MUST reject requests for these under the reference profile rather than accepting inert feature switches. They are research/prospective modules, not missing decisions inside this specification.

**Success criteria:** (i) numerical equivalence and convergence for matched PDEs; (ii) stable, conservative and divergence-compatible AMR/MPI evolution; (iii) physics that does not hide artificial resistivity, heating or floors; (iv) performance evaluated as time and memory to a specified accuracy; (v) no extra MHD memory or changes in results when MHD is disabled.

## 2. Source basis, current implementation and precedence

### 2.1 Existing CHUÍ contracts inspected

- `include/cosmosim/hydro/hydro_core_solver.hpp`: `HydroConservedState`, `HydroConservedStateSoa`, `HydroPatchGeometry`, flux-register sink, update/source contexts; current state includes gas mass, momentum, thermal-plus-kinetic total energy and passive metal mass but **no magnetic energy or face flux**.
- `include/cosmosim/hydro/hydro_riemann.hpp`: HLLC default and HLLE fallback for hydrodynamics; the MHD solver must use a separate MHD Riemann contract rather than reinterpret HLLC.
- `include/cosmosim/amr/amr_framework.hpp` and `include/cosmosim/amr/amr_hydro_orchestrator.hpp`: stable patch identity, cell transfer, ordinary finite-volume flux registers, reflux and AMR execution.
- `include/cosmosim/core/time_integration.hpp`: existing eight-stage integration contract; the MHD update is owned within `kHydroUpdate`, not by inserting a ninth global stage.
- `include/cosmosim/core/simulation_state.hpp`: canonical state, stable gas-cell IDs and sidecars. Magnetic face fluxes are patch-geometric state, not particle-indexed gas sidecars.
- `docs/architecture/amr_hydro_subcycling_and_flux_registers.md`: describes locally scoped AMR subcycling, temporal ghosts and outstanding MPI-distributed subcycling limitations at the date of inspection. MHD must not assume those limitations have already disappeared.
- `docs/architecture/adr_param_txt_configuration_format.md`: **authoritative simulation configuration is `.param.txt`**, overriding earlier project literature suggesting YAML/JSON input.
- `docs/output_schema.md`: science snapshot and restart are distinct; current document states native snapshot v6 and checkpoint v24. MHD extensions require explicit new schema versions, not mutation of old schemas in place.
- CosmoSim naming rules and CHUÍ memory-governor specification: units/frames must be explicit; magnetic state must have clear ownership, bounded communication and phase-local scratch.

**Repository naming resolution:** existing implementation uses `.cpp` and public `.hpp`, despite an older document preferring `.cc/.hh`. This specification uses `.cpp/.hpp` for new files matching the live tree, and requires an explicit naming-governance approval if the broader convention is changed. All new stems remain `lower_snake`. The existing `cosmosim::...` namespace remains authoritative; no parallel `cosmo::...` tree is introduced.

### 2.2 Literature reference basis

Athena++/AthenaK establish the modern structured-AMR, finite-volume and CT reference class [R1–R3]. IDEFIX supplies well-documented Ohmic, ambipolar, Hall and RKL2 designs [R4]. AREPO moving-mesh Ohmic/ambipolar/Hall work is the independent irregular-mesh and physical-stabilization comparator [R5–R7]. Balsara and Tóth–Roe specify divergence-preserving AMR transfers [R8–R9]. Gardiner–Stone and Mignone–Del Zanna support unsplit and upwind CT [R10–R11]. RKG and new vector-potential research are comparative directions, **not dependencies of the reference** [R12–R13]. See References.

## 3. Frozen numerical architecture: CHUI-MHD-REF-v1

| Design axis | Reference decision | Why |
|---|---|---|
| Gas discretization | Conservative second-order **method-of-lines unsplit finite volume** | Robust shocks; direct compatibility with existing cell averages and reflux |
| Time integrator for ideal subsystem | Two-stage **SSPRK2 (Heun)**; one common stage-time for gas and magnetic flux | Coherent CT and gas updates; no ambiguous MUSCL-Hancock predictor reuse |
| Reconstruction | Monotone piecewise linear MC limiter for mandatory minimal reference; **PPM is the certified accuracy profile** | A stable fallback and a competitive modern high-accuracy mode |
| Riemann solver | **HLLD** for ideal MHD; **HLLE** per-face fallback | Accurate Alfvén/contact structure with robust escape path |
| Magnetic truth | **Oriented face-integrated magnetic flux** `Phi_x/y/z` in double precision | Natural CT/AMR and exact discrete topological constraint |
| Discrete induction | Upwind **UCT-HLL** edge EMFs, with consistent one-dimensional reduction | Avoids simple arithmetic edge EMFs and supports multidimensional upwinding |
| Magnetic divergence | Geometric face-flux CT, no Powell or GLM terms in reference evolution | Divergence-preserving discrete topology |
| AMR refinement | Ratio two, divergence-preserving Tóth–Roe-style face prolongation; area-flux restriction | Compatible with block AMR |
| AMR synchronization | Gas **reflux** plus magnetic **reflux-curl** using time-integrated edge EMFs | Both conservation laws maintained |
| Non-ideal coupling | Symmetric Strang composition: `D(dt/2) → H(dt/2) → I(dt) → H(dt/2) → D(dt/2)` | Reproducible second-order split with separately testable operators |
| Diffusion `D` | Ohmic+ambipolar, explicit SSPRK2 if nonstiff; **RKL2** when stiffness threshold crossed | Reference method well documented in IDEFIX |
| Hall `H` | Whistler-speed-aware **HLL induction flux**, SSPRK2 with bounded subcycles | Stable reference without undocumented artificial resistivity |
| Stiffness defaults | `cfl_ideal=0.30`, `cfl_diffusion=0.15`, `cfl_hall=0.10`, `rkl_switch_ratio=4`, `rkl_max_stages=128` | Conservative, fixed starting safety policy; must pass gates before certification |
| Thermodynamics | Ideal-gas `gamma=5/3` in certification; compatible existing thermodynamic closure only after independent equivalence tests | Prevents effective-ISM EOS from obscuring MHD validation |
| Energy | One authoritative **gas+kinetic+magnetic total energy** in MHD mode | Supports conservative Poynting flux and resistive heating |
| Comoving magnetic variable | `B_comov = a² B_phys`; physical coefficient conversion `eta_comov = eta_phys/a²` for physical diffusivity | No artificial expansion term in magnetic flux evolution |
| Main execution | CPU/OpenMP+MPI first, SoA patch kernels; CUDA/GPU backend after reference gates | Desktop-first, accelerator-amenable without two physics codes |
| Failure behavior | Fail-closed on invalid divergence, stale topology, nonphysical transport, impossible timesteps, or unbudgeted memory | No silent physical changes |
| Reference comparison | Athena++/AthenaK (ideal/AMR), IDEFIX (non-ideal), AREPO (cosmological/moving-mesh) | Independent strengths |

**No runtime mix-and-match of geometric authorities:** MHD-REF-v1 has exactly one authoritative magnetic representation (face fluxes). The derived cell-centered field is a cache or output quantity. Research solvers with a different representation must use separately versioned state and controlled conversion, never reinterpret the reference restart bytes.

## 4. Governing physical model and units

### 4.1 Physical-coordinate equations (normalized Lorentz–Heaviside form)

For mass density `rho`, physical peculiar velocity `u`, thermal pressure `p`, magnetic field `B` (with magnetic permeability absorbed), internal energy per mass `e_int`, and total energy per volume `E = rho*e_int + 0.5*rho*|u|² + 0.5*|B|²`, the homogeneous ideal part is:

$$
\partial_t\rho+\nabla\cdot(\rho\mathbf u)=0,
$$
$$
\partial_t(\rho\mathbf u)+\nabla\cdot[\rho\mathbf u\mathbf u+(p+B^2/2)\mathbf I-\mathbf B\mathbf B]=\mathbf S_m,
$$
$$
\partial_t E+\nabla\cdot[(E+p+B^2/2)\mathbf u-(\mathbf u\cdot\mathbf B)\mathbf B+\mathbf E_{ni}\times\mathbf B]=S_E,
$$
$$
\partial_t\mathbf B=-\nabla\times\mathbf E,\qquad\nabla\cdot\mathbf B=0,
$$
$$
\mathbf E=-\mathbf u\times\mathbf B+\mathbf E_{ni},\quad \mathbf J=\nabla\times\mathbf B,\quad \hat{\mathbf b}=\mathbf B/|\mathbf B|,
$$
$$
\mathbf E_{ni}=\eta_O\mathbf J+\eta_H(\mathbf J\times\hat{\mathbf b})+\eta_A[\mathbf J-(\mathbf J\cdot\hat{\mathbf b})\hat{\mathbf b}].
$$

**Important sign and normalization convention:** `eta_H` is **signed**. `eta_O` and `eta_A` are nonnegative physical magnetic diffusivities, expressed in `length²/time`. In SI, a physical magnetic field in tesla is converted to the code's normalized magnetic field consistently with `sqrt(mu0)`; in Gaussian units use the corresponding `sqrt(4*pi)` convention. One shared unit conversion module owns these transforms. No other module introduces `4π`, `mu0`, or `c` factors into numerical kernels.

For `B=0`, `b_hat` is set to zero **only inside coefficients that must vanish**, with a specified regularization `B_reg = max(|B|, 1e-30 * B_unit_code)`; coefficients are required to approach their physically admissible zero-field limit. No artificial Hall direction is selected at a null field.

For calorically perfect gas: `p=(gamma-1)*rho*e_int`. The conserved passive metallicity lane uses `rho*Z` and is advected with gas as in the preexisting hydro contract.

### 4.2 Non-ideal energy and entropy

The non-ideal energy flux is `F_E,ni = E_ni × B`. With correctly coupled magnetic/thermal energy and nonnegative transport coefficients,

$$
Q_{ni}=\mathbf J\cdot\mathbf E_{ni}=\eta_O|\mathbf J|^2+\eta_A|\mathbf J_\perp|^2\geq0,\qquad\mathbf J\cdot\mathbf E_H=0.
$$

No second independent `+Q_ni` is added to the already conservative total-energy equation. Heating is recovered from the decrease in magnetic energy and the conservative non-ideal Poynting-energy transport. Diagnostics independently calculate `Q_ni` to detect mismatches. The reference solver MUST demonstrate that Hall-only evolution does not generate net physical heat and that periodic Ohmic/ambipolar systems preserve total energy within tolerance, with magnetic loss appearing as thermal gain.

### 4.3 Cosmological variables and exact expansion convention

Define `r=a*x`, `u = dr/dt - H*r` (physical peculiar velocity), `rho_c=a³*rho_phys`, `p_c=a³*p_phys`, and

$$
\mathbf B_c=a^2\mathbf B_{phys},\qquad
\mathbf m_c=\rho_c\mathbf u,\qquad
E_c=\rho_c(e_{int}+|\mathbf u|^2/2)+\frac{|\mathbf B_c|^2}{2a}.
$$

At each RK stage with stage scale factor `a`, define `B_eff=B_c/sqrt(a)`. Standard ideal-MHD Riemann formulas using `rho_c,p_c,u,B_eff,E_c` give the comoving-volume gas+magnetic flux; the spatial Godunov divergence carries the existing factor `1/a`. The `B_c`-CT EMF is `E_c,edge=-(u×B_c)/a` for ideal induction. This preserves `B_c=const` for a homogeneous flux-frozen field expanding with the Universe.

In comoving coordinates, the induction equation is

$$
\partial_t\mathbf B_c=-\nabla_x\times\left[-\frac{\mathbf u\times\mathbf B_c}{a}+\frac{\boldsymbol{\eta}_{phys}}{a^2}\mathcal{J}(\mathbf B_c)\right],\qquad \mathcal J=\nabla_x\times\mathbf B_c,
$$

where the tensor notation denotes applying the Ohmic/Hall/ambipolar constitutive operator to `mathcal J` after evaluating all coefficients from **physical** gas and field quantities. Each physical diffusivity `eta_phys` enters the comoving induction operator as `eta_phys/a²`; its actual code numeric value is divided by `L_unit²/T_unit` for the adopted simulation base units. Hall polarity and ambipolar projection use `B_phys` directions, equivalent to `B_c` directions for `a>0`.

The comoving momentum flux magnetic terms are `(B_c²/(2a))*I - (B_c B_c)/a`; the magnetic energy is `B_c²/(2a)`. Gas and magnetic energy fluxes MUST use the **same stage value of `a`** and compatible induction/energy quadratures. The MHD expansion source in a spatially uniform, force-free state is exactly:

$$
\dot{\mathbf m}_c=-H\mathbf m_c,\qquad
\dot E_c=-H\left(2K_c+3p_c+M_c\right),\quad
K_c=\frac{|\mathbf m_c|^2}{2\rho_c},\ M_c=\frac{B_c^2}{2a}.
$$

The expanded nonuniform formulation is the conservative comoving flux system plus this expansion source, gravitational work and other registered source terms. It MUST be derived and implemented using the same cosmological timeline's `a(t), H(t)` and midpoint/stage endpoints—not independently recomputed in MHD. The existing hydro expansion/source ownership is extended, not duplicated. In **isolated modes**, `a=1,H=0`, and the equations reduce exactly to the physical-coordinate system.

*Scope note:* these are Newtonian peculiar-velocity equations; relativistic corrections, relativistic electric energy and displacement current are outside REF-v1.

## 5. Discrete degrees of freedom and topology

### 5.1 Mesh/indices and canonical magnetic state

For Cartesian patch dimensions `nx,ny,nz` and spacings `dx,dy,dz` in comoving coordinates, store **oriented magnetic face flux** (not point values):

- `Phi_x(i,j,k)` with `0≤i≤nx`, `0≤j<ny`, `0≤k<nz`; shape `(nx+1,ny,nz)`.
- `Phi_y(i,j,k)` with `0≤i<nx`, `0≤j≤ny`, `0≤k<nz`; shape `(nx,ny+1,nz)`.
- `Phi_z(i,j,k)` with `0≤i<nx`, `0≤j<ny`, `0≤k≤nz`; shape `(nx,ny,nz+1)`.

The arrays use axis-fastest storage following existing Cartesian cell indexing: `i+Nx*(j+Ny*k)` for each array's own shape. Face orientation is the positive coordinate axis. A face value is `Phi_axis=B_axis,comov*area_axis,comov`. Fields for Riemann reconstruction and science output are derived from face values; for uniform rectangular cells the cell-centered field uses the arithmetic mean of its two bordering face field values on each axis.

Define edge-integrated EMFs on the staggered edge arrays with natural shapes `(nx,ny+1,nz+1)` for `E_x`, `(nx+1,ny,nz+1)` for `E_y`, `(nx+1,ny+1,nz)` for `E_z`. At a time stage, their values are EMF line integrals `E_e*length_e`. Use right-handed edge orientations and the incidence signs dictated by Stokes's theorem. The discrete CT update is:

$$
\Phi_f^{new}=\Phi_f^{old}-\Delta t\sum_{e\in\partial f}s_{fe}\mathcal E_e.
$$

The discrete divergence in a cell is `D(Phi)=Phi_x(i+1)-Phi_x(i)+Phi_y(j+1)-Phi_y(j)+Phi_z(k+1)-Phi_z(k)` divided by its volume if dimensional divergence is requested. Because `D*C=0`, a consistent update preserves the discrete flux sum. A rank or AMR edge may NEVER use independently rounded or independently computed EMFs on two sides of the same logical interface.

### 5.2 Face/edge ownership and representation

- One authoritative active flux value exists per *unique logical leaf face*; patch-local storage may mirror shared boundary values but must name one deterministic owner (lowest stable `(refinement_level, patch_id, face_key)` tuple). Every mirrored value is synchronized from that owner at each stage.
- One authoritative time-integrated edge EMF exists per unique logical edge. Coarse-fine edges are reconciled from covered fine contributions before the coarse correction is applied.
- A face/edge key is geometric and includes level, axis, integer global lattice coordinates and patch-hierarchy epoch. Never derive identity from changing dense gas-cell row order.
- Patch migration transfers face state by stable patch ID and geometry fingerprint; receiving rank rejects inconsistent shape, orientation, duplication or epoch.
- `B_cell` is recomputable and non-authoritative. No separate persistent vector potential or GLM scalar is stored by the reference.

### 5.3 Zero-field and large dynamic range

`B=0` is legal. The reference must recover ordinary hydro to roundoff in zero-field runs, within expected algorithmic differences. Zero-field handling uses well-defined regularized unit vectors, positivity checks and HLLD→HLLE degeneracy fallback. Magnetic floors are prohibited. The code may impose gas density/thermal floors only under the preexisting audited hydro policy; floor-induced energy is measured and reported.

## 6. Spatial solver and energy recovery

### 6.1 Two-stage unsplit update

Each ideal-MHD stage reconstructs left/right primitive states at every active face, obtains a single shared normal magnetic field from the face flux, calculates HLLD (or fallback HLLE) fluxes, constructs **upwind** edge EMFs from face electric fluxes, updates cell averages and face fluxes with the same SSPRK2 stage coefficient, and recovers primitives using total-energy semantics. Reconstruction is multidimensional in the method-of-lines sense: all x/y/z numerical fluxes are evaluated from the same stage state, not directionally split in successive sweeps.

For a fixed level with conservative state `U`, magnetic flux vector `Phi`, operator `L`, and level step `dt`:

$$
Y_1=Y_n+\Delta t L(Y_n,t_n),\qquad
Y_{n+1}=\tfrac12 Y_n+\tfrac12[Y_1+\Delta t L(Y_1,t_n+\Delta t)],
$$

with `Y=(U,Phi)`, including stage-expansion sources and boundary reconciliation. SSPRK2 is mandatory for ideal REF; no Heun/MUSCL-Hancock hybrid is allowed. The old hydro solver remains unchanged for `mhd.enabled=false` and does not share magnetic-stage scratch accidentally.

### 6.2 Reconstruction, wave speeds and fallbacks

- Default minimal reference: characteristic or primitive monotone linear reconstruction with MC limiter, verified second order in smooth flow. The certified `ppm` accuracy profile uses bounded characteristic PPM, including normal-field continuity and multidimensional magnetic compatibility. Ghost depth is **three** active cells in all directions for the certified `ppm` mode; `plm_mc` may use two. A run with inadequate ghost width fails before allocation.
- Normalize Riemann magnetic inputs as `B_eff=B_comov/sqrt(a)` for cosmological runs. The normal `B_eff` equals one shared face value on both sides. Derived sound/fast speeds use the active EOS and magnetic contributions consistently.
- HLLD is the default for finite positive input states. If its candidate states violate finite density, pressure, nonnegative thermal energy, or wave-speed sanity, recompute that face flux with HLLE. If HLLE is invalid, use a first-order positivity-preserving local flux, halve the candidate timestep and retry the full stage. Retry at most **eight** times; then fail with a diagnostic and preserve the last accepted synchronized state.
- Correction/floor triggers are counted by patch, level, rank and run. A reference-validation run with unexplained fallback rates above **0.1% of faces** or floor-added energy exceeding **1e-8** of the initial absolute total-energy scale fails the fidelity gate even if it finishes.
- Fluxes for gas, metals and magnetic energy use a consistent conservative operator; no independently reconstructed B authority is introduced.

### 6.3 Thermal-energy recovery

For a cell in comoving MHD mode:

$$
U_{int,c}=E_c-\frac{|\mathbf m_c|^2}{2\rho_c}-\frac{|\mathbf B_c|^2}{2a},\qquad p_c=(\gamma-1)U_{int,c}.
$$

Cell-centered `B_c` in this formula is derived from the authoritative face fluxes, never separately evolved. A positivity-aware face-flux correction and a timestep retry are attempted **before** invoking the existing audited thermal floor. A dual-energy/entropy tracer MAY be derived or introduced only under a separate versioned contract after tests; it is not an untracked energy authority in v1. At regrid, a change in the fine-level magnetic-energy sum may not silently alter aggregate total energy; coarse-to-fine thermal allocation is adjusted conservatively when positive, otherwise prolongation order is reduced or refinement is rejected.

## 7. Constrained transport and patch-AMR mathematics

### 7.1 Topological invariant

Let `C` denote the edge-to-face incidence matrix and `D` the face-to-cell incidence matrix. The global leaf complex MUST obey `D*C=0`, including orientation, periodic wrapping and coarse-fine interfaces. The reference numerical divergence uses **flux imbalance** rather than the divergence of a reconstructed cell-centered interpolation.

### 7.2 Refine and coarsen (ratio 2)

On refine, split each parent boundary face flux into its four child boundary faces with a limited second-order face reconstruction whose sum equals the parent flux exactly. Determine new internal child faces using the **Tóth–Roe compatible divergence-and-curl construction**, so each fine-cell flux imbalance equals the volume-split parent imbalance (zero in accepted states). Enforce exact parent/child flux sums. If the limited reconstruction yields invalid thermodynamic states after energy-consistent cell transfer, fall back to piecewise-constant magnetic split; if still invalid reject the regrid transaction before committing.

On coarsen, restrict face magnetic flux by oriented **sum**, not mean: `Phi_parent=sum(Phi_child_covering_parent)`. Restrict gas conserved quantities by volume integration. Reconcile total-energy and magnetic-energy budgets, applying no unreported heat. After successful coarsening the discarded children are deleted only after memory-governor transaction admission and successful global rank agreement.

### 7.3 Reflux-curl

At every coarse-fine synchronization interval `Δt_coarse`, accumulate coarse and fine **time-integrated** edge EMFs `I_e = ∫ E_edge dt` in the same orientation. Coarse edges covered by fine edges receive the line-integral sum of appropriately oriented fine contributions. Apply the difference as a coarse-face curl correction (`Phi_coarse ← Phi_coarse - C(I_fine-I_coarse)`) while ordinary gas reflux corrects cell-centered conserved quantities. Energy consistency after reflux-curl is checked from authoritative total energy and magnetic reconstruction; correction is accounted for through conservative reflux and thermal recovery, never a hidden global magnetic-energy rescale.

All contributions must be bound to `(patch_id, level, geometry_generation, interval, edge_key)` and rejected if incomplete, out of order, stale or inconsistent with the current hierarchy. The deterministic reduction order is global key, rank and substep ordinal (fixed order in deterministic mode).

### 7.4 Temporal boundary and AMR subcycling

Fine-level time advancement uses ratio-two recursive subcycling only if the **global** AMR level scheduler and MPI remote temporal ghost interpolation are implemented and pass equivalence tests. Because the inspected CHUÍ documentation explicitly states distributed subcycling is incomplete, the **initial released MHD reference configuration uses synchronized AMR level sweeps (`amr.mhd_subcycling=false`) for any multi-rank run**. The fully specified distributed subcycling variant is mandatory for the `mhd_ref_amr_subcycled` certification tier, not silently enabled by local-only logic.

For that tier, coarse boundary data are stored at accepted beginning/end times, fine ghost states interpolate conserved gas variables in time, magnetic boundary fluxes interpolate consistently within the divergence-preserving face complex, and the time-integrated coarse EMF history is supplied to fine boundary updates. Interpolation requests outside a strictly known interval are rejected. No refinement, derefinement or patch migration is permitted with unresolved EMF registers or active temporal boundary histories. Coarse and fine stages must have matching physical timestamps and scale factor histories.

### 7.5 MPI and physical boundaries

- **Periodic:** map corresponding face fluxes and tangential edge EMFs with exact orientation; global net magnetic flux is supported (it need not arise from a periodic vector potential).
- **Open/outflow:** use zero-gradient extrapolation for gas/tangential magnetic state, enforce one authoritative normal face field and compatible boundary EMFs, and account for physical boundary fluxes in conservation reports.
- **Reflecting perfectly conducting wall:** set wall normal gas velocity to zero, tangential electric field at the wall to zero, and preserve wall-normal magnetic flux; construct ghost states to satisfy those conditions rather than relying on inconsistent component-parity recipes.
- **Imported MPI ghost:** immutable for its stage, timestamped, ownership-checked, chunked and invalidated by migration/regrid/epoch change.
- **Shearing-periodic:** explicitly **rejected in REF-v1**.

Each unique interface face and edge obtains one ordered owner value; no halo exchange may depend on retained `O(N_local × stencil)` per-contribution messages. Collective error agreement occurs before any MPI rank applies irreversible topology or reflux changes.

## 8. Non-ideal induction and transport operators

### 8.1 Coefficient provider contract

`MhdTransportCoefficients` consists of `eta_ohmic_phys`, signed `eta_hall_phys`, `eta_ambipolar_phys`, a validity mask, coefficient provenance ID and a maximum physically valid local timestep hint; all coefficients are physical magnetic diffusivities and converted to comoving code units exactly once by the source/property adapter.

**Supported providers and mandatory behavior:**

1. `constant`: three finite values in physical diffusivity units, including a signed Hall value; for analytic tests and controlled numerical comparison. It is *not* automatically a claim of realistic galactic microphysics.
2. `table`: versioned HDF5 table with axes `(log_rho_phys_cgs,log_temperature_k,log_b_phys_gauss,composition_id,grain_model_id)` and datasets `EtaOhmicCgs`, `EtaHallCgs`, `EtaAmbipolarCgs`, with an explicit `BZeroLimit` slice for zero-field queries; log-space monotone interpolation for strictly positive Ohmic/ambipolar magnitudes, sign-aware linear interpolation for Hall with an exact zero crossing. `composition_id` and `grain_model_id` are categorical exact-match indices and are never interpolated across different physical models. All continuous axes must be sorted, unique and unit-declared. For exactly `B=0`, use the declared `BZeroLimit` dataset rather than a logarithmic-field extrapolation. Outside table coverage is **fatal** in certified physics runs; no silent extrapolation, clamping or coefficient replacement. The table content SHA-256 is part of run provenance and restart admission.
3. `chemistry`: receives electron, ion, neutral and charged-grain number densities, collision rates and field strength from an explicitly validated chemical/charging provider and derives conductivities and diffusivities. If the live network lacks the required charge carriers or collision metadata, the run **fails at configuration admission** rather than substituting an arbitrary ionization fraction. A validated table provider is the production alternative when live detailed charging is too expensive.

The reference implementation is considered **numerically certified** once `constant` and `table` pass the full non-ideal test suite. `chemistry` is a *separate physics certification*: it requires charge neutrality, thermochemical and grain benchmarks, and comparisons with independently computed conductivity tables. No publication using `chemistry` can claim its certification before those gates pass.

`eta_O >= 0`, `eta_A >= 0`, and all coefficients finite are mandatory. Negative Hall values are valid. `B=0` requires a well-defined limiting transport law; an undefined charge-state coefficient is rejected, not silently set to zero.

### 8.2 Magnetic current and edge EMFs

Compute current from the face magnetic fluxes using the compatible discrete curl and required dual-mesh metric factors. The edge-current stencil is second-order centered in smooth regions, uses the same stage-time magnetic state on both sides, and preserves sign/orientation. Construct edge non-ideal contributions:

`E_ohmic = eta_O * J`, `E_hall = eta_H * (J × b_hat)`, `E_ambi = eta_A * (J - (J·b_hat)*b_hat)`.

For ideal and non-ideal EMFs, the **edge authority is one combined time-integrated value**, even though the integration is operator-split. The sum of substep edge integrals is what AMR/MPI synchronizes. The energy operator uses matching face/edge averaging to construct the non-ideal Poynting flux and the heating audit.

### 8.3 Ohmic and ambipolar update `D`

`D` evolves magnetic face fluxes and the conservative total-energy equation with Ohmic+ambipolar EMFs and Poynting-energy fluxes. The coefficients are reevaluated at each explicit/RKL stage using the current magnetic field and physical thermochemical properties; they are never frozen for an entire dynamical step if a dependency changes.

Use two-stage SSPRK2 for `R = dt_requested/dt_diff_stable ≤ 4`; use second-order RKL2 with the least integer stage count whose stability interval contains `R` when `R>4`; cap at **128 stages per cycle**, and partition a longer requested interval into as many consecutively audited RKL2 cycles as needed. The stage count is deterministic from stage-local coefficient bounds and is included in diagnostics. No STS step may exceed its formally derived stability interval, and any positivity/energy failure retries with a smaller interval rather than raising artificial coefficients. Nonlinear RKL2 is certified separately at variable coefficients, not inferred solely from constant-coefficient tests.

The safe explicit diffusivity candidate for the diagonal Cartesian metric is

$$
\Delta t_D \leq \frac{0.15}{\max_{cell}(\eta_O+\eta_A)\,\big(\Delta x^{-2}+\Delta y^{-2}+\Delta z^{-2}\big)}.
$$

Here diffusivities and lengths are in the same *comoving code* units. This is deliberately conservative; the actual implementation also evaluates the relevant anisotropic operator spectral bounds. Zero coefficients impose no diffusion constraint.

### 8.4 Hall update `H`

`H` uses an IDEFIX-inspired whistler-speed-expanded **HLL induction flux** (not the ideal HLLD gas Riemann solver) and CT edge EMFs. The gas mass/momentum flux is zero during a pure Hall split, and the conservative total-energy Poynting-flux contribution is updated consistently. This separation prevents falsely applying a whistler-augmented flux to gas variables that are not transported by the Hall operator.

The reference Hall integrator is explicit SSPRK2 with bounded subcycling. Set the safety candidate

$$
\Delta t_H \leq \frac{0.10}{\max_{cell}|\eta_H|\,\big(\Delta x^{-2}+\Delta y^{-2}+\Delta z^{-2}\big)};
$$

and independently ensure the HLL Hall wave fan bounds resolved whistler characteristics at the grid Nyquist scale. If they disagree, the **smaller** step wins. Recompute `eta_H`, `B` and Hall wave fan every Hall substage. The Hall path must preserve the expected odd-in-polarity behavior and not add artificial Ohmic diffusion. If Hall subcycling is prohibitively small, the reference reports the predicted cost and may terminate with `MHD_STIFFNESS_BUDGET_EXCEEDED`; it **must not** silently switch to a non-certified implicit algorithm or modify `eta_O`.

### 8.5 Symmetric coupling, source order and timestep hierarchy

For one accepted level-local MHD interval `dt`, execute:

`D(dt/2) → H(dt/2) → I(dt) → H(dt/2) → D(dt/2)`.

`I` is the ideal SSPRK2 MHD update, including cell-centered conservative gas fluxes, magnetic CT, and cosmological stage sources. Each suboperator carries its own matching total-energy and edge-integral update. Cooling, ionization and feedback remain owned by CHUÍ's existing `kSourceTerms` scheduler stage. Where sources modify thermochemistry or momentum during an MHD interval, the orchestrator MUST compose their two half-interval updates as `S(dt/2) → [D/H/I/H/D](dt) → S(dt/2)` at legal source-ownership boundaries (without renumbering the eight global stage IDs). If the current workflow cannot supply both half-source callbacks, combined MHD+those sources fails preflight until that callback contract exists; when their values determine non-ideal coefficients, a **read-only consistent source snapshot at each substage** is required. Source-stage updates are applied only at legal synchronization points and never concurrently write a cell that an MHD substage reads. An explicit predictor at a substage may derive thermochemistry without mutating chemistry authority.

The operator coupling is targeted at global second-order accuracy; its convergence must be demonstrated for field-dependent coefficients and cosmological expansion, not assumed from the formal Strang composition alone. The global gas-cell scheduler MUST include the minimum of ideal fast-magnetosonic, Hall and diffusion candidate times, after eligible RKL acceleration and Hall subcycle cost bounds. The reference does not integrate any magnetic step across a forbidden gravity/AMR KDK boundary.

## 9. EOS, gravity, cooling and star/AGN feedback interfaces

- **Equation of state:** baseline `gamma=5/3` perfect gas; pressure from conserved thermal energy. Existing effective-multiphase ISM closures can be used only if their pressure, sound-speed derivative, internal-energy target and magnetic energy recovery are all evaluated consistently, and after a dedicated MHD-EOS verification gate. A reference profile with this closure disabled remains independently reproducible.
- **Gravity:** TreePM gravitational acceleration acts on mass and momentum exactly once under existing source ownership. MHD does not add a second gravitational potential, gravity solve or independent gas acceleration. Magnetic pressure/tension come exclusively through the MHD fluxes. Energy gravitational work uses the same source scheme as hydro.
- **Cooling and radiation:** cooling acts on thermal energy, never directly on total energy without recalculating the thermal/magnetic split. Cooling and chemistry property sampling must use physical `rho,T,B` at correct `a` and source timestamp. If an externally imposed UV/cosmic-ray ionization parameter influences `eta`, that parameter enters provenance.
- **Star formation and feedback:** source modules can convert gas into stars or inject mass/momentum/energy only at the legal source mutation boundary; any mass/energy removed from gas has an explicit budget. Magnetic fields remain face-owned continuum state; no star particle acquires a secret duplicate `B` authority. Changes to gas topology trigger a synchronized face-flux reconciliation before any subsequent MHD step.
- **Physical applicability:** the single-fluid approximation is formally invalid where charged/neutral drift timescales or charge-carrier magnetization invalidate the assumed closure. The provider reports `closure_valid`; such cells cannot be run as certified non-ideal physics. The solver may refuse the run or require an explicit different physics profile; no automatic arbitrary Hall suppression.

## 10. Initial conditions, boundaries and magnetic seeding

Accepted initial magnetic data are oriented face fluxes satisfying the discrete divergence norm in Section 20. For analytic initialization, construct `Phi = Phi_uniform + C*A_edge`: a uniform net-flux part plus the discrete curl of an edge vector potential, so periodic domains may represent nonzero mean magnetic flux. The converter records the `seed_model`, amplitude/unit, field orientation, random seed where relevant, code version and initial divergence diagnostic. Existing GADGET/AREPO cell-centered `MagneticField` datasets are **not** accepted as exact face-flux authority without an explicit, offline, versioned compatible conversion that satisfies the same numerical acceptance gates; no silent per-cell averaging or hidden global divergence projection is permitted.

Initial-condition presets supported by v1: uniform field, discrete vector-potential loop, force-free/Alfvénic wave, magnetized shock, divergence-free stochastic seed generated in Fourier space with Hermitian symmetry and projected perpendicular to `k`, and imported certified face-flux HDF5. Cosmological initial `B_phys = B_comov/a_initial²`; a uniform field with `B_comov` constant passes the homogeneous expansion check. Seed generation MUST be deterministic with seed and algorithm recorded.

## 11. Software module and API specification

### 11.1 Ownership decisions

- `SimulationState` owns patch-attached authoritative **persistent** magnetic face fluxes through a versioned MHD state component; gas-cell sidecars continue to own gas thermochemistry and identity. MHD state is not packed into `ParticleSoa`.
- `MhdPatchState` binds `patch_id`, `geometry_generation`, `refinement_level`, shapes and three face-flux arrays. It owns no duplicate cell-centered `B` lane.
- `MhdStageWorkspace` owns patch-local primitive/reconstruction arrays, edge EMFs, currents, temporary coefficients, RK stage buffers and temporary energy fluxes, all governed phase scratch.
- `MhdAmrRegister` owns only **in-flight** integrated edge EMFs plus interval/level identity; it is released immediately after accepted reconciliation and cannot be published in a science snapshot as authoritative field.
- `MhdCoefficientProvider` reads a time-consistent cell/source view and returns physical coefficient values; it cannot mutate gas, field, chemistry or scheduler state.
- `MhdSourceAudit` accumulates physically differentiated flux/source/floor/retry contributions per stage, patch and rank.

### 11.2 Required public C++ interfaces (semantic contract)

```cpp
namespace cosmosim::mhd {
struct MhdFaceFluxView;       // read-only/read-write x/y/z face flux spans + geometry epoch
struct MhdStageContext;       // dt, t, a, H, level, patch, stage, source epoch
struct MhdConservationBudget; // mass, momentum, E_gas, E_mag, E_total, fluxes, sources
struct MhdTransportCoefficients;
class MhdCoefficientProvider; // no mutation of physical authoritative state
class MhdRiemann;             // HLLD and HLLE, patch dispatch outside face loop
class MhdEmfBuilder;          // UCT-HLL and Hall-HLL edge EMFs
class MhdStepOperator;        // ideal or non-ideal accepted substep
class MhdAmrSynchronizer;     // prolong, restrict, reflux-curl, epoch validation
class MhdDiagnosticSink;
}
```

Concrete class layout MAY be optimized without breaking this semantic contract. Values are passed as span/SoA views. Algorithm selection occurs before patch sweeps; no virtual call or dynamic allocation is allowed per face/edge. All public functions use `camelCase`, classes `PascalCase`, locals `lower_snake`, private fields `m_lower_snake`, explicit `_comov/_phys/_code/_cgs` where ambiguity exists.

### 11.3 File/module layout (mandatory stems)

```text
include/cosmosim/mhd/mhd_state.hpp
include/cosmosim/mhd/mhd_geometry.hpp
include/cosmosim/mhd/mhd_riemann.hpp
include/cosmosim/mhd/mhd_induction.hpp
include/cosmosim/mhd/mhd_transport.hpp
include/cosmosim/mhd/mhd_amr.hpp
include/cosmosim/mhd/mhd_validation.hpp
src/mhd/mhd_state.cpp
src/mhd/mhd_geometry.cpp
src/mhd/mhd_riemann.cpp
src/mhd/mhd_induction.cpp
src/mhd/mhd_transport.cpp
src/mhd/mhd_amr.cpp
src/mhd/mhd_validation.cpp
tests/unit/mhd/test_mhd_discrete_complex.cpp
tests/integration/test_mhd_ideal_reference.cpp
tests/integration/test_mhd_nonideal_reference.cpp
tests/integration/test_mhd_amr_mpi_restarts.cpp
configs/cosmo_cube/cube_mhd_ref_v01.param.txt
configs/isolated_galaxy/disk_mhd_ref_v01.param.txt
docs/architecture/chui_mhd_architecture_specification_v01.md
```

These names are normative artifacts to create. Additional implementation-local source files MAY be added as required, but no parallel monolithic fluid engine or duplicative persistent magnetic authority may be introduced.

### 11.4 Integration into the current runtime

`kHydroUpdate` remains one top-level stage. Its handler chooses hydro or coupled MHD by validated run config and then executes the five-part MHD split internally. The `StageContract` required/mutated state masks explicitly include gas, mesh-face and ghost domains; existing `StageDataDomain` requires a backward-compatible new `kMagneticFaces` bit and schema bump, never aliasing `kMeshCells` to hide write permissions. Likewise, add a semantic `kMhd` subsystem owner while keeping existing stage numbering stable for checkpoint compatibility. Any update to stage enums requires serialization migration tests. Existing source/analysis/output stages retain their order.

A global AMR regrid transaction must present consistent gas and magnetic patch mappings, magnetic flux validation and memory admission before commit. The first implementation executes synchronized MPI AMR when MHD is enabled; attempted unvalidated distributed MHD level subcycling fails at preflight.

## 12. Configuration: authoritative `.param.txt` profile

The following is a *normative proposed config stanza*, to be integrated into the existing typed `.param.txt` parser with schema-validated keys. It is not YAML/JSON and is not claimed to parse in current CHUÍ binaries.

```ini
[mhd]
enabled = true
reference_profile = mhd_ref_nonideal_v01
magnetic_representation = face_flux_ct
spatial_method = finite_volume_unsplit
reconstruction = ppm
riemann_solver = hlld
riemann_fallback = hlle
emf_scheme = uct_hll
time_integrator = ssprk2
divergence_policy = strict_ct
precision = float64
cfl_ideal = 0.30
cfl_diffusion = 0.15
cfl_hall = 0.10
max_step_retries = 8
fallback_fraction_limit = 0.001
floor_energy_fraction_limit = 1.0e-8
deterministic = true

[mhd_nonideal]
enabled = true
ohmic_enabled = false
ambipolar_enabled = true
hall_enabled = true
integration = strang_d_h_i_h_d
diffusion_integrator = auto_rkl2
hall_integrator = explicit_hll_ssprk2
rkl_switch_ratio = 4
rkl_max_stages = 128
coefficient_provider = table
coefficient_table_path = tables/mhd_transport_reference_v01.hdf5

[mhd_amr]
refinement_ratio = 2
magnetic_transfer = toth_roe_compatible
emf_reflux = curl
distributed_subcycling = false
checkpoint_sync_only = true

[mhd_io]
science_b_field = derived_cell_center
checkpoint_face_flux = true
initial_field_format = face_flux
input_projection = reject
output_divergence_diagnostics = true
```

The above is the single normative example of the new key families; the authoritative existing parser must accept sectioned typed key/value syntax before a configuration is runnable. `ohmic_enabled` is the accepted spelling; any unknown key fails preflight. `coefficient_table_path` must name a real, validated HDF5 table whose SHA-256 is persisted. The example file name does not imply that such a table is already distributed.
| Config group | Canonical keys and defaults | Validation |
|---|---|---|
| `[mhd]` | `enabled=false`, `reference_profile=mhd_ref_ideal_v01`, `magnetic_representation=face_flux_ct`, `spatial_method=finite_volume_unsplit`, `reconstruction=ppm`, `riemann_solver=hlld`, `riemann_fallback=hlle`, `emf_scheme=uct_hll`, `time_integrator=ssprk2`, `precision=float64` | Reject unsupported enum; no implicit fallback |
| `[mhd]` stability | `cfl_ideal=0.30`, `cfl_diffusion=0.15`, `cfl_hall=0.10`, `max_step_retries=8`, `fallback_fraction_limit=0.001`, `floor_energy_fraction_limit=1e-8` | Positives and ranges validated |
| `[mhd_nonideal]` | `enabled=false`, `ohmic_enabled=false`, `ambipolar_enabled=false`, `hall_enabled=false`, `integration=strang_d_h_i_h_d`, `diffusion_integrator=auto_rkl2`, `hall_integrator=explicit_hll_ssprk2`, `rkl_switch_ratio=4`, `rkl_max_stages=128`, `coefficient_provider=constant` | Nonideal requires MHD and at least one enabled term |
| `[mhd_amr]` | `refinement_ratio=2`, `magnetic_transfer=toth_roe_compatible`, `emf_reflux=curl`, `distributed_subcycling=false`, `checkpoint_sync_only=true` | Other ratios/missing sync rejected |
| `[mhd_io]` | `science_b_field=derived_cell_center`, `checkpoint_face_flux=true`, `initial_field_format=face_flux`, `input_projection=reject`, `output_divergence_diagnostics=true` | Refuse noncertified IC conversion |

`mhd_ref_ideal_v01` requires nonideal disabled; `mhd_ref_nonideal_v01` requires nonideal enabled and the requested physical coefficients. The live parser MUST reject an inconsistent combination and must include the normalized parameters in provenance and restart equivalence checks. The reference has no `automatic_new_physics` toggle.

## 13. HDF5 data, snapshots, restart and provenance

### 13.1 Science output — new schema version

The current science schema v6 MUST remain readable unchanged. MHD science output creates **`chui_science_snapshot_v7`**, adding derived `PartType0/MagneticField` as a three-component **`B_comov` in normalized CHUÍ code units** dataset with `Frame=Comoving` and `MagneticUnit=Code` attributes; conversion to external dialect fields is performed only by an explicitly named dialect converter, optional `PartType0/MagneticPressure`, `PartType0/PlasmaBeta`, and a new `/MHD/Diagnostics` collection for flux-based divergence and magnetic-energy totals. For AMR scientific analysis, include `/MHD/FaceFlux` ragged arrays and patch metadata **only when** the configured analysis product requests topology-resolved magnetic field reconstruction. Standard science snapshots remain independently interpretable but are **not** exact restart authority.

Provenance attributes include `MhdSchemeVersion=1`, `MagneticUnits`, `MagneticFrame`, `InductionOperator`, `EmfScheme`, `Reconstruction`, `NonidealIntegrator`, coefficients provider/table SHA, dimensional diffusivity convention, `DivergenceMaxNorm`, floor/retry rates, and cosmological scale factor. A new science schema must pass existing snapshot completion-marker publication semantics.

### 13.2 Exact restart — new checkpoint schema

The current checkpoint v24 remains readable and explicitly represents an **MHD-disabled** state. MHD execution writes **checkpoint v25**, with authoritative HDF5:

```text
/state/mhd/patch_id                 uint64[P]
/state/mhd/geometry_generation     uint64[P]
/state/mhd/refinement_level        uint8[P]
/state/mhd/face_dims_xyz           uint32[P,3]
/state/mhd/face_offset_x           uint64[P+1]
/state/mhd/face_offset_y           uint64[P+1]
/state/mhd/face_offset_z           uint64[P+1]
/state/mhd/face_flux_x_comov       float64[sum_fx]
/state/mhd/face_flux_y_comov       float64[sum_fy]
/state/mhd/face_flux_z_comov       float64[sum_fz]
/state/mhd/reference_contract      string
/state/mhd/coefficient_table_sha256 string
/state/mhd/divergence_max_norm     float64
```

MHD checkpoint writes occur **only at globally agreed MHD-safe synchronization boundaries**: no active RK/STS/Hall partial stage, no unresolved gas/EMF flux registers, no active coarse temporal history, and no in-flight patch migration. A checkpoint request at an unsafe point is **deferred** until the next legal boundary; it never serializes an inconsistent partial magnetic state. This is a deliberate v1 decision; it avoids claiming exact mid-stage restart while preserving existing hydro pending-register behavior for non-MHD runs.

On restart, validate schema version, normalized config hash, source identity/provenance, topology/patch IDs, face shapes/offset monotonicity, unit/scale-factor consistency, divergence bound, provider table hash, and rank compatibility. Reference v1 requires the **same MPI rank topology** for exact restart unless an independently certified redistribution migration exists; unsupported rank-count changes fail before loading into authoritative state. Science snapshots must not be used as restart files.

### 13.3 Output transactionality

Use existing CHUÍ `.partial`/completion marker discipline, including collective agreement before publishing successful MPI HDF5. Streaming/chunked writing is mandatory: no full additional magnetic-state copy is allocated merely for serialization. Checksum and bounded readback follow the current schema's accepted verification rules. Schema migration is additive for non-MHD runs and explicit for MHD-enabled runs; v24 does not pretend to contain missing B data.

## 14. Memory governance, HPC and GPUs

### 14.1 Storage cost model

For one patch, persistent face flux capacity excluding ghost widths is

$$
M_{\Phi}=8[(n_x+1)n_yn_z+n_x(n_y+1)n_z+n_xn_y(n_z+1)]\ \mathrm{bytes}.
$$

For a uniform 512³ grid, the asymptotic field-state cost is approximately **3 GiB** globally in float64, excluding boundary faces, ghost layers, migration coexistence and scratch. Each active patch additionally uses bounded edge-EMF arrays, current/transport coefficient evaluation and RK scratch. No `B_cell` triplet, `J` triplet, gradients or `eta` triplet may be permanently allocated across the whole gas population unless a separate memory-review proves its time-to-accuracy benefit.

### 14.2 Memory classes and admission

- `CanonicalPersistent`: face fluxes, patch geometric keys and safe-boundary committed MHD metadata.
- `PhaseResident`: patch primitives and any repeatedly read face reconstruction cache for the active MHD stage.
- `ScratchArena`: temporary edge EMFs, currents, diffusion/Hall stage states, flux buffers, local sorting and deterministic reductions.
- `Communication`: bounded halo slabs and cross-rank EMF/reflux messages; must have a configured worst-case budget independent of global gas count.
- `Diagnostic`: optional derived B/current/spectra allocated only when scheduled safely.

The memory governor must reserve both old+new patch-face data during regridding and both send/receive buffers during migration. If a reservation would exceed the `M_safe` ceiling, reduce patch-batch size or serialize optional tasks. It cannot silently downgrade magnetic precision or disable a physical operator. Failure is deterministic with subsystem, bytes requested, committed/high-water and safe limit. 48-GiB workstation engineering defaults inherit the memory-spec 42-GiB safe ceiling and 36-GiB nominal target (hardware-dependent, not universal values).

### 14.3 Execution kernels and scaling

The CPU reference uses contiguous SoA patch sweeps and OpenMP over owned faces/edges/patches, with explicit anti-race ownership of CT edges. GPU is a backend for the same numerical operator, not a separate scientific implementation. MPI domain balancing includes gas-cell workload, AMR patch count, MHD faces/edges, diffusive stage count, memory pressure and communication estimates. Profile `bytes moved`, `per-cell stages`, `edge EMF updates`, `HLLD fallback`, `Hall subcycles`, `MPI critical path`, `rank-max RSS` and GPU transfer/scratch. Never optimize cells/sec at the expense of scientific time-to-accuracy.

## 15. Determinism, observability and failures

`mhd.deterministic=true` orders rank-local reductions and shared face/edge updates by stable geometric keys, uses fixed-stage scheduling for matched deterministic configurations, and records compiler arithmetic flags. Cross-platform bitwise equivalence is *not* promised; toleranced scientific equivalence is the default external contract. All rejected steps are atomic: a failed RK, Hall or RKL cycle does not partially commit authoritative state.

At each accepted level interval record: gas+magnetic total energy, kinetic/internal/magnetic split, flux boundary balance, gravitational and cooling source work, dissipation integrals, maximum and RMS normalized divergence, number of positivity/fallback occurrences, Hall/diffusion minimum timestep, subcycle/stage counts, active patch count, peak reserved and observed memory, and max/mean rank imbalance. Logging must distinguish expected boundary energy transport from true numerical nonconservation.

Fatal diagnostic codes: `MHD_INVALID_FIELD`, `MHD_DIVERGENCE_CONTRACT`, `MHD_STALE_GEOMETRY`, `MHD_BAD_COEFFICIENT`, `MHD_UNSUPPORTED_PHYSICS`, `MHD_STIFFNESS_BUDGET_EXCEEDED`, `MHD_NONPOSITIVE_STATE`, `MHD_AMR_REGISTER_INCOMPLETE`, `MHD_RESTART_SCHEMA_MISMATCH`, `MHD_MEMORY_BUDGET_EXCEEDED`, `MHD_MPI_EPOCH_MISMATCH`. Each failure preserves last synchronized state and emits problem identifiers, face/edge/cell keys, physical time/scale factor and relevant measured values.

## 16. Validation framework and measurable acceptance gates

All gates are **predeclared engineering targets**, not claims that current CHUÍ runs have passed. Every gate runs in serial and parallel where applicable, on both ideal/minimal and `ppm` reference when relevant. Benchmark code, initial conditions and reference tolerances are versioned and never changed to improve a result retroactively.

### 16.1 General deterministic gates

| ID | Test | Acceptance rule |
|---|---|---|
| G01 | Zero-field hydro limit | Identically zero initial B remains zero to roundoff; MHD and legacy hydro satisfy the same smooth analytic hydro benchmark with observed order ≥1.8 in comparable settings; MHD-disabled legacy regression remains unchanged |
| G02 | Discrete complex | Numerically verify `D*C=0` on periodic, open, reflective, nested AMR and MPI interface topologies; maximum normalized flux residual ≤1e-12 in small linear-operator test |
| G03 | Long-run divergence | `max_cells |Σoriented Phi| / max(Σ|Phi|, B_ref*A_ref, eps_scale) ≤ 1e-11` after 1,000 accepted steps, including dynamic AMR and MPI; record normalizer |
| G04 | Closed periodic conservation | Without physical source work or floors, relative volume-integrated mass and total energy budget residual ≤1e-9; momentum similarly for zero-source periodic problems |
| G05 | Floors/fallback | No unexplained nonfinite primitive states, total floor energy ≤1e-8 of initial absolute total-energy scale and fallback fraction ≤0.001 on certified smooth/wave problems |
| G06 | MPI rank equivalence | 1-, 2-, 4-rank runs satisfy all analytic tolerances and agree on integral diagnostics within 1e-8 relative; same-topology restart within 1e-10 on key science fields |
| G07 | Regrid transaction | Ratio-2 refine/restrict cycles preserve parent boundary flux exactly to floating-point tolerance; total energy budget within 1e-9; all patches remain positively admissible |
| G08 | I/O | v25 MHD checkpoint exact-safe roundtrip; v24 legacy non-MHD read; v7 science output and complete marker; invalid table SHA and missing field datasets rejected |
| G09 | Memory boundedness | Scratch/communication plateau after repeated steady iterations; peak within governor hard ceiling; disabled MHD adds no resident field/scratch arrays |

`eps_scale` in G03 is `1e-30 * B_unit_code * A_ref`; `B_ref` is maximum nonzero face-derived field on the current run. This definition prevents reporting artificial zero-divergence excellence from division by zero.

### 16.2 Ideal and AMR numerical tests

| ID | Benchmark | Acceptance rule |
|---|---|---|
| I01 | 1D and oblique circularly polarized Alfvén waves | L1 error on 32/64/128 cells per wavelength exhibits measured order ≥1.8 on the smooth second-order reference; phase and amplitude against analytic solution |
| I02 | Fast/slow magnetosonic and entropy eigenmodes | Observed second-order convergence in smooth amplitudes; correct wave-speed ordering and polarization |
| I03 | Brio–Wu shock tube | No nonphysical states, correct intermediate-wave ordering; refine sequence converges to documented high-resolution reference |
| I04 | Orszag–Tang vortex | Energy, shock morphology and magnetic spectra converge with refinement; conservation and divergence gates hold |
| I05 | Magnetic loop translation in stationary/translating frames | Quantify advection error, loop deformation and induced spurious out-of-plane field; stage order and Galilean behavior documented |
| I06 | Strong magnetized blast, low-beta rotor and current sheet | Positivity, HLLD→HLLE fallback accounting and robust conservation; no hidden diffusivity |
| I07 | MHD turbulent dynamo | Compare magnetic growth and saturation and kinetic/magnetic spectra across resolution and codes; no universal arbitrary absolute target claimed |
| I08 | AMR crossing and reflux-curl | Divergence G03, parent/child EMF interval completion and shock/loop continuity through regrid and MPI boundaries |
| I09 | Expansion-only magnetic field | Homogeneous `B_comov` invariant within 1e-12 relative and `B_phys∝a^-2`; homogeneous `E_mag,comov∝a^-1` |
| I10 | Cosmological waves, gravity coupling | Second-order convergence and no double application of Hubble/gravity sources; external benchmark comparisons |

### 16.3 Non-ideal transport tests

| ID | Benchmark | Acceptance rule |
|---|---|---|
| N01 | Ohmic sinusoidal field decay | Decay-rate relative error ≤1% at 128 cells per wavelength and spatial order ≥1.8 on smooth refinement sequence |
| N02 | Ambipolar diffusion wave/C-shock | Correct diffusion/heating rate and shock profile against analytic or high-resolution reference; convergent, nonnegative thermal heating |
| N03 | Hall whistler circular waves | Phase-speed relative error ≤2% at 128 cells per wavelength, convergent group/phase error and correct polarity reversal |
| N04 | Hall-only energy | No spurious physical Joule heat; conservative total-energy residual ≤1e-9 on closed periodic benchmark |
| N05 | Strongly variable eta/charged-grain table | No extrapolation outside table; correct sign-reversing eta_H, validated interpolation convergence, table hash restart admission |
| N06 | Combined D+H+I | Smooth manufactured solution order ≥1.8 with variable coefficients and nonzero expansion; consistency of split source work and EMF integrals |
| N07 | STS vs explicit | RKL2 and explicit converge to same reference solution under matched accuracy, no unphysical overshoot/floor energy beyond gate G05 |
| N08 | Collapsing magnetized cloud | Resolution study of magnetic braking, mass-to-flux ratio, thermodynamics and angular momentum; independent reference-code comparison, not a single arbitrary observational target |

### 16.4 Cross-code scientific and hardware evidence

Use matched physical equations and units first, with the same boundary conditions and compatible initial data. Compare (a) same effective resolution and (b) **time/memory to a fixed predefined accuracy target**. Required reference partners: Athena++/AthenaK for ideal/AMR; IDEFIX for non-ideal; AREPO for independent cosmological/moving-mesh behavior. Comparing an algorithm at the same cell count but with different PDEs or coefficient tables is **not** a valid numerical accuracy ranking.

For scientific observables use independently calibrated quantities where possible: dynamo spectra, Faraday rotation measures, field–density correlations, magnetized outflow properties and cluster/galaxy diagnostics. No single agreement statistic is treated as proof of correct physics; preserve error bars, model assumptions and numerical convergence. `REF-v1` certification requires classical numerics gates and cross-code agreement on at least four representative problems (Alfvén wave, Orszag–Tang, AMR loop and Ohmic/Hall case), with reference data and normalization published.

## 17. Pareto-frontier evaluation and algorithm promotion

The research comparison vector is `P=(error_vector, wall_seconds, peak_RSS_bytes, peak_device_bytes, MPI_bytes, energy_joules_if_measured, robustness_failures)`. Numerically admissible methods must pass conservation/divergence and physics-validity gates **before** Pareto ranking. A method is non-dominated if no alternative has no-worse measured objectives and strictly improves one. The time-to-accuracy target is the cheapest method configuration satisfying a fixed set of error tolerances on the same problem and hardware.

Three research candidates are formally authorized but **disabled in reference profiles**: (A) high-order entropy-stable FV/compatible DG magnetic representation; (B) vector-potential MHD with compatible curl and a separately specified gauge/AMR transfer; (C) semi-implicit/IMEX Hall and other stiff transport. Each gets a distinct schema ID, dedicated workspaces, identical physical coefficient provider and external benchmark control. Promotion to the reference requires: all reference gates; at least **10% improvement in time-to-fixed-accuracy** or **10% lower error at equal resource budget** in at least two declared target workloads, with the 95% confidence interval of the measured improvement remaining above zero across at least five independent timings; no unexplained >5% regression in peak memory on the 48-GiB class target; no severe regression on any mandatory safety benchmark; and a versioned migration plan. Novelty alone is not a promotion criterion.

## 18. Implementation order, dependency gates and deliverables

| Phase | Mandatory deliverable | Completion gate |
|---|---|---|
| M0: contracts | Freeze geometric face/edge keying, energy/units, `MhdPatchState`, config and versioned schemas; preflight | Compile-independent unit/provenance checks; disabled-MHD regression unchanged |
| M1: ideal single patch | SSPRK2, PLM, HLLD/HLLE, UCT-HLL, exact face CT, energy recovery | I01–I06, G01–G05 |
| M2: accuracy and AMR | PPM, ratio-2 prolong/restrict/reflux-curl, synchronized MPI exchange | G02–G09, I05/I08, MPI and restart gates |
| M3: cosmology & astrophysical coupling | `B_c` conversion, expansion source, source-state contracts, cooling/gravity | I09/I10, conservation and source equivalence |
| M4: non-ideal transport | constant/table providers, conservative D and H, Strang, RKL2, Hall-HLL | N01–N07, nonideal schema and restart tests |
| M5: science qualification | cloud collapse, dynamos, cross-code campaign; one-node and small-cluster profiles | N08, I07 and cross-code evidence; certify REF-I and REF-NI separately |
| M6: research and optimization | optimized UCT/GPU, compatible high-order candidate, IMEX Hall prototypes | Separate Pareto reports; promotion only by Section 17 |

**Workflow policy:** design/source completion may be batched for development convenience, but no phase earns a “certified” label before its tests run and the measurement reports are archived. MHD-disabled hydro remains production-supported throughout. Full distributed subcycling is a separately gated extension of M2 and is explicitly not enabled by default while its MPI temporal-ghost prerequisite remains unvalidated.

## 19. Security, reproducibility and source-compatibility rules

Every external table, IC source and benchmark reference must carry a content digest and documented units; fail rather than guess missing metadata. External paper/code algorithms may be reimplemented from published descriptions, subject to source-code licenses. No unlicensed wholesale code copying from AREPO or other research repositories. Record compiler, flags, git SHA/source archive manifest, CPU/GPU and reproducibility mode as currently expected by CHUÍ. Unknown configuration keys fail, and invalid optional features cannot be silently ignored. Hidden numerical dissipation, coefficient clamping, unaccounted energy injection and silent output schema conversion are prohibited.

## 20. Release decision and fixed acceptance interpretation

The v1 implementation has **two independent official badges of numerical readiness**: `REF-I CERTIFIED` after G01–G09/I01–I10 and at least the required ideal cross-code comparisons; and `REF-NI CERTIFIED` after REF-I plus N01–N08, table-provider acceptance and required non-ideal cross-code comparisons. A failed individual test means the profile is **uncertified**; it does not authorize changing tolerances without a review/version bump. There is no claim that either badge is currently attained.

This specification deliberately closes the questions of field authority, equations, defaults, divergence control, solver selection, AMR strategy, coefficient behavior, source ownership, memory, I/O, failure behavior and tests. It does **not** fabricate benchmark outcomes, validated diffusion tables or production-ready implementation status. Those are empirical deliverables governed by fixed acceptance gates.

## Appendix A — Full reference `.param.txt` example and parser rejection rules

The configuration shown in Section 12 is the complete normative MHD extension for the `mhd_ref_nonideal_v01` profile. It must be combined with the preexisting CHUÍ run, physics, cosmology, memory, output and parallel parameter sections. For `mhd_ref_ideal_v01`, use exactly these changes: `reference_profile=mhd_ref_ideal_v01`, `[mhd_nonideal].enabled=false`, and set all three operator flags to `false`. Unknown keys, accidental whitespace in key names, missing required table files, mismatched profile flags and unsupported units must fail with a path-specific typed configuration error. The parser never guesses from `.yaml` or `.json`.

To run analytic transport tests without a table, use `coefficient_provider=constant`, define `eta_ohmic_phys_cgs`, `eta_hall_phys_cgs`, `eta_ambipolar_phys_cgs` explicitly in `[mhd_nonideal]`, and set `coefficient_table_path` to an empty value. Values are magnetic diffusivities in cm²/s, nonnegative for Ohmic/ambipolar, and signed for Hall. A table provider refuses the constant keys and a constant provider refuses a nonempty table path, preventing mixed authority.

A production example is not considered executable or scientifically certified until its table exists, its metadata are validated and the preflight routine passes. No synthetic microphysics file is supplied or implied by this design document.

## Appendix B — Pseudocode for one accepted global MHD interval

```text
preflight(): validate schema, profile, EOS, transport tables, topology, memory
begin accepted global interval dt at t_n
  freeze legal gas/AMR scheduler active set and source epochs
  for each active synchronized AMR level interval in deterministic order:
    ensure stage-time face/edge MPI ghost consistency
    D(dt/2): SSPRK2 or bounded RKL2, including conservative E flux
    H(dt/2): whistler-HLL CT substeps + conservative E flux
    I(dt): SSPRK2 unsplit HLLD/HLLE gas + UCT-HLL induction
    H(dt/2): symmetric Hall update
    D(dt/2): symmetric Ohmic+ambipolar update
    append time-integrated edge EMFs with stable geometric keys
    accumulate ordinary cell-face fluid fluxes
    reconcile complete gas reflux and magnetic reflux-curl
    check D*C, positivity, energy, memory, time interval completion
  if any failure: restore last accepted interval, halve dt, retry <=8
  after collective acceptance: publish committed physics state
  permit source work, regrid and output only at legal lifecycle boundaries
  permit checkpoint only if no pending MHD stage/register/history/migration
end interval
```

The actual implementation must schedule AMR coarse/fine recursion according to the validated time-coordinator contracts; the pseudocode above shows the logical operations, not a license to apply coarse/fine corrections prematurely.

## Appendix C — Requirements traceability

| Requirement | Owners | Test evidence |
|---|---|---|
| Discrete magnetic divergence | `mhd_geometry`, `mhd_induction`, `mhd_amr` | G02, G03, I08 |
| Coupled total energy | `mhd_riemann`, `mhd_transport`, hydro source adapter | G04, N01–N04 |
| Cosmological units | `mhd_state`, timeline adapter | I09, I10, N06 |
| AMR/MPI/epoch/restart | `mhd_amr`, workflow, I/O | G06–G09 |
| Physical coefficient semantics | `mhd_transport`, chemistry adapter | N05–N08 |
| Numeric robustness | `mhd_riemann`, `mhd_validation` | G05, I03, I06 |
| Resource budget | `memory_governor`, stage workspace, MPI | G09, benchmarks |
| Research Pareto claims | external runner, archives | Section 17 reports |

## References and provenance

**Repository and supplied project references**

- [P1] CHUÍ repository: https://github.com/AlexandreKoz/CHUI-Cosmology-Simulator (public `main`, inspected 22 November 2026).
- [P2] `CosmoSim_Master_Design_Document_Revised-1.pdf`, 11 March 2026: modular Godunov/TreePM/patch-AMR, future MHD and staged verification.
- [P3] `CosmoSim_Variable_Naming_Conventions(1).md`, 6 November 2025: naming, units, HDF5 external compatibility.
- [P4] `CosmoSim_File_Naming_Conventions.md`, 6 November 2025: lower-snake stems; superseded on `.param.txt` by accepted repository ADR.
- [P5] `chui_memory_architecture_and_resource_governance_specification_v01(1).docx`, 23 August 2026: memory authority, budgets, scratch classes, chunked MPI, cost justification.
- [P6] `arepo lessons.pdf`: architectural comparisons, data motion, state separation and scientific fidelity.

**Independent literature**

- [R1] Stone et al. (2020), *The Athena++ Adaptive Mesh Refinement Framework: Design and Magnetohydrodynamic Solvers*, ApJS 249:4. https://doi.org/10.3847/1538-4365/ab929b
- [R2] AthenaK documentation, MHD methods and solvers, accessed 22 November 2026. https://ias-astrophysics.github.io/athenak-docs/mhd.html
- [R3] AthenaK framework documentation. https://ias-astrophysics.github.io/athenak-docs/
- [R4] Lesur et al. (2023), *IDEFIX: A versatile performance-portable Godunov code for astrophysical flows*, A&A 677:A9. https://doi.org/10.1051/0004-6361/202346005
- [R5] Zier, Springel & Mayer (2024), *Non-ideal MHD on a moving mesh I: Ohmic and ambipolar diffusion*, MNRAS 527, 1563. https://academic.oup.com/mnras/article/527/1/1563/7324843
- [R6] Zier, Mayer & Springel (2024), *Non-ideal MHD on a moving mesh II: Hall effect*, MNRAS 527, 8355. https://academic.oup.com/mnras/article/527/3/8355/7459936
- [R7] Springel, Mayer & Pakmor (2026), *Divergence-free magnetohydrodynamics on a moving mesh*, preprint. https://arxiv.org/abs/2609.35989
- [R8] Balsara (2001), *Divergence-free adaptive mesh refinement for magnetohydrodynamics*, JCP. https://www.sciencedirect.com/science/article/abs/pii/S0021999101969177
- [R9] Tóth & Roe (2002), *Divergence- and curl-preserving prolongation and restriction formulas*, JCP. https://www.sciencedirect.com/science/article/pii/S0021999102971202
- [R10] Gardiner & Stone (2008), *An unsplit Godunov method for ideal MHD via constrained transport in three dimensions*, JCP. https://www.sciencedirect.com/science/article/pii/S0021999107005669
- [R11] Mignone & Del Zanna (2020), *Systematic construction of upwind constrained transport schemes for MHD*, JCP. https://www.sciencedirect.com/science/article/pii/S0021999120305222
- [R12] Mattia et al. (2026), *A robust super-time-stepping scheme for Ohmic and ambipolar diffusion*, A&A. https://doi.org/10.1051/0004-6361/202659730
- [R13] Springel, Mayer & Pakmor (2026), vector-potential formulation; research alternative, not reference requirement. https://arxiv.org/abs/2609.35989
- [R14] Pakmor et al. (2013), *Simulations of magnetic fields in isolated disc galaxies*, cosmological magnetic-field conventions. https://doi.org/10.1093/mnras/stt428
