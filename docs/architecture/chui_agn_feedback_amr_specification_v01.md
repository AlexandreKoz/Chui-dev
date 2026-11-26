# CHUÍ — AGN Feedback, Conservative Injection and AMR Specification

**Version:** 0.1 · **Date:** 2026-11-26 · **Status:** proposed reference design  
**Parent:** [Master design](chui_rt_bh_agn_master_design_v01.md)

## 1. Scope and decisions

P0: continuous and accumulated thermal feedback, conservative kinetic kicks, wide-angle mass-loaded winds, bipolar volumetric jets, spectral emission, and resolution-aware injection/outflow refinement. P1: resolved nozzle, magnetic injection compatible with the MHD constrained-transport (CT) authority, spin-powered budgets and detailed radiation. Research: specialized wind carriers, relativistic propagation, cosmic rays and adaptive closure switching. Mechanical feedback is not limited to mutually exclusive high/low accretion modes; multiple channels may operate with a common budget. The inner engine sets energetics; AGN maps those budgets to resolved gas. Cooling/chemistry/ISM, gravity, and magnetic field evolution remain owned by their existing modules.

## 2. Budget, units and transaction model

The BH/disk engine issues integrated proper-frame quantities per source interval: radiative energy per band, thermal energy, mechanical energy (kinetic/thermal/magnetic partitions), ejected baryonic mass and species, momentum and angular momentum, emission interval, orientation, source velocity, model ID and source/event ID. The engine separates accretion-powered energy, spin extraction and previously accumulated reserves. The AGN module may redistribute only within the allowed budget; it cannot increase total energy or infer a second luminosity from BH mass. Pending energy remains restart truth until released or explicitly lost; radiation and mechanics cannot debit the same energy twice. Escape and unresolved absorption have named terms. The Newtonian ledger and physical mass-defect ledger described in the BH specification cannot be mixed. In a relativistic inner model, rest mass is not ordinary gas kinetic/thermal energy.

All weights act on integrated quantities, not density. At physical scale factor a, evaluate proper support volume, mass, and peculiar velocities; convert the resulting conserved increments once to solver code units and canonical cell averages. For nonrelativistic ejection in the BH rest frame, the simulation frame momentum increment is dp=dp_rest+dm v_bh and energy increment dE=dE_rest+v_bh·dp_rest+0.5 dm |v_bh|². Corresponding donor motion/recoil must be accounted for; radiation momentum uses the declared relativistic/order-v/c emission model, not a Galilean photon boost. Newtonian jets must be restricted to a declared subrelativistic domain, with initial proposed max v/c=0.1 (not a universal astrophysical boundary). Reject incompatible injection models rather than cap velocities silently.

## 3. Leaf support and MPI commit

A feedback kernel samples only leaf volumes intersecting physical support, excluding covered parents and ghost duplicates. Use cell-volume integral of nonnegative kernel K, with center sampling only when qualified by convergence. Normalized scalar weights w_i= integral_cell K dV / sum_j integral_cell K dV. Angular kernels may include lobes and disk orientation, resolved geometry, and target mass constraints. For vector momentum, normalized scalar weights alone do not ensure balanced lobes; compute each lobe vector sum explicitly. A bipolar source in the rest frame must have equal and opposite net momenta if the selected model promises zero recoil. Otherwise transfer reaction to the source. Record angular momentum about the source plus orbital term, including any unresolved torque.

Ownership protocol: source owner queries intersecting patch owners, gets bounded aggregate statistics, creates a deterministic plan for recipient stable cell IDs and topology epoch, all participants check current donor/cell feasibility and reserve memory, then commit at a legal common scheduler time. Rank participation and transaction sequence are deterministic in debug mode. No globally replicated recipient list, no unbounded event-by-cell records. Competing events use deterministic ordering or joint budgets. If no valid support, retain budget and diagnose; never deposit into a stale host. Strong events wake/synchronize affected gas and invalidate primitive, thermal and force caches as required. A rejected plan leaves source and targets unchanged. An interrupted distributed commit cannot resume from partial state; recover last fully committed checkpoint.

## 4. Mechanical channels

### 4.1 Thermal injection

Continuous thermal energy Q_th delivered by normalized weights, dU_i=w_i Q_th. Reservoir mode accumulates until Q_th >= M_heat * integral_T^(T+dT) c_v(T',Y) dT' (constant cv only for ideal gas). Select a target mass/heating increment consistent with actual EOS, cooling and resolution. Stochastic event targeting, if used, must preserve exact total debited budget and deterministic counter-based RNG keyed to BH ID and event ID; statistical spatial variation allowed, not stochastic creation of energy. Track waiting times, injection temperatures and radiative losses. No undocumented delayed cooling. If the desired event cannot fit in allowed mass or budget, keep reserve or use a user-selected explicitly distinct model, not arbitrary temperature changes.

### 4.2 Momentum-only / kinetic kicks

For cell with integrated mass m and momentum p, after adding dp at fixed mass, dK=v·dp+|dp|²/(2m). For a trial directional pattern dp_i=A q_i, sum energy change is A B + A² C, B=sum v_i·q_i, C=sum |q_i|²/(2m_i). Given allocated kinetic energy Q>=0, solve C A²+B A=Q with nonnegative root using numerically stable quadratic evaluation; handle C=0 and Q=0 degeneracies explicitly. Momentum pattern q_i includes recoil or equal-opposite constraints; if recoil changes source kinetic energy include that contribution in the same budget. Do not force both independently prescribed arbitrary dp and Q when inconsistent. Check internal energy remains admissible. Kinetic energy can be redistributed or dissipated only through resolved dynamics or explicit physical source terms, not hidden energy correction.

### 4.3 Mass-loaded winds and bipolar jets

For actual ejected mass dm, wind velocity v_w and thermal specific energy u_w in source frame, dp=dm v_w, dE=dm(u_w+0.5|v_w|²) plus other declared components. Wind mass comes from an explicitly debited disk or entrained gas reservoir; entrainment of ambient gas is a transfer, not new mass. Species and metal masses sum consistently with total mass. With prescribed nonrelativistic kinetic luminosity L_k and outflow speed v_w, mdot_w=2 L_k/v_w² and momentum flux=mdot_w v_w. Alternatively prescribe mass loading and energy, derive speed. Three independently chosen incompatible values must be rejected. For jets, choose cone angle, injection radius, axis from disk/spin or named model, and separate lobe weights; equal-power equal-mass opposite lobes by default. Direction changes follow physical source model and time; no random reorientation unless explicitly selected. Source orientation singular at zero spin handled by disk axis or named fallback.

Integrated cell update:

$$m'=m+dm,\quad p'=p+dp,\quad E'=E+dE,$$
$$\Delta K=|p+dp|^2/[2(m+dm)]-|p|^2/(2m),$$
$$\Delta U=\Delta E-\Delta K-\Delta E_B-\Delta E_{other}.$$

Canonical hydro/MHD total energy remains authority; update internal-energy auxiliary consistently with existing dual-energy/chemical conventions. Check positivity after transaction, including existing fluid velocity and magnetic terms. Mixing injected material at different velocity generates the corresponding thermal energy through conservative update. It is not permissible to deposit the requested kinetic energy again as heat. No double count with stellar feedback if both target same gas interval.

### 4.4 Resolved nozzle (P1)

Define internal control surface in proper coordinates with resolved face area and direction. Integrate mass, momentum, enthalpy, kinetic and magnetic energy flux through actual surface; distinguish volumetric source from boundary flux (only one applied). Source rest-frame normal and motion explicit. If embedded nozzle geometry cannot be faithfully represented by Cartesian AMR, qualify geometric quadrature or reject that setup. Ensure correct pressure matching/overpressure, boundary Riemann flux, no mass creation beyond budget. Validate against volumetric injection at converged scales; do not force nozzle model when its physical diameter is unresolved.

## 5. Radiation channel

Given allocated L_rad, use normalized spectral shape phi_nu so integral phi_nu dnu=1 and L_g=integral_group L_rad phi_nu dnu. Integrate over actual source interval and transfer budget to RT once. Spectra depend on accretion state, spin and obscuration only via versioned providers; no invented X-ray fraction or AGN SED. If a torus/unresolved absorption model is used, track absorbed energy, local matter coupling and IR re-emission separately. Angular distribution normalized over solid angle; direct ray solver for anisotropy P1. Isotropic M1 deposition is valid only within its angular model. In no-RT runs, a named local empirical coupling may consume a declared fraction and record escape; it cannot also operate on already transported luminosity. Background radiation and stellar luminosity are owned by their separate providers. No duplicate photoheating in thermochemistry. Photon counts, radiation source units and reduced-c normalization follow the RT specification.

## 6. Magnetic feedback (P1)

The existing MHD design uses Lorentz–Heaviside normalized magnetic fields, staggered oriented face-integrated flux, CT and reflux-curl. Magnetic jet injection MUST use that geometry, not independent cell-centered B updates. For normalized B, e_B=B²/2 and ideal Poynting flux S_B=B² v-(v·B)B. Gaussian units require factors 1/(8 pi) and 1/(4 pi), respectively; convert exactly once. Construct compact edge vector-potential increment A or time-integrated edge EMF such that discrete dPhi=-curl_edge(EMF) preserves discrete divergence (sign based on repository convention). Net boundary flux and return field must be compatible with selected topology. Use poloidal/toroidal profiles with regular axis, smooth support, declared helicity and no magnetic monopole. Update through the same CT authority at a synchronized stage, including AMR edge reflux-curl and compatible restriction/prolongation.

For B'=B+s b (b derived from a compatible unit-amplitude curl), the actual magnetic-energy change using MHD's canonical energy reconstruction is

$$\Delta E_B(s)=s\sum_i V_i B_i\cdot b_i+\tfrac12 s^2\sum_i V_i |b_i|^2.$$

Solve for amplitude using available magnetic budget and consistent energy quadrature, not norm of injected field alone. The cross term can be negative; no automatic deletion of existing field energy. If a feasible root and CT-consistent total-energy update do not exist, reject/defer event and report. Sum over owned leaf cells only. Account for thermal/kinetic changes in same transaction, use actual post-CT energy if discrete representation differs. A nozzle Poynting flux and volumetric curl must not double-inject energy. Verify divergence remains within MHD reference tolerances on regrid and migration, not just on a uniform isolated grid.

## 7. Refinement policy

Arkenstone BH motivates avoiding dilution of high-specific-energy outflows into large gas masses; it does not establish that a Cartesian AMR analogue is already validated. CHUÍ uses targeted refinement around injection and downstream structure rather than copying moving-mesh wind carriers.

| Indicator | Criterion / measurement | Guard |
|---|---|---|
| Injection geometry | nozzle diameter / dx_phys >= N_j | choose N_j by convergence (trial 8,16,32); mandatory for nozzle, not universal physical constant |
| Dilution | actual heated/ejected mass and e_inj=Q/m_coupled relative to requested specific energy | do not refine infinitely in dense gas; record unmet requirement |
| Cooling scale | l_cool≈v_postshock t_cool; l_cool/dx >= N_cool | only in identified radiative shocks; cooling service supplies time; unresolved diagnostic at max level |
| Jet tracer | optional conservative passive scalar for injected material | bounded [0,1], no new chemistry; avoid refining entire mixed halo |
| Pressure/enthalpy | dimensionless slope dx |grad q| / max(|q|,q_floor) | threshold and hysteresis explicit |
| Shear/magnetic gradients | vorticity/current indicators normalized by local characteristic scales | P1, noise and bandwidth guards |

Use priority classes: mandatory injection adequacy for selected resolved model, scientifically important shocks/interfaces, optional downstream refinement. Apply hysteresis, buffer cells/patches, minimum residence and safe derefinement cadence. Refine before deposition only at legal topology boundary; hold budget pending if required. Recompute recipient plan after regrid. For derefinement preserve conserved energy, momentum, species, magnetic face flux and radiation moments; report artificial mixing and resolution loss. A memory denial of mandatory refinement fails closed or defers only if latency is explicitly permitted; it never selects a different physical channel. Optional additional resolution within an explicitly configured error/resource policy may be skipped with diagnostics. At maximum level, either the explicitly chosen unresolved model remains valid or the requested resolved model is rejected. No hidden model switch. Refinement priorities feed memory reservations and work-weighted decomposition; 3D fixed-volume refinement multiplies leaf cells by eight per level and may double subcycle count, so runtime cost is not just memory.

## 8. Execution, observability and proposed files

Suggested `include/cosmosim/physics/agn/{agn_budget,agn_injection,agn_kernel,agn_refinement,agn_spectrum}.hpp`, matching `.cpp` under `src/physics/agn/`; magnetic coupling adapter owned jointly with existing MHD interface. `AgnFeedbackBudget` immutable request; `AgnInjectionPlan` bounded phase-local geometry and ownership; `AgnInjectionReport` actual applied ledger, failure/defer reason and refinement status. Reuse stellar feedback transactions where semantically identical; do not create a parallel feedback authority. Event history compact per source; full per-cell audit only sampled/debug. Include source ID, mode, actual mass/energy per channel, lobe imbalance, recoil, reservoir balance, temperature, injection support size, tracer distribution, radiated/escaped energy, energy/momentum residual, and resolution inadequate flag.

## 9. Qualification matrix

| Test | Required result |
|---|---|
| Single/multiple cell injection | exact ledger within fp64 tolerance, positive thermal state, moving-frame energy/momentum correctness |
| Opposing lobes | balance including irregular AMR/MPI support, recoil when asymmetric, angular momentum ledger |
| Reservoir duty cycle | dt refinement, no double release, restart exact event identity |
| Thermal blast | Sedov–Taylor limit without cooling; radiative blast compare qualified cooling module |
| Wind bubble | analytic energy-driven bubble scaling R proportional (L/rho)^(1/5)t^(3/5) where assumptions apply |
| Jet propagation | head speed, cocoon pressure, cavity volume/enthalpy, orientation dependence, resolution convergence |
| Magnetic injection | discrete divB, actual magnetic energy including cross term, AMR CT compatibility, jet topology |
| Radiation | spectral luminosity integral, escape/absorption/debit, no duplicate heating |
| MPI | uneven ranks, overlapping supports, simultaneous donors, rejection/rollback, migration and checkpoint |
| Isolated galaxy/cluster | outflow mass/energy loading, cavity sizes and entropy profiles, cooling-flow balance, comparison at matched model/resolution |
| Cosmological zoom | BH growth, AGN duty and luminosity distributions, quenched fractions, baryon fractions, resolution/calibration sensitivity |

Qualification gates in master §8 apply. Pre-register model-dependent observational calibration targets and held-out metrics; identical total injected energy does not imply equivalent coupling. No empirical efficiency or refinement threshold is certified by this document.

## References

- A1. Arkenstone BH (2026), high-specific-energy feedback, isolated-galaxy demonstration: https://doi.org/10.1093/mnras/stag1520
- A2. Huško et al. (2026), hybrid AGN winds/jets/spin: https://doi.org/10.1093/mnras/stag324
- A3. Davé et al. (2019), SIMBA feedback: https://doi.org/10.1093/MNRAS/STZ937
- A4. CHUÍ repository MHD architecture and CT contracts, `docs/architecture/chui_mhd_architecture_specification_v01.md`.

The injection weighting, bounded distributed commit and refinement policy are CHUÍ design proposals. They require engineering verification and physical convergence studies; they are not claimed to be published Arkenstone algorithms.
