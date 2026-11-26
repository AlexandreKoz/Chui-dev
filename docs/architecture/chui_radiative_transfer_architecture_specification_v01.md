# CHUÍ — Radiative Transfer Architecture Specification

**Version:** 0.1 · **Date:** 2026-11-26 · **Status:** proposed reference design  
**Parent:** [Master design](chui_rt_bh_agn_master_design_v01.md)

## 1. Scope and model selection

P0 SHALL provide time-dependent nonrelativistic multigroup two-moment transport with M1 closure, second-order finite volume, two-stage asymptotically correct IMEX, full-c and separately qualified reduced-c modes, source spectra, matter exchange, cosmological expansion and patch AMR/MPI. P1 adds a direct-ray backend, VET closure, hard-spectrum and dust IR coupling and Monte Carlo reference tools. FLD remains a comparator for diffusion; directional ordinates and adaptive angular/spectral refinement are research. No claim of general relativistic RHD or all-frequency scattering accuracy. The existing chemistry/cooling/ISM campaign owns rates and species. The RT solver must not implement a parallel chemical network.

| Representation | Why use it | Limitation / role |
|---|---|---|
| M1 | fixed 4 moments/group, local stencils, source-count-independent transport | cannot distinguish crossing beams; P0 |
| Direct adaptive rays | directional fidelity for sparse luminous sources | source-dependent cost and irregular exchange; P1 |
| VET | nonlocal angular information with moment evolution | closure solve cost and time dependence; P1 |
| Monte Carlo | independent spectra/scattering reference and observables | sampling variance; P1 offline |
| S_N | explicit directional transport | angular artifacts, memory, sweep dependencies; research |

## 2. Governing equations and frame

For each physical frequency group g, lab-frame radiation energy density E_g, energy flux F_g and pressure P_g:

$$\partial_t E_g+\nabla\cdot\mathbf F_g=j_g-cG_g^0,$$
$$\partial_t(\mathbf F_g/c^2)+\nabla\cdot\mathbb P_g=\mathbf s_{p,g}-\mathbf G_g.$$

`j_g` and `s_p,g` denote external emission sources (with corresponding emitter recoil if modeled). The four-force convention gives matter +c G^0 and +G, radiation the negative; do not assign sign independently in two modules. Full-c closed static-frame exchange preserves gas+rad energy and gas momentum+F/c², including chemical energy if tracked separately. Derive group-dependent opacities and mixed-frame O(v/c) terms from [R1,R2]; absorption/emission, scattering, work, frame transforms and group boundary terms must follow the same order of approximation. Static opacity formulas alone do not qualify dynamic diffusion. In static gas with elastic isotropic transport scattering, G=chi_tr F/c, where chi_tr is inverse length. Absorption and scattering have different energy terms.

M1:

$$f=|\mathbf F|/(cE),\quad \chi_M=(3+4f^2)/(5+2\sqrt{4-3f^2}),$$
$$\mathbb P=E[(1-\chi_M)\mathbb I/2+(3\chi_M-1)\hat{\mathbf F}\hat{\mathbf F}/2].$$

At F=0 set P=E I/3 directly; at E=0 require F=0. Physical admissibility E>=0 and |F|<=cE. M1 is an angular closure, not an estimator of its own missing angular information.

## 3. Canonical state and spectral contract

P0 stores `rad_energy_density` and three `rad_flux` lanes, fp64, per group per authoritative patch cell; exact stored comoving convention below. Covered parents may be derived except when required for time integration; ghost and stage copies are not independent truth. Persistent per-cell photon-number state is not P0. For ionizing groups use a declared reconstructed within-group shape to derive photon counts, cross sections and excess energy per absorption. The discrete group photon diagnostic N_g=E_g/<h nu> is closure-dependent and is not an independent conserved moment when the mean changes. Validate absorption against the same spectral quadrature used for chemistry; if a science case demands independent number and energy transport, use a separately qualified extension rather than promising both exact integrals from one scalar. Track unresolved within-group spectral error. For photons removed by absorption, enforce energy partition and stoichiometry in the local solver, including ionization thresholds and secondary-electron yields where enabled.

Proposed physical band edges in eV (profiles select subsets with explicit coverage; no hardcoded group count):

| Profile | Edges / coverage | Notes |
|---|---|---|
| grey | declared model bolometric range | RHD reference; not a general ionization model |
| ionizing3 | 13.6,24.6,54.4, upper bound from spectrum | H/He ionization; hard-tail validity limited |
| broad UV/IR | 0.01,1,5.8,11.2,13.6,24.6,54.4,100 | illustrative boundaries, dust/stellar groups depend on data |
| AGN hard | broad profile plus 100,2000,100000 | X-ray photoheating; Compton and secondary electron treatment need qualification |

An explicit upper/lower boundary policy SHALL account for photons redshifted out of the modeled range. No interpolation silently outside tabulated domains. Tabulate group integrals of emissivity, absorption, scattering, photon number and heating against declared spectra; version all table generators and source datasets. User-selectable group definitions require compatible chemistry and opacity metadata. No arbitrary grey X-ray opacity for production.

## 4. Cosmological representation

Use proper E_g and F_g as physics interface; proposed canonical moment densities `E_comov=a^4 E_phys`, `F_comov=a^4 F_phys` for radiation on a comoving volume, at physical frequency bin edges. For collisionless radiation in homogeneous FLRW with peculiar velocity zero, ignoring gravitational redshift and metric perturbations:

$$\partial_t E_\nu+3HE_\nu-H\nu\partial_\nu E_\nu+a^{-1}\nabla_x\cdot F_\nu=0.$$

Group integration yields

$$\partial_t E_g+4HE_g+a^{-1}\nabla_x\cdot F_g=H[\nu E_\nu]_{\nu_l}^{\nu_u}.$$

The same cosmological dilution/redshift structure applies to each flux component with its angular moment transport. Multiply by a^4 to remove the homogeneous dilution; use a conservative frequency-space flux toward lower physical frequency (upwind or positivity-preserving reconstruction of spectral edge state). Summing groups cancels internal edge transfers. Physical energy in proper volume scales as a^-1 for homogeneous freely propagating photons; comoving photon count remains constant in a closed full-spectrum calculation. Track boundary leakage and expansion work explicitly. Convert local physical interaction increments by a^4, and geometric divergence by 1/a; do not reuse hydro's a³ density conversion for radiation. Evaluate scale factor at stage time, integrate proper dt from the cosmology service. This proposed convention must be reconciled with the actual hydro integrator; no independent cosmological clock. RSLA applies to spatial transport and interactions as derived below, not arbitrarily to H or frequency redshift. A cosmological reduced-c equation and its invariants require separate validation before enablement.

## 5. Spatial discretization

Finite volume over volume V and faces A: dU_i/dt= -sum_faces A Fhat/V + S_i. Conserved radiation vector for transport may use (E,F/c²) in full-c form. Use unsplit second-order piecewise-linear MC reconstruction, limit slopes jointly to preserve E>=0 and |F|<=cE at every face; first-order cell-average fallback must still obey the actual CFL and realizability. PPM is an optional higher accuracy candidate, not assumed compatible automatically. HLL with valid extremal characteristic speed bounds; GLF/Rusanov uses global safe signal bound as per-face fallback. Equal-and-opposite face fluxes shared by neighbors, never separately recomputed with asymmetric rounding/limiting. Fluxes integrated with stage weights; conserve correction at AMR interfaces. Test rotational symmetry and oblique beams.

For isotropic thick static scattering, P≈E I/3 and flux relaxation yields F≈-c grad(E)/(3 chi_tr). Physical diffusion coefficient D=c/(3 chi_tr). First-order GLF diffusion scales c dx/2, potentially much larger than D. Adopt published He et al. PD-ARS scheme with consistent reconstruction and flux; do not add an ad hoc opacity wave-speed multiplier by default. The paper demonstrates its scheme without the empirical wave-speed correction; modified HLL/GLF fallback and AMR extensions must independently pass diffusion gates. If they fail, qualify an asymptotic-preserving spatial correction or alternative flux as an explicit amendment, not a hidden tuning knob.

## 6. Time integration

Let T(U,t) be spatial transport and S(U,t) local stiff exchange. Reference He et al. 2024 equations (28)–(29), with source terms evaluated at the implicit end stage:

$$U^*=U^n+\Delta t\,T(U^n,t_n)+\Delta t\,S(U^*,t_{n+1}),$$
$$U^{n+1}=U^n+\tfrac{\Delta t}{2}[T(U^n,t_n)+T(U^*,t_{n+1})]+\tfrac{\Delta t}{2}[S(U^*,t_{n+1})+S(U^{n+1},t_{n+1})].$$

The first intermediate state is evaluated at end time, not a physical half-step despite the source paper's superscript notation. Both stages retain stiff interactions. Streaming is second-order; do not claim second-order temporal convergence in asymptotic diffusion. Nonautonomous emission events need their own time-integral contract; do not inject full luminosity in both stages. An explicit or split comparator is permitted but not the certified thick-limit path. Reference transport CFL `dt <= C_rt min(dx_phys/c_eff)` (initial C_rt=0.3 provisional; verify multidimensional sum-of-speeds condition). Subcycle relative to hydro interval; adjust final substep to exact boundary. Never assume implicit local exchange removes the transport CFL. A faster optically thick step is research unless proven stable and accurate.

### Local nonlinear solve

The existing thermochemistry operator exposes reaction rates, heating/cooling, species and energy; extend it with a shared radiation exchange residual and Jacobian or derivative API. Unknowns include gas thermal energy, independent species fractions, group E and flux/momentum as needed. Use local elimination for linear flux relaxation where valid, block Jacobian/Schur complement for independent group couplings, otherwise damped Newton with bounds and residual control. Do not use stale opacity for a strongly coupled Newton solve unless convergence and energy error are controlled. Use positivity-preserving line search, elemental and charge constraints, finite iteration cap and bounded substeps. Define residual weights per conserved quantity; initial relative residual 1e-10, absolute tolerance per units/state scale. If no convergence, roll back, reduce dt and retry with bounded count (initial 8); fail with cell ID, regime, residual and state if exhausted. Do not repair with untracked energy clipping. Return a staged delta and ledger; gas and radiation commit together. Chemistry owns chemical binding and escaped cooling. If cooling photons are emitted into transported groups, subtract energy only once. Include radiation work in gas total energy; calculate kinetic update from actual momentum exchange.

### Reduced-speed-of-light approximation

For mixed-frame convention [R1], use physical flux F, physical c in opacity/closure, and reduced c_hat in radiation time evolution:

$$\partial_t E+(\hat c/c)\nabla\cdot F=-\hat c G^0,$$
$$\partial_t(F/(c\hat c))+\nabla\cdot P=-G.$$

Matter receives cG^0 and G. Closed source exchange invariant is E_gas+(c/c_hat)sum E_g and p_gas+sum F_g/(c c_hat). Physical full-c balance only holds at c_hat=c. Source injection, escape flux and emission must be scaled consistently to preserve approximation invariants; report both approximation and physical diagnostics. Do not mix a photon-number RSLA convention with this one. `c_hat/c` in (0,1], initially globally constant; full-c reference. Reduced-c selection requires convergence of target observables (ionization fronts, reionization history, equilibrium ion fractions, dynamic diffusion and AGN variability) at multiple c_hat. No automatic decrease under resource pressure. Level-dependent c and implicit transport are not P0.

## 7. AMR, parallelization and sources

- Use geometric patch storage, ratio-2 prolongation with conservative parent averages and joint realizability limiter; volume-weighted restriction preserves admissibility via convexity. Flux-register correction must preserve both group moments; if inadmissible, use coordinated limited flux correction or rollback, never independent post-reflux clipping.
- Time-weighted flux registers include both IMEX transport stages. Coarse temporal ghosts interpolate admissible states; if cross-rank temporal AMR unsupported, use synchronized updates. Parent cells covered by finer leaf cells cannot absorb/emit a second time. Regrid only at consistent boundaries and migrate moment state, ownership and history together.
- Spatial MPI ghost exchange per stage; rank boundary uses one flux authority or proven equivalent symmetric computation. Bounded exchange, no global radiation cell or source replication. Use patch-level work estimates reflecting RT subcycles and source Newton iterations.
- Stellar spectra supplied from the existing age/metallicity/IMF provider. AGN spectra from actual BH budget. Deposit luminosity over normalized leaf kernel; source timing integrated over active interval, no double count on restart. For isotropic sources net emission momentum zero in rest frame; anisotropic and moving-source recoil explicitly modeled or named approximation. External boundary profiles distinguish vacuum/outflow, incident intensity and periodic domains; reflecting boundary normal flux reverses. Cosmological radiation periodicity is a declared choice, not implied by gravity; avoid artificial background double counting.
- P1 direct rays: launch angular quadrature with adaptive splitting, carry integrated photon/energy budgets and time, attenuate using cell optical depth and traverse AMR leaves once per segment. Rank migration bounded; deterministic ray IDs. Direct and diffuse fields must be partitioned: photons absorbed, escaped, or converted once to diffuse group moment source with energy/momentum bookkeeping. Never inject direct luminosity in both ray and M1. Finite propagation speed or explicit steady-state limit required.
- P1 VET: compute normalized positive semidefinite tensor, trace 1, from angular solution; closure update cadence based on source/opacity changes, causality caveats documented. VET may need distinct implicit transport backend; do not assume a generic tensor plug-in preserves M1 hyperbolicity or realizability. Verify moment/closure consistency, nonlocal communication and accuracy; initially qualify grey VET before multigroup.

## 8. Interface and implementation proposal

Suggested paths (not promises of current code): `include/cosmosim/radiation/radiation_state.hpp`, `radiation_transport.hpp`, `radiation_closure.hpp`, `radiation_imex.hpp`, `radiation_spectrum.hpp`, `radiation_amr.hpp`, `radiation_matter_coupling.hpp`; corresponding `.cpp` files in `src/radiation/`. `RadiationPatchState` owns moment arrays, `RadiationTransportBackend` declares capabilities and state schema, `RadiationSourceProvider` returns time-integrated group emission, `RadiationMatterCoupler` uses existing chemistry service. Typed views avoid copying all cells. Reuse MHD/hydro patch metadata and memory reservations. No radiation allocation when disabled; prove disabled equivalence. Exact API must follow existing source runtime resource keys and stage contracts.

## 9. Tests and qualification

| Gate | Reference / measurement |
|---|---|
| Closure algebra | isotropic/free streaming, rotational invariance, PSD, trace and bounds |
| Vacuum transport | pulse speed, oblique beam, manufactured smooth convergence >=1.8 spatial/time streaming at full-c |
| Local absorption/emission | analytic exponential, LTE equilibrium, radiation/matter energy+momentum, charge/element conservation |
| Static and dynamic diffusion | Marshak wave, moving pulse, advecting uniform radiation, shock; compare [R1,R2], no spurious leading numerical diffusion; vary tau and dx |
| Ionization | Strömgren radius R_S=(3 Q/(4 pi alpha_B n_H²))^(1/3), time evolution in stated on-the-spot limit; compare RAMSES-RT test suite |
| Angular | crossing beams and shadow: document M1 failure envelope, P1 rays/VET improve matched metric |
| Spectral | group-edge redshift, absorption hardening, blackbody integration, hard-spectrum reference rates |
| Cosmological | homogeneous a^-4 physical energy density, conserved full-spectrum comoving photon number in reference spectrum, expansion ledger |
| AMR/MPI | moving front across patch/rank interfaces, reflux conservation, regrid/migrate, restart, empty ranks, source ownership |
| Coupling | radiation pressure acceleration and MHD radiation shock; no duplicate photoheating or cooling, source-only ledger |
| Performance | time-to-error, peak host/device, stage bytes and communication scaling; CPU/GPU equivalence |

Quantitative proposed tolerances are in master §8; specific test norms/reference solutions and refined-grid thresholds must be pinned before qualification. M1 angular errors are model errors, not acceptable hidden failures in a science case requiring directional fidelity. Full-c and RSLA are separate certified profiles.

## References

- R1. He, Wibking & Krumholz (2024), asymptotically correct IMEX RHD, MNRAS 531, 1228: https://arxiv.org/abs/2404.08247 (particularly equations 25–29, diffusion analysis).
- R2. He, Wibking & Krumholz (2024), mixed-frame multigroup QUOKKA: https://arxiv.org/abs/2407.18304
- R3. Wibking & Krumholz (2022), QUOKKA: https://arxiv.org/abs/2110.01792
- R4. Rosdahl et al. (2013), RAMSES-RT: https://arxiv.org/abs/1304.7126
- R5. Menon et al. (2022), VETTAM: https://doi.org/10.1093/mnras/stac485

The specific CHUÍ cosmological/group/chemistry integration and API are proposed engineering extensions; the cited results do not establish their correctness. Opacity tables, spectral datasets, exact current MHD/source interface and complete P1 closure choices are preflight/qualification dependencies.
