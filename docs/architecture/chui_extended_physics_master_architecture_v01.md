# CHUÍ — Extended Physics Master Architecture

**Version:** 0.1 · 27 November 2026  
**Status:** Accepted conceptual direction; proposed normative design, pending detailed subsystem specifications and qualification  
**Target:** Powerful desktops, workstations and small-to-medium HPC systems  
**Scope:** Massive neutrinos, radiation hydrodynamics, cosmic rays, Lyman-alpha, CMB, radio astrophysics and synthetic instruments  
**Repository basis:** `Chui-dev-main(8)(1).zip` (latest supplied snapshot); source-only review

## 1. Scientific mandate and design decisions

CHUÍ's ambition is to become the most physically faithful numerical representation of the Universe that can be constructed within a declared scientific and computational envelope. Advanced means fidelity, not the largest particle count, longest feature list or most elaborate mathematics. The code should first achieve feature, numerical and scientific parity with established methods in matched regimes, then pursue improvements that demonstrably reduce physical or numerical error. No single code is a universal parity reference.

The reference architecture is desktop/workstation-first but must scale to distributed HPC. Hardware constraints determine feasible discretizations and model choices, not silent omissions of relevant physics. The existing master design and resource governance specification remain in force; this document extends, rather than replaces, them.

The following requirements are proposed as normative. MUST/SHALL denotes an intended requirement, not a claim of implementation. Three method families SHALL coexist: **reference**, using transparent established algorithms and independent comparators; **production**, using qualified optimizations and scientifically declared approximations; and **experimental**, opt-in research methods with independent acceptance gates. A high-fidelity method may be a parity target without being an affordable default. A production dispatch path is not scientific qualification.

| Domain | Affordable/reference comparator | Preferred advanced parity target |
|---|---|---|
| Neutrinos | Background, transfer ratio, linear response, ordinary particles | Weighted particle delta-f |
| Radiation | Explicit multigroup M1, FLD for diffusion regime | M1 + IMEX; distinct VET backend; direct/diffuse rays |
| CR transport | Grey fluid, advection, diffusion/streaming | Spectral two-moment protons and electrons |
| CR spectra | Fixed momentum bins, steady-state electron model | Piecewise power-law groups; time-dependent spectra; tracer option |
| Ly-alpha | Voigt absorption; unaccelerated scattering reference | Resonant Monte Carlo, qualified acceleration and coupling |
| CMB | CLASS/CAMB, Born lensing, nonrelativistic SZ | Multi-plane lensing, relativistic tSZ, kSZ, ISW/RS |
| Radio | Thin HI, scalar transfer, analytic screens | Full optical-depth HI, full Stokes and evolving CR electrons |
| Instruments | External OSKAR/pyuvsim, direct summation | Native full RIME and accelerated transforms |

The user accepted this overall direction in the design discussion. Specific schemas, tolerance thresholds, implementation sequencing and resource coefficients remain proposed; they require review in the detailed specifications. No new research algorithm is claimed to be novel or superior without a dedicated literature and validation campaign.

## 2. Source baseline and blockers

This snapshot is not a certified release. The RT/BH/AGN pass-01 ledger explicitly marks all tests NOT_RUN; its 35 rows comprise 17 IMPLEMENTED_UNVALIDATED, 16 PARTIAL, one FAIL_CLOSED and one EXCLUDED_RESEARCH. The MHD pass-02 ledger and runtime capability guards likewise describe provisional, restricted dispatch. The source audit did not compile, run, benchmark or edit the archive.

| Evidence in latest snapshot | Consequence for extension |
|---|---|
| `radiation_closure`, `radiation_transport`, `radiation_runtime`, `radiation_cosmology`, `radiation_amr` | Reuse multigroup M1, HLL/GLF, PD-ARS and existing radiation moments; do not replace them wholesale. |
| `radiation_matter_coupling`, `TransportedPhotoRates` and source ledger | Photoionization reconciliation not wired; grey coupling only in admitted reference path. Stellar ionizing sources refused. |
| `radiation_angular.hpp` | Direct rays are steady-state; VET tensor construction is not a VET transport backend; Monte Carlo supports isotropic elastic scattering only. |
| RT capability admission | Single-rank, level-synchronous M1; MPI, dynamic regrid and non-M1 production transport remain refused. |
| MHD pass 02 and `runtime_capabilities.cpp` | Restricted synchronized CPU leaf-AMR dispatch, with halo/EMF exchange in source; dynamic magnetic regrid, recursive subcycling, GPU and chemistry-derived coefficients blocked. |
| BH/AGN source engine and transactions | Reuse reservoir/budget and source transaction contracts; source populations and exchange pathways require qualification. |
| Restart v27; radiation science export incomplete | Add versioned scientific fields and restart comparison before relying on downstream observables. |
| RT ledger PERF01 | RT work copies estimated ~5 times canonical radiation size; scratch not fully governed. This is a source estimate, not an RSS measurement. |
| Neutrinos, CR, radio, Ly-alpha, CMB | Dedicated production modules not established by inspected evidence; CR excluded from latest RT campaign. |

Foundation closure SHALL prioritize governed RT workspaces; bounded MPI halos; conservative regrid/migration and reflux; stellar source luminosities; radiation/chemistry reconciliation; source and energy ledgers; scalable topology construction; science-snapshot export; MHD synchronization and qualified CT transfer. Do not lift a capability refusal merely to expose a configuration option. The existing source-only campaign policy permits implementation and test-source authoring before the separate owner qualification campaign.

## 3. System architecture and authority boundaries

The core integrator owns time, geometry, canonical particle/cell state, gravity, hydro and MHD. Cosmology owns expansion and species accounting. Radiation and CR solvers own their own evolved moments; chemistry owns species abundances; source populations own their emission and injection budgets. Observables consume immutable state views or versioned lightcone records. Instruments consume skies, never mutate astrophysical state.

| Layer | Modules and data flow |
|---|---|
| Cosmological sector | Boltzmann tables → ICs, H(a), species response → gravity and expansion |
| Evolving matter | TreePM ↔ hydro/AMR/MHD ↔ chemistry, stars and BH/AGN |
| Transport | Radiation and CR moments ↔ conservative matter/source transactions |
| Geometry & spectral services | AMR leaf traversal, lightcones, atomic data, frequency/momentum grids |
| Observables | Ly-alpha, HI, continuum/Stokes, CMB, DM/RM and statistics |
| Instruments | Sky cube → beam/Jones model → visibilities or detector data → reconstruction |

One physical quantity must have one canonical authority. The same hydrogen and electron fields SHALL feed Ly-alpha, HI, SZ, DM and free-free calculations; the same magnetic field SHALL feed CR transport and Faraday synthesis; the same electron spectrum SHALL feed synchrotron and inverse Compton emission. Derived caches may not become second authorities. Share geometry, units, interpolation, spectral quadrature, RNG, transport primitives and resource management, not an inappropriate universal solver for all particles and photons.

Required common contracts: `CosmologicalBackground`, radiation-background authority, versioned atomic/rate tables, spectral discretization descriptors, emissivity/opacity providers, source transaction and conservation ledgers, leaf-geometry/ray interfaces, read-only science-state views, restart codecs, capability admission and memory leases. Names are illustrative until reconciled with repository interfaces. Use the repository's current `.cpp/.hpp` convention, lower_snake names, explicit frame/unit suffixes, and stable external HDF5 schema names.

## 4. Massive neutrinos and cosmological radiation

**Decision N01:** retain background-only, instantaneous transfer-function correction, time-dependent linear response, conventional particles and delta-f weighted particles as distinct models. Prefer linear response for efficient cosmology and delta-f for nonlinear halo/cluster science. A momentum-shell hierarchy and hybrid phase-space model are future research.

For scale factor a, the expansion includes photons, baryons, cold dark matter, curvature, dark energy and the neutrino energy density converted to mass-equivalent units:

$$H^2/H_0^2=\Omega_\gamma a^{-4}+(\Omega_b+\Omega_c)a^{-3}+\Omega_k a^{-2}+\Omega_\Lambda+\rho_\nu(a)/\rho_{\mathrm{crit},0}.$$

Integrate the Fermi–Dirac momentum distribution per mass eigenstate; tabulate H(a), densities and pressure over the relevant epoch. Treat massless radiation consistently and prevent double-counting with Omega_m or Omega_r. IC transfer functions must specify gauge, species, normalization, k units, scale factor and velocity convention. CLASS/CAMB adapters shall preserve table provenance and reject unsupported extrapolation.

For the nonlinear reference, use collisionless phase-space evolution with f = f_ref + delta-f. Analytic reference density plus signed particle perturbation weights reduce variance compared with sampling the full thermal distribution. Derive relativistically consistent drift, kick, weight and deposition formulas for the supported weak-field/Newtonian regime; do not advertise full relativistic perturbation evolution. Signed weights need special care: multipoles, centers of mass and opening criteria cannot assume positive masses. PM deposition and tree force splitting must not double-count or omit species contributions. A PM-only neutrino force model is allowed only as an explicit scale-limited approximation.

| Method | Strength | Limitation |
|---|---|---|
| Background | Negligible cost; correct homogeneous transition | No density perturbations |
| Instantaneous transfer ratio | Cheap, no sampling noise | Prescribed correlation with cold matter |
| Linear response | Responds to evolving potential, modest memory | History kernels; nonlinear phase-space limitations |
| Ordinary particles | Independent nonlinear comparator | Thermal shot noise and large population |
| Delta-f | Variance reduction, nonlinear clustering | Signed weights, weight evolution, specialized gravity |
| Shell/hybrid | Potentially richer phase space per byte | Representation interfaces and fast crossing streams |

Outputs: H(a), species densities, weights/noise estimates, P_cc, P_nunu, cross spectra and total matter P(k), halo neutrino profiles and response convergence. Validate against CLASS/CAMB, SWIFT and matched particle calculations; include massless/zero-mass limits, free-streaming, gauge consistency, momentum sampling, weight conservation where applicable, force accuracy and rank/restart equivalence. Do not infer nonlinear halo accuracy from linear power-spectrum agreement alone.

## 5. Multigroup radiation hydrodynamics

**Decision R01:** retain and finish the current multigroup M1 + PD-ARS IMEX implementation; add independent VET and adaptive direct-ray transport backends. FLD remains a specialized optically thick comparator; Monte Carlo is a separate reference and spectral transport engine. Never route a VET tensor through an M1-specific Riemann solver without a new derivation.

For each frequency group, schematic proper-frame transport is

$$\partial_t E_g+\nabla\cdot\mathbf F_g=S_{E,g},\qquad \partial_t\mathbf F_g+c^2\nabla\cdot\mathbb P_g=\mathbf S_{F,g}.$$

M1 sets the pressure tensor from local energy and reduced flux; VET obtains it from a nonlocal angular transport solution. The detailed specification SHALL derive the actual comoving variables, expansion work, redshifting across frequency bins, frame transformations, physical/reduced light speeds and source exchange. Energy and photon-number accounting must be consistent with group-integrated cross sections and the selected moment variables. An ionizing photon budget is not inferred from bolometric energy without the assumed spectrum.

Start with HI (13.6–24.6 eV), HeI (24.6–54.4 eV) and HeII-ionizing groups. Provide extensible FUV/Lyman-Werner, optical, IR, X-ray and hard AGN groups, with nonoverlapping edges, cross sections, spectral assumptions and redshift remapping. The ionizing baseline is not a complete spectrum; Ly-alpha resonant transport remains distinct. Stellar population/AGN emission tables, radiation pressure, photoheating, H/He non-equilibrium chemistry, shielding and dust extensions must share source accounting. Explicit UV backgrounds and transported photons cannot be counted twice.

| Method | Advantages | Risks/cost |
|---|---|---|
| FLD | Small state; useful in diffusion limit | Poor optically thin anisotropy/shadows |
| M1 | Local, source-count independent, AMR-friendly | Crossing beams merge artificially; closure error |
| VET | Better angular geometry and multiple sources | Nonlocal angular solve, iteration and MPI cost |
| Adaptive rays | Accurate direct illumination, shadows | Source count, migration, diffuse scattering; steady-state caveat |
| Discrete ordinates | Angularly explicit reference | Ray effects, storage and angular convergence |
| Monte Carlo | General scattering and spectra | Noise, thick-medium expense, irregular workloads |
| Hybrid direct/diffuse | Directional bright sources plus cheap diffuse field | Transfer accounting, interfaces, double counting |

Time integration: preserve explicit physical-c comparator, radiation subcycling and implicit stiff matter coupling. IMEX is preferred; validate order, stiff limit and positivity/realizability (E >= 0, |F| <= appropriate transport speed times E). Fully implicit and asymptotic-preserving alternatives are advanced methods, not presumed qualified by PD-ARS alone. Reduced c is explicit, regime-limited and convergence-tested; never a hidden optimization. Current cosmological reduced-c refusal remains until a consistent contract exists.

AMR requires realizable conservative prolongation/restriction, per-group reflux, synchronized source updates and bounded halo exchange. Regridding must reserve old/new coexistence, preserve radiation and CT magnetic authority and restart at a supported synchronization point. VET needs an angular solve, closure-update cadence, boundary conditions and a compatible transport integrator. Causal ray transport must retain photon travel time; the existing steady-state ray routine is only applicable when retardation is negligible. Hybrid direct/diffuse transport must ledger emission, absorption, scattering, escape and unresolved queues exactly once.

Validation: free streaming, crossing beams, shadows, diffusion limit, Strömgren sphere, moving ionization front, photoheating, radiation pressure, cosmological redshift, photon conservation, AMR interface and MPI equivalence, source-off regression and stiff-source convergence. Compare with RAMSES-RT, AREPO-RT/THESAN, Quokka and VETTAM under matched assumptions.

## 6. Cosmic-ray magnetohydrodynamics

**Decision CR01:** design for separate proton and electron spectra with two-moment transport from the outset; implement grey and one-moment comparators as first usable slices. The advanced parity target is conservative, spectrally resolved two-moment CR-MHD using piecewise power-law momentum groups. Do not equate a proton-pressure model with a predictive electron synchrotron model.

A schematic grey energy equation is

$$\partial_t e_{\mathrm{cr}}+\nabla\cdot(e_{\mathrm{cr}}\mathbf u)=-P_{\mathrm{cr}}\nabla\cdot\mathbf u-\nabla\cdot\mathbf F_{\mathrm{tr}}+Q_{\mathrm{cr}}-\Lambda_{\mathrm{cr}}.$$

The pressure force on gas, CR work, streaming heating, injection and losses must form a consistent exchange ledger. The ultrarelativistic grey closure P=e/3 is not valid for every spectral distribution. Compute pressure from the spectrum for general populations. Diffusion is anisotropic along the CT magnetic field, with parallel/perpendicular coefficients and a defined weak-field limit. A two-moment model evolves independent fluxes and finite-speed propagation under a declared scattering closure. It does not solve microscopic plasma uncertainties by itself.

Support advection, field-aligned diffusion, streaming, two-moment transport and energy-dependent coefficients. Explicit diffusion has a fine-grid timestep proportional to dx^2/kappa; super-time-stepping, implicit or IMEX backends must be checked for monotonicity, positivity, conservation and MPI scaling. Streaming must account for wave-mediated energy transfer without counting heating twice. Reduced signal speed must be documented and tested.

| Representation | Advantages | Limitations |
|---|---|---|
| Grey / separate grey species | Affordable dynamical backreaction | No detailed spectral ageing; closure assumptions |
| Fixed bins | Transparent reference | Large uniform momentum storage |
| Piecewise power-law | Rich spectra with fewer groups | Reconstruction, cutoffs, number/energy consistency |
| Tracer spectra | Affordable electron histories | Gas tracers alone do not model relative CR transport |
| Spectral two-moment | Energy-dependent propagation, pressure and emission | High storage, stiff coupling, closure uncertainty |
| Kinetic/PIC | Microscopic physics in resolved domains | Not viable as a global cosmological solver |

Use physically specified momentum bounds and distributions, conservative momentum-space fluxes and explicit boundary losses. Evolve number and energy per bin (and corresponding vector fluxes in the full two-moment variant); reconstruct slopes only from admissible moments. Injection providers for SNe, shocks and AGN must debit their budgets. Proton losses: Coulomb/hadronic; electron losses: Coulomb, synchrotron, inverse Compton and bremsstrahlung. Include secondary e+/e- production, CR heating/ionization and gamma-ray diagnostics as qualified models. Separate energy deposited in thermal matter from energy escaping as radiation/secondaries.

Reuse MHD field authority, AMR and source transactions; no CR state on DM-only runs. Reflux spatial transport and momentum-group interfaces, preserve nonnegative spectra, test energy/number accounting including open boundaries and losses. Compare grey and spectral solutions with AREPO CR, FIRE, RAMSES-MCR, CRESCENDO and CREST as applicable. Qualification cases: advection, anisotropic diffusion, streaming equilibrium/heating, CR shocks, injection/loss spectra, spectral ageing, secondary production, AMR/MPI and joint gas+CR energy closure. Spectral adaptivity is deferred until uniform-bin and piecewise power-law solutions qualify.

## 7. Lyman-alpha physics

**Decision L01:** deterministic Voigt-profile absorption and a distinct frequency-dependent resonant Monte Carlo engine. Three use cases are independently admitted: forest absorption, resonant emission/escape, and cosmic-dawn scattering-rate/feedback coupling.

For absorption, integrate optical depth tau(nu) = integral n_HI sigma_alpha(nu_local,T) ds_phys and F = exp(-tau). Use Voigt profiles, proper column density, thermal broadening, peculiar velocity, Hubble/redshift mapping and adaptive sightline integration. Use stable profile evaluation and quadrature with damping-wing coverage. Produce mock spectra, column densities, flux PDF, mean transmitted flux and 1D/3D flux spectra. Doppler-only and fluctuating Gunn–Peterson approximations are named cheap comparators, not defaults for damped systems. Diffuse IGM resolution is essential; halo-centric AMR alone does not guarantee forest convergence.

For resonant emission, sample photons from recombination and collisional excitation, then evolve frequency/direction through atom-frame redistribution, thermal velocities, recoil where required, dust absorption, cosmological redshifting and appropriate phase functions. Reuse packet IDs, geometry and RNG infrastructure but not the existing isotropic elastic scattering kernel as physical Ly-alpha scattering. Produce escaped spectra, maps, escape fractions, scattering rates and momentum deposition. Unaccelerated MC is a reference; core-skipping or diffusion acceleration requires observable-specific equivalence. A method valid for escape spectra may be invalid for scattering counts, Wouthuysen–Field coupling or force estimates. Packet truncation and unresolved energy must be reported with error bounds, not hidden.

Validation: analytic slab optical depths and equivalent widths, thermal profiles, damping wings, moving gas, homogeneous sphere/slab resonant escape, frequency redistribution, dust, photon budgets, packet-number convergence and rank reproducibility. Compare forest statistics to Nyx and resonant results to COLT/RASCAS. Coupled Ly-alpha force or spin-temperature calculations require separate estimator validation.

## 8. Cosmic microwave background

**Decision C01:** use CLASS and CAMB for primary perturbations, recombination and primary CMB spectra; native late-time physics and secondary observables. A native Einstein–Boltzmann engine is separate long-term research, not a modification to Newtonian TreePM.

A shared radiation-background authority supplies T_CMB(a)=T_CMB,0/a and u_CMB=a_rad T_CMB^4 under standard adiabatic evolution. It feeds Compton thermal exchange, inverse Compton CR losses, radiative excitation and HI spin coupling. Preserve the background versus perturbation distinction and avoid independent temperatures in modules.

| Observable | Baseline and preferred extension |
|---|---|
| Primary TT/TE/EE and lensing inputs | External Boltzmann tables, consistent cosmology and gauge |
| CMB lensing | Born projection comparator; multi-plane ray tracing with lens-lens coupling |
| Thermal SZ | Electron-pressure integral; relativistic spectral corrections in hot gas |
| Kinetic SZ | Optical-depth-weighted electron LOS momentum |
| ISW/Rees–Sciama | Time-varying potentials with controlled temporal interpolation |

The nonrelativistic Compton-y is sigma_T/(m_e c^2) times the proper-path integral of electron pressure. kSZ requires consistent LOS velocity sign, electron density, optical depth and observer frame. Multi-plane lensing must define distances, projected potential normalization, plane thickness, deflection interpolation and boundary treatment. Use consistent primary realizations and lightcone geometry for cross-correlations.

Small periodic volumes do not contain all relevant long modes, and replicated boxes can produce artificial correlations. Report volume, redshift coverage, replication strategy, missing-mode caveats and angular-scale validity. Qualification: homogeneous and analytic lens/cluster projections, Born versus multi-plane convergence, relativistic SZ kernels, uniform electron slabs, redshift/velocity signs, temporal interpolation and matched cross-spectra against Agora or other reference lightcones. Do not claim full-sky cosmological accuracy from a small repeated box.

## 9. Radio astrophysics

**Decision RA01:** native frequency-dependent, full-optical-depth, full-Stokes radiative transfer, with composable emission/absorption/propagation providers. Start with isolated-galaxy HI cubes, then add cosmic-dawn/reionization, synchrotron and polarization once physical dependencies qualify.

The HI rest frequency is approximately 1420.40575 MHz; observed frequency is rest frequency/(1+z). Under quasi-equilibrium, inverse spin temperature is the weighted sum of inverse radiation, kinetic and Ly-alpha color temperatures with collisional and Wouthuysen–Field couplings. A homogeneous segment against a background has observed differential brightness temperature (T_s−T_gamma)(1−exp(−tau_21))/(1+z). For nonuniform gas, perform ordered radiative transfer; do not apply the homogeneous formula to an entire arbitrary sightline. Saturated spin temperature and optically thin assumptions are named regime-limited shortcuts, not universal defaults.

For cosmic dawn, support Ly-alpha coupling, X-ray heating, collisions, neutral fraction and velocity effects. Semi-numerical reionization (21cmFAST-like) is a useful independently classified large-volume model; photon conservation, source assumptions and resolution dependence must be audited. Coupled RHD remains the higher-fidelity target for resolved source/gas interactions. Never imply a semi-numerical field is equivalent to full coupled evolution without evidence.

For polarized transfer evolve Stokes S=(I,Q,U,V) along the ray:

$$\frac{d\mathbf S_\nu}{ds}=\mathbf j_\nu-\mathbf K_\nu\mathbf S_\nu.$$

For constant coefficients, use a matrix-exponential formal solution and stable evaluation of the source integral even for singular K. DELO-type methods and adaptive quadrature are selectable alternatives. Define polarization basis, handedness, frequency, cosmological intensity transformation, Faraday rotation/conversion conventions, and physical realizability. Handle rapid low-frequency rotation and bandwidth depolarization. A scalar or screen solution is an independent comparator.

| Physics provider | Requirements and caveats |
|---|---|
| HI 21 cm | Neutral H, spin level populations, broadening, velocity and background |
| Synchrotron/SSA | Electron spectra and pitch-angle model, magnetic geometry; not total CR energy alone |
| Thermal free-free | Electron/ion populations and temperature, Gaunt factors |
| Faraday rotation/conversion | Magnetic LOS field, plasma distributions and polarization conventions |
| Dispersion/rotation measure | Proper path, redshift weighting, electron density and B; optional FRB source models |
| Recombination/molecular lines | Atomic/molecular level populations and excitation; later explicit providers |

Keep power-law and steady-state electron approximations available for diagnostics. Time-dependent tracer and Eulerian electron spectra are advanced parity paths. Generate brightness temperature, intensity, spectral cubes, Stokes maps, spectral-index maps, DM/RM and correlations. Molecular lines, masers, cyclotron and detailed plasma effects require separate physical specifications; they are not automatically implied by a full-Stokes formal solver.

Validate uniform slabs, HI thin/thick limits, velocity broadening, spin-temperature limiting cases, synchrotron spectra/ageing, free-free emission/absorption, analytic Faraday screens, constant-matrix solutions, spectral/channel convergence and comparison with POLARIS/21cmFAST where applicable. Separate numerical transport error from uncertain CR acceleration, magnetic and ionization models.

## 10. Shared lightcone and virtual observatory

**Decision O01:** implement one read-only lightcone/ray service for observables and a standalone instrument layer. The native reference is direct full-polarization RIME; production backends may use GPU direct summation, NUFFT or antenna-based matrix multiplication. Initially interoperate with OSKAR, pyuvsim, CASA and WSClean rather than duplicating every mature telescope tool.

The lightcone service SHALL define observer position, orientation, cosmology, emission time, radial shell boundaries, spatial/temporal interpolation, periodic wrapping, replication and missing-mode metadata. Traverse AMR leaves once per ray segment; never double count covered parents. Distinguish snapshot rays, evolving lightcones and gravitationally deflected paths. Provide streaming/tiled accumulation, bounded MPI exchange, shell checksums and versioned science products. Verify units, cosmological redshifting, LOS signs and interpolation convergence independently of the emissivity model.

A schematic full-polarization measurement equation is

$$\mathbf V_{pq}=\int \mathbf J_p(\hat s)\mathbf B(\hat s)\mathbf J_q^\dagger(\hat s)e^{-2\pi i\mathbf u_{pq}\cdot(\hat s-\hat s_0)}d\Omega.$$

Define sky coherency/Stokes normalization, phase center, baseline and Fourier signs, coordinate frames, polarization basis, full wide-field w-term, integration bandwidth/time, beam normalization and units. Add station positions, Earth rotation, analytic/measured/EM-derived beams, gains, phases, delays, bandpass, leakage, thermal noise, ionospheric phase/Faraday effects, flagging and calibration error. Foreground and instrumental chromaticity must be modelled separately from the desired cosmological 21 cm signal.

| Backend | Advantages and limitations |
|---|---|
| Direct source sum/RIME | Transparent discretized reference; cost O(N_baseline N_source N_time N_freq) |
| FFT / NUFFT | Fast; interpolation and direction-dependent corrections need quantified error |
| Antenna matrix methods | GPU-friendly; batch/working-set constraints |
| External OSKAR / pyuvsim | Mature independent reference and interoperability; external dependencies |
| Statistical receiver noise | Efficient integrated products; misses detailed waveform artifacts |
| Complex baseband | Quantization/filtering/correlator studies; enormous sample data, specialized opt-in |

Native products: FITS with WCS, spectral and Stokes axes; HEALPix where appropriate; UVH5 and Measurement Sets through tested adapters; mock spectra and calibrated maps. Include physical simulation and instrument configuration hashes, seeds, beam/calibration model and selection mask. For N_a antennas the cross-baseline count is N_a(N_a−1)/2. Channel/time/source batching and streamed output are mandatory for large arrays. Validate point-source phases, analytic beams, polarization leakage, noise variance, time/bandwidth smearing, wide-field geometry and flux recovery; compare complex visibilities with external reference tools before using accelerated modes in precision inference.

## 11. Coupling, scheduling and numerical invariants

The runtime registry must declare dependencies, read/write sets, frame, time interval, memory contract, retry behavior and capability status. Each enabled physics configuration must be a legal composition, not an arbitrary flag combination. Reuse `SourceTransactionBatch` and existing proper-frame/source ledgers where applicable. Specify coupling order in detailed subsystem designs; the following is a logical sequence, not a claim of second-order splitting: synchronize geometry and clocks; prepare gravity and matter state; evolve hydro/MHD and transport stages; reconcile stiff exchanges and sources; reflux/restrict and commit; then checkpoint or sample observables.

SHALL enforce: (1) no duplicated source injection or background absorption; (2) gas+CR+radiation+source budgets, with explicit expansion work, escape and losses; (3) baryon/charge/species accounting where applicable; (4) positivity and radiation/CR moment realizability; (5) magnetic divergence/CT authority; (6) conserved AMR transfers and MPI ownership; (7) deterministic or tolerance-defined restart/retry semantics. Global energy conservation in expanding cosmology is not a naive constant-energy invariant. Frame conversion and expansion terms belong in the ledger. Any rejected step must roll back all affected authoritative fields, reservoirs, counters and RNG state or fail safely to a known checkpoint.

## 12. Memory, performance and hardware profiles

Persistent optional state SHALL be lazy and owned by the relevant species, patch or module. Follow the MemoryGovernor doctrine: reserve large allocations before commit, report capacities/high-water and rank maxima, bound communication, reuse phase arenas, stream analysis, avoid duplicated state and fail before OOM. Memory pressure may reduce batches or postpone optional work; it may not silently reduce fidelity or disable required physics. Any validated alternate algorithm must be user-authorized and provenance-visible.

For N cubed cells and n FP64 scalar fields, canonical additional state is 8 n N^3 bytes. Illustrative global totals below exclude parent/ghost cells, transport stages, exchange, caches, metadata, old/new regridding, I/O, existing gravity/hydro/MHD and allocator overhead.

| Additional state | 128^3 | 256^3 | 512^3 |
|---|---:|---:|---:|
| One scalar | 16 MiB | 128 MiB | 1 GiB |
| Grey CR E + 3 flux components | 64 MiB | 0.5 GiB | 4 GiB |
| M1, 3 groups (12 scalars) | 192 MiB | 1.5 GiB | 12 GiB |
| M1, 6 groups (24 scalars) | 384 MiB | 3 GiB | 24 GiB |
| 2 CR species × 16 bins × 2 scalars | 1 GiB | 8 GiB | 64 GiB |
| 2 CR species × 16 bins × 8 scalars | 4 GiB | 32 GiB | 256 GiB |

The last row is the illustrative full two-moment number/energy plus two vector fluxes model, not a universal storage requirement. Avoid duplicating grey fields when derived from spectral moments. A 256^3 spectral two-moment example already costs 32 GiB before all other physics. The existing 48 GiB policy proposes a 42 GiB simulation ceiling and 36 GiB nominal target; actual hardware reservation is configurable and measured. Thus full uniform spectral CR+RHD is not a credible default on 48 GiB. Zooms, sparse physical state where justified, coarser spectral models and distributed storage are necessary choices, not silent degradation.

| Target tier | Intended use, not demonstrated capacity |
|---|---|
| 32–48 GiB desktop | Small reference decks, focused zooms, grey CR, streamed observables |
| 64–128 GiB workstation | Coupled MHD/RHD and CR at appropriate resolution |
| 128–512 GiB node | Richer groups, larger zooms and spectral studies |
| 0.5–8 TiB aggregate HPC | Distributed coupled volumes and ensemble campaigns |

Require module-specific byte/element and peak-live-set models; patch-local SoA, hot/cold split, bounded packet/ray buffers, scalable topology indices, shared read-only rate tables, CPU/OpenMP first and accelerator-specific host/device leases. GPU kernels must preserve reference semantics; FP32/mixed precision is a dedicated validation decision. Measure wall time, rank imbalance, memory bandwidth, bytes communicated, packet variance and accuracy per unit resource. No hardware-fit guarantee is made here.

## 13. Parity and qualification matrix

Feature presence, successful compilation, numerical verification, physical validation, calibration, performance qualification and production admission are distinct states. Maintain a per-feature ledger with exact source commit, config, test, evidence and failure reason. Source-only campaigns shall author all required tests and fail-closed guards but not report unexecuted tests as passed. Qualification can be staged by dependency during a dedicated execution campaign; a feature cannot be production-qualified merely because downstream source work is complete.

| Domain | Independent references | Core parity measurements |
|---|---|---|
| Cosmology/gravity | CLASS/CAMB, SWIFT, GADGET-4/AREPO | H(a), scale-dependent growth, species P(k), halo profiles |
| Hydro/AMR/MHD | RAMSES, Nyx, Athena++ | Shock/wave convergence, conservative fluxes, divergence |
| RHD | RAMSES-RT, AREPO-RT, Quokka, VETTAM | Fronts, shadows, energy/photon balance, stiff limits |
| CR | AREPO, FIRE, RAMSES-MCR, CRESCENDO/CREST | Transport, heating, spectral and pressure evolution |
| Ly-alpha | Nyx, RASCAS, COLT | Forest statistics, resonant escape, scattering estimators |
| Radio | 21cmFAST, POLARIS | HI brightness and spectra, Stokes, RM |
| CMB | CLASS/CAMB, Agora | Primary spectra, lensing/SZ, correlations |
| Instrument | OSKAR, pyuvsim | Complex visibility errors and calibrated recovery |

Proposed initial engineering targets for detailed gate design: relative background error 10^-6 on a declared range; neutrino-induced P(k) change to 1% on converged scales (define denominator and absolute floor); analytic ionization-front radius and CR profiles to 1% under adequate resolution; deterministic optical-depth and analytic CMB kernels to 0.1% where well-conditioned. Monte Carlo gates must report confidence intervals and systematic errors; visibility tolerance depends on science case/foreground leakage. These are discussion-derived candidate thresholds, not published universal standards or achieved metrics. Set each norm, reference, domain, resolution ladder, statistical confidence and near-zero handling before runs.

All modules require disabled-feature regression, conservation/accounting tests, analytic limiting cases, resolution/time/group convergence, MPI/rank and restart equivalence, fault/rollback tests, resource admission/high-water tests and matched cross-code comparisons. Scientific parity applies to matched problems, not an assertion of parity with every feature of every code.

## 14. Repository organization and configuration

Proposed directories (adapt to existing layout rather than move working code unnecessarily): `include/cosmosim/{neutrinos,radiation,cosmic_rays,lightcone,observables,instruments}/` and corresponding `src/` modules; observable subdirectories `lyman_alpha`, `cmb`, `radio`; instrument subdirectories `interferometry`, `single_dish`, `receiver`. Each module owns typed configuration, kernels, storage, runtime adapter, memory contract, restart/diagnostics, unit/integration tests and benchmarks. Preserve existing radiation, MHD, source and BH/AGN interfaces where possible.

Illustrative configuration only (not current accepted keys or a runnable deck):

```yaml
neutrinos:
  method: delta_f  # background, transfer_function, linear_response, particle
radiation:
  transport: m1   # fld, vet, ray_tracing, hybrid
  integrator: imex
cosmic_rays:
  transport: two_moment  # advection, diffusion_streaming
  spectrum: piecewise_power_law  # grey, fixed_bins, adaptive_bins
  species: [proton, electron]
observables:
  lyman_alpha: voigt_absorption
  cmb_lensing: multi_plane
  radio: full_stokes
instruments:
  backend: external_oskar  # direct_rime, nufft, matrix
```

Validate inactive keys, incompatible choices, supported modes and memory before allocation. Reject missing physics dependencies, e.g. polarized synchrotron with neither B nor a declared electron distribution; scattering-coupled spin temperature without a qualified coupling provider; unsupported neutrino force/gauge regimes; two-moment CR without its coupling contract. Store actual resolved config, data-table versions/hashes, commit, compiler, precision, rank/thread layout, RNG policy, assumptions and capability statuses. Restart schema changes must be versioned with explicit compatibility rules; no silent loss of canonical state.

## 15. Implementation roadmap and decision gates

| Phase | Deliverable and exit condition |
|---|---|
| 0 — Foundation closure | Finish RT/BH/AGN/MHD integration, chemistry, memory, MPI/regrid and export; preserve refusals for unsupported modes. |
| 1 — Shared foundations | Cosmology/Boltzmann tables, neutrino background and linear response, atomic/spectral services, lightcone/ray geometry and science schema. |
| 2 — First observables | Voigt forest, HI cubes, Born lensing/SZ, DM, external interferometry; analytic and cross-code test sources. |
| 3 — Coupled physics | Delta-f particles; complete M1/IMEX source coupling; grey CR, anisotropic diffusion/streaming and two-moment infrastructure; VET/ray backend. |
| 4 — Advanced parity | Spectral proton/electron transport; resonant Ly-alpha; cosmic dawn; relativistic SZ/multi-plane lensing; full Stokes and native RIME. |
| 5 — Qualification | Owner-controlled compile, tests, convergence, benchmark and production admission campaigns. Staged dependency qualification is permitted. |
| 6 — Frontier research | Compare new methods with retained qualified references; promote only after error/resource evidence. |

Order is a proposed dependency plan, not an estimate of delivery dates. Specifications can proceed in parallel; source work must not remove existing behavior or bypass blocked dependencies. First flagship campaigns: massive-neutrino matter suppression and halo profiles; isolated-galaxy HI cube to visibilities; IGM forest with matched neutrino cosmologies; CR-MHD galaxy winds and spectral radio maps; reionization/21 cm coupled lightcones; cluster lensing and SZ. Select volume/resolution to fit the question and budget.

## 16. Frontier research: adaptive multiphysics fidelity

After parity, investigate adapting the physical discretization itself, beyond spatial AMR: spatial and temporal resolution, momentum groups, frequency groups, angular closure/quadrature, chemistry and gravity tolerances. The goal is to minimize science error subject to hard memory and wall-time budgets and physical invariants. The error vector must distinguish discretization, closure/model uncertainty, sampling and finite-volume effects; a local indicator is not automatically a reliable global observable bound.

| Research direction | Potential benefit | Required proof / main risk |
|---|---|---|
| Spectral CR AMR | Resolve cooling breaks and injection only where necessary | Conservative number/energy/flux remapping, positivity, MPI interfaces |
| Adaptive radiation groups | More spectral fidelity per byte | Photon/energy redistribution, group-dependent cross sections, redshift |
| M1/VET/direct-diffuse hybrid | Correct angular complexity where needed | Causal closure, interface fluxes, no double counting |
| Neutrino phase-space hybrids | Low noise with selective nonlinear resolution | Fast crossing streams, signed perturbations, consistent force split |
| Ly-alpha MC/diffusion | Thick-region acceleration | Frequency redistribution, scattering counts, forces and escape errors |
| Adaptive Stokes spectral rays | Avoid oversampling smooth spectra | Narrow lines, rapid Faraday rotation, bandwidth error |
| Observable-aware AMR/error estimation | Better science return per resource | Bias to selected observables; independent holdout tests |

Start research with spectral adaptation and adaptive polarized quadrature, after their uniform references qualify. Require conservative common refinement at heterogeneous interfaces, hysteresis, error estimates, memory reservation for transitions, restart persistence and reproducible decisions. Compare against globally rich reference runs and independent observables; record cases where adaptivity worsens accuracy or wall time. No claim of originality until a dedicated literature audit. Native early-Universe GR/Einstein–Boltzmann evolution and microscopic plasma/kinetic solvers remain separately scoped research programs.

## 17. Specification register and source basis

This document is the umbrella design. Next detailed specifications (version 0.1 unless amended):

1. `chui_neutrino_cosmology_specification_v01.md`
2. `chui_radiation_hydrodynamics_specification_v01.md` (extension/amendment to existing RT design)
3. `chui_cosmic_ray_transport_specification_v01.md`
4. `chui_lightcone_ray_engine_specification_v01.md`
5. `chui_lyman_alpha_specification_v01.md`
6. `chui_cmb_physics_observables_specification_v01.md`
7. `chui_radio_astrophysics_specification_v01.md`
8. `chui_virtual_observatory_specification_v01.md`
9. `chui_extended_physics_memory_performance_specification_v01.md`
10. `chui_multiphysics_parity_validation_specification_v01.md`
11. `chui_adaptive_multiphysics_fidelity_research_dossier_v01.md`

Internal sources: CosmoSim Master Design Document (11 March 2026), §§2–4, 7–13; CHUI Memory Architecture and Resource Governance Specification v0.1 (23 August 2026), §§2, 4–17; AREPO Lessons technical summary; Master Variable and File Naming Conventions v1.0; latest supplied archive, especially `docs/architecture/chui_rt_bh_agn_implementation_pass01_status.md`, `chui_mhd_implementation_pass02_status.md`, `chui_rt_bh_agn_master_design_v01.md`, `src/workflows/runtime_capabilities.cpp` and `include/cosmosim/radiation/radiation_angular.hpp`; design discussions of 8–27 November 2026. The archive ledger's base commit is not asserted to be the current archive HEAD.

External comparison targets recorded from the discussion: CLASS, CAMB, SWIFT, GADGET-4, AREPO, RAMSES-RT, THESAN, Quokka, VETTAM, Nyx, RAMSES-MCR, FIRE, CRESCENDO, CREST, COLT, RASCAS, 21cmFAST, POLARIS, Agora, OSKAR, pyuvsim, CASA, WSClean, fftvis and matvis. These names are reference leads, not a fresh exhaustive literature review or an assertion that each supports all methods listed. Detailed specifications must cite exact papers, versions, equations and benchmark datasets before implementation approval.

**Open decisions:** exact species/force implementation, neutrino IC gauge and relativistic regime, radiation group boundaries and closures, CR scattering and spectral model, angular solver and coupling scheme, lightcone schema, supported instrument formats, memory admission coefficients, scientific tolerance norms and calibration priors. Equations above are schematic where noted; detailed comoving derivations and discrete conservation proofs are required before numerical implementation. No code was changed and no simulation executed for this document.

**Closing doctrine:** establish credible established physics first; preserve the comparators; optimize without corrupting scientific meaning; innovate only when a better physical answer can be demonstrated.
