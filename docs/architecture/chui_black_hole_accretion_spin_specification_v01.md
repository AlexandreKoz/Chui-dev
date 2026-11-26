# CHUÍ — Black-Hole Lifecycle, Accretion and Spin Specification

**Version:** 0.1 · **Date:** 2026-11-26 · **Status:** proposed reference design  
**Parent:** [Master design](chui_rt_bh_agn_master_design_v01.md)

## 1. Model architecture and compatibility

Separate seed provider, orbital dynamics, gas-supply estimator, conservative gas withdrawal, disk evolution, horizon evolution, merger policy and emission budget. Legacy `BlackHoleAgnModel` remains available as `legacy_host_cell` comparison with its existing Newtonian mass semantics. New models must not duplicate its callback. P0: halo seed, neighborhood hot/cold supply, friction, bound-pair merger, scalar disk and Eddington-limited engine; P1: local seed channels, resolved supply, vector disk/spin, three accretion states, super-Eddington, binary delays and recoil. Research: radial disk PDE, resolved PN or GRMHD.

## 2. Canonical BH state

Use species-specific SoA indexed by stable particle ID, with row-index lookup cache invalidated on mutation. Persistent core: physical horizon/gravitational mass M_bh, numerical gravitational particle mass and explicit meaning, disk baryonic mass M_d, disk constituent mass bookkeeping (metals/species as required by existing chemistry contract), pending mechanical and spectral energy/mass budgets, cumulative physical loss counters, source time and event sequence. P1 adds spin angular momentum J_bh, disk angular momentum J_d, model state and magnetic-flux proxy if explicitly modeled, binary state and merger lineage. Derived spin a*, luminosity, Eddington ratio, disk radius, accretion rates and kernel environment are caches/diagnostics unless required for restart. Never allocate BH fields per DM/gas particle. Disk and pending ejecta gravitational contributions must be specified to avoid double counting (normally M_dyn=M_bh+M_d+bound pending material, subject to declared mass-energy approximation). Any enhanced numerical mass is a separate labeled comparator, not hidden mass creation. Preserve old output field meaning via versioned schema and explicit export.

## 3. Seeding and lifecycle

### 3.1 P0 halo provider

Production seeding is currently fail-closed in source runtime. Integrate a distributed authoritative halo catalog/candidate provider on declared cadence. Candidate record includes stable halo ID or documented membership identity, total/eligible mass, center and gas donor, existing BH occupancy, owner and time. Select at most one new seed per eligible halo per event under the reference rule; reconcile duplicates across overlapping ranks, halo merger/split and restart. Halo identity persistence and occupied-halo policy must be explicit; never reseed because local index changed. Use existing particle-ID precommit registry. Convert gas to BH with exactly transferred mass, composition and momentum; seed mass cannot exceed available donor mass minus floors. A failed candidate is recorded, not a synthetic mass source. Schedule new BH at next legal tick, invalidate gravity/cache. P1 local seed providers: stellar-remnant, direct-collapse and cluster-collapse criteria depend on existing chemistry/stellar/radiation inputs and separate calibration; no invented astrophysical thresholds.

### 3.2 Orbit and unresolved drag

Use physical peculiar velocities, standard TreePM orbit and appropriate softening. Collisionless Chandrasekhar acceleration for an isotropic Maxwellian background:

$$\mathbf a_{df}=-4\pi G^2 M_{bh}\rho\ln\Lambda[\mathrm{erf}(X)-2Xe^{-X^2}/\sqrt\pi]\mathbf v/v^3,$$
$$X=v/(\sqrt2\sigma).$$

Use stable small-X series (bracket≈4X³/(3 sqrt(pi))) and zero limit at v=0. Separate DM and stellar density/dispersion; gas drag requires a different Mach-dependent closure, not substitution into this formula. Estimate unresolved impact-parameter range with b_min tied to physical strong scattering/softening, b_max limited to unresolved scales and local gradients; set correction zero when range vanishes. This is a proposed resolution-aware policy, not a universally calibrated friction law. Avoid double counting resolved wake. Transfer reaction momentum and work into represented background where feasible with a conservative kernel; if unresolved background reservoir is used, explicitly track its momentum/energy and limitations. No unconditional position snapping in reference. Retain repositioning and boosted dynamical mass as labeled nonphysical comparison options. Limit orbit/drag timestep by velocity change and local dynamical time, validate orbital decay across softening and particle mass ratios.

## 4. Gas supply estimators

Gather only owned leaf cells intersecting compact physical kernel. Use stable IDs, periodic minimum image when applicable, and conservative quadrature weights. Reduce sufficient statistics, not entire neighbor populations. Record kernel radius, resolved cell count/effective mass, uncertainty, phase fractions and inadequacy flags. Radius chosen by explicit target mass/min/max physical radius and resolution; reject invalid/empty support, do not manufacture density. Velocity estimates use relative peculiar velocities, physical density and temperature/sound speed from existing thermochemistry/EOS. Hot and cold components must be mutually exclusive or partition of unity; do not apply both to same mass. Smooth phase weighting is a candidate; version its thresholds/transition and validate.

### 4.1 Bondi and angular-suppressed reference

$$\dot M_{BHL}=4\pi\alpha G^2M_{bh}^2\rho_h/(c_{s,h}^2+v_{rel,h}^2)^{3/2}.$$

Use alpha=1 as unboosted comparator; any boost is an explicit calibrated closure. Denominator degeneracy requires physical EOS floor or rejected estimator, not a magic epsilon. A rotational suppression f_J in [0,1] may be selected via a named published prescription; compare to unmodified Bondi. Compute and record r_B=GM/(c_s²+v_rel²) relative to dx and softening. Do not claim resolved Bondi flow when unresolved. Rates are requested supply, not accepted transfer.

### 4.2 Cold gravitational torque supply

Reference Hopkins–Quataert / SIMBA-inspired model [B1]:

$$\dot M_T=\epsilon_T f_d^{5/2}(M_{bh}/10^8M_\odot)^{1/6}(M_{enc}/10^9M_\odot)(R_0/100\,pc)^{-3/2}(1+f_0/f_{gas})^{-1} M_\odot/yr,$$
$$f_0=0.31 f_d^2(M_d^{gal}/10^9M_\odot)^{-1/3}.$$

Here M_enc is enclosed baryonic mass, M_d^gal is resolved galactic disk mass (NOT the subgrid reservoir), f_d=M_d^gal/M_enc, f_gas is gas fraction of the disk, R_0 is physical estimator radius. Zero eligible disk/gas means zero torque supply. Estimate disk via declared kinematic decomposition of stars and cold gas using local angular momentum and circularity; require sufficient samples and reliable axis, record all measured inputs. Do not conflate thermal phase with rotational support. A precise disk classifier and radius policy must be selected and verified against the published implementation before calibration. SIMBA's normalization includes unresolved retention; CHUÍ must reinterpret or calibrate epsilon_T as *supply* at disk boundary, not reuse retention factor and then independently remove winds. Retain pure published instantaneous comparator with original definitions separately. Hot/cold supply estimator must have well-defined behavior in dispersion-supported nuclei; no silent substitution of a different prescription.

### 4.3 Resolved flux (P1)

$$\dot M_{in}=-\int_{S_r}\rho\min[(\mathbf v-\mathbf v_{surface})\cdot\hat n,0]dA.$$

Use moving control surface in proper coordinates, include Hubble flow if appropriate for the defined physical surface; use finite-volume face flux or qualified embedded quadrature. Distinguish net from inward flux and recrossing gas. Debit only mass actually withdrawn, not merely measured crossing rate. Require convergence versus surface radius, resolution and center motion; flux is nuclear supply, not horizon rate.

## 5. Conservative withdrawal and scalar reservoir

Requested donor mass `dm_req=mdot_supply*dt_phys`. Allocate weights based on eligible mass/phase kernel, clamp each to available mass after floors and other pending transactions; redistribute shortfall deterministically or record unmet supply. Remove mass, momentum, thermal/kinetic energy, species/metals in proportions consistent with canonical state. Do not remove magnetic face flux by mass fraction: MHD sink policy must preserve CT and account for retained magnetic energy. Enforce density/internal-energy floors and gas fraction limits. Update BH+disk center-of-mass motion and capture angular momentum relative to BH, with explicit orbital/internal partition; do not count the same angular momentum as both. Any unresolved binding/thermal energy liberated by withdrawal goes to named reservoir or loss model, not hidden energy.

Scalar reservoir: dM_d/dt=S-(1+eta_w)M_d/t_visc for nonnegative constant accepted supply S, wind loading eta_w and viscous time. Let tau=t_visc/(1+eta_w). Exact update:

$$M_d'=M_d e^{-dt/\tau}+S\tau(1-e^{-dt/\tau}).$$

Integrated drained mass D=M_d+S dt-M_d'; inner inflow D/(1+eta_w), wind eta_w D/(1+eta_w). Use expm1/series to avoid cancellation at small dt/tau; preserve nonnegativity and exact ledger. Eddington cap changes ODE: if horizon inflow capped, excess remains in disk unless an explicitly selected wind removes it. Use piecewise analytic integration across cap threshold or bounded positive implicit ODE, not clip after draining. Time-varying coefficients use substeps with error control. Disk wind energy/momentum require physical budget, not independent from mass loading.

P0 characteristic viscous time configured with explicit physical units. P1 alpha disk t_visc≈[alpha_SS (H/R)^2 Omega_K]^-1, Omega_K=sqrt(GM/R³). For a Keplerian characteristic radius R_c≈(|J_d|/M_d)^2/(GM), valid only when central gravity dominates. Disk self-gravity/Toomre fragmentation, thick-flow geometry, relativistic radius limits require named closure and validity diagnostics. No unbounded reservoir growth disguised as stable accretion.

## 6. Horizon mass, energy and Eddington convention

Two mutually exclusive versioned ledgers:

- `legacy_newtonian`: old `dm_bh=dm_acc` and feedback yield is phenomenological; not a physical rest-energy conservation claim.
- `relativistic_defect`: define dm_0 as baryonic mass reaching inner accretion engine AFTER disk winds. Energy budget includes binding energy released before crossing horizon and extracted spin energy. If epsilon_acc includes all outgoing accretion-powered radiation/mechanical energy relative to dm_0 c² and no other internal retained terms, dm_bh=(1-epsilon_acc)dm_0-E_spin_extracted/c². Define additional retained disk/jet budgets explicitly when needed. Do not subtract wind rest mass twice; its mass was never part of dm_0. Distinguish rest mass of ejecta from mechanical energy above rest mass. Do not double subtract BZ energy if an efficiency fit already includes it. Use actual accepted horizon transfer. Track escaping rest-energy equivalents and physical energy released without injecting all rest mass into Newtonian hydro energy lanes.

Eddington L=4 pi G M m_p c/sigma_T (ionized hydrogen Thomson approximation); define mdot_Edd,ref=L/(epsilon_ref c²), default epsilon_ref=0.1 proposed reference convention. Physical epsilon_r may differ; luminosity cap and mass inflow cap are distinct. P0 reference can cap rate at fixed mdot_Edd,ref, document actual L/L_Edd. P1 slim-disk model supports photon trapping and super-Eddington inflow with a qualified mass/spin-dependent radiative closure. Report both mdot/mdot_Edd,ref and L/L_Edd. No unqualified extrapolation of thin disk efficiency to other regimes.

## 7. Vector disk and spin (P1)

$$\dot{\mathbf J}_d=\dot{\mathbf J}_{supply}-\dot{\mathbf J}_{in}-\dot{\mathbf J}_{wind}-\mathbf T_{align},$$
$$\dot{\mathbf J}_{bh}=\dot{\mathbf J}_{in}+\mathbf T_{align}-\mathbf T_{jet}+\dot{\mathbf J}_{merge}.$$

Internal transfers cancel in the sum; unresolved radiation angular momentum and orbital transfer require their own terms if modeled. Evolve J and M, derive a*=cJ/(GM²); derivative da*/dt=c Jdot/(GM²)-2 a* Mdot/M. Zero-spin axis defined via disk or previous orientation only as documented fallback. Integrate coupled mass/J with adaptive error bounds, exactly equal/opposite alignment increments; use a vector rotation or conservative implicit torque update. Reject |a*|>1 beyond tolerance; a Thorne-like equilibrium cap near 0.998 is a model for photon capture, not numerical clipping. Excess angular momentum/energy must be accounted for in a physical saturation model. Fiacconi and Cenci thin/warped disk closures are references [B2,B3]; implement their specific alignment/counteralignment equations only after pinning model version and reproducing semi-analytic tests. Thick flows require precession/alignment model distinct from thin disks.

### Kerr thin disk

For signed equatorial a relative to disk axis, r_ISCO in GM/c²:

$$Z_1=1+(1-a^2)^{1/3}[(1+a)^{1/3}+(1-a)^{1/3}],\quad Z_2=\sqrt{3a^2+Z_1^2},$$
$$r_{ISCO}=3+Z_2-\mathrm{sgn}(a)\sqrt{(3-Z_1)(3+Z_1+2Z_2)},$$
$$E_{ISCO}=\sqrt{1-2/(3r_{ISCO})},\quad \epsilon_r=1-E_{ISCO}.$$

For equatorial aligned/counteraligned thin disk use corresponding Kerr specific angular momentum at ISCO, not Newtonian circular momentum. In units GM/c, along disk angular momentum direction with signed a, l=(r²-2a sqrt(r)+a²)/(r^(3/4) sqrt(r^(3/2)-3 sqrt(r)+2a)); handle near-extremal conditioning. Ideal thin disk dm_bh=E_ISCO dm_0, dJ=l GM/c dm_0 absent extra torques. Binding efficiency assumes zero stress at ISCO and neglects photon capture unless modeled separately. Verify Schwarzschild epsilon≈0.05719, prograde/retrograde limits and spin-up trajectories. Misaligned flow uses explicit warp model, not signed-spin formula blindly.

### Accretion state and spin-powered jets

P1 thick/ADAF, thin and slim disk regimes inspired by [B4], with state-dependent epsilon_r, epsilon_w, epsilon_j and H/R. Proposed mdot ratio boundaries near 0.01 and 1 are candidate defaults only; define exact normalization, transition smoothing/hysteresis and efficiency fit validity. A continuous state transition should not create energy or mass jumps. Use a versioned published model for photon trapping and magnetic flux; no arbitrary efficiency curve accepted as calibrated.

BZ-informed jet P≈(kappa/(4 pi c)) Phi_BH² Omega_H² f(Omega_H) in Gaussian convention, r_H=GM/c²(1+sqrt(1-a²)), Omega_H=a c/(2r_H). Choose flux convention and closure from reference (MAD/SANE assumptions explicit). E_rot=M c²[1-sqrt((1+sqrt(1-a²))/2)] is a theoretical upper extraction bound, not a second independent reservoir on top of M,J. For extraction dE/dt=P and dJ/dt=P/Omega_F along spin axis in a specified field-line angular velocity model; require correct zero-spin limit and 0<Omega_F/Omega_H<1 for extraction. Update gravitational mass and J simultaneously; check irreducible mass/area nondecrease for classical extraction and avoid negative reservoir. BZ power can exceed instantaneous accretion binding power, but only with actual spin depletion. Distinguish accretion-powered jets from BZ contribution; pass energy once to AGN ledger.

## 8. Mergers and binary evolution

P0: proximity relative to softening/kernel, relative velocity and pair binding using consistent softened potential and local tidal context. Require unresolved bound orbit and converging/physically eligible pair; proximity alone prohibited. Deterministic candidate ordering with disjoint pairs, distributed ID authority, owner migration and one commit. Merge dynamical mass/momentum, disk state and pending events with explicit energy/angular momentum ledger. Legacy Newtonian merger conserves tracked mass and momentum but cannot claim GW remnant accuracy. P1 remnant mass/spin and recoil use named numerical-relativity fits with domain checks, GW energy/momentum/angular momentum losses. Do not add orbital J arbitrarily to horizon J. Optional binary delay evolves separation/eccentricity via stellar/gas hardening and GW Peters rates when appropriate; report unresolved assumptions. Recoil affects orbit and occupation. Retain merger trees and event provenance.

## 9. Verification and interfaces

Suggested `src/physics/black_hole/` components `seed_provider`, `bh_neighborhood`, `bh_supply`, `bh_disk`, `bh_spin`, `bh_dynamics`, `bh_merger`, `bh_budget`, existing façade retained for compatibility. Public `.hpp` under matching `include/cosmosim/physics/black_hole/`; no competing namespace. `BhSupplyEstimate` is request+diagnostics; `BhAccretionTransaction` records accepted mass and momentum; `BhEmissionBudget` passes integrated radiative/mechanical budgets, source frame/time. Stable BH IDs and existing ID registry, migration and scheduler contracts mandatory.

| Gate | Measurement |
|---|---|
| Seed | unique across ranks/restart, gas mass/momentum/metal conservation, no duplicates in occupied halo |
| Bondi | analytic rates, unit conversion a^-3 density, spherical flow convergence and unresolved flags |
| Torque | reproduce reference equation and kinematic decomposition, phase exclusivity and resolution study |
| Disk | analytic constant/pulse supply and capped transition, exact integrated mass; stochastic/restart reproducibility |
| Dynamics | analytic drag limits, decay in controlled collisionless model, resolved-wake correction decreases |
| Spin | Kerr values, aligned/counteraligned and warped-disk ODE comparison, area/energy constraint |
| Merger | flyby rejection, bound pairs, distributed races, mass/spin/recoil fit benchmarks |
| Combined | nuclear inflow + disk + AGN self-regulation, dt and resolution convergence, comparison against published models |

Calibration parameters are versioned science inputs. Full published warp, remnant and super-Eddington coefficient sets remain model-selection/data dependencies; the present document defines interfaces and physical invariants, not an invented fit. No P1 certification until these are pinned.

## References

- B1. Davé et al. (2019), SIMBA: https://doi.org/10.1093/MNRAS/STZ937 (torque-limited + hot Bondi reference).
- B2. Fiacconi, Sijacki & Pringle (2018): https://arxiv.org/abs/1712.00023
- B3. Cenci et al. (2021): https://arxiv.org/abs/2011.06596
- B4. Huško et al. (2026), hybrid spin-dependent AGN: https://doi.org/10.1093/mnras/stag324
- B5. Super-Eddington spinning BH subgrid model (2026): https://doi.org/10.1093/mnras/stag003

These methods were developed for different resolutions and codes. The combination of supply and disk closures requires recalibration and cannot inherit their science validation automatically.
