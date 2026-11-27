# CHUI extended physics — equation-to-code numerical contracts, pass 01

**Campaign:** CHUI-EXT-PHYS-REF-v1 · source-only · **nothing here was compiled, executed or validated.**

Each section states: reference (verification status), equations, frame/gauge and units, closure,
discrete update, boundary/timestep treatment, invariants, domain of validity, and the independent
benchmark authored as test source. Reference status vocabulary:

* **derived** — follows from stated equations/defining constants in this document; no literature value used.
* **recalled lead** — a published source recalled from training knowledge. *External access to verify
  title/DOI/equation numbers was not used in this pass*; the lead must be verified in audit before
  qualification. No numerical coefficient in code depends on a recalled lead unless stated.
* **external data** — a coefficient set that must be supplied in a versioned, hashed table; the code
  refuses to run without it.

---

## S01 / N01 — homogeneous background authority (`core::LambdaCdmBackground`)

* **Equations (derived).** `E²(a) = (Ω_γ+Ω_ur)a⁻⁴ + Ω_cb a⁻³ + Σ_i Ω_ncdm,i(a) + Ω_k a⁻² + Ω_Λ`;
  `Ω_γ = a_rad T₀⁴/(c² ρ_crit,0)`, `a_rad = π²k⁴/(15ħ³c³)`; `Ω_ur = n_ur (7/8)(4/11)^{4/3} Ω_γ`;
  `Ω_ncdm,i(a) = (7/8) r⁴ Ω_γ a⁻⁴ I_ρ(y_i a)/I_ρ(0)`, `y_i,0 = m_i c²/(k r T₀)`; pressure with `I_P`.
  `I_ρ(y)=∫x²√(x²+y²)/(eˣ+1)dx`, `I_P(y)=∫x⁴/(3√(x²+y²))/(eˣ+1)dx`, `I_ρ(0)=7π⁴/120`.
  Each eigenstate has g = 2 (ν + ν̄), zero chemical potential, instantaneous decoupling.
* **Convention.** `cosmology.omega_matter` = Ω_b + Ω_cdm + Ω_ncdm,0 (CLASS / GADGET `Omega0` total
  clustering matter). `Ω_cb = Ω_m − Ω_ncdm,0` is what N-body particles must represent. Strict closure
  `Ω_γ+Ω_ur+Ω_m+Ω_k+Ω_Λ = 1` (1e-10); no renormalization. T_ncdm/T_cmb is an explicit input
  (`instantaneous_decoupling` = (4/11)^{1/3}, or a numeric value such as the CLASS default 0.71611,
  **recalled lead**: Lesgourgues & Tram 2011, CLASS IV).
* **Discretization.** Composite 16-node Gauss–Legendre (nodes by Newton iteration on P₁₆) on
  x ∈ [0, 80], 40 panels; 4001-point table in ln y ∈ [ln 1e-5, ln 1e6], 4-point Lagrange in (ln y, ln I);
  two-term asymptotic expansions outside (closed-form moments ∫xⁿ/(eˣ+1) = (1−2⁻ⁿ) n! ζ(n+1)).
  Distances: adaptive Simpson in ln a of `1/(aH)`; inversion by bisection in ln a.
* **Invariants.** E(1)=1; Ω_cb+Ω_ncdm,0 = Ω_m; P ≤ ρ/3; massless and nonrelativistic limits;
  legacy E(a) bit-identical when disabled; `cmbTemperatureK` throws when disabled (no second
  temperature authority).
* **Validity.** Homogeneous FLRW background only. a_begin must exceed matter–radiation equality
  (computed from the same background, radiation-like part = photons + ur + 3P_ncdm) because the
  Newtonian TreePM evolves no relativistic perturbations. No Einstein–Boltzmann evolution.
* **Single authority.** `cosmologyBackgroundConfigFromSimulation` is the only config→background path;
  the three base construction sites (`config.cpp` integration authority, `time_coordinator.cpp`
  cosmology dt, `reference_workflow.cpp` runtime) were switched to it. Inherited exception recorded in
  the ledger: `StellarEvolutionBookkeeper::evaluateStarAgeYears(a_f, a)` overload builds a default
  background (pre-existing; not used by the RT stellar provider of this pass).
* **Benchmarks (test source `unit_cosmological_background`, NOT_RUN).** EXT-N01-T01 FD limits and
  table-vs-quadrature 1e-9; EXT-S01-T01 Ω_γh² = 2.4728e-5 and n_γ = 410.7 cm⁻³ at 2.7255 K
  (**recalled lead** values, 2e-4/1e-3 tolerances); EXT-N01-T02 massless/NR limits from the number
  density formula (1e-5..1e-6), CLASS 93.14 eV relation (1e-3, recalled lead), no double counting
  (1e-14); EXT-S01-T02 legacy bit equivalence and refusals; EXT-S01-T03 EdS distances (1e-9);
  EXT-S01-T04 `.param.txt` admission/round trip/refusals.

---

## S02 — Boltzmann, atomic and spectral data (`cosmosim_spectral`)

### Provider status (`spectral/provider_status.hpp`)
Three states: `kAvailable`, `kPhysicalZero` (true vacuum/absent species), `kUnavailable` (unknown; carries a
reason; `require()` throws). Observable kernels use it wherever a physics input may be absent.

### Boltzmann transfer tables (`spectral/boltzmann_tables.hpp`)
* **Reference.** CLASS (Blas, Lesgourgues & Tram 2011) and CAMB (Lewis, Challinor & Lasenby 2000):
  **recalled leads** for file layouts. The CLASS `N:name` header parsing is generic; the CAMB 13-column
  layout (`k/h CDM baryon photon nu mass_nu total no_nu total_de Weyl v_CDM v_b v_b-v_c`) is a recalled
  lead and is verified against the header names when present.
* **Canonical epoch.** `delta_cb = (Ω_cdm δ_cdm + Ω_b δ_b)/Ω_cb` with Ω_cdm = Ω_cb − Ω_b from the S01
  authority; `delta_ncdm = Σ ρ_i(a) δ_i / Σ ρ_i(a)`; `theta_cb` mass-weighted, `theta_ncdm` weighted by
  (ρ_i + P_i) (relativistic momentum density); θ is w.r.t. conformal time in 1/Mpc (CLASS units).
  Only ratios and shapes are used (common arbitrary normalization).
* **Gauge.** Declared; chui_v1 header must match; CLASS `t_cdm` with declared synchronous gauge is
  refused; CAMB must be synchronous; velocity columns only in Newtonian gauge (linear-response
  initialization therefore requires Newtonian-gauge CLASS or chui_v1).
* **Interpolation.** Ratio linear in ln k on each epoch; linear in ln a between bracketing epochs; any
  k or a outside the common domain throws (no extrapolation). Eigenstate count/masses checked against the
  configured background (chui_v1 masses to 1e-6 relative; CLASS count only — masses are not recorded in
  CLASS files, a documented provenance limitation).
* **Tests (`unit_spectral_data`, NOT_RUN).** EXT-N02-T01 synthetic chui_v1/CLASS/CAMB parsing, exact ln-k
  interpolation for ratios linear in ln k, ln-a midpoint, out-of-domain throws, gauge/redshift refusals.

### Atomic lines and Voigt (`spectral/atomic_data.hpp`)
* Lyα: λ = 1215.67 Å, f = 0.4164, A = 6.265e8 s⁻¹, g = 2→6; 21 cm: ν = 1420.405751768 MHz,
  A₁₀ = 2.85e-15 s⁻¹ (f derived from A). **recalled leads (NIST ASD / standard reviews)**, flagged in the
  ledger for verification; A₁₀ literature spread ≈ 1 %.
* `σ(ν) = (e²/(4ε₀ m_e c)) f φ(ν)`, `Δν_D = (ν₀/c)√(2kT/m + b²)`, `a = A/(4πΔν_D)`,
  `φ = H(a,x)/(√π Δν_D)`.
* Faddeeva w(z), Im z ≥ 0: Weideman (1994) rational expansion with N = 32,
  `w = 2Σ_{n=1}^{N} a_n Zⁿ⁻¹/(L−iz)² + 1/(√π(L−iz))`, `Z = (L+iz)/(L−iz)`, `L = (N/√2)^{1/2}`,
  `a_n = (1/2M) Σ_{k=−M+1}^{M−1} e^{−t_k²}(L²+t_k²) cos(nθ_k)`, `θ_k = kπ/M`, `t_k = L tan(θ_k/2)`,
  M = 2N (**recalled lead**, implemented from the formula; coefficients computed, not transcribed).
* Independent reference: `H = (1/π)∫_{−π/2}^{π/2} exp(−(x + a tan φ)²) dφ` (exact substitution of the
  Voigt convolution, derived), adaptive 8-point Gauss–Legendre split at the peak.
* **Tests.** EXT-L01-T01: H(a,0) = e^{a²}erfc(a) (1e-12), Weideman vs quadrature for a ∈ [1e-4, 5],
  |x| ≤ 40 (1e-9 relative + 1e-14), Lorentz wing, ∫H dx = √π (1e-4). EXT-L01-T02/EXT-RA01-T01:
  A↔f consistency (5e-3), πe²/m_ec, T* = hν/k = 0.06817 K, λ₂₁ = 21.1061 cm, b(10⁴ K) = 12.85 km/s.

### Frequency axis (`spectral/frequency_axis.hpp`)
Explicit strictly increasing channel edges; per-channel Gauss–Legendre averaging with `subsamples` nodes
(bandwidth averaging is explicit). Test EXT-RA05-T01: cubic averaged exactly by 4 nodes.

---

## N02–N06 — massive-neutrino perturbation methods (`cosmosim_neutrinos`)

Models are distinct and separately admitted (`neutrinos.method`): `background_only` (S01 only),
`transfer_ratio`, `linear_response` (both integrated, unqualified), `particles` / `delta_f` (reference
kernels; admission FAIL_CLOSED). No model substitutes another.

### PM coupling contract (`gravity/pm_spectral_hook.hpp`, `pm_solver.cpp::applyDerivedPoissonOperator`)
`Ψ_k = −4πG W⁻² S_split M(|k|) δρ_cb,k / k²`, `M = 1 + (ρ̄_ν/ρ̄_cb) R(k)`, k ≠ 0 (k = 0 stays 0, so ρ̄_ν
enters only H(a)). Deposited density is cold matter only (Ω_cb particles). M commutes with window
deconvolution and the Gaussian long-range split; the short-range tree carries no neutrino force
(declared scale limit; `shortRangeOmissionBound(k) = M(k) − 1`). Observation of the cold amplitude is a
deterministic serial storage-order pass (`forEachLogicalSpectralModeSerial`), Hermitian multiplicity 1/2
for the r2c half spectrum, then a collective `allreduceSumDoublesInPlace`. The mesh |k| domain is checked
before any OpenMP loop (`requireWavenumberDomain`), so the noexcept multiplier never extrapolates. The
hook identity is folded into the PM collective fingerprint. CUDA and isolated/open PM refuse the hook.

### Transfer ratio (N03a)
`R(k,a) = δ_ncdm/δ_cb` from `BoltzmannTransferSet` (ln k / ln a interpolation, no extrapolation). Prescribed
linear-theory correlation with the current cold mode; does not respond to the simulated potential.

### Linear response (N03b) — derivation (derived)
Vlasov in comoving x, superconformal time `ds = dt/a²`, canonical velocity `u = dx/ds = q/m` (NR):
`∂_s f + u·∇f − a²∇φ·∂_u f = 0`, with `∇²φ = 4πG a² δρ` (ρ physical). Linearize, Fourier (`e^{ik·x}`):
`δf(s) = e^{−ik·uτ}δf_i + ∫ ds' e^{−ik·u(s−s')} a² i φ k·∂_u f̄`. Integrating over u with
`∫ e^{−ik·uτ} ∂_u f̄ = ikτ n̄ F(kτ)` gives
`δ_ν(s) = [δ_ν,i − θ_ν,i τ] F(k u₀τ) + ∫_{s_i}^{s} ds' (s−s') F(k u₀(s−s')) 4πG a⁴ (ρ̄_cb δ_cb + ρ̄_ν δ_ν)(s')`,
`F(y) = ∫x²f(x) sinc(xy)dx / ∫x²f`, `u₀ = kT_ν,0/(mc)`, `θ_s = i k·dx/ds = a θ_conformal`
(continuity `dδ/ds = −θ_s`). Cold limit F = 1 reproduces `d²δ/ds² = 4πGρ̄a⁴δ` (checked analytically and by
test EXT-N03-T04). History closure: phase-coherent amplitude approximation `δ_cb(k,s') = δ_cb(k,s)
A(k,s')/A(k,s)` with A the measured rms of the deposited cold density per |k| bin (**recalled lead**:
Ali-Haïmoud & Bird 2013, MNRAS 428, 3375). Quadrature: trapezoid over the PM-solve history, explicit
(weight (s − s') = 0 at s' = s). Bins without mesh modes interpolate ln A in ln k between observed bins.
F: 5001-point table on y ∈ [0, 50] (oscillatory Gauss–Legendre quadrature, ≥ 16 nodes per period) plus
the derived asymptotic series `∫x f sin(xy)dx = 1/(2y³)+1/(2y⁵)+3/(2y⁷)+…`.
* **Validity (admission).** Degenerate eigenstates; Newtonian-gauge tables with θ columns (CLASS
  Newtonian or chui_v1; CAMB refused); `<|p|>/(mc)` at a_begin below the user-declared
  `response_max_mean_velocity_c` (no default bound); periodic host PM; table epochs bracket z_begin.
* **State.** Bounded history (s, a, A_b, R_b) reserved through `MemoryGovernor`
  (`kCanonicalPersistent`, `capacity × (2 + 2 n_bins) × 8 B`); persisted as whole-block module sidecar
  `chui_neutrino_response_v1` (no restart schema bump; existing module-block persistence) after every
  refreshing PM solve and restored on first execute; replicated on all ranks after the reduction;
  deterministic replay truncates entries with a ≥ a_solve; a full history refuses.

### Particles and δf (N04, reference kernels only)
Lagrangian `L = −mca√((1+2Ψ)c² − (1−2Φ)v²)`, Φ = Ψ, first order: `dx/dt = q/(aε)`,
`dq/dt = −(ε/a)(1+v²/c²)∇Ψ`, `εc² = √(q²c²+m²c⁴a²)` (derived; dΦ/dt term neglected). δf weight
`w = 1 − f_ref(q)/f_i` with `f_i = f_ref(q_initial)` (**recalled lead**: Elbers et al. 2021, MNRAS 507,
2614); gravitational estimator mass `w × m_rest γ`; signed masses refused on the tree path; no clipping.
Thermal sampling by inverse CDF of `x²/(eˣ+1)`, counter-based SplitMix64 keyed by (seed, id, draw).
**Missing edges (FAIL_CLOSED):** neutrino particle species in `core::ParticleSpecies`, IC displacement and
bulk velocity from transfer tables (N02 IC adapter), signed PM deposition path, migration/restart lanes.

### Species statistics (N05)
CIC deposit against the BACKGROUND mean; `P_XY = V⟨Re δ_X δ_Y*⟩/N⁶`, optional `W⁻²` CIC deconvolution;
`P_m = f_cb²P_cc + 2f_cb f_ν P_cν + f_ν² P_νν`; shot noise `VΣm²/(Σm)²`; δf noise `VΣ(wm)²/M_bg²`;
response family `P_cν = R P_cc`, `P_νν = R² P_cc`. Reference serial radix-2 FFT (`spectral/reference_fft`).

### Tests (`unit_neutrino_methods`, NOT_RUN)
EXT-N03-T01 F(0)=1, small-y `1 − y²⟨x²⟩/6`, table vs quadrature 2e-7, asymptotic 1e-5; EXT-N03-T02 EdS
superconformal time 1e-12; EXT-N03-T03 transfer-ratio multiplier and domain refusal; EXT-N03-T04 cold
limit vs independent RK4 (|R−1| < 2e-3); EXT-N03-T05 sidecar round trip and replay; EXT-N03-T06 capacity
refusal; EXT-N04-T01 δf weight sign, signed-mass refusal, NR drift 1e-9, FD sample mean (6σ);
EXT-N05-T01 single-mode spectrum normalization, shot noise, ratio-family total. Missing (named for the
execution campaign): PM-integrated neutrino DMO run vs CLASS P(k) suppression, MPI/rank equivalence of
the reduction, restart equivalence through the real checkpoint writer.

---

## S03 — read-only leaf ray and lightcone service (`cosmosim_lightcone`)

* **Geometry (derived).** Leaves = populated patches (`cell_count > 0`); parents carry no cells, so covered
  parents are excluded by construction. Cell row = `first_cell + i + nx(j + ny k)` (canonical
  `amr::linearPatchCellIndex`). Point location: CSR uniform bucket index (≤ 256³ buckets, each leaf in
  every bucket it overlaps): O(1) expected, replacing linear patch scans. Half-open leaf boxes.
* **Traversal.** Comoving length s; at each crossing the next cell is found at `p + 10⁻⁶ Δx_min d`
  (ties on faces/edges/corners resolve toward the direction of travel); exit distance = nearest face of
  that cell measured in the probe image frame (periodic images handled by an explicit shift). Segments
  strictly increasing in s; Σ lengths = requested length. Periodic: unbounded unwrapped position with
  wrapped probes (box replicas). Open: stop at the face. A tiling hole throws; `max_segments` capacity
  throws (no silent truncation). Proper path length for a snapshot ray at fixed a: `dℓ = a ds`.
* **MPI.** Rank-local leaves; additive LOS integrals reduce per ray; ordered transfer needs
  `mergeRankSegments` (sorted by s then global patch id, overlap refused). The exchange of segment
  records between ranks is NOT implemented in pass 01 (single-rank ordered transfer only; ledger gap).
* **Lightcone.** `χ(a)` table on ln a (cumulative adaptive conformal-time integrals of the S01 authority,
  cubic Hermite interpolation, bisection inverse). Shell i owns χ between midpoints to its neighbours;
  replicas `⌈χ_max/L⌉`, missing modes k < 2π/L and repeated-structure flag recorded. Temporal weights linear
  in a between bracketing snapshots (no extrapolation). Observer frame `1+z_obs = (1+z)(1+v_los/c)` (first
  order, declared); snapshot Hubble velocity `v = aH dx`.
* **Not implemented (ledger):** HEALPix pixelization, distributed streaming/tiled lightcone accumulation
  with shell checksums during the run, gravitationally deflected rays (multi-plane lensing in C02 handles
  deflection on planes separately).
* **Tests (`unit_lightcone_geometry`, NOT_RUN).** EXT-S03-T01 axis ray segment lengths/rows; T02 oblique ray
  closure and midpoint containment, tie determinism; T03 AMR leaves via PatchSoa with metadata parent,
  open face, hole and capacity refusals; T04 rank merge and overlap refusal; T05 EdS χ(a) 1e-8, inverse,
  contiguous shells, replication count, temporal weights.

---

## F01 / F02 / R01–R06 — radiation foundations and reference backends

### F01 governed RT workspaces (`radiation_runtime`)
Live allocation formula (source estimate, not RSS): `radiationIntervalWorkspaceBytes = 8·(6·4+1)·n_values
+ coupled·3·n_cells·(4 + 6·8)` bytes, i.e. ≈ 6.25× the canonical 4-lane radiation state plus three gas
work copies (the prior ledger's "≈5×" undercounted the source lane and the stage copies). Reserved as
`kPhaseResident` through the run's `MemoryGovernor` (forwarded from `RuntimeServices` by
`SourceRuntimeImpl` → `RtBhAgnSourceBridge::setMemoryGovernor` → `RadiationRuntime`) BEFORE the work
copies are allocated and committed after; released at interval end (RAII). **Not done (ledger):** level
solver topology/ghost buffers are not in the formula; RT topology construction remains O(P²); MPI
owned-face halo exchange, regrid transfer and synchronized reflux across ranks remain refused (A14).

### F02 stellar luminosity provider (`radiation_stellar_sources`)
`E_g = M_* ∫_{τ0}^{τ1} L_g(τ, Z) dτ`, L_g from the external table `chui_ssp_group_luminosity_v1`
(group luminosity per unit BIRTH mass vs age and metallicity; source/IMF declared; edges must equal the RT
group edges to 1e-9); piecewise-linear in age (the trapezoid over the union of nodes is the exact
integral of the interpolant), linear in log₁₀Z; no extrapolation in age or Z (the age grid must start at 0).
Ages from `formation_scale_factor` via the S01 background → admitted only in cosmological modes.
Deposition into the containing radiation record cell (bucket-indexed locator; no O(P) scan per star),
isotropic, comoving a⁴ storage like the AGN sources; emission is a pure function of age, so disjoint
accepted intervals emit exactly once (no pending state, nothing to checkpoint). Photon number is not
inferred: the report keeps the table Q integral and the transported photon equivalent E_g/⟨hν⟩_g of the RT
reference shape (auditable mismatch; group cross sections stay those of the declared reference shape —
a named multigroup approximation). Admission now requires the table (refusal replaced by explicit
requirements); RT's existing single-rank/no-regrid/no-MHD restrictions still apply.

### R01 / R02 (existing, unchanged)
Multigroup M1 + PD-ARS IMEX, explicit comparator, HLL/GLF, realizability and the ionizing3 threshold
edges (13.598434 / 24.587389 / 54.417765 eV with explicit terminal edge) are EXISTING_UNQUALIFIED; no
change in this pass. Cosmological reduced-c remains refused. Photoionization reconciliation with the
thermochemistry operator remains refused (FAIL_CLOSED) — not closed in this pass.

### R03 reference backends on a uniform patch (`radiation_reference_transport`)
* **Discrete ordinates (derived).** Step (first-order upwind) discretization
  `I_c = (η + χ_s J + Σ_a |n_a| I_up,a/Δx_a)/(χ_t + Σ_a |n_a|/Δx_a)` per Fibonacci direction (equal
  weights 4π/N), upwind sweep order, source iteration for isotropic scattering, periodic inflow iteration,
  vacuum inflow at open faces; exact for uniform periodic media (I = η/χ). Reflecting boundaries refused.
* **VET.** `f_ij = Σ w I n_i n_j / Σ w I` (convex combination of projectors: PSD, trace 1); moment system
  `∂E/∂t + ∇·F = cχ_a(aT⁴ − E)`, `∂F/∂t + c²∇·(fE) = −cχ_t F`; LLF fluxes at speed c (eigenvalues of f
  along an axis in [0, 1] ⇒ stable for `Σ_a cΔt/Δx_a ≤ cfl`), backward-Euler local sources, closure from a
  steady formal solution of the current emissivity every `closure_cadence` steps (time-independent VET
  closure, declared), realizability projection |F| ≤ cE counted. Distinct backend; never routed through
  the M1 Riemann solver. Leads: VETTAM (Menon et al. 2022), Davis, Stone & Jiang 2012 (**recalled leads**).
* **FLD comparator.** `D = cλ(R)/χ_t`, Levermore–Pomraning `λ = (2+R)/(6+3R+R²)` (**recalled lead**:
  Levermore & Pomraning 1981; limits 1/3 and 1/R verified in test), backward Euler, harmonic face D,
  Jacobi-preconditioned CG (refuses non-convergence), periodic/reflecting/vacuum boundaries.

### R04 causal rays / R05 HG Monte Carlo (`radiation_causal_transport`)
* Causal packets at `x_e + c(t − t_e)n`, exact exponential attenuation per leaf segment, escape at open
  faces, explicit `unresolved` ledger below an energy floor, exactly-once closure
  `emitted = absorbed + escaped + unresolved + in_flight`, bounded queue refusing over-capacity emission
  before mutation, serialize/restore of packet state. Steady-state `traceDirectRays` keeps its
  infinite-speed label. Direct/diffuse hand-off into M1 is NOT implemented (no hybrid).
* Analog MC with Henyey–Greenstein sampling `μ = (1+g² − ((1−g²)/(1−g+2gξ))²)/(2g)` (**recalled lead**,
  property ⟨μ⟩ = g verified in test), counter-based RNG (`radiationCounterRandom`), batch standard errors,
  truncation ledger. Asymptotic-preserving/implicit MC (R05 P1) FAIL_CLOSED: no published discrete
  contract selected.

### Tests (`unit_rhd_reference_backends`, NOT_RUN)
EXT-R03-T01 FLD limits, Gaussian thick diffusion amplitude 1/√2 (±0.02), conservation 1e-10, LTE
relaxation; T02 DO exact uniform solution, isotropic f (±0.02), scattering J = η/χ_a; T03 VET budget
closure 1e-10, relaxation `1 − e^{−cχ_a t}` (±0.02); EXT-R04-T01 light travel, attenuation, ledger
closure 1e-14, restart, queue refusal; EXT-R05-T01 ⟨μ⟩ = g, cube escape vs exact angular average (5σ),
energy closure; EXT-F02-T01 table parse, exact age integral 1e-12, log-Z interpolation, domain refusals,
F01 formula. Missing (named): Strömgren sphere with stellar sources through RadiationRuntime, shadow and
crossing-beam M1 vs VET comparison, MPI/regrid RT, RT–MHD coupling.

---

## CR01–CR08 — cosmic rays (`cosmosim_cosmic_rays`, kernels on a uniform patch)

* **CR01 grey (recalled lead: Pfrommer et al. 2017 energy equation).**
  `∂e/∂t + ∇·(eu) = −P∇·u − ∇·F_diff − ∇·F_st + v_st·∇P + Q − Λ`, `P = (γ_cr − 1)e` (declared γ_cr; 4/3 only
  for the ultrarelativistic grey approximation). Gas: momentum `−∇P`, energy `−u·∇P` and `H_st = −v_st·∇P ≥ 0`.
  Donor-cell advection with face velocity `(u_L+u_R)/2`; central `∇·u`, `∇P` (the CR/gas pair is conservative
  only to truncation error for non-uniform u; the ledger reports the residual). No CR allocation unless a
  kernel object is created (lazy).
* **CR02 anisotropic diffusion (recalled lead: Sharma & Hammett 2007 limited scheme).**
  `F = −κ_∥ b(b·∇e) − κ_⊥(∇e − b(b·∇e))`; face normal gradient from the two cells, transverse gradients by a
  4-argument monotonized-central limiter of the one-sided differences of both cells; face b = normalized
  average; weak field `|B| < threshold` → declared isotropic κ (never an undefined direction). Explicit bound
  `Δt ≤ cfl/(2κ_max Σ Δx⁻²)`; RKL2 super-steps reuse `mhd::makeRkl2Coefficients` (Meyer, Balsara & Aslam
  2014, as already cited in the MHD module). Coefficients are inputs (no universal constants).
* **Streaming (recalled lead: Sharma, Colella & Martin 2010 regularization).**
  `v_st = −v_s b tanh(l_reg b·∇e / e)`, `F_st = (e + P)v_st`; heating `H = v_s tanh(·)(γ−1)b·∇e ≥ 0` removed
  from the CR and added to the gas in the same step (exact balance; no double counting); explicit bound
  `Δt ≤ cfl Δx/v_s · min(1, Δx/l_reg)`.
* **CR03 two-moment (derived).** `∂F/∂t = −(F − F_target)/τ`, `F_target = (e+P)(u+v_st) − K∇e`,
  `τ = 3K/c_r²` — algebraically identical to `∂F/∂t + (c_r²/3)∇e = −(c_r²/3)K⁻¹(F − F_eq)`. Staggered
  face-normal fluxes, exact exponential relaxation per face with the normal diffusivity
  `κ_n = b_a²κ_∥ + (1−b_a²)κ_⊥` (approximation for oblique fields: the full tensor coupling of the face
  normal and transverse flux components is reduced to κ_n; recorded), donor-cell enthalpy flux, CR02 limited
  diffusive flux; gas impulse = time integral of `f_gas = (1/3)K⁻¹(F − F_eq)` using the exact exponential
  average of F, energy `(u+v_st)·impulse` split between adjacent cells. Stiff limit = CR02 stencil
  (spatially asymptotic preserving); `Δt ≤ cfl/Σ c_r/Δx`. Periodic/reflecting only; positivity refusal.
  Not a copy of photon M1 (isotropic CR pressure, no M1 closure). Leads for later comparison: Jiang & Oh
  2018, Thomas & Pfrommer 2019 (**recalled leads**, formulations not transcribed).
* **CR04 spectra (recalled leads: Girichidis et al. 2020; Miniati 2001).** p̃ = p/(mc) bins; fixed-bin
  reference (number only, f constant per bin) and piecewise power law (number + energy conserved,
  slope by bisection on the monotone ⟨T⟩(q), bound slopes counted). Momentum-space fluxes
  `Φ_N = 4πp̃²f_up A`, `Φ_E = TΦ_N`, in-bin work `W = 4π∫p̃²f T' A dp̃` (24-node GL in ln p̃), outer-edge
  losses ledgered, zero inflow from outside, explicit with positivity refusal; pressure from the spectrum.
  Spectral two-moment (number/energy fluxes per bin): **NOT_STARTED** (no verified discrete contract).
* **CR05 sources.** Injection normalizes the declared p^{−α} shape to the caller's (already debited)
  energy exactly (fixed bins: normalized to the represented energy). Electron synchrotron + IC (Thomson,
  derived): `dp̃/dt = −C p̃ γ`, `C = (4/3)σ_T(U_B + U_rad)/(m_e c)`, `U_rad` = u_cmb(a) of S01 (caller passes
  it). Proton Coulomb/hadronic coefficients and thermalized hadronic fraction only from a versioned
  `chui_cr_loss_rates_v1` table with a mandatory reference line; missing → unavailable → refusal. Secondary
  e±, γ-ray emissivity, CR ionization: providers report UNAVAILABLE.
* **CR06/CR07.** Eulerian electron spectra = CR04 + CR05 radiative losses (tracer spectra NOT_STARTED).
  Memory: `CrSpectralState::bytesFor = lanes·bins·cells·8` reserved as `kCanonicalPersistent` before
  allocation (lanes 2 for piecewise power law, 1 for fixed bins); e.g. 2 species × 16 bins × 2 lanes at
  256³ = 8 GiB (matches master §12); the full 8-scalar two-moment variant would be 32 GiB at 256³ (not
  implemented). Groups are never dropped under pressure (reservation fails instead).
* **Runtime (I01/I02).** `cosmic_rays.model ≠ none` is refused at admission with the exact edge: no MHD
  level-stepper CR stage, no AMR/halo/regrid/migration/restart lanes, no feedback budget debit.
* **Tests (`unit_cosmic_ray_transport`, NOT_RUN).** EXT-CR01-T01 advection conservation/ledger; CR02-T01
  field-aligned rows and variance `σ₀² + 2κt`; CR02-T02 oblique-field monotonicity (property claimed by the
  recalled lead; to verify); CR02-T03 RKL2 vs explicit; CR02-T04 streaming balance; CR03-T01 diffusive
  limit vs one-moment (3 %); CR04-T01 slope recovery; CR04-T02 momentum-space adiabatic work
  `E₀e^{At}`, number conservation; CR06-T01 synchrotron cooling to p₀/2 after one cooling time (10 %);
  CR05-T01 loss refusal and gas/escape split.

---

## L01–L04, C01–C04, RA01–RA06 — observables (`cosmosim_observables`, read-only)

* **Science accessor.** One authority per quantity from `SimulationState`: n_X from composition species lanes
  (persistent `ChemicalSpecies` slot order) over the proper leaf-cell volume; n_e by charge neutrality of
  H/He ions (metal electrons excluded, stated); T from the kelvin cache; v peculiar; B only when a CT
  sample is attached. Missing inputs → `kUnavailable` → `require()` throws (no zero substitution).
* **L01 forest.** Velocity axis = proper Hubble velocity `u = aHs + v·n` (periodic `aHL`) for snapshot
  sightlines; segment-averaged Voigt/Doppler profile (3-point GL per b/4) preserves ∫τ du = Nσλ₀ for any
  Hubble spread; `τ = Σ n ℓ σ_int λ₀ ⟨φ_v⟩`, `φ_v = H(a,Δv/b)/(√πb)`, `a = Aλ₀/(4πb)` (derived from S02).
  Statistics: mean flux, flux PDF, 1-D `P_F(k) = L_v|δ̃_k|²/N²` (k in s/m), observed-frame EW, column.
  Doppler-only = named comparator. **NOT_STARTED:** FGPA comparator (needs a recombination-rate provider),
  3-D flux statistics, lightcone forest driver (the kernel supports log-λ velocity coordinates; no driver).
* **L02/L03 resonant MC (reference, unaccelerated).** Exact parallel-velocity sampler (derived two-component
  envelope: Gaussian on u < x − d accepted with d²/((x−u)²+a²); truncated Cauchy on u ≥ x − d accepted with
  e^{−u²}/max e^{−u²}; d = max(1,|x|/2)); coherent atom-frame redistribution `x' = x − u·k + u·k'`, optional
  recoil `g = hν₀/(c√(2kTm_H))` (derived, equal to Field's hΔν_D/2kT), isotropic or dipole phase functions
  only (no unverified mixtures), dust absorption, optional uniform Hubble flow, separate estimators (escape
  spectrum, escape/dust/truncation fractions, scattering rate, momentum deposition). Core skipping,
  diffusion hybrids, dynamical/WF coupling: **FAIL_CLOSED/NOT_STARTED**.
* **RA01 HI 21 cm.** `α = 3c³A₁₀h n_HI φ_v/(32πν₀²kT_s)` (Einstein relations, derived); ordered formal solution
  from the far end; background T_cmb(z) (S01) or 0; observer-frame δT_b = (T_b − T_γ)/(1+z); thin and
  saturated-spin shortcuts separately named. Spin temperature (Field 1958, **recalled lead**), x_c from a
  versioned κ₁₀ table only, x_α = 4P_αT_*/(27A₁₀T_γ) (**recalled lead**) available as a function but the
  pipeline passes x_α = 0 (no admitted P_α provider). **RA02 cosmic dawn / semi-numerical:** NOT_STARTED.
* **C01.** `chui_cmb_cl_v1` tables (units μK², dl or cl normalization, provenance); flat-sky Gaussian
  realization `⟨|a_k|²⟩ = A C_ℓ`, Hermitian symmetry, counter-based draws. Primary TE/EE realizations and
  lensing-potential inputs: parsed, not realized (TT only).
* **C02.** Plane Poisson `∇²_⊥ψ = (3H₀²Ω_m/c²)Σ`, deflection by spectral gradient, multi-plane recursion
  with lens-lens coupling, bilinear periodic interpolation, κ from the finite-difference Jacobian; Born
  comparator. Single plane: multi-plane ≡ Born analytically (test). Lightcone plane construction from
  the simulation (projected δ/a per shell) has NO driver yet (kernel only).
* **C03.** y, τ_e, kSZ (receding ⇒ decrement), nonrelativistic tSZ factor `x coth(x/2) − 4` (Kompaneets,
  derived); ISW/RS `ΔT/T = (2/c²)∫Φ̇ dt`, Φ̇ from two bracketing potential snapshots linear in a.
  Relativistic tSZ: **FAIL_CLOSED** (no verified kernel).
* **RA03–RA05 radio.** Free-free RL 5.14 with Draine 2011 Gaunt (**recalled leads**), Kirchhoff absorption;
  synchrotron single-particle kernel F, G by the derived single-integral representations
  (`F(x) = x∫cosh(5t/3)e^{−x cosh t}/cosh t dt`, `G(x) = x∫e^{−x cosh t}cosh(2t/3)dt`), emission/absorption by
  integration over an arbitrary electron spectrum (time-dependent CR04 spectra can feed it), RL power-law
  closed forms kept as an independent comparator (test asserts agreement); Faraday rotation
  `ρ_V = e³λ²n_eB_∥/(4π²ε₀m_e²c³)`; exact constant-coefficient Stokes step by the 5×5 augmented matrix
  exponential (Padé(6,6), scaling-squaring; valid for singular K); DM/RM with (1+z)⁻¹, (1+z)⁻² weights.
  Named omissions: circular emission, Faraday conversion, polarized absorption, DELO, recombination and
  molecular lines, FRB populations.
* **I05 pipeline.** `[observables]` → `ObservablePipeline` in `AnalysisRuntime::executeDiagnostics`
  (completed-step label, cadence) → leaf geometry + accessor → kernels → FITS (`.partial` + rename) and a
  `chui_observables_summary_v1` text summary. Single rank, Part 2 thermochemistry required.
* **Tests (`unit_observables_kernels`, NOT_RUN).** EXT-L01-T03 (Doppler core 2e-3, thin EW 1 %, damping
  wing 1 %, velocity-integral conservation 1e-3); EXT-RA01-T02 (1.823e18 column relation 1 %, two-layer
  ordered solution 1e-12, spin limits); EXT-C03-T01 (slab y 1e-14, kSZ sign, tSZ null 217.5 ± 1.5 GHz, ISW);
  EXT-C02-T01 (Born vs multi-plane 2 %); EXT-C01-T01 (white-noise variance 5 %); EXT-RA03/RA04-T01
  (Kirchhoff, RM 0.812, Faraday screen 1e-12, singular K); EXT-RA03-T02 (∫F = 8π/9√3 2e-3, kernel vs closed
  form 2–3 %, Π 1 %); EXT-L02-T01 (sampler moments 5σ, photon budget); EXT-O01-T01 FITS blocking.
