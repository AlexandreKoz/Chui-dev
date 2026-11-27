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
