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
