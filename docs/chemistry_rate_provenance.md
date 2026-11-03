# CHUI thermochemistry rate provenance (Part 2.5, 2026-11-03)

This file records, for every rate and cooling term used by the
`equilibrium_primordial`, `fixed_hybrid` and `fixed_detailed` backends:

- the publication and the equation or table it comes from;
- its coefficients, units and validity range;
- what was actually checked;
- independent reference values.

The reference values are computed **outside** the C++ code from the
published equations and are asserted by
`tests/unit/test_chemistry_rate_references.cpp`.

**Every test named here is NOT RUN (USER-RESERVED).** No compilation, test,
benchmark or simulation was executed in this pass. Production qualification
is **not** declared.

## 0. Status summary

| Status (`RateVerification`) | Meaning | Entries |
|---|---|---|
| `kVerifiedAgainstSource` | Coefficients, exponents, temperature units and piecewise ranges compared line by line with the **primary** text (arXiv copy) | R01–R12, R14, R16, R19; R07–R09 (data contract + Rahmati); cooling C01, C03, C11 |
| `kDerivedFromConstants` | Computed from CODATA constants | Compton (C02) |
| `kVerifiedAgainstSecondarySource` | Primary not accessible. Checked against a peer-reviewed compilation that reproduces it | R13, R15, R17 (Abel et al. 1997); R20, R21 (Glover & Clark 2012 restating NL97); H₂ LTE cooling C04, C05 (KROME = Grassi et al. 2014, eq. 10); fine structure C10 (KROME Table B) |
| `kUnverifiedTranscription` | Neither the primary nor a reproducing compilation was checked | R18 (Jura 1975 dust H₂ formation, 3e-17); C06 (0.4 eV per H₂ photodissociation, Black & Dalgarno 1977) |

The gate (`physics.thermochemistry_allow_unqualified_rates`) **remains in
force for all three backends**:

- `fixed_hybrid` and `fixed_detailed`: `NetworkManifest::qualified()` is
  false. Secondary or unverified entries never qualify
  (`rateVerificationQualifies`).
- `equilibrium_primordial`: its rate set (R01–R09, C01–C03) is
  primary-verified. Provenance alone never lifts the gate, though, and the
  numerical validation (USER-RESERVED) is not yet recorded.

Lifting the gate is an owner decision after the user-reserved tests pass.

## 1. Corrections made in Part 2.5

| Item | Before | After | Evidence |
|---|---|---|---|
| R16 H₂ + H → 3H | removed (corrective gate) | **Restored** from Glover & Abel 2008 (GA08) Table A1 reaction 9. v=0 limit k = 6.67e-12 T^½ exp[−(1 + 63593/T)] (Mac Low & Shull 1986). LTE limit k = 3.52e-9 exp(−43900/T) (Lepp & Shull 1983). Log blend, GA08 eq. 13: log k = [(n/n_cr)/(1+n/n_cr)] log k_LTE + [1/(1+n/n_cr)] log k_v=0. Critical density, eqs. 14–17: 1/n_cr = x_H/n_cr,H + x_H₂/n_cr,H₂ + x_He/n_cr,He, with x_H = n_HI/n, x_H₂ = 2n_H₂/n, x_He = n_He/n, where n is the number density of hydrogen nuclei. n_cr,H = dex[3.0 − 0.416 log T₄ − 0.327 (log T₄)²]; n_cr,H₂ = dex[4.845 − 1.3 log T₄ + 1.62 (log T₄)²]; n_cr,He = dex[5.0792 (1 − 1.23e-5 (T − 2000))]. Evaluated in log space (no underflow to log 0). The 4.48 eV binding energy is debited to `chemical_cooling`, as for R15. | The Abel et al. 1997 k13 text prints k13 = 1.067e-10 T^2.012 exp[−(4.463/T)(1 + 0.2472 T)^3.512] with T in eV, exactly the removed implementation. It decreases above ~1 eV, so the printed form is unusable and probably contains a typo. Dove & Mandy 1986 was not accessible. GA08 is the replacement. |
| H₂ self-shielding Doppler parameter | b5 = max(0.29 √(T/1e3), 0.01), i.e. 0.0917 at 100 K | b5 = √(kT/m_H)/1 km s⁻¹ = 0.0908 √T, i.e. 0.908 at 100 K (`h2ThermalDopplerB5`) | Thermal b = √(2kT/m_H₂). The old value was **10× too small**, over-shielding by up to two orders of magnitude in the (1 + x/b5)⁻² regime. Turbulent broadening is not added and is named in `omitted_channels`. DB96 eq. 37 itself was verified (primary). |
| Photoelectric heating | ψ = G₀ T^½ / n_e | ψ = G₀ T^½ / (n_e φ_PAH), φ_PAH = 0.5 | Wolfire et al. 2003 eqs. 23–24, standard model φ_PAH = 0.5 |
| He 2³S term 5.01e-27 T^-0.1687 f₅ e^{−55338/T} n_e² n_He⁺ | in `cooling::primordial` | **Removed**. Named in `omitted_channels` | Not in KWH96 Table 1 (verified). KROME lists Cen (1992) He triplet excitation, but Cen 1992 was not accessible. It contributes 8e-5 of the HeII cooling at 3e4 K and n_e = 1 cm⁻³. |
| O I 63 µm ΔE | 228 K | 227.7 K | hc/(kλ) for λ = 63.18 µm. KROME rounds to 230 K. |
| R12 3.2e4 K cap | labelled an "implementation choice" | GP98 Table 1 H8 states validity **1 ≤ T ≤ 32000 K**. The value at the limit is held beyond it, which is documented and is not an extrapolation of the cubic. | GP98 primary |
| Metal-table abundance scaling | X_e / X_e,ref | (X_e/X_H)/(X_e,ref/X_H,ref), new required table field `reference_hydrogen_mass_fraction` | WSS09 eq. 4 defines 10^[i/H] relative to hydrogen. The old form was 7.6% wrong for X_H = 0.76 against the CLOUDY X_H = 0.7065. |
| Metal-table UVB | undeclared | Required header field `uv_background`. Compatibility policy in §5. | WSS09 tables embed HM01 + CMB photoionization |

Network manifest versions bumped: `chui_hybrid_light` v2 and
`chui_detailed_hybrid_co_nl97` v2.

## 2. Per-rate record

The temperature T is in K unless stated. "Primary" means the arXiv text of
the cited paper was read in this pass.

### H/He ionization and recombination: KWH96 Table 2 (primary)

KWH96 takes these from Black (1981) and Cen (1992). Its assumptions are an
optically thin ionization balance and collisional ionization equilibrium
plus photoionization (KWH96 §3). f₅ = [1 + (T/1e5)^½]⁻¹ and
T₃ = T/1e3, T₆ = T/1e6.

| ID | Rate | Coefficients as printed (matched exactly) |
|---|---|---|
| R01 | Γ_eH0 | 5.85e-11 T^½ e^{−157809.1/T} f₅ |
| R02 | α_H+ | 8.4e-11 T^−½ T₃^−0.2 (1 + T₆^0.7)⁻¹ |
| R03 | Γ_eHe0 | 2.38e-11 T^½ e^{−285335.4/T} f₅ |
| R04 | α_He+ + α_d | 1.5e-10 T^−0.6353 + 1.9e-3 T^−1.5 e^{−470000/T}(1 + 0.3 e^{−94000/T}) |
| R05 | Γ_eHe+ | 5.68e-12 T^½ e^{−631515/T} f₅ |
| R06 | α_He++ | 3.36e-10 T^−½ T₃^−0.2 (1 + T₆^0.7)⁻¹ |

**Case A vs case B.** KWH96 does not state the case. Its α_H+ gives
5.10e-13 cm³ s⁻¹ at 1e4 K. The Ferland et al. (1992) fits printed in GA08
Table A1 give:

- case A: 1.269e-13 (315614/T)^1.503 [1 + (604625/T)^0.470]^−1.923, which is
  4.30e-13 at 1e4 K;
- case B: 2.753e-14 (315614/T)^1.5 [1 + (115188/T)^0.407]^−2.242, which is
  2.59e-13 at 1e4 K.

The KWH96 fit is therefore **case-A-like** (+19% versus case A). CHUI uses it
everywhere, which is correct for optically thin gas. **Case B for
self-shielded gas is not implemented** and is a named regime limitation.
Both Ferland fits are exposed only as cross-check functions.

**Electron bookkeeping.** n_e = Σ_i q_i n_i over the stored species, by
charge neutrality, clamped ≥ 0. Every reaction's net change balances H, He,
C, O nuclei and charge in integer arithmetic (`validateNetworkManifest`).
Collisional ionization debits χ × rate to `chemical_cooling`. The binding
ledger is state-based (`bindingEnergyDensity`). Metal electrons enter the
detailed network only through C⁺; other donors are omitted and named.

### Photo-rates R07–R09 (data + Rahmati 2013, primary)

Γ and the photoheating per ion come from the user TREECOOL table. The
column order and units were verified against the FG20 README; see
`tools/chui_data/data_sources.md`.

The local shielding factor is Rahmati et al. 2013 eqs. 13–14:
f = 0.98 [1 + (n_H/n₀)^1.64]^−2.28 + 0.02 (1 + n_H/n₀)^−0.84, with
n₀ = 6.73e-3 (σ̄/2.49e-18)^−2/3 T₄^0.17 Γ₋₁₂^2/3 (f_g/0.17)^−1/3 cm⁻³.
Two assumptions apply here:

- **σ̄ = 2.49e-18 cm² is the HM01 z = 3 value.** Rahmati Table 2 shows that
  σ̄ depends on the UVB model and redshift, so with FG20 it is an
  approximation. f_g = 0.17 is fixed.
- The same factor is applied to HeI and HeII rates and heating (an
  assumption; Rahmati fits HI only).

### Molecular hydrogen chemistry

| ID | Reaction | Form | Source and status | Validity |
|---|---|---|---|---|
| R10 | H + e → H⁻ + γ | 1.4e-18 T^0.928 e^{−T/16200} | GP98 Table 1 H3, primary | — |
| R11 | H⁻ + H → H₂ + e | 1.5e-9 (T ≤ 300), 4.0e-9 T^−0.17 | GP98 H5, primary | — |
| R12 | H + H⁺ → H₂⁺ + γ | dex[−19.38 − 1.523 lg T + 1.118 lg²T − 0.1269 lg³T] | GP98 H8, primary | **1–32000 K**, held beyond |
| R13 | H₂⁺ + H → H₂ + H⁺ | 6.4e-10 | Karpas et al. 1979 via Abel et al. 1997 k10 ((6.4 ± 1.2)e-10), secondary | constant |
| R14 | H⁻ + H⁺ → 2H | 5.7e-6 T^−½ + 6.3e-8 − 9.2e-11 T^½ + 4.4e-13 T | GP98 H7 ("fit by PAMS"), primary | No range printed. Grows above ~1e5 K, where H⁻ is negligible. Kept as published. |
| R15 | H₂ + e → 2H + e | 5.6e-11 T^½ e^{−102124/T} | Donahue & Shull 1991 via Abel 1997 k12, secondary | — |
| R16 | H₂ + H → 3H | GA08 reaction 9 blend (§1) | GA08 Table A1 + eqs. 13–17, primary | No range printed in GA08. H₂–H₂ and H₂–He dissociation (GA08 reactions 10–11) are **omitted and named**. |
| R17 | H₂⁺ + e → 2H | 1.0e-8 (T < 617), 1.32e-6 T^−0.76 | Schneider et al. 1994 via Abel 1997 k18, secondary | — |
| R18 | 2H(grain) → H₂ | 3e-17 (D/D_MW) n_H n_HI | Jura 1975, **unverified**. GML07 quotes only "R_H₂ ~ 1e-17 s⁻¹ (Jura 1974)". | No grain-temperature or sticking dependence |
| R19 | H₂ + γ → 2H | 3.3e-11 G₀ f_sh e^{−τ}, τ = 2e-21 N_H = 3.74 A_V | GML07 eqs. 32 and 34; DB96 eq. 37 (f_sh = 0.965/(1+x/b5)² + 0.035/√(1+x) e^{−8.5e-4 √(1+x)}, x = N_H₂/5e14), primary | Local Jeans-length column (not radiative transfer), thermal b only |

### Carbon/oxygen (detailed backend)

R20 and R21 follow the NL97 reduced network:

- k₀ = 5e-16 cm³ s⁻¹ and k₁ = 5e-10 cm³ s⁻¹;
- Γ_CHx = 5e-10 G₀ e^{−2.5 A_V} s⁻¹ and Γ_CO = 1e-10 G₀ e^{−2.5 A_V} s⁻¹;
- β = k₁x_O/(k₁x_O + Γ_CHx/n).

These were checked against Glover & Clark 2012 (MNRAS 421, 9) eqs. 2–5,
which also note typos in NL97 eqs. 18/20 (secondary; the NL97 text was not
accessible). NL97 assumes fully molecular hydrogen and treats all non-CO
carbon as C⁺. CO line self-shielding is omitted.

### Cooling and heating channels

| ID | Channel | Source and status |
|---|---|---|
| C01 | KWH96 Table 1 terms, matched line by line. Collisional excitation: HI 7.50e-19 f₅ e^{−118348/T}; HeII 5.54e-17 T^−0.397 f₅ e^{−473638/T}. Collisional ionization: 13.6/24.6/54.4 eV × Γ, debited to `chemical_cooling`. Recombination: 8.70e-27 T^½ T₃^−0.2 (1+T₆^0.7)⁻¹ (H⁺); 1.55e-26 T^0.3647 (He⁺); 3.48e-26 (He⁺⁺, same form). Dielectronic: 1.24e-13 T^−1.5 e^{−470000/T}(1+0.3e^{−94000/T}). Free-free: 1.42e-27 g_ff T^½ (n_H+ + n_He+ + 4 n_He++) n_e, g_ff = 1.1 + 0.34 exp[−(5.5 − lg T)²/3] (KWH96 eq. 23). | KWH96, primary |
| C02 | Compton: 4σ_T a k T_γ⁴/(m_e c) (T − T_γ) n_e, T_γ = 2.7255(1+z) | CODATA. Equals the KWH96 5.41e-36 n_e T (1+z)⁴ for T₀ = 2.7 K. Applied in cosmological modes only. |
| C03 | UVB photoheating Σ q_i n_i × shielding | data + Rahmati 2013 |
| C04 | H₂ rovibrational. Λ = Λ_LTE/(1 + n_cr/n_H) with n_cr/n_H = Λ_LTE/Λ(n→0). Λ_LTE (HM79): rot 9.5e-22 T₃^3.76/(1 + 0.12 T₃^2.1) e^{−(0.13/T₃)³} + 3e-24 e^{−0.51/T₃}; vib 6.7e-19 e^{−5.86/T₃} + 1.6e-18 e^{−11.7/T₃}. lg Λ(n→0) = −103.0 + 97.59 lg T − 48.05 lg²T + 10.80 lg³T − 0.9032 lg⁴T, valid 10–1e4 K (clamped). | GP98 eqs. A5–A7 primary. HM79 term secondary (KROME eq. 10; the PDF prints "1.6×10^18", an evident sign typo for 1.6e-18). |
| C05 | H₂ formation heating 4.48 eV × n/(n + n_cr), plus dissociation debits (R15, R16) | constants + C04, secondary |
| C06 | 0.4 eV per H₂ photodissociation | Black & Dalgarno 1977, **unverified** |
| C10 | C⁺ 158 µm: A = 2.4e-6 s⁻¹, ΔE = 91.2 K, g_u/g_l = 2, γ_H = 8e-10 T₂^0.07, γ_e = 2.8e-7 T₂^−½. O I 63 µm: A = 8.9e-5, ΔE = 227.7 K, g_u/g_l = 3/5, γ_H = 9.2e-11 T₂^0.67. Two-level, collisions with H (+e for C⁺). | HM89 via KROME Table B, secondary. O I 145 µm and the electron/H₂ colliders are omitted. |
| C11 | Photoelectric heating: 1.3e-24 ε G₀ n_H (D/D_MW), with ε = 4.9e-2/[1 + 4e-3 ψ^0.73] + 3.7e-2 (T/1e4)^0.7/[1 + 2e-4 ψ] and ψ = G₀T^½/(n_e φ_PAH) | Wolfire et al. 2003 eqs. 23–24, primary. Grain recombination cooling (eq. 25) is omitted. |

## 3. Independent reference values

These were computed from the published equations by a stand-alone script
(not from CHUI code) and asserted in
`tests/unit/test_chemistry_rate_references.cpp` (NOT RUN).

| Quantity | Value |
|---|---|
| α_H+ KWH96 (1e4 K) | 5.0971216268e-13 |
| α_He++ KWH96 (1e4 K) | 2.0388486507e-12 |
| α_He+ radiative KWH96 (1e4 K) | 4.3141104433e-13 |
| α_d dielectronic (1e5 K) | 6.1051443167e-13 |
| Γ_eH0 / Γ_eHe0 / Γ_eHe+ (2e4 K) | 2.1397210736e-12 / 1.4810844953e-15 / 1.0743770256e-23 |
| Γ_eH0 / Γ_eHe0 / Γ_eHe+ (1e5 K) | 1.9088405716e-09 / 2.1694570885e-10 / 1.6243646986e-12 |
| Ferland case A / case B (1e4 K) | 4.2969509300e-13 / 2.5918095507e-13 |
| GP98 H3 (1e3 K) | 8.0042332219e-16 |
| GP98 H5 (100 K / 1e3 K) | 1.5e-9 / 1.2361181730e-09 |
| GP98 H8 (1e3 K / 3.2e4 K, held above) | 4.8607132352e-18 / 7.0225570932e-16 |
| GP98 H7 (1e3 K) | 2.4078053118e-07 |
| Abel k12 (3e3 K) / k18 (1e3 K) | 5.0441202456e-24 / 6.9274584753e-09 |
| GA08 R9 k_v=0 / k_LTE (3e3 K) | 8.3630042370e-20 / 1.5536984947e-15 |
| GA08 n_cr,H / n_cr,H₂ / n_cr,He (3e3 K) | 1.3431207062e+03 / 9.2822637427e+05 / 1.0392638956e+05 |
| GA08 n_cr for x_H = 0.9, x_H₂ = 0.1, x_He = 0.08 (3e3 K) | 1.4904045741e+03 |
| GA08 blended k at n = 1e-2 / 1e3 / 1e8 (same mixture) | 8.3635558213e-20 / 4.3305534887e-18 / 1.5534708933e-15 |
| b5 thermal (100 K) | 9.0828328526e-01 |
| DB96 f_sh (N_H₂ = 1e14 / 1e18, 100 K) | 6.8006009867e-01 / 7.5343554122e-04 |
| Photoelectric heating (T = 100, n_H = 30, n_e = 0.03, G₀ = 1.7, D = 1) | 2.0147474478e-24 erg cm⁻³ s⁻¹ |
| O I 63 µm (T = 100, n_OI = 1e-3, n_H = 100) | 1.7801196802e-26 erg cm⁻³ s⁻¹ |
| C⁺ 158 µm (T = 100, n_C+ = 1e-4, n_H = 100, n_e = 0.01) | 7.8858398320e-26 erg cm⁻³ s⁻¹ |
| KWH96 HeII-only radiative cooling + Compton (T = 3e4, n_e = 1, n_He+ = 0.1, z = 0) | 1.0819664145e-25 erg cm⁻³ s⁻¹ (the removed Cen term would add 8.99e-30) |

## 4. Physical regimes that remain unsupported

- **Optically thick gas.**
  - Case-B recombination is not implemented.
  - Shielding is the local Rahmati or Jeans-length approximation, not
    radiative transfer.
  - σ̄ is fixed at the HM01 value.
- **H₂ destruction.**
  - H₂–H₂ and H₂–He collisional dissociation are absent.
  - Three-body H₂ formation is absent, so n_H > 1e8 cm⁻³ is out of range.
  - H⁻ photodetachment and H₂⁺ photodissociation are absent.
- **Dust H₂ formation rate (R18).** Unverified, with no grain-temperature
  dependence.
- **Missing species and cooling.**
  - No HD chemistry or HD cooling.
  - No CO cooling and no neutral atomic carbon.
  - No CO line self-shielding.
- **Metal ionization.**
  - Metal electron donors other than C⁺ are absent.
  - The metal-table ionization balance is that of the source table
    (WSS09: equilibrium, optically thin, HM01).
  - The table is limited to n_H ≤ 1 cm⁻³ and T ≥ 100 K. Outside that range
    CHUI refuses, or clamps density and redshift by policy, but never
    temperature.
- **H₂ self-shielding.** Uses the thermal Doppler parameter only.

## 5. UV background and metal-table compatibility policy

- **UVB declaration.** A `chui_metal_cooling_table_v1` table must declare
  `uv_background`, either `none` or a text identity. It may optionally
  declare `uv_background_treecool_sha256`.
- **Compatible pairings** (`MetalCoolingTable::checkUvBackgroundCompatibility`,
  called by `loadBaryonicDataProducts`):
  - the table declares `none` and CHUI has no UVB; or
  - the table's TREECOOL SHA-256 equals the configured TREECOOL file's
    SHA-256.
- **Mismatches.** Every other pairing is a mismatch, for example WSS09/HM01
  with FG20, or a photoionized table with no CHUI UVB. A mismatch is
  **refused** unless
  `physics.metal_cooling_accept_uv_background_mismatch = true`.
  - The key defaults to false and has typed validation.
  - Setting it without a table mode is refused as a no-op.
  - It is part of the normalized dump and the config hash, so an accepted
    mismatch is recorded in the run provenance.
- **No silent mixing.** The H/He photo-rates and heating always come from
  the configured TREECOOL. The metal-line increments always come from the
  table. CHUI never rescales one to the other.
- **No double counting.** The table's H/He (`Metal_free`) contribution and
  Compton term are never converted (the converter excludes them), and the
  reader refuses H/He columns.
- **Recommendation.** For FG20 runs, use a metal table computed with FG20,
  for example a future Ploeckinger & Schaye 2020 conversion (blocked; see
  `tools/chui_data/data_sources.md`). Alternatively, accept the WSS09/HM01
  mismatch explicitly and quantify it in validation.

## 6. Owner validation still required (USER-RESERVED)

- Compile, then run:
  - `unit_chemistry_rate_references`;
  - `unit_thermochemistry_backends`;
  - `unit_config_part2_baryonic_physics`;
  - the Part 2 suites.
- Equilibrium ionization fractions versus an independent CIE solver
  (e.g. KWH96 Fig. 1).
- H₂ formation and destruction against a GA08-like one-zone collapse.
- Shielded versus unshielded photodissociation.
- Metal-table interpolation against the WSS09 `Total_Metals` column at
  solar abundances. With the electron correction disabled, the sum over
  the nine elements should agree with `Total_Metals/Net_cooling` within
  the WSS09-stated accuracy.
- A decision on lifting the gate for `equilibrium_primordial`.
