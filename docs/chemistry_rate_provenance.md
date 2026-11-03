# CHUI thermochemistry rate provenance (corrective gate, 2026-11-03)

**Status: every coefficient used by the `equilibrium_primordial`,
`fixed_hybrid` and `fixed_detailed` backends is UNQUALIFIED.**
The coefficients were transcribed from the cited publications without a
documented line-by-line check. An attempt to check them against the
originals during the corrective gate failed because the cloud session's
network policy blocks `arxiv.org` (DNS ENOTFOUND). No coefficient was
changed to "look right" and none was invented.

Consequences enforced in source:
* `RateVerification::kUnverifiedTranscription` on every reaction and on the
  cooling fits; `NetworkManifest::qualified()` is false.
* Those backends refuse to run (config validation and solver construction)
  unless `physics.thermochemistry_allow_unqualified_rates = true`; runs with
  the opt-in are explicitly unqualified research runs, not
  production-validated science.
* `analytic_reference` (REFERENCE ONLY, closed-form) needs no rate data.
* Compton coefficient: computed from CODATA 2018 constants
  (4 sigma_T a k_B T_cmb^4 / (m_e c), T_cmb = 2.7255 (1+z) K, Fixsen 2009) —
  the only term derived without transcription.

## Owner verification procedure (required to qualify a backend)

For each row below: open the cited source, compare coefficient, exponents,
temperature units (K vs eV) and piecewise ranges; record the equation/table
identifier and the check (reviewer, date) in this file; then set the entry's
`verification` to `kVerifiedAgainstSource` in `src/physics/chemical_network.cpp`.
Use the numerical reference values (computed from the IMPLEMENTED formulas)
to compare against an independent evaluation of the published fit.

| ID | Reaction | Implemented form | Claimed source (to verify) | Regime notes / risks |
|---|---|---|---|---|
| R01 | H + e → H⁺ + 2e | 5.85e-11 T^½ e^{-157809.1/T} (1+√(T/1e5))⁻¹ | Katz, Weinberg & Hernquist 1996 (KWH96) Table 2, after Cen 1992 | fit for T ≳ 1e4 K; negligible below |
| R02 | H⁺ + e → H (case A) | 8.4e-11 T^-½ (T/1e3)^-0.2 [1+(T/1e6)^0.7]⁻¹ | KWH96 Table 2 | gives 5.10e-13 at 1e4 K vs 4.18e-13 (Osterbrock & Ferland 2006 case A): +22%, consistent with a fit, but unverified; case A vs B choice matters for optically thick gas |
| R03 | He + e → He⁺ + 2e | 2.38e-11 T^½ e^{-285335.4/T} (1+√(T/1e5))⁻¹ | KWH96 Table 2 | |
| R04 | He⁺ + e → He | 1.5e-10 T^-0.6353 + 1.9e-3 T^-1.5 e^{-470000/T}(1+0.3e^{-94000/T}) | KWH96 Table 2 (radiative + dielectronic) | |
| R05 | He⁺ + e → He⁺⁺ + 2e | 5.68e-12 T^½ e^{-631515/T} (1+√(T/1e5))⁻¹ | KWH96 Table 2 | |
| R06 | He⁺⁺ + e → He⁺ | 3.36e-10 T^-½ (T/1e3)^-0.2 [1+(T/1e6)^0.7]⁻¹ | KWH96 Table 2 | |
| R07–R09 | photoionization | Γ n from the user TREECOOL table × Rahmati et al. 2013 eq. 14 factor | TREECOOL provenance + Rahmati 2013 | shielding fit applied to He as an approximation |
| R10 | H + e → H⁻ + γ | 1.4e-18 T^0.928 e^{-T/16200} | Galli & Palla 1998 (GP98) | |
| R11 | H⁻ + H → H₂ + e | 1.5e-9 (T<300), 4.0e-9 T^-0.17 | GP98 | |
| R12 | H + H⁺ → H₂⁺ + γ | 10^(-19.38 - 1.523 lgT + 1.118 lgT² - 0.1269 lgT³), T clamped to ≤3.2e4 K | GP98 | the 3.2e4 K clamp is an implementation choice, not a verified limit |
| R13 | H₂⁺ + H → H₂ + H⁺ | 6.4e-10 | Karpas, Anicich & Huntress 1979 | constant |
| R14 | H⁻ + H⁺ → 2H | 5.7e-6 T^-½ + 6.3e-8 - 9.2e-11 T^½ + 4.4e-13 T | GP98 | grows above ~1e5 K (4e-6 at 1e7 K): fit used beyond its range; H⁻ is negligible there but the range must be verified |
| R15 | H₂ + e → 2H + e | 5.6e-11 T^½ e^{-102124/T} | Donahue & Shull 1991 via Abel et al. 1997 k12 | |
| R16 | H₂ + H → 3H | **REMOVED** | (Dove & Mandy 1986 via Abel et al. 1997 k13) | transcribed form decreased above 1e4 K and vanished by 1e6 K — defect; channel omitted, H₂ destruction in warm atomic gas underestimated |
| R17 | H₂⁺ + e → 2H | 1.0e-8 (T<617), 1.32e-6 T^-0.76 | Abel et al. 1997 k18 | |
| R18 | 2H (grain) → H₂ | 3e-17 (D/D_MW) n_H n_HI | Jura 1975 | no grain-temperature / sticking dependence |
| R19 | H₂ + γ → 2H | 3.3e-11 G0 f_sh(N_H₂) e^{-3.74 A_V} | Glover & Mac Low 2007; Draine & Bertoldi 1996 shielding | local Jeans-length column (approximation) |
| R20 | C⁺ + H₂ (+O,e) → CO + 2H | k0=5e-16, k1=5e-10, Γ_CHx=5e-10 G0 | Nelson & Langer 1997 | NL97 reduced network; all non-CO C as C⁺ |
| R21 | CO + γ → C⁺ + O + e | 1e-10 G0 e^{-2.5 A_V} | Nelson & Langer 1997 | no CO / H₂ line self-shielding |

Cooling fits (same UNQUALIFIED status): KWH96 Table 1 primordial terms
(collisional excitation, ionization, recombination, dielectronic,
free-free with Spitzer & Hart 1979 Gaunt fit); H₂: Hollenbach & McKee 1979
LTE + GP98 low-density fit (10–1e4 K); C⁺ 158 µm and O I 63 µm two-level
atoms (Hollenbach & McKee 1989 collision rates); photoelectric heating
(Bakes & Tielens 1994 / Wolfire et al. 2003).

## Numerical reference values of the implemented formulas [cm³ s⁻¹]

Computed by evaluating the formulas exactly as implemented (scratch
arithmetic; not a project build or test). Compare each against the
published fit.

| Rate (implemented formula) | T=1e+02 K | T=1e+03 K | T=1e+04 K | T=1e+05 K | T=1e+06 K | T=1e+07 K |
|---|---|---|---|---|---|---|
| R01 collisionalIonizationHI | 0.0000e+00 | 4.8994e-78 | 6.2268e-16 | 1.9088e-09 | 1.2003e-08 | 1.6554e-08 |
| R03 collisionalIonizationHeI | 0.0000e+00 | 8.2337e-134 | 7.3331e-22 | 2.1695e-10 | 4.2986e-09 | 6.6496e-09 |
| R05 collisionalIonizationHeII | 0.0000e+00 | 8.9018e-285 | 1.6168e-37 | 1.6244e-12 | 7.2569e-10 | 1.5330e-09 |
| R02 recombinationHII (A) | 1.3292e-11 | 2.6354e-12 | 5.0971e-13 | 8.8160e-14 | 1.0550e-14 | 7.0028e-16 |
| R04a recombinationHeII | 8.0444e-12 | 1.8629e-12 | 4.3141e-13 | 9.9906e-14 | 2.3136e-14 | 5.3579e-15 |
| R04b dielectronicHeII | 0.0000e+00 | 4.5745e-212 | 7.3608e-30 | 6.1051e-13 | 1.5118e-12 | 7.4361e-14 |
| R06 recombinationHeIII | 5.3168e-11 | 1.0542e-11 | 2.0388e-12 | 3.5264e-13 | 4.2200e-14 | 2.8011e-15 |
| R10 formationHminus | 9.9873e-17 | 8.0042e-16 | 3.8909e-15 | 1.2744e-16 | 8.0506e-40 | 3.6234e-280 |
| R11 associativeDetachment | 1.5000e-09 | 1.2361e-09 | 8.3572e-10 | 5.6502e-10 | 3.8200e-10 | 2.5826e-10 |
| R12 formationH2plus | 1.0735e-19 | 4.8607e-18 | 1.9697e-16 | 7.0226e-16 | 7.0226e-16 | 7.0226e-16 |
| R13 chargeTransferH2plus | 6.4000e-10 | 6.4000e-10 | 6.4000e-10 | 6.4000e-10 | 6.4000e-10 | 6.4000e-10 |
| R14 mutualNeutralization | 6.3212e-07 | 2.4078e-07 | 1.1520e-07 | 9.5932e-08 | 4.1670e-07 | 4.1739e-06 |
| R15 electronDissociationH2 | 0.0000e+00 | 7.8759e-54 | 2.0559e-13 | 6.3778e-09 | 5.0563e-08 | 1.7529e-07 |
| R17 dissociativeRecombinationH2plus | 1.0000e-08 | 6.9275e-09 | 1.2039e-09 | 2.0921e-10 | 3.6356e-11 | 6.3179e-12 |

## Missing channels and their consequences

Listed in each manifest's `omitted_channels`: H⁻ photodetachment and H₂⁺
photodissociation (overestimated H₂ formation under strong radiation), H⁻ +
e collisional detachment and H₂ + H⁺ (minor above 1e4 K), three-body H₂
formation (invalid above n_H ~ 1e8 cm⁻³), H₂ + H collisional dissociation
(removed, see R16), HD chemistry/cooling (primordial gas below ~200 K
cools too slowly), metal electron donors except C⁺, grain recombination
cooling, neutral C, CO cooling (dense molecular gas cools too slowly),
CO/H₂ line self-shielding (CO photodissociated too easily).
