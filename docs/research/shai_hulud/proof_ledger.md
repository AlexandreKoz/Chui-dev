# SHAI-HULUD proof ledger

Status of every mathematical statement in `manuscript/shai_hulud.pdf`. Numbering follows the manuscript.

Categories

- **Established** — known result; reference given in the manuscript; not re-proved (a short verification may be included).
- **New derivation** — proved in the manuscript. This is *not* a priority claim: the literature check for this note was limited.
- **Numerical observation** — supported only by computation in `prototype/`; no proof.
- **Conditional** — proved under a hypothesis that is itself not proved in general.
- **Open** — neither proved nor refuted here.

## A. Established results (used, not claimed)

| Ref. | Statement | Source |
|---|---|---|
| Prop. 2.2 | A flow-transported sheet `F_t#μ` is a weak solution of Vlasov | chain rule; standard |
| Prop. 2.4 | `F_t*ω` constant; Poincaré loop invariants; conserved Liouville class; coordinate form = Cauchy invariants | Cartan formula; Arnold; Frisch–Villone / Zheligovsky–Frisch for the Cauchy-invariant reading |
| Thm 2.6 (a)(b) | Action transport: `dΣ_t = F_t*θ`; the discrete action is the generating function of the step map | Hamilton–Jacobi theory; Marsden–West; Hairer–Lubich–Wanner |
| Thm 2.10 | The Lagrangian condition is the zero level of the momentum map of mass-preserving relabelling | Donaldson (1999); fluid analogue: Salmon |
| Nonlinear form of Prop. 2.12 | Nearby Lagrangians are graphs of closed one-forms; same Liouville class ⇔ exact | Weinstein (1971) |
| Thm 3.2 | Existence, uniqueness, gradient structure and regularity of the optimal map on the torus; Monge–Ampère equation | Brenier; Cordero-Erausquin; McCann; Caffarelli |
| Prop. 3.3 | d = 1: monotone rearrangement, one-parameter family of offsets | classical; Delon–Salomon–Sobolevski |
| Second-order accuracy of KDK | time discretisation error | standard |

## B. New derivations (proved in the manuscript)

| Ref. | Statement | Proof | Independent check |
|---|---|---|---|
| Thm 2.6 (c) | Generating functions of drift and kick: `A = D(p·∇K − K)`, `A = −κV`; exact per-vertex action update for any time-dependent KDK | 3 lines | Poincaré period at round-off in all runs, including time-dependent Hamiltonians |
| Ex. 2.8 | Zel'dovich data are exactly Lagrangian with explicit action; 2LPT data have defect `a²(D₁Ḋ₂ − D₂Ḋ₁)[∇∇χ₁, ∇∇χ₂]` | direct computation | none (algebra only) |
| Thm 2.9 | Relabelling invariance of measure, Lagrangian property and periods; ALE system for moving labels | elementary | — |
| Prop. 2.12 | First-order characterisation: a perturbation preserves the structure iff its normal error one-form `β_e = ω(e, dF)` is exact | Cartan formula for maps | consistent with all experiments |
| Prop. 2.13 | Cubic form `C_abc = ω(∂_a F, ∂_bc F)` is a totally symmetric tensor, affine-symplectic invariant | proved | symbolic (graph case) |
| Prop. 3.4 | Sampling is not mass: pull-back density keeps `f`; masses attached to moved samples give an O(1) error bounded by `Lip(F)‖R − id‖` | elementary | cell-mass ratios observed |
| Lemma 3.5 | Edge action defect `ε_e = −C(e,e,e)/12 + O(h⁴)` | Taylor expansion | sympy (`edge_defect_p1_leading_term`) |
| Eq. (3.9)–(3.10) | Monitor exponents; equi-affine sampling; complexity functional `𝒜`, gain `G ≥ 1` | variational argument; Hölder | measured `𝒜(t)`, `G(t)` |
| Prop. 4.2 | P1 interpolant of a Lagrangian sheet has flux `−¼[C(u,u,v) − C(u,v,v)] + O(h⁴)` per triangle | from Lemma 3.5 | sympy, generic quartic potential, and the explicit example |
| Prop. 4.3 | Total mass is the degree of the label map: exact for any conforming label interpolation | degree formula | 0 (d = 1) and ≤ 1.1e-16 (d = 2) in all runs |
| Thm 5.1 | Carrier transport: with a common-potential KDK the nodal data stay samples of an exactly transported carrier; all invariants exact between remaps | from Thm 2.6 | oracle construction; fixed-mesh runs |
| Thm 5.4 | Invariant budget of an interpolatory remap, with no loss of derivative | integration by parts | drift ∝ number of remaps for linear/cubic/quintic remaps |
| Cex. 5.5 | Midpoint remap of the regular N-gon multiplies the action by `cos²(π/N)` while conserving mass and momentum; any inscribed re-sampling of a strictly convex polygon loses action | elementary | numeric, N = 8, 64, 512 |
| Thm 5.6 | Rigidity: exactly isotropic P1 sheets over the Kuhn-triangulated lattice (d = 2) are `(f(i₁), g(i₂)) + k(i₁−i₂)(1,−1)`; dimension `3n − 1` of `2n²` | discrete Fourier analysis | SVD kernel dimension `3n` (incl. the Σ constant) for n = 5, 6, 8 |
| Thm 5.8 (i),(ii),(iv) | Action-corrected interpolation AC(m,1): continuity, **exact edge actions and exact period**, m = 2 is the cubic Hermite interpolant of the generating function (parabolic arc, sub-arc area `ε s³`) | area-variation identity (5.7), linear in the bubble amplitude | period error ≤ 6e-15 in all runs; sympy for the `s³` law |
| Thm 5.9 | Generating-function elements: any C¹ interpolant of `W` with nodal data `(W_i, DXᵀP_i)` gives a continuous, exactly Lagrangian carrier in a chart; error bound | two lines | d = 2 circulation at quadrature round-off |
| Eq. (5.8) | Quadratic form `Q_M` relating generating functions in two linear symplectic frames | direct computation | — |
| Thm 6.1 | Structure preservation of the coupled scheme (mass; d = 1 exact; d ≥ 2 exact chart-wise) | combination of the above | experiments |
| Lemma 6.2 | Cubic Hermite front: `‖W_H′‖∞ ≤ 2‖W′‖∞`, sharp | elementary | numeric maximum = 2 |
| Prop. 8.1 | Comoving Hamiltonian: same statements with `A_drift = D|p|²/2`, `A_kick = −κφ_n` | corollary | external-potential EdS test |

## C. Proved only to leading order, or conditional

| Ref. | Statement | What is missing |
|---|---|---|
| Thm 5.8 (iii) | AC(m,1): tangential error `O(h^m)`, normal error `O(h^{m+1})` with explicit leading constant (includes the hairpin factor `‖F″‖/min|F′|`) | higher-order terms not tracked with explicit constants. Measured: order 5.00 for m = 4 on exact data, ≈ 7 for m = 6 in a development test |
| Thm 6.4 | Convergence for a bounded number of remaps, `error ≤ Σ Λ^{N_r−k} e^{ℓ(t−t_k)} C h^r |F∘T_k|_r` | hypothesis (R2) (Lipschitz stability of the reconstruction with respect to the carrier) is proved only for plain Lagrange interpolation and for AC(2,1) in a frozen frame (Λ = 2); assumed for AC(m,1) with moving chord frames |

## D. Numerical observations (no proof)

| Ref. | Observation | Evidence |
|---|---|---|
| Obs. 5.7 | Generic (jittered) vertex positions admit only trivial exactly isotropic P1 sheets (kernel dimension 3 in d = 2, 4 in d = 3); regular d = 3 Kuhn lattice: `6n − 2` | SVD ranks, n = 5, 6, 8 (d = 2), n = 4, 5 (d = 3) |
| Obs. 6.3 | Spectral radius of the linearised AC(2,1), AC(4,1), AC(6,1) remap ≤ 1 + 9e-16 for all shifts and wavenumbers | 101 × 240 sample grid; necessary condition only |
| §7.1 | Errors of AC(4,1) nearly independent of the number of remaps (1–399) | `exp1d_accumulation_twist.json` |
| §7.1 | Quadratic-element edge action error is `O(h⁵)` | measured orders 4.93, 4.98, 5.00 |
| Rem. 5.2 | Bilinear interpolation of a centred-difference grid force has pointwise relative curl 0.43, 0.22, 0.11 at n = 16, 32, 64 | `verify_propositions.py` |
| Rem. 3.6 | Remap order collapses to ≈ 1 when the monitor is smoothed only at the mesh scale; restored by fixed-scale mollification | d = 2 development runs (described, final runs use the mollified monitor) |
| §7.2 | In d = 2 exact Lagrangian reconstruction costs about one order of momentum accuracy versus an unconstrained spline of equal degree | `exp2d.json` |

## E. Open problems and unverified items

1. **Global conforming exactly Lagrangian carrier across caustics, d ≥ 2** (Open problem 5.11). No construction; the method is exact chart-wise and falls back on the budget of Thm 5.4 at chart transitions.
2. **Uniform stability under unboundedly many remaps** (power-boundedness). Evidence only (Obs. 6.3, accumulation test).
3. **(R2) for AC(m,1)** with moving chord frames; explicit constants in Thm 5.8 (iii).
4. **Self-consistent Vlasov–Poisson**: existence/regularity of the cold sheet solution after shell crossing in d ≥ 2 (assumed, Assumption 2.3); stability and convergence of the coupled scheme; whether exact sheet invariants improve self-gravitating runs. No self-gravitating computation was made.
5. **AC(m, m−1)** (all stencil cell actions): formally of normal order `2m − 1`; observed irregular convergence at hairpin turns; not analysed.
6. **Anisotropy**: optimal transport equidistributes a scalar; alignment with the cubic form is not addressed.
7. **Novelty**: only the abstract of the ColDICE paper could be fetched in this session; details of ColDICE beyond its abstract, and the absence of an equivalent action-based remap in the geometric-optics (wavefront construction) and geometric-integration literature, are unverified.
