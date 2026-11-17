# Local refinement contract — locally h-adaptive high-order phase sheet (2D Vlasov–Poisson)

Scope: the research solver in `src/shadapt/` (this directory). Evidence tags as in the dossier:
PROVED (argument given here, elementary), MEASURED (number in `results/` or asserted by a test in
`tests/`), CONJECTURE, FAILED, OPEN. Nothing here is peer reviewed. Nothing here is part of CHUÍ.

## 0. Frozen conventions (unchanged from the Fable breakthrough package)

Label torus `ξ ∈ T² = [0,1)²`, sheet `F(ξ) = (x(ξ), p(ξ))`, `x(ξ) = ξ + d(ξ)` with `d` periodic
(degree-one map), specific momentum `p = dx/dt`, `θ = p·dx`, action `Σ` with periods `c`.
Box 1, total mass 1, `4πG = 1`: `Δφ = ρ − 1`, static background, time variable `t`.
KDK with one field solve per step after the drift, `dt = 0.01`, 300 steps to `t = 3`.
PM: PCS assignment, continuous Green function `−1/k²`, exact-gradient ("grad") force law
`a(x) = −∇V(x)`, `V(x) = Σ_g φ_g W(x − x_g)`. Initial data `GrowingModeIC(eps=(0.10,0.07,0.03), phase=0.9)`.
Mass measure: Lebesgue measure of the **label** torus; a label cell of area `A` has mass `A`.
The original code was re-run before any change: 30/30 tests; the 512² control at PM 64² and four
saved arms reproduce **bit-identically** (`results/exp00_baseline_repro.json`).

## 1. Topology

Level `l` has `N_l = n0·2^l` cells per dimension and `N_l/B` blocks per dimension; block `(l,I,J)`
covers cells `[IB,(I+1)B)×[JB,(J+1)B)` (`B = 8`). All level-0 blocks exist. A block may be *split*
into its four children. Blocks that exist and are not split are *leaves*.

Invariants (all enforced by `Forest.check()`, exercised on random forests in `tests/test_forest.py`):

* **nesting** — a level-`l≥1` block exists iff its parent exists and is split;
* **2:1 balance** including corners and the periodic seam — every 8-neighbour of a split block
  exists at the same level. A split request first splits coarser neighbours; those splits are
  recorded with reason `"balance"` and the id of the block that caused them;
* **exact tiling** — leaf areas (dyadic numbers) sum to exactly 1; every finest-level slot is
  covered by exactly one leaf (`leaf_level_map` asserts multiplicity 1).

**Nodes and ownership.** The active nodes of level `l` are the level-`l` grid points in the
*closed* region covered by level-`l` blocks that are not level-`(l−1)` grid points (at least one
odd index); level 0 owns its whole grid. Hence (PROVED, trivial; tested) every active node has
exactly one owner level and one integer key `i·N_l + j`; nodes on edges shared by blocks are stored
once; a node on a coarse–fine interface (a "hanging" node) is an ordinary active node with its own
degrees of freedom. **No hanging-node constraint is needed** because the carrier below is globally
smooth by construction (§2); this is the "documented alternative" to trace matching.

Lineage: every non-root block records parent, creation step/time and reason; every split/merge is
appended to `Forest.events` with the indicator values that caused it. `Forest.to_dict()` /
`from_dict()` round-trip the topology and the history (checkpoint format).

## 2. The composite hierarchical-surplus carrier

Let `β = β^k` be the centred cardinal B-spline of odd degree `k ∈ {3,5}` and
`β_{m,j}(ξ) = β(N_m ξ₁ − j₁) β(N_m ξ₂ − j₂)`. Let `Ω_m` be the closed region covered by level-`m`
blocks and `T_m` the set of level-`m` grid points in `Ω_m` dilated by `h_t` rings of blocks
(`h_t = 1` by default: a halo of `B` fine cells). For periodic nodal data `z` (each component of
`d`, `p`, `Σ − c·ξ` separately):

```
F = Σ_{m=0}^{top} S_m ,        S_m = Σ_{j∈T_m} c_{m,j} β_{m,j} ,
S_0 = periodic interpolating spline of the level-0 nodes              (the baseline carrier),
S_m (m≥1): the unique element of span{β_{m,j}: j∈T_m} with
      S_m(ξ_j) = z_j − Σ_{l<m} S_l(ξ_j)    at the active level-m nodes       (the "surplus"),
      S_m(ξ_j) = 0                          at every other grid point of T_m.
```

**P1 (well posed, SPD, conditioned) — PROVED.** The 1D periodic collocation matrix `A_N` of `β` at
the knots is a symmetric circulant with symbol `λ(θ) = Σ_j β(j) e^{ijθ} ∈ [λ_min, 1]`,
`λ_min = λ(π) = 1/3` (cubic), `2/15` (quintic). The matrix of the level-`m` problem is the
principal submatrix of `A_{N_m} ⊗ A_{N_m}` on `T_m`; by Cauchy interlacing it is symmetric positive
definite with spectrum in `[λ_min², 1]`: condition number `≤ 9` (cubic), `≤ 56.25` (quintic), for
**every** shape of `T_m`. `F` is a linear function of the nodal data.

**P2 (regularity) — PROVED and MEASURED.** Every `S_m` is a finite sum of `C^{k−1}` tensor
B-splines, so `F ∈ C^{k−1}(T²)`: **C⁴ for the quintic carrier, C² for the cubic**, across every
coarse–fine interface and everywhere else; it is *not* `C^k`. Measured two-sided jumps of the
normal derivatives across all 96 hanging-interface probes of a two-level patch (exact
polynomial-fit probe, no difference-quotient truncation; `diagnostics.continuity_jumps`,
`test_interface_continuity_is_k_minus_1`), relative to the size of the derivative:

| order of normal derivative | 0 | 1 | 2 | 3 | 4 | 5 |
|:--|--:|--:|--:|--:|--:|--:|
| quintic | 5e-15 | 4e-12 | 5e-10 | 5e-8 | 2e-6 | **0.5** |
| cubic | 2e-15 | 5e-13 | 2e-11 | **0.5** | | |

(The growth with the order is the conditioning of the probe, not a defect: it is the same on
interior block edges, `kind="all"`.) In the self-gravitating runs the same probe is applied at
`t = 3` to every locally refined arm (`interface_jumps` in the result files).

**P3 (interpolation) — PROVED up to the skirt leak; leak MEASURED.** Let `x` be an active node of
level `l`. For `m < l` the terms are in the definition of the surplus; `S_l(x)` is the surplus;
for `m > l`, `x` is a level-`m` grid point that is not an active level-`m` node, so `S_m(x) = 0`
if `x ∈ T_m` (a condition of the system) or if `x` is more than `(k−1)/2` grid steps outside `T_m`
(outside the support). In the remaining *skirt* (`≤ (k−1)/2` steps outside `T_m`) `S_m(x)` is a
combination of the outermost coefficients, which are `≥ h_t·B` steps away from every non-zero of
the right-hand side. By the Demko–Moss–Smith estimate for inverses of banded SPD matrices they are
bounded by `C q^{h_t B/h}·max|surplus|`, `q = (√κ−1)/(√κ+1)`, `h = (k−1)/2`: a guaranteed rate
per grid step of `0.5` (cubic) and `0.875` (quintic). The **measured** rate is the root of the
cardinal spline, `λ = 0.268` (cubic), `0.431` (quintic): the interpolation defect at any node is
`≤ 2 λ^{h_t B}·max|surplus|` (`test_interpolates_every_active_node`: `≤ 2·λ⁸` for `h_t = 1`,
`≤ 2·λ¹⁶` for `h_t = 2`, i.e. `1.2e-3` resp. `1.4e-6` of the surplus for quintics; in absolute
terms `1e-11…1e-15` in the tests). Nodes inside `Ω_m` are interpolated exactly by `Σ_{l≤m} S_l`.

**P4 (uniform limit) — PROVED and MEASURED.** If all blocks of levels `1..m` exist then
`T_l` is the whole torus and, by induction, `S_l = I_l[z] − Σ_{j<l} S_j` (`I_l` the periodic
interpolation operator of level `l`, a projector on the nested spline space), so
`F = I_m[z]`: the composite carrier **is** the baseline periodic spline of all nodes on the fine
grid. Measured: `2e-15` against `shsheet.pspline.PSpline2` two levels deep; the coupled solver with
everything split evolves like the baseline uniform sheet to `1e-12`
(`test_fully_refined_sheet_matches_baseline_fine_sheet`), and the unrefined solver reproduces the
baseline `VPRun` cubic and quintic arms to `1e-13` over 40 steps.

**P5 (transfer = identity) — PROVED up to the skirt leak; MEASURED.** A split creates nodes with
values `F_old(ξ_new)` for `(d, p, Σ)`; old nodes are kept bit for bit. The old corrections satisfy
the new systems except at grid points of `T_m^new \ T_m^old` in the skirt of the old `T_m`, so
`F_new = F_old` up to the leak of P3: the represented sheet, its density, momentum, kinetic and
potential energy and every loop action are unchanged by a split. Measured: `‖F_new − F_old‖∞ ≤
2 λ⁸·surplus` (1e-11 in the test), density change `< 1e-8 ρ_max`, `|ΔP|, |ΔK|, |ΔU| < 1e-10 K`
(`test_split_leaves_the_carrier_unchanged`, `test_mid_run_split_conserves_…`). Mass is exact
(§4). No exact/reference value is ever used: the only inputs are the old nodal state.
A merge drops the surplus of the removed nodes (an `O(tol)` change of the carrier by construction
of the hysteresis rule) and conserves mass exactly; it is **not** the identity.
*Re-basing*: when every level-0 block is split, levels 0 and 1 are relabelled as a new level 0
(`n0 → 2 n0`). By P4 this is exact (measured `1e-13`); it returns the cost of the carrier to that
of a uniform sheet.

**P6 (order) — MEASURED.** Same refinement pattern on `n0 = 16, 32, 64` (all `h` halve):

| region | cubic order | quintic order |
|:--|--:|--:|
| fine interior | 4.0 | 6.2 |
| ±¼ coarse cell around the coarse–fine interface | 4.0 | 6.2 |
| coarse region | 4.1 | 6.3 |

i.e. order `k+1` everywhere including the interface strip; near an interface the error is that of
the **coarse** level (the correction is forced to vanish outside `Ω_m`), decaying to the fine
level's error at the rate `λ` per fine cell into the patch.

**Locality — PROVED (structure) and MEASURED.** `S_m` is supported in `T_m` plus `(k+1)/2` cells;
stored coefficient tiles of level `m` cover `Ω_m` plus three block rings (not the level-`m` grid);
a unit bump on a fine node changes the carrier by `< 1e-4` twelve fine cells away and by exactly 0
beyond 24, and never changes a coarser level (`test_correction_support_is_local`). One split of a
level-0 block creates `(2B+1)² − (B+1)² = 208` nodes and `3B²` leaf cells (`test_one_split_is_local`).

### 2.1 Evaluation (accumulated representation)

`A_m = Σ_{l≤m} S_l` is written in level-`m` coefficients, `a_m = P_m a_{m−1} + c_m`, with the
two-scale mask of `β` (`[1,6,15,20,15,6,1]/32` for quintics). `a_m` is stored on tile sets `D_m`.
A leaf block of level `L` evaluates from `A_{ℓe}`, `ℓe ≥ L` the finest level whose correction can be
non-zero on the block (`ℓe ∈ {L, L+1}` for `h_t = 1` under 2:1 balance). **Lemma (PROVED):** for
every `m > ℓe` the support of `S_m` (`T_m` + one tile) does not meet the block, so `F = A_{ℓe}` on
the block; `D_m` is built finest-to-coarsest so that every evaluation window and every subdivision
stencil is stored (checked at layout time; violations raise). Evaluation on the tensor sample
lattice of a block is two dense matrix products `W₁ · A · W₂ᵀ` on a `(B·2^{ℓe−L}+k)²` window.

### 2.2 Linear solver

On each grid line of `T_m` the 1D collocation operator is inverted exactly (banded Cholesky,
segments rotated so that none straddles the periodic seam; complete periodic lines by FFT).
**Lemma (PROVED):** `(L_i L_j)_{pq} = M_{pq}` iff the corner point `(i_q, j_p)` belongs to `T_m`;
hence `L_i L_j = M` on rectangles, bands and the full torus, and the defect `M − L_i L_j` couples
only pairs of points diagonally across a re-entrant corner of `T_m`, which lie `≥ h_t B − h`
steps from every non-zero of the right-hand side. The default solver applies **one** sweep
`c = L_j⁻¹ L_i⁻¹ r` (a fixed linear operator; no iteration, no tolerance); a level whose `T_m` is
the whole torus uses the 2D FFT prefilter. MEASURED on random forests
(`test_fixed_sweep_solver_matches_exact_solves`): coefficients within `3e-3·surplus` of a sparse-LU
solve, and the node interpolation defect is unchanged (it is the skirt leak of P3 either way).
Verification solvers `adi-tol` (sweeps to `1e-14`, CG fallback never triggered) and `splu` agree to
`1e-12`. Sparse LU is not used in production runs: its fill is ~20 MB per 10⁴ unknowns here.

## 3. Mass measure and quadrature

Every leaf block deposits its label mass `4^{-L} B²/n0²` exactly once, as `nq₁×nq₂` samples of
equal mass located on a tensor midpoint lattice of the block and placed at `x = ξ + d_h(ξ)` with the
composite carrier. `nq` are powers of two, so all masses are dyadic and the **global ledger is
exactly 1.0 in floating point** (asserted with `==`); the grid mass is `1 ± 2e-13`.
Mass belongs to label cells: no transfer, re-base or regrid reweights it.

* **Lattice mode (default).** `nq = B·max(1, sub0/2^L)`: for `sub0 ≥ 2^{top}` the samples of all
  levels form **one uniform label lattice** `(n0·sub0)²`, identical before and after any split.
  The density is then the rectangle rule on the torus for the periodic integrand
  `ξ ↦ W(x_h(ξ) − x_g)`, whose error is governed by the smoothness of `W∘x_h`, not by `O(δ²)`.
  Consequences MEASURED: the deposit of a cold sheet at rest is exactly uniform on any refinement
  pattern (`1e-12`); a split leaves `ρ` unchanged (`1e-9 ρ_max`); translation through the periodic
  seam by whole PM cells rolls `ρ` (`1e-10`).
* **Spacing-adaptive mode — FAILED (recorded).** Choosing `nq` per block from the largest sample
  spacing of its image (target in PM cells) changes the sample spacing across block edges. The
  Euler–Maclaurin edge terms `(δ₁² − δ₂²)/24 ∮ ∂_n f` of neighbouring blocks no longer cancel, the
  rule drops to second order at those edges and the changes in time inject noise. At PM 256²,
  quintic `n = 256`: `x_rms(t=3) = 1.2e-5` with 6.7e5 samples/step against `3.2e-7` with 2.6e5
  lattice samples (§4 of `results/SUMMARY.md`).
* **Gauss panels — alternative tested, also loses.** Two Gauss–Legendre points per panel (equal
  positive weights) are fourth order with or without uniformity and repair most of the damage of
  the adaptive rule, but a uniform Gauss lattice leaves an alias at the *panel* frequency with
  relative amplitude `|cos(π/√3)| = 0.24`, about 4× the error of the midpoint lattice at equal
  sample count (measured 3–6×), and the adaptive Gauss rule never beats the uniform midpoint
  lattice at equal samples. **Decision: uniform label lattice; quadrature is not adapted.**

Transient memory is bounded: blocks are evaluated and deposited in chunks of `2^15` samples.

## 4. Field coupling and what is (not) conserved

All nodes of all levels are advanced by KDK with `a = −∇V` gathered at their own positions; kicks
of consecutive steps are merged (one gather per step; equal to the half-kick sequence to
round-off). Each node is therefore an exact characteristic of the one-particle Hamiltonian
`|p|²/2 + V(x,t)`, and its action `Σ` is transported exactly along its own trajectory.

**This nodal force is not the gradient of the discrete energy.** With
`U_h = ½⟨ρ_h − 1, φ_h⟩` and `ρ_h` assembled from the carrier, `∂U_h/∂X_a = Σ_s m_s ℓ_a(ξ_s) ∇V(x_s)`
(`ℓ_a` the cardinal function of node `a`), and the consistent kinetic energy `½∫|p_h|²dμ` has a
non-diagonal Gram matrix. The scheme used here is the consistent, non-variational "tracers of the
common potential" discretisation; `E_h = K_h + U_h` is therefore **not** an invariant even for
`dt → 0`. Measured residuals at `t = 3` are tabulated per arm in `results/SUMMARY.md`
(mass ledger exactly 1; momentum drift `1e-7…1e-6` — the "grad" law has a self-force; relative
energy drift; carrier loop-action deviation from the conserved periods). The adjoint
(energy-derived) force is the subject of Call 2 and was deliberately not attempted here.

## 5. State-only indicators

For a degree-`k` spline the `(k+1)`-th undivided difference of the coefficients is the jump of its
`k`-th derivative across a knot, `[∂^k s]·h^k = Δ^{k+1}c`, and for smooth data the interpolation
error is `κ_k h^{k+1} f^{(k+1)}·b_k(t)` with `max|b_k| = 1`, `κ₃ = 1/384`, `κ₅ = 1/15360`
(MEASURED to 0.2 %, `test_error_constants_of_the_periodic_splines`). Per leaf block, in PM cells:

* `eta_d6 = κ_k 2^{k(ℓe−L)} max|Δ^{k+1} a| / Δ_pm` for the displacement and `τ_p ×` momentum — the
  truncation ("d4/d6") monitor of the block at its own resolution. It predicts the measured
  carrier error within a factor 2 on smooth sheets and its block map correlates at 0.95 (log–log)
  with the true error map of a folded sheet, but it measures the *instantaneous* interpolation
  error: at `t = 3` the accumulated trajectory error was ~25× larger (offline analysis).
* `eta_pred`: the same functional applied to the look-ahead `x + τ p + ½τ² a(x)` with the PM
  acceleration already gathered at the nodes — the predictive strain/force-informed candidate.
  For ballistic motion it equals `eta_d6` of the state at `t + τ` exactly
  (`test_predictive_indicator_is_the_d6_indicator_of_the_future_state`).
* strain diagnostics: `det Dx` and the ballistic time-to-caustic `τ_c` (smallest positive root of
  `det(Dx + τ Dp) = c₀ + c₁τ + c₂τ²`) at the cell centres; used for the gates `precaustic`
  (`τ_c ≤ τ_lead`) and `crossing` (`det Dx ≤ 0`).

Controller contract (`controller.Controller`): threshold `tol`; hysteresis (`merge` only if
`2^{k+1}·eta_children < hyst·tol`); hard node budget (never exceeded; each denied request is
rolled back and counted); cooldown; resolution floor in PM cells; 2:1 balance; optional buffer
ring; every split carries `reason ∈ {indicator:<name>, balance, buffer}` and the indicator values.
`test_indicators_use_only_the_state` runs the controller with `open` and `numpy.load` disabled.

## 6. Invariants and where they are tested

| Invariant | Status | Test / evidence |
|:--|:--|:--|
| unique leaf ownership, exact tiling, 2:1 incl. corners and seam | PROVED + tested | `test_forest.py` |
| unique deterministic node ownership, hanging nodes counted once | PROVED + tested | `test_forest.py` |
| constants / affine sheets exact on any topology | tested (`2e-15`) | `test_carrier.py` |
| level 0 ≡ baseline carrier; uniform limit ≡ fine baseline spline | PROVED + tested | `test_carrier.py`, `test_solver.py` |
| `C^{k−1}` across interfaces; `C^k` fails | PROVED + MEASURED | `test_interface_continuity_is_k_minus_1` |
| order `k+1` in patch, at interface, outside | MEASURED | `test_convergence_order_in_patch_and_at_interface` |
| split = identity (state, density, momentum, energy, loop action) | PROVED up to leak + MEASURED | `test_split_…`, `test_mid_run_split_…` |
| mass ledger exactly 1; one deposit per leaf; no interface density | PROVED + tested | `test_deposit.py` |
| periodic seam / translation equivariance | tested | `test_seam_crossing_and_periodicity` |
| determinism; bitwise checkpoint/restart (both quadrature modes) | tested | `test_solver.py` |
| budget never exceeded; cooldown; reasons; no oscillation | tested | `test_solver.py`, `test_controller.py` |
| locality (would fail for a global regrid) | tested | `test_one_split_is_local`, `test_correction_support_is_local`, `test_adaptive_run_refines_locally…` |
| under-resolved folds: finite, mass exact, topology valid | tested (no accuracy claim) | `test_underresolved_fold_stress…` |
| energy, momentum, action periods | **not invariants**; trends MEASURED | `results/SUMMARY.md` |

## 7. Known failures and limits

1. **Spacing-adaptive quadrature** (midpoint and Gauss panels): worse than the uniform label
   lattice at equal sample count (§3). FAILED.
2. **Late triggers.** A trigger that waits for the first crossing is 30× worse at `t = 3` than one
   that fires by `t ≈ 1.5` (`exp04`); refinement cannot repair error already committed.
3. **Truncation indicators under-estimate accumulated error** (~25×): the tolerance has to be set
   well below the target accuracy (`tol ≈ 1e-6` PM cells for a `3e-8` position error at PM 256²).
4. **Refinement above the PM-cell scale does not converge at the formal order.** After shell
   crossing the sheet has structure at the PM cell scale; errors collapse only once the node
   spacing reaches about one PM cell (quintic `n = ng/2 → ng`: >200× at PM 256²).
5. **On the gate problem the refined region ends up being (almost) the whole sheet**, so local
   refinement brings no accuracy-per-cost advantage over a uniform sheet with `n = ng`
   (`reports/decision.md`).
6. Prototype overheads: regridding rebuilds the operators of the changed levels from scratch
   (`O(tiles)`), and scipy sparse matrices cost 12 bytes per non-zero; peak RSS of locally refined
   arms is 3–5× that of the uniform sheet although the evolving state is smaller.
7. Energy is not conserved by construction (§4); no Lagrangian/symplectic property of the nodal
   map is claimed; no claim about expanding backgrounds, 3D, or novelty.
8. Interpolation at nodes in the skirt of a finer level is exact only to `λ^{h_t B}·surplus`.

## 8. What a 3D extension (Codex Astra) needs

*Carries over verbatim (tensor structure):* octree of `B³` blocks with 26-neighbour 2:1 balance and
the same closed-region node ownership; hierarchical-surplus carrier with tri-quintic B-splines
(P1 with `κ ≤ λ_min^{-3}`: 27 cubic, 422 quintic; P2–P5 unchanged; `C⁴` across faces, edges and
corners with no hanging constraints); three exact line inverses per sweep (exact on boxes; defects
only at re-entrant edges/corners of `T_m`, `≥ h_t B` cells from the data); accumulated
representation with the 3D two-scale mask; evaluation of a block's sample lattice by three small
matrix products; uniform label lattice quadrature with dyadic masses; indicators per direction
(`Δ^{k+1}` along each axis already gives a directional estimate for anisotropic splitting).

*Must be redesigned:*
1. **Storage.** The block-resolution lookup tables are dense here; in 3D they must be hashed, and
   the per-block evaluation windows `(B+k)³` should be assembled from tile ghost layers rather than
   stored as index arrays (2197 indices per block for `B = 8`).
2. **Halo cost.** With `h_t = 1` each patch drags three rings of coefficient tiles: in 3D that is
   `(n+6)³/n³` tiles; a fixed-size skirt in cells (instead of whole tiles) or `B = 16` is needed.
3. **Node budget.** Finding 4 implies nodes at ≈ one PM cell spacing wherever the sheet has folded:
   `~ng³` nodes × 9 state numbers in folded regions, against `(4ng)³` particles × 6 for N-body of
   equal deposit lattice — a 40× state-memory advantage, but only if the folded volume fraction is
   small enough for locality to matter (it was not in the 2D gate problem).
4. **Deposit cost.** PCS in 3D touches 64 cells per sample; the sheet's per-sample advantage
   (no force gather per sample) shrinks or grows with the ratio of gather to deposit cost and has
   to be re-measured in compiled code.
5. **Force law.** The adjoint (energy-derived) nodal force and the consistent mass matrix (Call 2),
   before any claim about energy or action in 3D.
6. **Time.** Comoving KDK factors and Layzer–Irvine; block-wise time steps would need a contract
   for the surplus of nodes advanced on different clocks (the carrier mixes levels linearly).
7. **Parallel ownership.** Blocks as ownership units; coefficient halos of one tile suffice for
   evaluation; the line solves couple ranks along lines — replace by block-local solves with the
   halo as overlap (the one-sweep operator is already local up to `λ^{B}`).
8. **Reference solutions.** A 3D control at `(8ng)³` particles is out of reach; the scoring
   protocol needs sheet self-convergence plus N-body cross-checks instead of a single control.
