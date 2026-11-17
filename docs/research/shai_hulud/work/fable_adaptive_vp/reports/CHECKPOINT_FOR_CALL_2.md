# CHECKPOINT FOR CALL 2 (IX-GRAVITAS)

Compact hand-over. Details: `../FINAL_REPORT.md` (results), `../RESUME.md` (API map, how to run,
what not to repeat), `../research/local_refinement_contract.md` (mathematics, invariants).

## Verified (use as given)

* Composite hierarchical-surplus carrier `F = Σ_m S_m`: linear in the nodal data, SPD per level
  (κ ≤ 56.25 quintic), globally C⁴, interpolatory up to a 2λ⁸·surplus leak, exactly the baseline
  spline in the uniform limit, unchanged by splits. 71 tests.
* Mass ledger exactly 1; uniform label lattice is the right quadrature (adaptive / Gauss fail).
* Nodal dynamics = tracers of the common PM potential (exact-gradient law): action transported
  exactly per node; (E−E0)/K = −3.6e-5 at t = 3 for N-body and resolved sheets alike; momentum
  drift 2e-7 / 2.8e-6; carrier loop-action deviation 9e-7 (n = ng).
* Gate B not passed at t = 3 (37 % folded) on PM 256² and 512²; sheets win at equal wall time and
  memory up to ~20 % folded labels.
* After crossing: n ≈ N_p/2 nodes per dimension are needed, for cubic and quintic alike.

## Failed / do not retry as is

Spacing-adaptive quadrature; Gauss-panel quadrature; refinement gated on first crossing;
truncation-monitor refinement beyond n = ng under a node budget; hard node budgets below the
uniform count; bases coarser than ng with the `ms` trigger (mask empty before the first fold).

## What Call 2 inherits for the force question

The scheme is *not* the Hamiltonian flow of `H_h = ½∫|p_h|²dμ + U_h[ρ_h]`. With the carrier
`x_h(ξ) = Σ_a X_a ℓ_a(ξ)` (ℓ_a = cardinal functions, i.e. columns of the linear carrier operator),

    ∂U_h/∂X_a = Σ_s m_s ℓ_a(ξ_s) ∇V(x_s),        K_h = ½ PᵀG P,  G_ab = Σ_s m_s ℓ_a(ξ_s) ℓ_b(ξ_s),

so the energy-derived force is the **transpose of the evaluation operator** applied to the sample
forces, and the consistent momentum equation carries the Gram matrix `G`. All pieces exist:
`block_values_cf` (evaluation), `lv.P` (two-scale), `lv.Mn` and the line factors (level solves),
`prefilter_periodic` (level 0). The transpose of the one-sweep solve is the same two line solves in
reverse order.

## Exact next experiment

**E1 (first, one day): adjoint force on the uniform quintic sheet** (no refinement; removes every
confound). PM 256², n = 256, 1024² lattice, same IC and control. Arms: (a) current tracer force,
(b) adjoint force with lumped mass, (c) adjoint force with the Gram matrix (CG on G; G is SPD and
banded). Report: x_rms / rho / force error at t = 2, 2.5, 3; (E−E0)/K versus dt (0.01, 0.005);
momentum; loop actions; self-force on an isolated node. *Decisive question:* does the adjoint
force change finding "n ≈ N_p/2 after crossing" (5.74e-8 at PM 512² n = ng vs 5.16e-9 N-body)?
If yes, re-run `exp03_gate.py`; if no, the equal-cost verdict stands and the force work is about
invariants only.

**E2 (parallel, half a day): hybrid release** on the `ms` mask (see `RESUME.md` §6.3) — the
route that can actually deliver N-body accuracy at t = 3 with sheet economy in the single-stream
region.

Start commands: `RESUME.md` §2. Controls must be regenerated once (`exp01_reference.py 256 4096`,
≈ 20 min) because they are not in git.
