# PROJECT DUNE — SHAI-HULUD: Master Research Dossier v1.2

**Date:** 24 November 2026  
**Project:** CHUÍ Cosmological Simulator / Experimental Numerical Research Programme  
**Status:** archived research baseline; development deferred behind flagship physics.  
**Scope:** scientific architecture, numerical evidence, mathematical hypotheses, and a reproducible continuation plan.  
**Authority:** supersedes v1.1 as a roadmap, not as a replacement for the original evidence. This document is not peer reviewed and does not authorize changes to production CHUÍ.

## 0. Read this first

The Maker is an attempt to construct a mathematically controlled method in which the *representation of matter* can change during evolution. It is not an administrative scheduler, a collection of exotic algorithms, or a claim to have discovered new physical laws. CHUÍ's flagship remains conventional N-body TreePM, cosmological time integration, Godunov hydrodynamics and patch AMR, with its baryonic and magnetic physics developed through their own acceptance campaigns.

The first two SHAI-HULUD campaigns produced action-aware phase-space reconstruction, a conditional locally Lagrangian atlas, and a 2D self-gravitating Vlasov–Poisson prototype. The third campaign, The Adaptive Maker (Call 1, 17 November), implemented actual local block h-refinement and subjected the sheet to equal-wall-time and equal-memory comparisons. The outcome is unequivocal: **the original late-time all-sheet replacement gate failed**. The uniform high-order sheet is economical before severe folding, but after shell crossing its node and quadrature requirements approach those of particles. Locality disappears on the benchmark, and the prototype's regrid memory overhead is substantial.

Consequently:

- N-body remains the production collisionless representation; TreePM remains the production gravitational field method. Sheets and TreePM are not competing concepts: one represents matter, the other approximates gravity.
- Preserve the early-time sheet as a specialized experimental representation and geometric analysis instrument. Do not integrate the all-sheet prototype as a production replacement.
- The next DUNE experiment, when work resumes, is a **2D hybrid sheet-to-particle release** under shared gravity, with causal predictive triggers, exact label-mass ownership and a fair N-body baseline. This transition has not yet been implemented.
- Remove repeated optimal-transport node redistribution from the collisionless critical path. Preserve its negative evidence and the separate IX-FORGE gas-mesh research branch.
- Keep expensive exact atlases, generating families, cohomology, variational gravity and GENERIC as independent scientific research branches. Do not require all of them for the first hybrid experiment.
- The ultimate mathematical objective is a representation-independent, measure-valued Hamiltonian evolution with certified approximate transfer maps. The proposed unifying formulation is not a new physical equation, and no general hybrid error theorem or Lean formalization has been completed.
- **Development priority:** finish CHUÍ's flagship modules, particularly MHD, and their source/validation campaigns. DUNE is on the back burner, not abandoned. No speculative DUNE code is to enter production as a side effect of this dossier.

This file is intentionally readable independently, but v1.1 and the original research packages remain the detailed historical record.

## 1. Revision record and evidence policy

| Version | Interpretation |
|---|---|
| v0.x | Dune-inspired experimental family; Hamiltonian, phase-sheet, OT and GENERIC synthesis |
| v1.0, 16 Nov | First Fable manuscript and DeepSeek falsification audit; prescribed-field prototypes |
| v1.1, 16 Nov | Self-gravitating 2D breakthrough, conditional atlas, negative OT/atlas performance; local refinement and IX split |
| v1.2, 24 Nov | Call 1 closed; Gate B failed, early-sheet niche retained; hybrid release next; wider representation research, formal mathematics and deferred development |

Evidence vocabulary:

- **Established**: mathematical fact under stated assumptions, rather than a new SHAI-HULUD invention.
- **Author-derived**: a derivation in the research package, not an independently peer-reviewed or machine-checked theorem.
- **Recorded**: reported result in the archived source, log, JSON or report.
- **Previously reproduced**: explicitly identified reruns from earlier review sessions, not a new test in this update.
- **Proposed**: a design or candidate algorithm not yet demonstrated.
- **Open / negative**: unresolved problem or failed tested approach.

This update inspected the repository's Call 1 final report, resume instructions, README, refinement contract and gate-evaluation data, plus the v1.1 dossier and later exported discussion. It does **not** claim a fresh execution of the 71-test suite or the full benchmark campaign. The repository's current research tree is the primary source; the conversation is evidence of design decisions, not a scientific proof. A referenced pytest.log was not retrievable at its expected path during this update; 71/71 is therefore a reported count corroborated by the report and README, not a newly inspected log.

## 2. Scientific charter and complete DUNE family

The goal is to adapt the representation to the physics, with an explicit physical and numerical contract for every approximation or information-losing transformation. The thematic hierarchy reflects scientific ambition, not a rigid software dependency tree.

| Module | Intended mathematical role | Status and boundary |
|---|---|---|
| SHAI-HULUD — The Maker | Actual numerical method for adaptive Hamiltonian/geometric and eventually metriplectic evolution of heterogeneous representations | 2D sheet research exists; hybrid and generalized method open |
| SPICE — Structure-Preserving Interchange and Conservative Evolution | Reusable operators, discrete gradients, conservation/entropy contracts, error ledgers, differentiable interfaces | Proposed shared foundation; may serve ordinary CHUÍ independently. Name also honors electrical-engineering SPICE |
| ARRAKIS | Periodic domain geometry, topology, geometric queries, coordinate and ownership conventions | Research extension of existing infrastructure; not a second domain authority |
| LITTLE MAKER | Local h/p/hp adaptation, geometric refinement, subdivision, hierarchical compression | Call 1 implemented local h-refinement in research Python; p and anisotropy remain proposed |
| IX-GRAVITAS | Variational gravitational field discretization, adjoint deposition/gather, compatible operators, mass matrices and quadrature | Priority mathematical field branch; not yet implemented in Call 1 |
| IX-FORGE | Laguerre/power diagrams, semi-discrete OT, gas-side moving mesh and ALE geometry | Optional future research; no collisionless dependency |
| KWISATZ HADERACH | Controlled representation transfers, conservative projection, information-loss estimates | Hybrid transfer is its first decisive future application |
| MENTAT | Error/cost modeling and resource-constrained selection | Administrative policy, explicitly not The Maker's mathematical identity |
| GUILD NAVIGATOR | Predictive error growth, pre-caustic seeding, short-horizon adaptation | Call 1 tested state-only predictive triggers; general controller open |
| WATER OF LIFE | GENERIC/metriplectic thermodynamics, irreversible physical closures | Oscillator-level research only; open in cosmological gas |
| BENE GESSERIT | Discrete adjoints, sensitivity and differentiable evolution | Future; topology switches complicate differentiation |
| GOLDEN PATH | Alternative physical histories, uncertainty quantification, ensemble design | Future; not a current solver feature |
| GOM JABBAR | Adversarial tests, proof obligations, pathological geometry, Lean formalization | Mathematical audit role; formal certification remains open |

IX is one family with two research divisions. Thus this remains the original twelve named modules, with IX subdivided. The Little Makers reshape the computational environment; the Maker develops the physical numerical formulation; SPICE makes the reusable mathematical machinery available. A successful scientific result may be promoted into ordinary CHUÍ only after its own acceptance process.

## 3. Physical and geometric foundations

### 3.1 Collisionless evolution in an expanding periodic domain

Use the comoving configuration torus $X=\mathbb T^d$, typically $d=3$, and canonical phase space $T^*X$. For clarity use **specific canonical momentum** $p=a^2\dot x$ and a mass-weighted phase-space measure $\mu$ (the individual particle canonical momentum is $m p$). Let $\Phi$ be the peculiar Newtonian potential, $\rho_c$ the comoving mass density, and $a(t)$ the scale factor. In the usual Newtonian cosmological convention,

$$
H_\mu(x,p,t)=\frac{|p|^2}{2a^2(t)}+\Phi_\mu(x,t),\qquad
\dot x=\frac{p}{a^2},\quad \dot p=-\nabla_x\Phi_\mu .
$$

$$
\Delta_x\Phi_\mu=4\pi G a^{-1}(\rho_c-\bar\rho_c),\qquad
\rho_c=(\pi_x)_\#\mu.
$$

Here the projected measure is interpreted as a density where appropriate; in numerical work it is deposited or regularized. With physical density $\rho_{\rm phys}=a^{-3}\rho_c$, the Poisson factor is equivalently $4\pi G a^2(\rho_{\rm phys}-\bar\rho_{\rm phys})$. These conventions must be translated to the actual CHUÍ units before implementation. For any sufficiently smooth time-independent test function $\varphi$,

$$
\frac{d}{dt}\langle\mu_t,\varphi\rangle
=\langle\mu_t,\{\varphi,H_{\mu_t}\}\rangle,\qquad
\{f,g\}=\nabla_x f\cdot\nabla_p g-\nabla_p f\cdot\nabla_x g .
$$

This is the weak Vlasov evolution, not a new Maker equation. For time-dependent test functions include $\langle\mu,\partial_t\varphi\rangle$. Assume adequate well-posedness/force regularity; unregularized singular gravitating measures require separate analysis.

The total peculiar energy is $E=K+U$, where $K=\int |p|^2/(2a^2)\,d\mu$ and $U=\tfrac12\int\Phi\,d\rho_c$ with zero-mean potential. For the standard collisionless cosmological system with suitable regularity and no external exchanges, the Layzer–Irvine balance is $\dot E+(\dot a/a)(2K+U)=0$. Do not impose ordinary energy constancy in an expanding cube. The static 2D test normalizes box and mass to one, sets $4\pi G=1$, uses $p=\dot x$, and solves $\Delta\phi=\rho-1$.

### 3.2 Cold sheets and caustics

For initially cold matter let $Q=\mathbb T^d$, $\mu_0(dq)=\rho_0(q)dq$, and $F_t(q)=(x(q,t),p(q,t))$. Then $\mu_t=(F_t)_\#\mu_0$. Where the position projection has finitely many regular preimages, the spatial density is

$$
\rho_c(x,t)=\sum_{q:x(q,t)=x}\frac{\rho_0(q)}{|\det D_qx(q,t)|}.
$$

Shell crossing is a degeneracy of the projection $D_qx$, not necessarily a singularity of $F$ in full phase space. Negative determinant is a useful orientation/folding diagnostic, not a count of multistream mass or a universal release threshold. Multiple streams can coexist in a spatial cell even if particular labels retain positive orientation.

Use $\theta=p\cdot dx$ and $\omega=-d\theta=\sum dx_i\wedge dp_i$. Under smooth Hamiltonian evolution, $\omega$ is preserved even when the Hamiltonian explicitly depends on time. If $F_0^*\omega=0$ and the embedding remains regular, $F_t^*\omega=0$. For a periodic sheet, $F_t^*\theta$ can be closed without being globally exact. In unit-period label coordinates write

$$
F_t^*\theta=dS_t+\sum_{i=1}^d c_i\,dq_i.
$$

The $c_i$ encode periods around fundamental cycles and are invariant under exact Hamiltonian isotopies with the stated conventions. A one-dimensional sheet has a vacuous pulled-back two-form condition; it can still have nontrivial action periods. Do not infer multi-dimensional Lagrangian preservation from a 1D test. Exact continuous reparameterization $T:Q\to Q$ with the correctly pulled-back mass measure leaves the physical pushforward unchanged; finite-order resampling generally does not.

### 3.3 Retained geometric research

Fable's initial action-aware reconstruction evolves nodal $(x,p,S)$ and uses action-compatible reconstruction instead of unconstrained interpolation. Its later atlas constructs conditional cell-wise Lagrangian carriers through configuration-space folds in 2D. Smooth, compatible local charts and action data matter; fallback regions and chart degeneracies prevent a general global theorem. Generating families offer a possible extension:

$$
\partial_\alpha S(x,\alpha)=0,\qquad p=\partial_x S(x,\alpha).
$$

Multiple critical points may encode multiple streams. Local existence is not a guarantee of one global generating family, cheap evaluation or smooth chart gluing. Retain action periods, discrete differential forms/cohomology, generating families, compatible atlas transitions, multisymplectic integration and structure-aware multiresolution compression as research tools, not requirements in every cell. The 3D/6D conforming carrier and face/corner compatibility remain open.

## 4. Historical results through v1.1

The v1.0 and v1.1 continuity archives preserve the original papers, code, ledgers, data and negative findings; consult them for exact definitions and all tables.

First Fable mission: 31-page manuscript; 1D action-corrected interpolation maintained the tested action-period residual near roundoff and improved normal reconstruction error substantially in collapse and chaotic prescribed-field tests. This was not a corresponding improvement in full position/momentum or gravitational error. In smooth 2D tests the structure-aware element reduced circulation defects at a cost to ordinary momentum accuracy. The initial monitor derivation contains an apparent exponent inconsistency: minimizing $\int C h^{3-d}dq$ under $\int h^{-d}dq=N$ gives a stationary scaling $h^3\propto C^{-1}$ (with dimension-dependent qualifications), not the claimed $C^{d/(d+3)}$ sampling density. The latter follows for $\int C h^3dq$ under the stated count constraint. At $d=3$ the first objective is independent of $h$; no optimum can be inferred from it. The gravity branch no longer depends on this monitor, but any reuse requires re-derivation.

DeepSeek: a restricted linear averaging no-go theorem, not a no-go for all adaptive sheet methods. For nonnegative row-stochastic $D$ acting identically on position and specific momentum, positive diagonal mass matrices $M,M'$, and a weighted canonical two-form, preservation requires $D^\top M'D=M$. Off-diagonal entries are sums of nonnegative terms; each row can have at most one nonzero, so only copying/permutation with appropriate split masses can satisfy it. Genuine averaging/coarsening in this class is not symplectic. An exact kinetic-energy loss identity follows from weighted variance under the averaging operation; removed resolved energy is not automatically physical heat. DeepSeek's graded-mesh test at 128² gave RMS errors 0.004723 uniform, 0.001754 freshly sampled graded, and 0.005625 remapped graded: resampling erased the geometric benefit. Its 2D two-form test separated geometric defects from caustics; its GENERIC oscillator is a useful reduced verification problem.

Second Fable (breakthrough): 2D periodic static self-gravity with sheet quadrature and KDK, with conditional Lagrangian atlas through caustics. Earlier review ran 30/30 regression tests. For a representative t=3 comparison with 128² deposited samples, 4,096 cubic sheet nodes gave position RMS 2.53e-6 versus 3.68e-6 for 16,384 N-body particles, at comparable time. At 1,024 nodes the cubic sheet (1.65e-5) was more accurate than the action-aware (5.17e-5) and exact-atlas (6.01e-5) alternatives in that experiment. Single 64² remap timings were approximately 3.0 ms cubic, 16.6 ms quintic, 888.7 ms exact cubic atlas, 1,277.5 ms exact quintic atlas. Python timing does not predict optimized C++ costs. The exact atlas and repeated OT redistribution did not justify default use. Keep ordinary high-order splines, nodal action diagnostics and separate geometric verification.

These findings motivated the local-refinement campaign now recorded below. Neither the initial prescribed-potential cosmological example nor any static 2D test proves a full expanding 3D cosmological solver.

## 5. Call 1 — The Adaptive Maker: implementation and scientific evidence

Primary source: [FINAL_REPORT.md](work/fable_adaptive_vp/FINAL_REPORT.md), dated 17 November 2026; also [local refinement contract](work/fable_adaptive_vp/research/local_refinement_contract.md), [README](work/fable_adaptive_vp/README.md), [RESUME](work/fable_adaptive_vp/RESUME.md) and [gate evaluation](work/fable_adaptive_vp/results/gate_evaluation.json). The code is already under this repository's research directory. No new solver is delivered by this dossier.

### 5.1 Frozen reference and numerical contract

- Label torus [0,1)², periodic displacement $x(q)=q+d(q)$, action with global periods, total mass one, static background, $4\pi G=1$.
- Three-mode GrowingModeIC, amplitudes (0.10, 0.07, 0.03), phase 0.9; KDK, dt=0.01, 300 steps to t=3.
- PM with PCS assignment, continuous Green function $-1/k^2$, exact-gradient gathered potential. One field solve per step, merged kicks.
- Python/NumPy/SciPy/Numba; reported 71 regression tests. Baseline reproduction includes 30/30 previous tests and bit-identical control and four saved arms.
- N-body control uses the same field conventions and compiled kernels, in-place kick/drift and one gather per step. Its earlier slower version was rejected as unfair. Timed phases exclude scoring controls from memory; score offline.
- Same-PM high-resolution controls: 4096² particles. Field discretization error is materially larger than the inter-representation errors.

### 5.2 Topology, composite carrier and transfer

A periodic quadtree of B×B-cell blocks (B=8) has exact dyadic tiling, 2:1 balance including corners and periodic seams, deterministic node ownership and logged parent/child lineage. All base blocks exist. Every label mass belongs to exactly one leaf. Hanging nodes are ordinary owned degrees of freedom, not independent overlapping mass elements.

For cubic or quintic centered cardinal B-splines, define the tensor product basis $\beta_{m,j}$ on grid level $m$, and a fine correction support $T_m$ enlarged by a halo. The represented carrier is

$$
F_h=\sum_m S_m,\qquad
S_m=\sum_{j\in T_m}c_{m,j}\beta_{m,j}.
$$

The base is the periodic interpolant. Each fine correction interpolates the nodal surplus relative to coarser levels at active new nodes and is zero at other knots inside its support. The author-derived 2D principal-submatrix SPD argument bounds condition numbers by 9 for cubics and 56.25 for quintics. The sum is globally $C^{k-1}$ ($C^4$ for quintic), including coarse/fine interfaces, without hanging-node trace constraints. The exact full refinement limit equals the baseline fine spline. The prototype uses a fixed line-solve sweep; exactness on general re-entrant geometries is not claimed. The finite-support skirt causes a small interpolation leak. A measured rate $2\lambda^8$ of the surplus with $\lambda=0.431$ (quintic) or 0.268 (cubic) must not be conflated with a universal sharp bound; the contract distinguishes the weaker theoretical rate.

Splitting reads the new nodes from the old represented carrier, never exact reference data. It is near-identity up to the skirt/solver defect; the report records tiny density, momentum, energy and action differences, exact mass ledger. Rebase of a complete level is exact up to floating arithmetic. Coarsening drops surplus and is not identity; its error must be accounted for. Measured smooth-problem convergence was near fourth/sixth order, including interfaces, not a proof of that rate at caustics or in coupled dynamics.

Mass belongs to label cells. In default lattice mode the quadrature samples form a single uniform label lattice across refinement levels (when the base subdivision supports the highest level). The ledger is exactly 1.0 for the tested dyadic weights; PM grid deposition differs at roundoff. **Adapt the carrier, not the quadrature lattice** in the current implementation. Adaptive spacing-driven block quadrature was 30–100× worse at equal sample count in the tested regime, and Gauss panels 3–6× worse; boundary cancellation and aliasing explain the observed sensitivity. These are experiment-specific conclusions, not universal quadrature impossibility theorems.

### 5.3 Measured gate and verdict

The original Gate B required sheet position AND density errors no greater than N-body at equal wall time, and error no greater than half N-body at equal memory, at t=3 with at least 30% folded labels on PM 256² AND 512². **NOT PASSED.** No all-sheet replacement in the multistream regime is justified.

Selected reported comparisons (errors are relative to controls for the same PM field model):

| Regime / method | Evolving nodes/particles | Wall time | Position RMS | Interpretation |
|---|---:|---:|---:|---|
| t=2, PM 512², quintic n256, 2048² lattice | 65,536 | 57 s | 2.02e-9 | early-collapse sheet |
| t=2, PM 512², N-body 2048² | 4,194,304 | 97 s | 1.92e-9 | approximately matching position error; control floor matters |
| t=3, PM 256², uniform quintic n256 | 65,536 | 15.8 s | 2.57e-8 | wall position ratio 0.73, density 1.14; fails joint gate |
| t=3, PM 256², local 128→256 pred | 65,536 | 16.5 s | 2.85e-8 | wall position 0.89, density 1.13; peak RSS overhead |
| t=3, PM 512², uniform quintic n512, 2048² lattice | 262,144 | 70.6 s | 5.74e-8 | wall position 5.02, density 3.50; fails |
| t=3, PM 512², local 256→512 pred | 254,992 | 79.8 s | 5.74e-8 | wall position 6.84, density 3.45; fails |
| t=3, PM 512², N-body 2048² | 4,194,304 | 96.8 s | 5.16e-9 | reference ladder point |
| t=3, PM 512², uniform n1024 | 1,048,576 | 118.7 s | 6.20e-9 | comparable accuracy at greater wall cost |
| t=3, PM 512², local 512→1024 ms | 876,680 | 137.1 s | 6.18e-9 | similar accuracy, fewer node-steps but no speed advantage |

The last two timings are indicative because those arms ran concurrently; do not present them as clean matched timings. In the report's early t=2 comparison (5.8% folded labels), the sheet uses 64× fewer evolving nodes, approximately 1.7× less wall time and 17× lower peak incremental memory (24 vs 416 MB) than the 2048² N-body run at comparable position error. The equal-wall-time error ratios 0.28 position and 0.42 density are separately reported. At t=2.5 (19.5% folded labels), PM 512² uniform n512 has 0.94 position and 0.71 density ratios. These are **conditional early-regime advantages**, not speedups over production CHUÍ or physical error improvements. The reported position control resolution at PM 512² is about 5e-9, so very small early errors should not be over-interpreted.

By t=3 the folded-label fraction is approximately 37%. Half the blocks contain a fold; 78–86% lie in or adjacent to multistream regions; refined regions cover 85–100% of the sheet for relevant tolerances. Local refinement saves at most roughly 16% of terminal nodes in the deep comparison. The ms run has 117M versus 315M node-steps, but deposition samples dominate cost. Local arms fail measured RSS comparisons despite good evolving-state byte ratios because sparse operator storage and full regrid rebuilds are expensive. Uniform sheets retain an equal-memory advantage in tested arms; that is not equivalent to a combined memory/time/accuracy win.

The PM 256²→512² field change produces position differences around 5.7e-4 at t=3, far above the same-PM inter-representation differences. The tests therefore compare discretizations of a common numerical model, not convergence to physical Vlasov–Poisson at those error values. Density sampling ripple is nonmonotone; report treats near-unity density ratios within roughly ±20% as ties.

### 5.4 Trigger and force lessons

- First-crossing-triggered refinement is too late: uniform refinements at t=1.0/1.5/1.75/2.0 gave 2.58e-8 / 2.70e-8 / 5.12e-8 / 1.22e-6. A newly refined carrier cannot repair error committed during earlier evolution.
- Predictive monitor evaluates an approximate future state $x+\tau p+\tfrac12\tau^2 a$ in static coordinates. It fires around 1.4–1.6 in this case and beats the pure d6 truncation indicator at equal tolerance.
- Once h approaches the PM scale, truncation monitors lose correlation with actual error. Near-fold boundaries can be limiting even when the sheet looks smooth.
- A hard node budget can fail ungracefully: at 0.6 ng² nodes the reported error is 3.5e-6, about 120× the uncapped case. A budget must support a deliberate fallback, not silent degradation.
- The Eulerian multistream proximity mask reproduces the uniform fine solution with fewer nodes, but requires a pre-caustic seed because the mask is empty before the first fold. Its lead dependence is nonmonotonic (PM 512²: 6.2e-9, 7.9e-9, 1.5e-8 for lead 0.3, 0.6, 1.0). Test the frozen-momentum-error hypothesis directly; do not assume a larger lead is better.
- The nodal common-potential force is not the derivative of the quadrature-defined sheet energy with respect to its DOFs. The measured relative energy drift at t=3 is about -3.6e-5 for resolved arms, with momentum and loop-action defects. No semidiscrete self-consistent symplectic or exact energy theorem is established.

Call 1 is one large collapse in a static 2D box, one-thread Python/Numba, not a 3D or expanding cosmological study. The observed approximately N_p/2 nodes per dimension after folding is empirical, not a universal law.

## 6. Revised Maker: heterogeneous measure representation

### 6.1 One physical law, multiple numerical realizations

Use a partition of original label mass into active sheet material and material released into other representations. Let $\mathcal A_t$ denote active labels (possibly disconnected), and $\mu_0|_{\mathcal A_t}$ its restricted measure. A candidate representation is

$$
\mu_t^{\rm SH}=(F_t)_\#(\mu_0|_{\mathcal A_t})+
\sum_i w_i\,\mathcal K_i(t),
\qquad w_i\ge0,\quad \mathcal K_i\ge0,\quad \int d\mathcal K_i=1.
$$

Dirac $\mathcal K_i=\delta_{z_i}$ gives ordinary particles. Extended nonnegative kernels can encode coarse phase-space moments. The restriction and release weights must be disjoint in *mass ownership*, even where the images overlap in position or phase space. All components contribute to one gravitational density and one synchronized field. Numerical carriers are not different physical species; cold collisionless matter does not acquire physical pressure, collisions or entropy merely because it is represented by a packet.

The proposed evolution law is the weak Vlasov equation in §3 plus a **discrete adaptive approximation rule**. At an event time $\tau_j$ a transfer operator replaces a representation $\mu^-$ by $\mu^+=\mathcal P_j(\mu^-)$. Exact equality of measures is exceptional. The operative question is what observables and invariants remain controlled, and at what cost.

### 6.2 The immediate hybrid release hypothesis

Use the existing Call 1 leaf quadrature, state-only triggers and lineage. Begin with the simplest one-way deterministic sheet-to-Dirac release; no Gaussian or wave method is needed to test its value.

Proposed event transaction:

1. Synchronize all momenta to the same KDK event time (Call 1 merges kicks; reading unsynchronized values would produce a false transfer error).
2. Predict under-resolution from pre-caustic and multistream proximity diagnostics, local force/deposition error estimates and the cost/memory budget. Online triggers must never read scoring/reference arrays.
3. Select complete, uniquely owned label leaves. If a candidate would leave an invalid topology or violate the budget, reject or defer it with a recorded reason.
4. Evaluate the old carrier at the *same* uniform label quadrature points already used for PM deposition. Assign each particle its existing positive sample weight and synchronized position and momentum.
5. Construct a provisional new state and verify mass ownership, finite values, positivity, momentum/kinetic moments, deposition residual and field compatibility. Remove the released mass from sheet deposition exactly once. Keep any needed geometric ghost/lineage information explicitly non-mass-bearing.
6. Atomically commit the transfer, update shared density/force and ownership ledger, and record pre/post diagnostics. If the implementation cannot satisfy the contract, roll back.
7. Resume common gravitational evolution. Initially, no particle-to-sheet reconstruction, no automatic merging of distinct streams and no claim of reversibility.

If quadrature samples and weights are identical on both sides, the discrete deposited density is identical at the event up to summation/roundoff, and its quadrature-defined linear momentum and kinetic energy are unchanged. This elementary conditional result is not equality with the continuous curved sheet or a guarantee of future trajectory equality: particle forces are now gathered at samples, rather than interpolated from evolving sheet nodes. For arbitrary alternative sampling, potential-energy and force defects must be assessed anew.

Define the released mass fraction, active node count, particle count, field sample count, node/particle steps, transfer frequency and discarded geometry. Report all of them; counting only sheet nodes after release hides costs. The union of active sheet mass and released particle mass must equal initial total mass at every checkpoint, including across periodic seams, domain migration and restart. The topology of retained sheet and the status of formerly shared nodes must be resolved before calling this production-capable.

### 6.3 MENTAT / Navigator contract

Select admissible actions from retain, h-refine, eventually p-refine, release, or (only after validation) coarsen. Use hard physical and resource constraints first, then minimize a normalized objective such as

$$
J(r,K)=\eta_{\rm obs}(r,K)^2+\lambda_g\eta_{\rm geom}^2+
\lambda_t C_{\rm time}+\lambda_m C_{\rm memory}+\lambda_c C_{\rm comm}+
\lambda_P\eta_{\rm transfer}^2 .
$$

This is a policy hypothesis, not a universal variational law. Predictive trigger must act before committing coarse error; thresholds are calibrated on multiple held-out ICs, not inferred universally from the 20% folded-label observation. Use hysteresis/cooldown; release is one-way initially. If constrained choices cannot meet accuracy, report failure rather than quietly violate invariants. Assess cost with actual field-sampling work, not node counts alone.

## 7. Beyond N-body release: the research candidates

Particles are the reference release representation, not the final intellectual limit. Keep each alternative separately falsifiable. Candidate order below is a proposed research priority, not demonstrated performance.

### 7.1 Moment-enriched carriers and Gaussian clouds

A packet has mass $w$, center $m=\langle z\rangle$, covariance $\Sigma=\langle(z-m)(z-m)^\top\rangle\succeq0$ in $z=(x,p)$, potentially additional moments. Compute these from a positive quadrature on the released label region, lifting periodic positions consistently. A Gaussian or other nonnegative normalized kernel can match these coarse moments; for a singular covariance use a degenerate measure or an explicitly reported regularization, not a silently invertible Gaussian.

For a smooth Hamiltonian vector field $b(z)=J\nabla H$, the exact moment equations for one normalized transported packet are

$$
\dot m=\langle b(z)\rangle,\qquad
\dot\Sigma=\langle (b(z)-\dot m)(z-m)^\top+(z-m)(b(z)-\dot m)^\top\rangle.
$$

Under a local affine approximation, $A=Db(m)$,

$$
\dot m\approx b(m),\qquad \dot\Sigma\approx A\Sigma+\Sigma A^\top.
$$

The covariance equation is exact for affine dynamics, not for general self-gravity; the neglected hierarchy contains higher moments. Under affine Hamiltonian evolution $\Sigma(t)=S\Sigma_0S^\top$ with $S$ symplectic, preserving positive semidefiniteness and corresponding symplectic invariants. Nonlinear closure and packet split/merge need their own consistency tests. Spatial marginals determine density deposition; internal velocity dispersion contributes to kinetic energy. A packet is not a collisional fluid element.

A bound for replacing a packet by its center in a smooth observable with Hessian norm bounded by $B$ is $|\langle\varphi\rangle-\varphi(m)|\le (B/2)\operatorname{tr}\Sigma$, by Taylor expansion in a suitable local lift. Matching covariance cancels the quadratic moment difference between distributions; higher derivatives and higher moments then control the residual. These are restricted approximation results, not a universal closure theorem.

Compare center particles, moment-enriched particles and anisotropic clouds at matched accuracy and *total* bytes (a full 6×6 symmetric covariance needs 21 additional scalars per packet before metadata). Test harmonic transport, tidal shear, caustic crossings, splitting, covariance conditioning, density deposition, and force accuracy. Strong phase mixing may force so many packets that the method loses again. If so, archive the negative result. Existing cloudy Vlasov methods are prior art; novelty cannot rest on using Gaussian kernels alone.

### 7.2 Wave representation (Schrödinger method)

In the static specific-momentum convention, one possible normalization is

$$
i\hbar_{\rm eff}\partial_t\psi=-\frac{\hbar_{\rm eff}^2}{2}\Delta\psi+\Phi\psi,
\qquad \rho=|\psi|^2 .
$$

With comoving coordinates the kinetic term acquires $a^{-2}$ under corresponding conventions. Effective $\hbar$ sets approximation/regularization scales; the Wigner/Husimi interpretation and interference, phase coherence, positivity, initialization and matching at hybrid boundaries require explicit treatment. One complex wavefunction is not an exact arbitrary multistream classical distribution. Research this as an independent controlled comparison, not as a drop-in lossless sheet release. No claim about physical ultralight dark matter is implied.

### 7.3 Other admissible approximations

Waterbags evolve regions of approximately constant phase-space density; cold singular measures cannot be exactly replaced by finite-thickness bags. Adaptive Vlasov elements or spectral bases can retain continuum information but encounter six-dimensional memory costs. Multiresolution coefficients and generating-family compression may reduce sheet costs in smooth regions; independent position/momentum truncation can break geometric structure. Preserve all as optional laboratory tracks. No general switch among them is implemented.

## 8. IX-GRAVITAS: the mathematical field programme

This is a separate, high-value research branch that may change the observed all-sheet performance verdict, but must not be assumed to do so. Call 1's actual carrier is linear in nodal data; use that fact to derive the discrete force consistently.

For fixed topology, let $x_h(q)=\sum_a\ell_a(q)X_a$ in an appropriate lifted/periodic coordinate representation; a known affine baseline may be treated separately. At quadrature sample $s$ with mass $m_s$, $x_s=x_h(q_s)$, and similarly $v_s=\sum_a\ell_a(q_s)V_a$. Let $W_g(x)$ be the PM assignment kernel and $\rho_g$ the resulting grid mass/density in a consistent grid inner product. If the mean-subtracted Poisson operator is self-adjoint on the zero-mean space,

$$
U_h(X)=\frac12\langle\rho_h-\bar\rho_h,G_h(\rho_h-\bar\rho_h)\rangle_h.
$$

With $\Phi_h=G_h(\rho_h-\bar\rho_h)$, differentiating this quadratic functional gives

$$
\partial_{X_a}U_h=\sum_s m_s\,\ell_a(q_s)\,\nabla V_h(x_s),
$$

where $V_h$ is the potential interpolated using the adjoint assignment convention. The consistent mass matrix is

$$
M_{ab}=\sum_s m_s\ell_a(q_s)\ell_b(q_s),\qquad
P=MV,\qquad
K_h=\frac12 P^\top M^{-1}P
$$

for static geometry and $M$ positive definite (check independence, quadrature rank, constraints, and vector-component tensor factors). For fixed topology, a time-independent semidiscrete Hamiltonian $K_h+U_h$ yields $\dot X=M^{-1}P$, $\dot P=-\nabla_XU_h$. In comoving coordinates $K_h$ gains the appropriate $a^{-2}$ factor when $P$ is canonical. Do not mix nodal specific momentum, mass-weighted canonical generalized momentum and velocity. Changing topology changes basis and mass matrix and is not covered by fixed-topology conservation.

Candidate Call 2: implement the transpose of the actual linear carrier and quadrature evaluation, compare adjoint force to finite differences of $U_h$, check matrix conditioning and mass lumping error, then test field smoothness, self-force, translations, momentum balance, energy/timestep convergence and matched-cost hybrid benchmarks. Nodal tracer KDK and variational nodal KDK are different discretizations. A Hamiltonian derived from a grid energy does not automatically have continuous translation symmetry; zero self-force and exact total momentum need independent proof. Fourier-Galerkin, compatible/mimetic operators and high-order geometric quadrature are alternatives. Consider short-range TreePM compatibility only after PM validation. No competing production gravity authority.

## 9. The ultimate mathematical question and formal verification

The original ambition of a single polynomial differential equation is retained as a question, not a requirement. An algebraic differential equation is polynomial in unknowns/derivatives; differential algebra studies such structures. Vlasov–Poisson is an infinite-dimensional, nonlinear, nonlocal evolution, and changing representation introduces event maps and dimension changes. An ordinary finite polynomial ODE has not been shown to encode this adaptive system. Do not replace physical equivalence with a clever but unverified encoding. A measure-valued Hamiltonian evolution plus a hybrid event law and constrained approximation principle is currently the more natural candidate. Generating a unified differential/algebraic formulation remains a legitimate long-horizon mathematical experiment.

### 9.1 A precise proposed theorem target

Consider a regularized self-consistent model on an appropriate phase-space metric, with common total mass normalized to one and finite moments. Let $\mathcal E_{t,s}$ be an exact evolution on measures and let $d$ be a metric controlling desired observables. Assume the **explicit stability hypothesis**

$$
d(\mathcal E_{t,s}\mu,\mathcal E_{t,s}\nu)\le e^{L(t-s)}d(\mu,\nu).
$$

At transfer time $\tau_j$, assume mass/positivity preservation and $d(\mathcal P_j\mu,\mu)\le\epsilon_j$ on the admissible domain. Let numerical evolution steps have independently bounded local defects $\delta_k$. Then repeated use of the triangle inequality and stability gives a conditional global bound of the form

$$
d(\mu_T^{\rm SH},\mu_T)\le e^{LT}d(\mu_0^{\rm SH},\mu_0)+
\sum_j e^{L(T-\tau_j)}\epsilon_j+
\sum_k e^{L(T-t_{k+1})}\delta_k .
$$

Proof sketch: at a transfer, error increases by at most its certified defect; on an interval exact evolution multiplies prior error by at most the stability factor, and the numerical step adds its local defect. Iterate and telescope. The statement is elementary **conditional mathematics**, not a proof that SHAI-HULUD satisfies the assumptions. Uniform stability for singular Newtonian forces, certified event error bounds, variable quadrature, chaotic growth and long-time usefulness remain open. Start with smooth softened forces or prescribed Hamiltonians and an appropriate bounded-Lipschitz/Wasserstein metric; using a Wasserstein distance as an *error measure* does not reintroduce OT node redistribution into the solver.

A useful first restricted result may be mass ownership plus event-wise moment/observable error bounds for a fixed positive quadrature and regularized field. The eventual novelty would be a specific validated adaptive method with usable error estimates, not the recurrence or weak Vlasov law by itself.

### 9.2 GOM JABBAR / Lean research plan

Use Lean 4/Mathlib (or a later suitable version) once the definitions are settled. No Lean file was built, compiled or certified by this update.

| Proposed theorem | Initial formal scope |
|---|---|
| Nonnegative averaging no-go | Finite matrices, positive diagonal weights, row-stochastic map; explicit dimension and copying characterization |
| Mass ownership | Disjoint finite partitions, positive weights; split/release preserves total mass |
| Moment preservation | Same weighted sample list before/after a representation event; separate continuous quadrature error |
| Gaussian moment bound | Finite-dimensional affine transport / covariance PSD; Taylor remainder with clear domain |
| Event defect recurrence | Abstract metric space, Lipschitz evolution, triangle inequality and finite event schedule |
| Fixed-topology variational field | Discrete self-adjoint energy derivative, consistent Gram matrix, conditional invariants |
| GENERIC degeneracies | Finite-dimensional energy and entropy balances, then a specific discrete integrator |

Write definitions, hypotheses, theorem statements and counterexamples before coding proofs. Require a clean machine check without sorry/admit or unstated axioms. Track mathematical assumptions and the code-to-theorem correspondence separately. Lean verifies a statement, not runtime performance, floating-point behavior, scientific novelty or correctness of an unrelated solver. A theorem with unrealistic hypotheses is not sufficient evidence for production.

## 10. SPICE and WATER OF LIFE: retained grand vision

The mature Maker may couple collisionless Hamiltonian dynamics to appropriate dissipative gas physics, while SPICE exposes reusable mathematical contracts. Standard GENERIC for a closed system is

$$
\dot z=L(z)\nabla E+M(z)\nabla S,\quad
L^\top=-L,\quad M^\top=M\succeq0,\quad
L\nabla S=0,\quad M\nabla E=0.
$$

The reversible bracket must satisfy Jacobi where claimed. Degeneracy implies $\dot E=0$ and $\dot S=\nabla S^\top M\nabla S\ge0$ for autonomous closed systems. Dissipative evolution is not globally symplectic. Ideal collisionless Vlasov has no physical dissipative term. Numerical release/coarse-graining can lose fine-grained information but is not physical entropy production. Expanding cosmology, escaping radiation and other open exchanges require the appropriate energy/entropy fluxes, not naive closed-system conservation.

Reference oscillator: $z=(q,p,u)$, $E=p^2/(2m)+kq^2/2+u$, $S=C\ln(u/u_*)$, $T=u/C$, $v=(0,1,-p/m)^\top$, $M=\gamma T vv^\top$. The degeneracy $v^\top\nabla E=0$ and $M\nabla S$ give mechanical damping with reservoir heating. DeepSeek studied the oscillator and midpoint corrections; none of this establishes thermodynamically compatible astrophysical MHD or a universal adaptive metriplectic scheme.

SPICE should publish explicit operator contracts (domains, preconditions, positivity, conservation residuals, approximation error, derivatives when defined, and failure behavior). WATER OF LIFE owns physical irreversible closures and reservoir bookkeeping; MENTAT owns selection. The Maker must contain actual mathematics, not just an administrative control loop. Preserve the possibility of combining a Hamiltonian measure sector, suitable GENERIC gas sector and separately certified adaptation. A single global structure-preserving formulation remains an open research problem. Do not divert CHUÍ's ordinary MHD implementation to depend on this research.

## 11. Decision and proof ledger

| ID | Claim / decision | Evidence status |
|---|---|---|
| D1 | Production gravity remains N-body TreePM | Project decision |
| D2 | Label-space OT redistribution is off collisionless critical path; hydro OT independent | Negative tested arms / design decision |
| D3 | High-order spline is the practical sheet baseline; exact atlas is laboratory/selective | Recorded benchmark; no universal exclusion |
| D4 | Local block h-refinement with hierarchical surplus works as a geometric kernel | Author-derived restricted contract and recorded tests |
| D5 | All-sheet late multistream replacement meets original Gate B | FAILED |
| D6 | Early-collapse sheet can be more economical in the tested 2D regime | Recorded, conditional; not a general threshold |
| D7 | Same-lattice deposition and unique label ownership are key | Recorded for Call 1; hybrid extension proposed |
| D8 | Predictive refinement beats first-crossing gate in this test | Recorded |
| D9 | Releasing sheets gives the intended late-time cost/accuracy advantage | OPEN, not implemented |
| D10 | Positive averaging is symplectic | False for the specified class except trivial copying; restricted proof |
| D11 | Local atlas guarantees global adaptive Hamiltonian structure | Not established |
| D12 | Variational adjoint gravity resolves momentum/action/self-force tradeoff | OPEN |
| D13 | Gaussian/cloud/wave carriers improve over Dirac release | OPEN |
| D14 | Global hybrid error bound | Conditional abstract argument only; solver assumptions unverified |
| D15 | Lean certification, general 3D conforming geometry, expanding self-gravity and GENERIC coupling | OPEN |
| D16 | DUNE development deferred behind flagship physics/MHD | Project decision, 23 Nov branch discussion |
| D17 | Research only, no production integration or public feature claims | Retained charter |

Negative results are first-class deliverables. A future proposal may reopen a rejected arm only by explaining the changed mechanism and demonstrating fair improvements under controlled comparisons.

## 12. Resumption gates and acceptance strategy

The following is a proposed revised programme. It does not imply current work is scheduled or underway.

**Gate P — Pause and preserve (current).** Keep all source packages, hashes, reports, negative arms and continuation instructions. Complete flagship CHUÍ priorities independently. No DUNE production integration.

**Gate H0 — Reproduce and isolate.** Reproduce Call 1 smoke/regression tests and selected saved arms with the environment pinned. Resolve absent artifacts by regeneration, not invented results. Freeze the benchmark, synchronized momentum convention, label lattice and N-body control. Record source and environment hashes.

**Gate H1 — Predictive release, 2D.** Diagnose the ms lead anomaly using offline creation-momentum errors; add pre-caustic seeding. Implement deterministic one-way same-lattice release with a transaction/lineage ledger and shared PM gravity. Unit tests: zero/all/partial release; repeated requests; release at seam; ownership conflicts; no density double counting; conservation of sample mass and moments; rollback; restart; deterministic ordering. Compare initial exact sample release with alternative reduced sampling only as a separate arm. The release trigger is state-only. First show that the hybrid at t=3 can approach N-body accuracy with active sheet DOFs confined to the economical regime. Then require practical measured cost improvements, not just fewer sheet nodes.

**Gate H2 — Matched-cost scientific evaluation.** Compare lean N-body, uniform sheet, adaptive sheet, fixed-time release and causal predictive release on PM 256² and 512², multiple folding patterns and ICs. Include a sparse/local collapse rather than only the original worst case. Repeat one arm at a time and report variability. Use same mass, IC, field, timestep and reference; verify separate PM and dt convergence. Position and momentum comparisons must respect label correspondence; for general kernels use distribution/observable comparisons. Report density, force, power spectra, mass, momentum, energy balance, action defect where defined, released mass, cumulative transfer error, runtime/CPU/peak RSS/state bytes, quadrature count and actual total carrier count. Acceptance proposal: hybrid matches both position and density/force error targets at a specified reference accuracy with lower measured wall time or RSS, without unacceptable degradation in the other metric. A practical threshold should be fixed before runs (for example a 20% wall-time improvement at matched accuracy with no memory penalty, or a substantial measured memory reduction with no wall penalty). Do not silently change the original failed gate or use an extrapolated cost ratio as primary proof.

**Gate I — IX-GRAVITAS (parallel mathematical investigation when resumed).** Prove and finite-difference-check adjoint force and consistent mass matrix on fixed geometry; quantify self-force, momentum, action, energy and cost. Determine whether regularized/smoother field or variational force changes the late-time interpolation limitation. Apply to a hybrid only after standalone correctness. Compare the original nonvariational common-potential and momentum-conserving force arms fairly.

**Gate A — Advanced representations.** Only after ordinary release is validated, test moment/cloud representations under affine Hamiltonian flow and prescribed nonlinearity, then self-gravity. Include a cheap particle baseline and total memory overhead. Schrödinger, waterbag, higher-order Vlasov and hierarchical compression remain independent investigations. No automatic inverse conversion to a pristine sheet.

**Gate L — Formal foundations.** Formalize small restricted lemmas in Lean; refine the abstract error contract and transfer bounds. This can proceed independently in research but does not imply code validation. Revisit the single-law ambition once a concrete hybrid model exists.

**Gate 3D — Dimensional extension.** Conforming 3D/6D carrier, octree/face/edge/corner/periodic-seam ownership and compatible topology, 3D deposition, expanding Hamiltonian and Layzer–Irvine diagnostics, memory scaling and communication. Require 2D hybrid or specialized-sheet benefit in relevant regimes before expensive 3D work. For ordinary cosmological observables compare appropriate N-body references and resolution studies, not only a smooth toy problem.

**Gate G — Production and multiphysics, last.** Separate proposal, source audit, validation and explicit user decision. Respect production gravity and timestep ownership, patch-AMR/hydro, chemistry/cooling, MHD, MPI/OpenMP and snapshot/restart contracts. No silent enablement. GENERIC and gas OT belong in separate experiments until compatible coupling is derived and measured.

## 13. Implementation boundaries and operational contracts

Future C++ design must follow CHUÍ repository conventions (lower_snake_case paths, PascalCase types, camelCase methods, lower_snake_case variables). The Python work tree stays research-only. Use a single owner of gravity field and global timestep; an event scheduler cannot step physical time independently. Adopt concrete contiguous layouts and specialized hot kernels, not per-node virtual dispatch or unbounded allocations. Record geometry/node ownership separately from mass ownership.

The hybrid event journal should include event ID, step/time/frame, source label-region IDs and lineage, target representation and weights, physical and numerical pre/post residuals, trigger/threshold, field convention, source/config hashes, budget reservation, status/rollback and checkpoint schema. MPI eventually needs deterministic migration, global reductions and no double deposition across ranks; disabled features must leave ordinary CHUÍ behavior unchanged.

Memory planning must include simultaneous old/new state, rollback scratch, field buffers, ghost geometry, carrier operator storage, sample chunks, communication and checkpoint buffers, not just resident node arrays. Call 1's incremental RSS measurement is distinct from evolving-state memory. Rebuild/sparse storage overhead needs improvement before scaling. Maintain exact ownership and determinism contracts without claiming bitwise identity of parallel floating-point reductions unless actually tested.

## 14. Reproduction and source map

Repository source locations (relative to this file):

| Source | Purpose |
|---|---|
| [Original research README](README.md) | First manuscript, prescribed-Hamiltonian prototypes, 11-check reported verification |
| [Original proof ledger](proof_ledger.md) | Initial proof claims and restrictions |
| [v1.0 continuity pack](CHUI_DUNE_Research_Continuity_Pack_2026-11-16%281%29.zip) | Historical original sources and dossier |
| [v1.1 continuity pack](CHUI_DUNE_Research_Continuity_Pack_v1.1_2026-11-16.zip) | Prior master dossier, original Fable/DeepSeek packages and breakthrough data |
| [Call 1 final report](work/fable_adaptive_vp/FINAL_REPORT.md) | Measured verdict, implementation, known failures |
| [Call 1 local contract](work/fable_adaptive_vp/research/local_refinement_contract.md) | Geometry, mass, field, proof details |
| [Call 1 resume](work/fable_adaptive_vp/RESUME.md) | Exact continuation, APIs, caveats and experiments not to repeat |
| [Call 1 result summary](work/fable_adaptive_vp/results/SUMMARY.md) | Generated benchmark tables |
| [Call 1 gate JSON](work/fable_adaptive_vp/results/gate_evaluation.json) | Machine-readable comparison flags, metrics and provenance |
| [Call 1 code](work/fable_adaptive_vp/src/shadapt/) | forest, hcarrier, sheet, solver, controller, kernels, diagnostics |
| [Call 1 tests](work/fable_adaptive_vp/tests/) | Regression checks |
| [Call 1 experiment drivers](work/fable_adaptive_vp/experiments/) | exp00–exp08 and scoring |
| [Call 1 vendor](work/fable_adaptive_vp/vendor/) | Frozen breakthrough source and checksums |
| [Developer placement](../../developer_file_placement.md) | Repository conventions |

Reference blobs inspected for this update: FINAL_REPORT.md 980d04a098f50ec434f716e0272e201abf5d692c; RESUME.md 65e7d03c89aef44c342617c8dd917b40ec0411ec; local_refinement_contract.md 1e1309f6986a92ebc337807bc6f735d641d6a34f; gate_evaluation.json bcf6472845836b8515f90bf28aa958d508efcb39. The private conversation export is not uploaded; its 23 November decisions are incorporated as design history only.

Reproduction commands from repository root, adapted from Call 1 README/RESUME (not executed for v1.2):

    cd docs/research/shai_hulud/work/fable_adaptive_vp
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    python3 -m pytest -q tests -p no:cacheprovider
    cd experiments
    python3 exp01_reference.py 256 4096
    python3 exp01_reference.py 512 4096
    python3 exp03_gate.py 256 --repeat
    python3 exp03_gate.py 512 --repeat
    python3 analyze.py

Use the recorded package versions (Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, Numba 0.68.1, pytest 9.2.0 in the source campaign) when reproducing. Large reference snapshots are regenerable, git-ignored and not committed; the original breakthrough tree is inside v1.1's original sources for exp00. Preserve hygiene (no caches or reference binaries committed unintentionally). Read RESUME before any new experiment, especially synchronization before momentum inspection, sub0 supporting the highest refinement level, and floor_cells=0.0 for the specific ms deep-refinement comparison.

## 15. Final instruction to the next research session

Read this dossier, then FINAL_REPORT, local_refinement_contract, RESUME and the original v1.1 materials. Do not ask a new agent to reinvent the previous numerical tests. Start only when the project owner reactivates DUNE. First address the causal release experiment; develop the variational force and formal proof programme as separate research tasks. Produce a small falsifiable experiment, source/proof ledger, commands, raw results and an explicit verdict. Do not treat early success as evidence for all epochs or declare a general no-go from one 2D test. No unsupported novelty, no uncompiled Lean proof claims, no replacement of CHUÍ's flagship until a real validated advantage exists.

**The Maker's surviving question:** Can a numerical universe dynamically change the mathematical objects representing its matter, while controlling the physical information lost in those changes and retaining the useful structure of its evolution?

This remains an ambitious research programme, not an advertised production feature. The negative results are why it is now better defined.
