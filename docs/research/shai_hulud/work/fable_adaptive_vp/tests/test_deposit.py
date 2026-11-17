"""PM kernels and leaf-quadrature deposition: equivalence with the baseline
PM, exact mass ledger, single ownership of mass, periodic seams, invariance
of the deposit under refinement."""
import numpy as np
import pytest

from helpers import random_forest, smooth
from shadapt.kernels import PMGrid, deposit_pcs, deposit_pcs_w, gather_pcs_hess
from shadapt.sheet import AdaptiveSheet
from shadapt.solver import AdaptiveVP, Quadrature
from shsheet.pm import PM


def test_kernels_match_the_baseline_pm():
    rng = np.random.default_rng(0)
    x = rng.random((5000, 2)) * 3.0 - 1.0          # also outside [0,1): periodic wrap
    base = PM(32, 4, "grad")
    pm = PMGrid(32)
    acc = pm.new_density()
    deposit_pcs(acc, x, 1.0 / 5000)
    rho = pm.finalize_density(acc)
    rb = base.deposit(x, 1.0 / 5000)
    assert np.abs(rho - rb).max() < 1e-12 and abs(rho.mean() - 1.0) < 1e-13
    acc2 = pm.new_density()
    deposit_pcs_w(acc2, x, np.full(5000, 1.0 / 5000))
    assert np.abs(acc2 - acc).max() < 1e-18
    phi, phib = pm.potential(rho), base.potential(rb)
    assert np.abs(phi - phib).max() < 1e-14
    V, a = pm.field(phi, x)
    Vb, ab = base.field(phib, x)
    assert np.abs(V - Vb).max() < 1e-14 and np.abs(a - ab).max() < 1e-12
    assert abs(pm.energy_potential(rho, phi) - base.energy_potential(rb, phib)) < 1e-15


def test_hessian_kernel_is_the_derivative_of_the_force():
    rng = np.random.default_rng(1)
    pm = PMGrid(16)
    phi = rng.standard_normal((16, 16))
    x = rng.random((200, 2))
    H = np.empty((200, 3))
    gather_pcs_hess(phi, x, H)
    e = 1e-6
    for ax, cols in ((0, (0, 1)), (1, (1, 2))):
        d = np.zeros(2)
        d[ax] = e
        fd = -(pm.field(phi, x + d)[1] - pm.field(phi, x - d)[1]) / (2 * e)
        assert np.abs(fd - H[:, list(cols)]).max() < 1e-5 * np.abs(H).max()


def _sheet(f, amp=0.02):
    sh = AdaptiveSheet(f, k=5)
    for l in range(len(sh.keys)):
        z = smooth(sh.xi[l])
        sh.d[l], sh.p[l] = amp * z, 0.1 * z[:, ::-1]
    return sh


@pytest.mark.parametrize("mode", ("lattice", "adaptive"))
@pytest.mark.parametrize("seed", (0, 1, 2))
def test_mass_ledger_is_exact_on_any_topology(mode, seed):
    f = random_forest(32, 8, 3, nsplit=30, seed=seed)
    run = AdaptiveVP(_sheet(f), 64, Quadrature(mode, 4, target=0.5))
    for _ in range(3):                       # adaptive mode changes the sample counts between solves
        run.solve_field()
        assert run.mass_ledger == 1.0        # dyadic masses: exact, not just round-off
        assert abs(run.rho.mean() - 1.0) < 2e-13
        assert np.all(run.quad.nq >= run.sh.car.B) and np.all((run.quad.nq & (run.quad.nq - 1)) == 0)
    # every leaf block is visited by exactly one quadrature group
    seen = np.zeros(len(run.sh.car.leafL), dtype=int)
    for L, d, n1, n2, sel in run.quadrature_groups():
        seen[run.sh.car.groups[(L, d)][0][sel]] += 1
    assert np.all(seen == 1)


def test_uniform_density_of_the_identity_sheet_on_a_refined_mesh():
    """Cold sheet at rest: the deposit is exactly uniform whatever the
    refinement pattern (no spurious density at coarse-fine interfaces)."""
    f = random_forest(32, 8, 3, nsplit=30, seed=3)
    run = AdaptiveVP(AdaptiveSheet(f, k=5), 32, Quadrature("lattice", 8))   # 256^2 lattice on a 32^2 grid
    run.solve_field()
    assert np.abs(run.rho - 1.0).max() < 1e-12 and np.abs(run.phi).max() < 1e-13
    run.step(0.05)
    assert max(np.abs(d).max() for d in run.sh.d) < 1e-13     # no spurious force: nothing moves


def test_deposit_is_invariant_under_refinement():
    """Splitting blocks changes neither the carrier nor (in lattice mode) the
    sample lattice: the PM density is the same to round-off."""
    f = random_forest(32, 8, 3, nsplit=6, seed=4)
    run = AdaptiveVP(_sheet(f), 64, Quadrature("lattice", 8))
    run.solve_field()
    rho0, n0, s0 = run.rho.copy(), run.sh.n_nodes(), run.last_samples
    leaves = [b for b in run.sh.f.leaves() if b[0] < 2][:10]
    run.sh.refine(leaves, reason="test")
    run.solve_field()
    assert run.sh.n_nodes() > n0 and run.last_samples == s0
    assert np.abs(run.rho - rho0).max() < 1e-9 * rho0.max()
    assert run.mass_ledger == 1.0


def test_seam_crossing_and_periodicity():
    """Translating the whole sheet by an integer number of PM cells through
    the periodic seam rolls the density; a translation by a full box is the
    identity."""
    f = random_forest(32, 8, 2, nsplit=10, seed=5)
    run = AdaptiveVP(_sheet(f), 64, Quadrature("lattice", 4))
    run.solve_field()
    rho0 = run.rho.copy()
    for l in range(len(run.sh.keys)):
        run.sh.d[l] = run.sh.d[l] + [37.0 / 64, -5.0 / 64]
    run.solve_field()
    assert np.abs(run.rho - np.roll(np.roll(rho0, 37, 0), -5, 1)).max() < 1e-10
    for l in range(len(run.sh.keys)):
        run.sh.d[l] = run.sh.d[l] + [1.0 - 37.0 / 64, 2.0 + 5.0 / 64]
    run.solve_field()
    assert np.abs(run.rho - rho0).max() < 1e-10


def test_adaptive_quadrature_meets_its_spacing_target():
    f = random_forest(32, 8, 2, nsplit=10, seed=6)
    run = AdaptiveVP(_sheet(f, amp=0.08), 256, Quadrature("adaptive", 1, target=0.5, nq_max=256))
    for _ in range(4):
        run.solve_field()
    assert run._spacing_cells.max() <= 0.5 + 1e-12
    assert run.mass_ledger == 1.0
    assert len({tuple(v) for v in run.quad.nq.tolist()}) > 1      # genuinely non-uniform
