"""Self-consistent 2D periodic Vlasov-Poisson evolution of a cold sheet with
a locally h-adaptive high-order carrier, and the matched N-body PM control.

Frozen conventions (identical to the baseline `shsheet.vp.VPRun`):
  time variable t on a static background; units L = 1, total mass 1,
  4 pi G = 1 (lap phi = rho - 1); KDK with the field solved once per step
  after the drift; PCS assignment; exact-gradient ("grad") force law, i.e.
  every node and every particle is a tracer of ONE potential V(x, t).

Sheet arm.  The density is assembled by leaf-cell quadrature of the composite
carrier: every leaf block deposits its label mass exactly once, through a
tensor midpoint lattice with positive equal weights (`Quadrature`).  Nodes of
all levels are kicked with -grad V at their own positions.  This nodal force
is the gradient of the PM potential at the node; it is NOT the derivative of
the discrete potential energy with respect to the nodal position (the
adjoint of the quadrature), so the discrete energy is not exactly conserved.
"""
from __future__ import annotations

import json
import time

import numpy as np

from .forest import Forest
from .kernels import PMGrid, block_spacing, deposit_pcs, deposit_pcs_blocks, drift_inplace, eval_weights, kick_pcs
from .sheet import AdaptiveSheet


def _pow2_ceil(v):
    return np.exp2(np.ceil(np.log2(np.maximum(v, 1e-300))))


class Quadrature:
    """Per-leaf-block tensor midpoint quadrature.

    mode "lattice": nq = B * max(1, sub0 >> L) samples per block and
        dimension, i.e. the uniform deposit lattice (n0 sub0)^2 of the
        baseline when sub0 >= 2^L (the baseline rule; used for controls).
    mode "adaptive": nq1, nq2 (powers of two in [nq_min, nq_max]) chosen per
        block and per label direction from the CURRENT image of the block so
        that the largest spacing between neighbouring samples is at most
        `target` PM cells (max-norm), with hysteresis.  State-only.
    Masses: block mass / (nq1 nq2) per sample -- positive, dyadic, and the
    block masses sum to 1 exactly."""

    def __init__(self, mode="lattice", sub0=4, target=0.5, nq_min=None, nq_max=None, shrink=0.45, rule="midpoint"):
        if rule not in ("midpoint", "gauss2"):
            raise ValueError(rule)
        # rule "midpoint": nq equispaced midpoints across the block (a lattice rule: spectrally
        #   accurate for smooth periodic integrands when the lattice is GLOBALLY uniform, but only
        #   second order across an interface where the sample spacing changes);
        # rule "gauss2": nq/2 panels with the two Gauss-Legendre points each (equal positive
        #   weights, exact for cubics panel by panel: fourth order with or without uniformity).
        self.rule = rule
        self.mode, self.sub0, self.target, self.shrink = mode, int(sub0), float(target), float(shrink)
        self.nq_min, self.nq_max = nq_min, nq_max
        self.nq = None        # (nleaf, 2) int
        self.version = None

    def lattice_nq(self, car):
        B = car.B
        s = np.maximum(1, self.sub0 >> car.leafL)
        return np.stack([B * s, B * s], axis=1).astype(np.int64)

    def sync(self, car, old=None):
        """(Re)initialise per-block sample counts after a topology change.
        `old`: dict (L, I, J) -> (nq1, nq2) of the previous leaf set; children
        inherit half of the parent's counts (same sample lattice), parents of
        merged blocks twice the largest child count."""
        B = car.B
        lo = self.nq_min or B
        hi = self.nq_max or 8 * B
        if self.mode == "lattice":
            self.nq = self.lattice_nq(car)
        else:
            nq = np.empty((len(car.leafL), 2), dtype=np.int64)
            for q, (L, I, J) in enumerate(zip(car.leafL.tolist(), car.leafI.tolist(), car.leafJ.tolist())):
                if old is None:
                    v = (B * max(1, self.sub0 >> L),) * 2
                elif (L, I, J) in old:
                    v = old[(L, I, J)]
                elif (L - 1, I // 2, J // 2) in old:
                    par = old[(L - 1, I // 2, J // 2)]
                    v = (par[0] // 2, par[1] // 2)
                else:
                    kids = [old.get((L + 1, 2 * I + a, 2 * J + b)) for a in (0, 1) for b in (0, 1)]
                    kids = [c for c in kids if c is not None]
                    v = (2 * max(c[0] for c in kids), 2 * max(c[1] for c in kids)) if kids else (B, B)
                nq[q] = (min(max(v[0], lo), hi), min(max(v[1], lo), hi))
            self.nq = nq
        self.version = car.version

    def as_dict(self, car):
        return {(int(L), int(I), int(J)): (int(a), int(b)) for L, I, J, (a, b) in zip(car.leafL, car.leafI, car.leafJ, self.nq)}

    def adapt(self, car, spacing_cells):
        """spacing_cells (nleaf, 2): largest sample spacing of each block along
        each label direction, in PM cells, measured with the current counts."""
        if self.mode != "adaptive":
            return
        B = car.B
        lo = self.nq_min or B
        hi = self.nq_max or 8 * B
        ratio = spacing_cells / self.target
        up = ratio > 1.0
        down = ratio < self.shrink
        fac = np.where(up, _pow2_ceil(ratio), np.where(down, 0.5, 1.0))
        self.nq = np.clip(np.rint(self.nq * fac), lo, hi).astype(np.int64)


def sample_abscissae(nq, B, rule):
    """Positions (in own-level cell units, 0..B) of the nq samples of a block
    along one label direction.  All weights are equal (1/nq of the block)."""
    if rule == "midpoint":
        return (np.arange(nq) + 0.5) * (B / nq)
    if nq % 2:
        raise ValueError("gauss2 needs an even number of samples per block and direction")
    w = B / (nq // 2)
    c = (np.arange(nq // 2) + 0.5) * w
    off = 0.5 * w / np.sqrt(3.0)
    return np.stack([c - off, c + off], axis=1).ravel()


class AdaptiveVP:
    """KDK Vlasov-Poisson driver for an AdaptiveSheet.

    Kicks of consecutive steps use the same field and are merged (one force
    gather per step); `sync()` applies the pending half kick and must precede
    any use of the momenta (diagnostics, refinement decisions, output).  The
    result equals the plain half-kick sequence up to round-off."""

    def __init__(self, sheet, ng, quadrature=None, controller=None, chunk_points=1 << 15):
        self.sh = sheet
        self.pm = PMGrid(ng)
        self.quad = quadrature or Quadrature("lattice", 4)
        self.ctrl = controller
        self.chunk = int(chunk_points)
        self.t, self.step_count = 0.0, 0
        self.rho = self.phi = None
        self.acc = None                     # acceleration of all nodes at the last kick
        self._lag, self._dt_lag = False, 0.0
        self.cpu = dict(carrier=0.0, evaluate=0.0, deposit=0.0, poisson=0.0, gather=0.0, adapt=0.0)
        self.n_field_solves = 0
        self.deposit_samples = 0            # cumulative number of deposited samples (work measure)
        self.last_samples = 0
        self.mass_ledger = 0.0
        self.node_steps = 0
        self._spacing_cells = None
        self.quad.sync(sheet.car)

    # ------------------------------------------------------------------ density
    def _block_xi(self, L, ids, n1, n2):
        car = self.sh.car
        B = car.B
        N = float(self.sh.f.ncell(L))
        u1 = sample_abscissae(n1, B, self.quad.rule)
        u2 = sample_abscissae(n2, B, self.quad.rule)
        x1 = (car.leafI[ids][:, None] * B + u1[None, :]) / N
        x2 = (car.leafJ[ids][:, None] * B + u2[None, :]) / N
        return x1, x2

    def _weights(self, d, nq, deriv=False):
        """Evaluation matrix for nq midpoint samples across a block evaluated
        d levels finer; deriv=True: derivative w.r.t. the level-(L+d) GRID
        coordinate (multiply by N_{L+d} for d/dxi)."""
        car = self.sh.car
        key = ("q", d, int(nq), deriv, self.quad.rule)
        W = car._wcache.get(key)
        if W is None:
            B, k = car.B, car.k
            u = sample_abscissae(nq, B, self.quad.rule) * (1 << d)
            i0 = np.floor(u).astype(np.int64)
            _, ws = eval_weights(u - i0, k, 1 if deriv else 0)
            W = np.zeros((nq, (B << d) + k))
            rr = np.arange(nq)
            for q in range(k + 1):
                W[rr, i0 + q] = ws[1 if deriv else 0][q]
            car._wcache[key] = W
        return W

    def quadrature_groups(self):
        """Yield (L, d, n1, n2, sel) with `sel` indexing the (L, d) group."""
        car = self.sh.car
        if self.quad.version != car.version:
            self.quad.sync(car)
        for (L, d), (ids, _) in car.groups.items():
            nq = self.quad.nq[ids]
            code = nq[:, 0] * 100000 + nq[:, 1]
            for c in np.unique(code):
                sel = np.nonzero(code == c)[0]
                yield L, d, int(c // 100000), int(c % 100000), sel

    def deposit(self):
        """Assemble the PM density from the leaf quadrature.  Returns the
        accumulated mass grid (not yet normalised)."""
        sh, car, pm = self.sh, self.sh.car, self.pm
        t0 = time.perf_counter()
        a = sh.coeffs("x")
        t1 = time.perf_counter()
        self.cpu["carrier"] += t1 - t0
        if self.quad.version != car.version:
            self.quad.sync(car)
        acc = pm.new_density()
        ledger = 0.0
        nsamp = 0
        measure = self.quad.mode == "adaptive"
        spacing = np.zeros((len(car.leafL), 2)) if measure else None
        te = td = 0.0
        for L, d, n1, n2, sel in self.quadrature_groups():
            W1, W2 = self._weights(d, n1), self._weights(d, n2)
            mass = sh.f.cell_mass(L) * car.B ** 2 / (n1 * n2)
            per = max(1, self.chunk // (n1 * n2))
            for s0 in range(0, len(sel), per):
                ss = sel[s0:s0 + per]
                tA = time.perf_counter()
                ids, V = car.block_values_cf(a, L, d, W1, W2, ss)     # displacement (nblk, 2, n1, n2)
                x1, x2 = self._block_xi(L, ids, n1, n2)
                V[:, 0] += x1[:, :, None]
                V[:, 1] += x2[:, None, :]
                if measure:
                    sp = np.empty((len(ids), 2))
                    block_spacing(V, sp)
                    spacing[ids] = sp
                tB = time.perf_counter()
                deposit_pcs_blocks(acc, V, mass)
                tC = time.perf_counter()
                te += tB - tA
                td += tC - tB
                ledger += mass * V.shape[0] * n1 * n2
                nsamp += V.shape[0] * n1 * n2
        self.cpu["evaluate"] += te
        self.cpu["deposit"] += td
        self.mass_ledger = ledger
        self.last_samples = nsamp
        self.deposit_samples += nsamp
        if measure:
            self._spacing_cells = spacing * pm.ng
        return acc

    def solve_field(self):
        acc = self.deposit()
        t0 = time.perf_counter()
        self.rho = self.pm.finalize_density(acc)
        self.phi = self.pm.potential(self.rho)
        self.cpu["poisson"] += time.perf_counter() - t0
        self.n_field_solves += 1
        # state-only update of the sample counts for the NEXT field solve
        if self.quad.mode == "adaptive":
            self.quad.adapt(self.sh.car, self._spacing_cells)

    # ------------------------------------------------------------------ KDK
    def _kick(self, kappa):
        t0 = time.perf_counter()
        sh = self.sh
        V, a = self.pm.field(self.phi, sh.x_nodes())
        self.acc = a
        for l, sl in enumerate(sh.level_slices()):
            sh.p[l] = sh.p[l] + kappa * a[sl]
            sh.S[l] = sh.S[l] - kappa * V[sl]
        self.cpu["gather"] += time.perf_counter() - t0

    def _drift(self, dt):
        sh = self.sh
        for l in range(len(sh.keys)):
            dn = sh.d[l] + dt * sh.p[l]
            sh.S[l] = sh.S[l] + np.sum(sh.p[l] * (dn - sh.d[l]), axis=-1) - 0.5 * dt * np.sum(sh.p[l] ** 2, axis=-1)
            sh.d[l] = dn

    def sync(self):
        """Apply the pending half kick: momenta are at time t afterwards."""
        if self._lag:
            self._kick(0.5 * self._dt_lag)
            self._lag = False

    def step(self, dt):
        if self.phi is None:
            self.solve_field()
        if self._lag and self._dt_lag != dt:
            self.sync()
        self._kick(dt if self._lag else 0.5 * dt)
        self._drift(dt)
        self.solve_field()
        self._lag, self._dt_lag = True, dt
        self.t += dt
        self.step_count += 1
        self.node_steps += self.sh.n_nodes()
        if self.ctrl is not None and self.ctrl.due(self.step_count):
            t0 = time.perf_counter()
            self.sync()
            changed = self.ctrl.maybe_adapt(self, dt)
            if changed:
                self.solve_field()      # same carrier (identity transfer); refresh the field on the new mesh
            self.cpu["adapt"] += time.perf_counter() - t0

    # ------------------------------------------------------------------ diagnostics
    def sheet_integrals(self):
        """Label-measure integrals of the carrier by leaf quadrature:
        mass, momentum, kinetic energy, folded (det Dx < 0) mass fraction,
        minimum of det Dx.  Uses the current quadrature lattice."""
        self.sync()
        sh, car = self.sh, self.sh.car
        a = sh.coeffs("xp")
        if self.quad.version != car.version:
            self.quad.sync(car)
        mass = 0.0
        P = np.zeros(2)
        K = 0.0
        folded = 0.0
        detmin = np.inf
        for L, d, n1, n2, sel in self.quadrature_groups():
            W1, W2 = self._weights(d, n1), self._weights(d, n2)
            D1, D2 = self._weights(d, n1, True), self._weights(d, n2, True)
            m = sh.f.cell_mass(L) * car.B ** 2 / (n1 * n2)
            scale = float(sh.f.ncell(L + d))      # d/dxi = N_{L+d} d/du
            per = max(1, self.chunk // (n1 * n2))
            for s0 in range(0, len(sel), per):
                ss = sel[s0:s0 + per]
                _, V = car.block_values(a, L, d, W1, W2, ss)
                _, G1 = car.block_values(a, L, d, D1, W2, ss)
                _, G2 = car.block_values(a, L, d, W1, D2, ss)
                p = V[..., 2:4]
                x11 = 1.0 + scale * G1[..., 0]
                x21 = scale * G1[..., 1]
                x12 = scale * G2[..., 0]
                x22 = 1.0 + scale * G2[..., 1]
                det = x11 * x22 - x12 * x21
                mass += m * det.size
                P += m * p.reshape(-1, 2).sum(axis=0)
                K += 0.5 * m * float(np.sum(p * p))
                folded += m * float((det < 0).sum())
                detmin = min(detmin, float(det.min()))
        return dict(mass=mass, Px=float(P[0]), Py=float(P[1]), K=K, folded_fraction=folded, det_min=detmin)

    def conserved(self):
        if self.phi is None:
            self.solve_field()
        self.sync()
        s = self.sheet_integrals()
        U = self.pm.energy_potential(self.rho, self.phi)
        s.update(t=self.t, U=U, E=s["K"] + U, grid_mass=float(self.rho.mean()), ledger_mass=float(self.mass_ledger))
        return s

    def memory_report(self):
        sh, car = self.sh, self.sh.car
        mem = car.memory_floats(2)
        return dict(nodes=sh.n_nodes(), nodes_per_level=sh.nodes_per_level(), leaf_cells=sh.f.n_leaf_cells(),
                    state_floats=sh.state_floats(), carrier_coefficient_floats=mem["coefficients"],
                    carrier_operator_floats=mem["operators"], carrier_window_ints=mem["windows"],
                    pm_floats=4 * self.pm.ng ** 2, samples_last=self.last_samples)

    # ------------------------------------------------------------------ checkpoints
    def save_checkpoint(self, path):
        sh = self.sh
        arrs = dict(t=self.t, step_count=self.step_count, ng=self.pm.ng, k=sh.k, periods=sh.periods,
                    lag=np.array([1.0 if self._lag else 0.0, self._dt_lag]), chunk=self.chunk,
                    forest=json.dumps(sh.f.to_dict()), nlev=len(sh.keys), quad_nq=self.quad.nq,
                    quad=json.dumps(dict(mode=self.quad.mode, sub0=self.quad.sub0, target=self.quad.target,
                                         nq_min=self.quad.nq_min, nq_max=self.quad.nq_max, shrink=self.quad.shrink,
                                         rule=self.quad.rule)),
                    car=json.dumps(dict(halo_tiles=sh.car.ht, solver=sh.car.solver, sweeps=sh.car.sweeps)),
                    ctrl=json.dumps(self.ctrl.state_dict() if self.ctrl is not None else None),
                    counters=json.dumps(dict(n_field_solves=self.n_field_solves, deposit_samples=self.deposit_samples,
                                             node_steps=self.node_steps, cpu=self.cpu)))
        for l in range(len(sh.keys)):
            arrs[f"keys{l}"], arrs[f"d{l}"], arrs[f"p{l}"], arrs[f"S{l}"] = sh.keys[l], sh.d[l], sh.p[l], sh.S[l]
        if self.phi is not None:
            # the field actually used by the run (with adaptive quadrature it is NOT
            # recomputable bit-for-bit from the state after the sample counts moved on)
            arrs["rho"], arrs["phi"] = self.rho, self.phi
            arrs["mass_ledger"], arrs["last_samples"] = self.mass_ledger, self.last_samples
        if self.acc is not None:
            arrs["acc"] = self.acc
        np.savez(path, **arrs)

    @classmethod
    def load_checkpoint(cls, path, controller_factory=None):
        z = np.load(path, allow_pickle=False)
        f = Forest.from_dict(json.loads(str(z["forest"])))
        cc = json.loads(str(z["car"]))
        sh = AdaptiveSheet(f, k=int(z["k"]), periods=z["periods"], halo_tiles=cc["halo_tiles"], solver=cc["solver"], sweeps=cc["sweeps"])
        for l in range(int(z["nlev"])):
            assert np.array_equal(sh.keys[l], z[f"keys{l}"]), "checkpoint node set does not match the forest"
            sh.d[l], sh.p[l], sh.S[l] = z[f"d{l}"].copy(), z[f"p{l}"].copy(), z[f"S{l}"].copy()
        q = json.loads(str(z["quad"]))
        quad = Quadrature(q["mode"], q["sub0"], q["target"], q["nq_min"], q["nq_max"], q["shrink"], q.get("rule", "midpoint"))
        cs = json.loads(str(z["ctrl"]))
        ctrl = controller_factory(cs) if (controller_factory is not None and cs is not None) else None
        run = cls(sh, int(z["ng"]), quad, ctrl, chunk_points=int(z["chunk"]))
        run._lag, run._dt_lag = bool(z["lag"][0]), float(z["lag"][1])
        run.quad.nq = z["quad_nq"].copy()
        run.quad.version = sh.car.version
        run.t, run.step_count = float(z["t"]), int(z["step_count"])
        c = json.loads(str(z["counters"]))
        run.n_field_solves, run.deposit_samples, run.node_steps = c["n_field_solves"], c["deposit_samples"], c["node_steps"]
        run.cpu.update(c["cpu"])
        if "phi" in z.files:
            run.rho, run.phi = z["rho"].copy(), z["phi"].copy()
            run.mass_ledger, run.last_samples = float(z["mass_ledger"]), int(z["last_samples"])
        if "acc" in z.files:
            run.acc = z["acc"].copy()
        return run


class NBodyVP:
    """Conventional N-body PM with the SAME kernels, force law and KDK:
    n x n equal-mass particles on the label lattice (the baseline "nodes"
    deposit).  Lean implementation for fair cost comparisons: in-place
    compiled kick and drift, one force gather per step (merged kicks), two
    arrays (x, p) of state.  `track_action=True` reproduces the baseline's
    per-particle action transport (used only by equivalence tests)."""

    def __init__(self, ic, n, ng, track_action=False):
        g = np.arange(n) / n
        q = np.stack(np.meshgrid(g, g, indexing="ij"), axis=-1).reshape(-1, 2)
        x, p, S = ic(q)
        self.n, self.x, self.p = n, np.ascontiguousarray(x), np.ascontiguousarray(p)
        self.S = S.copy() if track_action else None
        self.pm = PMGrid(ng)
        self.t, self.step_count = 0.0, 0
        self.rho = self.phi = None
        self._lag, self._dt_lag = False, 0.0
        self.cpu = dict(deposit=0.0, poisson=0.0, gather=0.0, drift=0.0)
        self.n_field_solves = 0
        self.deposit_samples = 0

    def solve_field(self):
        t0 = time.perf_counter()
        acc = self.pm.new_density()
        deposit_pcs(acc, self.x, 1.0 / (self.n * self.n))
        t1 = time.perf_counter()
        self.rho = self.pm.finalize_density(acc)
        self.phi = self.pm.potential(self.rho)
        self.cpu["deposit"] += t1 - t0
        self.cpu["poisson"] += time.perf_counter() - t1
        self.n_field_solves += 1
        self.deposit_samples += self.n * self.n

    def _kick(self, kappa):
        t0 = time.perf_counter()
        if self.S is None:
            kick_pcs(self.phi, self.x, self.p, kappa)
        else:
            V, a = self.pm.field(self.phi, self.x)
            self.p += kappa * a
            self.S -= kappa * V
        self.cpu["gather"] += time.perf_counter() - t0

    def sync(self):
        if self._lag:
            self._kick(0.5 * self._dt_lag)
            self._lag = False

    def step(self, dt):
        if self.phi is None:
            self.solve_field()
        if self._lag and self._dt_lag != dt:
            self.sync()
        self._kick(dt if self._lag else 0.5 * dt)
        t0 = time.perf_counter()
        if self.S is None:
            drift_inplace(self.x, self.p, dt)
        else:
            xn = self.x + dt * self.p
            self.S += np.sum(self.p * (xn - self.x), axis=-1) - 0.5 * dt * np.sum(self.p * self.p, axis=-1)
            self.x = xn
        self.cpu["drift"] += time.perf_counter() - t0
        self.solve_field()
        self._lag, self._dt_lag = True, dt
        self.t += dt
        self.step_count += 1

    def conserved(self):
        if self.phi is None:
            self.solve_field()
        self.sync()
        m = 1.0 / (self.n * self.n)
        K = 0.5 * m * float(np.sum(self.p ** 2))
        U = self.pm.energy_potential(self.rho, self.phi)
        P = m * self.p.sum(axis=0)
        return dict(t=self.t, K=K, U=U, E=K + U, Px=float(P[0]), Py=float(P[1]), mass=float(self.rho.mean()))

    def folded_fraction(self):
        """Fraction of label cells with negative finite-difference det Dx."""
        n = self.n
        x = self.x.reshape(n, n, 2)
        a = np.roll(x, -1, 0) - x
        a[-1, :, 0] += 1.0
        b = np.roll(x, -1, 1) - x
        b[:, -1, 1] += 1.0
        det = (a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]) * n * n
        return float((det < 0).mean()), float(det.min())

    def memory_report(self):
        return dict(particles=self.n * self.n, state_floats=4 * self.n * self.n, pm_floats=4 * self.pm.ng ** 2)
