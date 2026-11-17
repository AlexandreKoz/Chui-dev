"""State-only refinement indicators and the refinement controller.

Everything here is a function of the evolving nodal state (x, p), of the PM
field computed from it, and of the clock.  No reference / oracle / future
trajectory enters any decision.

Indicators (one value per leaf block; all in units of PM cells):

  "d6"   truncation monitor.  For a degree-k spline the (k+1)-th undivided
         difference of the B-spline coefficients IS the jump of its k-th
         derivative across a knot, and  kappa_k |Delta^{k+1} c|  is the
         leading term of the interpolation error of the block at its own
         resolution (kappa_3 = 1/384, kappa_5 = 1/15360: periodic-spline
         error constants, measured in tests).  Evaluated for the
         displacement and for tau_p * momentum; the block value is the max.
         (For k = 3 this is the classical fourth-difference "d4" monitor.)
  "pred" predictive strain/force-informed monitor: the same functional
         applied to the second-order look-ahead of the node positions
              x_pred = x + tau p + tau^2/2 a(x),
         a = PM acceleration at the nodes (already gathered by the last
         kick).  It measures the truncation error the CURRENT mesh would
         have at time t + tau if the current strain rate and tidal field
         persisted, and therefore fires before the under-resolution exists.
  strain diagnostics (logged, and used for the trigger gates): det Dx and the
         ballistic time-to-caustic tau_c (smallest positive root of
         det(Dx + tau Dp) = 0) at the cell centres of every block.

Controller: refinement threshold `tol`, hysteresis (merge only if the
estimated parent-level error is below `hyst * tol`), hard node budget
`max_nodes` (never exceeded; denied requests are logged), cooldown (a block
younger than `cooldown` steps is neither split nor merged), resolution floor
(a block whose cells are already smaller than `floor_cells` PM cells is not
split), 2:1 balance, optional buffer ring, and a reason for every split.
Gates for the trigger-timing study: `t_enable` (clock) and `enable_on`
("always" | "precaustic" | "crossing").
"""
from __future__ import annotations

import copy

import numpy as np

KAPPA = {3: 1.0 / 384.0, 5: 1.0 / 15360.0}


class Controller:
    def __init__(self, indicator="pred", tol=1e-3, hyst=0.1, interval=5, cooldown=10, max_nodes=None,
                 tau=0.3, tau_p=0.3, floor_cells=0.5, t_enable=0.0, enable_on="always", tau_lead=0.3,
                 buffer=0, coarsen=False, max_splits_per_event=None):
        if indicator not in ("d6", "pred", "max"):
            raise ValueError(indicator)
        self.indicator, self.tol, self.hyst = indicator, float(tol), float(hyst)
        self.interval, self.cooldown, self.max_nodes = int(interval), int(cooldown), max_nodes
        self.tau, self.tau_p, self.floor_cells = float(tau), float(tau_p), float(floor_cells)
        self.t_enable, self.enable_on, self.tau_lead = float(t_enable), enable_on, float(tau_lead)
        self.buffer, self.coarsen = int(buffer), bool(coarsen)
        self.max_splits_per_event = max_splits_per_event
        self.latched = False            # gate latch for "precaustic" / "crossing"
        self.latch_step = None
        self.merged_at = {}             # (l, I, J) -> step of last merge
        self.history = []               # one record per evaluation
        self.first_split_t = None

    # ------------------------------------------------------------------ state
    def params(self):
        return dict(indicator=self.indicator, tol=self.tol, hyst=self.hyst, interval=self.interval, cooldown=self.cooldown,
                    max_nodes=self.max_nodes, tau=self.tau, tau_p=self.tau_p, floor_cells=self.floor_cells,
                    t_enable=self.t_enable, enable_on=self.enable_on, tau_lead=self.tau_lead, buffer=self.buffer,
                    coarsen=self.coarsen, max_splits_per_event=self.max_splits_per_event)

    def state_dict(self):
        return dict(params=self.params(), latched=self.latched, latch_step=self.latch_step,
                    merged_at=[[list(k), v] for k, v in sorted(self.merged_at.items())],
                    history=self.history, first_split_t=self.first_split_t)

    @classmethod
    def from_state(cls, st):
        c = cls(**st["params"])
        c.latched, c.latch_step = st["latched"], st["latch_step"]
        c.merged_at = {tuple(k): v for k, v in st["merged_at"]}
        c.history = list(st["history"])
        c.first_split_t = st["first_split_t"]
        return c

    # ------------------------------------------------------------------ indicators
    def indicators(self, run):
        """Per-leaf-block indicator arrays (state-only)."""
        sh, car, pm = run.sh, run.sh.car, run.pm
        k = car.k
        kap = KAPPA[k]
        acc = run.acc
        nlev = len(sh.keys)
        zs = []
        for l, sl in enumerate(sh.level_slices()):
            zs.append(np.concatenate([sh.d[l], sh.p[l], acc[sl]], axis=1))
        a = car.update(zs, tag="ind")
        nl = len(car.leafL)
        out = dict(eta_x=np.zeros(nl), eta_p=np.zeros(nl), eta_pred=np.zeros(nl), det_min=np.full(nl, np.inf),
                   tau_c=np.full(nl, np.inf), extent=np.zeros(nl))
        tau = self.tau
        B = car.B
        for (L, d), (ids, idx) in car.groups.items():
            A = a[L + d][idx]                                   # (nblk, mw, mw, 6)
            sc = kap * float(2 ** (k * d)) * pm.ng              # -> PM cells, at the block's own resolution

            def jump(F):
                d1 = np.abs(np.diff(F, n=k + 1, axis=1)).max(axis=(1, 2, 3))
                d2 = np.abs(np.diff(F, n=k + 1, axis=2)).max(axis=(1, 2, 3))
                return np.maximum(d1, d2)
            out["eta_x"][ids] = sc * jump(A[..., 0:2])
            out["eta_p"][ids] = sc * self.tau_p * jump(A[..., 2:4])
            out["eta_pred"][ids] = sc * jump(A[..., 0:2] + tau * A[..., 2:4] + 0.5 * tau * tau * A[..., 4:6])
            # strain at cell centres: Dx = 1 + grad d, Dp = grad p
            W, Dw = run._weights(d, B), run._weights(d, B, True)
            scale = float(1 << (L + d))
            t1 = np.tensordot(A[..., 0:4], Dw, axes=([1], [1]))  # (nblk, mw, 4, B) d/dxi1
            g1 = scale * np.tensordot(t1, W, axes=([1], [1]))    # (nblk, 4, B, B)
            t2 = np.tensordot(A[..., 0:4], W, axes=([1], [1]))
            g2 = scale * np.tensordot(t2, Dw, axes=([1], [1]))
            x11, x21, p11, p21 = 1.0 + g1[:, 0], g1[:, 1], g1[:, 2], g1[:, 3]
            x12, x22, p12, p22 = g2[:, 0], 1.0 + g2[:, 1], g2[:, 2], g2[:, 3]
            c0 = x11 * x22 - x12 * x21
            c1 = x11 * p22 + p11 * x22 - x12 * p21 - p12 * x21
            c2 = p11 * p22 - p12 * p21
            disc = c1 * c1 - 4.0 * c2 * c0
            with np.errstate(divide="ignore", invalid="ignore"):
                sq = np.sqrt(np.maximum(disc, 0.0))
                r1 = np.where(np.abs(c2) > 1e-300, (-c1 - sq) / (2 * c2), -c0 / c1)
                r2 = np.where(np.abs(c2) > 1e-300, (-c1 + sq) / (2 * c2), -c0 / c1)
            bad = (disc < 0) | ~np.isfinite(r1) | ~np.isfinite(r2)
            r1 = np.where(bad | (r1 <= 0), np.inf, r1)
            r2 = np.where(bad | (r2 <= 0), np.inf, r2)
            tc = np.where(c0 > 0, np.minimum(r1, r2), 0.0)       # already folded: 0
            out["tau_c"][ids] = tc.min(axis=(1, 2))
            out["det_min"][ids] = c0.min(axis=(1, 2))
            hL = 1.0 / sh.f.ncell(L)
            ext = np.maximum(np.maximum(np.abs(x11), np.abs(x21)), np.maximum(np.abs(x12), np.abs(x22))).max(axis=(1, 2))
            out["extent"][ids] = ext * hL * pm.ng                # largest cell edge in PM cells
        out["eta_d6"] = np.maximum(out["eta_x"], out["eta_p"])
        out["eta_max"] = np.maximum(out["eta_d6"], out["eta_pred"])
        return out

    def select(self, ind):
        return {"d6": ind["eta_d6"], "pred": ind["eta_pred"], "max": ind["eta_max"]}[self.indicator]

    # ------------------------------------------------------------------ decisions
    def maybe_adapt(self, run, dt):
        if self.interval <= 0 or run.step_count % self.interval != 0:
            return False
        sh, car, f = run.sh, run.sh.car, run.sh.f
        ind = self.indicators(run)
        eta = self.select(ind)
        step, t = run.step_count, run.t
        # gates (state-only latches)
        if not self.latched:
            if self.enable_on == "always":
                self.latched = True
            elif self.enable_on == "crossing" and float(ind["det_min"].min()) <= 0.0:
                self.latched = True
            elif self.enable_on == "precaustic" and float(ind["tau_c"].min()) <= self.tau_lead:
                self.latched = True
            if self.latched:
                self.latch_step = step
        enabled = self.latched and t >= self.t_enable - 1e-12
        rec = dict(step=step, t=t, eta_max=float(eta.max()), eta_d6_max=float(ind["eta_d6"].max()),
                   eta_pred_max=float(ind["eta_pred"].max()), det_min=float(ind["det_min"].min()),
                   tau_c_min=float(ind["tau_c"].min()), enabled=bool(enabled), nodes=sh.n_nodes(),
                   n_split=0, n_balance=0, n_buffer=0, n_merge=0, n_denied_budget=0, n_denied_cooldown=0,
                   n_denied_floor=0, n_denied_level=0)
        changed = False
        if enabled:
            L, I, J = car.leafL, car.leafI, car.leafJ
            age = np.array([step - max(f.lineage.get((int(l), int(i), int(j)), {}).get("step", -10 ** 9),
                                       self.merged_at.get((int(l), int(i), int(j)), -10 ** 9))
                            for l, i, j in zip(L, I, J)])
            over = eta > self.tol
            lvl_ok = L < f.lmax
            cool_ok = age >= self.cooldown
            floor_ok = ind["extent"] > self.floor_cells
            rec["n_denied_level"] = int((over & ~lvl_ok).sum())
            rec["n_denied_cooldown"] = int((over & lvl_ok & ~cool_ok).sum())
            rec["n_denied_floor"] = int((over & lvl_ok & cool_ok & ~floor_ok).sum())
            cand = np.nonzero(over & lvl_ok & cool_ok & floor_ok)[0]
            cand = cand[np.argsort(-eta[cand], kind="stable")]
            if self.max_splits_per_event is not None:
                cand = cand[: self.max_splits_per_event]
            requests = []                                   # (block, reason, info)
            seen = set()
            for q in cand.tolist():
                blk = (int(L[q]), int(I[q]), int(J[q]))
                info = dict(indicator=self.indicator, value=float(eta[q]), tol=self.tol, eta_x=float(ind["eta_x"][q]),
                            eta_p=float(ind["eta_p"][q]), eta_pred=float(ind["eta_pred"][q]),
                            det_min=float(ind["det_min"][q]), tau_c=float(ind["tau_c"][q]), extent_cells=float(ind["extent"][q]))
                requests.append((blk, "indicator:" + self.indicator, info))
                seen.add(blk)
            if self.buffer > 0:
                for (blk, _, _) in list(requests):
                    ring = [blk]
                    for _ in range(self.buffer):
                        ring = [(blk[0], a, b) for (_, i0, j0) in ring for (a, b) in f.neighbours(blk[0], i0, j0)]
                        for nb in ring:
                            if nb not in seen and f.is_leaf(*nb) and nb[0] < f.lmax:
                                requests.append((nb, "buffer", dict(caused_by=list(blk))))
                                seen.add(nb)
            # merges (hysteresis): parents whose four children are quiet leaves
            merges = []
            if self.coarsen:
                eta_of = {(int(l), int(i), int(j)): float(e) for l, i, j, e in zip(L, I, J, eta)}
                kfac = float(2 ** (car.k + 1))
                for l in range(f.lmax):
                    for (Ip, Jp) in sorted(f.refined[l]):
                        kids = [(l + 1, 2 * Ip + a, 2 * Jp + b) for a in (0, 1) for b in (0, 1)]
                        if not all(kd in eta_of for kd in kids):
                            continue
                        if any(step - f.lineage[kd]["step"] < self.cooldown for kd in kids):
                            continue
                        if kfac * max(eta_of[kd] for kd in kids) < self.hyst * self.tol and f.can_coarsen(l, Ip, Jp):
                            if not any(nb in seen for nb in [(l, a, b) for (a, b) in f.neighbours(l, Ip, Jp)] + kids):
                                merges.append(((l, Ip, Jp), dict(value=kfac * max(eta_of[kd] for kd in kids), tol=self.hyst * self.tol)))
            if requests or merges:
                old_nq = run.quad.as_dict(car) if run.quad.nq is not None else None
                budget = self.max_nodes
                n_before = len(f.events)

                def change(forest):
                    for (blk, _info) in merges:
                        if forest.can_coarsen(*blk):
                            forest.coarsen(*blk, reason="hysteresis", step=step, t=t, info=_info)
                            self.merged_at[blk] = step
                    for (blk, reason, info) in requests:
                        if not forest.is_leaf(*blk):
                            continue                      # already split by a balance cascade
                        if budget is not None:
                            trial = copy.deepcopy((forest.exists, forest.refined))
                            ev, lin, ver = len(forest.events), dict(forest.lineage), forest.version
                            forest.refine(*blk, reason=reason, step=step, t=t, info=info)
                            if forest.n_nodes() > budget:   # hard budget: roll this request back
                                forest.exists, forest.refined = trial
                                del forest.events[ev:]
                                forest.lineage, forest.version = lin, ver
                                rec["n_denied_budget"] += 1
                                continue
                        else:
                            forest.refine(*blk, reason=reason, step=step, t=t, info=info)
                res = sh.apply_topology_change(change)
                new_events = f.events[n_before:]
                rec["n_split"] = sum(1 for e in new_events if e["kind"] == "split" and e["reason"].startswith("indicator"))
                rec["n_balance"] = sum(1 for e in new_events if e["kind"] == "split" and e["reason"] == "balance")
                rec["n_buffer"] = sum(1 for e in new_events if e["kind"] == "split" and e["reason"] == "buffer")
                rec["n_merge"] = sum(1 for e in new_events if e["kind"] == "merge")
                changed = bool(new_events)
                if changed:
                    run.quad.sync(sh.car, old_nq)
                    rec["nodes_after"] = sh.n_nodes()
                    rec["created"], rec["removed"] = res["created"], res["removed"]
                    if self.first_split_t is None and (rec["n_split"] + rec["n_buffer"]) > 0:
                        self.first_split_t = t
        rec["levels"] = sh.nodes_per_level()
        self.history.append(rec)
        return changed
