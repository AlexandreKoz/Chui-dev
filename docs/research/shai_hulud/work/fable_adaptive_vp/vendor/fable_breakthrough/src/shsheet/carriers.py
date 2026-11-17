"""Uniform interface to every reconstruction ("carrier") compared in the
benchmarks.  A carrier maps reference coordinates xi (M, 2) to (x, p, Sigma)
using ONLY the nodal state of a Sheet (no oracle).

  lin, cub, quin    direct periodic spline interpolation of (x, p, Sigma)
  front3, front5    Fable's generating-function elements in the global
                    configuration chart, p = (Dx)^{-T} grad W  (vendored,
                    unmodified sheet2d.Carrier2; identity label map only)
  ac3, ac5          skeleton-exact carrier: spline + nodal Lagrangian tangent
                    fix + edge-action bubbles (exact edge actions, zero net
                    omega flux per cell; frame-free; not pointwise Lagrangian)
  csa3, csa5        conforming symplectic-atlas carrier (cell-wise exactly
                    Lagrangian, osculating unitary frame per cell)
  csa3x, csa3a      same with the configuration chart / best Arnold chart
  ac3g, csa3g, csa5g  position-preserving gauge: pure-momentum edge bubbles and
                    the configuration chart wherever the sheet is a well
                    conditioned graph over x, osculating frames near caustics
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from .atlas import AtlasCarrier
from .pspline import PSpline2, grid

_VENDOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "vendor", "fable_prototype")

ALL = ("lin", "cub", "quin", "front3", "front5", "ac3", "ac5", "csa3", "csa5", "csa3x", "csa3a", "ac3g", "csa3g", "csa5g")


class SplineCarrier:
    def __init__(self, sh, order):
        n = sh.n
        G = grid(n)
        Sper = sh.S - G @ sh.periods
        self.sp = PSpline2(np.concatenate([sh.x - G, sh.p, Sper[..., None]], axis=-1), order)
        self.periods = sh.periods

    def __call__(self, xi):
        v = self.sp(xi)
        return v[:, :2] + xi, v[:, 2:4], v[:, 4] + xi @ self.periods

    def memory_floats(self):
        return int(self.sp.c.size)


class FableFront:
    def __init__(self, sh, degree):
        if _VENDOR not in sys.path:
            sys.path.insert(0, _VENDOR)
        import sheet2d as F  # vendored, unmodified
        if np.any(sh.periods != 0):
            raise ValueError("Fable front carriers assume zero Liouville periods")
        self.car = F.Carrier2(F.Sheet2(grid(sh.n), sh.x, sh.p, sh.S), f"front{degree}")

    def __call__(self, xi):
        return self.car(xi)

    def memory_floats(self):
        n = self.car.n
        return int(2 * n * n + sum(a.size for row in self.car.Dn for a in row))


def make_carrier(name, sh, **kw):
    """Returns (carrier, setup_seconds)."""
    t0 = time.perf_counter()
    if name in ("lin", "cub", "quin"):
        car = SplineCarrier(sh, {"lin": 1, "cub": 3, "quin": 5}[name])
    elif name in ("front3", "front5"):
        car = FableFront(sh, int(name[-1]))
    elif name.startswith("ac"):
        g = dict(gauge="x") if name.endswith("g") else {}
        car = AtlasCarrier(sh.x, sh.p, sh.S, order=int(name[2]), periods=sh.periods, lagrangian=False, interior=0, **g, **kw)
    elif name.startswith("csa"):
        frame = {"": "osc", "x": "x", "a": "arnold", "g": "auto"}[name[4:]]
        g = dict(gauge="x") if name.endswith("g") else {}
        car = AtlasCarrier(sh.x, sh.p, sh.S, order=int(name[3]), periods=sh.periods, frame=frame, **g, **kw)
    else:
        raise ValueError(name)
    return car, time.perf_counter() - t0
