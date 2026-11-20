#!/usr/bin/env python3
"""KERNSYM20261120 B: the square is a symmetry-forced critical point of every spectral functional.

Heat sums H(tau, s) = sum exp(-s lambda^2) over the Dirac spectrum of the flat torus of unit area.
Every spectral functional of the original (zeta family, residue, finite part, determinant with a
mass) is an integral of H over s, so statements for all s hold for all of them, whatever the sign.
B1  antiperiodic structure: the gradient of H vanishes at tau = i for every s (fixed point of tau -> -1/tau)
B2  it does not vanish at the hexagonal point: for this spin structure the square is the only forced point
B3  periodic structure: forced at both the square and the hexagon
B4  along the mirror family H is largest at the square for every s, across it smallest: a saddle of
    fixed orientation for every functional with a positive weight in s; a sign flip keeps the saddle
B5  along the mirror family H falls monotonically away from the square for every s: no further
    critical point there
B6  the residue at s = 1 is the area term: it does not depend on tau at fixed area
"""
import numpy as np
from common import Checks

c = Checks("b_forced")
R = 60
mm, nn = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")

def H(x, y, s, a=0.5, b=0.5):
    lam2 = (2 * np.pi) ** 2 * (((mm + a) * x - (nn + b)) ** 2 + ((mm + a) * y) ** 2) / y      # unit area
    return float(2 * np.sum(np.exp(-s * lam2)))

def grad(x, y, s, a=0.5, b=0.5, h=1e-5):
    return np.array([(H(x + h, y, s, a, b) - H(x - h, y, s, a, b)) / (2 * h), (H(x, y + h, s, a, b) - H(x, y - h, s, a, b)) / (2 * h)])

ss = [0.01, 0.03, 0.1, 0.3, 1.0]
rho = (0.5, np.sqrt(3) / 2)
scale = {s: H(0, 1, s) for s in ss}
g_i = max(np.linalg.norm(grad(0, 1, s)) / scale[s] for s in ss)
g_rho = min(np.linalg.norm(grad(*rho, s)) / scale[s] for s in ss[1:])
c.ok("B1 antiperiodic structure: gradient zero at tau = i for every s", g_i < 1e-6, f"{g_i:.1e}")
c.ok("B2 antiperiodic structure: gradient not zero at the hexagonal point", g_rho > 1e-3, f"{g_rho:.1e}")
gp = max(max(np.linalg.norm(grad(0, 1, s, 0, 0)), np.linalg.norm(grad(*rho, s, 0, 0))) / H(0, 1, s, 0, 0) for s in ss)
c.ok("B3 periodic structure: gradient zero at both the square and the hexagon", gp < 1e-6, f"{gp:.1e}")

ts = np.exp(np.linspace(0, 1.2, 25))
xs = np.linspace(0, 0.3, 13)
along = all(H(0, t, s) <= H(0, 1, s) * (1 + 1e-12) for s in ss for t in ts)
across = all(H(x, 1, s) >= H(0, 1, s) * (1 - 1e-12) for s in ss for x in xs)
strict = all(H(0, 1.3, s) < H(0, 1, s) and H(0.15, 1, s) > H(0, 1, s) for s in ss[1:])
c.ok("B4 for every s: largest at the square along tau = i t, smallest across it (saddle of fixed orientation)", along and across and strict)
mono = all(np.all(np.diff([H(0, t, s) for t in ts]) <= 1e-12 * scale[s]) for s in ss)
c.ok("B5 for every s: monotone along the mirror family away from the square", mono)
weyl = [H(0, 1, 1e-3) * 1e-3 * 2 * np.pi, H(0.2, 1.7, 1e-3) * 1e-3 * 2 * np.pi]
c.ok("B6 small s: s H -> 1/(2 pi) = 4 c3 per unit area for every tau: the residue term carries no modulus",
     all(abs(w - 1) < 1e-6 for w in weyl), [round(w, 8) for w in weyl])
c.done()
