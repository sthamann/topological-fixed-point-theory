#!/usr/bin/env python3
"""Step 6 (KERNSIX20261120): spacetime and gravity after step 3.

A  the entropy of the face state is finite and proportional to the area of the seam: S = (9 zeta(3)/2) c3 A / ell^2
B  the same coefficient on the round sphere (it only uses the Weyl law)
C  written as S = A/(4 G): G = ell^2 / (18 zeta(3) c3); the width ell is free, so G is not predicted
D  width zero (the flat kernel state): log 2 per mode, the coefficient grows with the cutoff
E  the modular generator commutes with translations: an equilibrium at one temperature, no horizon
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import zeta
from common import Checks

c = Checks("6_gravity")
c3 = 1 / (8 * np.pi)

def s(x):
    x = np.minimum(np.abs(x), 700.0)
    return np.log1p(np.exp(-x)) + x / (np.exp(x) + 1.0)     # entropy of one mode with Fermi weight 1/(1 + e^x)

I, _ = quad(lambda x: x * s(x), 0, 80)
c.ok("A0 integral of x s(x) over the half line = (9/4) zeta(3)", abs(I - 2.25 * zeta(3)) < 1e-9, round(I, 9))
coeff = 4.5 * zeta(3) * c3

# ---------------------------------------------------------------- A
R = 400
mm, nn = np.meshgrid(np.arange(-R, R), np.arange(-R, R), indexing="ij")
kabs = 2 * np.pi * np.sqrt((mm + 0.5) ** 2 + (nn + 0.5) ** 2)       # torus of side 1, antiperiodic; eigenvalues +-|k|
ratios = []
for ell in (0.04, 0.02, 0.01):
    S = 2 * np.sum(s(2 * ell * kabs))
    ratios.append(S * ell ** 2 / coeff)
c.ok("A1 flat torus: S ell^2 / A tends to (9 zeta(3)/2) c3 = 0.2152 as the width gets small against the size",
     abs(ratios[-1] - 1) < 2e-3 and abs(ratios[-1] - 1) <= abs(ratios[0] - 1) + 1e-4, [round(float(r), 5) for r in ratios])

# ---------------------------------------------------------------- B
ls = np.arange(0, 20000)
ratios_s = []
for ell in (0.04, 0.02, 0.01):                                     # sphere of radius 1, area 4 pi
    S = np.sum(4 * (ls + 1) * s(2 * ell * (ls + 1)))
    ratios_s.append(S * ell ** 2 / (coeff * 4 * np.pi))
c.ok("B1 round sphere: the same coefficient (the entropy only uses the density of modes per area)",
     abs(ratios_s[-1] - 1) < 2e-3, [round(float(r), 5) for r in ratios_s])

# ---------------------------------------------------------------- C
G_over_ell2 = 1 / (4 * coeff)
c.ok("C1 as S = A/(4 G): G = ell^2 / (18 zeta(3) c3) = 1.1615 ell^2", abs(G_over_ell2 - 1 / (18 * zeta(3) * c3)) < 1e-12
     and abs(G_over_ell2 - 1.1615) < 1e-3, round(float(G_over_ell2), 5))
c.ok("C2 the entropy needs no cutoff: the width of the collar plays the role of the Planck length; "
     "the original datum fixes no width, so no value of G follows", True, "ell free")

# ---------------------------------------------------------------- D
counts = [int(np.sum(kabs < L)) * 2 for L in (200.0, 400.0)]
c.ok("D1 width zero: log 2 per mode; the number of modes below a cutoff L is A L^2/(2 pi) = 4 c3 A L^2 and grows without bound",
     all(abs(n / (L ** 2 / (2 * np.pi)) - 1) < 0.01 for n, L in zip(counts, (200.0, 400.0))), counts)

# ---------------------------------------------------------------- E
Nx = 64
shift = np.roll(np.eye(Nx), 1, axis=0)
p = (shift - shift.T) / (2j)
Z = np.zeros((Nx, Nx))
B = np.block([[Z, p], [p, Z]])
T = np.block([[shift, Z], [Z, shift]])
ell = 0.7
w, v = np.linalg.eigh(B)
C0 = (v / (1 + np.exp(-2 * ell * w))) @ v.conj().T
wc, vc = np.linalg.eigh(C0)
K = (vc * np.log((1 - wc) / wc)) @ vc.conj().T
c.ok("E1 the modular generator of the face state commutes with the translations of the seam: one temperature everywhere",
     np.allclose(K @ T, T @ K, atol=1e-8))
x = np.diag(np.tile(np.arange(Nx) - Nx / 2, 2))
boost = (x @ B + B @ x) / 2
c.ok("E2 it is not a boost generator (which would grow with the distance from an edge): no horizon, no Unruh reading",
     np.linalg.norm(K - (-2 * ell) * B) < 1e-8 and np.linalg.norm(K @ T - T @ K) < 1e-8 < np.linalg.norm(boost @ T - T @ boost))
c.ok("E3 what a derivation of gravity along Jacobson would still need: a width that varies from place to place (a lapse), "
     "a horizon reading of the flow, and a value of the width", True, "open")

c.done({"coefficient_S_ell2_over_A": coeff, "G_over_ell2": G_over_ell2, "torus": ratios, "sphere": ratios_s})
