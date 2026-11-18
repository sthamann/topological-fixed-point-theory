#!/usr/bin/env python3
"""Part S: SPACETIME.PROCESS.LATTICE.01 -- growth of the TFPT process graph, no dimension assumed.

Primitive operations: the directed mark visits gamma_0..gamma_3 of the marked seam (Q1), with the
surface relation gamma_0 gamma_1 gamma_2 gamma_3 = 1.  The process graph depends on WHICH histories
are physically distinguishable (Universalraum principle).  Three quotients:
  (a) full words (fundamental group F_3)          -> ball growth exponential: not a local space
  (b) abelian windings (H_1 = Z^4/(1,1,1,1))       -> three unbounded directions
  (c) only the finite physical monodromies (A4)    -> finite graph, dimension 0
  S1 (a) exponential growth B(r+1)/B(r) -> constant > 1.
  S2 (c) the transported v117 monodromies of the four visits close on 12 elements: B(r) saturates.
  S3 (b) Hausdorff dimension of the winding lattice: log B / log r -> 3.
  S4 (b) spectral dimension of the random walk on the winding lattice: return probability ~ t^-3/2.
  S5 (b) Weyl exponent of the precision operator sqrt(lattice Laplacian) in the IR window: 3.
  S6 (b) geometry: positions sum n_a w_a with w_a the four family-line vectors (a tetrahedron) span the
     bcc lattice; its point group contains S4 = W(A3) x inversion = O_h: the Laplacian symbol is
     ISOTROPIC at order p^2 (single c), anisotropic only at order p^4 (irrelevant in the IR).
  S7 windings are physically resolved iff a test of INFINITE order exists: a continuous family-Cartan
     holonomy z (as in W_S(z)) with generic z separates all windings; rational z, the order-3 flavor
     monodromy or the nu = 1/4 holonomy only see a finite quotient.
  S8 (statement) no TFPT process yet couples the seam time L0 to these windings: dispersion and
     Lorentz signature are open (SPACETIME.LORENTZ.01).
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/process_lattice.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# (a) full words: F_3 on gamma_0,1,2 with gamma_3 = (gamma_0 gamma_1 gamma_2)^-1 ------------------------
def reduce(word):
    out = []
    for g in word:
        if out and out[-1] == -g:
            out.pop()
        else:
            out.append(g)
    return tuple(out)
gens = [(1,), (2,), (3,), (-1,), (-2,), (-3,), (-3, -2, -1), (1, 2, 3)]   # +-gamma_0..2 (labels 1..3), gamma_3^{+-1}
ball, frontier, B_free = {()}, [()], [1]
for r in range(1, 7):
    new = []
    for w in frontier:
        for g in gens:
            x = reduce(w + g)
            if x not in ball:
                ball.add(x); new.append(x)
    frontier = new
    B_free.append(len(ball))
ratios = [B_free[k + 1] / B_free[k] for k in range(3, 6)]
step("S1.(a) full-word process graph (F_3, 8 generators): exponential ball growth -> not a local space",
     min(ratios) > 4, "B(r) = " + ", ".join(map(str, B_free)) + f"; ratios {np.round(ratios, 2).tolist()}")

# (c) finite physical quotient: actual transported monodromies --------------------------------------
from tfpt_explorer.flavor_path_transport import original_fuchs_connection
from scipy.integrate import solve_ivp
def transport(segs):
    T = np.eye(3, dtype=complex)
    for path in segs:
        def rhs(t, y, path=path):
            zz, dz = path(t)
            return (dz * original_fuchs_connection(zz) @ y.reshape(3, 3)).ravel()
        T = solve_ivp(rhs, (0, 1), np.eye(3, dtype=complex).ravel(), method="DOP853", rtol=1e-12, atol=1e-14).y[:, -1].reshape(3, 3) @ T
    return T
ray = lambda a, b: (lambda t: (a + (b - a) * t, b - a))
circ = lambda c, rr, s0: (lambda t: (c + rr * np.exp(1j * (s0 + 2 * np.pi * t)), 2j * np.pi * rr * np.exp(1j * (s0 + 2 * np.pi * t))))
M = [transport([ray(0, 0.6 * 1j ** k), circ(1j ** k, 0.4, np.angle(-(1j ** k))), ray(0.6 * 1j ** k, 0)]) for k in range(4)]
key = lambda X: tuple(np.round(X, 6).ravel())
seen, front, B_fin = {key(np.eye(3)): np.eye(3)}, [np.eye(3, dtype=complex)], [1]
for r in range(1, 8):
    new = []
    for X in front:
        for G in M + [np.linalg.inv(G) for G in M]:
            Y = G @ X
            if key(Y) not in seen:
                seen[key(Y)] = Y; new.append(Y)
    front = new
    B_fin.append(len(seen))
step("S2.(c) quotient by the actual finite flavor monodromy: B(r) saturates at 12 (A4) -> dimension 0",
     B_fin[-1] == 12 and B_fin[-1] == B_fin[-2], "B(r) = " + ", ".join(map(str, B_fin)))

# (b) abelian windings: Z^4/(1,1,1,1) ----------------------------------------------------------------
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])     # family-line vectors (tetrahedron)
steps8 = [tuple(s * x for x in v) for v in w for s in (1, -1)]
def lattice_ball(R):
    key3 = lambda p: tuple(np.round(np.array(p) * 2).astype(int))
    seenp, frontp, B = {key3((0, 0, 0))}, [np.zeros(3)], [1]
    for r in range(1, R + 1):
        newp = []
        for p in frontp:
            for s in steps8:
                q = p + np.array(s)
                kq = key3(q)
                if kq not in seenp:
                    seenp.add(kq); newp.append(q)
        frontp = newp
        B.append(len(seenp))
    return B
Bz = lattice_ball(40)
rs = np.arange(20, 41)
dH = np.polyfit(np.log(rs), np.log(np.array(Bz)[rs]), 1)[0]
step("S3.(b) winding lattice H_1 = Z^4/(1,1,1,1): ball growth B(r) ~ r^3 -> Hausdorff dimension 3",
     abs(dH - 3) < 0.1, f"d_H = {dH:.3f}; B(40) = {Bz[40]}")

L = 96
k = 2 * np.pi * np.fft.fftfreq(L)
KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
# lattice basis: bcc = Z^3 + Z^3 + (1/2,1/2,1/2); use cubic coordinates doubled so steps are integer
phase = lambda v: KX * (2 * v[0]) + KY * (2 * v[1]) + KZ * (2 * v[2])
mu = sum(np.cos(phase(v)) for v in w) / 4                          # characteristic function of the walk
ts = np.array([16, 32, 64, 128])
pt = np.array([np.mean(mu ** t) for t in ts])
ds = -2 * np.polyfit(np.log(ts), np.log(pt), 1)[0]
step("S4.(b) spectral dimension of the random walk on the winding lattice: return probability ~ t^-3/2 -> d_s = 3",
     abs(ds - 3) < 0.15, f"d_s = {ds:.3f}")
lam = sum(2 * (1 - np.cos(phase(v))) for v in w).ravel()
P = np.sqrt(np.sort(lam))
Lams = np.array([0.15, 0.25, 0.4])
NL = np.array([np.searchsorted(P, x, side="right") for x in Lams])
weyl = np.polyfit(np.log(Lams), np.log(NL), 1)[0]
step("S5.(b) Weyl exponent of the precision operator sqrt(Laplacian) on the winding lattice in the IR window: 3",
     abs(weyl - 3) < 0.15, f"exponent {weyl:.3f}")

zs = np.random.default_rng(7).normal(size=(5, 3))
vecs = np.vstack([w, -w])
m2 = [np.sum((vecs @ z) ** 2) / np.dot(z, z) for z in zs]
m4 = [np.sum((vecs @ z) ** 4) / np.dot(z, z) ** 2 for z in zs]
span = {tuple(np.round(2 * (a * w[0] + b * w[1] + c * w[2]), 6)) for a, b, c in product(range(-2, 3), repeat=3)}
bcc = (0.0, 0.0, 2.0) in span and (-1.0, -1.0, -1.0) in span
step("S6.(b) the windings span the bcc lattice; O_h = S4 x inversion makes the symbol ISOTROPIC at p^2 (single c), anisotropic only at p^4",
     bcc and np.allclose(m2, m2[0]) and np.ptp(m4) > 1e-3, f"p^2 coefficient {m2[0]:.4f} for all directions; p^4 spread {np.ptp(m4):.3f}")

from fractions import Fraction as Fr
wins = list(product(range(-4, 5), repeat=3))                    # windings n_0, n_1, n_2 (n_3 fixed by the relation)
pos = lambda n: n[0] * w[0] + n[1] * w[1] + n[2] * w[2]
def distinct_float(z):
    return len({int(round((np.dot(z, pos(n)) % 1.0) * 1e9)) % 10 ** 9 for n in wins})
wf = [[Fr(int(2 * x), 2) for x in v] for v in w]
def distinct_exact(z):
    phases = set()
    for n in wins:
        val = sum(z[i] * (n[0] * wf[0][i] + n[1] * wf[1][i] + n[2] * wf[2][i]) for i in range(3))
        phases.add(val - (val.numerator // val.denominator))
    return len(phases)
z_gen = np.array([np.sqrt(2), np.sqrt(3), np.sqrt(5)]) / 7
z_rat = [Fr(2, 3)] * 3                                          # the order-3 family holonomy z3 (Part C7)
step("S7.windings resolved iff an infinite-order test exists: a generic continuous family-Cartan holonomy separates all windings; the order-3 family holonomy sees a finite quotient",
     distinct_float(z_gen) == len(wins) and distinct_exact(z_rat) == 3,
     f"generic z: {distinct_float(z_gen)} of {len(wins)} distinct phases; z3 = (2/3)^3: {distinct_exact(z_rat)}")
step("S8.(statement) no TFPT process yet couples the seam time L0 to the winding lattice: dispersion omega(p) and Lorentz signature are open",
     True, "next: SPACETIME.LORENTZ.01")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "B_free": B_free, "B_finite": B_fin, "d_H": dH, "d_s": ds, "weyl": weyl}
(HERE / "process_lattice_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
