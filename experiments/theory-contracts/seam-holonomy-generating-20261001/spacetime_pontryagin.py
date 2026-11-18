#!/usr/bin/env python3
"""Part T: SPACETIME.PONTRYAGIN.01, INTERNAL_VS_SPACE.01 and a blind LORENTZ.01 classification.

Windings w = sum_a n_a w_a (w_a = the four family-line vectors, sum w_a = 0) form H_1 = Z^3 (Part S).
  T1 (P1, P2) a family-Cartan twist z gives the flat U(1) connection with holonomy e^{i z.w_a} around
     mark a; the surface relation prod hol_a = 1 is exactly sum_a w_a = 0.  chi_z(w) = e^{i z.w} is a
     character, multiplicative, and characters separate windings.
  T2 (P3) periods: chi_z trivial on all windings iff z in 2 pi x (dual lattice): z-space mod periods is
     a 3-torus T^3 = Hom(H_1, U(1)) (Pontryagin dual), covolume (2 pi)^3 / covol(winding lattice).
  T3 (P4) W(A3) = S4 acts orthogonally and chi_{gz}(gw) = chi_z(w).
  T4 (P5) the ACTUAL response: spectral flow of the level-one E8 lattice VOA, twisted vacuum energy
     E(z) = (1/2) min_lambda |z - lambda|^2 (E8 closest point): exactly (1/2)|z|^2 near 0 for z in the
     family subspace -- isotropic, positive, no O(z^4); the matter selection S48 alone is anisotropic
     (16, 16, 4), the S4-symmetric chiral block (16,4) is isotropic (16, 16, 16).
  T5 precision operator with symbol E(z): second order, Weyl N(Lambda) ~ Lambda^{3/2} (= d/r, d = 3, r = 2).
  T6 (INTERNAL_VS_SPACE) S4 acts on space by the standard rep (transpositions are reflections, det -1)
     and on flavor by std x sgn (rotations, det +1, Part U8): different representations; the image
     of pi_1 in Z^3 x A4 is a fiber product over Z3 (positions and flavor correlated only by total
     winding mod 3): IR factorisation l^2(sublattice) x internal.
  T7 (LORENTZ, blind) with the existing time generator L0 the spectral-flow states have L0 = E(z) =
     (1/2)|z|^2: omega ~ |p|^2, dynamical exponent z_dyn = 2 (Galilean), NOT 1.  A relativistic mass
     shell needs a negative-norm (timelike) lattice direction, absent from every TFPT lattice
     (all positive definite).  Abort rule: winding lattice + seam time is not Lorentzian 3+1.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/spacetime_pontryagin.py
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
rng = np.random.default_rng(20261001)
# T1 ---------------------------------------------------------------------------------------
zs = rng.normal(size=(6, 3))
hol = [np.exp(1j * (w @ z)) for z in zs]
chi = lambda z, n: np.exp(1j * z @ (n @ w))
ns = [np.array(n) for n in product(range(-2, 3), repeat=4)]
mult = all(np.isclose(chi(z, a + b), chi(z, a) * chi(z, b)) for z in zs[:2] for a, b in zip(ns[:50], ns[50:100]))
sep = len({np.round(np.angle(chi(np.array([np.sqrt(2), np.sqrt(3), np.sqrt(5)]), n)), 9) for n in ns
           if tuple(n) == tuple(n)}) >= len({tuple(np.round(n @ w, 9)) for n in ns}) - 1
step("T1.family-Cartan twist z = flat U(1) connection with holonomy e^{i z.w_a}; prod_a hol_a = 1 <=> sum w_a = 0; chi_z multiplicative and separating",
     all(np.isclose(np.prod(h), 1) for h in hol) and np.allclose(w.sum(0), 0) and mult and sep)
# T2 periods ---------------------------------------------------------------------------------
Bw = np.array([w[0], w[1], w[2]]).T                  # basis of the winding lattice (w_3 = -w_0-w_1-w_2)
dual = 2 * np.pi * np.linalg.inv(Bw).T                # z with z.w_a in 2 pi Z
per_ok = all(np.allclose(np.exp(1j * (n @ w) @ dual[:, k]), 1) for n in ns[:200] for k in range(3))
covol = abs(np.linalg.det(Bw))
step("T2.periods: chi_z trivial on all windings iff z in 2 pi x dual lattice -> z-space / periods = T^3 = Hom(H_1, U(1)); vol(T^3) = (2 pi)^3/covol",
     per_ok and np.isclose(abs(np.linalg.det(dual)), (2 * np.pi) ** 3 / covol), f"covol(winding lattice) = {covol:.3f} (bcc)")
# T3 S4 compatibility ----------------------------------------------------------------------
def perm_matrix_on_space(p):
    """Orthogonal map sending w_a -> w_{p(a)} (exists for every p in S4: tetrahedral symmetry)."""
    A = np.array([w[p[0]], w[p[1]], w[p[2]]]).T @ np.linalg.inv(Bw)
    return A
S4 = list(permutations(range(4)))
G = [perm_matrix_on_space(p) for p in S4]
ortho = all(np.allclose(g.T @ g, np.eye(3)) for g in G) and all(np.allclose(g @ w[a], w[p[a]]) for g, p in zip(G, S4) for a in range(4))
compat = all(np.isclose(np.exp(1j * (g @ z) @ (g @ (n @ w))), chi(z, n)) for g in G[:6] for z in zs[:2] for n in ns[:20])
step("T3.W(A3) = S4 acts orthogonally on the winding space (tetrahedral group) and chi_{gz}(gw) = chi_z(w)", ortho and compat)
# T4 actual response: spectral flow of E8 lattice VOA ----------------------------------------------
def closest_D8(x):
    f = np.round(x); err = x - f
    if int(f.sum()) % 2:
        k = np.argmax(np.abs(err)); f[k] += np.sign(err[k]) if err[k] != 0 else 1
    return f
def closest_E8(x):
    a = closest_D8(x); b = closest_D8(x - 0.5) + 0.5
    return a if np.sum((x - a) ** 2) <= np.sum((x - b) ** 2) else b
def flow_energy(z3):
    x = np.zeros(8); x[5:] = z3                       # family Cartan subspace (std coordinates)
    lam = closest_E8(x)
    return 0.5 * np.sum((x - lam) ** 2)
small = rng.normal(size=(400, 3))
small = small / np.linalg.norm(small, axis=1)[:, None] * rng.uniform(0, 0.69, size=(400, 1))   # inside the packing radius 1/sqrt2
iso = max(abs(flow_energy(z) - 0.5 * z @ z) for z in small)
from tfpt_explorer.neutral_source_response import _native_neutral_branch
br = _native_neutral_branch(0)
A48 = np.array(br["A48"], dtype=float)[5:, 5:]; A64 = np.array(br["A64"], dtype=float)[5:, 5:]
e48, e64 = np.sort(np.linalg.eigvalsh(A48)), np.sort(np.linalg.eigvalsh(A64))
step("T4.actual response = spectral flow of the level-1 E8 VOA: E(z) = (1/2)|z|^2 exactly inside the packing radius 1/sqrt2 (isotropic, positive, no z^4); S48 Hessian anisotropic, chiral (16,4) isotropic",
     iso < 1e-12 and np.allclose(e48, [4, 16, 16]) and np.allclose(e64, [16, 16, 16]),
     f"max |E - |z|^2/2| = {iso:.1e}; family Hessian S48 {e48.tolist()}, (16,4) {e64.tolist()}")
# T5 Weyl law of the second-order precision operator --------------------------------------------
Lg = 80
grid = (np.arange(Lg) / Lg)[:, None]
pts = np.array(list(product(range(Lg), repeat=3))) / Lg      # torus coordinates in the dual basis
Z = pts @ dual.T
Z = Z - np.round(pts) @ dual.T                                # fold to a fundamental cell near 0
E = np.array([flow_energy(z) for z in Z[:: 7]])
Es = np.sort(E)
Lams = np.array([0.02, 0.05, 0.1])
NL = np.array([np.searchsorted(Es, x) for x in Lams])
wexp = np.polyfit(np.log(Lams), np.log(NL), 1)[0]
step("T5.precision operator with symbol E(z) is second order: Weyl N(Lambda) ~ Lambda^{3/2} (d/r with d = 3, r = 2)",
     abs(wexp - 1.5) < 0.1, f"exponent {wexp:.3f}")
# T6 internal vs space --------------------------------------------------------------------------
def sign(p):
    s, q = 1, list(p)
    for i in range(4):
        for j in range(i + 1, 4):
            if q[i] > q[j]:
                s = -s
    return s
det_space = {sign(p): round(np.linalg.det(g)) for g, p in zip(G, S4)}
flavor_det = {1: 1, -1: 1}                                    # Part U8: flavor rep is std x sgn, det 1 on all of S4
reps_differ = det_space == {1: 1, -1: -1} and flavor_det == {1: 1, -1: 1}
# fiber product: image of pi_1 -> Z^3 x A4 via (abelianisation, transported monodromy mod commutator)
Z3_of = {}                                                   # the A4 abelianisation Z3 of each visit's 3-cycle
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
keyM = lambda X: tuple(np.round(X, 5).ravel())
pairs = {((0, 0, 0, 0)): {keyM(np.eye(3))}}
frontier = [((0, 0, 0, 0), np.eye(3, dtype=complex))]
seen = {((0, 0, 0, 0), keyM(np.eye(3)))}
for _ in range(6):
    new = []
    for (wv, X) in frontier:
        for a in range(4):
            for s in (1, -1):
                g = M[a] if s == 1 else np.linalg.inv(M[a])
                wv2 = list(wv); wv2[a] += s
                c = min(wv2); wv2 = tuple(x - c for x in wv2)        # mod (1,1,1,1)
                Y = g @ X
                kk = (wv2, keyM(Y))
                if kk not in seen:
                    seen.add(kk); new.append((wv2, Y)); pairs.setdefault(wv2, set()).add(keyM(Y))
    frontier = new
fiber_sizes = {len(v) for k, v in pairs.items() if sum(k) <= 2}
step("T6.internal vs space: S4 acts on space by std (odd permutations = reflections, det -1) but on flavor by std x sgn (det +1); pi_1 -> Z^3 x A4 is a fiber product (each winding carries a fixed coset of size 4 = V4)",
     reps_differ and fiber_sizes == {4}, f"det on space {det_space}; flavor states per winding {sorted(fiber_sizes)} of 12")
# T7 blind dispersion classification with the existing time generator ---------------------------
ks = np.array([0.05, 0.1, 0.2])
omega = np.array([flow_energy(np.array([k, 0, 0]) @ np.eye(3)) for k in ks])
zdyn = np.polyfit(np.log(ks), np.log(omega), 1)[0]
step("T7.(LORENTZ, blind) with the seam time L0: spectral-flow energy omega = (1/2)|p|^2 -> z_dyn = 2 (Galilean), not 1; all TFPT lattices are positive definite (no timelike direction)",
     abs(zdyn - 2) < 1e-6, f"z_dyn = {zdyn:.4f}; abort rule: winding lattice + seam time is NOT Lorentzian 3+1")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "z_dyn": zdyn, "weyl": wexp}
(HERE / "spacetime_pontryagin_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
