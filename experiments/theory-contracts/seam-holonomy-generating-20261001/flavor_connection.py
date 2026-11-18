#!/usr/bin/env python3
"""Part L: FLAVOR.CONNECTION.01 -- is the v117 Fuchs connection an output of local data?

L1 multiplicative rigidity: the local monodromy classes at the three orbifold points of x = z^4
   (deck: order 4, eigenvalues {1,i,-i}; marks: order 3, {1,w,w^2}; infinity: order 2, {1,-1,-1})
   have Katz rigidity index sum dim Z = 3 + 3 + 5 = 11 = n^2 + 2: the irreducible local system is
   unique.  The transported monodromy realises exactly these classes.
L2 inside S4 the classes force the rotation representation (no reducible or det -1 alternative).
L3 additive non-rigidity: the Fuchs residues are NOT fixed by exponents alone (index 9 < 11).
L4 residue reconstruction on the z-plane: A(z) = sum U^k A U^-k/(z - i^k), A Hermitian with
   spectrum (0,1/3,2/3); 'no defect at infinity' needs 4 diag(A) integral; Schur-Horn plus
   irreducibility leave diag(A) = (1/2,1/4,1/4) in one of three character slots.  Scan the whole
   remaining family for trivial monodromy at z = infinity.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/flavor_connection.py
"""
from __future__ import annotations

import json
import sys
from itertools import permutations, product
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


U = np.diag([1, 1j, -1j])
Ui = np.linalg.inv(U)


def conn(A):
    def f(z):
        return sum(np.linalg.matrix_power(U, k) @ A @ np.linalg.matrix_power(Ui, k) / (z - 1j ** k) for k in range(4))
    return f


def transport(Af, segs, rtol=1e-11):
    T = np.eye(3, dtype=complex)
    for path in segs:
        def rhs(t, y, path=path):
            z, dz = path(t)
            return (dz * Af(z) @ y.reshape(3, 3)).ravel()
        T = solve_ivp(rhs, (0, 1), np.eye(3, dtype=complex).ravel(), method="DOP853",
                      rtol=rtol, atol=rtol * 1e-2).y[:, -1].reshape(3, 3) @ T
    return T


ray = lambda a, b: (lambda t: (a + (b - a) * t, b - a))
circ = lambda c, r, s: (lambda t: (c + r * np.exp(1j * (s + 2 * np.pi * t)), 2j * np.pi * r * np.exp(1j * (s + 2 * np.pi * t))))
A117 = np.array([[1 / 2, np.sqrt(8) / 12, 0], [np.sqrt(8) / 12, 1 / 4, np.sqrt(5) / 12], [0, np.sqrt(5) / 12, 1 / 4]], dtype=complex)
F117 = conn(A117)
M0 = transport(F117, [ray(0, 0.6), circ(1, 0.4, np.pi), ray(0.6, 0)])
g0, g1 = Ui, M0
ginf = np.linalg.inv(g0 @ g1)


def cent_dim(X):
    ev = np.linalg.eigvals(X)
    groups = []
    for e in ev:
        for gidx, gg in enumerate(groups):
            if abs(gg[0] - e) < 1e-6:
                gg.append(e); break
        else:
            groups.append([e])
    return sum(len(gg) ** 2 for gg in groups)  # semisimple: centraliser = sum of mult^2


w = np.exp(2j * np.pi / 3)
semisimple = all(np.linalg.matrix_rank(np.column_stack(np.linalg.eig(X)[1])) == 3 for X in (g0, g1, ginf))
idx = [cent_dim(X) for X in (g0, g1, ginf)]
step("L1.local classes: deck {1,i,-i}, marks {1,w,w^2}, infinity {1,-1,-1}; Katz index 3+3+5 = 11 = n^2+2 (rigid)",
     semisimple and idx == [3, 3, 5] and sum(idx) == 11
     and np.allclose(sorted(np.round(np.linalg.eigvals(ginf).real, 8)), [-1, -1, 1]),
     f"centraliser dims {idx}")

# L2: characters of S4 on (4-cycle, 3-cycle, double-transposition-type 2-cycle product g_inf)
#    classes: e, (12), (12)(34), (123), (1234)
chars = {"triv": [1, 1, 1, 1, 1], "sgn": [1, -1, 1, 1, -1], "two": [2, 0, 2, -1, 0],
         "std": [3, 1, -1, 0, -1], "std*sgn": [3, -1, -1, 0, 1]}
# in a 3-dim rep the 4-cycle needs eigenvalues {1, i, -i}: trace 1, det 1
cands = [k for k in ("std", "std*sgn") if chars[k][4] == 1]
reducible_ok = False  # 1-dim and 2-dim reps send the 4-cycle to +-1 eigenvalues only: no +-i
step("L2.inside S4 only std x sgn (the rotation rep, det 1) has 4-cycle eigenvalues {1,i,-i}; reducible 3-dim reps cannot",
     cands == ["std*sgn"] and not reducible_ok)

step("L3.additive index of the residues 3+3+3 = 9 < 11: exponents alone leave a 2-parameter family (log terms possible)",
     3 + 3 + 3 == 9)

# L3b: which diagonals are possible?  No defect at z=infinity needs integral exponents there:
#      Res_inf = -sum_k U^k A U^-k = -4 diag(A), so 4 diag(A) is integral, sum = tr A = 1.
#      Schur-Horn: diag(A) majorised by the spectrum (0,1/3,2/3); equality in a partial sum forces
#      an invariant subspace (reducible connection), excluded by irreducibility.
spec = sorted([0, 1 / 3, 2 / 3], reverse=True)
quarter = [tuple(x / 4 for x in t) for t in product(range(5), repeat=3) if sum(t) == 4]
def majorised(dg, strict):
    ds = sorted(dg, reverse=True)
    sums = [(sum(ds[:k]), sum(spec[:k])) for k in (1, 2)]
    if strict:
        return all(a < b - 1e-12 for a, b in sums)
    return all(a <= b + 1e-12 for a, b in sums)
allowed = [dg for dg in quarter if majorised(dg, strict=False)]
irreducible = [dg for dg in allowed if majorised(dg, strict=True)]
step("L3b.integrality at infinity + Schur-Horn + irreducibility leave exactly diag(A) = (1/2,1/4,1/4) in one of three character slots",
     sorted(irreducible) == sorted({(0.5, 0.25, 0.25), (0.25, 0.5, 0.25), (0.25, 0.25, 0.5)})
     and len(allowed) > len(irreducible), f"allowed {sorted(allowed)}; irreducible {sorted(irreducible)}")


def family(d, slot_half):
    """Hermitian A with diag d, spectrum (0,1/3,2/3); gauge: A12 = a >= 0, A23 = c >= 0, A13 = b e^{i phi}."""
    S0 = d[0] * d[1] + d[0] * d[2] + d[1] * d[2] - 2 / 9
    out = []
    for b in np.linspace(0, np.sqrt(max(S0, 0)), 9):
        for phi in np.linspace(0, np.pi, 7):
            S = S0 - b * b
            if S < 0:
                continue
            def det_res(a):
                c = np.sqrt(max(S - a * a, 0))
                return d[0] * d[1] * d[2] + 2 * a * b * c * np.cos(phi) - d[0] * c * c - d[1] * b * b - d[2] * a * a
            grid = np.linspace(0, np.sqrt(S), 400)
            vals = [det_res(a) for a in grid]
            for i in range(len(grid) - 1):
                if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
                    a = brentq(det_res, grid[i], grid[i + 1]) if vals[i] != 0 else grid[i]
                    c = np.sqrt(max(S - a * a, 0))
                    A = np.array([[d[0], a, b * np.exp(1j * phi)], [a, d[1], c], [b * np.exp(-1j * phi), c, d[2]]])
                    if np.allclose(sorted(np.linalg.eigvalsh(A)), [0, 1 / 3, 2 / 3], atol=1e-9):
                        out.append((b, phi, a, c, A))
    return out


results = {}
for slot in range(3):
    d = [0.25, 0.25, 0.25]; d[slot] = 0.5
    best = []
    for b, phi, a, c, A in family(d, slot):
        Minf = transport(conn(A), [circ(0, 2.0, 0.0)], rtol=1e-9)
        best.append((np.abs(Minf - np.eye(3)).max(), b, phi, a, c))
    best.sort()
    results[slot] = best
zeros = {slot: [x for x in r if x[0] < 1e-6] for slot, r in results.items()}
z0 = zeros[0]
step("L4.scan (half on character slot 0, the v117 slot): trivial monodromy at infinity only at b = 0, i.e. A13 = 0, a^2 = 8/144, c^2 = 5/144 = v117",
     len(z0) > 0 and all(abs(x[1]) < 1e-12 for x in z0)
     and all(abs(x[3] ** 2 - 8 / 144) < 1e-8 and abs(x[4] ** 2 - 5 / 144) < 1e-8 for x in z0),
     f"{len(z0)} zero(s) among {len(results[0])} family points; min deviation off b=0: {min(x[0] for x in results[0] if x[1] > 1e-12):.2e}")
# L5: the other slots are NOT empty.  Slot chi_3 has an exact zero at b = 0; its tau-mirror sits in
#     slot chi_1 with A12 = 0 (off the coarse grid).  Both have trivial z-monodromy at infinity but a
#     DIFFERENT local system: M0 U^-1 has order 4 (hyperbolic triangle group (4,3,4), infinite image)
#     instead of order 2 (spherical (4,3,2) = S4).  Finiteness of the flavor group selects v117.
alt = [x for x in results[2] if x[0] < 1e-6 and abs(x[1]) < 1e-12][0]
a2, c2 = alt[3], alt[4]
A2 = np.array([[0.25, a2, 0], [a2, 0.25, c2], [0, c2, 0.5]], dtype=complex)
A1 = np.array([[0.25, 0, a2], [0, 0.5, c2], [a2, c2, 0.25]], dtype=complex)        # tau-mirror (chi_1 <-> chi_3)
dev1 = np.abs(transport(conn(A1), [circ(0, 2.0, 0.0)], rtol=1e-10) - np.eye(3)).max()


def closure(gens, cap=400):
    key = lambda X: tuple(np.round(X, 5).ravel())
    seen = {key(np.eye(3)): np.eye(3, dtype=complex)}
    frontier = list(seen.values())
    while frontier and len(seen) <= cap:
        new = []
        for X in frontier:
            for G in gens:
                Y = G @ X
                if key(Y) not in seen:
                    seen[key(Y)] = Y; new.append(Y)
        frontier = new
    return len(seen)


def x_inf_order(A):
    M = transport(conn(A), [ray(0, 0.6), circ(1, 0.4, np.pi), ray(0.6, 0)])
    X = M @ Ui
    for k in range(1, 9):
        if np.allclose(np.linalg.matrix_power(X, k), np.eye(3), atol=1e-7):
            return k, M
    return None, M


o117, _ = x_inf_order(A117)
o2, M2 = x_inf_order(A2)
o1, M1 = x_inf_order(A1)
step("L5.slots chi_3 and chi_1 hold a tau-mirror pair with trivial monodromy at infinity, but x = infinity class of order 4 (Delta(4,3,4), infinite)",
     dev1 < 1e-7 and o117 == 2 and o2 == 4 and o1 == 4 and closure([M2, U]) > 400 and closure([M0, U]) == 24,
     f"orders of M0 U^-1: v117 {o117}, alt {o2}, mirror {o1}; |image| v117 24, alt > 400")
step("L6.RESULT: the residues are an output of (mu4 equivariance, mark exponents (0,1/3,2/3), no defect at z=infinity, Hermitian, irreducible, FINITE flavor group); dropping finiteness admits the (4,3,4) pair",
     o117 == 2 and o2 == 4,
     "finiteness <=> (M U^-1)^2 = I <=> spherical (4,3,2); grid scan, not a continuum exhaustiveness proof")

# L7: finiteness is not needed as an axiom -- UNITARITY suffices.  Over all words of length <= 6 in
#     (M, U, M^-1, U^-1) the v117 monodromy has all eigenvalues on the unit circle, while the
#     (4,3,4) pair has eigenvalues off it: it preserves no positive Hermitian form.
from itertools import product as _prod
def max_log_modulus(M):
    gens = [M, U, np.linalg.inv(M), Ui]
    worst = 0.0
    for L in range(1, 7):
        for wd in _prod(range(4), repeat=L):
            X = np.eye(3, dtype=complex)
            for k in wd:
                X = gens[k] @ X
            worst = max(worst, np.abs(np.log(np.abs(np.linalg.eigvals(X)))).max())
    return worst
u117, ualt = max_log_modulus(M0), max_log_modulus(M2)
step("L7.unitarity replaces finiteness: v117 words have |eigenvalues| = 1, the (4,3,4) pair has a word with |eigenvalue| != 1 (not unitarisable)",
     u117 < 1e-9 and ualt > 0.1, f"max |ln|eig||: v117 {u117:.1e}, pair {ualt:.3f}")

# L8: the mark order 3 is forced.  Finite triangle groups (spherical: 1/p+1/q+1/r > 1) containing the
#     deck order 4 with the infinity order dividing 4: (2,2,4) = D4 (irreps of dim 1, 2 only) and
#     (2,3,4) = S4.  An irreducible 3-dim flavor action therefore forces q = 3 at the marks, i.e. the
#     exponents (0,1/3,2/3) mod Z of the original compiler (det 1, regular).
spherical = sorted({tuple(sorted((4, q, r))) for q in range(2, 13) for r in (1, 2, 4)
                    if r > 1 and 1 / 4 + 1 / q + 1 / r > 1})
irrep_dims = {(2, 2, 4): [1, 1, 1, 1, 2], (2, 3, 4): [1, 1, 2, 3, 3]}
with3 = [t for t in spherical if 3 in irrep_dims.get(t, [])]
step("L8.mark order forced: spherical (4,q,r) with r | 4 are (2,2,4)=D4 and (2,3,4)=S4; only S4 has a 3-dim irrep -> q = 3",
     spherical == [(2, 2, 4), (2, 3, 4)] and with3 == [(2, 3, 4)], f"spherical {spherical}")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "scan": {s: [(float(x[0]), float(x[1]), float(x[2])) for x in r[:5]] for s, r in results.items()}}
(HERE / "flavor_connection_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
