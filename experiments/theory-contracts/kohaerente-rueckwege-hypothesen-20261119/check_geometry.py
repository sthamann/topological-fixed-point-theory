"""MARKGEOMETRY20261119 (hypotheses D and 7).
D: do the three pairings of four marks give three spatial translations?
7: does a causal action on correlation operators F_w = -e_w^+ Sigma e_w see anything
   with the forms the source already has (positive ones), and what is the least extra structure?"""
from common import *
import sympy as sp
import numpy as np

c = Checks("geometry")
lam = sp.symbols("lambda")
# marks 0, 1, infinity, lambda; the three Moebius involutions realising the three pairings
P = {"0inf|1lam": sp.Matrix([[0, lam], [1, 0]]),          # z -> lam/z
     "01|inflam": sp.Matrix([[lam, -lam], [1, -lam]]),    # z -> lam (z - 1)/(z - lam)
     "0lam|1inf": sp.Matrix([[1, -lam], [1, -1]])}        # z -> (z - lam)/(z - 1)
pt = {"0": sp.Matrix([0, 1]), "1": sp.Matrix([1, 1]), "inf": sp.Matrix([1, 0]), "lam": sp.Matrix([lam, 1])}
sends = lambda M, u, w: sp.simplify((M * pt[u])[0] * pt[w][1] - (M * pt[u])[1] * pt[w][0]) == 0
c.ok("each matrix swaps its two pairs of marks (homogeneous coordinates)",
     all(sends(P["0inf|1lam"], u, w) for u, w in (("0", "inf"), ("inf", "0"), ("1", "lam"), ("lam", "1")))
     and all(sends(P["01|inflam"], u, w) for u, w in (("0", "1"), ("1", "0"), ("inf", "lam"), ("lam", "inf")))
     and all(sends(P["0lam|1inf"], u, w) for u, w in (("0", "lam"), ("lam", "0"), ("1", "inf"), ("inf", "1"))))
Ms = list(P.values())
c.ok("all three are traceless: rotations by pi about a geodesic of hyperbolic three-space", all(M.trace() == 0 for M in Ms))
c.ok("they anticommute pairwise for every cross ratio: three mutually orthogonal axes through one point",
     all(sp.simplify(Ms[i] * Ms[j] + Ms[j] * Ms[i]) == sp.zeros(2) for i in range(3) for j in range(i)))
prod = sp.simplify(Ms[0] * Ms[1])
c.ok("the product of two is proportional to the third: Klein four-group in PSL(2,C), quaternion frame in SL(2,C)",
     sp.simplify(prod[0, 0] * Ms[2][0, 1] - prod[0, 1] * Ms[2][0, 0]) == 0 and sp.simplify(prod[1, 0] * Ms[2][1, 1] - prod[1, 1] * Ms[2][1, 0]) == 0)
comm = sp.simplify(Ms[0] * Ms[1] - Ms[1] * Ms[0])
c.ok("the three generators do not commute: [A, B] = 2AB is the third one", comm != sp.zeros(2) and sp.simplify(comm - 2 * prod) == sp.zeros(2))
vec = sp.Matrix([list(M) for M in Ms])
c.ok("their complex span is sl(2,C) = so(3,1): three boosts and three rotations, no translation", vec.rank() == 3)
# canonical observer: positive Hermitian form fixed by the frame (numerical, generic and special cross ratios)
def tetrad(l):
    N = [np.array(M.subs(lam, l), dtype=complex) for M in Ms]
    N = [M / np.sqrt(np.linalg.det(M)) for M in N]
    grp = [np.eye(2), *N]
    Pm = sum(h.conj().T @ h for h in grp)
    return all(np.allclose(h.conj().T @ Pm @ h, Pm) for h in N) and np.linalg.eigvalsh(Pm).min() > 0
c.ok("a common fixed point exists (invariant positive Hermitian form) for generic, harmonic and equianharmonic marks",
     all(tetrad(l) for l in (0.3 + 0.7j, -1, np.exp(1j * np.pi / 3), 2.5)))
c.ok("verdict D: four marks give a Lorentz frame (three orthogonal axes and an observer) kinematically, for any four points; commuting translations do not arise", True)

# --- hypothesis 7 ---------------------------------------------------------------
rng = np.random.default_rng(7)
def closed_chain_spectrum(Fx, Fy):
    ev = np.linalg.eigvals(Fx @ Fy)
    return ev[np.argsort(-abs(ev))]
bad = 0
for _ in range(2000):
    n, r = 6, int(rng.integers(1, 4))
    ex, ey = (rng.normal(size=(r, n)) + 1j * rng.normal(size=(r, n)) for _ in range(2))
    S = rng.normal(size=(r, r)) + 1j * rng.normal(size=(r, r)); S = S @ S.conj().T + 0.1 * np.eye(r)   # positive spin form
    ev = closed_chain_spectrum(-ex.conj().T @ S @ ex, -ey.conj().T @ S @ ey)
    bad += np.abs(ev.imag).max() > 1e-8 or ev.real.min() < -1e-8
c.ok("positive spin form (reflection positive data): every closed chain F_x F_y has real nonnegative spectrum, 2000 random cases", bad == 0)
c.ok("then no pair of points is spacelike in the causal-fermion sense (that needs a non-real spectrum): no causal distinction", bad == 0)
bad1 = 0
for _ in range(2000):
    n = 6
    ex, ey = (rng.normal(size=(1, n)) + 1j * rng.normal(size=(1, n)) for _ in range(2))
    ev = closed_chain_spectrum(-rng.normal() * ex.conj().T @ ex, -rng.normal() * ey.conj().T @ ey)
    bad1 += np.abs(ev.imag).max() > 1e-8
c.ok("scalar evaluations (rank one, either sign): the closed chain has one real eigenvalue s_x s_y |<x|y>|^2, never a non-real pair", bad1 == 0)
# least extra structure: pair each point with its mirror, spin form = swap, signature (1,1)
from math import comb
def coh(w):                                               # degree-four sections on the sphere, <v_x, v_y> = (1 + conj(x) y)^4 / norms
    v = np.array([np.sqrt(comb(4, k)) * w**k for k in range(5)])
    return v / np.linalg.norm(v)
def F(w):
    a, b = coh(w), coh(1 / np.conj(w))                    # mirror through the equator |w| = 1
    return -(np.outer(a, b.conj()) + np.outer(b, a.conj()))
counts = {"spacelike": 0, "timelike": 0}
for _ in range(4000):
    x, y = (rng.uniform(0.05, 0.95) * np.exp(2j * np.pi * rng.uniform()) for _ in range(2))
    ev = closed_chain_spectrum(F(x), F(y))[:2]
    counts["spacelike" if abs(ev[0].imag) > 1e-9 else "timelike"] += 1
c.ok("with the mirrored pair and the swap form (signature (1,1)) both spacelike and timelike pairs occur", min(counts.values()) > 0, counts)
eq = np.linalg.matrix_rank(F(np.exp(0.3j)), tol=1e-10)
sig = np.sign(np.round(np.linalg.eigvalsh(F(0.4 + 0.2j)), 12))
c.ok("each F_w then has one positive and one negative eigenvalue; on the mirror line it degenerates to rank one",
     sorted(sig[sig != 0]) == [-1, 1] and eq == 1)
c.ok("verdict 7: with the source's positive forms the causal structure is trivial; the first nontrivial case needs the reflected copy, a chosen kernel and a chosen measure", True)
c.done({"causal_types_with_mirror_pair": counts})
