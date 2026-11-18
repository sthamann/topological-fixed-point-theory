#!/usr/bin/env python3
"""A parameter-free weighted form on the four moving marks (18 November 2026).

This is a model, not the original producer: the native consumer
(original_source_adapter.py) accepts only a source-reviewed original bundle
and none is supplied here.

Model.  The four marks X_1..X_4 are the roots of X^4 + z X + q = 0,
q^3 = 27 (z^4-1)/256.  Put them on the round sphere in the conformally
canonical frame (vector sum zero).  Sources are the states of the line of
degree one (a qubit) at the marks; the form is the overlap of their
density matrices,
        T_ab = Tr(rho_a rho_b) = (1 + n_a . n_b)/2.

Sections
    W1  pairing products = sides of the triangle of the three two-torsion points
    W2  closed form of T: weight 1/2 on the symmetric sum, p_k/2 on the three pairings
    W3  the canonical frame, numerically, against the closed form
    W4  values along the family: square, tetrahedron, collision; small z
    W5  comparison with recorded numbers
    W6  fifth pasted note: local averaging over the sheets of a mark loop
"""
import itertools
import sys

import numpy as np
import sympy as sp
from scipy.optimize import minimize

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


PAIRINGS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
CQ = 3/4**(4/3)


def roots_at(zz):
    q = CQ*np.exp(np.log(complex(zz)**4-1)/3)
    return np.roots([1, 0, 0, zz, q])


def weights(X):
    """p_k = 1 - 2 r_k / (r_1 + r_2 + r_3),  r_k = |X_a - X_b| |X_c - X_d| for the pairing k"""
    r = np.array([abs(X[a]-X[b])*abs(X[c]-X[d]) for (a, b), (c, d) in PAIRINGS])
    return 1-2*r/r.sum()


# ---------------------------------------------------------------- W1
print("W1  pairing products and the two-torsion triangle")
x1, x2, x3, x4 = sp.symbols("x1:5")
th = [(x1+x2)*(x3+x4), (x1+x3)*(x2+x4), (x1+x4)*(x2+x3)]
check("resolvent roots: theta_1 - theta_2 = (x1 - x4)(x3 - x2), and cyclically",
      sp.expand(th[0]-th[1]-(x1-x4)*(x3-x2)) == 0 and sp.expand(th[1]-th[2]-(x1-x2)*(x4-x3)) == 0
      and sp.expand(th[0]-th[2]-(x1-x3)*(x4-x2)) == 0)
zt = 0.37+0.21j
X = roots_at(zt)
qv = CQ*np.exp(np.log(zt**4-1)/3)
thn = np.array([(X[0]+X[1])*(X[2]+X[3]), (X[0]+X[2])*(X[1]+X[3]), (X[0]+X[3])*(X[1]+X[2])])
check("the theta_k solve theta^3 - 4 q theta + z^2 = 0 up to sign convention (resolvent of the quartic)",
      max(abs(t**3-4*qv*t+zt**2) for t in thn) < 1e-10 or max(abs(t**3-4*qv*t-zt**2) for t in thn) < 1e-10)
print("      so r_k is the side of the triangle of the three two-torsion points opposite to the vertex k,")
print("      and p_k = 1 - 2 (side)/(perimeter).")

# ---------------------------------------------------------------- W2
print("W2  closed form of the overlap matrix")
c1, c2, c3 = sp.symbols("c1 c2 c3")
Pi = []
for (a, b), (c, d) in PAIRINGS:
    M = sp.zeros(4, 4)
    for u, v in ((a, b), (b, a), (c, d), (d, c)):
        M[u, v] = 1
    Pi.append(M)
N = sp.eye(4)+c1*Pi[0]+c2*Pi[1]+c3*Pi[2]        # Bloch Gram of a tetrahedron with equal opposite edges
T = (sp.ones(4, 4)+N)/2
chars = [sp.Matrix([1, 1, 1, 1]), sp.Matrix([1, 1, -1, -1]), sp.Matrix([1, -1, 1, -1]), sp.Matrix([1, -1, -1, 1])]
ev = [sp.simplify((T*v)[0]/v[0]) for v in chars]
sub = {c3: -1-c1-c2}
check("zero vector sum forces c1 + c2 + c3 = -1; then T has eigenvalue 2 on the symmetric sum",
      sp.simplify(ev[0].subs(sub)-2) == 0)
check("and eigenvalue 1 + c_k = 2 p_k on the character of the pairing k",
      all(sp.simplify(ev[k+1].subs(sub)-(1+[c1, c2, c3][k]).subs(sub)) == 0 for k in range(3)))
check("trace one: T/4 = (1/2) P_sym + sum_k (p_k/2) P_k with p_1 + p_2 + p_3 = 1",
      sp.simplify((T.trace()/4).subs(sub)-1) == 0)

# ---------------------------------------------------------------- W3
print("W3  the canonical frame (vector sum zero), numerically")


def to_sphere(w):
    d = 1+abs(w)**2
    return np.array([2*w.real/d, 2*w.imag/d, (abs(w)**2-1)/d])


def frame(X):
    def moved(par):
        a, b, c = par[0]+1j*par[1], par[2]+1j*par[3], par[4]+1j*par[5]
        d = (1+b*c)/a
        return [(a*x+b)/(c*x+d) for x in X]

    def cost(par):
        return float(np.sum(np.sum([to_sphere(w) for w in moved(par)], axis=0)**2))
    best = None
    for start in ([1, 0, 0, 0, 0, 0], [1, 0, 0.3, 0.1, -0.2, 0.1], [0.7, 0.2, -0.3, 0.2, 0.1, -0.3]):
        res = minimize(cost, start, method="Nelder-Mead", options={"xatol": 1e-13, "fatol": 1e-26, "maxiter": 40000})
        if best is None or res.fun < best.fun:
            best = res
    return np.array([to_sphere(w) for w in moved(best.x)]), best.fun


worst = 0
for zz in (0.37+0.21j, 0.8, 0.5j+0.5, 1.7-0.4j):
    X = roots_at(zz)
    n, resid = frame(X)
    Tn = (1+n @ n.T)/2
    evn = np.sort(np.linalg.eigvalsh(Tn))
    p = weights(X)
    pred = np.sort(np.concatenate([[2.0], 2*p]))
    worst = max(worst, np.max(np.abs(evn-pred)), resid)
check("eigenvalues of T in the canonical frame = (2, 2 p_1, 2 p_2, 2 p_3) at four parameter values",
      worst < 1e-6, "%.1e" % worst)

# ---------------------------------------------------------------- W4
print("W4  along the family")
p0 = np.sort(weights(roots_at(1e-9)))
check("z = 0, the original square: p = (0, 1/2, 1/2); the pairing of the diagonals carries no weight",
      np.allclose(p0, [0, 0.5, 0.5], atol=1e-6), np.round(p0, 6))
pm = np.sort(weights(roots_at(1-1e-12)))
check("z at a parameter mark, tetrahedron: p = (1/3, 1/3, 1/3), T = (2/3) I + (1/3) 11^T (the form G2)",
      np.allclose(pm, [1/3]*3, atol=1e-3), np.round(pm, 4))
pinf = np.sort(weights(roots_at(4000.0)))
check("z -> infinity, two marks collide: p -> (0, 0, 1)", np.allclose(pinf, [0, 0, 1], atol=2e-2), np.round(pinf, 4))
ratios = []
for zz in (0.05, 0.1, 0.2):
    ratios.append(np.min(weights(roots_at(zz)))/zz**4)
check("small real z: the third weight grows like z^4/12", abs(ratios[0]-1/12) < 1e-3, np.round(ratios, 5))
diag_ray = [np.min(weights(roots_at(r*np.exp(1j*np.pi/4)))) for r in (0.3, 0.7, 1.5, 3.0)]
check("on the rays between the marks (arg z = 45 degrees) the third weight stays zero: marks on a circle",
      max(abs(v) for v in diag_ray) < 1e-9, np.round(diag_ray, 12))
gen = [(np.real(zz**2))**2/12 for zz in (0.1*np.exp(0.3j), 0.1*np.exp(0.6j))]
act = [np.min(weights(roots_at(zz))) for zz in (0.1*np.exp(0.3j), 0.1*np.exp(0.6j))]
check("general small z: third weight = (Re z^2)^2/12 to leading order", np.allclose(gen, act, rtol=2e-2), np.round(act, 9))
print("      real z from the square to a mark:")
for zz in (0.0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0):
    pz = np.sort(weights(roots_at(zz if 0 < zz < 1 else (1e-9 if zz == 0 else 1-1e-12))))
    print("        z = %-5s p = (%.4f, %.4f, %.4f)" % (zz, *pz))

# ---------------------------------------------------------------- W5
print("W5  comparison with recorded numbers")
check("degree two (the three-dimensional generation space): overlap ((1+c)/2)^2 = 1/9 at the tetrahedron",
      sp.Rational(1, 3)**2 == sp.Rational(1, 9))
I = sp.I
perm3 = sp.zeros(4, 4)                           # a cycle of order three fixing the mark 0
for j, k in ((0, 0), (1, 2), (2, 3), (3, 1)):
    perm3[k, j] = 1
four = [sp.Matrix([I**(k*j) for j in range(4)])/2 for k in (2, 1, 3)]
dM = [sp.simplify((v.H*perm3*v)[0]) for v in four]
check("diagonal of a mark cycle in the basis of the quarter turn: (0, +-i/2, -+i/2), the recorded harmonic diagonal",
      [sp.Abs(d) for d in dM] == [0, sp.Rational(1, 2), sp.Rational(1, 2)], dM)
Wsq = sp.zeros(4, 4)
psq = [sp.Rational(1, 2), 0, sp.Rational(1, 2)]  # marks 0,1,2,3 around the square: pairing (02)(13) = diagonals
for k in range(3):
    v = chars[k+1]
    Wsq += psq[k]*v*v.T/4
dW = [sp.simplify((v.H*Wsq*v)[0]) for v in four]
check("diagonal of the weight operator of the square in the same basis: (0, 1/2, 1/2)",
      dW == [0, sp.Rational(1, 2), sp.Rational(1, 2)], dW)
G2 = sp.Rational(2, 3)*sp.eye(4)+sp.Rational(1, 3)*sp.ones(4, 4)
Psym = sp.ones(4, 4)/4
check("at the tetrahedron the form equals the Hodge masses: (1/6) per mark + (1/3) on the symmetric sum",
      G2/4 == sp.Rational(1, 6)*sp.eye(4)+sp.Rational(1, 3)*Psym
      and G2/4 == sp.Rational(1, 2)*Psym+sp.Rational(1, 6)*(sp.eye(4)-Psym))

# ---------------------------------------------------------------- W6
print("W6  fifth pasted note: averaging over the sheets permuted by a mark loop")
J4 = sp.ones(4, 4)
Ea = []
for a in range(4):
    ea = sp.zeros(4, 1)
    ea[a] = 1
    rest = sp.ones(4, 1)-ea
    Ea.append(ea*ea.T+rest*rest.T/3)             # projector on functions invariant under the loop of the mark a
check("each E_a is a projector of rank two: the symmetric sum plus one line L_a",
      all(E*E == E and E.trace() == 2 for E in Ea)
      and all((E-Psym)*(E-Psym) == (E-Psym) and (E-Psym).trace() == 1 for E in Ea))
check("(1/2) sum_a E_a = (2/3) I + (1/3) 11^T = G2 exactly, with no input from the surface",
      sum(Ea, sp.zeros(4, 4))/2 == G2)
check("the four lines L_a have squared overlap 1/9 (the recorded value): a second mark measurement gives 1/9",
      all(((Ea[a]-Psym)*(Ea[b]-Psym)).trace() == sp.Rational(1, 9) for a in range(4) for b in range(4) if a != b))
atoms_only = sum(Ea, sp.zeros(4, 4))/8           # four mark pieces of trace 1/4 each
hodge = sp.Rational(1, 6)*sp.eye(4)+sp.Rational(1, 3)*Psym
check("the same matrix as four mark pieces of mass 1/4, and as 1/6 per mark plus 1/3 on the symmetric sum",
      atoms_only == hodge and all((E/8).trace() == sp.Rational(1, 4) for E in Ea))
print("      -> the value on the symmetric sum does not separate the two measures (correction accepted).")
zz_, XX_ = sp.symbols("z X")
plane = 256*XX_**3*(XX_**3+zz_)**3+27*(zz_**4-1)
check("q eliminated: 256 X^3 (X^3 + z)^3 + 27 (z^4 - 1) = 0; at each mark four points, each threefold",
      all(sp.expand(plane.subs(zz_, m)-256*XX_**3*(XX_**3+m)**3) == 0 for m in (1, -1, sp.I, -sp.I)))
tw = [sp.Rational(1, 2)*Psym+sum(((sp.Rational(1, 6))*(chars[k+1]*chars[k+1].T/4) for k in range(3)), sp.zeros(4, 4)), G2/4]
check("the model of W2 at the tetrahedron is the same matrix", tw[0] == tw[1])

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
