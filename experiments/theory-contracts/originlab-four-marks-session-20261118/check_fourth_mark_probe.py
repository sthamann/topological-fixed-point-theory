#!/usr/bin/env python3
"""Thirteenth pasted note: the probe at the fourth mark, and the single fermion at a mark (18 November 2026).

Sections
    P1  stress insertion at -i for a three-point word at 1, i, -1 (conformal Ward identity)
    P2  the same word at (1, i, -1) and at (1, -1, -i): equal modulus one-sided, doubled when paired
    P3  one fermion at a mark through the recorded reflection: norm square and mutual product
    P4  the tetrahedral torus has no mirror with two fixed circles
"""
import itertools
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 30
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- P1
print("P1  stress insertion at the fourth mark")


def three_point(zs, hs):
    """holomorphic three-point function of primaries with weights hs"""
    z1, z2, z3 = zs
    h1, h2, h3 = hs
    return (z1-z2)**(-(h1+h2-h3))*(z2-z3)**(-(h2+h3-h1))*(z1-z3)**(-(h1+h3-h2))


def ward_ratio(zs, hs, z4):
    """<T(z4) phi phi phi>/<phi phi phi> = sum_i [h_i/(z4-z_i)^2 + (1/(z4-z_i)) d_i log F]"""
    zsym = sp.symbols("z1:4")
    F = three_point(zsym, hs)
    tot = 0
    for i in range(3):
        tot += hs[i]/(z4-zsym[i])**2+sp.diff(sp.log(F), zsym[i])/(z4-zsym[i])
    return sp.simplify(tot.subs(dict(zip(zsym, zs))))


pts = (sp.Integer(1), sp.I, sp.Integer(-1))
left = ward_ratio(pts, (1, 1, 1), -sp.I)
check("three legs of weight one at 1, i, -1, stress at -i: ratio 3/4 (the recorded original value)",
      left == sp.Rational(3, 4), left)
right_chiral = ward_ratio(pts, (0, 0, 0), -sp.I)
right_paired = ward_ratio(pts, (1, 0, 1), -sp.I)
check("right-moving part: 0 if the outer legs are one-sided, 1 if they carry right weight one (paired)",
      right_chiral == 0 and sp.simplify(right_paired-1) == 0, (right_chiral, right_paired))

# ---------------------------------------------------------------- P2
print("P2  two placements of the same word")
mod = lambda zs, hs: sp.Abs(sp.simplify(three_point(zs, hs)))
a = (sp.Integer(1), sp.I, sp.Integer(-1))
b = (sp.Integer(1), sp.Integer(-1), -sp.I)
ch = sp.simplify(mod(b, (1, 1, 1))/mod(a, (1, 1, 1)))
pr = sp.simplify(mod(b, (1, 1, 1))*mod(b, (1, 0, 1))/(mod(a, (1, 1, 1))*mod(a, (1, 0, 1))))
check("moduli at (1, -1, -i) against (1, i, -1): equal for one-sided legs, doubled for paired legs",
      ch == 1 and pr == 2, (ch, pr))

# ---------------------------------------------------------------- P3
print("P3  one fermion at a mark")
tau = mp.mpc(0, 1)
N = 40
th1 = lambda z: 2*mp.nsum(lambda n: (-1)**int(n)*mp.exp(1j*PI*tau*(n+mp.mpf(1)/2)**2)*mp.sin((2*n+1)*z), [0, N])
th3 = lambda z: mp.nsum(lambda n: mp.exp(1j*PI*tau*n*n+2j*n*z), [-N, N])
d1 = mp.diff(th1, 0)
S = lambda u: PI*th3(PI*u)*d1/(th3(0)*th1(PI*u))
refl = lambda z: -mp.conj(z)+1                   # recorded rho_+ = h Q kappa
m1, m2 = mp.mpc("0.25", "0.25"), mp.mpc("0.25", "-0.25")
n1 = abs(S(refl(m1)-m1))/(2*PI)
n2 = abs(S(refl(m2)-m2))/(2*PI)
cross = abs(S(refl(m1)-m2))/(2*PI)
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI))
check("norm square of one fermion at a mark of the half annulus: varpi/(2 pi) = 0.417313420837...",
      abs(n1-mp.mpf("0.417313420837")) < mp.mpf(10)**(-11) and abs(n1-varpi/(2*PI)) < mp.mpf(10)**(-25)
      and abs(n2-n1) < mp.mpf(10)**(-25), mp.nstr(n1, 14))
check("the two marks of the half annulus give orthogonal states", cross < mp.mpf(10)**(-25))
print("      a single fermion is an admissible state of the free source; products of fermions change the")
print("      eight charges by integers and never reach the half-integer charges of a spinor current.")
spin = [sp.Rational(s, 2) for s in (1, 1, 1, 1, 1, 1, 1, 1)]
check("spinor weight (1/2, ..., 1/2): norm 2, charges half-integer, not in the integer lattice of fermion products",
      sum(x*x for x in spin) == 2 and all(x.q == 2 for x in spin))

# ---------------------------------------------------------------- P4
print("P4  mirrors of the tetrahedral torus")
om = sp.exp(2*sp.pi*sp.I/3)
basis = (sp.Integer(1), om)


def coords(z):
    """coordinates of z in the basis (1, omega)"""
    im = sp.im(z)/sp.im(om)
    re = sp.re(z)-im*sp.re(om)
    return sp.nsimplify(re), sp.nsimplify(im)


two_circle = 0
count = 0
for k in range(6):
    al = sp.exp(sp.I*sp.pi*k/3)                   # the six linear mirrors z -> al conj(z)
    M = sp.Matrix([[coords(sp.simplify(al*sp.conjugate(bv)))[i] for bv in basis] for i in range(2)])
    if all(x.is_integer for x in M):
        count += 1
        if all((M[i, j]-(1 if i == j else 0)) % 2 == 0 for i in range(2) for j in range(2)):
            two_circle += 1
check("all six linear mirrors preserve the lattice; none acts as the identity modulo two: one fixed circle each",
      count == 6 and two_circle == 0, (count, two_circle))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
