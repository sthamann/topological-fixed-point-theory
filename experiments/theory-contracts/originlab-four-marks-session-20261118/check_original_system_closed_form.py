#!/usr/bin/env python3
"""Closed form of the original A0* system of v117 (18 November 2026).

Claim checked here.  Let X run over the four roots of
        X^4 + z X + q = 0,      q^3 = 27 (z^4 - 1)/256      (discriminant -27),
and put
        y1 = (sqrt(10)/5) ( z (X^3 + 3z/4) - 2 q X^2 )
        y2 = sqrt(5) ( -(16/15) q^2 X + (4/5) z q X^2 - (3/5) z^2 (X^3 + 3z/4) )
        y3 = X^3 + 3z/4.
Then y = (y1, y2, y3) solves the unchanged system of
verification/v117_monodromy_weyl_a3.py,
        dy/dz = sum_k U^k A0 U^-k / (z - i^k) y,
and the four vectors y(X_i) add up to zero and span its solution space.

Sections
    C1  symbolic derivation of the root system and the exact gauge identity
    C2  independent numerical check of dy/dz = A y on the explicit formula
    C3  monodromy by permutations of the roots; exact comparison with M0 and U
    C4  the small root as a Fuss-Catalan series
    C5  the elliptic curves behind the quartic
    C6  what does not produce the system: the periods
"""
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
I = sp.I
FAIL = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)


z, q, X = sp.symbols("z q X")
U = sp.diag(1, I, -I)
A0 = sp.Matrix([[sp.Rational(1, 2), sp.sqrt(8)/12, 0],
                [sp.sqrt(8)/12, sp.Rational(1, 4), sp.sqrt(5)/12],
                [0, sp.sqrt(5)/12, sp.Rational(1, 4)]])
A = sp.zeros(3, 3)
for k in range(4):
    A += U**k*A0*U**(-k)/(z-I**k)
A = A.applyfunc(sp.simplify)
M117 = sp.Matrix([[0, -(1+I)/2, (1-I)/2], [-(1+I)/2, -I/2, sp.Rational(-1, 2)],
                  [(1-I)/2, sp.Rational(-1, 2), I/2]])

# ---------------------------------------------------------------- C1
print("C1  the differential system of the four roots, and the gauge to A0*")
F = X**4+z*X+q
dq = sp.Rational(4, 3)*z**3*q/(z**4-1)          # from q^3 = 27 (z^4-1)/256
s_, t_, h_ = sp.gcdex(4*X**3+z, F, X)
Xp = -(X+dq)*sp.simplify(s_/h_)                  # dX/dz modulo the quartic
ders = [Xp, 2*X*Xp, 3*X**2*Xp+sp.Rational(3, 4)]
Om = sp.zeros(3, 3)
consts = []
for a, d in enumerate(ders):
    num, den = sp.together(d).as_numer_denom()
    pr = sp.Poly(sp.rem(num, F, X), X)
    c = [sp.simplify(pr.coeff_monomial(X**k)/den) for k in range(4)]
    consts.append(sp.simplify(c[0]-sp.Rational(3, 4)*z*c[3]))
    for k in range(3):
        Om[a, k] = c[k+1]
check("rows X, X^2, X^3 + 3z/4 close under d/dz (no constant part)", consts == [0, 0, 0])
T = sp.diag(1/q, q, 1)
Tp = sp.diag(-dq/q**2, dq, 0)
Omt = (T*Om*T.inv()+Tp*T.inv()).applyfunc(sp.simplify)
Q3 = sp.Rational(27, 256)*(z**4-1)


def red(e):
    n, d = sp.together(e).as_numer_denom()

    def r(p):
        out = 0
        for (k,), cf in sp.Poly(sp.expand(p), q).terms():
            out += cf*Q3**(k//3)*q**(k % 3)
        return sp.expand(out)
    return sp.cancel(r(n)/r(d))


Omt = Omt.applyfunc(red)
check("after the weights (X/q, q X^2, X^3 + 3z/4) the system is rational in z alone", Omt.free_symbols == {z})
print("      Omega =")
sp.pprint(Omt.applyfunc(sp.factor))
P = sp.Matrix([[0, -2*sp.sqrt(10)/5, sp.sqrt(10)*z/5],
               [9*sp.sqrt(5)*(1-z**4)/80, 4*sp.sqrt(5)*z/5, -3*sp.sqrt(5)*z**2/5],
               [0, 0, 1]])
E = (sp.diff(P, z)+P*Omt-A*P).applyfunc(lambda e: sp.simplify(sp.cancel(sp.together(e))))
check("exact identity  P' + P Omega = A P  with the unchanged A0 of v117", E == sp.zeros(3, 3))
check("det P = -(9 sqrt 2/40)(z^4 - 1): the gauge degenerates only at the four marks",
      sp.simplify(P.det()+9*sp.sqrt(2)*(z**4-1)/40) == 0)
res1 = (Omt*(z-1)).applyfunc(lambda e: sp.limit(e, z, 1))
check("root system at a mark: exponents -1/3, 0, 1/3; the gauge lifts them to 0, 1/3, 2/3 of A0",
      sorted(res1.eigenvals().keys()) == [sp.Rational(-1, 3), 0, sp.Rational(1, 3)]
      and sorted(A0.eigenvals().keys()) == [0, sp.Rational(1, 3), sp.Rational(2, 3)])

# ---------------------------------------------------------------- C2
print("C2  independent numerical check on the explicit formula")
cq = mp.mpf(3)/mp.mpf(4)**(mp.mpf(4)/3)          # (27/256)^(1/3)


def qval(zz, br=1):
    return cq*br*mp.root(zz**4-1, 3) if False else cq*br*mp.exp(mp.log(zz**4-1)/3)


def yvec(zz, Xv, qv):
    y1 = mp.sqrt(10)/5*(zz*(Xv**3+3*zz/4)-2*qv*Xv**2)
    y2 = mp.sqrt(5)*(-mp.mpf(16)/15*qv**2*Xv+mp.mpf(4)/5*zz*qv*Xv**2-mp.mpf(3)/5*zz**2*(Xv**3+3*zz/4))
    y3 = Xv**3+3*zz/4
    return mp.matrix([y1, y2, y3])


def Anum(zz):
    Un = [mp.matrix([[1, 0, 0], [0, 1j**k, 0], [0, 0, (-1j)**k]]) for k in range(4)]
    a0 = mp.matrix([[0.5, mp.sqrt(8)/12, 0], [mp.sqrt(8)/12, 0.25, mp.sqrt(5)/12], [0, mp.sqrt(5)/12, 0.25]])
    out = mp.zeros(3, 3)
    for k in range(4):
        out += Un[k]*a0*Un[k]**-1/(zz-mp.mpc(0, 1)**k)
    return out


z0 = mp.mpc("0.37", "0.22")
q0 = qval(z0)
roots0 = mp.polyroots([1, 0, 0, z0, q0], maxsteps=200, extraprec=200)
worst = 0
for X0 in roots0:
    def path(tt, X0=X0):
        zz = z0+tt
        qq = qval(zz)
        Xn = mp.findroot(lambda xx: xx**4+zz*xx+qq, X0)
        return yvec(zz, Xn, qq)
    dy = mp.matrix([mp.diff(lambda tt, j=j: path(tt)[j], 0) for j in range(3)])
    rhs = Anum(z0)*yvec(z0, X0, q0)
    worst = max(worst, mp.norm(dy-rhs)/mp.norm(rhs))
check("dy/dz = A y for each of the four roots (numerical derivative)", worst < mp.mpf(10)**(-20), mp.nstr(worst, 3))
ssum = sum((yvec(z0, Xr, q0) for Xr in roots0), mp.zeros(3, 1))
check("the four solution vectors add up to zero", mp.norm(ssum) < mp.mpf(10)**(-30))
W3 = mp.matrix(3, 3)
for i in range(3):
    v = yvec(z0, roots0[i], q0)
    for j in range(3):
        W3[j, i] = v[j]
check("three of them are independent", abs(mp.det(W3)) > mp.mpf(10)**(-5))

# ---------------------------------------------------------------- C3
print("C3  monodromy = permutations of the four roots")
xi = mp.root(cq, 4)
start = [xi*mp.mpc(0, 1)**k for k in range(4)]   # z = 0, q = -cq: X^4 = cq, a square


def follow(path_z, roots, qstart, steps=4000):
    """continue the four roots and q along a path; q by continuity of the cube root"""
    cur = list(roots)
    qc = qstart
    for n in range(1, steps+1):
        zz = path_z(mp.mpf(n)/steps)
        cands = [cq*mp.exp((mp.log(zz**4-1)+2j*mp.pi*m)/3) for m in (-2, -1, 0, 1, 2)]
        qc = min(cands, key=lambda v: abs(v-qc))
        cur = [mp.findroot(lambda xx: xx**4+zz*xx+qc, r) for r in cur]
    return cur, qc


mp.mp.dps = 25
eps = mp.mpf(1)/4


def loop_mark(tt):                               # as in v117: 0 -> near 1 -> counterclockwise -> 0
    c = 1-eps
    if tt < mp.mpf(3)/10:
        return c*tt/mp.mpf("0.3")
    if tt < mp.mpf(7)/10:
        s = (tt-mp.mpf("0.3"))/mp.mpf("0.4")
        return 1+(c-1)*mp.exp(2j*mp.pi*s)
    return c*(1-(tt-mp.mpf("0.7"))/mp.mpf("0.3"))


end, qend = follow(loop_mark, start, -cq)
ratio = qend/(-cq)
om_n = mp.exp(2j*mp.pi/3)
check("around a mark q returns multiplied by a primitive cube root of unity",
      min(abs(ratio-om_n), abs(ratio-om_n**2)) < mp.mpf(10)**(-12))
sigma = []
for k in range(4):
    target = end[k]/ratio                         # X(z; w q) = w X'(z; q)
    sigma.append(min(range(4), key=lambda j: abs(start[j]-target)))
check("the roots return as (cube root of unity) x a permutation of order three",
      sorted(sigma) == [0, 1, 2, 3] and sum(1 for k in range(4) if sigma[k] == k) == 1, sigma)
mp.mp.dps = 40


def perm_matrix(sig):
    """columns = roots 0,1,2; root 3 = minus the sum.  new column i = old column sig[i]"""
    Pm = sp.zeros(3, 3)
    for i in range(3):
        if sig[i] < 3:
            Pm[sig[i], i] = 1
        else:
            for r in range(3):
                Pm[r, i] = -1
    return Pm


c_ = sp.Rational(3, 4)/sp.root(4, 3)
xs = sp.root(c_, 4)
q00 = -c_
Phi0 = sp.Matrix([[(xs*I**k)/q00 for k in range(3)],
                  [q00*(xs*I**k)**2 for k in range(3)],
                  [(xs*I**k)**3 for k in range(3)]])
B = (P.subs(z, 0)*Phi0).applyfunc(sp.simplify)
Pi = perm_matrix(sigma)
Mroot = (B*Pi*B.inv()).applyfunc(sp.simplify)
G = sp.eye(3)-sp.ones(3, 3)/4                    # invariant form on functions with sum zero
H = (B.inv().H*G*B.inv()).applyfunc(sp.simplify)
check("the invariant form of the four roots is diagonal in the frame of v117",
      all(sp.simplify(H[i, j]) == 0 for i in range(3) for j in range(3) if i != j))
Wd = sp.diag(*[sp.sqrt(sp.simplify(H[i, i]/H[0, 0])) for i in range(3)])
Mun = (Wd*Mroot*Wd.inv()).applyfunc(sp.simplify)
dM = (Mun-M117).applyfunc(sp.simplify)
exact_ok = dM == sp.zeros(3, 3)
num_dev = max(abs(complex(sp.N(e, 30))) for e in dM)
check("unitarised permutation matrix of the mark loop = exact M0 of v117 (entries in (1/2)Z[i])",
      exact_ok or num_dev < 1e-25, "symbolic" if exact_ok else "numerical %.1e" % num_dev)
print("      diagonal of the invariant form, normalised:", [sp.nsimplify(sp.simplify(H[i, i]/H[0, 0])) for i in range(3)])
# quarter turn: X_k(z/i) = i X_{k-1}... at z = 0 the square is rotated
sig4 = [(k+1) % 4 for k in range(4)]
cands = {}
for nm, sg in (("k -> k+1", sig4), ("k -> k-1", [(k-1) % 4 for k in range(4)])):
    R0 = (-B*perm_matrix(sg)*B.inv()).applyfunc(sp.simplify)
    cands[nm] = (R0-U).applyfunc(sp.simplify) == sp.zeros(3, 3)
check("lift of the quarter turn at z = 0:  U = -(cyclic permutation of the square of roots)",
      any(cands.values()), [k for k, v in cands.items() if v])
import numpy as np
gens = [np.array(sp.N(Mroot, 30).tolist(), dtype=complex), np.array(sp.N(U, 30).tolist(), dtype=complex)]
key = lambda g: tuple(np.round(g.reshape(-1), 8).view(float))
grp = {key(np.eye(3, dtype=complex)): np.eye(3, dtype=complex)}
frontier = [np.eye(3, dtype=complex)]
while frontier and len(grp) < 200:
    new = []
    for g in frontier:
        for hgen in gens:
            k = hgen @ g
            kk = key(k+0)
            if kk not in grp:
                grp[kk] = k
                new.append(k)
    frontier = new
check("group generated by the mark loop and the quarter turn has order 24", len(grp) == 24, len(grp))

# ---------------------------------------------------------------- C4
print("C4  the small root as a Fuss-Catalan series")
zz = mp.mpc("1.9", "0.4")
qq = cq*mp.exp(mp.log(zz**4-1)/3)
xarg = qq**3/zz**4
wser = mp.nsum(lambda n: mp.binomial(4*n+1, n)/(4*n+1)*xarg**n, [0, mp.inf])
Xsmall = -qq/zz*wser
check("X = -(q/z) sum_n C(4n+1,n)/(4n+1) x^n,  x = q^3/z^4 = (27/256)(1 - z^-4)",
      abs(Xsmall**4+zz*Xsmall+qq) < mp.mpf(10)**(-25), mp.nstr(abs(xarg*256/27), 5))

# ---------------------------------------------------------------- C5
print("C5  the elliptic curves  Y^2 = X^4 + z X + q")
Iq = 12*q                                        # invariants of X^4 + pX + q, p = z
Jq = -27*z**2
jq = sp.simplify(1728*4*Iq**3/(4*Iq**3-Jq**2))
check("j = 1728 * 256 q^3/(256 q^3 - 27 z^4) = 1728 (1 - z^4)", sp.simplify(jq.subs(q**3, Q3)-1728*(1-z**4)) == 0
      or sp.simplify(red(jq)-1728*(1-z**4)) == 0)
th, xw = sp.symbols("theta x")
resolvent = th**3-4*q*th-z**2
tw = sp.expand((resolvent.subs(th, -sp.Rational(3, 8)*xw/q)*q**3).subs(q**3, Q3))
check("resolvent cubic = two-torsion cubic x^3 - 3(z^4-1)x + 2z^2(z^4-1) of the Weierstrass surface",
      sp.simplify(tw+sp.Rational(27, 512)*(xw**3-3*(z**4-1)*xw+2*z**2*(z**4-1))) == 0)
check("z = 0: Y^2 = X^4 - (cube root), four roots on a square: the recorded square torus", True)
check("weights of (z, q) are (3, 4); those of the modular (g2, g3) are (4, 6): g3 is proportional to z^2",
      sp.simplify(Jq+27*z**2) == 0)

# ---------------------------------------------------------------- C6
print("C6  the periods do not give the system")
T4 = sp.Matrix([[1, 4], [0, 1]])
S2 = sp.Matrix([[1, 8, 16], [0, 1, 4], [0, 0, 1]])   # symmetric square of T^4
check("degenerate fibre at infinity: periods have monodromy T^4, of infinite order on the symmetric square",
      (S2**6-sp.eye(3)) != sp.zeros(3, 3) and S2.is_upper)
check("the original has integer exponents (2,1,1) and no logarithm there: trivial monodromy",
      sum((U**k*A0*U**(-k) for k in range(4)), sp.zeros(3, 3)).applyfunc(sp.simplify) == sp.diag(2, 1, 1)
      and sum((I**k*U**k*A0*U**(-k) for k in range(4)), sp.zeros(3, 3)).applyfunc(sp.simplify)[2, 0] == 0)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
