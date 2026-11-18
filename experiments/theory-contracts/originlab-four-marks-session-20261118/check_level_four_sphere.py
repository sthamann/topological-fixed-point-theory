#!/usr/bin/env python3
"""Session check of 18 November 2026: hypothesis H-X4.

The sphere that carries the four marks 1, i, -1, -i is read as the modular
curve X(4) = H / Gamma(4), with coordinate

    f(tau) = b/a,   a = theta_3(2 tau),  b = theta_2(2 tau)

(the two weight-1/2 theta functions of one lattice direction; a/eta, b/eta are
the two characters of SU(2) level one).  Checked here, numerically at 40
digits unless stated otherwise:

X1  modular action: T multiplies f by i, S acts as f -> (1-f)/(1+f);
    the six cusps are 0, infinity, 1, i, -1, -i
X2  hyperbolic area of X(4) is 8 pi; the weight-1/2 bundle has degree one
X3  Klein's octahedral invariants in (a, b) are 2 eta^6, E4, E6;
    the degree-8 one is the weight enumerator of the Hamming code
X4  projective image of the modular group on (a, b) has order 24; on the
    symmetric square it gives the recorded 1/9, 4/9, 4/9 table
X5  special points: tau = i, tau = omega
X6  the recorded base curve w^3 = (z^2-1)/(z^2+1) under z = f
X7  three-fold covers of the sphere branched in the octahedron vertices
"""
import itertools
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
PI = mp.pi


def thetas(tau):
    """theta_2, theta_3, theta_4 at argument tau, as explicit sums (branch of q^(1/4) fixed by tau)."""
    N = 40
    t2 = mp.nsum(lambda n: mp.exp(1j*PI*tau*(n+mp.mpf(1)/2)**2), [-N, N])
    t3 = mp.nsum(lambda n: mp.exp(1j*PI*tau*n**2), [-N, N])
    t4 = mp.nsum(lambda n: (-1)**int(n)*mp.exp(1j*PI*tau*n**2), [-N, N])
    return t2, t3, t4


def ab(tau):
    t2, t3, _ = thetas(2*tau)
    return t3, t2


def eta(tau):
    q = mp.exp(2j*PI*tau)
    return mp.exp(2j*PI*tau/24)*mp.qp(q)


def eis(k, tau):
    q = mp.exp(2j*PI*tau)
    c = {4: 240, 6: -504}[k]
    return 1+c*mp.nsum(lambda n: n**(k-1)*q**n/(1-q**n), [1, mp.inf])


f = lambda tau: ab(tau)[1]/ab(tau)[0]
close = lambda x, y, tol=mp.mpf(10)**-30: abs(x-y) < tol*(1+abs(y))
samples = [mp.mpc(0.13, 0.71), mp.mpc(-0.37, 1.2), mp.mpc(0.45, 0.33)]

print("X1. modular action on f = theta_2(2 tau) / theta_3(2 tau)")
for tau in samples:
    assert close(f(tau+1), 1j*f(tau))
    assert close(f(-1/tau), (1-f(tau))/(1+f(tau)))
print("    f(tau + 1) = i f(tau);  f(-1/tau) = (1 - f)/(1 + f): verified at three points")
print("    cusp tau = i*infinity: f ->", mp.nstr(abs(f(mp.mpc(0, 12))), 5), "| cusp tau = 0: f ->", mp.nstr(f(mp.mpc(0, mp.mpf(1)/60)), 12))
orbit = {sp.Integer(0)}
S_ = lambda z: sp.oo if z == -1 else (sp.Integer(-1) if z == sp.oo else sp.simplify((1-z)/(1+z)))
T_ = lambda z: sp.oo if z == sp.oo else sp.simplify(sp.I*z)
for _ in range(6):
    orbit |= {S_(z) for z in orbit} | {T_(z) for z in orbit}
print("    orbit of the cusp f = 0 under S and T:", sorted(map(str, orbit)))
assert len(orbit) == 6

print("X2. area and flux")
index = 24                                              # |PSL2(Z/4)|
print("    index of Gamma(4) in PSL2(Z):", index, "| hyperbolic area = index * pi/3 = 8 pi = 2 pi (6 cusps - 2)")
assert sp.Rational(index, 3) == 8 and 2*(6-2) == 8
print("    curvature form of a weight-k form: k/(4 pi) dx dy / y^2; k = 1/2 gives 1/(8 pi)")
print("    flux of the weight-1/2 bundle over X(4): (1/(8 pi)) * 8 pi = 1; degree = k * index / 12 =", sp.Rational(1, 2)*index/12)
print("    over H/Gamma(2) (area 2 pi): 1/4;  over H/PSL2(Z): 1/24")
print("    sections: weight k on Gamma(4) = polynomials of degree 2k in (a, b): dimensions",
      {str(k): 2*k+1 for k in (sp.Rational(1, 2), 1, 4)})

print("X3. octahedral invariants")
x, y = sp.symbols("x y")
V6 = x*y*(x**4-y**4)
W8 = x**8+14*x**4*y**4+y**8
E12 = x**12-33*x**8*y**4-33*x**4*y**8+y**12
assert sp.expand(W8**3-E12**2-108*V6**4) == 0
print("    syzygy  W8^3 - E12^2 = 108 V6^4  (exact)")
for tau in samples:
    a, b = ab(tau)
    assert close(a*b*(a**4-b**4), 2*eta(tau)**6)
    assert close(a**8+14*a**4*b**4+b**8, eis(4, tau))
    assert close(a**12-33*a**8*b**4-33*a**4*b**8+b**12, eis(6, tau))
print("    a b (a^4 - b^4) = 2 eta^6;  a^8 + 14 a^4 b^4 + b^8 = E4;  a^12 - 33(...) + b^12 = E6: verified")
print("    hence 1728 eta^24 = E4^3 - E6^2 is Klein's syzygy (108 * 16 = 1728)")
# Hamming code enumerator
import numpy as np
Gm = np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1], [1, 0, 1, 0, 1, 0, 1, 0]])
words = {tuple((np.array(c) @ Gm) % 2) for c in itertools.product((0, 1), repeat=4)}
wts = {}
for wd in words:
    wts[sum(wd)] = wts.get(sum(wd), 0)+1
assert wts == {0: 1, 4: 14, 8: 1}
print("    weight enumerator of the [8,4,4] code: x^8 + 14 x^4 y^4 + y^8 = W8: the E8 theta function is")
print("    the degree-8 invariant; its zeros on the sphere are the 8 cube vertices (f^8 + 14 f^4 + 1 = 0)")

print("X4. the modular group on the qubit (a, b) and on its symmetric square")
Tm = sp.diag(1, sp.I)
Sm = sp.Matrix([[1, 1], [1, -1]])/sp.sqrt(2)
def proj_key(Mx):
    Mx = sp.simplify(Mx)
    piv = next(v for v in Mx if v != 0)
    return tuple(sp.nsimplify(sp.simplify(v/piv)) for v in Mx)
grp = {proj_key(sp.eye(2)): sp.eye(2)}
frontier = [sp.eye(2)]
while frontier:
    nxt = []
    for g_ in frontier:
        for h_ in (Tm, Sm):
            n_ = sp.simplify(h_*g_)
            k_ = proj_key(n_)
            if k_ not in grp:
                grp[k_] = n_
                nxt.append(n_)
    frontier = nxt
print("    projective group generated by S (Hadamard) and T (diag(1, i)): order", len(grp))
assert len(grp) == 24
def sym2(Mx):
    p, q, r, s = Mx[0, 0], Mx[0, 1], Mx[1, 0], Mx[1, 1]
    # orthonormal basis x^2, sqrt2 xy, y^2
    return sp.Matrix([[p*p, sp.sqrt(2)*p*q, q*q], [sp.sqrt(2)*p*r, p*s+q*r, sp.sqrt(2)*q*s], [r*r, sp.sqrt(2)*r*s, s*s]])
ST = sym2(Sm*Tm)
T2 = sym2(Tm)
ev = ST.eigenvects()
vecs = [(sp.simplify(e[2][0]/sp.sqrt((e[2][0].H*e[2][0])[0, 0]))) for e in ev]
assert len(vecs) == 3
table = sp.Matrix(3, 3, lambda i_, j_: sp.nsimplify(sp.simplify(sp.Abs((vecs[i_].H*T2*vecs[j_])[0, 0])**2), rational=True))
print("    order of S T on the symmetric square:", next(n_ for n_ in range(1, 13) if sp.simplify(ST**n_ - (ST**n_)[0, 0]*sp.eye(3)) == sp.zeros(3)),
      "| transition table of T between its eigenlines:", sorted(map(str, table)))
assert sorted(map(str, table)) == sorted(["1/9"]*3+["4/9"]*6)

print("X5. special points")
fi = f(mp.mpc(0, 1))
print("    tau = i: f =", mp.nstr(fi, 20), "= sqrt(2) - 1 = tan(pi/8):", close(fi, mp.sqrt(2)-1),
      "| E6 = 0 there:", abs(eis(6, mp.mpc(0, 1))) < mp.mpf(10)**-30)
om = mp.exp(2j*PI/3)
fo = f(om)
print("    tau = omega: f^8 + 14 f^4 + 1 =", mp.nstr(abs(fo**8+14*fo**4+1), 5), "(a cube vertex); |f| =", mp.nstr(abs(fo), 12))
chi = lambda tau: eis(4, tau)/eta(tau)**8
print("    E8 character at tau = i:", mp.nstr(chi(mp.mpc(0, 1)), 15),
      "| as a function on the sphere: chi^3 = 16 (1 + 14 f^4 + f^8)^3 / (f^4 (1 - f^4)^4):",
      all(close(chi(t_)**3, 16*(1+14*f(t_)**4+f(t_)**8)**3/(f(t_)**4*(1-f(t_)**4)**4)) for t_ in samples))

print("X6. the recorded base curve under z = f")
for tau in samples:
    t2, t3, t4 = thetas(tau)
    assert close((f(tau)**2-1)/(f(tau)**2+1), -t4**2/t3**2)
print("    (f^2 - 1)/(f^2 + 1) = - theta_4(tau)^2 / theta_3(tau)^2: verified; so w^3 = - theta_4^2/theta_3^2")

print("X7. three-fold covers branched in octahedron vertices (local exponents mod 3 at 1, i, -1, -i, 0, infinity)")
covers = {"recorded base curve  w^3 = (z^2-1)/(z^2+1)": (1, 2, 1, 2, 0, 0),
          "special Picard curve y^3 = z^4 - 1": (1, 1, 1, 1, 0, 2),
          "character curve      t^3 = z - z^5": (1, 1, 1, 1, 1, 1)}
for nm, ex in covers.items():
    assert sum(ex) % 3 == 0
    nb = sum(1 for e in ex if e)
    genus = (2-6+2*nb)//2
    print("    %-44s exponents %s | branch points %d | genus %d | invariant under all 24 rotations: %s"
          % (nm, ex, nb, genus, len(set(ex)) == 1))
print("    Gamma(4)/Gamma(12) = PSL2(Z/3) = A4 has exactly two nontrivial characters of order 3; both give the")
print("    character curve (adjoining j^(1/3) = chi_E8). The other two covers have level 12 but do not contain")
print("    Gamma(12): they are not congruence covers of X(4).")

print("X8. signature of the hermitian form over Z[omega] on the first cohomology (Chevalley-Weil)")
from fractions import Fraction as Fr
def signature(ex):
    frac = lambda p, q: Fr(p % q, q)
    d1 = -1+sum(frac(e, 3) for e in ex if e)            # forms on which the deck map acts by omega^2
    d2 = -1+sum(frac(2*e, 3) for e in ex if e)          # forms on which it acts by omega
    return int(d2), int(d1)
for nm, ex in covers.items():
    print("    %-44s holomorphic forms by deck eigenvalue: %s" % (nm, signature(ex)))
assert signature((1, 1, 1, 1, 1, 1)) == (3, 1) and signature((1, 2, 1, 2, 0, 0)) == (1, 1) and signature((1, 1, 1, 1, 0, 2)) == (2, 1)
print("    n branch points with equal exponents (n divisible by 3): signature (2n/3 - 1, n/3 - 1):",
      {n_: signature((1,)*n_) for n_ in (3, 6, 9, 12)})
print("    flux of the weight-1/2 bundle on a sphere with n cusps: (n - 2)/4:", {n_: str(Fr(n_-2, 4)) for n_ in (3, 4, 6, 9)})
print("    n = 6 is the only case with unit flux, and the only equal-exponent case with exactly one reversed direction")
# explicit differentials on t^3 = z - z^5: z^k dz/t^2 (k = 0,1,2) and dz/t
print("    explicit basis on t^3 = z - z^5: dz/t^2, z dz/t^2, z^2 dz/t^2 (three) and dz/t (one)")
# Fermat-quotient reading: t^12 = x (1-x)^4 with x = z^4; characters k with <k/12>+<4k/12>+<7k/12> = 1
hol = [k for k in range(1, 12) if all((k*a_) % 12 for a_ in (1, 4, 7)) and
       Fr(k % 12, 12)+Fr((4*k) % 12, 12)+Fr((7*k) % 12, 12) == 1]
print("    as a cyclic cover of degree 12 of the t-sphere (t = z^4): holomorphic characters", hol, "of 12")
assert hol == [1, 2, 4, 7]
print("    types: {1, 7} is stable under multiplication by 7, which fixes omega: every factor of the Jacobian")
print("    has complex multiplication by Z[omega] (four hexagonal tori up to isogeny)")
print("ALL CHECKS PASSED")
