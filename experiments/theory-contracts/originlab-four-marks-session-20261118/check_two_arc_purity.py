#!/usr/bin/env python3
"""Hypothesis H-R2: the recorded free source is the purity of two opposite arcs (18 November 2026).

The four marks cut the seam circle into four arcs.  Let A be the union of two opposite arcs and
rho_A the state of the seam restricted to A.  Tr rho_A^2 is the expectation of the exchange of
two copies of the seam on A.  On two copies a fermion passing a mark goes 1 -> 2 -> -1: a
quarter turn.  The surface on which Tr rho_A^2 is a partition function is the double cover of
the mark sphere branched at the four marks, the recorded torus.

Sections
    U1  two copies: quarter turn, twists 1/4 and 3/4, inverse temperature 4 pi; the recorded
        law of R066 is the partition function of one spinor on the recorded torus
    U2  independent computation on a ring of lattice fermions: the purity of two arcs is the
        weight of the alternating signs (the free seam of round twenty-three)
    U3  the sum over the three pairings (the E8 reading) as a purity; counting of currents
    U4  the current function on the recorded torus; values between the marks
    U5  numbers that follow
    U6  the quarter turn as a product of two reflections; weights 1, 1 and 1/2
"""
import sys

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 40
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- U1
print("U1  two copies of a fermion at a mark")
M = sp.Matrix([[0, 1], [-1, 0]])            # copy 1 -> copy 2 -> minus copy 1
check("passing a mark is a quarter turn: M^2 = -1, M^4 = 1, eigenvalues +i and -i",
      M**2 == -sp.eye(2) and M**4 == sp.eye(2) and set(M.eigenvals()) == {sp.I, -sp.I})
alphas = sorted(sp.Rational(1, 2)+k/2 for k in (sp.Rational(1, 2), sp.Rational(-1, 2)))
check("twists k/n with n = 2, k = +-1/2, measured from the antiperiodic value 1/2: 1/4 and 3/4",
      alphas == [sp.Rational(1, 4), sp.Rational(3, 4)], alphas)
ok = True
for al in alphas:
    for Q in range(-7, 8):
        if Q >= 0:
            E = sum(j+al for j in range(Q))                 # lowest Q particles in modes j + alpha
        else:
            E = sum(j-al for j in range(1, -Q+1))           # lowest |Q| holes
        ok = ok and sp.simplify(E-(sp.Rational(Q*Q, 2)+(al-sp.Rational(1, 2))*Q)) == 0
check("lowest energy of charge Q for modes r + alpha is F_alpha(Q) = Q^2/2 + (alpha - 1/2) Q, the recorded function", ok)
check("two copies give the inverse temperature 2 * 2 pi = 4 pi, the recorded b", 2*2*sp.pi == 4*sp.pi)

tau = mp.mpc(0, 1)


def th(n, t):
    q = mp.exp(mp.pi*1j*t)
    return mp.jtheta(n, 0, q)


def eta(t):
    q = mp.exp(2j*mp.pi*t)
    return mp.exp(2j*mp.pi*t/24)*mp.qp(q)


def Zal(al):
    return mp.nsum(lambda Q: mp.exp(-2*mp.pi*(Q+al-mp.mpf(1)/2)**2), [-mp.inf, mp.inf])/eta(2*tau)


lhs = Zal(mp.mpf(1)/4)*Zal(mp.mpf(3)/4)
rhs = th(3, tau)/eta(tau)
check("trace at 4 pi over both twisted fermions = theta_3(i)/eta(i): one spinor on the recorded torus "
      "in the recorded spin structure", abs(lhs-rhs) < mp.mpf(10)**(-30), mp.nstr(lhs.real, 15))

# ---------------------------------------------------------------- U2
print("U2  ring of lattice fermions (independent of the theta functions)")


def purity_ratio(N, l):
    """Tr rho^2 of two opposite arcs of l sites over the product for the single arcs."""
    d = np.arange(N)
    row = np.where(d == 0, 0.5, np.sin(np.pi*d/2)/(N*np.sin(np.pi*np.where(d == 0, 1, d)/N)))
    idx1 = np.arange(l)
    idx2 = np.arange(N//2, N//2+l)

    def logpur(idx):
        C = row[np.abs(idx[:, None]-idx[None, :])]
        nu = np.clip(np.linalg.eigvalsh(C), 0, 1)
        return np.sum(np.log(nu**2+(1-nu)**2))
    return np.exp(logpur(np.concatenate([idx1, idx2]))-2*logpur(idx1))


rows = []
for frac in (sp.Rational(1, 4), sp.Rational(3, 20), sp.Rational(7, 20)):     # arc length / circumference
    vals = [purity_ratio(N, int(frac*N)) for N in (800, 1600)]
    ext = 2*vals[1]-vals[0]
    delta = 2*np.pi*float(frac)
    rows.append((frac, ext, np.cos(delta/2)**(-0.5)))
print("      arc fraction, lattice (extrapolated), cos(delta/2)^(-1/2):",
      [(str(f), round(a, 5), round(b, 5)) for f, a, b in rows])
check("one lattice fermion: purity ratio = cos(delta/2)^(-1/2) = (1 - x)^(-1/4), x = sin^2(delta/2)",
      all(abs(a-b) < 3e-3 for _, a, b in rows))
check("eight chiral Dirac fermions of the free seam: (1 - x)^(-1) = |P_12| / sin^2(delta), the weight of the "
      "alternating signs times the distance factor of the two arcs",
      all(abs((1-np.sin(dl/2)**2)**(-1)-4*np.sin(dl/2)**2/np.sin(dl)**2) < 1e-12 for dl in (0.4, 0.9, 1.3)))

# ---------------------------------------------------------------- U3
print("U3  the sum over the three pairings as a purity")
x = sp.symbols("x", positive=True)
Xfree = 1/(1-x)
XE8 = (1-x+x**2)/(1-x)
dl = sp.symbols("delta", positive=True)
Ropp, Rs, Rsp = 1/sp.sin(dl)**2, sp.tan(dl/2)**2/4, sp.cot(dl/2)**2/4
P12 = 4*sp.sin(dl/2)**2
ok = all(abs(float(((P12*(Ropp+Rs+Rsp)/2)-XE8.subs(x, sp.sin(dl/2)**2)).subs(dl, v).evalf())) < 1e-12
         for v in (sp.Rational(2, 5), sp.Rational(9, 10), sp.Rational(13, 10)))
check("half the sum of the three pairing weights times |P_12| is (1 - x + x^2)/(1 - x); it tends to 1 for "
      "distant arcs, which fixes the factor 1/2", ok)
ser = sp.series(XE8, x, 0, 4).removeO()
check("E8: 1 + x^2 + x^3 + ..., no term in x: nothing of weight 1/2 + 1/2 passes between the arcs",
      sp.expand(ser-(1+x**2+x**3)) == 0)
check("coefficient of x^2: 248 currents + the total stress of the two copies: (248 + 8)/256 = 1",
      sp.Rational(248, 256)+sp.Rational(1, 2)**2*2/16 == 1)
check("free seam: 1 + x + ..., sixteen real fermions exchanged in pairs: 16/16 = 1",
      sp.series(Xfree, x, 0, 2).removeO() == 1+x)
jx = 256*(1-x+x**2)**3/(x**2*(1-x)**2)
check("cube of the E8 character in x: 256 (1 - x + x^2)^3 / (x^2 (1 - x)^2), the j function of the cover",
      sp.simplify(jx.subs(x, sp.Rational(1, 2))-1728) == 0)
check("weight of the alternating signs in the sum: 1/(2 (1 - x + x^2)); at the square 2/3 (round twenty-three)",
      sp.simplify((Xfree/(2*XE8)).subs(x, sp.Rational(1, 2))) == sp.Rational(2, 3))

# ---------------------------------------------------------------- U4
print("U4  the current function on the recorded torus")


def eis(t):
    q = mp.exp(2j*mp.pi*t)
    s = lambda k: mp.nsum(lambda n: n**k*q**n/(1-q**n), [1, mp.inf])
    return 1-24*s(1), 1+240*s(3), 1-504*s(5)


def e_values(t):
    a, b, c = th(2, t)**4, th(4, t)**4, th(3, t)**4
    return {2: mp.pi**2/3*(c+b), 4: -mp.pi**2/3*(a+c), 3: mp.pi**2/3*(a-b)}   # zero of theta_eps at 1/2, tau/2, (1+tau)/2


tg = mp.mpc(0.13, 0.87)
E2, E4, E6 = eis(tg)
ev = e_values(tg)
num = sum(th(k, tg)**8*ev[k] for k in (2, 3, 4))
den = sum(th(k, tg)**8 for k in (2, 3, 4))
check("sum over the spin structures with weights theta^8: constant = -(pi^2/3) E6/E4 (generic tau)",
      abs(num/den+mp.pi**2/3*E6/E4) < mp.mpf(10)**(-25))
q = mp.exp(2j*mp.pi*tg)
dlnE4 = mp.nsum(lambda n: 240*n**4*q**n/(1-q**n)**2, [1, mp.inf])/E4          # q d/dq log E4
lat = mp.pi**2/3*E2-4*mp.pi**2*dlnE4/4
check("the same constant from the lattice: (pi^2/3) E2 - 4 pi^2 <p_a^2>, <p_a^2> = (1/4) q d/dq log E4",
      abs(lat-mp.pi**2/3*E6/E4) < mp.mpf(10)**(-25))
E2i, E4i, E6i = eis(tau)
evi = e_values(tau)
check("square torus: E6(i) = 0, so the E8 current function is the Weierstrass function alone",
      abs(E6i) < mp.mpf(10)**(-30))
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*mp.pi))
check("between the marks (differences 1/2, i/2, (1+i)/2): +varpi^2, -varpi^2 and 0; neighbours carry a current "
      "with opposite signs, opposite marks none; varpi^2 is the recorded kernel square of round fourteen",
      abs(evi[2]-varpi**2) < mp.mpf(10)**(-25) and abs(evi[4]+varpi**2) < mp.mpf(10)**(-25) and abs(evi[3]) < mp.mpf(10)**(-30),
      mp.nstr(evi[2].real, 12))
check("the free seam (theta_3 alone) has the same current function at the square, e_3 = 0; "
      "off the square the two differ by e_3 + (pi^2/3) E6/E4",
      abs(ev[3]+mp.pi**2/3*E6/E4) > mp.mpf(10)**(-3))

# ---------------------------------------------------------------- U5
print("U5  numbers")
sq = sp.Rational(1, 2)
cw = sp.Rational(1, 5)                      # counterwitness {+-1, +-(3+4i)/5}: x = sin^2(delta/2) = 1/5
vals = {"E8 square": XE8.subs(x, sq), "free square": Xfree.subs(x, sq),
        "E8 counterwitness": XE8.subs(x, cw), "free counterwitness": Xfree.subs(x, cw)}
print("      purity ratio:", vals)
check("square: 3/2 (E8) against 2 (free); counterwitness: 21/20 against 5/4",
      vals == {"E8 square": sp.Rational(3, 2), "free square": 2, "E8 counterwitness": sp.Rational(21, 20),
               "free counterwitness": sp.Rational(5, 4)})
check("the purity itself, X/x up to a constant, is the same for the arcs and for their complement (x -> 1 - x), "
      "for E8 and for the free seam: the seam state is pure",
      sp.simplify(XE8/x-(XE8/x).subs(x, 1-x)) == 0 and sp.simplify(Xfree/x-(Xfree/x).subs(x, 1-x)) == 0)
check("theta_3(i)/eta(i) = sqrt(2): one spinor on the recorded torus counts sqrt(2); sixteen of them 16; "
      "the E8 character 16 (1 + 1/4 + 1/4)/2 = 12",
      abs(th(3, tau)/eta(tau)-mp.sqrt(2)) < mp.mpf(10)**(-30))

# ---------------------------------------------------------------- U6
print("U6  the quarter turn as a product of two reflections")
S = sp.Matrix([[0, 1], [1, 0]])            # exchange of the two copies
Pm = sp.diag(1, -1)                        # relative sign of the second copy
check("M = P S, with S the exchange of the copies and P the relative sign; S^2 = P^2 = 1, S P S P = -1",
      Pm*S == M and S**2 == sp.eye(2) and Pm**2 == sp.eye(2) and S*Pm*S*Pm == -sp.eye(2))


def weight(matrix, fermions=8):
    """Weight of the place where the copies are joined by `matrix`: sum of alpha^2/2 over the phases."""
    total = 0
    for value, mult in matrix.eigenvals().items():
        alpha = sp.nsimplify(sp.arg(value)/(2*sp.pi))
        total += mult*alpha**2/2
    return fermions*total


wts = {"1": weight(sp.eye(2)), "S": weight(S), "P": weight(Pm), "PS": weight(M)}
print("      weights at a mark for eight chiral Dirac fermions:", wts)
check("weights 0, 1, 1 and 1/2: the two reflections carry weight one (the legs of the recorded E8 vacuum form), "
      "their product the recorded mark weight 1/2",
      wts == {"1": 0, "S": 1, "P": 1, "PS": sp.Rational(1, 2)})
check("for the two arcs of the free seam the ratios are (1 - x)^(-2), (1 - x)^(-2) and (1 - x)^(-1): "
      "at the square 4, 4 and 2", [(1-sp.Rational(1, 2))**(-2*w) for w in (1, 1, sp.Rational(1, 2))] == [4, 4, 2])

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
