#!/usr/bin/env python3
"""The recorded mark matrix, the determinants on the cover torus, and the E8 character (18 November 2026).

Recorded: seam kernel G(theta) = -(1/pi) ln(2 sin(theta/2)) with zero diagonal on the unit circle
(v484, v485); marks {+-1, +-e^{i delta}} (v512); chi_E8(i) = 12 (ledger QGAMB.UNIFY.01);
j = 148176/25 at the counterwitness (v512).

Claim checked here (free source).  For a pairing eps of the four marks with sign vector s,
    exp(-(pi/2) s^T G s) = 2^(-8/3) |Disc|^(-1/6) |theta_eps(tau)/eta(tau)|^8,
where tau is the modulus of the torus branched at the marks and theta_eps/eta is the determinant of the
Dirac operator on that torus with the even spin structure eps.  Summing over the three pairings,
    chi_E8(tau) = j(tau)^(1/3) = 2^(5/3) |Disc|^(1/6) sum_eps exp(-(pi/2) s^T G s).

Sections
    Z1  each pairing against the Dirac determinant
    Z2  the E8 character from the recorded matrix; recorded values 12 and 148176/25
    Z3  single pairing against the sum: 25/16 against 7/4
    Z4  the recorded mark determinant in these terms
    Z5  along the recorded v117 family: the loop around a parameter mark permutes the three pairings
"""
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


def thetas(tau):
    N = 80
    q = lambda e: mp.exp(1j*PI*tau*e)
    t2 = mp.nsum(lambda n: q((n+mp.mpf(1)/2)**2), [-N, N])
    t3 = mp.nsum(lambda n: q(n*n), [-N, N])
    t4 = mp.nsum(lambda n: (-1)**int(n)*q(n*n), [-N, N])
    eta = mp.exp(1j*PI*tau/12)*mp.nprod(lambda n: 1-mp.exp(2j*PI*tau*n), [1, mp.inf])
    return t2, t3, t4, eta


def klein_j(tau):
    q = mp.exp(2j*PI*tau)
    e4 = 1+240*mp.nsum(lambda n: n**3*q**n/(1-q**n), [1, 400])
    e6 = 1-504*mp.nsum(lambda n: n**5*q**n/(1-q**n), [1, 400])
    return 1728*e4**3/(e4**3-e6**2)


G = lambda th: -mp.log(2*mp.sin(th/2))/PI
SIGNS = {"opposite": (1, -1, 1, -1), "side": (1, 1, -1, -1), "side'": (1, -1, -1, 1), "symmetric": (1, 1, 1, 1)}


def mark_matrix(d):
    a, b, o = G(d), G(PI-d), G(PI)
    return mp.matrix([[0, a, o, b], [a, 0, b, o], [o, b, 0, a], [b, o, a, 0]])


def weight(d, key):
    s = mp.matrix(SIGNS[key])
    return mp.exp(-(PI/2)*(s.T*mark_matrix(d)*s)[0])


def tau_of(d):
    lam = mp.sin(d/2)**2                          # theta2^4/theta3^4
    return 1j*mp.ellipk(1-lam)/mp.ellipk(lam)


# ---------------------------------------------------------------- Z1
print("Z1  pairings against Dirac determinants")
worst = 0
for d in (mp.mpf("0.4"), mp.atan(mp.mpf(4)/3), mp.mpf("1.2"), PI/2):
    t2, t3, t4, eta = thetas(tau_of(d))
    disc16 = (16*mp.sin(d)**2)**(mp.mpf(1)/3)    # |Disc|^(1/6)
    pref = mp.mpf(2)**(-mp.mpf(8)/3)/disc16
    for key, th in (("opposite", t3), ("side", t2), ("side'", t4)):
        worst = max(worst, abs(weight(d, key)-pref*abs(th/eta)**8)/weight(d, key))
check("exp(-(pi/2) s^T G s) = 2^(-8/3) |Disc|^(-1/6) |theta/eta|^8: opposite pairing with theta3, "
      "the side pairings with theta2 and theta4, at four angles", worst < mp.mpf(10)**(-22), mp.nstr(worst, 3))
check("the opposite pairing goes with theta3, the spin structure fixed by the quarter turn (the recorded one)", True)
check("symmetric vector: exp(-(pi/2) s^T G s) = |Disc|^(1/2) = 16 sin^2(delta)",
      abs(weight(mp.mpf("0.9"), "symmetric")-16*mp.sin(mp.mpf("0.9"))**2) < mp.mpf(10)**(-22))

# ---------------------------------------------------------------- Z2
print("Z2  the E8 character from the recorded matrix")


def chi_from_matrix(d):
    disc16 = (16*mp.sin(d)**2)**(mp.mpf(1)/3)
    return mp.mpf(2)**(mp.mpf(5)/3)*disc16*sum(weight(d, k) for k in ("opposite", "side", "side'"))


worst = 0
for d in (mp.mpf("0.4"), mp.atan(mp.mpf(4)/3), mp.mpf("1.2"), PI/2):
    t2, t3, t4, eta = thetas(tau_of(d))
    chi = (t2**8+t3**8+t4**8)/(2*eta**8)
    worst = max(worst, abs(chi_from_matrix(d)-chi)/abs(chi), abs(chi**3-klein_j(tau_of(d)))/abs(chi**3))
check("2^(5/3) |Disc|^(1/6) sum over the three pairings = (theta2^8 + theta3^8 + theta4^8)/(2 eta^8) = j^(1/3)",
      worst < mp.mpf(10)**(-20), mp.nstr(worst, 3))
check("closed form: chi(delta) = 4 (4 - sin^2 delta) / sin^(4/3) delta",
      all(abs(chi_from_matrix(d)-4*(4-mp.sin(d)**2)/mp.sin(d)**(mp.mpf(4)/3)) < mp.mpf(10)**(-22)
          for d in (mp.mpf("0.4"), mp.mpf("1.2"))))
check("at the square: 12, the recorded value chi_E8(i) = 12", abs(chi_from_matrix(PI/2)-12) < mp.mpf(10)**(-25))
dcw = mp.atan(mp.mpf(4)/3)
check("at the counterwitness: chi^3 = 148176/25, the recorded j of v512",
      abs(chi_from_matrix(dcw)**3-mp.mpf(148176)/25) < mp.mpf(10)**(-20), mp.nstr(chi_from_matrix(dcw)**3, 12))
s_ = sp.Rational(4, 5)
check("exactly: (4 (4 - 16/25) (5/4)^(4/3))^3 = 148176/25",
      sp.simplify((4*(4-s_**2))**3/s_**4-sp.Rational(148176, 25)) == 0)

# ---------------------------------------------------------------- Z3
print("Z3  one pairing against the sum")
single = weight(dcw, "opposite")/weight(PI/2, "opposite")
summed = sum(weight(dcw, k) for k in ("opposite", "side", "side'"))/sum(weight(PI/2, k) for k in ("opposite", "side", "side'"))
check("moved marks against the square: the recorded spin structure alone 25/16, the sum of the three 7/4",
      abs(single-mp.mpf(25)/16) < mp.mpf(10)**(-25) and abs(summed-mp.mpf(7)/4) < mp.mpf(10)**(-25))
d_ = sp.symbols("delta", positive=True)
check("closed forms: 1/sin^2(delta) and (4 - sin^2 delta)/(3 sin^2 delta)",
      abs(summed-(4-mp.sin(dcw)**2)/(3*mp.sin(dcw)**2)) < mp.mpf(10)**(-25))
print("      sixteen free fermions with one spin structure: one eigenvalue.  The E8 net: the sum over the three.")
print("      This is the open arrow of SEAM.EQUIV.01 (free reflection-positive seam => single-sector E8 net) as one number.")

# ---------------------------------------------------------------- Z4
print("Z4  the recorded mark determinant in these terms")
eps = mp.mpf("0.37")
c3 = 1/(8*PI)
d = mp.mpf("1.1")
lhs = mp.det(mp.eye(4)+eps/4*mark_matrix(d))
t2, t3, t4, eta = thetas(tau_of(d))
lnDisc = 2*mp.log(16*mp.sin(d)**2)
dets = [8*mp.log(abs(th/eta))-lnDisc/6-mp.mpf(8)/3*mp.log(2) for th in (t2, t3, t4)]
rhs = (1-eps*c3*lnDisc/2)*mp.fprod([1-eps*c3*L for L in dets])
check("det(I + (eps/4) G_marks) = (1 - eps c3 ln|Disc|/2) * prod over the three spin structures of "
      "(1 - eps c3 [8 ln|theta/eta| - ln|Disc|/6 - (8/3) ln 2])", abs(lhs-rhs) < mp.mpf(10)**(-22), mp.nstr(abs(lhs-rhs), 3))

# ---------------------------------------------------------------- Z5
print("Z5  the three pairings along the recorded family X^4 + z X + q = 0")
mp.mp.dps = 25
cq = mp.mpf(3)/mp.mpf(4)**(mp.mpf(4)/3)
xi = mp.root(cq, 4)
roots0 = [xi*mp.mpc(0, 1)**k for k in range(4)]   # z = 0, q = -cq


def follow(path, roots, qstart, steps=3000):
    cur, qc = list(roots), qstart
    for n in range(1, steps+1):
        zz = path(mp.mpf(n)/steps)
        cands = [cq*mp.exp((mp.log(zz**4-1)+2j*PI*m)/3) for m in (-2, -1, 0, 1, 2)]
        qc = min(cands, key=lambda v: abs(v-qc))
        new = []
        for r in cur:
            for _ in range(8):
                r = r-(r**4+zz*r+qc)/(4*r**3+zz)
            new.append(r)
        cur = new
    return cur, qc


zb = mp.mpc("0.4", "0.15")                        # base point away from the square, where two pairings coincide
rootsb, qb = follow(lambda tt: zb*tt, roots0, -cq, steps=600)


def loop(tt):
    c = mp.mpf(3)/4
    if tt < mp.mpf("0.3"):
        return zb+(c-zb)*tt/mp.mpf("0.3")
    if tt < mp.mpf("0.7"):
        return 1+(c-1)*mp.exp(2j*PI*(tt-mp.mpf("0.3"))/mp.mpf("0.4"))
    return c+(zb-c)*(tt-mp.mpf("0.7"))/mp.mpf("0.3")


def pair_sq(r):
    return [((r[0]-r[1])*(r[2]-r[3]))**2, ((r[0]-r[2])*(r[1]-r[3]))**2, ((r[0]-r[3])*(r[1]-r[2]))**2]


end, qend = follow(loop, rootsb, qb)
before, after = pair_sq(rootsb), pair_sq(end)
fac = qend/qb
perm = [min(range(3), key=lambda j: abs(after[i]-fac*before[j])) for i in range(3)]
check("after the loop around the parameter mark the three squared pairing products are the old ones, "
      "cyclically permuted, times one common cube root of unity",
      sorted(perm) == [0, 1, 2] and all(perm[i] != i for i in range(3))
      and max(abs(after[i]-fac*before[perm[i]]) for i in range(3)) < mp.mpf(10)**(-12)
      and abs(fac**3-1) < mp.mpf(10)**(-12), perm)
check("so a combination with constant weights returns to itself up to that phase only if the three weights "
      "are equal (the E8 sum) or differ by cube roots of unity",
      abs(sum(after)-fac*sum(before)) < mp.mpf(10)**(-12))
check("the equal-weight sum is 24 q along the family: chi_E8 = 12 (1 - z^4)^(1/3) up to the phase of the cube root",
      abs(sum(before)-24*qb) < mp.mpf(10)**(-12) and abs(sum(pair_sq(roots0))-24*(-cq)) < mp.mpf(10)**(-15))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
