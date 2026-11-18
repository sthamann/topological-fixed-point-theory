#!/usr/bin/env python3
"""The recorded free elliptic source against the level-two Gaussian of the marks (18 November 2026).

Recorded, unchanged:
    experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/
        ELLIPTIC-LOCAL-INSERTION-PROOF.txt, section 4 (R066):
        b = 4 pi, independent quarter charges q, r, F_1/4(q) + F_3/4(r) = (p^2 + d^2 - d)/4,
        P(p) = exp(-pi p^2)/sum, Var(Q) = 1/(8 pi), means +-0.04134328958148
    Atlas finding actual-elliptic-gaussian-source-determines-full-quarter-current-functional:
        C_r = (1 + exp(4 pi (r + alpha)))^-1, rho = exp(-4 pi H_alpha)/Z, H_alpha = N + Q^2/2 + (alpha - 1/2) Q
    experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/MODULAR-C3-SOURCE-TEST.txt:
        E_i = C/(Z + iZ), Q: z -> iz, half periods a = 1/2, b = i/2, c = (1+i)/2, Q(a) = b, Q(c) = c

Sections
    R1  the recorded charge law is the level-two Gaussian with shifts -+1/4
    R2  recorded numbers reproduced (means, variances, independence, total law)
    R3  fermion occupation, charge law and the quarter blocks (Jacobi triple product)
    R4  the two tori: E (branched at the four marks) and the recorded E_i (branched at 0, inf, +-1 of z^2)
    R5  ledger combinations against the cumulants of the quarter charge
    R6  uniqueness over all offsets and widths
    R7  the recorded kernel of the free source at the marks
"""
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- R1
print("R1  recorded law = level-two Gaussian")
q, r, p, d = sp.symbols("q r p d")
F = lambda al, Q: Q**2/2+(al-sp.Rational(1, 2))*Q
rec = F(sp.Rational(1, 4), q)+F(sp.Rational(3, 4), r)
check("F_1/4(q) + F_3/4(r) = (p^2 + d^2 - d)/4 with p = q + r, d = q - r",
      sp.simplify(rec.subs({q: (p+d)/2, r: (p-d)/2})-(p**2+d**2-d)/4) == 0)
gauss = 2*((q-sp.Rational(1, 4))**2+(r+sp.Rational(1, 4))**2)
check("4 pi (F_1/4 + F_3/4) = 2 pi [(q - 1/4)^2 + (r + 1/4)^2] - pi/4: the weight exp(-2 pi |v|^2), shifts -+1/4",
      sp.simplify(4*rec-gauss+sp.Rational(1, 4)) == 0)

# ---------------------------------------------------------------- R2
print("R2  recorded numbers")


def law(al):
    w = lambda n: mp.exp(-4*PI*(n**2/2+(al-mp.mpf(1)/2)*n))
    Z = mp.nsum(w, [-mp.inf, mp.inf])
    m1 = mp.nsum(lambda n: n*w(n), [-mp.inf, mp.inf])/Z
    m2 = mp.nsum(lambda n: n*n*w(n), [-mp.inf, mp.inf])/Z
    return m1, m2-m1**2


mq, vq = law(mp.mpf(1)/4)
mr, vr = law(mp.mpf(3)/4)
check("means +-0.04134328958148 (recorded value)",
      abs(mq-mp.mpf("0.04134328958148")) < mp.mpf(10)**(-13) and abs(mr+mq) < mp.mpf(10)**(-30), mp.nstr(mq, 16))
check("both variances equal 1/(8 pi)", abs(vq-1/(8*PI)) < mp.mpf(10)**(-36) and abs(vr-vq) < mp.mpf(10)**(-36))
Zt = lambda t: mp.nsum(lambda n: mp.exp(-PI*t*n*n), [-mp.inf, mp.inf])
check("recorded proof: Z(t) = t^-1/2 Z(1/t), (log Z)'(1) = -1/4, so Var(p) = 1/(4 pi)",
      abs(Zt(mp.mpf("1.7"))-Zt(1/mp.mpf("1.7"))/mp.sqrt(mp.mpf("1.7"))) < mp.mpf(10)**(-30)
      and abs(mp.diff(lambda t: mp.log(Zt(t)), 1)+mp.mpf(1)/4) < mp.mpf(10)**(-25))
tot = {}
for N in (-2, -1, 0, 1, 2, 3):
    tot[N] = mp.nsum(lambda n: mp.exp(-4*PI*(n**2/2-n/4+(N-n)**2/2+(N-n)/4)), [-mp.inf, mp.inf])
check("total charge law exp(-pi p^2), discrete, for even and odd p alike",
      all(abs(tot[N]/tot[0]-mp.exp(-PI*N*N)) < mp.mpf(10)**(-30) for N in tot))

# ---------------------------------------------------------------- R3
print("R3  fermion blocks and charge law")
x = mp.exp(-4*PI)
for al in (mp.mpf(1)/4, mp.mpf(3)/4):
    yv = mp.mpf("0.83")
    lhs = mp.nprod(lambda n: (1+yv*x**(n+al))*(1+x**(n+1-al)/yv), [0, mp.inf])
    rhs = mp.nsum(lambda Q: yv**Q*x**(Q*Q/2+(al-mp.mpf(1)/2)*Q), [-mp.inf, mp.inf])/mp.nprod(lambda n: 1-x**n, [1, mp.inf])
    ok = abs(lhs-rhs) < mp.mpf(10)**(-30)
    check("modes r + %s with occupation (1 + exp(4 pi (r + alpha)))^-1 carry the charge law exp(-4 pi F_alpha)" % mp.nstr(al, 3), ok)
ks = [mp.mpf(n)+mp.mpf(1)/2 for n in range(-4, 4)]
eig = [mp.exp(1j*PI*k) for k in ks]
check("on the spinor modes k in Z + 1/2 of the unit square torus the translation by the half period 1/2 has "
      "eigenvalues +-i; k in 2Z + 1/2 is alpha = 1/4, k in 2Z + 3/2 is alpha = 3/4",
      all(abs(e-1j) < mp.mpf(10)**(-30) or abs(e+1j) < mp.mpf(10)**(-30) for e in eig)
      and all(abs(mp.exp(-2*PI*k)-mp.exp(-4*PI*(mp.floor(k/2)+(k/2-mp.floor(k/2))))) < mp.mpf(10)**(-30) for k in ks))

# ---------------------------------------------------------------- R4
print("R4  the two tori")
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI))
sl = lambda u: mp.ellipfun("sn", u, m=-1)
W = lambda u: sl(u)**2                              # function of z^2
u0 = mp.mpc("0.31", "0.17")
check("W = z^2 = sl^2 has the larger lattice varpi (Z + iZ): a square torus E_i, of which E is a double cover",
      abs(W(u0+varpi)-W(u0)) < mp.mpf(10)**(-25) and abs(W(u0+1j*varpi)-W(u0)) < mp.mpf(10)**(-25)
      and abs(sl(u0+varpi)+sl(u0)) < mp.mpf(10)**(-25))
check("on E_i (lattice Z + iZ): 0 is z = 0, the half periods a = 1/2 and b = i/2 are z^2 = 1 and z^2 = -1 "
      "(the marks +-1 and +-i), c = (1+i)/2 is z = infinity",
      abs(W(varpi/2)-1) < mp.mpf(10)**(-25) and abs(W(1j*varpi/2)+1) < mp.mpf(10)**(-25)
      and abs(1/W((1+1j)*varpi/2+mp.mpf(10)**(-15))) < mp.mpf(10)**(-20) and abs(W(mp.mpf(10)**(-15))) < mp.mpf(10)**(-20))
check("the quarter turn fixes 0 and c and exchanges a and b: the recorded action Q = S on the half periods",
      abs(W(1j*u0)+W(u0)) < mp.mpf(10)**(-25))
# level one on E_i at the half period versus level two on E at the mark
lvl1 = mp.nsum(lambda P, D: mp.exp(-PI*(P*P+(D-mp.mpf(1)/2)**2)) if (int(P)-int(D)) % 2 == 0 else 0,
               [-12, 12], [-12, 12])
lvl2 = mp.nsum(lambda a, b: mp.exp(-2*PI*((a-mp.mpf(1)/4)**2+(b+mp.mpf(1)/4)**2)), [-12, 12], [-12, 12])
check("self-dual Gaussian exp(-pi |w|^2) of E_i over the lifts of the half period, one sheet of E "
      "= Gaussian exp(-2 pi |v|^2) of E over the lifts of a mark (w = (1+i) v)",
      abs(lvl1-lvl2) < mp.mpf(10)**(-30))
t2 = mp.nsum(lambda n: mp.exp(-PI*(n+mp.mpf(1)/2)**2), [-mp.inf, mp.inf])
t3 = mp.nsum(lambda n: mp.exp(-PI*n*n), [-mp.inf, mp.inf])
check("without the sheet condition the sum is theta2(i) theta3(i); one sheet carries exactly half",
      abs(2*lvl1-t2*t3) < mp.mpf(10)**(-30))
print("      dictionary: p = coordinate of the lift along the direction in which the mark is a lattice point,")
print("      d - 1/2 = coordinate along the direction in which it is a half period; q, r = the diagonals.")

# ---------------------------------------------------------------- R5
print("R5  ledger combinations against the cumulants (a control, expected negative)")
c3 = 1/(8*PI)


def cumulants(al):
    w = lambda n: mp.exp(-4*PI*(n**2/2+(al-mp.mpf(1)/2)*n))
    Z = mp.nsum(w, [-mp.inf, mp.inf])
    m = [mp.nsum(lambda n, k=k: n**k*w(n), [-mp.inf, mp.inf])/Z for k in range(1, 5)]
    k1 = m[0]
    k2 = m[1]-m[0]**2
    k3 = m[2]-3*m[1]*m[0]+2*m[0]**3
    k4 = m[3]-4*m[2]*m[0]-3*m[1]**2+12*m[1]*m[0]**2-6*m[0]**4
    return k1, k2, k3, k4


k1, k2, k3, k4 = cumulants(mp.mpf(1)/4)
targets = {"c3 + 36 c3^4": c3+36*c3**4, "phi0 = (4/3)(c3 + 36 c3^4)": mp.mpf(4)/3*(c3+36*c3**4)}
cands = {"mean": abs(k1), "third cumulant": abs(k3), "fourth cumulant": abs(k4), "mean/variance": abs(k1)/k2,
         "variance + |fourth|": k2+abs(k4)}
print("      cumulants:", [mp.nstr(v, 10) for v in (k1, k2, k3, k4)])
hits = [(a, b) for a, va in cands.items() for b, vb in targets.items() if abs(va-vb) < mp.mpf(10)**(-6)]
check("no simple cumulant equals c3 + 36 c3^4 or phi0 (no further ledger number comes for free)", not hits, hits)

# ---------------------------------------------------------------- R6
print("R6  uniqueness over all offsets (a1, a2) and widths beta")


def parity_gap(delta, beta):
    """sum over D of (-1)^D exp(-(beta/4)(D + delta)^2): zero iff even and odd totals share one law"""
    return mp.nsum(lambda D: (-1)**int(D)*mp.exp(-beta/4*(D+delta)**2), [-mp.inf, mp.inf])


zeros = {}
for beta in (mp.mpf(2), 4*PI, mp.mpf(30)):
    grid = [mp.mpf(k)/200 for k in range(0, 201)]
    vals = [parity_gap(g, beta) for g in grid]
    zs = [grid[k] for k in range(200) if vals[k]*vals[k+1] < 0 or vals[k] == 0]
    zeros[mp.nstr(beta, 4)] = [mp.nstr(mp.findroot(lambda t: parity_gap(t, beta), z0+mp.mpf(1)/400), 10) for z0 in zs]
check("the parity gap vanishes exactly at a1 - a2 = 1/2 (mod 1) for every width: zeros of theta2",
      all(v == ["0.5"] for v in zeros.values()), zeros)
dual_dev = {}
thq = mp.mpf("0.3")
per = lambda t: mp.nsum(lambda m: mp.exp(-PI*(m-t)**2), [-mp.inf, mp.inf])
for sft in ("0", "0.1", "0.25", "0.5"):
    sv = mp.mpf(sft)
    Zs = mp.nsum(lambda N: mp.exp(-PI*(N+sv)**2), [-mp.inf, mp.inf])
    char = mp.nsum(lambda N: mp.exp(-PI*(N+sv)**2)*mp.exp(2j*PI*N*thq), [-mp.inf, mp.inf])/Zs
    dual_dev[sft] = abs(char-per(thq)/per(0))
check("the characteristic function of the total is the same Gaussian, periodised, only for a1 + a2 in Z",
      dual_dev["0"] < mp.mpf(10)**(-30) and all(dual_dev[k] > mp.mpf(10)**(-3) for k in ("0.1", "0.25", "0.5")),
      {k: mp.nstr(v, 3) for k, v in dual_dev.items()})
widths = {mp.nstr(b, 4): abs(mp.diff(lambda t: mp.log(mp.nsum(lambda N: mp.exp(-b/4*t*N*N), [-mp.inf, mp.inf])), 1)
                           + mp.mpf(1)/4) for b in (mp.mpf(8), 4*PI, mp.mpf(16))}
check("self-dual width only at beta = 4 pi", widths[mp.nstr(4*PI, 4)] < mp.mpf(10)**(-25)
      and all(v > mp.mpf(10)**(-4) for k, v in widths.items() if k != mp.nstr(4*PI, 4)), {k: mp.nstr(v, 3) for k, v in widths.items()})
print("      together: (a1, a2) = (+-1/4, -+1/4) and beta = 4 pi; the other two marks are the same statement")
print("      with total and relative charge exchanged.")

# ---------------------------------------------------------------- R7
print("R7  the recorded Szego kernel of the spin structure fixed by the quarter turn, at the marks")
tau = 1j


def th(z, kind):
    N = 40
    if kind == 1:
        return 2*mp.nsum(lambda n: (-1)**int(n)*mp.exp(1j*PI*tau*(n+mp.mpf(1)/2)**2)*mp.sin((2*n+1)*z), [0, N])
    return mp.nsum(lambda n: mp.exp(1j*PI*tau*n*n+2j*n*z), [-N, N])


th1p = mp.diff(lambda z: th(z, 1), 0)
S = lambda u: PI*th(PI*u, 3)*th1p/(th(0, 3)*th(PI*u, 1))    # Szego kernel, pole 1/u at 0, lattice Z + iZ
ua, ub, uc = mp.mpf(1)/2, 1j/2, (1+1j)/2
check("simple pole with residue one at the centre of the quarter turn", abs(S(mp.mpf(10)**(-8))*mp.mpf(10)**(-8)-1) < mp.mpf(10)**(-10))
check("S(i u) = -i S(u): from the centre, the two mark classes are seen with a relative quarter phase",
      abs(S(1j*mp.mpc("0.2", "0.1"))+1j*S(mp.mpc("0.2", "0.1"))) < mp.mpf(10)**(-25) and abs(S(ub)+1j*S(ua)) < mp.mpf(10)**(-25))
check("S vanishes at c = a - b: the kernel between neighbouring marks (+-1 against +-i) is exactly zero",
      abs(S(uc)) < mp.mpf(10)**(-25))
check("S(a)^2 = e_a - e_c = varpi^2 = Gamma(1/4)^4/(8 pi)", abs(S(ua)**2-varpi**2) < mp.mpf(10)**(-20), mp.nstr(S(ua)**2, 12))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
