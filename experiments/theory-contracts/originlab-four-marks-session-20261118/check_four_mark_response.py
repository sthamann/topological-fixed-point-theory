#!/usr/bin/env python3
"""Tenth pasted note: the four-mark response of the free source, and the phase step at the marks
(18 November 2026).

Sections
    F1  marks {1, -1, u, -u}: the three pairing contributions P^2/sqrt(Disc) and their half sum I/sqrt(Disc)
    F2  the table of the note at u = (3+4i)/5 against the square u = i
    F3  the phase step at the marks on the original circle: first mode, logarithmic growth
    F4  overlap with the unmarked circle: power law in the number of sites
    F5  the recorded seam covariance of v512 (H2.1) at the four marks against the three pairing responses
"""
import sys

import numpy as np
import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- F1
print("F1  three pairings and their half sum")
u, X = sp.symbols("u X")
marks = [1, -1, u, -u]
P = {"opposite {1,-1}{u,-u}": (marks[0]-marks[1])*(marks[2]-marks[3]),
     "{1,u}{-1,-u}": (marks[0]-marks[2])*(marks[1]-marks[3]),
     "{1,-u}{-1,u}": (marks[0]-marks[3])*(marks[1]-marks[2])}
poly = sp.Poly(sp.expand(sp.prod([X-m for m in marks])), X)
c4, c3, c2, c1, c0 = poly.all_coeffs()
I = sp.expand(12*c4*c0-3*c3*c1+c2**2)
check("I = half the sum of the squared pairing products", sp.expand(I-sum(p**2 for p in P.values())/2) == 0)
check("I = 1 + 14 u^2 + u^4 for the marks {1, -1, u, -u}", sp.expand(I-(1+14*u**2+u**4)) == 0)
disc = sp.prod([(marks[i]-marks[j])**2 for i in range(4) for j in range(i+1, 4)])
rootdisc = 4*u*(1-u**2)**2
check("Disc = (4 u (1 - u^2)^2)^2", sp.expand(disc-rootdisc**2) == 0)

# ---------------------------------------------------------------- F2
print("F2  the table")
uv = sp.Rational(3, 5)+sp.Rational(4, 5)*sp.I
mod = lambda e, val: sp.simplify(sp.Abs(e.subs(u, val)))
single = {k: sp.simplify(mod(p**2/rootdisc, uv)/mod(p**2/rootdisc, sp.I)) for k, p in P.items()}
for k, vv in single.items():
    print("      single pairing %-24s moved marks / square = %s" % (k, vv))
check("recorded spin structure (the pairing of opposite marks): 25/16", single["opposite {1,-1}{u,-u}"] == sp.Rational(25, 16))
tot = sp.simplify(mod(I/rootdisc, uv)/mod(I/rootdisc, sp.I))
check("half sum of the three: 7/4", tot == sp.Rational(7, 4), tot)
check("the other two single pairings give 1/4 and 4", sorted([single["{1,u}{-1,-u}"], single["{1,-u}{-1,u}"]]) == [sp.Rational(1, 4), 4])
check("at the square the half sum has modulus 3/4 of the recorded single pairing",
      sp.simplify(mod(I/rootdisc, sp.I)/mod(P["opposite {1,-1}{u,-u}"]**2/rootdisc, sp.I)) == sp.Rational(3, 4))

# ---------------------------------------------------------------- F3
print("F3  the phase step on the original circle")
# g = 1 on two opposite quarter arcs, i on the other two: g = a + b s, s = sign(sin 2 theta)
a_, b_ = (1+1j)/2, (1-1j)/2


def shat(n):
    if n % 2 or (n//2) % 2 == 0:
        return 0.0
    return 2/(1j*np.pi*(n//2))


first = sum(abs(b_*shat(n))**2 for n in range(1, 400001))
check("first outer mode: a quarter of its norm square ends up in the inner modes", abs(first-0.25) < 1e-5, round(first, 6))
tail = {}
for K in (64, 256, 1024, 4096):
    tail[K] = sum(min(n, K)*abs(b_*shat(n))**2 for n in range(1, 200001))
slopes = [(tail[4*K]-tail[K])/np.log(4) for K in (64, 256, 1024)]
check("summed over the outer modes the weight grows like (2/pi^2) log K: no bound, the step is not implementable",
      all(abs(s-2/np.pi**2) < 5e-3 for s in slopes), [round(s, 5) for s in slopes])

# ---------------------------------------------------------------- F4
print("F4  overlap with the unmarked circle on N sites")


def toeplitz_logdet(N):
    idx = np.arange(N)
    d = idx[:, None]-idx[None, :]
    T = np.zeros((N, N), complex)
    T[d == 0] = a_
    vals = {n: b_*shat(n) for n in range(-N, N+1) if n != 0}
    for n, vv in vals.items():
        if vv != 0:
            T[d == n] = vv
    return np.linalg.slogdet(T)[1]


Ns = [16, 32, 64, 128, 256, 512]
ld = {N: toeplitz_logdet(N) for N in Ns}
sl = [(ld[2*N]-ld[N])/np.log(2) for N in Ns[:-1]]
print("      log |D_N|:", {N: round(v, 5) for N, v in ld.items()})
check("the overlap falls like N^(-1/4): four steps of a quarter phase, 4 * (1/4)^2",
      abs(sl[-1]+0.25) < 5e-3, [round(s, 4) for s in sl])
check("both sectors together N^(-1/2) per spinor, N^(-4) for eight spinors: weight 1/2 per mark, "
      "as in the sixteenth round", 2*sp.Rational(1, 4) == sp.Rational(1, 2) and 8*sp.Rational(1, 2) == 4)
print("      at the sixteen sites of v519 the single-sector overlap is |D_16| = %.6f" % np.exp(ld[16]))

# ---------------------------------------------------------------- F5
print("F5  the recorded 4 x 4 seam covariance of v512")
# verification/v512_seam_tau_flag.py, H2.1: marks at angles 0, delta, pi, pi + delta;
# G(theta) = -(1/pi) log(2 sin(theta/2)), regularised diagonal c left free.
dl, c = sp.symbols("delta c", positive=True)
uu = -sp.log(2*sp.sin(dl/2))/sp.pi
vv = -sp.log(2*sp.cos(dl/2))/sp.pi
ww = -sp.log(2)/sp.pi
G4 = sp.Matrix([[c, uu, ww, vv], [uu, c, vv, ww], [ww, vv, c, uu], [vv, ww, uu, c]])
chars = {"symmetric": sp.Matrix([1, 1, 1, 1]), "opposite": sp.Matrix([1, -1, 1, -1]),
         "side A": sp.Matrix([1, -1, -1, 1]), "side B": sp.Matrix([1, 1, -1, -1])}
lam = {k: sp.simplify((G4*vec)[0]/vec[0]) for k, vec in chars.items()}
check("the four characters of the pairings are eigenvectors (as recorded)",
      all(sp.simplify(G4*vec-lam[k]*vec) == sp.zeros(4, 1) for k, vec in chars.items()))
resp = {"opposite": 1/sp.sin(dl)**2, "side A": sp.cot(dl/2)**2/4, "side B": sp.tan(dl/2)**2/4}
ud = sp.exp(sp.I*dl)
free = {"opposite": P["opposite {1,-1}{u,-u}"]**2/rootdisc, "side B": P["{1,u}{-1,-u}"]**2/rootdisc,
        "side A": P["{1,-u}{-1,u}"]**2/rootdisc}
num_ok = True
for k in resp:
    for dval in (0.4, 0.9273, 1.3):
        fv = abs(complex(free[k].subs(u, sp.exp(sp.I*dval)).evalf()))
        rv = float(resp[k].subs(dl, dval))
        num_ok &= abs(fv-rv) < 1e-9
check("moduli of the three single pairing responses: 1/sin^2(delta), cot^2(delta/2)/4, tan^2(delta/2)/4", num_ok)
ident = True
for k in resp:
    for dval in (0.4, 0.9273, 1.3):
        lhs = float(sp.exp(-2*sp.pi*lam[k].subs(c, 0)).subs(dl, dval).evalf())
        ident &= abs(lhs-float(resp[k].subs(dl, dval))) < 1e-9
check("exp(-2 pi (eigenvalue - c)) of the recorded covariance = the three free pairing responses, for every delta",
      ident)
symv = float(sp.exp(-2*sp.pi*lam["symmetric"].subs(c, 0)).subs(dl, 0.9273).evalf())
check("on the symmetric vector it is |Disc|^(1/2) = 16 sin^2(delta)",
      abs(symv-16*np.sin(0.9273)**2) < 1e-9)
split = sp.simplify(lam["side B"]-lam["side A"])
check("recorded odd-doublet split (2/pi) log cot(delta/2) = (1/2 pi) log of the ratio of the two side responses",
      abs(float((split-2*sp.log(sp.cot(dl/2))/sp.pi).subs(dl, 0.7).evalf())) < 1e-12)
dcw = sp.atan(sp.Rational(4, 3))
r_opp = sp.simplify(sp.exp(-2*sp.pi*(lam["opposite"].subs(dl, dcw)-lam["opposite"].subs(dl, sp.pi/2))))
check("recorded covariance, opposite pairing, counterwitness against square: exactly 25/16",
      abs(float(r_opp.evalf())-25/16) < 1e-12, sp.nsimplify(r_opp.evalf(), rational=True))
print("      the free diagonal c of the record drops out of every ratio; it is the missing absolute reference.")

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
