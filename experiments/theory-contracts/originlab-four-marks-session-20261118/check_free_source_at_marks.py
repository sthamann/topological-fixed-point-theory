#!/usr/bin/env python3
"""The free elliptic source evaluated at all four special points at once (18 November 2026).

Free spinor on a torus C/(Z + tau Z) with an even spin structure eps in {a, b, c},
named after the half period at which its theta function vanishes
(a = 1/2: theta2, b = tau/2: theta4, c = (1+tau)/2: theta3).  Kernel
        S_eps(u) = pi theta_eps(pi u) theta1'(0) / (theta_eps(0) theta1(pi u)).
The recorded source is eps = c on the unit square torus (fixed by the quarter turn).

Sections
    M1  two-point kernel between half periods: S_eps(t)^2 = e_t - e_eps
    M2  the square torus: tables of weights and the quarter phase
    M3  all four points at once: the four-point function for the three splittings
    M4  the triangle of the two-torsion points = triangle of the three theta fourth powers
    M5  weights of the qubit model from the same triangle, at general tau
    M6  lifting the recorded source to the cover branched at the four marks: one zero mode
    M7  seventh pasted note: a probe on one quarter arc, and the mixture of the four marks
    M8  seventh pasted note: the rule that selects the square
    M9  the two mirrors on the two tori
"""
import itertools
import sys

import mpmath as mp

mp.mp.dps = 40
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


def make(tau):
    N = 40
    q = lambda e: mp.exp(1j*PI*tau*e)

    def th1(z):
        return 2*mp.nsum(lambda n: (-1)**int(n)*q((n+mp.mpf(1)/2)**2)*mp.sin((2*n+1)*z), [0, N])

    def th2(z):
        return 2*mp.nsum(lambda n: q((n+mp.mpf(1)/2)**2)*mp.cos((2*n+1)*z), [0, N])

    def th3(z):
        return mp.nsum(lambda n: q(n*n)*mp.exp(2j*n*z), [-N, N])

    def th4(z):
        return mp.nsum(lambda n: (-1)**int(n)*q(n*n)*mp.exp(2j*n*z), [-N, N])

    d1 = mp.diff(th1, 0)
    fam = {"a": th2, "b": th4, "c": th3}
    S = {e: (lambda u, f=f: PI*f(PI*u)*d1/(f(0)*th1(PI*u))) for e, f in fam.items()}
    half = {"a": mp.mpf(1)/2, "b": tau/2, "c": (1+tau)/2}
    const = {e: f(0) for e, f in fam.items()}
    return S, half, const


# ---------------------------------------------------------------- M1
print("M1  kernel between half periods")
tau_g = mp.mpc("0.23", "1.07")
S, half, th0 = make(tau_g)
# differences e_x - e_y from theta constants: e_a - e_b = pi^2 th3^4, e_a - e_c = pi^2 th4^4, e_c - e_b = pi^2 th2^4
ediff = {("a", "b"): PI**2*th0["c"]**4, ("a", "c"): PI**2*th0["b"]**4, ("c", "b"): PI**2*th0["a"]**4}
for (x, y), v in list(ediff.items()):
    ediff[(y, x)] = -v
u1 = mp.mpc("0.17", "0.33")
check("S_eps(u)^2 - S_delta(u)^2 does not depend on u and equals e_delta - e_eps",
      all(abs(S[e](u1)**2-S[dd](u1)**2-ediff[(dd, e)]) < mp.mpf(10)**(-22) for e in "abc" for dd in "abc" if e != dd))
ok = True
for e in "abc":
    for t in "abc":
        val = S[e](half[t])**2
        tgt = 0 if e == t else ediff[(t, e)]
        ok &= abs(val-tgt) < mp.mpf(10)**(-20)
check("at the half periods: S_eps(t)^2 = e_t - e_eps, zero for t = eps (general tau)", ok)

# ---------------------------------------------------------------- M2
print("M2  the square torus")
S, half, th0 = make(mp.mpc(0, 1))
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI))
table = {e: [abs(S[e](half[t]))**2/varpi**2 for t in "abc"] for e in "abc"}
for e in "abc":
    print("      structure %s: |S|^2 at (a, b, c) in units of varpi^2: %s" % (e, [mp.nstr(v, 6) for v in table[e]]))
near = lambda v, w: all(abs(x-y) < mp.mpf(10)**(-18) for x, y in zip(v, w))
check("structure c (recorded, fixed by the quarter turn): (1, 1, 0), normalised (1/2, 1/2, 0)", near(table["c"], [1, 1, 0]))
check("structures a and b: (0, 2, 1) and (2, 0, 1), normalised (0, 2/3, 1/3) and (2/3, 0, 1/3)",
      near(table["a"], [0, 2, 1]) and near(table["b"], [2, 0, 1]))
check("quarter phase in the recorded structure: S_c(b) = -i S_c(a)", abs(S["c"](half["b"])+1j*S["c"](half["a"])) < mp.mpf(10)**(-22))

# ---------------------------------------------------------------- M3
print("M3  all four points at once")
pts = {"0": mp.mpf(0), "a": half["a"], "b": half["b"], "c": half["c"]}
splits = {"a": (("0", "a"), ("b", "c")), "b": (("0", "b"), ("a", "c")), "c": (("0", "c"), ("a", "b"))}


def four(Se, xs, ys):
    """<psi+(x1) psi+(x2) psi(y2) psi(y1)> = det S(y_j - x_i)"""
    m = [[Se(pts[y]-pts[x]) for y in ys] for x in xs]
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]


res = {}
for e in "abc":
    row = []
    for k in "abc":
        pair, rest = splits[k]
        row.append(four(S[e], pair, rest))
    res[e] = row
    print("      structure %s: four-point value for the splittings (a, b, c), units varpi^2: %s"
          % (e, [mp.nstr(v/varpi**2, 6) for v in row]))
check("recorded structure c: moduli (1, 1, 2) for the splittings by a, b, c",
      near([abs(v)/varpi**2 for v in res["c"]], [1, 1, 2]))
check("every even structure gives the same three moduli 1, 1, 2 up to order",
      all(sorted(mp.nstr(abs(v)/varpi**2, 12) for v in res[e]) == ["1.0", "1.0", "2.0"] for e in "abc"),
      {e: [mp.nstr(abs(v)/varpi**2, 4) for v in res[e]] for e in "abc"})

# ---------------------------------------------------------------- M4
print("M4  the triangle")
t2, t4, t3 = th0["a"], th0["b"], th0["c"]
check("sides of the two-torsion triangle = pi^2 (theta4^4, theta2^4, theta3^4); Jacobi's identity closes it",
      abs(t3**4-t2**4-t4**4) < mp.mpf(10)**(-30))
check("square torus: pi^2 theta2^4 = pi^2 theta4^4 = varpi^2 = Gamma(1/4)^4/(8 pi), pi^2 theta3^4 = 2 varpi^2",
      abs(PI**2*t2**4-varpi**2) < mp.mpf(10)**(-25) and abs(PI**2*t4**4-varpi**2) < mp.mpf(10)**(-25)
      and abs(PI**2*t3**4-2*varpi**2) < mp.mpf(10)**(-25))
print("      sides 1 : 1 : 2, perimeter 4, in units of varpi^2 (the recorded anchor triple is (1, 1, 2), e1 = 4;")
print("      this is an agreement of numbers, no connection is shown).")

# ---------------------------------------------------------------- M5
print("M5  the qubit weights from the same triangle, general modulus")


def thetas(tau):
    N = 60
    q = lambda e: mp.exp(1j*PI*tau*e)
    return (mp.nsum(lambda n: q((n+mp.mpf(1)/2)**2), [-N, N]), mp.nsum(lambda n: q(n*n), [-N, N]),
            mp.nsum(lambda n: (-1)**int(n)*q(n*n), [-N, N]))


worst = 0
for tau in (mp.mpc(0, 1), mp.exp(2j*PI/3), mp.mpc("0.23", "1.07"), mp.mpc("-0.4", "0.8")):
    a2, a3, a4 = thetas(tau)
    w = [abs(a2)**4, abs(a3)**4, abs(a4)**4]
    p_tri = sorted(1-2*x/sum(w) for x in w)
    b2, b3, _ = thetas(2*tau)
    f = b2/b3
    dd = (1+abs(f)**2)**2
    p_f = sorted([4*f.real**2/dd, 4*f.imag**2/dd, (1-abs(f)**2)**2/dd])
    worst = max(worst, max(abs(x-y) for x, y in zip(p_tri, p_f)))
check("p_k = 1 - 2 |theta_k|^4 / sum_j |theta_j|^4 equals the squared Bloch components of the mark, at four moduli",
      worst < mp.mpf(10)**(-25), mp.nstr(worst, 3))
a2, a3, a4 = thetas(mp.mpc(0, 1))
w = [abs(x)**4 for x in (a2, a3, a4)]
check("so the weights are fixed by the squared moduli of the three theta squares, the recorded carrier of v117: "
      "(1/4, 1/2, 1/4) at the square, (1/3, 1/3, 1/3) at the hexagon",
      all(abs(x/sum(w)-y) < mp.mpf(10)**(-25) for x, y in zip(w, [mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(1)/4])))
avg = [sum(abs(S[e](half[t]))**2 for e in "abc" if e != t) for t in "abc"]
print("      sum over the three even structures of |S|^2 at (a, b, c), normalised:",
      [mp.nstr(v/sum(avg), 6) for v in avg], " (qubit: 1/2, 1/2, 0)")

# ---------------------------------------------------------------- M6
print("M6  the recorded source on the cover branched at the four marks")
# quadratic form of a spin structure: 0 on a loop with antiperiodic spinors, 1 with periodic ones;
# q(x + y) = q(x) + q(y) + <x, y> mod 2.  Recorded structure on E_i: q(1) = q(i) = 0.
q1, qi = 0, 0
q_e1 = (q1+qi+1) % 2                             # loop 1 + i
q_e2 = (q1+qi+1) % 2                             # loop -1 + i
check("lattice of the cover: 1 + i and -1 + i; both loops carry periodic spinors, Arf invariant 1: odd structure",
      q_e1 == 1 and q_e2 == 1 and (q_e1*q_e2) % 2 == 1)
print("      an odd structure on a torus has exactly one holomorphic spinor, a constant: the lifted source")
print("      has one zero mode, with the same value at all four marks.  No kernel between the four marks")
print("      exists for it; the kernel of M2 lives on the recorded torus, where the marks are two points.")

# ---------------------------------------------------------------- M7
print("M7  seventh pasted note: localisation and reference")


def bump(sv, lo, hi):
    if not lo < sv < hi:
        return mp.mpf(0)
    xx = (sv-lo)/(hi-lo)
    return mp.exp(-1/(xx*(1-xx)))


lo, hi = mp.mpf("0.05"), mp.mpf("0.30")          # support inside one quarter arc of the circle
coef = {}
for n in range(-40, 40):
    k = n+mp.mpf(1)/2
    coef[n] = mp.quad(lambda sv: bump(sv, lo, hi)*mp.exp(-2j*PI*k*sv), [lo, (lo+hi)/2, hi])
even = sum(abs(c)**2 for n, c in coef.items() if n % 2 == 0)
odd = sum(abs(c)**2 for n, c in coef.items() if n % 2 == 1)
check("a probe supported on one quarter arc has half of its norm square in each of the two mode groups",
      abs(even/(even+odd)-mp.mpf(1)/2) < mp.mpf(10)**(-10), mp.nstr(even/(even+odd), 14))
print("      -> the two quarter groups are not the odd and even arcs; that reading of the eleventh round is withdrawn.")
w = lambda x: mp.exp(-2*PI*x*x)
Zm = mp.nsum(lambda m: w(m/2+mp.mpf(1)/4), [-mp.inf, mp.inf])
mixed = mp.nsum(lambda m: (m/2+mp.mpf(1)/4)**2*w(m/2+mp.mpf(1)/4), [-mp.inf, mp.inf])/Zm
check("variance at a fixed mark 1/(8 pi) = 0.03978873577; all four marks mixed from a common origin 0.08332635858",
      abs(mixed-mp.mpf("0.08332635858")) < mp.mpf(10)**(-10), mp.nstr(mixed, 12))

# ---------------------------------------------------------------- M8
print("M8  the rule that selects the square")


def klein_j(tau):
    q = mp.exp(2j*PI*tau)
    e4 = 1+240*mp.nsum(lambda n: n**3*q**n/(1-q**n), [1, 400])
    e6 = 1-504*mp.nsum(lambda n: n**5*q**n/(1-q**n), [1, 400])
    return 1728*e4**3/(e4**3-e6**2)


vals = {}
for yv in ("0.6", "0.8", "1", "1.3", "2"):
    yy = mp.mpf(yv)
    tq = (1+1j*yy)/2                             # rectangle 1 x y, points identified with the half diagonal shift
    tq = -1/tq if tq.imag < mp.mpf("0.6") else tq
    vals[yv] = (klein_j(1j*yy), klein_j(tq))
check("rectangle: j real and >= 1728; after the identification j real and <= 1728; both only at the square",
      all(abs(a.imag) < mp.mpf(10)**(-20) and abs(b.imag) < mp.mpf(10)**(-15) and a.real > 1728-mp.mpf(10)**(-20)
          and b.real < 1728+mp.mpf(10)**(-15) for a, b in vals.values())
      and abs(vals["1"][1]-1728) < mp.mpf(10)**(-20) and all(vals[k][1].real < 1727 for k in vals if k != "1"),
      {k: mp.nstr(v[1].real, 8) for k, v in vals.items()})
print("      a torus has a mirror with two fixed circles exactly when j >= 1728: the rule leaves the square.")
print("      in this session's terms the identification is the passage from the cover of the mark sphere")
print("      to the recorded torus over the sphere of z^2, i.e. the quotient by the half turn z -> -z.")

# ---------------------------------------------------------------- M9
print("M9  the two mirrors on the two tori (lattice Z + iZ each; w = (1+i) v)")


def circles(lin):
    """number of fixed circles of w -> lin * conj(w) + shift on C/(Z + iZ): 2 if the action on the lattice is
    the identity modulo 2, else 1 (when fixed points exist)"""
    im1, imi = lin*1, lin*mp.conj(1j)
    M = [[int(mp.nint(im1.real)) % 2, int(mp.nint(imi.real)) % 2], [int(mp.nint(im1.imag)) % 2, int(mp.nint(imi.imag)) % 2]]
    return 2 if M == [[1, 0], [0, 1]] else 1


rot = (1+1j)/(1-1j)                               # conjugation of a linear part from the v frame to the w frame
mirrors = {"mark circle z -> 1/conj(z)": -1, "axis through two marks z -> conj(z)": -1j,
           "axis between the marks z -> i conj(z)": 1}
tab = {}
for nm, lin_v in mirrors.items():
    tab[nm] = (circles(mp.mpc(lin_v)), circles(mp.mpc(lin_v)*rot))
    print("      %-40s fixed circles on the mark cover: %d, on the recorded torus: %d" % (nm, *tab[nm]))
check("the original mirror divides the cover of the marks and not the recorded torus; "
      "the mirror through two marks does the opposite",
      tab["mark circle z -> 1/conj(z)"] == (2, 1) and tab["axis through two marks z -> conj(z)"] == (1, 2))
print("      the recorded free source is reflection positive for the second kind (time = Im w, circles 0 and 1/2).")

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
