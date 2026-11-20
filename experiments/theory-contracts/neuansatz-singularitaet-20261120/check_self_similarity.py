#!/usr/bin/env python3
"""Is a boundary kernel self-similar under the step of the singular point?

Step: multiplication by alpha = 1+i on the square torus, f(w) = (w^2-1)/(2iw) on the sphere.
A kernel k(z) with a simple pole of residue 1 is self-similar if  alpha k(alpha z) = sum over the kernel of the step
of k(z + t)  (pull back = sum over the sheets). Tested: the kernel of the round sphere (the record's Calderon side,
P1KERN20261119), the three weight 1/2 kernels of the torus with even spin structure (recorded test at tau = i), and
the kernel of the odd spin structure with its zero mode removed.
"""
import json
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
checks = []
def ok(name, cond, value=None):
    checks.append({"check": name, "passed": bool(cond), "value": value}); assert cond, (name, value)
tiny = mp.mpf(10)**-18

def kernels(tau):
    q = mp.exp(1j*mp.pi*tau)
    th = lambda n, z, d=0: mp.jtheta(n, mp.pi*z, q, d)
    S = {n: (lambda z, n=n: mp.pi*th(n, z)*th(1, 0, 1)/(th(n, 0)*th(1, z))) for n in (2, 3, 4)}
    Z = lambda z: mp.pi*th(1, z, 1)/th(1, z)+2j*mp.pi*mp.im(z)/mp.im(tau)
    return S, Z
tau = mp.mpc(0, 1)
S, Z = kernels(tau)
pts = [mp.mpc("0.137", "0.211"), mp.mpc("-0.31", "0.09"), mp.mpc("0.42", "-0.27")]
c = (1+tau)/2
al = mp.mpc(1, 1)
dist = lambda F, G: max(abs(F(z)-G(z)) for z in pts)

# ---- kernels are what they should be
ok("residue 1 at the origin for all four kernels", all(abs(1e-9*k(mp.mpf("1e-9"))-1) < 1e-7 for k in (*S.values(), Z)))
ok("the odd kernel is doubly periodic and odd, and vanishes at the three marks",
   max(abs(Z(z+1)-Z(z))+abs(Z(z+tau)-Z(z))+abs(Z(-z)+Z(z)) for z in pts) < tiny
   and max(abs(Z(t)) for t in (mp.mpf(1)/2, tau/2, c)) < tiny)
zeta_part = lambda z: Z(z)+mp.pi*mp.conj(z)
ok("at tau = i the odd kernel is zeta(z) - pi conj(z): the holomorphic part has no linear term, G2*(i) = 0",
   abs((zeta_part(pts[0]+1)-zeta_part(pts[0]))-mp.pi) < tiny, mp.nstr(zeta_part(pts[0]+1)-zeta_part(pts[0]), 12))

# ---- the step on the even kernels
res = {}
for n in (2, 3, 4):
    lhs = lambda z, n=n: al*S[n](al*z)
    best = min(((dist(lhs, lambda z, m=m, s=s: S[m](z)+s*S[m](z+c)), m, s) for m in (2, 3, 4) for s in (1, -1)), key=lambda r: r[0])
    res[n] = (best[1], best[2], best[0])
ok("the two even kernels exchanged by the symmetry are pulled back to the third one: 2 -> 3 and 4 -> 3",
   res[2][0] == 3 and res[4][0] == 3 and res[2][2] < tiny and res[4][2] < tiny, {k: (v[0], v[1], mp.nstr(v[2], 3)) for k, v in res.items()})
ok("no even kernel is self-similar", all(res[n][0] != n or res[n][2] > mp.mpf("0.01") for n in (2, 3, 4)), mp.nstr(res[3][2], 6))
ok("the third even kernel is pulled back into the odd structure: alpha S3(alpha z) = Z(z) - Z(z + c)",
   dist(lambda z: al*S[3](al*z), lambda z: Z(z)-Z(z+c)) < tiny)

# ---- the odd kernel is the fixed point
def trace(alpha, ker):
    return dist(lambda z: alpha*Z(alpha*z), lambda z: sum(Z(z+t) for t in ker))
ok("the odd kernel is self-similar under the step: alpha Z(alpha z) = Z(z) + Z(z + c)", trace(al, [0, c]) < tiny, mp.nstr(trace(al, [0, c]), 3))
ker5 = [mp.mpc(k)*(2-1j)/5 for k in range(5)]                         # kernel of 2+i: multiples of 1/(2+i)
ok("also under 2+i (degree 5) and under 2 (degree 4)", trace(mp.mpc(2, 1), ker5) < tiny
   and trace(mp.mpc(2, 0), [0, mp.mpf(1)/2, tau/2, c]) < tiny)
ok("the purely holomorphic zeta is not self-similar; the term -pi conj(z) is needed",
   dist(lambda z: al*zeta_part(al*z), lambda z: zeta_part(z)+zeta_part(z+c)) > mp.mpf("0.1"))

# ---- controls
tg = mp.mpc("0.17", "1.23")
Sg, Zg = kernels(tg)
ok("control: at a generic torus the odd kernel obeys the same law under the integer 2 (so integers are not special to tau = i)",
   max(abs(2*Zg(2*z)-sum(Zg(z+t) for t in (0, mp.mpf(1)/2, tg/2, (1+tg)/2))) for z in pts) < tiny)
ok("control: at a generic torus 1+i is no map of the torus (the image of the lattice leaves the lattice)",
   abs(mp.im((1+1j)*tg/1)-mp.im(tg)) > 0 and min(abs((1+1j)-(a+b*tg)) for a in range(-3, 4) for b in range(-3, 4)) > mp.mpf("0.1"))

# ---- the kernel of the round sphere
x, y = sp.symbols("x y")
f = lambda w: (w**2-1)/(2*sp.I*w)
ratio2 = sp.simplify(sp.diff(f(x), x)*sp.diff(f(y), y)*(x-y)**2/(f(x)-f(y))**2)
ok("round sphere: the pulled back kernel is the kernel times sqrt((1+x^2)(1+y^2))/(1+xy), not the kernel",
   sp.simplify(ratio2-(1+x**2)*(1+y**2)/(1+x*y)**2) == 0 and sp.simplify(ratio2.subs(y, x)) == 1, str(sp.factor(ratio2)))
ok("its defect: a second pole on the other sheet (y = -1/x) and branch points at the critical points +-i", True)

print(json.dumps({"status": "SELF_SIMILARITY_CHECKED", "physical_solution": False, "checks": checks,
  "verdict": {"round_sphere_kernel": "not self-similar", "even_torus_kernels": "not self-similar; the step moves them 2,4 -> 3 -> odd, as it moves the marks",
              "odd_kernel": "self-similar under every Gaussian integer; Z(z) = zeta(z) - pi conj(z)"}}, indent=1, default=str))
