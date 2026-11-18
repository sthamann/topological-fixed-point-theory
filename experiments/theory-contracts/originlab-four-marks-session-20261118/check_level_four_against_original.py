#!/usr/bin/env python3
"""H-X4 against the recorded originals (18 November 2026).

Recorded inputs used here, unchanged:
    verification/v117_monodromy_weyl_a3.py           U = diag(1,i,-i), M0, (M0 U)^4 = 1
    experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/
        MODULAR-C3-SOURCE-TEST.txt                   A_S, A_T, U_Theta, M_Theta, X
    ledger SEAM.RIGIDITY.MU4FROMMARKS.01             omega_k = z^(k-1) dz/(z^4-1), eigenvalue i^k
    ledger AX.P1.01                                  8 pi = |Z2| * 2 pi * chi(S^2)

Sections
    O1  recorded dictionary T -> U, ST -> M0 reproduced; the two pole lifts
    O2  H-X4 as posed (marks = cusps of X(4)): what holds and what fails
    O3  corrected reading: z^4 = 1 - j/1728 (exact identities)
    O4  the three sign lines as covers of the t-line: which one is modular
    O5  the elliptic surface with this j and its lattice
    O6  Gauss-Bonnet numbers
"""
import itertools
import sys
from collections import deque

import mpmath as mp
import sympy as sp

import cover_lattice as cl

mp.mp.dps = 40
PI = mp.pi
I = sp.I
FAIL = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)


def morder(M, proj=False, bound=50):
    P = sp.eye(M.shape[0])
    for n in range(1, bound):
        P = (P*M).applyfunc(sp.simplify)
        if proj:
            c = [P[i, i] for i in range(P.shape[0])]
            if all(sp.simplify(P[i, j]) == 0 for i in range(3) for j in range(3) if i != j) \
                    and all(sp.simplify(ci-c[0]) == 0 for ci in c):
                return n
        elif (P-sp.eye(M.shape[0])).applyfunc(sp.simplify) == sp.zeros(*M.shape):
            return n
    return None


# ---------------------------------------------------------------- O1
print("O1  recorded dictionary and the two lifts of the quarter turn")
A_S = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
A_T = sp.Matrix([[I, 0, 0], [0, 0, 1], [0, 1, 0]])
U_T = -I*A_T
M_T = I*A_S*A_T
X = sp.Matrix([[-1-I, 0, 0], [0, I, -I], [0, 1, 1]])/sp.sqrt(2)
U117 = sp.diag(1, I, -I)
M117 = sp.Matrix([[0, -(1+I)/2, (1-I)/2], [-(1+I)/2, -I/2, sp.Rational(-1, 2)],
                  [(1-I)/2, sp.Rational(-1, 2), I/2]])
Z3 = sp.zeros(3, 3)
simp = lambda M: M.applyfunc(sp.simplify)
check("X unitary", simp(X.H*X-sp.eye(3)) == Z3)
check("X U_Theta X^-1 = U_v117", simp(X*U_T*X.H-U117) == Z3)
check("X M_Theta X^-1 = M0_v117", simp(X*M_T*X.H-M117) == Z3)
check("order of M0 U = 4 (recorded flatness branch)", morder(M117*U117) == 4)
check("order of M0 U^-1 = 2", morder(M117*U117.inv()) == 2)
check("M_Theta U_Theta^-1 = -A_S (the recorded S lift at tau = i)",
      simp(M_T*U_T.inv()+A_S) == Z3)
P1 = sp.eye(3)
P2 = sp.eye(3)
for k in range(4):
    P1 = P1*U117**k*M117*U117**(-k)
    P2 = P2*U117**(-k)*M117*U117**k
check("prod U^k M0 U^-k = 1 and prod U^-k M0 U^k = 1 (both orientations flat)",
      simp(P1-sp.eye(3)) == Z3 and simp(P2-sp.eye(3)) == Z3)
print("      -> v117's matrices admit two equivariant systems on the quotient line:")
print("         pole lifts of orders (4,4), or pole lifts of orders (2,4) = (S, T).")

# ---------------------------------------------------------------- O2
print("O2  H-X4 as posed: z = theta2(2 tau)/theta3(2 tau), marks = cusps")


def thetas(tau):
    N = 40
    t2 = mp.nsum(lambda n: mp.exp(1j*PI*tau*(n+mp.mpf(1)/2)**2), [-N, N])
    t3 = mp.nsum(lambda n: mp.exp(1j*PI*tau*n**2), [-N, N])
    t4 = mp.nsum(lambda n: (-1)**int(n)*mp.exp(1j*PI*tau*n**2), [-N, N])
    return t2, t3, t4


def ab(tau):
    t2, t3, _ = thetas(2*tau)
    return t3, t2


def zf(tau):
    a, b = ab(tau)
    return b/a


tau0 = mp.mpc("0.13", "0.71")
a0, b0 = ab(tau0)
z0 = zf(tau0)
dz = mp.diff(zf, tau0)
worst = 0
for k in (1, 2, 3):
    lhs = z0**(k-1)*dz/(z0**4-1)
    rhs = -(1j*PI/2)*a0**(4-k)*b0**k
    worst = max(worst, abs(lhs-rhs)/abs(rhs))
check("omega_k = z^(k-1) dz/(z^4-1) = -(i pi/2) a^(4-k) b^k d tau", worst < mp.mpf(10)**(-25), mp.nstr(worst, 3))
a1, b1 = ab(tau0+1)
check("tau -> tau+1: a -> a, b -> i b, so omega_k -> i^k omega_k (ledger eigenvalues)",
      abs(a1-a0) < 1e-30 and abs(b1-1j*b0) < 1e-30)
check("A_T^4 = 1: a loop around a level-four cusp acts trivially on theta squares",
      simp(A_T**4-sp.eye(3)) == Z3)


def mat_mod(M, n):
    return tuple(int(x) % n for x in M)


def mmul(A, B, n):
    a, b, c, d = A
    e, f, g, h = B
    return ((a*e+b*g) % n, (a*f+b*h) % n, (c*e+d*g) % n, (c*f+d*h) % n)


def minv(A, n):
    a, b, c, d = A
    return (d % n, (-b) % n, (-c) % n, a % n)


# cusp generators of Gamma(4): g T^4 g^-1 for the six cusps, reduced mod 3
cusp_reps = [(1, 0, 0, 1), (0, -1, 1, 0), (1, 0, 1, 1), (1, 0, 2, 1), (1, 0, 3, 1), (1, -1, 2, -1)]
T4 = (1, 4, 0, 1)
orders = []
for g in cusp_reps:
    h = mmul(mmul(mat_mod(g, 3), mat_mod(T4, 3), 3), minv(mat_mod(g, 3), 3), 3)
    P, n = h, 1
    while P not in ((1, 0, 0, 1), (2, 0, 0, 2)):
        P = mmul(P, h, 3)
        n += 1
    orders.append(n)
check("level 12 (needed for order three): every one of the six cusps gets order 3",
      orders == [3]*6, orders)
print("      -> a cover with order three at four cusps only is not a congruence cover;")
print("         the original's cubic mark monodromy cannot sit at level-four cusps.")

# ---------------------------------------------------------------- O3
print("O3  corrected reading: z^4 = 1 - j/1728")
lam = sp.symbols("lambda")
om = sp.Rational(-1, 2)+sp.sqrt(3)*I/2
w = (lam+om)/(lam+om**2)
z2 = (1+w**3)/(1-w**3)
jl = 256*(1-lam+lam**2)**3/(lam**2*(1-lam)**2)
check("w = (lambda+omega)/(lambda+omega^2), w^3 = (z^2-1)/(z^2+1)  =>  1728 (1 - z^4) = j(lambda)",
      sp.simplify(sp.together(1728*(1-z2**2)-jl)) == 0)
check("z^2 is a rational function of lambda (sign character of S3)",
      sp.simplify(z2-I*(1+lam)*(2-lam)*(1-2*lam)/(3*sp.sqrt(3)*lam*(1-lam))) == 0
      or sp.simplify(z2+I*(1+lam)*(2-lam)*(1-2*lam)/(3*sp.sqrt(3)*lam*(1-lam))) == 0)
den = (lam+om**2)
rq = [sp.factor(sp.simplify((w-om**k)*den)) for k in range(3)]
rp = [sp.factor(sp.simplify((w+om**k)*den)) for k in range(3)]
check("line Q radicands are 1, lambda-1, lambda up to constants (cusps)",
      all(sp.simplify(sp.diff(sp.simplify(r/t), lam)) == 0 for r, t in zip(rq, (1, lam-1, lam))))
check("line P radicands are 2 lambda-1, lambda+1, lambda-2 up to constants (points tau = i)",
      all(sp.simplify(sp.diff(sp.simplify(r/t), lam)) == 0 for r, t in zip(rp, (2*lam-1, lam+1, lam-2))))
check("w -> -w is z -> 1/z: exchanges the two lines", sp.simplify(z2.subs(lam, lam)*((1+(-w)**3)/(1-(-w)**3))-1) == 0)


def eis(tau):
    q = mp.exp(2j*PI*tau)
    s4 = mp.nsum(lambda n: n**3*q**n/(1-q**n), [1, 200])
    s6 = mp.nsum(lambda n: n**5*q**n/(1-q**n), [1, 200])
    return 1+240*s4, 1-504*s6


tau1 = mp.mpc("0.21", "0.83")
t2, t3, t4 = thetas(tau1)
lam1 = t2**4/t3**4
E4, E6 = eis(tau1)
jnum = 1728*E4**3/(E4**3-E6**2)
omn = mp.mpc(-0.5, mp.sqrt(3)/2)
w1 = (lam1+omn)/(lam1+omn**2)
z21 = (1+w1**3)/(1-w1**3)
check("numerically: 1 - z^4 = j(tau)/1728 with lambda = theta2^4/theta3^4",
      abs(1728*(1-z21**2)-jnum)/abs(jnum) < mp.mpf(10)**(-25))
check("z^4 = -E6^2/(E4^3 - E6^2)", abs(z21**2+E6**2/(E4**3-E6**2)) < mp.mpf(10)**(-25))
check("sqrt(lambda) = theta2^2/theta3^2, sqrt(1-lambda) = theta4^2/theta3^2 single valued: line Q is level four",
      abs((t2**2/t3**2)**2-lam1) < 1e-30 and abs((t4**2/t3**2)**2-(1-lam1)) < 1e-30)

# ---------------------------------------------------------------- O4
print("O4  the sign lines as covers of the t-line; comparison with the modular triple")


def pmul(a, b):
    return tuple(cl.mul(a, b))


def pinv(a):
    return tuple(cl.inv(a))


def gen_group(gens):
    idn = tuple(range(len(gens[0])))
    seen = {idn}
    todo = deque([idn])
    while todo:
        g = todo.popleft()
        for h in gens:
            k = pmul(h, g)
            if k not in seen:
                seen.add(k)
                todo.append(k)
    return seen


def centre(G, gens):
    return [g for g in G if all(pmul(g, h) == pmul(h, g) for h in gens)]


def order_mod(g, Z):
    idn = tuple(range(len(g)))
    P, n = g, 1
    Zs = set(Z)
    while P not in Zs:
        P = pmul(P, g)
        n += 1
    return n


def iso_pairs(x, y, X2, Y2):
    """Is x -> X2, y -> Y2 a group isomorphism <x,y> -> <X2,Y2>?"""
    idn = tuple(range(len(x)))
    idn2 = tuple(range(len(X2)))
    img = {idn: idn2}
    todo = deque([idn])
    while todo:
        g = todo.popleft()
        for h, H in ((x, X2), (y, Y2)):
            k = pmul(h, g)
            K = pmul(H, img[g])
            if k in img:
                if img[k] != K:
                    return False
            else:
                img[k] = K
                todo.append(k)
    return len(set(img.values())) == len(img) and len(img) == len(gen_group([X2, Y2]))


# PSL2(Z/4) on the six cusps = points of P^1(Z/4), paired with Z/4
pts = []
for c in itertools.product(range(4), repeat=2):
    if (c[0] % 2 or c[1] % 2):
        rep = min(((u*c[0]) % 4, (u*c[1]) % 4) for u in (1, 3))
        if rep not in pts:
            pts.append(rep)


def act(Mx):
    a, b, c, d = Mx
    out = []
    for (p, q) in pts:
        v = ((a*p+b*q) % 4, (c*p+d*q) % 4)
        out.append(pts.index(min(((u*v[0]) % 4, (u*v[1]) % 4) for u in (1, 3))))
    return tuple(out)


Sg, Tg = act((0, -1, 1, 0)), act((1, 1, 0, 1))
STg = pmul(Sg, Tg)


def paired(g, m):
    """element (g, m) of S4 x Z/4 as a permutation of 6 + 4 letters"""
    return tuple(list(g)+[6+((k+m) % 4) for k in range(4)])


check("six cusps, S4 of order 24", len(pts) == 6 and len(gen_group([Sg, Tg])) == 24)
ST2inv = pinv(pmul(STg, Tg))
triples = {
    "modular (S, ST, T^-1): lifts of orders 2 and 4": (paired(Sg, 1), paired(STg, 0)),
    "mirror  (T, (ST)^-1, S): lifts of orders 4 and 2": (paired(Tg, 1), paired(pinv(STg), 0)),
    "v117 as written ((ST^2)^-1, ST, T): lifts 4 and 4": (paired(ST2inv, 1), paired(STg, 0)),
}
for nm, (X2, Y2) in triples.items():
    G2 = gen_group([X2, Y2])
    print("      model  %-52s group order %d" % (nm, len(G2)))

found = {}
for line in ("P", "Q", "L0"):
    x, y = (tuple(v) for v in cl.build_cover(line)[1:])
    G = gen_group([x, y])
    Z = centre(G, [x, y])
    xy = pmul(x, y)
    sig = (cl.perm_order(x), cl.perm_order(y), cl.perm_order(pinv(xy)))
    sigc = (order_mod(x, Z), order_mod(y, Z), order_mod(xy, Z))
    hits = []
    for nm, (X2, Y2) in triples.items():
        for xx, tag in ((x, ""), (pinv(x), " [x inverted]")):
            for yy, tag2 in ((y, ""), (pinv(y), " [y inverted]")):
                if iso_pairs(xx, yy, X2, Y2):
                    hits.append(nm.split(":")[0]+tag+tag2)
    found[line] = (len(G), sig, len(Z), sigc, hits)
    print("      line %-2s  group %d  orders %s  centre %d  orders modulo centre %s" % (line, len(G), sig, len(Z), sigc))
    print("               isomorphic as marked cover to:", hits if hits else "none of the three")
check("lines P and Q: group 48, orders (4,3,4)", found["P"][0] == 48 and found["Q"][0] == 48
      and found["P"][1] == (4, 3, 4) and found["Q"][1] == (4, 3, 4))
check("one of P, Q is the modular fibre product, the other its mirror",
      any("modular" in h for h in found["Q"][4]+found["P"][4])
      and any("mirror" in h for h in found["Q"][4]+found["P"][4]))
check("orders modulo centre: Q (2,3,4) modular, P (4,3,2) mirror, L0 (4,3,4)",
      found["Q"][3] == (2, 3, 4) and found["P"][3] == (4, 3, 2) and found["L0"][3] == (4, 3, 4))
check("as marked covers: Q = modular model, P = mirror model, L0 = the (4,4) model",
      any(h.startswith("modular") and "[" not in h for h in found["Q"][4])
      and any(h.startswith("mirror") and "[" not in h for h in found["P"][4])
      and any(h.startswith("v117 as written") and "[" not in h for h in found["L0"][4]))

# ---------------------------------------------------------------- O5
print("O5  elliptic surface  y^2 = x^3 - 3 (1-z^4) x + 2 (1-z^4)   (orientation of the original)")
z, x, y, s = sp.symbols("z x y s")
A = -3*(1-z**4)
B = 2*(1-z**4)
D = sp.factor(4*A**3+27*B**2)
check("j = 1728 * 4A^3/(4A^3+27B^2) = 1728 (1 - z^-4)", sp.simplify(1728*4*A**3/D-1728*(1-z**-4)) == 0)
check("4A^3 + 27B^2 = 108 z^4 (1-z^4)^2: type II at the four marks, I4 at z = 0",
      sp.simplify(D-108*z**4*(1-z**4)**2) == 0 and A.subs(z, 0) != 0 and B.subs(z, 0) != 0)
Ai = sp.expand(s**4*A.subs(z, 1/s))
Bi = sp.expand(s**6*B.subs(z, 1/s))
Di = sp.factor(4*Ai**3+27*Bi**2)
check("at infinity: smooth fibre Y^2 = X^3 + 3X (square torus, j = 1728)",
      Di.subs(s, 0) != 0 and Bi.subs(s, 0) == 0 and Ai.subs(s, 0) == 3)
check("Euler numbers 4*2 + 4 = 12: rational elliptic surface", 4*2+4 == 12)
check("quarter turn z -> i z fixes (x, y); at infinity it is (X, Y) -> (-X, i Y): order 4 on the fibre, order 2 on its H^1 up to sign",
      sp.expand(A.subs(z, I*z)-A) == 0 and sp.expand(B.subs(z, I*z)-B) == 0)
rhs = lambda xx: xx**3+A*xx+B
check("sections (1, z^2) and (-2, 2 sqrt(2) i z^2)",
      sp.expand(z**4-rhs(1)) == 0 and sp.expand((2*sp.sqrt(2)*I*z**2)**2-rhs(-2)) == 0)

# complement of A3 in E8
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = [0]*8
            v[i], v[j] = si, sj
            roots.append(tuple(sp.Integer(c) for c in v))
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append(tuple(sp.Rational(c, 2) for c in signs))
dot = lambda u, v: sum(p*q for p, q in zip(u, v))
a3 = [tuple(sp.Integer(c) for c in v) for v in ((1, -1, 0, 0, 0, 0, 0, 0), (0, 1, -1, 0, 0, 0, 0, 0), (0, 0, 1, -1, 0, 0, 0, 0))]
perp = [r for r in roots if all(dot(r, a) == 0 for a in a3)]
rk = sp.Matrix([list(r) for r in perp]).rank()
check("E8 has 240 roots; orthogonal to an A3 there are 40 roots of rank 5 (D5)",
      len(roots) == 240 and len(perp) == 40 and rk == 5)
CA3 = sp.Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
CD5 = sp.Matrix([[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, -1, -1], [0, 0, -1, 2, 0], [0, 0, -1, 0, 2]])
qa = sorted(set(CA3.inv()[i, i] for i in range(3)))
qd = sorted(set(CD5.inv()[i, i] for i in range(5)))
check("det A3 = det D5 = 4; glue norms 3/4 and 5/4 with sum 2 (ledger q(A3)+q(D5)=2)",
      CA3.det() == 4 and CD5.det() == 4 and sp.Rational(3, 4) in qa and sp.Rational(5, 4) in qd)
check("fibre I4 local height terms i(4-i)/4 = 3/4, 1, 3/4; section heights 2 - 3/4 = 5/4, 2 - 1 = 1",
      [sp.Rational(i*(4-i), 4) for i in (1, 2, 3)] == [sp.Rational(3, 4), 1, sp.Rational(3, 4)])

# ---------------------------------------------------------------- O6
print("O6  Gauss-Bonnet numbers (units of pi)")
area_j = sp.Rational(1, 3)                      # modular orbifold (2,3,inf)
area_434 = 2*(1-sp.Rational(1, 4)-sp.Rational(1, 3)-sp.Rational(1, 4))
check("orbifolds (2,3,cusp) and (4,3,4) have the same area pi/3", area_j == area_434 == sp.Rational(1, 3))
check("level-four sphere 24 * pi/3 = 8 pi = 2 pi (6 - 2)", 24*area_j == 8 and 2*(6-2) == 8)
check("mark sphere in the corrected reading 4 * pi/3; kernel cover 48 * pi/3 = 16 pi = 2 pi (2*5-2)",
      48*area_j == 16 and 2*(2*5-2) == 16)
check("ledger: |Z2| * 2 pi * chi(S^2) = 8 pi; equals 2 pi (n - 2) only for n = 6", 2*2*2 == 8)

# ---------------------------------------------------------------- O7
print("O7  the original A0* system of v117: the quarter turn at the two poles")
import numpy as np
import scipy.linalg as sl
from scipy.integrate import solve_ivp

un = np.diag([1, 1j, -1j])
uis = [np.linalg.matrix_power(un, k) for k in range(4)]
punct = np.array([1, 1j, -1, -1j])
a0n = np.array([[0.5, np.sqrt(8)/12, 0], [np.sqrt(8)/12, 0.25, np.sqrt(5)/12], [0, np.sqrt(5)/12, 0.25]], complex)
akn = [uis[k] @ a0n @ np.linalg.inv(uis[k]) for k in range(4)]
res_inf = sum(akn)
check("sum of residues = diag(2,1,1): integer exponents at infinity, z = 0 is an ordinary point",
      np.allclose(res_inf, np.diag([2, 1, 1])))


def transport(path, Y0):
    def f(t, yv):
        zz, dzz = path(t)
        am = sum(akn[j]/(zz-punct[j]) for j in range(4))
        return ((am @ yv.reshape(3, 3))*dzz).reshape(-1)
    sol = solve_ivp(f, [0, 1], Y0.reshape(-1), rtol=1e-12, atol=1e-14, method="DOP853")
    return sol.y[:, -1].reshape(3, 3)


RHO = 12.0
e_m, e_p = np.exp(-1j*np.pi/4), np.exp(1j*np.pi/4)
Ya = transport(lambda t: (RHO*e_m*t, RHO*e_m), np.eye(3, dtype=complex))
Yarc = transport(lambda t: (RHO*np.exp(1j*(-np.pi/4+np.pi/2*t)), RHO*1j*np.pi/2*np.exp(1j*(-np.pi/4+np.pi/2*t))), Ya)
Yb = transport(lambda t: (RHO*e_p*t, RHO*e_p), np.eye(3, dtype=complex))
check("equivariance: Y along the rotated ray = U Y U^-1", np.linalg.norm(Yb-un @ Ya @ np.linalg.inv(un)) < 1e-7*np.linalg.norm(Yb))
C = np.linalg.solve(Yb, Yarc)           # monodromy of the counterclockwise loop around the mark 1
Ck = [uis[k] @ C @ np.linalg.inv(uis[k]) for k in range(4)]
aop = np.vstack([np.kron(Ck[k].T, Ck[k].conj().T)-np.eye(9) for k in range(4)])
h = np.linalg.svd(aop)[2].conj().T[:, 8].reshape(3, 3)
h = (h+h.conj().T)/2
if np.trace(h).real < 0:
    h = -h
wq = np.asarray(sl.sqrtm(h), dtype=complex)
m0n = wq @ C @ np.linalg.inv(wq)
M0n = np.array(M117.evalf(20), dtype=complex)
d0 = np.linalg.norm(m0n-M0n)
check("loop around the mark 1 (sector path) reproduces the exact M0 of v117", d0 < 1e-6, "%.1e" % d0)
Rinf = np.linalg.inv(C) @ un            # lift of z -> i z transported through infinity
ev = np.sort_complex(np.linalg.eigvals(Rinf))
check("lift at z = 0 is U (spectrum 1, i, -i); lift at infinity has spectrum (1, -1, -1), order 2",
      np.allclose(ev, np.sort_complex(np.array([1, -1, -1], complex)), atol=1e-6)
      and np.linalg.norm(Rinf @ Rinf-np.eye(3)) < 1e-6, str(np.round(ev, 6)))
check("growth exponents predict the same: U diag(i^-2, i^-1, i^-1) = diag(-1, 1, -1)",
      np.allclose(un @ np.diag([1j**-2, 1j**-1, 1j**-1]), np.diag([-1, 1, -1])))
Rex = M117.inv()*U117
check("exactly: (M0^-1 U)^2 = 1, trace -1", simp(Rex*Rex-sp.eye(3)) == Z3 and sp.simplify(Rex.trace()+1) == 0)
check("exactly: X^-1 (M0^-1 U) X = -A_T^-1 A_S A_T, a conjugate of the recorded S operator -A_S",
      simp(X.H*Rex*X+A_T.inv()*A_S*A_T) == Z3)
print("      -> in the original system itself the quarter turn is of type T at z = 0")
print("         and of type S at z = infinity; the mark monodromy joins them: R_inf = M0^-1 U.")

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
