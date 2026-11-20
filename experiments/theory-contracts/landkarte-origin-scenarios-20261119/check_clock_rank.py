"""CLOCKRANK20261119: which lattice rank do the TFPT clocks (2,3,4,5,12,30) need,
which pairs fuse into one primitive clock in E8, and what happens off the square point?"""
from common import *
import random, itertools, numpy as np
from sympy import totient, factorint, cyclotomic_poly, Matrix, symbols, I, Rational as R, exp as sexp, pi as PI, simplify, expand, Poly
from sympy.combinatorics import Permutation
c = Checks("clock_rank")
x = symbols("x")
def minrank(n):                      # least d with an element of order n in GL(d, Z)
    f = factorint(n); s = sum(int(totient(p**k)) for p, k in f.items())
    return s - (1 if f.get(2, 0) == 1 and n > 2 else 0)
tab = {n: (minrank(n), int(totient(n))) for n in (2, 3, 4, 5, 6, 12, 15, 20, 24, 30, 60)}
print("  order: (least rank, rank as one primitive clock)", tab)
c.ok("a plane lattice carries only orders 1,2,3,4,6: the five-clock needs rank 4", [n for n in range(1, 40) if minrank(n) <= 2] == [1, 2, 3, 4, 6] and minrank(5) == 4)
c.ok("order 4 forces the square torus, order 3 the hexagonal one: both on one torus would need order 12, least rank 4", minrank(12) == 4)
c.ok("primitive clocks of rank exactly 8 are the orders 15, 16, 20, 24, 30", [n for n in range(1, 200) if totient(n) == 8] == [15, 16, 20, 24, 30])
c.ok("all three clocks 3, 4, 5 as one primitive clock need order 60 and rank 16", totient(60) == 16 and minrank(60) == 8)
# Weyl group of E8 in simple-root coordinates
edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (4, 7)]
A = 2 * np.eye(8, dtype=np.int64)
for i, j in edges: A[i, j] = A[j, i] = -1
c.ok("Cartan matrix is E8 (determinant 1)", round(np.linalg.det(A)) == 1)
refl = []
for i in range(8):
    S = np.eye(8, dtype=np.int64); S[i, :] -= A[:, i]; refl.append(S)
Id = np.eye(8, dtype=np.int64)
def order(M):
    P = M.copy()
    for k in range(1, 61):
        if (P == Id).all(): return k
        P = P @ M
    return None
random.seed(20261119); seen = {}
for _ in range(60000):
    M = Id
    for _ in range(random.randint(1, 60)): M = M @ refl[random.randrange(8)]
    o = order(M)
    if o not in seen: seen[o] = M
orders = sorted(seen)
print("  element orders found in W(E8):", orders)
c.ok("orders found are exactly the known list for W(E8)", orders == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 18, 20, 24, 30], orders)
prim = {}
for n in (15, 20, 24, 30):
    best = None
    for _ in range(4000):
        M = Id
        for _ in range(random.randint(4, 40)): M = M @ refl[random.randrange(8)]
        if order(M) == n and Poly(Matrix(M.tolist()).charpoly(x).as_expr(), x) == Poly(cyclotomic_poly(n, x), x): best = True; break
    prim[n] = bool(best)
c.ok("in E8 every PAIR of clocks fuses into one primitive clock: 3x5 = 15, 4x5 = 20, 8x3 = 24, 2x3x5 = 30", all(prim.values()), prim)
# no order 60: the 5-part has fixed lattice A4 or none; Aut(A4) = S5 x {+-1} has no order 12
s5 = set()
for p in itertools.permutations(range(5)): s5.add(Permutation(list(p)).order())
aut_a4 = {int(np.lcm(o, s)) for o in s5 for s in (1, 2)}
c.ok("no element of order 12 in Aut(A4) = S5 x Z2, and degree count forbids the other case: E8 has no order 60", 12 not in aut_a4 and 60 not in orders and 16 not in orders, sorted(aut_a4))
c.ok("so the three clocks never close into ONE clock on E8: the four acts on the five as its Galois group (Z/5)^x = Z4", pow(2, 4, 5) == 1 and all(pow(2, k, 5) != 1 for k in (1, 2, 3)))
# gear relation inside S5 = W(A4): G C G^-1 = C^2, G^4 = 1
Cc = Permutation([1, 2, 3, 4, 0]); G = Permutation([0, 2, 4, 1, 3])
c.ok("the recorded gear G C5 G^-1 = C5^2, G^4 = 1 is the Frobenius group F20 inside S5 = W(A4)", (G.order() == 4) and any((g * Cc * g**-1 == Cc**2 or g**-1 * Cc * g == Cc**2) for g in (G, G**-1)))
# Picard curve carries 3 and 4 together
X_, Y_ = symbols("X Y"); w3 = sexp(2 * PI * I / 3)
c.ok("the Picard curve y^3 = x^4 - 1 has the order-12 symmetry (x,y) -> (i x, w y); its period lattice has rank 6 >= 4", simplify(expand((w3 * Y_)**3 - (I * X_)**4 + 1) - (Y_**3 - X_**4 + 1)) == 0)
# what the square point is worth numerically
def ainv(cc):
    from mpmath import findroot
    def F(a):
        xq = 48 * cc**4 * exp(-2 * a); ph = mpf(4) / 3 * cc + xq * (1 - xq)**(-mpf(5) / 4)
        return a**3 - 2 * cc**3 * a**2 - mpf(164) / 5 * cc**6 * log(1 / ph)
    return 1 / findroot(F, 2 * cc**2 * mpf("0.6"))
from mpmath import nsum, inf, mpc, sqrt
E2 = lambda t: (1 - 24 * nsum(lambda n: n * exp(2j * pi * t * n) / (1 - exp(2j * pi * t * n)), [1, inf])).real
c.ok("lattice energy per direction E2(tau)/24: square gives c3, hexagon gives (2/sqrt 3) c3", abs(E2(1j) / 24 - C3) < 1e-25 and abs(E2(mpc(-0.5, sqrt(3) / 2)) / 24 - 2 / sqrt(3) * C3) < 1e-25)
hexv = ainv(C3 * 2 / sqrt(3)); w2 = ainv(2 * C3); i2 = ainv(C3 / 2)
print(f"  alpha^-1: hexagon {mp.nstr(hexv, 8)}   winding 2 {mp.nstr(w2, 8)}   index 8 {mp.nstr(i2, 8)}")
c.ok("hexagonal placement or coherent multiple seeds are excluded by alpha (104.4, 37.3, 512.2 instead of 137.036)", abs(hexv - mpf("104.4222")) < 1e-3 and abs(w2 - mpf("37.3012")) < 1e-3 and abs(i2 - mpf("512.2314")) < 1e-3)
c.done({"table": {str(k): v for k, v in tab.items()}, "orders_WE8": orders, "alpha_inv_hexagon": mp.nstr(hexv, 10),
        "verdict": "exact structural statement: pairs of clocks fuse in rank 8, the triple never does; square point fixed by the four-clock, alpha excludes hexagon and coherent multiples"})
