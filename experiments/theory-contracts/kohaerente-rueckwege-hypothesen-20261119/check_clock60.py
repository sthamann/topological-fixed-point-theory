"""CLOCK60DOUBLE20261119 (hypotheses A and E): the order-60 clock on E8 + E8.
Claims under test: U = [[0, C], [1, 0]] with C a Coxeter element has chi_U = Phi_60,
its powers give the tacts 30, 5, 4, 3. Counter-checks named by the hypothesis:
the Galois action of order four on the five-tact, and the reflected (opposite) copy.
Additional control: compatibility with a fixed D5 + A3 frame."""
from common import *
import sympy as sp
from math import gcd

c = Checks("clock60")
x = sp.symbols("x")
edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (4, 7)]
A = 2 * sp.eye(8)
for i, j in edges: A[i, j] = A[j, i] = -1
c.ok("Cartan matrix is E8 (determinant 1)", A.det() == 1)
refl = []
for i in range(8):
    S = sp.eye(8); S[i, :] = S[i, :] - A[:, i].T; refl.append(S)
C = sp.eye(8)
for S in refl: C = C * S
I8, Z8 = sp.eye(8), sp.zeros(8)
c.ok("C is a lattice isometry of order 30 with characteristic polynomial Phi_30",
     C.T * A * C == A and C**30 == I8 and C**15 == -I8
     and sp.expand(C.charpoly(x).as_expr() - sp.cyclotomic_poly(30, x)) == 0)

blk = lambda a, b, cc, d: sp.Matrix(sp.BlockMatrix([[a, b], [cc, d]]))
U = blk(Z8, C, I8, Z8)
G2 = blk(A, Z8, Z8, A)
I16 = sp.eye(16)
c.ok("U preserves the Gram form of E8 + E8 (even, unimodular, positive, rank 16)",
     U.T * G2 * U == G2 and G2.det() == 1)
c.ok("U^2 = C + C, U^30 = -1, U^60 = 1, and no proper divisor of 60 is the order",
     U**2 == blk(C, Z8, Z8, C) and U**30 == -I16 and U**60 == I16
     and all(U**(60 // p) != I16 for p in (2, 3, 5)))
c.ok("chi_U(x) = Phi_30(x^2) = Phi_60(x): primitive order-60 clock of rank 16 = phi(60)",
     sp.expand(U.charpoly(x).as_expr() - sp.cyclotomic_poly(60, x)) == 0
     and sp.expand(sp.cyclotomic_poly(30, x**2) - sp.cyclotomic_poly(60, x)) == 0)


def order(M, cap=120):
    P = M
    for k in range(1, cap + 1):
        if P == sp.eye(M.rows): return k
        P = P * M
tacts = {k: order(U**k) for k in (2, 12, 15, 20)}
c.ok("tacts as powers: U^2 order 30, U^12 order 5, U^15 order 4, U^20 order 3",
     tacts == {2: 30, 12: 5, 15: 4, 20: 3}, tacts)

# --- counter-check 1: Galois action. Find g in W(E8) with g C g^-1 = C^k for every unit k.
roots, frontier = set(), [tuple(int(i == j) for j in range(8)) for i in range(8)]
while frontier:
    r = frontier.pop()
    if r in roots: continue
    roots.add(r)
    v = sp.Matrix(r)
    frontier += [tuple(S * v) for S in refl]
roots = sorted(roots)
c.ok("240 roots; C acts freely with eight orbits of length 30", len(roots) == 240 and
     len({frozenset(tuple(C**k * sp.Matrix(r)) for k in range(30)) for r in roots}) == 8)
Cp = [sp.eye(8)]
for _ in range(60): Cp.append(Cp[-1] * C)
v0 = sp.Matrix(roots[0])
X = sp.Matrix.hstack(*[Cp[k] * v0 for k in range(8)])
Xinv = X.inv()


def conjugator(k):
    """Isometry g of E8 with g C = C^k g, built from the image of one root."""
    for r in roots:
        Y = sp.Matrix.hstack(*[Cp[(k * m) % 30] * sp.Matrix(r) for m in range(8)])
        g = Y * Xinv
        if all(e.is_integer for e in g) and g.T * A * g == A and g * C == Cp[k % 30] * g:
            return g
units60 = [k for k in range(60) if gcd(k, 60) == 1]
Up = [I16]
for _ in range(60): Up.append(Up[-1] * U)
galois = {}
for k in units60:
    g = conjugator(k % 30)
    if g is None: break
    Gk = blk(g, Z8, Z8, Cp[((k - 1) // 2) % 30] * g)
    galois[k] = Gk.T * G2 * Gk == G2 and Gk * U == Up[k] * Gk
c.ok("every Galois conjugate: an isometry G_k of E8 + E8 with G_k U G_k^-1 = U^k for all 16 units k mod 60",
     len(galois) == 16 and all(galois.values()), sorted(galois))
g7 = conjugator(7)
G7 = blk(g7, Z8, Z8, Cp[3] * g7)
five = Up[12]
c.ok("order-four Galois action on the five-tact: G_7 U^12 G_7^-1 = (U^12)^2, and four steps return it",
     G7 * five == five**2 * G7 and G7**4 * five == five * G7**4 and G7**2 * five != five * G7**2)
c.ok("G_7 is not a power of U (powers of U commute with U): the Galois action is extra data",
     G7 * U != U * G7)

# --- counter-check 2: fixed D5 + A3 frame
sizes = {k: 30 // gcd(k, 30) for k in range(1, 30)}
keeps = [k for k, n in sizes.items() if 40 % n == 0 and 12 % n == 0]
c.ok("a D5 + A3 root frame (40 + 12 roots) can only be preserved by C^15 = -1: orbit lengths of C^k must divide 40 and 12",
     keeps == [15], keeps)
c.ok("hence the five-tact U^12 = C^6 + C^6 and the three-tact U^20 = C^10 + C^10 move every D5 + A3 frame",
     sizes[6] == 5 and sizes[10] == 3)

# --- counter-check 3 (hypothesis E): reflected copy with opposite sign, II_{8,8} = E8 + E8(-1)
Gm = blk(A, Z8, Z8, -A)
c.ok("on E8 + E8(-1) the same U is an anti-isometry (U^T G U = -G), its square an isometry",
     U.T * Gm * U == -Gm and (U**2).T * Gm * U**2 == Gm)
odd_inv = [k for k in range(1, 60, 2) if Up[k]**2 == I16]
c.ok("no odd power of U is an involution; U^15 squares to -1", odd_inv == [] and Up[15]**2 == -I16)
Sw = blk(Z8, I8, I8, Z8)
c.ok("the plain swap is the involutive anti-isometry; U = diag(C, 1) * swap is the swap followed by one-sided Coxeter transport",
     Sw**2 == I16 and Sw.T * Gm * Sw == -Gm and U == blk(C, Z8, Z8, I8) * Sw)
c.ok("signature of the reflected pairing <x, swap y> on E8 + E8 is (8, 8)",
     sorted(sp.sign(e) for e in (G2 * Sw).eigenvals(multiple=True)).count(1) == 8)

# --- how special is 60: the same doubling for the other rank-eight clocks
dbl = {n: sp.factor_list(sp.cyclotomic_poly(n, x**2))[1] for n in (15, 16, 20, 24, 30)}
prim = {n: len(f) == 1 for n, f in dbl.items()}
c.ok("doubling is generic: primitive for 16, 20, 24, 30 (orders 32, 40, 48, 60), not for 15", prim ==
     {15: False, 16: True, 20: True, 24: True, 30: True}, prim)
c.done({"tacts": tacts, "galois_units": units60, "frame_preserving_powers_of_C": keeps})
