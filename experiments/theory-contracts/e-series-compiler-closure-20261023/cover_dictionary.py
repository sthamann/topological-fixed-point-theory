#!/usr/bin/env python3
"""Exact common-cover and marked-clock certificate.

Pure mathematics in the active marked-seam program.  This checker does not
change the verification ledger.  It compares the already fixed objects

    C_n : u_n**n = x**4 - 1,

the original Gaussian E8 operators J, sigma and c = J sigma, and the covered
seam lift of v623.  The comparison is by actual morphisms and generator
relations, not by equality of orders or genera.
"""

from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd, lcm

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def cyclic_order(exponent, modulus):
    return modulus // gcd(exponent, modulus)


def pair_power(pair, power, deck_modulus):
    """Power in C4(base rotation) x C_n(deck)."""
    base, deck = pair
    return (power * base % 4, power * deck % deck_modulus)


def pair_order(pair, deck_modulus):
    return lcm(cyclic_order(pair[0], 4), cyclic_order(pair[1], deck_modulus))


def cover_record(n):
    infinity_points = gcd(n, 4)
    infinity_index = n // infinity_points
    genus = (3 * n - infinity_points - 2) // 2
    peripheral = (1, 1, 1, 1, (-4) % n)
    require(sum(peripheral) % n == 0, f"C_{n}: peripheral product")
    return {
        "degree": n,
        "genus": genus,
        "finite_ramification": n,
        "infinity_points": infinity_points,
        "infinity_ramification": infinity_index,
        "peripheral_exponents": peripheral,
    }


# ---------------------------------------------------------------- covers

COVERS = {n: cover_record(n) for n in (2, 3, 4, 6, 12)}
require([COVERS[n]["genus"] for n in (2, 3, 4, 6, 12)] == [1, 3, 3, 7, 15],
        "compact genera")
require([COVERS[n]["infinity_ramification"] for n in (2, 3, 4, 6, 12)]
        == [1, 3, 1, 3, 3], "infinity ramification")


def check_morphism(n, m):
    """C_n -> C_m, u_m = u_n**(n/m), for m dividing n."""
    require(n % m == 0, f"C_{n} does not dominate C_{m}")
    power = n // m
    require(power * m == n, "equation does not commute")
    # Delta_n maps to Delta_m because zeta_n**(n/m) = zeta_m.
    for exponent in range(n):
        image = exponent % m
        require((exponent * power) % n == (image * power) % n,
                f"deck exponent map C_{n}->C_{m}")
    # Peripheral vectors reduce under the same deck-generator map.
    require(tuple(a % m for a in COVERS[n]["peripheral_exponents"])
            == COVERS[m]["peripheral_exponents"],
            f"peripheral naturality C_{n}->C_{m}")


for source, target in ((6, 2), (6, 3), (12, 2), (12, 3),
                       (12, 4), (12, 6)):
    check_morphism(source, target)


def crt_image_size(a, b):
    modulus = lcm(a, b)
    return len({(k % a, k % b) for k in range(modulus)})


require(lcm(2, 3) == 6 and crt_image_size(2, 3) == 6,
        "C6 is the connected minimal common C2/C3 cover")
require(lcm(4, 3) == 12 and crt_image_size(4, 3) == 12,
        "C12 is the connected minimal common C4/C3 cover")


# ------------------------------------------------ actual E8 compiler clock

def e8_roots_doubled():
    roots = []
    for v in product(range(-1, 2), repeat=8):
        if sum(a * a for a in v) == 2 and sum(v) % 2 == 0:
            roots.append(tuple(2 * a for a in v))
    for y in product((0, -1), repeat=8):
        v = tuple(2 * a + 1 for a in y)
        if sum(a * a for a in v) == 8 and sum(v) % 4 == 0:
            roots.append(v)
    return tuple(sorted(roots))


def J_vec(x):
    out = []
    for k in range(0, 8, 2):
        out.extend((-x[k + 1], x[k]))
    return tuple(out)


SIGMA_IDX = (4, 5, 0, 1, 2, 3, 6, 7)


def sigma_vec(x):
    return tuple(x[k] for k in SIGMA_IDX)


def compose(left, right):
    return lambda x: left(right(x))


IDENTITY = lambda x: tuple(x)


def map_power(f, n):
    out = IDENTITY
    for _ in range(n):
        out = compose(f, out)
    return out


ROOTS = e8_roots_doubled()
ROOT_SET = set(ROOTS)
require(len(ROOTS) == 240, "E8 root shell")
require({J_vec(v) for v in ROOTS} == ROOT_SET, "J preserves E8 roots")
require({sigma_vec(v) for v in ROOTS} == ROOT_SET, "sigma preserves E8 roots")
require(all(J_vec(sigma_vec(v)) == sigma_vec(J_vec(v)) for v in ROOTS),
        "J and sigma commute")

C_VEC = compose(J_vec, sigma_vec)
require(all(map_power(J_vec, 4)(v) == v for v in ROOTS), "J^4=1")
require(all(map_power(sigma_vec, 3)(v) == v for v in ROOTS), "sigma^3=1")
require(all(map_power(C_VEC, 4)(v) == sigma_vec(v) for v in ROOTS),
        "c^4=sigma")
require(all(map_power(C_VEC, 9)(v) == J_vec(v) for v in ROOTS), "c^9=J")
require(all(map_power(C_VEC, 12)(v) == v for v in ROOTS), "c^12=1")
require(all(any(map_power(C_VEC, d)(v) != v for v in ROOTS)
            for d in divisors(12) if d < 12), "c has exact order12")

root_index = {v: i for i, v in enumerate(ROOTS)}
c_perm = [root_index[C_VEC(v)] for v in ROOTS]
seen = set()
cycle_census = Counter()
for start in range(240):
    if start in seen:
        continue
    j = start
    length = 0
    while j not in seen:
        seen.add(j)
        length += 1
        j = c_perm[j]
    cycle_census[length] += 1
require(cycle_census == Counter({12: 19, 4: 3}),
        "actual compiler c root-cycle census")


# -------------------------------------------- actual carrier representations

X = sp.symbols("X")
I = sp.I
PHI1 = X - 1
PHI3 = X**2 + X + 1
PHI4 = X**2 + 1
PHI6 = X**2 - X + 1
PHI12 = X**4 - X**2 + 1


def linear_matrix(f, n):
    """Matrix of an already fixed integer coordinate action."""
    cols = []
    for j in range(n):
        e = [0] * n
        e[j] = 1
        cols.append(sp.Matrix(f(tuple(e))))
    return sp.Matrix.hstack(*cols)


def matrix_polynomial(A, p):
    coeffs = sp.Poly(sp.expand(p), X).all_coeffs()
    out = sp.zeros(A.rows)
    for a in coeffs:
        out = out * A + a * sp.eye(A.rows)
    return sp.simplify(out)


def vanishes(A, p):
    return matrix_polynomial(A, p) == sp.zeros(A.rows)


J8 = linear_matrix(J_vec, 8)
SIGMA8 = linear_matrix(sigma_vec, 8)
C8 = J8 * SIGMA8
require(C8 == linear_matrix(C_VEC, 8), "actual compiler real8 matrix")
require(sp.factor(C8.charpoly(X).as_expr()) == PHI4**2 * PHI12,
        "real8 compiler characteristic polynomial")
require(vanishes(C8, PHI4 * PHI12)
        and not vanishes(C8, PHI4) and not vanishes(C8, PHI12),
        "real8 compiler minimal polynomial")

# In z_j=x_(2j)+i x_(2j+1), the fixed real chart is the realification of
# c=i*sigma on V4, with sigma cycling the first three coordinates.
P4 = sp.zeros(4)
for source, target in enumerate((1, 2, 0, 3)):
    P4[target, source] = 1
C4_COMPLEX = I * P4


def realify(A):
    out = sp.zeros(2 * A.rows)
    for a in range(A.rows):
        for b in range(A.cols):
            re, im = sp.expand_complex(A[a, b]).as_real_imag()
            out[2*a:2*a+2, 2*b:2*b+2] = sp.Matrix([[re, -im], [im, re]])
    return out


require(realify(C4_COMPLEX) == C8, "V4 realification is the actual E8 chart")
require(sp.factor(C4_COMPLEX.charpoly(X).as_expr(), extension=I)
        == (X - I)**2 * (X**2 + I*X - 1), "V4 complex characteristic polynomial")
require(vanishes(C4_COMPLEX, (X-I)*(X**2+I*X-1))
        and not vanishes(C4_COMPLEX, X-I)
        and not vanishes(C4_COMPLEX, X**2+I*X-1), "V4 minimal polynomial")


def wedge2_matrix(A):
    pairs = list(combinations(range(A.rows), 2))
    out = sp.zeros(len(pairs))
    for col, (i, j) in enumerate(pairs):
        for row, (a, b) in enumerate(pairs):
            out[row, col] = A[a, i]*A[b, j] - A[a, j]*A[b, i]
    return out


W6 = wedge2_matrix(C4_COMPLEX)
require(sp.factor(W6.charpoly(X).as_expr()) == (X + 1)**2 * PHI6**2,
        "actual exterior2 V4 characteristic polynomial")
require(vanishes(W6, (X+1)*PHI6)
        and not vanishes(W6, X+1) and not vanishes(W6, PHI6),
        "actual exterior2 V4 minimal polynomial")

# The five supports are the literal sparse basis of process._code_basis.
# Acting on every tensor word avoids replacing it by a fitted five-space.
def word_permutations(word):
    return set(permutations(word))


P_SUPPORTS = [set((j, j, j, j) for j in range(4))]
P_SUPPORTS += [word_permutations(a) | word_permutations(b) for a, b in (
    ((0, 0, 1, 1), (2, 2, 3, 3)),
    ((0, 0, 2, 2), (1, 1, 3, 3)),
    ((0, 0, 3, 3), (1, 1, 2, 2)),
)]
P_SUPPORTS += [word_permutations((0, 1, 2, 3))]
index_map = (1, 2, 0, 3)
P5 = sp.zeros(5)
for col, support in enumerate(P_SUPPORTS):
    image = {tuple(index_map[j] for j in word) for word in support}
    targets = [row for row, target in enumerate(P_SUPPORTS) if image == target]
    require(len(targets) == 1, "native quartic support is c-stable")
    P5[targets[0], col] = 1       # i**4=1
require(sp.factor(P5.charpoly(X).as_expr()) == PHI1**3 * PHI3,
        "actual quartic5 characteristic polynomial")
require(vanishes(P5, X**3-1)
        and not vanishes(P5, PHI1) and not vanishes(P5, PHI3),
        "actual quartic5 minimal polynomial")

# The actual erasure20 is (bar V4) tensor P5.  Independently, the twenty
# displayed derivatives of the same five quartics span all Sym3(V4), so this
# is not a dimension-only replacement.
z = sp.symbols("z0:4")
quartics = [sum(q**4 for q in z),
            z[0]**2*z[1]**2 + z[2]**2*z[3]**2,
            z[0]**2*z[2]**2 + z[1]**2*z[3]**2,
            z[0]**2*z[3]**2 + z[1]**2*z[2]**2,
            z[0]*z[1]*z[2]*z[3]]
cubic_monomials = [sp.prod(z[j] for j in word)
                   for word in combinations_with_replacement(range(4), 3)]
derivative_coefficients = sp.zeros(20)
for col, q in enumerate(sp.diff(p, za) for p in quartics for za in z):
    qp = sp.Poly(q, *z)
    for row, monomial in enumerate(cubic_monomials):
        derivative_coefficients[row, col] = qp.coeff_monomial(monomial)
require(derivative_coefficients.det() != 0, "actual twenty erasures span Sym3 V4")
ERASURE20 = sp.kronecker_product(P5, C4_COMPLEX.conjugate())
require(sp.factor(ERASURE20.charpoly(X).as_expr(), extension=I)
        == (X+I)**8 * (X**2-I*X-1)**6,
        "actual erasure20 characteristic polynomial")
require(vanishes(ERASURE20, (X+I)*(X**2-I*X-1))
        and not vanishes(ERASURE20, X+I)
        and not vanishes(ERASURE20, X**2-I*X-1),
        "actual erasure20 minimal polynomial")


# ------------------------------------------- actual mu3 Burau/H1 carrier

T = sp.symbols("T")


def burau_gen(i, n=4):
    """Literal reduced Burau generator used by v597/v610/v613."""
    m = n - 1
    M = sp.eye(m)
    if i - 2 >= 0:
        M[i - 2, i - 1] = T
    M[i - 1, i - 1] = -T
    if i < m:
        M[i, i - 1] = 1
    return M


BUR_T = sp.expand(burau_gen(1) * burau_gen(2) * burau_gen(3))
burau_char_reduced = sp.rem(sp.Poly(BUR_T.charpoly(X).as_expr(), T),
                            sp.Poly(T**2+T+1, T)).as_expr()
require(sp.expand(burau_char_reduced
                  - (X**3 + T*X**2 + (-T-1)*X + 1)) == 0,
        "actual rank3 Eisenstein Burau characteristic polynomial")


def eisenstein_block(a, b):
    """Multiplication by a+b*omega in the integral basis (1,omega)."""
    return sp.Matrix([[a, -b], [b, a-b]])


BUR6 = sp.zeros(6)
for a in range(3):
    for b in range(3):
        reduced = sp.rem(sp.Poly(BUR_T[a, b], T), sp.Poly(T**2+T+1, T))
        BUR6[2*a:2*a+2, 2*b:2*b+2] = eisenstein_block(reduced.nth(0), reduced.nth(1))
DECK6 = sp.diag(*([eisenstein_block(0, 1)] * 3))
PLAIN6 = BUR6**9
require(BUR6**4 == DECK6 and BUR6**12 == sp.eye(6),
        "actual Burau r fourth power is oriented omega deck")
require(sp.factor(BUR6.charpoly(X).as_expr()) == PHI6 * PHI12,
        "retained rational Burau characteristic polynomial")
require(vanishes(BUR6, PHI6*PHI12)
        and not vanishes(BUR6, PHI6) and not vanishes(BUR6, PHI12),
        "retained rational Burau minimal polynomial")
require(sp.gcd(X**3-1, PHI6*PHI12) == 1,
        "quartic5 has no common rational clock constituent with H1")
require(sp.gcd((X+1)*PHI6, PHI6*PHI12) == PHI6,
        "exterior6 shares exactly Phi6 with H1")
require(sp.gcd(PHI4*PHI12, PHI6*PHI12) == PHI12,
        "erasure20 realification shares exactly Phi12 with H1")


# ---------------------------- maximal canonical common transport component

# These are saturated integral cyclic bases of ker Phi12(c) and ker Phi12(r).
# In both bases the generator is the same companion matrix.  Thus the map
# B_NATIVE*a -> B_CURVE*a is an explicit integral Z[C12]-module isomorphism.
B_NATIVE = sp.Matrix([
    [-1, 0, -1, 0],
    [ 0, 0,  0, 1],
    [ 1, 0,  0, 0],
    [ 0,-1,  0,-1],
    [ 0, 0,  1, 0],
    [ 0, 1,  0, 0],
    [ 0, 0,  0, 0],
    [ 0, 0,  0, 0],
])
B_CURVE = sp.Matrix([
    [0,1,0,0],
    [0,1,0,1],
    [1,0,1,0],
    [0,0,1,0],
    [0,1,0,1],
    [0,0,0,1],
])
COMP12 = sp.Matrix([[0,0,0,-1], [1,0,0,0], [0,1,0,1], [0,0,1,0]])
require(matrix_polynomial(C8, PHI12)*B_NATIVE == sp.zeros(8, 4)
        and matrix_polynomial(BUR6, PHI12)*B_CURVE == sp.zeros(6, 4),
        "displayed bases lie in the canonical Phi12 kernels")
require(C8*B_NATIVE == B_NATIVE*COMP12
        and BUR6*B_CURVE == B_CURVE*COMP12,
        "explicit integral Phi12 intertwiner")


def content_of_maximal_minors(B):
    values = [abs(int(B[list(rows), :].det()))
              for rows in combinations(range(B.rows), B.cols)]
    return gcd(*values)


require(content_of_maximal_minors(B_NATIVE) == 1
        and content_of_maximal_minors(B_CURVE) == 1,
        "both Phi12 cyclic lattices are saturated")
# U_E has its last complex coordinate zero, so the half-integral E8 coset
# cannot meet it.  Its displayed integer generators all have even coordinate
# sum.  Thus the saturated Z8 intersection just checked is also U_E cap L_E8.
require(B_NATIVE[6:8, :] == sp.zeros(2, 4)
        and all(sum(B_NATIVE[:, j]) % 2 == 0 for j in range(4)),
        "native Phi12 lattice is the actual E8-lattice intersection")
J_COMMON = sp.Matrix([[0,1,0,1], [0,0,1,0], [0,-1,0,0], [-1,0,-1,0]])
SIGMA_COMMON = sp.Matrix([[-1,0,-1,0], [0,-1,0,-1], [1,0,0,0], [0,1,0,0]])
require(J8*B_NATIVE == B_NATIVE*J_COMMON
        and PLAIN6*B_CURVE == B_CURVE*J_COMMON,
        "Phi12 transport sends compiler J to plain quarter rotation")
require(SIGMA8*B_NATIVE == B_NATIVE*SIGMA_COMMON
        and DECK6*B_CURVE == B_CURVE*SIGMA_COMMON,
        "Phi12 transport sends compiler sigma to oriented deck")

# Retain the existing integral forms.  The Burau invariant Hermitian matrix
# called J in v597 gives E=Tr((1+2 omega)J_B)/6 on the rationalized
# Eisenstein lattice.
def emul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]-a[1]*b[1]


def econj(a):
    return a[0]-a[1], -a[1]


def etrace(a):
    return 2*a[0]-a[1]


H_BUR = [[(0,0), (2,2), (0,2)],
         [(0,-2), (2,0), (2,2)],
         [(-2,-2), (0,-2), (0,0)]]
E_BUR = sp.zeros(6)
ebasis = ((1,0), (0,1))
for aa in range(6):
    i0, i1 = divmod(aa, 2)
    for bb in range(6):
        j0, j1 = divmod(bb, 2)
        h = emul(emul(econj(ebasis[i1]), H_BUR[i0][j0]), ebasis[j1])
        E_BUR[aa, bb] = sp.Rational(etrace(emul((1,2), h)), 6)
require(E_BUR == -E_BUR.T and E_BUR.det() == 1
        and all(v.q == 1 for v in E_BUR)
        and BUR6.T*E_BUR*BUR6 == E_BUR,
        "primitive unimodular Burau alternating form")
E_NATIVE12 = B_NATIVE.T*J8*B_NATIVE
E_CURVE12 = B_CURVE.T*E_BUR*B_CURVE
require(E_NATIVE12.det() == 9 and E_CURVE12.det() == 4,
        "restricted retained alternating-form determinants")
require(gcd(*(abs(int(v)) for v in E_NATIVE12 if v)) == 1
        and gcd(*(abs(int(v)) for v in E_CURVE12 if v)) == 1,
        "restricted alternating-form first elementary divisor")
# Hence the types are (1,3) and (1,2): the module isomorphism is not a
# polarized integral isomorphism.  This obstruction does not erase the exact
# generator, J/plain and sigma/deck intertwining above.


# ------------------------------------------------------- marked clock lifts

Q3 = (1, 0)       # plain algebraic lift (x,w) -> (i x,w)
D3 = (0, 1)       # deck (x,w) -> (x,omega w)
R3 = (1, 1)       # original non-plain lift (x,w) -> (i x,omega w)
require(R3 == ((Q3[0] + D3[0]) % 4, (Q3[1] + D3[1]) % 3), "R3=Q3 D3")
require(pair_power(R3, 4, 3) == D3, "R3^4=D3")
require(pair_power(R3, 9, 3) == Q3, "R3^9=Q3")
require(pair_order(R3, 3) == 12 and pair_order(Q3, 3) == 4,
        "non-plain versus plain mu3 lift")


def lifts_descending_to_R3(n):
    """Lifts (q,Delta_n^k) whose C3 image is R3."""
    require(n % 3 == 0, "a C3 quotient is required")
    return [(1, k) for k in range(n) if k % 3 == 1]


def compiler_marked_lift(n):
    """The unique lift also satisfying ninth power = plain quarter lift."""
    candidates = lifts_descending_to_R3(n)
    selected = [r for r in candidates if pair_power(r, 9, n) == (1, 0)]
    require(len(selected) == 1, f"unique compiler-marked lift on C_{n}")
    return candidates, selected[0]


CANDIDATES6, CLOCK6 = compiler_marked_lift(6)
CANDIDATES12, CLOCK12 = compiler_marked_lift(12)
require(CANDIDATES6 == [(1, 1), (1, 4)] and CLOCK6 == (1, 4),
        "C6 lift census")
require(CANDIDATES12 == [(1, 1), (1, 4), (1, 7), (1, 10)]
        and CLOCK12 == (1, 4), "C12 lift census")
require(pair_power(CLOCK6, 4, 6) == (0, 4)
        and pair_power(CLOCK12, 4, 12) == (0, 4), "fourth power is mu3 deck")
require(pair_order(CLOCK6, 6) == pair_order(CLOCK12, 12) == 12,
        "common-cover clocks have exact order12")


def project_clock(pair, target_modulus):
    return (pair[0] % 4, pair[1] % target_modulus)


require(project_clock(CLOCK12, 6) == CLOCK6, "C12 clock descends to C6")
require(project_clock(CLOCK12, 4) == (1, 0), "C12 clock descends plainly to C4")
require(project_clock(CLOCK12, 3) == R3, "C12 clock descends to original C3 lift")
require(project_clock(CLOCK12, 2) == (1, 0), "C12 clock descends plainly to C2")
require(project_clock(CLOCK6, 3) == R3 and project_clock(CLOCK6, 2) == (1, 0),
        "C6 clock descents")

# The cover alone has the finite lift ambiguity recorded in the predecessor.
parity_lifts = [(1, k) for k in range(2)]
kernel_lifts = [(1, k) for k in range(4)]
require(all(pair_order(r, 2) == 4 for r in parity_lifts), "two parity C4 lifts")
require(all(pair_order(r, 4) == 4 for r in kernel_lifts), "four kernel C4 lifts")


# ---------------------------------------------------- v623 seam realization

N = 16
MARK_BONDS = (0, 4, 8, 12)


def successor(point):
    a, sheet = point
    a2 = (a + 1) % N
    if a2 in MARK_BONDS:
        sheet = (sheet + 1) % 3
    return a2, sheet


def seam_clock(point):
    for _ in range(4):
        point = successor(point)
    return point


def iterate(f, point, n):
    for _ in range(n):
        point = f(point)
    return point


seam_points = [(a, s) for a in range(16) for s in range(3)]
require(all(iterate(seam_clock, p, 4) == (p[0], (p[1] + 1) % 3)
            for p in seam_points), "v623 L^4=deck on all48 points")
require(all(iterate(seam_clock, p, 12) == p for p in seam_points),
        "v623 L^12=1 on all48 points")
require(all(any(iterate(seam_clock, p, d) != p for p in seam_points)
            for d in divisors(12) if d < 12), "v623 L exact order12")


def cycle_lengths_on(items, f):
    index = {item: k for k, item in enumerate(items)}
    perm = [index[f(item)] for item in items]
    used = set()
    counts = Counter()
    for start in range(len(items)):
        if start in used:
            continue
        j = start
        length = 0
        while j not in used:
            used.add(j)
            length += 1
            j = perm[j]
        counts[length] += 1
    return counts


seam_cycle_census = cycle_lengths_on(seam_points, seam_clock)
require(seam_cycle_census == Counter({12: 4}), "v623 is four regular C12 orbits")
# Root c: 19 regular C12 orbits plus three C4 orbits.  Therefore
#   chi_root=(X^12-1)^19 (X^4-1)^3,
#   chi_seam=(X^12-1)^4.
# Both minimal polynomials are X^12-1, but their permutation characters differ.
root_cyclotomic_multiplicity = {1: 22, 2: 22, 3: 19, 4: 22, 6: 19, 12: 19}
seam_cyclotomic_multiplicity = {d: 4 for d in (1, 2, 3, 4, 6, 12)}
require(sum(sp.totient(d)*m for d, m in root_cyclotomic_multiplicity.items()) == 240
        and sum(sp.totient(d)*m for d, m in seam_cyclotomic_multiplicity.items()) == 48,
        "root and seam permutation characteristic degrees")
require(root_cyclotomic_multiplicity != seam_cyclotomic_multiplicity,
        "root and seam permutation modules are not isomorphic")

# There are two exact but different uses of J.  Peripheral all-J lattice
# monodromy is the deck generator Delta4=(0,1) of the C4 kernel cover.  The
# compiler relation J=c^9 transports to the plain lift Q4=(1,0).  A ninth
# power of any lift covering q still has base component 1, so it cannot be a
# deck map, whose base component is 0.  No single cover action can silently
# identify these roles.
PERIPHERAL_J_DECK4 = (0, 1)
COMPILER_J_PLAIN4 = project_clock(CLOCK12, 4)
require(PERIPHERAL_J_DECK4 != COMPILER_J_PLAIN4
        and COMPILER_J_PLAIN4 == (1, 0), "peripheral-J versus compiler-J roles")
require(all(pair_power((1, k), 9, 4)[0] == 1 for k in range(4)),
        "ninth power of a q-lift cannot be a deck transformation")


print("PASS common covers: C6 is minimal for C2/C3; C12 is minimal for C4/C3")
for n in (2, 3, 4, 6, 12):
    r = COVERS[n]
    print("C%d genus=%d peripheral=%s infinity=%d point(s), e=%d" %
          (n, r["genus"], r["peripheral_exponents"],
           r["infinity_points"], r["infinity_ramification"]))
print("PASS compiler: c=J sigma, c^4=sigma, c^9=J, root cycles=%s" %
      dict(sorted(cycle_census.items())))
print("PASS clock dictionary: J<->Q, sigma<->deck, c<->Q deck")
print("PASS actual carriers: real8 Phi4^2 Phi12; retained Burau/H1 model Phi6 Phi12")
print("PASS maximal transport: saturated integral Phi12 modules, J<->plain, sigma<->deck")
print("PASS retained-form obstruction: native type (1,3), Burau type (1,2)")
print("PASS native descendants: quartic5 Phi1^3 Phi3; wedge6 (X+1)^2 Phi6^2")
print("PASS erasure20 realification: Phi4^8 Phi12^6")
print("PASS descendant overlap: quartic none; wedge Phi6 (2 copies); erasure Phi12 (6 copies)")
print("PASS unique marked lifts: C6 %s; C12 %s" % (CLOCK6, CLOCK12))
print("PASS v623: cycles=%s versus E8-root cycles=%s (characters differ)" %
      (dict(seam_cycle_census), dict(sorted(cycle_census.items()))))
print("PASS role boundary: peripheral J is Delta4; compiler J=c^9 is Q4")
print("ALL COVER DICTIONARY CHECKS PASSED")
