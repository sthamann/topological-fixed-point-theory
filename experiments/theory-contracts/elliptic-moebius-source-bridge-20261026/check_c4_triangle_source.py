"""Exact finite certificate for the recovered C4 triangle descent.

The analytic implication from a regular singular three-point connection to
its local monodromy classes is proved in C4-DESCENT-FINITE-MONODROMY.txt.
This file certifies the rational gauge calculation, residues, local
polynomials, marked finite model, source map, and the older D4 obstruction.
"""

import itertools
import json

import sympy as sp


I = sp.I
z, t, lam = sp.symbols("z t lam")
rt2, rt5 = sp.sqrt(2), sp.sqrt(5)
eye = sp.eye(3)
zero = sp.zeros(3)
checks = []


def require(name, condition):
    if not bool(condition):
        raise RuntimeError(name)
    checks.append(name)


def clean_matrix(matrix):
    return sp.simplify(sp.expand(matrix))


def freeze(matrix):
    return tuple(sp.simplify(sp.expand_complex(sp.expand(x))) for x in matrix)


def enumerate_matrix_group(generators):
    elements = {freeze(eye): eye}
    frontier = [eye]
    while frontier:
        following = []
        for element in frontier:
            for generator in generators:
                product = clean_matrix(generator * element)
                key = freeze(product)
                if key not in elements:
                    elements[key] = product
                    following.append(product)
        frontier = following
    return list(elements.values())


def matrix_order(matrix, bound=24):
    power = eye
    for order in range(1, bound + 1):
        power = clean_matrix(power * matrix)
        if power == eye:
            return order
    return None


# The actual v117 source.
U = sp.diag(1, I, -I)
A0 = sp.Matrix([
    [sp.Rational(1, 2), rt2 / 6, 0],
    [rt2 / 6, sp.Rational(1, 4), rt5 / 12],
    [0, rt5 / 12, sp.Rational(1, 4)],
])
A = zero
for k in range(4):
    A += U**k * A0 * U**(-k) / (z - I**k)
A = clean_matrix(A)
A_expected = sp.Matrix([
    [2*z**3, 2*rt2*z**2/3, 0],
    [2*rt2/3, z**3, rt5*z/3],
    [0, rt5*z/3, z**3],
]) / (z**4 - 1)
require("exact four-residue sum", clean_matrix(A - A_expected) == zero)
require("exact C4 covariance", clean_matrix(I*A.subs(z, I*z) - U*A*U**-1) == zero)


# Short meromorphic gauge: Y=z^R X(t), t=z^4, R=diag(0,1,-1).
R = sp.diag(0, 1, -1)
T = sp.diag(1, z, z**-1)
B_direct = clean_matrix((T**-1 * A * T - T**-1 * sp.diff(T, z)) / (4*z**3))
B0 = sp.Matrix([
    [0, 0, 0],
    [-rt2/6, -sp.Rational(1, 4), -rt5/12],
    [0, 0, sp.Rational(1, 4)],
])
B1 = A0
Binf = -B0 - B1
B_partial = B0/t + B1/(t - 1)
require("exact C4 quotient gauge", clean_matrix(B_direct - B_partial.subs(t, z**4)) == zero)
require("Fuchs residue sum", clean_matrix(B0 + B1 + Binf) == zero)

char0 = sp.factor(B0.charpoly(lam).as_expr())
char1 = sp.factor(B1.charpoly(lam).as_expr())
charinf = sp.factor(Binf.charpoly(lam).as_expr())
require("residue characteristic polynomials",
        char0 == lam*(4*lam - 1)*(4*lam + 1)/16
        and char1 == lam*(3*lam - 2)*(3*lam - 1)/9
        and charinf == lam*(2*lam + 1)**2/4)
require("infinity is semisimple in the short gauge",
        clean_matrix(Binf*(Binf + sp.Rational(1, 2)*eye)) == zero
        and Binf != zero and Binf != -sp.Rational(1, 2)*eye)

# The old source gauge R=diag(0,1,3) and the short gauge differ by
# diag(1,1,t) on the quotient.  This independently reconciles both formulas.
R_old = sp.diag(0, 1, 3)
T_old = sp.diag(1, z, z**3)
B_old_direct = clean_matrix((T_old**-1*A*T_old - T_old**-1*sp.diff(T_old, z))/(4*z**3))
B0_old = sp.Matrix([
    [0, 0, 0],
    [-rt2/6, -sp.Rational(1, 4), 0],
    [0, -rt5/12, -sp.Rational(3, 4)],
])
require("recovered September source gauge",
        clean_matrix(B_old_direct - (B0_old/t + A0/(t-1)).subs(t, z**4)) == zero)
H = sp.diag(1, 1, t)
B_from_old = clean_matrix(H*(B0_old/t + A0/(t-1))*H**-1 + sp.diff(H, t)*H**-1)
require("old and short quotient gauges agree", clean_matrix(B_from_old - B_partial) == zero)


# Exact marked S4 model from v117.
M = sp.Matrix([
    [0, -(1+I)/2, (1-I)/2],
    [-(1+I)/2, -I/2, -sp.Rational(1, 2)],
    [(1-I)/2, -sp.Rational(1, 2), I/2],
])
c0 = U**-1
m1 = M
binf = clean_matrix(m1*c0)
require("marked local orders 4 3 2",
        matrix_order(c0) == 4 and matrix_order(m1) == 3 and matrix_order(binf) == 2)
require("marked oriented product", clean_matrix(m1*c0*binf) == eye)
require("local zero class is inverse deck fibre",
        c0 == sp.diag(1, -I, I)
        and sp.factor(c0.charpoly(lam).as_expr()) == (lam-1)*(lam**2+1))

extended = enumerate_matrix_group([c0, m1])
peripherals = [clean_matrix(U**k*M*U**(-k)) for k in range(4)]
ordinary = enumerate_matrix_group(peripherals)
require("extended marked image has order 24", len(extended) == 24)
require("ordinary four-puncture image has order 12", len(ordinary) == 12)
product = eye
for peripheral in peripherals:
    product = clean_matrix(product*peripheral)
require("four peripheral product", product == eye)


# Unique ordered Nielsen class: a 4-cycle and a 3-cycle whose product is an
# involution form one simultaneous-conjugacy orbit of 24 pairs in S4.
permutations = list(itertools.permutations(range(4)))


def compose(p, q):
    return tuple(p[q[i]] for i in range(4))


def inverse(p):
    return tuple(p.index(i) for i in range(4))


def cycle_type(p):
    remaining = set(range(4))
    lengths = []
    while remaining:
        current = start = next(iter(remaining))
        length = 0
        while current in remaining:
            remaining.remove(current)
            length += 1
            current = p[current]
        lengths.append(length)
    return tuple(sorted(lengths))


cycles4 = [p for p in permutations if cycle_type(p) == (4,)]
cycles3 = [p for p in permutations if cycle_type(p) == (1, 3)]
pairs = [(c, m) for c in cycles4 for m in cycles3
         if cycle_type(compose(m, c)) == (1, 1, 2)]
first = pairs[0]
orbit = {
    (compose(compose(k, first[0]), inverse(k)),
     compose(compose(k, first[1]), inverse(k)))
    for k in permutations
}
require("unique marked 4 3 2 Nielsen orbit", len(pairs) == 24 and set(pairs) == orbit)


# Actual degree-four source coordinate for the spherical finite quotient.
J_actual = 1 - z**-4
require("actual map is C4 invariant", clean_matrix(sp.Matrix([[J_actual.subs(z, I*z)-J_actual]])) == sp.zeros(1, 1))
require("actual map is not D4 invariant", sp.simplify(J_actual.subs(z, 1/z)-J_actual) != 0)
require("actual map labels",
        sp.simplify(z**4*(J_actual-1)) == -1
        and sp.limit(J_actual, z, sp.oo) == 1
        and all(sp.simplify(J_actual.subs(z, I**k)) == 0 for k in range(4)))


# Concise exact replay of the already proved absence of a rational inversion
# lift.  Pole analysis in D4-CONNECTION-LIFT.txt proves this nine-parameter
# polynomial is exhaustive; here we certify that its coefficient system has
# full rank.
names = sp.symbols("a b c d e f g h j")
a, b, c, d, e, f, g, h, j = names
G = sp.Matrix([
    [a*z**2, b*z+c*z**2, d*z+e*z**2],
    [0, f*z, g*z],
    [0, h*z, j*z],
])
inversion_equation = clean_matrix(sp.diff(G, z) - A*G - z**-2*G*A.subs(z, 1/z))
numerators = []
for entry in inversion_equation:
    numerator, _ = sp.together(entry).as_numer_denom()
    polynomial = sp.Poly(sp.expand(numerator), z)
    numerators.extend(polynomial.all_coeffs())
linear_matrix, linear_rhs = sp.linear_eq_to_matrix(numerators, names)
require("old rational inversion obstruction has rank nine",
        linear_rhs == sp.zeros(linear_rhs.rows, 1) and linear_matrix.rank() == 9)


result = {
    "verdict": "PASS_EXACT_REUSED_SOURCE",
    "exact_checks": len(checks),
    "checks": checks,
    "B0": str(B0.tolist()),
    "B1": str(B1.tolist()),
    "Binf": str(Binf.tolist()),
    "residue_characteristic_polynomials": [str(char0), str(char1), str(charinf)],
    "extended_orbifold_group_order": len(extended),
    "ordinary_four_puncture_group_order": len(ordinary),
    "actual_source_map": "1 - z**(-4)",
    "recovered_source": {
        "id": "paper-source:10bdafe2930b29584f91db06",
        "sha256": "574058876ebe23d693a9fce061f68da6fdec63419af78243ac9d457d7858f74b",
        "lines": "203-327",
    },
    "analytic_scope": (
        "The report proves regular-singular reduction, triangle-group "
        "faithfulness, and exact marked conjugacy."
    ),
}
print(json.dumps(result, indent=2, sort_keys=True))
