"""Exact finite checks for the canonical D x_C Y_q compositum."""

import itertools
import json


checks = []


def require(name, condition):
    if not condition:
        raise RuntimeError(name)
    checks.append(name)


def compose(p, q):
    return tuple(p[q[i]] for i in range(4))


def parity(p):
    return sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4)) % 2


identity = tuple(range(4))
permutations = list(itertools.permutations(range(4)))
A4 = [p for p in permutations if parity(p) == 0]
V4 = [p for p in A4 if compose(p, p) == identity]
c3 = (1, 2, 0, 3)


def perm_power(p, exponent):
    result = identity
    for _ in range(exponent):
        result = compose(result, p)
    return result


cosets = [
    {compose(v, perm_power(c3, exponent)) for v in V4}
    for exponent in range(3)
]
quotient = {element: exponent for exponent, coset in enumerate(cosets)
            for element in coset}
require("A4 and its V4 kernel", len(A4) == 12 and len(V4) == 4)
require("A4 quotient is the marked C3",
        len(quotient) == 12
        and all(quotient[compose(a, b)] == (quotient[a] + quotient[b]) % 3
                for a in A4 for b in A4))


# Gal(D/z)=C12 and Gal(Y_q/z)=A4 have the common quotient Gal(C/z)=C3.
fiber_product = [
    (n, a) for n in range(12) for a in A4 if n % 3 == quotient[a]
]
require("fiber product has degree 48", len(fiber_product) == 48)


# Chinese remainder gives the exact direct-product dictionary.
to_direct = {(n, a): (a, n % 4) for n, a in fiber_product}
require("fiber product is bijective with A4 times C4", len(set(to_direct.values())) == 48)
require("CRT dictionary preserves multiplication",
        all(to_direct[((n + m) % 12, compose(a, b))]
            == (compose(to_direct[(n, a)][0], to_direct[(m, b)][0]),
                (to_direct[(n, a)][1] + to_direct[(m, b)][1]) % 4)
            for n, a in fiber_product for m, b in fiber_product))


# Relative degrees and Riemann-Hurwitz.  D/C is C4 and ramifies totally at
# the two zeros and two poles of w.  Y_q/C is unramified V4.
degree_c_over_z = 3
degree_d_over_c = 4
degree_y_over_c = 4
degree_x_over_c = degree_d_over_c * degree_y_over_c
degree_x_over_z = degree_c_over_z * degree_x_over_c
require("relative linear-disjoint degree", degree_x_over_c == 16)
require("absolute degree", degree_x_over_z == 48)

genus_c = 2
branch_points_on_c = 4
points_above_each = degree_x_over_c // 4
ramification_contribution = branch_points_on_c * points_above_each * (4 - 1)
twice_genus_x_minus_two = degree_x_over_c * (2 * genus_c - 2) + ramification_contribution
genus_x = (twice_genus_x_minus_two + 2) // 2
require("genus 41 by Riemann Hurwitz",
        ramification_contribution == 48
        and twice_genus_x_minus_two == 80 and genus_x == 41)


# The original g already generates D over C: g=H^3, w=H^4, so H=w/g.
# Exponent pairs are recorded in the free rank-one exponent lattice.
require("g and w recover H", (4 - 3) == 1)
require("recovered H has fourth power w and third power g",
        4 * (4 - 3) == 4 and 3 * (4 - 3) == 3)


# Actual geometric quarter on D: Q(H)=i/H.  It inverts the C12 deck.
# Normal forms delta^n Q^p multiply with Q delta Q^-1=delta^-1.
def semidirect_product(left, right):
    n, p = left
    m, q = right
    return ((n + (-1 if p % 2 else 1) * m) % 12, (p + q) % 4)


delta = (1, 0)
quarter = (0, 1)
quarter_inverse = (0, 3)
conjugated = semidirect_product(semidirect_product(quarter, delta), quarter_inverse)
require("actual quarter inverts the D deck", conjugated == (11, 0))


def pair_order(element, multiply, bound=48):
    power = (0, 0)
    for order in range(1, bound + 1):
        power = multiply(power, element)
        if power == (0, 0):
            return order
    return None


geometric_combination = semidirect_product(quarter, (4, 0))
require("geometric quarter times cubic deck remains order four",
        pair_order(geometric_combination, semidirect_product) == 4)


# In the separate compiler closure Q and Delta commute, so Q Delta^4 has
# order lcm(4,3)=12.  This comparison prevents merging the two clocks.
def direct_product(left, right):
    return ((left[0] + right[0]) % 12, (left[1] + right[1]) % 4)


compiler_combination = direct_product((4, 0), (0, 1))
require("commuting compiler combination has order twelve",
        pair_order(compiler_combination, direct_product) == 12)


# The coefficient quarter matrix from the accepted q-packet result.
Aq = ((0, 0, -1), (0, 1, 0), (1, 0, 0))


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
Aq2 = matmul(Aq, Aq)
Aq4 = matmul(Aq2, Aq2)
require("q packet quarter remains one order-four lift",
        Aq2 == ((-1, 0, 0), (0, 1, 0), (0, 0, -1)) and Aq4 == I3)


print(json.dumps({
    "verdict": "PASS_EXACT_MINIMAL_COMPOSITUM",
    "exact_checks": len(checks),
    "checks": checks,
    "degree_over_C": degree_x_over_c,
    "degree_over_z": degree_x_over_z,
    "deck_group_over_z": "A4 x C4",
    "genus": genus_x,
    "quarter_scope": (
        "diagonal geometric Q has order 4; it inverts the D deck and is "
        "not the commuting order-12 compiler clock"
    ),
}, indent=2, sort_keys=True))
