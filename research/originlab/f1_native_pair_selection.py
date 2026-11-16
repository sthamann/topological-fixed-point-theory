"""Absolute selection cost and four-letter word on the existing v746 NS ring."""

from fractions import Fraction
from itertools import combinations, product
from math import factorial
from pathlib import Path
import runpy

import numpy as np
import sympy as sp

from reconstruction import PAULI, require


RING_SOURCE = Path("verification/v746_phys_gnet_local_functor.py")


def absolute_weight(k):
    """One fermion in every antipodal pair of the half-filled 4k ring."""
    return Fraction(factorial(2*k), factorial(k)*(4*k)**k)


def exact_slater_sum(n):
    """Sum occupation determinants as integer polynomials at a primitive root."""
    pairs = n//2
    numerator = [0]*n
    for bits in product((0, 1), repeat=pairs):
        sites = [site+pairs*bit for site, bit in enumerate(bits)]
        term = [1]+[0]*(n-1)
        for first, second in combinations(sites, 2):
            shift = (first-second) % n
            term = [2*term[j]-term[(j-shift) % n]-term[(j+shift) % n] for j in range(n)]
        for j, value in enumerate(term):
            numerator[j] += value
    x = sp.Symbol("x")
    polynomial = sp.Poly.from_list(numerator[::-1], x)
    remainder = polynomial.rem(sp.Poly(sp.cyclotomic_poly(n, x), x))
    require(remainder.degree() <= 0, "Slater sum is not rational")
    return Fraction(int(remainder.nth(0)), n**pairs)


def controls(root):
    checks = []

    def check(name, condition):
        require(condition, name)
        checks.append(name)

    ring = runpy.run_path(str(root/RING_SOURCE))
    identity, x, y, z = PAULI
    eye4, tensor = sp.eye(4), sp.kronecker_product
    destroy = sp.Matrix([[0, 1], [0, 0]])
    c1, c2 = tensor(destroy, identity), tensor(z, destroy)
    n1, n2 = c1.H*c1, c2.H*c2
    current = sp.I*(c1.H*c2-c2.H*c1)
    odd = n1+n2-2*n1*n2
    a = eye4-2*n1
    b = current+eye4-odd
    check("actual opposite-pair current squares to the singly occupied support", current**2 == odd and odd.rank() == 2)
    check("the local parity anticommutes with the pair current", a*current == -current*a)
    check("padding gives two involutions and the exact unconditioned word", a*a == b*b == eye4
          and a*b*a*b == eye4-2*odd and a*current*a*current == -odd)
    check("conditional word is minus one while its unconditioned values are zero and minus one half", (
        sp.trace(odd*a*b*a*b)/sp.trace(odd), sp.trace(a*b*a*b)/4, sp.trace(a*current*a*current)/4
    ) == (-1, 0, -sp.Rational(1, 2)))

    # The joint unsharp comparison is a positive four-effect measurement.
    ax, bz = x/sp.sqrt(2), z/sp.sqrt(2)
    joint = [(identity+s*ax+t*bz)/4 for s, t in product((1, -1), repeat=2)]
    check("positive joint effects attain the compatible qubit four-word bound", all(
        effect.eigenvals() == {0: 1, sp.Rational(1, 2): 1} for effect in joint
    ) and sum(joint, sp.zeros(2)) == identity and (ax*bz)**2 == -identity/4)

    for n in (4, 8, 12, 16):
        check(f"complete exact Slater occupation sum on ring {n}", exact_slater_sum(n) == absolute_weight(n//4))

    weights = {}
    for n in (48, 96, 192):
        opposite = n//2
        covariance = ring["covariance"](n)
        precision = max(abs(covariance[0, 0]-0.5), abs(covariance[opposite, opposite]-0.5),
                        abs(covariance[0, opposite]))
        check(f"original v746 NS sea on the antipodal pair at N={n}", precision < 1e-12)
        shift = ring["shift"](n)
        pair = {0, opposite}
        step_image = {int(abs(shift[:, site]).argmax()) for site in pair}
        # On the one-particle sector the corresponding occupation projector
        # is just the coordinate projection. Its shifted support is disjoint.
        check(f"one microscopic shift moves the fixed support at N={n}", step_image == {1, opposite+1}
              and pair.isdisjoint(step_image))
        placements = tuple((site, site+opposite) for site in range(opposite))
        placement_sets = {frozenset(placement) for placement in placements}
        check(f"the collection of antipodal constraints is shift invariant at N={n}", all(
            frozenset((site+1) % n for site in placement) in placement_sets for placement in placements))
        half_turn = np.linalg.matrix_power(shift, opposite)[np.ix_([0, opposite], [0, opposite])]
        local_turn = sp.Matrix(half_turn.astype(int).tolist())
        check(f"moving support returns with an oriented half-turn at N={n}", local_turn == x*z
              and local_turn*z*local_turn.H == -z and local_turn*y*local_turn.H == y)
        weight = absolute_weight(n//4)
        weights[str(n)] = {"exact": str(weight), "approx": float(weight)}

    # On the full even sector an exterior one-body sign is parity, whereas
    # an actual even CAR operation with a shared Majorana can retain -I.
    global_parity = tensor(z, z)
    even = (eye4+global_parity)/2
    gamma = (tensor(x, identity), tensor(y, identity), tensor(z, x))
    m, n = (sp.I*gamma[i]*gamma[j] for i, j in ((0, 1), (1, 2)))
    check("even Majorana pair operations give minus one even on the empty cap", m*m == n*n == eye4
          and m*n*m*n == -eye4 and global_parity*m == m*global_parity and global_parity*n == n*global_parity
          and c1.H*c1*n != n*c1.H*c1 and (m*n*m*n)[0, 0] == -1)
    check("full exterior parity is plus one on the even support", global_parity*even == even)

    return {
        "status": "EXACT_SCOPED_CONTROLS_PASSED_F1_OPEN",
        "checks": checks,
        "old_source": str(RING_SOURCE),
        "preselection_local_single_occupation": "1/2",
        "conditional_binary_word": -1,
        "unconditioned_padded_binary_word": 0,
        "unconditioned_current_word": "-1/2",
        "jointly_measurable_four_word_bound_saturation_control": "-1/4",
        "all_opposite_pair_absolute_weights": weights,
        "physical_solution": False,
    }
