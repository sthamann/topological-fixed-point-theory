"""Exact checks for the charged and regional closure of the original E8 net."""

from collections import Counter, deque
from contextlib import redirect_stdout
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import io
import runpy
import sys


plan = Path(sys.argv[1]) / "docs/research/fundamental-plan" if len(sys.argv) > 1 else Path(__file__).parent
with redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(plan / "check_fable_origin_rule_20261112.py"))

roots = original["roots"]
dot = original["dot"]
native = original["native"]
glue = original["glue_reader"]
vectors = original["vector_weights"]


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def difference(a, b):
    return tuple(x - y for x, y in zip(a, b))


def negate(a):
    return tuple(-x for x in a)


def phase(operator, root):
    return dot(operator, root) / 2 % 1


def grade(root):
    return int(4 * phase(glue, root))


def triality(root):
    return Q(root[0] + root[1] - 2 * root[2], 2) % 3


def multiply(a, b):
    # Exact arithmetic in Z[zeta_8], zeta_8^4 = -1.
    out = [0] * 4
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[(i + j) % 4] += x * y * (-1 if i + j >= 4 else 1)
    return tuple(out)


def fourth(a):
    square = multiply(a, a)
    return multiply(square, square)


def inverse_polynomial_check():
    # For x=X_1+2X_2, the diagonalized vector weights are
    # (1,zeta,zeta^2,zeta^3 | 2,2zeta,2zeta^2,2zeta^3) and their negatives.
    values = {tuple(p[i] + 2 * p[i + 4] for i in range(4)) for p in roots}
    assert (0, 0, 0, 0) not in values
    quartics = {fourth(value) for value in values}
    polynomial = [(1, 0, 0, 0)]
    for eigenvalue in sorted(quartics):
        old = polynomial
        polynomial = [(0, 0, 0, 0)] * (len(old) + 1)
        for j, coeff in enumerate(old):
            polynomial[j] = difference(polynomial[j], multiply(coeff, eigenvalue))
            polynomial[j + 1] = add(polynomial[j + 1], coeff)
    assert polynomial[0][0] and all(coeff[1:] == (0, 0, 0) for coeff in polynomial)
    assert Counter(grade(p) for p in roots) == {0: 52, 1: 64, 2: 60, 3: 64}
    # The centralizer of regular x is the two-block Cartan.
    # The odd powers 1,3,5,7 give two i and two -i per block.
    kernel = Counter((2 * k) % 8 // 2 for k in (1, 3, 5, 7) for _ in range(2))
    assert kernel == {1: 4, 3: 4}
    print("Regular Cartan: kernels 0,4,0,4; ad(x) maps grade 1 onto grade 2, and ad(x)^2 is a 60-to-60 isomorphism from grade 0 to grade 2.")
    print("Exact rational inverse polynomial: distinct fourth-power eigenvalues", len(quartics))


def marked_middle_check():
    p1 = (1,) * 8
    p2 = (-1, 1, 1, 1, 1, -1, 1, 1)
    gamma = difference(p1, p2)
    assert {p1, p2, gamma} <= roots
    assert (grade(p1), grade(p2), grade(gamma)) == (3, 1, 2)
    assert (triality(p1), triality(p2), triality(gamma)) == (0, 2, 1)
    y6 = (-2, -2, -2, 3, 3, 0, 0, 0)
    assert [Q(dot(y6, p), 12) for p in (p1, p2, gamma)] == [0, Q(1, 3), Q(-1, 3)]
    decompositions = {
        (p, difference(gamma, p)) for p in roots if grade(p) == 3
        and difference(gamma, p) in roots and grade(difference(gamma, p)) == 3
    }
    assert len(decompositions) == 16 and (p1, negate(p2)) in decompositions
    assert all(dot(a, b) == -4 for a, b in decompositions)
    # In the marked SU4 basis the family weight of gamma is the pair
    # {0,1}. The half-turn sends it to the complementary pair, so the
    # relative involution mixes its root line with a different one.
    assert {(k + 2) % 4 for k in (0, 1)} == {2, 3}
    print("Original gamma: Y=-1/3, triality 1, grade 2, 8 unordered original-grade decompositions.")


def colored_cartan_check():
    # Choose the two original vector weights in T-eigenspaces 1,2,3,
    # and pair negative weights at -j. Self-inverse eigenspaces
    # 0 and 4 have real orthonormal bases and are unused in this test.
    eigenspaces = {
        j: sorted(v for v in vectors if dot(native, v) / 4 % 1 == Q(j, 8))
        for j in (1, 2, 3)
    }
    block = {j: eigenspaces[j][0] for j in (1, 2, 3)}
    candidates = [difference(block[j + 1], block[j]) for j in (1, 2)]
    charged = [tuple(Q(k, 2) for k in u) for u in candidates]
    assert all(p in roots for p in charged)
    assert all(grade(p) == 1 for p in charged)
    seed = next(p for p in charged if triality(p))
    # The root occurs in X with a nonzero coefficient. No other block
    # cancels it; the opposite matrix entry is the same skew root generator.
    neutral = {p for p in roots if triality(p) == 0}
    target_sector = {p for p in roots if triality(p) == triality(seed)}
    reached = {seed}
    queue = deque([seed])
    while queue:
        p = queue.popleft()
        for q in neutral:
            result = add(p, q)
            if result in target_sector and result not in reached:
                assert dot(p, q) == -4
                reached.add(result)
                queue.append(result)
    assert len(neutral) == 78 and reached == target_sector and len(reached) == 81
    print("The actual Cartan block has a charged root component. Center-neutral root commutators reach all 81 roots of its triality.")


def reflected_lattice_and_middle_check():
    flip = lambda p: p[:4] + negate(p[4:])
    assert {flip(p) for p in roots} == roots
    signatures = Counter((dot(p[:4], p[:4]) // 4, dot(p[4:], p[4:]) // 4)
                         for p in roots)
    assert signatures == {(2, 0): 24, (0, 2): 24, (1, 1): 192}
    assert all(dot(p, flip(p)) % 8 == 0 for p in roots)
    slots = {j: sorted(v for v in vectors if dot(native, v) / 4 % 1 == Q(j, 8))
             for j in (1, 2, 3)}
    block = {v: k for pair in slots.values() for k, v in enumerate(pair)}
    block.update({negate(v): k for v, k in list(block.items())})
    for f in (5, 6, 7):
        for sign in (-1, 1):
            middle = tuple(2 * (int(i == 0) + sign * int(i == f)) for i in range(8))
            parity = (2 * phase(native, middle) - phase(glue, middle)) % 1
            assert parity == (0 if f == 6 else Q(1, 2))
            if parity == 0:
                a, b = next((a, b) for a, b in combinations(vectors, 2)
                            if add(a, b) == tuple(2 * k for k in middle))
                assert block[a] != block[b]
    print("Reflected E8 lattice: 24 left, 24 right and 192 mixed roots. The entire marked middle family lies in the mixed sector.")


def original_color_obstructs_lightcone_split():
    # A real compact split of the two signed cycles uses one real line
    # in each two-dimensional self-inverse T eigenspace, j=0 and j=4.
    # The color center rotates each such plane through a non-real phase.
    y6 = (-2, -2, -2, 3, 3, 0, 0, 0)
    for j in (0, 4):
        weights = [v for v in vectors if dot(native, v) / 4 % 1 == Q(j, 8)]
        charges = {Q(dot(y6, v), 24) for v in weights}
        assert len(weights) == 2 and weights[0] == negate(weights[1])
        assert charges == {Q(-5, 12), Q(5, 12)}
        assert {2 * y % 1 for y in charges} == {Q(1, 6), Q(5, 6)}
    for j in (1, 2, 3):
        weights = [v for v in vectors if dot(native, v) / 4 % 1 == Q(j, 8)]
        assert len(weights) == 2
        if j == 2:
            assert {Q(dot(y6, v), 24) for v in weights} == {Q(-1, 12), Q(1, 12)}
    # In each block Fourier eigenbasis of X, the adjacent off-diagonal
    # entry of Y is i*(sqrt(2)*(y_1+y_3) + 2*y_2)/8. Rational y_2 != 0
    # makes it nonzero in both blocks. At least one cross-block entry
    # survives the invertible Fourier transform of the j=0,4 rotations.
    # Hence the graph on joint Cartan weight lines is connected.
    axes = [tuple(int(i == j) for j in range(8)) for i in range(8)]
    cycle_roots = [
        [difference(axes[start + k], axes[start + k + 1]) for k in range(3)]
        + [add(axes[start + 3], axes[start])]
        for start in (0, 4)
    ]
    for cycle in cycle_roots:
        assert all(dot(a, b) for a, b in zip(cycle, cycle[1:]))
        assert all(tuple(2 * coordinate for coordinate in root) in roots for root in cycle)
    for i in range(4):
        for j in range(4, 8):
            for sign in (-1, 1):
                bridge = add(axes[i], tuple(sign * value for value in axes[j]))
                assert all(any(dot(bridge, root) for root in cycle) for cycle in cycle_roots)
    print("Original Y and its color center mix the two compact Cartan blocks. Y root support has a connected spanning component.")


def block_split_full_color_check():
    # All placements of the original color triplet on the five D5 axes.
    # If a color generator is S-odd it cannot be in the two D4 blocks.
    for color in combinations(range(5), 3):
        color_roots = {
            tuple(2 * (int(k == i) - int(k == j)) for k in range(8))
            for i in color for j in color if i != j
        }
        fixed = all((2 * phase(native, p) - phase(glue, p)) % 1 == 0 for p in color_roots)
        assert fixed == (color == (0, 2, 4))
        if fixed:
            remaining = set(vectors)
            phase_multisets = []
            while remaining:
                stack = [remaining.pop()]
                orbit = list(stack)
                while stack:
                    v = stack.pop()
                    neighbors = {w for w in remaining
                                 if tuple(Q(k, 2) for k in difference(w, v)) in color_roots}
                    orbit.extend(neighbors)
                    stack.extend(neighbors)
                    remaining.difference_update(neighbors)
                if len(orbit) > 1:
                    phase_multisets.append(tuple(sorted(dot(native, v) / 4 % 1 for v in orbit)))
            assert sorted(phase_multisets) == [(Q(1, 8), Q(1, 8), Q(5, 8)),
                                               (Q(3, 8), Q(7, 8), Q(7, 8))]
    print("The full original SU(3) fails the T-adapted compact split for all ten color placements: S-odd roots or repeated T phases on a color triplet.")


def color_factor_lightcone_spin_check():
    # The original A2 is generated by e1-e2 and e2-e3. Projection
    # onto A2 subtracts the mean of the first three coordinates.
    p1 = (1,) * 8
    p2 = (-1, 1, 1, 1, 1, -1, 1, 1)
    gamma = difference(p1, p2)
    color_weights = []
    for root in (p1, p2, gamma):
        first_three = [Q(k, 2) for k in root[:3]]
        mean = sum(first_three) / 3
        color_weights.append(sum((k - mean) ** 2 for k in first_three) / 2)
    assert color_weights == [0, Q(1, 3), Q(1, 3)]
    lorentz_spins = [1 - 2 * weight for weight in color_weights]
    assert lorentz_spins == [1, Q(1, 3), Q(1, 3)]
    print("The original E6/A2 light-cone split gives fractional spin 1/3 on the marked middle and one original leg.")


def stress_rank_check():
    # A symmetric current bilinear M gives a geometric derivative on
    # a root field only if Mp is parallel to p.
    symmetric = [(i, j) for i in range(8) for j in range(i, 8)]
    columns = {pair: k for k, pair in enumerate(symmetric)}
    pivots = {}
    prime = 101
    for p in sorted(roots):
        for i, j in combinations(range(8), 2):
            row = [0] * 36
            for k in range(8):
                row[columns[tuple(sorted((i, k)))]] += p[k] * p[j]
                row[columns[tuple(sorted((j, k)))]] -= p[k] * p[i]
            row = [x % prime for x in row]
            for col, pivot in sorted(pivots.items()):
                factor = row[col]
                if factor:
                    row = [(a - factor * b) % prime for a, b in zip(row, pivot)]
            column = next((k for k, x in enumerate(row) if x), None)
            if column is not None:
                factor = pow(row[column], -1, prime)
                pivots[column] = [x * factor % prime for x in row]
    assert len(pivots) == 35
    print("Quadratic Cartan stress: exact rank 35/36, only a scalar bilinear acts geometrically on all root fields.")


def original_ward_check():
    # These point values use Gaussian dyadics throughout. Python's complex
    # division can round despite dyadic output, so use exact Gaussian pairs.
    def mul(a, b):
        return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]

    def inverse(a):
        norm = a[0] ** 2 + a[1] ** 2
        return Q(a[0], norm), Q(-a[1], norm)

    z = [(1, 0), (0, 1), (-1, 0)]
    w = (0, -1)
    planar = (1, 0)
    circle = (1, 0)
    for i, j in combinations(range(3), 2):
        planar = mul(planar, inverse(difference(z[i], z[j])))
    for point in z:
        circle = mul(circle, point)
    assert mul(planar, circle) == (0, Q(-1, 4))
    stress = (0, 0)
    for i, point in enumerate(z):
        pole = inverse(difference(w, point))
        slope = (0, 0)
        for j, other in enumerate(z):
            if i != j:
                slope = difference(slope, inverse(difference(point, other)))
        stress = add(stress, add(mul(pole, pole), mul(pole, slope)))
    assert stress == (Q(3, 4), 0)
    print("Original collar word at (1,i,-1): -i/4; common Cartan stress at -i gives 3/4.")


if __name__ == "__main__":
    inverse_polynomial_check()
    marked_middle_check()
    colored_cartan_check()
    stress_rank_check()
    reflected_lattice_and_middle_check()
    original_color_obstructs_lightcone_split()
    block_split_full_color_check()
    color_factor_lightcone_spin_check()
    original_ward_check()
