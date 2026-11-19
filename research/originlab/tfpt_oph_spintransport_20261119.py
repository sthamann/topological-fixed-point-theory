"""Exact finite audit of the OPH port action against the TFPT pair (D, T).

Run with Python 3.9 or newer. Uses only the standard library.

OPH inputs are three rows of Lean/Screen/A5PortAction.lean, blob
a542679c8c59685774fc35b06d694e8327392bfa. The TFPT torus is from
check_fable_origin_rule_20261112.py, Atlas source c318ca8e76e0dfc417784b3a.
The OPH edge Hamiltonian is from model/polynomial_sv/SPEC.md, section 3.
The neutral summand, selected modes, pulse schedule and extra memory model
are supplied. This audit does not evaluate native R113.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, prod


PORTS = tuple(range(12))
NEIGHBORS = (
    (1, 2, 3, 4, 6), (0, 2, 3, 5, 7), (0, 1, 4, 5, 8), (0, 1, 6, 7, 9),
    (0, 2, 6, 8, 10), (1, 2, 7, 8, 11), (0, 3, 4, 9, 10), (1, 3, 5, 9, 11),
    (2, 4, 5, 10, 11), (3, 6, 7, 10, 11), (4, 6, 8, 9, 11), (5, 7, 8, 9, 10),
)
NATIVE_EDGES = {
    (0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3),
    (2, 4), (2, 6), (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9),
    (4, 11), (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11),
}
# From the Lean port labels to model/polynomial_sv/carrier.py.
NATIVE_PORT = (0, 1, 5, 7, 11, 9, 10, 8, 4, 6, 2, 3)
HALF_TURN = (1, 0, 3, 2, 7, 6, 5, 4, 9, 8, 11, 10)
THIRD_TURN = (1, 2, 0, 5, 3, 4, 7, 8, 6, 11, 9, 10)
OTHER_HALF_TURN = (10, 11, 9, 8, 6, 7, 4, 5, 3, 2, 0, 1)
BASE_SIGNS = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def compose(left, right):
    return tuple(left[x] for x in right)


def order(perm):
    power = PORTS
    for step in range(1, 61):
        power = compose(perm, power)
        if power == PORTS:
            return step
    raise ValueError('Permutation order exceeded 60')


def closure():
    group = {PORTS}
    queue = [PORTS]
    for element in queue:
        for move in (HALF_TURN, THIRD_TURN):
            result = compose(move, element)
            if result not in group:
                group.add(result)
                queue.append(result)
    return group


def clifford(left, right):
    """Multiply signed basis monomials in the negative Clifford convention."""
    sign_left, indices_left = left
    sign_right, indices_right = right
    exponent = sum(i >= j for i in indices_left for j in indices_right)
    return sign_left * sign_right * (-1)**exponent, indices_left ^ indices_right


def basis_and_signs():
    a, b = HALF_TURN, OTHER_HALF_TURN
    ab = compose(a, b)
    assert ab == compose(b, a) and len({PORTS, a, b, ab}) == 4

    orbits = {frozenset(move[i] for move in (PORTS, a, b, ab)) for i in PORTS}
    assert sorted(map(len, orbits)) == [4, 4, 4]
    basis, signs = [], []
    for orbit in sorted(orbits, key=min):
        i = min(orbit)
        for a_sign, b_sign in BASE_SIGNS:
            vec = [0] * 16
            for slot, value in zip((i, a[i], b[i], ab[i]), (1, a_sign, b_sign, a_sign*b_sign)):
                vec[slot] = value
            basis.append(tuple(Fraction(v, 2) for v in vec))
            signs.append((a_sign, b_sign))
    for j in range(12, 16):
        basis.append(tuple(Fraction(i == j) for i in range(16)))
        signs.append((1, 1))

    for i, u in enumerate(basis):
        for j, v in enumerate(basis):
            assert sum(x*y for x, y in zip(u, v)) == (i == j)
        for move, eig in zip((a, b), signs[i]):
            assert tuple(u[j] for j in (*move, 12, 13, 14, 15)) == tuple(eig*x for x in u)
    return basis, signs


def profile(signs):
    counts = Counter(signs)
    return tuple(counts[pair] for pair in BASE_SIGNS)


def flip_second_sign(signs, indices):
    return [(a, -b if i in indices else b) for i, (a, b) in enumerate(signs)]


def inner(left, right):
    return sum(a*b for a, b in zip(left, right))


def reflect(vector, line):
    factor = 2 * inner(vector, line) / inner(line, line)
    return tuple(x - factor*y for x, y in zip(vector, line))


def port_and_clifford_checks():
    group = closure()
    assert Counter(map(order, group)) == {1: 1, 2: 15, 3: 20, 5: 24}
    half_turns = {p for p in group if order(p) == 2}
    assert all(not any(i == p[i] for i in PORTS) for p in half_turns)
    assert all(sum(compose(p, r) == compose(r, p) for p in group) == 4 for r in half_turns)
    inversion = PORTS[::-1]
    assert all(compose(inversion, move) == compose(move, inversion) for move in group)
    assert {order(compose(inversion, move)) for move in group} == {2, 6, 10}

    basis, signs = basis_and_signs()
    assert profile(signs) == (7, 3, 3, 3)
    target = (6, 4, 4, 2)
    assert all(profile(flip_second_sign(signs, {i})) != target for i in range(16))
    corrections = [pair for pair in combinations(range(16), 2)
                   if profile(flip_second_sign(signs, set(pair))) == target]
    assert len(corrections) == 21
    assert all({signs[i], signs[j]} == {(1, 1), (-1, -1)} for i, j in corrections)

    # This witness fixes the total-port sum. Use a difference of two orbit sums
    # as the common even line and a common odd line on the third orbit.
    f_even = tuple(a - b for a, b in zip(basis[0], basis[4]))
    f_odd = basis[11]
    total_load = (Fraction(1),)*12 + (Fraction(0),)*4
    assert reflect(reflect(total_load, f_even), f_odd) == total_load
    assert all(x == 0 for x in (*f_even[12:], *f_odd[12:]))
    # Use an orthogonal basis (normalization is unnecessary for sign checks).
    pair = {0, 11}
    basis[0], basis[4] = f_even, tuple(a + b for a, b in zip(basis[0], basis[4]))
    for i, vec in enumerate(basis):
        reflected = reflect(reflect(vec, f_even), f_odd)
        assert reflected == tuple((-x if i in pair else x) for x in vec)
    assert profile(flip_second_sign(signs, pair)) == target
    plus_one = (1, frozenset())
    minus_one = (-1, frozenset())
    d = (1, frozenset(i for i, (a, _) in enumerate(signs) if a == -1))
    old_b = (1, frozenset(i for i, (_, b) in enumerate(signs) if b == -1))
    c = (1, frozenset(pair))
    s = clifford(old_b, c)
    assert [clifford(x, x) for x in (d, old_b, c, s)] == [minus_one]*4
    assert clifford(clifford(d, old_b), clifford(d, old_b)) == minus_one
    assert clifford(d, old_b) != clifford(old_b, d)
    assert clifford(d, s) == clifford(s, d)
    assert clifford(clifford(d, s), clifford(d, s)) == plus_one
    return len(half_turns), len(corrections)


def e8_roots():
    for a, b in combinations(range(8), 2):
        for sa, sb in product((-1, 1), repeat=2):
            yield tuple(sa if j == a else sb if j == b else 0 for j in range(8))
    for signs in product((-1, 1), repeat=8):
        if signs.count(-1) % 2 == 0:
            yield tuple(Fraction(sign, 2) for sign in signs)


def spin_checks():
    # Torus angles in eighth-turn units. The four nontrivial corrected joint
    # eigenspaces have real dimensions 6, 4, 4, 2, in BASE_SIGNS order.
    # Put two quarter planes in (+,-) and one in (-,-). Put one half-turn
    # plane in (+,+), one in (-,+), and act trivially on the remaining planes.
    native_d = (0, 0, 0, 0, 0, 4, 4, 4)
    native_t = (0, 2, 4, 6, 8, 0, 2, 4)
    constructed_t = (4, 0, 0, 2, -2, 4, 0, 2)
    # Both an explicit plane relabelling and the lattice phase table agree.
    relabelling = (2, 0, 4, 1, 3, 7, 5, 6)
    shift = tuple(native_t[i] - t for i, t in zip(relabelling, constructed_t))
    assert all(delta % 8 == 0 for delta in shift) and sum(shift) % 16 == 0
    assert tuple(native_d[i] for i in relabelling) == native_d

    def phase(weight, angles):
        value = sum(x*y for x, y in zip(weight, angles))
        assert int(value) == value
        return int(value) % 8

    roots = [*e8_roots(), *([(0,)*8]*8)]
    joint = Counter((phase(weight, native_d), phase(weight, native_t)) for weight in roots)
    assert joint == Counter((phase(weight, native_d), phase(weight, constructed_t)) for weight in roots)
    degrees = [sum(n for (d, _), n in joint.items() if d == k) for k in (0, 2, 4, 6)]
    assert degrees == [60, 64, 60, 64]
    for k in (2, 6):
        assert sum(n for (d, t), n in joint.items() if d == k and 2*t % 8 == d) == 32
    assert all(4*t % 8 == 2*d % 8 for d, t in joint)
    vector_counts = Counter(((-1)**(d//4), (-1)**(t//2)) for d, t in zip(native_d, constructed_t))
    assert tuple(2*vector_counts[key] for key in BASE_SIGNS) == (6, 4, 4, 2)
    assert comb(10, 2) + comb(6, 2) == 60

    # The same reflection of a common fixed plane switches the D8 spinors.
    # Joint phase spectroscopy alone is blind to that change.
    spinors = [tuple(Fraction(x, 2) for x in signs) for signs in product((-1, 1), repeat=8)]
    plus = [weight for weight in spinors if prod(weight) > 0]
    minus = [weight for weight in spinors if prod(weight) < 0]
    mirror = lambda weight: (-weight[0], *weight[1:])
    def phase_pair(weight):
        return phase(weight, native_d), phase(weight, native_t)
    assert {mirror(weight) for weight in plus} == set(minus)
    assert all(phase_pair(weight) == phase_pair(mirror(weight)) for weight in plus)
    assert Counter(map(phase_pair, plus)) == Counter(map(phase_pair, minus))
    orientation = tuple(sum(prod(weight) for weight in sector) for sector in (plus, minus))
    assert orientation == (Fraction(1, 2), Fraction(-1, 2))
    return degrees, [(d, t, joint[d, t]) for d, t in sorted(joint)], orientation


def memory_checks():
    """A single memory coordinate mediates a quarter-turn on any named plane."""
    dimension = 17
    memory = dimension - 1
    identity = tuple(range(1, dimension + 1))

    def update(slot, return_sign):
        # A fresh write sends both (x, 0) to (0, x), with different return signs.
        mapping = list(identity)
        mapping[slot], mapping[memory] = memory + 1, return_sign * (slot + 1)
        return tuple(mapping)

    def multiply(left, right):
        return tuple(left[abs(v) - 1] if v > 0 else -left[abs(v) - 1] for v in right)

    def word(*operations):
        result = identity
        for move in operations:
            result = multiply(move, result)
        return result

    def inverse(move):
        return tuple(next((i + 1) if v > 0 else -(i + 1) for i, v in enumerate(move) if abs(v) == j)
                     for j in range(1, dimension + 1))

    unsigned = update(0, 1)
    oriented = update(0, -1)
    assert unsigned[:memory] == oriented[:memory]
    assert word(unsigned, unsigned) == identity
    assert word(oriented, oriented)[0] == -1
    assert word(oriented, oriented)[memory] == -(memory + 1)

    # Two repeated writes on the selected pair make the sign correction.
    companion = update(11, -1)
    correction = word(oriented, oriented, companion, companion)
    assert correction == tuple(-v if i in (0, 11) else v for i, v in enumerate(identity))

    def quarter(i, j):
        a, b = update(i, -1), update(j, -1)
        step = word(a, b, inverse(a))
        assert step[i] == -(j + 1) and step[j] == i + 1
        assert all(step[k] == identity[k] for k in range(dimension) if k not in (i, j))
        return inverse(step)

    # Realizations of the comparison torus, with one memory reused by every plane.
    main_t = word(*(quarter(2*i, 2*i + 1)
                    for i, turns in enumerate((2, 0, 0, 1, -1, 2, 0, 1))
                    for _ in range(turns % 4)))
    main_d = word(*(quarter(2*i, 2*i + 1) for i in (5, 6, 7) for _ in range(2)))
    assert main_t[memory] == main_d[memory] == memory + 1
    assert word(main_t, main_d) == word(main_d, main_t)
    assert word(main_t, main_t) != main_d
    assert word(main_t, main_t, main_t, main_t) == identity
    assert word(main_d, main_d) == identity


def native_edge_pulse_checks():
    # On one-load states the OPH blank-memory seam generator is
    # 8 I + (I-SWAP)/2 on the seam and zero off it. For t = pi/2 the
    # common phase is 1 and the antisymmetric phase is -i.
    def add(x, y):
        return x[0] + y[0], x[1] + y[1]

    def sub(x, y):
        return x[0] - y[0], x[1] - y[1]

    def inverse(pulses):
        return [(a, b, 4 - t) for a, b, t in reversed(pulses)]

    def preparation(pairs):
        if len(pairs) == 1:
            return []
        half = len(pairs)//2
        left, right = pairs[:half], pairs[half:]
        crossed = [(a, b, 1) for pair_a, pair_b in zip(left, right) for a, b in zip(pair_a, pair_b)]
        right_prep = preparation(right)
        return preparation(left) + crossed + inverse(right_prep) + [(*right[0], 1)] + right_prep

    def reflection(vector):
        plus = [i for i, x in enumerate(vector[:12]) if x > 0]
        minus = [i for i, x in enumerate(vector[:12]) if x < 0]
        pairs = list(zip(plus, minus))
        prep = preparation(pairs)
        return inverse(prep) + [(*pairs[0], 2)] + prep

    def path(a, b):
        paths = {a: (a,)}
        queue = [a]
        for i in queue:
            for j in NEIGHBORS[i]:
                if j not in paths:
                    paths[j] = (*paths[i], j)
                    queue.append(j)
        return paths[b]

    def localize(pulses):
        result = []
        for a, b, t in pulses:
            route = path(a, b)
            swaps = [(i, j, 2) for i, j in zip(route[:-2], route[1:-1])]
            result.extend(swaps + [(*route[-2:], t)] + list(reversed(swaps)))
        return result

    def run(pulses, vector):
        state = [(Fraction(x), Fraction(0)) for x in vector]
        for a, b, t in pulses:
            for _ in range(t):
                x, y = state[a], state[b]
                total, delta = add(x, y), sub(x, y)
                rotated = (delta[1], -delta[0])
                state[a] = tuple(v/2 for v in add(total, rotated))
                state[b] = tuple(v/2 for v in sub(total, rotated))
        return state

    basis, _ = basis_and_signs()
    even = tuple(a - b for a, b in zip(basis[0], basis[4]))
    odd = basis[11]
    first = localize(reflection(even))
    second = localize(reflection(odd))
    assert len(set(NATIVE_PORT)) == 12
    assert {tuple(sorted((NATIVE_PORT[a], NATIVE_PORT[b]))) for a in PORTS for b in NEIGHBORS[a]} == NATIVE_EDGES
    for pulses in (first, second):
        assert all(tuple(sorted((NATIVE_PORT[a], NATIVE_PORT[b]))) in NATIVE_EDGES for a, b, _ in pulses)
    columns = []
    for i in PORTS:
        test = tuple(Fraction(i == j) for j in PORTS)
        expected = reflect(reflect(test, even[:12]), odd[:12])
        assert run(first + second, test) == [(x, 0) for x in expected]
        columns.append(expected)

    # The homogeneous OPH Hamiltonian on the cell commutes with all port
    # rotations, including its identical reference terms. Our pulse word doesn't.
    symmetry_breaks = [move for move in closure()
                       if any(columns[move[i]][move[j]] != columns[i][j] for i in PORTS for j in PORTS)]
    assert THIRD_TURN in symmetry_breaks
    return len(first + second), sum(t for _, _, t in first + second)


if __name__ == '__main__':
    half_turn_count, basis_correction_count = port_and_clifford_checks()
    phase_dimensions, joint_table, orientation = spin_checks()
    memory_checks()
    edge_blocks, quarter_pulses = native_edge_pulse_checks()
    print('OPH half turns:', half_turn_count, '(6 disjoint transpositions each)')
    print('Joint real parity dimensions: (7, 3, 3, 3) -> (6, 4, 4, 2)')
    print('Minimal coordinate-pair sign corrections in the test basis:', basis_correction_count)
    print('Clifford lifts: anticommutation before the correction, commutation after it')
    print('Dimensions of the four D grades in E8:', phase_dimensions)
    print('Joint phase table (D eighths, T eighths, multiplicity):', joint_table)
    print('D8 half spinors are joint-isospectral; oriented sums:', *orientation)
    print('One signed memory mediates both comparison torus operations and returns unchanged.')
    print('Two completions agree on fresh writes; only the oriented return gives the sign correction.')
    print('OPH edge Hamiltonians implement the correction exactly with supplied local pulses:',
          edge_blocks, 'edge blocks,', quarter_pulses, 'pulses of duration pi/2')
    print('The homogeneous cell Hamiltonian cannot give this correction: the word breaks port symmetry.')
