"""Exact tests of Fable 5.2's proposed geometric quarter action."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import prod

weights = (0, 1, 2, 3, 4, 1, 2, 3)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


integer_roots = {
    tuple(2 * a * int(k == i) + 2 * b * int(k == j) for k in range(8))
    for i, j in combinations(range(8), 2)
    for a, b in product((-1, 1), repeat=2)
}
spinor_roots = {
    signs for signs in product((-1, 1), repeat=8) if signs.count(-1) % 2 == 0
}
roots = integer_roots | spinor_roots
assert len(roots) == 240


def root_components(system):
    remaining = set(system)
    sizes = []
    while remaining:
        component = {remaining.pop()}
        while True:
            neighbors = {v for v in remaining if any(dot(u, v) for u in component)}
            if not neighbors:
                break
            component |= neighbors
            remaining -= neighbors
        sizes.append(len(component))
    return sorted(sizes)


angles = {p: dot(weights, p) for p in roots}
finite = Counter(a % 8 for a in angles.values())
finite[0] += 8
assert dict(finite) == {0: 60, 2: 64, 4: 60, 6: 64}
assert root_components(p for p, a in angles.items() if a % 8 == 0) == [12, 40]
assert root_components(p for p, a in angles.items() if a % 4 == 0) == [112]
print("Fable's free-mode quarter: order 4, graded dimensions 60,64,60,64, fixed roots A3+D5.")
print("Fable's half-turn: fixed root system D8.")

alternative = (weights[0] + 8, *weights[1:])
alt_angles = {p: dot(alternative, p) for p in roots}
assert all((alt_angles[p] - angles[p]) % 8 == 0 for p in roots)
continuous = root_components(p for p, a in angles.items() if a == 0)
continuous_alt = root_components(p for p, a in alt_angles.items() if a == 0)
assert continuous != continuous_alt
print("Continuous representatives with identical full E8 quarter:", continuous, continuous_alt)
print("Signed integer weight histograms differ:", Counter(angles.values()) != Counter(alt_angles.values()))

fable_continuous = tuple(k - 2 for k in weights)
fable_other = (fable_continuous[0] + 8, *fable_continuous[1:])
fable_components = []
for representative in (fable_continuous, fable_other):
    actual = {p: dot(representative, p) for p in roots}
    assert all((actual[p] - angles[p]) % 8 == 0 for p in roots)
    fable_components.append(root_components(p for p, a in actual.items() if not a))
assert fable_components == [[6, 40], [2, 6, 12]]
print("Fable's selected continuous representative and an exact same-quarter alternative:", fable_components)

forms = (1, 2, 3)
functions = Counter(k % 4 for k in (0, 1, 2, 3, 4))
trace_free_square = Counter((a + b) % 4 for a, b in combinations_with_replacement(forms, 2))
trace_free_square[0] -= 1
assert trace_free_square != functions
assert Counter({(k + 2) % 4: v for k, v in trace_free_square.items()}) == functions
print("The trace-free symmetric square requires an extra sign character to match the actual pullback.")

# In the original D5+A3 chart the first five entries are twice the D5
# weights; the last four are four times the A3 weights (sum zero).
d5_roots = {
    tuple(2 * a * int(k == i) + 2 * b * int(k == j) for k in range(5))
    for i, j in combinations(range(5), 2)
    for a, b in product((-1, 1), repeat=2)
}
family_roots = {
    tuple(4 * (int(k == i) - int(k == j)) for k in range(4))
    for i, j in permutations(range(4), 2)
}
meridians = [tuple(4 * int(k == i) - 1 for k in range(4)) for i in range(4)]
spinor = {s for s in product((-1, 1), repeat=5) if s.count(-1) % 2 == 0}
d5_vectors = {tuple(2 * s * int(k == i) for k in range(5))
              for i, s in product(range(5), (-1, 1))}
family_six = {tuple(a + b for a, b in zip(u, v)) for u, v in combinations(meridians, 2)}
sectors = [
    [(r, (0,) * 4) for r in d5_roots] + [((0,) * 5, r) for r in family_roots],
    list(product(spinor, meridians)),
    list(product(d5_vectors, family_six)),
]
sectors.append([(tuple(-x for x in r), tuple(-x for x in s)) for r, s in sectors[1]])
assert [len(s) for s in sectors] == [52, 64, 60, 64]

carrier_angles = [Fraction(k, 4) for k in (0, 1, 2, -1, 0)]
family_angles = [Fraction(k, 8) for k in (3, 1, -1, -3)]
comparison = [
    (8 * (dot(carrier_angles, r) / 2 + dot(family_angles, s) / 4), grade)
    for grade, sector in enumerate(sectors) for r, s in sector
]
assert all(phase.denominator == 1 for phase, _ in comparison)
comparison = [(int(phase) % 8, grade) for phase, grade in comparison]
old_phases = Counter(phase for phase, _ in comparison)
old_phases[0] += 8
assert [old_phases[k] for k in range(8)] == [30, 32] * 4
assert all(4 * phase % 8 == 4 * grade % 8 for phase, grade in comparison)
assert any(phase and 2 * phase % 8 for phase, grade in comparison if grade == 0)

# Read both operators in eighth-phase units, and generate their joint action.
joint_powers = {
    tuple((m * phase + 2 * n * grade) % 8 for phase, grade in comparison)
    for m in range(8) for n in range(4)
}
assert len(joint_powers) == 16
print("Original quarter T and the glue character D commute, T^4=D^2, joint group C8×C2.")


def shortest_representatives(weight, order, power):
    shift = [Fraction(power * k, order) for k in weight]
    best = None
    vectors = set()
    for offset in (Fraction(0), Fraction(1, 2)):
        costs = []
        for coordinate in shift:
            parity_costs = []
            for parity in (0, 1):
                target = (-coordinate - offset - parity) / 2
                lower = target.numerator // target.denominator
                choices = [coordinate + offset + parity + 2 * j for j in (lower, lower + 1)]
                minimum = min(x ** 2 for x in choices)
                parity_costs.append((minimum, [x for x in choices if x ** 2 == minimum]))
            costs.append(parity_costs)
        for parities in product((0, 1), repeat=8):
            if sum(parities) % 2:
                continue
            selected = [cost[parity] for cost, parity in zip(costs, parities)]
            energy = sum(x[0] for x in selected) / 2
            if best is None or energy < best:
                best = energy
                vectors.clear()
            if energy == best:
                vectors.update(product(*(x[1] for x in selected)))
    return best, vectors


phase_twists = [shortest_representatives(weights, 4, j)[0] for j in (1, 2, 3)]
old_octant = (0, 2, 4, -2, 0, 4, 2, 0)
assert phase_twists == [Fraction(3, 8), Fraction(1, 2), Fraction(3, 8)]
assert 16 * phase_twists[0] % 4 == 2
assert 64 * shortest_representatives(old_octant, 8, 1)[0] % 8 == 6
print("Cyclic anomaly comparison: glue C4 type 2 mod 4; original C8 type 6 mod 8.")

energy, shortest = shortest_representatives(weights, 4, 1)
assert len(shortest) == 4 and energy == Fraction(3, 8)
assert tuple(Fraction(k, 4) for k in fable_continuous) in shortest
assert {tuple(root_components(p for p in roots if not dot(v, p))) for v in shortest} == {(6, 40)}
assert {dot(a, b) for a, b in combinations(shortest, 2)} == {Fraction(-1, 4)}
print("Four shortest logarithms of D: norm squared 3/4, mutual products -1/4, common D5+A2 root type.")

# Fable's second-round torus comparison, in one fixed D8+ chart.
native = tuple(Fraction(k, 4) for k in (0, 1, 2, 3, 4, 0, 1, 2))
form_number_shift = tuple(Fraction(k, 4) for k in (0, 0, 0, 0, 0, 1, 1, 1))
free = tuple(Fraction(k, 4) for k in weights)
glue_reader = tuple(2 * k for k in form_number_shift)
assert free == tuple(x + y for x, y in zip(native, form_number_shift))
assert all((2 * dot(native, p) - dot(glue_reader, p)) % 1 == 0 for p in roots)
matter = [p for p in roots if dot(glue_reader, p) / 2 % 1 == Fraction(3, 4)]
assert len(matter) == 64
relative = Counter((dot(native, p) - dot(glue_reader, p) / 2) % 1 for p in matter)
assert relative == {0: 32, Fraction(1, 2): 32}
print("Native square vs. glue reader on the actual 64-dimensional grade: 32 agree and 32 differ by sign.")

unit = [tuple(int(i == j) for i in range(8)) for j in range(8)]
simple = [
    tuple(Fraction(sign, 2) for sign in (1, -1, -1, -1, -1, -1, -1, 1)),
    tuple(a + b for a, b in zip(unit[0], unit[1])),
    *(tuple(a - b for a, b in zip(unit[i + 1], unit[i])) for i in range(6)),
]
root_space = [tuple(Fraction(a, 2) for a in p) for p in roots]
rho = (0, 1, 2, 3, 4, 5, 6, 23)
assert all(r in root_space and dot(r, rho) == 1 for r in simple)
highest = max(root_space, key=lambda r: dot(r, rho))
assert dot(highest, rho) == 29


def alcove(x):
    for _ in range(2000):
        wall = next((r for r in simple if dot(r, x) < 0), None)
        if wall is not None:
            amount = dot(wall, x)
        elif dot(highest, x) > 1:
            wall = highest
            amount = dot(wall, x) - 1
        else:
            return x
        x = tuple(v - amount * a for v, a in zip(x, wall))
    raise AssertionError("Affine Weyl reduction did not finish")


native8 = alcove(native)
principal8 = alcove(tuple(Fraction(k, 8) for k in rho))
assert native8 == principal8
assert (
    alcove(free)
    == alcove(glue_reader)
    == alcove(tuple(2 * k for k in native))
    == alcove(tuple(Fraction(k, 4) for k in rho))
)
assert alcove(tuple(4 * k for k in native)) == alcove(tuple(Fraction(k, 2) for k in rho))
print("Affine Weyl reduction confirms principal conjugacy classes at 1/8, 1/4 and 1/2.")

# Recover the vector-weight frame of the fixed D8 of T^2 D^{-1}.
joint_involution = tuple(2 * a - b for a, b in zip(native, glue_reader))
fixed_d8 = {p for p in roots if dot(joint_involution, p) / 2 % 1 == 0}
assert len(fixed_d8) == 112 and root_components(fixed_d8) == [112]
candidate_vectors = {
    tuple(a + b for a, b in zip(p, q))
    for p, q in combinations(fixed_d8, 2) if dot(p, q) == 0
}
vector_weights = sorted(v for v in candidate_vectors if all(dot(v, p) % 8 == 0 for p in fixed_d8))
assert len(vector_weights) == 16
assert all(dot(v, v) == 16 for v in vector_weights)
assert all(dot(v, w) in (0, -16) for v, w in combinations(vector_weights, 2))
assert Counter(dot(native, v) / 4 % 1 for v in vector_weights) == {
    Fraction(j, 8): 2 for j in range(8)
}
assert all((dot(glue_reader, v) / 4 - dot(native, v) / 2) % 1 == Fraction(1, 2)
           for v in vector_weights)

# On each of two 8-dimensional vector blocks pair index j with -j.
# X moves j to j+1 with sign a_j. X and X^5 are skew and semisimple.
signs = [1] * 4 + [-1] * 4
assert all(signs[(-j - 1) % 8] == -signs[j] for j in range(8))
assert prod(signs) == 1
assert all(power % 2 == 1 and 2 * power % 8 == 2 for power in (1, 5))

# X on the first block and 2X on the second give a regular semisimple
# combination. For a spinor root, the cyclotomic basis coefficients
# are +/-1 +/-2, so its evaluation is nonzero for either parity.
assert all(any(a + 2 * b for a, b in zip(first, second))
           for first, second in product(product((-1, 1), repeat=4), repeat=2))
print("A stable Cartan subspace in grade i: two 8-cycles, X and X^5 in each block, T^2=iI.")
print("T is regular there: one regular eigenvector is X in the first block plus 2X in the second.")
