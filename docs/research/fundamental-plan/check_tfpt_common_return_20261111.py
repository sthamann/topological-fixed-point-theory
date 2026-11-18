"""Exact charged-return checks in the two original simultaneous lattice frames."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from pathlib import Path
import sys


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def axis(index, sign=1):
    return tuple(sign * int(i == index) for i in range(8))


def apply(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def klein(left, right):
    return 1 if sum(left[i] * right[j] for i in range(8) for j in range(i)) % 2 == 0 else -1


def integer_part(charge):
    fractional = charge[0] % 1
    assert fractional in (0, Q(1, 2))
    assert all(value % 1 == fractional for value in charge)
    return tuple(int(value - fractional) for value in charge)


def phase_of_word(word, basis):
    prefix = (0,) * 8
    phase = 1
    for root in word:
        r = integer_part(apply(basis, root))
        phase *= klein(r, prefix)
        prefix = add(prefix, r)
    return phase


repo = Path(sys.argv[1])
fold_path = "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
fold = json.loads((repo / fold_path).read_text())
M = tuple(tuple(map(Q, row)) for row in fold["explicit_Weyl_matrix"])
MT = tuple(zip(*M))
identity = tuple(axis(i) for i in range(8))
assert tuple(apply(MT, col) for col in zip(*M)) == identity

roots = {add(axis(i, a), axis(j, b))
         for i, j in combinations(range(8), 2)
         for a, b in product((-1, 1), repeat=2)}
common = roots & {apply(M, root) for root in roots}
blocks = ((0, 1, 5, 6), (2, 3, 4, 7))
assert len(common) == 48
assert common == {root for root in roots if
                  any(all(not root[j] or j in block for j in range(8)) for block in blocks)}

w = (0, 2, 4, -2, 0, 4, 2, 0)
delta = tuple(x // 2 for x in w)
v = axis(7, -1)
dressing = sub(v, delta)
alpha = (Q(1, 2),) * 8
minus_beta = sub(alpha, add(axis(0), axis(5)))
gamma = sub(alpha, minus_beta)
starts = tuple(add(p, delta) for p in (alpha, minus_beta))
targets = tuple(add(p, v) for p in (alpha, minus_beta))
assert gamma in common
assert all(sum(dressing[j] for j in block) % 2 == 0 for block in blocks)
assert [sum(delta[j] for j in block) % 2 for block in blocks] == [0, 1]
assert [sum(v[j] for j in block) % 2 for block in blocks] == [0, 1]

candidates = [root for root in sorted(roots)
              if all(not x or x * y > 0 for x, y in zip(root, dressing))]
counts = Counter()
phases = defaultdict(Counter)
witnesses = {}
for multiset in combinations_with_replacement(candidates, 4):
    if tuple(map(sum, zip(*multiset))) != dressing:
        continue
    for word in sorted(set(permutations(multiset))):
        momenta = starts
        drops = []
        for root in word:
            drops.extend(1 + dot(p, root) for p in momenta)
            momenta = tuple(add(p, root) for p in momenta)
        assert momenta == targets
        if any(drop > 0 for drop in drops):
            continue
        populations = ["canonical"]
        if all(root in common for root in word):
            populations.append("common")
        for population in populations:
            counts[population] += 1
            for name, basis in (("compiler", identity), ("fold", MT)):
                if name == "fold" and population != "common":
                    continue
                coefficient = phase_of_word(word, basis)
                phases[population, name][coefficient] += 1
                witnesses.setdefault((population, name, multiset, coefficient), word)
                total = integer_part(apply(basis, dressing))
                for start in starts if name == "compiler" else ():
                    n = integer_part(apply(basis, start))
                    separate = 1
                    for root in word:
                        r = integer_part(apply(basis, root))
                        separate *= klein(r, n)
                        n = add(n, r)
                    assert separate == coefficient * klein(total, integer_part(apply(basis, start)))

print(f"Common current roots: {len(common)}, D4 plus D4. Both delta and v have block parities (0,1).")
print(f"Minimal nonincreasing words: {dict(counts)}.")
print("Signs versus one fused charge operator: " +
      ", ".join(f"{key} {dict(sorted(value.items()))}" for key, value in sorted(phases.items())))
mixed = [(key, word, witnesses[key[:-1] + (-key[-1],)])
         for key, word in witnesses.items()
         if key[:2] == ("common", "fold") and key[-1] == 1
         and key[:-1] + (-1,) in witnesses]
print(f"Common multisets with opposite signs: {len(mixed)}.")
if mixed:
    key, positive, negative = mixed[0]
    print(f"Original fold cocycle, same multiset, opposite signs: {positive}, {negative}.")
for label, basis in (("compiler", identity), ("fold", MT)):
    r = integer_part(apply(basis, dressing))
    g = integer_part(apply(basis, gamma))
    relative = klein(r, g)
    assert relative * klein(g, r) == 1 and dot(r, g) == -2
    print(f"{label}: (lambda,gamma)= -2; exchange sign +1; section-dependent ratio {relative}.")

def rotate(vector, steps):
    return tuple(vector[(j - steps) % 8] for j in range(8))


native_gamma = apply(MT, gamma)
gamma_periods = [k for k in range(8) if rotate(native_gamma, k) == native_gamma]
assert gamma_periods == [0]
native_v = apply(MT, v)
assert len(set(native_v)) > 1
assert len({value % 1 for value in native_v}) == 2
print(f"Original W1: native middle charge {native_gamma}, periods mod 8 {gamma_periods}.")
print(f"Original W1: desired v in native frame {native_v}. All powers fixing the middle charge shift by k*(1,...,1).")
print("Therefore no integer power of the original same-source W1 hits both marked endpoint charges.")
print(f"After four W1 steps, target middle charge is {apply(M, rotate(native_gamma, 4))}.")
assert counts["canonical"] == 468


def in_d8_dual(charge):
    return len({value % 1 for value in charge}) == 1


def in_e8(charge):
    return in_d8_dual(charge) and sum(charge) % 2 == 0


def in_common_dual(charge):
    return all(len({charge[j] % 1 for j in block}) == 1
               and charge[block[0]] % 1 in (0, Q(1, 2)) for block in blocks)


fold_roots = {apply(M, root) for root in roots}
assert all(in_e8(root) for root in fold_roots)
assert all(not in_d8_dual(apply(MT, p)) for p in (*starts, *targets))
assert all(in_d8_dual(apply(MT, p)) for p in (alpha, minus_beta))
assert all((dot(delta, root) - dot(v, root)) % 1 == 0 for root in fold_roots)
assert {root for root in fold_roots if dot(v, root) % 1 == 0} == common
assert {root for root in fold_roots if dot(v, root) % 1 == Q(1, 2)} == fold_roots - common
print("Both marked initial momenta are in the old fold; both bare and corrected interior momenta are outside its four sectors.")
print("The native half-clock and the corrected target give identical fold current monodromy: + on 48 common roots and - on 64 others.")

fold_vector = apply(M, axis(0))
assert all(in_common_dual(q) for q in (v, fold_vector))
cosets = ((0,) * 8, v, fold_vector, add(v, fold_vector))
assert all(in_e8(add(q, q)) for q in cosets)
assert all(not in_e8(sub(p, q)) for p, q in combinations(cosets, 2))
print("Minimal additive charge group containing both original four-sector charge groups: (D4 plus D4)*.")
print("It has 16 classes over D4 plus D4, four over E8, with E8 representatives 0, v, M e_1, v+M e_1.")

other_half_pattern = tuple((value + Q(1, 2)) % 1 for value in native_v)
half_patterns = {tuple(value % 1 for value in native_v), other_half_pattern}
assert all(sum(bool(x) for x in pattern) == 4 for pattern in half_patterns)
minima = set()
for pattern in half_patterns:
    for signs in product((-1, 1), repeat=4):
        assigned = iter(signs)
        charge = tuple(value * next(assigned) if value else value for value in pattern)
        assert dot(charge, charge) == 1 and in_d8_dual(sub(charge, native_v))
        minima.add(apply(M, charge))
assert len(minima) == 32
canonical_vectors = {axis(i, sign) for i, sign in product(range(8), (-1, 1))}
assert len(minima & canonical_vectors) == 16
assert len(minima - canonical_vectors) == 16
flavor_and_clock_neutral = {p for p in minima if all(x == 0 for x in p[:5])
                            and dot(w, p) % 8 == 0 and dot(p, gamma) == 0}
assert flavor_and_clock_neutral == {v, tuple(-x for x in v)}
print("The missing half-clock fold sector has 32 minimum-weight momenta of weight 1/2.")
print("Only v and -v are neutral for D5 Cartan, clock phase and middle current among these 32; the prescribed direction fixes v.")

assert phases["common", "compiler"] == phases["common", "fold"] == {-1: 12, 1: 12}
assert phases["canonical", "compiler"] == {-1: 232, 1: 236}
assert counts["common"] == 24 and len(mixed) == 1

# The only common multiset has all 24 orders. Two pairs of roots live
# in orthogonal blocks. The cocycle flips when either pair is reversed.
multiset = mixed[0][0][2]
a = sub((0,) * 8, add(axis(1), axis(5)))
d = sub((0,) * 8, add(axis(5), axis(6)))
b = sub((0,) * 8, add(axis(2), axis(7)))
c = sub(axis(3), axis(2))
assert set(multiset) == {a, d, b, c}
assert dot(a, d) == dot(b, c) == 1
assert all(dot(x, y) == 0 for x, y in product((a, d), (b, c)))
assert all(dot(p, sub(a, d)) == dot(p, sub(b, c)) == 0 for p in starts)


def conformal(charge):
    return dot(charge, charge) / 2


def coefficients(time, origin, gap_unit):
    result = Counter()
    for word in permutations((a, d, b, c)):
        p = origin
        exponent = 0
        for root in word[:-1]:
            p = add(p, root)
            exponent += time(p) * gap_unit
        assert exponent.denominator == 1
        result[int(exponent)] += phase_of_word(word, MT)
    return {power: value for power, value in sorted(result.items()) if value}


for p in starts:
    assert coefficients(conformal, p, 32) == {}
print("The 24-order sum cancels with conformal propagation at all three gap lengths: both internal reversals preserve each intermediate conformal weight.")
for quarter in (Q(1, 4), Q(3, 4)):
    d_native = tuple((j + quarter) / 8 - Q(1, 2) for j in range(8))
    d_target = apply(M, d_native)
    assert dot(d_target, v) == Q(1, 4)
    assert dot(d_target, sub(a, d)) == -Q(1, 4)
    assert dot(d_target, sub(b, c)) == Q(3, 4)

    def shifted_time(charge):
        return conformal(charge) + dot(d_target, charge)

    for initial in (alpha, minus_beta):
        assert conformal(add(initial, v)) == conformal(initial)
        assert shifted_time(add(initial, v)) - shifted_time(initial) == Q(1, 4)

    polynomials = [coefficients(shifted_time, p, 32) for p in starts]
    assert all(polynomials)
    for polynomial in polynomials:
        exact = sum(n * Q(1, 2)**power for power, n in polynomial.items())
        assert exact
    print(f"Old quarter {quarter}: transport raises K by 1/4 on both legs; internal order splits are -1/4 and +3/4.")
    print(f"Old quarter {quarter}: exact equal-gap coefficients are nonzero on both legs.")
