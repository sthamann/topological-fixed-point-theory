"""Exact original E8 pair, Fock norm, and neutral descendant controls."""

from contextlib import redirect_stdout
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import io
import runpy
import sys


plan = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
with redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(plan / "check_fable_origin_rule_20261112.py"))
roots, dot = original["roots"], original["dot"]
glue = original["glue_reader"]


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def bit_mask(weight):
    return sum(1 << i for i, x in enumerate(weight) if x == 1)


def create(i, state):
    bit = 1 << i
    if state & bit:
        return 0, state
    return (-1) ** (state & (bit - 1)).bit_count(), state | bit


def create_pair(i, j, state):
    sign_j, intermediate = create(j, state)
    sign_i, result = create(i, intermediate)
    return sign_j * sign_i, result


p1 = (1,) * 8
p2 = (-1, 1, 1, 1, 1, -1, 1, 1)
gamma = sub(p1, p2)
b = (0, 0, 0, 0, 0, 0, -2, -2)
ends = [(p1, p2), (add(p1, b), add(p2, b))]
for top, bottom in ends:
    assert top in roots and bottom in roots and sub(top, bottom) == gamma
    assert create_pair(0, 5, bit_mask(bottom)) == (1, bit_mask(top))
    assert create_pair(0, 5, bit_mask(top))[0] == 0
assert create_pair(0, 5, 0) == (1, 33)
assert create_pair(0, 5, 33)[0] == 0
print("D8 original pair (1,6): vacuum creates mask 33; repeat zero; both Ramond endpoint signs +1.")

def grade(p):
    return int(4 * (dot(glue, p) / 2 % 1))


for g in (1, 3):
    pairs = sorted({
        tuple(sorted((p, sub(gamma, p))))
        for p in roots if grade(p) == g and sub(gamma, p) in roots
        and grade(sub(gamma, p)) == g
    })
    assert len(pairs) == 8 and len({p for pair in pairs for p in pair}) == 16
    assert all(dot(p, q) == -4 for p, q in pairs)
    # Any original Chevalley signs give the same four-particle norm.
    for signs in product((-1, 1), repeat=8):
        amplitudes = {}
        for i, a in enumerate(signs):
            for j, c in enumerate(signs):
                initial = (1 << (2*j)) | (1 << (2*j+1))
                sign, final = create_pair(2*i, 2*i+1, initial)
                if sign:
                    amplitudes[final] = amplitudes.get(final, 0) + a*c*sign
        assert len(amplitudes) == 28
        assert sum(v*v for v in amplitudes.values()) == 112
    assert Q(112, 8**2) == Q(7, 4)
    print(f"D5/A3 grade {g}: eight disjoint pairs, rank 16; normalized repeat norm squared 7/4.")

# The exact neutral return in longitudinal grade two is the Cartan reflection.
root = tuple(Q(x, 2) for x in gamma)
basis = [tuple(int(i == j) for j in range(8)) for i in range(8)]
reflected = [tuple(x - dot(u, root)*r for x, r in zip(u, root)) for u in basis]
assert all(dot(reflected[i], reflected[j]) == int(i == j)
           for i in range(8) for j in range(8))
assert reflected[0] == tuple(-x for x in basis[5])
assert reflected[5] == tuple(-x for x in basis[0])
assert all(reflected[i] == basis[i] for i in (1, 2, 3, 4, 6, 7))
assert tuple(x - dot(root, root)*r for x, r in zip(root, root)) == tuple(-x for x in root)
print("Grade-two charged descendants: neutral return is the exact Weyl reflection of rank 8; circle derivative is the gamma line.")
