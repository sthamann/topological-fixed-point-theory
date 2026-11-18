"""Exact controls for the existing circle's relative CAR defect.

The analytic infinite-source proof is DIRECT-DEFECT-CLOSURE.txt. Finite CAR
windows test restrictions of its vector states, not a cyclic approximation
to an infinite chiral implementer. No physical source identification is made.
"""

from fractions import Fraction as F
from pathlib import Path
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def det(matrix):
    a = [[F(x) for x in row] for row in matrix]
    result = F(1)
    for col in range(len(a)):
        pivot = next((i for i in range(col, len(a)) if a[i][col]), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            result = -result
        result *= a[col][col]
        for row in range(col + 1, len(a)):
            scale = a[row][col] / a[col][col]
            for j in range(col + 1, len(a)):
                a[row][j] -= scale * a[col][j]
            a[row][col] = F(0)
    return result


def cofactor(matrix, i, j):
    minor = [[x for n, x in enumerate(row) if n != j]
             for m, row in enumerate(matrix) if m != i]
    return (-1)**(i+j) * det(minor)


def apply(state, vector, create):
    out = {}
    for bits, amplitude in state.items():
        for site, coefficient in enumerate(vector):
            occupied = bool(bits & (1 << site))
            if occupied == create:
                continue
            sign = (-1)**((bits & ((1 << site)-1)).bit_count())
            target = bits ^ (1 << site)
            out[target] = out.get(target, F(0)) + sign*coefficient*amplitude
    return out


def word_expectation(bits, fs, gs):
    operators = [(f, True) for f in fs] + [(g, False) for g in reversed(gs)]
    state = {bits: F(1)}
    for vector, create in reversed(operators):
        state = apply(state, vector, create)
    return state.get(bits, F(0))


for n in range(-20, 21):
    old = int(n < 0)
    shifted = int(n-1 < 0)
    reverse = int(n+1 < 0)
    need(shifted-old == int(n == 0), "forward polarization defect")
    need(reverse-old == -int(n == -1), "reverse polarization defect")

relative_values = []
for z in (F(0), F(1, 3), F(-1, 2), F(3, 4)):
    for radius in (1, 2, 7):
        ratio = F(1)
        for n in range(-radius, radius+1):
            ratio *= (1-z*int(n < 0))/(1-z*int(n-1 < 0))
        need(ratio == 1/(1-z), "directed relative Fredholm quotient")
    relative_values.append({"z": str(z), "inverse_relative": str(ratio)})

# A FINITE cyclic shift has a compensating edge. Its total defect trace is0.
labels = list(range(-2, 3))
cyclic_difference = [int(labels[(i-1) % len(labels)] < 0)-int(n < 0)
                     for i, n in enumerate(labels)]
need(cyclic_difference == [-1, 0, 1, 0, 0], "finite compensating edge")
need(sum(cyclic_difference) == 0, "finite cyclic index cannot be one")

# Actual finite restrictions of the two infinite sea VECTOR STATES.
labels = [-2, -1, 0, 1]
old_bits, new_bits = 0b0011, 0b0111
fs = [[F(v) for v in row] for row in
      [(1, 2, 1, 0), (0, 1, -1, 2), (2, 0, 1, 1), (1, -1, 0, 1)]]
gs = [[F(v) for v in row] for row in
      [(1, 0, 2, 1), (1, 1, 0, -1), (0, 2, 1, 1), (2, 1, -1, 0)]]
word_data = []
for k in range(1, 5):
    f, g = fs[:k], gs[:k]
    matrix = [[sum(f[i][n]*g[j][n] for n in (0, 1))
               for j in range(k)] for i in range(k)]
    changed = [[matrix[i][j]+f[i][2]*g[j][2]
                for j in range(k)] for i in range(k)]
    old_word = word_expectation(old_bits, f, g)
    new_word = word_expectation(new_bits, f, g)
    need(old_word == det(matrix), "old source CAR/Wick word")
    need(new_word == det(changed), "shifted source CAR/Wick word")
    correction = sum(cofactor(matrix, i, j)*f[i][2]*g[j][2]
                     for i in range(k) for j in range(k))
    need(new_word-old_word == correction, "rank-one full-word response")
    word_data.append({"pairs": k, "old": str(old_word), "new": str(new_word)})

energies = []
for a in (F(1, 4), F(3, 4)):
    old_energy = sum(F(n)+a for i, n in enumerate(labels) if old_bits & (1 << i))
    new_energy = sum(F(n)+a for i, n in enumerate(labels) if new_bits & (1 << i))
    need(new_energy-old_energy == a, "forward relative excitation energy")
    hole_bits = old_bits & ~(1 << labels.index(-1))
    hole_energy = sum(F(n)+a for i, n in enumerate(labels) if hole_bits & (1 << i))
    need(hole_energy-old_energy == 1-a, "reverse hole energy")
    energies.append({"twist": str(a), "particle": str(a), "hole": str(1-a)})

# z=e^-alpha=1/2: an alpha-dependent trivial sector changes the ratio.
weight = F(1, 2)
z_def, z_triv = weight**2, weight
need(z_def/z_triv == weight and z_def/z_triv != weight**2,
     "reweighting needs a constant residual quotient")
need(1+weight != 1-weight and 1/(1-weight) != 1+weight,
     "ordinary fermion trace, supertrace and inverse determinant differ")

result = {
    "verdict": "PASS",
    "directed_relative_determinants": relative_values,
    "full_normal_ordered_word_checks": word_data,
    "excitation_energies": energies,
    "finite_cyclic_compensating_edge": cyclic_difference,
    "nonconstant_trivial_sector_countercontrol": "PASS",
    "physical_EM_response_identification": "OPEN",
    "heat_to_EM_parameter_identification_required_for_direct_defect": False,
    "complete_TFPT_solution": False,
}
payload = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)+"\n"
Path(__file__).with_name("direct-seam-defect-results.json").write_text(payload)
print(payload, end="")
