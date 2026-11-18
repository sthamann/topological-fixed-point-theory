"""Test the zero mode of the forced charge correction on the old Fold packets."""
from collections import Counter
from contextlib import redirect_stdout
from functools import cache
from io import StringIO
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import json
import runpy
import sys

here = Path(__file__).resolve().parent
root = here.parents[2]
sys.argv[1:] = [str(root)]
source = root / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(StringIO()):
    old = runpy.run_path(str(source / "check_ramond_full_return_20261110.py"))
    half = runpy.run_path(str(here / "check_original_aps_halfshift.py"))

sp, Q, m = old["sp"], old["Q"], old["m"]
q, r, b, bad, kappa = (half[name] for name in ("q", "r", "b", "bad_b", "correction"))
axes = sp.Matrix.hstack(old["V"], q)
expansions = [
    [(i, Q.from_sympy(sp.simplify(axes[i, j]))) for i in range(8) if axes[i, j] != 0]
    for j in range(axes.cols)
]


def collect(terms):
    result = {}
    for word, value in terms:
        result[word] = result.get(word, Q.zero) + value
    return {word: value for word, value in result.items() if value}


def to_orthonormal(poly):
    terms = []
    for word, value in poly.items():
        partial = {(): value}
        for n, j in word:
            partial = collect((tuple(sorted(w + ((n, i),))), a*c)
                              for w, a in partial.items() for i, c in expansions[j])
        terms.extend(partial.items())
    return collect(terms)


@cache
def exponential(charge, degree):
    if degree == 0:
        return {(): Q.one}
    # Degree recurrence for exp(sum_n (charge . a[-n]) z^n/n).
    return collect((tuple(sorted(word + ((n, i),))), coeff*c/degree)
                   for n in range(1, degree+1)
                   for word, coeff in exponential(charge, degree-n).items()
                   for i, c in enumerate(charge) if c)


def vertex(poly, charge, momentum, l0_step=0):
    # H drift difference is handled separately. These are lattice L0 zero modes.
    shift = int(l0_step - charge.dot(charge)/2 - charge.dot(momentum))
    charge = tuple(map(int, charge))
    terms = []
    for word, coeff in poly.items():
        for length in range(len(word)+1):
            for gone in combinations(range(len(word)), length):
                a = coeff*prod(-charge[word[i][1]] for i in gone)
                if not a:
                    continue
                degree = sum(word[i][0] for i in gone)
                if degree + shift < 0:
                    continue
                kept = tuple(x for i, x in enumerate(word) if i not in gone)
                terms.extend((tuple(sorted(kept+new)), a*c)
                             for new, c in exponential(charge, degree + shift).items())
    return collect(terms)


def inner(left, right):
    return Q.to_sympy(sum((m.conj(value)*right[word]*
                          prod(n**count*factorial(count) for (n, _), count in Counter(word).items())
                          for word, value in left.items() if word in right), Q.zero))


def gram(left, right):
    return sp.Matrix(len(left), len(right), lambda i, j: inner(left[i], right[j]))


def rank(matrix):
    return matrix.to_DM().convert_to(Q).rank()


results = {}
packets, _ = m.code()
for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
    drift = sp.Matrix([j+alpha-4 for j in range(8)])
    assert drift.dot(kappa) == 0 and r.dot(r) == (q+bad).dot(q+bad)
    initial = [old["vertex_zero"](old["series_response"](
        f, [Q.from_sympy(x) for x in old["cl"](old["V"].T*drift/8)], 4)) for f in packets]
    ortho = list(map(to_orthonormal, initial))
    assert inner(ortho[0], ortho[0]) == old["inner"](initial[0], initial[0])
    correct = [vertex(v, b, q) for v in ortho]
    reversed_ = [vertex(v, bad, q) for v in ortho]
    corrected = [vertex(w, kappa, q+bad) for w in reversed_]
    # All four compositions have the same product cocycle since kappa.bad=-2.
    # The double-pole mode is [u_(1),v_(1)] - [u_(0),v_(2)].
    assert kappa.dot(bad) == -2
    sidebands = []
    for start, ab, expected in zip(ortho, corrected, correct, strict=True):
        ba = vertex(vertex(start, kappa, q), bad, q+kappa)
        ab_shift = vertex(vertex(start, bad, q, -1), kappa, q+bad, 1)
        ba_shift = vertex(vertex(start, kappa, q, 1), bad, q+kappa, -1)
        sideband = collect([*((w, -c) for w, c in ab_shift.items()), *ba_shift.items()])
        sidebands.append(sideband)
        pole = collect([*ab.items(), *((w, -c) for w, c in ba.items()), *sideband.items()])
        assert pole == expected
    assert inner(correct[0], correct[0])/inner(ortho[0], ortho[0]) == sp.Rational(19017, 4277)
    correct_gram = gram(correct, correct)
    corrected_gram = gram(corrected, corrected)
    cross_gram = gram(correct, corrected)
    ratios = []
    for i in range(5):
        overlap = cross_gram[i, i]
        ratio = sp.simplify(overlap*sp.conjugate(overlap)/(correct_gram[i, i]*corrected_gram[i, i]))
        assert 0 <= ratio < 1
        ratios.append(str(ratio))
    combined = gram(correct + corrected, correct + corrected)
    assert [rank(g) for g in (correct_gram, corrected_gram, cross_gram, combined)] == [5, 5, 5, 10]
    assert rank(gram(sidebands, sidebands)) == 5
    results[str(alpha)] = {
        "double_pole_contact_equals_correct_on_all_five_packets": True,
        "required_sideband_rank": rank(gram(sidebands, sidebands)),
        "correct_rank": rank(correct_gram),
        "corrected_rank": rank(corrected_gram),
        "cross_rank": rank(cross_gram),
        "combined_rank": rank(combined),
        "packet_squared_normalized_overlaps": ratios,
        "packet_squared_normalized_overlaps_numeric": [float(sp.sympify(x)) for x in ratios],
    }
print(json.dumps(results, indent=2))
