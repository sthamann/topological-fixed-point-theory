"""A single native adjacent-arc bond on the exact marked Fold energy shell."""

import json
import sys
from itertools import product
from pathlib import Path

import sympy as sp

repo = Path(sys.argv[1]).resolve()
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
prior = json.loads((contract / "direct-fold-word-results.json").read_text())
p0, spinor, _ = [tuple(map(sp.Rational, p)) for p in prior["source_orbit"]]
zeta = (1 + sp.I) / sp.sqrt(2)


def e0(charge):
    """8 K_alpha minus alpha Q on primary vectors."""
    return sum(4 * p**2 + (k - 4) * p for k, p in enumerate(charge))


def moments(charge):
    """4 times the Hermitian and anti-Hermitian readings of one arc bond."""
    value = sp.expand(sum(p * zeta**k for k, p in enumerate(charge)))
    return tuple(map(sp.simplify, (sp.re(value), sp.im(value))))


assert sum(spinor) == 2
assert e0(spinor) == 6
assert moments(spinor) == (sp.sqrt(2), 0)
assert e0(p0) == 4
even = [
    p for p in product((-1, 0, 1), repeat=8)
    if sum(p) == 2 and e0(p) == 6
]
assert len(even) == 7
assert all(moments(p) != moments(spinor) for p in even)
even_moments = [(p, *moments(p)) for p in even]

# Every integer charge cost is nonnegative, |p_k| >= 2 costs >= 8,
# and N >= 1 adds 8. Hence the above enumeration is exact.
for k in range(8):
    assert 4 - abs(k - 4) >= 0
    assert 16 - 2 * abs(k - 4) >= 8

# In the half-integer class the coordinate minima sum to zero. A
# coordinate of absolute value >= 5/2 costs >= 16 above its minimum.
assert sum(1 - sp.Rational(abs(k - 4), 2) for k in range(8)) == 0
assert all(24 - 2 * abs(k - 4) >= 16 for k in range(8))
half = [tuple(sp.Rational(x, 2) for x in p)
        for p in product((-3, -1, 1, 3), repeat=8)
        if sum(p) == 4 and sum(2 * x**2 + (k - 4) * x for k, x in enumerate(p)) == 12]
assert len(half) == 7
half_moments = [(p, *moments(p)) for p in half]
assert [p for p, x, _ in half_moments if x == sp.sqrt(2)] == [spinor]

# The bond's diagonal arc-to-residue entry:
# (U* A U)_kk = conj(U_0k) exp(-i*pi*alpha/4) U_1k.
# Use a formal quarter branch factor t, which cancels.
t = sp.symbols("t", nonzero=True)
for k in range(8):
    assert sp.simplify(t * zeta**k / (8 * t) - zeta**k / 8) == 0

# Test the explicit mixture reproducing both components of the bond mean.
p_a = (1, 0, 0, 0, 0, 0, 1, 0)
p_b = (1, 1, 1, 0, 0, -1, 0, 0)
weight_b = 1 / (2 + sp.sqrt(2))
assert p_a in even and p_b in even
mixed = tuple(sp.simplify((1 - weight_b) * a + weight_b * b)
              for a, b in zip(moments(p_a), moments(p_b)))
assert mixed == moments(spinor)
other_eigenvalues = {x for p, x, _ in even_moments + half_moments if p != spinor}
differences = [sp.simplify((x - sp.sqrt(2))**2) for x in other_eigenvalues]
variance_bound = 3 - 2 * sp.sqrt(2)
gap_from_bound = [sp.simplify(diff - variance_bound) for diff in differences]
assert all(diff.is_nonnegative for diff in gap_from_bound) and 0 in gap_from_bound

X = sp.symbols("X")
projection = sp.prod((X - x) / (sp.sqrt(2) - x) for x in other_eigenvalues)
assert sp.simplify(projection.subs(X, sp.sqrt(2))) == 1
assert all(sp.simplify(projection.subs(X, x)) == 0 for x in other_eigenvalues)

cases = []
for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
    k0 = (e0(p0) + alpha * sum(p0)) / 8
    k1 = (e0(spinor) + alpha * sum(spinor)) / 8
    assert k1 - k0 == sp.Rational(1, 4)
    phase = sp.exp(sp.I * sp.pi * e0(spinor) / 2)
    assert phase == -1
    cases.append({"quarter": str(alpha), "intermediate_energy": str(k1),
                  "relative_gap": str(k1 - k0), "intermediate_phase": int(phase)})

print(json.dumps({
    "arc_bond": [0, 1],
    "total_charge": 2,
    "intermediate_energy_8K_minus_alphaQ": 6,
    "spinor_bond_reading": list(map(str, moments(spinor))),
    "even_shell": [{"p": p, "real_mark": str(x), "imaginary_mark": str(y)}
                   for p, x, y in even_moments],
    "spinor_shell": [{"p": list(map(str, p)), "real_mark": str(x), "imaginary_mark": str(y)}
                     for p, x, y in half_moments],
    "matching_even_mixture": {"a": p_a, "b": p_b, "weight_b": str(weight_b)},
    "unique_spinor_on_whole_shell_with_target_real_mark": True,
    "whole_shell_min_squared_real_mark_difference_off_target": str(variance_bound),
    "quarter_cases": cases,
}, indent=2))
