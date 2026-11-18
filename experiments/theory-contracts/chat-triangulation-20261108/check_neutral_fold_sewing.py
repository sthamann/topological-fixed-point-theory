"""Independent neutral OPE contraction of the marked Fold pair."""

import json
import sys
from itertools import combinations
from pathlib import Path

import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo))
from tfpt_explorer.current_block_geometry import _actual_current_source

contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
previous = json.loads((contract / "direct-fold-word-results.json").read_text())
p0, p1, p2 = (sp.Matrix([sp.Rational(k) for k in row]) for row in previous["source_orbit"])
a, b = p1 - p0, p2 - p1
c = a + b
assert c == sp.Matrix([1, -1, 0, -1, -1, 0, 0, 0])
assert (b.dot(a), b.dot(p0), a.dot(p0), c.dot(c)) == (-2, 0, -2, 4)
assert sum(c) % 2 == 0

fold = repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = sp.Matrix(json.loads(fold.read_text())["explicit_Weyl_matrix"])
ch = _actual_current_source()["chevalley"]


def cocycle(left, right):
    coordinates = []
    for charge in (left, right):
        doubled = tuple(2 * n for n in M * charge)
        assert all(x.q == 1 for x in doubled)
        simple = ch.rs.alpha_coords(tuple(int(x) for x in doubled))
        assert all(x.denominator == 1 for x in simple)
        coordinates.append(simple)
    u, v = coordinates
    return (-1) ** (sum(int(u[i]) * ch.B[i][j] * int(v[j])
                        for i in range(8) for j in range(8)) % 2)


cocycle_sewn = cocycle(b, a) * cocycle(c, p0)
assert cocycle_sewn == -1
assert cocycle_sewn == cocycle(b, p1) * cocycle(a, p0)

# Oscillator contraction from the actual Heisenberg commutator:
# x F'(x) = (-b.a) (sum_{n>=1} x^n) F(x), F(0)=1.
coeffs = [sp.Integer(1)]
for n in range(1, 17):
    coeffs.append(-b.dot(a) * sum(coeffs) / n)
assert coeffs == list(range(1, 18))
x = sp.symbols("x")
kernel = (1 - x) ** b.dot(a)
assert list(reversed(sp.Poly(kernel.series(x, 0, 17).removeO(), x).all_coeffs())) == coeffs

# Exact relative sign in the original marked primary subspace even with an
# arbitrary endpoint compression, without choosing a physical projector.
b1 = sp.SparseMatrix(3, 3, {(1, 0): cocycle(a, p0)})
b2 = sp.SparseMatrix(3, 3, {(2, 1): cocycle(b, p1)})
D = sp.diag(1, -1, 1)
P = sp.Matrix(3, 3, sp.symbols("P:9"))
t1 = sp.diag(*sp.symbols("u:3"))
t2 = sp.diag(*sp.symbols("v:3"))
reader = (b2 * b1).T
first = P * reader * t2 * b2 * D * t1 * b1 * P
last = P * reader * t2 * D * b2 * t1 * b1 * P
assert first == -last

cases = []
even_lookalikes = [(i, j) for i, j in combinations(range(8), 2) if i + j == 6]
eye = sp.eye(8)
assert even_lookalikes == [(0, 6), (1, 5), (2, 4)]


def energy(p, delta):
    return p.dot(p) / 2 + delta.dot(p)


for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
    delta = sp.Matrix([(i + alpha - 4) / 8 for i in range(8)])
    k = [energy(p, delta) for p in (p0, p1, p2)]
    assert k[1] - k[0] == sp.Rational(1, 4)
    for i, j in even_lookalikes:
        root = eye[:, i] + eye[:, j]
        assert energy(root, delta) == k[1]
        assert sum(root) == sum(p1) == 2
        assert 8 * energy(root, delta) - alpha * sum(root) == 6
        assert root[0].q == 1 and p1[0].q == 2
    cases.append({
        "quarter": str(alpha),
        "relative_intermediate_gap_in_K_time": str(k[1] - k[0]),
        "equal_strips_zero_mode_exponent": str(k[1] + k[2] - 2 * k[0]),
    })

# Radial current insertion in |z1| < |w| < |z2|. The first term in
# the outer geometric series is regular at w=0. The inner one is not.
w, z1, z2 = sp.symbols("w z1 z2")
inner_series = p0[0] / w + a[0] * sp.series(1 / (w - z1), z1, 0, 3).removeO()
outer_series = b[0] * sp.series(1 / (w - z2), w, 0, 3).removeO()
cartan = (inner_series + outer_series).expand().coeff(w, -1)
assert cartan == p1[0] == sp.Rational(1, 2)

print(json.dumps({
    "even_vertex_charge": list(map(int, c)),
    "pair_inner_product": int(a.dot(b)),
    "sewn_cocycle": cocycle_sewn,
    "neutral_current_kernel": str(kernel),
    "oscillator_weights_first_17": list(map(int, coeffs)),
    "opposite_phase_for_arbitrary_endpoint_compression": True,
    "even_roots_with_same_intermediate_energy_charge_and_phase": even_lookalikes,
    "intermediate_first_cartan_ward_residue": str(cartan),
    "cases": cases,
}, indent=2))
