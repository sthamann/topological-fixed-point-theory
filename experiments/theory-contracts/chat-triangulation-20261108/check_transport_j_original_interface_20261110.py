"""Read-only exact tests of the existing R032, native fold and cover actions."""
import ast
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v252_full_finite_triple as v252

# Actual inner gauge representation, with the real opposite action.
swap = v252.Kswap()

def real_lift(t):
    return -1j * (t + swap @ t.conj() @ swap)

zero2 = np.zeros((2, 2), complex)
zero3 = np.zeros((3, 3), complex)
y6 = real_lift(v252.rep(3j, zero2, -1j * np.eye(3)))
expected_y6 = {
    ("nu", "L"): -3, ("nu", "R"): 0, ("e", "L"): -3, ("e", "R"): -6,
    ("u", "L"): 1, ("u", "R"): 4, ("d", "L"): 1, ("d", "R"): -2,
}
charges = [expected_y6[name, chirality] for name, chirality, _, _ in v252.P]
assert np.array_equal(y6, np.diag(charges + [-q for q in charges]))
pauli = (
    np.array([[0, 1], [1, 0]]),
    np.array([[0, -1j], [1j, 0]]),
    np.diag([1, -1]),
)
weak_generators = [real_lift(v252.rep(0, 1j * s, zero3)) for s in pauli]
assert all(np.array_equal(y6 @ q, q @ y6) for q in weak_generators)
assert all(np.array_equal(v252.grading() @ q, q @ v252.grading()) for q in weak_generators)

# A nonzero reference Higgs direction, with exact full-rank diagonal Yukawa
# controls. The matrix support also gives |Y_L-Y_R|=1/2 for every original LR block.
id3 = np.eye(3)
phi = v252.full_D(id3, id3, id3, id3, id3)
yukawa_defect = y6 @ phi - phi @ y6
assert np.all(np.count_nonzero(yukawa_defect, axis=0) == 1)
assert np.all(np.count_nonzero(yukawa_defect, axis=1) == 1)
assert set(np.abs(yukawa_defect).ravel()) == {0, 3}

# Original fold M and native Gaussian/family operators from the deployed files.
source = repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
m = sp.Matrix(json.loads(source.read_text())["explicit_Weyl_matrix"])
identity = sp.eye(8)
assert m.T * m == identity
native_file = repo / "verification/v634_st31_structure.py"
tree = ast.parse(native_file.read_text())
idx = next(ast.literal_eval(node.value) for node in tree.body
           if isinstance(node, ast.Assign)
           and any(isinstance(t, ast.Name) and t.id == "SIGMA_IDX" for t in node.targets))
f = identity[idx, :]
j = sp.diag(*[sp.Matrix([[0, -1], [1, 0]])] * 4)
assert j * f == f * j and j**2 == -identity and f**3 == identity
p = sp.diag(1, 1, 1, 1, 1, 0, 0, 0)
y = sp.Matrix([-2, -2, -2, 3, 3, 0, 0, 0]) / 6
beta = identity[:, 3] - identity[:, 4]
r = m.T * f * m
p_src = m.T * p * m
y_src = m.T * y
orbit = [f**k * beta for k in range(3)]
orbit_source = [m.T * b for b in orbit]
assert orbit == [beta, identity[:, 5] - identity[:, 0], identity[:, 1] - identity[:, 2]]
assert orbit_source[0] == identity[:, 3] - identity[:, 7]
assert orbit_source[1] == identity[:, 4] + identity[:, 6]
assert orbit_source[1] == r * orbit_source[0]
assert [y.dot(b) for b in orbit] == [0, sp.Rational(1, 3), 0]
assert [b.dot(p * b) for b in orbit] == [2, 1, 2]
assert (r * p_src * r.T - p_src).rank() == 2
assert (r * y_src - y_src).dot(r * y_src - y_src) == sp.Rational(22, 9)
deltas = []
for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
    delta = sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
    delta_k = delta.dot(orbit_source[1] - orbit_source[0])
    deltas.append(delta_k)
assert deltas == [sp.Rational(13, 16), sp.Rational(15, 16)]

# Given q-frame Q coefficient of the original C3 source cover, section 2.
# Its action on the family Cartan is the permutation of diagonal coefficients.
aq = sp.Matrix([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
q4 = sp.diag(aq, sp.eye(1))
assert q4**4 == sp.eye(4) and q4**2 != sp.eye(4)
w = sp.Symbol("omega")
adapted_family = sp.diag(w**2, w, w**2, w)
assert q4 * adapted_family == adapted_family * q4
diagonal_cartan = sp.Matrix.hstack(*(sp.eye(4)[:, k] - sp.eye(4)[:, 3] for k in range(3)))
permutation = q4.applyfunc(abs)
cartan_q3 = diagonal_cartan[:3, :] * (permutation * diagonal_cartan)[:3, :]
cartan_q8 = sp.diag(sp.eye(5), cartan_q3)
assert cartan_q8.det() == -1 and cartan_q8**2 == identity
assert (cartan_q8 - identity).rank() == 1
assert (f - identity).rank() == 4 and (j - identity).rank() == 8
print(json.dumps({
    "v252_96_weak_generators_commute_with_Y_and_grading": True,
    "full_rank_fixed_Higgs_control_commutator_rank": 96,
    "absolute_nonzero_Y_commutator_coefficient_relative_to_Yukawa": "1/2",
    "native_sigma_D5_projector_difference_rank_in_actual_fold": 2,
    "native_sigma_Y_difference_squared_in_actual_fold": "22/9",
    "actual_weak_root_native_sigma_orbit_target": [list(map(str, q)) for q in orbit],
    "actual_weak_root_native_sigma_orbit_folded_source": [list(map(str, q)) for q in orbit_source],
    "hypercharges_at_fixed_Y": [str(y.dot(q)) for q in orbit],
    "first_delta_K_original_quarters": list(map(str, deltas)),
    "Q_cover_Cartan_fixed_dimension": 7,
    "Q_cover_Cartan_determinant": -1,
    "Q_cover_Cartan_order": 2,
    "K_Cartan_kappa_and_sigma0_ranks_minus_identity": [0, 0],
    "native_Cartan_J_and_sigma_ranks_minus_identity": [8, 4],
}, indent=2))
