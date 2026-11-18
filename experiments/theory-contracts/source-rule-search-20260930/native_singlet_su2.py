"""Exact native E8 SM-singlet + M/U-clock centralizer check.

This checks the source-declared native frame only. It does not identify the
geometric P1 class with the family matrices or select a current state.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp

THEORY = "/Users/stefanhamann/Projekte/tfpt-theoryv4"
sys.path.insert(0, THEORY)

from tfpt_explorer.current_block_geometry import (  # noqa: E402
    _actual_current_source,
    _kappa,
)
from tfpt_explorer.flavor_path_transport import exact_clock_dictionary  # noqa: E402
from tfpt_explorer.source_charge_response import hypercharge_vector  # noqa: E402


source = _actual_current_source()
chev = source["chevalley"]
roots = chev.roots
root_index = chev.ridx
index_root = {index: root for root, index in root_index.items()}
clock = exact_clock_dictionary()
M = sp.Matrix(clock["M"])
U = sp.diag(1, sp.I, -sp.I)
I3, I4, I5 = sp.eye(3), sp.eye(4), sp.eye(5)

# Native SM simple roots used by the existing source implementation.
sm_simple_roots = [
    (2, -2, 0, 0, 0, 0, 0, 0),
    (0, 2, -2, 0, 0, 0, 0, 0),
    (0, 0, 0, 2, -2, 0, 0, 0),
]
sm_roots = sm_simple_roots + [tuple(-x for x in r) for r in sm_simple_roots]
Y = sp.Matrix(hypercharge_vector())

def centralizes_sm(root: tuple[int, ...]) -> bool:
    if sp.simplify((Y.T * sp.Matrix(root) / 2)[0]) != 0:
        return False
    i = root_index[root]
    return all(not chev.bracket(root_index[g], i) for g in sm_roots)


# Find the full SM-singlet root subsystem from native charge and Chevalley data.
sm_singlet_roots = {r for r in roots if centralizes_sm(r)}
a3_roots = {r for r in roots if not any(r[:5])}
c_roots = {
    tuple([-1] * 5 + list(tail))
    for tail in [(-1, -1, -1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]
}
expected_a4_roots = a3_roots | c_roots | {tuple(-x for x in r) for r in c_roots}
assert len(a3_roots) == 12 and len(c_roots) == 4
assert sm_singlet_roots == expected_a4_roots
assert len(sm_singlet_roots) == 20

# Exact root-span and coroot-span ranks are four; Y is nonzero and orthogonal
# to every A4 root, so it supplies the fifth Cartan direction of the full
# SM-centralizer u(5)=su(5)+u(1)_Y.
root_matrix = sp.Matrix(sorted(sm_singlet_roots))
coroot_matrix = sp.Matrix([
    [sp.Rational(chev.hvec_index(r).get(240 + j, 0).numerator,
                 chev.hvec_index(r).get(240 + j, 0).denominator)
     for j in range(8)]
    for r in sorted(sm_singlet_roots)
])
assert root_matrix.rank() == 4
assert coroot_matrix.rank() == 4
assert Y != sp.zeros(8, 1)
assert all(sp.simplify((Y.T * sp.Matrix(r) / 2)[0]) == 0 for r in sm_singlet_roots)
assert sp.Matrix.vstack(root_matrix, Y.T).rank() == 5

# Every root-valued Chevalley bracket of the 20-root subsystem remains in it.
# Cartan-valued brackets occur only in opposite-root pairs and lie in the
# rank-four span of its coroots.
root_bracket_outputs = 0
cartan_bracket_outputs = 0
bracket_leaks = []
for left in sm_singlet_roots:
    for right in sm_singlet_roots:
        for target, coefficient in chev.bracket(root_index[left], root_index[right]).items():
            if target < len(roots):
                root_bracket_outputs += 1
                if index_root[target] not in sm_singlet_roots:
                    bracket_leaks.append((left, right, index_root[target], str(coefficient)))
            else:
                cartan_bracket_outputs += 1
assert root_bracket_outputs == 120
assert not bracket_leaks
assert cartan_bracket_outputs > 0

# The actual source clock on the SU(5) fundamental is 1 + m4, with
# m4=1_(C0)+M on the retained three, plus the already-existing extra singlet.
# Thus its basis order is (external singlet, C0, retained 3): 2*1 + 3.
m4 = sp.diag(1, M)
u4 = sp.diag(1, U)
M5 = sp.diag(1, m4)
U5 = sp.diag(1, u4)
assert sp.simplify(M5.det() - 1) == 0
assert sp.simplify(U5.det() - 1) == 0
assert sp.simplify(M5.H * M5 - I5) == sp.zeros(5)
assert sp.simplify(U5.H * U5 - I5) == sp.zeros(5)
assert M5[0, 0] == U5[0, 0] == M5[1, 1] == U5[1, 1] == 1
assert sp.simplify(M**3 - I3) == sp.zeros(3)
assert U**4 == I3 and U**2 != I3

# M has three distinct eigenlines in the retained 3. A reducible unitary
# three-dimensional representation would have a common invariant line;
# none of the M-eigenlines is U-invariant, proving irreducibility exactly.
omega = -sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2
P_eigen = [
    sp.simplify((I3 + M + M**2) / 3),
    sp.simplify((I3 + omega**-1 * M + omega**-2 * M**2) / 3),
    sp.simplify((I3 + omega * M + omega**2 * M**2) / 3),
]
assert all(sp.simplify(P**2 - P) == sp.zeros(3) for P in P_eigen)
assert all(P.rank() == 1 for P in P_eigen)
assert all(sp.simplify(U * P - P * U) != sp.zeros(3) for P in P_eigen)

# Compute the simultaneous commutant in gl(5,C) exactly.
entries = sp.symbols("x0:25")
X = sp.Matrix(5, 5, entries)
equations = list(M5 * X - X * M5) + list(U5 * X - X * U5)
coefficient_matrix, _ = sp.linear_eq_to_matrix(equations, entries)
commutant_basis = coefficient_matrix.nullspace()
assert len(commutant_basis) == 5
for vector in commutant_basis:
    B = sp.Matrix(5, 5, vector)
    assert B[:2, 2:] == sp.zeros(2, 3)
    assert B[2:, :2] == sp.zeros(3, 2)
    assert B[2:, 2:] == sp.eye(3) * B[2, 2]

# The full SM centralizer is u(5)=su(5)+u(1)_Y. In its clock-fixed algebra,
# U(2) on the two trivial lines and a scalar U(1) on the irreducible 3 remain.
# Intersecting with su(5) imposes tr(A)+3z=0, leaving s(u(2)+u(1)) of dim 4.
traceless_constraint_rank = 1
su5_clock_commutant_dimension = len(commutant_basis) - traceless_constraint_rank
full_clock_commutant_dimension = su5_clock_commutant_dimension + 1  # add u(1)_Y
assert su5_clock_commutant_dimension == 4
assert full_clock_commutant_dimension == 5

# Native C0 root is lambda=(-1)^8/2. It commutes with the SM roots and Y;
# both trivial fundamental lines have clock eigenvalue 1, so E_lambda and
# E_-lambda are fixed. Their Chevalley triple is the unique semisimple su(2).
lambda_root = tuple([-1] * 8)
lambda_vec = sp.Matrix(lambda_root) / 2
assert lambda_root in sm_singlet_roots
assert lambda_root in c_roots
assert tuple(-x for x in lambda_root) in sm_singlet_roots
assert lambda_vec.dot(lambda_vec) == 2
assert sp.simplify((Y.T * lambda_vec)[0]) == 0
assert all(not chev.bracket(root_index[g], root_index[lambda_root]) for g in sm_roots)
lambda_cartan = chev.hvec_index(lambda_root)
kappa_h_lambda = _kappa(source, lambda_cartan, lambda_cartan)
assert kappa_h_lambda == 2
normalized_embedding_index = sp.Rational(kappa_h_lambda, lambda_vec.dot(lambda_vec))
assert normalized_embedding_index == 1

result = {
    "verdict": "UNIQUE_SEMISIMPLE_SU2_IN_THE_DECLARED_NATIVE_SM_CLOCK_CENTRALIZER",
    "scope": (
        "Exact result in the source-declared native E8 root frame, SM embedding, "
        "and actual M/U family matrices. It does not derive the family operator, "
        "the P1-to-E8 map, a current state, a coupling, or full TFPT."
    ),
    "native_sm_singlet_root_system": {
        "root_count": len(sm_singlet_roots),
        "decomposition": "12 A3 roots + 4 C + 4 conjugate C",
        "root_span_rank": root_matrix.rank(),
        "coroot_span_rank": coroot_matrix.rank(),
        "root_valued_chevalley_bracket_outputs": root_bracket_outputs,
        "outside_root_bracket_outputs": len(bracket_leaks),
        "cartan_rank_after_adding_Y": 5,
        "full_sm_centralizer_lie_algebra": "u(5) = su(5) + u(1)_Y",
    },
    "clock_action_on_su5_fundamental": {
        "basis_order": "external singlet, omitted C0, retained family 3",
        "M5": "diag(1,1,M)",
        "U5": "diag(1,1,U)",
        "decomposition": "2 trivial + 3 irreducible",
        "M_eigenline_tests_U_invariant": [False, False, False],
    },
    "clock_fixed_centralizer": {
        "gl5_complex_dimension": len(commutant_basis),
        "gl5_complex_form": "Mat_2(C) direct-sum C*I_3",
        "su5_dimension": su5_clock_commutant_dimension,
        "su5_form": "s(u(2) + u(1)) = su(2) + u(1)",
        "full_sm_centralizer_dimension": full_clock_commutant_dimension,
        "full_form": "u(2) + u(1)_Y",
        "unique_semisimple_factor": "su(2)_lambda",
    },
    "lambda_sl2": {
        "root": [str(x) for x in lambda_vec],
        "norm_squared": str(lambda_vec.dot(lambda_vec)),
        "SM_gauge_and_hypercharge_neutral": True,
        "fixed_by_M_and_U": True,
        "kappa_Hlambda_Hlambda": str(kappa_h_lambda),
        "normalized_embedding_index": str(normalized_embedding_index),
    },
}

Path(__file__).with_name("native_singlet_su2_results.json").write_text(
    json.dumps(result, indent=2, ensure_ascii=False) + "\n"
)
print(json.dumps(result, indent=2, ensure_ascii=False))
