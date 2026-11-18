"""Exact coefficient check for the operator-skeleton/algebra ambiguity.

This is a finite matrix identity certificate used in the analytic product proof.
It neither constructs a physical local QFT nor selects a TFPT source.
No Standard Model representation, Yukawa matrix, or fitted coupling is used.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    basis = []
    for i in range(2):
        for j in range(2):
            e = sp.zeros(2)
            e[i, j] = 1
            basis.append(e)

    def matrix_of(op):
        return sp.Matrix.hstack(*[
            sp.Matrix(list(op(x))) for x in basis
        ])

    left = [matrix_of(lambda x, a=a: a*x) for a in basis]
    right = [matrix_of(lambda x, a=a: x*a) for a in basis]
    transpose = matrix_of(lambda x: x.T)
    zero = sp.zeros(4)
    count = 0
    for a, lm in zip(basis, left):
        # J is complex conjugation followed by transpose:
        # J L_(a*) J^-1 = R_a. Real coefficient basis suffices
        # by complex bilinearity; the analytic proof treats adjoints.
        assert transpose*lm.T*transpose == matrix_of(lambda x: x*a)
        for rm in right:
            assert lm*rm-rm*lm == zero
            count += 1

    def vec(a):
        return sp.Matrix(list(a))

    def commutant_dimension(generators):
        # Basis-independent linear equations for XG-GX=0.
        unknowns = sp.symbols('x:16')
        x = sp.Matrix(4, 4, unknowns)
        equations = [entry for g in generators for entry in x*g-g*x]
        coeff, _ = sp.linear_eq_to_matrix(equations, unknowns)
        return 16-coeff.rank()

    scalar_dimension = sp.Matrix.hstack(vec(sp.eye(4))).rank()
    matrix_dimension = sp.Matrix.hstack(*map(vec, left)).rank()
    assert scalar_dimension == 1 and matrix_dimension == 4
    assert commutant_dimension([sp.eye(4)]) == 16
    assert commutant_dimension(left) == 4
    assert commutant_dimension(left+right) == 1
    # The center of the matrix insertion algebra is its intersection
    # with its commutant, computed in its own four-dimensional span.
    c = sp.symbols('c:4')
    z = sum((c[i]*left[i] for i in range(4)), sp.zeros(4))
    equations = [entry for lm in left for entry in z*lm-lm*z]
    coeff, _ = sp.linear_eq_to_matrix(equations, c)
    assert 4-coeff.rank() == 1

    root = Path(__file__).resolve().parents[3]
    source_paths = [
        '_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex',
        '_archive/tfpt-45/source_extracts/02_carrier_source.tex',
        'experiments/theory-contracts/raw-carrier-origin-gate-20260920/PROOF.txt',
        'experiments/theory-contracts/source-variation-origin-gate-20260922/PROOF.md',
    ]
    # An unmarked, supplied E=C5 gives an internal anti-linear Hodge
    # operator. This is NOT a physical-J or local-field identification.
    g = 5
    hodge = sp.zeros(2**g)
    grading = sp.diag(*[(-1)**i.bit_count() for i in range(2**g)])
    for i in range(2**g):
        subset = [k for k in range(g) if i & (1 << k)]
        complement = [k for k in range(g) if not i & (1 << k)]
        inversions = sum(a > b for a in subset for b in complement)
        hodge[(2**g-1)^i, i] = (-1)**inversions
    assert hodge*hodge == sp.eye(32)
    assert hodge*grading == -grading*hodge
    assert (grading*hodge)**2 == -sp.eye(32)
    # The full even internal endomorphism algebra is invariant under
    # Hodge conjugation, hence cannot satisfy the order-zero condition
    # with that J: exhibit a nonzero commutator in one parity block.
    a = sp.zeros(32); a[0, 3] = 1
    b = sp.zeros(32); b[3, 0] = 1
    opposite_seed = hodge*b.T*hodge
    opposite = hodge*opposite_seed.T*hodge
    assert opposite == b
    order_zero_defect = a*opposite-opposite*a
    assert order_zero_defect != sp.zeros(32)

    result = {
        'contract': 'UR.SEAM.JOINT_SOURCE_ORIGIN_AUDIT.01',
        'annex': 'local-algebra-identifiability-20261016',
        'verdict': 'OPERATOR_SKELETON_DOES_NOT_SELECT_INSERTION_ALGEBRA',
        'source_sha256': {
            p: hashlib.sha256((root/p).read_bytes()).hexdigest()
            for p in source_paths
        },
        'matrix_identity_check': {
            'internal_hilbert_space': 'HS(M2(C))',
            'left_right_basis_commutators_zero': count,
            'scalar_insertion_fiber_dimension': scalar_dimension,
            'matrix_insertion_fiber_dimension': matrix_dimension,
            'scalar_insertion_commutant_dimension': 16,
            'matrix_insertion_commutant_dimension': 4,
            'matrix_insertion_center_dimension': 1,
            'joint_left_right_commutant_dimension': 1,
            'opposite_action': 'J L_(a*) J^-1 = R_a',
        },
        'analytic_scope': 'Product Dirac skeleton with zero finite Dirac term; smooth left multipliers, common domain; original broad center/product-collar conditions only.',
        'unmarked_internal_hodge_control': {
            'carrier_rank_granted_not_derived': 5,
            'exterior_dimension': 32,
            'antilinear_hodge_square': 1,
            'hodge_anticommutes_exterior_grading': True,
            'grading_times_hodge_square': -1,
            'physical_real_structure_identified': False,
            'finite_ko6_triple_derived': False,
            'full_even_internal_algebra_order_zero_defect_rank': order_zero_defect.rank(),
            'full_even_internal_algebra_is_coordinate_algebra': False,
        },
        'same_full_original_tuple': False,
        'full_p1_p2_or_minimal_branch_countermodels': False,
        'physical_local_field_net_derived': False,
        'internal_clifford_identified_with_physical_car': False,
        'new_su5_no_selector_result_claimed': False,
        'standard_model_structure_inserted': False,
        'yukawa_coupling_inserted': False,
        'verification_or_ledger_modified': False,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
