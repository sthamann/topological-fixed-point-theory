#!/usr/bin/env python3
"""Exact tangent test on the actual quartic code; no group census.

The compact-group conclusion additionally uses Nebe--Rains--Sloane,
arXiv:math/0001038v2 Theorem 6.5, as explained in transversal_audit.md.
"""
from __future__ import annotations
from collections import defaultdict
import itertools as it
import json
from pathlib import Path
import sys
import sympy as s

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(REPO))
from tfpt_explorer.process import _code_basis
from tfpt_explorer.marker_selection import build_marker_data
from tfpt_explorer.charge_selection import _parse_matrix, _process_frame
from tfpt_explorer.marked_source_process import _source_operators


def clean(m):
    return m.applyfunc(s.simplify)


def main():
    _, supports = _code_basis()
    sizes = list(map(len, supports))
    assert sizes == [4, 12, 12, 12, 24]
    swap = s.zeros(16)
    for left, right in it.product(range(4), repeat=2):
        swap[4*right+left, 4*left+right] = 1
    pair_marginals = []
    for first, second in it.combinations(range(4), 2):
        others = [j for j in range(4) if j not in (first, second)]
        two_body = s.zeros(16)
        for support in supports:
            for left in support:
                for right in support:
                    if all(left[j] == right[j] for j in others):
                        two_body[4*left[first]+left[second], 4*right[first]+right[second]] += s.Rational(1, len(support))
        assert two_body == (s.eye(16)+swap)/4
        pair_marginals.append([first, second])
    I = s.eye(2)
    X = s.Matrix([[0, 1], [1, 0]])
    Y = s.Matrix([[0, -s.I], [s.I, 0]])
    Z = s.diag(1, -1)
    labels = [''.join(p) for p in it.product('IXYZ', repeat=2)]
    hs = [s.kronecker_product(a, b) for a, b in it.product([I, X, Y, Z], repeat=2)]

    # Columns use the exact, UNNORMALIZED support indicators W_j, with
    # W^*W=diag(sizes). Correct norms are restored in the Gram matrix.
    derivatives = []
    compressions = []
    for H in hs:
        cols = []
        for support in supports:
            out = defaultdict(lambda: s.Integer(0))
            for word in support:
                for site in range(4):
                    for target in range(4):
                        value = H[target, word[site]]
                        if value:
                            w = list(word)
                            w[site] = target
                            out[tuple(w)] += value
            cols.append({w: s.expand(v) for w, v in out.items() if s.expand(v) != 0})
        compression = s.Matrix(5, 5, lambda a, b: sum(cols[b].get(w, 0) for w in supports[a]))
        assert compression == s.trace(H) * s.diag(*sizes)
        derivatives.append(cols)
        compressions.append(compression)

    def inner(left, right):
        return s.simplify(sum(
            s.Rational(1, sizes[j]) * sum(s.conjugate(value) * right[j].get(word, 0) for word, value in left[j].items())
            for j in range(5)
        ))

    gram = s.Matrix(15, 15, lambda a, b: inner(derivatives[a+1], derivatives[b+1]))
    assert gram == 32 * s.eye(15)
    # Scalar tangent D(I)=4I has zero leakage; all traceless tangents leak.
    assert derivatives[0] == [{word: s.Integer(4) for word in support} for support in supports]

    marker = build_marker_data()['data']
    rotation = _parse_matrix(marker['basis_alignment']['rotation_matrix'])
    yn = _parse_matrix(marker['native_polar_transport']['native_charge_matrix'])
    y = clean(rotation.T * yn * rotation)
    p2 = clean((6*y + 2*s.eye(5))/5)
    p3 = s.eye(5)-p2
    frame = _process_frame()

    def perm(p):
        mat = s.zeros(6)
        for col, row in enumerate(p):
            mat[row, col] = 1
        return clean(frame.T * mat * frame)

    sigma = perm((0, 1, 3, 4, 2, 5))
    tau = perm((1, 0, 2, 3, 4, 5))
    fsigma = clean((s.eye(5)+sigma+sigma**2)/3)
    R = clean(p3*(s.eye(5)-2*fsigma)+s.I*p2*tau)

    # Original-source N labels and the process C labels must be transported
    # together.  The earlier R here is a conjugate of the source-consistent R.
    reverse5 = perm((5, 4, 3, 2, 1, 0))
    assert clean(rotation+reverse5) == s.zeros(5)
    sigma4 = _source_operators()[0]
    mapping = [next(row for row in range(4) if sigma4[row, col] == 1) for col in range(4)]
    sigma_actual = s.Matrix(5, 5, lambda a, b: s.simplify(
        sum(tuple(mapping[x] for x in word) in supports[a] for word in supports[b])
        / s.sqrt(sizes[a]*sizes[b])
    ))
    assert clean(sigma_actual-perm((0, 3, 1, 2, 4, 5))) == s.zeros(5)
    assert clean(sigma_actual-rotation*sigma*rotation.T) == s.zeros(5)
    tau_source = clean(rotation*tau*rotation.T)
    assert tau_source == perm((0, 1, 2, 3, 5, 4))
    p2_source = clean((6*yn+2*s.eye(5))/5)
    f_source = clean((s.eye(5)+sigma_actual+sigma_actual**2)/3)
    R_source = clean((s.eye(5)-p2_source)*(s.eye(5)-2*f_source)+s.I*p2_source*tau_source)
    assert clean(R_source-rotation*R*rotation.T) == s.zeros(5)
    assert s.simplify(s.trace(R_source)) == 1
    assert s.simplify(s.trace(((R_source-R_source.conjugate())/(2*s.I)).T*((R_source-R_source.conjugate())/(2*s.I)))) == 2
    assert all(entry.is_real is True for entry in p2)
    assert clean(p2**2-p2) == s.zeros(5) and s.trace(p2) == 2
    assert clean(R.conjugate().T*R-s.eye(5)) == s.zeros(5)
    assert s.simplify(s.trace(R)) == 1
    assert clean((R-R.conjugate())/(2*s.I)-p2*tau) == s.zeros(5)
    imaginary_hs2 = s.simplify(s.trace((p2*tau).T*(p2*tau)))
    assert imaginary_hs2 == 2
    charpoly = s.factor(R.charpoly().as_expr())
    lam = s.Symbol('lambda')
    assert s.expand(charpoly-(lam-1)**2*(lam+1)*(lam**2+1)) == 0
    assert clean(R**4-s.eye(5)) == s.zeros(5)
    out = {
        'status': 'PASS_EXACT_NATIVE_TRANSVERSAL_TANGENT_AND_PHASE_AUDIT',
        'actual_constructor': str(REPO/'tfpt_explorer/process.py')+'::_code_basis',
        'code_support_sizes': sizes,
        'two_register_marginal_of_code_projector': '(I_16+Swap)/4',
        'all_register_pairs_checked_exactly': pair_marginals,
        'traceless_hermitian_basis': labels[1:],
        'exact_leakage_gram': '32 I_15',
        'exact_tangent_leakage_rank': 15,
        'scalar_tangent_leakage_rank': 0,
        'lie_algebra_dimension_real': 1,
        'native_marked_R_trace': '1',
        'native_marked_R_characteristic_polynomial': str(charpoly),
        'native_marked_R_imaginary_part_HS_squared': str(imaginary_hs2),
        'source_frame_correction': {
            'O_equals': '-rho_C(reverse)',
            'unchanged_sigma4_C_quotient': '(1 3 2)',
            'unchanged_sigma4_N_quotient': '(2 3 4)',
            'weak_N01_equals_C': '(4 5)',
            'source_consistent_R_equals': 'O R_previous_C O^T',
            'effect_on_stabilizer_and_projective_nonlift': 'none; O in U(1)rho5(G31)',
        },
        'projective_S6_nonmembership_proof': 'All S6 process matrices are real. Trace R=1 forces any unit phase relating R to one of them to be +/-1; Im(R) has HS squared 2.',
        'compact_stabilizer_conclusion': 'K=U(1)G31, conditional on prior actual Clifford identification and NRS Theorem 6.5',
        'general_product_stabilizer_conclusion': 'Every code-preserving U1 tensor U2 tensor U3 tensor U4 has Ui=phase_i U, by pair-marginal invariance; hence logical action is phase times rho5(G31).',
        'classical_native_mixture_conclusion': 'No convex mixture of native group-word channels equals Ad_R5 on all logical inputs, by rank-one Choi extremality and projective R5 nonmembership.',
        'scope': 'All product unitaries preserving this actual code, including unequal local factors, and convex mixtures of native group-word logical channels. Collective operations, ancillary protocols and general LOCC are not excluded.',
        'literature': 'https://arxiv.org/pdf/math/0001038v2#page=18',
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
