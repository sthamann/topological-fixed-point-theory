"""Exact controls for ordered family loops and selection of a pre-P1 process."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import runpy
import sys

import sympy as sp


root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root / 'verification'))
from v117_monodromy_weyl_a3 import M0_EXACT as M, U_EXACT as U, enumerate_group
from v94_sheet_diamond import R, Q, K, SIG, ONE, A as MARK


I = sp.eye(3)
II = sp.I
S = sp.diag(-1, 1, 1)
S[1, 1] = S[2, 2] = 0
S[1, 2] = S[2, 1] = 1


def clean(value):
    return value.applyfunc(sp.expand)


def deck(value):
    return clean(sum((U**k * value * U**(-k) for k in range(4)), sp.zeros(3)) / 4)


def charge(value, q):
    return clean(sum((II**(-q*k)*U**k*value*U**(-k) for k in range(4)), sp.zeros(3)) / 4)


def ordered_words():
    words = [I, M, U, M.H, U.H, M*U, M*U.H, M.H*U, M.H*U.H]
    flattened = sp.Matrix.hstack(*(w.T.reshape(9, 1) for w in words))
    determinant = sp.simplify(flattened.det())
    defect = clean(deck(M*M.H) - deck(M)*deck(M.H))
    assert abs(determinant) == 16 and defect == sp.diag(1, sp.Rational(3, 4), sp.Rational(3, 4))
    sectors = [charge(M, q) for q in range(4)]
    assert sum(sectors, sp.zeros(3)) == M
    assert clean(sum((piece*piece.H for piece in sectors), sp.zeros(3))) == I
    loop = clean((M*U-U*M).H*(M*U-U*M))
    assert deck(loop) == 2*I
    assert deck(M)*U == U*deck(M)
    assert clean((M*U).H*(U*M)).trace() == 0
    return {'word_basis_det_column_major': str(determinant), 'word_basis_det_row_major': str(-determinant),
            'early_deck_loss': list(map(str, defect.diagonal())),
            'charge_return_diagonals': [list(map(str, clean(part*part.H).diagonal())) for part in sectors],
            'closed_order_contrast': str(loop.trace()/3), 'closed_order_contrast_after_early_projection': '0'}


def cp_trace(a, b):
    return sp.expand(sp.trace(clean(a*b-b*a)**3))


def orientation():
    assert S*S.conjugate() == I
    assert clean(S*M.conjugate()*S) == M and clean(S*U.conjugate()*S) == U
    a = 3*I + M + M.H
    v = M*U*M.H
    b = 3*I + v + v.H
    k = clean(a+b)
    signed = clean(II*(a*b-b*a)/2)
    zero = clean(U*k*U.H)
    plus = clean(U*(k+signed)*U.H)
    minus = clean(U*(k-signed)*U.H)
    assert k.is_positive_definite and plus.is_positive_definite and minus.is_positive_definite
    assert all(sp.discriminant(x.charpoly().as_expr(), x.charpoly().gen) != 0 for x in (k, zero, plus, minus))
    assert clean(S*plus.conjugate()*S) == minus
    assert plus.charpoly().as_expr() == minus.charpoly().as_expr()
    traces = [cp_trace(k, item) for item in (zero, plus, minus)]
    assert traces == [0, -24*II, 24*II]
    assert clean(k*zero-zero*k) != sp.zeros(3)
    return {'common_real_structure_matrix': S.tolist(), 'example_CP_traces': list(map(str, traces)),
            'opposite_orientation_spectra_equal': True, 'example_positive': True}


def marked_frame():
    change = I.copy()
    change[0, 2] = -1
    first = I[:, 0]
    assert change.det() == 1 and change*first == first
    changed_k, changed_r = change.inv()*K*change, change.inv()*R*change
    changed_q, changed_sig = change.inv()*Q*change, change.inv()*SIG*change
    assert changed_k == changed_r + changed_q*changed_sig
    def pairs(rows):
        x, y = rows
        return [x[i]*y[j]-x[j]*y[i] for i, j in ((0, 1), (0, 2), (1, 2))]

    rows = [mark.T*K for mark in (ONE, MARK)]
    changed_rows = [mark.T*change*changed_k for mark in (ONE, MARK)]
    assert changed_rows[0][0] == rows[0][0] == 13
    assert pairs(rows) == [-1, 6, 4] and pairs(changed_rows) == [-1, 6, 3]
    assert pairs([row*change.inv() for row in changed_rows]) == pairs(rows)
    return {'basis_change': change.tolist(), 'det_R': int(changed_r.det()), 'det_K': int(changed_k.det()),
            'denominator_first': 13, 'coordinate_readouts': ['55/117', '50/117'],
            'with_mark_vectors_transported': '55/117'}


def channel_selection():
    group = enumerate_group([M, U])
    assert len(group) == 24
    def orbit(g):
        return set(tuple(clean(h*g*h.H)) for h in group)

    def action(g):
        return clean(sp.kronecker_product(g.conjugate(), g))

    def mean(elements):
        return clean(sum((action(g) for g in elements), sp.zeros(9)) / len(elements))

    haar = mean(group)
    identity_vector = I.reshape(9, 1)
    assert haar == identity_vector*identity_vector.T/3
    pairs = []
    for g in (M, U):
        class_elements = [sp.Matrix(3, 3, w) for w in sorted(orbit(g), key=str)]
        average = mean(class_elements)
        assert average == average.H
        channel = clean(sp.eye(9)/2 + haar/4 + average/4)
        assert all(clean(action(h)*channel-channel*action(h)) == sp.zeros(9) for h in (M, U))
        eigenvalues = channel.eigenvals()
        assert eigenvalues[1] == 1 and all(x >= sp.Rational(1, 4) for x in eigenvalues)
        marks = {name: word.T.reshape(9, 1) for name, word in [('M', M), ('U', U), ('MU', M*U)]}
        pairs.append({'class_size': len(class_elements), 'eigenvalues_with_multiplicity':
                      {str(x): int(y) for x, y in sorted(eigenvalues.items())},
                      'marked_returns': {name: str(sp.expand((v.H*channel*v)[0]/3))
                                         for name, v in marks.items()}})
    assert pairs[0] != pairs[1]
    assert pairs[0]['marked_returns'] != pairs[1]['marked_returns']
    # Coherent reversible realization of the very same group average.
    # Environment records are projected onto the nine matrix coefficient vectors.
    coefficient = sp.Matrix(24, 9, lambda g, ij: group[g][ij//3, ij%3] / sp.sqrt(8))
    assert clean(coefficient.H*coefficient) == sp.eye(9)

    # Track the actual Weyl sign through the two native generators.
    signs = {tuple(I): 1}
    frontier = [I]
    while frontier:
        current = frontier.pop()
        for generator, parity in ((M, 1), (U, -1)):
            value = clean(generator * current)
            key, sign = tuple(value), signs[tuple(current)]*parity
            if key in signs:
                assert signs[key] == sign
            else:
                signs[key] = sign
                frontier.append(value)
    assert len(signs) == 24 and sum(signs.values()) == 0
    signed_part = clean(sum((signs[tuple(g)]*action(g) for g in group), sp.zeros(9))/24)
    assert signed_part == sp.zeros(9)
    # Two positive distributions for the same channel, distinct if the
    # recorded Weyl sign is retained. No channel-only observable sees it.
    bias = sp.Rational(1, 2)
    recorded_sign = [sum(s * (1 + epsilon*s)/24 for s in signs.values()) for epsilon in (bias, -bias)]
    assert recorded_sign == [bias, -bias]
    return {'group_size': 24, 'haar_channel': 'Tr(rho) I/3',
            'two_covariant_positive_channels': pairs,
            'coherent_environment_isometry_rank': 9, 'complement_in_matrix_coefficient_frame': 'I/3 tensor rho',
            'same_reduced_channel_two_records': {'sign_expectations': list(map(str, recorded_sign)),
                                                 'signed_adjoint_group_average': 'zero'}}


def affine_order():
    path = (root / 'experiments/theory-contracts/recursive-mirrored-source-20261103'
            / 'check_seam_yukawa_ordered_norms_20261111.py')
    with redirect_stdout(io.StringIO()):
        source = runpy.run_path(str(path))
    point = {**source['replay'], **source['one']}
    first, second, overlap = (source[key].subs(point)[:3, :3] for key in ('G_AB', 'G_BA', 'G_cross'))
    assert first*second-second*first == sp.zeros(3)
    loop = clean(first + second - overlap - overlap.H)
    assert loop.is_positive_definite
    assert clean(S*loop.conjugate()*S) == loop
    visibility = sp.cancel(loop.trace() / (first + second).trace())
    assert visibility > 0
    return {'original_two_orders_return_commutator': 'zero',
            'closed_order_difference_gram_diagonal': list(map(str, loop.diagonal())),
            'normalized_closed_order_difference': str(visibility),
            'normalized_closed_order_difference_numeric': float(visibility),
            'closed_order_difference_invariant_under_J': True}


print(json.dumps({'ordered_words': ordered_words(), 'CP': orientation(), 'marked_frame': marked_frame(),
                  'channel_selection': channel_selection(), 'actual_affine_order': affine_order()},
                 default=str, indent=2))
