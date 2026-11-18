"""Bounded exact search for one vs. two neutral source directions.

This imports and calls the TFPT source implementations; it does not alter them.
"""
from __future__ import annotations
import argparse
import os
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

REPO = Path(os.environ.get('TFPT_REPO', '/Users/stefanhamann/Projekte/tfpt-theoryv4'))
sys.path.insert(0, str(REPO))

from tfpt_explorer.source_charge_response import exact_charge_response
from tfpt_explorer.source_majorana_pairs import _native_pairs, _joint_neutral_pair




def s(x):
    return str(sp.factor(x))


def run():
    src = exact_charge_response()  # computes from the original 48 retained E8 roots
    R = src['R']
    RY = src['R_Y']
    SY = sp.simplify(RY - sp.Rational(5, 24) * R)
    y = src['hypercharge']
    gram_R_SY = sp.Matrix([
        [sp.trace(R*R)/2, sp.trace(R*SY)/2],
        [sp.trace(SY*R)/2, sp.trace(SY*SY)/2],
    ])
    gram_R_RY = sp.Matrix([
        [sp.trace(R*R)/2, sp.trace(R*RY)/2],
        [sp.trace(R*RY)/2, sp.trace(RY*RY)/2],
    ])
    family_coeff = sp.Rational(5, 24)
    family_SY = sp.simplify(sp.trace(R*SY)/sp.trace(R*R))
    ward_R = sp.simplify((y.T*R*y)[0])
    ward_RY = sp.simplify((y.T*RY*y)[0])
    ward_SY = sp.simplify((y.T*SY*y)[0])

    # General real coefficients a,b for one operator F=a R+b RY.
    a, b = sp.symbols('a b', real=True)
    family_F = sp.simplify(a + family_coeff*b)
    norm_F = sp.factor(sp.trace((a*R+b*RY)**2)/2)
    norm_in_family_frame = sp.factor(
        48*family_F**2 + sp.Rational(115, 12)*b**2
    )
    norm_identity = sp.simplify(norm_F - norm_in_family_frame)
    # Enforce family coefficient one, then norm 48; positive term forces b=0.
    norm_under_family_one = sp.factor(norm_in_family_frame.subs(a, 1-family_coeff*b))
    norm_excess = sp.factor(norm_under_family_one - 48)
    ward_F = sp.factor(a*ward_R + b*ward_RY)
    ward_under_family_one = sp.factor(ward_F.subs(a, 1-family_coeff*b))

    # Gamma-pair responses come directly from source-owned pair states and
    # original R/R_Y matrices, not from a result JSON/metadata.
    native_pairs = _native_pairs()
    pair_R = _joint_neutral_pair(native_pairs)
    response_R_pairs = {
        tuple(row['family_pair']): sp.Rational(row['R0_eigenvalue'])
        for row in pair_R['symmetric_retained_six']
    }
    response_RY_pairs = {
        tuple(row['families']): sp.Rational(row['total'])
        for row in src['pair_responses']
    }
    response_SY_pairs = {
        pair: sp.simplify(response_RY_pairs[pair] - family_coeff*response_R_pairs[pair])
        for pair in response_R_pairs
    }

    # Enumerate the explicitly bounded rational coefficient grid; exact checks
    # verify that no one-direction candidate passes all three hard conditions.
    denominators = (1, 2, 3, 4, 6, 12, 24)
    coeffs = sorted({sp.Rational(n, d) for d in denominators for n in range(-12, 13)})
    candidates = []
    for aa in coeffs:
        for bb in coeffs:
            fcomp = sp.simplify(aa + family_coeff*bb)
            nrm = sp.simplify(norm_in_family_frame.subs({a: aa, b: bb}))
            w = sp.simplify(ward_F.subs({a: aa, b: bb}))
            if fcomp == 1 and nrm == 48 and w != 0:
                candidates.append({'a': s(aa), 'b': s(bb), 'norm': s(nrm), 'ward': s(w)})

    checks = {
        'same_original_48_root_set_used': len(src['charges']) == 48,
        'RY_equals_5over24_R_plus_SY': sp.simplify(RY-sp.Rational(5,24)*R-SY) == sp.zeros(8),
        'R_SY_orthogonal': gram_R_SY[0,1] == 0,
        'exact_gram_R_SY': gram_R_SY == sp.diag(48, sp.Rational(115,12)),
        'exact_gram_R_RY': gram_R_RY == sp.Matrix([[48,10],[10,sp.Rational(35,3)]]),
        'single_vertex_norm_identity': norm_identity == 0,
        'single_vertex_constraints_force_b_zero': norm_excess == sp.Rational(115,12)*b**2,
        'single_vertex_ward_zero_under_constraints': ward_under_family_one.subs(b, 0) == 0,
        'bounded_grid_has_no_counterexample': len(candidates) == 0,
        'source_R_and_RY_pair_responses': all(
            response_R_pairs[pair] == 4 and response_RY_pairs[pair] == sp.Rational(-5,3)
            for pair in response_R_pairs
        ),
        'SY_pair_response_minus_5over2': all(
            value == sp.Rational(-5,2) for value in response_SY_pairs.values()
        ),
    }
    assert all(checks.values()), checks

    result = {
        'status': 'PASS_SCOPED_SINGLE_DIRECTION_EXCLUSION_AND_TWO_DIRECTION_SOURCE_CANDIDATE',
        'scope': 'Only one neutral Grade-2 vertex F=aR+bRY with fixed family component 1 and norm 48; not an all-source or all-TFPT no-go.',
        'source_calls': [
            'tfpt_explorer.source_charge_response.exact_charge_response()',
            'tfpt_explorer.source_majorana_pairs._native_pairs()',
            'tfpt_explorer.source_majorana_pairs._joint_neutral_pair(native_pairs)',
        ],
        'native_48_root_outputs': {
            'count': len(src['charges']),
            'relative_coefficient': s(src['relative_coefficient']),
            'Gram_R_RY': [[s(v) for v in row] for row in gram_R_RY.tolist()],
            'norm_R': s(gram_R_SY[0,0]),
            'norm_RY': s(gram_R_RY[1,1]),
            'inner_R_RY': s(gram_R_RY[0,1]),
            'SY_definition': 'SY=RY-(5/24)R',
            'Gram_R_SY': [[s(v) for v in row] for row in gram_R_SY.tolist()],
            'family_response_SY_relative_to_R': s(family_SY),
            'R_JY_JY': s(ward_R),
            'RY_JY_JY': s(ward_RY),
            'SY_JY_JY': s(ward_SY),
            'chiral_current_level': s(src['current_level']),
        },
        'single_vertex_F_aR_bRY': {
            'family_component': s(family_F),
            'norm_squared_expanded': s(norm_F),
            'norm_squared_orthogonal_frame': s(norm_in_family_frame),
            'with_family_component_1_norm_squared': s(norm_under_family_one),
            'norm_excess_over_48': s(norm_excess),
            'chiral_hypercharge_ward': s(ward_F),
            'ward_with_family_component_1': s(ward_under_family_one),
            'bounded_rational_grid': {
                'numerators': [-12,12], 'denominators': list(denominators),
                'unique_coefficients_each': len(coeffs), 'tested_pairs': len(coeffs)**2,
                'passing_candidates': candidates,
            },
            'global_real_coefficients_proof': 'For real a,b, family=1 and norm²=48 imply 48+115 b²/12=48, hence b=0; then the J_Y,J_Y Ward is 0. Thus no real pair passes all conditions.',
        },
        'two_direction_candidate': {
            'directions': ['R', 'SY'],
            'Gram': [[s(v) for v in row] for row in gram_R_SY.tolist()],
            'determinant': s(gram_R_SY.det()),
            'ward_vector': [s(ward_R), s(ward_SY)],
            'RY_reconstruction': 'RY=(5/24)R+SY',
            'coefficient_5over24': s(family_coeff),
            'source_derived_not_fit': sp.simplify(RY-sp.Rational(5,24)*R-SY) == sp.zeros(8),
        },
        'Gamma_pair_responses_direct_source': {
            'family_pairs': [[list(pair), s(value)] for pair,value in response_R_pairs.items()],
            'R': sorted({s(v) for v in response_R_pairs.values()}),
            'RY': sorted({s(v) for v in response_RY_pairs.values()}),
            'SY': sorted({s(v) for v in response_SY_pairs.values()}),
        },
        'checks': checks,
    }
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outputpath', type=Path)
    args = parser.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, ensure_ascii=False) + '\n'
    if args.outputpath:
        args.outputpath.parent.mkdir(parents=True, exist_ok=True)
        args.outputpath.write_text(payload)
    else:
        print(payload, end='')
