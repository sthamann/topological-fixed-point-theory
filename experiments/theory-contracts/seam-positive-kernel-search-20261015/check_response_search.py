"""Decisive exact full-representation controls, not new TFPT source inputs."""
from copy import deepcopy
import json
from pathlib import Path

import sympy as sp

from kernel_search import MomentSearch
from response_search import (evaluate, normalize_response, affine_certificate,
                             verify_affine_certificate)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def run():
    # Complete 2D CAR representation. Every word exists, with no truncation.
    t = sp.Symbol('t', real=True)
    g0, g1 = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]])
    transfer, omega = sp.diag(1, t), sp.Matrix([1, 0])
    quarter = sp.diag(1, sp.I)
    require(g0*g1 == -g1*g0 and g0*g0 == sp.eye(2) and g1*g1 == sp.eye(2), 'CAR')
    require(quarter*g0*quarter.H == g1 and quarter*g1*quarter.H == -g0, 'Quarter')
    require(quarter*omega == omega and quarter*transfer*quarter.H == transfer, 'Whole-state invariance')
    require(transfer*omega == omega, 'Reference eigenvector')
    a = (omega.H*g0*transfer*g0*omega)[0]
    b = (omega.H*g1*transfer*g1*omega)[0]
    require(a == t and b == t, 'Exact full witness responses')
    # 0<=T<=I iff 0<=t<=1. Hence these exact full witnesses attain both
    # global bounds for g1*T*g1; the finite SDP endpoints alone do not prove it.
    require(transfer.eigenvals() == {sp.Integer(1): 1, t: 1}, 'Exact positive-transfer interval')

    query = {'terms': [{'word': [1, -1, 1], 're': '1'}]}
    unpinned = MomentSearch(m=1, degree=2)
    free = evaluate(unpinned, query)
    require(not free['affine_certificate']['identity_proved'], 'Free response must remain free')
    require(len(free['numerical_endpoints']) == 2 and all('estimate' in x for x in free['numerical_endpoints']), 'Both searches must actually succeed')
    require(abs(free['numerical_endpoints'][0]['estimate']) < 1e-6 and abs(free['numerical_endpoints'][1]['estimate']-1) < 1e-6, 'Numerical range versus exact full witnesses')

    fixed = [{'word': [0, -1, 0], 're': '1/3'}]
    linked = MomentSearch(m=1, degree=2, fixed_moments=fixed)
    constrained = evaluate(linked, query)
    proof = constrained['affine_certificate']
    require(proof.get('value') == '1/3', 'Held-out response must follow from a different supplied response')
    polynomial = normalize_response(linked, query)
    require(verify_affine_certificate(linked, polynomial, proof), 'Independent exact replay')
    require(a.subs(t, sp.Rational(1, 3)) == b.subs(t, sp.Rational(1, 3)) == sp.Rational(1, 3), 'Exact nonempty full witness for linked control')

    altered = deepcopy(proof); altered['value'] = '2/3'
    try:
        verify_affine_certificate(linked, polynomial, altered)
    except ValueError:
        pass
    else:
        raise RuntimeError('A changed certificate value must fail')

    contradictory = MomentSearch(m=1, degree=2, fixed_moments=fixed + [{'word': [1, -1, 1], 're': '2/3'}])
    bad = evaluate(contradictory, query)
    require(bad['affine_certificate']['status'] == 'INCONSISTENT_AFFINE_INPUT' and not bad['numerical_endpoints'], 'No vacuous uniqueness from contradictory sources')

    complex_search = MomentSearch(m=1, degree=2, fixed_moments=[{'word': [0, 1], 're': '0', 'im': '1/2'}])
    complex_query = {'terms': [{'word': [0, 1], 're': '0', 'im': '1'}]}
    complex_proof = affine_certificate(complex_search, normalize_response(complex_search, complex_query))
    require(complex_proof.get('value') == '-1/2', 'Imaginary moment coefficient and adjoint sign')
    empty_positive = MomentSearch(m=1, degree=2, fixed_moments=[{'word': [0, 1], 're': '0', 'im': '2'}])
    empty_result = evaluate(empty_positive, complex_query)
    require(empty_result['affine_certificate'].get('value') == '-2' and not empty_result['positivity_feasibility_proved'], 'Affine identity must not imply existence')
    require(all(e['solver_status'] == 'infeasible' for e in empty_result['numerical_endpoints']), 'Positivity contradiction must remain visible')
    rounded = MomentSearch(m=1, degree=2, fixed_moments=[{'word': [], 're': '9007199254740993/9007199254740992'}])
    require(affine_certificate(rounded, normalize_response(rounded, query))['status'] == 'INCONSISTENT_AFFINE_INPUT', 'Rational contradiction must survive float rounding')
    try:
        normalize_response(unpinned, {'terms': [{'word': [0, 1], 're': '1'}]})
    except ValueError:
        pass
    else:
        raise RuntimeError('Anti-Hermitian response was accepted')
    try:
        normalize_response(unpinned, {**query, 'target': '1/3'})
    except ValueError:
        pass
    else:
        raise RuntimeError('A query target must not become a fitted constraint')

    result = {'status': 'PASS', 'role': 'software_and_exact_full_control_only',
              'exact_full_unpinned_range': ['0', '1'], 'exact_full_linked_response': '1/3',
              'unconstrained_search': free, 'linked_search': constrained,
              'contradictory_input': bad, 'complex_coefficient_certificate': complex_proof,
              'positive_empty_input': empty_result, 'float_rounding_contradiction_detected': True,
              'fundamental_solution': False, 'physical_P1_source_input': False}
    output = Path(__file__).with_name('runs')/'response-search-20261102'
    output.mkdir(parents=True, exist_ok=True)
    (output/'controls.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'role', 'exact_full_unpinned_range', 'exact_full_linked_response', 'fundamental_solution')}, indent=2))
    return result


if __name__ == '__main__':
    run()
