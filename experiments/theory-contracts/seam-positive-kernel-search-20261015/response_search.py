"""Withhold a Hermitian response, then search its range in the EXISTING class.

Exact affine certificates are conditional identities, not existence proofs.
SDP endpoints are numerical estimates, never certified outer bounds.
Nothing in this module identifies the declared CAR interface with physical P1.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from kernel_search import MomentSearch, load_source_control
from linear_relations import _q
from word_algebra import T


def add(row, key, value):
    row[key] = row.get(key, Q(0)) + value
    if not row[key]:
        del row[key]


def normalize_response(search, response):
    if not isinstance(response, dict) or set(response) != {'terms'}:
        raise ValueError('Response requires terms only; no fitted target value')
    if not isinstance(response['terms'], list) or not response['terms']:
        raise ValueError('Response terms must be a nonempty list')
    polynomial = {}
    for term in response['terms']:
        if not isinstance(term, dict) or set(term) - {'word', 're', 'im'} or not {'word', 're'} <= set(term):
            raise ValueError('Response terms require word/re and optional im')
        sign, word = search.algebra.reduce_word(term['word'])
        if word not in search.index:
            raise ValueError('Response exceeds this relaxation; explicitly increase degree')
        re, im = polynomial.get(word, (Q(0), Q(0)))
        polynomial[word] = (re + sign*_q(term['re']), im + sign*_q(term.get('im', 0)))
    polynomial = {w: c for w, c in polynomial.items() if c != (0, 0)}
    adjoint = {}
    for word, (re, im) in polynomial.items():
        sign, star = search.algebra.star(word)
        adjoint[star] = (sign*re, -sign*im)
    if polynomial != adjoint:
        raise ValueError('Response is not Hermitian in the exact original word algebra')
    return polynomial


def affine_rows(search):
    """Replay the ideal rational contract of MomentSearch, before float rounding."""
    rows = []
    def equation(label, terms, rhs=(Q(0), Q(0))):
        real, imag = {}, {}
        for word, re, im, conjugate in terms:
            i = 2*search.index[word]
            add(real, i, re); add(real, i+1, im if conjugate else -im)
            add(imag, i, im); add(imag, i+1, -re if conjugate else re)
        rows.extend([(label+':re', real, Q(rhs[0])), (label+':im', imag, Q(rhs[1]))])
    equation('normalized_conditional_state', [((), Q(1), Q(0), False)], (Q(1), Q(0)))
    lam = Q(str(search.reference_eigenvalue))
    for word in search.words:
        sign, star = search.algebra.star(word)
        equation('adjoint:'+str(word), [(star, Q(sign), Q(0), False), (word, Q(-1), Q(0), True)])
        sign, rotated = search.algebra.quarter(word)
        equation('declared_internal_quarter:'+str(word), [(rotated, Q(sign), Q(0), False), (word, Q(-1), Q(0), False)])
        for side, remainder in [('left', word[1:]), ('right', word[:-1])]:
            if word and (word[0] if side == 'left' else word[-1]) == T:
                equation('reference_'+side+':'+str(word), [(word, Q(1), Q(0), False), (remainder, -lam, Q(0), False)])
    for j, row in enumerate(search.fixed_moments):
        equation('supplied_fixed:'+str(j), [(tuple(row['word']), Q(1), Q(0), False)], (Q(row['re']), Q(row['im'])))
    for j, row in enumerate(search.linear_relations):
        equation('supplied_relation:'+str(j), [(tuple(t['word']), Q(t['re']), Q(t['im']), False) for t in row['terms']], (Q(row['re']), Q(row['im'])))
    return rows


def target_row(search, polynomial):
    row = {}
    for word, (re, im) in polynomial.items():
        add(row, 2*search.index[word], re)
        add(row, 2*search.index[word]+1, -im)
    return row


def subtract(row, other, factor):
    for key, value in other.items():
        add(row, key, -factor*value)


def affine_certificate(search, polynomial):
    rows, basis = affine_rows(search), {}
    for index, (_, initial, initial_rhs) in enumerate(rows):
        row, rhs, combination = dict(initial), initial_rhs, {index: Q(1)}
        while row and min(row) in basis:
            pivot = min(row); factor = row[pivot]
            old, old_rhs, old_combination = basis[pivot]
            subtract(row, old, factor); rhs -= factor*old_rhs
            subtract(combination, old_combination, factor)
        if not row:
            if rhs:
                return {'status': 'INCONSISTENT_AFFINE_INPUT', 'identity_proved': False}
            continue
        pivot = min(row); factor = row[pivot]
        basis[pivot] = ({k: v/factor for k, v in row.items()}, rhs/factor,
                        {k: v/factor for k, v in combination.items()})
    residual, value, combination = target_row(search, polynomial), Q(0), {}
    # Reduce every pivot; a free column must not hide later pivot columns.
    for pivot in sorted(basis):
        factor = residual.get(pivot, Q(0))
        if not factor:
            continue
        row, rhs, weights = basis[pivot]
        subtract(residual, row, factor); value += factor*rhs
        for key, weight in weights.items():
            add(combination, key, factor*weight)
    if residual:
        return {'status': 'NOT_FORCED_BY_AFFINE_EQUALITIES', 'identity_proved': False,
                'remaining_real_coordinate_count': len(residual),
                'note': 'Positivity may impose further restrictions; this is not a nonuniqueness proof.'}
    certificate = {'status': 'EXACT_IDENTITY_CONDITIONAL_ON_EXISTENCE', 'identity_proved': True,
                   'value': str(value),
                   'weights': [{'row': j, 'label': rows[j][0], 'weight': str(weight)} for j, weight in sorted(combination.items())],
                   'existence_proved': False}
    verify_affine_certificate(search, polynomial, certificate)
    return certificate


def verify_affine_certificate(search, polynomial, certificate):
    """Independent replay of a claimed linear-combination certificate."""
    if certificate.get('status') != 'EXACT_IDENTITY_CONDITIONAL_ON_EXISTENCE':
        raise ValueError('Not an identity certificate')
    rows, combined, value = affine_rows(search), {}, Q(0)
    for term in certificate['weights']:
        label, row, rhs = rows[term['row']]
        if label != term['label']:
            raise ValueError('Certificate row label mismatch')
        weight = Q(term['weight'])
        for key, coefficient in row.items():
            add(combined, key, weight*coefficient)
        value += weight*rhs
    if combined != target_row(search, polynomial) or value != Q(certificate['value']):
        raise ValueError('Exact response certificate replay failed')
    return True


def evaluate(search, response, solver='CLARABEL', max_iters=20000):
    polynomial = normalize_response(search, response)
    certificate = affine_certificate(search, polynomial)
    expression = cp.real(sum((complex(float(re), float(im))*search.moment(word)
                              for word, (re, im) in polynomial.items()), 0*search.y[0]))
    endpoints = []
    if certificate['status'] != 'INCONSISTENT_AFFINE_INPUT':
        for sense in ('min', 'max'):
            objective = cp.Minimize(expression) if sense == 'min' else cp.Maximize(expression)
            problem = cp.Problem(objective, search.constraints)
            options = {'eps': 1e-8, 'max_iters': max_iters} if solver == 'SCS' else {
                'tol_gap_abs': 1e-8, 'tol_feas': 1e-8, 'tol_gap_rel': 1e-8, 'max_iter': max_iters}
            endpoint = {'sense': sense}
            try:
                problem.solve(solver=solver, **options)
                endpoint['solver_status'] = problem.status
                if search.y.value is not None and problem.status in ('optimal', 'optimal_inaccurate'):
                    endpoint['estimate'] = float(expression.value)
                    endpoint['max_constraint_violation'] = max(float(np.max(c.violation())) for c in search.constraints)
            except cp.error.SolverError as error:
                endpoint.update(solver_status='solver_error', error=str(error))
            endpoints.append(endpoint)
    return {'schema': 'tfpt.withheld-response-search.v1', 'response': response,
            'affine_certificate': certificate, 'numerical_endpoints': endpoints,
            'numerical_endpoints_are_certified_bounds': False,
            'numerical_endpoints_are_exact_realizations': False,
            'positivity_feasibility_proved': False,
            'coefficient_semantics': {'certificate': 'Exact rational input contract before float conversion',
                                      'numerical_search': 'Binary floating approximation of that contract; tiny exact contradictions can disappear under rounding'},
            'source_control': {'fixed_moments': search.fixed_moments, 'linear_relations': search.linear_relations,
                               'source_observables': search.source_observables},
            'class': {'m': search.m, 'degree': search.degree, 'reference_eigenvalue': search.reference_eigenvalue,
                      'transfer_floor': search.transfer_floor, 'transfer_ceiling': search.transfer_ceiling,
                      'transfer_kind': search.transfer_kind},
            'scope': 'Existing truncated CAR/internal-quarter/reference-transfer class. A conditional affine identity applies to every full realization of these equalities. Numerical separation needs exact full witnesses before a nonuniqueness claim.',
            'physical_source_identified': False, 'absolute_Z0_recovered': False,
            'fundamental_solution': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--response', type=Path)
    parser.add_argument('--verify', type=Path, help='Replay a stored exact affine certificate without any optimizer call')
    parser.add_argument('--source-control', type=Path)
    parser.add_argument('--m', type=int, default=5)
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--transfer-floor', type=float, default=0)
    parser.add_argument('--transfer-ceiling', type=float, default=1)
    parser.add_argument('--reference-eigenvalue', type=float, default=1)
    parser.add_argument('--transfer-kind', choices=('vacuum_effect', 'thermal_gns', 'bounded_positive'), default='vacuum_effect')
    parser.add_argument('--solver', choices=('CLARABEL', 'SCS'), default='CLARABEL')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.verify:
        result = json.loads(args.verify.read_text())
        if result.get('schema') != 'tfpt.withheld-response-search.v1':
            parser.error('Unexpected response-result schema')
        config = result['class']
        search = MomentSearch(config['m'], config['degree'], config['transfer_floor'],
                              transfer_ceiling=config['transfer_ceiling'], transfer_kind=config['transfer_kind'],
                              reference_eigenvalue=config['reference_eigenvalue'], **result['source_control'])
        verify_affine_certificate(search, normalize_response(search, result['response']), result['affine_certificate'])
        print(json.dumps({'status': 'EXACT_CONDITIONAL_IDENTITY_REPLAYED',
                          'value': result['affine_certificate']['value'],
                          'existence_proved': False, 'physical_source_identified': False}))
        return
    if not args.response or not args.out:
        parser.error('--response and --out are required when not replaying a certificate')
    search = MomentSearch(args.m, args.degree, args.transfer_floor,
                          transfer_ceiling=args.transfer_ceiling, transfer_kind=args.transfer_kind,
                          reference_eigenvalue=args.reference_eigenvalue, **load_source_control(args.source_control))
    response_bytes = args.response.read_bytes()
    result = evaluate(search, json.loads(response_bytes), args.solver)
    result['response_sha256'] = hashlib.sha256(response_bytes).hexdigest()
    result['code_sha256'] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in ('response_search.py', 'kernel_search.py', 'word_algebra.py', 'linear_relations.py')}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation preserves previous runs, including late races.
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_control', 'code_sha256')}, indent=2))


if __name__ == '__main__':
    main()
