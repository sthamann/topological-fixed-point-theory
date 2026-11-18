"""Certify an entire response segment between stored exact finite processes.

Only succeeds when ALL supplied moment polynomials stay constant, both full
representations pass the existing independent checker, and state, fields and
quarter map agree. This is not physical P1 source selection.
"""
import argparse
from functools import lru_cache
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp

from proof_check import check_pair, _matrix


def certify(run):
    run = Path(run)
    payload = json.loads((run/'moments.json').read_text())
    low = json.loads((run/'min-trace.json').read_text())
    high = json.loads((run/'max-trace.json').read_text())
    question = json.loads((run/'next-source-question.json').read_text())
    checked = check_pair(low, high, payload, question)
    if not checked.get('ok'):
        raise ValueError('Stored endpoints failed existing exact full-process checks: '+str(checked))
    if sp.Rational(str(payload.get('cluster_gap', 0))) <= 0:
        raise ValueError('This segment certificate requires an explicitly positive common clustering gap')
    for key in ('Omega', 'gammas', 'quarter'):
        if low['exact_matrices'][key] != high['exact_matrices'][key]:
            raise ValueError('Segment requires the same original state/fields/quarter: '+key)
    for key in ('m', 'dimension', 'inferred_occupations', 'Omega_basis_index'):
        if low[key] != high[key]:
            raise ValueError('Endpoint representation mismatch: '+key)
    n, s = low['dimension'], sp.Symbol('s', real=True)
    a = _matrix(low['exact_matrices']['T'], n, n, 'Tmin')
    b = _matrix(high['exact_matrices']['T'], n, n, 'Tmax')
    omega = _matrix(low['exact_matrices']['Omega'], n, 1, 'Omega')
    operators = {-1: a+s*(b-a)}
    operators.update({j: _matrix(g, n, n, 'gamma') for j, g in enumerate(low['exact_matrices']['gammas'])})
    annihilates_vacuum = (b-a)*omega == sp.zeros(n, 1)
    annihilates_single_fields = [(b-a)*operators[j]*omega == sp.zeros(n, 1) for j in range(2*low['m'])]
    columns = {}
    for token, operator in operators.items():
        entries = [[] for _ in range(n)]
        for (row, col), coefficient in operator.todok().items():
            entries[col].append((row, coefficient))
        columns[token] = entries
    @lru_cache(None)
    def vector(word):
        if not word:
            return tuple((row, c) for (row, col), c in omega.todok().items())
        answer = {}
        for col, value in vector(word[1:]):
            for row, coefficient in columns[word[0]][col]:
                answer[row] = answer.get(row, 0)+coefficient*value
        return tuple((row, sp.expand(value)) for row, value in sorted(answer.items()) if sp.expand(value) != 0)
    def moment(word):
        values = dict(vector(tuple(word)))
        return sp.expand(sum(sp.conjugate(c)*values.get(row, 0) for (row, col), c in omega.todok().items()))
    def recover(value):
        q = Fraction(value).limit_denominator(64)
        if abs(float(q)-value) > 1e-6:
            raise ValueError('Moment rational recovery failed')
        return sp.Rational(q.numerator, q.denominator)
    for row in payload['moments']:
        expected = recover(row['re'])+sp.I*recover(row['im'])
        if sp.expand(moment(row['word'])-expected) != 0:
            raise ValueError('Intermediate processes change supplied moment: '+str(row['word']))
    responses = []
    for row in question['differences']:
        polynomial = moment(row['word'])
        derivative = sp.diff(polynomial, s)
        if derivative.has(s) or derivative == 0:
            raise ValueError('This response does not have a nonconstant affine segment')
        ends = [polynomial.subs(s, x) for x in (0, 1)]
        if not all(x.is_real for x in ends):
            raise ValueError('Response endpoints must be real')
        responses.append({'word': row['word'], 'polynomial': str(polynomial),
                          'attained_interval_on_this_segment': [str(min(ends)), str(max(ends))]})
    # Reference-vector relation, quarter commutation and both positivity/gap
    # Loewner inequalities are affine in T. Equal original endpoints satisfy
    # them, hence every convex combination does. No spectrum interpolation
    # identity is claimed. All nonlinear input moments were checked above.
    return {'schema': 'tfpt.exact-full-response-segment.v1', 'status': 'EXACT_FULL_PROCESS_SEGMENT',
            'parameter': '0 <= s <= 1', 'transfer': 'T(s)=(1-s)Tmin+sTmax',
            'dimension': n, 'all_replayed_moment_polynomials_constant': len(payload['moments']),
            'structural_annihilator': {'difference_annihilates_reference': annihilates_vacuum,
                                      'difference_annihilates_each_single_field': annihilates_single_fields},
            'retained_constraints': ['same full CAR representation and normalized vector',
                                     'same internal quarter action', 'selfadjoint positive bounded transfer',
                                     'same reference eigenvector', 'same declared positive clustering gap'],
            'responses': responses, 'endpoint_check': checked,
            'source_sha256': {name: hashlib.sha256((run/name).read_bytes()).hexdigest()
                              for name in ('moments.json', 'min-trace.json', 'max-trace.json', 'next-source-question.json')},
            'scope': 'Constructed full finite segment, not all admissible processes. Attained intervals are not global extrema. Original input moments were rationally recovered under the existing checker contract; they are not P1-derived data.',
            'physical_P1_source_identified': False, 'absolute_Z0_recovered': False,
            'fundamental_solution': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    result = certify(args.run)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps({k: result[k] for k in ('status', 'dimension', 'all_replayed_moment_polynomials_constant', 'responses', 'fundamental_solution')}, indent=2))


if __name__ == '__main__':
    main()
