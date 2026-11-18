"""Finite faithful CAR GNS recovery in the supplied state, with exact replay.

The metric is inferred from all supplied CAR moments. No pure state, Hamiltonian,
thermal formula or Wick rule is substituted. An unsuccessful recovery is
inconclusive. This is a finite source interface, not a P1 selection theorem.
"""
from fractions import Fraction
from itertools import combinations
from hashlib import sha256
import json

import cvxpy as cp
import numpy as np
import sympy as sp

from word_algebra import T, WordAlgebra
from linear_relations import check_linear_relations


def _number(value):
    if isinstance(value, bool):
        raise ValueError('Boolean is not an exact moment component')
    if isinstance(value, str):
        q = Fraction(value)
    else:
        q = Fraction(str(value)).limit_denominator(1000000)
        if not np.isfinite(float(value)) or abs(float(q)-float(value)) > 1e-8:
            raise ValueError('Moment has no supported exact recovery')
    return sp.Rational(q.numerator, q.denominator)


def _json_matrix(matrix):
    return [[{'re': str(sp.re(x)), 'im': str(sp.im(x))} for x in row]
            for row in matrix.tolist()]


def _table(payload, algebra):
    values = {}
    fixed = {}
    for row in payload.get('fixed_moments', []):
        phase, word = algebra.reduce_word(row['word'])
        target = sp.expand(phase*(_number(row['re'])+sp.I*_number(row.get('im', 0))))
        if word in fixed and fixed[word] != target:
            raise ValueError('Conflicting exact source controls')
        fixed[word] = target
    for row in payload['moments']:
        phase, word = algebra.reduce_word(row['word'])
        exact_input = isinstance(row['re'], str) and isinstance(row.get('im', '0'), str)
        if word in fixed and not exact_input:
            numeric = phase*complex(row['re'], row.get('im', 0))
            if not np.isfinite(numeric) or abs(numeric-complex(fixed[word])) > 1e-6:
                raise ValueError('Numerical candidate disagrees with exact source control')
            value = fixed[word]
        else:
            value = sp.expand(phase*(_number(row['re'])+sp.I*_number(row.get('im', 0))))
            if word in fixed and value != fixed[word]:
                raise ValueError('Exact candidate disagrees with exact source control')
        if word in values and values[word] != value:
            raise ValueError('Conflicting exact moments')
        values[word] = value
    # Direct completion may receive the newly queried word before a larger
    # moment prefix is generated. Preserve and replay that source answer too.
    for word, value in fixed.items():
        values.setdefault(word, value)
    return values


def first_missing_probe(payload):
    """Find a new Hermitian one-transfer query from the finite CAR basis."""
    m = payload['m']
    algebra = WordAlgebra(m)
    basis = [tuple(c) for k in range(2*m+1) for c in combinations(range(2*m), k)]
    available = {algebra.reduce_word(row['word'])[1]
                 for row in payload['moments']+payload.get('fixed_moments', [])}
    queries = set()
    for left in basis:
        for right in basis:
            _, word = algebra.reduce_word(left+(T,)+right)
            phase, adjoint = algebra.star(word)
            if word not in available and phase == 1 and adjoint == word:
                queries.add(word)
    return min(queries, key=lambda word: (len(word), word)) if queries else None


def recover(payload, objective='feasibility', probe_word=None):
    """Recover a bounded positive transfer on the faithful CAR cyclic space.

    Linear one-transfer equations construct candidates. Every supplied word,
    including nonlinear multi-transfer equations, is replayed before acceptance.
    Probe extrema explore ambiguity only and never select a physical source.
    """
    from finite_completion import exact_psd

    m = payload['m']
    if m > 2:
        raise ValueError('Faithful GNS recovery currently supports at most two modes')
    algebra = WordAlgebra(m)
    basis = [tuple(c) for k in range(2*m+1) for c in combinations(range(2*m), k)]
    index = {word: i for i, word in enumerate(basis)}
    dim = len(basis)
    values = _table(payload, algebra)
    if payload.get('linear_relations'):
        check_linear_relations(payload['linear_relations'],values,m,payload['degree'])
    if values.get(()) != 1:
        raise ValueError('Normalized identity moment must equal one')
    if any(word not in values for word in basis):
        raise ValueError('Complete CAR-only state moments are required')

    def state(word):
        phase, normal = algebra.reduce_word(word)
        return phase*values[normal]

    gram = sp.Matrix(dim, dim, lambda i, j: state(tuple(reversed(basis[i]))+basis[j]))
    gram_proof = exact_psd(gram)
    if gram_proof['rank'] != dim:
        raise ValueError('State is not faithful; a null-space quotient is required')
    gram_inverse = gram.inv()
    omega = sp.eye(dim)[:, 0]
    gammas = []
    for token in range(2*m):
        action = sp.zeros(dim)
        for j, word in enumerate(basis):
            phase, normal = algebra.reduce_word((token,)+word)
            action[index[normal], j] = phase
        gammas.append(action)
    quarter = sp.zeros(dim)
    for j, word in enumerate(basis):
        phase, normal = algebra.quarter(word)
        quarter[index[normal], j] = phase
    if quarter.adjoint()*gram*quarter != gram:
        raise ValueError('Supplied state is not quarter invariant')

    lower = _number(payload.get('transfer_floor', 0))
    upper = _number(payload.get('transfer_ceiling', 1))
    reference = _number(payload.get('reference_eigenvalue', 1))
    if not 0 <= lower <= reference <= upper or reference <= 0:
        raise ValueError('Invalid transfer/reference bounds')
    kind = payload.get('transfer_kind', 'bounded_positive')
    if kind not in ('thermal_gns', 'bounded_positive', 'vacuum_effect'):
        raise ValueError('Unknown transfer kind')
    if kind == 'vacuum_effect' and (upper != 1 or reference != 1):
        raise ValueError('vacuum_effect requires ceiling and reference one')
    if payload.get('cluster_gap', 0):
        raise ValueError('Pure-Fock cluster cap is not a thermal GNS condition')

    # A=G T is Hermitian in Euclidean coordinates; T is G-selfadjoint.
    a = cp.Variable((dim, dim), hermitian=True)
    gn = np.array(gram, dtype=complex)
    constraints = [a-float(lower)*gn >> 0, float(upper)*gn-a >> 0,
                   a[:, 0] == float(reference)*gn[:, 0]]
    qn = np.array(quarter, dtype=complex)
    constraints.append(qn.conj().T@a@qn == a)

    def linear_entry(word, matrix):
        phase, normal = algebra.reduce_word(word)
        if normal.count(T) != 1:
            raise ValueError('Probe must contain exactly one transfer')
        split = normal.index(T)
        left, right = normal[:split], normal[split+1:]
        # <Omega,left T right Omega> = phase(left*) A[left*,right].
        star_phase, left_star = algebra.star(left)
        return phase*star_phase*matrix[index[left_star], index[right]]

    determined_entries = set()
    for word, value in values.items():
        if word.count(T) == 1:
            constraints.append(linear_entry(word, a) == complex(value))
            split = word.index(T)
            _, left_star = algebra.star(word[:split])
            determined_entries.add((index[left_star], index[word[split+1:]]))
    if objective == 'feasibility':
        target = cp.Minimize(0)
    elif objective in ('min-probe', 'max-probe') and probe_word is not None:
        probe = cp.real(linear_entry(tuple(probe_word), a))
        target = cp.Minimize(probe) if objective == 'min-probe' else cp.Maximize(probe)
    else:
        raise ValueError('Choose feasibility or an explicitly labeled probe extremum')
    problem = cp.Problem(target, constraints)
    problem.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_gap_rel=1e-10,
                  tol_feas=1e-10, max_iter=300)
    if problem.status not in ('optimal', 'optimal_inaccurate') or a.value is None:
        raise ValueError('No numerical GNS completion: '+str(problem.status))
    exact_a = sp.Matrix(dim, dim, lambda i, j:
                        _number(float(a.value[i, j].real))
                        +sp.I*_number(float(a.value[i, j].imag)))
    checks = {'state_gram': gram_proof,
              'lower_bound': exact_psd(exact_a-lower*gram),
              'upper_bound': exact_psd(upper*gram-exact_a)}
    transfer = (gram_inverse*exact_a).applyfunc(sp.expand)
    if transfer*omega != reference*omega or transfer*quarter != quarter*transfer:
        raise ValueError('Exact reference or quarter replay failed')
    if transfer.adjoint()*gram != gram*transfer:
        raise ValueError('Exact GNS adjoint failed')
    identity = sp.eye(dim)
    for j, gamma in enumerate(gammas):
        if gamma.adjoint()*gram != gram*gamma:
            raise ValueError('Exact field adjoint failed')
        for k, other in enumerate(gammas):
            if gamma*other+other*gamma != (2*identity if j == k else sp.zeros(dim)):
                raise ValueError('Exact CAR failed')

    def moment(word):
        vector = omega
        for token in reversed(word):
            vector = (transfer if token == T else gammas[token])*vector
        return sp.expand((omega.adjoint()*gram*vector)[0])

    for word, value in values.items():
        if moment(word) != value:
            raise ValueError('Exact replay failed at '+str(word))
    probes = []
    for left in basis:
        for right in basis:
            word = left+(T,)+right
            phase, normal = algebra.reduce_word(word)
            value = moment(normal)
            probes.append({'word': list(normal), 'value': {
                're': str(sp.re(value)), 'im': str(sp.im(value))}})
    return {'schema': 'tfpt.faithful-car-gns-completion.v1', 'ok': True,
            'verdict': 'EXACT_FAITHFUL_FINITE_GNS_COMPLETION',
            'complete_finite_kernel_certificate': True,
            'representation': 'faithful_car_gns', 'dimension': dim,
            'objective': objective, 'objective_scope': 'ambiguity probe, not physical selection',
            'basis': [list(word) for word in basis], 'gram': _json_matrix(gram),
            'transfer': _json_matrix(transfer), 'exact_checks': checks,
            'all_input_moments_replayed': len(values),
            'multi_T_moments_replayed': sum(w.count(T) > 1 for w in values),
            'one_T_metric_entries_determined': len(determined_entries),
            'unique_in_fixed_faithful_representation': len(determined_entries) == dim*dim,
            'uniqueness_scope': 'Sufficient complete-entry certificate in the stipulated CAR-cyclic representation only; no larger-dilation or original-source uniqueness claim.',
            'candidate_predictions': probes, 'solver_status': problem.status,
            'input_sha256': sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'pure_state_assumed': False, 'wick_rule_assumed': False,
            'exactness_scope': 'Exact for supplied rational strings, or explicitly rational-recovered numeric moments; exact controls are retained after numerical residual checks. No exactness claim about unrounded floating input.',
            'original_source_derivation': False, 'absolute_partition_function_recovered': False,
            'scope': 'Declared finite faithful CAR cyclic source and bounded transfer; no original P1 map, physical time, continuum or T1-T8 closure.'}
