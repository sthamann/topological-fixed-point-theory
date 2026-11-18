"""Original two-mode exterior-transfer relations as source-word constraints.

Conditional finite alpha=1/4 particle/hole control, not a physical source law.
"""
import json
import sys
from pathlib import Path

import cvxpy as cp
import sympy as s

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parent
COLLAR = SEARCH.parent / 'seam-collar-multifield-20261016'
sys.path[:0] = [str(SEARCH), str(COLLAR)]

from charged_pair_control import RESPONSE, majoranas, moment, source_lift
from kernel_search import MomentSearch
from linear_relations import normalize_linear_relations
from source_sewing import creator


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def relation_terms(left, right, mode, a):
    """Twice <Omega| left (T c_i^* - a_i c_i^* T) right |Omega>."""
    terms = []
    for g, coefficient in ((2 * mode, s.Integer(1)), (2 * mode + 1, -s.I)):
        for middle, c in (((-1, g), coefficient), ((g, -1), -a * coefficient)):
            terms.append({'word': list(left + middle + right),
                          're': str(s.re(c)), 'im': str(s.im(c))})
    return {'terms': terms, 're': '0'}


def residual(row, transfer):
    return s.expand(sum((s.Rational(t['re']) + s.I*s.Rational(t['im'])) *
                        moment(t['word'], transfer) for t in row['terms']))


def main():
    # Relative filled-sea creators: gamma_0, gamma_2 generate the two
    # orthonormal one-particle vectors. The double state is gamma_2 gamma_0 Ω.
    ket_words = [(), (0,), (2,), (2, 0)]
    relations = [relation_terms(tuple(reversed(bra)), ket, i, RESPONSE[i, i])
                 for i in range(2) for bra in ket_words for ket in ket_words]
    relations = normalize_linear_relations(2, 3, relations)
    original = source_lift(RESPONSE)
    certificate = json.loads((HERE/'run-20261102-verified/exact-completion.json').read_text())
    def exact_entry(x):
        return s.Rational(*x['re']) + s.I*s.Rational(*x['im'])
    prior = s.Matrix([[exact_entry(x) for x in row] for row in certificate['exact_matrices']['T']])
    need(prior == s.diag(1, s.Rational(1, 4), s.Rational(1, 64), 0), 'stored prior candidate mismatch')
    need(original == s.diag(1, s.Rational(1, 4), s.Rational(1, 64), s.Rational(1, 256)), 'original source mismatch')
    source_residual = [residual(r, original) for r in relations]
    prior_residual = [residual(r, prior) for r in relations]
    need(all(v == 0 for v in source_residual), 'source violates relation')
    need(sum(v != 0 for v in prior_residual) == 2, 'prior violation count mismatch')

    # Completeness in the exact four-dimensional source representation:
    # the 32 basis sandwiches are exactly both full matrix equations.
    unknown = s.Matrix(4, 4, s.symbols('t0:16'))
    omega = s.Matrix([1, 0, 0, 0])
    equations = list(unknown*omega-omega)
    for i in range(2):
        c = creator(2, i)
        equations.extend(list(unknown*c - RESPONSE[i, i]*c*unknown))
    solution = s.linsolve(equations, list(unknown))
    need(solution == s.FiniteSet(tuple(original)), 'symbolic source uniqueness mismatch')

    # Same search as R085: only the two vacuum CAR moments and these operator
    # relations are supplied. Its internal quarter is not the charged clock.
    search = MomentSearch(m=2, degree=3, transfer_floor=0,
                          fixed_moments=[{'word':[0,1], 're':'0', 'im':'1'},
                                         {'word':[2,3], 're':'0', 'im':'1'}],
                          linear_relations=relations,
                          source_observables={'scope':'declared free alpha=1/4 particle/hole control'})
    probe = cp.real(search.moment((0,2,-1,2,0)))
    bounds = []
    for objective in (cp.Minimize(probe), cp.Maximize(probe)):
        problem = cp.Problem(objective, search.constraints)
        value = problem.solve(solver='SCS', eps=1e-8, max_iters=50000)
        bounds.append({'status': problem.status, 'value': value})

    output = {'scope':'Fixed two-mode free particle/hole control; not an interacting or integrated-H source; absolute Z0 and physical origin remain open.',
              'basis_kets': [list(x) for x in ket_words], 'relations': relations,
              'independent_matrix_equations': 2, 'basis_sandwiches':len(relations),
              'source_nonzero_residuals':sum(x != 0 for x in source_residual),
              'prior_nonzero_residuals':sum(x != 0 for x in prior_residual),
              'prior_residuals': [str(x) for x in prior_residual if x != 0],
              'relation_normalization': 'twice the creator equation; divide these residuals by two for creator-equation residuals',
              'prior_source': 'run-20261102-verified/exact-completion.json',
              'exact_matrix_solution': [str(original[i,i]) for i in range(4)],
              'exact_target_value': '1/256', 'numerical_min_max':bounds}
    (HERE/'generator-relations-result.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='relations'},indent=2))


if __name__ == '__main__':
    main()
