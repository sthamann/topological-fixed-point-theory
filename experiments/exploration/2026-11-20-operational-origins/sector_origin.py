#!/usr/bin/env python3
"""Exact test and ablations of a minimal fixed-observable algebra proposal.

Two nonisomorphic noncommutative simple sectors are ASSUMPTIONS. This test
exposes their role; it does not derive them or physical superselection.
Only the explicitly constructed bistochastic conditional expectation is used.
"""
from itertools import combinations
import json
import sympy as s


def matrix_unit(n, i, j):
    m = s.zeros(n)
    m[i, j] = 1
    return m


def sector_universe(max_carrier=12):
    return [dict(p=p, q=q, carrier=p+q, algebra_dimension=p*p+q*q,
                 clauses={*(('both_noncommutative',) if p >= 2 else ()),
                          *(('nonisomorphic',) if p != q else ())})
            for p in range(1, max_carrier)
            for q in range(p, max_carrier-p+1)]


def minima(objects, clauses):
    eligible = [o for o in objects if set(clauses) <= o['clauses']]
    n = min(o['carrier'] for o in eligible)
    return [(o['p'], o['q']) for o in eligible if o['carrier'] == n]


def conditional_expectation(p=2, q=3):
    n = p+q
    P = s.diag(*([1]*p+[0]*q))
    Q = s.eye(n)-P
    E = lambda x: P*x*P+Q*x*Q
    units = [matrix_unit(n,i,j) for i in range(n) for j in range(n)]
    assert E(s.eye(n)) == s.eye(n)
    assert all(E(E(x)) == E(x) and s.trace(E(x)) == s.trace(x) for x in units)
    superoperator = s.Matrix.hstack(*(s.Matrix(n*n,1,list(E(x))) for x in units))
    choi = sum((s.kronecker_product(x,E(x)) for x in units),s.zeros(n*n))
    eig = choi.eigenvals()
    assert all(v >= 0 for v in eig)
    # This central generator is normalized to determinant phases +/-theta.
    # It is NOT a measured electromagnetic normalization.
    Y = -P/s.Integer(p)+Q/s.Integer(q)
    assert s.trace(Y) == 0
    return {'carrier':n,'fixed_algebra_dimension':superoperator.rank(),
            'choi_spectrum':{str(k):int(v) for k,v in eig.items()},
            'idempotent_unital_trace_preserving':True,
            'central_traceless_generator_by_determinant_convention':list(map(str,Y.diagonal())),
            'physical_hypercharge_normalization_derived':False}


def unital_cp_counterexample():
    # Unital CP + idempotent is not sufficient for inherited multiplication.
    K0=s.diag(1,1,0)
    K1=matrix_unit(3,2,0)/s.sqrt(2)
    K2=matrix_unit(3,2,1)/s.sqrt(2)
    E=lambda x:sum((K*x*K.conjugate().T for K in (K0,K1,K2)),s.zeros(3))
    B=s.diag(1,-1,0)
    assert E(s.eye(3)) == s.eye(3)
    assert all(E(E(matrix_unit(3,i,j))) == E(matrix_unit(3,i,j))
               for i in range(3) for j in range(3))
    assert E(B) == B and E(B*B) != B*B
    return {'unital_and_idempotent':True,'fixed_B':str(B),
            'B_squared':str(B*B),'E_B_squared':str(E(B*B)),
            'trace_preserving':False,'inherited_subalgebra':False}


def partition_average():
    # Giving the five slots every 3+2 cut with equal weight does not mark one.
    marks=list(combinations(range(5),2))
    ys=[s.diag(*(s.Rational(1,2) if i in pair else -s.Rational(1,3)
                 for i in range(5))) for pair in marks]
    mean=sum(ys,s.zeros(5))/len(ys)
    second=sum((y*y for y in ys),s.zeros(5))/len(ys)
    assert mean == s.zeros(5) and second == s.eye(5)/6
    return {'marks':len(marks),'mean_Y_diagonal':list(map(str,mean.diagonal())),
            'mean_Y_squared_diagonal':list(map(str,second.diagonal())),
            'scope':'linear classical averaging of all cuts; not conditioned or coherent branch dynamics'}


def main():
    universe=sector_universe()
    clauses=['both_noncommutative','nonisomorphic']
    rules=[]
    for k in range(3):
        for rule in combinations(clauses,k):
            rules.append({'rule':list(rule),'minimum_pairs':minima(universe,rule)})
    unconditioned=[o for o in universe if o['p'] >= 2]
    unequal=[o for o in unconditioned if o['p'] < o['q']]
    selected=minima(universe,clauses)
    assert selected == [(2,3)]
    return {'principle':'smallest faithful carrier with two nonisomorphic noncommutative simple sectors',
            'maximum_carrier_in_census':12,'all_two_sector_algebras':len(universe),
            'both_noncommutative_count':len(unconditioned),
            'both_noncommutative_nonisomorphic_count':len(unequal),
            'ablation':rules,'conditional_expectation':conditional_expectation(),
            'unital_cp_negative_control':unital_cp_counterexample(),
            'uniform_mark_negative_control':partition_average(),
            'stage':1,'verdict':'OFFEN: exact conditional algebraic P2 pattern, no independent sector-selection mechanism',
            'first_free_choice':'two inequivalent genuinely quantum sectors, then minimal carrier'}


if __name__ == '__main__':
    print(json.dumps(main(),indent=2,ensure_ascii=False))
