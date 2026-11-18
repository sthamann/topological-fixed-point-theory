"""Construct the unique finite CAR lift from a FULL creator covariance equation.

Input R is the supplied one-particle source transfer, not an inferred physical
seam. No target many-body matrix or energy list is used. The recurrence solves
Q c_j* = sum_i R_ij c_i* Q, Q|0>=|0>. This is stronger than vacuum moments.
"""
from __future__ import annotations
from functools import lru_cache
import argparse
import json
from pathlib import Path
import sympy as sp


def creator(m, j):
    dim=1<<m
    result=sp.zeros(dim)
    for state in range(dim):
        if not state & (1<<j):
            result[state|(1<<j),state]=(-1)**((state&((1<<j)-1)).bit_count())
    return result


def source_lift(one_particle):
    """Solve creator covariance column by column, with exact arithmetic."""
    r=sp.Matrix(one_particle)
    m=r.rows
    if r.cols != m or not 1<=m<=6:
        raise ValueError('A square one-particle matrix of dimension 1..6 is required')
    dim=1<<m
    cs=[creator(m,j) for j in range(m)]
    transported=[sum((r[i,j]*cs[i] for i in range(m)),sp.zeros(dim)) for j in range(m)]
    vacuum=sp.eye(dim)[:,0]

    @lru_cache(None)
    def column(state):
        if state==0: return vacuum
        # c_j*|state without j> = phase |state>. Every state is reachable.
        j=state.bit_length()-1
        remainder=state^(1<<j)
        phase=(-1)**((remainder&((1<<j)-1)).bit_count())
        return (phase*transported[j]*column(remainder)).applyfunc(sp.simplify)

    q=sp.Matrix.hstack(*(column(state) for state in range(dim)))
    if q*vacuum != vacuum:
        raise ArithmeticError('Vacuum normalization failed')
    for j in range(m):
        if (q*cs[j]-transported[j]*q).applyfunc(sp.simplify)!=sp.zeros(dim):
            raise ArithmeticError('Full creator covariance failed')
    return q


def covariance_residual(q, r):
    q,r=sp.Matrix(q),sp.Matrix(r)
    m=r.rows; dim=1<<m
    if q.shape != (dim,dim): raise ValueError('Wrong many-body matrix shape')
    failures=[]
    for j in range(m):
        right=sum((r[i,j]*creator(m,i) for i in range(m)),sp.zeros(dim))
        residual=(q*creator(m,j)-right*q).applyfunc(sp.simplify)
        if residual!=sp.zeros(dim): failures.append(j)
    return failures


def determinant_lift(r):
    """Independent exterior-minor formula, used only to check recurrence output."""
    r=sp.Matrix(r);m=r.rows;dim=1<<m
    q=sp.zeros(dim)
    for row in range(dim):
        ii=[j for j in range(m) if row&(1<<j)]
        for col in range(dim):
            jj=[j for j in range(m) if col&(1<<j)]
            if len(ii)==len(jj):
                q[row,col]=r.extract(ii,jj).det() if ii else 1
    return q


def checks():
    a=sp.Rational(1,4)
    r=sp.diag(a,a)
    q=source_lift(r)
    assert q==sp.diag(1,a,a,a*a)
    assert q==determinant_lift(r)
    # Interacting correction invisible in vacuum and one-particle propagation.
    alt=sp.diag(1,a,a,sp.Rational(1,2))
    same_first_columns=q[:,:3]==alt[:,:3]
    bad=covariance_residual(alt,r)
    assert same_first_columns and bad==[0,1]
    # Noncommuting transfers test actual chronological composition, not spectra.
    r1=sp.Matrix([[sp.Rational(1,2),sp.Rational(1,6)],[sp.Rational(1,6),sp.Rational(1,2)]])
    r2=sp.diag(sp.Rational(1,3),sp.Rational(2,3))
    assert r1*r2!=r2*r1
    assert source_lift(r2*r1)==source_lift(r2)*source_lift(r1)
    # Quarter source lift: full particle words acquire all their field phases.
    u=sp.diag(sp.I,-sp.I)
    assert source_lift(u)==sp.diag(1,sp.I,-sp.I,1)
    return {
        'schema':'tfpt.source-field-covariance-sewing.v1',
        'status':'EXACT_CONDITIONAL_MULTIFIELD_CONSTRUCTION',
        'free_control':{'one_field':str(a),'two_field':str(a*a),'spectrum':[str(x) for x in q.diagonal()]},
        'operator_covariance_rejects_ambiguous_extension':{'same_vacuum_and_one_particle_columns':same_first_columns,'failed_generator_indices':bad},
        'noncommuting_transfer_composition_exact':True,
        'quarter_lift_exact':True,
        'full_field_covariance_assumed':True,
        'full_field_covariance_derived_from_P1_P2':False,
        'proved':'Full creator covariance and reference normalization uniquely determine all CAR product actions and chronological sewing.',
        'scope':'Finite CAR lift; no physical source selection, local net construction or event instrument.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path)
    args=p.parse_args();result=checks()
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
