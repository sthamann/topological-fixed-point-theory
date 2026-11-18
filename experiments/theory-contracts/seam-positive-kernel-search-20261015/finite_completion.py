#!/usr/bin/env python3
"""Conditional minimal-Fock completion of a discovered truncated word kernel."""
from __future__ import annotations
import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
import cvxpy as cp
import numpy as np
import sympy as sp
from word_algebra import WordAlgebra, T
from linear_relations import check_linear_relations

TOLERANCE=1e-6
MAX_DENOMINATOR=64


def rational(value):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
        raise ValueError('Moment components must be finite JSON numbers')
    q=Fraction(value).limit_denominator(MAX_DENOMINATOR)
    if abs(float(q)-value)>TOLERANCE:
        raise ValueError(f'No rational with denominator <=64 within1e-6 of {value}')
    return sp.Rational(q.numerator,q.denominator)


def scalar_json(value):
    value=sp.expand(value)
    return {'re':str(sp.re(value)),'im':str(sp.im(value))}


def matrix_json(matrix):
    def part(value):
        q=sp.Rational(value)
        return [int(q.p),int(q.q)]
    return [[{'re':part(sp.re(matrix[row,col])),'im':part(sp.im(matrix[row,col]))}
             for col in range(matrix.cols)] for row in range(matrix.rows)]


def gamma_action(generator,state):
    mode=generator//2
    phase=(-1)**((state&((1<<mode)-1)).bit_count())
    if generator%2: phase*=(-sp.I if state&(1<<mode) else sp.I)
    return state^(1<<mode),phase


def clifford_action(word,state):
    phase=sp.Integer(1)
    for token in reversed(word):
        state,factor=gamma_action(token,state);phase*=factor
    return state,phase


def exact_psd(matrix):
    if matrix!=matrix.adjoint(): raise ValueError('Exact PSD matrix is not Hermitian')
    work=sp.Matrix(matrix)
    pivots=[]
    while work.rows:
        diagonal=[work[j,j] for j in range(work.rows)]
        if any(value.is_real is not True or value<0 for value in diagonal):
            raise ValueError('Negative/nonreal exact Schur pivot')
        pivot=next((j for j,value in enumerate(diagonal) if value>0),None)
        if pivot is None:
            if any(work): raise ValueError('Zero diagonal with nonzero exact row')
            break
        if pivot:
            indices=[pivot]+[j for j in range(work.rows) if j!=pivot]
            work=work.extract(indices,indices)
        value=work[0,0];pivots.append(str(value))
        work=(work[1:,1:]-work[1:,0]*work[0,1:]/value).applyfunc(sp.expand)
    return {'method':'exact pivoted Schur/LDL','rank':len(pivots),'pivots':pivots}


class Representation:
    def __init__(self,m,omega,transfer):
        self.m=m;self.dim=1<<m;self.omega=omega;self.transfer=transfer
        self.columns=[[] for _ in range(self.dim)]
        for (row,col),value in transfer.todok().items(): self.columns[col].append((row,value))

    @lru_cache(None)
    def vector(self,word):
        if not word: return {self.omega:sp.Integer(1)}
        source=self.vector(word[1:]);out={}
        if word[0]==T:
            for col,value in source.items():
                for row,factor in self.columns[col]: out[row]=out.get(row,0)+factor*value
        else:
            for col,value in source.items():
                row,factor=gamma_action(word[0],col)
                out[row]=out.get(row,0)+factor*value
        return {row:sp.expand(value) for row,value in out.items() if sp.expand(value)!=0}

    def moment(self,word): return self.vector(tuple(word)).get(self.omega,sp.Integer(0))


def _complete(payload,objective):
    # This recovery ansatz proves T Omega=Omega and T<=I. Other source
    # classes must not acquire that certificate through the direct API/CLI.
    if (payload.get('transfer_kind','vacuum_effect')!='vacuum_effect'
            or float(payload.get('reference_eigenvalue',1))!=1.0
            or float(payload.get('transfer_ceiling',1))!=1.0):
        raise ValueError('Finite Fock completion supports only vacuum_effect with reference_eigenvalue=1 and transfer_ceiling=1')
    objective={'mintrace':'min-trace','maxtrace':'max-trace'}.get(objective,objective)
    if objective not in ('min-trace','max-trace'): raise ValueError('objective must be min-trace or max-trace')
    m=payload['m'];algebra=WordAlgebra(m);dim=1<<m
    rows=payload['moments']
    if not rows: raise ValueError('No supplied moments')
    raw={}
    for row in rows:
        word=tuple(row['word']);phase,normal=algebra.reduce_word(word)
        if word!=normal or phase!=1: raise ValueError('Input moments must use canonical basis words')
        if word in raw: raise ValueError('Duplicate input word')
        raw[word]=row
    if () not in raw: raise ValueError('Missing normalization moment')
    occupations=[]
    for j in range(m):
        moment=raw.get((2*j,2*j+1))
        if moment is None: raise ValueError('Missing occupation-defining Majorana moment')
        n=(1-float(moment['im']))/2
        if min(abs(n),abs(n-1))>TOLERANCE: raise ValueError(f'Occupation {j} is not pure: {n}')
        occupations.append(int(n>.5))
    omega=sum(bit<<j for j,bit in enumerate(occupations))
    # Recover all inputs before certification, including every multi-T moment.
    moments={}
    for word,row in raw.items():
        try: moments[word]=rational(row['re'])+sp.I*rational(row['im'])
        except ValueError as error: raise ValueError(f'Input rational recovery failed at word {word}: {error}') from error
    if moments[()]!=1: raise ValueError('Exact normalization is not one')
    check_linear_relations(payload.get('linear_relations',[]),moments,m,payload['degree'])
    floor=rational(payload.get('transfer_floor',0))
    if not 0<=floor<=1: raise ValueError('Transfer floor outside [0,1]')
    gap=rational(payload.get('cluster_gap',0))
    if not 0<=gap<1: raise ValueError('cluster_gap must lie in [0,1)')
    if gap>0 and floor>1-gap:
        raise ValueError('Infeasible bounds: transfer_floor > 1-cluster_gap on the nonempty Omega-orthogonal complement')
    fixed={}
    for word,value in moments.items():
        if word.count(T)==0:
            state,coefficient=clifford_action(word,omega)
            expected=coefficient if state==omega else 0
            if expected!=value: raise ValueError(f'Pure CAR state disagrees at word {word}')
        elif word.count(T)==1:
            position=word.index(T);left=word[:position];right=word[position+1:]
            row,left_phase=clifford_action(tuple(reversed(left)),omega)
            col,right_phase=clifford_action(right,omega)
            coefficient=sp.conjugate(left_phase)*right_phase
            entry=sp.expand(value/coefficient)
            if (row,col) in fixed and fixed[row,col]!=entry:
                raise ValueError(f'Conflicting one-T constraints at word {word}')
            fixed[row,col]=entry
    transfer=cp.Variable((dim,dim),hermitian=True)
    constraints=[transfer-float(floor)*np.eye(dim)>>0,np.eye(dim)-transfer>>0]
    target=np.zeros(dim);target[omega]=1
    constraints.append(transfer[:,omega]==target)
    cluster_cap=(1-float(gap))*np.eye(dim)+float(gap)*np.outer(target,target)
    if gap>0: constraints.append(cluster_cap-transfer>>0)
    phases=np.asarray([state.bit_count()%4 for state in range(dim)])
    constraints.append(cp.multiply(phases[:,None]!=phases[None,:],transfer)==0)
    for (row,col),value in fixed.items(): constraints.append(transfer[row,col]==complex(value))
    trace=cp.real(cp.trace(transfer))
    problem=cp.Problem(cp.Minimize(trace) if objective=='min-trace' else cp.Maximize(trace),constraints)
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=300)
    if problem.status not in ('optimal','optimal_inaccurate') or transfer.value is None:
        raise ValueError('Finite completion SDP failed: '+str(problem.status))
    entries={}
    for row in range(dim):
        for col in range(dim):
            value=transfer.value[row,col]
            exact=rational(float(value.real))+sp.I*rational(float(value.imag))
            if exact: entries[row,col]=exact
    exact_t=sp.SparseMatrix(dim,dim,entries)
    if exact_t!=exact_t.adjoint(): raise ValueError('Recovered T is not exactly Hermitian')
    identity=sp.eye(dim)
    t_psd=exact_psd(exact_t)
    floor_psd=exact_psd(exact_t-floor*identity)
    upper_psd=exact_psd(identity-exact_t)
    e=sp.zeros(dim,1);e[omega]=1
    cluster_cap_psd=exact_psd((1-gap)*identity+gap*e*e.adjoint()-exact_t)
    finite_clustering=upper_psd['rank']==dim-1
    if exact_t*e!=e: raise ValueError('Exact stationarity T Omega=Omega failed')
    internal_quarter=sp.diag(*(sp.I**int(p) for p in phases))
    gns_phases=(phases-omega.bit_count())%4
    quarter=sp.diag(*(sp.I**int(p) for p in gns_phases))
    if quarter**4!=identity or quarter*exact_t!=exact_t*quarter:
        raise ValueError('Exact quarter covariance/order failed')
    if quarter*e!=e: raise ValueError('Normalized GNS quarter implementer does not fix Omega')
    # Check the complete original CAR matrices, not numerical anticommutators.
    gammas=[]
    for generator in range(2*m):
        elements={}
        for state in range(dim):
            destination,phase=gamma_action(generator,state);elements[destination,state]=phase
        gammas.append(sp.SparseMatrix(dim,dim,elements))
    for j,gamma in enumerate(gammas):
        if gamma.adjoint()!=gamma or gamma*gamma!=identity: raise ValueError('Exact Majorana check failed')
        expected=gammas[j+1] if j%2==0 else -gammas[j-1]
        if quarter*gamma*quarter.adjoint()!=expected: raise ValueError('Exact quarter field covariance failed')
        for previous in gammas[:j]:
            if gamma*previous+previous*gamma!=sp.zeros(dim): raise ValueError('Exact CAR check failed')
    representation=Representation(m,omega,exact_t)
    for word,value in moments.items():
        actual=representation.moment(word)
        if actual!=value:
            return {'verdict':'INCONCLUSIVE_MULTI_T_REPLAY_FAILED','objective':objective,'m':m,
                'ok':False,'dimension':dim,'exact_matrices':None,
                'replay_failure':{'word':list(word),'expected':scalar_json(value),'actual':scalar_json(actual)},
                'complete_finite_kernel_certificate':False,'extra_realization_ansatz':'minimal irreducible CAR Fock representation of dimension2^m'}
    predictions=[]
    for j,k in __import__('itertools').combinations(range(m),2):
        word=(2*j,2*k,T,2*k,2*j)
        predictions.append({'word':list(word),'value':scalar_json(representation.moment(word))})
    return {'schema':'tfpt.positive-word-kernel.finite-completion.v1',
        'verdict':'EXACT_FINITE_KERNEL_COMPLETION','complete_finite_kernel_certificate':True,
        'ok':True,
        'm':m,'dimension':dim,'objective':objective,'objective_scope':'min/max trace only explore undetermined directions, not a physical selection principle',
        'inferred_occupations':occupations,'Omega_basis_index':omega,
        'T':{'format':'sparse Gaussian rationals','entries':[
            {'row':row,'column':col,**scalar_json(value)} for (row,col),value in sorted(entries.items())]},
        'gammas':{'format':'monomial','definitions':'gamma_2j=c_j+c_j*, gamma_(2j+1)=i(c_j*-c_j)',
            'actions':[{'permutation':[gamma_action(j,state)[0] for state in range(dim)],
                        'phases':[scalar_json(gamma_action(j,state)[1]) for state in range(dim)]} for j in range(2*m)]},
        'R':{'format':'diagonal roots of unity','i_exponents':gns_phases.tolist(),
             'meaning':'normalized GNS implementer i^(N-N_Omega), fixes Omega'},
        'original_internal_W':{'format':'diagonal roots of unity','i_exponents':phases.tolist(),
                              'meaning':'original internal rotation i^N; need not fix Omega'},
        'exact_matrices':{'T':matrix_json(exact_t),'Omega':matrix_json(e),
                          'gammas':[matrix_json(gamma) for gamma in gammas],
                          'quarter':matrix_json(quarter)},
        'exact_checks':{'all_input_moments_replayed':len(moments),'multi_T_moments_replayed':sum(w.count(T)>1 for w in moments),
            'one_T_linear_constraints':len(fixed),'CAR':True,'stationarity':True,'quarter_order_and_covariance':True,'quarter_fixes_Omega':True,
            'T_PSD':t_psd,'floor_PSD':floor_psd,'I_minus_T_PSD':upper_psd,'cluster_cap_minus_T_PSD':cluster_cap_psd},
        'finite_clustering':finite_clustering,
        'stationary_eigenspace_dimension':dim-upper_psd['rank'],
        'clustering':{'declared_gap_bound':str(gap),'gap_was_explicit_input':'cluster_gap' in payload,
            'certified_gap_lower_bound':str(gap) if gap>0 else None,
            'bound':'T <= (1-gap)I+gap|Omega><Omega|',
            'limit':'T^n -> |Omega><Omega|' if finite_clustering else 'Additional stationary directions prevent the unique-Omega limit',
            'scope':'Exact finite-dimensional clustering analogue. The numerical gap is a search restriction, not a native physical scale or an OS-net identification.'},
        'candidate_predictions':predictions,'solver_status':problem.status,
        'finite_selfadjoint_H_minus_log_T_exists':t_psd['rank']==dim,
        'finite_unitary_field_time_is_UCP_automorphism':t_psd['rank']==dim,
        'time_scope':'If T is strictly positive, H=-log(T) and Ad_(T^(it)) define finite unitary field time, a UCP automorphism. Not established: a Euclidean UCP channel with GNS implementer T, a physical source selection, or an event instrument.',
        'flatness_claimed':False,'physical_source_selection_proved':False,
        'extra_realization_ansatz':'minimal irreducible CAR Fock representation of dimension2^m with the inferred pure occupation vector; this ansatz is not inferred from truncated flatness',
        'input_sha256':sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()}


def complete(payload,objective='min-trace',*,representation='minimal_fock',probe_word=None):
    if representation=='faithful_car_gns':
        from gns_completion import recover
        try:
            return recover(payload,objective,probe_word)
        except (ValueError,KeyError,TypeError,OverflowError,cp.error.SolverError) as error:
            return {'verdict':'INCONCLUSIVE_FAITHFUL_GNS_COMPLETION','objective':objective,
                    'ok':False,'dimension':None,'reason':str(error),
                    'representation':representation,'physical_source_selection_proved':False}
    if representation!='minimal_fock':
        raise ValueError('Unknown finite completion representation')
    try: return _complete(payload,objective)
    except (ValueError,KeyError,TypeError,OverflowError,cp.error.SolverError) as error:
        return {'verdict':'INCONCLUSIVE_FINITE_FOCK_COMPLETION','objective':objective,
                'ok':False,'dimension':None,'exact_matrices':None,
                'complete_finite_kernel_certificate':False,'reason':str(error),
                'extra_realization_ansatz':'minimal irreducible CAR Fock representation of dimension2^m; not a physical source selection'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--objective',choices=('min-trace','max-trace','mintrace','maxtrace'),default='min-trace')
    parser.add_argument('--representation',choices=('minimal_fock','faithful_car_gns'),default='minimal_fock')
    parser.add_argument('--gns-objective',choices=('feasibility','min-probe','max-probe'),default='feasibility')
    parser.add_argument('--probe-word',type=int,nargs='+')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=complete(json.loads(args.input.read_text()),
        args.gns_objective if args.representation=='faithful_car_gns' else args.objective,
        representation=args.representation,probe_word=args.probe_word)
    if args.output: args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:result.get(key) for key in ('verdict','complete_finite_kernel_certificate','objective','representation','reason')},sort_keys=True))


if __name__=='__main__': main()
