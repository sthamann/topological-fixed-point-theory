"""Search unknown positive CAR/transfer word moments; never input a target state.

A finite SDP is a relaxation, not a positive functional on all words. Exact
reconstruction lives in exact_certificate.py and is a separate acceptance step.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
import math
from fractions import Fraction
from pathlib import Path
import cvxpy as cp
import numpy as np
import scipy.sparse as sparse
from word_algebra import WordAlgebra, T
from linear_relations import normalize_linear_relations


OBJECTIVES=('feasibility','diversity','identity-control','random')
TRANSFER_KINDS=('vacuum_effect','thermal_gns','bounded_positive')
OBJECTIVE_DESCRIPTIONS={
    'feasibility':'Constant zero objective; no occupation or response is optimized.',
    'diversity':'Calibration/boundary probe: maximize summed single-field response deficits plus a seeded occupation tie-break. This can saturate the global transfer floor; it is not a physical selector.',
    'identity-control':'Calibration: minimize single-field response deficits toward reference_eigenvalue times the identity transfer.',
    'random':'Calibration: maximize a seeded random real linear functional of the truncated moments.'}


def normalize_fixed_moments(m,degree,rows=None):
    """Safe exact Gaussian-rational controls; no executable input expressions."""
    if rows is None: return []
    if not isinstance(rows,list): raise ValueError('fixed_moments must be a list')
    algebra=WordAlgebra(m);fixed={}
    def component(value):
        if isinstance(value,bool) or not isinstance(value,(int,float,str,Fraction)):
            raise ValueError('Fixed moment components must be finite numbers or rational strings')
        try:
            q=Fraction(str(value)) if not isinstance(value,Fraction) else value
            if not math.isfinite(float(q)): raise ValueError('Nonfinite fixed moment')
            return q
        except (ValueError,OverflowError,ZeroDivisionError) as error:
            raise ValueError('Invalid finite rational fixed moment component') from error
    for row in rows:
        if not isinstance(row,dict) or not {'word','re'}<=set(row) or set(row)-{'word','re','im'}:
            raise ValueError('Each fixed moment requires word/re and optional im only')
        if not isinstance(row['word'],(list,tuple)): raise ValueError('Fixed word must be a token list')
        sign,word=algebra.reduce_word(row['word'])
        if len(word)>2*degree: raise ValueError('Fixed word exceeds available moment degree; increase degree')
        target=(sign*component(row['re']),sign*component(row.get('im',0)))
        if word in fixed and fixed[word]!=target: raise ValueError('Conflicting fixed moments after exact CAR reduction')
        fixed[word]=target
    return [{'word':list(word),'re':str(value[0]),'im':str(value[1])}
            for word,value in sorted(fixed.items(),key=lambda item:(len(item[0]),item[0]))]


def load_source_control(path):
    if path is None: return {}
    path=Path(path);data=path.read_bytes();control=json.loads(data)
    if not isinstance(control,dict) or set(control)-{'fixed_moments','linear_relations','source_observables'}:
        raise ValueError('Source control accepts only fixed_moments/linear_relations/source_observables')
    metadata=control.get('source_observables',{})
    if not isinstance(metadata,dict): raise ValueError('source_observables must be a metadata object')
    return {'fixed_moments':control.get('fixed_moments',[]),
            'linear_relations':control.get('linear_relations',[]),
            'source_observables':{**metadata,'input_path':str(path.resolve()),
                                  'input_sha256':hashlib.sha256(data).hexdigest()}}


class MomentSearch:
    def __init__(self, m=5, degree=2, transfer_floor=0.0,*,transfer_ceiling=1.0,transfer_kind='vacuum_effect',reference_eigenvalue=1.0,fixed_moments=None,linear_relations=None,source_observables=None):
        if degree < 2:
            raise ValueError('degree >= 2 is required for a transfer response')
        if not math.isfinite(reference_eigenvalue) or reference_eigenvalue <= 0:
            raise ValueError('reference_eigenvalue must be finite and positive')
        if not math.isfinite(transfer_floor) or not 0 <= transfer_floor <= reference_eigenvalue:
            raise ValueError('transfer_floor must be finite, nonnegative and <=reference_eigenvalue')
        if not math.isfinite(transfer_ceiling) or transfer_ceiling < reference_eigenvalue:
            raise ValueError('transfer_ceiling must be finite and >=reference_eigenvalue')
        if transfer_kind not in TRANSFER_KINDS:
            raise ValueError('Unknown transfer_kind')
        if transfer_kind=='vacuum_effect' and (transfer_ceiling!=1 or reference_eigenvalue!=1):
            raise ValueError('vacuum_effect requires ceiling1 and reference_eigenvalue1; declare bounded_positive for a scaled source')
        self.algebra = a = WordAlgebra(m)
        self.m, self.degree, self.transfer_floor = m, degree, transfer_floor
        self.transfer_ceiling,self.transfer_kind=transfer_ceiling,transfer_kind
        self.reference_eigenvalue=reference_eigenvalue
        self.fixed_moments=normalize_fixed_moments(m,degree,fixed_moments)
        self.linear_relations=normalize_linear_relations(m,degree,linear_relations)
        if source_observables is not None and not isinstance(source_observables,dict):
            raise ValueError('source_observables must be a metadata object')
        self.source_observables=json.loads(json.dumps(source_observables or {},allow_nan=False))
        self.words = a.words(2 * degree)
        self.index = {w: i for i, w in enumerate(self.words)}
        self.y = cp.Variable(len(self.words), complex=True, name='unknown_moments')
        self.basis = a.words(degree)
        self.lower = a.words(degree - 1)
        self.gram = self.matrix(self.basis)
        self.lower_gram = self.matrix(self.lower)
        self.transfer = self.matrix(self.lower, middle=(T,))
        self.transfer_squared = self.matrix(self.lower, middle=(T,T))
        # Separate lower and upper localizers; no unjustified global commutation.
        self.constraints = [self.gram >> 0,
                            self.transfer - transfer_floor * self.lower_gram >> 0,
                            transfer_ceiling*self.lower_gram - self.transfer >> 0,
                            (transfer_ceiling+transfer_floor)*self.transfer - self.transfer_squared
                            - transfer_floor*transfer_ceiling*self.lower_gram >> 0,
                            self.y[self.index[()]] == 1]
        star_index, star_sign, quarter_index, quarter_sign = [], [], [], []
        stationarity_left, stationarity_right = [], []
        for w in self.words:
            s, v = a.star(w)
            star_index.append(self.index[v]); star_sign.append(s)
            s, v = a.quarter(w)
            quarter_index.append(self.index[v]); quarter_sign.append(s)
            # T Omega = lambda Omega gives a factor lambda for exterior T letters.
            # It does NOT erase an internal T and does NOT impose T = I.
            if w and w[0] == T:
                stationarity_left.append(self.index[w])
                stationarity_right.append(self.index[w[1:]])
            if w and w[-1] == T:
                stationarity_left.append(self.index[w])
                stationarity_right.append(self.index[w[:-1]])
        self.constraints += [
            cp.multiply(star_sign, self.y[star_index]) == cp.conj(self.y),
            cp.multiply(quarter_sign, self.y[quarter_index]) == self.y,
            self.y[stationarity_left] == reference_eigenvalue*self.y[stationarity_right]]
        for row in self.fixed_moments:
            value=complex(float(Fraction(row['re'])),float(Fraction(row['im'])))
            self.constraints.append(self.moment(row['word'])==value)
        for row in self.linear_relations:
            expression=0*self.y[0]
            for term in row['terms']:
                coefficient=complex(float(Fraction(term['re'])),float(Fraction(term['im'])))
                expression += coefficient*self.moment(term['word'])
            target=complex(float(Fraction(row['re'])),float(Fraction(row['im'])))
            self.constraints.append(expression==target)

    def moment(self, word):
        s, w = self.algebra.reduce_word(word)
        return s * self.y[self.index[w]]

    def matrix(self, left, right=None, middle=()):
        right = left if right is None else right
        rows, cols, coeff = [], [], []
        for i, u in enumerate(left):
            s, us = self.algebra.star(u)
            for j, v in enumerate(right):
                q, w = self.algebra.reduce_word(us + middle + v)
                rows.append(i * len(right) + j)
                cols.append(self.index[w]); coeff.append(s * q)
        transform = sparse.csr_matrix((coeff, (rows, cols)),
                                      shape=(len(left)*len(right), len(self.words)))
        return cp.reshape(transform @ self.y, (len(left),len(right)), order='C')

    def solve(self, objective='feasibility', seed=0, solver='CLARABEL',
              max_iters=20000, tolerance=1e-8):
        started = time.monotonic()
        # All objectives select probes in a feasible family; none is asserted
        # to be a TFPT source-selection law. Supplied control moments stay inputs.
        response = cp.real(sum(self.reference_eigenvalue-self.moment((j,T,j)) for j in range(2*self.m)))
        rng = np.random.default_rng(seed)
        if objective == 'feasibility':
            target=cp.Minimize(0)
        elif objective == 'diversity':
            # Signed occupation probes break ties without specifying a state.
            weights = rng.integers(-7,8,size=self.m).astype(float)
            weights[weights == 0] = 1
            state_probe = cp.real(sum(weights[j]*1j*self.moment((2*j,2*j+1))
                                      for j in range(self.m))) / (20*self.m)
            target = cp.Maximize(response + state_probe)
        elif objective == 'identity-control':
            target = cp.Minimize(response)
        elif objective == 'random':
            w = rng.normal(size=len(self.words)) + 1j*rng.normal(size=len(self.words))
            target = cp.Maximize(cp.real(w @ self.y)/len(w))
        else:
            raise ValueError('unknown objective')
        problem = cp.Problem(target, self.constraints)
        options = ({'eps':tolerance,'max_iters':max_iters} if solver == 'SCS'
                   else {'tol_gap_abs':tolerance,'tol_feas':tolerance,
                         'tol_gap_rel':tolerance,'max_iter':max_iters})
        solver_error = None
        try:
            problem.solve(solver=solver, **options)
        except cp.error.SolverError as error:
            solver_error = str(error)
        result = {'schema':'tfpt.unknown-positive-word-kernel.v1',
                  'm':self.m,'degree':self.degree,'seed':seed,
                  'objective':objective,'transfer_floor':self.transfer_floor,
                  'transfer_ceiling':self.transfer_ceiling,'transfer_kind':self.transfer_kind,
                  'reference_eigenvalue':self.reference_eigenvalue,
                  'reference_eigenvalue_scope':'Explicit supplied eigenvalue in T Omega=lambda Omega, not a value derived from P1 or an absolute source normalization',
                  'conditional_normalized_state':True,
                  'absolute_partition_function_recovered':False,
                  'state_normalization_scope':'omega(1)=1 normalizes the conditional state only; it is not a recovered Z(0)',
                  'transfer_scope':'Declared bounded positive reference-eigenvector transfer class; ceiling1 is a vacuum/effect restriction, not universal thermal GNS time. A supplied source outside these bounds is a class mismatch, not physical exclusion.',
                  'objective_role':('feasibility_only_no_state_or_response_optimization' if objective=='feasibility'
                                    else 'explicit_calibration_or_boundary_probe_not_native_selector'),
                  'objective_definition':OBJECTIVE_DESCRIPTIONS[objective],
                  'fixed_moments':self.fixed_moments,'fixed_moment_count':len(self.fixed_moments),
                  'linear_relations':self.linear_relations,'linear_relation_count':len(self.linear_relations),
                  'source_observables':self.source_observables,
                  'source_observables_scope':'Supplied control data and provenance labels, not an inferred original source',
                  'solver':solver,'solver_status':problem.status or 'solver_error',
                  'solver_error':solver_error,
                  'elapsed_seconds':time.monotonic()-started,
                  'basis_size':len(self.basis),'moment_count':len(self.words),
                  'assumptions':['finite CAR interface with m modes',
                                 'selfadjoint extra transfer operation T',
                                 'T bounded between declared transfer_floor and transfer_ceiling',
                                 'T Omega = reference_eigenvalue Omega (explicit input)',
                                 'invariance under internal CAR quarter rotation'],
                  'not_assumed':['Hamiltonian matrix','density matrix','pure binary occupation state',
                                 '3+2 occupation split','physical chiral holonomy'],
                  'selection_scope':'objective, seed and global transfer_floor are search choices, not P1/P2 consequences; feasibility imposes no optimization of occupations or transfer responses',
                  'verdict':'NO_NUMERICAL_CANDIDATE',
                  'engine_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        if solver_error or self.y.value is None:
            return result
        y = self.y.value
        gram = np.asarray(self.gram.value)
        low = np.asarray(self.lower_gram.value)
        local = np.asarray(self.transfer.value)
        local_squared = np.asarray(self.transfer_squared.value)
        eigen = lambda x: np.linalg.eigvalsh((x+x.conj().T)/2)
        rank = lambda x: int(np.count_nonzero(eigen(x)>1e-6*max(1,float(np.linalg.norm(x,2)))))
        result.update({
            'verdict':'NUMERICAL_PSEUDOMOMENT_CANDIDATE_REQUIRES_EXACT_RECOVERY',
            'objective_value':float(problem.value),
            'response_sum':float(response.value),
            'numeric_rank_lower':rank(low),'numeric_rank_full':rank(gram),
            'numeric_flatness_candidate':rank(low)==rank(gram),
            'min_eigenvalues':{'gram':float(eigen(gram)[0]),
                              'lower_localizer':float(eigen(local-self.transfer_floor*low)[0]),
                              'upper_localizer':float(eigen(self.transfer_ceiling*low-local)[0])},
            'max_constraint_violation':max(float(np.max(c.violation())) for c in self.constraints),
            'occupations':[(1+float((1j*self.moment((2*j,2*j+1))).value.real))/2
                           for j in range(self.m)],
            'moments':[{'word':list(w),'re':float(y[i].real),'im':float(y[i].imag)}
                       for i,w in enumerate(self.words)]})
        result['min_eigenvalues']['contraction_localizer'] = float(eigen(
            (self.transfer_ceiling+self.transfer_floor)*local-local_squared
            -self.transfer_floor*self.transfer_ceiling*low)[0])
        return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--m',type=int,default=5)
    parser.add_argument('--degree',type=int,default=2)
    parser.add_argument('--seed',type=int,default=7)
    parser.add_argument('--objective',choices=OBJECTIVES,default='feasibility')
    parser.add_argument('--source-control',type=Path,help='JSON fixed_moments/linear_relations/source_observables; explicit supplied controls')
    parser.add_argument('--transfer-floor',type=float,default=0)
    parser.add_argument('--transfer-ceiling',type=float,default=1)
    parser.add_argument('--reference-eigenvalue',type=float,default=1)
    parser.add_argument('--transfer-kind',choices=TRANSFER_KINDS,default='vacuum_effect')
    parser.add_argument('--solver',choices=['CLARABEL','SCS'],default='CLARABEL')
    parser.add_argument('--max-iters',type=int,default=20000)
    parser.add_argument('--tolerance',type=float,default=1e-8)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists(): parser.error('Output already exists; choose a fresh file to preserve historical runs')
    search=MomentSearch(args.m,args.degree,args.transfer_floor,transfer_ceiling=args.transfer_ceiling,
                        transfer_kind=args.transfer_kind,reference_eigenvalue=args.reference_eigenvalue,
                        **load_source_control(args.source_control))
    out=search.solve(args.objective,args.seed,args.solver,args.max_iters,args.tolerance)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='moments'},indent=2))


if __name__=='__main__': main()
