"""Whole R104 Gamma(r) boundary words by original renewal and exterior minors.

No inverse, no new measure, no new Hamiltonian. The fixed original dimension
is 16. Null-P traces use at most 2^16 original bit states, not a model search.
"""
from functools import lru_cache
from pathlib import Path
import json
import math
import sys

import numpy as np
from scipy import linalg,sparse
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0,str(Path(__file__).parent))
import twined_original_det16 as source
original=source.original


class ExteriorBoundary:
    def __init__(self,r):
        self.r=np.asarray(r,dtype=complex)
        if self.r.shape!=(16,16):raise ValueError('original sixteen modes required')

    @lru_cache(maxsize=None)
    def minor(self,row_bits,column_bits):
        if row_bits.bit_count()!=column_bits.bit_count():return 0j
        if not row_bits:return 1+0j
        rows=[i for i in range(16) if row_bits>>i&1]
        cols=[i for i in range(16) if column_bits>>i&1]
        return complex(np.linalg.det(self.r[np.ix_(rows,cols)]))

    def split(self,letters):
        positions=[i for i,x in enumerate(letters) if x==('G',self)]
        if len(positions)>1:raise ValueError('one original boundary Gamma(r) per word')
        if not positions:return None
        j=positions[0]
        return letters[:j],letters[j+1:]

    def phi(self,letters):
        pieces=self.split(letters)
        if pieces is None:return original.word_phi(letters)
        left,right=pieces
        left_adjoint=[(i,not dagger) for i,dagger in reversed(left)]
        endpoints=[(0,1/math.sqrt(2)),(original.D-1,np.exp(1j*original.phi)/math.sqrt(2))]
        answer=0j
        for x,cx in endpoints:
            target,al=original.apply_word(left_adjoint,x)
            if target is None:continue
            for y,cy in endpoints:
                state,ar=original.apply_word(right,y)
                if state is not None:
                    answer+=np.conjugate(cx*al)*cy*ar*self.minor(target,state)
        return answer

    def trace(self,letters):
        pieces=self.split(letters)
        if pieces is None:return original.word_fock_trace(letters)
        left,right=pieces
        word=right+left # ordinary trace cyclicity; no CAR letter is commuted.
        if not word:return complex(np.linalg.det(np.eye(16)+self.r))
        if sum(1 if dagger else -1 for _,dagger in word)!=0:return 0j
        bits={}
        for mode,dagger in reversed(word):
            if mode not in bits:bits[mode]=0 if dagger else 1
        base=sum(bit<<mode for mode,bit in bits.items())
        free=[i for i in range(16) if i not in bits]
        total=0j
        for choice in range(1<<len(free)):
            state=base+sum(((choice>>j)&1)<<mode for j,mode in enumerate(free))
            target,amplitude=original.apply_word(word,state)
            if target is not None:total+=amplitude*self.minor(state,target)
        return total

    def word(self,letters,spans):
        return source.generalized_renewal(list(letters)+[('G',self)],list(spans)+[0.],phi_eval=self.phi,trace_eval=self.trace)


def native_word_matrix(letters,spans):
    result=original.I
    for index,t in enumerate(spans):
        result=result@(math.exp(-2*t)*original.I-math.expm1(-2*t)*original.P)
        if index<len(letters):
            i,dagger=letters[index]
            result=result@(original.cre[i] if dagger else original.ann[i])
    result.eliminate_zeros()
    return result.tocoo()


def native_trace_by_actual_generator(Q,many_body_generator,angle):
    # Direct canonical Fock unitary from the ACTUAL native dGamma(G), using
    # sparse exponential action. No exterior minors are used in this check.
    rows=sorted(set(int(x) for x in Q.row))
    basis=np.zeros((original.D,len(rows)),dtype=complex)
    basis[rows,np.arange(len(rows))]=1
    evolved=expm_multiply(1j*angle*many_body_generator,basis,traceA=0)
    lookup={row:i for i,row in enumerate(rows)}
    return sum(value*evolved[col,lookup[int(row)]] for row,col,value in zip(Q.row,Q.col,Q.data)),len(rows)


def main():
    original.initialize_native_controls()
    generators,pairs=original.native.positive_weyl_generators()
    # Same original v1002 diagnostic angle; a gauge insertion, no new coupling.
    angle=.2718281828459045
    index=pairs.index((0,2))
    G=generators[index]
    r=linalg.expm(1j*angle*G)
    many=original.native.second_quantize(G)
    rows,cols=np.nonzero(np.abs(G-np.diag(np.diag(G)))>1e-12)
    i,j=int(rows[0]),int(cols[0])
    off=[(i,True),(j,False)]+[(m,d) for m in range(16) if m not in (i,j) for d in (True,False)]
    cases=[('density32',[(m,d) for m in range(16) for d in (True,False)]),
           ('offdiagonal30',off),('omega16',[(m,False) for m in range(16)])]
    boundary=ExteriorBoundary(r)
    checks={
        'actual_noncentral_spin10_unitary':bool(np.max(np.abs(r.conj().T@r-np.eye(16)))<2e-14),
        'determinant_one':abs(np.linalg.det(r)-1)<2e-14,
        'noncentral_actual_generator':bool(np.max(np.abs(r-np.eye(16)*r[0,0]))>.01),
    }
    output=[];maximum=0.
    for name,letters in cases:
        spans=[(k+1)/sum(range(1,len(letters)+2)) for k in range(len(letters)+1)]
        Q=native_word_matrix(letters,spans)
        direct,native_columns=native_trace_by_actual_generator(Q,many,angle)
        computed=boundary.word(letters,spans)
        error=float(abs(direct-computed));maximum=max(maximum,error)
        checks[name+'_native_sparse_exponential']=error<2e-11*max(1,abs(direct))
        output.append({'word':name,'length':len(letters),'native':source.zcomplex(direct),'minor_renewal':source.zcomplex(computed),'error':error,'native_word_nnz':int(Q.nnz),'native_exponential_columns':native_columns})
    # General algorithm also covers original singular det(I+r) center exactly.
    singular=ExteriorBoundary(-np.eye(16,dtype=complex))
    checks['singular_empty_character_zero']=abs(singular.word([],[0.]))<1e-14
    checks['singular_marked_full_density_one']=abs(singular.word(cases[0][1],[0.]*(len(cases[0][1])+1))-1)<1e-14
    # Compare complete generic-minor evaluation to independently derived q-reader.
    for label,q in source.TWISTS.items():
        central=ExteriorBoundary(q*np.eye(16,dtype=complex))
        for name,letters in cases:
            spans=[1/(len(letters)+1)]*(len(letters)+1)
            checks[label+'_'+name+'_matches_original_center_algorithm']=abs(central.word(letters,spans)-source.twined_word(letters,spans,q))<2e-11
    checks={k:bool(v) for k,v in checks.items()};assert all(checks.values()),{k:v for k,v in checks.items() if not v}
    report={'status':'pass','source':'Actual original R104 rank-one source and native v1002 positive-Weyl Spin10 generator pair(0,2); original diagnostic angle .2718281828459045.',
            'scope':'Constructive complete finite R104 boundary Gamma(r) word reader for actually given canonical r. Not the full v1013 gauge embedding, primitive S_adm or a new selected measure.',
            'whole_word_formula':'Same original cyclic renewal. Each P-block with Gamma(r) is a four-endpoint sum of det(r[T,S]); no-P term is the finite original trace sum over CAR-admissible bit states. No determinant inverse.',
            'scaling':'For fixed original16 modes, at most2^16 native bit states per no-P trace, polynomial work in word length in the renewal. No scaling claim uniform in mode number or exact bit complexity.',
            'checked':sum(checks.values()),'total':len(checks),'max_absolute_error':maximum,'checks':checks,'values':output}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
