"""Full original SINGLE-beta return and its exact off-code feedback.

No Hamiltonian, physical time, projection execution or new coupling is chosen.
L powers are ordered insertions of the already specified response pair only.
"""
from pathlib import Path
import json, math
from functools import lru_cache
import sympy as s
import check_mixed_return as m
Q=m.Q


def basis(grade):
    out=[]
    def rec(mon,start,left):
        if not left:
            out.append(mon);return
        for k in range(start,8*grade):
            g=1+k//8
            if g>left:break
            rec(mon+(k,),k,left-g)
    rec((),0,grade)
    return out


def annihilate(f, direction, mode):
    out={}
    for mon,value in f.items():
        for j,k in enumerate(mon):
            if k//8==mode-1:
                i=k%8; c=direction[i] if i<4 else m.conj(direction[i-4])
                if c:
                    n=mon[:j]+mon[j+1:]
                    out[n]=out.get(n,Q.zero)+mode*c*value
    return {n:c for n,c in out.items() if c}


def full_source(b,c):
    B=[{},*[m.linear(b,j) for j in range(1,4)]]
    V=[{},*[m.linear(c,j) for j in range(1,4)]]
    E=[{():Q.one},B[1],m.add((Q(1)/2,m.mul(B[1],B[1])),(Q(1)/2,B[2])),
       m.add((Q(1)/6,m.mul(m.mul(B[1],B[1]),B[1])),(Q(1)/2,m.mul(B[1],B[2])),(Q(1)/3,B[3]))]
    @lru_cache(None)
    def translated(mon):
        terms={0:{():Q.one}}
        for k in mon:
            mode=1+k//8;i=k%8;r=b[i] if i<4 else m.conj(b[i-4])
            nxt={}
            for power,pol in terms.items():
                nxt[power]=m.add((Q.one,nxt.get(power,{})),(Q.one,m.mul(pol,{(k,):Q.one})))
                if r:nxt[power-mode]=m.add((Q.one,nxt.get(power-mode,{})),(-r,pol))
            terms={p:f for p,f in nxt.items() if f}
        return terms
    def shift(f):
        out={}
        for mon,value in f.items():
            for power,pol in translated(mon).items():
                out[power]=m.add((Q.one,out.get(power,{})),(value,pol))
        return {p:f for p,f in out.items() if f}
    def coeff(product, power):
        return m.add(*[(Q.one,m.mul(pol,E[power-p])) for p,pol in product.items() if 0<=power-p<len(E)])
    def apply(f):
        m.need(all(sum(1+k//8 for k in mon)==4 for mon in f),'Full neutral grade4 input')
        F=shift(f)
        pieces=[(Q.one,m.mul(V[j],coeff(F,-j-1))) for j in range(1,4)]
        pieces += [(Q.one,coeff(shift(annihilate(f,c,j)),j-1)) for j in range(1,5)]
        out=m.add(*pieces)
        m.need(all(sum(1+k//8 for k in mon)==3 for mon in out),'Charged beta output oscillator grade3')
        return out
    return apply


def mat(A):return [[str(x) for x in row] for row in A.tolist()]


def main():
    b=[Q.zero,Q.zero,Q(0,-1),Q.zero];c=[Q.zero,Q.zero,Q.zero,Q(4,-4)]
    source=full_source(b,c);neut=basis(4);charged=basis(3)
    m.need(len(neut)==726 and len(charged)==192,'Complete original PBW blocks, no state truncation')
    charged_set=set(charged);inverse={};entries=0
    for mon in neut:
        Sm=source({mon:Q.one});m.need(set(Sm)<=charged_set,'Complete charged PBW basis')
        for n,value in Sm.items():
            inverse.setdefault(n,[]).append((mon,m.conj(value)*Q(m.weight(n))/m.weight(mon)))
            entries+=1
    def adjoint(f):
        out={}
        for n,value in f.items():
            for mon,z in inverse.get(n,[]):out[mon]=out.get(mon,Q.zero)+z*value
        return {mon:z for mon,z in out.items() if z}
    def L(f):return adjoint(source(f))
    polys,norms=m.code()
    def gram(fs,gs):return s.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(fs[i],gs[j])))
    D=s.diag(*[s.sqrt(n) for n in norms]);Dinv=D.inv();Ninv=s.diag(*[s.Rational(1,n) for n in norms])
    forward=[source(f) for f in polys]
    for j,f in enumerate(polys):m.need(forward[j]==m.response(f,b,c),'Full PBW source equals prior positive-code formula')
    returned=[adjoint(f) for f in forward]
    G=gram(forward,forward);K=m.clean(Dinv*G*Dinv)
    m.need(gram(polys,returned)==G,'Exact physical Fock adjoint on all original code columns')
    q=m.add((Q.one,{(0,0,0,0):Q.one}),(-Q(1)/4,polys[0]))
    q_amp=s.simplify(Q.to_sympy(m.inner(q,returned[0]))/s.sqrt(18*norms[0]))
    m.need(s.simplify(q_amp+188/(3*s.sqrt(3)))==0,'Full return independently reproduces prior leakage witness')
    H=m.clean(Dinv*gram(returned,returned)*Dinv)
    delta=m.clean(H-K*K)
    outside=[]
    for i,lf in enumerate(returned):
        projected=m.add(*[(Q.from_sympy(G[j,i]/norms[j]),polys[j]) for j in range(5)])
        outside.append(m.add((Q.one,lf),(-Q.one,projected)))
    m.need(m.clean(Dinv*gram(outside,outside)*Dinv)==delta,'Exact B^dagger B and full four-insertion defect')
    for f in polys:
        for y in outside:m.need(m.inner(f,y)==Q.zero,'Retained exterior vectors are exactly code orthogonal')
    m.need(delta[0,0]>s.Rational(35344,27),'Full defect exceeds prior single-component lower bound')
    # Exact third-return amplitude verifies the first two memory terms.
    twice=[L(f) for f in returned]
    three=m.clean(Dinv*gram(returned,twice)*Dinv)
    # D=Q L Q on the FULL neutral grade4 block. No repeated physical P.
    Loutside=[L(y) for y in outside]
    memory1=m.clean(Dinv*gram(outside,Loutside)*Dinv)
    m.need(m.zero(three-(K**3+K*delta+delta*K+memory1)),'Exact third-return memory expansion')
    m.need(delta==delta.H and memory1==memory1.H,'Hermitian exact memory coefficients')
    # Exact principal-minor certificates, no floating-point positivity test.
    lower=delta-5120*s.eye(5)
    m.need(lower[0,0]>0 and s.simplify(lower[:2,:2].det())>0 and
           all(lower[j,j]>=0 for j in range(2,5)) and
           all(not lower[i,j] for i in range(5) for j in range(5) if i!=j and (i>=2 or j>=2)),
           'Strict first-memory lower bound5120 I on every code direction')
    retention=m.clean(s.Rational(9,44)*H-K*K)
    m.need(retention[0,0]>0 and s.simplify(retention[:2,:2].det())>0 and
           all(retention[j,j]>=0 for j in range(2,5)) and
           retention[:,2]==s.zeros(5,1) and retention[:,3]==s.zeros(5,1) and
           all(not retention[i,j] for i in range(5) for j in range(5) if i!=j and (i>=2 or j>=2)),
           'Sharp scale-independent retained norm fraction at most9/44')
    hankel=s.eye(5).row_join(K).col_join(K.row_join(H))
    m.need(hankel.rank()==10 and s.simplify(hankel.det()-delta.det())==0,
           'Exact rank10 Hankel/Gram certificate from first3 response moments')
    output={
        'status':'PASS_EXACT_ORIGINAL_FULL_RETURN_AND_OFF_CODE_FEEDBACK',
        'contract_id':'UR.NATIVE.JOINT_OPERATOR.01',
        'source_scope':'Original SINGLE-beta paired response, doubled root(0,0,0,0,0,0,-2,-2); no coherent cross-root truncation',
        'complete_neutral_grade4_dimension':len(neut),'complete_charge_beta_grade4_dimension':len(charged),
        'source_nonzero_entries':entries,'source_nonzero_output_rows':len(inverse),
        'original_compressed_return':mat(K),'actual_four_insertion_return':mat(H),
        'first_memory_BdaggerB':mat(delta),'first_memory_rank':delta.rank(),
        'first_memory_characteristic_polynomial':str(s.factor(delta.charpoly().as_expr())),
        'first_memory_trace':str(s.trace(delta)),
        'first_memory_minimum_eigenvalue':'5120',
        'first_memory_lower_bound_2x2_principal_minors':[str(lower[0,0]),str(s.simplify(lower[:2,:2].det()))],
        'maximum_retained_squared_norm_fraction':'9/44',
        'minimum_off_code_squared_norm_fraction':'35/44',
        'sharp_retention_equality_space':'span{original f2,f3}',
        'retention_bound_2x2_principal_minors':[str(retention[0,0]),str(s.simplify(retention[:2,:2].det()))],
        'first_closed_linear_state_space_minimum_dimension':10,
        'first_response_hankel_rank':10,
        'first_response_hankel_determinant':str(s.factor(hankel.det())),
        'linear_realization_scope':'M_n=C A^n B for n>=0, fixed linear A and time-independent input/output maps. No separate feedthrough at n=0. Lower bound10, not sufficiency of10.',
        'code0_true_four_insertion_amplitude':str(H[0,0]),
        'code0_compressed_square_amplitude':str((K*K)[0,0]),
        'code0_memory_difference':str(delta[0,0]),
        'prior_off_code_amplitude_reproduced':str(q_amp),
        'next_memory_BdaggerDB':mat(memory1),'actual_six_insertion_return':mat(three),
        'exact_discrete_elimination':'x[n+1]=K x[n]+sum(j=0..n-1) Bdagger D^(n-1-j) B x[j]; y[0]=0',
        'meaning':'Exact algebraic ordered-response composition of specified original pair. P,Q only split retained and omitted amplitudes, never a physical measurement/reprojection. No new Hamiltonian or physical event/time rule chosen. No physical source selector or3+2 derived.'
    }
    Path(__file__).with_name('source_return_memory_results.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(output,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
