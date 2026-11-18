"""Actual source W4 and full quartic state responses, no closure assumption."""
import json
from pathlib import Path
import sympy as s
import check_mixed_return as m
from check_fold_quartic_transport import cl,zero,series_response

def main():
    ld=m._lattice_data();T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE));O=cl(T*ld['basis'].inv()/s.sqrt(2))
    M=s.Matrix(json.loads((m.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
    P=s.zeros(8)
    for j in range(8):P[(j+4)%8,j]=1
    G=cl(O.T*M*P*M.T*O);U=s.zeros(8,4)
    for j in range(4):U[2*j,j]=1/s.sqrt(2);U[2*j+1,j]=-s.I/s.sqrt(2)
    J=s.zeros(8)
    for j in range(4):J[2*j,2*j+1]=-1;J[2*j+1,2*j]=1
    m.need(G*G==s.eye(8) and G.T*G==s.eye(8),'Actual oscillator part of original W4')
    A=cl(U.H*G*U);B=cl(U.T*G*U)
    m.need(A==s.diag(0,-1,1,0) and B==s.diag(s.I,0,0,-s.I),'Literal actual polarization mixing, not fitted U4')
    def rotate(f):
        out={}
        for mon,c in f.items():
            phase=m.Q.one;word=[]
            for k in mon:
                n,a=divmod(k,8);m.need(a<4,'Initial shifted quartic positive-mode input')
                aa=(a+4 if a in (0,3) else a)
                phase*=m.Q.from_sympy(s.sympify([s.I,-1,1,-s.I][a]));word.append(8*n+aa)
            w=tuple(sorted(word));out[w]=out.get(w,m.Q.zero)+phase*c
        return {w:c for w,c in out.items() if c}
    ps,ns=m.code();scale=s.diag(*[s.sqrt(n) for n in ns]);a=s.Matrix([1]*4+[0]*4);ac=M*a
    control=s.Matrix(5,5,lambda i,j:m.Q.to_sympy(m.inner(ps[i],rotate(ps[j]))))
    cn=cl(scale.inv()*control*scale.inv())
    m.need(cn==s.diag(s.Rational(1,2),0,0,s.Rational(1,2),0),'Unshifted homogeneous chart has apparent rank2 overlap')
    rows={}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        d=M*s.Matrix([(j+alpha)/8-s.Rational(1,2) for j in range(8)])
        b=cl(U.T*O.T*d);bout=cl(U.T*O.T*(d+ac))
        incoming=[series_response(f,[m.Q.from_sympy(v) for v in b],4) for f in ps]
        outgoing=[series_response(f,[m.Q.from_sympy(v) for v in bout],4) for f in ps]
        actual=[rotate(f) for f in incoming]
        originalGram=s.Matrix(5,5,lambda i,j:m.Q.to_sympy(m.inner(incoming[i],incoming[j])))
        actualGram=s.Matrix(5,5,lambda i,j:m.Q.to_sympy(m.inner(actual[i],actual[j])))
        m.need(actualGram==originalGram,'Actual W4 preserves full state Gram, including exterior polarization/modes')
        cross=s.Matrix(5,5,lambda i,j:m.Q.to_sympy(m.inner(outgoing[i],actual[j])))
        m.need(cross.rank()==5,'Actual source state restores full rank overlap; no inherited rank2 kernel')
        increment=s.simplify((ac.T*ac)[0]/2+(d.T*ac)[0])
        m.need(increment==s.Rational(3,4)+alpha/2,'Original source W4 Ward time increment retained')
        rows[str(alpha)]={'incoming_shift':[str(v) for v in b],'outgoing_shift':[str(v) for v in bout],
             'cross_gram_raw':m.mat(cross),'cross_rank':cross.rank(),'cross_determinant':str(s.factor(cross.det())),
             'source_time_increment':str(8*increment),'fold_time_increment':str(increment),
             'full_W4_output_gram_equals_input':True}
    out={'status':'PASS_EXACT_ORIGINAL_WINDING_FULL_STATE_QUARTIC_RESPONSE','date':'2026-10-17',
         'native_real_oscillator_map':m.mat(G),'positive_to_positive':m.mat(A),'positive_to_negative':m.mat(B),
         'commutes_with_native_J':G*J==J*G,'unshifted_quartic_overlap':m.mat(cn),
         'unshifted_overlap_rank':2,'actual_charge_after_W4':[str(v) for v in ac],'cases':rows,
         'scope':'ActualoriginalglobalW4,existingevenCAR/E8factorization,existingM/OandfixedNSfielddictionary. CrossGramofactualfullY_d(F)[z0]inputandliteralY_d(F)[z0]outgoingcharge-a fieldresponses;notanexecutedprojection/closedcode/instrumentorautonomousclock. Outgoingcharge/modes/polarizationretained. Theapparent2+3rankofhomogeneousvacuumchartisnotactualstate-response rank. Noallhistory/physical3+2no-go.'}
    Path(__file__).with_name('fold_primitive_winding_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'Jcommutes':out['commutes_with_native_J'],
       'unshifted_rank':2,'actual_cases':{a:{k:r[k] for k in ('cross_rank','cross_determinant','source_time_increment')} for a,r in rows.items()}},indent=2))
if __name__=='__main__':main()
