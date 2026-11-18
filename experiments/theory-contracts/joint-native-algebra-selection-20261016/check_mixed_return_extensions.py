"""All original family omissions and the actual trimer erasure return test."""
import check_mixed_return as m
from pathlib import Path
import json
import sympy as s

def main():
    Q=m.Q;need=m.need;clean=m.clean;zero=m.zero;mat=m.mat
    L=m._lattice_data();T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE));O=clean(T*L['basis'].inv()/s.sqrt(2))
    native=m._native_neutral_branch();roots=m._actual_current_source()['chevalley'].roots
    ri=clean(s.sqrt(2)*O.T/2)
    native_roots=[ri*s.Matrix(r) for r in roots];need(all(all(t.is_Integer for t in r) for r in native_roots),'All actual native integer roots')
    polys,norms=m.code();D=s.diag(*[s.sqrt(n) for n in norms]);sums=[]
    prior=json.loads(Path(__file__).with_name('mixed_return_results.json').read_text())
    for a in range(4):
        subset=[r for r in native['spinor'] if tuple(r[5:])==tuple(native['family'][a])]
        A16=sum((s.Matrix(r)*s.Matrix(r).T/4 for r in subset),s.zeros(8));AN=clean(O.T*(4*s.eye(8)-A16)*O)
        total=s.zeros(5)
        for r in native_roots:
            vv=-AN*r+((r.T*AN*r)[0]/4)*r
            b=[Q(int(r[2*j]),-int(r[2*j+1]))/2 for j in range(4)]
            c=[Q.from_sympy(vv[2*j])-Q(0,1)*Q.from_sympy(vv[2*j+1]) for j in range(4)];c=[v/2 for v in c]
            W=[m.response(f,b,c) for f in polys]
            for i in range(5):
                for j in range(i,5):
                    v=Q.to_sympy(m.inner(W[i],W[j]));total[i,j]+=v
                    if i!=j:total[j,i]+=s.conjugate(v)
        K=clean(D.inv()*total*D.inv());sums.append(K);print('family',a,'sum',mat(K),flush=True)
    need(sums[0]==s.Matrix(prior['ordered_return_sum']),'Generic quartic expression reproduces full prior sum')
    # Actual original decoder columns are the annihilated normalized code
    # vectors. Check the whole20-column Fock isometry, not a dimension match.
    derivs=[m.deriv(f,[Q(int(j==a)) for j in range(4)]) for f in polys for a in range(4)]
    for i in range(20):
        for j in range(20):
            need(m.inner(derivs[i],derivs[j])==(Q(norms[i//4]) if i==j else Q.zero),'Actual derivative/erasure isometry')
    # Original simple source-root witness, transported by the literal O.
    source_root=(0,0,0,0,0,0,-2,-2);need(source_root in roots,'Original witness root')
    b=[Q.zero,Q.zero,Q(0,-1),Q.zero];c=[Q.zero,Q.zero,Q.zero,Q(4,-4)]
    W5=[m.response(f,b,c) for f in polys];G5=s.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(W5[i],W5[j])));K5=clean(D.inv()*G5*D.inv())
    old=next(v for v in prior['individual_return_grams'] if tuple(v['root_doubled'])==source_root)
    need(G5==s.Matrix(old['raw_gram']),'Independent original source witness agrees with all-root calculation')
    W20=[m.response(f,b,c) for f in derivs];G20=s.Matrix(20,20,lambda i,j:Q.to_sympy(m.inner(W20[i],W20[j])))
    D20=s.diag(*[s.sqrt(norms[i//4]) for i in range(20)]);K20=clean(D20.inv()*G20*D20.inv())
    need(not W20[0] and K20[:,0]==s.zeros(20,1),'Actual erased |0,0> trimer is annihilated')
    expected=s.kronecker_product(K5,s.eye(4));res=clean(K20-expected);hs=s.simplify(s.trace(res.H*res))
    M20=s.diag(*[s.Rational(1,norms[i//4]) for i in range(20)])*G20
    reshuffle=s.Matrix(25,16,lambda row,col:M20[4*(row//5)+col//4,4*(row%5)+col%4]);schmidt=reshuffle.rank()
    out={'status':'PASS_EXACT_ALL_FAMILY_RETURNS_AND_NATIVE_ERASURE_TEST','all_four_family_return_sums':[mat(A) for A in sums],
         'all_four_sums_equal':all(A==sums[0] for A in sums),'witness_root_doubled':list(source_root),
         'actual_erasure_isometry':'g_(i,a)=partial_a f_i; all20 normalized columns orthonormal; annihilator=sqrt4 times original tensor contraction',
         'witness_quartic_return':mat(K5),'witness_trimer_return':mat(K20),'quartic_rank':K5.rank(),'trimer_rank':K20.rank(),
         'trimer_operator_schmidt_rank':schmidt,'direct_decoder_lift_HS_residual_squared':str(hs),
         'zero_witness':'g_(0,0)=x0^3/sqrt6 is killed by actual source response, hence K20 g00=0',
         'decoder_lift_same_vector_diagonal':str(K5[0,0]),
         'scope':'Actual native neutral quartic/trimer Fock realizations, not the differently charged E8 F20 block. Excludes these single ordered response branches as a common full-erasure factorized lift. Does not exclude other combined histories, instruments or source-selected channels; no physical projection inserted.'}
    Path(__file__).with_name('mixed_return_extensions_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('all_four_family_return_sums','witness_quartic_return','witness_trimer_return')},indent=2))
if __name__=='__main__':main()
