"""Literal quarter-fold current shift and complete quartic response.

Fixed original conformal point-splitting reference and existing M/O map.
No selected clock, rank2 cut, coupling, optimization or profile.
"""
from pathlib import Path
from collections import Counter
import json,itertools as it
import sympy as s
import check_mixed_return as m

Q=m.Q
def cl(A):return A.applyfunc(s.simplify)
def zero(A):return cl(A)==s.zeros(*A.shape)
def series_response(f,b,grade):
    """[t^grade] f(b+sum a(-n)t^n), actual full shifted field on sea."""
    result={}
    for mon,coef in f.items():
        p={():coef}
        for a in mon:
            out={}
            for word,c in p.items():
                out[word]=out.get(word,Q.zero)+c*b[a]
                weight=sum(1+k//8 for k in word)
                for n in range(1,grade-weight+1):
                    key=tuple(sorted(word+((n-1)*8+a,)))
                    out[key]=out.get(key,Q.zero)+c
            p={k:c for k,c in out.items() if c}
        for word,c in p.items():
            if sum(1+k//8 for k in word)==grade:
                result[word]=result.get(word,Q.zero)+c
    return {k:c for k,c in result.items() if c}

def main():
    ld=m._lattice_data();T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE))
    O=cl(T*ld['basis'].inv()/s.sqrt(2));U=s.zeros(8,4)
    for j in range(4):U[2*j,j]=1/s.sqrt(2);U[2*j+1,j]=-s.I/s.sqrt(2)
    m.need(zero(O.T*O-s.eye(8)) and zero(U.H*U-s.eye(4)),'Original O and actual positive polarization')
    fold=json.loads((m.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())
    M=s.Matrix(fold['explicit_Weyl_matrix']);m.need(M.T*M==s.eye(8),'Fixed original M')
    ps,ns=m.code();x=s.symbols('x0:4')
    polys=[sum(Q.to_sympy(c)*s.prod(x[j] for j in mon) for mon,c in f.items()) for f in ps]
    D=s.diag(*[s.sqrt(n) for n in ns])
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],
                  [0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48))
    # Exact finite Fourier matrix; alpha contributes a diagonal unitary.
    zeta=(1+s.I)/s.sqrt(2)
    fourier=s.Matrix(8,8,lambda r,k:s.expand(zeta**(r*k))/s.sqrt(8))
    m.need(zero(fourier.H*fourier-s.eye(8)),'Exact8channel DFT, alpha diagonal cancels')
    rows={}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        mu=s.Matrix([(k+alpha)/8 for k in range(8)]);delta=mu-s.ones(8,1)/2
        # All integer mode labels: source n=8l+k, negative iff l<=-1.
        for ell in range(-4,5):
            for k in range(8):
                m.need((8*ell+k+alpha<0)==(ell<=-1),'Negative-D chiral orientation: transported sea, no new copies/state')
                m.need(8*(ell+mu[k])==8*ell+k+alpha,'Actual physical source time factor8 retained')
        m.need(all(s.Rational(1,2)+delta[k]==mu[k] and s.Rational(1,2)-delta[k]==1-mu[k] for k in range(8)),
               'Original particle/hole gaps fix current-reference shift sign')
        d=M*delta; bn=cl(U.T*O.T*d);b=[Q.from_sympy(v) for v in bn]
        # Twisted conformal field one-point values, not vacuum charge in VE8.
        one=cl(D.inv()*s.Matrix([p.subs(dict(zip(x,bn)),simultaneous=True) for p in polys]))
        matches=[q for q in range(6) if s.Matrix.hstack(one,N[:,q]).rank()==1]
        responses=[series_response(f,b,4) for f in ps]
        m.need(all({mon:c for mon,c in f.items() if len(mon)==4}==ps[i] for i,f in enumerate(responses)),
               'Full native quartic component retained, all additional modes included')
        B=cl(s.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(responses[i],responses[j]))/s.sqrt(ns[i]*ns[j])))
        m.need(B==B.H,'Complete GNS response Gram is Hermitian')
        parts={tuple(sorted((1+k//8 for k in word),reverse=True)) for f in responses for word in f}
        m.need(parts=={(4,),(3,1),(2,2),(2,1,1),(1,1,1,1)},'All level4 oscillator partitions, no code compression')
        comm=[]
        for q in range(6):
            for pair in it.combinations([j for j in range(6) if j!=q],2):
                uv=[N[:,j]-(s.sqrt(6)-1)*N[:,q]/5 for j in pair]
                P=cl(sum([u*u.T for u in uv],s.zeros(5)))
                C=cl(B*P-P*B);dist=s.simplify(s.trace(C.H*C))
                comm.append({'q':q,'pair':list(pair),'commutes':dist==0,'commutator_norm_squared':str(dist)})
        # Exact characteristic polynomial, not a fitted eigenvalue ordering.
        y=s.Symbol('y');char=s.factor(B.charpoly(y).as_expr())
        rows[str(alpha)]={'mu':[str(v) for v in mu],'delta':[str(v) for v in delta],
             'd_alpha':[str(v) for v in d],'native_positive_shift_b':[str(v) for v in bn],
             'quartic_onepoint':m.mat(one),'onepoint_simplex_ray_matches':matches,
             'full_level4_response_gram':m.mat(B),'characteristic_polynomial':str(char),
             'response_monomial_counts':[len(f) for f in responses],
             'all_level4_partitions':[list(p) for p in sorted(parts)],
             'commuting_polar_flags':[v for v in comm if v['commutes']],
             'all_polar_commutators':comm,
             'shift_ground_conformal_weight':str((d.T*d)[0]/2)}
    out={'status':'PASS_EXACT_FOLD_CURRENT_TO_QUARTIC_RESPONSE','contract_id':'UR.NATIVE.JOINT_OPERATOR.01',
         'date':'2026-10-17','cases':rows,
         'scope':'Negative-DchiralorientationofexistingCARfold/sea,fixedconformalpointsplitreference,existingM/Oandnativequartics. OriginalQWZpositivehalflinerequires explicitcurrent/particle-holeorientationdictionary;oppositeshiftb->-bkeepsallhomogeneousquarticGrams. Current/neutralquarticpullback,notrawP1selectionornewoddCAR/E8spinextension. GNSmodeGramofliteralY_d(F)[z^0]state,noparticularlocalprofile/instrument. State-adaptedsubtractedfieldsareDIFFERENTlabels;normalorderingcountertermsmustbeexplicit. Directonepoint/spectralresponsecomparisononly,nogeneralstate/P1selectorno-go.'}
    Path(__file__).with_name('fold_quartic_transport_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'cases':{a:{'b':r['native_positive_shift_b'],'onepoint':r['quartic_onepoint'],
       'onepoint_matches':r['onepoint_simplex_ray_matches'],'commuting_flags':r['commuting_polar_flags'],
       'characteristic_polynomial':r['characteristic_polynomial'],'response_terms':r['response_monomial_counts']} for a,r in rows.items()}},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
