"""Full native local quartic OPE; intrinsic10-ray and6-simplex reconstruction.

The original simplex/mark is loaded ONLY for comparison after reconstruction.
No root selection, family omission, intermediate projection or fitted map.
"""
from pathlib import Path
from collections import Counter
import itertools as it,math,json
import sympy as s
import check_mixed_return as m
from check_conformal_erasure_carry import pd
Q=m.Q

def clean(A):return A.applyfunc(s.simplify)
def factor(mon):return math.prod(math.factorial(v) for v in Counter(mon).values())
def multi_diff(f,mon):
    for k in mon:f=pd(f,k)
    return f
def shift_to_negative(f):return {tuple(k+4 for k in mon):m.conj(c) for mon,c in f.items()}

def main():
    ps,ns=m.code();fields={};grams={};spec={}
    for r in range(5):
        fs=[]
        for i in range(5):
            for j in range(5):
                fs.append(m.add(*[(Q.one/factor(a),m.mul(multi_diff(ps[i],a),shift_to_negative(multi_diff(ps[j],a))))
                                  for a in it.combinations_with_replacement(range(4),r)]))
        gram=s.Matrix(25,25,lambda k,l:Q.to_sympy(m.inner(fs[k],fs[l]))/s.sqrt(ns[k//5]*ns[k%5]*ns[l//5]*ns[l%5]))
        fields[r]=fs;grams[r]=gram;spec[str(r)]={str(v):mult for v,mult in gram.eigenvals().items()}
    m.need(grams[0]==s.eye(25),'Full uncontracted local quartic-pair component encodes all25logical matrix components')
    m.need(grams[1]==4*s.eye(25),'Original erasure20 full matched-index return has Gram4I25')
    m.need(grams[2].eigenvals()=={s.Integer(18):1,s.Integer(8):9,s.Integer(0):15},'First nontrivial2body source OPE spectrum')
    # Three contractions leave exactly canonicalT times delta_ij.
    T=m.add(*[(Q.one,{(a,a+4):Q.one}) for a in range(4)])
    for i in range(5):
        for j in range(5):
            m.need(fields[3][5*i+j]==(m.add((Q(ns[i]),T)) if i==j else {}),'Actual native OPE gives originalT, not imported')
            m.need(fields[4][5*i+j]==({():Q(ns[i])} if i==j else {}),'Leading original OPE vacuum normalization')
    erased=[[pd(f,a) for a in range(4)] for f in ps];decoder_checks=0
    for i in range(5):
        for j in range(5):
            for k in range(5):
                for b in range(4):
                    terms=[]
                    for a in range(4):
                        d=m.add(*[(m.conj(c),multi_diff(erased[k][b],mon)) for mon,c in erased[j][a].items()])
                        terms.append((Q.one/ns[j],m.mul(erased[i][a],d)))
                    m.need(m.add(*terms)==(erased[i][b] if j==k else {}),'Literal r1 matched mode component gives actual decoderEij tensorI4')
                    decoder_checks+=1
    # Flatten the ACTUAL code at2+2; no simplex, target algebra or mark input.
    pairs=list(it.combinations_with_replacement(range(4),2))
    hs=[]
    for i in range(5):
        hs.append(s.Matrix(10,10,lambda a,b:Q.to_sympy(m.inner(ps[i],{tuple(sorted(pairs[a]+pairs[b])):Q.one}))/s.sqrt(6*ns[i]*factor(pairs[a])*factor(pairs[b]))))
    for a in range(5):
        for b in range(5):
            m.need(clean(hs[a]*hs[b]-hs[b]*hs[a])==s.zeros(10),'All actual pair flattenings commute')
            m.need(s.simplify(s.trace(hs[a]*hs[b]))==int(a==b),'Source flattening HS normalization')
    blocks=[s.eye(10)]
    for mat in hs:
        new=[]
        for B in blocks:
            small=clean(B.T*mat*B)
            for val,mult,vecs in small.eigenvects():
                on=s.GramSchmidt(vecs,orthonormal=True)
                new.append(clean(B*s.Matrix.hstack(*on)))
        blocks=new
    m.need(len(blocks)==10 and all(B.cols==1 for B in blocks),'Ten unique simultaneous source eigenlines')
    U=clean(s.Matrix.hstack(*blocks));m.need(clean(U.T*U)==s.eye(10),'No fitted transformation in simultaneous diagonalization')
    H=clean(s.Matrix(10,5,lambda t,i:(U.T*hs[i]*U)[t,t]))
    m.need(clean(H.T*H)==s.eye(5),'Actual eigenvalue rows form tight frame')
    gram=clean(H*H.T);sixgram=clean(6*gram)
    m.need(all(gram[i,i]==s.Rational(1,2) and all(abs(gram[i,j])==s.Rational(1,6) for j in range(10) if j!=i) for i in range(10)),'Derived ten rays have fixed1/3 absolute normalized overlaps')
    # Verify full C2 source-field spectral expansion in its ON basis.
    Bfields=[m.add(*[(Q.from_sympy(s.simplify(U[a,t]*U[b,t]/s.sqrt(factor(pairs[a])*factor(pairs[b])))),
                      {tuple(sorted(pairs[a]+tuple(k+4 for k in pairs[b]))):Q.one})
                     for a in range(10) for b in range(10) if U[a,t] and U[b,t]]) for t in range(10)]
    for i in range(5):
        for j in range(5):
            reconstructed=m.add(*[(Q.from_sympy(s.simplify(6*s.sqrt(ns[i]*ns[j])*H[t,i]*H[t,j])),Bfields[t]) for t in range(10)])
            m.need(reconstructed==fields[2][5*i+j],'Literal full local pair-OPE ten-ray decomposition')
    # Intrinsic4th moment: >=2/5, equality is flat squared overlaps.
    # Exhaust ALL sign choices; no numerical extremum search/choice of potential.
    solutions=[]
    intgram=[[int(x) for x in sixgram.row(i)] for i in range(10)]
    for values in it.product((1,-1),repeat=10):
        if all(sum(intgram[i][j]*values[j] for j in range(10))==6*values[i] for i in range(10)):
            solutions.append(values)
    m.need(len(solutions)==12,'Exhaustive1024 signs give exactly12 equality points')
    signs=[s.Matrix(v) for v in solutions if v[0]==1]
    W=[clean(H.T*v/s.sqrt(10)) for v in signs]
    m.need(len(W)==6,'Six intrinsic antipodal simplex axes')
    oriented=[W[0]]
    for v in W[1:]:
        overlap=s.simplify((W[0].T*v)[0]);m.need(abs(overlap)==s.Rational(1,5),'Simplex axis overlap')
        oriented.append(v if overlap<0 else -v)
    S=s.Matrix.hstack(*oriented)
    m.need(clean(S.T*S)==s.Rational(6,5)*s.eye(6)-s.ones(6)/5,'Recovered regular simplex in actual native coordinates')
    m.need(clean(S*s.ones(6,1))==s.zeros(5,1),'No chosen center/direction for recovered simplex')
    # Essential control: these are minima on the intrinsic projective REAL
    # locus, not all complex minima or a physical real-state restriction.
    complex_w=clean((W[0]+s.I*W[1])/s.sqrt(2))
    overlaps=clean(s.sqrt(2)*H*complex_w)
    m.need(s.simplify((complex_w.H*complex_w)[0])==1,'Complex control is unit norm')
    m.need(all(s.simplify(v*s.conjugate(v))==s.Rational(1,5) for v in overlaps),'Additional complex flat-overlap minimizer is retained')
    m.need(all(clean(complex_w*complex_w.H-w*w.T)!=s.zeros(5) for w in W),'Complex minimum is not falsely identified with six real rays')
    # COMPARISON ONLY starts here. Nothing above uses originalN or P2.
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=clean(s.diag(*[s.sqrt(v) for v in (4,12,12,12,24)])*raw/s.sqrt(48))
    matching=[]
    for w in W:
        found=[j for j in range(6) if clean(w*w.T-s.Rational(6,5)*N[:,j]*N[:,j].T)==s.zeros(5)]
        m.need(len(found)==1,'OriginalW6 simplex axes match LITERALLY; no rotation fit')
        matching.append(found[0])
    m.need(sorted(matching)==list(range(6)),'All six original rays recovered, not selected subset')
    old=json.loads((m.BASE/'polar_fibre_results.json').read_text());e0=s.Matrix(old['fixed_basis']['e0']);e1=s.Matrix(old['fixed_basis']['e1']);a=s.Matrix(old['fixed_basis']['a'])
    cp=-s.sqrt(s.Rational(3,5))*e0+s.sqrt(s.Rational(2,5))*e1
    P=clean(a*a.H+e0*e0.H+e1*e1.H-cp*cp.H)
    distances=[]
    for i,j in it.combinations(range(6),2):
        pair=N[:,[i,j]];plane=clean(pair*(pair.T*pair).inv()*pair.T)
        distance=s.simplify(s.trace((P-plane)**2));m.need(distance>0,'Original polarP2 is no source simplex pair plane')
        distances.append(distance)
    minimum=min(distances,key=lambda x:float(x))
    m.need(minimum==(21-6*s.sqrt(6))/25,'Exact original polar-mark comparison')
    out={'status':'PASS_EXACT_FULL_QUARTIC_LOCAL_OPE_AND_INTRINSIC_SIMPLEX','date':'2026-10-17',
      'contract_id':'UR.NATIVE.JOINT_OPERATOR.01',
      'full_wick_identity':'F_i(z)F_j†(w)=sum_r=0^4 (z-w)^(-2r) sum_|alpha|=r :partial_alpha f_i(u(z)) partial_alpha conjugate(f_j)(ubar(w)): /alpha!, with all Taylor terms retained',
      'primitive_contraction_component_weights':{'0':8,'1':6,'2':4,'3':2,'4':0},
      'source_gram_spectra_by_contraction_number':spec,
      'erasure_component':'r1= sum_a g_ia(u) conjugate(g_ja)(ubar); Gram4 I25; its3creator3annihilator mode part actsEij tensorI4 on actual20 at grade3',
      'literal_erasure_decoder_mode_checks':decoder_checks,
      'stress_component':'r3 gives delta_ij T_E8, r4 gives delta_ij vacuum; no assumed separateT',
      'pair_component':'r2 has source-rank10; actual symmetric pair flattenings commute, with ten distinct simultaneous eigenlines',
      'derived_pair_eigenvalue_rows':[[str(v) for v in H.row(t)] for t in range(10)],
      'equiangular_unit_ray_overlap_squared':'1/9','tight_frame_sum':'sum_t v_t v_t^T=2I5',
      'intrinsic_fourth_moment_bound':'sum_t (v_t.w)^4 >=2/5 on real unitw; equality iff all squared overlaps1/5',
      'exhaustive_sign_candidates':1024,'exact_equality_points':len(solutions),'antipodal_axes':len(W),
      'recovered_simplex_gram':'6/5 I6-1/5 ones6','reconstruction_uses_original_simplex':False,
      'intrinsic_real_form':'Ten source rays spanC5 and have pairwise nonzero overlaps. Antiunitary fixing every ray exists and is unique up to global phase, hence defines intrinsic projective real locus. No physical real-state restriction.',
      'complex_minimum_control':'PASS: (w1+i w2)/sqrt2 is an additional complex minimum; continuous complex minimum family is NOT excluded or selected away',
      'literal_original_simplex_axis_permutation':matching,'fitted_rotation':False,
      'source_derived_rank2_simplex_pair_planes':15,
      'original_polar_P2_distance_squared_minimum':str(minimum),
      'original_polar_P2_pair_plane_distance_multiplicities':dict(Counter(str(x) for x in distances)),
      'positive_success':'Previously separately supplied unmarked nativeW6 simplex AXES are reconstructed from full actual local quartic OPE on its intrinsically defined projective real locus, then match literal original coordinates. Overallorientation/permutation remain; not all complex minima, physical real-state selection or positive polar transport/mark.',
      'scope':'Within existing E8 current realization and actual quartic Fock embedding. Equality-set minimization is an algebraic reconstruction, not a postulated physical potential or state selection. Wick contraction is field algebra, not measurement or state projection. Individual Fourier-mode components are not locally supported interventions. No3+2 selected, no rawP1 realization, no common physical time, no fermion/Yukawa assignment.'}
    Path(__file__).with_name('full_quartic_ope_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
