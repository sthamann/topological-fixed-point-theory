"""Clock-free actual quartic product geometry, state and composition controls.

Actual source polynomials only until explicit six-simplex comparison.
No target P2/clock/real-state input, no physically assumed selection rule.
"""
from pathlib import Path
from collections import Counter
import itertools as it
import math,json
import sympy as s
import check_mixed_return as m
from check_full_quartic_ope import multi_diff,shift_to_negative,factor

Q=m.Q
def cl(A):return A.applyfunc(s.simplify)
def zero(A):return cl(A)==s.zeros(*A.shape)
def spectrum(A):return {str(k):v for k,v in A.eigenvals().items()}
def symvec(v,ij):return s.Matrix([(1 if i==j else s.sqrt(2))*v[i]*v[j] for i,j in ij])
def symrep(T,J):return cl(J.T*s.kronecker_product(T,T)*J)
def coeffmat(polys,vs):
    mons=[vs[i]*vs[j] for i,j in it.combinations_with_replacement(range(len(vs)),2)]
    return s.Matrix([[s.Poly(p,*vs).coeff_monomial(mon) for mon in mons] for p in polys])

def main():
    ps,ns=m.code();ij=list(it.combinations_with_replacement(range(5),2))
    pp=[m.mul(ps[i],ps[j]) for i,j in ij]
    c=[s.sqrt(s.Rational(1 if i==j else 2,ns[i]*ns[j])) for i,j in ij]
    H=s.Matrix(15,15,lambda a,b:c[a]*c[b]*Q.to_sympy(m.inner(pp[a],pp[b])))
    m.need(spectrum(H)=={'28':1,'10':5,'18':9},'Actual full holomorphic source multiplication Gram, no clock input')
    I=s.eye(15)
    P1=cl((H-10*I)*(H-18*I)/180)
    P5=cl((H-28*I)*(H-18*I)/144)
    P9=cl((H-28*I)*(H-10*I)/(-80))
    for P,rank in ((P1,1),(P5,5),(P9,9)):
        m.need(zero(P*P-P) and P==P.H and s.trace(P)==rank,'Canonical spectral projectors of actual multiplication')
    m.need(P1+P5+P9==I,'Complete intrinsic1+5+9 decomposition')
    J=s.zeros(25,15)
    for col,(i,j) in enumerate(ij):
        if i==j:J[5*i+j,col]=1
        else:J[5*i+j,col]=J[5*j+i,col]=1/s.sqrt(2)
    sing=s.Matrix([int(i==j) for i,j in ij])/s.sqrt(5)
    m.need(zero(P1-sing*sing.T),'Actual singlet is nondegenerate symmetric bilinear form')
    # Independent entire Wick-r2 source tensor, not its imported spectrum.
    ww=[]
    for i in range(5):
        for j in range(5):
            ww.append(m.add(*[(Q.one/factor(a),m.mul(multi_diff(ps[i],a),shift_to_negative(multi_diff(ps[j],a))))
                              for a in it.combinations_with_replacement(range(4),2)]))
    G2=s.Matrix(25,25,lambda a,b:Q.to_sympy(m.inner(ww[a],ww[b]))/s.sqrt(ns[a//5]*ns[a%5]*ns[b//5]*ns[b%5]))
    m.need(zero(H-10*I-J.T*G2*J),'Literal full15x15 relation between holomorphic product and conjugate pair tensor')
    # Solve intrinsic P9(v tensor v)=0 OVER C before importing original N.
    a,b,c0,d,e=s.symbols('a b c d e')
    vs=(a,b,c0,d,e);v=s.Matrix([a,s.sqrt(3)*b,s.sqrt(3)*c0,s.sqrt(3)*d,s.sqrt(6)*e])
    actual=list(P9*symvec(v,ij))
    quadrics=[a*a-d*d-4*e*e,a*b+c0*d,a*c0+b*d,a*d+b*c0,
              b*b-d*d,c0*c0-d*d,b*e,c0*e,d*e]
    CA=coeffmat(actual,vs);CQ=coeffmat(quadrics,vs)
    m.need(CA.rank()==CQ.rank()==CA.col_join(CQ).rank()==9,'Exact nine intrinsic quadratic relations, no real-form assumption')
    gb=s.groebner([p.subs(a,1) for p in quadrics],b,c0,d,e)
    solutions=s.solve_poly_system(list(gb),b,c0,d,e)
    m.need(len(solutions)==6,'All affine complex roots found in a=1 chart')
    for var in (b,c0,d,e):
        rest=[x for x in (b,c0,d,e) if x!=var]
        excluded=s.groebner([p.subs({a:0,var:1}) for p in quadrics],*rest)
        m.need(list(excluded)==[1],'No projective solutions missed on a=0 hyperplane')
    rays=[]
    for sol in solutions:
        raw=s.Matrix([1,s.sqrt(3)*sol[0],s.sqrt(3)*sol[1],s.sqrt(3)*sol[2],s.sqrt(6)*sol[3]])
        rays.append(cl(raw/s.sqrt((raw.H*raw)[0])))
    m.need(all(s.simplify(abs((x.H*y)[0])**2)==s.Rational(1,25) for i,x in enumerate(rays) for j,y in enumerate(rays) if i!=j),'Six derived COMPLEX projective rays are the exact regular simplex axes')
    oriented=[rays[0]]+[w if (rays[0].T*w)[0]<0 else -w for w in rays[1:]]
    W=s.Matrix.hstack(*oriented)
    m.need(zero(W*s.ones(6,1)) and zero(W.T*W-s.Rational(6,5)*s.eye(6)+s.ones(6)/5),'Intrinsic oriented simplex, only common sign/permutation remains')
    # Original target is introduced ONLY NOW, for literal unmarked comparison.
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48))
    match=[]
    for w in rays:
        found=[q for q in range(6) if zero(w*w.H-s.Rational(6,5)*N[:,q]*N[:,q].T)]
        m.need(len(found)==1,'Literal original simplex match after full complex reconstruction');match.append(found[0])
    # The 5-sector is not merely another dimension five: actual equivariance.
    uq=[symvec(w,ij)-sing/s.sqrt(5) for w in oriented]
    I5=cl(s.sqrt(5)*s.Matrix.hstack(*uq[:5])*W[:,:5].inv()/2)
    m.need(zero(I5.H*I5-s.eye(5)) and zero(H*I5-10*I5),'Isometric copy of actual native five in full source weight8')
    event_checks=0
    for p,q in it.combinations(range(6),2):
        perm=s.eye(6);perm.row_swap(p,q);T=cl(N*perm*N.T)
        m.need(zero(symrep(T,J)*I5-I5*T),'Actual S6 event action, not dimension matching')
        event_checks+=1
    # Clock-free all-edge source relations and relative effect family.
    edges=list(it.combinations(range(6),2));B=s.zeros(6,15)
    for k,(p,q) in enumerate(edges):B[p,k]=1;B[q,k]=-1
    Dall=cl(N*B)
    m.need(zero(Dall*Dall.T-6*s.eye(5)) and spectrum(B.T*B)=={'6':5,'0':10},'Intrinsic full edge Gram and cycle kernel')
    effects=[cl(Dall[:,k]*Dall[:,k].T/6) for k in range(15)]
    m.need(zero(sum(effects,s.zeros(5))-s.eye(5)) and all(zero(E*E-E/3) for E in effects),'All-edge effects are one-third rank1 projectors, not projectors')
    # All rooted relative bases, without any preselected clock/2subset.
    all_rank2=[];alpha=(s.sqrt(6)-1)/5;starchecks=0
    for q in range(6):
        slots=[j for j in range(6) if j!=q]
        U=cl(s.Matrix.hstack(*[N[:,j]-alpha*N[:,q] for j in slots]))
        m.need(zero(U.T*U-s.eye(5)),'Every marked star gives a complete relative ON basis')
        for i,j in it.combinations(range(5),2):all_rank2.append(cl(U[:,i]*U[:,i].T+U[:,j]*U[:,j].T))
        starchecks+=1
    # Compare equality in exact6D source coordinates, avoiding a fitted frame.
    keys={tuple(s.expand(x) for x in cl(N.T*P*N)) for P in all_rank2}
    m.need(len(keys)==60,'No clock chooses one of60 source-covariant rooted rank2 candidates')
    # Complete logical event-pair classification: no chosen original clock.
    perms=[]
    for aa,bb in edges:
        p0=list(range(6));p0[aa],p0[bb]=p0[bb],p0[aa];perms.append(tuple(p0))
    paircounts=Counter();cycles=set()
    for p0,p1 in it.product(perms,repeat=2):
        prod=tuple(p0[p1[k]] for k in range(6));reverse=tuple(p1[p0[k]] for k in range(6))
        moved=sum(k!=prod[k] for k in range(6))
        if moved==0:paircounts['identity']+=1
        elif prod==reverse:paircounts['double_transposition']+=1
        else:
            paircounts['three_cycle']+=1;cycles.add(prod)
            m.need(moved==3 and all(prod[prod[prod[k]]]==k for k in range(6)),'Every logical noncommuting event pair is a three-cycle')
    m.need(paircounts=={'identity':15,'double_transposition':90,'three_cycle':120} and len(cycles)==40,'Complete original logical clock-class census')
    clockmarks=set()
    for cyc in cycles:
        fixed=[k for k in range(6) if cyc[k]==k]
        for q in fixed:
            two=[k for k in fixed if k!=q]
            P=cl(sum([(N[:,j]-alpha*N[:,q])*(N[:,j]-alpha*N[:,q]).T for j in two],s.zeros(5)))
            clockmarks.add(tuple(s.expand(x) for x in cl(N.T*P*N)))
    m.need(clockmarks==keys,'Derived complete three-cycle/fixed-point family gives exactly same60 source polar flags')
    # Reference-state controls in the SAME original full source.
    stateGrams=[]
    for k in range(5):
        out=[m.mul(ps[k],f) for f in ps]
        stateGrams.append(s.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(out[i],out[j]))/(ns[k]*s.sqrt(ns[i]*ns[j]))))
    Bf0=s.diag(19,8,8,8,5);Bf4=s.diag(5,9,9,9,16)
    m.need(stateGrams[0]==Bf0 and stateGrams[4]==Bf4,'Actual normalized source states have different full regular response spectra')
    # Full marked q response in field coordinates and fixed-vacuum relative metric.
    Aq=m.add(*[(Q.from_sympy(raw[i,5]*[s.Rational(1,24),s.Rational(1,4),s.Rational(1,4),s.Rational(1,4),1][i]),ps[i]) for i in range(5)])
    out=[m.mul(Aq,f) for f in ps]
    Bq=s.Matrix(5,5,lambda i,j:s.Rational(3,5)*Q.to_sympy(m.inner(out[i],out[j]))/s.sqrt(ns[i]*ns[j]))
    wq=cl(s.sqrt(s.Rational(6,5))*N[:,5])
    m.need(zero(Bq-s.Rational(43,5)*s.eye(5)-5*wq*wq.T),'Marked source state selects1+4 spectral response, no clock required')
    # Actual event product, retaining all weight8 monomials.
    A0=m.add(*[(Q.from_sympy(raw[i,0]*[s.Rational(1,24),s.Rational(1,4),s.Rational(1,4),s.Rational(1,4),1][i]),ps[i]) for i in range(5)])
    delta=m.add((Q.one,A0),(-Q.one,Aq))
    m.need(m.inner(delta,delta)/2==Q(2) and m.inner(m.mul(delta,delta),m.mul(delta,delta))/4==Q(64),'Same-source actual event pair norm64 versus independent4')
    # Actual erasure derivation law on full local products.
    from check_conformal_erasure_carry import pd
    for i,j,aa in it.product(range(5),range(5),range(4)):
        left=pd(m.mul(ps[i],ps[j]),aa)
        right=m.add((Q.one,m.mul(pd(ps[i],aa),ps[j])),(Q.one,m.mul(ps[i],pd(ps[j],aa))))
        m.need(left==right,'Actual source erasure/product Leibniz identity, full grade7 output')
    # Independent copies vs correlated state with identical marginals.
    Bmean=(s.eye(5)+Bq)/2
    jointcorr=(s.eye(25)+s.kronecker_product(Bq,Bq))/2
    jointprod=s.kronecker_product(Bmean,Bmean)
    m.need(zero(jointcorr-jointprod-s.kronecker_product(Bq-s.eye(5),Bq-s.eye(5))/4),'Explicit original-state correlation contribution, no added interaction')
    m.need(jointcorr!=jointprod,'Identical marginals do not determine joint event geometry')
    # General invertible label-coordinate change WITH its transformed metric.
    D=cl(N[:,:5]-N[:,5]*s.ones(1,5));A=s.eye(5);A[0,1]=1
    metric=A.T*A;cut=s.diag(1,0,0,0,0);cutprime=A.inv()*cut*A;Tprime=D*A
    m.need(zero(Tprime*metric.inv()*Tprime.T-D*D.T),'Nonunitary coordinate covariance of total response with metric retained')
    m.need(zero(Tprime*cutprime*metric.inv()*Tprime.T-D*cut*D.T),'Nonunitary label covariance, no hidden resetting metric toI')
    # The exact30 rank1 and60rank2 configurations are not physical selectors.
    result={'status':'PASS_EXACT_CLOCK_FREE_INTRINSIC_EVENT_GEOMETRY','date':'2026-10-17','contract_id':'UR.NATIVE.JOINT_OPERATOR.01',
      'source_product':'Actual holomorphic local quartic multiplication Sym2(C5)->full weight8; no selected sigma/P2/real-state input',
      'clock_scope':'No selected sigma/Dreierclock or target event partition; accepted compiler-compatible complex structureJ and unitaryE8 source normalization remain premises.',
      'multiplication_gram_spectrum':spectrum(H),'multiplication_gram':m.mat(H),
      'crossing_gram_identity':'H_holo=10 I15+ H_r2_conjugate_restricted_to_Sym2; literal matrices, not spectrum only',
      'intrinsic_quadrics_coordinates':'v=(a,sqrt3 b,sqrt3 c,sqrt3 d,sqrt6 e)',
      'nine_complex_quadrics':[str(x) for x in quadrics],
      'complex_affine_chart_solutions_a1':[[str(t) for t in sol] for sol in solutions],
      'a_zero_projective_charts_excluded':4,'complex_projective_rays':6,'real_state_restriction':False,
      'six_ray_original_comparison_only_permutation':match,'target_polar_or_clock_input':False,
      'native_five_intertwiner_event_checks':event_checks,'new_native_five_source_weight':8,
      'all_edge_gram_spectrum':spectrum(B.T*B),'all_edge_kernel_dimension':10,
      'all_edge_effects':'15 rank1 effects E=P_root/3, sum I5 on response support; not identity on full H_E8',
      'logical_event_pair_census':dict(paircounts),'derived_logical_three_cycles':len(cycles),'derived_clock_fixed_mark_choices':120,'marked_star_ON_bases':starchecks,'unselected_distinct_rank2_candidates':len(keys),
      'state_f0_regular_response_spectrum':spectrum(Bf0),'state_f4_regular_response_spectrum':spectrum(Bf4),
      'marked_state_regular_response':'43/5 I+5 |wq><wq|, spectra68/5 and43/5^4',
      'same_source_actual_delta0_pair_squared_norm':64,'independent_delta0_pair_squared_norm':4,
      'actual_erasure_product_checks':100,'ordered_holomorphic_product_rank':15,'antisymmetric_label_kernel_dimension':10,
      'correlated_joint_actual_state_spectrum':spectrum(jointcorr),'independent_same_marginal_joint_spectrum':spectrum(jointprod),
      'general_coordinate_metric_correction':'PASS exact invertible nonunitary original-label shear',
      'reference_scope':'Central leading OPE coefficients are state-independent algebraic data for normalized states on same OPE algebra with defined coefficients. Full finite-response Grams are state-dependent. Positive finite-energy control states are not automatically admissible original P1 seeds.',
      'composition_scope':'Independent tensor copies/product states factorize polar geometry. Same-source common-profile holomorphic local products kill antisymmetric labels for every smooth common profile; explicit complete weight8 coefficients give15 rank. General fields/adjoints remain full OPE, not all histories classified.',
      'selection_scope':'Unique source-covariant geometry and1+5+9 spectral decomposition; six complex rays and actual native-five intertwiner. No selected marker, clock,2subset or physical3+2. Rank2 normalization inherits supplied label projector; no physical collapse/processing rule or P1 state/time selection.'}
    Path(__file__).with_name('intrinsic_event_geometry_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='multiplication_gram'},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
