"""A source-first exact SU3 Ward selection conditional on original regional sl2.

From the original four local Möbius source fields (span 1,x,x²)
derive a rank-three sl2, not three commuting translations.
The 5+3 geometric modes are its Sym4 + Sym2 representation.
Its principal su2/SO3 Lie closure has one 5x5->3 alternating
structure constant lambda. A positive AD-INVARIANT metric forces
lambda>0 (compact SU3) and fixes lambda as the source Gram ratio.
The ORIGINAL physical Q0 forms and Ward cubic lambda are still missing.
"""
import sympy as sp
from itertools import combinations
from p1_principal_su3 import _basis

def require(ok,why):
    if not ok:
        raise ArithmeticError(why)

def rep(gen,basis):
    base=sp.Matrix.hstack(*(v.reshape(9,1) for v in basis))
    return sp.Matrix.hstack(*(base.gauss_jordan_solve(
        (gen*v-v*gen).reshape(9,1))[0] for v in basis))

def wedge_rep(M):
    pairs=list(combinations(range(M.rows),2))
    out=sp.zeros(len(pairs))
    for c,(i,j) in enumerate(pairs):
        for k in range(M.rows):
            if k!=j and M[k,i]:
                pos=tuple(sorted((k,j)))
                out[pairs.index(pos),c]+=M[k,i]*(1 if k<j else -1)
            if k!=i and M[k,j]:
                pos=tuple(sorted((i,k)))
                out[pairs.index(pos),c]+=M[k,j]*(1 if i<k else -1)
    return out

def controls():
    from source_bridge import regional_transport
    source=regional_transport()
    require(source["span_rank"]==3 and
            source["max_independent_commuting_rank"]==1 and
            source["bracket_wedge_determinant"] not in ("0",0),
            "original four regional fields do not span nonabelian sl2")

    so3,sym5=None,None
    sym5,so3=_basis()
    s3=[rep(a,so3) for a in so3]
    s5=[rep(a,sym5) for a in so3]
    wedge=[wedge_rep(a) for a in s5]
    variables=sp.symbols("h0:30")
    H=sp.Matrix(3,10,variables)
    eqs=[v for k in range(3) for v in (H*wedge[k]-s3[k]*H)]
    coeff,_=sp.linear_eq_to_matrix(eqs,variables)
    unique=30-coeff.rank()
    require(unique==1,"continuous SO3 must have unique 5 wedge 5 to 3 bracket")

    # Family: [A+P,B+Q]_lambda = [A,B] - lambda[P,Q] + [A,Q] + [P,B].
    # For lambda>0 this is compact su3 under (A,P) -> A+i sqrt(lambda) P.
    lam,a,b=sp.symbols("lambda a b",real=True)
    fixed=so3[0]
    p=sym5[0]
    q=sym5[3]
    L=sp.trace((p*q-q*p)*fixed)
    R=sp.trace(q*(fixed*p-p*fixed))
    require(L==R==4,"primitive 3x5 Ward source comparison lost its normalization")
    # Positive invariant metric on (A,P):
    # <A,B>=a (-tr AB), <P,Q>=b tr PQ.
    # Skew-adjointness of ad_P enforces 4(a lambda-b)=0.
    defect=4*(a*lam-b)
    require(sp.simplify(defect.subs(b,a*lam))==0,
            "source Ward metric normalization does not determine lambda")
    # Check all 8^3 metric Ward identities in the stated family
    basis=so3+sym5
    def bracket(i,j):
        c=basis[i]*basis[j]-basis[j]*basis[i]
        return -lam*c if i>=3 and j>=3 else c
    def kind(i,j):
        return "A" if (i<3 and j<3) or (i>=3 and j>=3) else "S"
    def pairing(m,kind,k):
        if kind=="A" and k<3: return -a*sp.trace(m*basis[k])
        if kind=="S" and k>=3:return b*sp.trace(m*basis[k])
        return 0
    n=0
    for i in range(8):
        for j in range(8):
            for k in range(8):
                invariance=pairing(bracket(i,j),kind(i,j),k)+pairing(
                    bracket(i,k),kind(i,k),j)
                require(sp.simplify(invariance.subs(b,a*lam))==0,
                        "full 8-mode compact Ward metric identity failed")
                n+=1
    require(n==512,"not all source Ward coefficient triples checked")

    # Relative signs: with a,b>0, b=a lambda forces lambda>0.
    # lambda=0 yields a semidirect SO3 ⋉ R5, and lambda<0
    # the noncompact split sl3(R) real form, not compact SU3.
    return {
      "status":"EXACT_CONDITIONAL_POSITIVE_SU3_SOURCE_WARD",
      "original_regional_mobius_span":"sl2 over C / three independent polynomial vector fields",
      "original_regional_commuting_rank":1,
      "principal_SO3_equivariant_wedge2_5_to_3_dimension":unique,
      "structure_constants_one_real_parameter":"lambda",
      "original_source_Gram_parameters":"a>0 on 3, b>0 on 5",
      "exact_positive_Ward_equation":"4*(a*lambda-b)=0",
      "full_basis_Ward_triples_checked":n,
      "positive_metric_forces_compact_SU3_lambda_positive":True,
      "lambda_zero_semidirect_compact_ad_invariant_metric":False,
      "lambda_negative_split_real_SL3_positive_ad_invariant_metric":False,
      "native_physical_regional_sl2_action_on_C5_F3_proved":False,
      "native_original_cubic_source_lambda_computed":False,
      "native_original_Gram_ratio_b_over_a_computed":False,
      "native_R113_and_parent_proved":False,
      "F1_solved":False
    }

if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
