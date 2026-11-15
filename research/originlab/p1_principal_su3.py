"""Exact 5+3 Riemann-Roch principal SO3->SU3 Lie closure with D4 twist.

Two INEQUIVALENT physical models may share the original geometric C4
character: eight chiral complex fermion fields (conditional E8) OR an
eight-generator compact SU3 adjoint algebra (also conditional).
This tests the representations and Lie brackets, NOT the original Q0,
a P1 gauge-source field map, physical current OPEs or F1.
"""
from collections import Counter
from itertools import combinations, product
import sympy as sp

def require(ok,msg):
    if not ok:
        raise ArithmeticError(msg)

def _E(i,j):
    return sp.eye(3)[:,i]*sp.eye(3)[j,:]

def _basis():
    anti=[_E(0,1)-_E(1,0),_E(0,2)-_E(2,0),_E(1,2)-_E(2,1)]
    sym=[_E(0,1)+_E(1,0),_E(0,2)+_E(2,0),_E(1,2)+_E(2,1),
         sp.diag(1,-1,0),sp.diag(1,1,-2)]
    return sym,anti

def _rep(basis,rotation):
    f=sp.Matrix.hstack(*(sp.Matrix(x).reshape(9,1) for x in basis))
    return sp.Matrix.hstack(*(f.gauss_jordan_solve(
        (rotation*x*rotation.T).reshape(9,1))[0] for x in basis))

def _commutant_rank(actions):
    n=actions[0].rows
    vars=sp.symbols(f"a0:{n*n}")
    x=sp.Matrix(n,n,vars)
    equations=[entry for a in actions for entry in (x*a-a*x)]
    co,_=sp.linear_eq_to_matrix(equations,vars)
    return n*n-co.rank()

def _a2_subsystems():
    """Full native E8 root chart and its C4 fixed root census."""
    from affine_operators import ROOTS, charge_signature
    allroots=set(ROOTS)
    sub=set()
    for a in ROOTS:
        for b in ROOTS:
            if sum(x*y for x,y in zip(a,b))!=-4:
                continue
            c=tuple(x+y for x,y in zip(a,b))
            require(c in allroots,"E8 root sum missing")
            sub.add(frozenset((a,b,c,tuple(-x for x in a),
                               tuple(-x for x in b),tuple(-x for x in c))))
    fixed={r for r in ROOTS if charge_signature(r)[0]==0}
    nfixed=sum(s<=fixed for s in sub)
    require(len(sub)==1120 and nfixed==44 and len(fixed)==52,
            "E8 A2 gauge embeddings/quarter-fixed embeddings changed")
    return len(sub),nfixed

def controls():
    sym,anti=_basis()
    basis=[sp.I*v for v in sym]+anti
    flat=sp.Matrix.hstack(*(v.reshape(9,1) for v in basis))
    require(flat.rank()==8 and
            all(v.H==-v and sp.trace(v)==0 for v in basis),
            "eight compact su3 generators must be skew-Hermitian traceless")
    metric=sp.Matrix(8,8,lambda i,j:-sp.trace(basis[i]*basis[j]))
    require(metric.eigenvals()=={2:7,6:1},
            "positive ad-invariant compact su3 metric changed")

    commutators={}
    for i,j in combinations(range(8),2):
        c=basis[i]*basis[j]-basis[j]*basis[i]
        coeff=flat.gauss_jordan_solve(c.reshape(9,1))[0]
        require(all(sp.im(v)==0 for v in coeff),
                "su3 real structure constants are not real")
        commutators[(i,j)]=coeff
    require(sp.Matrix.hstack(*commutators.values()).rank()==8,
            "geometric 3+5 compact Lie algebra not perfect")
    require(all(commutators[i,j][:5,0]==sp.zeros(5,1)
                for i,j in combinations(range(5),2)) and
            all(commutators[i,j][5:8,0]==sp.zeros(3,1)
                for i in range(5) for j in range(5,8)),
            "the principal su2 bracket 5x5->3 or 3x5->5 is missing")
    # Matrix Jacobi identity is automatic, but check the 56 generator triples.
    for i,j,k in combinations(range(8),3):
        cycle=(
          basis[i]*(basis[j]*basis[k]-basis[k]*basis[j])
             -(basis[j]*basis[k]-basis[k]*basis[j])*basis[i]
          +basis[j]*(basis[k]*basis[i]-basis[i]*basis[k])
             -(basis[k]*basis[i]-basis[i]*basis[k])*basis[j]
          +basis[k]*(basis[i]*basis[j]-basis[j]*basis[i])
             -(basis[i]*basis[j]-basis[j]*basis[i])*basis[k]
        )
        require(cycle==sp.zeros(3),"principal su3 Jacobi failed")

    quarter=sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    inversion=sp.diag(1,-1,-1)
    rp=sp.diag(_rep(sym,quarter),_rep(anti,quarter))
    spm=sp.diag(_rep(sym,inversion),_rep(anti,inversion))
    rg5=sp.diag(1,sp.I,-1,-sp.I,1)
    sg5=-sp.Matrix(5,5,lambda i,j:int(i+j==4))
    rg3=sp.diag(sp.I,-1,-sp.I)
    sg3=sp.Matrix(3,3,lambda i,j:int(i+j==2))
    rg=sp.diag(rg5,rg3); sg=sp.diag(sg5,sg3)

    # Explicit D4 intertwiners between principal spin 2/1 and the
    # original Riemann-Roch degree 4/2 representations TWISTED by
    # the single order-two D4 character chi(r)=chi(s)=-1.
    t5=sp.Matrix([[-4,0,0,0,4],[0,8*sp.I,0,-8*sp.I,0],
                  [0,-8,0,-8,0],[4,0,0,0,4],[0,0,8,0,0]])
    t3=sp.Matrix([[0,-16,0],[-4,0,-4],[-4*sp.I,0,4*sp.I]])
    tr=sp.diag(t5,t3)
    require(tr.det()!=0 and
            tr*(-rg)==rp*tr and tr*(-sg)==spm*tr and
            rg**4==sp.eye(8) and sg**2==sp.eye(8) and
            sg*rg*sg==rg.inv(),
            "geometric D4 and principal SU3 characters not exactly sign twisted")
    require([_commutant_rank([_rep(sym,quarter),_rep(sym,inversion)]),
             _commutant_rank([_rep(anti,quarter),_rep(anti,inversion)])]==[4,2],
            "D4 commutants of spin 2 and spin 1 changed")

    # The original meromorphic central section f2=z²/(z⁴-1)
    # carries EXACTLY chi as a 1D line under geometric r,s.
    z=sp.symbols("z")
    f2=z*z/(z**4-1)
    require(rg5[2,2]==sg5[2,2]==-1 and
            sp.simplify(f2.subs(z,sp.I*z)+f2)==0 and
            sp.simplify(f2.subs(z,1/z)+f2)==0 and
            sp.limit(f2/z**2,z,0)==-1 and
            sp.limit(f2*z**2,z,sp.oo)==1,
            "original f2 cannot supply chi twist with double zero at both fixed points")

    all_a2,quarter_a2=_a2_subsystems()
    # Kostant's rank-two principal sl2 exponents from the standard
    # classification: A2=(1,2), B2=(1,3), G2=(1,5).
    # The two original spaces are (n-1,n+1) for n marks.
    rank_two={"A2":(1,2),"B2":(1,3),"G2":(1,5)}
    matches=[]
    for name,exps in rank_two.items():
        dimensions=(2*exps[0]+1,2*exps[1]+1)
        n=sum(dimensions)//2
        if dimensions==(n-1,n+1):
            matches.append((name,n))
    require(matches==[("A2",4)],"principal rank-two mark compatibility not unique")

    return {
      "status":"EXACT_RR_5PLUS3_PRINCIPAL_SU3_D4_TWIST",
      "compact_eight_generator_lie_algebra":"su3 (real) / sl3(C) (complex)",
      "full_commutator_pairs":len(commutators),
      "derived_lie_rank":sp.Matrix.hstack(*commutators.values()).rank(),
      "jacobi_generator_triples_checked":56,
      "positive_compact_metric_eigenvalues":{"2":7,"6":1},
      "principal_so3_spin1_dimension":3,
      "principal_so3_spin2_dimension":5,
      "geometric_five_intertwiner_determinant":str(t5.det()),
      "geometric_three_intertwiner_determinant":str(t3.det()),
      "compensator_character":"chi(r)=chi(s)=-1",
      "compensator_geometry":"f2(z)=z²/(z⁴-1), double zero at 0 and infinity",
      "untwisted_geometric_D4_is_an_SU3_bracket_automorphism":False,
      "twisted_geometric_D4_can_act_as_principal_SU3_automorphism":True,
      "principal_rank_two_simple_n_mark_matches":matches,
      "all_E8_A2_root_subalgebras":all_a2,
      "quarter_fixed_E8_A2_root_subalgebras":quarter_a2,
      "source_bracket_and_cubic_OPE_identified":False,
      "physical_colour_SU3_embedding_selected":False,
      "8_complex_free_fermion_interpretation_excluded":False,
      "original_weighted_Q0_and_R113_produced":False,
      "F1_solved":False,
    }

if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
