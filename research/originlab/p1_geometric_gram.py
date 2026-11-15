"""Exact D4-symmetric Gram-rank source gate for geometric C5 plus F3.

Shows the ORIGINAL geometric reflection leaves only one mixed complex
two-point parameter. Positivity permits an 8-mode CAR source (|z|<1)
but ALSO a 6-mode null quotient (|z|=1); it does not create Q0.
"""
import sympy as sp


def require(ok,why):
    if not ok:
        raise ArithmeticError(why)


def controls():
    r5=sp.diag(1,sp.I,-1,-sp.I,1)
    rf=sp.diag(sp.I,-1,-sp.I)
    s5=-sp.Matrix(5,5,lambda i,j:int(i+j==4))
    sf=sp.Matrix(3,3,lambda i,j:int(i+j==2))
    xvars=sp.symbols("x0:15")
    unknown=sp.Matrix(5,3,xvars)
    equations=[a for a in r5.H*unknown*rf-unknown]
    equations += [a for a in s5*unknown*sf-unknown]
    linear,rhs=sp.linear_eq_to_matrix(equations,xvars)
    kernel=linear.nullspace()
    require(len(kernel)==1, "physical D4 cross-block space not dimension one")
    basic=sp.Matrix(5,3,lambda i,j:int(i==1 and j==0)-
                    int(i==3 and j==2))
    require(linear*sp.Matrix(list(basic))==sp.zeros(linear.rows,1),
            "canonical mixed source Gram not D4-equivariant")

    z=sp.symbols("z")
    M=z*basic
    G=sp.eye(5).row_join(M).col_join(M.H.row_join(sp.eye(3)))
    r8=sp.diag(r5,rf)
    s8=sp.diag(s5,sf)
    require(r8.H*G*r8==G and s8.H*G*s8==G,
            "full positive source matrix not geometric D4 invariant")
    require(sp.factor(G.det()-(1-z*sp.conjugate(z))**2)==0,
            "D4 source Gram determinant not the predicted one-parameter rank gate")

    real=sp.symbols("a",real=True)
    Greal=G.subs(z,real)
    lam=sp.symbols("lambda")
    require(sp.factor(Greal.charpoly(lam).as_expr()-
              (lam-1)**4*((lam-1)**2-real**2)**2)==0,
            "physical eight-mode Gram eigenvalues changed")
    witness={}
    for value in (0,sp.Rational(1,2),1,sp.Rational(3,2)):
        matrix=Greal.subs(real,value)
        eigen=matrix.eigenvals()
        witness[str(value)]={
            "rank":matrix.rank(),
            "positive_semidefinite":all(e>=0 for e in eigen),
            "positive_definite":all(e>0 for e in eigen),
        }
    require(witness=={
        "0":{"rank":8,"positive_semidefinite":True,"positive_definite":True},
        "1/2":{"rank":8,"positive_semidefinite":True,"positive_definite":True},
        "1":{"rank":6,"positive_semidefinite":True,"positive_definite":False},
        "3/2":{"rank":8,"positive_semidefinite":False,"positive_definite":False},
    }, "D4 Gram source positivity and null witness changed")
    # Exact Cholesky/Schur-complement normalization on the ENTIRE
    # physical-positive range |z|<1. This eliminates the apparent
    # unknown z from the D4 representation, not from source dynamics.
    h=sp.symbols("h",positive=True)
    Kinverse=sp.diag(1/h,1,1/h)
    T=sp.eye(5).row_join(-M*Kinverse).col_join(
        sp.zeros(3,5).row_join(Kinverse))
    defect=(T.H*G*T-sp.eye(8))
    structural=sp.zeros(8)
    structural[5,5]=-(h*h+z*sp.conjugate(z)-1)/(h*h)
    structural[7,7]=structural[5,5]
    require((defect-structural).applyfunc(sp.simplify)==sp.zeros(8)
            and T*r8==r8*T and T*s8==s8*T,
            "D4-equivariant physical Gram whitening failed")
    for signed_complex_z in (sp.Rational(1,2), (3+4*sp.I)/10):
        positive_h=sp.sqrt(sp.simplify(1-signed_complex_z*
                                      sp.conjugate(signed_complex_z)))
        Tsample=T.subs({z:signed_complex_z,h:positive_h})
        Gsample=G.subs(z,signed_complex_z)
        require((Tsample.H*Gsample*Tsample-sp.eye(8)).applyfunc(
                    sp.simplify)==sp.zeros(8),
                "D4 physical-positive CAR frame cannot be whitened")

    return {
      "status":"EXACT_D4_ORIGINAL_5PLUS3_GRAM_SOURCE_GATE",
      "mixed_complex_intertwiner_dimension":len(kernel),
      "normalized_mixed_matrix_nonzero_entries":[[1,0,"z"],[3,2,"-z"]],
      "reflection_forbidden_middle_grade_mixing":True,
      "rank_5_double_counting_impossible_under_given_D4_and_positive_marginals":True,
      "full_gram_determinant":"(1-|z|^2)^2",
      "full_gram_spectrum_for_real_z":"1 (x4), 1+|z| (x2), 1-|z| (x2)",
      "full_8_CAR_rank_condition":"|z|<1",
      "source_D4_equivariant_whitening_exists_for_every_strictly_positive_z":True,
      "whitening_T":"[[I5,-M*diag(1/h,1,1/h)],[0,diag(1/h,1,1/h)]], h=sqrt(1-|z|^2)",
      "whitening_preserves_geometric_C4_D4_operators":True,
      "source_parameter_z_does_not_affect_abstract_D4_CAR_character_within_positive_rank_8":True,
      "source_parameter_z_still_needed_for_weighted_Q0_and_time":True,
      "saturated_6_mode_null_quotient_condition":"|z|=1",
      "source_positivity_failure":"|z|>1",
      "regression_witnesses":witness,
      "actual_original_weighted_mixed_gram_z_computed":False,
      "source_actual_fermion_CAR_assignments_proved":False,
      "original_R113_Q0_and_time_computed":False,
      "F1_solved":False
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
