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

    # New sharp **source-level** sufficient condition: the given
    # spectral triple has a grading Gamma. IF its actual OS/CAR source
    # realization makes this Gamma act by +1 on the five sections and
    # -1 on the three logarithmic differentials, and preserves the
    # positive Hermitian form, the mixed modulus z MUST vanish.
    # This identification with the original physical Gamma is not
    # assumed in any prior calculation.
    gamma=sp.diag(*([1]*5+[-1]*3))
    ward=gamma.H*G*gamma-G
    expected_ward=sp.zeros(5).row_join(-2*M).col_join(
        (-2*M.H).row_join(sp.zeros(3)))
    require(ward==expected_ward and
            (ward.subs(z,0)==sp.zeros(8)) and
            ward.subs(z,sp.Rational(1,2))!=sp.zeros(8),
            "chiral source-orthogonality Ward gate changed")
    # Even the additional demand that the candidate geometric
    # odd operator D be physically G-selfadjoint only fixes Im z=0.
    # An actual Gamma source isometry is strictly stronger: z=0.
    embedding=sp.Matrix(5,3,lambda i,j:int(i==j+1))
    candidateD=sp.zeros(5).row_join(embedding).col_join(
        embedding.T.row_join(sp.zeros(3)))
    d_ward=(candidateD.H*G-G*candidateD)
    real_z=sp.symbols("r",real=True)
    require(d_ward.subs(z,real_z)==sp.zeros(8) and
            d_ward.subs(z,(3+4*sp.I)/10)!=sp.zeros(8),
            "physical Hermiticity not isolated to a real overlap")
    require(gamma*r8==r8*gamma and gamma*s8==s8*gamma and
            gamma*candidateD==-candidateD*gamma,
            "original D4 and geometric form-parity candidate incompatible")

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
      "form_parity_ward_orthogonality_forces_z_zero":True,
      "physical_Gamma_must_act_plus_on_C5_minus_on_F3":"conditional source-map premise",
      "source_graded_eight_mode_positive_rank_when_Gamma_isometric":8,
      "candidate_D_physical_G_selfadjoint_only_for_real_z":True,
      "original_spectral_triple_Gamma_identified_with_form_degree":False,
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
