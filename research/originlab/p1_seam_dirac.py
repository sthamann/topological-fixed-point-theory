"""A finite geometric seam Dirac operator on the original 5+3 typed modes.

This is an exact differential-complex / spectral-flow CANDIDATE from the
four-mark evaluation sequence. The actual original geometric Dirac
operator, source normalization, physical time, CAR field assignment,
selector and R113 values have not been identified.
"""
import sympy as sp


def require(ok,reason):
    if not ok:
        raise ArithmeticError(reason)


def controls():
    I=sp.Matrix(5,3,lambda row,col:int(row==col+1))
    R5=sp.diag(1,sp.I,-1,-sp.I,1)
    RF=sp.diag(sp.I,-1,-sp.I)
    S5=-sp.Matrix(5,5,lambda row,col:int(row+col==4))
    SF=sp.Matrix(3,3,lambda row,col:int(row+col==2))
    R=sp.diag(R5,RF)
    S=sp.diag(S5,SF)
    chirality=sp.diag(*([1]*5+[-1]*3))
    D=sp.zeros(5,5).row_join(I).col_join(
      I.T.row_join(sp.zeros(3,3)))
    m=sp.symbols("m",real=True)
    H=D+m*chirality

    require(I.rank()==3 and (I.T*I)==sp.eye(3),
            "geometric logarithmic inclusion lost injectivity")
    require(R5.det()==RF.det()==-1 and S5.det()==SF.det()==-1 and
            R.det()==S.det()==1,
            "five-plus-three determinant orientation fails SU8 lift")
    require(R.H*R==sp.eye(8) and S.H*S==sp.eye(8) and
            S*R*S==R.inv() and
            R*D==D*R and S*D*S==-D,
            "D4 chiral seam Dirac equivariance changed")
    require(D==D.H and D*chirality==-chirality*D and
            D.rank()==6 and D.nullspace()==[
                sp.eye(8)[:,0],sp.eye(8)[:,4]],
            "D4 Dirac complex does not have 3 paired and 2 endpoint modes")
    require(sp.factor(H.charpoly().as_expr())==
            sp.factor((sp.Symbol("lambda")-m)**2 *
                      (sp.Symbol("lambda")**2-(m*m+1))**3),
            "Dirac spectral flow polynomial changed")
    require(H.det()==-m*m*(1+m*m)**3 and
            S*H*S==-H.subs(m,-m) and R*H==H*R,
            "seam reflection and signed mass flow do not agree")

    e0=sp.eye(8)[:,0]
    einfty=sp.eye(8)[:,4]
    require(S*e0==-einfty and S*einfty==-e0 and
            (e0.H*H*e0)[0]==(einfty.H*H*einfty)[0]==m,
            "two fixed endpoints not exchanged by D4 reflection")
    # The full geometric index is +2. Restricting to just one
    # endpoint gives a single line crossing m=0. The one-sided
    # cap selection is NOT derived from the original physical source.
    return {
      "status":"EXACT_GEOMETRIC_D4_CHIRAL_COMPLEX_CANDIDATE",
      "original_five_plus_three_mode_rank":8,
      "canonical_inclusion_rank":3,
      "geometric_D4_rotation_determinants":["-1","-1","1"],
      "geometric_D4_reflection_determinants":["-1","-1","1"],
      "complex_D4_action_lifts_to_SU8":True,
      "Dirac_nonzero_modes":"3 positive and 3 negative for m=0",
      "Dirac_zero_modes":"two at geometric quarter fixed points 0 and infinity",
      "unreduced_geometric_index":2,
      "one_selected_endpoint_index":1,
      "Dirac_mass_family_spectrum":"m (x2), plus/minus sqrt(1+m^2) (x3 each)",
      "Dirac_mass_determinant":"-m^2(1+m^2)^3",
      "rotation_invariant":True,
      "reflection_odd_mass_intertwiner":"S H(m) S = -H(-m)",
      "single_endpoint_selection_from_original_P1_proved":False,
      "native_Dirac_selfadjoint_identified":False,
      "physical_nonnegative_time_generator_derived":False,
      "physically_independent_8_CAR_sources_proved":False,
      "full_R113_weighted_values_produced":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
