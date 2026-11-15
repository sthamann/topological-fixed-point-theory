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
    require(sp.factor(H.det())==-m*m*(1+m*m)**3 and
            S*H*S==-H.subs(m,-m) and R*H==H*R,
            "seam reflection and signed mass flow do not agree")

    e0=sp.eye(8)[:,0]
    einfty=sp.eye(8)[:,4]
    require(S*e0==-einfty and S*einfty==-e0 and
            (e0.H*H*e0)[0]==(einfty.H*H*einfty)[0]==m,
            "two fixed endpoints not exchanged by D4 reflection")
    # A finite "select one side, then run positive time" construction:
    # the two fixed endpoint states are the only zero modes of D^2.
    # Remove ONE endpoint at the cap (still a new physical assumption).
    # Both one-sided transfers have a unique vacuum, normalized gap 1,
    # and are exchanged by the original geometric reflection.
    p0=sp.eye(8)-einfty*einfty.H
    pinfty=sp.eye(8)-e0*e0.H
    h2=D*D
    require(p0*h2==h2*p0 and pinfty*h2==h2*pinfty and
            S*p0*S==pinfty and S*h2*S==h2,
            "one-sided finite time does not respect the seam double")
    keeper0=[0,1,2,3,5,6,7]
    keeper_infty=[1,2,3,4,5,6,7]
    h0=h2.extract(keeper0,keeper0)
    hinfty=h2.extract(keeper_infty,keeper_infty)
    require(h0.eigenvals()=={0:1,1:6} and
            hinfty.eigenvals()=={0:1,1:6} and
            sp.trace(h0)==sp.trace(hinfty)==6,
            "one-sided positive time lacks the expected unique vacuum")
    # On each one-sided Hilbert space T_t=P_vac+exp(-t)(I-P_vac).
    # The parameter t is only a normalized clock variable here:
    # physical elapsed time and absolute scale are NOT supplied.
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
      "one_sided_positive_time_H_squared_spectrum":["0 (x1)","1 (x6)"],
      "one_sided_transfer_formula":"T(t)=P_vac+exp(-t)*(I-P_vac)",
      "one_sided_vacuum_rank":1,
      "one_sided_two_branch_intertwinement":"S P_0 S=P_infinity",
      "one_sided_cap_selection_and_time_scale_derived":False,
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
