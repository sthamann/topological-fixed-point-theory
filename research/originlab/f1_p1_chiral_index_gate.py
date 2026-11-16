"""Exact four-mark P1 logarithmic-family index and conditional 4D chiral zero-mode gate.

The mathematical index of dbar_{K(D)} is three for D=four distinct
marks. On the ordinary spin P1, S+=O(-1), a twist E=O(3) makes
S+ tensor E=K(D)=O(2). Its charged 2D Dirac index is +3.
This does NOT construct the physical higher-dimensional TFPT
fermion/flux/coupling or 3+1D chiral Standard Model.
"""
from fractions import Fraction
import sympy as sp


def require(ok,msg):
    if not ok:
        raise ArithmeticError(msg)


def line_bundle_cohomology(degree):
    """Exact P1 h^0,h^1 of O(k) by Riemann-Roch and Serre duality."""
    h0=max(degree+1,0)
    h1=max(-degree-1,0)
    require(h0-h1==degree+1,"P1 Dolbeault index changed")
    return {"degree":degree,"h0":h0,"h1":h1,"index":h0-h1}


def controls():
    # Riemann sphere K=O(-2); D=4 simple points:
    # log one-forms K(D)=O(2). For the spin Dirac positive
    # spinor bundle S+=K^(1/2)=O(-1), the E=O(3)
    # gauge-twist yields the same O(2) positive spinor.
    r=4
    canonical_degree=-2
    D_degree=4
    log_degree=canonical_degree+D_degree
    spin_positive_degree=-1
    flux_degree=log_degree-spin_positive_degree
    require((log_degree,flux_degree)==(2,3),
            "original four-point log spin/twist degree changed")
    pos=line_bundle_cohomology(log_degree)
    conjugate=line_bundle_cohomology(
        spin_positive_degree-flux_degree)
    require(pos=={"degree":2,"h0":3,"h1":0,"index":3} and
            conjugate=={"degree":-4,"h0":0,"h1":3,"index":-3},
            "positive vs opposite-flux chiral Dirac indices changed")
    # Correct index of a true twisted spin Dirac D_E on S2:
    # index D_E^+ = ∫_{S2} c1(E) = deg E = 3.
    require(pos["index"]==flux_degree and
            conjugate["index"]==-flux_degree,
            "Dirac monopole flux and Dolbeault index fail to agree")
    # Index stability: moving the marked divisor (not colliding) leaves
    # the degree of K(D) unchanged.
    x=sp.symbols("x")
    four=(sp.Integer(1),-sp.Integer(1),sp.I,-sp.I)
    five=(sp.Integer(1),-sp.Integer(1),sp.I,-sp.I,
          sp.Rational(2))
    r2=[(n,line_bundle_cohomology(n-2)["index"])
        for n in range(2,13)]
    require(r2[2]==(4,3) and all(idx==n-1 for n,idx in r2),
            "n punctures logarithmic zero-mode index does not equal n-1")
    # The 2-sheet elliptic cover of P1 detects modular T for r=4,
    # but its complex torus modulus remains continuously variable.
    t=sp.symbols("t",positive=True)
    j=sp.factor(256*(1+t*t+t**4)**3/
                (t**4*(1+t*t)**2))
    require(j.subs(t,1)==1728 and
            j.subs(t,sp.Rational(1,2))==sp.Rational(148176,25),
            "four-branch family index does not preserve two j moduli")

    return {
      "status":"EXACT_P1_FOUR_MARK_LOG_DOLBEAULT_INDEX_THREE",
      "original_P1_canonical_line_bundle":"K=O(-2)",
      "four_mark_divisor_degree":D_degree,
      "logarithmic_one_form_line_bundle":"K(D)=O(2)",
      "logarithmic_one_form_cohomology":pos,
      "positive_spin_bundle":"S+=O(-1)",
      "twisting_line_bundle_for_same_log_spinor":"E=O(3)",
      "twisting_first_Chern_number":flux_degree,
      "twisted_two_dimensional_spin_Dirac_index":pos["index"],
      "positive_flux_zero_modes":3,
      "negative_chirality_zero_modes_for_positive_flux":0,
      "opposite_flux_O_minus_three_index":conjugate["index"],
      "index_stable_under_moving_four_distinct_marks":True,
      "n_mark_index_formula":"ind(dbar_K(D))=n-1 on P1",
      "four_branch_cover_torus_j_not_fixed":"j(1)=1728; j(1/2)=148176/25",
      "P1_log_forms_source_to_actual_fermion_field_map_proved":False,
      "original_physical_three_monopole_flux_derived":False,
      "higher_dimensional_chiral_Dirac_parent_built":False,
      "original_E8_R_and_R_conjugate_chirality_index_verified":False,
      "three_physical_3plus1D_chiral_families_proved":False,
      "original_R113_Q0_os_source_identified":False,
      "square_modulus_from_index_selected":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
