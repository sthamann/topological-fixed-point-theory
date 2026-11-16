"""Necessary gravitational-anomaly/inflow gate for physical chiral E8_1.

This completes the STANDARD mathematical E8 invertible 2+1D bulk/1+1D
edge anomaly matching from the same native E8 lattice. It does NOT
derive that bulk, its framing/eta phase, the 2D E8 worldsheet, 3+1D
chirality or the original all-field TFPT Q0/R113 physical source.

Original QFT003 'I6_tot=0' checks 4D SM gauge anomaly data; that
does not, by itself, cancel the E8_1 1+1D gravitational I4.
"""
from fractions import Fraction
import sympy as sp
from f1_ramond_e8_grade_bridge import SIMPLE_DOUBLED


def require(condition,why):
    if not condition:
        raise ArithmeticError(why)


def standard_model_left_weyl_anomalies():
    # One SM family including right-handed neutrino, written as LEFT
    # Weyl fields (Q,u^c,d^c,L,e^c,nu^c).
    one=Fraction(1,1)
    multiplets=(
      (6,Fraction(1,6)),(3,Fraction(-2,3)),
      (3,Fraction(1,3)),(2,Fraction(-1,2)),
      (1,one),(1,Fraction(0)))
    grav=sum(m*y for m,y in multiplets)
    cubic=sum(m*y**3 for m,y in multiplets)
    su3sq=2*Fraction(1,6)+Fraction(-2,3)+Fraction(1,3)
    su2sq=3*Fraction(1,6)+Fraction(-1,2)
    su3cub=2-1-1
    doublets=3+1
    require((grav,cubic,su3sq,su2sq,su3cub)==(0,0,0,0,0)
            and doublets%2==0,
            "one ordinary SM family has not cancelled gauge anomalies")
    return {
       "SM_U1_mixed_gravitational_trace":"0",
       "SM_U1_cubed_trace":"0",
       "SM_SU3_squared_U1_coefficient":"0",
       "SM_SU2_squared_U1_coefficient":"0",
       "SM_SU3_cubic_coefficient":"0",
       "SM_SU2_Witten_number_of_doublets":doublets,
       "SM_I6_gauge_and_mixed_anomaly_zero":True,
       "this_proves_chiral_E8_1_gravitational_I4_zero":False,
    }


def controls():
    B=SIMPLE_DOUBLED
    require(len(B)==8 and all(len(b)==8 for b in B),
            "native E8 simple-root frame missing")
    matrix=sp.Matrix(8,8,lambda i,j:
             sum(x*y for x,y in zip(B[i],B[j]))//4)
    minors=[matrix[:k,:k].det() for k in range(1,9)]
    require(minors==[2,4,6,5,4,3,2,1]
            and matrix.is_positive_definite and matrix.det()==1
            and all(entry==2 for entry in matrix.diagonal())
            and all(v.q==1 for v in matrix.inv()),
            "E8 Cartan K-matrix not positive even unimodular")
    # Taking the mirror E8 current sector cancels c_-, but it also
    # admits a rank-eight primitive null sublattice. A doubled boundary
    # CAN be gapped by appropriate interactions; mere orientation
    # doubling does not protect the desired chiral E8 currents.
    doubled=sp.diag(matrix,-matrix)
    null=sp.eye(8).col_join(sp.eye(8))
    require(doubled.det()==1 and
            null.rank()==8 and
            null.T*doubled*null==sp.zeros(8)
            and null[:8,:]==sp.eye(8),
            "doubled E8 opposite chirality null gapping lattice changed")
    c_left=8
    c_right=0
    T=sp.exp(-2*sp.pi*sp.I*sp.Rational(c_left-c_right,24))
    inverse=sp.conjugate(T)
    require(sp.simplify(T**3-1)==0 and
            sp.simplify(T-1)!=0 and
            sp.simplify(T*inverse-1)==0,
            "chiral E8 modular T/paired cancellation changed")

    # Tensor product with the opposite chirality:
    # modular T phase cancels, signed chiral charge c_- = 0.
    # Three SAME-handed E8 layers instead have a trivial scalar T
    # multiplier but still retain perturbative c_-=24 anomaly.
    require(sp.simplify(T**3-1)==0 and c_left*3!=0,
            "triple modular multiplier conflated with local anomaly")
    I4_magnitude=Fraction(c_left-c_right,24)
    require(I4_magnitude==Fraction(1,3),
            "local 1+1D gravitational anomaly coefficient changed")
    # Exact Riemann-Hurwitz for a degree-two cover of genus-zero P1:
    # 2g-2=2*(-2)+r, g=(r-2)/2 for r simple branch points.
    # The MINIMAL positive genus is g=1, r=4. This is conditional
    # on physical source selecting a genus-one double cover as the
    # minimal modular-anomaly diagnostic worldsheet.
    covers={r:(r-2)//2 for r in range(2,14,2)}
    require(covers=={2:0,4:1,6:2,8:3,10:4,12:5},
            "Riemann Hurwitz two-sheet branch formula changed")
    t=sp.symbols("t",positive=True)
    j=sp.factor(256*(1+t*t+t**4)**3/(t**4*(1+t*t)**2))
    require(j.subs(t,1)==1728 and j.subs(t,sp.Rational(1,2))==
            sp.Rational(148176,25),
            "four branch points do not retain distinct elliptic moduli")
    anomaly_sm=standard_model_left_weyl_anomalies()
    return {
       "status":"EXACT_CHIRAL_E8_1_GRAVITATIONAL_ANOMALY_AND_INFLOW_GATE",
       "chiral_E8_lattice_VOA_cL":8,
       "chiral_E8_lattice_VOA_cR":0,
       "chiral_c_minus":c_left-c_right,
       "E8_1_torus_character":"chi=E4/eta^8",
       "E8_1_modular_S_phase":"1",
       "E8_1_modular_T_phase":"exp(-2*pi*i/3)",
       "E8_1_modular_T_is_nontrivial":True,
       "local_two_dimensional_gravitational_I4_coefficient_magnitude":"p1/3",
       "original_TFPT_I6_gauge_anomaly_cancels_I4_grav_by_itself":False,
       "opposite_chirality_E8_double_partition":"abs(chi_E8)^2",
       "opposite_chirality_E8_double_modular_T_phase":"1",
       "opposite_chirality_E8_double_chiral_c_minus":0,
       "opposite_chirality_doubled_K_determinant":int(doubled.det()),
       "opposite_chirality_doubled_K_signature":0,
       "opposite_chirality_doubled_K_has_8_mutually_null_gapping_vectors":True,
       "opposite_chirality_doubled_chiral_edge_protected_by_anomaly":False,
       "three_same_chirality_E8_central_charge":24,
       "three_same_chirality_E8_modular_T_phase":"1",
       "three_same_chirality_E8_local_gravitational_anomaly_vanishes":False,
       "minimum_nontrivial_even_positive_unimodular_rank_by_known_lattice_classification":8,
       "unique_even_positive_unimodular_rank_eight_lattice_by_known_classification":"E8",
       "bosonic_chiral_invertible_minimality_axioms_from_original_Q0_derived":False,
       "riemann_hurwitz_degree_two_cover":"2g-2=-4+r, hence g=(r-2)/2",
       "simple_branch_points_to_genus":{str(r):g for r,g in covers.items()},
       "minimal_genus_one_double_cover_branch_marks":4,
       "modular_torus_cycle_required_for_diagnostic_not_chiral_CFT_existence":True,
       "physical_double_cover_minimality_axiom_derived":False,
       "four_marks_do_not_select_square_j":"j(t=1)=1728, j(t=1/2)=148176/25",
       "square_modulus_from_E8_anomaly_inflow_selected":False,
       "minimal_known_invertible_bulk":"oriented 2+1D E8 abelian Chern-Simons K theory",
       "bulk_E8_K_Cartan_matrix":matrix.tolist(),
       "bulk_K_principal_minors":[int(n) for n in minors],
       "bulk_K_determinant":int(matrix.det()),
       "bulk_K_rank":matrix.rank(),
       "bulk_K_signature":8,
       "bulk_K_even_integral_unimodular":True,
       "bulk_anyonic_ground_state_deg_on_torus":1,
       "bulk_boundary_chiral_c_minus":8,
       "bulk_boundary_matching_gravitational_anomaly_inflow_known":True,
       "source_required_inverse_T_holonomy":"exp(+2*pi*i/3)",
       "source_actual_eta_Pfaffian_holonomy_measured":False,
       "original_2D_chiral_E8_worldsheet_identified":False,
       "original_2plus1D_E8_bulk_or_opposite_sector_derived":False,
       "four_dimensional_chiral_SM_index_from_this_double_proved":False,
       "actual_original_R113_absolute_all_field_Q0_computed":False,
       "ordinary_SM_one_family_check":anomaly_sm,
       "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
