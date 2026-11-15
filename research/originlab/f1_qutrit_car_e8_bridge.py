"""Exact constraints for joining the original 1+3+5 qutrit algebra to CAR8 and E8.

This is a scoped mathematical source-architecture audit. The eight
independent complex CAR fields and physical original Q0/R113 are NOT derived.
No finite Hilbert space represents a nonzero-level affine current algebra.
"""
from collections import Counter
from math import lcm
import sympy as sp
from affine_operators import ROOTS, ROOT_SET, dot, components
from p1_e8_geometric_lift import (
    GEOM8, GEOM_TORUS, NATIVE_GLUE_TORUS, SIMPLE_DOUBLED,
    WEYL_WORD, reflect, character)


def require(ok, why):
    if not ok:
        raise ArithmeticError(why)


def _car_action(mode, mask, creation):
    occupied=(mask >> mode) & 1
    if (creation and occupied) or (not creation and not occupied):
        return None
    sign=-1 if (mask & ((1<<mode)-1)).bit_count()%2 else 1
    return (sign,mask^(1<<mode))


def _car_pair(first,second,mask):
    right=_car_action(second[0],mask,second[1])
    if right is None:
        return {}
    left=_car_action(first[0],right[1],first[1])
    return {} if left is None else {left[1]:right[0]*left[0]}


def _sum_actions(x,y):
    z={}
    for action in (x,y):
        for k,v in action.items():
            z[k]=z.get(k,0)+v
    return {k:v for k,v in z.items() if v}


def verify_eight_mode_fock_and_qutrit_tensor_host():
    """Check every Fock basis mask and all 8x8 CAR anticommutators.

    A commuting qutrit action is defined on a SEPARATE tensor factor.
    It is not a physical TFPT field identification or E8_1 VOA.
    """
    checks=0
    for i in range(8):
        for j in range(8):
            for mask in range(256):
                aa=_sum_actions(
                    _car_pair((i,False),(j,False),mask),
                    _car_pair((j,False),(i,False),mask))
                mixed=_sum_actions(
                    _car_pair((i,False),(j,True),mask),
                    _car_pair((j,True),(i,False),mask))
                creates=_sum_actions(
                    _car_pair((i,True),(j,True),mask),
                    _car_pair((j,True),(i,True),mask))
                require(aa=={} and creates=={} and
                        mixed==({mask:1} if i==j else {}),
                        "actual 8-mode Fock CAR Clifford relations failed")
                checks+=3
    # Each E_ab acts on a three-state index, and CAR8 on the separate
    # 256-state occupation index. Verify the commutation independently
    # on all 768 product basis states and all 9 matrix units.
    commuting_checks=0
    for a in range(3):
        for b in range(3):
            for mode in range(8):
                for state in range(3):
                    for mask in range(256):
                        def Eab(label):
                            if label is None or label[1]!=b:
                                return None
                            return (label[0],a,label[2])
                        def car_operator(label):
                            if label is None:
                                return None
                            new=_car_action(mode,label[2],False)
                            return None if new is None else (
                                label[0]*new[0],label[1],new[1])
                        vector=(1,state,mask)
                        require(Eab(car_operator(vector))==
                                car_operator(Eab(vector)),
                                "separate qutrit and fermion tensor actions do not commute")
                        commuting_checks+=1
    require(checks==49152 and commuting_checks==55296,
            "the complete product basis action audit was not exhaustive")
    return checks,commuting_checks


def matrix_algebra_type_gate():
    car_checks, commuting_checks=verify_eight_mode_fock_and_qutrit_tensor_host()
    qutrit_dim=3
    modes=8
    fock=1 << modes
    carrier=lcm(qutrit_dim, fock)
    require(fock==256 and carrier==768 and
            carrier%qutrit_dim==carrier%fock==0 and
            fock%qutrit_dim!=0 and (fock//2)%qutrit_dim!=0,
            "unital finite matrix/CAR representation divisibility changed")

    # A compressed three-dimensional code P in the CAR Fock space is
    # NOT a unital *-homomorphism on the full CAR operator products:
    # C(A* A) != C(A*) C(A) for an operator leaving the code.
    p=sp.diag(1,1,1,0)
    a=sp.eye(4)[:,3]*sp.eye(4)[0,:]
    cap_lhs=p*a.T*a*p
    cap_rhs=(p*a.T*p)*(p*a*p)
    require(cap_lhs==sp.diag(1,0,0,0) and
            cap_rhs==sp.zeros(4),
            "finite code compression incorrectly preserves all CAR products")

    # Nonzero-level affine E8 current relation has a c-number
    # [H_1,H_-1]=k I. Trace of any finite-dimensional matrix
    # commutator is zero, while k*dim(H)>0 for k=1.
    x=sp.Matrix([[0,1,0],[0,0,0],[0,0,0]])
    y=x.T
    require(sp.trace(x*y-y*x)==0 and sp.trace(sp.eye(3))==3,
            "finite-dimensional affine central-trace obstruction changed")
    return {
      "status":"EXACT_UNITAL_M3_VERSUS_EIGHT_COMPLEX_CAR_AND_LEVEL_ONE_GATE",
      "original_quantum_qutrit_matrix_algebra":"M3(C) on a 3-dimensional state space",
      "eight_complex_CAR_algebra":"CAR(C8)=M256(C)",
      "irreducible_CAR8_Fock_dimension":fock,
      "exhaustive_CAR_anticommutator_basis_checks":car_checks,
      "exhaustive_commuting_qutrit_CAR8_basis_checks":commuting_checks,
      "unital_M3_inside_irreducible_CAR8_Fock_exists":False,
      "unital_M3_inside_each_CAR8_even_parity_128_block_exists":False,
      "minimum_common_finite_unital_representation_dimension":carrier,
      "explicit_minimal_commuting_host":"C3 tensor wedge(C8), dimension 768",
      "exact_full_product_algebra":"M3 tensor M256=M768",
      "code_rank_three_compression_preserves_full_CAR_products":False,
      "exact_code_compression_product_defect":"C(A* A)=|0><0|, C(A*)C(A)=0",
      "finite_dimension_nonzero_affine_Kac_Moody_level_possible":False,
      "affine_central_trace_no_go":"Tr[H1,H(-1)]=0 != k*dim(H) for k=1",
      "infinite_mode_extension_required_for_E8_level_one":True,
      "original_physical_qutrit_and_8_CAR_source_identified":False,
      "actual_E8_current_VOA_on_host_proved":False,
    }


def _root_at(i,j):
    return tuple(2*(int(k==i)-int(k==j)) for k in range(8))


def _weyl(root, inverse=False):
    word=reversed(WEYL_WORD) if inverse else WEYL_WORD
    for k in word:
        root=reflect(root,SIMPLE_DOUBLED[k-1])
    return root


def colour_branch():
    # Pinned native E8 color A2 is coordinate SU3 on axes 5,6,7
    # in the doubled D8 root chart, fully neutral for native C4.
    a,b=_root_at(5,6),_root_at(6,7)
    ab=tuple(x+y for x,y in zip(a,b))
    require(a in ROOT_SET and b in ROOT_SET and ab in ROOT_SET
            and dot(a,b)==-1,"native color root A2 chart changed")
    roots=(a,b,ab)
    qn=tuple(2*x for x in NATIVE_GLUE_TORUS)
    qg=tuple(2*x for x in GEOM_TORUS)
    require(all(character(z,qn)==0 for z in roots),
            "native pinned color SU3 is no longer quarter neutral")

    pre=tuple(_weyl(z,inverse=True) for z in roots)
    require(pre==(
        (0,0,0,0,-2,-2,0,0),
        (1,1,1,1,1,1,-1,-1),
        (1,1,1,1,-1,-1,-1,-1)
    ),"explicit 31-Weyl root inverse source preimages changed")
    require(all(character(z,qg)==0 for z in pre),
            "geometric E8 preimage not quarter neutral")
    def typename(z):
        return "bilinear_D8" if all(x%2==0 for x in z) else "spinor"
    require([typename(z) for z in pre]==
            ["bilinear_D8","spinor","spinor"],
            "geometric color source does not require E8 twisted spinors")

    commute=[r for r in ROOTS if dot(r,a)==0 and dot(r,b)==0]
    require(len(commute)==72 and components(commute)==["E6"],
            "native A2 color E6 centralizer not closed")
    labels=Counter((dot(r,a),dot(r,b)) for r in ROOTS)
    adjoint_labels={(2,-1),(-2,1),(-1,2),(1,-2),(1,1),(-1,-1)}
    triplets={(1,0),(-1,1),(0,-1),(-1,0),(1,-1),(0,1)}
    require(labels[(0,0)]==72 and
            all(labels[label]==1 for label in adjoint_labels) and
            all(labels[label]==27 for label in triplets) and
            sum(labels.values())==240 and len(labels)==13,
            "248=(8,1)+(1,78)+(3,27)+(3bar,27bar) root branching failed")

    # IMPORTANT: the grouped GEOM_TORUS ordering interleaves
    # C5 and F3: the three F3 axes are 3,5,7, NOT 5,6,7.
    source_a,source_b=_root_at(3,5),_root_at(5,7)
    source_ab=tuple(x+y for x,y in zip(source_a,source_b))
    rawA2=(source_a,source_b,source_ab)
    profile=Counter(character(z,qg) for z in rawA2)
    profile+=Counter(character(tuple(-x for x in z),qg) for z in rawA2)
    profile[0]+=2  # two color Cartan generators
    require([profile[k] for k in range(4)]==[2,2,2,2],
            "natural three F3 complex-fermion SU3 is not quarter balanced")
    # Any common scalar twist of F3 changes all fundamental phases by
    # the same amount but leaves their root differences unchanged.
    # Thus its SU3 adjoint cannot become grade-zero simply by multiplying
    # the entire triplet by the f2 orientation character.
    source_phases=[GEOM_TORUS[k] for k in (3,5,7)]
    require(source_phases==[1,2,3] and
            all(len({(x+shift)%4 for x in source_phases})==3
                for shift in (0,2)),
            "scalar orientation f2 twist unexpectedly neutralizes SU3 roots")
    # The native Weyl map preserves the C4 multiset if BOTH clock and
    # currents are conjugated. It cannot turn this balanced source SU3
    # into the fixed all-grade-zero native SU3 simultaneously.
    images=tuple(_weyl(z) for z in rawA2)
    require(all(character(images[k],qn)==character(rawA2[k],qg)
                for k in range(3)) and
            [typename(z) for z in images].count("spinor")==2,
            "geometric F3 current/clock simultaneous Weyl audit changed")

    # For a genuinely physical 3 complex chiral color fundamental plus
    # five spectators, the free-fermion SU3 current has level 1.
    # For 8 complex fermions in the SU3 adjoint it has level 6.
    h=(1,-1,0)
    trfund=sum(x*x for x in h)
    tradj=sum((h[i]-h[j])**2 for i in range(3)
              for j in range(3) if i!=j)
    require(trfund==2 and tradj==12 and tradj//trfund==6,
            "fundamental-versus-adjoint SU3 Dynkin index changed")
    return {
      "status":"EXACT_NATIVE_E8_COLOR_PREIMAGE_REQUIRES_SPINOR_CURRENTS",
      "native_color_simple_roots_doubled":[list(a),list(b)],
      "geometric_E8_inverse_weyl_color_simple_roots":[list(pre[0]),list(pre[1])],
      "geometric_E8_inverse_weyl_color_positive_roots":[list(z) for z in pre],
      "geometric_E8_inverse_weyl_colour_positive_types":
          ["bilinear_D8","spinor","spinor"],
      "all_six_colour_root_types":{"bilinear_D8":2,"spinor":4},
      "native_colour_E6_commutant_roots":72,
      "native_colour_E6_commutant_dimension_including_Cartan":78,
      "charged_E8_color_weights_each_multiplicity":27,
      "charged_E8_color_weight_count":6,
      "exact_E8_branching":"248=(8,1)+(1,78)+(3,27)+(3bar,27bar)",
      "source_GEOM_TORUS_in_interleaved_C5_and_F3_mode_order":list(GEOM_TORUS),
      "source_F3_complex_mode_axes":[3,5,7],
      "source_F3_complex_mode_quarter_phases":source_phases,
      "natural_F3_fundamental_SU3_quarter_current_profile":[2,2,2,2],
      "native_pinned_color_SU3_quarter_current_profile":[8,0,0,0],
      "scalar_f2_orientation_character_removes_nontrivial_colour_root_grades":False,
      "natural_F3_SU3_recovered_as_native_colour_by_simultaneous_Weyl_conjugacy":False,
      "index_of_3_fundamental_plus_5_singlets_SU3_complex_fermions":1,
      "index_of_8_adjoint_SU3_complex_fermions":6,
      "physical_A2_color_fermion_source_identified":False,
      "native_original_colour_Weyl_intertwiner_sourced":False,
    }


def continuum_time_gate():
    # The exact one-step qutrit CP channel from f1_functorial_source
    # has q fixed. For microscopic a -> 0 and macroscopic t>0,
    # q^(t/a) -> 0 instantly. The limit is discontinuous at t=0.
    # A nontrivial strongly continuous limit requires q(a)=1-O(a).
    a,t,rate=sp.symbols("a t rate",positive=True)
    fixed=sp.Rational(1,2)
    require(sp.limit(fixed**(t/a),a,0,dir="+")==0,
            "scale-independent CP depolarization did not collapse")
    source_q=sp.exp(-rate*a)
    recovered=sp.simplify(source_q**(t/a))
    require(sp.simplify(recovered-sp.exp(-rate*t))==0 and
            sp.limit(source_q,a,0,dir="+")==1,
            "critical microscopic q(a) scaling failed")
    q=sp.symbols("q")
    require(sp.solve(q*q-q,q)==[0,1],
            "fixed exact blocked one-step kernel unexpectedly nontrivial")
    return {
      "status":"EXACT_CONTINUUM_CLOCK_SCALING_NECESSITY",
      "fixed_microscopic_q_half_at_any_lattice_spacing":"q(a)=1/2",
      "fixed_q_at_any_positive_time_continuum_limit":"0 (instant depolarization)",
      "fixed_q_continuum_lim_t_to_zero_of_Tt":"not identity",
      "strongly_continuous_nontrivial_time_from_fixed_q":False,
      "microscopic_scaling_for_finite_continuum_rate":"q(a)=exp(-rate*a)",
      "emergent_continuum_transfer":"q(t)=exp(-rate*t)",
      "physical_time_rate_and_absolute_clock_derived":False,
      "rate_zero_gapless_limit_needed_for_critical_affine_current":"conditional; not sufficient",
      "emergent_3plus1D_spatial_locality_proved":False
    }


def controls():
    return {
      "unital_qutrit_CAR_source_embedding":matrix_algebra_type_gate(),
      "E8_native_color_and_E6_backtransport":colour_branch(),
      "continuum_microscopic_clock_gate":continuum_time_gate(),
      "physical_original_Q0_R113_produced":False,
      "physical_E8_level_one_VOA_produced":False,
      "physical_chiral_SM_family_count_produced":False,
      "physical_3plus1D_gravity_produced":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
