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


def matrix_algebra_type_gate():
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
            len({(x+shift)%4 for x in source_phases})==3
            for shift in (0,2)),
    require(False, "unreachable")
