"""P1 physical-origin experiment: local gauge constraints and twisted seams.

Starts with the independently established Z2 toric-code/double-semion
topological source. This is CONDITIONAL on choosing a 2+1D square
lattice gauge model, NOT the original pregeometric TFPT Q0.

Exact tests:
- four-edge commuting Gauss/flux constraints retain four states;
- the nontrivial Z2 3-cocycle cannot be trivialized on full K=Z2
  by a normalized 2-cochain (boundary Frobenius obstruction);
- a trivial K=1 boundary restriction is unobstructed;
- double-semion four topological sectors have spins (1,i,-i,1);
- E8 invertible stack shifts central charge by 8 without changing
  the four anyons; only the absolute framing/vacuum phase distinguishes.
"""
import json
from itertools import product
import sympy as sp


def require(ok, msg):
    if not ok:
        raise ArithmeticError(msg)


def pauli_gate():
    X=sp.Matrix([[0,1],[1,0]])
    Z=sp.diag(1,-1)
    I=sp.eye(2)
    def tensor(*ops):
        out=sp.eye(1)
        for op in ops:
            out=sp.kronecker_product(out,op)
        return out
    star=tensor(X,X,X,X)
    flux_even=tensor(Z,Z,I,I)  # two shared edges, local overlap model
    flux_odd=tensor(Z,I,I,I)   # one shared edge, forbidden overlap
    require(star*flux_even==flux_even*star
            and star*flux_odd==-flux_odd*star
            and star*star==sp.eye(16)
            and flux_even*flux_even==sp.eye(16),
            "local star/flux commutators no longer obey intersection parity")
    p=(sp.eye(16)+star)*(sp.eye(16)+flux_even)/4
    require(p*p==p and p.rank()==4,
            "two local commuting gauge checks do not select fourfold sector")
    return {
        "four_edge_star_Gauss_law_commutes_with_even_flux_overlap":True,
        "odd_overlap_anticommutes":True,
        "two_commuting_checks_joint_rank":4,
        "unique_P1_from_local_Gauss_plus_one_flux":False,
    }


def cocycle_gate():
    g=(0,1)
    omega=lambda a,b,c:(-1)**(a*b*c)
    for a,b,c,d in product(g,repeat=4):
        left=omega(b,c,d)*omega(a,(b+c)%2,d)*omega(a,b,c)
        right=omega((a+b)%2,c,d)*omega(a,b,(c+d)%2)
        require(left==right,"nontrivial Z2 3-cocycle associator failed")
    require(all(omega(a,b,c)==1
                for a,b,c in product(g,repeat=3) if 0 in (a,b,c))
            and omega(1,1,1)==-1,
            "normalized nontrivial Z2 cocycle source altered")

    # A normalized 2-cochain beta has beta(0,a)=beta(a,0)=1.
    # Its arbitrary beta(1,1)=z cancels out from the coboundary
    # at (1,1,1), leaving delta beta=+1, not omega=-1.
    z=sp.symbols("z", nonzero=True)
    beta=lambda a,b: z if a==b==1 else sp.Integer(1)
    cob=lambda a,b,c:sp.cancel(beta(b,c)*beta(a,(b+c)%2)/
                                   (beta((a+b)%2,c)*beta(a,b)))
    require(all(cob(a,b,c)==1 for a,b,c in product(g,repeat=3)),
            "a normalized full-boundary 2-cochain unexpectedly cancels twist")
    return {
      "normalized_three_cocycle_source":"omega(a,b,c)=(-1)^(abc)",
      "Z2_cocycle_identity_checked":16,
      "omega_111":-1,
      "normalized_two_cochain_coboundary_111":"+1",
      "full_subgroup_K_Z2_cancels_omega_with_boundary_two_cochain":False,
      "trivial_subgroup_K_identity_restriction_unobstructed":True,
      "native_TFPT_boundary_uses_this_topological_G_or_K":False,
    }


def topological_spin_gate():
    Kds=sp.diag(2,-2)
    Ktc=sp.Matrix([[0,2],[2,0]])
    E8=2*sp.eye(8)
    for i,j in ((0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)):
        E8[i,j]=E8[j,i]=-1
    require(E8.det()==1 and all(E8[:i,:i].det()>0 for i in range(1,9))
            and Kds.det()==-4 and Ktc.det()==-4,
            "double-semion/Kitaev/E8 integral lattice source changed")
    labels=((0,0),(1,0),(0,1),(1,1))
    theta=[]
    for a,b in labels:
        v=sp.Matrix([a,b])
        q=(v.T*Kds.inv()*v)[0]
        theta.append(sp.simplify(sp.exp(sp.pi*sp.I*q)))
    require(theta==[1,sp.I,-sp.I,1],
            "double-semion anyon fourth roots changed")
    gauss=sum(theta)/2
    require(gauss==1,
            "Gauss-Milgram sum must detect chiral charge only modulo 8")
    composite=sp.diag(Kds,E8)
    require(composite.det()==-4 and E8.rank()==8,
            "E8 invertible stack added anyon sectors")
    # A rank-eight positive-definite unimodular even block has c=8,
    # but has zero discriminant group and no anyon labels.
    t8=sp.simplify(sp.exp(-2*sp.pi*sp.I*sp.Rational(8,24)))
    require(t8!=1 and sp.simplify(t8**3)==1,
            "gravitational framing phase cannot distinguish E8 stack")
    return {
      "double_semion_topological_spins":["1","i","-i","1"],
      "double_semion_four_anyon_types":4,
      "double_semion_chiral_c":0,
      "Gauss_Milgram_C_mod_eight":0,
      "E8_even_unimodular_rank":8,
      "E8_determinant":1,
      "double_semion_stacked_E8_anyon_count":4,
      "double_semion_stacked_E8_chiral_c":8,
      "anyon_spins_distinguish_E8_stack":False,
      "absolute_modular_T_vacuum_phase_for_c8":"exp(-2*pi*i/3)",
      "requires_original_framing_eta_or_absolute_reference":True,
    }


def controls():
    return {
      "status":"EXACT_CONDITIONAL_TOPOLOGICAL_P1_QUARTER_PHASE_GATE",
      "local_gauge_seam":pauli_gate(),
      "Z2_nontrivial_associator":cocycle_gate(),
      "double_semion_and_E8":topological_spin_gate(),
      "source_type":"preassigned 2+1D topological lattice gauge / Chern-Simons",
      "physical_P1_from_pregeometric_original_Q0":False,
      "original_R113_outer_Padm_and_all_fields_evaluated":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    print(json.dumps(controls(),indent=2,sort_keys=True))
