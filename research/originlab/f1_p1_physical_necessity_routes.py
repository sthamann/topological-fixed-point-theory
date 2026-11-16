"""Independent P1 origin tests against established physical constraints.

These tests are finite exact countermodels to common extrapolations:
(a) no proper half-sided modular inclusion with common cyclic/separating
    vacuum exists in finite-dimensional operator algebras;
(b) gauge-interface matching does not determine coherent cross-amplitudes;
(c) stationary four-clock constraint histories can yield relational time
    without determining a time arrow, modular time, or a unique history;
(d) repeated Zeno projection by itself does not activate dissipative time;
(e) four-branch elliptic geometry does not select the square modulus;
(f) a geometric bounce is an IR gravitational scenario, not a pregeometric
    derivation. None is the actual original weighted R113 Q0.
"""
from pathlib import Path
import json
import sympy as sp


def require(ok, msg):
    if not ok:
        raise ArithmeticError(msg)


def projective_origin_test():
    # Four branch points of a degree-two CP1 cover give g=(4-2)/2=1.
    # This is true for every noncolliding cross ratio, not just the square.
    a=sp.symbols("lambda")
    j=256*(1-a+a*a)**3/(a*a*(1-a)**2)
    square=sp.simplify(j.subs(a,sp.Rational(1,2)))
    hexagonal=sp.simplify(j.subs(a,(1+sp.sqrt(3)*sp.I)/2))
    require(square==1728 and hexagonal==0,
            "the two inequivalent elliptic enhanced-symmetry moduli changed")
    require((4-2)//2==1 and (6-2)//2==2,
            "Riemann-Hurwitz four-branch genus changed")
    return {
        "four_simple_branch_points_genus":1,
        "square_harmonic_cross_ratio_j":int(square),
        "hexagonal_equi_anharmonic_cross_ratio_j":int(hexagonal),
        "genus_one_and_modular_symmetry_uniquely_force_square":False
    }


def gauge_edge_mode_test():
    # A finite four-charge matching projector models a cut with the same
    # edge label on either side. NO continuum gravity is claimed.
    P=sp.diag(*[int(i//4 == i%4) for i in range(16)])
    a=sp.eye(16)[:,0]   # |0,0>
    b=sp.eye(16)[:,5]   # |1,1>
    plus=(a+b)/sp.sqrt(2)
    minus=(a-b)/sp.sqrt(2)
    flip=sp.Matrix([[0,1,0,0],[1,0,0,0],
                    [0,0,0,1],[0,0,1,0]])
    O=sp.kronecker_product(flip,flip)
    def reduced(psi, left=True):
        rho=psi*psi.H
        if left:
            return sp.Matrix(4,4,lambda i,j:
                    sum(rho[4*i+k,4*j+k] for k in range(4)))
        return sp.Matrix(4,4,lambda i,j:
                    sum(rho[4*k+i,4*k+j] for k in range(4)))
    require(P.rank()==4 and P*plus==plus and P*minus==minus,
            "physical boundary charge matching failed")
    require(reduced(plus)==reduced(minus)
            and reduced(plus,False)==reduced(minus,False),
            "edge-matched histories unexpectedly distinguish marginals")
    p=(plus.H*O*plus)[0]
    m=(minus.H*O*minus)[0]
    require(p==1 and m==-1, "coherent interfacial mixed mark changed")
    return {
        "edge_matching_projection_rank":4,
        "both_history_states_admissible":True,
        "both_one_side_reduced_states_identical":True,
        "cross_mark_values":["+1","-1"],
        "gauss_matching_fixes_original_coherent_amplitudes":False
    }


def stationary_relational_clock_test():
    S=sp.zeros(4)
    for t in range(4):
        S[(t+1)%4,t]=1
    psi=sp.zeros(16,1)
    for t in range(4):
        psi[5*t]=sp.Rational(1,2)
    G=sp.kronecker_product(S,S)
    rho=psi*psi.H
    clock=sp.Matrix(4,4,lambda i,j:
        sum(rho[4*i+k,4*j+k] for k in range(4)))
    require(G*psi==psi
            and len((G-sp.eye(16)).nullspace())==4,
            "period-four stationary source fixed space changed")
    require(clock==sp.eye(4)/4 and S**4==sp.eye(4)
            and S != sp.eye(4),
            "stationary relative clock must have trivial reduced modular flow")
    # Page-Wootters history also has a reversed solution: cycle direction
    # alone cannot select an irreversible arrow or absolute time scale.
    return {
        "four_clock_system_constraint":"(S tensor S) Psi=Psi",
        "invariant_global_state":True,
        "constraint_fixed_subspace_dimension":4,
        "clock_reduced_state":"I4/4",
        "reduced_clock_modular_hamiltonian":"scalar log(4)*I4",
        "conditionally_nontrivial_system_evolution":True,
        "unique_p1_or_arrow_from_stationarity":False,
    }



def joint_constraint_selection_test():
    # Conditional physical construction: match four edge sectors across two
    # sides and impose a single simultaneous cyclic *global* constraint.
    # Each separately admits rank four. Their intersection has rank ONE.
    S=sp.zeros(4)
    for j in range(4):
        S[(j+1)%4,j]=1
    G=sp.kronecker_product(S,S)
    match=sp.diag(*[int(i//4==i%4) for i in range(16)])
    invariant=sum((G**k for k in range(4)),sp.zeros(16))/4
    joint=match*invariant
    psi=sp.zeros(16,1)
    for j in range(4):
        psi[5*j]=sp.Rational(1,2)
    require(match.rank()==4 and invariant.rank()==4
            and match*invariant==invariant*match
            and joint.rank()==1 and joint*joint==joint
            and joint==psi*psi.H,
            "matched-edge plus cyclic group average must select unique Bell seam")
    reverse=sum((G**(-k) for k in range(4)),sp.zeros(16))/4
    require(invariant==reverse,
            "matched-edge consistency does not choose cyclic time orientation")
    return {
        "matching_projector_rank":4,
        "global_cycle_constraint_projector_rank":4,
        "commuting_constraint_intersection_rank":1,
        "unique_normalized_state_up_to_global_phase":
          "(|00>+|11>+|22>+|33>)/2",
        "relative_four_amplitudes_unique_under_assumed_constraints":True,
        "cyclic_orientation_selected":False,
        "original_six_D4_history_orbit_weight_selected":False,
        "physical_original_R113_projectors_identified":False
    }


def modular_finite_no_go():
    # Under the half-sided modular inclusion theorem both N and M
    # require the same cyclic + separating vector Omega. In finite
    # dimensions A -> A Omega is an injective surjection A -> H,
    # so dim(A)=dim(H). Proper inclusion N < M is impossible.
    #
    # Independent obstruction: if unitary conjugation scales a
    # finite-dimensional P>=0 as e^{-2 pi s} P, Hilbert-Schmidt norm
    # invariance forces P=0 for any s!=0.
    P=sp.diag(1,2,3,4)
    U=sp.diag(1,sp.I,-1,-sp.I)
    hs=sp.trace(P.H*P)
    require(sp.trace((U*P*U.H).H*(U*P*U.H))==hs
            and hs!=0 and sp.exp(-4*sp.pi)!=1,
            "finite Hilbert-Schmidt modular scaling obstruction failed")
    return {
        "finite_common_cyclic_separating_algebras_can_be_properly_nested":False,
        "finite_positive_nonzero_translation_satisfies_exact_modular_dilation":False,
        "standard_modular_inclusion_requires_nontrivial_infinite_source":True
    }


def gravitational_bounce_and_zeno_test():
    # Flat classical FRW GR: Hdot=-4piG(rho+p). With NEC, Hdot<=0,
    # so a contraction -> expansion bounce requires changing premises.
    G,rho,rc,w=sp.symbols("G rho rc w",positive=True)
    GR=-4*sp.pi*G*rho*(1+w)
    LQC=-4*sp.pi*G*rho*(1+w)*(1-2*rho/rc)
    require(sp.simplify(LQC.subs(rho,rc) + GR.subs(rho,rc))==0
            and sp.simplify(LQC.subs(rho,rc).subs({w:sp.Rational(1,3),
                                       G:1,rc:1}))==sp.Rational(16,3)*sp.pi,
            "effective bounce sign gate changed")
    # Rank-one Zeno projection with H=sigma_x:
    # P exp(-it H/n) P = cos(t/n) P,
    # [P exp(-it H/n) P]^n -> P.
    t=sp.symbols("t",real=True)
    n=sp.symbols("n",positive=True,integer=True)
    require(sp.limit(n*(sp.cos(t/n)-1),n,sp.oo)==0,
            "repeated rank-one Zeno projection limit changed")
    return {
        "flat_GR_bounce_under_NEC":False,
        "effective_LQC_bounce_sign_at_rho_crit":"Hdot=+4piG(rho+p)",
        "effective_bounce_is_pregeometric_P1_derivation":False,
        "repeated_rank_one_projection_generates_nontrivial_clock":False
    }


def original_source_gate():
    root=Path(__file__).resolve().parents[2]
    source=root/("experiments/theory-contracts/"
        "original-paper-source-reconciliation-20261102/"
        "contribution_20261106Q/domain/native-diagnostic-v1.json")
    j=json.loads(source.read_text(encoding="utf-8"))
    require(j["whole_source_bundle"] is None
            and j["coordinate_role"]=="spatial_tangential_circle_NOT_OS_time",
            "native R113 source gate changed: review physical source now")
    return {
        "native_R113_source_bundle":None,
        "diagnostic_circle_is_physical_OS_time":False,
        "fully_derived_p1_and_f1":False,
    }


def controls():
    return {
        "status":"INDEPENDENT_P1_ORIGIN_PHYSICS_CONSTRAINT_AUDIT",
        "finite_modular_translation":modular_finite_no_go(),
        "gauge_interface":gauge_edge_mode_test(),
        "relational_time":stationary_relational_clock_test(),
        "constraint_intersection":joint_constraint_selection_test(),
        "elliptic_modulus":projective_origin_test(),
        "bounce_and_zeno":gravitational_bounce_and_zeno_test(),
        "original_source":original_source_gate(),
        "F1_solved":False
    }


if __name__=="__main__":
    print(json.dumps(controls(),indent=2,sort_keys=True))
