"""Source-adjacent charged P1 modulus tests, all exact, NO native R113 values.

Classifies the four-mark modulus dependence of the *free* level-one E8
lattice root-word and a SPECULATIVE paired relative-polarization
determinant insertion. These are not the complete TFPT source.
"""
import sympy as sp
from affine_operators import ROOTS, ROOT_SET, MARKED_P2, MARKED_GAMMA, MARKED_P1, dot, add, neg


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def neutral_chord_exponents(charge_four):
    """Exponents of 2 sin(delta/2), 2 cos(delta/2), and 2, in that order."""
    require(len(charge_four) == 4, "exactly four Euclidean charge vectors needed")
    q = [sp.Matrix(row) for row in charge_four]
    require(all(v.rows==q[0].rows for v in q), "charges must share a positive Euclidean carrier")
    require(sum(q,sp.zeros(q[0].rows,1)) == sp.zeros(q[0].rows,1),
            "charge word is not neutral")
    a = sp.simplify(q[0].dot(q[1])+q[2].dot(q[3]))
    b = sp.simplify(q[1].dot(q[2])+q[3].dot(q[0]))
    opposite = sp.simplify(q[0].dot(q[2])+q[1].dot(q[3]))
    require(sp.simplify(a+b+(q[0]+q[2]).dot(q[0]+q[2]))==0,
            "neutrality must force A+B=-||q0+q2||^2")
    return a,b,opposite


def root_triple_no_max():
    """Exact root pairs, then a no-maximum theorem for positive boson metric."""
    p2,g,p1=MARKED_P2,MARKED_GAMMA,MARKED_P1
    charge=(p2,g,neg(p1))
    require(add(p2,g)==p1 and all(dot(q,q)==2 for q in charge) and
            all(dot(charge[i],charge[j])==-1
                for i in range(3) for j in range(i+1,3)),
            "the native E8 marked root triangle was altered")
    count=0
    by_root=[]
    for root in ROOTS:
        count_root=0
        for other in ROOTS:
            if dot(root,other)!=-1:
                continue
            third=neg(add(root,other))
            require(third in ROOT_SET and dot(third,third)==2 and
                    dot(root,third)==dot(other,third)==-1,
                    "non-A2 E8 root triangle found")
            count+=1
            count_root+=1
        by_root.append(count_root)
    require(len(ROOTS)==240 and count==13440 and
            set(by_root)=={56}, "E8 root-pair census is not 240 times 56")

    theta=sp.symbols("delta",real=True,positive=True)
    free_abs=1/(4*sp.sin(theta))
    require(sp.diff(free_abs,theta).subs(theta,sp.pi/2)==0 and
            sp.simplify(sp.diff(free_abs,theta,2).subs(theta,sp.pi/2))==sp.Rational(1,4),
            "free charged return not minimized at square")
    # This theorem applies to ANY fixed neutral real Euclidean charge
    # word: A+B=-||q0+q2||² <=0. In x=sin²(delta/2), the log
    # amplitude is A/2 log(x)+B/2 log(1-x). A stationary point needs
    # A and B to have the same sign. Both cannot be positive.
    # If they are negative, the log second derivative is positive:
    # every nonconstant interior stationary point is a strict minimum.
    return {
        "status":"EXACT_NATIVE_E8_FREE_VERTEX_NO_SQUARE_SELECTION",
        "charged_step_degrees_mod_12":[1,2,9],
        "original_marked_word":["p2","gamma","minus_p1","vacuum"],
        "e8_roots":240,
        "ordered_root_pairs_with_dot_minus_one":count,
        "partners_per_root":56,
        "root_triangle_pairings":[-1,-1,-1],
        "free_vertex_modulus":"1/(4*sin(delta))",
        "square_stationary_type":"strict minimum",
        "general_neutral_charge_theorem":"one fixed positive-metric neutral vertex word cannot have a strict interior maximum",
        "native_full_weighted_source_computed":False
    }


def paired_determinant_candidate():
    """Exact *conditional* compensation of the marked root-word poles.

    The model uses U_t from the proposed relative polarization, not the
    original Calderon/deck physical operator. It only specifies a
    source-binding experiment for the actual Pfaffian/eta factor.
    """
    t=sp.symbols("t",positive=True)
    a=(1-t*t)/(1+t*t)
    b=2*t/(1+t*t)
    R=sp.diag(1,-1)
    X=sp.Matrix([[0,1],[1,0]])
    U=(a*R+b*X)*R
    dm=sp.factor((sp.eye(2)-U).det())
    dp=sp.factor((sp.eye(2)+U).det())
    free_vertex=sp.factor((1+t*t)/(8*t))  # 1/(4 sin(delta)); sin(delta)=2t/(1+t²)
    dressed=sp.factor(dm*dp*free_vertex)
    expected=2*t/(1+t*t)
    require(sp.simplify(dressed-expected)==0 and
            sp.simplify(dm*dp-16*t*t/(1+t*t)**2)==0,
            "minimal paired-determinant source compensation failed")
    derivative=sp.factor(sp.diff(dressed,t))
    require(sp.simplify(derivative-2*(1-t*t)/(1+t*t)**2)==0 and
            sp.simplify(sp.diff(dressed,t,2).subs(t,1))==-1,
            "candidate square is not a strict maximum")
    withheld={str(p):str(dressed.subs(t,p)) for p in
              (sp.Rational(1,2),sp.Rational(1,3),sp.Rational(2,3),sp.S.One)}
    require(withheld=={"1/2":"4/5","1/3":"3/5","2/3":"12/13","1":"1"},
            "candidate withheld rational predictions changed")
    # If relative two-line determinant factors enter to a net exponent
    # m, the full modulus shape is (sin(delta))**(2*m-1).
    # A strict square maximum requires m > 1/2.
    return {
        "status":"EXACT_DECLARED_PAIRED_DETERMINANT_CANDIDATE",
        "free_E8_root_vertex":"1/(4*sin(delta))",
        "relative_even_determinant":"4/(1+t*t)",
        "relative_odd_determinant":"4*t*t/(1+t*t)",
        "paired_determinant":"4*sin(delta)**2",
        "product_times_root_word":"sin(delta)=2*t/(1+t*t)",
        "square_is_unique_maximum_with_this_pairing":True,
        "source_net_determinant_exponent_threshold":"m>1/2 in (det_even*det_odd)^m",
        "m_equal_half":"flat in delta",
        "m_equal_zero":"square strict minimum",
        "withheld_ratios_relative_to_square":withheld,
        "source_Pfaffian_and_ghost_exponent_derived":False,
        "original_Calderon_deck_operators_identified":False,
        "original_R113_values_computed":False
    }



def _pfaffian(antisymmetric):
    """Recursive exact Pfaffian for a small skew matrix (test only)."""
    n=antisymmetric.rows
    require(n==antisymmetric.cols and n%2==0 and
            antisymmetric == -antisymmetric.T,
            "Pfaffian needs an even skew symmetric matrix")
    if n==0:
        return sp.S.One
    return sp.expand(sum(
        (-1)**(j+1)*antisymmetric[0,j]*_pfaffian(
            antisymmetric.extract([k for k in range(1,n) if k!=j],
                                  [k for k in range(1,n) if k!=j]))
        for j in range(1,n)))


def two_face_fermion_sewing():
    """Four exact two-face CAR orientations with a relative unitary U.

    This is a *declared* quadratic fermion comparison, not the actual
    full TFPT Dirac/Pfaffian source or a proof of OS positivity.
    """
    t,g=sp.symbols("t g",positive=True)
    a=(1-t*t)/(1+t*t)
    b=2*t/(1+t*t)
    R=sp.diag(1,-1)
    X=sp.Matrix([[0,1],[1,0]])
    U=(a*R+b*X)*R
    identity=sp.eye(2)
    sin_delta=2*t/(1+t*t)
    forward=(1-g*g)**2+4*g*g*sin_delta**2
    antioriented=(1+g*g)**2-4*g*g*sin_delta**2
    expected={
      "same_forward_sign":forward,
      "same_opposite_sign":antioriented,
      "reverse_forward_sign":(1-g*g)**2,
      "reverse_opposite_sign":(1+g*g)**2
    }
    results={}
    for reverse in (False,True):
        for sign in (-1,1):
            row_top=identity.row_join(-g*U)
            second=U.inv() if reverse else U
            K=row_top.col_join((sign*g*second).row_join(identity))
            name=("reverse_" if reverse else "same_")+(
                "forward_sign" if sign==-1 else "opposite_sign")
            residue=sp.factor(K.det()-expected[name])
            require(residue==0,
                    "two-face determinant changed in the orientation/sign class")
            results[name]=str(sp.factor(K.det().subs(g,1)))
            if name=="same_forward_sign":
                for p in (sp.Rational(1,2),sp.S.One):
                    trial=K.subs({t:p,g:1})
                    block=sp.zeros(4).row_join(trial).col_join(
                        (-trial.T).row_join(sp.zeros(4)))
                    pf=_pfaffian(block)
                    require(sp.simplify(pf-trial.det())==0,
                            "eight-Majorana Grassmann Pfaffian failed its determinant identity")
    require(sp.simplify(forward.subs(g,1)-4*sin_delta**2)==0 and
            sp.simplify(antioriented.subs(g,1)-4*(1-sin_delta**2))==0 and
            sp.simplify(expected["reverse_forward_sign"].subs(g,1))==0,
            "unit-coupling orientation controls do not separate the shapes")
    # For g<1, the root-pole 1/(4sin delta) is NOT canceled:
    # W_g = (1-g²)²/(4sin delta)+g² sin delta.
    square_second=sp.factor((1-g*g)**2/4-g*g)
    require(sp.factor(square_second.subs(g,sp.sqrt(2)-1))==0 and
            sp.factor(square_second.subs(g,1))==-1,
            "finite sewing coupling has the wrong local-max threshold")
    # Holomorphic root vertex pole and full real paired-determinant
    # cancel to one oriented inverse-coordinate factor. The charged
    # cocycle remains an independent phase in the original E8 VOA.
    u=sp.symbols("u",nonzero=True)
    root_chiral=1/(2*(1-u*u))
    paired=(1-u*u)*(1-u**(-2))
    dressed=sp.factor(root_chiral*paired)
    on_unit_circle=dressed.subs(u,(1+sp.I*t)/(1-sp.I*t))
    require(sp.factor(dressed-(1-u**(-2))/2)==0 and
            sp.simplify(on_unit_circle*sp.conjugate(on_unit_circle)-
                4*t*t/(1+t*t)**2)==0,
            "oriented E8 vertex-pole cancellation failed")

    return {
      "status":"EXACT_FINITE_TWO_FACE_CAR_PFAFFIAN",
      "two_face_block":"K=[[I,-g U],[sign*g U_or_inverse,I]]",
      "unit_coupling_sewing_determinants":results,
      "same_orientation_forward_det":"(1-g*g)^2+4*g*g*sin(delta)^2",
      "reverse_orientation_forward_det":"(1-g*g)^2 (delta-independent)",
      "opposite_sign_same_orientation_det":"(1+g*g)^2-4*g*g*sin(delta)^2",
      "same_forward_unit_coupling_pfaffian":"Pf([[0,K],[-K.T,0]])=det(K)=4*sin(delta)^2",
      "oriented_chiral_pole_cancellation":"A_root(u)*Pf(K)=(1-u^(-2))/2 (up to root cocycle)",
      "fermionic_variables_in_this_control":8,
      "charged_root_dressed_at_g_one":"sin(delta), unique global maximum square",
      "g_below_one_pole_uncanceled":True,
      "local_maximum_square_threshold":"g>sqrt(2)-1 (does not give global maximum if g<1)",
      "physical_unit_coupling_and_CAR_orientation_proved":False,
      "actual_original_Dirac_eta_ghost_source_identified":False,
      "OS_reflection_positivity_proved":False,
      "original_R113_values_computed":False
    }



def original_r113_exterior_contract_gate():
    """Audit the ACTUAL R113 empty exterior state versus a NEW fermion trace.

    Owning text: experiments/theory-contracts/chat-triangulation-20261108/
    R113-ORIGINAL-UND-REFERENZ-20261109.txt, section 1.
    Q[b]=s[b] Gamma(R[b]); <Omega_F,Gamma(R)Omega_F>=1.
    Any determinant of a full Fock trace is therefore a DIFFERENT
    contraction unless independently bound inside original s[b].
    """
    t=sp.symbols("t",positive=True)
    a=(1-t*t)/(1+t*t)
    b=2*t/(1+t*t)
    U=(a*sp.diag(1,-1)+b*sp.Matrix([[0,1],[1,0]]))*sp.diag(1,-1)
    require(sp.factor(U.det()-1)==0, "relative finite word changed")
    # Lambda^0, Lambda^1, Lambda^2 of the two-mode carrier.
    # A full trace sums 0, 1, and 2 particle sectors; the original
    # R113 empty cap reads only the first 1x1 block.
    def gamma(r):
        return sp.Matrix([[1,0,0,0],
                          [0,r[0,0],r[0,1],0],
                          [0,r[1,0],r[1,1],0],
                          [0,0,0,r.det()]])
    vac=sp.Matrix([1,0,0,0])
    gu=gamma(U)
    gm=gamma(-U)
    require((vac.T*gu*vac)[0] == (vac.T*gm*vac)[0] == 1 and
            sp.factor(sp.trace(gu)-(sp.eye(2)+U).det())==0 and
            sp.factor(sp.trace(gm)-(sp.eye(2)-U).det())==0,
            "R113 exterior vacuum and full trace were conflated")
    # The very same two-mode operator U has an ordinary trace and a
    # parity-graded supertrace. This no longer needs a second -U
    # evolution, but STILL differs from the original R113 empty cap.
    parity=sp.diag(1,-1,-1,1)
    even_projector=(sp.eye(4)+parity)/2
    odd_projector=(sp.eye(4)-parity)/2
    even_trace=sp.factor(sp.trace(even_projector*gu))
    odd_trace=sp.factor(sp.trace(odd_projector*gu))
    supertrace=sp.factor(sp.trace(parity*gu))
    require(even_trace==2 and odd_trace==sp.factor(sp.trace(U)) and
            sp.factor(supertrace-(sp.eye(2)-U).det())==0 and
            sp.factor(sp.trace(gu)*supertrace-16*t*t/(1+t*t)**2)==0,
            "native parity grading failed to reproduce determinant pair")
    full_trace=sp.factor(sp.trace(gu)*sp.trace(gm))
    require(sp.factor(full_trace-16*t*t/(1+t*t)**2)==0 and
            (vac.T*gu*vac)[0]*(vac.T*gm*vac)[0]==1,
            "paired determinant cannot be inferred from empty Fock sources")
    return {
      "status":"NO_BARE_R113_VACUUM_DETERMINANT",
      "original_fermion_word":"Q[b] = s[b] Gamma(R[b])",
      "native_zero_fock":"<Omega_F,Gamma(R[b])Omega_F>=1",
      "proposed_full_trace":"Tr Gamma(U)=det(I+U)",
      "opposite_channel_full_trace":"Tr Gamma(-U)=det(I-U)",
      "same_operator_graded_trace":"Str((-1)^F Gamma(U))=det(I-U)",
      "same_operator_even_sector_trace":"2",
      "same_operator_odd_sector_trace":"2*(1-t*t)/(1+t*t)",
      "ordinary_times_supertrace":"16*t*t/(1+t*t)^2",
      "paired_full_trace":"16*t*t/(1+t*t)^2",
      "paired_empty_contraction":"1",
      "source_scalar_s_and_selection_could_contain_pair":True,
      "source_scalar_s_or_selector_proved_to_contain_pair":False,
      "direct_pfaffian_substitution_in_native_empty_Fock_permitted":False,
      "original_full_R113_produced":False
    }


def controls():
    return {
        "status":"EXACT_P1_CHARGED_MODULUS_CONTROLS",
        "free_charge_triangle":root_triple_no_max(),
        "paired_relative_determinant":paired_determinant_candidate(),
        "two_face_fermionic_sewing":two_face_fermion_sewing(),
        "original_R113_exterior_contract_gate":original_r113_exterior_contract_gate(),
        "physical_source_identity_proved":False,
        "F1_solved":False
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
