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


def controls():
    return {
        "free_charge_triangle":root_triple_no_max(),
        "paired_relative_determinant":paired_determinant_candidate(),
        "physical_source_identity_proved":False,
        "F1_solved":False
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
