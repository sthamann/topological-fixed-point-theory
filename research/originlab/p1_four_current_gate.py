"""Exact E8 genuine four-current channel and conformal frame gates.

The four charges in both lanes are actual norm-two E8 roots.
These are chiral level-one lattice VOA correlators, not original
fully weighted R113 U/Y marked source forms.
"""
import sympy as sp
from affine_operators import MARKED_P2, MARKED_GAMMA, ROOT_SET, dot, neg


def require(ok, reason):
    if not ok:
        raise ArithmeticError(reason)


def free_vertex(z, charge):
    """Free lattice vacuum correlator without the constant cocycle phase."""
    require(len(z)==len(charge)==4, "need four actual field insertions")
    require(all(q in ROOT_SET for q in charge), "all four fields must be native E8 roots")
    require(all(sum(q[k] for q in charge)==0 for k in range(8)),
            "root currents must be charge neutral")
    return sp.factor(sp.prod(
        (z[i]-z[j])**dot(charge[i],charge[j])
        for i in range(4) for j in range(i+1,4)))


def controls():
    a,b=MARKED_P2,MARKED_GAMMA
    require(dot(a,a)==dot(b,b)==2 and dot(a,b)==-1, "native marked E8 roots changed")
    qA=(a,b,neg(a),neg(b))
    qB=(a,neg(b),neg(a),b)
    u=sp.symbols("u",nonzero=True)
    z=(sp.Integer(1),u,sp.Integer(-1),-u)
    GA,GB=free_vertex(z,qA),free_vertex(z,qB)
    x=sp.factor((z[0]-z[1])*(z[2]-z[3])/
                ((z[0]-z[2])*(z[1]-z[3])))
    ratio=sp.factor(GA/GB)
    require(sp.factor(ratio-((1-x)/x)**2)==0,
            "charged four-current crossing ratio changed")
    require(sp.factor(free_vertex((z[0],z[3],z[2],z[1]),qA)-GB)==0,
            "crossing must compare swapped coordinates, not two same-frame values")

    t=sp.symbols("t",positive=True)
    uparam=(1+sp.I*t)/(1-sp.I*t)
    require(sp.factor(x.subs(u,uparam)-t*t/(1+t*t))==0 and
            sp.factor(ratio.subs(u,uparam)-t**(-4))==0,
            "genuine four-current ratio not exact tangent fourth power")
    odd_split=(2/sp.pi)*sp.log(1/t)
    require(sp.simplify(sp.exp(2*sp.pi*odd_split)-t**(-4))==0,
            "charged ratio not tied to original v512 odd DtN splitting")
    require(ratio.subs(u,uparam).subs(t,sp.Rational(1,2))==16 and
            ratio.subs(u,uparam).subs(t,1)==1,
            "original nonsquare and square charged comparisons changed")

    # A single neutral three-current word has no intrinsic four-point
    # cross ratio: only a three-point conformal structure and an identity.
    z3=z[:3]
    def three_point(points):
        return sp.factor(1/sp.prod((points[i]-points[j])
                           for i in range(3) for j in range(i+1,3)))
    def mobius(v):
        return (2*v+3)/(v+2)
    def deriv(v):
        return 1/(v+2)**2
    require(sp.factor(three_point(tuple(mobius(v) for v in z3))*
                      sp.prod(deriv(v) for v in z3)-three_point(z3))==0,
            "native E8 weight-one three-point Möbius covariance failed")

    # The literal even-fermion projector does NOT multiply full trace
    # and supertrace. It yields the constant 2, not 4sin(delta)^2.
    aa=(1-t*t)/(1+t*t)
    bb=2*t/(1+t*t)
    U=sp.Matrix([[aa,-bb],[bb,aa]])
    gamma=sp.Matrix([[1,0,0,0],[0,U[0,0],U[0,1],0],
                     [0,U[1,0],U[1,1],0],[0,0,0,U.det()]])
    parity=sp.diag(1,-1,-1,1)
    even=(sp.eye(4)+parity)/2
    odd=(sp.eye(4)-parity)/2
    even_trace=sp.factor(sp.trace(even*gamma))
    odd_trace=sp.factor(sp.trace(odd*gamma))
    product=sp.factor(sp.trace(gamma)*sp.trace(parity*gamma))
    require(even_trace==2 and odd_trace==2*(1-t*t)/(1+t*t) and
            product==16*t*t/(1+t*t)**2,
            "fermionic even projection confused with trace/supertrace product")

    return {
        "status":"EXACT_FOUR_E8_CURRENT_CHANNEL_GATE",
        "native_E8_nonzero_currents_in_each_word":4,
        "four_current_channel_A":"(p2,gamma,-p2,-gamma)",
        "four_current_channel_B":"(p2,-gamma,-p2,gamma)",
        "cross_ratio":"x=sin(delta/2)^2=t^2/(1+t^2)",
        "free_channel_ratio":"G_A/G_B=((1-x)/x)^2=t^(-4)",
        "non_square_original_v512_t_half_channel_ratio":"16",
        "square_t_one_channel_ratio":"1",
        "original_v512_odd_split":"Delta_odd=(2/pi)*log(1/t)",
        "new_charged_neutral_relation":"log|G_A/G_B|=2*pi*Delta_odd",
        "crossing_relation":"G_B(z)=G_A(z0,z3,z2,z1), for all moduli",
        "crossing_symmetry_alone_selects_square":False,
        "equal_amplitudes_at_same_fixed_frame_would_select_square":True,
        "same_frame_channel_equality_from_original_source_proved":False,
        "three_root_plus_identity_is_only_three_point":True,
        "three_point_carries_fourth_mark_invariant":False,
        "even_fermion_projection_trace":"2 (all moduli)",
        "ordinary_times_supertrace":"16*t*t/(1+t*t)^2",
        "original_Padm_equals_literal_even_Fock_projector_proved":False,
        "original_R113_UY_sources_computed":False,
        "physical_3plus1_parent_derived":False,
        "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
