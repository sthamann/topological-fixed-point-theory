"""Resolve the charged two-face collar jet from joint covariance, time and charge data."""

from itertools import product
from pathlib import Path

import sympy as sp

from affine_operators import dot
from folded_time import SOURCES, marked_clock
from reconstruction import matrix, require, scalar
from verify import digest
from weak_clock import _root_transport
from weak_source import WEAK_ROOT, triangle_components


TANGENT = Path("experiments/theory-contracts/seam-collar-multifield-20261016/P1-INSERTION-TANGENT.txt")
INPUTS = (*SOURCES,TANGENT)


def _comm(a, b):
    return a*b-b*a


def _chart(energies, charges, base, scale):
    size = len(energies)
    require(size==len(charges) and size>1 and all(x.is_Rational for x in (*energies,*charges)),
            "supply a common rational energy and hypercharge diagram")
    require(type(scale) is int and scale>0 and base.is_Rational and 0<base<1 and
            all((scale*x).is_Integer for x in energies),
            "the radial step must give exact rational heat factors")
    require(all(a!=b or q==r for (a,q),(b,r) in product(zip(energies,charges),repeat=2)),
            "the declared charged bands must have disjoint energies")
    k = sp.diag(*energies)
    u = sp.diag(*(base**(scale*x) for x in energies))
    a = (sp.eye(size)+u**2).inv()
    n,s = sp.eye(size)-a,a*u
    p = a.row_join(s).col_join(s.row_join(n))
    require(p**2==p,"the original two-face formula must give a projection")
    return k,n,s,p,sp.diag(u,u),sp.diag(*charges,*charges)


def _divide_charged(value, diagonal):
    return sp.Matrix(value.rows,value.cols,
                     lambda i,j:value[i,j]/(diagonal[i]-diagonal[j]) if diagonal[i]!=diagonal[j] else 0)


def _resolve(chart, covariance_jet, charge_jet, heat_jet):
    k,n,_,p,u,y = chart
    size = k.rows
    require(all(m.shape==p.shape and m==m.H for m in (covariance_jet,charge_jet,heat_jet)),
            "supply three whole Hermitian jets in the same two-face order")
    # In a moving frame the commutator combination cancels the ambient
    # frame derivative. On nonzero hypercharge differences it is invertible.
    omega = _divide_charged(charge_jet,y.diagonal())
    intrinsic = _divide_charged(_comm(y,covariance_jet)-_comm(p,charge_jet),y.diagonal())
    rotation = _divide_charged(intrinsic[size:,size:],n.diagonal())
    field = _comm(k,rotation)
    diagonal_frame = sp.diag(rotation,rotation)
    recovered = (_comm(p,diagonal_frame+omega),_comm(y,omega),_comm(u,diagonal_frame+omega))
    checks = [expected==actual for expected,actual in zip(recovered,(covariance_jet,charge_jet,heat_jet))]
    if not all(checks):
        status = "EXCLUDED_IN_DECLARED_CHARGED_JET_CLASS"
    elif field.is_zero_matrix:
        status = "CARRIER_MOTION_SUFFICES_IN_DECLARED_CLASS"
    else:
        status = "INTRINSIC_CHARGED_VARIATION_REQUIRED_IN_DECLARED_CLASS"
    return {"status":status,
            "rank_of_intrinsic_field":field.rank() if all(checks) else None,
            "intrinsic_field":[list(map(str,row)) for row in field.tolist()] if all(checks) else None,
            "same_clock_motion_possible":None if not all(checks) else _comm(u,omega).is_zero_matrix,
            "replay":dict(zip(("covariance","hypercharge","heat"),checks)),
            "physical_source_identity_proved":False}


def from_input(data):
    chart = _chart(list(map(scalar,data["energies"])),list(map(scalar,data["hypercharges"])),
                   scalar(data["radial_base"]),data["exponent_scale"])
    return {**_resolve(chart,*(matrix(data[name]) for name in
                             ("two_face_covariance_tangent","two_face_hypercharge_tangent","two_face_heat_tangent"))),
            "scope":"exact charged first variation in the declared two-face collar with independently anchored hypercharge"}


def _controls_for_field(chart, field):
    k,n,s,p,u,y = chart
    size = k.rows
    rotation = _divide_charged(field,k.diagonal())
    diagonal_frame = sp.diag(rotation,rotation)
    dp,du = (_comm(operator,diagonal_frame) for operator in (p,u))
    dn = dp[size:,size:]
    z = sp.Matrix(size,size,lambda i,j:dn[i,j]/(s[i,i]+s[j,j]))
    zero = sp.zeros(size)
    cross_frame = zero.row_join(z).col_join((-z.H).row_join(zero))
    cross_heat = _comm(u,cross_frame)
    require(_comm(k,rotation)==field and rotation==-rotation.H and
            _comm(p,cross_frame)==dp and cross_heat[:size,size:]!=zero and du[:size,size:]==zero,
            "original collar must agree on the covariance and separate the two cross-face time responses")

    intrinsic = _resolve(chart,dp,sp.zeros(2*size),du)
    pure_motion = _resolve(chart,dp,_comm(y,diagonal_frame),du)
    mixed = _resolve(chart,dp+_comm(p,cross_frame),_comm(y,cross_frame),du+cross_heat)
    mismatch = _resolve(chart,dp,sp.zeros(2*size),cross_heat)
    require(intrinsic["intrinsic_field"]==mixed["intrinsic_field"]==
            [list(map(str,row)) for row in field.tolist()] and intrinsic["rank_of_intrinsic_field"]==2 and
            pure_motion["rank_of_intrinsic_field"]==0 and pure_motion["same_clock_motion_possible"] is False and
            mismatch["status"]=="EXCLUDED_IN_DECLARED_CHARGED_JET_CLASS",
            "the whole charged jet must distinguish an intrinsic insertion from a moving carrier")
    require(_comm(y,dp)==_comm(p,_comm(y,diagonal_frame)) and not _comm(y,dp).is_zero_matrix,
            "the independent hypercharge anchor must detect the transported field")
    return {"same_covariance_under_cross_face_motion":True,
            "same_covariance_and_heat_under_within_face_motion":True,
            "pure_motion_preserving_fixed_hypercharge_possible":False,
            "intrinsic_replayed":intrinsic["status"],"moving_replayed":pure_motion["status"],
            "mixed_replayed":mixed["status"],"independent_heat_mutant":mismatch["status"]}


def controls(root):
    folds = marked_clock(root)
    words = triangle_components()
    lifts,_ = _root_transport(words)
    old = sorted((right for _,_,right in words["spinor"]),key=lambda p:p.index(-1))
    weak = sorted({right for _,_,right in words["weak_spinor_minus"]},key=lambda p:-dot(p,WEAK_ROOT))
    require(len(old)==3 and len(weak)==2,"the original early collar needs three marked and two weak lines")
    cases = []
    for case in folds["cases"]:
        d = list(map(sp.Rational,case["folded_cartan"]))
        energies = [1+sum(x*y/2 for x,y in zip(d,r)) for r in old+weak]
        chart = _chart(energies,[sp.Rational(1,3)]*3+[sp.S.Half]*2,sp.S.Half,16)
        results = []
        for lift in lifts:
            colour = lift["colour_slot"]
            weak_index = 0 if lift["removed_weak_slot"]==5 else 1
            field = sp.zeros(5)
            field[3+weak_index,colour-1] = field[colour-1,3+weak_index] = lift["relative_phase"]
            results.append({"colour":colour,"weak_weight":dot(weak[weak_index],WEAK_ROOT),
                            **_controls_for_field(chart,field)})
        cases.append({"quarter":case["quarter"],"energies":list(map(str,energies)),"channels":results})
    return {"status":"EXACT_CHARGED_COLLAR_JET_SEARCH",
            "source_refs":{**folds["source_refs"],str(TANGENT):digest(root/TANGENT)},
            "cases":cases,"charged_components_checked":sum(len(case["channels"]) for case in cases),
            "next_original_observation":(
                "Bind the original two-face covariance, transported hypercharge and heat jets in one absolute frame. "
                "Keep the cross-face hypercharge blocks and common weighted empty reference."),
            "scope":"original folded energies and affine root maps in the specified finite two-face collar construction",
            "physical_source_identity_proved":False,"original_r113_values":None}
