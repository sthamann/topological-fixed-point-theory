"""Ward test for Cartan preparations on the original three-leg affine word."""
from itertools import product

import sympy as sp

from affine_time import GLUE, HYPERCHARGES, TRIALITY, charge_signature, marked_charges, radial_word
from reconstruction import require, scalar


def root_matrix():
    return sp.Matrix.hstack(*(sp.Matrix(row["root"]).applyfunc(sp.Rational)
                               for row in marked_charges()))


def connected_response(roots, left, right, x, y):
    # [h(n), V_alpha(z)] = <h, alpha> z^n V_alpha(z).
    # Moving the left annihilator and the right creator to the vacuum
    # gives the bare word times the metric minus the Ward outer product.
    positions = sp.Matrix([1, sp.I*x, -x*y])
    positive = left.T*roots*positions
    negative = right.T*roots*positions.applyfunc(lambda value: 1/value)
    return -(positive*negative.T).applyfunc(sp.expand)


def matched_charge_words(original):
    gamma = original[:,1]
    expected = list(zip(GLUE,TRIALITY,HYPERCHARGES))
    words = []
    for signs in product((1,-1),repeat=8):
        if signs.count(-1)%2:
            continue
        right = sp.Matrix(signs)/2
        if gamma.dot(right)!=-1:
            continue
        legs = (-right-gamma,gamma,right)
        if [charge_signature(leg) for leg in legs]==expected:
            words.append(sp.Matrix.hstack(*legs))
    require(len(words)==2 and words[0]==original, "fixed middle root and charge-signature census changed")
    basis = sp.eye(8)
    delta = basis[:,6]+basis[:,7]
    reflection = basis-delta*delta.T
    require(reflection*original==words[1] and reflection*gamma==gamma,
            "the alternative word must be related by the marked root reflection")
    pair = basis[:,[1,6]]
    projections = sp.Matrix.hstack(*(pair.T*word[:,2] for word in words))
    require(pair.T*gamma==sp.zeros(2,1) and projections.det()==-sp.Rational(1,2),
            "the two charge signatures must distinguish independent gamma-orthogonal preparations")
    replies = [connected_response(word,pair,pair,sp.S.Half,sp.S.Half) for word in words]
    expected_offdiagonal = [-sp.Rational(25,16),sp.Rational(25,16)]
    require([reply[0,1] for reply in replies]==expected_offdiagonal,
            "the matched charge words must reverse the off-diagonal Ward response")
    return {"matching_words":len(words), "candidate_preparation_coordinates_1_based":[2,7],
            "right_roots":[list(map(str,word[:,2])) for word in words],
            "weyl_reflection_root_1_based":[7,8],
            "joint_projection_rank":projections.rank(),
            "connected_offdiagonal_at_half":list(map(str,expected_offdiagonal)),
            "scope":"The two words are Weyl-related. The sign contrast uses fixed polarization marks; "
                    "co-transporting them removes it. The second word requires a native R113 binding."}


def controls():
    roots = root_matrix()
    active = roots[:, :2]
    circle = sp.eye(8)
    old_pair = circle[:, 1:3]
    silent_pair = (circle[:, 1]-circle[:, 2]).row_join(circle[:, 3]-circle[:, 4])
    require(roots.rank()==2 and active.T*sp.Matrix.hstack(*roots.T.nullspace())==sp.zeros(2,6),
            "marked charge plane changed")
    gamma_orthogonal = sp.Matrix.hstack(*roots[:,1].T.nullspace())
    require((gamma_orthogonal.T*roots).rank()==1,
            "the restriction to gamma-orthogonal preparations must leave only one charge direction")
    require(old_pair[:, 0].dot(roots[:, 1])==0 and old_pair[:, 1].dot(roots[:, 1])==0
            and roots.T*(old_pair[:, 0]-old_pair[:, 1])==sp.zeros(3,1),
            "previous preparations must coincide on the marked charges")
    require(roots.T*silent_pair==sp.zeros(3,2), "orthogonal negative control changed")
    samples = [(x,y) for x in (sp.S.Half,sp.Rational(1,4)) for y in (sp.S.Half,sp.Rational(1,4))]
    active_responses = [connected_response(roots,active,active,x,y) for x,y in samples]
    time_response = sp.Matrix.hstack(*(response.reshape(4,1) for response in active_responses))
    require(time_response.det()==(66339-428652*sp.I)/256,
            "the marked time rectangle lost a charge-plane matrix direction")
    old = connected_response(roots,old_pair,old_pair,*samples[0])
    require(old==(-sp.Rational(25,16))*sp.ones(2), "unexpected old pair response")
    require(all(connected_response(roots,silent_pair,active,x,y)==sp.zeros(2)
                and connected_response(roots,active,silent_pair,x,y)==sp.zeros(2)
                for x,y in samples), "silent preparations failed the cross control")
    return {
        "status":"EXACT_ORIGINAL_AFFINE_POLARIZATION",
        "marked_charge_plane_rank":2, "silent_cartan_dimension":6,
        "gamma_orthogonal_charge_projection_rank":1,
        "single_time_pair_connected_rank":1, "sampled_charge_plane_matrix_span_rank":4,
        "time_rectangle_determinant":str(time_response.det()),
        "previous_pair_connected_at_half":str(old.tolist()),
        "matched_charge_words":matched_charge_words(roots),
        "negative_control":"Any polarization in the common orthogonal complement of all three roots.",
        "physical_source_identity_proved":False, "weighted_r113_values":None,
    }


def checked_polarization(data, source_id, time_id, contract):
    """Check marked endpoint forms against the same supplied bare-word normalization."""
    require(data.get("source_id")==source_id and data.get("time_id")==time_id and
            all(isinstance(data.get(key),str) and data[key].strip() for key in (
                "native_preparation_and_adjoint_binding", "physical_to_radial_time_binding",
                "absolute_reference_binding")),
            "polarization comparison needs the same source, time and absolute reference")
    original = root_matrix()
    roots = original
    if "right_root" in data:
        right_root = sp.Matrix([scalar(value) for value in data["right_root"]])
        require(right_root.rows==8, "the right root must have eight D8 coordinates")
        gamma = roots[:,1]
        require(all(value in (-sp.S.Half,sp.S.Half) for value in right_root) and
                sum(1 for value in right_root if value<0)%2==0 and gamma.dot(right_root)==-1,
                "the alternative word must use a root paired with the same middle current")
        legs = (-right_root-gamma,gamma,right_root)
        require([charge_signature(leg) for leg in legs]==list(zip(GLUE,TRIALITY,HYPERCHARGES)),
                "the alternative word must preserve the original discrete charge signature")
        roots = sp.Matrix.hstack(*legs)
    same_word = roots==original
    bases = [sp.Matrix([[scalar(value) for value in row] for row in data[key]])
             for key in ("left_polarizations","right_polarizations")]
    require(all(basis.rows==8 and basis.cols>0 and
                all(value.is_Rational for value in basis) and basis.rank()==basis.cols for basis in bases),
            "polarizations need independent real rational columns in the original D8 coordinates")
    left,right = bases
    seen = set()
    rows = []
    for sample in data["samples"]:
        x,y = (scalar(sample[key]) for key in ("left_radial_factor","right_radial_factor"))
        require(all(value.is_Rational and 0<value<1 for value in (x,y)) and (x,y) not in seen,
                "time samples need distinct rational pairs strictly between zero and one")
        seen.add((x,y))
        forms = sample["prepared_word_forms"]
        require(len(forms)==left.cols and all(len(row)==right.cols for row in forms),
                "prepared forms must cover both polarization bases")
        bare = sp.cancel(contract(sample["bare_word_form"]),extension=sp.I)
        expected_bare = sp.cancel(radial_word(x,y),extension=sp.I)
        matrix = sp.Matrix([[contract(form) for form in row] for row in forms])
        ward = connected_response(roots,left,right,x,y)
        residual = (matrix-bare*(left.T*right+ward)).applyfunc(sp.simplify)
        rows.append({"radial_factors":[str(x),str(y)], "absolute_bare_word":str(bare),
                     "affine_bare_word":str(expected_bare) if same_word else None,
                     "absolute_affine_bare_match":bare==expected_bare if same_word else None,
                     "normalization_degenerate":bare==0,
                     "prepared_residual_rank":residual.rank()})
    require(len(seen)>=2, "compare at least two ordered time pairs")
    distinguished = any(row["prepared_residual_rank"] for row in rows)
    degenerate = any(row["normalization_degenerate"] for row in rows)
    if distinguished:
        status = "DISTINGUISHED_SUPPLIED_POLARIZATION"
    elif degenerate:
        status = "INCONCLUSIVE_ZERO_BARE_WORD"
    else:
        status = "EXACT_SUPPLIED_WARD_SHAPE"
    return {"status":status, "samples":rows, "physical_source_identity_proved":False}
