"""Original affine return anchors and optional same-source endpoint checks."""
from pathlib import Path
from math import gcd

import sympy as sp

from language import positive
from reconstruction import require
from original_paths import AFFINE
from sewing_dynamics import return_forms
from verify import digest

REGIONAL = Path("docs/research/fundamental-plan/TFPT-REGIONALE-KOMPOSITION-20261112.txt")
R102 = Path("experiments/theory-contracts/seam-positive-kernel-search-20261015/"
            "PRIMITIVE-INTEGRAL-COMPLETION-20261102.txt")
CALLER = Path("experiments/theory-contracts/original-paper-source-reconciliation-20261102/"
              "triangulation_20261106O/caller/ACTUAL-PRIMITIVE-REGULATOR-CALLER-TRACE.txt")
ORIGINAL_FILES = (AFFINE,REGIONAL,R102,CALLER)
ANCHOR_HASHES = {
    REGIONAL:"6888980fcb80f292c4da087a5bc3f486395097be914794cb22b6ce876f40c64d",
    R102:"695a3249781bc7f6c46e091a0faec0c8eda384a9cf64c6d67f501d5ddf10c0e6",
}


def affine_return(root):
    from original_paths import original_word
    original = original_word(root)
    refs = {str(path):digest(root/path) for path in ORIGINAL_FILES}
    require(all(refs[str(path)]==sha for path,sha in ANCHOR_HASHES.items()),
            "review changed original affine mode or charge derivation")
    # From the original level-one affine form on the vacuum/current pair.
    # The two dimensional matrices are consequences of the affine
    # commutators at these grades, not a truncation of all R113 fields.
    e = sp.Matrix([[0,0],[1,0]])
    f = e.H
    l0 = sp.diag(0,1)
    p = sp.diag(1,0)
    require((f*e-e*f)*p==p and l0*e-e*l0==e and
            not any(p*e*p) and p*f*(sp.eye(2)-p)*e*p==p,"original affine return identities failed")
    family_weight = sp.eye(4)[:,0]-sp.ones(4,1)/4
    carrier_weight = sp.ones(5,1)/2
    norm = carrier_weight.dot(carrier_weight)+family_weight.dot(family_weight)
    # Same marked D5+A3 root coordinates as the cited affine source.
    require(carrier_weight.dot(carrier_weight)==5*sp.Rational(1,4) and
            family_weight.dot(family_weight)==3*sp.Rational(1,4) and norm==2,
            "wrong positive primitive root pairing")
    q = sp.I
    family_clock = sp.diag(1,q)
    require(family_clock*e==q*e*family_clock,"original primitive phase failed")
    return {"status":"EXACT_ORIGINAL_AFFINE_RETURN",
            "source_refs":refs,
            "source_identity":"E8 level-one chiral affine vacuum, original circle radial grading",
            "two_vectors":"vacuum and a normalized class-one root current at grade one",
            "direct_vacuum_insertion":"0",
            "adjoint_ordered_return":"1",
            "quarter_marked_return":"I",
            "radial_return":"exp(-t)",
            "primitive_family_class_mod4":1,
            "positive_root_norm":str(norm),
            "positive_norm_parts":{"carrier":str(carrier_weight.dot(carrier_weight)),
                                   "family":str(family_weight.dot(family_weight))},
            "neutral_word_on_three_original_marks":original["vacuum_word"],
            "r113_values":None,
            "first_transport_test":{
                "from":"original charged current and its adjoint in the E8 affine vacuum",
                "into":"whole original weighted R113 word quotient before outer selection",
                "constraints":["same full gauge and family charge", "adjoint and orientation",
                               "same physical product domain", "the original time map",
                               "same bound absolute empty reference"],
                "native_middle_charge":original["r113_binding"]["charge_of_original_middle"],
                "rejected_weak_charge":original["r113_binding"]["charge_of_r032_weak_comparison"],
                "required_operations":["charged insertion", "charge projector", "one positive time step",
                                       "adjoint return", "fourth-root action", "absolute empty"],
                "transport_proved":False}}


def affine_transport_constraint(root,original_readout):
    anchor = affine_return(root)
    supplied = original_readout.get("affine_return_anchor")
    if supplied is not None:
        return {"status":"AFFINE_ANCHOR_CHECKED_ON_SUPPLIED_FORMS","source_refs":anchor["source_refs"],
                "original_ordered_value":None,"supplied_anchor":supplied}
    # Rank-one compressed data from the affine vacuum. A positive
    # completion fixes the vacuum but leaves the charged step free.
    p = sp.diag(1,0)
    charged = sp.Matrix([[0,0],[1,0]])
    clock = sp.diag(1,sp.I)
    trials = []
    for s in (sp.Rational(1,2),sp.Rational(3,4)):
        time_step = sp.diag(1,s)
        whole,direct,outside,_ = return_forms(sp.eye(2),p,charged,charged.H,time_step)
        require(clock*time_step==time_step*clock and time_step*p==p and
                whole==outside==s*p and direct==sp.zeros(2),"positive completion comparison failed")
        trials.append({"charged_time_eigenvalue":str(s),"normalized_ordered_return":str(whole[0,0])})
    return {"status":"TRANSPORT_IDENTITY_UNCOMPUTED",
            "source_refs":anchor["source_refs"],
            "native_source_status":original_readout.get("status"),
            "affine_input":{"class":1,"radial_grade":1,"norm":"1","phase":"I",
                            "radial_time_known":True},
            "weaker_data_test":{"status":"UNDERDETERMINED_ON_DECLARED_INTERFACE",
                "fixed":["vacuum step", "affine instantaneous norm", "order-four action", "adjoint", "positivity"],
                "missing":["original physical time intertwining on the charged line"],
                "finite_completions":trials,
                "scope":"algebraic affine return data with the R113 time map omitted"},
            "exact_observation":"P A* Q T(t) Q A P with the original E8 charged insertion A",
            "compare":"V* T_R113(t) V on the affine vacuum and charged current line, in the same time and absolute normalization",
            "original_ordered_value":None}


def checked_affine_anchor(data,source_id,time_id,gram,pivots,selector,read):
    """Synthesize the marked vacuum/current isometry in a supplied quotient."""
    require(data.get("source_id")==source_id and data.get("time_id")==time_id and
            all(isinstance(data.get(key),str) and data[key].strip() for key in (
                "native_charge_dictionary_binding","affine_adjoint_binding","reference_binding","time_calibration_binding")),
            "affine anchor needs the original charge, adjoint, time calibration and reference")
    metric = gram.extract(pivots,pivots)
    from reconstruction import scalar
    omega = metric.inv()*gram.extract(pivots,[data["identity_index"]])
    require(omega.H*metric*omega==sp.Matrix([[1]]),"the original empty word must represent the normalized vacuum")
    a,adjoint,clock,step = (read(data[k]) for k in ("first","last","quarter_action","step"))
    theta = scalar(data["expected_charged_step"])
    require(theta.is_Rational and 0<theta<1,"a calibrated charged step in (0,1) is required")
    metric_adjoint = metric.inv()*a.H*metric
    require(adjoint==metric_adjoint and positive(metric*step) and positive(metric*(sp.eye(metric.rows)-step)),
            "anchor requires a metric adjoint and one positive contractive step")
    charged = a*omega
    v = omega.row_join(charged)
    residuals = {
        "metric":v.H*metric*v-sp.eye(2),
        "outer_selection":selector*v-v*sp.diag(1,0),
        "primitive_quarter":clock*v-v*sp.diag(1,sp.I),
        "original_time":step*v-v*sp.diag(1,theta),
        "adjoint_return":adjoint*charged-omega,
        "adjoint_annihilates_vacuum":adjoint*omega,
    }
    defects = {k:matrix.applyfunc(sp.expand).rank() for k,matrix in residuals.items()}
    norm = (omega.H*metric*omega)[0]
    value = None
    if norm:
        value = sp.simplify((omega.H*metric*adjoint*step*a*omega)[0]/norm)
    return {"status":"EXACT_MARKED_AFFINE_ANCHOR_ON_SUPPLIED_FORMS" if not any(defects.values())
                     else "REFUTED_SUPPLIED_AFFINE_ANCHOR",
            "defect_ranks":defects,"candidate_column_rank":v.rank(),
            "ordered_return":str(value) if value is not None else None,
            "physical_source_identity_proved":False,"tested_only_on_vacuum_and_first_current":True}


def checked_return(data,source_id,time_id,metric,selector,read):
    require(data.get("source_id")==source_id and data.get("time_id")==time_id and
            all(isinstance(data.get(key),str) and data[key].strip()
                for key in ("field_binding","charge_binding","time_step_binding","reference_binding")),
            "ordered return requires the same original source, time, reference and charge bindings")
    first,last,step = (read(data[k]) for k in ("first","last","step"))
    require(first.H*metric==metric*last,"return operator is not the metric adjoint of the first insertion")
    whole,direct,outside,cross = return_forms(metric,selector,first,last,step)
    values = {"direct_insertion_zero":not any(selector*first*selector),
              "outside_return_rank":outside.rank(),"full_return_rank":whole.rank(),
              "projected_at_every_step_rank":direct.rank(),"mixed_rank":cross.rank(),
              "outside_return_positive":positive(metric*outside),
              "whole_return_positive":positive(metric*whole),
              "physical_state_or_continuum_derived":False}
    if "ordered_return_form" in data:
        observed = read(data["ordered_return_form"])
        values["full_return_matches_independent_form"] = whole==selector*observed*selector
    if "quarter_action" in data:
        action = read(data["quarter_action"])
        eye = sp.eye(metric.rows)
        phase_power = data["character_power"]
        require(type(phase_power) is int and 0<=phase_power<4,"declare a fourth-root character index")
        preserves_selector = selector*action==action*selector
        unitary = action.H*metric*action==metric and action**4==eye
        charged = action*first==sp.I**phase_power*first*action
        values["quarter"] = {"unitary_order_divides_four":unitary,
                             "preserves_selector":preserves_selector,
                             "preserves_time_step":action*step==step*action,
                             "declared_phase_matches":charged,
                             "primitive_character":charged and gcd(phase_power,4)==1 and bool(first.rank())}
    return values,whole,direct,outside


def checked_temporal_modes(data,read,metric,selector,step):
    # A scalar rate ratio is only inferred when both accessible excitation
    # columns are nonzero and each lies in a single time eigenspace.
    state = read(data["state_projector"])
    n = metric.rows
    eye = sp.eye(n)
    require(positive(metric*step) and positive(metric*(eye-step)),"time comparison needs a positive contraction")
    require(state*state==state and state.H*metric==metric*state and state.rank()==1 and
            selector*state==state and step*state==state,
            "need one normalized invariant accessible state ray")
    rates = {}
    for label in ("reference","relational"):
        probe = read(data[label])
        vector = selector*probe*state
        if not any(vector):
            rates[label] = {"status":"INVISIBLE_ON_SELECTED_ALGEBRA"}
            continue
        adjoint = metric.inv()*vector.H*metric
        t = sp.simplify((adjoint*step*vector).trace()/(adjoint*vector).trace())
        single = not any((step*vector-t*vector).applyfunc(sp.expand))
        rates[label] = {"status":"ISOLATED_EIGENMODE" if single else "MULTIPLE_TIME_COMPONENTS",
                        "step_eigenvalue":str(t) if single else None}
    observed = [v.get("step_eigenvalue") for v in rates.values()]
    ratio = None
    if all(x is not None and 0<sp.Rational(x)<1 for x in observed):
        ratio = f"log({observed[0]})/log({observed[1]})"
    return {"modes":rates,"dimensionless_rate_ratio":ratio}
