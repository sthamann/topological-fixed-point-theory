import json
from pathlib import Path
from tempfile import TemporaryDirectory

import sympy as sp

from root_doublet import checked_doublet, root_doublet
from source import COMPONENTS, MARKS, evaluate
from verify import digest


def operators():
    vacuum, lower, upper = [sp.eye(3)[:,i] for i in range(3)]
    create_lower = lower*vacuum.T
    create_upper = upper*vacuum.T
    raising = upper*lower.T
    return (vacuum,create_lower,create_upper,create_lower.T,create_upper.T,
            raising,raising.T,sp.diag(0,-1,1),sp.diag(1,(-1-sp.I*sp.sqrt(3))/2,1),
            sp.diag(0,sp.Rational(1,3),0))


def test_root_string_and_temporal_split():
    root = root_doublet()
    assert root["normalized_middle_abs_squared"]=="1"
    assert root["cyclic_horizontal_test"]["same_color"]==2
    steps = [sp.diag(1,sp.Rational(1,2),sp.Rational(1,2)),
             sp.diag(1,sp.Rational(1,4),sp.Rational(1,2))]
    result = checked_doublet(sp.eye(3),*operators(),steps)
    assert result["status"]=="EXACT_SUPPLIED_ROOT_DOUBLET"
    assert [p["middle_time_commutator_rank"] for p in result["time_steps"]]==[0,1]
    assert [p["relative_ordered_return"] for p in result["time_steps"]]==["1/4","1/8"]
    bad = list(operators())
    bad[6] = -bad[6]
    assert checked_doublet(sp.eye(3),*bad,steps)["kinematic_defect_ranks"]["middle_adjoint"]==1


def test_endpoint_reads_exact_colour_action():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        original = root/"native.json"
        original.write_text(json.dumps({"domain_id":"control","source_originals":[
            {"path":"evidence.txt","sha256":digest(Path(__file__))}]}))
        (root/"evidence.txt").write_bytes(Path(__file__).read_bytes())
        binding = {"path":"native.json","sha256":digest(original)}
        domain = json.loads(original.read_text())
        _,a,b,ad,bd,e,f,h,_,y = operators()

        def form(operator):
            return [[[str(sp.re(v)),str(sp.im(v))] for v in row] for row in operator.tolist()]

        unit = form(sp.eye(3))
        data = {"source_id":"control:"+binding["sha256"],"time_id":"T","identity_index":0,
                **dict.fromkeys(("field_and_charge_binding","horizontal_mode_binding","time_steps_binding"),"control"),
                **dict(zip(("create_lower","create_upper","lower_end","upper_end","raising","lowering","root_weight"),
                           map(form,(a,b,ad,bd,e,f,h)))),
                "color":{"rational":form(sp.diag(1,-sp.S.Half,1)),"sqrt3":form(sp.diag(0,-sp.I/2,0))},
                "hypercharge":form(y),"steps":[unit,form(sp.diag(1,sp.S.Half,sp.S.Half))]}
        bundle = {"schema":"tfpt.original-blocked-source.v1","domain_id":"control",
                  "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
                  "matrix_convention":"sesquilinear_forms_in_original_word_basis",
                  "ensemble_scope":"whole_original_ensemble","source_refs":domain["source_originals"],
                  "representation":"original_exact_blocked_endpoint_source",
                  "component_bindings":dict.fromkeys(COMPONENTS,"control"),
                  "retained_words":[[],["lower"],["upper"]],
                  "time":{"id":"T","kind":"original_euclidean_ordered","T":1,"t_U":0,"t_Y":1},
                  "gram_form":unit,"physical_selector_form":form(sp.diag(1,0,0)),
                  "boundary_pairing_form":unit,"scalar_factor":[1,0],"reference_scalar":[7,0],
                  "ordered_source_forms":dict.fromkeys(MARKS,unit),"marked_root_doublet":data}
        attached = root/"attached.json"
        attached.write_text(json.dumps({**domain,"original_domain_binding":binding,"whole_source_bundle":bundle}))
        result = evaluate(root,attached)
        assert result["marked_root_doublet"]["status"]=="EXACT_SUPPLIED_ROOT_DOUBLET"
        assert result["marked_root_doublet"]["time_steps"][1]["relative_ordered_return"]=="1/4"
        assert result["values"]["Z0"]=="1/7"
        assert result["physical_source_identity_proved"] is False
