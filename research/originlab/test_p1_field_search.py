import json
from pathlib import Path
import subprocess
import sys

import pytest
import sympy as sp

from p1_field_search import controls
from source import COMPONENTS, MARKS, NATIVE, evaluate
from verify import digest


ROOT = Path(__file__).resolve().parents[2]


def test_original_folds_replay_all_charged_jets():
    report = controls(ROOT)
    assert report["charged_components_checked"] == 12
    assert [case["energies"] for case in report["cases"]] == [
        ["7/8","3/8","5/16","17/16","25/16"],["7/8","3/8","3/16","15/16","23/16"]]
    for case in report["cases"]:
        assert {row["weak_weight"] for row in case["channels"]} == {-1,1}
        for row in case["channels"]:
            assert row["same_covariance_and_heat_under_within_face_motion"]
            assert not row["pure_motion_preserving_fixed_hypercharge_possible"]
            assert row["independent_heat_mutant"] == "EXCLUDED_IN_DECLARED_CHARGED_JET_CLASS"
    assert report["original_r113_values"] is None


def test_supplied_joint_source_recovers_complex_intrinsic_field(tmp_path):
    def form(matrix):
        return [[[str(sp.re(value)),str(sp.im(value))] for value in row] for row in matrix.tolist()]

    # Two anchored charges and rational heat factors.
    transfer=sp.diag(sp.S.Half,sp.Rational(1,4))
    even=sp.diag(sp.Rational(4,5),sp.Rational(16,17))
    odd=sp.eye(2)-even
    off=even*transfer
    p=even.row_join(off).col_join(off.row_join(odd))
    h=sp.diag(transfer,transfer)
    y=sp.diag(sp.Rational(1,3),sp.S.Half,sp.Rational(1,3),sp.S.Half)
    rotation=sp.Matrix([[0,-sp.I],[-sp.I,0]])
    body=sp.diag(rotation,rotation)
    anchor=sp.zeros(4)
    anchor[0,3]=2
    anchor[3,0]=-2
    derivative=body+anchor
    field={"energies":[1,2],"hypercharges":["1/3","1/2"],"radial_base":"1/2","exponent_scale":1,
           "two_face_covariance_tangent":form(p*derivative-derivative*p),
           "two_face_heat_tangent":form(h*derivative-derivative*h),
           "two_face_hypercharge_tangent":form(y*anchor-anchor*y)}
    evidence=tmp_path/"evidence.txt"
    evidence.write_text("Exact joint collar arithmetic control")
    ref={"path":evidence.name,"sha256":digest(evidence)}
    native=tmp_path/"native.json"
    domain={"domain_id":"charged-collar-control","source_originals":[ref]}
    native.write_text(json.dumps(domain))
    binding={"path":native.name,"sha256":digest(native)}
    field.update(source_id=domain["domain_id"]+":"+binding["sha256"],time_id="T",
                 **dict.fromkeys(("two_face_frame_binding","independently_anchored_hypercharge_binding",
                                  "covariance_and_heat_tangents_binding","absolute_reference_binding"),"control"))
    unit=form(sp.eye(1))
    bundle={"schema":"tfpt.original-blocked-source.v1","domain_id":domain["domain_id"],
            "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
            "matrix_convention":"sesquilinear_forms_in_original_word_basis",
            "ensemble_scope":"whole_original_ensemble","representation":"original_exact_blocked_endpoint_source",
            "retained_words":[[]],"source_refs":[ref],"component_bindings":dict.fromkeys(COMPONENTS,"control"),
            "time":{"id":"T","kind":"original_euclidean_ordered","T":1,"t_U":0,"t_Y":1},
            "gram_form":unit,"physical_selector_form":unit,"boundary_pairing_form":unit,
            "scalar_factor":[1,0],"reference_scalar":[2,0],"ordered_source_forms":dict.fromkeys(MARKS,unit),
            "p1_charged_jet":field}
    attached=tmp_path/"attached.json"

    def read():
        attached.write_text(json.dumps({**domain,"original_domain_binding":binding,"whole_source_bundle":bundle}))
        return evaluate(tmp_path,attached)["p1_charged_jet"]

    answer=read()
    assert answer["status"] == "INTRINSIC_CHARGED_VARIATION_REQUIRED_IN_DECLARED_CLASS"
    assert answer["intrinsic_field"] == [["0","I"],["-I","0"]]
    assert not answer["physical_source_identity_proved"]
    field["two_face_heat_tangent"]=form(sp.zeros(4))
    assert read()["status"] == "EXCLUDED_IN_DECLARED_CHARGED_JET_CLASS"
    field["time_id"]="other time"
    with pytest.raises(ValueError,match="same source"):
        read()


def test_public_runner_continues_without_local_originals(tmp_path):
    subprocess.run(["git","init","-q",str(tmp_path)],check=True)
    subprocess.run(["git","-C",str(tmp_path),"-c","user.name=Test","-c","user.email=test@example.invalid",
                    "commit","-q","--allow-empty","-m","test"],check=True)
    out=tmp_path/"runs"
    output=subprocess.run([sys.executable,"-B",str(ROOT/"research/originlab/main.py"),"--repo",str(tmp_path),
                           "--out",str(out),"--campaign","r113","--campaign","p0","--campaign","f1",
                           "--jobs","1","--size","2","--mutations","2","--word-depth","2"],
                          capture_output=True,text=True)
    assert output.returncode == 0, output.stderr
    first=out/"round-0001"
    report=json.loads((first/"summary.json").read_text())
    assert report["context"]["atlas_provenance"] == "local_atlas_unavailable"
    assert report["outcomes"]["p0"] == "finite_search"
    assert report["outcomes"]["f1"] == "OPEN"
    r113=json.loads((first/"r113.json").read_text())
    assert r113["missing"] == [str(NATIVE)] and r113["original_values"] == dict.fromkeys(MARKS)


def test_public_r006_replays_the_byte_identical_archived_evidence(tmp_path):
    from verify import CAR_ORIGINAL_HASHES, car_replay

    report=car_replay(tmp_path,5)
    assert report["evidence_inputs"] == CAR_ORIGINAL_HASHES
    assert set(report["archived_inputs"]) == set(CAR_ORIGINAL_HASHES)
    assert report["shared_moments"] == 2245
    assert report["discriminator"]["shortest_length"] == 5
