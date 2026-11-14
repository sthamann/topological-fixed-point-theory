import json
from pathlib import Path

import sympy as sp

from affine_polarization import connected_response, controls, root_matrix
from affine_time import radial_word
from source import COMPONENTS, MARKS, evaluate
from verify import digest


def test_original_polarization_receipt():
    result = controls()
    assert result["marked_charge_plane_rank"] == 2
    assert result["silent_cartan_dimension"] == 6
    assert result["sampled_charge_plane_matrix_span_rank"] == 4
    assert result["weighted_r113_values"] is None
    assert result["matched_charge_words"]["matching_words"] == 2
    assert result["matched_charge_words"]["joint_projection_rank"] == 2
    assert result["matched_charge_words"]["connected_offdiagonal_at_half"] == ["-25/16","25/16"]


def test_original_consumer_compares_supplied_forms(tmp_path):
    def form(value):
        value = sp.expand(value)
        return [[[str(sp.re(value)),str(sp.im(value))]]]

    original = tmp_path/"native.json"
    evidence = {"path":str(Path(__file__).resolve()),"sha256":digest(Path(__file__))}
    domain = {"domain_id":"polarization-control","source_originals":[evidence]}
    original.write_text(json.dumps(domain))
    binding = {"path":"native.json","sha256":digest(original)}
    roots = root_matrix()
    basis = roots[:, :2]
    basis_input = [[str(value) for value in row] for row in basis.tolist()]
    samples = []
    for x,y in [(sp.S.Half,sp.S.Half),(sp.Rational(1,4),sp.Rational(1,4))]:
        # The artificial common scalar must not be mistaken for absolute agreement.
        bare = 2*sp.cancel(radial_word(x,y),extension=sp.I)
        expected = bare*(basis.T*basis+connected_response(roots,basis,basis,x,y))
        samples.append({"left_radial_factor":str(x),"right_radial_factor":str(y),
                        "bare_word_form":form(bare),
                        "prepared_word_forms":[[form(v) for v in row] for row in expected.tolist()]})
    data = {"source_id":domain["domain_id"]+":"+binding["sha256"],"time_id":"T",
            "native_preparation_and_adjoint_binding":"control",
            "physical_to_radial_time_binding":"control","absolute_reference_binding":"control",
            "left_polarizations":basis_input,"right_polarizations":basis_input,"samples":samples}
    packet = {"schema":"tfpt.original-blocked-source.v1","domain_id":domain["domain_id"],
              "input_role":"arithmetic_control_only","native_inputs_sha256":binding["sha256"],
              "matrix_convention":"sesquilinear_forms_in_original_word_basis",
              "ensemble_scope":"whole_original_ensemble",
              "representation":"original_exact_blocked_endpoint_source",
              "retained_words":[[]],"source_refs":[evidence],
              "component_bindings":dict.fromkeys(COMPONENTS,"control"),
              "time":{"id":"T","kind":"original_euclidean_ordered","T":1,"t_U":0,"t_Y":1},
              "gram_form":form(1),"physical_selector_form":form(1),"boundary_pairing_form":form(1),
              "scalar_factor":[1,0],"reference_scalar":[3,0],
              "ordered_source_forms":dict.fromkeys(MARKS,form(1)),"affine_polarization":data}
    attached = tmp_path/"attached.json"

    def read():
        attached.write_text(json.dumps({**domain,"original_domain_binding":binding,
                                        "whole_source_bundle":packet}))
        return evaluate(tmp_path,attached)["affine_polarization"]

    result = read()
    assert result["status"] == "EXACT_SUPPLIED_WARD_SHAPE"
    assert all(not row["absolute_affine_bare_match"] for row in result["samples"])
    assert not result["physical_source_identity_proved"]
    alternative_root = controls()["matched_charge_words"]["right_roots"][1]
    alternative = sp.Matrix(alternative_root).applyfunc(sp.Rational)
    alt_roots = sp.Matrix.hstack(-alternative-roots[:,1],roots[:,1],alternative)
    alt_samples = []
    for sample in samples:
        x,y = (sp.Rational(sample[k]) for k in ("left_radial_factor","right_radial_factor"))
        bare = 2*sp.cancel(radial_word(x,y),extension=sp.I)
        expected = bare*(basis.T*basis+connected_response(alt_roots,basis,basis,x,y))
        alt_samples.append({**sample,"prepared_word_forms":[[form(v) for v in row] for row in expected.tolist()]})
    data["alternative"] = {"right_root":alternative_root,"samples":alt_samples,"native_field_binding":"control"}
    result = read()["alternative"]
    assert result["status"] == "EXACT_SUPPLIED_WARD_SHAPE"
    assert all(row["affine_bare_word"] is None for row in result["samples"])
    del data["alternative"]
    saved = samples[0]["prepared_word_forms"][0][0]
    samples[0]["prepared_word_forms"][0][0] = form(0)
    result = read()
    assert result["status"] == "DISTINGUISHED_SUPPLIED_POLARIZATION"
    assert [row["prepared_residual_rank"] for row in result["samples"]] == [1,0]
    samples[0]["prepared_word_forms"][0][0] = saved
    samples[0]["bare_word_form"] = form(0)
    samples[0]["prepared_word_forms"] = [[form(0),form(0)],[form(0),form(0)]]
    assert read()["status"] == "INCONCLUSIVE_ZERO_BARE_WORD"
