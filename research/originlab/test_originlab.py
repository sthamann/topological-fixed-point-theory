from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest

import sympy as sp

from language import Evaluator, Term, TypeMismatch, positive
from source import COMPONENTS, MARKS, compression_defect, evaluate, source_input_key
from verify import digest, r102, separator


class OriginLabTests(unittest.TestCase):
    def test_psd_with_null_pivot(self):
        self.assertFalse(positive(sp.Matrix([[0,1],[1,0]])))
        self.assertTrue(positive(sp.Matrix([[1,1],[1,1]])))

    def test_typed_process(self):
        evaluator = Evaluator({"x":sp.ImmutableMatrix([[1,2]])})
        with self.assertRaises(TypeMismatch):
            evaluator.evaluate(Term.parse(["compose",["atom","x"],["atom","x"]]))
        operator = evaluator.evaluate(Term.parse(["glue",["atom","x"],["mirror",["atom","x"]]]))
        self.assertEqual(operator,sp.Matrix([[5]]))

    def test_exact_separator_and_hidden_step(self):
        state = sp.Matrix([1,0])
        first = sp.Matrix([[sp.Rational(1,2),0],[0,sp.Rational(1,2)]])
        second = first+sp.Matrix([[0,sp.Rational(1,4)],[sp.Rational(1,4),0]])
        result = separator((state,{0:first}),(state,{0:second}),3)
        self.assertEqual(result["shortest_length"],2)
        self.assertEqual(compression_defect(sp.eye(1),state.H*second*state,state.H*second**2*state),
                         sp.Matrix([[sp.Rational(1,16)]]))

    def test_rank_search_parameter_is_variable(self):
        self.assertFalse(r102(5,16)["solutions"])
        self.assertEqual({row["n"] for row in r102(3,24)["solutions"]},{4})

    def test_external_packet_change_rekeys_domain(self):
        with TemporaryDirectory() as folder:
            root = Path(folder)
            domain = root/"domain.json"
            packet = root/"packet.json"
            domain.write_text(json.dumps({"whole_source_bundle":packet.name}))
            missing = source_input_key(domain)
            packet.write_text("{}")
            first = source_input_key(domain)
            packet.write_text('{"updated":true}')
            self.assertEqual(len({missing,first,source_input_key(domain)}),3)

    def test_packet_evidence_change_rekeys_input(self):
        with TemporaryDirectory() as folder:
            root = Path(folder)
            domain, packet, evidence = (root/name for name in ("domain.json","packet.json","evidence.txt"))
            domain.write_text(json.dumps({"whole_source_bundle":packet.name}))
            packet.write_text(json.dumps({"source_refs":[{"path":evidence.name}]}))
            missing = source_input_key(domain,root)
            evidence.write_text("A")
            first = source_input_key(domain,root)
            evidence.write_text("B")
            self.assertEqual(len({missing,first,source_input_key(domain,root)}),3)

    def test_original_contraction_uses_same_packet_and_reference(self):
        with TemporaryDirectory() as folder:
            root = Path(folder)
            evidence = root/"test-input.json"
            evidence.write_text("{}")
            ref = {"path":evidence.name,"sha256":digest(evidence)}
            native = root/"native.json"
            domain = {"domain_id":"arithmetic-test","source_originals":[ref]}
            native.write_text(json.dumps(domain))
            result = evaluate(root,native)
            self.assertEqual(result["status"],"BLOCKED")
            self.assertTrue(all(v is None for v in result["original_values"].values()))

            def form(entry):
                return [[[entry,0],[entry,0]],[[entry,0],[entry,0]]]

            forms = {k:form(n) for k,n in zip(MARKS,(4,8,12,28))}
            bundle = {
                "schema":"tfpt.original-blocked-source.v1","domain_id":domain["domain_id"],
                "native_inputs_sha256":digest(native),"input_role":"arithmetic_control_only",
                "matrix_convention":"sesquilinear_forms_in_original_word_basis",
                "ensemble_scope":"whole_original_ensemble",
                "representation":"original_exact_blocked_endpoint_source",
                "source_refs":[ref],"component_bindings":dict.fromkeys(COMPONENTS,"test"),
                "retained_words":[[],["test"]],"time":{"id":"test","kind":"original_euclidean_ordered",
                                                       "T":"1","t_U":"1/4","t_Y":"3/4"},
                "gram_form":form(1),"physical_selector_form":form(1),"boundary_pairing_form":form(1),
                "scalar_factor":[1,0],"reference_scalar":[4,0],"ordered_source_forms":forms,
            }
            path = root/"packet-domain.json"
            domain["original_domain_binding"] = {"path":native.name,"sha256":digest(native)}
            domain["whole_source_bundle"] = bundle
            path.write_text(json.dumps(domain))
            result = evaluate(root,path)
            self.assertEqual(result["input_role"],"arithmetic_control_only")
            self.assertEqual(result["quotient_rank"],1)
            self.assertEqual(result["values"],dict(zip(MARKS,("1","2","3","7"))))
            self.assertEqual(result["connected_UY"],"1")
            self.assertFalse(result["physical_source_identity_proved"])

            bundle.update(representation="original_exact_finite_fiber_sum",
                          configuration_unit={"kind":"whole_original_history_with_both_bosonic_endpoints",
                                              "native_map_binding":"test"},
                          measure_rule={"kind":"original_finite_atomic_measure","source_binding":"test"})
            bundle["fiber_records"] = [dict(original_configuration_id=f"atom-{i}",
                                            original_measure_binding="test",whole_fiber_binding="test",
                                            original_measure_mass=["1/2",0],unnormalized_scalar_weight=[1,0],
                                            ordered_source_forms=forms) for i in range(2)]
            path.write_text(json.dumps(domain))
            self.assertEqual(evaluate(root,path)["values"],result["values"])
            bundle["fiber_records"][1]["original_configuration_id"] = "atom-0"
            path.write_text(json.dumps(domain))
            with self.assertRaises(ValueError):
                evaluate(root,path)


if __name__ == "__main__":
    unittest.main()
