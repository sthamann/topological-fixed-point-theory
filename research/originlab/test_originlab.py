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

            bundle["normalization_scope"] = "relative_mark_contraction"
            del bundle["reference_scalar"]
            path.write_text(json.dumps(domain))
            relative = evaluate(root,path)
            self.assertEqual(relative["status"],"EXACT_RELATIVE_CONTRACTION")
            self.assertEqual(relative["connected_UY"],"1")
            self.assertEqual(relative["normalized_by_common_empty"],result["values"])
            self.assertTrue(all(v is None for v in relative["original_values"].values()))

            from r113_reflected_radial import RADIAL_DATA, RADIAL_PROFILE_SHA256

            profile = root/RADIAL_DATA
            profile.parent.mkdir(parents=True)
            profile.write_bytes((Path(__file__).resolve().parents[2]/RADIAL_DATA).read_bytes())
            bundle["ordered_source_forms"] = {**forms,"ZU":form("73/100"),"ZUY":form("211/100")}
            bundle["radial_moment_probe"] = {
                "source_id":f"{domain['domain_id']}:{digest(native)}","time_id":"test",
                "radial_profile_sha256":RADIAL_PROFILE_SHA256,
                "positive_radial_law_binding":"control","normalized_late_profile_binding":"control",
                "paired_radial_operator_binding":"control","late_radial_form":form("73/100"),
                "early_radial_square_form":form("7633/40000"),"late_radial_square_form":form("7633/40000"),
                "cross_radial_form":form("5329/40000"),
            }
            path.write_text(json.dumps(domain))
            moment = evaluate(root,path)["radial_moment_probe"]
            self.assertEqual(moment["conditional_connected_UY_upper_bound"],"-63/6250")
            self.assertTrue(moment["strict_negative_bound"])
            self.assertFalse(moment["supplied_response_exceeds_bound"])
            bundle["ordered_source_forms"]["ZUY"] = form("219/100")
            path.write_text(json.dumps(domain))
            self.assertTrue(evaluate(root,path)["radial_moment_probe"]["supplied_response_exceeds_bound"])
            bundle["ordered_source_forms"] = {k:bundle["ordered_source_forms"][k] for k in ("Z0","ZU")}
            bundle["normalization_scope"] = "relative_radial_moment_bound"
            path.write_text(json.dumps(domain))
            bound_only = evaluate(root,path)
            self.assertEqual(bound_only["status"],"EXACT_RADIAL_MOMENT_BOUND")
            self.assertIsNone(bound_only["connected_UY"])
            self.assertTrue(bound_only["radial_moment_probe"]["strict_negative_bound"])
            del bundle["radial_moment_probe"]
            bundle["ordered_source_forms"] = forms

            del bundle["normalization_scope"]
            bundle["reference_scalar"] = [4,0]

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

    def test_autonomy_discovers_then_replays_a_domain_bound_packet(self):
        from autonomy import discovered_source_evidence, execute, research_followups, specifications
        from evidence_gates import assess
        from original_discovery import inventory

        with TemporaryDirectory() as folder:
            root = Path(folder)
            sources = root/"research/originlab"
            sources.mkdir(parents=True)
            evidence = sources/"evidence.txt"
            evidence.write_text("Synthetic source discovery test")
            ref = {"path":str(evidence.relative_to(root)),"sha256":digest(evidence)}
            native = root/"domain.json"
            native.write_text(json.dumps({"domain_id":"discovery-fixture","source_originals":[ref]}))
            def form(n):
                return [[[n,0]]]
            packet = {
                "schema":"tfpt.original-blocked-source.v1","domain_id":"discovery-fixture",
                "native_inputs_sha256":digest(native),"input_role":"source_reviewed_original",
                "matrix_convention":"sesquilinear_forms_in_original_word_basis",
                "ensemble_scope":"whole_original_ensemble","representation":"original_exact_blocked_endpoint_source",
                "component_bindings":dict.fromkeys(COMPONENTS,"synthetic fixture"),
                "source_refs":[ref],"retained_words":[[]],
                "time":{"id":"test","kind":"original_euclidean_ordered","T":"1","t_U":"1/4","t_Y":"3/4"},
                "gram_form":form(1),"physical_selector_form":form(1),"boundary_pairing_form":form(1),
                "scalar_factor":[1,0],"reference_scalar":[1,0],
                "ordered_source_forms":{k:form(v) for k,v in zip(MARKS,(1,2,3,7))},
            }
            (sources/"reviewed.json").write_text(json.dumps({"whole_source_bundle":packet}))
            (sources/"control.json").write_text(json.dumps({**packet,"input_role":"arithmetic_control_only"}))
            detected = inventory(root,native)
            self.assertEqual(len(detected["matching"]),1)
            self.assertEqual(detected["rejected"][0]["reason"],"CONTROL_OR_UNREVIEWED_ROLE")
            specs = [job for _,job in specifications(4,native,None,[],False,detected["matching"])]
            full_search = next(job for job in specs if job["lane"]=="f1_solution")
            self.assertEqual(full_search["stage"],3)
            self.assertEqual(full_search["domain"],str(native))
            spec = next(job for job in specs if job["lane"]=="discovered_original_packet")
            result = execute(spec,root)
            self.assertEqual(result["connected_UY"],"1")
            self.assertFalse(result["physical_solution"])
            completed = {"spec":spec,"result":result,"queue_status":"done"}
            jobs = research_followups([completed],native)
            self.assertEqual(jobs[0][1]["lane"],"discovered_q0_comparison")
            original = {"status":"BLOCKED"}
            received,selection = discovered_source_evidence(original,[completed])
            self.assertEqual(selection["status"],"REPLAYED_SOURCE_REQUIRES_REVIEW")
            physical = assess(received)
            self.assertEqual(physical["gates"]["G3"]["original_values"],result["values"])
            self.assertFalse(physical["physical_solution"])
            packet["ordered_source_forms"]["ZUY"] = form(8)
            (sources/"reviewed-second.json").write_text(json.dumps(packet))
            other_entry = next(entry for entry in inventory(root,native)["matching"]
                               if entry["path"].endswith("reviewed-second.json"))
            other_spec = {**spec,"entry":other_entry}
            other = {"spec":other_spec,"result":execute(other_spec,root),"queue_status":"done"}
            received,selection = discovered_source_evidence(original,[completed,other])
            self.assertEqual(received,original)
            self.assertEqual(selection["status"],"MULTIPLE_SOURCE_PACKETS_REQUIRE_REVIEW")


if __name__ == "__main__":
    unittest.main()
