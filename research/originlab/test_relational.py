import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import sympy as sp

from autonomy import execute, extra_specs
from coherence import compare_memory, det_word, r104
from observability import compiler_control
from original_paths import original_word
from relational import controls, flux_inference, network
from search import profile_job
from source import COMPONENTS, MARKS, NATIVE, evaluate, fibre_forms
from tomography import controls as tomography_controls, solve as tomography_solve
from verify import digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
# The tests also run against an isolated staging copy.
if not (ROOT/"semantic_explorer").exists():
    import os
    ROOT = Path(os.environ["TFPT_REPO_ROOT"])


class RelationalTests(unittest.TestCase):
    def test_absolute_timed_word_matches_dense_selected_trace(self):
        n = 3
        word = [[3,"annihilate"],[1,"create"],[1,"annihilate"],[2,"annihilate"],[1,"annihilate"]]
        times = ["1/7","1/7",0,"2/7","3/7","1/7"]
        q,z = sp.symbols("q z",nonzero=True)
        cap = sp.eye(8)[:,0]+z*sp.eye(8)[:,7]
        dual = sp.eye(8)[0,:]+sp.eye(8)[7,:]/z
        p = cap*dual/2

        def dense_letter(mode, op):
            x = sp.Matrix([[0,1],[0,0]])
            if op=="create":
                x = x.T
            return sp.kronecker_product(*(sp.diag(1,-1) if i<mode else x if i==mode else sp.eye(2)
                                          for i in range(n,0,-1)))

        for modulus in (None,2,3):
            def heat(t):
                t = q**(7*sp.Rational(t))
                return t*sp.eye(8)+(1-t)*p
            dense = sp.diag(*(int(modulus is None or i.bit_count()%modulus==0) for i in range(8)))
            for (i,op),t in zip(word,times):
                dense = (dense*heat(t)*dense_letter(i,op)).applyfunc(sp.expand)
            expected = sp.expand((dense*heat(times[-1])).trace())
            sparse = det_word(n,"1/7",word,times,modulus)
            self.assertEqual(sp.sympify(sparse["polynomial_in_q_z"],locals={"q":q,"z":z}),expected)

    def test_literal_r104_and_relational_blindness(self):
        result = r104(ROOT)
        self.assertEqual(result["coherent_gram"]["ranks"][-2:],[52666,65536])
        self.assertEqual(result["locality"]["four_mode_proper_support_pauli_checks"],175)
        self.assertEqual(result["markov"]["petz_trace_distance"],"1/2")
        self.assertEqual(result["number_sector"]["absolute_filtered_trace"],"1 + 16511*exp(-2*beta)")
        odd = r104(ROOT,5)
        self.assertEqual(odd["coherent_gram"]["first_full_depth"],2)
        self.assertEqual(odd["coherent_gram"]["first_phase_sensitive_gram_depth"],3)
        result = controls()
        self.assertEqual(result["curved"]["path_comparisons"][0]["trace"],"exp(I*pi/3)")
        self.assertFalse(result["curved"]["flat_on_supplied_faces"])
        self.assertEqual([x["regional"][0]["distinguished_by_this_charge"]
                          for x in result["charged_path_combinations"]],[False,True])
        self.assertEqual([x["indistinguishable_from_zero"] for x in result["charge_visibility"]],[True,False,False])
        self.assertEqual([x["covariant"] for x in result["stationary_symmetry"]],[True,False,True])
        self.assertEqual(result["mixed_gaussian"]["mixed_responses"][0]["connected"],"-3/8")
        self.assertEqual(len(result["inverse_flux"]["flux_residue_candidates"]),16)
        self.assertEqual(result["inverse_flux"]["next_probes"][0]["charge"],1)
        narrowed = flux_inference({"fundamental_loops":["square"],"phase_order":16,
            "observations":[{"charge":16,"winding":[1],"residue":0},
                            {"charge":1,"winding":[1],"residue":5}],"allowed_probes":[]})
        self.assertEqual(narrowed["flux_residue_candidates"],[[5]])

    def test_nonabelian_order_and_flat_face_scope(self):
        x = [[0,1],[1,0]]
        z = [[1,0],[0,-1]]
        i = [[1,0],[0,1]]
        data = {"vertices":["A","B","C"],"dimension":2,
                "edges":{"ab":{"from":"A","to":"B","matrix":x},
                         "bc":{"from":"B","to":"C","matrix":z},
                         "ac":{"from":"A","to":"C","matrix":i}},
                "path_pairs":[{"from":"A","to":"C","routes":[[["ab",1],["bc",1]],[["ac",1]]]}],
                "face_boundaries":[{"id":"contractible","base":"A","boundary":[["ab",1],["ab",-1]]}]}
        r = network(data)
        self.assertEqual(r["path_comparisons"][0]["holonomy"],[["0","1"],["-1","0"]])
        self.assertTrue(r["flat_on_supplied_faces"])
        self.assertFalse(r["all_graph_loops_trivial"])

    def test_reconstruction_and_exact_cap_readout(self):
        r = tomography_controls()
        self.assertTrue(r["restricted"]["full_reconstruction_on_accessible_operations"])
        self.assertEqual(r["restricted"]["positive_pair_extension_responses"],["3/8","1/8"])
        self.assertTrue(r["complete"]["full_linear_reconstruction"])
        channel = [[1,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,1]]
        effect = [[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        intervention = {"input_state":[[1,0],[0,0]],"effect":effect}
        data = {"dimensions":[2,2],"process":channel,"interventions":[{"m":intervention}],
                "observations":[{"word":["m"],"value":1}]}
        self.assertTrue(tomography_solve(data)["causal_process_exists"])
        intervention["effect"] = channel
        with self.assertRaisesRegex(ValueError,"normalization"):
            tomography_solve(data)
        left = {"interface":{"reference":"absolute"},"cap":[[1]],"operations":{"x":[[1]]}}
        right = {**left,"cap":[[0]]}
        self.assertEqual(compare_memory(left,right)["word"],[])
        compiler = compiler_control(ROOT)
        self.assertEqual(compiler["accessible_ranks_by_depth"],[1,11,26])
        self.assertEqual(compiler["extended_ranks_by_depth"],[1,13,34,36])

    def test_typed_search_reuses_counterexample(self):
        spec = next(case for case in extra_specs(HERE/"relational-cases.json") if case["lane"]=="origin_dsl")
        spec = {**spec,"observations":[{"dimension":2,"probe":"east-readout",
                                      "word":["candidate","candidate"],"value":"1/4","evidence_id":"control"}]}
        first = profile_job(spec)
        self.assertIn("combined",first["grammar_families"])
        self.assertTrue(any(row["term"][0]=="glue" for row in first["survivors"]))
        self.assertTrue(any(row["term"][0]=="compose" and row["status"]=="ill_typed" for row in first["excluded"]))
        self.assertTrue(any(row["status"]=="counterexample" for row in first["excluded"]))
        archived = {"excluded":first["excluded"],"questions":first["questions"]}
        second = profile_job(spec,archived)
        self.assertGreater(second["feedback"]["reused_negative_certificates"],0)
        self.assertEqual(first["survivors"],second["survivors"])

    def test_source_combines_kernels_and_checks_regional_holdout(self):
        def form(re,im=0):
            return [[[re,im]]]
        with TemporaryDirectory() as folder:
            root = Path(folder)
            original = root/"original.json"
            source = root/"source.txt"
            source.write_text("arithmetic control")
            ref = {"path":source.name,"sha256":digest(source)}
            domain = {"domain_id":"test","source_originals":[ref]}
            original.write_text(json.dumps(domain))
            binding = {"path":original.name,"sha256":digest(original)}
            source_id = "test:"+binding["sha256"]
            regional = {"source_id":source_id,"time_id":"clock","endpoint_charge":2,
                        "unit_charge_path_comparison_binding":"control unit charge",
                        **dict.fromkeys(("charge_spin_binding","endpoint_field_binding","transport_binding",
                                         "admissibility_binding"),"arithmetic control"),
                        "connection":{"vertices":["A","B"],"edges":{"a":{"from":"A","to":"B","phase_pi":0},
                                     "b":{"from":"A","to":"B","phase_pi":"1/2"}},
                                     "path_pairs":[{"from":"A","to":"B","routes":[[["a",1]],[["b",1]]]}]},
                        "paths":[{"id":"a","route":[["a",1]],"word":["a"],"form":form(2)},
                                 {"id":"b","route":[["b",1]],"word":["b"],"form":form(0,2)}]}
            memory = {"interface":{"source_id":source_id,"time_id":"clock",
                                   "reference":"same_bound_absolute_empty"},
                      "cap":[[1]],"operations":{"U":[[2]],"Y":[[3]],"A":[[1]],"B":[[1]]},
                      "mark_words":{"Z0":[],"ZU":["U"],"ZY":["Y"],"ZUY":["U","Y"]},
                      "regional_mark_words":{"a":["A"],"b":["B"]}}
            packet = {"schema":"tfpt.original-blocked-source.v1","domain_id":"test",
                      "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
                      "matrix_convention":"sesquilinear_forms_in_original_word_basis",
                      "ensemble_scope":"whole_original_ensemble","source_refs":[ref],
                      "representation":"original_exact_finite_fiber_sum",
                      "component_bindings":dict.fromkeys(COMPONENTS,"control"),
                      "configuration_unit":{"kind":"whole_original_history_with_both_bosonic_endpoints",
                                            "native_map_binding":"control"},
                      "measure_rule":{"kind":"original_finite_atomic_measure","source_binding":"control"},
                      "retained_words":[[],["input_basis_1"],["input_basis_2"]],
                      "time":{"id":"clock","kind":"original_euclidean_ordered","T":"1","t_U":"1/4","t_Y":"3/4"},
                      "scalar_factor":[1,0],"reference_scalar":[2,0],
                      "relational_path_forms":regional,"coherent_memory":memory}
            # A rank-one metric on three labelled words, including the
            # nontrivial null quotient shared by every endpoint form.
            def extend(q):
                return [[q[0][0]]*3 for _ in range(3)]
            for field in ("gram_form","physical_selector_form","boundary_pairing_form"):
                packet[field] = extend(form(1))
            for path in regional["paths"]:
                path["form"] = extend(path["form"])
            packet["joint_source_checks"] = {"source_id":source_id,"time_id":"clock","ward_binding":"control",
                "extra_forms":{"stress":extend(form(4)),"stressY":extend(form(20))},
                "mixed":[{"first":"stress","second":"ZY","joint":"stressY"}],
                "ward_identities":[{"id":"probe","coefficients":{"stress":1,"ZU":-1},"rhs":0}]}
            packet["fiber_records"] = [{"original_configuration_id":"path",
                 "whole_fiber_binding":"test","original_measure_binding":"test","original_measure_mass":[1,0],
                 "unnormalized_scalar_weight":[1,0],"ordered_window_kernels":{
                      "kind":"whole_history_fixed_bosonic_scalar_U","original_scalar_mark_binding":"test",
                      "U":2,"K0":extend(form(2)),"KY":extend(form(6))}}]
            attached = root/"attached.json"
            attached.write_text(json.dumps({**domain,"original_domain_binding":binding,"whole_source_bundle":packet}))
            result = evaluate(root,attached)
            self.assertEqual(result["values"],dict(zip(MARKS,["1","2","3","6"])))
            self.assertEqual(result["relational_original"]["absolute_responses"],["1","I"])
            self.assertEqual(result["relational_original"]["response_holonomy_check"]["status"],"DISTINGUISHED")
            self.assertEqual(result["coherent_memory"]["differences"],{"regional:b":"1 - I"})
            self.assertEqual(result["joint_source_checks"]["mixed"][0]["connected"],"4")
            self.assertTrue(result["joint_source_checks"]["linear_identities"]["probe"]["satisfied"])
            atom = packet["fiber_records"][0]
            atom["ordered_window_kernels"]["K0"] = extend(form(0))
            self.assertEqual(fibre_forms(atom,3,sp.S.Zero,sp.S.One)["ZUY"][0,0],12)

    def test_default_native_witness_has_no_promoted_source_values(self):
        r = execute({"lane":"r113_orientation"},ROOT)
        self.assertEqual(r["status"],"EXCLUDED_SCOPED")
        self.assertEqual(r["physical_original_values"],dict.fromkeys(MARKS))
        self.assertEqual(evaluate(ROOT)["relational_original"]["status"],"BLOCKED")
        native = execute({"lane":"r113","domain":str(ROOT/NATIVE)},ROOT)["native_endpoint_readout"]
        if native.get("status")=="UNAVAILABLE":
            self.skipTest("original runtime requires NumPy")
        self.assertEqual(native["values"],dict.fromkeys(MARKS))
        regional = original_word(ROOT)
        self.assertEqual(regional["vacuum_word"],"-I/4")
        self.assertEqual(regional["r113_binding"]["absolute_values"],[])
        self.assertFalse(regional["r113_binding"]["interchangeable_by_the_supplied_charge_dictionary"])


if __name__ == "__main__":
    unittest.main()
