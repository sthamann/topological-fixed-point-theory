from pathlib import Path
from tempfile import TemporaryDirectory
import json
import os
import unittest

import sympy as sp

from autonomy import execute
from empirical import audit
from selection import clifford_matrices, matrix_candidate, operator_checks
from source import COMPONENTS, MARKS, evaluate
from verify import digest

ROOT = Path(os.environ.get("TFPT_REPO_ROOT",Path(__file__).resolve().parents[2]))


class SelectionTests(unittest.TestCase):
    def test_full_and_split_invariance_and_original_carrier(self):
        result = execute({"lane":"four_mark_selection"},ROOT)
        self.assertEqual(result["geometry"]["cycle_lattice_determinant"],16)
        self.assertTrue(result["geometry"]["s4_mod2_actions_intertwine"])
        self.assertEqual(result["geometry"]["primitive_classes"],[1]*4)
        self.assertEqual(result["split"]["invariant_forms"]["S4"]["nondegenerate"],0)
        self.assertEqual(result["split"]["invariant_forms"]["C3"]["nondegenerate"],1)
        self.assertEqual(result["split"]["clock_and_anchor_preserving_intertwiners"],6)
        self.assertEqual(result["unsplit"]["single_face_values"],[0,1])
        self.assertEqual(result["unsplit"]["readout_branches"],
                         [{"value":0,"remaining_candidates":1},{"value":1,"remaining_candidates":1}])
        charged = result["clifford"]["same_four_state_carrier_for_s4_primitive_loops"]
        self.assertEqual(charged["minimum_charged_carrier_dimension"],15)
        self.assertEqual(charged["minimum_joint_carrier_with_clifford"],16)
        self.assertEqual(charged["conditional_joint_car_control"]["coherent_return"],"1/2")
        oriented = charged["conditional_oriented_s4_car_control"]
        self.assertEqual(oriented["verified_permutations"],24)
        self.assertEqual(oriented["invariant_hermitian_return_expectations"],["1","-1"])
        self.assertEqual(result["native_fixed_quark_block"]["status"],"REFUTED_SCOPED")

    def test_commuting_standins_fail_the_pair_relation(self):
        ops = [sp.Matrix([[0,-1,0],[1,0,0],[0,0,0]])]*4
        raw = [a.tolist() for a in ops]
        raw = [[[int(x) for x in row] for row in a] for a in raw]
        r = matrix_candidate({"operators":raw})
        self.assertEqual(r["status"],"REFUTED_SCOPED")
        self.assertEqual(r["defects"]["pair_1_2"]["state_norm_squared"],"4")

    def test_charged_operations_can_expose_a_hidden_clifford_failure(self):
        ops = [sp.diag(a,sp.zeros(3)) for a in clifford_matrices()]
        eye = sp.eye(7)
        loops = [eye[:,target]*eye[:,start].T for start,target in ((6,0),(5,6),(4,5),(0,4))]
        cartan = [-eye[:,i]*eye[:,i].T for i in (6,5,4)]
        isolated = operator_checks(ops,eye,eye[:,0],[-1]*4)
        joint = operator_checks(ops,eye,eye[:,0],[-1]*4,(loops,cartan))
        self.assertEqual(isolated["status"],"EXACT_CLIFFORD_ON_CYCLIC_SPACE")
        self.assertEqual(joint["primitive_loops"]["status"],"EXACT_SUPPLIED_WEIGHT_ACTION")
        self.assertEqual(joint["cyclic_rank"],7)
        self.assertEqual(joint["status"],"REFUTED_SCOPED")
        self.assertEqual(joint["defects"]["square_1"]["state_norm_squared"],"0")
        self.assertEqual(joint["defects"]["square_1"]["counterexample_word"],["loop1*"])

    def test_exact_packet_exercises_common_original_consumer(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root/"input.txt"
            evidence.write_text("arithmetic control")
            ref = {"path":evidence.name,"sha256":digest(evidence)}
            original = root/"original.json"
            native = {"domain_id":"selection-test","source_originals":[ref]}
            original.write_text(json.dumps(native))
            binding = {"path":original.name,"sha256":digest(original)}
            source_id = native["domain_id"]+":"+binding["sha256"]

            def encode(m):
                return [[[int(sp.re(x)),int(sp.im(x))] for x in row] for row in m.tolist()]

            eye = sp.eye(4)
            identity = encode(eye)
            operations = [encode(sp.I*op) for op in clifford_matrices()]
            loops = [eye[:,target]*eye[:,start].T for start,target in ((3,0),(2,3),(1,2),(0,1))]
            cartan = [-eye[:,i]*eye[:,i].T for i in (3,2,1)]
            bundle = {
                "schema":"tfpt.original-blocked-source.v1","domain_id":native["domain_id"],
                "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
                "matrix_convention":"sesquilinear_forms_in_original_word_basis",
                "ensemble_scope":"whole_original_ensemble","source_refs":[ref],
                "representation":"original_exact_blocked_endpoint_source",
                "component_bindings":dict.fromkeys(COMPONENTS,"declared arithmetic test"),
                "retained_words":[[],["x"],["y"],["z"]],
                "time":{"id":"T","kind":"original_euclidean_ordered","T":1,"t_U":"1/4","t_Y":"3/4"},
                "gram_form":identity,"physical_selector_form":identity,"boundary_pairing_form":identity,
                "scalar_factor":[1,0],"reference_scalar":[4,0],"ordered_source_forms":dict.fromkeys(MARKS,identity),
                "four_mark_selection":{"source_id":source_id,"time_id":"T","identity_index":0,
                    **dict.fromkeys(("orientation_binding","adjoint_binding","state_binding",
                                     "charge_dictionary_binding","face_code_binding"),"declared"),
                    "operations":{str(i):op for i,op in enumerate(operations)},"square_phases":[1]*4,
                    "primitive_loop_forms":[encode(op) for op in loops],
                    "cartan_forms":[encode(op) for op in cartan],
                    "transport":{"binding":"declared","step":identity,"transported_operations":operations},
                    "neutral_ordered_return_form":encode(eye[:,0]*eye[:,0].T)},
            }
            path = root/"packet.json"
            path.write_text(json.dumps({**native,"original_domain_binding":binding,"whole_source_bundle":bundle}))
            result = evaluate(root,path)
            selection = result["four_mark_original"]
            self.assertEqual(selection["status"],"EXACT_CLIFFORD_ON_CYCLIC_SPACE")
            self.assertEqual(selection["observed_cyclic_square_phases"],[1]*4)
            self.assertEqual(selection["face_refinement"]["compatible_refinements"][0]["squares"],[0]*4)
            self.assertEqual(selection["primitive_loops"]["status"],"EXACT_SUPPLIED_WEIGHT_ACTION")
            self.assertEqual(selection["declared_transport_intertwines"],[True]*4)
            self.assertEqual(selection["absolute_neutral_return_of_supplied_operations"],"1/4")
            self.assertTrue(selection["matches_ordered_return_operator_on_quotient"])
            self.assertFalse(selection["physical_source_identity_proved"])
            self.assertEqual(set(result["values"]),set(MARKS))
            bundle["four_mark_selection"]["neutral_ordered_return_form"] = encode(eye[:,1]*eye[:,1].T)
            path.write_text(json.dumps({**native,"original_domain_binding":binding,"whole_source_bundle":bundle}))
            selection = evaluate(root,path)["four_mark_original"]
            self.assertTrue(selection["matches_original_ordered_return"])
            self.assertFalse(selection["matches_ordered_return_operator_on_quotient"])

    def test_spectral_source_rejects_reflection_counterexample(self):
        r = execute({"lane":"spectral_source_controls"},ROOT)
        self.assertTrue(r["small"]["connected_excludes_zero"])
        self.assertEqual(r["small"]["reflection"]["combined_kernel"]["status"],"POSITIVE_DEFINITE")
        for name in ("unweighted","stronger_profile"):
            self.assertEqual(r[name]["reflection"]["combined_kernel"]["status"],"POSITIVE_DEFINITE")
        failed = r["reflection_counterexample"]["reflection"]["combined_kernel"]
        self.assertEqual(failed["status"],"REFUTED_SCOPED")
        self.assertEqual(failed["witness_indices"],list(range(5)))
        self.assertEqual(len(failed["witness_coefficients"]),5)
        self.assertEqual(r["small"]["unbounded_lebesgue_normalization"]["status"],"REFUTED_SCOPED")
        self.assertEqual(r["small"]["physical_original_values"],dict.fromkeys(MARKS))

    def test_frozen_comparison_preserves_scheme(self):
        r = audit(ROOT)
        self.assertEqual(len(r["alpha"]["stages"]),5)
        self.assertEqual(r["alpha"]["open_target"],"ALPHA.QUILLEN.EXACT.01")
        self.assertEqual(r["physical_gate_G7"],"OPEN")
        self.assertTrue(all(row["scheme_matched_prediction"]=="BLOCKED_MISSING_SOURCE_TO_POLE_MAP"
                            for row in r["leptons"]["rows"]))


if __name__ == "__main__":
    unittest.main()
