import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
import json
import sympy as sp

from reconstruction import compare_documents, controls, energy_moment_check, matrix, reconstruct, regional_check
from space import finite_space, finite_process_space
from source import COMPONENTS, evaluate
from verify import digest


class ReconstructionTests(unittest.TestCase):
    def test_global_phase_and_rank(self):
        result = controls()
        self.assertEqual(result["ghz"]["probe"],"XXX")
        self.assertEqual(result["ghz"]["values"],["1","-1"])
        self.assertEqual([x["gram_rank"] for x in result["rank_controls"]],[6,12,18,24,30,36])
        self.assertEqual(result["all_word_comparison"]["status"],"ALL_WORDS_EQUAL")
        self.assertEqual(result["distinguishing_comparison"]["word"],["X","X"])

    def test_null_quotient_rejects_inconsistent_operator(self):
        g = sp.ones(2)
        with self.assertRaisesRegex(ValueError,"null space"):
            reconstruct(g,{"X":sp.eye(2)},0)
        report,_ = reconstruct(g,{"X":g},0,step=g/2,twice=g/3)
        self.assertEqual(report["gram_rank"],1)
        self.assertEqual(report["two_step_compression"]["defect_rank"],1)

    def test_energy_moments_separate_equal_first_rates(self):
        gram = sp.ones(2)
        self.assertEqual(energy_moment_check(gram,gram,gram)["status"],"GENERATOR_INVARIANT")
        spread = energy_moment_check(gram,gram,5*gram/4)
        self.assertEqual(spread["minimum_additional_directions"],1)
        self.assertEqual(spread["defect_form"],[["1/4"]])
        witness = spread["finite_extension_witness"]
        metric = matrix(witness["gram"])
        for form in witness["energy_forms"]:
            energy = matrix(form)
            self.assertTrue(energy.is_positive_semidefinite)
            self.assertEqual(energy[0,0],1)
            self.assertEqual((energy*metric.inv()*energy)[0,0],sp.Rational(5,4))
        self.assertEqual(witness["third_moment_difference"],[["1/4"]])
        with self.assertRaisesRegex(ValueError,"negative energy moment"):
            energy_moment_check(gram,gram,gram/2)
        with self.assertRaisesRegex(ValueError,"zero-energy direction"):
            energy_moment_check(gram,sp.zeros(2),gram)

    def test_exact_operational_equivalence_requires_matching_interface(self):
        base = dict(interface={"region":"seam","reference":"A"},gram=[[1,0],[0,1]],
                    forms={"X":[[0,1],[1,0]]},identity_index=0,adjoints={"X":"X"})
        changed = {**base,"gram":[[1,0],[0,4]],"forms":{"X":[[0,2],[2,0]]}}
        self.assertTrue(compare_documents(base,changed)["cyclic_star_equivalence"])
        with self.assertRaisesRegex(ValueError,"marked interface"):
            compare_documents(base,{**changed,"interface":{"region":"different"}})

    def test_full_regional_algebra_moves(self):
        metric = sp.eye(2)
        p = sp.diag(1,0)
        x = sp.Matrix([[0,1],[1,0]])
        report = regional_check(metric,[metric,sp.zeros(2),sp.zeros(2),sp.zeros(2)],
            {"A":[metric,p],"B":[metric,metric-p]},[{"from":"A","to":"B","operator":x}])
        self.assertTrue(report["positive_joint_spinor_form"])
        self.assertTrue(report["transports"][0]["transports_full_algebra"])
        self.assertFalse(report["transports"][0]["source_and_target_algebras_distinct"])
        self.assertEqual(report["spatial_rank_modulo_identity_and_time"],0)

    def test_regional_motion_distinguishes_internal_actions_and_time(self):
        identity = sp.eye(3)
        spatial = [sp.diag(1,0,0),sp.diag(0,1,0),sp.zeros(3)]
        internal = {"A":[identity,*spatial[:2]]}
        report = regional_check(identity,[identity,*spatial],internal,[])
        self.assertEqual(report["spatial_rank_modulo_identity_and_time"],2)
        self.assertEqual(report["regional_motion"]["spatial_rank_modulo_time_and_internal_actions"],0)
        uniform = sp.ones(3)/3
        report = regional_check(identity,[identity,*spatial],{"A":[identity,uniform]},[])
        self.assertEqual(report["regional_motion"]["spatial_rank_modulo_time_and_internal_actions"],2)
        report = regional_check(identity,[spatial[0],*spatial],{"A":[identity,uniform]},[])
        self.assertEqual(report["regional_motion"]["spatial_rank_modulo_time_and_internal_actions"],1)

    def test_reconstruction_through_original_packet_entry(self):
        with TemporaryDirectory() as name:
            root = Path(name)
            domain = {"domain_id":"control","source_originals":[{"path":"proof","sha256":""}]}
            (root/"proof").write_text("arithmetic control")
            domain["source_originals"][0]["sha256"] = digest(root/"proof")
            (root/"native.json").write_text(json.dumps(domain))
            binding = {"path":"native.json","sha256":digest(root/"native.json")}
            unit = [[[1,0]]]
            packet = dict(schema="tfpt.original-blocked-source.v1",domain_id="control",
                          matrix_convention="sesquilinear_forms_in_original_word_basis",
                          ensemble_scope="whole_original_ensemble",
                          input_role="arithmetic_control_only",source_refs=domain["source_originals"],
                          native_inputs_sha256=binding["sha256"],component_bindings=dict.fromkeys(COMPONENTS,"control"),
                          retained_words=[[]],time={"T":1,"t_U":0,"t_Y":1,"id":"control","kind":"original_euclidean_ordered"},
                          gram_form=unit,physical_selector_form=unit,boundary_pairing_form=unit,
                          scalar_factor=[1,0],reference_scalar=[1,0],
                          representation="original_exact_blocked_endpoint_source",
                          ordered_source_forms=dict.fromkeys(("Z0","ZU","ZY","ZUY"),unit),
                          energy_moments={"source_id":f"control:{binding['sha256']}","time_id":"control",
                                          "generator_binding":"control","ground_energy_binding":"control",
                                          "word_domain_binding":"control","first":unit,"second":[[["5/4",0]]]},
                          physical_reconstruction={"identity_index":0,"forms":{"X":unit},
                                                   "adjoints":{"X":"X"},
                                                   "density_model":{"density":[[1]],"word_matrices":[[[1]]]}})
            path = root/"input.json"
            path.write_text(json.dumps({**domain,"original_domain_binding":binding,"whole_source_bundle":packet}))
            result = evaluate(root,path)
            self.assertTrue(result["physical_reconstruction"]["density_model"]["rank_formula_applies"])
            self.assertEqual(result["energy_moments"]["minimum_additional_directions"],1)
            self.assertFalse(result["physical_source_identity_proved"])
            packet["energy_moments"]["time_id"] = "different"
            path.write_text(json.dumps({**domain,"original_domain_binding":binding,"whole_source_bundle":packet}))
            with self.assertRaisesRegex(ValueError,"same original source, time"):
                evaluate(root,path)

    def test_exhaustive_bounded_spaces(self):
        self.assertEqual(finite_space(4,3)["count"],114)
        r = finite_process_space(2,1)
        self.assertEqual(r["tested"],81)
        self.assertEqual(r["positive_contractions"],4)


if __name__ == "__main__":
    unittest.main()
