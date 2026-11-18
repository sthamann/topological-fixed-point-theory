"""Regression checks for affine search and the bound endpoint comparison."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from p1_origin_hypotheses import search as hypotheses
from q0_rule_search import run, search
from source import COMPONENTS, MARKS
from verify import digest


class SourceRuleTests(unittest.TestCase):
    def test_affine_rule_ablation_across_grades(self):
        results = [search(depth) for depth in (1, 2, 3)]
        self.assertEqual([r["stages"][-1]["linear_dimension"] for r in results], [0, 2, 6])
        self.assertEqual(results[-1]["grade_dimensions"], [1, 3, 4, 7])
        self.assertEqual(len(results[-1]["independent_one_window_probes"]), 6)
        separator = results[-1]["separating_word"]
        self.assertEqual(separator["left"], [["E", -2]])
        self.assertEqual(separator["unit_word_candidates"], ["1/16", "9/16"])
        self.assertTrue(all(r["original_values"] == dict.fromkeys(MARKS) for r in results))

    def test_probe_uses_original_selector_time_and_reference(self):
        def diagonal(first, second):
            return [[[first, 0], [0, 0]], [[0, 0], [second, 0]]]

        with TemporaryDirectory() as folder:
            root = Path(folder)
            native = root/"native.json"
            reference_file = root/"declared-control.txt"
            reference_file.write_text("Synthetic arithmetic control\n")
            reference = {"path": reference_file.name, "sha256": digest(reference_file)}
            base = {"domain_id": "q0-control", "source_originals": [reference]}
            native.write_text(json.dumps(base))
            anchor = {"path": native.name, "sha256": digest(native)}
            probe = {
                "source_id": f"q0-control:{anchor['sha256']}", "time_id": "two-windows",
                "operator_dictionary_binding": "fixture", "two_window_time_binding": "fixture",
                "absolute_reference_binding": "fixture", "step": "1",
                "left_index": 1, "right_index": 1,
                "left_current_word": [["E", -2]], "right_current_word": [["E", -2]],
            }
            bundle = {
                "schema": "tfpt.original-blocked-source.v1", "domain_id": "q0-control",
                "input_role": "arithmetic_control_only", "source_refs": [reference],
                "native_inputs_sha256": anchor["sha256"],
                "matrix_convention": "sesquilinear_forms_in_original_word_basis",
                "ensemble_scope": "whole_original_ensemble",
                "representation": "original_exact_blocked_endpoint_source",
                "component_bindings": dict.fromkeys(COMPONENTS, "fixture"),
                "retained_words": [[], ["native-fixture-letter"]],
                "time": {"id": "two-windows", "kind": "original_euclidean_ordered",
                         "T": "2", "t_U": "1/2", "t_Y": "3/2"},
                "gram_form": diagonal(1, 2), "physical_selector_form": diagonal(1, 2),
                "boundary_pairing_form": diagonal(1, 0), "scalar_factor": [1, 0],
                "reference_scalar": [2, 0],
                "ordered_source_forms": dict.fromkeys(MARKS, diagonal("1/4", "1/4")),
                "q0_rule_probe": probe,
            }
            path = root/"source.json"
            data = dict(base, original_domain_binding=anchor, whole_source_bundle=bundle)
            path.write_text(json.dumps(data))
            result = run(root, path)["native_comparison"]
            self.assertEqual(result["surviving_pair_indices"], [0])
            self.assertEqual(result["probe"]["input_role"], "arithmetic_control_only")
            self.assertEqual(result["probe"]["absolute_numerator"], "1/4")
            self.assertEqual(result["probe"]["absolute_reference"], "2")
            self.assertFalse(result["physical_source_identity_proved"])
            bundle["physical_selector_form"] = diagonal(1, 0)
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "outer selector"):
                run(root, path)
            bundle["physical_selector_form"] = diagonal(1, 2)
            probe["step"] = "1/2"
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "two windows"):
                run(root, path)

    def test_degree_four_retains_nonunique_affine_forms(self):
        result = search(4)
        self.assertEqual(result["grade_dimensions"], [1, 3, 4, 7, 13])
        self.assertEqual(result["word_count"], 192)
        self.assertEqual(result["stages"][-1]["linear_dimension"], 13)
        self.assertEqual(result["separating_word"]["unit_word_candidates"], ["1/16", "9/16"])

    def test_report_hypotheses_and_negative_controls_run_by_default(self):
        result = hypotheses(max_cone_order=4)
        self.assertTrue({"minimal_self_dual_relations", "directed_anchor_groupoid",
                         "fourier_reciprocal_sources", "independent_dyadic_composition"}
                        <= result["hypotheses"].keys())
        evidence = result["evidence"]
        self.assertEqual(evidence["self_dual_relations"]["smith"], [1, 4, 4])
        self.assertFalse(evidence["self_dual_relations"]["scalar_reflection_positive"])
        self.assertEqual(evidence["orientation_flags"]["flags"], 30)
        self.assertEqual(evidence["reciprocal_sources"]["discrete_cap_probabilities"], ["7/10", "13/14"])
        self.assertFalse(result["physical_solution"])

    def test_singularity_and_modular_cases_are_selectable(self):
        result = hypotheses(["symmetric_singularity", "modular_index_kms"], 4)
        self.assertEqual(set(result["evidence"]), {"cone", "modular"})
        self.assertEqual(result["evidence"]["cone"]["order_four"]["fixed_point_coefficients"],
                         ["5/16", "-1/16", "-3/16", "-1/16"])
        self.assertIsNone(result["original_source_binding"])
        with self.assertRaises(ValueError):
            hypotheses(["unregistered-claim"])
