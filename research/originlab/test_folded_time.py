from pathlib import Path
from tempfile import TemporaryDirectory
import json
import os
import unittest

import sympy as sp

from affine_time import GLUE, HYPERCHARGES, LEGS, TRIALITY
from autonomy import execute
from folded_time import SOURCES
from source import COMPONENTS, MARKS, evaluate
from verify import digest


ROOT = Path(os.environ.get("TFPT_REPO_ROOT",Path(__file__).resolve().parents[2]))


class FoldedTimeTests(unittest.TestCase):
    def test_both_original_quarters_and_gauge_projection(self):
        data = execute({"lane":"original_folded_time"},ROOT)
        first, second = data["cases"]
        self.assertEqual(first["marked_energies_in_H8_time"],["3/2","7"])
        self.assertEqual(second["marked_energies_in_H8_time"],["5/2","7"])
        self.assertEqual([row["positive_gauge_average"]["coefficient_of_Y"] for row in data["cases"]],
                         ["-1/8","-7/40"])
        self.assertEqual(first["twisted_mirror_gaps_relative_to_minus_e8_reference_in_H8_time"],
                         {"e1_control":"-6","reflected_p1":"-5/2","reflected_p2":"3"})
        self.assertEqual(second["fixed_background_translated_lattice_bottom_by_coset"]["vector"][
                            "ground_momenta_doubled"],[[0,0,0,2,0,0,0,0],[2,0,0,0,0,0,0,0]])
        self.assertEqual([q["norm_squared"] for q in second["clock_preserving_mirror_backgrounds"]],
                         ["149/128","277/128"])
        full_mirror = second["whole_lattice_mirror"]
        self.assertEqual(full_mirror["target_empty_impulse"],"0")
        self.assertEqual(full_mirror["fixed_background_marked_energies_in_H8_time"],["9/2","9"])
        self.assertEqual(full_mirror["transported_background_marked_energies_in_H8_time"],["5/2","7"])
        self.assertEqual(full_mirror["relative_character_of_pulled_back_fixed_generator_at_K_angle_2pi"],"I")
        self.assertEqual(data["clock_choice_on_fixed_word"]["symmetric_exact_reduced_radial_probe_ratio"],"81/85")
        self.assertIsNone(data["physical_r113_endpoint_values"])

    def test_weighted_grid_accepts_fold_time_and_tests_held_out_entry(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            refs = []
            for bound in SOURCES:
                (root/bound).parent.mkdir(parents=True,exist_ok=True)
                (root/bound).write_bytes((ROOT/bound).read_bytes())
                refs.append({"path":str(bound),"sha256":digest(root/bound)})
            domain = {"domain_id":"weighted-clock-control","source_originals":refs}
            original = root/"base.json"
            original.write_text(json.dumps(domain))
            binding = {"path":original.name,"sha256":digest(original)}
            source_id = domain["domain_id"]+":"+binding["sha256"]

            def form(value):
                real,imag = sp.expand_complex(value).as_real_imag()
                return [[[str(real),str(imag)]]]

            # r=2^-16. Directly multiply the marked point separations and
            # the original quarter-one leg factors. Use a nonunit reference.
            width = sp.Integer(2)**-16
            axis = [{"physical_time":j,"radial_factor":str(width**j)} for j in (1,2)]
            grid = []
            for m in (1,2):
                row = []
                for n in (1,2):
                    a,b,c = 1,sp.I*width**m,-width**(m+n)
                    value = sp.Integer(2)**(13*m+2*n) * a*b*c/((a-b)*(b-c)*(a-c))
                    row.append(form(5*value))
                grid.append(row)
            time_grid = {"source_id":source_id,"time_id":"T","native_circle_marks":["1","I","-1"],
                         **dict.fromkeys(("native_field_and_charge_binding","physical_to_radial_time_binding",
                                          "absolute_reference_binding"),"arithmetic control"),
                         "ordered_fields":list(LEGS),"internal_glue_classes":list(GLUE),
                         "color_trialities":list(TRIALITY),"hypercharges":list(HYPERCHARGES),
                         "left_axis":axis,"right_axis":[dict(r) for r in axis],"ordered_word_forms":grid}
            bundle = {"schema":"tfpt.original-blocked-source.v1","domain_id":domain["domain_id"],
                      "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
                      "matrix_convention":"sesquilinear_forms_in_original_word_basis",
                      "ensemble_scope":"whole_original_ensemble","source_refs":refs,
                      "representation":"original_exact_blocked_endpoint_source",
                      "component_bindings":dict.fromkeys(COMPONENTS,"arithmetic control"),
                      "retained_words":[[]],"time":{"id":"T","kind":"original_euclidean_ordered",
                                                   "T":4,"t_U":1,"t_Y":3},
                      "gram_form":form(1),"physical_selector_form":form(1),"boundary_pairing_form":form(1),
                      "scalar_factor":[1,0],"reference_scalar":[5,0],
                      "ordered_source_forms":dict.fromkeys(MARKS,form(1)),"marked_affine_time_grid":time_grid}
            attached = root/"attached.json"

            def run():
                attached.write_text(json.dumps({**domain,"original_domain_binding":binding,
                                                "whole_source_bundle":bundle}))
                return evaluate(root,attached)["marked_affine_time_grid"]

            result = run()
            quarters = result["folded_clock_comparison"]["quarters"]
            self.assertEqual(result["status"],"REFUTED_SUPPLIED_RADIAL_TIME_GRID")
            self.assertEqual(result["mixed_time_ratios"]["measured"],result["mixed_time_ratios"]["radial"])
            self.assertEqual(result["separable_deviation_mismatch_indices"],[])
            self.assertTrue(quarters[0]["absolute_matches_supplied_grid"])
            self.assertFalse(quarters[1]["matches_up_to_one_global_normalization"])
            self.assertTrue(result["folded_clock_comparison"]["stationary_mixture_of_original_quarters"][
                            "positive_normalized_mixture"])

            time_grid["ordered_word_forms"][1][1] = form(1)
            changed = run()
            self.assertEqual(changed["separable_deviation_mismatch_indices"],[(1,1)])
            self.assertEqual(changed["folded_clock_comparison"]["stationary_mixture_of_original_quarters"][
                             "held_out_mismatch_indices"],[(1,1)])
            time_grid["right_axis"][1]["physical_time"] = 3
            self.assertEqual(run()["folded_clock_comparison"]["status"],
                             "SUPPLIED_RADIAL_TIME_MAP_NOT_A_COMMON_EXPONENTIAL")


if __name__=="__main__":
    unittest.main()
