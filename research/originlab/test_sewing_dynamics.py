from pathlib import Path
from tempfile import TemporaryDirectory
import json
import os
import unittest

import sympy as sp

from autonomy import execute, research_followups
from return_search import completion_search, HADAMARD
from return_moments import reconstruct_return
from sewing_dynamics import complete_graphs, finite_gauge_link, link_sewing, two_cells
from source import COMPONENTS, MARKS, NATIVE, evaluate
from verify import digest

ROOT = Path(os.environ.get("TFPT_REPO_ROOT",Path(__file__).resolve().parents[2]))


class SewingDynamicsTests(unittest.TestCase):
    def test_rank_selection_has_an_orientation_condition(self):
        report = complete_graphs()
        self.assertEqual(report["positive_rank_selection"],[4])
        self.assertEqual(report["k4_transposition_traces"],{"vertex_relations":1,"cycles":-1})
        self.assertEqual([r["determinant"] for r in report["rows"]],[3,16,125,1296,16807])

    def test_gauge_and_local_parity_change_the_observable_rates(self):
        zero = two_cells({"reference_weight":0})
        normal = two_cells({"reference_weight":2})
        opposite = two_cells({"parity_signs":[1,-1]})
        self.assertEqual(zero["ungauged_graph"]["ground_multiplicity"],2)
        self.assertEqual(normal["ungauged_equal"]["energy_multiplicities"],
                         [["0",1],["2",4],["4",6],["6",4],["8",1]])
        self.assertEqual(normal["unprojected_time_rates"]["ratio"],"1/2")
        self.assertEqual(normal["binary_link"]["projected_time_rates"],
                         {"single_reference":None,"reference_graph_pair":"12","graph_graph_pair":"16"})
        self.assertEqual(opposite["binary_link"]["equal"]["ground_multiplicity"],4)
        self.assertEqual(normal["local_projection"]["second_order_spectrum"],[["-2",1],["0",3]])

    def test_fourfold_link_can_be_trivial_after_closed_gauss_projection(self):
        report = finite_gauge_link(4)
        self.assertEqual(report["physical_dimension"],2)
        self.assertEqual(report["closed_primitive_hop_rank"],0)
        self.assertEqual(report["boundary_sectors"][1]["primitive_hop_rank"],1)
        self.assertEqual(report["boundary_sectors"][1]["right_boundary_character"],"I")
        sew = link_sewing(4)
        self.assertTrue(sew["associative_isometry"])
        self.assertEqual([r["neutral_flux_sectors"] for r in sew["topologies"]],[1,4,16])
        custom = link_sewing(4,[{"name":"reversed","vertices":3,"edges":[[1,0],[1,2],[2,0]]},
                                {"name":"point","vertices":1,"edges":[]}])
        self.assertEqual([r["neutral_flux_sectors"] for r in custom["topologies"]],[4,1])

    def test_synthesis_reacts_to_an_independently_chosen_answer(self):
        # The inverse algorithm sees only these oracle answers, not weights.
        hidden = sp.Matrix([2,3,5,8])
        values = [str(x*x/8) for x in HADAMARD*hidden]
        data = {"gap":8,"input_role":"declared_finite_oracle",
                "projected_return_observations":values,"require_positive_weights":True}
        spec = {"lane":"return_generator_search","input":data}
        report = execute(spec,ROOT)
        self.assertEqual(report["status"],"UNDERDETERMINED_ON_DECLARED_INTERFACE")
        self.assertEqual(report["same_interface"]["status"],"ALL_WORDS_EQUAL")
        self.assertEqual(len(report["counterexamples"]),2)
        self.assertEqual(report["number_of_candidates_up_to_overall_sign"],8)
        followup = research_followups([{"queue_status":"done","spec":spec,"result":report}],ROOT/NATIVE)
        _,job = followup[0]
        branches = execute(job,ROOT)["hypothetical_branches"]
        self.assertEqual(len(branches),8)
        signs = job["input"]["signs"]
        observed = sum(w*s for w,s in zip(hidden,signs))**2/8
        result = completion_search({**data,"opposite_parity_observations":[{"signs":signs,"value":str(observed)}]})
        self.assertEqual(result["candidates"],[["2","3","5","8"]])

    def test_positive_reconstruction_and_blind_moment(self):
        moments = [sp.Matrix([[(sp.Rational(1,4)**i+sp.Rational(3,4)**i)/2]]) for i in range(8)]
        found = reconstruct_return(moments[:6],moments[6:])
        self.assertEqual(found["rank"],2)
        self.assertEqual(found["matched_held_out_moments"],2)
        self.assertEqual(reconstruct_return(moments[:4])["status"],"WORD_SPACE_STILL_GROWING")
        self.assertEqual(reconstruct_return(moments[:6],[moments[6]+sp.eye(1)])["first_difference"],6)
        self.assertEqual(reconstruct_return([sp.Matrix([[x]]) for x in (1,2,4,8)])["status"],
                         "REFUTED_POSITIVE_CONTRACTION")

    def test_complement_return_and_metric_on_actual_source_reader(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root/"evidence.txt"
            evidence.write_text("arithmetic control")
            ref = {"path":evidence.name,"sha256":digest(evidence)}
            original = root/"native.json"
            domain = {"domain_id":"return-test","source_originals":[ref]}
            original.write_text(json.dumps(domain))
            binding = {"path":original.name,"sha256":digest(original)}
            source_id = "return-test:"+binding["sha256"]
            gram = sp.diag(2,3)
            selector = sp.diag(1,0)
            step = sp.diag(1,sp.Rational(1,2))
            first = sp.Matrix([[0,0],[1,0]])
            last = gram.inv()*first.H*gram
            def form(op):
                return [[[str(sp.re(v)),str(sp.im(v))] for v in row] for row in (gram*op).tolist()]
            bundle = {"schema":"tfpt.original-blocked-source.v1","domain_id":domain["domain_id"],
                      "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
                      "matrix_convention":"sesquilinear_forms_in_original_word_basis",
                      "ensemble_scope":"whole_original_ensemble","source_refs":[ref],
                      "representation":"original_exact_blocked_endpoint_source",
                      "component_bindings":dict.fromkeys(COMPONENTS,"arithmetic"),
                      "retained_words":[[],["a"]],"time":{
                          "id":"T","kind":"original_euclidean_ordered","T":1,"t_U":"1/4","t_Y":"3/4"},
                      "gram_form":form(sp.eye(2)),"physical_selector_form":form(selector),
                      "boundary_pairing_form":form(sp.eye(2)),"scalar_factor":[1,0],"reference_scalar":[2,0],
                      "ordered_source_forms":dict.fromkeys(MARKS,form(sp.eye(2))),
                      "ordered_return":{"source_id":source_id,"time_id":"T",
                          **dict.fromkeys(("field_binding","charge_binding","time_step_binding","reference_binding"),"arithmetic"),
                          "first":form(first),"last":form(last),"step":form(step),
                          "ordered_return_form":form(sp.Rational(3,4)*selector),
                          "quarter_action":form(sp.diag(1,sp.I)),"character_power":1}}
            path = root/"packet.json"
            def run():
                path.write_text(json.dumps({**domain,"original_domain_binding":binding,"whole_source_bundle":bundle}))
                return evaluate(root,path)
            report = run()["ordered_return"]
            self.assertEqual(report["absolute_readouts"],
                             {"full":"3/8","project_after_each_event":"0","via_complement":"3/8"})
            self.assertTrue(report["full_return_matches_independent_form"])
            self.assertTrue(report["quarter"]["primitive_character"])
            self.assertEqual(report["input_role"],"arithmetic_control_only")
            # Rebind the normalized vacuum and current in a second exact
            # test; the charged transfer must still match the same anchor.
            gram = sp.eye(2)
            last = first.H
            bundle["gram_form"] = form(sp.eye(2))
            bundle["physical_selector_form"] = form(selector)
            bundle["boundary_pairing_form"] = form(sp.eye(2))
            bundle["ordered_source_forms"] = dict.fromkeys(MARKS,form(sp.eye(2)))
            bundle["ordered_return"].update(first=form(first),last=form(last),step=form(step),
                quarter_action=form(sp.diag(1,sp.I)),ordered_return_form=form(selector/2))
            bundle["affine_return_anchor"] = {
                "source_id":source_id,"time_id":"T","identity_index":0,
                **dict.fromkeys(("native_charge_dictionary_binding","affine_adjoint_binding","reference_binding",
                                 "time_calibration_binding"),"arithmetic"),
                **{k:bundle["ordered_return"][k] for k in ("first","last","step","quarter_action")},
                "expected_charged_step":"1/2"}
            bundle["ordered_return_moments"] = {
                "source_id":source_id,"time_id":"T",
                **dict.fromkeys(("charged_alphabet_binding","physical_time_step_binding",
                                 "absolute_reference_binding"),"arithmetic"),
                "moments":[form(selector/2**i) for i in range(4)],
                "held_out":[form(selector/16)]}
            actual = run()
            self.assertEqual(actual["affine_return_anchor"]["ordered_return"],"1/2")
            self.assertEqual(actual["affine_return_anchor"]["status"],"EXACT_MARKED_AFFINE_ANCHOR_ON_SUPPLIED_FORMS")
            self.assertEqual(actual["ordered_return_moment_reconstruction"]["matched_held_out_moments"],1)
            # Exact arithmetic control with time unit log(4), using the
            # original signed fields and a held-out mixed response.
            axis = [{"physical_time":1,"radial_factor":"1/4"},
                    {"physical_time":2,"radial_factor":"1/16"}]
            grid = [[-sp.Rational(240,4913)-sp.I*sp.Rational(128,4913),
                     -sp.Rational(4032,283985)-sp.I*sp.Rational(256,56797)],
                    [-sp.Rational(4032,283985)-sp.I*sp.Rational(256,56797),
                     -sp.Rational(65280,16974593)-sp.I*sp.Rational(8192,16974593)]]
            bundle["marked_affine_time_grid"] = {
                "source_id":source_id,"time_id":"T","native_circle_marks":["1","I","-1"],
                **dict.fromkeys(("native_field_and_charge_binding","physical_to_radial_time_binding",
                                 "absolute_reference_binding"),"arithmetic control"),
                "ordered_fields":["J^(-p1)","J^gamma","J^p2"],"internal_glue_classes":[1,2,1],
                "color_trialities":[0,1,2],"hypercharges":["0","-1/3","1/3"],
                "left_axis":axis,"right_axis":axis,
                "ordered_word_forms":[[form(2*value*selector) for value in row] for row in grid]}
            self.assertEqual(run()["marked_affine_time_grid"]["status"],"EXACT_SUPPLIED_RADIAL_TIME_GRID")
            bundle["marked_affine_time_grid"]["ordered_word_forms"][1][1] = form(selector)
            self.assertEqual(run()["marked_affine_time_grid"]["mismatch_indices"],[[1,1]])
            bundle["marked_affine_time_grid"]["internal_glue_classes"] = [1,1,1]
            with self.assertRaisesRegex(ValueError,"signed original three legs"):
                run()
            del bundle["marked_affine_time_grid"]
            bundle["affine_return_anchor"]["expected_charged_step"] = "3/4"
            self.assertEqual(run()["affine_return_anchor"]["defect_ranks"]["original_time"],1)
            bundle["ordered_return"]["source_id"] = "different source"
            with self.assertRaisesRegex(ValueError,"same original source"):
                run()

    def test_original_anchor_is_affine_and_followup_is_uncomputed(self):
        original = execute({"lane":"original_affine_return"},ROOT)
        self.assertEqual(original["adjoint_ordered_return"],"1")
        self.assertEqual(original["quarter_marked_return"],"I")
        self.assertIsNone(original["r113_values"])
        rows = [{"queue_status":"done","spec":{"lane":"r113"},"result":{"status":"BLOCKED"}},
                {"queue_status":"done","spec":{"lane":"original_affine_return"},"result":original}]
        _,job = research_followups(rows,ROOT/NATIVE)[0]
        self.assertIsNone(execute(job,ROOT)["original_ordered_value"])

    def test_original_marked_time_and_single_fibre_obstruction(self):
        result = execute({"lane":"original_affine_time"},ROOT)
        self.assertEqual([leg["D_phase"] for leg in result["marked_charges"]],["I","-1","I"])
        self.assertEqual(result["radial_spectrum"]["analytic_separation_rank"],"infinite")
        self.assertEqual(result["radial_spectrum"]["distinct_positive_time_grid_rank"],4)
        self.assertEqual(result["clocks"]["parabolic_on_original_marks"]["-1"]["first_derivatives"],
                         ["-1/4","0"])
        self.assertEqual(result["clocks"]["parabolic_on_normalized_current"]["held_out_rank_one_defect"],
                         "17/324")
        self.assertEqual(result["clocks"]["balanced_positive_radial_composition"][
                         "radial_factor_1_over_4_steps"],["1/3","3/4","1/3"])
        self.assertEqual(result["finite_degree_two_restriction"][
                         "lower_bound_complex_car_modes_for_scalar_quadratic_embedding"],124)
        self.assertEqual(result["finite_degree_two_restriction"][
                         "minimum_level_for_linear_majorana_action_and_quadratic_free_currents"],30)
        self.assertIsNone(result["weighted_r113_values"])


if __name__=="__main__":
    unittest.main()
