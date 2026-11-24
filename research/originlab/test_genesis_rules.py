import unittest

import genesis_rules

from autonomy import execute, specifications
from genesis_rules import (clock_lattices, codes, cuts, double_cover, folds, norm_two_vectors, run, search, sign_kernel,
                           six, spinor_check, symmetric)


class GenesisRules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run()

    def selected(self, universe):
        return {row["selected"]: row for row in self.result["universes"][universe]["unique_selections"]}

    def test_rule_search_reports_ambiguity_and_empty_rules(self):
        objects = [{"n": 1, "clauses": {"a"}}, {"n": 1, "clauses": {"a", "b"}}, {"n": 2, "clauses": {"b"}}]
        found = search(objects, ["a", "b", "c"], lambda o: o["n"], lambda o: f"{o['n']}{sorted(o['clauses'])}")
        self.assertEqual(found["verdicts"], {"EMPTY": 4, "UNIQUE": 2, "AMBIGUOUS": 2})

    def test_minimum_is_not_mistaken_for_unique_satisfying_object(self):
        objects = [{"n": 1, "clauses": {"a"}}, {"n": 2, "clauses": {"a"}}]
        result = genesis_rules.search_certified(objects, ["a"], lambda o: o["n"], lambda o: str(o["n"]))
        row = result["unique_selections"][0]
        support = row["selection_support"][0]
        self.assertEqual(support["matching_labels"], 2)
        self.assertEqual(support["minimal_labels"], 1)
        self.assertTrue(support["requires_size_minimization"])


    def test_self_dual_doubly_even_selects_one_code_of_length_eight(self):
        row = self.selected("codes")["n=8 k=4 weights (1, 0, 0, 0, 14, 0, 0, 0, 1)"]
        self.assertIn(["doubly_even", "self_dual"], row["minimal_rules"])
        derived = self.result["universes"]["codes"]["derived"]
        self.assertEqual(derived["self_dual_doubly_even_codes_of_length_8"], 30)
        self.assertTrue(derived["one_orbit_under_relabelling"])
        self.assertEqual(derived["construction_a"]["norm_two_vectors"], 240)

    def test_dropping_a_clause_changes_the_code(self):
        row = self.selected("codes")["n=8 k=4 weights (1, 0, 0, 0, 14, 0, 0, 0, 1)"]
        self.assertEqual(row["ablation"]["doubly_even+self_dual without doubly_even"][1], ["n=2 k=1 weights (1, 0, 1)"])

    def test_cut_needs_both_glue_clauses_for_d5_a3(self):
        rows = self.selected("cuts")
        self.assertEqual(rows["A3 | D5"]["minimal_rules"], [["glue_has_conjugate_pair", "glue_has_neutral_class"]])
        self.assertEqual(rows["A2 | E6"]["minimal_rules"], [["glue_has_conjugate_pair"]])
        self.assertEqual(len(cuts(8)), 11)
        self.assertFalse(any(cuts(r) for r in range(2, 8)))

    def test_clocks_generate_an_even_self_dual_lattice_only_in_rank_eight(self):
        derived = self.result["universes"]["clocks"]["derived"]
        self.assertEqual(derived["even_self_dual_below_rank_8"], 0)
        self.assertEqual(derived["clocks_without_subclock"], [15, 20, 24, 30])
        self.assertEqual(derived["sample"]["norm_two_vectors"], 240)

    def test_clock_gram_of_a_non_unimodular_control_has_fewer_roots(self):
        control = next(c for c in clock_lattices(8, norms=(2,)) if c["det"] > 1)
        gram = [[control["seq"][abs(i-j)] for j in range(8)] for i in range(8)]
        self.assertLess(norm_two_vectors(gram), 240)

    def test_fold_with_a_shape_parameter_has_four_marks(self):
        rows = self.selected("folds")
        self.assertIn(["has_a_shape_parameter"], rows["(2, 2, 2, 2)"]["minimal_rules"])
        self.assertIn("(2, 3, 5)", rows)
        only_four = [d["marks"] for d in map(double_cover, range(2, 13, 2))
                     if d["modes_are_twice_marks"] and d["odd_modes_are_the_loops"] and d["cover_is_flat"]]
        self.assertEqual(only_four, [4])
        self.assertEqual(sum("flat" in f["clauses"] for f in folds()), 4)

    def test_five_signs_are_the_smallest_set_blind_to_third_order(self):
        rows = self.selected("signs")
        self.assertIn(["blind_to_third_order"], rows["n=5 degree 5"]["minimal_rules"])
        self.assertEqual([sign_kernel(n) for n in range(1, 6)], [0, 0, 0, 0, 1])
        five = spinor_check()
        self.assertEqual(five["states"], 16)
        self.assertTrue(five["third_moments_all_zero"])
        self.assertEqual(five["sum_of_cubed_charges"], "0")
        self.assertEqual(five["charge_3_plus_2_multiplicities"],
                         {"-1": 1, "-1/3": 3, "-1/6": 6, "0": 1, "1/2": 2, "2/3": 3})

    def test_blind_signs_completed_self_dually_need_a_partner_of_rank_three(self):
        found = self.result["composed_rule"]["result"]
        self.assertEqual(found["lattice_spanned_by_the_words"]["determinant"], 4)
        self.assertEqual(found["lattice_spanned_by_the_words"]["vectors_of_norm_two"], 40)
        self.assertEqual(found["partner_rank_mod_8"], 3)
        self.assertEqual([len(found["partners_by_rank"][r]) for r in (1, 2, 3)], [0, 0, 1])
        self.assertEqual(found["partners_by_rank"][3][0]["vectors_of_norm_two"], 12)

    def test_words_equal_flips_only_for_eight_signs(self):
        found = self.result["one_rule_variant"]["result"]
        self.assertEqual(found["n_with_words_of_flip_norm"], [8])
        self.assertEqual(found["n_integral_and_even"], [8, 16, 24])
        self.assertEqual(found["eight_signs"]["together"], 240)
        self.assertTrue(found["rest_of_three_signs"]["closed_under_product"])

    def test_completion_keeps_no_blind_chiral_label_once_partner_directions_are_probed(self):
        found = self.result["red_team_test"]["result"]
        self.assertEqual(found["labels_blind_in_all_eight_directions"], 0)
        self.assertGreater(found["labels_blind_in_the_five_directions_only"], 0)
        blind = lambda name: [found["family_labels"][name]["blind_with_k_partner_directions_probed"][k] for k in range(4)]
        self.assertEqual(blind("four_families"), [True, True, True, False])
        self.assertEqual(blind("three_families"), [True, False, False, False])
        self.assertEqual(found["family_labels"]["three_families"]["words_marked_plus"], 48)

    def test_six_is_the_only_number_with_a_second_description(self):
        self.assertEqual([o["n"] for o in symmetric() if o["clauses"]], [6])
        found = six()
        self.assertTrue(found["faithful"] and found["generalised_quadrangle_15_points_15_lines"])
        self.assertEqual(found["cycle_types_exchanged"]["(3, 2, 1)"], ["(6,)"])

    def test_nothing_is_promoted_and_the_lane_is_registered(self):
        self.assertFalse(self.result["physical_solution"])
        self.assertIn({"lane": "genesis_rules"}, [spec for _, spec in specifications(2, ".", None, None)])
        self.assertEqual(execute({"lane": "genesis_rules"}, ".")["status"], "GENESIS_RULES_SEARCHED")


class CurrentProbes(unittest.TestCase):
    def test_family_labels_under_root_currents(self):
        r = genesis_rules.current_probes()
        three, four = r["three_families"], r["four_families"]
        self.assertEqual(three["mixed"]["seen_at_first_order"], 60)
        self.assertEqual(three["mixed"]["first_order_values"], [-16, 16])
        self.assertEqual(three["partner"]["seen_at_first_order"], 6)
        self.assertEqual(three["blind_part"]["seen_at_third_order"], 0)
        self.assertEqual([k for k, v in three["sets_of_currents"].items() if v["label_blind"]], ["blind_part_only"])
        self.assertEqual(four["mixed"]["seen_at_first_order"]+four["mixed"]["seen_at_third_order"], 0)
        self.assertEqual(four["mixed"]["marked_words_kept_flipped_removed"], [(0, 32, 0)])
        self.assertTrue(four["sets_of_currents"]["blind_part_and_mixed_into_two_partner_directions"]["label_blind"])
        self.assertFalse(four["sets_of_currents"]["blind_part_and_partner"]["label_blind"])
        self.assertEqual(four["mixed_cubic_of_three_currents_into_three_partner_directions"], 128)


if __name__ == "__main__":
    unittest.main()
