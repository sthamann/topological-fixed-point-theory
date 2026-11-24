import unittest
from fractions import Fraction as Fr

from autonomy import execute, specifications
import numpy as np

from process_rules import (BETA_DET, BETA_V, Q_DIRAC, Q_ORIGINAL, edge_form_beta, flavour_groups, flavour_pairs, hardy_rect, metric_star, ladder_occupation, matrix_kernels, original_current_norm,
                           q_rules, readouts, run, schur_coefficient, seam_layers, spin2_kernels, transfer_processes)


class ProcessRules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run()
        cls.h6h23 = run(transfer=False)
        cls.h6h24 = run(state=False)
        cls.h6h26 = run(metric=False)

    def selected(self, universe):
        return {row["selected"]: row for row in self.result["universes"][universe]["unique_selections"]}

    def test_schur_coefficients_follow_the_wallach_set(self):
        self.assertEqual(schur_coefficient(Fr(1), (3, 0)), 1)
        self.assertEqual(schur_coefficient(Fr(1), (1, 1)), 0)
        self.assertEqual(schur_coefficient(Fr(1, 2), (1, 1)), Fr(-1, 8))
        self.assertEqual(schur_coefficient(Fr(2), (1, 1)), 1)
        positive = {o["nu"] for o in matrix_kernels() if o["eps"] == 0 and "positive" in o["clauses"]}
        self.assertEqual(positive, {Fr(0), Fr(1), Fr(3, 2), Fr(2), Fr(3)})

    def test_hardy_slice_alone_leaves_the_det_term_open(self):
        rows = self.selected("matrix_kernels")["nu=1 eps=0"]
        self.assertIn(["hardy_slice", "wave_equation"], rows["minimal_rules"])
        self.assertIn(["conformally_covariant", "hardy_slice"], rows["minimal_rules"])
        hardy = [o for o in matrix_kernels() if "hardy_slice" in o["clauses"]]
        self.assertEqual({(o["nu"], o["eps"]) for o in hardy}, {(Fr(1), 0), (Fr(1), 1)})

    def test_spectral_action_kernels_are_never_positive(self):
        kernels = {o["name"]: o for o in spin2_kernels()}
        for name in ("spectral action e^-u", "spectral action step", "spectral action (1+u)^-3"):
            self.assertNotIn("no_euclidean_zero", kernels[name]["clauses"])
        self.assertNotIn("kallen_lehmann", kernels["ledger ansatz H = x e^x"]["clauses"])
        self.assertIn("kallen_lehmann", kernels["Schur w=1/2 m=5/3"]["clauses"])

    def test_positive_dirac_kernel_needs_the_mirror_weight(self):
        rows = self.selected("spin2_kernels")
        self.assertIn(["dirac_low_energy", "mirror_weight"], rows["Schur w=1/2 m=5/3"]["minimal_rules"])
        ambiguous = [o["name"] for o in spin2_kernels(False) if {"dirac_low_energy", "kallen_lehmann"} <= o["clauses"]]
        self.assertEqual(len(ambiguous), 3)
        ambiguous = [o["name"] for o in spin2_kernels() if {"dirac_low_energy", "kallen_lehmann"} <= o["clauses"]]
        self.assertEqual(len(ambiguous), 4)

    def test_transfer_memory_keeps_the_dirac_correction_only_at_q_dirac(self):
        self.assertAlmostEqual(Q_DIRAC**2/(1 - Q_DIRAC**2)**2, 0.3, places=12)
        memories = {o["name"]: o for o in spin2_kernels(state=False) if "transfer_memory" in o["clauses"]}
        self.assertEqual(len(memories), 5)
        self.assertEqual(len([o for o in spin2_kernels() if "transfer_memory" in o["clauses"]]), 7)
        self.assertTrue(all("kallen_lehmann" in o["clauses"] for o in memories.values()))
        self.assertEqual([n for n, o in memories.items() if "dirac_low_energy" in o["clauses"]], ["transfer memory q=q_dirac"])
        rows = self.selected("spin2_kernels")
        self.assertIn(["dirac_low_energy", "transfer_memory"], rows["transfer memory q=q_dirac"]["minimal_rules"])

    def test_transfer_structure_fixes_the_quarter_turn_but_not_q(self):
        tp = transfer_processes()
        quarter = {(o["d"], o["k"]) for o in tp if {"quarter_covariant", "kernel_one_grade"} <= o["clauses"]}
        self.assertEqual(quarter, {(1, 4), (1, 8)})
        self.assertTrue(all(o["kernel_dim"] == o["d"] for o in tp))
        one_grade = {(o["d"], o["k"]) for o in tp if "kernel_one_grade" in o["clauses"]}
        self.assertEqual(one_grade, {(1, 4), (1, 8)})
        self.assertEqual(self.result["universes"]["transfer_processes"]["verdicts"]["UNIQUE"], 0)

    def test_no_rank_one_readout_is_massless_and_connected_or_universal(self):
        self.assertEqual(len(readouts()), 10)
        empty = self.h6h26["universes"]["readouts"]["minimal_empty_rules"]   # scoped statement of z45, corrected in z50
        self.assertEqual(empty, [["connected", "massless_pole"], ["universal_residue"]])

    def test_real_seam_needs_spin10_level_one(self):
        rows = self.selected("seam_layers")
        self.assertIn(["c_eight", "sheet_neutral", "spin10_level_one"], rows["n=16 real reversed"]["minimal_rules"])
        self.assertIn(["c_eight", "sheet_neutral"], rows["n=8 complex reversed"]["minimal_rules"])
        self.assertEqual(len(seam_layers()), 128)

    def test_tri_bimaximal_and_family_clock_patterns(self):
        patterns = {(o["pair"], o["pattern"]) for o in flavour_pairs()}
        tbm = (("0", "1/3", "2/3"), ("1/2", "1/3", "1/6"), ("1/2", "1/3", "1/6"))
        self.assertIn(("Z3 x V4'", tbm), patterns)
        clock = (("1/3 - sqrt(3)/6", "1/3", "sqrt(3)/6 + 1/3"), ("1/3", "1/3", "1/3"), ("sqrt(3)/6 + 1/3", "1/3", "1/3 - sqrt(3)/6"))
        self.assertIn(("Z3 x Z4", clock), patterns)           # (2 - sqrt 3)/6, printed as 2/45 before z46

    def test_whole_process_needs_eight_load_bearing_principles(self):
        comp = self.h6h23["composed_process"]
        self.assertEqual(comp["longest_minimal_rule"], 8)
        chosen = {r["selected"] for r in comp["selected_by_the_longest_rules"]}
        self.assertEqual(len(chosen), 1)
        label = chosen.pop()
        for part in ("kernel nu=1 eps=0", "spin2 Schur w=1/2 m=5/3", "seam n=16 real reversed", "flavour Z3 x Z4"):
            self.assertIn(part, label)
        row = comp["selected_by_the_longest_rules"][0]
        self.assertTrue(all(v[1] != [label] for v in row["ablation"].values()))

    def test_one_transfer_process_needs_nine_load_bearing_principles(self):
        comp = self.result["composed_process"]
        self.assertEqual(comp["longest_minimal_rule"], 9)
        chosen = {r["selected"] for r in comp["selected_by_the_longest_rules"]}
        self.assertEqual(len(chosen), 1)
        label = chosen.pop()
        for part in ("kernel nu=1 eps=0", "spin2 transfer memory q=q_dirac", "seam n=16 real reversed", "flavour Z3 x Z4"):
            self.assertIn(part, label)
        row = comp["selected_by_the_longest_rules"][0]
        self.assertIn("one_transfer", row["minimal_rules"][0])
        self.assertNotIn("mirror_sheet", row["minimal_rules"][0])
        self.assertTrue(all(v[1] != [label] for v in row["ablation"].values()))

    def test_earlier_stage_is_reproduced(self):
        self.assertEqual(self.h6h24["composed_process"]["objects"], 67200)
        self.assertEqual(self.h6h24["composed_process"]["longest_minimal_rule"], 9)
        self.assertNotIn("q_rules", self.h6h24["universes"])

    def test_original_state_fixes_q(self):
        omega = original_current_norm()
        self.assertAlmostEqual(omega/(1/(np.exp(4*np.pi) - 1)), 1.0, places=10)
        self.assertAlmostEqual(original_current_norm(alpha=0.75)/omega, 1.0, places=10)
        self.assertAlmostEqual(ladder_occupation(Q_ORIGINAL)/omega, 1.0, places=8)
        rows = {r["selected"]: r for r in self.result["universes"]["q_rules"]["unique_selections"]}
        self.assertEqual(rows["e^-2pi"]["minimal_rules"], [["original_state"]])
        self.assertEqual(rows["0 (Hardy state)"]["minimal_rules"], [["hardy_state"]])
        empty = self.result["universes"]["q_rules"]["minimal_empty_rules"]
        for pair in (["dirac_same_scale", "original_state"], ["equal_weight", "original_state"], ["hardy_state", "memory"]):
            self.assertIn(pair, empty)
        self.assertEqual(len(q_rules()), 7)

    def test_scale_free_search_selects_the_original_state_process(self):
        free = self.result["composed_process"]["without_scale_convention"]
        self.assertEqual(free["longest_minimal_rule"], 8)
        chosen = {r["selected"] for r in free["selected_by_the_longest_rules"]}
        self.assertEqual(len(chosen), 1)
        label = chosen.pop()
        for part in ("kernel nu=1 eps=0", "spin2 transfer memory q=e^-2pi", "seam n=16 real reversed", "flavour Z3 x Z4"):
            self.assertIn(part, label)
        ties = free["selected_by_the_longest_rules"][0]["tie_sets_without_size"]
        sizes = sorted(len(v) for v in ties.values())
        self.assertEqual(sizes, [1, 2])                   # with e8_seam no size tie; with positivity n = 16 and 32 tie
        two = [v for v in ties.values() if len(v) == 2][0]
        self.assertTrue(any("seam n=32 real reversed" in t for t in two))

    def test_delta_groups_give_only_tm2_and_no_data_fit(self):
        fg = flavour_groups(ns=(2, 4, 8))
        self.assertFalse(any("clock_pattern" in o["clauses"] for o in fg if o["n"] == 2))
        self.assertTrue(any("clock_pattern" in o["clauses"] for o in fg if o["n"] == 4))
        mixed = [o for o in fg if "mixes" in o["clauses"]]
        self.assertTrue(all("tm2" in o["clauses"] and "tm1" not in o["clauses"] for o in mixed))
        self.assertTrue(any("0.011358058" in o["pattern"][0] for o in mixed))   # no spurious exact label
        empty = self.result["universes"]["flavour_groups"]["minimal_empty_rules"]
        for rule in (["fits_data"], ["theta13_window"], ["mixes", "tm1"]):
            self.assertIn(rule, empty)

    def test_metric_star_is_the_quarter_transfer(self):
        for qv in (0.5, Q_ORIGINAL):
            st, A = metric_star(qv, 5), hardy_rect(qv, 5)
            self.assertTrue(np.allclose(np.sort(st["sv"]), np.sort(np.linalg.svd(A, compute_uv=False)), atol=1e-12))
            self.assertEqual(st["nullity"], 1)
            ratio = st["eta"]/st["M"]
            self.assertLess(np.ptp(ratio)/np.mean(ratio), 1e-12)
            self.assertAlmostEqual(st["svd_kernel_overlap"], 1.0, places=10)
        self.assertTrue(all(v["singular_values_agree"] for v in self.result["metric_star"].values()))

    def test_only_the_quartic_number_is_free_with_undirected_edges(self):
        pots = self.result["universes"]["metric_potentials"]
        self.assertEqual(pots["tie_flat_mass_undirected"], ["V", "V + det(I-X)/2", "V - det(I-X)/2"])
        d = pots["derived"]
        self.assertIn("cubic_as_V", d["V + det(I-X)/2"]["clauses"])
        self.assertNotIn("quartic_as_V", d["V + det(I-X)/2"]["clauses"])
        self.assertNotIn("cubic_as_V", d["V + cubic/2"]["clauses"])
        self.assertNotIn("undirected_edge", d["V + cubic/2"]["clauses"])

    def test_field_covariance_is_connected_and_massless(self):
        empty = self.result["universes"]["readouts"]["minimal_empty_rules"]
        self.assertNotIn(["field_connected", "massless_pole"], empty)
        self.assertNotIn(["universal_residue"], empty)
        self.assertTrue(all("vacuum_subtracted" in r for r in empty))
        self.assertNotIn("metric_potentials", self.h6h26["universes"])
        self.assertEqual(self.h6h26["universes"]["readouts"]["minimal_empty_rules"],
                         self.h6h24["universes"]["readouts"]["minimal_empty_rules"])

    def test_quartic_number_is_the_shape_of_the_endpoint_form(self):
        for x, y in ((1, 1), (1, 2), (0, 1), (1, 0), (0.3, -1.1)):
            beta = np.array(edge_form_beta(x, y))
            t = 2*(beta[0] - BETA_V[0])
            self.assertTrue(np.allclose(beta, np.array(BETA_V) + t/2*np.array(BETA_DET), atol=1e-12))
        forms = self.result["universes"]["edge_forms"]["derived"]
        ts = {k.split(" (")[0]: round(v["t"], 9) for k, v in forms.items()}
        self.assertEqual(sorted(ts.values()), [-2.0, 0.0, 1.0, 4.0])
        low = [k for k, v in forms.items() if v["edge_degree"] <= 3]
        self.assertEqual(low, ["Q = e_p e_k (one vielbein of each endpoint)"])

    def test_dirac_residuals_give_no_tm1_either(self):
        fg = flavour_groups(ns=(2, 4), neutrinos="cyclic")
        self.assertFalse(any("tm1" in o["clauses"] for o in fg if "mixes" in o["clauses"]))
        self.assertFalse(any("fits_data" in o["clauses"] for o in fg))

    def test_nothing_is_promoted_and_the_lane_is_registered(self):
        self.assertFalse(self.result["physical_solution"])
        self.assertIn({"lane": "process_rules"}, [spec for _, spec in specifications(2, ".", None, None)])
        self.assertEqual(execute({"lane": "process_rules"}, ".")["status"], "PROCESS_RULES_SEARCHED")


if __name__ == "__main__":
    unittest.main()
