import unittest
from fractions import Fraction as Fr

from autonomy import execute, specifications
from process_rules import flavour_pairs, matrix_kernels, run, schur_coefficient, seam_layers, spin2_kernels


class ProcessRules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run()

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
        ambiguous = [o["name"] for o in spin2_kernels() if {"dirac_low_energy", "kallen_lehmann"} <= o["clauses"]]
        self.assertEqual(len(ambiguous), 3)

    def test_real_seam_needs_spin10_level_one(self):
        rows = self.selected("seam_layers")
        self.assertIn(["c_eight", "sheet_neutral", "spin10_level_one"], rows["n=16 real reversed"]["minimal_rules"])
        self.assertIn(["c_eight", "sheet_neutral"], rows["n=8 complex reversed"]["minimal_rules"])
        self.assertEqual(len(seam_layers()), 128)

    def test_tri_bimaximal_and_family_clock_patterns(self):
        patterns = {(o["pair"], o["pattern"]) for o in flavour_pairs()}
        tbm = (("0", "1/3", "2/3"), ("1/2", "1/3", "1/6"), ("1/2", "1/3", "1/6"))
        self.assertIn(("Z3 x V4'", tbm), patterns)
        clock = (("2/45", "1/3", "28/45"), ("1/3", "1/3", "1/3"), ("28/45", "1/3", "2/45"))
        self.assertIn(("Z3 x Z4", clock), patterns)

    def test_whole_process_needs_eight_load_bearing_principles(self):
        comp = self.result["composed_process"]
        self.assertEqual(comp["longest_minimal_rule"], 8)
        chosen = {r["selected"] for r in comp["selected_by_the_longest_rules"]}
        self.assertEqual(len(chosen), 1)
        label = chosen.pop()
        for part in ("kernel nu=1 eps=0", "spin2 Schur w=1/2 m=5/3", "seam n=16 real reversed", "flavour Z3 x Z4"):
            self.assertIn(part, label)
        row = comp["selected_by_the_longest_rules"][0]
        self.assertTrue(all(v[1] != [label] for v in row["ablation"].values()))

    def test_nothing_is_promoted_and_the_lane_is_registered(self):
        self.assertFalse(self.result["physical_solution"])
        self.assertIn({"lane": "process_rules"}, [spec for _, spec in specifications(2, ".", None, None)])
        self.assertEqual(execute({"lane": "process_rules"}, ".")["status"], "PROCESS_RULES_SEARCHED")


if __name__ == "__main__":
    unittest.main()
