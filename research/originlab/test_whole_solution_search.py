"""Tests for the whole-solution search (H6h50, z171)."""
import unittest

import whole_solution_search as w


class WholeSolutionSearch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = w.run_whole_solution_search()

    def test_ledger_root(self):
        self.assertAlmostEqual(w.alpha_inverse(48, 1.25, 4.1), 137.0359992168407, places=9)

    def test_family_group_is_s4(self):
        g = w.family_group()
        self.assertEqual(len(g), 24)
        self.assertEqual(sorted(len([e for e in g if e[2] == c]) for c in set(e[2] for e in g)), [1, 3, 6, 6, 8])

    def test_shells(self):
        sh = w.shells(2)
        self.assertEqual((len(sh[1]), len(sh[2])), (240, 2160))

    def test_world_hits_are_one_object(self):
        self.assertEqual(self.r["world_hits"], ["1 x fredholm | class 1 | family standard | kernel D5",
                                                "1 x trace | class 1 | family standard | kernel D5"])
        self.assertEqual(self.r["all_clauses"], ["1 x fredholm | class 1 | family standard | kernel D5"])
        d = self.r["derived"]["1 x fredholm | class 1 | family standard | kernel D5"]
        self.assertEqual((d["n_top"], d["b1"], d["base"], d["power"]), ("48", "41/10", "4/3", "5/4"))
        self.assertEqual(d["series"], ["0", "48", "360", "1072", "2484"])
        self.assertLess(abs(d["alpha_pull_sigma"]), 2)

    def test_nearest_miss_is_far(self):
        self.assertTrue(all(abs(v["alpha_pull_sigma"]) > 50 for v in self.r["nearest_miss"].values()))

    def test_no_solution_claim(self):
        self.assertFalse(self.r["physical_solution"])
        self.assertEqual(len(self.r["ledger_inputs_not_derived"]), 5)


class LawFormsAndCompositions(unittest.TestCase):
    def test_pairing_action(self):
        self.assertEqual(w.pairing_action((0, 2, 3, 1)), [[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        self.assertEqual(w.pairing_action((1, 2, 3, 0)), [[0, 0, 1], [0, -1, 0], [-1, 0, 0]])

    def test_law_form_is_unique(self):
        r = w.run_law_form_search()
        self.assertEqual(r["objects"], 1092)
        self.assertEqual(r["hits"], ["c3 on opposite | 3-cycle | b1 ln(1/phi_seam) | alpha = 2 c3^2 t"])
        self.assertTrue(all(abs(v) > 1e4 for v in r["nearest_miss"].values()))
        self.assertFalse(r["physical_solution"])

    def test_compositions(self):
        c = w.run_composition_search()
        self.assertEqual(c["lorentzian"], ["1 partner | hermitian (2, 2)"])
        self.assertEqual(c["euclidean"], ["1 partner | quaternionic"])
        self.assertEqual(c["signature_without_real_structure"], {"(4,0)": (6, 0), "(3,1)": (0, 0)})
        self.assertEqual(c["common_slice"], (1, 4))
        self.assertFalse(c["physical_solution"])


if __name__ == "__main__":
    unittest.main()
