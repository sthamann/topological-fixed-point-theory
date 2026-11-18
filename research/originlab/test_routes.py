import unittest
from copy import deepcopy

from routes import examples, examine, ROUTES


class SixRoutesTests(unittest.TestCase):
    def test_scope_combinations_and_variation(self):
        candidates = list(examples())
        self.assertEqual(len(candidates),32)
        good = examine(candidates[0])
        self.assertEqual(set(good["routes"]),set(ROUTES))
        self.assertEqual(len(good["combinations"]),57)
        self.assertEqual(good["joint_markov_moments"]["status"],"PASS_SCOPED")
        self.assertFalse(good["physical_solution"])
        other = examine(candidates[16])
        self.assertEqual(other["routes"]["markov"]["status"],"REFUTED_SCOPED")
        self.assertEqual(other["routes"]["moments"]["status"],"PASS_SCOPED")

    def test_cross_source_and_wrong_products(self):
        candidate = next(examples())
        bad = deepcopy(candidate)
        bad["routes"]["regional"]["source_id"] = "different"
        with self.assertRaisesRegex(ValueError,"cross-source"):
            examine(bad)
        bad = deepcopy(candidate)
        bad["markov_moment_binding"]["generators"]["X"] = bad["markov_moment_binding"]["word_matrices"][0]
        self.assertEqual(examine(bad)["joint_markov_moments"]["status"],"REFUTED_SCOPED")

    def test_absent_source_data_remain_open(self):
        report = examine({"source_id":"original","routes":{}})
        self.assertTrue(all(row["status"]=="SOURCE_REQUIRED" for row in report["routes"].values()))


if __name__ == "__main__":
    unittest.main()
