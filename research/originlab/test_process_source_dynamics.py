import json
import unittest

from process_source_dynamics import parent_dynamics_universe, source_dynamics_audit


class SourceDynamicsTests(unittest.TestCase):
    def test_existing_engine_needs_the_explicit_quadratic_rule(self):
        result = parent_dynamics_universe()
        row = next(r for r in result["unique_selections"] if r["selected"] == "quadratic CCR parent")
        self.assertEqual(row["minimal_rules"], [["quadratic_in_the_two_current_modes"]])
        common = {"positive_response", "same_KMS_purification", "one_ground_state", "quarter_symmetry"}
        matches = [o for o in result["candidates"] if common <= set(o["clauses"])]
        self.assertEqual(len(matches), 4)
        self.assertFalse(result["unqualified_unique_process"])
        self.assertTrue(all(not o["physical_solution"] for o in matches))

    def test_combined_certificate_is_strict_json_and_has_no_promotion(self):
        result = source_dynamics_audit()
        self.assertTrue(result["all_checks_pass"], result["checks"])
        json.dumps(result, allow_nan=False)
        self.assertFalse(result["physical_solution"])
        self.assertEqual(result["ledger_promotions"], [])


if __name__ == "__main__":
    unittest.main()
