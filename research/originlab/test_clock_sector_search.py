import unittest

import clock_sector_search as cs


class ClockSectorSearch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = cs.run_clock_sector_search()
        cls.uni = cls.res["universes"]

    def rule(self, universe, clauses):
        rows = [r for r in self._rules(universe) if set(r["minimal_rules"][0]) and
                any(set(m) == set(clauses) for m in r["minimal_rules"])]
        return rows[0]["selected"] if rows else None

    def _rules(self, universe):
        return self.uni[universe]["search"]["unique_selections"]

    def test_no_physical_solution_claimed(self):
        self.assertEqual(self.res["status"], "CLOCK_SECTORS_SEARCHED")
        self.assertFalse(self.res["physical_solution"])

    def test_clock_classes(self):
        c = self.uni["clocks"]
        self.assertEqual(c["classes_per_order"], {2: 2, 3: 4, 4: 7, 5: 14, 6: 20, 7: 38, 8: 53})
        zero = c["trace_zero_classes"]
        self.assertEqual(len(zero), 5)
        self.assertEqual([z for z in zero if z.startswith(("Z2", "Z3", "Z4", "Z5"))], [c["h6_clock"]])
        self.assertEqual(c["derived"][c["h6_clock"]], {"multiplicities": [60, 64, 60, 64], "twist_ground": "3/8",
                                                       "anomaly_type": 2})

    def test_clock_is_pinned_by_trace_zero(self):
        h6 = self.uni["clocks"]["h6_clock"]
        self.assertEqual(self.rule("clocks", ["order_at_most_4", "trace_zero"]), h6)
        self.assertEqual(self.rule("clocks", ["semisimple_fixed", "trace_zero"]), h6)
        clauses = self.uni["clocks"]["clauses"][h6]
        self.assertIn("square_is_fermion_parity", clauses)
        self.assertNotIn("anomaly_free", clauses)

    def test_completions(self):
        c = self.uni["completions"]
        self.assertEqual(c["candidates"], 27)
        self.assertEqual(c["counts"]["bosonic"], 10)
        self.assertEqual(c["counts"]["spin"], 17)
        graphs = [l for l, cl in c["clauses"].items() if {"graph", "bosonic"} <= set(cl)]
        self.assertEqual(len(graphs), 4)                          # the four sign isometries: two bits
        self.assertTrue(all({"scalar_16", "scalar_4", "mirror_symmetric"} <= set(c["clauses"][l]) for l in graphs))

    def test_e8_currents_exclude_a_local_16(self):
        c = self.uni["completions"]
        self.assertEqual(c["e8_currents_and_local_16"], [])
        for label, cl in c["clauses"].items():
            if "factorised" in cl and "left [(0, 0), (1, " in label:      # left chiral algebra is E8
                self.assertNotIn("clock_sectors_present", cl)

    def test_response_power(self):
        self.assertTrue(self.rule("responses", ["insertion_matches", "native"]).startswith("coherent kernel of the 16"))
        cl = self.uni["responses"]["clauses"]
        self.assertNotIn("native", cl["quartic measure on five real directions (analysis 22), F = q K"])
        self.assertNotIn("insertion_matches", cl["coherent kernel of the 16, F = q d log K/dq"])


if __name__ == "__main__":
    unittest.main()
