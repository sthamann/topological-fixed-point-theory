import unittest

from autonomy import execute, specifications
from foreign_universes import TARGETS, automata, graphs, hits, polytopes, run, score


class ForeignUniverses(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run()

    def test_enumerations_have_the_known_sizes(self):
        self.assertEqual([len(graphs(n)) for n in (3, 4, 5)], [3, 9, 30])
        self.assertEqual(len(graphs(6)), 142)
        self.assertEqual(len(automata()), 88)
        self.assertEqual(len(polytopes()), 33)

    def test_petersen_like_checks_on_small_graphs(self):
        five_cycle = [g for g in graphs(5) if g["edges"] == 5 and "regular" in g["clauses"]]
        self.assertEqual(len(five_cycle), 1)
        self.assertIn("self_complementary", five_cycle[0]["clauses"])
        self.assertIn("vertex_transitive", five_cycle[0]["clauses"])
        self.assertEqual(five_cycle[0]["invariants"]["automorphisms"], 10)
        self.assertEqual(five_cycle[0]["invariants"]["spanning_trees"], 5)

    def test_automata_clauses(self):
        rules = {a["rule"]: a for a in automata()}
        self.assertIn("additive", rules[90]["clauses"])
        self.assertIn("additive", rules[150]["clauses"])
        self.assertIn("conserves_number", rules[170]["clauses"])
        self.assertIn("reversible_on_rings", rules[15]["clauses"])
        self.assertNotIn("reversible_on_rings", rules[30]["clauses"])

    def test_polytopes_have_euler_characteristic(self):
        for p in polytopes():
            if p["dim"] == 3:
                v, e, f = (p["invariants"][k] for k in ("vertices", "edges", "facets"))
                self.assertEqual(v-e+f, 2, p["name"])
        names = {p["name"]: p for p in polytopes()}
        self.assertEqual(names["24-cell"]["invariants"]["symmetry_order"], 1152)
        self.assertEqual(names["5-orthoplex"]["invariants"]["vertices"], 10)

    def test_trivial_target_does_not_beat_its_baseline(self):
        rows = score(polytopes(), lambda o: o["name"])
        self.assertGreater(rows["koide_q"]["objects_hit"], 0)
        self.assertFalse(rows["koide_q"]["beats_baseline"])

    def test_hits_use_only_the_grammar(self):
        self.assertIn("koide_q", hits({"a": 2, "b": 3}))
        self.assertEqual(hits({"a": 1}), {})
        self.assertEqual(set(self.result["targets"]), set(TARGETS))

    def test_run_is_scoped_and_registered(self):
        self.assertFalse(self.result["physical_solution"])
        self.assertEqual(self.result["objects_enumerated"], 142+88+33)
        lanes = [spec["lane"] for _, spec in specifications(3, ".", "x", None)]
        self.assertIn("genesis_foreign", lanes)
        self.assertEqual(execute({"lane": "genesis_foreign"}, None)["status"], "FOREIGN_UNIVERSES_SEARCHED")


if __name__ == "__main__":
    unittest.main()
