import unittest

import structure_objects

from atlas_scenarios import assess, knowledge
from autonomy import HERE, execute, specifications
from structure_scenarios import Search, atlas_fingerprint, run, stated_numbers, universe


def finding(problem, text, status="established_scoped", premises=("a",)):
    return {"id": text, "status": status, "premise_ids": list(premises), "problem": problem,
            "title": {"de": ""}, "statement": {"de": text}, "sources": []}


class StructureScenarios(unittest.TestCase):
    def test_numbers_exclude_identifiers_dates_and_decimals(self):
        found = stated_numbers("v724 und R113 am 2026-11-19: Rang 16, Anteil 5/4, 6/3, Wert 3,14 und 12 % sowie 240.")
        self.assertEqual(found, {("int", 16), ("frac", "5/4"), ("int", 240)})

    def test_bookkeeping_of_the_record_is_not_structure(self):
        found = stated_numbers("Residuen<3.2e−13, Radien 0 und .35, 26 exakte Kontrollen bestehen, mod 65537, "
                               "v28:55–75, auf 75 Dezimalstellen, 33/33 gezielte Kontrollen, Rang 16, 5/4.")
        self.assertEqual(found, {("int", 16), ("frac", "5/4")})

    def test_glue_classes_of_e8_are_computed_from_the_roots(self):
        from structure_objects import ROOTS, _e8_classes, subsystems
        classes = _e8_classes()
        full = {kind for kind in subsystems(("E", 8)) if sum(n for _, n in kind) == 8}
        self.assertEqual(set(classes) & full, full)
        self.assertTrue(all(sum(classes[kind]) == 240-sum(ROOTS(t, n) for t, n in kind) for kind in full))
        self.assertEqual(classes[(("A", 3), ("D", 5))], [60, 64, 64])
        self.assertEqual(classes[(("A", 2), ("E", 6))], [81, 81])
        self.assertEqual(classes[(("D", 4), ("D", 4))], [64, 64, 64])

    def test_nested_split_carries_charges_and_operators(self):
        pool = {row["id"]: row for row in universe()}
        feats = pool["nested:3+2"]["feats"]
        self.assertTrue({("frac", "1/6", "charge of an even Fock state"), ("frac", "2/3", "charge of an even Fock state"),
                         ("int", 6, "even Fock states of one charge"), ("int", 13, "operators that keep the split"),
                         ("int", 12, "traceless operators that keep the split")} <= feats)
        self.assertTrue({("int", 5, "two pairs: four distinct places"), ("int", 45, "two pairs: mixed symmetry"),
                         ("int", 50, "two pairs: window symmetry"), ("int", 24, "operators on pairs: adjoint part"),
                         ("int", 75, "operators on pairs: remaining part")} <= pool["graph:K(5,2)"]["feats"])

    def test_move_rule_gives_stay_hop_loss_wall_and_hand(self):
        rows = {row["id"]: row for row in structure_objects.gen_moves()}
        cube = rows["moves:3x2"]["feats"]
        for feature in (("int", 6, "moves"), ("int", 8, "states"), ("int", 4, "levels with the wall"), ("frac", "1/2", "probability to stay"),
                        ("frac", "1/6", "probability per neighbour"), ("frac", "1/3", "loss from a pair"), ("frac", "2/3", "rate of a level"),
                        ("frac", "64/729", "rate of a level after one hand"), ("frac", "1/729", "rate of a level after one hand"),
                        ("int", 729, "histories of one hand"), ("int", 6, "Majorana operators"), ("int", 15, "dimension of so(2d)")):
            self.assertIn(feature, cube)
        self.assertNotIn(("frac", "64/729", "rate of a level after one hand"), rows["moves:3x3"]["feats"])   # its hand has nine steps
        fermions = {row["id"]: row for row in structure_objects.gen_fermions()}
        self.assertEqual(fermions["fermion:3"]["parts"], ["moves:3x2"])

    def test_fingerprint_needs_independent_established_problems(self):
        kb = {"findings": [finding("p1", "Rang 16 und 5/4"), finding("p2", "wieder 16, dazu 5/4", premises=("b",)),
                           finding("p3", "nochmals 16"), finding("p3", "und 16"), finding("p4", "5/4", "excluded_scoped")]}
        self.assertEqual([(row["value"], row["problems"]) for row in atlas_fingerprint(kb)], [(16, 3)])
        self.assertEqual({row["value"] for row in atlas_fingerprint(kb, 2)}, {16, "5/4"})

    def test_pool_reproduces_known_facts(self):
        pool = {row["id"]: row for row in universe()}
        has = lambda name, *feature: any(f[:2] == feature for f in pool[name]["feats"])
        full = [row for name, row in pool.items() if name.startswith("lattice:") and name.endswith("<E8")
                and ("int", 8, "rank") in row["feats"]]
        self.assertEqual(len(full), 15)
        self.assertTrue(has("lattice:D5 + A3<E8", "type", "Z4") and has("lattice:D5 + A3<E8", "frac", "5/4"))
        self.assertTrue(("int", 46080, "order") in pool["reflection:G31"]["feats"])
        self.assertTrue(all(("int", 12, "Euler number") in row["feats"] for row in pool.values()
                            if row["fam"] == "Elliptic surface"))
        petersen = pool["graph:K(5,2)"]["feats"]
        self.assertTrue({("frac", "1/3", "normalised eigenvalue"), ("frac", "2/3", "normalised eigenvalue"),
                         ("int", 4, "multiplicity"), ("int", 5, "multiplicity"), ("int", 15, "edges")} <= petersen)
        self.assertTrue(has("fermion:8", "int", 240) and has("fermion:8", "int", 128) and not has("fermion:7", "int", 240))
        self.assertTrue(has("nested:3+2", "frac", "5/6") and has("nested:3+2", "frac", "10/3")
                        and has("nested:3+2", "frac", "5/3") and has("nested:(3+2)+3", "frac", "41/10"))

    def test_a_closure_holds_its_parts_and_their_parts(self):
        pool = {row["id"]: row for row in universe()}
        closure = pool["closure:nested:(3+2)+3"]
        self.assertTrue({"fermion:5+3", "fermion:3+2", "graph:K(5,2)", "lattice:D5 + A3<E8"} <= set(closure["parts"]))
        self.assertTrue(pool["lattice:D5 + A3<E8"]["feats"] <= closure["feats"])

    def test_ratios_are_optional_and_labelled(self):
        plain = {row["id"]: row for row in universe()}["rep:A4 = SU(5)"]["feats"]
        ratio = {row["id"]: row for row in universe(ratios=True)}["rep:A4 = SU(5)"]["feats"]
        self.assertTrue(plain < ratio and all(role.startswith("ratio ") for _, _, role in ratio-plain))

    def test_random_fingerprints_stay_inside_the_pool(self):
        import random
        search = Search(universe())
        fake = search.fake([("int", 240), ("frac", "5/4"), ("int", 16)], random.Random(3))
        self.assertTrue(all(key in search.have for key in fake) and len(set(fake)) == 3)
        self.assertNotIn(("int", 240), fake)

    def test_live_atlas_run_is_repeatable_and_its_scenarios_are_admissible(self):
        root = HERE.parents[1]
        self.assertIn({"lane": "structure_scenarios"}, [job for _, job in specifications(3, root, None, [])])
        first, again = run(root, fakes=20, splits=4, seed=5), run(root, fakes=20, splits=4, seed=5)
        self.assertEqual(first, again)
        live = execute({"lane": "structure_scenarios"}, root)
        self.assertEqual(live["status"], "STRUCTURE_SCENARIOS_RANKED")
        self.assertFalse(live["physical_solution"])
        self.assertTrue(live["proposals"] and all(0 < row["p_search_wide"] <= 1 for row in live["singles"]))
        self.assertTrue(0 <= live["margin_over_runner_up"]["first_ahead_share"] <= 1)
        kb = knowledge(root)
        for scenario in live["proposals"]:
            verdict = assess(kb, scenario)
            self.assertEqual(verdict["status"], "ADMISSIBLE_FOR_CHECKS")
            self.assertFalse(verdict["physical_solution"])

    def test_missing_atlas_is_reported(self):
        self.assertEqual(run(HERE)["status"], "ATLAS_UNAVAILABLE")


if __name__ == "__main__":
    unittest.main()
