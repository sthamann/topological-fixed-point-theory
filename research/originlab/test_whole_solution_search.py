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


class SeedWorlds(unittest.TestCase):
    """Seed worlds, readings and assemblies (exploration 2026-11-27)."""
    @classmethod
    def setUpClass(cls):
        cls.objs = w.worlds()
        cls.rules = w.world_rule_search(cls.objs)
        cls.assemblies = w.run_assembly_search(cls.objs)

    def test_universe_and_ledger_numbers(self):
        self.assertEqual(len(self.objs), 16464)
        world = next(o for o in self.objs if o["seed"] == w.TFPT_WORLD)
        self.assertEqual(world["clauses"], set(w.PRINCIPLES))                 # every principle, no rival
        self.assertAlmostEqual(world["numbers"]["alpha_inverse"], w.AINV, places=9)
        self.assertAlmostEqual(world["numbers"]["b1"], 4.1, places=12)
        self.assertEqual(w.positive_splits(5), [(3, 2)])
        self.assertEqual([N for N in range(2, 30) if len(w.positive_splits(N)) == 1], [5])

    def test_rule_search_can_select_other_worlds(self):
        r = self.rules
        self.assertEqual(r["rules"], 83682)
        self.assertEqual(r["verdicts"], {"EMPTY": 53660, "UNIQUE": 29587, "AMBIGUOUS": 435})
        self.assertEqual(r["rules_selecting_target"], 6692)
        self.assertEqual(r["minimal_rules_by_length"], {3: 60, 4: 287, 5: 52})
        self.assertEqual(r["uniquely_selected_worlds"][0], ("k2 n4 N5 f3 d4", 6692))
        self.assertEqual(r["uniquely_selected_worlds"][1][0], "k2 n4 N5 f3 d2")
        self.assertIn({"rule": ["flat_cover", "spin_sheets", "unimodular_glue"],
                       "without": {"flat_cover": "k1 n2 N5 f3 d2", "spin_sheets": "AMBIGUOUS",
                                   "unimodular_glue": "k2 n4 N2 f1 d4"}}, r["minimal_rules"])

    def test_bitmask_search_agrees_with_genesis_search(self):
        from genesis_rules import search
        small = w.worlds({"k": (1, 3), "n": (2, 4), "N": (4, 5), "f": (2, 3), "d": (3, 4)})
        pool = ["flat_cover", "spin_sheets", "unimodular_glue", "three_sheets", "anyonic_space", "odd_marks"]
        slow = search(small, pool, lambda o: o["size"], lambda o: o["label"])
        fast = w.world_rule_search(small, pool, max_len=len(pool))
        self.assertEqual(fast["verdicts"], slow["verdicts"])
        chosen = {tuple(m) for row in slow["unique_selections"] if row["selected"] == "k2 n4 N5 f3 d4"
                  for m in row["minimal_rules"]}
        self.assertEqual({tuple(row["rule"]) for row in fast["minimal_rules"]}, chosen)

    def test_only_one_seed_closes(self):
        c = w.coherent_worlds(self.objs)
        self.assertEqual(c["closing_seeds_without_dimension"], [[2, 4, 5, 3]])
        self.assertEqual(c["closing_with_spacetime_clauses"], ["k2 n4 N5 f3 d4"])
        self.assertEqual(c["dimensions_allowed_by_closure_alone"], [2, 3, 4, 5, 6, 7, 8])
        self.assertEqual(c["flat_divisor_worlds_and_failures"]["k3 n3 N4 f2"],
                         ["positive_split", "type_a_glue", "unimodular_glue", "z4_anomaly"])

    def test_number_baseline(self):
        n = w.number_baseline(self.objs)
        self.assertEqual((n["seeds"], n["with_numbers"]), (336, 210))
        self.assertEqual(n["hits_histogram"], {0: 276, 1: 25, 2: 30, 3: 4, 4: 1})
        self.assertEqual(n["per_target"]["alpha_measured"], ["k2 N5 f3"])
        self.assertEqual(n["neighbours"]["k3 N5 f3"]["hits"], ["baryon", "cabibbo", "solar"])   # alpha decides
        self.assertAlmostEqual(n["neighbours"]["k3 N5 f3"]["alpha_inverse"], 308.866873, places=5)

    def test_assemblies(self):
        a = self.assemblies
        self.assertEqual((a["readings"], a["assemblies"], a["contradictions"]), (10, 386, 0))
        self.assertEqual(a["complete"], [])
        self.assertEqual(a["best"][0]["members"], ["anomaly_minimal", "background_free", "real_form"])
        self.assertEqual(a["best"][0]["distance"], 7)
        self.assertEqual(a["roles_no_reading_delivers"],
                         ["rp", "time", "continuum", "einstein", "area", "flavour_rest"])
        self.assertEqual([m for m, row in a["single_readings"].items() if row["seed_forced"]], ["sky_spin2"])
        self.assertEqual(a["smallest_seed_forcing"][0], ["sky_spin2"])
        self.assertTrue(all(v == "forced" for k, v in a["all_readings"]["roles"].items() if k in w.COMPUTED_ROLES))

    def test_complete_solution_is_representable(self):
        audit = w.completeness_audit(self.objs)
        self.assertEqual(audit["status"], "COMPLETE_SOLUTION_REPRESENTABLE")
        self.assertTrue(audit["planted_complete"])
        self.assertEqual(audit["planted_without_time"], ["time"])
        self.assertIn("marks4", audit["planted_without_marks_clause"])
        self.assertEqual(audit["planted_with_unbridged_partner"], 1)
        self.assertEqual((audit["roles"], audit["present_complete"], audit["present_best_distance"]), (24, [], 7))

    def test_ladder_gates(self):
        gates = w.ladder_gates()
        row = gates["so(4,2)"]
        self.assertTrue(row["kernel_positive"] and row["time_positive"] and row["time_continuous"])
        self.assertEqual([gates[f"so({d},2)"]["spacetime_dimension"] for d in range(3, 8)], [3, 4, 5, 6, 7])
        self.assertFalse(gates["so(6) control"]["time_continuous"])

    def test_leverage(self):
        lev = w.run_leverage(self.objs)
        self.assertEqual((lev["distance_now"], lev["lowest_distance"], lev["statements_needed"]), (7, 1, 3))
        self.assertEqual(lev["single"]["ladder"]["gain"], 3)
        self.assertEqual(lev["smallest_sets"], [{"statements": ["equilibrium", "irreducible", "ladder"],
                                                 "members": ["background_free", "real_form"],
                                                 "open": ["flavour_rest"]}])


class ShapeSearch(unittest.TestCase):
    """Prediction roles, ledger reading of the twistor bridge and Z4 gates (second pass of 2026-11-27)."""

    @classmethod
    def setUpClass(cls):
        cls.result = w.run_shape_search()

    def test_gates(self):
        g = self.result["gates"]
        self.assertEqual((g["roots"], g["classes"], g["compact_dimension"]),
                         (240, {0: 52, 1: 64, 2: 60, 3: 64}, 136))
        for key in ("centres_identified", "half_period_is_theta_times_clock", "full_period_is_fermion_parity",
                    "one_helicity_per_16", "mirror_excluded", "orientation_swaps", "z4_trivial_on_photons"):
            self.assertTrue(g[key], key)
        self.assertFalse(g["helicity_u1_in_e8"])
        self.assertEqual(g["rule"][1], {0: ["-1", "1"], 1: ["-1/2"], 2: ["0"], 3: ["1/2"]})

    def test_family_clock_gates(self):
        f = self.result["family_gates"]
        self.assertTrue(f["clock_is_coxeter_element"] and f["squares_are_the_three_pairings"])
        self.assertEqual((f["clock_square_signature"], f["real_form_selected"], f["real_rank"]),
                         ((1, 2), ["su(2,2)"], 2))
        self.assertFalse(f["outer_forms_in_weyl_group"])
        self.assertEqual(f["centraliser_order"], 8)
        self.assertTrue(f["restricted_system_is_c2"])
        self.assertEqual((f["ad_spectrum_twistor_clock"], f["ad_spectrum_family_clock"]), ((7, 4, 0, 4), (3, 4, 4, 4)))
        self.assertFalse(f["clocks_conjugate"])
        self.assertTrue(f["clocks_share_square"])
        self.assertEqual(f["e8_compact_dimension"], 136)

    def test_contract_has_prediction_roles(self):
        r = self.result
        self.assertEqual((r["roles"], r["assemblies"]["readings"], r["assemblies"]["complete"]), (26, 11, []))
        self.assertEqual(r["assemblies"]["best"][0]["distance"], 9)
        self.assertEqual(r["assemblies"]["roles_no_reading_delivers"][-2:], ["chirality", "parity_response"])
        self.assertTrue(r["planted_complete"])
        self.assertEqual(r["planted_without_parity_response"], ["parity_response"])

    def test_leverage_names_five_statements(self):
        lev = self.result["leverage"]
        self.assertEqual((lev["distance_now"], lev["lowest_distance"], lev["statements_needed"]), (9, 1, 5))
        self.assertEqual(lev["smallest_sets"][0]["statements"],
                         ["equilibrium", "family_is_twistor", "irreducible", "ladder", "oriented_state"])
        self.assertEqual(lev["smallest_sets"][0]["open"], ["flavour_rest"])
        self.assertEqual([t["statements"] for t in self.result["tensions"]],
                         [["ladder", "oriented_state"], ["family_is_twistor"]])

    def test_section_d_contract_is_unchanged(self):
        self.assertEqual((len(w.ROLES), len(w.READINGS)), (24, 10))


if __name__ == "__main__":
    unittest.main()
