import json
import time
import unittest
from fractions import Fraction as Fr
from math import gcd

import sympy as sp

from reconstruction_search import (ETA_SPINOR, QUARTER_PAIR, SEAM_DOUBLET, SPIN_BRIDGE, _transport,
                                   PREMISE_MEANING, scenario_candidates,
                                   CHARACTERS, DOUBLE_COVER, M_VALUES, QUARTERS, U, Z, cover_genus, cov,
                                   doublet_state_candidates, exceptional_orbits, fock_multiplicities, fock_patterns,
                                   invariants, mark_symmetries, normal_bundle_candidates, original_symmetries, pushforward,
                                   run_reconstruction, spin1_matrix, spin2_family_candidates, symmetry_structure, v4_pencil,
                                   _weight)

TOMITA, EULER = "O(1) x Vbar (Tomita partner)", "O(1) x H^0(O(1)) (Euler bundle)"
POLES = "Hecke O(-1) -> O(-1)(D), D = {0,oo}"
RHO = {1: "diag(1, exp(-pi))/trace", 2: "diag(1, exp(-2*pi))/trace", 3: "diag(1, exp(-3*pi))/trace",
       4: "diag(1, exp(-4*pi))/trace"}


class ReconstructionSearch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        start = time.perf_counter()
        cls.result = run_reconstruction()
        cls.seconds = time.perf_counter() - start
        cls.bundles = {o["name"]: o for o in normal_bundle_candidates()}
        cls.states = doublet_state_candidates()

    def universe(self, name):
        return self.result["universes"][name]

    # ------------------------------------------------------------ doublet identifications (H6h40, z79)
    def test_transport_on_the_eta_doublet_is_exact(self):
        self.assertEqual(_transport(0, -1), sp.Matrix([[0, -1], [1, 0]]))       # spinor lift: square -1, det 1
        self.assertEqual(_transport(0, 1), sp.Matrix([[0, 1], [1, 0]]))
        self.assertEqual(_transport(0, 0), sp.diag(1, -1))                       # deck (z64: not scalar)

    def test_quarter_pair_is_the_eta_doublet_and_not_the_seam_doublet(self):
        q = self.universe("doublet_identifications")["questions"]
        self.assertEqual(q["r_equivalent_to_seam"]["survivors"], [SEAM_DOUBLET])
        self.assertEqual(sorted(q["spinor+r_equivalent_to_quarter_pair"]["survivors"]), sorted([ETA_SPINOR, QUARTER_PAIR]))
        d = self.universe("doublet_identifications")["derived"]
        self.assertEqual(d[SEAM_DOUBLET]["turns_of_16"], [0, 4])                 # diag(1, i): projective order 4
        self.assertIn("projective_order_1", d["quarter oscillators (a_1, a_2), every mode m"]["clauses"])

    def test_spin_bridge_positions_are_euclidean_and_no_native_marked_space_is_lorentzian(self):
        u = self.universe("position_spaces")
        self.assertEqual(u["derived"][SPIN_BRIDGE]["turns_of_16"], [0, 4, 8, 12])  # regular representation of mu_4
        self.assertTrue(u["derived"][SPIN_BRIDGE]["factorises_as_vstar_x_eta"])
        self.assertEqual(u["questions"]["native_r+lorentz_compatible+marked"]["verdict"], "EMPTY")
        self.assertEqual(u["questions"]["native_r+euclidean_compatible+marked"]["survivors"], [SPIN_BRIDGE])
        self.assertTrue(self.result["answers"]["quarter_pair_is_the_eta_doublet_not_the_seam_doublet"])

    # ------------------------------------------------------------ whole-theory scenarios (H6h42, S7)
    def test_scenarios_enumerate_all_axes_and_respect_recorded_exclusions(self):
        u = self.universe("scenarios")
        self.assertEqual(u["objects"], 4*3*4*2*10*2*5*6)                 # + pair-summed glue (H6h50, z174)
        sc, info = scenario_candidates()
        for o in sc:
            if (o["scenario"]["time"] == "modular flow" or o["scenario"]["spinor_index"] == "quarter oscillators"
                    or o["scenario"]["gravity"] in ("composite spin 2", "Virasoro T(z)",
                                                    "entanglement: original|partner cut", "moved algebra at fixed state")
                    or o["scenario"]["interaction"] == "transfer exchange (twirl)"):
                self.assertNotIn("admissible", o["clauses"])
        self.assertEqual(info["r_order"], {"quarter oscillators": 1, "step type": 1, "eta doublet": 2, "seam doublet V_1": 4})
        self.assertEqual(info["native_lorentz"], {"quarter oscillators": True, "step type": True, "eta doublet": False,
                                                  "seam doublet V_1": False})

    def test_scenarios_leave_graviton_open_and_modular_statistics_to_the_pair_sum(self):
        u = self.universe("scenarios")
        # H6h47 left graviton and chirality open; with the twistor-line index (z114) chirality is conditional under
        # P_GSO, and the blocker of H6h49 was the modular consistency of the statistics rule (z113). z174 ran that test
        # for the Z4 glue: anomalous per seam, modular and positive summed over a pair of seams (P_paar)
        self.assertEqual(u["minimal_empty_rules"], [["graviton"]])
        matter = "admissible+coherent+lorentzian_native+positive_time+one_field+spin_statistics+interacting"
        for name in (matter, matter + "+chiral_matter"):
            q = u["questions"][name]
            self.assertEqual((q["verdict"], len(q["survivors"])), ("AMBIGUOUS", 162))
        self.assertEqual(u["questions"]["admissible+chiral_matter"]["verdict"], "AMBIGUOUS")
        self.assertEqual(u["questions"]["admissible+modular_statistics"]["verdict"], "AMBIGUOUS")
        self.assertEqual(u["questions"]["admissible+graviton"]["verdict"], "EMPTY")
        vals = u["derived"]["values_in_base_survivors"]
        self.assertEqual(vals["statistics"], ["Z4 glue on the spinor module (P_glue)",       # z117, one bit
                                              "Z4 glue summed over a pair of seams (P_glue, P_paar)",                                       # z174
                                              "clock-twisted lift (P_GSO)"])                 # z97/z113, conditional
        self.assertEqual(vals["interaction"], ["chosen two-body term", "local coherent transfer gate",
                                               "twistor action from seam currents (T)"])    # z99, z120, z110/z111
        self.assertEqual(vals["gravity"], ["MacDowell-Mansouri: modular Cartan connection",  # z98, flat vacuum only
                                           "P0-graded two-seam composition",
                                           "double copy of seam amplitudes (P_diag)"])      # z112
        sc, info = scenario_candidates()
        excluded_interaction = ("fusion over even or linear algebras", "vertex operators = level-1 currents",
                                "transfer exchange (twirl)")
        for o in sc:
            if o["scenario"]["statistics"] != "untwisted lift" and o["scenario"]["spinor_index"] != "step type":
                self.assertNotIn("coherent", o["clauses"])
            if o["scenario"]["statistics"] == "untwisted lift" or o["scenario"]["interaction"] in excluded_interaction:
                self.assertNotIn("admissible", o["clauses"])
            if o["scenario"]["statistics"] in ("clock-twisted lift (P_GSO)", "Z4 glue on the spinor module (P_glue)"):
                self.assertIn("chiral_matter", o["clauses"])                               # z114 (T, K2, P_RJ), z117
                self.assertNotIn("modular_statistics", o["clauses"])                       # z113; z174: anomalous per seam
            if o["scenario"]["statistics"] == "Z4 glue summed over a pair of seams (P_glue, P_paar)":
                self.assertTrue({"chiral_matter", "modular_statistics"} <= o["clauses"])    # z117, z174
                self.assertTrue({"P_IF", "P_paar"} <= o["premises"])
        ev = info["evidence"]
        self.assertEqual(ev["statistics: E8 glue with A3 = conformal (P_A3)"]["status"], "excluded")
        self.assertEqual(ev["statistics: clock-twisted lift (P_GSO)"]["source"], "z113_ergebnis.json")
        self.assertEqual(ev["interaction: vertex operators = level-1 currents"]["status"], "excluded")
        self.assertEqual(ev["interaction: transfer exchange (twirl)"]["status"], "excluded")
        self.assertEqual(ev["interaction: local coherent transfer gate"]["status"], "conditional")
        self.assertEqual(ev["gravity: MacDowell-Mansouri: modular Cartan connection"]["status"], "conditional")
        self.assertEqual(ev["gravity: MacDowell-Mansouri: modular Cartan connection"]["source"], "z136_ergebnis.json")
        self.assertEqual(ev["interaction: local coherent transfer gate"]["source"], "z144_ergebnis.json")
        self.assertEqual(ev["statistics: Z4 glue on the spinor module (P_glue)"]["status"], "conditional")
        self.assertEqual((ev["statistics: Z4 glue summed over a pair of seams (P_glue, P_paar)"]["status"], ev["statistics: Z4 glue summed over a pair of seams (P_glue, P_paar)"]["source"]),
                         ("conditional", "z174_ergebnis.json"))
        self.assertEqual(ev["gravity: moved algebra at fixed state"]["status"], "excluded")
        self.assertEqual({k: v["status"] for k, v in info["clause_evidence"].items()},
                         {"modular_statistics": "conditional", "twistor_line_index": "conditional"})

    def test_premise_closure_names_what_the_best_candidate_still_assumes(self):
        u = self.universe("scenarios")
        pc = u["premise_closure"]
        # free parameters weigh first (K14): the chosen two-body term loses to the coherent gate, which trades P_sym for
        # P_univ (z125) and needs P_ord (z121) and P_uhr, which replaces P_null and P_viertel (z144); chirality costs one
        # bit P_IF on both routes (z117, z122). With the Tomita partner as second seam P_univ (z130) and P_uhr (z174)
        # merge into the one premise P_std; the pair-summed glue adds P_paar and is therefore not the minimal set here
        self.assertEqual((pc["survivors"], pc["least_free_parameters"], pc["least_premises"]), (162, 1, 8))
        self.assertEqual(pc["minimal_premise_sets"],
                         [("A*", "I_step", "P_G", "P_IF", "P_curv", "P_ord", "P_std", "P_zentral")])
        self.assertEqual(pc["clauses_no_survivor_has"], ["graviton"])
        self.assertEqual({i["replaces"]: i["holds"] for i in u["derived"]["implications"]},
                         {"P_RJ": True, "P_sym": True, "P_null": True, "P_viertel": True, "P_univ": True, "P_uhr": True,
                          "P_lift": True, "P_R": True, "P_pair": True})
        mm = [o for o in scenario_candidates()[0]
              if o["scenario"]["gravity"] == "MacDowell-Mansouri: modular Cartan connection"]
        self.assertTrue(all({"P_mod", "P_idx", "P_IF"} <= o["premises"] and not {"P_R", "P_pair"} & o["premises"] for o in mm))
        for p in pc["premise_frequency"]:
            self.assertIn(p, PREMISE_MEANING)                       # every premise is named, none dropped silently
        sc, _ = scenario_candidates()
        for o in sc:
            if o["scenario"]["statistics"] == "Z4 glue on the spinor module (P_glue)":
                self.assertIn("P_IF", o["premises"])

    def test_scenario_base_survivors_place_the_marks_only_where_coherence_allows(self):
        sc, _ = scenario_candidates()
        base = {"admissible", "coherent", "lorentzian_native", "positive_time"}
        for o in sc:
            if base <= o["clauses"] and o["scenario"]["marks"] == "celestial sphere":
                self.assertEqual(o["scenario"]["spinor_index"], "seam doublet V_1")
                self.assertEqual(o["scenario"]["time"], "OS reconstruction")
            if base <= o["clauses"] and o["scenario"]["marks"] == "internal spin-Z4":
                self.assertIn(o["scenario"]["spinor_index"], ("quarter oscillators", "step type"))

    # ------------------------------------------------------------ marks
    def test_mark_symmetries_are_d4_with_the_regular_klein_group_of_h6(self):
        info = symmetry_structure()
        self.assertEqual((info["group"], info["order"], info["abelian"]), ("D4", 8, False))
        self.assertEqual(info["element_orders"], [1, 2, 2, 2, 2, 2, 4, 4])
        self.assertTrue(info["regular_klein_is_group"])
        self.assertEqual(sorted(info["regular_klein_group"]), ["z -> -I/z", "z -> -z", "z -> I/z", "z -> z"])   # z36's V4
        self.assertIn("z -> I*z", [e["name"] for e in mark_symmetries()])

    def test_the_clock_fixed_poles_are_the_only_exceptional_orbit_of_size_two(self):
        orbits = exceptional_orbits()
        self.assertEqual([len(o) for o in orbits], [2, 4, 4])
        self.assertEqual(orbits[0], ["0", "oo"])
        self.assertIn(["-1", "-I", "1", "I"], orbits)

    # ------------------------------------------------------------ (B) normal bundles
    def test_double_cover_bundle_splits_as_o1_plus_o_minus1(self):
        self.assertEqual(pushforward(DOUBLE_COVER, 0), [0, -2])          # r_*O_E = O + O(-2)
        self.assertEqual(pushforward(DOUBLE_COVER), [1, -1])             # O(1) x r_*O_E = O(1) + O(-1)
        self.assertEqual(self.bundles["O(1) x r_*O_E (double cover)"]["degrees"], [1, -1])
        self.assertEqual((cover_genus(DOUBLE_COVER), cover_genus(CHARACTERS)), (1, 1))
        self.assertEqual(self.universe("normal_bundles")["cover_genus"]["h1(r_*O_E)"], 1)          # = genus
        self.assertEqual(self.universe("normal_bundles")["cover_genus"]["h1(r_*O_C)"], 1)

    def test_z2_cubed_cover_and_its_halves(self):
        self.assertEqual(pushforward(CHARACTERS), [1, 0, 0, 0, 0, 0, 0, -1])        # O(1) + 6 O + O(-1) (z34)
        plus = self.bundles["E+ under the point involution at mark 1"]
        minus = self.bundles["E- under the point involution at mark 1"]
        self.assertEqual((plus["degrees"], plus["h0"]), ([1, 0, 0, 0], 5))
        self.assertEqual((minus["degrees"], minus["h0"], minus["h1"]), ([0, 0, 0, -1], 3, 0))
        self.assertEqual(sorted(-a for a in plus["degrees"]), sorted(minus["degrees"]))   # E- = E+^dual
        self.assertEqual(plus["orbit_size"], 4)                    # the deck involution is one of four point involutions
        self.assertNotIn("choice_free", plus["clauses"])
        self.assertIn("mark_symmetric", plus["clauses"])           # no modification divisor
        central = self.bundles["E+ under the half period {1,-1}|{i,-i}"]
        self.assertEqual((central["degrees"], central["orbit_size"]), ([1, 0, 0, -1], 1))
        self.assertEqual(self.bundles["E+ under the half period {1,i}|{-1,-i}"]["orbit_size"], 2)

    def test_only_bundles_of_type_o1_plus_o1_are_twistor_lines(self):
        for name, b in self.bundles.items():
            self.assertEqual("twistor_line" in b["clauses"], b["degrees"] == [1, 1], name)
            if b["degrees"] == [1, 1]:
                self.assertEqual((b["h0"], b["h1"]), (4, 0))              # complete four-parameter family (Kodaira)
        hecke = {n: b for n, b in self.bundles.items() if n.startswith("Hecke")}
        self.assertEqual(len(hecke), 18)                                  # 15 mark subsets and 3 pole subsets
        for name, b in hecke.items():
            self.assertEqual("twistor_line" in b["clauses"], name.split("{")[1].count(",") == 1, name)

    def test_choice_free_twistor_lines(self):
        q = self.universe("normal_bundles")["questions"]
        self.assertEqual(q["twistor_line+choice_free"]["survivors"], [POLES, EULER, TOMITA])
        self.assertEqual(q["twistor_line+choice_free"]["verdict"], "AMBIGUOUS")
        self.assertEqual(q["twistor_line+mark_symmetric"]["survivors"], [POLES, EULER, TOMITA])
        self.assertEqual(len(q["twistor_line"]["survivors"]), 9)
        orbits = self.result["summary"]["normal_bundles"]["cover_derived_twistor_lines"]
        self.assertEqual(sorted(orbits.values()), [1, 2, 2, 4, 4, 4, 4])  # two marks: orbits 2 and 4; the poles: 1
        answers = self.result["answers"]
        self.assertEqual(answers["cover_derived_choice_free_twistor_lines"], [POLES])
        self.assertFalse(answers["choice_free_twistor_line_from_a_mark_divisor"])
        self.assertEqual(self.bundles["Hecke O(-1) -> O(-1)(D), D = {0}"]["orbit_size"], 2)     # one pole is a choice

    def test_lorentz_cone_needs_the_conjugate_doublet(self):
        self.assertEqual(self.result["summary"]["normal_bundles"]["lorentz_cone"], [TOMITA])
        self.assertEqual(self.universe("normal_bundles")["selected"], {TOMITA: [["lorentz_cone"]]})
        for name in (EULER, POLES):
            sig = {(r["type"], tuple(sorted(r["signature"]))) for r in self.bundles[name]["reality"]}
            self.assertEqual(sig, {("real", (2, 2)), ("quaternionic", (0, 4))}, name)   # Kleinian or Euclidean only
        tomita = sorted(tuple(r["signature"]) for r in self.bundles[TOMITA]["reality"] if r["type"] == "hermitian")
        self.assertEqual(tomita, [(1, 3), (3, 1)])
        self.assertEqual(self.result["summary"]["normal_bundles"]["seam_reality_signatures"], [[0, 4], [2, 2]])
        self.assertTrue(self.result["answers"]["lorentzian_real_form_only_with_conjugate_doublet"])

    # ------------------------------------------------------------ (A) doublet states
    def test_original_covariance_symmetries(self):
        sym = original_symmetries()
        self.assertTrue(sym["mirror x -> -x with particle-hole"])
        self.assertTrue(sym["quarter rotation z -> iz"])
        self.assertFalse(sym["reflection x -> -x"])
        self.assertFalse(sym["shift x -> x + 1"])
        self.assertFalse(sym["sector exchange x -> x + 1/2"])
        for x in (Fr(1, 4), Fr(-3, 4), Fr(7, 4)):
            self.assertEqual(sp.cancel(cov(x)/(1 - cov(x)) - U**int(-4*x)), 0)     # c/(1 - c) = e^{-pi k}, k = 4x

    def test_fock_state_of_each_sector_is_boltzmann_in_the_mark(self):
        for a in QUARTERS:
            counted, formula = fock_multiplicities(a)
            self.assertEqual(counted, formula)                                     # bosonization (Jacobi)
            self.assertEqual(counted, [1, 1, 0, 1, 1, 1, 1, 1, 2])                 # mark 2 is empty in each sector
            for s in fock_patterns(a):
                self.assertEqual(sp.cancel(_weight(s["flips"]) - U**-s["M"]), 0)   # weight e^{-pi M} relative to Omega

    def test_rho_m_is_faithful_kms_and_base_point_independent_for_every_m(self):
        derived = self.universe("doublet_states")["derived"]["rho_m"]
        for m in M_VALUES:
            corners = [s for s in self.states if s["construction"] == f"Fock corner, m = {m}"]
            self.assertTrue(corners)
            for s in corners:
                self.assertEqual(s["label"], RHO[m])                                 # diag(1, e^{-pi m})/(1 + e^{-pi m})
                self.assertTrue({"faithful", "kms", "base_point_independent", "symmetric"} <= s["clauses"])
                self.assertEqual(s["r_order"], 4//gcd(m, 4))
                self.assertEqual("neutral" in s["clauses"], m == 4)
                self.assertEqual("marks_separated" in s["clauses"], m % 2 == 1)
            self.assertEqual(derived[str(m)]["state"], RHO[m])
        self.assertEqual(self.universe("doublet_states")["derived"]["one_particle_pairs_exist_for_m"], [2, 4])
        for m in (2, 4):                                                             # N=1 corners of mode pairs agree
            pairs = [s for s in self.states if s["construction"] == f"N=1 corner of a mode pair, m = {m}"]
            self.assertEqual({s["label"] for s in pairs}, {RHO[m]})
        self.assertTrue(self.result["answers"]["rho_m_faithful_kms_base_point_independent_for_every_m"])
        rho1 = next(s for s in self.states if s["construction"] == "Fock corner, m = 1")
        self.assertEqual(rho1["base"], ["alpha = 1/4", "Omega (M=0, Q=0)", "particle 1/4 (M=1, Q=1)"])   # twist state

    def test_one_particle_restrictions_are_never_kms(self):
        rows = [s for s in self.states if s["construction"].startswith("one-particle")]
        self.assertEqual(len(rows), 16)
        self.assertTrue(all("faithful" in s["clauses"] and "kms" not in s["clauses"] for s in rows))
        self.assertEqual([s["base"] for s in rows if "symmetric" in s["clauses"]], [["-1/4", "1/4"]])
        fermi = next(s for s in rows if s["base"] == ["-1/4", "1/4"])
        self.assertEqual(fermi["label"], RHO[1])          # c(1/4) + c(-1/4) = 1: the matrix of rho_1, not KMS for gap 2 pi
        hardy = [s for s in self.states if "Hardy" in s["construction"]]
        self.assertTrue(all("faithful" not in s["clauses"] and "kms" not in s["clauses"] for s in hardy))

    def test_doublet_search_leaves_m_open_and_separated_marks_are_charged(self):
        u = self.universe("doublet_states")
        q = u["questions"]
        tie = q["faithful+kms+symmetric+base_point_independent"]
        self.assertEqual((tie["verdict"], tie["survivors"]), ("AMBIGUOUS", [RHO[2], RHO[3], RHO[4], RHO[1]]))
        self.assertEqual(q["kms+neutral"]["survivors"], [RHO[4]])
        self.assertEqual(q["kms+mirror_closed"]["survivors"], [RHO[2]])
        self.assertEqual(q["kms+marks_separated"]["survivors"], [RHO[3], RHO[1]])
        self.assertEqual(u["minimal_empty_rules"], [["marks_separated", "mirror_closed"], ["marks_separated", "neutral"],
                                                    ["mirror_closed", "neutral"]])
        self.assertIn(["kms", "neutral"], u["selected"][RHO[4]])
        self.assertIn(["kms", "mirror_closed"], u["selected"][RHO[2]])
        self.assertNotIn(RHO[1], u["selected"])
        answers = self.result["answers"]
        self.assertFalse(answers["doublet_state_fixed_by_faithful_kms_symmetric_base_point_independent"])
        self.assertEqual(answers["doublet_states_left_by_faithful_kms_symmetric_base_point_independent"], 4)
        self.assertFalse(answers["neutral_doublet_with_separated_marks_exists"])

    # ------------------------------------------------------------ spin-2 pencil (minimal)
    def test_spin2_pencil_minimal(self):
        a = sp.symbols("a0:5")
        psi = sum(c*Z**k for k, c in enumerate(a))
        M = spin1_matrix(psi)
        I_, J_, _ = invariants(psi)
        lam = sp.Symbol("lam")
        self.assertEqual(sp.expand(M.charpoly(lam).as_expr() - (lam**3 - I_*lam - 2*J_)), 0)   # tr 0, 2I, 2J
        self.assertEqual(sp.Matrix([[sp.expand(p).coeff(Z, k) for k in range(5)] for p in v4_pencil()]
                                   + [[1, 0, 0, 0, -1], [0, 0, 1, 0, 0]]).rank(), 2)         # span{1 - z^4, z^2}
        h6 = next(o for o in spin2_family_candidates()[0] if o["t"] == "0")
        self.assertEqual((h6["quartic"], h6["I"], h6["J"], h6["eigenvalues"]), ("1 - z**4", "-1", "0", ["-I", "0", "I"]))
        self.assertEqual(h6["kernel_phi00_phi01_phi11"], [["0", "1", "0"]])
        self.assertEqual(self.universe("spin2_family")["equivalent_clauses"], [["j_zero", "zero_eigenvalue"]])

    # ------------------------------------------------------------ run
    def test_run_is_fast_deterministic_and_promotes_nothing(self):
        self.assertLess(self.seconds, 20)
        self.assertFalse(self.result["physical_solution"])
        self.assertEqual(self.result["status"], "RECONSTRUCTION_SEARCHED")
        text = json.dumps(self.result, sort_keys=True, allow_nan=False)
        self.assertEqual(text, json.dumps(run_reconstruction(), sort_keys=True, allow_nan=False))
        self.assertEqual(self.result["objects_enumerated"], 54 + 66 + 10 + 6 + 4 + 57600)  # + doublets and position spaces (H6h40)


class CompleteSolutionDistance(unittest.TestCase):
    def test_distance_and_frontier(self):
        import reconstruction_search as rs
        d = rs.complete_solution_distance()
        self.assertFalse(d["complete"])
        # z173: (2, 1, 7) and eleven items; z174: the pair-summed glue carries modular_statistics under P_paar, and
        # P_univ, P_uhr merge into P_std
        self.assertEqual(d["distance"], {"missing_clauses": 1, "free_numbers": 1, "premises_beyond_axiom": 7})
        self.assertTrue(d["conventions"]["P_IF"]["holds"])
        # z178, z180: K19 follows from P_std, the GNS rule and A*, K20 from P_std and P_zentral with the geometric
        # partner; neither is an independent item any more, and neither is closed. Eight independent items remain
        self.assertEqual(d["open_total"], 8)
        self.assertEqual(d["by_question"], {"origin": [], "composition": ["P_ord", "P_std"],
                                            "gravity": ["P_G", "P_curv", "graviton"],
                                            "identification": ["I_step", "P_zentral"],
                                            "statistics": ["P_paar"]})
        self.assertEqual(sorted(d["bench_derived"]), ["K19", "K20"])
        self.assertEqual(d["bench_derived"]["K19"]["rests_on"], ["P_std"])
        self.assertEqual(d["bench_derived"]["K20"]["rests_on"], ["P_std", "P_zentral"])
        self.assertEqual(d["bench_open"], {})
        self.assertTrue(all("pair of seams" in label for label in d["best"]))

    def test_unification_candidates_close_nothing_by_derivation(self):
        import reconstruction_search as rs
        u = rs.unification_candidates()
        # z175 to z177: eleven candidates for one object behind several open items; z178 to z182: eight more;
        # z183: the attack on P_std adds one equivalence and one exclusion
        self.assertEqual(len(u["candidates"]), 21)
        self.assertEqual(u["counts"], {"excluded": 12, "holds": 6, "reading": 1, "derived": 2})
        self.assertEqual(u["closed_by_derivation"], [])
        self.assertEqual(u["derived_under_premises"], {"K19": ["P_std"], "K20": ["P_std", "P_zentral"]})
        self.assertEqual(u["open_total"], 8)
        verdict = {r["candidate"]: r["verdict"] for r in u["candidates"]}
        self.assertEqual(verdict["pair of seams = the two sheets of P1"], "excluded")
        self.assertEqual(verdict["number of seams from anomaly, positivity and normal statistics"], "holds")
        self.assertEqual(verdict["superselection of the pair clock = gauge sum"], "reading")
        self.assertEqual(verdict["graviton as an excitation of one glued object"], "excluded")
        self.assertEqual(verdict["gauge completion of the pair clock from superselection, the GNS rule and A*"], "derived")
        self.assertEqual(verdict["clock = deck parity of the sections extended to the seam (P_deck)"], "holds")
        self.assertEqual(verdict["graviton as a collective excitation of several seams"], "excluded")
        self.assertEqual(verdict["pair-clock part of P_std as locality: no one-sided current mixes the deck halves"],
                         "holds")
        self.assertEqual(verdict["P_std from A*, time neutrality or the covering of the four marks"], "excluded")
        # the z176 reading left eight items; z178 and z180 reach the same eight without the reading
        self.assertEqual(u["open_under_readings"],
                         ["I_step", "P_G", "P_curv", "P_ord", "P_paar", "P_std", "P_zentral", "graviton"])

    def test_derived_bench_row_needs_its_premise_and_its_result(self):
        import pathlib
        import tempfile
        import reconstruction_search as rs
        with tempfile.TemporaryDirectory() as tmp:                      # without the result files nothing is derived
            d = rs.complete_solution_distance(pathlib.Path(tmp))
            self.assertEqual(d["bench_derived"], {})
            self.assertEqual(sorted(d["bench_open"]), ["K19", "K20"])

    def test_counter_models(self):
        import reconstruction_search as rs
        c = rs.counter_models()
        rep = {r["conditions"]: r for r in c["completions_of_the_clock_sectors"]}
        self.assertEqual(rep["base changes S alone (fixed spin structure)"]["models"], 27)
        self.assertEqual(rep["A*: S and T"]["models"], 10)
        self.assertEqual(rep["A*: S and T"]["verdict"], "ambiguous")
        # A* and superselection alone do not gauge: several completions without charged local field remain
        self.assertEqual(rep["A* and no local field with pair-clock charge"]["verdict"], "ambiguous")
        self.assertEqual((rep["A*, superselection and the GNS rule"]["models"],
                          rep["A*, superselection and the GNS rule"]["verdict"]), (1, "unique"))
        # without T a spin form survives: flux 1 with charge 2 and spin 1/2
        without_t = rep["S, superselection and the GNS rule, without T"]
        self.assertEqual((without_t["models"], without_t["verdict"]), (2, "ambiguous"))
        # z183: conditions without a mirror bit, counted by spectrum. Time neutrality and orientation reversal of the
        # section bases leave counter-models to P_std; one local field with the deck parity as monodromy gives the gauge sum
        self.assertEqual({r["conditions"]: r["spectra"] for r in c["completions_of_the_clock_sectors"]}["A*: S and T"], 3)
        ex = rep["A* and time neutrality (exchange of seam and partner)"]
        self.assertEqual((ex["models"], ex["spectra"], ex["verdict"]), (8, 3, "ambiguous"))
        refl = rep["A* and orientation reversal of the section bases"]
        self.assertEqual((refl["models"], refl["spectra"], refl["verdict"]), (6, 2, "ambiguous"))
        mark = rep["A* and a local field with the deck parity as monodromy of the parents"]
        self.assertEqual((mark["models"], mark["spectra"], mark["verdict"]), (4, 1, "unique"))
        self.assertEqual(c["ambiguous"], ["value of G/a^2 under the area rule", "gate phase under P_uhr",
                                          "orientation of the circle through the four marks",
                                          "spectrum of the pair of seam and partner under A*"])
        self.assertEqual([r["verdict"] for r in c["recorded"]],
                         ["ambiguous", "ambiguous", "ambiguous", "unique", "ambiguous"])


if __name__ == "__main__":
    unittest.main()
