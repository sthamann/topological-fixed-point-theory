"""Negative controls protect the scope of the flavour obstructions."""
import json
import unittest
import sympy as sp

from process_flavour_obstructions import (
    algebraic_integer_certificate, delta96_majorana_certificate, direct_reactor_fractions, direct_reactor_scan,
    flavour_obstruction_audit, full_klein_certificate, minimal_doublet_alignment_data, rational_weight_angle_certificate,
    semidirect_potential_data, tm1_mu_tau_matrix, Z49_REACTOR_RANGE,
)


class FlavourObstructionTests(unittest.TestCase):
    def test_finite_order_trace_has_correct_integrality_control(self):
        self.assertTrue(algebraic_integer_certificate(2*sp.cos(sp.pi/5))["is_algebraic_integer"])
        self.assertFalse(algebraic_integer_certificate(-(2+sp.sqrt(3))/4)["is_algebraic_integer"])
        with self.assertRaises(ValueError):
            algebraic_integer_certificate(0.5)

    def test_full_klein_proof_does_not_gate_semidirect_model(self):
        d = full_klein_certificate()
        self.assertTrue(all(d["checks"].values()))
        self.assertFalse(d["semidirect_Z2_CP_excluded"])
        self.assertEqual(d["quadratic_field_norm_of_nonintegral_traces"], "1/16")

    def test_named_angle_is_an_input_and_other_angles_are_allowed(self):
        for angle in (0, sp.pi/13, sp.pi/12):
            U = tm1_mu_tau_matrix(angle)
            self.assertEqual(sp.simplify(U.conjugate().T*U), sp.eye(3))
            self.assertEqual(sp.simplify(U[0, 2]*sp.conjugate(U[0, 2])-sp.sin(angle)**2/3), 0)

    def test_Delta96_match_changes_the_neutral_symmetry(self):
        self.assertTrue(all(delta96_majorana_certificate()["checks"].values()))

    def test_rational_thermal_weights_include_degenerate_escape(self):
        t = sp.Symbol("t")
        witness = rational_weight_angle_certificate((1+t)/(1-t), t/(1-t), 0, t)
        self.assertTrue(witness["pi_over_12_excluded_at_transcendental_t"])
        degenerate = rational_weight_angle_certificate(1/(1-t), 0, 1/(1-t), t)
        self.assertTrue(degenerate["degenerate_angle_unselected"])
        self.assertFalse(degenerate["pi_over_12_excluded_at_transcendental_t"])

    def test_new_algebraic_or_nonrational_weights_are_not_falsely_excluded(self):
        t = sp.Symbol("t")
        for A, B, C in [(1, 1/sp.sqrt(24), 0), (sp.sin(t), 1, 0), (0.2, 1, 0)]:
            with self.assertRaises(ValueError):
                rational_weight_angle_certificate(A, B, C, t)

    def test_three_corner_near_miss_is_exactly_not_fifteen_degrees(self):
        t = sp.Symbol("t")
        d = rational_weight_angle_certificate(sp.Rational(7, 3), -sp.Rational(1, 3), sp.Rational(2, 3), t)
        self.assertEqual(sp.sympify(d["angle_equation"]), sp.Rational(1, 9))

    def test_larger_groups_cross_reactor_windows_but_not_solar_gate(self):
        self.assertEqual(direct_reactor_fractions(17), ())
        old = direct_reactor_scan(22, Z49_REACTOR_RANGE)
        new = direct_reactor_scan(22)
        self.assertEqual(old["hits"][0]["n"], 18)
        self.assertEqual(new["hits"][0]["n"], 22)
        self.assertEqual(new["hits"][0]["angle_over_pi"], "2/33")
        self.assertGreater(sp.Rational(1, 3), sp.Rational("0.3295"))

    def test_stable_semidirect_potential_keeps_direction_free(self):
        d = semidirect_potential_data()
        self.assertTrue(all(d["checks"].values()))
        self.assertFalse(d["angle_dynamically_derived_from_TFPT"])
        self.assertEqual(d["tan_2theta"], "y/x=b/a")

    def test_angular_invariant_does_not_silently_become_PMNS_prediction(self):
        d = minimal_doublet_alignment_data()
        self.assertTrue(all(d["checks"].values()))
        self.assertFalse(d["PMNS_angle_derived"])
        self.assertIn("V4", d["minimal_S4_limit"])

    def test_json_audit_preserves_scope(self):
        d = flavour_obstruction_audit()
        json.dumps(d)
        self.assertTrue(d["all_exact_checks_pass"])
        self.assertFalse(d["physical_solution"])
        self.assertFalse(d["thermal_rational_weights"]["excludes_all_S4_models"])
        self.assertFalse(d["z49_scope"]["all_integers_through_16_scanned"])
        self.assertTrue(d["larger_direct_family"]["intervals_are_not_joint_likelihood"])


if __name__ == "__main__":
    unittest.main()
