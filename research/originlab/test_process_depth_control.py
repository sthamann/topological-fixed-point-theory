"""Checks of actual normalizations, interval decisions, and scope."""
from fractions import Fraction
import json
import unittest

import sympy as sp

from process_depth_control import (
    depth_certificate, depth_control_audit, edge_depth_data, native_q_bounds,
    truncated_preparation_data,
)


class DepthControlTests(unittest.TestCase):
    def test_native_exponential_is_bounded_without_float_predicates(self):
        result = native_q_bounds()
        self.assertTrue(result["lower_proved"])
        self.assertTrue(result["upper_proved"])
        self.assertEqual(result["lower"], Fraction(1867442, 10**9))

    def test_first_edge_and_bulk_normalizations(self):
        q, M, mu = sp.symbols("q M mu", positive=True)
        first = edge_depth_data(1, q, M, mu)
        second = edge_depth_data(2, q, M, mu, weighted=True)
        self.assertEqual(first["rho"], q/2)
        self.assertEqual(second["rho"], q)
        self.assertEqual(first["isolated_edge_mass_squared"], mu**2*(1+q*q/4))
        self.assertEqual(second["isolated_edge_mass_squared"], 2*mu**2*(1+q*q))
        for data in (first, second):
            self.assertEqual(sp.simplify(data["edge_stiffness"]/data["M_relative"]**2
                                         - data["isolated_edge_mass_squared"]), 0)

    def test_five_layers_and_zero_based_index_are_distinct(self):
        q = native_q_bounds()
        for weighted in (False, True):
            result = depth_certificate(q["lower"], q["upper"], Fraction(1, 10**16), weighted=weighted)
            self.assertEqual(result["status"], "CERTIFIED_LOCAL_DEPTH")
            self.assertEqual(result["certified_layers"], 5)
            self.assertEqual(result["last_zero_based_source_index"], 4)
            self.assertEqual(result["first_failed_depth"], 6)

    def test_scale_separation_margin_changes_depth(self):
        q = native_q_bounds()
        ratio = Fraction(1000, 2435*10**15)
        for weighted in (False, True):
            counts = [depth_certificate(q["lower"], q["upper"], ratio, weighted=weighted,
                                        margin=margin)["certified_layers"]
                      for margin in (1, 10, 100)]
            self.assertEqual(counts, [5, 4, 3])

    def test_relative_mass_correction_matters_away_from_tiny_q(self):
        result = depth_certificate(Fraction(1, 2), Fraction(1, 2), Fraction(1, 10**10))
        self.assertEqual(result["certified_layers"], 31)
        # The source's approximate node check M_child>=mu instead gives32.
        self.assertGreater(Fraction(1, 2)**33, Fraction(1, 10**10))

    def test_quartic_weight_does_not_double_the_q_exponent(self):
        q, M, mu = sp.symbols("q M mu", positive=True)
        flat2 = edge_depth_data(2, q, M, mu)
        flat3 = edge_depth_data(3, q, M, mu)
        weighted3 = edge_depth_data(3, q, M, mu, weighted=True)
        self.assertEqual(sp.simplify(flat3["quartic_relative_field_coefficient_per_tau"]
                                    / flat2["quartic_relative_field_coefficient_per_tau"]), q**-2)
        self.assertEqual(sp.simplify(weighted3["quartic_transfer_row_coefficient_per_tau"]
                                    / flat3["quartic_transfer_row_coefficient_per_tau"]), 3)

    def test_uncertain_interval_is_not_relabelled_as_proof(self):
        result = depth_certificate(Fraction(1, 10), Fraction(9, 10), Fraction(1, 10))
        self.assertEqual(result["status"], "UNDECIDED_INTERVAL")
        self.assertIsNone(result["first_failed_depth"])
        with self.assertRaises(TypeError):
            depth_certificate(0.1, Fraction(1, 2), Fraction(1, 1000))

    def test_audit_is_serializable_and_not_a_global_cutoff_claim(self):
        result = depth_control_audit()
        self.assertTrue(result["all_checks_pass"], result["checks"])
        self.assertFalse(result["physical_solution"])
        self.assertIn("collective multigravity cutoff", result["not_derived"])
        json.dumps(result)

    def test_finite_D2_full_parent_and_thermal_marginal(self):
        q = sp.Rational(2, 3)
        result = truncated_preparation_data(2, q, scale=sp.Rational(7, 3))
        self.assertEqual(result["paired_parent"], result["physical_radial_parent"])
        self.assertEqual(result["paired_parent"][-1, -1], sp.Rational(14, 3))
        self.assertEqual(result["full_parent"].rank(), 8)
        self.assertEqual(result["full_parent"]*result["unnormalized_tfd"], sp.zeros(9, 1))
        self.assertEqual(sum(result["pair_ground_weights"]), 1)
        for n in (0, 1):
            self.assertEqual(result["pair_ground_weights"][n+1]/result["pair_ground_weights"][n], q*q)
        self.assertEqual(result["single_mode_commutator"], sp.diag(1, 1, -2))

    def test_symbolic_finite_pair_coefficients_and_boundary_defect(self):
        q, scale = sp.symbols("q scale", positive=True)
        result = truncated_preparation_data(2, q, scale)
        expected = scale*sp.Matrix([[q*q, -q, 0], [-q, 1+2*q*q, -2*q], [0, -2*q, 2]])
        self.assertEqual(sp.simplify(result["paired_parent"]-expected), sp.zeros(3))
        n, D = sp.symbols("n D", integer=True, nonnegative=True)
        # Each annihilator contributes n+q^2*(n+1) below the top.
        interior_from_two_modes = scale*(2*n+2*q*q*(n+1))/2
        self.assertEqual(sp.expand(interior_from_two_modes-scale*(n+q*q*(n+1))), 0)
        # At the top, both raising norms vanish by aa^dagger|D>=0.
        terminal_from_two_modes = scale*(2*D+q*q*(0+0))/2
        self.assertEqual(terminal_from_two_modes, scale*D)


if __name__ == "__main__":
    unittest.main()
