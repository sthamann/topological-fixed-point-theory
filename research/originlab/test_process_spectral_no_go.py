import unittest

import sympy as sp

from process_spectral_no_go import spectral_no_go_audit


class SpectralNoGo(unittest.TestCase):
    def test_profile_fubini_factors_cancel_with_integration_by_parts(self):
        result = spectral_no_go_audit()["analytic_identities"]
        self.assertEqual(sp.sympify(result["first_beta_factor"]), 3*sp.pi)
        self.assertEqual(sp.sympify(result["second_beta_factor"]), 3*sp.pi/2)
        self.assertEqual(sp.sympify(result["sum_rule_factor"]), 0)

    def test_positive_schur_control_fails_the_zero_sum_rule(self):
        x, M = sp.symbols("x M", positive=True)
        C = 1/x+1/(x+M)
        H = sp.factor(1/C)
        self.assertEqual(sp.limit(H/x, x, sp.oo), sp.Rational(1, 2))
        h = sp.factor(H-x/2)
        self.assertEqual(sp.simplify(h-M*x/(2*(2*x+M))), 0)
        integral = sp.integrate(h/x**sp.Rational(3, 2), (x, 0, sp.oo))
        self.assertEqual(sp.simplify(integral-sp.pi*sp.sqrt(M)/(2*sp.sqrt(2))), 0)

    def test_using_ir_slope_as_uv_slope_is_an_invalid_control(self):
        x, M = sp.symbols("x M", positive=True)
        H = 1/(1/x+1/(x+M))
        self.assertEqual(sp.limit((H-x)/x, x, sp.oo), -sp.Rational(1, 2))

    def test_scope_does_not_confuse_axis_positivity_with_spectral_positivity(self):
        result = spectral_no_go_audit()
        self.assertTrue(result["einstein_addition_cannot_restore_positive_spectral_measure"])
        self.assertFalse(result["every_einstein_augmented_kernel_has_a_euclidean_zero"])
        self.assertIn("full interacting quantum positivity", result["not_proved"])


if __name__ == "__main__":
    unittest.main()
