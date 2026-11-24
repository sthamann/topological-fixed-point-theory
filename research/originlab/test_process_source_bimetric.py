"""Meaningful exact tests for source scope, spectra and physical amplitudes."""
from fractions import Fraction
import json
import unittest

import sympy as sp

from process_source_bimetric import (
    Q_ORIGINAL_EXACT, bimetric_source_data, eikonal_family_data,
    massless_harmonic_eotvos, massive_yukawa_eotvos,
    quartic_forward_data, source_bimetric_audit, transfer_band_vs_bimetric,
)


class SourceBimetricTests(unittest.TestCase):
    def assertExact(self, left, right):
        self.assertEqual(sp.simplify(left-right), 0)

    def test_complete_analytic_audit_is_json(self):
        audit = source_bimetric_audit()
        self.assertTrue(audit["all_checks_pass"], audit["checks"])
        json.dumps(audit)
        self.assertFalse(audit["native_massless_EP_claim"]["Ward_identity_derived_by_Fourier_average"])
        self.assertFalse(audit["physical_solution"])

    def test_physical_source_residues_and_cross_channel(self):
        data = bimetric_source_data(2, 1, 3)
        self.assertEqual(data["massless_residue_matrix"], sp.ones(2)/5)
        self.assertEqual(data["massive_residue_matrix"], sp.Matrix([[sp.Rational(1, 4), -1], [-1, 4]])/5)
        self.assertEqual(data["massive_residue_matrix"].det(), 0)
        self.assertGreater(sp.trace(data["massive_residue_matrix"]), 0)
        self.assertEqual(data["mass_matrix"]*data["zero_mode"], sp.zeros(2, 1))

    def test_source_covariance_inverts_physical_quadratic_matrix(self):
        x = sp.Symbol("x", positive=True)
        data = bimetric_source_data(2, 1, 3)
        D = sp.diag(2, 1)
        parent = D*(x*sp.eye(2)+data["mass_matrix"])*D
        self.assertEqual(sp.simplify(parent*data["physical_TT_covariance"]), sp.eye(2))

    def test_massless_native_formula_and_nonzero_universal_counterexample(self):
        self.assertEqual(massless_harmonic_eotvos(Fraction(1, 2), 0, 1), -sp.Rational(2, 5))
        self.assertEqual(massless_harmonic_eotvos(Fraction(1, 2), 3, 3), 0)
        h = sp.Symbol("h", real=True)
        expr = massless_harmonic_eotvos(Q_ORIGINAL_EXACT, h, 0)
        self.assertExact(sp.diff(expr, h).subs(h, 0), Q_ORIGINAL_EXACT)

    def test_massive_common_metric_has_no_WEP_difference(self):
        z, alpha = sp.symbols("z alpha", positive=True)
        self.assertEqual(massive_yukawa_eotvos(-alpha, -alpha, -alpha, z), 0)

    def test_massive_range_is_separate_from_massless_source_variation(self):
        z = sp.Symbol("z", positive=True)
        eta = massive_yukawa_eotvos(1, 1, 0, z)
        self.assertEqual(sp.limit(eta, z, sp.oo), 0)
        self.assertEqual(massive_yukawa_eotvos(1, 1, 0, 0), sp.Rational(4, 5))
        self.assertNotEqual(massless_harmonic_eotvos(Fraction(1, 2), 1, 0), 0)

    def test_complete_band_cannot_match_one_massive_atom(self):
        comparison = transfer_band_vs_bimetric(Fraction(1, 2))
        self.assertEqual(comparison["massive_weight"], sp.Rational(1, 3))
        self.assertEqual(comparison["matching_pole_mass_squared"], sp.Rational(3, 4))
        self.assertExact(comparison["b_difference"], comparison["a"]**2)
        self.assertGreater(comparison["KL_hankel_MP"], 0)
        self.assertEqual(comparison["KL_hankel_one_massive_atom"], 0)
        self.assertFalse(comparison["full_response_match"])

    def test_original_nome_weight_is_exact(self):
        comparison = transfer_band_vs_bimetric(Q_ORIGINAL_EXACT)
        self.assertExact(comparison["massive_weight"], 1/(sp.exp(4*sp.pi)-1))
        self.assertFalse(comparison["b_difference"].has(sp.Float))

    def test_moment_scaling_is_preserved(self):
        base = transfer_band_vs_bimetric(Fraction(1, 2))
        scaled = transfer_band_vs_bimetric(Fraction(1, 2), mu_squared=3)
        for key, degree in (("a", 1), ("b_MP", 2), ("c_MP", 3), ("KL_hankel_MP", 2)):
            self.assertExact(scaled[key], base[key]/3**degree)

    def test_quartic_difference_and_forward_coefficient(self):
        data = quartic_forward_data(3, 4, 5)
        s = data["s"]
        self.assertExact(data["delta_amplitude_0000_forward"], s*(s-16)/60)
        self.assertEqual(data["half_second_s_derivative"], sp.Rational(1, 60))
        self.assertEqual(data["delta_amplitude_0000_forward"].subs(s, 16), 0)
        self.assertEqual(quartic_forward_data(0, 4, 5)["delta_amplitude_0000_forward"], 0)

    def test_eikonal_family_stays_flat_and_FP_normalized(self):
        tau = sp.Symbol("tau", real=True)
        data = eikonal_family_data(tau)
        b = data["beta"]
        self.assertExact(b[0]+3*b[1]+3*b[2]+b[3], 0)
        self.assertExact(b[4]+3*b[3]+3*b[2]+b[1], 0)
        self.assertExact(b[1]+2*b[2]+b[3], 1)
        self.assertEqual(data["c3"], sp.Rational(1, 4))
        self.assertExact(data["d5"], (tau-4)/96)

    def test_negative_input_controls(self):
        for q in (0, 1, -1, 0.5):
            with self.assertRaises((TypeError, ValueError)):
                transfer_band_vs_bimetric(q)
        with self.assertRaises(ValueError):
            bimetric_source_data(1, 0, 1)
        with self.assertRaises(TypeError):
            quartic_forward_data(0.1, 1, 1)
        with self.assertRaises(ValueError):
            massive_yukawa_eotvos(1, 1, 0, -1)


if __name__ == "__main__":
    unittest.main()
