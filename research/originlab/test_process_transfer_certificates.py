"""Exact identities and negative controls for Hardy transfer certificates."""
from fractions import Fraction
import json
import unittest

import sympy as sp

from process_transfer_certificates import (
    EXACT_Q_DIRAC, EXACT_R_DIRAC, EXACT_TRANSFER_Q_VALUES,
    dirac_matching_scale, transfer_certificate,
)


class TransferCertificateTests(unittest.TestCase):
    def assertExact(self, left, right):
        self.assertEqual(sp.simplify(left-right), 0)

    def test_half_coefficients(self):
        cert = transfer_certificate(Fraction(1, 2))
        self.assertEqual(cert["a"], sp.Rational(4, 9))
        self.assertEqual(cert["b"], sp.Rational(80, 81))
        self.assertEqual(cert["c"], sp.Rational(1856, 729))
        self.assertEqual(cert["raw_massless_residue"], sp.Rational(3, 4))
        self.assertEqual(cert["total_massive_weight"], sp.Rational(1, 3))
        self.assertEqual(cert["spectral_band"], (sp.Rational(1, 4), sp.Rational(9, 4)))
        self.assertEqual(cert["invariant_b_over_a2"], 5)
        self.assertEqual(cert["clauses"], {"no_euclidean_zero", "kallen_lehmann", "transfer_memory"})

    def test_exact_algebraic_dirac_match(self):
        cert = transfer_certificate(EXACT_Q_DIRAC)
        self.assertExact(cert["r"], EXACT_R_DIRAC)
        self.assertEqual(cert["a"], sp.Rational(3, 10))
        self.assertExact(cert["b"], 3*(11+sp.sqrt(55))/100)
        self.assertExact(cert["c"], 3*(8+sp.sqrt(55))/40)
        self.assertExact(cert["total_massive_weight"], (sp.sqrt(55)-5)/10)
        self.assertTrue(cert["checks"]["dirac_low_energy"])
        self.assertFalse(cert["checks"]["mirror_weight"])

    def test_coefficients_against_implicit_resolvent_equation(self):
        # The exact resolvent solves x*Craw^2+(x-Z)*Craw-1=0.
        # Substitute the certified inverse series into its independent
        # polynomial equation H^2-Z*(x-Z)*H-Z^2*x=0.
        x = sp.Symbol("x")
        for q in (sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4), EXACT_Q_DIRAC):
            cert = transfer_certificate(q)
            H = x-cert["a"]*x**2+cert["b"]*x**3-cert["c"]*x**4
            Z = cert["Z"]
            equation = sp.Poly(sp.expand(H**2-Z*(x-Z)*H-Z**2*x), x)
            for degree in range(5):
                self.assertExact(equation.nth(degree), 0)

    def test_scale_law_and_invariant(self):
        base = transfer_certificate("1/2")
        scaled = transfer_certificate("1/2", scale=Fraction(7, 3))
        self.assertExact(scaled["a"], base["a"]*sp.Rational(3, 7))
        self.assertExact(scaled["b"], base["b"]*sp.Rational(3, 7)**2)
        self.assertExact(scaled["c"], base["c"]*sp.Rational(3, 7)**3)
        self.assertEqual(scaled["total_massive_weight"], base["total_massive_weight"])
        self.assertEqual(scaled["invariant_b_over_a2"], base["invariant_b_over_a2"])

    def test_every_supported_nonzero_q_can_match_dirac(self):
        for _, q in EXACT_TRANSFER_Q_VALUES:
            with self.subTest(q=q):
                cert = transfer_certificate(q, scale=dirac_matching_scale(q))
                self.assertExact(cert["a"], sp.Rational(3, 10))
                self.assertTrue(cert["checks"]["dirac_low_energy"])

    def test_mirror_and_dirac_simultaneously_after_scaling(self):
        cert = transfer_certificate(sp.sqrt(sp.Rational(1, 2)), scale=Fraction(20, 3))
        self.assertEqual(cert["a"], sp.Rational(3, 10))
        self.assertEqual(cert["b"], sp.Rational(27, 100))
        self.assertEqual(cert["c"], sp.Rational(297, 1000))
        self.assertEqual(cert["total_massive_weight"], 1)
        self.assertTrue(cert["checks"]["dirac_low_energy"])
        self.assertTrue(cert["checks"]["mirror_weight"])
        self.assertFalse(cert["checks"]["rational"])

    def test_nome_is_exact_and_not_rounded_to_a_match(self):
        for q in ("e^-pi", "e^-2pi", sp.exp(-sp.pi/3)):
            cert = transfer_certificate(q)
            self.assertEqual(cert["parameter_class"], "nome")
            self.assertFalse(cert["checks"]["dirac_low_energy"])
            self.assertFalse(cert["checks"]["mirror_weight"])
            self.assertTrue(cert["checks"]["kallen_lehmann"])
            self.assertFalse(cert["a"].has(sp.Float))

    def test_massless_limit_is_explicit(self):
        cert = transfer_certificate(0, scale="7/2")
        self.assertEqual((cert["a"], cert["b"], cert["c"]), (0, 0, 0))
        self.assertEqual(cert["total_massive_weight"], 0)
        self.assertIsNone(cert["spectral_band"])
        self.assertIsNone(cert["invariant_b_over_a2"])
        self.assertTrue(cert["checks"]["rational"])
        self.assertFalse(cert["nontrivial_memory"])
        with self.assertRaises(ValueError):
            dirac_matching_scale(0)

    def test_rejects_float_and_outside_class_instead_of_sampling(self):
        for q in (0.5, sp.Float("0.5"), True):
            with self.assertRaises((TypeError, ValueError)):
                transfer_certificate(q)
        for q in (-1, "-1/2", 1, 2, sp.I/2, sp.Symbol("q")):
            with self.assertRaises(ValueError):
                transfer_certificate(q)
        for scale in (0, -1, 1.0, sp.Symbol("lambda")):
            with self.assertRaises((TypeError, ValueError)):
                transfer_certificate("1/2", scale=scale)

    def test_positive_hankel_witnesses(self):
        for q in (sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4)):
            cert = transfer_certificate(q)
            self.assertGreater(cert["witnesses"]["moment_hankel"]["b_minus_a_squared"], 0)
            self.assertGreater(cert["witnesses"]["moment_hankel"]["a_c_minus_b_squared"], 0)

    def test_json_conversion_contract(self):
        def exact_default(value):
            if isinstance(value, sp.Basic):
                return str(value)
            if isinstance(value, set):
                return sorted(value)
            raise TypeError(type(value).__name__)
        encoded = json.dumps(transfer_certificate("q_dirac"), default=exact_default)
        loaded = json.loads(encoded)
        self.assertEqual(loaded["a"], "3/10")
        self.assertIn("kallen_lehmann", loaded["clauses"])

    def test_connected_source_scope_separates_particle_and_vacuum(self):
        scope = transfer_certificate("1/2")["source_scope"]
        self.assertIn("retains its massless", scope["one_particle_mass_squared"])
        self.assertIn("no vacuum pole", scope["vacuum_hamiltonian"])


if __name__ == "__main__":
    unittest.main()
