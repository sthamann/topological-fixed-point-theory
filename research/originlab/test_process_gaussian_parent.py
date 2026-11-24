"""Exact current, KMS and Gaussian-parent checks with no sampled proof gate."""

from fractions import Fraction
import json
import unittest

import sympy as sp

from process_gaussian_parent import (
    EXACT_Q_ORIGINAL, current_kms_certificate, gaussian_eigenpolynomial,
    gaussian_parent_audit, gaussian_parent_certificate, gaussian_resolvent_partial,
)


class GaussianParentTests(unittest.TestCase):
    def test_current_commutator_shift_is_not_a_state_space_shift(self):
        audit = gaussian_parent_audit()
        witness = audit["current_operator"]["finite_CAR_witness"]
        self.assertTrue(witness["CAR_commutator_shifts_each_fermion_creation_mode"])
        self.assertEqual(witness["current_powers_norm_squared"], [1, 1, 2])
        self.assertEqual(witness["unilateral_shift_powers_norm_squared"], [1, 1, 1])
        self.assertEqual(audit["current_operator"]["full_current_level_dimensions_0_through_5"],
                         [1, 1, 2, 3, 5, 7])

    def test_full_KMS_fixes_original_q_without_a_separate_one_mode_Gaussian_assumption(self):
        cert = current_kms_certificate()
        self.assertEqual(sp.sympify(cert["canonical_tfd_q"]), EXACT_Q_ORIGINAL)
        self.assertEqual(sp.simplify(sp.sympify(cert["occupation"])-1/(sp.exp(4*sp.pi)-1)), 0)
        self.assertTrue(cert["gaussianity_of_this_mode_is_a_consequence"])
        # Independent bounded matrix-unit KMS relation.
        r = sp.Rational(1, 4)
        rho = sp.diag(*[(1-r)*r**n for n in range(5)])
        for n in range(4):
            A = sp.zeros(5); A[n, n+1] = 1
            B = A.T
            self.assertEqual(sp.trace(rho*A*(r*B)), sp.trace(rho*B*A))

    def test_equal_occupation_and_two_point_relation_do_not_imply_KMS(self):
        cert = current_kms_certificate(sp.log(4))
        example = cert["same_mean_counterexample"]
        p = [sp.sympify(value) for value in example["first_three_probabilities"]]
        r, Z = sp.Rational(1, 4), sp.Rational(3, 4)
        delta = [p[n]-Z*r**n for n in range(3)]
        self.assertTrue(all(value > 0 for value in p))
        self.assertEqual(sum(delta), 0)
        self.assertEqual(sum(n*delta[n] for n in range(3)), 0)
        self.assertEqual(sum(n*(n-1)*delta[n] for n in range(3)),
                         sp.sympify(example["fourth_normal_ordered_moment_excess"]))
        self.assertNotEqual(p[1]-r*p[0], 0)
        self.assertEqual(p[1]-r*p[0], sp.sympify(example["full_KMS_defect_p1_minus_r_p0"]))

    def test_single_mode_geometric_coefficients_fail_a_Gaussian_annihilator(self):
        q, u = sp.symbols("q u", nonzero=True)
        # n=0 fixes alpha=u*q; n=1 fixes v=u*q^2*(1-sqrt(2)).
        alpha = u*q
        v = u*q*q*(1-sp.sqrt(2))
        residual = sp.simplify(u*q*sp.sqrt(3)+v/q*sp.sqrt(2)-alpha)
        self.assertEqual(residual, u*q*(-3+sp.sqrt(2)+sp.sqrt(3)))
        self.assertNotEqual(sp.sqrt(2)+sp.sqrt(3)-3, 0)
        self.assertFalse(gaussian_parent_audit()["one_mode_geometric_state"]["pure_Gaussian_for_nonzero_q"])

    @staticmethod
    def compressed_jacobi(q, scale, count, weighted=True):
        A = sp.diag(*([q]*count))
        for n in range(count-1):
            A[n, n+1] = -1
        W = sp.diag(*range(1, count+1)) if weighted else sp.eye(count)
        return scale*A.T*W*A

    def test_weighted_factorization_and_first_distinguishing_mass_moment(self):
        q, scale, count = sp.Rational(1, 2), sp.Integer(3), 5
        r = q*q
        weighted = self.compressed_jacobi(q, scale, count)
        flat = self.compressed_jacobi(q, scale, count, weighted=False)
        for n in range(count):
            self.assertEqual(weighted[n, n], scale*((1+r)*n+r))
        for n in range(count-1):
            self.assertEqual(weighted[n, n+1], -scale*q*(n+1))
        cert = gaussian_parent_certificate(q, scale)
        moments = [sp.sympify(value) for value in cert["raw_mass_moments"]]
        for degree in range(4):
            self.assertEqual((weighted**degree)[0, 0], moments[degree])
        for degree in range(3):
            self.assertEqual((weighted**degree)[0, 0], (flat**degree)[0, 0])
        self.assertEqual((weighted**3-flat**3)[0, 0], scale**3*r*r)

    def test_geometric_spectral_weights_reproduce_independent_mass_moments(self):
        r = sp.Symbol("r", positive=True)
        generating = 1/(1-r)
        moment = generating
        cert = gaussian_parent_certificate(Fraction(1, 2), scale=3)
        delta = sp.sympify(cert["pair_spacing"])
        for degree in range(4):
            expected = ((1-r)*delta**degree*moment).subs(r, sp.Rational(1, 4))
            self.assertEqual(sp.simplify(expected-sp.sympify(cert["raw_mass_moments"][degree])), 0)
            moment = sp.diff(moment, r)*r
        self.assertEqual(sp.sympify(cert["full_two_mode_gap"])*2, delta)

    def test_meixner_recurrence_certifies_eigenvalue_spacing(self):
        q = sp.Rational(1, 2)
        m = sp.Symbol("m", integer=True, nonnegative=True)
        for n in range(5):
            f = gaussian_eigenpolynomial(n, m, q)
            up = gaussian_eigenpolynomial(n+1, m, q)
            down = gaussian_eigenpolynomial(n-1, m, q) if n else 0
            result = ((1+q*q)*n+q*q)*f-q*(n+1)*up-q*n*down
            self.assertEqual(sp.expand(result-(1-q*q)*m*f), 0)

    def test_meixner_spectral_transform_has_the_stated_norm(self):
        q = sp.Rational(1, 2)
        m = sp.Symbol("m", integer=True, nonnegative=True)
        r = sp.Symbol("r", positive=True)
        geometric_moments = [1/(1-r)]
        for _ in range(6):
            geometric_moments.append(r*sp.diff(geometric_moments[-1], r))
        for n in range(3):
            for j in range(3):
                polynomial = sp.Poly(gaussian_eigenpolynomial(n, m, q)*gaussian_eigenpolynomial(j, m, q), m)
                inner = sum(coefficient*geometric_moments[powers[0]]
                            for powers, coefficient in polynomial.terms())
                value = sp.simplify(((1-r)*inner).subs(r, q*q))
                self.assertEqual(value, int(n == j))

    def test_resolvent_remainder_has_exact_positive_bounds(self):
        short = gaussian_resolvent_partial(Fraction(1, 2), x=2, terms=2, scale=3)
        long = gaussian_resolvent_partial(Fraction(1, 2), x=2, terms=8, scale=3)
        self.assertGreater(long["partial"], short["partial"])
        self.assertLess(long["partial"]-short["partial"], short["tail_bound"])
        self.assertLess(long["tail_bound"], short["tail_bound"])
        self.assertLessEqual(long["partial"]+long["tail_bound"], sp.Rational(1, 2))
        self.assertGreaterEqual(long["partial"], sp.Rational(3, 8))

    def test_inverse_coefficients_follow_the_positive_pole_tower(self):
        cert = gaussian_parent_certificate(Fraction(1, 2), scale=3)
        x = sp.Symbol("x")
        r, delta = sp.Rational(1, 4), sp.sympify(cert["pair_spacing"])
        H = x
        for degree, name, sign in ((2, "a", -1), (3, "b", 1), (4, "c", -1)):
            H += sign*sp.sympify(cert["normalized_inverse_coefficients"][name])*x**degree
        C = 1/x-sp.log(1-r)/delta-sp.polylog(2, r)*x/delta**2+sp.polylog(3, r)*x*x/delta**3
        product = sp.expand(H*C)
        self.assertEqual(product.coeff(x, 0), 1)
        for degree in range(1, 4):
            self.assertEqual(sp.simplify(product.coeff(x, degree)), 0)

    def test_named_symmetries_select_equal_parent_frequencies(self):
        a, b, c, d = sp.symbols("a b c d")
        G = sp.Matrix([[a, b], [c, d]])
        charge = sp.diag(1, -1)
        exchange = sp.Matrix([[0, 1], [1, 0]])
        constraints = list(charge*G-G*charge)+list(exchange*G-G*exchange)
        solution = sp.solve(constraints, [b, c, d], dict=True)
        self.assertEqual(solution, [{b: 0, c: 0, d: a}])

    def test_zero_mode_limit_and_invalid_parameters(self):
        cert = gaussian_parent_certificate(0)
        self.assertEqual(cert["raw_massless_residue"], "1")
        self.assertEqual(cert["normalized_inverse_coefficients"], {"a": "0", "b": "0", "c": "0"})
        self.assertFalse(cert["boundary_is_cyclic"])
        self.assertEqual(gaussian_resolvent_partial(0, x=2)["partial"], sp.Rational(1, 2))
        for q in (0.5, True, sp.Float("0.5"), -1, 1, 2, sp.I/2, sp.Symbol("q")):
            with self.subTest(q=q), self.assertRaises((TypeError, ValueError)):
                gaussian_parent_certificate(q)
        with self.assertRaises(ValueError):
            gaussian_parent_certificate(Fraction(1, 2), scale=0)
        with self.assertRaises(ValueError):
            current_kms_certificate(0)

    def test_full_report_is_JSON_and_physical_identifications_remain_explicit(self):
        report = gaussian_parent_audit()
        json.dumps(report)
        self.assertTrue(report["all_checks_pass"])
        self.assertFalse(report["physical_solution"])
        self.assertFalse(report["original_q_gaussian_parent"]["flat_hardy_comparison"]["same_process_operator"])
        self.assertTrue(report["original_q_gaussian_parent"]["preparation_parent_unique_under_named_symmetries"]["remaining_scale_free"])


if __name__ == "__main__":
    unittest.main()
