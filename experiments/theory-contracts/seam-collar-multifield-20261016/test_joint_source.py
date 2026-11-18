"""Independent identities for source integration, scalar and phase provenance."""
import unittest
import sympy as s
from joint_source import (berezin_source, determinant_source, schur_source,
                          profile_obstruction, sea_energy_zeta, circle_source,
                          spectral_state_ward, control)


class JointSourceTests(unittest.TestCase):
    def test_direct_grassmann_integral_matches_independent_matrix_formula(self):
        # This nondiagonal block is an algebraic sign control, not a new source.
        for k in (s.Matrix([[3]]), s.Matrix([[2, s.I], [-s.I, 3]])):
            self.assertEqual(berezin_source(k), determinant_source(k))

    def test_schur_keeps_absolute_interior_determinant(self):
        k = s.Matrix([[2, 1], [1, 3]])
        eff, inner = schur_source(k, [0])
        self.assertEqual(eff, s.Matrix([[s.Rational(5, 3)]]))
        self.assertEqual(inner*eff.det(), k.det())
        self.assertNotEqual(eff.det(), k.det())
        self.assertEqual(eff.inv(), k.inv().extract([0], [0]))
        with self.assertRaises(ValueError):
            schur_source(s.Matrix([[1, 1], [1, 0]]), [0])

    def test_profile_identity_is_not_a_spectral_determinant_relabeling(self):
        x = s.Symbol('x', positive=True)
        self.assertEqual(profile_obstruction(-s.log(x)+7, x), 0)
        self.assertEqual(profile_obstruction(s.exp(-x*x), x).subs(x, 1), 1-2/s.E)

    def test_abel_counterterm_and_conjugate_domain(self):
        eps = s.Symbol('eps', positive=True)
        alpha = s.Rational(1, 4)
        regulated_count = s.exp(-(1-alpha)*eps)/(1-s.exp(-eps))
        expansion = s.series(s.diff(regulated_count, eps), eps, 0, 1).removeO()
        self.assertEqual(s.expand(expansion).coeff(eps, 0), sea_energy_zeta(alpha))
        self.assertEqual(sea_energy_zeta(alpha), -s.Rational(1, 96))
        self.assertEqual(sea_energy_zeta(alpha), sea_energy_zeta(1-alpha))

    def test_actual_collar_joint_source_and_scope(self):
        report = control()
        self.assertEqual(report['joint_outputs']['determinant_scalar'], '69649/65536')
        self.assertEqual(report['joint_outputs']['normal_ordered_occupation_diag'], ['1/17','1/4097'])
        self.assertFalse(report['absolute_original_TFPT_source_recovered'])
        self.assertFalse(report['conditional_absolute_zeta_control']['P1_P2_selection_of_prescription_proved'])

    def test_actual_angular_dirac_green_phase_and_normalizer(self):
        for alpha in (s.Rational(1,4),s.Rational(3,4),s.Rational(1,3)):
            a = circle_source(alpha)
            self.assertEqual(s.simplify(a['determinant']-(1-a['holonomy'])), 0)
            self.assertEqual(s.simplify(a['green_after']-a['green_before']), s.I)
            self.assertEqual(s.simplify(a['green_after']-a['holonomy']*a['green_before']), 0)
        with self.assertRaises(ValueError):
            circle_source(0)

    def test_circle_ward_contact_term_is_retained(self):
        alpha = s.Symbol('alpha', real=True)
        h = s.exp(2*s.pi*s.I*alpha)
        deriv = s.diff(s.log(1-h), alpha).subs(alpha, s.Rational(1,4))
        actual = circle_source(s.Rational(1,4))
        self.assertEqual(s.simplify(deriv-actual['zeta_logdet_derivative']), 0)
        self.assertEqual(actual['zeta_logdet_derivative']-actual['pv_inverse_trace'], s.I*s.pi)

    def test_same_spectral_source_ward_response_matches_independent_fock_trace(self):
        beta = 8*s.log(4)
        for alpha in (s.Rational(1,4),s.Rational(3,4),s.Rational(1,2)):
            ward = spectral_state_ward(alpha, beta)
            # Independent trace over |0>, |particle>, |hole>, |particle,hole>.
            energy = sea_energy_zeta(alpha)
            weights = [s.exp(-beta*(energy+exc)) for exc in (0,alpha,1-alpha,1)]
            # The common shift is obtained from the separate eta invariant.
            qsea = -circle_source(alpha)['eta']/2
            trace_charge = s.simplify(sum(w*(qsea+q) for w,q in zip(weights,(0,1,-1,0)))/sum(weights))
            self.assertEqual(ward['charge_from_log_partition'], trace_charge)
            for key in ('energy_eta_residual','sea_derivative_eta_residual','thermal_ward_residual'):
                self.assertEqual(ward[key], 0)
        actual = spectral_state_ward(s.Rational(1,4), beta)
        self.assertEqual(actual['charge_from_log_partition'], -s.Rational(53329,278596))
        with self.assertRaises(ValueError):
            spectral_state_ward(0, beta)
        with self.assertRaises(ValueError):
            spectral_state_ward(s.Rational(1,4), 0)


if __name__ == '__main__':
    unittest.main()
