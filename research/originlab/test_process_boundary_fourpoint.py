import json
import unittest
from fractions import Fraction

import sympy as sp

from process_boundary_fourpoint import (
    boundary_fourpoint_audit, determinant_vertex, edge_source_response,
    edge_source_response_zero, flat_depth_contact_prefactor,
    flat_infinite_contact_data, fourpoint_prefactor, physical_two_point,
    star_boundary_data, tree_source_data,
)
from process_star_metric import star_metric_data
from process_weighted_star import weighted_star_data


class BoundaryFourpointTests(unittest.TestCase):
    def assertExact(self, actual, expected):
        self.assertEqual(sp.simplify(actual-expected), 0)

    def test_single_edge_has_physical_source_factors(self):
        Mg, Mf, k, x = sp.symbols("Mg Mf k x", positive=True)
        data = tree_source_data((Mg, Mf), ((0, 1),), (k,))
        mass = k*(Mg**-2+Mf**-2)
        self.assertExact(edge_source_response(data, x, 0)[0], -1/(Mg**2*(x+mass)))
        self.assertExact(edge_source_response(data, x, 1)[0], 1/(Mf**2*(x+mass)))
        self.assertExact(edge_source_response_zero(data, 0)[0], -Mf**2/(k*(Mg**2+Mf**2)))

    def test_single_edge_full_generic_momentum_vertex(self):
        data = tree_source_data((2, 3), ((0, 1),), (5,))
        xs = sp.symbols("x1 x2 x3 x4", positive=True)
        tau = sp.Symbol("tau", real=True)
        mass = sp.Rational(65, 36)
        result = fourpoint_prefactor(data, xs, taus=tau)
        self.assertExact(result, -5*tau/(2*2**8*sp.prod(x+mass for x in xs)))
        mixed = fourpoint_prefactor(data, xs, (0, 0, 0, 1), taus=tau)
        self.assertExact(mixed/result, -sp.Rational(4, 9))

    def test_full_finite_stars_match_existing_canonical_operators(self):
        for weighted in (False, True):
            weights = (1, 2) if weighted else (1, 1)
            data = star_boundary_data(2, Fraction(1, 2), 3, 5, weights)
            old = (weighted_star_data(2, sp.Rational(1, 2), 3, 5) if weighted
                   else star_metric_data(2, sp.Rational(1, 2), 3, 5))
            self.assertEqual(data["canonical_mass"], old["mass_matrix"])
            x = sp.Rational(7, 3)
            canonical = (x*sp.eye(9)+old["mass_matrix"]).inv()
            self.assertExact(physical_two_point(data, x), canonical[0, 0]/9)

    def test_single_layer_closed_contact_and_soft_limit(self):
        q, Mc, mu2 = sp.symbols("q Mc mu2", positive=True)
        x = sp.Symbol("x", positive=True)
        data = star_boundary_data(1, q, Mc, mu2)
        profile = edge_source_response(data, x)
        expected = -1/(Mc**2*(x+mu2*(1+q*q)))
        self.assertTrue(all(sp.simplify(p-expected) == 0 for p in profile))
        self.assertExact(fourpoint_prefactor(data, (0,)*4),
                         -q*q/(2*mu2**3*Mc**6*(1+q*q)**4))

    def test_general_subtree_flux_for_noncentral_sources(self):
        data = tree_source_data((2, 3, 5, 7), ((0, 1), (1, 2), (1, 3)), (11, 13, 17))
        x = sp.Symbol("x", positive=True)
        source = (1, -2, 3, -4)
        finite = edge_source_response(data, x, source)
        soft = edge_source_response_zero(data, source)
        self.assertTrue(all(sp.limit(a, x, 0) == b for a, b in zip(finite, soft)))
        # A central source gives the same sign on every outward edge.
        self.assertTrue(all(value < 0 for value in edge_source_response_zero(data, 0)))

    def test_massless_pole_cancels_and_any_pure_massless_leg_is_invisible(self):
        data = star_boundary_data(2, sp.Rational(1, 2), 3, 5)
        x = sp.Symbol("x", positive=True)
        self.assertExact(sp.limit(x*physical_two_point(data, x), x, 0),
                         1/data["M_eff_squared"])
        null_source = tuple(data["mass_squared_weights"])
        self.assertEqual(edge_source_response(data, 3, null_source), sp.zeros(8, 1))
        for leg in range(4):
            sources = [0]*4
            sources[leg] = null_source
            self.assertEqual(fourpoint_prefactor(data, (1, 2, 3, 4), sources), 0)

    def test_polarized_determinant_normalization_from_polynomial(self):
        A, B = sp.diag(0, 1, 1, -2), sp.diag(1, 0, 1, -2)
        z = sp.symbols("z1:5")
        polynomial = sp.Poly(((z[0]+z[2])*A+(z[1]+z[3])*B).det(), *z)
        self.assertEqual(determinant_vertex(A, B, A, B),
                         polynomial.coeff_monomial(sp.prod(z)))
        self.assertEqual(determinant_vertex(A, B, A, B), -16)
        T = sp.diag(1, 1, -1, -1)
        self.assertEqual(determinant_vertex(T, T, T, T), 24*T.det())

    def test_physical_TT_visibility_and_collinear_exception(self):
        momenta = [sp.Matrix(p) for p in ((1, 0, 0, 0), (0, 1, 0, 0),
                                         (-1, 0, 0, 0), (0, -1, 0, 0))]
        A, B = sp.diag(0, 1, 1, -2), sp.diag(1, 0, 1, -2)
        tensors = (A, B, A, B)
        self.assertEqual(sum(momenta, sp.zeros(4, 1)), sp.zeros(4, 1))
        for p, T in zip(momenta, tensors):
            self.assertEqual(T*p, sp.zeros(4, 1))
            self.assertEqual(sp.trace(T), 0)
        data = star_boundary_data(1, sp.Rational(1, 2))
        self.assertNotEqual(fourpoint_prefactor(data, (1,)*4)*determinant_vertex(*tensors), 0)
        self.assertEqual(determinant_vertex(A, A, A, A), 0)

    def test_boundary_amputation_is_eighth_order_in_common_soft_momentum(self):
        data = star_boundary_data(1, sp.Rational(1, 2), 3, 5)
        x = sp.Symbol("x", positive=True)
        response = fourpoint_prefactor(data, (x,)*4)
        C = physical_two_point(data, x)
        constant = fourpoint_prefactor(data, (0,)*4)
        self.assertNotEqual(constant, 0)
        self.assertExact(sp.limit(response/C**4/x**4, x, 0),
                         constant*data["M_eff_squared"]**4)
        self.assertEqual(sp.limit(response/C**4, x, 0), 0)

    def test_independent_edge_couplings_and_arm_averaging(self):
        data = star_boundary_data(2, sp.Rational(1, 2))
        a = fourpoint_prefactor(data, (1, 2, 3, 4), taus=(1, 2, 3, 4, 5, 6, 7, 8))
        b = fourpoint_prefactor(data, (1, 2, 3, 4),
                               taus=(sp.Rational(5, 2),)*4+(sp.Rational(13, 2),)*4)
        self.assertExact(a, b)
        self.assertNotEqual(fourpoint_prefactor(data, (0,)*4, taus=(1,)*4+(0,)*4),
                            fourpoint_prefactor(data, (0,)*4, taus=(0,)*4+(1,)*4))

    def test_infinite_flat_profile_and_closed_geometric_sum(self):
        q, mu2 = sp.Rational(1, 2), sp.Integer(3)
        result = flat_infinite_contact_data(q, (2, 2, 2, 2), 5, mu2, tau=7)
        r, z, x = q*q, result["z"][0], sp.Rational(2, 3)
        self.assertExact(r*z*z-(x+1+r)*z+1, 0)
        raw = (sp.sqrt((x+(1-q)**2)*(x+(1+q)**2))-x+1-r)/(2*x)
        self.assertExact(raw*(z-1), -z)
        P = result["P"]
        finite = flat_depth_contact_prefactor(q, (2,)*4, (7,)*5, 5, mu2)
        tail = -7*(r*P)**6/(2*mu2**3*5**6*(1-r*P))
        self.assertExact(result["euclidean_prefactor"]-finite, tail)

    def test_IR_sum_does_not_identify_depth_couplings(self):
        q, r = sp.Rational(1, 2), sp.Rational(1, 4)
        self.assertEqual(flat_depth_contact_prefactor(q, (0,)*4, (r, -1)), 0)
        self.assertNotEqual(flat_depth_contact_prefactor(q, (1,)*4, (r, -1)), 0)

    def test_input_domains(self):
        with self.assertRaises(ValueError):
            tree_source_data((1, 2, 3), ((1, 2), (2, 1)), (1, 1))
        with self.assertRaises(ValueError):
            flat_infinite_contact_data(1, (1,)*4)
        with self.assertRaises(TypeError):
            star_boundary_data(1, 0.5)
        data = star_boundary_data(1, Fraction(1, 2))
        with self.assertRaises(ValueError):
            fourpoint_prefactor(data, (1,)*4, taus=(1, 2))

    def test_audit_is_exact_json_and_has_no_dynamics_promotion(self):
        report = boundary_fourpoint_audit()
        self.assertTrue(report["all_checks_pass"], report["checks"])
        self.assertFalse(report["physical_solution"])
        self.assertFalse(report["coupling_selection"]["tau_predicted"])
        self.assertIn("Gaussian", report["coupling_selection"]["Gaussian_state_is_not_a_zero_gravity_cumulant"])
        json.dumps(report, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
