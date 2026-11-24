"""Exact metric, transfer, symmetry, source, and boundary-condition checks."""
import json
import unittest

import sympy as sp

from process_star_metric import star_metric_audit, star_metric_data


class StarMetricTests(unittest.TestCase):
    def test_symbolic_full_audit_and_json(self):
        audit = star_metric_audit()
        self.assertTrue(audit["all_checks_pass"], audit["checks"])
        json.dumps(audit)
        self.assertFalse(audit["physical_solution"])
        self.assertFalse(audit["quantum_or_infinite_interacting_completion_proved"])

    def test_native_full_operator_and_first_order_for_multiple_lengths(self):
        q = sp.Rational(1, 2)
        for m in (1, 2, 4):
            d = star_metric_data(m, q, M_c=3, mu_squared=sp.Rational(7, 2))
            U, V, A, B = (d[k] for k in ("U_vertex_fourier", "V_edge_fourier", "native_A_rect", "canonical_incidence"))
            self.assertEqual(B*U-V*A, sp.zeros(4*m, 4*m+1))
            self.assertEqual(sp.simplify(U.conjugate().T*d["mass_matrix"]*U-sp.Rational(7, 2)*A.T*A), sp.zeros(4*m+1))

    def test_all_four_blocks_and_physical_endpoint(self):
        m, q = 3, sp.Rational(2, 3)
        d = star_metric_data(m, q)
        native_mass = d["native_A_rect"].T*d["native_A_rect"]
        for r in range(4):
            indices = list(range(r, 4*m+1, 4))
            block = native_mass.extract(indices, indices)
            expected = sp.diag(*([q*q]+[1+q*q]*(m-1)+[1] if r == 0 else [1+q*q]*(m-1)+[1]))
            for n in range(len(indices)-1):
                expected[n, n+1] = expected[n+1, n] = -q
            self.assertEqual(block, expected)
        # A square quadratic compression would keep 1+q^2 at these terminal
        # entries.  The actual graph has 1 and an exact single zero mode.
        self.assertEqual(native_mass[-1, -1], 1)
        self.assertEqual(len(native_mass.nullspace()), 1)

    def test_universal_source_normalization(self):
        d = star_metric_data(2, sp.Rational(1, 2), M_c=4)
        self.assertEqual(d["M_eff_squared"], 21)
        self.assertEqual(d["physical_pole_matrix"], sp.ones(9)/21)
        self.assertEqual(d["boundary_zero_weight"], sp.Rational(16, 21))
        self.assertEqual(d["source_couplings"], sp.ones(9, 1)/sp.sqrt(21))

    def test_quarter_rotation_without_dropping_three_sectors(self):
        d = star_metric_data(2, sp.Rational(3, 4))
        U, R, Rn = (d[k] for k in ("U_vertex_fourier", "rotation", "native_rotation"))
        self.assertEqual(R**4, sp.eye(9))
        self.assertEqual(R*U, U*Rn)
        self.assertEqual(R*d["mass_matrix"], d["mass_matrix"]*R)
        self.assertEqual(R.eigenvals(), {sp.Integer(1): 3, sp.Integer(-1): 2, sp.I: 2, -sp.I: 2})

    def test_universal_zero_heavy_heavy_overlap(self):
        d = star_metric_data(2, sp.Rational(1, 2))
        coeff = d["source_couplings"][0]
        self.assertEqual(d["cubic_zero_weight_matrix"], coeff*sp.eye(9))
        # This identity makes every orthonormal change of massive basis
        # yield delta_ab/M_eff without choosing numeric eigenvectors.

    def test_parameter_negative_controls(self):
        for length in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                star_metric_data(length, sp.Rational(1, 2))
        for q in (0, -1, 0.5, sp.Symbol("q")):
            with self.assertRaises((ValueError, TypeError)):
                star_metric_data(2, q)


if __name__ == "__main__":
    unittest.main()
