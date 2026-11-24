"""Independent checks of the weighted metric parent and its limiting response."""
import json
import unittest

import numpy as np
import sympy as sp
from scipy.linalg import solve_banded

from process_weighted_star import radial_finite_matrix, weighted_star_audit, weighted_star_data


class WeightedStarTests(unittest.TestCase):
    def test_symbolic_audit_is_exact_and_scoped(self):
        result = weighted_star_audit()
        self.assertTrue(result["all_checks_pass"], result["checks"])
        self.assertFalse(result["physical_solution"])
        self.assertFalse(result["quantum_completion_proved"])
        self.assertFalse(result["G_is_predicted"])
        json.dumps(result)

    def test_full_operator_and_all_sectors_at_other_lengths(self):
        for length in (1, 2, 4):
            q = sp.Rational(2, 3)
            data = weighted_star_data(length, q, M_c=3, scale=sp.Rational(7, 2))
            self.assertEqual(data["fourier_mass"], data["native_weighted_mass"])
            radial = radial_finite_matrix(length, q, sp.Rational(7, 2))
            self.assertEqual(data["blocks"][0], radial)
            for block in data["blocks"][1:]:
                self.assertEqual(block, radial[1:, 1:])
            self.assertEqual(len(data["mass_matrix"].nullspace()), 1)

    def test_boundary_resolvent_matches_independent_discrete_measure(self):
        # Rectangular graph truncation and its exact infinite geometric measure.
        # q is deliberately large enough that more than one massive atom matters.
        q, scale, length = 0.6, 1.7, 180
        n = np.arange(length+1)
        diagonal = scale*(n+q*q*(n+1))
        diagonal[-1] = scale*length
        off = -scale*q*np.arange(1, length+1)
        rhs = np.zeros(length+1)
        rhs[0] = 1
        for x in (0.04, 0.8, 3.0):
            banded = np.zeros((3, length+1))
            banded[1] = diagonal+x
            banded[0, 1:] = off
            banded[2, :-1] = off
            actual = solve_banded((1, 1), banded, rhs)[0]
            levels = np.arange(500)
            expected = np.sum((1-q*q)*q**(2*levels)/(x+scale*(1-q*q)*levels))
            self.assertAlmostEqual(actual, expected, places=11)

    def test_other_sectors_are_positive_and_interlace(self):
        q, length = 0.55, 90
        n = np.arange(length+1)
        diagonal = n+q*q*(n+1)
        diagonal[-1] = length
        off = -q*np.arange(1, length+1)
        radial = np.diag(diagonal)+np.diag(off, 1)+np.diag(off, -1)
        eigenvalues = np.linalg.eigvalsh(radial[1:, 1:])
        gap = 1-q*q
        for k in range(8):
            self.assertGreater(eigenvalues[k], gap*k)
            self.assertLess(eigenvalues[k], gap*(k+1))
        self.assertGreater(eigenvalues[0], 0)

    def test_state_does_not_fix_unweighted_vs_weighted_dynamics(self):
        q, length = sp.Rational(1, 2), 3
        data = weighted_star_data(length, q)
        A = data["native_A_rect"]
        unweighted = A.T*A
        weighted = data["native_weighted_mass"]
        self.assertNotEqual(weighted, unweighted)
        vector = data["native_zero"]
        self.assertEqual(weighted*vector, sp.zeros(13, 1))
        self.assertEqual(unweighted*vector, sp.zeros(13, 1))


if __name__ == "__main__":
    unittest.main()
