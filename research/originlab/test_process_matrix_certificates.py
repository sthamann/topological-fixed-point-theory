import json
import unittest
from fractions import Fraction as Fr

import sympy

from process_matrix_certificates import matrix_certificate


class MatrixCertificates(unittest.TestCase):
    def test_hardy_positive_deformation_has_exact_boost_defect(self):
        result = matrix_certificate(1, 1)
        self.assertIn("hardy_slice", result["clauses"])
        self.assertIn("positive", result["clauses"])
        self.assertNotIn("wave_equation", result["clauses"])
        self.assertNotIn("conformally_covariant", result["clauses"])
        self.assertEqual(result["witnesses"]["conformally_covariant"]["defect"], "1/16")

    def test_repairing_degree_two_does_not_repair_positivity(self):
        result = matrix_certificate(Fr(1, 2), 100)
        self.assertEqual(result["witnesses"]["positive"], {
            "kind": "negative_schur_coefficient", "partition": [2, 1],
            "coefficient": "-1/8",
        })

    def test_allowed_negative_deformation_boundary(self):
        self.assertIn("positive", matrix_certificate(2, -1)["clauses"])
        witness = matrix_certificate(2, Fr(-101, 100))["witnesses"]["positive"]
        self.assertEqual(witness["partition"], [1, 1])
        self.assertEqual(witness["coefficient"], "-1/100")

    def test_single_zero_of_wave_residual_is_not_an_identity(self):
        result = matrix_certificate(2, -1)
        self.assertNotIn("wave_equation", result["clauses"])
        witness = result["witnesses"]["wave_equation"]
        self.assertEqual(witness["Z"], [["1/2", "0"], ["0", "1/2"]])
        self.assertIn("nonzero_reason", witness)

    def test_wave_identity_includes_nonpositive_exception(self):
        result = matrix_certificate(-1, -1)  # K=1-tr(W*Z)
        self.assertIn("wave_equation", result["clauses"])
        self.assertNotIn("positive", result["clauses"])
        self.assertEqual(result["witnesses"]["positive"]["partition"], [1, 0])

    def test_cayley_formula_and_boost_witness_are_exact(self):
        a, b, c, d = sympy.symbols("a b c d")
        q = (1-a)*(1-d)-b*c
        self.assertEqual(sympy.expand(sympy.diff(q,a)*sympy.diff(q,d)
                                     - sympy.diff(q,b)*sympy.diff(q,c)), q.expand())
        self.assertEqual(sympy.diff(q,a,d)-sympy.diff(q,b,c), 2)
        r = sympy.Rational(1, 2)
        for nu in (0, 1, 2, sympy.Rational(3, 2)):
            lhs = (1-r*r)**(-2*nu) + r**4
            rhs = (sympy.Rational(4, 3))**(2*nu)
            self.assertEqual(sympy.simplify(lhs-rhs), sympy.Rational(1, 16))

    def test_exact_rational_inputs_and_serialization(self):
        result = matrix_certificate(sympy.Rational(3, 2), sympy.Integer(0))
        self.assertEqual(result["nu"], "3/2")
        self.assertFalse(result["integer_determinant_multiplier"])
        self.assertIn("covering group", result["evidence"]["multiplier_scope"])
        json.dumps(result)
        with self.assertRaises(TypeError):
            matrix_certificate(1.0, 0)


if __name__ == "__main__":
    unittest.main()
