import json
import unittest
from fractions import Fraction as Fr

import sympy as sp

from process_spin2_certificates import (
    CLAUSE_NAMES, EXACT_COUNTEREXAMPLES, STANDARD_KINDS,
    candidate_certificate, rational_certificate,
)


def _json_ready(value):
    if isinstance(value, Fr):
        return str(value)
    if isinstance(value, set):
        return sorted(value)
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


def _sym(value):
    return sp.Rational(value.numerator, value.denominator)


class Spin2Certificates(unittest.TestCase):
    def test_selected_single_channel_has_exact_equal_residues(self):
        result = rational_certificate(["1/2"], ["5/3"])
        self.assertEqual(result["clauses"], set(CLAUSE_NAMES))
        self.assertEqual((result["a"], result["b"], result["c"]),
                         (Fr(3, 10), Fr(9, 50), Fr(27, 250)))
        self.assertEqual(result["total_massive_weight"], 1)
        x = sp.symbols("x", positive=True)
        H = x-sp.Rational(1,2)*x*x/(x+sp.Rational(5,3))
        self.assertEqual(sp.factor(1/H-1/x-1/(x+sp.Rational(10,3))), 0)

    def test_two_counterexamples_survive_all_scalar_clauses(self):
        next_coefficients = []
        for case in EXACT_COUNTEREXAMPLES:
            result = rational_certificate(case["weights"], case["masses"])
            self.assertEqual(result["clauses"], set(CLAUSE_NAMES))
            self.assertEqual(result["a"], Fr(3, 10))
            self.assertEqual(result["total_massive_weight"], 1)
            next_coefficients.append(result["b"])
        self.assertEqual(next_coefficients, [Fr(6, 25), Fr(13, 50)])
        self.assertNotIn(Fr(9, 50), next_coefficients)

    def test_parent_matrix_independently_reproduces_counterexample_response(self):
        x = sp.symbols("x", positive=True)
        for case in EXACT_COUNTEREXAMPLES:
            ws = [_sym(w) for w in case["weights"]]
            ms = [_sym(m) for m in case["masses"]]
            B = sp.eye(3)
            for j, w in enumerate(ws, 1):
                B[0,j] = B[j,0] = sp.sqrt(w)
            A = sp.diag(0, *ms)
            parent = A+x*B
            C = sp.factor(parent.inv()[0,0])
            H = x-sum(w*x*x/(x+m) for w,m in zip(ws,ms))
            self.assertEqual(sp.factor(C*H-1), 0)
            self.assertEqual(B.det(), sp.Rational(1,2))
            self.assertEqual(sp.limit(x*C,x,0), 1)
            self.assertEqual(sp.limit(x*C,x,sp.oo), 2)
            series = sp.series(H,x,0,5).removeO().expand()
            result = rational_certificate(case["weights"],case["masses"])
            self.assertEqual(series.coeff(x,2), -_sym(result["a"]))
            self.assertEqual(series.coeff(x,3), _sym(result["b"]))
            self.assertEqual(series.coeff(x,4), -_sym(result["c"]))

    def test_general_weight_is_not_silently_called_mirror(self):
        result = rational_certificate(["1/4"], ["5/6"])
        self.assertIn("dirac_low_energy", result["clauses"])
        self.assertIn("kallen_lehmann", result["clauses"])
        self.assertEqual(result["total_massive_weight"], Fr(1,3))
        self.assertNotIn("mirror_weight", result["clauses"])

    def test_inputs_outside_positive_pencil_class_are_rejected_not_classified(self):
        for weights,masses in (([0],[1]), ([-1],[1]), (["1/2"],[0]),
                               ([1],[1]), (["3/4","1/4"],[1,2]), (["1/2"],[])):
            with self.subTest(weights=weights,masses=masses), self.assertRaises(ValueError):
                rational_certificate(weights,masses)
        with self.assertRaises(TypeError):
            rational_certificate([0.5],["5/3"])
        with self.assertRaises(TypeError):
            rational_certificate(["1/2"],[True])

    def test_empty_auxiliary_set_is_einstein(self):
        empty = rational_certificate([],[])
        named = candidate_certificate("einstein")
        self.assertEqual(empty["clauses"],named["clauses"])
        self.assertEqual(empty["total_massive_weight"],0)
        self.assertEqual((empty["a"],empty["b"],empty["c"]),(0,0,0))

    def test_standard_positive_and_negative_controls(self):
        positive_euclidean = {"einstein","stelle","ledger_entire","exp_plus_einstein"}
        for kind in STANDARD_KINDS:
            result = candidate_certificate(kind)
            self.assertEqual(result["checks"]["no_euclidean_zero"],kind in positive_euclidean)
            self.assertEqual(result["checks"]["kallen_lehmann"],kind == "einstein")
        self.assertEqual(candidate_certificate("stelle")["witnesses"]["kallen_lehmann"]["residue"],-1)
        self.assertEqual(candidate_certificate("dirac_local")["witnesses"]["no_euclidean_zero"]["x"],Fr(10,3))

    def test_exponential_taylor_data_from_beta_integrals(self):
        x = sp.symbols("x")
        F = sum((-x)**n*sp.factorial(n)/sp.factorial(2*n+1) for n in range(6))
        H = sp.expand((3*x+12)*F-12)
        result = candidate_certificate("exp")
        self.assertEqual(H.coeff(x,1),1)
        self.assertEqual(H.coeff(x,2),-_sym(result["a"]))
        self.assertEqual(H.coeff(x,3),_sym(result["b"]))
        self.assertEqual(H.coeff(x,4),-_sym(result["c"]))
        self.assertEqual(result["witnesses"]["kallen_lehmann"]["b_minus_a_squared"],-Fr(19,350))

    def test_step_zero_is_within_exact_polynomial_region(self):
        x,t = sp.symbols("x t")
        v=x*t*(1-t)
        H=sp.integrate(3*x*(1-v)+6*(1-v)**2-6,(t,0,1))
        self.assertEqual(sp.expand(H),x-sp.Rational(3,10)*x*x)
        result=candidate_certificate("step")
        zero=result["witnesses"]["no_euclidean_zero"]["x"]
        self.assertTrue(0<zero<4)
        self.assertEqual(H.subs(x,_sym(zero)),0)

    def test_laplace3_normalization_and_taylor_coefficients(self):
        x,t=sp.symbols("x t")
        v=x*t*(1-t)
        integrand=sp.series(3*x/(1+v)**2+12/(1+v)-12,x,0,5).removeO()
        H=sp.integrate(integrand,(t,0,1)).expand()
        result=candidate_certificate("laplace3")
        self.assertEqual(H.coeff(x,1),1)
        self.assertEqual(H.coeff(x,2),-_sym(result["a"]))
        self.assertEqual(H.coeff(x,3),_sym(result["b"]))
        self.assertEqual(H.coeff(x,4),-_sym(result["c"]))
        self.assertEqual(result["witnesses"]["no_euclidean_zero"]["H_at_infinity"],-6)

    def test_einstein_addition_is_positive_but_fails_next_hankel_test(self):
        result=candidate_certificate("exp_plus_einstein",kappa="2")
        self.assertIn("no_euclidean_zero",result["clauses"])
        self.assertNotIn("kallen_lehmann",result["clauses"])
        self.assertEqual(result["infrared_slope"],3)
        self.assertEqual(result["massless_residue"],Fr(1,3))
        a,b,c=result["a"],result["b"],result["c"]
        self.assertGreater(b,a*a)
        self.assertEqual(a*c-b*b,-Fr(13,264600))
        self.assertEqual(result["witnesses"]["kallen_lehmann"]["a_c_minus_b_squared"],a*c-b*b)
        for kappa in ("5/2",10):
            other=candidate_certificate("exp_plus_einstein",kappa=kappa)
            self.assertTrue(other["checks"]["no_euclidean_zero"])
            self.assertFalse(other["checks"]["kallen_lehmann"])
        with self.assertRaises(ValueError):
            candidate_certificate("exp_plus_einstein",kappa="1")

    def test_exact_results_can_be_serialized_by_the_lane(self):
        cases=[candidate_certificate(kind) for kind in STANDARD_KINDS]
        cases += [rational_certificate(case["weights"],case["masses"]) for case in EXACT_COUNTEREXAMPLES]
        for result in cases:
            json.dumps(_json_ready(result))
        with self.assertRaises(ValueError):
            candidate_certificate("unknown")
        with self.assertRaises(ValueError):
            candidate_certificate("einstein",kappa=2)


if __name__ == "__main__":
    unittest.main()
