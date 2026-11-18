"""Independent controls for exact positive-kernel recovery."""

from __future__ import annotations

import copy
import unittest

import sympy as sp

from exact_certificate import recover
from word_algebra import T, WordAlgebra


def _positive_payload(*, noise: float = 0.0) -> dict:
    """One-mode exact model: quarter-covariant vacuum and nontrivial T."""
    algebra = WordAlgebra(1)
    gamma0 = sp.Matrix([[0, 1], [1, 0]])
    gamma1 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    transfer = sp.diag(1, sp.Rational(1, 2))
    matrices = {0: gamma0, 1: gamma1, T: transfer}
    omega = sp.Matrix([1, 0])

    def moment(word):
        operator = sp.eye(2)
        for token in word:
            operator *= matrices[token]
        return sp.expand((omega.adjoint() * operator * omega)[0])

    moments = []
    # degree=2 certificate also verifies the next Gram level d+1=3, hence
    # moments through degree 6 are supplied.
    for index, word in enumerate(algebra.words(6)):
        value = moment(word)
        perturb = noise if index == 7 else 0.0
        moments.append({
            "word": list(word),
            "re": float(sp.re(value)) + perturb,
            "im": float(sp.im(value)),
        })
    return {
        "m": 1,
        "degree": 2,
        "transfer_floor": 0.0,
        "selection_scope": "test fixture only",
        "moments": moments,
    }


class ExactCertificateTests(unittest.TestCase):
    def test_recovers_nontrivial_positive_flat_model(self):
        result = recover(_positive_payload(noise=2e-8))
        self.assertTrue(result["ok"], result.get("error"))
        self.assertEqual(result["verdict"], "EXACT_FINITE_POSITIVE_GNS_KERNEL_WITH_EFFECT_T")
        self.assertEqual(result["basis"]["rank"], 2)
        names = {row["name"] for row in result["checks"] if row["passed"]}
        self.assertIn("T_fixes_Omega", names)
        self.assertIn("quarter_order_four", names)
        self.assertIn("all_supplied_moments_replay", names)
        transfer = result["exact_matrices"]["generators"]["T"]
        self.assertNotEqual(transfer, result["exact_matrices"]["metric_G"])

    def test_rejects_missing_next_degree_moment(self):
        payload = _positive_payload()
        payload["moments"].pop()
        result = recover(payload)
        self.assertFalse(result["ok"])
        self.assertIn("missing recovered moment", result["error"])

    def test_rejects_exact_negative_T_localizer(self):
        payload = _positive_payload()
        # omega(gamma0 T gamma0)=1/2 is the excited-state T expectation.
        target = next(row for row in payload["moments"] if row["word"] == [0, T, 0])
        target["re"] = -0.5
        result = recover(payload)
        self.assertFalse(result["ok"])
        self.assertTrue(
            "positive" in result["error"] or "represented by the recovered flat basis" in result["error"],
            result["error"],
        )

    def test_rejects_floor_above_excited_T_eigenvalue(self):
        payload = _positive_payload()
        payload["transfer_floor"] = 0.75
        result = recover(payload)
        self.assertFalse(result["ok"])
        self.assertIn("T-transfer_floor", result["error"])

    def test_rejects_float_without_allowed_rational_recovery(self):
        payload = copy.deepcopy(_positive_payload())
        payload["moments"][3]["re"] = 0.123456789
        result = recover(payload, max_denominator=8, tolerance=1e-9)
        self.assertFalse(result["ok"])
        self.assertIn("no rational", result["error"])

    def test_denominator_64_control_rejects_one_over_256(self):
        payload = copy.deepcopy(_positive_payload())
        payload["moments"][0]["re"] = 1 / 256
        result = recover(payload, max_denominator=64, tolerance=1e-6)
        self.assertFalse(result["ok"])
        self.assertIn("denominator <= 64", result["error"])


if __name__ == "__main__":
    unittest.main()
