"""Regression tests for the actual charged alpha=1/4 collar pair."""
from __future__ import annotations

import sympy as sp
import unittest

from charged_pair_control import (
    Q_BASE,
    RAW_RESPONSE,
    RESPONSE,
    charged_fixed_moments,
    run_moment_search,
    verify_charged_pair,
)


class ChargedPairControlTests(unittest.TestCase):
    def test_actual_charged_pair_reconstructs_full_operator_ray(self):
        report = verify_charged_pair()
        self.assertTrue(report["ok"])
        self.assertEqual(RESPONSE, sp.diag(sp.Rational(1, 4), sp.Rational(1, 64)))
        self.assertEqual(RAW_RESPONSE, sp.diag(sp.Rational(1, 4), 64))
        self.assertEqual(
            report["exact_results"]["raw_signed_Q"],
            ["1", "1/4", "64", "16"],
        )
        self.assertEqual(report["exact_results"]["filled_sea_reference_eigenvalue"], "64")
        self.assertTrue(report["raw_signed_route_receipt"]["state_ray_unique_from_CAR_projection"])
        self.assertEqual(report["exact_results"]["relative_charge_diag"], ["0", "1", "-1", "0"])
        self.assertEqual(report["exact_results"]["relative_H_diag"], ["0", "1/4", "3/4", "1"])
        self.assertEqual(
            report["exact_results"]["raw_filled_sea"],
            {"charge": "1", "energy": "-3/4", "quarter_phase": "i"},
        )
        self.assertEqual(
            report["exact_results"]["Q_ray_spectrum"],
            ["1", "1/4", "1/64", "1/256"],
        )
        self.assertEqual(
            report["exact_results"]["thermal_rho_diag"],
            ["65536/69649", "4096/69649", "16/69649", "1/69649"],
        )
        self.assertFalse(report["declared_source_inputs"]["absolute_partition_scalar_fixed"])
        self.assertFalse(report["declared_source_inputs"]["source_choice_from_P1_P2"])

    def test_charged_quarter_and_chronological_sewing_are_exact(self):
        report = verify_charged_pair()
        self.assertEqual(report["exact_results"]["charged_quarter_W"], ["1", "i", "-i", "1"])
        self.assertTrue(report["exact_results"]["W_equals_exp_2pi_i_H_on_source_Fock_words"])
        self.assertEqual(
            report["exact_results"]["chronological_two_collar_Q"],
            ["1", "1/16", "1/4096", "1/65536"],
        )
        self.assertTrue(report["route_receipt"]["state_ray_unique_from_CAR_projection"])

    def test_existing_search_accepts_zero_floor_and_rejects_quarter_floor(self):
        fixed = charged_fixed_moments()
        self.assertIn(
            {"word": [2, -1, 2], "re": "1/64", "im": "0"},
            fixed,
        )
        self.assertIn(
            {"word": [0, 2, -1, 2, 0], "re": "1/256", "im": "0"},
            fixed,
        )
        accepted = run_moment_search(0)
        rejected = run_moment_search(Q_BASE)
        self.assertIn(accepted["solver_status"], ("optimal", "optimal_inaccurate"))
        self.assertIn(rejected["solver_status"], ("infeasible", "infeasible_inaccurate"))


if __name__ == "__main__":
    unittest.main()
