"""Exact controls for the two-face collar multifield witness."""
from __future__ import annotations

import sympy as sp
import unittest

from collar_control import (
    B,
    T_SCALAR,
    adversarial_interacting_transfer,
    build_control,
    build_source_choi_control,
    exterior_transfer,
    run_moment_search,
    run_thermal_moment_search,
    source_fixed_moments,
    thermal_fixed_moments,
    thermal_gns_transfer,
    verify_control,
)


class CollarMultifieldControlTests(unittest.TestCase):
    def test_exact_collar_graph_state_time_and_moments(self):
        control = build_control()
        report = verify_control()
        self.assertTrue(report["ok"])
        self.assertEqual(control.graph, sp.Matrix([[16, 4], [4, 1]]) / 17)
        self.assertEqual(control.graph ** 2, control.graph)
        self.assertEqual(control.one_particle_h, B * (sp.eye(2) - 2 * control.graph))
        self.assertEqual(
            control.transfer.eigenvals(),
            {sp.Integer(1): 1, sp.Rational(1, 4): 2, sp.Rational(1, 16): 1},
        )
        self.assertEqual(control.moment((0, -1, 0)), sp.Rational(1, 4))
        self.assertEqual(control.moment((0, 2, -1, 2, 0)), sp.Rational(1, 16))
        self.assertEqual(
            exterior_transfer(control.graph, sp.Rational(1, 2)) ** 2,
            exterior_transfer(control.graph, T_SCALAR),
        )
        self.assertFalse(report["quarter_scope"]["identified"])

    def test_source_graph_state_unvectorizes_to_two_field_choi_ray(self):
        control = build_source_choi_control()
        self.assertEqual(control.response, sp.eye(2) / 4)
        self.assertEqual(
            control.transfer,
            sp.diag(1, sp.Rational(1, 4), sp.Rational(1, 4), sp.Rational(1, 16)),
        )
        self.assertEqual(
            control.thermal_density,
            sp.diag(256, 16, 16, 1) / 289,
        )
        self.assertEqual(control.boundary_graph[:2, :2], sp.eye(2) * sp.Rational(16, 17))
        self.assertEqual(control.boundary_graph[2:, 2:], sp.eye(2) * sp.Rational(1, 17))
        self.assertTrue(control.route_receipt["state_ray_unique_from_CAR_projection"])
        self.assertFalse(control.route_receipt["operator_scalar_recovered"])

    def test_report_preserves_source_faces_and_parent_provenance(self):
        report = verify_control()
        self.assertEqual(report["schema"], "tfpt.collar-multifield-control.v2")
        self.assertIn("two declared identical lowest-mode copies",
                      report["declared_inputs"]["control_mode_selection"])
        faces = report["exact_results"]["source_boundary_face_covariance"]
        self.assertEqual(faces, {
            "input_dual_boundary_CAR": "(16/17) I_2",
            "output_boundary_CAR": "(1/17) I_2",
        })
        matrices = report["exact_matrices"]
        self.assertNotIn("exterior_transfer_Q", matrices)
        self.assertNotIn("one_particle_h", matrices)
        self.assertNotIn("H_spectrum", report["exact_results"])
        source_q = matrices["source_choi_ray_Q_Fock4"]
        self.assertEqual([source_q[i][i]["re"] for i in range(4)],
                         ["1", "1/4", "1/4", "1/16"])
        parent = matrices["optional_pure_graph_parent"]
        self.assertEqual(parent["source_derived"], False)
        self.assertNotEqual(parent["transfer_Q"], source_q)

    def test_thermal_gns_transfer_is_positive_but_not_an_effect(self):
        _, liouville, transfer, omega = thermal_gns_transfer()
        self.assertEqual(
            liouville.eigenvals(),
            {sp.Integer(0): 2, sp.Rational(1, 4): 1, sp.Rational(-1, 4): 1},
        )
        self.assertEqual(
            transfer.eigenvals(),
            {sp.Integer(1): 2, sp.Rational(1, 4): 1, sp.Integer(4): 1},
        )
        self.assertEqual(transfer * omega, omega)
        self.assertEqual(liouville * omega, sp.zeros(4, 1))
        self.assertEqual(transfer[1, 1], 4)
        report = verify_control()
        self.assertFalse(
            report["time_type_audit"]["thermal_GNS_Liouville"]["is_positive_contraction"]
        )

    def test_full_field_covariance_is_the_extra_multifield_input(self):
        control = build_control()
        first, second = control.quasiparticle_creators
        for creator in (first, second):
            self.assertEqual(control.transfer * creator, T_SCALAR * creator * control.transfer)
        self.assertEqual(
            control.transfer * second * first,
            T_SCALAR ** 2 * second * first * control.transfer,
        )
        adversarial = adversarial_interacting_transfer(control)
        grade_zero_one = [
            control.omega,
            first * control.omega,
            second * control.omega,
        ]
        for vector in grade_zero_one:
            self.assertEqual(adversarial * vector, control.transfer * vector)
        two_particle = second * first * control.omega
        self.assertEqual(control.transfer * two_particle, sp.Rational(1, 16) * two_particle)
        self.assertEqual(adversarial * two_particle, sp.Rational(1, 8) * two_particle)
        self.assertTrue(any(
            adversarial * creator != T_SCALAR * creator * adversarial
            for creator in (first, second)
        ))

    def test_source_choi_search_accepts_zero_floor_and_rejects_quarter_floor(self):
        fixed = source_fixed_moments()
        self.assertIn(
            {"word": [0, 2, -1, 2, 0], "re": "1/16", "im": "0"},
            fixed,
        )
        accepted = run_moment_search(0.0)
        rejected = run_moment_search(0.25)
        self.assertIn(accepted["solver_status"], ("optimal", "optimal_inaccurate"))
        self.assertIn(rejected["solver_status"], ("infeasible", "infeasible_inaccurate"))

    def test_thermal_gns_search_requires_its_declared_ceiling(self):
        self.assertIn(
            {"word": [0, -1, -1, 0], "re": "1", "im": "0"},
            thermal_fixed_moments(),
        )
        accepted = run_thermal_moment_search(4)
        wrong_effect_class = run_thermal_moment_search(1)
        self.assertIn(accepted["solver_status"], ("optimal", "optimal_inaccurate"))
        self.assertIn(
            wrong_effect_class["solver_status"],
            ("infeasible", "infeasible_inaccurate"),
        )


if __name__ == "__main__":
    unittest.main()
