import json
import tempfile
import unittest
from pathlib import Path

import sympy as sp

from kernel_functional import (
    action_split,
    completion_report,
    coherent_trace,
    controls,
    exponential,
    kernel,
    matrix_from_kernel,
    quarter_operator,
    sew,
    variable,
)
from joint_source import add, multiply
from origin_operator_checks import fixed_field_ward_check
from source_sewing import source_lift


class CoherentKernelTests(unittest.TestCase):
    def test_identity_kernel_is_coherent_resolution_of_identity(self):
        for modes in (1, 2, 3):
            identity = sp.eye(1 << modes)
            quadratic = {}
            for j in range(modes):
                quadratic = add(quadratic, multiply(variable(j), variable(modes + j)))
            self.assertEqual(kernel(identity), exponential(quadratic, 2 * modes))

    def test_generic_even_two_mode_matrix_roundtrips(self):
        t = sp.Matrix([[2, 0, 0, 3], [0, 5, 7, 0],
                       [0, 11, 13, 0], [17, 0, 0, 19]])
        self.assertEqual(matrix_from_kernel(kernel(t), 2), t)
        self.assertEqual(coherent_trace(t), sp.trace(t))
        w = quarter_operator(2)
        self.assertEqual(coherent_trace(t, w), sp.trace(t * w))
        with self.assertRaises(ValueError):
            kernel(sp.Matrix([[0, 1], [0, 0]]))

    def test_noncommuting_chronological_sewing_matches_matrix_product(self):
        left = sp.eye(4)
        right = sp.eye(4)
        left[0, 3] = 1
        right[3, 0] = 1
        self.assertNotEqual(left * right, right * left)
        self.assertEqual(sew(left, right), kernel(left * right))
        self.assertEqual(sew(sp.eye(4), sp.eye(4)), kernel(sp.eye(4)))

    def test_gaussian_and_interacting_actions_traces_and_quarter_time(self):
        transfer = sp.diag(sp.Rational(1, 4), sp.Rational(1, 64))
        gaussian = source_lift(transfer)
        _, gaussian_action = action_split(gaussian)
        self.assertTrue(all(mask.bit_count() <= 2 for mask in gaussian_action))

        interacting = sp.diag(1, sp.Rational(1, 4), sp.Rational(1, 64), sp.Rational(1, 65536))
        scalar, action = action_split(interacting)
        quartic = [value for mask, value in action.items() if mask.bit_count() == 4]
        self.assertEqual(scalar, 1)
        self.assertEqual(quartic, [-sp.Rational(255, 65536)])

        quarter = quarter_operator(2, charges=(1, -1))
        h_alt = sp.diag(0, sp.Rational(1, 4), sp.Rational(3, 4), 2)
        ell = 8 * sp.log(2)
        self.assertEqual(interacting, sp.diag(*(sp.exp(-ell * h_alt[j, j]) for j in range(4))))
        self.assertEqual(sp.diag(*(sp.exp(2 * sp.pi * sp.I * h_alt[j, j]) for j in range(4))), quarter)
        self.assertEqual(sp.trace(interacting * quarter), sp.Rational(65537, 65536) + sp.I * sp.Rational(15, 64))
        sea_scalar = sp.Pow(2, sp.Rational(1, 12))
        scaled_body, scaled_action = action_split(sea_scalar * interacting)
        self.assertEqual(scaled_body, sea_scalar)
        self.assertEqual(scaled_action, action)
        self.assertEqual(sp.simplify(sea_scalar * sp.exp(-ell*h_alt[1, 1])
                                     - sp.exp(-ell*(h_alt[1, 1]-sp.Rational(1, 96)))), 0)
        self.assertTrue(controls()["quarter_time_equation_exact"])
        self.assertTrue(controls()["full_hamiltonian_loop_phase_exact"])

    def test_equal_vacuum_and_one_particle_columns_do_not_fix_extension(self):
        gaussian = source_lift(sp.diag(sp.Rational(1, 4), sp.Rational(1, 64)))
        interacting = sp.diag(1, sp.Rational(1, 4), sp.Rational(1, 64), sp.Rational(1, 65536))
        self.assertEqual(gaussian[:, :3], interacting[:, :3])
        self.assertNotEqual(kernel(gaussian), kernel(interacting))

    def test_singular_transfer_can_have_a_finite_kernel_log_chart(self):
        singular = sp.diag(1, 0)
        scalar, action = action_split(singular)
        self.assertEqual(scalar, 1)
        self.assertEqual(matrix_from_kernel(kernel(singular), 1), singular)
        self.assertFalse(action)
        with self.assertRaises(ValueError):
            action_split(sp.diag(0, 1))

    def test_source_clustering_certificates_are_reconstructed_from_exact_entries(self):
        root = Path(__file__).parents[1] / "seam-positive-kernel-search-20261015" / "runs" / "source-clustering"
        for name in ("min-trace.json", "max-trace.json"):
            source = root / name
            report = completion_report(source)
            payload = json.loads(source.read_text())
            original = sp.diag(*(sp.I**q for q in payload["original_internal_W"]["i_exponents"]))
            normalized = sp.diag(*(sp.I**q for q in payload["R"]["i_exponents"]))
            matrix = sp.diag(*(sp.Rational(entry["re"]) + sp.I*sp.Rational(entry["im"])
                               for entry in payload["T"]["entries"]))
            self.assertTrue(report["finite_matrix_checks"]["diagonal_hermitian_psd_checked_exactly"])
            self.assertTrue(report["finite_matrix_checks"]["sparse_dense_matrix_reconstruction"])
            self.assertTrue(report["roundtrip_exact"])
            self.assertFalse(report["physical_source_selection_proved"])
            self.assertEqual(len(report["kernel_coefficients"]), 32)
            self.assertEqual(len(report["input_file_sha256"]), 64)
            self.assertEqual(report["original_internal_quarter_trace"], str(sp.trace(matrix*original)))
            self.assertEqual(report["normalized_GNS_quarter_trace"], str(sp.trace(matrix*normalized)))
            self.assertEqual(report["quarter_inserted_trace_implementer"], "original_internal_W = i^N")
            if name == "min-trace.json":
                self.assertEqual(report["min_trace_closed_form"]["matrix_matches"], True)
                self.assertEqual(report["min_trace_closed_form"]["kernel_matches"], True)
                self.assertEqual(report["min_trace_closed_form"]["nilpotent_action_matches"], True)
                self.assertIn("not physical vertex selection", report["min_trace_closed_form"]["label"])

    def test_completion_adapter_rejects_normalized_quarter_mismatch(self):
        source = Path(__file__).parents[1] / "seam-positive-kernel-search-20261015" / "runs" / "source-clustering" / "min-trace.json"
        payload = json.loads(source.read_text())
        payload["R"]["i_exponents"][0] = (payload["R"]["i_exponents"][0] + 1) % 4
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad-quarter.json"
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "normalized R"):
                completion_report(path)

    def test_completion_adapter_rejects_bad_flags_types_and_parity(self):
        root = Path(__file__).parents[1] / "seam-positive-kernel-search-20261015" / "runs" / "source-clustering"
        payload = json.loads((root / "min-trace.json").read_text())
        cases = []
        bad_flag = json.loads(json.dumps(payload))
        bad_flag["ok"] = False
        cases.append(bad_flag)
        bad_scalar = json.loads(json.dumps(payload))
        bad_scalar["T"]["entries"][0]["re"] = "__import__('os')"
        cases.append(bad_scalar)
        odd_matrix = json.loads(json.dumps(payload))
        odd_matrix["T"]["entries"].append({"row": 0, "column": 1, "re": "1", "im": "0"})
        cases.append(odd_matrix)
        for case in cases:
            with self.subTest(case=case["ok"]):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "bad.json"
                    path.write_text(json.dumps(case))
                    with self.assertRaises(ValueError):
                        completion_report(path)

    def test_existing_h_h_prime_control_keeps_its_scope(self):
        result = fixed_field_ward_check()
        self.assertFalse(result["geometric_locality_of_the_many_body_conjugation_proved"])
        self.assertFalse(result["P1_P2_source_selection_proved"])


if __name__ == "__main__":
    unittest.main()
