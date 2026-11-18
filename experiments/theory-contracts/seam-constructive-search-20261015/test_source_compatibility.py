import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from source_compatibility import (
    DOMAIN_SPEC,
    SOURCE_PATH,
    SourceEvidenceError,
    assess_source_compatibility,
)


class SourceCompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.candidate = {
            "schema": "tfpt.crossed_defect.v1",
            "chi_exp": [3] * 8,
            "zeta_exp": 0,
            "kappa_exp": 2,
        }
        candidate_payload = json.dumps(self.candidate, sort_keys=True, separators=(",", ":"))
        self.receipt = {
            "ok": True,
            "candidate_sha256": hashlib.sha256(candidate_payload.encode()).hexdigest(),
        }

    def test_current_fixed_class_is_excluded_only_as_unchanged_source(self):
        result = assess_source_compatibility(self.candidate, self.receipt)
        self.assertTrue(result["algebra_certificate_ok"])
        self.assertEqual(result["source_equivalence_status"], "EXCLUDED_IN_UNCHANGED_SOURCE")
        self.assertEqual(result["half_charge_operator_obstruction"]["status"], "APPLIES")
        self.assertEqual(result["ground_space"]["original_source_multiplicity"], 1)
        self.assertEqual(result["ground_space"]["stipulated_e8_extension_multiplicity"], 2)
        for branch in result["ground_space"]["rational_connection_check"]:
            self.assertEqual(branch["K_zero"], "0")
            self.assertEqual(branch["K_r0"], "0")
            self.assertEqual(branch["delta_dot_r0"], "-1")

    def test_changed_dimension_or_time_is_not_given_the_source_verdict(self):
        base = json.loads(DOMAIN_SPEC.read_text())
        mutations = []
        changed_dimension = copy.deepcopy(base)
        changed_dimension["fixed_source_data"]["dimension"] = 10
        mutations.append(changed_dimension)
        changed_time = copy.deepcopy(base)
        changed_time["fixed_source_data"]["time"]["K_alpha"] = "K_prime(p,N)"
        mutations.append(changed_time)

        with tempfile.TemporaryDirectory() as directory:
            for index, spec in enumerate(mutations):
                with self.subTest(index=index):
                    path = Path(directory) / f"domain-{index}.json"
                    path.write_text(json.dumps(spec))
                    result = assess_source_compatibility(
                        self.candidate, self.receipt, domain_spec_path=path
                    )
                    self.assertEqual(result["source_equivalence_status"], "UNKNOWN")
                    self.assertEqual(
                        result["half_charge_operator_obstruction"]["status"], "NOT_TRANSFERRED"
                    )
                    self.assertEqual(result["ground_space"]["rational_connection_check"], [])
                    self.assertIsNone(result["ground_space"]["original_source_multiplicity"])
                    self.assertIsNone(result["ground_space"]["stipulated_e8_extension_multiplicity"])

    def test_changed_obstruction_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            changed = Path(directory) / "PROOF.txt"
            changed.write_bytes(SOURCE_PATH.read_bytes() + b"\nchanged\n")
            with self.assertRaises(SourceEvidenceError):
                assess_source_compatibility(self.candidate, self.receipt, source_path=changed)


if __name__ == "__main__":
    unittest.main()
