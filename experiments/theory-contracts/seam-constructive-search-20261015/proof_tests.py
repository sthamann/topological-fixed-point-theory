"""Independent positive and mutated-negative controls for proof_check.py."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import proof_check as checker


GOOD = {
    "schema": "tfpt.crossed_defect.v1",
    "chi_exp": [3] * 8,
    "zeta_exp": 2,
    "kappa_exp": 2,
}


class CertificateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_positive_candidate_reaches_kernel_and_has_no_extra_axioms(self):
        receipt = checker.check_candidate(GOOD, self.root / "positive")
        self.assertTrue(receipt["ok"], receipt)
        self.assertEqual(receipt["stage"], "complete")
        self.assertTrue(receipt["signature_seen"])
        self.assertLessEqual(set(receipt["axioms"]), checker.ALLOWED_AXIOMS)
        self.assertTrue(all(row["passed"] for row in receipt["checks"]))

    def test_ordinary_local_character_mutation_fails_before_lean(self):
        mutant = {**GOOD, "chi_exp": [0] * 8}
        receipt = checker.check_candidate(mutant, self.root / "local-mutant")
        self.assertFalse(receipt["ok"])
        self.assertEqual(receipt["stage"], "exact_cpu")
        self.assertIn("character_square_generators", receipt["error"])

    def test_zeta_rephasing_is_accepted_as_gauge(self):
        mutant = {**GOOD, "zeta_exp": 0}
        receipt = checker.check_candidate(mutant, self.root / "zeta-mutant")
        self.assertTrue(receipt["ok"], receipt)
        self.assertEqual(receipt["candidate_class"]["zeta_gauge_exp"], 0)

    def test_all_four_orientation_and_p_lift_classes_reach_lean(self):
        cases = [
            (1, 0, 0, "plus_i"),
            (1, 1, 2, "plus_i"),
            (3, 2, 0, "minus_i"),
            (3, 3, 2, "minus_i"),
        ]
        for orientation, zeta, kappa, label in cases:
            with self.subTest(orientation=orientation, zeta=zeta, kappa=kappa):
                candidate = {
                    "schema": GOOD["schema"],
                    "chi_exp": [orientation] * 8,
                    "zeta_exp": zeta,
                    "kappa_exp": kappa,
                }
                receipt = checker.check_candidate(
                    candidate, self.root / f"class-{orientation}-{zeta}-{kappa}")
                self.assertTrue(receipt["ok"], receipt)
                self.assertEqual(receipt["candidate_class"]["orientation"], label)
                self.assertEqual(receipt["candidate_class"]["p_lift_exp"], kappa)

    def test_invalid_p_lift_phase_fails_square(self):
        mutant = {**GOOD, "kappa_exp": 1}
        receipt = checker.check_candidate(mutant, self.root / "kappa-mutant")
        self.assertFalse(receipt["ok"])
        self.assertIn("p_lift_square", receipt["error"])

    def test_schema_rejects_text_injection_and_extra_fields(self):
        mutant = {**GOOD, "theorem": "axiom target : False"}
        receipt = checker.check_candidate(mutant, self.root / "schema-mutant")
        self.assertFalse(receipt["ok"])
        self.assertEqual(receipt["stage"], "schema")
        self.assertIn("extra=['theorem']", receipt["error"])

    def test_generated_file_is_literals_and_fixed_proof_is_byte_identical(self):
        receipt = checker.check_candidate(GOOD, self.root / "isolation")
        self.assertTrue(receipt["ok"], receipt)
        output = self.root / "isolation"
        self.assertEqual((output / "FixedCertificate.lean").read_bytes(),
                         checker.FIXED_CERTIFICATE.read_bytes())
        data = (output / "CandidateData.lean").read_text()
        self.assertNotRegex(data, r"\b(axiom|constant|sorry|admit|unsafe|opaque|partial)\b")
        self.assertNotIn("theorem", data)

    def test_mutated_fixed_target_is_rejected_by_lean(self):
        output = self.root / "mutated-target"
        output.mkdir()
        (output / "CandidateData.lean").write_text(checker.render_candidate_data(GOOD))
        source = checker.FIXED_CERTIFICATE.read_text()
        old = "epsExp q n - epsExp n q = 2 * (totalCharge n : Phase)"
        new = "epsExp q n - epsExp n q = 2 * (totalCharge n : Phase) + 1"
        self.assertIn(old, source)
        (output / "FixedCertificate.lean").write_text(source.replace(old, new, 1))
        env = checker.lean_environment(output)
        data_run = subprocess.run([str(checker.LEAN), "-o", "CandidateData.olean",
                                   "CandidateData.lean"],
                                  cwd=output, env=env, text=True,
                                  capture_output=True, timeout=120, check=False)
        self.assertEqual(data_run.returncode, 0, data_run.stdout + data_run.stderr)
        run = subprocess.run([str(checker.LEAN), "FixedCertificate.lean"],
                             cwd=output, env=env, text=True,
                             capture_output=True, timeout=120, check=False)
        self.assertNotEqual(run.returncode, 0, run.stdout + run.stderr)

    def test_cli_writes_machine_readable_failure_receipt(self):
        candidate = self.root / "candidate.json"
        candidate.write_text(json.dumps({**GOOD, "chi_exp": [0] * 8}))
        output = self.root / "cli"
        run = subprocess.run(["python3", str(Path(checker.__file__)),
                              "--candidate", str(candidate), "--output", str(output)],
                             text=True, capture_output=True, timeout=120, check=False)
        self.assertEqual(run.returncode, 1)
        receipt = json.loads((output / "receipt.json").read_text())
        self.assertFalse(receipt["ok"])
        self.assertEqual(json.loads(run.stdout), receipt)


if __name__ == "__main__":
    unittest.main()
