"""Independent controls for the optimizer-free exact m=5 proof checker."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from proof_check import check_certificate, check_pair


ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "runs/certified-positive"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text())


class OfflineProofCheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.moments = _load("moments.json")
        cls.minimum = _load("min-trace.json")
        cls.maximum = _load("max-trace.json")
        cls.next_source_question = _load("next-source-question.json")

    def test_both_real_completions_replay_and_are_separated(self):
        result = check_pair(
            self.minimum, self.maximum, self.moments, self.next_source_question
        )
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["minimum"]["moment_count"], 2245)
        self.assertEqual(result["maximum"]["moment_count"], 2245)
        self.assertEqual(result["pair_check"]["separating_degree5_words"], 10)
        self.assertTrue(result["pair_check"]["all_declared_degree5_values_differ"])
        self.assertTrue(result["pair_check"]["next_source_question_replayed"])
        self.assertTrue(result["pair_check"]["minimum_finite_clustering"])
        self.assertFalse(result["pair_check"]["maximum_finite_clustering"])

    def test_false_clustering_claim_on_old_maximum_is_rejected(self):
        mutated = copy.deepcopy(self.maximum)
        mutated["finite_clustering"] = True
        result = check_certificate(mutated, self.moments)
        self.assertFalse(result["ok"])
        self.assertIn("rank(I-T)", result["error"])

    def test_source_clustering_run_has_exact_gap_and_one_stationary_direction(self):
        source_run = ROOT / "runs/source-clustering"
        moments = json.loads((source_run / "moments.json").read_text())
        minimum = json.loads((source_run / "min-trace.json").read_text())
        maximum = json.loads((source_run / "max-trace.json").read_text())
        question = json.loads((source_run / "next-source-question.json").read_text())
        result = check_pair(minimum, maximum, moments, question)
        self.assertTrue(result["ok"], result)
        for completion in (result["minimum"], result["maximum"]):
            self.assertEqual(completion["cluster_gap"], "1/4")
            self.assertTrue(completion["finite_clustering"])
            self.assertEqual(completion["stationary_dimension"], 1)
            self.assertIsNotNone(
                completion["positivity"]["cluster_upper_bound_minus_T"]
            )

    def test_mutated_transfer_fails_exactly(self):
        mutated = copy.deepcopy(self.minimum)
        mutated["exact_matrices"]["T"][0][0] = {"re": [0, 1], "im": [0, 1]}
        result = check_certificate(mutated, self.moments)
        self.assertFalse(result["ok"])
        self.assertIn("positive", result["error"])

    def test_mutated_omega_fails_exactly(self):
        mutated = copy.deepcopy(self.minimum)
        nonzero = next(
            row for row in mutated["exact_matrices"]["Omega"]
            if row[0]["re"] != [0, 1] or row[0]["im"] != [0, 1]
        )
        nonzero[0] = {"re": [0, 1], "im": [0, 1]}
        result = check_certificate(mutated, self.moments)
        self.assertFalse(result["ok"])
        self.assertIn("normalized", result["error"])


if __name__ == "__main__":
    unittest.main()
