import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import manage
import problem


class CollectorTests(unittest.TestCase):
    def test_exhaustive_only_survivor_reaches_certificate(self):
        # A candidate missed by genetics must still reach Lean after enumeration.
        row = {"key": "exhaustive-only", "genotype": problem.exhaustive_generator_classification()[0],
               "result": {"feasible": True}}
        snap = {"0-status": {"reason": "classification_complete", "total_evaluated_here": 128},
                "1-status": {"reason": "classification_complete", "total_evaluated_here": 128},
                "classification": {"complete": True, "classified_unique": 262144, "survivors": [row]}}
        with tempfile.TemporaryDirectory() as directory:
            run = {"hours": 1, "local_run": directory}
            with patch.object(manage, "verify_inputs"), patch.object(manage, "snapshot", return_value=snap), \
                 patch("proof_check.check_candidate", return_value={"ok": True}) as proof:
                manage.collect(run)
            state = json.loads((Path(directory) / "proof_state.json").read_text())
            self.assertEqual(state["status"], "finite_class_finished")
            self.assertEqual(state["certified_candidate_count"], 1)
            self.assertEqual(state["classified_unique"], 262144)
            self.assertFalse(state["physical_solution"])
            proof.assert_called_once()

    def test_changed_source_stops_owned_run_before_proof(self):
        with tempfile.TemporaryDirectory() as directory:
            run = {"hours": 1, "local_run": directory, "remote_run": "/owned/test"}
            with patch.object(manage, "verify_inputs", side_effect=ValueError("source changed")), \
                 patch.object(manage, "ssh_script") as stop, patch("proof_check.check_candidate") as proof:
                manage.collect(run)
            state = json.loads((Path(directory) / "proof_state.json").read_text())
            self.assertEqual(state["status"], "source_changed_stopped")
            proof.assert_not_called()
            self.assertIn("/owned/test", stop.call_args.args[0])

    def test_dead_workers_do_not_look_like_ongoing_search(self):
        snap = {"0-status": {"reason": "running"}, "1-status": {"reason": "running"},
                "launch": {"workers": [{"process_exists": False}, {"process_exists": False}]}}
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(manage, "verify_inputs"), patch.object(manage, "snapshot", return_value=snap):
                manage.collect({"hours": 1, "local_run": directory})
            state = json.loads((Path(directory) / "proof_state.json").read_text())
            self.assertEqual(state["status"], "workers_exited_without_complete_classification")


if __name__ == "__main__":
    unittest.main()
