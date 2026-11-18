"""Regressions for model proposal boundaries, not tests of physical truth."""
import unittest
from unittest.mock import patch
from origin_hypothesis_search import premise_ids, run_checks, CHECK_SCRIPT, scoped_predictions_match


class ProposalBoundaries(unittest.TestCase):
    def test_unusable_expression_does_not_count_as_matching_candidate(self):
        entry = {'prediction_checks': [{'matches': True}],
                 'source_candidate_assessment': {'verdict': 'typing_error'}}
        self.assertFalse(scoped_predictions_match(entry))
        entry['source_candidate_assessment']['verdict'] = 'needs_original_source'
        self.assertTrue(scoped_predictions_match(entry))

    def test_findings_are_not_premise_ids(self):
        data = {'program': {'premises': [{'id': 'premise:given'}]},
                'ledger': [{'claim_id': 'AX.P1.01'}],
                'findings': [{'id': 'finding:known'}],
                'registry': [{'script': 'v1_e8_glue'}]}
        ids = premise_ids(data)
        self.assertIn('premise:given', ids)
        self.assertIn('AX.P1.01', ids)
        self.assertNotIn('finding:known', ids)
        self.assertNotIn('v1_e8_glue', ids)

    def test_curated_file_cannot_replace_fixed_command(self):
        catalog = {'quillen-relative-normalization': {'command': ['python3', '-c', 'raise SystemExit(0)']}}
        with patch('origin_hypothesis_search.subprocess.run') as process:
            with self.assertRaises(ValueError):
                run_checks({'tests': ['quillen-relative-normalization']}, catalog)
            process.assert_not_called()

    def test_mismatched_checker_response_rejected(self):
        key = 'quillen-relative-normalization'
        catalog = {key: {'command': ['python3', '-B', '-OO', CHECK_SCRIPT, '--check', key]}}
        with patch('origin_hypothesis_search.subprocess.run') as process:
            process.return_value.returncode = 0
            process.return_value.stdout = '{"id":"other","status":"pass"}'
            process.return_value.stderr = ''
            with self.assertRaises(ValueError):
                run_checks({'tests': [key]}, catalog)


if __name__ == '__main__':
    unittest.main()
