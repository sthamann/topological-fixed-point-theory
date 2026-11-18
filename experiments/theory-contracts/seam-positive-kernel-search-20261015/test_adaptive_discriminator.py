import unittest
from adaptive_discriminator import import_hypotheses, rank_responses


class AdaptiveDiscriminatorTests(unittest.TestCase):
    def test_affine_duplicates_are_one_signal_and_cheapest_is_selected(self):
        rows = [{'word': [0, 1, -1, 1, 0], 'values': ['1/4', '3/4'], 'cost': 5},
                {'word': [2], 'values': ['3/2', '5/2'], 'cost': 2},
                {'word': [3], 'values': ['7', '7'], 'cost': 1}]
        result = rank_responses(rows, 2)
        self.assertEqual(result['affine_signal_rank'], 1)
        self.assertEqual(result['affine_equivalence_groups'], 1)
        self.assertEqual(result['ranked_questions'][0]['word'], [2])
        self.assertEqual(len(result['ranked_questions'][0]['affine_equivalent_words']), 2)
        self.assertEqual(result['candidate_pairs_unseparated'], 0)

    def test_nonduplicate_linear_dependence_is_not_counted_as_rank(self):
        rows = [{'word': [0], 'values': [0, 1, 0]},
                {'word': [1], 'values': [0, 0, 1]},
                {'word': [2], 'values': [0, 1, 1]}]
        result = rank_responses(rows, 3)
        self.assertEqual(result['affine_equivalence_groups'], 3)
        self.assertEqual(result['affine_signal_rank'], 2)
        self.assertEqual(result['ranked_questions'][0]['new_candidate_pairs_separated'], 2)
        self.assertEqual(result['candidate_pairs_unseparated'], 0)

    def test_nonpositive_cost_and_mismatched_carrier_count_rejected(self):
        with self.assertRaises(ValueError):
            rank_responses([{'word': [0], 'values': [0, 1], 'cost': 0}], 2)
        with self.assertRaises(ValueError):
            rank_responses([{'word': [0], 'values': [0, 1]}], 3)

    def test_ai_proposals_are_unverified_and_never_constraints(self):
        proposals = [{'hypothesis': 'Source predicts a relation', 'source_ids': ['claim:AX.P1.01'],
                      'premises': ['marked source'], 'additional_assumptions': [],
                      'tests': ['compare source words']}]
        result = import_hypotheses(proposals)
        self.assertFalse(result['added_as_constraints'])
        self.assertFalse(result['source_ids_verified'])
        with self.assertRaises(ValueError):
            import_hypotheses([dict(proposals[0], source_ids=[])])


if __name__ == '__main__':
    unittest.main()
