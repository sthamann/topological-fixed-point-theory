import json
import unittest
from unittest.mock import patch

from integral_origin_checks import primitive_integral_completion, function_differential_typing
from origin_hypothesis_search import (CHECK_SCRIPTS, run_checks,
                                      validate_proposal_ids, ALLOWED)


class IntegralOriginTests(unittest.TestCase):
    def test_primitive_integral_and_actual_negative_controls(self):
        result = primitive_integral_completion()
        facts = result['facts']
        self.assertEqual(facts['representatives_checked'], 11216)
        self.assertEqual({r['n'] for r in facts['survivors']}, {4})
        self.assertEqual(facts['n4_inverse_paired_glue_count'], 2)
        self.assertEqual(facts['n16_nonprimitive_norm'], '8')
        self.assertEqual(facts['n5_original_norm'], '23/10')
        self.assertEqual(facts['n5_scaled_norm'], '46')
        self.assertEqual(facts['n5_doubled_braiding_norm'], '0')
        self.assertEqual(facts['n5_doubled_positive_norm'], '23/5')
        self.assertFalse(result['inference_claims'][0]['holds'])

    def test_divisor_and_real_residue_computation(self):
        result = function_differential_typing()
        self.assertEqual(result['facts']['principal_parts'], 4)
        self.assertEqual(result['facts']['h0_K_D'], 3)
        self.assertEqual(result['facts']['residue_at_1'], 1)
        self.assertEqual(result['facts']['residue_at_infinity'], -1)
        self.assertFalse(result['inference_claims'][0]['holds'])

    def test_fixed_command_injection_and_wrong_response_id_rejected(self):
        test = 'primitive-integral-completion'
        with patch('origin_hypothesis_search.subprocess.run') as run:
            with self.assertRaises(ValueError):
                run_checks({'tests':[test]}, {test:{'command':['python3','-c','print(1)']}})
            run.assert_not_called()
            command = ['python3','-B','-OO',CHECK_SCRIPTS[test],'--check',test]
            run.return_value.returncode = 0
            run.return_value.stdout = json.dumps({'id':'wrong','status':'pass'})
            run.return_value.stderr = ''
            with self.assertRaises(ValueError):
                run_checks({'tests':[test]}, {test:{'command':command}})

    def test_wrong_source_and_premise_and_test_ids_rejected(self):
        ident = 'function-differential-typing'
        catalog = {ident:{'premise_ids':['premise:p']}}
        valid = {'source_ids':['script:v189'], 'premises':['premise:p'], 'tests':[ident]}
        validate_proposal_ids(valid, {'script:v189'}, {'premise:p'}, catalog)
        for update in ({'source_ids':['script:wrong']}, {'premises':['finding:wrong']},
                       {'tests':['wrong']}, {'premises':[]}):
            with self.assertRaises(ValueError):
                validate_proposal_ids({**valid, **update}, {'script:v189'}, {'premise:p'}, catalog)

    def test_existing_and_new_claim_mappings_are_retained(self):
        self.assertTrue({'quillen-relative-normalization', 'theta-heat-shape-response',
                         'doubled-source-counting', 'variable-mark-completion',
                         'primitive-integral-completion', 'function-differential-typing'} <= set(ALLOWED))
        self.assertIn('quillen-relative-normalization', ALLOWED)
        self.assertEqual(ALLOWED['primitive-integral-completion'],
                         'primitive_integral_projection_allows_n_other_than_4')


if __name__ == '__main__':
    unittest.main()
