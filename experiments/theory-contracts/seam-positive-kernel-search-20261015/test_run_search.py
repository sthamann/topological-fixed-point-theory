import unittest
from pathlib import Path
import tempfile
from run_search import discriminate, campaign
from word_algebra import WordAlgebra


class DiscriminatorTests(unittest.TestCase):
    def test_scaled_reference_campaign_is_inconclusive_and_cache_binds_lambda(self):
        from unittest.mock import patch
        controls=[{'word':[-1],'re':'1/16'},{'word':[0,1],'re':'0','im':'1'},
                  {'word':[0,-1,0],'re':'1/64'}]
        with tempfile.TemporaryDirectory() as tmp,patch('run_search.complete',side_effect=AssertionError('cannot recover scaled source')):
            out=Path(tmp)/'scaled'
            result=campaign(out,m=1,degree=2,transfer_kind='bounded_positive',
                            transfer_ceiling=1/16,reference_eigenvalue=1/16,fixed_moments=controls)
            self.assertEqual(result['configuration']['reference_eigenvalue'],1/16)
            self.assertEqual(result['certified_finite_completions'],0)
            self.assertTrue(result['conditional_normalized_state'])
            self.assertFalse(result['absolute_partition_function_recovered'])
            self.assertTrue(all('reference_eigenvalue1' in c['reason'] for c in result['completions']))
            with self.assertRaisesRegex(ValueError,'different code/configuration'):
                campaign(out,m=1,degree=2,transfer_kind='bounded_positive',
                         transfer_ceiling=1/16,reference_eigenvalue=1/32,fixed_moments=controls)

    def test_thermal_campaign_never_calls_pure_Fock_completion(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp,patch('run_search.complete',side_effect=AssertionError('wrong completion type')):
            result=campaign(Path(tmp)/'thermal',m=1,degree=2,transfer_kind='thermal_gns',transfer_ceiling=4)
            self.assertEqual(result['configuration']['transfer_ceiling'],4)
            self.assertEqual(result['certified_finite_completions'],0)
            self.assertTrue(all(item['verdict']=='INCONCLUSIVE_UNSUPPORTED_TRANSFER_TYPE' for item in result['completions']))
            self.assertTrue(all('thermal GNS' in item['reason'] for item in result['completions']))
            self.assertTrue(result['inconclusive_completion_is_not_source_exclusion'])

    def test_default_campaign_is_feasibility_and_mixed_completion_stays_inconclusive(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'fresh'
            result=campaign(out,m=1,degree=2)
            self.assertEqual(result['configuration']['objective'],'feasibility')
            self.assertEqual(result['configuration']['transfer_floor'],0)
            self.assertEqual(result['configuration']['cluster_gap'],0)
            self.assertEqual(result['certified_finite_completions'],0)
            self.assertEqual(result['verdict'],'SEARCH_INCONCLUSIVE_NO_CERTIFIED_COMPLETION')
            self.assertTrue(all('not pure' in item['reason'] for item in result['completions']))
            self.assertIn('minimal irreducible',result['completion_ansatz']['representation'])
            self.assertTrue(result['inconclusive_completion_is_not_source_exclusion'])
            self.assertEqual(campaign(out,m=1,degree=2)['campaign_identity'],result['campaign_identity'])
            with self.assertRaises(ValueError):campaign(out,m=1,degree=2,objective='diversity')

    def test_orphaned_historical_artifacts_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);historical=out/'moments.json';historical.write_text('historical')
            with self.assertRaisesRegex(ValueError,'fresh output directory'):campaign(out,m=1,degree=2)
            self.assertEqual(historical.read_text(),'historical')

    def candidates(self):
        word=[0,2,-1,2,0]
        return [{'ok':True,'input_sha256':'same',
                 'candidate_predictions':[{'word':word,'value':{'re':value,'im':'0'}}]}
                for value in ['1/4','1']]

    def test_shortest_requires_complete_lower_word_data(self):
        first,second=self.candidates()
        payload={'m':2,'degree':2,'moments':[{'word':list(w)} for w in WordAlgebra(2).words(4)]}
        result=discriminate(first,second,payload)
        self.assertTrue(result['shortest_certified'])
        payload['moments'].pop()
        self.assertFalse(discriminate(first,second,payload)['shortest_certified'])

    def test_different_input_data_not_compared(self):
        first,second=self.candidates()
        second['input_sha256']='different'
        with self.assertRaises(ValueError): discriminate(first,second,{})

    def test_failed_completion_never_becomes_discriminator(self):
        first,second=self.candidates();second['ok']=False
        self.assertEqual(discriminate(first,second,{})['status'],'NO_TWO_CERTIFIED_COMPLETIONS')


if __name__=='__main__': unittest.main()
