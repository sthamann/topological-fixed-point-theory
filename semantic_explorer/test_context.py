import json
import tempfile
import unittest
from pathlib import Path

from semantic_explorer import build_context as bc
from build_context import assess,topic_role,merge_original_contract

class ContextTests(unittest.TestCase):
    def test_reviewed_classification_is_bound_to_original_row(self):
        row={'claim_id':'E8.FIXTURE.01','claim':'Finite identity plus failed identification',
             'status':'Exact [E] identities + [E, KILL] false dictionary','canonical_status':'Finite result; narrow kill','active':'true'}
        review={'row_sha256':bc.assessment_fingerprint(row),'category':'mixed','reason_de':'Getrennte Aussagen'}
        result=bc.assess(row,{row['claim_id']:review})
        self.assertEqual(result['category'],'mixed')
        self.assertEqual(result['scope'],row['canonical_status'])
        self.assertEqual(result['classification_method'],'source_reviewed_row_mapping')
        changed={**row,'claim':'Changed source proposition'}
        stale=bc.assess(changed,{row['claim_id']:review})
        self.assertFalse(stale['classification_review']['current'])
        self.assertTrue(stale['classification_uncertainty'])
        self.assertNotEqual(stale['classification_method'],'source_reviewed_row_mapping')
        self.assertEqual(bc.assess({**row,'active':'false'},{row['claim_id']:review})['category'],'historical')

    def row(self,status,canon='',active='true'):
        return {'status':status,'canonical_status':canon,'active':active}
    def test_partial_result_does_not_close_parent(self):
        a=assess(self.row('Identity (analytic theorem; fixed regulator)','[E] narrow norm theorem only; all T1-T8 remain [O]'))
        self.assertEqual(a['category'],'recorded_established')
        self.assertIn('T1-T8',a['scope'])
    def test_active_firewall_is_a_methodological_rule_not_an_unknown_result(self):
        status='Firewall (standing discipline rule, in force; audit census recorded)'
        a=assess(self.row(status))
        self.assertEqual(a['category'],'methodological_rule')
        self.assertFalse(a['classification_uncertainty'])
        self.assertEqual(assess(self.row(status,active='false'))['category'],'historical')
        self.assertEqual(assess(self.row('Firewall check remains open'))['category'],'open')
    def test_mixed_result_keeps_open_identification(self):
        self.assertEqual(assess(self.row('Identity [E]; physical identification [O]'))['category'],'mixed')
    def test_unproved_remainder_keeps_status_mixed(self):
        for status in ['Measured finite PSD + exact [E] rule; persistence unproven', 'Measured martingale hypotheses; normality step not proven']:
            self.assertIn(assess(self.row(status))['category'], ['mixed', 'open'])
    def test_axiom_and_inactive_not_solved(self):
        self.assertEqual(assess(self.row('Axiom (declared input)'))['category'],'assumed')
        self.assertEqual(assess(self.row('Identity [E]',active='false'))['category'],'historical')
    def test_numeric_is_scoped_reported_not_physical_completion(self):
        a=assess(self.row('Numerical'));self.assertEqual(a['category'],'recorded_established')
        self.assertEqual(a['authority'],'ledger');self.assertIn('keine neue',a['notice'])
    def test_real_ledger_marker_variants_and_measurements(self):
        for status in ['Theorem-grade [E neu] with reading [C neu]', 'Free-fermion invertibility [E/C]', 'Measured evidence; two control FAILs']:
            self.assertEqual(assess(self.row(status))['category'],'mixed')
        self.assertEqual(assess(self.row('Conceptual (verified)'))['category'],'recorded_established')
        self.assertEqual(assess(self.row('Not proven'))['category'],'open')
        self.assertEqual(assess(self.row('Kill-test [X]'))['category'],'unknown')
    def test_topic_focus_does_not_erase_cross_topic_connection(self):
        row={'claim_id':'PRIME.FLOOR.RATIO.01','claim':'Prime floor contract; related E8 theta coefficients'}
        self.assertEqual(topic_role(row,'topic:e8'),'related')
        self.assertEqual(topic_role(row,'topic:rh'),'primary')
        self.assertEqual(topic_role({'claim_id':'E8.SHELL.GLOBAL.01','claim':'Shell law'},'topic:e8'),'primary')

    def test_compound_markers_keep_every_qualification(self):
        for status in ['[E/O]', '[O/E]', '[E, C]', '[O | N]']:
            self.assertEqual(assess(self.row(status))['category'],'mixed')

    def test_canonical_result_scope_and_scoped_exclusions(self):
        self.assertEqual(assess(self.row(
            'Identity (quadratic target only; global zero-mode prescription required)',
            '[C] radiative target after global zero-mode removal; T7 remains [O]'))['category'],'conditional')
        self.assertEqual(assess(self.row(
            'Honest [E, kill] dynamic-cone dead verdict + [MEASURED] leaf-functional separation side finding',
            'CONE-DYNAMICS-DEAD; static cone untouched'))['category'],'mixed')
        self.assertEqual(assess(self.row('Honest [E, kill] dynamic-cone dead verdict'))['category'],'excluded')
        self.assertEqual(assess(self.row(
            'Proof-near [E] sector lemma with named residues',
            'L1-SECTOR-PROVEN-MODULO-ELEMENTARY; (b)/(c) stay open'))['category'],'mixed')

    def test_original_contract_replaces_stale_parsed_findings(self):
        old={'leaves_open':[], 'kills_free':['truncated...'], 'summary_de':'old', 'findings':[{'text':'historical indexed finding'}]}
        original={'verdict':'PARTIAL','leaves_open':['same-source geometry'], 'kills':['only the selected class is excluded'], 'summary_de':'current source summary', 'scope':'finite only'}
        merged=merge_original_contract(old,original)
        self.assertEqual(merged['leaves_open'],original['leaves_open'])
        self.assertEqual(merged['kills_free'],original['kills'])
        self.assertEqual(merged['summary_de'],original['summary_de'])
        self.assertEqual(merged['scope'],'finite only')
        self.assertEqual(merged['findings'],old['findings'])
        self.assertEqual(old['leaves_open'],[])

    def test_context_builder_uses_native_inventory_without_legacy_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); atlas=root/'atlas'; (atlas/'data').mkdir(parents=True)
            (root/'verification').mkdir(); (root/'verification/status_ledger.csv').write_text(
                'claim_id,claim,status,canonical_status,active\nCLAIM.ONE.01,First claim,Open,[O],true\n'
            )
            (atlas/'orientation.json').write_text(json.dumps({'topics':[
                {'id':'topic:test','label':'Test','summary':'Question','keywords':[]},
                {'id':'topic:e8','label':'E8','summary':'E8 topic','keywords':[]}
            ]}))
            contract=root/'experiments/theory-contracts/bridge'; contract.mkdir(parents=True)
            (contract/'contract_index.json').write_text(json.dumps({
                'contract':'bridge','question':'A scoped test','verdict':'PARTIAL',
                'topics':['topic:e8'],
                'claims_tested':['CLAIM.ONE.01'],
                'leaves_open':['The local derivation remains open.']
            }))
            (contract/'PROOF.txt').write_text('Original source excerpt.')
            (atlas/'data/knowledge.json').write_text(json.dumps({'nodes':[{
                'id':'claim:CLAIM.ONE.01','type':'claim','label':'CLAIM.ONE.01','summary':'',
                'attrs':{'semantic':{'topic_ids':['topic:test']}}
            }],'edges':[]}))
            old_here=bc.HERE; bc.HERE=atlas
            try:
                context=bc.build(root)
            finally:
                bc.HERE=old_here
            by_id={node['id']:node for node in context['nodes']}
            self.assertIn('contract:bridge',by_id)
            self.assertEqual(by_id['contract:bridge']['attrs']['leaves_open'],[
                'The local derivation remains open.'
            ])
            self.assertTrue(any(edge['src']=='contract:bridge' and edge['dst']=='claim:CLAIM.ONE.01' and edge['type']=='tests' for edge in context['edges']))
            self.assertTrue(any(edge['src']=='contract:bridge' and edge['dst']=='topic:e8' and edge['type']=='about' and edge['attrs'].get('source')=='contract_index_topics' for edge in context['edges']))
            self.assertTrue(context['sources'])
            self.assertFalse((root/'verification/theory_graph.json').exists())

if __name__=='__main__':unittest.main()
