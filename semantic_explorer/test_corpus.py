"""Focused tests for deterministic, provenance-preserving corpus extraction."""
import json
import hashlib
import tempfile
import unittest
from pathlib import Path

from semantic_explorer import build_corpus as bc
from semantic_explorer import source_inventory


class CorpusTests(unittest.TestCase):
    def test_declared_checkers_with_equal_basenames_keep_distinct_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory).resolve();old=bc.ROOT;bc.ROOT=root
            try:
                sources=[];nodes=[];identities=[]
                for name in ['a','b']:
                    folder=root/'experiments/theory-contracts'/name;folder.mkdir(parents=True)
                    checker=folder/'checker.py';checker.write_text('raise RuntimeError("must never execute")\n')
                    resolved=bc.resolve_declared_script('checker.py',folder/'contract_index.json',root/'experiments')
                    self.assertEqual(resolved,checker)
                    self.assertEqual(bc.resolve_declared_script(checker.relative_to(root).as_posix(),folder/'contract_index.json',folder),checker)
                    identities.append(bc.add_declared_script(resolved,sources,nodes))
                self.assertNotEqual(*identities)
                self.assertEqual({n['attrs']['path'] for n in nodes},{'experiments/theory-contracts/a/checker.py','experiments/theory-contracts/b/checker.py'})
                self.assertTrue(all(n['attrs']['executed'] is False for n in nodes))
                self.assertIsNone(bc.resolve_declared_script('../../../../../outside.py',folder/'contract_index.json',root))
            finally:bc.ROOT=old

    def test_unregistered_verification_helper_gets_its_own_source_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory).resolve();folder=root/'verification';folder.mkdir()
            (folder/'script_registry.csv').write_text('script,cluster\nv1_registered,core\n')
            helper=folder/'checker.py';helper.write_text('raise RuntimeError("do not execute")\n')
            registered=folder/'v1_registered.py';registered.write_text('pass\n')
            old=bc.ROOT;bc.ROOT=root
            try:
                sources=[];nodes=[]
                identity=bc.add_declared_script(helper,sources,nodes)
                self.assertTrue(identity.startswith('experiment_script:'))
                self.assertEqual(nodes[0]['attrs']['path'],'verification/checker.py')
                self.assertEqual(bc.add_declared_script(registered,sources,nodes),'script:v1_registered')
            finally:bc.ROOT=old

    def test_contract_alias_resolves_only_from_exact_source_ownership(self):
        sources=[{'id':'s1','path':'experiments/theory-contracts/a/contract_index.json'},
                 {'id':'s2','path':'experiments/theory-contracts/b/results.json'},
                 {'id':'s3','path':'unrelated/results.json'}]
        nodes=[{'id':'contract:a'},{'id':'contract:b'}]
        edges=[{'src':'record:'+source['id'],'dst':'contract:SAME.SCIENTIFIC.ID','type':'references_contract','attrs':{'source_id':source['id']}} for source in sources]
        gaps=[];bc.resolve_contract_references(sources,nodes,edges,gaps)
        self.assertEqual([edge['dst'] for edge in edges],['contract:a','contract:b'])
        self.assertEqual(len(gaps),1)
        self.assertTrue(all(edge['attrs']['semantic'] is False for edge in edges))
        self.assertTrue(all(edge['attrs']['declared_contract_id']=='SAME.SCIENTIFIC.ID' for edge in edges))

    def test_tex_parser_preserves_structure_without_asserting_proof(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); (root/'main.tex').write_text('\\section{Alpha}\\label{sec:a}\n\\begin{theorem}\\label{thm:x}\nA = B.\\end{theorem}\n\\begin{equation}x=1\\end{equation}\n\\ref{sec:a}\\n')
            old=bc.ROOT; bc.ROOT=root
            try:
                sources=[]; nodes=[]; edges=[]; refs=[]; labels={}
                bc.parse_tex(root/'main.tex',sources,nodes,edges,refs,labels,set(),{})
                types={n['type'] for n in nodes}
                self.assertIn('paper_section',types); self.assertIn('theorem_environment',types); self.assertIn('formula_context',types)
                thm=next(n for n in nodes if n['type']=='theorem_environment')
                self.assertEqual(thm['attrs']['epistemic'],'source_text_only')
                self.assertFalse(any('proves' in e['type'] or 'supports' in e['type'] for e in edges))
                self.assertTrue(sources[0]['sha256']); self.assertGreaterEqual(sources[0]['line_count'],1)
            finally: bc.ROOT=old

    def test_include_closure_and_escape_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); (root/'a.tex').write_text('\\input{sub/b}\n\\input{../../outside}\n'); (root/'sub').mkdir(); (root/'sub/b.tex').write_text('content')
            old=bc.ROOT; bc.ROOT=root
            try:
                files,gaps=bc.tex_files(['a.tex'])
                self.assertEqual([p.name for p in files],['a.tex','b.tex'])
                self.assertTrue(any('outside' in x['reason'] for x in gaps))
            finally: bc.ROOT=old

    def test_deterministic_ids_and_language_heuristic(self):
        self.assertEqual(bc.sid('a.tex','abc'),bc.sid('a.tex','abc'))
        self.assertNotEqual(bc.sid('a.tex','abc'),bc.sid('b.tex','abc'))
        self.assertEqual(bc.detect_lang('notes_de.tex','The and this'), 'de')
        self.assertEqual(bc.detect_lang('notes.tex','der die und das'), 'de')

    def test_unique_id_validator_rejects_empty_and_duplicate_ids(self):
        with self.assertRaisesRegex(ValueError,'duplicate id'):
            bc.validate_unique_ids([{'id':'source:a'},{'id':'source:a'}],'sources')
        with self.assertRaisesRegex(ValueError,'empty'):
            bc.validate_unique_ids([{'id':''}],'nodes')

    def test_artifact_allowlist_excludes_logs(self):
        self.assertIn('results.json',bc.ALLOWED_EXPERIMENT_FILES)
        self.assertNotIn('run.log',bc.ALLOWED_EXPERIMENT_FILES)
        self.assertLessEqual(bc.MAX_ARTIFACT_BYTES,262144)

    def test_reviewed_mapping_edges_are_relevance_only(self):
        doc={'reviewed_scope':'fixture only','aspects':[{
            'id':'aspect:fixture','type':'research_question','label':'Open question',
            'summary_de':'Frage','summary_en':'Question','epistemic':'open','scope':'fixture',
            'topic_ids':['topic:origin'],'claim_ids':['CLAIM.ONE.01'],
            'mapping_rationale':'Related open claim only.',
            'source_refs':[{'path':'paper.tex','line_start':1,'line_end':2,'review_note':'fixture'}]
        }]}
        sources=[{'id':'source:one','path':'paper.tex','line_count':2}]
        nodes=[{'id':'section:one','type':'paper_section','attrs':{'path':'paper.tex','line_start':1,'line_end':2}},
               {'id':'source:one','type':'tex_source','attrs':{}},
               {'id':'claim:CLAIM.ONE.01','type':'claim','attrs':{}},
               {'id':'topic:origin','type':'topic','attrs':{}}]
        edges=[]
        bc.integrate_reviewed_aspects(doc,sources,nodes,edges,{'CLAIM.ONE.01'},{'topic:origin'})
        aspect=next(n for n in nodes if n['id']=='aspect:fixture')
        self.assertEqual(aspect['type'],'research_question')
        self.assertEqual(aspect['attrs']['extraction_method'],'source_reviewed_mapping')
        claim_edge=next(e for e in edges if e['type']=='associated_claim')
        self.assertFalse(claim_edge['attrs']['supports']); self.assertFalse(claim_edge['attrs']['proves'])
        self.assertFalse(claim_edge['attrs']['status_changed'])
        self.assertTrue(any(e['type']=='grounded_in_section' for e in edges))

    def test_finding_classifier_keeps_completed_old_questions_distinct_from_open_ones(self):
        completed = [
            'Damit ist diese zuvor offengelassene Möglichkeit konkret ausgewertet.',
            'Er entstand vor der hier ausgeführten kanonischen Ergänzung; seine damals offene algebraische Teilfrage wird in der Herleitung beantwortet.',
        ]
        still_open = 'Die erste offene Herkunftsfrage lautet jetzt genauer: Welche Regel wählt den Ursprung aus?'
        self.assertEqual([source_inventory._finding_kind(text) for text in completed], ['supports', 'supports'])
        self.assertEqual(source_inventory._finding_kind(still_open), 'open')
        extracted = source_inventory._extract_findings('\n\n'.join(completed + [still_open]), None)
        self.assertEqual([item['kind'] for item in extracted], ['supports', 'supports', 'open'])

    def test_readme_finding_keeps_negated_proof_open(self):
        # Actual README excerpt: seam-full-collar-reconstruction-20261015.
        excerpt = 'Die physische Auswahl dieser Quelle aus P1/P2 ist damit nicht bewiesen.'
        findings = source_inventory._extract_findings(excerpt, None)
        self.assertEqual([(item['kind'], item['text']) for item in findings], [('open', excerpt)])
        negatives = [
            'Die Quelle ist nicht aus P1/P2 hergeleitet.',
            'Die Quelle ist nicht vollständig bewiesen.',
            'The physical source is not proved.',
            'The physical source is not yet derived.',
            'The physical source has not yet been rigorously proved.',
        ]
        for text in negatives:
            with self.subTest(text=text):
                self.assertEqual(source_inventory._finding_kind(text), 'open')
                self.assertEqual(source_inventory._extract_findings(text, None)[0]['kind'], 'open')
        for text in [
            'Der bedingte Satz ist bewiesen.',
            'The conditional theorem is proved.',
            'Der Satz ist nicht nur bewiesen, sondern exakt bestätigt.',
            'The theorem is not only proved, but confirmed.',
            'This does not refute the theorem, which is proved conditionally.',
        ]:
            with self.subTest(text=text):
                self.assertEqual(source_inventory._finding_kind(text), 'supports')
        self.assertEqual(source_inventory._finding_kind('This does not refute the conditional theorem.'), 'other')

    def test_native_inventory_reads_sources_without_legacy_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve()
            contract=root/'experiments/theory-contracts/bridge'; contract.mkdir(parents=True)
            (contract/'contract_index.json').write_text(json.dumps({
                'contract':'bridge','question':'A scoped test','verdict':'PARTIAL',
                'scope':'Finite conditional fixture only.',
                'claims_tested':['CLAIM.ONE.01'],
                'gates':[],
                'topics':['topic:quantum','topic:clocks','topic:universalraum'],
                'leaves_open':['The local geometric action is not derived.','T1'],
                'kills':['The broad closure claim is excluded by this bounded test.'],
                'relations':[{'type':'followup_of','target':'prior-check'},
                             {'type':'uses','target':'prior-check'},
                             {'type':'extends','target':'prior-check'},
                             {'type':'refines','target':'prior-check'}]
            }))
            (root/'experiments/theory-contracts/prior-check').mkdir(parents=True)
            (root/'experiments/next.txt').write_text(
                '# 2026-01-02 bridge: follow up from prior-check\n'
            )
            (contract/'PROOF.txt').write_text('Source excerpt.')
            (contract/'README.md').write_text(
                '# Bridge\n\nThe local geometric action is not derived; T1 remains open.'
            )
            experiment=root/'experiments/finite-check'; experiment.mkdir(parents=True)
            (experiment/'README.md').write_text('# Finite check\nTests CLAIM.ONE.01 only.')
            (root/'experiments').mkdir(exist_ok=True)
            (root/'experiments/evidence_scorecard.json').write_text(json.dumps({'rows':[{
                'domain':'D','observable':'Observable','status':'consistent',
                'source':'finite-check/checker.py'
            }]}))
            rh=root/'rh/catalog'; rh.mkdir(parents=True)
            (rh/'rh_semantic_catalog.json').write_text(json.dumps({'records':[{
                'path':'rh/sample.py','outcome':'KILLED','mechanism':'finite counterexample',
                'ledger_ids':['CLAIM.ONE.01']
            }]}))
            sidecar=root/'verification/graph_v2'; sidecar.mkdir(parents=True)
            (sidecar/'failure_scopes.json').write_text(json.dumps({'records':[{
                'id':'FAIL.SAMPLE.01','target':'one bounded adapter class',
                'evidence':['contract:bridge'],'excluded_under':['finite fixture']
            }]}))
            inventory=source_inventory.build(root,{'CLAIM.ONE.01'})
            by_id={n['id']:n for n in inventory['nodes']}
            self.assertEqual(by_id['contract:bridge']['attrs']['leaves_open'],[
                'The local geometric action is not derived.', 'T1'
            ])
            self.assertEqual(by_id['contract:bridge']['attrs']['kills_free'],[
                'The broad closure claim is excluded by this bounded test.'
            ])
            self.assertEqual(by_id['contract:bridge']['attrs']['original_contract']['kills'],[
                'The broad closure claim is excluded by this bounded test.'
            ])
            self.assertEqual(by_id['contract:bridge']['attrs']['topics'],[
                'topic:quantum','topic:clocks','topic:universalraum'
            ])
            self.assertEqual(by_id['contract:bridge']['attrs']['scope'],'Finite conditional fixture only.')
            self.assertEqual(by_id['contract:bridge']['attrs']['prose_gate_mentions'],['T1'])
            self.assertTrue(any(e['src']=='contract:bridge' and e['dst']=='goal:T1' and e['type']=='leaves_open' for e in inventory['edges']))
            self.assertTrue(all(any(e['src']=='contract:bridge' and e['dst']==topic and e['type']=='about' for e in inventory['edges'])
                                for topic in ('topic:quantum','topic:clocks','topic:universalraum')))
            self.assertTrue(any(e['src']=='contract:bridge' and e['dst']=='contract:prior-check' and e['type']=='references' for e in inventory['edges']))
            self.assertTrue(any(e['src']=='contract:bridge' and e['dst']=='contract:prior-check' and e['type']=='followup_of' for e in inventory['edges']))
            self.assertTrue(all(any(e['src']=='contract:bridge' and e['dst']=='contract:prior-check' and e['type']==kind for e in inventory['edges']) for kind in ('uses','extends','refines')))
            self.assertTrue(any(e['src']=='contract:bridge' and e['dst']=='contract:prior-check' and e['type']=='references' and e['attrs'].get('source')=='next_txt' for e in inventory['edges']))
            self.assertEqual(by_id['failure:FAIL.SAMPLE.01']['attrs']['target'],'one bounded adapter class')
            self.assertIn('rh:rh/sample.py',by_id)
            self.assertIn('scorecard:D:observable',by_id)
            self.assertTrue(any(e['src']=='contract:bridge' and e['dst']=='claim:CLAIM.ONE.01' and e['type']=='tests' for e in inventory['edges']))
            self.assertTrue(any(e['src']=='failure:FAIL.SAMPLE.01' and e['dst']=='contract:bridge' and e['type']=='excludes' for e in inventory['edges']))
            self.assertTrue(any(e['src']=='scorecard:D:observable' and e['dst']=='experiment:finite-check' and e['type']=='sourced_from' for e in inventory['edges']))
            self.assertTrue(inventory['sources'])

    def test_contract_prose_cannot_expand_empty_gate_declaration(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve()
            contract=root/'experiments/theory-contracts/native-sector-calibration-20261013'; contract.mkdir(parents=True)
            (contract/'contract_index.json').write_text(json.dumps({
                'contract':'native-sector-calibration-20261013','verdict':'PARTIAL',
                'gates':[],'topics':['topic:quantum','topic:clocks','topic:universalraum']
            }))
            (contract/'README.md').write_text('# Calibration\n\nT1–T8 remain open.\n')
            inventory=source_inventory.build(root,set())
            node=next(n for n in inventory['nodes'] if n['id']=='contract:native-sector-calibration-20261013')
            self.assertEqual(node['attrs']['gates'],[])
            self.assertEqual(node['attrs']['leaves_open'],[])
            self.assertEqual(node['attrs']['prose_gate_mentions'],[f'T{i}' for i in range(1,9)])
            self.assertFalse(any(e['src']==node['id'] and e['type'] in {'leaves_open','closes'} and e['dst'].startswith('goal:')
                                 for e in inventory['edges']))
            self.assertTrue(all(any(e['src']==node['id'] and e['dst']==topic and e['type']=='about'
                                    for e in inventory['edges'])
                                for topic in ('topic:quantum','topic:clocks','topic:universalraum')))
            self.assertFalse((root/'verification/theory_graph.json').exists())

    def test_declared_proof_and_checker_are_registered_as_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); contract=root/'experiments/theory-contracts/register'; contract.mkdir(parents=True)
            (contract/'contract_index.json').write_text(json.dumps({
                'contract':'register','question':'A scoped check','verdict':'PARTIAL',
                'proof':'evidence/PROOF.md','checker':'checker.py'
            }))
            (contract/'evidence').mkdir()
            (contract/'evidence/PROOF.md').write_text('# Proof record\nA bounded derivation record.\n')
            (contract/'checker.py').write_text('def check():\n    return True\n')
            inventory=source_inventory.build(root,set())
            contract_node=next(n for n in inventory['nodes'] if n['id']=='contract:register')
            registered={s['path']:s['id'] for s in inventory['sources']}
            self.assertIn('experiments/theory-contracts/register/evidence/PROOF.md',registered)
            self.assertIn('experiments/theory-contracts/register/checker.py',registered)
            self.assertTrue(set(contract_node['attrs']['source_refs']).issuperset({
                registered['experiments/theory-contracts/register/evidence/PROOF.md'],
                registered['experiments/theory-contracts/register/checker.py'],
            }))

    def test_source_manifest_registers_original_dependency_and_hash_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); contract=root/'experiments/theory-contracts/native'; contract.mkdir(parents=True)
            dependency=root/'tfpt_explorer/neutral_source_response.py'; dependency.parent.mkdir(parents=True)
            dependency.write_text('def source_response():\n    return "native"\n')
            digest=hashlib.sha256(dependency.read_bytes()).hexdigest()
            (contract/'source_manifest.json').write_text(json.dumps({'sources':[{
                'path':'tfpt_explorer/neutral_source_response.py','sha256':digest
            }]}))
            (contract/'contract_index.json').write_text(json.dumps({
                'contract':'native','question':'A source-dependent check','source_manifest':'source_manifest.json'
            }))
            inventory=source_inventory.build(root,set())
            node=next(n for n in inventory['nodes'] if n['id']=='contract:native')
            dep_source=next(s for s in inventory['sources'] if s['path']=='tfpt_explorer/neutral_source_response.py')
            self.assertIn(dep_source['id'],node['attrs']['source_refs'])
            self.assertEqual(node['attrs']['source_manifest_entries'],[{
                'path':'tfpt_explorer/neutral_source_response.py','expected_sha256':digest,
                'actual_sha256':digest,'hash_matches':True,'source_ref':dep_source['id']
            }])

    def test_inventory_check_detects_changed_source_without_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); folder=root/'experiments/theory-contracts/bridge'; folder.mkdir(parents=True)
            source=folder/'contract_index.json'
            source.write_text(json.dumps({'contract':'bridge','question':'A test','leaves_open':['gap']}))
            inventory=source_inventory.build(root,set())
            corpus_path=root/'semantic_explorer/data/corpus.json'; corpus_path.parent.mkdir(parents=True)
            corpus_path.write_text(json.dumps({'sources':inventory['sources'],'nodes':[
                {'id':node['id']} for node in inventory['nodes']
            ]}))
            self.assertTrue(source_inventory.check_against_corpus(root,inventory,corpus_path)['ok'])
            source.write_text(json.dumps({'contract':'bridge','question':'Changed source','leaves_open':['gap']}))
            changed=source_inventory.build(root,set())
            result=source_inventory.check_against_corpus(root,changed,corpus_path)
            self.assertFalse(result['ok'])
            self.assertIn('experiments/theory-contracts/bridge/contract_index.json',result['stale_sources'])

    def test_manifest_can_reference_an_already_indexed_small_paper(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve()
            paper=root/'origin_theory.tex'
            paper.write_text('\\section{Closure}\nA source shared by a paper and a contract.\n')
            atlas=root/'semantic_explorer'; atlas.mkdir()
            (atlas/'orientation.json').write_text(json.dumps({'topics':[]}))
            folder=root/'experiments/theory-contracts/shared-paper'; folder.mkdir(parents=True)
            (folder/'contract_index.json').write_text(json.dumps({
                'contract':'shared-paper','question':'Reuses an original paper',
                'source_manifest':'source_manifest.json'
            }))
            (folder/'source_manifest.json').write_text(json.dumps({'sources':[{
                'path':'origin_theory.tex','sha256':hashlib.sha256(paper.read_bytes()).hexdigest()
            }]}))
            old_root,old_out=bc.ROOT,bc.OUT
            bc.ROOT=root; bc.OUT=atlas/'data/corpus.json'
            try:
                corpus=bc.main()
            finally:
                bc.ROOT,bc.OUT=old_root,old_out
            shared=[s for s in corpus['sources'] if s['path']=='origin_theory.tex']
            self.assertEqual(len(shared),1)
            self.assertEqual(shared[0]['kind'],'tex')
            records=[n for n in corpus['nodes'] if n.get('attrs',{}).get('source_id')==shared[0]['id']]
            self.assertTrue(any(n['type']=='paper_source' for n in records))
            self.assertTrue(any(n['type']=='experiment_record' for n in records))

    def test_corpus_builder_merges_source_inventory_without_legacy_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); (root/'verification').mkdir()
            (root/'verification/status_ledger.csv').write_text(
                'claim_id,claim,status,canonical_status,active\nCLAIM.ONE.01,First claim,Open,[O],true\n'
            )
            atlas=root/'semantic_explorer'; atlas.mkdir()
            (atlas/'orientation.json').write_text(json.dumps({'topics':[{
                'id':'topic:test','label':'Test','summary':'Question','keywords':[]
            }]}))
            folder=root/'experiments/theory-contracts/bridge'; folder.mkdir(parents=True)
            (folder/'contract_index.json').write_text(json.dumps({
                'contract':'bridge','question':'A scoped test','verdict':'PARTIAL',
                'claims_tested':['CLAIM.ONE.01'],
                'proof':'PROOF.md','checker':'checker.py',
                'leaves_open':['The local derivation remains open.']
            }))
            (folder/'PROOF.md').write_text('# Proof\nScoped result.\n')
            (folder/'checker.py').write_text('def check():\n    return True\n')
            dep=root/'tfpt_explorer/neutral_source_response.py'; dep.parent.mkdir(parents=True)
            dep.write_text('def respond():\n    return "yes"\n')
            dep_hash=hashlib.sha256(dep.read_bytes()).hexdigest()
            (folder/'source_manifest.json').write_text(json.dumps({'sources':[{
                'path':'tfpt_explorer/neutral_source_response.py','sha256':dep_hash
            }]}))
            index=json.loads((folder/'contract_index.json').read_text())
            index['source_manifest']='source_manifest.json'
            (folder/'contract_index.json').write_text(json.dumps(index))
            old_root,old_out=bc.ROOT,bc.OUT
            bc.ROOT=root; bc.OUT=atlas/'data/corpus.json'
            try:
                corpus=bc.main()
            finally:
                bc.ROOT,bc.OUT=old_root,old_out
            contract=next(n for n in corpus['nodes'] if n['id']=='contract:bridge')
            self.assertEqual(contract['type'],'attempt')
            self.assertEqual(contract['attrs']['leaves_open'],['The local derivation remains open.'])
            self.assertTrue(any(e['src']=='contract:bridge' and e['dst']=='claim:CLAIM.ONE.01' and e['type']=='tests' for e in corpus['edges']))
            registered={s['path'] for s in corpus['sources']}
            self.assertIn('experiments/theory-contracts/bridge/PROOF.md',registered)
            self.assertIn('experiments/theory-contracts/bridge/checker.py',registered)
            self.assertIn('experiments/theory-contracts/bridge/source_manifest.json',registered)
            self.assertIn('tfpt_explorer/neutral_source_response.py',registered)
            node_ids=[n['id'] for n in corpus['nodes']]
            self.assertEqual(len(node_ids),len(set(node_ids)))
            self.assertFalse((root/'verification/theory_graph.json').exists())


if __name__=='__main__': unittest.main()
