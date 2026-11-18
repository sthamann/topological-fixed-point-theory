import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest
import extract

FIELDS=['claim_id','claim','status','location','dependencies','script','external_data','supersedes','canonical_status','active']

class ExtractionTests(unittest.TestCase):
    def test_frozen_static_candidate_inventory_has_exact_source_provenance(self):
        path=extract.HERE/'data/static-candidate-resolutions.json'
        data=json.loads(path.read_text())
        self.assertEqual((data['baseline_group_count'],data['case_count'],data['claim_count'],data['script_count']),
                         (89,195,69,67))
        self.assertEqual(data['unique_source_expression_count'],194)
        self.assertEqual(len({c['id'] for c in data['cases']}),195)
        corpus=extract.Corpus(extract.REPO)
        for case in data['cases']:
            self.assertIn(case['interpretation_status'],{'source_origin_resolved_role_scoped','source_reviewed_role_and_scope'})
            source=extract.REPO/case['source']['path']
            if source.is_file() and extract.sha(source.read_bytes())==case['source']['sha256']:
                self.assertEqual(case['source']['line_text'].strip(),source.read_text().splitlines()[case['source']['line']-1].strip())
            else:
                # Frozen historical/stale bytes are not a current mapping.
                packet=corpus.packet(corpus.by_id[case['claim_id']]) if case['claim_id'] in corpus.by_id else {'sources':[],'scripts':[]}
                record={'semantic':{}}
                extract.apply_static_candidate_origins(record,case['claim_id'],packet,corpus)
                accepted=[item for norm in record.get('normalizations',[]) for item in norm.get('cases',[])]
                self.assertNotIn(case['id'],{item['id'] for item in accepted})
            self.assertIn(case['source_origin'],{'function_default_in_source','literal_in_source','runtime_or_external_read','self_source_file_read','expression_in_source','helper_expression_in_source'})

    def test_static_candidate_origins_attach_only_when_source_and_expression_match(self):
        corpus=extract.Corpus(extract.REPO)
        data=json.loads((extract.HERE/'data/static-candidate-resolutions.json').read_text())
        case=data['cases'][0]
        row=corpus.by_id[case['claim_id']]
        packet=corpus.packet(row)
        record={'semantic':{},'source_refs':[],'source_records':[]}
        count=extract.apply_static_candidate_origins(record,case['claim_id'],packet,corpus)
        self.assertGreaterEqual(count,1)
        norm=next(n for n in record['normalizations'] if n['kind']=='static_candidate_source_origins')
        mapped=next(c for c in norm['cases'] if c['candidate_id']==case['candidate_id'])
        sources={s['id']:s for s in record['source_records']}
        self.assertEqual(sources[mapped['source']['source_ref']]['text'],case['source']['line_text'])
        self.assertEqual(sources[mapped['source']['context_ref']]['sha256'],case['source']['sha256'])

    def test_static_candidate_mapping_preserves_legacy_source_current_reviews(self):
        corpus=extract.Corpus(extract.REPO)
        cid='DHR.MODULAR.01'
        packet=corpus.packet(corpus.by_id[cid])
        record={'semantic':{'candidate_review':[
                    {'script_id':'script:v242_spin_statistics_modular','lines':[38],
                     'classification':'technical','reason':'N=16 follows from |G|=4x4.',
                     'input_names':['N']},
                    {'script_id':'script:v242_spin_statistics_modular','lines':[78],
                     'classification':'unknown','reason':'unresolved before source review.',
                     'input_names':['Cexp[IDX[a], IDX[neg(a)]]']}],
                    'unknowns':[
                    'Static candidate origin unresolved: script:v242_spin_statistics_modular:38:N',
                    'Static candidate origin unresolved: script:v242_spin_statistics_modular:78:Cexp[IDX[a], IDX[neg(a)]]']},
                'source_refs':[],'source_records':[]}
        count=extract.apply_static_candidate_origins(record,cid,packet,corpus)
        self.assertEqual(count,1)
        groups=record['semantic']['candidate_review']
        n=next(g for g in groups if '38:N' in g.get('candidate_ids',[]))
        self.assertEqual(n['classification'],'technical')
        self.assertEqual(n['reason'],'N=16 follows from |G|=4x4.')
        self.assertNotIn('Static candidate origin unresolved: script:v242_spin_statistics_modular:38:N',
                         record['semantic']['unknowns'])
        self.assertFalse(any('78:Cexp[IDX[a], IDX[neg(a)]]' in u for u in record['semantic']['unknowns']))

    def test_source_reviewed_supplement_maps_all_33_exact_candidates(self):
        corpus=extract.Corpus(extract.REPO)
        supplement=json.loads((extract.HERE/'data/static-candidate-supplement.json').read_text())
        self.assertEqual(supplement['case_count'],33)
        self.assertEqual(supplement['classification_counts'],
                         {'technical':11,'intermediate':14,'load_bearing':8})
        total=0;historical=0
        supplement_ids={c['id'] for c in supplement['cases']}
        for claim_id in sorted({c['claim_id'] for c in supplement['cases']}):
            if claim_id not in corpus.by_id:
                self.assertTrue(claim_id.startswith(('PRIME.','RH.')))
                historical+=sum(c['claim_id']==claim_id for c in supplement['cases'])
                self.assertEqual(extract.apply_static_candidate_origins({'semantic':{}},claim_id,{'sources':[],'scripts':[]},corpus),0)
                continue
            packet=corpus.packet(corpus.by_id[claim_id])
            record={'semantic':{},'source_refs':[],'source_records':[]}
            extract.apply_static_candidate_origins(record,claim_id,packet,corpus)
            norm=next(n for n in record['normalizations'] if n['kind']=='static_candidate_source_origins')
            cases=[c for c in norm['cases'] if c.get('id') in supplement_ids]
            total+=len(cases)
            self.assertTrue(all(c.get('physical_derivation_claimed') is False for c in cases))
            self.assertTrue(all(c.get('purpose') and c.get('limits') for c in cases))
            self.assertTrue(all(c['source']['source_ref'] in record['source_refs'] for c in cases))
        self.assertEqual(total+historical,33)
        self.assertGreater(total,0)

    def test_stale_static_origin_mapping_is_removed_and_candidate_returns_to_unknown(self):
        corpus=extract.Corpus(extract.REPO)
        data=json.loads((extract.HERE/'data/static-candidate-resolutions.json').read_text())
        case=next(c for c in data['cases'] if c['interpretation_status']=='source_reviewed_role_and_scope')
        packet=corpus.packet(corpus.by_id[case['claim_id']])
        changed=copy.deepcopy(case); changed['source']={**case['source'],'sha256':'0'*64}
        record={'semantic':{'candidate_review':[{'script_id':case['script_id'],'candidate_ids':[case['candidate_id']],
                   'lines':[case['source']['line']],'classification':'technical','reason':'prior mapping','input_names':[case['candidate_name']]}],
                   'unknowns':[]},'source_refs':[],'source_records':[],
                 'normalizations':[{'kind':'static_candidate_source_origins','claim_id':case['claim_id'],'cases':[case]}]}
        with tempfile.TemporaryDirectory() as td:
            resolution=Path(td)/'resolutions.json'
            resolution.write_text(json.dumps({'cases':[changed]}))
            self.assertEqual(extract.apply_static_candidate_origins(record,case['claim_id'],packet,corpus,resolution),0)
        self.assertFalse(any(n.get('kind')=='static_candidate_source_origins' for n in record['normalizations']))
        self.assertEqual(record['semantic']['candidate_review'][0]['classification'],'unknown')
        self.assertTrue(record['semantic']['unknowns'])

    def test_targeted_run_preserves_unselected_source_current_annotation_without_revalidation(self):
        record={'semantic':{'claim_id':'C.01','relationships':[{'type':'uses_claim_under_scope','target_id':'OLD.ID'}]},'input_sha256':'old-prompt'}
        self.assertTrue(extract.preserve_unselected_cached_record('C.01',record,{'TARGET.01'}))
        self.assertFalse(extract.preserve_unselected_cached_record('TARGET.01',record,{'TARGET.01'}))
        self.assertFalse(extract.preserve_unselected_cached_record('C.01',{'error':'previous error'},{'TARGET.01'}))

    def test_self_relationship_normalization_preserves_audit_record_and_other_edges(self):
        self_edge={'target_claim_id':'C.01','relation':'mentions','source_refs':['source:a']}
        other={'target_claim_id':'D.01','relation':'uses_result','source_refs':['source:b']}
        rec={'semantic':{'relationships':[self_edge,other]}}
        self.assertEqual(extract.remove_self_relationships('C.01',rec),1)
        self.assertEqual(rec['semantic']['relationships'],[other])
        self.assertEqual(rec['normalizations'][0]['removed'],[self_edge])
        self.assertEqual(extract.remove_self_relationships('C.01',rec),0)

    def test_self_claim_origin_is_cleared_as_local_ownership_not_upstream_evidence(self):
        rec={'semantic':{'premises':[{'id':'P1','origin_claim_id':'C.01','text':'Ledger-owned premise','source_refs':['source:x']}],
                         'input_audit':[{'name':'a','upstream_claim_id':'C.01','value':'1','source_refs':['source:y']}],
                         'derivations':[{'id':'D1','upstream_claim_ids':['C.01','UPSTREAM.01']}]}}
        self.assertEqual(extract.remove_self_claim_provenance('C.01',rec),3)
        self.assertEqual(rec['semantic']['premises'][0]['origin_claim_id'],'')
        self.assertEqual(rec['semantic']['input_audit'][0]['upstream_claim_id'],'')
        self.assertEqual(rec['semantic']['derivations'][0]['upstream_claim_ids'],['UPSTREAM.01'])
        self.assertEqual(rec['semantic']['premises'][0]['source_refs'],['source:x'])
        self.assertEqual(rec['semantic']['input_audit'][0]['source_refs'],['source:y'])
        normalization=rec['normalizations'][-1]
        self.assertEqual(normalization['kind'],'remove_self_claim_provenance')
        self.assertEqual(len(normalization['items']),3)

    def test_relationship_without_target_is_preserved_as_local_source_context(self):
        rel={'target_claim_id':'','relation':'mentions','reason':'A nearby claim is discussed.',
             'conditions':['textual only'],'source_refs':['source:x']}
        rec={'semantic':{'relationships':[rel],'objects':[]}}
        self.assertEqual(extract.preserve_unresolved_relationship_notes('C.01',rec),1)
        self.assertEqual(rec['semantic']['relationships'],[])
        note=rec['semantic']['objects'][0]
        self.assertEqual(note['kind'],'source_context')
        self.assertIn('A nearby claim is discussed.',note['description'])
        self.assertEqual(note['source_refs'],['source:x'])
        self.assertEqual(rec['normalizations'][-1]['items'],[rel])

    def test_source_registry_merge_preserves_historical_exact_ids_and_current_wins(self):
        current={'id':'source:current','path':'a.py','sha256':'new','line_start':1,'line_end':1}
        oldsame={'id':'source:current','path':'old.py','sha256':'old','line_start':2,'line_end':2}
        oldmissing={'id':'source:historical','path':'b.py','sha256':'hash','line_start':4,'line_end':8}
        graph={'sources':[current]}
        cache={'C.01':{'source_records':[oldsame,{'id':'source:historical','path':'b.py','sha256':'hash','line_start':4,'line_end':8,'text':'old excerpt'},
                                            {'id':'source:cached','path':'c.py','sha256':'c','line_start':1,'line_end':2,'text':'excerpt'}]}}
        extract.merge_source_registry(graph,cache)
        sources={s['id']:s for s in graph['sources']}
        self.assertEqual(sources['source:current'],current)
        self.assertEqual(sources['source:historical']['text'],'old excerpt')
        self.assertIn('text',sources['source:cached'])

    def test_cache_source_hydration_reconstructs_only_hash_current_exact_ranges(self):
        class FakeCorpus:
            def __init__(self, root): self.root=root; self.sources={}
            def read(self,path):
                data=path.read_bytes(); return data.decode(),extract.sha(data)
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); src=root/'paper.txt'; src.write_text('one\ntwo\nthree\n')
            digest=extract.sha(src.read_bytes());
            exact={'id':'source:exact','path':'paper.txt','sha256':digest,'line_start':2,'line_end':3}
            stale={'id':'source:stale','path':'paper.txt','sha256':'0'*64,'line_start':1,'line_end':1}
            corpus=FakeCorpus(root);corpus.sources={'source:exact':exact,'source:stale':stale}
            cache={'C.01':{'source_refs':['source:exact','source:stale']}}
            count=extract.hydrate_cache_source_records(cache,corpus)
        self.assertEqual(count,2)
        records={x['id']:x for x in cache['C.01']['source_records']}
        self.assertEqual(records['source:exact']['text'],'two\nthree')
        self.assertNotIn('text',records['source:stale'])

    def test_unregistered_upstream_becomes_source_scoped_premise_or_is_removed(self):
        cited={'id':'source:x','text':'Contract explicitly names LEGACY.RESULT.01 and no active result.'}
        uncited={'id':'source:y','text':'Only arithmetic is described.'}
        record={'source_records':[cited,uncited],'semantic':{
            'premises':[{'id':'P1','text':'A','kind':'input','origin_claim_id':'','source_refs':['source:x']}],
            'derivations':[{'id':'D1','premise_ids':['P1'],'upstream_claim_ids':['LEGACY.RESULT.01','HALLUCINATED.01'],
                'result_ids':['R1'],'argument':'Uses source and local premise','scope':'conditional','source_refs':['source:x','source:y']}],
            'results':[{'id':'R1'}],'relationships':[
                {'target_claim_id':'LEGACY.MENTION.01','relation':'mentions','reason':'Source mention','conditions':[], 'source_refs':['source:x']},
                {'target_claim_id':'HALLUCINATED.02','relation':'mentions','reason':'Unsupported target','conditions':[], 'source_refs':['source:y']}],
            'input_audit':[]}}
        count=extract.normalize_unregistered_claim_references('C.01',record,{'C.01'})
        sem=record['semantic'];d=sem['derivations'][0]
        self.assertEqual(count,4)
        self.assertEqual(d['upstream_claim_ids'],[])
        self.assertIn('PCTX1',d['premise_ids'])
        self.assertEqual(sem['premises'][-1]['kind'],'unregistered_reference')
        self.assertIn('not imported',sem['premises'][-1]['text'])
        self.assertEqual(sem['relationships'],[])
        self.assertIn('source_context',[o['kind'] for o in sem['objects']])
        self.assertTrue(any('not named' in d['scope'] for d in sem['derivations']))

    def test_semantic_source_freshness_uses_persisted_citations_not_packet_subset(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); f=root/'source.tex';f.write_text('evidence\n')
            digest=extract.sha(f.read_bytes())
            record={'semantic':{'results':[{'source_refs':['source:full']} ]},
                'source_records':[{'id':'source:full','path':'source.tex','sha256':digest,'line_start':1,'line_end':1,'text':'evidence'}]}
            fresh=extract.semantic_source_freshness(record,root)
            self.assertEqual(fresh['status'],'verified_current')
            f.write_text('changed\n')
            stale=extract.semantic_source_freshness(record,root)
            self.assertEqual(stale['status'],'stale')

    def test_static_candidate_classification_separates_derived_checks_and_selected_literals(self):
        derived={'kind':'numeric_assignment','value':'Nstar**2 * c3 / 24','name':'A_s'}
        check={'kind':'numeric_assignment','value':'all(x == y for x, y in pairs)','name':'ok'}
        literal={'kind':'numeric_assignment','value':'0.0073','name':'alpha0'}
        external={'kind':'external_or_state_input','value':'os.environ["MODEL"]','name':'model'}
        self.assertEqual(extract.classify_static_candidate(derived)[0],'intermediate')
        self.assertEqual(extract.classify_static_candidate(check)[0],'technical')
        self.assertEqual(extract.classify_static_candidate(literal)[0],'unknown')
        self.assertEqual(extract.classify_static_candidate(external)[0],'unknown')

    def test_static_candidate_coverage_is_exact_and_preserves_unknown_origins(self):
        data={'candidate_review':[],'unknowns':[]}
        packet={'scripts':[{'script_id':'script:v1','input_candidates':[
            {'candidate_id':'script:v1:L10:x','line':10,'name':'x','kind':'numeric_assignment','value':'a+b'},
            {'candidate_id':'script:v1:L11:seed','line':11,'name':'seed','kind':'numeric_assignment','value':'0.1'}]}]}
        additions=extract.mark_omitted_candidates_unknown(data,packet)
        self.assertEqual({x['classification'] for x in additions},{'intermediate','unknown'})
        self.assertEqual(sum(len(x['candidate_ids']) for x in data['candidate_review']),2)
        self.assertTrue(any('origin unresolved' in x for x in data['unknowns']))
        self.assertEqual(extract.mark_omitted_candidates_unknown(data,packet),[])

    def test_input_origin_correction_is_exact_and_audited(self):
        refs=['source:ledger']
        rec={'semantic':{'input_audit':[{'name':'status','origin':'declared_axiom','selection_status':'assumed_or_selected','impact':'conditional','explanation':'old','source_refs':refs}]}}
        self.assertTrue(extract.correct_input_audit_type(rec,'FLAV.CKM.01','status','declared_axiom','external_data',selection_status='not_applicable',impact='benign',reason='Ledger status is metadata, not a physical input.',source_refs=refs))
        item=rec['semantic']['input_audit'][0]
        self.assertEqual((item['origin'],item['selection_status'],item['impact']),('external_data','not_applicable','benign'))
        self.assertEqual(rec['normalizations'][0]['before']['origin'],'declared_axiom')
        self.assertFalse(extract.correct_input_audit_type(rec,'FLAV.CKM.01','status','declared_axiom','external_data',reason='wrong source',source_refs=['source:wrong']))

    def test_incomplete_derivation_retyping_preserves_text_without_derives_edge(self):
        d={'id':'D1','label':'Open route','premise_ids':['P1'],'upstream_claim_ids':[],'result_ids':[],
           'argument':'The source marks the route open.','scope':'open only','source_refs':['source:a']}
        rec={'semantic':{'derivations':[d],'proof_gaps':[]}}
        self.assertTrue(extract.reclassify_unfinished_derivation('C.01',rec,'D1','open_question'))
        self.assertFalse(rec['semantic']['derivations'])
        self.assertEqual(rec['semantic']['proof_gaps'][0]['source_refs'],['source:a'])
        self.assertIn('The source marks the route open.',rec['normalizations'][0]['original']['argument'])
        self.assertFalse(extract.reclassify_unfinished_derivation('C.01',rec,'D1','open_question'))

    def test_packet_schema_enumerates_only_supplied_external_ids(self):
        packet={'claim_id':'C.01','sources':[{'id':'source:one','path':'a.py','text':'Related claim D.X.02 is cited.'}],
                'scripts':[{'script_id':'script:a'}],'lean_candidates':[{'theorem_id':'A.b'}],
                'related_claims':[{'claim_id':'B.03'}],'topics':{'topic:origin':{}}}
        p=extract.schema_for_packet(packet)['properties']
        self.assertEqual(p['relationships']['items']['properties']['target_claim_id']['enum'],['B.03','C.01','D.X.02'])
        self.assertEqual(p['input_audit']['items']['properties']['upstream_claim_id']['enum'],['','B.03','C.01','D.X.02'])
        self.assertEqual(p['script_analysis']['items']['properties']['script_id']['enum'],['script:a'])
        self.assertEqual(p['lean_links']['items']['properties']['theorem_id']['enum'],['A.b'])
        self.assertEqual(p['relationships']['items']['properties']['source_refs']['items']['enum'],['source:one'])

    def test_cache_freshness_distinguishes_changed_inputs_from_stale_sources(self):
        rec={'semantic':{'title_de':'x'},'input_sha256':'same'}
        self.assertEqual(extract.cache_freshness(rec,'same',True)['state'],'current')
        self.assertEqual(extract.cache_freshness(rec,'new',True)['state'],'source_current_prompt_or_packet_changed')
        self.assertEqual(extract.cache_freshness(rec,'same',False)['state'],'source_unverified')
        self.assertEqual(extract.cache_freshness({},'new',False)['state'],'missing_annotation')

    def fixture(self,root):
        (root/'verification').mkdir()
        rows=[dict(zip(FIELDS,['TEST.ONE.01','under A and B: D','Conditional','','','v1_test.py','-','','Conditional','true']))]
        with (root/'verification/status_ledger.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
        (root/'verification/v1_test.py').write_text('import os\neta = 1/4\nmodel = os.environ["MODEL"]\ndef run(rate=0.1):\n    x = 0\n    x += 2\n')
        return extract.Corpus(root)

    def test_hidden_input_forms_and_original_preservation(self):
        with tempfile.TemporaryDirectory() as d:
            c=self.fixture(Path(d));p=c.packet(c.rows[0]);inputs=p['scripts'][0]['input_candidates']
            by_name={x['name']:x for x in inputs}
            self.assertEqual(by_name['eta']['value'],'1/4')
            self.assertEqual(by_name['model']['kind'],'external_or_state_input')
            self.assertEqual(by_name['run.rate']['kind'],'function_default')
            self.assertTrue(any(x['line']==6 for x in inputs))
            g=extract.build_graph(c,{},'fixture');claim=next(n for n in g['nodes'] if n['id']=='claim:TEST.ONE.01')
            self.assertEqual(claim['attrs']['ledger'],c.rows[0]);self.assertEqual(claim['attrs']['enrichment_state'],'pending')

    def test_corpus_read_and_static_script_inventory_refresh_after_edit(self):
        with tempfile.TemporaryDirectory() as d:
            c=self.fixture(Path(d));p=Path(d)/'verification/v1_test.py'
            sid,first,first_candidates,_=c.script(p)
            p.write_text(p.read_text().replace('eta = 1/4','eta = 1/5'))
            sid2,second,second_candidates,_=c.script(p)
            self.assertEqual(sid,sid2)
            self.assertNotEqual(first['sha256'],second['sha256'])
            self.assertIn('eta = 1/5',second['text'])
            self.assertEqual(next(c for c in first_candidates if c['name']=='eta')['value'],'1/4')
            self.assertEqual(next(c for c in second_candidates if c['name']=='eta')['value'],'1/5')

    def test_cached_embeddings_hydrate_checkpoints_by_exact_text_digest_and_contract(self):
        with tempfile.TemporaryDirectory() as d:
            cache_path=Path(d)/'embeddings-cache.json'
            exact={'id':'claim:A','type':'claim','label':'Normalisierte Materieantwort',
                'summary':'Flächennormierung','attrs':{}}
            changed={'id':'claim:B','type':'claim','label':'Normalisierte Materieantwort',
                'summary':'Andere Flächennormierung','attrs':{}}
            excluded={'id':'aspect:X','type':'theory_aspect','label':'normalisierte Materieantwort',
                'summary':'Flächennormierung','attrs':{}}
            digest=extract.sha((exact['label']+'\n'+exact['summary']).encode())
            vector=[0.125]*extract.EMBEDDING_DIMENSIONS
            cache_path.write_text(json.dumps({digest:vector}))  # legacy fixed-model cache
            graph={'meta':{},'nodes':[exact,changed,excluded]}
            extract.hydrate_cached_embeddings(graph,cache_path)
            self.assertEqual(exact['attrs']['embedding'],vector)
            self.assertEqual(exact['attrs']['embedding_model'],extract.EMBEDDING_MODEL)
            self.assertNotIn('embedding',changed['attrs'])
            self.assertNotIn('embedding',excluded['attrs'])
            self.assertEqual(graph['meta']['embedding_coverage'],{
                'eligible_nodes':2,'embedded_nodes':1,'cache_hits':1,'missing_nodes':1,
                'invalid_cached_vectors':0,'cache_contract_valid':True})

            cache_path.write_text(json.dumps({extract.EMBEDDING_CACHE_META_KEY:{
                'model':'wrong-model','dimensions':extract.EMBEDDING_DIMENSIONS},digest:vector}))
            graph={'meta':{},'nodes':[{'id':'claim:A','type':'claim','label':exact['label'],
                'summary':exact['summary'],'attrs':{}}]}
            extract.hydrate_cached_embeddings(graph,cache_path)
            self.assertNotIn('embedding',graph['nodes'][0]['attrs'])
            self.assertFalse(graph['meta']['embedding_coverage']['cache_contract_valid'])

    def test_ckm_claim_location_adds_narrow_formula_windows_to_primary_packet(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);c=self.fixture(root)
            row=c.rows[0]
            row.update({'claim_id':'FLAV.CKM.01','claim':'CKM core s12=lambda_C; s23; s13; delta=pi/3+3lC^2',
                'location':'tfpt_1;tfpt_2','script':'-'})
            first=root/'tfpt_1_architecture_e8.tex'
            second=root/'tfpt_2_standard_model.tex'
            first.write_text('before\nCKM CP phase & $\\delta_{\\rm CKM}=\\tfrac\\pi3+3\\lC^2$\nafter\n')
            second.write_text('before\n$s_{12}^{\\mathrm{CKM}}$ & $\\lamC$\n$s_{23}^{\\mathrm{CKM}}$ & formula\n')
            c.docs['tfpt_1']=first;c.docs['tfpt_2']=second
            packet=c.packet(row)
            tex_sources=[s for s in packet['sources'] if s['path'].endswith('.tex')]
            self.assertEqual({s['path'] for s in tex_sources},{'tfpt_1_architecture_e8.tex','tfpt_2_standard_model.tex'})
            self.assertTrue(any('delta_{\\rm CKM}' in s['text'] for s in tex_sources))
            self.assertTrue(any('s_{12}^{\\mathrm{CKM}}' in s['text'] for s in tex_sources))
            for source in tex_sources:
                self.assertIn(source['id'],c.sources)
                self.assertLessEqual(source['line_end']-source['line_start'],22)

    def test_reject_missing_evidence_and_silently_omitted_inputs(self):
        data={k:[] for k in ['premises','results','derivations','script_analysis','lean_links','input_audit','proof_gaps','candidate_review']}
        packet={'sources':[{'id':'s'}],'scripts':[{'script_id':'script:v1','input_candidates':[{'line':4}]}],'lean_candidates':[]}
        with self.assertRaisesRegex(ValueError,'unreviewed'):extract.validate(data,packet)
        data['candidate_review']=[{'script_id':'script:v1','lines':[4]}]
        data['results']=[{'id':'R1','source_refs':[]}]
        with self.assertRaisesRegex(ValueError,'missing evidence'):extract.validate(data,packet)
        data['results'][0]['source_refs']=['invented']
        with self.assertRaisesRegex(ValueError,'unknown source'):extract.validate(data,packet)

    def test_all_premises_required_and_reference_is_not_proof(self):
        with tempfile.TemporaryDirectory() as d:
            c=self.fixture(Path(d));p=c.packet(c.rows[0]);s=p['sources'][0]['id']
            sem={'premises':[{'id':'P1','text':'A','origin_claim_id':''},{'id':'P2','text':'B','origin_claim_id':''}],
                 'results':[{'id':'R1','label':'D','statement':'D'}],
                 'derivations':[{'id':'D1','label':'A and B => D','argument':'conjunction','premise_ids':['P1','P2'],'upstream_claim_ids':[],'result_ids':['R1'],'scope':'conditional'}]}
            g=extract.build_graph(c,{'TEST.ONE.01':{'semantic':sem}},'fixture')
            dnode=next(n for n in g['nodes'] if n['type']=='derivation')
            self.assertTrue(dnode['attrs']['all_premises_required'])
            self.assertEqual(sum(e['type']=='required_together' for e in g['edges']),2)
            self.assertFalse(any(e['type']=='proves' for e in g['edges']))

    def test_unused_parameter_and_candidate_identity(self):
        with tempfile.TemporaryDirectory() as d:
            c=self.fixture(Path(d));p=Path(d)/'verification/v2_unused.py'
            p.write_text('def root(c3, Nfac=8):\n    return 1/Nfac\n')
            sid,_,_,_=c.script(p)
            self.assertEqual(c.script_nodes[sid]['attrs']['static_inputs']['callable_audit'][0]['unused_parameters'],['c3'])
        script={'script_id':'s','input_candidates':[{'line':3,'name':'a, b','candidate_id':'3:a, b'}]}
        groups=[{'script_id':'s','lines':[3],'input_names':['a']}]
        self.assertNotIn('3:a, b',extract.reviewed_candidate_ids(groups,script))
        groups[0]['input_names'].append('b')
        self.assertIn('3:a, b',extract.reviewed_candidate_ids(groups,script))

    def test_self_dependency_and_result_cycle_rejected(self):
        data={k:[] for k in ['premises','results','derivations','script_analysis','lean_links','input_audit','proof_gaps','candidate_review']}
        packet={'claim_id':'T.ONE.01','sources':[{'id':'s'}],'scripts':[],'lean_candidates':[]}
        data['results']=[{'id':'R1','source_refs':['s']}]
        data['derivations']=[{'id':'D1','premise_ids':['R1'],'result_ids':['R1'],'upstream_claim_ids':['T.ONE.01'],'source_refs':['s']}]
        with self.assertRaisesRegex(ValueError,'self claim|cyclic'):extract.validate(data,packet)

    def test_key_file_never_copied_into_graph(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'key.env';p.write_text('OPENAI_API_KEY="sk-test-placeholder" # note\n')
            self.assertEqual(extract.resolve_key(p),'sk-test-placeholder')

    def test_omitted_candidates_are_explicitly_unknown(self):
        packet={'scripts':[{'script_id':'script:v1','input_candidates':[{'line':4},{'line':5}]}]}
        data={'candidate_review':[{'script_id':'script:v1','lines':[4]}],'unknowns':[]}
        extract.mark_omitted_candidates_unknown(data,packet)
        self.assertEqual(data['candidate_review'][-1]['lines'],[5])
        self.assertEqual(data['candidate_review'][-1]['classification'],'unknown')
        self.assertIn('not inferred',data['candidate_review'][-1]['reason'])

    def test_alias_normalization_requires_exact_source_match(self):
        packet={'scripts':[{'script_id':'script:v1'}],'sources':[{'path':'verification/helper.py'}],'lean_candidates':[]}
        data={'script_analysis':[{'script_id':'v1.py'},{'script_id':'invented.py'}],'candidate_review':[],'lean_links':[]}
        extract.normalize_ids(data,packet)
        self.assertEqual(data['script_analysis'][0]['script_id'],'script:v1')
        self.assertEqual(data['script_analysis'][1]['script_id'],'invented.py')

    def test_benign_selected_status_normalizes_only_for_exact_technical_candidate(self):
        packet={'claim_id':'TEST.ONE.01','sources':[{'id':'s'}],
            'scripts':[{'script_id':'script:precision','source_id':'s','input_candidates':[{'candidate_id':'10:mp.mp.dps','line':10,'name':'mp.mp.dps'}]}],
            'lean_candidates':[],'topics':{}}
        data={k:[] for k in ['premises','results','derivations','script_analysis','lean_links','input_audit','proof_gaps','candidate_review']}
        data['candidate_review']=[{'script_id':'script:precision','candidate_ids':['10:mp.mp.dps'],'lines':[10],'classification':'technical','reason':'numeric precision only','input_names':['mp.mp.dps']}]
        data['input_audit']=[{'name':'mp.mp.dps','selection_status':'assumed_or_selected','impact':'benign','source_refs':['s'],'explanation':'Code assigns precision.'}]
        extract.normalize_technical_input_audits(data,packet)
        self.assertEqual(data['input_audit'][0]['selection_status'],'not_applicable')
        extract.validate(data,packet)

        data['candidate_review'][0]['classification']='load_bearing'
        data['input_audit'][0]['selection_status']='assumed_or_selected'
        extract.normalize_technical_input_audits(data,packet)
        self.assertEqual(data['input_audit'][0]['selection_status'],'assumed_or_selected')
        with self.assertRaisesRegex(ValueError,'selected load-bearing input cannot be benign'):
            extract.validate(data,packet)

    def test_duplicate_input_names_are_scoped_to_script_source_and_candidate_identity(self):
        packet={'scripts':[
            {'script_id':'script:technical','source_id':'source:technical','input_candidates':[{'candidate_id':'10:x','line':10,'name':'x'}]},
            {'script_id':'script:physical','source_id':'source:physical','input_candidates':[{'candidate_id':'22:x','line':22,'name':'x'}]},
        ]}
        data={'candidate_review':[
            {'script_id':'script:technical','candidate_ids':['10:x'],'classification':'technical'},
            {'script_id':'script:physical','candidate_ids':['22:x'],'classification':'load_bearing'},
        ],'input_audit':[
            {'name':'x','source_refs':['source:technical'],'selection_status':'assumed_or_selected','impact':'benign','explanation':''},
            {'name':'x','source_refs':['source:physical'],'selection_status':'assumed_or_selected','impact':'benign','explanation':''},
        ]}
        extract.normalize_technical_input_audits(data,packet)
        self.assertEqual(data['input_audit'][0]['selection_status'],'not_applicable')
        self.assertEqual(data['input_audit'][1]['selection_status'],'assumed_or_selected')

    def test_candidate_with_conflicting_classifications_is_not_normalized(self):
        packet={'scripts':[{'script_id':'script:precision','source_id':'source:precision',
            'input_candidates':[{'candidate_id':'10:x','line':10,'name':'x'}]}]}
        data={'candidate_review':[
            {'script_id':'script:precision','candidate_ids':['10:x'],'classification':'technical'},
            {'script_id':'script:precision','candidate_ids':['10:x'],'classification':'unknown'},
        ],'input_audit':[{'name':'x','source_refs':['source:precision'],
            'selection_status':'assumed_or_selected','impact':'benign','explanation':''}]}
        extract.normalize_technical_input_audits(data,packet)
        self.assertEqual(data['input_audit'][0]['selection_status'],'assumed_or_selected')


    def test_reconciled_run_unknown_requires_current_receipt_and_exact_script(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);script=root/'check.py';script.write_text('print("ok")\n')
            receipt=root/'run.json';receipt.write_text(json.dumps({'command':['python3','check.py'],'exit_code':0,'script_sha256':extract.sha(script.read_bytes())}))
            record={'semantic':{'unknowns':['Ein Ausführungsresultat des Skripts ist nicht geliefert.','Physical selection remains open.']},
                'source_records':[{'id':'run','path':'run.json','sha256':extract.sha(receipt.read_bytes())}],
                'normalizations':[{'kind':'run_evidence_gap_closed','claim_id':'C.01','source_refs':['run']}]}
            stale=copy.deepcopy(record);script.write_text('different source')
            self.assertEqual(extract.reconcile_run_evidence('C.01',stale,root),0)
            script.write_text('print("ok")\n')
            self.assertEqual(extract.reconcile_run_evidence('C.01',record,root),1)
            self.assertEqual(record['semantic']['unknowns'],['Physical selection remains open.'])
            self.assertEqual(extract.reconcile_run_evidence('C.01',record,root),0)

    def test_text_repair_recovers_encoding_and_preserves_ambiguous_damage(self):
        record={'semantic':{'overview_de':'f\x00fcr \x03bc und \x1b[31mJ\x1b[0m','scope':'x\x10y'}}
        extract.normalize_annotation_text('C.01',record,Path('.'),{})
        self.assertEqual(record['semantic']['overview_de'],'für μ und J')
        self.assertIn('[unresolved encoding U+0010]',record['semantic']['scope'])
        self.assertIn('no missing mathematical symbol was guessed',record['semantic']['unknowns'][0])
        self.assertEqual(extract.annotation_text_issues(record['semantic']),[])
        self.assertTrue(extract.annotation_requires_repair(record))
        self.assertTrue(extract.annotation_text_issues({'text':'bad\x00data'}))
        self.assertEqual(extract.normalize_annotation_text('C.01',record,Path('.'),{}),0)

    def test_external_rh_links_remain_scoped_context_without_fabricated_identity(self):
        moved={'theorem_id':'lean_theorem:rh/lean/RH/CoxeterCompletion.lean:t','relation':'conditional','explanation':'Requires the original operator.','result_ids':['R1'],'source_refs':['s']}
        local={**moved,'theorem_id':'lean_theorem:experiments/local.lean:t'}
        record={'semantic':{'lean_links':[moved,local]}}
        self.assertEqual(extract.preserve_external_lean_context('C.01',record,set()),1)
        self.assertEqual(record['semantic']['lean_links'],[local])
        context=record['semantic']['objects'][0]
        self.assertIn('Requires the original operator.',context['description'])
        self.assertIn('unresolved',context['description'])
        self.assertEqual(record['normalizations'][0]['items'],[moved])

    def test_extra_lean_sources_keep_namespace_roles_and_unchecked_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'verification').mkdir()
            with (root/'verification/status_ledger.csv').open('w',newline='') as stream:csv.DictWriter(stream,fieldnames=FIELDS).writeheader()
            path=root/'experiments/theory-contracts/test/Proof.lean';path.parent.mkdir(parents=True)
            path.write_text('/- theorem fake : False := sorry -/\nnamespace A\nsection local\naxiom chosen : True\ntheorem same : True := by trivial\nend local\nend A\nnamespace B\ntheorem same : True := by sorry\nend B\n')
            run=path.parent/'runs/run1/Probe.lean';run.parent.mkdir(parents=True);run.write_text('theorem probe : True := by trivial\n')
            corpus=extract.Corpus(root);decls=list(corpus.lean_nodes.values())
            self.assertEqual(len(decls),4)
            self.assertEqual({n['label'] for n in decls},{'A.chosen','A.same','B.same','probe'})
            self.assertTrue(all(n['attrs']['kernel_checked'] is False for n in decls))
            self.assertTrue(all(n['attrs']['scope']=='local_private' for n in decls))
            self.assertEqual(next(n for n in decls if n['label']=='probe')['attrs']['source_role'],'generated_run')
            self.assertTrue(next(n for n in decls if n['label']=='B.same')['attrs']['contains_sorry'])
            self.assertTrue(next(n for n in decls if n['label']=='A.chosen')['attrs']['declared_axiom'])
            self.assertTrue(any(source['path'].endswith('Proof.lean') and source['line_start']==1 for source in corpus.sources.values()))

    def test_lean_named_types_instances_and_multiple_attributes_are_inventory_not_proofs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'verification').mkdir()
            with (root/'verification/status_ledger.csv').open('w',newline='') as stream:
                csv.DictWriter(stream,fieldnames=FIELDS).writeheader()
            path=root/'experiments/theory-contracts/test/Types.lean';path.parent.mkdir(parents=True)
            path.write_text('/- structure Fake /- class Nested -/ -/\nnamespace Types\n'
                            'structure Record where\n  n : Nat\n'
                            'class Marker where\n  ok : True\n'
                            'inductive Flag where\n  | selected\n'
                            'instance namedInstance : Inhabited Nat := ⟨0⟩\n'
                            '@[simp] @[local simp] private theorem decorated : True := by sorry\n'
                            'end Types\n')
            corpus=extract.Corpus(root);decls=list(corpus.lean_nodes.values())
            self.assertEqual({n['label']:n['attrs']['declaration_kind']for n in decls},
                             {'Types.Record':'structure','Types.Marker':'class','Types.Flag':'inductive',
                              'Types.namedInstance':'instance','Types.decorated':'theorem'})
            self.assertTrue(all(n['attrs']['kernel_checked']is False for n in decls))
            self.assertTrue(all(n['type']=='lean_declaration'for n in decls if n['attrs']['declaration_kind']!='theorem'))
            self.assertTrue(next(n for n in decls if n['label']=='Types.decorated')['attrs']['contains_sorry'])


    def test_modular_paper_packet_uses_physical_source_owner(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'verification').mkdir();(root/'documents').mkdir()
            row={key:'' for key in FIELDS};row.update(claim_id='TEST.MODULE.01',claim='Conditional source',location='paper',script='v1_actual.py',active='true')
            with (root/'verification/status_ledger.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=FIELDS);writer.writeheader();writer.writerow(row)
            (root/'paper.tex').write_text('\\input{documents/part}\n')
            (root/'documents/part.tex').write_text('\\section{Actual owner}\n\\veri{v1\\_actual.py} requires a selected source.\n')
            corpus=extract.Corpus(root);packet=corpus.packet(corpus.by_id['TEST.MODULE.01'])
            physical=[source for source in packet['sources'] if source['path']=='documents/part.tex']
            self.assertEqual(len(physical),1);self.assertIn('requires a selected source',physical[0]['text'])
            self.assertEqual(packet['coverage']['document_source_errors'],[])
            (root/'documents/part.tex').unlink()
            packet=corpus.packet(corpus.by_id['TEST.MODULE.01'])
            self.assertTrue(packet['coverage']['document_source_errors'])

    def test_source_relocation_requires_complete_unique_unchanged_excerpt(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'verification').mkdir();(root/'documents').mkdir()
            with (root/'verification/status_ledger.csv').open('w',newline='') as stream:csv.DictWriter(stream,fieldnames=FIELDS).writeheader()
            (root/'paper.tex').write_text('\\input{documents/part}\n')
            part=root/'documents/part.tex';part.write_text('header\nExact old premise\n')
            corpus=extract.Corpus(root)
            source={'id':'old','path':'paper.tex','sha256':'historical full-file hash','line_start':100,'line_end':100,'text':'Exact old premise'}
            record={'semantic':{'premises':[{'source_refs':['old']}]},'source_refs':['old'],'source_records':[source]}
            original=copy.deepcopy(record)
            preview=copy.deepcopy(original);preview['source_records'][0]['line_end']=120
            self.assertEqual(extract.refresh_equivalent_source_citations('C.01',preview,corpus),0)
            self.assertEqual(preview['semantic']['premises'][0]['source_refs'],['old'])
            self.assertEqual(extract.refresh_equivalent_source_citations('C.01',record,corpus),1)
            self.assertEqual(record['source_records'][0]['path'],'documents/part.tex')
            self.assertEqual(record['normalizations'][0]['items'][0]['original_source'],source)
            self.assertNotEqual(record['semantic']['premises'][0]['source_refs'],['old'])
            part.write_text('Exact old premise\nExact old premise\n')
            self.assertEqual(extract.refresh_equivalent_source_citations('C.01',copy.deepcopy(original),corpus),0)
            part.write_text('Changed premise\n')
            self.assertEqual(extract.refresh_equivalent_source_citations('C.01',copy.deepcopy(original),corpus),0)


    def test_unregistered_premise_origin_retains_full_premise_and_conjunction(self):
        premise={'id':'P3','text':'LEGACY.01 supplies the conditional homology statement.','kind':'ledger_claim','origin_claim_id':'LEGACY.01','source_refs':['source:premise']}
        record={'source_records':[{'id':'source:premise','text':'LEGACY.01 is named historically.'}],
                'semantic':{'premises':[copy.deepcopy(premise)],'derivations':[{'id':'D1','premise_ids':['P3'],'upstream_claim_ids':[],'result_ids':['R1'],'scope':'Only the original source model.'}]}}
        self.assertEqual(extract.normalize_unregistered_claim_references('CURRENT.01',record,{'CURRENT.01'}),1)
        expected={**premise,'origin_claim_id':''}
        self.assertEqual(record['semantic']['premises'],[expected])
        self.assertEqual(record['semantic']['derivations'][0]['premise_ids'],['P3'])
        self.assertEqual(record['semantic']['derivations'][0]['scope'],'Only the original source model.')
        change=record['normalizations'][0]['items'][0]
        self.assertEqual(change['original_premise'],premise)
        self.assertTrue(change['source_named_in_cited_excerpt'])
        self.assertIn(premise['text'],record['semantic']['objects'][0]['description'])
        self.assertIn('absent from the current ledger',record['semantic']['unknowns'][0])
        self.assertEqual(extract.normalize_unregistered_claim_references('CURRENT.01',record,{'CURRENT.01'}),0)

    def test_uncited_historical_origin_is_distinguished_from_own_premise_evidence(self):
        record={'source_records':[{'id':'own','text':'The actual arithmetic statement.'},{'id':'other','text':'LEGACY.01 occurs only in the ledger.'}],
                'semantic':{'premises':[{'id':'P2','text':'Actual arithmetic remains a premise.','kind':'arithmetic','origin_claim_id':'LEGACY.01','source_refs':['own']}]}}
        extract.normalize_unregistered_claim_references('CURRENT.01',record,{'CURRENT.01'})
        change=record['normalizations'][0]['items'][0]
        self.assertFalse(change['source_named_in_cited_excerpt'])
        self.assertTrue(change['source_named_in_supplied_records'])
        self.assertEqual(record['semantic']['objects'][0]['source_refs'],['own','other'])
        self.assertIn('not named in this premise',record['semantic']['objects'][0]['description'])

    def test_new_origin_validator_and_schema_require_actual_canonical_ledger_id(self):
        packet={'claim_id':'CURRENT.01','sources':[{'id':'s','text':'Historical LEGACY.01 and canonical KNOWN.01.'}],
                'scripts':[],'lean_candidates':[],'topics':{},'related_claims':[{'claim_id':'KNOWN.01'}]}
        data={key:[] for key in ['premises','results','derivations','script_analysis','lean_links','input_audit','proof_gaps','candidate_review','relationships']}
        data['premises']=[{'id':'P1','text':'LEGACY.01 is a historical source label.','kind':'context','origin_claim_id':'LEGACY.01','source_refs':['s']}]
        known={'CURRENT.01','KNOWN.01'}
        with self.assertRaisesRegex(ValueError,'unregistered premise origin'):
            extract.validate(data,packet,known)
        enum=extract.schema_for_packet(packet,known)['properties']['premises']['items']['properties']['origin_claim_id']['enum']
        self.assertNotIn('LEGACY.01',enum);self.assertIn('KNOWN.01',enum)
        data['premises'][0]['origin_claim_id']='';extract.validate(data,packet,known)
        data['premises'][0]['origin_claim_id']='KNOWN.01';extract.validate(data,packet,known)


class ControlFeedbackTests(unittest.TestCase):
    def packet(self):
        return {'claim_id':'TEST.TEXT.01','sources':[{'id':'s','path':'original.py','text':'The original source states mu4; physical selection remains open.'}],
                'scripts':[],'lean_candidates':[],'topics':{},'related_claims':[]}

    def annotation(self,text):
        result={key:[] for key in ['premises','results','derivations','script_analysis','lean_links','input_audit','proof_gaps','candidate_review','relationships']}
        result['overview_de']=text
        return result

    def response(self,data):
        return {'status':'completed','output':[{'content':[{'type':'output_text','text':json.dumps(data)}]}]}

    def test_retry_gets_field_paths_and_original_source_without_sanitizing(self):
        from unittest.mock import patch
        packet=self.packet();invalid=self.annotation('bad\x00notation');valid=self.annotation('mu4; physical selection remains open.')
        calls=[];responses=iter([self.response(invalid),self.response(valid)])
        def api(path,request,key):
            calls.append(copy.deepcopy(request));return next(responses)
        with patch.object(extract,'api',side_effect=api):record=extract.extract_one(packet,'unused','gpt-5.4-mini')
        self.assertEqual(record['semantic'],valid)
        self.assertEqual(len(calls),2)
        self.assertTrue(all(call['instructions']==extract.PROMPT for call in calls))
        feedback=json.loads(calls[1]['input'])
        self.assertEqual(feedback['source_packet'],packet)
        self.assertEqual(feedback['invalid_extraction']['overview_de'],'bad\x00notation')
        self.assertEqual(feedback['annotation_text_issues'],[{'path':'/overview_de','codepoints':['U+0000']}])
        self.assertIn('Do not strip controls',feedback['repair_instruction'])
        self.assertEqual(calls[1]['model'],'gpt-5.6-terra')

    def test_repeated_damaged_output_remains_a_validation_error(self):
        from unittest.mock import patch
        damaged=self.response(self.annotation('bad\x1bnotation'))
        with patch.object(extract,'api',return_value=damaged) as api:
            with self.assertRaisesRegex(ValueError,'annotation contains damaged control characters'):
                extract.extract_one(self.packet(),'unused','gpt-5.6-terra')
        self.assertEqual(api.call_count,2)

if __name__=='__main__':unittest.main()
