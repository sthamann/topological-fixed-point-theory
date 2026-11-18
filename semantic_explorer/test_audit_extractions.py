import copy
import tempfile
import unittest
from pathlib import Path

import audit_extractions as audit
from extract import sha,validate,schema_for_packet


class AuditProposalSourceTests(unittest.TestCase):
    def test_proposal_issue_refs_are_registered_in_quality_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "paper.tex"
            raw = b"intro\n\\delta_{CKM}=...\nconclusion\n"
            source_path.write_bytes(raw)
            digest = sha(raw)
            proposal = {
                "source_hashes": [{"source_id": "source:paper", "path": "paper.tex", "sha256": digest}],
                "sources": [{"source_id": "source:paper", "path": "paper.tex", "sha256": digest,
                             "line_start": 2, "line_end": 2}],
                "review": {"issues": [{"field_path": "semantic.overview_de", "finding": "Check source wording.",
                    "evidence": [{"source_id": "source:paper", "quote": "\\delta_{CKM}=..."}]}]},
            }
            registered = audit.verified_proposal_sources(proposal, root)
            issues, sources = audit.proposal_quality_findings(proposal, registered)
            self.assertEqual(issues[0]["source_refs"], ["source:paper"])
            self.assertEqual(sources["source:paper"]["line_start"], 2)
            self.assertEqual(sources["source:paper"]["line_end"], 2)
            self.assertEqual(sources["source:paper"]["sha256"], digest)

            del proposal["sources"]
            registered = audit.verified_proposal_sources(proposal, root)
            issues, sources = audit.proposal_quality_findings(proposal, registered)
            self.assertEqual(issues[0]["code"], "second_pass_source_unregistered")
            self.assertEqual(sources, {})

    def test_proposal_source_metadata_fails_closed_when_stale_or_out_of_range(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "paper.tex").write_text("one\ntwo\n", encoding="utf-8")
            proposal = {"source_hashes": [{"source_id": "source:x", "path": "paper.tex", "sha256": "0" * 64}],
                "sources": [{"source_id": "source:x", "path": "paper.tex", "sha256": "0" * 64,
                             "line_start": 1, "line_end": 3}]}
            self.assertEqual(audit.verified_proposal_sources(proposal, root), {})

    def test_annotation_refs_validate_against_full_registry_not_packet_subset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "paper.tex"
            source_path.write_text("source text\n", encoding="utf-8")
            record = {"id":"source:paper", "path":"paper.tex", "sha256":sha(source_path.read_bytes()),
                      "line_start":1, "line_end":1}
            issues = audit.annotation_source_issues({"source_refs":["source:paper"]}, {"source:paper":record}, root)
            self.assertEqual(issues, [])
            source_path.write_text("changed text\n", encoding="utf-8")
            stale = audit.annotation_source_issues({"source_refs":["source:paper"]}, {"source:paper":record}, root)
            self.assertEqual(stale[0]["code"], "source_ref_stale")

    def test_quality_source_records_rebuild_even_when_proposal_claim_is_cached(self):
        proposal = {"review":{"issues":[{"field_path":"overview_de","finding":"Finding",
            "evidence":[{"source_id":"source:paper"}]}]}}
        record = {"source:paper":{"id":"source:paper","path":"paper.tex","sha256":"digest",
                                  "line_start":1,"line_end":2}}
        issues, sources = audit.proposal_quality_findings(proposal, record)
        self.assertEqual(len(issues), 1)
        self.assertIn("source:paper", sources)

    def test_stale_self_upstream_review_is_suppressed_only_after_link_is_gone(self):
        finding={'code':'self_upstream','message':'old self link'}
        self.assertFalse(audit.manual_issue_still_applies(finding,{'derivations':[{'upstream_claim_ids':[]} ]},'C.01'))
        self.assertTrue(audit.manual_issue_still_applies(finding,{'derivations':[{'upstream_claim_ids':['C.01']} ]},'C.01'))
        self.assertTrue(audit.manual_issue_still_applies({'code':'annotation_integrity'}, {}, 'C.01'))

    def test_validation_adds_text_from_recorded_annotation_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source_path=root/'paper.tex'
            raw=b'QGEO.SYM.01 is cited here.\n'
            source_path.write_bytes(raw)
            packet={'sources':[]}
            semantic={'relationships':[{'target_claim_id':'QGEO.SYM.01','source_refs':['source:paper']}]}
            registry={'source:paper':{'id':'source:paper','path':'paper.tex','sha256':sha(raw),
                                      'line_start':1,'line_end':1}}
            sources=audit.validation_source_records(packet,semantic,registry,root)
            self.assertEqual(sources[0]['id'],'source:paper')
            self.assertIn('QGEO.SYM.01',sources[0]['text'])


class AuditLeanRegistryTests(unittest.TestCase):
    def fixture(self,root):
        raw=b'namespace Example\ntheorem beyondPacket : True := by trivial\nend Example\n'
        (root/'Example.lean').write_bytes(raw)
        tid='lean_theorem:Example.lean:beyondPacket'
        source={'id':'source:lean','path':'Example.lean','sha256':sha(raw),'line_start':1,'line_end':3}
        node={'id':tid,'type':'lean_theorem','label':'Example.beyondPacket',
              'summary':'theorem beyondPacket : True',
              'attrs':{'path':'Example.lean','line':2,'source_refs':['source:lean'],
                       'declaration_kind':'theorem','namespace':'Example',
                       'source_role':'exploratory','scope':'local_private',
                       'kernel_checked':False,'declared_axiom':False,'contains_sorry':False}}
        packet={'claim_id':'FORM.EXAMPLE.01','related_claims':[],'scripts':[],
                'sources':[{**source,'text':raw.decode()}],'lean_candidates':[]}
        semantic={key:[] for key in ['results','premises','derivations','input_audit',
                  'proof_gaps','objects','script_analysis','relationships']}
        semantic['lean_links']=[{'theorem_id':tid,'relation':'candidate','result_ids':[],
                                 'source_refs':['source:lean']}]
        return packet,semantic,{tid:node},{'source:lean':source}

    def test_audit_accepts_cited_current_full_registry_declaration_without_widening_api_packet(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);packet,semantic,lean,sources=self.fixture(root)
            original=copy.deepcopy(packet)
            with self.assertRaisesRegex(ValueError,'unknown theorem id'):validate(semantic,packet)
            candidates=audit.validation_lean_candidates(packet,semantic,lean,sources,root)
            validate(semantic,{**packet,'lean_candidates':candidates})
            self.assertEqual(packet,original)
            self.assertEqual(candidates[0]['namespace'],'Example')
            self.assertEqual(candidates[0]['source_role'],'exploratory')
            self.assertEqual(candidates[0]['scope'],'local_private')
            self.assertFalse(candidates[0]['kernel_checked'])
            self.assertEqual(candidates[0]['source_refs'],['source:lean'])
            self.assertNotIn(semantic['lean_links'][0]['theorem_id'],schema_for_packet(packet)['properties']['lean_links']['items']['properties']['theorem_id'].get('enum',[]))
            self.assertEqual(audit.validation_lean_candidates(packet,{'lean_links':[]},lean,sources,root),[])

    def test_absent_wrong_type_stale_and_out_of_scope_declarations_still_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);packet,semantic,lean,sources=self.fixture(root)
            tid=semantic['lean_links'][0]['theorem_id']
            cases=[]
            cases.append(('absent',{},sources))
            wrong=copy.deepcopy(lean);wrong[tid]['type']='claim';cases.append(('wrong type',wrong,sources))
            stale=copy.deepcopy(sources);stale['source:lean']['sha256']='0'*64;cases.append(('stale',lean,stale))
            outside=copy.deepcopy(lean);outside[tid]['attrs']['line']=4;cases.append(('outside range',outside,sources))
            moved=copy.deepcopy(sources);moved['source:lean']['path']='Other.lean';cases.append(('different owner',lean,moved))
            for name,index,registry in cases:
                with self.subTest(name=name):
                    candidates=audit.validation_lean_candidates(packet,semantic,index,registry,root)
                    with self.assertRaisesRegex(ValueError,'unknown theorem id'):
                        validate(semantic,{**packet,'lean_candidates':candidates})

    def test_file_edit_after_first_check_is_not_accepted_from_a_previous_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);packet,semantic,lean,sources=self.fixture(root)
            self.assertEqual(len(audit.validation_lean_candidates(packet,semantic,lean,sources,root)),1)
            (root/'Example.lean').write_text('changed declaration\n',encoding='utf-8')
            self.assertEqual(audit.validation_lean_candidates(packet,semantic,lean,sources,root),[])


if __name__ == "__main__":
    unittest.main()
