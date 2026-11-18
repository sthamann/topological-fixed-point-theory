import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import extract


class GraphIntegrityTests(unittest.TestCase):
    def test_unresolved_semantic_edges_and_self_relations_are_visible_not_proof_links(self):
        row = {key: '' for key in ['claim_id', 'claim', 'status', 'location',
            'dependencies', 'script', 'external_data', 'supersedes', 'canonical_status', 'active']}
        row.update(claim_id='TEST.ONE.01', claim='A conditional calculation', active='true')
        sem = {
            'premises': [{'id': 'P1', 'text': 'Chosen input', 'origin_claim_id': ''}],
            'results': [{'id': 'R1', 'label': 'Output', 'statement': 'Conditional output'}],
            'derivations': [{'id': 'D1', 'label': 'Calculation', 'argument': 'Conditional',
                'premise_ids': ['P1'], 'upstream_claim_ids': ['MISSING.01'],
                'result_ids': ['R1'], 'scope': 'under supplied assumptions'}],
            'relationships': [{'target_claim_id': 'TEST.ONE.01', 'relation': 'uses_result'}],
        }
        with tempfile.TemporaryDirectory() as directory:
            corpus = SimpleNamespace(rows=[row], root=Path(directory), script_nodes={}, lean_nodes={}, sources={})
            graph = extract.build_graph(corpus, {'TEST.ONE.01': {'semantic': sem}}, 'fixture')
        issues = graph['meta']['graph_integrity']['issues']
        self.assertEqual({x['reason'] for x in issues}, {'unresolved_semantic_reference', 'self_relation'})
        self.assertEqual(graph['meta']['graph_integrity']['unresolved_semantic_links'], 2)
        self.assertTrue(any(x['missing_endpoints'] == ['claim:MISSING.01'] for x in issues))
        self.assertFalse(any(e['src'] == e['dst'] for e in graph['edges']))
        self.assertEqual(sum(e['type'] == 'required_together' for e in graph['edges']), 1)
        self.assertEqual(sum(e['type'] == 'derives' for e in graph['edges']), 1)
        claim = next(n for n in graph['nodes'] if n['id'] == 'claim:TEST.ONE.01')
        self.assertEqual(claim['attrs']['semantic'], sem)


if __name__ == '__main__':
    unittest.main()
