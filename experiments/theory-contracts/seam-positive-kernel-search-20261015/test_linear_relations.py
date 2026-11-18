"""Synthetic relation controls; these are not TFPT source relations."""
import json
import tempfile
import unittest
from pathlib import Path

import sympy as sp
from kernel_search import MomentSearch, load_source_control
from linear_relations import normalize_linear_relations, check_linear_relations
from run_search import campaign
from word_algebra import T
from finite_completion import complete
from test_finite_completion import fixture


def relation(target):
    return {'terms': [{'word': [0,T,0], 're': '1'}, {'word': [1,T,1], 're': '1'}],
            're': str(target)}


class LinearRelationTests(unittest.TestCase):
    def test_canonicalization_preserves_car_sign_and_combines(self):
        relation_input = [{'terms': [{'word': [1,0], 're': '2/3', 'im': '1/3'},
                                     {'word': [0,1], 're': '1/6'}], 're': '0'}]
        result = normalize_linear_relations(1, 2, relation_input)
        self.assertEqual(result[0]['terms'], [{'word': [0,1], 're': '-1/2', 'im': '-1/3'}])
        for invalid in ('__import__("os")', 0.5, True, '1/0'):
            with self.assertRaises(ValueError):
                normalize_linear_relations(1,2,[{'terms':[{'word':[],'re':invalid}],'re':0}])

    def test_numeric_feasible_and_conflicting_conditions(self):
        good = MomentSearch(1,2,linear_relations=[relation(1)]).solve(solver='SCS')
        self.assertIn(good['solver_status'], ('optimal', 'optimal_inaccurate'))
        self.assertEqual(good['linear_relation_count'], 1)
        self.assertLess(good['max_constraint_violation'], 1e-6)
        bad = MomentSearch(1,2,linear_relations=[relation(1),relation(3)]).solve(solver='SCS')
        self.assertEqual(bad['solver_status'], 'infeasible')
        self.assertNotIn('moments', bad)
        self.assertEqual(bad['verdict'], 'NO_NUMERICAL_CANDIDATE')

    def test_exact_replay_checks_nontrivial_complex_relations(self):
        rows=[{'terms':[{'word':[0], 're':'1/3', 'im':'1/2'},
                        {'word':[1], 're':'-1/2'}],
               're':'-2/3', 'im':'1/2'}]
        self.assertEqual(check_linear_relations(rows,{(0,):sp.Integer(1),(1,):sp.Integer(2)},1,2),1)
        with self.assertRaisesRegex(ValueError,'fails exact replay'):
            check_linear_relations(rows,{(0,):sp.Integer(1),(1,):sp.Integer(3)},1,2)
        with self.assertRaisesRegex(ValueError,'missing exact word'):
            check_linear_relations(rows,{(0,):sp.Integer(1)},1,2)

    def test_campaign_reads_and_preserves_relations_and_legacy_input(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp)
            source=base/'input.json'
            source.write_text(json.dumps({'linear_relations':[relation(1)],
                                          'source_observables':{'label':'synthetic control'}}))
            control=load_source_control(source)
            result=campaign(base/'run',m=1,degree=2,**control)
            self.assertEqual(len(result['configuration']['linear_relations']),1)
            payload=json.loads((base/'run'/'moments.json').read_text())
            self.assertEqual(payload['linear_relations'],result['configuration']['linear_relations'])
            self.assertIn(result['solver_status'], ('optimal','optimal_inaccurate'))
            legacy=base/'old.json';legacy.write_text(json.dumps({'fixed_moments':[]}))
            self.assertEqual(load_source_control(legacy)['linear_relations'],[])

    def test_exact_completion_refuses_a_violated_relation(self):
        payload=fixture()
        payload['linear_relations']=[relation(1)]
        self.assertTrue(complete(payload)['ok'])
        payload['linear_relations']=[relation(3)]
        result=complete(payload)
        self.assertFalse(result['ok'])
        self.assertIn('fails exact replay',result['reason'])


if __name__ == '__main__':
    unittest.main()
