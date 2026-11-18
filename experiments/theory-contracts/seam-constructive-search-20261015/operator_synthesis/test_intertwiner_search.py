import json
from pathlib import Path
import tempfile
import unittest

import sympy as s

from intertwiner_search import dump_matrix, equation, load_problem, run, scalar, solve_branch


class InverseSynthesisTests(unittest.TestCase):
    def test_vectorization_and_antilinearity(self):
        a = s.Matrix([[1, s.I], [-s.I, 2]])
        b = s.Matrix([[2, 1], [0, 1]])
        x = s.Matrix([[1, 2], [s.I, -1]])
        for anti in (False, True):
            expected = x * (a.conjugate() if anti else a) - b*x
            self.assertEqual(equation(a, b, anti) * s.Matrix(list(x)), s.Matrix(list(expected)))

    def test_recovers_unseeded_unitary_from_operator_pairs(self):
        h = s.Matrix([[s.Rational(3,5), -s.Rational(4,5)], [s.Rational(4,5), s.Rational(3,5)]])
        z = s.diag(1, -1)
        x = s.Matrix([[0,1],[1,0]])
        y = s.Matrix([[0,-s.I],[s.I,0]])
        relations = [(name, a, h*a*h.H) for name,a in [('z',z),('x',x),('y',y)]]
        raw = {'premises':['synthetic exact recovery control']}
        result = solve_branch(raw, relations, 2, False, 10)
        self.assertEqual(result['status'], 'exact_intertwiner')
        self.assertEqual(result['intertwiner_dimension'], 1)
        self.assertTrue(result['all_relation_residuals_zero'])
        self.assertGreaterEqual(len(result['rounds']), 1)
        anti = solve_branch(raw, relations, 2, True, 10)
        self.assertEqual(anti['status'], 'no_nonzero_intertwiner')
        self.assertEqual(set(anti['core']), {'z','x','y'})

    def test_missing_constraints_do_not_inflate_final_family_dimension(self):
        # Identity immediately satisfies everything, but the whole commutant
        # is scalar, not the initial four-dimensional space of all matrices.
        x, z = s.Matrix([[0,1],[1,0]]), s.diag(1,-1)
        r = solve_branch({}, [('x',x,x),('z',z,z)], 2, False, 0)
        self.assertEqual(r['intertwiner_dimension'], 1)

    def test_rank_search_never_claims_nonexistence(self):
        a, b = s.diag(1,2), s.diag(1,3)
        r = solve_branch({}, [('spectral',a,b)], 2, False, 0)
        self.assertEqual(r['status'], 'full_rank_search_inconclusive')
        self.assertEqual(r['intertwiner_dimension'], 1)

    def test_zero_intertwiner_core(self):
        r = solve_branch({}, [('zero',s.zeros(2),s.eye(2))], 2, False, 0)
        self.assertEqual(r['status'], 'no_nonzero_intertwiner')
        self.assertEqual(r['core'], ['zero'])

    def test_cache_key_contains_changed_specification(self):
        z, x = s.diag(1,-1), s.Matrix([[0,1],[1,0]])
        a = solve_branch({'version':1}, [('z',z,x)], 2, False, 2)
        b = solve_branch({'version':2}, [('z',z,x)], 2, False, 2)
        self.assertNotEqual(a['rounds'][0]['key'], b['rounds'][0]['key'])
        c = solve_branch({'version':1}, [('z',z,x)], 2, False, 3)
        self.assertNotEqual(a['rounds'][0]['key'], c['rounds'][0]['key'])
        d = solve_branch({'version':1}, [('z',z,-x)], 2, False, 2)
        self.assertNotEqual(a['rounds'][0]['key'], d['rounds'][0]['key'])

    def test_json_never_evaluates_expressions_or_accepts_float(self):
        for value in [0.5, 'sqrt(2)', '__import__("os")', True, '1/0']:
            with self.assertRaises(ValueError): scalar(value)
        self.assertEqual(scalar(['1/2','-3/4']), s.Rational(1,2)-s.I*s.Rational(3,4))

    def test_end_to_end_json_and_scientific_scope(self):
        z, x = s.diag(1,-1), s.Matrix([[0,1],[1,0]])
        data = {'schema':'tfpt.operator_intertwiner.v1','id':'control',
                'premises':['synthetic test only'],
                'relations':[{'id':'z','source':dump_matrix(z),'target':dump_matrix(x)}]}
        with tempfile.TemporaryDirectory() as tmp:
            src, out = Path(tmp)/'input.json', Path(tmp)/'output.json'
            src.write_text(json.dumps(data))
            r = run(src,out,5)
            self.assertTrue(out.exists())
            self.assertFalse(r['physical_solution'])
            self.assertFalse(r['source_selection_proved'])
            self.assertEqual(r['branches'][0]['status'],'exact_intertwiner')
            data['target_map'] = dump_matrix(s.eye(2))
            src.write_text(json.dumps(data))
            with self.assertRaises(ValueError): load_problem(src)

    def test_similarity_is_not_unitary_equivalence(self):
        a = s.Matrix([[0,1],[0,0]])
        data = {'schema':'tfpt.operator_intertwiner.v1','id':'nonnormal-control',
                'premises':['Similarity alone is not enough'],
                'relations':[{'id':'nilpotent','source':dump_matrix(a),'target':dump_matrix(2*a)}]}
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'operators.json'
            p.write_text(json.dumps(data))
            r = run(p, trials=5)
            self.assertTrue(all(b['status']=='no_nonzero_intertwiner' for b in r['branches']))


if __name__ == '__main__':
    unittest.main()
