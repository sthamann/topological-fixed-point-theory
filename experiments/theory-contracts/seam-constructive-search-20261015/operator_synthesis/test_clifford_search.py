"""Exact source-map recovery, discriminating counterexamples and cache regressions."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('clifford_search',HERE/'clifford_search.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)


def equal(a,b):
    return all(np.array_equal(x,y) for x,y in zip(a,b))


def majorana(j):
    annihilator=c.annihilation(j//2)
    creator=c.adjoint(annihilator)
    if j%2: return c.phase(c.subtract(creator,annihilator),1)
    return annihilator[0]+creator[0],annihilator[1]+creator[1]


class CliffordSearchTests(unittest.TestCase):
    def test_original_CAR_and_majorana_algebra(self):
        identity=c.pair(np.eye(c.DIM,dtype=np.int64))
        zeros=c.pair(np.zeros((c.DIM,c.DIM),dtype=np.int64))
        for j in range(10):
            g=majorana(j)
            self.assertTrue(equal(g,c.adjoint(g)))
            self.assertTrue(equal(c.multiply(g,g),identity))
            for k in range(j):
                a,b=c.multiply(g,majorana(k)),c.multiply(majorana(k),g)
                self.assertTrue(equal((a[0]+b[0],a[1]+b[1]),zeros))

    def test_word_canonicalization_sign_and_cancellation(self):
        self.assertEqual(c.canonical_word([2,1]),(-1,(1<<1)|(1<<2)))
        self.assertEqual(c.canonical_word([2,1,2]),(-1,1<<1))
        self.assertEqual(c.canonical_word([3,3]),(1,0))
        self.assertEqual(c.canonical_word([1,2,1,2]),(-1,0))
        with self.assertRaises(ValueError): c.canonical_word([10])

    def test_monomial_action_against_independent_dense_exact_products(self):
        identity=c.pair(np.eye(c.DIM,dtype=np.int64))
        matrix=c.fixture_matrices(c.fixture_key({'kind':'one_body','family':'imaginary_antisymmetric','j':0,'k':2}),False)[0]
        for word in ([],[1],[2,1,2],[9,2,0,5,7],[3,3,4,1]):
            sign,support=c.canonical_word(word)
            u=identity
            for j in word: u=c.multiply(u,majorana(j))
            for anti in (False,True):
                source=(matrix[0],-matrix[1]) if anti else matrix
                expected=c.multiply(c.multiply(u,source),c.adjoint(u))
                self.assertTrue(equal(c.transform({'support':support,'antilinear':anti},matrix),expected))

    def test_generation_has_all_structures_and_both_parities(self):
        candidates=list(c.generate_candidates(19))
        self.assertEqual(len(candidates),2048)
        self.assertEqual(len({c.candidate_key(g) for g in candidates}),2048)
        self.assertEqual(sum(g['support'].bit_count()%2 for g in candidates),1024)
        self.assertEqual(sum(g['antilinear'] for g in candidates),1024)
        self.assertEqual(len(list(c.one_body_basis())),25)

    def test_cegis_generates_then_learns_imaginary_counterexample(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=c.search(tmp,seed=13)
            self.assertEqual(len(result['rounds']),2)
            first=result['rounds'][0]['survivors']
            self.assertEqual(len(first),2)
            self.assertEqual({g['antilinear'] for g in first},{False,True})
            self.assertEqual(len(result['survivors']),1)
            survivor=result['survivors'][0]
            self.assertTrue(survivor['antilinear'])
            self.assertIsNone(c.verify_complete(survivor))
            self.assertEqual(len(result['counterexample_trace']),1)
            failure=result['counterexample_trace'][0]
            self.assertEqual(failure['counterexample']['fixture']['family'],'imaginary_antisymmetric')
            self.assertNotEqual(failure['counterexample']['actual'],failure['counterexample']['expected'])
            self.assertNotEqual(failure['old_fixture_hash'],failure['new_fixture_hash'])
            self.assertEqual(result['stale_cache_entries_invalidated'],2048)
            self.assertFalse(result['source_selection_proved'])

    def test_erased_anti_flag_and_wrong_parity_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            winner=c.search(tmp)['survivors'][0]
        linear={**winner,'antilinear':False}
        self.assertTrue(all(c.counterexample(linear,f) is None for f in c.initial_fixtures()))
        self.assertEqual(c.verify_complete(linear)['fixture']['family'],'imaginary_antisymmetric')
        even={**winner,'support':winner['support']^1}
        self.assertEqual(even['support'].bit_count()%2,0)
        self.assertIsNotNone(c.verify_complete(even))

    def test_quarter_coefficients_are_conjugated_for_antiunitarity(self):
        with tempfile.TemporaryDirectory() as tmp:
            winner=c.search(tmp)['survivors'][0]
        quarter={'kind':'quarter'}
        for anti in (False,True):
            candidate={**winner,'antilinear':anti}
            self.assertIsNone(c.counterexample(candidate,quarter))
        matrix,linear_expected=c.fixture_matrices(c.fixture_key(quarter),False)
        anti_actual=c.transform(winner,matrix)
        self.assertFalse(equal(anti_actual,linear_expected))

    def test_persistent_constraints_and_stale_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            first=c.search(tmp)
            second=c.search(tmp)
            self.assertEqual(second['survivors'],first['survivors'])
            self.assertEqual(second['persisted_constraints_reused'],10)
            self.assertEqual(second['cache_hits'],2048)
            self.assertEqual(second['counterexample_trace'],first['counterexample_trace'])
            path=Path(tmp)/'separate-cache.json'
            cache=c.ConstraintCache(path)
            linear={**first['survivors'][0],'antilinear':False}
            old=c.initial_fixtures()
            self.assertIsNone(cache.check(linear,old));cache.save()
            loaded=c.ConstraintCache(path)
            updated=old+[first['counterexample_trace'][0]['counterexample']['fixture']]
            self.assertIsNotNone(loaded.check(linear,updated))
            self.assertEqual(loaded.hits,0)
            self.assertEqual(loaded.invalidated,1)

    def test_fresh_process_recovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            process=subprocess.run([sys.executable,str(HERE/'clifford_search.py'),'--run-dir',tmp,'--seed','7'],capture_output=True,text=True,check=True)
            result=json.loads(process.stdout)
            self.assertEqual(result['candidate_count'],2048)
            self.assertEqual(len(result['survivors']),1)
            self.assertTrue(result['survivors'][0]['antilinear'])
            stored=json.loads((Path(tmp)/'results.json').read_text())
            self.assertEqual(stored['fixture_sha256'],c.digest(stored['active_fixtures']))


if __name__=='__main__': unittest.main()
