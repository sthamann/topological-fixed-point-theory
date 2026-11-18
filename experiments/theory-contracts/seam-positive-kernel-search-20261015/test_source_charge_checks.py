"""Regression tests for the scoped source charge argument and its controls."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import source_charge_checks as subject


class SourceChargeChecks(unittest.TestCase):
    def test_expression_source_typing_and_tags(self):
        p = lambda name: {'op': 'primitive', 'name': name}
        scalar = subject.assess_expression(p('scalar_quarter_current'))
        self.assertEqual(scalar['verdict'], 'excluded_scoped')
        bilinear = {'op': 'product', 'args': [p('matter_spinor'),
                    {'op': 'adjoint', 'arg': p('matter_spinor')}]}
        self.assertEqual(subject.assess_expression(bilinear)['tags'], [[0, 0]])
        self.assertEqual(subject.assess_expression(p('clifford_vector'))['verdict'], 'typing_error')
        mixed = {'op': 'sum', 'args': [p('matter_spinor'), p('tensor_boson')]}
        self.assertEqual(subject.assess_expression(mixed)['tags'], [[0, 0], [1, 1]])
        sector = {'op': 'proposed_sector_insertion', 'source_ref': 'proposal:disorder'}
        candidate = subject.assess_expression({'op': 'derivative', 'arg': sector})
        self.assertEqual(candidate['verdict'], 'needs_original_source')
        self.assertTrue(candidate['possible_target_tag'])
        self.assertFalse(candidate['source_binding_verified'])
        self.assertFalse(candidate['candidate_physical_proof'])
        self.assertEqual(candidate['sector_insertion_status'], 'unverified_original_binding')
        self.assertFalse(scalar['source_binding_verified'])
        # An untyped new field might have tag (1,0); multiplying by a
        # matter spinor could then yield (0,1). Never falsely exclude it.
        composite = subject.assess_expression({'op': 'product',
                                              'args': [sector, p('matter_spinor')]})
        self.assertEqual(composite['verdict'], 'needs_original_source')
        self.assertTrue(composite['possible_target_tag'])
        self.assertEqual(composite['tags'], [[0, 0], [0, 1], [1, 0], [1, 1]])

    def test_expression_malformed_and_complexity_limits(self):
        bad = [None, True, {'op': 'unknown'}, {'op': True},
               {'op': 'primitive', 'name': True},
               {'op': 'proposed_sector_insertion', 'source_ref': True},
               {'op': 'sum', 'args': []}, {'op': 'product', 'args': [False]},
               {'op': 'primitive', 'name': 'tensor_boson', 'extra': 0}]
        for expr in bad:
            with self.subTest(expr=expr):
                self.assertEqual(subject.assess_expression(expr)['verdict'], 'typing_error')
        atom = {'op': 'primitive', 'name': 'tensor_boson'}
        self.assertEqual(subject.assess_expression({'op': 'sum', 'args': [atom]*100})['verdict'], 'typing_error')
        for _ in range(24):
            atom = {'op': 'adjoint', 'arg': atom}
        self.assertEqual(subject.assess_expression(atom)['verdict'], 'typing_error')

    def test_original_controls_and_scope(self):
        result = subject.check()
        self.assertEqual(result['status'], 'pass')
        self.assertEqual(len(result['checks']), 19)
        self.assertEqual(result['tag_class'], [[0, 0], [1, 1]])
        self.assertEqual(result['target_tag'], [0, 1])
        self.assertFalse(result['physical_source_identified'])
        self.assertFalse(result['fundamental_solution'])
        self.assertIn('supplied', result['source_scope'])
        self.assertEqual([x['path'] for x in result['source_paths']], [
            'verification/v252_full_finite_triple.py', 'verification/v1_e8_glue.py'])

    def test_disorder_escape_and_neutral_pair(self):
        diagonal = {(0, 0), (1, 1)}
        self.assertEqual(subject.add((1, 1), (1, 1)), (0, 0))
        self.assertNotIn((0, 1), diagonal)
        escaped = diagonal | {subject.add(tag, (0, 1)) for tag in diagonal}
        self.assertEqual(escaped, {(0, 0), (0, 1), (1, 0), (1, 1)})

    def test_failed_condition_raises_under_optimized_python(self):
        code = ('import source_charge_checks as s; '
                's.need(False, "negative control", [])')
        proc = subprocess.run([sys.executable, '-B', '-OO', '-c', code],
                              cwd=Path(subject.__file__).parent, capture_output=True, text=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn('Failed: negative control', proc.stderr)

    def test_cli_result_and_exclusive_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'result.json'
            command = [sys.executable, '-B', '-OO', subject.__file__, '--out', str(output)]
            first = subprocess.run(command, capture_output=True, text=True, check=True)
            self.assertEqual(json.loads(first.stdout), json.loads(output.read_text()))
            second = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(second.returncode, 0)
            self.assertEqual(json.loads(output.read_text())['status'], 'pass')

    def test_three_source_controls_run_through_existing_search(self):
        from origin_hypothesis_search import run_checks, ALLOWED
        path = Path(subject.__file__).with_name('source_genesis_checks.json')
        catalog = {row['id']: row for row in json.loads(path.read_text())['checks']}
        rows = run_checks({'tests': list(catalog)}, catalog)
        self.assertEqual(len(rows), 3)
        for row in rows:
            self.assertEqual(row['output']['inference_claims'],
                             [{'id': ALLOWED[row['test_id']], 'holds': False}])
            self.assertFalse(row['output']['physical_source_identified'])


if __name__ == '__main__':
    unittest.main()
