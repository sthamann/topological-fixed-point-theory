"""Controller controls only. Fake backend is never the research backend."""
from pathlib import Path
import json
import tempfile
import unittest
import hashlib
import sys

from backend import command
from context import write_json
from launch import run, rows, ADMISSION, target_assessment, followup, decision_history
from prompts import response_schema, GOAL_ANCHOR


class FakeControlBackend:
    calls = []
    fail_critic_once = False

    def __init__(self, stop, update):
        self.stop, self.update = stop, update

    def terminate(self):
        pass

    def __call__(self, role, folder, manifest, instruction):
        self.calls.append(role)
        folder.mkdir(parents=True)
        if role == 'critic' and self.fail_critic_once:
            type(self).fail_critic_once = False
            raise RuntimeError('CONTROL_ONLY simulated provider failure')
        schema = response_schema(role)
        result = {k: [] if v['type'] == 'array' else 'CONTROL_ONLY' for k, v in schema['properties'].items()}
        result['proposed_status' if role == 'proposer' else 'verdict'] = (
            'candidate_complete_solution' if role == 'proposer' else 'candidate_requires_proof_review')
        metrics = {'role': role, 'call_dir': str(folder), 'completed': True,
                   'tool_events': 1, 'completed_commands': 1, 'failed_commands': 0,
                   'receipt': {'verified': True, 'control_only': True}}
        write_json(folder / 'metrics.json', metrics)
        return result, metrics


def control_context(repo, out):
    return {'context_version': 'CONTROL_ONLY', 'python': '/usr/bin/python3',
            'scientific_readonly_roots': [], 'exact_source_files': []}


class Controls(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.root = Path(self.temp.name)
        self.repo = self.root / 'repo'; self.repo.mkdir()
        self.out = self.root / 'out'
        FakeControlBackend.calls = []; FakeControlBackend.fail_critic_once = False

    def tearDown(self):
        self.temp.cleanup()

    def execute(self, count):
        return run(self.repo, self.out, count, backend_factory=FakeControlBackend,
                   context_loader=control_context)

    def test_pinned_model_scoped_permissions_no_fallback(self):
        folder = self.root / 'out/rounds/id/proposer-001'
        argv, policy = command(folder, control_context(self.repo, self.out))
        self.assertEqual(argv[argv.index('--model') + 1], 'gpt-6.1-sol')
        for flag in ('--ignore-user-config', '--ignore-rules', '--ephemeral', '--json'):
            self.assertIn(flag, argv)
        self.assertIn('approval_policy="never"', argv)
        self.assertIn('model_provider="openai"', argv)
        self.assertIn('mcp_servers={}', argv)
        self.assertFalse(policy['network']['enabled'])
        self.assertNotIn(str(Path.home()), policy['filesystem'])
        self.assertEqual(policy['filesystem'][str(folder / 'scratch')], 'write')
        self.assertEqual(policy['filesystem']['/opt/homebrew/bin/ollama'], 'deny')
        self.assertEqual(policy['filesystem'][str(Path('/opt/homebrew/bin/ollama').resolve())], 'deny')
        with self.assertRaises(ValueError):
            command(folder, control_context(self.repo, self.out), model='qwen')

    def test_complete_candidate_never_auto_physics_pass_and_resume(self):
        self.assertEqual(self.execute(1), 0)
        self.assertEqual(FakeControlBackend.calls, ['proposer', 'critic'])
        first = rows(self.out / 'research-rounds.jsonl')[0]
        self.assertEqual(first['physics_admission'], ADMISSION)
        self.assertIn('Independently verify', first['next_task'])
        self.assertEqual(self.execute(2), 0)
        self.assertEqual(FakeControlBackend.calls, ['proposer', 'critic', 'proposer', 'critic'])
        self.assertEqual(len(rows(self.out / 'research-rounds.jsonl')), 2)

    def test_provider_failure_preserves_proposal_and_resumes_critic(self):
        FakeControlBackend.fail_critic_once = True
        self.assertEqual(self.execute(1), 1)
        worker = json.loads((self.out / 'worker.json').read_text())
        self.assertEqual(worker['status'], 'error')
        self.assertEqual(self.execute(1), 0)
        self.assertEqual(FakeControlBackend.calls, ['proposer', 'critic', 'critic'])
        self.assertEqual(len(rows(self.out / 'research-rounds.jsonl')), 1)

    def test_partial_last_line_preserves_all_complete_history(self):
        self.out.mkdir()
        p = self.out / 'research-rounds.jsonl'; p.write_bytes(b'{"id":1}\n{"id":')
        self.assertEqual(rows(p), [{'id': 1}])
        self.assertEqual(p.read_text(), '{"id":1}\n')

    def test_rephrased_conditional_math_cannot_continue_unchanged_blocker(self):
        proposal = {'proposed_status': 'new_derivation',
                    'new_information': 'Neue korrekte lokale Formel im anderen Chart.',
                    'next_task': 'Noch eine lokale Formel'}
        review = {'verdict': 'scoped_result_survives', 'target_progress': 'none',
                  'blocker_class': 'JH_test_pairing_domain',
                  'next_task': 'Weitere interne Identität berechnen'}
        task = followup(proposal, review, [], {'effective_progress': 'none'})
        self.assertTrue(task.startswith(GOAL_ANCHOR))
        self.assertIn('Bearbeite den benannten fehlenden ursprünglichen Input selbst', task)
        self.assertNotEqual(task, review['next_task'])

    def test_presupposed_input_is_not_counted_as_derived(self):
        scratch = self.root / 'call/scratch'; scratch.mkdir(parents=True)
        (scratch / 'proof.md').write_text('CONTROL ONLY conditional lemma')
        proposal = {'progress_evidence': ['proof.md']}
        review = {'target_progress': 'input_derived', 'verdict': 'scoped_result_survives',
                  'premises_presupposing_target': ['CONTROL_ONLY assumed target'],
                  'progress_evidence': ['proof.md']}
        call = {'call_dir': str(scratch.parent)}
        result = target_assessment(proposal, review, call, call)
        self.assertEqual(result['effective_progress'], 'none')
        self.assertEqual(result['physics_admission'], ADMISSION)

    def test_candidate_claim_needs_files_values_and_observed_execution(self):
        scratch = self.root / 'candidate/scratch'; scratch.mkdir(parents=True)
        proposal = {'candidate_id': 'CONTROL', 'candidate_bundle': 'bundle.json',
                    'progress_evidence': ['bundle.json']}
        review = {'target_progress': 'candidate_computed', 'verdict': 'scoped_result_survives'}
        call = {'call_dir': str(scratch.parent), 'commands': []}
        self.assertEqual(target_assessment(proposal, review, call, call)['effective_progress'], 'none')
        write_json(scratch / 'bundle.json', {'candidate_id': 'CONTROL',
                   'producer_file': 'producer.py', 'result_file': 'values.json',
                   'slots': {k: {'status': 'hypothesis', 'description': 'CONTROL ONLY',
                                  'new_assumption': 'Not physical evidence'} for k in (
                       'regulator_family', 'field_operator', 'reference_measure',
                       'marked_time_dictionary', 'absolute_reference', 'whole_contributions')}})
        (scratch / 'producer.py').write_text('# CONTROL ONLY, no physics')
        write_json(scratch / 'values.json', {'unnormalized_values': {k: 1 for k in ['Z0','ZU','ZY','ZUY']}})
        self.assertEqual(target_assessment(proposal, review, call, call)['effective_progress'], 'none')
        token = 'RESULT_SHA256=' + hashlib.sha256((scratch / 'values.json').read_bytes()).hexdigest()
        call['commands'] = [{'command': sys.executable + ' -B producer.py', 'exit_code': 0, 'output': token}]
        result = target_assessment(proposal, review, call, call)
        self.assertEqual(result['effective_progress'], 'candidate_computed')
        self.assertEqual(result['physics_admission'], ADMISSION)
        for command in ('cat producer.py', 'echo python producer.py', "/bin/zsh -lc 'cat producer.py'"):
            call['commands'][0]['command'] = command
            self.assertEqual(target_assessment(proposal, review, call, call)['effective_progress'], 'none')
        call['commands'][0]['command'] = sys.executable + ' -B producer.py'
        (scratch / 'values.json').write_text((scratch / 'values.json').read_text() + '\n')
        self.assertEqual(target_assessment(proposal, review, call, call)['effective_progress'], 'none')
        proposal['candidate_bundle'] = '../outside.json'
        self.assertEqual(target_assessment(proposal, review, call, call)['effective_progress'], 'none')

    def test_old_history_is_explicitly_unknown_and_in_prompt(self):
        row = {'round_id': 'research-legacy', 'proposal': {'question': 'Old question'},
               'review': {'first_extra_assumption': 'Old open input'}}
        table = decision_history([row])
        self.assertEqual(table[0]['blocker_class'], 'unknown')
        self.assertEqual(table[0]['reported_target_progress'], 'unknown')
        from prompts import prompt
        rendered = prompt('proposer', {}, 'next', self.out / 'journal', decision_history=table)
        self.assertIn('Old open input', rendered)

    def test_empty_legacy_checkpoint_keeps_old_question_and_uses_new_goal(self):
        p = self.out / 'rounds/research-000001/round-state.json'
        p.parent.mkdir(parents=True)
        write_json(p, {'round_id': 'research-000001', 'question': 'Old local chart task',
                       'started_utc': 'CONTROL_ONLY', 'context_version': 'OLD'})
        self.assertEqual(self.execute(1), 0)
        state = json.loads(p.read_text())
        self.assertTrue(state['question'].startswith(GOAL_ANCHOR))
        self.assertEqual(state['prior_questions'][0]['question'], 'Old local chart task')

    def test_known_r055_rediscovery_does_not_count_as_new_target_progress(self):
        scratch = self.root / 'call/scratch'; scratch.mkdir(parents=True)
        (scratch / 'proof.md').write_text('CONTROL ONLY reuse of old proof')
        proposal = {'proposed_status': 'scoped_exclusion', 'new_information': 'Rediscovery',
                    'next_task': 'More local algebra', 'candidate_id': 'CONTROL',
                    'progress_evidence': ['proof.md']}
        review = {'verdict': 'scoped_result_survives', 'target_progress': 'candidate_excluded',
                  'blocker_class': 'nu_Vn_base_law', 'next_task': 'More local algebra',
                  'already_implied_by': ['finding:bounded-fixed-rank-spectral-action-does-not-normalize-noncompact-higgs-volume']}
        call = {'call_dir': str(scratch.parent)}
        result = target_assessment(proposal, review, call, call)
        self.assertEqual(result['effective_progress'], 'none')
        self.assertNotIn('More local algebra', followup(proposal, review, [], result))

    def test_required_guardrail_absence_is_explicit(self):
        from context import prepare
        from unittest.mock import patch
        for relative in ('semantic_explorer/research/program.json', 'docs/research/fundamental-plan/progress.json'):
            p = self.repo / relative; p.parent.mkdir(parents=True, exist_ok=True); p.write_text('{}')
        with patch('context.build_knowledge_context', return_value={}):
            with self.assertRaisesRegex(FileNotFoundError, 'Required target-source findings missing'):
                prepare(self.repo, self.out)

    def test_round_boundary_stop_request_is_consumed_once(self):
        self.out.mkdir(); (self.out / 'stop-after-round').touch()
        self.assertEqual(self.execute(2), 0)
        self.assertEqual(len(rows(self.out / 'research-rounds.jsonl')), 1)
        self.assertFalse((self.out / 'stop-after-round').exists())
        self.assertEqual(self.execute(2), 0)
        self.assertEqual(len(rows(self.out / 'research-rounds.jsonl')), 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
