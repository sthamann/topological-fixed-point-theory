"""Autonomous original-source proposer -> independent critic -> next question.

No model call happens on import. No canonical scientific status is written.
"""
from pathlib import Path
import argparse
import fcntl
import hashlib
import json
import os
import shlex
import re
import math
import signal
import sys

from context import prepare, write_json
from backend import Backend, utc
from prompts import INITIAL_TASK, prompt, GOAL_ANCHOR, DIRECT_TARGET_TASK

ADMISSION = 'NOT_ADMITTED_BY_RUNTIME'


def rows(path):
    if not path.exists():
        return []
    # A crash can leave an incomplete final append. Repair only that line; never
    # discard completed prior science or silently skip interior corruption.
    result, valid_bytes = [], 0
    with path.open('rb') as stream:
        for line in stream:
            try:
                result.append(json.loads(line)); valid_bytes += len(line)
            except json.JSONDecodeError:
                if stream.read():
                    raise ValueError('Interior journal corruption: ' + str(path))
                with path.open('r+b') as repair:
                    repair.truncate(valid_bytes)
                break
    return result


def decision_history(history):
    """Small explicit continuity table; the complete journal remains available."""
    return [{'round_id': h['round_id'], 'question': h['proposal'].get('question'),
             'first_extra_assumption': h['review'].get('first_extra_assumption'),
             'blocker_class': h['review'].get('blocker_class', 'unknown'),
             'reported_target_progress': h['review'].get('target_progress', 'unknown'),
             'host_progress': h.get('target_assessment', {}).get('effective_progress', 'unknown'),
             'definition_decision': h['review'].get('definition_decision', 'not_assessed'),
             'route': h['review'].get('route', 'not_identified'),
             'candidate_id': h['proposal'].get('candidate_id', ''),
             'already_implied_by': h['review'].get('already_implied_by', [])}
            for h in history]


def producer_invocation(command, producer, scratch, python):
    """Recognize a direct Python script invocation, not cat/echo or shell text."""
    try:
        args = shlex.split(command)
        for _ in range(3):
            if len(args) == 3 and Path(args[0]).name in ('sh', 'bash', 'zsh') and args[1] in ('-c', '-lc'):
                args = shlex.split(args[2])
            else:
                break
        while args and re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*=.*', args[0]):
            if args[0].split('=', 1)[0] in ('PATH', 'PYTHONHOME'):
                return False
            args.pop(0)
        if not args or any(token in (';', '&&', '||', '|', '&', '>', '>>', '<') for token in args):
            return False
        executable = Path(args.pop(0))
        if not executable.is_absolute() and not re.fullmatch(r'python(?:[0-9.]+)?', executable.name):
            return False
        actual_python = executable if executable.is_absolute() else Path(python).parent / executable
        if actual_python.resolve() != Path(python).resolve():
            return False
        while args and args[0] in ('-B', '-u', '-I', '-E', '-s', '-S', '-O', '-OO', '-q'):
            args.pop(0)
        if not args or args[0].startswith('-'):
            return False
        return (scratch / args[0]).resolve() == producer.resolve()
    except (ValueError, TypeError, OSError):
        return False


def computed_value(value):
    if isinstance(value, bool):
        return False
    if isinstance(value, int):
        return True
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, list) and len(value) == 2:
        return all(isinstance(x, (int, float)) and computed_value(x) for x in value)
    return False


def target_assessment(proposal, review, proposal_call, review_call, python=sys.executable):
    """Evidence of task relevance, never a scientific admission or proof checker."""
    claimed = review.get('target_progress', 'none')
    effective, reasons = claimed, []
    if review.get('verdict') in ('falsified', 'unresolved'):
        effective = 'none'; reasons.append('Critic did not establish the proposed result')
    if claimed == 'input_derived' and review.get('premises_presupposing_target'):
        effective = 'none'; reasons.append('The claimed original input is still presupposed')
    if claimed in ('input_derived', 'candidate_excluded', 'definition_decided') and review.get('already_implied_by'):
        effective = 'none'; reasons.append('Already implied by: ' + ', '.join(review['already_implied_by']))
    evidence = []
    for response, call in ((proposal, proposal_call), (review, review_call)):
        scratch = Path(call.get('call_dir', '/nonexistent')) / 'scratch'
        for name in response.get('progress_evidence', []):
            p = (scratch / name).resolve()
            if p.is_relative_to(scratch.resolve()) and p.is_file():
                evidence.append(str(p))
    if claimed != 'none' and not evidence:
        effective = 'none'; reasons.append('No referenced local progress artifact was found')
    if claimed == 'definition_decided' and review.get('definition_decision') not in (
            'defined', 'not_found_in_scoped_sources'):
        effective = 'none'; reasons.append('No scoped definition decision was supplied')
    if claimed == 'candidate_excluded' and not (
            proposal.get('candidate_id') or review.get('candidate_id')):
        effective = 'none'; reasons.append('No concrete excluded candidate was identified')
    candidate = {'status': 'not_supplied'}
    if claimed == 'candidate_computed':
        scratch = Path(proposal_call.get('call_dir', '/nonexistent')) / 'scratch'
        def local(name):
            p = (scratch / name).resolve()
            if not name or not p.is_relative_to(scratch.resolve()) or not p.is_file():
                raise ValueError('Candidate artifact must be an existing file in proposer scratch')
            return p
        try:
            bundle_path = local(proposal.get('candidate_bundle', ''))
            bundle = json.loads(bundle_path.read_text())
            if not proposal.get('candidate_id') or bundle.get('candidate_id') != proposal['candidate_id']:
                raise ValueError('Candidate identity does not bind the supplied artifact')
            slots = bundle.get('slots', {})
            for key in ('regulator_family', 'field_operator', 'reference_measure',
                        'marked_time_dictionary', 'absolute_reference', 'whole_contributions'):
                slot = slots.get(key, {})
                if slot.get('status') not in ('original', 'hypothesis') or not slot.get('description'):
                    raise ValueError('Candidate has an open or undescribed source slot: ' + key)
                if slot['status'] == 'original' and not slot.get('source_references'):
                    raise ValueError('Claimed original slot lacks source references: ' + key)
                if slot['status'] == 'hypothesis' and not slot.get('new_assumption'):
                    raise ValueError('Hypothesis slot hides its extra assumption: ' + key)
            producer, result_file = local(bundle['producer_file']), local(bundle['result_file'])
            result = json.loads(result_file.read_text())
            values = result.get('unnormalized_values', {})
            if any(k not in values or not computed_value(values[k]) for k in ('Z0', 'ZU', 'ZY', 'ZUY')):
                raise ValueError('Candidate has no calculated four unnormalized contributions')
            commands = proposal_call.get('commands', [])
            digest = hashlib.sha256(result_file.read_bytes()).hexdigest()
            if not any(c.get('exit_code') == 0 and
                       producer_invocation(c.get('command', ''), producer, scratch, python) and
                       ('RESULT_SHA256=' + digest) in (c.get('output') or '') for c in commands):
                raise ValueError('No successful direct producer invocation emitting the current result hash')
            candidate = {'status': 'candidate_artifacts_and_execution_observed',
                         'candidate_id': proposal['candidate_id'], 'bundle': str(bundle_path),
                         'producer': str(producer), 'result_file': str(result_file),
                         'scope': 'File and execution evidence only; source identity, calculations and physics require independent review'}
        except (ValueError, KeyError, TypeError, OSError) as exc:
            effective = 'none'; candidate = {'status': 'incomplete_candidate', 'reason': str(exc)}
            reasons.append(str(exc))
    return {'reported_progress': claimed, 'effective_progress': effective,
            'blocker_class': review.get('blocker_class', 'unknown'),
            'evidence_files_found': sorted(set(evidence)), 'candidate': candidate,
            'reasons': reasons, 'physics_admission': ADMISSION,
            'scope': 'Task steering and observed artifacts, not verification of the scientific claim'}


def followup(proposal, review, history, assessment=None):
    task = review['next_task'].strip() or proposal['next_task'].strip() or INITIAL_TASK
    if proposal['proposed_status'] == 'candidate_complete_solution':
        return ('Independently verify the candidate complete solution against the actual whole-source '
                'and T1..T8 obligations. Actively seek a counterexample to its first unproved '
                'identification, rerun the actual original check or Lean statement. ' + task)
    if assessment is not None:
        progress = assessment['effective_progress']
        blocker = review.get('blocker_class', 'unknown')
        if review.get('verdict') == 'falsified':
            return (GOAL_ANCHOR + 'Die konkrete Kandidatenidentifikation wurde widerlegt. '
                    'Erhalte den Gegenbefund und ändere genau die widerlegte Annahme. ' + task)
        if progress == 'none':
            previous = [h['round_id'] for h in history
                        if h['review'].get('blocker_class') == blocker]
            return (GOAL_ANCHOR + 'Unverkleinerter Zielblocker: ' + blocker + '. '
                    + ('Bereits behandelt in: ' + ', '.join(previous) + '. ' if previous else '')
                    + DIRECT_TARGET_TASK)
        if progress == 'definition_decided' and review.get('definition_decision') == 'not_found_in_scoped_sources':
            return (GOAL_ANCHOR + 'Der protokollierte Quellenumfang liefert diesen Input nicht. '
                    'Erhalte diesen begrenzten Suchbefund und wiederhole ihn nur bei neuer Quelle. '
                    'Entwickle jetzt die begründete primitive Kandidatenregel und prüfe zuerst '
                    'Normalisierbarkeit, Originaltyp und ihre erste neue Annahme. ' + task)
        return GOAL_ANCHOR + task
    if proposal['proposed_status'] == 'no_new_information' or any(
        h['proposal']['new_information'].strip() == proposal['new_information'].strip() for h in history):
        return ('The previous route produced no distinct new information. Read the complete journal '
                'and choose a different concrete original equation, source dependency or countertest '
                'upstream of the missing nu_Vn/D_Vn input or identical marked IR caller. '
                'Name the changed question and run its decisive test; do not repeat a renamed kernel. ' + task)
    return task


def run(repo, out, requested_rounds=0, *, backend_factory=Backend, context_loader=prepare):
    repo, out = Path(repo).resolve(), Path(out).resolve()
    if out.is_relative_to(repo):
        raise ValueError('Run output must be outside the read-only canonical repository')
    out.mkdir(parents=True, exist_ok=True)
    lock = (out / 'controller.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    stopped = False
    def stop_handler(signum, frame):
        nonlocal stopped
        stopped = True
    old_handlers = {s: signal.signal(s, stop_handler) for s in (signal.SIGTERM, signal.SIGINT)}
    journal = out / 'research-rounds.jsonl'
    journal.touch(exist_ok=True)
    history = rows(journal)
    previous = json.loads((out / 'worker.json').read_text()) if (out / 'worker.json').exists() else {}
    totals = {'rounds': len(history), 'model_calls': 0, 'tool_events': 0,
              'completed_commands': 0, 'failed_commands': 0}
    # Persist real calls even if a process stopped between proposal and critic.
    for p in (out / 'rounds').glob('*/*/metrics.json'):
        m = json.loads(p.read_text())
        totals['model_calls'] += 1
        for key in ('tool_events', 'completed_commands', 'failed_commands'):
            totals[key] += m.get(key, 0)
    worker = {'kind': 'original_source_research', 'pid': os.getpid(), 'status': 'running',
        'started_utc': previous.get('started_utc', utc()), 'last_update_utc': utc(),
        'requested_rounds': requested_rounds, 'finished_rounds': len(history),
        'backend_identity': {'provider': 'openai_codex', 'model': 'gpt-6.1-sol',
                             'receipt_verified': bool(previous.get('backend_receipt'))},
        'backend_receipt': previous.get('backend_receipt'), 'trial_summary': totals,
        'physics_admission': ADMISSION}
    def persist():
        worker['last_update_utc'] = utc(); write_json(out / 'worker.json', worker)
    def update(call):
        worker['active_call'] = call
        worker['phase']['description'] = 'Originalquellen und Rechnungen im laufenden ' + call['role'] + ' Aufruf'
        persist()
    backend = backend_factory(lambda: stopped, update)
    task = (history[-1]['next_task'] if history and 'target_assessment' in history[-1]
            else GOAL_ANCHOR + DIRECT_TARGET_TASK + INITIAL_TASK)
    try:
        while not stopped and (requested_rounds == 0 or len(history) < requested_rounds):
            manifest = context_loader(repo, out)
            worker['knowledge_context'] = {k: v for k, v in manifest.items() if k != 'exact_source_files'}
            round_id = 'research-' + format(len(history) + 1, '06d')
            directory = out / 'rounds' / round_id; directory.mkdir(parents=True, exist_ok=True)
            state_path = directory / 'round-state.json'
            state = json.loads(state_path.read_text()) if state_path.exists() else {
                'round_id': round_id, 'question': task, 'started_utc': utc(),
                'context_version': manifest['context_version'], 'steering_version': 2}
            if 'steering_version' not in state and not any(k in state for k in ('proposer', 'critic')):
                # An old controller can have opened the next directory before
                # the upgrade signal. Preserve that question and all attempt
                # files, while applying the user's new task before any cached
                # completed role exists. Completed proposals are never replaced.
                state.setdefault('prior_questions', []).append({
                    'question': state['question'], 'context_version': state.get('context_version'),
                    'reason': 'User-authorized source-candidate steering repair', 'at': utc()})
                state.update(question=task, context_version=manifest['context_version'], steering_version=2)
                write_json(state_path, state)
            task = state['question']
            for role in ('proposer', 'critic'):
                if role in state:
                    continue
                if stopped:
                    break
                worker['active_question'] = task
                worker['phase'] = {'id': role, 'label': 'Herleitung' if role == 'proposer' else 'Unabhängige Gegenprüfung',
                    'description': task, 'current': totals['model_calls'],
                    'total': requested_rounds * 2 if requested_rounds else None, 'started_at': utc()}
                persist()
                attempt = len(list(directory.glob(role + '-*'))) + 1
                folder = directory / (role + '-' + format(attempt, '03d'))
                compact_manifest = {k: v for k, v in manifest.items()
                                    if k not in {'exact_source_files', 'scientific_readonly_roots'}}
                operator_notes = out / 'operator-notes.txt'
                if operator_notes.is_file():
                    compact_manifest['operator_notes_path'] = str(operator_notes)
                candidate = {'response': state['proposer'],
                    'call': {k: v for k, v in state['proposer_call'].items() if k != 'commands'},
                    'proposal_context_version': state['context_version'],
                    'context_changed_since_proposal': state['context_version'] != manifest['context_version']} if role == 'critic' else None
                totals['model_calls'] += 1
                persist()
                result, metrics = backend(role, folder, manifest,
                    prompt(role, compact_manifest, task, journal, candidate,
                           decision_history=decision_history(history)))
                state[role], state[role + '_call'] = result, metrics
                write_json(state_path, state)
                for key in ('tool_events', 'completed_commands', 'failed_commands'):
                    totals[key] += metrics.get(key, 0)
                worker['backend_receipt'] = metrics['receipt']
                worker['backend_identity']['receipt_verified'] = True
                persist()
            if stopped:
                break
            proposal, review = state['proposer'], state['critic']
            assessment = target_assessment(proposal, review, state['proposer_call'], state['critic_call'], manifest['python'])
            next_task = followup(proposal, review, history, assessment)
            row = {'round_id': round_id, 'question': task, 'started_utc': state['started_utc'],
                'completed_utc': utc(), 'proposal': proposal, 'review': review,
                'proposal_call': state['proposer_call'], 'review_call': state['critic_call'],
                'proposed_status': proposal['proposed_status'], 'critic_verdict': review['verdict'],
                'actual_result': {'observed_commands': state['proposer_call'].get('commands', []) + state['critic_call'].get('commands', []),
                    'scope': 'Host observed CLI command outputs and exits; scientific interpretation remains unadmitted'},
                'physics_admission': ADMISSION, 'target_assessment': assessment,
                'next_task': next_task}
            with journal.open('a') as stream:
                stream.write(json.dumps(row, ensure_ascii=False) + '\n'); stream.flush(); os.fsync(stream.fileno())
            history.append(row); task = next_task
            worker['finished_rounds'] = len(history); totals['rounds'] = len(history)
            worker['active_call'] = None; persist()
            if (out / 'stop-after-round').exists():
                (out / 'stop-after-round').unlink()
                stopped = True
        worker['status'] = 'stopped' if stopped else 'completed'
    except InterruptedError:
        worker['status'] = 'stopped'
    except Exception as exc:
        worker['status'] = 'error'; worker['error'] = {'type': type(exc).__name__, 'message': str(exc)}
    finally:
        backend.terminate()
        worker['phase'] = {'id': worker['status'], 'label': worker['status'],
            'description': worker.get('error', {}).get('message', task), 'current': len(history),
            'total': requested_rounds or None, 'started_at': utc()}
        persist()
        for s, handler in old_handlers.items():
            signal.signal(s, handler)
        lock.close()
    return 1 if worker['status'] == 'error' else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--rounds', type=int, default=0)
    args = parser.parse_args()
    if args.rounds < 0:
        parser.error('--rounds must be 0 (unbounded) or positive')
    sys.exit(run(args.repo, args.out, args.rounds))
