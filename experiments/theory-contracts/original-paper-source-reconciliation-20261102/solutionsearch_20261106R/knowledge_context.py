"""Lossless German reviewed Atlas packet for the existing proposal path.

This reads the owning records, not generated summaries.  Source hashes certify
file freshness only; they do not certify proofs.  No model calls or writes.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path


GLOBAL_CONTEXT = Path('/Users/stefanhamann/.codex/context/tfpt-gesamtkontext.md')
# These are communications/administration, not source-reviewed science.
# Their exact local source and hashes remain in the attachment manifest.
LOCAL_PLAN_FIELDS = frozenset({
    'counterparts', 'communication_protocol', 'user_directive',
    'latest_session_thread', 'entrypoints', 'canonical_plan',
})


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def german(value):
    """Remove the duplicate EN translation, retaining every German character.

    Nonlanguage metadata and EN-only records remain explicit, never silently
    dropped. English aliases and technical identifiers are not translations.
    """
    if isinstance(value, list):
        return [german(x) for x in value]
    if not isinstance(value, dict):
        return value
    if 'de' in value and set(value).issubset({'de', 'en'}):
        return german(value['de'])
    if set(value) == {'en'}:
        return {'translation_state': 'GERMAN_NOT_SUPPLIED', 'original_en': value['en']}
    return {k: german(v) for k, v in value.items()
            if not (k == 'en' and 'de' in value)}


def build_knowledge_context(repo, *, global_context_path=None, verify_sources=True):
    """Return all reviewed problems/premises/strategies/tasks plus the live plan.

    Sources are deduplicated by their complete original handle, including line
    ranges. File IO is separately deduplicated by resolved path. No last-N or
    character limits are applied anywhere in this function.
    """
    repo = Path(repo).resolve()
    atlas_path = repo / 'semantic_explorer/research/program.json'
    plan_path = repo / 'docs/research/fundamental-plan/progress.json'
    context_path = Path(global_context_path) if global_context_path else GLOBAL_CONTEXT
    program = json.loads(atlas_path.read_text())
    plan = json.loads(plan_path.read_text())
    pool, file_cache, source_files, handle_keys, file_keys = {}, {}, {}, {}, {}

    def intern(source):
        normalized = german(source)
        encoded = json.dumps(normalized, sort_keys=True, ensure_ascii=False,
                             separators=(',', ':')).encode()
        handle_key = hashlib.sha256(encoded).hexdigest()
        if handle_key not in handle_keys:
            sid = 's' + str(len(pool) + 1)
            handle_keys[handle_key] = sid
            original_path = Path(source['path'])
            path = (repo / original_path).resolve() if not original_path.is_absolute() else original_path.resolve()
            key = str(path)
            if key not in file_cache:
                if not verify_sources:
                    file_cache[key] = {'state': 'NOT_CHECKED'}
                else:
                    try:
                        file_cache[key] = {'state': 'READ', 'actual_sha256': file_sha256(path)}
                    except (OSError, ValueError) as exc:
                        file_cache[key] = {'state': 'UNREADABLE', 'error': type(exc).__name__ + ': ' + str(exc)}
            actual = file_cache[key]
            expected = source.get('sha256')
            state = ('NOT_CHECKED' if not verify_sources else
                     'MISSING_OR_UNREADABLE' if actual['state'] != 'READ' else
                     'HASH_NOT_SUPPLIED' if not expected else
                     'CURRENT' if actual['actual_sha256'] == expected else 'STALE')
            file_key = (source['path'], expected)
            if file_key not in file_keys:
                fid = 'f' + str(len(source_files) + 1)
                file_keys[file_key] = fid
                source_files[fid] = {'path': source['path'], 'sha256': expected,
                                     'freshness': state,
                                     **{k: v for k, v in actual.items()
                                        if k != 'state' and not (k == 'actual_sha256' and state == 'CURRENT')}}
            pool[sid] = {'file_ref': file_keys[file_key],
                         **{k: v for k, v in normalized.items() if k not in {'path', 'sha256'}}}
        return handle_keys[handle_key]

    def convert(value):
        if isinstance(value, list):
            return [convert(v) for v in value]
        if not isinstance(value, dict):
            return value
        if 'path' in value and 'sha256' in value:
            return {'source_handle': intern(value)}
        return {k: convert(v) for k, v in value.items()}

    # Preserve every field, including status, accepted/required premises,
    # conditional scope, source handles, attempts, tasks and exclusion changes.
    problems = convert(german(program['problems']))
    premises = convert(german(program['premises']))
    live_plan = convert(german({k: v for k, v in plan.items() if k not in LOCAL_PLAN_FIELDS}))
    attachments = [{
        'source_path': str(plan_path), 'json_pointer': '/' + key,
        'classification': 'unreviewed_communication_or_administration_NOT_scientific_evidence',
        'sha256': hashlib.sha256(json.dumps(plan[key], sort_keys=True, ensure_ascii=False,
                                          separators=(',', ':')).encode()).hexdigest(),
        'source_file_sha256': file_sha256(plan_path),
    } for key in sorted(LOCAL_PLAN_FIELDS) if key in plan]
    findings = [f for p in problems for f in p.get('findings', [])]
    tasks = [t for p in problems for t in p.get('tasks', [])]
    strategies = [s for p in problems for s in p.get('strategies', [])]
    completed = [{'problem_id': p['id'], 'task_id': t['id']}
                 for p in problems for t in p.get('tasks', [])
                 if t.get('status') in {'done', 'resolved', 'completed', 'closed'}]
    plan_completed = [{'parent_task_id': t['id'], 'subtask_id': s['id']}
                      for t in live_plan.get('tasks', [])
                      for s in t.get('completed_subtasks', [])]
    manifest = {
        'atlas_revision': program['revision'], 'atlas_sha256': file_sha256(atlas_path),
        'plan_revision': plan.get('revision'), 'plan_sha256': file_sha256(plan_path),
        'plan_atlas_revision': plan.get('current_atlas_revision'),
        'plan_atlas_revision_matches': plan.get('current_atlas_revision') == program['revision'],
        'global_context_path': str(context_path), 'global_context_sha256': file_sha256(context_path),
        'counts': {'problems': len(problems), 'findings': len(findings),
                   'premises': len(premises), 'strategies': len(strategies),
                   'problem_tasks': len(tasks), 'completed_problem_tasks': len(completed),
                   'plan_tasks': len(live_plan.get('tasks', [])),
                   'completed_plan_subtasks': len(plan_completed),
                   'unique_source_handles': len(pool), 'unique_source_files': len(file_cache)},
        'finding_status_counts': dict(Counter(f.get('status', 'UNSPECIFIED') for f in findings)),
        'source_freshness_counts': dict(Counter(source_files[s['file_ref']]['freshness'] for s in pool.values())),
        'language': 'de', 'arbitrary_truncation': False,
        'source_hash_check_is_proof': False,
    }
    return {
        'schema': 'tfpt.reviewed-model-context.v1', 'manifest': manifest,
        'interpretation': ('Alle Aussagen gelten genau mit ihren ursprünglichen Voraussetzungen und Geltungsbereichen. '
                           'Status, Drafts, historische Befunde und veraltete Quellen sind keine neuen Beweise. '
                           'CURRENT bezeichnet nur einen passenden Dateihash. Vorhandene Abschlüsse werden erhalten; '
                           'Ausschlüsse gelten nur für die benannte Klasse. Änderungen benötigen eine ausdrücklich geänderte Voraussetzung. '
                           'Der lebende Plan entscheidet über ältere Prioritätsangaben im Gesamtkontext.'),
        'global_context': context_path.read_text(),
        'premises': premises, 'problems': problems, 'current_plan': live_plan,
        'completed_problem_tasks': completed,
        'completed_plan_subtasks': plan_completed,
        'source_pool': pool, 'source_files': source_files,
        'local_attachments': attachments,
    }


def evidence(repo):
    """Drop-in name for the campaign's existing evidence loader."""
    return build_knowledge_context(repo)


def expand_source_handles(value, packet):
    """Recover original source dictionaries exactly, for consumers and controls."""
    if isinstance(value, list):
        return [expand_source_handles(x, packet) for x in value]
    if not isinstance(value, dict):
        return value
    if set(value) == {'source_handle'}:
        source = packet['source_pool'][value['source_handle']]
        file = packet['source_files'][source['file_ref']]
        return {'path': file['path'], 'sha256': file['sha256'],
                **{k: v for k, v in source.items() if k != 'file_ref'}}
    return {k: expand_source_handles(v, packet) for k, v in value.items()}
