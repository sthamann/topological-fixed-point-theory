"""Durable complete context with read-only, scientifically scoped originals."""
from pathlib import Path
import hashlib
import json
import shutil
import sys

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE))
from knowledge_context import build_knowledge_context
from proposal_context import all_trial_history

TARGET_GUARDRAILS = {
    'finding:bounded-fixed-rank-spectral-action-does-not-normalize-noncompact-higgs-volume',
    'finding:rec20261103e-high-gap-retains-coupled-higgs-radial-tail',
    'finding:positive-spectral-reweighting-does-not-establish-original-source-rp',
    'finding:literal-archived-global-rg-kernel-does-not-preserve-finite-measures',
    'finding:whole-undecorated-transfer-does-not-determine-original-marked-field-source',
    'finding:actual-disjoint-local-weyls-acquire-spurious-commutator-under-bare-native-compression',
    'finding:solutionsearch20261106r3-output-search-no-physical-signal',
}


def target_guardrails(program):
    found = {}
    def visit(value):
        if isinstance(value, dict):
            if value.get('id') in TARGET_GUARDRAILS and 'sources' in value:
                found[value['id']] = value
            for v in value.values():
                visit(v)
        elif isinstance(value, list):
            for v in value:
                visit(v)
    visit(program)
    return {'findings': list(found.values()),
            'requested_ids_not_found': sorted(TARGET_GUARDRAILS - found.keys()),
            'scope': 'Selected source-target guardrails; the full Atlas remains binding, no claim of complete topical coverage'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temp.replace(path)


def prepare(repo, out):
    paths = [repo / 'semantic_explorer/research/program.json',
             repo / 'docs/research/fundamental-plan/progress.json',
             Path(__file__), Path(__file__).with_name('prompts.py')]
    reviews = [p for p in (BASE / 'CLAUDE-LIVE-RUN-REVIEW-20261106R5.txt',
                          BASE / 'CLAUDE-CANDIDATE-STRATEGY-20261106R5.txt') if p.is_file()]
    paths += reviews
    version = hashlib.sha256(''.join(digest(p) for p in paths).encode()).hexdigest()[:20]
    folder = out / 'context' / version
    manifest_path = folder / 'manifest.json'
    if manifest_path.exists():
        return json.loads(manifest_path.read_text())
    folder.mkdir(parents=True, exist_ok=True)
    pack = build_knowledge_context(repo)
    write_json(folder / 'full-knowledge-context.json', pack)
    guardrails = target_guardrails(json.loads(paths[0].read_text()))
    if guardrails['requested_ids_not_found']:
        raise FileNotFoundError('Required target-source findings missing: ' + ', '.join(guardrails['requested_ids_not_found']))
    write_json(folder / 'target-source-guardrails.json', guardrails)
    history = None
    history_origin = None
    for runs in dict.fromkeys((out.parent, BASE / 'runs')):
        if runs.is_dir() and any(runs.glob('*/results.jsonl')):
            candidate = all_trial_history(runs)
            if candidate['whole_programs']:
                history, history_origin = candidate, str(runs)
                break
    if history is None:
        # An installed package need not carry mutable runs. Reuse the complete
        # snapshot already archived by this same output, including all versions.
        archived = [(p, json.loads(p.read_text())) for p in
                    (out / 'context').glob('*/all-finite-trial-history.json')]
        archived = [(p, h) for p, h in archived if h.get('whole_programs')]
        if not archived:
            raise FileNotFoundError('No complete prior finite trial history: neither original runs nor an archived context snapshot is available')
        path, history = max(archived, key=lambda pair: len(pair[1]['whole_programs']))
        history_origin = str(path)
    write_json(folder / 'all-finite-trial-history.json', history)
    # Known scientific roots plus the actual trusted Q caller, never the user's
    # home, email, browser, credentials, website or unrelated projects.
    roots = [repo / x for x in ('documents', 'docs', 'verification', 'experiments',
             'verify', 'wolfram', 'notes', 'research-20261027-source-unification',
             'work/contribution20261106Q', '.venv') if (repo / x).exists()]
    for q in (out.parent.parent.parent / 'contribution20261106Q',
              BASE.parent / 'contribution20261106Q'):
        if q.is_dir():
            roots.append(q)
            break
    files = [p for p in repo.iterdir() if p.is_file() and p.suffix == '.tex']
    for row in pack['source_files'].values():
        p = Path(row['path']); p = p if p.is_absolute() else repo / p
        if p.is_file() and p.is_relative_to(repo) and not any(x.startswith('.') for x in p.relative_to(repo).parts):
            files.append(p.resolve())
    index = {'scientific_readonly_roots': list(map(str, roots)),
             'exact_source_files': sorted(set(map(str, files))),
             'reviewed_source_files': pack['source_files'], 'source_handles': pack['source_pool']}
    write_json(folder / 'source-index.json', index)
    for name, src in [('claude-correction.md', BASE / 'CLAUDE-REVIEW-AND-SEARCH-CORRECTION-20261106R3.md'),
                      ('original-integrand-review.txt', BASE / 'source_review/ORIGINAL-FINITE-INTEGRAND-SOURCE-REVIEW.txt'),
                      ('original-provenance.json', BASE / 'source_review/original-integrand-provenance.json')]:
        shutil.copyfile(src, folder / name)
    review_files = []
    for src in reviews:
        target = folder / src.name
        shutil.copyfile(src, target)
        review_files.append(str(target))
    manifest = {**pack['manifest'], 'context_version': version, 'repo': str(repo),
                'context_dir': str(folder), 'full_context_path': str(folder / 'full-knowledge-context.json'),
                'full_trial_history_path': str(folder / 'all-finite-trial-history.json'),
                'full_trial_history_origin': history_origin,
                'source_index_path': str(folder / 'source-index.json'),
                'scientific_readonly_roots': index['scientific_readonly_roots'],
                'exact_source_files': index['exact_source_files'],
                'claude_correction': str(folder / 'claude-correction.md'),
                'original_integrand_review': str(folder / 'original-integrand-review.txt'),
                'current_claude_reviews': review_files,
                'target_guardrails_path': str(folder / 'target-source-guardrails.json'),
                'fixed_target_guardrails': guardrails,
                'python': str(repo / 'experiments/tfpt-discovery/.venv/bin/python'),
                'global_context_embedded_in_full_context': True}
    write_json(manifest_path, manifest)
    return manifest
