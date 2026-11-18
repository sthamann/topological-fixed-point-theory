#!/usr/bin/env python3
"""Refresh/check the single Atlas without LLM calls or research execution."""
import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def refresh_papers(check=False):
    """Preserve the explicit local scope and vectors without any model request."""
    lock = HERE / 'data/.papers.lock'
    if lock.exists() and not check:
        try:
            pid = int(lock.read_text())
            if pid <= 0:
                raise ValueError('invalid paper importer PID')
            os.kill(pid, 0)
        except ProcessLookupError:
            print(json.dumps({'paper_refresh': 'blocked_stale_lock', 'pid': pid,
                              'notice': 'Inspect the stale importer lock; no paper data was changed.'}))
            return 2
        except (ValueError, OSError) as error:
            print(json.dumps({'paper_refresh': 'blocked_lock_error', 'error': str(error)}))
            return 2
        print(json.dumps({'paper_refresh': 'active_import_preserved', 'pid': pid,
                          'notice': 'The active importer owns publication; offline refresh did not overwrite it.'}))
        return 0
    return subprocess.run([sys.executable, str(HERE / 'papers.py'),
                           '--check' if check else '--offline', '--refresh-registry'], cwd=REPO).returncode


def check_ledger():
    with (REPO / 'verification/status_ledger.csv').open(encoding='utf-8-sig', newline='') as handle:
        rows = {row['claim_id']: row for row in csv.DictReader(handle)}
    graph = json.loads((HERE / 'data/knowledge.json').read_text())
    recorded = {n['id'].removeprefix('claim:'): n.get('attrs', {}).get('ledger')
                for n in graph['nodes'] if n['type'] == 'claim'}
    changed = sorted(cid for cid in rows.keys() | recorded.keys() if rows.get(cid) != recorded.get(cid))
    with (REPO / 'verification/script_registry.csv').open(encoding='utf-8-sig', newline='') as handle:
        scripts = list(csv.DictReader(handle))
    script_nodes = {n['id']: n['attrs'] for n in graph['nodes'] if n['type'] == 'script'}
    linked = {edge['dst'] for edge in graph['edges'] if edge['type'] == 'implemented_check'}
    script_issues = []
    for row in scripts:
        script_id = 'script:' + row['script']
        path = f"verification/{row['script']}.py"
        attrs = script_nodes.get(script_id, {})
        source = REPO / path
        if attrs.get('path') != path or not source.is_file():
            script_issues.append({'id': script_id, 'issue': 'missing_script_or_source'})
        elif attrs.get('sha256') != sha256(source.read_bytes()).hexdigest():
            script_issues.append({'id': script_id, 'issue': 'stale_source'})
        if attrs.get('cluster') != row['cluster']:
            script_issues.append({'id': script_id, 'issue': 'registry_cluster_mismatch'})
        if script_id not in linked:
            script_issues.append({'id': script_id, 'issue': 'no_ledger_link'})
    print(json.dumps({'atlas_ledger_current': not changed, 'changed_claims': changed,
                      'registered_scripts_checked': len(scripts), 'script_issues': script_issues,
                      'coverage': graph.get('meta', {}).get('coverage'),
                      'notice': 'Freshness is not a scientific or semantic correctness verdict.'}))
    return not changed and not script_issues


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='read-only source and ledger freshness check')
    args = parser.parse_args()
    if args.check:
        inventory = subprocess.run([sys.executable, str(HERE / 'source_inventory.py'), '--check'], cwd=REPO)
        papers = refresh_papers(check=True)
        return 0 if check_ledger() and inventory.returncode == 0 and papers == 0 else 1
    # The extractor owns its process lock. If an extraction is active, stop here
    # instead of racing its publication or running a second model pass.
    for script, flags in [('extract.py', ['--offline', '--preserve-cache']),
                          ('build_corpus.py', []), ('build_context.py', []),
                          ('audit_extractions.py', [])]:
        result = subprocess.run([sys.executable, str(HERE / script), *flags], cwd=REPO)
        if result.returncode:
            return result.returncode
    return refresh_papers()


if __name__ == '__main__':
    raise SystemExit(main())
