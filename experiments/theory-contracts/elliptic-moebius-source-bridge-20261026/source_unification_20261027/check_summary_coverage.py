"""Read-only continuity audit of the four user-supplied Atlas summaries.

Verifies recorded source identity and extent, not mathematical correctness.
No calculation or old long-running experiment is rerun.
"""
import hashlib
import json
from pathlib import Path

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
BASE = REPO / 'experiments/theory-contracts'
OUT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


program = read(REPO / 'semantic_explorer/research/program.json')
records = {
    f['id']: (p['id'], f)
    for p in program['problems'] for f in p.get('findings', [])
}
e = read(BASE / 'elliptic-moebius-source-bridge-20261026/SESSION-COVERAGE-20261026.json')
a = read(BASE / 'native-arc-charge-source-interface-20261026/SESSION-ATLAS-MANIFEST.json')
d = read(Path('/Users/stefanhamann/Documents/Codex/2026-10-26/wrkl-r-mir-die-wrgwbnisse-der/outputs/Atlas-Abgleich-2026-10-26.json'))
k = read(BASE / 'original-admissible-source-closure-20261025/kohonen_session_coverage/coverage.json')
groups = {
    'early_geometry_and_spin_theta': [f['id'] for g in e['investigation_navigation'] for f in g['findings']],
    'algorithmic_and_joint_boundary_source': a['required_atlas_finding_ids'],
    'w8_single_seam_mode_and_d5_response': [f['finding_id'] for g in d['coverage'] for f in g['entries']],
    'recursive_machine_and_original_source': [f['id'] for g in k['groups'] for f in g['findings']],
}
errors = []
cache = {}


def inspect(s, owner):
    path = s.get('path', s.get('archived_path'))
    target = REPO / path
    if path not in cache:
        if target.is_file():
            raw = target.read_bytes()
            cache[path] = {'sha256': hashlib.sha256(raw).hexdigest(),
                           'bytes': len(raw), 'lines': len(raw.splitlines())}
        else:
            cache[path] = None
    actual = cache[path]
    if actual is None:
        errors.append({'owner': owner, 'path': path, 'error': 'missing'})
    elif actual['sha256'] != s['sha256']:
        errors.append({'owner': owner, 'path': path, 'error': 'hash_mismatch'})
    elif 'line_start' in s and not 1 <= s['line_start'] <= s['line_end'] <= actual['lines']:
        errors.append({'owner': owner, 'path': path, 'error': 'invalid_extent'})


coverage = []
for name, ids in groups.items():
    missing = [fid for fid in ids if fid not in records]
    sources = [s for fid in ids if fid in records for s in records[fid][1]['sources']]
    for fid in ids:
        if fid in records:
            for s in records[fid][1]['sources']:
                inspect(s, fid)
    coverage.append({'summary': name, 'count': len(ids), 'unique_findings': len(set(ids)),
                     'missing_ids': missing, 'source_references': len(sources),
                     'source_files': len({s['path'] for s in sources}), 'finding_ids': ids})
    errors.extend({'owner': name, 'error': 'missing_finding', 'id': fid} for fid in missing)

for s in a['figure_and_source_assets']:
    inspect(s, 'algorithmic_assets')
for s in k['reports']:
    inspect(s, 'machine_archived_report')

new_ids = [fid for p in program['problems'] if p['id'] in {
    'problem:marked-seam-monodromy-and-lattice-lifts',
    'problem:atlas-clock-defect-and-operational-reconstruction'}
    for fid in [f['id'] for f in p['findings']]
    if fid not in set(groups['early_geometry_and_spin_theta'])]
for fid in new_ids:
    for s in records[fid][1]['sources']:
        inspect(s, fid)

result = {'date': '2026-10-27', 'atlas_revision': program['revision'],
          'scope': 'Four supplied historical coverage claims checked against current canonical Atlas. Source identities and ranges only; not a new proof replay.',
          'current_total_findings': len(records), 'groups': coverage,
          'unique_findings_in_four_summaries': len(set().union(*map(set, groups.values()))),
          'early_geometry_and_spin_theta_current_count': len(groups['early_geometry_and_spin_theta']) + len(new_ids),
          'findings_added_since_197_snapshot': new_ids,
          'files_checked': len(cache),
          'explanatory_png_count': sum(s['path'].endswith('.png') for s in a['figure_and_source_assets']),
          'archived_machine_reports': len(k['reports']),
          'errors': errors, 'proofs_replayed': False, 'physical_gates_promoted': False}
(OUT / 'SUMMARY-COVERAGE-CHECK.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
print(json.dumps({key: value for key, value in result.items() if key not in {'groups', 'findings_added_since_197_snapshot'}}, ensure_ascii=False, indent=2))
if errors:
    raise SystemExit(1)
