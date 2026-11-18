"""Recompute the bounded search against the original TFPT source."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

REPO = Path(os.environ.get('TFPT_REPO', '/Users/stefanhamann/Projekte/tfpt-theoryv4'))
sys.path.insert(0, str(REPO))
import sympy as sp
from tfpt_explorer.neutral_source_response import _native_neutral_branch
import source_filter_search as filters
import source_filter_joint_gate as ports
import source_joint_vertex_search as vertices


def run(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    native = _native_neutral_branch()
    Rstate = native['states']['R']
    Rnorm = (Rstate.T*native['gram']*Rstate)[0]
    assert Rnorm == 48 and sp.trace(native['R_matrix']) == 0
    # Compare the source-owned state norm with the independent Laurent expansion.
    for n in range(2, 21):
        assert filters.standard_primary_mode_norm(int(Rnorm), n) == 8*n*(n*n-1)
    results = {
        'filter': filters.main(),
        'ports': ports.main(),
        'vertices': vertices.run(),
    }
    for name, data in results.items():
        (outdir / (name + '.json')).write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
    source_paths = [
        'tfpt_explorer/neutral_source_response.py',
        'tfpt_explorer/source_charge_response.py',
        'tfpt_explorer/source_majorana_pairs.py',
        'tfpt_explorer/source_flavor_correlator.py',
    ]
    summary = {
        'contract': 'source-rule-search-20260930',
        'research_id': 'UR.SOURCE.RULE_SEARCH.01',
        'verdict': 'PARTIAL',
        'full_underlying_rule_found': False,
        'native_source_R_norm': str(Rnorm),
        'source_mode_norm': '8*n*(n*n-1), n>=2',
        'tested': {'shifted_monomial_filters': 37349, 'pure_power_filters': 221,
                   'rational_vertex_pairs': 8649},
        'exact_filter_hits': results['filter']['casimir_shift_class']['exact_hits'],
        'single_vertex_joint_hits': results['vertices']['single_vertex_F_aR_bRY']['bounded_rational_grid']['passing_candidates'],
        'source_files_sha256': {p: hashlib.sha256((REPO/p).read_bytes()).hexdigest() for p in source_paths},
        'scope': 'Bounded declared filter and neutral-vertex classes, fixed original source branch. Multiport consistency is conditional; no P1/P2 origin or physical EM derivation.',
        'all_executable_assertions_passed': True,
        'physical_gates_closed': [],
    }
    (outdir / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: summary[k] for k in ('research_id','verdict','full_underlying_rule_found','tested','all_executable_assertions_passed')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parent/'results')
    args = parser.parse_args()
    run(args.outdir)
