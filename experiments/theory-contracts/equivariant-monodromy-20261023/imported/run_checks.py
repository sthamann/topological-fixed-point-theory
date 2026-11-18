"""Run each exact checker in a fresh workspace; compare all enumerated results.
Usage: python run_checks.py [--optimized]
Dependencies: Python 3, numpy, sympy. No network calls or repository access.
"""
from __future__ import annotations
import argparse,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path

def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--optimized',action='store_true',help='also execute all checks with Python -OO')
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    expected=json.loads((here/'expected_results.json').read_text())
    scripts=['01_monodromy.py','02_clocks.py','03_bundles.py','04_integral.py','05_cocycle.py']
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    outputs=[]
    modes=[False,True] if args.optimized else [False]
    for optimized in modes:
        with tempfile.TemporaryDirectory(prefix='tfpt_math_') as t:
            target=Path(t)
            for name in scripts:shutil.copyfile(here/name,target/name)
            for name in scripts:
                cmd=[sys.executable]+(['-OO'] if optimized else [])+[str(target/name)]
                result=subprocess.run(cmd,cwd=target,env=env,check=True,text=True,capture_output=True)
                print(f'PASS {name}; optimized={optimized}')
            actual=json.loads((target/'work/results.json').read_text())
            if actual != expected:
                raise RuntimeError('Computed classification differs from expected_results.json')
            cocycle=(target/'work/cocycle_results.json').read_bytes()
            outputs.append((actual,cocycle))
            (here/'verified_results.json').write_text(json.dumps(actual,indent=2)+'\n')
            (here/'verified_cocycle.json').write_bytes(cocycle)
    if any(o!=outputs[0] for o in outputs[1:]):
        raise RuntimeError('Normal and optimized execution differ')
    summary={
      'group_order':expected['group_order'],
      'ordered_reflection_quadruples':expected['quadruple_count'],
      'conjugacy_orbits':len(expected['orbits']),
      'hurwitz_orbits':len(expected['hurwitz_orbits']),
      'strict_quarterclock_orbits':sum(sum(c['order']==4 for c in o.get('clock_classes',[])) for o in expected['orbits']),
      'strict_clocked_pairs':sum(o['size']*o.get('clock_candidates_by_order',{}).get('4',0) for o in expected['orbits']),
      'modes':['normal','optimized'] if args.optimized else ['normal'],
      'all_passed':True
    }
    (here/'verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
