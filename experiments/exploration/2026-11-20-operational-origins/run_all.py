#!/usr/bin/env python3
"""Reproduce every shipped exact output. Python >=3.10 and SymPy required."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--originlab-root',type=Path)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    scripts=[('qec_origin','qec_output'),('syndrome_spread','syndrome_spread_output'),
             ('causal_origin','causal_output'),('qec_bridge_audit','qec_bridge_output'),
             ('qec_lift_check','qec_lift_output'),
             ('sector_origin','sector_origin_output'),('arithmetic_baseline','arithmetic_baseline_output'),
             ('reversible_control','reversible_control_output'),('originlab_adapter','originlab_adapter_output')]
    for script,output in scripts:
        command=[sys.executable,str(here/(script+'.py'))]
        if script == 'originlab_adapter' and args.originlab_root:
            command += ['--originlab-root',str(args.originlab_root)]
        result=subprocess.run(command,capture_output=True,text=True,check=True,cwd=here)
        payload=json.loads(result.stdout)
        (here/(output+'.json')).write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        print(script+': PASS',flush=True)


if __name__ == '__main__':
    main()
