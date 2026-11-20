#!/usr/bin/env python3
"""Single entrypoint. Emits machine-readable results after exact assertions."""
import json
from pathlib import Path
from pachner_minimal import test as pachner
from five_qubit_qec import test as qec
from pentagon_coherence import test as pentagon
from bridge_isometry import test as bridge
from modular_code_e8 import test as modular
from flat_sign_baseline import test as flat_null
from syndrome_geometry import test as syndrome
from information_time_nogos import test as nogos


def main():
    results={'pachner':pachner(), 'pentagon_coherence':pentagon(),
             'five_qubit_qec':qec(), 'conditional_tensor_bridge':bridge(),
             'flat_sign_null':flat_null(), 'syndrome_geometry_e8_bridge':syndrome(),
             'modular_code_e8_control':modular(),
             'information_time_no_gos':nogos()}
    dst=Path(__file__).with_name('RESULTS.json')
    dst.write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
    for name, r in results.items():
        print('===',name,'===')
        print(json.dumps(r,ensure_ascii=False,indent=2))
    print('PASS: all exact assertions; wrote',dst)


if __name__=='__main__':main()
