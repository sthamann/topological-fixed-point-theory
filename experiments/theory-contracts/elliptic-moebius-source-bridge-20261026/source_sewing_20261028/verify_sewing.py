"""Replay the narrow certificates, preserving optimized runtime checks."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

HERE = Path(__file__).resolve().parent
CASES = [
    ('current/check_native_gamma_current.py', 'current/native-gamma-current-results.json'),
    ('clock/check_native_six_ovoid_clock_bundle.py', 'clock/native-six-ovoid-clock-bundle-results.json'),
    ('clock/check_e8_common_clock_carrier.py', 'clock/e8-common-clock-results.json'),
    ('polarization/check_spin5_quartic_polarization.py', 'polarization/spin5-quartic-polarization-results.json'),
    ('sewing/check_ward_charged_module_selection.py', 'sewing/ward-charged-module-selection-results.json'),
    ('polarization/check_charge_typed_hodge_bridge.py', 'polarization/charge-typed-hodge-bridge-results.json'),
]

def replay(case):
    checker, result = case
    observed = []
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    for flags in ([], ['-OO']):
        process = subprocess.run([sys.executable, *flags, str(HERE/checker)],
                                 cwd=HERE, env=env, capture_output=True,
                                 timeout=480, check=False)
        if process.returncode:
            raise RuntimeError(checker + ': ' + process.stderr.decode()[-3000:])
        observed.append(((HERE/result).read_bytes(), process.stdout))
    if observed[0] != observed[1]:
        raise RuntimeError('Normal and optimized results differ: '+checker)
    return {'checker':checker, 'result':result, 'normal_and_optimized':'PASS',
            'result_sha256':hashlib.sha256(observed[1][0]).hexdigest(),
            'stdout_sha256':hashlib.sha256(observed[1][1]).hexdigest()}

def main():
    # Charge pins the current result. Run it only after that result is final.
    with ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(replay, CASES[:-1]))
    rows.append(replay(CASES[-1]))
    result = {'status':'PASS', 'scope':'Six exact scoped certificates only; no full physical TFPT proof.',
              'python':sys.version.split()[0], 'certificates':rows}
    (HERE/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
