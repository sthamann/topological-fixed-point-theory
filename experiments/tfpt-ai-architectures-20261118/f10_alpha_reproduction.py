"""F10 - constants: reproduction group only.

Runs the existing TFPT alpha computation (verification/v3_em_alpha.py, unchanged) and compares
its positive root with CODATA 2022 and with the three most precise individual determinations.
All four reference values were known before this run, so this is a consistency check of the
evaluation, NOT a prediction test.  The reference numbers are typed in from the published
values and should be re-checked against the sources before any use beyond this folder.
"""
import json, re, subprocess, sys, os
root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
out = subprocess.run([sys.executable, "verification/v3_em_alpha.py"], cwd=root, capture_output=True, text=True).stdout
tfpt = float(re.search(r"alpha\^-1 = ([0-9.]+)", out).group(1))
refs = {"CODATA 2022 adjustment": (137.035999177, 21e-9),
        "Rb recoil (Morel et al. 2020)": (137.035999206, 11e-9),
        "electron g-2 (Fan et al. 2023)": (137.035999166, 15e-9),
        "Cs recoil (Parker et al. 2018)": (137.035999046, 27e-9)}
res = {"tfpt_alpha_inverse": tfpt, "script_checks": re.findall(r"(\d+ passed, \d+ failed)", out)}
for k, (v, s) in refs.items():
    res[k] = dict(value=v, sigma=s, tfpt_minus_value=round(tfpt - v, 12), pull_sigma=round((tfpt - v) / s, 2))
res["spread_between_measurements"] = "Rb and Cs recoil differ by %.1f sigma from each other" % (
    (137.035999206 - 137.035999046) / (11e-9 ** 2 + 27e-9 ** 2) ** 0.5)
res["status"] = ("Reproduction only.  A prospective test needs a frozen prediction vector with joint covariance "
                 "(ledger PRED.JOINTLIKELIHOOD.01) and target values released after the freeze; neither exists here.")
json.dump(res, open("out/f10_alpha_reproduction.json", "w"), indent=1)
print(json.dumps(res, indent=1))
