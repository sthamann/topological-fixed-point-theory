#!/usr/bin/env python3
"""Session control of 17 November 2026: the unchanged native whole-source
consumer (original_source_adapter.evaluate_source) fed with the E8 word of
E8WORD20261116 as an ARITHMETIC CONTROL.

The bundle carries input_role = arithmetic_control_only. It supplies no
original source: the seven component bindings name the control, the time
values are placeholders, cap = identity, scalar = reference = 1. The run
shows what the consumer's own contraction Tr(cap P q P) returns for this
word. Argument: repository root. Prints JSON, writes nothing.
"""
import hashlib, json, sys
from fractions import Fraction as Fr
from itertools import combinations, product
from pathlib import Path
import numpy as np

repo = Path(sys.argv[1])
base = repo/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
import original_source_adapter as adapter

R = []
for i, j in combinations(range(8), 2):
    for a, b in product((2, -2), repeat=2):
        v = [0]*8; v[i], v[j] = a, b; R.append(tuple(v))
matter = [s for s in product((1, -1), repeat=8) if s.count(-1) % 2 == 0]
RS = set(R) | set(matter)
bif = [r for r in R if sum(1 for i in range(5) if r[i]) == 1]
add = lambda a, b: tuple(x+y for x, y in zip(a, b))
Y = (Fr(-1, 3), Fr(-1, 3), Fr(-1, 3), Fr(1, 2), Fr(1, 2), 0, 0, 0)
dot = lambda h, r: sum(Fr(x)*y for x, y in zip(h, r))/2
fam = lambda r: tuple(r[i] for i in range(5, 8))
fam4 = sorted({fam(r) for r in matter if r[:5].count(-1) % 2 == 0})
n = len(matter)
cplx = lambda d: np.stack([np.diag(d).astype(float), np.zeros((n, n))], axis=-1).tolist()
eye = cplx(np.ones(n))
word = repo/"experiments/theory-contracts/origin-diagnosis-new-hypotheses-20261116/check_word.py"
native = base/"domain/native-diagnostic-v1.json"
dom_id = json.loads(native.read_text()).get("domain_id", "native-diagnostic-v1")
sha = hashlib.sha256(native.read_bytes()).hexdigest()
control = "arithmetic control: E8 word of E8WORD20261116, not an original binding"

def run(keep, g):
    steps = np.array([float(add(a, g) in RS) for a in matter])
    zu = np.array([float(dot(Y, a)**2) for a in matter])
    zuy = np.array([float(dot(Y, a)*dot(Y, add(a, g))) if add(a, g) in RS else 0.0 for a in matter])
    bundle = {"schema": adapter.BUNDLE, "input_role": "arithmetic_control_only",
              "representation": "original_exact_blocked_endpoint_source",
              "matrix_convention": "sesquilinear_forms_in_original_word_basis",
              "ensemble_scope": "whole_original_ensemble",
              "native_inputs_sha256": sha, "domain_id": dom_id,
              "source_refs": [{"id": "E8WORD20261116/check_word.py", "path": str(word),
                               "sha256": hashlib.sha256(word.read_bytes()).hexdigest()}],
              "component_bindings": {k: control for k in adapter.REQUIRED_COMPONENTS},
              "scalar_factor": [1, 0], "reference_scalar": [1, 0],
              "retained_words": [list(a) for a in matter],
              "time": {"id": "control-placeholder", "kind": "original_euclidean_ordered",
                       "T": 1.0, "t_U": 0.25, "t_Y": 0.5},
              "gram_form": eye, "boundary_pairing_form": eye,
              "physical_selector_form": cplx(np.array([float(fam(a) in keep or tuple(-x for x in fam(a)) in keep)
                                                       for a in matter])),
              "ordered_source_forms": {"Z0": eye, "ZU": cplx(zu), "ZY": cplx(steps), "ZUY": cplx(zuy)}}
    return adapter.evaluate_source({"candidate_id": "e8-word-control"},
                                   {"native_inputs_sha256": sha, "domain_id": dom_id},
                                   {"whole_source_bundle": bundle, "_sha256": sha, "domain_id": dom_id})

frac = lambda z: str(Fr(z[0]).limit_denominator(10000))
table, statuses = {}, set()
for omit in [None]+fam4:
    keep = [f for f in fam4 if f != omit]
    for g in bif:
        r = run(keep, g)
        statuses.add((r["status"], r["physical_source_complete"]))
        nz = r["normalized"]
        kind = "weak" if abs(dot(Y, g)) == Fr(1, 2) else "colour"
        key = ("128" if omit is None else "96", kind, frac(r["Z0"]), frac(nz["U"]), frac(nz["Y"]),
               frac(nz["UY"]), frac(nz["connected_UY"]))
        table[key] = table.get(key, 0)+1
print(json.dumps({"consumer_status": sorted(map(list, statuses)),
                  "rows[space,current,Z0,U,Y,UY,connected]": [list(k)+[v] for k, v in sorted(table.items())],
                  "without_bundle": adapter.evaluate_source({}, {"native_inputs_sha256": sha, "domain_id": dom_id},
                                                            {"_sha256": sha, "domain_id": dom_id})["status"]}, indent=1))
