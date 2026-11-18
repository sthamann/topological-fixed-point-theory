#!/usr/bin/env python3
"""Session control of 17 November 2026: scalar return marker U = J^2 against
the Cartan marker, through the unchanged native consumer (arithmetic control).

J = X_gamma + X_-gamma is the Hermitian middle current; J^2 is the projector
on the states that gamma moves. U = J^2 is scalar and nonnegative on the
fibre of gamma. Compared: normalized U, Y, UY and the connected mixed
response for U = J^2 and for U = hypercharge Cartan (recorded contract).
Also: average over the 24 real weak directions of the compressed exchange.
Argument: repository root. Prints JSON, writes nothing.
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

matter = [s for s in product((1, -1), repeat=8) if s.count(-1) % 2 == 0]
index = {w: n for n, w in enumerate(matter)}
n = len(matter)
Y = (Fr(-1, 3), Fr(-1, 3), Fr(-1, 3), Fr(1, 2), Fr(1, 2), 0, 0, 0)
fam = lambda r: tuple(r[i] for i in range(5, 8))
fam4 = sorted({fam(r) for r in matter if r[:5].count(-1) % 2 == 0})
def moved(w, i, s, j, t):
    if w[i] == -s and w[j] == -t:
        v = list(w); v[i], v[j] = s, t
        return tuple(v)
cplx = lambda m: np.stack([np.asarray(m, float), np.zeros((n, n))], axis=-1).tolist()
eye = cplx(np.eye(n))
native = base/"domain/native-diagnostic-v1.json"
dom_id = json.loads(native.read_text()).get("domain_id", "native-diagnostic-v1")
sha = hashlib.sha256(native.read_bytes()).hexdigest()
word = repo/"experiments/theory-contracts/origin-diagnosis-new-hypotheses-20261116/check_word.py"
control = "arithmetic control: scalar return marker on the E8 word, not an original binding"

def consumer(keep, forms):
    sel = np.diag([float(fam(a) in keep or tuple(-x for x in fam(a)) in keep) for a in matter])
    bundle = {"schema": adapter.BUNDLE, "input_role": "arithmetic_control_only",
              "representation": "original_exact_blocked_endpoint_source",
              "matrix_convention": "sesquilinear_forms_in_original_word_basis",
              "ensemble_scope": "whole_original_ensemble", "native_inputs_sha256": sha, "domain_id": dom_id,
              "source_refs": [{"id": "E8WORD20261116/check_word.py", "path": str(word),
                               "sha256": hashlib.sha256(word.read_bytes()).hexdigest()}],
              "component_bindings": {k: control for k in adapter.REQUIRED_COMPONENTS},
              "scalar_factor": [1, 0], "reference_scalar": [1, 0],
              "retained_words": [list(a) for a in matter],
              "time": {"id": "control-placeholder", "kind": "original_euclidean_ordered", "T": 1.0, "t_U": 0.25, "t_Y": 0.5},
              "gram_form": eye, "boundary_pairing_form": eye, "physical_selector_form": cplx(sel),
              "ordered_source_forms": {k: cplx(v) for k, v in forms.items()}}
    r = adapter.evaluate_source({}, {"native_inputs_sha256": sha, "domain_id": dom_id},
                                {"whole_source_bundle": bundle, "_sha256": sha, "domain_id": dom_id})
    assert r["status"] == "evaluated_arithmetic_control_only", r
    return tuple(str(Fr(r["normalized"][k][0]).limit_denominator(10000)) for k in ("U", "Y", "UY", "connected_UY"))

rows = {}
for omit in [None]+fam4:
    keep = [f for f in fam4 if f != omit]
    for i in range(5):
        for j in (5, 6, 7):
            for s, t in product((1, -1), repeat=2):
                X = np.zeros((n, n))
                for w in matter:
                    v = moved(w, i, s, j, t)
                    if v: X[index[v], index[w]] = 1.0          # signs cancel in every mark below
                J2 = np.diag(((X+X.T) @ (X+X.T)).diagonal())     # projector on moved states
                ret = X.T @ X                                     # X_-gamma X_gamma
                forms = {"Z0": np.eye(n), "ZU": J2, "ZY": ret, "ZUY": J2 @ ret @ J2}
                key = ("128" if omit is None else "96", "weak" if i > 2 else "colour")+consumer(keep, forms)
                rows[key] = rows.get(key, 0)+1
print(json.dumps({"scalar_marker[space,current,U,Y,UY,connected,count]": [list(k)+[v] for k, v in sorted(rows.items())]}, indent=1))
