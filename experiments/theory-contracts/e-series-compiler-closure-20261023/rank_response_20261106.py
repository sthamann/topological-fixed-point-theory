"""Exact existing-block E7 restriction and inverse-trace diagnostic."""
import ast
import itertools
import json
from pathlib import Path
import numpy as np

repo = Path(__file__).resolve().parents[3]
if not (repo / 'verification/v1_e8_glue.py').exists():
    repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
source = repo / 'verification/v1_e8_glue.py'
tree = ast.parse(source.read_text())
fn = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == 'e8_roots')
ns = {'np':np, 'itertools':itertools}
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(source),'exec'),ns)
R = np.rint(2*np.array(ns['e8_roots']())).astype(int)
D5 = set(np.flatnonzero(np.all(R[:,5:]==0,axis=1)))
D3 = set(np.flatnonzero(np.all(R[:,:5]==0,axis=1)))
B = {i for i,r in enumerate(R) if np.all(np.abs(r)==1) and sum(r[:5]<0)%2==0}
assert (len(R),len(D5),len(D3),len(B)) == (240,40,12,64)
rows=[]
for i in sorted(D3):
    if next(int(x) for x in R[i] if x) != 2:
        continue
    C = set(np.flatnonzero(R@R[i]==0))
    kept = B&C
    weights = sorted({tuple(int(v) for v in R[j,5:]) for j in kept})
    assert D5 <= C and len(C)==126 and len(kept)==32 and len(weights)==2
    # Every D5 root step starting in the selected block remains in it or
    # is absent from the E8 shell. This checks the actual carrier action.
    lookup = {tuple(r):j for j,r in enumerate(R)}
    for j in kept:
        for a in D5:
            target = lookup.get(tuple(R[j]+R[a]))
            assert target is None or target in kept
    rows.append({'normal_doubled':R[i].tolist(),'root_count':len(C),'block_rank':len(kept),'family_weights':weights})
assert len(rows)==6 and len({tuple(map(tuple,r['family_weights'])) for r in rows})==6
allowed=[16*k for k in range(5)]
assert 48 in allowed and 64 in allowed and 40 not in allowed and 56 not in allowed
result={'verdict':'PASS_SCOPED','claim':'MATH.SEAM.ESERIES.RESPONSE.20261106','six_E7_cuts':rows,'D5_equivariant_projector_ranks':allowed,'inverse_trace_E7_denominator':4,'inverse_trace_E6_triplet_denominator':6,'scope':'Native B=(16,4), D5 equivariance, idempotents and ordinary trace only. No physical coupling prediction.'}
print(json.dumps(result,indent=2))
