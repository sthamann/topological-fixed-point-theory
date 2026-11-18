"""Exact finite verification for TFPT mathematical topology, 2026-10-23.
Reconstructed from the pinned Gaussian-root formula; no upstream suite is invoked.
All matrix keys and group equalities use integer arithmetic.
"""

from pathlib import Path
WORK = Path(__file__).resolve().parent / 'work'
WORK.mkdir(exist_ok=True)
def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)

import pickle, json, itertools as it, time
from collections import deque, Counter
import numpy as np
import sympy as sp
start = time.time()
root = str(WORK) + '/'
refs, perms, pg, chosen, orbits, full, orbit_id, raykeys = pickle.load(open(root + 'state.pkl', 'rb'))
d = json.load(open(root + 'clocks.json'))

def key(x):
    return tuple(map(int, x[0].ravel())) + tuple(map(int, x[1].ravel()))

def mul(x, y):
    a, b = x
    c, e = y
    r = a @ c - b @ e
    s = a @ e + b @ c
    require(not np.any(r % 2) and (not np.any(s % 2)), 'Exact check, original line 14: not np.any(r % 2) and (not np.any(s % 2))')
    return (r // 2, s // 2)

def adj(x):
    return (x[0].T, -x[1].T)
I = (2 * np.eye(4, dtype=np.int64), np.zeros((4, 4), dtype=np.int64))
ki = key(I)

def cp(x, j):
    return key(mul(mul(x, refs[j]), adj(x)))
kr = [key(r) for r in refs]
x = sp.symbols('x0:4')
f = [sum((z ** 4 for z in x)), 6 * (x[0] ** 2 * x[1] ** 2 + x[2] ** 2 * x[3] ** 2), 6 * (x[0] ** 2 * x[2] ** 2 + x[1] ** 2 * x[3] ** 2), 6 * (x[0] ** 2 * x[3] ** 2 + x[1] ** 2 * x[2] ** 2), 24 * x[0] * x[1] * x[2] * x[3]]
mons = [(4, 0, 0, 0), (2, 2, 0, 0), (2, 0, 2, 0), (2, 0, 0, 2), (1, 1, 1, 1)]
mult = [1, 6, 6, 6, 24]
qs = {}
ws = {}
Vs = {}
for j in set(chosen) | {v for o in d['orbits'] for v in o['representative']}:
    M = sp.Matrix(refs[j][0].tolist()) / 2 + sp.I * sp.Matrix(refs[j][1].tolist()) / 2
    subs = {x[i]: sum((M[k, i] * x[k] for k in range(4))) for i in range(4)}
    cols = []
    for v in f:
        transformed = sp.Poly(sp.expand(v.subs(subs, simultaneous=True)), *x)
        col = [sp.cancel(transformed.coeff_monomial(mons[a]) / mult[a]) for a in range(5)]
        require(sp.Poly(transformed.as_expr() - sum((col[a] * f[a] for a in range(5))), *x).is_zero, 'Exact check, original line 36: sp.Poly(transformed.as_expr() - sum((col[a] * f[a] for a in range(5))), *x).is_zero')
        cols.append(sp.Matrix(col))
    Q = sp.Matrix.hstack(*cols)
    require(Q.det() == -1 and Q.trace() == 3, 'Exact check, original line 39: Q.det() == -1 and Q.trace() == 3')
    D = sp.diag(4, 12, 12, 12, 24)
    require(Q.conjugate().T * D * Q == D, 'Exact check, original line 41: Q.conjugate().T * D * Q == D')
    pairs = list(it.combinations(range(4), 2))
    W = sp.Matrix([[M[a, c] * M[b, e] - M[a, e] * M[b, c] for c, e in pairs] for a, b in pairs])
    qs[j] = Q
    ws[j] = W
    Vs[j] = M

def invariant_dim(mats, n):
    return n - sp.Matrix.vstack(*[m - sp.eye(n) for m in mats]).rank()
for o in d['orbits']:
    ts = o['representative']
    res = {}
    for name, bank, n in [('V4', Vs, 4), ('Q5', qs, 5), ('Lambda2V4', ws, 6)]:
        h0 = invariant_dim([bank[j] for j in ts], n)
        res[name] = {'h0': h0, 'h1': 2 * n + h0}
    o['cohomology'] = res
print('cohomology done', time.time() - start, flush=True)
I5 = sp.eye(5)
I6 = sp.eye(6)
eq = sp.Matrix.vstack(*[sp.kronecker_product(I5, ws[j]) - sp.kronecker_product(qs[j].T, I6) for j in chosen])
eqc = sp.Matrix.vstack(*[sp.kronecker_product(I5, ws[j]) - sp.kronecker_product((-qs[j]).T, I6) for j in chosen])
E6 = sp.Matrix.vstack(*[sp.kronecker_product(I6, ws[j]) - sp.kronecker_product(ws[j].T, I6) for j in chosen])
d['bundle_tests'] = {'hom_Q5_to_Lambda2V4_dimension': 30 - eq.rank(), 'hom_detQ5_to_Lambda2V4_dimension': 30 - eqc.rank(), 'commutant_Lambda2V4_dimension': 36 - E6.rank()}
seen = set()
bc = []
for i in range(len(d['orbits'])):
    if i in seen:
        continue
    orb = {i}
    todo = [i]
    seen.add(i)
    while todo:
        j = todo.pop()
        for k in d['orbits'][j]['hurwitz_ids']:
            if k not in seen:
                orb.add(k)
                seen.add(k)
                todo.append(k)
    bc.append(sorted(orb))
d['hurwitz_orbits'] = bc
open(root + 'results.json', 'w').write(json.dumps(d, indent=2))
print(json.dumps(d['bundle_tests']), flush=True)
print('FINAL elapsed', time.time() - start, flush=True)
