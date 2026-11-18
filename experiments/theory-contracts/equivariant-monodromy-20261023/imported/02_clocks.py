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
d = json.load(open(root + 'interim.json'))

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
stable = [i for i, o in enumerate(d['orbits']) if i == o['rotation_orbit_id']]
print('stable tuple orbits', stable, flush=True)

def closure(gens):
    seen = {ki: I}
    queue = deque([I])
    while queue:
        x = queue.popleft()
        for g in gens:
            y = mul(g, x)
            k = key(y)
            if k not in seen:
                seen[k] = y
                queue.append(y)
    return seen
for oid in stable:
    o = d['orbits'][oid]
    t = o['representative']
    constraints = list(dict.fromkeys(zip(t, t[1:] + t[:1])))
    cent = []
    candidates = []
    byorder = Counter()
    for x in full.values():
        if all((cp(x, j) == kr[j] for j in set(t))):
            cent.append(x)
        if not all((cp(x, j) == kr[k] for j, k in constraints)):
            continue
        x2 = mul(x, x)
        x4 = mul(x2, x2)
        if key(x4) != ki:
            continue
        order = 1 if key(x) == ki else 2 if key(x2) == ki else 4
        byorder[order] += 1
        candidates.append((x, order))
    print('case', oid, 'centralizer', len(cent), 'clocks', byorder, flush=True)
    gen = []
    cg = {ki: I}
    for x in cent:
        if key(x) not in cg:
            gen.append(x)
            cg = closure(gen)
        if len(cg) == len(cent):
            break
    require(len(cg) == len(cent), 'Exact check, original line 53: len(cg) == len(cent)')
    remaining = {key(x): (x, ordr) for x, ordr in candidates}
    clockorbits = []
    while remaining:
        k = min(remaining)
        x, order = remaining.pop(k)
        orb = {k}
        todo = [x]
        while todo:
            y = todo.pop()
            for g in gen:
                z = mul(mul(g, y), adj(g))
                kz = key(z)
                if kz not in orb:
                    orb.add(kz)
                    todo.append(z)
                    require(kz in remaining, 'Exact check, original line 64: kz in remaining')
                    del remaining[kz]
        M = sp.Matrix(x[0].tolist()) / 2 + sp.I * sp.Matrix(x[1].tolist()) / 2
        clockorbits.append({'order': order, 'fiber_orbit_size': len(orb), 'automorphism_order': len(cent) // len(orb), 'matrix_real_twice': x[0].tolist(), 'matrix_imag_twice': x[1].tolist(), 'characteristic_polynomial': str(M.charpoly().as_expr())})
    o['clock_candidates_by_order'] = {str(k): v for k, v in byorder.items()}
    o['clock_classes'] = clockorbits
    open(root + 'clocks.json', 'w').write(json.dumps(d, indent=2))
    print('clock classes', len(clockorbits), Counter((z['order'] for z in clockorbits)), 'elapsed', time.time() - start, flush=True)
print('DONE CLOCKS', time.time() - start)
