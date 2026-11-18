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

import itertools as it
import numpy as np
from collections import defaultdict, Counter, deque
import json, time
start = time.time()

def gkey(re, im):
    return tuple(map(int, re.ravel())) + tuple(map(int, im.ravel()))

def mul(x, y):
    a, b = x
    c, d = y
    return (a @ c - b @ d, a @ d + b @ c)

def canon_ray(z):
    a, b = z
    variants = [(a, b), (-b, a), (-a, -b), (b, -a)]
    return min((tuple(((int(x), int(y)) for x, y in zip(c, d))) for c, d in variants))
supports = [tuple(sorted((2 * i, 2 * i + 1, 2 * j, 2 * j + 1))) for i, j in it.combinations(range(4), 2)]
supports += [tuple((2 * k + d[k] for k in range(4))) for d in it.product((0, 1), repeat=4) if sum(d) % 2 == 0]
roots = []
for k in range(8):
    for sign in (-1, 1):
        x = [0] * 8
        x[k] = 2 * sign
        roots.append(x)
for sup in supports:
    for signs in it.product((-1, 1), repeat=4):
        x = [0] * 8
        for k, sg in zip(sup, signs):
            x[k] = sg
        roots.append(x)
raykeys = sorted({canon_ray((np.array(x[::2]), np.array(x[1::2]))) for x in roots})
refs = []
for key in raykeys:
    a = np.array([v[0] for v in key], dtype=np.int64)
    b = np.array([v[1] for v in key], dtype=np.int64)
    refs.append((2 * np.eye(4, dtype=np.int64) - np.outer(a, a) - np.outer(b, b), np.outer(a, b) - np.outer(b, a)))
require(len(roots) == 240 and len(refs) == 60, 'Exact check, original line 33: len(roots) == 240 and len(refs) == 60')
lookup = {gkey(*a): i for i, a in enumerate(refs)}
perms = []
for r in refs:
    perm = []
    for s in refs:
        re, im = mul(mul(r, s), r)
        require(not np.any(re % 4) and (not np.any(im % 4)), 'Exact check, original line 40: not np.any(re % 4) and (not np.any(im % 4))')
        perm.append(lookup[gkey(re // 4, im // 4)])
    perms.append(tuple(perm))

def pcompose(a, b):
    return tuple((a[j] for j in b))

def pclosure(gens):
    ident = tuple(range(60))
    seen = {ident}
    queue = deque([ident])
    while queue:
        x = queue.popleft()
        for g in gens:
            y = pcompose(g, x)
            if y not in seen:
                seen.add(y)
                queue.append(y)
    return seen
chosen = []
pg = {tuple(range(60))}
for i, p in enumerate(perms):
    if p not in pg:
        chosen.append(i)
        pg = pclosure([perms[j] for j in chosen])
        print('gen', chosen, 'order', len(pg), flush=True)
    if len(pg) == 11520:
        break
require(all((p in pg for p in perms)), 'Exact check, original line 58: all((p in pg for p in perms))')
print('permutation group', len(pg), 'elapsed', time.time() - start, flush=True)
buckets = defaultdict(list)
for i, a in enumerate(refs):
    for j, b in enumerate(refs):
        buckets[gkey(*mul(a, b))].append((i, j))
tuples = set()
for pairs in buckets.values():
    for i, j in pairs:
        for l, k in pairs:
            tuples.add((i, j, k, l))
print('quadruples', len(tuples), 'product buckets', len(buckets), 'bucket sizes', Counter(map(len, buckets.values())), flush=True)
unseen = tuples.copy()
orbits = []
orbit_id = {}
gens = [perms[i] for i in chosen]
while unseen:
    root = min(unseen)
    orbit = {root}
    queue = deque([root])
    unseen.remove(root)
    while queue:
        t = queue.popleft()
        for g in gens:
            u = tuple((g[i] for i in t))
            if u not in orbit:
                require(u in tuples, 'Exact check, original line 79: u in tuples')
                orbit.add(u)
                queue.append(u)
                unseen.remove(u)
    oid = len(orbits)
    for t in orbit:
        orbit_id[t] = oid
    orbits.append({'representative': root, 'size': len(orbit), 'stabilizer_projective': len(pg) // len(orbit), 'orbit': orbit})
print('orbits', len(orbits), 'sizes', Counter((o['size'] for o in orbits)), 'elapsed', time.time() - start, flush=True)
Id = (2 * np.eye(4, dtype=np.int64), np.zeros((4, 4), dtype=np.int64))

def mm(a, b):
    r, i = mul(a, b)
    require(not np.any(r % 2) and (not np.any(i % 2)), 'Exact check, original line 89: not np.any(r % 2) and (not np.any(i % 2))')
    return (r // 2, i // 2)
full = {gkey(*Id): Id}
queue = deque([Id])
while queue:
    x = queue.popleft()
    for j in chosen:
        y = mm(refs[j], x)
        k = gkey(*y)
        if k not in full:
            full[k] = y
            queue.append(y)
print('full group', len(full), 'elapsed', time.time() - start, flush=True)
require(len(full) == 46080, 'Exact check, original line 98: len(full) == 46080')
for o in orbits:
    sub = {gkey(*Id): Id}
    queue = deque([Id])
    gs = list(dict.fromkeys(o['representative']))
    while queue:
        x = queue.popleft()
        for j in gs:
            y = mm(refs[j], x)
            k = gkey(*y)
            if k not in sub:
                sub[k] = y
                queue.append(y)
    o['image_order'] = len(sub)
    o['stabilizer_full'] = 46080 // o['size']
    print({k: v for k, v in o.items() if k != 'orbit'}, flush=True)
for o in orbits:
    t = o['representative']
    o['rotation_orbit_id'] = orbit_id[t[1:] + t[:1]]

def hurwitz(t, j):
    u = list(t)
    a, b = u[j:j + 2]
    u[j] = perms[a][b]
    u[j + 1] = a
    return tuple(u)
for o in orbits:
    o['hurwitz_ids'] = [orbit_id[hurwitz(o['representative'], j)] for j in range(3)]
out = {'source_commit': '7e7fefa42f0e26cc4229c12b35543fe34ad56f18', 'raykeys': raykeys, 'generator_indices': chosen, 'group_order': len(full), 'projective_order': len(pg), 'quadruple_count': len(tuples), 'orbits': [{k: v for k, v in o.items() if k != 'orbit'} for o in orbits]}
open(str(WORK / 'interim.json'), 'w').write(json.dumps(out, indent=2))
import pickle
with open(str(WORK / 'state.pkl'), 'wb') as f:
    pickle.dump((refs, perms, pg, chosen, orbits, full, orbit_id, raykeys), f)
print('DONE', time.time() - start, flush=True)

# Second orbit count: Burnside, without using the conjugacy partition.
require(all(gkey(*r) in full for r in refs), 'All original reflections occur in full group')
II=np.eye(4,dtype=np.int64);ZZ=np.zeros((4,4),dtype=np.int64)
require(all(gkey(*u) in full for u in [(2*II,ZZ),(-2*II,ZZ),(ZZ,2*II),(ZZ,-2*II)]), 'Four scalar central elements')
fixed_subsets=Counter(tuple(i for i in range(60) if p[i]==i) for p in pg)
TT=np.array(sorted(tuples),dtype=np.int16)
burnside_sum=0
for subset,multiplicity in fixed_subsets.items():
    mask=np.zeros(60,dtype=bool);mask[list(subset)]=True
    burnside_sum+=multiplicity*int(np.all(mask[TT],axis=1).sum())
require(burnside_sum==195840, 'Exact Burnside sum')
require(burnside_sum//len(pg)==17, 'Independent count: 17 classes')
print('Independent Burnside orbit count:',burnside_sum//len(pg),flush=True)
