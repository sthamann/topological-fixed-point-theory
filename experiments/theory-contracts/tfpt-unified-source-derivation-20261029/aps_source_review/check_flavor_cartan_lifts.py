"""Exact fixed D5+D3 Cartan comparison of Claude19; no physical identification."""
from itertools import combinations, permutations, product
from collections import Counter
import json
from pathlib import Path


def roots():
    result = set()
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            v = [0] * 8
            v[i], v[j] = a, b
            result.add(tuple(v))
    result.update(v for v in product((-1, 1), repeat=8) if v.count(-1) % 2 == 0)
    return result


def transform(v, permutation, signs, square=False):
    first = tuple(signs[i] * v[permutation[i]] for i in range(5))
    last = (-v[5], -v[6], v[7]) if square else (-v[6], v[5], v[7])
    return first + last


root_set = roots()
w64 = {v for v in root_set if all(abs(t) == 1 for t in v)
       and v[:5].count(-1) % 2 == 1 and v[5:].count(-1) % 2 == 1}
w64_conjugate = {tuple(-t for t in v) for v in w64}
charge = (-2, -2, -2, 3, 3, 0, 0, 0)
counts = Counter()
sm_signs = Counter()
for permutation in permutations(range(5)):
    for signs in product((-1, 1), repeat=5):
        counts['signed_D5_candidates'] += 1
        image = {transform(v, permutation, signs) for v in root_set}
        preserves = image == root_set
        odd = signs.count(-1) % 2 == 1
        if preserves != odd:
            raise RuntimeError('Root test disagrees with sign parity theorem')
        if not preserves:
            continue
        counts['E8_lifts'] += 1
        if {transform(v, permutation, signs) for v in w64} != w64_conjugate:
            raise RuntimeError('Fixed W64 not exchanged with conjugate')
        counts['W64_to_conjugate'] += 1
        block_preserving = set(permutation[:3]) == {0, 1, 2}
        if block_preserving:
            charge_image = transform(charge, permutation, signs)
            for sign in (-1, 1):
                if charge_image == tuple(sign * t for t in charge):
                    counts['SM_root_blocks_and_Y_line'] += 1
                    sm_signs[str(sign)] += 1

identity = tuple(range(5))
positive = (1,) * 5
square_preserves = {transform(v, identity, positive, True) for v in root_set} == root_set
square_preserves_w64 = {transform(v, identity, positive, True) for v in w64} == w64
if counts != Counter(signed_D5_candidates=3840, E8_lifts=1920,
                     W64_to_conjugate=1920, SM_root_blocks_and_Y_line=12):
    raise RuntimeError(dict(counts))
if sm_signs != Counter({'-1': 12}) or not square_preserves or not square_preserves_w64:
    raise RuntimeError('Control failed')
result = {
    'roots': len(root_set), 'W64': len(w64), 'counts': dict(counts),
    'SM_X_line_signs': dict(sm_signs),
    'square_with_identity_D5_preserves_E8': square_preserves,
    'square_with_identity_D5_preserves_W64': square_preserves_w64,
    'scope': 'Fixed orthogonal Cartan identification; quarter is (-y,x,z). '
             'Not the separate slot SU3 embedding, original flavor action, '
             'physical chirality, OS reflection, or native full-clock identity.'
}
path = Path(__file__).with_name('FLAVOR-CARTAN-LIFT-CHECK.json')
path.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
