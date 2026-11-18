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

import itertools as it, json, pickle
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ
root = str(WORK) + '/'
refs, perms, pg, chosen, orbits, full, orbit_id, raykeys = pickle.load(open(root + 'state.pkl', 'rb'))
d = json.load(open(root + 'results.json'))
supports = [tuple(sorted((2 * i, 2 * i + 1, 2 * j, 2 * j + 1))) for i, j in it.combinations(range(4), 2)]
supports += [tuple((2 * k + t[k] for k in range(4))) for t in it.product((0, 1), repeat=4) if sum(t) % 2 == 0]
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
B = hermite_normal_form(s.Matrix.hstack(*[s.Matrix(x) for x in roots]))
Bi = B.inv()
Gram = B.T * B / 2
require(Gram.det() == 1 and all((Gram[i, i] % 2 == 0 for i in range(8))), 'Exact check, original line 21: Gram.det() == 1 and all((Gram[i, i] % 2 == 0 for i in range(8)))')
ints = []
for A, C in refs:
    M = s.zeros(8)
    for i in range(4):
        for j in range(4):
            M[2 * i, 2 * j] = s.Rational(int(A[i, j]), 2)
            M[2 * i, 2 * j + 1] = -s.Rational(int(C[i, j]), 2)
            M[2 * i + 1, 2 * j] = s.Rational(int(C[i, j]), 2)
            M[2 * i + 1, 2 * j + 1] = s.Rational(int(A[i, j]), 2)
    R = Bi * M * B
    require(all((x.q == 1 for x in R)) and R.T * Gram * R == Gram, 'Exact check, original line 30: all((x.q == 1 for x in R)) and R.T * Gram * R == Gram')
    ints.append(R)
for o in d['orbits']:
    D = s.Matrix.vstack(*[ints[i] - s.eye(8) for i in o['representative'][:3]])
    SN = smith_normal_form(D, domain=ZZ)
    nonzero = [abs(int(SN[i, i])) for i in range(8) if SN[i, i] != 0]
    o['integral_H1'] = {'free_rank': 24 - len(nonzero), 'torsion_invariant_factors': [x for x in nonzero if x > 1], 'd0_smith_nonzero': nonzero}
    print(o['representative'], o['integral_H1'])
d['integral_lattice_basis'] = list(map(list, B.tolist()))
d['integral_gram'] = Gram.tolist()
open(root + 'results.json', 'w').write(json.dumps(d, indent=2, default=int))
