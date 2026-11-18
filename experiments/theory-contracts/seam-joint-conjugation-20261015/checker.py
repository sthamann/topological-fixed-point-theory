"""Exact, scoped checks for UR.SEAM.JOINT_CONJUGATION.01. No origin claim."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def require(ok, label):
    if not ok:
        raise AssertionError(label)


def creation(j, n=5):
    out = np.zeros((1 << n, 1 << n), dtype=np.int64)
    for m in range(1 << n):
        if not (m >> j) & 1:
            out[m | (1 << j), m] = (-1) ** (m & ((1 << j) - 1)).bit_count()
    return out


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def act(matrix, v):
    return tuple(dot(row, v) for row in matrix)


def run():
    counts = Counter()

    def equal(x, y, group):
        require(np.array_equal(x, y), group)
        counts[group] += 1

    e = [creation(i) for i in range(5)]
    a = [v.T for v in e]
    I = np.eye(32, dtype=np.int64)
    N = [e[i] @ a[i] for i in range(5)]
    B = I.copy()
    for i in range(5):
        B = B @ (e[i] + a[i])
    equal(B @ B, I, 'spin_conjugation_square')
    chirality = np.diag([(-1) ** n.bit_count() for n in range(32)])
    equal(B @ chirality @ B, -chirality, 'spin_chirality')
    generators = [2 * n - I for n in N]  # doubled diagonal Cartans
    for i, j in combinations(range(5), 2):
        g = e[i] @ a[j]
        generators += [g + g.T, 1j * (g - g.T)]
        p = e[i] @ e[j]
        generators += [p + p.T, 1j * (p - p.T)]
    require(len(generators) == 45, 'all Spin10 directions')
    for g in generators:
        equal(B @ g.conjugate() @ B, -g, 'Spin10_antiunitary_generators')
        equal(g.T @ B + B @ g, np.zeros((32, 32)), 'Chevalley_bilinear_invariance')
        # Hermitian same-sheet invariant form: skew generator i*g satisfies this.
        ix = np.ix_([m for m in range(32) if m.bit_count() % 2 == 0],
                    [m for m in range(32) if m.bit_count() % 2 == 0])
        x = (1j * g)[ix]
        equal(x.conjugate().T + x, np.zeros((16, 16)), 'same_sheet_Hermitian_scalar')
    for m in range(32):
        indices = [i + 1 for i in range(5) if m >> i & 1]
        equal(B[:, m], (-1) ** sum(5 - i for i in indices) * I[:, 31 - m],
              'Chevalley_complement_signs')
    even = [m for m in range(32) if m.bit_count() % 2 == 0]
    require(F(sum(int((2 * N[0] - I)[m, m]) ** 2 for m in even), 4) == 4,
            'complex spinor current coefficient is 4')
    counts['species_and_level_guard'] += 1

    # Actual finite CAR lift, on excitation Fock space. The two internal labels
    # here test a general algebraic law; they are not added source species.
    fcreate = [creation(j, 4) for j in range(4)]
    fnumber = [v @ v.T for v in fcreate]
    energies = [-3, -1, 1, 3]
    fields = [fcreate[j] if energies[j] < 0 else fcreate[j].T for j in range(4)]
    U = np.zeros((16, 16), dtype=np.int64)
    for mask in range(16):
        js = [j for j in range(4) if mask >> j & 1]
        inversions = len(js) * (len(js)-1) // 2  # reversal of all occupied labels
        U[sum(1 << (3-j) for j in js), mask] = (-1) ** inversions
    equal(U @ U.T, np.eye(16), 'polarized_Fock_unitarity')
    equal(U[:,0], np.eye(16)[:,0], 'polarized_vacuum_transport')
    for j in range(4):
        equal(U @ fields[j] @ U.T, fields[3-j].T, 'polarized_CAR_field_transport')
    def normal_charge(q):
        return (sum(q[i,j] * fields[i].T @ fields[j] for i in range(4) for j in range(4))
                - sum(q[i,i] for i in (0,1)) * np.eye(16))
    qbasis = [np.diag([int(i==j) for i in range(4)]) for j in range(4)]
    for i,j in combinations(range(4),2):
        for z in (1,1j):
            qtest = np.zeros((4,4), dtype=np.complex128)
            qtest[i,j], qtest[j,i] = z, np.conjugate(z)
            qbasis.append(qtest)
    for qtest in qbasis:
        jqj = qtest.conjugate()[::-1,::-1]
        equal(U @ normal_charge(qtest) @ U.T, -normal_charge(jqj), 'general_normal_ordering_law')
    for qtest, sign in [(np.eye(4),-1),(np.diag([-2,3,-3,2]),1),(np.diag(energies),1)]:
        equal(U @ normal_charge(qtest) @ U.T, sign * normal_charge(qtest), 'number_charge_time_signs')
    equal(normal_charge(np.diag(energies)), sum(abs(energies[j])*fnumber[j] for j in range(4)),
          'positive_hole_energy')

    # Reuse the actual M, not an invented replacement or a fitted permutation.
    fold = ROOT / 'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json'
    M = [[F(x) for x in row] for row in json.loads(fold.read_text())['explicit_Weyl_matrix']]
    O = [[-sum(M[i][k] * M[j][7-k] for k in range(8)) for j in range(8)] for i in range(8)]
    expected_columns = [(3, 1), (4, -1), (5, -1), (0, 1), (1, -1), (2, -1), (7, 1), (6, 1)]
    for j, (row, sign) in enumerate(expected_columns):
        require([O[i][j] for i in range(8)] == [sign * (i == row) for i in range(8)], 'actual O')
    roots = set()
    for i, j in combinations(range(8), 2):
        for si, sj in product([-1, 1], repeat=2):
            v = [F(0)] * 8
            v[i], v[j] = F(si), F(sj)
            roots.add(tuple(v))
    for signs in product([-1, 1], repeat=8):
        if sum(s == -1 for s in signs) % 2 == 0:
            roots.add(tuple(F(s, 2) for s in signs))
    require(len(roots) == 240 and {act(O, r) for r in roots} == roots, 'E8 root preservation')
    counts['E8_root_images'] = 240
    d5 = {r for r in roots if all(x == 0 for x in r[5:])}
    d3 = {r for r in roots if all(x == 0 for x in r[:5])}
    od5, od3 = {act(O, r) for r in d5}, {act(O, r) for r in d3}
    require((len(d5), len(d3), len(d5 & od5), len(d3 & od3)) == (40, 12, 24, 4), 'marking intersections')
    counts['marking_intersection_census'] = 1
    moment = [[sum(r[i] * r[j] for r in od5 | od3) for j in range(8)] for i in range(8)]
    transported_axes = {0, 1, 3, 4, 5}
    require(moment == [[(16 if i in transported_axes else 8) * (i == j) for j in range(8)] for i in range(8)], 'covariant marking recovered from neutral roots')
    counts['transported_marking_reconstruction'] = 1

    y = tuple(map(F, ['-1/3','-1/3','-1/3','1/2','1/2','0','0','0']))
    q = tuple(map(F, [0,0,2,0,0,2,0,0]))
    require(act(O, q) == tuple(-v for v in q), 'source charge reverses')
    require(act(O, y) == tuple(map(F, ['1/2','-1/2','0','-1/3','1/3','1/3','0','0'])), 'actual hypercharge image')
    counts['charge_images'] = 2
    quarter_data = []
    for alpha in (F(1,4), F(3,4)):
        beta = 1 - alpha
        delta = tuple((j + alpha) / 8 - F(1,2) for j in range(8))
        other = tuple((j + beta) / 8 - F(1,2) for j in range(8))
        d, db = act(M, delta), act(M, other)
        require(act(O, d) == db, 'quarter time transported')
        for n in range(-32, 32):
            l, j = divmod(n, 8)
            require(divmod(-n-1, 8) == (-l-1, 7-j), 'actual mode fold')
            require((-n-1)+beta == -(n+alpha), 'quarter domain')
            counts['original_quarter_modes'] += 1
        for r in roots:
            require(dot(d, r) == dot(db, act(O, r)), 'all current Ward coefficients')
            counts['root_Ward_transport'] += 1
        fixed = sum(db[i] ** 2 for i in range(5))
        moved = sum(db[i] ** 2 for i in transported_axes)
        require(fixed != moved and abs(fixed - moved) == F(1,32), 'no time-stabilizing orthogonal repair')
        quarter_data.append({'source_alpha':str(alpha),'target_beta':str(beta),
                             'fixed_D5_time_norm':str(fixed),'transported_D5_time_norm':str(moved)})
    counts['orthogonal_obstruction_cases'] = 2

    # Concrete all-integer cocycle lift on an E8 lattice basis. Identity extends
    # to all integer coordinates by bilinearity/polarization, not by this census.
    cols = [[F(1,2)] + [F(-1,2)] * 6 + [F(1,2)], [1,1,0,0,0,0,0,0]]
    for i in range(6):
        v = [0] * 8
        v[i], v[i+1] = -1, 1
        cols.append(v)
    L = sp.Matrix.hstack(*(sp.Matrix(v) for v in cols))
    require(abs(L.det()) == 1, 'E8 basis determinant')
    G = L.T * L
    C = sp.zeros(8)
    for i in range(8):
        C[i,i] = G[i,i] / 2
        for j in range(i):
            C[i,j] = G[i,j]
    W = L.inv() * sp.Matrix(O) * L
    require(all(v.q == 1 for v in W), 'integral E8 action')
    A = W.T * C * W - C
    require(all(int(A[i,j]-A[j,i]) % 2 == 0 for i in range(8) for j in range(8)), 'cocycle symmetry')
    ai = [[int(A[i,j]) for j in range(8)] for i in range(8)]
    def phase(x):
        return (sum(ai[i][i] * x[i]*(x[i]-1)//2 for i in range(8))
                + sum(ai[i][j]*x[i]*x[j] for i in range(8) for j in range(i+1,8))) % 2
    samples = list(product([0,1], repeat=8))
    # Basis/negative/large-coordinate tests supplement the analytic polarization.
    probes = [tuple(int(i==j) for i in range(8)) for j in range(8)]
    probes += [tuple(range(-4,4)), (1000001,-3,2,0,1,-8,9,2)]
    for x in samples:
        for z in probes:
            require((phase(tuple(a+b for a,b in zip(x,z)))-phase(x)-phase(z)
                     -sum(x[i]*ai[i][j]*z[j] for i in range(8) for j in range(8))) % 2 == 0,
                    'vertex cocycle lift')
            counts['cocycle_polarization_cases'] += 1

    return {'research_id':'UR.SEAM.JOINT_CONJUGATION.01','verdict':'EXACT_COVARIANT_QUARTER_TRANSPORT',
            'checks':dict(counts),'checks_total':sum(counts.values()),
            'O':[[int(v) for v in r] for r in O], 'quarters':quarter_data,
            'D5_roots_retained':24,'D5_roots_total':40,'A3_roots_retained':4,'A3_roots_total':12,
            'first_marking_given':True,'opposite_marking_derived_by_transport':True,
            'same_fixed_marking_preserved':False,'quarter_phase_erased':False,
            'raw_P1P2_derivation':False,'E8_extension_from_raw_source':False,
            'complete_TFPT_solution':False,
            'fold_sha256':hashlib.sha256(fold.read_bytes()).hexdigest()}


if __name__ == '__main__':
    result = run()
    (HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
