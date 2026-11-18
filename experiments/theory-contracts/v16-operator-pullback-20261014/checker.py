"""Exact local spinor-to-vector connection pullback; not a raw field map."""
import argparse
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--repo', type=Path, default=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4'))
args = parser.parse_args()
checks = 0
names = []
odd = [x for x in range(256) if x.bit_count() % 2]
packet = [x for x in odd if (x & 31).bit_count() % 2 == 0]
pairs = list(itertools.combinations(range(16), 2))
phase = (1, 1j, -1, -1j)  # Gaussian integers, never approximate phases


def ck(ok, name):
    global checks
    if not ok:
        raise AssertionError(name)
    checks += 1
    if name not in names:
        names.append(name)


def gamma(a, mask):
    j, second = divmod(a, 2)
    lower = (mask & ((1 << j)-1)).bit_count()
    p = 2*lower + (1+2*((mask >> j)&1) if second else 0)
    return mask ^ (1 << j), p % 4


def word(indices, mask):
    p = 0
    for a in reversed(indices):
        mask, q = gamma(a, mask)
        p += q
    return mask, p % 4


def trace_product(left, right):
    total = 0j
    for mask in odd:
        m, p = word(right, mask)
        k, q = word(left, m)
        if k == mask:
            total += phase[(p+q) % 4]
    return total


for a in range(16):
    for b in range(a, 16):
        for mask in (0, 1, 37, 128, 255):
            ab = word((a, b), mask)
            ba = word((b, a), mask)
            ck(ab[0] == ba[0] and phase[ab[1]]+phase[ba[1]] == (2 if a == b else 0),
               'actual Clifford anticommutator')
for ab in pairs:
    for cd in pairs:
        # T_ab=2 Sigma_ab, so tr(T_ab T_cd)=-128 delta.
        ck(trace_product(ab, cd) == (-128 if ab == cd else 0),
           'all 120 spin generators: exact trace Gram')


def clean(matrix):
    return {k: v for k, v in matrix.items() if v}


def decode(omega):
    coefficients = []
    for ab in pairs:
        tr = 0j
        for col in odd:
            row, p = word(ab, col)
            tr += phase[p] * omega.get((col, row), 0)
        q = -tr/64  # -tr(Sigma Omega)/32
        ck(q.imag == 0 and q.real == int(q.real), 'source inverse coefficients are exact integers')
        coefficients.append(int(q.real))
    reconstructed = defaultdict(complex)
    for ab, q in zip(pairs, coefficients):
        if q:
            for col in odd:
                row, p = word(ab, col)
                reconstructed[row, col] += q*phase[p]/2
    residual = clean({key: omega.get(key, 0)-reconstructed.get(key, 0)
                      for key in omega.keys() | reconstructed.keys()})
    return coefficients, residual


def add(target, row, col, value):
    target[row, col] += value


def dg5(h, mask):
    result = defaultdict(complex)
    for b in range(5):
        if not mask & (1 << b):
            continue
        s1 = (-1)**((mask & ((1 << b)-1)).bit_count())
        removed = mask ^ (1 << b)
        for a in range(5):
            if removed & (1 << a):
                continue
            s2 = (-1)**((removed & ((1 << a)-1)).bit_count())
            result[removed | (1 << a)] += s1*s2*h[a][b]
    return clean(result)


def charge_conjugation(mask):
    # J=gamma_1 gamma_3 ... gamma_15 conjugation; J²=1, preserves chirality.
    return word(tuple(range(1, 16, 2)), mask)


family_index = {1: 0, 2: 1, 4: 2, 7: 3}


def source_packet_connection(h, m):
    t = sum(h[j][j] for j in range(5))
    omega = defaultdict(complex)
    for col in packet:
        for row, v in dg5(h, col).items():
            add(omega, row, col, 1j*v)
        add(omega, col, col, 1j*m[family_index[col >> 5]]*t)
    # Realification of the actual complex64 source packet, not 128 new species.
    for (row, col), value in list(omega.items()):
        rr, pr = charge_conjugation(row)
        cc, pc = charge_conjugation(col)
        add(omega, rr, cc, phase[(pr-pc) % 4]*value.conjugate())
    return clean(omega)


def realify_u5(h):
    # real(i h), interleaved real/imaginary coordinates in each original E line
    a = [[0]*16 for _ in range(16)]
    for i in range(5):
        for j in range(5):
            z = 1j*h[i][j]
            a[2*i][2*j], a[2*i][2*j+1] = int(z.real), -int(z.imag)
            a[2*i+1][2*j], a[2*i+1][2*j+1] = int(z.imag), int(z.real)
    return a


def vector_matrix(coefficients):
    a = [[0]*16 for _ in range(16)]
    for (i, j), q in zip(pairs, coefficients):
        a[i][j], a[j][i] = q, -q
    return a


basis = []
for i in range(5):
    h = [[0j]*5 for _ in range(5)]
    h[i][i] = 1
    basis.append(h)
for i, j in itertools.combinations(range(5), 2):
    for z in (1, 1j):
        h = [[0j]*5 for _ in range(5)]
        h[i][j], h[j][i] = z, z.conjugate() if isinstance(z, complex) else z
        basis.append(h)
ck(len(basis) == 25, 'complete original U5 connection basis')

for m in ((0, 0, -1, -1), (0, 0, 0, -2), (1, 1, 1, -5)):
    ck(sum(m) == -2, 'original global lift condition')
    family_charges = [m[i]+m[3]+1 for i in range(3)]
    for h in basis:
        trace_h = complex(sum(h[j][j] for j in range(5)))
        ck(trace_h.imag == 0, 'original Hermitian carrier trace is real')
        t = int(trace_h.real)
        omega = source_packet_connection(h, m)
        coeff, residual = decode(omega)
        ck(not residual, 'original packet connection lies in Spin16 image')
        actual = vector_matrix(coeff)
        expected = realify_u5(h)
        for i, a in enumerate(family_charges):
            expected[10+2*i][11+2*i] = -a*t
            expected[11+2*i][10+2*i] = a*t
        ck(actual == expected, 'independent source-packet inverse equals natural E_R plus Lambda2F transport')
        ck(all(actual[i][j] == 0 for i in range(10) for j in range(10, 16)),
           'induced original connection contains no inferred cross term')

# Join the fold's actual transported charge 2(e3+e6) to the original
# determinant-lift formula. The E8/plane dictionary is a premise; m is
# then solved from its charge matching, not freely fitted.
charge_a = (1, 0, 0)
selected_m4 = (sum(charge_a)-1)//2
selected_m = tuple(a-selected_m4-1 for a in charge_a)+(selected_m4,)
ck(selected_m == (0, -1, -1, 0) and sum(selected_m) == -2,
   'transported charge uniquely fixes original determinant lift within dictionary')
q_e = [[0]*5 for _ in range(5)]
q_e[2][2] = 2
source_charge_omega = source_packet_connection(q_e, selected_m)
q_coeff, q_residual = decode(source_charge_omega)
q_terms = {ab: q for ab, q in zip(pairs, q_coeff) if q}
ck(not q_residual and q_terms == {(4, 5): -2, (10, 11): -2},
   'fold charge and source packet give identical two-plane V16 action')
ck(all(row == col for row, col in source_charge_omega)
   and {-1j*v for v in source_charge_omega.values()} == {-2, 2},
   'transported packet charge has exact Cartan weights')

# Every actual mixed Lie generator is reconstructed without postulating its
# coefficient. This is an action/dictionary test, not a source field proof.
for a in range(10):
    for b in range(10, 16):
        omega = {}
        for col in odd:
            row, p = word((a, b), col)
            omega[row, col] = phase[p]  # 2 Sigma_ab
        coeff, residual = decode(omega)
        ck(not residual and coeff[pairs.index((a, b))] == 2
           and sum(abs(x) for x in coeff) == 2, 'all sixty mixed actions uniquely recovered')

# Use the existing original C/X implementation as the actual algebraic
# input, keeping its source-identification guard intact.
sys.path.insert(0, str(args.repo.resolve()))
from tfpt_explorer.source_realization import reconstruct_marked_fermion_source
witness = reconstruct_marked_fermion_source(0)['original_mixed_current_witness']
term, = witness['terms']
ck(term['doubled_root'] == [-2, 0, 0, 0, 0, -2, 0, 0]
   and term['coefficient'] == '1', 'original C0 X01 bracket fixes the mixed root')
# Its two real compact generators are A-A*, i(A+A*) for A=c_0 c_5.
# With the stated Clifford convention these are Sigma_0,10-Sigma_1,11
# and -Sigma_0,11-Sigma_1,10, respectively. Decode their *actual* CAR
# entries independently, rather than constructing them from that formula.
annih = {}
for col in odd:
    if not (col & 1 and col & (1 << 5)):
        continue
    after5 = col ^ (1 << 5)
    sign5 = (-1)**((col & 31).bit_count())
    after0 = after5 ^ 1
    annih[after0, col] = sign5
adjoint = {(col, row): value.conjugate() if isinstance(value, complex) else value
           for (row, col), value in annih.items()}
real_generators = [clean({key: annih.get(key, 0)-adjoint.get(key, 0)
                          for key in annih.keys() | adjoint.keys()}),
                   clean({key: 1j*(annih.get(key, 0)+adjoint.get(key, 0))
                          for key in annih.keys() | adjoint.keys()})]
for omega, expected_terms in zip(real_generators,
                                ({(0, 10): 1, (1, 11): -1},
                                 {(0, 11): -1, (1, 10): -1})):
    coeff, residual = decode(omega)
    observed = {ab: q for ab, q in zip(pairs, coeff) if q}
    ck(not residual and observed == expected_terms, 'original mixed root has explicit V16 cross action')

# A skewadjoint quartic Clifford term is genuine source data outside spin16;
# the inverse must expose it instead of silently erasing it.
quartic = {}
for col in odd:
    row, p = word((0, 2, 4, 6), col)
    quartic[row, col] = 1j*phase[p]
coeff, residual = decode(quartic)
ck(not any(coeff) and residual == quartic, 'outside-image term retained as full residual')
ck(all(word(tuple(range(16)), mask)[0] == mask
       and phase[word(tuple(range(16)), mask)[1]] == -1 for mask in odd),
   'central negative-volume is invisible on chosen half-spin module')
ck(all(word(tuple(range(16)), mask)[0] == mask
       and phase[word(tuple(range(16)), mask)[1]] == 1 for mask in range(256)
       if mask.bit_count() % 2 == 0), 'opposite chirality sees the different central sign')

result = {'research_id':'UR.SEAM.V16_OPERATOR_PULLBACK.01', 'verdict':'PARTIAL',
          'check_verdict':'PASS', 'checks':checks, 'check_categories':names,
          'trace_gram':'tr(Sigma_ab Sigma_cd)=-32 delta',
          'inverse':'omega_ab=-tr(Sigma_ab Omega_W)/32',
          'source_connection_cases':75, 'mixed_actions_checked':60,
          'source_charge_matching_case':1,
          'charge_selected_lift_within_E8_dictionary':list(selected_m),
          'charge_matching_vector_action':'-2 L_4,5 -2 L_10,11 = i Qprime in the interleaved vector convention',
          'original_mixed_root_witness':witness,
          'original_mixed_vector_generators':['L_0,10-L_1,11', '-L_0,11-L_1,10'],
          'raw_Clifford_multiplicity_map_derived':False,
          'complete_TFPT_solution':False,
          'scope':'Exact representation/connection inverse on marked V16 and original supplied-connection lift family; no new physical source or field multiplicity.'}
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({'check_verdict':'PASS','checks':checks,'source_connection_cases':75,
                  'mixed_actions_checked':60,'raw_operator_derived':False}))
