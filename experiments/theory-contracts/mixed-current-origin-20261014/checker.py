"""Original-source and mixed-current identification check; no source copies added."""
import argparse
import hashlib
import importlib.util
import itertools
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import sympy as S

parser = argparse.ArgumentParser()
parser.add_argument('--repo', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
root = args.repo.resolve()
checks = []


def ck(condition, name):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def load(path):
    spec = importlib.util.spec_from_file_location('mixed_origin_car', root / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Execute the ORIGINAL transitive adapter, including its existing source pins.
car = load('experiments/theory-contracts/microscopic-charged-car-limit/checker.py')
_, source = car.inherited()


def exact(a):
    result = S.Matrix([[S.Rational(float(z.real)) + S.I*S.Rational(float(z.imag))
                        for z in row] for row in a])
    ck(all(S.denom(S.re(z)) <= 2 and S.denom(S.im(z)) <= 2 for z in result),
       'original source Gaussian dyadic coefficients')
    return result


h0 = exact(source.strip_at_momentum(0, 8))
tx, sz = exact(source.TX), exact(source.SZ)
derivative = S.diag(*([-S.I*tx + S.I*tx.conjugate().T]*8))
x = S.diag(*([S.Matrix([[0, 1], [1, 0]])]*8))
z = S.diag(*([sz]*8))
rt, rb = S.zeros(16, 1), S.zeros(16, 1)
rt[14] = rt[15] = 1/S.sqrt(2)
rb[0], rb[1] = 1/S.sqrt(2), -1/S.sqrt(2)
v = S.Matrix.hstack(rt, rb)
ck(h0.rank() == 14 and h0**3 == h0, 'original two edge zero modes, fourteen massive modes')
ck(v.T*v == S.eye(2) and h0*v == S.zeros(16, 2), 'complete original zero mode basis')
ck(derivative == -x and v.T*derivative*v == S.diag(-1, 1), 'one original channel per orientation')
ck(h0*x+x*h0 == S.zeros(16) and z*x+x*z == S.zeros(16), 'original global crossing identity')
ck((h0+2*z).det() != 0, 'no further zero crossing at pi')
b = [S.conjugate(rt[i])*rt[j] for i in range(16) for j in range(16)]
support = [value for value in b if value]
gram = S.Matrix(support)*S.Matrix(support).conjugate().T
ck(len(support) == 4 and gram.rank() == 1, 'all 256 transverse bilinears: leading current rank one')

# In a one-channel U(1) Fock module, rational quarter-twist charge m/4 has
# weight m²/32+N, N a nonnegative integer. Exact inequality bounds |m|<=5
# if weight<=1. Thus this is exhaustive for the requested weight, not a
# finite scan claimed to prove an unbounded theorem.
one_channel_h1 = [(m, n) for m in range(-5, 6) for n in range(2)
                  if F(m*m, 32)+n == 1]
ck(one_channel_h1 == [(0, 1)], 'every quarter-twist fusion/descendant at weight one is neutral')
ck(F(1, 32)*4 == F(1, 8), 'doubled one-channel quarter twist has weight one eighth')
ck(8*F(1, 8) == 1, 'eight independent half-charge channels are necessary for stated spinor')
ns_weight_one = {n: (2*n)*(2*n-1)//2 for n in (1, 8)}
ck(ns_weight_one == {1: 1, 8: 120}, 'exact even-NS weight-one Fock dimensions')

# Independently enumerate the original rank-eight root coordinates. Doubled
# coordinates make all arithmetic integer; no Lie coefficient is inferred
# from an arbitrary chosen finite matrix representation.
integer_roots = set()
for i, j in itertools.combinations(range(8), 2):
    for a, c in itertools.product((-2, 2), repeat=2):
        r = [0]*8
        r[i], r[j] = a, c
        integer_roots.add(tuple(r))
spinor_roots = {r for r in itertools.product((-1, 1), repeat=8)
                if sum(value < 0 for value in r) % 2 == 0}
e8 = integer_roots | spinor_roots
base = {r for r in integer_roots
        if not any(r[:5]) or not any(r[5:])}
ck((len(base), len(integer_roots), len(spinor_roots), len(e8)) == (52, 112, 128, 240),
   'original D5+D3, D8 and E8 root census')
ck(all(sum(value*value for value in r) == 8 for r in e8), 'every target root has conformal weight one')
root_grades = Counter(sum(r[:5]) % 4 for r in e8)
full_grades = root_grades.copy()
full_grades[0] += 8
ck([root_grades[i] for i in range(4)] == [52, 64, 60, 64], '240 roots by original glue grade')
ck([full_grades[i] for i in range(4)] == [60, 64, 60, 64], '248 full currents by original glue grade')
ck(all(sum(r[:5]) % 4 in (1, 3) for r in spinor_roots)
   and all(sum(r[:5]) % 4 in (0, 2) for r in integer_roots),
   'actual odd glue grade already selects spinor type; vector type is even')


def closure(seed):
    current = set(seed)
    while True:
        added = {tuple(a+b for a, b in zip(r, t)) for r in current for t in current}
        added &= e8
        if added <= current:
            return current
        current |= added


lam = (1,)*8
vector = (2, 0, 0, 0, 0, 2, 0, 0)
ck(closure(base | {lam, tuple(-a for a in lam)}) == e8,
   'base currents plus actual spinor and adjoint generate all E8 roots')
ck(closure(base | {vector, tuple(-a for a in vector)}) == integer_roots,
   'base currents plus mixed vector and adjoint generate only D8 roots')
ck(all(all(value % 2 == 0 for value in r) for r in integer_roots),
   'no spinor charge in even-NS current sector')

# Audit Claude's extended selector in the stipulated root-lattice class.
# For d>9 every nontrivial mixed D-spinor is already too long; pure A
# solutions are bounded by k<=4 after k<->d-k, and vector equation gives
# d=k+1+1/(k-1). The manuscript supplies the analytic all-d argument.
spin_solutions, vector_solutions, pure_a = [], [], []
for d in range(2, 201):
    for k in range(d):
        a = F(k*(d-k), d)
        if a+F(d+1, 4) == 2:
            spin_solutions.append([d, k])
        if k and a+1 == 2:
            vector_solutions.append([d, k])
        if a == 2:
            pure_a.append([d, k])
ck(spin_solutions == [[4, 1], [4, 3], [7, 0]], 'spinor selector including pure d7 exception')
ck(vector_solutions == [[4, 2]], 'mixed vector selector')
ck(pure_a == [[8, 4], [9, 3], [9, 6]], 'pure A exceptions')
ck(F(16, 4**2) == 1 and F(16, 2**2) == 4, 'spinor/vector extension determinants')

# Same scalar two-point, different Ward charge: exact counterexample to the
# pasted identification. Coordinates are separated; these are target
# analytic formulas, not measured microscopic charged field responses.
z0, w, u = F(0), F(1), F(3)
pair = 1/(w-u)**2
ward = [F(1, 2)*pair*(1/(z0-w)-1/(z0-u)) for _ in range(8)]
ck(pair == F(1, 4) and ward == [F(-1, 12)]*8,
   'spinor Cartan Ward fingerprint distinct from neutral abelian current zero fingerprint')
coords = [F(0), F(1), F(3), F(7)]
charges = [1, -1, 1, -1]
vertex_four = F(1)
for i, j in itertools.combinations(range(4), 2):
    vertex_four *= (coords[i]-coords[j])**(2*charges[i]*charges[j])
neutral_four = sum(1/((coords[a]-coords[b])**2*(coords[c]-coords[d])**2)
                   for a, b, c, d in [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)])
ck(vertex_four == F(81, 784) and neutral_four != vertex_four,
   'identical weight-one scalar two-point is not identical field response')
ck(sum(F(1, 4)**2 for _ in range(8))/2 == F(1, 4)
   and sum(F(1, 2)**2 for _ in range(8))/2 == 1
   and sum(F(1, 4)**2 for _ in range(8)) == F(1, 2),
   'conditional double-quarter-twist OPE weight and exponent')
ck(F(1, 2) not in [F(n) for n in range(-2, 3)], 'spinor is outside each elementary CAR integer charge lattice')

result = {
    'research_id': 'UR.SEAM.MIXED_CURRENT_ORIGIN.01',
    'verdict': 'PARTIAL', 'check_verdict': 'PASS', 'checks': len(checks),
    'check_names': checks,
    'original_source': {'copies': 1, 'width': 8, 'top_complex_channels': 1,
                        'transverse_current_rank': 1, 'weight_one_dimension': 1,
                        'quarter_fusion_weight_one_solutions': one_channel_h1},
    'conditional_eight_copy_target': {'even_NS_weight_one': 120, 'spinor_extension': 128,
                                    'E8_weight_one': 248, 'source_multiplicity_derived': False},
    'selector': {'spinor': spin_solutions, 'vector': vector_solutions, 'pure_A': pure_a},
    'pair_counterexample': {'pair': str(pair), 'spinor_ward': [str(x) for x in ward],
                            'neutral_ward': ['0']*8, 'vertex_four': str(vertex_four),
                            'neutral_four': str(neutral_four)},
    'complete_TFPT_solution': False, 'local_mixed_spinor_from_primitive_seam': False,
    'source_sha256': hashlib.sha256((root/'verification/v1033_charged_disorder.py').read_bytes()).hexdigest()
}
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({'check_verdict': 'PASS', 'checks': len(checks),
                  'source_current_rank': 1, 'local_origin_derived': False}))
