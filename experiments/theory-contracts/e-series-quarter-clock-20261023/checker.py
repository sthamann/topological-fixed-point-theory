"""Exact finite probe of the original E8 roots and Gaussian quarter clock.

Exploratory mathematics; no change to TFPT's verification ledger.
Native functions are extracted without executing their module-level suites.
"""
import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
HERE = Path(__file__).resolve().parent
PINS = {}


def native(relative, name):
    path = REPO / relative
    raw = path.read_bytes()
    PINS[relative] = sha256(raw).hexdigest()
    tree = ast.parse(raw)
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    namespace = {'np': np, 'itertools': itertools}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace[name]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


roots = np.array(native('verification/v1_e8_glue.py', 'e8_roots')())
require(np.array_equal(2 * roots, np.rint(2 * roots)), 'half-integral roots')
R = np.rint(2 * roots).astype(np.int64)
Jvec = native('verification/v634_st31_structure.py', 'J_vec')
J = np.column_stack([Jvec(v) for v in np.eye(8, dtype=np.int64)])
I = np.eye(8, dtype=np.int64)
lookup = {tuple(v): i for i, v in enumerate(R)}
jp = np.array([lookup[tuple(Jvec(v))] for v in R])
G = R @ R.T  # four times the actual inner products
require(len(lookup) == 240 and np.all(np.diag(G) == 8), 'E8 shell')
require(np.array_equal(J @ J, -I) and np.array_equal(J.T @ J, I), 'native J')

e7_pairs = set()
reflection_matrices = set()
reflection_permutations = {}
intersection_counts = Counter()
shells = Counter()
for a in range(240):
    b = int(jp[a])
    require(G[a, b] == 0, 'alpha perpendicular Jalpha')
    e7a = frozenset(np.flatnonzero(G[a] == 0))
    e7b = frozenset(np.flatnonzero(G[b] == 0))
    require(len(e7a) == len(e7b) == 126 and e7a != e7b, 'two E7 slices')
    common = e7a & e7b
    require(len(common) == 60, 'D6 common roots')
    intersection_counts[len(common)] += 1
    shells[tuple(sorted(Counter(int(x // 4) for x in G[a]).items()))] += 1
    e7_pairs.add(tuple(sorted((tuple(sorted(e7a)), tuple(sorted(e7b))))))
    # 4*s_alpha = 4I - (2alpha)(2alpha)^T, exact integer arithmetic.
    sa4 = 4 * I - np.outer(R[a], R[a])
    sb4 = 4 * I - np.outer(R[b], R[b])
    numerator = sa4 @ sb4
    require(np.all(numerator % 4 == 0), 'half-integral double-reflection')
    refl4 = numerator // 4  # four times s_alpha*s_Jalpha
    require(np.array_equal(refl4 @ J, J @ refl4), 'complex linearity')
    require(np.array_equal(refl4 @ refl4, 16 * I), 'involution')
    require(np.trace(refl4) == 16, 'real trace4, complex trace2')
    images4 = R @ refl4.T
    require(np.all(images4 % 4 == 0), 'root image integrality in doubled chart')
    require(all(tuple(v) in lookup for v in images4 // 4), 'root preservation')
    reflection_matrices.add(tuple(int(v) for v in refl4.flat))
    reflection_permutations[tuple(int(v) for v in refl4.flat)] = [lookup[tuple(v)] for v in images4 // 4]
require(len(e7_pairs) == len(reflection_matrices) == 60, '60 E7 pairs / reflections')
reflection_group_order = PermutationGroup([
    Permutation(p) for p in reflection_permutations.values()
]).order()
require(reflection_group_order == 46080, 'original G31 order')

# The original group acts transitively on the entire60-chart family.
lines = sorted({tuple(sorted((i, int(jp[i]), int(jp[jp[i]]), int(jp[jp[jp[i]]]))))
                for i in range(240)})
line_id = {root: k for k, line in enumerate(lines) for root in line}
line_actions = [[line_id[p[line[0]]] for line in lines]
                for p in reflection_permutations.values()]
orbit, frontier = {0}, [0]
while frontier:
    nxt = {p[x] for x in frontier for p in line_actions} - orbit
    orbit.update(nxt)
    frontier = list(nxt)
require(len(orbit) == 60, 'transitive60-chart family')

# Root brackets of the two E7 slices recover every initially missing root.
bracket_closures = Counter()
root_sum_table = np.full((240, 240), -1, dtype=np.int64)
for u, v in zip(*np.where(G == -4)):
    root_sum_table[u, v] = lookup[tuple(R[u] + R[v])]
for a in range(240):
    b = int(jp[a])
    initial = set(np.flatnonzero((G[a] == 0) | (G[b] == 0)))
    require(len(initial) == 192, 'E7 pair union192 roots')
    generated = initial | set(int(x) for x in np.unique(root_sum_table[np.ix_(sorted(initial), sorted(initial))]) if x >= 0)
    require(len(generated) == 240, 'one bracket round generates all E8 roots')
    bracket_closures[(len(initial), len(generated))] += 1

# Integral quotient groups retain the actual gluing type, not only its index.
B2 = sp.Matrix(np.rint(2 * np.array(native('verification/v1_e8_glue.py', 'e8_simple_roots')())).astype(int))
root_coordinates = sp.Matrix(R) * B2.inv()
require(all(x.q == 1 for x in root_coordinates), 'integral E8 root coordinates')


def quotient_group(indices):
    generators = sp.Matrix([[int(root_coordinates[i, j]) for j in range(8)] for i in sorted(indices)])
    basis = hermite_normal_form(generators.T)
    require(basis.shape == (8, 8), 'full-rank sublattice')
    sf = smith_normal_form(basis, domain=sp.ZZ)
    return [abs(int(sf[i, i])) for i in range(8) if abs(sf[i, i]) != 1]


qalpha = lookup[(2, 2, 0, 0, 0, 0, 0, 0)]
qbeta = int(jp[qalpha])
normal4 = {qalpha, qbeta, lookup[tuple(-R[qalpha])], lookup[tuple(-R[qbeta])]}
d6 = set(np.flatnonzero((G[qalpha] == 0) & (G[qbeta] == 0)))
e7q = set(np.flatnonzero(G[qalpha] == 0))
e7j = set(np.flatnonzero(G[qbeta] == 0))
d5d3 = set(np.flatnonzero(np.all(R[:, 5:] == 0, axis=1) | np.all(R[:, :5] == 0, axis=1)))
glue = {'D6+A1+A1': quotient_group(d6 | normal4),
        'E7+A1': quotient_group(e7q | {qalpha, lookup[tuple(-R[qalpha])]}),
        'D5+A3': quotient_group(d5d3),
        'sum_of_two_E7_root_lattices': quotient_group(e7q | e7j)}
require(glue == {'D6+A1+A1': [2, 2], 'E7+A1': [2], 'D5+A3': [4],
                 'sum_of_two_E7_root_lattices': []}, 'distinct actual glue groups')

# Every regular A2 subsystem consists of ±alpha, ±beta, ±(alpha+beta).
a2s = set()
for a, b in zip(*np.where(G == -4)):
    s = R[a] + R[b]
    a2s.add(tuple(sorted(lookup[tuple(v)] for v in (R[a], -R[a], R[b], -R[b], s, -s))))
require(len(a2s) == 1120, 'all A2 subsystems')
a2_fixed = 0
a2_pairs = set()
e6_intersections = Counter()
for sub in a2s:
    jsub = tuple(sorted(int(jp[i]) for i in sub))
    require(jsub in a2s, 'J transports A2')
    a2_fixed += int(jsub == sub)
    a2_pairs.add(tuple(sorted((sub, jsub))))
    e6a = set(np.flatnonzero(np.all(G[list(sub)] == 0, axis=0)))
    e6b = set(np.flatnonzero(np.all(G[list(jsub)] == 0, axis=0)))
    require(len(e6a) == len(e6b) == 72, 'E6 orthogonal complement')
    e6_intersections[len(e6a & e6b)] += 1
require(a2_fixed == 0 and len(a2_pairs) == 560, 'no J-fixed regular E6 slice')

# A compatible D5 -> E6 -> E7 -> E8 flag in the actual standard D5+D3 chart.
def idx(v):
    return lookup[tuple(2 * x for x in v)]

alpha = idx((0, 0, 0, 0, 0, 1, -1, 0))
beta = idx((0, 0, 0, 0, 0, 0, 1, -1))
e7flag = set(np.flatnonzero(G[alpha] == 0))
e6flag = set(np.flatnonzero((G[alpha] == 0) & (G[beta] == 0)))
d5flag = set(np.flatnonzero(np.all(R[:, 5:] == 0, axis=1)))
require(len(d5flag) == 40 and d5flag < e6flag < e7flag, 'compatible nested flag')
a2flag = {alpha, beta, lookup[tuple(-R[alpha])], lookup[tuple(-R[beta])],
          lookup[tuple(R[alpha] + R[beta])], lookup[tuple(-R[alpha] - R[beta])]}
glue['E6+A2'] = quotient_group(e6flag | a2flag)
require(glue['E6+A2'] == [3], 'actual index3 E6 glue')

# Analytic covering statements are checked by their exact finite formulas.
cover_data = {}
for n in (2, 3, 4):
    g = Fraction((n - 1) * 3 - sp.gcd(n, 4) + 1, 2)
    cover_data[str(n)] = {'equation': f'w**{n}=x**4-1', 'compact_genus': int(g),
                          'ramification_index_at_infinity': int(n // sp.gcd(n, 4))}
require([cover_data[str(n)]['compact_genus'] for n in (2, 3, 4)] == [1, 3, 3], 'Riemann-Hurwitz')
lam = Fraction(-1)
jcurve = 256 * (1 - lam + lam ** 2) ** 3 / (lam ** 2 * (1 - lam) ** 2)
require(jcurve == 1728, 'harmonic elliptic j')
quarter_lifts = [list(bits) for bits in itertools.product((1, 3), repeat=4) if sum(bits) % 4 == 0]
cyclic_lifts = [bits for bits in quarter_lifts if bits[1:] + bits[:1] == bits]
require(len(quarter_lifts) == 8 and len(cyclic_lifts) == 2, 'cover does not select quarter orientation')
normalization = {name: {'rank': rank, 'lattice_extension_index': index,
                         'inverse_rank_coefficient_times_pi': str(Fraction(1, rank)),
                         'inverse_twice_index_coefficient_times_pi': str(Fraction(1, 2 * index))}
                 for name, rank, index in [('E6+A2', 6, 3), ('E7+A1', 7, 2), ('D5+A3', 8, 4)]}
require(normalization['E7+A1']['inverse_rank_coefficient_times_pi'] !=
        normalization['E7+A1']['inverse_twice_index_coefficient_times_pi'], 'E7 distinguishes normalization prescriptions')

levels = {name: str(Fraction(dim, 1 + h)) for name, dim, h in
          [('E6', 78, 12), ('E7', 133, 18), ('E8', 248, 30)]}
result = {
    'scope': 'standard regular root subsystems with fixed original Gaussian J; no all-embedding claim',
    'source_sha256': PINS,
    'E7': {'root_slices': 120, 'J_fixed_slices': 0, 'J_paired_slices': len(e7_pairs),
           'common_root_count': 60, 'common_root_type': 'D6',
           'distinct_complex_root_reflections': len(reflection_matrices),
           'generated_reflection_group_order': int(reflection_group_order),
           'orbit_size_on_60_pairs': len(orbit),
           'one_root_bracket_round': [{'before': k[0], 'after': k[1], 'roots_tested': v}
                                      for k, v in bracket_closures.items()],
           'inner_product_shells': [{'multiplicities': dict(k), 'roots_tested': v} for k, v in shells.items()]},
    'E6': {'A2_subsystems': len(a2s), 'J_fixed_A2': a2_fixed, 'J_pairs': len(a2_pairs),
           'paired_E6_intersection_root_counts_over_all_1120': dict(e6_intersections)},
    'compatible_D5_E6_E7_E8_flag': {'root_counts': [40, 72, 126, 240],
                                   'alpha': [0, 0, 0, 0, 0, 1, -1, 0],
                                   'beta': [0, 0, 0, 0, 0, 0, 1, -1],
                                   'marking_chosen': True},
    'level_one_sugawara_central_charges': levels,
    'actual_integral_quotient_invariant_factors': glue,
    'derived_cyclic_curve_data': cover_data,
    'double_cover': {'compact_j': int(jcurve), 'restricted_holonomy_exponents_mod4': [0, 2],
                     'branch_puncture_fibre_holonomy_exponent': 2,
                     'extends_as_ordinary_flat_bundle_over_compact_curve': False,
                     'quarter_lifts_of_parity_character': quarter_lifts,
                     'cyclically_invariant_quarter_lifts': cyclic_lifts},
    'normalization_diagnostic': normalization,
    'seed': {'c3': '1/(8*pi)', 'phi_base': '1/(6*pi)', 'ratio': '4/3',
             'phi0': '1/(6*pi)+3/(256*pi**4)', 'E7_window_historically_present': True},
}
for relative, digest in PINS.items():
    require(sha256((REPO / relative).read_bytes()).hexdigest() == digest, 'source unchanged')
(HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
