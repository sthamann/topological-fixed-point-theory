"""Exact test of the original five source states and their J8 response."""
import json
from pathlib import Path
import sympy as s

root = Path(__file__).resolve().parents[4]
source = root / 'experiments/theory-contracts/joint-native-algebra-selection-20261016'
car = json.loads((source / 'original_density_current_results.json').read_text())
hall = json.loads((source / 'original_current_hall_kernel_results.json').read_text())
det = s.Rational(1607041858632363297, 1073741824)
expected_r32 = {'1/4': s.Rational(1142672, 328635469),
                '3/4': s.Rational(1560464, 331641421)}
for alpha in ('1/4', '3/4'):
    case = car['cases'][alpha]
    assert case['source_input_energy'] == 32
    assert hall['cases'][alpha]['weighted_degree'] == 32
    g = s.Matrix([[s.sympify(z) for z in row] for row in
        case['all_positive_density_mode_responses']['8']['raw_gram']])
    assert g == g.conjugate().T
    assert all(s.det(g[:k, :k]) > 0 for k in range(1, 6))
    assert s.factor(g.det()) == det
    g0 = s.sympify(hall['cases'][alpha]['normal_ordered_first_mode_jet'][0]['raw_gram'][0][0])
    g32 = s.sympify(case['all_positive_density_mode_responses']['32']['raw_gram'][0][0])
    assert s.factor(g32 / g0) == expected_r32[alpha]
    print(alpha, 'J8 Gram positive definite, determinant:', det,
          'normalized J32 vacuum weight:', expected_r32[alpha])
print('PASS: exact Gram and original state data; no physical P1 identification')
