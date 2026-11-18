from pathlib import Path
from itertools import product
import json
import sys

import mpmath as mp
import numpy as np

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(repo / 'verification'))
from v17_hexagonal_resolvent import U6_power
from v330_qgamb_admissible_measure import symmetric_transfer

transfer = symmetric_transfer()[0]
shift = np.array(U6_power(1)).astype(float)
np.testing.assert_allclose(
    transfer, np.array([[1651, 1267, 1456], [1267, 1651, 1456], [1456, 1456, 1462]]) / 4374,
    rtol=0, atol=3e-16,
)
iv = mp.iv
iv.dps = 60
delta = ((794 - 7 * iv.sqrt(9961)) / 2187) ** (iv.mpf(1) / 6)
cos6 = [iv.mpf(x) / 2 for x in (2, 1, -1, -2, -1, 1)]
target = iv.mpf(1456) ** 2 / (1651 * 1462)
separation = iv.mpf('0.00000001')
entries = []
tails = []
mode_sets = []
normalized_matrices = []
for numerator, denominator in [(1, 1), (1, 3), (2, 3)]:
    y = iv.mpf(numerator) / denominator
    a, b = y * y + delta * delta, y * delta
    gram_modes = [a - 2 * b * c for c in cos6]
    # Work with G^(-n/2) / ||G^(-n/2)||; the quotient cancels from the invariant.
    normalized_modes = [iv.mpf(1)] + [iv.sqrt(gram_modes[0] / x) for x in gram_modes[1:]]
    q = normalized_modes[1]
    assert all(x.b <= q.b for x in normalized_modes[1:])
    assert q.b < 1
    n = 1
    # All nonconstant Fourier modes are bounded above by q^n.
    # Then K_d/K_0 >= 1 - 10*q^n uniformly in d.
    while True:
        lower_bound = 1 - 10 * q ** n
        if lower_bound.a > 0 and (lower_bound ** 2 - target).a > separation:
            break
        n += 1
    tails.append({'y': f'{numerator}/{denominator}', 'first_tail_power_n': n,
                  'uniform_invariant_lower_bound': float((lower_bound ** 2).a)})
    gram_float = float(a.mid) * np.eye(6) - float(b.mid) * (shift + shift.T)
    evals, evecs = np.linalg.eigh(gram_float)
    mode_sets.append(normalized_modes)
    normalized_matrices.append((evals[0] / evals, evecs))

power_ranges = [range(row['first_tail_power_n']) for row in tails]
for powers in product(*power_ranges):
    if not any(powers):
        continue
    spectrum = [iv.mpf(1)] * 6
    matrix = np.eye(6)
    for modes, (eigenvalues, eigenvectors), power in zip(mode_sets, normalized_matrices, powers):
        spectrum = [value * mode ** power for value, mode in zip(spectrum, modes)]
        matrix = matrix @ ((eigenvectors * (eigenvalues ** (power / 2))) @ eigenvectors.T)
    circulant = [sum(spectrum[k] * cos6[(k * d) % 6] for k in range(6)) / 6 for d in range(4)]
    for distance in (1, 2, 3):
        ratio = (circulant[distance] / circulant[0]) ** 2
        numerical_ratio = matrix[0, distance] ** 2 / (matrix[0, 0] * matrix[distance, distance])
        assert abs(float(ratio.mid) - numerical_ratio) < 5e-13
        gap = abs(ratio - target)
        assert gap.a > separation
        entries.append({'powers_by_y_1_1over3_2over3': powers, 'distance': distance,
                        'invariant': float(ratio.mid), 'separation_lower_bound': float(gap.a)})
closest = min(entries, key=lambda x: x['separation_lower_bound'])
# Check the other direct interpretation: an arbitrary reducing subspace of
# the product, with a common scalar and with the v330 spectral ratios.
# There are four distinct Fourier eigenvalues with multiplicities 1,2,2,1.
spectral_target = iv.mpf(64) / 729
adjacent_ratios = [[modes[j + 1] / modes[j] for j in range(3)] for modes in mode_sets]
spectral_limits = []
for ratios in adjacent_ratios:
    worst_ratio = max(ratios, key=lambda x: float(x.b)).b
    assert worst_ratio < 1
    exponent = 1
    while not (worst_ratio ** exponent).b < spectral_target.a - separation:
        exponent += 1
    spectral_limits.append(exponent)
spectral_checks = []
for powers in product(*(range(limit) for limit in spectral_limits)):
    if not any(powers):
        continue
    spectrum = [iv.mpf(1)] * 4
    for modes, power in zip(mode_sets, powers):
        spectrum = [value * mode ** power for value, mode in zip(spectrum, modes[:4])]
    for i in range(3):
        for j in range(i+1, 4):
            ratio = spectrum[j] / spectrum[i]
            gap = abs(ratio - spectral_target)
            assert gap.a > separation
            spectral_checks.append(float(gap.a))

result = {
    'scope': 'The original v330 transfer cannot be a three-coordinate principal submatrix of any nonempty product of original zero-phase C6 Green and structural inverse-root kernels at the three original y, with arbitrary positive diagonal endpoint factors and scalar. The three kernels share a Fourier basis, so order does not matter before endpoint weighting. Interleaved weights and other marked transfers are outside this test.',
    'v330_index_pair_one_based': [1, 3], 'target_exact': '1456^2/(1651*1462)', 'target': float(target.mid),
    'finite_checks': entries, 'analytic_tails': tails, 'closest_finite_check': closest,
    'reducing_subspace_test': {'target_ratio': '64/729', 'finite_pair_checks': len(spectral_checks),
                               'smallest_finite_separation': min(spectral_checks),
                               'first_tail_powers': spectral_limits},
    'original_whole_Z0': None, 'original_whole_ZU': None,
}
Path(__file__).with_name('c6-v330-diagonal-bridge.json').write_text(json.dumps(result, indent=2) + '\n')
print('target', float(target.mid), 'finite checks', len(entries), 'closest', closest)
print('tails', tails)
print('reducing subspace', result['reducing_subspace_test'])
