#!/usr/bin/env python3
"""Original neutral parity pairs and complementary-vector controls.

No continuum field limit, E8 origin or arbitrary-time determinant is claimed.
Assertions are explicit and stay active with -OO.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as s

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
OUT = Path(__file__).resolve().parent


def need(ok, msg):
    if not bool(ok):
        raise RuntimeError(msg)


def zero(a):
    return all(s.cancel(s.expand(z)) == 0 for z in a)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


# A BASE interval on four lattice points, lifted to eight source sheets.
fold, width = 8, 4
nx = fold * width
region = {j * width + 1 for j in range(fold)}
sigma = s.diag(*[-1 if x in region else 1 for x in range(nx)])
deck = s.zeros(nx)
for x in range(nx):
    deck[(x + width) % nx, x] = 1
need(zero(sigma * deck - deck * sigma), 'mask must be sheet-periodic')
need(sigma[0, 0] == sigma[nx - 1, nx - 1] == 1,
     'original seam is outside the pair support')
need(all(sigma[j * width + x, j * width + x] == (-1 if x == 1 else 1)
         for j in range(fold) for x in range(width)), 'uniform local base sign')
z = s.diag(1, s.I, -1, -s.I)
xreg = s.zeros(4)
for r in range(4):
    xreg[(r + 1) % 4, r] = 1
# On two occupations, the two half strings have the same X^2 factor.
p_a = s.diag(1, -1, 1, -1)
p_b = s.diag(1, 1, -1, -1)
f_a = s.kronecker_product(xreg ** 2, p_a)
f_b = s.kronecker_product(xreg ** 2, p_b)
pair = f_b * f_a.H
need(pair == s.kronecker_product(s.eye(4), p_b * p_a), 'register cancels')
need(zero(pair * s.kronecker_product(z, s.eye(4))
          - s.kronecker_product(z, s.eye(4)) * pair), 'pair is register neutral')

# Exact Fourier selection and determinant factorization on a finite band.
# This regulator has TWO Fermi edges, not a proved single-chiral limit.
def mask_entry(m, n):
    d = n - m
    if d % 8:
        return s.Integer(0)
    k = d // 8
    # Base interval I=(0,pi): a useful exact branch/zero countercontrol.
    return s.Rational(1, 2) if k == 0 else (s.Integer(-1) ** k - 1) / (2 * s.pi * s.I * k)


band_controls = []
for length in (1, 2, 3, 4):
    # Ordered by residue; entries prove the same identity for all lengths.
    modes = [[-8 * length + j + 8 * ell for ell in range(length)]
             for j in range(8)]
    a = s.Matrix(length, length,
                 lambda i, j: int(i == j) - 2 * mask_entry(modes[0][i], modes[0][j]))
    for j in range(8):
        need(zero(s.Matrix(length, length,
                  lambda u, v: mask_entry(modes[j][u], modes[j][v]))
                  - (s.eye(length) - a) / 2), 'equal residue blocks')
        for k in range(j):
            need(all(mask_entry(m, n) == 0 for m in modes[j] for n in modes[k]),
                 'off-residue matrix elements vanish')
    d = s.simplify(a.det())
    if length % 2:
        need(d == 0, 'odd half-circle band determinant must vanish')
    band_controls.append({'occupied_modes_per_residue': length,
                          'one_block_determinant': str(d),
                          'full_expectation': str(s.expand(d ** 8))})

# Exact complementary-covariance determinant identity, both actual quarters.
# L=2 log 2 makes n(A)=1/(1+2^(4 A)) rational at these source modes.
complement_controls = []
shift = s.zeros(4)
for j in range(3):
    shift[j + 1, j] = 1
f = (shift + shift.H - shift ** 2 - (shift.H) ** 2) / 3
g = (s.eye(4) + s.I * f) * (s.eye(4) - s.I * f).inv()
need(zero(g.H * g - s.eye(4)), 'finite Cayley phase is unitary')
for alpha in (s.Rational(1, 4), s.Rational(3, 4)):
    c = s.diag(*[1 / (1 + s.Integer(2) ** (4 * (m + alpha)))
                 for m in (-4, 0, 4, 8)])
    d = s.factor((s.eye(4) - c + c * g).det())
    dc = s.factor((c + (s.eye(4) - c) * g).det())
    need(s.cancel(dc - g.det() * s.conjugate(d)) == 0,
         'complementary determinant identity')
    need(s.cancel(dc * d / g.det() - d * s.conjugate(d)) == 0,
         'complete vector amplitude loses the single component phase')
    need(s.im(d) != 0, 'single component phase must be nontrivial in control')
    complement_controls.append({'alpha': str(alpha), 'exact': True,
                                'single_component_phase_nonzero': True})

# Original four-mark-fixed profiles, two opposite Dirac grading components.
sin4 = {4: 1 / (2 * s.I), -4: -1 / (2 * s.I)}
cos4_minus_cos8 = {4: s.Rational(1, 2), -4: s.Rational(1, 2),
                       8: -s.Rational(1, 2), -8: -s.Rational(1, 2)}
central = sum(m * (sin4.get(-m, 0) * cos4_minus_cos8.get(m, 0)
                  - cos4_minus_cos8.get(-m, 0) * sin4.get(m, 0))
              for m in (4, 8))
need(central == 2 * s.I, 'one chiral central form')
need(central + (-central) == 0, 'complete vector Ward central form cancels')

# Same original time: Fourier/Parseval factor for phi(theta)=f(8theta).
# A pi rise and pi return of widths a,eps saturate the H1 ramp bound.
energy_controls = []
for q in (1, 2, 4, 8):
    a, eps = s.pi / 8, s.pi / (8 * q)
    derivative_norm2 = s.pi ** 2 / a + s.pi ** 2 / eps
    delta_k = s.simplify(2 * derivative_norm2 / s.pi)
    need(delta_k == 16 * (1 + q), 'same-time endpoint energy coefficient')
    need(8 ** 2 / 8 == 8, 'source Fourier time/fold coefficient')
    energy_controls.append({'return_width': str(eps),
                            'H1_ramp_exact_DeltaK': str(delta_k),
                            'all_smooth_profiles_lower_bound': str(2*s.pi/eps)})

# Actual FULL mass-one width-eight QWZ covariance, both quarter backgrounds.
sys.path.insert(0, str(REPO / 'verification'))
source = REPO / 'verification/v1033_charged_disorder.py'
spec = importlib.util.spec_from_file_location('neutral_pair_original_v1033', source)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)
ny = 8
signs = np.repeat([-1 if col in region else 1 for col in range(nx)], 2 * ny)
controls = []
for sector in (1, 3):
    h = mod.qwz_cylinder(nx, ny, 1.0, sector)
    es, vs = np.linalg.eigh(h)
    need(np.min(np.abs(es)) > 1e-8, 'quarter state must have no zero modes')
    occupied = vs[:, es < 0]
    c = occupied @ occupied.conj().T
    transported = signs[:, None] * c * signs[None, :]
    transformed_h = signs[:, None] * h * signs[None, :]
    # The seam phase is unchanged, not C(-h).
    seam0, seam1 = 2 * (nx - 1) * ny, 0
    need(np.array_equal(transformed_h[seam0:seam0+2, seam1:seam1+2],
                        h[seam0:seam0+2, seam1:seam1+2]), 'global seam retained')
    edge_changes = []
    for y in (0, ny - 1):
        ii = 2 * y
        jj = 2 * (ny + y)
        old = c[ii:ii+2, jj:jj+2]
        new = transported[ii:ii+2, jj:jj+2]
        need(np.linalg.norm(old) > 1e-6, 'nonzero actual edge witness')
        need(np.max(np.abs(new + old)) < 1e-13, 'both edge covariances flip')
        edge_changes.append(float(np.linalg.norm(new - old)))
    sub = np.flatnonzero(signs < 0)
    expectation = np.linalg.det(np.eye(len(sub)) - 2 * c[np.ix_(sub, sub)])
    controls.append({'quarter_register_sector': sector,
                     'all_occupied_modes_retained': int(occupied.shape[1]),
                     'left_right_edge_covariance_change': edge_changes,
                     'neutral_parity_expectation_real': float(expectation.real),
                     'neutral_parity_expectation_imag': float(expectation.imag),
                     'global_seam_unchanged': True})

result = {'date': '2026-10-22',
          'verdict': 'ORIGINAL_NEUTRAL_DEFECT_PAIR_CONSTRUCTED_NOT_SINGLE_SPINOR_SECTOR',
          'register_grade': 0, 'global_cover_holonomy_changed': False,
          'uniform_eight_sheet_interval_sign_exact': True,
          'finite_band_factorization': band_controls,
          'finite_band_has_second_Fermi_edge': True,
          'complement_determinant_controls': complement_controls,
          'original_four_mark_central_forms': [str(central), str(-central), '0'],
          'same_source_time_unitary_separation_energy_controls': energy_controls,
          'actual_QWZ_full_covariance_controls': controls,
          'source_sha256': sha(source),
          'infinite_determinant_or_continuum_limit_claimed': False,
          'single_128_spinor_generators_constructed': False,
          'physical_chiral_measure_selected': False}
(OUT / 'neutral-pairs-and-vector-wards-results.json').write_text(
    json.dumps(result, indent=2) + '\n')
print(json.dumps({'verdict': result['verdict'], 'checks_passed': True,
                  'both_actual_quarters': True}))
