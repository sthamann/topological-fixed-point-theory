#!/usr/bin/env python3
"""Run the requested finite-matrix check on the original TFPT QWZ builder."""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parent
REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
SOURCE = REPO / 'verification/v1033_charged_disorder.py'
FOLD_RESULTS = REPO / 'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json'

# Import and invoke the source's actual matrix builder; no Hamiltonian is rebuilt here.
sys.path.insert(0, str(REPO / 'verification'))
spec = importlib.util.spec_from_file_location('tfpt_original_v1033', SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError(f'cannot import {SOURCE}')
source = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = source
spec.loader.exec_module(source)

N_VALUES = (8, 16, 32)
NY = 8
MASS = 1.0
SECTORS = (1, 3)
ZERO_TOL = 1e-10
NUMERIC_BOUND_TOL = 1e-9

spectra = []
for sector in SECTORS:
    sector_results = []
    for nx in N_VALUES:
        h = source.qwz_cylinder(nx=nx, ny=NY, mass=MASS, sector=sector)
        eigenvalues = np.linalg.eigvalsh(h)  # full 2*N*8 source matrix
        minimum = float(np.min(np.abs(eigenvalues)))
        scaled = nx * minimum / (2 * math.pi)
        analytic_finite_bound = nx * math.sin(math.pi / (2 * nx)) / (2 * math.pi)
        negative_count = int(np.count_nonzero(eigenvalues < -ZERO_TOL))
        zero_count = int(np.count_nonzero(np.abs(eigenvalues) <= ZERO_TOL))
        passed = minimum + NUMERIC_BOUND_TOL >= math.sin(math.pi / (2 * nx))
        if not passed or zero_count != 0 or negative_count != nx * NY:
            raise AssertionError(
                f'failed sector={sector}, N={nx}: min={minimum}, bound='
                f'{math.sin(math.pi/(2*nx))}, neg={negative_count}, zero={zero_count}'
            )
        sector_results.append({
            'N': nx,
            'matrix_dimension': int(h.shape[0]),
            'min_abs_eigenvalue': minimum,
            'scaled_N_over_2pi_min_abs_eigenvalue': scaled,
            'analytic_scaled_finite_bound_N_sin_pi_over_2N_over_2pi': analytic_finite_bound,
            'negative_eigenvalue_count': negative_count,
            'zero_count_at_abs_eigenvalue_le_1e-10': zero_count,
            'numeric_minimum_ge_analytic_bound_with_tolerance_1e-9': passed,
        })
    spectra.append({'sector': sector, 'seam_phase_i_to_sector': ['1', 'i', '-1', '-i'][sector % 4], 'cases': sector_results})

fold = json.loads(FOLD_RESULTS.read_text(encoding='utf-8'))
M = [[Fraction(x) for x in row] for row in fold['explicit_Weyl_matrix']]
p_star = [Fraction(1, 2)] * 4 + [Fraction(-1, 2)] * 4
image = [sum(M[i][j] * p_star[j] for j in range(8)) for i in range(8)]
source_charge = [Fraction(x) for x in fold['transported_source_U1_charge']]
dot = sum(a * b for a, b in zip(image, source_charge))
expected_image = [Fraction(1), Fraction(0), Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(0), Fraction(0)]
weyl_pass = image == expected_image and dot == 0
if not weyl_pass:
    raise AssertionError(f'Weyl transport mismatch: image={image}, dot={dot}')

result = {
    'verdict': 'PASS: requested microscopic checks pass; this is not a complete P1/P2 derivation.',
    'assumptions': [
        'Call the original verification/v1033_charged_disorder.py qwz_cylinder with mass=1, ny=8, sectors=1 and 3.',
        'Diagonalize the full original 2*N*ny matrix for N=8,16,32; numerical zero count uses |eigenvalue| <= 1e-10.',
        'The numerical comparison tolerance against the exact analytic finite-N lower bound is 1e-9.',
        'The analytic gap statement is scoped to the stated noninteracting source QWZ cylinder and quarter-twist sectors.',
        'Weyl transport uses the stipulated folded E8 matrix and source Cartan charge already recorded in fold-results.json.'
    ],
    'categories': {
        'original_source_operator_spectrum': 'numerically_checked_full_finite_matrices',
        'source_noninteracting_gap': 'analytic_from_exact_block_identity; numeric matrices are consistency checks',
        'folded_Weyl_transport': 'exact_rational_check_within_stipulated_folded_extension',
        'complete_TFPT_solution': False,
        'raw_P1P2_derivation': False,
        'raw_quotient_or_E8_extension_selected_by_original_axioms': False
    },
    'checks': {
        'original_QWZ_full_matrix_spectra': {
            'source_function': 'qwz_cylinder',
            'source_path': str(SOURCE),
            'source_line_range': '145-166',
            'matrix_dimension_formula': '2*N*8',
            'all_cases_pass': True,
            'cases': spectra,
        },
        'analytic_source_gap': {
            'source_identity_path': str(REPO / 'experiments/theory-contracts/source-rg-clock-bridge-20260920/checker.py'),
            'source_identity_lines': '86-104',
            'identity': 'H(k)=[[-s I,B*],[B,s I]], H(k)^2=diag(s^2 I+B*B,s^2 I+BB*) >= sin(k)^2 I',
            'quarter_twist_momenta': 'k=2*pi*(j +/- 1/4)/N',
            'minimum_absolute_sine': 'sin(pi/(2N))',
            'analytic_finite_N_gap_bound': 'min|lambda| >= sin(pi/(2N))',
            'scaled_lower_bound': 'N/(2*pi)*sin(pi/(2N))',
            'uniform_scaled_lower_bound': '>= 1/(2*pi)',
            'limit': '1/4',
            'scope': 'the specified source noninteracting one-particle/cylinder operator; not an all-P1/P2 no-go or full interacting theory result',
            'numeric_matrix_check_passes_bound': True,
            'zero_absence_is_analytical': True,
        },
        'folded_E8_Weyl_transport': {
            'fold_results_path': str(FOLD_RESULTS),
            'p_star': ['1/2', '1/2', '1/2', '1/2', '-1/2', '-1/2', '-1/2', '-1/2'],
            'M_times_p_star': [str(x) for x in image],
            'expected_image': ['1', '0', '0', '1', '0', '0', '0', '0'],
            'transported_source_U1_charge': [str(x) for x in source_charge],
            'image_dot_source_charge': str(dot),
            'norm_squared_of_image': str(sum(x*x for x in image)),
            'D5_carrier_root_reason': 'image=e1+e4 has two nonzero integer coordinates within the first five coordinates and norm squared 2; hence it is an ordinary D5 root in the stipulated coordinate embedding.',
            'exact_rational_check_passes': weyl_pass,
            'scope': 'conditional on the stipulated folded E8 Weyl map and its embedding; does not derive that extension from raw P1/P2.'
        }
    }
}
(WORK / 'microscopic-results.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps({
    'verdict': result['verdict'],
    'spectra': [{
        'sector': block['sector'],
        'cases': [{k: case[k] for k in ('N', 'matrix_dimension', 'min_abs_eigenvalue', 'scaled_N_over_2pi_min_abs_eigenvalue', 'negative_eigenvalue_count', 'zero_count_at_abs_eigenvalue_le_1e-10')} for case in block['cases']]
    } for block in spectra],
    'weyl': result['checks']['folded_E8_Weyl_transport'],
}, indent=2))
