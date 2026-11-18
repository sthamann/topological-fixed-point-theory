"""Direct current-source spectral controls for the written all-energy argument.

Finite numerical spectra corroborate, but do not prove, the uniform theorem.
No historical source pin is overwritten or bypassed with a fabricated PASS.
"""
import argparse
import hashlib
import importlib
import itertools
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import sympy as S


def need(condition, message):
    if not condition:
        raise ValueError(message)


def exact(matrix):
    return S.Matrix([[S.Rational(float(z.real)) + S.I*S.Rational(float(z.imag))
                      for z in row] for row in matrix])


def count_fock(energies, cutoff):
    totals = [0.0]
    for energy in energies:
        totals += [total + energy for total in totals if total + energy < cutoff]
    return len(totals)


def run(repo):
    sys.path.insert(0, str(repo / 'verification'))
    source = importlib.import_module('v1033_charged_disorder')
    h0n = source.qwz_cylinder(1, 8, 1.0, 0)
    h0 = exact(h0n)
    x = S.diag(*([exact(source.SX)]*8))
    z = S.diag(*([exact(source.SZ)]*8))
    ident = S.eye(16)
    top, bottom = S.zeros(16,1), S.zeros(16,1)
    top[14] = top[15] = 1/S.sqrt(2)
    bottom[0], bottom[1] = 1/S.sqrt(2), -1/S.sqrt(2)
    edge = S.Matrix.hstack(top,bottom)
    identities = {
        'H0_cubic': h0**3 == h0,
        'H0_rank14': h0.rank() == 14,
        'X_square': x*x == ident,
        'Z_square': z*z == ident,
        'X_H0_anticommutator': x*h0+h0*x == S.zeros(16),
        'X_Z_anticommutator': x*z+z*x == S.zeros(16),
        'pi_gap': (h0+2*z).det() != 0,
        'edge_orthonormality': edge.T*edge == S.eye(2),
        'edge_zero_vectors': h0*edge == S.zeros(16,2),
        'opposite_edge_slopes': edge.T*(-x)*edge == S.diag(-1,1),
    }
    need(all(identities.values()), 'actual source does not satisfy theorem premises')
    tx = np.kron(np.eye(8), source.TX)

    def bloch(p):
        return h0n + (np.exp(-1j*p)-1)*tx + (np.exp(1j*p)-1)*tx.conj().T

    # Compare with the actual full-cylinder builder, not just our Bloch expression.
    full_controls = []
    for sector in (1, 3):
        n = 8
        blocks = np.sort(np.concatenate([
            np.linalg.eigvalsh(bloch(2*np.pi*(j-sector/4)/n)) for j in range(n)]))
        direct = np.linalg.eigvalsh(source.qwz_cylinder(n, 8, 1.0, sector))
        error = float(np.max(abs(blocks-direct)))
        need(error < 1e-12, 'Bloch blocks disagree with current full source')
        full_controls.append({'sector': sector, 'N': n, 'max_error': error})

    # A grid is a diagnostic of the analytical three-region bound, not its proof.
    bulk_min = 10.0
    max_subhalf = 0
    for p in np.linspace(-np.pi, np.pi, 513):
        values = np.sort(abs(np.linalg.eigvalsh(bloch(p))))
        bulk_min = min(bulk_min, float(values[2]))
        max_subhalf = max(max_subhalf, int(np.sum(values < .5-1e-12)))
    need(bulk_min >= .5-1e-12 and max_subhalf <= 2, 'uniform-gap grid control')

    spectral = []
    for sector in (1, 3):
        alpha = sector/4
        for n in (64, 128, 256):
            actual = np.sort(np.concatenate([
                abs(n*np.linalg.eigvalsh(bloch(2*np.pi*(j-alpha)/n))/(2*np.pi))
                for j in range(n)]))
            for cutoff in (2.13, 3.13):
                low = actual[actual < cutoff]
                limit = sorted(abs(j-alpha) for j in range(-10, 11)
                               for _ in range(2) if abs(j-alpha) < cutoff)
                need(len(low) == len(limit), 'one-body spectral rank mismatch')
                finite_fock = count_fock(low, cutoff)
                limiting_fock = count_fock(limit, cutoff)
                need(finite_fock == limiting_fock, 'many-body spectral rank mismatch')
                spectral.append({'sector': sector, 'N': n, 'cutoff': cutoff,
                                 'one_body_rank': len(low),
                                 'fock_spectral_rank': finite_fock,
                                 'max_energy_error': float(max(abs(low-limit)))})

    # Exhaustive norm-two census: |m_i+1/4| <= sqrt(2) permits only -1,0,1.
    quarter = []
    for m in itertools.product((-1, 0, 1), repeat=8):
        q = tuple(F(1, 4)+a for a in m)
        if sum(a*a for a in q) == 2:
            quarter.append((m, q))
    need(len(quarter) == 64, 'quarter norm-two count')
    need(all(sum(m) % 2 == 1 for m, q in quarter), 'relative odd dressing')
    need(all((q[0]+q[1]).denominator == 2 for m, q in quarter),
         'nontrivial monodromy against the fixed D8 root e1+e2')
    # All four scalar Z4 phases, tested against all D8 roots: only 0 and 1/2 local.
    scalar_locality = {str(F(r, 4)): r % 2 == 0 for r in range(4)}
    for r in range(4):
        pairings = [F(r,4)*(a+b) for a,b in itertools.product((-1,1),repeat=2)]
        need(all(v.denominator == 1 for v in pairings) == scalar_locality[str(F(r,4))],
             'scalar Z4 locality classification')

    half_twist_time = []
    for alpha in (F(1,4), F(3,4)):
        delta = [(F(k)+alpha)/8-F(1,2) for k in range(8)]
        for chirality in (0,1):
            values = [1+sum(F(s,2)*d for s,d in zip(signs,delta))
                      for signs in itertools.product((-1,1),repeat=8)
                      if sum(s<0 for s in signs)%2 == chirality]
            minimum = min(values)
            need(minimum == (F(0) if chirality==0 else F(1,32)),
                 'actual common-time spinor-chirality minimum')
            need(values.count(minimum)==1, 'spinor minimum multiplicity')
            half_twist_time.append({'alpha':str(alpha),'minus_parity':chirality,
                                   'minimum_fold_energy':str(minimum),'minimum_multiplicity':1})

    return {
        'date': '2026-10-22',
        'verdict': 'PASS_TARGETED_CONTROLS',
        'analytical_theorem': 'FULL-FINITE-ENERGY-AND-SECTOR-MENU.txt',
        'source_sha256': hashlib.sha256((repo/'verification/v1033_charged_disorder.py').read_bytes()).hexdigest(),
        'exact_source_identities': identities,
        'full_builder_controls': full_controls,
        'gap_grid_diagnostic': {'minimum_third_absolute_eigenvalue': bulk_min,
                                'maximum_subhalf_eigenvalue_count': max_subhalf},
        'spectral_controls': spectral,
        'quarter_h1_count': len(quarter),
        'quarter_h1_all_odd_relative_dressing': True,
        'quarter_h1_all_nonlocal_to_full_D8': True,
        'D8_relative_locality_of_uniform_scalar_Z4_classes': scalar_locality,
        'half_twist_common_time_by_chirality': half_twist_time,
        'source_vacuum_and_time_retained': True,
        'new_P1_spinor_sector_constructed': False,
        'full_P1_or_TFPT_exclusion': False,
        'historical_pinned_checker_repaired_or_replayed': False,
        'verification_boundary': 'Exact matrix/arithmetic checks plus finite floating spectral controls. The all-N and compactness claims rest on the written proof, not this finite run.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'verdict': result['verdict'], 'source_identities': len(result['exact_source_identities']),
                      'full_builder_controls': len(result['full_builder_controls']),
                      'spectral_cases': len(result['spectral_controls']), 'quarter_fields': result['quarter_h1_count']}))
