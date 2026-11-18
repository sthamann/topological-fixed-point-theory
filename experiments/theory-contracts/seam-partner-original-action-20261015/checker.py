"""Exact finite carrier audit; no physical vacuum/time identification is assumed."""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def creation(j, modes=5):
    out = np.zeros((2**modes, 2**modes), dtype=np.int64)
    for n in range(2**modes):
        if not (n >> j) & 1:
            out[n | (1 << j), n] = (-1) ** ((n & ((1 << j) - 1)).bit_count())
    return out


def check():
    identity = np.eye(32, dtype=np.int64)
    eps = [creation(j) for j in range(5)]
    iota = [e.T for e in eps]
    number = [eps[j] @ iota[j] for j in range(5)]
    even = [n for n in range(32) if n.bit_count() % 2 == 0]
    odd = [n for n in range(32) if n.bit_count() % 2 == 1]
    checks = {}

    for i in range(5):
        for j in range(5):
            assert np.array_equal(iota[i] @ eps[j] + eps[j] @ iota[i], identity * (i == j))
    checks['CAR_matrix_relations'] = 25

    Y6 = sum(q * n for q, n in zip([-2, -2, -2, 3, 3], number))
    T3_2 = number[3] - number[4]
    Qem6 = Y6 + 3 * T3_2
    Qcarrier = 2 * number[2]
    pair = eps[0] @ eps[3]
    empty = identity[:, 0]
    result = pair @ empty
    assert np.array_equal(result, identity[:, 9])
    for charge, multiple in [(Y6, 1), (Qem6, 4), (Qcarrier, 0)]:
        assert np.array_equal(charge @ pair - pair @ charge, multiple * pair)
    checks['pair_charge_commutators'] = 3

    # Real and imaginary integer entries are exactly represented here; equality
    # uses no floating tolerance. All intermediate magnitudes are far below 2**53.
    hermitian_basis = []
    for i in range(5):
        h = np.zeros((5, 5), dtype=np.complex128)
        h[i, i] = 1
        hermitian_basis.append(h)
    for i, j in combinations(range(5), 2):
        for z in (1, 1j):
            h = np.zeros((5, 5), dtype=np.complex128)
            h[i, j], h[j, i] = z, np.conjugate(z)
            hermitian_basis.append(h)

    def lift(h):
        return sum(h[j, k] * eps[j] @ iota[k] for j in range(5) for k in range(5))

    B = identity.copy()
    for j in range(5):
        B = B @ (eps[j] + iota[j])
    assert np.array_equal(B @ B, identity)
    parity = np.diag([(-1) ** n.bit_count() for n in range(32)])
    assert np.array_equal(B @ parity @ B, -parity)
    assert np.array_equal(B @ Y6 @ B, -Y6)
    assert np.array_equal(B @ empty, identity[:, 31])
    for h in hermitian_basis:
        H = lift(h)
        assert np.array_equal(B @ H.conjugate() @ B, np.trace(h) * identity - H)
        for i, a in combinations(range(5), 2):
            A = eps[i] @ eps[a]
            rhs = sum(h[j, i] * eps[j] @ eps[a] + h[j, a] * eps[i] @ eps[j] for j in range(5))
            assert np.array_equal(H @ A - A @ H, rhs)
    checks['time_commutator_basis_cases'] = 250
    checks['antiunitary_charge_conjugation_basis_cases'] = 25

    zero_even = [n for n in even if Y6[n, n] == 0]
    zero_odd = [n for n in odd if Y6[n, n] == 0]
    assert zero_even == [0] and zero_odd == [31]
    # The carrier involution is an explicit representation-level construction;
    # equality with the geometric seam reflection is NOT checked or assumed.

    fold_path = ROOT / 'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json'
    recorded = json.loads(fold_path.read_text())
    M = [[F(v) for v in row] for row in recorded['explicit_Weyl_matrix']]
    assert all(sum(M[i][k] * M[j][k] for k in range(8)) == (i == j) for i in range(8) for j in range(8))

    def transform(v):
        return [sum(a * b for a, b in zip(row, v)) for row in M]

    r0 = [F(1, 2)] * 4 + [F(-1, 2)] * 4
    root = transform(r0)
    assert root == [1, 0, 0, 1, 0, 0, 0, 0]
    assert transform([F(1)] * 8) == [0, 0, 2, 0, 0, 2, 0, 0]
    spectra = []
    for alpha in (F(1, 4), F(3, 4)):
        delta = [(j + alpha) / 8 - F(1, 2) for j in range(8)]
        transported = transform(delta)
        expected = [F(-1, 2), 0, alpha / 4, F(-1, 2), 0, (alpha - 1) / 4, F(-1, 4), F(-1, 4)]
        assert transported == expected
        energies = {f'e{i+1}+e{a+1}': 1 + transported[i] + transported[a] for i in range(3) for a in (3, 4)}
        assert list(energies.values()) == [0, F(1, 2), F(1, 2), 1, F(1, 2) + alpha / 4, 1 + alpha / 4]
        # The root operators form (3,2), but Q_source is not constant on it.
        q_values = {f'e{i+1}+e{a+1}': 2 * (i == 2) for i in range(3) for a in (3, 4)}
        ward = [transported[i] - transported[j] for i, j in [(0, 1), (1, 2), (3, 4)]]
        assert ward == [F(-1, 2), -alpha / 4, F(-1, 2)]
        spectra.append({'alpha': str(alpha), 'transported_delta': list(map(str, transported)),
                        'six_root_energies': {k: str(v) for k, v in energies.items()},
                        'six_source_charges': q_values,
                        'fixed_SM_simple_root_time_commutators': list(map(str, ward))})
    checks['quarter_spectra'] = 2

    # Actual one-channel source modes; this compares opposite quarter domains,
    # not two vacua of the same domain and not a new lattice Hamiltonian.
    quarter_pairs = []
    # A finite exact Fock implementation checks the ordering signs of the same
    # mode map, rather than inferring CAR intertwining from matched energies.
    ns = [-2, -1, 0, 1]
    fplus = [creation(j, 4) for j in range(4)]
    fnumber = [a @ a.T for a in fplus]
    U = np.zeros((16, 16), dtype=np.int64)
    for mask in range(16):
        indices = [j for j in range(4) if (mask >> j) & 1]
        mapped = [3 - j for j in indices]
        inversions = sum(mapped[j] > mapped[k] for j in range(len(mapped)) for k in range(j + 1, len(mapped)))
        U[sum(1 << j for j in mapped), mask] = (-1)**inversions
    assert np.array_equal(U @ U.T, np.eye(16, dtype=np.int64))
    assert U[0, 0] == 1
    for alpha in (F(1, 4), F(3, 4)):
        other = 1 - alpha
        for n in range(-16, 17):
            partner_n = -n - 1
            assert partner_n + other == -(n + alpha)
            assert abs(partner_n + other) == abs(n + alpha)
            assert ((n + alpha) < 0) == ((partner_n + other) > 0)
        field = [fplus[j].T if n + alpha > 0 else fplus[j] for j, n in enumerate(ns)]
        other_field = [fplus[j].T if n + other > 0 else fplus[j] for j, n in enumerate(ns)]
        for j in range(4):
            assert np.array_equal(U @ field[j] @ U.T, other_field[3-j].T)
        H4 = sum(int(4 * abs(n + alpha)) * fnumber[j] for j, n in enumerate(ns))
        other_H4 = sum(int(4 * abs(n + other)) * fnumber[j] for j, n in enumerate(ns))
        charge = sum((1 if n + alpha > 0 else -1) * fnumber[j] for j, n in enumerate(ns))
        other_charge = sum((1 if n + other > 0 else -1) * fnumber[j] for j, n in enumerate(ns))
        assert np.array_equal(U @ H4 @ U.T, other_H4)
        assert np.array_equal(U @ charge @ U.T, -other_charge)
        quarter_pairs.append({'alpha': str(alpha), 'conjugate_alpha': str(other),
                              'mode_map': 'n -> -n-1', 'physical_gap': '1/4',
                              'new_same_sector_zero_state': False})
    checks['original_chiral_quarter_mode_pairs'] = 66
    checks['finite_polarized_Fock_CAR_intertwining'] = 8
    checks['finite_polarized_Fock_charge_and_energy_intertwining'] = 4

    return {
        'research_id': 'UR.SEAM.PARTNER_ORIGINAL_ACTION.01',
        'verdict': 'EXACT_TYPED_CARRIER_AND_FOLDED_TIME_AUDIT',
        'checks': checks,
        'carrier': {'dimension_even': len(even), 'dimension_odd': len(odd),
                    'root_action': 'epsilon_1 epsilon_4 . 1 = u_1 wedge u_4',
                    'pair_hypercharge': '1/6', 'pair_electric_charge_T3_plus_convention': '2/3',
                    'pair_transported_source_charge_shift': 0,
                    'only_hypercharge_zero_even_occupation': zero_even,
                    'only_hypercharge_zero_odd_occupation': zero_odd,
                    'antiunitary': 'J=(product_j(epsilon_j+iota_j)) complex_conjugation',
                    'time_identity': 'J dGamma(h) J^-1 = tr(h) I - dGamma(h)'},
        'folded_current_module': spectra,
        'original_chiral_conjugation': quarter_pairs,
        'source_fold_sha256': hashlib.sha256(fold_path.read_bytes()).hexdigest(),
        'raw_P1P2_derivation': False,
        'physical_vacuum_identified_with_exterior_unit': False,
        'source_seam_reflection_identified_with_carrier_J': False,
        'quarter_phase_erased': False,
        'each_quarter_domain_retained': True,
        'conjugation_exchanges_quarter_domains': True,
        'physical_time_selected': False,
        'complete_TFPT_solution': False,
        'scope': 'Exact finite carrier action plus already stipulated E8 fold. No full source countermodel, origin proof, physical state selection or gauge-invariant-time assumption is inferred.'
    }


if __name__ == '__main__':
    result = check()
    (HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
