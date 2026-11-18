#!/usr/bin/env python3
"""Reconstruct known native TFPT objects against a fixed quantum-information map.

This is a source-comparison probe, not a physical runtime or a new model.
All identities below use bounded Gaussian/integer or rational arithmetic.
"""
import hashlib
import itertools as it
import json
import math
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(REPO))
from tfpt_explorer.process import _code_basis
from tfpt_explorer.cartan_source import _roots, _rays


def need(ok, message):
    if not bool(ok):
        raise RuntimeError(message)


def main():
    _, supports = _code_basis()
    words = list(it.product(range(4), repeat=4))
    index = {word: i for i, word in enumerate(words)}
    pc24 = np.zeros((256, 256), dtype=np.int64)
    for support in supports:
        ids = [index[word] for word in support]
        pc24[np.ix_(ids, ids)] = 24 // len(support)
    ps24 = np.zeros_like(pc24)
    orbits = {}
    for i, word in enumerate(words):
        orbits.setdefault(tuple(sorted(word)), []).append(i)
    for ids in orbits.values():
        ps24[np.ix_(ids, ids)] = 24 // len(ids)

    # Q = (1/16) sum_{a,b in F2^2} W_ab tensor4, in the ORIGINAL frame.
    # W_ab = i^(a.b) X^a Z^b. Its phase disappears in four tensor factors.
    # This signed-permutation formula uses int64 only, no complex tolerance.
    q16 = np.zeros_like(pc24)
    for a, b in it.product(range(4), repeat=2):
        for col, word in enumerate(words):
            row = index[tuple(x ^ a for x in word)]
            parity = sum((b & x).bit_count() for x in word) % 2
            q16[row, col] += 1 if parity == 0 else -1
    need(np.array_equal(q16.T, q16), 'Pauli average is Hermitian')
    need(np.array_equal(q16 @ q16, 16*q16), 'Full Pauli projector')
    need(np.array_equal(q16 @ ps24, 16*pc24), 'Literal symmetric Pauli section is native code')
    need(np.array_equal(ps24 @ q16, 16*pc24), 'Reverse ordered product agrees')
    need(np.trace(q16) == 256, 'Full Pauli fixed rank is 16, not 5')
    need(np.trace(pc24) == 120, 'Symmetric code rank is 5')
    need(np.trace(ps24) == 840, 'Symmetric fourth tensor rank is 35')
    wrong_pc24 = pc24.copy()
    wrong_pc24[0, 0] += 1
    need(not np.array_equal(q16 @ ps24, 16*wrong_pc24), 'Changed code negative control')

    # A standard external non-Clifford algorithmic gate is NOT declared physical.
    # U=diag(1,1,exp(i*pi/4),exp(i*pi/4)); in four slots every code support
    # has even first-bit occupation, so phases are exactly 1,i,-1,-i.
    # The unchanged source code has literal compression diag(0,0,i,i,i).
    nonclifford_compression = []
    for support in supports:
        real = imag = 0
        for word in support:
            count = sum((x >> 1) & 1 for x in word)
            need(count % 2 == 0, 'Even code occupation')
            dr, di = [(1,0),(0,1),(-1,0),(0,-1)][(count//2) % 4]
            real += dr
            imag += di
        nonclifford_compression.append((Fraction(real,len(support)),Fraction(imag,len(support))))
    need(nonclifford_compression == [(0,0),(0,0),(0,1),(0,1),(0,1)], 'Fixed non-Clifford control leaks exactly two code directions')

    # Direct erased-register error detection, independently of a dimension match.
    # Pi_C (|a><b| on slot s) Pi_C = delta_ab Pi_C / 4.
    for slot, a, b in it.product(range(4), range(4), range(4)):
        op = np.zeros_like(pc24)
        for col, word in enumerate(words):
            if word[slot] == b:
                out = list(word)
                out[slot] = a
                op[index[tuple(out)], col] = 1
        target = 6*pc24 if a == b else np.zeros_like(pc24)
        need(np.array_equal(pc24 @ op @ pc24, target), 'Exact one-register Knill Laflamme identity')

    rays = _rays(_roots())
    gram_squared = []
    for x, y in it.product(rays, repeat=2):
        re = sum(x[2*j]*y[2*j] + x[2*j+1]*y[2*j+1] for j in range(4))
        im = sum(x[2*j]*y[2*j+1] - x[2*j+1]*y[2*j] for j in range(4))
        gram_squared.append(re*re + im*im)
    fp = {}
    for t in range(1, 5):
        actual = Fraction(sum(g**t for g in gram_squared), 60**2 * 16**t)
        haar = Fraction(1, math.comb(t+3, 3))
        fp[str(t)] = {'native': str(actual), 'haar': str(haar), 'ratio': str(actual/haar)}
        need((actual == haar) == (t <= 3), 'Exact projective design order')
    need(fp['4']['native'] == '1/32', 'Fourth native frame potential')
    need(fp['4']['ratio'] == '35/32', 'Fourth excess is 35/32')
    # Squared HS norm of M4 - Pi_sym4/35 follows from orthogonal projectors.
    defect = 5*(Fraction(1,20)-Fraction(1,35))**2 + 30*(Fraction(1,40)-Fraction(1,35))**2
    need(defect == Fraction(3,1120), 'Exact fourth-moment defect norm')

    # Equal three-register marginals do NOT imply full-state recoverability.
    # For the known pair rho_C, rho_sym: D before = log(7), after = 0.
    # These exact symbolic values follow from uniform projector spectra.
    output = {
        'verdict': 'SOURCE_COMPARISON_EXACT; KNOWN_STRUCTURE_RECONSTRUCTED; NO_PHYSICAL_LAW_DERIVED',
        'full_pauli_fixed_rank': 16,
        'symmetric_pauli_fixed_rank': 5,
        'symmetric_fourth_tensor_rank': 35,
        'operator_identity': 'Pi_native = Q_Pauli Pi_sym4 = Pi_sym4 Q_Pauli',
        'one_register_erasure': {'matrix_unit_identities': 64, 'all_exact': True, 'meaning': 'Known finite quantum error correction property; no physical noise or recovery operation selected.'},
        'projective_frame_potentials': fp,
        'fourth_moment_defect_squared_hs_norm': str(defect),
        'negative_controls': {'changed_code_rejected': True, 'external_T_tensor4_compression': ['0','0','i','i','i'], 'external_T_leakage_rank': 2, 'physical_permission': 'not claimed; external standard algorithm comparison only'},
        'known_pair_compression': {'D_before': 'log(7)', 'D_after_any_proper_marginal': '0', 'global_trace_distance': '6/7', 'implication': 'This marginal map is not sufficient for this pair; not a no-go for other observations.'},
        'source_sha256': {p: hashlib.sha256((REPO/p).read_bytes()).hexdigest() for p in ['tfpt_explorer/process.py','tfpt_explorer/cartan_source.py']},
        'not_derived': ['simulation ontology','physical source selection','3+1D','native physical 3+2 mark','physical time generator or coupling','Born measurement rule','48*c3^4 or alpha input selection'],
    }
    target = Path(__file__).with_name('algorithmic_source_results.json')
    target.write_text(json.dumps(output, indent=2, sort_keys=True)+'\n')
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
