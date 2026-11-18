#!/usr/bin/env python3
"""Two processes on the actual six TFPT channel words, not two full TFPT models.

Exact arithmetic certifies a scoped inference counterexample. The all-time
statement is proved in PROOF.md; finite checks are independent diagnostics.
"""
import importlib.util
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    import sympy as s
    path = ROOT / 'experiments/theory-contracts/compiler-clifford-bridge/checker.py'
    spec = importlib.util.spec_from_file_location('history_bridge', path)
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    src, _ = bridge.inherited(ROOT)
    finite = bridge.finite_data(src)
    zero = (0, 0, 0, 0)
    five = [v for v in src.W16 if v != zero and finite['q'][v] == 0]
    five = [next(v for v in five if src.iota(v)[j] == 0) for j in range(5)]
    labels = [zero] + five
    words = [bridge.monomial(v, bridge.generators()) for v in labels]
    p = [F(7, 12)] + [F(1, 12)] * 5
    eps = F(1, 288)
    triples = list(itertools.permutations(range(1, 6), 3))

    def xor(a, b):
        return tuple(x ^ y for x, y in zip(a, b))

    def delta(triple):
        d = [[0] * 6 for _ in range(6)]
        for a, b in zip(triple, triple[1:] + triple[:1]):
            d[a][b] = 1
            d[b][a] = -1
        return d

    Ds = [delta(t) for t in triples]
    Js = [[[p[a] * p[b] + eps*d[a][b] for b in range(6)]
           for a in range(6)] for d in Ds]
    q2 = {v: F(0) for v in src.W16}
    for a in range(6):
        for b in range(6):
            q2[xor(labels[a], labels[b])] += p[a]*p[b]

    for J in Js:
        require(all(x >= 0 for row in J for x in row), 'negative probability')
        require([sum(row) for row in J] == p, 'row marginal')
        require([sum(J[a][b] for a in range(6)) for b in range(6)] == p, 'column marginal')
        q = {v: F(0) for v in src.W16}
        for a in range(6):
            for b in range(6):
                q[xor(labels[a], labels[b])] += J[a][b]
                U = words[b]*words[a]
                W = bridge.monomial(xor(labels[a], labels[b]), bridge.generators())
                phase = s.simplify(s.trace(W.H*U)/4)
                require(U == phase*W and phase in (1, -1, s.I, -s.I), 'word action')
        require(q == q2, 'different endpoint distribution')
    require(all(sum(d[a][b] for d in Ds) == 0 for a in range(6) for b in range(6)), 'nonzero mean perturbation')

    def probability(history):
        """Stationary mixture of the two dimer offsets and shared hidden frame."""
        out = F(0)
        for J in Js:
            for offset in (0, 1):
                weight = F(1)
                i = 0
                if offset and history:
                    weight *= p[history[0]]
                    i = 1
                while i + 1 < len(history):
                    weight *= J[history[i]][history[i+1]]
                    i += 2
                if i < len(history):
                    weight *= p[history[i]]
                out += weight
        return out / (2*len(Js))

    # All three-step histories agree; the four-step witness does not.
    from math import prod
    for n in (1, 2, 3):
        for h in itertools.product(range(6), repeat=n):
            require(probability(h) == prod(p[a] for a in h), 'short-time mismatch')
    h = (1, 2, 1, 2)
    iid = prod(p[a] for a in h)
    correlated = probability(h)
    require(iid == F(1,20736), 'iid witness')
    require(correlated - iid == F(1,552960), 'correlated witness')
    # Exact extension consistency around the witness is a check of stationarity.
    require(sum(probability((a,) + h) for a in range(6)) == correlated, 'left extension')
    require(sum(probability(h + (a,)) for a in range(6)) == correlated, 'right extension')
    # Full S5 invariance follows by permuting the complete uniform triple set.
    for perm in itertools.permutations(range(1,6)):
        trans = (0,) + perm
        require(probability(tuple(trans[a] for a in h)) == correlated, 'S5 witness covariance')

    # Bell-state probes identify each of the six original unitary operations.
    gram = s.Matrix(6,6, lambda a,b: s.trace(words[a].H*words[b])/4)
    require(gram == s.eye(6), 'Bell label observability')

    # Stronger version: random block lengths 1/2, and mixing hidden frames.
    # This removes both the global conserved frame and periodic dimer phase.
    kappa = F(1,2)
    def mixing_probability(history):
        R = [F(1,90)] * 60
        T = [[p[a]/180 for a in range(6)] for _ in range(60)]
        for a in history:
            endings = [R[h]*p[a]/2 + sum(T[h][b]*Js[h][b][a]/p[b] for b in range(6)) for h in range(60)]
            total = sum(endings)
            newT = [[F(0)]*6 for _ in range(60)]
            for h in range(60):
                newT[h][a] = R[h]*p[a]/2
            R = [kappa*x + (1-kappa)*total/60 for x in endings]
            T = newT
        return sum(R) + sum(sum(row) for row in T)

    require(mixing_probability(()) == 1, 'mixing initial normalization')
    # Sum of all emitted symbols preserves the exact stationary hidden masses.
    Rstar = F(1,90)
    for hidx in range(60):
        ending_mass = Rstar/2 + sum(p[a]/180 for a in range(6))
        require(ending_mass == Rstar, 'mixing ready stationary mass')
        require(all(Rstar*p[a]/2 == p[a]/180 for a in range(6)), 'mixing pending stationary mass')
    for n in (1,2,3):
        for hist in itertools.product(range(6), repeat=n):
            require(mixing_probability(hist) == prod(p[a] for a in hist), 'mixing short-word mismatch')
    mixing_witness = mixing_probability(h)
    require(mixing_witness - iid == F(1,3317760), 'mixing four-step witness')
    require(sum(mixing_probability((a,)+h) for a in range(6)) == mixing_witness, 'mixing left extension')
    require(sum(mixing_probability(h+(a,)) for a in range(6)) == mixing_witness, 'mixing right extension')
    result = {
        'scope': 'Counterexample to identifying a multitime process from Phi^n, stationary I4/4 and signed marked S5 covariance; not two complete TFPT realizations.',
        'same_original_five_words': [list(v) for v in five],
        'p': [str(x) for x in p], 'epsilon': str(eps),
        'hidden_frames': len(triples),
        'all_frame_pair_laws_nonnegative_with_same_marginals_and_q2': True,
        'iid_and_correlated_all_unperturbed_interval_channels': 'exactly Phi^n for every n>=0; analytic block proof',
        'both_stationary_and_marked_S5_covariant': True,
        'all_histories_of_lengths_1_2_3_identical': True,
        'four_step_probe': {'labels': list(h), 'iid':str(iid), 'correlated':str(correlated), 'difference':str(correlated-iid)},
        'four_step_probability_left_right_extension_consistent': True,
        'six_Bell_error_states_are_orthogonal': True,
        'mixing_ergodic_version': {
            'hidden_states': 420,
            'block_lengths': {'1':'1/2','2':'1/2'},
            'frame_transition': 'I/2 + Uniform60/2',
            'stationary_hidden_masses': {'each_ready':'1/90','pending_h_a':'p[a]/180'},
            'irreducible_aperiodic': 'all frame resets and emissions positive; each ready state has a positive self-loop',
            'all_histories_of_lengths_1_2_3_identical': True,
            'all_unperturbed_interval_channels': 'Phi^n for every n by the same conditional block proof',
            'four_step_iid': str(iid),
            'four_step_correlated': str(mixing_witness),
            'four_step_difference': str(mixing_witness-iid),
            'left_right_extension_consistent': True,
        },
        'physical_gates_closed': [],
    }
    Path(__file__).with_name('history-identifiability-result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
