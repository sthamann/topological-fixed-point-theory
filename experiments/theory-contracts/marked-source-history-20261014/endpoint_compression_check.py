#!/usr/bin/env python3
"""Exact endpoint compression and open-slot obstruction for repeated Phi."""
import importlib.util
import json
from pathlib import Path

import sympy as s

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def xor(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def ad(u):
    return s.kronecker_product(u.conjugate(), u)


def vec(u):
    return s.Matrix([u[i, j] / 2 for j in range(4) for i in range(4)])


def main():
    bridge = load('bridge_endpoint', ROOT / 'experiments/theory-contracts/compiler-clifford-bridge/checker.py')
    src, _ = bridge.inherited(ROOT)
    finite = bridge.finite_data(src)
    words = {v: bridge.monomial(v, bridge.generators()) for v in src.W16}
    zero = (0, 0, 0, 0)
    o5 = [v for v in src.W16 if v != zero and finite['q'][v] == 0]
    o10 = [v for v in src.W16 if finite['q'][v] == 1]
    five = [next(v for v in o5 if src.iota(v)[j] == 0) for j in range(5)]

    p = {v: s.Integer(0) for v in src.W16}
    p[zero] = s.Rational(7, 12)
    for v in five:
        p[v] = s.Rational(1, 12)

    def conv(a, b):
        return {v: s.simplify(sum(a[u] * b[xor(u, v)] for u in src.W16)) for v in src.W16}

    def formula(n, v):
        a = s.Rational(1, 3) ** n
        b = s.Rational(2, 3) ** n
        if v == zero:
            return s.simplify((1 + 5*a + 10*b) / 16)
        if v in o5:
            return s.simplify((1 - 3*a + 2*b) / 16)
        return s.simplify((1 + a - 2*b) / 16)

    distributions = {1: p}
    cur = p
    for n in range(2, 9):
        cur = conv(cur, p)
        distributions[n] = cur
    formula_exact = all(distributions[n][v] == formula(n, v) for n in range(1, 9) for v in src.W16)

    ranks = {}
    choi_rank_direct = {}
    for n in range(1, 9):
        qn = distributions[n]
        ranks[n] = sum(1 for x in qn.values() if x > 0)
        C = sum((qn[v] * vec(words[v]) * vec(words[v]).H for v in src.W16), s.zeros(16))
        choi_rank_direct[n] = C.rank()

    # Direct two-step collapse, including the projective Pauli phase.
    labels = [zero] + five
    histories_by_endpoint = {v: [] for v in src.W16}
    for first in labels:
        for second in labels:
            endpoint = xor(first, second)
            product = s.simplify(words[second] * words[first])
            phase = s.simplify(s.trace(words[endpoint].H * product) / 4)
            assert phase in (1, -1, s.I, -s.I)
            assert product == phase * words[endpoint]
            histories_by_endpoint[endpoint].append({
                'history': [list(first), list(second)],
                'weight': str(p[first] * p[second]),
                'phase': str(phase),
            })
    two_step_endpoint_norms = {
        str(v): s.simplify(sum(p[tuple(h['history'][0])] * p[tuple(h['history'][1])] for h in hs))
        for v, hs in histories_by_endpoint.items()
    }
    collapse_norms_exact = all(two_step_endpoint_norms[str(v)] == distributions[2][v] for v in src.W16)
    direct_two_step = sum((p[a] * p[b] * ad(words[b] * words[a]) for a in labels for b in labels), s.zeros(16))
    endpoint_two_step = sum((distributions[2][v] * ad(words[v]) for v in src.W16), s.zeros(16))

    # Same endpoint without an insertion, distinct open-slot histories with one.
    vi, vk = five[0], five[1]
    gi, gk = words[vi], words[vk]
    assert gi * gk == -gk * gi
    A = (s.eye(4) + s.I * gk) / s.sqrt(2)
    assert s.simplify(A.H * A - s.eye(4)) == s.zeros(4)
    left = s.simplify(gi * A)   # history (0,i): second gi after insertion
    right = s.simplify(A * gi)  # history (i,0): insertion after first gi
    hs_inner = s.simplify(s.trace(left.H * right) / 4)
    open_slot_maps_distinct = s.simplify(ad(left) - ad(right)) != s.zeros(16)

    q2 = distributions[2]
    result = {
        'premise': 'identical independently repeated uses of the existing one-step channel Phi',
        'actual_five_labels': [list(v) for v in five],
        'general_endpoint_distribution': {
            'identity': '(1 + 5*(1/3)^n + 10*(2/3)^n)/16',
            'each_of_five': '(1 - 3*(1/3)^n + 2*(2/3)^n)/16',
            'each_of_ten': '(1 + (1/3)^n - 2*(2/3)^n)/16',
        },
        'formula_verified_by_exact_convolution_n_1_to_8': bool(formula_exact),
        'q2': {
            'identity': str(q2[zero]),
            'each_of_five': str(q2[o5[0]]),
            'each_of_ten': str(q2[o10[0]]),
            'normalization': str(sum(q2.values())),
        },
        'choi_ranks_n_1_to_8_from_positive_endpoints': ranks,
        'choi_ranks_n_1_to_8_direct': choi_rank_direct,
        'two_step_projective_history_collapse_exact': bool(collapse_norms_exact and direct_two_step == endpoint_two_step),
        'two_step_number_of_histories': sum(map(len, histories_by_endpoint.values())),
        'two_step_number_of_endpoint_classes': sum(bool(x) for x in histories_by_endpoint.values()),
        'endpoint_coherent_dilation': 'V_n psi = sum_v sqrt(q_n(v)) W_v psi tensor |e_v^(n)>',
        'endpoint_vector_definition': '|e_v^(n)> = q_n(v)^(-1/2) sum_{h:end(h)=v} sqrt(p_h) omega(h)|h>',
        'endpoint_vectors_orthogonal_reason': 'disjoint history support for distinct v',
        'mark_transport_identity': 'P_v=W_v P W_v^*, hence (sum_v P_v tensor |e_v><e_v|) V_n P = V_n P',
        'open_slot_witness': {
            'endpoint_label': list(vi),
            'anticommuting_insert_label': list(vk),
            'A': '(I + i*g_k)/sqrt(2)',
            'history_0_i_operator': 'g_i A',
            'history_i_0_operator': 'A g_i',
            'normalized_HS_inner_product': str(hs_inner),
            'superoperators_distinct': bool(open_slot_maps_distinct),
        },
        'scope': [
            'Endpoint compression is exact for the terminal, unperturbed composite channel Phi^n.',
            'It does not retain the open slots of the repeated process and cannot predict arbitrary intermediate interventions by itself.',
            'Independent repetitions are an explicit extra operational premise, not a P1/P2 derivation of Phi or of source execution.',
            'The repeated process is Markov under that premise; the correct process object factorizes into one-step channels even though its terminal channel has a 16-dimensional minimal dilation.',
        ],
    }
    out = Path(__file__).with_name('endpoint-compression-result.json')
    out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
