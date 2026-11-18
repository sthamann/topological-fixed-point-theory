#!/usr/bin/env python3
"""Both D8 Ramond modules and same-family regular OPE coefficients.

Uses the existing E8 lattice/cocycle/response reader. No module is selected as
physical matter, and no coefficient of B[Phi] is supplied. A representative shift
for the second module is a basis convention, not a coupling or sector choice.
"""
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path
import sys

import e8_boundary_response as e8

sp = e8.sp


def chirality(p):
    return sum(x < 0 for x in p) % 2


def both_ramond_weights():
    return tuple(tuple(sp.Rational(x, 2) for x in signs)
                 for signs in product((-1, 1), repeat=8))


def module_phase_argument(p):
    # S weights are in E8. For C, adding a representative of the vector coset
    # puts the phase argument in E8. Other representatives give equivalent bases.
    return p if chirality(p) == 0 else e8.add(p, (0,)*7+(1,))


def ground_action(nu, weight, weights, index):
    action = {}
    for k, p in enumerate(weights):
        if e8.dot(nu, p) != -weight:
            continue
        target = e8.add(nu, p)
        assert target in index
        assert chirality(target) == chirality(p)
        action[k] = (index[target], e8.epsilon(nu, module_phase_argument(p)))
    return action


def root_module_check(weights, index):
    roots = tuple(r for r in e8.roots() if all(x.is_Integer for x in r))
    actions = {r: ground_action(r, 1, weights, index) for r in roots}
    for r in roots:
        for s in roots:
            nu = e8.add(r, s)
            expected = actions.get(nu)
            factor = e8.epsilon(r, s) if expected is not None else None
            ar, ass = actions[r], actions[s]
            for k, p in enumerate(weights):
                terms = defaultdict(int)
                if k in ass:
                    middle, x = ass[k]
                    if middle in ar:
                        target, y = ar[middle]
                        terms[target] += x*y
                if k in ar:
                    middle, x = ar[k]
                    if middle in ass:
                        target, y = ass[middle]
                        terms[target] -= x*y
                if expected is not None and k in expected:
                    target, x = expected[k]
                    terms[target] -= factor*x
                elif s == e8.neg(r):
                    terms[k] -= e8.epsilon(r, s)*e8.dot(r, p)
                assert not any(terms.values()), (r, s, p, terms)
        adj = e8.adjoint_sign(r)
        for k, (target, coefficient) in actions[r].items():
            assert actions[e8.neg(r)][target] == (k, adj*coefficient)
    return {"roots": len(roots), "ramond_module_dimensions": [128, 128],
            "root_pair_commutators": len(roots)**2,
            "scope": "Both stipulated D8 spin modules; not physical endpoint selection"}


def composite_check(weights, index):
    spin_s = tuple(p for p in weights if chirality(p) == 0)
    factorizations = defaultdict(list)
    for p in spin_s:
        for q in spin_s:
            if p[5:] != q[5:]:
                continue
            nu = e8.add(p, q)
            if sum(x != 0 for x in nu[:5]) != 1 or sum(x*x for x in nu[:5]) != 1:
                continue
            assert e8.dot(p, q) == 0 and e8.dot(nu, nu) == 4
            e8.coordinates(nu)
            factorizations[nu].append((p, q))
    assert len(factorizations) == 80
    assert {len(pairs) for pairs in factorizations.values()} == {8}
    per_state = Counter()
    response_values = Counter()
    parity_edges = Counter()
    sample = None
    for nu, pairs in factorizations.items():
        assert e8.grade(nu) == 2
        # A regular coefficient is already inside the existing vertex algebra.
        p, q = pairs[0]
        assert e8.epsilon(p, q) == e8.epsilon(q, p)
        action = ground_action(nu, 2, weights, index)
        assert len(action) == 16
        assert Counter(chirality(weights[k]) for k in action) == {0: 8, 1: 8}
        reverse = ground_action(e8.neg(nu), 2, weights, index)
        adj = e8.adjoint_sign(nu)
        for k, (target, coefficient) in action.items():
            incoming, outgoing = weights[k], weights[target]
            assert outgoing[5:] == e8.neg(incoming)[5:]
            assert reverse[target] == (k, adj*coefficient)
            assert e8.grade(outgoing) == (e8.grade(incoming)+2) % 4
            per_state[k] += 1
            parity_edges[chirality(incoming)] += 1
            if chirality(incoming) == 0:
                value = e8.current_response((nu, incoming, e8.neg(outgoing)), (3, 2, 1))
                norm = e8.current_response((outgoing, e8.neg(outgoing)), (2, 1))
                assert 4*value/norm == coefficient
                response_values[str(value)] += 1
                if sample is None:
                    sample = {"composite_charge": list(map(str, nu)),
                              "incoming": list(map(str, incoming)),
                              "outgoing": list(map(str, outgoing)),
                              "three_point_3_2_1": str(value),
                              "weight_two_zero_mode_coefficient": coefficient}
    assert len(per_state) == 256 and set(per_state.values()) == {5}
    assert dict(parity_edges) == {0: 640, 1: 640}
    return {"diagonal_family_composite_momenta": len(factorizations),
            "boundary_conformal_weight": 2,
            "ordered_spin_field_factorizations": sum(map(len, factorizations.values())),
            "factorizations_per_momentum": 8,
            "ground_state_actions_per_ramond_module": dict(parity_edges),
            "actions_per_spinor_weight": 5, "rank_per_composite_on_both_modules": 16,
            "S_module_existing_reader_values": dict(response_values), "sample": sample,
            "scope": "Bosonic OPE coefficients and Ramond-state actions, not physical Higgs/fermion identification"}


def discriminant_check():
    zero = (sp.S.Zero,)*8
    v = (sp.S.One,)+(sp.S.Zero,)*7
    s = (sp.Rational(1, 2),)*8
    c = e8.add(s, e8.neg(v))
    reps = [zero, v, s, c]
    monodromy = [[int((-1)**int(2*e8.dot(a, b))) for b in reps] for a in reps]
    subgroups = [[0], [0, 1], [0, 2], [0, 3], [0, 1, 2, 3]]
    local = [h for h in subgroups if all(monodromy[a][b] == 1 for a in h for b in h)]
    assert local == subgroups[:-1]
    return {"order": ["0", "v", "s", "c"],
            "ground_conformal_weights": [str(e8.dot(x, x)/2) for x in reps],
            "monodromy": monodromy,
            "local_extension_supports": [[['0','v','s','c'][i] for i in h] for h in local],
            "spinorial_cosets_have_integral_conformal_spin": True,
            "scope": "Fixed D8 chiral lattice extensions only; Ramond states are not declared local NS fields"}


def main():
    weights = both_ramond_weights()
    index = {p: i for i, p in enumerate(weights)}
    result = {"verdict": "RAMOND_ACTION_AND_WEIGHT_TWO_COMPOSITE_DERIVED_WITHIN_GIVEN_LATTICE",
              "discriminant": discriminant_check(),
              "ramond_action": root_module_check(weights, index),
              "composites": composite_check(weights, index),
              "physical_sector_selected": False, "physical_Phi_identified": False,
              "P1_field_dependent_boundary_operator_derived": False,
              "physical_gates_closed": []}
    if len(sys.argv)>1:
        Path(sys.argv[1]).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
