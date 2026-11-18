#!/usr/bin/env python3
"""Exact native six-ovoid clock-bundle test.

The earlier four-dimensional obstruction required one symplectic form to be
fixed by both the native family clock and a selected native F20.  Here the
six *existing* antisymmetric two-qubit Pauli forms are retained as fibers.
The test uses the actual frozen v783 sigma, the actual compatible theta alpha,
and the exact native F20 matrices already intertwined with v419.
"""

from collections import Counter, deque
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
RESULT = HERE / "native-six-ovoid-clock-bundle-results.json"


def canonical(x):
    return sp.simplify(sp.expand_complex(x))


def is_zero(m):
    return all(canonical(x) == 0 for x in m)


def matrix_order(m, bound=120):
    ident = sp.eye(m.rows)
    p = ident
    for n in range(1, bound + 1):
        p = (p * m).applyfunc(canonical)
        if p == ident:
            return n
    raise RuntimeError("REQUIRE failed: matrix order exceeds bound")


def require(condition, message):
    """Runtime certificate check retained by python -O and -OO."""
    if not condition:
        raise RuntimeError("REQUIRE failed: " + message)


def compose(p, q):
    """Permutation p after q."""
    return tuple(p[q[i]] for i in range(len(p)))


def generated_group(generators):
    ident = tuple(range(len(generators[0])))
    seen = {ident}
    queue = deque([ident])
    while queue:
        x = queue.popleft()
        for g in generators:
            y = compose(g, x)
            if y not in seen:
                seen.add(y)
                queue.append(y)
    return seen


def permutation_order(p):
    ident = tuple(range(len(p)))
    x = ident
    for n in range(1, 61):
        x = compose(p, x)
        if x == ident:
            return n
    raise RuntimeError("REQUIRE failed: permutation order exceeds bound")


def conjugate_permutation(p, relabel):
    inv = [0] * len(relabel)
    for old, new in enumerate(relabel):
        inv[new] = old
    return tuple(relabel[p[inv[i]]] for i in range(len(p)))


def pgl2_five_group():
    """Exact PGL(2,5) action on P^1(F5), labeled 0..4,infinity=5."""
    perms = set()
    for a, b, c, d in itertools.product(range(5), repeat=4):
        if (a * d - b * c) % 5 == 0:
            continue
        image = []
        for x in range(6):
            if x == 5:
                image.append(5 if c == 0 else (a * pow(c, -1, 5)) % 5)
            else:
                den = (c * x + d) % 5
                image.append(5 if den == 0 else
                             ((a * x + b) * pow(den, -1, 5)) % 5)
        perms.add(tuple(image))
    require(len(perms) == 120, "PGL(2,5) projective action has order 120")
    return perms


def main():
    ii = sp.I
    sqrt2 = sp.sqrt(2)
    I2 = sp.eye(2)
    X = sp.Matrix(((0, 1), (1, 0)))
    Y = sp.Matrix(((0, -ii), (ii, 0)))
    Z = sp.diag(1, -1)
    one = {"I": I2, "X": X, "Y": Y, "Z": Z}
    pauli = {
        a + b: sp.kronecker_product(one[a], one[b])
        for a in "IXYZ" for b in "IXYZ"
    }

    # Exactly the six invertible alternating Pauli forms on the frozen C^4.
    forms = sorted(name for name, m in pauli.items() if m.T == -m)
    require(forms == ["IY", "XY", "YI", "YX", "YZ", "ZY"],
            "six alternating Pauli forms")

    # Actual native F20 selected in the earlier v783 subgroup census.
    C = sp.Matrix((
        (-ii/2, -ii/2, -ii/2, -ii/2),
        (-ii/2, ii/2, -ii/2, ii/2),
        (sp.Rational(1, 2), sp.Rational(1, 2), -sp.Rational(1, 2), -sp.Rational(1, 2)),
        (-sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), -sp.Rational(1, 2)),
    ))
    G = sp.Matrix((
        (0, 0, 0, ii),
        (0, 0, -1, 0),
        (ii, 0, 0, 0),
        (0, 1, 0, 0),
    ))
    require(matrix_order(C) == 5 and matrix_order(G) == 4,
            "native F20 generators have orders 5 and 4")
    require(is_zero(G * C * G.inv() - C**2),
            "native F20 Frobenius relation")

    # Original frozen v783 family clock and the reviewed theta alpha lift.
    sigma = sp.Matrix(((0, 0, 1, 0), (1, 0, 0, 0),
                       (0, 1, 0, 0), (0, 0, 0, 1)))
    U = sp.Matrix(((1, 0, 0, 0), (0, 0, 0, -1),
                   (0, 0, 1, 0), (0, 1, 0, 0)))
    zeta8 = (1 + ii) / sqrt2
    alpha = (zeta8 * sp.kronecker_product(X, I2) * U).applyfunc(canonical)
    J = ii * sp.eye(4)
    require(matrix_order(sigma) == 3, "native sigma has order 3")
    require(sp.trace(sigma) == 1,
            "native sigma has complex trace 1 on C^4")
    require(is_zero(alpha**2 - ii * sp.eye(4)), "theta alpha squares to iI")
    require(is_zero(alpha * sigma * alpha.conjugate().T - sigma**2),
            "theta alpha inverts native sigma")
    require(J * sigma == sigma * J and J * C == C * J and J * G == G * J,
            "central J commutes with sigma and F20")

    phase_set = (sp.Integer(1), -sp.Integer(1), ii, -ii)

    def form_action(operator):
        perm = []
        phases = []
        for name in forms:
            image = (operator * pauli[name] * operator.T).applyfunc(canonical)
            hits = [
                (phase, target)
                for target in forms for phase in phase_set
                if is_zero(image - phase * pauli[target])
            ]
            require(len(hits) == 1,
                    "operator has a unique projective image on each form")
            phase, target = hits[0]
            perm.append(forms.index(target))
            phases.append(str(phase))
        return tuple(perm), phases

    actions = {}
    for name, op in (("C", C), ("G", G), ("sigma", sigma),
                     ("alpha", alpha), ("J", J)):
        perm, phases = form_action(op)
        actions[name] = {
            "images": {forms[i]: forms[perm[i]] for i in range(6)},
            "phases": {forms[i]: phases[i] for i in range(6)},
            "permutation": list(perm),
        }

    pC = tuple(actions["C"]["permutation"])
    pG = tuple(actions["G"]["permutation"])
    pS = tuple(actions["sigma"]["permutation"])
    pA = tuple(actions["alpha"]["permutation"])
    pJ = tuple(actions["J"]["permutation"])

    group_f20 = generated_group([pC, pG])
    group_s5 = generated_group([pC, pG, pS])
    group_s6 = generated_group([pC, pG, pS, pA])
    require(len(group_f20) == 20, "C,G induce F20 of order 20")
    require(len(group_s5) == 120, "C,G,sigma group has order 120")
    require(len(group_s6) == 720, "C,G,sigma,alpha group has order 720")
    s5_order_census = {
        str(k): v for k, v in sorted(
            Counter(permutation_order(p) for p in group_s5).items()
        )
    }
    require(s5_order_census == {
        "1": 1, "2": 25, "3": 20, "4": 30, "5": 24, "6": 20
    }, "order census of the 120-element group")
    require(pJ == tuple(range(6)), "J acts projectively trivially on fibers")
    require(len({p[0] for p in group_s5}) == 6,
            "C,G,sigma group is transitive on six fibers")

    # Prove the named identification, rather than infer it from the census:
    # find an explicit relabeling under which the computed permutation group
    # is exactly PGL(2,5) on P^1(F5).  The classical faithful degree-six
    # action gives PGL(2,5) ~= S5.
    pgl2 = pgl2_five_group()
    pgl2_relabel = next((r for r in itertools.permutations(range(6))
                         if {conjugate_permutation(p, r) for p in group_s5}
                         == pgl2), None)
    require(pgl2_relabel is not None,
            "computed 120-group is conjugate to PGL(2,5) in S6")

    # Each form canonically determines the five Hodge/Clifford currents.
    # They are the six point-ovoids of the Pauli doily, not the raw spread
    # labels of v783; those are related by v783's point-line duality.
    ovoids = {}
    nontrivial = [name for name in sorted(pauli) if name != "II"]
    for omega_name in forms:
        omega = pauli[omega_name]
        currents = [
            name for name in nontrivial
            if is_zero(pauli[name].T * omega - omega * pauli[name])
        ]
        require(len(currents) == 5, "each alternating form selects five currents")
        require(all(is_zero(pauli[a] * pauli[b] + pauli[b] * pauli[a])
                    for i, a in enumerate(currents) for b in currents[i+1:]),
                "each selected five-current set pairwise anticommutes")
        ovoids[omega_name] = currents
    require(ovoids["IY"] == ["XI", "YX", "YY", "YZ", "ZI"],
            "IY fiber is the identified Hodge K5")
    require(len({tuple(v) for v in ovoids.values()}) == 6,
            "six current ovoids are distinct")
    require(all(sum(name in v for v in ovoids.values()) == 2
                for name in nontrivial),
            "every nontrivial Pauli lies in two ovoids")
    require(all(len(set(ovoids[a]) & set(ovoids[b])) == 1
                for i, a in enumerate(forms) for b in forms[i+1:]),
            "each pair of ovoids intersects in one Pauli")

    # Direct covariance of both the form and its five-current fiber.
    for op_name, op in (("C", C), ("G", G), ("sigma", sigma),
                        ("alpha", alpha), ("J", J)):
        image_map = actions[op_name]["images"]
        for omega_name in forms:
            target = image_map[omega_name]
            transported = set()
            for current in ovoids[omega_name]:
                image = (op * pauli[current] * op.conjugate().T).applyfunc(canonical)
                hit = next(name for name in nontrivial
                           if any(is_zero(image - phase * pauli[name])
                                  for phase in phase_set))
                transported.add(hit)
            require(transported == set(ovoids[target]),
                    "full five-current fibers transport covariantly")

    result = {
        "verdict": "NATIVE_SIX_FIBER_CLOCK_BUNDLE_EXACT",
        "fiber_labels": forms,
        "fiber_ovoids": ovoids,
        "actions": actions,
        "permutation_group_orders": {
            "F20_C_G": len(group_f20),
            "F20_plus_sigma": len(group_s5),
            "F20_plus_sigma_plus_alpha": len(group_s6),
        },
        "S5_element_order_census": s5_order_census,
        "S5_identification_certificate": {
            "model": "PGL(2,5) on P^1(F5), classically isomorphic to S5",
            "fiber_to_projective_line_relabeling": list(pgl2_relabel),
            "exact_group_equality_after_relabeling": True,
        },
        "trace_typing": {
            "native_sigma_complex_C4": 1,
            "native_sigma_realification_R8": 2,
        },
        "transitive_with_sigma": True,
        "named_group_identifications": {
            "F20_C_G": "F20 stabilizer of the marked IY fiber",
            "F20_plus_sigma": "transitive S5 action on six fibers",
            "F20_plus_sigma_plus_alpha": "S6 action on six fibers",
            "J": "projectively trivial on fibers",
        },
        "matrix_relations": {
            "C5": "I",
            "G4": "I",
            "G_C_Ginv": "C^2",
            "sigma3": "I",
            "alpha2": "iI",
            "alpha_sigma_alpha_inv": "sigma^-1",
        },
        "source_scope": [
            "C and G are the selected native F20 word pair already exactly intertwined with v419",
            "sigma is the frozen v783 native family clock",
            "alpha is the reviewed scalar-square theta lift",
            "the bundle removes the false requirement of one sigma-fixed symplectic form",
            "selection of this conjugate native F20 as the physical v419 subgroup remains a source marking premise",
            "no P1/P2 state, time, or source instrument follows from the finite bundle alone",
        ],
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    RESULT.write_text(payload)
    print(payload, end="")
    print("RESULT_SHA256", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
