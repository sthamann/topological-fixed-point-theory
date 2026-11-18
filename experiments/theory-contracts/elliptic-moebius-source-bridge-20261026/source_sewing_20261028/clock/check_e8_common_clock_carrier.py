#!/usr/bin/env python3
"""Exact E8 common-clock carrier test.

This checker asks one bounded question.  Does the *existing* E8 lattice carry
both the Coxeter/Galois clock of v223/v315/v419 and the native compiler clock
of v634, and what exact relations hold?  It does not construct a tensor
product or identify the two order-three hands.

All calculations use exact SymPy arithmetic.  Completeness of the G-lift
census follows from the unimodular cyclic Krylov basis K=(r,W r,...,W^7 r):
an integral isometry satisfying G W G^-1=W^7 is uniquely determined by G r,
and G r must be one of the 240 E8 roots.
"""

from collections import Counter, deque
import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
RESULT = HERE / "e8-common-clock-results.json"


def matrix_key(m):
    return tuple(int(x) for x in m)


def matrix_order(m, bound=240):
    ident = sp.eye(m.rows)
    p = ident
    for n in range(1, bound + 1):
        p = p * m
        if p == ident:
            return n
    raise RuntimeError("REQUIRE failed: matrix order exceeds bound")


def require(condition, message):
    """Runtime certificate check retained by python -O and -OO."""
    if not condition:
        raise RuntimeError("REQUIRE failed: " + message)


def charpoly(m):
    x = sp.symbols("x")
    return str(sp.factor(m.charpoly(x).as_expr()))


def main():
    I8 = sp.eye(8)
    e = sp.eye(8)
    h = sp.Rational(1, 2)

    # Bourbaki E8 simple roots, verbatim convention used in v781/v1.
    a1 = sp.Matrix([
        h * (e[j, 0] - e[j, 1] - e[j, 2] - e[j, 3]
             - e[j, 4] - e[j, 5] - e[j, 6] + e[j, 7])
        for j in range(8)
    ])
    simple = [
        a1, e[:, 0] + e[:, 1], e[:, 1] - e[:, 0],
        e[:, 2] - e[:, 1], e[:, 3] - e[:, 2],
        e[:, 4] - e[:, 3], e[:, 5] - e[:, 4],
        e[:, 6] - e[:, 5],
    ]
    B = sp.Matrix.hstack(*simple)
    gram = B.T * B
    require(B.det() == -1 and gram.det() == 1,
            "Bourbaki basis and E8 Gram determinants")

    def reflection_std(a):
        return I8 - a * a.T

    refl_std = [reflection_std(a) for a in simple]
    W_std = I8
    for r in refl_std:
        W_std = W_std * r
    W = sp.simplify(B.inv() * W_std * B)
    require(matrix_order(W) == 30 and W**15 == -I8,
            "Coxeter element has order 30 and W^15=-I")

    # Native v634 operators in the same real E8 coordinates.
    J_std = sp.zeros(8)
    for k in range(0, 8, 2):
        J_std[k, k + 1] = -1
        J_std[k + 1, k] = 1
    sigma_std = sp.zeros(8)
    for row, col in enumerate([4, 5, 0, 1, 2, 3, 6, 7]):
        sigma_std[row, col] = 1
    J = sp.simplify(B.inv() * J_std * B)
    sigma = sp.simplify(B.inv() * sigma_std * B)
    compiler = J * sigma
    require(all(J_std * e[:, k] == e[:, k + 1] and
                J_std * e[:, k + 1] == -e[:, k]
                for k in range(0, 8, 2)),
            "v634 J is multiplication by +i in the x+iy convention")
    require(J.T * gram * J == gram, "native J preserves the E8 Gram form")
    require(sigma.T * gram * sigma == gram,
            "native sigma preserves the E8 Gram form")
    require(J * sigma == sigma * J, "native J and sigma commute")
    require(matrix_order(J) == 4 and matrix_order(sigma) == 3,
            "native J and sigma have orders 4 and 3")
    require(matrix_order(compiler) == 12, "native compiler has order 12")
    require(compiler**4 == sigma and compiler**9 == J,
            "compiler powers recover sigma and J")

    # Generate all 240 E8 roots in the simple-root coordinates.
    reflections = [sp.simplify(B.inv() * r * B) for r in refl_std]
    roots = set()
    queue = deque()
    for i in range(8):
        for sign in (1, -1):
            v = tuple(sign * int(x) for x in e[:, i])
            roots.add(v)
            queue.append(v)
    while queue:
        v = sp.Matrix(queue.popleft())
        for r in reflections:
            z = tuple(int(x) for x in r * v)
            if z not in roots:
                roots.add(z)
                queue.append(z)
    require(len(roots) == 240, "E8 root orbit has 240 roots")

    # Complete census of integral lifts of exponent multiplication m -> 7m.
    root0 = e[:, 0]
    krylov = sp.Matrix.hstack(*[W**k * root0 for k in range(8)])
    require(abs(krylov.det()) == 1, "Coxeter Krylov basis is unimodular")
    lifts = []
    for image in sorted(roots):
        r = sp.Matrix(image)
        candidate = sp.simplify(
            sp.Matrix.hstack(*[W**(7 * k) * r for k in range(8)])
            * krylov.inv()
        )
        if candidate.T * gram * candidate == gram:
            require(candidate * W * candidate.inv() == W**7,
                    "enumerated isometry conjugates W to W^7")
            lifts.append(candidate)
    require(len(lifts) == 30, "complete Galois-lift census has size 30")

    fourth_power = {}
    order_counts = Counter()
    for g in lifts:
        o = matrix_order(g)
        order_counts[o] += 1
        fourth_power.setdefault(str(o), Counter())
        hit = next(k for k in range(30) if g**4 == W**k)
        fourth_power[str(o)][str(hit)] += 1
    require(order_counts == Counter({8: 10, 24: 20}),
            "Galois-lift order census")
    require(fourth_power == {
        "8": Counter({"15": 10}),
        "24": Counter({"5": 10, "25": 10}),
    }, "Galois-lift fourth-power census")

    order8 = sorted((g for g in lifts if matrix_order(g) == 8), key=matrix_key)
    G = order8[0]  # certificate representative only; no source-canonical choice.
    require(G**4 == W**15 and W**15 == -I8,
            "order-eight lift has fourth power -I")
    require(G * W**6 * G.inv() == W**12,
            "order-eight lift acts by Frobenius on C5")
    require(G * W**10 * G.inv() == W**10,
            "order-eight lift centralizes Coxeter C3")
    require(G**2 not in (I8, -I8),
            "order-eight lift square is not scalar")

    # The lift does not identify or normalize the native compiler hands.
    native_normalizer_counts = {
        "J_to_J": 0,
        "J_to_J_inverse": 0,
        "sigma_to_sigma": 0,
        "sigma_to_sigma_inverse": 0,
        "sigma_to_coxeter_family": 0,
        "sigma_to_coxeter_family_inverse": 0,
    }
    for g in lifts:
        gj = g * J * g.inv()
        gs = g * sigma * g.inv()
        native_normalizer_counts["J_to_J"] += int(gj == J)
        native_normalizer_counts["J_to_J_inverse"] += int(gj == J.inv())
        native_normalizer_counts["sigma_to_sigma"] += int(gs == sigma)
        native_normalizer_counts["sigma_to_sigma_inverse"] += int(gs == sigma.inv())
        native_normalizer_counts["sigma_to_coxeter_family"] += int(gs == W**10)
        native_normalizer_counts["sigma_to_coxeter_family_inverse"] += int(gs == W**20)
    require(all(v == 0 for v in native_normalizer_counts.values()),
            "none of the 30 lifts identifies the native compiler hands")

    def orbit_census(generators):
        unseen = set(roots)
        census = Counter()
        while unseen:
            start = next(iter(unseen))
            unseen.remove(start)
            orb = {start}
            q = deque([start])
            while q:
                x = sp.Matrix(q.popleft())
                for gen in generators:
                    y = tuple(int(v) for v in gen * x)
                    if y not in orb:
                        orb.add(y)
                        unseen.remove(y)
                        q.append(y)
            census[len(orb)] += 1
        return dict(sorted(census.items()))

    spectra = {}
    for name, matrix in [
        ("native_J", J),
        ("native_sigma", sigma),
        ("native_compiler", compiler),
        ("coxeter_W", W),
        ("coxeter_carrier_C5", W**6),
        ("coxeter_family_C3", W**10),
        ("galois_lift_G8", G),
        ("galois_lift_G8_square", G**2),
    ]:
        spectra[name] = {
            "order": matrix_order(matrix),
            "trace": int(sp.trace(matrix)),
            "charpoly": charpoly(matrix),
        }

    # This character mismatch is the decisive typed obstruction.
    require(spectra["native_sigma"] != spectra["coxeter_family_C3"],
            "native and Coxeter C3 spectra differ")
    require(sp.trace(sigma) == 2 and sp.trace(W**10) == -4,
            "native and Coxeter C3 real traces are 2 and -4")

    result = {
        "verdict": "COMMON_E8_CARRIER_WITH_DISTINCT_CLOCK_SUBALGEBRAS",
        "exact_arithmetic": True,
        "root_count": len(roots),
        "gram_det": int(gram.det()),
        "krylov_det": int(krylov.det()),
        "galois_lift_census": {
            "total": len(lifts),
            "orders": {str(k): v for k, v in sorted(order_counts.items())},
            "fourth_power_W_exponents": {
                o: dict(sorted(c.items())) for o, c in fourth_power.items()
            },
        },
        "canonical_selection": False,
        "representative_G8_matrix": [
            [int(G[r, c]) for c in range(8)] for r in range(8)
        ],
        "relations": {
            "G_W_Ginv": "W^7",
            "G_C5_Ginv": "C5^2",
            "G_C3_Ginv": "C3",
            "G4": "W^15=-I",
            "compiler": "c=J*sigma, c^4=sigma, c^9=J",
            "native_J_coordinate_sign": "J(x,y)=(-y,x), multiplication by +i",
        },
        "native_normalizer_counts_among_30_lifts": native_normalizer_counts,
        "spectra": spectra,
        "root_orbits": {
            "compiler_and_coxeter": orbit_census([compiler, W]),
            "all_named_operators": orbit_census([J, sigma, W, G]),
        },
        "scope": [
            "common integral E8 lattice and 240-root carrier",
            "no equality of native and cyclotomic order-three clocks",
            "no source-derived P1/P2 sewing functor",
            "no identification of the scalar-square theta alpha with G8",
        ],
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    RESULT.write_text(payload)
    print(payload, end="")
    print("RESULT_SHA256", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
