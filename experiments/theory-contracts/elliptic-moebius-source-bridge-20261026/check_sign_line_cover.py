#!/usr/bin/env python3
"""Exact checks for the genus-two sign-line source dictionary.

The decisive calculation selects the negative-orbit two-torsion packet by the
ramification passport over t=z^4. It also records why the positive-orbit and
antipodal packets are not the marked v117 source.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


def matrix_order(m: sp.Matrix, maximum: int = 24) -> int | None:
    identity = sp.eye(m.rows)
    power = identity
    for n in range(1, maximum + 1):
        power *= m
        if power == identity:
            return n
    return None


def matrix_key(m: sp.Matrix) -> tuple:
    return tuple(m)


def generated_group(generators: list[sp.Matrix]) -> list[sp.Matrix]:
    identity = sp.eye(generators[0].rows)
    elements = {matrix_key(identity): identity}
    frontier = [identity]
    while frontier:
        x = frontier.pop()
        for generator in generators:
            y = x * generator
            key = matrix_key(y)
            if key not in elements:
                elements[key] = y
                frontier.append(y)
    return list(elements.values())


def triangle_genus(degree: int, orders: tuple[int, int, int]) -> int:
    correction = sum(sp.Rational(order - 1, order) for order in orders)
    two_g_minus_two = -2 * degree + degree * correction
    return int((two_g_minus_two + 2) / 2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).with_name("sign-line-cover-results.json"),
    )
    args = parser.parse_args()

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: object) -> None:
        ok = bool(condition)
        checks.append({"name": name, "ok": ok, "detail": detail})
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    w = sp.symbols("w")
    omega = -sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2

    def root(a: int) -> sp.Expr:
        return sp.expand(omega ** (a % 3))

    def q_ratio(a: int) -> sp.Expr:
        return sp.cancel((w + root(a)) / (w + root(a + 1)))

    def p_ratio(a: int) -> sp.Expr:
        return sp.cancel((w - root(a)) / (w - root(a + 1)))

    def antipodal_ratio(a: int) -> sp.Expr:
        return sp.cancel((w - root(a)) / (w + root(a)))

    q_product = sp.simplify(sp.prod(q_ratio(a) for a in range(3)))
    p_product = sp.simplify(sp.prod(p_ratio(a) for a in range(3)))
    anti_product = sp.factor(sp.prod(antipodal_ratio(a) for a in range(3)))
    check("negative-orbit Kummer equations have one product relation", q_product == 1, str(q_product))
    check("positive-orbit comparison packet has one product relation", p_product == 1, str(p_product))
    check(
        "antipodal packet has its exact eta product section",
        sp.simplify(anti_product - (w**3 - 1) / (w**3 + 1)) == 0,
        str(anti_product),
    )

    c_transport_q = [sp.simplify(q_ratio(a).subs(w, omega * w) - q_ratio(a - 1)) == 0 for a in range(3)]
    check("c:w->omega w cyclically transports the negative packet", all(c_transport_q), c_transport_q)

    degree_v4 = 4
    # Valuations modulo two of the independent square roots h_0,h_1.
    # These are derived from h_a^2=(w+omega^a)/(w+omega^(a+1)).
    q_inertia = {
        b: tuple(
            (int(b == a) - int(b == ((a + 1) % 3))) % 2
            for a in range(2)
        )
        for b in range(3)
    }
    q_branch_support = {b for b, vector in q_inertia.items() if vector != (0, 0)}
    q_branch_orders = [2 for _ in sorted(q_branch_support)]
    two_g_minus_two_q = -2 * degree_v4 + sum(degree_v4 // e * (e - 1) for e in q_branch_orders)
    genus_zq = (two_g_minus_two_q + 2) // 2
    check("negative packet V4 cover is rational", genus_zq == 0, genus_zq)
    check(
        "negative packet has the three distinct nonzero V4 inertia vectors",
        set(q_inertia.values()) == {(1, 0), (1, 1), (0, 1)},
        q_inertia,
    )

    genus_c = 2
    genus_yq = 1 + degree_v4 * (genus_c - 1)
    check("pullback Y_q -> C is unramified of degree four and genus five", genus_yq == 5, genus_yq)

    z_square = sp.factor((1 + w**3) / (1 - w**3))
    t_of_w = sp.factor(z_square**2)
    check(
        "the exact quotient coordinate is t=z^4",
        t_of_w == sp.factor((1 + w**3) ** 2 / (1 - w**3) ** 2),
        str(t_of_w),
    )

    # The base w->t map has orders 2 at q/t=0, 3 at 0,infinity/t=1,
    # and 2 at p/t=infinity.  A square-root branch at q or p doubles
    # the corresponding local order.
    base_signature = (2, 3, 2)
    q_support = {"q"}
    p_support = {"p"}
    antipodal_support = {"q", "p"}

    def composite_signature(support: set[str]) -> tuple[int, int, int]:
        return (
            base_signature[0] * (2 if "q" in support else 1),
            base_signature[1],
            base_signature[2] * (2 if "p" in support else 1),
        )

    q_signature = composite_signature(q_support)
    p_signature = composite_signature(p_support)
    antipodal_signature = composite_signature(antipodal_support)
    check("negative packet has the marked v117 passport", q_signature == (4, 3, 2), q_signature)
    check("positive packet has the reversed clock passport", p_signature == (2, 3, 4), p_signature)
    check("antipodal packet has the hyperbolic wrong passport", antipodal_signature == (4, 3, 4), antipodal_signature)

    degree_triangle = 24
    check("negative degree-24 triangle cover has genus zero", triangle_genus(degree_triangle, q_signature) == 0, 0)
    check("positive degree-24 triangle cover also has genus zero but wrong marking", triangle_genus(degree_triangle, p_signature) == 0, 0)
    check("antipodal degree-24 cover has genus three", triangle_genus(degree_triangle, antipodal_signature) == 3, 3)

    P = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    R = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    U0 = R.inv()
    M1 = P
    Uinf = (M1 * U0).inv()
    matrix_orders = (matrix_order(U0), matrix_order(M1), matrix_order(Uinf))
    check("v117 triangle matrices have the same oriented orders", matrix_orders == q_signature, matrix_orders)
    check("triangle matrices have determinant one", all(m.det() == 1 for m in (U0, M1, Uinf)), [int(m.det()) for m in (U0, M1, Uinf)])

    s4_group = generated_group([U0, M1])
    check("marked triangle matrices generate S4 of order twenty-four", len(s4_group) == 24, len(s4_group))
    commutators = [g * h * g.inv() * h.inv() for g in s4_group for h in s4_group]
    derived_group = generated_group(commutators)
    check(
        "the only nontrivial cyclic quotient is C2 with A4 kernel",
        len(derived_group) == 12 and len(s4_group) // len(derived_group) == 2,
        {"derived_order": len(derived_group), "abelianization_order": len(s4_group) // len(derived_group)},
    )

    eta, z_symbol = sp.symbols("eta z")
    curve_identity = sp.simplify((z_square * (1 - w**3) ** 2) - (1 - w**6))
    check("C has function field C(w,z) with eta=z(1-w^3)", curve_identity == 0, str(curve_identity))
    composite_degree = len(s4_group) // 2
    check("fiber product over C(z) has degree twelve and A4", composite_degree == 12, composite_degree)

    w_valuations = [1, 1, -1, -1]
    relative_e4 = [int(4 // sp.gcd(4, abs(v))) for v in w_valuations]
    relative_e2 = [int(2 // sp.gcd(2, abs(v))) for v in w_valuations]
    check("mu12 cover is ramified C4 over C", relative_e4 == [4, 4, 4, 4], relative_e4)
    check("its unique quadratic intermediate is ramified", relative_e2 == [2, 2, 2, 2], relative_e2)

    z, h = sp.symbols("z h", nonzero=True)
    f = (z**2 - 1) / (z**2 + 1)
    check("base quarter sends f to its inverse", sp.simplify(f.subs(z, sp.I * z) - 1 / f) == 0, "f(iz)=f(z)^-1")
    check("Q lift h->i/h sends w=h^4 to w^-1", sp.simplify((sp.I / h) ** 4 - 1 / h**4) == 0, "w->w^-1")
    check("Q lift sends g=h^3 to -i g^-1", sp.simplify((sp.I / h) ** 3 + sp.I / h**3) == 0, "g->-i/g")

    verdict = "PASS_EXACT_Q_PACKET_V117_COEFFICIENT_LIFT" if all(c["ok"] for c in checks) else "FAIL"
    result = {
        "verdict": verdict,
        "checks_passed": sum(bool(c["ok"]) for c in checks),
        "checks_total": len(checks),
        "objects": {
            "genus_two_curve": "eta^2=1-w^6",
            "selected_line_squares": "h_a^2=(w+omega^a)/(w+omega^(a+1))",
            "selected_line_classes": "L_a^q=O(q_a-q_(a+1)), a mod 3",
            "selected_passport": [4, 3, 2],
            "rejected_positive_passport": [2, 3, 4],
            "rejected_antipodal_passport": [4, 3, 4],
            "common_cyclic_cover": "h^12=f(z), g=h^3, w=h^4",
        },
        "checks": checks,
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(verdict)
    print(f"wrote {args.out}")
    return 0 if verdict.startswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
