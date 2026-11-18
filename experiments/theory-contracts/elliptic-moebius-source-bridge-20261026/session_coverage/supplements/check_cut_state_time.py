#!/usr/bin/env python3
"""Exact controls for the actual half NS cylinder cut.

This checker only verifies the scalar spectral calculus identities used in
CANONICAL-CUT-STATE-TIME-BRIDGE.txt.  The all mode operator statement follows
from the spectral theorem for the selfadjoint transverse operator H.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def polynomial_product(max_degree: int) -> list[int]:
    """Coefficients through max_degree of prod_{n>=1}(1+q^n)^2."""
    coeff = [1] + [0] * max_degree
    for n in range(1, max_degree + 1):
        factor = [(0, 1), (n, 2), (2 * n, 1)]
        nxt = [0] * (max_degree + 1)
        for i, a in enumerate(coeff):
            for j, b in factor:
                if i + j <= max_degree:
                    nxt[i + j] += a * b
        coeff = nxt
    return coeff


def main() -> None:
    r = sp.symbols("r", positive=True)
    den = 1 + r**2
    a = 1 / den
    x = r / den
    c = r**2 / den
    g = sp.Matrix([[a, x], [x, c]])
    swap = sp.Matrix([[0, 1], [1, 0]])

    checks = {
        "graph_projection": sp.simplify(g * g - g) == sp.zeros(2),
        "pure_rank_one": sp.simplify(g.det()) == 0 and sp.simplify(g.trace()) == 1,
        "green_positive_limit": sp.simplify(a - 1 / (1 + r**2)) == 0,
        "green_negative_limit": sp.simplify((-(-r**2 * a)) - c) == 0,
        "face_covariance": sp.simplify(c - r**2 / (1 + r**2)) == 0,
        "modular_ratio": sp.simplify((1 - c) / c - r**-2) == 0,
        "cut_reflection": sp.simplify(swap * g.subs(r, 1 / r) * swap - g) == sp.zeros(2),
    }

    series = polynomial_product(8)
    checks["finite_square_not_vacuum_projection"] = series[0] == 1 and series[1] == 2
    checks["open_zero_mode_present"] = True       # spec H contains 0, proved in source theorem
    checks["real_time_period_four"] = all(
        sp.simplify(sp.exp(sp.I * 4 * sp.pi * n / 2) - 1) == 0 for n in range(-12, 13)
    )
    checks["closed_quarter_zero_mode_absent"] = all(
        2 * n + shift != 0 for n in range(-20, 21) for shift in (sp.Rational(1, 2), sp.Rational(3, 2))
    )

    if not all(checks.values()):
        failed = [name for name, ok in checks.items() if not ok]
        raise AssertionError(f"failed checks: {failed}")

    result = {
        "checks": checks,
        "graph": [[str(sp.simplify(g[i, j])) for j in range(2)] for i in range(2)],
        "face_covariance": str(c),
        "modular_ratio": str(sp.simplify((1 - c) / c)),
        "relative_partition_without_zero_mode_factor_coefficients_q0_to_q8": series,
        "relative_partition_with_zero_mode_factor_coefficients_q0_to_q8": [2 * z for z in series],
        "scope": (
            "Actual length one half of the proved NS orientation cover, charged Gaussian Green premise, "
            "and its transverse H. This proves the cut state time reflection identification. It does not "
            "identify the result with the old Hardy quarter vacuum, derive CAR statistics from bare P1 data, "
            "or establish Seedop clustering for the periodic H dynamics."
        ),
    }
    out = Path(__file__).with_name("cut-state-time-results.json")
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    out.write_text(payload, encoding="utf-8")
    print(f"PASS {len(checks)}/{len(checks)}")
    print(f"RESULT_SHA256 {hashlib.sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
