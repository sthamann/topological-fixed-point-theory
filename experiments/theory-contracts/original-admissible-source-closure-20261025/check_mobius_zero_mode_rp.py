#!/usr/bin/env python3
"""Exact zero-transverse-mode reflection-positivity test.

This is not a model construction.  It restricts the stated flat Moebius
Dirac Green operator to the invariant transverse zero branch from Theorem
8.2 of arXiv:2610.03745v1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def encode(x: sp.Expr) -> str:
    return str(sp.simplify(x))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    for label, q in (("plus_i", sp.I), ("minus_i", -sp.I)):
        # Green kernel of partial_x on q-twisted functions.  On (-1,0) and
        # (0,1) it is constant.  The jump is one and G(s+1)=q G(s).
        g_minus = sp.simplify(1 / (q - 1))
        g_plus = sp.simplify(q / (q - 1))
        occupation = sp.simplify(-g_minus)

        assert sp.simplify(g_plus - g_minus) == 1
        assert sp.simplify(g_plus - q * g_minus) == 0
        assert sp.im(g_minus) != 0
        assert sp.im(occupation) != 0

        # The same functional is the normalized one-mode trace with the
        # temporal insertion z=-q^{-1}.
        z = sp.simplify(-1 / q)
        c = sp.Matrix([[0, 1], [0, 0]])
        cdag = c.T
        number = cdag * c
        rho_num = sp.diag(1, z)
        omega_number = sp.simplify(sp.trace(rho_num * number) /
                                   sp.trace(rho_num))
        assert number * number == number
        assert sp.simplify(omega_number - occupation) == 0

        # Every unit-modulus phase in the standard charge-exchanging,
        # order-reversing fermionic reflection cancels on c^* c.  Thus its
        # reflected square has the same non-real value omega(N^2)=omega(N).
        reflected_number_norm = omega_number

        rows.append({
            "sector": label,
            "q": encode(q),
            "green_0_minus": encode(g_minus),
            "green_0_plus": encode(g_plus),
            "twisted_trace_insertion_z": encode(z),
            "occupation": encode(occupation),
            "occupation_real": bool(sp.im(occupation) == 0),
            "neutral_number_reflected_norm": encode(reflected_number_norm),
            "neutral_number_norm_real": bool(sp.im(reflected_number_norm) == 0),
            "standard_linear_reflection_target_q": encode(1 / q),
            "fixed_sector_under_standard_linear_reflection": bool(
                sp.simplify(1 / q - q) == 0),
        })

    ratio = sp.simplify(
        sp.sympify(rows[0]["green_0_minus"]) /
        sp.sympify(rows[1]["green_0_minus"])
    )
    assert ratio == sp.I

    result = {
        "claim": "fixed Pin-minus zero branch fails full CAR OS positivity",
        "convention": (
            "Euclidean skew Dirac is slashD=iD; on f0 it is b*partial_x. "
            "The Clifford b cancels in the reflected scalar kernel."
        ),
        "half_interval": "0 < x < 1/2",
        "sectors": rows,
        "ratio_of_two_standard_reflected_constants": encode(ratio),
        "conclusion": (
            "The standard linear Clifford reflection exchanges q and q^-1. "
            "A constant field phase may rotate a charged one-particle rank-one "
            "form, but cancels on the neutral number word; its reflected norm "
            "is non-real in both quarter sectors."
        ),
    }
    body = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.write_text(body, encoding="utf-8")
    print(hashlib.sha256(body.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
