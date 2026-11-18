#!/usr/bin/env python3
"""Exact conformal-symplectic structure behind the v419 F20 bridge."""

from __future__ import annotations

from itertools import combinations
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
PIN = "cc6152ad7ea68c67725020e386e49eec32810ab25bb483c76bae99964399c2fd"


def require(value: bool, label: str, checks: list[str]) -> None:
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)


def permutation(images: list[int]) -> sp.Matrix:
    result = sp.zeros(len(images))
    for source, target in enumerate(images):
        result[target, source] = 1
    return result


def exterior_square(matrix: sp.Matrix) -> sp.Matrix:
    pairs = list(combinations(range(matrix.rows), 2))
    return sp.Matrix.hstack(*[
        sp.Matrix([
            matrix[i, a] * matrix[j, b] - matrix[j, a] * matrix[i, b]
            for i, j in pairs
        ])
        for a, b in pairs
    ])


def sym2_character(matrix: sp.Matrix) -> sp.Expr:
    return sp.simplify((sp.trace(matrix) ** 2 + sp.trace(matrix**2)) / 2)


def ext2_character(matrix: sp.Matrix) -> sp.Expr:
    return sp.simplify((sp.trace(matrix) ** 2 - sp.trace(matrix**2)) / 2)


def encode(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix.tolist()]


def main() -> None:
    checks: list[str] = []
    source = ROOT / "verification/v419_seam_galois_carrier.py"
    require(hashlib.sha256(source.read_bytes()).hexdigest() == PIN,
            "v419 source pin", checks)

    x = sp.symbols("x")
    phi5 = x**4 + x**3 + x**2 + x + 1
    c5 = sp.Matrix([
        [0, 0, 0, -1], [1, 0, 0, -1],
        [0, 1, 0, -1], [0, 0, 1, -1],
    ])

    def column(exponent: int) -> sp.Matrix:
        poly = sp.Poly(sp.rem(x**exponent, phi5, x), x)
        return sp.Matrix([poly.coeff_monomial(x**j) for j in range(4)])

    galois = sp.Matrix.hstack(*(column(2 * k) for k in range(4)))

    # The complete two-dimensional C5-invariant alternating pencil.
    omega_a = sp.Matrix([
        [0,0,1,-1], [0,0,0,1], [-1,0,0,0], [1,-1,0,0],
    ])
    omega_b = sp.Matrix([
        [0,1,0,0], [-1,0,1,0], [0,-1,0,1], [0,0,-1,0],
    ])
    omega_plus = sp.I * omega_a + omega_b
    omega_minus = -sp.I * omega_a + omega_b

    require(c5.T * omega_a * c5 == omega_a and
            c5.T * omega_b * c5 == omega_b,
            "C5 invariant alternating pencil", checks)
    require(galois.T * omega_plus * galois == sp.I * omega_plus,
            "Omega plus multiplier i", checks)
    require(galois.T * omega_minus * galois == -sp.I * omega_minus,
            "Omega minus multiplier minus i", checks)
    require(omega_plus.rank() == omega_minus.rank() == 4,
            "both alternating eigenforms rank four", checks)
    require(omega_plus.det() == 3 - 4*sp.I and
            omega_minus.det() == 3 + 4*sp.I,
            "nonzero determinants 3 minusplus 4i", checks)

    zeta8 = (1 + sp.I) / sp.sqrt(2)
    lifted_g = zeta8**-1 * galois
    require(sp.simplify(lifted_g.T * omega_plus * lifted_g - omega_plus) == sp.zeros(4),
            "mu8 normalized G is symplectic", checks)
    require(sp.simplify(lifted_g**4 + sp.eye(4)) == sp.zeros(4) and
            sp.simplify(lifted_g**8 - sp.eye(4)) == sp.zeros(4),
            "mu8 normalized G has order eight", checks)
    require(sp.simplify(lifted_g * c5 * lifted_g.inv() - c5**2) == sp.zeros(4),
            "central normalization preserves Frobenius relation", checks)

    pc = permutation([(j + 1) % 5 for j in range(5)])
    pg = permutation([(2 * j) % 5 for j in range(5)])

    vector_rows = []
    adjoint_rows = []
    for b in range(4):
        for a in range(5):
            raw = c5**a * galois**b
            multiplier = sp.I**b
            # The inverse symplectic bivector is the one-dimensional line of
            # Lambda2(V), with the same conformal multiplier.
            primitive_normalized = sp.simplify(
                multiplier**-1 * (ext2_character(raw) - multiplier)
            )
            sign_permutation = sp.simplify((-1)**b * sp.trace(pc**a * pg**b))
            require(primitive_normalized == sign_permutation,
                    f"Spin5 vector character a={a} b={b}", checks)
            vector_rows.append({
                "a": a, "b": b,
                "primitive_normalized": str(primitive_normalized),
                "sign_permutation": str(sign_permutation),
            })

            normalized_spin_adjoint = sp.simplify(multiplier**-1 * sym2_character(raw))
            # Lambda2(sign tensor Perm5)=Lambda2(Perm5).
            so5_adjoint = ext2_character(pc**a * pg**b)
            require(normalized_spin_adjoint == so5_adjoint,
                    f"Spin5 adjoint character a={a} b={b}", checks)
            adjoint_rows.append({
                "a": a, "b": b,
                "sym2_normalized": str(normalized_spin_adjoint),
                "lambda2_vector5": str(so5_adjoint),
            })

    output = {
        "verdict": "EXACT_CONFORMAL_SP4_AND_SPIN5_STRUCTURE",
        "invariant_pencil_basis": [encode(omega_a), encode(omega_b)],
        "Omega_plus": encode(omega_plus),
        "Omega_minus": encode(omega_minus),
        "multipliers": {"Omega_plus_G": "I", "Omega_minus_G": "-I"},
        "determinants": {"Omega_plus": "3 - 4*I", "Omega_minus": "3 + 4*I"},
        "rank": 4,
        "normalized_G": "zeta8^(-1) G for Omega_plus",
        "normalized_G_order": 8,
        "vector5": "(Lambda2(V)/Omega_inverse) tensor chi^(-1) = sign tensor Perm5",
        "adjoint10": "Sym2(V) tensor chi^(-1) = Lambda2(vector5)",
        "vector_characters": vector_rows,
        "adjoint_characters": adjoint_rows,
        "checks": checks,
        "scope": (
            "Exact structure on the v419 cyclotomic four. The original v419 "
            "source does not itself choose the central zeta8 scalar and this "
            "does not identify the separate NS/Pauli Bell four."
        ),
    }
    out = Path(__file__).with_name("f20-conformal-symplectic-results.json")
    out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "verdict": output["verdict"],
        "checks": len(checks), "rank": output["rank"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
