#!/usr/bin/env python3
"""Exact F20 comparison of the original cyclotomic four and the two tens."""

from __future__ import annotations

from itertools import combinations, combinations_with_replacement
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
SOURCES = {
    "verification/v419_seam_galois_carrier.py":
        "cc6152ad7ea68c67725020e386e49eec32810ab25bb483c76bae99964399c2fd",
    "experiments/theory-contracts/compiler-current-descendant-20260919/checker.py":
        "b342a8b102d9b1dedf150ccbd4812ccb26fa1b9a4b7cb70401b6cdd8201b4420",
}


def require(value: bool, label: str, checks: list[str]) -> None:
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)


def symmetric_square(matrix: sp.Matrix) -> sp.Matrix:
    pairs = list(combinations_with_replacement(range(matrix.rows), 2))
    columns = []
    for a, b in pairs:
        tensor = sp.zeros(matrix.rows)
        tensor[a, b] = 1
        tensor[b, a] = 1
        image = matrix * tensor * matrix.T
        columns.append(sp.Matrix([image[i, j] for i, j in pairs]))
    return sp.Matrix.hstack(*columns)


def permutation_matrix(images: list[int]) -> sp.Matrix:
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


def encode(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix.tolist()]


def main() -> None:
    checks: list[str] = []
    for relative, expected in SOURCES.items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, "source pin " + relative, checks)

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
    require(galois * c5 * galois.inv() == c5**2,
            "original Frobenius relation", checks)
    require(c5**5 == sp.eye(4) and galois**4 == sp.eye(4),
            "original F20 orders", checks)

    # Affine action x -> x+1 and x -> 2x on the five carrier slots.
    p5 = permutation_matrix([(j + 1) % 5 for j in range(5)])
    pg = permutation_matrix([(2 * j) % 5 for j in range(5)])
    require(pg * p5 * pg.inv() == p5**2,
            "permutation F20 relation", checks)

    native_c = symmetric_square(c5)
    native_g = symmetric_square(galois)
    hodge_c = exterior_square(p5)
    hodge_g = exterior_square(pg)

    # An exact invertible solution of
    # T (chi tensor Sym2)(a) = Lambda2(Perm5)(a) T,
    # with chi(C5)=1 and chi(G)=i.
    intertwiner = sp.Matrix([
        [0,0,0,0,0,0,0,0,1,-1],
        [0,0,0,0,sp.I,0,-sp.I,0,0,0],
        [0,0,0,0,0,-sp.I,0,0,0,0],
        [0,0,0,0,0,-1,1,1,-1,0],
        [0,0,0,-1,0,0,0,0,0,0],
        [0,0,-sp.I,0,0,0,0,sp.I,0,0],
        [-sp.I,0,sp.I,sp.I,0,0,0,0,-sp.I,0],
        [-1,1,0,0,0,0,0,0,0,0],
        [0,sp.I,0,-sp.I,0,0,-sp.I,0,0,sp.I],
        [0,1,-1,0,-1,1,0,0,0,0],
    ])
    require(intertwiner.det() == 3 * sp.I,
            "intertwiner invertible determinant 3i", checks)
    require(intertwiner * native_c == hodge_c * intertwiner,
            "C5 intertwining", checks)
    require(intertwiner * (sp.I * native_g) == hodge_g * intertwiner,
            "quarter-character G intertwining", checks)

    character_rows = []
    for b in range(4):
        for a in range(5):
            native = sp.trace(native_c**a * native_g**b)
            twisted = sp.I**b * native
            hodge = sp.trace(hodge_c**a * hodge_g**b)
            require(sp.simplify(twisted - hodge) == 0,
                    f"character equality a={a} b={b}", checks)
            character_rows.append({
                "a": a, "b": b, "native": str(native),
                "chi": str(sp.I**b), "twisted": str(twisted),
                "hodge": str(hodge),
            })

    require(sp.trace(native_g**2) == 2 and sp.trace(hodge_g**2) == -2,
            "untwisted involution character mismatch plus2 minus2", checks)
    require(sp.I**2 * sp.trace(native_g**2) == sp.trace(hodge_g**2),
            "quarter character repairs involution class", checks)

    output = {
        "verdict": "EXACT_F20_TWISTED_MODULE_EQUIVALENCE",
        "group": "F20=C5 semidirect C4, G C5 G^-1=C5^2",
        "native": "Sym2(Q(zeta5) cyclotomic four)",
        "target": "Lambda2(permutation five)",
        "character": {"chi_C5": "1", "chi_G": "I"},
        "basis_native": [str(p) for p in combinations_with_replacement(range(4), 2)],
        "basis_target": [str(p) for p in combinations(range(5), 2)],
        "C5": encode(c5), "G": encode(galois),
        "intertwiner": encode(intertwiner),
        "intertwiner_determinant": str(intertwiner.det()),
        "characters": character_rows,
        "checks": checks,
        "scope": (
            "Exact finite F20 representation dictionary. It does not identify "
            "the old Bell channel with the Hodge channel and does not supply "
            "the still absent physical identification of the Bell fundamental "
            "C4 with the cyclotomic carrier four."
        ),
    }
    out = Path(__file__).with_name("f20-twisted-intertwiner-results.json")
    out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "verdict": output["verdict"],
        "checks": len(checks),
        "det": output["intertwiner_determinant"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
