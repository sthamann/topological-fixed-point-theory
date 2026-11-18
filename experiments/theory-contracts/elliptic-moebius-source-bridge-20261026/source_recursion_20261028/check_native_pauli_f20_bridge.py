#!/usr/bin/env python3
"""Exact native Pauli F20 subgroup and unitary bridge to the v419 four."""

from __future__ import annotations

from itertools import combinations_with_replacement, product
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
V419_PIN = "cc6152ad7ea68c67725020e386e49eec32810ab25bb483c76bae99964399c2fd"
V783_PIN = "8f4851634b83f61671b04f3a6211059c40758d21caaf1f779302c799e1e9d6c4"


def require(value: bool, label: str, checks: list[str]) -> None:
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)


def is_zero(matrix: sp.Matrix) -> bool:
    return sp.simplify(matrix) == sp.zeros(*matrix.shape)


def encode(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(value)) for value in row] for row in matrix.tolist()]


def sym2_character(matrix: sp.Matrix) -> sp.Expr:
    return sp.simplify((sp.trace(matrix) ** 2 + sp.trace(matrix**2)) / 2)


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


def induced_symmetric_gram(gram: sp.Matrix) -> sp.Matrix:
    pairs = list(combinations_with_replacement(range(gram.rows), 2))
    vectors = []
    for a, b in pairs:
        vector = sp.zeros(gram.rows**2, 1)
        vector[a * gram.rows + b] = 1
        if a != b:
            vector[b * gram.rows + a] = 1
        vectors.append(vector)
    basis = sp.Matrix.hstack(*vectors)
    return sp.simplify(basis.conjugate().T * sp.kronecker_product(gram, gram) * basis)


def main() -> None:
    checks: list[str] = []
    v419_path = ROOT / "verification/v419_seam_galois_carrier.py"
    v783_path = ROOT / "verification/v783_two_qubit_clifford.py"
    require(hashlib.sha256(v419_path.read_bytes()).hexdigest() == V419_PIN,
            "v419 source pin", checks)
    require(hashlib.sha256(v783_path.read_bytes()).hexdigest() == V783_PIN,
            "v783 source pin", checks)

    ii = sp.I
    sqrt2 = sp.sqrt(2)
    identity2 = sp.eye(2)
    identity4 = sp.eye(4)
    X = sp.Matrix([[0, 1], [1, 0]])
    Y = sp.Matrix([[0, -ii], [ii, 0]])
    Z = sp.diag(1, -1)
    H = sp.Matrix([[1, 1], [1, -1]]) / sqrt2
    S = sp.diag(1, ii)
    CNOT = sp.Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
    ])
    H1, H2 = sp.kronecker_product(H, identity2), sp.kronecker_product(identity2, H)
    S1, S2 = sp.kronecker_product(S, identity2), sp.kronecker_product(identity2, S)

    # Rightmost gate acts first. These are words in the frozen v783 source
    # generators, followed only by source-proved Pauli/central corrections.
    c_raw = sp.simplify(H2 * CNOT * S1 * H1)
    g_raw = sp.simplify(CNOT * S2 * S1)
    c_native = sp.simplify(-ii * c_raw)
    g_native = sp.simplify(sp.kronecker_product(Y, Z) * g_raw)

    require(is_zero(c_raw**5 - ii * identity4), "raw Singer word fifth power i", checks)
    require(is_zero(c_native**5 - identity4), "normalized native C has order five", checks)
    require(is_zero(g_native**4 - identity4), "Pauli corrected native G has order four", checks)
    require(is_zero(g_native * c_native * g_native.inv() - c_native**2),
            "exact native Frobenius relation", checks)
    require(sp.trace(c_native) == -1 and sp.trace(g_native) == 0 and
            sp.trace(g_native**2) == 0,
            "native four has F20 character generators", checks)

    paulis = {name1 + name2: sp.kronecker_product(p1, p2)
              for name1, p1 in (("I", identity2), ("X", X), ("Y", Y), ("Z", Z))
              for name2, p2 in (("I", identity2), ("X", X), ("Y", Y), ("Z", Z))}
    zeta8 = (1 + ii) / sqrt2
    exact_lifts: list[dict[str, object]] = []
    for name, pauli in paulis.items():
        for k in range(8):
            candidate = sp.simplify(zeta8**k * pauli * g_raw)
            if is_zero(candidate**4 - identity4) and is_zero(
                    candidate * c_native * candidate.inv() - c_native**2):
                exact_lifts.append({"pauli": name, "zeta8_power": k})
    require(exact_lifts == [
        {"pauli": "YZ", "zeta8_power": 0},
        {"pauli": "YZ", "zeta8_power": 2},
        {"pauli": "YZ", "zeta8_power": 4},
        {"pauli": "YZ", "zeta8_power": 6},
    ], "all exact source lift corrections are the mu4 orbit of YZ", checks)

    # The actual v419 cyclotomic carrier.
    x = sp.symbols("x")
    phi5 = x**4 + x**3 + x**2 + x + 1
    c_v419 = sp.Matrix([
        [0, 0, 0, -1], [1, 0, 0, -1],
        [0, 1, 0, -1], [0, 0, 1, -1],
    ])

    def galois_column(exponent: int) -> sp.Matrix:
        polynomial = sp.Poly(sp.rem(x**exponent, phi5, x), x)
        return sp.Matrix([polynomial.coeff_monomial(x**j) for j in range(4)])

    g_v419 = sp.Matrix.hstack(*(galois_column(2 * k) for k in range(4)))
    field_trace_gram = 5 * identity4 - sp.ones(4)
    require(c_v419.T * field_trace_gram * c_v419 == field_trace_gram and
            g_v419.T * field_trace_gram * g_v419 == field_trace_gram,
            "v419 field trace Gram is F20 invariant", checks)
    require(field_trace_gram.det() == 125 and
            field_trace_gram.eigenvals() == {sp.Integer(1): 1, sp.Integer(5): 3},
            "v419 Gram is positive with eigenvalues one five five five", checks)

    intertwiner = sp.Matrix([
        [-1-ii, ii, ii, ii],
        [-1+ii, -2-ii, 1, 1],
        [1-ii, -1, -1, 2+ii],
        [-1+ii, 1, -2-ii, 1],
    ])
    require(is_zero(intertwiner * c_v419 - c_native * intertwiner) and
            is_zero(intertwiner * g_v419 - g_native * intertwiner),
            "exact four dimensional F20 intertwiner", checks)
    require(intertwiner.det() == -8 + 44*ii, "intertwiner determinant", checks)
    require(is_zero(intertwiner.conjugate().T * intertwiner - 2 * field_trace_gram),
            "intertwiner pulls standard metric back to twice field trace Gram", checks)

    # The solution space is one-dimensional, so the relative unitary is unique
    # up to a phase. U=T/sqrt(2) obeys U^dagger U=H, the correct notion of
    # unitarity from (V,H) to standard C^4.
    variables = sp.symbols("t0:16")
    unknown = sp.Matrix(4, 4, variables)
    equations = list(unknown * c_v419 - c_native * unknown)
    equations += list(unknown * g_v419 - g_native * unknown)
    coefficient, _ = sp.linear_eq_to_matrix(equations, variables)
    require(16 - coefficient.rank() == 1,
            "four dimensional intertwiner space has dimension one", checks)
    relative_unitary = sp.simplify(intertwiner / sqrt2)
    require(is_zero(relative_unitary.conjugate().T * relative_unitary - field_trace_gram),
            "T over sqrt two is a metric preserving isomorphism", checks)
    require(sp.simplify(abs(relative_unitary.det())**2) == 125,
            "relative unitary determinant matches Gram determinant", checks)

    character_rows: list[dict[str, object]] = []
    for b, a in product(range(4), range(5)):
        matrix = sp.simplify(c_native**a * g_native**b)
        value = sp.simplify(sp.trace(matrix))
        expected = 4 if (a, b) == (0, 0) else (-1 if b == 0 else 0)
        require(value == expected, f"native F20 character a={a} b={b}", checks)
        character_rows.append({"a": a, "b": b, "trace": str(value)})

    # The native source realizes the same conformal symplectic multiplier.
    standard_symplectic = sp.diag(1, 1, 1, 1)
    standard_symplectic[0, 0] = 0
    standard_symplectic[1, 1] = 0
    standard_symplectic[2, 2] = 0
    standard_symplectic[3, 3] = 0
    standard_symplectic[0, 1], standard_symplectic[1, 0] = 1, -1
    standard_symplectic[2, 3], standard_symplectic[3, 2] = 1, -1
    require(c_native.T * standard_symplectic * c_native == standard_symplectic,
            "native C preserves standard alternating form", checks)
    require(g_native.T * standard_symplectic * g_native == ii * standard_symplectic,
            "native G has conformal multiplier i", checks)
    g_spin = sp.simplify(zeta8**-1 * g_native)
    require(is_zero(g_spin.T * standard_symplectic * g_spin - standard_symplectic) and
            is_zero(g_spin**4 + identity4) and is_zero(g_spin**8 - identity4),
            "mu8 normalized native G is symplectic of order eight", checks)

    omega_a = sp.Matrix([
        [0, 0, 1, -1], [0, 0, 0, 1],
        [-1, 0, 0, 0], [1, -1, 0, 0],
    ])
    omega_b = sp.Matrix([
        [0, 1, 0, 0], [-1, 0, 1, 0],
        [0, -1, 0, 1], [0, 0, -1, 0],
    ])
    omega_plus = ii * omega_a + omega_b
    transported = sp.simplify(intertwiner.inv().T * omega_plus * intertwiner.inv())
    require(is_zero(transported - (1 - 2*ii) * standard_symplectic / 10),
            "intertwiner transports v419 Omega plus to standard form", checks)

    # On chi^(-1) Sym2(V), the ten decomposes as 2V4 + chi_1 + chi_3.
    rep10_character: dict[tuple[int, int], sp.Expr] = {}
    for b, a in product(range(4), range(5)):
        matrix = sp.simplify(c_native**a * g_native**b)
        rep10_character[(a, b)] = sp.simplify(ii**(-b) * sym2_character(matrix))
    one_dimensional_multiplicities: dict[str, str] = {}
    for k in range(4):
        multiplicity = sp.simplify(sum(
            sp.conjugate(ii**(k*b)) * rep10_character[(a, b)]
            for b, a in product(range(4), range(5))) / 20)
        one_dimensional_multiplicities[f"chi_{k}"] = str(multiplicity)
    four_multiplicity = sp.simplify(sum(
        sp.conjugate(sp.trace(c_native**a * g_native**b)) *
        rep10_character[(a, b)]
        for b, a in product(range(4), range(5))) / 20)
    require(one_dimensional_multiplicities == {
        "chi_0": "0", "chi_1": "1", "chi_2": "0", "chi_3": "1"
    } and four_multiplicity == 2,
            "ten decomposes as two fours plus chi one plus chi three", checks)
    commutant_dimension = int(four_multiplicity**2 + 1 + 1)
    require(commutant_dimension == 6,
            "ten dimensional commutant has complex dimension six", checks)

    # The earlier determinant-3i rational ten-intertwiner is algebraic, but it
    # is not already an isometry for the field-trace metric. Its canonical
    # metric polar correction is governed by a positive commutant element.
    old_ten_intertwiner = sp.Matrix([
        [0,0,0,0,0,0,0,0,1,-1],
        [0,0,0,0,ii,0,-ii,0,0,0],
        [0,0,0,0,0,-ii,0,0,0,0],
        [0,0,0,0,0,-1,1,1,-1,0],
        [0,0,0,-1,0,0,0,0,0,0],
        [0,0,-ii,0,0,0,0,ii,0,0],
        [-ii,0,ii,ii,0,0,0,0,-ii,0],
        [-1,1,0,0,0,0,0,0,0,0],
        [0,ii,0,-ii,0,0,-ii,0,0,ii],
        [0,1,-1,0,-1,1,0,0,0,0],
    ])
    sym_gram = induced_symmetric_gram(field_trace_gram)
    ten_c = symmetric_square(c_v419)
    ten_g = ii * symmetric_square(g_v419)
    require(is_zero(ten_c.conjugate().T * sym_gram * ten_c - sym_gram) and
            is_zero(ten_g.conjugate().T * sym_gram * ten_g - sym_gram),
            "correct induced field trace metric is invariant on twisted ten", checks)
    require(not is_zero(old_ten_intertwiner.conjugate().T * old_ten_intertwiner - sym_gram),
            "old determinant 3i ten intertwiner is not metric preserving", checks)
    polar_square = sp.simplify(sym_gram.inv() *
                               old_ten_intertwiner.conjugate().T * old_ten_intertwiner)
    require(is_zero(polar_square * ten_c - ten_c * polar_square) and
            is_zero(polar_square * ten_g - ten_g * polar_square),
            "ten polar square lies in the F20 commutant", checks)

    output = {
        "verdict": "EXACT_NATIVE_SOURCE_F20_AND_METRIC_BRIDGE",
        "source_words_rightmost_first": {
            "C_raw": "(I tensor H) CNOT (S tensor I) (H tensor I)",
            "C": "-i C_raw",
            "G_raw": "CNOT (I tensor S) (S tensor I)",
            "G": "(Y tensor Z) G_raw",
        },
        "exact_G_lifts": exact_lifts,
        "C_native": encode(c_native),
        "G_native": encode(g_native),
        "v419_Gram": encode(field_trace_gram),
        "intertwiner_T": encode(intertwiner),
        "intertwiner_determinant": str(intertwiner.det()),
        "metric_identity": "T^dagger T = 2 (5 I - J)",
        "relative_unitary": "U=T/sqrt(2), U^dagger U=5 I-J",
        "intertwiner_space_dimension": 1,
        "characters": character_rows,
        "conformal_multiplier_G": "i",
        "spin_lift": "zeta8^(-1) G, order 8",
        "transported_alternating_form": "(1-2i)/10 times J_standard",
        "ten_decomposition": "2 V4 + chi_1 + chi_3",
        "ten_multiplicities": {
            **one_dimensional_multiplicities,
            "V4": str(four_multiplicity),
        },
        "ten_commutant": "M2(C) plus C plus C",
        "ten_commutant_dimension": commutant_dimension,
        "old_ten_intertwiner": "determinant 3i, not metric preserving",
        "ten_polar_repair": (
            "P^2=H_sym^(-1) A^dagger A; U=A P^(-1), with P the positive "
            "H_sym square root"
        ),
        "unitary_ten_residual_gauge": "U(2) times U(1) times U(1)",
        "checks": checks,
        "scope": (
            "The F20 is constructed from the frozen v783 native Clifford "
            "generators and source-proved Pauli/mu8 kernel. v783 does not "
            "designate this Singer normalizer as the v419 Galois clock. The "
            "result is a carrier representation bridge, not a physical "
            "four-dimensional field or channel identification."
        ),
    }
    out = Path(__file__).with_name("native-pauli-f20-results.json")
    out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "verdict": output["verdict"],
        "checks": len(checks),
        "ten_commutant_dimension": commutant_dimension,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
