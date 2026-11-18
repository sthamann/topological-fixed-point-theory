#!/usr/bin/env python3
"""Exact joint-clock typecheck for the native F20 four and v634 sigma.

This checker keeps three representations separate:

* the selected native Pauli F20 four;
* the v634 native/theta family clock sigma on the same complex four;
* the original arithmetic Q(zeta30) module, where carrier and family clocks
  coexist on rational rank eight.
"""

from __future__ import annotations

import hashlib
import json
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
HERE = Path(__file__).resolve().parent
PINS = {
    "verification/v315_coxeter_coupling.py":
        "490969af0827cf8829f236e3999710bb526396dc9eb2ee0672d8bdd5c6af7b0c",
    "verification/v417_eisenstein_cp_operator.py":
        "c5f1fe26b9e6fd30bf2ad45354de6819a184c4a86b491f19ab21cb53b21c2b7e",
    "verification/v419_seam_galois_carrier.py":
        "cc6152ad7ea68c67725020e386e49eec32810ab25bb483c76bae99964399c2fd",
    "verification/v634_st31_structure.py":
        "5ff53ed031ca03103187ac22edf723aabce88bee3df1085c52486ecb0044d37c",
    "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/"
    "PAULI-SPIN-CENTRAL-LIFT.txt":
        "e5ef50c899eba80dca03e1b7023edf816df3b45a669534fccf2fef215b0ae8e6",
}


def encode(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(value)) for value in row] for row in matrix.tolist()]


def main() -> None:
    checks: list[str] = []

    def require(condition: bool, label: str) -> None:
        if not bool(condition):
            raise AssertionError(label)
        checks.append(label)

    for relative, expected in PINS.items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, f"source pin {relative}")

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
    H1 = sp.kronecker_product(H, identity2)
    H2 = sp.kronecker_product(identity2, H)
    S1 = sp.kronecker_product(S, identity2)
    S2 = sp.kronecker_product(identity2, S)

    c_native = sp.simplify(-ii * H2 * CNOT * S1 * H1)
    g_native = sp.simplify(sp.kronecker_product(Y, Z) * CNOT * S2 * S1)
    sigma = sp.Matrix([
        [0, 0, 1, 0],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
    ])
    require(c_native**5 == identity4, "selected native carrier has order five")
    require(g_native**4 == identity4 and
            sp.simplify(g_native * c_native * g_native.inv() - c_native**2) ==
            sp.zeros(4), "selected native F20 relations")
    require(sigma**3 == identity4 and sp.trace(sigma) == 1,
            "actual v634 sigma has order three and trace one")
    require(sigma.charpoly().as_expr().factor() ==
            (sigma.charpoly().gen - 1)**2 *
            (sigma.charpoly().gen**2 + sigma.charpoly().gen + 1),
            "actual v634 sigma characteristic polynomial")

    commutator = sp.simplify(c_native * sigma - sigma * c_native)
    g_plus = sp.simplify(g_native * sigma * g_native.inv() - sigma)
    g_minus = sp.simplify(g_native * sigma * g_native.inv() - sigma.inv())
    require(commutator.rank() == 4, "actual selected C5 does not commute with sigma")
    require(g_plus.rank() == 4, "actual selected G does not fix sigma")
    require(g_minus.rank() == 4, "actual selected G does not invert sigma")

    # The two order-eight quarter lifts have different projective type.
    zeta8 = (1 + ii) / sqrt2
    g_hat = sp.simplify(g_native / zeta8)
    g_hat_square = sp.simplify(g_hat**2)
    require(g_hat**4 == -identity4 and g_hat**8 == identity4,
            "F20 symplectic normalization has order eight")
    require(g_hat_square != ii * identity4 and
            g_hat_square != -ii * identity4 and
            sp.trace(g_hat_square) == 0,
            "F20 quarter lift has nonscalar square and projective order four")

    # Abstract four-dimensional obstruction under either extra joint relation.
    # In the C5 eigenbasis, G cycles the four one-dimensional eigenspaces.
    # A commuting order-three sigma is diag(omega^e_j).  G sigma G^-1=sigma^s
    # becomes e_{j+1}=s e_j mod 3, with s=+1 or -1.
    allowed: dict[str, list[dict[str, object]]] = {"plus": [], "minus": []}
    for label, sign in (("plus", 1), ("minus", -1)):
        for exponents in product(range(3), repeat=4):
            if all(exponents[(j + 1) % 4] ==
                   (sign * exponents[j]) % 3 for j in range(4)):
                counts = [exponents.count(k) for k in range(3)]
                # Trace = n0 + n1 omega + n2 omega^2; omega^2=-1-omega.
                trace_pair = (counts[0] - counts[2], counts[1] - counts[2])
                allowed[label].append({
                    "exponents": list(exponents),
                    "trace_in_basis_1_omega": list(trace_pair),
                })
    require(all(item["trace_in_basis_1_omega"] != [1, 0]
                for rows in allowed.values() for item in rows),
            "no commuting F20 four sigma has order three trace one under plus or minus")
    require({tuple(item["trace_in_basis_1_omega"])
             for item in allowed["plus"]} ==
            {(4, 0), (0, 4), (-4, -4)},
            "plus relation forces scalar sigma")
    require({tuple(item["trace_in_basis_1_omega"])
             for item in allowed["minus"]} == {(4, 0), (-2, 0)},
            "minus relation gives trace four or minus two")

    # The already original arithmetic joint module Q(zeta30), rational rank 8.
    x = sp.symbols("x")
    phi30 = sp.cyclotomic_poly(30, x)

    def column(polynomial: sp.Expr) -> sp.Matrix:
        remainder = sp.Poly(sp.rem(polynomial, phi30, x), x)
        return sp.Matrix([remainder.coeff_monomial(x**j) for j in range(8)])

    def multiplication(exponent: int) -> sp.Matrix:
        return sp.Matrix.hstack(*(column(x**(exponent + j)) for j in range(8)))

    c5_q30 = multiplication(6)       # zeta30^6=zeta5
    sigma_q30 = multiplication(10)  # zeta30^10=zeta3
    g_q30 = sp.Matrix.hstack(*(column(x**(7 * j)) for j in range(8)))
    identity8 = sp.eye(8)
    require(c5_q30**5 == identity8 and sigma_q30**3 == identity8 and
            g_q30**4 == identity8, "Q(zeta30) carrier family and Galois orders")
    require(c5_q30 * sigma_q30 == sigma_q30 * c5_q30,
            "Q(zeta30) cyclic carrier and family clocks commute")
    require(sp.simplify(g_q30 * c5_q30 * g_q30.inv() - c5_q30**2) ==
            sp.zeros(8), "Q(zeta30) carrier Frobenius relation")
    require(g_q30 * sigma_q30 == sigma_q30 * g_q30,
            "Q(zeta30) carrier Galois generator fixes the family clock")
    require(sp.trace(sigma_q30) == -4,
            "Q(zeta30) family trace is minus four not native trace one")

    result = {
        "checks_pass": True,
        "checks": checks,
        "source_pins": PINS,
        "actual_native_four": {
            "sigma": encode(sigma),
            "sigma_trace": str(sp.trace(sigma)),
            "sigma_charpoly": str(sp.factor(sigma.charpoly().as_expr())),
            "rank_C5_sigma_commutator": commutator.rank(),
            "rank_G_sigma_Ginv_minus_sigma": g_plus.rank(),
            "rank_G_sigma_Ginv_minus_sigma_inv": g_minus.rank(),
        },
        "quarter_lift_distinction": {
            "old_theta_alpha_square": "+i I4",
            "old_theta_alpha_projective_order": 2,
            "new_F20_Ghat_square": encode(g_hat_square),
            "new_F20_Ghat_square_trace": str(sp.trace(g_hat_square)),
            "new_F20_Ghat_projective_order": 4,
        },
        "abstract_four_dimensional_joint_relation": allowed,
        "original_Qzeta30_rank_eight": {
            "minimal_polynomial": str(phi30),
            "carrier_family_commute": True,
            "G_C5_Ginv_equals_C5_squared": True,
            "G_sigma_Ginv_equals_sigma": True,
            "sigma_trace": str(sp.trace(sigma_q30)),
        },
        "verdict": "JOINT_CLOCK_GAP_ON_SELECTED_NATIVE_F20_FOUR",
        "scope": (
            "The gap is for the selected native F20 four with the actual v634 "
            "sigma and for any irreducible F20 four under either stated extra "
            "joint relation. The original Q(zeta30) arithmetic rank-eight module "
            "does carry the commuting carrier and family clocks. No physical "
            "no-go or identification of the theta alpha with the F20 Ghat follows."
        ),
    }
    output = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    (HERE / "native-f20-joint-clock-results.json").write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
