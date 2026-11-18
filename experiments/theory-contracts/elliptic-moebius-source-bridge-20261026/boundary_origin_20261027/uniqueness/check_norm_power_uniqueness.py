#!/usr/bin/env python3
"""Exact residue check for the marked q norm-power family.

All phases are exponents of zeta_8=exp(pi i/4), so the calculation is
integer arithmetic modulo eight.  No floating-point comparison is used.
"""

from __future__ import annotations

import json
from pathlib import Path


OUT = Path(__file__).with_name("norm-power-uniqueness-results.json")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def phase_table(alpha_on_norm: int) -> list[dict[str, int | bool]]:
    """Return the W_k central data for one of the two spin-lift signs."""
    require(alpha_on_norm % 2 == 1, "alpha phase must have odd zeta_8 exponent")
    require(
        (2 * alpha_on_norm) % 8 == 6,
        "alpha_N^2 must equal -i=zeta_8^6",
    )
    rows = []
    for k in range(8):
        beta_exp = (2 * k) % 8  # beta on N^{-k} is i^k
        alpha_relative_exp = (-alpha_on_norm * k) % 8
        rows.append(
            {
                "k_mod_8": k,
                "beta_exp_mod_8": beta_exp,
                "alpha_relative_exp_mod_8": alpha_relative_exp,
                "native_beta_plus_i": beta_exp == 2,
                "alpha_fourth_is_minus_one": (2 * beta_exp) % 8 == 4,
            }
        )
    return rows


def main() -> None:
    # The two choices differ by alpha -> -alpha.
    tables = {str(a): phase_table(a) for a in (3, 7)}

    for rows in tables.values():
        require(
            [r["k_mod_8"] for r in rows if r["native_beta_plus_i"]] == [1, 5],
            "native beta phase must select k=1,5 modulo eight",
        )
        require(
            all(r["alpha_fourth_is_minus_one"] for r in rows if r["native_beta_plus_i"]),
            "selected alpha must have fourth power -1",
        )

        by_k = {int(r["k_mod_8"]): r for r in rows}
        for k in range(4):
            # Adding four norm copies flips alpha and leaves beta unchanged.
            require(
                by_k[k + 4]["beta_exp_mod_8"] == by_k[k]["beta_exp_mod_8"],
                "adding four copies must preserve beta",
            )
            require(
                (
                    int(by_k[k + 4]["alpha_relative_exp_mod_8"])
                    - int(by_k[k]["alpha_relative_exp_mod_8"])
                )
                % 8
                == 4,
                "adding four copies must flip alpha",
            )

    # eta is a local coordinate at each q point.  Q^*deta=i deta there,
    # hence Omega=product_q deta has phase i^3=-i=zeta_8^6.
    omega_q_exp = 6
    require((4 * omega_q_exp) % 8 == 0, "Omega^4 must be Q invariant")

    # The sign eta -> -eta changes Omega by -1, but not Omega^4.
    eta_sign_on_omega_exp = 4
    require(
        (4 * eta_sign_on_omega_exp) % 8 == 0,
        "Omega^4 must be independent of eta sign",
    )

    # alpha_N has exact order eight for either spin-lift sign.  Thus no
    # smaller positive tensor power of N is equivariantly trivial.
    for a in (3, 7):
        require(
            next(m for m in range(1, 9) if (m * a) % 8 == 0) == 8,
            "the norm character must have exact order eight",
        )

    result = {
        "phase_convention": "zeta_8 exponent modulo 8",
        "alpha_on_norm_choices": [3, 7],
        "alpha_on_norm_square": 6,
        "beta_on_inverse_norm_power": "2*k mod 8",
        "native_allowed_k_mod_8": [1, 5],
        "native_allowed_k_mod_4": [1],
        "k_plus_4_effect": "same beta and gamma, alpha changes sign",
        "k_plus_8_effect": "same full finite action",
        "normK_omega_Q_exp": omega_q_exp,
        "normK_omega_fourth_Q_exp": (4 * omega_q_exp) % 8,
        "eta_sign_on_omega_fourth_exp": (4 * eta_sign_on_omega_exp) % 8,
        "minimal_equivariantly_trivial_norm_power": 8,
        "gamma_on_norm_power": 1,
        "tables": tables,
        "verdict": "FULL_RELATIONS_FIX_K_MOD_4_NOT_LITERAL_ONE_COPY",
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    OUT.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
