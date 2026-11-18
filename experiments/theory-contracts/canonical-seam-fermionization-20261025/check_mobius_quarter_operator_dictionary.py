#!/usr/bin/env python3
"""Exact reduced-longitudinal checks for the Moebius NS-cover quarters.

The full one-step Pin-minus monodromy is cJ with J including transverse
reflection. This checker begins after projection to the two J eigenspaces;
it does not replace the full pinor bundle action by a scalar translation.
"""

from fractions import Fraction


def need(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"PASS {label}")


def cover_label(epsilon: int, m: int) -> int:
    return 2 * m + epsilon


def cover_eigenvalue(n: int) -> Fraction:
    return Fraction(n, 1) + Fraction(1, 2)


def quarter_eigenvalue(epsilon: int, m: int) -> Fraction:
    a = Fraction(1, 4) if epsilon == 0 else Fraction(3, 4)
    return Fraction(m, 1) + a


def u_action(epsilon: int, m: int) -> tuple[int, int]:
    # [[0,u_base],[I,0]] in the ordered (1/4,3/4) coefficient basis.
    return (1, m) if epsilon == 0 else (0, m + 1)


def crossing_count(values: list[Fraction], shift: Fraction) -> int:
    return sum(1 for value in values if value < 0 < value + shift)


def run() -> None:
    labels = [(epsilon, m) for epsilon in (0, 1) for m in range(-8, 9)]

    need(
        all(
            cover_eigenvalue(cover_label(epsilon, m))
            == 2 * quarter_eigenvalue(epsilon, m)
            for epsilon, m in labels
        ),
        "B_cov = 2(B_1/4 direct-sum B_3/4)",
    )

    need(
        all(
            ((1 if cover_label(epsilon, m) % 2 == 0 else -1)
             == (1 if epsilon == 0 else -1))
            for epsilon, m in labels
        ),
        "half-deck eigenphases are +i on 1/4 and -i on 3/4",
    )

    need(
        all(
            cover_label(*u_action(epsilon, m))
            == cover_label(epsilon, m) + 1
            for epsilon, m in labels
        ),
        "u_cov has block action [[0,u_base],[I,0]]",
    )

    need(
        all(
            u_action(*u_action(epsilon, m)) == (epsilon, m + 1)
            for epsilon, m in labels
        ),
        "u_cov squared is diag(u_base,u_base)",
    )

    # Each W_epsilon maps an orthonormal Fourier basis to a disjoint parity
    # subset of the NS-cover Fourier basis.  A coefficient norm is therefore
    # preserved exactly, and the two images are orthogonal.
    coeffs = {
        (0, -2): Fraction(3, 5),
        (0, 1): Fraction(-7, 11),
        (1, -1): Fraction(5, 13),
        (1, 3): Fraction(2, 7),
    }
    source_norm = sum(value * value for value in coeffs.values())
    cover_coeffs = {
        cover_label(epsilon, m): value
        for (epsilon, m), value in coeffs.items()
    }
    target_norm = sum(value * value for value in cover_coeffs.values())
    need(
        len(cover_coeffs) == len(coeffs) and source_norm == target_norm,
        "quarter-to-cover map is isometric and sector images are orthogonal",
    )

    cover_values = [cover_eigenvalue(n) for n in range(-8, 9)]
    need(
        crossing_count(cover_values, Fraction(1, 1)) == 1,
        "u_cov winding one has full-cover spectral flow one",
    )
    need(
        crossing_count(cover_values, Fraction(2, 1)) == 2,
        "u_cov squared has full-cover spectral flow two",
    )

    for epsilon, name in ((0, "1/4"), (1, "3/4")):
        quarter_values = [quarter_eigenvalue(epsilon, m) for m in range(-8, 9)]
        need(
            crossing_count(quarter_values, Fraction(1, 1)) == 1,
            f"u_base has spectral flow one on the {name} domain",
        )

    print("RESULT exact Moebius quarter operator dictionary verified")


if __name__ == "__main__":
    run()
