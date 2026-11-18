#!/usr/bin/env python3
"""Exact E8 clock-twisted sectors for the native quarter turn.

This checker uses only integer and rational arithmetic.  In doubled standard
coordinates the E8 lattice consists of vectors y whose coordinates have one
parity and whose coordinate sum is divisible by four.  The inner clock

    kappa = exp(2*pi*i*v(0)),  v=(1/2^5,0^3),

has kappa^a-twisted charges E8-a*v in the Li Delta-operator convention.
The already proved full-field inner conjugacy kappa ~ Jhat is an input; this
file does not reconstruct that conjugacy or a physical clock.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import comb, gcd, isqrt
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "clock-twisted-results.json"
MAX_CHARACTER_DEGREE = 24  # degree in x=q^(1/4), after removing x^(1/6)


def e8_shift_vectors(a: int, max_square: int) -> Iterable[tuple[int, ...]]:
    """Yield z=2*(beta-a*v) with beta in E8 and sum(z_i^2)<=max_square."""

    radius = isqrt(max_square)

    def rec(prefix: tuple[int, ...], square: int) -> Iterable[tuple[int, ...]]:
        if len(prefix) == 8:
            y = tuple(value + (a if i < 5 else 0) for i, value in enumerate(prefix))
            if len({value & 1 for value in y}) == 1 and sum(y) % 4 == 0:
                yield prefix
            return
        remaining = max_square - square
        local_radius = min(radius, isqrt(remaining))
        for value in range(-local_radius, local_radius + 1):
            yield from rec(prefix + (value,), square + value * value)

    yield from rec((), 0)


def t_chart(last_three: tuple[int, int, int]) -> tuple[Fraction, ...]:
    """D3 to A3 isometry on doubled input coordinates."""

    x, y, z = (Fraction(value, 2) for value in last_three)
    return (
        (x + y + z) / 2,
        (x - y - z) / 2,
        (-x + y - z) / 2,
        (-x - y + z) / 2,
    )


def frac(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def shell_data(a: int) -> dict[str, object]:
    candidates = list(e8_shift_vectors(a, 12))
    minimum = min(sum(value * value for value in vector) for vector in candidates)
    shell = sorted(vector for vector in candidates if sum(value * value for value in vector) == minimum)
    split = Counter(
        (sum(value * value for value in vector[:5]), sum(value * value for value in vector[5:]))
        for vector in shell
    )
    a3_weights = sorted(t_chart(vector[5:]) for vector in shell)

    expected_representations = {
        0: "(1,1)",
        1: "(1,4)",
        2: "(10,1)+(1,6)",
        3: "(1,bar4)",
    }
    expected_splits = {
        0: {(0, 0): 1},
        1: {(0, 3): 4},
        2: {(4, 0): 10, (0, 4): 6},
        3: {(0, 3): 4},
    }
    expected_dimensions = (1, 4, 16, 4)
    expected_minima = (0, 3, 4, 3)
    zero5 = (0,) * 5
    expected_shells = {
        0: {(0,) * 8},
        1: {
            zero5 + last
            for last in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
        },
        2: {
            tuple(sign * 2 if coordinate == index else 0 for coordinate in range(8))
            for index in range(8)
            for sign in (-1, 1)
        },
        3: {
            zero5 + last
            for last in ((-1, -1, -1), (-1, 1, 1), (1, -1, 1), (1, 1, -1))
        },
    }

    return {
        "a": a,
        "minimum_norm_squared": frac(Fraction(minimum, 4)),
        "ground_weight": frac(Fraction(minimum, 8)),
        "ground_dimension": len(shell),
        "D5_A3_representation": expected_representations[a],
        "doubled_norm_split_counts": {
            f"{left}+{right}": count for (left, right), count in sorted(split.items())
        },
        "shifted_shell_doubled": [list(vector) for vector in shell],
        "A3_weights": [[frac(value) for value in weight] for weight in a3_weights],
        "checks": {
            "minimum": minimum == expected_minima[a],
            "dimension": len(shell) == expected_dimensions[a],
            "representation_split": dict(split) == expected_splits[a],
            "exact_D5_D3_weight_set": set(shell) == expected_shells[a],
        },
    }


def multiply_inverse_power(series: list[int], degree: int, power: int) -> list[int]:
    result = series[:]
    for total in range(len(series)):
        if series[total] == 0:
            continue
        for copies in range(1, (len(series) - 1 - total) // degree + 1):
            result[total + copies * degree] += series[total] * comb(copies + power - 1, power - 1)
    return result


def oscillator_character(max_degree: int) -> list[int]:
    series = [0] * (max_degree + 1)
    series[0] = 4
    for degree in range(1, max_degree + 1, 4):
        series = multiply_inverse_power(series, degree, 4)
    for degree in range(3, max_degree + 1, 4):
        series = multiply_inverse_power(series, degree, 4)
    return series


def shifted_theta_eta_character(max_degree: int) -> tuple[list[int], dict[int, int]]:
    # Relative x-degree is (sum(z_i^2)-3)/2 for x=q^(1/4).
    max_square = 3 + 2 * max_degree
    theta: Counter[int] = Counter()
    for vector in e8_shift_vectors(1, max_square):
        square = sum(value * value for value in vector)
        difference = square - 3
        if difference >= 0 and difference % 2 == 0:
            degree = difference // 2
            if degree <= max_degree:
                theta[degree] += 1

    series = [theta.get(degree, 0) for degree in range(max_degree + 1)]
    for degree in range(4, max_degree + 1, 4):
        series = multiply_inverse_power(series, degree, 8)
    return series, dict(sorted(theta.items()))


def phase_data() -> dict[str, object]:
    rho = Fraction(3, 8)
    central_shift = Fraction(8, 24)
    exponents = [((rho + Fraction(k, 4) - central_shift) % 1) for k in range(4)]
    orders = [value.denominator // gcd(value.numerator, value.denominator) for value in exponents]
    orbifold_type = (16 * rho.numerator // rho.denominator) % 4
    discriminant_level = 16 // gcd(orbifold_type, 4)
    modular_level = 3 * discriminant_level  # gcd(3,8)=1
    return {
        "central_charge": 8,
        "twisted_ground_weight": "3/8",
        "ground_rotation_exponent": "3/8",
        "ground_rotation_order": 8,
        "normalized_T_exponents_mod_1": [frac(value) for value in exponents],
        "normalized_T_orders": orders,
        "T_fourth_power_exponent": "1/6",
        "orbifold_type": "4{2}",
        "discriminant_form_level": discriminant_level,
        "twisted_trace_congruence_level": modular_level,
        "checks": {
            "type_is_4_2": orbifold_type == 2,
            "all_normalized_T_phases_have_order_24": orders == [24, 24, 24, 24],
            "modular_level_is_24": modular_level == 24,
        },
    }


def build_result() -> dict[str, object]:
    sectors = [shell_data(a) for a in range(4)]
    oscillator = oscillator_character(MAX_CHARACTER_DEGREE)
    theta_eta, theta = shifted_theta_eta_character(MAX_CHARACTER_DEGREE)
    character_equal = oscillator == theta_eta
    phase = phase_data()
    checks = {
        "all_shell_checks": all(all(item["checks"].values()) for item in sectors),
        "character_coefficients_equal": character_equal,
        "phase_checks": all(phase["checks"].values()),
    }
    return {
        "id": "native-clock-twisted-sectors",
        "status": "exact" if all(checks.values()) else "failed",
        "assumption": "existing full-field inner conjugacy kappa~Jhat preserving the conformal vector",
        "sectors": sectors,
        "phases": phase,
        "character_identity": {
            "variable": "x=q^(1/4), common prefactor x^(1/6)=q^(1/24) removed",
            "through_x_degree": MAX_CHARACTER_DEGREE,
            "oscillator_coefficients": oscillator,
            "theta_eta_coefficients": theta_eta,
            "shifted_theta_relative_shell_counts": {str(key): value for key, value in theta.items()},
            "equal": character_equal,
        },
        "checks": checks,
        "scope": [
            "The level 24 is a modular congruence level and an order of normalized T phases.",
            "No C12 or C24 base cover and no physical time are inferred.",
            "The A3 labels belong to the toral kappa chart and transport with the conjugator to Jhat.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build_result()
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")
    for name, ok in result["checks"].items():
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    print(f"wrote {args.output}")
    return 0 if result["status"] == "exact" else 1


if __name__ == "__main__":
    raise SystemExit(main())
