#!/usr/bin/env python3
"""Exact color-center test for the two H8 ground charges.

The infinite-lattice ground classification is proved coordinatewise.  The
finite enumeration below is used only for the first positive center-neutral
energy.  Its cutoff is certified by explicit coordinate lower bounds.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json


HERE = Path(__file__).resolve().parent
REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
FOLD_RESULTS = REPO / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"

ALPHAS = (F(1, 4), F(3, 4))
DIM = 8
PSTAR = (F(1, 2),) * 4 + (F(-1, 2),) * 4
ZERO = (F(0),) * DIM


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fstr(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def delta(alpha: F) -> tuple[F, ...]:
    return tuple((F(k) + alpha) / 8 - F(1, 2) for k in range(DIM))


def energy_k(p: tuple[F, ...], alpha: F, oscillator_level: int = 0) -> F:
    d = delta(alpha)
    return F(oscillator_level) + sum(x * x for x in p) / 2 + sum(x * y for x, y in zip(d, p))


def matvec(matrix: tuple[tuple[F, ...], ...], vector: tuple[F, ...]) -> tuple[F, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def color_triality(carrier_weight: tuple[F, ...]) -> int:
    """SU(3) center exponent in the fixed first-three-slot carrier frame.

    The coweight pairing is q1+q2-2*q3 modulo 3.  It is integral on the
    restricted E8 weights used below.  On an integral root it agrees with the
    sum of the first three coordinates modulo 3.
    """

    numerator = carrier_weight[0] + carrier_weight[1] - 2 * carrier_weight[2]
    need(numerator.denominator == 1, f"nonintegral SU(3) coweight pairing: {numerator}")
    return int(numerator) % 3


def psing_eigenvalue(triality: int) -> int:
    """Exact cyclotomic average (1+omega^t+omega^(2t))/3."""

    return 1 if triality % 3 == 0 else 0


def e8_roots() -> list[tuple[F, ...]]:
    roots: list[tuple[F, ...]] = []
    for i, j in combinations(range(DIM), 2):
        for si, sj in product((-1, 1), repeat=2):
            root = [F(0)] * DIM
            root[i], root[j] = F(si), F(sj)
            roots.append(tuple(root))
    for signs in product((-1, 1), repeat=DIM):
        if sum(s < 0 for s in signs) % 2 == 0:
            roots.append(tuple(F(s, 2) for s in signs))
    need(len(roots) == 240 and len(set(roots)) == 240, "E8 root construction failed")
    need(all(sum(x * x for x in root) == 2 for root in roots), "non-root in E8 root list")
    return roots


def main() -> None:
    fold = json.loads(FOLD_RESULTS.read_text(encoding="utf-8"))
    matrix = tuple(tuple(F(x) for x in row) for row in fold["explicit_Weyl_matrix"])
    need(len(matrix) == DIM and all(len(row) == DIM for row in matrix), "Weyl matrix is not 8x8")
    need(all(sum(matrix[i][k] * matrix[j][k] for k in range(DIM)) == (i == j)
             for i in range(DIM) for j in range(DIM)), "Weyl matrix is not orthogonal")

    pstar_image = matvec(matrix, PSTAR)
    expected_image = (F(1), F(0), F(0), F(1), F(0), F(0), F(0), F(0))
    need(pstar_image == expected_image, f"unexpected M p*: {pstar_image}")

    # Original 3+2 carrier mark: X=6Y=(-2,-2,-2,3,3).
    xmark = (F(-2), F(-2), F(-2), F(3), F(3), F(0), F(0), F(0))
    xcharge = sum(a * b for a, b in zip(xmark, pstar_image))
    source_charge = tuple(F(x) for x in fold["transported_source_U1_charge"])
    source_u1 = sum(a * b for a, b in zip(source_charge, pstar_image))
    pstar_triality = color_triality(pstar_image)
    need(xcharge == 1, f"unexpected scaled hypercharge X: {xcharge}")
    need(source_u1 == 0, f"unexpected transported source U1 charge: {source_u1}")
    need(pstar_triality == 1, f"unexpected p* color triality: {pstar_triality}")
    need(psing_eigenvalue(pstar_triality) == 0, "P_sing failed to remove colored p*")
    need(psing_eigenvalue(color_triality(matvec(matrix, ZERO))) == 1, "P_sing removed vacuum")

    roots = e8_roots()
    cases: dict[str, object] = {}
    for alpha in ALPHAS:
        d = delta(alpha)
        need(all(abs(x) < F(1, 2) for x in d), "integer coordinate convexity premise failed")
        need(all(x < 0 for x in d[:4]) and all(x > 0 for x in d[4:]),
             "half-coordinate minimizer sign pattern changed")

        # Infinite lattice proof: integer coordinates have the unique minimum
        # zero; half coordinates have the unique minimum PSTAR.  Both obey the
        # required D8 parity condition and both total energies are zero.
        integer_unique = all(
            (F(n * n, 2) + dk * n > 0)
            for dk in d for n in (-3, -2, -1, 1, 2, 3)
        )
        # The displayed finite n test is accompanied by the sign formula in the
        # JSON/report: n(n/2+d)>0 for every nonzero integer n when |d|<1/2.
        half_minimizer = tuple(F(1, 2) if x < 0 else F(-1, 2) for x in d)
        need(integer_unique and half_minimizer == PSTAR, "ground minimizer classification failed")
        need(energy_k(ZERO, alpha) == 0 and energy_k(PSTAR, alpha) == 0,
             "ground energies are not zero")

        # E8 is even, so every nonzero non-root lattice vector has norm^2>=4.
        # Here ||d||^2=85/128<49/64, hence ||d||<7/8.  Cauchy gives, for
        # ||p||>=2, K(p)>=||p||^2/2-||d||||p||>2-7/4=1/4.
        # Thus the exact gap <=1/4 is decided by the 240 roots alone.
        delta_norm_squared = sum(x * x for x in d)
        need(delta_norm_squared == F(85, 128), "delta norm changed")
        need(delta_norm_squared < F(49, 64), "radial non-root bound failed")

        root_states = []
        for p in roots:
            q = matvec(matrix, p)
            root_states.append((energy_k(p, alpha), color_triality(q), "root", p, q))
        root_states.sort(key=lambda row: (row[0], row[1], row[3]))
        zero_states = [(F(0), 0, "vacuum", ZERO, matvec(matrix, ZERO))]
        zero_states += [row for row in root_states if row[0] == 0]
        neutral_states = [row for row in zero_states + root_states if row[1] == 0]
        neutral_zeros = [row for row in zero_states if row[1] == 0]
        positive_neutral = [row for row in root_states if row[1] == 0 and row[0] > 0]
        need(len(zero_states) == 2 and {row[3] for row in zero_states} == {ZERO, PSTAR},
             f"zero-state classification changed: {zero_states}")
        need(len(neutral_zeros) == 1 and neutral_zeros[0][3] == ZERO,
             "center-neutral zero space is not the vacuum line")
        need(positive_neutral, "no positive neutral state in certified window")
        gap_k = positive_neutral[0][0]
        need(gap_k <= F(1, 4), "root census does not contain the certified neutral gap")

        cases[fstr(alpha)] = {
            "delta": [fstr(x) for x in d],
            "all_H8_zero_charges": [
                {
                    "source_charge": [fstr(x) for x in row[3]],
                    "carrier_weight_Mp": [fstr(x) for x in row[4]],
                    "color_triality": row[1],
                    "P_sing_eigenvalue": psing_eigenvalue(row[1]),
                }
                for row in zero_states
            ],
            "full_ground_dimension": len(zero_states),
            "center_neutral_ground_dimension": len(neutral_zeros),
            "center_neutral_gap_K": fstr(gap_k),
            "center_neutral_gap_H8_equals_8K": fstr(8 * gap_k),
            "first_positive_center_neutral_charge": [fstr(x) for x in positive_neutral[0][3]],
            "first_positive_center_neutral_carrier_weight": [fstr(x) for x in positive_neutral[0][4]],
            "root_census_size": len(roots),
            "delta_norm_squared": fstr(delta_norm_squared),
            "nonroot_radial_bound": "For norm(p)>=2, K(p)>=norm(p)^2/2-||delta|| norm(p)>1/4 because ||delta||^2=85/128<49/64.",
            "oscillator_levels_in_gap_window": [0],
        }

    result = {
        "verdict": "PASS_CENTER_PROJECTOR_REMOVES_THE_ONLY_EXTRA_H8_ZERO_STATE",
        "sources": {
            "fold_results": str(FOLD_RESULTS),
            "fold_results_sha256": hashlib.sha256(FOLD_RESULTS.read_bytes()).hexdigest(),
            "old_ground_proof": str(REPO / "experiments/theory-contracts/quotient-locality-e8-20261014/checker.py"),
            "carrier_action_proof": str(REPO / "experiments/theory-contracts/seam-partner-original-action-20261015/checker.py"),
            "original_admissibility_source": str(REPO / "documents/archive--tfpt-45--source_extracts--04_qft_source/sections/002-section-closure-theorems-from-the-admissibility-complex.tex"),
        },
        "fixed_dictionary": {
            "p_star": [fstr(x) for x in PSTAR],
            "M_p_star": [fstr(x) for x in pstar_image],
            "carrier_name": "e1+e4",
            "carrier_representation": "(3,2)_(1/6)",
            "scaled_hypercharge_X": fstr(xcharge),
            "hypercharge_Y": fstr(xcharge / 6),
            "transported_source_U1_charge": fstr(source_u1),
            "color_triality": pstar_triality,
            "color_center_phase": "omega",
            "P_sing_eigenvalue": psing_eigenvalue(pstar_triality),
            "vacuum_P_sing_eigenvalue": 1,
        },
        "exact_projector": "P_sing=(1+Z_c+Z_c^2)/3; on triality t it is delta_(t mod 3,0)",
        "commutation": "H8 and Z_c are diagonal Cartan actions on lattice charge states, hence [H8,P_sing]=0 on the stipulated extension.",
        "ground_classification_proof": {
            "integer_coset": "For nonzero integer n, n^2/2+d*n=n(n/2+d)>0 when |d|<1/2; the unique minimizer is p=0.",
            "half_coset": "Coordinatewise strict convexity selects +1/2 where d<0 and -1/2 where d>0, hence uniquely p_star; its total energy is zero.",
            "all_zero_charges": ["0", "p_star"],
        },
        "cases": cases,
        "scope": {
            "proved": "On the stipulated E8 lattice extension with the fixed Weyl/carrier dictionary, the original color-center projector removes p_star and the center-neutral restriction has a one-dimensional H8 kernel and a positive exact gap.",
            "not_proved": [
                "that raw P1 alone already contains the later carrier projector",
                "that P_prim and P_Theta retain the vacuum line",
                "that H8 is the reconstructed physical OS transfer Hamiltonian",
                "clustering of the full graded field net or a 4D completion",
            ],
        },
        "checks_pass": True,
    }

    out = HERE / "clustering-admissible-sector-results.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
