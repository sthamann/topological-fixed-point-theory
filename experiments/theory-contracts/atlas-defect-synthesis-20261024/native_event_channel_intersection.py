#!/usr/bin/env python3
"""Exact intersection of the six-effect channels with native random unitaries.

The native H5 action factors through the actual S6 quotient of G31.  A fully
covariant random-unitary channel can therefore be conjugation averaged.  This
checker enumerates the eleven S6 classes, reconstructs their actual channel
matrices from native reflection words, and solves the remaining rational
three-equality polytope exactly.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path
from typing import Any

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DOILY_SOURCE = HERE / "doily_quartic_bridge.py"
INSTRUMENT_SOURCE = HERE / "covariant_instruments.py"
DEFAULT_OUTPUT = HERE / "native_event_channel_intersection_results.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def clean(matrix: sp.MatrixBase) -> sp.Matrix:
    return sp.Matrix(matrix).applyfunc(sp.factor)


def expr_json(value: sp.Expr) -> int | str:
    value = sp.factor(value)
    return int(value) if value.is_Integer else str(value)


def matrix_json(matrix: sp.MatrixBase) -> list[list[int | str]]:
    return [
        [expr_json(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def matrix_hash(matrix: sp.MatrixBase) -> str:
    payload = json.dumps(matrix_json(matrix), separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_doily() -> Any:
    spec = importlib.util.spec_from_file_location("doily_event_intersection", DOILY_SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load doily bridge")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(6))


def cycle_type(permutation: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(6):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def partitions(number: int, maximum: int | None = None):
    if number == 0:
        yield ()
        return
    maximum = min(number, maximum if maximum is not None else number)
    for first in range(maximum, 0, -1):
        for remainder in partitions(number - first, first):
            yield (first,) + remainder


def characters(cycles: tuple[int, ...]) -> tuple[sp.Integer, sp.Integer, sp.Integer]:
    multiplicities = Counter(cycles)
    m1 = multiplicities[1]
    m2 = multiplicities[2]
    chi5 = sp.Integer(m1 - 1)
    chi_square = sp.Integer(m1 + 2 * m2 - 1)
    chi9 = (chi5**2 + chi_square) / 2 - chi5 - 1
    chi10 = (chi5**2 - chi_square) / 2
    require(all(value.is_Integer for value in (chi5, chi9, chi10)), "character left Z")
    return chi5, sp.Integer(chi9), sp.Integer(chi10)


def column_vectorize(matrix: sp.MatrixBase) -> sp.Matrix:
    return sp.Matrix(
        [
            matrix[row, column]
            for column in range(matrix.cols)
            for row in range(matrix.rows)
        ]
    )


def conjugation_superoperator(matrix: sp.Matrix) -> sp.Matrix:
    return sp.kronecker_product(matrix.inv().T, matrix)


def group_closure(
    generators: list[tuple[int, tuple[int, ...], sp.Matrix]],
) -> dict[tuple[int, ...], tuple[sp.Matrix, list[int]]]:
    identity = tuple(range(6))
    group: dict[tuple[int, ...], tuple[sp.Matrix, list[int]]] = {
        identity: (sp.eye(5), [])
    }
    queue: deque[tuple[int, ...]] = deque([identity])
    while queue:
        permutation = queue.popleft()
        matrix, word = group[permutation]
        for reflection_index, generator_permutation, generator_matrix in generators:
            next_permutation = compose(permutation, generator_permutation)
            next_matrix = clean(matrix * generator_matrix)
            if next_permutation in group:
                require(group[next_permutation][0] == next_matrix, "S6 action is inconsistent")
                continue
            group[next_permutation] = (next_matrix, word + [reflection_index])
            queue.append(next_permutation)
    return group


def instrument_channel(
    parameter: sp.Rational,
    effects: list[sp.Matrix],
    projectors: list[sp.Matrix],
) -> sp.Matrix:
    identity = sp.eye(5)
    channel = sp.zeros(25)
    for effect, projector in zip(effects, projectors):
        output = clean(parameter * projector + (1 - parameter) * (identity - projector) / 4)
        channel += column_vectorize(output) * column_vectorize(effect.T).T
    return clean(channel)


def exact_lp_vertices(class_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    vertices: dict[tuple[sp.Expr, ...], dict[str, Any]] = {}
    for support_size in range(1, 4):
        for support in it.combinations(range(len(class_rows)), support_size):
            constraint = sp.Matrix(
                [
                    [1 for _ in support],
                    [class_rows[index]["chi9"] for index in support],
                    [class_rows[index]["chi10"] for index in support],
                ]
            )
            # A basic feasible point must have independent active columns.
            if constraint.rank() != support_size:
                continue
            solutions = list(sp.linsolve((constraint, sp.Matrix([1, 0, 0]))))
            if len(solutions) != 1:
                continue
            solution = tuple(sp.factor(value) for value in solutions[0])
            if any(value <= 0 for value in solution):
                continue
            full = tuple(
                solution[support.index(index)] if index in support else sp.Integer(0)
                for index in range(len(class_rows))
            )
            average_chi5 = sp.factor(
                sum(full[index] * class_rows[index]["chi5"] for index in range(len(class_rows)))
            )
            lambda5 = sp.factor(average_chi5 / 5)
            parameter = sp.factor(lambda5 + sp.Rational(1, 5))
            vertices[full] = {
                "p": parameter,
                "lambda5": lambda5,
                "average_chi5": average_chi5,
                "weights": [
                    {
                        "cycle_type": list(class_rows[index]["cycle_type"]),
                        "total_class_weight": solution[position],
                        "class_size": class_rows[index]["size"],
                        "per_element_weight": sp.factor(
                            solution[position] / class_rows[index]["size"]
                        ),
                    }
                    for position, index in enumerate(support)
                ],
                "full_weight_vector": full,
            }
    return sorted(
        vertices.values(),
        key=lambda vertex: (
            vertex["p"],
            len(vertex["weights"]),
            tuple(tuple(weight["cycle_type"]) for weight in vertex["weights"]),
        ),
    )


def serialise_vertex(vertex: dict[str, Any]) -> dict[str, Any]:
    return {
        "p": expr_json(vertex["p"]),
        "lambda5": expr_json(vertex["lambda5"]),
        "average_chi5": expr_json(vertex["average_chi5"]),
        "weights": [
            {
                "cycle_type": weight["cycle_type"],
                "total_class_weight": expr_json(weight["total_class_weight"]),
                "class_size": weight["class_size"],
                "per_element_weight": expr_json(weight["per_element_weight"]),
            }
            for weight in vertex["weights"]
        ],
        "full_weight_vector": [expr_json(value) for value in vertex["full_weight_vector"]],
    }


def endpoint_word_certificate(
    weights: dict[tuple[int, ...], sp.Rational],
    classes: dict[tuple[int, ...], list[tuple[tuple[int, ...], sp.Matrix, list[int]]]],
) -> list[dict[str, Any]]:
    certificate = []
    for cycles, total_weight in weights.items():
        members = sorted(classes[cycles], key=lambda item: item[0])
        certificate.append(
            {
                "cycle_type": list(cycles),
                "class_size": len(members),
                "total_class_weight": expr_json(total_weight),
                "per_element_weight": expr_json(total_weight / len(members)),
                "actual_native_words": [
                    {"permutation": list(permutation), "reflection_indices": word}
                    for permutation, _, word in members
                ],
            }
        )
    return certificate


def main(output: Path) -> dict[str, Any]:
    doily = load_doily()
    bridge = doily.build_report()
    rows = bridge["quartic_bridge"]["per_reflection"]
    require(len(rows) == 60, "native reflection count changed")

    unique_generators: dict[tuple[int, ...], tuple[int, sp.Matrix]] = {}
    for row in rows:
        permutation = tuple(row["ovoid_permutation"])
        matrix = sp.Matrix(row["quartic_numerator"]) / 16
        if permutation in unique_generators:
            require(unique_generators[permutation][1] == matrix, "duplicate action mismatch")
        else:
            unique_generators[permutation] = (row["index"], matrix)
    require(len(unique_generators) == 15, "expected fifteen transposition actions")
    generators = [
        (index, permutation, matrix)
        for permutation, (index, matrix) in sorted(unique_generators.items())
    ]
    group = group_closure(generators)
    require(len(group) == 720, "native words did not generate S6")

    cycle_types = list(partitions(6))
    require(len(cycle_types) == 11, "S6 class count changed")
    classes: dict[tuple[int, ...], list[tuple[tuple[int, ...], sp.Matrix, list[int]]]] = {
        cycles: [] for cycles in cycle_types
    }
    for permutation, (matrix, word) in group.items():
        classes[cycle_type(permutation)].append((permutation, matrix, word))
    require(sum(map(len, classes.values())) == 720, "class partition lost group elements")

    class_rows: list[dict[str, Any]] = []
    class_averages: dict[tuple[int, ...], sp.Matrix] = {}
    for cycles in cycle_types:
        chi5, chi9, chi10 = characters(cycles)
        members = classes[cycles]
        average = clean(
            sum((conjugation_superoperator(matrix) for _, matrix, _ in members), sp.zeros(25))
            / len(members)
        )
        class_averages[cycles] = average
        class_rows.append(
            {
                "cycle_type": cycles,
                "size": len(members),
                "m1": cycles.count(1),
                "m2": cycles.count(2),
                "chi5": chi5,
                "chi9": chi9,
                "chi10": chi10,
                "lambda5": sp.factor(chi5 / 5),
                "lambda9": sp.factor(chi9 / 9),
                "lambda10": sp.factor(chi10 / 10),
                "actual_class_channel_sha256": matrix_hash(average),
            }
        )

    # Reconstruct the four isotypic projectors on End(H5) from the actual
    # class channel matrices and the exact S6 characters.
    dimensions = {"1": 1, "5": 5, "9": 9, "10": 10}
    character_key = {"1": None, "5": "chi5", "9": "chi9", "10": "chi10"}
    isotypic: dict[str, sp.Matrix] = {}
    for label, dimension in dimensions.items():
        projector = sp.zeros(25)
        for class_row in class_rows:
            character = 1 if label == "1" else class_row[character_key[label]]
            projector += class_row["size"] * character * class_averages[class_row["cycle_type"]]
        isotypic[label] = clean(sp.Rational(dimension, 720) * projector)
    require(
        all(isotypic[label].rank() == dimensions[label] for label in dimensions),
        "End(H5) isotypic ranks changed",
    )
    require(
        all(clean(projector**2 - projector) == sp.zeros(25) for projector in isotypic.values()),
        "central character operator is not a projector",
    )
    require(
        all(
            clean(isotypic[left] * isotypic[right]) == sp.zeros(25)
            for left in dimensions
            for right in dimensions
            if left != right
        ),
        "isotypic projectors overlap",
    )
    require(clean(sum(isotypic.values(), sp.zeros(25))) == sp.eye(25), "isotypic sum is not I25")
    for class_row in class_rows:
        average = class_averages[class_row["cycle_type"]]
        for label, dimension in dimensions.items():
            character = 1 if label == "1" else class_row[character_key[label]]
            require(
                clean(average * isotypic[label] - sp.Rational(character, dimension) * isotypic[label])
                == sp.zeros(25),
                "actual class channel disagrees with its character ratio",
            )

    vertices = exact_lp_vertices(class_rows)
    require(len(vertices) == 29, "LP vertex census changed")
    minimum_p = min(vertex["p"] for vertex in vertices)
    maximum_p = max(vertex["p"] for vertex in vertices)
    require((minimum_p, maximum_p) == (sp.Rational(1, 10), sp.Rational(1, 2)), "p interval changed")

    # Analytic dual certificate for the complete interval.  Put chi=chi5.
    # The two zero-eigenvalue constraints imply E[chi^2]=1+E[chi].
    # Since chi is integral, chi(chi+1)>=0 and (chi-1)(chi-2)>=0 pointwise.
    # Averaging and substituting the moment identity yields
    # -1/2 <= E[chi] <= 3/2, hence 1/10 <= p=(1+E[chi])/5 <= 1/2.
    possible_chi5 = sorted({row["chi5"] for row in class_rows})
    require(
        all(value * (value + 1) >= 0 for value in possible_chi5),
        "lower dual polynomial lost positivity",
    )
    require(
        all((value - 1) * (value - 2) >= 0 for value in possible_chi5),
        "upper dual polynomial lost positivity",
    )

    lower_weights = {(6,): sp.Rational(1, 2), (3, 2, 1): sp.Rational(1, 2)}
    upper_weights = {
        (4, 1, 1): sp.Rational(1, 4),
        (3, 1, 1, 1): sp.Rational(1, 2),
        (2, 2, 1, 1): sp.Rational(1, 4),
    }

    metric = sp.diag(1, 3, 3, 3, 6)
    rays = [sp.Matrix(ray) for ray in bridge["selection_data"]["native_ovoid_rays"]]
    effects = [clean(ray * ray.T * metric / 12) for ray in rays]
    projectors = [clean(sp.Rational(6, 5) * effect) for effect in effects]

    def class_mixture(weights: dict[tuple[int, ...], sp.Rational]) -> sp.Matrix:
        return clean(sum((weight * class_averages[cycles] for cycles, weight in weights.items()), sp.zeros(25)))

    endpoint_channels = {
        "p=1/10": class_mixture(lower_weights),
        "p=1/2": class_mixture(upper_weights),
    }
    require(endpoint_channels["p=1/10"] == instrument_channel(sp.Rational(1, 10), effects, projectors), "lower endpoint is not Phi_p")
    require(endpoint_channels["p=1/2"] == instrument_channel(sp.Rational(1, 2), effects, projectors), "upper endpoint is not Phi_p")

    # Every enumerated vertex is also checked against the actual operator, not
    # merely against its four eigenvalue labels.
    for vertex in vertices:
        weights = {
            tuple(weight["cycle_type"]): weight["total_class_weight"]
            for weight in vertex["weights"]
        }
        require(
            class_mixture(weights) == instrument_channel(vertex["p"], effects, projectors),
            "LP vertex character solution failed as an actual H5 channel",
        )

    # Convex interpolation between the two explicit endpoint distributions
    # covers the complete interval.  The interpolation coefficient is
    # t=(10p-1)/4.
    parameter = sp.Symbol("p", real=True)
    interpolation = sp.factor((10 * parameter - 1) / 4)
    require(interpolation.subs(parameter, minimum_p) == 0, "lower interpolation changed")
    require(interpolation.subs(parameter, maximum_p) == 1, "upper interpolation changed")

    # The Lueders member p=1 lies strictly outside the random-unitary section.
    lueders = instrument_channel(sp.Integer(1), effects, projectors)
    require(lueders.rank() == 6, "Lueders Liouville rank changed")
    require(maximum_p < 1, "Lueders channel unexpectedly entered feasible interval")

    p_histogram = Counter(vertex["p"] for vertex in vertices)
    report = {
        "claim": "intersection of native central random-unitary channels with the six-effect Phi_p family",
        "status": "PASS",
        "arithmetic": "exact integers and rationals; actual native 5 by 5 word matrices",
        "sources": {
            str(DOILY_SOURCE.relative_to(REPO)): file_hash(DOILY_SOURCE),
            str(INSTRUMENT_SOURCE.relative_to(REPO)): file_hash(INSTRUMENT_SOURCE),
        },
        "actual_native_group": {
            "quotient": "S6 on the six ovoids",
            "order": len(group),
            "generator_count": len(generators),
            "generators": [
                {"reflection_index": index, "permutation": list(permutation)}
                for index, permutation, _ in generators
            ],
            "every_element_has_native_reflection_word": True,
        },
        "end_H5_decomposition": {
            "dimensions": [1, 5, 9, 10],
            "sum": 25,
            "projector_ranks": {label: isotypic[label].rank() for label in dimensions},
            "projector_sha256": {label: matrix_hash(isotypic[label]) for label in dimensions},
            "all_projectors_exact_idempotents": True,
            "all_class_eigenvalues_match_character_ratios": True,
        },
        "class_table": [
            {
                key: (list(value) if key == "cycle_type" else expr_json(value) if isinstance(value, sp.Expr) else value)
                for key, value in class_row.items()
            }
            for class_row in class_rows
        ],
        "linear_program": {
            "constraints": [
                "sum_C w_C = 1",
                "sum_C w_C chi9(C) = 0",
                "sum_C w_C chi10(C) = 0",
                "w_C >= 0",
            ],
            "objective_relation": "lambda5=sum_C w_C chi5(C)/5=p-1/5",
            "class_count": len(class_rows),
            "constraint_rank": sp.Matrix(
                [
                    [1] * len(class_rows),
                    [row["chi9"] for row in class_rows],
                    [row["chi10"] for row in class_rows],
                ]
            ).rank(),
            "vertex_count": len(vertices),
            "vertices": [serialise_vertex(vertex) for vertex in vertices],
            "vertex_p_histogram": {str(value): count for value, count in sorted(p_histogram.items())},
            "reachable_lambda5_interval": ["-1/10", "3/10"],
            "reachable_p_interval": ["1/10", "1/2"],
            "all_vertex_channels_checked_as_actual_25_by_25_operators": True,
        },
        "analytic_dual_certificate": {
            "chi": "chi5=m1-1 is integral on every S6 class",
            "possible_chi_values": [int(value) for value in possible_chi5],
            "moment_identity_from_lambda9_equals_lambda10_equals_zero": "E[chi^2]=1+E[chi]",
            "lower_pointwise_polynomial": "chi(chi+1)>=0",
            "lower_consequence": "E[chi]>=-1/2",
            "upper_pointwise_polynomial": "(chi-1)(chi-2)>=0",
            "upper_consequence": "E[chi]<=3/2",
            "parameter_relation": "p=(1+E[chi])/5",
            "certified_interval": ["1/10", "1/2"],
        },
        "explicit_endpoint_mixtures": {
            "p=1/10": {
                "weights": {"6": "1/2", "3,2,1": "1/2"},
                "actual_word_certificate": endpoint_word_certificate(lower_weights, classes),
                "channel_sha256": matrix_hash(endpoint_channels["p=1/10"]),
            },
            "p=1/2": {
                "weights": {"4,1,1": "1/4", "3,1,1,1": "1/2", "2,2,1,1": "1/4"},
                "actual_word_certificate": endpoint_word_certificate(upper_weights, classes),
                "channel_sha256": matrix_hash(endpoint_channels["p=1/2"]),
            },
            "all_intermediate_p": "mix lower and upper endpoints with t=(10p-1)/4",
        },
        "lueders_test": {
            "p": 1,
            "lambda5": "4/5",
            "liouville_rank": lueders.rank(),
            "native_random_unitary_feasible": False,
            "reason": "maximum feasible lambda5 is 3/10, equivalently maximum p is 1/2",
        },
        "interpretation": {
            "positive": "Phi_p is an exact mixture of actual native event conjugations exactly for 1/10 <= p <= 1/2",
            "negative": "the p=1 Lueders channel and all p outside that interval are not such mixtures",
            "physical_weight_selection_derived": False,
        },
    }
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = main(args.output)
    print(
        "PASS: 11 S6 classes, 29 exact LP vertices, "
        "native random-unitary intersection p in [1/10,1/2]; p=1 excluded."
    )
