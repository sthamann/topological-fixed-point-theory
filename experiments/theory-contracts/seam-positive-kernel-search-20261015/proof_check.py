#!/usr/bin/env python3
"""Offline exact checker for the two stored m=5 finite completions.

The checker never calls an optimizer.  It reads the Gaussian-rational matrices
from a completion, checks their algebraic and positivity relations, and replays
the numerical discovery moments only after independent rational recovery.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp

from exact_certificate import _psd_ldl
from word_algebra import T, WordAlgebra
from linear_relations import check_linear_relations


SCHEMA = "tfpt.positive-word-kernel.offline-proof-check.v1"
COMPLETION_SCHEMA = "tfpt.positive-word-kernel.finite-completion.v1"
MOMENT_SCHEMA = "tfpt.unknown-positive-word-kernel.v1"
MAX_DENOMINATOR = 64
TOLERANCE = 1e-6


class ProofCheckError(ValueError):
    """The stored data do not pass the exact offline proof checks."""


def _rational_pair(value: Any, label: str) -> sp.Rational:
    if not isinstance(value, list) or len(value) != 2:
        raise ProofCheckError(f"{label} must be an exact [numerator, denominator] pair")
    numerator, denominator = value
    if type(numerator) is not int or type(denominator) is not int or denominator <= 0:
        raise ProofCheckError(f"{label} is not a normalized rational pair")
    return sp.Rational(numerator, denominator)


def _matrix_scalar(value: Any, label: str) -> sp.Expr:
    if not isinstance(value, dict) or set(value) != {"re", "im"}:
        raise ProofCheckError(f"{label} must contain exactly re and im")
    return _rational_pair(value["re"], f"{label}.re") + sp.I * _rational_pair(
        value["im"], f"{label}.im"
    )


def _matrix(value: Any, rows: int, columns: int, label: str) -> sp.SparseMatrix:
    if not isinstance(value, list) or len(value) != rows:
        raise ProofCheckError(f"{label} must have {rows} rows")
    entries: dict[tuple[int, int], sp.Expr] = {}
    for row, source_row in enumerate(value):
        if not isinstance(source_row, list) or len(source_row) != columns:
            raise ProofCheckError(f"{label} row {row} must have {columns} columns")
        for column, source_value in enumerate(source_row):
            exact = sp.expand(_matrix_scalar(source_value, f"{label}[{row},{column}]"))
            if exact != 0:
                entries[row, column] = exact
    return sp.SparseMatrix(rows, columns, entries)


def _declared_scalar(value: Any, label: str) -> sp.Expr:
    if not isinstance(value, dict) or set(value) != {"re", "im"}:
        raise ProofCheckError(f"{label} must contain exactly re and im")
    try:
        real = sp.Rational(value["re"])
        imaginary = sp.Rational(value["im"])
    except (TypeError, ValueError) as error:
        raise ProofCheckError(f"{label} is not a Gaussian rational") from error
    return sp.expand(real + sp.I * imaginary)


def _recover_float(value: Any, label: str) -> tuple[sp.Rational, float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProofCheckError(f"{label} must be a finite JSON number")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ProofCheckError(f"{label} must be finite")
    candidate = Fraction(numeric).limit_denominator(MAX_DENOMINATOR)
    error = abs(numeric - float(candidate))
    if error > TOLERANCE:
        raise ProofCheckError(
            f"{label} has no rational denominator <= {MAX_DENOMINATOR} "
            f"within {TOLERANCE}: error={error}"
        )
    return sp.Rational(candidate.numerator, candidate.denominator), error


def _canonical_json_hash(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


class _ExactRepresentation:
    def __init__(self, omega: sp.SparseMatrix, transfer: sp.SparseMatrix,
                 gammas: list[sp.SparseMatrix]):
        self.omega = omega
        self.operators = {T: transfer, **dict(enumerate(gammas))}
        self.columns: dict[int, tuple[tuple[tuple[int, sp.Expr], ...], ...]] = {}
        for token, operator in self.operators.items():
            by_column: list[list[tuple[int, sp.Expr]]] = [
                [] for _ in range(operator.cols)
            ]
            for (row, column), exact in operator.todok().items():
                by_column[column].append((row, exact))
            self.columns[token] = tuple(tuple(column) for column in by_column)
        self.omega_sparse = {
            row: value for (row, column), value in omega.todok().items()
            if column == 0 and value != 0
        }

    @lru_cache(None)
    def vector(self, word: tuple[int, ...]) -> tuple[tuple[int, sp.Expr], ...]:
        if not word:
            return tuple(sorted(self.omega_sparse.items()))
        token = word[0]
        if token not in self.columns:
            raise ProofCheckError(f"unknown word token {token}")
        result: dict[int, sp.Expr] = {}
        for column, coefficient in self.vector(word[1:]):
            for row, factor in self.columns[token][column]:
                result[row] = result.get(row, 0) + factor * coefficient
        return tuple(sorted(
            (row, sp.expand(value)) for row, value in result.items()
            if sp.expand(value) != 0
        ))

    def moment(self, word: tuple[int, ...]) -> sp.Expr:
        vector = dict(self.vector(word))
        return sp.expand(sum(
            sp.conjugate(coefficient) * vector.get(row, 0)
            for row, coefficient in self.omega_sparse.items()
        ))


def _load_moments(payload: dict[str, Any], algebra: WordAlgebra) -> tuple[dict[tuple[int, ...], sp.Expr], float]:
    if payload.get("schema") != MOMENT_SCHEMA or payload.get("m") != algebra.m:
        raise ProofCheckError("moment source has the wrong schema or mode count")
    rows = payload.get("moments")
    if not isinstance(rows, list) or len(rows) != 2245:
        raise ProofCheckError("the m=5 source must contain exactly 2245 moments")
    result: dict[tuple[int, ...], sp.Expr] = {}
    max_error = 0.0
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or not {"word", "re", "im"}.issubset(row):
            raise ProofCheckError(f"moment {index} lacks word/re/im")
        if not isinstance(row["word"], list):
            raise ProofCheckError(f"moment {index} word is not a list")
        raw_word = tuple(row["word"])
        try:
            phase, word = algebra.reduce_word(raw_word)
        except ValueError as error:
            raise ProofCheckError(f"invalid moment word {raw_word}: {error}") from error
        if phase != 1 or word != raw_word:
            raise ProofCheckError(f"moment word {raw_word} is not canonical")
        if word in result:
            raise ProofCheckError(f"duplicate moment word {word}")
        real, real_error = _recover_float(row["re"], f"moment {index}.re")
        imaginary, imaginary_error = _recover_float(row["im"], f"moment {index}.im")
        result[word] = sp.expand(real + sp.I * imaginary)
        max_error = max(max_error, real_error, imaginary_error)
    expected_words = set(algebra.words(4))
    if set(result) != expected_words:
        missing = len(expected_words - set(result))
        extra = len(set(result) - expected_words)
        raise ProofCheckError(
            f"moment source is not the complete degree<=4 basis: missing={missing}, extra={extra}"
        )
    return result, max_error


def _check_certificate_impl(certificate: dict[str, Any], moments_payload: dict[str, Any]) -> dict[str, Any]:
    if certificate.get("schema") != COMPLETION_SCHEMA:
        raise ProofCheckError("wrong finite-completion schema")
    if certificate.get("m") != 5 or certificate.get("dimension") != 32:
        raise ProofCheckError("this checker requires the stored m=5, dimension-32 completion")
    if not certificate.get("ok") or not certificate.get("complete_finite_kernel_certificate"):
        raise ProofCheckError("input does not declare a completed finite certificate")
    expected_hash = _canonical_json_hash(moments_payload)
    if certificate.get("input_sha256") != expected_hash:
        raise ProofCheckError("completion input hash does not match the supplied moment source")

    algebra = WordAlgebra(5)
    moments, max_recovery_error = _load_moments(moments_payload, algebra)
    exact = certificate.get("exact_matrices")
    if not isinstance(exact, dict):
        raise ProofCheckError("missing exact_matrices")
    transfer = _matrix(exact.get("T"), 32, 32, "T")
    omega = _matrix(exact.get("Omega"), 32, 1, "Omega")
    gamma_json = exact.get("gammas")
    if not isinstance(gamma_json, list) or len(gamma_json) != 10:
        raise ProofCheckError("exact_matrices.gammas must contain ten matrices")
    gammas = [_matrix(value, 32, 32, f"gamma_{index}")
              for index, value in enumerate(gamma_json)]
    quarter = _matrix(exact.get("quarter"), 32, 32, "quarter")
    identity = sp.eye(32)
    zero = sp.zeros(32)

    if transfer != transfer.adjoint():
        raise ProofCheckError("T is not exactly selfadjoint")
    if (omega.adjoint() * omega)[0] != 1:
        raise ProofCheckError("Omega is not exactly normalized")
    if transfer * omega != omega:
        raise ProofCheckError("T Omega != Omega")
    for index, gamma in enumerate(gammas):
        if gamma != gamma.adjoint():
            raise ProofCheckError(f"gamma_{index} is not selfadjoint")
        for previous in range(index + 1):
            expected = 2 * identity if previous == index else zero
            if gamma * gammas[previous] + gammas[previous] * gamma != expected:
                raise ProofCheckError(f"CAR failed for gamma_{index}, gamma_{previous}")

    if quarter.adjoint() * quarter != identity:
        raise ProofCheckError("R is not exactly unitary")
    if quarter ** 4 != identity:
        raise ProofCheckError("R^4 != I")
    if quarter * omega != omega:
        raise ProofCheckError("R Omega != Omega")
    if quarter * transfer != transfer * quarter:
        raise ProofCheckError("R T != T R")
    for index, gamma in enumerate(gammas):
        expected = gammas[index + 1] if index % 2 == 0 else -gammas[index - 1]
        if quarter * gamma * quarter.adjoint() != expected:
            raise ProofCheckError(f"quarter covariance failed for gamma_{index}")

    floor, floor_error = _recover_float(
        moments_payload.get("transfer_floor"), "transfer_floor"
    )
    if not 0 <= floor <= 1:
        raise ProofCheckError("transfer_floor is outside [0,1]")
    transfer_psd = _psd_ldl(transfer, "T")
    floor_psd = _psd_ldl(transfer - floor * identity, "T-transfer_floor*I")
    complement_psd = _psd_ldl(identity - transfer, "I-T")
    finite_clustering = complement_psd["rank"] == 31
    stationary_dimension = 32 - complement_psd["rank"]
    declared_clustering = certificate.get("finite_clustering")
    if declared_clustering is not None:
        if type(declared_clustering) is not bool:
            raise ProofCheckError("finite_clustering must be a boolean when present")
        if declared_clustering is not finite_clustering:
            raise ProofCheckError(
                "declared finite_clustering disagrees with exact rank(I-T): "
                f"declared={declared_clustering}, computed={finite_clustering}"
            )
    cluster_gap, cluster_gap_error = _recover_float(
        moments_payload.get("cluster_gap", 0), "cluster_gap"
    )
    if not 0 <= cluster_gap <= 1:
        raise ProofCheckError("cluster_gap is outside [0,1]")
    cluster_bound_psd = None
    if cluster_gap > 0:
        omega_projection = omega * omega.adjoint()
        cluster_bound = ((1 - cluster_gap) * identity
                         + cluster_gap * omega_projection - transfer)
        cluster_bound_psd = _psd_ldl(cluster_bound, "cluster upper bound minus T")
        if not finite_clustering:
            raise ProofCheckError(
                "positive cluster_gap was declared but rank(I-T) does not prove finite clustering"
            )

    representation = _ExactRepresentation(omega, transfer, gammas)
    failures = []
    for word, expected in moments.items():
        actual = representation.moment(word)
        if actual != expected:
            failures.append({"word": list(word), "expected": str(expected), "actual": str(actual)})
            if len(failures) == 3:
                break
    if failures:
        raise ProofCheckError(f"supplied moment replay failed: {failures}")
    try:
        check_linear_relations(moments_payload.get("linear_relations", []), moments, 5, moments_payload["degree"])
    except ValueError as error:
        raise ProofCheckError(str(error)) from error

    predictions = certificate.get("candidate_predictions")
    if not isinstance(predictions, list) or not predictions:
        raise ProofCheckError("missing declared degree-5 discriminator")
    discriminator: dict[tuple[int, ...], sp.Expr] = {}
    for index, row in enumerate(predictions):
        if not isinstance(row, dict) or set(row) != {"word", "value"}:
            raise ProofCheckError(f"discriminator {index} has the wrong shape")
        word = tuple(row["word"])
        try:
            algebra.reduce_word(word)
        except ValueError as error:
            raise ProofCheckError(f"invalid discriminator word {word}: {error}") from error
        if len(word) != 5 or word.count(T) != 1:
            raise ProofCheckError(f"discriminator word {word} is not degree 5 with one T")
        if word in moments or word in discriminator:
            raise ProofCheckError(f"discriminator word {word} is supplied or duplicated")
        declared = _declared_scalar(row["value"], f"discriminator {index}.value")
        actual = representation.moment(word)
        if actual != declared:
            raise ProofCheckError(
                f"degree-5 discriminator replay failed at {word}: {actual} != {declared}"
            )
        discriminator[word] = actual

    return {
        "schema": SCHEMA,
        "ok": True,
        "verdict": "EXACT_OFFLINE_FINITE_COMPLETION_REPLAY",
        "objective": certificate.get("objective"),
        "dimension": 32,
        "moment_count": len(moments),
        "max_moment_recovery_error": max(
            max_recovery_error, floor_error, cluster_gap_error
        ),
        "transfer_floor": str(floor),
        "cluster_gap": str(cluster_gap),
        "finite_clustering": finite_clustering,
        "stationary_dimension": stationary_dimension,
        "degree5_discriminator": {
            "word_count": len(discriminator),
            "values": {json.dumps(list(word)): str(value)
                       for word, value in discriminator.items()},
        },
        "positivity": {
            "T": transfer_psd,
            "T_minus_floor": floor_psd,
            "I_minus_T": complement_psd,
            "cluster_upper_bound_minus_T": cluster_bound_psd,
        },
        "checks": [
            "input_hash",
            "exact_matrix_shapes_and_scalars",
            "CAR",
            "Omega_norm_and_T_stationarity",
            "quarter_unitarity_order_vacuum_and_covariance",
            "exact_transfer_interval",
            "finite_clustering_by_rank_I_minus_T",
            "cluster_gap_upper_bound_when_declared",
            "all_2245_source_moments",
            "declared_degree5_discriminator",
        ],
        "scope": {
            "proved": "exact replay of this explicit 32-dimensional completion",
            "not_proved": [
                "flatness or uniqueness of the truncated moment problem",
                "selection of either completion by P1/P2 or the source theory",
                "a Euclidean UCP time channel or event instrument",
                "infinite-volume or OS clustering beyond this finite analogue",
            ],
        },
    }


def check_certificate(certificate: dict[str, Any], moments_payload: dict[str, Any]) -> dict[str, Any]:
    """Check one completion and return a compact success or honest failure."""
    try:
        return _check_certificate_impl(certificate, moments_payload)
    except (ProofCheckError, ValueError, TypeError, KeyError, ZeroDivisionError) as error:
        return {
            "schema": SCHEMA,
            "ok": False,
            "verdict": "OFFLINE_EXACT_CHECK_FAILED",
            "error": str(error),
        }


def check_pair(minimum: dict[str, Any], maximum: dict[str, Any],
               moments_payload: dict[str, Any],
               next_source_question: dict[str, Any] | None = None) -> dict[str, Any]:
    """Check both stored completions and their declared degree-5 separation."""
    minimum_result = check_certificate(minimum, moments_payload)
    maximum_result = check_certificate(maximum, moments_payload)
    if not minimum_result["ok"] or not maximum_result["ok"]:
        return {
            "schema": SCHEMA,
            "ok": False,
            "verdict": "OFFLINE_EXACT_CHECK_FAILED",
            "minimum": minimum_result,
            "maximum": maximum_result,
        }
    if minimum_result["objective"] != "min-trace" or maximum_result["objective"] != "max-trace":
        return {"schema": SCHEMA, "ok": False, "verdict": "OFFLINE_EXACT_CHECK_FAILED",
                "error": "completion objectives are not min-trace/max-trace"}
    low = minimum_result["degree5_discriminator"]["values"]
    high = maximum_result["degree5_discriminator"]["values"]
    if set(low) != set(high) or not low or any(low[word] == high[word] for word in low):
        return {"schema": SCHEMA, "ok": False, "verdict": "OFFLINE_EXACT_CHECK_FAILED",
                "error": "declared degree-5 words do not exactly separate both completions"}
    if next_source_question is not None:
        try:
            if (next_source_question.get("status") != "EXACT_DISTINGUISHING_WORDS_FOUND"
                    or next_source_question.get("same_input_moments_replayed") != 2245
                    or next_source_question.get("all_words_through_degree_replayed") != 4
                    or next_source_question.get("shortest_found_degree") != 5
                    or next_source_question.get("shortest_certified") is not True):
                raise ProofCheckError("next-source summary has inconsistent degree/replay claims")
            declared_differences = next_source_question.get("differences")
            if not isinstance(declared_differences, list):
                raise ProofCheckError("next-source differences must be a list")
            reconstructed = {}
            for index, row in enumerate(declared_differences):
                if not isinstance(row, dict) or set(row) != {"word", "first", "second"}:
                    raise ProofCheckError(f"next-source difference {index} has the wrong shape")
                key = json.dumps(row["word"])
                reconstructed[key] = (
                    str(_declared_scalar(row["first"], f"difference {index}.first")),
                    str(_declared_scalar(row["second"], f"difference {index}.second")),
                )
            expected = {word: (low[word], high[word]) for word in low}
            if reconstructed != expected:
                raise ProofCheckError("next-source differences do not match exact matrix replay")
        except (ProofCheckError, TypeError, KeyError, ValueError) as error:
            return {"schema": SCHEMA, "ok": False, "verdict": "OFFLINE_EXACT_CHECK_FAILED",
                    "error": str(error)}
    return {
        "schema": SCHEMA,
        "ok": True,
        "verdict": "TWO_EXACT_FINITE_COMPLETIONS_WITH_DEGREE5_SEPARATOR",
        "minimum": minimum_result,
        "maximum": maximum_result,
        "pair_check": {
            "shared_source_moments": 2245,
            "all_shared_words_through_degree": 4,
            "separating_degree5_words": len(low),
            "all_declared_degree5_values_differ": True,
            "next_source_question_replayed": next_source_question is not None,
            "minimum_finite_clustering": minimum_result["finite_clustering"],
            "maximum_finite_clustering": maximum_result["finite_clustering"],
        },
        "scope": minimum_result["scope"],
    }


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", nargs="?", type=Path,
                        default=root / "runs/certified-positive")
    parser.add_argument("--out", "--output", dest="output", type=Path)
    args = parser.parse_args()
    run_directory = args.run_directory.resolve()
    minimum_path = run_directory / "min-trace.json"
    maximum_path = run_directory / "max-trace.json"
    moments_path = run_directory / "moments.json"
    question_path = run_directory / "next-source-question.json"
    minimum = json.loads(minimum_path.read_text())
    maximum = json.loads(maximum_path.read_text())
    moments_payload = json.loads(moments_path.read_text())
    next_source_question = json.loads(question_path.read_text())
    result = check_pair(minimum, maximum, moments_payload, next_source_question)
    result["provenance"] = {
        "run_directory": str(run_directory),
        "min_sha256": sha256(minimum_path.read_bytes()).hexdigest(),
        "max_sha256": sha256(maximum_path.read_bytes()).hexdigest(),
        "moments_sha256": sha256(moments_path.read_bytes()).hexdigest(),
        "next_source_question_sha256": sha256(question_path.read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    # Keep stdout compact even though the source certificates are about 2 MB.
    print(json.dumps({
        "ok": result["ok"],
        "verdict": result["verdict"],
        "min_moments": result.get("minimum", {}).get("moment_count"),
        "max_moments": result.get("maximum", {}).get("moment_count"),
        "min_finite_clustering": result.get("minimum", {}).get("finite_clustering"),
        "max_finite_clustering": result.get("maximum", {}).get("finite_clustering"),
        "degree5_separators": result.get("pair_check", {}).get("separating_degree5_words"),
    }, sort_keys=True))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
