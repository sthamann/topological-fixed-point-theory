#!/usr/bin/env python3
"""Exact recovery and finite-GNS certificate for positive word kernels.

The SDP is discovery only.  This module reconstructs every supplied complex
moment as a Gaussian rational, rebuilds the moment/localizing data using the
shared CAR word algebra, extracts a flat finite representation, and replays
all supplied moments exactly.  It proves neither a UCP time channel nor the
physical selection of the kernel.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp

from word_algebra import T, WordAlgebra, SCHEMA_ID


CERTIFICATE_SCHEMA = "tfpt.positive_word_kernel.exact_certificate.v1"


class CertificateError(ValueError):
    """A numerical candidate did not recover to an exact valid certificate."""


def _recover_rational(value: Any, max_denominator: int, tolerance: float) -> tuple[sp.Rational, float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CertificateError("moment components must be finite JSON numbers")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise CertificateError("moment components must be finite")
    candidate = Fraction(numeric).limit_denominator(max_denominator)
    error = abs(numeric - float(candidate))
    if error > tolerance:
        raise CertificateError(
            f"no rational with denominator <= {max_denominator} within tolerance: "
            f"value={numeric!r}, best={candidate}, error={error}"
        )
    return sp.Rational(candidate.numerator, candidate.denominator), error


def _as_gaussian_rational(row: dict[str, Any], max_denominator: int, tolerance: float) -> tuple[sp.Expr, float]:
    if not isinstance(row, dict) or not {"word", "re", "im"}.issubset(row):
        raise CertificateError("each moment requires word, re, and im")
    re_part, re_error = _recover_rational(row["re"], max_denominator, tolerance)
    im_part, im_error = _recover_rational(row["im"], max_denominator, tolerance)
    return re_part + sp.I * im_part, max(re_error, im_error)


def _conj(value: sp.Expr) -> sp.Expr:
    return sp.expand(sp.conjugate(value))


def _is_zero(value: sp.Expr) -> bool:
    return sp.expand(value) == 0


def _matrix_zero(matrix: sp.MatrixBase) -> bool:
    return all(_is_zero(value) for value in matrix)


def _matrix_equal(left: sp.MatrixBase, right: sp.MatrixBase) -> bool:
    return left.shape == right.shape and _matrix_zero(left - right)


def _hermitian(matrix: sp.MatrixBase) -> bool:
    return _matrix_equal(matrix, matrix.adjoint())


def _rational_real(value: sp.Expr, label: str) -> sp.Rational:
    value = sp.expand(value)
    if sp.im(value) != 0 or not sp.re(value).is_Rational:
        raise CertificateError(f"{label} is not an exact rational real number: {value}")
    return sp.Rational(sp.re(value))


def _psd_ldl(matrix: sp.MatrixBase, label: str, *, positive_definite: bool = False) -> dict[str, Any]:
    """Exact pivoted Schur/LDL test for a Gaussian-rational Hermitian matrix."""
    if matrix.rows != matrix.cols or not _hermitian(matrix):
        raise CertificateError(f"{label} is not Hermitian square")
    work = sp.Matrix(matrix)
    pivots: list[sp.Rational] = []
    zero_dimensions = 0
    while work.rows:
        diagonal = [_rational_real(work[i, i], f"{label} diagonal") for i in range(work.rows)]
        if any(value < 0 for value in diagonal):
            raise CertificateError(f"{label} is not positive semidefinite: negative diagonal/Schur pivot")
        pivot_index = next((i for i, value in enumerate(diagonal) if value > 0), None)
        if pivot_index is None:
            if not _matrix_zero(work):
                raise CertificateError(f"{label} has a zero diagonal with a nonzero row")
            zero_dimensions += work.rows
            break
        if pivot_index:
            order = [pivot_index] + [i for i in range(work.rows) if i != pivot_index]
            work = work.extract(order, order)
        pivot = _rational_real(work[0, 0], f"{label} pivot")
        pivots.append(pivot)
        if work.rows == 1:
            work = sp.zeros(0)
        else:
            column = work[1:, 0]
            row = work[0, 1:]
            work = (work[1:, 1:] - column * row / pivot).applyfunc(sp.expand)
            if not _hermitian(work):
                raise CertificateError(f"{label} Schur complement lost Hermiticity")
    if positive_definite and (zero_dimensions or len(pivots) != matrix.rows):
        raise CertificateError(f"{label} is not positive definite")
    return {
        "rank": len(pivots),
        "positive_pivots": [_json_rational(value) for value in pivots],
        "zero_dimensions": zero_dimensions,
    }


def _json_rational(value: sp.Rational) -> list[int]:
    value = sp.Rational(value)
    return [int(value.p), int(value.q)]


def _json_scalar(value: sp.Expr) -> dict[str, list[int]]:
    value = sp.expand(value)
    real = _rational_real(sp.re(value), "matrix real part")
    imag = _rational_real(sp.im(value), "matrix imaginary part")
    return {"re": _json_rational(real), "im": _json_rational(imag)}


def _json_matrix(matrix: sp.MatrixBase) -> list[list[dict[str, list[int]]]]:
    return [[_json_scalar(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _canonical_word_list(value: Any, algebra: WordAlgebra) -> tuple[int, ...]:
    if not isinstance(value, list):
        raise CertificateError("moment word must be a JSON list")
    try:
        return tuple(value)
    except TypeError as exc:
        raise CertificateError("invalid moment word") from exc


class _MomentTable:
    def __init__(self, payload: dict[str, Any], algebra: WordAlgebra,
                 max_denominator: int, tolerance: float):
        rows = payload.get("moments")
        if not isinstance(rows, list) or not rows:
            raise CertificateError("moments must be a nonempty list")
        self.algebra = algebra
        self.canonical: dict[tuple[int, ...], sp.Expr] = {}
        self.supplied: list[tuple[tuple[int, ...], sp.Expr]] = []
        self.max_error = 0.0
        for row in rows:
            value, error = _as_gaussian_rational(row, max_denominator, tolerance)
            raw_word = _canonical_word_list(row["word"], algebra)
            try:
                sign, canonical = algebra.reduce_word(raw_word)
            except ValueError as exc:
                raise CertificateError(str(exc)) from exc
            canonical_value = sp.expand(sign * value)
            if canonical in self.canonical and self.canonical[canonical] != canonical_value:
                raise CertificateError(f"inconsistent duplicate moment for canonical word {canonical}")
            self.canonical[canonical] = canonical_value
            self.supplied.append((raw_word, value))
            self.max_error = max(self.max_error, error)

    def moment(self, word: tuple[int, ...]) -> sp.Expr:
        sign, canonical = self.algebra.reduce_word(word)
        if canonical not in self.canonical:
            raise CertificateError(f"missing recovered moment for canonical word {list(canonical)}")
        return sp.expand(sign * self.canonical[canonical])

    def inner(self, left: tuple[int, ...], right: tuple[int, ...]) -> sp.Expr:
        star_sign, left_star = self.algebra.star(left)
        mul_sign, product = self.algebra.multiply(left_star, right)
        return sp.expand(star_sign * mul_sign * self.moment(product))


def _moment_matrix(table: _MomentTable, word_list: tuple[tuple[int, ...], ...]) -> sp.Matrix:
    return sp.Matrix([[table.inner(left, right) for right in word_list] for left in word_list])


def _independent_basis(table: _MomentTable, candidates: tuple[tuple[int, ...], ...]) -> tuple[list[tuple[int, ...]], sp.Matrix, dict[str, Any]]:
    matrix = _moment_matrix(table, candidates)
    psd = _psd_ldl(matrix, "M_(d-1)")
    _, pivots = matrix.rref()
    basis = [candidates[index] for index in pivots]
    if not basis:
        raise CertificateError("normalization produced a zero moment space")
    metric = sp.Matrix([[table.inner(left, right) for right in basis] for left in basis])
    metric_ldl = _psd_ldl(metric, "basis metric G", positive_definite=True)
    if metric_ldl["rank"] != psd["rank"]:
        raise CertificateError("basis rank does not match M_(d-1) rank")
    return basis, metric, psd


def _recover_impl(payload: dict[str, Any], max_denominator: int, tolerance: float) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise CertificateError("payload must be a JSON object")
    m = payload.get("m")
    degree = payload.get("degree")
    if type(m) is not int or not 1 <= m <= 5:
        raise CertificateError("m must be an integer from 1 through 5")
    if type(degree) is not int or degree < 1:
        raise CertificateError("degree must be a positive integer")
    if type(max_denominator) is not int or max_denominator < 1:
        raise CertificateError("max_denominator must be positive")
    if tolerance < 0:
        raise CertificateError("tolerance must be nonnegative")

    algebra = WordAlgebra(m)
    table = _MomentTable(payload, algebra, max_denominator, tolerance)
    if "transfer_floor" not in payload:
        raise CertificateError("transfer_floor is required for exact lower-bound replay")
    transfer_floor, floor_error = _recover_rational(
        payload["transfer_floor"], max_denominator, tolerance
    )
    if transfer_floor < 0 or transfer_floor > 1:
        raise CertificateError("transfer_floor must recover exactly inside [0,1]")
    table.max_error = max(table.max_error, floor_error)
    if table.moment(()) != 1:
        raise CertificateError("normalization failed: omega(1) is not exactly 1")

    levels = {level: algebra.words(level) for level in (degree - 1, degree, degree + 1)}
    basis, metric, previous_psd = _independent_basis(table, levels[degree - 1])
    metric_inverse = metric.inv()
    rank = len(basis)

    coordinate_cache: dict[tuple[int, ...], sp.Matrix] = {}

    def coordinates(word: tuple[int, ...]) -> sp.Matrix:
        sign, canonical = algebra.reduce_word(word)
        if canonical not in coordinate_cache:
            column = sp.Matrix([table.inner(left, canonical) for left in basis])
            coordinate_cache[canonical] = (metric_inverse * column).applyfunc(sp.expand)
        return sign * coordinate_cache[canonical]

    flat_checks: dict[str, Any] = {}
    for level in (degree - 1, degree, degree + 1):
        word_list = levels[level]
        coords = sp.Matrix.hstack(*(coordinates(word) for word in word_list))
        actual = _moment_matrix(table, word_list)
        predicted = (coords.adjoint() * metric * coords).applyfunc(sp.expand)
        if not _matrix_equal(actual, predicted):
            raise CertificateError(f"M_{level} is not represented by the recovered flat basis")
        flat_checks[str(level)] = {
            "word_count": len(word_list),
            "rank": rank,
            "rank_proof": "exact Gram factorization through positive-definite basis metric",
        }
    if previous_psd["rank"] != rank:
        raise CertificateError("internal previous-rank mismatch")

    generator_tokens = tuple(range(2 * m)) + (T,)
    generator_matrices: dict[int, sp.Matrix] = {}
    for generator in generator_tokens:
        columns = []
        for word in basis:
            sign, product = algebra.multiply((generator,), word)
            columns.append(sign * coordinates(product))
        generator_matrices[generator] = sp.Matrix.hstack(*columns).applyfunc(sp.expand)

    omega = coordinates(())
    quarter_columns = []
    for word in basis:
        sign, image = algebra.quarter(word)
        quarter_columns.append(sign * coordinates(image))
    quarter_matrix = sp.Matrix.hstack(*quarter_columns).applyfunc(sp.expand)
    identity = sp.eye(rank)

    checks: list[dict[str, Any]] = []

    def require(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        if not passed:
            raise CertificateError(f"{name} failed: {detail}")

    require("flat_rank_d_minus_1_d", flat_checks[str(degree - 1)]["rank"] == flat_checks[str(degree)]["rank"],
            f"rank={rank}")
    require("flat_rank_next_degree", flat_checks[str(degree)]["rank"] == flat_checks[str(degree + 1)]["rank"],
            f"rank={rank}")
    metric_ldl = _psd_ldl(metric, "G", positive_definite=True)
    require("positive_metric", metric_ldl["rank"] == rank, f"rank={rank}")

    for generator, matrix in generator_matrices.items():
        require(f"metric_selfadjoint_{generator}",
                _matrix_equal(metric * matrix, matrix.adjoint() * metric),
                "G*A=A^*G")

    for left in range(2 * m):
        for right in range(left, 2 * m):
            anticommutator = (generator_matrices[left] * generator_matrices[right]
                              + generator_matrices[right] * generator_matrices[left])
            expected = 2 * identity if left == right else sp.zeros(rank)
            require(f"CAR_{left}_{right}", _matrix_equal(anticommutator, expected),
                    "{gamma_i,gamma_j}=2 delta_ij")

    transfer = generator_matrices[T]
    transfer_form = (metric * transfer).applyfunc(sp.expand)
    floor_form = (metric * (transfer - transfer_floor * identity)).applyfunc(sp.expand)
    complement_form = (metric * (identity - transfer)).applyfunc(sp.expand)
    transfer_ldl = _psd_ldl(transfer_form, "G*T")
    floor_ldl = _psd_ldl(floor_form, "G*(T-transfer_floor*I)")
    complement_ldl = _psd_ldl(complement_form, "G*(I-T)")
    require("T_positive", transfer_ldl["rank"] >= 0, "G*T is exact PSD")
    require("T_respects_transfer_floor", floor_ldl["rank"] >= 0,
            f"G*(T-{transfer_floor}*I) is exact PSD")
    require("I_minus_T_positive", complement_ldl["rank"] >= 0, "G*(I-T) is exact PSD")
    require("T_fixes_Omega", _matrix_equal(transfer * omega, omega), "T*Omega=Omega")

    require("quarter_metric_unitary",
            _matrix_equal(quarter_matrix.adjoint() * metric * quarter_matrix, metric),
            "R^* G R=G")
    require("quarter_order_four", _matrix_equal(quarter_matrix ** 4, identity), "R^4=I")
    require("quarter_fixes_Omega", _matrix_equal(quarter_matrix * omega, omega), "R*Omega=Omega")
    for generator in generator_tokens:
        sign, image = algebra.quarter((generator,))
        if len(image) != 1:
            raise CertificateError("quarter image of a generator is not a generator")
        require(f"quarter_action_{generator}",
                _matrix_equal(quarter_matrix * generator_matrices[generator],
                              sign * generator_matrices[image[0]] * quarter_matrix),
                "R*A_g=phase*A_quarter(g)*R")

    for generator in generator_tokens:
        matrix = generator_matrices[generator]
        for word in levels[degree]:
            sign, product = algebra.multiply((generator,), word)
            require_name = f"closure_g{generator}_w{list(word)}"
            if not _matrix_equal(matrix * coordinates(word), sign * coordinates(product)):
                raise CertificateError(f"{require_name} failed: generator closure/null stability")
    checks.append({
        "name": "all_generator_closure_through_next_degree",
        "passed": True,
        "detail": f"{len(generator_tokens)} generators on {len(levels[degree])} words",
    })

    vector_cache: dict[tuple[int, ...], sp.Matrix] = {(): omega}

    def represented_vector(raw_word: tuple[int, ...]) -> sp.Matrix:
        if raw_word in vector_cache:
            return vector_cache[raw_word]
        token = raw_word[0]
        suffix = raw_word[1:]
        vector_cache[raw_word] = generator_matrices[token] * represented_vector(suffix)
        return vector_cache[raw_word]

    replay_failures = []
    omega_bra = omega.adjoint() * metric
    for raw_word, expected in table.supplied:
        actual = sp.expand((omega_bra * represented_vector(raw_word))[0])
        if actual != expected:
            replay_failures.append({"word": list(raw_word), "expected": str(expected), "actual": str(actual)})
            if len(replay_failures) == 5:
                break
    require("all_supplied_moments_replay", not replay_failures,
            "all recovered moments replay exactly" if not replay_failures else str(replay_failures))

    named_matrices = {
        ("T" if token == T else f"gamma_{token}"): _json_matrix(matrix)
        for token, matrix in generator_matrices.items()
    }
    result = {
        "schema": CERTIFICATE_SCHEMA,
        "ok": True,
        "verdict": "EXACT_FINITE_POSITIVE_GNS_KERNEL_WITH_EFFECT_T",
        "input": {
            "word_algebra_schema": SCHEMA_ID,
            "m": m,
            "degree": degree,
            "moment_count": len(table.supplied),
            "max_denominator": max_denominator,
            "tolerance": tolerance,
            "max_recovery_error": table.max_error,
            "transfer_floor": _json_rational(transfer_floor),
            "original_sdp_degree": payload.get("original_sdp_degree", payload.get("sdp_degree", degree)),
            "selection_scope": payload.get("selection_scope"),
            "search_choices": {
                key: payload[key] for key in ("objective", "seed", "solver", "solver_status")
                if key in payload
            },
            "declared_assumptions": payload.get("assumptions", []),
            "declared_not_assumed": payload.get("not_assumed", []),
        },
        "basis": {
            "rank": rank,
            "words": [list(word) for word in basis],
            "flat_levels": flat_checks,
        },
        "positivity": {
            "M_d_minus_1": previous_psd,
            "metric": metric_ldl,
            "T": transfer_ldl,
            "T_minus_transfer_floor": floor_ldl,
            "I_minus_T": complement_ldl,
        },
        "checks": checks,
        "exact_matrices": {
            "metric_G": _json_matrix(metric),
            "generators": named_matrices,
            "Omega": _json_matrix(omega),
            "quarter_R": _json_matrix(quarter_matrix),
        },
        "scope": {
            "proved": [
                "exact Gaussian-rational recovery of every supplied moment",
                "flat positive finite GNS representation through degree d+1",
                "CAR, generator closure/null stability, the supplied exact lower bound transfer_floor <= T <= 1, fixed cyclic vector, and quarter covariance",
                "exact replay of every supplied recovered moment",
            ],
            "not_proved": [
                "a unital completely positive time channel or semigroup",
                "physical state/time interpretation or uniqueness",
                "selection of this kernel from P1/P2 or the original seam source",
                "derivation of the objective, seed, transfer_floor, or other recorded search choices",
                "continuum, locality, or T1-T8 closure",
            ],
        },
    }
    return result


def recover(payload: dict[str, Any], max_denominator: int = 64,
            tolerance: float = 1e-6) -> dict[str, Any]:
    """Recover and verify an exact certificate, returning an honest failure object."""
    try:
        return _recover_impl(payload, max_denominator, tolerance)
    except (CertificateError, ValueError, ZeroDivisionError) as exc:
        return {
            "schema": CERTIFICATE_SCHEMA,
            "ok": False,
            "verdict": "EXACT_RECOVERY_FAILED",
            "error": str(exc),
            "scope": {"proved": [], "not_proved": ["the supplied moments define an exact flat positive kernel"]},
        }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-denominator", type=int, default=64)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    result = recover(payload, args.max_denominator, args.tolerance)
    result["provenance"] = {
        "input_path": str(args.input.resolve()),
        "input_sha256": _sha256(args.input),
        "checker_sha256": _sha256(Path(__file__).resolve()),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
