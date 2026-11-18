#!/usr/bin/env python3
"""Exact inverse synthesis of simultaneous linear/antilinear operator maps.

This solves a declared finite representation problem, not the physical choice
of that representation. No desired map or fitted Hamiltonian enters the input.
Counterexamples add whole operator relations. Every round has its own key.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import re

import sympy as s

RATIONAL = re.compile(r"^[+-]?\d+(?:/[1-9]\d*)?$")


def scalar(value):
    if type(value) is int:
        return s.Integer(value)
    if isinstance(value, str) and RATIONAL.fullmatch(value):
        return s.Rational(value)
    if isinstance(value, list) and len(value) == 2:
        if any(isinstance(x, list) for x in value):
            raise ValueError("Gaussian entries need two rational scalars")
        return scalar(value[0]) + s.I * scalar(value[1])
    raise ValueError("Use exact integers, rational strings, or [real, imaginary]; no floats/code")


def matrix(rows):
    if not rows or not all(isinstance(row, list) and len(row) == len(rows) for row in rows):
        raise ValueError("Operator must be a nonempty square matrix")
    return s.Matrix([[scalar(v) for v in row] for row in rows])


def dump_matrix(m):
    return [[[str(s.re(v)), str(s.im(v))] for v in m.row(i)] for i in range(m.rows)]


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_problem(path):
    raw = json.loads(Path(path).read_text())
    if raw.get("schema") != "tfpt.operator_intertwiner.v1":
        raise ValueError("Unsupported problem schema")
    if not raw.get("premises") or not raw.get("relations"):
        raise ValueError("Explicit premises and operator relations are required")
    if "target_map" in raw:
        raise ValueError("A target map is not a synthesis input")
    ids, relations, dimension = set(), [], None
    for row in raw["relations"]:
        name = row["id"]
        if not isinstance(name, str) or name in ids or "::" in name:
            raise ValueError("Relation IDs must be unique strings")
        ids.add(name)
        a, b = matrix(row["source"]), matrix(row["target"])
        dimension = dimension or a.rows
        if a.shape != b.shape or a.rows != dimension:
            raise ValueError("This solver handles equal finite dimensions only")
        # Adjoint closure is necessary for the polar-normalization theorem.
        relations.append((name, a, b))
        relations.append((name + "::adjoint", a.H, b.H))
    return raw, relations, dimension


def equation(a, b, anti):
    """Row-major vectorization of X conjugate(a) - b X, exact over Q(i)."""
    if anti:
        a = a.conjugate()
    n = a.rows
    return s.kronecker_product(s.eye(n), a.T) - s.kronecker_product(b, s.eye(n))


def restrict(basis, equations):
    for eq in equations:
        kernel = (eq * basis).nullspace()
        if not kernel:
            return s.zeros(basis.rows, 0)
        basis = basis * s.Matrix.hstack(*kernel)
    return basis


def choose_full_rank(basis, n, trials):
    """Bounded deterministic search; failure is NEVER a no-isomorphism claim."""
    columns = basis.cols
    choices = [tuple(int(i == j) for i in range(columns)) for j in range(columns)]
    choices += [tuple(1 for _ in range(columns))]
    rng = random.Random(971)
    for _ in range(trials):
        choices.append(tuple(rng.choice((-2, -1, 0, 1, 2)) for _ in range(columns)))
    for coeff in choices:
        flat = basis * s.Matrix(coeff)
        candidate = s.Matrix(n, n, list(flat))
        if candidate.det() != 0:
            return candidate, list(coeff)
    return None, None


def residual(candidate, a, b, anti):
    return s.simplify(candidate * (a.conjugate() if anti else a) - b * candidate)


def conflict_core(relations, active, anti, n):
    """Deletion-minimal exact zero-intertwiner core, within the active list."""
    core = list(active)
    for idx in list(core):
        rest = [j for j in core if j != idx]
        basis = restrict(s.eye(n*n), [equation(relations[j][1], relations[j][2], anti) for j in rest])
        if basis.cols == 0:
            core = rest
    return [relations[j][0] for j in core]


def solve_branch(raw, relations, n, anti=False, trials=128):
    if n > 16:
        raise ValueError("Generic exact solver limit is 16; use the CAR word solver for larger sparse sources")
    active, rounds, basis = [], [], s.eye(n*n)
    code_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    actual_relations = [{"id": name, "source": dump_matrix(a), "target": dump_matrix(b)}
                        for name, a, b in relations]
    for _ in range(len(relations) + 1):
        round_key = canonical_hash({"problem": raw, "relations": actual_relations, "anti": anti,
                                    "active": active, "code": code_sha, "coefficient_trials": trials})
        if not basis.cols:
            return {"status": "no_nonzero_intertwiner", "antilinear": anti,
                    "core": conflict_core(relations, active, anti, n), "rounds": rounds,
                    "scope": "only the supplied finite representation pair and declared map type"}
        candidate, coeff = choose_full_rank(basis, n, trials)
        if candidate is None:
            # Add unused equations first: this can expose a true zero kernel.
            missing = next((j for j in range(len(relations)) if j not in active), None)
            if missing is None:
                return {"status": "full_rank_search_inconclusive", "antilinear": anti,
                        "intertwiner_dimension": basis.cols, "rounds": rounds,
                        "scope": "bounded coefficient attempts; not a proof of absence"}
            violation = missing
            witness = {"kind": "add_relation_after_bounded_rank_search"}
        else:
            violation, witness = None, None
            for j, (name, a, b) in enumerate(relations):
                diff = residual(candidate, a, b, anti)
                nonzero = next(((r, c) for r in range(n) for c in range(n) if diff[r, c] != 0), None)
                if nonzero is not None:
                    violation = j
                    r, c = nonzero
                    witness = {"kind": "exact_matrix_entry", "relation": name, "row": r,
                               "column": c, "residual": str(diff[r, c])}
                    break
            if violation is None:
                # One successful candidate does not determine the dimension of
                # the full solution family. Finish all equations for that claim.
                basis = restrict(basis, [equation(a, b, anti) for _, a, b in relations])
                gram = s.simplify(candidate.H * candidate)
                gram_commutes = all(s.simplify(gram * (a.conjugate() if anti else a)
                                               - (a.conjugate() if anti else a) * gram) == s.zeros(n)
                                    for _, a, _ in relations)
                if not gram_commutes:
                    raise AssertionError("Adjoint-closed intertwiner lost its polar commutant")
                c = gram[0, 0]
                scalar_gram = gram == c * s.eye(n)
                unitary = candidate / s.sqrt(c) if scalar_gram else None
                if unitary is not None:
                    assert s.simplify(unitary.H * unitary) == s.eye(n)
                    assert all(residual(unitary, a, b, anti) == s.zeros(n) for _, a, b in relations)
                return {"status": "exact_intertwiner", "antilinear": anti,
                        "rounds": rounds, "final_round_key": round_key,
                        "intertwiner_dimension": basis.cols, "coefficients": coeff,
                        "rational_intertwiner": dump_matrix(candidate), "gram": dump_matrix(gram),
                        "normalization": "U=S/sqrt(c)" if scalar_gram else "U=S*(S†S)^(-1/2)",
                        "normalization_scalar": str(c) if scalar_gram else None,
                        "unitary_matrix": dump_matrix(unitary) if unitary is not None else None,
                        "gram_positive_definite": candidate.det() != 0,
                        "gram_commutes_with_full_source_star_algebra": gram_commutes,
                        "all_relation_residuals_zero": True,
                        "scope": "simultaneous finite star-algebra transport; no source-selection proof"}
        if violation in active:
            raise AssertionError("Already learned exact constraint failed")
        rounds.append({"key": round_key, "dimension_before": basis.cols,
                       "added_relation": relations[violation][0], "counterexample": witness})
        active.append(violation)
        basis = restrict(basis, [equation(relations[violation][1], relations[violation][2], anti)])
    raise AssertionError("Finite CEGIS failed to terminate")


def run(path, output=None, trials=128):
    raw, relations, n = load_problem(path)
    result = {"schema": "tfpt.operator_intertwiner.result.v1", "problem": raw.get("id"),
              "problem_sha256": canonical_hash(raw), "premises": raw["premises"],
              "source_selection_proved": False, "physical_solution": False,
              "branches": [solve_branch(raw, relations, n, anti, trials) for anti in (False, True)]}
    if output:
        Path(output).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coefficient-trials", type=int, default=128)
    args = parser.parse_args()
    if args.coefficient_trials < 0:
        parser.error("coefficient trials must be nonnegative")
    data = run(args.input, args.output, args.coefficient_trials)
    print(json.dumps({"problem": data["problem"], "branches": [
        {"antilinear": b["antilinear"], "status": b["status"], "rounds": len(b["rounds"])}
        for b in data["branches"]], "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
