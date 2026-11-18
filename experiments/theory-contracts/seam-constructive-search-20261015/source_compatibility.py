"""Assess a proved algebraic candidate against the unchanged TFPT source.

This is deliberately separate from ``proof_check.py``.  The Lean receipt says
that the crossed-defect identities hold.  This module asks the later, physical
question whether that algebra can be realized on the *unchanged* source
Hilbert space with its integer Cartans and its fixed quadratic time.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DOMAIN_SPEC = HERE / "domain_spec.json"
SOURCE_RELATIVE = Path("experiments/theory-contracts/quotient-locality-e8-20261014/PROOF.txt")
SOURCE_PATH = REPO / SOURCE_RELATIVE
SOURCE_SHA256 = "dfa4b094f660d47bc28cf2729e1c12eb29ae0b79cdfda36e72cb40f9406f5714"

EXPECTED_DIMENSION = 8
EXPECTED_S_TWICE = [1] * 8
EXPECTED_A = [1, 1, 1, 1, 0, 0, 0, 0]
EXPECTED_B = [0, 0, 0, 0, 1, 1, 1, 1]
EXPECTED_R0_TWICE = [1, 1, 1, 1, -1, -1, -1, -1]
EXPECTED_TIME = {
    "kind": "fixed_input_extended_to_all_charges",
    "charges": "p=n+c*s, n in Z^8, c in {0,1}",
    "K_alpha": "K_alpha(p,N)=N+|p|^2/2+delta_alpha dot p",
    "delta_alpha": "delta_alpha[j]=(j+alpha)/8-1/2",
    "alpha_values": ["1/4", "3/4"],
    "W4_ward": "[K_alpha,W4]=W4*(Q/2+3/4+alpha/2)",
    "quarter_character_mu8_exp": {
        "alpha_1_4": "7+4*(Q mod 2) mod 8",
        "alpha_3_4": "1+4*(Q mod 2) mod 8",
    },
}


class SourceEvidenceError(ValueError):
    """The pinned obstruction source is absent or has changed."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _candidate_sha256(candidate: dict[str, Any]) -> str:
    payload = json.dumps(candidate, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode()).hexdigest()


def _scope_checks(spec: dict[str, Any]) -> list[dict[str, Any]]:
    fixed = spec.get("fixed_source_data", {})
    vectors = fixed.get("vectors", {})
    time_data = fixed.get("time", {})
    definitions = fixed.get("definitions", {})
    checks: list[tuple[str, Any, Any]] = [
        ("dimension", fixed.get("dimension"), EXPECTED_DIMENSION),
        ("basis", fixed.get("basis"), [f"e{j}" for j in range(EXPECTED_DIMENSION)]),
        ("integer_charge", definitions.get("charge"), "Q(n)=sum_j n_j"),
        ("s_twice", vectors.get("s_twice"), EXPECTED_S_TWICE),
        ("a", vectors.get("a"), EXPECTED_A),
        ("b", vectors.get("b"), EXPECTED_B),
        ("r0_twice", vectors.get("r0_twice"), EXPECTED_R0_TWICE),
    ]
    checks.extend((f"time.{key}", time_data.get(key), expected) for key, expected in EXPECTED_TIME.items())

    s_twice = vectors.get("s_twice")
    b = vectors.get("b")
    r0_twice = vectors.get("r0_twice")
    a = vectors.get("a")
    dimension = fixed.get("dimension")
    vector_shape_ok = all(
        isinstance(value, list) and isinstance(dimension, int) and len(value) == dimension
        for value in (s_twice, b, r0_twice, a)
    )
    r0_definition = vector_shape_ok and all(r == s - 2 * bv for r, s, bv in zip(r0_twice, s_twice, b))
    a_definition = vector_shape_ok and all(2 * av == s + r for av, s, r in zip(a, s_twice, r0_twice))
    checks.extend([
        ("r0=s-b", r0_definition, True),
        ("a=s+r0", a_definition, True),
    ])
    return [
        {"name": name, "passed": actual == expected, "actual": actual, "expected": expected}
        for name, actual, expected in checks
    ]


def _rational_k_checks(spec: dict[str, Any]) -> list[dict[str, str]]:
    fixed = spec["fixed_source_data"]
    r0 = [Fraction(value, 2) for value in fixed["vectors"]["r0_twice"]]
    dimension = fixed["dimension"]
    norm_squared = sum(value * value for value in r0)
    results = []
    for alpha_text in fixed["time"]["alpha_values"]:
        alpha = Fraction(alpha_text)
        delta = [Fraction(j, 1) / dimension + alpha / dimension - Fraction(1, 2)
                 for j in range(dimension)]
        delta_dot_r0 = sum(d * p for d, p in zip(delta, r0))
        k_zero = Fraction(0)
        k_r0 = norm_squared / 2 + delta_dot_r0
        results.append({
            "alpha": alpha_text,
            "norm_r0_squared": _fraction_text(norm_squared),
            "quadratic_term": _fraction_text(norm_squared / 2),
            "delta_dot_r0": _fraction_text(delta_dot_r0),
            "K_zero": _fraction_text(k_zero),
            "K_r0": _fraction_text(k_r0),
        })
    return results


def assess_source_compatibility(
    candidate: dict[str, Any],
    algebra_receipt: dict[str, Any],
    *,
    domain_spec_path: Path = DOMAIN_SPEC,
    source_path: Path = SOURCE_PATH,
) -> dict[str, Any]:
    """Return a source assessment without changing the algebraic receipt."""
    source_path = Path(source_path)
    if not source_path.is_file():
        raise SourceEvidenceError(f"Pinned source evidence is missing: {source_path}")
    actual_source_sha = _sha256(source_path)
    if actual_source_sha != SOURCE_SHA256:
        raise SourceEvidenceError(
            f"Pinned source evidence changed: expected {SOURCE_SHA256}, observed {actual_source_sha}"
        )

    spec = json.loads(Path(domain_spec_path).read_text())
    scope_checks = _scope_checks(spec)
    same_fixed_scope = all(row["passed"] for row in scope_checks)
    candidate_sha = _candidate_sha256(candidate)
    receipt_candidate_sha = algebra_receipt.get("candidate_sha256")
    receipt_matches_candidate = receipt_candidate_sha == candidate_sha
    algebra_ok = algebra_receipt.get("ok") is True and receipt_matches_candidate
    status = "EXCLUDED_IN_UNCHANGED_SOURCE" if algebra_ok and same_fixed_scope else "UNKNOWN"
    k_checks = _rational_k_checks(spec) if same_fixed_scope else []

    return {
        "schema": "tfpt.source_compatibility.v1",
        "candidate_sha256": candidate_sha,
        "algebra_certificate_ok": algebra_ok,
        "algebra_receipt_candidate_match": receipt_matches_candidate,
        "source_equivalence_status": status,
        "status_reason": (
            "The proved crossed-defect algebra requires a nonzero half-Cartan-charge operator and an "
            "additional K_alpha ground state, both excluded on the unchanged source Hilbert space."
            if status == "EXCLUDED_IN_UNCHANGED_SOURCE"
            else "The pinned unchanged-source exclusion was not transferred because the algebra receipt "
                 "or at least one fixed source/time premise does not match."
        ),
        "half_charge_operator_obstruction": {
            "status": "APPLIES" if status == "EXCLUDED_IN_UNCHANGED_SOURCE" else "NOT_TRANSFERRED",
            "statement": (
                "Integer source Cartans satisfy Ad(exp(2*pi*i*N_k))=identity; a nonzero operator with "
                "Cartan charge +/-1/2 would instead acquire phase -1 and therefore must vanish."
            ),
            "candidate_interpretation": (
                "The formal defect generator X can define an algebraic extension, but it is not an "
                "operator realization on the unchanged source Hilbert space."
            ),
        },
        "ground_space": {
            "original_source_multiplicity": 1 if same_fixed_scope else None,
            "stipulated_e8_extension_multiplicity": 2 if same_fixed_scope else None,
            "rational_connection_check": k_checks,
            "connection_check_scope": (
                "The displayed K(0)=K(r0)=0 identities only connect the candidate's fixed r0 and quarter "
                "branches to the source theorem. The full-lattice exclusion comes from the pinned analytic "
                "coordinate-minimization and microscopic-gap proof, not from this finite calculation."
            ),
        },
        "scope_checks": scope_checks,
        "source": {
            "research_id": "UR.SEAM.QUOTIENT_LOCALITY_E8.01",
            "path": SOURCE_RELATIVE.as_posix(),
            "sha256": SOURCE_SHA256,
            "citations": [
                {"section": "3", "lines": "82-103", "claim": "integer-Cartan half-charge operator obstruction"},
                {"section": "4", "lines": "105-147", "claim": "same-K ground-space multiplicities 1 versus 2"},
                {"section": "5", "lines": "149-172", "claim": "microscopic source gap and unique sea"},
            ],
        },
        "conditional_scope": {
            "requires": [
                "the same microscopic source Hilbert space",
                "the same eight integer occupation Cartans",
                "the same normal ordering and source vacuum constant",
                "the fixed K_alpha quadratic time with alpha in {1/4,3/4}",
                "the fixed eight-dimensional charges s and r0",
            ],
            "does_not_claim": [
                "a universal no-go for other representations, sources, scalings, or added sectors",
                "a derivation or refutation of P1/P2 source selection",
                "that failure of unchanged-source equivalence invalidates the algebraic Lean certificate",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--domain-spec", type=Path, default=DOMAIN_SPEC)
    parser.add_argument("--source", type=Path, default=SOURCE_PATH)
    args = parser.parse_args()
    assessment = assess_source_compatibility(
        json.loads(args.candidate.read_text()),
        json.loads(args.receipt.read_text()),
        domain_spec_path=args.domain_spec,
        source_path=args.source,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(assessment, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
