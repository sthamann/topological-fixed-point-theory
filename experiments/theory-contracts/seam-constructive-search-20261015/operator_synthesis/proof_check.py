#!/usr/bin/env python3
"""Kernel-check the fixed generic operator-transport certificate.

This checker deliberately does not certify a synthesized matrix candidate.
The exact solver must first verify that candidate's matrix equations; this
file certifies only the universal Lean implications used after that step.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
LEAN_SOURCE = HERE / "lean" / "OperatorTransport.lean"
LEAN_PROJECT = REPO / "experiments" / "lean4-carrier-rigidity"
LAKE = Path.home() / ".elan" / "bin" / "lake"
DEFAULT_RECEIPT = HERE / "lean-receipt.json"

EXPECTED_THEOREMS = [
    "TFPT.OperatorTransport.RingExpr.map_eval",
    "TFPT.OperatorTransport.intertwines_add",
    "TFPT.OperatorTransport.intertwines_neg",
    "TFPT.OperatorTransport.intertwines_mul",
    "TFPT.OperatorTransport.intertwines_pow",
    "TFPT.OperatorTransport.intertwines_eval",
    "TFPT.OperatorTransport.intertwines_target_eq_conjugate",
    "TFPT.OperatorTransport.intertwines_star_of_unitary",
    "TFPT.OperatorTransport.conjugateBy_star",
    "TFPT.OperatorTransport.conjugateBy_positive_square",
    "TFPT.OperatorTransport.conjugate_intertwines_eval",
    "TFPT.OperatorTransport.imaginary_generator_changes_sign",
]
ALLOWED_AXIOMS = {"propext", "Quot.sound"}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def _parse_axiom_audit(output: str) -> tuple[dict[str, list[str]], list[str]]:
    dependencies: dict[str, list[str]] = {}
    for line in output.splitlines():
        with_axioms = re.fullmatch(r"'([^']+)' depends on axioms: \[(.*)\]", line.strip())
        if with_axioms:
            name, body = with_axioms.groups()
            dependencies[name] = [item.strip() for item in body.split(",") if item.strip()]
            continue
        without_axioms = re.fullmatch(r"'([^']+)' does not depend on any axioms", line.strip())
        if without_axioms:
            dependencies[without_axioms.group(1)] = []
    errors: list[str] = []
    missing = [name for name in EXPECTED_THEOREMS if name not in dependencies]
    if missing:
        errors.append(f"missing #print axioms output: {missing}")
    unexpected = sorted(
        {axiom for values in dependencies.values() for axiom in values} - ALLOWED_AXIOMS
    )
    if unexpected:
        errors.append(f"unexpected axioms: {unexpected}")
    if "sorryAx" in output:
        errors.append("sorryAx occurred in Lean output")
    return dependencies, errors


def _version(command: list[str]) -> str:
    run = subprocess.run(
        command, cwd=LEAN_PROJECT, text=True, capture_output=True,
        timeout=30, check=False,
    )
    return (run.stdout or run.stderr).strip()


def check_transport_kernel(receipt_path: Path = DEFAULT_RECEIPT) -> dict[str, Any]:
    """Run Lean on the immutable generic theorem file and write an audit receipt."""
    command = [str(LAKE), "env", "lean", str(LEAN_SOURCE)]
    run = subprocess.run(
        command, cwd=LEAN_PROJECT, text=True, capture_output=True,
        timeout=120, check=False,
    )
    combined = run.stdout + "\n" + run.stderr
    dependencies, audit_errors = _parse_axiom_audit(combined)
    source_text = LEAN_SOURCE.read_text()
    forbidden_source = [
        token for token in ("sorry", "admit", "axiom ") if token in source_text
    ]
    if forbidden_source:
        audit_errors.append(f"forbidden source tokens: {forbidden_source}")
    if run.returncode != 0:
        audit_errors.append(f"Lean exited with {run.returncode}")

    receipt: dict[str, Any] = {
        "schema": "tfpt.operator_transport.lean_receipt.v1",
        "ok": not audit_errors,
        "verdict": "KERNEL_CHECKED_GENERIC_OPERATOR_TRANSPORT" if not audit_errors
                   else "LEAN_CERTIFICATE_FAILED",
        "command": command,
        "exit_code": run.returncode,
        "source": {
            "path": str(LEAN_SOURCE.relative_to(REPO)),
            "sha256": sha256_file(LEAN_SOURCE),
        },
        "checker_sha256": sha256_file(Path(__file__).resolve()),
        "environment": {
            "lean": _version([str(LAKE), "env", "lean", "--version"]),
            "lake": _version([str(LAKE), "--version"]),
            "lean_toolchain_sha256": sha256_file(LEAN_PROJECT / "lean-toolchain"),
            "lakefile_sha256": sha256_file(LEAN_PROJECT / "lakefile.lean"),
            "lake_manifest_sha256": sha256_file(LEAN_PROJECT / "lake-manifest.json"),
        },
        "theorems": EXPECTED_THEOREMS,
        "axiom_dependencies": dependencies,
        "allowed_axioms": sorted(ALLOWED_AXIOMS),
        "sorryAx_present": "sorryAx" in combined,
        "audit_errors": audit_errors,
        "scope": {
            "proved": [
                "simultaneous intertwining of checked generators extends to every finite noncommutative integer-polynomial expression",
                "one checked intertwining relation extends to every natural power",
                "unitary-like transport identifies the target by conjugation and preserves adjoints",
                "unitary-like conjugation transports algebraic positive squares star(a)*a",
                "an explicit conjugation automorphism extends generator checks to every generated expression",
                "a conjugation-odd imaginary generator changes sign in the conjugate-linear branch",
            ],
            "not_proved": [
                "any concrete solver candidate satisfies the theorem hypotheses",
                "existence or uniqueness of a linear or antilinear intertwiner",
                "analytic positivity, logarithms, operator domains, or C-star completion",
                "selection of the source representation, state, dynamics, or P1/P2 origin",
            ],
        },
        "integration_rule": (
            "A candidate receipt may invoke these theorems only after an independent exact "
            "matrix checker has verified the listed generator, unitary, star, and conjugation equations."
        ),
        "lean_stdout": run.stdout,
        "lean_stderr": run.stderr,
    }
    _atomic_json(Path(receipt_path), receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    receipt = check_transport_kernel(args.receipt.resolve())
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if receipt["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
