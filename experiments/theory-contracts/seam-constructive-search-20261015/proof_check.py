#!/usr/bin/env python3
"""Exact CPU + Lean certificate gate for crossed-defect phase candidates.

The candidate controls only ten small integers.  It cannot supply Lean
source, imports, theorem statements, assumptions, or proof terms.  The fixed
certificate proves universal all-charge identities over ``Fin 8 -> Int``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Any


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
LEAN_PROJECT = HERE.parents[1] / "lean4-carrier-rigidity"
FIXED_CERTIFICATE = HERE / "lean" / "FixedCertificate.lean"
DOMAIN_SPEC = HERE / "domain_spec.json"
LAKE = Path("/Users/stefanhamann/.elan/bin/lake")
LEAN = Path("/Users/stefanhamann/.elan/toolchains/leanprover--lean4---v4.29.1/bin/lean")
SCHEMA = "tfpt.crossed_defect.v1"
EXPECTED_KEYS = {"schema", "chi_exp", "zeta_exp", "kappa_exp"}
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def atomic_write(path: Path, data: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=path.name + ".",
                                     delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _phase(value: Any, label: str) -> int:
    if type(value) is not int or not 0 <= value <= 3:
        raise ValueError(f"{label} must be an integer in 0..3")
    return value


def validate_schema(candidate: Any) -> dict[str, Any]:
    if not isinstance(candidate, dict):
        raise ValueError("candidate must be a JSON object")
    keys = set(candidate)
    if keys != EXPECTED_KEYS:
        missing = sorted(EXPECTED_KEYS - keys)
        extra = sorted(keys - EXPECTED_KEYS)
        raise ValueError(f"candidate keys differ from fixed schema; missing={missing}, extra={extra}")
    if candidate["schema"] != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA!r}")
    values = candidate["chi_exp"]
    if not isinstance(values, list) or len(values) != 8:
        raise ValueError("chi_exp must contain exactly eight entries")
    chi = [_phase(value, f"chi_exp[{index}]") for index, value in enumerate(values)]
    return {"schema": SCHEMA, "chi_exp": chi,
            "zeta_exp": _phase(candidate["zeta_exp"], "zeta_exp"),
            "kappa_exp": _phase(candidate["kappa_exp"], "kappa_exp")}


def exact_checks(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    """Check the finite coefficient conditions that imply the Lean theorem."""
    chi = candidate["chi_exp"]
    zeta = candidate["zeta_exp"]
    kappa = candidate["kappa_exp"]
    checks: list[tuple[str, bool, Any]] = []
    checks.append(("character_square_generators",
                   all((2 * value - 2) % 4 == 0 for value in chi),
                   [(2 * value - 2) % 4 for value in chi]))
    checks.append(("chi_q", sum(chi) % 4 == 0, sum(chi) % 4))
    difference_residuals = [(chi[i] - chi[j]) % 4 for i in range(8) for j in range(i + 1, 8)]
    pair_residuals = [(chi[i] + chi[j] - 2) % 4 for i in range(8) for j in range(i + 1, 8)]
    checks.append(("d8_difference_roots", not any(difference_residuals), difference_residuals))
    checks.append(("d8_pair_roots", not any(pair_residuals), pair_residuals))
    permutation = [4, 5, 6, 7, 0, 1, 2, 3]
    p_residuals = [(chi[j] - chi[permutation[j]]) % 4 for j in range(8)]
    checks.append(("half_swap_character", not any(p_residuals), p_residuals))
    checks.append(("p_lift_square", 2 * kappa % 4 == 0, 2 * kappa % 4))
    l_minus_b = -sum(chi[4:])
    checks.append(("us_rc_equals_w4_phase", l_minus_b % 4 == 0, l_minus_b % 4))
    # Coefficient form of eps(q,n)/eps(n,q)=(-1)^Q for every integer n.
    comm_coefficients = [(2 * ((7 - j) - j)) % 4 for j in range(8)]
    checks.append(("all_charge_cocycle_coefficients", comm_coefficients == [2] * 8,
                   comm_coefficients))
    return [{"name": name, "passed": passed, "residual": residual}
            for name, passed, residual in checks]


def fixed_source_checks() -> list[dict[str, Any]]:
    spec = json.loads(DOMAIN_SPEC.read_text())
    checks: list[dict[str, Any]] = []
    checks.append({"name": "domain_id", "passed": spec.get("domain_id") == SCHEMA,
                   "residual": spec.get("domain_id")})
    genotype = spec.get("genotype_schema", {})
    checks.append({"name": "strict_genotype_schema",
                   "passed": genotype.get("additionalProperties") is False and
                             set(genotype.get("required", [])) == EXPECTED_KEYS,
                   "residual": genotype.get("required")})
    seen: set[tuple[str, str]] = set()
    for entry in spec.get("source_manifest", []):
        relative = entry.get("path")
        expected = entry.get("sha256")
        key = (str(relative), str(expected))
        if key in seen:
            continue
        seen.add(key)
        path = REPO_ROOT / str(relative)
        actual = sha256_file(path) if path.is_file() else None
        checks.append({"name": f"source_pin:{relative}", "passed": actual == expected,
                       "residual": {"expected": expected, "actual": actual}})
    return checks


def render_candidate_data(candidate: dict[str, Any]) -> str:
    values = ", ".join(str(value) for value in candidate["chi_exp"])
    return f"""/- Generated data literals only. -/
import Mathlib

namespace CandidateData

abbrev Phase := ZMod 4

def chiExpTable : Fin 8 → Phase := ![{values}]
def zetaExp : Phase := {candidate['zeta_exp']}
def kappaExp : Phase := {candidate['kappa_exp']}

end CandidateData
"""


def _parse_axioms(output: str) -> tuple[bool, list[str], str | None]:
    matches = re.findall(r"depends on axioms:\s*\[([^]]*)\]", output)
    if not matches:
        # Lean prints "does not depend on any axioms" for a completely constructive proof.
        if "does not depend on any axioms" in output:
            return True, [], None
        return False, [], "missing #print axioms receipt"
    axioms = sorted({item.strip() for match in matches for item in match.split(",") if item.strip()})
    unexpected = sorted(set(axioms) - ALLOWED_AXIOMS)
    if unexpected:
        return False, axioms, f"unexpected axiom dependencies: {unexpected}"
    return True, axioms, None


def lean_environment(module_dir: Path) -> dict[str, str]:
    """Use the pinned Lake dependency path while compiling isolated sources."""
    query = subprocess.run([str(LAKE), "env", "printenv", "LEAN_PATH"],
                           cwd=LEAN_PROJECT, text=True, capture_output=True,
                           timeout=30, check=False)
    if query.returncode != 0:
        raise RuntimeError(f"cannot resolve Lake LEAN_PATH: {query.stderr.strip()}")
    env = os.environ.copy()
    env["LEAN_PATH"] = str(module_dir) + os.pathsep + query.stdout.strip()
    return env


def _failure(candidate_input: Any, output_dir: Path, stage: str, error: str,
             checks: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    receipt = {"schema": "tfpt.crossed_defect.receipt.v1", "ok": False,
               "stage": stage, "error": error, "checks": checks or []}
    try:
        receipt["candidate_sha256"] = sha256_bytes(canonical_json(candidate_input).encode())
    except (TypeError, ValueError):
        receipt["candidate_sha256"] = None
    atomic_write(output_dir / "receipt.json", json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def check_candidate(candidate: dict, output_dir: Path) -> dict:
    """Validate, render and kernel-check one candidate; always write a receipt."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        validated = validate_schema(candidate)
    except (TypeError, ValueError) as exc:
        return _failure(candidate, output_dir, "schema", str(exc))
    try:
        source_checks = fixed_source_checks()
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return _failure(validated, output_dir, "fixed_sources", str(exc))
    source_failures = [row["name"] for row in source_checks if not row["passed"]]
    if source_failures:
        return _failure(validated, output_dir, "fixed_sources",
                        f"fixed source validation failed: {source_failures}", source_checks)
    checks = source_checks + exact_checks(validated)
    failed = [row["name"] for row in checks if not row["passed"]]
    if failed:
        return _failure(validated, output_dir, "exact_cpu",
                        f"exact conditions failed: {failed}", checks)

    data_source = render_candidate_data(validated)
    forbidden = re.compile(r"\b(axiom|constant|sorry|admit|unsafe|opaque|partial)\b")
    if forbidden.search(data_source):
        return _failure(validated, output_dir, "renderer",
                        "generated literal source contains a forbidden declaration token", checks)
    data_path = output_dir / "CandidateData.lean"
    certificate_path = output_dir / "FixedCertificate.lean"
    atomic_write(data_path, data_source)
    fixed_bytes = FIXED_CERTIFICATE.read_bytes()
    atomic_write(certificate_path, fixed_bytes.decode())
    if certificate_path.read_bytes() != fixed_bytes:
        return _failure(validated, output_dir, "renderer", "fixed certificate copy changed", checks)

    try:
        env = lean_environment(output_dir)
        data_command = [str(LEAN), "-o", "CandidateData.olean", "CandidateData.lean"]
        command = [str(LEAN), "FixedCertificate.lean"]
        data_run = subprocess.run(data_command, cwd=output_dir, env=env, text=True,
                                  capture_output=True, timeout=120, check=False)
        if data_run.returncode == 0:
            run = subprocess.run(command, cwd=output_dir, env=env, text=True,
                                 capture_output=True, timeout=120, check=False)
        else:
            run = data_run
    except (OSError, subprocess.TimeoutExpired) as exc:
        return _failure(validated, output_dir, "lean", f"Lean invocation failed: {exc}", checks)
    combined = run.stdout + "\n" + run.stderr
    axiom_ok, axioms, axiom_error = _parse_axioms(combined)
    signature_seen = "SeamConstructiveCertificate.candidate_certificate" in combined
    ok = run.returncode == 0 and axiom_ok and signature_seen
    version = subprocess.run([str(LEAN), "--version"], cwd=output_dir,
                             text=True, capture_output=True, timeout=30, check=False)
    receipt = {
        "schema": "tfpt.crossed_defect.receipt.v1",
        "ok": ok,
        "stage": "complete" if ok else "lean",
        "candidate_sha256": sha256_bytes(canonical_json(validated).encode()),
        "candidate_class": {
            "orientation": "minus_i" if validated["chi_exp"][0] == 3 else "plus_i",
            "zeta_gauge_exp": validated["zeta_exp"],
            "p_lift_exp": validated["kappa_exp"],
        },
        "checks": checks,
        "fixed_theorem": "SeamConstructiveCertificate.candidate_certificate",
        "signature_seen": signature_seen,
        "axioms": axioms,
        "axiom_error": axiom_error,
        "lean": {"returncode": run.returncode, "candidate_data_returncode": data_run.returncode,
                 "version": version.stdout.strip(),
                 "stdout": run.stdout, "stderr": run.stderr},
        "source_hashes": {
            "proof_check.py": sha256_file(Path(__file__)),
            "domain_spec.json": sha256_file(DOMAIN_SPEC),
            "lean/FixedCertificate.lean": sha256_file(FIXED_CERTIFICATE),
            "CandidateData.lean": sha256_file(data_path),
            "FixedCertificate.copy.lean": sha256_file(certificate_path),
            "lean-toolchain": sha256_file(LEAN_PROJECT / "lean-toolchain"),
            "lakefile.lean": sha256_file(LEAN_PROJECT / "lakefile.lean"),
            "lake-manifest.json": sha256_file(LEAN_PROJECT / "lake-manifest.json"),
        },
        "scope": {
            "proved": [
                "C1-C7 exact phase, character, D8, P-lift, and UsRc relations",
                "C8 universal normal-form associativity for every n,m,l in Z^8 and c,d,f in Bool",
                "C9 partial: Q(Pp+a)-Q(p)=4 for all rational charges, including the half-lattice",
                "C10 fixed rational K_alpha W4 Ward identity, source-channel clock pullback, and mu8 quarter table",
                "fixed four-sector toggle, T-sign, and Arf-coefficient fixture only",
            ],
            "not_proved": [
                "C9 remaining parity-operator transport through U_s, R_c, and every sector",
                "C11 finite-word expectation SOS and positivity",
                "C12 full P_hat/U_s/R_c/W4 action and S/T/Arf compatibility",
                "P1/P2 source selection", "local vertex existence", "event law",
                "C-star completion", "physical state selection", "T1-T8 closure",
            ],
        },
    }
    if not ok:
        receipt["error"] = axiom_error or f"Lean exit={run.returncode}, signature_seen={signature_seen}"
    atomic_write(output_dir / "receipt.json", json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--candidate", required=True, type=Path)
    result.add_argument("--output", required=True, type=Path)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        candidate = json.loads(args.candidate.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        receipt = _failure({}, args.output.resolve(), "input", str(exc))
    else:
        receipt = check_candidate(candidate, args.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0 if receipt["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
