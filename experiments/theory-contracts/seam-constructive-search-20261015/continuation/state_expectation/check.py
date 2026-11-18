#!/usr/bin/env python3
"""Compile the isolated C11 theorem and audit its axiom receipt."""
from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import re
import subprocess
import sys


HERE = Path(__file__).resolve().parent
LEAN_PROJECT = HERE.parents[3] / "lean4-carrier-rigidity"
LAKE = Path("/Users/stefanhamann/.elan/bin/lake")
LEAN = Path("/Users/stefanhamann/.elan/toolchains/leanprover--lean4---v4.29.1/bin/lean")
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
REQUIRED_THEOREMS = {
    "SeamStateExpectation.basisMul_assoc",
    "SeamStateExpectation.basis_unitary_left",
    "SeamStateExpectation.basis_unitary_right",
    "SeamStateExpectation.crossedMul_assoc",
    "SeamStateExpectation.crossedStar_involutive",
    "SeamStateExpectation.crossedStar_mul",
    "SeamStateExpectation.crossedX_unitary_left",
    "SeamStateExpectation.crossedX_unitary_right",
    "SeamStateExpectation.expectation_star_mul",
    "SeamStateExpectation.induced_positive",
    "SeamStateExpectation.dual_invariant_extension_unique",
    "SeamStateExpectation.induced_time_invariant",
}


def main() -> int:
    query = subprocess.run([str(LAKE), "env", "printenv", "LEAN_PATH"],
                           cwd=LEAN_PROJECT, text=True, capture_output=True,
                           timeout=30, check=False)
    if query.returncode:
        raise RuntimeError(query.stderr)
    env = os.environ.copy()
    env["LEAN_PATH"] = query.stdout.strip()
    run = subprocess.run([str(LEAN), "StateExpectation.lean"], cwd=HERE, env=env,
                         text=True, capture_output=True, timeout=120, check=False)
    output = run.stdout + run.stderr
    matches = re.findall(r"depends on axioms:\s*\[([^]]*)\]", output)
    axioms = sorted({a.strip() for match in matches for a in match.split(",") if a.strip()})
    missing_theorems = sorted(name for name in REQUIRED_THEOREMS if name not in output)
    receipt = {
        "ok": run.returncode == 0 and not missing_theorems and
              set(axioms) <= ALLOWED_AXIOMS and "sorryAx" not in output,
        "returncode": run.returncode,
        "axioms": axioms,
        "missing_theorems": missing_theorems,
        "stdout": run.stdout,
        "stderr": run.stderr,
        "scope": {
            "lean_proved": ["actual ordered Z^8 cocycle associativity and unitary basis monomials",
                       "general CrossedData normal-form associativity, star laws, and unitary X",
                       "general CrossedData finite-word SOS identity and induced positivity",
                       "dual-Z2-invariant extension uniqueness",
                       "conditional grade-preserving time invariance"],
            "analytical_bridge": [
                "construct A as finite complex sums of ordered-cocycle E_n",
                "derive every CrossedData field from production chi_add, chi_q, chi-square/epsilon and basis unitarity",
                "extend the verified basis identities by finite bilinearity and conjugate-linearity",
            ],
            "not_proved": ["C-star completion", "state selection", "KMS/Gibbs property",
                           "kernel-level concrete Finsupp CrossedData instantiation",
                           "unconditional K_alpha invariance", "P1/P2 origin"],
        },
        "source_sha256": hashlib.sha256((HERE / "StateExpectation.lean").read_bytes()).hexdigest(),
        "bridge_sha256": hashlib.sha256((HERE / "SOURCE_BRIDGE.md").read_bytes()).hexdigest(),
    }
    (HERE / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0 if receipt["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
