#!/usr/bin/env python3
"""Checker for UR.MINIMAL_RUNTIME.RECONSTRUCTION.01 (minimal-runtime-reconstruction-20261005).

Runs the five parts, requires every step to pass, compares the load-bearing numbers with frozen pins and writes the
aggregate verdict.  Nothing here edits verification/, the ledger, papers or the website.
    python3 -B experiments/theory-contracts/minimal-runtime-reconstruction-20261005/checker.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARTS = ["future_equivalence", "mark_clifford", "filter_uniqueness", "kinetic_reconstruction", "quartic_vertex"]
PINS = {
    "mark_clifford": {"A_gram_diag": 0.75, "A_gram_off": -0.25, "A_rank": 8, "A4_order": 12,
                      "B_ranks": [1, 9, 55, 223, 535, 991, 1645, 2533, 3691]},
    "filter_uniqueness": {"b_max": 1.0},
    "kinetic_reconstruction": {"T_fixed_points": ["-P", "Gamma", "H", "P"], "theta_c_deg": 34.06},
    "quartic_vertex": {"lambda2_multiplicities": [10, 126, 360], "lambda4_singlets": 1, "hermitian_invariant_dim": 5,
                       "physical_couplings": 4},
}

out, ok = {}, True
for part in PARTS:
    r = subprocess.run([sys.executable, "-B", str(HERE / f"{part}.py")], cwd=ROOT, capture_output=True, text=True)
    res = json.loads((HERE / f"{part}_results.json").read_text())
    passed = r.returncode == 0 and res["passed"] == res["total"]
    pin_ok = all(res.get("pins", {}).get(k) == v for k, v in PINS.get(part, {}).items())
    ok &= passed and pin_ok
    out[part] = {"passed": res["passed"], "total": res["total"], "pins_ok": pin_ok,
                 "verdict": res.get("verdict"), "sha256": hashlib.sha256((HERE / f"{part}.py").read_bytes()).hexdigest()[:16]}
    print(f"{'PASS' if passed and pin_ok else 'FAIL'} {part}: {res['passed']}/{res['total']} steps, pins {'ok' if pin_ok else 'MISMATCH'}")

summary = {
    "research_id": "UR.MINIMAL_RUNTIME.RECONSTRUCTION.01",
    "verdict": "NOT_DERIVED" if ok else "CHECKER_FAILED",
    "sub_verdicts": {
        "UR.FUTURE_EQUIVALENCE.FINITERANK.01": "PASS_MATHEMATICAL",
        "UR.CONTINUATION_KERNEL.01": "CONSTRUCTED (three source channels, kernel-only)",
        "UR.RESIDUAL_RANK.01": "STABLE_INTERNAL (A: 8, B-fiber: 9/3, C: 4/16); ADDRESS_UNBOUNDED (B: ~L^3)",
        "SOURCE.MARK_ACTION.01": "CONSTRUCTED",
        "SOURCE.MARK_CLIFFORD.01": "NO_GO_FROM_CURRENT_SOURCE",
        "SOURCE.FILTER.UNIQUENESS.01": "CONDITIONAL_UNIQUE",
        "SOURCE.KINETIC.RECONSTRUCTION.01": "NOT_DERIVED",
        "SOURCE.QUARTIC.VERTEX.01": "NOT_SELECTED",
    },
    "minimal_missing_inputs": [
        "spin: one Z2 projective class bit for a T-symmetric Cl(4) coin (AX.MARKSIGN, g = I); for the chiral Theta_1 only the "
        "two-component internal-singlet coin (AX.WEYLCOIN; spin-1/2 class IR-automatic, mark-projective-lift-20261005)",
        "kinetic: T-breaking (forced by minimal doubling) + maximal-residual-symmetry premise (one mark) + fundamental (non-IR) "
        "conditions (ii)-(vi), (viii)",
        "filter: reality / reversible Markov averaging premise (b = 0)",
        "quartic: 4 physical couplings, Theta_1 = one ray; source response zero",
    ],
    "parts": out,
}
(HERE / "results.json").write_text(json.dumps(summary, indent=2))
print(f"{summary['research_id']}: {summary['verdict']}")
sys.exit(0 if ok else 1)
