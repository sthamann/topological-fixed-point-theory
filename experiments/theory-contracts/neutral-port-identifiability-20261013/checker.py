#!/usr/bin/env python3
"""Exact scoped rejection of the selected-neutral-R route to D5's 5/4.

This imports the existing TFPT source construction read-only.  It does not
test the whole TFPT claim.  It tests what information survives the declared
rank-one P0 port built from the selected neutral state R.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import sys
from pathlib import Path

import sympy as sp


REPO = Path(__file__).resolve().parents[3]
SOURCE = REPO / "tfpt_explorer/neutral_source_response.py"
sys.path.insert(0, str(REPO))
neutral = importlib.import_module("tfpt_explorer.neutral_source_response")

native = neutral._native_neutral_branch()
p5, p2, p1 = (native[key] for key in ("P5", "P2", "P1"))
gram = native["gram"]
states = native["states"]
norms = native["norms"]

# Exact sector coefficients in B48 and in its stress-orthogonal remainder R.
b48 = native["A48"]
r_matrix = native["R_matrix"]
projectors = [p5, p2, p1]
ranks = [5, 2, 1]


def sector_coefficients(matrix: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [sp.simplify(sp.trace(projector * matrix) / rank)
         for projector, rank in zip(projectors, ranks)]
    )


def state_norm_squared(vector: sp.Matrix) -> sp.Expr:
    return sp.simplify((vector.T * gram * vector)[0])


b48_coefficients = sector_coefficients(b48)
r_coefficients = sector_coefficients(r_matrix)

# A distinct Virasoro-primary sector contrast with the same exact norm.
# It is a reflection of R inside 5a+2b+c=0.  This is an information-loss
# witness for the port; it is not claimed to be another selected 48-root branch.
alternate_coefficients = sp.Matrix(
    [sp.Rational(16, 7), -sp.Rational(12, 7), -8]
)
alternate_matrix = (
    alternate_coefficients[0] * p5
    + alternate_coefficients[1] * p2
    + alternate_coefficients[2] * p1
)
alternate_state = neutral._quadratic_state(alternate_matrix)

p0_vector = sp.ones(4, 1) / 2
p0 = p0_vector * p0_vector.T
r_port = states["R"] * p0_vector.T
alternate_port = alternate_state * p0_vector.T
r_port_compression = sp.simplify(r_port.T * gram * r_port)
alternate_port_compression = sp.simplify(
    alternate_port.T * gram * alternate_port
)

checks = {
    "B48_coefficients_are_D5_12_A2_16_U1_4":
        b48_coefficients == sp.Matrix([12, 16, 4]),
    "selected_R_coefficients_are_D5_0_A2_4_U1_minus8":
        r_coefficients == sp.Matrix([0, 4, -8]),
    "selected_R_D5_coefficient_is_exactly_zero": r_coefficients[0] == 0,
    "selected_R_norm_squared_is_48": norms["R"] == 48,
    "B48_norm_squared_is_624": norms["B48"] == 624,
    "selected_R_is_stress_orthogonal":
        5 * r_coefficients[0] + 2 * r_coefficients[1] + r_coefficients[2] == 0,
    "alternate_is_stress_orthogonal":
        5 * alternate_coefficients[0]
        + 2 * alternate_coefficients[1]
        + alternate_coefficients[2] == 0,
    "alternate_has_nonzero_D5_coefficient": alternate_coefficients[0] != 0,
    "alternate_norm_squared_is_48": state_norm_squared(alternate_state) == 48,
    "rank_one_P0_compressions_are_equal":
        r_port_compression == alternate_port_compression == 48 * p0,
}

result = {
    "research_id": "UR.SOURCE.P0_D5_IDENTIFIABILITY.COLD_CHECK.01",
    "verdict": "SCOPED_ROUTE_REJECTED",
    "scope": (
        "The declared selected R -> rank-one P0 port cannot carry or identify "
        "the D5 5/4 factor. This does not reject other carrier-resolved, charged, "
        "or defect responses and is not a whole-claim counterexample."
    ),
    "source": {
        "path": str(SOURCE.relative_to(REPO)),
        "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "lines": {
            "native_48_branch_and_remainder": "72-110",
            "selected_port": "114-132",
            "registered_checks": "149-177",
            "reported_operator_identities": "211-225",
        },
        "atlas_registration": (
            "During the initial cold-agent test, no direct source record was returned after "
            "exact-path, function-name, and formula searches."
        ),
    },
    "registered_context_sources": [
        "source:45d734c654de8fb46a8a",
        "source:ad8a961941e631e265d7",
        "source:8cbcc53b7f027d3211a6",
    ],
    "B48_sector_coefficients_D5_A2_U1": list(map(str, b48_coefficients)),
    "R_sector_coefficients_D5_A2_U1": list(map(str, r_coefficients)),
    "R_norm_squared": str(norms["R"]),
    "B48_norm_squared": str(norms["B48"]),
    "equal_norm_primary_witness_coefficients_D5_A2_U1":
        list(map(str, alternate_coefficients)),
    "equal_norm_primary_witness_norm_squared":
        str(state_norm_squared(alternate_state)),
    "P0_compression": [[str(x) for x in row]
                       for row in r_port_compression.tolist()],
    "checks": checks,
    "passed": sum(bool(value) for value in checks.values()),
    "failed": sum(not bool(value) for value in checks.values()),
    "physical_gates_closed": [],
}

print(json.dumps(result, indent=2))
assert all(checks.values())
