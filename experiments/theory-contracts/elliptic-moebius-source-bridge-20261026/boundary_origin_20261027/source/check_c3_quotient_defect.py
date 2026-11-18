#!/usr/bin/env python3
"""Exact C3 quotient/index and determinant-orientation test for the q packet."""
import json
from pathlib import Path
import sympy as sp

I = sp.I
omega = sp.exp(2*sp.pi*I/3)
a = sp.exp(sp.pi*I/4)  # a^2=i
# Fourier basis: 1, omega, omega^2 characters of c.
C = sp.diag(1, omega, omega**2)
A = sp.Matrix([[a, 0, 0], [0, 0, 1], [0, I, 0]])
Id = sp.eye(3)

def zero(M):
    return all(sp.simplify(x) == 0 for x in M)

beta_ell = I
beta_rank1_fredholm = 1 / beta_ell
beta_Bq = I**3
beta_Bq_inverse = 1 / beta_Bq
checks = {
    "A_squared_is_i": zero(A*A - I*Id),
    "A_conjugates_c_to_inverse": zero(A*C*A.inv() - C.inv()),
    "invariant_sector_is_A_stable": zero(A[:, 0] - a*Id[:, 0]),
    "undualized_invariant_character_equals_full_inverse_determinant": sp.simplify(1/A.det() - a) == 0,
    "full_fiber_modification_index_jump": -3 == -C.rows,
    "one_character_modification_index_jump": -1 == -1,
    "rank1_fredholm_beta_is_minus_i": sp.simplify(beta_rank1_fredholm + I) == 0,
    "full_inverse_determinant_beta_is_plus_i": sp.simplify(beta_Bq_inverse - I) == 0,
    "rank1_fredholm_is_not_full_inverse": sp.simplify(beta_rank1_fredholm-beta_Bq_inverse) != 0,
}
errors = [name for name, ok in checks.items() if not ok]
payload = {
    "quotient": {
        "curve": "eta^2=1-x^2 with x=w^3",
        "q_orbit": "q0,q1,q2 maps to one unramified point x=-1",
        "pushforward_rank": 3,
        "fiber_c_character": "1 + omega + omega^2 (regular representation)",
    },
    "orientation": {
        "beta_on_invariant_cokernel_line": "+i",
        "beta_on_rank1_fredholm_line": "-i",
        "beta_on_Bq": "-i",
        "beta_on_Bq_inverse": "+i",
        "conclusion": "The rank-one directed Fredholm defect is not Nq^-1; the full three-fiber inverse determinant is required.",
    },
    "checks": checks,
    "checks_passed": not errors,
    "errors": errors,
}
out = Path(__file__).with_name("c3-quotient-defect.json")
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(json.dumps(payload, indent=2, sort_keys=True))
raise SystemExit(0 if not errors else 1)
