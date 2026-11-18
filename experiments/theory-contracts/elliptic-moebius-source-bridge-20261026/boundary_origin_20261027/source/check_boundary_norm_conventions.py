#!/usr/bin/env python3
"""Exact exponent and beta-character audit for the proposed q-boundary insertion."""
from dataclasses import dataclass
import json
from pathlib import Path

@dataclass(frozen=True)
class Line:
    theta: int
    nq: int
    np: int
    beta_exp: int  # i**beta_exp

    def __mul__(self, other):
        return Line(self.theta + other.theta, self.nq + other.nq,
                    self.np + other.np, (self.beta_exp + other.beta_exp) % 4)

    def inv(self):
        return Line(-self.theta, -self.nq, -self.np, (-self.beta_exp) % 4)

    def as_dict(self):
        return {"theta_degree": self.theta, "Nq_power": self.nq,
                "Np_power": self.np, "beta_phase": ["1", "i", "-1", "-i"][self.beta_exp]}

ONE = Line(0, 0, 0, 0)
THETA = Line(1, 0, 0, 0)          # D(E)^(-1)
NQ = Line(0, 1, 0, 3)             # beta on Nq is -i
NP = Line(0, 0, 1, 1)             # beta on Np is +i

cases = {
    "target_q_inverse": THETA * THETA * NQ.inv(),
    "mixed_q_vanishing": THETA * (THETA * NQ),
    "both_q_vanishing": (THETA * NQ) * (THETA * NQ),
    "mixed_p_vanishing": THETA * (THETA * NP),
}

expected = {
    "target_q_inverse": {"theta_degree": 2, "Nq_power": -1, "Np_power": 0, "beta_phase": "i"},
    "mixed_q_vanishing": {"theta_degree": 2, "Nq_power": 1, "Np_power": 0, "beta_phase": "-i"},
    "both_q_vanishing": {"theta_degree": 2, "Nq_power": 2, "Np_power": 0, "beta_phase": "-1"},
    "mixed_p_vanishing": {"theta_degree": 2, "Nq_power": 0, "Np_power": 1, "beta_phase": "i"},
}

actual = {name: line.as_dict() for name, line in cases.items()}
errors = []
for name, want in expected.items():
    if actual.get(name) != want:
        errors.append({"case": name, "expected": want, "actual": actual.get(name)})

# Logical consequences checked without assert, so python -OO is identical.
if actual["mixed_q_vanishing"] == actual["target_q_inverse"]:
    errors.append({"case": "separation", "error": "q-vanishing was incorrectly identified with target"})
if actual["mixed_p_vanishing"]["beta_phase"] != actual["target_q_inverse"]["beta_phase"]:
    errors.append({"case": "p_phase", "error": "p direct phase should match target phase"})
if actual["mixed_p_vanishing"]["Np_power"] != 1 or actual["mixed_p_vanishing"]["Nq_power"] != 0:
    errors.append({"case": "p_type", "error": "phase match must not erase q/p line type"})

payload = {
    "conventions": {
        "D(E)": "det H0(E) tensor det H1(E)^(-1)",
        "O(Theta)": "D(E)^(-1)",
        "D(E(-Dq))": "D(E) tensor Nq^(-1)",
        "D(E(-Dp))": "D(E) tensor Np^(-1)",
    },
    "cases": actual,
    "checks_passed": not errors,
    "errors": errors,
}
out = Path(__file__).with_name("boundary-norm-conventions.json")
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(json.dumps(payload, indent=2, sort_keys=True))
raise SystemExit(0 if not errors else 1)
