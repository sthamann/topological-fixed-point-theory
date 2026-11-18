"""Targeted audit of the finite graph particle-hole block in Claude48."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "verification"))

from verification import v252_full_finite_triple as v252  # noqa: E402


a = sp.symbols("a", positive=True)
c = sp.sqrt(a)
d = sp.sqrt(1 - a)
V = sp.Matrix([c, d])
Vh = sp.Matrix([-d, c])
Swap = sp.Matrix([[0, 1], [1, 0]])

face = {
    "pp": sp.simplify((V.T * Swap * V)[0]),
    "hp": sp.simplify((Vh.T * Swap * V)[0]),
    "ph": sp.simplify((V.T * Swap * Vh)[0]),
    "hh": sp.simplify((Vh.T * Swap * Vh)[0]),
}
assert face == {
    "pp": 2 * sp.sqrt(a) * sp.sqrt(1 - a),
    "hp": 2 * a - 1,
    "ph": 2 * a - 1,
    "hh": -2 * sp.sqrt(a) * sp.sqrt(1 - a),
}

ell = sp.symbols("ell", positive=True)
aell = 1 / (1 + sp.exp(-ell))
def exp_simplify(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.together(expr.rewrite(sp.exp)))

assert exp_simplify(face["pp"].subs(a, aell) ** 2
                    - 1 / sp.cosh(ell / 2) ** 2) == 0
assert exp_simplify(face["ph"].subs(a, aell) - sp.tanh(ell / 2)) == 0
assert exp_simplify((face["ph"].subs(a, aell) / face["pp"].subs(a, aell)) ** 2
                    - sp.sinh(ell / 2) ** 2) == 0

u, v = sp.symbols("I02 I21", nonzero=True)
ii = sp.I
def E(r: int, s: int) -> sp.Matrix:
    out = sp.zeros(3)
    out[r, s] = 1
    return out

mF = u * E(0, 2) + v * E(2, 1)
mN = -ii * (u * E(1, 2) + v * E(2, 0))
traces = {
    "FF": sp.simplify(sp.trace(mF * mF)),
    "NN": sp.simplify(sp.trace(mN * mN)),
    "FN": sp.simplify(sp.trace(mF * mN)),
    "NF": sp.simplify(sp.trace(mN * mF)),
}
assert traces == {"FF": 0, "NN": 0, "FN": -2 * ii * u * v,
                  "NF": -2 * ii * u * v}

gamma = v252.grading()
rank_plus_full_96 = int(((gamma.diagonal() + 1) / 2).sum())
assert gamma.shape == (96, 96)
assert rank_plus_full_96 == 48
assert v252.NG == 3
rank_plus_per_family = rank_plus_full_96 // v252.NG
assert rank_plus_per_family == 16

out = {
    "scope": (
        "Finite Pi/Pi_h contribution for the existing Phi=0 NS graph source. "
        "No selection of this compression as the whole original kernel."
    ),
    "face_blocks": {k: str(v) for k, v in face.items()},
    "at_half_mode": {
        "a": "1/(1+exp(-ell_NS))",
        "pp": "sech(ell_NS/2)",
        "ph": "hp=tanh(ell_NS/2)",
        "ratio_ph_to_pp": "sinh(ell_NS/2)",
    },
    "family_traces": {k: str(v) for k, v in traces.items()},
    "actual_v252": {
        "dimension": int(gamma.shape[0]),
        "families": int(v252.NG),
        "rank_P_gamma_plus_full": rank_plus_full_96,
        "rank_P_gamma_plus_per_family": rank_plus_per_family,
    },
    "normalization": {
        "per_Pgamma_component": "(2/pi)^2 tanh(ell_NS/2)^2 Tr_family(m_mu m_nu)",
        "ordinary_full_finite_trace": "16 times the per-component expression",
        "note": "A selfdual second-quantization convention may add its stated global 1/2, but cannot remove the carrier multiplicity without an explicit normalized trace.",
    },
    "passed": True,
}
target = HERE / "claude48-pauli-ph-face-block-results.json"
target.write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
