"""Use the actual v252 constructor to identify all fixed free quarter channels.

These are exact controls of a GIVEN source, not derived physical masses.
The analytic all-word/R/time certificate is THERMAL-RIGIDITY.txt.
"""
from pathlib import Path
import importlib.util
import json
import sys

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
helper = ROOT / "experiments/theory-contracts/tfpt-unified-source-derivation-20261029/original_insertions/check_original_insertions.py"
spec = importlib.util.spec_from_file_location("existing_exact_constructor", helper)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
v, exact = old.v, old.exact


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run():
    yu, yd, ye = [np.diag(x).astype(complex) for x in
                  ([11, 13, 17], [19, 23, 29], [31, 37, 41])]
    y0, mr0 = np.diag([1, 2, 3]).astype(complex), np.diag([1, 4, 9]).astype(complex)
    g = exact(v.grading())
    plus = [i for i in range(v.NH) if g[i, i] == 1]
    eplus = sp.eye(v.NH)[:, plus]
    # Actual gauge neutral positive grade right antineutrino family module.
    neutral_ids = [v.NP + v.IDX[("nu", "R", 0, j)] for j in range(3)]
    eneutral = sp.eye(v.NH)[:, neutral_ids]
    rows = []
    for name, yn, mr, expected in [
        ("existing_diagonal_control", y0, mr0, 0),
        # Rank-deficient controls verify discovery of ALL channels and exclusion
        # when Majorana input restores the coupling. No physical input proposed.
        ("one_neutral_channel_control", np.diag([0, 2, 3]).astype(complex),
         np.diag([0, 4, 9]).astype(complex), 2),
        ("same_dirac_kernel_majorana_restored", np.diag([0, 2, 3]).astype(complex),
         mr0, 1),
    ]:
        phi = exact(v.full_D(yu, yd, ye, yn, mr))
        require(phi == phi.H, "source adjoint")
        require(phi*g + g*phi == sp.zeros(v.NH), "source oddness")
        kp = (phi*eplus).nullspace()
        kn = (phi*eneutral).nullspace()
        stacked = exact(np.vstack([mr, np.conjugate(yn)]))
        require(len(kn) == len(stacked.nullspace()), "neutral stacked criterion")
        require(len(kp) == expected, "complete positive channel dimension")
        images = [eplus*x for x in kp]
        require(all(phi*x == sp.zeros(v.NH, 1) and g*x == x for x in images),
                "all identified source channels")
        # One nonzero p intertwines. The theorem, not a numerical fit, then
        # supplies all p, H1 domains, both faces, products, reflection and time.
        p = sp.Rational(1, 4)
        h = p*g-sp.I*g*phi
        require(all(h*x == p*x for x in images), "same original quarter momentum")
        rows.append({"control": name, "positive_grade_dimension": len(plus),
                     "all_fixed_positive_channels": len(kp),
                     "gauge_neutral_channels": len(kn),
                     "neutral_stacked_rank": stacked.rank(),
                     "support_of_positive_kernel_basis": [
                         [i for i in range(v.NH) if x[i] != 0] for x in images]})
    # Nonreal check prevents replacing entrywise conjugation by adjoint.
    yn = np.array([[1, 1j, 0], [0, 0, 2], [0, 0, 0]], complex)
    mr = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 3]], complex)
    phi = exact(v.full_D(yu, yd, ye, yn, mr))
    direct = phi*eneutral
    stack = exact(np.vstack([mr, np.conjugate(yn)]))
    require(direct.nullspace() == stack.nullspace(), "nonreal neutral kernel")
    wrong = exact(np.vstack([mr, yn.conj().T]))
    require(stack.nullspace() != wrong.nullspace(), "adjoint substitution rejected")
    out = {"scope": "actual finite constructor, given constant source only",
           "arithmetic": "exact Gaussian integers and rational quarter p",
           "target_not_fitted": True, "controls": rows,
           "nonreal_entrywise_conjugation_control": True,
           "physical_source_origin_derived": False}
    Path(__file__).with_name("exact-source-channels.json").write_text(
        json.dumps(out, indent=2)+"\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    run()
