"""Small exact gate on Claude7: Lie-bracket tensor is not yet a Weyl Yukawa.

Class tested: one skew family form and arbitrary polynomials of the SAME
anchor/triplet mark on its two legs. No full TFPT exclusion is asserted.
"""
from pathlib import Path
import json
import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sym(a):
    return (a+a.T)/2


def run():
    x = sp.diag(3, -1, -1, -1)
    require(x*x == 2*x+3*sp.eye(4), "anchor polynomial closure")
    skew = []
    for i in range(4):
        for j in range(i+1, 4):
            b = sp.zeros(4)
            b[i, j], b[j, i] = 1, -1
            skew.append(b)
    eps = sp.Matrix([[0, 1], [-1, 0]])
    require(all(sym(b) == sp.zeros(4) for b in skew), "bare symmetric part")
    # Grassmann quadratic coefficients use the antisymmetric TOTAL matrix.
    # With skew internal b and skew Lorentz epsilon, b tensor epsilon is
    # symmetric and its Grassmann quadratic is identically zero.
    require(all(sp.kronecker_product(b, eps) ==
                sp.kronecker_product(b, eps).T for b in skew), "Weyl statistics")
    candidates = [sym(l*b*r) for b in skew
                  for l in (sp.eye(4), x) for r in (sp.eye(4), x)]
    require(all(a[1:, 1:] == sp.zeros(3) for a in candidates),
            "whole triplet mass block vanishes")
    span_dim = sp.Matrix.hstack(*(a.reshape(16, 1) for a in candidates)).rank()
    require(span_dim == 3, "complete anchor-only symmetric image")
    # Entire symmetric image is the star form [[0,v^T],[v,0_3]]. For complex v,
    # AA* has eigenvalues ||v||^2 twice and zero twice; rank at most two.
    v1, v2, v3 = sp.symbols("v1 v2 v3", real=True)
    a = sp.Matrix([[0, v1, v2, v3], [v1, 0, 0, 0],
                   [v2, 0, 0, 0], [v3, 0, 0, 0]])
    t = sp.symbols("t")
    want = t*t*(t-(v1*v1+v2*v2+v3*v3))**2
    require(sp.expand((a*a.T).charpoly(t).as_expr()-want) == 0,
            "real star singular-square certificate")
    # A genuinely distinguishing family generator changes the candidate class.
    z = sp.diag(3, 1, -1, -3)
    broken = [sym(b*z) for b in skew]
    broken_dim = sp.Matrix.hstack(*(b.reshape(16, 1) for b in broken)).rank()
    require(broken_dim == 6, "distinct-family generator enlarges image")
    out = {"scope": "one Lie-bracket family form, same anchor mark only",
           "bare_Weyl_Yukawa_zero": True,
           "anchor_polynomial_basis": ["I", "X"],
           "skew_basis_size": len(skew), "all_reduced_words": len(candidates),
           "symmetric_image_dimension": span_dim,
           "candidate_symmetrized_triplet_block_zero": True,
           "single_4_by_4_family_channel_rank_bound": 2,
           "distinct_family_control_image_dimension": broken_dim,
           "other_original_family_operations_excluded": False,
           "E8_source_is_physical_4D_Yukawa_proved": False}
    Path(__file__).with_name("e8-anchor-yukawa.json").write_text(
        json.dumps(out, indent=2)+"\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    run()
