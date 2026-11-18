"""Exact algebra checks for the infinite-circle pair-index proof.

All finite supports below sit inside the original bilateral circle. There is
no periodic wrap identification at the cutoff.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def run() -> dict:
    z = sp.Symbol("z")
    checks = []

    def check(name, ok, **detail):
        checks.append({"name": name, "ok": bool(ok), **detail})
        if not ok:
            raise AssertionError(name)

    c, s = sp.Rational(3, 5), sp.Rational(4, 5)
    rotation = sp.Matrix([[c, -s], [s, c]])
    p = sp.diag(0, 1, 0)
    q = sp.diag(1, 0, 0)
    q[1:3, 1:3] = rotation * sp.diag(1, 0) * rotation.T
    identity = sp.eye(3)
    relative = (identity-z*q)*(identity-z*p).inv()
    check("noncommuting_pair_is_orthogonal_projections",
          p*p == p and q*q == q and p.T == p and q.T == q
          and p*q != q*p)
    check("indefinite_rank_three_index_one",
          (q-p).rank() == 3 and sp.trace(q-p) == 1
          and set((q-p).eigenvals()) == {sp.Integer(1), s, -s},
          eigenvalues=[str(x) for x in (q-p).eigenvals()])
    check("pair_relative_determinant_is_one_minus_z",
          sp.factor(relative.det()) == 1-z,
          determinant=str(sp.factor(relative.det())))
    derivative = relative.diff(z)*relative.inv()
    expected = (identity-z*q)*(p-q)*(identity-z*q).inv()/(1-z)
    check("exact_operator_log_derivative",
          (derivative-expected).applyfunc(sp.factor) == sp.zeros(3))

    # P is the negative sea on the whole bilateral circle. Only this window
    # supports C-P and UCU*-C. The shift takes -1->0 and 1->2 without wrap.
    modes = [-2, -1, 0, 1, 2]
    idx = {mode: j for j, mode in enumerate(modes)}
    sea = sp.diag(*[int(mode < 0) for mode in modes])
    cov = sea.copy()
    block = rotation*sp.diag(sp.Rational(2, 3), sp.Rational(1, 4))*rotation.T
    for i, ni in enumerate([-1, 1]):
        for j, nj in enumerate([-1, 1]):
            cov[idx[ni], idx[nj]] = block[i, j]
    delta = cov-sea
    shifted_delta = sp.zeros(5)
    for i, ni in enumerate(modes):
        for j, nj in enumerate(modes):
            if delta[i, j] != 0:
                shifted_delta[idx[ni+1], idx[nj+1]] = delta[i, j]
    mix = sp.eye(5)
    for i, ni in enumerate([-2, 2]):
        for j, nj in enumerate([-2, 2]):
            mix[idx[ni], idx[nj]] = rotation[i, j]
    pi0 = sp.zeros(5)
    pi0[idx[0], idx[0]] = 1
    shifted_cov = mix*(sea+pi0+shifted_delta)*mix.T
    # Outside this window the covariances equal the same sea, so this is
    # precisely the full relative Fredholm determinant of their difference.
    ratio = sp.factor((sp.eye(5)-z*shifted_cov).det()
                      / (sp.eye(5)-z*cov).det())
    check("nonprojective_covariance_is_positive_contraction",
          cov*cov != cov
          and all(0 <= value <= 1 for value in cov.eigenvals())
          and all(0 <= value <= 1 for value in shifted_cov.eigenvals()))
    check("same_bilateral_source_nonprojective_factor", ratio == 1-z,
          determinant=str(ratio),
          shifted_eigenvalues={str(k): v for k, v in shifted_cov.eigenvals().items()})

    weak = sp.diag(z, 1)
    scalar_weak = z*sp.eye(2)
    full_special = sp.diag(1/z, 1, 1, z, 1)
    check("carrier_clutch_degree_is_not_rank",
          weak.det() == z and scalar_weak.det() == z*z
          and full_special.det() == 1,
          weak_determinant=str(weak.det()),
          scalar_determinant=str(scalar_weak.det()),
          full_special_determinant=str(full_special.det()))

    root = Path(__file__).resolve().parents[3]
    source_paths = [
        Path(__file__).relative_to(root),
        Path(__file__).with_name("PAIR-INDEX-DETERMINANT.txt").relative_to(root),
        Path("documents/archive--paper-latex--tfpt-42/sections/005-subsection-internal-carrier-stabilizer-and-one-family-pa.tex"),
        Path("documents/archive--paper-latex--tfpt-42/sections/007-subsection-carrier-form-of-the-electromagnetic-closure.tex"),
        Path("documents/archive--paper-latex--tfpt-42/sections/012-subsection-bosonic-relative-index-and-higgs-selection.tex"),
        Path("documents/archive--paper-latex--tfpt-31/sections/007-subsection-occupancy-and-boundary-to-bulk-bridge.tex"),
        Path("documents/archive--paper-latex--tfpt-32/sections/007-subsection-occupancy-and-boundary-to-bulk-bridge.tex"),
    ]
    sources = [{"path": str(path),
                "sha256": hashlib.sha256((root/path).read_bytes()).hexdigest()}
               for path in source_paths]
    return {"claim": "MATH.ORIGINAL-SEAM.PAIR-INDEX.01",
            "verdict": "PASS_EXACT_ALGEBRA",
            "scope": "Infinite analytic theorem plus exact finite-support checks; no physical source selection or EM Ward proof.",
            "checks": checks, "sources": sources}


if __name__ == "__main__":
    result = run()
    target = Path(__file__).with_name("pair-index-results.json")
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps({"verdict": result["verdict"],
                      "passed": len(result["checks"]),
                      "result": str(target)}, ensure_ascii=False))
