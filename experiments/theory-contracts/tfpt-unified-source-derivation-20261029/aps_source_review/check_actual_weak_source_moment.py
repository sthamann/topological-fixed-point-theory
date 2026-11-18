"""Exact source moment from the actual v252 weak lift and finite Dirac.

Diagnostic coefficients test orientation sensitivity, not mass selection.
No original source, verification suite or physical parameters are changed.
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "verification"))
import v252_full_finite_triple as v


def exact(a):
    return s.MutableSparseMatrix(a.shape[0], a.shape[1], {
        (i, j): s.Rational(str(a[i, j].real)) + s.I*s.Rational(str(a[i, j].imag))
        for i, j in zip(*np.nonzero(a))
    })


def main():
    z = np.zeros((3, 3), complex)
    sx = np.array([[0, 1], [1, 0]], complex)
    # Real gauge lift acts on the antihermitian generator before multiplying
    # by -i. Antilinearity therefore gives MINUS on the antiparticle Q block.
    T = v.rep(0, 1j*sx, z)
    K = v.Kswap()
    Q = exact(-1j*(T + K@T.conj()@K))
    gamma = exact(v.grading())
    yu = np.diag([1, 2, 3]).astype(complex)
    yd = np.diag([4, 5, 6]).astype(complex)
    perm = np.eye(3)[[2, 1, 0]]
    # Same, nonzero lepton/Majorana data in both controls; these do not
    # interact with quark flavour in this original weak representation.
    ye, yn = np.diag([2, 3, 4]), np.diag([1, 3, 5])
    mr = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])
    outputs = []
    checks = {"Q_selfadjoint": Q == Q.H,
              "Q_even": Q*gamma == gamma*Q}
    for label, down in (("aligned", yd), ("left_rotated_down", perm@yd)):
        D = exact(v.full_D(yu, down, ye, yn, mr))
        checks[label+"_D_odd"] = gamma*D + D*gamma == s.zeros(96)
        H0 = gamma - 4*s.I*gamma*D
        H1 = 4*Q*gamma
        checks[label+"_h_selfadjoint"] = H0 == H0.H and H1 == H1.H
        # H0+tH1 = 4h at the existing quarter mode p=1/4.
        poly = [s.eye(96)]
        moments = {}
        for n in range(1, 7):
            new = [s.zeros(96) for _ in range(len(poly)+1)]
            for k, coefficient in enumerate(poly):
                new[k] += coefficient*H0
                new[k+1] += coefficient*H1
            poly = new
            if n in (2, 4, 6):
                moments[n] = [s.expand(s.trace(c)/4**n) for c in poly]
        outputs.append(moments)
    differences = {n: [s.expand(a-b) for a, b in zip(outputs[0][n], outputs[1][n])]
                   for n in (2, 4, 6)}
    checks["second_and_fourth_moments_agree"] = all(x == 0 for n in (2, 4) for x in differences[n])
    checks["sixth_moment_difference_exact"] = differences[6] == [0, 0, 5760, 0, 0, 0, 0]
    checks["same_down_singular_data"] = exact(yd).H*exact(yd) == exact(perm@yd).H*exact(perm@yd)
    result = {
        "verdict": "ACTUAL_GAUGE_COUPLED_SOURCE_RETAINS_RELATIVE_FAMILY_ORIENTATION",
        "operator": "h(t)=(I/4+tQ)gamma-i gamma D_F, Q=-i(T+J_F T J_F^-1), T=rep(0,i sigma_x,0)",
        "moments_t_coefficients": [{str(n): [str(c) for c in rows] for n, rows in out.items()} for out in outputs],
        "difference_t_coefficients": {str(n): [str(c) for c in rows] for n, rows in differences.items()},
        "cross_face_trace_difference_leading_ell6": "-244*t**2",
        "cross_face_logdet_difference_leading_ell6": "-128*t**2",
        "checks": {k: bool(vv) for k, vv in checks.items()},
        "all_checks_pass": all(checks.values()),
        "scope": "Finite p=1/4 block of the existing source. Arbitrary formal weak field t; fixed diagnostic Yukawas and nonzero fixed MR. No selected gauge background, Yukawa variation, sum over circle modes, or complete Gamma_ren coefficient is asserted."
    }
    raw = json.dumps(result, indent=2, sort_keys=True)+"\n"
    (HERE/"actual-weak-source-moment-results.json").write_text(raw)
    print(raw, end="")
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
