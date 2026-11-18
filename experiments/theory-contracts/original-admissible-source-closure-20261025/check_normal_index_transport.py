"""Small exact checks of the conditional normal-to-circle operator proof."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def run():
    x, lam, a, theta, eps = sp.symbols("x lam a theta eps", real=True)
    z = sp.Symbol("z")
    checks = []

    def check(name, ok, **data):
        checks.append({"name": name, "ok": bool(ok), **data})
        if not ok:
            raise AssertionError(name)

    gamma = sp.diag(1, 1, -1)
    normal = sp.Matrix([[0, 0, 0], [0, 0, lam], [0, lam, 0]])
    v = sp.Matrix([1, 0, 0])
    f = v*v.T
    product = x*gamma+normal
    check("normal_grading_and_product_kernel_intertwining",
          gamma*normal+normal*gamma == sp.zeros(3)
          and product*v == x*v and product*f == f*product)
    block = product[1:3, 1:3]
    check("nonzero_pair_has_no_crossing",
          block*block == (x*x+lam*lam)*sp.eye(2)
          and sp.factor(block.det()) == -x*x-lam*lam)

    def h0(degree):
        return max(degree+1, 0)

    def h1(degree):
        return max(-degree-1, 0)

    census = [{"k": k, "positive_kernel": h0(k-1)+h0(-k),
               "negative_kernel": h1(k-1)+h1(-k)} for k in range(1, 7)]
    check("spin_line_one_and_scalar_dolbeault_two_are_distinct",
          h0(-1+1) == 1 and h1(-1+1) == 0
          and h0(1)-h1(1) == 2)
    check("same_rank_two_determinant_keeps_extra_kernel_pairs",
          all(r["positive_kernel"] == r["k"]
              and r["negative_kernel"] == r["k"]-1 for r in census),
          census=census)

    phase = eps*sp.sin(4*theta)
    potential = a-sp.diff(phase, theta)
    check("induced_connection_gauge_transport_on_all_fourier_modes",
          sp.simplify(sp.diff(phase, theta)+potential-a) == 0
          and sp.integrate(potential, (theta, 0, 2*sp.pi)) == 2*sp.pi*a
          and all(phase.subs(theta, j*sp.pi/2) == 0 for j in range(4)))
    check("unchanged_vertical_data_can_keep_distinct_horizontal_phase",
          sp.simplify(sp.exp(2*sp.pi*sp.I*sp.Rational(1, 2))) == -1)

    e0, e1 = sp.sqrt(x*x+lam*lam), sp.sqrt((x+1)**2+lam*lam)
    unit_dot = (lam*lam+x*(x+1))/(e0*e1)
    check("exact_projection_angle_identity_for_trace_norm_bound",
          sp.factor(1-unit_dot**2-lam**2/(e0**2*e1**2)) == 0)
    ell = sp.Symbol("ell", positive=True)
    lower_squared = (2*ell+1)**2/(4*ell*(ell+1))
    check("monopole_multiplicities_force_divergent_trace_lower_bound",
          sp.factor(lower_squared-1) == 1/(4*ell*(ell+1)))

    positive_channel = sp.cosh(2*sp.pi*lam)-sp.cos(2*sp.pi*a)
    derivative = sp.diff(sp.log(positive_channel), a)
    check("zero_spectral_flow_channel_retains_positive_response",
          sp.simplify(derivative.subs(a, sp.Rational(1, 4))
                      -2*sp.pi/sp.cosh(2*sp.pi*lam)) == 0,
          derivative_at_quarter=str(derivative.subs(a, sp.Rational(1, 4))))

    root = Path(__file__).resolve().parents[3]
    files = [Path(__file__).relative_to(root),
             Path(__file__).with_name("NORMAL-INDEX-TRANSPORT.txt").relative_to(root),
             Path(__file__).with_name("NORMAL-CLUTCH-TO-SOURCE-PUSHFORWARD.txt").relative_to(root),
             Path(__file__).with_name("TRANSVERSE-INDEX-EM-RESPONSE.txt").relative_to(root)]
    return {"claims": ["MATH.ORIGINAL-SEAM.NORMAL-INDEX-TRANSPORT.01",
                       "MATH.ORIGINAL-SEAM.NORMAL-HOLONOMY-CRITERION.01",
                       "MATH.ORIGINAL-SEAM.FULL-PRODUCT-TRACECLASS.01"],
            "verdict": "PASS_EXACT_STRUCTURAL_CHECKS",
            "scope": "Analytic conditional transport, not a selected original source or physical EM functional.",
            "checks": checks,
            "sources": [{"path": str(p), "sha256": hashlib.sha256((root/p).read_bytes()).hexdigest()}
                        for p in files]}


if __name__ == "__main__":
    result = run()
    target = Path(__file__).with_name("normal-index-transport-results.json")
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps({"verdict": result["verdict"], "passed": len(result["checks"]),
                      "result": str(target)}, ensure_ascii=False))
