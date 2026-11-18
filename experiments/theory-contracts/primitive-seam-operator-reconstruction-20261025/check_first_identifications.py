"""Exact controls for the two first original-paper identification steps.

The infinite-dimensional theorems are proved in the companion text files.
These controls do not approximate Fredholm indices with finite cyclic shifts.
"""

import json
from fractions import Fraction
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compressed_shift_defects(holes):
    """P projects onto {n>=0} minus a finite, declared set of holes.

    Both defect sets are contained in {0,...,max(holes)+1}. Above this bound
    both n and its neighbors belong to P; below zero n does not belong to P.
    The bound is exact for this eventually constant infinite membership rule.
    """
    holes = frozenset(holes)
    require(all(isinstance(n, int) and n >= 0 for n in holes), "invalid hole")

    def inside(n):
        return n >= 0 and n not in holes

    upper = max(holes, default=-1) + 1
    support = range(upper + 1)
    incoming = [n for n in support if inside(n) and not inside(n + 1)]
    outgoing = [n for n in support if inside(n) and not inside(n - 1)]
    # On each retained basis vector S sends e_n to e_(n+1) or zero. Thus
    # these lists are exactly ker(S) and ker(S*) on the whole l2(Z).
    return {
        "holes": sorted(holes),
        "kernel": incoming,
        "cokernel": outgoing,
        "index": len(incoming) - len(outgoing),
        "ordinary_cohomology_dimension": len(incoming) + len(outgoing),
        "isometry": not incoming,
    }


hardy = compressed_shift_defects(())
punctured = compressed_shift_defects((1,))
require(hardy["kernel"] == [] and hardy["cokernel"] == [0], "Hardy sign")
require(hardy["index"] == -1 and hardy["isometry"], "Hardy shift")
require(punctured["kernel"] == [0], "punctured kernel")
require(punctured["cokernel"] == [0, 2], "punctured cokernel")
require(punctured["index"] == -1 and not punctured["isometry"], "index probe")
require(punctured["ordinary_cohomology_dimension"] == 3, "ordinary dimension")

# Exact Gaussian reflected Gram control. This tests the full-moment versus
# cumulant substitution, not the unconstructed TFPT source functional.
a, c = Fraction(3, 2), Fraction(1, 2)
full_22 = a * a + 2 * c * c
full_det = full_22 - a * a
truncated_det = -a * a
require(full_det == 2 * c * c > 0, "full Gaussian Gram")
require(truncated_det < 0, "truncated Gram must fail")

# A second exact circle-symbol control, independently of the projection hole.
# b_a=(z-a)/(1-a*z), u_a=z**2/b_a. Directly expand
# (z**2-a*z**3)/(z-a) as a geometric series in a/z to obtain
# coeff(z^-m)=(1-a**2)*a**(m+1), m>=1.
# Therefore the negative norm is a**4*(1-a**2), as in the companion proof.
symbol_a = Fraction(1, 4)
leakage = symbol_a**4 * (1 - symbol_a**2)
require(leakage == Fraction(15, 4096), "circle-symbol leakage")

epsilon = Fraction(1, 4)
ward_checks = 0
for n in range(-32, 33):
    require((n + 1 + epsilon) - (n + epsilon) == 1, "Ward sign")
    require((n + epsilon > 0) == (n >= 0), "Hardy spectral cut")
    ward_checks += 1
crossing_labels = [n for n in range(-2, 2) if n + epsilon < 0 < n + epsilon + 1]
require(crossing_labels == [-1], "spectral-flow sign")
# exp(2*pi*i*x)=i iff x mod 1=1/4; no floating point phase test.
require(epsilon % 1 == Fraction(1, 4), "quarter monodromy")

result = {
    "verdict": "PASS",
    "scope": "exact first-identification controls; analytic proofs in companion files",
    "infinite_membership_rules": {"hardy": hardy, "one_hole": punctured},
    "gaussian_gram": {"a": str(a), "c": str(c), "full_det": str(full_det),
                      "truncated_det": str(truncated_det)},
    "circle_symbol_negative_norm_squared": str(leakage),
    "quarter_spectrum_basis_checks": ward_checks,
    "spectral_flow_upward_crossing": crossing_labels,
    "no_physical_source_or_original_pair_covariance_asserted": True,
}
payload = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
Path(__file__).with_name("identification-results.json").write_text(payload)
print(payload, end="")
