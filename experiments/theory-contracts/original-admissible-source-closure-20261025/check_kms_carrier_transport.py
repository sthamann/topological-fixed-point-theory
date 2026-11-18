"""Exact controls for KMS-AND-CARRIER-TRANSPORT.txt.

All finite covariances below are restrictions of the existing circle or of
its pure limit. They test boundary factors and block identities; the analytic
trace-class proof, not a finite truncation, proves the infinite determinant.
"""

from fractions import Fraction as F
import json
from pathlib import Path


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def power2(exponent):
    return F(2) ** exponent


def occupation(n, quarter):
    # beta=8 log(2), a=quarter/4. Both quarter sectors and full NS are exact.
    return 1 / (1 + power2(8 * n + 2 * quarter))


def window_ratio(radius, winding, z, quarter):
    value = F(1)
    for n in range(-radius, radius + 1):
        value *= (1 - z * occupation(n - winding, quarter))
        value /= (1 - z * occupation(n, quarter))
    return value


def two_det(a):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def a_matrix(c, z):
    return [[int(i == j) - z * c[i][j] for j in range(2)]
            for i in range(2)]


def rotate_diagonal(x, y):
    # An exact change of basis inside two existing circle modes.
    c, s = F(3, 5), F(4, 5)
    return [[c*c*x+s*s*y, c*s*(x-y)],
            [c*s*(x-y), s*s*x+c*c*y]]


def schur(a):
    need(a[1][1] != 0, "QQ block must be invertible")
    return a[0][0] - a[0][1] * a[1][0] / a[1][1]


records = []
for quarter in (1, 2, 3):
    for z in (F(-1, 2), F(0), F(1, 3), F(3, 4)):
        errors = []
        for radius in (1, 2, 5, 10):
            value = window_ratio(radius, 1, z, quarter)
            left = occupation(-radius-1, quarter)
            right = occupation(radius, quarter)
            boundary = (1-z*left)/(1-z*right)
            need(value == boundary, "actual open-window boundary quotient")
            defect_trace = left-right
            need(0 < defect_trace < 1, "finite carrier misses thermal tail")
            error = abs(value-(1-z))
            bound = abs(z)*((1-left)+(1+abs(z))*right)/(1-abs(z))
            need(error <= bound, "explicit Fredholm tail error bound")
            errors.append(error)
        if z:
            need(all(a > b for a, b in zip(errors, errors[1:])),
                 "finite boundary error decreases without cyclic closure")
        records.append({"quarter": quarter, "z": str(z),
                        "radius10_error": str(errors[-1])})

# Cover versus base winding: retain the boundary modes of the actual shift.
cover_checks = 0
radius = 5
for winding in (-2, -1, 1, 2):
    for quarter in (1, 2, 3):
        z = F(1, 3)
        direct = window_ratio(radius, winding, z, quarter)
        if winding > 0:
            boundary = F(1)
            for n in range(-radius-winding, -radius):
                boundary *= 1-z*occupation(n, quarter)
            for n in range(radius-winding+1, radius+1):
                boundary /= 1-z*occupation(n, quarter)
        else:
            boundary = F(1)
            for n in range(radius+1, radius-winding+1):
                boundary *= 1-z*occupation(n, quarter)
            for n in range(-radius, -radius-winding):
                boundary /= 1-z*occupation(n, quarter)
        need(direct == boundary, "index-m winding boundary identity")
        cover_checks += 1

# A finite CYCLIC unitary shift has the compensating edge and total det one.
values = [occupation(n, 1) for n in range(-2, 3)]
cyclic = F(1)
z = F(1, 3)
for i, value in enumerate(values):
    cyclic *= (1-z*values[(i-1) % len(values)])/(1-z*value)
need(cyclic == 1, "finite cyclic shift cannot carry the infinite index")

# The actual pure two-face Cauchy graph has rank one in every Fourier block.
graph_blocks = 0
for n in range(-4, 5):
    for quarter in (1, 2, 3):
        t = power2(-4*n-quarter)
        denom = 1+t*t
        g = [[1/denom, t/denom], [t/denom, t*t/denom]]
        for z in (F(1, 3), F(-1, 2)):
            need(two_det(a_matrix(g, z)) == 1-z,
                 "full Cauchy graph block determinant")
        graph_blocks += 1

# Nonreducing compression: check the complement, not only the diagonal entry.
schur_records = []
for quarter in (1, 2, 3):
    old = rotate_diagonal(occupation(-1, quarter), occupation(0, quarter))
    new = rotate_diagonal(occupation(-2, quarter), occupation(-1, quarter))
    for z in (F(1, 3), F(3, 4)):
        a, a_new = a_matrix(old, z), a_matrix(new, z)
        full = two_det(a_new)/two_det(a)
        carrier = schur(a_new)/schur(a)
        complement = a_new[1][1]/a[1][1]
        naive = a_new[0][0]/a[0][0]
        need(full == carrier*complement, "exact Schur relative determinant")
        need(naive != full and naive != carrier,
             "naive covariance compression loses the actual Schur term")
        schur_records.append({"quarter": quarter, "z": str(z),
                              "full": str(full), "naive": str(naive)})

# Pure circle restriction: retaining only part of the actual defect is not
# the same as retaining its normalized rank-one projector.
old = rotate_diagonal(F(1), F(0))
new = rotate_diagonal(F(1), F(1))
for z in (F(1, 3), F(3, 4)):
    a, a_new = a_matrix(old, z), a_matrix(new, z)
    w = F(16, 25)
    naive = a_new[0][0]/a[0][0]
    expected = 1-z*w/(1-z*F(9, 25))
    need(naive == expected and naive != 1-z,
         "rank-one weight and compressed resolvent are retained")
    need(two_det(a_new)/two_det(a) == 1-z,
         "full pure circle restriction retains the entire defect")
    need((schur(a_new)/schur(a))*(a_new[1][1]/a[1][1]) == 1-z,
         "carrier plus complement restores the pure factor")

# Reference alpha=0 is not a seam-free reference q=0.
delta = F(1, 16)
q_alpha = delta/2
relative = (1-q_alpha)/(1-delta)
need(relative == F(31, 30) and relative != 1-q_alpha,
     "EM reference normalization must retain q(alpha_*)")

result = {
    "verdict": "PASS",
    "analytic_proof": "KMS-AND-CARRIER-TRANSPORT.txt",
    "open_window_exact_boundary_controls": records,
    "cover_and_reverse_winding_controls": cover_checks,
    "finite_cyclic_compensation": "PASS",
    "full_cauchy_graph_blocks_checked": graph_blocks,
    "nonreducing_actual_two_mode_schur_controls": schur_records,
    "pure_partial_retention_control": "PASS",
    "EM_reference_normalization_control": "PASS",
    "active_v258_forward_source_map": "NOT_CONSTRUCTED_BY_V258",
    "physical_EM_measure_identity": "OPEN",
    "complete_TFPT_solution": False,
}
payload = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)+"\n"
Path(__file__).with_name("kms-carrier-transport-results.json").write_text(payload)
print(payload, end="")
