"""Exact controls for recovered original circle and full heat determinant.

Infinite-dimensional statements are proved in HEAT-TRANSFER-IR-CLOSURE.txt.
No chosen physical source, parameter identification or finite Fredholm
approximation is asserted by these algebraic controls.
"""

import json
from fractions import Fraction as F
from pathlib import Path


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def mul(left, right, degree):
    out = [F(0) for _ in range(degree + 1)]
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if i + j <= degree:
                out[i + j] += x * y
    return out


degree = 20
b = F(5, 4)
direct = [F(1)] + [F(0)] * degree
for odd in range(1, degree + 1, 2):
    factor = [F(0)] * (degree + 1)
    coefficient = F(1)
    factor[0] = coefficient
    for k in range(1, degree // odd + 1):
        coefficient *= (b + k - 1) / k
        factor[odd * k] = coefficient
    direct = mul(direct, factor, degree)

sigma = [0] + [sum(d for d in range(1, n + 1, 2) if n % d == 0)
               for n in range(1, degree + 1)]
recursive = [F(1)] + [F(0)] * degree
for n in range(1, degree + 1):
    recursive[n] = b * sum(F(sigma[k]) * recursive[n - k]
                           for k in range(1, n + 1)) / n
need(direct == recursive, "odd product and all-order recurrence disagree")

rank_one = [F(1)]
for n in range(1, degree + 1):
    rank_one.append(rank_one[-1] * (b + n - 1) / n)
first_difference = next(n + 1 for n in range(degree + 1)
                        if direct[n] != rank_one[n])
need(first_difference == 4, "first opening difference must be fourth order")
need(direct[:4] == [F(1), F(5, 4), F(45, 32), F(355, 128)],
     "full opening coefficients")
need(rank_one[3] == F(195, 128), "rank-one fourth opening coefficient")
need(direct[3] - rank_one[3] == F(5, 4), "fourth opening correction")

halves = []
for residue in (1, 3):
    half = [F(1)] + [F(0)] * degree
    for exponent in range(residue, degree + 1, 4):
        factor = [F(0)] * (degree + 1)
        coefficient = F(1)
        factor[0] = coefficient
        for k in range(1, degree // exponent + 1):
            coefficient *= (b + k - 1) / k
            factor[exponent * k] = coefficient
        half = mul(half, factor, degree)
    half_recursive = [F(1)] + [F(0)] * degree
    half_sigma = [0] + [sum(d for d in range(1, n + 1)
                            if n % d == 0 and d % 4 == residue)
                        for n in range(1, degree + 1)]
    for n in range(1, degree + 1):
        half_recursive[n] = b * sum(F(half_sigma[k]) * half_recursive[n-k]
                                    for k in range(1, n+1)) / n
    need(half == half_recursive, "half product and recurrence disagree")
    halves.append(half)
need(mul(halves[0], halves[1], degree) == direct,
     "positive/negative determinant factorization")
positive_first_difference = next(n+1 for n in range(degree+1)
                                 if halves[0][n] != rank_one[n])
need(positive_first_difference == 6, "positive-half opening degree")
need(halves[0][5] - rank_one[5] == b, "positive-half sixth-order correction")
for a in (F(1, 10), F(1, 4), F(1, 2), F(3, 4), F(9, 10)):
    need(all((n+a)-a == n for n in range(degree)),
         "positive-half normalized spectrum for every interior twist")

quarter_data = []
for a in (F(1, 4), F(3, 4)):
    pairs = sorted((abs(n + a), n) for n in range(-degree, degree + 1))
    need([p[0] for p in pairs[:degree]] == [F(2 * m + 1, 4)
                                           for m in range(degree)],
         "quarter absolute spectrum")
    gap_label = pairs[0][1]
    need(gap_label == (0 if a == F(1, 4) else -1), "gap branch")
    for n in range(-degree, degree + 1):
        need((n + 1 + a) - (n + a) == 1, "oriented Ward sign")
        need((-(n + 1 + a)) - (-(n + a)) == -1, "printed i-partial sign")
    quarter_data.append({"twist": str(a), "gap_label": gap_label,
                         "positive_toeplitz_defect_label": 0})

# Exact geometric-series control of the trace norm at one rational r.
r, Q = F(1, 3), F(1, 2)
trace_tail = r / (1 - r)
for cutoff in (1, 2, 7):
    partial = sum(r**m for m in range(1, cutoff + 1))
    tail = r**(cutoff + 1) / (1 - r)
    need(partial + tail == trace_tail, "infinite geometric tail identity")
log_bound = Q * r / ((1 - r) * (1 - Q * r))
need(trace_tail == F(1, 2) and log_bound == F(3, 10), "IR bounds")

result = {
    "verdict": "PASS",
    "scope": "exact spectral signs, gap lines, formal coefficients and geometric tails",
    "coefficient_degree_compared": degree,
    "full_opening_coefficients_q1_to_q8": [str(x) for x in direct[:8]],
    "rank_one_opening_coefficients_q1_to_q8": [str(x) for x in rank_one[:8]],
    "first_different_opening_degree": first_difference,
    "first_difference": str(direct[3] - rank_one[3]),
    "positive_half_opening_coefficients_q1_to_q8": [str(x) for x in halves[0][:8]],
    "positive_half_first_different_opening_degree": positive_first_difference,
    "positive_half_first_difference": str(halves[0][5]-rank_one[5]),
    "positive_negative_factorization": "PASS",
    "quarters": quarter_data,
    "trace_norm_tail_at_r_one_third": str(trace_tail),
    "log_bound_at_Q_one_half_r_one_third": str(log_bound),
    "original_admissible_sector_and_EM_parameter_identification": "OPEN",
    "physical_TFPT_completion_claimed": False,
}
payload = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
Path(__file__).with_name("older-circle-heat-results.json").write_text(payload)
print(payload, end="")
