"""Identify all independently chosen two-step partners in the finite full96 menu."""

from collections import Counter, defaultdict
from itertools import combinations, product

import numpy as np

from reconstruction import require


def series_census(profiles, mirrors, coefficients, independent_pairs):
    """Match first derivatives and second skew derivatives at common width zero."""
    rank = len(profiles[0]) + 1
    require(independent_pairs == rank*(rank-1)//2, "a full affine commutator certificate is required")
    vectors = tuple((1, *(coefficients[g] for g in profile)) for profile in profiles)
    histories = tuple(product(range(len(profiles)), repeat=2))
    field_pairs = tuple(combinations(range(rank), 2))
    by_sum, by_jet = defaultdict(list), defaultdict(list)
    for i, j in histories:
        a, b = vectors[i], vectors[j]
        first = tuple(x+y for x,y in zip(a,b))
        skew = tuple(b[u]*a[v]-b[v]*a[u] for u,v in field_pairs)
        by_sum[first].append((i,j))
        by_jet[first,skew].append((i,j))
    require(all(len(bucket)==1 for bucket in by_jet.values()), "two-step jets do not identify the ordered history")
    survivor_map = [{"left":[list(profiles[i]),list(profiles[j])],
                     "right":[list(profiles[mirrors[i]]),list(profiles[mirrors[j]])]}
                    for i,j in histories]
    census = {
        "candidate_pairs":len(histories)**2,
        "after_first_derivative":sum(len(bucket)**2 for bucket in by_sum.values()),
        "first_derivative_fibre_sizes":dict(sorted(Counter(map(len,by_sum.values())).items())),
        "after_second_skew_derivative":len(survivor_map),
        "matched_histories":survivor_map,
        "scope":"Equality as a function of a common variable step width near zero; numerical fixed-width recheck reported separately.",
    }
    return histories, census


def numeric_fixed_width(transfers, histories, entries):
    """Independent finite-width replay using the commutator witness entries."""
    dimension = transfers[0].shape[0]
    positions = list(dict.fromkeys([divmod(k,dimension) for k in entries]+[(0,0),(1,1),(2,2),(3,3)]))
    rows, cols = zip(*positions)
    left = np.stack([matrix[list(rows),:] for matrix in transfers])
    right = np.stack([matrix[:,list(cols)].T for matrix in transfers])
    signatures = np.stack([np.einsum("ij,ij->i",left[j],right[i]) for i,j in histories])
    closest = float("inf")
    below_tolerance = 0
    tolerance = 1e-9
    for i, signature in enumerate(signatures[:-1]):
        separations = np.max(np.abs(signatures[i+1:]-signature),axis=1)
        closest = min(closest,float(np.min(separations)))
        below_tolerance += int(np.sum(separations <= tolerance))
    return {"smallest_distinct_history_gap":closest,"tested_matrix_entries":[list(ij) for ij in positions],
            "unresolved_off_diagonal_pairs":below_tolerance,"tolerance":tolerance,
            "is_interval_enclosed":False}
