"""Test the candidate involution against the pre-existing v252 weak source."""
import json
import runpy
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

repo = Path(sys.argv[1]).resolve()
source = repo / "experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/check_actual_weak_source_moment.py"
original = runpy.run_path(str(source))
v, exact = original["v"], original["exact"]
swap = v.Kswap()
gamma = exact(v.grading())
tx = v.rep(0, 1j * np.array([[0, 1], [1, 0]]), np.zeros((3, 3)))
q = exact(-1j * (tx + swap @ tx.conj() @ swap))

# q(H)=(Htilde,H) is the identity frame here. Keep H (second column)
# even, Htilde odd, and use the color factor of the carrier lift.
a = np.diag([-1, 1])
linear = v.rep(1j, -1j * a, 1j * np.eye(3))
g = exact(linear @ swap @ linear.conj() @ swap)
identity = sp.eye(96)
zero = sp.zeros(96)
assert g == g.H and g * g == identity and g * gamma == gamma * g
assert g * q * g == -q and q * gamma == gamma * q
assert (q * q) ** 2 == q * q
left_count = sp.trace(q * q)
assert left_count == 48
assert g[v.IDX[("e", "R", 0, 0)], v.IDX[("e", "R", 0, 0)]] == -1
assert g[v.IDX[("nu", "R", 0, 0)], v.IDX[("nu", "R", 0, 0)]] == 1

data = []
for down in (np.diag([4, 5, 6]), np.eye(3)[[2, 1, 0]] @ np.diag([4, 5, 6])):
    d = exact(v.full_D(np.diag([1, 2, 3]), down, np.diag([2, 3, 4]),
                       np.diag([1, 3, 5]), np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])))
    h0 = gamma / 4 - sp.I * gamma * d
    h1 = q * gamma
    assert g * h0 * g == h0
    assert g * h1 * g == -h1
    assert (g * h1 - h1 * g).H * (g * h1 - h1 * g) == 4 * q * q
    arr = lambda x: np.array(x.tolist(), dtype=complex)
    h0, h1, mark = arr(h0), arr(h1), arr(g)
    positive = (np.eye(96) + mark) / 2
    negative = (np.eye(96) - mark) / 2
    row = []
    for t in (-.7, 0., .7):
        hp, hm = h0 + t * h1, h0 - t * h1
        assert abs(np.linalg.norm(mark @ hp - hp @ mark, 2) - 2 * abs(t)) < 1e-12
        half = expm(-hp / 8)
        other_half = expm(-hm / 8)
        assert np.linalg.norm(mark @ half @ mark - other_half) < 1e-12
        mix = float(np.linalg.norm(negative @ half @ positive, "fro") ** 2)
        assert (mix == 0. if t == 0. else mix > 0.)
        row.append({"t": t, "reflected_odd_weight_from_even": mix})
    data.append(row)
print(json.dumps({
    "scope": "Actual existing quarter-mode 96 component formal weak probe, fixed diagnostic Yukawas; before physical P_adm",
    "exact_covariance": "G h(t) G = h(-t)",
    "exact_generator_defect_operator_norm": "2*abs(t)",
    "exact_generator_defect_frobenius_norm_squared": "192*t^2",
    "weak_left_lift_rank": int(left_count),
    "all_positive_times": "If t != 0, positive heat has a nonzero even-to-odd block by selfadjointness and unique Hermitian logarithm",
    "source_ordered_pair_example": data
}, indent=2))
