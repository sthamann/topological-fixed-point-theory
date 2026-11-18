"""Use the existing v252 weak lift to evaluate the Fold arc transporter."""

import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v252_full_finite_triple as v

sys.path.insert(0, str(repo / "experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review"))
from check_actual_weak_source_moment import exact

zero = np.zeros((3, 3), complex)
sx = np.array([[0, 1], [1, 0]], complex)
T = v.rep(0, 1j * sx, zero)
J = v.Kswap()
Q = exact(-1j * (T + J @ T.conj() @ J))
gamma = exact(v.grading())
I = sp.eye(96)
assert Q.H == Q and Q**3 == Q and Q * gamma == gamma * Q

P_zero = I - Q**2
P_plus = (Q**2 + Q) / 2
P_minus = (Q**2 - Q) / 2
sectors = {0: P_zero, 1: P_plus, -1: P_minus}
for q, proj in sectors.items():
    assert proj**2 == proj and Q * proj == q * proj
assert sum(sectors.values(), sp.zeros(96)) == I
dimensions = {q: int(sp.trace(proj)) for q, proj in sectors.items()}

sz = np.diag([1, -1])
Tz = v.rep(0, 1j * sz, zero)
Qz = exact(-1j * (Tz + J @ Tz.conj() @ J))
weak_flip = I - Qz**2 + sp.I * Qz  # exp(i*pi*Qz/2) in the actual weak lift
assert weak_flip * weak_flip.H == I
assert weak_flip * Q * weak_flip.H == -Q

# Same original R072/R032-era weak diagnostic branch and family data as
# ACTUAL-WEAK-CONNECTION-AND-JOINT-SOURCE; no new field values selected.
yu = np.diag([1, 2, 3])
yd = np.diag([4, 5, 6])
ye = np.diag([2, 3, 4])
yn = np.diag([1, 3, 5])
mr = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])
D = exact(v.full_D(yu, yd, ye, yn, mr))

x = sp.symbols("x", real=True)
W = P_zero + sp.exp(sp.I * x) * P_plus + sp.exp(-sp.I * x) * P_minus
assert (W * W.H - I).applyfunc(sp.simplify) == sp.zeros(96)
assert W.subs(x, 0) == I
assert (W.subs(x, sp.pi / 4)**8 - I).applyfunc(sp.simplify) == sp.zeros(96)
assert W.subs(x, sp.pi / 4) != I

def squared_norm(a):
    return sp.expand(sp.trace(a.H * a))

blocks = {
    (q, r): squared_norm(Pq * D * Pr)
    for q, Pq in sectors.items()
    for r, Pr in sectors.items()
}
assert blocks[(1, -1)] == blocks[(-1, 1)] == 0
coefficient = squared_norm(Q * D - D * Q)
assert coefficient > 0
assert squared_norm(P_zero * D - D * P_zero) == coefficient
# The exact sector expression is equivalent to 4*sin(x/2)^2 times the
# commutator norm of Q since every nonzero cross-sector block has |q-r|=1.
assert sum((q - r)**2 * n for (q, r), n in blocks.items()) == coefficient
defect = squared_norm(W * D - D * W) - (2 - 2 * sp.cos(x)) * coefficient
assert sp.trigsimp(sp.expand_complex(defect)) == 0

# At half a weak turn the holonomy is central under the actual weak group.
# The root continued from the geometric scalar branch is not.
weak_half_turn = P_zero - Q**2
assert weak_flip * weak_half_turn * weak_flip.H == weak_half_turn
root_limit = W.subs(x, sp.pi / 8)
root_non_covariance = squared_norm(weak_flip * root_limit * weak_flip.H - root_limit)
assert sp.simplify(sp.expand_complex(root_non_covariance)) == 96 - 48 * sp.sqrt(2)

# At a=1 the total gauge holonomy is unchanged, but the older correction
# keeps the arc eigenvalue zeta**(m+q); verify its shifted passband.
zeta = (1 + sp.I) / sp.sqrt(2)
a1_selected = {}
for q in sectors:
    a1_selected[q] = [
        m for m in range(8)
        if sp.simplify(sum(zeta**(j * (m + q)) for j in range(8)) / 8) == 1
    ]
    assert a1_selected[q] == [(-q) % 8]

print(json.dumps({
    "operator_dimension": 96,
    "weak_sector_dimensions": dimensions,
    "arc_length": "pi/4",
    "fixed_quarter_average_passbands_at_a_1": a1_selected,
    "at_a_1_weak_loop_is_identity_but_arc_is_nontrivial": True,
    "squared_weak_dirac_commutator": int(coefficient),
    "squared_weak_singlet_projector_dirac_commutator": int(coefficient),
    "squared_arc_dirac_commutator_at_a_1": str(sp.simplify((2 - sp.sqrt(2)) * coefficient)),
    "squared_root_limit_covariance_defect_at_central_half_turn": "96 - 48*sqrt(2)",
    "cross_sector_dirac_block_squared_norms": {
        f"{q},{r}": int(n) for (q, r), n in blocks.items() if q != r
    },
}, indent=2))
