"""Exact E8 weak-root and native normal-equator CAR checks."""
from collections import Counter
from fractions import Fraction
from math import exp, fsum
from pathlib import Path
import sys

import numpy as np

root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(root / "verification"))
from v1_e8_glue import e8_roots
import v252_full_finite_triple as native
from v330_qgamb_admissible_measure import symmetric_transfer

sys.path.insert(0, str(root / "experiments/theory-contracts/original-admissible-source-closure-20261025"))
from check_joint_lift_native_clock import cocharacter, roots_doubled

roots = set(roots_doubled())
assert roots == {tuple(int(2*x) for x in r) for r in e8_roots()}
assert len(roots) == 240

hypercharge6 = (-2, -2, -2, 3, 3, 0, 0, 0)
weak_cartan2 = (0, 0, 0, 1, -1, 0, 0, 0)
weak_root2 = tuple(2*x for x in weak_cartan2)

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

neutral = {r for r in roots if dot(r, hypercharge6) == 0}
assert len(neutral) == 28
assert {r for r in neutral if dot(r, weak_cartan2) == 4} == {weak_root2}
assert {r for r in neutral if dot(r, weak_cartan2) == -4} == {tuple(-x for x in weak_root2)}
assert sum(dot(r, weak_cartan2) == 0 for r in neutral) == 26
assert len({r for r in neutral if r[3] or r[4]}) == 10
for a in range(-8, 9):
    h = cocharacter((0, 0, 0, a, 1-a), (-1, -1, 0, 0))
    assert dot(h, weak_root2) == 2*(2*a-1)

swap = native.Kswap()
def lift(matrix):
    left = native.rep(0, matrix, np.zeros((3, 3)))
    return left - swap @ left.conj() @ swap

plus = (lift(np.array([[0, 1], [1, 0]])) +
        1j*lift(np.array([[0, -1j], [1j, 0]])))/2
minus = plus.conj().T
cartan2 = lift(np.diag([1, -1]))
grading = np.diag(native.grading()).astype(int)
assert np.array_equal(plus @ minus - minus @ plus, cartan2)
assert np.array_equal(cartan2 @ plus - plus @ cartan2, 2*plus)
zero = np.zeros((3, 3))
majorana = native.full_D(zero, zero, zero, zero, np.eye(3))
assert not np.any(majorana @ plus) and not np.any(plus @ majorana)

# Integrate the native hypercharge of the real bimodule in its own basis.
# d/dt pi(exp(i t/2),I2,exp(-i t/6) I3) at zero equals i times this diagonal.
derivative = native.rep(3j, np.zeros((2, 2)), -1j*np.eye(3))
hypercharge6_native = np.diag((derivative - swap @ derivative @ swap)/1j).real.astype(int)
edges = list(zip(*np.nonzero(plus)))
assert len(edges) == 24 and set(abs(plus[i,j]) for i,j in edges) == {1.0}
assert all(cartan2[i,i] == 1 and cartan2[j,j] == -1 and grading[i] == grading[j]
           and hypercharge6_native[i] == hypercharge6_native[j] for i,j in edges)
sector_sizes = Counter((int(grading[i]), int(hypercharge6_native[i])) for i,j in edges)
assert sector_sizes == {(1, -3): 3, (1, 1): 9, (-1, 3): 3, (-1, -1): 9}

def signed_sea(alpha, co_transport):
    """Crossings need only n=-1 and n=0, checked with an extra halo."""
    crossings = Counter()
    contact = Counter()
    for i,j in edges:
        g = int(grading[i])
        # Only the difference of the two clutch weights enters. The
        # common integer changes the origin, not these crossings.
        # Use actual m_L=-1 and its conjugate: (0,-1) and (1,0).
        lower_degree = -1 if g == 1 else 0
        upper_degree = lower_degree + 1
        for n in range(-4, 5):
            p_in = n + alpha - co_transport*lower_degree
            p_out = n + 1 + alpha - co_transport*upper_degree
            occupied_in, occupied_out = g*p_in < 0, g*p_out < 0
            if occupied_in and not occupied_out:
                crossings[(g, "plus", g*(p_out-p_in))] += 1
                contact[g] -= 1
            if occupied_out and not occupied_in:
                crossings[(g, "minus", g*(p_in-p_out))] += 1
                contact[g] += 1
    return crossings, contact

for alpha in (Fraction(1, 4), Fraction(3, 4)):
    fixed, cocycle = signed_sea(alpha, False)
    carried, carried_cocycle = signed_sea(alpha, True)
    assert fixed == {(1, "plus", 1): 12, (-1, "minus", 1): 12}
    assert cocycle == {1: -12, -1: 12}
    assert not carried and not carried_cocycle

    beta, time = 2.0, 0.7
    def occupation(energy):
        return 1/(1+exp(beta*energy))
    for orientation in (-1, 1):
        direct = 12*fsum(
            occupation(g*(n+float(alpha))) *
            (1-occupation(g*(n+orientation+float(alpha)))) *
            exp(-time*g*orientation)
            for g in (-1, 1) for n in range(-40, 41)
        )
        expected = 12*(exp(-time)+exp(-(beta-time)))/(1-exp(-beta))
        assert abs(direct-expected) < 1e-13

# Round Dirac towers for adjacent degrees use the exact Casimir values.
# Interlacing is proved for all levels in the derivation, these checks
# additionally compare the two charge conventions and finite spectra.
for k in range(-12, 13):
    degrees = (k, k+1)
    bands = [{n*(n+abs(d)) for n in range(1, 129)} for d in degrees]
    assert not (bands[0] & bands[1])
    assert sum(degrees) % 2 == 1

# Exact chiral boundary intersection and cokernel. The finite range
# contains all discrepant Fourier modes for the tested line degrees.
halo = set(range(-16, 17))
north = {n for n in halo if n >= 0}
for d in range(-12, 13):
    south = {n for n in halo if n <= d-1}
    next_south = {n for n in halo if n <= d}
    assert next_south - south == {d}
    assert len(north & south) == max(d, 0)
    assert len(halo - (north | south)) == max(-d, 0)
    assert {n+1 for n in south if n+1 in halo} ^ next_south == {-16}
assert 0 in north and -1 not in north

# The declared scalar transfer has no positive energy of multiplicity two.
transfer, _, lam2, lam3 = symmetric_transfer()
assert 0 < lam3 < lam2 < 1
assert np.allclose(np.linalg.eigvalsh(transfer), [lam3, lam2, 1], rtol=0, atol=1e-14)

print("E8: 240 roots, 28 Y-neutral, precisely one T3=+1 and one T3=-1")
print("v252: weak positive-root ranks by (gamma_F,6Y):", dict(sorted(sector_sizes.items())))
print("both original quarters: fixed sea gives 12 unit-gap returns per order")
print("signed central contacts: -12 and +12; co-transported sea gives no returns")
print("positive-temperature fixed-sea KMS return agrees with telescoping formula")
print("adjacent round Dirac bands: no shared positive eigenvalue in checked range")
print("normal cap Fredholm indices differ by one; v330 excited eigenvalues are simple")
