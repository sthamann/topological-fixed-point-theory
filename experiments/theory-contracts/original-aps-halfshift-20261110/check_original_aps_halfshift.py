"""Fixed-label half shifts on the original folded R074 two-face collar."""
from fractions import Fraction as F
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

repo = Path(__file__).resolve().parents[3]
fold = json.loads((repo / "experiments/theory-contracts/v16-operator-pullback-20261014/"
                   "fold-results.json").read_text())
M = sp.Matrix(fold["explicit_Weyl_matrix"])
I = sp.eye(8)
q = M.T * (I[:, 7] - I[:, 6])
r = M.T * (I[:, 4] + I[:, 5])
assert M.T * M == I and q.dot(q) == r.dot(r) == 2 and q.dot(r) == 0
b = r - q
swap = [4, 5, 6, 7, 0, 1, 2, 3]
half_residue = lambda x: x[swap, :]
assert half_residue(q) == q
bad_b = half_residue(r) - q
correction = b - bad_b
Y = sp.Matrix([-2, -2, -2, 3, 3, 0, 0, 0]) / 6
T3 = (I[:, 3] - I[:, 4]) / 2
assert M * correction == 2 * I[:, 4] and correction.dot(correction) == 4
assert [Y.dot(M * a) for a in (b, bad_b, correction)] == [sp.Rational(1, 2), -sp.Rational(1, 2), 1]
assert T3.dot(M * correction) == -1
assert r.dot(-half_residue(r)) == 0 and r - half_residue(r) == correction


def shift(a, n):
    return n + int(8 * a[n % 8])


def charge_from_sea(a):
    charges = [0] * 8
    for k in range(8):
        final = shift(a, k)
        target = final % 8
        charges[target] = (final - target) // 8
    return charges


for n in range(-32, 32):
    assert shift(-q, shift(r, n)) == shift(b, n)
    assert shift(r, shift(-q, n)) == shift(bad_b, n)
    assert shift(correction, shift(r, shift(-q, n))) == shift(b, n)
    assert shift(half_residue(r), shift(-q, n)) == shift(b, n)
    assert shift(-q, shift(q, n)) == n
    assert shift(-half_residue(r), shift(r, n)) == n

for a in (q, r, b, bad_b, correction):
    assert len({shift(a, k) % 8 for k in range(8)}) == 8
    charges = charge_from_sea(a)
    assert all(isinstance(x, int) for x in charges)
    assert sum(charges) == sum(a)
    if all(v.q == 1 for v in a):
        assert charges == list(a)

# The required charge correction has two creation and two annihilation
# legs in the actual folded Fourier chart. Exhibit a separated-arc term.
zeta = sp.sqrt(2) * (1 + sp.I) / 2
def minor(channels, arcs):
    entries = sp.Matrix([[zeta**((j*k) % 8) for j in arcs] for k in channels])
    return sp.simplify(sp.expand_complex(entries.det() / 8))

creations = [j for j in range(8) if correction[j] == 1]
annihilations = [j for j in range(8) if correction[j] == -1]
coefficient = minor(creations, [0, 1]) * sp.conjugate(minor(annihilations, [2, 3]))
assert sp.simplify(coefficient * sp.conjugate(coefficient)) == F(1, 1024)


def energy_from_sea(a, alpha):
    result = F(0)
    for k, excess in enumerate(charge_from_sea(a)):
        occupied = range(0, excess) if excess >= 0 else range(excess, 0)
        result += (1 if excess >= 0 else -1) * sum(8*n + k + alpha for n in occupied)
    return result


mp.mp.dps = 70
beta = 4 * mp.pi


def graph(p):
    t = mp.exp(-beta * p / 2)
    return mp.matrix([[1, t], [t, t*t]]) / (1+t*t)


def frame(p):
    t = mp.exp(-beta * p / 2)
    return mp.matrix([[1, -t], [t, 1]]) / mp.sqrt(1+t*t)


def square_norm(matrix):
    return mp.fsum(abs(entry)**2 for entry in matrix)


rows = {}
for quarter in (F(1, 4), F(3, 4)):
    alpha = mp.mpf(quarter.numerator) / quarter.denominator
    rows[str(quarter)] = {}
    for name, a in (("correct_order", b), ("reversed_fixed_labels", bad_b), ("charge_correction", correction)):
        energy = energy_from_sea(a, quarter)
        assert energy == sum(4*x*x + (k+quarter-4)*x for k, x in enumerate(a))
        defect = mp.mpf(0)
        covariance_trace = [mp.mpf(0) for _ in range(8)]
        worst_graph_residual = mp.mpf(0)
        for k in range(8):
            for n in range(-7, 8):
                p = 8*n + k + alpha
                p_new = p + int(8*a[k])
                g, g_new = graph(p), graph(p_new)
                defect += square_norm(g_new - g)
                covariance_trace[k] += (g - g_new)[1, 1]
                rotation = frame(p_new) * frame(p).T
                worst_graph_residual = max(worst_graph_residual, square_norm(g_new * rotation - rotation * g))
        expected = 2 * sum(abs(x) for x in a) * mp.tanh(2 * beta)
        assert abs(defect - expected) < mp.mpf("1e-55")
        assert max(abs(x-y) for x, y in zip(covariance_trace, a, strict=True)) < mp.mpf("1e-55")
        assert worst_graph_residual < mp.mpf("1e-125")
        rows[str(quarter)][name] = {
            "filled_sea_energy": str(energy),
            "bare_two_face_hilbert_schmidt_squared": mp.nstr(defect, 45),
            "half_frame_graph_transport_max_squared_residual": mp.nstr(worst_graph_residual, 5),
        }

print(json.dumps({
    "fold_weights": {name: [str(x) for x in v] for name, v in
                     (("q", q), ("r", r), ("B", b), ("wrong_B", bad_b), ("correction", correction))},
    "half_shift_sea_charges": {"q": charge_from_sea(q), "r": charge_from_sea(r)},
    "separated_arc_coefficient_abs": "1/32",
    "quarters": rows,
    "physical_R113_weight_computed": False
}, indent=2))
