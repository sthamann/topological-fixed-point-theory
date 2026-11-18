"""Test the native Pin lift and reflected charged Fold readouts."""

import json
import sys
from collections import Counter
from pathlib import Path

import sympy as s

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo))
from tfpt_explorer.current_block_geometry import _actual_current_source

base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
data = json.loads((base / "one-original-fold-link-results.json").read_text())
orbit = json.loads((base / "direct-fold-word-results.json").read_text())["source_orbit"]
orbit = [tuple(map(s.Rational, p)) for p in orbit]
zeta = (1 + s.I) / s.sqrt(2)
sqrt2 = s.sqrt(2)

sx = s.Matrix([[0, 1], [1, 0]])
sy = s.Matrix([[0, -s.I], [s.I, 0]])
sz = s.diag(1, -1)
choices = {"I": s.eye(2), "sigma1": sx, "sigma2": sy, "sigma3": sz}
spin_actions = []
for rx in (0, 1):
    for rw in (0, 1):
        for label, u in choices.items():
            # Coefficients of the actual D=-i(sigma1 dx+sigma2 dw).
            if u * sx * u.H != (-1)**rx * sx or u * sy * u.H != -(-1)**rw * sy:
                continue
            sign_x = (-1)**rx
            c_before_u = -s.I if not rx else s.I  # original c=i, after K and rx
            c_after = c_before_u * sign_x
            assert c_after == -s.I
            spin_actions.append({"reflect_x": bool(rx), "reflect_w": bool(rw),
                                 "spin_factor": label, "c_image": "-i",
                                 "preserves_boundary_sign_b": sign_x == 1})
assert len(spin_actions) == 4


def partner(p):
    return tuple(-x for x in reversed(p))


def energy(alpha, p):
    return sum(x**2 / 2 + s.Rational(k + alpha - 4, 8) * x for k, x in enumerate(p))


def edge(p):
    return s.expand(sum(x * zeta**k for k, x in enumerate(p)))


def transform_edge(j):
    return s.simplify(-s.conjugate(zeta) * s.conjugate(j))


def charge_oriented_edge(p):
    j = edge(p)
    return s.simplify(s.re(j if sum(p) >= 0 else transform_edge(j)))


def clock_index(alpha, p, oriented=False):
    q = sum(p)
    e = 8 * energy(alpha, p) - alpha * q
    return e + (q - abs(q)) / 2 if oriented else e


first, middle, last = orbit
p_partner = partner(middle)
assert edge(middle) == sqrt2 and edge(p_partner) == -1 + s.I
assert s.expand((s.re(edge(p_partner)) - sqrt2)**2) == 3 + 2 * sqrt2

cases = []
for alpha in (s.Rational(1, 4), s.Rational(3, 4)):
    companion = 1 - alpha
    for p in orbit:
        assert energy(alpha, p) == energy(companion, partner(p))
        assert s.simplify(edge(partner(p)) - transform_edge(edge(p))) == 0
        assert clock_index(alpha, p, True) == clock_index(companion, partner(p), True)
    raw = [clock_index(companion, partner(p)) for p in orbit]
    fixed = [clock_index(companion, partner(p), True) for p in orbit]
    assert raw == [6, 8, 4] and fixed == [4, 6, 4]
    cases.append({"quarter": str(alpha), "paired_quarter": str(companion),
                  "paired_energies": [str(energy(companion, partner(p))) for p in orbit],
                  "native_partner_phase_word": [int(s.I**int(n)) for n in raw],
                  "transported_phase_word": [int(s.I**int(n)) for n in fixed],
                  "equal_strip_reference_normalized_exponent": str(
                      energy(alpha, middle) + energy(alpha, last) - 2 * energy(alpha, first))})

def parse(key):
    return [tuple(map(s.Rational, row["p"])) for row in data[key]]


even, spinor = parse("even_shell"), parse("spinor_shell")
assert Counter(s.simplify(edge(p)) for p in even) == Counter(transform_edge(edge(p)) for p in spinor)
assert [p for p in even if s.simplify(edge(p) - transform_edge(edge(middle))) == 0] == [
    (0, 0, 1, 0, 1, 0, 0, 0)]
assert all(charge_oriented_edge(p) == charge_oriented_edge(partner(p)) for p in even + spinor)
assert [p for p in even + spinor if charge_oriented_edge(p) == sqrt2] == [middle]
gap = 3 - 2 * sqrt2
outside = [s.simplify((charge_oriented_edge(p) - sqrt2)**2 - gap)
           for p in even + spinor if p != middle]
assert 0 in outside and all(x.is_nonnegative for x in outside)

# Check the charge reflection with the actual Fold and lattice cocycle.
fold = repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = s.Matrix(json.loads(fold.read_text())["explicit_Weyl_matrix"])
reflect = -s.eye(8)[list(reversed(range(8))), :]
ch = _actual_current_source()["chevalley"]


def integral_coords(double_charge):
    assert all(x.q == 1 for x in double_charge)
    c = ch.rs.alpha_coords(tuple(map(int, double_charge)))
    assert all(x.denominator == 1 for x in c)
    return s.Matrix(list(map(int, c)))


U = s.Matrix.hstack(*(integral_coords(M * reflect * M.T * s.Matrix(root))
                      for root in ch.rs.simple))
B = s.Matrix(ch.B)
assert U**2 == s.eye(8)
delta = (U.T * B * U - B).applyfunc(lambda k: k % 2)
assert delta == delta.T


def eta_exponent(v):
    return int(sum(delta[i, i] * v[i] * (v[i] - 1) / 2 for i in range(8)) +
               sum(delta[i, j] * v[i] * v[j] for i in range(8) for j in range(i + 1, 8))) % 2


unit = [s.eye(8)[:, i] for i in range(8)]
lift_square = [((eta_exponent(e) + eta_exponent(U * e)) % 2) for e in unit]
assert lift_square == [0] * 8

print(json.dumps({
    "actual_flat_pin_sign_reversing_antiunitaries": spin_actions,
    "all_seven_even_edge_values_equal_reflected_spinor_multiset": True,
    "spinor_partner_native_edge": str(edge(p_partner)),
    "native_partner_real_second_moment_about_sqrt2": str(3 + 2 * sqrt2),
    "charge_oriented_partner_real_mark": str(charge_oriented_edge(p_partner)),
    "charge_oriented_other_state_squared_gap": str(gap),
    "actual_cocycle_ratio_is_symmetric": True,
    "actual_cocycle_lift_square_exponents": lift_square,
    "cases": cases,
}, indent=2))
