"""Source-scale e/u checks for the original flavor coefficient and transport formulas."""
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(sys.argv[1]).resolve() / "verification"))
from v18_quark_yukawa import TAB, lY, phi0, yuk

mp.mp.dps = 30
rat = sp.Rational
c_up = (rat(1, 2), rat(34, 41), rat(4, 13))
c_down = (rat(117, 110), rat(47, 41), rat(8, 3))
c_lep = (rat(16, 7), rat(4, 3), rat(7, 6))
s = sp.symbols("s", positive=True)
pair_scale = (s, 1 / s, 1)
up_prime = tuple(c * scale for c, scale in zip(c_up, pair_scale))
down_prime = tuple(c * scale for c, scale in zip(c_down, pair_scale))
assert tuple(u / d for u, d in zip(up_prime, down_prime)) == (
    rat(55, 117), rat(34, 47), rat(3, 26)
)
assert sp.prod(up_prime) == sp.prod(c_up)
assert sp.prod(down_prime) == sp.prod(c_down)
assert sp.prod(c_lep) == rat(32, 9)
assert up_prime[-1] == c_up[-1]
assert (c_lep[0] / up_prime[0]) / (c_lep[0] / c_up[0]) == 1 / s
assert (c_lep[0] / (s * c_up[0])) / (c_lep[0] / c_up[0]) == 1 / s
print("Exact: reciprocal pair rescaling preserves all three ratios, both products, leptons and top.")
print("Exact: uniform quark rescaling preserves every quark c-ratio and the determinant powers.")
print("Both change e/u by the reciprocal rescaling factor.")

lam_e, lam_u = (mp.mpf(TAB[f][1]) for f in ("e", "u"))
c_u_tab = yuk("u") / (mp.pi * phi0**4)
lam_e_exact = mp.pi * mp.mpf(16) / 7 * phi0**5 / lY**8
table_ratio = yuk("e") / yuk("u")
rational_ratio = 32 * phi0 / 7
mixed_ratio = 16 * phi0 / (7 * c_u_tab)
half_unit = mp.mpf("0.0005")
lo = lY * (lam_e - half_unit) / (lam_u + half_unit)
hi = lY * (lam_e + half_unit) / (lam_u - half_unit)
assert lo < table_ratio < hi < min(mixed_ratio, rational_ratio) < mp.mpf(1) / 3
assert abs(lam_e_exact - lam_e) > half_unit
print("phi0, lambdaY:", mp.nstr(phi0, 14), mp.nstr(lY, 14))
print("e/u table, exact electron + table up, rational normalization gauge:",
      *(mp.nstr(value, 14) for value in (table_ratio, mixed_ratio, rational_ratio)))
print("e/u conditional table rounding interval:", mp.nstr(lo, 14), mp.nstr(hi, 14))
print("Electron Lambda from exact 16/7:", mp.nstr(lam_e_exact, 14))
print("Up c from v18 table:", mp.nstr(c_u_tab, 14))

# A common factor is compatible with the three printed lepton entries.
normalization_intervals = []
for field, coefficient, exponent in zip(("e", "mu", "tau"), c_lep, (5, 3, 2)):
    length, amplitude = TAB[field]
    rational = mp.mpf(int(coefficient.p)) / int(coefficient.q)
    expected = mp.pi * rational * phi0**exponent / lY**length
    midpoint = mp.mpf(amplitude)
    normalization_intervals.append(((midpoint - half_unit) / expected,
                                    (midpoint + half_unit) / expected))
norm_lo = max(interval[0] for interval in normalization_intervals)
norm_hi = min(interval[1] for interval in normalization_intervals)
assert norm_lo < mp.mpf("0.98") < norm_hi < 1
# The same factor cannot account for the quark table relative to the rational gauge.
q_u = mp.pi * phi0**4 / (2 * lY**7)
assert (lam_u - half_unit) / q_u > norm_hi
print("All three lepton table entries are consistent with one factor in:",
      mp.nstr(norm_lo, 12), mp.nstr(norm_hi, 12))
print("The up table and the rational quark gauge are incompatible with that same factor.")

# A single common matrix with species scalar prefactors requires constant e_i/u_i.
family_ratios = [yuk(e) / yuk(u) for e, u in zip(("e", "mu", "tau"), ("u", "c", "t"))]
def ratio_bounds(lepton, quark):
    l_lep, a_lep = TAB[lepton]
    l_quark, a_quark = TAB[quark]
    a_lep, a_quark = mp.mpf(a_lep), mp.mpf(a_quark)
    power = lY**(l_lep - l_quark)
    return (power * (a_lep - half_unit) / (a_quark + half_unit),
            power * (a_lep + half_unit) / (a_quark - half_unit))

bounds = [ratio_bounds(e, u) for e, u in zip(("e", "mu", "tau"), ("u", "c", "t"))]
assert bounds[0][0] > bounds[1][1] and bounds[1][0] > bounds[2][1]
print("Same-source family ratios e/u, mu/c, tau/t:", *(mp.nstr(x, 12) for x in family_ratios))

# Full H K with a unitary H cancels H from singular values. A common P and
# ||V_f|| < epsilon/2 force sigma_max ratios above (a-epsilon/2)/(a+epsilon/2).
epsilon, excess = sp.symbols("epsilon excess", positive=True)
a = epsilon + excess
difference = (a - epsilon / 2) / (a + epsilon / 2) - rat(1, 3)
assert sp.cancel(difference - 4 * excess / (3 * (3 * epsilon + 2 * excess))) == 0
eta = epsilon / 4
assert (epsilon / 2 - eta).is_positive
assert (2 * eta * epsilon + eta**2).is_positive
assert max(rational_ratio, mixed_ratio, hi) < mp.mpf(1) / 3
print("Exact: a single full unitary times the common gap-stable positive kernel cannot give the source readouts.")

# Work in Q[omega]/(omega^2 + omega + 1) for the archival generator matrices.
w = sp.symbols("omega")
def reduce_cubic_root(matrix):
    return matrix.applyfunc(lambda x: sp.rem(sp.expand(x), w**2 + w + 1, w))

D = sp.diag(1, w, w**2)
P = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
third = P**2 * D * P**-2
forced_fourth = reduce_cubic_root(D * (P * D * P.inv()) * third)
assert forced_fourth == sp.eye(3)
assert reduce_cubic_root(P * third * P.inv() - D) == sp.zeros(3)
assert forced_fourth.trace() == 3 and sp.rem(D.trace(), w**2 + w + 1, w) == 0
print("Exact: the archival D, PDP^-1, P^2DP^-2 force fourth generator I; the stated rotation forces D.")
