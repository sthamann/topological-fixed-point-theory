"""Do the Dirac eta values of the Z_4 and Z_5 cuts match the glue data of the E8 charts D5+A3 and A4+A4?

Pattern test only.  eta from the spectral formula (checked in check_eta_phase_and_period.py):
  eta_eps(Z_p) = (1/p) sum_j eps^j / (2 sin^2(pi j/p)),  eps = +1 / -1 (second spin structure, p even).
Glue data: spinor class of D_n has norm n/4 (conformal weight n/8); A_{p-1} class k has norm k(p-k)/p.
"""
from fractions import Fraction as Fr
from mpmath import mpf, pi, sin
from common import Checks

c = Checks("cut_glue_match")
def eta(p, eps=1):
    v = sum(eps ** j / (2 * sin(pi * mpf(j) / p) ** 2) for j in range(1, p)) / p
    return Fr(str(v)).limit_denominator(10 ** 6)

c.ok("closed forms: eta_+ = (p^2-1)/(6p); eta_- = -(p^2+2)/(12p) for even p (p <= 40)",
     all(eta(p) == Fr(p * p - 1, 6 * p) for p in range(2, 41)) and all(eta(p, -1) == -Fr(p * p + 2, 12 * p) for p in range(2, 41, 2)))
e4 = (eta(4), abs(eta(4, -1)))
c.ok("Z_4 cut: the two spin structures give 5/8 and 3/8 = conformal weights of the D5 and A3 = D3 spin fields (norms 5/4, 3/4)",
     e4 == (Fr(5, 8), Fr(3, 8)) and Fr(5, 4) / 2 == e4[0] and Fr(3, 4) / 2 == e4[1], e4)
sums = {p: eta(p) + abs(eta(p, -1)) for p in range(2, 41, 2)}
c.ok("eta_+ + |eta_-| = p/4: the two values add up to weight 1 (a norm 2 glue vector) only for p = 4",
     all(v == Fr(p, 4) for p, v in sums.items()) and [p for p, v in sums.items() if v == 1] == [4])
split = [n for n in range(0, 9) if Fr(n, 8) == e4[0] and Fr(8 - n, 8) == e4[1]]
c.ok("as a split D_n + D_(8-n) of D8 inside E8 the pair (n/8, (8-n)/8) selects n = 5, i.e. D5 + D3", split == [5])
# look-elsewhere: for which even p is eta_+ a spinor weight n/8 at all?
mult8 = [p for p in range(2, 41, 2) if (eta(p) * 8).denominator == 1]
c.ok("eta_+ is a multiple of 1/8 only for p = 2 and p = 4 among even p <= 40; p = 2 gives 1/4 + 1/4, not weight 1", mult8 == [2, 4], mult8)
# Z_5: eta equals the discriminant norm (p-1)/p of A_(p-1) only for p = 5
eq = [p for p in range(2, 200) if Fr(p * p - 1, 6 * p) == Fr(p - 1, p)]
c.ok("eta(Z_p) = (p-1)/p (norm of the A_(p-1) glue class 1) holds only for p = 5: eta = 4/5", eq == [5])
tw5 = {s: Fr(-(4 - s * (5 - s)), 5) for s in range(5)}
c.ok("Z_5 twisted sectors contain 2/5 = weight of the A4 glue class 1 (norm 4/5), but not 3/5 (class 2): partial match only",
     Fr(2, 5) in tw5.values() and Fr(3, 5) not in [abs(v) for v in tw5.values()], {s: str(v) for s, v in tw5.items()})
# Z_4 twisted sectors: 5/8, -1/8, -3/8, -1/8
tw4 = {s: Fr(5, 8) - Fr(s * (4 - s), 4) for s in range(4)}
c.ok("Z_4 twisted sectors (first spin structure): 5/8, -1/8, -3/8, -1/8, all spinor weights n/8",
     [tw4[s] for s in range(4)] == [Fr(5, 8), Fr(-1, 8), Fr(-3, 8), Fr(-1, 8)])
c.ok("inconsistency recorded: Z_4 matches at the level 'eta = weight = norm/2', Z_5 at the level 'eta = norm'; no common rule",
     e4[0] == Fr(5, 4) / 2 and eta(5) == Fr(4, 5) and eta(5) != Fr(4, 5) / 2)
c.done({"eta4": [str(x) for x in e4], "eta5": str(eta(5)), "tw4": {s: str(v) for s, v in tw4.items()}})
