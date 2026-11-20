#!/usr/bin/env python3
"""Step 4 (KERNSIX20261120): where do the four marks sit?

Four marks on the seam circle are the branch points of a rectangular torus tau = i t. The action
of the original is a spectral sum of an operator divided by its lowest eigenvalue (the seed scale).
Model: S_f(t) = sum f(lambda / lambda_1) over the Dirac spectrum of the flat torus with sides 1, t.

A  antiperiodic spin structure: S(t) = S(1/t); the square t = 1 is the only minimum
B  the two other types of spin structure: no smooth minimum at the square
C  off the mirror family (twist) the square is not a minimum: the mirror premise is needed
D  the quarter turn is a symmetry of the spectrum only at the square
E  marks on the circle <-> t: the square is the marks at the fourth roots of unity
"""
import numpy as np
from scipy.special import ellipk
from common import Checks

c = Checks("4_marks")
R = 90
mm, nn = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")

def levels(t, a, b, x=0.0):
    """lambda^2 / (2 pi)^2 for the torus with periods 1 and tau = x + i t, spin structure (a, b)."""
    return ((mm + a) * x - (nn + b)) ** 2 / t ** 2 + (mm + a) ** 2

def action(f, t, a, b, x=0.0):
    u = levels(t, a, b, x)
    u = u[u > 1e-14]
    return float(2 * np.sum(f(u / u.min())))       # two spinor components; f as a function of (lambda/lambda_1)^2

profiles = {
    "exp(-x^2/2)": lambda u: np.exp(-u / 2.0),
    "exp(-x^2/8)": lambda u: np.exp(-u / 8.0),
    "exp(-x^2/30)": lambda u: np.exp(-u / 30.0),
    "exp(-x)": lambda u: np.exp(-np.sqrt(u)),
}
ts = np.exp(np.linspace(-1.2, 1.2, 97))            # symmetric grid in log t, t = 1 in the middle

# ---------------------------------------------------------------- A
dual, mins, curv = True, True, True
for name, f in profiles.items():
    S = np.array([action(f, t, 0.5, 0.5) for t in ts])
    dual &= np.allclose(S, S[::-1], rtol=1e-10)
    mins &= int(np.argmin(S)) == len(ts) // 2 and np.all(np.diff(S[len(ts) // 2:]) > 0)
    k = len(ts) // 2
    curv &= S[k + 1] + S[k - 1] - 2 * S[k] > 0
c.ok("A1 antiperiodic on both cycles: S(t) = S(1/t) for every profile, so the square is stationary", dual)
c.ok("A2 the square is the only minimum on the rectangular family (S rises monotonically away from t = 1), four profiles", mins)
c.ok("A3 the minimum is non-degenerate (positive curvature in log t)", curv)
weyl = np.array([t * levels(t, 0.5, 0.5).min() for t in ts])
c.ok("A4 leading term: area times lambda_1^2 = pi^2 (t + 1/t), smallest at the square", np.allclose(weyl, (ts + 1 / ts) / 4))

# ---------------------------------------------------------------- B
f = profiles["exp(-x^2/8)"]
Sp = np.array([action(f, t, 0.0, 0.0) for t in ts])
k = len(ts) // 2
c.ok("B1 periodic on both cycles: the square is a maximum with a kink (both one-sided slopes point down)",
     Sp[k] > Sp[k + 1] and Sp[k] > Sp[k - 1] and abs((Sp[k + 1] - Sp[k]) - (Sp[k] - Sp[k - 1])) > 0.5 * abs(Sp[k + 1] - Sp[k]))
Sm = np.array([action(f, t, 0.5, 0.0) for t in ts])
c.ok("B2 antiperiodic on one cycle only: no symmetry t -> 1/t and no stationary point at the square",
     not np.allclose(Sm, Sm[::-1], rtol=1e-3) and abs(Sm[k + 1] - Sm[k - 1]) > 1e-3 * Sm[k])

# ---------------------------------------------------------------- C
xs = np.array([-0.06, -0.03, 0.0, 0.03, 0.06])
Sx = np.array([action(f, 1.0, 0.5, 0.5, x) for x in xs])
c.ok("C1 twisting the square (tau = x + i) lowers the action on both sides, with a kink: the square is a saddle among all tori",
     Sx[2] > Sx[1] > Sx[0] and Sx[2] > Sx[3] > Sx[4], [round(float(v), 4) for v in Sx - Sx[2]])

# ---------------------------------------------------------------- D
def turn_invariant(t, r=6):
    """Is the set of momentum vectors ((m + 1/2), (n + 1/2)/t) mapped to itself by a quarter turn?"""
    pts = {(round(m + 0.5, 9), round((n + 0.5) / t, 9)) for m in range(-40, 40) for n in range(-40, 40)}
    small = [p for p in pts if p[0] ** 2 + p[1] ** 2 < r ** 2]
    return all((round(-p[1], 9), round(p[0], 9)) in pts for p in small)
c.ok("D1 the quarter turn maps the momenta of the modes to themselves at the square and not at t = 1.1, 2, 3",
     turn_invariant(1.0) and not any(turn_invariant(t) for t in (1.1, 2.0, 3.0)))

# ---------------------------------------------------------------- E
def modulus(theta):
    z = np.exp(1j * np.array([theta, np.pi - theta, np.pi + theta, -theta]))   # four marks, symmetric under both mirrors
    lam = ((z[0] - z[1]) * (z[2] - z[3]) / ((z[0] - z[2]) * (z[1] - z[3]))).real
    return lam, ellipk(1 - lam) / ellipk(lam)
thetas = np.linspace(0.15, np.pi / 2 - 0.15, 41)
tt = np.array([modulus(th)[1] for th in thetas])
lam4, t4 = modulus(np.pi / 4)
c.ok("E1 four marks on the circle have a real cross ratio and give tau = i t; the map angle -> t is monotone",
     np.all(np.diff(tt) != 0) and (np.all(np.diff(tt) > 0) or np.all(np.diff(tt) < 0)))
c.ok("E2 t = 1 exactly for the marks at the fourth roots of unity (cross ratio 1/2)", abs(lam4 - 0.5) < 1e-12 and abs(t4 - 1) < 1e-12)

c.done({"S_minus_S1": {n: [round(action(f_, t, 0.5, 0.5) - action(f_, 1.0, 0.5, 0.5), 6) for t in (0.5, 0.8, 1.25, 2.0)]
                       for n, f_ in profiles.items()}})
