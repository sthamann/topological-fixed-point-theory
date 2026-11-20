"""Sphere <-> torus bridge and the 'three routes are one' test.

Z+(beta) = Tr_{B>0} exp(-beta B) for the Dirac boundary operator on the round unit S^2
(eigenvalue k, multiplicity 2k, k = 1, 2, ...).  Winding sum over m thermal windings versus
the Eisenstein series E2 on the square torus tau = i (q = exp(-2 pi), i.e. beta = 2 pi).
"""
from mpmath import mpf, pi, exp, sqrt, nsum, inf, log, diff, gamma
from common import Checks, C3, DTOP, PHI_TREE, alpha_inv, z

c = Checks("torus_bridge")
tol = mpf(10) ** -30

def Zplus(beta):                       # Dirac half heat trace on S^2, closed form 2q/(1-q)^2
    q = exp(-beta)
    return 2 * q / (1 - q) ** 2

def E2(q):
    return 1 - 24 * nsum(lambda n: n * q ** n / (1 - q ** n), [1, inf])

# 1 closed form of the half heat trace
direct = nsum(lambda k: 2 * k * exp(-mpf(3) * k), [1, inf])
c.ok("Dirac half heat trace on S^2: sum 2k e^{-beta k} = 2q/(1-q)^2", abs(direct - Zplus(3)) < tol)

# 2 it is the beta-derivative of the circle spectrum (chiral boson modes n = 1, 2, ...)
circ = lambda b: nsum(lambda n: exp(-b * n), [1, inf])
c.ok("S^2 Dirac half trace = -2 d/dbeta of the unit-circle mode trace (sphere channel k <-> circle mode k)",
     abs(Zplus(3) + 2 * diff(circ, 3)) < mpf(10) ** -25)

# 3 Eisenstein series as winding sum of the sphere trace
for b in (mpf(1), 2 * pi, mpf(7)):
    lhs = E2(exp(-b)); rhs = 1 - 12 * nsum(lambda m: Zplus(m * b), [1, inf])
    if abs(lhs - rhs) > tol: break
c.ok("E2(q = e^{-beta}) = 1 - 12 * sum over windings m of Z+_{S^2}(m beta) (beta = 1, 2 pi, 7)", abs(lhs - rhs) < tol)

# 4 square torus: E2(i) = 3/pi, so the winding sum at beta = 2 pi is 1/12 - 2 c3
E2i = E2(exp(-2 * pi)); W = nsum(lambda m: Zplus(2 * pi * m), [1, inf])
c.ok("at the Unruh/KMS value beta = 2 pi (square torus tau = i): E2 = 3/pi and winding sum = 1/12 - 2 c3",
     abs(E2i - 3 / pi) < tol and abs(W - (mpf(1) / 12 - 2 * C3)) < tol, f"W = {W}")

# 5 energy of a chiral boson on the unit circle at beta = 2 pi: thermal 1/24 - c3, with Casimir -1/24 -> -c3
Eth = nsum(lambda n: n / (exp(2 * pi * n) - 1), [1, inf])
c.ok("chiral boson on the square torus: oscillator energy incl. Casimir is exactly -c3 (thermal part 1/24 - c3)",
     abs(Eth - (mpf(1) / 24 - C3)) < tol and abs(Eth - mpf(1) / 24 + C3) < tol, f"E_osc = {Eth - mpf(1)/24}")

# 6 the zero mode adds T/4 = 1/(4 beta) (chiral half of Z ~ beta^{-1/2}); at beta = 2 pi this is +c3: total 0
zero_mode = lambda b: -diff(lambda x: log(x ** (-mpf(1) / 4)), b)
c.ok("zero-mode equipartition energy of the chiral half is T/4; at beta = 2 pi it is +c3 and cancels the oscillators",
     abs(zero_mode(2 * pi) - C3) < mpf(10) ** -25 and abs(zero_mode(2 * pi) + Eth - mpf(1) / 24) < mpf(10) ** -25)

# 7 the cancellation is forced: E2* = E2 - 3/(pi Im tau) has weight 2, so it vanishes at every S-fixed point
def E2star(b):                         # purely imaginary tau = i b/(2 pi)
    return E2(exp(-b)) - 3 / (pi * b / (2 * pi))
ok7 = abs(E2star(2 * pi)) < tol and all(abs(E2star(b)) > mpf(10) ** -3 for b in (mpf(5), mpf(6), mpf(7), mpf(8)))
# weight-2 behaviour under S: E2*(1/b') = -(b')^2 E2*(b') with b' = Im tau
t = mpf("1.3"); s_rel = E2star(2 * pi / t) + t ** 2 * E2star(2 * pi * t)
c.ok("E2* transforms with weight 2 under the quarter turn S, hence vanishes at tau = i only (on the imaginary axis)",
     ok7 and abs(s_rel) < mpf(10) ** -25, f"S-relation residual {s_rel}")

# 8 hexagonal torus tau = rho: also anomaly-only, value (2/sqrt3) c3 -> not c3
def E2_complex(tau):
    from mpmath import mpc, expj
    q = expj(2 * pi * tau)
    return 1 - 24 * nsum(lambda n: n * q ** n / (1 - q ** n), [1, inf])
from mpmath import mpc
rho = mpc(-mpf(1) / 2, sqrt(3) / 2)
E2r = E2_complex(rho)
c.ok("hexagon tau = rho: E2/24 = (2/sqrt3) c3 (anomaly only as well); c3 itself needs the square",
     abs(E2r / 24 - 2 / sqrt(3) * C3) < mpf(10) ** -25, f"E2(rho)/24 = {E2r/24}")

# 9 the three routes are one number: Unruh T (hbar=1), hbar at h=1, anomaly 1/(Im tau) -> all 1/(2 pi), c3 = that / 4
T_unruh, hbar_h1, anomaly = 1 / (2 * pi), 1 / (2 * pi), 3 / pi / 6
c.ok("KMS route 1/(I beta), vacuum route (h = 1), torus route E2(i)/24: all equal (1/(2 pi))/4, not independent",
     abs(T_unruh / 4 - C3) < tol and abs(hbar_h1 / 4 - C3) < tol and abs(anomaly / 4 - C3) < tol)

# 10 I = 4 decomposes as 2 (one Gaussian zero-mode integral) x 2 (chiral half)
full = lambda b: -diff(lambda x: log(x ** (-mpf(1) / 2)), b)
c.ok("I = 4 = 2 x 2: full boson zero mode gives T/2 (equipartition), its chiral half T/4",
     abs(full(2 * pi) - 2 * C3) < mpf(10) ** -25)

# 11 damping as thermal Debye-Waller factor: classical zero mode <x^2> = T at unit stiffness, e^2 = 4 pi alpha, |b|^2 = 2
a = mpf(1) / 137
dw = exp(-(4 * pi * a) * 2 * T_unruh / 2)
c.ok("exp(-2 alpha) = thermal Debye-Waller factor of one unit-stiffness zero mode at T = 1/(2 pi) (one quadrature)",
     abs(dw - exp(-2 * a)) < tol)

# 12 is (1-q)^(-5/4) the first factor of a full mode product prod_n (1-q^n)^(-5/4)?  alpha cannot tell
single = alpha_inv(lambda al: PHI_TREE + (lambda q: q * (1 - q) ** (-mpf(5) / 4))(DTOP * exp(-2 * al)))
def prod_phi(al):
    q = DTOP * exp(-2 * al); p = mpf(1)
    for n in range(1, 30): p *= (1 - q ** n) ** (-mpf(5) / 4)
    return PHI_TREE + q * p
full_prod = alpha_inv(prod_phi)
dz = (full_prod - single) / mpf("0.000000021")
c.ok("full mode product versus single factor shifts alpha^-1 by far less than 1 sigma: undecidable by alpha",
     abs(dz) < mpf("0.1"), f"shift = {float(dz):.3e} sigma")

# 13 exponent 5/4 as two-point exponent 2h of the spin field of 2n = 10 real fermions (h = n/8, n = 5)
h = {n: mpf(n) / 8 for n in (4, 5, 6)}
zs = {n: z(alpha_inv(lambda al, n=n: PHI_TREE + (lambda q: q * (1 - q) ** (-2 * h[n]))(DTOP * exp(-2 * al)))) for n in h}
c.ok("(1-q)^(-5/4) = two-point function exponent 2h of the D5 spin field (h = 5/8); D4 (2h = 1) and D6 (2h = 3/2) fail",
     2 * h[5] == mpf(5) / 4 and abs(zs[5]) < 3 and abs(zs[4]) > 5 and abs(zs[6]) > 5, {n: float(v) for n, v in zs.items()})

c.done({"E2_i": str(E2i), "winding_sum": str(W), "z_spin_field": {n: str(v) for n, v in zs.items()}, "prod_shift_sigma": str(dz)})
