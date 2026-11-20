"""Is the square torus (tau = i) dynamically selected?  f(tau) = log( sqrt(Im tau) |eta(tau)|^2 ) = -log Z_boson."""
from mpmath import mpf, mpc, pi, exp, expj, log, sqrt, diff, nprod, inf, fabs
from common import Checks, C3

c = Checks("square_stability")

def f(x, y):
    tau = mpc(x, y); q = expj(2 * pi * tau)
    eta = expj(pi * tau / 12) * nprod(lambda n: 1 - q ** n, [1, inf])
    return log(sqrt(y) * fabs(eta) ** 2)

h = mpf("1e-6")
d2 = lambda g, t: (g(t + h) - 2 * g(t) + g(t - h)) / h ** 2
d1 = lambda g, t: (g(t + h) - g(t - h)) / (2 * h)

fy_i, fx_i = d1(lambda y: f(0, y), mpf(1)), d1(lambda x: f(x, 1), mpf(0))
c.ok("tau = i is a critical point of the modular-invariant free energy (both first derivatives vanish)",
     abs(fy_i) < mpf(10) ** -10 and abs(fx_i) < mpf(10) ** -10, f"{fy_i}, {fx_i} (central differences, h = 1e-6)")
fyy, fxx = d2(lambda y: f(0, y), mpf(1)), d2(lambda x: f(x, 1), mpf(0))
c.ok("Laplacian of f at tau = i equals -1/2 (only the zero mode sqrt(Im tau) is non-harmonic)",
     abs(fxx + fyy + mpf(1) / 2) < mpf(10) ** -9, f"f_yy = {fyy}, f_xx = {fxx}")
c.ok("the square is a SADDLE: along temperature (Im tau) Z has a minimum, along twist (Re tau) a maximum",
     fyy < 0 and fxx > 0)
r3 = sqrt(3) / 2
gy, gx = d2(lambda y: f(mpf(1) / 2, y), r3), d2(lambda x: f(x, r3), mpf(1) / 2)
c.ok("the hexagon tau = rho is the true extremum (both curvatures negative, f maximal, Z minimal)",
     gy < 0 and gx < 0 and f(mpf(1) / 2, r3) > f(0, 1), f"f(rho) = {f(mpf(1)/2, r3)}, f(i) = {f(0, 1)}")
# along the reflection-symmetric line Re tau = 0 (rectangular tori) the square is the unique maximum of f
ys = [mpf(k) / 20 for k in range(8, 61)]
vals = [f(0, y) for y in ys]
c.ok("with twist frozen by the reflection (rectangular tori, Re tau = 0) the square is the unique extremum",
     max(vals) == f(0, 1) and all(abs(f(0, y) - f(0, 1 / y)) < mpf(10) ** -25 for y in ys[:6]))
# energy readings at the two fixed points: oscillator energy -E2/24 -> -c3 (square), -(2/sqrt3) c3 (hexagon)
c.ok("if the twist were free the surface would roll to the hexagon and c3 would become (2/sqrt3) c3 = 0.045944",
     abs(2 / sqrt(3) * C3 - mpf("0.0459440746")) < mpf(10) ** -9)
# first-law reading: dE = T dS with T = kappa/(2 pi), S = A/4  ->  dE = (kappa/(8 pi)) dA = c3 kappa dA
c.ok("c3 = T_Unruh * (1/4): the same '/4' as entropy per horizon area in dM = (kappa / 8 pi) dA (arithmetic reading only)",
     abs(1 / (2 * pi) * mpf(1) / 4 - C3) < mpf(10) ** -35)
c.done({"f_yy": str(fyy), "f_xx": str(fxx), "g_yy": str(gy), "g_xx": str(gx)})
