"""Does the original force 'round S^2' and 'positive half'?

Source facts used (read, not computed):
  * boundary kernel source: Sigma = boundary of M_+ is a general compact oriented boundary; no roundness stated;
    the polarization is the Calderon projector of the APS-realized problem; chi_seed uses lambda_1^+ ;
  * paper section 013: flat spherical Z_2 branch class, four marks by Gauss-Bonnet, common divisor mu_4 in P^1.
Computed here: what these data give instead of a round sphere.
"""
from mpmath import mpf, mpc, pi, sqrt, ellipk, expj, log, fabs, nprod, inf
from common import Checks, C3

c = Checks("original_forcing")
tol = mpf(10) ** -25

# 2 double cover of P^1 branched over mu_4 = {1, i, -1, -i}: cross-ratio and j-invariant
z1, z2, z3, z4 = mpc(1), mpc(0, 1), mpc(-1), mpc(0, -1)
lam = ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))
j = 256 * (lam ** 2 - lam + 1) ** 3 / (lam ** 2 * (lam - 1) ** 2)
c.ok("marks mu_4: cross-ratio 2 (orbit {-1, 1/2, 2}), j = 1728: the double cover is exactly the square torus tau = i",
     abs(lam - 2) < tol and abs(j - 1728) < tol, f"lambda = {lam}, j = {j}")
tau = lambda m: 1j * ellipk(1 - m) / ellipk(m)
c.ok("tau(lambda = 1/2) = i from the period ratio i K(1-lambda)/K(lambda)", abs(tau(mpf(1) / 2) - 1j) < tol)
# 3 concyclic marks <-> real cross-ratio in (0,1) <-> rectangular torus; free energy along that line
def f(y):
    q = expj(2 * pi * mpc(0, y))
    return log(sqrt(y) * fabs(expj(pi * mpc(0, y) / 12) * nprod(lambda n: 1 - q ** n, [1, inf])) ** 2)
ms = [mpf(k) / 40 for k in range(2, 39)]
vals = [f(tau(m).imag) for m in ms]
best = ms[vals.index(max(vals))]
c.ok("four marks on a circle (real cross-ratio) give rectangular tori; the free energy has its only extremum at lambda = 1/2, i.e. mu_4",
     best == mpf(1) / 2 and all(abs(f(tau(m).imag) - f(tau(1 - m).imag)) < mpf(10) ** -20 for m in ms[:8]), f"argmax lambda = {best}")
# 4 other position of the marks changes the number: lambda = 1/4 -> tau = i*1.2792..., oscillator energy != -c3
from mpmath import nsum, exp
def E_osc(y):   # chiral boson on rectangular torus, units of the space circle 2 pi
    b = 2 * pi * y
    return nsum(lambda n: n / (exp(b * n) - 1), [1, inf]) - mpf(1) / 24
y4 = tau(mpf(1) / 4).imag
c.ok("moving the marks moves the value: lambda = 1/4 gives oscillator energy -0.0413, not -c3 = -0.0398 (marks matter)",
     abs(E_osc(1) + C3) < tol and abs(E_osc(y4) + C3) > mpf("1e-4"), f"Im tau = {y4}, E = {E_osc(y4)}")
# 5 'round' is not needed on average: Weyl law in two dimensions gives 1/(2 pi) per unit of spectrum on ANY 2D cut
#   flat square torus, Dirac (rank 2), positive half: modes k in (2 pi/L) Z^2, eigenvalue |k|, one positive mode per k
L = mpf(1); K = 400
cnt = mpf(0)
R2 = (K * float(L) / (2 * 3.141592653589793)) ** 2
N = int(R2 ** 0.5) + 1
for m in range(-N, N + 1):
    for n in range(-N, N + 1):
        r2 = m * m + n * n
        if 0 < r2 <= R2: cnt += 1 / (2 * pi / L * sqrt(mpf(r2)))
weight = cnt / L ** 2 / K
c.ok("flat pillow/torus: sum of 1/|lambda| per area per unit spectrum tends to 1/(2 pi) = 4 c3 (Weyl law); round S^2 has it exactly per channel",
     abs(weight - 1 / (2 * pi)) < mpf("2e-3"), f"{weight} vs {1/(2*pi)}")
# 7 flat pillow has no integer spectrum: the quarter clock there is the geometric rotation by i on Z[i], not exp(i pi B/2)
spec = sorted({m * m + n * n for m in range(-3, 4) for n in range(-3, 4)} - {0})[:6]
c.ok("square torus spectrum |m + i n| = 1, sqrt2, 2, sqrt5, ...: not integer spaced; order 4 acts geometrically (multiplication by i)",
     spec == [1, 2, 4, 5, 8, 9])
c.done({"lambda_mu4": str(lam), "j": str(j), "weyl_weight": str(weight)})
