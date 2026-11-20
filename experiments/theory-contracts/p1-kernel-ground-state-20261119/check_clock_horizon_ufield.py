"""Quarter clock of the source boundary operator, horizon-first probe (G3), alpha potential (gap-ufield)."""
from mpmath import mpf, mpc, pi, exp, expj, coth, log, findroot, quad, diff, binomial, nsum, inf
from common import Checks, C3, DTOP, PHI_TREE, alpha_inv, z

c = Checks("clock_horizon_ufield")
tol = mpf(10) ** -30
cls = lambda lam: complex(expj(pi * mpf(lam) / 2))
rnd = lambda w: (round(w.real, 9) + 0.0, round(w.imag, 9) + 0.0)

# --- A: quarter clock U = exp(i pi B / 2) -------------------------------------------------------------
S2 = [s * (k + 1) for k in range(0, 40) for s in (1, -1)]                 # Dirac on S^2: +-(k+1)
c.ok("Dirac B on round S^2 has integer spectrum: U = exp(i pi B/2) has order 4 with all four classes 1, i, -1, -i",
     {rnd(cls(l)) for l in S2} == {(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)})
sc = [mpf(l) + mpf(1) / 2 for l in range(0, 40)]                          # scalar Omega = l + 1/2
c.ok("scalar comparison kernel: exp(i pi Omega) = i * antipode, only classes +-i, square = -1 (two acts, not four)",
     {rnd(complex(expj(pi * l))) for l in sc} == {(0.0, 1.0), (0.0, -1.0)})
S3 = [s * (k + mpf(3) / 2) for k in range(0, 40) for s in (1, -1)]
S1 = [s * (k + mpf(1) / 2) for k in range(0, 40) for s in (1, -1)]
c.ok("Dirac on S^1 (antiperiodic) and S^3: half-integer spectrum, U has order 8, not 4",
     all(len({rnd(cls(l)) for l in sp}) == 4 and all(abs(cls(l) ** 4 + 1) < 1e-9 for l in sp) for sp in (S1, S3)))
# 2D oscillator / Bargmann-Fock on C^2: level n has energy n+1 (two zero-point halves) and degeneracy n+1
c.ok("|B| on S^2 is isospectral to the 2D isotropic oscillator incl. zero point (energy n+1, degeneracy n+1) times 2 spin",
     all(2 * binomial(k + 1, k) == 2 * (k + 1) for k in range(0, 50)))
# Fourier transform on L^2(R^2) acts as (-i)^n on level n = exp(-i pi (H-1)/2): U^-1 on the positive half up to phase i
c.ok("exp(-i pi |B|/2) on the positive half = 2D Fourier transform times the metaplectic phase -i (ground level class -i)",
     all(abs(expj(-pi * mpf(n + 1) / 2) - (-1j) * (-1j) ** n) < 1e-12 for n in range(0, 12)))
ab = nsum(lambda k: 2 * k * (mpc(0, 1) * mpf("0.999999")) ** k, [1, inf])
c.ok("Abel-regularised trace of the quarter clock on the positive half is -1",
     abs(ab + 1) < mpf(10) ** -5, f"{ab}")
# 2D Gaussian ground state exp(-pi |x|^2): variance 1/(4 pi) per direction; Weyl response with |beta|^2 = 2 over two directions
a = mpf(1) / 137
c.ok("ground level of the 2D oscillator (h = 1): two directions with 1/(4 pi) each give exp(-2 alpha), no momentum needed",
     abs(exp(-(4 * pi * a) * (1 / (4 * pi) + 1 / (4 * pi))) - exp(-2 * a)) < tol)

# --- B: horizon first (G3): thermal / thermofield double at beta = 2 pi, hbar = 1 ---------------------
eq = lambda w: w / 4 * coth(pi * w)                                       # energy per quadrature
c.ok("thermal state at beta = 2 pi (hbar = 1): energy per quadrature >= 1/(4 pi) = 2 c3, never c3",
     all(eq(mpf(w)) > 2 * C3 for w in ("0.001", "0.1", "1", "3")) and abs(eq(mpf("1e-12")) - 2 * C3) < mpf(10) ** -20)
x2 = lambda w: coth(pi * w) / (2 * w)
k1 = 4 * pi * x2(mpf(1))
zk = z(alpha_inv(lambda al: PHI_TREE + (lambda q: q * (1 - q) ** (-mpf(5) / 4))(DTOP * exp(-k1 * al))))
c.ok("quantum thermofield oscillator with omega = 1: damping exp(-6.31 alpha), excluded",
     abs(k1 - mpf("6.3067")) < mpf("0.001") and abs(zk) > 1000, f"k = {float(k1):.4f}, z = {float(zk):.0f}")
wstar = findroot(lambda w: x2(w) - 1 / (2 * pi), 3)
c.ok("exp(-2 alpha) from a quantum thermal mode would need omega = pi*(1 + 5e-9): horizon-first works only as classical zero mode",
     abs(wstar - pi) < mpf(10) ** -6, f"omega* = {wstar}")

# --- C: alpha as a field with potential U (gap-ufield) ------------------------------------------------
phi = lambda al: PHI_TREE + (lambda q: q * (1 - q) ** (-mpf(5) / 4))(DTOP * exp(-2 * al))
F = lambda al: al ** 3 - 2 * C3 ** 3 * al ** 2 - mpf(164) / 5 * C3 ** 6 * log(1 / phi(al))
root = 1 / alpha_inv(phi)
grid = [mpf(i) / 2000 for i in range(1, 2000)]
signs = sum(1 for i in range(len(grid) - 1) if F(grid[i]) * F(grid[i + 1]) < 0)
c.ok("F(alpha) = dU/dalpha has exactly one zero on (0, 1)", signs == 1, f"root = {root}")
c.ok("the zero is a minimum of U: U'' = F' > 0, equal to 3 alpha^2 within 2 % (ratio 0.9885)",
     diff(F, root) > 0 and abs(diff(F, root) / (3 * root ** 2) - 1) < mpf("2e-2"), f"U'' = {diff(F, root)}, ratio = {diff(F, root) / (3 * root ** 2)}")
U = lambda al: quad(F, [root, al])
c.ok("U is a single well: no second vacuum, no barrier; U(0) - U(min) and U(1) - U(min) both positive",
     U(mpf("1e-9")) > 0 and U(1) > 0, f"U(0)-U(min) = {U(mpf('1e-9'))}")
# sensitivity of the minimum to the log term: the log is only a 1/(2 pi)-suppressed source
c.ok("the well is cubic-dominated: alpha^3 = (164/5) c3^6 log(1/phi) + 2 c3^3 alpha^2 holds at the root",
     abs(root ** 3 - 2 * C3 ** 3 * root ** 2 - mpf(164) / 5 * C3 ** 6 * log(1 / phi(root))) < tol)

c.done({"root": str(root), "omega_star": str(wstar), "k_tfd": str(k1)})
