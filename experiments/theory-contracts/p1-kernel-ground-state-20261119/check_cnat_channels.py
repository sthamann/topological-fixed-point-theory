"""C_nat test: does the P1 boundary kernel carry the Fourier-self-dual ground state weight 1/(4 pi)?

Two kernels on the round unit sphere S^{d-1}:
  scalar comparison kernel  C = 1/(2 kappa Omega),  Omega = l + (d-2)/2   (F1-ORIGIN-SCENARIOS section 5)
  source seam response      |B_Sigma|^{-1}, Dirac on S^n, eigenvalues +-(k + n/2)  (boundary kernel source)
Channel weight = (degeneracy / sphere volume) * (kernel eigenvalue): the coincidence value contributed
by one angular channel.
"""
from mpmath import mpf, pi, gamma, binomial, exp, legendre, nsum, inf, zeta, altzeta
from common import Checks, C3, DTOP, PHI_TREE, alpha_inv, z

c = Checks("cnat_channels")

def vol_sphere(n):            # volume of unit S^n
    return 2 * pi ** (mpf(n + 1) / 2) / gamma(mpf(n + 1) / 2)

def deg_scalar(l, d):         # harmonic polynomials of degree l in R^d
    if d == 2:
        return 1 if l == 0 else 2
    return binomial(l + d - 1, d - 1) - (binomial(l + d - 3, d - 1) if l >= 2 else 0)

def w_scalar(l, d, kappa=1):  # channel weight of 1/(2 kappa Omega) on S^{d-1}
    om = l + mpf(d - 2) / 2
    return deg_scalar(l, d) / vol_sphere(d - 1) / (2 * kappa * om)

def deg_dirac(k, n):          # multiplicity of eigenvalue +(k+n/2) of Dirac on S^n
    return 2 ** (n // 2) * binomial(k + n - 1, k)

def w_dirac(k, n):            # channel weight of |B|^{-1} on the positive spectral half (Calderon side) of S^n
    return deg_dirac(k, n) / vol_sphere(n) / (k + mpf(n) / 2)

# 1 scalar kernel on S^2: every channel carries exactly 1/(4 pi kappa)
ws = [w_scalar(l, 3) for l in range(0, 40)]
c.ok("scalar D-to-N kernel on S^2: every angular channel weighs exactly 1/(4 pi) = 2 c3 at kappa=1",
     all(abs(w - 2 * C3) < mpf(10) ** -35 for w in ws), f"w_l = {ws[0]}")

# 2 flatness selects d = 3
flat = {d: max(abs(w_scalar(l, d) - w_scalar(1, d)) for l in range(1, 30)) < mpf(10) ** -30 for d in range(2, 9)}
c.ok("flat scalar channel spectrum holds only for S^2 in R^3 (d = 2..8 tested)",
     flat == {d: d == 3 for d in range(2, 9)}, flat)

# 3 Dirac source kernel on S^2: flat, value 1/(2 pi) summed over both signs = 1/(4 pi) per sign
wd = [w_dirac(k, 2) for k in range(0, 40)]
c.ok("source seam response |B_Sigma|^-1 on round S^2, positive spectral half: every channel weighs 1/(2 pi) = 4 c3 (2 c3 per spinor component; both signs together 1/pi)",
     all(abs(w - 4 * C3) < mpf(10) ** -35 for w in wd), f"w_k = {wd[0]}")

# 4 Dirac flatness selects n = 2
flatd = {n: max(abs(w_dirac(k, n) - w_dirac(0, n)) for k in range(0, 30)) < mpf(10) ** -30 for n in range(1, 8)}
c.ok("flat Dirac channel spectrum holds only on S^2 (n = 1..7 tested)",
     flatd == {n: n == 2 for n in range(1, 8)}, flatd)

# 5 chi_seed = 1/lambda_1^+ on round unit S^2 is 1; times lowest-channel density 2*2/(4 pi) -> 1/(2 pi)
lam1 = mpf(1)
c.ok("chi_seed = 1/lambda_1^+(|B|) = 1 on the unit S^2; times lowest positive channel density 2/(4 pi) gives 1/(2 pi) = 4 c3",
     abs(deg_dirac(0, 2) / vol_sphere(2) / lam1 - 4 * C3) < mpf(10) ** -35)

# 6 ground state: psi0 = 2^(1/4) exp(-pi x^2), variance 1/(4 pi) with h = 1.
#   channel weight = hbar_eff / 2  ->  hbar_eff = 1/(2 pi): the sphere area supplies h = 1
hbar_eff = 2 * ws[0]
c.ok("each S^2 channel is an h = 1 oscillator quadrature: weight = hbar_eff/2 with hbar_eff = 1/(2 pi)",
     abs(hbar_eff - 1 / (2 * pi)) < mpf(10) ** -35, f"h_eff = {2 * pi * hbar_eff}")

# 7 damping: <exp(i e beta.phi)> = exp(-e^2 |beta|^2 C / 2), e^2 = 4 pi alpha, |beta|^2 = 2
def damp_exponent(C):         # coefficient k in exp(-k alpha)
    return 4 * pi * C
c.ok("one positive Dirac channel (weight 1/(2 pi)) gives exp(-2 alpha) exactly; one scalar channel gives exp(-alpha); both Dirac signs would give exp(-4 alpha)",
     abs(damp_exponent(wd[0]) - 2) < mpf(10) ** -35 and abs(damp_exponent(ws[0]) - 1) < mpf(10) ** -35)

# 8 kappa is not selected by the scalar kernel: weight 1/(4 pi kappa)
c.ok("scalar weight scales as 1/kappa (kappa = 1/2, 1, 2 -> 1/(2 pi), 1/(4 pi), 1/(8 pi)): kappa stays an input",
     all(abs(w_scalar(3, 3, k) - 1 / (4 * pi * k)) < mpf(10) ** -35 for k in (mpf(1) / 2, 1, 2)))

# 9 antipodal readout: sum_l (-1)^l * 1/(4 pi) -> Abel value 1/(8 pi) = c3
abel = lambda r: nsum(lambda l: (-r) ** l, [0, inf])
ab = [abel(mpf(1) - mpf(10) ** -k) for k in (2, 4, 6)]
c.ok("antipodal sum of the flat spectrum has Abel value 1/2: R*C(x,-x) = c3 (equals closed form 1/(8 pi))",
     abs(ab[-1] - mpf(1) / 2) < mpf(10) ** -6 and abs(altzeta(0) - mpf(1) / 2) < mpf(10) ** -30,
     f"Abel -> {ab[-1]}")

# 10 closed form check of the full kernel 1/(4 pi sqrt(2-2cos g)) against the channel sum at g = 2
g = mpf(2)
from mpmath import cos, sqrt
part = nsum(lambda l: 2 * C3 * legendre(int(l), cos(g)) * mpf("0.9") ** l, [0, 600])
closed = 2 * C3 / sqrt(1 - 2 * mpf("0.9") * cos(g) + mpf("0.9") ** 2)
c.ok("channel sum reproduces the Coulomb kernel 1/(4 pi |x-y|) (generating function, r = 0.9)",
     abs(part - closed) < mpf(10) ** -20, f"{part} vs {closed}")

# 11 zeta-regularised coincidence value is convention dependent: zeta(0) = -1/2, Hurwitz zeta(0, a) = 1/2 - a
hz = {a: zeta(0, a) for a in (mpf(1) / 2, 1, mpf(3) / 2)}
c.ok("regularised coincidence sum is NOT fixed: Hurwitz zeta(0,a) = 1/2 - a gives 0, -1/2, -1 for a = 1/2, 1, 3/2",
     all(abs(hz[a] - (mpf(1) / 2 - a)) < mpf(10) ** -30 for a in hz), hz)

# 12 alpha under alternative channel counts (damping exp(-k alpha)); only k = 2 survives
def ainv(k):
    return alpha_inv(lambda a: PHI_TREE + (lambda q: q * (1 - q) ** (-mpf(5) / 4))(DTOP * exp(-k * a)))
zs = {k: z(ainv(k)) for k in (0, 1, 2, 3, 4)}
c.ok("only two units of channel weight (k = 2) fit: k = 1 (one scalar channel) and k = 3, 4 are excluded",
     abs(zs[2]) < 3 and all(abs(zs[k]) > 5 for k in (0, 1, 3, 4)), {k: float(v) for k, v in zs.items()})

c.done({"z_by_k": {k: str(v) for k, v in zs.items()}, "w_scalar_S2": str(ws[0]), "w_dirac_S2": str(wd[0])})
